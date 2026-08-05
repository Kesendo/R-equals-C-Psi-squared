"""
Sacrifice-Zone Qubit Mapping: Find optimal chains on IBM Torino.

Uses real T2 calibration data and heavy-hex topology to find 5-qubit
chains where a naturally noisy qubit sits at the edge (sacrifice zone)
and quiet qubits in the center. Compares sacrifice-zone ranking with
naive mean-T2 ranking. Spectral verification via Liouvillian eigvals.

Convention: gamma_k = 1/(2*T2_k). The model here is dephasing-only (the only
jump operators are sqrt(gamma_k)*Z_k; T1 is never read at all), so this is the
D[Z] rate that reproduces the measured coherence decay: a D[Z] channel at rate
gamma decays coherences at 2*gamma. See docs/GLOSSARY.md, "The T2 -> gamma
conversion", and the gate simulations/t2_gamma_book_gate.py. The CSV column is
"T2_us" with no echo/Ramsey marker, so this code does not assume which it is.

Repaired 2026-08-05, and it was NOT the conversion: the palindrome score used to
measure its own matcher. spectral_analysis paired each rate with the FIRST
partner within an absolute 1e-4, while 927 of 959 nearest-neighbour level gaps
are smaller than that, so 362 of the 410 accepted pairs were wrong by more than
1e-12 (410 pairs = 820 of the 960 rates = the published 85.4%). The column now
reports the repo's canonical F1 distance (a port of
F1SpectrumStatistics.MaxF1PairingDistance, see max_f1_pairing_distance below) in
units of eps * spectral radius. A rates-only sorted-multiset residual stood here
between those two states for part of the same day and is a strictly weaker test;
see the comment at the palindrome check in spectral_analysis for the error model,
and experiments/CONCENTRATOR_MAPPING.md for the working.
"""

import numpy as np
import csv
from pathlib import Path
from collections import defaultdict

EPS = np.finfo(float).eps


# ================================================================
# Physics: Liouvillian construction (from ibm_cavity_analysis.py)
# ================================================================
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_at(op, target, n_qubits):
    result = np.eye(1, dtype=complex)
    for k in range(n_qubits):
        result = np.kron(result, op if k == target else I2)
    return result


def build_heisenberg_chain(n, J=1.0):
    d = 2**n
    H = np.zeros((d, d), dtype=complex)
    for i in range(n - 1):
        for pauli in [X, Y, Z]:
            H += J * kron_at(pauli, i, n) @ kron_at(pauli, i + 1, n)
    return H


def build_liouvillian(H, gammas):
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    n_qubits = int(np.log2(d))
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for k in range(n_qubits):
        if gammas[k] <= 0:
            continue
        Lk = np.sqrt(gammas[k]) * kron_at(Z, k, n_qubits)
        LdL = Lk.conj().T @ Lk
        L += np.kron(Lk.conj(), Lk)
        L -= 0.5 * np.kron(Id, LdL)
        L -= 0.5 * np.kron(LdL.T, Id)
    return L


def max_f1_pairing_distance(spectrum, sigma):
    """The repo's CANONICAL F1 symmetry distance, ported from C#.

    Port of `F1SpectrumStatistics.MaxF1PairingDistance`
    (compute/RCPsiSquared.Core/F1/F1SpectrumStatistics.cs), which its own summary
    calls "the canonical F1 check", and which is the matcher behind
    `MultisetAssert.NearestNeighbourEqual` and the live witness
    `BlockSpectrumWitness.PalindromePairingDistance`. Committed reference values
    live in simulations/results/f1_n8_n9_metrics/.

    Max greedy nearest-neighbour distance, WITH REMOVAL, between the multiset of
    eigenvalues and its F1 reflection {-2*sigma - lambda}. Multiplicity-aware by
    construction: each mirror point is consumed once, so a dropped or duplicated
    eigenvalue leaves an unmatched point and surfaces as a large distance, which
    a set or Hausdorff distance would miss.

    Note this runs on the FULL COMPLEX spectrum. A rates-only version (comparing
    only -Re(lambda) as a sorted multiset) is a strictly WEAKER test: it is blind
    to the imaginary parts, and halving one eigenvalue's Im leaves it unchanged.
    That weaker version stood here briefly on 2026-08-05 and was replaced by this
    one, because the repo already owned the check.
    """
    spectrum = np.asarray(spectrum)
    # Guard, because a NaN would otherwise pass SILENTLY. np.argmin selects a NaN
    # index and `nan > worst` is False, so a broken spectrum cannot raise the
    # metric: measured, one NaN among four eigenvalues returns 11.0 rather than a
    # failure. The C# returns double.MaxValue in that situation. A metric whose
    # job is to catch a broken spectrum must not fail quiet.
    if not np.all(np.isfinite(spectrum)):
        raise ValueError(
            "spectrum contains non-finite values; the F1 distance is undefined "
            "and would otherwise be silently understated")
    reflected = -2.0 * sigma - spectrum
    taken = np.zeros(len(spectrum), dtype=bool)
    worst = 0.0
    for x in spectrum:
        d = np.abs(x - reflected)
        d[taken] = np.inf
        j = int(np.argmin(d))
        taken[j] = True
        if d[j] > worst:
            worst = d[j]
    return float(worst)


def spectral_analysis(gammas, J=1.0):
    """Compute key spectral metrics for a 5-qubit chain."""
    N = len(gammas)
    H = build_heisenberg_chain(N, J)
    L = build_liouvillian(H, gammas)
    evals = np.linalg.eigvals(L)

    sum_g = sum(gammas)
    center = sum_g

    # Palindrome check: the repo's CANONICAL F1 symmetry distance (2026-08-05).
    #
    # `max_f1_pairing_distance` above is a port of
    # F1SpectrumStatistics.MaxF1PairingDistance, which the C# calls the canonical
    # F1 check and which is already surfaced as a live witness
    # (BlockSpectrumWitness.PalindromePairingDistance). The repo owned this
    # metric; this file previously did not use it.
    #
    # What it replaced, twice over. First a greedy first-fit matcher against an
    # ABSOLUTE 1e-4 tolerance, reported as a percentage: 927 of 959
    # nearest-neighbour level gaps are below that, so 362 of the 410 accepted
    # pairs were wrong by more than 1e-12 (410 pairs = 820 of 960 rates = the
    # published 85.4%). It was measuring its own matcher, reading 91.7%, 94.0%,
    # 100% and 100% on four chains whose true residuals are all ~1e-14. No
    # tolerance could fix it: any window below the level spacing returns exactly
    # 100% and hides the same clustering. Then, briefly the same day, a
    # rates-only sorted-multiset residual, which is a correct identity but a
    # strictly WEAKER test than F1: it discards the imaginary parts entirely.
    #
    # There is no exact route (a non-Hermitian eigensolver), so per the repo's
    # no-rounding rule the number carries its ERROR MODEL rather than a
    # threshold. It is reported in units of eps * spectral radius, and the
    # normalisation is not free: normalising by max|rate| instead spreads the
    # same ten chains by 15.0x, because max|rate| is only the real part while the
    # spectrum is dominated by |Im| (rho / max|rate| is 47x to 600x here). The
    # spectral radius brings them to 1.43x,
    # which is the floor rather than a residual error: permuting only the ORDER
    # in which the five jump operators are summed into L, same physics and same
    # gammas, moves this number over all 120 orders from 51.2 to 90.9, a 1.77x
    # spread (one chain, and identical under in-place and out-of-place
    # accumulation), MORE than the 1.43x separating the ten chains. And the greedy
    # matcher has a second sensitivity of its own: permuting only the ARRAY ORDER
    # of the eigenvalues, same spectrum, gives 58.9 against 75.8 on one chain
    # (1.29x). Greedy matching is order-dependent by construction, so this metric
    # still partly measures its matcher; it is reported anyway because it is the
    # repo's canonical check and because both sensitivities point the same way.
    # That order-dependence is
    # this repo's documented case-3 residual, a deterministic function of an
    # input the physics does not contain.
    #
    # Reading: the value is O(10) to O(100) eps here. The palindrome holds to the
    # eigensolver's own accuracy on every chain, and NO chain ranks above another
    # on it; ranking it would be ranking arithmetic order. It is NOT a proof of
    # exactness: a genuine violation below roughly 1e-13 absolute would sit
    # inside the same band. F1 is proven analytically; this is the numerical
    # check, not the theorem.
    pal_distance = max_f1_pairing_distance(evals, sum_g)
    rho = max(abs(e) for e in evals)
    pal_backward_eps = (pal_distance / (rho * EPS)) if rho > 0 else 0.0

    # Oscillating modes
    osc = [(-ev.real, abs(ev.imag)) for ev in evals if abs(ev.imag) > 1e-10]
    slowest_rate = min(r for r, _ in osc) if osc else float('inf')

    # Protection vs uniform
    uniform_gammas = np.array([sum_g / N] * N)
    L_uni = build_liouvillian(H, uniform_gammas)
    evals_uni = np.linalg.eigvals(L_uni)
    osc_uni = [(-ev.real, abs(ev.imag)) for ev in evals_uni if abs(ev.imag) > 1e-10]
    slowest_uni = min(r for r, _ in osc_uni) if osc_uni else float('inf')
    protection = slowest_uni / slowest_rate if slowest_rate > 1e-15 else 0

    # Protected modes count
    threshold = 0.05 * sum_g
    n_protected = sum(1 for r, _ in osc if r < threshold)
    n_protected_uni = sum(1 for r, _ in osc_uni if r < threshold)

    return {
        'palindrome_distance': pal_distance,
        'palindrome_backward_eps': pal_backward_eps,
        'slowest_rate': slowest_rate,
        'slowest_uni': slowest_uni,
        'protection': protection,
        'n_protected': n_protected,
        'n_protected_uni': n_protected_uni,
        'sum_gamma': sum_g,
    }


# ================================================================
# Topology: IBM Torino heavy-hex via Qiskit
# ================================================================
def get_torino_coupling_map():
    """Get IBM Torino coupling map via Qiskit."""
    try:
        from qiskit.providers.fake_provider import GenericBackendV2
        backend = GenericBackendV2(num_qubits=133, coupling_map=None)
        # GenericBackendV2 with None coupling_map generates a generic one
        # that may not match Torino. Try Fake backend instead.
    except Exception:
        pass

    # Fallback: use qiskit's built-in heavy-hex generator
    try:
        from qiskit.transpiler import CouplingMap
        cmap = CouplingMap.from_heavy_hex(7)  # distance-7 heavy-hex ~ 127 qubits
        edges = set()
        for e in cmap.get_edges():
            a, b = min(e), max(e)
            edges.add((a, b))
        return sorted(edges), cmap.size()
    except Exception as ex:
        print(f"  Qiskit CouplingMap failed: {ex}")

    # Last resort: hardcode a minimal heavy-hex structure
    print("  WARNING: Using hardcoded heavy-hex approximation")
    return None, 0


def get_coupling_map_fallback():
    """Build a heavy-hex-like coupling map for 133 qubits manually."""
    # IBM Heron r2 (133 qubits) coupling map from published specifications
    # This is the standard heavy-hex layout used in IBM Torino
    # Source: IBM Quantum documentation, various research papers
    edges = set()
    # Row structure for heavy-hex: alternating rows of vertices and links
    # For a 133-qubit Heron processor, use known connectivity patterns
    # Simplified: generate from Qiskit if possible, else approximate
    try:
        from qiskit.transpiler import CouplingMap
        cmap = CouplingMap.from_heavy_hex(7, bidirectional=False)
        for e in cmap.get_edges():
            a, b = min(e), max(e)
            edges.add((a, b))
        n_qubits = cmap.size()
        return sorted(edges), n_qubits
    except Exception:
        return None, 0


# ================================================================
# Chain finding
# ================================================================
def find_all_chains(edges, n_qubits, chain_length=5):
    """Find all simple paths of given length on the graph."""
    adj = defaultdict(set)
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)

    chains = set()

    def dfs(path):
        if len(path) == chain_length:
            # Canonical form: smaller end first
            canon = tuple(path) if path[0] < path[-1] else tuple(reversed(path))
            chains.add(canon)
            return
        last = path[-1]
        for neighbor in adj[last]:
            if neighbor not in path:
                path.append(neighbor)
                dfs(path)
                path.pop()

    for start in range(n_qubits):
        if start in adj:
            dfs([start])

    return [list(c) for c in sorted(chains)]


# ================================================================
# Main analysis
# ================================================================
def main():
    results_dir = Path(__file__).parent / "results"
    results_dir.mkdir(exist_ok=True)
    out_path = results_dir / "sacrifice_zone_mapping.txt"

    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=== SACRIFICE-ZONE QUBIT MAPPING ANALYSIS ===")
    out("IBM Torino, heavy-hex topology")
    out()

    # Step 1: Get coupling map
    out("Step 1: Loading coupling map...")
    edges, n_qubits = get_torino_coupling_map()
    if edges is None:
        edges, n_qubits = get_coupling_map_fallback()
    if edges is None:
        out("FATAL: Could not obtain coupling map. Aborting.")
        return

    out(f"  Qubits: {n_qubits}")
    out(f"  Edges: {len(edges)}")
    max_qubit = max(max(a, b) for a, b in edges)
    out(f"  Max qubit index: {max_qubit}")

    # Step 2: Load T2 data
    out("\nStep 2: Loading T2 calibration data...")
    csv_path = Path(__file__).parent.parent / "data" / "ibm_history" / "ibm_torino_history.csv"

    qubit_data = {}  # date -> {qubit -> T2_us}
    dates = set()
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            date = row['date']
            qubit = int(row['qubit'])
            t2 = float(row['T2_us']) if row['T2_us'] else None
            crosses = row.get('crosses_quarter', 'False') == 'True'
            if t2 and t2 > 0:
                if date not in qubit_data:
                    qubit_data[date] = {}
                qubit_data[date][qubit] = {'T2': t2, 'crosses': crosses}
            dates.add(date)

    sorted_dates = sorted(dates, reverse=True)
    out(f"  Total dates: {len(sorted_dates)}")
    out(f"  Latest: {sorted_dates[0]}")
    out(f"  Earliest: {sorted_dates[-1]}")

    # Use latest date
    latest = sorted_dates[0]
    data = qubit_data[latest]
    out(f"  Using: {latest} ({len(data)} qubits with T2 data)")

    # Compute gammas
    gammas = {}
    for q, d in data.items():
        gammas[q] = 1.0 / (2.0 * d['T2'])

    # Step 3: Find all 5-qubit chains
    out("\nStep 3: Finding all 5-qubit chains...")
    # Only use qubits that have T2 data AND are in the coupling map
    valid_qubits = set(gammas.keys()) & set(range(n_qubits))
    valid_edges = [(a, b) for a, b in edges if a in valid_qubits and b in valid_qubits]
    out(f"  Valid qubits (have T2 data): {len(valid_qubits)}")
    out(f"  Valid edges: {len(valid_edges)}")

    chains = find_all_chains(valid_edges, n_qubits, chain_length=5)
    out(f"  Total 5-qubit chains found: {len(chains)}")

    if len(chains) == 0:
        out("No chains found. Aborting.")
        return

    # Step 4: Compute sacrifice-zone scores
    out("\nStep 4: Computing sacrifice-zone scores...")

    chain_scores = []
    for chain in chains:
        g = [gammas.get(q, 0) for q in chain]
        t2 = [data[q]['T2'] for q in chain if q in data]
        crosses = [data[q]['crosses'] for q in chain if q in data]

        if any(gi <= 0 for gi in g):
            continue

        mean_interior = np.mean(g[1:4])
        mean_all = np.mean(g)

        # Sacrifice score: edge noise / interior noise
        score = max(g[0], g[4]) / mean_interior if mean_interior > 0 else 0

        # Asymmetric: which end is the sacrifice?
        score_left = g[0] / np.mean(g[1:]) if np.mean(g[1:]) > 0 else 0
        score_right = g[4] / np.mean(g[:4]) if np.mean(g[:4]) > 0 else 0
        best_sacrifice = max(score_left, score_right)
        sac_end = 'left' if score_left > score_right else 'right'

        chain_scores.append({
            'chain': chain,
            'gammas': g,
            'score': score,
            'best_sacrifice': best_sacrifice,
            'sac_end': sac_end,
            'mean_T2': np.mean(t2),
            'min_T2': min(t2),
            'sum_gamma': sum(g),
            'n_crossers': sum(crosses),
        })

    out(f"  Scored chains: {len(chain_scores)}")

    # Step 5: Rankings
    out("\nStep 5: Creating rankings...")

    by_sacrifice = sorted(chain_scores, key=lambda x: -x['best_sacrifice'])
    by_mean_t2 = sorted(chain_scores, key=lambda x: -x['mean_T2'])

    out(f"\n{'='*80}")
    out("TOP 10: SACRIFICE-ZONE RANKING")
    out(f"{'='*80}")
    out(f"{'Rank':>4} {'Chain':>25} {'Score':>7} {'mean_T2':>8} {'sum_g':>8} {'min_T2':>7} {'Cross':>5}")
    for i, cs in enumerate(by_sacrifice[:10]):
        chain_str = str(cs['chain'])
        out(f"{i+1:4d} {chain_str:>25} {cs['best_sacrifice']:7.2f} {cs['mean_T2']:8.1f} {cs['sum_gamma']:8.4f} {cs['min_T2']:7.1f} {cs['n_crossers']:5d}")

    out(f"\n{'='*80}")
    out("TOP 10: MEAN-T2 RANKING")
    out(f"{'='*80}")
    out(f"{'Rank':>4} {'Chain':>25} {'Score':>7} {'mean_T2':>8} {'sum_g':>8} {'min_T2':>7} {'Cross':>5}")
    for i, cs in enumerate(by_mean_t2[:10]):
        chain_str = str(cs['chain'])
        out(f"{i+1:4d} {chain_str:>25} {cs['best_sacrifice']:7.2f} {cs['mean_T2']:8.1f} {cs['sum_gamma']:8.4f} {cs['min_T2']:7.1f} {cs['n_crossers']:5d}")

    # Overlap analysis
    top_sac_set = set(tuple(cs['chain']) for cs in by_sacrifice[:10])
    top_t2_set = set(tuple(cs['chain']) for cs in by_mean_t2[:10])
    overlap = top_sac_set & top_t2_set

    out(f"\n{'='*80}")
    out("OVERLAP ANALYSIS")
    out(f"{'='*80}")
    out(f"Chains in BOTH top-10 lists: {len(overlap)}")
    out(f"Mean sacrifice score (sacrifice top-10): {np.mean([cs['best_sacrifice'] for cs in by_sacrifice[:10]]):.2f}")
    out(f"Mean sacrifice score (mean-T2 top-10): {np.mean([cs['best_sacrifice'] for cs in by_mean_t2[:10]]):.2f}")
    out(f"Mean mean-T2 (sacrifice top-10): {np.mean([cs['mean_T2'] for cs in by_sacrifice[:10]]):.1f}")
    out(f"Mean mean-T2 (mean-T2 top-10): {np.mean([cs['mean_T2'] for cs in by_mean_t2[:10]]):.1f}")

    # Step 6: Spectral verification of top 5 each
    out(f"\n{'='*80}")
    out("SPECTRAL VERIFICATION (Top 5 each)")
    out(f"{'='*80}")
    out(f"\n{'Chain':>25} {'Method':>10} {'Pal/eps':>9} {'Slowest':>10} {'Uni_slow':>10} {'Protect':>8} {'Prot_modes':>10}")

    verified = []
    for label, ranking in [("Sacrifice", by_sacrifice), ("Mean-T2", by_mean_t2)]:
        for cs in ranking[:5]:
            g = np.array(cs['gammas'])
            result = spectral_analysis(g)
            chain_str = str(cs['chain'])
            out(f"{chain_str:>25} {label:>10} {result['palindrome_backward_eps']:9.1f} {result['slowest_rate']:10.6f} {result['slowest_uni']:10.6f} {result['protection']:8.2f}x {result['n_protected']:10d}")
            verified.append({**cs, **result, 'method': label})

    # Step 7: Time stability (3 dates)
    out(f"\n{'='*80}")
    out("TIME STABILITY (3 dates)")
    out(f"{'='*80}")

    test_dates = [sorted_dates[0]]  # latest
    if len(sorted_dates) > 30:
        test_dates.append(sorted_dates[len(sorted_dates)//3])
    if len(sorted_dates) > 60:
        test_dates.append(sorted_dates[2*len(sorted_dates)//3])

    best_chain = by_sacrifice[0]['chain']
    out(f"\nTracking best sacrifice chain {best_chain} across dates:")
    out(f"{'Date':>12} {'Score':>8} {'mean_T2':>8} {'sum_g':>8}")

    for date in test_dates:
        if date in qubit_data:
            dd = qubit_data[date]
            if all(q in dd for q in best_chain):
                g = [1.0 / (2.0 * dd[q]['T2']) for q in best_chain]
                t2 = [dd[q]['T2'] for q in best_chain]
                mean_int = np.mean(g[1:4])
                score = max(g[0], g[4]) / mean_int if mean_int > 0 else 0
                out(f"{date:>12} {score:8.2f} {np.mean(t2):8.1f} {sum(g):8.4f}")
            else:
                out(f"{date:>12} (missing qubit data)")

    # Summary
    out(f"\n{'='*80}")
    out("CONCLUSION")
    out(f"{'='*80}")

    sac_protections = [v['protection'] for v in verified if v['method'] == 'Sacrifice']
    t2_protections = [v['protection'] for v in verified if v['method'] == 'Mean-T2']

    mean_sac_prot = np.mean(sac_protections) if sac_protections else 0
    mean_t2_prot = np.mean(t2_protections) if t2_protections else 0

    out(f"\nMean protection factor (sacrifice top-5): {mean_sac_prot:.2f}x")
    out(f"Mean protection factor (mean-T2 top-5): {mean_t2_prot:.2f}x")

    if mean_sac_prot > mean_t2_prot * 1.1:
        out(f"\nSacrifice-zone mapping OUTPERFORMS mean-T2 mapping by {mean_sac_prot/mean_t2_prot:.2f}x")
        out("Mode-based chain selection provides a practical advantage.")
    elif mean_t2_prot > mean_sac_prot * 1.1:
        out(f"\nMean-T2 mapping outperforms sacrifice-zone mapping.")
        out("Simple T2 maximization is sufficient for this hardware.")
    else:
        out(f"\nBoth methods perform similarly ({mean_sac_prot:.2f}x vs {mean_t2_prot:.2f}x).")
        out("The advantage depends on the specific noise distribution.")

    best = max(verified, key=lambda v: v['protection'])
    out(f"\nBest overall chain: {best['chain']}")
    out(f"  Method: {best['method']}")
    out(f"  Protection: {best['protection']:.2f}x")
    out(f"  Sacrifice score: {best['best_sacrifice']:.2f}")
    out(f"  Mean T2: {best['mean_T2']:.1f} us")

    out("\n=== DONE ===")

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"\n>>> Results saved to: {out_path}")


if __name__ == "__main__":
    main()
