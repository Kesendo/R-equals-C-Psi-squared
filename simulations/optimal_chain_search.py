"""
Optimal Chain Search: Find the sweet spot between contrast and total noise.

The Chain Selection Test showed:
  - High sacrifice score + high total noise (Chain A): wins for |+>^5 only
  - Low total noise + no contrast (Chain B): wins for |01010>
  - We need BOTH: enough contrast for mode protection, low enough total noise

This script searches all 330 chains on IBM Torino heavy-hex for the
optimal combined score, then validates the best candidate with full
spectral analysis and time evolution.
"""

import numpy as np
from scipy.linalg import expm, eig
from scipy.stats import pearsonr
import csv
import itertools
from pathlib import Path
from collections import defaultdict


# === Pauli matrices ===
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}


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


def partial_trace(rho, keep, n_qubits):
    dims = [2] * n_qubits
    rho_tensor = rho.reshape(dims + dims)
    trace_over = sorted([q for q in range(n_qubits) if q not in keep],
                        reverse=True)
    nq = n_qubits
    for q in trace_over:
        rho_tensor = np.trace(rho_tensor, axis1=q, axis2=q + nq)
        nq -= 1
    return rho_tensor.reshape(2**len(keep), 2**len(keep))


def von_neumann_entropy(rho):
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-15]
    return -np.sum(evals * np.log2(evals))


def mutual_information(rho, qi, qj, n_qubits):
    rho_ij = partial_trace(rho, [qi, qj], n_qubits)
    rho_i = partial_trace(rho, [qi], n_qubits)
    rho_j = partial_trace(rho, [qj], n_qubits)
    return (von_neumann_entropy(rho_i) + von_neumann_entropy(rho_j)
            - von_neumann_entropy(rho_ij))


def compute_observables(rho, n_qubits):
    purity = np.real(np.trace(rho @ rho))
    pairs = [(i, i+1) for i in range(n_qubits - 1)]
    mi_values = {}
    for qi, qj in pairs:
        mi_values[f"MI_{qi}{qj}"] = mutual_information(rho, qi, qj, n_qubits)
    return {'purity': purity, 'sum_mi': sum(mi_values.values()), **mi_values}


def make_plus_state(n):
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    psi = plus
    for _ in range(n - 1):
        psi = np.kron(psi, plus)
    return np.outer(psi, psi.conj())


def make_neel_state(n):
    d = 2**n
    idx = 0
    for k in range(n):
        bit = k % 2
        idx += bit * (2**(n - 1 - k))
    psi = np.zeros(d, dtype=complex)
    psi[idx] = 1.0
    return np.outer(psi, psi.conj())


def spectral_analysis(gammas, J=1.0):
    N = len(gammas)
    H = build_heisenberg_chain(N, J)
    L = build_liouvillian(H, gammas)
    evals = np.linalg.eigvals(L)

    center = sum(gammas)
    rates = sorted(-ev.real for ev in evals if abs(ev.imag) > 1e-10)
    paired = 0
    used = [False] * len(rates)
    for i in range(len(rates)):
        if used[i]:
            continue
        partner = 2 * center - rates[i]
        for j in range(len(rates)):
            if not used[j] and j != i and abs(rates[j] - partner) < 1e-4:
                paired += 2
                used[i] = True
                used[j] = True
                break
    pal = paired / max(len(rates), 1)

    osc = [(-ev.real, abs(ev.imag)) for ev in evals if abs(ev.imag) > 1e-10]
    slowest = min(r for r, _ in osc) if osc else float('inf')

    uniform_gammas = np.array([center / N] * N)
    L_uni = build_liouvillian(H, uniform_gammas)
    evals_uni = np.linalg.eigvals(L_uni)
    osc_uni = [(-ev.real, abs(ev.imag)) for ev in evals_uni
               if abs(ev.imag) > 1e-10]
    slowest_uni = min(r for r, _ in osc_uni) if osc_uni else float('inf')
    protection = slowest_uni / slowest if slowest > 1e-15 else 0

    return {'palindrome': pal, 'slowest': slowest, 'slowest_uni': slowest_uni,
            'protection': protection, 'sum_gamma': center}


def time_evolution(gammas, rho0, times, J=1.0):
    N = len(gammas)
    d = 2**N
    H = build_heisenberg_chain(N, J)
    L = build_liouvillian(H, gammas)
    rho0_vec = rho0.flatten()
    results = []
    for t in times:
        if t == 0:
            rho_vec = rho0_vec.copy()
        else:
            rho_vec = expm(L * t) @ rho0_vec
        rho = rho_vec.reshape(d, d)
        rho = (rho + rho.conj().T) / 2
        obs = compute_observables(rho, N)
        obs['t'] = t
        results.append(obs)
    return results


def get_coupling_map():
    """Get heavy-hex coupling map via Qiskit."""
    try:
        from qiskit.transpiler import CouplingMap
        cmap = CouplingMap.from_heavy_hex(7, bidirectional=False)
        edges = set()
        for e in cmap.get_edges():
            a, b = min(e), max(e)
            edges.add((a, b))
        return sorted(edges), cmap.size()
    except Exception as ex:
        print(f"Qiskit failed: {ex}")
        return None, 0


def find_all_chains(edges, n_qubits, chain_length=5):
    adj = defaultdict(set)
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    chains = set()

    def dfs(path):
        if len(path) == chain_length:
            canon = tuple(path) if path[0] < path[-1] else tuple(reversed(path))
            chains.add(canon)
            return
        for neighbor in adj[path[-1]]:
            if neighbor not in path:
                path.append(neighbor)
                dfs(path)
                path.pop()

    for start in range(n_qubits):
        if start in adj:
            dfs([start])
    return [list(c) for c in sorted(chains)]


def load_all_t2(csv_path):
    """Load T2 data for all qubits across all dates."""
    data = {}  # date -> {qubit -> {T2, T1}}
    dates = set()
    with open(csv_path, 'r') as f:
        for row in csv.DictReader(f):
            date = row['date']
            q = int(row['qubit'])
            t2 = float(row['T2_us']) if row['T2_us'] else None
            t1 = float(row['T1_us']) if row['T1_us'] else None
            dates.add(date)
            if t2 and t2 > 0 and t1 and t1 > 0:
                if date not in data:
                    data[date] = {}
                data[date][q] = {'T2': t2, 'T1': t1}
    return data, sorted(dates, reverse=True)


def main():
    results_dir = Path(__file__).parent / "results"
    results_dir.mkdir(exist_ok=True)
    out_path = results_dir / "optimal_chain_search.txt"

    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 70)
    out("OPTIMAL CHAIN SEARCH: Sweet spot between contrast and total noise")
    out("=" * 70)

    N = 5
    J = 1.0
    H = build_heisenberg_chain(N, J)

    # ================================================================
    # Step 1: Load data and enumerate chains
    # ================================================================
    out("\n--- STEP 1: Data and topology ---")
    csv_path = (Path(__file__).parent.parent / "data" / "ibm_history"
                / "ibm_torino_history.csv")
    all_data, sorted_dates = load_all_t2(csv_path)
    latest = sorted_dates[0]
    out(f"Latest date: {latest}, {len(sorted_dates)} dates total")

    edges, n_qubits = get_coupling_map()
    if edges is None:
        out("FATAL: No coupling map. Aborting.")
        return
    out(f"Topology: {n_qubits} qubits, {len(edges)} edges")

    data = all_data[latest]
    valid_qubits = set(data.keys()) & set(range(n_qubits))
    valid_edges = [(a, b) for a, b in edges
                   if a in valid_qubits and b in valid_qubits]
    chains = find_all_chains(valid_edges, n_qubits, 5)
    out(f"Chains found: {len(chains)}")

    # ================================================================
    # Step 2: Compute proxy scores for all chains
    # ================================================================
    out("\n--- STEP 2: Proxy scores for all chains ---")

    scored = []
    for chain in chains:
        g = np.array([1.0 / data[q]['T2'] for q in chain])
        t2 = np.array([data[q]['T2'] for q in chain])

        # Identify sacrifice end
        if g[0] >= g[4]:
            sac_idx, sac_gamma = 0, g[0]
            prot_gammas = g[1:]
            prot_t2 = t2[1:]
        else:
            sac_idx, sac_gamma = 4, g[4]
            prot_gammas = g[:4]
            prot_t2 = t2[:4]

        sacrifice_score = sac_gamma / np.mean(prot_gammas)
        sum_gamma = np.sum(g)
        mean_t2_prot = np.mean(prot_t2)
        contrast = np.max(g) / np.min(g)

        # Combined scores
        # A: protection per unit noise (higher = better)
        score_A = sacrifice_score / sum_gamma
        # B: protection * protected coherence (higher = better)
        score_B = sacrifice_score * mean_t2_prot
        # C: slowest mode estimate (1/sum_gamma as proxy, weighted by contrast)
        score_C = contrast / sum_gamma

        scored.append({
            'chain': chain,
            'gammas': g,
            'sacrifice_score': sacrifice_score,
            'sum_gamma': sum_gamma,
            'mean_t2': np.mean(t2),
            'mean_t2_prot': mean_t2_prot,
            'contrast': contrast,
            'sac_idx': sac_idx,
            'sac_t2': t2[sac_idx],
            'score_A': score_A,
            'score_B': score_B,
            'score_C': score_C,
        })

    # ================================================================
    # Step 3: Rankings
    # ================================================================
    out("\n--- STEP 3: Rankings by combined scores ---")

    rankings = {
        'A: sacrifice/sum_gamma': sorted(scored, key=lambda x: -x['score_A']),
        'B: sacrifice*mean_T2_prot': sorted(scored, key=lambda x: -x['score_B']),
        'C: contrast/sum_gamma': sorted(scored, key=lambda x: -x['score_C']),
    }

    for name, ranking in rankings.items():
        out(f"\n{'='*70}")
        out(f"TOP 10: {name}")
        out(f"{'='*70}")
        out(f"{'#':>3} {'Chain':>25} {'Score':>8} {'Sac':>6}"
            f" {'SumG':>8} {'MeanT2p':>8} {'Contr':>6} {'SacT2':>6}")
        for i, s in enumerate(ranking[:10]):
            sc = s['score_A'] if 'A' in name else (
                s['score_B'] if 'B' in name else s['score_C'])
            out(f"{i+1:3d} {str(s['chain']):>25} {sc:8.1f}"
                f" {s['sacrifice_score']:6.1f} {s['sum_gamma']:8.4f}"
                f" {s['mean_t2_prot']:8.1f} {s['contrast']:6.1f}"
                f" {s['sac_t2']:6.1f}")

    # ================================================================
    # Step 4: Full spectral analysis on top candidates
    # ================================================================
    out("\n--- STEP 4: Spectral verification of top candidates ---")

    # Collect unique top-10 chains from all rankings
    candidates = {}
    for name, ranking in rankings.items():
        for s in ranking[:10]:
            key = tuple(s['chain'])
            if key not in candidates:
                candidates[key] = s

    # Add reference chains
    ref_chains = {
        'Chain A (Sac)': [80, 8, 79, 53, 85],
        'Chain B (T2)': [18, 89, 19, 90, 60],
    }
    for label, chain in ref_chains.items():
        key = tuple(chain)
        if key not in candidates:
            g = np.array([1.0 / data[q]['T2'] for q in chain])
            candidates[key] = {
                'chain': chain, 'gammas': g,
                'sum_gamma': np.sum(g), 'mean_t2': np.mean(
                    [data[q]['T2'] for q in chain]),
                'sacrifice_score': 0, 'label': label,
            }

    out(f"\nComputing spectral analysis for {len(candidates)} unique chains...")
    verified = []
    for key, s in candidates.items():
        g = np.array([1.0 / data[q]['T2'] for q in s['chain']])
        spec = spectral_analysis(g, J)
        verified.append({**s, **spec})
        chain_str = str(s['chain'])[:25]
        out(f"  {chain_str:>25} prot={spec['protection']:.2f}x"
            f" slow={spec['slowest']:.6f} sumG={spec['sum_gamma']:.4f}")

    # Rank by absolute slowest rate (the actual thing we want to minimize)
    verified.sort(key=lambda x: x['slowest'])

    out(f"\n{'='*70}")
    out("RANKED BY ABSOLUTE SLOWEST MODE (lower = better)")
    out(f"{'='*70}")
    out(f"{'#':>3} {'Chain':>25} {'Slowest':>10} {'Prot':>6}"
        f" {'SumG':>8} {'Contr':>6} {'MeanT2':>7}")
    for i, v in enumerate(verified[:15]):
        out(f"{i+1:3d} {str(v['chain']):>25} {v['slowest']:10.6f}"
            f" {v['protection']:6.2f}x {v['sum_gamma']:8.4f}"
            f" {v.get('contrast', 0):6.1f} {v['mean_t2']:7.1f}")

    # Rank by protection / sum_gamma (combined efficiency)
    by_efficiency = sorted(verified,
                           key=lambda x: -x['protection'] / x['sum_gamma'])

    out(f"\n{'='*70}")
    out("RANKED BY PROTECTION EFFICIENCY (protection / sum_gamma)")
    out(f"{'='*70}")
    out(f"{'#':>3} {'Chain':>25} {'Effic':>8} {'Prot':>6}"
        f" {'SumG':>8} {'Slowest':>10} {'MeanT2':>7}")
    for i, v in enumerate(by_efficiency[:15]):
        eff = v['protection'] / v['sum_gamma']
        out(f"{i+1:3d} {str(v['chain']):>25} {eff:8.1f}"
            f" {v['protection']:6.2f}x {v['sum_gamma']:8.4f}"
            f" {v['slowest']:10.6f} {v['mean_t2']:7.1f}")

    # ================================================================
    # Step 5: Pick the winner and compare
    # ================================================================
    out("\n--- STEP 5: Winner selection ---")

    # The optimal chain minimizes the absolute slowest rate
    # while maintaining meaningful protection (> 1.5x)
    good_prot = [v for v in verified if v['protection'] > 1.3]
    if good_prot:
        winner = min(good_prot, key=lambda x: x['slowest'])
    else:
        winner = verified[0]

    chain_A_spec = next(v for v in verified
                        if tuple(v['chain']) == tuple([80, 8, 79, 53, 85]))
    chain_B_spec = next(v for v in verified
                        if tuple(v['chain']) == tuple([18, 89, 19, 90, 60]))

    out(f"\n{'Metric':>25} {'Winner':>12} {'Chain A':>12} {'Chain B':>12}")
    out(f"{'-'*25} {'-'*12} {'-'*12} {'-'*12}")
    for name, w, a, b in [
        ('Chain', str(winner['chain'])[:12], str(chain_A_spec['chain'])[:12],
         str(chain_B_spec['chain'])[:12]),
        ('Protection', f"{winner['protection']:.2f}x",
         f"{chain_A_spec['protection']:.2f}x",
         f"{chain_B_spec['protection']:.2f}x"),
        ('Sum(gamma)', f"{winner['sum_gamma']:.4f}",
         f"{chain_A_spec['sum_gamma']:.4f}",
         f"{chain_B_spec['sum_gamma']:.4f}"),
        ('Slowest rate', f"{winner['slowest']:.6f}",
         f"{chain_A_spec['slowest']:.6f}",
         f"{chain_B_spec['slowest']:.6f}"),
        ('Mean T2 (us)', f"{winner['mean_t2']:.1f}",
         f"{chain_A_spec['mean_t2']:.1f}",
         f"{chain_B_spec['mean_t2']:.1f}"),
    ]:
        out(f"{name:>25} {w:>12} {a:>12} {b:>12}")

    # ================================================================
    # Step 6: Time evolution comparison
    # ================================================================
    out("\n--- STEP 6: Time evolution ---")

    times = np.linspace(0, 10, 51)
    rho_plus = make_plus_state(N)
    rho_neel = make_neel_state(N)

    gammas_W = np.array([1.0 / data[q]['T2'] for q in winner['chain']])
    gammas_A = np.array([1.0 / data[q]['T2'] for q in [80, 8, 79, 53, 85]])
    gammas_B = np.array([1.0 / data[q]['T2'] for q in [18, 89, 19, 90, 60]])

    out(f"\nWinner: {winner['chain']}")
    out(f"gammas: [{', '.join(f'{g:.6f}' for g in gammas_W)}]")

    # |+>^5
    out("\nComputing |+>^5 evolution (3 chains)...")
    res_W_plus = time_evolution(gammas_W, rho_plus, times, J)
    res_A_plus = time_evolution(gammas_A, rho_plus, times, J)
    res_B_plus = time_evolution(gammas_B, rho_plus, times, J)

    out(f"\n--- SumMI(t): |+>^5 ---")
    out(f"{'t':>5} {'Winner':>10} {'ChainA':>10} {'ChainB':>10}"
        f" {'W/A':>6} {'W/B':>6}")
    for rW, rA, rB in zip(res_W_plus, res_A_plus, res_B_plus):
        t = rW['t']
        if t % 2.0 < 0.01 or t < 0.5:
            wa = (rW['sum_mi'] / rA['sum_mi']
                  if rA['sum_mi'] > 1e-8 else float('inf'))
            wb = (rW['sum_mi'] / rB['sum_mi']
                  if rB['sum_mi'] > 1e-8 else float('inf'))
            out(f"{t:5.1f} {rW['sum_mi']:10.6f} {rA['sum_mi']:10.6f}"
                f" {rB['sum_mi']:10.6f} {wa:6.1f}x {wb:6.1f}x")

    # |01010>
    out("\nComputing |01010> evolution (3 chains)...")
    res_W_neel = time_evolution(gammas_W, rho_neel, times, J)
    res_A_neel = time_evolution(gammas_A, rho_neel, times, J)
    res_B_neel = time_evolution(gammas_B, rho_neel, times, J)

    out(f"\n--- SumMI(t): |01010> ---")
    out(f"{'t':>5} {'Winner':>10} {'ChainA':>10} {'ChainB':>10}"
        f" {'W/A':>6} {'W/B':>6}")
    for rW, rA, rB in zip(res_W_neel, res_A_neel, res_B_neel):
        t = rW['t']
        if t % 2.0 < 0.01 or t < 0.5:
            wa = (rW['sum_mi'] / rA['sum_mi']
                  if rA['sum_mi'] > 1e-8 else float('inf'))
            wb = (rW['sum_mi'] / rB['sum_mi']
                  if rB['sum_mi'] > 1e-8 else float('inf'))
            out(f"{t:5.1f} {rW['sum_mi']:10.6f} {rA['sum_mi']:10.6f}"
                f" {rB['sum_mi']:10.6f} {wa:6.1f}x {wb:6.1f}x")

    # ================================================================
    # Step 7: Selective DD on the winner
    # ================================================================
    out("\n--- STEP 7: Selective DD on winner ---")

    # DD reduces effective gamma to ~1/T2echo on protected qubits
    # Without DD: gamma = 1/T2echo (already our convention)
    # With DD: gamma_protected ~ 1/(2*T2echo) (factor 2 reduction)
    # Sacrifice qubit: no DD, stays at 1/T2echo
    sac_end = 0 if gammas_W[0] > gammas_W[4] else 4

    gammas_W_dd = gammas_W.copy()
    for k in range(N):
        if k != sac_end:
            gammas_W_dd[k] *= 0.5  # DD reduces gamma by ~2x on protected qubits

    out(f"Sacrifice end: position {sac_end} (Q{winner['chain'][sac_end]})")
    out(f"DD gammas: [{', '.join(f'{g:.6f}' for g in gammas_W_dd)}]")

    spec_dd = spectral_analysis(gammas_W_dd, J)
    out(f"Protection with DD: {spec_dd['protection']:.2f}x"
        f" (was {winner['protection']:.2f}x without)")
    out(f"Slowest rate with DD: {spec_dd['slowest']:.6f}"
        f" (was {winner['slowest']:.6f} without)")

    res_W_dd_plus = time_evolution(gammas_W_dd, rho_plus, times, J)
    res_W_dd_neel = time_evolution(gammas_W_dd, rho_neel, times, J)

    peak_W_plus = max(res_W_plus, key=lambda r: r['sum_mi'])
    peak_W_neel = max(res_W_neel, key=lambda r: r['sum_mi'])
    peak_W_dd_plus = max(res_W_dd_plus, key=lambda r: r['sum_mi'])
    peak_W_dd_neel = max(res_W_dd_neel, key=lambda r: r['sum_mi'])
    peak_A_neel = max(res_A_neel, key=lambda r: r['sum_mi'])
    peak_B_neel = max(res_B_neel, key=lambda r: r['sum_mi'])

    out(f"\n--- Peak SumMI comparison ---")
    out(f"{'Config':>25} {'|+>^5':>10} {'|01010>':>10}")
    out(f"{'Winner (no DD)':>25} {peak_W_plus['sum_mi']:10.6f}"
        f" {peak_W_neel['sum_mi']:10.6f}")
    out(f"{'Winner (sel. DD)':>25} {peak_W_dd_plus['sum_mi']:10.6f}"
        f" {peak_W_dd_neel['sum_mi']:10.6f}")
    out(f"{'Chain A (no DD)':>25}"
        f" {max(res_A_plus, key=lambda r: r['sum_mi'])['sum_mi']:10.6f}"
        f" {peak_A_neel['sum_mi']:10.6f}")
    out(f"{'Chain B (no DD)':>25}"
        f" {max(res_B_plus, key=lambda r: r['sum_mi'])['sum_mi']:10.6f}"
        f" {peak_B_neel['sum_mi']:10.6f}")

    # ================================================================
    # Step 8: Temporal stability
    # ================================================================
    out("\n--- STEP 8: Temporal stability ---")

    test_dates = [sorted_dates[0]]
    n_dates = len(sorted_dates)
    for frac in [0.25, 0.5, 0.75]:
        idx = int(frac * n_dates)
        if idx < n_dates:
            test_dates.append(sorted_dates[idx])
    if n_dates > 1:
        test_dates.append(sorted_dates[-1])

    out(f"\nTracking winner {winner['chain']} across {len(test_dates)} dates:")
    out(f"{'Date':>12} {'SumG':>8} {'Contr':>6} {'SacT2':>7}"
        f" {'ProtT2':>7}")

    for date in test_dates:
        dd = all_data.get(date, {})
        if all(q in dd for q in winner['chain']):
            g = [1.0 / dd[q]['T2'] for q in winner['chain']]
            t2 = [dd[q]['T2'] for q in winner['chain']]
            sum_g = sum(g)
            contrast = max(g) / min(g)
            sac_t2 = t2[sac_end]
            prot_t2 = np.mean([t2[k] for k in range(N) if k != sac_end])
            out(f"{date:>12} {sum_g:8.4f} {contrast:6.1f} {sac_t2:7.1f}"
                f" {prot_t2:7.1f}")
        else:
            out(f"{date:>12} (missing data)")

    # ================================================================
    # Step 9: Final recommendation
    # ================================================================
    out("\n" + "=" * 70)
    out("RECOMMENDATION")
    out("=" * 70)

    out(f"\nOptimal chain: {winner['chain']}")
    out(f"  Protection factor: {winner['protection']:.2f}x")
    out(f"  Sum(gamma): {winner['sum_gamma']:.4f}")
    out(f"  Slowest mode rate: {winner['slowest']:.6f}")
    out(f"  Mean T2: {winner['mean_t2']:.1f} us")
    out(f"  Contrast: {winner.get('contrast', 0):.1f}x")

    out(f"\nWith selective DD:")
    out(f"  Protection factor: {spec_dd['protection']:.2f}x")
    out(f"  Slowest mode rate: {spec_dd['slowest']:.6f}")

    out(f"\nComparison (peak SumMI, |01010>):")
    out(f"  Winner:  {peak_W_neel['sum_mi']:.4f}")
    out(f"  +DD:     {peak_W_dd_neel['sum_mi']:.4f}")
    out(f"  Chain A: {peak_A_neel['sum_mi']:.4f}")
    out(f"  Chain B: {peak_B_neel['sum_mi']:.4f}")

    # Determine the actual best configuration
    configs = [
        ('Winner+DD', peak_W_dd_neel['sum_mi']),
        ('Winner', peak_W_neel['sum_mi']),
        ('Chain A', peak_A_neel['sum_mi']),
        ('Chain B', peak_B_neel['sum_mi']),
    ]
    configs.sort(key=lambda x: -x[1])
    out(f"\nBest config for |01010>: {configs[0][0]} ({configs[0][1]:.4f})")

    out("\n=== DONE ===")

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"\n>>> Results saved to: {out_path}")


if __name__ == "__main__":
    main()
