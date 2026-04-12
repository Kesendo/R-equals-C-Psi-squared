"""
Symmetry Census: enumerating the block structure of the Liouvillian
====================================================================
For N-qubit Heisenberg chain under Z-dephasing, systematically enumerate:
  1. All symmetries that commute with the Liouvillian superoperator
  2. The sector decomposition (simultaneous eigenspaces of commuting symmetries)
  3. Asymptotic attractors per sector (steady states, limit cycles)
  4. Which sectors are reachable from physical initial states
  5. How the picture changes across topologies (chain, ring, star, complete)

Symmetries checked:
  - U(1) excitation number conservation: quantum numbers (w_bra, w_ket)
  - n_XY parity: popcount(i XOR j) mod 2 for basis element |i><j|
  - Spin-flip X^{otimes N}: maps |i><j| to |~i><~j|, pairs (w,w') <-> (N-w,N-w')
  - Spatial reflection (uniform gamma only): maps site k to N-1-k
  - Pi conjugation: spectral mirror (L Pi = -Pi (L + 2 Sigma_gamma I)),
    documented in docs/proofs/MIRROR_SYMMETRY_PROOF.md. Not a block-
    diagonalizing symmetry but a spectral pairing.

The finest block structure from commuting symmetries is:
  Non-uniform gamma: sectors = (w_bra, w_ket, nxy_parity)
  Uniform gamma: above + spatial_reflection_eigenvalue within each sector

Authors: Thomas Wicht, Claude (Opus 4.6)
Date: April 12, 2026
"""
import os
import sys
sys.stdout.reconfigure(encoding='utf-8')
sys.stderr.reconfigure(encoding='utf-8')
import json
import time as _time
from collections import defaultdict
import numpy as np
from scipy import linalg

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ibm_april_predictions import heisenberg_H, build_liouvillian

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "results", "symmetry_census")
os.makedirs(OUT_DIR, exist_ok=True)


# ===================================================================
# Quantum number assignment for the operator basis |i><j|
# ===================================================================
def popcount(x):
    return bin(x).count('1')


def bit_reverse(x, N):
    """Reverse the bit order of x (N-bit integer)."""
    result = 0
    for k in range(N):
        if x & (1 << k):
            result |= 1 << (N - 1 - k)
    return result


def assign_quantum_numbers(N):
    """For each operator basis element |i><j|, compute quantum numbers.

    Returns dict mapping (i, j) -> dict of quantum numbers.
    Also returns sector_map: (w_bra, w_ket, nxy_par) -> list of (i,j).
    """
    D = 2 ** N
    mask = D - 1
    sector_map = defaultdict(list)
    qn_map = {}

    for i in range(D):
        for j in range(D):
            w_bra = popcount(i)
            w_ket = popcount(j)
            nxy_par = popcount(i ^ j) % 2
            flip_i = (~i) & mask
            flip_j = (~j) & mask
            ref_i = bit_reverse(i, N)
            ref_j = bit_reverse(j, N)

            qn = dict(
                w_bra=w_bra, w_ket=w_ket, nxy_par=nxy_par,
                flip_image=(flip_i, flip_j),
                ref_image=(ref_i, ref_j),
            )
            qn_map[(i, j)] = qn
            sector_map[(w_bra, w_ket, nxy_par)].append((i, j))

    return qn_map, sector_map


# ===================================================================
# Sector enumeration
# ===================================================================
def enumerate_sectors(N):
    """Count sectors and their dimensions for chain of length N.

    Returns list of dicts with sector info.
    """
    _, sector_map = assign_quantum_numbers(N)
    sectors = []
    for (w, wp, p), elements in sorted(sector_map.items()):
        sectors.append(dict(
            w_bra=w, w_ket=wp, nxy_par=p,
            dim=len(elements),
            flip_partner=(N - w, N - wp, p),
            is_self_conjugate=(w == N - w and wp == N - wp),
        ))
    return sectors


def sector_summary(N):
    """Print sector count summary for one N."""
    sectors = enumerate_sectors(N)
    total_dim = sum(s['dim'] for s in sectors)
    n_nonzero = len([s for s in sectors if s['dim'] > 0])

    # Spin-flip pairing
    paired = set()
    n_sf_blocks = 0
    for s in sectors:
        key = (s['w_bra'], s['w_ket'], s['nxy_par'])
        partner = s['flip_partner']
        if key in paired:
            continue
        if key == partner:
            # Self-conjugate: splits into + and -
            n_sf_blocks += 2
            paired.add(key)
        else:
            # Paired: combined space splits into + and -
            n_sf_blocks += 2
            paired.add(key)
            paired.add(partner)

    return dict(
        N=N, D=2**N, d2=(2**N)**2,
        n_basic_sectors=n_nonzero,
        n_sf_blocks=n_sf_blocks,
        total_dim=total_dim,
        sectors=sectors,
    )


# ===================================================================
# Build Liouvillian restricted to a sector
# ===================================================================
def build_sector_L(N, gammas, sector_elements, J=1.0):
    """Build the Liouvillian restricted to the given sector.

    sector_elements: list of (i, j) basis element pairs.
    Returns (L_sector, element_list).
    """
    D = 2 ** N
    # Build full Liouvillian (only for N <= 6 to keep tractable)
    H = heisenberg_H(N, J)
    L_full = build_liouvillian(H, gammas)

    M = len(sector_elements)
    # Map (i, j) -> flat index in full Liouvillian
    elem_to_flat = {(i, j): i * D + j for i, j in sector_elements}
    flat_indices = [elem_to_flat[e] for e in sector_elements]

    L_sector = L_full[np.ix_(flat_indices, flat_indices)]
    return L_sector, sector_elements


# ===================================================================
# Attractor analysis per sector
# ===================================================================
def analyze_sector_attractors(N, gammas, J=1.0):
    """For each sector at given N, find steady states and characterize dynamics.

    Returns list of sector results.
    """
    _, sector_map = assign_quantum_numbers(N)
    results = []

    for (w, wp, p), elements in sorted(sector_map.items()):
        dim = len(elements)
        if dim == 0:
            continue

        L_sec, _ = build_sector_L(N, gammas, elements, J)
        eigvals = linalg.eigvals(L_sec)

        # Steady states: Re(lambda) ≈ 0
        n_steady = int(np.sum(np.abs(eigvals.real) < 1e-10))
        # Oscillatory modes: nonzero Im
        n_osc = int(np.sum((np.abs(eigvals.imag) > 1e-10) &
                           (np.abs(eigvals.real) > 1e-10)))
        # Purely decaying
        n_decay = dim - n_steady - n_osc

        # Slowest non-stationary rate
        nonstat = np.abs(eigvals.real) > 1e-10
        if nonstat.any():
            slowest = float(np.min(np.abs(eigvals[nonstat].real)))
        else:
            slowest = None

        results.append(dict(
            w_bra=w, w_ket=wp, nxy_par=p, dim=dim,
            n_steady=n_steady, n_osc=n_osc, n_decay=n_decay,
            slowest_rate=slowest,
        ))

    return results


# ===================================================================
# Degeneracy check (are known symmetries sufficient?)
# ===================================================================
def check_degeneracies(N, gammas, J=1.0):
    """Check if all eigenvalue degeneracies are explained by known symmetries.

    Build the full Liouvillian, compute its spectrum, group eigenvalues,
    and check if degeneracies match the sector pairing predictions.
    """
    D = 2 ** N
    H = heisenberg_H(N, J)
    L = build_liouvillian(H, gammas)
    eigvals = linalg.eigvals(L)

    # Sort by real part, then imaginary
    order = np.lexsort((eigvals.imag, eigvals.real))
    eigvals = eigvals[order]

    # Group nearly-equal eigenvalues
    tol = 1e-8
    groups = []
    used = np.zeros(len(eigvals), dtype=bool)
    for i in range(len(eigvals)):
        if used[i]:
            continue
        group = [i]
        for j in range(i + 1, len(eigvals)):
            if used[j]:
                continue
            if abs(eigvals[i] - eigvals[j]) < tol:
                group.append(j)
                used[j] = True
        used[i] = True
        groups.append(dict(
            eigenvalue=complex(eigvals[i]),
            multiplicity=len(group),
        ))

    # Check palindromic pairing (Pi conjugation)
    Sg = float(np.sum(gammas))
    center = -Sg
    rates = -eigvals.real
    paired = 0
    unpaired = 0
    for g in groups:
        ev = g['eigenvalue']
        mirror = complex(-ev.real - 2 * Sg, -ev.imag)
        found = False
        for g2 in groups:
            if abs(g2['eigenvalue'] - mirror) < tol:
                found = True
                break
        if found:
            paired += g['multiplicity']
        else:
            unpaired += g['multiplicity']

    return dict(
        n_eigenvalues=len(eigvals),
        n_distinct=len(groups),
        max_multiplicity=max(g['multiplicity'] for g in groups),
        palindromic_paired=paired,
        unpaired=unpaired,
        groups=groups,
    )


# ===================================================================
# Reachability from physical initial states
# ===================================================================
def state_sector_population(psi, N):
    """For a pure state |psi>, compute which sectors rho_0 = |psi><psi| populates.

    Returns dict of (w_bra, w_ket, nxy_par) -> weight.
    """
    D = 2 ** N
    rho = np.outer(psi, psi.conj())
    pop = defaultdict(float)
    for i in range(D):
        for j in range(D):
            val = abs(rho[i, j]) ** 2
            if val < 1e-30:
                continue
            w = popcount(i)
            wp = popcount(j)
            p = popcount(i ^ j) % 2
            pop[(w, wp, p)] += val
    # Normalize
    total = sum(pop.values())
    return {k: v / total for k, v in pop.items()}


def reachability_table(N):
    """Compute sector populations for common physical initial states."""
    D = 2 ** N
    states = {}

    # Product state |00...0>
    psi = np.zeros(D, dtype=complex)
    psi[0] = 1.0
    states['|0>^N'] = psi

    # Product state |11...1>
    psi = np.zeros(D, dtype=complex)
    psi[D - 1] = 1.0
    states['|1>^N'] = psi

    # |+>^N
    psi = np.ones(D, dtype=complex) / np.sqrt(D)
    states['|+>^N'] = psi

    # W_N state
    psi = np.zeros(D, dtype=complex)
    for k in range(N):
        psi[1 << (N - 1 - k)] = 1.0 / np.sqrt(N)
    states['W_N'] = psi

    # GHZ state
    psi = np.zeros(D, dtype=complex)
    psi[0] = 1.0 / np.sqrt(2)
    psi[D - 1] = 1.0 / np.sqrt(2)
    states['GHZ'] = psi

    # Bell+ on central pair, |0> elsewhere
    c1 = (N - 1) // 2
    c2 = c1 + 1
    psi = np.zeros(D, dtype=complex)
    psi[0] = 1.0 / np.sqrt(2)
    idx_11 = (1 << (N - 1 - c1)) | (1 << (N - 1 - c2))
    psi[idx_11] = 1.0 / np.sqrt(2)
    states[f'Bell+({c1},{c2})|0>'] = psi

    # Neel state |01010...>
    neel_idx = 0
    for k in range(N):
        if k % 2 == 1:
            neel_idx |= 1 << (N - 1 - k)
    psi = np.zeros(D, dtype=complex)
    psi[neel_idx] = 1.0
    states['Neel'] = psi

    results = {}
    for name, psi in states.items():
        pop = state_sector_population(psi, N)
        n_sectors = len(pop)
        # Collect unique (w_bra, w_ket) pairs
        ww_pairs = set((k[0], k[1]) for k in pop.keys())
        results[name] = dict(
            n_sectors=n_sectors,
            n_ww_pairs=len(ww_pairs),
            sectors=pop,
        )

    return results


# ===================================================================
# Topology comparison
# ===================================================================
def topology_bonds(N, topo):
    """Return bond list for given topology."""
    if topo == 'chain':
        return [(i, i + 1) for i in range(N - 1)]
    elif topo == 'ring':
        return [(i, (i + 1) % N) for i in range(N)]
    elif topo == 'star':
        return [(0, i) for i in range(1, N)]
    elif topo == 'complete':
        return [(i, j) for i in range(N) for j in range(i + 1, N)]
    raise ValueError(f"Unknown topology: {topo}")


def heisenberg_H_topo(N, bonds, J=1.0):
    """Build Heisenberg Hamiltonian for arbitrary topology."""
    I2 = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    D = 2 ** N
    H = np.zeros((D, D), dtype=complex)
    for (i, j) in bonds:
        for pauli in [X, Y, Z]:
            ops = [I2] * N
            ops[i] = pauli
            ops[j] = pauli
            term = ops[0]
            for op in ops[1:]:
                term = np.kron(term, op)
            H += J * term
    return H


def topology_comparison(N):
    """Compare sector structure and spectrum across topologies."""
    gammas = np.ones(N) * 0.1  # uniform for comparison
    results = {}

    for topo in ['chain', 'ring', 'star', 'complete']:
        bonds = topology_bonds(N, topo)
        H = heisenberg_H_topo(N, bonds)
        L = build_liouvillian(H, gammas)
        eigvals = linalg.eigvals(L)

        # Spectrum statistics
        rates = -eigvals.real
        Sg = float(gammas.sum())
        n_stat = int(np.sum(np.abs(rates) < 1e-10))
        n_osc = int(np.sum(np.abs(eigvals.imag) > 1e-10))

        # Group distinct eigenvalues
        tol = 1e-8
        distinct = []
        used = np.zeros(len(eigvals), dtype=bool)
        for i in range(len(eigvals)):
            if used[i]:
                continue
            mult = 1
            for j in range(i + 1, len(eigvals)):
                if not used[j] and abs(eigvals[i] - eigvals[j]) < tol:
                    mult += 1
                    used[j] = True
            used[i] = True
            distinct.append(mult)

        max_deg = max(distinct) if distinct else 0

        # Additional symmetries: ring has cyclic C_N
        extra_sym = []
        if topo == 'ring':
            extra_sym.append(f'C_{N} cyclic')
        elif topo == 'star':
            extra_sym.append(f'S_{N-1} leaf permutation')
        elif topo == 'complete':
            extra_sym.append(f'S_{N} full permutation')

        results[topo] = dict(
            n_bonds=len(bonds),
            n_eigenvalues=len(eigvals),
            n_stationary=n_stat,
            n_oscillatory=n_osc,
            n_distinct=len(distinct),
            max_degeneracy=max_deg,
            extra_symmetries=extra_sym,
        )

    return results


# ===================================================================
# Main
# ===================================================================
if __name__ == "__main__":
    print("Symmetry Census: Liouvillian Block Structure")
    print("=" * 60)

    # ============================================================
    # Q1 + Q2: Sector enumeration for N=3-7
    # ============================================================
    print("\n--- SECTOR ENUMERATION ---")
    all_sector_data = {}
    for N in [3, 4, 5, 6, 7]:
        info = sector_summary(N)
        all_sector_data[N] = info
        print(f"  N={N}: D={info['D']}, d²={info['d2']}, "
              f"basic sectors={info['n_basic_sectors']}, "
              f"with spin-flip={info['n_sf_blocks']}")

    # Detailed sector table for N=5
    print("\n--- N=5 SECTOR TABLE ---")
    sectors_5 = all_sector_data[5]['sectors']
    print(f"  {"(w,w')":<8} {'par':>3} {'dim':>5} {'flip':>10} {'self?':>5}")
    print(f"  {'-' * 35}")
    for s in sectors_5:
        if s['dim'] > 0:
            fp = f"({s['flip_partner'][0]},{s['flip_partner'][1]})"
            print(f"  ({s['w_bra']},{s['w_ket']}){'':<3} {s['nxy_par']:>3} "
                  f"{s['dim']:>5} {fp:>10} "
                  f"{'*' if s['is_self_conjugate'] else '':>5}")

    # ============================================================
    # Q3: Attractor analysis for N=5 (uniform gamma)
    # ============================================================
    print("\n--- ATTRACTOR ANALYSIS (N=5, uniform gamma=0.1) ---")
    N_detail = 5
    gamma_uniform = np.ones(N_detail) * 0.1
    t0 = _time.time()
    attractors = analyze_sector_attractors(N_detail, gamma_uniform)
    print(f"  Computed in {_time.time() - t0:.1f}s")
    print(f"\n  {"(w,w')":>8} {'p':>2} {'dim':>5} {'steady':>6} "
          f"{'osc':>5} {'decay':>5} {'slowest':>10}")
    print(f"  {'-' * 50}")
    total_steady = 0
    for a in attractors:
        total_steady += a['n_steady']
        sr = f"{a['slowest_rate']:.4f}" if a['slowest_rate'] else "n/a"
        print(f"  ({a['w_bra']},{a['w_ket']}){'':<3} {a['nxy_par']:>2} "
              f"{a['dim']:>5} {a['n_steady']:>6} "
              f"{a['n_osc']:>5} {a['n_decay']:>5} {sr:>10}")
    print(f"\n  Total steady states: {total_steady}")
    print(f"  (Expected from sector conservation: one per diagonal sector "
          f"(w,w,even) = {N_detail + 1})")

    # ============================================================
    # Q1 continued: Degeneracy check (N=5, uniform)
    # ============================================================
    print("\n--- DEGENERACY CHECK (N=5, uniform gamma=0.1) ---")
    deg_info = check_degeneracies(N_detail, gamma_uniform)
    print(f"  Total eigenvalues: {deg_info['n_eigenvalues']}")
    print(f"  Distinct eigenvalues: {deg_info['n_distinct']}")
    print(f"  Max multiplicity: {deg_info['max_multiplicity']}")
    print(f"  Palindromic paired: {deg_info['palindromic_paired']}")
    print(f"  Unpaired: {deg_info['unpaired']}")

    # High-multiplicity groups
    high_deg = [g for g in deg_info['groups'] if g['multiplicity'] > 2]
    if high_deg:
        print(f"\n  High-multiplicity eigenvalues (>2):")
        for g in sorted(high_deg, key=lambda x: -x['multiplicity'])[:15]:
            ev = g['eigenvalue']
            print(f"    Re={ev.real:+.6f} Im={ev.imag:+.6f}  "
                  f"mult={g['multiplicity']}")

    # Also check sacrifice profile
    print("\n--- DEGENERACY CHECK (N=5, IBM sacrifice) ---")
    gamma_sac = np.array([2.336, 0.099, 0.050, 0.072, 0.051])
    # Rescale to match sum
    gamma_sac = gamma_sac / gamma_sac.sum() * gamma_uniform.sum()
    deg_sac = check_degeneracies(N_detail, gamma_sac)
    print(f"  Distinct eigenvalues: {deg_sac['n_distinct']}")
    print(f"  Max multiplicity: {deg_sac['max_multiplicity']}")
    high_sac = [g for g in deg_sac['groups'] if g['multiplicity'] > 1]
    print(f"  Eigenvalues with multiplicity > 1: {len(high_sac)}")

    # ============================================================
    # Q4: Reachability
    # ============================================================
    print("\n--- REACHABILITY TABLE (N=5) ---")
    reach = reachability_table(N_detail)
    print(f"  {'State':<20} {'#sectors':>8} {'#(w,w)':>8}")
    print(f"  {'-' * 40}")
    for name, info in reach.items():
        print(f"  {name:<20} {info['n_sectors']:>8} {info['n_ww_pairs']:>8}")

    # Detailed sector populations for a few states
    for name in ['W_N', 'GHZ', '|+>^N']:
        print(f"\n  {name} sector populations:")
        for (w, wp, p), weight in sorted(reach[name]['sectors'].items(),
                                         key=lambda x: -x[1]):
            if weight > 0.001:
                print(f"    ({w},{wp}) par={p}: {weight:.4f}")

    # ============================================================
    # Q5: Topology comparison (N=4)
    # ============================================================
    print("\n--- TOPOLOGY COMPARISON (N=4, uniform gamma=0.1) ---")
    topo_results = topology_comparison(4)
    print(f"  {'Topo':<10} {'bonds':>5} {'stationary':>10} "
          f"{'oscillatory':>11} {'distinct':>8} {'max_deg':>7} {'extra':>20}")
    print(f"  {'-' * 75}")
    for topo, info in topo_results.items():
        extra = ', '.join(info['extra_symmetries']) if info['extra_symmetries'] else '-'
        print(f"  {topo:<10} {info['n_bonds']:>5} {info['n_stationary']:>10} "
              f"{info['n_oscillatory']:>11} {info['n_distinct']:>8} "
              f"{info['max_degeneracy']:>7} {extra:>20}")

    # Also N=5 comparison
    print("\n--- TOPOLOGY COMPARISON (N=5, uniform gamma=0.1) ---")
    topo_results_5 = topology_comparison(5)
    print(f"  {'Topo':<10} {'bonds':>5} {'stationary':>10} "
          f"{'distinct':>8} {'max_deg':>7}")
    print(f"  {'-' * 45}")
    for topo, info in topo_results_5.items():
        print(f"  {topo:<10} {info['n_bonds']:>5} {info['n_stationary']:>10} "
              f"{info['n_distinct']:>8} {info['max_degeneracy']:>7}")

    # ============================================================
    # Save results
    # ============================================================
    output = dict(
        sector_counts={str(N): {
            'n_basic': d['n_basic_sectors'],
            'n_with_flip': d['n_sf_blocks'],
        } for N, d in all_sector_data.items()},
        attractors_N5=[{k: v for k, v in a.items()}
                       for a in attractors],
        degeneracy_uniform=dict(
            n_distinct=deg_info['n_distinct'],
            max_mult=deg_info['max_multiplicity'],
        ),
        degeneracy_sacrifice=dict(
            n_distinct=deg_sac['n_distinct'],
            max_mult=deg_sac['max_multiplicity'],
        ),
        reachability={name: {
            'n_sectors': info['n_sectors'],
            'sectors': {f"({k[0]},{k[1]},{k[2]})": round(v, 6)
                        for k, v in info['sectors'].items()}
        } for name, info in reach.items()},
        topology_N4={t: {k: v for k, v in i.items() if k != 'extra_symmetries'}
                     for t, i in topo_results.items()},
        topology_N5={t: {k: v for k, v in i.items() if k != 'extra_symmetries'}
                     for t, i in topo_results_5.items()},
    )

    with open(os.path.join(OUT_DIR, 'census_results.json'), 'w',
              encoding='utf-8') as f:
        json.dump(output, f, indent=2, ensure_ascii=False, default=str)

    print(f"\n  Results saved to {OUT_DIR}/census_results.json")
    print(f"\nDone.")
