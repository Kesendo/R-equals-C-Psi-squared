"""
IBM April 2026 Pre-Registration: Quantitative Predictions
=========================================================
Computes predictions BEFORE the hardware experiment.

Priority 1: Absorption Theorem Verification (N=5)
Priority 3: Quantitative Heartbeat / Dwell Prefactor (N=3)

Authors: Thomas Wicht, Claude (Opus 4.6)
Date: April 9, 2026
"""

import numpy as np
from scipy import linalg
from itertools import product as iter_product
import os

# === Pauli infrastructure ===
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}

def kron_list(mats):
    result = mats[0]
    for m in mats[1:]:
        result = np.kron(result, m)
    return result

def pauli_basis(n):
    labels = ['I', 'X', 'Y', 'Z']
    strings = []
    matrices = []
    for combo in iter_product(labels, repeat=n):
        strings.append(''.join(combo))
        matrices.append(kron_list([PAULIS[c] for c in combo]))
    return strings, matrices

def n_xy_string(s):
    return sum(1 for c in s if c in ('X', 'Y'))

def n_xy_weighted(s, gammas):
    return sum(gammas[k] for k, c in enumerate(s) if c in ('X', 'Y'))

# === Hamiltonian and Liouvillian ===
def heisenberg_H(n, J=1.0):
    d = 2**n
    H = np.zeros((d, d), dtype=complex)
    for i in range(n - 1):
        for pauli in [X, Y, Z]:
            ops = [I2] * n
            ops[i] = pauli
            ops[i + 1] = pauli
            H += J * kron_list(ops)
    return H

def build_liouvillian(H, gammas):
    d = H.shape[0]
    n = int(np.log2(d))
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(n):
        ops = [I2] * n
        ops[k] = Z
        Zk = kron_list(ops)
        gk = gammas[k]
        Lk = np.sqrt(gk) * Zk
        LdL = Lk.conj().T @ Lk
        L += np.kron(Lk.conj(), Lk)
        L -= 0.5 * (np.kron(LdL.T, Id) + np.kron(Id, LdL))
    return L

# === Initial states ===
def bell_plus_center(n):
    c1 = (n - 1) // 2
    bell = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    rho_parts = []
    i = 0
    while i < n:
        if i == c1:
            rho_parts.append(np.outer(bell, bell.conj()))
            i += 2
        else:
            rho_parts.append(np.outer(plus, plus.conj()))
            i += 1
    rho = rho_parts[0]
    for p in rho_parts[1:]:
        rho = np.kron(rho, p)
    return rho


def w_full(n):
    """Fully delocalized W_n state on all n qubits.

    |W_n> = (|10...0> + |01...0> + ... + |00...1>) / sqrt(n)
    """
    d = 2**n
    psi = np.zeros(d, dtype=complex)
    for k in range(n):
        psi[1 << (n - 1 - k)] = 1.0 / np.sqrt(n)
    return np.outer(psi, psi.conj())


def w_center3_plus(n):
    """W3 on the three center qubits, |+> on the rest.

    Direct W-analog of bell_plus_center: for n=5, sites {1,2,3} host W3,
    sites {0,4} host |+>. Same spatial geometry as Bell+|+>, different
    initial-state Pauli content.
    """
    if n < 3:
        raise ValueError("w_center3_plus requires n >= 3")
    w3 = np.zeros(8, dtype=complex)
    # |100> = index 4, |010> = index 2, |001> = index 1
    w3[4] = w3[2] = w3[1] = 1.0 / np.sqrt(3)
    rho_w3 = np.outer(w3, w3.conj())
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    rho_plus = np.outer(plus, plus.conj())

    c = (n - 1) // 2  # center index
    # W3 block occupies sites c-1, c, c+1
    rho_parts = []
    i = 0
    while i < n:
        if i == c - 1:
            rho_parts.append(rho_w3)
            i += 3
        else:
            rho_parts.append(rho_plus)
            i += 1
    rho = rho_parts[0]
    for p in rho_parts[1:]:
        rho = np.kron(rho, p)
    return rho


# Ordered list of initial states to test in priority1.
# Each entry: (display_name, constructor(n) -> rho).
INITIAL_STATES = [
    ("Bell+|+⟩",     bell_plus_center),
    ("W3 center|+⟩", w_center3_plus),
    ("W5 full",      w_full),
]

# =============================================================================
# Priority 1: Absorption Theorem (N=5 eigenvalue analysis only)
# =============================================================================
def priority1():
    print("=" * 70)
    print("PRIORITY 1: ABSORPTION THEOREM VERIFICATION")
    print("=" * 70)

    N = 5
    J = 1.0

    # IBM calibration: T2 for chain [80, 8, 79, 53, 85]
    T2_us = np.array([5.22, 122.70, 243.85, 169.97, 237.57])
    gamma_phys = 1.0 / (2.0 * T2_us)

    # Scale to dimensionless units with γ_min = 0.05
    gamma_min = 0.05
    gamma_sacrifice = gamma_phys / gamma_phys.min() * gamma_min
    gamma_uniform = np.ones(N) * np.mean(gamma_sacrifice)

    print(f"\nChain [80, 8, 79, 53, 85]:")
    print(f"  T2 (μs):       {T2_us}")
    print(f"  γ sacrifice:   {np.round(gamma_sacrifice, 5)}")
    print(f"  γ uniform:     {np.round(gamma_uniform, 5)}")
    print(f"  Sacrifice ratio Q80/Q79: {gamma_sacrifice[0]/gamma_sacrifice[2]:.1f}x")
    print(f"  Σγ sacrifice:  {np.sum(gamma_sacrifice):.5f}")
    print(f"  Σγ uniform:    {np.sum(gamma_uniform):.5f}")

    H = heisenberg_H(N, J)
    strings, P_mats = pauli_basis(N)

    results = {}

    for label, gammas in [("sacrifice", gamma_sacrifice), ("uniform", gamma_uniform)]:
        print(f"\n--- {label.upper()} ---")

        L = build_liouvillian(H, gammas)

        # Full eigendecomposition (eigenvectors needed for Pauli analysis)
        eigvals_full, R = linalg.eig(L)
        order = np.argsort(-eigvals_full.real)
        eigvals_full = eigvals_full[order]
        R = R[:, order]

        # Pauli decomposition for ALL modes (sandbox cap removed, full 4^N = 1024 at N=5)
        n_analyze = len(eigvals_full)
        P_flat = np.array([P.flatten() for P in P_mats])  # 4^N x d²
        d2 = 2**(2*N)

        mode_nxy_w = np.zeros(n_analyze)
        mode_nxy_c = np.zeros(n_analyze)

        for k in range(n_analyze):
            v = R[:, k]
            # Pauli decomposition: c_P = Tr(P† · V) / 2^N
            coeffs = P_flat.conj() @ v / (2**N)
            probs = np.abs(coeffs)**2
            p_sum = probs.sum()
            if p_sum < 1e-30:
                continue
            w = probs / p_sum
            mode_nxy_w[k] = sum(wi * n_xy_weighted(s, gammas) for wi, s in zip(w, strings))
            mode_nxy_c[k] = sum(wi * n_xy_string(s) for wi, s in zip(w, strings))

        # Verify Absorption Theorem
        predicted = -2.0 * mode_nxy_w
        actual = eigvals_full[:n_analyze].real
        mask = np.abs(actual) > 1e-10
        if mask.sum() > 0:
            dev = np.abs(predicted[mask] - actual[mask])
            print(f"  Absorption Theorem: max dev = {dev.max():.2e}, mean = {dev.mean():.2e}")

        # Slowest non-stationary modes
        slow_mask = np.abs(eigvals_full[:n_analyze].real) > 1e-10
        slow_idx = np.where(slow_mask)[0][:5]
        print(f"  Slowest 5 non-stationary:")
        for i in slow_idx:
            print(f"    Re(λ) = {eigvals_full[i].real:+.6f}"
                  f"  Im(λ) = {eigvals_full[i].imag:+.6f}"
                  f"  ⟨n_XY⟩ = {mode_nxy_c[i]:.3f}")

        # Initialize per-profile results container
        results[label] = {'sum_gamma': float(np.sum(gammas))}

        # Thresholds for rate measures
        OVERLAP_SLOW_THRESHOLD = 0.001
        OVERLAP_EFF_THRESHOLD = 1e-6

        # Loop over initial states, reusing the eigendecomposition
        for state_name, state_maker in INITIAL_STATES:
            print(f"\n  === Initial state: {state_name} ===")
            rho0 = state_maker(N)
            rho0_vec = rho0.flatten()
            try:
                c0 = linalg.solve(R, rho0_vec)
            except Exception:
                c0 = np.linalg.lstsq(R, rho0_vec, rcond=None)[0]
            overlap = np.abs(c0[:n_analyze])**2

            # Top-5 overlap modes (diagnostic)
            dom = np.argsort(-overlap)[:5]
            print(f"    Top-5 overlap modes:")
            for i in dom:
                if overlap[i] > 1e-6:
                    print(f"      overlap={overlap[i]:.4f}"
                          f"  Re(lambda)={eigvals_full[i].real:+.6f}"
                          f"  <n_XY>={mode_nxy_c[i]:.3f}")

            ev_real = eigvals_full[:n_analyze].real
            nonstat = np.abs(ev_real) > 1e-10
            ov = overlap[:n_analyze]

            # (a) SLOWEST: langsamster Mode mit overlap > threshold
            slow_mask = nonstat & (ov > OVERLAP_SLOW_THRESHOLD)
            if slow_mask.sum() > 0:
                slow_idx_arr = np.where(slow_mask)[0]
                slow_i = slow_idx_arr[np.argmax(ev_real[slow_idx_arr])]
                slowest_rate = float(ev_real[slow_i])
                slowest_overlap = float(ov[slow_i])
            else:
                slowest_rate = slowest_overlap = float('nan')

            # (b) DOMINANT: Mode mit höchstem overlap
            if nonstat.sum() > 0:
                nonstat_idx = np.where(nonstat)[0]
                dom_i = nonstat_idx[np.argmax(ov[nonstat_idx])]
                dominant_rate = float(ev_real[dom_i])
                dominant_overlap = float(ov[dom_i])
                dominant_nxy = float(mode_nxy_c[dom_i])
            else:
                dominant_rate = dominant_overlap = dominant_nxy = float('nan')

            # (c) EFFECTIVE: overlap-gewichtetes |Re(lambda)|
            eff_mask = nonstat & (ov > OVERLAP_EFF_THRESHOLD)
            if eff_mask.sum() > 0:
                w = ov[eff_mask]
                rates = np.abs(ev_real[eff_mask])
                effective_rate = -float(np.sum(w * rates) / np.sum(w))
                effective_overlap_sum = float(np.sum(w))
            else:
                effective_rate = float('nan')
                effective_overlap_sum = 0.0

            results[label][state_name] = {
                'slowest_bell_rate': slowest_rate,
                'slowest_bell_overlap': slowest_overlap,
                'dominant_bell_rate': dominant_rate,
                'dominant_bell_overlap': dominant_overlap,
                'dominant_bell_nxy': dominant_nxy,
                'effective_bell_rate': effective_rate,
                'effective_overlap_sum': effective_overlap_sum,
            }
            print(f"    Slowest:   Re(lambda)={slowest_rate:+.6f}"
                  f"  overlap={slowest_overlap:.4f}")
            print(f"    Dominant:  Re(lambda)={dominant_rate:+.6f}"
                  f"  overlap={dominant_overlap:.4f}"
                  f"  <n_XY>={dominant_nxy:.3f}")
            print(f"    Effective: |Re(lambda)|={-effective_rate:.6f}"
                  f"  (sum overlap={effective_overlap_sum:.4f})")

    # Protection factors per initial state (uniform rate / sacrifice rate)
    pf_by_state = {}
    if len(results) == 2:
        s = results['sacrifice']
        u = results['uniform']

        print(f"\n{'=' * 70}")
        print("PROTECTION FACTORS (uniform rate / sacrifice rate)")
        print(f"{'=' * 70}")
        for state_name, _ in INITIAL_STATES:
            if state_name not in s or state_name not in u:
                continue
            ss = s[state_name]
            uu = u[state_name]
            pf_slow = uu['slowest_bell_rate'] / ss['slowest_bell_rate']
            pf_dom = uu['dominant_bell_rate'] / ss['dominant_bell_rate']
            pf_eff = uu['effective_bell_rate'] / ss['effective_bell_rate']
            pf_by_state[state_name] = {
                'slow': pf_slow,
                'dom': pf_dom,
                'eff': pf_eff,
            }
            print(f"\n  [{state_name}]")
            print(f"    Slowest-mode:   {pf_slow:.3f}x"
                  f"  (s={ss['slowest_bell_rate']:.4f},"
                  f" u={uu['slowest_bell_rate']:.4f},"
                  f" overlap_s={ss['slowest_bell_overlap']:.4f})")
            print(f"    Dominant-mode:  {pf_dom:.3f}x"
                  f"  (s={ss['dominant_bell_rate']:.4f},"
                  f" u={uu['dominant_bell_rate']:.4f},"
                  f" overlap_s={ss['dominant_bell_overlap']:.4f})")
            print(f"    Effective:      {pf_eff:.3f}x"
                  f"  (overlap-weighted mean |Re(lambda)|)")
        print(f"\n{'=' * 70}")
        return pf_by_state, results
    return None, results

# =============================================================================
# Priority 3: CΨ Heartbeat (N=3)
# =============================================================================
def partial_trace_to_pair(rho, n, i, j):
    rho_r = rho.reshape([2]*n + [2]*n)
    keep = sorted([i, j])
    trace_over = [k for k in range(n) if k not in keep]
    in_idx = list(range(2*n))
    out_idx = [keep[0], keep[1], keep[0]+n, keep[1]+n]
    for k in trace_over:
        in_idx[k + n] = in_idx[k]
    result = np.einsum(rho_r, in_idx, out_idx)
    return result.reshape(4, 4)

def wootters_concurrence(rho):
    sy = np.kron(Y, Y)
    rho_tilde = sy @ rho.conj() @ sy
    sqrt_rho = linalg.sqrtm(rho)
    M = sqrt_rho @ rho_tilde @ sqrt_rho
    eigs = np.sort(np.real(np.sqrt(np.maximum(linalg.eigvals(M).real, 0))))[::-1]
    return max(0, eigs[0] - eigs[1] - eigs[2] - eigs[3])

def priority3():
    print(f"\n{'=' * 70}")
    print("PRIORITY 3: QUANTITATIVE HEARTBEAT (N=3)")
    print("=" * 70)

    N = 3
    J = 1.0
    H = heisenberg_H(N, J)
    pairs = [(0, 1), (1, 2), (0, 2)]

    # States to test
    bell3 = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    rho_bell = np.kron(np.outer(bell3, bell3.conj()), np.outer(plus, plus.conj()))

    d = 2**N
    w_psi = np.zeros(d, dtype=complex)
    for k in range(N):
        w_psi[1 << (N - 1 - k)] = 1.0 / np.sqrt(N)
    rho_w = np.outer(w_psi, w_psi.conj())

    ghz_psi = np.zeros(d, dtype=complex)
    ghz_psi[0] = ghz_psi[-1] = 1.0 / np.sqrt(2)
    rho_ghz = np.outer(ghz_psi, ghz_psi.conj())

    gamma = 0.05
    gammas = np.ones(N) * gamma

    L = build_liouvillian(H, gammas)
    eigvals, R = linalg.eig(L)
    order = np.argsort(-eigvals.real)
    eigvals = eigvals[order]
    R = R[:, order]
    R_inv = linalg.inv(R)

    times = np.linspace(0, 60, 601)
    state_results = {}

    for state_name, rho0 in [("Bell+|+⟩", rho_bell), ("W₃", rho_w), ("GHZ₃", rho_ghz)]:
        print(f"\n--- {state_name} at γ = {gamma} ---")

        c0 = R_inv @ rho0.flatten()

        cpsi_traj = []
        for t in times:
            rho_vec = R @ (c0 * np.exp(eigvals * t))
            rho_t = rho_vec.reshape(d, d)
            rho_t = (rho_t + rho_t.conj().T) / 2

            concs = []
            for (i, j) in pairs:
                rho2 = partial_trace_to_pair(rho_t, N, i, j)
                concs.append(wootters_concurrence(rho2))
            cpsi_traj.append(np.mean(concs))

        cpsi = np.array(cpsi_traj)
        above = cpsi > 0.25
        crossings = np.sum(np.diff(above.astype(int)) != 0)

        eps = 0.01
        near = np.abs(cpsi - 0.25) < eps
        dwell = np.sum(near) * (times[1] - times[0])

        cross_times = []
        cross_idx = np.where(np.diff(above.astype(int)) != 0)[0]
        for ci in cross_idx[:10]:
            cross_times.append(times[ci])

        print(f"  CΨ(0) = {cpsi[0]:.6f}")
        print(f"  CΨ_max = {cpsi.max():.6f}")
        print(f"  Crossings of 1/4: {crossings}")
        print(f"  Dwell time near 1/4 (±0.01): {dwell:.2f}")
        if cross_times:
            print(f"  First 10 crossing times: {[f'{t:.2f}' for t in cross_times]}")
        if cpsi.max() < 0.25:
            print(f"  ** NEVER reaches CΨ = 1/4 **")

        state_results[state_name] = {
            'cpsi_0': float(cpsi[0]),
            'cpsi_max': float(cpsi.max()),
            'crossings': int(crossings),
            'dwell_near_quarter': float(dwell),
            'crossing_times': [float(t) for t in cross_times],
            'reaches_quarter': bool(cpsi.max() >= 0.25),
            'gamma': float(gamma),
        }

    # Analytical prefactor comparison
    print(f"\n--- Analytical Dwell Prefactors ---")
    print(f"  Bell+ predicted: 1.080088 (from DWELL_PREFACTOR_GENERALIZED)")
    print(f"  W₃ predicted:    0.876832")
    print(f"  Bell+ lingers {1.080088/0.876832:.2f}x longer at the boundary")
    print(f"  GHZ₃: structurally excluded (born below 1/4)")

    return state_results

# =============================================================================
if __name__ == "__main__":
    out = os.path.join(os.path.dirname(__file__), "results", "ibm_april_predictions")
    os.makedirs(out, exist_ok=True)

    print("IBM April 2026 Pre-Registration")
    print(f"Date: {np.datetime64('today')}\n")

    pf, p1 = priority1()
    p3 = priority3()

    # Save summary
    with open(os.path.join(out, "predictions.txt"), 'w', encoding='utf-8') as f:
        f.write("IBM April 2026 Pre-Registered Predictions\n")
        f.write(f"Computed: {np.datetime64('today')}\n")
        f.write("\n")
        f.write("METHODOLOGICAL NOTE\n")
        f.write("-------------------\n")
        f.write("Priority 1 uses a SHAPE COMPARISON at matched Sum(gamma).\n")
        f.write("The sacrifice and uniform profiles have IDENTICAL total noise,\n")
        f.write("only the spatial distribution differs. This isolates the shape\n")
        f.write("effect from total-noise differences. On IBM hardware, selective\n")
        f.write("DD vs uniform DD do NOT have matched Sum(gamma); the numbers\n")
        f.write("below are an upper bound on the shape component of any measured\n")
        f.write("hardware advantage, not a direct hardware prediction.\n")
        f.write("\n")
        f.write("=" * 70 + "\n")
        f.write("PRIORITY 1: Absorption Theorem (N=5, chain [80, 8, 79, 53, 85])\n")
        f.write("=" * 70 + "\n")
        if pf:
            f.write("\nProtection factors per initial state"
                    " (uniform rate / sacrifice rate):\n\n")
            for state_name, d in pf.items():
                f.write(f"  [{state_name}]\n")
                f.write(f"    Slowest-mode:   {d['slow']:.3f}x\n")
                f.write(f"    Dominant-mode:  {d['dom']:.3f}x\n")
                f.write(f"    Effective:      {d['eff']:.3f}x"
                        f"  (overlap-weighted mean |Re(lambda)|)\n")
            f.write("\nProfile details:\n")
            for label in ('sacrifice', 'uniform'):
                profile = p1[label]
                f.write(f"\n  [{label}]  Sum(gamma) = {profile['sum_gamma']:.5f}\n")
                for state_name, _ in INITIAL_STATES:
                    if state_name not in profile:
                        continue
                    v = profile[state_name]
                    f.write(f"    --- {state_name} ---\n")
                    f.write(f"      slowest rate   = {v['slowest_bell_rate']:.6f}"
                            f"  (overlap {v['slowest_bell_overlap']:.4f})\n")
                    f.write(f"      dominant rate  = {v['dominant_bell_rate']:.6f}"
                            f"  (overlap {v['dominant_bell_overlap']:.4f},"
                            f" <n_XY> {v['dominant_bell_nxy']:.3f})\n")
                    f.write(f"      effective rate = {v['effective_bell_rate']:.6f}"
                            f"  (sum overlap {v['effective_overlap_sum']:.4f})\n")

        f.write("\n")
        f.write("=" * 70 + "\n")
        f.write("PRIORITY 3: Heartbeat (N=3)\n")
        f.write("=" * 70 + "\n")
        for state, d in p3.items():
            f.write(f"\n  {state}  (gamma = {d['gamma']})\n")
            f.write(f"    CPsi(0)           = {d['cpsi_0']:.6f}\n")
            f.write(f"    CPsi_max          = {d['cpsi_max']:.6f}\n")
            f.write(f"    reaches 1/4       = {d['reaches_quarter']}\n")
            f.write(f"    crossings of 1/4  = {d['crossings']}\n")
            f.write(f"    dwell (+-0.01)    = {d['dwell_near_quarter']:.3f}\n")
            if d['crossing_times']:
                ct_str = ', '.join(f"{t:.3f}" for t in d['crossing_times'])
                f.write(f"    crossing times    = [{ct_str}]\n")
            else:
                f.write("    crossing times    = []\n")
        f.write("\nAnalytical dwell prefactors (from DWELL_PREFACTOR_GENERALIZED):\n")
        f.write("  Bell+ = 1.080088\n")
        f.write("  W3    = 0.876832\n")
        f.write("  ratio = 1.232x  (Bell+ lingers longer at the boundary)\n")

    print(f"\nSaved to {out}/predictions.txt")
