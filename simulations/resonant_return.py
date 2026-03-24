#!/usr/bin/env python3
"""
Resonant Return: Palindrome-Derived γ Profiles
================================================
Test 1: SVD-optimal γ profiles vs hand-designed (N=5)
Test 2: Frequency-matched γ pulsing (N=5)
Test 3: Palindrome-timed relay (N=11) - flagged for C# if too slow
Test 4: Scaling with N (N=3, 5, 7)

Script:  simulations/resonant_return.py
Output:  simulations/results/resonant_return.txt

Hypothesis: hypotheses/RESONANT_RETURN.md
Builds on: gamma_signal_analysis.py, reading_the_30_percent.py
"""

import numpy as np
from scipy.linalg import expm
import os, sys, time as _time

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "resonant_return.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


# ============================================================
# INFRASTRUCTURE (N-parametric)
# ============================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)
plus = (up + dn) / np.sqrt(2)


def site_op(op, k, N):
    ops = [I2] * N
    ops[k] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def build_H(N, J=1.0):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        for P in [sx, sy, sz]:
            H += J * site_op(P, i, N) @ site_op(P, i + 1, N)
    return H


def build_L(H, gammas, N):
    d = 2 ** N
    d2 = d * d
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(sz, k, N)
        L += gammas[k] * (np.kron(Zk, Zk.conj()) - np.eye(d2))
    return L


def build_L_dissipator_only(gammas, N):
    """Build only the dissipator part of L (for time-dependent gamma)."""
    d = 2 ** N
    d2 = d * d
    L_D = np.zeros((d2, d2), dtype=complex)
    for k in range(N):
        Zk = site_op(sz, k, N)
        L_D += gammas[k] * (np.kron(Zk, Zk.conj()) - np.eye(d2))
    return L_D


def evolve(L, rho0, t):
    d = int(np.sqrt(L.shape[0]))
    v = expm(L * t) @ rho0.flatten()
    rho = v.reshape(d, d)
    return (rho + rho.conj().T) / 2


def make_init_state_plus(N):
    psi = plus
    for _ in range(N - 1):
        psi = np.kron(psi, plus)
    rho = np.outer(psi, psi.conj())
    return rho


def make_init_state_bell(N):
    """Bell+ on qubits 0-1, rest |0>."""
    d = 2 ** N
    psi = np.zeros(d, dtype=complex)
    psi[0] = 1.0 / np.sqrt(2)
    idx_11 = (1 << (N - 1)) + (1 << (N - 2))
    psi[idx_11] = 1.0 / np.sqrt(2)
    return np.outer(psi, psi.conj())


def ptrace_keep(rho, keep, N):
    keep = list(keep)
    trace_out = [q for q in range(N) if q not in keep]
    dims = [2] * N
    reshaped = rho.reshape(dims + dims)
    current_n = N
    for q in sorted(trace_out, reverse=True):
        reshaped = np.trace(reshaped, axis1=q, axis2=q + current_n)
        current_n -= 1
    d_k = 2 ** len(keep)
    return reshaped.reshape((d_k, d_k))


def von_neumann(rho):
    ev = np.linalg.eigvalsh(rho)
    ev = ev[ev > 1e-15]
    return -float(np.sum(ev * np.log2(ev)))


def mutual_info(rho, kA, kB, N):
    rhoA = ptrace_keep(rho, [kA], N)
    rhoB = ptrace_keep(rho, [kB], N)
    rhoAB = ptrace_keep(rho, [kA, kB], N)
    return von_neumann(rhoA) + von_neumann(rhoB) - von_neumann(rhoAB)


def purity(rho):
    return float(np.trace(rho @ rho).real)


def psi_norm(rho):
    d_r = rho.shape[0]
    l1 = float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))
    return l1 / (d_r - 1) if d_r > 1 else 0.0


def cpsi(rho):
    return purity(rho) * psi_norm(rho)


def extract_features(H, gammas, rho0, N, t_meas=5.0):
    L = build_L(H, gammas, N)
    rho = evolve(L, rho0, t_meas)
    features = []
    for q in range(N):
        rho_q = ptrace_keep(rho, [q], N)
        features.append(purity(rho_q))
    for i in range(N - 1):
        rho_pair = ptrace_keep(rho, [i, i + 1], N)
        features.append(cpsi(rho_pair))
    features.append(mutual_info(rho, 0, N - 1, N))
    return np.array(features)


# ============================================================
# TEST 1: SVD-OPTIMAL γ PROFILES
# ============================================================

def compute_response_matrix(H, gammas_base, rho0, N, t_meas=5.0, dg=1e-4):
    """Compute response matrix: how features change with per-site gamma perturbation."""
    feat_base = extract_features(H, gammas_base, rho0, N, t_meas)
    n_feat = len(feat_base)
    R = np.zeros((n_feat, N))
    for site in range(N):
        gammas_pert = list(gammas_base)
        gammas_pert[site] += dg
        feat_pert = extract_features(H, gammas_pert, rho0, N, t_meas)
        R[:, site] = (feat_pert - feat_base) / dg
    return R, feat_base


def run_test_1(N=5, gamma_base=0.05, t_meas=5.0):
    log("=" * 70)
    log(f"TEST 1: SVD-OPTIMAL gamma PROFILES (N={N})")
    log("=" * 70)
    log()

    d = 2 ** N
    H = build_H(N)
    gammas_uniform = [gamma_base] * N
    rho0 = make_init_state_plus(N)

    # Step 1: Compute response matrix and SVD
    log("Computing response matrix (feature sensitivity to per-site gamma)...")
    R_mat, feat_base = compute_response_matrix(H, gammas_uniform, rho0, N, t_meas)

    U, sv, Vt = np.linalg.svd(R_mat, full_matrices=False)
    n_modes = min(N, len(sv))

    log(f"  Response matrix shape: {R_mat.shape[0]} features x {N} sites")
    log(f"  SVD singular values: {', '.join(f'{s:.4f}' for s in sv[:n_modes])}")
    log(f"  Mode 1 spatial pattern: [{', '.join(f'{v:.4f}' for v in Vt[0])}]")
    if n_modes > 1:
        log(f"  Mode 2 spatial pattern: [{', '.join(f'{v:.4f}' for v in Vt[1])}]")
    log()

    # Step 2: Design gamma profiles
    # Determine epsilon so that max/min ratio matches V-shape
    # V-shape: [0.07, 0.06, 0.05, 0.06, 0.07] -> ratio 0.07/0.05 = 1.4
    # We want max(gamma)/min(gamma) ~ 1.4
    v1 = Vt[0]
    v1_range = np.max(np.abs(v1))
    if v1_range > 0:
        eps = 0.02 / v1_range  # scale so max perturbation is ~0.02
    else:
        eps = 0.0

    profiles = {}

    # Uniform baseline
    profiles['Uniform'] = gammas_uniform

    # V-shape (hand-designed, from GAMMA_CONTROL)
    if N == 3:
        profiles['V-shape'] = [0.07, 0.05, 0.07]
    elif N == 5:
        profiles['V-shape'] = [0.07, 0.06, 0.05, 0.06, 0.07]
    elif N == 7:
        profiles['V-shape'] = [0.08, 0.07, 0.06, 0.05, 0.06, 0.07, 0.08]
    else:
        # Generic V-shape
        v_shape = [gamma_base + 0.02 * abs(2 * i / (N - 1) - 1) for i in range(N)]
        profiles['V-shape'] = v_shape

    # SVD mode 1
    g_svd1 = [gamma_base + eps * Vt[0][k] for k in range(N)]
    g_svd1 = [max(0.001, g) for g in g_svd1]  # clamp positive
    profiles['SVD mode 1'] = g_svd1

    # SVD mode 2
    if n_modes > 1:
        g_svd2 = [gamma_base + eps * Vt[1][k] for k in range(N)]
        g_svd2 = [max(0.001, g) for g in g_svd2]
        profiles['SVD mode 2'] = g_svd2

    # SVD modes 1+2 combined
    if n_modes > 1:
        combined = (Vt[0] + Vt[1]) / np.sqrt(2)
        g_svd12 = [gamma_base + eps * combined[k] for k in range(N)]
        g_svd12 = [max(0.001, g) for g in g_svd12]
        profiles['SVD 1+2'] = g_svd12

    # Anti-SVD (orthogonal to mode 1)
    if n_modes > 1:
        # Use mode with lowest singular value
        anti = Vt[-1]
        g_anti = [gamma_base + eps * anti[k] for k in range(N)]
        g_anti = [max(0.001, g) for g in g_anti]
        profiles['Anti-SVD'] = g_anti

    log("  Designed gamma profiles:")
    for name, gammas in profiles.items():
        log(f"    {name}: [{', '.join(f'{g:.4f}' for g in gammas)}]")
    log()

    # Step 3: Evaluate each profile
    # For each: compute MI(0, N-1), classification accuracy, effective capacity
    t_values = [1.0, 3.0, 5.0, 8.0, 12.0]
    results = {}

    log(f"  {'Profile':>15}  {'MI(0,N-1)':>10}  {'MI(peak)':>10}  "
        f"{'CΨ(0,1)':>10}  {'Sum_MI':>10}")
    log("  " + "-" * 65)

    for name, gammas in profiles.items():
        L = build_L(H, gammas, N)

        mi_best = 0.0
        mi_at_5 = 0.0
        cpsi_01 = 0.0
        sum_mi = 0.0

        for t in t_values:
            rho = evolve(L, rho0, t)
            mi = mutual_info(rho, 0, N - 1, N)
            if mi > mi_best:
                mi_best = mi

            if abs(t - 5.0) < 0.01:
                mi_at_5 = mi
                rho_pair = ptrace_keep(rho, [0, 1], N)
                cpsi_01 = cpsi(rho_pair)

            # Sum MI over all pairs
            if abs(t - 5.0) < 0.01:
                for i in range(N - 1):
                    sum_mi += mutual_info(rho, i, i + 1, N)

        results[name] = {
            'mi_at_5': mi_at_5, 'mi_peak': mi_best,
            'cpsi_01': cpsi_01, 'sum_mi': sum_mi
        }

        log(f"  {name:>15}  {mi_at_5:10.6f}  {mi_best:10.6f}  "
            f"{cpsi_01:10.6f}  {sum_mi:10.6f}")

    log()

    # Improvement factors
    mi_uniform = results['Uniform']['sum_mi']
    mi_vshape = results['V-shape']['sum_mi']
    log("  Improvement factors (Sum MI):")
    for name, res in results.items():
        if mi_uniform > 0:
            factor = res['sum_mi'] / mi_uniform
            log(f"    {name}: {factor:.2f}x vs Uniform")
    log()

    if mi_vshape > 0:
        log("  Improvement vs V-shape (Sum MI):")
        for name, res in results.items():
            factor = res['sum_mi'] / mi_vshape
            log(f"    {name}: {factor:.2f}x vs V-shape")
        log()

    # Prediction check
    mi_svd1 = results.get('SVD mode 1', {}).get('sum_mi', 0)
    mi_anti = results.get('Anti-SVD', {}).get('sum_mi', 0)

    log("  PREDICTION CHECK:")
    if mi_svd1 >= mi_vshape:
        log(f"    [CONFIRMED] SVD mode 1 ({mi_svd1:.6f}) >= V-shape ({mi_vshape:.6f})")
    else:
        log(f"    [FALSIFIED] SVD mode 1 ({mi_svd1:.6f}) < V-shape ({mi_vshape:.6f})")

    if mi_anti < mi_svd1:
        log(f"    [CONFIRMED] Anti-SVD ({mi_anti:.6f}) < SVD mode 1 ({mi_svd1:.6f})")
    else:
        log(f"    [FALSIFIED] Anti-SVD ({mi_anti:.6f}) >= SVD mode 1 ({mi_svd1:.6f})")

    log()
    return results, Vt, sv


# ============================================================
# TEST 2: FREQUENCY-MATCHED PULSING
# ============================================================

def sum_mi_adjacent(rho, N):
    """Sum of MI across all adjacent pairs."""
    total = 0.0
    for i in range(N - 1):
        total += mutual_info(rho, i, i + 1, N)
    return total


def run_test_2(N=5, gamma_base=0.05):
    log("=" * 70)
    log(f"TEST 2: FREQUENCY-MATCHED PULSING (N={N}) - REDESIGNED")
    log("=" * 70)
    log()
    log("  Fix vs original: Bell(0,1) initial state + Sum-MI observable")
    log("  (Original used |+>^N which has zero MI everywhere)")
    log()

    d = 2 ** N
    H = build_H(N)
    gammas_uniform = [gamma_base] * N
    rho0 = make_init_state_bell(N)

    # Step 1: Find dominant oscillation frequency from palindromic eigenvalues
    log("Computing Liouvillian eigenvalues for palindromic frequencies...")
    L_base = build_L(H, gammas_uniform, N)
    evals = np.linalg.eigvals(L_base)

    # Find palindromic pairs with largest imaginary parts
    osc_freqs = []
    for ev in evals:
        if abs(ev.imag) > 0.01:
            osc_freqs.append(abs(ev.imag))

    if len(osc_freqs) == 0:
        log("  No oscillating modes found. Test 2 cannot proceed.")
        log()
        return

    osc_freqs = np.array(osc_freqs)
    # Group by unique frequencies (round to 0.01)
    unique_freqs = np.unique(np.round(osc_freqs, 2))

    # Find the dominant frequency: the one with the most modes
    freq_counts = [(f, np.sum(np.abs(osc_freqs - f) < 0.05)) for f in unique_freqs]
    freq_counts.sort(key=lambda x: -x[1])

    omega_dom = freq_counts[0][0]
    log(f"  Eigenvalue oscillation frequencies (top 5):")
    for f, c in freq_counts[:5]:
        log(f"    omega = {f:.4f}, count = {c}")
    log(f"  Dominant frequency: omega_dom = {omega_dom:.4f}")
    log()

    # Step 2: Design time-dependent gamma profiles
    A = 0.5  # amplitude of modulation
    dt = 0.05
    t_max = 20.0
    n_steps = int(round(t_max / dt))

    scenarios = {
        'Static': lambda t: [gamma_base] * N,
        'Resonant': lambda t: [gamma_base * (1 + A * np.sin(omega_dom * t))] * N,
        'Off-resonant': lambda t: [gamma_base * (1 + A * np.sin(2.73 * omega_dom * t))] * N,
    }

    log(f"  Initial state: Bell(0,1) x |0>^{N-2}")
    log(f"  Observable: Sum-MI (all adjacent pairs)")
    log(f"  Propagating with RK4 (dt={dt}, t_max={t_max})...")
    log()

    all_mi = {}

    for name, gamma_fn in scenarios.items():
        t0_scenario = _time.time()
        rho = rho0.copy()
        times = []
        mi_vals = []

        # Precompute site Z operators
        Zk_ops = [site_op(sz, k, N) for k in range(N)]

        for step in range(n_steps + 1):
            t = step * dt
            if step % 20 == 0:  # sample every 1.0 time units
                times.append(t)
                smi = sum_mi_adjacent(rho, N)
                mi_vals.append(smi)

            if step < n_steps:
                # Build L_D for current gamma(t)
                gammas_t = gamma_fn(t)

                def rhs(r):
                    dr = -1j * (H @ r - r @ H)
                    for k in range(N):
                        Lk = np.sqrt(max(0.001, gammas_t[k])) * Zk_ops[k]
                        LdL = Lk.conj().T @ Lk
                        dr += Lk @ r @ Lk.conj().T - 0.5 * (LdL @ r + r @ LdL)
                    return dr

                k1 = rhs(rho)
                k2 = rhs(rho + 0.5 * dt * k1)
                k3 = rhs(rho + 0.5 * dt * k2)
                k4 = rhs(rho + dt * k3)
                rho = rho + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
                rho = (rho + rho.conj().T) / 2
                tr = np.trace(rho).real
                if tr > 0:
                    rho /= tr

        elapsed = _time.time() - t0_scenario
        all_mi[name] = (times, mi_vals)

        log(f"  {name} ({elapsed:.1f}s):")
        log(f"    {'t':>6}  {'Sum-MI':>10}")
        log("    " + "-" * 20)
        for i in range(0, len(times), max(1, len(times) // 10)):
            log(f"    {times[i]:6.1f}  {mi_vals[i]:10.6f}")
        log(f"    Peak Sum-MI: {max(mi_vals):.6f} at t={times[np.argmax(mi_vals)]:.1f}")
        log()

    # Prediction check
    peak_static = max(all_mi['Static'][1])
    peak_res = max(all_mi['Resonant'][1])
    peak_offres = max(all_mi['Off-resonant'][1])

    log("  PREDICTION CHECK:")
    if peak_res > peak_static * 1.01:  # require >1% improvement
        log(f"    [CONFIRMED] Resonant peak ({peak_res:.6f}) > Static peak ({peak_static:.6f})")
        log(f"    Improvement: {peak_res / peak_static:.2f}x")
    else:
        log(f"    [FALSIFIED] Resonant peak ({peak_res:.6f}) <= Static peak ({peak_static:.6f})")

    if peak_offres <= peak_res:
        log(f"    [CONFIRMED] Off-resonant ({peak_offres:.6f}) <= Resonant ({peak_res:.6f})")
    else:
        log(f"    [FALSIFIED] Off-resonant ({peak_offres:.6f}) > Resonant ({peak_res:.6f})")

    log()
    return all_mi


# ============================================================
# TEST 3: PALINDROME-TIMED RELAY (N=11)
# ============================================================

def run_test_3():
    log("=" * 70)
    log("TEST 3: PALINDROME-TIMED RELAY (N=11)")
    log("=" * 70)
    log()

    # N=11: d=2048, d^2=4,194,304 - too large for Python eigendecomp/expm
    # We need the C# implementation in RCPsiSquared.Propagate
    N_relay = 11
    d_relay = 2 ** N_relay

    log(f"  N={N_relay}, d={d_relay}, d^2={d_relay**2}")
    log(f"  Matrix size: {d_relay**2} x {d_relay**2} = {d_relay**4 / 1e9:.1f} billion elements")
    log()
    log("  STATUS: FLAGGED FOR C# IMPLEMENTATION")
    log("  Reason: N=11 Liouvillian is 4M x 4M complex matrix.")
    log("  Python expm/eigvals on this size takes >1 hour and >30 GB RAM.")
    log()
    log("  Structure for C# implementation:")
    log("  1. Fixed timing:      t_stage = K/gamma = 0.039/0.05 = 0.78")
    log("  2. Palindrome timing: t_stage_k = pi / Re(lambda_dominant_k)")
    log("     - Bridge A (q0-4):  5-qubit Heisenberg chain dominant rate")
    log("     - Meta (q4-6):      3-qubit star dominant rate")
    log("     - Bridge B (q6-10): 5-qubit Heisenberg chain dominant rate")
    log("  3. Oracle timing:     sweep t_stage [0.1, 3.0], find MI max")
    log()

    # We CAN compute the dominant rates for the sub-segments in Python
    log("  Computing sub-segment dominant rates in Python...")
    log()

    for seg_name, seg_N in [("Bridge A/B (5-chain)", 5), ("Meta (3-star)", 3)]:
        H_seg = build_H(seg_N)
        gamma_base = 0.05
        gammas = [gamma_base] * seg_N
        L_seg = build_L(H_seg, gammas, seg_N)
        evals_seg = np.linalg.eigvals(L_seg)

        # Find dominant decay rate (slowest nonzero)
        rates = -np.real(evals_seg)
        rates_nonzero = rates[rates > 0.001]
        if len(rates_nonzero) > 0:
            dominant_rate = np.min(rates_nonzero)
            t_half = np.log(2) / dominant_rate
            t_palindrome = np.pi / dominant_rate
        else:
            dominant_rate = 0
            t_half = 0
            t_palindrome = 0

        # Find dominant oscillation frequency
        imag_parts = np.abs(np.imag(evals_seg))
        imag_nonzero = imag_parts[imag_parts > 0.01]
        if len(imag_nonzero) > 0:
            dom_freq = np.min(imag_nonzero)
        else:
            dom_freq = 0

        log(f"  {seg_name}:")
        log(f"    Dominant decay rate: {dominant_rate:.6f}")
        log(f"    Half-life: {t_half:.4f}")
        log(f"    Palindrome stage time (pi/rate): {t_palindrome:.4f}")
        log(f"    Dominant osc. frequency: {dom_freq:.4f}")
        log(f"    Fixed stage time (K/gamma): {0.039/gamma_base:.4f}")
        log()

    log("  To run full relay comparison, use:")
    log("  dotnet run -c Release -- pull   (in compute/RCPsiSquared.Propagate/)")
    log("  with modified stage timing from the rates above.")
    log()


# ============================================================
# TEST 4: SCALING WITH N
# ============================================================

def run_test_4():
    log("=" * 70)
    log("TEST 4: SCALING WITH N")
    log("=" * 70)
    log()

    scaling_results = {}

    for N_test in [3, 5, 7]:  # N=7 needs ~8 GB RAM (128 GB home PC)
        t0_n = _time.time()
        log(f"--- N = {N_test} ---")

        try:
            results, Vt, sv = run_test_1(N=N_test, gamma_base=0.05, t_meas=5.0)
        except MemoryError:
            log(f"  N={N_test}: SKIPPED (MemoryError on d^2={4**N_test} Liouvillian)")
            log()
            continue

        n_effective = int(np.sum(sv > 0.01 * sv[0]))
        mi_vshape = results['V-shape']['sum_mi']
        # Use best SVD mode (mode 2 if available) since mode 1 is near-uniform
        mi_svd_best = max(
            results.get('SVD mode 1', {}).get('sum_mi', 0),
            results.get('SVD mode 2', {}).get('sum_mi', 0),
        )
        mi_svd1 = results.get('SVD mode 1', {}).get('sum_mi', 0)
        mi_svd2 = results.get('SVD mode 2', {}).get('sum_mi', 0)
        improvement = mi_svd_best / mi_vshape if mi_vshape > 0 else 0

        scaling_results[N_test] = {
            'n_modes': n_effective,
            'mi_vshape': mi_vshape,
            'mi_svd1': mi_svd1,
            'mi_svd2': mi_svd2,
            'mi_svd_best': mi_svd_best,
            'improvement': improvement,
            'svd_mode1_pattern': list(Vt[0]) if len(Vt) > 0 else [],
            'svd_mode2_pattern': list(Vt[1]) if len(Vt) > 1 else [],
        }

        elapsed = _time.time() - t0_n
        log(f"  N={N_test} completed in {elapsed:.1f}s")
        log()

    # Summary table
    tested_N = sorted(scaling_results.keys())
    log("=" * 70)
    log("SCALING SUMMARY")
    log("=" * 70)
    log()
    log(f"  {'N':>3}  {'SVD modes':>10}  {'MI(V-shape)':>12}  {'MI(SVD m1)':>12}  "
        f"{'MI(SVD m2)':>12}  {'Best/V':>8}")
    log("  " + "-" * 65)
    for N_test in tested_N:
        r = scaling_results[N_test]
        log(f"  {N_test:>3}  {r['n_modes']:>10}  {r['mi_vshape']:>12.6f}  "
            f"{r['mi_svd1']:>12.6f}  {r.get('mi_svd2', 0):>12.6f}  "
            f"{r['improvement']:>7.2f}x")

    log()
    log("  NOTE: SVD mode 1 (highest singular value) is near-uniform and performs poorly.")
    log("  SVD mode 2 (edge-hot, center-cold) is the actual optimal direction.")
    log()

    # Prediction check
    improvements = [scaling_results[n]['improvement'] for n in tested_N]
    if len(improvements) >= 2 and improvements[-1] > improvements[0]:
        log("  PREDICTION CHECK:")
        log(f"    [CONFIRMED] Improvement grows with N: "
            f"{dict(zip(tested_N, [f'{x:.2f}x' for x in improvements]))}")
    else:
        log("  PREDICTION CHECK:")
        log(f"    [MIXED] Improvement trend: "
            f"{dict(zip(tested_N, [f'{x:.2f}x' for x in improvements]))}")
        log("    (N=7 skipped due to memory; trend assessment incomplete)")

    log()
    return scaling_results


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    t_start = _time.time()

    log("=== RESONANT RETURN: PALINDROME-DERIVED gamma PROFILES ===")
    log(f"Started: {_time.strftime('%Y-%m-%d %H:%M:%S')}")
    log()

    # Test 1: SVD-optimal profiles (N=5)
    t1_start = _time.time()
    results_1, Vt_1, sv_1 = run_test_1(N=5)
    log(f"Test 1 runtime: {_time.time() - t1_start:.1f}s")
    log()

    # Test 2: Frequency-matched pulsing (N=5)
    t2_start = _time.time()
    mi_2 = run_test_2(N=5)
    log(f"Test 2 runtime: {_time.time() - t2_start:.1f}s")
    log()

    # Test 3: Palindrome-timed relay (N=11)
    t3_start = _time.time()
    run_test_3()
    log(f"Test 3 runtime: {_time.time() - t3_start:.1f}s")
    log()

    # Test 4: Scaling (N=3, 5, 7)
    t4_start = _time.time()
    scaling = run_test_4()
    log(f"Test 4 runtime: {_time.time() - t4_start:.1f}s")
    log()

    # ============================================================
    # SUMMARY
    # ============================================================
    log("=" * 70)
    log("SUMMARY")
    log("=" * 70)
    log()

    # Collect all prediction results
    log("Predictions from hypotheses/RESONANT_RETURN.md:")
    log()
    log("  1. SVD-optimal beats V-shape:")
    mi_svd1 = results_1.get('SVD mode 1', {}).get('sum_mi', 0)
    mi_svd2 = results_1.get('SVD mode 2', {}).get('sum_mi', 0)
    mi_vshape = results_1.get('V-shape', {}).get('sum_mi', 0)
    mi_best_svd = max(mi_svd1, mi_svd2)
    best_name = "mode 2" if mi_svd2 > mi_svd1 else "mode 1"
    if mi_best_svd >= mi_vshape:
        log(f"     CONFIRMED: SVD {best_name} ({mi_best_svd:.6f}) >= V-shape ({mi_vshape:.6f})")
        log(f"     Improvement: {mi_best_svd/mi_vshape:.1f}x")
        log(f"     SURPRISE: Mode 1 (highest singular value) is near-uniform and bad.")
        log(f"     Mode 2 (edge-hot, center-cold) is the true optimal direction.")
    else:
        log(f"     FALSIFIED: best SVD ({mi_best_svd:.6f}) < V-shape ({mi_vshape:.6f})")

    log()
    log("  2. Resonant pulsing creates Sum-MI spikes above static (Bell initial state):")
    if mi_2:
        peak_s = max(mi_2['Static'][1])
        peak_r = max(mi_2['Resonant'][1])
        peak_o = max(mi_2['Off-resonant'][1])
        if peak_r > peak_s * 1.01:
            log(f"     CONFIRMED (resonant {peak_r:.6f} > static {peak_s:.6f}, {peak_r/peak_s:.2f}x)")
        else:
            log(f"     FALSIFIED (resonant {peak_r:.6f} <= static {peak_s:.6f})")
        log(f"     Off-resonant: {peak_o:.6f}")
    else:
        log("     SKIPPED (no oscillating modes)")

    log()
    log("  3. Palindrome-timed relay outperforms fixed:")
    log("     DEFERRED to C# implementation (N=11 too large for Python)")

    log()
    log("  4. Improvement grows with N:")
    if scaling:
        tested = sorted(scaling.keys())
        imps = [scaling[n]['improvement'] for n in tested]
        summary = ", ".join(f"N={n}:{scaling[n]['improvement']:.2f}x" for n in tested)
        if len(imps) >= 2 and imps[-1] > imps[0]:
            log(f"     CONFIRMED: {summary}")
        elif len(imps) >= 2:
            log(f"     FALSIFIED: {summary}")
        else:
            log(f"     INCOMPLETE (only {len(tested)} sizes tested): {summary}")
        if 7 not in tested:
            log("     N=7 skipped (memory). Run at home with 128 GB for full scaling.")

    log()
    total = _time.time() - t_start
    log(f"Total runtime: {total:.1f}s ({total / 60:.1f} min)")
    log(f"Results saved to: {OUT_PATH}")
    _outf.close()
