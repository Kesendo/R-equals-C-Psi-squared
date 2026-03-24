#!/usr/bin/env python3
"""
Resonant Return V2: Multi-Mode Optimization + Spatially Structured Pulsing
==========================================================================
Test A: Multi-mode optimization (N=5, N=7)
  A1: Individual SVD modes (2, 3, 4)
  A2: Pairwise combinations of modes 2+3
  A3: Numerical optimization (N=5 only, Nelder-Mead)
Test B: Spatially structured pulsing (N=5)
  B1: Mode 2 spatial x resonant frequency (omega_dom)
  B2: Mode 2 spatial x slow frequency (omega_dom/10)
  B3: Mode 2 spatial x double frequency (2*omega_dom)

Script:  simulations/resonant_return_v2.py
Output:  simulations/results/resonant_return_v2.txt
Builds on: resonant_return.py
"""

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize
import os, sys, time as _time

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "resonant_return_v2.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


# ============================================================
# INFRASTRUCTURE (from resonant_return.py)
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


def evolve_expm(L, rho0, t):
    d = int(np.sqrt(L.shape[0]))
    v = expm(L * t) @ rho0.flatten()
    rho = v.reshape(d, d)
    return (rho + rho.conj().T) / 2


def make_init_state_plus(N):
    psi = plus
    for _ in range(N - 1):
        psi = np.kron(psi, plus)
    return np.outer(psi, psi.conj())


def make_init_state_bell(N):
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


def sum_mi_adjacent(rho, N):
    total = 0.0
    for i in range(N - 1):
        total += mutual_info(rho, i, i + 1, N)
    return total


# ============================================================
# RK4 EVOLUTION (density matrix space - fast for large N)
# ============================================================

def evolve_rk4(H, gammas, rho0, N, t_target, dt=0.05, Zk_ops=None):
    """Evolve density matrix to t_target using RK4.

    Works in d x d space instead of d^2 x d^2 Liouville space.
    At N=7: ~5s per evaluation vs ~8 min with expm.
    """
    if Zk_ops is None:
        Zk_ops = [site_op(sz, k, N) for k in range(N)]

    rho = rho0.astype(complex).copy()
    n_steps = max(1, int(round(t_target / dt)))
    h = t_target / n_steps
    g_clamped = [max(0.001, g) for g in gammas]

    def rhs(r):
        dr = -1j * (H @ r - r @ H)
        for k in range(N):
            dr += g_clamped[k] * (Zk_ops[k] @ r @ Zk_ops[k] - r)
        return dr

    for _ in range(n_steps):
        k1 = rhs(rho)
        k2 = rhs(rho + 0.5 * h * k1)
        k3 = rhs(rho + 0.5 * h * k2)
        k4 = rhs(rho + h * k3)
        rho = rho + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        rho = (rho + rho.conj().T) / 2
        tr = np.trace(rho).real
        if tr > 0:
            rho /= tr

    return rho


# ============================================================
# EVALUATION HELPERS
# ============================================================

def eval_sum_mi(gammas, H, rho0, N, t_meas=5.0, method='expm', Zk_ops=None):
    """Compute Sum_MI at t_meas for a given gamma profile."""
    if method == 'expm':
        L = build_L(H, gammas, N)
        rho = evolve_expm(L, rho0, t_meas)
    else:
        rho = evolve_rk4(H, gammas, rho0, N, t_meas, Zk_ops=Zk_ops)
    return sum_mi_adjacent(rho, N)


def extract_features(H, gammas, rho0, N, t_meas=5.0, method='expm', Zk_ops=None):
    """Extract feature vector (purities, cpsi, MI) for response matrix."""
    if method == 'expm':
        L = build_L(H, gammas, N)
        rho = evolve_expm(L, rho0, t_meas)
    else:
        rho = evolve_rk4(H, gammas, rho0, N, t_meas, Zk_ops=Zk_ops)

    features = []
    for q in range(N):
        features.append(purity(ptrace_keep(rho, [q], N)))
    for i in range(N - 1):
        features.append(cpsi(ptrace_keep(rho, [i, i + 1], N)))
    features.append(mutual_info(rho, 0, N - 1, N))
    return np.array(features)


def compute_response_svd(N, gamma_base, t_meas, method, Zk_ops=None):
    """Compute response matrix SVD. Returns H, rho0, Vt, sv."""
    H = build_H(N)
    rho0 = make_init_state_plus(N)
    gammas_base = [gamma_base] * N
    dg = 1e-4

    feat_base = extract_features(H, gammas_base, rho0, N, t_meas, method, Zk_ops)
    n_feat = len(feat_base)
    R = np.zeros((n_feat, N))

    for site in range(N):
        gammas_pert = list(gammas_base)
        gammas_pert[site] += dg
        feat_pert = extract_features(H, gammas_pert, rho0, N, t_meas, method, Zk_ops)
        R[:, site] = (feat_pert - feat_base) / dg

    _, sv, Vt = np.linalg.svd(R, full_matrices=False)
    return H, rho0, Vt, sv


def compute_eps_for_ratio(v, gamma_base, target_ratio):
    """Compute epsilon so max(gamma)/min(gamma) ~ target_ratio."""
    v_max, v_min = np.max(v), np.min(v)
    if v_max - v_min < 1e-10:
        return 0.0

    denom = v_max - target_ratio * v_min
    if denom > 1e-10:
        eps = gamma_base * (target_ratio - 1) / denom
    else:
        # Fallback: scale max perturbation to 0.02
        v_range = max(abs(v_max), abs(v_min))
        eps = 0.02 / v_range if v_range > 0 else 0.0

    # Ensure all gammas > 0.001
    for val in v:
        if gamma_base + eps * val < 0.001:
            if abs(val) > 1e-10:
                eps = min(eps, (gamma_base - 0.001) / abs(val))

    return max(0.0, eps)


def make_profile(gamma_base, eps, v, N):
    """Create clamped gamma profile."""
    return [max(0.001, gamma_base + eps * v[k]) for k in range(N)]


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    t_start = _time.time()

    log("=== RESONANT RETURN V2: MULTI-MODE OPTIMIZATION + SPATIALLY STRUCTURED PULSING ===")
    log(f"Started: {_time.strftime('%Y-%m-%d %H:%M:%S')}")
    log()

    gamma_base = 0.05
    t_meas = 5.0

    # ========================================================
    # TEST A: MULTI-MODE OPTIMIZATION
    # ========================================================
    log("=" * 70)
    log("TEST A: MULTI-MODE OPTIMIZATION")
    log("=" * 70)
    log()

    saved = {}

    for N in [5, 7]:
        method = 'expm' if N <= 5 else 'rk4'
        t0_N = _time.time()

        log(f"{'=' * 50}")
        log(f"N = {N} (method = {method})")
        log(f"{'=' * 50}")
        log()

        # Precompute Z operators for RK4
        Zk_ops = [site_op(sz, k, N) for k in range(N)]

        # --- Compute SVD ---
        log("Computing response matrix SVD...")
        H, rho0, Vt, sv = compute_response_svd(N, gamma_base, t_meas, method, Zk_ops)
        n_modes = min(N, len(sv))

        log(f"  Singular values: {', '.join(f'{s:.4f}' for s in sv[:n_modes])}")
        for m in range(min(4, n_modes)):
            log(f"  Mode {m+1}: [{', '.join(f'{v:.4f}' for v in Vt[m])}]")
        log()

        # --- V-shape baseline ---
        if N == 5:
            vshape = [0.07, 0.06, 0.05, 0.06, 0.07]
        elif N == 7:
            vshape = [0.08, 0.07, 0.06, 0.05, 0.06, 0.07, 0.08]
        else:
            vshape = [gamma_base + 0.02 * abs(2 * i / (N - 1) - 1) for i in range(N)]

        mi_vshape = eval_sum_mi(vshape, H, rho0, N, t_meas, method, Zk_ops)

        # --- Reference: mode 2 profile with original eps ---
        v1_range = np.max(np.abs(Vt[0]))
        eps_orig = 0.02 / v1_range if v1_range > 0 else 0.0
        g_mode2_orig = make_profile(gamma_base, eps_orig, Vt[1], N)
        target_ratio = max(g_mode2_orig) / min(g_mode2_orig) if min(g_mode2_orig) > 0 else 4.0
        mi_mode2 = eval_sum_mi(g_mode2_orig, H, rho0, N, t_meas, method, Zk_ops)

        log(f"  V-shape:  Sum_MI = {mi_vshape:.6f}")
        log(f"  Mode 2:   Sum_MI = {mi_mode2:.6f} ({mi_mode2/mi_vshape:.2f}x vs V-shape)")
        log(f"  Mode 2 profile: [{', '.join(f'{g:.4f}' for g in g_mode2_orig)}]")
        log(f"  Target max/min ratio: {target_ratio:.2f}")
        log()

        # --- RK4 cross-validation (N=5 only) ---
        if N == 5:
            mi_rk4_check = eval_sum_mi(g_mode2_orig, H, rho0, N, t_meas, 'rk4', Zk_ops)
            log(f"  RK4 validation: expm={mi_mode2:.6f}, rk4={mi_rk4_check:.6f}, "
                f"diff={abs(mi_mode2 - mi_rk4_check):.2e}")
            log()

        # -----------------------------------------------
        # A1: Individual Modes
        # -----------------------------------------------
        log("--- A1: Individual Modes ---")
        a1 = {}

        for mode_idx in [1, 2, 3]:  # SVD modes 2, 3, 4
            mode_name = f"Mode {mode_idx + 1}"
            if mode_idx >= n_modes:
                log(f"  {mode_name}: N/A (only {n_modes} modes)")
                continue

            if mode_idx == 1:
                # Mode 2: use pre-computed result for exact baseline match
                gammas_m = g_mode2_orig
                mi_m = mi_mode2
                eps_m = eps_orig
            else:
                # Modes 3, 4: scale to match mode 2 contrast
                v_m = Vt[mode_idx]
                eps_m = compute_eps_for_ratio(v_m, gamma_base, target_ratio)
                gammas_m = make_profile(gamma_base, eps_m, v_m, N)
                mi_m = eval_sum_mi(gammas_m, H, rho0, N, t_meas, method, Zk_ops)

            factor_m = mi_m / mi_vshape if mi_vshape > 0 else float('inf')
            ratio_m = max(gammas_m) / min(gammas_m) if min(gammas_m) > 0 else float('inf')

            a1[mode_idx] = {'mi': mi_m, 'factor': factor_m, 'gammas': gammas_m}

            log(f"  {mode_name}: Sum_MI = {mi_m:.6f} ({factor_m:.2f}x vs V-shape)")
            log(f"    Profile: [{', '.join(f'{g:.4f}' for g in gammas_m)}]")
            log(f"    max/min = {ratio_m:.2f}, sigma = {sv[mode_idx]:.4f}")

        log()

        # -----------------------------------------------
        # A2: Mode 2+3 Combinations
        # -----------------------------------------------
        log("--- A2: Mode 2+3 Combinations ---")
        a2 = {}

        if n_modes >= 3:
            weights = [(1.00, 0.00), (0.75, 0.25), (0.50, 0.50),
                       (0.25, 0.75), (0.00, 1.00)]

            best_mi_a2 = 0.0
            best_w = None

            for w2, w3 in weights:
                v_combo = w2 * Vt[1] + w3 * Vt[2]
                norm = np.linalg.norm(v_combo)
                if norm > 1e-10:
                    v_combo = v_combo / norm

                eps_c = compute_eps_for_ratio(v_combo, gamma_base, target_ratio)
                gammas_c = make_profile(gamma_base, eps_c, v_combo, N)
                mi_c = eval_sum_mi(gammas_c, H, rho0, N, t_meas, method, Zk_ops)
                factor_c = mi_c / mi_vshape if mi_vshape > 0 else float('inf')

                a2[(w2, w3)] = {'mi': mi_c, 'factor': factor_c, 'gammas': gammas_c}

                if mi_c > best_mi_a2:
                    best_mi_a2 = mi_c
                    best_w = (w2, w3)

                log(f"  w2={w2:.2f}, w3={w3:.2f}: Sum_MI = {mi_c:.6f} ({factor_c:.2f}x vs V-shape)")

            log()
            if best_w:
                log(f"  Best: w2={best_w[0]:.2f}, w3={best_w[1]:.2f} "
                    f"-> Sum_MI = {best_mi_a2:.6f} ({best_mi_a2/mi_vshape:.2f}x)")
        else:
            log("  N/A (fewer than 3 SVD modes)")

        log()

        elapsed_N = _time.time() - t0_N
        log(f"N={N} completed in {elapsed_N:.1f}s")
        log()

        saved[N] = {
            'a1': a1, 'a2': a2, 'Vt': Vt, 'sv': sv,
            'mi_vshape': mi_vshape, 'mi_mode2': mi_mode2,
            'mode2_gammas': g_mode2_orig, 'eps_orig': eps_orig,
            'H': H, 'rho0': rho0, 'Zk_ops': Zk_ops,
            'target_ratio': target_ratio,
        }

    # ========================================================
    # A3: NUMERICAL OPTIMIZATION (N=5)
    # ========================================================
    log("=" * 50)
    log("A3: NUMERICAL OPTIMIZATION (N=5)")
    log("=" * 50)
    log()

    N = 5
    s5 = saved[5]
    H5, rho0_5 = s5['H'], s5['rho0']
    Vt5, sv5 = s5['Vt'], s5['sv']
    mi_mode2_5 = s5['mi_mode2']
    g_mode2_5 = s5['mode2_gammas']

    log(f"  Starting point (mode 2): Sum_MI = {mi_mode2_5:.6f}")
    log(f"  Profile: [{', '.join(f'{g:.4f}' for g in g_mode2_5)}]")
    log(f"  Mean: {np.mean(g_mode2_5):.6f}")
    log()

    # Parametrize: N-1 free deviations, Nth maintains mean = gamma_base
    delta_start = np.array([g_mode2_5[k] - gamma_base for k in range(N)])
    delta_start -= np.mean(delta_start)  # center to zero-mean
    x0 = delta_start[:N - 1]

    eval_count = [0]
    best_found = [0.0]

    def objective(x):
        eval_count[0] += 1
        delta_last = -np.sum(x)
        gammas = [gamma_base + x[k] for k in range(N - 1)] + [gamma_base + delta_last]

        # Penalty for gammas below threshold
        penalty = 0.0
        for g in gammas:
            if g < 0.001:
                penalty += 10000 * (0.001 - g) ** 2
        if penalty > 0:
            return penalty

        mi = eval_sum_mi(gammas, H5, rho0_5, N, t_meas, 'expm')

        if mi > best_found[0]:
            best_found[0] = mi
        if eval_count[0] % 50 == 0:
            log(f"    [{eval_count[0]} evals] best = {best_found[0]:.6f}")

        return -mi

    log("  Running Nelder-Mead optimizer (max 500 evals)...")
    t0_opt = _time.time()

    result = minimize(objective, x0, method='Nelder-Mead',
                      options={'maxfev': 500, 'xatol': 1e-5, 'fatol': 1e-8,
                               'adaptive': True})

    opt_time = _time.time() - t0_opt

    # Reconstruct optimal profile
    x_opt = result.x
    delta_last = -np.sum(x_opt)
    opt_gammas = [gamma_base + x_opt[k] for k in range(N - 1)] + [gamma_base + delta_last]
    opt_mi = eval_sum_mi(opt_gammas, H5, rho0_5, N, t_meas, 'expm')

    log()
    log(f"  Optimizer finished: {eval_count[0]} evaluations in {opt_time:.1f}s")
    log(f"  Converged: {result.success} ({result.message})")
    log(f"  Optimized profile: [{', '.join(f'{g:.4f}' for g in opt_gammas)}]")
    log(f"  Mean: {np.mean(opt_gammas):.6f}")
    log(f"  Optimized Sum_MI: {opt_mi:.6f}")
    log(f"  Improvement over mode 2: {opt_mi / mi_mode2_5:.2f}x")
    log(f"  SVD efficiency: {mi_mode2_5 / opt_mi * 100:.1f}%")
    log()

    # Decompose into SVD modes
    delta_opt = np.array([opt_gammas[k] - gamma_base for k in range(N)])
    svd_weights = Vt5 @ delta_opt

    log("  SVD decomposition of optimal profile:")
    for m in range(N):
        pct = abs(svd_weights[m]) / np.linalg.norm(delta_opt) * 100 if np.linalg.norm(delta_opt) > 0 else 0
        log(f"    Mode {m+1}: weight = {svd_weights[m]:+.6f} ({pct:.1f}% of ||delta||, sigma = {sv5[m]:.4f})")

    residual = delta_opt - Vt5.T @ svd_weights
    log(f"  Residual: {np.linalg.norm(residual):.2e} (should be ~0)")
    log()

    # ========================================================
    # TEST B: SPATIALLY STRUCTURED PULSING (N=5)
    # ========================================================
    log("=" * 70)
    log("TEST B: SPATIALLY STRUCTURED PULSING (N=5)")
    log("=" * 70)
    log()

    N = 5
    H_B = saved[5]['H']
    Vt_B = saved[5]['Vt']
    eps_B = saved[5]['eps_orig']
    v_mode2 = Vt_B[1]
    Zk_ops_B = saved[5]['Zk_ops']

    # Bell initial state for pulsing test
    rho0_bell = make_init_state_bell(N)

    # Find dominant frequency from Liouvillian eigenvalues
    L_uniform = build_L(H_B, [gamma_base] * N, N)
    evals = np.linalg.eigvals(L_uniform)
    osc_freqs = np.abs(evals.imag)
    osc_freqs = osc_freqs[osc_freqs > 0.01]
    unique_freqs = np.unique(np.round(osc_freqs, 2))
    freq_counts = [(f, np.sum(np.abs(osc_freqs - f) < 0.05)) for f in unique_freqs]
    freq_counts.sort(key=lambda x: -x[1])
    omega_dom = freq_counts[0][0]

    log(f"  Dominant frequency: omega_dom = {omega_dom:.4f}")
    log(f"  Mode 2 pattern: [{', '.join(f'{v:.2f}' for v in v_mode2)}]")
    log(f"  Mode 2 eps: {eps_B:.6f}")

    static_gammas = make_profile(gamma_base, eps_B, v_mode2, N)
    log(f"  Static mode 2 profile: [{', '.join(f'{g:.4f}' for g in static_gammas)}]")
    log()

    # Pulsing scenarios
    dt = 0.05
    t_max = 20.0
    n_steps = int(round(t_max / dt))

    def gamma_static(t):
        return static_gammas

    def gamma_b1(t):
        s = np.sin(omega_dom * t)
        return [max(0.001, gamma_base + eps_B * v_mode2[k] * s) for k in range(N)]

    def gamma_b2(t):
        s = np.sin(omega_dom / 10 * t)
        return [max(0.001, gamma_base + eps_B * v_mode2[k] * s) for k in range(N)]

    def gamma_b3(t):
        s = np.sin(2 * omega_dom * t)
        return [max(0.001, gamma_base + eps_B * v_mode2[k] * s) for k in range(N)]

    scenarios = [
        ("Static mode 2", gamma_static),
        ("B1: Mode 2 x resonant (w_dom)", gamma_b1),
        ("B2: Mode 2 x slow (w_dom/10)", gamma_b2),
        ("B3: Mode 2 x 2*w_dom", gamma_b3),
    ]

    all_b = {}

    for name, gamma_fn in scenarios:
        t0_sc = _time.time()
        rho = rho0_bell.copy()
        times = []
        mi_vals = []

        for step in range(n_steps + 1):
            t = step * dt

            if step % 20 == 0:  # sample every 1.0 time units
                times.append(t)
                mi_vals.append(sum_mi_adjacent(rho, N))

            if step < n_steps:
                gammas_t = gamma_fn(t)
                g_cl = [max(0.001, g) for g in gammas_t]

                def rhs(r, _gc=g_cl):
                    dr = -1j * (H_B @ r - r @ H_B)
                    for k in range(N):
                        dr += _gc[k] * (Zk_ops_B[k] @ r @ Zk_ops_B[k] - r)
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

        elapsed_sc = _time.time() - t0_sc
        all_b[name] = (times, mi_vals)

        log(f"  --- {name} ({elapsed_sc:.1f}s) ---")
        log(f"    {'t':>6}  {'Sum-MI':>10}")
        log("    " + "-" * 20)
        for i in range(0, len(times), max(1, len(times) // 10)):
            log(f"    {times[i]:6.1f}  {mi_vals[i]:10.6f}")
        peak_mi = max(mi_vals)
        peak_t = times[int(np.argmax(mi_vals))]
        log(f"    Peak Sum-MI: {peak_mi:.6f} at t={peak_t:.1f}")
        log()

    # Comparison
    log("  COMPARISON:")
    peak_static = max(all_b["Static mode 2"][1])
    for name, (_, mi_v) in all_b.items():
        pk = max(mi_v)
        if "Static" in name:
            log(f"    {name:35s}: peak Sum_MI = {pk:.6f}")
        else:
            ratio = pk / peak_static if peak_static > 0 else float('inf')
            log(f"    {name:35s}: peak Sum_MI = {pk:.6f} ({ratio:.2f}x vs static)")

    log()

    # Prediction check
    peak_b1 = max(all_b["B1: Mode 2 x resonant (w_dom)"][1])
    peak_b2 = max(all_b["B2: Mode 2 x slow (w_dom/10)"][1])
    peak_b3 = max(all_b["B3: Mode 2 x 2*w_dom"][1])

    log("  PREDICTION CHECK:")

    def check(pred, cond):
        tag = "CONFIRMED" if cond else "FALSIFIED"
        log(f"    [{tag}] {pred}")

    check(f"B1 > static (temporal resonance helps): {peak_b1:.6f} vs {peak_static:.6f}",
          peak_b1 > peak_static * 1.01)
    check(f"B1 > B2 (resonant > slow): {peak_b1:.6f} vs {peak_b2:.6f}",
          peak_b1 > peak_b2 * 1.01)
    check(f"B1 > B3 (fundamental > 2nd harmonic): {peak_b1:.6f} vs {peak_b3:.6f}",
          peak_b1 > peak_b3 * 1.01)
    check(f"B2 ~ static (slow pulsing averages out): |{peak_b2:.6f} - {peak_static:.6f}| / {peak_static:.6f} = {abs(peak_b2-peak_static)/peak_static:.4f}",
          abs(peak_b2 - peak_static) / peak_static < 0.05)

    log()

    # ========================================================
    # SUMMARY
    # ========================================================
    log("=" * 70)
    log("SUMMARY")
    log("=" * 70)
    log()

    log("  TEST A - Multi-Mode Optimization:")
    for N_s in [5, 7]:
        s = saved[N_s]
        m = 'expm' if N_s <= 5 else 'rk4'
        log(f"    N={N_s} ({m}):")
        log(f"      V-shape baseline:  {s['mi_vshape']:.6f}")
        for m_idx in sorted(s['a1'].keys()):
            d = s['a1'][m_idx]
            log(f"      Mode {m_idx+1}:           {d['mi']:.6f} ({d['factor']:.2f}x)")
        if s['a2']:
            best_k = max(s['a2'], key=lambda k: s['a2'][k]['mi'])
            bd = s['a2'][best_k]
            log(f"      Best combo (w2={best_k[0]:.2f},w3={best_k[1]:.2f}): "
                f"{bd['mi']:.6f} ({bd['factor']:.2f}x)")
        log()

    log("  TEST A3 - Numerical Optimization (N=5):")
    log(f"    Mode 2 Sum_MI:     {mi_mode2_5:.6f}")
    log(f"    Optimized Sum_MI:  {opt_mi:.6f}")
    log(f"    Improvement:       {opt_mi/mi_mode2_5:.2f}x")
    log(f"    SVD efficiency:    {mi_mode2_5/opt_mi*100:.1f}%")
    log()

    log("  TEST B - Spatially Structured Pulsing (N=5):")
    log(f"    Static mode 2:     {peak_static:.6f}")
    log(f"    B1 (resonant):     {peak_b1:.6f} ({peak_b1/peak_static:.2f}x)")
    log(f"    B2 (slow):         {peak_b2:.6f} ({peak_b2/peak_static:.2f}x)")
    log(f"    B3 (2x freq):      {peak_b3:.6f} ({peak_b3/peak_static:.2f}x)")
    log()

    # Key question: does N=7 multi-mode recover the N=5 improvement?
    if 7 in saved and 5 in saved:
        best_n7 = max(
            max((d['mi'] for d in saved[7]['a1'].values()), default=0),
            max((d['mi'] for d in saved[7]['a2'].values()), default=0),
        )
        best_n5 = max(
            max((d['mi'] for d in saved[5]['a1'].values()), default=0),
            max((d['mi'] for d in saved[5]['a2'].values()), default=0),
        )
        n7_best_factor = best_n7 / saved[7]['mi_vshape'] if saved[7]['mi_vshape'] > 0 else 0
        n5_factor = saved[5]['mi_mode2'] / saved[5]['mi_vshape'] if saved[5]['mi_vshape'] > 0 else 0
        log(f"  KEY QUESTION: Does N=7 multi-mode recover the N=5 improvement?")
        log(f"    N=5 best: {n5_factor:.2f}x vs V-shape")
        log(f"    N=7 best: {n7_best_factor:.2f}x vs V-shape")
        if n7_best_factor > n5_factor * 0.95:
            log(f"    -> YES: multi-mode combination recovers the lost performance")
        else:
            log(f"    -> NO: multi-mode does not fully recover (gap remains)")
        log()

    total_time = _time.time() - t_start
    log(f"Total runtime: {total_time:.1f}s ({total_time / 60:.1f} min)")
    log(f"Results saved to: {OUT_PATH}")
    _outf.close()
