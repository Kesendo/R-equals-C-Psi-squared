#!/usr/bin/env python3
"""Pi-pair closure investigation.

Testing whether palindromic pair rate conservation
(alpha_fast + alpha_slow = 2 * Sigma_gamma = 2 * N * gamma_0 for uniform gamma_0)
connects to the PTF closure law (Sum_i ln(alpha_i) = 0).

See ClaudeTasks/TASK_PI_PAIR_CLOSURE_INVESTIGATION.md for the full briefing.

Phases:
  1. Build XY chain, Liouvillian with uniform Z-dephasing, eigendecompose.
  2. Verify palindromic pairing and Tr(L) formula at N=3.
  3. Apply J-defect, verify palindromic pairing is preserved.
  4. Propagate rho(t), compute per-site purity, fit alpha_i, compute closure.
  5. Idea A: determinant / trace invariants of exp(L t).
  6. Idea B: voice vs memory mode contribution to alpha_i.
  7. Idea D: Liouvillian resolvent, per-site response functions.
  8. Scale to N=5 if N=3 results are consistent.

Rules from task:
  - XY coupling (XX + YY), NOT Heisenberg.
  - Numbers only from script output.
  - UTF-8 stdout on Windows.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d
from scipy.linalg import eig, expm
from scipy.optimize import minimize_scalar

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
GAMMA_0 = 0.05
J_UNIFORM = 1.0
DELTA_J = 0.1         # perturbation on bond (0,1)
DEFECT_BOND = (0, 1)
T_FINAL = 80.0
N_STEPS = 200          # propagation resolution; fit window uses a subset
T_FIT_MAX = 20.0       # fit alpha_i on t in [0, T_FIT_MAX]
ALPHA_BOUNDS = (0.1, 10.0)
TOL_PALINDROME = 1e-10

RESULTS_DIR = Path(__file__).parent / "results" / "pi_pair_closure_investigation"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Pauli operators and chain builders
# ---------------------------------------------------------------------------
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_chain(*ops):
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def site_op(op, site, N):
    """Op at site (big-endian), identity elsewhere."""
    factors = [I2] * N
    factors[site] = op
    return kron_chain(*factors)


def build_H_XY(J_list, N):
    """H = sum_b (J_b / 2) (X_b X_{b+1} + Y_b Y_{b+1})."""
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        Jb = J_list[b]
        H += (Jb / 2.0) * (site_op(X, b, N) @ site_op(X, b + 1, N) +
                           site_op(Y, b, N) @ site_op(Y, b + 1, N))
    return H


def build_liouvillian_matrix(H, gamma_0, N):
    """Vectorised L such that dvec(rho)/dt = L @ vec(rho).

    vec: column-stacking convention, so X rho Y -> (Y^T kron X) vec(rho).
    Dissipator: D[rho] = gamma_0 sum_i (Z_i rho Z_i - rho).
    """
    d = 2**N
    I_d = np.eye(d, dtype=complex)
    # Hamiltonian part: -i (H rho - rho H) -> -i (I kron H - H^T kron I) vec(rho)
    L = -1j * (np.kron(I_d, H) - np.kron(H.T, I_d))
    # Dissipator: for each site i, Z_i rho Z_i - rho
    #   (Z_i^T kron Z_i) - I kron I  [complex conj not needed: Z is real diagonal]
    for i in range(N):
        Zi = site_op(Z, i, N)
        L += gamma_0 * (np.kron(Zi.T, Zi) - np.kron(I_d, I_d))
    return L


def eig_LR(L):
    """Right and left eigendecomposition.
    Returns (eigvals, V_R, V_L) with L V_R = V_R diag(eigvals)
    and V_L^dag V_R = I (biorthonormal)."""
    eigvals, V_R = eig(L)
    # Left eigenvectors via inverse. Rows of V_L_mat^dag give left eigvecs.
    V_L_mat = np.linalg.inv(V_R)
    # Normalise: rows V_L[k] satisfy V_L[k] V_R[:,k] = 1 already after inv.
    return eigvals, V_R, V_L_mat


# ---------------------------------------------------------------------------
# Palindrome and Tr(L) checks
# ---------------------------------------------------------------------------
def verify_tr_L(L, N, gamma_0):
    """Tr(L) should equal -d^2 * N * gamma_0."""
    d = 2**N
    expected = -(d**2) * N * gamma_0
    actual = float(np.trace(L).real)
    return actual, expected, abs(actual - expected)


def verify_palindromic_pairs(eigvals, N, gamma_0, tol=TOL_PALINDROME):
    """For uniform gamma_0, every eigenmode pairs with a partner such that
    alpha_s + alpha_partner = 2*N*gamma_0, where alpha = -Re(lambda).
    Also Im(lambda_s) + Im(lambda_partner) = 0 (up to sign for standing waves).

    Returns (n_paired, n_unpaired, max_residual).
    """
    expected_sum = 2.0 * N * gamma_0
    alphas = -eigvals.real
    imags = eigvals.imag
    n = len(eigvals)
    used = np.zeros(n, dtype=bool)
    residuals = []
    n_paired = 0
    for i in range(n):
        if used[i]:
            continue
        # Look for partner j with alpha_i + alpha_j ~ expected_sum
        target_a = expected_sum - alphas[i]
        target_im = -imags[i]
        best_j = -1
        best_r = np.inf
        for j in range(n):
            if used[j] or j == i:
                continue
            r = abs(alphas[j] - target_a) + abs(imags[j] - target_im)
            if r < best_r:
                best_r = r
                best_j = j
        if best_j >= 0 and best_r < 1e-6:
            used[i] = True
            used[best_j] = True
            n_paired += 2
            residuals.append(best_r)
    n_unpaired = n - n_paired
    max_res = max(residuals) if residuals else 0.0
    return n_paired, n_unpaired, max_res


# ---------------------------------------------------------------------------
# Initial states
# ---------------------------------------------------------------------------
def vacuum_ket(N):
    v = np.zeros(2**N, dtype=complex)
    v[0] = 1.0
    return v


def single_excitation_mode(N, k=1):
    """F65 bonding mode: psi_k = sqrt(2/(N+1)) sum_i sin(pi k (i+1)/(N+1)) |1_i>.
    Big-endian: qubit i excited -> state index 2^(N-1-i)."""
    psi = np.zeros(2**N, dtype=complex)
    norm = np.sqrt(2.0 / (N + 1))
    for i in range(N):
        amp = norm * np.sin(np.pi * k * (i + 1) / (N + 1))
        idx = 2**(N - 1 - i)
        psi[idx] = amp
    return psi


def bonding_plus_vacuum(N, k=1):
    """phi = (|vac> + psi_k)/sqrt(2), the PTF-standard initial state."""
    v = vacuum_ket(N)
    psi = single_excitation_mode(N, k)
    return (v + psi) / np.sqrt(2.0)


def density_matrix(ket):
    return np.outer(ket, ket.conj())


# ---------------------------------------------------------------------------
# Propagation
# ---------------------------------------------------------------------------
def propagate_vectorised(L, rho0, times):
    """Return rho(t) for each t in times, using L (d^2 x d^2) as superoperator."""
    d = rho0.shape[0]
    rho0_vec = rho0.flatten(order='F')  # column-stacking (matches kron convention)
    out = np.empty((len(times), d, d), dtype=complex)
    # Eigendecompose L once (dense), then use exponential via eigvals
    eigvals, V_R = eig(L)
    V_Linv = np.linalg.inv(V_R)
    c0 = V_Linv @ rho0_vec
    for k, t in enumerate(times):
        rho_vec_t = V_R @ (np.exp(eigvals * t) * c0)
        out[k] = rho_vec_t.reshape(d, d, order='F')
    return out


# ---------------------------------------------------------------------------
# Per-site purity
# ---------------------------------------------------------------------------
def partial_trace_keep_site(rho, i, N):
    """Reduced 2x2 density matrix on site i (big-endian), tracing out others."""
    d = 2**N
    rho_reshaped = rho.reshape([2] * (2 * N))
    # After reshape, axes are (q0, q1, ..., q_{N-1}, q0', q1', ..., q_{N-1}')
    # Trace out all except site i.
    keep = i
    # Permute so site i axes are at positions 0 and N
    # Trace out j != i
    # Implement by iterative einsum trick
    rho_i = rho.copy()
    # Work with the standard (d,d) form; compute using bit tricks
    result = np.zeros((2, 2), dtype=complex)
    for a in range(2):
        for b in range(2):
            # sum over all configs where site i is a (row) and b (column)
            s = 0.0
            for state in range(d):
                bit_i = (state >> (N - 1 - i)) & 1
                if bit_i != a:
                    continue
                state_b = state ^ ((a ^ b) << (N - 1 - i))
                # Sum over all "other" bits: iterate states consistent
                # Actually we need: rho[state_row, state_col] where row has bit i = a,
                # col has bit i = b, and all other bits match.
                for other_mask in range(2**(N - 1)):
                    # Build row and col indices by inserting a (resp. b) at position i
                    # and using other_mask for remaining bits.
                    def insert_bit(mask, bit, pos_from_left, N):
                        pos = N - 1 - pos_from_left
                        hi = (mask >> pos) << (pos + 1)
                        lo = mask & ((1 << pos) - 1)
                        return hi | (bit << pos) | lo
                    row = insert_bit(other_mask, a, i, N)
                    col = insert_bit(other_mask, b, i, N)
                    s += rho[row, col]
                break  # only need one pass with the bit trick
            result[a, b] = s
    return result


def partial_trace_keep_site_fast(rho, i, N):
    """Partial trace to 2x2 reduced state on site i (big-endian convention).

    rho has shape (2**N, 2**N) with site 0 as outermost (MSB) kron factor.
    Reshape to [2]*2N, axes 0..N-1 = ket bits, axes N..2N-1 = bra bits.
    Trace out every pair (j, N+j) with j != i by np.trace.
    """
    shape_2N = [2] * (2 * N)
    out = rho.reshape(shape_2N)
    # Build list of (ket_axis, bra_axis) pairs to trace; note: indexing
    # shifts as axes disappear. Trace in reverse order to keep lower
    # indices valid.
    to_trace = [(j, N + j) for j in range(N) if j != i]
    # Trace pairs with HIGHER bra axis first so that remaining axis indices
    # stay consistent. After trace of (j, N+j) with j < N+j, any axis
    # index > N+j shifts by 2, any axis between j and N+j shifts by 1.
    # Easiest: use np.trace iteratively, recomputing axis positions.
    ket_axes = list(range(N))       # current axis of site j ket
    bra_axes = list(range(N, 2 * N))  # current axis of site j bra
    for j in range(N - 1, -1, -1):
        if j == i:
            continue
        a_k = ket_axes[j]
        a_b = bra_axes[j]
        out = np.trace(out, axis1=a_k, axis2=a_b)
        # Every axis index > a_b shifts -1; every axis index > a_k shifts -1 more.
        lo, hi = sorted((a_k, a_b))
        for k in range(N):
            if k == j:
                continue
            if ket_axes[k] > hi:
                ket_axes[k] -= 2
            elif ket_axes[k] > lo:
                ket_axes[k] -= 1
            if bra_axes[k] > hi:
                bra_axes[k] -= 2
            elif bra_axes[k] > lo:
                bra_axes[k] -= 1
    # Now out has shape (2, 2). If remaining site-i ket axis is 1 and bra is 0,
    # transpose. Typically site-i ket ended up at axis 0, bra at axis 1.
    a_k = ket_axes[i]
    a_b = bra_axes[i]
    if a_k == 1 and a_b == 0:
        out = out.T
    return out


def per_site_purity(rho_traj, N):
    """Compute Tr(rho_i^2) for each site i and each time step."""
    T = rho_traj.shape[0]
    P = np.zeros((T, N))
    for t in range(T):
        rho_t = rho_traj[t]
        for i in range(N):
            rho_i = partial_trace_keep_site_fast(rho_t, i, N)
            P[t, i] = float(np.trace(rho_i @ rho_i).real)
    return P


# ---------------------------------------------------------------------------
# Alpha fit
# ---------------------------------------------------------------------------
def fit_alpha(t, pA_i, pB_i, t_max=T_FIT_MAX, bounds=ALPHA_BOUNDS):
    interp = interp1d(t, pA_i, bounds_error=False,
                      fill_value=(float(pA_i[0]), float(pA_i[-1])),
                      kind='cubic')
    mask = t <= t_max
    te = t[mask]
    b = pB_i[mask]

    def mse(a):
        d = interp(a * te) - b
        return float(np.mean(d * d))

    res = minimize_scalar(mse, bounds=bounds, method='bounded',
                          options={'xatol': 1e-6})
    return float(res.x), float(np.sqrt(res.fun))


# ---------------------------------------------------------------------------
# Main investigation
# ---------------------------------------------------------------------------
def run_for_N(N, J_mod_values=(1.1,), label=""):
    print(f"\n{'='*60}")
    print(f"N = {N}  (d = {2**N}, L is {4**N} x {4**N})")
    print(f"{'='*60}")
    result = {"N": N, "gamma_0": GAMMA_0, "T_final": T_FINAL,
              "N_steps": N_STEPS, "T_fit_max": T_FIT_MAX,
              "defect_bond": list(DEFECT_BOND), "J_uniform": J_UNIFORM,
              "J_mod_values": list(J_mod_values), "label": label,
              "phases": {}}

    d = 2**N
    # --- Phase 1-2: unperturbed chain ---
    J_A = [J_UNIFORM] * (N - 1)
    t0 = time.time()
    H_A = build_H_XY(J_A, N)
    L_A = build_liouvillian_matrix(H_A, GAMMA_0, N)
    eigvals_A, V_R_A, V_L_A = eig_LR(L_A)
    tr_A, tr_A_exp, tr_A_err = verify_tr_L(L_A, N, GAMMA_0)
    n_paired_A, n_unp_A, max_res_A = verify_palindromic_pairs(eigvals_A, N, GAMMA_0)
    print(f"\nPhase 1-2 (unperturbed):")
    print(f"  Tr(L_A) = {tr_A:.6f}  expected {tr_A_exp:.6f}  err {tr_A_err:.2e}")
    print(f"  Palindromic pairing: {n_paired_A}/{4**N} paired, "
          f"{n_unp_A} unpaired, max residual {max_res_A:.2e}")
    print(f"  Build+eigendecomp time: {time.time()-t0:.2f} s")
    result["phases"]["unperturbed"] = {
        "tr_L": tr_A, "tr_L_expected": tr_A_exp, "tr_L_err": tr_A_err,
        "n_paired": n_paired_A, "n_unpaired": n_unp_A,
        "max_palindrome_residual": max_res_A}

    # --- Phase 3: perturbations ---
    alphas_by_Jmod = {}
    closure_by_Jmod = {}
    for J_mod in J_mod_values:
        print(f"\nPhase 3-4 (J_mod = {J_mod}):")
        J_B = list(J_A)
        J_B[DEFECT_BOND[0]] = J_mod
        H_B = build_H_XY(J_B, N)
        L_B = build_liouvillian_matrix(H_B, GAMMA_0, N)
        eigvals_B, V_R_B, V_L_B = eig_LR(L_B)
        tr_B, tr_B_exp, tr_B_err = verify_tr_L(L_B, N, GAMMA_0)
        n_paired_B, n_unp_B, max_res_B = verify_palindromic_pairs(eigvals_B, N, GAMMA_0)
        print(f"  Tr(L_B) = {tr_B:.6f}  err {tr_B_err:.2e}  "
              f"(Tr preserved under J-perturbation: diff {abs(tr_A-tr_B):.2e})")
        print(f"  Palindromic pairing (B): {n_paired_B}/{4**N} paired, "
              f"max residual {max_res_B:.2e}")

        # --- Phase 4: propagate and fit alpha_i ---
        phi = bonding_plus_vacuum(N, k=1)
        rho_0 = density_matrix(phi)
        times = np.linspace(0.0, T_FINAL, N_STEPS + 1)
        t0 = time.time()
        rho_A = propagate_vectorised(L_A, rho_0, times)
        rho_B = propagate_vectorised(L_B, rho_0, times)
        P_A = per_site_purity(rho_A, N)
        P_B = per_site_purity(rho_B, N)
        alpha = np.zeros(N)
        rmse = np.zeros(N)
        for i in range(N):
            a, r = fit_alpha(times, P_A[:, i], P_B[:, i])
            alpha[i] = a
            rmse[i] = r
        closure = float(np.sum(np.log(alpha)))
        print(f"  Propagation + fit time: {time.time()-t0:.2f} s")
        print(f"  alpha_i:     " + "  ".join(f"{a:.4f}" for a in alpha))
        print(f"  fit RMSE:    " + "  ".join(f"{r:.1e}" for r in rmse))
        print(f"  Sum_i ln(alpha_i) = {closure:+.6f}")
        print(f"  Prod_i alpha_i    = {np.prod(alpha):.6f}")
        alphas_by_Jmod[J_mod] = alpha.tolist()
        closure_by_Jmod[J_mod] = closure

    result["phases"]["alpha_fit"] = {
        "alphas_by_Jmod": {f"{j}": a for j, a in alphas_by_Jmod.items()},
        "closure_by_Jmod": {f"{j}": c for j, c in closure_by_Jmod.items()}}

    # --- Phase 5 (Idea A): determinant / trace invariants ---
    print(f"\nPhase 5 (Idea A: determinant/trace invariants):")
    # Tr(L) is J-independent (Hamiltonian commutator traceless).
    # det(exp(L t)) = exp(Tr(L) t) is J-independent.
    # Question: does this give a constraint on Prod_i alpha_i?
    # The alpha_i come from per-site reduced-state dynamics, not the full propagator.
    # Compute the "per-site" eigenvalue contribution: for each site i, sum of
    # eigenvalues weighted by their coupling to P_i at leading order.
    # Leading-order: dP_i/dt|_0 = -2 <rho_i(0)| Tr_{not i}(L rho_0)> (a scalar).
    # For the uniform chain we compute this and compare to the perturbed.
    # This is the "derivative" alpha_i: (dP_B/dt)/(dP_A/dt) at t=0.
    # If this ratio matches the fitted alpha_i, the fit captures a local rate.
    for J_mod in J_mod_values:
        J_B = list(J_A); J_B[DEFECT_BOND[0]] = J_mod
        H_B = build_H_XY(J_B, N)
        L_B = build_liouvillian_matrix(H_B, GAMMA_0, N)
        phi = bonding_plus_vacuum(N, k=1)
        rho_0 = density_matrix(phi)
        rho_0_vec = rho_0.flatten(order='F')
        dr_A = (L_A @ rho_0_vec).reshape(2**N, 2**N, order='F')
        dr_B = (L_B @ rho_0_vec).reshape(2**N, 2**N, order='F')
        # dP_i/dt|_0 = 2 Re Tr(rho_i(0) * dr_i(0))
        dP_A = np.zeros(N); dP_B = np.zeros(N)
        for i in range(N):
            rho_i0 = partial_trace_keep_site_fast(rho_0, i, N)
            drA_i = partial_trace_keep_site_fast(dr_A, i, N)
            drB_i = partial_trace_keep_site_fast(dr_B, i, N)
            dP_A[i] = 2.0 * float(np.trace(rho_i0 @ drA_i).real)
            dP_B[i] = 2.0 * float(np.trace(rho_i0 @ drB_i).real)
        # "Derivative alpha_i": ratio of initial purity derivatives.
        # Well-defined only where dP_A[i] is nonzero.
        alpha_deriv = np.where(np.abs(dP_A) > 1e-12, dP_B / dP_A, np.nan)
        print(f"  J_mod = {J_mod}: dP_A/dt|_0 = " +
              "  ".join(f"{v:+.3e}" for v in dP_A))
        print(f"                   dP_B/dt|_0 = " +
              "  ".join(f"{v:+.3e}" for v in dP_B))
        print(f"                   alpha_deriv (dP_B/dP_A) = " +
              "  ".join(f"{v:+.4f}" if np.isfinite(v) else "   NaN"
                       for v in alpha_deriv))
        if np.all(np.isfinite(alpha_deriv)):
            closure_deriv = float(np.sum(np.log(np.abs(alpha_deriv))))
            print(f"                   Sum_i ln|alpha_deriv_i| = {closure_deriv:+.4f}")

    # --- Phase 6 (Idea B): voice vs memory split ---
    # Decompose the Liouvillian eigenmodes into voice (alpha > 0) and memory
    # (alpha = 0, i.e., stationary). Check what fraction of the initial state
    # overlaps with voice modes, and how that fraction shifts under J-perturb.
    print(f"\nPhase 6 (Idea B: voice vs memory modes):")
    alpha_A = -eigvals_A.real
    voice_mask_A = alpha_A > 1e-10
    n_voice_A = int(np.sum(voice_mask_A))
    n_mem_A = len(eigvals_A) - n_voice_A
    phi = bonding_plus_vacuum(N, k=1)
    rho_0 = density_matrix(phi)
    rho_0_vec = rho_0.flatten(order='F')
    c_A = V_L_A @ rho_0_vec  # left-eigenvector overlaps
    voice_weight_A = float(np.sum(np.abs(c_A[voice_mask_A])**2).real)
    mem_weight_A = float(np.sum(np.abs(c_A[~voice_mask_A])**2).real)
    total = voice_weight_A + mem_weight_A
    print(f"  Unperturbed: {n_voice_A} voice modes, {n_mem_A} memory modes")
    print(f"  |c|^2 on voice: {voice_weight_A:.4f}   on memory: "
          f"{mem_weight_A:.4f}   (total {total:.4f})")
    # F4 stationary count prediction: for XY chain with Z-dephasing, there are
    # N+1 strict stationary modes (excitation-sector projectors).
    # But "memory" (alpha = 0) includes I operator contributions. Let us check.
    n_strict_zero_A = int(np.sum(np.abs(alpha_A) < 1e-10))
    print(f"  Strict zero modes (alpha = 0 exactly): {n_strict_zero_A}  "
          f"(F4 predicts N+1 = {N+1})")

    # --- Phase 7 (Idea D): resolvent trace per site ---
    # The Liouvillian resolvent G(z) = (z I - L)^{-1}. Its trace
    # Tr(G(z)) = Sum 1/(z - lambda_s) has poles at eigenvalues.
    # "Per-site" projection: for the operator-space basis, project onto
    # operators supported on site i. For a 2x2 per-site operator,
    # the relevant projector is P_i = sum over 3 non-identity Paulis:
    # Pauli_k at site i otimes (I/2)^(N-1) elsewhere.
    # Tr[G(z) P_i] tracks how site-i local observables respond to drive.
    # Test: does the product over i of some resolvent-derived quantity give
    # a J-independent constant?
    print(f"\nPhase 7 (Idea D: resolvent trace per site):")
    # Use a few probe frequencies. Since L is purely dissipative (Re(lambda) <= 0),
    # use z in the upper half plane just above the real axis.
    probe_z = [0.01 + 0j, 0.01 + 0.1j, 0.01 - 0.1j]
    d = 2**N
    for z in probe_z:
        GA = np.linalg.inv(z * np.eye(4**N) - L_A)
        trGA = complex(np.trace(GA))
        # Build per-site Pauli projectors in vectorised operator space.
        # Operator basis: |E_k> = vec(E_k) where E_k is a d x d matrix.
        # We work with site-i Pauli operators as vectors.
        site_traces = []
        for i in range(N):
            # Site-i Paulis (excluding I_i): X_i, Y_i, Z_i. Use vec.
            s_sum = 0.0 + 0j
            for P_local in [X, Y, Z]:
                Pi = site_op(P_local, i, N)
                vPi = Pi.flatten(order='F')
                s_sum += vPi.conj() @ GA @ vPi
            site_traces.append(s_sum)
        # For each J_mod compare
        for J_mod in J_mod_values:
            J_B = list(J_A); J_B[DEFECT_BOND[0]] = J_mod
            L_B_z = build_liouvillian_matrix(build_H_XY(J_B, N), GAMMA_0, N)
            GB = np.linalg.inv(z * np.eye(4**N) - L_B_z)
            trGB = complex(np.trace(GB))
            site_traces_B = []
            for i in range(N):
                s_sum = 0.0 + 0j
                for P_local in [X, Y, Z]:
                    Pi = site_op(P_local, i, N)
                    vPi = Pi.flatten(order='F')
                    s_sum += vPi.conj() @ GB @ vPi
                site_traces_B.append(s_sum)
            ratio = [complex(b) / complex(a) if abs(a) > 1e-12 else float('nan')
                     for a, b in zip(site_traces, site_traces_B)]
            print(f"  z = {z}: Tr(G_A) = {trGA:.4f}  Tr(G_B) = {trGB:.4f}  "
                  f"diff = {abs(trGA-trGB):.2e}")
            print(f"    per-site |Tr_i(G)|: A = " +
                  "  ".join(f"{abs(s):.3f}" for s in site_traces))
            print(f"    per-site |Tr_i(G)|: B = " +
                  "  ".join(f"{abs(s):.3f}" for s in site_traces_B))
            # Log-ratio (real part) - would be the resolvent analog of ln(alpha_i)
            log_ratios = [np.log(abs(b) / abs(a)) if abs(a) > 1e-12 else np.nan
                          for a, b in zip(site_traces, site_traces_B)]
            if all(np.isfinite(r) for r in log_ratios):
                sum_log = sum(log_ratios)
                print(f"    Sum_i ln|Tr_i(G_B)/Tr_i(G_A)| = {sum_log:+.4f}")
            break  # one J_mod per z is enough

    # Save
    out_path = RESULTS_DIR / f"n{N}_closure_investigation.json"
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2, default=str)
    print(f"\nSaved: {out_path}")
    return result


def run_fine_scan(N):
    """Fine scan of delta_J around zero to separate 1st and 2nd order closure
    contributions."""
    print(f"\n{'='*60}")
    print(f"Fine delta_J scan at N={N}")
    print(f"{'='*60}")
    d = 2**N
    J_A = [J_UNIFORM] * (N - 1)
    H_A = build_H_XY(J_A, N)
    L_A = build_liouvillian_matrix(H_A, GAMMA_0, N)
    phi = bonding_plus_vacuum(N, k=1)
    rho_0 = density_matrix(phi)
    times = np.linspace(0.0, T_FINAL, N_STEPS + 1)
    rho_A = propagate_vectorised(L_A, rho_0, times)
    P_A = per_site_purity(rho_A, N)
    dJ_values = [-0.2, -0.15, -0.1, -0.05, -0.02, -0.01, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2]
    print(f"  dJ       closure Σln(α)    per-site α_i")
    rows = []
    for dJ in dJ_values:
        J_B = list(J_A); J_B[DEFECT_BOND[0]] = J_UNIFORM + dJ
        L_B = build_liouvillian_matrix(build_H_XY(J_B, N), GAMMA_0, N)
        rho_B = propagate_vectorised(L_B, rho_0, times)
        P_B = per_site_purity(rho_B, N)
        alpha = np.zeros(N)
        for i in range(N):
            a, _ = fit_alpha(times, P_A[:, i], P_B[:, i])
            alpha[i] = a
        closure = float(np.sum(np.log(alpha)))
        print(f"  {dJ:+.3f}   {closure:+.6f}   " +
              " ".join(f"{a:.4f}" for a in alpha))
        rows.append({"dJ": dJ, "closure": closure, "alpha": alpha.tolist()})
    # Fit closure ~ a*dJ + b*dJ^2
    dJ_arr = np.array([r["dJ"] for r in rows])
    cl_arr = np.array([r["closure"] for r in rows])
    # Polynomial fit of degree 3
    coeffs = np.polyfit(dJ_arr, cl_arr, 3)
    print(f"  Fit closure(dJ) = {coeffs[3]:+.5f} + {coeffs[2]:+.4f}·dJ"
          f" + {coeffs[1]:+.4f}·dJ² + {coeffs[0]:+.4f}·dJ³")
    print(f"  Order-of-magnitude 1st/2nd order ratio at dJ=0.1: "
          f"|lin| = {abs(coeffs[2]*0.1):.4f}   |quad| = {abs(coeffs[1]*0.01):.4f}")
    return {"dJ": dJ_arr.tolist(), "closure": cl_arr.tolist(),
            "alpha": [r["alpha"] for r in rows],
            "polyfit_cubic": coeffs.tolist()}


def run_resolvent_scan(N):
    """Scan z along the real axis and check if Sum_i ln|Tr_i(G_B)/Tr_i(G_A)| ~ 0."""
    print(f"\n{'='*60}")
    print(f"Resolvent z-scan at N={N}")
    print(f"{'='*60}")
    d = 2**N
    J_A = [J_UNIFORM] * (N - 1)
    L_A = build_liouvillian_matrix(build_H_XY(J_A, N), GAMMA_0, N)
    J_B = list(J_A); J_B[DEFECT_BOND[0]] = 1.0 + 0.1
    L_B = build_liouvillian_matrix(build_H_XY(J_B, N), GAMMA_0, N)

    # Precompute site-Pauli vectors
    def site_pauli_vecs(N):
        vecs = []
        for i in range(N):
            for P_local in [X, Y, Z]:
                Pi = site_op(P_local, i, N)
                vecs.append((i, Pi.flatten(order='F')))
        return vecs
    site_vecs = site_pauli_vecs(N)

    z_values = [0.001, 0.005, 0.01, 0.03, 0.05, 0.1, 0.2, 0.5, 1.0]
    rows = []
    for z in z_values:
        zz = z + 0j
        I_big = np.eye(4**N)
        GA = np.linalg.solve(zz * I_big - L_A, I_big)
        GB = np.linalg.solve(zz * I_big - L_B, I_big)
        # Per-site |Tr_i(G)| = sum over X,Y,Z vectors at site i
        per_site_A = np.zeros(N)
        per_site_B = np.zeros(N)
        for i, vPi in site_vecs:
            per_site_A[i] += float(abs(vPi.conj() @ GA @ vPi))
            per_site_B[i] += float(abs(vPi.conj() @ GB @ vPi))
        log_ratios = np.log(per_site_B / per_site_A)
        sum_log = float(np.sum(log_ratios))
        trGA = complex(np.trace(GA))
        trGB = complex(np.trace(GB))
        print(f"  z = {z:6.3f}  Σ_i ln(|Tr_i(G_B)|/|Tr_i(G_A)|) = {sum_log:+.5f}"
              f"   Tr(G_B)-Tr(G_A) = {abs(trGB-trGA):.2e}")
        rows.append({"z": z, "sum_log": sum_log,
                     "per_site_A": per_site_A.tolist(),
                     "per_site_B": per_site_B.tolist(),
                     "log_ratios": log_ratios.tolist()})
    return rows


def run_multiple_initial_states(N):
    """Check closure at N=3 or N=5 with multiple initial states, including
    single-excitation ψ_k for k=1,2 and |+>^N."""
    print(f"\n{'='*60}")
    print(f"Initial state variation at N={N}")
    print(f"{'='*60}")
    d = 2**N
    J_A = [J_UNIFORM] * (N - 1)
    L_A = build_liouvillian_matrix(build_H_XY(J_A, N), GAMMA_0, N)
    J_B = list(J_A); J_B[DEFECT_BOND[0]] = 1.1
    L_B = build_liouvillian_matrix(build_H_XY(J_B, N), GAMMA_0, N)
    times = np.linspace(0.0, T_FINAL, N_STEPS + 1)

    states = {}
    # ψ_1 + vac (PTF standard)
    phi1 = bonding_plus_vacuum(N, k=1)
    states["phi_1_plus_vac"] = phi1
    # ψ_2 + vac (PTF falsification state)
    if N >= 2:
        phi2 = (vacuum_ket(N) + single_excitation_mode(N, k=2)) / np.sqrt(2.0)
        states["phi_2_plus_vac"] = phi2
    # ψ_1 only (pure single-excitation)
    states["psi_1_only"] = single_excitation_mode(N, k=1)
    # |+>^N (multi-sector state)
    plus_N = np.ones(d, dtype=complex) / np.sqrt(d)
    states["plus_N"] = plus_N
    # Bell-like: (|00..0> + |11..1>)/sqrt(2)
    ghz = np.zeros(d, dtype=complex); ghz[0] = 1; ghz[-1] = 1
    ghz /= np.sqrt(2.0)
    states["ghz"] = ghz

    results = {}
    for name, ket in states.items():
        rho_0 = density_matrix(ket)
        rho_A = propagate_vectorised(L_A, rho_0, times)
        rho_B = propagate_vectorised(L_B, rho_0, times)
        P_A = per_site_purity(rho_A, N)
        P_B = per_site_purity(rho_B, N)
        alpha = np.zeros(N); rmse = np.zeros(N)
        for i in range(N):
            a, r = fit_alpha(times, P_A[:, i], P_B[:, i])
            alpha[i] = a; rmse[i] = r
        closure = float(np.sum(np.log(alpha)))
        print(f"  {name:20s}  Σ ln α = {closure:+.5f}   "
              f"max RMSE = {rmse.max():.1e}   α = " +
              " ".join(f"{a:.4f}" for a in alpha))
        results[name] = {"alpha": alpha.tolist(), "rmse": rmse.tolist(),
                         "closure": closure}
    return results


def main():
    print("="*60)
    print("Pi-pair closure investigation")
    print(f"  gamma_0 = {GAMMA_0}")
    print(f"  defect bond = {DEFECT_BOND}, delta_J = +{DELTA_J}")
    print(f"  T_final = {T_FINAL}, T_fit_max = {T_FIT_MAX}")
    print(f"  Results -> {RESULTS_DIR}")

    all_results = {}
    # N=3 probe
    all_results["n3"] = run_for_N(3, J_mod_values=(0.9, 1.0, 1.1), label="N=3 probe")
    # Fine dJ scan at N=3
    all_results["n3_fine"] = run_fine_scan(3)
    # Resolvent scan at N=3
    all_results["n3_resolvent"] = run_resolvent_scan(3)
    # Multiple initial states at N=3
    all_results["n3_states"] = run_multiple_initial_states(3)

    # N=5 for validation
    all_results["n5"] = run_for_N(5, J_mod_values=(0.9, 1.0, 1.1), label="N=5 validation")
    # Fine dJ scan at N=5
    all_results["n5_fine"] = run_fine_scan(5)
    # Resolvent scan at N=5
    all_results["n5_resolvent"] = run_resolvent_scan(5)
    # Multiple initial states at N=5
    all_results["n5_states"] = run_multiple_initial_states(5)

    # Save combined results
    out_path = RESULTS_DIR / "summary.json"
    with open(out_path, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\n\nSummary saved: {out_path}")


if __name__ == "__main__":
    main()
