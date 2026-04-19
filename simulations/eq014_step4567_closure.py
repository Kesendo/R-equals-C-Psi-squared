#!/usr/bin/env python3
"""eq014_step4567_closure.py

EQ-014 Steps 4-7: V_L construction, first-order mixing, bilinear purity,
alpha_i extraction, closure law test.

Dependencies: eq014_step23_biorth.py must have been run first to produce
- eq014_biorth_fix.npz       (c_diag + degenerate clusters + slow indices)
- eq014_biorth_metadata.json

Input files (from ptf mode):
- eq014_eigvals_n7.bin        (16384 complex)
- eq014_right_eigvecs_n7.bin  (~4 GB, column-major)
- eq014_left_eigvecs_n7.bin   (~4 GB, column-major)

Pipeline:
  4. Build V_L = -i [H_pert, .] for H_pert = (1/2)(X_b X_{b+1} + Y_b Y_{b+1})
     for each target bond (default (0,1) and (3,4)).
  5. Compute mixing coefficients A_{s', s} = <W_{s'} | V_L | M_s>
     for every slow s and every mode s'.
  6. Bilinear purity expansion:
     P_i(t) = Σ_{s,s'} c_s c*_{s'} e^{(λ_s + λ*_{s'}) t} Tr(m_{s,i} m_{s',i}^†)
     Apply first-order perturbation (δc_s, δm_{s,i}) from the mixing, then
     compute perturbed P_B^i(t). Extract α_i from P_B(t) ≈ P_A(α_i t).
  7. Closure law test: Σ_i ln(α_i) per initial state and bond.

Outputs:
- eq014_mixing_coefficients_bond{0_1,3_4}.npz    (A matrix + slow eigenvalues)
- eq014_alpha_prediction.json                    (α_i per state, bond, site)
- eq014_closure_test.txt                         (Σ_i ln α_i + verdict)

Date: 2026-04-19
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sps

sys.stdout.reconfigure(encoding='utf-8')

RESULTS_DIR = Path(__file__).parent / "results"
DELTA_J = 0.1                  # perturbation strength δJ on the target bond
TIME_WINDOW = np.linspace(0.05, 20.0, 80)  # matches PTF's T_FIT = 20
BONDS = [(0, 1), (3, 4)]       # target bonds for the defect
CLOSURE_TARGET = 1e-6          # success if Σ_i ln(α_i) < this
CLOSURE_PARTIAL = 1e-2         # partial support threshold
DEGEN_EXCLUDE = 1e-10          # exclude mixing to eigenvalues within this of λ_s
                               # (within-cluster contributions handled by degenerate PT:
                               # step23's SVD biorth aligns the subspace so first-order
                               # mixing only happens across clusters)

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def load_bin(path, shape):
    data = np.fromfile(path, dtype=np.complex128)
    assert data.size == np.prod(shape), f"{path.name}: size {data.size} != {np.prod(shape)}"
    return data.reshape(shape, order='F')


def site_operator_sparse(op, site, N):
    """Sparse N-qubit site operator I⊗...⊗op⊗...⊗I with op on `site`
    (site 0 = MSB). Matches PauliOps.At / BuildDirectRaw convention."""
    op_s = sps.csr_matrix(op, dtype=complex)
    I2_s = sps.eye(2, dtype=complex, format='csr')
    full = sps.eye(1, dtype=complex, format='csr')
    for k in range(N):
        full = sps.kron(full, op_s if k == site else I2_s, format='csr')
    return full


def two_site_bond_H(bond, N):
    """H_bond = (1/2) (X_b X_{b+1} + Y_b Y_{b+1}) as sparse d x d matrix."""
    b, b_plus = bond
    XX = site_operator_sparse(X, b, N) @ site_operator_sparse(X, b_plus, N)
    YY = site_operator_sparse(Y, b, N) @ site_operator_sparse(Y, b_plus, N)
    return 0.5 * (XX + YY)


def build_V_L(bond, N):
    """Sparse V_L = -i (H_bond ⊗ I - I ⊗ H_bond^T) as d^2 x d^2.
    vec(ρ)[a*d + b] = ρ[a, b], so L * vec(ρ) = vec(-i[H, ρ])."""
    d = 2 ** N
    Hb = two_site_bond_H(bond, N)
    Id = sps.eye(d, dtype=complex, format='csr')
    return (-1j * (sps.kron(Hb, Id, format='csr')
                   - sps.kron(Id, Hb.T, format='csr'))).tocsr()


def biorthogonalize_cluster(W_sub, V_sub):
    """SVD-based biorthogonalization; same transform as in step23."""
    M = W_sub.conj().T @ V_sub
    U, s, Vh = np.linalg.svd(M)
    sqrt_s_inv = np.diag(1.0 / np.sqrt(s))
    W_new = W_sub @ (U @ sqrt_s_inv)
    V_new = V_sub @ (Vh.conj().T @ sqrt_s_inv)
    return W_new, V_new


def iter_clusters(cluster_flat, cluster_offsets):
    for k in range(len(cluster_offsets) - 1):
        yield cluster_flat[cluster_offsets[k]:cluster_offsets[k + 1]]


def apply_biorth_fix(R, W, c_diag, cluster_flat, cluster_offsets):
    """Apply normalization R /= c_diag + in-cluster SVD biorth fix, in place."""
    safe_c = np.where(np.abs(c_diag) < 1e-14, 1.0, c_diag)
    R /= safe_c[np.newaxis, :]
    for cl in iter_clusters(cluster_flat, cluster_offsets):
        if len(cl) < 2:
            continue
        W_sub = W[:, cl].copy()
        R_sub = R[:, cl].copy()
        W_new, R_new = biorthogonalize_cluster(W_sub, R_sub)
        W[:, cl] = W_new
        R[:, cl] = R_new


def mode_to_matrix(mode_vec, d):
    """vec(ρ)[a*d + b] = ρ[a, b] -> reshape to (d, d) row-major."""
    return mode_vec.reshape((d, d))


def partial_trace_all_sites(M, d, N):
    """Single-site marginals of d x d matrix M, returned as (N, 2, 2).
    Big-endian site convention: site 0 is MSB of the basis index.
    """
    out = np.zeros((N, 2, 2), dtype=complex)
    shape = [2] * (2 * N)
    T = M.reshape(shape)
    letters = "abcdefghijklmnop"
    for i in range(N):
        row_labels = list(letters[:N])
        col_labels = list(letters[N:2 * N])
        for j in range(N):
            if j != i:
                col_labels[j] = row_labels[j]  # trace j
        in_spec = "".join(row_labels) + "".join(col_labels)
        out_spec = row_labels[i] + col_labels[i]
        out[i] = np.einsum(f"{in_spec}->{out_spec}", T)
    return out


def build_initial_states(N, d):
    """Five initial states matching the PTF empirical data.
    Big-endian: site i corresponds to bit (N-1-i) of the basis index."""
    vac = np.zeros(d, dtype=complex)
    vac[0] = 1.0

    def single_excitation(k):
        psi = np.zeros(d, dtype=complex)
        norm = np.sqrt(2.0 / (N + 1))
        for i in range(N):
            amp = norm * np.sin(np.pi * k * (i + 1) / (N + 1))
            idx = 1 << (N - 1 - i)
            psi[idx] = amp
        return psi

    states = {}
    for k in range(1, 5):
        psi_k = single_excitation(k)
        phi = (vac + psi_k) / np.sqrt(2)
        rho0 = np.outer(phi, phi.conj())
        states[f"psi_{k}_bonding"] = rho0.flatten()

    plus = np.full(d, 1.0 / np.sqrt(2 ** N), dtype=complex)
    rho_plus = np.outer(plus, plus.conj())
    states["plus_7"] = rho_plus.flatten()
    return states


def compute_unperturbed_purity_exact(c_all, all_vals, R, d, N, times):
    """EXACT unperturbed purity P_A_i(t) via full-mode direct evolution.
    Uses ALL 16384 modes (not slow-only) to preserve accuracy at early times.
    """
    P = np.zeros((N, len(times)), dtype=float)
    for ti, t in enumerate(times):
        rho_vec = R @ (c_all * np.exp(all_vals * t))
        rho_mat = rho_vec.reshape((d, d))
        margs = partial_trace_all_sites(rho_mat, d, N)
        for i in range(N):
            P[i, ti] = float(np.real(np.trace(margs[i] @ margs[i].conj().T)))
    return P


def compute_perturbed_purity_exact(c_all, delta_c_slow, slow_indices, R,
                                   delta_R_slow, all_vals, slow_vals,
                                   d, N, times):
    """First-order δP_i(t) / δJ with FULL ρ_A (all modes) on one side
    and SLOW δρ (first-order correction) on the other:

        δP_i = 2 Re Tr[ρ_{A,i}(t) · (δρ_i(t))^†]

    where ρ_A,i is exact and δρ = R_slow @ (δc ⊙ exp) + δR_slow @ (c ⊙ exp).

    The bilinear-only-on-slow expansion used before was incorrect for
    bonding modes whose slow-subspace weight is ~70 %.
    """
    c_slow = c_all[slow_indices]
    R_slow = R[:, slow_indices]
    dP = np.zeros((N, len(times)), dtype=float)
    for ti, t in enumerate(times):
        rho_A_vec = R @ (c_all * np.exp(all_vals * t))
        rho_A_mat = rho_A_vec.reshape((d, d))
        exp_slow = np.exp(slow_vals * t)
        drho_from_dc = R_slow @ (delta_c_slow * exp_slow)
        drho_from_dv = delta_R_slow @ (c_slow * exp_slow)
        drho_vec = drho_from_dc + drho_from_dv
        drho_mat = drho_vec.reshape((d, d))
        margs_A = partial_trace_all_sites(rho_A_mat, d, N)
        margs_d = partial_trace_all_sites(drho_mat, d, N)
        for i in range(N):
            dP[i, ti] = 2.0 * float(np.real(
                np.trace(margs_A[i] @ margs_d[i].conj().T)))
    return dP


def fit_alpha(P_A, P_B, times, zero_dynamics_tol=1e-8,
              alpha_bounds=(0.1, 10.0)):
    """Bounded scalar fit α_i minimising ⟨[P_A(α t) − P_B(t)]²⟩ over the
    time window. Matches the PTF-style one-parameter time-rescale fit in
    observer_time_rescale.py. Falls back to α = 1 for sites with essentially-
    constant P_A (e.g. node-sites of ψ_2 / ψ_4) where the rescaling is
    operationally undefined.
    """
    from scipy.interpolate import interp1d
    from scipy.optimize import minimize_scalar

    N, T = P_A.shape
    alpha = np.zeros(N)
    for i in range(N):
        pa_range = np.max(P_A[i]) - np.min(P_A[i])
        pab_range = np.max(np.abs(P_B[i] - P_A[i]))
        if pa_range < zero_dynamics_tol and pab_range < zero_dynamics_tol:
            alpha[i] = 1.0
            continue
        # interp1d over expanded time axis: we may evaluate P_A at α*t
        # which can exceed times[-1]. Extrapolate with the last value.
        interp = interp1d(times, P_A[i], kind='cubic', bounds_error=False,
                          fill_value=(P_A[i, 0], P_A[i, -1]))

        def mse(a):
            d = interp(a * times) - P_B[i]
            return float(np.mean(d * d))

        try:
            res = minimize_scalar(mse, bounds=alpha_bounds, method='bounded')
            alpha[i] = float(res.x)
        except Exception:
            alpha[i] = np.nan
    return alpha


def main():
    t_total = time.time()
    print("=== EQ-014 Steps 4-7: closure law ===")

    with open(RESULTS_DIR / "eq014_biorth_metadata.json") as f:
        meta = json.load(f)
    N = meta["N"]
    d = meta["d"]
    d2 = meta["d2"]
    gamma_0 = meta["gamma"]
    slow_cutoff = meta["slow_cutoff"]
    print(f"N={N}, d={d}, d2={d2}, γ={gamma_0}, slow cutoff={slow_cutoff}")

    fix = np.load(RESULTS_DIR / "eq014_biorth_fix.npz")
    c_diag = fix["c_diag"]
    cluster_flat = fix["cluster_flat"]
    cluster_offsets = fix["cluster_offsets"]
    slow_indices = fix["slow_indices"]
    n_clusters = len(cluster_offsets) - 1
    print(f"Slow modes: {len(slow_indices)}, clusters: {n_clusters}")
    print()

    vals = np.fromfile(RESULTS_DIR / "eq014_eigvals_n7.bin", dtype=np.complex128)
    print("Loading right eigenvectors (4 GB)...")
    t0 = time.time()
    R = load_bin(RESULTS_DIR / "eq014_right_eigvecs_n7.bin", (d2, d2))
    print(f"  {time.time() - t0:.1f} s")
    print("Loading left eigenvectors (4 GB)...")
    t0 = time.time()
    W = load_bin(RESULTS_DIR / "eq014_left_eigvecs_n7.bin", (d2, d2))
    print(f"  {time.time() - t0:.1f} s")
    print()

    print("Applying biorth normalization + cluster fix...")
    t0 = time.time()
    apply_biorth_fix(R, W, c_diag, cluster_flat, cluster_offsets)
    print(f"  {time.time() - t0:.1f} s")
    print()

    R_slow = R[:, slow_indices]
    W_slow = W[:, slow_indices]
    slow_vals = vals[slow_indices]
    n_slow = len(slow_indices)
    print(f"R_slow {R_slow.shape}, W_slow {W_slow.shape}")

    print("Computing slow-mode single-site marginals...")
    t0 = time.time()
    slow_marginals = np.zeros((n_slow, N, 2, 2), dtype=complex)
    for idx in range(n_slow):
        M = mode_to_matrix(R_slow[:, idx], d)
        slow_marginals[idx] = partial_trace_all_sites(M, d, N)
    print(f"  {time.time() - t0:.1f} s -> shape {slow_marginals.shape}")
    print()

    states = build_initial_states(N, d)
    print(f"Initial states: {list(states.keys())}")

    # c_all = <u_s | ρ_0> for ALL modes (needed for exact ρ_A(t))
    # c_slow = c_all[slow_indices] (used for δρ)
    c_all_by_state = {name: W.conj().T @ rho for name, rho in states.items()}
    c_by_state = {name: c_all_by_state[name][slow_indices]
                  for name in states}

    print("Computing unperturbed P_A_i(t) via EXACT direct evolution (all modes)...")
    t0 = time.time()
    PA_by_state = {name: compute_unperturbed_purity_exact(
        c_all_by_state[name], vals, R, d, N, TIME_WINDOW)
        for name in states}
    print(f"  {time.time() - t0:.1f} s")
    print()

    results = {
        "N": N, "d": d, "d2": d2, "gamma": gamma_0,
        "slow_cutoff": slow_cutoff, "delta_J": DELTA_J,
        "slow_mode_count": int(n_slow),
        "initial_states": list(states.keys()),
        "bonds": [],
    }

    closure_log_lines = [
        "EQ-014 Closure Law Test (Steps 4-7)",
        f"Date: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        f"N={N}, d2={d2}, gamma={gamma_0}",
        f"Slow modes: {n_slow} (|Re λ| ≤ {slow_cutoff})",
        f"δJ = {DELTA_J}",
        "",
    ]

    for bond in BONDS:
        bond_tag = f"bond_{bond[0]}_{bond[1]}"
        print(f"--- Bond {bond} ---")

        print(f"Building V_L for bond {bond}...")
        t0 = time.time()
        V_L = build_V_L(bond, N)
        frob = float(np.sqrt((np.abs(V_L.data) ** 2).sum()))
        print(f"  V_L shape {V_L.shape} in {time.time() - t0:.1f} s, "
              f"nnz={V_L.nnz}, |V_L|_F = {frob:.3f}")

        print(f"Computing V_L @ R_slow...")
        t0 = time.time()
        VR = V_L @ R_slow
        print(f"  {time.time() - t0:.1f} s")

        print(f"A = W^H @ VR  (projecting onto all W)...")
        t0 = time.time()
        A = W.conj().T @ VR
        print(f"  {time.time() - t0:.1f} s -> A shape {A.shape}")

        mix_path = RESULTS_DIR / f"eq014_mixing_coefficients_{bond_tag}.npz"
        np.savez(mix_path,
                 A=A,
                 slow_indices=slow_indices,
                 slow_eigenvalues=slow_vals,
                 all_eigenvalues=vals,
                 bond=np.array(bond))
        print(f"  Wrote {mix_path.name}")

        # Degenerate-cluster handling: mask out all s' whose eigenvalue is
        # within DEGEN_EXCLUDE of λ_s. The SVD biorth fix in step23 already
        # aligned the subspace within each cluster; within-cluster mixing is
        # absorbed there and should not appear in first-order PT.
        print("Computing δM_s eigenvector corrections...")
        t0 = time.time()
        delta_R_slow = np.zeros((d2, n_slow), dtype=complex)
        n_excluded_total = 0
        for s_pos, s_glob in enumerate(slow_indices):
            denom = slow_vals[s_pos] - vals
            mask_tiny = np.abs(denom) < DEGEN_EXCLUDE
            n_excluded_total += int(np.sum(mask_tiny))
            denom_safe = np.where(mask_tiny, 1.0, denom)
            coeffs = A[:, s_pos] / denom_safe
            coeffs[mask_tiny] = 0.0
            delta_R_slow[:, s_pos] = R @ coeffs
        avg_excluded = n_excluded_total / max(1, n_slow)
        print(f"  {time.time() - t0:.1f} s (avg {avg_excluded:.1f} modes excluded per slow s)")

        print("Computing δ single-site marginals...")
        t0 = time.time()
        delta_slow_marginals = np.zeros((n_slow, N, 2, 2), dtype=complex)
        for idx in range(n_slow):
            dM = mode_to_matrix(delta_R_slow[:, idx], d)
            delta_slow_marginals[idx] = partial_trace_all_sites(dM, d, N)
        print(f"  {time.time() - t0:.1f} s")

        print("Computing B = W_slow^H V_L R for left corrections...")
        t0 = time.time()
        WV = W_slow.conj().T @ V_L
        B = WV @ R
        print(f"  {time.time() - t0:.1f} s -> B shape {B.shape}")

        # δc_s = Σ_{s'≠s} [u_s^H V_L v_{s'} / (λ_s - λ_{s'})] · c_all[s']
        # Same DEGEN_EXCLUDE masking as δM_s.
        print("Computing δc_s per initial state...")
        t0 = time.time()
        delta_c_by_state = {}
        for name, rho_flat in states.items():
            c_all = W.conj().T @ rho_flat
            dc = np.zeros(n_slow, dtype=complex)
            for s_pos, s_glob in enumerate(slow_indices):
                denom = slow_vals[s_pos] - vals
                mask_tiny = np.abs(denom) < DEGEN_EXCLUDE
                denom_safe = np.where(mask_tiny, 1.0, denom)
                term = B[s_pos, :] / denom_safe * c_all
                term[mask_tiny] = 0.0
                dc[s_pos] = np.sum(term)
            delta_c_by_state[name] = dc
        print(f"  {time.time() - t0:.1f} s")

        print("Computing δP_i(t) per state (FULL ρ_A ⊗ slow δρ)...")
        t0 = time.time()
        dP_by_state = {
            name: compute_perturbed_purity_exact(
                c_all_by_state[name], delta_c_by_state[name],
                slow_indices, R, delta_R_slow,
                vals, slow_vals, d, N, TIME_WINDOW)
            for name in states
        }
        print(f"  {time.time() - t0:.1f} s")

        print("Extracting α_i and testing closure...")
        bond_results = {"bond": list(bond), "states": {}}
        closure_log_lines.append(f"Bond {bond}:")
        closure_log_lines.append(f"  {'state':<20} {'α_i per site':<80} Σ_i ln(α_i)")
        for name in states:
            PA = PA_by_state[name]
            dP = dP_by_state[name]
            PB = PA + DELTA_J * dP
            alpha = fit_alpha(PA, PB, TIME_WINDOW)
            valid = np.isfinite(alpha) & (alpha > 0)
            sum_ln = float(np.sum(np.log(alpha))) if np.all(valid) else float('nan')
            alpha_str = "[" + ", ".join(f"{a:+.4f}" for a in alpha) + "]"
            closure_log_lines.append(f"  {name:<20} {alpha_str:<80} {sum_ln:+.6e}")
            bond_results["states"][name] = {
                "alpha_per_site": [float(a) for a in alpha],
                "sum_ln_alpha": sum_ln,
                "valid": bool(np.all(valid)),
            }
        results["bonds"].append(bond_results)
        closure_log_lines.append("")
        print()

    # Overall verdict
    closure_log_lines.append("=" * 60)
    closure_log_lines.append("Verdict thresholds:")
    closure_log_lines.append(f"  PASS (numerically verified): |Σ_i ln α_i| < {CLOSURE_TARGET}")
    closure_log_lines.append(f"  PARTIAL (approximate):       |Σ_i ln α_i| < {CLOSURE_PARTIAL}")
    closure_log_lines.append("")

    max_abs_closure = 0.0
    n_valid = 0
    n_invalid = 0
    for bond_res in results["bonds"]:
        for state_name, state_res in bond_res["states"].items():
            s = state_res["sum_ln_alpha"]
            if np.isfinite(s):
                max_abs_closure = max(max_abs_closure, abs(s))
                n_valid += 1
            else:
                n_invalid += 1

    closure_log_lines.append(f"Valid (finite α_i all sites): {n_valid}")
    closure_log_lines.append(f"Invalid (NaN in α_i):         {n_invalid}")

    if n_valid == 0:
        verdict = "INVALID (no states produced finite α_i; likely a bug)"
    elif max_abs_closure < CLOSURE_TARGET:
        verdict = f"PASS (closure law numerically verified on {n_valid}/{n_valid + n_invalid} states)"
    elif max_abs_closure < CLOSURE_PARTIAL:
        verdict = f"PARTIAL (approximate closure on {n_valid}/{n_valid + n_invalid} states)"
    else:
        verdict = f"FAIL (closure law not verified at this precision; {n_valid}/{n_valid + n_invalid} states valid)"
    closure_log_lines.append(f"Max |Σ_i ln α_i|: {max_abs_closure:.3e}")
    closure_log_lines.append(f"Overall: {verdict}")
    results["verdict"] = verdict
    results["max_abs_closure"] = float(max_abs_closure)

    with open(RESULTS_DIR / "eq014_alpha_prediction.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    print("Wrote eq014_alpha_prediction.json")

    with open(RESULTS_DIR / "eq014_closure_test.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(closure_log_lines) + "\n")
    print("Wrote eq014_closure_test.txt")

    print()
    print(f"=== Steps 4-7 complete in {time.time() - t_total:.1f} s ===")
    print(f"Max |Σ_i ln α_i|: {max_abs_closure:.3e}")
    print(f"Verdict: {verdict}")


if __name__ == "__main__":
    main()
