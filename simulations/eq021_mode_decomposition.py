#!/usr/bin/env python3
"""eq021_mode_decomposition.py

EQ-021 Phase 1: Decompose c_1 into per-mode contributions.

c_1 is the first-order coefficient of the PTF closure sum
Sigma_i ln(alpha_i) under a bond perturbation delta_J. Empirically
c_1 ~ 0.5 * V(N) for the initial state (|vac> + |psi_1>)/sqrt(2) at
bond (0,1), for N >= 4. The 0.5 factor has no derivation yet.

This script diagonalises the Liouvillian L_A (uniform XY chain + Z
dephasing) at small N, projects the initial state onto L_A eigenmodes,
and expresses P_i(t) = Tr(rho_i^2) as a sum of mode-pair contributions:

    P_i(t) = Sum_{s, s'} C_{s, s', i} * exp((lambda_s + lambda_{s'}*) t)

where C_{s, s', i} = c_s * c_{s'}^* * Tr(M_{s, i} * M_{s', i}^dag)
and M_{s, i} = Tr_{not i}(M_s).

Under perturbation the mode-pair contributions shift; the shift's
decomposition into C, lambda, and partial-trace corrections tells
us which modes drive c_1.

At N = 4, 5, 6 the Liouvillian is small enough that we can enumerate
all mode pairs and identify the dominant contributions.

Rules from task:
  - XY Hamiltonian H = (J/2)(XX+YY), NOT Heisenberg. Critical.
  - Numbers from script output only.
  - UTF-8 stdout.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eig

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).parent))
from pi_pair_closure_investigation import (
    GAMMA_0, J_UNIFORM, T_FINAL, N_STEPS,
    build_H_XY, build_liouvillian_matrix,
    vacuum_ket, single_excitation_mode, density_matrix,
    per_site_purity, fit_alpha, partial_trace_keep_site_fast,
)

RESULTS_DIR = Path(__file__).parent / "results" / "eq021_mode_decomposition"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DJ_EXTRACT = 0.01
T_FIT_MAX = 20.0


def v_effect(N):
    return 1.0 + np.cos(np.pi / N)


def bonding_plus_vacuum(N, k=1):
    v = vacuum_ket(N)
    psi = single_excitation_mode(N, k=k)
    return (v + psi) / np.sqrt(2.0)


def eig_biorthogonal(L):
    """Return (eigvals, V_R, W_dag) with W_dag @ V_R = I (biorthonormal).

    Here V_R columns are right eigenvectors and W_dag rows are left
    eigenvectors, both normalised so W_dag V_R = I.
    """
    eigvals, V_R = eig(L)
    W_dag = np.linalg.inv(V_R)  # rows of W_dag are left eigenvectors
    return eigvals, V_R, W_dag


def compute_mode_partial_traces(V_R, N):
    """For each mode (column of V_R reshaped to d x d), compute its
    partial trace to each of N sites. Returns array of shape
    (d^2, N, 2, 2) giving M_{s, i} as a 2x2 matrix."""
    d = 2**N
    n_modes = V_R.shape[1]
    out = np.zeros((n_modes, N, 2, 2), dtype=complex)
    for s in range(n_modes):
        M_s = V_R[:, s].reshape(d, d, order='F')
        for i in range(N):
            out[s, i] = partial_trace_keep_site_fast(M_s, i, N)
    return out


def compute_P_i_from_modes(c_s, eigvals, mode_margs, times, i):
    """P_i(t) = Sum_{s, s'} c_s c_{s'}^* exp((lambda_s + lambda_{s'}*) t)
                 * Tr(M_{s, i} @ M_{s', i}^dagger)"""
    n_modes = len(c_s)
    T = len(times)
    out = np.zeros(T, dtype=complex)
    # Precompute Tr(M_s M_s'^dag) matrix
    trace_mat = np.zeros((n_modes, n_modes), dtype=complex)
    for s in range(n_modes):
        for sp in range(n_modes):
            trace_mat[s, sp] = np.trace(mode_margs[s, i] @ mode_margs[sp, i].conj().T)
    # Coefficient matrix c_s c_{s'}^*
    coeff = np.outer(c_s, np.conj(c_s))
    # Time evolution: exp((lambda_s + lambda_{s'}*) t) outer sum
    # omega[s, s'] = lambda_s + conj(lambda_{s'})
    omega = eigvals[:, None] + np.conj(eigvals[None, :])
    # Summation over (s, s')
    for k_t, t in enumerate(times):
        weights = coeff * np.exp(omega * t)
        out[k_t] = np.sum(weights * trace_mat)
    return out.real


def compute_all_per_site_P(c_s, eigvals, mode_margs, times, N):
    out = np.zeros((len(times), N))
    for i in range(N):
        out[:, i] = compute_P_i_from_modes(c_s, eigvals, mode_margs, times, i)
    return out


def project_rho0_onto_modes(rho_0, W_dag):
    """c_s = <W_s|rho_0> = (W_dag @ vec(rho_0))_s"""
    return W_dag @ rho_0.flatten(order='F')


def dominant_mode_pairs(c_s, eigvals, mode_margs, i, top_k=10):
    """Identify the top_k mode pairs (s, s') contributing to P_i at
    t = 0 (initial weight). Returns list of (contribution, s, s') sorted
    by abs contribution descending."""
    n_modes = len(c_s)
    pairs = []
    for s in range(n_modes):
        for sp in range(n_modes):
            contrib = (c_s[s] * np.conj(c_s[sp]) *
                       np.trace(mode_margs[s, i] @ mode_margs[sp, i].conj().T))
            if abs(contrib) > 1e-10:
                pairs.append((abs(contrib), contrib, s, sp,
                              eigvals[s], eigvals[sp]))
    pairs.sort(reverse=True)
    return pairs[:top_k]


def run_N(N):
    print(f"\n{'=' * 70}")
    print(f"N = {N}")
    print(f"{'=' * 70}")
    d = 2**N
    V_N = v_effect(N)
    print(f"  V(N) = {V_N:.4f}, 0.5 * V(N) = {0.5 * V_N:.4f}")

    # Build L_A and L_B with small delta_J
    J_A = [J_UNIFORM] * (N - 1)
    J_Bp = list(J_A); J_Bp[0] = J_UNIFORM + DJ_EXTRACT
    J_Bm = list(J_A); J_Bm[0] = J_UNIFORM - DJ_EXTRACT

    t0 = time.time()
    L_A = build_liouvillian_matrix(build_H_XY(J_A, N), GAMMA_0, N)
    ev_A, V_R_A, W_dag_A = eig_biorthogonal(L_A)
    print(f"  L_A diagonalised in {time.time() - t0:.2f} s (d^2 = {d**2})")

    t0 = time.time()
    L_Bp = build_liouvillian_matrix(build_H_XY(J_Bp, N), GAMMA_0, N)
    ev_Bp, V_R_Bp, W_dag_Bp = eig_biorthogonal(L_Bp)
    L_Bm = build_liouvillian_matrix(build_H_XY(J_Bm, N), GAMMA_0, N)
    ev_Bm, V_R_Bm, W_dag_Bm = eig_biorthogonal(L_Bm)
    print(f"  L_B+/- diagonalised in {time.time() - t0:.2f} s")

    # Initial state and its mode-projection
    phi = bonding_plus_vacuum(N, k=1)
    rho_0 = density_matrix(phi)
    c_s_A = project_rho0_onto_modes(rho_0, W_dag_A)
    print(f"  |c_s_A|^2 range: min={np.abs(c_s_A**2).min():.2e}, "
          f"max={np.abs(c_s_A**2).max():.4f}")
    print(f"  sum |c_s_A| (weight): {np.sum(np.abs(c_s_A)):.4f}")
    n_significant_A = int(np.sum(np.abs(c_s_A) > 1e-8))
    print(f"  Modes with |c_s| > 1e-8: {n_significant_A} (out of {d**2})")

    # Compute partial-trace marginals for each mode at each site
    t0 = time.time()
    mode_margs_A = compute_mode_partial_traces(V_R_A, N)
    print(f"  Mode marginals (A): {time.time() - t0:.2f} s")

    # Direct time-domain P_i(t) from mode expansion
    times = np.linspace(0, T_FINAL, N_STEPS + 1)
    t0 = time.time()
    P_A_modes = compute_all_per_site_P(c_s_A, ev_A, mode_margs_A, times, N)
    print(f"  P_A(i, t) from modes: {time.time() - t0:.2f} s")

    # Cross-check with direct propagation
    rho0_vec = rho_0.flatten(order='F')
    c0 = W_dag_A @ rho0_vec
    P_A_direct = np.zeros((len(times), N))
    for k_t, t in enumerate(times):
        rho_vec_t = V_R_A @ (np.exp(ev_A * t) * c0)
        rho_t = rho_vec_t.reshape(d, d, order='F')
        for i in range(N):
            rho_i = partial_trace_keep_site_fast(rho_t, i, N)
            P_A_direct[k_t, i] = float(np.trace(rho_i @ rho_i).real)
    direct_vs_modes_diff = np.max(np.abs(P_A_modes - P_A_direct))
    print(f"  max |P_modes - P_direct|: {direct_vs_modes_diff:.2e}")
    if direct_vs_modes_diff > 1e-8:
        print(f"  WARNING: mode expansion and direct propagation disagree!")

    # Same for L_B+
    c_s_Bp = project_rho0_onto_modes(rho_0, W_dag_Bp)
    c_s_Bm = project_rho0_onto_modes(rho_0, W_dag_Bm)
    mode_margs_Bp = compute_mode_partial_traces(V_R_Bp, N)
    mode_margs_Bm = compute_mode_partial_traces(V_R_Bm, N)
    P_Bp_modes = compute_all_per_site_P(c_s_Bp, ev_Bp, mode_margs_Bp, times, N)
    P_Bm_modes = compute_all_per_site_P(c_s_Bm, ev_Bm, mode_margs_Bm, times, N)

    # Fit alpha_i and closure
    alpha_p = np.zeros(N); alpha_m = np.zeros(N)
    for i in range(N):
        a, _ = fit_alpha(times, P_A_modes[:, i], P_Bp_modes[:, i])
        alpha_p[i] = a
        a, _ = fit_alpha(times, P_A_modes[:, i], P_Bm_modes[:, i])
        alpha_m[i] = a
    closure_p = float(np.sum(np.log(alpha_p)))
    closure_m = float(np.sum(np.log(alpha_m)))
    c_1 = (closure_p - closure_m) / (2 * DJ_EXTRACT)
    ratio_to_05V = c_1 / (0.5 * V_N)
    print(f"\n  ALPHA FIT: c_1 = {c_1:+.5f}")
    print(f"  0.5 * V(N) = {0.5 * V_N:.5f}, ratio c_1 / (0.5 V(N)) = {ratio_to_05V:.4f}")

    # Per-mode-pair contribution to P_i at several key times
    # Focus on one representative site (middle and endpoint)
    test_sites = [0, N // 2]
    test_times_idx = [0, len(times) // 4, len(times) // 2, len(times) - 1]

    print(f"\n  Top mode pairs contributing to P_i(t=0) at site 0:")
    top_pairs_0 = dominant_mode_pairs(c_s_A, ev_A, mode_margs_A, i=0, top_k=10)
    for rank, (_abs, contrib, s, sp, lam_s, lam_sp) in enumerate(top_pairs_0):
        print(f"    #{rank+1}: contrib={contrib.real:+.5f} + {contrib.imag:+.5f}i | "
              f"s={s} (lam={lam_s.real:+.3f}{lam_s.imag:+.3f}i), "
              f"s'={sp} (lam={lam_sp.real:+.3f}{lam_sp.imag:+.3f}i)")

    # Save
    out = {
        "N": N, "gamma_0": GAMMA_0, "J": J_UNIFORM, "dJ": DJ_EXTRACT,
        "d_squared": d**2, "V_N": V_N,
        "c_1_alpha_fit": c_1,
        "c_1_over_05V": ratio_to_05V,
        "closure_plus": closure_p,
        "closure_minus": closure_m,
        "alpha_plus": alpha_p.tolist(),
        "alpha_minus": alpha_m.tolist(),
        "n_significant_modes": n_significant_A,
        "top_mode_pairs_site0": [
            {"rank": rank + 1, "s": int(s), "sp": int(sp),
             "abs_contrib": float(abs_),
             "contrib_real": float(contrib.real),
             "contrib_imag": float(contrib.imag),
             "lambda_s_real": float(lam_s.real),
             "lambda_s_imag": float(lam_s.imag),
             "lambda_sp_real": float(lam_sp.real),
             "lambda_sp_imag": float(lam_sp.imag)}
            for rank, (abs_, contrib, s, sp, lam_s, lam_sp) in enumerate(top_pairs_0)
        ],
    }
    path = RESULTS_DIR / f"mode_contributions_N{N}.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"  Saved: {path}")
    return out


def main():
    print("=" * 70)
    print("EQ-021 Phase 1: Liouvillian mode decomposition of c_1")
    print("=" * 70)
    print(f"  gamma_0 = {GAMMA_0}, J = {J_UNIFORM}, dJ = +/- {DJ_EXTRACT}")
    print(f"  Initial state: (|vac> + |psi_1>)/sqrt(2), bond (0,1) perturbation")

    results = {}
    for N in [4, 5]:
        results[f"N_{N}"] = run_N(N)

    summary_path = RESULTS_DIR / "decomposition_summary.txt"
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("EQ-021 Phase 1 Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"For (|vac> + |psi_1>)/sqrt(2) at bond (0,1) with delta_J = "
                f"+/- {DJ_EXTRACT}:\n\n")
        f.write(f"{'N':>3} {'V(N)':>8} {'0.5*V(N)':>10} {'c_1':>10} "
                f"{'ratio':>10}\n")
        for N in [4, 5]:
            r = results[f"N_{N}"]
            f.write(f"{N:>3} {r['V_N']:>8.4f} {0.5 * r['V_N']:>10.4f} "
                    f"{r['c_1_alpha_fit']:>+10.4f} "
                    f"{r['c_1_over_05V']:>10.4f}\n")
        f.write("\nPer-site mode-decomposition details in "
                "mode_contributions_N*.json.\n")
    print(f"\nSaved: {summary_path}")


if __name__ == "__main__":
    main()
