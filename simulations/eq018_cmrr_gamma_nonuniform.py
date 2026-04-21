#!/usr/bin/env python3
"""eq018_cmrr_gamma_nonuniform.py

Prediction (d) test from ORTHOGONALITY_SELECTION_FAMILY.md:

    Non-uniform gamma_i breaks the CMRR (Common-Mode Rejection Ratio) of
    the spatial-sum purity detector, making K_CC[0,1]_pr finite instead
    of exactly zero.

Meta-theorem claim:
    Under uniform gamma_0, the (vac, S_1) coherence has decay rate
    2*gamma_0 uniformly across all |1_i> basis states (d_H=1 to vacuum).
    Sine-basis orthogonality + this uniformity force
        sum_i 2|rho_{coh,i,01}(t)|^2 = (1/2) * exp(-4*gamma_0*t)
    Hamiltonian-independent, hence K_CC[0,1]_pr = 0 exactly.

    Break the uniform-gamma assumption:
        decay rate for |vac><1_i| becomes 2*gamma_i (site-dependent)
    The completeness-collapse argument fails; the sum no longer has a
    universal exp(-4*gamma t) form; delta_dJ of the sum is non-zero.

Setup: N=5 chain, bond 0 perturbed, gamma profiles:
    1. uniform gamma_0 = 0.05 (baseline; expected K_CC[0,1]_pr = 0)
    2. single-site bump: gamma_0 at site 0 = 0.05 + delta_g; rest = 0.05
    3. linear gradient: gamma_i = 0.05 + alpha * (i - (N-1)/2)
    4. random profile (fixed seed)

For profile 2, scan delta_g in {0, 0.01, 0.05, 0.1, 0.2}. Measure
K_CC[0,1]_pr at t_0 = 20. Expected: grows with delta_g.

Rules: UTF-8 stdout, no em-dashes, hyphens only.
"""
from __future__ import annotations

import json
import sys
import time
from math import comb
from pathlib import Path

import numpy as np
from scipy.linalg import eig

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).parent))
from pi_pair_closure_investigation import (
    GAMMA_0, J_UNIFORM, T_FINAL, N_STEPS,
    build_H_XY, X, Y, Z, I2, site_op, kron_chain,
    vacuum_ket, density_matrix,
    per_site_purity,
)
from c1_bilinearity_test import dicke_state

RESULTS_DIR = Path(__file__).parent / "results" / "eq018_cmrr_gamma_nonuniform"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DJ_EXTRACT = 0.01


def build_liouvillian_nonuniform(H, gamma_list, N):
    """L[rho] = -i[H,rho] + sum_i gamma_i (Z_i rho Z_i - rho)

    gamma_list: length-N array of per-site dephasing rates.
    """
    d = 2**N
    I_d = np.eye(d, dtype=complex)
    L = -1j * (np.kron(I_d, H) - np.kron(H.T, I_d))
    for i in range(N):
        Zi = site_op(Z, i, N)
        L += gamma_list[i] * (np.kron(Zi.T, Zi) - np.kron(I_d, I_d))
    return L


def eig_and_inv(L):
    eigvals, V_R = eig(L)
    V_Linv = np.linalg.inv(V_R)
    return eigvals, V_R, V_Linv


def propagate(eigvals, V_R, V_Linv, rho_0, times):
    d = rho_0.shape[0]
    rho0_vec = rho_0.flatten(order='F')
    c0 = V_Linv @ rho0_vec
    out = np.empty((len(times), d, d), dtype=complex)
    for k, t in enumerate(times):
        rho_vec_t = V_R @ (np.exp(eigvals * t) * c0)
        out[k] = rho_vec_t.reshape(d, d, order='F')
    return out


def build_decomps_nonuniform(N, bond, gamma_list, dj=DJ_EXTRACT, J=J_UNIFORM):
    """Return (decomp_A, decomp_B+, decomp_B-) using site-dependent gamma."""
    J_A = [J] * (N - 1)
    J_Bp = list(J_A); J_Bp[bond] = J + dj
    J_Bm = list(J_A); J_Bm[bond] = J - dj
    L_A = build_liouvillian_nonuniform(build_H_XY(J_A, N), gamma_list, N)
    L_Bp = build_liouvillian_nonuniform(build_H_XY(J_Bp, N), gamma_list, N)
    L_Bm = build_liouvillian_nonuniform(build_H_XY(J_Bm, N), gamma_list, N)
    return eig_and_inv(L_A), eig_and_inv(L_Bp), eig_and_inv(L_Bm)


def measure_c1_pr(rho_0, decomps, times, t_0, N, dj=DJ_EXTRACT):
    """Purity-response c_1 at time t_0."""
    (evA, VRA, VLA), (evBp, VRBp, VLBp), (evBm, VRBm, VLBm) = decomps
    rho_Bp = propagate(evBp, VRBp, VLBp, rho_0, times)
    rho_Bm = propagate(evBm, VRBm, VLBm, rho_0, times)
    P_Bp = per_site_purity(rho_Bp, N)
    P_Bm = per_site_purity(rho_Bm, N)
    t_0_idx = int(np.argmin(np.abs(times - t_0)))
    dP_B_dJ = (P_Bp[t_0_idx] - P_Bm[t_0_idx]) / (2 * dj)
    return float(np.sum(dP_B_dJ))


def measure_K_CC_01_pr(N, bond, decomps, times, t_0):
    """Extract K_CC[0, 1]_pr from coh vs mix probes."""
    S0 = dicke_state(N, 0)
    S1 = dicke_state(N, 1)
    rho_mix = 0.5 * (density_matrix(S0) + density_matrix(S1))
    rho_coh = density_matrix((S0 + S1) / np.sqrt(2.0))
    c_mix = measure_c1_pr(rho_mix, decomps, times, t_0, N)
    c_coh = measure_c1_pr(rho_coh, decomps, times, t_0, N)
    return 2.0 * (c_coh - c_mix), c_mix, c_coh


def measure_sum_coh_purity_time_series(rho_0_coh_op, decomps_A, times, N):
    """Measure sum_i 2|rho_coh_i(t)_01|^2 as a time series to verify
    whether it's still (1/2) exp(-4 <gamma> t) (uniform case) or
    H-dependent (non-uniform case)."""
    evA, VRA, VLA = decomps_A
    rho_A = propagate(evA, VRA, VLA, rho_0_coh_op, times)
    from pi_pair_closure_investigation import partial_trace_keep_site_fast
    T_len = rho_A.shape[0]
    S = np.zeros(T_len)
    for t_idx in range(T_len):
        acc = 0.0
        for i in range(N):
            ri = partial_trace_keep_site_fast(rho_A[t_idx], i, N)
            acc += 2 * abs(ri[0, 1])**2
        S[t_idx] = acc
    return S


def run_profile(N, bond, t_0, gamma_list, label):
    """Run one gamma profile and report K_CC[0,1]_pr."""
    times = np.linspace(0.0, T_FINAL, N_STEPS + 1)
    decomps = build_decomps_nonuniform(N, bond, gamma_list)
    K_CC, c_mix, c_coh = measure_K_CC_01_pr(N, bond, decomps, times, t_0)
    # Also measure the time series of the purity-sum to visualize the break
    S0 = dicke_state(N, 0)
    S1 = dicke_state(N, 1)
    rho_coh_only = 0.5 * (np.outer(S0, S1.conj()) + np.outer(S1, S0.conj()))
    S_t = measure_sum_coh_purity_time_series(rho_coh_only, decomps[0], times, N)
    # Gamma mean
    gamma_mean = float(np.mean(gamma_list))
    # Expected uniform-case value at t_0
    t_0_idx = int(np.argmin(np.abs(times - t_0)))
    S_at_t0 = float(S_t[t_0_idx])
    S_uniform_pred = 0.5 * np.exp(-4 * gamma_mean * t_0)
    return {
        "label": label,
        "gamma_list": list(gamma_list),
        "gamma_mean": gamma_mean,
        "gamma_variance": float(np.var(gamma_list)),
        "K_CC_01_pr": K_CC,
        "c_1_mix_pr": c_mix,
        "c_1_coh_pr": c_coh,
        "S_coh_sum_at_t0": S_at_t0,
        "S_uniform_pred": float(S_uniform_pred),
        "S_deviation_from_uniform": S_at_t0 - float(S_uniform_pred),
        "S_time_series": [float(x) for x in S_t],
    }


def main():
    start = time.time()
    N = 5
    bond = 0
    t_0 = 1.0 / GAMMA_0  # = 20.0
    print("=" * 78)
    print(f"EQ-018 Prediction (d): CMRR break under non-uniform gamma_i at N={N}, bond={bond}")
    print("=" * 78)
    print(f"  Baseline gamma_0 = {GAMMA_0}, J = {J_UNIFORM}, dJ = +/- {DJ_EXTRACT}, t_0 = {t_0}")

    # Profile 1: uniform (baseline)
    print(f"\n--- Profile 1: UNIFORM gamma = {GAMMA_0} (baseline) ---")
    gamma_unif = [GAMMA_0] * N
    r1 = run_profile(N, bond, t_0, gamma_unif, "uniform baseline")
    print(f"  K_CC[0,1]_pr = {r1['K_CC_01_pr']:+.5e}  (should be ~0)")
    print(f"  Sum_i 2|rho_coh_i(t_0)_01|^2 = {r1['S_coh_sum_at_t0']:.6e}")
    print(f"  Predicted (1/2)*exp(-4*gamma_mean*t_0) = {r1['S_uniform_pred']:.6e}")
    print(f"  Deviation from (1/2)*exp(-4*gamma*t_0): {r1['S_deviation_from_uniform']:.2e}")

    # Profile 2: single-site bump at site 0 with varying delta
    print(f"\n--- Profile 2: SINGLE-SITE BUMP at site 0, scan delta_g ---")
    delta_g_values = [0.0, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2]
    bump_results = []
    for dg in delta_g_values:
        gamma_prof = [GAMMA_0] * N
        gamma_prof[0] = GAMMA_0 + dg
        r = run_profile(N, bond, t_0, gamma_prof, f"bump site0 +{dg}")
        bump_results.append(r)
        print(f"  delta_g = {dg:>6.3f}: K_CC[0,1]_pr = {r['K_CC_01_pr']:+.5e}"
              f"  Var(gamma) = {r['gamma_variance']:.3e}")

    # Profile 3: linear gradient across the chain
    print(f"\n--- Profile 3: LINEAR GRADIENT, scan slope alpha ---")
    alpha_values = [0.0, 0.002, 0.005, 0.01, 0.02]
    grad_results = []
    for alpha in alpha_values:
        gamma_prof = [GAMMA_0 + alpha * (i - (N - 1) / 2) for i in range(N)]
        r = run_profile(N, bond, t_0, gamma_prof, f"gradient alpha={alpha}")
        grad_results.append(r)
        print(f"  alpha = {alpha:>6.3f}: gamma = {[f'{g:.3f}' for g in gamma_prof]}"
              f"  K_CC[0,1]_pr = {r['K_CC_01_pr']:+.5e}")

    # Profile 4: random profile (reproducible)
    print(f"\n--- Profile 4: RANDOM gamma profile (seed=42) ---")
    rng = np.random.default_rng(42)
    random_results = []
    for trial in range(3):
        perturb = rng.uniform(-0.02, 0.02, N)
        gamma_prof = [GAMMA_0 + float(p) for p in perturb]
        r = run_profile(N, bond, t_0, gamma_prof, f"random trial {trial}")
        random_results.append(r)
        print(f"  trial {trial}: gamma = {[f'{g:.3f}' for g in gamma_prof]}"
              f"  Var(gamma) = {r['gamma_variance']:.3e}"
              f"  K_CC[0,1]_pr = {r['K_CC_01_pr']:+.5e}")

    # Analytical comparison: for small delta_g, is K_CC linear or quadratic?
    print(f"\n--- Scaling analysis: K_CC[0,1]_pr vs gamma variance ---")
    print(f"  {'profile':>25} {'Var(g)':>12} {'K_CC':>14} {'K_CC/Var(g)':>14}")
    for r in bump_results:
        v = r["gamma_variance"]
        K = r["K_CC_01_pr"]
        ratio = K / v if v > 1e-12 else float('nan')
        print(f"  {r['label']:>25} {v:>12.3e} {K:>+14.5e} {ratio:>+14.5e}")

    # Save
    out = {
        "N": N, "bond": bond, "t_0": t_0,
        "gamma_0_baseline": GAMMA_0,
        "J": J_UNIFORM, "dJ": DJ_EXTRACT,
        "profile_1_uniform": r1,
        "profile_2_bump_scan": bump_results,
        "profile_3_gradient_scan": grad_results,
        "profile_4_random": random_results,
    }
    # Drop time series for compact JSON
    def _strip_ts(d):
        return {k: (v if k != "S_time_series" else "omitted_in_json") for k, v in d.items()}
    out_compact = {
        "N": N, "bond": bond, "t_0": t_0,
        "gamma_0_baseline": GAMMA_0,
        "J": J_UNIFORM, "dJ": DJ_EXTRACT,
        "profile_1_uniform": _strip_ts(r1),
        "profile_2_bump_scan": [_strip_ts(r) for r in bump_results],
        "profile_3_gradient_scan": [_strip_ts(r) for r in grad_results],
        "profile_4_random": [_strip_ts(r) for r in random_results],
    }
    path = RESULTS_DIR / "cmrr_gamma_nonuniform.json"
    with open(path, "w") as f:
        json.dump(out_compact, f, indent=2, default=str)
    # Save time series separately
    ts_path = RESULTS_DIR / "coh_purity_time_series.json"
    ts = {
        "times": list(np.linspace(0.0, T_FINAL, N_STEPS + 1)),
        "uniform_baseline": r1["S_time_series"],
        "bump_delta_0.1": next((r["S_time_series"] for r in bump_results
                                 if abs(r["gamma_list"][0] - (GAMMA_0 + 0.1)) < 1e-9), None),
        "bump_delta_0.2": next((r["S_time_series"] for r in bump_results
                                 if abs(r["gamma_list"][0] - (GAMMA_0 + 0.2)) < 1e-9), None),
    }
    with open(ts_path, "w") as f:
        json.dump(ts, f, indent=2, default=str)
    print(f"\nSaved: {path}")
    print(f"Saved: {ts_path}")
    print(f"\nTotal walltime: {time.time() - start:.1f} s")

    # Summary verdict
    print(f"\n{'=' * 78}")
    print(f"VERDICT")
    print(f"{'=' * 78}")
    K_baseline = r1["K_CC_01_pr"]
    K_strongest_bump = bump_results[-1]["K_CC_01_pr"]
    K_ratio = K_strongest_bump / (K_baseline + 1e-15) if abs(K_baseline) > 1e-15 else float('inf')
    print(f"  Baseline K_CC[0,1]_pr (uniform)       = {K_baseline:+.3e}")
    print(f"  Strongest non-uniform (delta=0.2)    = {K_strongest_bump:+.3e}")
    print(f"  Ratio                                  = {K_ratio:.3g}")
    print(f"  If K grows with gamma variance, meta-theorem prediction (d) CONFIRMED.")
    print(f"  If K stays at ~1e-12 regardless, prediction FALSIFIED.")


if __name__ == "__main__":
    main()
