#!/usr/bin/env python3
"""c1_sector_kernel.py

Isolate the sector-coherence contribution to c_1 by comparing mixed
vs pure superposition initial states.

For two states |A> and |B>:
  rho_mix  = (1/2)(|A><A| + |B><B|)          classical mixture, no coherence
  rho_coh  = (1/2)(|A>+|B>)(<A| + <B|)       coherent superposition, has coherence

Difference rho_coh - rho_mix = (1/2)(|A><B| + |B><A|)
                            = pure off-diagonal (A, B) block

Bilinearity: c_1(coh) = c_1(mix) + (coherence-only contribution).

The coherence-only contribution is the bilinear form's value on the
off-diagonal block paired with itself and with the diagonal blocks.
Testing across Dicke-state pairs (S_n, S_m) for various |n - m| gives
the sector-kernel K as a function of Delta N and sector position.

N=5, same setup. ~2 minutes compute.
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
    GAMMA_0, J_UNIFORM, T_FINAL, N_STEPS, T_FIT_MAX,
    build_H_XY, build_liouvillian_matrix,
    vacuum_ket, single_excitation_mode, density_matrix,
    per_site_purity, fit_alpha,
)
from c1_bilinearity_test import dicke_state

RESULTS_DIR = Path(__file__).parent / "results" / "c1_sector_kernel"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

N = 5
DJ_EXTRACT = 0.01


def mixed_density(ket_A, ket_B):
    """Classical mixture (|A><A| + |B><B|) / 2."""
    return 0.5 * (density_matrix(ket_A) + density_matrix(ket_B))


def coherent_density(ket_A, ket_B):
    """Pure superposition state density matrix from (|A> + |B>)/sqrt(2)."""
    v = ket_A + ket_B
    v = v / np.linalg.norm(v)
    return density_matrix(v)


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


def measure_c1(rho_0, decomps, times):
    """Measure c_1 for a given initial density matrix."""
    (evA, VRA, VLA), (evBp, VRBp, VLBp), (evBm, VRBm, VLBm) = decomps
    rho_A = propagate(evA, VRA, VLA, rho_0, times)
    rho_Bp = propagate(evBp, VRBp, VLBp, rho_0, times)
    rho_Bm = propagate(evBm, VRBm, VLBm, rho_0, times)
    P_A = per_site_purity(rho_A, N)
    P_Bp = per_site_purity(rho_Bp, N)
    P_Bm = per_site_purity(rho_Bm, N)
    alpha_p = np.zeros(N); alpha_m = np.zeros(N)
    rmse_p = 0.0
    for i in range(N):
        a, r = fit_alpha(times, P_A[:, i], P_Bp[:, i])
        alpha_p[i] = a; rmse_p = max(rmse_p, r)
        a, _ = fit_alpha(times, P_A[:, i], P_Bm[:, i])
        alpha_m[i] = a
    closure_p = float(np.sum(np.log(alpha_p)))
    closure_m = float(np.sum(np.log(alpha_m)))
    c_1 = (closure_p - closure_m) / (2 * DJ_EXTRACT)
    return {"c_1": c_1, "rmse": rmse_p,
            "alpha_plus": alpha_p.tolist(),
            "alpha_minus": alpha_m.tolist()}


def main():
    print("=" * 70)
    print(f"Sector-kernel extraction at N = {N}")
    print("=" * 70)
    print(f"  gamma_0 = {GAMMA_0}, J = {J_UNIFORM}, dJ = +/- {DJ_EXTRACT}\n")

    # Build Liouvillians
    J_A = [J_UNIFORM] * (N - 1)
    J_B_plus = list(J_A); J_B_plus[0] = J_UNIFORM + DJ_EXTRACT
    J_B_minus = list(J_A); J_B_minus[0] = J_UNIFORM - DJ_EXTRACT

    print("  Building Liouvillians...")
    t0 = time.time()
    L_A = build_liouvillian_matrix(build_H_XY(J_A, N), GAMMA_0, N)
    decA = eig_and_inv(L_A); del L_A
    L_Bp = build_liouvillian_matrix(build_H_XY(J_B_plus, N), GAMMA_0, N)
    decBp = eig_and_inv(L_Bp); del L_Bp
    L_Bm = build_liouvillian_matrix(build_H_XY(J_B_minus, N), GAMMA_0, N)
    decBm = eig_and_inv(L_Bm); del L_Bm
    decomps = (decA, decBp, decBm)
    print(f"  Built in {time.time()-t0:.1f} s")

    times_arr = np.linspace(0.0, T_FINAL, N_STEPS + 1)

    # Dicke states and vacuum
    S = {n: dicke_state(N, n) for n in range(N + 1)}

    # Pure c_1 for each Dicke state
    print(f"\n  Pure Dicke-state c_1 (diagonal-only initial rho):")
    c1_pure = {}
    for n in range(N + 1):
        rho_0 = density_matrix(S[n])
        r = measure_c1(rho_0, decomps, times_arr)
        c1_pure[n] = r["c_1"]
        print(f"    |S_{n}>: c_1 = {r['c_1']:+.5f}")

    # For each pair (n, m) with n < m, compute c_1(mixed) and c_1(coherent)
    # and isolate coherence contribution
    print(f"\n  (n, m) pairs: c_1(mixed), c_1(coh), coherence contribution")
    print(f"  {'n':>2} {'m':>2} {'Delta N':>8} {'c_1(mix)':>12} "
          f"{'c_1(coh)':>12} {'coh_contrib':>14}")
    kernel_data = {}
    for n in range(N + 1):
        for m in range(n + 1, N + 1):
            dN = m - n
            rho_mix = mixed_density(S[n], S[m])
            rho_coh = coherent_density(S[n], S[m])
            r_mix = measure_c1(rho_mix, decomps, times_arr)
            r_coh = measure_c1(rho_coh, decomps, times_arr)
            coh_contrib = r_coh["c_1"] - r_mix["c_1"]
            print(f"  {n:>2d} {m:>2d} {dN:>8d} "
                  f"{r_mix['c_1']:>+12.5f} "
                  f"{r_coh['c_1']:>+12.5f} "
                  f"{coh_contrib:>+14.5f}")
            kernel_data[f"{n}_{m}"] = {
                "n": n, "m": m, "delta_N": dN,
                "c_1_mixed": r_mix["c_1"],
                "c_1_coherent": r_coh["c_1"],
                "coherence_contribution": coh_contrib,
            }

    # Group by Delta N and report averages
    print(f"\n  Coherence contribution grouped by |Delta N|:")
    print(f"  {'Delta N':>8} {'n,m pairs':>30} {'coh_contrib values':>30}")
    by_dN = {}
    for key, d in kernel_data.items():
        dN = d["delta_N"]
        by_dN.setdefault(dN, []).append(d)
    for dN in sorted(by_dN.keys()):
        entries = by_dN[dN]
        pair_labels = ", ".join(f"({e['n']},{e['m']})" for e in entries)
        values = [e["coherence_contribution"] for e in entries]
        values_str = ", ".join(f"{v:+.3f}" for v in values)
        mean = sum(values) / len(values)
        print(f"  {dN:>8d} {pair_labels:>30} {values_str:>30}")
        print(f"           mean: {mean:+.4f}")

    # Mixed-state diagonal test: c_1(mix) = sum/combinations of diagonals?
    # For rho_mix = (|S_n><S_n| + |S_m><S_m|)/2:
    # Prediction if c_1 purely diagonal-bilinear: c_1(mix) = (1/4)[c_1(S_n) + c_1(S_m)] + 2*(cross-diag-bilinear)
    print(f"\n  Mixed-state diagonal-bilinear decomposition check:")
    print(f"  (n,m)  c_1(mix)  pred_0.25*(c_n+c_m)  cross-diag-bilinear")
    for key, d in kernel_data.items():
        n, m = d["n"], d["m"]
        pred = 0.25 * (c1_pure[n] + c1_pure[m])
        cross_diag = d["c_1_mixed"] - pred
        print(f"  ({n},{m})  {d['c_1_mixed']:+.5f}  {pred:+.5f}  {cross_diag:+.5f}")

    # Save
    out = {
        "N": N, "gamma_0": GAMMA_0, "J": J_UNIFORM, "dJ": DJ_EXTRACT,
        "defect_bond": [0, 1],
        "c1_pure_Sn": c1_pure,
        "pair_results": kernel_data,
    }
    path = RESULTS_DIR / "sector_kernel.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nSaved: {path}")


if __name__ == "__main__":
    main()
