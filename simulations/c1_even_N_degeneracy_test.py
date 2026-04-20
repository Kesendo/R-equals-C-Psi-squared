#!/usr/bin/env python3
"""c1_even_N_degeneracy_test.py

Falsification test for the "antisymmetric + self-Pi-partner -> amplification"
hypothesis, evaluated at N=6.

Observation pattern from prior investigations:
  N=3: psi_2 is self-Pi (k=2, 4-k=2), antisymmetric, c_1 = 0.576 (amplified ~2.2x)
  N=5: psi_3 is self-Pi (k=3, 6-k=3), symmetric,    c_1 = 0.677 (not amplified)
  N=7: psi_4 is self-Pi (k=4, 8-k=4), antisymmetric, c_1 = 2.14  (amplified ~2.2x)

Pattern: amplification happens only when the self-Pi-partner mode is
ALSO antisymmetric, which occurs at N mod 4 == 3 (N = 3, 7, 11, ...).

At N=6 (even), there is NO self-Pi-partner in the single-excitation sector.
Pi-partner of k is (N+1-k) = 7-k, giving pairs (1,6), (2,5), (3,4); none fixed.
Also: Pi-pair partners have OPPOSITE reflection parity at even N (since N+1 is
odd, so k and N+1-k differ in parity).

Two simultaneous predictions at N=6:
  (P1) No mode shows c_1 amplification (no self-Pi-antisym mode exists).
  (P2) Pi-pair partners have opposite parity; if "Pi-pair c_1 identity" wins,
       antisym-small behaviour breaks for some partners; if "parity" wins,
       Pi-pair identity breaks.

Setup identical to c1_past_future_test.py: PTF-standard bonding-plus-vacuum
initial states, bond (0,1) perturbation, dJ = +/- 0.01 symmetric difference.

N=6: d=64, d^2 = 4096, L is 4096x4096. Eigendecomp ~1 min each. Three
decomps (L_A, L_B+, L_B-) + 6 propagations. Total ~5-10 minutes.
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
    GAMMA_0, J_UNIFORM, T_FINAL, N_STEPS, T_FIT_MAX,
    build_H_XY, build_liouvillian_matrix,
    vacuum_ket, single_excitation_mode, density_matrix,
    per_site_purity, fit_alpha,
)

RESULTS_DIR = Path(__file__).parent / "results" / "c1_even_N_degeneracy_test"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

N = 6
DJ_EXTRACT = 0.01


def v_effect(N):
    return 1.0 + np.cos(np.pi / N)


def bonding_plus_vacuum(N, k):
    v = vacuum_ket(N)
    psi = single_excitation_mode(N, k=k)
    return (v + psi) / np.sqrt(2.0)


def single_excitation_H_spectrum(N, J):
    return np.array([2.0 * J * np.cos(np.pi * k / (N + 1)) for k in range(1, N + 1)])


def reflection_parity_label(k):
    """(-1)^(k+1) under i <-> N-1-i reflection. +1 = symmetric, -1 = antisymmetric."""
    return "symmetric" if (k + 1) % 2 == 0 else "antisymmetric"


def pi_partner(k, N):
    return N + 1 - k


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


def main():
    print("=" * 70)
    print(f"c_1 even-N degeneracy test at N = {N}")
    print("=" * 70)
    print(f"  gamma_0 = {GAMMA_0}, J = {J_UNIFORM}, dJ = +/- {DJ_EXTRACT}")
    print(f"  V(N) = {v_effect(N):.4f}")

    Es = single_excitation_H_spectrum(N, J_UNIFORM)
    print(f"\n  Single-excitation H-spectrum at N={N}:")
    for k in range(1, N + 1):
        parity = reflection_parity_label(k)
        partner = pi_partner(k, N)
        self_pi = "SELF-Pi" if k == partner else f"Pi-partner of psi_{partner}"
        print(f"    psi_{k}: E = {Es[k-1]:+.4f}  [{parity}, {self_pi}]")
    print(f"\n  Predictions under the 'antisym + self-Pi = amplification' hypothesis:")
    print(f"    No self-Pi-partner at even N -> no amplification expected.")
    print(f"    At N={N}, all 3 Pi-pairs are distinct, cross-parity (sym <-> antisym).")
    print(f"    Will test: does Pi-pair c_1 identity hold across parity boundary?")

    # Build Liouvillians
    J_A = [J_UNIFORM] * (N - 1)
    J_B_plus = list(J_A); J_B_plus[0] = J_UNIFORM + DJ_EXTRACT
    J_B_minus = list(J_A); J_B_minus[0] = J_UNIFORM - DJ_EXTRACT

    print(f"\n  Building Liouvillians (d^2 = {4**N})...")
    t0 = time.time()
    L_A = build_liouvillian_matrix(build_H_XY(J_A, N), GAMMA_0, N)
    ev_A, VR_A, VLinv_A = eig_and_inv(L_A)
    print(f"  L_A eigendecomp: {time.time()-t0:.1f} s")
    del L_A

    t0 = time.time()
    L_B_plus = build_liouvillian_matrix(build_H_XY(J_B_plus, N), GAMMA_0, N)
    ev_Bp, VR_Bp, VLinv_Bp = eig_and_inv(L_B_plus)
    del L_B_plus
    print(f"  L_B+ eigendecomp: {time.time()-t0:.1f} s")

    t0 = time.time()
    L_B_minus = build_liouvillian_matrix(build_H_XY(J_B_minus, N), GAMMA_0, N)
    ev_Bm, VR_Bm, VLinv_Bm = eig_and_inv(L_B_minus)
    del L_B_minus
    print(f"  L_B- eigendecomp: {time.time()-t0:.1f} s")

    # Per-state c_1
    times_arr = np.linspace(0.0, T_FINAL, N_STEPS + 1)
    results = {}
    print(f"\n  Per-state c_1 (bond (0,1), dJ = +/- {DJ_EXTRACT}):")
    print(f"  {'k':>2} {'E_k':>10} {'parity':>15} {'partner':>8} {'c_1':>10}")
    for k in range(1, N + 1):
        phi = bonding_plus_vacuum(N, k)
        rho_0 = density_matrix(phi)
        rho_A = propagate(ev_A, VR_A, VLinv_A, rho_0, times_arr)
        rho_Bp = propagate(ev_Bp, VR_Bp, VLinv_Bp, rho_0, times_arr)
        rho_Bm = propagate(ev_Bm, VR_Bm, VLinv_Bm, rho_0, times_arr)
        P_A = per_site_purity(rho_A, N)
        P_Bp = per_site_purity(rho_Bp, N)
        P_Bm = per_site_purity(rho_Bm, N)
        alpha_p = np.zeros(N); alpha_m = np.zeros(N)
        for i in range(N):
            a, _ = fit_alpha(times_arr, P_A[:, i], P_Bp[:, i])
            alpha_p[i] = a
            a, _ = fit_alpha(times_arr, P_A[:, i], P_Bm[:, i])
            alpha_m[i] = a
        closure_p = float(np.sum(np.log(alpha_p)))
        closure_m = float(np.sum(np.log(alpha_m)))
        c_1 = (closure_p - closure_m) / (2 * DJ_EXTRACT)
        parity = reflection_parity_label(k)
        partner = pi_partner(k, N)
        print(f"  {k:>2d} {Es[k-1]:>+10.4f} {parity:>15} {partner:>8d} {c_1:>+10.4f}")
        results[f"psi_{k}"] = {
            "k": k, "E_k": float(Es[k-1]), "parity": parity,
            "pi_partner": partner, "c_1": float(c_1),
            "closure_plus": closure_p, "closure_minus": closure_m,
            "alpha_plus": alpha_p.tolist(),
            "alpha_minus": alpha_m.tolist(),
        }

    # Pi-pair analysis
    print(f"\n  Pi-pair c_1 identity check:")
    print(f"  {'pair':>12} {'k-parity':>15} {'partner-parity':>18} "
          f"{'c_1(k)':>10} {'c_1(partner)':>14} {'diff':>10}")
    for k in range(1, N + 1):
        partner = pi_partner(k, N)
        if k >= partner:
            continue
        c_k = results[f"psi_{k}"]["c_1"]
        c_p = results[f"psi_{partner}"]["c_1"]
        p_k = reflection_parity_label(k)
        p_p = reflection_parity_label(partner)
        print(f"  {k:>3d}<->{partner:<3d}  {p_k:>15} {p_p:>18} "
              f"{c_k:>+10.4f} {c_p:>+14.4f} {abs(c_k-c_p):>10.2e}")

    # Compare with prior N=7 data
    print(f"\n  Comparison across N (all bond (0,1), dJ = +/- {DJ_EXTRACT}):")
    print(f"  ----------------------------------------------------------")
    print(f"  N=3, psi_2 (antisym, self-Pi):           c_1 = +0.5756")
    print(f"  N=5, psi_3 (sym,     self-Pi):           c_1 = +0.6765")
    print(f"  N=6, no self-Pi:                         see above")
    print(f"  N=7, psi_4 (antisym, self-Pi):           c_1 = +2.1357")
    print(f"  ")
    print(f"  Hypothesis 'antisym + self-Pi -> amplification':")
    print(f"    -> predicts N=6 has no single outlier c_1 value")
    print(f"    -> predicts Pi-pair c_1 identity test will decide if parity or")
    print(f"       Pi-symmetry is the governing symmetry")

    # Save
    out = {
        "N": N, "gamma_0": GAMMA_0, "J": J_UNIFORM, "dJ": DJ_EXTRACT,
        "V_N": v_effect(N),
        "defect_bond": [0, 1],
        "H_spectrum": Es.tolist(),
        "results_by_state": results,
        "prior_comparison": {
            "N3_psi2_self_pi_antisym": 0.5756,
            "N5_psi3_self_pi_sym":     0.6765,
            "N7_psi4_self_pi_antisym": 2.1357,
        },
    }
    path = RESULTS_DIR / "c1_even_N_test.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nSaved: {path}")


if __name__ == "__main__":
    main()
