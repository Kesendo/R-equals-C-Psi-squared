#!/usr/bin/env python3
"""Verify superposition of per-bond c_1 at N=5.

From pi_pair_closure_followup.py at N=5, psi_1+vac, delta_J = +0.01:
  bond (0,1): c_1 = +0.928
  bond (1,2): c_1 = -0.213
  bond (2,3): c_1 = -0.213
  bond (3,4): c_1 = +0.928
  Sum = +1.430

Linear-order superposition predicts: for uniform delta_J = 0.01 on ALL
four bonds, closure = (Sum_b c_1^(b)) * delta_J = 0.01430.

This script applies uniform delta_J and checks if the observed closure
matches the superposed prediction. A match confirms linear superposition.

Also test a specific linear combination that should CANCEL c_1: for instance,
delta_J^(0,1) + x * delta_J^(1,2) = 0 with x = c_1^(0,1) / c_1^(1,2) = 0.928 / 0.213
= 4.36. So delta_J on bond (0,1) = 0.01, delta_J on bond (1,2) = 0.01 * 4.36
should give closure ~ 0 at linear order.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from pi_pair_closure_investigation import (
    GAMMA_0, J_UNIFORM, T_FINAL, N_STEPS, T_FIT_MAX,
    build_H_XY, build_liouvillian_matrix,
    bonding_plus_vacuum, density_matrix, propagate_vectorised,
    per_site_purity, fit_alpha,
)

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def run_multi_bond(N, dJ_per_bond, state, label=""):
    """dJ_per_bond: list of delta_J for each bond (length N-1)."""
    J_A = [J_UNIFORM] * (N - 1)
    J_B = [J_UNIFORM + dj for dj in dJ_per_bond]
    H_A = build_H_XY(J_A, N); L_A = build_liouvillian_matrix(H_A, GAMMA_0, N)
    H_B = build_H_XY(J_B, N); L_B = build_liouvillian_matrix(H_B, GAMMA_0, N)
    times = np.linspace(0.0, T_FINAL, N_STEPS + 1)
    rho_0 = density_matrix(state)
    rho_A = propagate_vectorised(L_A, rho_0, times)
    rho_B = propagate_vectorised(L_B, rho_0, times)
    P_A = per_site_purity(rho_A, N)
    P_B = per_site_purity(rho_B, N)
    alpha = np.zeros(N)
    for i in range(N):
        a, _ = fit_alpha(times, P_A[:, i], P_B[:, i])
        alpha[i] = a
    closure = float(np.sum(np.log(alpha)))
    print(f"  {label}: dJ = " + "  ".join(f"{d:+.4f}" for d in dJ_per_bond))
    print(f"    alpha =   " + "  ".join(f"{a:.4f}" for a in alpha))
    print(f"    Sum ln = {closure:+.6f}")
    return closure, alpha.tolist()


def main():
    N = 5
    print("="*60)
    print(f"Superposition test at N={N}, state = psi_1 + vac")
    print("="*60)
    phi = bonding_plus_vacuum(N, k=1)

    # Baseline: no perturbation (sanity)
    run_multi_bond(N, [0, 0, 0, 0], phi, label="No perturbation (sanity)")

    # Per-bond c_1 values from earlier fit
    c1_bond = [+0.9284, -0.2133, -0.2133, +0.9284]
    sum_c1 = sum(c1_bond)
    print(f"\n  Per-bond c_1 from earlier fit: {c1_bond}")
    print(f"  Sum c_1 over bonds = {sum_c1:+.4f}")
    print(f"  Predicted closure at uniform dJ = 0.01: {sum_c1 * 0.01:+.5f}")

    # Test 1: uniform dJ on all bonds
    print(f"\n  Test 1: uniform dJ = 0.01 on all 4 bonds")
    cl1, _ = run_multi_bond(N, [0.01]*4, phi, label="Uniform +0.01")

    # Test 2: same but smaller
    print(f"\n  Test 2: uniform dJ = 0.001 on all 4 bonds (closer to linear)")
    cl2, _ = run_multi_bond(N, [0.001]*4, phi, label="Uniform +0.001")

    # Test 3: cancellation combination
    # c_1^(0,1) + x * c_1^(1,2) = 0 -> x = -c_1^(0,1) / c_1^(1,2) = -(0.9284 / -0.2133) = 4.35
    # With dJ_base = 0.001 on (0,1) and dJ_scale * dJ_base on (1,2):
    x = -c1_bond[0] / c1_bond[1]
    print(f"\n  Test 3: cancellation combination dJ^(0,1) : dJ^(1,2) = 1 : {x:.2f}")
    print(f"    Predicted closure = c_1^(0,1) * 0.001 + c_1^(1,2) * ({x}*0.001) * dJ = 0")
    run_multi_bond(N, [0.001, x*0.001, 0, 0], phi, label="Cancel (0,1)+(1,2)")

    # Test 4: same cancellation scaled up
    print(f"\n  Test 4: same cancellation at base 0.01")
    run_multi_bond(N, [0.01, x*0.01, 0, 0], phi, label="Cancel (0,1)+(1,2) at base 0.01")

    # Test 5: two-bond combination that sums c_1 to zero using all 4 bonds
    # Symmetric perturbation: dJ on (0,1) and (3,4), opposite dJ on (1,2), (2,3).
    # Ratio: 4 * 0.928 + 4 * (-0.213) * y = 0 won't work cleanly. Try dJ_outer = 1,
    # dJ_inner = (0.928 + 0.928) / (0.213 + 0.213) = 4.36.
    # So perturbation: +1 on (0,1), (3,4); +4.36 on (1,2), (2,3) -- wrong sign.
    # Correction: want c_1_total = 0.928 * 2 * u + (-0.213) * 2 * v = 0 with u, v same sign.
    # u = 1, v = 0.928 / 0.213 = 4.36 (same sign cancels)
    u = 0.001; v = 0.928 / 0.213 * 0.001
    print(f"\n  Test 5: symmetric cancellation: outer +{u:.5f}, inner +{v:.5f}")
    run_multi_bond(N, [u, v, v, u], phi, label="Symmetric cancel")


if __name__ == "__main__":
    main()
