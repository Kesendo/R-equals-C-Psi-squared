#!/usr/bin/env python3
"""c1_bond_scan_multi_N.py

Parameterised bond-scan of c_1 for Dicke states at N in {3, 4, 5, 6}.

Generalises c1_bond_scan.py (hardcoded N=5). For each N, for each bond
b in {0, ..., N-2}, for each Dicke state |S_n> with n in {0, ..., N},
measures c_1 from alpha-fit on per-site purity (standard EQ-018 c_1
metric). Produces the full (bond, n) matrix.

Task context: EQ-019 bond-position dependence of K_diag = c_1(|S_n>).
The acoustics lens asks whether this matrix shows clean node/antinode
structure, comparable to sine-basis psi_k amplitudes at each bond's
sites.

N=6 is expensive (d^2 = 4096, ~1 min per eigendecomposition, 11 decomps
per bond gives ~5 min per bond, 5 bonds = 25 min total). Parameterise
the N list via CLI.

Usage:
  python c1_bond_scan_multi_N.py --N 3,4,5,6
"""
from __future__ import annotations

import argparse
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
    density_matrix, per_site_purity, fit_alpha,
)
from c1_bilinearity_test import dicke_state

RESULTS_DIR = Path(__file__).parent / "results" / "c1_bond_scan_multi_N"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DJ_EXTRACT = 0.01


def eig_and_inv(L):
    eigvals, V_R = eig(L)
    V_Linv = np.linalg.inv(V_R)
    return eigvals, V_R, V_Linv


def propagate(dec, rho_0, times):
    eigvals, V_R, V_Linv = dec
    rho0_vec = rho_0.flatten(order='F')
    c0 = V_Linv @ rho0_vec
    d = rho_0.shape[0]
    out = np.empty((len(times), d, d), dtype=complex)
    for k, t in enumerate(times):
        rho_vec_t = V_R @ (np.exp(eigvals * t) * c0)
        out[k] = rho_vec_t.reshape(d, d, order='F')
    return out


def measure_c1(rho_0, decA, decBp, decBm, N, times):
    rho_A = propagate(decA, rho_0, times)
    rho_Bp = propagate(decBp, rho_0, times)
    rho_Bm = propagate(decBm, rho_0, times)
    P_A = per_site_purity(rho_A, N)
    P_Bp = per_site_purity(rho_Bp, N)
    P_Bm = per_site_purity(rho_Bm, N)
    alpha_p = np.zeros(N); alpha_m = np.zeros(N)
    rmse = 0.0
    for i in range(N):
        a, r = fit_alpha(times, P_A[:, i], P_Bp[:, i])
        alpha_p[i] = a; rmse = max(rmse, r)
        a, r = fit_alpha(times, P_A[:, i], P_Bm[:, i])
        alpha_m[i] = a
    closure_p = float(np.sum(np.log(alpha_p)))
    closure_m = float(np.sum(np.log(alpha_m)))
    c_1 = (closure_p - closure_m) / (2 * DJ_EXTRACT)
    return c_1, alpha_p.tolist(), alpha_m.tolist(), rmse


def sine_amp(n_qubits, k, site_index):
    """Sine-basis amplitude: psi_k(i) = sqrt(2/(N+1)) sin(pi k (i+1)/(N+1))."""
    return np.sqrt(2.0 / (n_qubits + 1)) * np.sin(
        np.pi * k * (site_index + 1) / (n_qubits + 1)
    )


def run_one_N(N):
    print(f"\n{'=' * 70}")
    print(f"N = {N}   (d^2 = {4**N})")
    print(f"{'=' * 70}")
    J_A = [J_UNIFORM] * (N - 1)
    times_arr = np.linspace(0.0, T_FINAL, N_STEPS + 1)

    t0 = time.time()
    L_A = build_liouvillian_matrix(build_H_XY(J_A, N), GAMMA_0, N)
    decA = eig_and_inv(L_A); del L_A
    print(f"  L_A built + eigendecomposed in {time.time() - t0:.1f} s")

    S = {n: dicke_state(N, n) for n in range(N + 1)}
    state_densities = {n: density_matrix(S[n]) for n in range(N + 1)}

    bond_results = {}
    for b in range(N - 1):
        print(f"\n  Bond b = {b} ({b}, {b+1}):")
        J_Bp = list(J_A); J_Bp[b] = J_UNIFORM + DJ_EXTRACT
        J_Bm = list(J_A); J_Bm[b] = J_UNIFORM - DJ_EXTRACT
        t0 = time.time()
        L_Bp = build_liouvillian_matrix(build_H_XY(J_Bp, N), GAMMA_0, N)
        decBp = eig_and_inv(L_Bp); del L_Bp
        L_Bm = build_liouvillian_matrix(build_H_XY(J_Bm, N), GAMMA_0, N)
        decBm = eig_and_inv(L_Bm); del L_Bm
        print(f"    L_B+/- built + eigendecomposed in {time.time() - t0:.1f} s")

        c1_by_n = {}
        for n in range(N + 1):
            c1, _, _, rmse = measure_c1(state_densities[n], decA, decBp, decBm,
                                          N, times_arr)
            c1_by_n[n] = {"c_1": c1, "rmse": rmse}
        bond_results[b] = c1_by_n

        # Print row
        labels = "  ".join(f"|S_{n}>" for n in range(N + 1))
        values = "  ".join(f"{c1_by_n[n]['c_1']:+.4f}" for n in range(N + 1))
        print(f"    state:     {labels}")
        print(f"    c_1:       {values}")

    # F71 mirror check across bonds
    print(f"\n  F71 bond-mirror check:")
    for b in range((N - 1) // 2 + 1):
        b_mirror = (N - 2) - b
        if b > b_mirror:
            continue
        if b == b_mirror:
            print(f"    bond {b} is self-mirror (center)")
            continue
        diffs = []
        for n in range(N + 1):
            diff = abs(bond_results[b][n]["c_1"] - bond_results[b_mirror][n]["c_1"])
            diffs.append(diff)
        max_diff = max(diffs)
        print(f"    bond {b} <-> bond {b_mirror}: max |c_1 diff| over n = {max_diff:.2e}")

    # Sine-basis comparison (acoustics lens): for each bond b, compare
    # c_1(|S_n>) to "product of sine amplitudes at sites b and b+1" for mode k = n.
    print(f"\n  Acoustics lens: sine-basis psi_k(site) amplitudes at each bond.")
    print(f"  {'k':>3}   " + "  ".join(f"bond{b}:psi({b})psi({b+1})" for b in range(N - 1)))
    for k in range(1, N + 1):
        amps = [sine_amp(N, k, b) * sine_amp(N, k, b + 1) for b in range(N - 1)]
        print(f"  {k:>3}   " + "  ".join(f"{a:>14.4f}" for a in amps))

    # Export full bond × n matrix
    out = {
        "N": N, "gamma_0": GAMMA_0, "J": J_UNIFORM, "dJ": DJ_EXTRACT,
        "bonds_tested": list(range(N - 1)),
        "bond_results": {
            f"bond_{b}_{b+1}": {
                f"S_{n}": bond_results[b][n]["c_1"] for n in range(N + 1)
            }
            for b in range(N - 1)
        },
        "bond_rmse": {
            f"bond_{b}_{b+1}": {
                f"S_{n}": bond_results[b][n]["rmse"] for n in range(N + 1)
            }
            for b in range(N - 1)
        },
        "sine_basis_bond_products": {
            f"k_{k}": [sine_amp(N, k, b) * sine_amp(N, k, b + 1) for b in range(N - 1)]
            for k in range(1, N + 1)
        },
    }
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=str, default="3,4,5,6",
                        help="Comma-separated N values")
    args = parser.parse_args()
    N_list = [int(x) for x in args.N.split(",")]

    t_start = time.time()
    all_N = {}
    for N in N_list:
        all_N[str(N)] = run_one_N(N)

    out_path = RESULTS_DIR / f"bond_scan_multi_N_{'_'.join(str(n) for n in N_list)}.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(all_N, f, indent=2)
    print(f"\nSaved: {out_path}")
    print(f"Total walltime: {time.time() - t_start:.1f} s")


if __name__ == "__main__":
    main()
