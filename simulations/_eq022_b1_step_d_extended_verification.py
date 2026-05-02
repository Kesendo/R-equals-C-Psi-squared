#!/usr/bin/env python3
"""_eq022_b1_step_d_extended_verification.py — verify F86 closed forms
beyond the original N=7 anchor.

Two scans:

(1) c=3 Interior Q_peak at N=8, N=9 — to test whether csc(π/5) = 1.7013
    really is the c=3 asymptote, or whether the N-trend (1.59 → 1.67 →
    1.70) continues past 1.7013.

(2) Endpoint Q_peak at N=5, 6, 7, 8 — to test csc(π/(N+1)) at multiple
    chain lengths. Currently only N=7 is verified.

Both use the analytical bond-derivative (Duhamel formula in eigenbasis,
one eigendecomp per J) for tractability. Block-L is the (n, n+1) coherence
block; Dicke probe; spatial-sum coherence kernel from the framework.
"""
from __future__ import annotations

import sys
import time
from math import pi, sin
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import framework as fw  # noqa: E402


def per_bond_Q_peak_scan(N, n, gamma_0, J_grid, t_grid):
    """Returns (best_Q, best_t, best_K) per bond via analytical bond derivative."""
    D, M_H_per_bond, _, _ = fw.block_L_split_xy(N, n, gamma_0)
    rho0 = fw.dicke_block_probe(N, n)
    S_kernel = fw.spatial_sum_coherence_kernel(N, n)
    n_bonds = N - 1
    M_H_total = sum(M_H_per_bond)

    best_K = np.zeros(n_bonds)
    best_Q = np.full(n_bonds, np.nan)
    best_t = np.full(n_bonds, np.nan)

    for J in J_grid:
        L = D + J * M_H_total
        evals, R = np.linalg.eig(L)
        R_inv = np.linalg.inv(R)
        c0 = R_inv @ rho0
        X_b_list = [R_inv @ Mb @ R for Mb in M_H_per_bond]
        for t in t_grid:
            e = np.exp(evals * t)
            lam_j = evals[:, None]
            lam_k = evals[None, :]
            with np.errstate(divide='ignore', invalid='ignore'):
                I_mat = np.where(np.abs(lam_k - lam_j) > 1e-10,
                                 (e[None, :] - e[:, None]) / (lam_k - lam_j),
                                 t * e[:, None])
            rho_t = R @ (e * c0)
            for b in range(n_bonds):
                F_b = X_b_list[b] * I_mat
                drho_b = R @ (F_b @ c0)
                K = 2.0 * float(np.real(np.vdot(rho_t, S_kernel @ drho_b)))
                K_abs = abs(K)
                if K_abs > best_K[b]:
                    best_K[b] = K_abs
                    best_Q[b] = J / gamma_0
                    best_t[b] = t
    return best_Q, best_t, best_K


def main():
    gamma_0 = 0.05
    csc_pi5 = 1.0 / sin(pi / 5)

    # --- Scan 1: c=3 Interior at N=8, 9 ---
    print(f"# Scan 1: c=3 Interior Q_peak at higher N (test csc(π/5) = {csc_pi5:.6f} asymptote)")
    print(f"# c=3 N=5 → 1.590, N=6 → 1.667, N=7 → 1.695 (prior data, dQ=0.01)")
    print()

    Q_grid = np.linspace(1.55, 1.85, 31)  # dQ = 0.01
    J_grid = Q_grid * gamma_0
    t_grid = np.linspace(4.0, 7.0, 16)  # dt = 0.2

    print(f"  Q grid: [{Q_grid[0]:.2f}, {Q_grid[-1]:.2f}], dQ = {Q_grid[1]-Q_grid[0]:.3f}")
    print(f"  t grid: [{t_grid[0]:.2f}, {t_grid[-1]:.2f}], dt = {t_grid[1]-t_grid[0]:.3f}")
    print()

    for N in (8, 9):
        n = 2  # c=3 at N=5,6,7,8,9 with n=2
        c = fw.chromaticity(N, n)
        if c != 3:
            print(f"  Skipped: (N={N}, n={n}) gives c={c}, not 3.")
            continue
        block_dim = len(fw.popcount_states(N, n)) * len(fw.popcount_states(N, n + 1))
        print(f"  Scanning c={c}, N={N}, n={n}, block dim = {block_dim}...", flush=True)
        t0 = time.time()
        best_Q, best_t, best_K = per_bond_Q_peak_scan(N, n, gamma_0, J_grid, t_grid)
        elapsed = time.time() - t0

        n_bonds = N - 1
        interior_Qs = best_Q[1:-1]  # exclude endpoints
        endpoint_Qs = [best_Q[0], best_Q[-1]]
        print(f"  ({elapsed:.1f}s) c={c} N={N} n={n}:")
        print(f"    Interior Q* = {interior_Qs}, mean = {interior_Qs.mean():.4f}")
        print(f"    csc(π/5) = {csc_pi5:.4f}, deviation = {abs(interior_Qs.mean() - csc_pi5):.4f}")
        print(f"    Endpoint Q* = {endpoint_Qs}")
        csc_endpoint = 1.0 / sin(pi / (N + 1))
        print(f"    csc(π/{N+1}) = {csc_endpoint:.4f}")
        print()

    # --- Scan 2: Endpoint at multiple N ---
    print(f"# Scan 2: Endpoint Q_peak at N=5, 6, 7, 8 (test csc(π/(N+1)) anchor)")
    print()

    for (N, n) in [(5, 2), (6, 2), (7, 2), (7, 3), (8, 3)]:
        c = fw.chromaticity(N, n)
        block_dim = len(fw.popcount_states(N, n)) * len(fw.popcount_states(N, n + 1))
        print(f"  Scanning c={c}, N={N}, n={n}, block dim = {block_dim}...", flush=True)

        # csc(π/6)=2, csc(π/7)=2.305, csc(π/8)=2.613, csc(π/9)=2.924
        Q_grid_ep = np.linspace(1.8, 3.0, 25)  # dQ = 0.05
        J_grid_ep = Q_grid_ep * gamma_0
        t_grid_ep = np.linspace(3.0, 7.0, 9)

        t0 = time.time()
        best_Q, best_t, best_K = per_bond_Q_peak_scan(N, n, gamma_0, J_grid_ep, t_grid_ep)
        elapsed = time.time() - t0

        endpoint_Qs = [best_Q[0], best_Q[-1]]
        endpoint_K = [best_K[0], best_K[-1]]
        csc_predict = 1.0 / sin(pi / (N + 1))
        print(f"  ({elapsed:.1f}s) c={c} N={N} n={n}:")
        print(f"    Endpoint Q* = {endpoint_Qs}, |K|max = {endpoint_K}")
        print(f"    csc(π/{N+1}) = {csc_predict:.4f}, deviation = {abs(np.mean(endpoint_Qs) - csc_predict):.4f}")
        print()


if __name__ == "__main__":
    main()
