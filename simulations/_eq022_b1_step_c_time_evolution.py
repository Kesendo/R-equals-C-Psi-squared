#!/usr/bin/env python3
"""_eq022_b1_step_c_time_evolution.py — EQ-022 (b1) Step (c).

Compute the actual Q_SCALE observable via time evolution on block-L:

    S(t, J) = Σᵢ 2 |(ρᵢ(t))_{0,1}|²,   ρ(t) = exp(L(J) · t) · ρ₀_Dicke

Find Q_peak as the J/γ₀ at which |∂S/∂J| is maximized over (J, t).

Question Tom raised: does this connect to PTF (Perspectival Time Field)?
PTF's machinery is bilinear-in-ρ J-perturbation observables on per-site
purity. K_CC_pr is the same machinery on the (n, n+1)-block spatial-sum
coherence — same kind of observable, different chromaticity. If Q_peak
falls out of an eigenmode-structure pattern at the peak, it parallels
PTF's α_i framework.

After finding Q_peak numerically, examine the block-L eigenstructure
at Q_peak: which modes carry probe weight, what are their pure-rate
projections, and is there a chromaticity-c-specific signature?
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import framework as fw  # noqa: E402


def evolve_S_grid(N, n, gamma_0, J_grid, t_grid):
    """Build S(t, J) over the (J, t) grid. Returns S 2D array shape (len(J), len(t))."""
    D, M_H_per_bond, _, _ = fw.block_L_split_xy(N, n, gamma_0)
    M_H = sum(M_H_per_bond)
    rho0 = fw.dicke_block_probe(N, n)
    S_kernel = fw.spatial_sum_coherence_kernel(N, n)

    S = np.zeros((len(J_grid), len(t_grid)))
    for i, J in enumerate(J_grid):
        L = D + J * M_H
        evals, R = np.linalg.eig(L)
        c0 = np.linalg.solve(R, rho0)
        for j, t in enumerate(t_grid):
            rho_t = R @ (np.exp(evals * t) * c0)
            S[i, j] = float(np.real(rho_t.conj() @ S_kernel @ rho_t))
    return S


def find_Q_peak(J_grid, t_grid, S, gamma_0):
    """Find (Q*, t*) where |∂S/∂J| is maximized over the grid."""
    dJ = J_grid[1] - J_grid[0]
    dS_dJ = np.gradient(S, dJ, axis=0)
    abs_K = np.abs(dS_dJ)
    flat_idx = int(np.argmax(abs_K))
    i, j = np.unravel_index(flat_idx, abs_K.shape)
    return J_grid[i] / gamma_0, t_grid[j], abs_K[i, j], abs_K


def eigenstructure_at_J(N, n, gamma_0, J):
    """Return (evals, R, c0_dicke, rho0) at given (N, n, J)."""
    D, M_H_per_bond, _, _ = fw.block_L_split_xy(N, n, gamma_0)
    rho0 = fw.dicke_block_probe(N, n)
    L = D + J * sum(M_H_per_bond)
    evals, R = np.linalg.eig(L)
    R_inv = np.linalg.inv(R)
    c0 = R_inv @ rho0
    return evals, R, c0, rho0


def main():
    gamma_0 = 0.05

    print(f"# EQ-022 (b1) Step (c): time evolution + Q_peak via |∂S/∂J|_max")
    print(f"# gamma_0 = {gamma_0}")
    print()

    # Refined Q-grid + t-grid for sharper Q_peak resolution
    Q_grid = np.linspace(1.0, 3.0, 41)  # dQ = 0.05, focus on peak region
    J_grid = Q_grid * gamma_0
    t_grid = np.linspace(1.0, 15.0, 71)  # dt = 0.2, finer resolution around t*

    print(f"  Q grid: [{Q_grid[0]:.2f}, {Q_grid[-1]:.2f}], dQ = {Q_grid[1]-Q_grid[0]:.3f}")
    print(f"  t grid: [{t_grid[0]:.2f}, {t_grid[-1]:.2f}], dt = {t_grid[1]-t_grid[0]:.3f}")
    print()

    test_cases = [
        # focus on small N for speed; multiple c values
        (4, 1), (5, 1), (6, 1),  # c=2
        (5, 2), (6, 2), (7, 2),  # c=3
        (7, 3),                  # c=4
    ]

    print("## Q_peak from time evolution: |∂S/∂J| max over (J, t)")
    print()
    print(f"  {'c':>3} {'N':>3} {'n':>3} {'Q*':>7} {'t*':>7} {'|K|max':>10}  F86")
    results = {}
    for (N, n) in test_cases:
        c = fw.chromaticity(N, n)
        t0 = time.time()
        S = evolve_S_grid(N, n, gamma_0, J_grid, t_grid)
        Q_star, t_star, K_max, _ = find_Q_peak(J_grid, t_grid, S, gamma_0)
        elapsed = time.time() - t0
        f86 = {2: "1.4-1.6 wobble", 3: "1.6", 4: "1.8", 5: "1.8"}[c]
        results[(N, n)] = (c, Q_star, t_star, K_max)
        print(f"  {c:>3} {N:>3} {n:>3} {Q_star:>7.3f} {t_star:>7.2f} {K_max:>10.4f}  {f86}  ({elapsed:.1f} s)")
    print()

    # Now examine eigenstructure at Q_peak for one case: c=4 N=7
    print("## Eigenstructure at Q_peak for c=4 N=7")
    N, n = 7, 3
    c, Q_star, t_star, _ = results[(N, n)]
    J_star = Q_star * gamma_0
    print(f"  At Q*={Q_star:.2f}, J*={J_star:.4f}, t*={t_star:.2f}")
    evals, R, c0, rho0 = eigenstructure_at_J(N, n, gamma_0, J_star)
    print(f"  block dim = {R.shape[0]}, eigenvalues = {len(evals)}")

    # Find modes with largest |c0| (Dicke probe weight)
    abs_c0 = np.abs(c0)
    order = np.argsort(-abs_c0)
    print(f"  Top 12 modes by Dicke probe weight |c0_j|:")
    print(f"  {'rank':>4} {'idx':>5} {'|c0|':>10} {'Re(λ)/γ₀':>10} {'Im(λ)':>12} {'pure-rate?':>12}")
    pure_rates = [-2 * gamma_0 * (2*k - 1) for k in range(1, c + 1)]
    for rank in range(12):
        j = order[rank]
        lj = evals[j]
        re_per_g = lj.real / gamma_0
        # Is this mode close to a pure rate?
        pr_distances = [abs(lj.real - pr) for pr in pure_rates]
        nearest_pr_idx = int(np.argmin(pr_distances))
        nearest_pr = pure_rates[nearest_pr_idx]
        is_pure = abs(lj.real - nearest_pr) / gamma_0 < 0.05  # within 5% of γ₀
        is_pure_str = f"HD={2*nearest_pr_idx+1}" if is_pure else "dressed"
        print(f"  {rank+1:>4} {j:>5} {abs_c0[j]:>10.5f} {re_per_g:>10.4f} {lj.imag:>+12.4f}  {is_pure_str:>12}")

    # Cumulative weight in pure vs dressed modes
    pure_weight = 0.0
    dressed_weight = 0.0
    for j in range(len(evals)):
        re = evals[j].real
        pr_distances = [abs(re - pr) for pr in pure_rates]
        is_pure = min(pr_distances) / gamma_0 < 0.05
        w = abs(c0[j]) ** 2
        if is_pure:
            pure_weight += w
        else:
            dressed_weight += w
    total = pure_weight + dressed_weight
    print()
    print(f"  Total weight in pure-rate modes: {pure_weight:.4f} ({pure_weight/total:.2%})")
    print(f"  Total weight in dressed modes:   {dressed_weight:.4f} ({dressed_weight/total:.2%})")
    print()

    # ---- Per-bond Q_peak verification (analytical bond derivative) ----
    print("## Per-bond Q_peak: ∂S/∂J_b at each bond via eigenbasis derivative")
    print()

    # Strategy: for each J, eigendecompose L(J) ONCE.
    # Use the duhamel formula in eigenbasis to compute ∂ρ/∂J_b cheaply.
    Q_grid_pb = np.linspace(1.0, 3.0, 21)  # dQ = 0.1
    J_grid_pb = Q_grid_pb * gamma_0
    t_grid_pb = np.linspace(1.0, 13.0, 25)
    print(f"  Q grid: [{Q_grid_pb[0]:.2f}, {Q_grid_pb[-1]:.2f}], dQ = {Q_grid_pb[1]-Q_grid_pb[0]:.3f}")
    print(f"  t grid: [{t_grid_pb[0]:.2f}, {t_grid_pb[-1]:.2f}], dt = {t_grid_pb[1]-t_grid_pb[0]:.3f}")
    print()

    def compute_per_bond_Q_peak(N, n, gamma_0, J_grid, t_grid):
        c = fw.chromaticity(N, n)
        D, M_H_per_bond, _, _ = fw.block_L_split_xy(N, n, gamma_0)
        rho0 = fw.dicke_block_probe(N, n)
        S_kernel = fw.spatial_sum_coherence_kernel(N, n)
        n_bonds = N - 1
        M_H_total = sum(M_H_per_bond)

        # Track best (Q, t) per bond
        best_K = np.zeros(n_bonds)
        best_Q = np.full(n_bonds, np.nan)
        best_t = np.full(n_bonds, np.nan)

        for J in J_grid:
            L = D + J * M_H_total
            evals, R = np.linalg.eig(L)
            R_inv = np.linalg.inv(R)
            c0 = R_inv @ rho0
            # Precompute X_b = R^{-1} M_H_b R for each bond
            X_b_list = [R_inv @ Mb @ R for Mb in M_H_per_bond]
            # Precompute Y = R^† S_kernel R (for inner products)
            Y = R.conj().T @ S_kernel @ R
            # Note: we need ∂S/∂J_b = 2·Re[ρ^† S_kernel · ∂ρ/∂J_b]
            #       In eigenbasis: ρ(t) = R · (exp(λ t) ⊙ c0), so ρ^† S_kernel ρ
            #       = (exp(conj(λ)t) ⊙ c0^*)^T · R^† S_kernel R · (exp(λt) ⊙ c0)
            #       = (e^* ⊙ c0^*)^T Y (e ⊙ c0)
            #       where e[k] = exp(λ_k t).
            # And ∂ρ(t)/∂J_b in eigenbasis:
            #   ∂_b ρ(t) = R · F_b(t) · c0
            #   F_b(t)[j, k] = X_b[j,k] · I_jk(t)
            #   I_jj(t) = t · exp(λ_j t)
            #   I_jk(t) = (exp(λ_k t) - exp(λ_j t))/(λ_k - λ_j)  (j ≠ k)
            for t in t_grid:
                e = np.exp(evals * t)
                # diff matrix for I_jk(t)
                # avoid division by zero by setting diagonal correctly
                lam_j = evals[:, None]
                lam_k = evals[None, :]
                with np.errstate(divide='ignore', invalid='ignore'):
                    I_mat = np.where(
                        np.abs(lam_k - lam_j) > 1e-10,
                        (e[None, :] - e[:, None]) / (lam_k - lam_j),
                        t * e[:, None]
                    )
                # ρ(t) in eigenbasis: e ⊙ c0
                ec = e * c0
                rho_t = R @ ec
                Srho = S_kernel @ rho_t
                # ∂ρ/∂J_b = R · F_b(t) · c0  where F_b[j,k] = X_b[j,k] · I_mat[j,k]
                for b in range(n_bonds):
                    F_b = X_b_list[b] * I_mat  # element-wise
                    drho_b = R @ (F_b @ c0)
                    K = 2.0 * float(np.real(np.vdot(rho_t, S_kernel @ drho_b)))
                    K_abs = abs(K)
                    if K_abs > best_K[b]:
                        best_K[b] = K_abs
                        best_Q[b] = J / gamma_0
                        best_t[b] = t
        return c, best_Q, best_K, best_t

    for (N, n) in [(5, 2), (6, 2), (7, 2), (7, 3), (5, 1), (6, 1)]:
        t0 = time.time()
        c, q_per_bond, k_per_bond, t_per_bond = compute_per_bond_Q_peak(
            N, n, gamma_0, J_grid_pb, t_grid_pb)
        elapsed = time.time() - t0
        max_b = int(np.argmax(k_per_bond))
        print(f"  c={c}, N={N}, n={n}: ({elapsed:.1f} s)")
        print(f"    bond  Q*     t*     |K|max")
        for b in range(N - 1):
            mark = " ←" if b == max_b else ""
            print(f"    {b}     {q_per_bond[b]:.2f}   {t_per_bond[b]:.2f}   {k_per_bond[b]:.4f}{mark}")
        f86 = {2: "1.4-1.6 wobble", 3: "1.6", 4: "1.8"}[c]
        print(f"    F86 reference: {f86}")
        print()

    # Also at Q=20 (plateau) for comparison
    print("## Eigenstructure at Q=20 (plateau) for c=4 N=7 (comparison)")
    J_plateau = 20 * gamma_0
    evals_p, R_p, c0_p, _ = eigenstructure_at_J(N, n, gamma_0, J_plateau)
    pure_w_p = 0.0
    dressed_w_p = 0.0
    for j in range(len(evals_p)):
        re = evals_p[j].real
        pr_distances = [abs(re - pr) for pr in pure_rates]
        is_pure = min(pr_distances) / gamma_0 < 0.05
        w = abs(c0_p[j]) ** 2
        if is_pure:
            pure_w_p += w
        else:
            dressed_w_p += w
    total_p = pure_w_p + dressed_w_p
    print(f"  Total weight in pure-rate modes: {pure_w_p:.4f} ({pure_w_p/total_p:.2%})")
    print(f"  Total weight in dressed modes:   {dressed_w_p:.4f} ({dressed_w_p/total_p:.2%})")


if __name__ == "__main__":
    main()
