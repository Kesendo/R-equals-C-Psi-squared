#!/usr/bin/env python3
"""_eq022_b1_step_e_resonance_shape.py — universality of the abs(K_CC_pr)(Q)
curve shape around Q_peak.

Hypothesis (2026-05-02): the position Q_peak is chain-specific (refuted
closed-form attempts), but the SHAPE of abs(K_CC_pr)(Q) around the peak
should reflect the universal EP-resonance physics: Q_EP = 2/g_eff is
the maximum-mixing point where the discriminant 4γ₀² − J²·g_eff² vanishes,
and the resonance form should follow from that.

Test: for multiple (c, N, Interior) cases, compute the full K(Q) curve at
fine grid (dQ = 0.025), find Q_peak and |K|max, normalize curves to
y = K(Q)/|K|max, x = Q − Q_peak. Compare.

If the normalized curves overlap: universal EP-resonance shape, Tier-1
candidate observation that survives the closed-form rollback.

If they don't: shape is also chain-specific. Then test alternative
normalizations (x = (Q − Q_peak)/HWHM, x = (Q − Q_peak)/Q_peak, etc.)
to find what universalises.

Output: JSON with per-case K(Q) curves + analysis summary.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import framework as fw  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = REPO_ROOT / "simulations" / "results" / "eq022_resonance_shape"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def per_bond_K_curve(N, n, gamma_0, J_grid, t_grid):
    """For each bond, compute the full |∂S/∂J_b|(Q) curve via the analytical
    bond-derivative (Duhamel formula, one eigendecomp per J).

    Returns dict per-bond:
      'Q_grid': Q values
      'K_curve': max |∂S/∂J_b| over t at each Q (peak-over-t per Q)
      'best_t_curve': t at which |K| peaks at each Q
    """
    D, M_H_per_bond, _, _ = fw.block_L_split_xy(N, n, gamma_0)
    rho0 = fw.dicke_block_probe(N, n)
    S_kernel = fw.spatial_sum_coherence_kernel(N, n)
    n_bonds = N - 1
    M_H_total = sum(M_H_per_bond)

    K_curves = np.zeros((n_bonds, len(J_grid)))
    t_curves = np.zeros((n_bonds, len(J_grid)))

    for i_J, J in enumerate(J_grid):
        L = D + J * M_H_total
        evals, R = np.linalg.eig(L)
        R_inv = np.linalg.inv(R)
        c0 = R_inv @ rho0
        X_b_list = [R_inv @ Mb @ R for Mb in M_H_per_bond]

        bond_K_at_J = np.zeros(n_bonds)
        bond_t_at_J = np.zeros(n_bonds)
        bond_K_max = np.zeros(n_bonds)
        bond_t_at_max = np.zeros(n_bonds)

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
                if K_abs > bond_K_max[b]:
                    bond_K_max[b] = K_abs
                    bond_t_at_max[b] = t

        K_curves[:, i_J] = bond_K_max
        t_curves[:, i_J] = bond_t_at_max

    return K_curves, t_curves


def find_peak_with_interp(Q_grid, K_curve):
    """Find peak Q* and |K|max with parabolic interpolation around the
    maximum grid point, plus HWHM on each side via linear interpolation.

    Returns:
      Q_star, K_max, hwhm_left, hwhm_right
    """
    i_max = int(np.argmax(K_curve))
    K_max = K_curve[i_max]

    # Parabolic interpolation around peak (if interior point)
    if 0 < i_max < len(Q_grid) - 1:
        x = Q_grid[i_max - 1: i_max + 2]
        y = K_curve[i_max - 1: i_max + 2]
        coefs = np.polyfit(x, y, 2)
        if coefs[0] < 0:  # downward parabola
            Q_star = -coefs[1] / (2 * coefs[0])
            K_max_interp = coefs[2] - coefs[1] ** 2 / (4 * coefs[0])
            if abs(Q_star - Q_grid[i_max]) <= (Q_grid[i_max + 1] - Q_grid[i_max - 1]):
                K_max = K_max_interp
            else:
                Q_star = Q_grid[i_max]
        else:
            Q_star = Q_grid[i_max]
    else:
        Q_star = Q_grid[i_max]

    # HWHM via linear interpolation
    half = K_max / 2.0
    hwhm_left = None
    hwhm_right = None
    # Left side: scan from i_max backwards until below half
    for i in range(i_max, -1, -1):
        if K_curve[i] < half:
            # Linear interp between i and i+1
            x0, x1 = Q_grid[i], Q_grid[i + 1]
            y0, y1 = K_curve[i], K_curve[i + 1]
            x_half = x0 + (half - y0) * (x1 - x0) / (y1 - y0)
            hwhm_left = Q_star - x_half
            break
    # Right side
    for i in range(i_max, len(Q_grid)):
        if K_curve[i] < half:
            x0, x1 = Q_grid[i - 1], Q_grid[i]
            y0, y1 = K_curve[i - 1], K_curve[i]
            x_half = x0 + (half - y0) * (x1 - x0) / (y1 - y0)
            hwhm_right = x_half - Q_star
            break

    return Q_star, K_max, hwhm_left, hwhm_right


def main():
    gamma_0 = 0.05

    # Wide Q grid covering pre-onset, peak, plateau
    Q_grid = np.arange(0.20, 4.01, 0.025)  # dQ = 0.025
    J_grid = Q_grid * gamma_0
    t_grid = np.linspace(3.0, 8.0, 21)  # dt = 0.25

    print(f"# EQ-022 (b1) Step (e): resonance shape universality test")
    print(f"# Q grid: [{Q_grid[0]:.3f}, {Q_grid[-1]:.3f}], dQ = {Q_grid[1]-Q_grid[0]:.3f}, {len(Q_grid)} points")
    print(f"# t grid: [{t_grid[0]:.1f}, {t_grid[-1]:.1f}], dt = {t_grid[1]-t_grid[0]:.2f}, {len(t_grid)} points")
    print(f"# γ₀ = {gamma_0}")
    print()

    # Cases to scan: (N, n) → (c, block_dim)
    cases = [
        (5, 2),  # c=3, dim=100
        (6, 2),  # c=3, dim=300
        (7, 2),  # c=3, dim=735
        (8, 2),  # c=3, dim=1568
        (7, 3),  # c=4, dim=1225
        (8, 3),  # c=4, dim=3920
    ]

    summary = []
    all_curves = {}

    for (N, n) in cases:
        c = fw.chromaticity(N, n)
        block_dim = len(fw.popcount_states(N, n)) * len(fw.popcount_states(N, n + 1))
        n_bonds = N - 1
        print(f"## c={c}, N={N}, n={n}, block dim = {block_dim}, n_bonds = {n_bonds}", flush=True)
        t0 = time.time()
        K_curves, t_curves = per_bond_K_curve(N, n, gamma_0, J_grid, t_grid)
        elapsed = time.time() - t0
        print(f"   ({elapsed:.1f}s)", flush=True)

        # For Interior bonds: average over interior bonds (exclude endpoints)
        interior_bonds = list(range(1, n_bonds - 1))
        endpoint_bonds = [0, n_bonds - 1]

        K_interior = K_curves[interior_bonds].mean(axis=0)
        K_endpoint = K_curves[endpoint_bonds].mean(axis=0)
        t_interior = t_curves[interior_bonds].mean(axis=0)
        t_endpoint = t_curves[endpoint_bonds].mean(axis=0)

        # Peak finding for interior and endpoint
        Q_int, K_int, hwhm_l_int, hwhm_r_int = find_peak_with_interp(Q_grid, K_interior)
        Q_ep, K_ep, hwhm_l_ep, hwhm_r_ep = find_peak_with_interp(Q_grid, K_endpoint)

        case_key = f"c{c}_N{N}_n{n}"
        all_curves[case_key] = {
            'c': c, 'N': N, 'n': n, 'block_dim': block_dim,
            'Q_grid': Q_grid.tolist(),
            'K_interior': K_interior.tolist(),
            'K_endpoint': K_endpoint.tolist(),
            't_interior': t_interior.tolist(),
            't_endpoint': t_endpoint.tolist(),
        }
        summary.append({
            'c': c, 'N': N, 'n': n, 'block_dim': block_dim,
            'Q_peak_interior': float(Q_int) if Q_int is not None else None,
            'K_max_interior': float(K_int) if K_int is not None else None,
            'hwhm_left_interior': float(hwhm_l_int) if hwhm_l_int is not None else None,
            'hwhm_right_interior': float(hwhm_r_int) if hwhm_r_int is not None else None,
            'Q_peak_endpoint': float(Q_ep) if Q_ep is not None else None,
            'K_max_endpoint': float(K_ep) if K_ep is not None else None,
            'hwhm_left_endpoint': float(hwhm_l_ep) if hwhm_l_ep is not None else None,
            'hwhm_right_endpoint': float(hwhm_r_ep) if hwhm_r_ep is not None else None,
            'elapsed_s': float(elapsed),
        })

        print(f"   Interior:  Q* = {Q_int:.4f}, |K|max = {K_int:.5f}", end="")
        if hwhm_l_int and hwhm_r_int:
            print(f", HWHM = ({hwhm_l_int:.3f}, {hwhm_r_int:.3f}), asym = {hwhm_r_int/hwhm_l_int:.3f}")
        else:
            print()
        print(f"   Endpoint:  Q* = {Q_ep:.4f}, |K|max = {K_ep:.5f}", end="")
        if hwhm_l_ep and hwhm_r_ep:
            print(f", HWHM = ({hwhm_l_ep:.3f}, {hwhm_r_ep:.3f}), asym = {hwhm_r_ep/hwhm_l_ep:.3f}")
        else:
            print()
        print()

    # Save
    out_summary = RESULTS_DIR / "summary.json"
    out_curves = RESULTS_DIR / "curves.json"
    with open(out_summary, 'w') as f:
        json.dump(summary, f, indent=2)
    with open(out_curves, 'w') as f:
        json.dump(all_curves, f, indent=2)
    print(f"# Saved summary to {out_summary}")
    print(f"# Saved curves to {out_curves}")

    # Universality analysis: print HWHM/Q_peak ratios per case
    print()
    print("# Universality candidates")
    print()
    print(f"{'case':<12} {'Q_peak':>8} {'|K|max':>8} {'HWHM-':>7} {'HWHM+':>7} {'HWHM-/Qp':>10} {'HWHM+/Qp':>10} {'asym':>6}")
    for s in summary:
        for kind in ('interior', 'endpoint'):
            Q = s.get(f'Q_peak_{kind}')
            K = s.get(f'K_max_{kind}')
            hl = s.get(f'hwhm_left_{kind}')
            hr = s.get(f'hwhm_right_{kind}')
            if Q is None:
                continue
            tag = f"c{s['c']}N{s['N']}{kind[:3]}"
            print(f"{tag:<12} {Q:>8.4f} {K:>8.5f}", end="")
            if hl and hr:
                print(f" {hl:>7.4f} {hr:>7.4f} {hl/Q:>10.4f} {hr/Q:>10.4f} {hr/hl:>6.3f}")
            else:
                print()


if __name__ == "__main__":
    main()
