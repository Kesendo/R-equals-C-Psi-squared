#!/usr/bin/env python3
"""_eq022_b1_step_g_two_level_decomposition.py — derive f_class(x) from the
2-level effective Liouvillian.

Stage 2 of Item 1 (analytical f_class derivation): compute the per-bond
V_b matrices in the slowest-pair (HD=1, HD=3) 2-level subspace of c=2
chains, average within each bond class, and test whether the 2-level
analytical prediction reproduces the full-block-L K_b(Q) shape and the
bond-class HWHM/Q* split.

Workflow:
  1. For each c=2 chain (N=5..8):
     - Build full block-L via fw.block_L_split_xy
     - Build channel-uniform basis P (2D for c=2)
     - For each bond b, compute V_b = P† · M_H_per_bond[b] · P (2×2)
     - Bond-class averages: ⟨V⟩_int (mean over interior bonds), ⟨V⟩_end
       (mean over the two endpoint bonds)
     - Sanity: ∑_b V_b should equal the heuristic L_eff off-diagonal
       structure (same-sign +ig_eff)
  2. For each chain and bond class, run the 2-level analytical prediction:
     - L_eff_class = diag(−2γ₀, −6γ₀) + J · ⟨V⟩_class
     - Compute K_class(Q, t) directly in the 2-level model with the
       channel-uniform-projected probe and S_kernel
     - Time-peak K_class(Q) = max_t |K_class(Q, t)|
     - Find Q_peak, HWHM_left/Q_peak
  3. Compare 2-level prediction to full-block-L data from step_f:
     - Same Q_peak ± ?
     - Same HWHM_left/Q_peak ± ?
     - Same shape f_class(x)?

If the 2-level prediction matches: f_class(x) is derivable analytically from
⟨V⟩_class and the universal HWHM-/Q* values follow from the 2×2 algebra. If
it deviates: the heuristic 2-level reduction misses something, escalate to
Item 4 (full block-L derivation).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
import framework as fw  # noqa: E402

RESULTS_DIR = REPO_ROOT / "simulations" / "results" / "eq022_two_level_decomposition"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def per_bond_V_in_two_level(N, n, gamma_0):
    """For c=2 chain (n=1, N≥5): return V_b matrices in the (HD=1, HD=3)
    2-level subspace and the channel-uniform-projected probe + S_kernel.

    Returns dict with:
      V_b_list: list of (N-1) 2x2 complex matrices
      D_2lvl: diagonal pure-rate part of L (2x2)
      probe_2lvl: channel-uniform-projected Dicke probe (2-vector)
      S_kernel_2lvl: spatial-sum kernel projected onto the 2-level subspace
                     (2x2 quadratic form on probe coefficients)
      g_eff_total: the heuristic g_eff = abs(off-diagonal of sum_b V_b)
      P_channel: channel-uniform basis (full block dim × 2)
    """
    D_full, M_H_per_bond, P_n, P_np1 = fw.block_L_split_xy(N, n, gamma_0)
    P_channel, HDs = fw.hd_channel_basis(N, n)
    assert HDs == [1, 3], f"expected c=2 with HD=[1,3], got {HDs}"
    Mtot = D_full.shape[0]

    # Project D onto the 2-level subspace (should be diagonal with -2γ₀, -6γ₀)
    D_2lvl = P_channel.conj().T @ D_full @ P_channel

    # Per-bond V_b in 2-level subspace
    V_b_list = [P_channel.conj().T @ Mb @ P_channel for Mb in M_H_per_bond]

    # Sanity: ∑_b V_b should be off-diagonal with ±ig_eff
    V_total = sum(V_b_list)
    g_eff_total = abs(V_total[0, 1])

    # Probe in 2-level subspace
    rho0 = fw.dicke_block_probe(N, n)
    probe_2lvl = P_channel.conj().T @ rho0  # 2-vector

    # S_kernel projected onto 2-level
    S_full = fw.spatial_sum_coherence_kernel(N, n)
    S_kernel_2lvl = P_channel.conj().T @ S_full @ P_channel

    return dict(
        V_b_list=V_b_list, D_2lvl=D_2lvl, probe_2lvl=probe_2lvl,
        S_kernel_2lvl=S_kernel_2lvl, g_eff_total=g_eff_total,
        P_channel=P_channel, V_total=V_total,
    )


def two_level_K_curve(D_2lvl, V_b, probe_2lvl, S_kernel_2lvl, gamma_0,
                      Q_grid, t_grid):
    """Compute K_b(Q) curve in the 2-level model for a single V_b matrix.

    L_eff(J) = D_2lvl + J · V_b   (2x2)
    ρ(t) = e^(L_eff t) · probe_2lvl
    ∂ρ/∂J = ∫_0^t e^(L_eff(t-s)) · V_b · e^(L_eff s) · probe_2lvl ds
    K(Q, t) = 2 · Re(⟨ρ(t), S_kernel_2lvl · ∂ρ/∂J⟩)
    K(Q) = max_t |K(Q, t)|

    Same Duhamel + spectral-contour formulation as step_e per_bond_K_curve,
    but in the 2x2 subspace.
    """
    K_curve = np.zeros(len(Q_grid))
    for i_Q, Q in enumerate(Q_grid):
        J = Q * gamma_0
        L = D_2lvl + J * V_b
        evals, R = np.linalg.eig(L)
        R_inv = np.linalg.inv(R)
        c0 = R_inv @ probe_2lvl
        Vb_eig = R_inv @ V_b @ R

        K_max = 0.0
        for t in t_grid:
            e = np.exp(evals * t)
            lam_j = evals[:, None]
            lam_k = evals[None, :]
            with np.errstate(divide='ignore', invalid='ignore'):
                I_mat = np.where(np.abs(lam_k - lam_j) > 1e-10,
                                 (e[None, :] - e[:, None]) / (lam_k - lam_j),
                                 t * e[:, None])
            rho_t = R @ (e * c0)
            F = Vb_eig * I_mat
            drho = R @ (F @ c0)
            K = 2.0 * float(np.real(np.vdot(rho_t, S_kernel_2lvl @ drho)))
            if abs(K) > K_max:
                K_max = abs(K)
        K_curve[i_Q] = K_max
    return K_curve


def find_peak_with_interp(Q_grid, K_curve):
    """Same routine as step_e: parabolic peak interpolation + linear HWHM."""
    i_max = int(np.argmax(K_curve))
    K_max = K_curve[i_max]
    if 0 < i_max < len(Q_grid) - 1:
        x = Q_grid[i_max - 1: i_max + 2]
        y = K_curve[i_max - 1: i_max + 2]
        coefs = np.polyfit(x, y, 2)
        if coefs[0] < 0:
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
    half = K_max / 2.0
    hwhm_l = hwhm_r = None
    for i in range(i_max, -1, -1):
        if K_curve[i] < half:
            x0, x1 = Q_grid[i], Q_grid[i + 1]
            y0, y1 = K_curve[i], K_curve[i + 1]
            x_half = x0 + (half - y0) * (x1 - x0) / (y1 - y0)
            hwhm_l = Q_star - x_half
            break
    for i in range(i_max, len(Q_grid)):
        if K_curve[i] < half:
            x0, x1 = Q_grid[i - 1], Q_grid[i]
            y0, y1 = K_curve[i - 1], K_curve[i]
            x_half = x0 + (half - y0) * (x1 - x0) / (y1 - y0)
            hwhm_r = x_half - Q_star
            break
    return Q_star, K_max, hwhm_l, hwhm_r


def fmt_complex_matrix(M, prec=4):
    """Format 2×2 complex matrix for printing."""
    lines = []
    for i in range(M.shape[0]):
        row = []
        for j in range(M.shape[1]):
            z = M[i, j]
            if abs(z.real) < 10**(-prec) and abs(z.imag) < 10**(-prec):
                row.append(f"{'0':>{2*prec+5}}")
            else:
                row.append(f"{z.real:+.{prec}f}{z.imag:+.{prec}f}j")
        lines.append("  [" + "  ".join(row) + "]")
    return "\n".join(lines)


def main():
    gamma_0 = 0.05
    Q_grid = np.arange(0.20, 4.01, 0.025)
    t_peak_universal = 1.0 / (4.0 * gamma_0)
    t_grid = np.linspace(0.6 * t_peak_universal, 1.6 * t_peak_universal, 21)

    print("# EQ-022 (b1) Step (g): 2-level decomposition for c=2 chains")
    print("# Goal: extract V_b per-bond, examine bond-class structure,")
    print("# test 2-level analytical prediction against full-block-L data.")
    print(f"# γ₀ = {gamma_0}, Q ∈ [{Q_grid[0]}, {Q_grid[-1]}], dQ = 0.025")
    print()

    cases = [(5, 1), (6, 1), (7, 1), (8, 1)]  # c=2 N=5..8
    all_data = {}

    for (N, n) in cases:
        c = fw.chromaticity(N, n)
        n_bonds = N - 1
        print("=" * 72)
        print(f"# c={c}, N={N}, n={n}, n_bonds={n_bonds}")
        print("=" * 72)
        t0 = time.time()
        decomp = per_bond_V_in_two_level(N, n, gamma_0)
        print(f"   Decomposition computed in {time.time()-t0:.2f}s")
        print(f"   D_2lvl (diag entries should be -2γ₀, -6γ₀ = {-2*gamma_0}, {-6*gamma_0}):")
        print(fmt_complex_matrix(decomp['D_2lvl']))
        print(f"   probe_2lvl: {decomp['probe_2lvl']}")
        print(f"   sum(V_b) (should be off-diagonal with ±ig_eff_total):")
        print(fmt_complex_matrix(decomp['V_total']))
        print(f"   g_eff_total = |V_total[0,1]| = {decomp['g_eff_total']:.6f}")
        print(f"   S_kernel_2lvl:")
        print(fmt_complex_matrix(decomp['S_kernel_2lvl']))
        print()
        print("   Per-bond V_b matrices:")
        for b, V_b in enumerate(decomp['V_b_list']):
            label = "endpoint" if b in (0, n_bonds - 1) else "interior"
            print(f"   bond {b} ({label}):")
            print(fmt_complex_matrix(V_b))
            print(f"     diag entries: {V_b[0,0]:+.4f}, {V_b[1,1]:+.4f}")
            print(f"     off-diag |V_b[0,1]| = {abs(V_b[0,1]):.4f}, "
                  f"phase = {np.angle(V_b[0,1])/np.pi:.3f}π")

        # Bond-class averages
        endpoint_bonds = [0, n_bonds - 1]
        interior_bonds = list(range(1, n_bonds - 1))
        V_int_avg = sum(decomp['V_b_list'][b] for b in interior_bonds) / max(len(interior_bonds), 1)
        V_end_avg = sum(decomp['V_b_list'][b] for b in endpoint_bonds) / len(endpoint_bonds)
        print()
        print("   Bond-class averages:")
        print(f"   ⟨V⟩_int (mean of {len(interior_bonds)} interior bonds):")
        print(fmt_complex_matrix(V_int_avg))
        print(f"   ⟨V⟩_end (mean of 2 endpoint bonds):")
        print(fmt_complex_matrix(V_end_avg))

        # 2-level prediction with bond-class-averaged V
        print()
        print("   2-level prediction with ⟨V⟩_class:")
        K_int_pred = two_level_K_curve(decomp['D_2lvl'], V_int_avg,
                                       decomp['probe_2lvl'],
                                       decomp['S_kernel_2lvl'],
                                       gamma_0, Q_grid, t_grid)
        K_end_pred = two_level_K_curve(decomp['D_2lvl'], V_end_avg,
                                       decomp['probe_2lvl'],
                                       decomp['S_kernel_2lvl'],
                                       gamma_0, Q_grid, t_grid)
        Q_int, K_int, hl_int, hr_int = find_peak_with_interp(Q_grid, K_int_pred)
        Q_end, K_end, hl_end, hr_end = find_peak_with_interp(Q_grid, K_end_pred)
        print(f"     Interior pred:  Q*={Q_int:.4f}, |K|={K_int:.5f}, "
              f"HWHM-/Q*={hl_int/Q_int:.4f}" if hl_int else f"   Q*={Q_int:.4f}")
        print(f"     Endpoint pred:  Q*={Q_end:.4f}, |K|={K_end:.5f}, "
              f"HWHM-/Q*={hl_end/Q_end:.4f}" if hl_end else f"   Q*={Q_end:.4f}")

        all_data[f"c2_N{N}"] = dict(
            N=N, n=n, n_bonds=n_bonds,
            D_2lvl_diag=[float(decomp['D_2lvl'][i,i].real) for i in range(2)],
            probe_2lvl=[(z.real, z.imag) for z in decomp['probe_2lvl']],
            V_total_offdiag=(decomp['V_total'][0,1].real, decomp['V_total'][0,1].imag),
            g_eff_total=float(decomp['g_eff_total']),
            S_kernel_2lvl=[[(z.real, z.imag) for z in row] for row in decomp['S_kernel_2lvl']],
            V_b_per_bond=[[[(z.real, z.imag) for z in row] for row in V] for V in decomp['V_b_list']],
            V_int_avg=[[(z.real, z.imag) for z in row] for row in V_int_avg],
            V_end_avg=[[(z.real, z.imag) for z in row] for row in V_end_avg],
            pred_Q_peak_int=float(Q_int), pred_HWHM_left_int=float(hl_int) if hl_int else None,
            pred_Q_peak_end=float(Q_end), pred_HWHM_left_end=float(hl_end) if hl_end else None,
            pred_HWHM_l_over_Q_int=float(hl_int/Q_int) if hl_int else None,
            pred_HWHM_l_over_Q_end=float(hl_end/Q_end) if hl_end else None,
            K_int_pred=K_int_pred.tolist(),
            K_end_pred=K_end_pred.tolist(),
        )
        print()

    out = RESULTS_DIR / "summary.json"
    out.write_text(json.dumps({
        'gamma_0': gamma_0, 'Q_grid': Q_grid.tolist(),
        'cases': all_data,
    }, indent=2))
    print(f"# Saved to {out}")

    # Comparison summary
    print()
    print("=" * 72)
    print("# Comparison: 2-level prediction vs full-block-L step_f data")
    print("=" * 72)
    full_block_data = {
        # from step_f log
        'c2_N5': dict(int_Q=1.4821, int_h_over_Q=0.7455, end_Q=2.5008, end_h_over_Q=0.7700),
        'c2_N6': dict(int_Q=1.5801, int_h_over_Q=0.7529, end_Q=2.5470, end_h_over_Q=0.7738),
        'c2_N7': dict(int_Q=1.5831, int_h_over_Q=0.7507, end_Q=2.5299, end_h_over_Q=0.7738),
        'c2_N8': dict(int_Q=1.6049, int_h_over_Q=0.7531, end_Q=2.5145, end_h_over_Q=0.7734),
    }
    print(f"{'case':<8} {'class':<10} "
          f"{'pred Q*':>10} {'full Q*':>10} {'ΔQ%':>7}  "
          f"{'pred H-/Q':>11} {'full H-/Q':>11} {'ΔH%':>7}")
    print('-' * 90)
    for case_key in ['c2_N5', 'c2_N6', 'c2_N7', 'c2_N8']:
        d = all_data[case_key]
        f = full_block_data[case_key]
        for kind in ('int', 'end'):
            pq = d[f'pred_Q_peak_{"int" if kind=="int" else "end"}']
            fq = f[f'{kind}_Q']
            pho = d[f'pred_HWHM_l_over_Q_{"int" if kind=="int" else "end"}']
            fho = f[f'{kind}_h_over_Q']
            dq_pct = (pq - fq) / fq * 100 if fq else None
            dh_pct = (pho - fho) / fho * 100 if pho and fho else None
            print(f"{case_key:<8} {kind:<10} "
                  f"{pq:>10.4f} {fq:>10.4f} {dq_pct:>+6.2f}%  "
                  f"{pho:>11.4f} {fho:>11.4f} "
                  f"{dh_pct:>+6.2f}%" if dh_pct is not None else "")


if __name__ == "__main__":
    main()
