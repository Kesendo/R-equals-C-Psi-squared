"""3D visualisation of the slow-mode trajectory across 11 decades of Q.

Two views, side by side:

  Left  (standard 3D):
      X = log₁₀(Q),  Y = ⟨n_XY⟩_slow,  Z = log₁₀|Re(λ)|_slow.
      Trajectory connects the slow-mode position at each Q.
      Both endpoints (Q→0 and Q→∞) sit at the F50 vertex (⟨n_XY⟩=1, |Re|=2γ₀).
      The interior dips into the Lebensader valley.

  Right (polar cylinder = '360° + alle Z-Schichten'):
      θ = log₁₀(Q) wrapped to [0°, 360°].
      r = ⟨n_XY⟩_slow (F50 sits at r=1, Lebensader at r near 0).
      z = dominant joint-popcount sum (p_row + p_col).
      The trajectory is a closed loop: same starting and ending point,
      one full revolution through the Z-layers along the way.

Output: simulations/results/q_scaling_visual/q_sweep_3d_trajectory_N5.png
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, "simulations")

from framework.lindblad import lindbladian_z_dephasing  # noqa: E402

from _f1_chain_heisenberg_small_n_anchor import (  # noqa: E402
    chain_bonds,
    heisenberg_graph_h,
)
from _q_scaling_as_gamma_distribution import (  # noqa: E402
    all_pauli_matrices_flat,
    all_pauli_strings_n_xy,
)


def popcount(i: int) -> int:
    return bin(i).count("1")


def main() -> None:
    N = 5
    gamma = 0.05
    n_points = 40
    J_min, J_max = 1e-10, 20.0

    print("=" * 100)
    print(f"3D trajectory sweep: γ₀={gamma}, N={N} chain, J ∈ [{J_min:.0e}, {J_max}]")
    print(f"                     {n_points} log-spaced points")
    print("=" * 100)

    J_values = np.logspace(math.log10(J_min), math.log10(J_max), n_points)
    Q_values = J_values / gamma

    d = 1 << N
    pauli_mats = all_pauli_matrices_flat(N)
    n_xy_pauli = all_pauli_strings_n_xy(N).astype(float)

    pops = np.array([popcount(i) for i in range(d)], dtype=np.int32)
    row_pop = np.repeat(pops, d)
    col_pop = np.tile(pops, d)
    sec_flat = row_pop * (N + 1) + col_pop

    slow_re = np.zeros(n_points)
    slow_im = np.zeros(n_points)
    slow_nxy = np.zeros(n_points)
    slow_p_row = np.zeros(n_points, dtype=int)
    slow_p_col = np.zeros(n_points, dtype=int)

    t_total = time.time()
    for i, J in enumerate(J_values):
        t0 = time.time()
        bonds = chain_bonds(N)
        H = heisenberg_graph_h(N, bonds, J=J)
        L = lindbladian_z_dephasing(H, [gamma] * N)
        evals, evecs = np.linalg.eig(L)
        re_abs = np.abs(evals.real)
        non_kernel = re_abs > 1e-13
        if not np.any(non_kernel):
            slow_re[i] = np.nan
            slow_nxy[i] = np.nan
            continue
        idx = int(np.argmin(np.where(non_kernel, re_abs, np.inf)))
        vec = evecs[:, idx]

        amps = (pauli_mats.conj() @ vec) * math.sqrt(1.0 / d)
        ovr = np.abs(amps) ** 2
        ovr_sum = ovr.sum()
        nxy = (ovr * n_xy_pauli).sum() / max(ovr_sum, 1e-15)

        weights = np.abs(vec) ** 2
        hist = np.bincount(sec_flat, weights=weights, minlength=(N + 1) ** 2)
        dom_flat = int(hist.argmax())
        p_row, p_col = divmod(dom_flat, N + 1)

        slow_re[i] = re_abs[idx]
        slow_im[i] = evals[idx].imag
        slow_nxy[i] = nxy
        slow_p_row[i] = p_row
        slow_p_col[i] = p_col

        if i % 5 == 0 or i == n_points - 1:
            print(f"  [{i+1:2d}/{n_points}] J={J:.2e} Q={Q_values[i]:.2e}  "
                  f"|Re|={re_abs[idx]:.3e}  ⟨n_XY⟩={nxy:.5f}  "
                  f"sector=({p_row},{p_col})  ({time.time()-t0:.1f}s)")
    print(f"  Total: {time.time()-t_total:.1f}s")

    fig = plt.figure(figsize=(20, 9))
    log_q = np.log10(Q_values)
    sec_sum = slow_p_row + slow_p_col

    eps = 1e-16
    log_re = np.log10(np.maximum(slow_re, eps))

    ax1 = fig.add_subplot(121, projection="3d")
    sc1 = ax1.scatter(log_q, slow_nxy, log_re,
                       c=sec_sum, cmap="plasma", s=80,
                       vmin=0, vmax=2 * N)
    ax1.plot(log_q, slow_nxy, log_re, color="gray", alpha=0.4, lw=1.5)
    ax1.scatter([log_q[0], log_q[-1]],
                [slow_nxy[0], slow_nxy[-1]],
                [log_re[0], log_re[-1]],
                color="red", s=250, marker="X",
                edgecolors="black", linewidths=1.5,
                label="F50 endpoints (Q→0 and Q→∞)", zorder=5)
    ax1.set_xlabel("log₁₀(Q)")
    ax1.set_ylabel("⟨n_XY⟩_slow")
    ax1.set_zlabel("log₁₀|Re(λ)|_slow")
    ax1.set_title(f"Slow-mode trajectory across 11 decades of Q\n"
                  f"(N={N} chain, γ₀={gamma})\n"
                  f"Colour = dominant joint-popcount sum")
    cb1 = plt.colorbar(sc1, ax=ax1, shrink=0.55, pad=0.10)
    cb1.set_label("p_row + p_col")
    ax1.legend(loc="upper left", fontsize=9)

    ax2 = fig.add_subplot(122, projection="3d")
    theta = (log_q - log_q.min()) / (log_q.max() - log_q.min()) * 2 * math.pi
    r = slow_nxy
    z = sec_sum.astype(float)
    x_polar = r * np.cos(theta)
    y_polar = r * np.sin(theta)
    sc2 = ax2.scatter(x_polar, y_polar, z,
                       c=log_q, cmap="viridis", s=80)
    ax2.plot(x_polar, y_polar, z, color="gray", alpha=0.5, lw=1.5)
    ax2.scatter([x_polar[0], x_polar[-1]],
                [y_polar[0], y_polar[-1]],
                [z[0], z[-1]],
                color="red", s=250, marker="X",
                edgecolors="black", linewidths=1.5,
                label="F50 endpoints", zorder=5)

    for z_ref in np.unique(z):
        theta_ref = np.linspace(0, 2 * math.pi, 100)
        ax2.plot(np.cos(theta_ref), np.sin(theta_ref),
                 np.full_like(theta_ref, z_ref),
                 color="black", alpha=0.12, lw=0.6)
        ax2.plot(0.5 * np.cos(theta_ref), 0.5 * np.sin(theta_ref),
                 np.full_like(theta_ref, z_ref),
                 color="black", alpha=0.08, lw=0.4)

    ax2.set_xlabel("cos(θ),  θ = log₁₀(Q) wrapped to 360°")
    ax2.set_ylabel("sin(θ)")
    ax2.set_zlabel("dominant joint-popcount sum  p_row + p_col")
    ax2.set_title("360° rotation through Z layers\n"
                  "θ = Q-angle around cylinder, r = ⟨n_XY⟩\n"
                  "Closed loop: same start = end at F50 vertex")
    cb2 = plt.colorbar(sc2, ax=ax2, shrink=0.55, pad=0.10)
    cb2.set_label("log₁₀(Q)")
    ax2.legend(loc="upper left", fontsize=9)

    plt.tight_layout()
    out_dir = Path("simulations/results/q_scaling_visual")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "q_sweep_3d_trajectory_N5.png"
    plt.savefig(out_file, dpi=120, bbox_inches="tight")
    print(f"\nPlot saved: {out_file}")

    print()
    print(f"Sector trajectory across Q (N={N}):")
    print(f"{'Q':>14}  {'|Re|':>12}  {'⟨n_XY⟩':>10}  {'sector':>12}  {'p_sum':>6}")
    sample_idx = list(range(0, n_points, max(1, n_points // 14)))
    if (n_points - 1) not in sample_idx:
        sample_idx.append(n_points - 1)
    for i in sample_idx:
        Q = Q_values[i]
        sec_label = f"({slow_p_row[i]},{slow_p_col[i]})"
        psum = slow_p_row[i] + slow_p_col[i]
        print(f"{Q:14.3e}  {slow_re[i]:12.3e}  {slow_nxy[i]:10.5f}  "
              f"{sec_label:>12}  {psum:6d}")


if __name__ == "__main__":
    main()
