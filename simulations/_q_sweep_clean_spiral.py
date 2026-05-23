"""Clean spiral: focus on the Q² Lebensader regime, no wrap, two revolutions.

The previous polar-cylinder version wrapped θ over the full 11-decade
J range and only showed half a rotation of actual spiral (Lebensader),
the rest being numerical F50-fallback and saturation plateau.

This version isolates the Q² Lebensader band: Q ∈ [0.01, 5], 40 points,
θ unwrapped with two full revolutions across the sweep. The trajectory
expands radially as it climbs in z (conical helix), because r = ⟨n_XY⟩
grows as Q² while z = log|Re| grows as 2·log(Q).

The geometry is a logarithmic spiral on a cone, slope 2 in (log Q, log|Re|).

Output:
  simulations/results/q_scaling_visual/q_sweep_clean_spiral.gif
  simulations/results/q_scaling_visual/q_sweep_clean_spiral_final.png
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, "simulations")

from _q_scaling_as_gamma_distribution import (  # noqa: E402
    all_pauli_matrices_flat,
    all_pauli_strings_n_xy,
    compute_eigenmode_n_xy,
)


def main() -> None:
    N = 5
    gamma = 0.05
    n_points = 40
    Q_min, Q_max = 0.01, 5.0
    n_revolutions = 2.0

    out_dir = Path("simulations/results/q_scaling_visual")
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_file = out_dir / f"q_sweep_clean_data_N{N}.npz"

    if cache_file.exists():
        print(f"Loading cached sweep: {cache_file}")
        data = np.load(cache_file)
        Q_values = data["Q_values"]
        slow_re = data["slow_re"]
        slow_nxy = data["slow_nxy"]
    else:
        print(f"Computing slow mode at {n_points} Q values "
              f"(Q ∈ [{Q_min}, {Q_max}], γ₀={gamma}, N={N})...")
        Q_values = np.logspace(math.log10(Q_min), math.log10(Q_max), n_points)
        J_values = Q_values * gamma
        pauli_mats = all_pauli_matrices_flat(N)
        n_xy_pauli = all_pauli_strings_n_xy(N).astype(float)
        slow_re = np.zeros(n_points)
        slow_nxy = np.zeros(n_points)
        for i, J in enumerate(J_values):
            t0 = time.time()
            re_abs, n_xy = compute_eigenmode_n_xy(
                N, J=J, gamma=gamma,
                pauli_mats=pauli_mats, n_xy_pauli=n_xy_pauli)
            non_kernel = re_abs > 1e-13
            idx = int(np.argmin(np.where(non_kernel, re_abs, np.inf)))
            slow_re[i] = re_abs[idx]
            slow_nxy[i] = n_xy[idx]
            if i % 5 == 0:
                print(f"  [{i+1:2d}/{n_points}] Q={Q_values[i]:.3f}  "
                      f"|Re|={slow_re[i]:.3e}  ⟨n_XY⟩={slow_nxy[i]:.4f}  "
                      f"({time.time()-t0:.1f}s)")
        np.savez(cache_file, Q_values=Q_values,
                 slow_re=slow_re, slow_nxy=slow_nxy)
        print(f"Cached: {cache_file}")

    log_Q = np.log10(Q_values)
    theta = ((log_Q - log_Q.min()) / (log_Q.max() - log_Q.min())
             * 2 * math.pi * n_revolutions)
    r = slow_nxy
    z = np.log10(np.maximum(slow_re, 1e-16))
    x_polar = r * np.cos(theta)
    y_polar = r * np.sin(theta)

    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection="3d")

    theta_ref = np.linspace(0, 2 * math.pi, 120)
    for z_ref in [-1, -2, -3, -4, -5, -6]:
        ax.plot(np.cos(theta_ref), np.sin(theta_ref),
                np.full_like(theta_ref, z_ref),
                color="black", alpha=0.10, lw=0.4)
        ax.text(1.10, 0.0, z_ref, f"{z_ref}", fontsize=8, color="gray")
    ax.plot([0, 0], [0, 0], [z.min() - 0.5, 0],
            color="gray", alpha=0.2, lw=0.5, ls=":")

    ax.scatter([1], [0], [-1], color="red", s=300, marker="X",
                edgecolors="black", linewidths=1.5,
                label=f"F50 vertex (Q ≈ {math.sqrt(N*N/0.55):.2f} ceiling)",
                zorder=10)

    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_zlim(z.min() - 0.5, 0)
    ax.set_xlabel("r · cos(θ)")
    ax.set_ylabel("r · sin(θ)")
    ax.set_zlabel("log₁₀|Re(λ)|")
    ax.set_title(f"Clean Q² Lebensader spiral on a cone\n"
                 f"θ unwrapped: {n_revolutions} revolutions across "
                 f"Q ∈ [{Q_min}, {Q_max}], slope 2 in (log Q, log|Re|)\n"
                 f"(N={N} chain, γ₀={gamma} constant)")
    ax.view_init(elev=22, azim=-58)

    sc = ax.scatter(x_polar[:1], y_polar[:1], z[:1],
                    c=log_Q[:1], cmap="viridis", s=90,
                    vmin=log_Q.min(), vmax=log_Q.max(),
                    edgecolors="black", linewidths=0.5, zorder=5)
    line, = ax.plot([], [], [], color="gray", alpha=0.55, lw=1.4, zorder=4)
    text = ax.text2D(0.02, 0.97, "", transform=ax.transAxes, fontsize=10,
                     verticalalignment="top", family="monospace",
                     bbox=dict(facecolor="white", alpha=0.85,
                                edgecolor="gray"))
    cb = plt.colorbar(sc, ax=ax, shrink=0.55, pad=0.10,
                      label="log₁₀(Q)  [stroke order]")
    ax.legend(loc="upper right", fontsize=9)

    def update(frame: int):
        sc._offsets3d = (x_polar[:frame + 1],
                         y_polar[:frame + 1],
                         z[:frame + 1])
        sc.set_array(log_Q[:frame + 1])
        line.set_data_3d(x_polar[:frame + 1],
                         y_polar[:frame + 1],
                         z[:frame + 1])
        text.set_text(
            f"stroke {frame+1:>2d}/{n_points}\n"
            f"Q  = {Q_values[frame]:7.3f}\n"
            f"θ  = {math.degrees(theta[frame]):6.1f}°\n"
            f"r  = ⟨n_XY⟩ = {slow_nxy[frame]:.4f}\n"
            f"z  = log|Re|  = {z[frame]:6.2f}"
        )
        return sc, line, text

    print(f"Rendering helix GIF ({n_points} frames at 5 fps)...")
    anim = FuncAnimation(fig, update, frames=n_points, interval=200,
                         blit=False, repeat=True)
    gif_file = out_dir / "q_sweep_clean_spiral.gif"
    anim.save(gif_file, writer=PillowWriter(fps=5))
    print(f"GIF saved: {gif_file}")

    update(n_points - 1)
    static_file = out_dir / "q_sweep_clean_spiral_final.png"
    plt.savefig(static_file, dpi=120, bbox_inches="tight")
    print(f"Final static: {static_file}")

    plt.close(fig)


if __name__ == "__main__":
    main()
