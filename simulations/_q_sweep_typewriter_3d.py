"""3D typewriter: polar cylinder, one stroke per J, closed loop visible.

Same data as the 2D typewriter, lifted into 3D:
  θ = log₁₀(J) wrapped to [0°, 360°]      (the parametric axis)
  r = ⟨n_XY⟩_slow                          (radial offset from cylinder axis)
  z = log₁₀|Re(λ)|_slow                    (height)

Both ends of the J sweep (J → 0 and J → 20) sit at the F50 vertex
(⟨n_XY⟩=1, |Re|=2γ₀=0.1), so θ=0 and θ=360° map to the SAME 3D point.
The trajectory is a literal closed loop in 3D space.

Output:
  simulations/results/q_scaling_visual/q_sweep_typewriter_3d.gif
  simulations/results/q_scaling_visual/q_sweep_typewriter_3d_final.png
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def main() -> None:
    N = 5
    gamma = 0.05
    out_dir = Path("simulations/results/q_scaling_visual")
    cache_file = out_dir / f"q_sweep_data_N{N}.npz"
    if not cache_file.exists():
        raise FileNotFoundError(
            f"Cache missing: {cache_file}. Run _q_sweep_typewriter.py first.")

    data = np.load(cache_file)
    J_values = data["J_values"]
    slow_re = np.maximum(data["slow_re"], 1e-16)
    slow_nxy = np.maximum(data["slow_nxy"], 1e-16)
    n_points = len(J_values)
    Q_values = J_values / gamma

    log_J = np.log10(J_values)
    theta = (log_J - log_J.min()) / (log_J.max() - log_J.min()) * 2 * math.pi
    r = slow_nxy
    z = np.log10(slow_re)
    x_polar = r * np.cos(theta)
    y_polar = r * np.sin(theta)

    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection="3d")

    theta_ref = np.linspace(0, 2 * math.pi, 120)
    for z_ref, label in [
        (-1, "F50 floor (2γ₀=0.1)"),
        (-3, "10⁻³"),
        (-6, "10⁻⁶"),
        (-9, "10⁻⁹"),
        (-12, "10⁻¹²"),
        (-15, "10⁻¹⁵"),
    ]:
        ax.plot(np.cos(theta_ref), np.sin(theta_ref),
                np.full_like(theta_ref, z_ref),
                color="black", alpha=0.18, lw=0.5)
        ax.text(1.15, 0, z_ref, label, fontsize=8, color="gray")

    ax.plot([1, 1], [0, 0], [z.min() - 0.5, 0],
            color="red", alpha=0.3, lw=1, ls="--")
    ax.scatter([1], [0], [-1], color="red", s=300, marker="X",
                edgecolors="black", linewidths=1.5,
                label="F50 vertex (start = end of typewriter)", zorder=10)

    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.4, 1.4)
    ax.set_zlim(-16.5, 0.2)
    ax.set_xlabel("r · cos(θ)")
    ax.set_ylabel("r · sin(θ)")
    ax.set_zlabel("log₁₀|Re(λ)|")
    ax.set_title(f"3D typewriter on a polar cylinder\n"
                 f"θ = log₁₀(J) wrapped to 360°,  r = ⟨n_XY⟩,  z = log₁₀|Re(λ)|\n"
                 f"(N={N} chain, γ₀={gamma} constant)")
    ax.view_init(elev=22, azim=-55)

    sc = ax.scatter(x_polar[:1], y_polar[:1], z[:1],
                    c=log_J[:1], cmap="viridis", s=110,
                    vmin=log_J.min(), vmax=log_J.max(),
                    edgecolors="black", linewidths=0.5, zorder=5)
    line, = ax.plot([], [], [], color="gray", alpha=0.55, lw=1.4, zorder=4)
    text = ax.text2D(0.02, 0.97, "", transform=ax.transAxes, fontsize=10,
                     verticalalignment="top", family="monospace",
                     bbox=dict(facecolor="white", alpha=0.85, edgecolor="gray"))
    cb = plt.colorbar(sc, ax=ax, shrink=0.55, pad=0.10,
                      label="log₁₀(J)  [stroke order]")
    ax.legend(loc="upper right", fontsize=9)

    def update(frame: int):
        sc._offsets3d = (x_polar[:frame + 1],
                         y_polar[:frame + 1],
                         z[:frame + 1])
        sc.set_array(log_J[:frame + 1])
        line.set_data_3d(x_polar[:frame + 1],
                         y_polar[:frame + 1],
                         z[:frame + 1])
        text.set_text(
            f"stroke {frame+1:>2d}/{n_points}\n"
            f"J  = {J_values[frame]:.2e}\n"
            f"Q  = {Q_values[frame]:.2e}\n"
            f"θ  = {math.degrees(theta[frame]):6.1f}°\n"
            f"r  = ⟨n_XY⟩ = {slow_nxy[frame]:.2e}\n"
            f"z  = log|Re|  = {z[frame]:6.2f}"
        )
        return sc, line, text

    print(f"Rendering 3D typewriter GIF ({n_points} frames at 4 fps)...")
    anim = FuncAnimation(fig, update, frames=n_points, interval=250,
                         blit=False, repeat=True)
    gif_file = out_dir / "q_sweep_typewriter_3d.gif"
    anim.save(gif_file, writer=PillowWriter(fps=4))
    print(f"GIF saved: {gif_file}")

    update(n_points - 1)
    static_file = out_dir / "q_sweep_typewriter_3d_final.png"
    plt.savefig(static_file, dpi=120, bbox_inches="tight")
    print(f"Final static: {static_file}")

    plt.close(fig)


if __name__ == "__main__":
    main()
