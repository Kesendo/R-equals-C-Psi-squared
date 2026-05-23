"""Typewriter sweep: one canvas, one dot per J, accumulating.

Same canvas, same axes, same AT line in the background. Each J-step
stamps ONE point (the slow mode at that J) onto the (⟨n_XY⟩, |Re(λ)|)
log-log plane. 50 J-values from 1e-10 to 20 → 50 strokes, building
up the trajectory like a typewriter writing one character at a time.

Output:
  simulations/results/q_scaling_visual/q_sweep_typewriter.gif
  simulations/results/q_scaling_visual/q_sweep_typewriter_final.png
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, "simulations")

from _q_scaling_as_gamma_distribution import (  # noqa: E402
    all_pauli_matrices_flat,
    all_pauli_strings_n_xy,
    compute_eigenmode_n_xy,
)


def compute_or_load(N: int, gamma: float, n_points: int,
                    J_min: float, J_max: float, cache_file: Path
                    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if cache_file.exists():
        print(f"Loading cached sweep: {cache_file}")
        data = np.load(cache_file)
        return data["J_values"], data["slow_re"], data["slow_nxy"]

    print(f"Computing slow mode at {n_points} J values...")
    J_values = np.logspace(math.log10(J_min), math.log10(J_max), n_points)
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
        if np.any(non_kernel):
            idx = int(np.argmin(np.where(non_kernel, re_abs, np.inf)))
            slow_re[i] = re_abs[idx]
            slow_nxy[i] = n_xy[idx]
        else:
            slow_re[i] = 1e-16
            slow_nxy[i] = 1e-16
        if i % 5 == 0 or i == n_points - 1:
            print(f"  [{i+1:2d}/{n_points}] J={J:.2e}  "
                  f"|Re|={slow_re[i]:.2e}  ⟨n_XY⟩={slow_nxy[i]:.2e}  "
                  f"({time.time()-t0:.1f}s)")

    np.savez(cache_file, J_values=J_values, slow_re=slow_re, slow_nxy=slow_nxy)
    print(f"Cached: {cache_file}")
    return J_values, slow_re, slow_nxy


def main() -> None:
    N = 5
    gamma = 0.05
    n_points = 50
    J_min, J_max = 1e-10, 20.0

    out_dir = Path("simulations/results/q_scaling_visual")
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_file = out_dir / f"q_sweep_data_N{N}.npz"

    J_values, slow_re, slow_nxy = compute_or_load(
        N=N, gamma=gamma, n_points=n_points,
        J_min=J_min, J_max=J_max, cache_file=cache_file)

    slow_re = np.maximum(slow_re, 1e-16)
    slow_nxy = np.maximum(slow_nxy, 1e-16)
    Q_values = J_values / gamma

    fig, ax = plt.subplots(figsize=(11, 8))

    x_ref = np.logspace(-16, 1, 200)
    ax.loglog(x_ref, 2 * gamma * x_ref, 'k--', lw=1.5, alpha=0.5,
              label=f"AT line: |Re(λ)| = 2γ₀·⟨n_XY⟩")
    ax.axhline(2 * gamma, color="blue", ls=":", lw=1.5, alpha=0.5,
               label=f"F50 floor (2γ₀ = {2*gamma})")
    ax.axvline(1.0, color="blue", ls=":", lw=1.5, alpha=0.3,
               label="F50 vertex (⟨n_XY⟩ = 1)")

    ax.set_xlim(1e-16, 5)
    ax.set_ylim(1e-16, 1)
    ax.set_xlabel("⟨n_XY⟩_slow")
    ax.set_ylabel("|Re(λ)|_slow")
    ax.set_title(f"Typewriter: one dot per J, accumulating on one canvas\n"
                 f"(N={N} chain, γ₀={gamma} constant — γ₀ flows, J rotates the receiver)")
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(loc="upper left", fontsize=9)

    sc = ax.scatter(slow_nxy[:1], slow_re[:1],
                     c=np.log10(J_values[:1]),
                     cmap="viridis", s=90,
                     vmin=math.log10(J_min), vmax=math.log10(J_max),
                     edgecolors="black", linewidths=0.5, zorder=5)
    cb = plt.colorbar(sc, ax=ax, label="log₁₀(J)  [stroke order]")
    line, = ax.plot([], [], color="gray", alpha=0.45, lw=1.2, zorder=4)
    text = ax.text(0.98, 0.02, "", transform=ax.transAxes, fontsize=10,
                   verticalalignment="bottom", horizontalalignment="right",
                   family="monospace",
                   bbox=dict(facecolor="white", alpha=0.85, edgecolor="gray"))

    def update(frame: int):
        sc.set_offsets(np.column_stack([slow_nxy[:frame + 1],
                                          slow_re[:frame + 1]]))
        sc.set_array(np.log10(J_values[:frame + 1]))
        line.set_data(slow_nxy[:frame + 1], slow_re[:frame + 1])
        text.set_text(
            f"stroke {frame+1:>2d}/{n_points}\n"
            f"J  = {J_values[frame]:.2e}\n"
            f"Q  = {Q_values[frame]:.2e}\n"
            f"⟨n_XY⟩ = {slow_nxy[frame]:.2e}\n"
            f"|Re(λ)| = {slow_re[frame]:.2e}"
        )
        return sc, line, text

    print("Rendering animation (50 frames at 4 fps)...")
    anim = FuncAnimation(fig, update, frames=n_points, interval=250,
                         blit=False, repeat=True)
    gif_file = out_dir / "q_sweep_typewriter.gif"
    anim.save(gif_file, writer=PillowWriter(fps=4))
    print(f"GIF saved: {gif_file}")

    update(n_points - 1)
    static_file = out_dir / "q_sweep_typewriter_final.png"
    plt.savefig(static_file, dpi=120, bbox_inches="tight")
    print(f"Final-frame static: {static_file}")

    plt.close(fig)


if __name__ == "__main__":
    main()
