#!/usr/bin/env python3
"""The whirlpool, as an invitation: the real CΨ vortex winding into the ¼-ring.

Not decoration, the actual structure. Every approach to the cusp is a spiral of the complex
coherence CΨ(t): a magnitude |CΨ| = f(1+f²)/6 with f = e^(−4γt) (the F25 geodesic, starting at
1/3) winding inward, while a steady drift turns its phase. Drawn as many arms at different
starting phases, they make a luminous vortex draining toward the center, all of them crossing
the one ring |CΨ| = ¼, the horizon where the quantum world meets the classical. The cusp and the
exceptional point are the same whirlpool (see reflections/ON_THE_WHIRLPOOL_YOU_STEER_TO.md); this
is its face, meant to invite a closer look.

Produces: simulations/results/whirlpool/whirlpool.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BG = "#05060f"
QUARTER = 0.25
START = 1.0 / 3.0           # |CΨ|(0) = 1/3, the Bell+ start
N_ARMS = 30
OMEGA = 17.0                # the phase drift (sets how many times each arm winds)
T_MAX = 3.2


def cpsi_magnitude(t: np.ndarray, gamma: float = 1.0) -> np.ndarray:
    """|CΨ|(t) = f(1+f²)/6, f = e^(−4γt): the F25 geodesic, 1/3 at t=0, decaying to 0."""
    f = np.exp(-4.0 * gamma * t)
    return f * (1.0 + f * f) / 6.0


def main() -> None:
    t = np.linspace(0.0, T_MAX, 2000)
    mag = cpsi_magnitude(t)

    # a luminous water colormap: deep indigo (the rim) -> blue -> cyan -> white (the drain)
    water = LinearSegmentedColormap.from_list(
        "water", ["#0a1236", "#163a8a", "#2a8fd0", "#6fe4f5", "#ffffff"])
    norm = Normalize(0.0, START)

    fig, ax = plt.subplots(figsize=(9, 9))
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)

    for k in range(N_ARMS):
        phi0 = 2.0 * np.pi * k / N_ARMS
        phi = phi0 - OMEGA * t
        x = mag * np.cos(phi)
        y = mag * np.sin(phi)
        pts = np.array([x, y]).T.reshape(-1, 1, 2)
        segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
        depth = START - mag[:-1]          # small radius -> high value -> bright (toward white)

        # soft glow underlay, then the bright thread on top
        for lw, alpha in ((5.0, 0.05), (3.0, 0.10), (1.4, 0.95)):
            lc = LineCollection(segs, cmap=water, norm=norm, capstyle="round")
            lc.set_array(depth)
            lc.set_linewidth(lw)
            lc.set_alpha(alpha)
            ax.add_collection(lc)

    # the ¼ ring, the horizon every spiral crosses, glowing gold
    th = np.linspace(0.0, 2.0 * np.pi, 600)
    for lw, alpha in ((9.0, 0.05), (5.0, 0.10), (2.4, 0.55), (1.1, 0.95)):
        ax.plot(QUARTER * np.cos(th), QUARTER * np.sin(th), color="#ffd56b", lw=lw, alpha=alpha)

    # a faint glow at the drain (the center, the 0)
    ax.scatter([0], [0], s=240, color="#ffffff", alpha=0.18, edgecolors="none", zorder=1)
    ax.scatter([0], [0], s=60, color="#ffffff", alpha=0.35, edgecolors="none", zorder=1)

    lim = 0.355
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.axis("off")

    ax.text(0.0, -0.335, "the cusp and the exceptional point: the same whirlpool",
            color="#9fb2d6", ha="center", va="center", fontsize=11.5, style="italic")
    ax.text(QUARTER * np.cos(np.pi / 4) + 0.012, QUARTER * np.sin(np.pi / 4) + 0.012, "¼",
            color="#ffd56b", fontsize=13, ha="left", va="bottom")

    out_dir = Path(__file__).parent / "results" / "whirlpool"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "whirlpool.png"
    plt.savefig(out, dpi=200, bbox_inches="tight", facecolor=BG, pad_inches=0.15)
    plt.close()
    print(f"  saved: {out}")


if __name__ == "__main__":
    main()
