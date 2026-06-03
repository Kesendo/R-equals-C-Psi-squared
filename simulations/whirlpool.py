#!/usr/bin/env python3
"""The whirlpool, as an invitation: the real CΨ vortex winding into the ¼-ring, and its mirror.

Not decoration, the actual structure. Every approach to the cusp is a spiral of the complex
coherence CΨ(t): a magnitude |CΨ| = f(1+f²)/6 with f = e^(−4γt) (the F25 geodesic, starting at
1/3) winding inward, while a steady drift turns its phase. Drawn as many arms at different
starting phases, they make a luminous vortex draining toward the center, all of them crossing
the one ring |CΨ| = ¼, the horizon where the quantum world meets the classical.

Two faces, because the cusp and the exceptional point are the same whirlpool seen from two sides
(reflections/ON_THE_WHIRLPOOL_YOU_STEER_TO.md): a cool one winding one way (the cusp, the one you
fall into) and its warm mirror winding the other (the exceptional point, the one you steer to).

Produces: simulations/results/whirlpool/whirlpool.png and whirlpool_mirror.png
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
from matplotlib.patches import Circle

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


def render(out_name: str, colors, ring_color: str, subtitle: str, handedness: float,
           center: str, eye_color: str = "#ff7a3c") -> Path:
    """One whirlpool: arms colored by depth (rim dim -> bright), the glowing ¼ ring.
    handedness = -1 winds one way (the cusp), +1 the other (the exceptional point).
    center = 'drain' lights the middle (the cusp, light you fall into); center = 'eye' punches a
    dark pupil (the EP, the defective pinch where the eigenvectors collapse, the one you steer to)."""
    t = np.linspace(0.0, T_MAX, 2000)
    mag = cpsi_magnitude(t)
    cmap = LinearSegmentedColormap.from_list(out_name, colors)
    norm = Normalize(0.0, START)

    fig, ax = plt.subplots(figsize=(9, 9))
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)

    for k in range(N_ARMS):
        phi0 = 2.0 * np.pi * k / N_ARMS
        phi = phi0 + handedness * OMEGA * t
        x = mag * np.cos(phi)
        y = mag * np.sin(phi)
        pts = np.array([x, y]).T.reshape(-1, 1, 2)
        segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
        depth = START - mag[:-1]          # small radius -> high value -> bright (toward white)
        for lw, alpha in ((5.0, 0.05), (3.0, 0.10), (1.4, 0.95)):
            lc = LineCollection(segs, cmap=cmap, norm=norm, capstyle="round")
            lc.set_array(depth)
            lc.set_linewidth(lw)
            lc.set_alpha(alpha)
            ax.add_collection(lc)

    th = np.linspace(0.0, 2.0 * np.pi, 600)
    for lw, alpha in ((9.0, 0.05), (5.0, 0.10), (2.4, 0.55), (1.1, 0.95)):
        ax.plot(QUARTER * np.cos(th), QUARTER * np.sin(th), color=ring_color, lw=lw, alpha=alpha)

    if center == "drain":
        ax.scatter([0], [0], s=240, color="#ffffff", alpha=0.18, edgecolors="none", zorder=4)
        ax.scatter([0], [0], s=60, color="#ffffff", alpha=0.35, edgecolors="none", zorder=4)
    else:  # 'eye': a dark pupil punched into the swirl (the defective pinch)
        for r, a in ((0.085, 1.0), (0.10, 0.55), (0.12, 0.25)):
            ax.add_patch(Circle((0.0, 0.0), r, color=BG, alpha=a, zorder=4, linewidth=0))
        ax.add_patch(Circle((0.0, 0.0), 0.085, fill=False, edgecolor=eye_color, lw=1.4,
                            alpha=0.8, zorder=5))

    lim = 0.355
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.text(0.0, -0.335, subtitle, color="#9fb2d6", ha="center", va="center",
            fontsize=11.5, style="italic")
    lab_color = ring_color if ring_color != "#79e6ff" else "#bfefff"
    ax.text(QUARTER * np.cos(np.pi / 4) + 0.012, QUARTER * np.sin(np.pi / 4) + 0.012, "¼",
            color=lab_color, fontsize=13, ha="left", va="bottom")

    out_dir = Path(__file__).parent / "results" / "whirlpool"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / out_name
    plt.savefig(out, dpi=200, bbox_inches="tight", facecolor=BG, pad_inches=0.15)
    plt.close()
    return out


def main() -> None:
    # the cusp: cool water, winding one way into a bright drain (the light you fall into)
    a = render(
        "whirlpool.png",
        ["#0a1236", "#163a8a", "#2a8fd0", "#6fe4f5", "#ffffff"],
        "#ffd56b",
        "the cusp and the exceptional point: the same whirlpool",
        handedness=-1.0,
        center="drain")
    # the exceptional point: warm embers, turned the other way, around a dark pinch-eye
    b = render(
        "whirlpool_mirror.png",
        ["#150418", "#5e1230", "#b32c1e", "#ef8a2b", "#ffd45e", "#fff6e0"],
        "#79e6ff",
        "the same whirlpool, turned: the dark pinch you steer to",
        handedness=1.0,
        center="eye")
    print(f"  saved: {a}")
    print(f"  saved: {b}")


if __name__ == "__main__":
    main()
