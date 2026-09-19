#!/usr/bin/env python3
"""Two abstract whirlpool drawings anchored to one exact F25 scalar curve.

The dynamical input is only the Bell+ local-Z-dephasing radius
``r=f(1+f^2)/6`` with ``f=exp(-4 gamma t)``, ``gamma>0`` and ``t>=0``.
The imposed phase, 30 arms, colours, glow, handedness and eye are drawing
choices. Since ``dr/df=(1+3f^2)/6>0``, ``r=1/4`` has one finite crossing:
``f=0.8612240997395736...``. Bell+ is still entangled there because its
concurrence is ``C=f>0``. The control ``f=1/2`` gives
``r=5/48<1/4`` and ``C=1/2``. The radial quarter ring is a drawing reference.
It marks neither separability nor a horizon, and it is not the cardioid boundary.

The recurrence cusp at ``c=+1/4`` and the specified toy EP in
``experiments/F86_EP_THROUGH_THE_CLOCK.md`` are separate mathematical
objects. Their resemblance to the drawings is an invitation, not an
identity, and no F25 geodesic is inferred.

Produces exactly ``whirlpool.png`` and ``whirlpool_mirror.png`` in the
caller-selected output directory.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Circle
import numpy as np

BG = "#05060f"
QUARTER = 0.25
START = 1.0 / 3.0
N_ARMS = 30
OMEGA = 17.0
T_MAX = 3.2


def cpsi_magnitude(t: np.ndarray, gamma: float = 1.0) -> np.ndarray:
    """Return the F25 Bell+ scalar radius for local Z dephasing."""
    f = np.exp(-4.0 * gamma * t)
    return f * (1.0 + f * f) / 6.0


def make_canvas(colors, ring_color: str, handedness: float, eye: bool):
    """Build one deterministic drawing; none of its styling is measured."""
    t = np.linspace(0.0, T_MAX, 2000)
    radius = cpsi_magnitude(t)
    cmap = LinearSegmentedColormap.from_list("whirlpool-colors", colors)
    norm = Normalize(0.0, START)
    fig, ax = plt.subplots(figsize=(9, 9))
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)

    for arm in range(N_ARMS):
        phase = 2.0 * np.pi * arm / N_ARMS + handedness * OMEGA * t
        x = radius * np.cos(phase)
        y = radius * np.sin(phase)
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        for width, alpha in ((5.0, 0.05), (3.0, 0.10), (1.4, 0.95)):
            line = LineCollection(segments, cmap=cmap, norm=norm, capstyle="round")
            line.set_array(START - radius[:-1])
            line.set_linewidth(width)
            line.set_alpha(alpha)
            ax.add_collection(line)

    angle = np.linspace(0.0, 2.0 * np.pi, 600)
    for width, alpha in ((9.0, 0.05), (5.0, 0.10), (2.4, 0.55), (1.1, 0.95)):
        ax.plot(QUARTER * np.cos(angle), QUARTER * np.sin(angle),
                color=ring_color, lw=width, alpha=alpha)

    if eye:
        for eye_radius, alpha in ((0.085, 1.0), (0.10, 0.55), (0.12, 0.25)):
            ax.add_patch(Circle((0.0, 0.0), eye_radius, color=BG, alpha=alpha,
                                zorder=4, linewidth=0))
        ax.add_patch(Circle((0.0, 0.0), 0.085, fill=False,
                            edgecolor="#ff7a3c", lw=1.4, alpha=0.8, zorder=5))
    else:
        ax.scatter([0], [0], s=240, color="#ffffff", alpha=0.18,
                   edgecolors="none", zorder=4)
        ax.scatter([0], [0], s=60, color="#ffffff", alpha=0.35,
                   edgecolors="none", zorder=4)

    ax.text(QUARTER * np.cos(np.pi / 4) + 0.012,
            QUARTER * np.sin(np.pi / 4) + 0.012, "1/4 radial reference",
            color=ring_color, fontsize=11, ha="left", va="bottom")
    ax.set_xlim(-0.355, 0.355)
    ax.set_ylim(-0.355, 0.355)
    ax.set_aspect("equal")
    ax.axis("off")
    return fig, ax


def render_cool(output_path: Path) -> None:
    fig, ax = make_canvas(
        ["#0a1236", "#163a8a", "#2a8fd0", "#6fe4f5", "#ffffff"],
        "#ffd56b", -1.0, False)
    ax.text(0.0, -0.323, "cool inward drawing; phase and arms are chosen",
            color="#9fb2d6", ha="center", fontsize=10, style="italic")
    fig.text(0.5, 0.065, "drawing: Bell+ local-Z-dephasing radius",
             color="#cfe0ff", ha="center", fontsize=11)
    fig.text(0.5, 0.035,
             "F25 controls: unique r=1/4 at f*=0.8612240997395736 with C=f>0; "
             "f=1/2 gives r=5/48<1/4 and C=1/2",
             color="#9fb2d6", ha="center", fontsize=9.5)
    fig.savefig(output_path, dpi=200, bbox_inches="tight", facecolor=BG,
                pad_inches=0.15)
    plt.close(fig)


def render_warm(output_path: Path) -> None:
    fig, ax = make_canvas(
        ["#150418", "#5e1230", "#b32c1e", "#ef8a2b", "#ffd45e", "#fff6e0"],
        "#79e6ff", 1.0, True)
    ax.text(0.0, -0.323, "warm reflected drawing; the eye is interpretation",
            color="#9fb2d6", ha="center", fontsize=10, style="italic")
    fig.text(0.5, 0.065,
             "radial quarter ring; cardioid cusp and F86 toy EP are separate",
             color="#cfe0ff", ha="center", fontsize=10.5)
    fig.text(0.5, 0.035,
             "F25 controls: unique r=1/4 at f*=0.8612240997395736 with C=f>0; "
             "f=1/2 gives r=5/48<1/4 and C=1/2",
             color="#9fb2d6", ha="center", fontsize=9.5)
    fig.savefig(output_path, dpi=200, bbox_inches="tight", facecolor=BG,
                pad_inches=0.15)
    plt.close(fig)


def main(output_dir) -> None:
    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)
    render_cool(output_root / "whirlpool.png")
    render_warm(output_root / "whirlpool_mirror.png")
    print("saved: whirlpool.png")
    print("saved: whirlpool_mirror.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Render the two abstract Whirlpool drawings.")
    parser.add_argument("--output-dir", required=True)
    arguments = parser.parse_args()
    main(arguments.output_dir)
