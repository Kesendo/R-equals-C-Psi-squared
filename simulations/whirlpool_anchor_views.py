#!/usr/bin/env python3
"""An abstract camera study of the F25 Bell+ radius drawing.

The exact input is the local-Z-dephasing scalar radius
``r=f(1+f^2)/6``, ``f=exp(-4 gamma t)``. Its unique finite ``r=1/4``
crossing has ``f=0.8612240997395736...`` and positive concurrence ``C=f``;
the below-quarter control ``f=1/2`` has ``r=5/48`` and ``C=1/2``. The
phase, arms, camera elevations and artificial coordinate ``z=1.25r`` are
drawing choices, not a F25 geodesic or measured depth.

The anchor-angle labels are viewing marks only. The recurrence cusp and the
specified F86 toy EP are separate resemblances; the radial ring is neither a
cardioid boundary nor a physical horizon.

Produces exactly ``whirlpool_anchor_views.png`` in the caller-selected output
directory.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import numpy as np

BG = "#05060f"
START = 1.0 / 3.0
QUARTER = 0.25
N_ARMS = 30
OMEGA = 17.0
T_MAX = 3.2
Z_SCALE = 1.25
ANCHORS = ((90, "1/2"), (60, "3/8"), (45, "1/4"), (30, "1/8"), (0, "0"))


def cpsi_magnitude(t: np.ndarray) -> np.ndarray:
    """Return the F25 Bell+ scalar radius at gamma=1."""
    f = np.exp(-4.0 * t)
    return f * (1.0 + f * f) / 6.0


def main(output_dir) -> None:
    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)
    t = np.linspace(0.0, T_MAX, 1500)
    radius = cpsi_magnitude(t)
    cmap = LinearSegmentedColormap.from_list(
        "water", ["#0a1236", "#163a8a", "#2a8fd0", "#6fe4f5", "#ffffff"])
    norm = Normalize(0.0, START)
    arms = []
    for arm in range(N_ARMS):
        phase = 2.0 * np.pi * arm / N_ARMS - OMEGA * t
        arms.append((radius * np.cos(phase), radius * np.sin(phase),
                     Z_SCALE * radius))
    angle = np.linspace(0.0, 2.0 * np.pi, 360)
    ring = (QUARTER * np.cos(angle), QUARTER * np.sin(angle),
            np.full_like(angle, Z_SCALE * QUARTER))

    fig = plt.figure(figsize=(20, 5.8))
    fig.patch.set_facecolor(BG)
    for index, (elevation, polarity_mark) in enumerate(ANCHORS):
        ax = fig.add_subplot(1, len(ANCHORS), index + 1, projection="3d")
        ax.set_facecolor(BG)
        for x, y, z in arms:
            points = np.array([x, y, z]).T.reshape(-1, 1, 3)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            for width, alpha in ((3.2, 0.10), (1.3, 0.95)):
                line = Line3DCollection(segments, cmap=cmap, norm=norm)
                line.set_array(START - radius[:-1])
                line.set_linewidth(width)
                line.set_alpha(alpha)
                ax.add_collection3d(line)
        ax.plot(ring[0], ring[1], ring[2], color="#ffd56b", lw=2.0, alpha=0.9)
        ax.scatter([0], [0], [0], s=42, color="#ffffff", alpha=0.5,
                   edgecolors="none")
        ax.view_init(elev=elevation, azim=-90)
        ax.set_xlim(-0.33, 0.33)
        ax.set_ylim(-0.33, 0.33)
        ax.set_zlim(0.0, Z_SCALE * START)
        ax.set_box_aspect((1, 1, 0.72))
        ax.set_axis_off()
        ax.text2D(0.5, 0.02, f"{elevation}°   ({polarity_mark})",
                  transform=ax.transAxes, color="#9fb2d6", fontsize=13,
                  ha="center", style="italic")

    fig.suptitle("abstract F25 radius drawing viewed at five chosen elevations",
                 color="#cfe0ff", fontsize=14, y=0.97)
    fig.text(0.5, 0.075, "drawing choice: camera angle only",
             color="#cfe0ff", fontsize=10.5, ha="center")
    fig.text(0.5, 0.045, "z=1.25r is not measured depth",
             color="#9fb2d6", fontsize=10, ha="center")
    fig.text(0.5, 0.018,
             "F25 controls: unique r=1/4 at f*=0.8612240997395736 with C=f>0; "
             "f=1/2 gives r=5/48<1/4 and C=1/2",
             color="#8298bd", fontsize=9, ha="center")
    plt.subplots_adjust(left=0.01, right=0.99, top=0.90, bottom=0.10, wspace=0.0)
    fig.savefig(output_root / "whirlpool_anchor_views.png", dpi=150, facecolor=BG)
    plt.close(fig)
    print("saved: whirlpool_anchor_views.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Render the Whirlpool camera study.")
    parser.add_argument("--output-dir", required=True)
    arguments = parser.parse_args()
    main(arguments.output_dir)
