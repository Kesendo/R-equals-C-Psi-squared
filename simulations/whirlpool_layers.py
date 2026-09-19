#!/usr/bin/env python3
"""An abstract layered drawing around the exact F25 Bell+ radius.

For local Z dephasing the only dynamical input is
``r=f(1+f^2)/6``, ``f=exp(-4 gamma t)``. The finite ``r=1/4`` crossing
is unique and remains entangled because concurrence is ``C=f>0``. The
control ``f=1/2`` gives ``r=5/48<1/4`` and ``C=1/2``. Horizontal plates,
their heights, the 30 arms and the funnel are drawing choices; they are not a
F25 geodesic, physical strata, or a classicality horizon.

The recurrence cusp and the specified F86 toy EP are separate mathematical
resemblances. The radial quarter ring is not the cardioid boundary.

Produces exactly ``whirlpool_layers.png`` in the caller-selected output
directory.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection
import numpy as np

BG = "#05060f"
START = 1.0 / 3.0
QUARTER = 0.25
N_ARMS = 30
OMEGA = 17.0
T_MAX = 3.2
Z_SCALE = 1.25
LAYERS = (1.0 / 3.0, 0.295, QUARTER, 0.205, 0.165, 0.125,
          0.092, 0.063, 0.038, 0.018)


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
    angle = np.linspace(0.0, 2.0 * np.pi, 240)

    fig = plt.figure(figsize=(10, 9.2))
    fig.patch.set_facecolor(BG)
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor(BG)
    for layer_radius in LAYERS:
        z = Z_SCALE * layer_radius
        is_quarter = abs(layer_radius - QUARTER) < 1e-6
        ring_color = "#ffd56b" if is_quarter else cmap(norm(START - layer_radius))
        vertices = [list(zip(layer_radius * np.cos(angle),
                             layer_radius * np.sin(angle),
                             np.full_like(angle, z)))]
        plate = Poly3DCollection(
            vertices, facecolor="#ffd56b" if is_quarter else "#2a8fd0",
            alpha=0.10 if is_quarter else 0.05, edgecolor="none")
        ax.add_collection3d(plate)
        ax.plot(layer_radius * np.cos(angle), layer_radius * np.sin(angle),
                np.full_like(angle, z), color=ring_color,
                lw=3.0 if is_quarter else 1.4,
                alpha=0.95 if is_quarter else 0.7)

    for arm in range(N_ARMS):
        phase = 2.0 * np.pi * arm / N_ARMS - OMEGA * t
        x = radius * np.cos(phase)
        y = radius * np.sin(phase)
        z = Z_SCALE * radius
        points = np.array([x, y, z]).T.reshape(-1, 1, 3)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        line = Line3DCollection(segments, cmap=cmap, norm=norm)
        line.set_array(START - radius[:-1])
        line.set_linewidth(1.0)
        line.set_alpha(0.55)
        ax.add_collection3d(line)

    ax.scatter([0], [0], [0], s=60, color="#ffffff", alpha=0.6,
               edgecolors="none")
    ax.text(QUARTER, 0.0, Z_SCALE * QUARTER, "  1/4  radial reference",
            color="#ffd56b", fontsize=11, ha="left", va="center")
    ax.text(0.125, 0.0, Z_SCALE * 0.125, "  1/8 drawn layer",
            color="#9fd0ff", fontsize=10, ha="left", va="center")
    ax.view_init(elev=24, azim=-72)
    ax.set_xlim(-0.33, 0.33)
    ax.set_ylim(-0.33, 0.33)
    ax.set_zlim(0.0, Z_SCALE * START)
    ax.set_box_aspect((1, 1, 0.82))
    ax.set_axis_off()

    fig.suptitle("layered abstract view of the Bell+ radius drawing",
                 color="#cfe0ff", fontsize=13, y=0.96)
    fig.text(0.5, 0.045,
             "drawing layers; 1/4 is a radial reference, not a horizon",
             color="#cfe0ff", fontsize=10.5, ha="center")
    fig.text(0.5, 0.018,
             "F25 controls: unique r=1/4 at f*=0.8612240997395736 with C=f>0; "
             "f=1/2 gives r=5/48<1/4 and C=1/2",
             color="#8298bd", fontsize=9, ha="center")
    plt.subplots_adjust(left=0.0, right=1.0, top=0.92, bottom=0.08)
    fig.savefig(output_root / "whirlpool_layers.png", dpi=170, facecolor=BG)
    plt.close(fig)
    print("saved: whirlpool_layers.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Render the Whirlpool layer drawing.")
    parser.add_argument("--output-dir", required=True)
    arguments = parser.parse_args()
    main(arguments.output_dir)
