#!/usr/bin/env python3
"""A 4D projection study of the abstract F25 Bell+ radius drawing.

The exact dynamical input remains the local-Z-dephasing scalar radius
``r=f(1+f^2)/6`` with ``f=exp(-4 gamma t)``. Its unique finite ``r=1/4``
crossing is entangled (concurrence ``C=f>0``), while ``f=1/2`` supplies the
below-quarter control ``r=5/48`` and ``C=1/2``. Phase, arms, height, colour
and the 4D rotation are drawing choices, not an F25 geodesic.

The construction embeds ``(x,y,z)`` at ``w=0``, rotates the ``(x,w)`` plane,
then drops ``w``. At beta=90 degrees x-width vanishes, but y and z remain: the
shadow is a yz-plane family whose fixed-height rings are segments, not a
single line. The projection pinch is only a visual metaphor and computes no
exceptional point. The recurrence cusp and the specified F86 toy EP are
separate mathematical resemblances. The radial circle also contains
``c=-1/4``, where ``1-4c=2``; it is not the cardioid boundary or a horizon.

Produces exactly ``whirlpool_4d.png`` in the caller-selected output directory.
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
FRAMES = (
    (0, "0°  original shadow", ["#0a1236", "#2a8fd0", "#ffffff"]),
    (45, "45°", ["#0c1430", "#3f9fb8", "#ffe6b0"]),
    (90, "90°  yz-plane segments", ["#1a1408", "#caa23a", "#fff4d6"]),
    (135, "135°", ["#1a0a10", "#d2662a", "#ffd45e"]),
    (180, "180°  reflected shadow", ["#150418", "#b32c1e", "#fff0e0"]),
)


def cpsi_magnitude(t: np.ndarray) -> np.ndarray:
    """Return the F25 Bell+ scalar radius at gamma=1."""
    f = np.exp(-4.0 * t)
    return f * (1.0 + f * f) / 6.0


def main(output_dir) -> None:
    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)
    t = np.linspace(0.0, T_MAX, 1400)
    radius = cpsi_magnitude(t)
    norm = Normalize(0.0, START)
    arms = []
    for arm in range(N_ARMS):
        phase = 2.0 * np.pi * arm / N_ARMS - OMEGA * t
        arms.append((radius * np.cos(phase), radius * np.sin(phase),
                     Z_SCALE * radius))
    angle = np.linspace(0.0, 2.0 * np.pi, 240)
    ring_x = QUARTER * np.cos(angle)
    ring_y = QUARTER * np.sin(angle)
    ring_z = np.full_like(angle, Z_SCALE * QUARTER)

    fig = plt.figure(figsize=(20, 5.8))
    fig.patch.set_facecolor(BG)
    for index, (degrees, label, palette) in enumerate(FRAMES):
        ax = fig.add_subplot(1, len(FRAMES), index + 1, projection="3d")
        ax.set_facecolor(BG)
        cmap = LinearSegmentedColormap.from_list(f"frame-{degrees}", palette)
        c = np.cos(np.radians(degrees))
        for x, y, z in arms:
            xr = x * c
            points = np.array([xr, y, z]).T.reshape(-1, 1, 3)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            for width, alpha in ((3.0, 0.10), (1.25, 0.95)):
                line = Line3DCollection(segments, cmap=cmap, norm=norm)
                line.set_array(START - radius[:-1])
                line.set_linewidth(width)
                line.set_alpha(alpha)
                ax.add_collection3d(line)
        ax.plot(ring_x * c, ring_y, ring_z, color="#ffd56b", lw=2.0, alpha=0.9)
        ax.scatter([0], [0], [0], s=36, color="#ffffff", alpha=0.5,
                   edgecolors="none")
        ax.view_init(elev=55, azim=-90)
        ax.set_xlim(-0.33, 0.33)
        ax.set_ylim(-0.33, 0.33)
        ax.set_zlim(0.0, Z_SCALE * START)
        ax.set_box_aspect((1, 1, 0.72))
        ax.set_axis_off()
        ax.text2D(0.5, 0.02, label, transform=ax.transAxes, color="#9fb2d6",
                  fontsize=12, ha="center", style="italic")

    fig.suptitle("drawn 4D rotation and orthographic shadow of the F25 radius",
                 color="#cfe0ff", fontsize=13, y=0.97)
    fig.text(0.5, 0.072, "90° projection: yz-plane segments",
             color="#cfe0ff", fontsize=10.5, ha="center")
    fig.text(0.5, 0.043,
             "projection pinch is a visual metaphor; y and z survive; no EP is computed",
             color="#9fb2d6", fontsize=9.5, ha="center")
    fig.text(0.5, 0.016,
             "F25 controls: unique r=1/4 at f*=0.8612240997395736 with C=f>0; "
             "f=1/2 gives r=5/48<1/4 and C=1/2",
             color="#8298bd", fontsize=9, ha="center")
    plt.subplots_adjust(left=0.01, right=0.99, top=0.88, bottom=0.10, wspace=0.0)
    fig.savefig(output_root / "whirlpool_4d.png", dpi=150, facecolor=BG)
    plt.close(fig)
    print("saved: whirlpool_4d.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Render the Whirlpool 4D projection study.")
    parser.add_argument("--output-dir", required=True)
    arguments = parser.parse_args()
    main(arguments.output_dir)
