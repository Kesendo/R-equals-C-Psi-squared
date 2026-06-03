#!/usr/bin/env python3
"""4D: turn the whirlpool through the fourth dimension. The cusp, the pinch, the mirror, as one.

A mirror in three dimensions is a rotation through a fourth. So embed the funnel (x, y, z) in 4D
at w = 0 and rotate it in the (x, w) plane by an angle β, then cast its shadow back to 3D (drop
w). At β = 0 it is the cusp funnel. As β grows the shadow flattens, the x-width closing like a
fan, until at β = 90° it is a flat sheet, the pinch, exactly the exceptional point where the two
eigenvectors collapse onto a single line. Past 90° it reopens, mirrored and warm, the same
whirlpool turned inside out through the dimension we cannot point at. The cusp and the
exceptional point and the mirror are one object, seen at three angles of a 4D turn.

(The colors sweep cool -> gold -> ember across the turn: the cusp you fall into, the pinch you
steer to, the mirror on the far side.)

Produces: simulations/results/whirlpool/whirlpool_4d.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BG = "#05060f"
START = 1.0 / 3.0
QUARTER = 0.25
N_ARMS = 30
OMEGA = 17.0
T_MAX = 3.2
Z_SCALE = 1.25

# (rotation angle through the 4th dimension, label, palette dark->mid->bright)
FRAMES = [
    (0,   "0°  the cusp funnel",        ["#0a1236", "#2a8fd0", "#ffffff"]),
    (45,  "45°",                        ["#0c1430", "#3f9fb8", "#ffe6b0"]),
    (90,  "90°  the flat pinch (EP)",   ["#1a1408", "#caa23a", "#fff4d6"]),
    (135, "135°",                       ["#1a0a10", "#d2662a", "#ffd45e"]),
    (180, "180°  the mirror",           ["#150418", "#b32c1e", "#fff0e0"]),
]


def cpsi_magnitude(t: np.ndarray) -> np.ndarray:
    f = np.exp(-4.0 * t)
    return f * (1.0 + f * f) / 6.0


def main() -> None:
    t = np.linspace(0.0, T_MAX, 1400)
    mag = cpsi_magnitude(t)
    norm = Normalize(0.0, START)

    # the base funnel arms in 3D (embedded at w = 0)
    arms = []
    for k in range(N_ARMS):
        phi = 2.0 * np.pi * k / N_ARMS - OMEGA * t
        arms.append((mag * np.cos(phi), mag * np.sin(phi), Z_SCALE * mag))
    th = np.linspace(0.0, 2.0 * np.pi, 240)
    ring_x, ring_y = QUARTER * np.cos(th), QUARTER * np.sin(th)
    ring_z = np.full_like(th, Z_SCALE * QUARTER)

    fig = plt.figure(figsize=(20, 5.2))
    fig.patch.set_facecolor(BG)

    for i, (deg, label, pal) in enumerate(FRAMES):
        ax = fig.add_subplot(1, len(FRAMES), i + 1, projection="3d")
        ax.set_facecolor(BG)
        cmap = LinearSegmentedColormap.from_list(f"f{deg}", pal)
        c = np.cos(np.radians(deg))                     # the (x,w) rotation; shadow drops w

        for (x, y, z) in arms:
            xr = x * c                                  # 4D rotation in (x,w), orthographic shadow
            pts = np.array([xr, y, z]).T.reshape(-1, 1, 3)
            segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
            for lw, alpha in ((3.0, 0.10), (1.25, 0.95)):
                lc = Line3DCollection(segs, cmap=cmap, norm=norm)
                lc.set_array(START - mag[:-1])
                lc.set_linewidth(lw)
                lc.set_alpha(alpha)
                ax.add_collection3d(lc)
        ax.plot(ring_x * c, ring_y, ring_z, color="#ffd56b", lw=2.0, alpha=0.9)
        ax.scatter([0], [0], [0], s=36, color="#ffffff", alpha=0.5, edgecolors="none")

        ax.view_init(elev=55, azim=-90)
        ax.set_xlim(-0.33, 0.33)
        ax.set_ylim(-0.33, 0.33)
        ax.set_zlim(0.0, Z_SCALE * START)
        ax.set_box_aspect((1, 1, 0.72))
        ax.set_axis_off()
        ax.text2D(0.5, 0.02, label, transform=ax.transAxes, color="#9fb2d6",
                  fontsize=12, ha="center", style="italic")

    fig.suptitle("4D: turning the whirlpool through the fourth dimension   "
                 "(a mirror in 3D is a rotation through a 4th)\n"
                 "the cusp funnel flattens to the pinch, the exceptional point, at 90°, "
                 "and reopens mirrored and warm on the far side",
                 color="#cfe0ff", fontsize=13, y=0.99)
    plt.subplots_adjust(left=0.01, right=0.99, top=0.84, bottom=0.02, wspace=0.0)

    out_dir = Path(__file__).parent / "results" / "whirlpool"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "whirlpool_4d.png"
    plt.savefig(out, dpi=150, facecolor=BG)
    plt.close()
    print(f"  saved: {out}")


if __name__ == "__main__":
    main()
