#!/usr/bin/env python3
"""The whirlpool was a top view. Tilt our viewpoint through the anchor angles, in 3D.

The luminous 2D whirlpool (simulations/whirlpool.py) is the bird's-eye view, looking straight
down, of a 3D object: a funnel. Each spiral arm of the complex coherence CΨ winds down it, the
magnitude |CΨ| = f(1+f²)/6 (the F25 geodesic) serving as both the radius and the depth, so the
arms descend a cone toward the drain at the center, all crossing the ring |CΨ| = ¼.

The point is not to spin the picture; it is to move OUR VIEWPOINT. We tilt the camera through
the framework's anchor angles, the canonical heading marks where the polarity α = sin²θ/2 lands
on the dyadic ladder: 90° (½), 60° (⅜), 45° (¼), 30° (⅛), 0° (0). At 90° we look straight down
the funnel and recover the flat whirlpool; as we tilt toward 0° we see it edge-on, a cone. The
marks are not things in the whirlpool; they are the angles we look from.

Produces: simulations/results/whirlpool/whirlpool_anchor_views.png
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
Z_SCALE = 1.25          # depth = Z_SCALE * |CΨ|: the funnel's height per unit radius

# the anchor angles: camera elevation, with the dyadic polarity α = sin²θ/2 it marks
ANCHORS = [(90, "½"), (60, "⅜"), (45, "¼"), (30, "⅛"), (0, "0")]


def cpsi_magnitude(t: np.ndarray) -> np.ndarray:
    f = np.exp(-4.0 * t)
    return f * (1.0 + f * f) / 6.0


def main() -> None:
    t = np.linspace(0.0, T_MAX, 1500)
    mag = cpsi_magnitude(t)
    cmap = LinearSegmentedColormap.from_list("water",
                                             ["#0a1236", "#163a8a", "#2a8fd0", "#6fe4f5", "#ffffff"])
    norm = Normalize(0.0, START)

    arms = []
    for k in range(N_ARMS):
        phi = 2.0 * np.pi * k / N_ARMS - OMEGA * t
        arms.append((mag * np.cos(phi), mag * np.sin(phi), Z_SCALE * mag))

    th = np.linspace(0.0, 2.0 * np.pi, 360)
    ring = (QUARTER * np.cos(th), QUARTER * np.sin(th), np.full_like(th, Z_SCALE * QUARTER))

    fig = plt.figure(figsize=(20, 5.4))
    fig.patch.set_facecolor(BG)

    for i, (elev, pol) in enumerate(ANCHORS):
        ax = fig.add_subplot(1, len(ANCHORS), i + 1, projection="3d")
        ax.set_facecolor(BG)
        for (x, y, z) in arms:
            pts = np.array([x, y, z]).T.reshape(-1, 1, 3)
            segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
            for lw, alpha in ((3.2, 0.10), (1.3, 0.95)):
                lc = Line3DCollection(segs, cmap=cmap, norm=norm)
                lc.set_array(START - mag[:-1])
                lc.set_linewidth(lw)
                lc.set_alpha(alpha)
                ax.add_collection3d(lc)
        ax.plot(ring[0], ring[1], ring[2], color="#ffd56b", lw=2.0, alpha=0.9)
        ax.scatter([0], [0], [0], s=42, color="#ffffff", alpha=0.5, edgecolors="none")

        ax.view_init(elev=elev, azim=-90)
        ax.set_xlim(-0.33, 0.33)
        ax.set_ylim(-0.33, 0.33)
        ax.set_zlim(0.0, Z_SCALE * START)
        ax.set_box_aspect((1, 1, 0.72))
        ax.set_axis_off()
        ax.text2D(0.5, 0.02, f"{elev}°   ({pol})", transform=ax.transAxes,
                  color="#9fb2d6", fontsize=13, ha="center", style="italic")

    fig.suptitle("the same whirlpool, our viewpoint tilted through the anchor angles   "
                 "(90° straight down the funnel  →  0° edge-on)",
                 color="#cfe0ff", fontsize=14, y=0.96)
    plt.subplots_adjust(left=0.01, right=0.99, top=0.90, bottom=0.02, wspace=0.0)

    out_dir = Path(__file__).parent / "results" / "whirlpool"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "whirlpool_anchor_views.png"
    plt.savefig(out, dpi=150, facecolor=BG)
    plt.close()
    print(f"  saved: {out}")


if __name__ == "__main__":
    main()
