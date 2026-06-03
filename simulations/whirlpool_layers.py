#!/usr/bin/env python3
"""The whirlpool is layers: the funnel sliced into strata, the ¼ cusp one layer among them.

Tilting the top-view whirlpool into 3D (simulations/whirlpool_anchor_views.py) stood it up as a
funnel. Look closer and the depth is layered: the spiral of the complex coherence CΨ descends
through a stack of levels, each a horizontal ring at one value of |CΨ|. The cusp |CΨ| = ¼, the
line where the quantum world meets the classical, is one of those layers, a horizon part-way
down. The drain at the center is the bottom; the mouth |CΨ| = 1/3 is the rim.

Produces: simulations/results/whirlpool/whirlpool_layers.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection
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

# the layers: rings at these |CΨ| levels (the cusp ¼ and the dyadic ⅛ are named)
LAYERS = [1.0 / 3.0, 0.295, QUARTER, 0.205, 0.165, 0.125, 0.092, 0.063, 0.038, 0.018]


def cpsi_magnitude(t: np.ndarray) -> np.ndarray:
    f = np.exp(-4.0 * t)
    return f * (1.0 + f * f) / 6.0


def main() -> None:
    t = np.linspace(0.0, T_MAX, 1500)
    mag = cpsi_magnitude(t)
    cmap = LinearSegmentedColormap.from_list("water",
                                             ["#0a1236", "#163a8a", "#2a8fd0", "#6fe4f5", "#ffffff"])
    norm = Normalize(0.0, START)
    th = np.linspace(0.0, 2.0 * np.pi, 240)

    fig = plt.figure(figsize=(10, 9.2))
    fig.patch.set_facecolor(BG)
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor(BG)

    # the layers: a translucent disk + a bright ring at each level (the strata of the funnel)
    for r in LAYERS:
        z = Z_SCALE * r
        is_cusp = abs(r - QUARTER) < 1e-6
        ring_c = "#ffd56b" if is_cusp else cmap(norm(START - r))
        # the plate (a faint filled disk: the layer you could rest on)
        verts = [list(zip(r * np.cos(th), r * np.sin(th), np.full_like(th, z)))]
        plate = Poly3DCollection(verts, facecolor=("#ffd56b" if is_cusp else "#2a8fd0"),
                                 alpha=0.10 if is_cusp else 0.05, edgecolor="none")
        ax.add_collection3d(plate)
        # the rim of the layer
        ax.plot(r * np.cos(th), r * np.sin(th), np.full_like(th, z),
                color=ring_c, lw=3.0 if is_cusp else 1.4, alpha=0.95 if is_cusp else 0.7)

    # the spiral arms threading down through the layers (dimmed, so the strata read)
    for k in range(N_ARMS):
        phi = 2.0 * np.pi * k / N_ARMS - OMEGA * t
        x, y, z = mag * np.cos(phi), mag * np.sin(phi), Z_SCALE * mag
        pts = np.array([x, y, z]).T.reshape(-1, 1, 3)
        segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
        lc = Line3DCollection(segs, cmap=cmap, norm=norm)
        lc.set_array(START - mag[:-1])
        lc.set_linewidth(1.0)
        lc.set_alpha(0.55)
        ax.add_collection3d(lc)

    ax.scatter([0], [0], [0], s=60, color="#ffffff", alpha=0.6, edgecolors="none")

    # labels on two named layers
    ax.text(QUARTER, 0.0, Z_SCALE * QUARTER, "  ¼  the cusp", color="#ffd56b", fontsize=12,
            ha="left", va="center")
    ax.text(0.125, 0.0, Z_SCALE * 0.125, "  ⅛", color="#9fd0ff", fontsize=11, ha="left", va="center")
    ax.text(0.0, 0.0, Z_SCALE * START + 0.03, "the mouth  1/3", color="#9fb2d6", fontsize=10,
            ha="center", va="bottom")

    ax.view_init(elev=24, azim=-72)
    ax.set_xlim(-0.33, 0.33)
    ax.set_ylim(-0.33, 0.33)
    ax.set_zlim(0.0, Z_SCALE * START)
    ax.set_box_aspect((1, 1, 0.82))
    ax.set_axis_off()

    fig.suptitle("the whirlpool is layers: the funnel sliced into strata, the spiral descending through them,\n"
                 "and the ¼ cusp one horizon-layer among many",
                 color="#cfe0ff", fontsize=13, y=0.95)
    plt.subplots_adjust(left=0.0, right=1.0, top=0.92, bottom=0.0)

    out_dir = Path(__file__).parent / "results" / "whirlpool"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "whirlpool_layers.png"
    plt.savefig(out, dpi=170, facecolor=BG)
    plt.close()
    print(f"  saved: {out}")


if __name__ == "__main__":
    main()
