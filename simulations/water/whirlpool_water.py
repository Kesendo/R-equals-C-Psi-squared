#!/usr/bin/env python3
"""Water-language imagery built from the F25 Bell+ drawing curve.

The plotted radius is the abstract local-Z-dephasing Bell+ scalar
``r=f(1+f^2)/6`` with ``f=exp(-4 gamma t)``. Its ``r=1/4`` crossing is
unique and entangled (concurrence ``C=f>0``); the control ``f=1/2`` has
``r=5/48<1/4`` and ``C=1/2``. The spiral phase, arms, radial ring and inset
double well are illustrative drawing choices, not an F25 geodesic. The ring
marks no physical horizon and is not the cardioid boundary.

The O-H...O panel asks a water-flavoured question. This script contains no
proton Hamiltonian, bath model, populations, time conversion, or measurement.
No proton decay law is computed, and no population equality or localization
result follows. The recurrence cusp and
specified F86 toy EP remain separate mathematical resemblances.

Produces exactly ``whirlpool_water.png`` in the caller-selected output
directory.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np

BG = "#05060f"
START = 1.0 / 3.0
QUARTER = 0.25
N_ARMS = 30
OMEGA = 17.0
T_MAX = 3.2


def coherence(t: np.ndarray) -> np.ndarray:
    """Return the F25 Bell+ scalar radius at gamma=1, not proton dynamics."""
    f = np.exp(-4.0 * t)
    return f * (1.0 + f * f) / 6.0


def double_well(ax) -> None:
    """Draw a schematic quartic well; it is not a fitted material potential."""
    x = np.linspace(-1.6, 1.6, 400)
    potential = (x * x - 1.0) ** 2
    ax.plot(x, potential, color="#6fe4f5", lw=1.8, alpha=0.9)
    ax.scatter([-1.0, 1.0], [0.0, 0.0], s=120, color="#ff9d5c",
               edgecolors="white", linewidths=0.6, zorder=4)
    ax.text(-1.0, -0.42, "O", color="#ffd0a0", fontsize=11, ha="center")
    ax.text(1.0, -0.42, "O", color="#ffd0a0", fontsize=11, ha="center")
    ax.scatter([0.0], [1.02], s=46, color="#ffffff", edgecolors="#9fd0ff",
               linewidths=0.8, zorder=5)
    ax.annotate("H (drawn)", (0.0, 1.02), (0.0, 1.5), color="#ffffff",
                fontsize=9, ha="center",
                arrowprops=dict(arrowstyle="-", color="#ffffff", lw=0.6))
    ax.set_xlim(-1.9, 1.9)
    ax.set_ylim(-0.6, 1.9)
    ax.axis("off")
    ax.set_facecolor(BG)


def main(output_dir) -> None:
    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)
    t = np.linspace(0.0, T_MAX, 2000)
    radius = coherence(t)
    water = LinearSegmentedColormap.from_list(
        "water", ["#0a1236", "#163a8a", "#2a8fd0", "#6fe4f5", "#ffffff"])
    norm = Normalize(0.0, START)
    fig, ax = plt.subplots(figsize=(9.2, 10.0))
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)

    for arm in range(N_ARMS):
        phase = 2.0 * np.pi * arm / N_ARMS - OMEGA * t
        x = radius * np.cos(phase)
        y = radius * np.sin(phase)
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        for width, alpha in ((5.0, 0.05), (3.0, 0.10), (1.4, 0.95)):
            line = LineCollection(segments, cmap=water, norm=norm, capstyle="round")
            line.set_array(START - radius[:-1])
            line.set_linewidth(width)
            line.set_alpha(alpha)
            ax.add_collection(line)

    angle = np.linspace(0.0, 2.0 * np.pi, 600)
    for width, alpha in ((9.0, 0.05), (5.0, 0.10), (2.4, 0.55), (1.1, 0.95)):
        ax.plot(QUARTER * np.cos(angle), QUARTER * np.sin(angle),
                color="#ffd56b", lw=width, alpha=alpha)
    ax.scatter([0], [0], s=240, color="#ffffff", alpha=0.18, edgecolors="none")
    ax.scatter([0], [0], s=60, color="#ffffff", alpha=0.35, edgecolors="none")
    ax.set_xlim(-0.36, 0.36)
    ax.set_ylim(-0.36, 0.36)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.annotate("borrowed Bell+ drawing curve", (0.0, 0.305), fontsize=9.5,
                color="#bcd0f0", ha="center", va="center")
    ax.annotate("r=1/4 radial reference\nnot a proton population balance",
                (0.255, -0.30), xytext=(0.10, -0.345), fontsize=9,
                color="#ffd56b", ha="center",
                arrowprops=dict(arrowstyle="->", color="#ffd56b", lw=0.8))
    ax.annotate("drawing limit\nnot a localization result", (0.028, 0.0),
                xytext=(-0.205, -0.085), fontsize=9, color="#eaf2ff",
                ha="center",
                arrowprops=dict(arrowstyle="->", color="#9fb2d6", lw=0.7))

    inset = fig.add_axes([0.04, 0.76, 0.27, 0.19])
    double_well(inset)
    inset.set_title("schematic O-H ... O", color="#9fb2d6", fontsize=10, pad=1)
    fig.suptitle("water imagery around an abstract Bell+ radius",
                 color="#cfe0ff", fontsize=12.5, y=0.985)
    fig.text(0.5, 0.075, "water imagery: illustrative O-H...O question",
             color="#cfe0ff", fontsize=10.5, ha="center")
    fig.text(0.5, 0.047, "no proton Hamiltonian, bath model, or measurement",
             color="#ffce8a", fontsize=10, ha="center")
    fig.text(0.5, 0.020,
             "F25 controls: unique r=1/4 at f*=0.8612240997395736 with C=f>0; "
             "f=1/2 gives r=5/48<1/4 and C=1/2",
             color="#8298bd", fontsize=9, ha="center")
    plt.subplots_adjust(left=0.02, right=0.98, top=0.92, bottom=0.11)
    fig.savefig(output_root / "whirlpool_water.png", dpi=180, facecolor=BG,
                bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)
    print("saved: whirlpool_water.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Render the water-language Whirlpool illustration.")
    parser.add_argument("--output-dir", required=True)
    arguments = parser.parse_args()
    main(arguments.output_dir)
