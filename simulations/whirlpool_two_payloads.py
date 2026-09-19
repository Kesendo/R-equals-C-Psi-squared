#!/usr/bin/env python3
"""One drawing recipe placed beside two deliberately different readings.

The left panel uses the exact F25 Bell+ local-Z-dephasing radius
``r=f(1+f^2)/6``. Its ``r=1/4`` crossing is unique and entangled because
concurrence is ``C=f>0``; ``f=1/2`` gives the below-quarter control
``r=5/48`` and ``C=1/2``. The phase, arms and radial ring remain drawing
choices, not an F25 geodesic. The ring marks no horizon and is not the
cardioid boundary.

The right panel asks how this imagery might look in the language of a proton
and two oxygens. Reusing pixels and arrays proves only code reuse. There is no
proton Hamiltonian, bath model, population dynamics, time conversion, or
measurement here, so the illustration cannot establish equal sharing,
localization in a well, or physics that transfers between substrates. The
recurrence cusp and specified F86 toy EP are separate mathematical
resemblances, not identities with either panel.

Produces exactly ``whirlpool_two_payloads.png`` in the caller-selected output
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
WATER_COLORS = ("#0a1236", "#163a8a", "#2a8fd0", "#6fe4f5", "#ffffff")


def coherence(t: np.ndarray) -> np.ndarray:
    """Return the F25 Bell+ scalar radius at gamma=1."""
    f = np.exp(-4.0 * t)
    return f * (1.0 + f * f) / 6.0


def draw_vortex(ax, title: str, rim: str, ring_label: str, center_label: str) -> None:
    """Draw the shared visual scaffold without transferring physical meaning."""
    t = np.linspace(0.0, T_MAX, 1600)
    radius = coherence(t)
    water = LinearSegmentedColormap.from_list("water", WATER_COLORS)
    norm = Normalize(0.0, START)
    ax.set_facecolor(BG)
    for arm in range(N_ARMS):
        phase = 2.0 * np.pi * arm / N_ARMS - OMEGA * t
        x = radius * np.cos(phase)
        y = radius * np.sin(phase)
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        for width, alpha in ((4.5, 0.06), (2.8, 0.10), (1.3, 0.95)):
            line = LineCollection(segments, cmap=water, norm=norm, capstyle="round")
            line.set_array(START - radius[:-1])
            line.set_linewidth(width)
            line.set_alpha(alpha)
            ax.add_collection(line)
    angle = np.linspace(0.0, 2.0 * np.pi, 600)
    for width, alpha in ((8.0, 0.05), (4.5, 0.10), (2.2, 0.55), (1.0, 0.95)):
        ax.plot(QUARTER * np.cos(angle), QUARTER * np.sin(angle),
                color="#ffd56b", lw=width, alpha=alpha)
    ax.scatter([0], [0], s=210, color="#ffffff", alpha=0.18, edgecolors="none")
    ax.scatter([0], [0], s=55, color="#ffffff", alpha=0.35, edgecolors="none")
    ax.annotate(rim, (0.0, 0.305), fontsize=9, color="#bcd0f0",
                ha="center", va="center")
    ax.annotate(ring_label, (0.252, -0.30), xytext=(0.06, -0.345),
                fontsize=8.8, color="#ffd56b", ha="center",
                arrowprops=dict(arrowstyle="->", color="#ffd56b", lw=0.8))
    ax.annotate(center_label, (0.028, 0.0), xytext=(-0.20, -0.085),
                fontsize=8.8, color="#eaf2ff", ha="center",
                arrowprops=dict(arrowstyle="->", color="#9fb2d6", lw=0.7))
    ax.set_xlim(-0.36, 0.36)
    ax.set_ylim(-0.36, 0.36)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, color="#cfe0ff", fontsize=13, pad=6)


def main(output_dir) -> None:
    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)
    fig, (left, right) = plt.subplots(1, 2, figsize=(15.5, 9.2))
    fig.patch.set_facecolor(BG)
    draw_vortex(left, "abstract Bell+ model drawing",
                "F25 radius at t=0", "r=1/4 radial reference",
                "finite-time limit approaches zero")
    draw_vortex(right, "hydrogen-bond imagery: a question",
                "borrowed visual scaffold", "not a population-balance claim",
                "not a localization result")
    right.text(0.0, -0.315, "O-H ... O", fontsize=10, color="#9fb2d6",
               ha="center", va="center")
    fig.suptitle("one drawing recipe, two readings kept visibly apart",
                 color="#dfeaff", fontsize=14, y=0.975)
    fig.text(0.5, 0.080,
             "illustrative translation; shared arrays do not prove substrate physics",
             color="#9fd0ff", fontsize=10.5, ha="center")
    fig.text(0.5, 0.050,
             "no proton Hamiltonian, bath model, or measurement",
             color="#ffce8a", fontsize=10.5, ha="center")
    fig.text(0.5, 0.020,
             "F25 controls: unique r=1/4 at f*=0.8612240997395736 with C=f>0; "
             "f=1/2 gives r=5/48<1/4 and C=1/2",
             color="#8298bd", fontsize=9, ha="center")
    plt.subplots_adjust(left=0.02, right=0.98, top=0.91, bottom=0.12, wspace=0.04)
    fig.savefig(output_root / "whirlpool_two_payloads.png", dpi=170, facecolor=BG)
    plt.close(fig)
    print("saved: whirlpool_two_payloads.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Render the two-reading Whirlpool illustration.")
    parser.add_argument("--output-dir", required=True)
    arguments = parser.parse_args()
    main(arguments.output_dir)
