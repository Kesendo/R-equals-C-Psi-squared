#!/usr/bin/env python3
"""The whirlpool, adapted to water: a proton crossing a hydrogen bond, seen as a spiral.

A hydrogen bond O-H...O is a double well: the proton can sit by the left oxygen (the donor) or
the right one (the acceptor). The picture asks what the proton's coherence would look like if it
drained the way our watched qubits do. It is the same vortex we drew for the abstract coherence,
the same numbers under water's names: the radius is the Bell+ fall |CΨ| = f(1+f²)/6, which
crosses the gold ring at ¼ in one clean pass, and the winding is ours. No proton Hamiltonian is
computed here, and the drawn curve balances no populations. The model in
docs/water/HYDROGEN_BOND_QUBIT.md tells a different story: one coordinate that starts in one well
with no coherence, rises through ¼ within picoseconds in its selected runs and relaxes to the even
mixture; whether a proton would draw this whirlpool is the open question the picture carries.
This is the most
common event in chemistry: a proton crossing a hydrogen bond, which happens countless times a
second in every drop of water and every enzyme.

(Translator's note: this speaks water's language. The framework that draws the spiral stays the
invisible microscope; its names are not on the picture.)

Produces: simulations/results/whirlpool/whirlpool_water.png
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

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BG = "#05060f"
START = 1.0 / 3.0
SHARE = 0.25                  # the ring at ¼, the gold ring
N_ARMS = 30
OMEGA = 17.0
T_MAX = 3.2


def coherence(t: np.ndarray) -> np.ndarray:
    """The drawn radius: the Bell+ fall |CΨ| = f(1+f²)/6 under water's names."""
    f = np.exp(-4.0 * t)
    return f * (1.0 + f * f) / 6.0


def double_well(ax) -> None:
    """A small O-H...O double well with the proton, to ground the picture in the bond."""
    x = np.linspace(-1.6, 1.6, 400)
    v = (x * x - 1.0) ** 2                       # two minima: donor (left), acceptor (right)
    ax.plot(x, v, color="#6fe4f5", lw=1.8, alpha=0.9)
    ax.scatter([-1.0, 1.0], [0.0, 0.0], s=120, color="#ff9d5c", edgecolors="white",
               linewidths=0.6, zorder=4)          # the two oxygens
    ax.text(-1.0, -0.42, "O", color="#ffd0a0", fontsize=11, ha="center")
    ax.text(1.0, -0.42, "O", color="#ffd0a0", fontsize=11, ha="center")
    ax.scatter([-0.7], [0.26], s=46, color="#ffffff", edgecolors="#9fd0ff", linewidths=0.8,
               zorder=5)                           # the proton, in the donor's well, where the model starts
    ax.annotate("H", (-0.7, 0.26), (-0.7, 0.9), color="#ffffff", fontsize=10, ha="center",
                arrowprops=dict(arrowstyle="-", color="#ffffff", lw=0.6))
    ax.annotate("", (1.25, 0.18), (-1.25, 0.18),
                arrowprops=dict(arrowstyle="<->", color="#9fb2d6", lw=0.9, alpha=0.7))
    ax.set_xlim(-1.9, 1.9)
    ax.set_ylim(-0.6, 1.9)
    ax.axis("off")
    ax.set_facecolor(BG)


def main() -> None:
    t = np.linspace(0.0, T_MAX, 2000)
    mag = coherence(t)
    water = LinearSegmentedColormap.from_list("water",
                                              ["#0a1236", "#163a8a", "#2a8fd0", "#6fe4f5", "#ffffff"])
    norm = Normalize(0.0, START)

    fig, ax = plt.subplots(figsize=(9.2, 9.6))
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)

    for k in range(N_ARMS):
        phi = 2.0 * np.pi * k / N_ARMS - OMEGA * t
        x, y = mag * np.cos(phi), mag * np.sin(phi)
        pts = np.array([x, y]).T.reshape(-1, 1, 2)
        segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
        for lw, alpha in ((5.0, 0.05), (3.0, 0.10), (1.4, 0.95)):
            lc = LineCollection(segs, cmap=water, norm=norm, capstyle="round")
            lc.set_array(START - mag[:-1])
            lc.set_linewidth(lw)
            lc.set_alpha(alpha)
            ax.add_collection(lc)

    th = np.linspace(0.0, 2.0 * np.pi, 600)
    for lw, alpha in ((9.0, 0.05), (5.0, 0.10), (2.4, 0.55), (1.1, 0.95)):
        ax.plot(SHARE * np.cos(th), SHARE * np.sin(th), color="#ffd56b", lw=lw, alpha=alpha)
    ax.scatter([0], [0], s=240, color="#ffffff", alpha=0.18, edgecolors="none")
    ax.scatter([0], [0], s=60, color="#ffffff", alpha=0.35, edgecolors="none")

    lim = 0.36
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.axis("off")

    # water-language labels (the tool-names stay off the picture)
    ax.annotate("the mouth: the same numbers,\nwater's names", (0.0, 0.305), fontsize=9.5,
                color="#bcd0f0", ha="center", va="center")
    ax.annotate("the ring at ¼:\nwould a proton meet it here?", (0.18, -0.17),
                xytext=(0.10, -0.345), fontsize=9, color="#ffd56b", ha="center",
                arrowprops=dict(arrowstyle="->", color="#ffd56b", lw=0.8))
    ax.annotate("the drain", (0.028, 0.0), xytext=(-0.205, -0.085),
                fontsize=9, color="#eaf2ff", ha="center",
                arrowprops=dict(arrowstyle="->", color="#9fb2d6", lw=0.7))

    inset = fig.add_axes([0.03, 0.70, 0.24, 0.17])
    double_well(inset)
    inset.set_title("O–H ⋯ O", color="#9fb2d6", fontsize=11, pad=1)

    fig.suptitle("A proton crossing a hydrogen bond, drawn as a whirlpool:\n"
                 "the same numbers under water's names, a question rather than a model",
                 color="#cfe0ff", fontsize=12.5, y=0.985)
    plt.subplots_adjust(left=0.02, right=0.98, top=0.90, bottom=0.02)

    out_dir = Path(__file__).parent.parent / "results" / "whirlpool"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "whirlpool_water.png"
    plt.savefig(out, dpi=180, facecolor=BG, bbox_inches="tight", pad_inches=0.12)
    plt.close()
    print(f"  saved: {out}")


if __name__ == "__main__":
    main()
