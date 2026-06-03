#!/usr/bin/env python3
"""The same whirlpool, two payloads: one structure read in two worlds, like a QR code twice.

The abstract qubit and a proton in a hydrogen bond draw the SAME whirlpool, to the pixel,
because the structure is the same: the mirror symmetry, the fold where the quantum gives way to
the classical, the carrier that never stops, are provably substrate-independent. What differs is
only the detail, the payload: what the spiral MEANS, and the clock it keeps. Like a QR code read
twice, one format, different data. The shared frame is the robust, error-correcting part; the
substrate is what you scan off the fine detail.

Produces: simulations/results/whirlpool/whirlpool_two_payloads.png
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
QUARTER = 0.25
N_ARMS = 30
OMEGA = 17.0
T_MAX = 3.2

WATER = LinearSegmentedColormap.from_list("water",
                                          ["#0a1236", "#163a8a", "#2a8fd0", "#6fe4f5", "#ffffff"])


def coherence(t: np.ndarray) -> np.ndarray:
    f = np.exp(-4.0 * t)
    return f * (1.0 + f * f) / 6.0


def draw_vortex(ax, title: str, rim: str, ring_label: str, drain: str, motif: str = "") -> None:
    t = np.linspace(0.0, T_MAX, 1600)
    mag = coherence(t)
    norm = Normalize(0.0, START)
    ax.set_facecolor(BG)
    for k in range(N_ARMS):
        phi = 2.0 * np.pi * k / N_ARMS - OMEGA * t
        x, y = mag * np.cos(phi), mag * np.sin(phi)
        pts = np.array([x, y]).T.reshape(-1, 1, 2)
        segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
        for lw, alpha in ((4.5, 0.06), (2.8, 0.10), (1.3, 0.95)):
            lc = LineCollection(segs, cmap=WATER, norm=norm, capstyle="round")
            lc.set_array(START - mag[:-1])
            lc.set_linewidth(lw)
            lc.set_alpha(alpha)
            ax.add_collection(lc)
    th = np.linspace(0.0, 2.0 * np.pi, 600)
    for lw, alpha in ((8.0, 0.05), (4.5, 0.10), (2.2, 0.55), (1.0, 0.95)):
        ax.plot(QUARTER * np.cos(th), QUARTER * np.sin(th), color="#ffd56b", lw=lw, alpha=alpha)
    ax.scatter([0], [0], s=210, color="#ffffff", alpha=0.18, edgecolors="none")
    ax.scatter([0], [0], s=55, color="#ffffff", alpha=0.35, edgecolors="none")

    ax.annotate(rim, (0.0, 0.305), fontsize=9, color="#bcd0f0", ha="center", va="center")
    ax.annotate(ring_label, (0.252, -0.30), xytext=(0.06, -0.345), fontsize=8.8, color="#ffd56b",
                ha="center", arrowprops=dict(arrowstyle="->", color="#ffd56b", lw=0.8))
    ax.annotate(drain, (0.028, 0.0), xytext=(-0.20, -0.085), fontsize=8.8, color="#eaf2ff",
                ha="center", arrowprops=dict(arrowstyle="->", color="#9fb2d6", lw=0.7))
    if motif:
        ax.text(0.0, -0.335, motif, fontsize=10, color="#9fb2d6", ha="center", va="center")

    lim = 0.36
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, color="#cfe0ff", fontsize=13, pad=6)


def main() -> None:
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(15.5, 8.8))
    fig.patch.set_facecolor(BG)

    draw_vortex(axL, "the abstract qubit",
                "a coherent superposition",
                "the fold: the quantum gives way to the classical",
                "decohered, classical")
    draw_vortex(axR, "a proton in a hydrogen bond",
                "shared across both oxygens",
                "the transfer point: the proton balanced",
                "settled in one well",
                motif="O–H ⋯ O")

    fig.suptitle("The same whirlpool, two payloads: one structure read in two worlds, "
                 "like a QR code scanned twice (one format, different data)",
                 color="#dfeaff", fontsize=14, y=0.975)
    fig.text(0.5, 0.075,
             "the frame, shared and robust:  a fold where the quantum gives way to the classical  ·  "
             "a carrier that never stops  ·  the cusp and its mirror",
             color="#9fd0ff", fontsize=10.5, ha="center")
    fig.text(0.5, 0.038,
             "the payload, the detail that is the substrate:   left, a coherence on a clock you set;   "
             "right, a proton between two oxygens, ticking in femtoseconds, parked at the transfer fold",
             color="#ffce8a", fontsize=10.5, ha="center")
    plt.subplots_adjust(left=0.02, right=0.98, top=0.92, bottom=0.12, wspace=0.04)

    out_dir = Path(__file__).parent / "results" / "whirlpool"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "whirlpool_two_payloads.png"
    plt.savefig(out, dpi=170, facecolor=BG)
    plt.close()
    print(f"  saved: {out}")


if __name__ == "__main__":
    main()
