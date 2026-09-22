"""Plot the exact formulas in PROOF_MISSING_PHASE_SLOW_READOUT.md, section 7.1."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


theta = np.linspace(0, 2*np.pi, 801)
zz = 1 - (1 + np.cos(theta))**2 / 2
triplet = np.sin(theta) * (1 + np.cos(theta)) / (2*np.sqrt(2))
marks = [np.pi/2, 3*np.pi/2]
colors = ["#c65f26", "#27788e"]

with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 12}):
    fig, axes = plt.subplots(2, 1, figsize=(9.2, 7.0), sharex=True,
                             gridspec_kw={"hspace": .28})
    fig.patch.set_facecolor("#faf9f5")
    for ax in axes:
        ax.set_facecolor("#faf9f5")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="#d5d8d9", linewidth=.7)
        for x, color in zip(marks, colors):
            ax.axvline(x, color=color, linewidth=1.1, linestyle="--", alpha=.65)
    axes[0].plot(theta, zz, color="#283d55", linewidth=2.3)
    axes[0].scatter(marks, [.5, .5], color=colors, s=70, zorder=3)
    axes[0].set_title("Two-spin endpoint correlation: same value, opposite slopes", loc="left", fontsize=13)
    axes[0].set_ylabel(r"$\langle Z_0 Z_6\rangle$")
    axes[0].set_ylim(-1.18, 1.18)
    axes[0].annotate("1/2", (marks[0], .5), xytext=(12, 10), textcoords="offset points", color=colors[0])
    axes[0].annotate("1/2", (marks[1], .5), xytext=(12, 10), textcoords="offset points", color=colors[1])
    axes[1].plot(theta, triplet, color="#283d55", linewidth=2.3)
    axes[1].scatter(marks, [np.sqrt(2)/4, -np.sqrt(2)/4], color=colors, s=70, zorder=3)
    axes[1].set_title("Three-spin current correlation: the sign tells them apart", loc="left", fontsize=13)
    axes[1].set_ylabel(r"$\langle Y_0 X_1 Z_3\rangle$")
    axes[1].set_xlabel(r"Motion phase $\theta = 2\sqrt{2}\,t$")
    axes[1].set_ylim(-.56, .56)
    axes[1].set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
                       ["0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])
    axes[1].annotate(r"$+\sqrt{2}/4$", (marks[0], np.sqrt(2)/4), xytext=(12, 8), textcoords="offset points", color=colors[0])
    axes[1].annotate(r"$-\sqrt{2}/4$", (marks[1], -np.sqrt(2)/4), xytext=(12, 15), textcoords="offset points", color=colors[1])
    fig.suptitle("A snapshot and a motion reading", x=.12, ha="left", fontsize=20, fontweight="bold")
    fig.text(.12, .915, "N = 7  |  equal XY bonds  |  J = 1  |  light at the centre only", color="#55616b", fontsize=11)
    fig.text(.12, .025, "At the marked times, all 21 pair density matrices agree (exact partial-trace check).\n"
             r"Time derivative: $d\langle Z_0 Z_6\rangle/dt = 8\langle Y_0 X_1 Z_3\rangle$. Curves show exact formulas.",
             fontsize=10, color="#55616b")
    fig.subplots_adjust(top=.85, bottom=.16, left=.12, right=.96)
    destination = Path(__file__).resolve().parent / "results" / "missing_phase_three_spin_readout.png"
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=180, facecolor=fig.get_facecolor())
    print(destination)
