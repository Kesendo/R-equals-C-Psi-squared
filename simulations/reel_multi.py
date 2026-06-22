"""Many reels overlaid: does the XOR pattern stay the same as more molecules join?

Tom's conjecture: throw in more chemical systems and they form the SAME pattern as the
inner three (the self-similar bands the eye reads). This is the test. We auto-discover
every system exported under results/symphony_reel/ (root = chain, plus each subdir), and
overlay all their centred spectra in distinct neon colours. If the conjecture holds, the
shared structure (the palindrome symmetry, the popcount-sector bands) is reinforced, not
broken, as systems are added; a new system that breaks it would refute the universality.

All systems are live C# Symphony exports at Q=1.5 (J=0.075, γ=0.05). Centred on each
system's own −σ. No physics computed here.

  python simulations/reel_multi.py
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BG, PANEL, FG, GRID = "#080b12", "#0c1018", "#7ef9ff", "#15303a"
plt.rcParams.update({
    "figure.facecolor": BG, "axes.facecolor": PANEL, "savefig.facecolor": BG,
    "axes.edgecolor": FG, "axes.labelcolor": FG, "axes.titlecolor": "#39ff14",
    "text.color": FG, "xtick.color": FG, "ytick.color": FG, "font.family": "monospace",
    "legend.facecolor": PANEL, "legend.edgecolor": GRID,
})

NEON = ["#08f7fe", "#fe53bb", "#39ff14", "#ffe700", "#b14aed", "#ff6a00",
        "#00ff9f", "#ff2d55"]

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "results", "symphony_reel")

# Friendly labels for known export names; anything else uses the dir name.
LABELS = {
    "_root_": "chain (XY, N5)",
    "water_heisenberg": "water (Heis chain, N5)",
    "benzene_xyring": "benzene (XY ring, N6)",
    "butadiene_xy_chain": "butadiene (XY chain, N4)",
    "cyclobutadiene_xy_ring": "cyclobutadiene (XY ring, N4)",
    "hexatriene_xy_chain": "hexatriene (XY chain, N6)",
}


def read(path):
    meta, ri = {}, []
    for ln in open(path, encoding="utf-8"):
        ln = ln.strip()
        if ln.startswith("#"):
            for t in ln[1:].split():
                if "=" in t:
                    k, v = t.split("=", 1); meta[k] = float(v)
            continue
        if not ln or ln.startswith("Re"):
            continue
        a, b = ln.split(","); ri.append((float(a), float(b)))
    a = np.array(ri)
    return meta, a[:, 0], a[:, 1]


# Discover systems: root first, then subdirs alphabetically.
systems = []
root_csv = os.path.join(ROOT, "symphony_eigenvalues.csv")
if os.path.exists(root_csv):
    systems.append(("_root_", root_csv))
for d in sorted(os.listdir(ROOT)):
    p = os.path.join(ROOT, d, "symphony_eigenvalues.csv")
    if os.path.isdir(os.path.join(ROOT, d)) and os.path.exists(p):
        systems.append((d, p))

fig, ax = plt.subplots(figsize=(9.5, 8.6))
ax.axvline(0.0, color=GRID, lw=0.8, alpha=0.6)
ax.axhline(0.0, color=GRID, lw=0.8, alpha=0.6)

loaded = []
for i, (name, path) in enumerate(systems):
    meta, re, im = read(path)
    sigma = meta.get("sigma", meta.get("N", 5) * meta.get("gamma", 0.05))
    loaded.append((name, re, im, sigma, len(re)))

# Draw heaviest first so smaller sets stay visible on top.
for name, re, im, sigma, n in sorted(loaded, key=lambda t: -t[4]):
    i = [s[0] for s in systems].index(name)
    col = NEON[i % len(NEON)]
    ax.scatter(re + sigma, im, s=5, color=col, alpha=0.45, edgecolors="none")
for i, (name, path) in enumerate(systems):
    col = NEON[i % len(NEON)]
    lab = LABELS.get(name, name)
    ax.scatter([], [], color=col, s=22, label=lab)
    print(f"{lab:32s} modes={dict((x[0],x[4]) for x in loaded)[name]:5d}")

ax.set_aspect("auto")
ax.set_xlabel(r"Re $\lambda + \sigma$  (centred on the palindrome)")
ax.set_ylabel(r"Im $\lambda$")
ax.set_title(f"MANY REELS OVERLAID: {len(systems)} chemical systems superimposed (live C# Symphony, Q=1.5)\n"
             "does the pattern hold as more molecules join?")
ax.legend(loc="upper left", fontsize=8, labelcolor=FG, ncol=2)
fig.tight_layout()
OUT = os.path.join(ROOT, f"multi_{len(systems)}_systems.png")
fig.savefig(OUT, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"\nwrote: {OUT}")
