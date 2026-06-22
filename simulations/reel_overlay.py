"""Overlay the reels: three molecular Liouvillian spectra on one canvas.

The companion to reel_and_projector.py. Where that drew one system, this LAYS THREE
OVER EACH OTHER (Tom's "uebereinander legen"): the time-erased spectrum, the QR code, of
three systems superimposed, so the shared skeleton (the F1 palindrome, the absorption
rungs, the one object) shows through, and where each molecule departs from it stands out.

Data is exported by the live C# lab (the Symphony witness, `inspect --root symphony
--export --export-name <name>`), all at the same operating point Q=1.5, then read here.
No physics is computed (the deprecated framework is not imported); cyberpunk palette.

  Feed the lab first (already done if the subdirs exist):
    inspect --root symphony --N 5 --J 0.075 --gamma 0.05 --export                          (chain, XY)
    inspect --root symphony --N 5 --J 0.075 --gamma 0.05 --htype heisenberg \
            --export --export-name water_heisenberg                                        (water)
    inspect --root symphony --N 6 --topology ring --htype xy --J 0.075 --gamma 0.05 \
            --export --export-name benzene_xyring                                          (benzene)
  Then:
    python simulations/reel_overlay.py
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
    "text.color": FG, "xtick.color": FG, "ytick.color": FG, "grid.color": GRID,
    "grid.alpha": 0.35, "font.family": "monospace", "legend.facecolor": PANEL,
    "legend.edgecolor": GRID, "axes.grid": True,
})


def glow_scatter(ax, x, y, color, s=10, alpha=0.85, zorder=3, layers=4, spread=6.0):
    for i in range(layers, 0, -1):
        ax.scatter(x, y, s=s * (1 + spread * i / layers), color=color, alpha=0.04,
                   edgecolors="none", zorder=zorder - 0.1)
    ax.scatter(x, y, s=s, color=color, alpha=alpha, edgecolors="none", zorder=zorder)


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "results", "symphony_reel")

# (label, subdir-or-None-for-root, neon color)
SYSTEMS = [
    ("chain  (XY, N=5)",            None,                "#08f7fe"),   # cyan
    ("water  (Heisenberg chain N=5)", "water_heisenberg", "#fe53bb"),  # pink
    ("benzene (XY ring, N=6)",      "benzene_xyring",    "#39ff14"),   # green
]


def read(subdir):
    path = os.path.join(ROOT if subdir is None else os.path.join(ROOT, subdir),
                        "symphony_eigenvalues.csv")
    meta, ri = {}, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                for tok in line[1:].split():
                    if "=" in tok:
                        k, v = tok.split("=", 1)
                        meta[k] = float(v)
                continue
            if line.startswith("Re"):
                continue
            r, i = line.split(",")
            ri.append((float(r), float(i)))
    a = np.array(ri)
    return meta, a[:, 0], a[:, 1]


loaded = []
for label, subdir, color in SYSTEMS:
    meta, re, im = read(subdir)
    sigma = meta.get("sigma", meta["N"] * meta["gamma"])
    loaded.append((label, color, re, im, sigma, int(meta.get("N", 0))))
    print(f"{label:32s} modes={len(re):5d}  sigma={sigma:.3g}  "
          f"Re[{re.min():.3f},{re.max():.3f}]  |Im|max={np.abs(im).max():.3f}")

OUT = os.path.join(ROOT, "overlay_three_systems.png")
fig, (axA, axB) = plt.subplots(1, 2, figsize=(15.5, 6.6))

# Panel A: raw, all at Q=1.5 on the true (Re, Im) plane. Heaviest point-count drawn first.
for label, color, re, im, sigma, N in sorted(loaded, key=lambda t: -len(t[2])):
    glow_scatter(axA, re, im, color, s=7, alpha=0.55)
for label, color, re, im, sigma, N in loaded:
    axA.axvline(-sigma, color=color, ls="--", lw=1.0, alpha=0.7)
    axA.scatter([], [], color=color, s=22, label=f"{label}  (−σ={-sigma:.2g})")
axA.set_xlabel(r"Re $\lambda$  (fade rate; γ shared, so rungs $-2\gamma n$ coincide)")
axA.set_ylabel(r"Im $\lambda$  (turn frequency)")
axA.set_title("raw: three systems at Q = 1.5, one canvas")
axA.legend(loc="upper left", fontsize=8.5, labelcolor=FG)

# Panel B: centered on −σ (Re + σ), so every palindrome center sits at 0 and the shared
# skeleton overlays; the molecular departures stand out.
for label, color, re, im, sigma, N in sorted(loaded, key=lambda t: -len(t[2])):
    glow_scatter(axB, re + sigma, im, color, s=7, alpha=0.6)
for label, color, re, im, sigma, N in loaded:
    axB.scatter([], [], color=color, s=22, label=label)
axB.axvline(0.0, color="#39ff14", ls="--", lw=1.2, alpha=0.85, label="shared center (−σ → 0)")
axB.set_xlabel(r"Re $\lambda + \sigma$  (centered on the palindrome)")
axB.set_ylabel(r"Im $\lambda$")
axB.set_title("centered on −σ: the shared skeleton, the one object")
axB.legend(loc="upper left", fontsize=8.5, labelcolor=FG)

fig.suptitle("THE REELS, OVERLAID: three molecular Liouvillian spectra superimposed (live C# Symphony, Q=1.5)",
             color="#39ff14", fontsize=13, y=1.02)
fig.tight_layout()
fig.savefig(OUT, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"\nwrote: {OUT}")
