"""Sector fingerprint: colour every mode by its disagreement-count sector, globally.

Tom's idea: colour by sector could be THE fingerprint, a global RGB shared by all systems.

The gift: we need NO eigenvectors. The absorption law holds exactly for every eigenmode,
    Re(λ) = −2γ·⟨n_XY⟩,
(because L r = λ r gives Re(λ)·‖r‖² = ⟨r|½(L+L†)|r⟩ = ⟨r|(−2γ·n_XY diagonal)|r⟩, the
Hamiltonian part being anti-Hermitian and contributing only to Im). So a mode's sector is
just ⟨n_XY⟩ = −Re(λ)/(2γ), read straight off the eigenvalue we already exported. The data
confirm it: max|Re| = 0.5 at N=5 and 0.6 at N=6, exactly n_XY = N.

So the horizontal axis IS the sector axis. Here we colour every mode of every molecule by
that one global coordinate (the disagreement count the light reads). It is the same scheme
for all systems, a global RGB, and it tests Tom's "the horizontal rectangles are sectors":
single-colour rectangles would BE sectors; rainbow ones would cross sectors.

Live C# Symphony exports, Q=1.5. No physics computed here.
  python simulations/reel_sector.py
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
})

GAMMA = 0.05
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "results", "symphony_reel")


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


systems = []
root_csv = os.path.join(ROOT, "symphony_eigenvalues.csv")
if os.path.exists(root_csv):
    systems.append(("chain", root_csv))
for d in sorted(os.listdir(ROOT)):
    p = os.path.join(ROOT, d, "symphony_eigenvalues.csv")
    if os.path.isdir(os.path.join(ROOT, d)) and os.path.exists(p):
        systems.append((d, p))

allx, ally, allsec = [], [], []
for name, path in systems:
    meta, re, im = read(path)
    sigma = meta.get("sigma", meta.get("N", 5) * GAMMA)
    sector = -re / (2.0 * GAMMA)                 # <n_XY> exactly, via the absorption law
    allx.append(re + sigma); ally.append(im); allsec.append(sector)
    print(f"{name:24s} sector range [{sector.min():.2f}, {sector.max():.2f}]  (N={meta.get('N')})")
x = np.concatenate(allx); y = np.concatenate(ally); sec = np.concatenate(allsec)

# Draw faded points first, brighter on top, coloured by sector (the global coordinate).
order = np.argsort(np.abs(y))[::-1]              # off-axis first, real-axis on top
fig, ax = plt.subplots(figsize=(10, 8.4))
scat = ax.scatter(x[order], y[order], c=sec[order], s=7, cmap="turbo",
                  alpha=0.7, edgecolors="none", vmin=0, vmax=6)
cb = fig.colorbar(scat, ax=ax, pad=0.01)
cb.set_label("sector  ⟨n_XY⟩ = −Re λ / 2γ  (the disagreement count the light reads)", color=FG)
cb.ax.yaxis.set_tick_params(color=FG)
plt.setp(plt.getp(cb.ax.axes, "yticklabels"), color=FG)
ax.set_xlabel(r"Re $\lambda + \sigma$  (centred; this axis IS the sector axis)")
ax.set_ylabel(r"Im $\lambda$  (the Hamiltonian's turn)")
ax.set_title("SECTOR FINGERPRINT: every mode of every molecule, coloured by its global sector\n"
             "(six systems superimposed; sector = the disagreement count, free from the absorption law)")
fig.tight_layout()
OUT = os.path.join(ROOT, "sector_fingerprint.png")
fig.savefig(OUT, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"\nwrote: {OUT}")
