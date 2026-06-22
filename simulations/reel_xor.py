"""The XOR of the reels: additive RGB density of three molecular spectra.

The second image (Tom's idea). The overlay (reel_overlay.py) draws the three spectra as
points, but a point plot has a z-order: something is drawn "on top", and "who is first?"
is a question with no real answer. So here we DROP the order and blend with transparency:
each system is one additive colour CHANNEL on a binned density of the centred spectrum.

  * where all three coincide, the channels add to WHITE: the shared skeleton, the one
    object, where the molecules agree;
  * where they part, a channel survives alone (or two mix): the colour that remains IS the
    XOR, the symmetric difference, where each molecule departs from the shared object.

XOR is the dephasing's native operation: the light content the environment reads is
popcount(i XOR j) (see experiments/XOR_SPACE.md), so asking "what is the XOR of these
reels" lands on the physics, not just the picture.

Caveat (honest): the binning resolution is a representation choice; the qualitative split
(white shared core, coloured difference fringe) is robust, the exact pixels are not.

Data: the live C# Symphony export at Q=1.5 (same three systems as reel_overlay.py).
  python simulations/reel_xor.py
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


BG, FG = "#080b12", "#7ef9ff"
plt.rcParams.update({
    "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
    "axes.edgecolor": FG, "axes.labelcolor": FG, "axes.titlecolor": "#39ff14",
    "text.color": FG, "xtick.color": FG, "ytick.color": FG, "font.family": "monospace",
})

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "results", "symphony_reel")

# (label, subdir-or-None, RGB channel index 0=R 1=G 2=B)
SYSTEMS = [
    ("chain (XY, N=5)               -> red",   None,               0),
    ("water (Heisenberg chain N=5)  -> green", "water_heisenberg", 1),
    ("benzene (XY ring, N=6)        -> blue",  "benzene_xyring",   2),
]


def read(subdir):
    path = os.path.join(ROOT if subdir is None else os.path.join(ROOT, subdir),
                        "symphony_eigenvalues.csv")
    meta, ri = {}, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("Re"):
                continue
            if line.startswith("#"):
                for tok in line[1:].split():
                    if "=" in tok:
                        k, v = tok.split("=", 1)
                        meta[k] = float(v)
                continue
            r, i = line.split(",")
            ri.append((float(r), float(i)))
    a = np.array(ri)
    return meta, a[:, 0], a[:, 1]


# Common centred grid (Re + sigma so every palindrome centre sits at 0).
NX, NY = 240, 320
RE_LIM, IM_LIM = 0.34, 0.60
xedges = np.linspace(-RE_LIM, RE_LIM, NX + 1)
yedges = np.linspace(-IM_LIM, IM_LIM, NY + 1)

rgb = np.zeros((NY, NX, 3), dtype=float)
for label, subdir, ch in SYSTEMS:
    meta, re, im = read(subdir)
    sigma = meta.get("sigma", meta["N"] * meta["gamma"])
    H, _, _ = np.histogram2d(im, re + sigma, bins=[yedges, xedges])   # rows=Im, cols=Re
    dens = np.sqrt(H)                                                 # sqrt for visibility
    hi = np.percentile(dens[dens > 0], 99.0) if np.any(dens > 0) else 1.0
    rgb[:, :, ch] = np.clip(dens / hi, 0.0, 1.0)
    print(f"{label:42s} modes={len(re):5d} sigma={sigma:.3g} maxcount={int(H.max())}")

# A faint glow: blend in a lightly blurred copy so single points still register.
def blur(a, k=1):
    out = a.copy()
    for _ in range(k):
        out = (out
               + np.roll(out, 1, 0) + np.roll(out, -1, 0)
               + np.roll(out, 1, 1) + np.roll(out, -1, 1)) / 5.0
    return out

rgb_show = np.clip(0.78 * rgb + 0.55 * blur(rgb, 2), 0.0, 1.0)

fig, ax = plt.subplots(figsize=(9.2, 8.4))
ax.imshow(rgb_show, origin="lower", extent=[-RE_LIM, RE_LIM, -IM_LIM, IM_LIM],
          aspect="auto", interpolation="nearest")
ax.axvline(0.0, color="#ffffff", ls=":", lw=0.8, alpha=0.5)
ax.set_xlabel(r"Re $\lambda + \sigma$  (centred on the palindrome)")
ax.set_ylabel(r"Im $\lambda$")
ax.set_title("THE XOR OF THE REELS\n"
             "additive: no system is 'first'. WHITE = all three agree (the one object); "
             "colour = where a molecule departs (the XOR)")
# legend as coloured text
for i, (label, _, ch) in enumerate(SYSTEMS):
    col = [0, 0, 0]; col[ch] = 1.0
    ax.text(0.015, 0.97 - 0.045 * i, label, transform=ax.transAxes,
            color=tuple(col), fontsize=9, va="top",
            bbox=dict(facecolor="#0c1018", edgecolor="none", alpha=0.6, pad=2))
ax.text(0.015, 0.97 - 0.045 * 3, "white = R+G+B = shared skeleton", transform=ax.transAxes,
        color="#ffffff", fontsize=9, va="top",
        bbox=dict(facecolor="#0c1018", edgecolor="none", alpha=0.6, pad=2))

fig.tight_layout()
OUT = os.path.join(ROOT, "xor_three_systems.png")
fig.savefig(OUT, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"\nwrote: {OUT}")
