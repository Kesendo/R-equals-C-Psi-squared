"""Draw MirrorWorld's Atmosphere: its price over the block lattice and its dark classes.

Every number is read off the C# object (compute/MirrorWorld/Atmosphere.cs) through its gallery export;
nothing physical is computed here, only drawn. The readings in the panel titles ("the GHZ pair",
"the Néel pair") are checked against the exported classes before anything is drawn.

  dotnet run --project compute/MirrorWorld -- atmosphere gallery simulations/results/atmosphere_gallery.json
  python simulations/atmosphere_gallery.py

Writes visualizations/atmosphere_price_landscape.png and visualizations/atmosphere_dark_classes.png.
"""
from __future__ import annotations

import json
import os
from math import comb

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgba

BG, PANEL, FG, GRID = "#080b12", "#0c1018", "#7ef9ff", "#15303a"
CYAN, PINK, GREEN, AMBER, VIOLET = "#08f7fe", "#fe53bb", "#39ff14", "#ffe700", "#b14aed"
plt.rcParams.update({
    "figure.facecolor": BG, "axes.facecolor": PANEL, "savefig.facecolor": BG,
    "axes.edgecolor": FG, "axes.labelcolor": FG, "axes.titlecolor": GREEN,
    "text.color": FG, "xtick.color": FG, "ytick.color": FG, "font.family": "monospace",
})

DEPTH = LinearSegmentedColormap.from_list("depth", [PANEL, "#3b1f6e", VIOLET, PINK, AMBER])
LIFE = LinearSegmentedColormap.from_list("life", [PANEL, "#1a1440", "#4a2a8a", PINK])

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "results", "atmosphere_gallery.json")
OUT = os.path.join(HERE, "..", "visualizations")


def glow_line(ax, x, y, color, lw=1.2, zorder=4, layers=6, spread=7.0):
    for i in range(layers, 0, -1):
        ax.plot(x, y, color=color, lw=lw + spread * i / layers, alpha=0.05, solid_capstyle="round", zorder=zorder - 0.1)
    ax.plot(x, y, color=color, lw=lw, solid_capstyle="round", zorder=zorder)


def bits(i, n):
    """A configuration as a string, site 0 first."""
    return "".join(str((i >> k) & 1) for k in range(n))


def bowl_panel(ax, panel, vmax, orbit=False):
    n = panel["N"]
    depth = {(p, q): price / 2 for p, q, price in panel["prices"]}
    for p in range(n + 1):
        for q in range(n + 1):
            x, y = q, n - p                                   # q to the right, p downward
            if (p, q) in depth:
                v = depth[(p, q)]
                ax.add_patch(plt.Rectangle((x - 0.46, y - 0.46), 0.92, 0.92, fc=DEPTH(v / vmax), ec=CYAN, lw=0.6, zorder=2))
                ax.text(x, y, f"{v:g}", ha="center", va="center", fontsize=6.5, zorder=6, color=BG if v / vmax > 0.55 else FG)
            else:
                ax.add_patch(plt.Rectangle((x - 0.46, y - 0.46), 0.92, 0.92, fill=False, ec="#2f6270", lw=0.9,
                                           ls=(0, (2, 2)), zorder=2))
    c = n / 2
    if len(depth) == (n + 1) ** 2:
        radii = sorted({(p - c) ** 2 + (q - c) ** 2 for p in range(n + 1) for q in range(n + 1)})
        t = np.linspace(0, 2 * np.pi, 400)
        for r2 in radii[1:]:
            ax.plot(c + np.sqrt(r2) * np.cos(t), c + np.sqrt(r2) * np.sin(t), color=CYAN, lw=0.5, alpha=0.22, zorder=3)
    if orbit:
        # one quarter-turn orbit, (p, q) -> (q, N - p), in the panel's frame x = q, y = N - p
        p, q, pts = 0, 1, []
        for _ in range(4):
            pts.append((q, n - p))
            p, q = q, n - p
        for x, y in pts:
            ax.add_patch(plt.Rectangle((x - 0.46, y - 0.46), 0.92, 0.92, fill=False, ec=PINK, lw=1.6, zorder=5))
        pts.append(pts[0])
        for (x1, y1), (x2, y2) in zip(pts, pts[1:]):          # each side from tile edge to tile edge
            dx, dy = x2 - x1, y2 - y1
            t = 0.5 / max(abs(dx), abs(dy))
            glow_line(ax, [x1 + t * dx, x2 - t * dx], [y1 + t * dy, y2 - t * dy], PINK, lw=0.8, layers=4, spread=4.0)
    ax.set_xlim(-0.6, n + 0.6); ax.set_ylim(-0.6, n + 0.6); ax.set_aspect("equal")
    ax.set_xticks(range(n + 1)); ax.set_yticks(range(n + 1)); ax.set_yticklabels([str(n - k) for k in range(n + 1)])
    ax.tick_params(labelsize=7, length=0)
    ax.set_xlabel("q (bra popcount)", fontsize=7.5); ax.set_ylabel("p (ket popcount)", fontsize=7.5)
    for s in ax.spines.values():
        s.set_color(GRID)


def landscape(data, out):
    bowl = {panel["name"]: panel for panel in data["bowl"]}
    anti, flat = bowl["anti-correlated end"], bowl["featureless"]
    assert all(anti["covariance"][k][l] * (anti["N"] - 1) == -anti["gamma"]
               for k in range(anti["N"]) for l in range(anti["N"]) if k != l), "the anti-correlated end"
    assert not any(any(row) for row in flat["covariance"]), "featureless"
    n, gamma = bowl["shared field"]["N"], bowl["shared field"]["gamma"]
    vmax = max(price / 2 for panel in data["bowl"] for _, _, price in panel["prices"])
    kappa = bowl["shared field"]["covariance"][0][1]
    fig = plt.figure(figsize=(15, 9.6))
    top = [("anti-correlated end", "the anti-correlated end, κ = −γ/(N−1)"),
           ("featureless", "featureless, κ = 0: the Mirror's one centre"),
           ("shared field", f"a shared field, κ = {kappa}")]
    for k, (name, title) in enumerate(top):
        ax = fig.add_axes([0.04 + k * 0.31, 0.475, 0.27, 0.39])
        bowl_panel(ax, bowl[name], vmax, orbit=True)
        ax.set_title(title, fontsize=9.5, pad=6)
    bottom = [("integer pattern", "an integer pattern"), ("ring", "a ring of covariances"),
              ("star", "a star of covariances"), ("shared field", "a shared field")]
    for k, (name, title) in enumerate(bottom):
        ax = fig.add_axes([0.035 + k * 0.24, 0.05, 0.20, 0.33])
        bowl_panel(ax, bowl[name], vmax)
        qual = bowl[name]["qualifying"]
        ax.set_title(f"{title}\none K on popcounts {{{','.join(map(str, qual))}}}: {len(qual) ** 2} blocks",
                     fontsize=8.5, pad=5, color=AMBER)
    fig.text(0.5, 0.965, f"The atmosphere's price over the block lattice, N = {n}, γ = {gamma}", ha="center",
             fontsize=13, color=GREEN)
    fig.text(0.5, 0.938, "each tile: (K_p + K_q)/2, the depth of the centre about which the spectra of block (p, q) "
             "and of its mirror image (q, N − p) pair", ha="center", fontsize=8.5)
    fig.text(0.5, 0.918, "(for a popcount-conserving H that commutes with flipping every spin); "
             "pink: one quarter-turn orbit; circles: equal depth", ha="center", fontsize=8.5)
    fig.text(0.5, 0.415, "four weathers, four shapes: the blocks whose two popcounts each see one value of K, each pairing "
             "about a centre of its own (dashed: K takes several values on at least one of the two)", ha="center", fontsize=8.5)
    fig.savefig(out, dpi=170)
    plt.close(fig)


def checked_reading(panel):
    """The reading in a panel's title, checked against the exported classes."""
    n, name = panel["N"], panel["name"]
    big = [sorted(cl) for cl in panel["classes"] if len(cl) > 1]
    full = (1 << n) - 1
    if name == "anti-correlated end":
        assert big == [[0, full]], big
        return "the GHZ pair"
    if name == "featureless":
        assert big == [], big
        return "the populations alone"
    if name == "full correlation":
        assert sorted(sorted(cl) for cl in panel["classes"]) == sorted(
            sorted(i for i in range(full + 1) if bin(i).count("1") == p) for p in range(n + 1)), "one class per popcount"
        assert panel["dark"] == comb(2 * n, n)
        return "the popcount sectors"
    if name == "three anti-correlated pairs":
        alike = sorted(i for i in range(full + 1) if all(((i >> k) & 1) == ((i >> (k + 1)) & 1) for k in range(0, n, 2)))
        assert sorted(max(panel["classes"], key=len)) == alike and panel["dark"] == 6 ** (n // 2)
        return "each pair cannot tell 00 from 11"
    if name == "ring, positive covariance":
        neel = sum(1 << k for k in range(0, n, 2))
        assert big == [sorted([neel, neel ^ full])], big
        return "the Néel pair"
    if name == "two triangles, positive covariance":
        assert big == [], big
        return "frustrated"
    raise KeyError(name)


def dark_classes(data, out):
    panels = data["cells"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 10.4))
    for ax, panel in zip(axes.flat, panels):
        n = panel["N"]
        d = 1 << n
        reading = checked_reading(panel)
        order = sorted(range(d), key=lambda i: (bin(i).count("1"), i))
        rates = np.array(panel["rates"])[np.ix_(order, order)]
        mag = -rates.astype(float)
        assert int((mag == 0).sum()) == panel["dark"]
        life = np.where(mag == 0, 1.0, 1.0 - mag / mag.max())
        img = LIFE(np.where(mag == 0, 0.0, 0.12 + life * 0.73))       # the fastest cell stays visible
        img[mag == 0] = to_rgba(GREEN)
        ax.imshow(img, interpolation="nearest", origin="upper", extent=(-0.5, d - 0.5, d - 0.5, -0.5), zorder=1)
        ys, xs = np.nonzero(mag == 0)
        for i in range(4, 0, -1):
            ax.scatter(xs, ys, s=2.5 * (1 + 5.0 * i / 4), color=GREEN, alpha=0.04, edgecolors="none", zorder=3)
        off = [(x, y) for x, y in zip(xs, ys) if x != y]
        if 0 < len(off) <= 4:
            ox, oy = zip(*off)
            ax.scatter(ox, oy, s=110, facecolors="none", edgecolors=AMBER, linewidths=1.1, zorder=6, clip_on=False)
        sizes = np.cumsum([0] + [comb(n, p) for p in range(n + 1)])
        for e in sizes[1:-1]:
            ax.axhline(e - 0.5, color=CYAN, lw=0.35, alpha=0.35, zorder=2)
            ax.axvline(e - 0.5, color=CYAN, lw=0.35, alpha=0.35, zorder=2)
        ticks = [(sizes[p] + sizes[p + 1] - 1) / 2 for p in range(n + 1)]
        ax.set_xlim(-0.5, d - 0.5); ax.set_ylim(d - 0.5, -0.5)
        ax.set_xticks(ticks); ax.set_xticklabels(range(n + 1), fontsize=6.5)
        ax.set_yticks(ticks); ax.set_yticklabels(range(n + 1), fontsize=6.5)
        ax.tick_params(length=0, pad=8)
        ax.set_xlabel("bra, by popcount", fontsize=7.5); ax.set_ylabel("ket, by popcount", fontsize=7.5)
        big = sorted((len(cl) for cl in panel["classes"] if len(cl) > 1), reverse=True)
        parts = [f"{s}" if big.count(s) == 1 else f"{s}×{big.count(s)}" for s in sorted(set(big), reverse=True)]
        class_text = (("classes " + ", ".join(parts) + ", the rest single") if parts
                      else "no class of more than one configuration")
        ax.set_title(f"{panel['name']}: {reading}\n{panel['dark']} dark cells; {class_text}", fontsize=8.5, pad=5)
        for s in ax.spines.values():
            s.set_color(GRID)
    fig.text(0.5, 0.975, f"What the atmosphere cannot tell apart: the {4 ** panels[0]['N']} cells |ket⟩⟨bra| at N = {panels[0]['N']}",
             ha="center", fontsize=13, color=GREEN)
    fig.text(0.5, 0.955, "green: a dark cell (rate exactly 0), ringed where few lie off the diagonal; "
             "lines: the popcount blocks", ha="center", fontsize=8.5)
    fig.text(0.5, 0.937, "the rest shaded by the noise's rate alone, brighter = slower, the colour position "
             "linear in the rate (per panel)", ha="center", fontsize=8.5)
    fig.subplots_adjust(left=0.04, right=0.99, bottom=0.05, top=0.885, wspace=0.18, hspace=0.32)
    fig.savefig(out, dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    with open(DATA, encoding="utf-8") as f:
        data = json.load(f)
    landscape(data, os.path.join(OUT, "atmosphere_price_landscape.png"))
    dark_classes(data, os.path.join(OUT, "atmosphere_dark_classes.png"))
    print("wrote visualizations/atmosphere_price_landscape.png and visualizations/atmosphere_dark_classes.png")
