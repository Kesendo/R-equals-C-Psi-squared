"""F127 residue-collapse figures for PROOF_F127_RESIDUE_COLLAPSE.md.

Renders two PNGs into visualizations/:
  f127_zero_valley.png   -- |T| on the constraint surface; the sheet curve is a
                            machine-zero valley (the SS3 core identity, seen).
  f127_sheet_lattice.png -- the 32 sheets x 9 events incidence, one event per
                            (i,j) block, colored by atom class (the SS2 lattice).

The zero-valley figure self-gates: it re-solves points ON the drawn sheet curve
and asserts |T| < 1e-12 there (the drawn curve must BE the theorem, not decoration).
Run from the repo root:  python simulations/f127_viz.py
"""
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

sys.path.insert(0, "simulations")
import f127_sheet_lattice as SL

OUT_DIR = "visualizations"
BG = "#0d1117"
NEON = ["#00e5ff", "#ff2ec4", "#adff2f", "#ffb347"]


def alpha(ang, i):
    j, l = [t for t in range(3) if t != i]
    return math.sin(ang[l]) - math.sin(ang[j])


def T(a, b):
    tot = 0.0
    for i in range(3):
        for j in range(3):
            d = math.tan((a[i] + b[j]) / 2)
            if abs(d) < 1e-9:
                return float("nan")
            tot += ((-1) ** (i + j)) * alpha(a, i) * alpha(b, j) / d
    return tot


def zero_valley():
    a1, a2 = 0.9, 2.1
    a3 = math.acos(-math.cos(a1) - math.cos(a2))
    A = [a1, a2, a3]
    Sa = sum(A)

    n = 601
    b1s = np.linspace(-math.pi, math.pi, n)
    b2s = np.linspace(-math.pi, math.pi, n)
    Z = np.full((n, n), np.nan)
    SHEET = np.full((n, n), np.nan)      # sin(v/2): 0 iff v = 0 mod 2pi, no wrap
    for iy, b2 in enumerate(b2s):
        for ix, b1 in enumerate(b1s):
            c3 = -math.cos(b1) - math.cos(b2)
            if abs(c3) > 1.0:
                continue
            b3 = math.acos(c3)           # the + branch of Cb
            Z[iy, ix] = abs(T(A, [b1, b2, b3]))
            SHEET[iy, ix] = math.sin((Sa + b1 + b2 + b3) / 2)

    # the gate: points solved ON the drawn curve must be machine zeros of T
    checked, worst = 0, 0.0
    for b1 in np.linspace(-3.0, 3.0, 61):
        def g(b2, b1=b1):
            c3 = -math.cos(b1) - math.cos(b2)
            if abs(c3) > 1.0:
                return None
            return math.sin((Sa + b1 + b2 + math.acos(c3)) / 2)
        prev = None
        for b2 in np.linspace(-3.1, 3.1, 400):
            cur = g(b2)
            if prev is not None and cur is not None and (prev < 0) != (cur < 0):
                lo, hi = b2 - 6.2 / 399, b2
                for _ in range(60):
                    mid = 0.5 * (lo + hi)
                    if (g(lo) < 0) != (g(mid) < 0):
                        hi = mid
                    else:
                        lo = mid
                b2c = 0.5 * (lo + hi)
                b3c = math.acos(-math.cos(b1) - math.cos(b2c))
                t = abs(T(A, [b1, b2c, b3c]))
                if not math.isnan(t):
                    checked += 1
                    worst = max(worst, t)
            prev = cur
    assert checked > 10 and worst < 1e-12, f"sheet curve NOT a zero valley: {worst:.2e}"
    print(f"[gate] {checked} points on the drawn sheet curve, worst |T| = {worst:.2e}")

    fig, ax = plt.subplots(figsize=(8.6, 7), facecolor=BG)
    ax.set_facecolor(BG)
    im = ax.pcolormesh(b1s, b2s, np.log10(np.clip(Z, 1e-18, None)),
                       cmap="magma", vmin=-16, vmax=2, shading="auto")
    ax.contour(b1s, b2s, SHEET, levels=[0.0], colors=[NEON[0]], linewidths=1.4)
    cb = fig.colorbar(im, ax=ax, fraction=0.046)
    cb.set_label("log10 |T|", color="w")
    plt.setp(cb.ax.get_yticklabels(), color="w")
    cb.ax.yaxis.set_tick_params(color="w")
    ax.set_xlabel("b1", color="w")
    ax.set_ylabel("b2", color="w")
    ax.tick_params(colors="w")
    ax.set_title("The core identity, seen: |T| on the constraint surface\n"
                 "(a on Ca fixed, b3 from Cb; cyan = the sheet, a machine-zero valley;\n"
                 "other dark lines are ordinary sign changes, allowed anywhere)",
                 color="w", fontsize=11)
    fig.tight_layout()
    out = os.path.join(OUT_DIR, "f127_zero_valley.png")
    fig.savefig(out, dpi=140, facecolor=BG)
    print("saved", out)


def sheet_lattice():
    sheets = {}
    for atom, s, tau, L in SL.EVENTS:
        key, eps = SL.canon(L)
        sheets.setdefault(key, []).append(atom)
    keys = sorted(sheets.keys())
    classes = {("psum", "psum"): 1, ("psum", "pdif"): 2,
               ("pdif", "psum"): 3, ("pdif", "pdif"): 4}
    grid = np.zeros((32, 9))
    for r, key in enumerate(keys):
        for (i, j, xi, up, e) in sheets[key]:
            grid[r, 3 * i + j] = classes[(xi, up)]
    assert (grid > 0).all(), "a (sheet, block) cell has no event"

    fig, ax = plt.subplots(figsize=(8.6, 7), facecolor=BG)
    ax.set_facecolor(BG)
    ax.pcolormesh(np.arange(10), np.arange(33), grid,
                  cmap=ListedColormap([BG] + NEON), vmin=0, vmax=4,
                  edgecolors="#222933", linewidth=0.4, shading="flat")
    ax.set_xticks(np.arange(9) + 0.5)
    ax.set_xticklabels([f"({i},{j})" for i in range(3) for j in range(3)],
                       color="w", fontsize=8)
    ax.set_yticks(np.arange(32) + 0.5)
    ax.set_yticklabels(["".join("+" if c > 0 else "-" for c in k) for k in keys],
                       color="w", fontsize=6, family="monospace")
    ax.invert_yaxis()
    ax.tick_params(colors="w", length=0)
    ax.set_xlabel("(i, j) block", color="w")
    ax.set_ylabel("sheet (sign pattern over a1 a2 a3 b1 b2 b3)", color="w")
    ax.set_title("The sheet lattice: 32 sheets x 9 events, one per (i,j) block\n"
                 "color = atom class (xi, upsilon) in {psum, pdif}^2",
                 color="w", fontsize=11)
    handles = [plt.Rectangle((0, 0), 1, 1, color=NEON[k]) for k in range(4)]
    ax.legend(handles, ["(psum,psum)", "(psum,pdif)", "(pdif,psum)", "(pdif,pdif)"],
              loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=8,
              facecolor=BG, edgecolor="#444", labelcolor="w")
    fig.tight_layout()
    out = os.path.join(OUT_DIR, "f127_sheet_lattice.png")
    fig.savefig(out, dpi=140, facecolor=BG)
    print("saved", out)


if __name__ == "__main__":
    zero_valley()
    sheet_lattice()
