#!/usr/bin/env python3
"""Quantify the formation map: ridge line, transition width, k- and N-stability.

Closes the formation-possibility thread (f90fc89 → 5ab050b → 41f2ef5 → a8ba1d5) with numbers
instead of eyeballed contours. Four readings, all on the band-edge formation order parameter
of formation_possibility.py:

  1. Ridge fit. Δ_ridge(j2), the 0.5-crossing of the formation order, fitted as a line
     a + b·j2 per topology (chain, ring), over the full 21-row j2 grid of the map.
  2. Width fit. W(j2) = Δ(0.65-crossing) − Δ(0.35-crossing), the width of the MARGINAL
     window (the scanner's own regime boundaries). Tests "integrability breaking does not
     only shift the threshold, it softens it" as a quantitative statement.
  3. k-independence. The ridge at k=2,3,4 (chain): same line, and the sharpening with k
     (the width falling) made explicit.
  4. N-stability. The ridge at N=10,12,14 (k=3, chain) to bound finite-size drift.

Crossings use first-upward-crossing interpolation (robust to small non-monotonic wiggles
from band-edge degeneracies on the ring; the maximum wiggle is reported).

Writes simulations/results/formation_ridge_fit.tsv and formation_ridge_fit.png.

  python simulations/formation_ridge_fit.py
"""
import os
import sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "simulations")
from formation_possibility import formation_order  # noqa: E402

DELTAS = np.linspace(0.05, 5.0, 100)
LEVELS = (0.35, 0.5, 0.65)


def first_crossing(row, deltas, level):
    """Δ of the first upward crossing of `level`; linear interpolation between grid points.
    NaN if the row never reaches the level."""
    above = row >= level
    if not above.any() or above[0]:
        return float("nan")
    i = int(np.argmax(above))
    d0, d1, r0, r1 = deltas[i - 1], deltas[i], row[i - 1], row[i]
    return float(d0 + (level - r0) * (d1 - d0) / (r1 - r0))


def scan_row(N, k, j2, ring):
    row = np.array([formation_order(N, k, float(d), float(j2), ring) for d in DELTAS])
    wiggle = float(np.max(np.maximum(0.0, -np.diff(row))))
    d35, d50, d65 = (first_crossing(row, DELTAS, lv) for lv in LEVELS)
    return d35, d50, d65, wiggle


def linfit(x, y):
    """Least-squares line y = a + b·x; returns a, b, rms residual."""
    x, y = np.asarray(x), np.asarray(y)
    ok = np.isfinite(y)
    b, a = np.polyfit(x[ok], y[ok], 1)
    rms = float(np.sqrt(np.mean((a + b * x[ok] - y[ok]) ** 2)))
    return float(a), float(b), rms


def main():
    rows = []          # (section, topo, N, k, j2, d35, d50, d65, width, wiggle)
    j2s = np.linspace(0.0, 1.0, 21)

    print("Formation ridge quantification (band-edge order parameter, first-crossing readout)\n")

    # 1+2: ridge + width over the full map grid, chain and ring (N=12, k=3)
    main_fits = {}
    for ring in (False, True):
        topo = "ring" if ring else "chain"
        d50s, widths = [], []
        for j2 in j2s:
            d35, d50, d65, wig = scan_row(12, 3, j2, ring)
            w = d65 - d35
            d50s.append(d50)
            widths.append(w)
            rows.append(("map", topo, 12, 3, j2, d35, d50, d65, w, wig))
        a, b, rms = linfit(j2s, d50s)
        wa, wb, wrms = linfit(j2s, widths)
        main_fits[topo] = (a, b, rms, wa, wb, wrms, list(d50s), list(widths))
        print(f"  {topo:5s} N=12 k=3:  Δ_ridge = {a:.3f} + {b:.3f}·j2   (rms {rms:.3f})")
        print(f"  {topo:5s} N=12 k=3:  width   = {wa:.3f} + {wb:.3f}·j2   (rms {wrms:.3f})")

    # 3: k-independence of the ridge, sharpening of the width (chain)
    print("\n  k-dependence (chain, N=12):")
    print(f"    {'k':>3} {'j2':>5} {'Δ_ridge':>8} {'width':>7}")
    for k in (2, 3, 4):
        for j2 in (0.0, 0.5, 1.0):
            d35, d50, d65, wig = scan_row(12, k, j2, False)
            w = d65 - d35
            rows.append(("kdep", "chain", 12, k, j2, d35, d50, d65, w, wig))
            print(f"    {k:>3} {j2:>5.2f} {d50:>8.3f} {w:>7.3f}")

    # 4: N-stability of the ridge (chain, k=3)
    print("\n  N-stability (chain, k=3):")
    print(f"    {'N':>3} {'j2':>5} {'Δ_ridge':>8} {'width':>7}")
    for N in (10, 12, 14):
        for j2 in (0.0, 0.5, 1.0):
            d35, d50, d65, wig = scan_row(N, 3, j2, False)
            w = d65 - d35
            rows.append(("ndep", "chain", N, 3, j2, d35, d50, d65, w, wig))
            print(f"    {N:>3} {j2:>5.2f} {d50:>8.3f} {w:>7.3f}")

    max_wiggle = max(r[9] for r in rows)
    print(f"\n  max non-monotonic wiggle of any row: {max_wiggle:.2e}")

    os.makedirs("simulations/results", exist_ok=True)
    tsv = "simulations/results/formation_ridge_fit.tsv"
    with open(tsv, "w", encoding="utf-8", newline="\n") as f:
        f.write("section\ttopo\tN\tk\tj2\td35\td50\td65\twidth\twiggle\n")
        for sec, topo, N, k, j2, d35, d50, d65, w, wig in rows:
            f.write(f"{sec}\t{topo}\t{N}\t{k}\t{j2:.2f}\t{d35:.4f}\t{d50:.4f}\t{d65:.4f}\t{w:.4f}\t{wig:.2e}\n")
    print(f"  wrote {tsv}")

    fig, (axl, axr) = plt.subplots(1, 2, figsize=(11.0, 4.4))
    for topo, color in (("chain", "tab:blue"), ("ring", "tab:orange")):
        a, b, rms, wa, wb, wrms, d50s, widths = main_fits[topo]
        axl.plot(j2s, d50s, "o", ms=4, color=color, label=f"{topo}: {a:.2f} + {b:.2f}·j2")
        axl.plot(j2s, a + b * j2s, "-", lw=1.2, color=color, alpha=0.7)
        axr.plot(j2s, widths, "o", ms=4, color=color, label=f"{topo}: {wa:.2f} + {wb:.2f}·j2")
        axr.plot(j2s, wa + wb * j2s, "-", lw=1.2, color=color, alpha=0.7)
    axl.set_xlabel("integrability breaking  j2")
    axl.set_ylabel("Δ_ridge  (formation = 0.5)")
    axl.set_title("The marginal ridge is a line")
    axl.legend()
    axr.set_xlabel("integrability breaking  j2")
    axr.set_ylabel("width  Δ(0.65) − Δ(0.35)")
    axr.set_title("… and the transition broadens linearly")
    axr.legend()
    fig.suptitle("Formation map quantified: 3-body complex on N=12, chain vs ring", y=1.0)
    fig.tight_layout()
    out = "simulations/results/formation_ridge_fit.png"
    fig.savefig(out, dpi=130)
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()
