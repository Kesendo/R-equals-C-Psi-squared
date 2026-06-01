#!/usr/bin/env python3
"""The clean formation map: the band-edge formation order parameter over (Delta, j2), as a heatmap.

Reproduces the real artifact (the image) rather than scalar thresholds, so the marginal
near-threshold ridge is seen, not inferred from a few crossings. The ridge (formation ~ 0.5)
should appear as a sloped band, Delta_threshold ~ 1.15 + 1.3*j2.

  python simulations/formation_map.py [N] [k]
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "simulations")
from formation_possibility import formation_order  # noqa: E402


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    k = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    deltas = np.linspace(0.1, 3.0, 36)
    j2s = np.linspace(0.0, 1.0, 21)
    Z = np.zeros((len(j2s), len(deltas)))
    for i, j2 in enumerate(j2s):
        for jx, d in enumerate(deltas):
            Z[i, jx] = formation_order(N, k, float(d), float(j2))

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    im = ax.pcolormesh(deltas, j2s, Z, cmap="viridis", shading="auto", vmin=0.0, vmax=1.0)
    ax.contour(deltas, j2s, Z, levels=[0.5], colors="white", linewidths=2.0, linestyles="--")
    ax.set_xlabel("binding  Delta")
    ax.set_ylabel("integrability breaking  j2")
    ax.set_title(f"Formation map: {k}-body complex on an N={N} chain\n"
                 f"color = band-edge clustering (0 unformed, 1 deeply bound); white dashed = marginal ridge")
    fig.colorbar(im, ax=ax, label="formation order")
    os.makedirs("simulations/results", exist_ok=True)
    out = "simulations/results/formation_map.png"
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")

    print("\nmarginal ridge (Delta where formation crosses 0.5) vs j2:")
    for i in range(0, len(j2s), 4):
        row = Z[i]
        cross = float(np.interp(0.5, row, deltas)) if row.max() > 0.5 > row.min() else float("nan")
        print(f"  j2={j2s[i]:.2f}   Delta_ridge ~ {cross:.2f}")


if __name__ == "__main__":
    main()
