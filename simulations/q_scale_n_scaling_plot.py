#!/usr/bin/env python3
"""q_scale_n_scaling_plot.py

Plot W(Q) and abs(K)_max(Q) curves across blocks at a given N. Also a
summary plot: Q_onset, Q_peak, W_peak, W_plateau vs (N, n).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

RESULTS_DIR = Path(__file__).parent / "results" / "q_scale_n_scaling"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tag", type=str, required=True)
    args = parser.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data = json.load(open(RESULTS_DIR / f"q_scan_{args.tag}.json"))
    cfg = data["config"]
    N = cfg["N"]
    gamma_0 = cfg["gamma_0"]

    # Figure 1: W(Q) and abs(K)_max(Q) per block
    blocks_sorted = sorted(data["blocks"].keys(), key=lambda k: int(k.split("_")[1]))
    n_blocks = len(blocks_sorted)
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, n_blocks))

    fig, (axW, axK) = plt.subplots(1, 2, figsize=(12, 5))
    for i, bk in enumerate(blocks_sorted):
        n_val = int(bk.split("_")[1])
        c = data["blocks"][bk]["chromaticity"]
        q_data = data["blocks"][bk]["q_data"]
        Qs = sorted(float(q) for q in q_data.keys())
        Ws = [q_data[str(q)]["W"] for q in Qs]
        Ks = [q_data[str(q)]["abs_K_max"] for q in Qs]

        label = f"n={n_val} (c={c})"
        axW.plot(Qs, Ws, color=colors[i], label=label, marker="o", markersize=4)
        axK.plot(Qs, Ks, color=colors[i], label=label, marker="o", markersize=4)

    axW.set_xscale("log")
    axW.set_xlabel("Q = J/γ₀")
    axW.set_ylabel("W (dressed mode weight)")
    axW.set_title(f"W(Q) at N={N}, γ₀={gamma_0}")
    axW.axvline(1.5, color="black", linewidth=0.5, alpha=0.5, linestyle=":")
    axW.axvline(0.3, color="gray", linewidth=0.5, alpha=0.5, linestyle=":")
    axW.grid(True, which="both", alpha=0.3)
    axW.legend(fontsize=8)
    axW.set_ylim(-0.05, 1.05)

    axK.set_xscale("log")
    axK.set_yscale("log")
    axK.set_xlabel("Q = J/γ₀")
    axK.set_ylabel("abs(K_CC_pr)_max over t")
    axK.set_title(f"abs(K)_max(Q) at N={N}, γ₀={gamma_0}")
    axK.axvline(1.5, color="black", linewidth=0.5, alpha=0.5, linestyle=":")
    axK.grid(True, which="both", alpha=0.3)
    axK.legend(fontsize=8)

    fig.suptitle(f"N={N} Q-scan across (n, n+1) blocks", fontsize=11)
    fig.tight_layout()
    out_path = RESULTS_DIR / f"q_scan_{args.tag}_curves.png"
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
