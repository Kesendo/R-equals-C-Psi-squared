#!/usr/bin/env python3
"""q_scale_n_scaling_summary.py

Cross-N summary plots: W_plateau(c) by N, abs(K)_peak(c) by N, Q_peak by c.
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
    parser.add_argument("--tags", type=str, default="N4_blockL,N5_blockL,N6_blockL,N7_gamma0_0p05",
                        help="Comma-separated tags to include")
    parser.add_argument("--out", type=str, default="summary_N4_N7")
    args = parser.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    N_DATA = {}
    for tag in args.tags.split(","):
        tag = tag.strip()
        path = RESULTS_DIR / f"q_scan_{tag}.json"
        if not path.exists():
            print(f"Skipping missing {path}")
            continue
        data = json.load(open(path))
        N = data["config"]["N"]
        N_DATA[N] = data

    Q_CAP = 3.0  # primary peak window

    # Compute aggregates per (N, c)
    agg = {}  # (N, c) -> list of blocks with that N, c
    for N, data in N_DATA.items():
        for bk, bd in data["blocks"].items():
            c = bd["chromaticity"]
            if c == 1:
                continue
            qdata = bd["q_data"]
            Qs = sorted(float(q) for q in qdata.keys())
            Ws = [qdata[str(q)]["W"] for q in Qs]
            Ks = [qdata[str(q)]["abs_K_max"] for q in Qs]
            # Primary peak in Q <= Q_CAP
            primary_Qs = [q for q in Qs if q <= Q_CAP]
            primary_Ws = Ws[:len(primary_Qs)]
            primary_Ks = Ks[:len(primary_Qs)]
            idx_W = int(np.argmax(primary_Ws))
            idx_K = int(np.argmax(primary_Ks))
            Q_onset = next((q for q, w in zip(Qs, Ws) if w > 0.05), None)
            entry = {
                "N": N,
                "n": int(bk.split("_")[1]),
                "c": c,
                "Q_onset": Q_onset,
                "Q_peak_W": primary_Qs[idx_W],
                "W_peak": primary_Ws[idx_W],
                "Q_peak_K": primary_Qs[idx_K],
                "K_peak": primary_Ks[idx_K],
                "W_plateau": Ws[-1],
                "Qs": Qs,
                "Ws": Ws,
                "Ks": Ks,
            }
            agg.setdefault((N, c), []).append(entry)

    # --- Figure 1: W_plateau(c) for each N ---
    fig, ax = plt.subplots(figsize=(7, 5))
    Ns_sorted = sorted(N_DATA.keys())
    for N in Ns_sorted:
        cs = sorted({key[1] for key in agg.keys() if key[0] == N})
        W_plat_vals = [np.mean([e["W_plateau"] for e in agg[(N, c)]]) for c in cs]
        ax.plot(cs, W_plat_vals, "o-", label=f"N={N}", markersize=8)
    ax.set_xlabel("chromaticity c")
    ax.set_ylabel("W_plateau (Q=50)")
    ax.set_title("W_plateau vs chromaticity, by N")
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_xticks(range(1, 6))
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / f"{args.out}_W_plateau_vs_c.png", dpi=120)
    plt.close(fig)

    # --- Figure 2: abs(K)_peak vs c ---
    fig, ax = plt.subplots(figsize=(7, 5))
    for N in Ns_sorted:
        cs = sorted({key[1] for key in agg.keys() if key[0] == N})
        K_peak_vals = [np.mean([e["K_peak"] for e in agg[(N, c)]]) for c in cs]
        ax.plot(cs, K_peak_vals, "o-", label=f"N={N}", markersize=8)
    ax.set_xlabel("chromaticity c")
    ax.set_ylabel("abs(K_CC_pr)_peak (primary, Q<=3)")
    ax.set_title("abs(K)_peak vs chromaticity, by N")
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_xticks(range(1, 6))
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / f"{args.out}_K_peak_vs_c.png", dpi=120)
    plt.close(fig)

    # --- Figure 3: Q_peak(W) and Q_peak(K) vs c, per N ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    for N in Ns_sorted:
        cs = sorted({key[1] for key in agg.keys() if key[0] == N})
        Q_W_vals = [np.mean([e["Q_peak_W"] for e in agg[(N, c)]]) for c in cs]
        Q_K_vals = [np.mean([e["Q_peak_K"] for e in agg[(N, c)]]) for c in cs]
        ax1.plot(cs, Q_W_vals, "o-", label=f"N={N}", markersize=8)
        ax2.plot(cs, Q_K_vals, "o-", label=f"N={N}", markersize=8)
    ax1.axhline(1.5, color="black", linewidth=0.5, alpha=0.5, linestyle=":")
    ax1.set_xlabel("c")
    ax1.set_ylabel("Q_peak(W) primary (Q<=3)")
    ax1.set_title("Q_peak(W) vs c")
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_xticks(range(1, 6))
    ax2.axhline(1.5, color="black", linewidth=0.5, alpha=0.5, linestyle=":")
    ax2.set_xlabel("c")
    ax2.set_ylabel("Q_peak(|K|) primary (Q<=3)")
    ax2.set_title("Q_peak(|K|) vs c")
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.set_xticks(range(1, 6))
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / f"{args.out}_Q_peak_vs_c.png", dpi=120)
    plt.close(fig)

    # Summary text
    print(f"{'N':>3}  {'c':>3}  {'#blocks':>8}  {'Q_onset_mean':>14}  {'Q_peak_W_mean':>14}  {'W_peak_mean':>12}  {'Q_peak_K_mean':>14}  {'K_peak_mean':>12}  {'W_plat_mean':>12}")
    for N in Ns_sorted:
        cs_here = sorted({key[1] for key in agg.keys() if key[0] == N})
        for c in cs_here:
            entries = agg[(N, c)]
            n_blocks = len(entries)
            qo_mean = np.mean([e["Q_onset"] for e in entries if e["Q_onset"] is not None])
            qpw_mean = np.mean([e["Q_peak_W"] for e in entries])
            wp_mean = np.mean([e["W_peak"] for e in entries])
            qpk_mean = np.mean([e["Q_peak_K"] for e in entries])
            kp_mean = np.mean([e["K_peak"] for e in entries])
            wpl_mean = np.mean([e["W_plateau"] for e in entries])
            print(f"{N:>3}  {c:>3}  {n_blocks:>8}  {qo_mean:>14.3f}  {qpw_mean:>14.3f}  {wp_mean:>12.4f}  {qpk_mean:>14.3f}  {kp_mean:>12.4e}  {wpl_mean:>12.4f}")

    print(f"\nSaved: {args.out}_W_plateau_vs_c.png, _K_peak_vs_c.png, _Q_peak_vs_c.png")


if __name__ == "__main__":
    main()
