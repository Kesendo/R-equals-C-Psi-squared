#!/usr/bin/env python3
"""eq018_kcc_pr_plot.py

Plot S_A(t) vs S_ref(t) = S(0) exp(-4 gamma_0 t) and K_CC_pr(t) from
kcc_pr_{tag}.json. Saves one figure per bond showing all five (n, n+1)
curves.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


RESULTS_DIR = Path(__file__).parent / "results" / "eq018_kcc_pr_extension"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tag", type=str, default="n5_full")
    args = parser.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data = json.load(open(RESULTS_DIR / f"kcc_pr_{args.tag}.json"))
    cfg = data["config"]
    N = cfg["N"]
    gamma_0 = cfg["gamma_0"]

    # Colors per (n, n+1)
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, N))

    # One figure with 2 rows per bond
    bonds = sorted(data["bonds"].keys(), key=lambda s: int(s.replace("bond_", "")))
    n_bonds = len(bonds)
    fig, axes = plt.subplots(2, n_bonds, figsize=(4.5 * n_bonds, 8), sharex=True)
    if n_bonds == 1:
        axes = axes[:, None]

    for b_idx, bond_key in enumerate(bonds):
        bond_num = int(bond_key.replace("bond_", ""))
        ax_S = axes[0, b_idx]
        ax_K = axes[1, b_idx]
        for n_coh in range(N):
            n_key = f"n_{n_coh}"
            if n_key not in data["bonds"][bond_key]:
                continue
            d = data["bonds"][bond_key][n_key]
            times = np.array(d["times"])
            S_A = np.array(d["S_A"])
            S_ref = np.array(d["S_ref"])
            K = np.array(d["K"])
            label = f"(n,n+1)=({n_coh},{n_coh+1})"
            ax_S.semilogy(times, np.maximum(S_A, 1e-20), color=colors[n_coh],
                          linewidth=1.5, label=label)
            ax_S.semilogy(times, np.maximum(S_ref, 1e-20), color=colors[n_coh],
                          linewidth=1.0, linestyle=":")
            ax_K.plot(times, K, color=colors[n_coh], linewidth=1.5, label=label)

        ax_S.set_title(f"S_A(t) (solid) vs S_ref=S(0)exp(-4gamma t) (dotted), bond {bond_num}")
        ax_S.set_ylabel("S(t)")
        ax_S.grid(True, which="both", alpha=0.3)
        ax_S.legend(fontsize=7, loc="upper right")

        ax_K.set_title(f"K_CC[n,n+1]_pr(t), bond {bond_num}")
        ax_K.set_xlabel("t")
        ax_K.set_ylabel("K_CC_pr(t)")
        ax_K.axhline(0.0, color="black", linewidth=0.6, alpha=0.5)
        ax_K.grid(True, which="both", alpha=0.3)
        ax_K.legend(fontsize=7, loc="upper right")

    fig.suptitle(
        f"K_CC[n, n+1]_pr extension at N={N}, gamma_0={gamma_0}, "
        f"J={cfg['J']}, dJ={cfg['dJ']}",
        fontsize=11,
    )
    fig.tight_layout()
    out_path = RESULTS_DIR / f"kcc_pr_{args.tag}_trajectories.png"
    fig.savefig(out_path, dpi=120)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
