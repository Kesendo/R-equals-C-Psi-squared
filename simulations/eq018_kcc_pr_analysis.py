#!/usr/bin/env python3
"""eq018_kcc_pr_analysis.py

Analysis layer on top of kcc_pr_{tag}.json:
  - Effective decay rate gamma_eff(t) = -d(log S_A)/dt, compared to 4 gamma_0.
  - Deviation rate D(t) = gamma_eff(t) - 4 gamma_0 = "extra decay" from
    H-induced leak out of the 1-site-differing slice.
  - Ratio R(t) = S_A(t) / S_ref(t), bounded above by 1 when amplitude
    leaks out.

Also produces a compact table of D(t=20) per (bond, n) for the RESULT.
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

    data = json.load(open(RESULTS_DIR / f"kcc_pr_{args.tag}.json"))
    cfg = data["config"]
    N = cfg["N"]
    gamma_0 = cfg["gamma_0"]
    t_pw = cfg["T_pointwise"]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    bonds = sorted(data["bonds"].keys(), key=lambda s: int(s.replace("bond_", "")))
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, N))

    fig, axes = plt.subplots(2, len(bonds), figsize=(4.2 * len(bonds), 7), sharex=True)
    if len(bonds) == 1:
        axes = axes[:, None]

    # Deviation table
    print("=" * 78)
    print(f"Deviation table at N = {N}, t = {t_pw}:")
    print(f"  R(t) = S_A(t) / S_ref(t)  where S_ref = S(0) exp(-4 gamma_0 t)")
    print(f"  D(t) = gamma_eff(t) - 4 gamma_0  [extra decay rate from H-leak]")
    print(f"    gamma_eff(t) = -d(log S_A)/dt, estimated from adjacent time points")
    print("=" * 78)

    for b_idx, bond_key in enumerate(bonds):
        bond_num = int(bond_key.replace("bond_", ""))
        ax_R = axes[0, b_idx]
        ax_D = axes[1, b_idx]

        print(f"\nBond {bond_num}:")
        print(f"  {'(n, n+1)':>10}  {'S_A(t=' + str(t_pw) + ')':>16}  "
              f"{'S_ref(t=' + str(t_pw) + ')':>16}  "
              f"{'R':>10}  {'D [rate]':>12}")

        for n_coh in range(N):
            n_key = f"n_{n_coh}"
            if n_key not in data["bonds"][bond_key]:
                continue
            d = data["bonds"][bond_key][n_key]
            times = np.array(d["times"])
            S_A = np.array(d["S_A"])
            S_ref = np.array(d["S_ref"])
            # Avoid divide-by-zero
            with np.errstate(divide="ignore", invalid="ignore"):
                R = np.where(S_ref > 1e-20, S_A / S_ref, 0.0)
                # gamma_eff(t) = -d(log S_A)/dt via central difference
                log_S = np.log(np.maximum(S_A, 1e-30))
                gamma_eff = np.zeros_like(times)
                for k in range(1, len(times) - 1):
                    gamma_eff[k] = -(log_S[k + 1] - log_S[k - 1]) / (times[k + 1] - times[k - 1])
                D = gamma_eff - 4.0 * gamma_0

            # Extract at t = t_pw
            idx_pw = int(np.argmin(np.abs(times - t_pw)))
            R_pw = R[idx_pw]
            D_pw = D[idx_pw]
            print(f"  ({n_coh}, {n_coh+1})     "
                  f"{S_A[idx_pw]:>16.4e}  {S_ref[idx_pw]:>16.4e}  "
                  f"{R_pw:>10.5f}  {D_pw:>+12.4e}")

            # Plot R and D over time
            label = f"({n_coh},{n_coh+1})"
            ax_R.plot(times, R, color=colors[n_coh], linewidth=1.3, label=label)
            # D is noisy at endpoints; plot from k=1 to k=N-2
            ax_D.plot(times[1:-1], D[1:-1], color=colors[n_coh],
                      linewidth=1.3, label=label)

        ax_R.set_title(f"R(t)=S_A/S_ref, bond {bond_num}")
        ax_R.set_ylabel("R(t)")
        ax_R.axhline(1.0, color="black", linewidth=0.6, alpha=0.5)
        ax_R.grid(True, which="both", alpha=0.3)
        ax_R.legend(fontsize=7, loc="upper right")

        ax_D.set_title(f"D(t) = gamma_eff - 4 gamma_0, bond {bond_num}")
        ax_D.set_xlabel("t")
        ax_D.set_ylabel("D [extra decay rate]")
        ax_D.axhline(0.0, color="black", linewidth=0.6, alpha=0.5)
        ax_D.grid(True, which="both", alpha=0.3)
        ax_D.legend(fontsize=7, loc="upper right")

    fig.suptitle(
        f"Effective-rate analysis at N={N}, gamma_0={gamma_0}",
        fontsize=11,
    )
    fig.tight_layout()
    out_path = RESULTS_DIR / f"kcc_pr_{args.tag}_rate_analysis.png"
    fig.savefig(out_path, dpi=120)
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
