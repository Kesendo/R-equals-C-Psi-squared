"""Compare bond-isolate S(t) trajectories across all bonds at N=7.

Reads the 6 CSVs produced by `dotnet run ... bond-isolate --N 7 --bond {0..5}`
with default ρ_cc probe at J=0.075, γ=0.05, tmax=30, dt=0.1.

Question: do F71-mirror-pair bonds (0↔5, 1↔4, 2↔3) produce identical S(t)?
And how does the orbit (Endpoint vs mid vs Center-near) shape decay?

Outputs:
  - Markdown table to stdout: bond, S(0), τ_half, decay-rate, F71-mirror-residual
  - 2-panel plot saved to simulations/results/bond_isolate/comparison_N7_S_overlay.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO_ROOT = Path(__file__).resolve().parents[1]
CSV_DIR = REPO_ROOT / "simulations" / "results" / "bond_isolate"
PARAMS = "J0.0750_gamma0.0500_probe-coherence"
N = 7
BONDS = list(range(N - 1))  # 0..5 for N=7

# F71 spatial mirror: bond b ↔ bond (N-2 - b)
ORBIT = {0: 0, 5: 0, 1: 1, 4: 1, 2: 2, 3: 2}  # orbit label per bond
ORBIT_NAME = {0: "Endpoint", 1: "mid-flank", 2: "Center-near"}
COLORS = {0: "tab:red", 1: "tab:blue", 2: "tab:green"}


def load_bond(b: int) -> tuple[np.ndarray, np.ndarray]:
    path = CSV_DIR / f"N{N}_b{b}_{PARAMS}.csv"
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    t = data[:, 0]
    S = data[:, -1]
    return t, S


def half_life(t: np.ndarray, S: np.ndarray) -> float:
    """First time t* with S(t*) ≤ S(0)/2. Linearly interpolated."""
    target = S[0] / 2.0
    below = np.where(S <= target)[0]
    if below.size == 0:
        return float("nan")
    i = int(below[0])
    if i == 0:
        return 0.0
    # linear interp between i-1 and i
    s_lo, s_hi = S[i - 1], S[i]
    t_lo, t_hi = t[i - 1], t[i]
    frac = (s_lo - target) / (s_lo - s_hi) if s_lo != s_hi else 0.0
    return float(t_lo + frac * (t_hi - t_lo))


def decay_rate(t: np.ndarray, S: np.ndarray, t_min: float = 5.0, t_max: float = 25.0) -> float:
    """Slope of -log(S) vs t over a clean late-time window. Units: 1/time."""
    mask = (t >= t_min) & (t <= t_max) & (S > 0)
    if mask.sum() < 3:
        return float("nan")
    slope, _ = np.polyfit(t[mask], -np.log(S[mask]), 1)
    return float(slope)


def main() -> None:
    series = {b: load_bond(b) for b in BONDS}

    # F71 mirror residual: max |S_b(t) - S_{N-2-b}(t)| / S(0)
    pairs = [(0, 5), (1, 4), (2, 3)]
    mirror_resid = {}
    for a, b in pairs:
        t_a, S_a = series[a]
        t_b, S_b = series[b]
        diff = np.max(np.abs(S_a - S_b))
        mirror_resid[(a, b)] = diff / S_a[0]

    # PAIRWISE check: is S(t) bit-identical across ALL bonds (universal),
    # or only within F71-mirror pairs (geometric)?
    print()
    print("# Pairwise max |S_b(t) - S_{b'}(t)| / S(0)  —  bond-pair universality check")
    print()
    header = "| b\\b' | " + " | ".join(f"{b}" for b in BONDS) + " |"
    print(header)
    print("|" + "------|" * (len(BONDS) + 1))
    for a in BONDS:
        _, Sa = series[a]
        row = [f"  {a}  "]
        for b in BONDS:
            _, Sb = series[b]
            d = np.max(np.abs(Sa - Sb)) / Sa[0]
            row.append(f"{d:.2e}")
        print("| " + " | ".join(row) + " |")
    # Also: max relative diff at any time across all (a,b) pairs
    max_rel_diff = 0.0
    for a in BONDS:
        _, Sa = series[a]
        for b in BONDS:
            if b <= a:
                continue
            _, Sb = series[b]
            mask = Sa > 1e-15
            rel = np.max(np.abs(Sa[mask] - Sb[mask]) / np.abs(Sa[mask]))
            if rel > max_rel_diff:
                max_rel_diff = rel
    print()
    print(f"Maximum POINTWISE relative diff across all pairs: {max_rel_diff:.2e}")

    # Table
    print()
    print("# Bond-isolate N=7 comparison (ρ_cc probe, J=0.075, γ=0.05)")
    print()
    print("| Bond | Orbit       | F71-mirror | S(0)     | τ_half   | Decay rate Γ | Mirror residual |")
    print("|------|-------------|------------|----------|----------|--------------|-----------------|")
    for b in BONDS:
        t, S = series[b]
        mirror = N - 2 - b
        orbit_name = ORBIT_NAME[ORBIT[b]]
        S0 = S[0]
        tau = half_life(t, S)
        gamma_eff = decay_rate(t, S)
        # which pair this bond belongs to
        pair = tuple(sorted((b, mirror)))
        resid = mirror_resid.get(pair, float("nan"))
        print(
            f"| {b}    | {orbit_name:11s} | b={mirror}        | {S0:.6f} | {tau:7.3f}  | {gamma_eff:.6f}     | {resid:.2e}        |"
        )
    print()
    print(f"S(0) closed-form: (N-1)/N = 6/7 = {6/7:.6f}")
    print()

    # Plot: linear S(t) and log S(t)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    for b in BONDS:
        t, S = series[b]
        orbit_id = ORBIT[b]
        color = COLORS[orbit_id]
        # solid for left half (b in 0..2), dashed for right half (b in 3..5)
        linestyle = "-" if b <= 2 else "--"
        label = f"b={b} ({ORBIT_NAME[orbit_id]})"
        ax1.plot(t, S, color=color, linestyle=linestyle, label=label, alpha=0.85)
        ax2.semilogy(t, S, color=color, linestyle=linestyle, label=label, alpha=0.85)
    ax1.set_ylabel("S(t) (spatial-sum coherence)")
    ax1.set_title(f"Bond-isolate N={N} comparison: S(t) per active bond (ρ_cc probe)")
    ax1.axhline(6 / 7, color="black", linestyle=":", alpha=0.4, label="S(0) = 6/7")
    ax1.legend(loc="upper right", fontsize=9)
    ax1.grid(alpha=0.3)
    ax2.set_xlabel("t")
    ax2.set_ylabel("log S(t)")
    ax2.set_title("Log-scale: linear region exposes decay rate")
    ax2.grid(alpha=0.3, which="both")
    plt.tight_layout()
    out = CSV_DIR / "comparison_N7_S_overlay.png"
    plt.savefig(out, dpi=120)
    print(f"Plot saved: {out.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
