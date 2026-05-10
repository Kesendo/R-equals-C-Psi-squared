"""F89 -> F86 bridge test at N=11: numerical ∂S/∂J_b at one bond around uniform-J full chain.

N=11 c=2 has 5 distinct F71 orbits:
  Orbit 0 (Endpoint):     b=0, b=9
  Orbit 1 (next-flank):   b=1, b=8
  Orbit 2 (next-flank):   b=2, b=7
  Orbit 3 (mid-flank):    b=3, b=6
  Orbit 4 (Center):       b=4, b=5

F86 K_b at N=11 c=2 (Dicke probe, extended grid, from c2hwhm_N11_q32.txt):
  Endpoint  b=0: |K|max = 0.0766 at Q_peak = 2.50
  flank     b=1: |K|max = 0.1037 at Q_peak = 8.79
  flank     b=2: |K|max = 0.1076 at Q_peak = 13.61
  mid-flank b=3: |K|max = 0.0502 at Q_peak = 1.59
  Center    b=4: |K|max = 0.0622 at Q_peak = 21.94

bond-isolate samples at Q = J0/γ = 0.075/0.05 = 1.5 (single Q point); F86 |K|max
at each bond's own Q_peak (different Q per orbit). Comparing orbit-fan ranking.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "simulations" / "results" / "bond_isolate_perturbation"

DELTA_J = 0.01

ORBITS = [
    ("Endpoint b=0", "N11_full_pert_b0_dJ0.01.csv", 0.0766, 2.50),
    ("flank    b=1", "N11_full_pert_b1_dJ0.01.csv", 0.1037, 8.79),
    ("flank    b=2", "N11_full_pert_b2_dJ0.01.csv", 0.1076, 13.61),
    ("mid-flank b=3", "N11_full_pert_b3_dJ0.01.csv", 0.0502, 1.59),
    ("Center   b=4", "N11_full_pert_b4_dJ0.01.csv", 0.0622, 21.94),
]


def load_S(path: Path) -> tuple[np.ndarray, np.ndarray]:
    arr = np.loadtxt(path, delimiter=",", skiprows=1)
    return arr[:, 0], arr[:, -1]


def main() -> None:
    t, S_base = load_S(DATA / "N11_full_baseline.csv")

    print(f"\n# F89 → F86 bridge at N=11: numerical ∂S/∂J_b on ρ_cc probe (full chain, J0=0.075, γ=0.05, δJ={DELTA_J})\n")
    print(f"# Baseline at J0 uniform: S(0) = {S_base[0]:.6f}, S(t=10) = {S_base[100]:.6f}, S(t=20) = {S_base[200]:.6f}, S(t=30) = {S_base[300]:.6f}")
    print()
    print("| Orbit | K_b(t=2) | K_b(t=5) | K_b(t=10) | K_b(t=15) | K_b(t=20) | K_b(t=30) | peak |K_b(t)| | t at peak | F86 |K|max | F86 Q_peak |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")

    results = []
    for label, fname, f86_kmax, f86_qpeak in ORBITS:
        _, S_pert = load_S(DATA / fname)
        K_t = (S_pert - S_base) / DELTA_J
        peak_idx = int(np.argmax(np.abs(K_t)))
        results.append((label, K_t, peak_idx, f86_kmax, f86_qpeak))
        sample_idx = [20, 50, 100, 150, 200, 300]
        sample_K = " | ".join(f"{K_t[i]:+.5f}" for i in sample_idx)
        print(f"| {label} | {sample_K} | {abs(K_t[peak_idx]):.5f} | {t[peak_idx]:.2f} | {f86_kmax:.4f} | {f86_qpeak:.2f} |")

    print()
    print("# Orbit ranking comparison (largest first)")
    print()
    sorted_bi = sorted(results, key=lambda r: -abs(r[1][r[2]]))
    sorted_f86 = sorted(results, key=lambda r: -r[3])

    print("bond-isolate ρ_cc K_b orbit ranking at Q=1.5:")
    for r in sorted_bi:
        print(f"  {r[0]}: peak |K_b| = {abs(r[1][r[2]]):.5f} at t = {t[r[2]]:.2f}")

    print()
    print("F86 Dicke |K|max orbit ranking (each at own Q_peak):")
    for r in sorted_f86:
        print(f"  {r[0]}: F86 |K|max = {r[3]:.4f} at Q_peak = {r[4]:.2f}")

    print()
    print("# Order match? bond-isolate vs F86:")
    bi_order = [r[0] for r in sorted_bi]
    f86_order = [r[0] for r in sorted_f86]
    for i, (a, b) in enumerate(zip(bi_order, f86_order)):
        marker = "==" if a == b else "≠"
        print(f"  rank {i+1}: bond-isolate={a:15s} {marker} F86={b}")


if __name__ == "__main__":
    main()
