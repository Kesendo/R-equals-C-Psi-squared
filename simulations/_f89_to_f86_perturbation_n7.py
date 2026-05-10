"""F89 -> F86 bridge test: numerical ∂S/∂J_b at one bond around uniform-J full chain.

Setup:
  - Baseline: N=7 full chain, all bonds at J = J0 = 0.075, ρ_cc probe.
  - Perturbed: same as baseline EXCEPT one bond at J0 + δJ, δJ = 0.01.
  - Compute K_b(t) ≈ (S_perturbed(t) − S_baseline(t)) / δJ.
  - Peak |K_b(t)| over t per orbit representative (b = 0 Endpoint, b = 1 flank, b = 2 Center-near).

Comparison target: F86's K_b(Q=1.5, t_peak) per orbit at N=7 c=2 extended grid.
F86 orbit K_max values (Dicke probe, ResonanceScan):
  Endpoint  b=0: |K|max = 0.168 at Q_peak = 2.57
  flank     b=1: |K|max = 0.149 at Q_peak = 7.29
  Center    b=2: |K|max = 0.100 at Q_peak = 1.51

We are sampling at Q = J0/γ = 0.075/0.05 = 1.5 for ALL bonds (single Q point);
F86's |K|max ranking is for each bond at its OWN Q_peak (different Q per orbit).
The bond-isolate ranking we compute here is K_b at one shared Q (= 1.5), comparable to
F86's K_b(Q=1.5) per orbit — NOT to F86's |K|max-per-orbit.

What we test: does the orbit-ranking of K_b(Q=1.5) on ρ_cc match F86's K_b(Q=1.5) on Dicke?
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
    ("Endpoint  b=0", "N7_full_pert_b0_dJ0.01.csv"),
    ("flank     b=1", "N7_full_pert_b1_dJ0.01.csv"),
    ("Center    b=2", "N7_full_pert_b2_dJ0.01.csv"),
]


def load_S(path: Path) -> tuple[np.ndarray, np.ndarray]:
    arr = np.loadtxt(path, delimiter=",", skiprows=1)
    return arr[:, 0], arr[:, -1]


def main() -> None:
    t, S_base = load_S(DATA / "N7_full_baseline.csv")

    print(f"\n# F89 → F86 bridge: numerical ∂S/∂J_b on ρ_cc probe (N=7 full chain, J0=0.075, γ=0.05, δJ={DELTA_J})\n")
    print(f"# Baseline at J0 uniform: S(0) = {S_base[0]:.6f}, S(t=10) = {S_base[100]:.6f}, S(t=20) = {S_base[200]:.6f}")
    print()
    print("| Orbit | K_b(t=0) | K_b(t=5) | K_b(t=10) | K_b(t=15) | K_b(t=20) | K_b(t=30) | peak |K_b(t)| | t at peak |")
    print("|---|---|---|---|---|---|---|---|---|")

    results = []
    for label, fname in ORBITS:
        _, S_pert = load_S(DATA / fname)
        K_t = (S_pert - S_base) / DELTA_J
        peak_idx = int(np.argmax(np.abs(K_t)))
        results.append((label, K_t, peak_idx))
        sample_idx = [0, 50, 100, 150, 200, 300]
        sample_K = " | ".join(f"{K_t[i]:+.4f}" for i in sample_idx)
        print(f"| {label} | {sample_K} | {abs(K_t[peak_idx]):.4f} | {t[peak_idx]:.2f} |")

    print()
    print("# F86 K_b at N=7 c=2 (Dicke probe, extended grid) for comparison:")
    print("|  Orbit  | F86 |K|max | F86 Q_peak |")
    print("|---|---|---|")
    print("| Endpoint b=0 | 0.1684 | 2.57 |")
    print("| flank b=1    | 0.1489 | 7.29 |")
    print("| Center b=2   | 0.1000 | 1.51 |")
    print()
    print("# Orbit ranking comparison")
    print()
    sorted_bond_isolate = sorted(results, key=lambda r: -abs(r[1][r[2]]))
    print("bond-isolate ρ_cc K_b orbit ranking (largest peak |K_b(t)| first):")
    for r in sorted_bond_isolate:
        print(f"  {r[0]}: peak |K_b| = {abs(r[1][r[2]]):.4f}")

    print()
    print("F86 Dicke |K|max orbit ranking:")
    print("  Endpoint  b=0: |K|max = 0.1684")
    print("  flank     b=1: |K|max = 0.1489")
    print("  Center    b=2: |K|max = 0.1000")


if __name__ == "__main__":
    main()
