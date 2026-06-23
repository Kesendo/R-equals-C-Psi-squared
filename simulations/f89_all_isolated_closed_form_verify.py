"""Verify the F89 all-isolated subclass closed form against bond-isolate CSVs.

Theorem (F89 Tier 1 derived candidate, late-tail subclass):

  For the all-isolated topology (1)^m on N qubits (m disjoint NN-bonds,
  N − 2m bare sites), uniform J coupling, uniform Z-dephasing γ₀, and
  ρ_cc initial state, the spatial-sum coherence has the EXACT closed form

    S_(1)^m, N(t) = [(N−1)/N + 4m(N−2)(cos(4Jt) − 1)/(N²(N−1))] · exp(−4γ₀ t)

  The asymptotic decay rate is 4γ₀ universal across all m (matching F73's
  vac-SE rate). The m-dependence enters only through a periodic correction
  with frequency 4J that vanishes at cos(4Jt) = 1.

This script loads the (1), (1,1), (1,1,1) CSVs from the bond-isolate
suite and overlays the closed-form prediction. Bit-exact match expected
modulo RK4 numerical noise.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO = Path(__file__).resolve().parents[1]
CSV_DIR = REPO / "simulations" / "results" / "bond_isolate"
PARAMS = "J0.0750_gamma0.0500_probe-coherence"

N = 7
J = 0.075
GAMMA = 0.05


def closed_form(n: int, m: int, j: float, gamma: float, t: np.ndarray) -> np.ndarray:
    base = (n - 1) / n
    correction = 4 * m * (n - 2) * (np.cos(4 * j * t) - 1) / (n**2 * (n - 1))
    return (base + correction) * np.exp(-4 * gamma * t)


def load_S(bonds: list[int]) -> tuple[np.ndarray, np.ndarray]:
    bs = sorted(bonds)
    label = f"b{bs[0]}" if len(bs) == 1 else "b" + "-".join(str(b) for b in bs)
    p = CSV_DIR / f"N{N}_{label}_{PARAMS}.csv"
    data = np.loadtxt(p, delimiter=",", skiprows=1)
    return data[:, 0], data[:, -1]


def main() -> None:
    cases = [
        ("(1) m=1", 1, [0]),
        ("(1,1) m=2", 2, [0, 2]),
        ("(1,1,1) m=3", 3, [0, 2, 4]),
    ]

    print()
    print(f"# Closed-form verification: S_(1)^m, N={N} (J={J}, γ={GAMMA})")
    print()

    for label, m, bonds in cases:
        t_data, S_data = load_S(bonds)
        S_pred = closed_form(N, m, J, GAMMA, t_data)
        max_abs = float(np.max(np.abs(S_data - S_pred)))
        max_rel = float(np.max(np.abs(S_data - S_pred) / np.maximum(np.abs(S_data), 1e-15)))
        print(f"  {label}: max |S_data − S_pred| = {max_abs:.3e}, max rel = {max_rel:.3e}")

    print()
    print("# Spot-check at sample times")
    print()
    sample_t = np.array([0.0, 5.0, 10.0, 20.0, 21.0, 30.0])
    print(f"  cos(4Jt) at sample times: {np.cos(4*J*sample_t)}")
    print()
    for label, m, _ in cases:
        S_pred = closed_form(N, m, J, GAMMA, sample_t)
        print(f"  {label}: S(t) = {[f'{s:.6f}' for s in S_pred]}")

    print()
    print(f"Asymptotic rate Γ_∞ = 4γ₀ = {4*GAMMA:.4f}")
    print(f"Period of m-correction = π/(2J) = {np.pi/(2*J):.3f}")
    print(f"In-phase moments (cos(4Jt)=1) at t = k·π/(2J), k = 0, 1, 2, ...")


if __name__ == "__main__":
    main()
