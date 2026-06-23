"""Slow-mode sector sweep: ring + star + Q-sweep + N=7 chain.

Batched runner for tasks #69-#72:
  - Ring N=4..6 at Q=2 (Marrakesh convention, γ=0.5, J=1).
  - Star N=3..6 at Q=2.
  - Chain N=5 across 6 canonical Q-anchors {0.5, 1.0, 1.5, √3, 2.0, 2.5}.
  - Chain N=7 at Q=2 (last N before block-spectrum bridge required).

For each (topology, N, Q): slow-mode sector, ⟨n_XY⟩, weight distribution,
Absorption-Theorem cross-check. Predictions to test:

  - Ring: same central-popcount-block + magnon-admixture story as chain,
    4× larger gap prefactor (cyclic vs open k_min²).
  - Star: different sector geometry; slow mode may live in non-central
    popcount sector or have hub-localized weight rather than n_XY=2 admixture.
  - Q-sweep at fixed N: 0.55 coefficient stability test (Q-invariant or
    Q-dependent drift?).
  - N=7 chain: confirm central-popcount-block result at scale frontier.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, "simulations")

from slow_mode_sector_diagnostic import run  # noqa: E402


def main() -> None:
    print("=" * 100)
    print("Slow-mode sector sweep: ring + star + Q-sweep + N=7 chain")
    print("=" * 100)

    results = {}

    # Task #69: Ring N=4..6 at Q=2 (Marrakesh convention)
    print("\n[Task 69] RING N=4..6 at Q=2")
    print("-" * 100)
    for N in (4, 5, 6):
        results[("ring", N, 2.0)] = run("ring", N=N, J=1.0, gamma=0.5)

    # Task #70: Star N=3..6 at Q=2
    print("\n[Task 70] STAR N=3..6 at Q=2")
    print("-" * 100)
    for N in (3, 4, 5, 6):
        results[("star", N, 2.0)] = run("star", N=N, J=1.0, gamma=0.5)

    # Task #71: Q-sweep at chain N=5, 6 canonical Q anchors, γ₀=0.05
    print("\n[Task 71] CHAIN N=5 Q-sweep at γ₀=0.05")
    print("-" * 100)
    sqrt3 = math.sqrt(3.0)
    q_anchors = [
        (0.5, "sub-balance"),
        (1.0, "Balance"),
        (1.5, "F86 Q_peak c=2"),
        (sqrt3, "canonical θ=60°"),
        (2.0, "Q_EP idealized"),
        (2.5, "Endpoint orbit"),
    ]
    gamma_substrate = 0.05
    for Q, _label in q_anchors:
        J = Q * gamma_substrate
        results[("chain", 5, Q)] = run("chain", N=5, J=J, gamma=gamma_substrate)

    # Task #72: N=7 chain at Q=2 (Marrakesh convention)
    print("\n[Task 72] CHAIN N=7 at Q=2 (will take ~5-10s for Pauli projection)")
    print("-" * 100)
    results[("chain", 7, 2.0)] = run("chain", N=7, J=1.0, gamma=0.5)

    # Comparison tables
    print("\n" + "=" * 100)
    print("Summary: ring + star at Q=2 (γ=0.5, J=1)")
    print(f"{'topo':>5} {'N':>3}  {'gap':>10}  {'sector':>10}  {'⟨n_XY⟩':>10}  "
          f"{'w_0':>7}  {'w_2':>7}  {'w_4':>7}")
    for key, r in results.items():
        topology, N, Q = key
        if Q == 2.0 and topology in ("ring", "star"):
            wd = r["weight_distribution"]
            w0 = wd[0] if len(wd) > 0 else 0
            w2 = wd[2] if len(wd) > 2 else 0
            w4 = wd[4] if len(wd) > 4 else 0
            print(f"{topology:>5} {N:3d}  {r['gap']:10.6e}  {str(r['slow_sector']):>10}  "
                  f"{r['n_xy_avg']:10.6f}  {w0:7.4f}  {w2:7.4f}  {w4:7.4f}")

    print()
    print("Summary: Q-sweep at chain N=5 (γ=0.05)")
    print(f"{'Q':>8}  {'gap':>10}  {'sector':>10}  {'⟨n_XY⟩':>10}  "
          f"{'predicted 0.55·Q²/N²':>22}  {'ratio':>8}")
    for Q, _label in q_anchors:
        r = results.get(("chain", 5, Q))
        if r is None:
            continue
        pred = 0.55 * Q * Q / 25
        ratio = r["n_xy_avg"] / pred if pred > 0 else float("nan")
        print(f"{Q:8.4f}  {r['gap']:10.6e}  {str(r['slow_sector']):>10}  "
              f"{r['n_xy_avg']:10.6f}  {pred:22.6f}  {ratio:8.4f}")

    print()
    print("N=7 chain at Q=2:")
    r = results[("chain", 7, 2.0)]
    print(f"  sector = {r['slow_sector']}, gap = {r['gap']:.6e}, ⟨n_XY⟩ = {r['n_xy_avg']:.6f}")
    print(f"  predicted ⟨n_XY⟩ = 0.55·Q²/N² = {0.55*4/49:.6f}")
    print(f"  weight distribution: {[f'{w:.4f}' for w in r['weight_distribution']]}")


if __name__ == "__main__":
    main()
