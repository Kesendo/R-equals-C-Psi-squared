"""Q=1.5 vs Q=2.0 split test at the canonical γ₀=0.05 substrate scale.

The 2026-05-19 N=9 chain sprint ran at the deep-quantum convention (γ=0.5, J=1,
Q = J/γ = 2) which happens to sit exactly on the Q_EP idealized anchor (per
`docs/Q_REGIME_ANCHORS.md` row Q=2.0). At that anchor, the star topology
saturates the bound `MaxImag ≤ J·N/2` because J = 2γ ↔ J·N/2 = σ = N·γ.

If we move the same four topologies to the framework's clean substrate
convention γ₀ = 0.05 with two J-values:
  Q = 1.5: J = 0.075 (F86 K_CC_pr Peak, c=2 saturated anchor)
  Q = 2.0: J = 0.100 (Q_EP idealized; reproduces today's Q=2 results,
           rescaled by γ ratio 1/10)

we get the following structural separation:
  - Star Im/σ at Q=2 should saturate to 1.0 exactly (J·N/2 = σ).
  - Star Im/σ at Q=1.5 should drop to ≤ 0.75 (J·N/2 = 0.75·σ).
  - The Q-invariant content (palindrome pairing, kernel dim, F1 identity) is
    identical between the two γ-conventions at the same Q.
  - The gap × N² constant of the chain is Q-dependent: at Q=2 with γ=0.5 it
    sat at 2.20 (today's anchor); at γ=0.05 it should rescale by γ² to 0.022,
    and at Q=1.5 the absolute number is whatever it is (we measure it).

Hamiltonian: H = (J/4) Σ_b (X_i X_j + Y_i Y_j + Z_i Z_j) on chain / ring / star,
uniform Z-dephasing γ₀ = 0.05 per site.

Output: simulations/results/q_split_anchor/{topology}_N{N}_Q{Q}.json plus a
comparison-table dump on stdout.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, "simulations")

from _f1_topology_heisenberg_small_n_anchor import (  # noqa: E402
    chain_bonds,
    ring_bonds,
    run,
    star_bonds,
)


GAMMA_SUBSTRATE = 0.05
Q_VALUES = (1.5, 2.0)
N_VALUES = (3, 4, 5, 6)
TOPOLOGIES = (
    ("chain", chain_bonds),
    ("ring", ring_bonds),
    ("star", star_bonds),
)


def main() -> None:
    print("F1 Q-split anchor: Q=1.5 vs Q=2.0 at γ₀ = 0.05")
    print("=" * 100)
    print(f"Substrate γ₀ = {GAMMA_SUBSTRATE}; J = Q · γ₀")
    print(f"Q=1.5 → J={1.5*GAMMA_SUBSTRATE} (F86 K_CC_pr Peak c=2 anchor)")
    print(f"Q=2.0 → J={2.0*GAMMA_SUBSTRATE} (Q_EP idealized; star-saturation anchor)")
    print()

    all_results: dict[tuple[str, int, float], dict] = {}
    out_dir = Path("simulations/results/q_split_anchor")
    out_dir.mkdir(parents=True, exist_ok=True)

    for Q in Q_VALUES:
        J = Q * GAMMA_SUBSTRATE
        print(f"\n--- Q = {Q} (J = {J}, γ = {GAMMA_SUBSTRATE}) ---")
        for topology, bonds_fn in TOPOLOGIES:
            for N in N_VALUES:
                if topology == "ring" and N < 3:
                    continue
                bonds = bonds_fn(N)
                print(f"  {topology} N={N}...", end=" ", flush=True)
                t0 = time.time()
                metrics = run(topology, bonds, N, J=J, gamma=GAMMA_SUBSTRATE)
                metrics["Q"] = Q
                print(f"{time.time()-t0:.2f}s, "
                      f"Im/σ={metrics['MaxImag']/(N*GAMMA_SUBSTRATE):.4f}, "
                      f"gap·N²={metrics['DissipationGap']*N*N:.6f}")
                all_results[(topology, N, Q)] = metrics
                out_path = out_dir / f"{topology}_N{N}_Q{Q:.1f}.json"
                with open(out_path, "w", encoding="utf-8") as f:
                    json.dump(metrics, f, indent=2)

    print("\n" + "=" * 100)
    print("Comparison table: dimensionless metrics at the same Q across γ-conventions")
    print()
    print(f"{'Topology':>8} {'N':>3} {'Q':>4} | "
          f"{'MinReal':>12} {'-2σ':>10} | "
          f"{'MaxImag':>10} {'Im/σ':>8} | "
          f"{'gap':>12} {'gap·N²':>10} {'gap·N²/γ²':>12} | "
          f"{'KernDim':>8} {'MaxPair':>10}")
    print("-" * 150)
    for topology, _ in TOPOLOGIES:
        for N in N_VALUES:
            for Q in Q_VALUES:
                if (topology, N, Q) not in all_results:
                    continue
                m = all_results[(topology, N, Q)]
                sigma = N * GAMMA_SUBSTRATE
                gap = m["DissipationGap"]
                print(f"{topology:>8} {N:3d} {Q:4.2f} | "
                      f"{m['MinReal']:12.6f} {-2*sigma:10.4f} | "
                      f"{m['MaxImag']:10.6f} {m['MaxImag']/sigma:8.4f} | "
                      f"{gap:12.7f} {gap*N*N:10.6f} {gap*N*N/(GAMMA_SUBSTRATE**2):12.4f} | "
                      f"{m['KernelDimension']:8d} {m['MaxPairingDistance']:10.3e}")
        print()

    print("=" * 100)
    print("Q-invariance check (Q=2 rescaled vs today's γ=0.5/J=1 anchors):")
    print("    Today's chain gap·N² at γ=0.5, J=1, Q=2: 2.18 (N=4), 2.21 (N=5), 2.18 (N=6)")
    print("    Q-invariant statement: gap·N² / γ² is the same dimensionless ratio.")
    print(f"    At γ={GAMMA_SUBSTRATE}: predicted gap·N²/γ² ≈ 2.18/0.25 ≈ 8.72")
    print()
    print("STAR Im/σ saturation check (at Q=2 only):")
    print("    Q=2: star Im/σ should = 1.000 (saturation, J·N/2 = σ)")
    print("    Q=1.5: star Im/σ should ≤ 0.750 (J·N/2 = 0.75·σ)")


if __name__ == "__main__":
    main()
