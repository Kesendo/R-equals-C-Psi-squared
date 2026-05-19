"""Q-sweep over canonical anchors {0.5, 1.0, 1.5, √3, 2.0, 2.5} at γ₀=0.05.

Extends `_f1_q_split_anchor.py` (which only ran Q ∈ {1.5, 2.0}) across the six
table-anchored Q values from `docs/Q_REGIME_ANCHORS.md`:

  Q = 0.5  (sub-balance, J/γ < 1)
  Q = 1.0  (Balance, J = γ exactly)
  Q = 1.5  (F86 Q_peak c=2)
  Q = √3 ≈ 1.732  (canonical θ=60°, LindbladAbsorptionMatch)
  Q = 2.0  (Q_EP idealized)
  Q = 2.5  (Endpoint orbit candidate)

Goal: pin down three structural hypotheses surfaced by the Q=1.5 vs Q=2.0 split:

  H1 (star saturation, universal): Im_max(star, N, J) = J·N/2 for all Q.
      → Im/σ = Q/2 should hold flat across all 6 anchors, not just Q=2.

  H2 (ring N=4 dihedral lock): Im_max(ring, N=4) = (3/4)·J·N = 3·J → Im/σ = 3Q/4.
      → Test 6 Q-values to confirm the 3Q/4 line; if it holds, the bound is
        a typed claim candidate. N=6 ring (if it also locks) would give a
        second dihedral-symmetric witness.

  H3 (chain gap·N²/γ as a Q-function): The Q=2 plateau ≈ 4.36 and Q=1.5
      plateau ≈ 2.55 are two points on a smooth curve f(Q). 6 anchors give
      enough shape to guess whether it's linear, quadratic, or has a
      cusp/bend.

Output: simulations/results/q_sweep_anchor/{topology}_N{N}_Q{Q}.json
+ a comparison-table dump showing Im/σ vs Q/2 (for star), 3Q/4 (for ring N=4),
and chain gap·N²/γ vs Q.
"""
from __future__ import annotations

import json
import math
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
Q_ANCHORS = (
    (0.5, "sub-balance"),
    (1.0, "Balance"),
    (1.5, "F86 Q_peak c=2"),
    (math.sqrt(3.0), "canonical θ=60°"),
    (2.0, "Q_EP idealized"),
    (2.5, "Endpoint orbit"),
)
N_VALUES = (3, 4, 5, 6)
TOPOLOGIES = (
    ("chain", chain_bonds),
    ("ring", ring_bonds),
    ("star", star_bonds),
)


def main() -> None:
    print("F1 Q-sweep anchor: 6 Q-values × 3 topologies × N=3..6 at γ₀ = 0.05")
    print("=" * 110)
    print(f"Substrate γ₀ = {GAMMA_SUBSTRATE}; J = Q · γ₀")
    for Q, label in Q_ANCHORS:
        print(f"  Q = {Q:.4f}  J = {Q*GAMMA_SUBSTRATE:.6f}  ({label})")
    print()

    all_results: dict[tuple[str, int, float], dict] = {}
    out_dir = Path("simulations/results/q_sweep_anchor")
    out_dir.mkdir(parents=True, exist_ok=True)

    t_start = time.time()
    for Q, label in Q_ANCHORS:
        J = Q * GAMMA_SUBSTRATE
        print(f"\n--- Q = {Q:.4f} ({label}, J = {J:.6f}) ---")
        for topology, bonds_fn in TOPOLOGIES:
            for N in N_VALUES:
                if topology == "ring" and N < 3:
                    continue
                bonds = bonds_fn(N)
                print(f"  {topology} N={N}...", end=" ", flush=True)
                t0 = time.time()
                metrics = run(topology, bonds, N, J=J, gamma=GAMMA_SUBSTRATE)
                metrics["Q"] = Q
                metrics["QLabel"] = label
                sigma = N * GAMMA_SUBSTRATE
                print(f"{time.time()-t0:.2f}s, "
                      f"Im/σ={metrics['MaxImag']/sigma:.4f}, "
                      f"gap·N²/γ={metrics['DissipationGap']*N*N/GAMMA_SUBSTRATE:.4f}")
                all_results[(topology, N, Q)] = metrics
                out_path = out_dir / f"{topology}_N{N}_Q{Q:.4f}.json"
                with open(out_path, "w", encoding="utf-8") as f:
                    json.dump(metrics, f, indent=2)
    print(f"\nTotal sweep wall: {time.time()-t_start:.1f}s")

    # ----------------------------------------------------------------------
    # H1: star Im/σ should equal Q/2 universally
    # ----------------------------------------------------------------------
    print("\n" + "=" * 110)
    print("H1: STAR Im/σ vs Q/2 (predicted Im_max(star) = J·N/2 ∀ Q)")
    print("-" * 110)
    print(f"{'N':>3} | " + " | ".join(f"Q={Q:.4f}" for Q, _ in Q_ANCHORS))
    print(f"{'pred':>3} | " + " | ".join(f"{Q/2:.4f}  " for Q, _ in Q_ANCHORS))
    for N in N_VALUES:
        cells = []
        for Q, _ in Q_ANCHORS:
            if ("star", N, Q) in all_results:
                im_sigma = all_results[("star", N, Q)]["MaxImag"] / (N * GAMMA_SUBSTRATE)
                cells.append(f"{im_sigma:.4f}  ")
            else:
                cells.append("  --    ")
        print(f"{N:3d} | " + " | ".join(cells))

    # ----------------------------------------------------------------------
    # H2: ring N=4 Im/σ should equal 3Q/4
    # ----------------------------------------------------------------------
    print("\n" + "=" * 110)
    print("H2: RING N=4 Im/σ vs 3Q/4 (predicted Im_max(ring, N=4) = (3/4)·J·N = 3·J)")
    print("-" * 110)
    print(f"{'Q':>8} {'predicted 3Q/4':>15} {'observed Im/σ':>15} {'relative error':>16}")
    for Q, _ in Q_ANCHORS:
        if ("ring", 4, Q) in all_results:
            pred = 3.0 * Q / 4.0
            obs = all_results[("ring", 4, Q)]["MaxImag"] / (4 * GAMMA_SUBSTRATE)
            rel = abs(obs - pred) / pred if pred > 0 else float("nan")
            print(f"{Q:8.4f} {pred:15.6f} {obs:15.6f} {rel:16.3e}")

    # Bonus: ring N=6 dihedral lock check
    print("\n     Ring N=6 (testing for dihedral lock at higher even N):")
    print(f"{'Q':>8} {'Im/σ':>10} {'Im/σ / Q':>12}")
    for Q, _ in Q_ANCHORS:
        if ("ring", 6, Q) in all_results:
            obs = all_results[("ring", 6, Q)]["MaxImag"] / (6 * GAMMA_SUBSTRATE)
            print(f"{Q:8.4f} {obs:10.4f} {obs/Q:12.6f}")

    # ----------------------------------------------------------------------
    # H3: chain gap·N²/γ as a Q-function
    # ----------------------------------------------------------------------
    print("\n" + "=" * 110)
    print("H3: CHAIN gap·N²/γ vs Q (looking for f(Q) shape)")
    print("-" * 110)
    print(f"{'Q':>8} | " + " | ".join(f"N={N}" for N in N_VALUES))
    for Q, label in Q_ANCHORS:
        cells = []
        for N in N_VALUES:
            if ("chain", N, Q) in all_results:
                gap_n2_g = all_results[("chain", N, Q)]["DissipationGap"] * N * N / GAMMA_SUBSTRATE
                cells.append(f"{gap_n2_g:7.4f}")
            else:
                cells.append("  --   ")
        print(f"{Q:8.4f} | " + " | ".join(cells) + f"  ({label})")

    print("\n     Q=1.5 plateau (N=4..6 mean) and Q=2.0 plateau (N=4..6 mean):")
    for Q, _ in Q_ANCHORS:
        vals = [all_results[("chain", N, Q)]["DissipationGap"] * N * N / GAMMA_SUBSTRATE
                for N in (4, 5, 6) if ("chain", N, Q) in all_results]
        if vals:
            mean = sum(vals) / len(vals)
            print(f"     Q={Q:.4f}: plateau (N=4..6) mean = {mean:.4f}")

    # ----------------------------------------------------------------------
    # F1 palindrome Q-invariance check (sanity)
    # ----------------------------------------------------------------------
    print("\n" + "=" * 110)
    print("F1 palindrome Q-invariance: MinReal vs −2σ (predicted bit-exact)")
    print("-" * 110)
    bad = []
    for (topology, N, Q), m in all_results.items():
        sigma = N * GAMMA_SUBSTRATE
        if abs(m["MinReal"] - (-2 * sigma)) > 1e-10:
            bad.append((topology, N, Q, m["MinReal"], -2 * sigma))
    if not bad:
        print(f"     All {len(all_results)} (topology, N, Q) combinations: MinReal = −2σ to < 1e-10 ✓")
    else:
        print(f"     {len(bad)} violations:")
        for t, N, Q, mr, pred in bad:
            print(f"       {t} N={N} Q={Q:.4f}: MinReal={mr:.6f} vs −2σ={pred:.6f}")

    # MaxPairingDistance summary
    print()
    print("MaxPairingDistance summary by topology (worst across all (N, Q)):")
    for topology, _ in TOPOLOGIES:
        worst = max(((N, Q, m["MaxPairingDistance"])
                     for (t, N, Q), m in all_results.items() if t == topology),
                    key=lambda x: x[2])
        N_worst, Q_worst, d_worst = worst
        print(f"     {topology}: worst at N={N_worst} Q={Q_worst:.4f}: {d_worst:.3e}")


if __name__ == "__main__":
    main()
