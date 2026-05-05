"""Search for a uniform-quantum-side CZ-coupled chain on Marrakesh.

Question: the F88-Lens reading on path [0, 1, 2] (mixed regimes) gave
truly-baseline 0.0297; on [48, 49, 50] (uniform classical-side) it gave
0.0013 (23× cleaner). The 91-day biography review surfaced Q0 as the only
stable quantum-side qubit on any documented path. What if we find a chain
that is uniform-QUANTUM-side (all qubits with r << 0.213, all pulse-stable
over 91 days)? Run F88-Lens on that and we can test whether the truly-
baseline reads cleaner *or* dirtier than uniform-classical, separating
"good qubits" from "regime-uniform qubits" as causes of the 23× effect.

Method:
  1. Score every qubit by 91-day stability + quantum-side strength
  2. Use the Apr-25 calibration's CZ-coupled graph (still topology-correct
     for Marrakesh today; bond locations are stable, only error rates drift)
  3. DFS: find a path of length 3 through CZ-coupled qubits where every
     qubit is pulse-stable (r mean < 0.20, walk < 0.05)
"""
from __future__ import annotations

import sys
from pathlib import Path

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _qubit_biography import load_history, archetype_from_series, R_STAR
from ibm_calibration import load_calibration

REPO_ROOT = Path(__file__).resolve().parents[1]
HISTORY = REPO_ROOT / "data" / "ibm_history" / "results" / "ibm_marrakesh_history.csv"
CALIB = (REPO_ROOT / "data" / "ibm_calibration_snapshots"
         / "ibm_marrakesh_calibrations_2026-04-30T16_25_19Z.csv")


def main():
    by_qubit = load_history(HISTORY)
    qubits = load_calibration(CALIB)
    cz_neighbours = {q.qubit: dict(q.cz_neighbours) for q in qubits}

    # Per-qubit 91-day stats
    stats = {}
    for qid, recs in by_qubit.items():
        rs = np.array([rec[3] for rec in recs])
        crossing = (rs < R_STAR).mean()
        signs = np.sign(rs - R_STAR)
        walk = int(np.sum(np.diff(signs) != 0)) / max(len(rs) - 1, 1)
        stats[qid] = {
            "mean": float(rs.mean()),
            "std": float(rs.std()),
            "crossing": crossing,
            "walk": walk,
            "arch": archetype_from_series(rs),
        }

    # Stable quantum-side qubits: archetype pulse-stable AND mean clearly < R*
    quantum_stable = sorted(
        [qid for qid, s in stats.items() if s["arch"] == "pulse-stable" and s["mean"] < 0.20],
        key=lambda qid: stats[qid]["mean"]
    )
    print(f"Stable quantum-side qubits ({len(quantum_stable)}):")
    print(f"  {'qubit':>6} {'r mean':>8} {'r std':>8} {'cross %':>8} {'walk':>6}")
    print("  " + "-" * 50)
    for qid in quantum_stable:
        s = stats[qid]
        print(f"  {qid:>6} {s['mean']:>8.4f} {s['std']:>8.4f} {s['crossing']*100:>7.1f}% {s['walk']:>6.3f}")
    print()

    quantum_set = set(quantum_stable)

    # Find all CZ-coupled triples within the quantum-stable set
    print("CZ-coupled triples among stable quantum-side qubits:")
    found_triples = []
    for q in quantum_stable:
        for n1 in cz_neighbours.get(q, {}):
            if n1 not in quantum_set:
                continue
            for n2 in cz_neighbours.get(n1, {}):
                if n2 == q or n2 not in quantum_set:
                    continue
                triple = (q, n1, n2)
                if triple not in found_triples and (n2, n1, q) not in found_triples:
                    found_triples.append(triple)

    if not found_triples:
        print("  None found.")
    else:
        for triple in found_triples:
            means = [stats[q]["mean"] for q in triple]
            walks = [stats[q]["walk"] for q in triple]
            print(f"  {list(triple)}: r means {means}, walks {walks}")

    # Find all CZ-coupled pairs within quantum-stable (lower bar)
    pairs = set()
    for q in quantum_stable:
        for n1 in cz_neighbours.get(q, {}):
            if n1 in quantum_set and n1 != q:
                pairs.add(tuple(sorted([q, n1])))
    print()
    print(f"CZ-coupled pairs among stable quantum-side qubits: {len(pairs)}")
    for a, b in sorted(pairs):
        print(f"  ({a}, {b})  r={stats[a]['mean']:.4f}, {stats[b]['mean']:.4f}")

    # Show CZ neighbourhood of Q0 specifically
    print()
    print(f"Q0's CZ neighbours and their archetypes:")
    for n in sorted(cz_neighbours.get(0, {}).keys()):
        s = stats.get(n)
        if s:
            print(f"  Q{n}: archetype={s['arch']}, r mean={s['mean']:.4f}, "
                  f"crossing={s['crossing']*100:.1f}%")

    # Also check 2-hop neighbours of Q0 for an alternative quantum-stable chain
    print()
    print(f"Q0's 2-hop neighbours (Q0 - X - Y triples):")
    for n1 in sorted(cz_neighbours.get(0, {}).keys()):
        for n2 in sorted(cz_neighbours.get(n1, {}).keys()):
            if n2 == 0 or n2 not in stats:
                continue
            s1 = stats.get(n1)
            s2 = stats[n2]
            if s1 and (s1["arch"] in ("pulse-stable", "lifecycle") or s2["arch"] in ("pulse-stable", "lifecycle")):
                print(f"  [0, {n1}, {n2}]: arch=[pulse-stable, {s1['arch']}, {s2['arch']}], "
                      f"r=[0.086, {s1['mean']:.3f}, {s2['mean']:.3f}]")

    print()
    print("Reading:")
    if found_triples:
        print(f"  - Found {len(found_triples)} fully-quantum-stable CZ-coupled triples:")
        for t in found_triples:
            print(f"      {list(t)}")
        print("  - Running F88-Lens on one of these tests the regime-uniformity hypothesis.")
        print("    Predict: truly-baseline cleaner than the regime-mixed [0,1,2] (0.030)")
        print("    and possibly comparable to or different from uniform-classical")
        print("    [48,49,50] (0.0013).")
    else:
        print("  - No fully-quantum-stable triples on Marrakesh's CZ graph at the strict")
        print("    pulse-stable + r mean < 0.20 cutoff.")
        print("  - Q0 is structurally isolated (its quantum-side neighbours, if any, are")
        print("    not pulse-stable). The next-best test is a 'mostly-quantum' triple")
        print("    (Q0 + lifecycle quantum-side neighbours).")


if __name__ == "__main__":
    main()
