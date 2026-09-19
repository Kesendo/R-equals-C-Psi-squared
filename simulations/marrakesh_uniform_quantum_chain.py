"""Search for a CZ-coupled chain whose calibration rows are all below R*.

Question: the F88b-Lens reading on path [0, 1, 2] (mixed R* bands) gave
truly-baseline 0.0297; on [48, 49, 50] (all at-or-above-R*) it gave 0.0013.
The 91-day proxy biography surfaced Q0 as the only stable-below history on a
documented path. Finding an all-below-R* triple permits another finite
comparison, but backend, path, date, and calibration all remain confounded;
the 23× association is not a causal band effect.

Method:
  1. Score every qubit by 91-day stability + distance below R*
  2. Use the Apr-25 calibration's CZ-coupled graph (still topology-correct
     for Marrakesh today; bond locations are stable, only error rates drift)
  3. DFS: find a path of length 3 through CZ-coupled qubits where every
     qubit is stable-below (r mean < 0.20, walk < 0.05)
"""
from __future__ import annotations

import sys
from pathlib import Path

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qubit_biography import load_history, archetype_from_series, R_STAR
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
        below_flags = rs < R_STAR
        below_fraction = below_flags.mean()
        walk = int(np.sum(below_flags[1:] != below_flags[:-1])) / max(len(rs) - 1, 1)
        stats[qid] = {
            "mean": float(rs.mean()),
            "std": float(rs.std()),
            "below_fraction": below_fraction,
            "walk": walk,
            "arch": archetype_from_series(rs),
        }

    # Stable below-R* histories with a conservative margin.
    below_rstar_stable = sorted(
        [qid for qid, s in stats.items() if s["arch"] == "stable-below" and s["mean"] < 0.20],
        key=lambda qid: stats[qid]["mean"]
    )
    print(f"Stable below-R* histories ({len(below_rstar_stable)}):")
    print(f"  {'qubit':>6} {'r mean':>8} {'r std':>8} {'below %':>8} {'walk':>6}")
    print("  " + "-" * 50)
    for qid in below_rstar_stable:
        s = stats[qid]
        print(f"  {qid:>6} {s['mean']:>8.4f} {s['std']:>8.4f} {s['below_fraction']*100:>7.1f}% {s['walk']:>6.3f}")
    print()

    below_rstar_set = set(below_rstar_stable)

    # Find all CZ-coupled triples within the stable below-R* set.
    print("CZ-coupled triples among stable below-R* histories:")
    found_triples = []
    for q in below_rstar_stable:
        for n1 in cz_neighbours.get(q, {}):
            if n1 not in below_rstar_set:
                continue
            for n2 in cz_neighbours.get(n1, {}):
                if n2 == q or n2 not in below_rstar_set:
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

    # Find all CZ-coupled pairs within the stable below-R* set (lower bar).
    pairs = set()
    for q in below_rstar_stable:
        for n1 in cz_neighbours.get(q, {}):
            if n1 in below_rstar_set and n1 != q:
                pairs.add(tuple(sorted([q, n1])))
    print()
    print(f"CZ-coupled pairs among stable below-R* histories: {len(pairs)}")
    for a, b in sorted(pairs):
        print(f"  ({a}, {b})  r={stats[a]['mean']:.4f}, {stats[b]['mean']:.4f}")

    # Show CZ neighbourhood of Q0 specifically
    print()
    print(f"Q0's CZ neighbours and their archetypes:")
    for n in sorted(cz_neighbours.get(0, {}).keys()):
        s = stats.get(n)
        if s:
            print(f"  Q{n}: archetype={s['arch']}, r mean={s['mean']:.4f}, "
                  f"below_fraction={s['below_fraction']*100:.1f}%")

    # Also check two-hop neighbours of Q0 for an alternative below-R* chain.
    print()
    print(f"Q0's 2-hop neighbours (Q0 - X - Y triples):")
    for n1 in sorted(cz_neighbours.get(0, {}).keys()):
        for n2 in sorted(cz_neighbours.get(n1, {}).keys()):
            if n2 == 0 or n2 not in stats:
                continue
            s1 = stats.get(n1)
            s2 = stats[n2]
            if s1 and (s1["arch"] in ("stable-below", "multi-band") or s2["arch"] in ("stable-below", "multi-band")):
                print(f"  [0, {n1}, {n2}]: arch=[stable-below, {s1['arch']}, {s2['arch']}], "
                      f"r=[0.086, {s1['mean']:.3f}, {s2['mean']:.3f}]")

    print()
    print("Reading:")
    if found_triples:
        print(f"  - Found {len(found_triples)} all-below-R* CZ-coupled triples:")
        for t in found_triples:
            print(f"      {list(t)}")
        print("  - These enable a finite comparison with mixed-band [0,1,2] (0.030)")
        print("    and all-at-or-above-R* [48,49,50] (0.0013).")
        print("    Backend/path/date differences keep the association confounded.")
    else:
        print("  - No all-below-R* triples on Marrakesh's CZ graph at the strict")
        print("    stable-below + r mean < 0.20 cutoff.")
        print("  - Q0 is structurally isolated from other stable-below histories.")
        print("    A mixed-band triple is a different, explicitly confounded comparison.")


if __name__ == "__main__":
    main()
