"""Pre-flight audit for the next Marrakesh hardware run, using the 2026-05-05
calibration snapshot (ibm_marrakesh_calibrations_2026-05-05T04_18_17Z.csv) plus
the 91-day history through 2026-05-05.

Threefold view:

1. **Today's best chains** by composite score on the May-5 snapshot. Compares
   to the Apr-25 and Apr-30 winners to flag whether the same chain wins
   today, or whether the topography has shifted.

2. **Per-named-path summary** across all three snapshots: framework_snapshots
   [0, 1, 2], soft_break/zn_mirror [48, 49, 50], the Apr-25 / Apr-30 best
   chains, plus today's best. Score, regime composition, F87 confirmations
   already on the path.

3. **Drift verdict from 91-day history** for each candidate path. Combines
   snapshot regime (RegimeSummary equivalent) with multi-day archetype
   (LifecycleSummary equivalent) so we know which candidates are not just
   high-scoring today but also stable over the experiment-window timescale.

Goal: pick one path for the next QPU run. Both score and stability matter;
the May-5 snapshot alone can mis-rank a recently-flipped cluster.
"""
from __future__ import annotations

import sys
from pathlib import Path

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ibm_calibration import (
    load_calibration, score_qubit, best_chain, chain_score,
)
from _qubit_biography import load_history, archetype_from_series, R_STAR

REPO_ROOT = Path(__file__).resolve().parents[1]
SNAP_DIR = REPO_ROOT / "data" / "ibm_calibration_snapshots"
SNAP_APR25 = SNAP_DIR / "ibm_marrakesh_calibrations_2026-04-25T11_28_00Z.csv"
SNAP_APR30 = SNAP_DIR / "ibm_marrakesh_calibrations_2026-04-30T16_25_19Z.csv"
SNAP_MAY05 = SNAP_DIR / "ibm_marrakesh_calibrations_2026-05-05T04_18_17Z.csv"
HISTORY = REPO_ROOT / "data" / "ibm_history" / "results" / "ibm_marrakesh_history.csv"

NAMED_PATHS = [
    ("framework_snapshots [0, 1, 2]", [0, 1, 2]),
    ("soft_break / zn_mirror [48, 49, 50]", [48, 49, 50]),
    ("Apr-25 best 3-chain [4, 3, 2]", [4, 3, 2]),
    ("Apr-25 best 5-chain [1, 2, 3, 4, 5]", [1, 2, 3, 4, 5]),
    ("Apr-30 best 3-chain [93, 94, 95]", [93, 94, 95]),
    ("Apr-30 best 5-chain [95, 94, 93, 92, 91]", [95, 94, 93, 92, 91]),
]


def regime_count(qubits, path):
    by_id = {q.qubit: q for q in qubits}
    q = sum(1 for qid in path if by_id[qid].t2_us / (2.0 * by_id[qid].t1_us) < R_STAR)
    return q, len(path) - q


def archetype_for_path(history, path):
    archs = []
    for qid in path:
        if qid not in history:
            archs.append("missing")
            continue
        rs = np.array([rec[3] for rec in history[qid]])
        archs.append(archetype_for_series_str(rs))
    return archs


def archetype_for_series_str(rs):
    return archetype_from_series(rs)


def main():
    print("Marrakesh pre-flight audit, 2026-05-05")
    print("=" * 78)

    snap25 = load_calibration(SNAP_APR25)
    snap30 = load_calibration(SNAP_APR30)
    snap05 = load_calibration(SNAP_MAY05)
    history = load_history(HISTORY)

    print()
    print(f"Snapshots: Apr-25 ({len(snap25)} qubits), Apr-30 ({len(snap30)}), May-05 ({len(snap05)})")
    print(f"History:   {len(history)} qubits, up to {max(len(v) for v in history.values())} days each")

    # --- Section 1: today's best chains ---
    print()
    print("=" * 78)
    print("1. TODAY'S BEST CHAINS (by composite score on May-5 snapshot)")
    print()
    p3, s3 = best_chain(snap05, 3)
    p5, s5 = best_chain(snap05, 5)
    print(f"  Best 3-chain:  {p3}, score = {s3:.2f}")
    print(f"  Best 5-chain:  {p5}, score = {s5:.2f}")
    print()
    print(f"  Apr-25 had:    [4, 3, 2] (score 867.07), [1, 2, 3, 4, 5] (1246.58)")
    print(f"  Apr-30 had:    [93, 94, 95] (917.06), [95, 94, 93, 92, 91] (1469.23)")

    # Add today's winners to the named paths if not already present
    today_paths = []
    if p3 not in [path for _, path in NAMED_PATHS]:
        today_paths.append((f"May-05 best 3-chain {p3}", list(p3)))
    if p5 not in [path for _, path in NAMED_PATHS]:
        today_paths.append((f"May-05 best 5-chain {p5}", list(p5)))
    all_paths = NAMED_PATHS + today_paths

    # --- Section 2: per-path scores across snapshots ---
    print()
    print("=" * 78)
    print("2. PER-PATH SCORES ACROSS THREE SNAPSHOTS")
    print()
    print(f"  {'path':<45} {'Apr-25':>10} {'Apr-30':>10} {'May-05':>10}  trend")
    print("  " + "-" * 88)
    for label, path in all_paths:
        s_apr25 = chain_score(snap25, path)
        s_apr30 = chain_score(snap30, path)
        s_may05 = chain_score(snap05, path)
        trend = ""
        if s_may05 > max(s_apr25, s_apr30):
            trend = "↑ rising"
        elif s_may05 < min(s_apr25, s_apr30):
            trend = "↓ falling"
        else:
            trend = "≈"
        print(f"  {label:<45} {s_apr25:>10.2f} {s_apr30:>10.2f} {s_may05:>10.2f}  {trend}")

    # --- Section 3: regime + lifecycle audit per path ---
    print()
    print("=" * 78)
    print("3. REGIME (May-5 snapshot) + LIFECYCLE (91-day history) PER PATH")
    print()
    print(f"  {'path':<45} {'regime':>22} {'archetypes (per qubit)'}")
    print("  " + "-" * 100)
    for label, path in all_paths:
        q, c = regime_count(snap05, path)
        regime_label = (f"uniform-quantum ({q}/{len(path)})" if q == len(path) else
                        f"uniform-classical ({c}/{len(path)})" if c == len(path) else
                        f"mixed ({q}q/{c}c)")
        archs = archetype_for_path(history, path)
        archs_short = [a.replace("-stable", "").replace("drifty-", "drifty") for a in archs]
        print(f"  {label:<45} {regime_label:>22}  {', '.join(archs_short)}")

    # --- Section 4: recommendation ---
    print()
    print("=" * 78)
    print("4. RECOMMENDATION")
    print()
    # Score + stability composite
    best_combo = None
    best_combo_score = -1
    for label, path in all_paths:
        s = chain_score(snap05, path)
        archs = archetype_for_path(history, path)
        n_stable = sum(1 for a in archs if a in ("pulse-stable", "silent-stable", "classic-stable"))
        n_twitch = sum(1 for a in archs if a == "twitch")
        # heuristic: prefer addressable-on-history (no missing), penalize twitch heavily,
        # weigh stable archetypes positively, then score as tiebreaker
        if "missing" in archs:
            continue
        composite = (n_stable * 100) - (n_twitch * 200) + s
        if composite > best_combo_score:
            best_combo_score = composite
            best_combo = (label, path, s, archs)

    if best_combo:
        label, path, s, archs = best_combo
        print(f"  → {label}")
        print(f"    May-05 score: {s:.2f}")
        print(f"    91-day archetypes: {', '.join(archs)}")
        print(f"    Composite (score + stability bonus - twitch penalty): {best_combo_score:.2f}")
    print()
    print("  Reading: a high May-5 score on a recently-flipped cluster (Apr-30 winners)")
    print("  may not be the right pick if the 91-day archetypes show twitch or lifecycle")
    print("  drift. Stable + addressable + decent score beats bleeding-edge score on a")
    print("  drift-volatile path.")


if __name__ == "__main__":
    main()
