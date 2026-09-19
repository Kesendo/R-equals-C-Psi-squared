"""Pull the 91-day biography for each qubit on our F88b-Lens experimental paths,
applying qubit_biography classifier to the freshly fetched Marrakesh history.

Cross-references:
- `marrakesh_quarter_boundary_review.py`: 5-day-window flipping (Apr-25 vs Apr-30)
- `qubit_biography.py`: full classifier (used here)
- `docs/BOTH_SIDES_VISIBLE.md`: original Torino 180-day analysis

Question: which path qubits have stable versus variable histories under the
free-single-transmon `|+>` normalized-purity proxy? The labels describe neutral
R* bands in finite calibration records, not physical regimes.
"""
from __future__ import annotations

import sys
from pathlib import Path

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qubit_biography import (
    load_history, biography, archetype_from_series, fingerprint, R_STAR
)

REPO_ROOT = Path(__file__).resolve().parents[1]
HISTORY = REPO_ROOT / "data" / "ibm_history" / "results" / "ibm_marrakesh_history.csv"

NAMED_PATHS = [
    ("framework_snapshots [0, 1, 2]", [0, 1, 2]),
    ("soft_break / zn_mirror [48, 49, 50]", [48, 49, 50]),
    ("Apr-25 best 3-chain [4, 3, 2]", [4, 3, 2]),
    ("Apr-25 best 5-chain [1, 2, 3, 4, 5]", [1, 2, 3, 4, 5]),
    ("Apr-30 best 3-chain [93, 94, 95]", [93, 94, 95]),
    ("Apr-30 best 5-chain [95, 94, 93, 92, 91]", [95, 94, 93, 92, 91]),
]


def main():
    print(f"Marrakesh 91-day biography for F88b-Lens experimental paths")
    print("=" * 78)
    by_qubit = load_history(HISTORY)

    all_path_qubits = sorted({qid for _, p in NAMED_PATHS for qid in p})

    print()
    print(f"  {'qubit':>6} {'archetype':>16} {'r mean':>8} {'r std':>8} "
          f"{'below %':>8} {'walk':>6}")
    print("  " + "-" * 70)
    qstats = {}
    for qid in all_path_qubits:
        if qid not in by_qubit:
            print(f"  Q{qid:<5} (no calibration data)")
            continue
        rs = np.array([rec[3] for rec in by_qubit[qid]])
        below_flags = rs < R_STAR
        below_fraction = below_flags.mean() * 100
        walk = int(np.sum(below_flags[1:] != below_flags[:-1])) / max(len(rs) - 1, 1)
        arch = archetype_from_series(rs)
        qstats[qid] = (arch, rs.mean(), rs.std(), below_fraction, walk)
        print(f"  {qid:>6} {arch:>16} {rs.mean():>8.4f} {rs.std():>8.4f} "
              f"{below_fraction:>7.1f}% {walk:>6.3f}")
    print()

    print("Path-level summary:")
    print(f"  {'path':<42} {'archetypes (n_qubits)':>40}")
    print("  " + "-" * 85)
    for label, path in NAMED_PATHS:
        archs = [qstats[q][0] for q in path if q in qstats]
        from collections import Counter
        counts = Counter(archs)
        summary = ", ".join(f"{a}={n}" for a, n in counts.most_common())
        print(f"  {label:<42} {summary:>40}")
    print()

    print("Per-qubit fingerprints (91-day strings):")
    print("  legend: . well-above  _ at-or-above  X below  ~ variable-near")
    print("          / rising-through  \\ falling-through  ! outlier")
    print()
    for qid in all_path_qubits:
        if qid not in by_qubit:
            continue
        rs = np.array([rec[3] for rec in by_qubit[qid]])
        phases = biography(rs)
        fp = fingerprint(phases)
        arch = qstats[qid][0]
        print(f"  Q{qid:<3} [{arch:<14}]")
        chunks = [fp[i:i+30] for i in range(0, len(fp), 30)]
        for chunk in chunks:
            print(f"    {chunk}")
        print()

    print("Reading:")
    print("  Path-quality on F88b-Lens timescale = (a) per-qubit archetype + (b) co-stability.")
    print("  Stable-below and stable-at-or-above are the lowest-walk R* histories.")
    print("  Variable-near histories warn that a later run may see another proxy band;")
    print("  this is a scheduling caveat, not a submission verdict.")


if __name__ == "__main__":
    main()
