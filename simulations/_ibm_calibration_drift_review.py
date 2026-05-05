"""Side-by-side review of the two Marrakesh calibration snapshots that bracket
the 2026-04-26 soft_break and 2026-04-29 zn_mirror runs:

- ibm_marrakesh_calibrations_2026-04-25T11_28_00Z.csv  (1 day before soft_break)
- ibm_marrakesh_calibrations_2026-04-30T16_25_19Z.csv  (1 day after zn_mirror)

Reports for each:
  - Top-10 qubits by composite score
  - Best 3-chain and 5-chain (DFS over CZ graph)
  - Per-named-path score: [0, 1, 2], [48, 49, 50], [4, 3, 2], [1, 2, 3, 4, 5]

Then prints the per-path drift across the four-day window: T2 and score
deltas show how much the calibration moved between the two experiment runs,
which directly bears on the F88-Lens truly-baseline noise floor reading
(0.0013 on path [48, 49, 50] for soft_break, 0.0102 for zn_mirror Heisenberg
on the same path; the path-quality didn't change much between runs, so the
~10× soft→Heisenberg gap is in the F88 reading, not the substrate).

Reading: with the corrected docstring (see commit b588e2d), the framing is
that calibration score is the *input*, F88-Lens state-level number is the
*consequence*. The score gap between paths is ~14% (linear in qubit
quality); the F88 truly-baseline gap is 23× because the state-level reading
amplifies per-qubit dephasing nonlinearly.
"""
from __future__ import annotations

import sys
from pathlib import Path

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ibm_calibration import (
    load_calibration,
    score_qubit,
    best_qubits,
    best_chain,
    chain_score,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
CALIB_DIR = REPO_ROOT / "data" / "ibm_calibration_snapshots"
CALIB_APR25 = CALIB_DIR / "ibm_marrakesh_calibrations_2026-04-25T11_28_00Z.csv"
CALIB_APR30 = CALIB_DIR / "ibm_marrakesh_calibrations_2026-04-30T16_25_19Z.csv"

NAMED_PATHS = [
    ("framework_snapshots [0, 1, 2]", [0, 1, 2]),
    ("soft_break / zn_mirror [48, 49, 50]", [48, 49, 50]),
    ("best 3-chain (Apr-25) [4, 3, 2]", [4, 3, 2]),
    ("best 5-chain (Apr-25) [1, 2, 3, 4, 5]", [1, 2, 3, 4, 5]),
]


def summarize(label: str, csv_path: Path):
    qubits = load_calibration(csv_path)
    n_op = sum(1 for q in qubits if q.operational)
    print(f"=== {label} ({csv_path.name}) ===")
    print(f"  {len(qubits)} qubits, {n_op} operational")
    print()
    print(f"  {'qubit':>6} {'T1 (us)':>9} {'T2 (us)':>9} {'sx err':>10} {'readout':>9} {'score':>10}")
    for q in best_qubits(qubits, 10):
        print(f"  {q.qubit:>6} {q.t1_us:>9.1f} {q.t2_us:>9.1f} {q.sx_error:>10.5f} {q.readout_error:>9.4f} {score_qubit(q):>10.2f}")

    p3, s3 = best_chain(qubits, 3)
    p5, s5 = best_chain(qubits, 5)
    print()
    print(f"  best 3-chain: {p3}, score = {s3:.2f}")
    print(f"  best 5-chain: {p5}, score = {s5:.2f}")
    return qubits, (p3, s3), (p5, s5)


def main():
    print("Two-snapshot drift review of ibm_marrakesh calibrations")
    print("=" * 70)
    print()

    qubits_25, _, _ = summarize("Apr 25 (1 day before soft_break)", CALIB_APR25)
    print()
    qubits_30, _, _ = summarize("Apr 30 (1 day after zn_mirror)", CALIB_APR30)

    print()
    print("=" * 70)
    print("Per-named-path scores across the snapshots")
    print()
    by_id_25 = {q.qubit: q for q in qubits_25}
    by_id_30 = {q.qubit: q for q in qubits_30}

    print(f"  {'path':<40} {'Apr 25':>10} {'Apr 30':>10} {'Δ':>8} {'%':>7}")
    print("  " + "-" * 78)
    for label, path in NAMED_PATHS:
        s25 = chain_score(qubits_25, path)
        s30 = chain_score(qubits_30, path)
        delta = s30 - s25
        pct = 100.0 * delta / s25 if s25 != 0 else 0.0
        print(f"  {label:<40} {s25:>10.2f} {s30:>10.2f} {delta:>+8.2f} {pct:>+6.1f}%")

    print()
    print("Per-qubit T2 drift on the path qubits (Apr-30 − Apr-25)")
    print()
    print(f"  {'qubit':>6} {'T2 Apr-25':>11} {'T2 Apr-30':>11} {'Δ T2':>8} {'%':>7} {'score Δ':>10}")
    print("  " + "-" * 60)
    for q_id in sorted({q for _, p in NAMED_PATHS for q in p}):
        if q_id not in by_id_25 or q_id not in by_id_30:
            continue
        q25 = by_id_25[q_id]
        q30 = by_id_30[q_id]
        dT2 = q30.t2_us - q25.t2_us
        pct = 100.0 * dT2 / q25.t2_us if q25.t2_us > 0 else 0.0
        ds = score_qubit(q30) - score_qubit(q25)
        print(f"  {q_id:>6} {q25.t2_us:>11.1f} {q30.t2_us:>11.1f} {dT2:>+8.1f} {pct:>+6.1f}% {ds:>+10.2f}")

    print()
    print("Reading:")
    print("  - Calibration score gap between paths (~14% [48,49,50] vs [0,1,2]) is")
    print("    the linear per-qubit-quality input.")
    print("  - F88-Lens state-level truly-baseline gap (~23× on same paths) is the")
    print("    consequence: dephasing accumulated through the circuit amplifies the")
    print("    per-qubit difference nonlinearly. Score and F88 reading move in the")
    print("    same direction but at different magnitudes; the score predicts the")
    print("    ordering, not the absolute F88 ratio.")
    print("  - Best chains [4,3,2] and [1,2,3,4,5] outscore the documented soft_break")
    print("    path [48,49,50]; running F88-Lens on those would be a clean follow-up")
    print("    to test whether the truly-baseline goes below 0.0013.")


if __name__ == "__main__":
    main()
