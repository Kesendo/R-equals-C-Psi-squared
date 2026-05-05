"""Marrakesh r = T2/(2T1) crossing review across the two April 2026 snapshots,
applying the BOTH_SIDES_VISIBLE / ibm_torino_history framework to the 5-day
window that brackets the soft_break (2026-04-26) and zn_mirror (2026-04-29)
runs.

The 6-month Torino history (`docs/BOTH_SIDES_VISIBLE.md`) found 16 of 133 Eagle-r3
qubits oscillating around r = 0.213 (the CΨ = ¼ boundary) on the daily-calibration
cadence. The same plot for Marrakesh on its 5-day window (Heron r2, 156 qubits):

- Which qubits sit BELOW r = 0.213 on each calibration (CΨ < ¼ regime: pure-quantum
  side decays faster than the palindromic mirror, T2 << 2·T1)
- Which qubits FLIP across the boundary between Apr-25 and Apr-30 (a 5-day-window
  oscillator candidate, comparable to Torino's Q98/Q72)
- Where the F88-Lens experimental paths sit on this map

Reading: a qubit that flips between Apr-25 and Apr-30 is one whose CΨ regime
identity is unstable on F88-Lens-experiment timescale. Picking such a qubit means
the r-regime in which the Hamiltonian is run is not the same regime the
calibration suggested. The C# `IbmCalibration` score does not currently include
the r-boundary signal; this review tests whether it should.

Note on framing: the March 2026 update to BOTH_SIDES_VISIBLE clarified that the
X/. complement reading is logically a bit-flip (tautological). The boundary
crossing itself is not tautological: it is observed daily T2 drift around a
specific physical threshold (r = 0.213, where CΨ touches ¼). This script
reports the crossing structure only; no complement claim is made.
"""
from __future__ import annotations

import sys
from pathlib import Path

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ibm_calibration import load_calibration, score_qubit

REPO_ROOT = Path(__file__).resolve().parents[1]
CALIB_DIR = REPO_ROOT / "data" / "ibm_calibration_snapshots"
CALIB_APR25 = CALIB_DIR / "ibm_marrakesh_calibrations_2026-04-25T11_28_00Z.csv"
CALIB_APR30 = CALIB_DIR / "ibm_marrakesh_calibrations_2026-04-30T16_25_19Z.csv"

R_THRESHOLD = 0.213  # CΨ = ¼ boundary; r = T2/(2T1) below means CΨ_min crosses ¼

NAMED_PATHS = [
    ("framework_snapshots [0, 1, 2]", [0, 1, 2]),
    ("soft_break / zn_mirror [48, 49, 50]", [48, 49, 50]),
    ("Apr-25 best 3-chain [4, 3, 2]", [4, 3, 2]),
    ("Apr-25 best 5-chain [1, 2, 3, 4, 5]", [1, 2, 3, 4, 5]),
    ("Apr-30 best 3-chain [93, 94, 95]", [93, 94, 95]),
    ("Apr-30 best 5-chain [95, 94, 93, 92, 91]", [95, 94, 93, 92, 91]),
]


def r_param(q):
    return q.t2_us / (2.0 * q.t1_us) if q.t1_us > 0 else 0.0


def crosses(q):
    return r_param(q) < R_THRESHOLD


def main():
    print("Marrakesh r = T2 / (2·T1) boundary review across 5-day window")
    print(f"Boundary: r = {R_THRESHOLD} (CΨ_min = ¼; below = quantum-side regime)")
    print("=" * 78)
    print()

    q25 = {q.qubit: q for q in load_calibration(CALIB_APR25)}
    q30 = {q.qubit: q for q in load_calibration(CALIB_APR30)}

    # Boundary-crossing tally
    cross_25 = {qid for qid, q in q25.items() if crosses(q)}
    cross_30 = {qid for qid, q in q30.items() if crosses(q)}
    flip_to_quantum = cross_30 - cross_25  # was above → went below
    flip_to_classical = cross_25 - cross_30  # was below → went above
    stable_quantum = cross_25 & cross_30
    stable_classical = (set(q25) - cross_25) & (set(q30) - cross_30)

    print(f"Apr-25: {len(cross_25):3d} / {len(q25)} qubits below r={R_THRESHOLD} (quantum-side regime)")
    print(f"Apr-30: {len(cross_30):3d} / {len(q30)} qubits below r={R_THRESHOLD}")
    print()
    print(f"Stable quantum side (both):    {len(stable_quantum):3d}")
    print(f"Stable classical side (both):  {len(stable_classical):3d}")
    print(f"Flipped → quantum (Apr-30):    {len(flip_to_quantum):3d}  ← oscillator candidates")
    print(f"Flipped → classical (Apr-30):  {len(flip_to_classical):3d}  ← oscillator candidates")
    print()

    # Show flippers in detail
    flippers = sorted(flip_to_quantum | flip_to_classical)
    if flippers:
        print(f"Qubits that flipped sides in the 5-day window ({len(flippers)} total):")
        print(f"  {'qubit':>6} {'r Apr-25':>10} {'r Apr-30':>10} {'Δr':>9} {'side25':>8} {'side30':>8}")
        print("  " + "-" * 60)
        for qid in flippers:
            r25 = r_param(q25[qid])
            r30 = r_param(q30[qid])
            s25 = "quantum" if crosses(q25[qid]) else "classic"
            s30 = "quantum" if crosses(q30[qid]) else "classic"
            print(f"  {qid:>6} {r25:>10.4f} {r30:>10.4f} {r30 - r25:>+9.4f} {s25:>8} {s30:>8}")
        print()

    # Closest-to-boundary qubits on each calibration (oscillator candidates extending
    # the 16/133 Torino result; these are the ones that could flip with another day's drift)
    by_dist_25 = sorted(q25.values(), key=lambda q: abs(r_param(q) - R_THRESHOLD))
    by_dist_30 = sorted(q30.values(), key=lambda q: abs(r_param(q) - R_THRESHOLD))
    print(f"Top-10 closest-to-boundary on Apr-25 (potential oscillators):")
    print(f"  {'qubit':>6} {'r':>8} {'|Δ|':>8} {'side':>8}")
    for q in by_dist_25[:10]:
        r = r_param(q)
        print(f"  {q.qubit:>6} {r:>8.4f} {abs(r - R_THRESHOLD):>8.4f} {'quantum' if r < R_THRESHOLD else 'classic':>8}")
    print()
    print(f"Top-10 closest-to-boundary on Apr-30:")
    print(f"  {'qubit':>6} {'r':>8} {'|Δ|':>8} {'side':>8}")
    for q in by_dist_30[:10]:
        r = r_param(q)
        print(f"  {q.qubit:>6} {r:>8.4f} {abs(r - R_THRESHOLD):>8.4f} {'quantum' if r < R_THRESHOLD else 'classic':>8}")
    print()

    # Where do experimental paths sit?
    print("Experimental paths through this lens:")
    print(f"  {'path':<42} {'side counts':>15}")
    print("  " + "-" * 70)
    for label, path in NAMED_PATHS:
        a25 = sum(1 for qid in path if qid in cross_25)
        a30 = sum(1 for qid in path if qid in cross_30)
        flips = sum(1 for qid in path if qid in flippers)
        n = len(path)
        print(f"  {label:<42} Apr-25 {a25}/{n} q  Apr-30 {a30}/{n} q  flips={flips}")
    print()
    print("Per-qubit r values on the experimental paths:")
    path_qubits = sorted({qid for _, p in NAMED_PATHS for qid in p})
    print(f"  {'qubit':>6} {'r Apr-25':>10} {'r Apr-30':>10} {'Δr':>9} {'flipped?':>10}")
    for qid in path_qubits:
        if qid in q25 and qid in q30:
            r25 = r_param(q25[qid])
            r30 = r_param(q30[qid])
            flag = "FLIP" if qid in flippers else ""
            print(f"  {qid:>6} {r25:>10.4f} {r30:>10.4f} {r30 - r25:>+9.4f} {flag:>10}")

    print()
    print("Reading:")
    print(f"  - Marrakesh has a meaningful sub-population of qubits in the CΨ-quantum")
    print(f"    regime (r < {R_THRESHOLD}). Daily-T2 drift can flip individual qubits across")
    print(f"    the boundary on the same timescale as our QPU experiments.")
    print(f"  - Path qubits sitting near the boundary are calibration-fragile for F88-Lens:")
    print(f"    the r-regime at submission time may not match the regime at execution.")
    print(f"  - The Torino 6-month-history finding (16/133 ≈ 12% oscillators) gives the")
    print(f"    long-window prior; the Marrakesh 5-day-window count is the consistency")
    print(f"    check that the same drift mechanism is active here.")


if __name__ == "__main__":
    main()
