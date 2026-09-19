"""Marrakesh r = T2/(2T1) proxy-band review across two April 2026 snapshots.
applying the BOTH_SIDES_VISIBLE / ibm_torino_history framework to the 5-day
window that brackets the soft_break (2026-04-26) and zn_mirror (2026-04-29)
runs.

The 6-month Torino history (`docs/BOTH_SIDES_VISIBLE.md`) found 16 of 133
Eagle-r3 qubits changing sides of a model-derived R* line on the daily
calibration cadence. The same finite tally is made for Marrakesh over five days:

- Which qubits sit below R* on each calibration
- Which qubits change proxy bands between Apr-25 and Apr-30
- Where the F88b-Lens experimental paths sit on this map

Reading: a band change records T1/T2 calibration drift relative to a
free-single-transmon `|+>` normalized-purity proxy. It is a scheduling covariate,
not a physical phase identity or a property of the later Hamiltonian run.

Equality belongs to the at-or-above band. The script reports raw-T1/T2 band
membership and changes only; it makes no quantum/classical or causal claim.
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

R_STAR = 0.21275477982200533

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


def is_below_rstar(q):
    return r_param(q) < R_STAR


def main():
    print("Marrakesh normalized-purity proxy bands across a 5-day window")
    print(f"R* = {R_STAR:.17g}: labels are below-R* or at-or-above-R*")
    print("=" * 78)
    print()

    q25 = {q.qubit: q for q in load_calibration(CALIB_APR25)}
    q30 = {q.qubit: q for q in load_calibration(CALIB_APR30)}

    below_25 = {qid for qid, q in q25.items() if is_below_rstar(q)}
    below_30 = {qid for qid, q in q30.items() if is_below_rstar(q)}
    changed_to_below = below_30 - below_25
    changed_to_at_or_above = below_25 - below_30
    stable_below = below_25 & below_30
    stable_at_or_above = (set(q25) - below_25) & (set(q30) - below_30)

    print(f"Apr-25: {len(below_25):3d} / {len(q25)} rows below-R*")
    print(f"Apr-30: {len(below_30):3d} / {len(q30)} rows below-R*")
    print()
    print(f"Below-R* on both dates:         {len(stable_below):3d}")
    print(f"At-or-above-R* on both dates:  {len(stable_at_or_above):3d}")
    print(f"Changed to below-R*:           {len(changed_to_below):3d}")
    print(f"Changed to at-or-above-R*:     {len(changed_to_at_or_above):3d}")
    print()

    # Show flippers in detail
    flippers = sorted(changed_to_below | changed_to_at_or_above)
    if flippers:
        print(f"Qubits that flipped sides in the 5-day window ({len(flippers)} total):")
        print(f"  {'qubit':>6} {'r Apr-25':>10} {'r Apr-30':>10} {'Δr':>9} {'side25':>8} {'side30':>8}")
        print("  " + "-" * 60)
        for qid in flippers:
            r25 = r_param(q25[qid])
            r30 = r_param(q30[qid])
            s25 = "below-R*" if is_below_rstar(q25[qid]) else "at/above"
            s30 = "below-R*" if is_below_rstar(q30[qid]) else "at/above"
            print(f"  {qid:>6} {r25:>10.4f} {r30:>10.4f} {r30 - r25:>+9.4f} {s25:>8} {s30:>8}")
        print()

    # Closest-to-boundary qubits on each calibration (oscillator candidates extending
    # the 16/133 Torino result; these are the ones that could flip with another day's drift)
    by_dist_25 = sorted(q25.values(), key=lambda q: abs(r_param(q) - R_STAR))
    by_dist_30 = sorted(q30.values(), key=lambda q: abs(r_param(q) - R_STAR))
    print(f"Top-10 closest-to-boundary on Apr-25 (potential oscillators):")
    print(f"  {'qubit':>6} {'r':>8} {'|Δ|':>8} {'side':>8}")
    for q in by_dist_25[:10]:
        r = r_param(q)
        print(f"  {q.qubit:>6} {r:>8.4f} {abs(r - R_STAR):>8.4f} {'below' if r < R_STAR else 'at/above':>8}")
    print()
    print(f"Top-10 closest-to-boundary on Apr-30:")
    print(f"  {'qubit':>6} {'r':>8} {'|Δ|':>8} {'side':>8}")
    for q in by_dist_30[:10]:
        r = r_param(q)
        print(f"  {q.qubit:>6} {r:>8.4f} {abs(r - R_STAR):>8.4f} {'below' if r < R_STAR else 'at/above':>8}")
    print()

    # Where do experimental paths sit?
    print("Experimental paths through this lens:")
    print(f"  {'path':<42} {'side counts':>15}")
    print("  " + "-" * 70)
    for label, path in NAMED_PATHS:
        a25 = sum(1 for qid in path if qid in below_25)
        a30 = sum(1 for qid in path if qid in below_30)
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
    print(f"  - The rows form below-R* and at-or-above-R* proxy bands.")
    print("  - A band change is a calibration-stability warning for later scheduling;")
    print("    it neither classifies the hardware state nor identifies a drift mechanism.")
    print("  - Torino and Marrakesh counts use different backends and windows, so their")
    print("    comparison is a confounded association, not a mechanism check.")


if __name__ == "__main__":
    main()
