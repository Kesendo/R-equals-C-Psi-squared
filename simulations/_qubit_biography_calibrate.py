"""Calibrate the biography classifier against the named qubits in BOTH_SIDES_VISIBLE.md:
Q98 (lifecycle), Q72 (rhythmic / twitch), Q105 (long active then silent),
Q80 (consistent crossing), Q70 (mostly silent). Print their raw r-trajectory
statistics so we can pick thresholds that match the doc's labels.
"""
from __future__ import annotations

import sys
from pathlib import Path

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _qubit_biography import load_history, R_STAR, DEFAULT_HISTORY, archetype_from_series


def main():
    by_qubit = load_history(DEFAULT_HISTORY)

    named = {
        72: "rhythmic / twitch (66.9% crossing)",
        80: "consistent crosser (1.9% deviation)",
        98: "lifecycle: tune → pulse → fade",
        105: "long active then silent (56.9%)",
        70: "mostly silent, brief pulses (26.0%)",
        68: "mostly silent (23.2%)",
    }

    print(f"R* = {R_STAR}")
    print()
    print(f"  {'qubit':>5} {'days':>5} {'r mean':>8} {'r std':>8} {'crossing %':>11} "
          f"{'walk':>6} {'archetype':>16} {'doc label'}")
    print("  " + "-" * 110)
    for qid, label in named.items():
        if qid not in by_qubit:
            continue
        recs = by_qubit[qid]
        rs = np.array([rec[3] for rec in recs])
        crossings = (rs < R_STAR).mean() * 100
        signs = np.sign(rs - R_STAR)
        flips = int(np.sum(np.diff(signs) != 0))
        walk = flips / max(len(rs) - 1, 1)
        arch = archetype_from_series(rs)
        print(f"  Q{qid:<4} {len(rs):>5} {rs.mean():>8.4f} {rs.std():>8.4f} {crossings:>10.1f}% "
              f"{walk:>6.3f} {arch:>16} {label}")
    print()
    print("Interpretation:")
    print(f"  - 'walk' = fraction of days where r flipped to the other side (boundary crossings).")
    print(f"  - 'crossing %' = total fraction of days below R*.")
    print(f"  - High walk + crossing ≈ 50% → twitch / rhythmic (boundary-living)")
    print(f"  - High walk + multiple long runs of pulse vs silent → lifecycle")
    print(f"  - Low walk + high crossing → pulse-stable")
    print(f"  - Low walk + low crossing → silent-stable")


if __name__ == "__main__":
    main()
