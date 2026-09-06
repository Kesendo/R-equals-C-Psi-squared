#!/usr/bin/env python3
"""Locate changes in a numerical oscillating-eigenvalue count.

A count change is only a transition bracket.  This script does not call it an
exceptional point: that requires an independent defectiveness/Jordan gate.
"""

import argparse
from pathlib import Path
import sys
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
from thermal_blackbody import GAMMA, TOL_FREQ, build_thermal_liouvillian


def oscillating_count(values, tol=TOL_FREQ):
    values = np.asarray(values, dtype=complex)
    return int(np.count_nonzero(np.abs(values.imag) > tol))


def count_at(n, n_bar):
    return oscillating_count(
        np.linalg.eigvals(build_thermal_liouvillian(n, [GAMMA] * n, n_bar, [GAMMA] * n))
    )


def transition_brackets(n, grid):
    rows = [(float(value), count_at(n, float(value))) for value in grid]
    return rows, [
        (left[0], right[0], left[1], right[1])
        for left, right in zip(rows, rows[1:]) if left[1] != right[1]
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=4)
    args = parser.parse_args()
    grid = np.unique(np.concatenate((np.linspace(0, 2, 81), [5, 10, 50])))
    print("THERMAL OSCILLATING-COUNT TRANSITIONS")
    print(f"imag tolerance={TOL_FREQ:.0e}; sigma- emission is present at nbar=0")
    for n in range(3, args.N + 1):
        rows, brackets = transition_brackets(n, grid)
        print(f"N={n}: {len(brackets)} count-change brackets on the stated grid")
        for lo, hi, before, after in brackets:
            print(f"  nbar in [{lo:.6g},{hi:.6g}]: {before}->{after}")
        print("  probes: " + ", ".join(f"{x:g}:{count}" for x, count in rows if x in (0, 1, 2, 5, 10, 50)))
    print("Verdict: these are numerical count transitions, not EP certificates.")
    print("A defectiveness/Jordan-chain test is required before using the EP label.")


if __name__ == "__main__":
    main()
