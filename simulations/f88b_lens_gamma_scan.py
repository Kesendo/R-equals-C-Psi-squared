"""Parameter scan: how do (γ_Z, γ_T1) affect F88b-Lens Π²-odd/mem for the three
F87 categories under |+−+⟩ at t=0.8 with h_y=0.05?

Motivation: Trotter-asymmetry hypothesis was falsified — Trotter pushes Π²-odd
UP, not down. Hard's hardware reading (0.276) sits BELOW continuous-Lindbladian
ideal (0.381) with the calibration-derived γ_Z=0.1, γ_T1=0.046 baseline.

Two new candidate explanations from calibration analysis:
(a) T1/T2 drift: experiment was 2026-04-26, calibration was 2026-04-25; over 10
    days the T1 values drift ±15% and q0's T2 by 30%. γ_T1 range 0.035-0.047.
(b) Hahn-echo vs T2*: IBM's "T2" column is T2_echo. The gate-time effective
    dephasing during Trotter execution corresponds to T2* (Ramsey), not T2_echo.
    Per IBM_ABSORPTION_THEOREM.md: γ*/γ_echo ≈ 6.2 on Heron-class hardware.
    So the effective γ_Z during gates might be up to ~6× our framework value.

This scan sweeps γ_Z ∈ {0.1, 0.2, 0.3, 0.5, 0.7} (calibration → ~7× boost) and
γ_T1 ∈ {0.046, 0.1, 0.2} for each of truly/soft/hard, reading Π²-odd/mem from
the F88b-Lens. Goal: identify which (γ_Z, γ_T1) region brings hard down to the
hardware value 0.276 — and check whether truly and soft remain consistent with
their hardware values (0.030, 0.744) in the same region.

If a single (γ_Z, γ_T1) point reproduces all three categories simultaneously,
then the hardware-vs-calibration mismatch is consistent and quantifiable.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from f88b_lens_ibm_framework_snapshots import (  # noqa: E402
    f88b_lens_2qubit,
    lindbladian_plus_minus_plus_rho_q0q2,
)

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


HARD_H_TERMS = [("X", "X", 1.0), ("X", "Y", 1.0)]
SOFT_H_TERMS = [("X", "Y", 1.0), ("Y", "X", 1.0)]
TRULY_H_TERMS = [("X", "X", 1.0), ("Y", "Y", 1.0)]
H_Y = 0.05
T = 0.8

HW_VALUES = {"truly": 0.0297, "soft": 0.7444, "hard": 0.2763}

CATEGORIES = [("truly", TRULY_H_TERMS), ("soft", SOFT_H_TERMS), ("hard", HARD_H_TERMS)]

GAMMA_Z_SWEEP = [0.1, 0.2, 0.3, 0.5, 0.7]
GAMMA_T1_SWEEP = [0.046, 0.1, 0.2]


def lens_value(H_terms, gamma_z, gamma_t1):
    rho = lindbladian_plus_minus_plus_rho_q0q2(
        t=T, H_terms=H_terms, gamma=gamma_z, h_y=H_Y, gamma_t1=gamma_t1,
    )
    return f88b_lens_2qubit(rho)["pi2_odd_in_memory"]


def main() -> None:
    print("F88b-Lens (γ_Z, γ_T1) parameter scan")
    print("=" * 78)
    print(f"Setup: N=3, |+−+⟩, t={T}, h_y={H_Y} (Marrakesh 2026-04-29 detected leak)")
    print(f"Hardware targets (Marrakesh 2026-04-26):")
    print(f"  truly={HW_VALUES['truly']:.4f}  soft={HW_VALUES['soft']:.4f}  hard={HW_VALUES['hard']:.4f}")
    print()
    print(f"γ_Z range  {GAMMA_Z_SWEEP}  (calibration baseline = 0.1)")
    print(f"γ_T1 range {GAMMA_T1_SWEEP}  (calibration baseline = 0.046)")
    print()

    for cat_label, H_terms in CATEGORIES:
        print(f"=== {cat_label} === (hardware {HW_VALUES[cat_label]:.4f})")
        # Header
        print(f"  γ_T1 \\ γ_Z " + "".join(f"  {gz:>6.2f}" for gz in GAMMA_Z_SWEEP))
        print("  " + "-" * (12 + 8 * len(GAMMA_Z_SWEEP)))
        for gt1 in GAMMA_T1_SWEEP:
            row = []
            for gz in GAMMA_Z_SWEEP:
                v = lens_value(H_terms, gz, gt1)
                # Mark the cell with * if within 0.02 of hardware
                marker = "*" if abs(v - HW_VALUES[cat_label]) < 0.02 else " "
                row.append(f"{v:>6.4f}{marker}")
            print(f"  {gt1:>10.3f}  " + "  ".join(row))
        print()

    # Find the best matching (γ_Z, γ_T1) — minimize total squared deviation
    # from all three hardware values simultaneously.
    print("=== Global fit: which (γ_Z, γ_T1) matches ALL three best? ===")
    best = None
    for gz in GAMMA_Z_SWEEP:
        for gt1 in GAMMA_T1_SWEEP:
            sse = 0.0
            vals = {}
            for cat_label, H_terms in CATEGORIES:
                v = lens_value(H_terms, gz, gt1)
                vals[cat_label] = v
                sse += (v - HW_VALUES[cat_label]) ** 2
            if best is None or sse < best[0]:
                best = (sse, gz, gt1, vals)
    sse_best, gz_best, gt1_best, vals_best = best
    print(f"  Best: γ_Z={gz_best}, γ_T1={gt1_best}, total SSE={sse_best:.4f}")
    for cat in ["truly", "soft", "hard"]:
        v = vals_best[cat]
        hw = HW_VALUES[cat]
        print(f"    {cat:>8s}: ideal {v:.4f}  hw {hw:.4f}  Δ {v - hw:+.4f}")


if __name__ == "__main__":
    main()
