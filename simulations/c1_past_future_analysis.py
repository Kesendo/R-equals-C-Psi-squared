#!/usr/bin/env python3
"""c1_past_future_analysis.py

Recovery script: the c1_past_future_test.py run completed all 3 heavy
eigendecomps and all 7 per-state fits (6 hours of CPU), then crashed
on a format-string bug in the Π-pair table header. The c_1 values
were printed to the log but not saved to JSON before the crash.

This script hardcodes the printed values, redoes the Π-pair analysis
(which is trivial scalar arithmetic), and writes the JSON so the data
is preserved.

Source log: simulations/results/c1_past_future_test_fullrun.log
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np


N = 7
GAMMA_0 = 0.05
J = 1.0
DJ = 0.01
V_N = 1.0 + np.cos(np.pi / N)

# From simulations/results/c1_past_future_test_fullrun.log, per-state output:
# The log shows E_k, role, c_1, and 7 alpha_i values per state at dJ=+0.01.
raw_data = [
    # (k, E_k, role, c_1, alpha_i list at dJ=+0.01)
    (1, +1.8478, "symmetric", +0.9700, [1.008, 1.014, 1.006, 1.000, 0.990, 0.992, 1.000]),
    (2, +1.4142, "antisymmetric", +0.0371, [1.011, 1.002, 0.983, 1.000, 1.010, 0.998, 0.996]),
    (3, +0.7654, "symmetric", +0.3574, [0.992, 0.988, 1.007, 1.006, 0.998, 1.010, 1.004]),
    (4, +0.0000, "antisymmetric (Pi-self-partner, E=0 standing wave)",
                  +2.1357, [0.989, 1.000, 1.008, 1.001, 1.010, 1.001, 1.013]),
    (5, -0.7654, "symmetric", +0.3574, [0.992, 0.988, 1.007, 1.006, 0.998, 1.010, 1.004]),
    (6, -1.4142, "antisymmetric", +0.0371, [1.011, 1.002, 0.983, 1.000, 1.010, 0.998, 0.996]),
    (7, -1.8478, "symmetric", +0.9700, [1.008, 1.014, 1.006, 1.000, 0.990, 0.992, 1.000]),
]


def main():
    print("=" * 70)
    print(f"c_1 past/future analysis (recovery) at N = {N}")
    print("=" * 70)
    print(f"  gamma_0 = {GAMMA_0}, J = {J}, dJ = +/- {DJ}, V(N) = {V_N:.4f}")
    print(f"  Prediction for psi_1+vac: c_1 ~ 0.5 * V(N) = {0.5 * V_N:.4f}")
    print(f"  EQ-014 reports c_1(psi_1) = 0.97 at N=7 (parallel session)")

    results_by_state = {}
    print(f"\n  per-state c_1 values:")
    print(f"  {'k':>2} {'E_k':>10} {'role':>50} {'c_1':>10}")
    for k, Ek, role, c1, _alpha in raw_data:
        print(f"  {k:>2d} {Ek:>+10.4f} {role:>50} {c1:>+10.4f}")
        results_by_state[f"psi_{k}"] = {
            "k": k, "E_k": Ek, "role": role, "c_1": c1,
            "alpha_plus": _alpha,
        }

    # Pi-pair analysis
    print(f"\n  Pi-pair analysis: c_1(psi_k) vs c_1(psi_{{N+1-k}})")
    print(f"  {'k':>2} {'mirror':>8} {'c_1(k)':>10} {'c_1(N+1-k)':>12} "
          f"{'sum':>10} {'diff':>10} {'identical?':>12}")
    pi_pair_data = {}
    for (k, _, _, c_k, _), (km, _, _, c_m, _) in zip(raw_data, reversed(raw_data)):
        if k > km:
            continue
        s = c_k + c_m
        d = c_k - c_m
        identical = abs(d) < 1e-4
        note = ""
        if k == km:
            note = "  (self-mirror)"
        print(f"  {k:>2d} {km:>8d} {c_k:>+10.4f} {c_m:>+12.4f} "
              f"{s:>+10.4f} {d:>+10.4f} {str(identical):>12}{note}")
        pi_pair_data[f"pair_{k}_{km}"] = {
            "c_1_k": c_k, "c_1_mirror": c_m, "sum": s, "diff": d,
            "identical_to_1e-4": identical,
        }

    # Symmetric vs antisymmetric summary
    print(f"\n  Summary by reflection behavior under i <-> N-1-i:")
    symmetric = [r for r in raw_data if r[2].startswith("symmetric")]
    antisym = [r for r in raw_data if r[2].startswith("antisym")]
    print(f"  Symmetric states (ψ_1, ψ_3, ψ_5, ψ_7):")
    for k, Ek, _, c1, _ in symmetric:
        print(f"    psi_{k}: c_1 = {c1:+.4f},  c_1/V(N) = {c1/V_N:+.4f}")
    print(f"  Antisymmetric states (ψ_2, ψ_4, ψ_6):")
    for k, Ek, _, c1, _ in antisym:
        note = "  [self-Pi-partner, NOT suppressed]" if k == 4 else ""
        print(f"    psi_{k}: c_1 = {c1:+.4f},  c_1/V(N) = {c1/V_N:+.4f}{note}")

    # Verdict
    print(f"\n  Key findings:")
    c4 = next(c for k, _, _, c, _ in raw_data if k == 4)
    c1_val = next(c for k, _, _, c, _ in raw_data if k == 1)
    print(f"    1. Pi-pair c_1 values are IDENTICAL to 4 decimal places:")
    print(f"       c_1(psi_1) = c_1(psi_7) = +0.9700")
    print(f"       c_1(psi_2) = c_1(psi_6) = +0.0371")
    print(f"       c_1(psi_3) = c_1(psi_5) = +0.3574")
    print(f"    2. psi_4 (E=0 standing wave, Pi-self-partner) has c_1 = {c4:+.4f},")
    print(f"       the LARGEST value, more than 2x c_1(psi_1) = {c1_val:+.4f}.")
    print(f"       Not closure-protected; maximally closure-breaking.")
    print(f"    3. Antisymmetric non-self-partner modes (psi_2, psi_6) are")
    print(f"       suppressed; antisymmetric self-partner (psi_4) is amplified.")
    print(f"       The suppression rule needs the Pi-partner qualifier.")
    print(f"    4. Tom's parallel session reported c_1(psi_4) ~ -0.12 via")
    print(f"       first-order slow-mode mixing; direct RK4 gives +2.14. Likely")
    print(f"       a degenerate-PT artifact (|vac> and psi_4 both at E=0).")

    # Save
    results_dir = Path(__file__).parent / "results" / "c1_past_future_test"
    results_dir.mkdir(parents=True, exist_ok=True)
    out = {
        "N": N, "gamma_0": GAMMA_0, "J": J, "dJ": DJ, "V_N": V_N,
        "defect_bond": [0, 1],
        "initial_state_family": "(|vac> + |psi_k>) / sqrt(2) for k = 1..7",
        "results_by_state": results_by_state,
        "pi_pair_analysis": pi_pair_data,
        "prediction_psi_1": 0.5 * V_N,
        "eq014_psi_1_match": 0.97,
        "source_log": "simulations/results/c1_past_future_test_fullrun.log",
        "note": ("c_1 values recovered from console log after the main "
                 "script crashed on a format-string bug in the Pi-pair "
                 "table header. All heavy computation (3 dense "
                 "eigendecomps + 7 state propagations at N=7) completed "
                 "successfully."),
    }
    path = results_dir / "past_future_test.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nSaved: {path}")


if __name__ == "__main__":
    main()
