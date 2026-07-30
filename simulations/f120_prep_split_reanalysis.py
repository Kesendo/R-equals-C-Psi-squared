#!/usr/bin/env python3
"""F120 Kingston run, prep-conditioned re-analysis: the same-day correction.

The 2026-06-11 moment-tower run's first reading ("q13 violates pump <= Gamma
at 4-6 sigma") compared the run's pump slopes against a standard-T1 arbiter
measured 16 minutes later. This script splits the run's own counts by each
qubit's preparation bit, which the 8-basis-state mixed preparation contains
for free:

    s_b = d<Z>/dt from preparation |b>;  in any two-level Lindblad model
    s0 = -2*gamma_up <= 0  (from the ground state <Z> can only fall),
    pump = (s1 + s0)/2,    Gamma = (s1 - s0)/2,
    so pump <= Gamma  <=>  s0 <= 0,  measured IN-SITU, epoch-matched.

Result (the correction): every qubit satisfies the bound in-situ
(pump/Gamma = 0.965..0.996), the 0.4-3.5% margins, whose half read the per-qubit thermal
excitation gamma_up = -s0/2, and the apparent violation was an epoch
artifact: q13's T1 telegraphs between ~315 us (during the run) and ~430 us
(at the arbiter); q9 switches the other way (~172 us in-run vs ~75-100 us
at the arbiter); q149 looked stable across these two epochs, a reading the 10:02Z telegraph chase retired. Two-level Lindblad holds within epochs;
the epoch was the hidden variable.

BIT-ORDERING GOTCHA (for reproducers): both the preparation keys ('x000'..)
and the measured bitstrings in the saved counts are little-endian relative
to the chain order, so chain site s addresses string position (2 - s). That spot-check
below cannot settle the ordering by itself: it is chain site 1, the fixed
point of (2 - s), so both candidate rules give it identically. What settles
it is the cross-arm assert in main(): the two arms swap chain sites 1 and 2,
so the SAME physical qubit must return the same in-situ Gamma in both. It
agrees to ~5% under this mapping and to 28-59% under the reversed one. The
spot-check still pins the conditioning itself against the unconditioned
curves of the main analysis (q13 mixed <Z_1> at tau = 50:
(0.9947 + (-0.6737))/2 = 0.1605, the main table's value exactly).

Data: data/ibm_moment_tower_june2026/. Writeup: the Correction section of
experiments/F120_MOMENT_TOWER_KINGSTON.md.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

DATA = Path(__file__).resolve().parent.parent / "data" / "ibm_moment_tower_june2026"
RUN = json.load(open(DATA / "moment_tower_hardware_ibm_kingston_20260611T073908Z.json",
                     encoding="utf-8"))

CAL = {}
for r in csv.DictReader(open(DATA / "ibm_kingston_calibrations_2026-06-11T06_33_42Z.csv",
                             encoding="utf-8-sig")):
    try:
        CAL[int(r["Qubit"])] = (float(r["Prob meas1 prep0"]), float(r["Prob meas0 prep1"]))
    except (ValueError, KeyError):
        continue

LAYOUT = {a: RUN["arms"][a]["chain_to_physical"] for a in ("A", "B")}
TAUS = [int(x) for x in RUN["tau_grid_us"]]   # count keys are "tau_0", "tau_25", ...
EARLY = RUN["linear_fit_window_us"]
ARBITER_GAMMA = {149: 1 / 424.6, 13: 1 / 430.3, 9: 1 / 99.9}   # 07:55Z, 16 min later


def corrected_p1(p1_raw, qubit):
    p01, p10 = CAL[qubit]
    return (p1_raw - p01) / (1.0 - p10 - p01)


def prep_split(arm, site):
    """<Z>(tau) conditioned on this chain site's preparation bit (0 and 1)."""
    qubit = LAYOUT[arm][site]
    pos = 2 - site                              # little-endian string position
    curves = {0: [], 1: []}
    for tau in TAUS:
        counts = RUN["counts"][arm][f"tau_{tau}"]["ZZZ"]
        agg = {0: [0, 0], 1: [0, 0]}
        for prep_key, cdict in counts.items():
            pb = int(prep_key.lstrip("x")[pos])
            for bits, n in cdict.items():
                agg[pb][int(bits[pos])] += n
        for pb in (0, 1):
            p1 = agg[pb][1] / (agg[pb][0] + agg[pb][1])
            curves[pb].append(1.0 - 2.0 * corrected_p1(p1, qubit))
    return curves


def main():
    t = np.array(TAUS, dtype=float)
    m = t <= EARLY
    print("F120 prep-conditioned re-analysis (early window tau <= 75 us, readout-corrected)")
    print("=" * 96)
    worst_ratio = 0.0
    any_s0_positive = False
    gamma_by_qubit = {}
    for arm in ("A", "B"):
        print(f"Arm {arm}: chain (0,1,2) = physical {LAYOUT[arm]}")
        for site in range(3):
            qubit = LAYOUT[arm][site]
            cv = prep_split(arm, site)
            # sanity: at tau = 0 the conditioning must reproduce the preparation
            assert cv[0][0] > 0.9 and cv[1][0] < -0.9, "prep/readout orientation broken"
            c0, v0 = np.polyfit(t[m], np.array(cv[0])[m], 1, cov=True)
            c1, v1 = np.polyfit(t[m], np.array(cv[1])[m], 1, cov=True)
            s0, e0 = c0[0], np.sqrt(v0[0, 0])
            s1, e1 = c1[0], np.sqrt(v1[0, 0])
            pump, gam = (s1 + s0) / 2, (s1 - s0) / 2
            ratio = pump / gam
            worst_ratio = max(worst_ratio, ratio)
            gamma_by_qubit.setdefault(qubit, []).append(gam)
            if s0 > 2 * e0:
                any_s0_positive = True
            pth = -s0 / 2 / gam * 100
            print(f"  site {site} (q{qubit:>3}): s0 = {s0:+.2e}±{e0:.0e}  s1 = {s1:+.2e}±{e1:.0e}"
                  f"  pump = {pump:.3e}  Γ_insitu = {gam:.3e}  pump/Γ = {ratio:.3f}"
                  f"  p_th ≈ {pth:.1f}%  (arbiter Γ = {ARBITER_GAMMA[qubit]:.3e})")
    print("=" * 96)
    # The site -> bitstring-position mapping is the one thing no other check here can see:
    # the orientation assert above applies the SAME pos to prep and measurement, so it
    # tests that the two agree, not that either lands on the right chain site. The arms
    # swap chain sites 1 and 2, so a per-qubit Gamma has to survive the swap. Under the
    # reversed mapping every qubit is reported against another qubit's curve.
    for qubit, gammas in sorted(gamma_by_qubit.items()):
        spread = abs(gammas[0] - gammas[1]) / max(abs(g) for g in gammas)
        print(f"  cross-arm Γ agreement q{qubit:>3}: {spread * 100:.1f}%")
        assert spread < 0.15, (
            f"q{qubit} disagrees across arms by {spread*100:.1f}%: the chain-site to "
            f"bitstring-position mapping is wrong ({gammas})")
    assert not any_s0_positive, "a qubit shows <Z> rising from |0>: genuine two-level violation"
    assert worst_ratio < 1.0, f"in-situ pump/Gamma reached {worst_ratio:.3f} >= 1"
    print(f"IN-SITU BOUND HOLDS EVERYWHERE (worst pump/Γ = {worst_ratio:.3f}); the cross-epoch")
    print("'violation' was T1 telegraphing between the run and the arbiter (q13: ~315 vs 430 us;")
    print("q9: ~172 vs ~75-100 us; q149 steady across THESE two epochs, a reading the")
    print("10:02Z telegraph chase later retired). The protocol self-arbitrates by prep-splitting.")


if __name__ == "__main__":
    main()
