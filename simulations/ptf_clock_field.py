"""ptf_clock_field.py: the real clock field (four IBM machines, three pre-registered questions).

"Let us find PTF in real" (Tom, 2026-07-18) -> the clock reading: the repo already
OWNS a measured per-site clock field: data/ibm_history/ (per-qubit daily T1/T2/
r_param over months, four machines). This script reads that field in clock
currency and asks three PRE-REGISTERED questions. Vocabulary note: mapping this
onto PTF is a READING (the in-model strong closure-law claim was falsified,
EQ-014); what is tested here is structure of the real data itself.

Findings are written up in experiments/CLOCK_FIELD_SITE_OWNED.md; the committed
output copy lives in simulations/results/clock_field/ptf_clock_field_out.txt.

Tick rate used: rate = 1/T2 (a genuine additive rate; the dephasing-dominated
clock). Control: 1/T1.

Q1 (is the clock field site-owned?): variance decomposition of ln(rate):
    between-qubit vs within-qubit-across-days; ICC = Vb/(Vb+Vw).
    PRE-REGISTERED verdict: ICC > 0.5 -> the field is a site property.
Q2 (total fixed, shares rotate? THE WEAVE SIGNATURE, labeled reading):
    on the coverage-complete matrix R[q,d]: CV of the chip total T(d)=sum_q R
    vs the independence null (each qubit's day-series independently permuted,
    500 perms). PRE-REGISTERED verdicts:
      CV_total < 2.5th pct of null  -> COMPENSATION (shares rotate, total held)
      CV_total > 97.5th pct of null -> COMMON-MODE chip drift
      else                          -> independent drift (no total conservation)
Q3 (the quarter-boundary field): fraction of records with crosses_quarter,
    and #qubits that sit on both sides of the 1/4 boundary across days
    (extends docs/BOTH_SIDES_VISIBLE.md).
Also reported: day-to-day telegraphing |delta ln rate| (median, p95).

Caveats (stated up front): calibration values carry per-measurement fit error
(inflates within-qubit variance -> biases Q1 ICC DOWN, biases Q2 toward
"independent"); chip-wide temperature/TLS drift is a physical common-mode.
"""

import csv
import numpy as np
from pathlib import Path

RNG = np.random.default_rng(135)
BASE = Path(__file__).resolve().parent.parent / "data" / "ibm_history"
FILES = {
    "torino": BASE / "ibm_torino_history.csv",
    "marrakesh": BASE / "results" / "ibm_marrakesh_history.csv",
    "kingston": BASE / "results" / "ibm_kingston_history.csv",
    "fez": BASE / "results" / "ibm_fez_history.csv",
}
N_PERM = 500


def load(fp):
    rows = []
    with open(fp, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                t1 = float(r["T1_us"])
                t2 = float(r["T2_us"])
            except (ValueError, KeyError):
                continue
            if t1 <= 0 or t2 <= 0:
                continue
            rows.append((r["date"], int(r["qubit"]), t1, t2,
                         r.get("crosses_quarter", ""), r.get("distance_from_quarter", "")))
    return rows


def coverage_matrix(rows, which):
    """R[q, d] of rate = 1/T2 (or 1/T1), qubits with >=90% date coverage,
    dates where all selected qubits are present."""
    dates = sorted({r[0] for r in rows})
    qubits = sorted({r[1] for r in rows})
    val = {}
    for (d, q, t1, t2, _, _) in rows:
        val[(q, d)] = 1.0 / (t2 if which == "T2" else t1)
    keep_q = [q for q in qubits if sum((q, d) in val for d in dates) >= 0.9 * len(dates)]
    keep_d = [d for d in dates if all((q, d) in val for q in keep_q)]
    R = np.array([[val[(q, d)] for d in keep_d] for q in keep_q])
    return R, keep_q, keep_d


def q1_icc(R):
    ln = np.log(R)
    grand = ln.mean()
    vb = ((ln.mean(axis=1) - grand) ** 2).mean()
    vw = ln.var(axis=1).mean()
    return vb / (vb + vw), vb, vw


def q2_total_vs_shares(R):
    total = R.sum(axis=0)
    cv_total = total.std() / total.mean()
    null = []
    for _ in range(N_PERM):
        P = np.array([RNG.permutation(row) for row in R])
        t = P.sum(axis=0)
        null.append(t.std() / t.mean())
    null = np.array(null)
    pct = float((null < cv_total).mean() * 100)
    return cv_total, null, pct


def main():
    print("PTF clock field: the real per-site clock (data/ibm_history)")
    print(f"rate = 1/T2 unless noted; permutations N={N_PERM}, seed 135")
    for name, fp in FILES.items():
        if not fp.exists():
            print(f"\n== {name}: FILE MISSING {fp}")
            continue
        rows = load(fp)
        R, qs, ds = coverage_matrix(rows, "T2")
        print(f"\n== {name}: {len(rows)} records; coverage matrix {R.shape[0]} qubits x {R.shape[1]} days")
        if R.shape[1] < 10:
            print("   too few complete days, skipping stats")
            continue

        icc, vb, vw = q1_icc(R)
        verdict1 = "SITE-OWNED (ICC > 0.5)" if icc > 0.5 else "day-noise dominated"
        print(f"  Q1 ln-rate variance: between-qubit {vb:.4f}  within-qubit {vw:.4f}  ICC={icc:.3f}  -> {verdict1}")

        cv, null, pct = q2_total_vs_shares(R)
        lo, hi = np.percentile(null, 2.5), np.percentile(null, 97.5)
        if cv < lo:
            verdict2 = "COMPENSATION: total held, shares rotate (the weave signature)"
        elif cv > hi:
            verdict2 = "COMMON-MODE chip drift (total varies MORE than independent)"
        else:
            verdict2 = "independent drift (no total conservation)"
        print(f"  Q2 CV(chip total)={cv:.4f}  independence null [{lo:.4f}, {hi:.4f}] (pct rank {pct:.0f}%)  -> {verdict2}")

        # control with 1/T1
        R1, _, _ = coverage_matrix(rows, "T1")
        cv1, null1, _ = q2_total_vs_shares(R1)
        lo1, hi1 = np.percentile(null1, 2.5), np.percentile(null1, 97.5)
        tag1 = "compensation" if cv1 < lo1 else ("common-mode" if cv1 > hi1 else "independent")
        print(f"  Q2 control 1/T1: CV={cv1:.4f} null [{lo1:.4f}, {hi1:.4f}] -> {tag1}")

        crossers = {q for (d, q, t1, t2, c, _) in rows if str(c).strip().lower() == "true"}
        sides = {}
        for (d, q, t1, t2, _, dist) in rows:
            try:
                s = float(dist)
            except ValueError:
                continue
            sides.setdefault(q, set()).add(s > 0)
        both = sum(1 for v in sides.values() if len(v) == 2)
        print(f"  Q3 quarter boundary: {len(crossers)} qubits ever cross in-window; {both} qubits sit on BOTH sides across days")

        dln = np.abs(np.diff(np.log(R), axis=1)).ravel()
        print(f"  telegraph |d ln rate|/day: median {np.median(dln):.3f}  p95 {np.percentile(dln, 95):.3f}")


if __name__ == "__main__":
    main()
