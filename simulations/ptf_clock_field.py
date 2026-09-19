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
Q3 (free-|+> normalized-purity proxy): carry the archive's raw-derived flag
    through a visibly marked legacy adapter and count qubits observed in both
    proxy bands, below-R* and at-or-above-R*.  Rounded display columns are not
    reclassified, and this is not a physical phase assignment.
Also reported: day-to-day telegraphing |delta ln rate| (median, p95).

Caveats (stated up front): calibration values carry per-measurement fit error
(inflates within-qubit variance -> biases Q1 ICC DOWN, biases Q2 toward
"independent"); chip-wide temperature/TLS drift is a physical common-mode.
"""

import argparse
import csv
import os
import numpy as np
from pathlib import Path

SCRIPT_PATH = Path(os.path.abspath(__file__))
BASE = SCRIPT_PATH.parent.parent / "data" / "ibm_history"
FILES = {
    "torino": BASE / "ibm_torino_history.csv",
    "marrakesh": BASE / "results" / "ibm_marrakesh_history.csv",
    "kingston": BASE / "results" / "ibm_kingston_history.csv",
    "fez": BASE / "results" / "ibm_fez_history.csv",
}
N_PERM = 500
RNG_SEED = 135
R_STAR = 0.21275477982200533
DEFAULT_OUTPUT_PATH = (
    SCRIPT_PATH.parent / "results" / "clock_field" /
    "ptf_clock_field_out.txt"
)


def _adapt_archived_proxy_band(row):
    """LEGACY-CALIBRATION-SCHEMA: carry the archive's raw-derived flag."""
    value = str(row.get("crosses_quarter", "")).strip().lower()
    if value not in ("true", "false"):
        raise ValueError("archive row lacks its raw-derived proxy flag")
    return "below-R*" if value == "true" else "at-or-above-R*"


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
            try:
                band = _adapt_archived_proxy_band(r)
            except ValueError:
                continue
            rows.append((r["date"], int(r["qubit"]), t1, t2, band))
    return rows


def coverage_matrix(rows, which):
    """R[q, d] of rate = 1/T2 (or 1/T1), qubits with >=90% date coverage,
    dates where all selected qubits are present."""
    dates = sorted({r[0] for r in rows})
    qubits = sorted({r[1] for r in rows})
    val = {}
    for (d, q, t1, t2, _) in rows:
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


def q2_total_vs_shares(R, rng):
    total = R.sum(axis=0)
    cv_total = total.std() / total.mean()
    null = []
    for _ in range(N_PERM):
        P = np.array([rng.permutation(row) for row in R])
        t = P.sum(axis=0)
        null.append(t.std() / t.mean())
    null = np.array(null)
    pct = float((null < cv_total).mean() * 100)
    return cv_total, null, pct


def build_report():
    """Build one deterministic UTF-8/LF payload from the immutable CSV inputs."""
    lines = []
    emit = lines.append
    rng = np.random.default_rng(RNG_SEED)
    emit("QUARTER-CURRENT")
    emit("Current reading: Q3 counts an archived normalized-purity-proxy flag, not active C_Psi crossings.")
    emit("")
    emit("PTF clock field: per-site calibration-rate reading (data/ibm_history)")
    emit(f"rate = 1/T2 unless noted; permutations N={N_PERM}, seed {RNG_SEED}")
    emit(f"Q3 carries the archived raw-derived proxy flag at numerical R*={R_STAR:.17g}")
    for name, fp in FILES.items():
        if not fp.exists():
            emit(f"\n== {name}: FILE MISSING {fp.name}")
            continue
        rows = load(fp)
        R, qs, ds = coverage_matrix(rows, "T2")
        emit(f"\n== {name}: {len(rows)} records; coverage matrix {R.shape[0]} qubits x {R.shape[1]} days")
        if R.shape[1] < 10:
            emit("   too few complete days, skipping stats")
            continue

        icc, vb, vw = q1_icc(R)
        verdict1 = "SITE-OWNED (ICC > 0.5)" if icc > 0.5 else "day-noise dominated"
        emit(f"  Q1 ln-rate variance: between-qubit {vb:.4f}  within-qubit {vw:.4f}  ICC={icc:.3f}  -> {verdict1}")

        cv, null, pct = q2_total_vs_shares(R, rng)
        lo, hi = np.percentile(null, 2.5), np.percentile(null, 97.5)
        if cv < lo:
            verdict2 = "COMPENSATION: total held, shares rotate (the weave signature)"
        elif cv > hi:
            verdict2 = "COMMON-MODE chip drift (total varies MORE than independent)"
        else:
            verdict2 = "independent drift (no total conservation)"
        emit(f"  Q2 CV(chip total)={cv:.4f}  independence null [{lo:.4f}, {hi:.4f}] (pct rank {pct:.0f}%)  -> {verdict2}")

        # control with 1/T1
        R1, _, _ = coverage_matrix(rows, "T1")
        cv1, null1, _ = q2_total_vs_shares(R1, rng)
        lo1, hi1 = np.percentile(null1, 2.5), np.percentile(null1, 97.5)
        tag1 = "compensation" if cv1 < lo1 else ("common-mode" if cv1 > hi1 else "independent")
        emit(f"  Q2 control 1/T1: CV={cv1:.4f} null [{lo1:.4f}, {hi1:.4f}] -> {tag1}")

        below = {q for (_, q, _, _, band) in rows if band == "below-R*"}
        sides = {}
        for (_, q, _, _, band) in rows:
            sides.setdefault(q, set()).add(band)
        both = sum(1 for bands in sides.values() if len(bands) == 2)
        emit(f"  Q3 normalized-purity proxy: {len(below)} qubits have a below-R* row; "
             f"{both} qubits occupy both proxy bands across days")

        dln = np.abs(np.diff(np.log(R), axis=1)).ravel()
        emit(f"  telegraph |d ln rate|/day: median {np.median(dln):.3f}  p95 {np.percentile(dln, 95):.3f}")
    return "\n".join(lines) + "\n"


def main(output_path=DEFAULT_OUTPUT_PATH):
    """Write the deterministic report only to the caller-selected destination."""
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = build_report()
    destination.write_bytes(payload.encode("utf-8"))
    print(payload, end="")
    return payload


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT_PATH)
    args = parser.parse_args()
    main(output_path=args.output)
