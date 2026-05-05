"""Qubit-biography analysis: per-qubit lifecycle phases over an IBM-history CSV.

Reads `data/ibm_history/ibm_<backend>_history.csv` (one row per qubit-day) and
classifies each qubit's daily r = T2/(2T1) trajectory into phases. Companion to
`docs/BOTH_SIDES_VISIBLE.md` which sketched Q98's "tune → pulse → fade" lifecycle
visually; this script makes the segmentation explicit.

Phase classification (sliding window, default 7 days):

- **silent**:     mean r ≥ R* + 0.10, std small  → stably classical-side, far from boundary
- **classic**:    mean r ≥ R*, std small         → classical-side, near or at boundary
- **pulse**:      mean r < R*, std small         → quantum-side, sustained crossing
- **twitch**:     std large, near boundary       → boundary-fluctuating
- **fade**:       trend crossing R* upward       → leaving quantum side
- **tune**:       trend crossing R* downward     → entering quantum side
- **anomaly**:    sudden T2 drop beyond ±3σ      → TLS hit / spike

Outputs:
  1. Per-qubit phase sequence (compact "fingerprint" character string)
  2. Aggregate phase-fraction table (which qubits live where)
  3. Co-movement clusters: pairs whose r time series correlate at ρ > 0.5
     (suggests shared physical mechanism, e.g. TLS defect, cooling line)
  4. Lifecycle archetype tally: how many "pulse-only", "lifecycle (tune→pulse→fade)",
     "twitch-dominant", etc.

Reading: this is the per-qubit version of the Marrakesh 5-day flipping count
(33/156 = 21% in `_marrakesh_quarter_boundary_review.py`). At 6-month depth on
Torino we expect Q98's lifecycle to surface as "tune → pulse → fade → silent",
and Q72's rhythm as "pulse → twitch → pulse" alternation.
"""
from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_HISTORY = REPO_ROOT / "data" / "ibm_history" / "ibm_torino_history.csv"

R_STAR = 0.212755  # canonical CΨ=¼ boundary, matching QubitRegime.R_STAR in C# Calibration cockpit
SILENT_DELTA = 0.10        # mean r ≥ R*+δ_silent → silent (far above boundary)
WINDOW = 7                 # days per sliding window
STD_HIGH = 0.05            # window std above this = "twitch" / volatile
TREND_THRESHOLD = 0.02     # window-end minus window-start > this = trending
ANOMALY_SIGMA = 3.0        # |Δr| > σ·trajectory_std → anomaly day

PHASE_CHARS = {
    "silent":  ".",
    "classic": "_",
    "pulse":   "X",
    "twitch":  "~",
    "fade":    "/",
    "tune":    "\\",
    "anomaly": "!",
    "missing": " ",
}


def load_history(csv_path: Path) -> Dict[int, List[Tuple[str, float, float, float]]]:
    """Return {qubit_id: [(date, T1, T2, r), ...]} sorted by date."""
    by_qubit: Dict[int, List[Tuple[str, float, float, float]]] = defaultdict(list)
    with open(csv_path, "r", encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            try:
                qid = int(row["qubit"])
                t1 = float(row["T1_us"]) if row.get("T1_us") else 0.0
                t2 = float(row["T2_us"]) if row.get("T2_us") else 0.0
                r = float(row["r_param"]) if row.get("r_param") else 0.0
                if t1 > 0 and t2 > 0:
                    by_qubit[qid].append((row["date"], t1, t2, r))
            except (ValueError, KeyError):
                continue
    for qid in by_qubit:
        by_qubit[qid].sort(key=lambda rec: rec[0])
    return by_qubit


def classify_window(window_r: np.ndarray) -> str:
    """Classify a sliding window of r values into a phase tag."""
    if len(window_r) < 3:
        return "missing"
    mean = window_r.mean()
    std = window_r.std()
    trend = window_r[-1] - window_r[0]
    if std > STD_HIGH:
        if abs(mean - R_STAR) < 0.05:
            return "twitch"
    if std > STD_HIGH and trend > TREND_THRESHOLD and window_r[0] < R_STAR <= window_r[-1]:
        return "fade"
    if std > STD_HIGH and trend < -TREND_THRESHOLD and window_r[0] > R_STAR >= window_r[-1]:
        return "tune"
    if mean < R_STAR:
        return "pulse"
    if mean >= R_STAR + SILENT_DELTA:
        return "silent"
    return "classic"


def biography(rs: np.ndarray) -> List[str]:
    """Per-day phase sequence (one tag per day, using sliding window centered at i)."""
    n = len(rs)
    half = WINDOW // 2
    phases: List[str] = []
    overall_std = rs.std()
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        win = rs[lo:hi]
        # anomaly check: this day's |delta from local mean| > 3*overall_std
        if i > 0 and abs(rs[i] - win.mean()) > ANOMALY_SIGMA * overall_std and overall_std > 0.01:
            phases.append("anomaly")
        else:
            phases.append(classify_window(win))
    return phases


def fingerprint(phases: List[str]) -> str:
    """One character per day; compresses to a fingerprint string."""
    return "".join(PHASE_CHARS[p] for p in phases)


def archetype_from_series(rs: np.ndarray) -> str:
    """Reduce a raw r-series to a coarse archetype label using boundary-walk
    rate, crossing fraction, and trajectory volatility. Calibrated against the
    BOTH_SIDES_VISIBLE.md named examples (Q72/Q98 twitch, Q80 pulse-stable,
    Q105 lifecycle, Q70/Q68 mostly silent with crossings)."""
    n = len(rs)
    if n < 2:
        return "empty"
    crossing = (rs < R_STAR).mean()
    signs = np.sign(rs - R_STAR)
    walk = int(np.sum(np.diff(signs) != 0)) / (n - 1)
    std = float(rs.std())

    if walk > 0.20:
        return "twitch"
    if walk > 0.05:
        return "lifecycle"
    if crossing > 0.7:
        return "pulse-stable"
    if crossing < 0.1:
        return "silent-stable" if std < 0.10 else "drifty-silent"
    return "classic-stable"


def co_movement(by_qubit: Dict[int, List], min_days: int = 30, threshold: float = 0.5):
    """Pairs of qubits whose r time series correlate above `threshold`. Restricted
    to pairs sharing at least `min_days` of overlapping calibration days."""
    qids = sorted(by_qubit.keys())
    series = {qid: {rec[0]: rec[3] for rec in by_qubit[qid]} for qid in qids}
    pairs = []
    for i, a in enumerate(qids):
        for b in qids[i + 1:]:
            common = sorted(set(series[a]) & set(series[b]))
            if len(common) < min_days:
                continue
            ra = np.array([series[a][d] for d in common])
            rb = np.array([series[b][d] for d in common])
            if ra.std() < 0.01 or rb.std() < 0.01:
                continue
            rho = np.corrcoef(ra, rb)[0, 1]
            if abs(rho) > threshold:
                pairs.append((a, b, rho, len(common)))
    pairs.sort(key=lambda p: -abs(p[2]))
    return pairs


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", nargs="?", type=str, default=str(DEFAULT_HISTORY))
    parser.add_argument("--show", type=int, default=15, help="how many qubit fingerprints to display")
    parser.add_argument("--corr-threshold", type=float, default=0.5)
    args = parser.parse_args()

    csv_path = Path(args.csv)
    print(f"Loading {csv_path.name}")
    by_qubit = load_history(csv_path)
    n_q = len(by_qubit)
    n_days = max(len(v) for v in by_qubit.values())
    print(f"  {n_q} qubits, up to {n_days} days each")
    print()

    # Compute biographies
    bios: Dict[int, Tuple[List[str], str]] = {}
    for qid, recs in by_qubit.items():
        rs = np.array([rec[3] for rec in recs])
        phases = biography(rs)
        bios[qid] = (phases, archetype_from_series(rs))

    # Archetype tally
    arch_count: Dict[str, int] = defaultdict(int)
    for _, arch in bios.values():
        arch_count[arch] += 1
    print("Archetype tally:")
    for arch, n in sorted(arch_count.items(), key=lambda kv: -kv[1]):
        print(f"  {arch:<20} {n:3d} qubits")
    print()

    # Show fingerprints for the most "interesting" qubits (lifecycle / twitch / drifty)
    interesting_archs = ("lifecycle", "twitch", "drifty-silent")
    interesting = [(qid, ph, ar) for qid, (ph, ar) in bios.items() if ar in interesting_archs]
    print(f"Interesting qubits ({len(interesting)} total):")
    print(f"  legend: . silent  _ classic  X pulse  ~ twitch  / fade  \\ tune  ! anomaly")
    print()
    for qid, phases, arch in interesting[:args.show]:
        fp = fingerprint(phases)
        chunks = [fp[i:i+30] for i in range(0, len(fp), 30)]
        print(f"  Q{qid:<3} [{arch:<14}] ({len(phases)} days)")
        for chunk in chunks:
            print(f"    {chunk}")
    print()

    stable_pulse = [(qid, ph) for qid, (ph, ar) in bios.items() if ar == "pulse-stable"]
    print(f"Stable pulse qubits (always quantum-side): {len(stable_pulse)}")
    for qid, phases in stable_pulse[:5]:
        fp = fingerprint(phases)
        print(f"  Q{qid:<3}: {fp[:90]}{'...' if len(fp) > 90 else ''}")
    print()

    # Co-movement clusters
    print(f"Co-movement (r-trajectory correlation > {args.corr_threshold}):")
    pairs = co_movement(by_qubit, threshold=args.corr_threshold)
    print(f"  {len(pairs)} significant pairs out of {n_q * (n_q - 1) // 2} possible")
    print()
    if pairs:
        print(f"  {'qubit A':>8} {'qubit B':>8} {'ρ':>8} {'n days':>8}")
        for a, b, rho, n in pairs[:20]:
            print(f"  {a:>8} {b:>8} {rho:>+8.3f} {n:>8}")
    print()

    print("Reading:")
    print(f"  - {arch_count.get('silent-stable', 0)} qubits never approach CΨ=¼ ({100*arch_count.get('silent-stable',0)/n_q:.0f}%)")
    print(f"  - {arch_count.get('pulse-stable', 0)} qubits live permanently in the quantum-side regime ({100*arch_count.get('pulse-stable',0)/n_q:.0f}%)")
    print(f"  - {arch_count.get('lifecycle', 0)} qubits show a slow tune-pulse-fade arc (Q105 archetype)")
    print(f"  - {arch_count.get('twitch', 0)} qubits are boundary-twitchers (Q72/Q98 archetype)")
    print(f"  - {arch_count.get('drifty-silent', 0)} qubits are mostly above the boundary but with large excursions")
    print(f"  - Co-movement pairs are candidates for shared physical mechanism")
    print(f"    (TLS defects coupling adjacent qubits, common drive-line drift).")


if __name__ == "__main__":
    main()
