"""Wave-breaking N-scan: chain OBC spectrum constructibility → Liouvillian Im=0 cancellations.

EQ-018 (2026, earlier) observed that N=4 OBC XY chain under Z-dephasing has 18
modes with exact Im(λ)=0 in the n_XY=2 sector, while N=6 has 0. The N=4
explanation: H_SE eigenvalues are 2cos(πk/5) = ±φ, ±1/φ, golden-ratio family
— a constructible algebraic field Q[√5] where exact cancellations
φ+(−φ)=0, 1/φ+(−1/φ)=0 happen. N=6 H_SE eigenvalues are 2cos(πk/7), a
NON-constructible heptagonal field where no such cancellations exist.

Tonight (2026-05-17) we derived F99 + Niven-completeness + off-Niven=wave-
breaking. The structural reading predicts which N have wave-breaking and which
do not, based ENTIRELY on whether π/(N+1) is constructible by compass+straightedge:

  N=2: π/3 = 60°       constructible Niven       → many Im=0 (integer ring)
  N=3: π/4 = 45°       constructible Niven       → many Im=0 (Q[√2])
  N=4: π/5 = 36°       constructible off-Niven   → many Im=0 (Q[√5], golden)  ← VERIFIED
  N=5: π/6 = 30°       constructible Niven       → many Im=0 (Q[√3])
  N=6: π/7 ≈ 25.7°     NOT constructible (heptagon) → 0 Im=0  ← VERIFIED
  N=7: π/8 = 22.5°     constructible off-Niven   → many Im=0 (Q[√2, √(2+√2)]) (compute-prohibitive here)
  N=8: π/9 = 20°       NOT constructible (nonagon) → 0 Im=0 (compute-prohibitive here)

This script tests the prediction at N=2, 3, 5 (which EQ-018 didn't cover) to
mature the connection from one anchor (N=4) to a pattern. N=4 and N=6 are re-
run for self-consistency.

Tom's structural reading (2026-05-17): "an dem Punkt wo der Goldene Schnitt
im Repo erscheint, wurde vorher auf irgendeine Art eine Welle gebrochen, ein
Mechanismus?" This script verifies whether the mechanism is constructibility
of the chain spectrum's underlying angle field.

Run:
  PYTHONIOENCODING=utf-8 python simulations/wave_breaking_chain_n_scan.py
"""
from __future__ import annotations

import json
import sys
import time
from math import cos, pi
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

sys.path.insert(0, str(Path(__file__).parent))
from pi_pair_closure_investigation import (
    GAMMA_0, J_UNIFORM,
    build_H_XY, build_liouvillian_matrix,
)

RESULTS_DIR = Path(__file__).parent / "results" / "wave_breaking_chain_n_scan"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

PHI = (1 + np.sqrt(5)) / 2


CHAIN_ANGLE_CLASS = {
    # N: (angle π/(N+1) in deg, constructibility class, expected field, prediction)
    2: (60.0,    "Niven",                "Z (integer)",          "many"),
    3: (45.0,    "Niven",                "Q[√2]",                "many"),
    4: (36.0,    "off-Niven golden",     "Q[√5]",                "many"),
    5: (30.0,    "Niven",                "Q[√3]",                "many"),
    6: (180/7,   "NOT constructible",    "Q[2cos(2π/7)] septagon", "zero"),
}


def closed_form_h_eigenvalues(N: int) -> list[float]:
    """OBC tridiagonal H_SE = J·(off-diag) eigenvalues = 2·J·cos(πk/(N+1))."""
    return sorted([2 * J_UNIFORM * cos(pi * k / (N + 1)) for k in range(1, N + 1)])


def identify_h_spectrum(N: int) -> str:
    """Closed-form algebraic label for the H_SE spectrum at this N."""
    if N == 2:
        return "{−1, +1}"
    if N == 3:
        return "{−√2, 0, +√2}"
    if N == 4:
        return "{−φ, −1/φ, +1/φ, +φ}  (golden)"
    if N == 5:
        return "{−√3, −1, 0, +1, +√3}"
    if N == 6:
        return "{2·cos(πk/7) for k=1..6}  (heptagonal, NOT constructible)"
    return "?"


def scan_one_N(N: int, im_tol: float = 1e-9) -> dict:
    """Build L_A, diagonalize, count Im=0 modes per Re sector."""
    print(f"\n{'='*86}")
    print(f"  N = {N}  (chain angle π/(N+1) = {180/(N+1):.4f}°, "
          f"class: {CHAIN_ANGLE_CLASS[N][1]})")
    print(f"{'='*86}")

    h_evals = closed_form_h_eigenvalues(N)
    print(f"  H_SE closed-form eigenvalues: {[f'{v:+.6f}' for v in h_evals]}")
    print(f"  Algebraic identity:           {identify_h_spectrum(N)}")
    print(f"  Predicted Im=0 cancellations: {CHAIN_ANGLE_CLASS[N][3].upper()}")

    t0 = time.time()
    J_list = [J_UNIFORM] * (N - 1)
    H = build_H_XY(J_list, N)
    L = build_liouvillian_matrix(H, GAMMA_0, N)
    print(f"  L_A size: {L.shape[0]}×{L.shape[0]} "
          f"({L.shape[0]**2 * 16 / 1e9:.3f} GB)")

    ev = np.linalg.eigvals(L)
    walltime = time.time() - t0
    print(f"  L_A diagonalised in {walltime:.2f} s ({len(ev)} eigenvalues)")

    im_abs = np.abs(np.imag(ev))
    n_zero_im = int(np.sum(im_abs < im_tol))
    n_total = len(ev)
    frac_zero = n_zero_im / n_total

    print(f"  Im(λ)=0 modes: {n_zero_im} of {n_total} total ({100*frac_zero:.2f}%)")

    # Bucket by Re(λ) at sigma_gamma scale precision
    re_vals = np.real(ev)
    re_rounded = np.round(re_vals / GAMMA_0).astype(int)
    sectors = {}
    for r, im in zip(re_rounded, im_abs):
        if r not in sectors:
            sectors[r] = {"count": 0, "zero_im": 0}
        sectors[r]["count"] += 1
        if im < im_tol:
            sectors[r]["zero_im"] += 1

    print(f"  Per Re-sector breakdown (Re ≈ k·γ_0 = {GAMMA_0}):")
    for r in sorted(sectors.keys()):
        s = sectors[r]
        if s["zero_im"] > 0 or s["count"] > 16:
            print(f"    Re ≈ {r}·γ_0:  {s['count']:5d} modes,  "
                  f"{s['zero_im']:4d} with Im=0  "
                  f"({100*s['zero_im']/s['count']:5.1f}%)")

    return {
        "N": N,
        "angle_deg": 180 / (N + 1),
        "class": CHAIN_ANGLE_CLASS[N][1],
        "prediction": CHAIN_ANGLE_CLASS[N][3],
        "h_evals": h_evals,
        "h_identity": identify_h_spectrum(N),
        "n_total_modes": n_total,
        "n_zero_im": n_zero_im,
        "fraction_zero_im": frac_zero,
        "walltime_s": walltime,
        "sectors": {int(r): s for r, s in sectors.items()},
    }


def print_summary_table(results: list[dict]):
    print()
    print("=" * 86)
    print("  SUMMARY: chain angle constructibility → wave-breaking (Im=0) prediction test")
    print("=" * 86)
    print()
    print(f"  {'N':>3} {'π/(N+1)':>10} {'class':<22} {'predict':<8} "
          f"{'#Im=0':>8} {'total':>8} {'fraction':>10} {'match?'}")
    print(f"  {'-'*3} {'-'*10} {'-'*22} {'-'*8} {'-'*8} {'-'*8} {'-'*10} {'-'*7}")
    for r in results:
        N = r["N"]
        frac = r["fraction_zero_im"]
        predicted = r["prediction"]
        if predicted == "many" and frac > 0.01:
            match = "✓ YES"
        elif predicted == "zero" and frac < 0.001:
            match = "✓ YES"
        elif predicted == "many" and frac <= 0.01:
            match = "✗ no"
        elif predicted == "zero" and frac >= 0.001:
            match = "✗ no"
        else:
            match = "?"
        print(f"  {N:>3} {r['angle_deg']:>9.4f}° {r['class']:<22} {predicted:<8} "
              f"{r['n_zero_im']:>8d} {r['n_total_modes']:>8d} "
              f"{frac:>10.4f} {match}")
    print()
    print("  Reading:")
    print("    'many' prediction = chain angle π/(N+1) is constructible → Liouvillian")
    print("    Im(λ)=0 cancellations are algebraically supported (rational sum to 0)")
    print()
    print("    'zero' prediction = chain angle is NOT constructible (heptagon, nonagon)")
    print("    → cyclotomic field has no rational-coefficient kernel → generic Im≠0")
    print()
    print("  If all rows match, the wave-breaking mechanism is confirmed structurally:")
    print("  off-Niven angles produce wave-breaking iff they are still CONSTRUCTIBLE.")
    print("  Niven angles are a strict subset; the broader criterion is constructibility")
    print("  of π·k/(N+1) (= Gauss-Wantzel: N+1 product of distinct Fermat primes × 2^k).")


def main():
    print()
    print("=" * 86)
    print("  Wave-breaking N-scan: testing the constructibility → Im=0 prediction")
    print("=" * 86)
    print()
    print("  Hypothesis (Tom 2026-05-17): OFF-Niven constructible angles populate the")
    print("  same wave-breaking landscape as Niven angles, just via algebraic-irrational")
    print("  α values (golden / silver / √3 families). Test: does the prediction hold")
    print("  on the EQ-018 chain spectrum scan at N = 2, 3, 5 (uncovered) and 4, 6 (anchors)?")

    results = []
    for N in [2, 3, 4, 5, 6]:
        r = scan_one_N(N)
        results.append(r)

    print_summary_table(results)

    # Save
    out_path = RESULTS_DIR / "chain_n_scan.json"
    serializable = []
    for r in results:
        s = dict(r)
        s["h_evals"] = [float(v) for v in s["h_evals"]]
        serializable.append(s)
    with open(out_path, "w") as f:
        json.dump(serializable, f, indent=2)
    print(f"\n  Saved: {out_path}")


if __name__ == "__main__":
    main()
