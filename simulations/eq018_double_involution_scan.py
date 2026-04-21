#!/usr/bin/env python3
"""eq018_double_involution_scan.py

Test the mirror-axis principle from ORTHOGONALITY_SELECTION_FAMILY §5a:

    Self-Pi modes in the Liouvillian at n_XY = N/2 sector appear when the
    XY chain's single-excitation H-eigenvalues {E_k = 2*cos(pi*k/(N+1))}
    carry a "double involution":
      - multiplicative: x * y = 1 for some pair (x, y) in the spectrum
      - additive: x - y = 1 for some pair (x, y) in the spectrum
    AND N is even (so the midpoint N/2 is an integer, enabling integer
    n_XY midpoint modes).

The unique real numbers satisfying x*(1/x) = 1 AND x - (1/x) = 1 are
x = (1 +/- sqrt(5))/2, i.e., phi and -1/phi. So the "double involution"
is exactly the Golden Ratio structure.

When does 2*cos(pi*k/(N+1)) = phi? Answer: when pi*k/(N+1) = pi/5, so
k/(N+1) = 1/5, so N+1 is a multiple of 5. That gives N+1 in {5, 10, 15, ...}
i.e., N in {4, 9, 14, 19, ...}.

Intersect with N even: N in {4, 14, 24, 34, ...}.

Prediction: self-Pi modes exist ONLY at N in {4, 14, 24, ...} in the
XY chain OBC setting. All other even N (6, 8, 10, 12, 16, 18, 20, 22, ...)
have no self-Pi modes at the midpoint sector.

This script verifies the prediction by scanning N = 3..30, listing for
each N:
  - H-eigenvalues
  - whether multiplicative pairs (product = 1) exist
  - whether additive pairs (difference = 1) exist
  - whether double-involution holds
  - whether N is even

Then cross-checks against the parity-scan data for N = 3, 4, 5, 6.
"""
from __future__ import annotations

import json
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

RESULTS_DIR = Path(__file__).parent / "results" / "eq018_double_involution_scan"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


PHI = (1 + np.sqrt(5)) / 2
ONE_OVER_PHI = 1 / PHI  # = phi - 1 = 0.618...


def h_single_excitation_eigenvalues(N):
    return [2 * np.cos(np.pi * k / (N + 1)) for k in range(1, N + 1)]


def find_multiplicative_pairs(eigvals, target=1.0, tol=1e-9):
    pairs = []
    for i, j in combinations(range(len(eigvals)), 2):
        if abs(eigvals[i] * eigvals[j] - target) < tol:
            pairs.append((i + 1, j + 1, float(eigvals[i]), float(eigvals[j])))
    return pairs


def find_additive_pairs(eigvals, diff=1.0, tol=1e-9):
    pairs = []
    n = len(eigvals)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if abs(eigvals[i] - eigvals[j] - diff) < tol:
                pairs.append((i + 1, j + 1, float(eigvals[i]), float(eigvals[j])))
    return pairs


def check_golden_ratio_presence(eigvals, tol=1e-9):
    """Does {phi, 1/phi} subset the spectrum (up to sign)?"""
    has_phi = any(abs(e - PHI) < tol or abs(e + PHI) < tol for e in eigvals)
    has_inv_phi = any(abs(e - ONE_OVER_PHI) < tol
                      or abs(e + ONE_OVER_PHI) < tol for e in eigvals)
    return has_phi, has_inv_phi


def main():
    start = time.time()
    print("=" * 78)
    print("Double-involution scan for XY-chain single-excitation eigenvalues")
    print("=" * 78)
    print(f"  phi = {PHI:.10f}")
    print(f"  1/phi = {ONE_OVER_PHI:.10f}")
    print()

    results = []
    for N in range(3, 31):
        eigvals = h_single_excitation_eigenvalues(N)
        mult_pairs = find_multiplicative_pairs(eigvals)
        add_pairs = find_additive_pairs(eigvals)
        has_phi, has_inv_phi = check_golden_ratio_presence(eigvals)
        even = (N % 2 == 0)
        double_invol = len(mult_pairs) > 0 and len(add_pairs) > 0
        predict_self_pi = double_invol and even

        results.append({
            "N": N,
            "N_plus_1": N + 1,
            "N_plus_1_mod_5": (N + 1) % 5,
            "even": even,
            "eigvals": [float(e) for e in eigvals],
            "has_phi": has_phi,
            "has_inv_phi": has_inv_phi,
            "multiplicative_pairs_count": len(mult_pairs),
            "additive_pairs_count": len(add_pairs),
            "double_involution": double_invol,
            "predict_self_pi": predict_self_pi,
        })

    # Cross-check with parity scan data
    parity_scan_known = {3: 0, 4: 18, 5: 0, 6: 0}  # from eq018_pi_parity_scan

    print(f"  {'N':>3} {'even?':>6} {'(N+1)%5':>9} {'has phi?':>10} "
          f"{'has 1/phi?':>11} {'mult pairs':>12} {'add pairs':>10} "
          f"{'predict':>9} {'actual':>8}")
    for r in results:
        actual = parity_scan_known.get(r["N"], "?")
        actual_s = str(actual) if actual != "?" else "?"
        print(f"  {r['N']:>3d} {str(r['even']):>6} {r['N_plus_1_mod_5']:>9d} "
              f"{str(r['has_phi']):>10} {str(r['has_inv_phi']):>11} "
              f"{r['multiplicative_pairs_count']:>12d} "
              f"{r['additive_pairs_count']:>10d} "
              f"{str(r['predict_self_pi']):>9} {actual_s:>8}")

    # Summary
    predicted_yes = [r["N"] for r in results if r["predict_self_pi"]]
    print(f"\n  N values predicted to have self-Pi modes (double involution + even):")
    print(f"    {predicted_yes}")
    print(f"  This should be {{4, 14, 24, ...}} if N+1 mod 5 = 0 and N even.")
    expected = [N for N in range(3, 31) if (N + 1) % 5 == 0 and N % 2 == 0]
    print(f"  Expected (N even AND N+1 divisible by 5): {expected}")

    # Full details for the first few N with double involution
    print(f"\n  Detailed pair listing for predicted-self-Pi N values:")
    for r in results:
        if not r["predict_self_pi"]:
            continue
        N = r["N"]
        eigvals = r["eigvals"]
        print(f"\n    N = {N}:")
        print(f"      eigvals = {[f'{e:+.4f}' for e in eigvals]}")
        mults = find_multiplicative_pairs(eigvals)
        adds = find_additive_pairs(eigvals)
        for (i, j, e_i, e_j) in mults[:3]:
            print(f"      E_{i}·E_{j} = {e_i:+.4f} * {e_j:+.4f} = 1")
        for (i, j, e_i, e_j) in adds[:3]:
            print(f"      E_{i} - E_{j} = {e_i:+.4f} - {e_j:+.4f} = 1")

    # Save
    out = {
        "scan_range": [3, 30],
        "phi": PHI,
        "one_over_phi": ONE_OVER_PHI,
        "results": results,
        "predicted_self_pi_N": predicted_yes,
        "parity_scan_known": parity_scan_known,
    }
    path = RESULTS_DIR / "double_involution_scan.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nSaved: {path}")

    # Final verdict
    print(f"\n{'=' * 78}")
    print(f"VERDICT")
    print(f"{'=' * 78}")
    all_match = True
    for N, actual in parity_scan_known.items():
        r = next(rr for rr in results if rr["N"] == N)
        predicted = r["predict_self_pi"]
        actual_has = actual > 0
        match = predicted == actual_has
        if not match:
            all_match = False
        print(f"  N={N}: predicted_self_pi = {predicted}, actual_self_pi_count = {actual}, "
              f"match = {match}")
    print(f"\n  All parity-scan data matches the double-involution prediction: "
          f"{all_match}")
    print(f"\n  Testable: N = 14 (even, N+1=15 divisible by 5). Prediction: "
          f"self-Pi modes present.")
    print(f"  Not testable here (d^2 = 2^28 = 2.7e8, too large for Python dense eig).")
    print(f"  N = 8 (even, N+1=9 NOT divisible by 5): prediction: NO self-Pi modes.")
    print(f"  d^2 = 65536 at N=8, also too large for Python but feasible for C# engine.")

    print(f"\nWalltime: {time.time() - start:.3f} s")


if __name__ == "__main__":
    main()
