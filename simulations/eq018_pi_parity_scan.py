#!/usr/bin/env python3
"""eq018_pi_parity_scan.py

Question: does the Pi-pair structure differ between odd N and even N?

Hypothesis: at odd N, the XY-weight midpoint N/2 is non-integer, so no
Liouvillian mode can sit at Re(lam) = -Sigma*gamma with Im(lam) = 0.
Hence zero self-Pi-modes. At even N, the midpoint is integer, so modes
with <n_XY> = N/2 exist and COULD be self-Pi if Im(lam) = 0.

Direct test: scan N in {3, 4, 5, 6}, count (pairs, self-pairs, unpaired)
and examine <n_XY> distribution.

Output: concise table. Short runtime (N=6 slowest at ~150s).
"""
from __future__ import annotations

import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.linalg import eig

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).parent))
from pi_pair_closure_investigation import (
    GAMMA_0, J_UNIFORM,
    build_H_XY, build_liouvillian_matrix,
)
from eq018_pi_pair_flow import find_pi_pairs

RESULTS_DIR = Path(__file__).parent / "results" / "eq018_pi_parity_scan"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def scan_N(N, J=J_UNIFORM, gamma_0=GAMMA_0):
    sigma_gamma = N * gamma_0
    target_sum_re = -2 * sigma_gamma

    J_list = [J] * (N - 1)
    t0 = time.time()
    L = build_liouvillian_matrix(build_H_XY(J_list, N), gamma_0, N)
    evals, _V = eig(L)
    t_diag = time.time() - t0

    pairs, self_pairs, used = find_pi_pairs(evals, target_sum_re, tol=1e-7)

    # n_XY distribution
    n_XY_values = [-float(lam.real) / (2 * gamma_0) for lam in evals]
    # Round to nearest integer (should be exactly integer up to numerical noise)
    n_XY_int = [int(round(v)) for v in n_XY_values]
    n_XY_hist = Counter(n_XY_int)

    # Self-pair details (Im(lam) of each self-Pi mode)
    self_pair_details = []
    for s in self_pairs:
        lam = evals[s]
        self_pair_details.append({
            "s": int(s),
            "Re": float(lam.real),
            "Im": float(lam.imag),
            "n_XY": float(n_XY_values[s]),
        })

    return {
        "N": N,
        "d_squared": 2 ** (2 * N),
        "sigma_gamma": sigma_gamma,
        "target_sum_re": target_sum_re,
        "n_modes": len(evals),
        "n_pairs": len(pairs),
        "n_self_pairs": len(self_pairs),
        "n_unpaired": int(np.sum(~used)),
        "diag_time_s": t_diag,
        "n_XY_histogram": {str(k): v for k, v in sorted(n_XY_hist.items())},
        "self_pair_count": len(self_pairs),
        "self_pair_details": self_pair_details[:20],  # cap for compactness
    }


def main():
    print("=" * 78)
    print("EQ-018 Pi-parity scan: even vs odd N, self-Pi-modes?")
    print("=" * 78)

    results = {}
    for N in [3, 4, 5, 6]:
        print(f"\n--- N = {N} ---")
        r = scan_N(N)
        results[N] = r
        print(f"  d^2 = {r['d_squared']},  diag in {r['diag_time_s']:.1f} s")
        print(f"  {r['n_pairs']} pairs ({r['n_pairs'] * 2} modes), "
              f"{r['n_self_pairs']} self-pairs, {r['n_unpaired']} unpaired")
        print(f"  <n_XY> histogram (value: count):")
        for k, v in r["n_XY_histogram"].items():
            print(f"    n_XY = {k}: {v}")
        if r["n_self_pairs"] > 0:
            print(f"  Self-Pi modes (first 5):")
            for sp in r["self_pair_details"][:5]:
                print(f"    index {sp['s']}: Re={sp['Re']:+.5f} "
                      f"Im={sp['Im']:+.5f} <n_XY>={sp['n_XY']:.4f}")

    # Cross-N summary
    print(f"\n{'=' * 78}")
    print(f"CROSS-N SUMMARY")
    print(f"{'=' * 78}")
    print(f"  {'N':>3} {'d^2':>6} {'parity':>8} {'N/2 integer?':>14} "
          f"{'pairs':>7} {'self':>6} {'unpaired':>10}")
    for N, r in results.items():
        parity = "odd" if N % 2 else "even"
        n_half_int = "yes" if N % 2 == 0 else "no"
        print(f"  {N:>3d} {r['d_squared']:>6d} {parity:>8} {n_half_int:>14} "
              f"{r['n_pairs']:>7d} {r['n_self_pairs']:>6d} "
              f"{r['n_unpaired']:>10d}")

    # Save
    path = RESULTS_DIR / "pi_parity_scan.json"
    with open(path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nSaved: {path}")


if __name__ == "__main__":
    main()
