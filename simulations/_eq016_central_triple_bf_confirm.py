#!/usr/bin/env python3
"""F69 central-triple: confirm finite-N -> 0.4312363 (the N->inf limit) with the VALIDATED
brute-force reduced_rho_AB (my analytic X-state reconstruction was wrong: it dropped the dN=1
coherences that the full reduced matrix carries). Symmetric (a,b,a) triple, central k=(N+1)//2.
"""
from __future__ import annotations

import sys
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from _eq016_n4_full_landscape import pair_cpsi, reduced_rho_AB, dicke_state_vector

CPSI_INF = 0.431236324950


def triple_max_bruteforce(N):
    k = (N + 1) // 2
    Dkm, Dk, Dkp = dicke_state_vector(N, k - 1), dicke_state_vector(N, k), dicke_state_vector(N, k + 1)

    def neg(t):
        a = math.sin(t) / math.sqrt(2.0); b = math.cos(t)
        psi = a * Dkm + b * Dk + a * Dkp
        psi = psi / np.linalg.norm(psi)
        return -pair_cpsi(reduced_rho_AB(psi, N))

    ts = np.linspace(0.0, math.pi, 400)
    t0 = ts[int(np.argmin([neg(t) for t in ts]))]
    res = minimize_scalar(neg, bracket=(t0 - 0.02, t0, t0 + 0.02), method="brent",
                          options={"xtol": 1e-12})
    return -res.fun


def main():
    print(f"  F69 central-triple, brute-force, toward cpsi_inf = {CPSI_INF:.10f}")
    print(f"  {'N':>4} {'cpsi(N)':>12} {'cpsi-inf':>12} {'(cpsi-inf)*N':>14}   (doc N=3..6: .8011 .7136 .6492 .6163)")
    for N in (3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13):
        cp = triple_max_bruteforce(N)
        d = cp - CPSI_INF
        print(f"  {N:>4} {cp:>12.6f} {d:>12.4e} {d * N:>14.4f}")
    print("\n  reading: if cpsi(N) decreases toward 0.43124 and (cpsi-inf)*N tends to a constant,")
    print("  the N->inf limit 0.43124 (sextic root) is confirmed with a ~1/N approach.")


if __name__ == "__main__":
    main()
