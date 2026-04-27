#!/usr/bin/env python3
"""EQ-016 sub-question: closed form for central-Dicke-triple max pair-CPsi.

Numerical sequence (max pair-CPsi on |D_{k-1}>+|D_k>+|D_{k+1}> slice):
  N=3 (k=2): 0.8011
  N=4 (k=2 or 3): 0.7136
  N=5 (k=3): 0.6492
  N=6 (k=3 or 4): 0.6163
  N=7 (k=4): ?
  N=8 (k=4): ?

Optimal coefficients consistently look like (a, b, a) with c_{k-1} = c_{k+1}
and a^2 + b^2/2 + a^2 = 1, i.e., 2a^2 + b^2 = 1. Approximate values (a, b) ≈
(0.51, 0.69) but exact optimum unclear.

Strategy:
  1. Numerically compute pair_AB(rho) for |psi(a,b,a)> in the central-Dicke
     triple, parametrise (a, b) on unit hemisphere.
  2. Find max with high precision.
  3. Verify the (a, b) form is exact (i.e. c_{k-1} = c_{k+1}) by comparing
     to fully-free (c_{k-1}, c_k, c_{k+1}) optimum.
  4. Look for a closed form: try (a, b) = (1/2, 1/sqrt(2)), then small
     deviations; try exact rational forms.
  5. Extend numerically to N = 7, 8, 9, 10 and try to fit cpsi(N) =
     f(N) as a rational, asymptotic, or algebraic expression.
"""
from __future__ import annotations

import math
import sys
from math import comb
from pathlib import Path

import numpy as np
from scipy.optimize import minimize, brentq

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from _eq016_n4_full_landscape import (
    dicke_state_vector, reduced_rho_AB, pair_cpsi, single_qubit_purity,
)


def cpsi_central_symmetric(N, k, a, b):
    """pair-CPsi of |psi> = a |D_{k-1}> + b |D_k> + a |D_{k+1}> normalized."""
    Dkm = dicke_state_vector(N, k - 1)
    Dk = dicke_state_vector(N, k)
    Dkp = dicke_state_vector(N, k + 1)
    psi = a * Dkm + b * Dk + a * Dkp
    norm = np.linalg.norm(psi)
    if norm < 1e-12:
        return 0.0
    psi /= norm
    return pair_cpsi(reduced_rho_AB(psi, N))


def cpsi_central_free(N, k, c1, c2, c3):
    """pair-CPsi of c1 |D_{k-1}> + c2 |D_k> + c3 |D_{k+1}>."""
    psi = (c1 * dicke_state_vector(N, k - 1)
           + c2 * dicke_state_vector(N, k)
           + c3 * dicke_state_vector(N, k + 1))
    norm = np.linalg.norm(psi)
    if norm < 1e-12:
        return 0.0
    psi /= norm
    return pair_cpsi(reduced_rho_AB(psi, N))


def find_max_symmetric(N, k):
    """Find max cpsi on (a, b, a) slice. Param t: a = sin(t)/sqrt(2), b = cos(t)."""
    def neg(t):
        a = math.sin(t) / math.sqrt(2)
        b = math.cos(t)
        return -cpsi_central_symmetric(N, k, a, b)

    # grid search first
    ts = np.linspace(0, math.pi, 1001)
    vals = [neg(t) for t in ts]
    best_t = ts[np.argmin(vals)]
    # refine
    res = minimize(lambda x: neg(x[0]), [best_t], method="Nelder-Mead",
                    options={"xatol": 1e-12, "fatol": 1e-15})
    t_opt = res.x[0]
    a_opt = math.sin(t_opt) / math.sqrt(2)
    b_opt = math.cos(t_opt)
    cpsi_opt = -res.fun
    return {
        "N": N, "k": k,
        "t_opt": t_opt, "a_opt": a_opt, "b_opt": b_opt,
        "cpsi_opt": cpsi_opt,
        "2a²+b²": 2 * a_opt ** 2 + b_opt ** 2,
        "a²/b²": a_opt ** 2 / b_opt ** 2 if abs(b_opt) > 1e-12 else float('inf'),
    }


def find_max_free(N, k):
    """Find max cpsi on full (c_{k-1}, c_k, c_{k+1}) slice (3D, real)."""
    def neg(c):
        return -cpsi_central_free(N, k, c[0], c[1], c[2])

    # multi-start
    best = (-1, None)
    for seed in range(50):
        rng = np.random.default_rng(seed)
        c0 = rng.normal(size=3)
        c0 /= np.linalg.norm(c0)
        res = minimize(neg, c0, method="Nelder-Mead",
                        options={"xatol": 1e-12, "fatol": 1e-15})
        if -res.fun > best[0]:
            best = (-res.fun, res.x / np.linalg.norm(res.x))
    return {
        "N": N, "k": k,
        "c": best[1].tolist(),
        "cpsi_opt": best[0],
    }


def main():
    print("=" * 80)
    print("Central-Dicke-triple maximum pair-CPsi: closed form search")
    print("=" * 80)
    print()
    print("Param: |psi> = c1 |D_{k-1}> + c2 |D_k> + c3 |D_{k+1}>, k around N/2.")
    print()

    # Step 1: verify symmetric (a, b, a) is correct ansatz
    print("Step 1: verify c_{k-1} = c_{k+1} ansatz vs fully free")
    print("-" * 80)
    print(f"  {'N':>3} {'k':>3} {'free max':>12} {'symmetric max':>14} {'gap':>10}")
    print("  " + "-" * 50)
    for N in [3, 4, 5, 6, 7, 8]:
        # Use central k
        k = (N + 1) // 2
        if k < 1 or k > N - 1:
            continue
        sym = find_max_symmetric(N, k)
        free = find_max_free(N, k)
        gap = abs(free["cpsi_opt"] - sym["cpsi_opt"])
        print(f"  {N:>3d} {k:>3d} {free['cpsi_opt']:>12.8f} "
              f"{sym['cpsi_opt']:>14.8f} {gap:>10.2e}")
    print()

    # Step 2: high-precision symmetric optimum
    print("Step 2: high-precision symmetric (a, b, a) optimum")
    print("-" * 80)
    print(f"  {'N':>3} {'k':>3} {'a':>12} {'b':>12} {'a/b':>10} "
          f"{'cpsi':>14}")
    print("  " + "-" * 70)
    results = []
    for N in [3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 20]:
        k = (N + 1) // 2
        if k < 1 or k > N - 1:
            continue
        sym = find_max_symmetric(N, k)
        results.append(sym)
        print(f"  {N:>3d} {k:>3d} {sym['a_opt']:>12.8f} {sym['b_opt']:>12.8f} "
              f"{sym['a_opt']/sym['b_opt']:>10.6f} {sym['cpsi_opt']:>14.8f}")
    print()

    # Step 3: pattern in (a/b) ratio and cpsi
    print("Step 3: scaling patterns")
    print("-" * 80)
    print(f"  {'N':>3} {'cpsi':>12} {'1/cpsi':>10} {'cpsi*N':>10} "
          f"{'1-cpsi':>10} {'(1-cpsi)*N':>12} {'a²/b²':>10}")
    print("  " + "-" * 75)
    for r in results:
        cp = r["cpsi_opt"]
        N = r["N"]
        ab2 = r["a_opt"] ** 2 / r["b_opt"] ** 2
        print(f"  {N:>3d} {cp:>12.6f} {1/cp:>10.4f} {cp*N:>10.4f} "
              f"{1-cp:>10.4f} {(1-cp)*N:>12.4f} {ab2:>10.4f}")
    print()

    # Try common closed forms and check residual
    print("Step 4: try candidate closed forms")
    print("-" * 80)
    candidates = [
        ("3/(N+1)", lambda N: 3.0 / (N + 1)),
        ("(N-1)/(N+1)", lambda N: (N - 1) / (N + 1)),
        ("4/(N+2)", lambda N: 4.0 / (N + 2)),
        ("1 - 1/N", lambda N: 1 - 1.0 / N),
        ("1/(1+log(N)/log(2))", lambda N: 1 / (1 + math.log(N, 2))),
        ("(N-1)^2/N(N+1)", lambda N: (N - 1) ** 2 / (N * (N + 1))),
    ]
    for name, formula in candidates:
        max_err = max(abs(r["cpsi_opt"] - formula(r["N"])) for r in results)
        print(f"  cpsi(N) = {name:<25}: max_err over tested N = {max_err:.4e}")
    print()


if __name__ == "__main__":
    main()
