#!/usr/bin/env python3
"""EQ-016 sub-question: asymptotic form of central-Dicke-triple max pair-CPsi.

Fast evaluator: compute ρ_AB analytically from binomial coefficients, no
need to instantiate the 2^N state vector. Allows N up to ~1000.

Setup:
  |psi> = a |D_{k-1}> + b |D_k> + a |D_{k+1}>,  k = (N+1)//2 (central),
  norm² = 2a² + b² = 1 (when normalized).

Pair-AB on qubits 0, 1 has structure (real symmetric):
  d_0 = rho[00,00] = a² · D[k-1, 0] + b² · D[k, 0] + a² · D[k+1, 0]
  d_1 = rho[01,01] = rho[10,10] = rho[01,10] = a² · D[k-1, 1] + ... (s=1 entry)
  d_2 = rho[11,11] = a² · D[k-1, 2] + ...
  e_01 = rho[00,01] = rho[00,10] = rho[01,11] = rho[10,11]
       = ab · [E[k-1, k] + E[k, k+1]]    (s=0↔1 sector)
  e_02 = rho[00,11] = a² · E[k-1, k+1]    (s=0↔2 sector)

where D[n, s] = C(N-2, n-s) / C(N, n) (sector s diagonal weight on Dicke n)
      E[n, m] = C(N-2, n-(s)) / sqrt(C(N,n)·C(N,m)) for matching s sector.

Tr(ρ²) = d_0² + 4·d_1² + d_2² + 8·e_01² + 2·e_02²
L1_off = 2·d_1 + 8·e_01 + 2·e_02     (all positive for real (a, b) > 0)
cpsi = Tr(ρ²) · L1_off / 3
"""
from __future__ import annotations

import math
import sys
from math import comb, log, exp
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log_binom(n, r):
    """Log of binomial coefficient C(n, r). For r out of range returns -inf."""
    if r < 0 or r > n:
        return float('-inf')
    return (math.lgamma(n + 1) - math.lgamma(r + 1) - math.lgamma(n - r + 1))


def cpsi_central_triple_fast(N, k, a, b):
    """pair-CPsi of |psi(a,b)> = a|D_{k-1}> + b|D_k> + a|D_{k+1}>.

    Direct evaluation in log space to handle large N without overflow.
    """
    def lC(n, r):
        return log_binom(n, r)

    def C(n, r):
        # Used only in flagging (r in [0, n] range check)
        if r < 0 or r > n:
            return 0
        return 1  # caller doesn't need value, just nonzero

    # Coefficients: c_{k-1} = a, c_k = b, c_{k+1} = a (symmetric ansatz).
    c = {k - 1: a, k: b, k + 1: a}

    def rho(s, sp):
        """ρ_AB[(s_pair), (s'_pair)] for s, s' ∈ {0, 1, 2}.

        Computes: Σ_n c_n c_m · C(N-2, n-s) / sqrt(C(N,n) C(N,m))
        with m = n + (sp - s), constraint n-s ∈ [0, N-2].
        Uses log-space binomials for large-N stability.
        """
        result = 0.0
        for n in c:
            m = n + (sp - s)
            if m not in c:
                continue
            res_pop = n - s
            if res_pop < 0 or res_pop > N - 2:
                continue
            log_factor = lC(N - 2, res_pop) - 0.5 * (lC(N, n) + lC(N, m))
            result += c[n] * c[m] * math.exp(log_factor)
        return result

    d0 = rho(0, 0)
    d1 = rho(1, 1)
    d2 = rho(2, 2)
    e_01 = rho(0, 1)   # ρ[|00>, |01>] = ρ[|00>, |10>]
    e_02 = rho(0, 2)   # ρ[|00>, |11>]
    e_12 = rho(1, 2)   # ρ[|01>, |11>] = ρ[|10>, |11>]   (≠ e_01 in general)

    # 4x4 ρ_AB in basis |00>, |01>, |10>, |11>:
    #   [d0,    e01,   e01,   e02]
    #   [e01,   d1,    d1,    e12]
    #   [e01,   d1,    d1,    e12]
    #   [e02,   e12,   e12,   d2]
    # Tr(ρ²) = sum of all |ρ_ij|² (real, Hermitian)
    #        = d0² + 2·d1² + d2² (diagonal)
    #        + 4·e01² + 2·e02² + 2·d1² + 4·e12²    (off-diagonal entries:
    #          e01 at (0,1)(1,0)(0,2)(2,0); e02 at (0,3)(3,0); d1 at (1,2)(2,1);
    #          e12 at (1,3)(3,1)(2,3)(3,2))
    trace_rho2 = (d0 ** 2 + 4 * d1 ** 2 + d2 ** 2
                  + 4 * e_01 ** 2 + 2 * e_02 ** 2 + 4 * e_12 ** 2)
    # L1_off = sum of |ρ_ij| for i ≠ j
    L1_off = (4 * abs(e_01) + 2 * abs(e_02) + 2 * abs(d1)
              + 4 * abs(e_12))
    return trace_rho2 * L1_off / 3.0


def central_triple_optimum_fast(N, n_seeds=20):
    k = (N + 1) // 2
    if k < 1 or k > N - 1:
        return None

    def neg(t):
        a = math.sin(t) / math.sqrt(2)
        b = math.cos(t)
        return -cpsi_central_triple_fast(N, k, a, b)

    best = (1.0, None)
    rng = np.random.default_rng(42)
    for seed in range(n_seeds):
        t0 = math.pi * rng.uniform()
        res = minimize(lambda x: neg(x[0]), [t0], method="Nelder-Mead",
                        options={"xatol": 1e-14, "fatol": 1e-16})
        if res.fun < best[0]:
            best = (res.fun, res.x[0])
    t_opt = best[1]
    a_opt = math.sin(t_opt) / math.sqrt(2)
    b_opt = math.cos(t_opt)
    if b_opt < 0:
        a_opt, b_opt = -a_opt, -b_opt
    if a_opt < 0:
        a_opt = -a_opt
    cpsi_opt = -best[0]
    return {"N": N, "k": k, "a": a_opt, "b": b_opt, "cpsi": cpsi_opt}


def main():
    print("=" * 80)
    print("Central-Dicke-triple asymptotics (fast binomial pair-CPsi)")
    print("=" * 80)

    Ns = [3, 4, 6, 8, 12, 16, 20, 30, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
    print()
    print(f"  {'N':>5} {'a':>14} {'b':>14} {'a²/b²':>12} {'cpsi':>14}")
    print("  " + "-" * 70)
    data = []
    for N in Ns:
        r = central_triple_optimum_fast(N, n_seeds=15)
        if r is None:
            continue
        data.append(r)
        print(f"  {r['N']:>5d} {r['a']:>14.10f} {r['b']:>14.10f} "
              f"{r['a']**2/r['b']**2:>12.6f} {r['cpsi']:>14.10f}")
    print()

    # Asymptotic
    a_inf = data[-1]['a']
    b_inf = data[-1]['b']
    cp_inf = data[-1]['cpsi']
    print(f"  Asymptotic at N={data[-1]['N']}:")
    print(f"    a = {a_inf:.10f}   (a² = {a_inf**2:.10f})")
    print(f"    b = {b_inf:.10f}   (b² = {b_inf**2:.10f})")
    print(f"    cpsi = {cp_inf:.10f}")
    print(f"    a²/b² = {a_inf**2/b_inf**2:.10f}")
    print(f"    a/b   = {a_inf/b_inf:.10f}")
    print()

    # Look at trends
    Ns_arr = np.array([r['N'] for r in data])
    cpsi_arr = np.array([r['cpsi'] for r in data])
    a_arr = np.array([r['a'] for r in data])
    b_arr = np.array([r['b'] for r in data])

    # Test if cpsi → 0 with various rates
    print("  Asymptotic decay tests (large N):")
    print(f"    {'N':>5} {'cpsi':>10} {'cpsi·log(N)':>13} {'cpsi·N^0.3':>12} "
          f"{'cpsi·sqrt(log)':>15}")
    for r, cp in zip(data, cpsi_arr):
        N = r['N']
        l = math.log(N)
        print(f"    {N:>5d} {cp:>10.6f} {cp*l:>13.6f} {cp*N**0.3:>12.6f} "
              f"{cp*math.sqrt(l):>15.6f}")
    print()

    # 1/cpsi vs ln(N): linear?
    inv_cp = 1.0 / cpsi_arr
    log_N = np.log(Ns_arr)
    coeffs = np.polyfit(log_N[-6:], inv_cp[-6:], 1)
    print(f"  Linear fit at large N: 1/cpsi ≈ {coeffs[0]:.6f} · ln(N) + "
          f"{coeffs[1]:.6f}")
    pred = np.polyval(coeffs, log_N)
    err = np.max(np.abs(pred - inv_cp))
    print(f"  Max abs error in 1/cpsi over all N: {err:.4e}")
    print()

    # Check if cpsi → 0 (it does)
    print(f"  cpsi at N={Ns[-1]} = {cpsi_arr[-1]:.6e}")
    print()


if __name__ == "__main__":
    main()
