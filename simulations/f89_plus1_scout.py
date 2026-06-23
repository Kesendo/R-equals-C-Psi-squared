#!/usr/bin/env python3
"""F89 '+1' lemma scout (2026-06-04): reduce val2(L~) to an explicit trig-sum trace.

The last open input in the cap proof is val2(L~) = (FA-1)+1 for odd k and 4|k, where
L~ = leading coeff of (G mod Phi_c), G = (1-c^2)^2 A(c)^2 B(c) the integer numerator,
Phi_c the orbit nodal polynomial (roots c_j = cos(2 pi j / m), j=1..FA, m=k+2).

This scout establishes (numerically, high precision) the analytic reformulation that
turns the '+1' into a concrete trigonometric trace whose 2-adic carrier is explicit:

  (H1) at a node theta = 2 pi j / m, the sin-denominators of U_j(cos theta)=sin((j+1)theta)/sin theta
       cancel against (1-c^2)^2 = sin^4 theta, giving  G(c_j) = S1(theta_j)^2 * S2(theta_j),
       S1(theta) = sum_{i=1}^{m-1} (m-2i) sin(i theta),  S2(theta) = sum_{i=1}^{m-1} (m-2i)^2 sin^2(i theta).
  (H2) S1(theta_j) = m * cot(pi j / m)   (closed form; Sum i sin(i th) = -(m/2) cot).
  (H3) L~ = sum_j G(c_j) / Phi'(c_j)      (leading coeff of the Lagrange interpolant = codifferent trace).
  (H4) for odd m:  T_m(c) - 1 = 2^{m-1} (c-1) Phi(c)^2,  hence
       Phi'(c_j)^2 = -m^2 / (2^m (c_j - 1)(1 - c_j^2)).   <-- the 2^m is the 2-adic carrier.

If all four hold, the '+1' is the statement that the trace sum_j m^2 cot^2(pi j/m) S2(theta_j)/Phi'(c_j),
an explicit algebraic number, has 2-adic valuation exactly (FA-1)+1. Odd-k case here; 4|k analogous.
"""
from __future__ import annotations

import sys
from mpmath import mp, mpf, sin, cos, cot, pi, sqrt, fabs, nint

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

mp.dps = 60


def cheb_u(j, x):
    """U_j(x) by recurrence."""
    um1, u0 = mpf(0), mpf(1)
    if j == 0:
        return u0
    for _ in range(j):
        um1, u0 = u0, 2 * x * u0 - um1
    return u0


def v2_int(n: int) -> int:
    n = abs(int(n))
    if n == 0:
        return 0
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def main():
    print("=" * 100)
    print("F89 '+1' scout: G(node) = S1^2 S2, S1 = m cot, Phi'^2 = -m^2/(2^m (c-1)(1-c^2)), L~ = trace")
    print("=" * 100)
    print()
    print(f"  {'k':>3} {'m':>3} {'FA':>3} | {'H1 maxerr':>10} {'H2 maxerr':>10} {'H4 maxerr':>10} | "
          f"{'L~ (trace)':>16} {'v2(L~)':>7} {'(FA-1)+1':>9} {'ok':>4}")
    for k in (5, 7, 9, 11, 13, 15):
        m = k + 2
        FA = (m - 1) // 2
        nodes = [cos(2 * pi * j / m) for j in range(1, FA + 1)]
        h1err = h2err = h4err = mpf(0)
        Ltr = mpf(0)
        for idx, j in enumerate(range(1, FA + 1)):
            c = nodes[idx]
            th = 2 * pi * j / m
            # G(c) from the integer-poly side (independent of the trig sums)
            A = sum(cheb_u(jj, c) * (k - 2 * jj) for jj in range(k + 1))
            B = sum(cheb_u(jj, c) ** 2 * (k - 2 * jj) ** 2 for jj in range(k + 1))
            Gc = (1 - c * c) ** 2 * A * A * B
            # trig-sum side
            S1 = sum((m - 2 * i) * sin(i * th) for i in range(1, m))
            S2 = sum((m - 2 * i) ** 2 * sin(i * th) ** 2 for i in range(1, m))
            h1err = max(h1err, fabs(Gc - S1 * S1 * S2))
            h2err = max(h2err, fabs(S1 - m * cot(pi * j / m)))
            # Phi'(c_j) = prod_{l != j}(c_j - c_l)
            Phip = mpf(1)
            for l in range(FA):
                if l != idx:
                    Phip *= (c - nodes[l])
            rhs = -mpf(m) ** 2 / (mpf(2) ** m * (c - 1) * (1 - c * c))
            h4err = max(h4err, fabs(Phip * Phip - rhs))
            Ltr += Gc / Phip
        Ln = int(nint(Ltr))
        ok = (v2_int(Ln) == (FA - 1) + 1) and (fabs(Ltr - Ln) < mpf(10) ** (-20))
        print(f"  {k:>3} {m:>3} {FA:>3} | {mp.nstr(h1err,2):>10} {mp.nstr(h2err,2):>10} {mp.nstr(h4err,2):>10} | "
              f"{Ln:>16} {v2_int(Ln):>7} {(FA-1)+1:>9} {('OK' if ok else 'FAIL'):>4}")
    print()
    print("  H1: G(c_j) == S1^2 S2   H2: S1 == m cot(pi j/m)   H4: Phi'^2 == -m^2/(2^m (c-1)(1-c^2))")
    print("  L~ (trace) integer, v2 = (FA-1)+1: the '+1' lemma in the trig-trace form.")
    print()
    print("  => If errors ~1e-50, the '+1' = [v2 of  sum_j m^2 cot^2(pi j/m) S2(theta_j) / Phi'(c_j)] - (FA-1).")
    print("     Phi'(c_j) carries 2^{-m/2} (from H4); the trace is an integer, so the 2-powers recombine")
    print("     across the Galois orbit to leave exactly +1 above the degree-growth (FA-1). That recombination")
    print("     is the explicit target: a 2-adic valuation of one trig sum, no polynomial reduction left.")


if __name__ == "__main__":
    main()
