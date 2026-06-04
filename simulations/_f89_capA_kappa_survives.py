#!/usr/bin/env python3
"""F89 capA probe 8: rigorize PIECE A -- the WHOLE kappa^2 survives the c-reduction (not just 2-adically).

Cleaner rigorous statement than '2-content = 2 v2(kappa)':
   the c-reduced denominator equals  kappa^2 / gcd-stuff  with the (kappa+1)^2 FULLY cancelled
   and kappa^2 surviving INTACT (odd part and 2-part both).

We prove it by exhibiting the reduction as an operation over Z[1/(kappa+1)] -- i.e. show that
after clearing ONLY the (kappa+1)^2 and the lone 2 (the m-tied content), the remaining object
   q(c) := (kappa^2 / 2) * p_n(c) = G(c)/(kappa+1)^2
reduces mod Phi_c to a polynomial with denominator COPRIME to kappa (so kappa^2 cannot be
diminished), while p_n itself reduces to (2/kappa^2)*[that], giving denominator exactly kappa^2
times a unit at 2 (the 2-content 2 v2(kappa) is then immediate, and so is the odd part: odd(kappa)^2).

Checks, full (not just v2):
  (a) full c-reduced denominator of p_n(c)               -> call it Den_c
  (b) Den_c / kappa^2   is an integer? (kappa^2 | Den_c)  -> kappa^2 divides the denom
  (c) Den_c == kappa^2 * (odd unit)?  i.e. Den_c / kappa^2 has the SAME odd part as ... we just
      report Den_c factorised and Den_c/kappa^2.
  (d) reduce q(c) = G/(kappa+1)^2 mod Phi_c; show its denominator is coprime to kappa.
  (e) reduce r(c) = G mod Phi_c; show (kappa+1)^2 divides its denominator's ... i.e. the (kappa+1)
      cancellation is what clears, leaving the kappa-coprime core.
"""
from __future__ import annotations

import sys
from math import gcd
from pathlib import Path

import sympy as sp

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))

c = sp.Symbol("c")


def Atil(k):
    kap = k // 2
    return sp.expand(sum(sp.chebyshevu(j, c) * (kap - j) for j in range(k + 1)))


def Btil(k):
    kap = k // 2
    return sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (kap - j) ** 2 for j in range(k + 1)))


def G_poly(k):
    return sp.expand((1 - c**2) ** 2 * Atil(k) ** 2 * Btil(k))


def orbit_poly_c(k):
    m = k + 2
    orbit = list(range(2, k + 2, 2))
    roots = [sp.cos(sp.pi * n / m) for n in orbit]
    mps, seen = [], set()
    for r in roots:
        mp_r = sp.minimal_polynomial(r, c)
        s = str(mp_r)
        if s not in seen:
            seen.add(s)
            mps.append(mp_r)
    return sp.expand(sp.prod(mps))


def denom_lcm(poly_expr):
    poly = sp.Poly(poly_expr, c, domain="QQ")
    d = sp.Integer(1)
    for coef in poly.all_coeffs():
        d = sp.lcm(d, sp.denom(sp.Rational(coef)))
    return int(d)


def main():
    print("=" * 110)
    print("F89 capA probe 8: PIECE A rigor -- whole kappa^2 survives; (kappa+1)^2 cancels; kappa-coprime core")
    print("=" * 110)
    print()
    hdr = (f"  {'k':>3} {'kap':>4} {'kap^2':>6} | {'Den_c':>10} {'kap2|Den':>8} {'Den/kap2':>9} "
           f"{'gcd(Den/kap2,kap)':>17} | {'denom(q)':>9} {'gcd(dq,kap)':>11}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))

    reps = [6, 10, 12, 20, 8, 24, 16, 18, 36, 40]
    for k in reps:
        kap = k // 2
        G = G_poly(k)
        op = orbit_poly_c(k)

        p_c = sp.Rational(2) / ((kap + 1) ** 2 * kap**2) * G
        red = sp.rem(sp.Poly(p_c, c, domain="QQ"), sp.Poly(op, c, domain="QQ")).as_expr()
        Den_c = denom_lcm(red)

        kap2_divides = (Den_c % (kap**2) == 0)
        Den_over = Den_c // (kap**2) if kap2_divides else -1
        g1 = gcd(Den_over, kap) if kap2_divides else -1

        # q = G/(kap+1)^2 reduced; denom should be coprime to kappa
        q = sp.Rational(1, (kap + 1) ** 2) * G
        redq = sp.rem(sp.Poly(q, c, domain="QQ"), sp.Poly(op, c, domain="QQ")).as_expr()
        dq = denom_lcm(redq)
        g2 = gcd(dq, kap)

        print(f"  {k:>3} {kap:>4} {kap**2:>6} | {Den_c:>10} {str(kap2_divides):>8} {Den_over:>9} "
              f"{g1:>17} | {dq:>9} {g2:>11}")

    print()
    print("  kap2|Den=True : kappa^2 divides the c-reduced denominator (the WHOLE kappa^2 survives).")
    print("  Den/kap2 with gcd(.,kappa)=1 : the surviving denominator is EXACTLY kappa^2 times a")
    print("     kappa-coprime cofactor, so v2(Den_c)=2 v2(kappa) AND odd-part contains odd(kappa)^2.")
    print("  denom(q) coprime to kappa : reducing G/(kappa+1)^2 (the m-tied part) yields a")
    print("     kappa-coprime denominator => kappa^2 enters ONLY through the explicit 1/kappa^2,")
    print("     never cancelled by the reduction. This is the mechanism, fully localized.")


if __name__ == "__main__":
    main()
