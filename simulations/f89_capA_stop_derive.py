#!/usr/bin/env python3
"""F89 capA probe 13: derive s_top analytically -- split L into Gtilde-leading and the 1/kappa^2.

After the even-k pull-out: p_n(c) = 8/(m^2 kappa^2) * Gtil(c),  Gtil = (1-c^2)^2 Atil^2 Btil,
INTEGER coeffs, m = 2(kappa+1).  So
   p_n(c) = 2/((kappa+1)^2 kappa^2) * Gtil(c).
Reduction is LINEAR in p_n, so
   L = [c^{FA-1}]( p_n mod Phi_c ) = 2/((kappa+1)^2 kappa^2) * Ltil,
   Ltil := [c^{FA-1}]( Gtil mod Phi_c )   -- an algebraic number with den | (orbit content).
Hence
   val2(L) = 1 - 2 v2(kappa+1) - 2 v2(kappa) + val2(Ltil).               (*)

For EVEN k, exactly one of kappa, kappa+1 is even:
   - kappa odd  (k=2*odd, v2k=1): v2(kappa)=0, so val2(L) = 1 - 2 v2(kappa+1) + val2(Ltil).
   - kappa even (4|k, v2k>=2):    v2(kappa+1)=0, so val2(L) = 1 - 2 v2(kappa) + val2(Ltil).

So s_top = val2(L). The claim s_top = (FA-1) + offset(v2k) becomes a claim about val2(Ltil),
the reduced leading coeff of the INTEGER polynomial Gtil (no 1/kappa^2, no 1/(kappa+1)^2).
This probe computes Ltil and val2(Ltil) and checks the identity (*), then reads what val2(Ltil)
must be:  val2(Ltil) = s_top - 1 + 2 v2(kappa) + 2 v2(kappa+1).

We tabulate val2(Ltil) and compare to (FA-1) + correction, isolating whether Ltil's 2-content
is the 'generic degree-growth' part (FA-1 + const) with the v2-dependence living ENTIRELY in the
explicit 1/(kappa^2 (kappa+1)^2) prefactor. THAT would be the clean analytic statement:
   v2(Ltil) = (FA-1) + 1   (a pure degree-growth constant, v2-INDEPENDENT)  ??? -- check.
If yes: s_top = (FA-1)+1 + 1 - 2v2(kap) - 2v2(kap+1) ... must reconcile to offset. Check numerically.
"""
from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))

c = sp.Symbol("c")


def v2_int(n):
    n = abs(int(n)); v = 0
    while n % 2 == 0:
        n //= 2; v += 1
    return v


def val2(r):
    r = sp.Rational(r)
    if r == 0:
        return None
    num, den = sp.fraction(r)
    return v2_int(num) - v2_int(den)


def Gtil(k):
    kap = k // 2
    At = sp.expand(sum(sp.chebyshevu(j, c) * (kap - j) for j in range(k + 1)))
    Bt = sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (kap - j) ** 2 for j in range(k + 1)))
    return sp.expand((1 - c**2) ** 2 * At**2 * Bt)


def Gfull(k):
    """non-pulled integer poly G = (1-c^2)^2 A^2 B (odd k uses this; even k Gtil*4 = ...)."""
    A = sp.expand(sum(sp.chebyshevu(j, c) * (k - 2 * j) for j in range(k + 1)))
    B = sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (k - 2 * j) ** 2 for j in range(k + 1)))
    return sp.expand((1 - c**2) ** 2 * A**2 * B)


def orbit_poly_c(k):
    m = k + 2
    orbit = list(range(2, k + 2, 2))
    roots = [sp.cos(sp.pi * n / m) for n in orbit]
    mps, seen = [], set()
    for r in roots:
        mp_r = sp.minimal_polynomial(r, c)
        s = str(mp_r)
        if s not in seen:
            seen.add(s); mps.append(mp_r)
    return sp.expand(sp.prod(mps))


def lead_reduced(poly_expr, k):
    FA = (k + 1) // 2
    P = sp.rem(sp.Poly(poly_expr, c, domain="QQ"),
               sp.Poly(orbit_poly_c(k), c, domain="QQ")).as_expr()
    return sp.Rational(sp.Poly(P, c, domain="QQ").coeff_monomial(c**(FA - 1)))


def offset(v):
    return 2 - v - max(0, v - 2)


def main():
    print("=" * 120)
    print("F89 capA probe 13: s_top via integer Gtil; v2(Ltil) -- is it v2-independent degree-growth?")
    print("=" * 120)
    print()
    hdr = (f"  {'k':>3} {'v2':>2} {'FA-1':>4} {'kap':>4} {'v2kap':>5} {'v2(kap+1)':>9} | "
           f"{'s_top':>5} {'v2(Ltil)':>8} {'v2(Lfull)':>9} | "
           f"{'identity(*)':>11} | {'v2(Ltil)-(FA-1)':>15}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))

    reps = [5, 7, 9, 11, 6, 10, 14, 22, 12, 20, 36, 8, 24, 40, 16, 48, 32, 64]
    for k in reps:
        FA = (k + 1) // 2
        v2k = v2_int(k)
        kap = k // 2

        # s_top from FULL p_n (ground truth path)
        m = k + 2
        A = sp.expand(sum(sp.chebyshevu(j, c) * (k - 2 * j) for j in range(k + 1)))
        B = sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (k - 2 * j) ** 2 for j in range(k + 1)))
        p_c = sp.expand(sp.Rational(2) / (m**2 * k**2) * (1 - c**2) ** 2 * A**2 * B)
        L = lead_reduced(p_c, k)
        s_top = val2(L)

        if k % 2 == 0:
            Lt = lead_reduced(Gtil(k), k)
            v2Lt = val2(Lt)
            v2kap = v2_int(kap)
            v2kap1 = v2_int(kap + 1)
            # identity (*): val2(L) =? 1 - 2 v2(kap+1) - 2 v2(kap) + v2(Ltil)
            rhs = 1 - 2 * v2kap1 - 2 * v2kap + v2Lt
            id_ok = (s_top == rhs)
            v2Lt_minus = v2Lt - (FA - 1)
        else:
            # odd k: no pull-out; use Gfull, p_n = 2/(m^2 k^2) Gfull
            Lt = lead_reduced(Gfull(k), k)
            v2Lt = val2(Lt)
            v2kap = -1; v2kap1 = -1
            # identity: val2(L) = 1 - 2 v2(m) - 2 v2(k) + v2(Lfull) ; k odd => v2(k)=0,v2(m)=v2(k+2)=0
            rhs = 1 - 2 * v2_int(m) - 2 * v2_int(k) + v2Lt
            id_ok = (s_top == rhs)
            v2Lt_minus = v2Lt - (FA - 1)

        print(f"  {k:>3} {v2k:>2} {FA-1:>4} {kap:>4} {v2kap:>5} {v2kap1:>9} | "
              f"{s_top:>5} {v2Lt:>8} {'-':>9} | {str(id_ok):>11} | {v2Lt_minus:>15}")

    print()
    print("  identity(*) True everywhere => s_top = 1 - 2 v2(kap) - 2 v2(kap+1) + v2(Ltil)  [even k].")
    print("  Watch column v2(Ltil)-(FA-1): if CONSTANT across v2 (e.g. always 1 or 2), the whole")
    print("  v2-dependence of s_top is carried by the explicit 2/((kap+1)^2 kap^2) prefactor and the")
    print("  reduction of the INTEGER Gtil contributes only generic degree-growth -> cap PROVEN.")
    print("  If v2(Ltil)-(FA-1) varies with v2, the cap hides partly inside Gtil's reduction (localize).")


if __name__ == "__main__":
    main()
