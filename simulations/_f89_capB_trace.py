#!/usr/bin/env python3
"""F89 capB: Rm_top as a trace over the orbit; 2-adic via the different (2026-06-04).

Rm_top = [u^{FA-1}](Nint mod Phi_u) = sum_{j=1}^{FA} Nint(u_j)/Phi_u'(u_j)   (top divided diff),
u_j = 2cos(2 pi j/m).  This is a Q-rational integer (Rm integer poly), and equals
   sum over distinct minimal-poly factors P of Phi_u of  Tr_{K_P/Q}( Nint(alpha)/P'(alpha) * [Phi_u/P](alpha)^{-1} ... )
-- but since Phi_u = prod of DISTINCT irreducible factors (squarefree: the u_j are distinct),
Phi_u'(u_j) = P'(u_j) * prod_{Q != P}(Q(u_j)) for the factor P with P(u_j)=0.

We verify Rm_top = sum_j Nint(u_j)/Phi_u'(u_j) numerically/exactly, then study v2.

KEY 2-adic question: what is v2 of the orbit codifferent contribution?  We compute, for each
distinct factor P (a real-cyclotomic minimal poly of 2cos(2 pi/d) for d|m-ish), the 2-adic
valuation data: v2(disc(P)), v2(P(0)), and whether 2 ramifies.  The trace's 2-content is
bounded by -v2(different) = -v2 of the codifferent.  We tabulate to see the cap emerge.

This probe is EXPLORATORY: it tests whether v2(Rm_top) tracks the orbit different's 2-part
or whether (as Angle B suggests) it is a numerator-specific cancellation that the different
only bounds.  Decisive for whether C2/C3 is 'algebraic-number-theory clean' or 'needs the
explicit Chebyshev numerator'.
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
u = sp.Symbol("u")


def v2(n) -> int:
    n = abs(int(n))
    if n == 0:
        return 10**9
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def build_Nint_s2(k):
    A = sp.expand(sum(sp.chebyshevu(j, c) * (k - 2 * j) for j in range(k + 1)))
    B = sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (k - 2 * j) ** 2 for j in range(k + 1)))
    Nc = sp.expand((1 - c**2) ** 2 * A**2 * B)
    Nu = sp.expand(Nc.subs(c, u / 2))
    poly = sp.Poly(Nu, u)
    d = sp.Integer(1)
    for co in poly.all_coeffs():
        d = sp.lcm(d, sp.denom(sp.Rational(co)))
    s2 = v2(int(d))
    return sp.expand(Nu * (2 ** s2)), s2


def distinct_factors(k):
    m = k + 2
    FA = (k + 1) // 2
    facs, seen = [], set()
    for j in range(1, FA + 1):
        r = 2 * sp.cos(2 * sp.pi * j / m)
        mp = sp.minimal_polynomial(r, u)
        s = str(mp)
        if s not in seen:
            seen.add(s)
            facs.append(mp)
    return facs


def main():
    print("=" * 110)
    print("F89 capB: orbit different 2-content vs v2(Rm_top) -- is the cap NT-clean?")
    print("=" * 110)
    print()
    reps = [5, 7, 9, 11, 6, 10, 14, 18, 22, 12, 20, 28, 8, 24, 16, 32]
    print(f"  {'k':>3} {'v2k':>3} {'m':>3} {'#factors':>8} {'v2(disc Phi)':>12} "
          f"{'sum v2(P(0))':>12} {'v2(Rm_top)':>10}")
    print("  " + "-" * 70)
    for k in reps:
        m = k + 2
        facs = distinct_factors(k)
        Phi = sp.expand(sp.prod(facs))
        disc = sp.discriminant(sp.Poly(Phi, u).as_expr(), u)
        v2disc = v2(int(disc)) if disc != 0 else -1
        # sum of v2 of each factor's constant term
        sum_v2P0 = 0
        for P in facs:
            P0 = sp.Poly(P, u).all_coeffs()[-1]
            sum_v2P0 += (v2(int(P0)) if P0 != 0 else 10**6)
        # Rm_top
        Nint, s2 = build_Nint_s2(k)
        Rm = sp.rem(sp.Poly(Nint, u), sp.Poly(Phi, u))
        Rt = int(sp.Poly(Rm.as_expr(), u).all_coeffs()[0])
        print(f"  {k:>3} {v2(k):>3} {m:>3} {len(facs):>8} {v2disc:>12} "
              f"{str(sum_v2P0) if sum_v2P0<10**6 else 'inf':>12} {v2(Rt):>10}")
    print()
    print("  Reading: if v2(Rm_top) << v2(disc Phi) (the different bound), the cap is NOT the")
    print("  full different -- it is a numerator-specific cancellation, consistent with Angle B")
    print("  (Washington: D_k not a Galois invariant).  The different only UPPER-bounds the denom.")
    print()
    # verify Rm_top = sum_j Nint(u_j)/Phi'(u_j) exactly for a couple small k (numeric high-prec)
    print("-" * 110)
    print("  numeric check Rm_top = sum_j Nint(u_j)/Phi_u'(u_j):")
    import mpmath as mp
    mp.mp.dps = 60
    for k in [5, 6, 8, 12]:
        m = k + 2
        FA = (k + 1) // 2
        Nint, s2 = build_Nint_s2(k)
        facs = distinct_factors(k)
        Phi = sp.expand(sp.prod(facs))
        dPhi = sp.diff(Phi, u)
        Nfun = sp.lambdify(u, Nint, "mpmath")
        dPhifun = sp.lambdify(u, dPhi, "mpmath")
        tot = mp.mpf(0)
        for j in range(1, FA + 1):
            uj = 2 * mp.cos(2 * mp.pi * j / m)
            tot += Nfun(uj) / dPhifun(uj)
        Rm = sp.rem(sp.Poly(Nint, u), sp.Poly(Phi, u))
        Rt = int(sp.Poly(Rm.as_expr(), u).all_coeffs()[0])
        print(f"    k={k}: sum={mp.nstr(tot, 20)}  Rm_top(exact)={Rt}  match={abs(tot-Rt)<1e-20}")
    print()


if __name__ == "__main__":
    main()
