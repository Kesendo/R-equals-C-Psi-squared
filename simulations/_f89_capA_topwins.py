#!/usr/bin/env python3
"""F89 capA probe 16: make 'top-degree wins' a Newton-polygon edge statement.

We use v2(D_k) = max_d (2 d - s_d), s_d = val2([c^d] P_c). Top-wins means d=FA-1 attains the max:
   for all d < FA-1:   2 d - s_d  <=  2(FA-1) - s_top
   <=>  s_top - s_d  <=  2 (FA-1 - d).                                  (TW)
i.e. going DOWN from the top by t steps, the signed valuation s drops by at most 2t:
   s_{FA-1} - s_{FA-1-t} <= 2t   for all t>=1.

Equivalently the upper boundary of the points (d, s_d) (the lower convex hull in the
2-adic-Newton sense for 'max 2d - s_d') has its supporting line of slope 2 touching at the top.

This probe:
  (1) prints, for each k, the per-step drop  delta_t := s_{FA-1} - s_{FA-1-t}  for t=1..min(6,FA-1),
      and the bound 2t, flagging any violation (none expected).
  (2) reports the MARGIN  min_t (2t - delta_t)  (how far from tight) -- if always >=0, TW holds;
      if it is frequently 0, TW is tight (interesting), if always >0 strictly, TW has slack.
  (3) reports the denominator v2 of each near-top coeff (must divide kappa^2 => den-v2 <= 2 v2(kap)),
      bounding s_d from below: s_d >= -2 v2(kappa). Combined with s_top closed form this gives a
      structural (non-fitted) proof of TW IF s_d >= s_top - 2(FA-1-d) follows from den<=kappa^2 and
      num>=0. We check num(>=0) i.e. s_d>= -2v2(kap) and whether that ALONE forces TW.
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


def p_n_in_c_expr(k):
    m = k + 2
    A = sp.expand(sum(sp.chebyshevu(j, c) * (k - 2 * j) for j in range(k + 1)))
    B = sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (k - 2 * j) ** 2 for j in range(k + 1)))
    return sp.expand(sp.Rational(2) / (m**2 * k**2) * (1 - c**2) ** 2 * A**2 * B)


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


def main():
    print("=" * 110)
    print("F89 capA probe 16: top-wins as Newton edge  s_{FA-1}-s_{FA-1-t} <= 2t")
    print("=" * 110)
    print()
    reps = [7, 9, 11, 10, 14, 22, 12, 20, 8, 24, 16, 32]
    global_min_margin = 10**9
    any_violation = False
    for k in reps:
        FA = (k + 1) // 2
        deg_red = FA - 1
        v2kap = v2_int(k // 2) if k % 2 == 0 else 0

        P_c = sp.rem(sp.Poly(p_n_in_c_expr(k), c, domain="QQ"),
                     sp.Poly(orbit_poly_c(k), c, domain="QQ")).as_expr()
        Pc = sp.Poly(P_c, c, domain="QQ")
        s = {d: val2(sp.Rational(Pc.coeff_monomial(c**d))) for d in range(deg_red + 1)}
        s_top = s[deg_red]

        T = min(6, deg_red)
        drops = []
        margins = []
        for t in range(1, T + 1):
            sd = s[deg_red - t]
            if sd is None:
                drops.append("inf"); continue
            delta = s_top - sd
            drops.append(delta)
            margins.append(2 * t - delta)
        mm = min(margins) if margins else 99
        if any(isinstance(d, int) and d > 2 * (i + 1) for i, d in enumerate(drops)):
            any_violation = True
        global_min_margin = min(global_min_margin, mm)

        # lower bound on s_d from den | kappa^2:  s_d >= -2 v2(kappa); does num>=0 hold (s_d>=-2v2kap)?
        min_s = min(v for v in s.values() if v is not None)
        den_bound_ok = (min_s >= -2 * v2kap)

        print(f"k={k:>3} FA-1={deg_red} v2kap={v2kap} s_top={s_top} | drops(t=1..{T})={drops} | "
              f"2t bound -> min margin={mm} | min s_d={min_s} >=-2v2kap({-2*v2kap}):{den_bound_ok}")
    print()
    print(f"any_violation={any_violation}  (False => TW holds on all tested k)")
    print(f"global min margin = {global_min_margin}  (>=0 => top is a supporting vertex of slope-2 line)")
    print()
    print("Note: TW is here a VERIFIED Newton-edge inequality. A fully rigorous TW needs a bound")
    print("s_d >= s_top - 2(FA-1-d). The den|kappa^2 bound gives s_d >= -2 v2(kap); that ALONE does")
    print("not yet force TW (it bounds s_d below by a constant, not by the slope-2 line). The slope-2")
    print("edge is the precise residual analytic statement (left localized, see report).")


if __name__ == "__main__":
    main()
