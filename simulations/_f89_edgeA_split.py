#!/usr/bin/env python3
"""F89 edge A, step 5: the clean split  (A) slope-1 for d>=1 all k  +  (B) the d=0 corner.

Step 4 finding: SLOPE1  val2(R_d) >= d + delta  holds for EVERY even k and every d, and for
odd k it holds for all d>=1 but FAILS at d=0 where val2(R_0) = -1 = delta - 2 exactly (delta=1).
The d=0 failure is harmless for dominance because TARGET-R only needs the slope-(-2) floor.

This script makes that precise and bit-exact over a WIDE k range (5..40):

  (A)  for all d in 1..FA-1:   val2(R_d) >= d + delta(k).          [SLOPE1>=1]
  (B0) at d=0:                 val2(R_0) = -1 (odd k),  >= delta (even k);
       and in BOTH cases val2(R_0) >= val2(R_top) - 2(FA-1)  [TARGET-R at d=0].

Then assembles the implication chain to TARGET-R for all d, hence top dominance:
  - d>=1: val2(R_d) >= d+delta = val2(R_top) - (FA-1-d) >= val2(R_top) - 2(FA-1-d).  [margin >=0]
  - d=0:  checked directly in (B0).
We print, per k: the min slack of (A) over d>=1, the d=0 values, and a single PASS/FAIL.
Also dumps R_0 itself (odd part) to expose that the odd-k constant term is 2*(odd)/(... ) with
exactly one 2 in the denominator (val2 = -1), the universal corner value.
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


def v2_int(n) -> int:
    n = abs(int(n))
    if n == 0:
        return 10**9
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def signed_val2(q):
    q = sp.Rational(q)
    if q == 0:
        return None
    num, den = sp.fraction(q)
    return v2_int(int(num)) - v2_int(int(den))


def reduced_R_c(k):
    if k % 2 == 1:
        A = sp.expand(sum(sp.chebyshevu(j, c) * (k - 2 * j) for j in range(k + 1)))
        B = sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (k - 2 * j) ** 2 for j in range(k + 1)))
    else:
        kap = k // 2
        A = sp.expand(sum(sp.chebyshevu(j, c) * (kap - j) for j in range(k + 1)))
        B = sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (kap - j) ** 2 for j in range(k + 1)))
    G = sp.expand((1 - c**2) ** 2 * A**2 * B)
    m = k + 2
    FA = (k + 1) // 2
    roots = [sp.cos(2 * sp.pi * j / m) for j in range(1, FA + 1)]
    mps, seen = [], set()
    for r in roots:
        mp_r = sp.minimal_polynomial(r, c)
        s = str(mp_r)
        if s not in seen:
            seen.add(s)
            mps.append(mp_r)
    Phi = sp.expand(sp.prod(mps))
    Phi_monic = sp.expand(Phi / sp.Rational(sp.Poly(Phi, c).LC()))
    R = sp.rem(sp.Poly(G, c, domain="QQ"), sp.Poly(Phi_monic, c, domain="QQ"))
    deg = R.degree()
    return [sp.Rational(R.coeff_monomial(c**d)) for d in range(deg + 1)], FA, m


def delta_closed(k):
    if k % 2 == 1:
        return 1
    if k % 4 == 0:
        return 1
    return 2 * v2_int(k // 2 + 1)


def cls_of(k):
    if k % 2 == 1:
        return "odd"
    if k % 4 == 0:
        return "4|k"
    return "2mod4"


def main():
    print("=" * 112)
    print("F89 edge A: SPLIT proof skeleton  (A) val2(R_d)>=d+delta for d>=1 all k  +  (B0) d=0 corner")
    print("=" * 112)
    print()
    print(f"  {'k':>3} {'cls':>6} {'FA':>3} {'delta':>5} {'minSlack(d>=1)':>14} "
          f"{'val2(R0)':>9} {'R0 odd-part':>20} {'TARGET-R@0_slack':>16} {'A_ok':>5} {'all_ok':>6}")
    print("  " + "-" * 108)
    allA = True
    all0 = True
    allDom = True
    for k in range(5, 41):
        R_d, FA, m = reduced_R_c(k)
        vR = [signed_val2(r) for r in R_d]
        top = FA - 1
        delta = delta_closed(k)
        vRtop = vR[top]
        # (A) d>=1 slope-1
        slacksA = [vR[d] - (d + delta) for d in range(1, FA) if vR[d] is not None]
        minA = min(slacksA) if slacksA else None
        A_ok = (minA is not None and minA >= 0)
        allA = allA and A_ok
        # (B0) d=0
        v0 = vR[0]
        # odd part of R_0 numerator
        num0, den0 = sp.fraction(R_d[0])
        odd_num0 = int(num0) // (2 ** v2_int(int(num0))) if num0 != 0 else 0
        tr0_floor = vRtop - 2 * (top - 0)   # slope-2 floor at d=0
        slack0 = (v0 - tr0_floor) if v0 is not None else None
        b0_ok = (slack0 is not None and slack0 >= 0)
        all0 = all0 and b0_ok
        # whole dominance via the chain
        dom_ok = A_ok and b0_ok
        allDom = allDom and dom_ok
        print(f"  {k:>3} {cls_of(k):>6} {FA:>3} {delta:>5} {str(minA):>14} "
              f"{str(v0):>9} {str(odd_num0):>20} {str(slack0):>16} {str(A_ok):>5} {str(dom_ok):>6}")
    print()
    print(f"  (A)  val2(R_d) >= d + delta(k) for all d>=1, all k=5..40:  {allA}")
    print(f"  (B0) d=0 satisfies the slope-(-2) floor (TARGET-R@0), all k:  {all0}")
    print(f"  => top dominance via the split, all k=5..40:  {allDom}")
    print()
    print("  Note: for odd k, val2(R0) = -1 universally (R0 = 2*odd / (2*...): one 2 in denom).")
    print("        That is below the slope-1 line (d+delta=1) by 2, but ABOVE the slope-2 floor.")
    print("        For even k, val2(R0) >= delta (slope-1 holds at 0 too).")


if __name__ == "__main__":
    main()
