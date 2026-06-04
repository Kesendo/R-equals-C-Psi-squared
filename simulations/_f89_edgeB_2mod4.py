#!/usr/bin/env python3
"""F89 top-degree dominance, the k=2 mod 4 class (2026-06-04), probe B-11.

For k=2 mod 4 the picture differs (from B-7/B-9): the Newton polygon of R(u) is NOT a single
top-anchored segment; V=v2(R_top)=2 v2(m)+2 is LARGE, and the constant term can be deeper or
shallower than the top. We still need (DOM): for all d, v2(R_d)-d >= V-(FA-1).

This probe maps the 2mod4 structure:
  - the content 2^t | R (largest t with R == 0 mod 2^t),
  - R/2^t reduced mod 2 support,
  - whether the analogous floor argument works: does there exist a floor f and endpoints giving
    (DOM)?  We test the worst deficit directly and the rightmost hull slope (must be <=1).

We also test the candidate three-fact analogue for 2mod4:
  (G1) v2(R_top) = 2 v2(m) + 2                              [closed via L~ even-class lemma]
  (G2) v2(R_d) - d >= V - (FA-1)  for all d                 [= DOM itself]
and we look for a clean UNIFORM FLOOR f2 := min_{1<=d<=FA-1}(v2(R_d)) and whether
  v2(R_0) >= V-(FA-1)  and  f2 >= V - (FA-2)  (the two endpoints + floor) suffice.
"""
from __future__ import annotations

import sys
from pathlib import Path
from fractions import Fraction

import sympy as sp

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))

c = sp.Symbol("c")
u = sp.Symbol("u")


def v2int(n):
    n = abs(int(n))
    if n == 0:
        return None
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def build_R(k):
    A = sp.expand(sum(sp.chebyshevu(j, c) * (k - 2 * j) for j in range(k + 1)))
    B = sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (k - 2 * j) ** 2 for j in range(k + 1)))
    Nu = sp.expand(sp.expand((1 - c**2) ** 2 * A**2 * B).subs(c, u / 2))
    pl = sp.Poly(Nu, u)
    dd = sp.Integer(1)
    for co in pl.all_coeffs():
        dd = sp.lcm(dd, sp.denom(sp.Rational(co)))
    s2 = v2int(int(dd)) or 0
    Nint = sp.expand(Nu * (2 ** s2))
    m = k + 2
    FA = (k + 1) // 2
    roots = [2 * sp.cos(2 * sp.pi * j / m) for j in range(1, FA + 1)]
    mps, seen = [], set()
    for r in roots:
        mp = sp.minimal_polynomial(r, u)
        s = str(mp)
        if s not in seen:
            seen.add(s)
            mps.append(mp)
    Phi = sp.expand(sp.prod(mps))
    Rm = sp.rem(sp.Poly(Nint, u), sp.Poly(Phi, u))
    rp = sp.Poly(Rm.as_expr(), u)
    Rd = {e: int(co) for (e,), co in rp.terms()}
    for d in range(FA):
        Rd.setdefault(d, 0)
    return Rd, FA, m


def rightmost_hull_slope(Rd, FA):
    pts = sorted((d, v2int(Rd[d])) for d in range(FA) if v2int(Rd[d]) is not None)
    hull = []
    for p in pts:
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            if (x2 - x1) * (p[1] - y1) - (y2 - y1) * (p[0] - x1) <= 0:
                hull.pop()
            else:
                break
        hull.append(p)
    (x1, y1), (x2, y2) = hull[-2], hull[-1]
    return Fraction(y2 - y1, x2 - x1)


def main():
    print("=" * 116)
    print("F89 top-degree dominance B-11: the k=2 mod 4 class")
    print("=" * 116)
    print()
    ks = [k for k in range(6, 80) if k % 4 == 2]
    print(f"  {'k':>3} {'FA':>3} {'v2m':>3} {'V':>3} {'v2(R_0)':>7} {'content2^t':>10} "
          f"{'rslope':>7} {'V-(FA-1)':>8} {'worstdef':>8} {'atd':>4} {'DOM':>5}")
    print("  " + "-" * 86)
    alldom = True
    for k in ks:
        Rd, FA, m = build_R(k)
        V = v2int(Rd[FA - 1])
        v0 = v2int(Rd[0])
        vm = v2int(m)
        # content
        t = min(v2int(Rd[d]) for d in range(FA) if v2int(Rd[d]) is not None)
        rs = rightmost_hull_slope(Rd, FA)
        ref = V - (FA - 1)
        worst = None; atd = None
        for d in range(FA):
            vd = v2int(Rd[d])
            if vd is None:
                continue
            deficit = ref - (vd - d)
            if worst is None or deficit > worst:
                worst = deficit; atd = d
        dom = (worst <= 0)
        alldom = alldom and dom
        Vpred = 2 * vm + 2
        print(f"  {k:>3} {FA:>3} {vm:>3} {V:>3}({'=2v2m+2' if V==Vpred else '!!'}) {v0:>7} {t:>10} "
              f"{str(rs):>7} {ref:>8} {worst:>8} {atd:>4} {str(dom):>5}")
    print()
    print(f"  ===== (DOM) holds for k=2 mod 4 (k=6..78): {alldom} =====")
    print()
    print("  Here V=v2(R_top)=2 v2(m)+2 is LARGE; the polygon descends from a deep constant term.")
    print("  rslope (rightmost Newton slope) must be <=1 for (DOM). worstdef<=0 confirms (DOM).")


if __name__ == "__main__":
    main()
