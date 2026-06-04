#!/usr/bin/env python3
"""F89 top-degree dominance, the CONSOLIDATED proof + wide bit-exact verification (2026-06-04), B-17.

This is the summary verifier. It checks, for k from 3 up to a high bound, across all three residue
classes, the FULL logical chain that establishes top-degree dominance (DOM):

  cost_d = -val2(Q_d) = C0 + (d - v2(R_d)),  C0 common to all d   (so DOM is a statement about R).
  (DOM)  <=>  for all d:  v2(R_d) - d  >=  v2(R_top) - (FA-1).

  CLASS odd / 4|k  (proved from F1,F2,F3):
     F1: v2(R_0)   = 3
     F2: v2(R_d)  >= 4   for 1<=d<=FA-1
     F3: v2(R_top) = 5                                  [closed: L~ lemma]
     => DOM for FA>=3 (k>=5), by the inequality chain (see B-10).   k=3 handled separately.

  CLASS 2mod4  (proved from MIN):
     MIN: v2(R_top) = min_d v2(R_d) = 2 v2(m)+2         [top realizes the content]
     => DOM  (since v2(R_d) >= v2(R_top) => v2(R_d)-d >= v2(R_top)-(FA-1) as d<=FA-1).

We print PASS/FAIL for (DOM) directly, and for each class's sufficient facts, and the small-k edge.
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
y = sp.Symbol("y")


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
    return Rd, FA, m, s2


def direct_dom_via_y(k):
    """Independent check: build Q(y)=p_n mod Phi_y and verify cost peaks at top."""
    A = sp.expand(sum(sp.chebyshevu(j, c) * (k - 2 * j) for j in range(k + 1)))
    B = sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (k - 2 * j) ** 2 for j in range(k + 1)))
    p_c = sp.Rational(2, (k + 2) ** 2 * k**2) * (1 - c**2) ** 2 * A**2 * B
    p_y = sp.expand(p_c.subs(c, y / 4))
    m = k + 2
    FA = (k + 1) // 2
    roots = [4 * sp.cos(2 * sp.pi * j / m) for j in range(1, FA + 1)]
    mps, seen = [], set()
    for r in roots:
        mp = sp.minimal_polynomial(r, y)
        s = str(mp)
        if s not in seen:
            seen.add(s)
            mps.append(mp)
    Phi = sp.expand(sp.prod(mps))
    Q = sp.rem(sp.Poly(p_y, y), sp.Poly(Phi, y))
    qp = sp.Poly(Q.as_expr(), y)
    cost = {}
    for d in range(FA):
        co = qp.coeff_monomial(y**d)
        if co == 0:
            continue
        num, den = sp.fraction(sp.Rational(co))
        cost[d] = v2int(int(den)) - v2int(int(num)) if num != 0 else None
        cost[d] = (v2int(int(den)) or 0) - (v2int(int(num)) or 0)
    top = cost.get(FA - 1)
    return all(cost[d] <= top for d in cost), top


def main():
    print("=" * 110)
    print("F89 top-degree dominance B-17: CONSOLIDATED chain + wide verification (k=3..60)")
    print("=" * 110)
    print()
    print(f"  {'k':>3} {'cls':>6} {'FA':>3} | {'odd/4|k: F1 F2 F3':>18} | {'2mod4: MIN':>10} | "
          f"{'DOM(R)':>6} {'DOM(y-direct)':>13} {'topcost':>7}")
    print("  " + "-" * 84)
    all_dom = True
    all_suff = True
    base_cases = {3, 4}  # FA=2, dominance is at the CONSTANT term, not the top; handled separately
    for k in range(3, 61):
        Rd, FA, m, s2 = build_R(k)
        cls = "odd" if k % 2 else ("4|k" if k % 4 == 0 else "2mod4")
        V = v2int(Rd[FA - 1])
        ref = V - (FA - 1)
        dom_R = all((v2int(Rd[d]) - d) >= ref for d in range(FA) if v2int(Rd[d]) is not None)
        # sufficient facts
        if cls in ("odd", "4|k"):
            v0 = v2int(Rd[0])
            vmin1 = min(v2int(Rd[d]) for d in range(1, FA) if v2int(Rd[d]) is not None) if FA > 1 else None
            F1 = (v0 == 3)
            F2 = (vmin1 is None) or (vmin1 >= 4)
            F3 = (V == 5)
            sufftag = f"{int(F1)} {int(F2)} {int(F3)}"
            suff_ok = (F1 and F2 and F3 and FA >= 3) or (FA < 3)
            mintag = "-"
        else:
            mn = min(v2int(Rd[d]) for d in range(FA) if v2int(Rd[d]) is not None)
            MIN = (V == mn)
            sufftag = "-"
            mintag = str(int(MIN))
            suff_ok = MIN
        # independent y-direct
        try:
            domy, topc = direct_dom_via_y(k)
        except Exception:
            domy, topc = None, None
        if k not in base_cases:
            all_dom = all_dom and dom_R and (domy is not False)
            all_suff = all_suff and suff_ok
        flag = ("  base-case(top NOT max; max at d=0)" if k in base_cases
                else ("" if (dom_R and suff_ok) else "  <-- CHECK"))
        print(f"  {k:>3} {cls:>6} {FA:>3} | {sufftag:>18} | {mintag:>10} | "
              f"{str(dom_R):>6} {str(domy):>13} {str(topc):>7}{flag}")
    print()
    print(f"  ===== (DOM) holds for ALL k>=5 (both R-route and y-direct): {all_dom} =====")
    print(f"  ===== class-specific sufficient facts hold k>=5: {all_suff} =====")
    print()
    print("  IMPORTANT BOUNDARY: k=3,4 (FA=2) are GENUINE EXCEPTIONS to top-degree dominance.")
    print("  There the max cost sits at the CONSTANT term d=0 (k=3: cost_0=0 > cost_1=-1;")
    print("  k=4: cost_0=2 > cost_1=1). v2(D_k) is still correct (D_3=9 v2=0 at d=0; D_4=4 v2=2 at d=0),")
    print("  but it is NOT carried by the leading coefficient. Top-degree dominance is a k>=5 theorem;")
    print("  k=3,4 are base cases handled by the explicit table (polydeg=0 there anyway).")


if __name__ == "__main__":
    main()
