#!/usr/bin/env python3
"""F89 top-degree dominance, the THREE-FACT sufficient set (2026-06-04), probe B-10.

For odd k and 4|k, (DOM) follows from THREE integer facts about R(u) = Nint mod Phi_u:

  (F1)  v2(R_0)   = 3
  (F2)  v2(R_d)   >= 4   for all d >= 1            (equiv: R_d == 0 mod 16 for d>=1; R_0==8*odd)
  (F3)  v2(R_top) = 5

PROOF that (F1)+(F2)+(F3)+(FA>=3) => (DOM):
  (DOM) <=> for all d:  v2(R_d) - d >= v2(R_top) - (FA-1) = 5 - (FA-1) = 6 - FA.
   - d = FA-1:   v2(R_top) - (FA-1) = 6-FA.  EQUALITY (reference).            [from F3]
   - d = 0:      v2(R_0) - 0 = 3 >= 6-FA  <=>  FA >= 3.                        [from F1]
   - 1<=d<=FA-2: v2(R_d) - d >= 4 - d.  Need 4-d >= 6-FA <=> FA >= d+2,
                 which holds since d <= FA-2.                                  [from F2]
  Hence min_d (v2(R_d)-d) attained at d=FA-1 (ties only possible at d=0 when FA=3, and at
  d=FA-2 when... 4-(FA-2)=6-FA, i.e. ALWAYS a tie candidate at d=FA-2 IF v2(R_{FA-2})=4).
  So the top is a (the) minimizer; (DOM) holds, ties allowed. QED (modulo F1,F2,F3).

This probe verifies F1,F2,F3 bit-exact (the inputs to the elementary argument), and re-derives
(DOM) from them via the inequality chain (independently of the direct profile check), for odd/4|k.
F3 is the CLOSED L~ lemma (v2(L~)=5). F1 and F2 are the NEW inputs this edge needs.
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


def main():
    print("=" * 110)
    print("F89 top-degree dominance B-10: the THREE-FACT set (F1) v2(R_0)=3 (F2) v2(R_d>=1)>=4 (F3) v2(R_top)=5")
    print("=" * 110)
    print()
    ks = [k for k in range(5, 50) if (k % 2 == 1 or k % 4 == 0)]
    print(f"  {'k':>3} {'cls':>4} {'FA':>3} {'v2(R_0)':>7} {'min_{{d>=1}}v2(R_d)':>17} {'v2(R_top)':>9} "
          f"{'F1':>4} {'F2':>4} {'F3':>4} {'=>DOM':>6} {'directDOM':>9}")
    print("  " + "-" * 92)
    all_facts = True
    all_impl = True
    all_direct = True
    for k in ks:
        Rd, FA, m = build_R(k)
        cls = "odd" if k % 2 else "4|k"
        v0 = v2int(Rd[0])
        vtop = v2int(Rd[FA - 1])
        vmin1 = min(v2int(Rd[d]) for d in range(1, FA) if v2int(Rd[d]) is not None)
        F1 = (v0 == 3)
        F2 = (vmin1 >= 4)
        F3 = (vtop == 5)
        facts = F1 and F2 and F3
        # implication: derive DOM from facts via the inequality chain (treat F2 as a floor of 4)
        impl_dom = True
        if facts and FA >= 3:
            ref = vtop - (FA - 1)  # = 6-FA
            # d=0
            if not (v0 - 0 >= ref):
                impl_dom = False
            # d>=1 using floor 4
            for d in range(1, FA):
                lb = 4 - d if d <= FA - 2 else (vtop - (FA - 1))  # top uses exact
                if d == FA - 1:
                    lb = vtop - (FA - 1)
                if not (lb >= ref):
                    impl_dom = False
                    break
        else:
            impl_dom = (FA >= 3) and facts  # if facts fail, implication vacuous-but-mark
        # direct DOM check from actual coeffs
        ref = vtop - (FA - 1)
        direct = all((v2int(Rd[d]) - d) >= ref for d in range(FA) if v2int(Rd[d]) is not None)
        all_facts = all_facts and facts
        all_impl = all_impl and (impl_dom if facts else True)
        all_direct = all_direct and direct
        print(f"  {k:>3} {cls:>4} {FA:>3} {v0:>7} {vmin1:>17} {vtop:>9} "
              f"{str(F1):>4} {str(F2):>4} {str(F3):>4} {str(impl_dom):>6} {str(direct):>9}")
    print()
    print(f"  ===== F1,F2,F3 all hold (odd/4|k, k=5..49): {all_facts} =====")
    print(f"  ===== (F1+F2+F3+FA>=3 => DOM) verified: {all_impl} ;  direct DOM: {all_direct} =====")
    print()
    print("  CONCLUSION: for odd k and 4|k, top-degree dominance reduces to the three integer facts")
    print("  F1 (constant term = 8*odd), F2 (every higher coeff divisible by 16), F3 (leading = 32*odd).")
    print("  F3 = the closed L~ lemma. F1, F2 are the remaining inputs. All bit-exact here.")


if __name__ == "__main__":
    main()
