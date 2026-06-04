#!/usr/bin/env python3
"""F89 capA probe 10: the TRUE object -- signed 2-adic Newton polygon of P_c (reduced-in-c poly).

Corrected understanding (probe 9): P_y(y) = P_c(y/4) EXACTLY (residue equality), so
   [y^d] P_y = ([c^d] P_c) / 4^d   ->   v2([y^d]P_y) = val2([c^d]P_c) - 2d,
where val2(.) is the SIGNED 2-adic valuation of a rational (num v2 minus den v2). Then
   v2(D_k) = - min_d v2([y^d]P_y) = max_d ( 2d - val2([c^d]P_c) ).

probe 4/9 only saw DENOMINATOR v2; the NUMERATOR 2-content of [c^d]P_c is what tames the 4^d.
This probe computes the full signed profile  s_d := val2([c^d]P_c)  and verifies:

  (i)   v2(D_k) = max_d (2d - s_d)
  (ii)  the max is attained at d = FA-1 (top), i.e. 2(FA-1) - s_{FA-1} = v2(D_k)
  (iii) s_{FA-1} (signed val2 of the LEADING coeff of P_c) has the closed form
            s_{FA-1} = 2(FA-1) - v2(D_k)  ==>  to be matched to 2 v2(kap) - R ... we just REPORT
        s_{FA-1} and 2(FA-1)-v2(D_k) and check equal, and tabulate s_{FA-1} vs (v2k, FA).

Then PIECE B becomes: characterise val2 of the leading coefficient of P_c. Since denom(P_c)=kappa^2
(PIECE A) and the leading coeff's DENOMINATOR v2 is (from probe 9 c-profile top) exactly... we
report whether top-coeff denom = kappa^2's full 2-part or less, plus its numerator 2-content.
"""
from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import f89_pathk_symbolic_derivation as fpk  # noqa: E402

c = sp.Symbol("c")
y = sp.Symbol("y")


def val2(r):
    """Signed 2-adic valuation of a rational. +inf for 0."""
    r = sp.Rational(r)
    if r == 0:
        return None
    num, den = sp.fraction(r)
    def v2i(n):
        n = abs(int(n)); v = 0
        while n % 2 == 0:
            n //= 2; v += 1
        return v
    return v2i(num) - v2i(den)


def p_n_in_c(k):
    m = k + 2
    A_poly = sp.expand(sum(sp.chebyshevu(j, c) * (k - 2 * j) for j in range(k + 1)))
    B_poly = sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (k - 2 * j) ** 2 for j in range(k + 1)))
    return sp.expand(sp.Rational(2) / (m**2 * k**2) * (1 - c**2) ** 2 * A_poly**2 * B_poly)


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


def v2_int(n):
    n = abs(int(n)); v = 0
    while n % 2 == 0:
        n //= 2; v += 1
    return v


def main():
    print("=" * 116)
    print("F89 capA probe 10: signed 2-adic Newton polygon of P_c; v2(D_k)=max_d(2d - s_d), attained at top")
    print("=" * 116)
    print()

    reps = [5, 7, 9, 11, 6, 10, 14, 22, 12, 20, 8, 24, 16, 32]
    allok = True
    for k in reps:
        FA = (k + 1) // 2
        deg_red = FA - 1
        polydeg = max(0, (k - 5) // 2)
        v2k = v2_int(k)

        P_c = sp.rem(sp.Poly(p_n_in_c(k), c, domain="QQ"),
                     sp.Poly(orbit_poly_c(k), c, domain="QQ")).as_expr()
        Pc = sp.Poly(P_c, c, domain="QQ")
        s = {}
        for d in range(deg_red + 1):
            s[d] = val2(sp.Rational(Pc.coeff_monomial(c**d)))

        # v2(D_k) = max_d (2d - s_d)
        vals = [(2 * d - s[d], d) for d in range(deg_red + 1) if s[d] is not None]
        v2D_pred, argmax_d = max(vals)
        argmaxes = sorted([d for (v, d) in vals if v == v2D_pred])

        _, D_k, _ = fpk.extract_path_polynomial(k)
        v2D = v2_int(D_k)

        top_at_FAm1 = (2 * deg_red - s[deg_red] == v2D) if s[deg_red] is not None else False
        s_top = s[deg_red]
        ok = (v2D_pred == v2D) and top_at_FAm1
        allok = allok and ok

        sprof = [("." if s[d] is None else s[d]) for d in range(deg_red + 1)]
        print(f"k={k:>3} v2k={v2k} FA={FA} polydeg={polydeg} | v2D={v2D} pred=max_d(2d-s_d)={v2D_pred} "
              f"at d={argmaxes} | top d={deg_red} wins: {top_at_FAm1} | ok={ok}")
        print(f"   s_d (signed val2 of [c^d]P_c): {sprof}")
        print(f"   s_top = s_{deg_red} = {s_top}   2*deg_red - s_top = {2*deg_red - s_top if s_top is not None else 'n/a'} (=v2D)")
        print()

    print(f"ALL OK: {allok}")
    print()
    print("If top always wins: v2(D_k) = 2(FA-1) - s_top, with s_top = signed val2 of the LEADING")
    print("coefficient of the c-reduced polynomial. PIECE B reduces to a closed form for s_top.")


if __name__ == "__main__":
    main()
