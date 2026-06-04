#!/usr/bin/env python3
"""F89 capA probe 11 (lean): closed form for s_top = val2(leading coeff of P_c).

Master reduction (probe 10, ALL OK):  v2(D_k) = 2(FA-1) - s_top, s_top = val2 of the LEADING
coefficient L of P_c = p_n(c) mod Phi_c, the maximum of (2d - s_d) ALWAYS attained at d=FA-1.

This probe confirms the closed form
    s_top = (FA-1) + offset(v2k),   offset(v) = 2 - v - max(0, v-2),
across v2 = 0..6 (k chosen to span the cap and beyond), and hence
    v2(D_k) = 2(FA-1) - s_top = (FA-1) - offset(v2k) = polydeg + a(k),
    a(k) = 2 v2(k) - min(v2(k),2),     using FA-1 = polydeg + 2 (k>=5).

It also reports the FULL leading coefficient L (exact rational) so the analytic attack on
val2(L) has the concrete object. L is computed as the c^{FA-1} coefficient of rem(p_n, Phi_c)
(fast, exact over QQ).
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


def val2(r):
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


def v2_int(n):
    n = abs(int(n)); v = 0
    while n % 2 == 0:
        n //= 2; v += 1
    return v


def p_n_in_c_expr(k):
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


def main():
    print("=" * 116)
    print("F89 capA probe 11 (lean): s_top = val2(leading coeff L of P_c) = (FA-1) + offset(v2k)")
    print("=" * 116)
    print()

    def offset(v):
        return 2 - v - max(0, v - 2)

    reps = [5, 7, 9, 11, 13, 6, 10, 14, 22, 26, 12, 20, 36, 8, 24, 40, 16, 48, 32, 64]
    allok = True
    for k in reps:
        FA = (k + 1) // 2
        deg_red = FA - 1
        polydeg = max(0, (k - 5) // 2)
        v2k = v2_int(k)

        P_c = sp.rem(sp.Poly(p_n_in_c_expr(k), c, domain="QQ"),
                     sp.Poly(orbit_poly_c(k), c, domain="QQ")).as_expr()
        L = sp.Rational(sp.Poly(P_c, c, domain="QQ").coeff_monomial(c**deg_red))

        s_top = val2(L)
        s_top_pred = deg_red + offset(v2k)
        v2D = 2 * deg_red - s_top
        target = polydeg + (2 * v2k - min(v2k, 2))
        ok = (s_top == s_top_pred) and (v2D == target)
        allok = allok and ok

        # show L compactly: odd-part factored
        num, den = sp.fraction(L)
        Lstr = f"{int(num)}/{int(den)}"
        if len(Lstr) > 26:
            Lstr = f"v2num={v2_int(num)} v2den={v2_int(den)}"
        print(f"k={k:>3} v2k={v2k} FA-1={deg_red} polydeg={polydeg} | s_top={s_top:>3} pred={s_top_pred:>3} "
              f"| v2D={v2D:>3} target={target:>3} | {'OK' if ok else 'X':>2} | L={Lstr}")

    print()
    print(f"ALL OK: {allok}")
    print()
    print("offset(v) = 2 - v - max(0,v-2):  {0:2, 1:1, 2:0, 3:-2, 4:-4, 5:-6, 6:-8}.")
    print("=> v2(D_k) = (FA-1) - offset = polydeg + 2 - offset = polydeg + (2 v2 - min(v2,2)) = polydeg + a(k).")


if __name__ == "__main__":
    main()
