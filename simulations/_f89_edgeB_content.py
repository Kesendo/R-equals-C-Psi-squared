#!/usr/bin/env python3
"""F89 top-degree dominance, deriving the content + F1/F2 via G mod Phi over Z (2026-06-04), probe B-12.

We attack the THREE integer facts (odd/4|k):  F1 v2(R_0)=3, F2 v2(R_d>=1)>=4, F3 v2(R_top)=5,
and the 2mod4 (MIN): content(R)=v2(R_top)=2 v2(m)+2.

Approach: R = Nint mod Phi over Z. The content / per-coeff v2 come from (i) Nint's own 2-content and
(ii) the reduction by Phi. We expose:
  - the 2-content of Nint (cont(Nint) = min v2 of its coeffs) and its top/bottom,
  - Phi mod 2 (its shape controls how reduction moves 2-adic mass),
  - DIRECT factorization R(u) = 2^t * Rprime(u), and Rprime mod 2 (support) and mod 4,
  - whether Rprime is, mod 2, equal to a UNIT*(Phi mod 2)-coprime monomial pattern.

The cleanest derivable handle we test:  Nint(u) = 2^{s2} (1-(u/2)^2)^2 A(u/2)^2 B(u/2).
The factor (1-(u/2)^2)^2 = ((4-u^2)/4)^2 = (4-u^2)^2 / 16 contributes 2^{s2-4}*(4-u^2)^2 ... and the
'4 from clearing (4-u^2)^2/16' is the proof's stated '4 of the 5'. So we separately track:
   Nint = (4-u^2)^2 * Atil(u)^2 * Btil(u) / 2^{w}    with Atil=2^a A(u/2) etc integerized,
and read where the residual single '+1' (the 5th factor of 2 in R_top, and the 8=2^3 content) is born.
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


def content_poly(poly_u):
    cs = [int(co) for co in sp.Poly(poly_u, u).all_coeffs()]
    vs = [v2int(x) for x in cs if x != 0]
    return min(vs) if vs else None


def build(k):
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
    Phi = sp.Poly(sp.expand(sp.prod(mps)), u)
    Rm = sp.rem(sp.Poly(Nint, u), Phi)
    rp = sp.Poly(Rm.as_expr(), u)
    Rd = {e: int(co) for (e,), co in rp.terms()}
    for d in range(FA):
        Rd.setdefault(d, 0)
    return Rd, FA, m, Phi, sp.Poly(Nint, u), s2


def main():
    print("=" * 116)
    print("F89 top-degree dominance B-12: content of Nint, of R, Phi(0), and the birth of 8 and 32")
    print("=" * 116)
    print()
    print(f"  {'k':>3} {'cls':>6} {'FA':>3} {'s2':>3} {'cont(Nint)':>10} {'cont(R)':>8} "
          f"{'v2(R_0)':>7} {'v2(R_top)':>9} {'v2(Phi0)':>8} {'v2(Phitop=1)':>12} {'v2(disc)':>8}")
    print("  " + "-" * 100)
    for k in [5, 7, 9, 11, 13, 15, 17, 8, 12, 16, 20, 24, 6, 10, 14, 18, 22]:
        Rd, FA, m, Phi, Nint, s2 = build(k)
        cls = "odd" if k % 2 else ("4|k" if k % 4 == 0 else "2mod4")
        cN = content_poly(Nint.as_expr())
        cR = content_poly(sum(Rd[d] * u**d for d in range(FA)))
        v0 = v2int(Rd[0]); vtop = v2int(Rd[FA - 1])
        phi0 = v2int(int(Phi.eval(0)))
        disc = v2int(int(sp.discriminant(Phi)))
        print(f"  {k:>3} {cls:>6} {FA:>3} {s2:>3} {str(cN):>10} {str(cR):>8} "
              f"{str(v0):>7} {str(vtop):>9} {str(phi0):>8} {'0':>12} {str(disc):>8}")
    print()
    print("  Reading: cont(Nint) is the 2-content the numerator carries before reduction. cont(R) is")
    print("  what survives mod Phi. For odd/4|k cont(R)=3 at the constant term; for 2mod4 cont(R)=V at top.")
    print("  v2(Phi0)=v2(Phi(0))=v2(product of nodes)=v2(N(nodes)); ties to the constant-term mass.")
    print()

    # The (4-u^2)^2/16 split: show s2 and the 'first 4' explicitly
    print("  --- the (4-u^2)^2 split: Nint = 2^{s2}*( (4-u^2)^2/16 )*A(u/2)^2*B(u/2). Track 2-powers. ---")
    print(f"  {'k':>3} {'cls':>6} {'s2':>3} {'v2(A(u/2) int-clear)':>20} {'v2(B intclear)':>15}")
    for k in [5, 7, 9, 8, 12, 16, 6, 10, 14]:
        m = k + 2
        Apoly = sp.expand(sum(sp.chebyshevu(j, c) * (k - 2 * j) for j in range(k + 1)).subs(c, u/2))
        Bpoly = sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (k - 2 * j) ** 2 for j in range(k + 1)).subs(c, u/2))
        cls = "odd" if k % 2 else ("4|k" if k % 4 == 0 else "2mod4")
        # integerize A(u/2): lcm of denominators
        def clearpow(P):
            dd = sp.Integer(1)
            for co in sp.Poly(P, u).all_coeffs():
                dd = sp.lcm(dd, sp.denom(sp.Rational(co)))
            return v2int(int(dd))
        print(f"  {k:>3} {cls:>6} {'-':>3} {str(clearpow(Apoly)):>20} {str(clearpow(Bpoly)):>15}")
    print()


if __name__ == "__main__":
    main()
