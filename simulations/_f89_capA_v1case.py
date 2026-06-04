#!/usr/bin/env python3
"""F89 capA probe 14: pin the v2(k)=1 (k=2*odd) residual -- the only case where v2(Ltil) != (FA-1)+1.

From probe 13: for v2(k)=0 and v2(k)>=2, v2(Ltil) = (FA-1)+1 (a v2-independent degree-growth
constant), and the cap is PROVEN. For v2(k)=1 (kappa odd, kappa+1 even), v2(Ltil)-(FA-1)
VARIED (4,2,6,4). We must characterise it.

Identity (rigorous): s_top = 1 - 2 v2(kappa) - 2 v2(kappa+1) + v2(Ltil), and for v2(k)=1,
v2(kappa)=0, so
    s_top = 1 - 2 v2(kappa+1) + v2(Ltil),    and we KNOW s_top = (FA-1)+1 (offset(1)=1).
Therefore the EXACT requirement is
    v2(Ltil) = (FA-1) + 2 v2(kappa+1).
Probe 13 hinted v2(Ltil)=(FA-1)+1+2 v2(kap+1) (off by one at k=6). Let's just MEASURE
v2(Ltil)-(FA-1) and 2 v2(kappa+1) side by side over MANY v2=1 cases and see which holds:
    candidate (i):  v2(Ltil)-(FA-1) = 2 v2(kappa+1)        [needed for s_top=(FA-1)+1]
    candidate (ii): v2(Ltil)-(FA-1) = 1 + 2 v2(kappa+1)    [probe-13 guess]
The one that holds bit-exact across all v2=1 k is the closure for this case.

kappa = k/2 odd; kappa+1 = k/2+1; v2(kappa+1) = v2(k/2+1) = v2(k+2)-1 = v2(m)-1 (since m=k+2,
k=2*odd => m=2*(odd+1), v2(m) = 1+v2(odd+1)... let's just compute v2(kappa+1) directly).
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


def main():
    print("=" * 110)
    print("F89 capA probe 14: v2(k)=1 case -- v2(Ltil)-(FA-1) vs 2 v2(kappa+1)")
    print("=" * 110)
    print()
    hdr = (f"  {'k':>3} {'kap':>4} {'kap+1':>5} {'v2(kap+1)':>9} {'FA-1':>4} | "
           f"{'v2(Ltil)':>8} {'-(FA-1)':>7} | {'2v2(kap+1)':>10} {'cand(i)':>7} | "
           f"{'1+2v2(kap+1)':>12} {'cand(ii)':>8} | {'s_top':>5} {'expect':>6}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))

    # all k = 2*odd up to ~70
    reps = [k for k in range(6, 72, 4)]  # 6,10,14,...,70 are 2*odd (k/2 odd)
    for k in reps:
        kap = k // 2
        if kap % 2 == 0:
            continue
        FA = (k + 1) // 2
        v2kap1 = v2_int(kap + 1)

        # s_top from full p_n
        m = k + 2
        A = sp.expand(sum(sp.chebyshevu(j, c) * (k - 2 * j) for j in range(k + 1)))
        B = sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (k - 2 * j) ** 2 for j in range(k + 1)))
        p_c = sp.expand(sp.Rational(2) / (m**2 * k**2) * (1 - c**2) ** 2 * A**2 * B)
        s_top = val2(lead_reduced(p_c, k))

        v2Lt = val2(lead_reduced(Gtil(k), k))
        delta = v2Lt - (FA - 1)
        cand_i = (delta == 2 * v2kap1)
        cand_ii = (delta == 1 + 2 * v2kap1)

        print(f"  {k:>3} {kap:>4} {kap+1:>5} {v2kap1:>9} {FA-1:>4} | "
              f"{v2Lt:>8} {delta:>7} | {2*v2kap1:>10} {str(cand_i):>7} | "
              f"{1+2*v2kap1:>12} {str(cand_ii):>8} | {s_top:>5} {(FA-1)+1:>6}")

    print()
    print("  Need: s_top = (FA-1)+1 (offset(1)=1). Via identity s_top=1-2v2(kap+1)+v2(Ltil),")
    print("  this REQUIRES v2(Ltil) = (FA-1) + 2 v2(kap+1)  i.e. cand(i)=True for all rows.")
    print("  If cand(i) holds bit-exact: v2=1 case CLOSED (v2(Ltil)=(FA-1)+2v2(kap+1)).")


if __name__ == "__main__":
    main()
