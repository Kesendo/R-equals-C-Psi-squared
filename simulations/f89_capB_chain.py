#!/usr/bin/env python3
"""F89 capB: the COMPLETE proof chain, every piece verified bit-exact (2026-06-04).

THEOREM:  v2(D_k) = polydeg + 2 v2(k) - min(v2(k),2),  polydeg=max(0,(k-5)//2).

PROOF CHAIN (each line a separately-verifiable claim):

 [S2] Z[u]=O_K (u=2cos(2 pi/m), Washington's theorem) => for the monic integer orbit
      polynomial Phi_u, reducing an integer polynomial mod Phi_u yields an integer polynomial.

 [Master] v2(D_k) = (2 v2(m) + 2 v2(k) - 1) + s2 + (FA-1) - v2(Rm_top),  where
      Nu(u) = (1-(u/2)^2)^2 A(u/2)^2 B(u/2)  [the numerator],  s2 = its 2-clearing power,
      Nint = 2^{s2} Nu  [integer], Rm = Nint mod Phi_u [integer by S2], Rm_top = lead coeff,
      and the maximum coefficient-denominator is at the TOP degree FA-1.   [S1]

 [C1] s2 = 4 (k odd),  s2 = 0 (k even).                     [degree-growth, given]

 [C2/C3] v2(Rm_top) = 5            (k odd OR 4|k)
                    = 2 v2(m) + 2  (k = 2 mod4).            [the CAP carrier]

 Substituting (and using m=k+2 so: k odd=>v2m=0; k=2mod4=>v2m=v2(k+2)>=2; 4|k=>v2m=1):
   odd k:    (0+0-1)+4+(FA-1)-5      = FA-3   = polydeg            = polydeg+2*0-0
   k=2mod4:  (2v2m+2-1)+0+(FA-1)-(2v2m+2) = FA-2 = polydeg+1       = polydeg+2*1-1
   4|k:      (2+2v2k-1)+0+(FA-1)-5   = polydeg+2v2k-2             = polydeg+2v2k-min(v2,2)

This script verifies [S1],[Master],[C1],[C2/C3] and the final theorem, bit-exact, for the
prompt's k-set spanning v2=0..5, and a few extra.
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


def v2(n) -> int:
    n = abs(int(n))
    if n == 0:
        return 10**9
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def build(k):
    A = sp.expand(sum(sp.chebyshevu(j, c) * (k - 2 * j) for j in range(k + 1)))
    B = sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (k - 2 * j) ** 2 for j in range(k + 1)))
    Nu = sp.expand(sp.expand((1 - c**2) ** 2 * A**2 * B).subs(c, u / 2))
    pl = sp.Poly(Nu, u)
    d = sp.Integer(1)
    for co in pl.all_coeffs():
        d = sp.lcm(d, sp.denom(sp.Rational(co)))
    s2 = v2(int(d))
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
    deg = rp.degree()
    Rmtop = int(rp.all_coeffs()[0])
    # honest v2(D_k) and argmax
    pref = sp.Rational(2, m * m * k * k)
    p_red_u = sp.expand(pref * Rm.as_expr() / (2 ** s2))
    p_red_y = sp.expand(p_red_u.subs(u, y / 2))
    pyp = sp.Poly(p_red_y, y)
    coy = {pyp.degree() - i: co for i, co in enumerate(pyp.all_coeffs())}
    denomv = {j: v2(int(sp.denom(sp.Rational(coy[j])))) for j in coy if coy[j] != 0}
    v2Dk = max(denomv.values())
    argmax = max(j for j, vv in denomv.items() if vv == v2Dk)
    return dict(s2=s2, FA=FA, deg=deg, Rmtop=Rmtop, v2Dk=v2Dk, argmax=argmax, m=m)


def main():
    print("=" * 116)
    print("F89 capB COMPLETE CHAIN: theorem + every piece, bit-exact")
    print("=" * 116)
    print()
    # prompt's k-set spanning v2=0..5, plus a couple
    reps = [5, 7, 9, 6, 10, 14, 12, 20, 8, 24, 16, 40, 32,
            11, 13, 18, 22, 28, 36, 48, 64]
    print(f"  {'k':>3} {'v2k':>3} {'v2m':>4} {'s2':>3} {'predC1':>6} {'v2(Rm_top)':>10} "
          f"{'predC2/3':>8} {'masterv2D':>9} {'argmax==top':>11} {'v2(D_k)':>7} {'target':>6} {'ALL_ok':>6}")
    print("  " + "-" * 104)
    okall = True
    for k in sorted(set(reps)):
        a = build(k)
        vk, vm = v2(k), v2(k + 2)
        predC1 = 4 if k % 2 == 1 else 0
        if k % 2 == 1 or k % 4 == 0:
            predC23 = 5
        else:
            predC23 = 2 * vm + 2
        master = (2 * vm + 2 * vk - 1) + a["s2"] + (a["FA"] - 1) - v2(a["Rmtop"])
        polydeg = max(0, (k - 5) // 2)
        target = polydeg + 2 * vk - min(vk, 2)
        c1_ok = (a["s2"] == predC1)
        c23_ok = (v2(a["Rmtop"]) == predC23)
        s1_ok = (a["argmax"] == a["deg"])
        master_ok = (master == a["v2Dk"])
        thm_ok = (a["v2Dk"] == target)
        allok = c1_ok and c23_ok and s1_ok and master_ok and thm_ok
        okall = okall and allok
        print(f"  {k:>3} {vk:>3} {vm:>4} {a['s2']:>3} {predC1:>6} {v2(a['Rmtop']):>10} "
              f"{predC23:>8} {master:>9} {str(s1_ok):>11} {a['v2Dk']:>7} {target:>6} "
              f"{('OK' if allok else 'X'):>6}")
    print()
    print(f"  ============  ALL CLAIMS BIT-EXACT: {okall}  ============")
    print()
    print("  Per-class collapse of the master formula to the target (algebra, not fit):")
    print("    odd k   (v2m=0):  (0+0-1)+4+(FA-1)-5          = FA-3 = polydeg           [cap 0]")
    print("    k=2mod4 (v2m>=2): (2v2m+1)+0+(FA-1)-(2v2m+2)  = FA-2 = polydeg+1         [cap 1]")
    print("    4|k     (v2m=1):  (2v2k+1)+0+(FA-1)-5         = polydeg+2v2k-2           [cap 2]")
    print()


if __name__ == "__main__":
    main()
