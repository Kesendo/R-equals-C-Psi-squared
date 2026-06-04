#!/usr/bin/env python3
"""F89 capB: the odd-k baseline v2(D_k)=polydeg, proven via 2-integrality of T_u (2026-06-04).

For odd k, m=k+2 odd, v2(k)=0.  Target: v2(D_k)=polydeg=max(0,(k-5)//2), cap=0.

Cleanest object: work in u=y/2.  Phi_u monic integer, Z[u]=O_K.  The numerator in u,
  Nu(u) = (1-(u/2)^2)^2 A(u/2)^2 B(u/2),
has 2-power denominators only (the c=u/2 substitution).  Write Nu(u) = M(u)/2^{s2},
M integer-primitive (odd content).  Reduction: M mod Phi_u =: Rm (integer poly).
p_n in u: p_red_u = (2/(m^2 k^2)) * Rm / 2^{s2}.  For ODD k, m,k odd, so v2 of every coeff
of p_red_u is  v2(Rm_j) + 1 - s2.  The y-denominator: a_j = [u^j]p_red_u / 2^j, so
  v2(denom a_j) = max(0, s2 - 1 - v2(Rm_j) + j).
v2(D_k) = max_j of that.  We verify it equals polydeg, AND that the max is at top j=FA-1
with v2(Rm_{FA-1}) such that s2-1-v2(Rm_top)+(FA-1) = polydeg.

We tabulate s2, v2(Rm_top), FA-1, and confirm  s2 - 1 - v2(Rm_top) + (FA-1) = polydeg
for odd k -- a CLEAN integer identity, the baseline.  We also confirm the interior
coefficients never exceed it (so top is the true max).
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


def v2(n) -> int:
    n = abs(int(n))
    if n == 0:
        return 10**9
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def build_M_u(k):
    A = sp.expand(sum(sp.chebyshevu(j, c) * (k - 2 * j) for j in range(k + 1)))
    B = sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (k - 2 * j) ** 2 for j in range(k + 1)))
    Nc = sp.expand((1 - c**2) ** 2 * A**2 * B)
    Nu = sp.expand(Nc.subs(c, u / 2))
    poly = sp.Poly(Nu, u)
    d = sp.Integer(1)
    for co in poly.all_coeffs():
        d = sp.lcm(d, sp.denom(sp.Rational(co)))
    s2 = v2(int(d))
    M = sp.expand(Nu * (2 ** s2))
    return M, s2


def build_phi_u(k):
    m = k + 2
    FA = (k + 1) // 2
    roots = [2 * sp.cos(2 * sp.pi * j / m) for j in range(1, FA + 1)]
    mps, seen = [], set()
    for r in roots:
        mp_r = sp.minimal_polynomial(r, u)
        s = str(mp_r)
        if s not in seen:
            seen.add(s)
            mps.append(mp_r)
    return sp.expand(sp.prod(mps))


def Rm_coeffs(k):
    M, s2 = build_M_u(k)
    phi = build_phi_u(k)
    rem = sp.rem(sp.Poly(M, u), sp.Poly(phi, u))
    rp = sp.Poly(rem.as_expr(), u)
    deg = rp.degree()
    co = {deg - i: int(c0) for i, c0 in enumerate(rp.all_coeffs())}
    return co, deg, s2


def main():
    print("=" * 110)
    print("F89 capB ODD-k baseline: s2 - 1 - v2(Rm_top) + (FA-1) = polydeg  (cap = 0)")
    print("=" * 110)
    print()
    odds = [5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25]
    print(f"  {'k':>3} {'FA':>3} {'s2':>4} {'v2(Rm_top)':>10} {'FA-1':>5} "
          f"{'s2-1-v2top+(FA-1)':>17} {'polydeg':>7} {'top is max?':>11} {'ok':>3}")
    print("  " + "-" * 86)
    okall = True
    for k in odds:
        co, deg, s2 = Rm_coeffs(k)
        FA = (k + 1) // 2
        v2top = v2(co[deg]) if co.get(deg, 0) != 0 else None
        # per-j  denom v2 = max(0, s2-1-v2(Rm_j)+j)
        denomv = {j: max(0, s2 - 1 - v2(co[j]) + j) for j in co if co[j] != 0}
        vmax = max(denomv.values())
        argmax = max(j for j, vv in denomv.items() if vv == vmax)
        polydeg = max(0, (k - 5) // 2)
        top_val = s2 - 1 - v2top + deg
        ok = (top_val == polydeg and argmax == deg and vmax == polydeg)
        okall = okall and ok
        print(f"  {k:>3} {FA:>3} {s2:>4} {v2top:>10} {deg:>5} "
              f"{top_val:>17} {polydeg:>7} {str(argmax==deg):>11} {('OK' if ok else 'X'):>3}")
    print()
    print(f"  ALL odd-k baseline (top_val=polydeg, top is the max) bit-exact: {okall}")
    print()
    print("  So for odd k, v2(D_k)=polydeg with NO cap, and the driver is the leading reduced")
    print("  coefficient. The cap min(v2,2) is the EXTRA that even k adds on top of this baseline.")
    print()
    # show s2 and v2top closed forms for odd k
    print("  closed-form pieces (odd k):")
    print(f"  {'k':>3} {'s2':>4} {'v2(Rm_top)':>10} {'s2-v2top':>8} {'expected s2=2*deg? ':>20}")
    for k in odds:
        co, deg, s2 = Rm_coeffs(k)
        v2top = v2(co[deg])
        print(f"  {k:>3} {s2:>4} {v2top:>10} {s2 - v2top:>8} {2*deg:>20}")
    print()


if __name__ == "__main__":
    main()
