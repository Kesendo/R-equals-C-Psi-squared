#!/usr/bin/env python3
"""F89 edge A, step 6: recast SLOPE1 in the integer u-variable as a clean divisibility.

c-route reduced coeffs R_d are RATIONAL (denominators powers of 2). Moving to u = 2c makes
the orbit polynomial monic-integer (Washington) and the reduced numerator INTEGER. Relation:

  G(c) = Nu(u)|_{u=2c},   Nint(u) = 2^{s2} Nu(u) in Z[u],   Phi_u(u) = 2^{FA} Phi_c(u/2) monic,
  Rm(u) = Nint(u) mod Phi_u(u) in Z[u],   Rm_d = [u^d] Rm.

Claim (verify): R_d = Rm_d * 2^{d-s2}, i.e.  val2(R_d) = val2(Rm_d) + d - s2.   (REL)

Then SLOPE1 (val2(R_d) >= d + delta) becomes, for d>=1:

      val2(Rm_d) >= s2 + delta.                                                  (DIV)

s2 is the integer 2-clearing power of Nu (rigorous: s2 = 4 odd k... but here Nu uses the
even-pull-out kappa-j for even k; we recompute s2 honestly).  So DIV says: every reduced
INTEGER coefficient Rm_d (d>=1) is divisible by 2^{s2+delta}.  This is a *uniform divisibility*
of the integer reduced polynomial by a fixed 2-power, the cleanest possible target.

We tabulate s2, delta, the threshold s2+delta, min_{d>=1} val2(Rm_d), and the d=0 value,
verifying (REL) and (DIV) bit-exact, k=5..36.
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


def build(k):
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
    # c-route
    roots_c = [sp.cos(2 * sp.pi * j / m) for j in range(1, FA + 1)]
    mps, seen = [], set()
    for r in roots_c:
        mp_r = sp.minimal_polynomial(r, c)
        s = str(mp_r)
        if s not in seen:
            seen.add(s); mps.append(mp_r)
    Phi_c = sp.expand(sp.prod(mps))
    Phi_c = sp.expand(Phi_c / sp.Rational(sp.Poly(Phi_c, c).LC()))
    Rc = sp.rem(sp.Poly(G, c, domain="QQ"), sp.Poly(Phi_c, c, domain="QQ"))
    R_d = [sp.Rational(Rc.coeff_monomial(c**d)) for d in range(Rc.degree() + 1)]
    # u-route
    Nu = sp.expand(G.subs(c, u / 2))
    dlcm = sp.Integer(1)
    for co in sp.Poly(Nu, u).all_coeffs():
        dlcm = sp.lcm(dlcm, sp.denom(sp.Rational(co)))
    s2 = v2_int(int(dlcm))
    Nint = sp.expand(Nu * 2 ** s2)
    roots_u = [2 * sp.cos(2 * sp.pi * j / m) for j in range(1, FA + 1)]
    mps2, seen2 = [], set()
    for r in roots_u:
        mp_r = sp.minimal_polynomial(r, u)
        s = str(mp_r)
        if s not in seen2:
            seen2.add(s); mps2.append(mp_r)
    Phi_u = sp.expand(sp.prod(mps2))
    Rm = sp.rem(sp.Poly(Nint, u, domain="ZZ"), sp.Poly(Phi_u, u, domain="ZZ"))
    Rm_d = [int(sp.Poly(Rm.as_expr(), u).coeff_monomial(u**d)) for d in range(sp.Poly(Rm.as_expr(), u).degree() + 1)]
    return R_d, Rm_d, s2, FA, m


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
    print("F89 edge A: u-route recast.  REL: val2(R_d)=val2(Rm_d)+d-s2.  DIV: val2(Rm_d)>=s2+delta (d>=1)")
    print("=" * 112)
    print()
    print(f"  {'k':>3} {'cls':>6} {'FA':>3} {'s2':>3} {'delta':>5} {'s2+delta':>8} "
          f"{'min v2(Rm_d>=1)':>15} {'v2(Rm_0)':>9} {'REL_ok':>6} {'DIV_ok':>6}")
    print("  " + "-" * 96)
    allrel = True
    alldiv = True
    for k in range(5, 37):
        R_d, Rm_d, s2, FA, m = build(k)
        # REL check
        rel_ok = True
        for d in range(min(len(R_d), len(Rm_d))):
            vrm = v2_int(Rm_d[d]) if Rm_d[d] != 0 else None
            vr = signed_val2(R_d[d])
            if vrm is None and vr is None:
                continue
            if vrm is None or vr is None:
                rel_ok = False
                continue
            if vr != vrm + d - s2:
                rel_ok = False
        allrel = allrel and rel_ok
        delta = delta_closed(k)
        thr = s2 + delta
        v2Rm = [v2_int(r) if r != 0 else None for r in Rm_d]
        d1 = [v2Rm[d] for d in range(1, len(Rm_d)) if v2Rm[d] is not None]
        minRm1 = min(d1) if d1 else None
        v0 = v2Rm[0]
        div_ok = (minRm1 is not None and minRm1 >= thr)
        alldiv = alldiv and div_ok
        print(f"  {k:>3} {cls_of(k):>6} {FA:>3} {s2:>3} {delta:>5} {thr:>8} "
              f"{str(minRm1):>15} {str(v0):>9} {str(rel_ok):>6} {str(div_ok):>6}")
    print()
    print(f"  REL  val2(R_d) = val2(Rm_d) + d - s2  (all k):  {allrel}")
    print(f"  DIV  val2(Rm_d) >= s2 + delta  for d>=1  (all k):  {alldiv}")
    print()
    print("  s2: for odd k s2=4 always (clearing (1-c^2)^2 A^2 B at half-integers); even k s2 also 4 here")
    print("  (the pull-out leaves (1-c^2)^2 A_til^2 B_til, same clearing). DIV is the integer")
    print("  divisibility: reduced numerator Rm(u) is divisible by 2^{s2+delta} away from the constant.")


if __name__ == "__main__":
    main()
