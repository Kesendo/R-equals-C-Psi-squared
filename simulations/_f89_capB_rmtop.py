#!/usr/bin/env python3
"""F89 capB: closed form of Rm_top = lead coeff of (Nint mod Phi_u), the cap carrier (2026-06-04).

Master (bit-exact): v2(D_k) = (2v2m+2v2k-1) + s2 + (FA-1) - v2(Rm_top).
Explicit: s2=4 (odd k), 0 (even k); FA-1=deg.  Atomic: v2(Rm_top).  Observed:
  odd k:        v2(Rm_top)=5
  k=2 mod4:     v2(Rm_top)=2 v2(m)+2
  4|k:          v2(Rm_top)=5
We (a) re-confirm at larger k, (b) compute Rm_top EXACTLY (the integer, not just v2) and
look for its full closed form, (c) test an ALTERNATIVE clean characterization:

Rm_top is the leading coeff of the remainder of Nint (deg 2k+4 in u) mod Phi_u (deg FA, monic).
By polynomial division the remainder's top coeff = Nint reduced; equivalently it is the
(FA-1)-th elementary-symmetric / divided-difference of Nint over the orbit roots.

We use Newton's identity for the leading remainder coeff: if Nint(u) = sum_i n_i u^i, then
reducing u^p mod Phi_u for p>=FA via the recurrence u^FA = -(sum of lower Phi_u coeffs)*u^...,
the leading (u^{FA-1}) coeff of the remainder is  sum_{p} n_p * h_{p-(FA-1)}  where h are the
"reduction weights" (complete homogeneous-like in the roots).  Concretely it equals the
(FA-1) coeff which we just read off sp.rem.  We instead VERIFY the v2 closed form holds and
factor Rm_top to see the 2-part structure vs the odd part.
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


def build_Nint_s2(k):
    A = sp.expand(sum(sp.chebyshevu(j, c) * (k - 2 * j) for j in range(k + 1)))
    B = sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (k - 2 * j) ** 2 for j in range(k + 1)))
    Nc = sp.expand((1 - c**2) ** 2 * A**2 * B)
    Nu = sp.expand(Nc.subs(c, u / 2))
    poly = sp.Poly(Nu, u)
    d = sp.Integer(1)
    for co in poly.all_coeffs():
        d = sp.lcm(d, sp.denom(sp.Rational(co)))
    s2 = v2(int(d))
    return sp.expand(Nu * (2 ** s2)), s2


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


def rm_top(k):
    Nint, s2 = build_Nint_s2(k)
    phi = build_phi_u(k)
    Rm = sp.rem(sp.Poly(Nint, u), sp.Poly(phi, u))
    rp = sp.Poly(Rm.as_expr(), u)
    return int(rp.all_coeffs()[0]), s2, rp.degree()


def main():
    print("=" * 110)
    print("F89 capB: Rm_top closed form (cap carrier).  factor out 2-part, inspect odd part.")
    print("=" * 110)
    print()
    reps = [5, 7, 9, 11, 13, 15, 17, 19,
            6, 10, 14, 18, 22, 26, 30, 34, 38,
            12, 20, 28, 36, 44,
            8, 24, 40, 56,
            16, 48, 32, 64]
    print(f"  {'k':>3} {'v2k':>3} {'kmod4':>5} {'v2m':>4} {'v2(Rm_top)':>10} "
          f"{'predicted':>9} {'odd(Rm_top)':>12} {'ok':>3}")
    print("  " + "-" * 70)
    okall = True
    for k in reps:
        Rt, s2, deg = rm_top(k)
        vk, vm = v2(k), v2(k + 2)
        if k % 2 == 1:
            pred = 5
        elif k % 4 == 2:
            pred = 2 * vm + 2
        else:  # 4|k
            pred = 5
        oddpart = abs(Rt) // (2 ** v2(Rt)) if Rt != 0 else 0
        ok = (v2(Rt) == pred)
        okall = okall and ok
        # show odd part mod small to look for pattern, but mostly just confirm v2
        odstr = str(oddpart)
        if len(odstr) > 12:
            odstr = odstr[:9] + "..."
        print(f"  {k:>3} {vk:>3} {k%4:>5} {vm:>4} {v2(Rt):>10} "
              f"{pred:>9} {odstr:>12} {('OK' if ok else 'X'):>3}")
    print()
    print(f"  ALL v2(Rm_top) closed form bit-exact: {okall}")
    print()
    print("  Summary of v2(Rm_top):")
    print("    odd k:      5")
    print("    k=2 mod4:   2 v2(m) + 2    (m=k+2, v2(m)>=2 here)")
    print("    4|k:        5")
    print()


if __name__ == "__main__":
    main()
