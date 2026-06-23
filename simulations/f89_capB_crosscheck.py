#!/usr/bin/env python3
"""F89 capB: definitive cross-check vs the established pipeline + constant-5 source (2026-06-04).

Cross-check the capB master-formula v2(D_k) against the ACTUAL D_k from the established
pipeline f89_pathk_symbolic_derivation.extract_path_polynomial (the ground-truth reduction
mod Phi_y in y), for the prompt's full k-set spanning v2=0..5.  This proves the capB
reformulation computes the SAME D_k, not a different object.

Also probe the constant '5' (= v2(Rm_top) for odd k and 4|k) for a clean source:
Rm_top for odd k is an ODD multiple of 32 (v2=5).  We factor 2^5 out and report the odd
cofactor's small structure; and we report the leading coeff of Nint (=k^4*2^{s2-4}) vs Rm_top
to see the '5' = '4 (from (1-c^2)^2 /16 cleared) + 1 (one residual reduction factor)'.
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
    Nu = sp.expand(sp.expand((1 - c**2) ** 2 * A**2 * B).subs(c, u / 2))
    pl = sp.Poly(Nu, u)
    d = sp.Integer(1)
    for co in pl.all_coeffs():
        d = sp.lcm(d, sp.denom(sp.Rational(co)))
    s2 = v2(int(d))
    return sp.expand(Nu * (2 ** s2)), s2


def build_phi_u(k):
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
    return sp.expand(sp.prod(mps))


def capB_v2Dk(k):
    Nint, s2 = build_Nint_s2(k)
    Phi = build_phi_u(k)
    Rm = sp.rem(sp.Poly(Nint, u), sp.Poly(Phi, u))
    Rmtop = int(sp.Poly(Rm.as_expr(), u).all_coeffs()[0])
    m = k + 2
    FA = (k + 1) // 2
    return (2 * v2(m) + 2 * v2(k) - 1) + s2 + (FA - 1) - v2(Rmtop), Rmtop, s2


def main():
    print("=" * 104)
    print("F89 capB: DEFINITIVE cross-check vs established pipeline D_k (prompt k-set)")
    print("=" * 104)
    print()
    # prompt's exact k-set
    kset = [5, 7, 9, 6, 10, 14, 12, 20, 8, 24, 16, 40, 32]
    print(f"  {'k':>3} {'v2k':>3} {'D_k(pipeline)':>14} {'v2(D_k)pipe':>11} "
          f"{'v2(D_k)capB':>11} {'target':>6} {'match':>6}")
    print("  " + "-" * 70)
    okall = True
    for k in kset:
        _, D_pipe, _ = fpk.extract_path_polynomial(k)
        v2_pipe = v2(D_pipe)
        v2_capB, Rmtop, s2 = capB_v2Dk(k)
        polydeg = max(0, (k - 5) // 2)
        target = polydeg + 2 * v2(k) - min(v2(k), 2)
        match = (v2_pipe == v2_capB == target)
        okall = okall and match
        print(f"  {k:>3} {v2(k):>3} {D_pipe:>14} {v2_pipe:>11} "
              f"{v2_capB:>11} {target:>6} {('OK' if match else 'X'):>6}")
    print()
    print(f"  ===== capB master == established-pipeline D_k == closed-form target: {okall} =====")
    print()

    # constant-5 source
    print("-" * 104)
    print("constant '5' = v2(Rm_top) for odd k and 4|k: factor 2^5 out, inspect odd cofactor")
    print("-" * 104)
    print(f"  {'k':>3} {'class':>8} {'Rm_top':>12} {'Rm_top/2^5 (odd)':>16}")
    for k in [5, 7, 9, 11, 13, 12, 20, 8, 24, 16, 32]:
        _, Rmtop, s2 = capB_v2Dk(k)
        cls = "odd" if k % 2 == 1 else "4|k"
        odd = Rmtop // (2 ** v2(Rmtop))
        print(f"  {k:>3} {cls:>8} {Rmtop:>12} {odd:>16}")
    print()
    print("  The '5' = 4 (from (1-c^2)^2 = (4-u^2)^2/16 contributing /16 cleared to numerator")
    print("  giving 2^4 in Nint for odd k via s2=4, OR the explicit 2^4 even-k pull-out) PLUS a")
    print("  single residual factor 2 from the leading reduction step. The exact derivation of")
    print("  this '+1' (and that it is EXACTLY 1, capping the cancellation) is the localized gap.")
    print()


if __name__ == "__main__":
    main()
