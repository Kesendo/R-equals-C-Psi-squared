#!/usr/bin/env python3
"""F89 Path-D: where does the cancellation cap 2^min(v2(k),2) come from? (scout, 2026-06-04)

D_k = 2^max(0, floor((k-5)/2)) * k^2 / 2^min(v2(k), 2). The poly-degree 2-power is the
Chebyshev degree growth (understood). This scout traces the v2-specific part a(k) :=
v2(D_k) - polydeg = 2*v2(k) - min(v2(k),2), to localise the cap.

Hypothesis (from the (k-2j) structure): in
    p_n(c) = (2 / (m^2 k^2)) (1-c^2)^2 A(c)^2 B(c),
    A(c) = sum_j U_j(c)(k-2j),  B(c) = sum_j U_j(c)^2 (k-2j)^2,
every (k-2j) is even iff k is even, = 2(kappa - j) with kappa = k/2. So for even k:
    A = 2*Atil,  B = 4*Btil   =>   p_n = (8 / (m^2 kappa^2)) (1-c^2)^2 Atil^2 Btil.
That converts k^2 -> kappa^2 = (k/2)^2 in the denominator and injects a FIXED 2^2 into the
numerator (the 8 vs 2), independent of how divisible k is. The pull-out works exactly ONCE
because (kappa - j) alternates parity (no uniform second factor). So:
  - odd k:        no pull-out, D-amplitude-2power a(k) = 0.
  - k = 2*odd:    one pull-out, kappa odd, a(k) = 1.
  - 4 | k:        one pull-out, kappa even, a(k) = 2, and it STAYS 2 for any higher v2.
This file checks the pull-out identity symbolically and traces v2 at each stage.
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


def v2(n: int) -> int:
    n = abs(int(n))
    if n == 0:
        return 0
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def denom_lcm_v2(expr_in_y: sp.Expr) -> int:
    """v2 of the LCM of coefficient denominators of a polynomial in y."""
    poly = sp.Poly(sp.expand(expr_in_y), y)
    d = sp.Integer(1)
    for coef in poly.all_coeffs():
        d = sp.lcm(d, sp.denom(sp.nsimplify(coef)))
    return v2(int(d))


def A_poly(k, var, shift=0):
    """A(c) with (k - 2j) -> (k - 2j) - 2*shift, i.e. pulling kappa = k/2 uses shift via scale."""
    return sp.expand(sum(sp.chebyshevu(j, var) * (k - 2 * j) for j in range(k + 1)))


def B_poly(k, var):
    return sp.expand(sum(sp.chebyshevu(j, var) ** 2 * (k - 2 * j) ** 2 for j in range(k + 1)))


def Atil_poly(k, var):
    kap = k // 2
    return sp.expand(sum(sp.chebyshevu(j, var) * (kap - j) for j in range(k + 1)))


def Btil_poly(k, var):
    kap = k // 2
    return sp.expand(sum(sp.chebyshevu(j, var) ** 2 * (kap - j) ** 2 for j in range(k + 1)))


def main():
    print("=" * 96)
    print("F89 Path-D cancellation cap 2^min(v2,2): scout trace")
    print("=" * 96)
    print()

    # --- 1. The (k-2j) = 2(kappa-j) pull-out, checked symbolically for even k -------------
    print("-" * 96)
    print("1. Pull-out identity for even k:  A = 2*Atil,  B = 4*Btil   (=> p_n gains 2^3 in num,")
    print("   denominator k^2 -> kappa^2 = (k/2)^2). Works once; (kappa-j) alternates parity.")
    print("-" * 96)
    print(f"  {'k':>3} {'v2':>3} {'A==2*Atil':>10} {'B==4*Btil':>10}  {'kappa=k/2':>10} {'v2(kappa)':>10}")
    for k in (4, 6, 8, 10, 12, 16, 20, 24):
        a_ok = sp.simplify(A_poly(k, c) - 2 * Atil_poly(k, c)) == 0
        b_ok = sp.simplify(B_poly(k, c) - 4 * Btil_poly(k, c)) == 0
        kap = k // 2
        print(f"  {k:>3} {v2(k):>3} {str(a_ok):>10} {str(b_ok):>10}  {kap:>10} {v2(kap):>10}")
    print()
    print("  => the factor 2 comes out of EVERY (k-2j) exactly once for even k; the residual")
    print("     (kappa - j) is a mix of even and odd, so no second uniform factor of 2 exists.")
    print("     This is why the even-specific 2-content is FIXED (the 2^3 in '8'), not v2-growing.")
    print()

    # --- 2. v2 trace through the reduction, representative k per v2 class -----------------
    print("-" * 96)
    print("2. v2 at each stage:  unreduced p_n(y) denominator  ->  reduced D_k.")
    print("   a(k) := v2(D_k) - polydeg should equal 2*v2(k) - min(v2(k),2) (the capped part).")
    print("-" * 96)
    header = f"  {'k':>3} {'v2':>3} {'polydeg':>8} {'v2(unred)':>10} {'v2(D_k)':>8} {'a(k)':>5} {'2v2-min(v2,2)':>14} {'ok':>4}"
    print(header)
    reps = [5, 7, 9, 6, 10, 14, 12, 20, 8, 24, 16]
    for k in reps:
        p_y = fpk.derive_p_n_in_y(k)
        v_unred = denom_lcm_v2(p_y)
        _, D_k, _ = fpk.extract_path_polynomial(k)
        polydeg = max(0, (k - 5) // 2)
        a = v2(D_k) - polydeg
        target = 2 * v2(k) - min(v2(k), 2)
        ok = (a == target)
        print(f"  {k:>3} {v2(k):>3} {polydeg:>8} {v_unred:>10} {v2(D_k):>8} {a:>5} {target:>14} {('OK' if ok else 'FAIL'):>4}")
    print()
    print("  reading: a(k) is 0 (odd), 1 (k=2*odd), 2 (4|k) and then STAYS 2. The unreduced")
    print("  denominator carries a large 2-power (from U_j(y/4) ~ 1/2^j and (1-c^2)^2 ~ 1/2^8);")
    print("  the reduction mod the orbit minimal polynomial cancels it down to exactly the cap.")
    print()
    print("  Open hinge for a full proof: the residual a(k) = min(v2,2) is the parity content of")
    print("  the single pull-out (1 for k=2*odd via odd kappa, 2 for 4|k via even kappa, capped")
    print("  because the pull-out happens once). The reduction itself contributes no further")
    print("  v2-dependence -- to be pinned down on the orbit-polynomial side.")


if __name__ == "__main__":
    main()
