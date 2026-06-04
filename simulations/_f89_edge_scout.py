#!/usr/bin/env python3
"""F89 cap proof, the last edge: top-degree dominance of the denominator profile (scout, 2026-06-04).

The two-route chain reduces v2(D_k) to v2(D_k) = max_d ( -val2(Q_d) ), where Q(y) = p_n(y) mod Phi_y
is the reduced amplitude and val2 is the signed 2-adic valuation (the c-variable's +2d cancels:
[c^d]P_c = Q_d 4^d so 2d - val2([c^d]) = -val2(Q_d)). 'Top-degree dominance' is the claim that
this max is attained at d = FA-1, i.e. the LEADING coefficient Q_{FA-1} carries the largest
denominator 2-power, with margin >= 1 over every lower coefficient.

This scout lays out the full per-coefficient profile (d, val2(Q_d), cost_d = -val2(Q_d)) to expose
the structure: is the cost a concave/monotone profile peaking at the top? What is the margin?
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

y = sp.Symbol("y")


def v2_int(n: int) -> int:
    n = abs(int(n))
    if n == 0:
        return 0
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def signed_val2(q) -> object:
    """signed 2-adic valuation of a rational; None for 0."""
    q = sp.nsimplify(q)
    if q == 0:
        return None
    num, den = sp.fraction(q)
    return v2_int(int(num)) - v2_int(int(den))


def reduced_coeffs(k):
    """coefficients Q_d (low->high) of Q(y) = p_n(y) mod Phi_y."""
    p_y = fpk.derive_p_n_in_y(k)
    mp, _ = fpk.get_orbit_polynomial(k)
    R = sp.rem(sp.Poly(p_y, y), sp.Poly(mp, y))
    poly = sp.Poly(R.as_expr(), y)
    deg = poly.degree()
    coeffs = {}
    for (e,), co in poly.terms():
        coeffs[e] = co
    return [coeffs.get(d, sp.Integer(0)) for d in range(deg + 1)]


def main():
    print("=" * 100)
    print("F89 last edge: top-degree dominance of the denominator profile  v2(D_k)=max_d(-val2(Q_d))")
    print("=" * 100)
    print()
    for k in (5, 7, 9, 11, 6, 8, 12, 16):
        coeffs = reduced_coeffs(k)
        FA = len(coeffs)
        vals = [signed_val2(c) for c in coeffs]
        # cost_d = -val2(Q_d): the denominator 2-power that coefficient d contributes to D_k
        costs = [(-v if v is not None else None) for v in vals]
        present = [(d, costs[d]) for d in range(FA) if costs[d] is not None]
        top_cost = costs[FA - 1]
        others = [c for (d, c) in present if d != FA - 1]
        vD = max(c for (d, c) in present)
        argmax_top = (top_cost == vD)
        margin = (top_cost - max(others)) if others else None
        print(f"  k={k:>2}  FA={FA:>2}  v2(D_k)={vD:>3}  top_cost={top_cost}  "
              f"argmax_at_top={argmax_top}  margin={margin}")
        # full profile
        prof = "  ".join(f"{d}:{costs[d]}" for d in range(FA))
        print(f"       cost_d (d:cost): {prof}")
        # slope-2 bound from the top: cost_d <= top_cost - 2*(FA-1-d)  would be the tight edge;
        # margin>=1 means cost_d <= top_cost - 2*(FA-1-d) - 1 for d<FA-1
        slack = [None if costs[d] is None else (top_cost - 2 * (FA - 1 - d)) - costs[d] for d in range(FA)]
        print(f"       slack to slope-2 edge (top - 2*(FA-1-d) - cost_d, >=0 means below edge): "
              + "  ".join(f"{d}:{slack[d]}" for d in range(FA)))
        print()
    print("  reading: if argmax_at_top is True with margin>=1 for all k, the leading coefficient")
    print("  dominates and v2(D_k)=2(FA-1)-val2(L~) holds. The slack row shows how far each lower")
    print("  coefficient sits below the slope-2 edge anchored at the top (looking for a clean pattern).")


if __name__ == "__main__":
    main()
