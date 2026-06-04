#!/usr/bin/env python3
"""F89 edge A, step 1b: cross-check the cost_d profile against the bit-exact y-pipeline.

The _f89_edge_scout / _f89_edgeA_profile read coefficients of sp.rem(p_y, Phi_y) then
nsimplify them. For high FA (k>=22) nsimplify on a degree-(FA-1) rational coefficient with
a huge orbit polynomial can MISREAD the 2-adic valuation (k=22 showed dom=False, top_cost=0,
which contradicts the proven v2(D_22)=9). This script recomputes cost_d directly from the
exact rational coefficients of the reduced polynomial, with NO nsimplify, using domain='QQ'
arithmetic, and confirms:
  (a) max_d cost_d == v2(D_k) from the closed form (sanity on the whole profile), and
  (b) the argmax set, top dominance, and margin, exactly.

If the exact route gives dom=True at k=22 while the nsimplify route gave dom=False, the
failure was a precision artifact, and we keep the EXACT route for all downstream tables.
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
        return 10**9
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def signed_val2_exact(q):
    """q a sympy Rational (exact). signed v2; None for 0."""
    q = sp.Rational(q)
    if q == 0:
        return None
    num, den = sp.fraction(q)
    return v2_int(int(num)) - v2_int(int(den))


def reduced_coeffs_exact(k):
    """Exact rational coeffs Q_d (low->high) of p_n(y) mod Phi_y, all in QQ, no nsimplify."""
    p_y = fpk.derive_p_n_in_y(k)
    mp, _ = fpk.get_orbit_polynomial(k)
    # force exact rational domain
    pp = sp.Poly(sp.expand(p_y), y, domain="QQ")
    phi = sp.Poly(sp.expand(mp), y, domain="QQ")
    R = sp.rem(pp, phi)
    deg = R.degree()
    out = []
    for d in range(deg + 1):
        out.append(sp.Rational(R.coeff_monomial(y**d)))
    return out


def cls_of(k):
    if k % 2 == 1:
        return "odd"
    if k % 4 == 0:
        return "4|k"
    return "2mod4"


def closed_form_v2D(k):
    polydeg = max(0, (k - 5) // 2)
    vk = v2_int(k)
    return polydeg + 2 * vk - min(vk, 2)


def main():
    print("=" * 110)
    print("F89 edge A: EXACT cost_d profile (domain=QQ, no nsimplify) vs closed-form v2(D_k)")
    print("=" * 110)
    print()
    ks = [5, 7, 9, 11, 13, 15, 17, 19, 21, 23,
          6, 10, 14, 18, 22,
          8, 12, 16, 20, 24]
    allok = True
    for k in sorted(ks):
        coeffs = reduced_coeffs_exact(k)
        FA = len(coeffs)
        # cost_d = -val2(Q_d): denominator 2-power contributed by coeff d
        vals = [signed_val2_exact(c) for c in coeffs]
        costs = [(-v if v is not None else None) for v in vals]
        present = [d for d in range(FA) if costs[d] is not None]
        top = FA - 1
        top_cost = costs[top]
        vD = max(costs[d] for d in present)
        argmaxes = [d for d in present if costs[d] == vD]
        dom = (top in argmaxes)
        below = [costs[d] for d in present if d < top]
        margin = (top_cost - max(below)) if below else None
        C_top = top_cost - top
        cf = closed_form_v2D(k)
        match_vD = (vD == cf)
        ok = dom and match_vD
        allok = allok and ok
        prof = " ".join(f"{d}:{('.' if costs[d] is None else costs[d])}" for d in range(FA))
        flag = "" if ok else "   <<< CHECK"
        print(f"  k={k:>2} [{cls_of(k):>5}] FA={FA:>2} dom={str(dom):>5} top_cost={top_cost:>3} "
              f"C={C_top:>3} maxcost={vD:>3} v2D_cf={cf:>3} match={str(match_vD):>5} "
              f"margin={str(margin):>4}{flag}")
        print(f"        cost_d: {prof}")
    print()
    print(f"  ===== EXACT route: all dom True AND maxcost==closed-form v2(D_k): {allok} =====")


if __name__ == "__main__":
    main()
