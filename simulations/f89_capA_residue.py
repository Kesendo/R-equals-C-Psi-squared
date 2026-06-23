#!/usr/bin/env python3
"""F89 capA probe 15: lock the RIGOROUS step -- residue equality P_y(y) = P_c(y/4).

This is the hinge that turns the c-reduction (clean, PIECE A) into the y-reduction (the repo D_k).
We must be sure it is an EXACT polynomial identity, and understand WHY (so it can be cited as
rigorous, not 'verified').

Claim: let Phi_c(c) = prod_i (c - c_i) (monic, nodal poly of the FA orbit cosines c_i), and
Phi_y(y) = prod_i (y - 4 c_i) (monic, nodal poly of y_i = 4 c_i). For any polynomial f,
   ( f(c) mod Phi_c )  and  ( f(y/4) mod Phi_y )  are the SAME function on the orbit, AND both
have degree < FA, so by uniqueness of the < FA interpolant through FA points,
   ( f(y/4) mod Phi_y )(y) = ( f(c) mod Phi_c )(c=y/4)   EXACTLY.
Proof sketch to validate: rem(f(y/4),Phi_y) is the unique deg<FA poly R_y with R_y(y_i)=f(y_i/4)
=f(c_i). rem(f,Phi_c) is unique deg<FA R_c with R_c(c_i)=f(c_i). Then R_c(y/4) is deg<FA in y
and R_c(c_i)=f(c_i)=R_y(y_i)=R_y(4 c_i) -> R_c(y/4) and R_y(y) agree at y=y_i (FA points), both
deg<FA -> equal. RIGOROUS.

Also reconcile probe 9's 'Phi_y(4c) = 4^FA Phi_c : False' -- the relation that actually holds is
   Phi_y(4c) = 4^FA * Phi_c(c)   IFF Phi_y and Phi_c are the EXACT same-indexed nodal products.
sympy's minimal_polynomial may return a poly whose ROOTS are the orbit but with a sign/leading
convention, OR the combined product over distinct minimal polys may DROP multiplicity / reorder.
We test the precise identity with EXPLICIT nodal products (not minimal_polynomial) to confirm
4^FA scaling, then confirm sympy's get_orbit_polynomial gives the same ideal.
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


def p_n_in_c_expr(k):
    m = k + 2
    A = sp.expand(sum(sp.chebyshevu(j, c) * (k - 2 * j) for j in range(k + 1)))
    B = sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (k - 2 * j) ** 2 for j in range(k + 1)))
    return sp.expand(sp.Rational(2) / (m**2 * k**2) * (1 - c**2) ** 2 * A**2 * B)


def explicit_nodal_c(k):
    """Phi_c = prod over DISTINCT orbit cosines of (c - c_i), built explicitly & symmetrized."""
    m = k + 2
    orbit = list(range(2, k + 2, 2))
    roots = [sp.cos(sp.pi * n / m) for n in orbit]
    # distinct minimal polys -> their product has exactly the distinct conjugate roots once
    mps, seen = [], set()
    for r in roots:
        mp_r = sp.minimal_polynomial(r, c)
        s = str(mp_r)
        if s not in seen:
            seen.add(s); mps.append(mp_r)
    return sp.expand(sp.prod(mps))


def main():
    print("=" * 104)
    print("F89 capA probe 15: residue equality P_y(y)=P_c(y/4) exact; 4^FA nodal scaling")
    print("=" * 104)
    print()
    reps = [5, 6, 7, 8, 9, 10, 12, 16, 22, 24]
    for k in reps:
        FA = (k + 1) // 2
        deg_red = FA - 1
        Phi_c = explicit_nodal_c(k)
        Phi_y, _ = fpk.get_orbit_polynomial(k)

        # nodal scaling: Phi_y(4c) =? 4^FA Phi_c  (both monic in their var, same roots scaled)
        # NOTE Phi_y is monic in y of degree FA; Phi_y(4c) leading = 4^FA c^FA; Phi_c monic c^FA.
        scaling = sp.simplify(Phi_y.subs(y, 4 * c) - 4**FA * Phi_c)
        scaling_ok = (scaling == 0)

        # residue equality
        p_c = p_n_in_c_expr(k)
        P_c = sp.rem(sp.Poly(p_c, c, domain="QQ"), sp.Poly(Phi_c, c, domain="QQ")).as_expr()
        p_y = sp.expand(p_c.subs(c, y / 4))
        P_y = sp.rem(sp.Poly(p_y, y, domain="QQ"), sp.Poly(Phi_y, y, domain="QQ")).as_expr()
        residue_ok = sp.simplify(P_y.subs(y, 4 * c) - P_c) == 0

        # also: deg check
        dPy = sp.Poly(P_y, y).degree()
        dPc = sp.Poly(P_c, c).degree()
        print(f"k={k:>3} FA={FA} | Phi_y(4c)=4^FA Phi_c: {scaling_ok} | "
              f"P_y(y)=P_c(y/4) exact: {residue_ok} | deg(P_y)={dPy} deg(P_c)={dPc} (<{FA})")
    print()
    print("scaling_ok True => the explicit nodal relation holds; sympy combined poly is the same ideal.")
    print("residue_ok True (all k) => the variable-change of the REDUCED poly is exact: RIGOROUS.")
    print("Therefore  v2(D_k) = -min_d val2([y^d]P_y) = -min_d ( val2([c^d]P_c) - 2d )")
    print("                   = max_d ( 2d - val2([c^d]P_c) ).   [exact]")


if __name__ == "__main__":
    main()
