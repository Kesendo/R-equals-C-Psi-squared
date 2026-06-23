#!/usr/bin/env python3
"""F89 capA: MASTER verification table for the cap decomposition (spans v2 = 0..6).

Consolidates the proof chain into ONE table, every intermediate quantity bit-exact:

  Stage 0 (given):   even-k pull-out p_n = 8/(m^2 kappa^2)(1-c^2)^2 Atil^2 Btil  [spec, symbolic-OK]
  Stage 1 (RIGOROUS identity):  P_y(y) = P_c(y/4)  [residue equality]
        => v2(D_k) = max_d ( 2d - val2([c^d] P_c) )
  Stage 2 (VERIFIED, top-wins):  max attained at d = FA-1
        => v2(D_k) = 2(FA-1) - s_top,  s_top = val2(leading coeff L of P_c)
  Stage 3 (VERIFIED closed form): s_top = (FA-1) + offset(v2k), offset(v)=2-v-max(0,v-2)
  Stage 4 (RIGOROUS algebra): FA-1 = polydeg + 2  =>  v2(D_k) = polydeg + a(k),
        a(k) = 2 v2(k) - min(v2(k), 2).

Columns prove each '=>'. We also recompute PIECE A (c-reduced denom = kappa^2) and the
ground-truth v2(D_k) from the sympy y-pipeline; all must agree.
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


def p_n_in_c_expr(k):
    m = k + 2
    A_poly = sp.expand(sum(sp.chebyshevu(j, c) * (k - 2 * j) for j in range(k + 1)))
    B_poly = sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (k - 2 * j) ** 2 for j in range(k + 1)))
    return sp.expand(sp.Rational(2) / (m**2 * k**2) * (1 - c**2) ** 2 * A_poly**2 * B_poly)


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


def denom_lcm(poly_expr, var):
    poly = sp.Poly(poly_expr, var, domain="QQ")
    d = sp.Integer(1)
    for coef in poly.all_coeffs():
        d = sp.lcm(d, sp.denom(sp.Rational(coef)))
    return int(d)


def offset(v):
    return 2 - v - max(0, v - 2)


def main():
    print("=" * 124)
    print("F89 capA MASTER TABLE: v2(D_k) = polydeg + a(k) via the c-reduction / residue-equality chain")
    print("=" * 124)
    print()
    hdr = (f"  {'k':>3} {'v2':>2} {'FA-1':>4} {'pdeg':>4} | "
           f"{'Den_c':>7} {'=kap^2':>6} | "
           f"{'top_wins':>8} {'s_top':>5} {'(FA-1)+off':>10} | "
           f"{'v2D_chain':>9} {'v2D_pipe':>8} {'pdeg+a(k)':>9} {'ALL':>4}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))

    reps = [5, 7, 9, 11, 13, 6, 10, 14, 22, 26, 12, 20, 36, 8, 24, 40, 16, 48, 32, 64]
    allok = True
    for k in reps:
        FA = (k + 1) // 2
        deg_red = FA - 1
        polydeg = max(0, (k - 5) // 2)
        v2k = v2_int(k)
        a_k = 2 * v2k - min(v2k, 2)

        p_c = p_n_in_c_expr(k)
        Phi_c = orbit_poly_c(k)
        P_c = sp.rem(sp.Poly(p_c, c, domain="QQ"), sp.Poly(Phi_c, c, domain="QQ")).as_expr()
        Pc = sp.Poly(P_c, c, domain="QQ")

        # PIECE A: c-reduced denom
        Den_c = denom_lcm(P_c, c)
        kap2 = (k // 2) ** 2 if k % 2 == 0 else None
        pieceA_ok = (Den_c == kap2) if kap2 is not None else (v2_int(Den_c) == 0)

        # signed profile + top-wins + s_top
        s = {d: val2(sp.Rational(Pc.coeff_monomial(c**d))) for d in range(deg_red + 1)}
        vals = [(2 * d - s[d], d) for d in range(deg_red + 1) if s[d] is not None]
        v2D_chain, _ = max(vals)
        argmaxes = [d for (v, d) in vals if v == v2D_chain]
        top_wins = (deg_red in argmaxes) and (2 * deg_red - s[deg_red] == v2D_chain)
        s_top = s[deg_red]
        s_top_pred = deg_red + offset(v2k)

        # ground truth from y-pipeline
        _, D_k, _ = fpk.extract_path_polynomial(k)
        v2D_pipe = v2_int(D_k)

        formula = polydeg + a_k
        ok = (pieceA_ok and top_wins and s_top == s_top_pred
              and v2D_chain == v2D_pipe == formula)
        allok = allok and ok

        kapstr = str(kap2) if kap2 is not None else "odd"
        print(f"  {k:>3} {v2k:>2} {deg_red:>4} {polydeg:>4} | "
              f"{Den_c:>7} {kapstr:>6} | "
              f"{str(top_wins):>8} {s_top:>5} {s_top_pred:>10} | "
              f"{v2D_chain:>9} {v2D_pipe:>8} {formula:>9} {('OK' if ok else 'FAIL'):>4}")

    print()
    print(f"  ALL STAGES OK (every k): {allok}")
    print()
    print("  Legend:")
    print("   Den_c = kap^2     : PIECE A, c-reduced denominator is exactly (k/2)^2 (even k), odd-denom-free (odd k)")
    print("   top_wins          : v2(D_k)=max_d(2d - val2 [c^d]P_c) attained at top degree FA-1  [Stage 2]")
    print("   s_top=(FA-1)+off  : Stage 3 closed form, offset(v)=2-v-max(0,v-2)")
    print("   v2D_chain=pipe=formula : chain value = y-pipeline ground truth = polydeg + (2 v2 - min(v2,2))")


if __name__ == "__main__":
    main()
