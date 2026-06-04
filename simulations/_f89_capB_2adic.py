#!/usr/bin/env python3
"""F89 capB: honest 2-adic accounting, odd/2 parts separated correctly (2026-06-04).

LESSON from the failed u-probe: p(u) is NOT an integer polynomial (the prefactor
2/(m^2 k^2) carries ODD denominators m_odd^2, k_odd^2 too), so 'R = 2^P p_red_u' was
not integer and the int() casts were garbage.  Work HONESTLY in y, localized at 2.

D_k = odd(k)^2 * 2^{v2(D_k)} is established.  We only want v2(D_k).  Compute it as the
honest max over j of v2(denom([y^j] p_red_y)).  We then localize the WHOLE pipeline at 2:
every rational number x = 2^{v2(x)} * (odd/odd); we track only v2.

Key clean facts to establish (bit-exact, v2=0..6):
  (1) v2(D_k) = polydeg + 2 v2(k) - min(v2(k),2),   polydeg = max(0,(k-5)//2).  [the goal]
  (2) The orbit polynomial Phi_y(y) (monic integer in y) has a known 2-adic structure:
        Phi_y(y) = 2^FA * Phi_u(y/2), Phi_u monic integer in u=y/2.
      We read v2 of Phi_y's coefficients and its discriminant.
  (3) The numerator p(y)*(m^2 k^2 / 2) = (1-c^2)^2 A^2 B with c=y/4 -- call it Ntot(y) --
      is an integer-coefficient? NO (c=y/4 gives /4^deg). Track v2 of its coeff denoms.

THE ACTUAL DECOMPOSITION we verify (all in y, honest):
  Let q(y) := p_red_y(y) = sum_j a_j y^j.
  Define for the TOP coeff a_{FA-1}: v2(denom a_{FA-1}).
  Empirically the max over j of v2(denom a_j) is the relevant v2(D_k); we record argmax,
  and the value, and confirm = target.  We ALSO record v2(denom a_j) for ALL j to see the
  shape (whether top-degree, or some interior j, drives it -- the u-probe wrongly forced top).
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


def build_p_y(k):
    m = k + 2
    A = sp.expand(sum(sp.chebyshevu(j, c) * (k - 2 * j) for j in range(k + 1)))
    B = sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (k - 2 * j) ** 2 for j in range(k + 1)))
    p = sp.Rational(2, m * m * k * k) * (1 - c**2) ** 2 * A**2 * B
    return sp.expand(sp.expand(p).subs(c, y / 4))


def build_phi_y(k):
    m = k + 2
    orbit = list(range(2, k + 2, 2))
    roots = [4 * sp.cos(sp.pi * n / m) for n in orbit]
    mps, seen = [], set()
    for r in roots:
        mp_r = sp.minimal_polynomial(r, y)
        s = str(mp_r)
        if s not in seen:
            seen.add(s)
            mps.append(mp_r)
    return sp.expand(sp.prod(mps))


def denom_v2_by_degree(expr, var):
    poly = sp.Poly(sp.expand(expr), var)
    deg = poly.degree()
    d = {deg - i: co for i, co in enumerate(poly.all_coeffs())}
    out = {}
    for j in range(deg + 1):
        co = d.get(j, sp.Integer(0))
        out[j] = None if co == 0 else v2(int(sp.denom(sp.Rational(co))))
    return out, deg


def main():
    print("=" * 118)
    print("F89 capB 2-adic honest: v2(denom a_j) shape, argmax, and target")
    print("=" * 118)
    print()
    reps = [5, 7, 9, 11, 13, 15, 6, 10, 14, 22, 26, 12, 20, 28, 8, 24, 40, 16, 48, 32, 64]
    print(f"  {'k':>3} {'v2k':>3} {'FA':>3} {'v2(D_k)':>7} {'target':>6} "
          f"{'argmax_j':>8} {'ok':>3}  v2(denom a_j) j=0..FA-1")
    print("  " + "-" * 100)
    okall = True
    for k in reps:
        p_y = build_p_y(k)
        phi_y = build_phi_y(k)
        q = sp.rem(sp.Poly(p_y, y), sp.Poly(phi_y, y)).as_expr()
        dv, deg = denom_v2_by_degree(q, y)
        vmax = max(x for x in dv.values() if x is not None)
        argmax = max(j for j, x in dv.items() if x == vmax)
        polydeg = max(0, (k - 5) // 2)
        target = polydeg + 2 * v2(k) - min(v2(k), 2)
        ok = (vmax == target)
        okall = okall and ok
        shape = " ".join(("." if dv[j] is None else str(dv[j])) for j in range(deg + 1))
        print(f"  {k:>3} {v2(k):>3} {(k+1)//2:>3} {vmax:>7} {target:>6} "
              f"{argmax:>8} {('OK' if ok else 'X'):>3}  [{shape}]")
    print()
    print(f"  ALL v2(D_k) == target bit-exact: {okall}")
    print()
    print("  Reading argmax: where in the degree ladder the maximal-denominator coefficient sits.")
    print("  For odd k it climbs interiorly; the top coefficient is NOT always the driver.")
    print()


if __name__ == "__main__":
    main()
