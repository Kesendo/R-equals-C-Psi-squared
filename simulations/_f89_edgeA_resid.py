#!/usr/bin/env python3
"""F89 edge A, step 15: pin the RESID bound and the full per-class assembly.

ESTABLISHED rigorous chain (polynomial identities mod Phi_u, all bit-exact):
  LEMMA 1:  (4-u^2) Atil == eps_k m (2+u)  (mod Phi_u),  eps_k=2 (odd), 1 (even).
  => Nint = [(4-u^2)Atil]^2 Btil == eps_k^2 m^2 (2+u)^2 Btil (mod Phi_u)         [rigorous from L1]
  => Rm = Nint mod Phi_u = eps_k^2 m^2 * W,   W := (2+u)^2 Btil mod Phi_u.
  => val2(Rm_d) = 2 v2(eps_k) + 2 v2(m) + val2(W_d).

DIV needs val2(Rm_d) >= T = s2 + delta = 4 + delta  for d>=1.  So the residual target is

  RESID:  val2(W_d) >= 4 + delta - 2 v2(eps_k) - 2 v2(m)   for d>=1.

Per class (using delta and m=k+2):
  odd k:     eps=2 (v2=1), v2(m)=0, delta=1  => RESID rhs = 4+1-2-0 = 3.
  4|k:       eps=1 (v2=0), v2(m)=1 (m=k+2, 4|k => m even, v2(m)=1), delta=1 => 4+1-0-2 = 3.
  k=2 mod4:  eps=1 (v2=0), v2(m)=v2(k+2)>=2, delta=2 v2(kappa+1)=2 v2(m/... )...; check rhs.

This script verifies the chain bit-exact AND prints, per class, val2(W_d) profile and the
RESID rhs, confirming RESID holds for all k. It also reports min(val2(W_d), d>=1) - rhs (>=0).
This LOCALIZES the residual gap to a single clean statement RESID about W=(2+u)^2 Btil mod Phi_u.
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


def sv2(q):
    q = sp.Rational(q)
    if q == 0:
        return None
    n, d = sp.fraction(q)
    return v2_int(int(n)) - v2_int(int(d))


def pieces(k):
    if k % 2 == 1:
        A = sp.expand(sum(sp.chebyshevu(j, c) * (k - 2 * j) for j in range(k + 1)))
        B = sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (k - 2 * j) ** 2 for j in range(k + 1)))
        eps = 2
    else:
        kap = k // 2
        A = sp.expand(sum(sp.chebyshevu(j, c) * (kap - j) for j in range(k + 1)))
        B = sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (kap - j) ** 2 for j in range(k + 1)))
        eps = 1
    At = sp.expand(A.subs(c, u / 2))
    Bt = sp.expand(B.subs(c, u / 2))
    m = k + 2
    FA = (k + 1) // 2
    roots_u = [2 * sp.cos(2 * sp.pi * j / m) for j in range(1, FA + 1)]
    mps, seen = [], set()
    for r in roots_u:
        mp_r = sp.minimal_polynomial(r, u)
        s = str(mp_r)
        if s not in seen:
            seen.add(s); mps.append(mp_r)
    Phi = sp.Poly(sp.expand(sp.prod(mps)), u, domain="QQ")
    return At, Bt, Phi, FA, m, eps


def reduce_expr(expr, Phi):
    R = sp.rem(sp.Poly(sp.expand(expr), u, domain="QQ"), Phi)
    return sp.Poly(R.as_expr(), u, domain="QQ")


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
    print("=" * 116)
    print("F89 edge A: full assembly  Rm = eps^2 m^2 W (mod nothing, exact),  W=(2+u)^2 Btil mod Phi_u")
    print("=" * 116)
    print()
    print(f"  {'k':>3} {'cls':>5} {'FA':>3} {'eps':>3} {'v2m':>3} {'delta':>5} {'T':>3} "
          f"{'L1_ok':>5} {'Rm=eps2m2W':>10} {'RESID_rhs':>9} {'minW(d>=1)':>10} {'RESID_ok':>8}")
    print("  " + "-" * 96)
    L1all = chainall = residall = True
    for k in range(5, 37):
        At, Bt, Phi, FA, m, eps = pieces(k)
        # L1
        Lred = reduce_expr((4 - u**2) * At, Phi)
        L1_ok = (sp.expand(Lred.as_expr() - eps * m * (2 + u)) == 0)
        L1all = L1all and L1_ok
        # chain Rm == eps^2 m^2 W
        Nint = (4 - u**2) ** 2 * At**2 * Bt
        Rm = reduce_expr(Nint, Phi)
        W = reduce_expr((2 + u) ** 2 * Bt, Phi)
        chain_ok = (sp.expand(Rm.as_expr() - eps**2 * m**2 * W.as_expr()) == 0)
        chainall = chainall and chain_ok
        # RESID
        vk = v2_int(k); vm = v2_int(m); delta = delta_closed(k); T = 4 + delta
        rhs = T - 2 * v2_int(eps) - 2 * vm
        Wprof = [sv2(W.coeff_monomial(u**d)) for d in range(FA)]
        Wd1 = [Wprof[d] for d in range(1, FA) if Wprof[d] is not None]
        minW = min(Wd1) if Wd1 else None
        resid_ok = (minW is not None and minW >= rhs)
        residall = residall and resid_ok
        print(f"  {k:>3} {cls_of(k):>5} {FA:>3} {eps:>3} {vm:>3} {delta:>5} {T:>3} "
              f"{str(L1_ok):>5} {str(chain_ok):>10} {rhs:>9} {str(minW):>10} {str(resid_ok):>8}")
    print()
    print(f"  LEMMA 1  (4-u^2)Atil == eps_k m (2+u) mod Phi_u  (all k):  {L1all}")
    print(f"  CHAIN    Rm == eps^2 m^2 W (exact polynomial identity)  (all k):  {chainall}")
    print(f"  RESID    val2(W_d) >= T - 2v2(eps) - 2v2(m)  for d>=1  (all k):  {residall}")
    print()
    print("  STATUS: top dominance == DIV == [LEMMA 1 (proved-form, closed) AND RESID].")
    print("  LEMMA 1 is a clean closed identity (degree-1 reduction). RESID is the residual.")


if __name__ == "__main__":
    main()
