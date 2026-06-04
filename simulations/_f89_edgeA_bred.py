#!/usr/bin/env python3
"""F89 edge A, step 14: hunt the second clean reduction for the Btil block, to close RESID.

Established (step 13, bit-exact all k):
   LINEAR LEMMA:  (4-u^2)*Atil(u) == eps_k * m * (2+u)  (mod Phi_u),  eps_k = 2 (odd k), 1 (even k).
   => Nint mod Phi_u = [(4-u^2)Atil]^2 Btil mod Phi = eps_k^2 m^2 (2+u)^2 Btil mod Phi_u.
   => val2(Rm_d) = 2*v2(eps_k) + 2*v2(m) + val2( [(2+u)^2 Btil mod Phi_u]_d ).

So DIV reduces to a clean target on the SINGLE object  W(u) := (2+u)^2 * Btil(u)  mod Phi_u:
   need  val2(W_d) >= T - 2 v2(eps_k) - 2 v2(m)   for d>=1.
For odd k: eps=2 => 2v2(eps)=2, v2(m)=0, RHS = T-2 = 5-2 = 3.  (matches RESID minv2=3.)

We look for closed/low-degree reductions of the Btil block:
   (a)  (4-u^2) * Btil  mod Phi_u  -- degree? closed form?
   (b)  (2+u)   * Btil  mod Phi_u
   (c)  (4-u^2)^2 * Btil mod Phi_u
   (d)  Btil itself mod Phi_u (full degree FA-1, but maybe 2-structured)
and print their reduced forms / v2 profiles, seeking a 'linear lemma' analogue or a clean
2-adic generating structure. Also test the SE/SE-anti partner: maybe (1-c^2)Btil ~ S2-type.

If (a) (4-u^2)Btil == (low-deg, 2-structured) mod Phi_u with the right 2-content, RESID closes
and odd-k DIV is fully proven.
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
        B = sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (k - 2 * j) ** 2 for j in range(k + 1)))
    else:
        kap = k // 2
        B = sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (kap - j) ** 2 for j in range(k + 1)))
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
    return Bt, Phi, FA, m


def reduce_expr(expr, Phi):
    R = sp.rem(sp.Poly(sp.expand(expr), u, domain="QQ"), Phi)
    return sp.Poly(R.as_expr(), u, domain="QQ")


def profile(P, FA):
    return [sv2(P.coeff_monomial(u**d)) for d in range(FA)]


def main():
    print("=" * 116)
    print("F89 edge A: Btil-block reductions (seek the second linear lemma to close RESID)")
    print("=" * 116)
    print()
    for k in (5, 7, 9, 11, 13):
        Bt, Phi, FA, m = pieces(k)
        Wa = reduce_expr((4 - u**2) * Bt, Phi)       # (a)
        Wb = reduce_expr((2 + u) * Bt, Phi)          # (b)
        Wc = reduce_expr((4 - u**2) ** 2 * Bt, Phi)  # (c)
        Wd = reduce_expr(Bt, Phi)                    # (d)
        Wresid = reduce_expr((2 + u) ** 2 * Bt, Phi)
        print(f"  k={k} (odd) FA={FA} m={m}:")
        print(f"    (4-u^2)Btil mod Phi  deg={Wa.degree():>2}  v2prof={profile(Wa,FA)}")
        print(f"      expr: {str(sp.expand(Wa.as_expr()))[:90]}")
        print(f"    (2+u)Btil   mod Phi  deg={Wb.degree():>2}  v2prof={profile(Wb,FA)}")
        print(f"    (4-u^2)^2Bt mod Phi  deg={Wc.degree():>2}  v2prof={profile(Wc,FA)}")
        print(f"    Btil        mod Phi  deg={Wd.degree():>2}  v2prof={profile(Wd,FA)}")
        print(f"    (2+u)^2 Btil mod Phi (=RESID W) v2prof={profile(Wresid,FA)}  minv2(d>=1)="
              f"{min([x for x in profile(Wresid,FA)[1:] if x is not None])}")
        print()
    print("  Seeking: a factor*(Btil) that reduces to LOW degree (like (4-u^2)Atil did -> linear).")
    print("  If (4-u^2)Btil is NOT low-degree, the Btil block has no single linear collapse; then")
    print("  RESID must be argued differently (e.g. Btil == Atil-like-square structure, or a")
    print("  generating 2-adic bound). The v2 profiles reveal the residual 2-content to bound.")


if __name__ == "__main__":
    main()
