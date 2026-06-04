#!/usr/bin/env python3
"""F89 edge A SUMMARY: complete status of the top-degree-dominance proof (2026-06-04).

TOP-DEGREE DOMINANCE (the last open input of the F89 cap proof):
   max_d cost_d  is attained at d = FA-1,   cost_d = -val2(Q_d),  Q = p_n(y) mod Phi_y.

This script is the single end-to-end verification + status report.  It (i) confirms dominance
bit-exact from the ground-truth y-pipeline for k=5..40 across all three residue classes, and
(ii) confirms every link of the rigorous reduction chain, marking which classes are CLOSED and
where the single residual sits.

PROOF (established this session; full details in _f89_edgeA_*.py, each link bit-exact):

  Reduction to the integer u-world (RIGOROUS):
    cost_d = 2d - val2(pref) - val2(R_d)  [residue equality];  s2 = 4 exactly (only denom is 16
    from (1-u^2/4)^2); Nint(u) = (4-u^2)^2 Atil(u)^2 Btil(u) with Atil=A(u/2),Btil=B(u/2) in Z[u].
    Top dominance <=> DIV: val2(Rm_d) >= T:=4+delta for d>=1, Rm = Nint mod Phi_u, plus the d=0
    corner (separate, satisfied).

  LEMMA 1 (RIGOROUS, node identity):  (4-u^2)Atil == eps m (2+u) (mod Phi_u), eps=2(odd)/1(even),
    from 4 sin th * S1(th) = eps m (2+2cos th), S1 = m cot(pi j/m).
    => Rm = eps^2 m^2 W,  W := (2+u)^2 Btil mod Phi_u; DIV <=> RESID: val2(W_d) >= rhs (d>=1),
    rhs = 3 (odd, 4|k) / 2 (2mod4).

  LEMMA 2 (RIGOROUS, node identity):  (4-u^2)^2 Btil == b0 + b2 u^2 (mod Phi_u) [degree 2],
    from 16 sin^2 th * S2(th).  Bridge: u^2 W == (4-u^2)^2 Btil == b0+b2 u^2 (mod (4,Phi_u))
    [since (2-u)^2 == u^2 (mod 4)].

  CLOSURE odd k & 4|k (RIGOROUS, COMPLETE):  Phi_u(0) = +/-1 (odd) => u is a 2-adic unit mod Phi_u
    => cancel u^2:  W == [u^2]^{-1}(b0+b2 u^2) (mod (4,Phi_u)).  With v2(b0)=3, v2(b2)=1 and
    u^2 [(2-u)^2]^{-1} == 1 (mod (4,Phi_u)), every d>=1 coeff of W == 0 (mod 8) => val2(W_d)>=3.
    => RESID, hence TOP DOMINANCE, PROVED for all odd k and all 4|k.

  RESIDUAL k = 2 mod4 (LOCALIZED GAP):  here Phi_u(0)=0 (u a zero-divisor), v2(b0)>=3, v2(b2)>=2 so
    u^2 W == b0+b2 u^2 == 0 (mod (4,Phi_u)).  Need to deduce val2(W_d)>=2 (d>=1).  The cancellation
    of u^2 is NOT available (u is a zero divisor; Phi_u = u*Psi with Psi(0) sometimes even).
    Verified bit-exact for k=6..42 that RESID holds, but a clean derivation of W==0 (mod(4,Phi))
    from u^2 W==0 (mod(4,Phi)) in this non-unit class is the one residual item.
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
y = sp.Symbol("y")


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
    return sp.Poly(sp.rem(sp.Poly(sp.expand(expr), u, domain="QQ"), Phi).as_expr(), u, domain="QQ")


def delta_closed(k):
    if k % 2 == 1:
        return 1
    if k % 4 == 0:
        return 1
    return 2 * v2_int(k // 2 + 1)


def cost_profile_y(k):
    p_y = fpk.derive_p_n_in_y(k)
    mp_y, _ = fpk.get_orbit_polynomial(k)
    R = sp.rem(sp.Poly(sp.expand(p_y), y, domain="QQ"), sp.Poly(sp.expand(mp_y), y, domain="QQ"))
    deg = R.degree()
    out = []
    for d in range(deg + 1):
        v = sv2(sp.Rational(R.coeff_monomial(y**d)))
        out.append(None if v is None else -v)
    return out


def cls_of(k):
    if k % 2 == 1:
        return "odd"
    if k % 4 == 0:
        return "4|k"
    return "2mod4"


def main():
    print("=" * 122)
    print("F89 TOP-DEGREE DOMINANCE -- complete status, end-to-end verification k=5..40 (all classes)")
    print("=" * 122)
    print()
    print(f"  {'k':>3} {'cls':>5} {'FA':>3} | {'v2Dk(=maxcost)':>13} {'cost_argmax':>11} {'DOM':>4} | "
          f"{'L1':>3} {'L2':>3} {'Phi(2)':>6} {'Phi(0)':>6} {'RESID':>5} {'class CLOSED?':>13}")
    print("  " + "-" * 110)
    domall = True
    closed_oddfour = True
    resid_2mod4_holds = True
    for k in range(5, 41):
        At, Bt, Phi, FA, m, eps = pieces(k)
        # ground truth dominance
        cy = cost_profile_y(k)
        present = [d for d in range(len(cy)) if cy[d] is not None]
        top = len(cy) - 1
        maxc = max(cy[d] for d in present)
        dom = (cy[top] == maxc)
        argmax = max(d for d in present if cy[d] == maxc)
        domall &= dom
        # chain links
        Lred = reduce_expr((4 - u**2) * At, Phi)
        L1 = (sp.expand(Lred.as_expr() - eps * m * (2 + u)) == 0)
        Q = reduce_expr((4 - u**2) ** 2 * Bt, Phi)
        L2 = (Q.degree() <= 2)
        Phi2 = int(Phi.as_expr().subs(u, 2))
        Phi0 = int(Phi.as_expr().subs(u, 0))
        W = reduce_expr((2 + u) ** 2 * Bt, Phi)
        rhs = 3 if (k % 2 == 1 or k % 4 == 0) else 2
        Wd1 = [sv2(W.coeff_monomial(u**d)) for d in range(1, FA)]
        Wd1 = [x for x in Wd1 if x is not None]
        resid = (min(Wd1) >= rhs) if Wd1 else True
        cl = cls_of(k)
        if cl in ("odd", "4|k"):
            cls_closed = "YES(proved)"
            closed_oddfour &= (L1 and L2 and (Phi2 % 2 == 1) and resid)
        else:
            cls_closed = "gap(2mod4)"
            resid_2mod4_holds &= resid
        print(f"  {k:>3} {cl:>5} {FA:>3} | {maxc:>13} {argmax:>11} {str(dom):>4} | "
              f"{str(L1):>3} {str(L2):>3} {Phi2:>6} {Phi0:>6} {str(resid):>5} {cls_closed:>13}")
    print()
    print(f"  END-TO-END dominance (y-pipeline ground truth), k=5..40:  {domall}")
    print(f"  odd k & 4|k FULLY PROVED (L1+L2+unit+RESID all hold):     {closed_oddfour}")
    print(f"  k=2mod4 RESID holds bit-exact (closure argument is the gap): {resid_2mod4_holds}")
    print()
    print("  ===========================================================================================")
    print("  STATUS:  CLOSED for odd k and 4|k (rigorous).")
    print("           CLOSED-MODULO-one-residual for k = 2 mod4:")
    print("             from u^2 W == 0 (mod(4,Phi_u)) deduce W == 0 (mod(4,Phi_u)) [val2(W_d)>=2, d>=1],")
    print("             where u is a 2-adic ZERO-DIVISOR (Phi_u(0)=0).  Verified bit-exact k=6..40.")
    print("  ===========================================================================================")


if __name__ == "__main__":
    main()
