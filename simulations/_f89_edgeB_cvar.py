#!/usr/bin/env python3
"""F89 top-degree dominance, the c-variable reduced numerator L_c and its 2-structure (2026-06-04), B-18.

We work in the c-variable directly: Phi_c = monic-up-to-lead orbit poly of c_j=cos(2 pi j/m), and
L_c(c) = G(c) mod Phi_c, where G=(1-c^2)^2 A^2 B integer in c.

WARNING: in c, the nodes c_j=cos(...) are roots of the MINIMAL polynomial of 2cos scaled; 2c_j=u_j is
the algebraic integer. So Phi_c is NOT monic-integer; Phi_u (in u=2c) is. Hence the u-variable is the
right home (Washington). We nonetheless tabulate the c-route to see if the proof's stated
'G~ mod Phi' (c-variable) has a clean mod-2 law, since the doc's L~ closed form is c-variable.

KEY TEST: the doc says val2(L~)=(FA-1)+1 and the cap '+1'. We connect F1/F2 (u-variable) to the
c-variable by  R(u)=2^{s2} G(u/2) mod Phi_u.  We verify the dictionary:
   R_d (u-coeff)  vs  Lc-coeff, and confirm the u-route facts F1,F2,F3 are the faithful carriers.

Then we test the SHARPEST derivable lemma for F2 in u:
   define D16(u) = (R(u) - R(0)) / 16  and check it is an INTEGER polynomial (this IS F2).
   We then test whether 16 | (R(u)-R(0)) follows from  Nint(u) - Nint(0) being in the ideal
   (16, Phi(u)-Phi(0))... i.e. a congruence Nint(u) ≡ Nint(0) mod (16, u*unit). We probe the
   ideal membership numerically: is (R(u)-R(0)) in 16*Z[u]?  (yes = F2), AND is (Nint(u)-Nint(0))
   already ≡ 0 mod (16) before reduction? (to see if 16-divisibility is pre- or post-reduction).
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


def v2int(n):
    n = abs(int(n))
    if n == 0:
        return None
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def build(k):
    A = sp.expand(sum(sp.chebyshevu(j, c) * (k - 2 * j) for j in range(k + 1)))
    B = sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (k - 2 * j) ** 2 for j in range(k + 1)))
    Gc = sp.expand((1 - c**2) ** 2 * A**2 * B)
    Nu = sp.expand(Gc.subs(c, u / 2))
    pl = sp.Poly(Nu, u)
    dd = sp.Integer(1)
    for co in pl.all_coeffs():
        dd = sp.lcm(dd, sp.denom(sp.Rational(co)))
    s2 = v2int(int(dd)) or 0
    Nint = sp.Poly(sp.expand(Nu * (2 ** s2)), u)
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
    Phi = sp.Poly(sp.expand(sp.prod(mps)), u)
    Rm = sp.rem(Nint, Phi)
    rp = sp.Poly(Rm.as_expr(), u)
    Rd = {e: int(co) for (e,), co in rp.terms()}
    for d in range(FA):
        Rd.setdefault(d, 0)
    return Rd, FA, m, Phi, Nint, s2


def main():
    print("=" * 110)
    print("F89 top-degree dominance B-18: F2 as 16 | (R(u)-R(0)); pre- vs post-reduction 16-divisibility")
    print("=" * 110)
    print()
    ks = [k for k in range(5, 40) if (k % 2 == 1 or k % 4 == 0)]
    print(f"  {'k':>3} {'cls':>4} {'FA':>3} {'16|(R-R0)':>10} {'16|(Nint-Nint0)?':>16} "
          f"{'(Nint-Nint0)/(4-u^2) 16div?':>26}")
    print("  " + "-" * 70)
    allF2 = True
    all_pre = True
    all_factor = True
    for k in ks:
        Rd, FA, m, Phi, Nint, s2 = build(k)
        cls = "odd" if k % 2 else "4|k"
        # F2: 16 | (R - R0)
        F2 = all((Rd[d] % 16) == 0 for d in range(1, FA))
        # pre-reduction: is 16 | (Nint(u)-Nint(0))?
        N0 = int(Nint.eval(0))
        diffN = sp.Poly(Nint.as_expr() - N0, u)
        pre16 = all((int(co) % 16) == 0 for co in diffN.all_coeffs() if co != 0) if diffN.as_expr() != 0 else True
        # factor out (4-u^2): Nint - Nint0 is divisible by u^2 (even powers); check (Nint-Nint0)/(u^2)?
        # actually test divisibility of (Nint - Nint0) by (4-u^2) and the 16-content of the quotient
        q, r = sp.div(sp.Poly(Nint.as_expr() - N0, u), sp.Poly(4 - u**2, u))
        if r.as_expr() == 0:
            fac16 = all((int(co) % 16) == 0 for co in q.all_coeffs() if co != 0)
        else:
            fac16 = None
        allF2 = allF2 and F2
        all_pre = all_pre and pre16
        if fac16 is not None:
            all_factor = all_factor and fac16
        print(f"  {k:>3} {cls:>4} {FA:>3} {str(F2):>10} {str(pre16):>16} {str(fac16):>26}")
    print()
    print(f"  ===== F2 (16|(R-R0)) post-reduction: {allF2} =====")
    print(f"  ===== 16|(Nint-Nint0) PRE-reduction: {all_pre}  (if False, the 16 is born in reduction) =====")
    print()
    print("  Interpretation: if PRE-reduction 16-divisibility is FALSE but POST is TRUE, the divisibility")
    print("  by 16 is created by the orbit reduction (the Phi-multiple subtracted carries the 2-mass).")
    print("  That is the crux that resists a purely-Nint derivation; it needs Phi's 2-adic structure.")


if __name__ == "__main__":
    main()
