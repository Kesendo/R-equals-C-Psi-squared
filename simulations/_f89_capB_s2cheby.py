#!/usr/bin/env python3
"""F89 capB: prove the explicit pieces s2 and the Nint leading coeff (2026-06-04).

C1 (s2): Nu(u) = (1-(u/2)^2)^2 A(u/2)^2 B(u/2), c=u/2.  U_j(c) has degree j leading coeff 2^j.
A(c)=sum_j U_j(c)(k-2j) has degree k, leading coeff 2^k*(k-2k)=2^k*(-k).  Hmm (k-2j) at j=k
is k-2k=-k.  So lead(A)=2^k*(-k)... wait top degree j=k contributes U_k leading 2^k times
(k-2k)=-k => lead(A)=-k*2^k.  A(u/2): replace c=u/2, lead in u = -k*2^k*(1/2)^k = -k.  Integer!
Similarly B(c)=sum U_j^2 (k-2j)^2, degree 2k, lead = (2^k)^2 * k^2 = 4^k k^2; B(u/2) lead in u
= 4^k k^2 (1/2)^{2k} = k^2.  (1-(u/2)^2)^2 lead in u = ( -1/4 u^2)^2 = u^4/16, lead coeff 1/16.
So Nu lead coeff in u = (1/16)*( -k... wait A^2 lead = k^2) ... lead(Nu) = (1/16)*(k^2)*(k^2)=k^4/16.
=> s2 from leading term alone needs 2^4 to clear /16; if the WHOLE Nu has no worse denom, s2=4.
For even k: A(u/2),B(u/2) have integer-with-/2^? structure giving s2=0 because the pull-out /
parity makes lower coeffs integer and /16 is cancelled by k^4 even part... we MEASURE & confirm
the LEADING-term reasoning vs actual s2.

We verify:
  lead_u(A(u/2)) = -k         (integer)
  lead_u(B(u/2)) = k^2        (integer)
  lead_u(Nu)     = k^4 / 16
  s2: odd k -> 4 (the 1/16 uncancelled since k odd), even k -> measure (k^4 even cancels part).
Also Nint leading coeff = 2^{s2} * k^4/16 = k^4 * 2^{s2-4}; for odd k = k^4 (odd), v2=0;
for even k (s2=0) = k^4/16 *... must be integer => k^4/16 integer iff 16|k^4 iff v2(k)>=1 (v2(k^4)=4v2k>=4
iff v2k>=1). For even k v2k>=1 so 4v2k>=4, k^4/16 integer with v2 = 4v2k-4. CHECK vs Nint lead.
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


def v2(n) -> int:
    n = abs(int(n))
    if n == 0:
        return 10**9
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def main():
    print("=" * 100)
    print("F89 capB: explicit pieces s2 and Nint structure (the given degree-growth part)")
    print("=" * 100)
    print()
    reps = [5, 7, 9, 11, 6, 10, 14, 12, 20, 8, 24, 16, 32]
    print(f"  {'k':>3} {'v2k':>3} {'lead_u(A/2)':>11} {'lead_u(B/2)':>11} "
          f"{'lead_u(Nu)':>12} {'s2':>3} {'predS2':>6} {'v2(Nint_lead)':>13} {'pred':>5} {'ok':>3}")
    print("  " + "-" * 84)
    okall = True
    for k in reps:
        A = sp.expand(sum(sp.chebyshevu(j, c) * (k - 2 * j) for j in range(k + 1)))
        B = sp.expand(sum(sp.chebyshevu(j, c) ** 2 * (k - 2 * j) ** 2 for j in range(k + 1)))
        Au = sp.expand(A.subs(c, u / 2))
        Bu = sp.expand(B.subs(c, u / 2))
        Nu = sp.expand((1 - (u / 2) ** 2) ** 2 * Au ** 2 * Bu)
        leadA = sp.Poly(Au, u).all_coeffs()[0]
        leadB = sp.Poly(Bu, u).all_coeffs()[0]
        leadNu = sp.Rational(sp.Poly(Nu, u).all_coeffs()[0])
        pl = sp.Poly(Nu, u)
        d = sp.Integer(1)
        for co in pl.all_coeffs():
            d = sp.lcm(d, sp.denom(sp.Rational(co)))
        s2 = v2(int(d))
        # predicted s2: odd k -> 4; even k -> 0 (claim)
        predS2 = 4 if k % 2 == 1 else 0
        Nint_lead = sp.Rational(leadNu) * (2 ** s2)
        v2Nl = v2(int(Nint_lead)) if Nint_lead == int(Nint_lead) else -99
        # predicted v2(Nint_lead): odd -> 0; even -> 4 v2k - 4 (from k^4/16)
        predNl = 0 if k % 2 == 1 else 4 * v2(k) - 4
        ok = (s2 == predS2 and v2Nl == predNl and leadA == -k and leadB == k * k)
        okall = okall and ok
        print(f"  {k:>3} {v2(k):>3} {str(leadA):>11} {str(leadB):>11} "
              f"{str(leadNu):>12} {s2:>3} {predS2:>6} {v2Nl:>13} {predNl:>5} {('OK' if ok else 'X'):>3}")
    print()
    print(f"  ALL (lead_u(A)=-k, lead_u(B)=k^2, s2 closed form, Nint_lead v2) bit-exact: {okall}")
    print()
    print("  Rigorous derivation of these:")
    print("   - U_j(c) leading coeff = 2^j (standard).  lead_c(A)=2^k*(k-2k)=-k*2^k; A(u/2) lead = -k.")
    print("   - lead_c(B)=4^k*k^2; B(u/2) lead = k^2.")
    print("   - (1-(u/2)^2)^2 lead = (u^2/4)^2 = u^4/16.  lead(Nu)=(1/16)*(-k)^2*(k^2)=k^4/16.")
    print("   - odd k: 16 has no common factor with k^4 (odd) => need 2^4 to clear => s2>=4;")
    print("     the rest of Nu has no worse 2-denominator (claim, verified) => s2=4 exactly.")
    print("   - even k: k^4 has v2=4v2k>=4 so k^4/16 is already an integer (v2=4v2k-4); and the")
    print("     full Nu is integer (verified) => s2=0.")
    print()


if __name__ == "__main__":
    main()
