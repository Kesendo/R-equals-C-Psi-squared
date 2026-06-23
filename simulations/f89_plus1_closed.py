#!/usr/bin/env python3
"""F89 '+1' lemma CLOSED: explicit closed form for the leading trace L~, and the cap made
elementary on every residue class (2026-06-04).

The cap proof (PROOF_F89_PATH_D_CLOSED_FORM.md, two-route chain) reduced everything to one
underived input: val2(L~), where L~ is the leading coefficient of (numerator mod Phi_c), and
the '+1' is the claim that the leading reduction step adds exactly one factor of 2.

This file closes that val2(L~) part. Via the node identities of f89_plus1_scout.py
(G(c_j) = S1^2 S2; S1 = m cot for odd k, (m/2) cot for even k after the pull-out;
L~ = sum_j G(c_j)/Phi'(c_j), the codifferent trace), L~ has a closed form:

    odd  k:  L~ = 2^{FA-2} * m^2     * (m^2 + 3)            m = k+2,  FA = (k+1)//2
    even k:  L~ = 2^{FA-1} * (m/2)^2 * ((m/2)^2 + 1)

and the 2-adic valuation reads off ELEMENTARILY, with the cap on every class:

    odd k:        v2(L~) = (FA-2) + v2(m^2+3)      = (FA-2)+2 = FA = (FA-1)+1
                  because m odd => m^2 == 1 (mod 8) => m^2+3 == 4 (mod 8), v2 = 2 exactly.
    4|k:          v2(L~) = (FA-1) + v2((m/2)^2+1)  = (FA-1)+1
                  because m/2 odd => (m/2)^2+1 == 2 (mod 8), v2 = 1 exactly.
    k=2 mod 4:    v2(L~) = (FA-1) + v2((m/2)^2)    = (FA-1) + 2 v2(kappa+1)
                  because m/2 even => (m/2)^2+1 is ODD, v2 = 0; and m/2 = kappa+1.

The bracket parities (m^2+3 == 4 mod 8; (m/2)^2+1 == 2 mod 8 for m/2 odd; odd for m/2 even)
are rigorous for ALL k. The closed forms for L~ are verified bit-exact here over k=3..51 and
match the two-route chain's val2(L~) (and hence reproduce v2(D_k) = polydeg + a(k)).
"""
from __future__ import annotations

import sys
from mpmath import mp, mpf, cos, pi, fabs, nint

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

mp.dps = 90


def cheb_u(j, x):
    um1, u0 = mpf(0), mpf(1)
    if j == 0:
        return u0
    for _ in range(j):
        um1, u0 = u0, 2 * x * u0 - um1
    return u0


def v2_int(n: int) -> int:
    n = abs(int(n))
    if n == 0:
        return 0
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def Ltilde(k: int) -> int:
    """Leading reduced coefficient L~ = sum_j G(c_j)/Phi'(c_j), correct numerator per parity."""
    m = k + 2
    FA = (k + 1) // 2
    even = (k % 2 == 0)
    kap = k // 2
    nodes = [cos(2 * pi * j / m) for j in range(1, FA + 1)]
    tot = mpf(0)
    for idx in range(FA):
        c = nodes[idx]
        if even:
            A = sum(cheb_u(jj, c) * (kap - jj) for jj in range(k + 1))
            B = sum(cheb_u(jj, c) ** 2 * (kap - jj) ** 2 for jj in range(k + 1))
        else:
            A = sum(cheb_u(jj, c) * (k - 2 * jj) for jj in range(k + 1))
            B = sum(cheb_u(jj, c) ** 2 * (k - 2 * jj) ** 2 for jj in range(k + 1))
        Gc = (1 - c * c) ** 2 * A * A * B
        Phip = mpf(1)
        for l in range(FA):
            if l != idx:
                Phip *= (c - nodes[l])
        tot += Gc / Phip
    n = int(nint(tot))
    assert fabs(tot - n) < mpf(10) ** (-25), f"L~ not integer at k={k}: {tot}"
    return n


def closed_form(k: int) -> int:
    m = k + 2
    FA = (k + 1) // 2
    if k % 2 == 1:
        return 2 ** (FA - 2) * m * m * (m * m + 3)
    h = m // 2
    return 2 ** (FA - 1) * h * h * (h * h + 1)


def target_v2(k: int) -> int:
    """The per-class val2(L~) the two-route chain needs (= cap a(k) restated as 2(FA-1)-v2(D_k)... )."""
    FA = (k + 1) // 2
    if k % 2 == 1:
        return FA  # = (FA-1)+1
    if k % 4 == 0:
        return (FA - 1) + 1
    return (FA - 1) + 2 * v2_int(k // 2 + 1)  # kappa+1


def main():
    print("=" * 100)
    print("F89 '+1' lemma CLOSED: L~ closed form + elementary 2-adic valuation, all residue classes")
    print("=" * 100)
    print()
    print(f"  {'k':>4} {'m':>4} {'FA':>4} {'class':>9} {'L~ (trace)':>22} {'closed form':>22} "
          f"{'eq':>4} {'v2(L~)':>7} {'target':>7} {'ok':>4}")
    all_ok = True
    for k in range(3, 52):
        L = Ltilde(k)
        cf = closed_form(k)
        m, FA = k + 2, (k + 1) // 2
        cls = "odd" if k % 2 else ("4|k" if k % 4 == 0 else "k=2mod4")
        tv = target_v2(k)
        ok = (L == cf) and (v2_int(L) == tv)
        all_ok = all_ok and ok
        print(f"  {k:>4} {m:>4} {FA:>4} {cls:>9} {L:>22} {cf:>22} {str(L == cf):>4} "
              f"{v2_int(L):>7} {tv:>7} {'OK' if ok else 'FAIL':>4}")

    # the elementary bracket-parity facts, checked for ALL m in range (rigorous, not a fit)
    print()
    bad = []
    for k in range(3, 200):
        m = k + 2
        if k % 2 == 1:
            if v2_int(m * m + 3) != 2:
                bad.append(("odd", k))
        else:
            h = m // 2
            want = 1 if (h % 2 == 1) else 0
            if v2_int(h * h + 1) != want:
                bad.append(("even", k))
    print(f"  elementary brackets (m²+3≡4 mod8; (m/2)²+1≡2 mod8 if m/2 odd, odd if even) for k=3..199: "
          f"{'ALL HOLD' if not bad else f'FAIL {bad[:5]}'}")
    print()
    print(f"  ===== closed forms bit-exact + elementary v2 = chain target, k=3..51: {all_ok} =====")
    print()
    print("  The '+1' (val2(L~) part of the cap lemma) is therefore CLOSED modulo the two L~ closed forms,")
    print("  whose 2-adic content is elementary and rigorous for all k. Remaining: top-degree dominance only.")
    assert all_ok and not bad


if __name__ == "__main__":
    main()
