#!/usr/bin/env python3
"""Derive the windowed F87 obstruction-size law max = min(2W-1, 2k-3) (scout, 2026-06-04).

A k-body diagonal-cell term's X/Y mask becomes a GF(2) polynomial p (bit j <-> x^j); placing it at
window w is multiplication by x^w. A pair's odd 𝔽₂-relation is therefore
    sum_{w in A} x^w p1  +  sum_{w in B} x^w p2 = 0   in GF(2)[x],
i.e. q_A p1 = q_B p2 with q_A = sum_{w in A} x^w (degree <= W-1, popcount = |A|), and the cycle size
is |A|+|B| = popcount(q_A)+popcount(q_B). With g = gcd(p1,p2), p1 = g*a, p2 = g*b (gcd(a,b)=1), the
relation is q_A a = q_B b, forcing q_A = b*s, q_B = a*s; the minimal (s=1) gives
    minimal cycle = popcount(p1/g) + popcount(p2/g)     (when its parity is odd, i.e. v1 != v2),
available once W-1 >= deg of the quotients (saturation). This scout checks (1) that formula against
the actual minimal odd cycle from the mask shifts, and (2) that the saturated max over pairs is
2k-3 = (k-1)+(k-2). Every diagonal-cell term mask has EVEN popcount, so (x+1) | p (p(1)=0).
"""
from __future__ import annotations

import sys
from itertools import product, combinations, combinations_with_replacement
from pathlib import Path

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw

DIAG = {"I", "Z"}


# ---- GF(2)[x] arithmetic on int-encoded polynomials (bit j = coeff of x^j) ----
def deg(p):
    return p.bit_length() - 1


def gf2_mul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
    return r


def gf2_divmod(a, b):
    """a = q*b + r over GF(2)[x]."""
    q = 0
    db = deg(b)
    while a and deg(a) >= db:
        shift = deg(a) - db
        q ^= 1 << shift
        a ^= b << shift
    return q, a


def gf2_gcd(a, b):
    while b:
        _, r = gf2_divmod(a, b)
        a, b = b, r
    return a


def popcount(x):
    return bin(x).count("1")


def min_odd_cycle(masks):
    masks = sorted(set(masks))
    for size in range(3, len(masks) + 1, 2):
        for combo in combinations(masks, size):
            x = 0
            for m in combo:
                x ^= m
            if x == 0:
                return size
    return 0


def main():
    print("=" * 96)
    print("F87 obstruction-size law: minimal odd cycle = popcount(p1/g)+popcount(p2/g), max = 2k-3")
    print("=" * 96)
    for k in (4, 5, 6):
        N = 2 * k          # plenty of windows (W = N-k+1 = k+1 >= k-1): saturated regime
        W = N - k + 1
        terms = [t for t in product("IXYZ", repeat=k)
                 if not all(L == "I" for L in t) and fw.klein_index(t) == (0, 1)
                 and any(L not in DIAG for L in t)]

        def mask(t):
            m = 0
            for i, L in enumerate(t):
                if L in ("X", "Y"):
                    m |= 1 << i
            return m

        formula_ok = True
        max_cyc = 0
        worst = None
        for t1, t2 in combinations_with_replacement(terms, 2):
            if (sum(c == "Y" for c in t1) % 2) != (sum(c == "Y" for c in t2) % 2):
                continue
            p1, p2 = mask(t1), mask(t2)
            S = {p1 << w for w in range(W)} | {p2 << w for w in range(W)}
            actual = min_odd_cycle(S)
            if actual == 0:
                continue  # soft
            # gcd-formula prediction (the v1 != v2 / odd branch)
            g = gf2_gcd(p1, p2)
            a = gf2_divmod(p1, g)[0]
            b = gf2_divmod(p2, g)[0]
            pred = popcount(a) + popcount(b)
            pred_is_odd = (pred % 2 == 1)
            # formula holds for the odd branch; for the even branch (v1==v2) the minimal ODD cycle
            # is larger, so we only assert the formula on the odd-branch pairs.
            if pred_is_odd and pred != actual:
                formula_ok = False
            if actual > max_cyc:
                max_cyc, worst = actual, (("".join(t1), "".join(t2)), p1, p2, g, a, b)
        print(f"\n  k={k}  (N={N}, W={W} windows, saturated):")
        print(f"    gcd-formula popcount(p1/g)+popcount(p2/g) == minimal odd cycle on odd-branch pairs: {formula_ok}")
        print(f"    max minimal-odd-cycle over hard pairs = {max_cyc}   (2k-3 = {2*k-3})   "
              f"{'MATCH' if max_cyc == 2*k-3 else 'MISMATCH'}")
        (lbl, p1, p2, g, a, b) = worst
        print(f"    extremal pair {lbl}: p1/g popcount {popcount(a)} (deg {deg(a)}), "
              f"p2/g popcount {popcount(b)} (deg {deg(b)}), gcd deg {deg(g)}; "
              f"sum {popcount(a)+popcount(b)} = (k-1)+(k-2)")
    print()
    print("  reading: if the formula holds and max=2k-3=(k-1)+(k-2), the saturated bound is derived:")
    print("  p1/g, p2/g have degree <= k-2 (g >= x+1 since masks are even-popcount), so popcount <= k-1;")
    print("  a coprime pair with odd popcount-sum maxes at (k-1)+(k-2)=2k-3. The window bound 2W-1 is")
    print("  |S| <= 2W (an odd subset has size <= 2W-1). So max = min(2W-1, 2k-3).")


if __name__ == "__main__":
    main()
