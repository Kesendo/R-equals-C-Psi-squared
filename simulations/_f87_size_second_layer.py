#!/usr/bin/env python3
"""Door 3 finding: the obstruction size has a SECOND valuation layer (2026-06-05).

(1+x) is the unique hardness prime (verified in _f87_beyond_x1_scout.py: no other irreducible governs
hard/soft). But the OTHER shared factors are not inert. Writing g = gcd(p1,p2) = (1+x)^m · g_rest with
g_rest coprime to (1+x) (the shared "other-frequency" content), the obstruction size obeys a refined
body bound:

  max obstruction over hard pairs with deg(g_rest) = d   =   2k - 3 - 2d.

Derivation: the gcd-formula size is popcount(p1/g) + popcount(p2/g); deg(p_i/g) <= (k-1) - deg(g),
deg(g) = m + d >= 1 + d, so popcount(p_i/g) <= k - 1 - d, and an odd sum of two such maxes at 2k-3-2d.
The §7.7 bound 2k-3 is the d=0 (coprime-apart-from-1+x) face. This script checks (1) the ACTUAL minimal
odd obstruction realises the layered cap, k=4,5,6, and (2) the gcd-formula cap = 2k-3-2d, k=4..10.
"""
from __future__ import annotations

import sys
from itertools import combinations
from collections import Counter

def deg(p): return p.bit_length() - 1
def popcount(x): return bin(x).count("1")

def gf2_divmod(a, b):
    q, db = 0, deg(b)
    while a and deg(a) >= db:
        sh = deg(a) - db; q ^= 1 << sh; a ^= b << sh
    return q, a

def gf2_gcd(a, b):
    while b: a, b = b, gf2_divmod(a, b)[1]
    return a

ONEPX = 0b11

def g_rest_degree(p1, p2):
    """deg of the non-(1+x) part of gcd(p1,p2)."""
    g = gf2_gcd(p1, p2)
    while g != 1 and gf2_divmod(g, ONEPX)[1] == 0:
        g = gf2_divmod(g, ONEPX)[0]
    return deg(g) if g > 1 else 0

def valuation(p, phi=ONEPX):
    v = 0
    while p != 0 and gf2_divmod(p, phi)[1] == 0:
        p = gf2_divmod(p, phi)[0]; v += 1
    return v

def gcd_formula_size(p1, p2):
    g = gf2_gcd(p1, p2)
    return popcount(gf2_divmod(p1, g)[0]) + popcount(gf2_divmod(p2, g)[0])

def even_masks(k):
    return [m for m in range(1, 1 << k) if popcount(m) % 2 == 0]

def min_odd_dependence(masks):
    masks = sorted(set(masks))
    for size in range(3, len(masks) + 1, 2):
        for combo in combinations(masks, size):
            x = 0
            for m in combo: x ^= m
            if x == 0: return size
    return 0


def main():
    print("=" * 84)
    print("Door 3: the obstruction-size second layer  max(size | deg g_rest = d) = 2k-3-2d")
    print("=" * 84)

    ok = True
    print("\n  (1) ACTUAL minimal odd obstruction, layered cap (k=4,5,6, N=2k saturated):")
    for k in (4, 5, 6):
        N = 2 * k; W = N - k + 1
        masks = even_masks(k)
        cap = {}
        for p1, p2 in combinations(masks, 2):
            if valuation(p1) == valuation(p2): continue  # soft
            d = g_rest_degree(p1, p2)
            S = {p1 << w for w in range(W)} | {p2 << w for w in range(W)}
            s = min_odd_dependence(S)
            cap[d] = max(cap.get(d, 0), s)
        line = []
        for d in sorted(cap):
            pred = 2 * k - 3 - 2 * d
            hit = (cap[d] == pred)
            ok = ok and hit
            line.append(f"d={d}:max={cap[d]}(={pred}{'' if hit else ' MISMATCH'})")
        print(f"    k={k}: " + "  ".join(line))

    print("\n  (2) gcd-formula cap = 2k-3-2d (the derivable bound), k=4..10:")
    for k in range(4, 11):
        masks = even_masks(k)
        cap = {}
        for p1, p2 in combinations(masks, 2):
            if valuation(p1) == valuation(p2): continue
            d = g_rest_degree(p1, p2)
            cap[d] = max(cap.get(d, 0), gcd_formula_size(p1, p2))
        allhit = all(cap[d] == 2 * k - 3 - 2 * d for d in cap)
        ok = ok and allhit
        caps = ", ".join(f"d={d}:{cap[d]}" for d in sorted(cap))
        print(f"    k={k:>2}: 2k-3={2*k-3:>2}  layered caps [{caps}]  {'OK' if allhit else 'MISMATCH'}")

    print(f"\n  ALL second-layer checks: {'OK' if ok else 'SOME MISMATCH'}")
    print("  reading: (1+x) decides hard/soft; the shared non-(1+x) factor degree d shrinks the")
    print("  obstruction by 2 per degree. The size law is layered: max = min(2W-1, 2k-3-2d).")


if __name__ == "__main__":
    main()
