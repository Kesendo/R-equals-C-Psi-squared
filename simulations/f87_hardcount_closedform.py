#!/usr/bin/env python3
"""Verify the closed-form count of F87 hard pairs from the (1+x)-valuation classes (2026-06-05).

The coding-theory agent surfaced (and the §7.7 valuation criterion explains) a clean closed form for
HOW MANY diagonal-cell pairs are hard at body count k. Hardness ⟺ different (1+x)-adic valuation, so:

  - even-popcount nonzero k-bit masks split into valuation classes v = 1..k-1 of size c_v = 2^(k-1-v);
  - a pair is hard iff the two masks lie in DIFFERENT classes, so
        #hard mask-pairs = e₂(c_1,…,c_{k-1}) = e₂(2^0,…,2^(k-2)) = (4^(k-1) − 3·2^(k-1) + 2)/3   (OEIS A203241);
  - dressed by the 2^(2k-3) Klein / y-parity constant, this is the WindowedObstructionScan hard count.

This script checks every step bit-exact against direct enumeration, plus the agent's d=3 closed form and
the convolutional free-distance (constant 4) vs minimum-odd-weight (2k-3) distinction.
"""
from __future__ import annotations

import sys
from itertools import combinations
from math import comb

# ---- GF(2)[x] ----
def deg(p): return p.bit_length() - 1
def popcount(x): return bin(x).count("1")

def gf2_mul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
    return r

def gf2_divq(a, b):
    q, db = 0, deg(b)
    while a and deg(a) >= db:
        sh = deg(a) - db; q ^= 1 << sh; a ^= b << sh
    return q

def gf2_gcd(a, b):
    while b:
        db = deg(b); r = a
        while r and deg(r) >= db: r ^= b << (deg(r) - db)
        a, b = b, r
    return a

def valuation_1px(p):
    v = 0
    while p and popcount(p) % 2 == 0:
        p = gf2_divq(p, 0b11); v += 1
    return v

def even_masks(k):
    """Nonzero even-popcount k-bit masks (the diagonal-cell X/Y mask space)."""
    return [m for m in range(1, 1 << k) if popcount(m) % 2 == 0]

def min_odd_dependence(masks):
    masks = sorted(set(masks))
    for size in range(3, len(masks) + 1, 2):
        for combo in combinations(masks, size):
            x = 0
            for m in combo: x ^= m
            if x == 0: return size
    return 0

def conv_free_distance(a, b, max_s_deg=8):
    """Minimum weight (ANY parity) of the convolutional code {(a*s, b*s)}, s up to max_s_deg."""
    best = None
    for s in range(1, 1 << (max_s_deg + 1)):
        w = popcount(gf2_mul(a, s)) + popcount(gf2_mul(b, s))
        if w > 0 and (best is None or w < best): best = w
    return best


def main():
    print("=" * 92)
    print("Closed-form count of F87 hard pairs from the (1+x)-valuation classes")
    print("=" * 92)

    print(f"\n{'k':>2} {'class sizes c_v (v=1..k-1)':<34} {'=2^(k-1-v)?':<11} "
          f"{'#hard=e2':<9} {'closed':<8} {'2^(2k-3)*hard':<13}")
    ok_all = True
    for k in range(3, 11):
        masks = even_masks(k)
        # class sizes by valuation
        cls = {}
        for m in masks:
            cls[valuation_1px(m)] = cls.get(valuation_1px(m), 0) + 1
        cvec = [cls.get(v, 0) for v in range(1, k)]
        cv_pred = [2 ** (k - 1 - v) for v in range(1, k)]
        cv_ok = (cvec == cv_pred)
        # #hard mask-pairs = pairs in different classes = e2(class sizes)
        hard_direct = sum(1 for a, b in combinations(masks, 2)
                          if valuation_1px(a) != valuation_1px(b))
        e2 = sum(cvec[i] * cvec[j] for i in range(len(cvec)) for j in range(i + 1, len(cvec)))
        closed = (4 ** (k - 1) - 3 * 2 ** (k - 1) + 2) // 3
        ok = cv_ok and hard_direct == e2 == closed
        ok_all = ok_all and ok
        print(f"{k:>2} {str(cvec):<34} {str(cv_ok):<11} {hard_direct:<9} {closed:<8} "
              f"{2**(2*k-3) * closed:<13} {'OK' if ok else 'MISMATCH'}")

    # cross-check the dressed totals against the known WindowedObstructionScan hard counts
    print("\n  dressed total 2^(2k-3)*closed  vs  C# WindowedObstructionScan saturated hard count:")
    scan_hard = {4: 448, 5: 8960, 6: 158720}
    for k, h in scan_hard.items():
        closed = (4 ** (k - 1) - 3 * 2 ** (k - 1) + 2) // 3
        dressed = 2 ** (2 * k - 3) * closed
        print(f"    k={k}: closed={closed}  dressed={dressed}  scan={h}  {'MATCH' if dressed == h else 'MISMATCH'}")

    # the d=3 mask-class closed form (agent: 5*2^(k-1) - (3k^2+k)/2 - 3), checked at MASK level (N=2k)
    print("\n  d=3 obstruction class, MASK count vs closed form 5*2^(k-1) - (3k^2+k)/2 - 3 (N=2k):")
    for k in range(3, 8):
        N = 2 * k; W = N - k + 1
        masks = even_masks(k)
        cnt = 0
        for p1, p2 in combinations(masks, 2):
            if valuation_1px(p1) == valuation_1px(p2): continue
            S = {p1 << w for w in range(W)} | {p2 << w for w in range(W)}
            if min_odd_dependence(S) == 3: cnt += 1
        pred = 5 * 2 ** (k - 1) - (3 * k * k + k) // 2 - 3
        print(f"    k={k}: d=3 mask-count={cnt}  closed={pred}  {'MATCH' if cnt == pred else 'MISMATCH'}")

    # convolutional free distance (any parity) vs minimum odd weight, extremal family
    print("\n  extremal family: convolutional free distance (any parity) vs min ODD weight (= 2k-3):")
    for k in range(4, 9):
        p1 = (1 << 1) | (1 << (k - 1)); p2 = 1 | (1 << (k - 1))
        g = gf2_gcd(p1, p2); a = gf2_divq(p2, g); b = gf2_divq(p1, g)
        fd = conv_free_distance(a, b)
        odd = popcount(a) + popcount(b)
        print(f"    k={k}: a={a:b} b={b:b}  free_dist(any parity)={fd}  min_odd(s=1)={odd}=2k-3={2*k-3}")

    print(f"\n  ALL valuation-class / closed-form checks: {'OK' if ok_all else 'SOME MISMATCH'}")


if __name__ == "__main__":
    main()
