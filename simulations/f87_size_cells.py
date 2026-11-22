#!/usr/bin/env python3
"""F115 MacWilliams kernel: the saturated per-SIZE distribution of the d=0 obstruction. 2026-06-08.

saturated_dist(k)[s] = Sum_D c(D,s) * max(0, k-1-D),
  c(D,s) = #{coprime reduced (a,b): v(a)=0, v(b)>=1, max(deg a,deg b)=D, minweight(a,b)=s}  (over all Delta),
  minweight(a,b) = min ODD weight of the convolutional code <(a,b)> = min_t odd[ wt(a t)+wt(b t) ]
                 = the saturated obstruction (s-decoupled; computed at cap = D+4, and the D = 7 figures
                   re-checked below at larger caps).

This script locates where the per-size distribution closes, in which layer, and where it turns
number-theoretic:

  FLOOR (size 3) -- CLOSED.  minweight==3 <=> the v=0 generator is a monomial x^j and the v=Delta generator
    has popcount 2 (so Delta is a power of two). Hence c(D,3)=3D-1 and the d0 size-3 total
    T3(k) = Sum_{D=1}^{k-2}(3D-1)(k-1-D) = (k-1)^2 (k-2)/2.

  MONOMIAL COLUMN -- CLOSED in the gcd-formula layer.  For a fixed even weight beta,
    #{(x^j, weight-beta b) coprime, max deg D} = 2*C(D,beta-1) + (D-1)*C(D-1,beta-2), a polynomial in D
    of degree beta-1: x^j is coprime to b exactly when j = 0 or b(0) = 1.

  CEILING (max size 2D+1) -- computed.  Count == 2 at D = 4..7: the repunit pair, R_D = 1+x+...+x^D with
    R_{D-1} or x*R_{D-1} (densest, no cancellation). D=2,3 are small edge cases (3, 4).

  THE HARD CORE -- LOCATED, not closed.  The irregularity enters exactly at the first popcount split with
    BOTH popcounts >= 2, namely (3,2) at size 5: a weight-3 a with a(0) = 1 is coprime to x^p(1+x^r) <=>
    coprime to 1+x^r, which depends on the factorization of 1+x^r (GF(2)[x] weighted coprimality). A
    polynomial count of total weight s has degree <= s-2 (coprimality puts a bit at degree 0, the maximal
    degree one at D); the order-(s-1) differences do not vanish, so the cell is not a polynomial in D. The
    convolutional cancellation (k>=6, the -2 cascade) sits on top.
"""
from __future__ import annotations
from itertools import combinations
from math import comb


def deg(p): return p.bit_length() - 1
def popcount(x): return bin(x).count("1")

def gf2_mul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
    return r

def gf2_divmod(a, b):
    q, db = 0, deg(b)
    while a and deg(a) >= db:
        sh = deg(a) - db; q ^= 1 << sh; a ^= b << sh
    return q, a

def gf2_gcd(a, b):
    while b: a, b = b, gf2_divmod(a, b)[1]
    return a

ONEPX = 0b11

def valuation(p):
    v = 0
    while p and gf2_divmod(p, ONEPX)[1] == 0:
        p = gf2_divmod(p, ONEPX)[0]; v += 1
    return v

def is_pow2(n): return n >= 1 and (n & (n - 1)) == 0
def is_monomial(p): return popcount(p) == 1
def repunit(D): return (1 << (D + 1)) - 1   # 1 + x + ... + x^D (all ones)


def minweight(a, b, cap=None):
    if cap is None:
        cap = max(deg(a), deg(b)) + 4
    best = 0
    for t in range(1, 1 << (cap + 1)):
        w = popcount(gf2_mul(a, t)) + popcount(gf2_mul(b, t))
        if w & 1 and (best == 0 or w < best):
            best = w
    return best


# ---- FLOOR: size 3 ----------------------------------------------------------
def floor_size3():
    print("=" * 96)
    print("FLOOR (size 3): characterization, c(D,3)=3D-1, and T3(k)=(k-1)^2(k-2)/2")
    print("=" * 96)
    # characterization both directions on a small grid
    bad = 0
    for Delta in (1, 2, 3, 4):
        us = [u for u in range(1, 1 << 6) if valuation(u) == 0]
        ws = [w for w in range(1, 1 << 6) if valuation(w) == Delta]
        for u in us:
            for w in ws:
                if gf2_gcd(u, w) != 1:
                    continue
                if (minweight(u, w) == 3) != (is_monomial(u) and popcount(w) == 2):
                    bad += 1
    print(f"  characterization minweight==3 <=> (mono u)&(pc(w)==2): {'HOLDS' if bad == 0 else f'FAILS({bad})'}")

    def cD3(D):
        c = 0
        for p in range(D + 1):
            for q in range(p + 1, D + 1):
                b = (1 << p) | (1 << q)
                for j in range(D + 1):
                    if (j == 0 or p == 0) and max(j, q) == D and gf2_gcd(1 << j, b) == 1:
                        c += 1
        return c
    cs = [cD3(D) for D in range(1, 9)]
    print(f"  c(D,3) for D=1..8: {cs}  (== 3D-1? {cs == [3 * D - 1 for D in range(1, 9)]})")
    ok = True
    for k in range(4, 13):
        t = sum((3 * D - 1) * max(0, k - 1 - D) for D in range(1, k - 1))
        closed = (k - 1) ** 2 * (k - 2) // 2
        ok = ok and t == closed
    print(f"  T3(k) = Sum(3D-1)(k-1-D) == (k-1)^2(k-2)/2 for k=4..12: {ok}\n")


# ---- MONOMIAL COLUMN --------------------------------------------------------
def monomial_column():
    print("=" * 96)
    print("MONOMIAL COLUMN: #{(x^j, weight-beta b) coprime, max deg D} = 2*C(D,beta-1) + (D-1)*C(D-1,beta-2)")
    print("=" * 96)

    def col(beta, D):
        c = 0
        Bs = [b for b in range(1, 1 << (D + 1)) if popcount(b) == beta and valuation(b) >= 1]
        for j in range(D + 1):
            a = 1 << j
            for b in Bs:
                if max(j, deg(b)) == D and gf2_gcd(a, b) == 1:
                    c += 1
        return c
    # x^j is coprime to b exactly when j = 0 or b(0) = 1: j = D and j = 0 each give C(D, beta-1),
    # each of the D-1 exponents in between gives C(D-1, beta-2).
    for beta in (2, 4, 6, 8):
        seq = [col(beta, D) for D in range(1, 14)]
        exact = [2 * comb(D, beta - 1) + (D - 1) * comb(D - 1, beta - 2) for D in range(1, 14)]
        assert seq == exact, f"monomial column beta={beta}: counted {seq}, closed form {exact}"
        print(f"  beta={beta}, D = 1..13: {seq}  == the closed form (a polynomial of degree {beta - 1} in D)")

    # The column counts by the generators' weights (the gcd-formula layer). In the saturated distribution a
    # monomial pair can still lose weight to a multiplier: b = 1+x+x^3+x^4+x^6+x^7 times 1+x+x^2 is 1+x^9.
    from collections import Counter
    D, beta = 7, 6
    pairs = [(1 << j, b) for j in range(D + 1) for b in range(1, 1 << (D + 1))
             if popcount(b) == beta and valuation(b) >= 1 and max(j, deg(b)) == D and gf2_gcd(1 << j, b) == 1]
    sat = dict(sorted(Counter(minweight(a, b) for a, b in pairs).items()))
    assert len(pairs) == 132 and sat == {5: 8, 7: 124}, (len(pairs), sat)
    sat14 = dict(sorted(Counter(minweight(a, b, cap=14) for a, b in pairs).items()))
    assert sat14 == sat, (sat14, sat)   # a larger cap finds no further drop
    print(f"  saturated, beta=6, D=7: {len(pairs)} monomial pairs, obstruction sizes {sat}")
    print()


# ---- CEILING: max size 2D+1 -------------------------------------------------
def ceiling():
    print("=" * 96)
    print("CEILING (max size 2D+1): count and the repunit pair")
    print("=" * 96)
    for D in range(2, 8):
        A = [a for a in range(1, 1 << (D + 1)) if valuation(a) == 0]
        Bs = [b for b in range(1, 1 << (D + 1)) if valuation(b) >= 1]
        tops = [(a, b) for a in A for b in Bs
                if max(deg(a), deg(b)) == D and gf2_gcd(a, b) == 1 and minweight(a, b) == 2 * D + 1]
        assert all(minweight(a, b, cap=D + 10) == 2 * D + 1 for a, b in tops)   # still 2D+1 at a larger cap
        rep = {repunit(D), repunit(D - 1), repunit(D - 1) << 1}
        on_repunits = all(a in rep and b in rep for a, b in tops)
        print(f"  D={D}: count={len(tops)}  pairs={[(f'{a:b}', f'{b:b}') for a, b in tops]}"
              f"  (repunit family? {on_repunits})")
    print()


# ---- HARD CORE: where polynomiality breaks ----------------------------------
def hard_core(Dmax=6):
    print("=" * 96)
    print("HARD CORE: gcd-formula layer per size sigma, and the (3,2) split = first irregular cell")
    print("=" * 96)

    def f_of_D(D):
        from collections import Counter
        out = Counter()
        A = [a for a in range(1, 1 << (D + 1)) if valuation(a) == 0]
        Bs = [b for b in range(1, 1 << (D + 1)) if valuation(b) >= 1]
        for a in A:
            for b in Bs:
                if max(deg(a), deg(b)) != D or gf2_gcd(a, b) != 1:
                    continue
                out[popcount(a) + popcount(b)] += 1
        return out
    rows = {D: f_of_D(D) for D in range(1, Dmax + 1)}

    # A count of coprime pairs of total weight s is O(D^(s-2)): coprimality puts a bit at degree 0 and the
    # maximal degree puts one at D. A polynomial count has degree <= s-2, so nonzero differences of order
    # s-1 rule out one polynomial from the first value on.
    def verdict(seq, s):
        if len(seq) < s:
            return f"too few points for order {s - 1}"
        d = list(seq)
        for _ in range(s - 1):
            d = [d[i + 1] - d[i] for i in range(len(d) - 1)]
        if any(d):
            return f"order-{s - 1} differences {d}: not a polynomial in D"
        return f"order-{s - 1} differences vanish: consistent with a polynomial of degree <= {s - 2}"
    for s in (3, 5, 7):
        ys = [rows[D].get(s, 0) for D in range(1, Dmax + 1)]
        ys = ys[next((i for i, y in enumerate(ys) if y > 0), 0):]
        print(f"  sigma={s}: {ys}  -> {verdict(ys, s)}")

    # the (3,2) split is the first 'both popcounts >= 2' cell, and it is the irregular one
    def split32(D):
        c = 0
        A = [a for a in range(1, 1 << (D + 1)) if valuation(a) == 0 and popcount(a) == 3]
        Bs = [b for b in range(1, 1 << (D + 1)) if valuation(b) >= 1 and popcount(b) == 2]
        for a in A:
            for b in Bs:
                if max(deg(a), deg(b)) == D and gf2_gcd(a, b) == 1:
                    c += 1
        return c
    seq = [split32(D) for D in range(2, Dmax + 5)]
    print(f"  (3,2) weight-3 x weight-2 split: {seq}  -> {verdict(seq, 5)}")
    print()


def main():
    floor_size3()
    monomial_column()
    ceiling()
    hard_core(10)   # gcd-formula layer is cheap; D up to 10 gives enough points for the poly test


if __name__ == "__main__":
    main()
