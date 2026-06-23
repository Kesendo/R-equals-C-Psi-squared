#!/usr/bin/env python3
"""Door 3 deepening (combinatorial): is the hard COUNT layered by d like the size CAP? (2026-06-05).

The total hard count is A203241 = e2(2^0..2^(k-2)). The obstruction-size CAP is layered by d=deg(g_rest)
(the shared non-(1+x) factor degree): max size = 2k-3-2d. Question: does the hard count ALSO split
cleanly by d, i.e. does #{hard mask-pairs with deg(g_rest)=d} have a closed form in (k,d)? If yes the
layering is total (count and size both clean by d); if not, only the cap is clean.
"""
from __future__ import annotations
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

def valuation(p, phi=ONEPX):
    v = 0
    while p != 0 and gf2_divmod(p, phi)[1] == 0:
        p = gf2_divmod(p, phi)[0]; v += 1
    return v

def g_rest_degree(p1, p2):
    g = gf2_gcd(p1, p2)
    while g != 1 and gf2_divmod(g, ONEPX)[1] == 0:
        g = gf2_divmod(g, ONEPX)[0]
    return deg(g) if g > 1 else 0

def even_masks(k):
    return [m for m in range(1, 1 << k) if popcount(m) % 2 == 0]


def main():
    print("=" * 78)
    print("Hard mask-pair count, split by d = deg(g_rest)")
    print("=" * 78)
    table = {}
    for k in range(3, 9):
        masks = even_masks(k)
        by_d = Counter()
        total = 0
        for p1, p2 in combinations(masks, 2):
            if valuation(p1) == valuation(p2):
                continue
            by_d[g_rest_degree(p1, p2)] += 1
            total += 1
        table[k] = by_d
        closed = (4 ** (k - 1) - 3 * 2 ** (k - 1) + 2) // 3
        row = "  ".join(f"d={d}:{by_d[d]}" for d in sorted(by_d))
        print(f"  k={k}: total={total} (A203241={closed}{'' if total == closed else ' !!'})   [{row}]")

    # CLOSED FORMS (verify, do not just fit):
    #   d=0 base B(k) = (4^k - 12k + 8) / 18   (from B(k) = 4 B(k-1) + 2(k-2), B(3)=2)
    #   #hard(deg g_rest = d, body k) = 2^(d-1) * B(k-d)   for d >= 1
    #   total = B(k) + sum_{d>=1} 2^(d-1) B(k-d) = A203241(k-1)
    def B(k): return (4 ** k - 12 * k + 8) // 18
    print("\n  closed forms (verified against the enumerated table):")
    base_ok = all(table[k][0] == B(k) for k in range(3, 9))
    print(f"    d=0 base B(k) = (4^k - 12k + 8)/18 :  {base_ok}   "
          f"({[B(k) for k in range(3,9)]})")
    layer_ok = all(table[k][d] == (1 << (d - 1)) * B(k - d)
                   for k in range(3, 9) for d in table[k] if d >= 1)
    print(f"    layer #hard(d,k) = 2^(d-1)*B(k-d) (d>=1) :  {layer_ok}")
    total_ok = all(
        sum(table[k].values()) == B(k) + sum((1 << (d - 1)) * B(k - d)
                                             for d in range(1, k - 2))
        for k in range(3, 9))
    print(f"    total = B(k) + sum_d 2^(d-1) B(k-d) = A203241 :  {total_ok}")
    print(f"\n    => the d-layering is TOTAL: both the size cap (2k-3-2d) and the COUNT are closed-form")
    print(f"       layered by the shared non-(1+x) factor degree. {'ALL OK' if base_ok and layer_ok and total_ok else 'MISMATCH'}")


if __name__ == "__main__":
    main()
