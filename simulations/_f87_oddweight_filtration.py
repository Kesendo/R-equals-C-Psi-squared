#!/usr/bin/env python3
"""F115 / "zero is the mirror" probe: the d=0 obstruction as an ODD-WEIGHT enumeration of the quasi-cyclic
code <(a,b)>, read through the (1+x)-filtration (x=1 = the mirror / DC / Perron point). 2026-06-08.

Per the d-reduction (verified in _f87_middle_distribution.py) the per-size distribution collapses to the
d=0 layer. For a d=0 hard pair (p1,p2): g=gcd=(1+x)^m, the reduced generators a=p2/g, b=p1/g are COPRIME,
and the windowed obstruction = the min ODD weight of {(a*s, b*s) : deg(s) <= W-1-max(deg a,deg b)}; the s=1
term is the gcd-formula size popcount(a)+popcount(b), and s!=1 is the cancellation. Because a,b are coprime
AND the pair is hard, {v(a), v(b)} = {0, Delta} with Delta = |v(p1)-v(p2)| >= 1 = the (1+x)-valuation GAP,
the two masks' "distance from the mirror" (Tom's lens: the zero/mirror mode is x=1, the (1+x) prime).

Parts:
  A. validate the s-multiplier (QC-code) view == the direct min-odd-dependence over the windowed shifts;
  B. bucket d=0 hard pairs by Delta; obstruction-size distribution per Delta, across W;
  C. cancellation (gcd-formula size vs actual obstruction) per Delta -- where does s!=1 help, and is the
     ODD-weight condition the clean x=1 statement a(1)+b(1)=1 & s(1)=1 ?
"""
from __future__ import annotations
from itertools import combinations
from collections import Counter, defaultdict


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

def g_rest_degree(p1, p2):
    g = gf2_gcd(p1, p2)
    while g != 1 and gf2_divmod(g, ONEPX)[1] == 0:
        g = gf2_divmod(g, ONEPX)[0]
    return deg(g) if g > 1 else 0

def even_masks(k):
    return [m for m in range(1, 1 << k) if popcount(m) % 2 == 0]

def min_odd_dependence(masks):
    masks = sorted(set(masks))
    n = len(masks)
    for size in range(3, n + 1, 2):
        for combo in combinations(masks, size):
            x = 0
            for m in combo: x ^= m
            if x == 0:
                return size
    return 0

def min_odd_weight_via_s(a, b, W):
    """min over s (deg(s) <= W-1-max(deg a,deg b)) of popcount(a*s)+popcount(b*s), odd only; 0 if none."""
    cap = W - 1 - max(deg(a), deg(b))
    if cap < 0:
        return 0
    best = 0
    for s in range(1, 1 << (cap + 1)):
        w = popcount(gf2_mul(a, s)) + popcount(gf2_mul(b, s))
        if w % 2 == 1 and (best == 0 or w < best):
            best = w
    return best


def d0_hard_pairs(k):
    """(p1,p2,a,b,Delta) for d=0 hard pairs at body k."""
    out = []
    masks = even_masks(k)
    for p1, p2 in combinations(masks, 2):
        if valuation(p1) == valuation(p2):
            continue
        if g_rest_degree(p1, p2) != 0:
            continue
        g = gf2_gcd(p1, p2)
        a = gf2_divmod(p2, g)[0]
        b = gf2_divmod(p1, g)[0]
        Delta = abs(valuation(p1) - valuation(p2))
        out.append((p1, p2, a, b, Delta))
    return out


def main():
    print("=" * 92)
    print("A. validate: s-multiplier QC view  ==  direct windowed min-odd-dependence  (d=0 hard pairs)")
    print("=" * 92)
    ok = True
    for k in (4, 5, 6):
        for W in (k - 1, k, 2 * k):
            N = W + k - 1
            mism = 0
            for (p1, p2, a, b, Delta) in d0_hard_pairs(k):
                S = {p1 << w for w in range(W)} | {p2 << w for w in range(W)}
                direct = min_odd_dependence(S)
                via_s = min_odd_weight_via_s(a, b, W)
                if direct != via_s:
                    mism += 1
            ok = ok and mism == 0
            print(f"  k={k} W={W} (N={N}): mismatches={mism}  {'OK' if mism == 0 else 'MISMATCH'}")
    print(f"  ==> s-multiplier view {'VALIDATED' if ok else 'BROKEN'}\n")

    print("=" * 92)
    print("B. d=0 obstruction-size distribution bucketed by Delta = |v(p1)-v(p2)| (the mirror gap)")
    print("=" * 92)
    for k in (4, 5, 6, 7):
        pairs = d0_hard_pairs(k)
        for W in (k, k + 1, 2 * k - 1, 2 * k):
            by_delta = defaultdict(Counter)
            for (p1, p2, a, b, Delta) in pairs:
                s = min_odd_weight_via_s(a, b, W)
                by_delta[Delta][s] += 1
            parts = []
            for Delta in sorted(by_delta):
                dist = ", ".join(f"{sz}:{c}" for sz, c in sorted(by_delta[Delta].items()))
                parts.append(f"Delta={Delta}: {{{dist}}}")
            print(f"  k={k} W={W}:  " + "   ".join(parts))
        print()

    print("=" * 92)
    print("C. cancellation per Delta: gcd-formula size (s=1) vs actual obstruction, at W=2k")
    print("   also: is the ODD condition the x=1 statement  a(1)+b(1)=1  (different parity at the mirror)?")
    print("=" * 92)
    for k in (5, 6, 7):
        W = 2 * k
        canc = Counter()           # (formula, actual) -> count, where actual < formula
        parity_ok = True
        per_delta_formula = defaultdict(Counter)
        per_delta_actual = defaultdict(Counter)
        for (p1, p2, a, b, Delta) in d0_hard_pairs(k):
            formula = popcount(a) + popcount(b)
            actual = min_odd_weight_via_s(a, b, W)
            per_delta_formula[Delta][formula] += 1
            per_delta_actual[Delta][actual] += 1
            if actual != 0 and actual < formula:
                canc[(formula, actual)] += 1
            # x=1 parity: a(1)=popcount(a)%2, b(1)=popcount(b)%2; ODD weight needs a(1)+b(1)=1
            a1 = popcount(a) % 2
            b1 = popcount(b) % 2
            if (a1 + b1) % 2 != 1:
                parity_ok = False
        print(f"  k={k} W={W}: a(1)+b(1)=1 for every d=0 hard pair? {parity_ok}")
        print(f"    cancellation events (formula->actual : count): "
              + (", ".join(f"{f}->{ac}:{c}" for (f, ac), c in sorted(canc.items())) or "none"))
        for Delta in sorted(per_delta_formula):
            fdist = dict(sorted(per_delta_formula[Delta].items()))
            adist = dict(sorted(per_delta_actual[Delta].items()))
            print(f"    Delta={Delta}: gcd-formula {fdist}  -->  actual {adist}")
        print()


if __name__ == "__main__":
    main()
