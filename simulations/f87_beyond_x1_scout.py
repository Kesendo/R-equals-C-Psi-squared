#!/usr/bin/env python3
"""Door 3: is there structure beyond the (1+x) valuation? (scout, 2026-06-05).

Hardness is governed entirely by v_{1+x} (the x=1 / DC / Perron point). But a mask is a GF(2)[x]
polynomial with a FULL factorization into irreducibles (1+x, 1+x+x², 1+x²+x³, 1+x+x³, …), each an
"other frequency" (a root of unity in an extension field). This scout asks, empirically and without
fishing, two grounded questions:

  (a) does the valuation difference at ANOTHER irreducible φ ≠ (1+x) ever govern hard/soft?
      Expectation NO — that would confirm (1+x) is the unique hardness prime.
  (b) what, if anything, do the non-(1+x) factors index (the obstruction SIZE within hard pairs?
      a different pair property? nothing)?

If (a) is a clean NO and (b) finds the non-(1+x) factors govern the obstruction SIZE, then the
valuation picture is layered: (1+x) decides IF, the other factors shape HOW BIG. If (b) is null, the
whole story lives at (1+x) and Door 3 closes — itself a clean finding.
"""
from __future__ import annotations

import sys
from itertools import combinations
from collections import Counter

# ---- GF(2)[x] ----
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

def irreducibles(max_deg):
    """All monic irreducible polynomials over GF(2) of degree 1..max_deg (constant term forced 1
    except x itself; masks have constant term 0 or not, but factors of interest are the cyclotomic
    ones with nonzero constant term, plus x)."""
    irr = []
    for d in range(1, max_deg + 1):
        for p in range(1 << d, 1 << (d + 1)):  # degree exactly d
            if all(gf2_divmod(p, q)[1] != 0 for q in irr if deg(q) <= d // 2 + 1):
                # trial division by all smaller-degree polys (not just irr) for correctness
                if all(gf2_divmod(p, c)[1] != 0 for c in range(2, 1 << (d // 2 + 1)) if 1 <= deg(c) <= d // 2):
                    irr.append(p)
    return irr

def factor(p, irr):
    """Factor p over GF(2)[x] as {irreducible: multiplicity}."""
    fac = Counter()
    # strip x = 0b10
    while p % 2 == 0 and p != 0:
        fac[0b10] += 1; p >>= 1
    for q in irr:
        if q == 0b10: continue
        while p != 1 and gf2_divmod(p, q)[1] == 0:
            fac[q] += 1; p = gf2_divmod(p, q)[0]
    return fac, p  # p = remaining cofactor (1 if fully factored over irr)

def valuation(p, phi):
    v = 0
    while p != 0 and gf2_divmod(p, phi)[1] == 0:
        p = gf2_divmod(p, phi)[0]; v += 1
    return v

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
    print("=" * 92)
    print("Door 3: structure beyond the (1+x) valuation?")
    print("=" * 92)

    irr = irreducibles(5)
    print(f"\n  irreducibles up to deg 5: {[format(q, 'b') for q in irr]}")
    ONEPX = 0b11  # 1+x

    for k in (4, 5, 6):
        N = 2 * k; W = N - k + 1
        masks = even_masks(k)
        # (a) does v_phi difference at phi != (1+x) ever govern hard/soft?
        # hard truth = v_{1+x} differs; test each other irreducible as a would-be predictor.
        phis = [q for q in irr if q != ONEPX and q != 0b10]
        predictor_hits = {q: 0 for q in phis}        # times "v_phi differs" == "hard"
        total = 0
        # (b) does the non-(1+x) gcd structure predict the obstruction SIZE?
        size_by_cofactor = {}                          # (non-1+x gcd cofactor deg) -> set of sizes
        for p1, p2 in combinations(masks, 2):
            hard = valuation(p1, ONEPX) != valuation(p2, ONEPX)
            total += 1
            for q in phis:
                if (valuation(p1, q) != valuation(p2, q)) == hard:
                    predictor_hits[q] += 1
            if hard:
                g = gf2_gcd(p1, p2)
                # strip the (1+x) part of g, keep the rest (the "other-frequency" shared content)
                g_rest = g
                while gf2_divmod(g_rest, ONEPX)[1] == 0 and g_rest != 1:
                    g_rest = gf2_divmod(g_rest, ONEPX)[0]
                S = {p1 << w for w in range(W)} | {p2 << w for w in range(W)}
                size = min_odd_dependence(S)
                key = deg(g_rest) if g_rest > 1 else 0
                size_by_cofactor.setdefault(key, Counter())[size] += 1

        print(f"\n  k={k} (N={N}):  {total} mask-pairs")
        print(f"    (a) only (1+x) governs hard/soft? other-irreducible predictors "
              f"(agreement-with-hard / total):")
        for q in phis:
            frac = predictor_hits[q] / total
            tag = "= (1+x) alone is the hardness prime" if abs(frac - 0.5) < 0.499 else "PERFECT PREDICTOR!"
            print(f"        phi={format(q,'b'):<7} agree={frac:.3f}  {tag if frac > 0.95 or frac < 0.05 else ''}")
        print(f"    (b) obstruction-size distribution by non-(1+x) shared-factor degree of gcd:")
        for key in sorted(size_by_cofactor):
            dist = dict(sorted(size_by_cofactor[key].items()))
            print(f"        non-(1+x) gcd cofactor deg {key}: sizes {dist}")

    print("\n  reading: (a) tells whether (1+x) is the unique hardness prime; (b) tells whether the")
    print("  other factors shape the obstruction SIZE (a second valuation layer) or are inert.")


if __name__ == "__main__":
    main()
