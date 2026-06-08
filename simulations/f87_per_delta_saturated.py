#!/usr/bin/env python3
"""F115: the Delta-bucket COUNT closed form (proven) + the saturated per-(Delta, size) distribution
(the open MacWilliams kernel). Builds on _f87_oddweight_filtration.py (2026-06-08).

A d=0 hard pair = {p1,p2} with gcd=(1+x)^m (pure power); reduced generators a=p2/g, b=p1/g are coprime
with {v(a),v(b)}={0,Delta}, Delta=|v(p1)-v(p2)| the (1+x)-valuation gap (distance from the mirror x=1).

CLOSED (proven here, verify_*):
  Delta-bucket count  n(Delta,k) = 2^(Delta-1) * a(k-Delta+1),  a(k)=B(k)-2B(k-1)=(4^(k-1)+6k-16)/9.
  It rests on the doubling  n(Delta+1,k+1)=2 n(Delta,k), which (linear independence in k) reduces to
  c_{Delta+1}(D+1)=2 c_Delta(D), where c_Delta(D)=#coprime reduced gens (u,w), v(u)=0, v(w)=Delta,
  max(deg u,deg w)=D. That lemma is proven by three ingredients (verify_doubling_lemma):
    (M) bijection (u,w)->(u,(1+x)w) [gap Delta<->Delta+1],  (R) reflection (u,w)->(u+w,w) [tie<->w-dom,
    so E=W],  (P) per-w doubling [#u of given degree coprime to w doubles in the degree, and gcd(u,w)=1
    forces v(u)=0 since (1+x)|w].

OPEN (the MacWilliams kernel): the saturated per-(Delta, SIZE) split. The saturated obstruction = min ODD
weight of the QC code <(a,b)> = min over ALL s of popcount(a*s)+popcount(b*s) (odd); s-decoupled, cap
stability-checked, split into the s=1 "gcd-formula" layer popcount(a)+popcount(b) (pre-cancellation) and
the saturated layer (after the -2 cascade), bucketed by Delta.
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

def min_odd_weight_capped(a, b, scap):
    best = 0
    for s in range(1, 1 << (scap + 1)):
        w = popcount(gf2_mul(a, s)) + popcount(gf2_mul(b, s))
        if w % 2 == 1 and (best == 0 or w < best):
            best = w
    return best

def d0_pairs(k):
    """(a, b, Delta) for d=0 hard pairs at body k (a,b coprime reduced generators)."""
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
        out.append((a, b, Delta))
    return out


def fmt(counter):
    return "{" + ", ".join(f"{s}:{c}" for s, c in sorted(counter.items())) + "}"


def B(k): return (4 ** k - 12 * k + 8) // 18
def a_seed(m): return (4 ** (m - 1) + 6 * m - 16) // 9   # = B(m) - 2*B(m-1) = n(Delta=1, body m)
def n_bucket(Delta, k): return (1 << (Delta - 1)) * a_seed(k - Delta + 1)   # closed-form Delta-bucket count


def c_split(Delta, D):
    """(E, U, W) for coprime reduced gens (u, w), v(u)=0, v(w)=Delta, max(deg u, deg w)=D,
    split by tie (deg u=deg w=D) / u-dominant (deg u=D>deg w) / w-dominant (deg w=D>deg u)."""
    top = 1 << (D + 1)
    us = [u for u in range(1, top) if valuation(u) == 0 and deg(u) <= D]
    ws = [w for w in range(1, top) if valuation(w) == Delta and deg(w) <= D]
    E = U = W = 0
    for u in us:
        for w in ws:
            if max(deg(u), deg(w)) != D or gf2_gcd(u, w) != 1:
                continue
            if deg(u) == deg(w):
                E += 1
            elif deg(u) > deg(w):
                U += 1
            else:
                W += 1
    return E, U, W


def verify_doubling_lemma(Dmax=6):
    """Derivation of the doubling n(Delta+1,k+1)=2 n(Delta,k). Linear independence in k reduces it to
    c_{Delta+1}(D+1) = 2 c_Delta(D); the three proof ingredients are each checked here:
      (M) bijection (u,w)->(u,(1+x)w):  c_{Delta+1}(D+1) = U_Delta(D+1) + E_Delta(D) + W_Delta(D)
      (R) reflection (u,w)->(u+w,w):    E_Delta(D) = W_Delta(D)
      (P) per-w doubling:               U_Delta(D+1) = 2 (E_Delta(D) + U_Delta(D))
    Combine: c_{Delta+1}(D+1) = U(D+1)+E(D)+W(D) = 2E(D)+2(E(D)+U(D)) = 4E(D)+2U(D) = 2 c_Delta(D)."""
    print("=" * 100)
    print("Doubling-lemma derivation:  c_{Delta+1}(D+1) = 2 c_Delta(D)  via (M) bijection / (R) reflection / (P) per-w")
    print("=" * 100)
    ok = True
    for Delta in (1, 2, 3):
        for D in range(Delta, Dmax):
            E, U, W = c_split(Delta, D)
            cD1 = sum(c_split(Delta + 1, D + 1))         # c_{Delta+1}(D+1)
            Un1 = c_split(Delta, D + 1)[1]               # U_Delta(D+1)
            iM = cD1 == (Un1 + E + W)                    # (M) M-pullback identity
            iR = E == W                                  # (R) reflection
            iP = Un1 == 2 * (E + U)                      # (P) per-w doubling
            lemma = cD1 == 2 * (E + U + W)               # the lemma
            ok = ok and iM and iR and iP and lemma
            if not (iM and iR and iP and lemma):
                print(f"  !! Delta={Delta} D={D}: M={iM} R={iR} P={iP} lemma={lemma}")
    print(f"  ==> ingredients (M),(R),(P) and the lemma all hold, Delta=1..3, D up to {Dmax - 1}: {ok}\n")
    return ok


def verify_bucket_counts(kmax=9):
    """The Delta-bucket COUNT is closed: n(Delta,k) = 2^(Delta-1) * a(k-Delta+1),
    a(m) = B(m)-2B(m-1) = (4^(m-1)+6m-16)/9.  Falls out of the doubling recurrence
    n(Delta+1,k+1)=2 n(Delta,k) + the proven Sum_Delta n = B(k) (telescoping). Verify k=4..kmax."""
    print("=" * 100)
    print("Delta-bucket COUNT closed form:  n(Delta,k) = 2^(Delta-1) * a(k-Delta+1),  a(m)=B(m)-2B(m-1)")
    print("=" * 100)
    assert all(B(m) - 2 * B(m - 1) == a_seed(m) for m in range(2, kmax + 2)), "a(m) closed form broken"
    all_ok = True
    for k in range(4, kmax + 1):
        cnt = Counter(D for (_a, _b, D) in d0_pairs(k))
        row, ok_k = [], True
        for D in sorted(cnt):
            pred = n_bucket(D, k)
            ok = pred == cnt[D]
            ok_k = ok_k and ok
            row.append(f"D{D}:{cnt[D]}" + ("" if ok else f"!=PRED{pred}"))
        tot_ok = sum(cnt.values()) == B(k)
        all_ok = all_ok and ok_k and tot_ok
        print(f"  k={k}: {' '.join(row)}   total={sum(cnt.values())} B(k)={B(k)} "
              f"{'OK' if ok_k and tot_ok else 'BAD'}")
    print(f"  ==> bucket-count closed form {'VERIFIED' if all_ok else 'BROKEN'} k=4..{kmax}\n")
    return all_ok


def main():
    import sys
    SCAP = {4: (6, 9), 5: (6, 9), 6: (8, 11), 7: (8, 11)}
    verify_doubling_lemma(6)
    verify_bucket_counts(9)

    print("=" * 100)
    print("SATURATED per-(Delta, size) distribution of the d=0 obstruction  (s-decoupled, stability-checked)")
    print("  formula = popcount(a)+popcount(b) (s=1, pre-cancellation);  saturated = true code min odd weight")
    print("=" * 100)
    for k in (4, 5, 6, 7):
        pairs = d0_pairs(k)
        lo, hi = SCAP[k]
        # stability check on the global saturated total
        tot_lo, tot_hi = Counter(), Counter()
        for (a, b, Delta) in pairs:
            tot_lo[min_odd_weight_capped(a, b, lo)] += 1
            tot_hi[min_odd_weight_capped(a, b, hi)] += 1
        stable = tot_lo == tot_hi
        print(f"\nk={k}  (saturation {'STABLE' if stable else '!! STILL MOVING !!'} "
              f"scap {lo} vs {hi})   total={fmt(tot_hi)}  B(k)={sum(tot_hi.values())}")
        by_delta_f = defaultdict(Counter)   # gcd-formula (s=1)
        by_delta_s = defaultdict(Counter)   # saturated
        for (a, b, Delta) in pairs:
            by_delta_f[Delta][popcount(a) + popcount(b)] += 1
            by_delta_s[Delta][min_odd_weight_capped(a, b, hi)] += 1
        for Delta in sorted(by_delta_f):
            nf = sum(by_delta_f[Delta].values())
            print(f"   Delta={Delta} (n={nf:3d}):  formula {fmt(by_delta_f[Delta])}"
                  f"  -->  saturated {fmt(by_delta_s[Delta])}")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
