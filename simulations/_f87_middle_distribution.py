#!/usr/bin/env python3
"""F115 open thread: the per-size DISTRIBUTION of the hardness obstruction (2026-06-08).

Closed so far: max size = min(2W-1, 2k-3-2d); total count A203241 = (4^(k-1)-3*2^(k-1)+2)/3;
d-layer count 2^(d-1)*B(k-d), B(k)=(4^k-12k+8)/18; triangle (size-3) sub-count. OPEN: the full per-size
distribution #{hard mask-pairs : obstruction size = s} for the MIDDLE s, "window-dependent" (the obstruction
is the min ODD weight of the quasi-cyclic code <(a,b)>, a=p2/g, b=p1/g, and a sparse multiple can cancel
below the gcd-formula popcount).

This scout does three things:
  (1) print the empirical per-size distribution across a (k,N) grid -> see the window-dependence;
  (2) TEST the d-reduction hypothesis: layer-d distribution =?= 2^(d-1) * (d=0 distribution at body k-d),
      same window count W. If it holds, the open problem collapses to the d=0 layer;
  (3) print the d=0 distribution as W grows -> isolate the core window-dependence to attack with MacWilliams.
Self-validates: total == A203241, per-layer max == 2k-3-2d.
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
            for m in combo:
                x ^= m
            if x == 0:
                return size
    return 0

def B(k): return (4 ** k - 12 * k + 8) // 18
def A203241(k): return (4 ** (k - 1) - 3 * 2 ** (k - 1) + 2) // 3


def distribution(k, N):
    """over hard mask-pairs at (k,N): Counter[(d,size)], Counter[size], W."""
    W = N - k + 1
    masks = even_masks(k)
    by_ds, by_s = Counter(), Counter()
    for p1, p2 in combinations(masks, 2):
        if valuation(p1) == valuation(p2):
            continue  # soft
        d = g_rest_degree(p1, p2)
        S = {p1 << w for w in range(W)} | {p2 << w for w in range(W)}
        s = min_odd_dependence(S)
        by_ds[(d, s)] += 1
        by_s[s] += 1
    return by_ds, by_s, W


def main():
    print("=" * 88)
    print("(1) per-size obstruction distribution by (k,N) [size:count], total vs A203241")
    print("=" * 88)
    for k in range(3, 7):
        for N in range(k + 1, k + 5):
            by_ds, by_s, W = distribution(k, N)
            total = sum(by_s.values())
            sizes = ", ".join(f"{s}:{by_s[s]}" for s in sorted(by_s))
            tag = "" if total == A203241(k) else f"  (A203241={A203241(k)} !!)"
            print(f"  k={k} N={N:2d} W={W}: total={total}{tag}   [{sizes}]")
        print()

    print("=" * 88)
    print("(2) d-reduction test: dist(layer d, body k) =?= 2^(d-1) * dist(d=0, body k-d), same W (N=2k)")
    print("=" * 88)
    all_ok = True
    for k in range(4, 8):
        N = 2 * k
        by_ds, _, W = distribution(k, N)
        ds = sorted({d for (d, s) in by_ds})
        for d in ds:
            if d == 0:
                continue
            kp = k - d
            if kp < 3:
                continue
            Np = W + kp - 1  # same window count W
            base_ds, _, Wp = distribution(kp, Np)
            assert Wp == W
            layer = {s: c for (dd, s), c in by_ds.items() if dd == d}
            base0 = {s: c for (dd, s), c in base_ds.items() if dd == 0}
            scaled = {s: (1 << (d - 1)) * c for s, c in base0.items()}
            match = layer == scaled
            all_ok = all_ok and match
            print(f"  k={k} d={d} (vs body {kp}, W={W}): "
                  f"layer={dict(sorted(layer.items()))}  2^(d-1)*base0={dict(sorted(scaled.items()))}  "
                  f"{'MATCH' if match else 'MISMATCH'}")
        print()
    print(f"  ==> d-reduction hypothesis: {'HOLDS on this grid' if all_ok else 'FAILS (see mismatches)'}\n")

    print("=" * 88)
    print("(3) the d=0 layer distribution as W grows [size:count] (the core window-dependence)")
    print("=" * 88)
    for k in range(4, 7):
        print(f"  body k={k}, d=0:")
        for N in range(k + 1, 2 * k + 2):
            by_ds, _, W = distribution(k, N)
            d0 = {s: c for (dd, s), c in by_ds.items() if dd == 0}
            cap = 2 * k - 3
            wcap = 2 * W - 1
            sizes = ", ".join(f"{s}:{d0[s]}" for s in sorted(d0))
            print(f"    N={N:2d} W={W} (min(2W-1,2k-3)=min({wcap},{cap})={min(wcap,cap)}): [{sizes}]")
        print()


if __name__ == "__main__":
    main()
