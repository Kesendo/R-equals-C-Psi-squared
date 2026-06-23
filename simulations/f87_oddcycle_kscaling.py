#!/usr/bin/env python3
"""How does the windowed odd 𝔽₂-cycle (the F87 hardness obstruction) scale with k? (2026-06-04)

The §7.5/§7.6 derivation proved soft ⟺ bipartite for any k, so bipartiteness of H's hopping graph
IS the ground-truth class (no Liouvillian needed). The ONE k-specific piece is the SHAPE of the odd
𝔽₂-relation in the edge-mask set S: at k=3 it is the K3 triangle (3 popcount-2 masks on 3 consecutive
sites). This probe characterises that shape for k=3,4,5 windowed (k<N), purely combinatorially
(2^N bipartite check + GF(2) mask relations, no eigendecomposition), to see whether it generalises.
"""
from __future__ import annotations

import sys
from itertools import product, combinations, combinations_with_replacement
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw
from framework.pauli import _build_kbody_chain
from f87_flip_generators import edge_generators, is_bipartite

DIAG = {"I", "Z"}


def min_odd_cycle(S):
    S = sorted(S)
    for size in (3, 5, 7, 9):
        for combo in combinations(S, size):
            x = 0
            for m in combo:
                x ^= m
            if x == 0:
                return combo
    return None


def popcount(x):
    return bin(x).count("1")


def site_span(masks):
    """The set of chain sites touched by the cycle's masks, and whether they are consecutive."""
    bits = set()
    for m in masks:
        i = 0
        while m:
            if m & 1:
                bits.add(i)
            m >>= 1
            i += 1
    sites = sorted(bits)
    consecutive = (len(sites) == sites[-1] - sites[0] + 1)
    return len(sites), consecutive


def main():
    print("=" * 90)
    print("F87 windowed odd-cycle (hardness obstruction) vs k  (combinatorial, no Liouvillian)")
    print("=" * 90)
    print()
    for (k, N) in [(3, 4), (4, 5), (5, 6)]:
        terms = [t for t in product("IXYZ", repeat=k)
                 if not all(L == "I" for L in t) and fw.klein_index(t) == (0, 1)
                 and any(L not in DIAG for L in t)]
        size_dist = {}
        pc_seen = set()
        span_seen = set()
        n_hard = 0
        n_pairs = 0
        for t1, t2 in combinations_with_replacement(terms, 2):
            if (sum(c == "Y" for c in t1) % 2) != (sum(c == "Y" for c in t2) % 2):
                continue
            n_pairs += 1
            H = _build_kbody_chain(N, [tuple(t1) + (1.0,), tuple(t2) + (1.0,)])
            if is_bipartite(H):
                continue  # soft
            n_hard += 1
            S = edge_generators(H)
            cyc = min_odd_cycle(S)
            if cyc is None:
                size_dist["none?"] = size_dist.get("none?", 0) + 1
                continue
            size_dist[len(cyc)] = size_dist.get(len(cyc), 0) + 1
            pc_seen.add(tuple(sorted(popcount(m) for m in cyc)))
            nsites, consec = site_span(cyc)
            span_seen.add((nsites, consec))
        print(f"  k={k}, N={N} (windowed, k<N):  {n_pairs} y_par-homog Mixed pairs, {n_hard} hard")
        print(f"     minimal odd-cycle SIZE distribution: {dict(sorted((str(x), c) for x, c in size_dist.items()))}")
        print(f"     popcounts of the cycle masks (sorted tuples seen): {sorted(pc_seen)}")
        print(f"     (#sites touched, consecutive?) seen: {sorted(span_seen)}")
        print()
    print("  reading: if the size stays 3 and masks stay popcount-2 on consecutive sites, the K3")
    print("  triangle IS the general windowed obstruction (k-independent shape). If size/popcount/")
    print("  span grow with k, the obstruction is a k-dependent family and the general form is what")
    print("  a full general-k proof must name (replacing §7.2's triangle).")


if __name__ == "__main__":
    main()
