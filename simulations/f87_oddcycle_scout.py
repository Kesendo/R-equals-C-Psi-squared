#!/usr/bin/env python3
"""F87 windowed converse (task B): odd-cycle ⟺ hard ⟺ block-asymmetry (scout, 2026-06-04).

The one open edge of the 42:8 question is the windowed (k<N) converse: non-bipartite ⟹ hard.
Three readings of the SAME diagonal-cell pairs live in the repo:
  GF(2):  S = edge-XOR mask set has an odd 𝔽₂-relation (odd-size subset XOR-ing to 0)
          ⟺ non-bipartite  (phi_exists_gf2, f87_flip_generators.py)
  class:  the pair is F87-hard  (classify_pauli_pair)
  block:  some degenerate first-order D̂-block is asymmetric about −N  (f87_block_localize.py)

"bipartite ⟺ no odd relation" is rigorous graph theory; "bipartite ⟹ soft" is derived (§7.1).
The lone open converse is "odd relation ⟹ hard", equivalently "odd relation ⟹ some block
asymmetric". This scout (a) confirms the three readings coincide bit-exact for the windowed
(N=4, k=3) Z-deph diagonal cell, and (b) extracts the minimal odd cycle for each hard pair,
the explicit 𝔽₂ object a direct proof of the converse must consume.
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
from f87_flip_generators import edge_generators, phi_exists_gf2
from f87_block_localize import blocks, asym_about

DIAG = {"I", "Z"}  # Z-deph diagonal letters


def is_mixed(t):
    return any(L not in DIAG for L in t)


def min_odd_cycle(S):
    """Smallest odd-size subset of S that XORs to 0 (the minimal odd 𝔽₂-relation); None if bipartite."""
    S = sorted(S)
    for size in (3, 5, 7):
        for combo in combinations(S, size):
            x = 0
            for m in combo:
                x ^= m
            if x == 0:
                return combo
    return None


def block_is_asymmetric(pair, N=4, tol=1e-6):
    """True iff some degenerate first-order D̂-block is asymmetric about −N (the 'hard' block signature)."""
    blk, _ = blocks(pair, N)
    worst = 0.0
    for _omega, (modes, M) in blk.items():
        s = np.linalg.eigvals(M).real
        worst = max(worst, asym_about(s, -N))
    return worst > tol, worst


def main():
    N, k = 4, 3
    chain = fw.ChainSystem(N=N)
    terms = [t for t in product("IXYZ", repeat=k)
             if not all(L == "I" for L in t) and fw.klein_index(t) == (0, 1)]
    mixed = [t for t in terms if is_mixed(t)]

    print("=" * 92)
    print("F87 windowed converse scout: odd-cycle ⟺ hard ⟺ block-asymmetry  (N=4, k=3, Z-deph (0,1))")
    print("=" * 92)
    print()

    rows = []
    for t1, t2 in combinations_with_replacement(mixed, 2):
        if (sum(c == "Y" for c in t1) % 2) != (sum(c == "Y" for c in t2) % 2):
            continue  # y_par-homogeneous only
        H = _build_kbody_chain(N, [tuple(t1) + (1.0,), tuple(t2) + (1.0,)])
        S = edge_generators(H)
        odd = not phi_exists_gf2(S, N)                       # GF(2): has odd relation
        cls = fw.classify_pauli_pair(chain, [t1, t2], dephase_letter="Z")
        hard = (cls == "hard")
        asym, worst = block_is_asymmetric([t1, t2], N)        # block reading
        rows.append((t1, t2, S, odd, hard, asym, worst))

    # three-way agreement
    n = len(rows)
    agree_oh = sum(1 for r in rows if r[3] == r[4])            # odd == hard
    agree_hb = sum(1 for r in rows if r[4] == r[5])            # hard == block-asym
    n_hard = sum(1 for r in rows if r[4])
    n_odd = sum(1 for r in rows if r[3])
    n_asym = sum(1 for r in rows if r[5])
    print(f"  pairs (y_par-homogeneous Mixed+Mixed): {n}")
    print(f"  odd-relation: {n_odd}    hard: {n_hard}    block-asymmetric: {n_asym}")
    print(f"  odd ⟺ hard agreement:            {agree_oh}/{n}  {'ALL' if agree_oh == n else 'MISMATCH'}")
    print(f"  hard ⟺ block-asym agreement:     {agree_hb}/{n}  {'ALL' if agree_hb == n else 'MISMATCH'}")
    assert agree_oh == n and agree_hb == n, "three-way equivalence broken"
    print(f"  => odd-relation ⟺ hard ⟺ block-asymmetric, bit-exact ({n} pairs).")
    print()

    # the minimal odd cycle for the hard pairs: the explicit 𝔽₂ object
    print("  Minimal odd 𝔽₂-cycle for the hard pairs (mask = X/Y bit positions in the N-bit string):")
    print(f"  {'t1':>10} {'t2':>10} {'|S|':>4} {'min odd cycle (masks, binary)':>42}")
    cyc_sizes = {}
    for (t1, t2, S, odd, hard, asym, worst) in rows:
        if not hard:
            continue
        cyc = min_odd_cycle(S)
        size = len(cyc) if cyc else 0
        cyc_sizes[size] = cyc_sizes.get(size, 0) + 1
        masks = "  ".join(format(m, f"0{N}b") for m in cyc) if cyc else "none"
        print(f"  {''.join(t1):>10} {''.join(t2):>10} {len(S):>4}   {masks}")
    print()
    print(f"  minimal-odd-cycle size distribution over the {n_hard} hard pairs: {dict(sorted(cyc_sizes.items()))}")
    print()
    print("  reading: every hard pair carries a minimal odd cycle (size 3 = a triangle in the Cayley")
    print("  graph of S). The converse to derive: this odd cycle FORCES some degenerate D̂-block")
    print("  asymmetric about −N (the block entries are built from the same XOR masks via the")
    print("  Absorption Theorem D̂|i⟩⟨j| = −2·popcount(i⊕j)). That 𝔽₂ ⟹ spectral implication is the gap.")


if __name__ == "__main__":
    main()
