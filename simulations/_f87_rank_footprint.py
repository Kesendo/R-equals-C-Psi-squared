#!/usr/bin/env python3
"""The rank's spectral footprint, corrected (2026-06-05).

CORRECTION of an earlier overclaim (dim ker L = 2^(N-r)): that held only for the HARD pairs first
tested. The right, universal statement is:

  span-rank r  ->  #hopping-components = 2^(N-r) = dim(diagonal commutant) = the # DIAGONAL conserved
                   sectors = the gain-channel +N Perron multiplicity n_plusN  (letter-independent).

These diagonal conserved sectors always sit inside ker L. For HARD pairs they are all of ker L, so
dim ker L = 2^(N-r). For SOFT (and partial-support) pairs there are ADDITIONAL steady states (the chiral
symmetry / free sites), so dim ker L > 2^(N-r). So dim ker L is NOT the clean rank footprint;
#components / the diagonal conserved sectors / n_plusN is. This script shows the distinction for soft
and hard pairs side by side.
"""
from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try: sys.stdout.reconfigure(encoding="utf-8")
    except Exception: pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw
from framework.pauli import _build_kbody_chain
from framework.lindblad import lindbladian_pauli_dephasing

def gf2_rank(vectors):
    basis = []
    for v in vectors:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v); basis.sort(reverse=True)
    return len(basis)

def windowed(mask, k, N): return [mask << w for w in range(N - k + 1)]
def mask_of(letters): return sum(1 << i for i, L in enumerate(letters) if L in ("X", "Y"))

def n_components(masks, N):
    parent = list(range(1 << N))
    def find(a):
        while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for i in range(1 << N):
        for m in masks:
            ri, rj = find(i), find(i ^ m)
            if ri != rj: parent[ri] = rj
    return len({find(i) for i in range(1 << N)})

def diagonal_commutant_dim(H, N):
    """# of diagonal operators D (D = diag(d)) with [H,D]=0 = # values d may take = # components."""
    # d_i = d_j whenever H_ij != 0; count connected components of the |H|>tol graph.
    parent = list(range(1 << N))
    def find(a):
        while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
        return a
    A = np.abs(H) > 1e-9
    for i in range(1 << N):
        for j in range(i + 1, 1 << N):
            if A[i, j]:
                ri, rj = find(i), find(j)
                if ri != rj: parent[ri] = rj
    return len({find(i) for i in range(1 << N)})

def dim_ker_L(letter_terms, N, gamma=0.05):
    H = _build_kbody_chain(N, [tuple(t) + (1.0,) for t in letter_terms])
    L = lindbladian_pauli_dephasing(H, [gamma] * N, dephase_letter='Z')
    return int(np.sum(np.abs(np.linalg.eigvals(L)) < 1e-8)), H


def main():
    print("=" * 96)
    print("The rank footprint, corrected: #components = diagonal sectors = 2^(N-r); dim ker L only for hard")
    print("=" * 96)
    N = 4
    pairs = [
        ("hard  XYI+YIX", [('X','Y','I'), ('Y','I','X')]),
        ("hard  XIY+IXY", [('X','I','Y'), ('I','X','Y')]),
        ("soft  XXI+YYI", [('X','X','I'), ('Y','Y','I')]),
        ("soft  XYZ+YXZ", [('X','Y','Z'), ('Y','X','Z')]),
        ("soft  XXI+IXX", [('X','X','I'), ('I','X','X')]),
    ]
    print(f"\n  {'pair':<16} {'r':<3} {'2^(N-r)':<8} {'#comp':<6} {'diag.commutant':<15} {'dim ker L':<10} hard?")
    print("  " + "-" * 78)
    for label, terms in pairs:
        k = len(terms[0])
        flips = []
        for t in terms: flips += windowed(mask_of(t), k, N)
        r = gf2_rank(flips)
        comps = n_components(set(flips), N)
        kerd, H = dim_ker_L(terms, N)
        dc = diagonal_commutant_dim(H, N)
        is_hard = "hard" in label
        foot = 1 << (N - r)
        note = "ker L = footprint" if kerd == foot else f"ker L EXCEEDS footprint by {kerd - foot}"
        print(f"  {label:<16} {r:<3} {foot:<8} {comps:<6} {dc:<15} {kerd:<10} "
              f"{'Y' if is_hard else 'n'}  ({note})")

    print("\n  reading: #components = diagonal-commutant dim = 2^(N-r) is the UNIVERSAL rank footprint")
    print("  (= the +N Perron multiplicity, the # diagonal conserved sectors, letter-independent).")
    print("  dim ker L equals it for HARD pairs but EXCEEDS it for soft / partial-support pairs (extra")
    print("  steady states). So the earlier 'dim ker L = 2^(N-r)' was a hard-pair coincidence; corrected.")


if __name__ == "__main__":
    main()
