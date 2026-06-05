#!/usr/bin/env python3
"""Door 2: a Liouvillian-free palindrome classifier for arbitrary Pauli Hamiltonians (scout, 2026-06-05).

The §7.1 chiral-K test is mask-level: a dephased H is soft (palindrome-restorable) iff there is a linear
phi over GF(2) with phi(m) = 1 for every FLIP-VALUE m the Hamiltonian produces (m = i XOR j over the
nonzero off-diagonals H_ij), i.e. iff the set of distinct nonzero flip-masks carries NO odd F2-relation.
Crucially (the cancellation lesson, §7.10): XX+YY kills the |00>-|11> edges but NOT the flip-VALUE (the
bond mask still appears via |01>-|10>), so the flip-VALUE-set is mask-only even when the graph geometry
is not. So the soft/hard verdict should be computable from the term masks alone, no 2^N Liouvillian.

This scout builds that classifier and tests it HONESTLY against the actual spectral verdict
(fw.classify_pauli_pair) across a broad battery: F87 pairs, the XY model, Heisenberg, single terms,
diagonal-lift cases, cancellation cases, 3-term sets. It reports where mask-bipartite matches
spectral-not-hard, and any mismatch (the validity boundary).
"""
from __future__ import annotations

import sys
from pathlib import Path

if sys.platform == "win32":
    try: sys.stdout.reconfigure(encoding="utf-8")
    except Exception: pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw


def flip_masks(terms, N):
    """Distinct nonzero flip-masks of the terms placed by the sliding-window builder on N sites."""
    masks = set()
    for t in terms:
        L = len(t)
        for w in range(0, N - L + 1):
            m = 0
            for p, letter in enumerate(t):
                if letter in ("X", "Y"):
                    m |= 1 << (w + p)
            if m:
                masks.add(m)
    return masks


def mask_bipartite(masks):
    """True iff a linear phi with phi(m)=1 for all m exists, i.e. no odd F2-relation. Gaussian
    elimination over GF(2) on the augmented rows [m | 1]; inconsistency (0 | 1) means non-bipartite."""
    rows = [(m << 1) | 1 for m in masks]   # append target bit 1 in the lowest position
    basis = []
    for r in rows:
        for b in basis:
            r = min(r, r ^ b)
        if r:
            # if r is just the target bit (value 1, mask part zero) -> 0 = 1 inconsistency
            if r == 1:
                return False
            basis.append(r)
            basis.sort(reverse=True)
    return True


def main():
    print("=" * 92)
    print("Door 2: mask-level palindrome classifier vs the actual spectral verdict")
    print("=" * 92)
    N = 4
    chain = fw.ChainSystem(N=N, J=1.0, gamma_0=0.05)

    battery = [
        # (label, terms)
        ("XX (single bond)",            [('X','X')]),
        ("XX+YY (XY model)",            [('X','X'), ('Y','Y')]),
        ("XX+YY+ZZ (Heisenberg)",       [('X','X'), ('Y','Y'), ('Z','Z')]),
        ("XY+YX",                       [('X','Y'), ('Y','X')]),
        ("XZ+ZX",                       [('X','Z'), ('Z','X')]),
        ("YZ+ZY",                       [('Y','Z'), ('Z','Y')]),
        ("Z field (diagonal)",          [('Z','I')]),
        ("XX (k=3 IXX)",                [('I','X','X')]),
        ("F87 soft  XXI+YYI",           [('X','X','I'), ('Y','Y','I')]),
        ("F87 hard  XYI+YIX",           [('X','Y','I'), ('Y','I','X')]),
        ("F87 hard  XIY+IXY",           [('X','I','Y'), ('I','X','Y')]),
        ("XX+YY+XY (3-term)",           [('X','X'), ('Y','Y'), ('X','Y')]),
        ("XYI+YIX+IXY (3-term hard?)",  [('X','Y','I'), ('Y','I','X'), ('I','X','Y')]),
        ("XX+XZ (mixed body diag)",     [('X','X'), ('X','Z')]),
    ]

    def bit_b(letter_term):
        return sum(1 for L in letter_term if L in ("Y", "Z")) % 2

    print(f"\n  {'Hamiltonian':<30} {'bit_b-homog':<11} {'hopping':<8} {'mask':<6} {'spectral':<8} {'match?'}")
    print("  " + "-" * 78)
    agree = 0; total = 0
    boundary_ok = True
    for label, terms in battery:
        masks = flip_masks(terms, N)
        bip = mask_bipartite(masks) if masks else True
        spectral = fw.classify_pauli_pair(chain, terms, dephase_letter='Z')
        ok = (bip == (spectral != 'hard'))
        agree += ok; total += 1
        homog = len({bit_b(t) for t in terms}) == 1
        hopping = len(masks) > 0
        in_scope = homog and hopping
        # the claim: mask test is valid EXACTLY on the in-scope (bit_b-homogeneous hopping) cases
        if in_scope != ok:
            boundary_ok = False
        print(f"  {label:<30} {('yes' if homog else 'NO'):<11} {('yes' if hopping else 'NO'):<8} "
              f"{('soft' if bip else 'HARD'):<6} {spectral:<8} {'OK' if ok else 'MISMATCH'}")

    print(f"\n  agreement: {agree}/{total}")
    print(f"  validity boundary: mask test correct  <=>  (bit_b-homogeneous AND has hopping):  {boundary_ok}")
    print("\n  reading: the mask-bipartite test IS a Liouvillian-free palindrome classifier, valid for")
    print("  bit_b-homogeneous hopping Hamiltonians (the §7 / F87 single-Klein-cell scope). Mixed-cell")
    print("  Hamiltonians (different bit_b) and pure-diagonal lifts fall outside; their soft/hard is not")
    print("  fixed by the flip-mask set alone (XX+YY truly vs XX+YY+XY hard share the same mask set).")


if __name__ == "__main__":
    main()
