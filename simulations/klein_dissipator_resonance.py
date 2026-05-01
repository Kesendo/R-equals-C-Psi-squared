#!/usr/bin/env python3
"""Verify the dissipator-resonance law for F77 hardness.

F77-hardness lives in the Klein cell that matches the dephasing letter:
  Z-dephasing → Klein (0, 1)  [Z's Klein index]
  X-dephasing → Klein (1, 0)  [X's Klein index]
  Y-dephasing → Klein (1, 1)  [Y's Klein index]

Method: full enumeration of Klein-homogeneous + Y-par-homogeneous k=3 pairs
at N=4, classified under each Pauli-letter dephasing channel via
`fw.classify_pauli_pair(chain, terms, dephase_letter=...)`. See
hypotheses/THE_POLARITY_LAYER.md for the conceptual reading.
"""
from __future__ import annotations

import sys
import time
from itertools import product
from pathlib import Path

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw


# Klein-cell representatives (k=3 Klein-homogeneous pairs, Y-par 0).
TEST_PAIRS = {
    'Klein (0,0)': [('X', 'X', 'I'), ('Y', 'Y', 'I')],
    'Klein (0,1)': [('I', 'I', 'Z'), ('Z', 'Z', 'Z')],
    'Klein (1,0)': [('I', 'I', 'X'), ('X', 'X', 'X')],
    'Klein (1,1)': [('I', 'I', 'Y'), ('Y', 'Y', 'Y')],
}


def enumerate_z2_homogeneous_k3(letters_alphabet=('I', 'X', 'Y', 'Z')):
    """Unordered Klein-homogeneous + Y-par-homogeneous k=3 pairs."""
    seen = set()
    pairs = []
    for term1 in product(letters_alphabet, repeat=3):
        for term2 in product(letters_alphabet, repeat=3):
            if all(L == 'I' for L in term1) or all(L == 'I' for L in term2):
                continue
            if fw.klein_index(term1) != fw.klein_index(term2):
                continue
            y1 = sum(1 for L in term1 if L == 'Y') % 2
            y2 = sum(1 for L in term2 if L == 'Y') % 2
            if y1 != y2:
                continue
            sorted_terms = tuple(sorted([term1, term2]))
            if sorted_terms in seen:
                continue
            seen.add(sorted_terms)
            pairs.append([term1, term2])
    return pairs


def main():
    N = 4
    chain = fw.ChainSystem(N=N)
    print(f"=== Dissipator-Resonance: F77-Hardness vs Klein cell ===")
    print(f"N = {N}, k = 3 (Z₂³-homogeneous pairs), γ = {chain.gamma_0} per site, J = 1.0")
    print()

    # 1) Targeted 4×3 sanity check
    print("--- Part 1: Targeted 4×3 (one representative per Klein cell) ---")
    print(f"{'Hamiltonian':<24s} {'Z-deph':>10s} {'X-deph':>10s} {'Y-deph':>10s}")
    print("-" * 60)
    for klein_label, terms in TEST_PAIRS.items():
        h_label = '+'.join([''.join(t) for t in terms])
        h_short = f"{klein_label} {h_label}"
        cells = [
            fw.classify_pauli_pair(chain, terms, dephase_letter=letter)
            for letter in ('Z', 'X', 'Y')
        ]
        print(f"{h_short:<24s} {cells[0]:>10s} {cells[1]:>10s} {cells[2]:>10s}")
    print()

    # 2) Full Z₂³-homogeneous sweep × 3 dephasing letters
    print("--- Part 2: Full sweep over all Z₂³-homogeneous pairs ---")
    pairs = enumerate_z2_homogeneous_k3()
    print(f"Enumerated {len(pairs)} Klein-homogeneous + Y-par-homogeneous k=3 pairs.")
    print()

    counts = {
        (klein, letter): {'truly': 0, 'soft': 0, 'hard': 0}
        for klein in [(0, 0), (0, 1), (1, 0), (1, 1)]
        for letter in ('Z', 'X', 'Y')
    }
    start = time.time()
    for terms in pairs:
        klein = fw.klein_index(terms[0])
        for letter in ('Z', 'X', 'Y'):
            cls = fw.classify_pauli_pair(chain, terms, dephase_letter=letter)
            counts[(klein, letter)][cls] += 1
    elapsed = time.time() - start

    print(f"Hardness counts by (Klein cell, dephasing letter)  [elapsed {elapsed:.1f}s]:")
    print()
    print(f"{'Klein cell':<14s} {'Z-deph':>14s} {'X-deph':>14s} {'Y-deph':>14s}  (#hard / #total)")
    print("-" * 75)
    for klein in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        cells = []
        for letter in ('Z', 'X', 'Y'):
            c = counts[(klein, letter)]
            total = c['truly'] + c['soft'] + c['hard']
            cells.append(f"{c['hard']:>3d} / {total:>3d}")
        print(f"{str(klein):<14s} {cells[0]:>14s} {cells[1]:>14s} {cells[2]:>14s}")

    print()
    print("Expected diagonal (conjecture):")
    print("  Z-deph: hardness only in Klein (0, 1)")
    print("  X-deph: hardness only in Klein (1, 0)")
    print("  Y-deph: hardness only in Klein (1, 1)")
    print("  Klein (0, 0): never hard (Mother sector)")


if __name__ == "__main__":
    main()
