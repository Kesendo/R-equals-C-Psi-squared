#!/usr/bin/env python3
"""Stage 3c empirical anchor: does y_par split the F87 trichotomy at k=3?

Re-enumerates the 294 Klein-homogeneous + Y-par-homogeneous k=3 pairs at N=4
(same enumeration as klein_dissipator_resonance.py) and buckets each pair by
(Klein cell, dephasing letter, y_par, trichotomy class). If y_par is a
trichotomy-independent axis, the y_par=0 and y_par=1 counts within each
(Klein, dephase, class) cell are distributed by combinatorial weight only; if
y_par actually splits the trichotomy, we should see distinct count patterns.

Output: a 4 (Klein) × 3 (letter) table per trichotomy class, with each cell
showing y_par=0 / y_par=1 / total.
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


def enumerate_z2_homogeneous_k3(letters=('I', 'X', 'Y', 'Z')):
    """Unordered Klein-homogeneous + Y-par-homogeneous k=3 pairs; returns
    list of (pair_terms, klein, y_par)."""
    seen = set()
    items = []
    for term1 in product(letters, repeat=3):
        for term2 in product(letters, repeat=3):
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
            items.append(([term1, term2], fw.klein_index(term1), y1))
    return items


def main():
    N = 4
    chain = fw.ChainSystem(N=N)
    items = enumerate_z2_homogeneous_k3()
    print(f"=== Stage 3c: F87 trichotomy × y_par at N={N}, k=3 ===")
    print(f"Enumerated {len(items)} Klein-homogeneous + Y-par-homogeneous pairs.")
    print()

    kleins = [(0, 0), (0, 1), (1, 0), (1, 1)]
    letters = ('Z', 'X', 'Y')
    classes = ('truly', 'soft', 'hard')

    # counts[klein][letter][y_par][class] = int
    counts = {
        klein: {
            letter: {y: {cls: 0 for cls in classes} for y in (0, 1)}
            for letter in letters
        }
        for klein in kleins
    }

    start = time.time()
    for terms, klein, y_par in items:
        for letter in letters:
            cls = fw.classify_pauli_pair(chain, terms, dephase_letter=letter)
            counts[klein][letter][y_par][cls] += 1
    elapsed = time.time() - start
    print(f"Classification complete in {elapsed:.1f}s.")
    print()

    # Per-trichotomy class tables
    for cls in classes:
        print(f"--- trichotomy class: {cls} ---")
        print(f"{'Klein':<8s}  {'Z (y0/y1/tot)':>16s}  {'X (y0/y1/tot)':>16s}  {'Y (y0/y1/tot)':>16s}")
        print("-" * 64)
        for klein in kleins:
            cells = []
            for letter in letters:
                c0 = counts[klein][letter][0][cls]
                c1 = counts[klein][letter][1][cls]
                cells.append(f"{c0:>3d}/{c1:>3d}/{c0+c1:>3d}")
            print(f"{str(klein):<8s}  {cells[0]:>16s}  {cells[1]:>16s}  {cells[2]:>16s}")
        print()

    # Pivoted summary: total pairs by (klein, y_par)
    print("--- enumeration breakdown by (Klein, y_par) ---")
    print(f"{'Klein':<8s}  {'y_par=0':>10s}  {'y_par=1':>10s}  {'total':>10s}")
    print("-" * 44)
    enum_counts = {(k, y): 0 for k in kleins for y in (0, 1)}
    for _, k, y in items:
        enum_counts[(k, y)] += 1
    for klein in kleins:
        c0 = enum_counts[(klein, 0)]
        c1 = enum_counts[(klein, 1)]
        print(f"{str(klein):<8s}  {c0:>10d}  {c1:>10d}  {c0+c1:>10d}")
    print()

    # Splitting verdict per (klein, letter, class): does y_par=0 differ from y_par=1?
    print("--- Splitting verdict: y_par splits each (Klein, letter, class)? ---")
    print(f"{'Klein':<8s}  {'letter':<8s}  {'class':<8s}  {'y0':>4s}  {'y1':>4s}  splits?")
    print("-" * 56)
    any_split = False
    for klein in kleins:
        for letter in letters:
            for cls in classes:
                c0 = counts[klein][letter][0][cls]
                c1 = counts[klein][letter][1][cls]
                if c0 == 0 and c1 == 0:
                    continue
                # Relative test: if the y_par=0/y_par=1 ratio matches the enumeration
                # breakdown (proportional), no split; otherwise split.
                enum_y0 = enum_counts[(klein, 0)]
                enum_y1 = enum_counts[(klein, 1)]
                if enum_y0 == 0 or enum_y1 == 0:
                    verdict = "(degenerate enum)"
                else:
                    # expected if no split: c0/enum_y0 == c1/enum_y1
                    p0 = c0 / enum_y0
                    p1 = c1 / enum_y1
                    diff = abs(p0 - p1)
                    if diff < 1e-9:
                        verdict = "no split (proportional)"
                    else:
                        verdict = f"SPLIT  p0={p0:.3f} p1={p1:.3f}"
                        any_split = True
                print(f"{str(klein):<8s}  {letter:<8s}  {cls:<8s}  {c0:>4d}  {c1:>4d}  {verdict}")
    print()
    if any_split:
        print(">>> y_par SPLITS at least one (Klein, letter, class) cell -> 3c is real <<<")
    else:
        print(">>> y_par is trichotomy-blind across all cells -> 3c is a non-splitting claim <<<")


if __name__ == "__main__":
    main()
