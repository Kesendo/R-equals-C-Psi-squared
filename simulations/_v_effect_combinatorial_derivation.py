#!/usr/bin/env python3
"""V-Effect 14/19/3 split — full combinatorial derivation at N=3.

For 36 unordered cross-pairs of bilinears T_k ∈ {X,Y,Z}^2 (the V-Effect
36-enum), the truly/soft/hard classification can be derived purely from
algebraic features of the pair (T1, T2) — no L computation needed.

Verified: all 36 cases match the empirical framework classification.

  TRULY (3):      both T1, T2 ∈ {XX, YY, ZZ}
                   (cross-pairs of both-parity-even bilinears)

  1-BPE (18):     exactly one of T1, T2 in {XX, YY, ZZ}
    HARD (4):     BPE = (a,a) with a ∈ {X, Y};
                   non-BPE shares 'a' at one site AND has bit_a-partner
                   ('Y' if a='X', 'X' if a='Y') at the other site
    SOFT (14):    everything else in 1-BPE (includes all 6 ZZ-pairs)

  0-BPE (15):     neither T1 nor T2 in {XX, YY, ZZ}
    SOFT (5):     bond-flipped (T2 = reverse(T1)) OR
                   Z-aligned (Z at same site in both T1 and T2)
    HARD (10):    everything else

  Total = 3 + 19 + 14 = 36 ✓

Mechanism: the V-Effect at N=3 emerges from the compatibility
requirements between the two bond bilinears. Compatibility means:
both-parity-even (truly), or bond-symmetric/Z-aligned (soft), or
bit_a-conflict-avoidance (1-BPE soft). When compatibility breaks
in any of these forms, the case becomes hard.
"""
import math
import sys
from itertools import combinations
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


BPE = {('X', 'X'), ('Y', 'Y'), ('Z', 'Z')}
PAULI_XYZ = ['X', 'Y', 'Z']
BILINEARS = [(a, b) for a in PAULI_XYZ for b in PAULI_XYZ]


def classify_combinatorial(T1, T2):
    """Combinatorial classification — pure algebra, no L computation."""
    in_BPE_1 = T1 in BPE
    in_BPE_2 = T2 in BPE

    # TRULY
    if in_BPE_1 and in_BPE_2:
        return 'truly'

    # 1-BPE
    if in_BPE_1 ^ in_BPE_2:
        bpe = T1 if in_BPE_1 else T2
        non_bpe = T2 if in_BPE_1 else T1
        a = bpe[0]  # bpe = (a, a)
        if a == 'Z':
            return 'soft'  # ZZ has no bit_a-partner
        partner = 'Y' if a == 'X' else 'X'
        if (non_bpe[0] == a and non_bpe[1] == partner) or \
           (non_bpe[1] == a and non_bpe[0] == partner):
            return 'hard'
        return 'soft'

    # 0-BPE
    if T1[::-1] == T2:
        return 'soft'  # bond-flipped
    if (T1[0] == 'Z' and T2[0] == 'Z') or \
       (T1[1] == 'Z' and T2[1] == 'Z'):
        return 'soft'  # Z-aligned
    return 'hard'


def classify_empirical(H, N, gamma):
    """Empirical classification via framework primitives."""
    L = fw.lindbladian_z_dephasing(H, [gamma] * N)
    M = fw.palindrome_residual(L, N * gamma, N)
    if float(np.linalg.norm(M)) < 1e-10:
        return 'truly'
    evals = np.linalg.eigvals(L)
    used = np.zeros(len(evals), dtype=bool)
    max_err = 0.0
    sigma_g = N * gamma
    for i in range(len(evals)):
        if used[i]:
            continue
        target = -evals[i] - 2 * sigma_g
        dists = np.abs(evals - target)
        for j in range(len(evals)):
            if used[j]:
                dists[j] = np.inf
        bj = int(np.argmin(dists))
        used[i] = True
        if bj != i:
            used[bj] = True
        max_err = max(max_err, float(dists[bj]))
    return 'soft' if max_err < 1e-6 else 'hard'


def main():
    N = 3
    GAMMA, J = 0.1, 1.0
    bonds = [(i, i + 1) for i in range(N - 1)]

    print(f"V-Effect 14/19/3 combinatorial derivation at N={N}")
    print(f"  γ_dephasing = {GAMMA}")
    print(f"  Total pairs: 36 (cross-pairs of bilinears in {{X,Y,Z}}²)")
    print()

    # Verify all 36 cases
    n_correct = 0
    n_total = 0
    mismatches = []
    counts_combinatorial = {'truly': 0, 'soft': 0, 'hard': 0}
    for t1, t2 in combinations(BILINEARS, 2):
        pair = tuple(sorted([t1, t2]))
        terms = [(pair[0][0], pair[0][1], J), (pair[1][0], pair[1][1], J)]
        H = fw._build_bilinear(N, bonds, terms)
        cat_emp = classify_empirical(H, N, GAMMA)
        cat_comb = classify_combinatorial(pair[0], pair[1])
        n_total += 1
        counts_combinatorial[cat_comb] = counts_combinatorial.get(cat_comb, 0) + 1
        if cat_emp == cat_comb:
            n_correct += 1
        else:
            label = f"{pair[0][0]}{pair[0][1]}+{pair[1][0]}{pair[1][1]}"
            mismatches.append((label, cat_emp, cat_comb))

    print(f"Verification: {n_correct}/{n_total} cases match")
    print(f"Combinatorial counts: {counts_combinatorial}")
    print()

    if mismatches:
        print("Mismatches:")
        for label, emp, comb in mismatches:
            print(f"  {label}: empirical={emp}, combinatorial={comb}")
    else:
        print("✓ Full algebraic derivation of V-Effect 14/19/3 split achieved.")
        print()
        print("The combinatorial rules (no Lindbladian computation needed):")
        print()
        print("  Rule TRULY:  both T1 AND T2 ∈ {XX, YY, ZZ}")
        print("  Rule 1-BPE:  exactly one of T1, T2 ∈ {XX, YY, ZZ}")
        print("    HARD: BPE = (a,a), a ∈ {X,Y}, non-BPE has 'a' at one site")
        print("           AND bit_a-partner (Y if a=X, X if a=Y) at other")
        print("    SOFT: otherwise (includes all ZZ-1-BPE pairs)")
        print("  Rule 0-BPE:  neither T1 nor T2 in {XX, YY, ZZ}")
        print("    SOFT: bond-flipped (T2 = reverse(T1)) OR Z-aligned")
        print("           (Z at same site in both T1 and T2)")
        print("    HARD: otherwise")
        print()
        print("  Sum: 3 truly + 19 soft + 14 hard = 36 ✓")


if __name__ == "__main__":
    main()
