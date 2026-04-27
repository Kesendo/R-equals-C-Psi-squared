#!/usr/bin/env python3
"""V-Effect triple-bilinear-per-bond test at chain N=3.

H per bond = T1 + T2 + T3 (three bilinears) instead of pair T1 + T2.

Hypothesis: is the triple's classification derivable from the three
pairwise compatibilities (T1,T2), (T1,T3), (T2,T3)?

  Naive composition rule:
    All three pairs truly  → triple truly
    At least one pair hard → triple hard
    Otherwise (some soft)  → triple soft

If this rule holds for all triples, no new triple-level compatibility
layer. If it fails for some triples, we have a new structural layer
that emerges only when 3+ bilinears co-act per bond.
"""
import sys
from itertools import combinations_with_replacement, combinations, product
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


GAMMA, J = 0.1, 1.0
BPE = {('I', 'I'), ('X', 'X'), ('Y', 'Y'), ('Z', 'Z')}


def n_I(T):
    return T.count('I')


def active_letters(T):
    return [c for c in T if c != 'I']


def classify_combinatorial_pair(T1, T2):
    """V-Effect classifier for pairs (commit 079c7ce)."""
    n1, n2 = n_I(T1), n_I(T2)
    if n1 == 1 and n2 == 1:
        letters = active_letters(T1) + active_letters(T2)
        unique_letters = set(letters)
        if 'Z' in unique_letters:
            return 'hard'
        if unique_letters == {'X'}:
            return 'truly'
        return 'soft'
    if (n1 == 1) != (n2 == 1):
        single = T1 if n1 == 1 else T2
        double = T2 if n1 == 1 else T1
        single_letter = active_letters(single)[0]
        if single_letter == 'Z':
            return 'hard'
        if single_letter == 'X':
            if double in BPE:
                return 'truly'
            if double in {('Y', 'Z'), ('Z', 'Y')}:
                return 'soft'
            return 'hard'
        if single_letter == 'Y':
            if double in BPE or double in {('X', 'Z'), ('Z', 'X')}:
                return 'soft'
            return 'hard'
    if n1 == 0 and n2 == 0:
        in_BPE_1 = T1 in BPE
        in_BPE_2 = T2 in BPE
        if in_BPE_1 and in_BPE_2:
            return 'truly'
        if T1 == T2:
            if T1 in BPE:
                return 'truly'
            return 'soft'
        if in_BPE_1 ^ in_BPE_2:
            bpe = T1 if in_BPE_1 else T2
            non_bpe = T2 if in_BPE_1 else T1
            a = bpe[0]
            if a == 'Z':
                return 'soft'
            partner = 'Y' if a == 'X' else 'X'
            if (non_bpe[0] == a and non_bpe[1] == partner) or \
               (non_bpe[1] == a and non_bpe[0] == partner):
                return 'hard'
            return 'soft'
        if T1[::-1] == T2:
            return 'soft'
        if (T1[0] == 'Z' and T2[0] == 'Z') or \
           (T1[1] == 'Z' and T2[1] == 'Z'):
            return 'soft'
        return 'hard'
    return 'soft'


def classify_empirical(H, N, gamma):
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


def naive_triple_classification(T1, T2, T3):
    """Compose pairwise classifications via:
      all truly → truly; any hard → hard; else soft."""
    cats = [classify_combinatorial_pair(T1, T2),
            classify_combinatorial_pair(T1, T3),
            classify_combinatorial_pair(T2, T3)]
    if all(c == 'truly' for c in cats):
        return 'truly'
    if 'hard' in cats:
        return 'hard'
    return 'soft'


def main():
    N = 3
    bonds = [(i, i + 1) for i in range(N - 1)]

    # 15 non-II bilinears
    paulis = ['I', 'X', 'Y', 'Z']
    bilinears = [t for t in product(paulis, repeat=2) if t != ('I', 'I')]
    assert len(bilinears) == 15

    # Enumerate triples with replacement, unordered
    triples = list(combinations_with_replacement(bilinears, 3))
    assert len(triples) == (17 * 16 * 15) // 6  # = 680
    print(f"Enumerated {len(triples)} unordered triples with replacement at N={N}")
    print()

    # Build H per bond from triple
    n_correct = 0
    cnt_emp = {'truly': 0, 'soft': 0, 'hard': 0}
    cnt_naive = {'truly': 0, 'soft': 0, 'hard': 0}
    mismatches_by_class = {'truly→soft': [], 'truly→hard': [],
                            'soft→truly': [], 'soft→hard': [],
                            'hard→truly': [], 'hard→soft': []}

    for triple in triples:
        T1, T2, T3 = triple
        terms = [(T1[0], T1[1], J), (T2[0], T2[1], J), (T3[0], T3[1], J)]
        H = fw._build_bilinear(N, bonds, terms)
        cat_emp = classify_empirical(H, N, GAMMA)
        cat_naive = naive_triple_classification(T1, T2, T3)
        cnt_emp[cat_emp] += 1
        cnt_naive[cat_naive] += 1
        if cat_emp == cat_naive:
            n_correct += 1
        else:
            label = (f"{T1[0]}{T1[1]}+{T2[0]}{T2[1]}+{T3[0]}{T3[1]}")
            key = f"{cat_naive}→{cat_emp}"
            if key in mismatches_by_class:
                mismatches_by_class[key].append(label)

    print(f"Empirical:    truly={cnt_emp['truly']}, soft={cnt_emp['soft']}, "
          f"hard={cnt_emp['hard']}")
    print(f"Naive (paired): truly={cnt_naive['truly']}, soft={cnt_naive['soft']}, "
          f"hard={cnt_naive['hard']}")
    print(f"Match: {n_correct}/{len(triples)} ({100*n_correct/len(triples):.0f}%)")
    print()
    print(f"Mismatch breakdown (naive→empirical):")
    for key, lst in mismatches_by_class.items():
        if lst:
            print(f"  {key}: {len(lst)} cases")
            for label in lst[:5]:
                print(f"    {label}")
            if len(lst) > 5:
                print(f"    ... ({len(lst)-5} more)")


if __name__ == "__main__":
    main()
