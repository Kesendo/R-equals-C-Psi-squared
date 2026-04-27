#!/usr/bin/env python3
"""V-Effect 120-enum combinatorial derivation (N≥3, including I-bilinears).

The 36-enum (no I) was fully derived in commit 81caf67. Extending to the
120-enum (with I allowed) at N=4: 15 truly / 46 soft / 59 hard. The
N=3-36-enum rule gives only 68% match in 120-enum — new rules needed for
single-letter bilinears.

This script tests an extended rule set:
  Layer 1: 36-enum rule (BPE + bit_a-partner + bond-flip + Z-align)
  Layer 2: both-single-letter rule (Pauli-letter-Π hierarchy)
  Layer 3: mixed (single + double) rule
"""
import sys
from itertools import product
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
BPE_extended = {('I', 'I'), ('X', 'X'), ('Y', 'Y'), ('Z', 'Z')}


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


def n_I(T):
    return T.count('I')


def active_letters(T):
    return [c for c in T if c != 'I']


def classify_combinatorial_v2(T1, T2):
    """Extended rule covering 120-enum (with I bilinears)."""
    n1, n2 = n_I(T1), n_I(T2)

    # Both single-letter (each has exactly one I)
    if n1 == 1 and n2 == 1:
        letters = active_letters(T1) + active_letters(T2)
        unique_letters = set(letters)
        if 'Z' in unique_letters:
            return 'hard'
        if unique_letters == {'X'}:
            return 'truly'
        return 'soft'  # has Y, no Z

    # Mixed: one is single-letter, the other is double-letter
    if (n1 == 1) != (n2 == 1):
        single = T1 if n1 == 1 else T2
        double = T2 if n1 == 1 else T1
        single_letter = active_letters(single)[0]

        # Pauli-letter-Π hierarchy:
        # X-single: Π-trivial (X → I), can saturate truly
        # Y-single: Π-half (Y → iZ), maxes at soft
        # Z-single: Π-full-conflict, always hard
        if single_letter == 'Z':
            return 'hard'

        if single_letter == 'X':
            # + BPE → truly
            if double in BPE_extended:
                return 'truly'
            # + Y-Z pair (commute with Z-dephasing structure) → soft
            if double in {('Y', 'Z'), ('Z', 'Y')}:
                return 'soft'
            # everything else → hard
            return 'hard'

        if single_letter == 'Y':
            # Y + BPE → soft
            if double in BPE_extended:
                return 'soft'
            # Y + X-Z pair → soft
            if double in {('X', 'Z'), ('Z', 'X')}:
                return 'soft'
            # everything else (X-Y, Y-Z pairs) → hard
            return 'hard'

    # Both double-letter
    if n1 == 0 and n2 == 0:
        in_BPE_1 = T1 in BPE_extended
        in_BPE_2 = T2 in BPE_extended
        if in_BPE_1 and in_BPE_2:
            return 'truly'
        if T1 == T2:
            # Self-pair: BPE → truly, non-BPE → soft
            if T1 in BPE_extended:
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
        # Both non-BPE
        if T1[::-1] == T2:
            return 'soft'
        if (T1[0] == 'Z' and T2[0] == 'Z') or \
           (T1[1] == 'Z' and T2[1] == 'Z'):
            return 'soft'
        return 'hard'

    return 'soft'  # placeholder


def enumerate_120():
    paulis = ['I', 'X', 'Y', 'Z']
    seen = set()
    out = []
    for t1 in product(paulis, repeat=2):
        for t2 in product(paulis, repeat=2):
            if t1 == ('I', 'I') or t2 == ('I', 'I'):
                continue
            sorted_t = tuple(sorted([t1, t2]))
            if sorted_t in seen:
                continue
            seen.add(sorted_t)
            out.append(sorted_t)
    return out


def main():
    pairs = enumerate_120()
    print(f"V-Effect 120-enum combinatorial test")
    print(f"  γ = {GAMMA}, total pairs: {len(pairs)}")
    print()

    for N in [3, 4]:
        bonds = [(i, i + 1) for i in range(N - 1)]
        n_correct = 0
        cnt_emp = {'truly': 0, 'soft': 0, 'hard': 0}
        cnt_comb = {'truly': 0, 'soft': 0, 'hard': 0}
        mismatches_by_class = {'both_single': [], 'mixed': [], 'both_double': []}

        for pair in pairs:
            t1, t2 = pair
            terms = [(t1[0], t1[1], J), (t2[0], t2[1], J)]
            H = fw._build_bilinear(N, bonds, terms)
            cat_emp = classify_empirical(H, N, GAMMA)
            cat_comb = classify_combinatorial_v2(t1, t2)
            cnt_emp[cat_emp] += 1
            cnt_comb[cat_comb] += 1
            if cat_emp == cat_comb:
                n_correct += 1
            else:
                # Classify the pair structure
                n1, n2 = n_I(t1), n_I(t2)
                if n1 == 1 and n2 == 1:
                    cls = 'both_single'
                elif n1 == 0 and n2 == 0:
                    cls = 'both_double'
                else:
                    cls = 'mixed'
                label = f"{t1[0]}{t1[1]}+{t2[0]}{t2[1]}"
                mismatches_by_class[cls].append((label, cat_emp, cat_comb))

        print(f"N={N}, 120-enum:")
        print(f"  Empirical:    truly={cnt_emp['truly']}, soft={cnt_emp['soft']}, hard={cnt_emp['hard']}")
        print(f"  Combinatorial: truly={cnt_comb['truly']}, soft={cnt_comb['soft']}, hard={cnt_comb['hard']}")
        print(f"  Match: {n_correct}/{len(pairs)} ({100 * n_correct / len(pairs):.0f}%)")
        for cls in ['both_single', 'mixed', 'both_double']:
            n_mis = len(mismatches_by_class[cls])
            print(f"    {cls} mismatches: {n_mis}")
        print()
        for cls in ['both_single', 'mixed', 'both_double']:
            mis = mismatches_by_class[cls]
            if mis:
                print(f"  {cls} mismatches (first 8):")
                for label, emp, comb in mis[:8]:
                    print(f"    {label}: empirical={emp}, combinatorial={comb}")
                print()


if __name__ == "__main__":
    main()
