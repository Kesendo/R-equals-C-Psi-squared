#!/usr/bin/env python3
"""V-Effect topology sensitivity: chain vs ring vs star at N=4.

The combinatorial classifier (commit 079c7ce) gives 120/120 match at chain
N≥3. Testing other topologies at N=4 reveals the V-Effect's activation is
topology-sensitive:

  chain N=4: 15/46/59 (120/120 match — 2-bond chain adjacency activates)
  ring N=4:  15/48/57 (118/120 — chain + 2 cycle-closure extras)
  star N=4:  15/72/33 (94/120 = IDENTICAL to N=2 single-bond split!)

Star has 3 bonds but no chain-adjacency: all bonds share the center qubit,
not a 'pass-through' boundary qubit like chain. So the 26 V-Effect-trigger
cases stay soft on star.

This refines V-Effect's emergence: trigger requires CHAIN-ADJACENT bond
pairs (one bond ending at a qubit, another beginning there), not just two
bonds.
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
BPE = {('I', 'I'), ('X', 'X'), ('Y', 'Y'), ('Z', 'Z')}


def n_I(T):
    return T.count('I')


def active_letters(T):
    return [c for c in T if c != 'I']


def classify_combinatorial(T1, T2):
    """V-Effect classifier (chain-derived, commit 079c7ce)."""
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
    topologies = {
        'chain N=4': (4, [(0, 1), (1, 2), (2, 3)]),
        'ring N=4':  (4, [(0, 1), (1, 2), (2, 3), (3, 0)]),
        'star N=4':  (4, [(0, 1), (0, 2), (0, 3)]),
        'chain N=2 (baseline)': (2, [(0, 1)]),
    }

    print("V-Effect topology sensitivity at N=4 vs chain N=2 baseline")
    print(f"  Combinatorial classifier (chain-derived) gives 15/46/59 always.")
    print(f"  Empirical varies by topology — reveals V-Effect's activation rules.")
    print()
    print(f"  {'Topology':<25s}  {'truly':>5s}  {'soft':>5s}  {'hard':>5s}  {'match':>5s}")
    print('-' * 60)

    for top_name, (N, bonds) in topologies.items():
        n_correct = 0
        cnt = {'truly': 0, 'soft': 0, 'hard': 0}
        for pair in pairs:
            t1, t2 = pair
            terms = [(t1[0], t1[1], J), (t2[0], t2[1], J)]
            H = fw._build_bilinear(N, bonds, terms)
            cat_emp = classify_empirical(H, N, GAMMA)
            cat_comb = classify_combinatorial(t1, t2)
            cnt[cat_emp] += 1
            if cat_emp == cat_comb:
                n_correct += 1
        print(f"  {top_name:<25s}  {cnt['truly']:>5d}  {cnt['soft']:>5d}  "
              f"{cnt['hard']:>5d}  {n_correct:>3d}/120")

    print()
    print("Reading:")
    print("  - Chain N=4: full V-Effect activation (chain classifier match perfect).")
    print("  - Ring N=4: chain + 2 cycle-closure extras.")
    print("  - Star N=4: NO V-Effect activation despite 3 bonds — same as N=2.")
    print()
    print("  V-Effect requires CHAIN-ADJACENT bond pairs (one bond's endpoint =")
    print("  another bond's startpoint, NOT a shared center). Star's bonds all")
    print("  radiate from a hub, which lacks this structure.")
    print()
    print("  The 'atmosphere complexity' grows with TOPOLOGICAL ADJACENCY,")
    print("  not just bond count or N.")


if __name__ == "__main__":
    main()
