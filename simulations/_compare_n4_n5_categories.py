#!/usr/bin/env python3
"""Compare which Pauli-pair Hamiltonians are soft / truly / hard at N=4 vs N=5.

Tom's hypothesis: the half-integer mirror at odd N (N=3, N=5) introduces
breaking that the integer mirror at even N (N=2, N=4) does not. If true,
the SPECIFIC Hamiltonians in the soft category should differ between
N=4 and N=5, even though the COUNTS happen to match (15/46/59 in both).

This script enumerates 120 unordered Pauli-pair Hamiltonians, classifies
each at N=4 and N=5 (truly / soft / hard), and reports:
  - intersection: same category at both N
  - symmetric difference: category changes between N=4 and N=5

If the difference is small, the trichotomy is mirror-parity-stable.
If the difference is large, the mirror parity matters.
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


def spectrum_pair_max_err(L, sigma_gamma):
    evals = np.linalg.eigvals(L)
    used = np.zeros(len(evals), dtype=bool)
    max_err = 0.0
    for i in range(len(evals)):
        if used[i]:
            continue
        target = -evals[i] - 2 * sigma_gamma
        dists = np.abs(evals - target)
        for j in range(len(evals)):
            if used[j]:
                dists[j] = np.inf
        best_j = int(np.argmin(dists))
        if best_j != i:
            used[i] = True
            used[best_j] = True
        else:
            used[i] = True
        max_err = max(max_err, float(dists[best_j]))
    return max_err


def classify_at_N(N, terms, gamma=0.1, J=1.0):
    bonds = [(i, i + 1) for i in range(N - 1)]
    bilinear = [(t[0], t[1], J) for t in terms]
    H = fw._build_bilinear(N, bonds, bilinear)
    L = fw.lindbladian_z_dephasing(H, [gamma] * N)
    M = fw.palindrome_residual(L, N * gamma, N)
    op_norm = float(np.linalg.norm(M))
    spec_err = spectrum_pair_max_err(L, N * gamma)
    op_ok = op_norm < 1e-10
    spec_ok = spec_err < 1e-6
    if op_ok:
        return 'truly'
    if spec_ok:
        return 'soft'
    return 'hard'


def enumerate_pairs():
    """Same enumeration as _pi_protected_test_n4.py: 120 unordered."""
    paulis = ['I', 'X', 'Y', 'Z']
    seen = set()
    pairs = []
    for term1 in product(paulis, repeat=2):
        for term2 in product(paulis, repeat=2):
            if term1 == ('I', 'I') or term2 == ('I', 'I'):
                continue
            sorted_terms = tuple(sorted([term1, term2]))
            if sorted_terms in seen:
                continue
            seen.add(sorted_terms)
            label = f"{term1[0]}{term1[1]}+{term2[0]}{term2[1]}"
            pairs.append((sorted_terms, label, [term1, term2]))
    return pairs


def main():
    pairs = enumerate_pairs()
    print(f"Classifying {len(pairs)} Hamiltonians at N=4 (integer mirror) and N=5 (half-integer)")
    print()

    classifications = {}  # label -> {'N=4': cat, 'N=5': cat}
    for sorted_terms, label, terms in pairs:
        cat4 = classify_at_N(4, terms)
        cat5 = classify_at_N(5, terms)
        classifications[label] = {'N=4': cat4, 'N=5': cat5}

    # Tabulate transitions
    transitions = {}  # (cat_n4, cat_n5) -> list of labels
    for label, c in classifications.items():
        key = (c['N=4'], c['N=5'])
        transitions.setdefault(key, []).append(label)

    print(f"  {'N=4':>10s} → {'N=5':>10s}    {'count':>5s}    examples")
    print(f"  {'-' * 10}   {'-' * 10}    {'-' * 5}    {'-' * 50}")
    for key in sorted(transitions.keys()):
        labels = transitions[key]
        same = "(same)" if key[0] == key[1] else "(SHIFT)"
        examples = ", ".join(labels[:5])
        if len(labels) > 5:
            examples += f", ... ({len(labels) - 5} more)"
        print(f"  {key[0]:>10s} → {key[1]:>10s}    {len(labels):>5d} {same}   {examples}")

    print()
    print(f"Summary:")
    n_stable = sum(len(v) for k, v in transitions.items() if k[0] == k[1])
    n_shifted = sum(len(v) for k, v in transitions.items() if k[0] != k[1])
    print(f"  stable across N=4 → N=5:  {n_stable}/{len(pairs)}")
    print(f"  category shifted:          {n_shifted}/{len(pairs)}")

    # Specifically: which soft at N=4 became NOT soft at N=5? And vice versa?
    soft_n4 = {l for l, c in classifications.items() if c['N=4'] == 'soft'}
    soft_n5 = {l for l, c in classifications.items() if c['N=5'] == 'soft'}
    only_soft_n4 = soft_n4 - soft_n5
    only_soft_n5 = soft_n5 - soft_n4
    soft_both = soft_n4 & soft_n5
    print()
    print(f"Soft-category overlap:")
    print(f"  soft at both N=4 and N=5:  {len(soft_both)}")
    print(f"  soft only at N=4:           {len(only_soft_n4)}  ({sorted(only_soft_n4)})")
    print(f"  soft only at N=5:           {len(only_soft_n5)}  ({sorted(only_soft_n5)})")


if __name__ == "__main__":
    main()
