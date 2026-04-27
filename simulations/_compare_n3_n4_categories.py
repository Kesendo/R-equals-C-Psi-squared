#!/usr/bin/env python3
"""EQ-031: N=3 → N=4 categorical mapping across the cusp.

For each of the 120 unordered two-term Pauli-pair Hamiltonians (same
enumeration as _pi_protected_test_n4.py and _compare_n4_n5_categories.py),
classify at N=3 (half-integer mirror w_XY = 1.5, no modes on it) and at
N=4 (integer mirror w_XY = 2, modes on it), then tabulate transitions.

If categories are preserved across the cusp: the cusp is structure-fixing
in identity (each Hamiltonian keeps its category through the transition).
If they shuffle: the cusp permutes categories.

This is the operator-level analogue of "what does crossing CΨ = 1/4 do?"
asked at the level of the parity-class structure.
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
    print(f"Classifying {len(pairs)} Hamiltonians at N=3 (cusp entry, half-int mirror) "
          f"and N=4 (cusp exit, int mirror)")
    print()

    classifications = {}
    for sorted_terms, label, terms in pairs:
        cat3 = classify_at_N(3, terms)
        cat4 = classify_at_N(4, terms)
        classifications[label] = {'N=3': cat3, 'N=4': cat4}

    # Tabulate transitions
    transitions = {}
    for label, c in classifications.items():
        key = (c['N=3'], c['N=4'])
        transitions.setdefault(key, []).append(label)

    print(f"  {'N=3':>10s} → {'N=4':>10s}    {'count':>5s}    examples")
    print(f"  {'-' * 10}   {'-' * 10}    {'-' * 5}    {'-' * 50}")
    for key in sorted(transitions.keys()):
        labels = transitions[key]
        same = "(same)" if key[0] == key[1] else "(SHIFT)"
        examples = ", ".join(labels[:6])
        if len(labels) > 6:
            examples += f", ... ({len(labels) - 6} more)"
        print(f"  {key[0]:>10s} → {key[1]:>10s}    {len(labels):>5d} {same:>8s}   {examples}")

    # Per-category summaries
    print()
    counts_n3 = {'truly': 0, 'soft': 0, 'hard': 0}
    counts_n4 = {'truly': 0, 'soft': 0, 'hard': 0}
    for label, c in classifications.items():
        counts_n3[c['N=3']] += 1
        counts_n4[c['N=4']] += 1
    print(f"At N=3 (with full 120 enumeration): "
          f"{counts_n3['truly']} truly, {counts_n3['soft']} soft, {counts_n3['hard']} hard")
    print(f"At N=4: "
          f"{counts_n4['truly']} truly, {counts_n4['soft']} soft, {counts_n4['hard']} hard")
    print()

    # Stable vs shifted
    n_stable = sum(len(v) for k, v in transitions.items() if k[0] == k[1])
    n_shifted = sum(len(v) for k, v in transitions.items() if k[0] != k[1])
    print(f"stable across N=3 → N=4:  {n_stable}/{len(pairs)}")
    print(f"category shifted:          {n_shifted}/{len(pairs)}")

    # Detailed shifts
    print()
    print("Specific category shifts at the cusp:")
    for key in sorted(transitions.keys()):
        if key[0] == key[1]:
            continue
        labels = transitions[key]
        print(f"  {key[0]:>10s} → {key[1]:<10s}  ({len(labels):>3d}): {labels}")


if __name__ == "__main__":
    main()
