#!/usr/bin/env python3
"""Where the F103/F105 k=3 sub-cell record starts: the lower bound in N.

F103 anchored the five Z2^3 sub-cell records at N=4, F105 reproduced them
bit-exactly at N=5, and both then carried the reading "N-stable" with no
lower bound. F103's own §7.4 supplies one: rule (b)'s odd cycle needs a term
placed at more than one window, i.e. k < N. At k=3, N=3 is full support, so
the adjacency half of the §6 counting rule is empty there.

This gate measures the boundary instead of arguing it. Same enumeration as
`f87_z2cubed_split_n4_k3.py` (the same
294 unordered PAIRS of length-3 templates over {I, X, Y, Z}, Klein- and
y_par-homogeneous, self-pairs included), classified on the N=3 and N=4 chain:

  * three records are already in place at N=3: truly stays y_par=0-pure,
    mother soft stays 0:21, and all six off-diagonal Pattern B / Pattern C
    cells are unchanged;
  * the two diagonal-cell records are not: N=3 reads 34:0 hard and 21:21
    soft against N=4's 42:8 and 13:13. The 34 is F103 §6's own template
    count, the split with the adjacency contribution removed.

The N=4 half is a live control, not decoration: it fails if the pipeline
stops computing, which no assertion on the N=3 differences alone would catch.

Above N=5 the spectral classifier gets expensive (it eigendecomposes a
4^N x 4^N Liouvillian), so the upper reach is measured through F103 §7's
criterion instead: a diagonal-cell pair is soft iff H's hopping graph is
bipartite in the dephase letter's eigenbasis, which reads the 2^N x 2^N
Hamiltonian only. Two guards keep that from being a free pass. The criterion
is checked against the spectral classifier PAIR BY PAIR over all three
diagonal cells at N=3 and N=4, where both are affordable; and at N=6 four
representative pairs are classified spectrally and have to match it.

Runtime ~5 min on an idle machine.
"""
from __future__ import annotations

import sys
import time
from collections import defaultdict, deque
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw  # noqa: E402

from f87_z2cubed_split_n4_k3 import enumerate_z2_homogeneous_k3  # noqa: E402
from framework.pauli import _build_kbody_chain  # noqa: E402

# How far F103 §7's criterion is swept, and the four pairs cross-checked
# spectrally at N=6: one for rule (a), two for rule (b) (one built from X
# letters and one from Y, so the check is not blind to the Y content), and
# one that neither rule fires on.
CRITERION_NMAX = 8
N6_CROSS_CHECK = [
    (('Z', 'I', 'I'), ('I', 'Z', 'I'), 'hard'),   # rule (a), a pure-Z template
    (('Z', 'X', 'X'), ('X', 'Z', 'X'), 'hard'),   # rule (b), |Δpos| = 1
    (('Z', 'X', 'X'), ('X', 'X', 'Z'), 'soft'),   # neither, |Δpos| = 2
    (('Y', 'Z', 'Y'), ('Y', 'Y', 'Z'), 'hard'),   # rule (b) again, with Y content
]

KLEINS = [(0, 0), (0, 1), (1, 0), (1, 1)]
LETTERS = ('Z', 'X', 'Y')
CLASSES = ('truly', 'soft', 'hard')

# The diagonal Klein cell of each dephase letter (F103 Notation block).
DIAGONAL_CELL = {'Z': (0, 1), 'X': (1, 0), 'Y': (1, 1)}

# F103's frozen N=4 records, per dephase letter, as (y_par=0, y_par=1).
N4_DIAGONAL_HARD = {'Z': (42, 8), 'X': (42, 8), 'Y': (8, 42)}
N4_DIAGONAL_SOFT = {'Z': (13, 13), 'X': (13, 13), 'Y': (13, 13)}

# The same two cells at N=3, with the adjacency contribution absent.
N3_DIAGONAL_HARD = {'Z': (34, 0), 'X': (34, 0), 'Y': (0, 34)}
N3_DIAGONAL_SOFT = {'Z': (21, 21), 'X': (21, 21), 'Y': (21, 21)}

# The three records that are already in place at N=3, shared by both chains.
MOTHER_SOFT = (0, 21)
PATTERN_B_C = {
    ((0, 1), 'Y'): (55, 21),
    ((1, 1), 'Z'): (21, 55),
    ((1, 1), 'X'): (21, 55),
    ((0, 1), 'X'): (0, 21),
    ((1, 0), 'Z'): (0, 21),
    ((1, 0), 'Y'): (0, 21),
}


# Single-qubit V with V Z V† = D, so V† (·) V reads an operator in D's
# eigenbasis. §7's criterion is stated in that basis, not the computational
# one; skipping the rotation makes every X- and Y-dephasing cell read soft.
_TO_DEPHASING_BASIS = {
    'Z': np.eye(2, dtype=complex),
    'X': np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2),
    'Y': np.array([[1, 1], [1j, -1j]], dtype=complex) / np.sqrt(2),
}


def in_dephasing_basis(H, N, letter):
    """V^⊗N† H V^⊗N: H read in the dephase letter's eigenbasis."""
    V = _TO_DEPHASING_BASIS[letter]
    V_all = np.array([[1.0 + 0j]])
    for _ in range(N):
        V_all = np.kron(V_all, V)
    return V_all.conj().T @ H @ V_all


def hopping_graph_is_bipartite(H, tol=1e-12):
    """F103 §7's criterion, read on H: soft iff G_H is 2-colourable.

    A nonzero diagonal is the rule (a) obstruction (no diagonal K can flip a
    diagonal entry's sign); an odd cycle is the rule (b) one. H must already
    be expressed in the dephase letter's eigenbasis (`in_dephasing_basis`).
    """
    if np.any(np.abs(np.diag(H)) > tol):
        return False
    adjacency = csr_matrix(np.abs(H) > tol)
    indptr, indices = adjacency.indptr, adjacency.indices
    colour = np.full(adjacency.shape[0], -1, dtype=np.int8)
    for source in range(adjacency.shape[0]):
        if colour[source] >= 0:
            continue
        colour[source] = 0
        queue = deque([source])
        while queue:
            u = queue.popleft()
            for v in indices[indptr[u]:indptr[u + 1]]:
                if colour[v] < 0:
                    colour[v] = 1 - colour[u]
                    queue.append(v)
                elif colour[v] == colour[u]:
                    return False
    return True


def criterion_verdict(N, terms, letter):
    """§7's soft/hard verdict for one pair on the N-chain."""
    H = _build_kbody_chain(N, [tuple(t) + (1.0,) for t in terms])
    return 'soft' if hopping_graph_is_bipartite(
        in_dephasing_basis(H, N, letter)) else 'hard'


def criterion_cell(N, pairs, letter):
    """(hard, soft) counts by y_par in one diagonal cell, via the criterion."""
    counts = defaultdict(int)
    for terms, y_par in pairs:
        counts[(criterion_verdict(N, terms, letter), y_par)] += 1
    return counts


def classify_grid(N, items):
    """counts[(klein, letter, y_par, class)] over the whole enumeration."""
    chain = fw.ChainSystem(N=N)
    counts = defaultdict(int)
    for terms, klein, y_par in items:
        for letter in LETTERS:
            cls = fw.classify_pauli_pair(chain, terms, dephase_letter=letter)
            counts[(klein, letter, y_par, cls)] += 1
    return counts


def cell(counts, klein, letter, cls):
    return (counts[(klein, letter, 0, cls)], counts[(klein, letter, 1, cls)])


def report(N, counts):
    print(f"--- N={N} ({N - 2} sliding window(s) for a length-3 template) ---")
    for cls in CLASSES:
        print(f"  {cls}:")
        for klein in KLEINS:
            cells = "  ".join(
                f"{letter}: {cell(counts, klein, letter, cls)}" for letter in LETTERS
            )
            print(f"    Klein {klein}   {cells}")


def main():
    items = enumerate_z2_homogeneous_k3()
    assert len(items) == 294, f"enumeration drifted: {len(items)} pairs, expected 294"

    start = time.time()
    grids = {N: classify_grid(N, items) for N in (3, 4)}
    print(f"=== F103/F105 k=3 sub-cell records at N=3 and N=4 "
          f"({time.time() - start:.0f}s) ===\n")
    for N in (3, 4):
        report(N, grids[N])
        print()

    checks = 0

    # The control: N=4 still reproduces F103's frozen anchor.
    for letter in LETTERS:
        klein = DIAGONAL_CELL[letter]
        got = cell(grids[4], klein, letter, 'hard')
        assert got == N4_DIAGONAL_HARD[letter], (letter, 'N=4 hard', got)
        got = cell(grids[4], klein, letter, 'soft')
        assert got == N4_DIAGONAL_SOFT[letter], (letter, 'N=4 soft', got)
        checks += 2

    # The boundary: the two diagonal-cell records read differently at N=3.
    for letter in LETTERS:
        klein = DIAGONAL_CELL[letter]
        got = cell(grids[3], klein, letter, 'hard')
        assert got == N3_DIAGONAL_HARD[letter], (letter, 'N=3 hard', got)
        assert got != N4_DIAGONAL_HARD[letter], (letter, 'N=3 hard did not move')
        got = cell(grids[3], klein, letter, 'soft')
        assert got == N3_DIAGONAL_SOFT[letter], (letter, 'N=3 soft', got)
        assert got != N4_DIAGONAL_SOFT[letter], (letter, 'N=3 soft did not move')
        checks += 4

    # The three records that are already in place at N=3, at both chain lengths.
    for N in (3, 4):
        truly_y1 = sum(counts for (_, _, y_par, cls), counts in grids[N].items()
                       if y_par == 1 and cls == 'truly')
        truly_total = sum(counts for (_, _, _, cls), counts in grids[N].items()
                          if cls == 'truly')
        assert truly_y1 == 0, (N, 'truly is not y_par=0-pure', truly_y1)
        assert truly_total == 300, (N, 'truly total', truly_total)
        checks += 2

        for letter in LETTERS:
            got = cell(grids[N], (0, 0), letter, 'soft')
            assert got == MOTHER_SOFT, (N, letter, 'mother soft', got)
            checks += 1

        for (klein, letter), expected in PATTERN_B_C.items():
            got = cell(grids[N], klein, letter, 'soft')
            assert got == expected, (N, klein, letter, 'off-diagonal soft', got)
            checks += 1

    # The upper reach, through §7's criterion. It has to agree with the
    # spectral verdicts above at BOTH N, or it is not measuring the same thing.
    print(f"=== F103 §7 criterion over the diagonal cell, N=3..{CRITERION_NMAX} ===")
    start = time.time()
    for letter in LETTERS:
        klein = DIAGONAL_CELL[letter]
        pairs = [(terms, y_par) for terms, k, y_par in items if k == klein]
        assert len(pairs) == 76, (letter, 'diagonal cell size', len(pairs))
        for N in range(3, CRITERION_NMAX + 1):
            counts = criterion_cell(N, pairs, letter)
            hard = (counts[('hard', 0)], counts[('hard', 1)])
            soft = (counts[('soft', 0)], counts[('soft', 1)])
            print(f"  {letter}-deph  N={N}  hard {hard}  soft {soft}")
            expected_hard = N3_DIAGONAL_HARD if N == 3 else N4_DIAGONAL_HARD
            expected_soft = N3_DIAGONAL_SOFT if N == 3 else N4_DIAGONAL_SOFT
            assert hard == expected_hard[letter], (letter, N, 'hard', hard)
            assert soft == expected_soft[letter], (letter, N, 'soft', soft)
            checks += 2
    print(f"  ({time.time() - start:.0f}s)\n")

    # Counts agreeing is weaker than the criterion agreeing: two verdicts could
    # swap and leave the totals alone. Pin it pair by pair where both are cheap.
    print("=== criterion vs spectral, pair by pair, over the diagonal cells ===")
    start = time.time()
    for N in (3, 4):
        chain = fw.ChainSystem(N=N)
        mismatches = 0
        for letter in LETTERS:
            klein = DIAGONAL_CELL[letter]
            for terms, k, _ in items:
                if k != klein:
                    continue
                spectral = fw.classify_pauli_pair(chain, terms, dephase_letter=letter)
                if criterion_verdict(N, terms, letter) != spectral:
                    mismatches += 1
        print(f"  N={N}  228 pairs (76 per letter)  mismatches = {mismatches}")
        assert mismatches == 0, (N, 'criterion disagrees with spectral', mismatches)
        checks += 1
    print(f"  ({time.time() - start:.0f}s)\n")

    # And the criterion is pinned to the spectral classifier once above N=5.
    print("=== N=6 spectral cross-check of the criterion ===")
    start = time.time()
    chain6 = fw.ChainSystem(N=6)
    for term1, term2, expected in N6_CROSS_CHECK:
        by_criterion = criterion_verdict(6, [term1, term2], 'Z')
        spectral = fw.classify_pauli_pair(chain6, [term1, term2], dephase_letter='Z')
        name = f"{''.join(term1)}+{''.join(term2)}"
        print(f"  {name}  criterion {by_criterion}  spectral {spectral}")
        assert by_criterion == expected, (name, 'criterion', by_criterion)
        assert spectral == expected, (name, 'spectral', spectral)
        checks += 2
    print(f"  ({time.time() - start:.0f}s)\n")

    print(f"{checks} checks pass.")
    print("The three alphabet-driven records hold at N=3 and N=4; the two")
    print("diagonal-cell records hold from N=4 (F103 §7.4: rule (b) needs k < N)")
    print(f"and, by F103 §7's criterion, on through N={CRITERION_NMAX}.")


if __name__ == "__main__":
    main()
