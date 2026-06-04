#!/usr/bin/env python3
"""F87 task B(b): unify the §6 counting rule with the bipartite obstruction (2026-06-04).

The §6 diagonal-cell hardness rule (PROOF_F103 §6): a k=3 pair is hard iff
  (a) at least one term is a pure-D template (all-diagonal), OR
  (b) both terms are single-diagonal with their {I,D} letter at chain-adjacent positions.
The §7 mechanism: hard iff the pair is NON-bipartite in the dephasing basis. A graph is
non-bipartite two ways: a DIAGONAL LIFT (a diagonal H-entry, no chiral K possible) or an
ODD CYCLE in the hopping graph. This script verifies the exact correspondence

  rule (a)  ⟺  diagonal lift            (template term ⟹ diagonal H, has_diag)
  rule (b)  ⟺  odd cycle (K3 triangle)  (adjacent single-diagonal ⟹ triangle in S)
  §6-hard (a∨b)  ⟺  non-bipartite  ⟺  F87-hard

over ALL diagonal-cell (Z-deph, Klein (0,1)) y_par-homogeneous pairs (templates included).
This closes the combinatorial half of the windowed converse: the §6 rule IS the bipartite
obstruction, split into its two mechanisms. Only the spectral half (odd cycle ⟹ the
degenerate D̂-block is asymmetric ⟹ hard) then remains.
"""
from __future__ import annotations

import sys
from itertools import product, combinations_with_replacement
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
from framework.pauli import _build_kbody_chain
from f87_flip_generators import is_bipartite

DIAG = {"I", "Z"}


def n_diag(t):
    return sum(1 for L in t if L in DIAG)


def is_template(t):
    return n_diag(t) == len(t) and not all(L == "I" for L in t)  # all-diagonal, not all-I


def is_single_diag(t):
    return n_diag(t) == 1


def diag_pos(t):
    return [i for i, L in enumerate(t) if L in DIAG][0]


def main():
    N, k = 4, 3
    chain = fw.ChainSystem(N=N)
    terms = [t for t in product("IXYZ", repeat=k)
             if not all(L == "I" for L in t) and fw.klein_index(t) == (0, 1)]

    print("=" * 96)
    print("F87 B(b): §6 rule  ⟺  bipartite obstruction  ⟺  F87-hard   (N=4, k=3, Z-deph (0,1))")
    print("=" * 96)
    print()

    tol = 1e-9
    rows = []
    for t1, t2 in combinations_with_replacement(terms, 2):
        if (sum(c == "Y" for c in t1) % 2) != (sum(c == "Y" for c in t2) % 2):
            continue  # y_par-homogeneous
        H = _build_kbody_chain(N, [tuple(t1) + (1.0,), tuple(t2) + (1.0,)])
        has_diag = bool(np.max(np.abs(np.diag(H))) > tol)
        non_bip = not is_bipartite(H)                 # is_bipartite returns False if has_diag OR odd cycle
        has_odd_cycle = non_bip and not has_diag
        # §6 rule
        rule_a = is_template(t1) or is_template(t2)
        rule_b = (is_single_diag(t1) and is_single_diag(t2)
                  and abs(diag_pos(t1) - diag_pos(t2)) == 1)
        rule_hard = rule_a or rule_b
        # ground truth
        cls = fw.classify_pauli_pair(chain, [t1, t2], dephase_letter="Z")
        hard = (cls == "hard")
        rows.append((rule_a, rule_b, rule_hard, has_diag, has_odd_cycle, non_bip, hard))

    n = len(rows)
    eq_a = sum(1 for r in rows if r[0] == r[3])              # rule_a ⟺ has_diag
    eq_b = sum(1 for r in rows if r[1] == r[4])              # rule_b ⟺ has_odd_cycle
    eq_hard1 = sum(1 for r in rows if r[2] == r[5])          # §6-hard ⟺ non_bip
    eq_hard2 = sum(1 for r in rows if r[2] == r[6])          # §6-hard ⟺ F87-hard
    eq_hard3 = sum(1 for r in rows if r[5] == r[6])          # non_bip ⟺ F87-hard

    print(f"  diagonal-cell y_par-homogeneous pairs: {n}")
    print(f"  rule (a)  ⟺  diagonal lift:        {eq_a}/{n}  {'ALL' if eq_a == n else 'MISMATCH'}")
    print(f"  rule (b)  ⟺  odd cycle (triangle): {eq_b}/{n}  {'ALL' if eq_b == n else 'MISMATCH'}")
    print(f"  §6-hard (a∨b)  ⟺  non-bipartite:   {eq_hard1}/{n}  {'ALL' if eq_hard1 == n else 'MISMATCH'}")
    print(f"  §6-hard (a∨b)  ⟺  F87-hard:        {eq_hard2}/{n}  {'ALL' if eq_hard2 == n else 'MISMATCH'}")
    print(f"  non-bipartite  ⟺  F87-hard:        {eq_hard3}/{n}  {'ALL' if eq_hard3 == n else 'MISMATCH'}")
    assert eq_a == n and eq_b == n and eq_hard1 == n and eq_hard2 == n and eq_hard3 == n

    # mechanism breakdown of the hard pairs (reproduce the §6 counting)
    n_hard = sum(1 for r in rows if r[6])
    via_lift = sum(1 for r in rows if r[6] and r[3])         # hard via diagonal lift
    via_cycle = sum(1 for r in rows if r[6] and r[4])        # hard via odd cycle
    via_both = sum(1 for r in rows if r[6] and r[3] and r[4])
    print()
    print(f"  hard pairs: {n_hard}   via diagonal lift (rule a): {via_lift}   "
          f"via odd cycle (rule b): {via_cycle}   (both: {via_both})")
    print()
    print("  => the §6 rule IS the bipartite obstruction, split into its two mechanisms:")
    print("     rule (a) template = diagonal lift, rule (b) adjacency = odd-cycle triangle.")
    print("     Combinatorial half of the windowed converse: CLOSED (bit-exact, this cell).")
    print("     Remaining: the spectral half, odd cycle ⟹ degenerate D̂-block asymmetric ⟹ hard.")


if __name__ == "__main__":
    main()
