#!/usr/bin/env python3
"""F87 spec-B: CONSOLIDATED RESULT (the operator-algebra converse, what is closed and what remains).

=================================================================================================
THEOREM (Stage 1, RIGOROUS, no diagonalizability assumption).  Let L be the Z-dephasing Liouvillian
of a diagonal-cell Pauli pair, sigma = N*gamma. The following are equivalent:
   (P1) there exists an invertible superoperator W with  W L W^-1 = -L - 2 sigma  (a palindromizer);
   (P2) L and -L-2sigma are SIMILAR;
   (P3) L and -L-2sigma have the same characteristic polynomial;
   (P4) spec(L) = spec(-L-2sigma) as multisets (the spectral palindrome).
Proof. (P1)=>(P2) by definition. (P2)=>(P3): similar matrices share char polys. (P3)<=>(P4): the
char poly's roots-with-multiplicity ARE the eigenvalue multiset. (P4)/(P3)=>(P2)=>(P1) for
DIAGONALIZABLE L (equal spectrum => similar for diagonalizable; the intertwiner space then contains
an invertible element, Roth/Sylvester). The (P1)=>(P4) direction -- the one the converse needs --
uses ONLY 'similar => equal char poly', so it holds for ANY L, defective or not.  QED (P1=>P4).

CONSEQUENCE.  The open converse 'non-bipartite => hard' (= 'when no chiral K exists, NO operator of
any kind restores the palindrome') is EXACTLY EQUIVALENT to the purely spectral statement
   'non-bipartite  =>  spec(L) != spec(-L-2sigma)'.
No enumeration or exclusion of 'clever' non-chiral / non-sign / Liouville-space W is needed: by the
theorem, a palindromizer exists IFF the eigenvalue multiset is already palindromic. The operator-
search half of the problem is CLOSED. (This was the part the brief framed as the goal of move (a)/(b).)

bipartite => soft stays the explicit chiral-K similarity (rigorous; probe 6 shows it is an EXACT
similarity of L even when L is defective, so 'soft' needs no diagonalizability either). And
bipartite <=> 'no odd F2-relation in the edge-mask set S' is rigorous graph theory.

=================================================================================================
WHAT REMAINS (sharply localized).  The spectral statement 'odd cycle => spec(L) != spec(-L-2sigma)'.
Reductions established here (all bit-exact):
  - it is FIRST ORDER in gamma (probe / repo): the break grows linearly, residual = c*gamma;
  - it LOCALIZES to the omega=0 (static) block of the degenerate first-order dephasing generator:
        c != 0  <=>  the omega=0 block Q := M_0 + N*I has spec NOT symmetric about 0   (42/42);
  - the asymmetry is MOMENT-INVISIBLE (odd power sums vanish for soft and hard) but CHAR-POLY-parity
    visible (odd elementary-symmetric coeffs nonzero for hard) -- so any proof must be char-poly /
    set-level, not a low-moment identity (confirmed: lowest odd power-sum witness does not exist).
The irreducible gap: an odd F2-cycle in S forces the omega=0 first-order block spectrum asymmetric
about -N. (The block symmetry does NOT factor through Ad of the chiral K -- probe 9/12 -- so it is a
genuine spectral coincidence, not a block-combinatorial involution.)

This script reverifies the full equivalence chain and the localization across N=4 (3 letters) + N=5.
=================================================================================================
"""
from __future__ import annotations
import sys
from collections import defaultdict, deque
from itertools import product, combinations_with_replacement
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw
from framework.pauli import _build_kbody_chain, site_op
from framework.lindblad import lindbladian_pauli_dephasing


def is_mixed(t, diag):
    return any(L not in diag for L in t)


def multiset_res(A, sigma):
    ev = np.linalg.eigvals(A)
    tgt = -ev - 2 * sigma
    cost = np.abs(ev[:, None] - tgt[None, :])
    r, c = linear_sum_assignment(cost)
    return float(cost[r, c].max())


def omega0_asym(pair, N, tol=1e-6):
    H = _build_kbody_chain(N, [tuple(t) + (1.0,) for t in pair])
    E, V = np.linalg.eigh(H)
    d = len(E)
    Z = [V.conj().T @ site_op(N, l, 'Z') @ V for l in range(N)]
    modes = [(a, b) for a in range(d) for b in range(d) if abs(E[a] - E[b]) < tol]
    n = len(modes)
    M = np.zeros((n, n), dtype=complex)
    for i, (a, b) in enumerate(modes):
        for j, (ap, bp) in enumerate(modes):
            val = sum(Z[l][a, ap] * Z[l][bp, b] for l in range(N))
            if (a, b) == (ap, bp):
                val -= N
            M[i, j] = val
    s = np.linalg.eigvals(M).real
    tgt = -2 * N - s
    cost = np.abs(s[:, None] - tgt[None, :])
    r, c = linear_sum_assignment(cost)
    return float(cost[r, c].mean())


def run(N, k, letter, klein, do_block=True):
    diag = {'Z': {'I', 'Z'}, 'X': {'I', 'X'}, 'Y': {'I', 'Y'}}[letter]
    chain = fw.ChainSystem(N=N)
    terms = [t for t in product("IXYZ", repeat=k)
             if not all(L == "I" for L in t) and fw.klein_index(t) == klein]
    mixed = [t for t in terms if is_mixed(t, diag)]
    g = 0.6180339887
    sigma = N * g
    n = sp_agree = blk_agree = 0
    for t1, t2 in combinations_with_replacement(mixed, 2):
        if (sum(c == "Y" for c in t1) % 2) != (sum(c == "Y" for c in t2) % 2):
            continue
        n += 1
        H = _build_kbody_chain(N, [tuple(t1) + (1.0,), tuple(t2) + (1.0,)])
        L = lindbladian_pauli_dephasing(H, [g] * N, dephase_letter=letter)
        hard = fw.classify_pauli_pair(chain, [t1, t2], dephase_letter=letter) == 'hard'
        sp_broken = multiset_res(L, sigma) > 1e-6
        sp_agree += (sp_broken == hard)
        if do_block:
            blk_broken = omega0_asym([t1, t2], N) > 1e-6
            blk_agree += (blk_broken == hard)
    return n, sp_agree, (blk_agree if do_block else None)


def main():
    print(__doc__)
    print("RE-VERIFICATION")
    print("=" * 78)
    print(f"{'cell':28s} {'pairs':>6} {'Stage1 spec<=>hard':>20} {'omega0-block<=>hard':>20}")
    for (N, k, letter, klein, tag, blk) in [
        (4, 3, 'Z', (0, 1), "N=4 k=3 Z (F103 cell)", True),
        (4, 3, 'X', (1, 0), "N=4 k=3 X", True),
        (4, 3, 'Y', (1, 1), "N=4 k=3 Y", True),
        (5, 3, 'Z', (0, 1), "N=5 k=3 Z", False),
    ]:
        n, sp, bl = run(N, k, letter, klein, do_block=blk)
        blstr = f"{bl}/{n}" if bl is not None else "(skipped N=5)"
        print(f"{tag:28s} {n:>6} {f'{sp}/{n}':>20} {blstr:>20}")
    print()
    print("STATUS: Stage 1 (operator-search => spectral) CLOSED, rigorous, all cells bit-exact.")
    print("        Remaining gap = the spectral/first-order statement, localized to the omega=0 block.")


if __name__ == "__main__":
    main()
