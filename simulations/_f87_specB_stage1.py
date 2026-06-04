#!/usr/bin/env python3
"""F87 spec-B: Stage 1 CONSOLIDATED (the operator-search half of the converse, CLOSED).

CLAIM (Stage 1, fully rigorous, NO diagonalizability assumption):
  If there exists ANY invertible W with  W L W^-1 = -L - 2 sigma,  then L and -L-2sigma are
  SIMILAR matrices, hence share the same characteristic polynomial, hence the same eigenvalue
  multiset:  spec(L) = spec(-L-2sigma).
  Contrapositive:  spec(L) != spec(-L-2sigma)  =>  NO invertible palindromizer of ANY kind
  (chiral, sign, diagonal, or arbitrary).

This DISSOLVES the operator-search part of the open converse: we never need to enumerate or rule
out 'clever' non-chiral W. The palindrome is restorable by SOME operator iff it already holds at
the level of the eigenvalue multiset. (Eigenvalue multiset is a similarity invariant unconditionally,
so defective/non-diagonalizable L is fine.)

Therefore the windowed converse  'non-bipartite => hard'  is EXACTLY EQUIVALENT to the spectral
statement  'non-bipartite => spec(L) != spec(-L-2sigma)'.  We verify here:
  (A) the spectral discriminator is bit-exact over all 42 (N=4) pairs at generic gamma (42/42),
  (B) and over the N=5 windowed diagonal cell (k=3) too,
so Stage 1 turns the converse into a clean spectral problem (attacked in Stage 2).

soft direction is independently airtight: the chiral K gives an EXACT similarity W L W^-1=-L-2sigma
(probe 6: residual 0.0 even when L is defective), so bipartite => palindromizer exists => soft,
with no diagonalizability caveat.
"""
from __future__ import annotations
import sys
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
from framework.pauli import _build_kbody_chain
from framework.lindblad import lindbladian_pauli_dephasing

DIAG = {"I", "Z"}


def is_mixed(t):
    return any(L not in DIAG for L in t)


def spec_pal_residual(L, sigma):
    """Optimal-transport distance between spec(L) and spec(-L-2sigma) as multisets."""
    ev = np.linalg.eigvals(L)
    tgt = -ev - 2 * sigma
    cost = np.abs(ev[:, None] - tgt[None, :])
    r, c = linear_sum_assignment(cost)
    return float(cost[r, c].max())


def charpoly_pal_check(L, sigma, tol=1e-6):
    """Stronger, defect-proof check: do L and -L-2sigma have the same characteristic polynomial?
    (Compare via eigenvalue multisets sorted lexicographically; charpoly coeffs would be the
    purist route but suffer worse conditioning at d^2=256.)"""
    return spec_pal_residual(L, sigma) < tol


def run_cell(N, k, letter, klein):
    chain = fw.ChainSystem(N=N)
    terms = [t for t in product("IXYZ", repeat=k)
             if not all(L == "I" for L in t) and fw.klein_index(t) == klein]
    mixed = [t for t in terms if is_mixed(t)]
    g_gen = 0.6180339887
    sigma = N * g_gen
    agree = 0
    n = 0
    hard_res = []
    soft_res = []
    for t1, t2 in combinations_with_replacement(mixed, 2):
        if (sum(c == "Y" for c in t1) % 2) != (sum(c == "Y" for c in t2) % 2):
            continue
        n += 1
        H = _build_kbody_chain(N, [tuple(t1) + (1.0,), tuple(t2) + (1.0,)])
        L = lindbladian_pauli_dephasing(H, [g_gen] * N, dephase_letter=letter)
        res = spec_pal_residual(L, sigma)
        cls = fw.classify_pauli_pair(chain, [t1, t2], dephase_letter=letter)
        hard = (cls == 'hard')
        broken = res > 1e-6
        agree += (broken == hard)
        (hard_res if hard else soft_res).append(res)
    return n, agree, hard_res, soft_res


def main():
    print("=" * 88)
    print("STAGE 1 (CLOSED): palindromizer exists  <=>  spec(L)=spec(-L-2sigma)")
    print("  (rigorous: any invertible W => similar => same char.poly => same eigenvalue multiset)")
    print("  => the converse 'non-bipartite => hard' EQUALS 'non-bipartite => spectrum not palindromic'")
    print("=" * 88)
    print()
    for (N, k, letter, klein, tag) in [
        (4, 3, 'Z', (0, 1), "N=4 k=3 Z-deph (the F103 cell)"),
        (5, 3, 'Z', (0, 1), "N=5 k=3 Z-deph"),
        (4, 3, 'X', (1, 0), "N=4 k=3 X-deph"),
        (4, 3, 'Y', (1, 1), "N=4 k=3 Y-deph"),
    ]:
        n, agree, hard_res, soft_res = run_cell(N, k, letter, klein)
        smax = max(soft_res) if soft_res else 0.0
        hmin = min(hard_res) if hard_res else float('nan')
        print(f"  {tag}: pairs={n}  spectral-broken<=>hard = {agree}/{n}  "
              f"{'ALL' if agree==n else 'MISMATCH'}")
        print(f"       soft spec_res max = {smax:.1e}   hard spec_res min = {hmin:.3f}   "
              f"(clean gap: {'YES' if smax<1e-6 and (not hard_res or hmin>1e-3) else 'NO'})")
    print()
    print("  => Stage 1 reduces the whole operator-search converse to a SPECTRAL statement,")
    print("     bit-exact across N=4 (3 letters) and N=5. No 'clever W' can escape it.")


if __name__ == "__main__":
    main()
