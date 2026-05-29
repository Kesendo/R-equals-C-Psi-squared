#!/usr/bin/env python3
"""Closed-form rule for the F87 42:8 diagonal-cell hard split (F103/F105/F106/F107/F110).

In the diagonal Klein cell for dephase letter D (the cell whose (bit_a,bit_b) = D's),
call {I, D} the DIAGONAL letters (they commute with D-dephasing). A k=3 cell term has
either 3 diagonal letters (an all-diagonal "pure-D template", n_X=n_Y=0) or exactly 1
(two off-diagonal letters {X,Y} plus one {I,D}); the cell forces n_diagonal in {1,3}.

    Hardness rule (verified bit-exact, 0 mismatches, at N=4 and N=5, for D in {Z,X,Y}):
      a k=3 pair is F87-hard  iff
        (a) at least one term is all-diagonal (pure-D template), OR
        (b) both terms are single-diagonal AND their lone {I,D} letter sits at
            chain-ADJACENT window positions (|Δpos| = 1);
      otherwise soft.

Consequences (the open 42:8 closed form):
  * single-diagonal pairs: 8 adjacent (hard) + 13 non-adjacent (soft), symmetric in y_par.
  * all-diagonal templates (4 of them: D at one site, or DDD) -> 34 hard pairs
    (= 55 - 21 = pairs in that y_par minus the 21 single-diagonal pairs).
  * pure-D templates are Y-free for D in {Z,X} (y_par=0) but ARE Y's for D=Y (y_par=1),
    so the 34 extra hard pairs land in y_par=0 for Z/X and y_par=1 for Y.
  => hard = (34 + 8):(0 + 8) = 42:8 for Z,X and 8:42 (Y-inversion) for Y; soft 13:13.

The two atomic sub-rules (a),(b) are the F87 mechanism at the term level; (a) is the
k=3 face of F111's pure-D template rule. Verified here; a palindrome-level proof of the
atomic rules is the remaining layer.
"""
from __future__ import annotations
import sys
from itertools import product, combinations_with_replacement
from pathlib import Path

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw

DIAG_CELL = {'Z': (0, 1), 'X': (1, 0), 'Y': (1, 1)}


def y_par(t):
    return sum(1 for L in t if L == 'Y') % 2


def check(N):
    chain = fw.ChainSystem(N=N)
    print(f"=== N={N} ===")
    all_ok = True
    for D in 'ZXY':
        cell = DIAG_CELL[D]
        diag = {'I', D}
        terms = [t for t in product('IXYZ', repeat=3)
                 if not all(L == 'I' for L in t) and fw.klein_index(t) == cell]

        def ndiag(t):
            return sum(1 for L in t if L in diag)

        def dpos(t):
            return next(i for i, L in enumerate(t) if L in diag)

        def predict(t1, t2):
            if ndiag(t1) == 3 or ndiag(t2) == 3:
                return 'hard'                                   # (a) all-diagonal template
            return 'hard' if abs(dpos(t1) - dpos(t2)) == 1 else 'soft'  # (b) adjacent

        counts = {(c, y): 0 for c in ('hard', 'soft', 'truly') for y in (0, 1)}
        mismatch = 0
        for t1, t2 in combinations_with_replacement(terms, 2):
            if y_par(t1) != y_par(t2):
                continue
            actual = fw.classify_pauli_pair(chain, [t1, t2], dephase_letter=D)
            counts[(actual, y_par(t1))] += 1
            if actual != predict(t1, t2):
                mismatch += 1
        if mismatch:
            all_ok = False
        h = (counts[('hard', 0)], counts[('hard', 1)])
        s = (counts[('soft', 0)], counts[('soft', 1)])
        print(f"  D={D} cell={cell}  hard={h}  soft={s}  mismatches={mismatch}")
    return all_ok


def main():
    ok = all(check(N) for N in (4, 5))
    print()
    print("RULE CONFIRMED (0 mismatches everywhere)." if ok else "RULE FAILED.")
    print("42:8 = (34 all-diagonal + 8 adjacent):(0 + 8); Y-inversion = templates carry y_par=1.")


if __name__ == "__main__":
    main()
