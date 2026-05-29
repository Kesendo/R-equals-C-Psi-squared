#!/usr/bin/env python3
"""Does the F103 §7 bipartite criterion extend to k=4 (the F111 / F106 regime)?

k=4 body terms on N=4 = full-support (k=N), one placement, no sliding window. F111: a k=4
diagonal-cell pair is hard iff >=1 term is a pure-D template; the blocked direction was
Mixed+Mixed = soft. If the bipartite criterion holds at k=4 too (soft <=> hopping graph
bipartite in the dephase basis), it would EXPLAIN F111: a template lifts the diagonal -> no
chiral K -> hard, and Mixed+Mixed keeps the diagonal zero -> (if bipartite) soft. Test it.
"""
from __future__ import annotations
import sys
from collections import deque
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

KLEIN = {'Z': (0, 1), 'X': (1, 0), 'Y': (1, 1)}
DIAGSET = {'Z': {'I', 'Z'}, 'X': {'I', 'X'}, 'Y': {'I', 'Y'}}
_R2 = 1.0 / np.sqrt(2.0)
ROT = {
    'Z': np.eye(2, dtype=complex),
    'X': _R2 * np.array([[1, 1], [1, -1]], dtype=complex),
    'Y': _R2 * np.array([[1, -1j], [1, 1j]], dtype=complex),
}


def rotate(H, N, dl):
    U = ROT[dl]
    for _ in range(N - 1):
        U = np.kron(U, ROT[dl])
    return U @ H @ U.conj().T


def is_bipartite(H, tol=1e-9):
    d = H.shape[0]
    if np.max(np.abs(np.diag(H))) > tol:
        return False
    color = {}
    for s in range(d):
        if s in color:
            continue
        color[s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            for v in np.where(np.abs(H[u]) > tol)[0]:
                if v == u:
                    continue
                if v not in color:
                    color[v] = color[u] ^ 1
                    q.append(v)
                elif color[v] == color[u]:
                    return False
    return True


def y_par(t):
    return sum(1 for L in t if L == 'Y') % 2


def is_pure_d_template(t, diag):
    return all(L in diag for L in t)


def main():
    N = 4
    k = 4
    chain = fw.ChainSystem(N=N)
    for dl in ('Z', 'X', 'Y'):
        cell = KLEIN[dl]
        diag = DIAGSET[dl]
        terms = [t for t in product('IXYZ', repeat=k)
                 if not all(L == 'I' for L in t) and fw.klein_index(t) == cell]
        mism = 0
        nhard = nsoft = ntruly = 0
        # cross-check: among non-bipartite hard pairs, how many are template-driven vs odd-cycle
        oddcycle_hard = 0
        for t1, t2 in combinations_with_replacement(terms, 2):
            if y_par(t1) != y_par(t2):
                continue
            H = _build_kbody_chain(N, [tuple(t1) + (1.0,), tuple(t2) + (1.0,)])
            bip = is_bipartite(rotate(H, N, dl))
            actual = fw.classify_pauli_pair(chain, [t1, t2], dephase_letter=dl)
            if actual == 'hard':
                nhard += 1
                has_template = is_pure_d_template(t1, diag) or is_pure_d_template(t2, diag)
                if not has_template:
                    oddcycle_hard += 1
            elif actual == 'soft':
                nsoft += 1
            else:
                ntruly += 1
            if actual in ('soft', 'hard'):
                pred = 'soft' if bip else 'hard'
                if actual != pred:
                    mism += 1
                    if mism <= 6:
                        print(f"  [{dl}] MISMATCH {''.join(t1)}+{''.join(t2)}: actual={actual} bipartite={bip}")
        print(f"dephase {dl}: hard={nhard} soft={nsoft} truly={ntruly}  "
              f"mismatches(soft<=>bipartite)={mism}  hard-without-template={oddcycle_hard}")


if __name__ == "__main__":
    main()
