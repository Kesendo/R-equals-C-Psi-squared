#!/usr/bin/env python3
"""Full diagonal-cell test: F87 soft <=> H's hopping graph bipartite, for ALL cell pairs.

Covers templates (pure-D, all-diagonal), single-diagonal, and mixed terms, both y-parities,
all three dephase letters. If soft<=>bipartite holds cell-wide with zero mismatches, the
bipartiteness criterion IS the F87 diagonal-cell rule, and the atomic rules (a) template->hard
and (b) adjacency->hard are its two combinatorial faces:
  (a) a pure-D template is diagonal -> nonzero H-diagonal -> not bipartite (no K with KHK=-H)
  (b) opposite position-parity -> odd hopping cycle -> not bipartite
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

# Single-qubit rotation U_P with U_P * P * U_P^dag = Z (maps the dephase letter to the
# computational axis, so the chiral K becomes diagonal in the comp basis there).
_R2 = 1.0 / np.sqrt(2.0)
ROT = {
    'Z': np.eye(2, dtype=complex),
    'X': _R2 * np.array([[1, 1], [1, -1]], dtype=complex),          # Hadamard: H X H = Z
    'Y': _R2 * np.array([[1, -1j], [1, 1j]], dtype=complex),        # U Y U^dag = Z
}


def rotate_to_dephase_basis(H, N, dl):
    U1 = ROT[dl]
    U = U1
    for _ in range(N - 1):
        U = np.kron(U, U1)
    return U @ H @ U.conj().T


def is_bipartite(H, tol=1e-9):
    d = H.shape[0]
    if np.max(np.abs(np.diag(H))) > tol:
        return False
    adj = [np.where(np.abs(H[a]) > tol)[0] for a in range(d)]
    color = {}
    for start in range(d):
        if start in color:
            continue
        color[start] = 0
        q = deque([start])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if v != u and v not in color:
                    color[v] = color[u] ^ 1
                    q.append(v)
                elif v != u and color[v] == color[u]:
                    return False
    return True


def y_par(t):
    return sum(1 for L in t if L == 'Y') % 2


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    letters = tuple(sys.argv[2]) if len(sys.argv) > 2 else ('Z', 'X', 'Y')
    chain = fw.ChainSystem(N=N)
    print(f"--- N={N} ---")
    for dl in letters:
        cell = KLEIN[dl]
        diag = DIAGSET[dl]
        terms = [t for t in product('IXYZ', repeat=3)
                 if not all(L == 'I' for L in t) and fw.klein_index(t) == cell]
        mism = 0
        nhard = nsoft = 0
        for t1, t2 in combinations_with_replacement(terms, 2):
            if y_par(t1) != y_par(t2):
                continue
            H = _build_kbody_chain(N, [tuple(t1) + (1.0,), tuple(t2) + (1.0,)])
            bip = is_bipartite(rotate_to_dephase_basis(H, N, dl))
            actual = fw.classify_pauli_pair(chain, [t1, t2], dephase_letter=dl)
            if actual == 'hard':
                nhard += 1
            elif actual == 'soft':
                nsoft += 1
            pred = 'soft' if bip else 'hard'
            # 'truly' pairs (||M||~0) are neither soft nor hard; bipartite trivially holds
            if actual in ('soft', 'hard') and actual != pred:
                mism += 1
                if mism <= 6:
                    print(f"  [{dl}] MISMATCH {''.join(t1)}+{''.join(t2)}: actual={actual} bipartite={bip}")
        print(f"dephase {dl}: hard={nhard} soft={nsoft}  mismatches(soft<=>bipartite) = {mism}")


if __name__ == "__main__":
    main()
