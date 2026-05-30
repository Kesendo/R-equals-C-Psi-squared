#!/usr/bin/env python3
"""Why no odd cycle at full support? The flip-generator count.

In the Z-basis a Pauli term acts as a single bit-flip mask (its X/Y positions). The hopping
graph G_H is the Cayley graph of the set S of distinct edge-XOR masks {a^b : H[a,b]!=0}. G_H
is bipartite iff there is a LINEAR phi: Z2^N -> Z2 with phi(s)=1 for every s in S (then
color(a)=phi(a) 2-colours it; this phi IS the chiral K = diag((-1)^phi)). Such phi exists iff
S has no odd linear relation.

At FULL SUPPORT k=N each term places once -> one mask per term -> |S| <= 2 -> phi always
exists -> bipartite -> soft. Odd cycles (rule b) need |S| >= 3, which requires k<N (sliding
windows give one mask per (term, window)). N=3 is special only as the smallest k=3 full-support
case. Test: |S| distribution for zero-diagonal (Mixed+Mixed) pairs, full support vs windows.
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

DIAG = {'I', 'Z'}  # Z-deph diagonal letters


def edge_generators(H, tol=1e-9):
    d = H.shape[0]
    S = set()
    for a in range(d):
        row = np.where(np.abs(H[a]) > tol)[0]
        for b in row:
            if b > a:
                S.add(int(a) ^ int(b))
    return S


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


def phi_exists_gf2(S, nbits):
    """Is the GF(2) system { sum_i phi_i s_i = 1 : s in S } consistent? (phi linear, no const)."""
    rows = [[(s >> i) & 1 for i in range(nbits)] + [1] for s in S]
    rows = [r[:] for r in rows]
    piv_col = 0
    r = 0
    R = len(rows)
    for c in range(nbits):
        pivot = next((i for i in range(r, R) if rows[i][c]), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        for i in range(R):
            if i != r and rows[i][c]:
                rows[i] = [(x ^ y) for x, y in zip(rows[i], rows[r])]
        r += 1
    # inconsistent iff a row is all-zero in coeffs but 1 in augmented col
    for row in rows:
        if not any(row[:nbits]) and row[nbits]:
            return False
    return True


def main():
    for (N, k, tag) in [(3, 3, 'FULL k=N=3 (N=3 special)'),
                        (4, 4, 'FULL k=N=4 (F111)'),
                        (4, 3, 'WINDOWS k=3<N=4 (F103)')]:
        chain = fw.ChainSystem(N=N)
        terms = [t for t in product('IXYZ', repeat=k)
                 if not all(L == 'I' for L in t) and fw.klein_index(t) == (0, 1)]
        # Mixed+Mixed = both terms have an off-diagonal letter (zero diagonal)
        def is_mixed(t):
            return any(L not in DIAG for L in t)
        mixed = [t for t in terms if is_mixed(t)]
        smax = 0
        bip_all = True
        soft = hard = 0
        phi_bip_mismatch = 0
        for t1, t2 in combinations_with_replacement(mixed, 2):
            if (sum(c == 'Y' for c in t1) % 2) != (sum(c == 'Y' for c in t2) % 2):
                continue
            H = _build_kbody_chain(N, [tuple(t1) + (1.0,), tuple(t2) + (1.0,)])
            S = edge_generators(H)
            smax = max(smax, len(S))
            bip = is_bipartite(H)
            bip_all = bip_all and bip
            if phi_exists_gf2(S, N) != bip:
                phi_bip_mismatch += 1
            cls = fw.classify_pauli_pair(chain, [t1, t2], dephase_letter='Z')
            if cls == 'soft':
                soft += 1
            elif cls == 'hard':
                hard += 1
        print(f"{tag}: Mixed+Mixed pairs  max|S|={smax}  all-bipartite={bip_all}  "
              f"soft={soft} hard={hard}  (phi==bipartite mismatches={phi_bip_mismatch})")


if __name__ == "__main__":
    main()
