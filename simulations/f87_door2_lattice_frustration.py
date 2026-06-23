#!/usr/bin/env python3
"""Door 2 (a): port the palindrome mask-test to 2D lattices -- is hardness = geometric frustration?
(scout, 2026-06-05).

On the chain the mask-test for bond hopping was: exists phi with phi(e_i + e_j) = 1 for every bond
(i,j). That is phi_i + phi_j = 1 on every edge, i.e. a proper 2-COLOURING of the lattice graph. So the
classifier predicts soft <=> the lattice is BIPARTITE. On a 1D chain that is always true (no odd cycle),
so this face was invisible; in 2D a frustrated lattice (triangle, triangular patch) is non-bipartite and
should be HARD. This scout checks that prediction against the actual spectral verdict on small graphs.
"""
from __future__ import annotations

import sys
from collections import deque
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try: sys.stdout.reconfigure(encoding="utf-8")
    except Exception: pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from framework.pauli import _build_bilinear
from framework.lindblad import lindbladian_z_dephasing, palindrome_residual


def graph_bipartite(N, bonds):
    """2-colour the graph; True iff bipartite (no odd cycle)."""
    color = {}
    adj = {i: [] for i in range(N)}
    for i, j in bonds:
        adj[i].append(j); adj[j].append(i)
    for s in range(N):
        if s in color: continue
        color[s] = 0; q = deque([s])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if v not in color: color[v] = color[u] ^ 1; q.append(v)
                elif color[v] == color[u]: return False
    return True


def spectral_verdict(N, bonds, bilinear, gamma=0.05, op_tol=1e-10, spec_tol=1e-6):
    """truly / soft / hard for H = sum_bonds bilinear, under uniform Z-dephasing."""
    H = _build_bilinear(N, bonds, bilinear)
    L = lindbladian_z_dephasing(H, [gamma] * N)
    sigma = N * gamma
    M = palindrome_residual(L, sigma, N)
    if np.linalg.norm(M) < op_tol:
        return 'truly'
    ev = np.linalg.eigvals(L)
    used = np.zeros(len(ev), dtype=bool); max_err = 0.0
    for i in range(len(ev)):
        if used[i]: continue
        target = -ev[i] - 2 * sigma
        d = np.abs(ev - target); d[used] = np.inf
        j = int(np.argmin(d))
        used[i] = True
        if j != i: used[j] = True
        max_err = max(max_err, float(d[j]))
    return 'soft' if max_err < spec_tol else 'hard'


def mask_bipartite(N, bonds):
    """The flip-mask bipartite test for single-flip-per-bond hopping = the graph 2-colouring."""
    masks = {(1 << i) | (1 << j) for i, j in bonds}
    basis = []
    for m in masks:
        r = (m << 1) | 1
        for b in basis: r = min(r, r ^ b)
        if r == 1: return False
        if r: basis.append(r); basis.sort(reverse=True)
    return True


def main():
    print("=" * 86)
    print("Door 2 (a): palindrome hardness vs lattice frustration (bipartiteness)")
    print("=" * 86)
    graphs = [
        ("chain-3 (bipartite)",    3, [(0, 1), (1, 2)]),
        ("triangle K3 (FRUSTRATED)", 3, [(0, 1), (1, 2), (0, 2)]),
        ("square-4 (bipartite)",   4, [(0, 1), (1, 2), (2, 3), (3, 0)]),
        ("4-path (bipartite)",     4, [(0, 1), (1, 2), (2, 3)]),
        ("5-cycle (FRUSTRATED)",   5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]),
        ("6-cycle (bipartite)",    6, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]),
        ("K4 (FRUSTRATED)",        4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]),
    ]
    hoppings = [
        ("XX",     [('X', 'X', 1.0)]),
        ("XX+YY",  [('X', 'X', 1.0), ('Y', 'Y', 1.0)]),
        ("XY+YX",  [('X', 'Y', 1.0), ('Y', 'X', 1.0)]),
    ]
    for hlabel, bilinear in hoppings:
        print(f"\n  hopping = {hlabel}:")
        print(f"    {'graph':<26} {'bipartite':<11} {'mask-test':<11} {'spectral':<9} hard==frustrated?")
        for glabel, N, bonds in graphs:
            bip = graph_bipartite(N, bonds)
            mt = mask_bipartite(N, bonds)
            sp = spectral_verdict(N, bonds, bilinear)
            # prediction: hard <=> non-bipartite (frustrated)
            consistent = (sp == 'hard') == (not bip)
            print(f"    {glabel:<26} {('yes' if bip else 'NO'):<11} "
                  f"{('soft' if mt else 'HARD'):<11} {sp:<9} {'OK' if consistent else 'MISMATCH'}")
    print("\n  reading: if hard appears EXACTLY on the non-bipartite (frustrated) graphs, then")
    print("  palindrome hardness = geometric frustration -- a face invisible on the always-bipartite")
    print("  1D chain, surfacing only when we leave it. (truly = the M=0 bucket, counts as not-hard.)")


if __name__ == "__main__":
    main()
