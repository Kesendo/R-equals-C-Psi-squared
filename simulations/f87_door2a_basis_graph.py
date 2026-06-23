#!/usr/bin/env python3
"""Door 2a resolution: the soft criterion is the BASIS-STATE graph bipartiteness, not the site graph
(scout, 2026-06-05).

The §7.5 -N-mode criterion soft <=> exists diagonal D with {H,D}=0. For an off-diagonal H and diagonal
D, {H,D}_ij = H_ij (D_i + D_j), so the condition is D_i = -D_j on every edge of the BASIS-STATE hopping
graph (2^N nodes, edges = nonzero off-diagonals of H). That is a proper signed 2-colouring, i.e. the
basis-state graph being BIPARTITE. This is the actual (letter-dependent) criterion; the mask/site
bipartite test (the chiral K on the LATTICE graph) is a chain-specific proxy that coincides with it on
the chain but not in general.

Hypothesis: XY+YX (pure pair term, Δn=±2) has a BIPARTITE basis-state graph (colour by n mod 4) even on
a frustrated lattice, so it is soft; XY / XY-YX (with Δn=0 hopping) do not. This scout builds the
basis-state graph for each and checks basis-bipartite <=> not-hard against the spectral verdict.
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


def spectral_verdict(N, bonds, bilinear, gamma=0.05, op_tol=1e-10, spec_tol=1e-6):
    H = _build_bilinear(N, bonds, bilinear)
    L = lindbladian_z_dephasing(H, [gamma] * N)
    sigma = N * gamma
    if np.linalg.norm(palindrome_residual(L, sigma, N)) < op_tol:
        return 'truly'
    ev = np.linalg.eigvals(L)
    used = np.zeros(len(ev), dtype=bool); mx = 0.0
    for i in range(len(ev)):
        if used[i]: continue
        d = np.abs(ev - (-ev[i] - 2 * sigma)); d[used] = np.inf
        j = int(np.argmin(d)); used[i] = True
        if j != i: used[j] = True
        mx = max(mx, float(d[j]))
    return 'soft' if mx < spec_tol else 'hard'


def basis_graph_bipartite(N, bonds, bilinear):
    """Build H, take its basis-state hopping graph (nonzero off-diagonals), test bipartiteness.
    Returns (is_bipartite, colouring_or_None)."""
    H = _build_bilinear(N, bonds, bilinear)
    d = 1 << N
    adj = {i: [] for i in range(d)}
    A = np.abs(H) > 1e-12
    for i in range(d):
        for j in range(i + 1, d):
            if A[i, j]:
                adj[i].append(j); adj[j].append(i)
    color = {}
    for s in range(d):
        if s in color: continue
        color[s] = 0; q = deque([s])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if v not in color: color[v] = color[u] ^ 1; q.append(v)
                elif color[v] == color[u]: return False, None
    return True, color


def main():
    print("=" * 84)
    print("Door 2a resolution: soft <=> BASIS-STATE graph bipartite (not the site graph)")
    print("=" * 84)
    triangle = (3, [(0, 1), (1, 2), (0, 2)])
    chain3 = (3, [(0, 1), (1, 2)])
    cases = [
        ("XY+YX  (pure pair)",   [('X', 'Y', 1.0), ('Y', 'X', 1.0)]),
        ("XY     (pair+hop)",    [('X', 'Y', 1.0)]),
        ("XY-YX  (pure hop)",    [('X', 'Y', 1.0), ('Y', 'X', -1.0)]),
        ("XX+YY  (sym hop)",     [('X', 'X', 1.0), ('Y', 'Y', 1.0)]),
    ]
    for glabel, (N, bonds) in [("triangle K3", triangle), ("chain-3", chain3)]:
        print(f"\n  {glabel}:")
        print(f"    {'H':<22} {'basis-bipartite':<16} {'spectral':<9} consistent?")
        for clabel, bilinear in cases:
            bip, col = basis_graph_bipartite(N, bonds, bilinear)
            sp = spectral_verdict(N, bonds, bilinear)
            ok = bip == (sp != 'hard')
            print(f"    {clabel:<22} {('yes' if bip else 'NO'):<16} {sp:<9} {'OK' if ok else 'MISMATCH'}")

    # confirm the XY+YX colouring is the excitation-number (n mod 4) structure
    print("\n  XY+YX basis-state 2-colouring vs n mod 4 (triangle):")
    N, bonds = triangle
    bip, col = basis_graph_bipartite(N, bonds, [('X', 'Y', 1.0), ('Y', 'X', 1.0)])
    if bip:
        agree = all(col[s] == ((bin(s).count('1') // 2) % 2) for s in col)
        print(f"    colour[s] == (popcount(s)//2) mod 2 for all s: {agree}  "
              f"(pairing changes n by +-2, so n mod 4 two-colours it)")

    print("\n  reading: if basis-bipartite == not-hard everywhere, the soft criterion is the BASIS-STATE")
    print("  graph 2-colouring (letter-dependent), and the mask/site test is its chain proxy. XY+YX is")
    print("  soft because its pairing keeps the basis graph bipartite (n mod 4) even when the site is")
    print("  frustrated -- the 'other mechanism' is the excitation-number structure, not a new symmetry.")


if __name__ == "__main__":
    main()
