#!/usr/bin/env python3
"""Door 2: is the 'residual' soft case really uncatchable by a scalable colouring? (scout, 2026-06-05).

Before accepting the ceiling (some soft H reachable by no structured colouring, needs the 2^N graph), make
SURE. For the residual witness XY+YX+XZ+ZX on a frustrated graph: (1) re-verify it is soft and its
basis-state graph bipartite (the Z-sign cancellations make the actual edge set sparser than the naive
flip-mask set, so this needs checking); (2) extract the 2-colouring c(s) and compute its GF(2) algebraic
normal form (ANF) DEGREE -- degree 1 = linear (already tried), 2 = quadratic, 3 = cubic, ... A low,
N-stable degree would mean a richer scalable strategy (a degree-d polynomial colouring) CAN catch it and
the ceiling lifts; a high or N-growing degree means it is genuinely non-scalable.
"""
from __future__ import annotations

import sys
from collections import deque
from itertools import combinations
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try: sys.stdout.reconfigure(encoding="utf-8")
    except Exception: pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from framework.pauli import _build_bilinear
from framework.lindblad import lindbladian_z_dephasing, palindrome_residual

def popcount(x): return bin(x).count("1")

def spectral_verdict(N, bonds, bil, gamma=0.05, op_tol=1e-10, spec_tol=1e-6):
    H = _build_bilinear(N, bonds, bil); L = lindbladian_z_dephasing(H, [gamma]*N); sg = N*gamma
    if np.linalg.norm(palindrome_residual(L, sg, N)) < op_tol: return 'truly'
    ev = np.linalg.eigvals(L); used = np.zeros(len(ev), bool); mx = 0.
    for i in range(len(ev)):
        if used[i]: continue
        d = np.abs(ev-(-ev[i]-2*sg)); d[used]=np.inf; j=int(np.argmin(d)); used[i]=True
        if j!=i: used[j]=True
        mx=max(mx, float(d[j]))
    return 'soft' if mx<spec_tol else 'hard'

def basis_adj(N, bonds, bil):
    H = _build_bilinear(N, bonds, bil); d = 1<<N
    A = np.abs(H) > 1e-12
    adj = {i: [j for j in range(d) if j!=i and A[i,j]] for i in range(d)}
    return adj, d

def colour(adj, d):
    """BFS 2-colour; return (bipartite, colour dict, #components)."""
    col = {}; comps = 0
    for s in range(d):
        if s in col: continue
        comps += 1; col[s]=0; q=deque([s])
        while q:
            u=q.popleft()
            for v in adj[u]:
                if v not in col: col[v]=col[u]^1; q.append(v)
                elif col[v]==col[u]: return False, None, comps
    return True, col, comps

def anf_degree(N, col):
    """GF(2) ANF degree of the colouring c (defined on 2^N states)."""
    d = 1<<N
    c = [col[s] for s in range(d)]
    # Mobius (zeta) transform over the subset lattice: a_T = XOR_{u subset of T} c[u]
    a = c[:]
    for i in range(N):
        for s in range(d):
            if s & (1<<i):
                a[s] ^= a[s ^ (1<<i)]
    deg = max((popcount(s) for s in range(d) if a[s]), default=0)
    nterms = sum(a)
    return deg, nterms


def main():
    print("="*84)
    print("Door 2: is the residual catchable? re-verify + ANF degree of its 2-colouring")
    print("="*84)
    cases = [
        ("triangle K3",      3, [(0,1),(1,2),(0,2)]),
        ("K4",               4, [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]),
        ("5-cycle",          5, [(0,1),(1,2),(2,3),(3,4),(4,0)]),
    ]
    bil = [('X','Y',1.), ('Y','X',1.), ('X','Z',1.), ('Z','X',1.)]   # XY+YX + XZ+ZX
    for glabel, N, bonds in cases:
        sp = spectral_verdict(N, bonds, bil)
        adj, d = basis_adj(N, bonds, bil)
        nedges = sum(len(v) for v in adj.values())//2
        bip, col, comps = colour(adj, d)
        print(f"\n  {glabel} (N={N}):  spectral={sp}  basis-edges={nedges}  components={comps}  bipartite={bip}")
        if bip and sp == 'soft':
            deg, nterms = anf_degree(N, col)
            print(f"    2-colouring ANF: degree={deg}  (#monomials={nterms})  "
                  f"[1=linear, 2=quadratic, 3=cubic, ...; N={N} so max possible deg={N}]")
            if comps > 1:
                print(f"    NOTE: {comps} components -> the colouring has free per-component swaps; the ANF")
                print(f"          degree shown is for one BFS choice, not necessarily the minimal one.")

    print("\n  reading: a low, N-STABLE ANF degree (e.g. always 2) means a degree-d polynomial colouring")
    print("  strategy could catch these, lifting the ceiling. A degree that GROWS with N (toward N) means")
    print("  no fixed-degree scalable colouring works -- the ceiling is real. (The #components matters: a")
    print("  disconnected basis-graph has colouring freedom that a single global polynomial need not have.)")


if __name__ == "__main__":
    main()
