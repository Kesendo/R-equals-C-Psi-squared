#!/usr/bin/env python3
"""Door 2 deeper: the family of scalable basis-graph 2-colourings (scout, 2026-06-05).

§7.12: among non-truly H, soft <=> the basis-state hopping graph is bipartite. That is the 2^N truth;
a SCALABLE soft-certifier needs structured 2-colourings c(s) checkable per-term. Two are in the
certifier: LINEAR c(s)=<phi,s> (the chiral K) and EXCITATION-2 c(s)=floor(n/2) mod 2 (pure Δn=±2
pairing). This scout maps the whole family on small N: for each soft H it builds the actual basis-graph
and tests which structured colourings 2-colour it --
    LIN     : some linear <phi,s>            (handles bipartite flip-value graphs / the chiral K)
    PAR     : n mod 2                         (handles all-edges-odd-Δn, i.e. every term odd k_xy)
    PAIR    : floor(n/2) mod 2                (handles pure Δn=±2 pairings)
    LIN+PAR : <phi,s> XOR (n mod 2)
    LIN+PAIR: <phi,s> XOR (floor(n/2) mod 2)
-- and flags any soft case where the basis-graph is bipartite but NO structured colouring works (the
residual that needs the full 2^N graph). Also confirms no colouring 2-colours a HARD case (safety).
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

def popcount(x): return bin(x).count("1")

def spectral_verdict(N, bonds, bilinear, gamma=0.05, op_tol=1e-10, spec_tol=1e-6):
    H = _build_bilinear(N, bonds, bilinear)
    L = lindbladian_z_dephasing(H, [gamma] * N); sigma = N * gamma
    if np.linalg.norm(palindrome_residual(L, sigma, N)) < op_tol: return 'truly'
    ev = np.linalg.eigvals(L); used = np.zeros(len(ev), bool); mx = 0.0
    for i in range(len(ev)):
        if used[i]: continue
        d = np.abs(ev - (-ev[i] - 2 * sigma)); d[used] = np.inf
        j = int(np.argmin(d)); used[i] = True
        if j != i: used[j] = True
        mx = max(mx, float(d[j]))
    return 'soft' if mx < spec_tol else 'hard'

def basis_edges(N, bonds, bilinear):
    H = _build_bilinear(N, bonds, bilinear); d = 1 << N
    A = np.abs(H) > 1e-12
    return [(i, j) for i in range(d) for j in range(i + 1, d) if A[i, j]]

def is_bipartite(N, edges):
    adj = {i: [] for i in range(1 << N)}
    for i, j in edges: adj[i].append(j); adj[j].append(i)
    col = {}
    for s in range(1 << N):
        if s in col: continue
        col[s] = 0; q = deque([s])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if v not in col: col[v] = col[u] ^ 1; q.append(v)
                elif col[v] == col[u]: return False
    return True

def colours_with(N, edges, c):
    return all(c(i) != c(j) for i, j in edges)

def which_structured(N, edges):
    """Return the set of structured colourings that 2-colour these edges."""
    out = []
    lin = [phi for phi in range(1 << N) if colours_with(N, edges, lambda s, p=phi: popcount(p & s) & 1)]
    if lin: out.append("LIN")
    if colours_with(N, edges, lambda s: popcount(s) & 1): out.append("PAR")
    if colours_with(N, edges, lambda s: (popcount(s) // 2) & 1): out.append("PAIR")
    if any(colours_with(N, edges, lambda s, p=phi: (popcount(p & s) & 1) ^ (popcount(s) & 1)) for phi in range(1 << N)):
        out.append("LIN+PAR")
    if any(colours_with(N, edges, lambda s, p=phi: (popcount(p & s) & 1) ^ ((popcount(s) // 2) & 1)) for phi in range(1 << N)):
        out.append("LIN+PAIR")
    return out


def main():
    print("=" * 90)
    print("Door 2: which structured basis-graph 2-colourings cover which soft Hamiltonians?")
    print("=" * 90)
    N = 4
    chain = [(0, 1), (1, 2), (2, 3)]
    triangle = [(0, 1), (1, 2), (0, 2)]
    battery = [
        ("chain", chain, "XY+YX  pair",   [('X','Y',1.), ('Y','X',1.)]),
        ("chain", chain, "XZ+ZX  odd-flip",[('X','Z',1.), ('Z','X',1.)]),
        ("chain", chain, "XY-YX  hop",     [('X','Y',1.), ('Y','X',-1.)]),
        ("chain", chain, "XY     pair+hop",[('X','Y',1.)]),
        ("chain", chain, "XZ     1-flip",  [('X','Z',1.)]),
        ("chain", chain, "all-odd 4-term", [('X','Z',1.), ('Z','X',1.), ('Y','Z',1.), ('Z','Y',1.)]),
        ("chain", chain, "pair+odd mix",   [('X','Y',1.), ('Y','X',1.), ('X','Z',1.), ('Z','X',1.)]),
        ("tri",   triangle, "XY+YX pair",  [('X','Y',1.), ('Y','X',1.)]),
        ("tri",   triangle, "XZ+ZX odd",   [('X','Z',1.), ('Z','X',1.)]),
        ("tri",   triangle, "XY-YX hop",   [('X','Y',1.), ('Y','X',-1.)]),
        ("tri",   triangle, "XY    hard",  [('X','Y',1.)]),
        ("tri",   triangle, "all-odd 4-term",[('X','Z',1.), ('Z','X',1.), ('Y','Z',1.), ('Z','Y',1.)]),
        ("tri",   triangle, "pair+odd mix", [('X','Y',1.), ('Y','X',1.), ('X','Z',1.), ('Z','X',1.)]),
    ]
    print(f"\n  {'graph':<5} {'H':<16} {'spectral':<8} {'bipartite':<10} {'structured colourings'}")
    print("  " + "-" * 80)
    safety_ok = True
    for glabel, bonds, hlabel, bil in battery:
        sp = spectral_verdict(N, bonds, bil)
        edges = basis_edges(N, bonds, bil)
        bip = is_bipartite(N, edges)
        cols = which_structured(N, edges) if bip else []
        if sp == 'hard' and cols:  # a structured colouring on a hard graph = false certificate
            safety_ok = False
        residual = (sp == 'soft' and bip and not cols)
        tag = "  <- RESIDUAL (bipartite, no structured colouring)" if residual else ""
        print(f"  {glabel:<5} {hlabel:<16} {sp:<8} {('yes' if bip else 'NO'):<10} {','.join(cols) if cols else '-'}{tag}")

    print(f"\n  safety (no structured colouring on a hard graph): {safety_ok}")
    print("  reading: PAR (n mod 2) should cover the all-odd-k_xy cases (XZ+ZX, XZ) even on the frustrated")
    print("  triangle, where LIN fails -- the third scalable strategy. Any RESIDUAL row is a soft case no")
    print("  structured colouring reaches (would need the 2^N graph); that bounds the certifier's reach.")


if __name__ == "__main__":
    main()
