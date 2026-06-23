#!/usr/bin/env python3
"""EQ-031 step (b): test ‖M(N)‖² scaling for non-chain topologies.

The chain absolute formulas are:
  main class         ‖M(N)‖² = c_H · (N − 1) · 4^(N−2)
  single-body class  ‖M(N)‖² = c_H · (2N − 3) · 4^(N−2)

Conjecture for general topology G (graph with N sites):
  main class         ‖M(N)‖² = c_H · B(G) · 4^(N−2)
  single-body class  ‖M(N)‖² = c_H · D2(G)/2 · 4^(N−2)

where B(G) = bond count and D2(G) = Σ_i deg_G(i)². The chain case
recovers B = N − 1 and D2/2 = (4N − 6)/2 = 2N − 3 from the path graph
P_N's degree sequence (1, 2, 2, ..., 2, 1).

Anchor c_H is topology-independent: at the smallest setting (N=3) chain
and 3-site path-star are graph-isomorphic, so c_H matches between them.
For ring and complete graph at N=3, the same per-bond-bilinear definition
applies.

This script tests at N=4 and N=5:
  Chain:  bonds (0,1), (1,2), ..., (N-2, N-1)            B = N-1
  Ring:   chain + closing bond (N-1, 0)                  B = N
  Star:   bonds (0, 1), (0, 2), ..., (0, N-1)            B = N-1, hub deg = N-1
  K_N:    all (i,j), i<j                                  B = N(N-1)/2
"""
import math
import sys
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


GAMMA = 0.1
J = 1.0


def topology_bonds(name, N):
    if name == 'chain':
        return [(i, i + 1) for i in range(N - 1)]
    if name == 'ring':
        return [(i, (i + 1) % N) for i in range(N)]
    if name == 'star':
        return [(0, i) for i in range(1, N)]
    if name == 'complete':
        return [(i, j) for i in range(N) for j in range(i + 1, N)]
    raise ValueError(f"unknown topology: {name}")


def graph_invariants(bonds, N):
    """Return (B, D2) where B is bond count and D2 = Σ deg²."""
    deg = [0] * N
    for i, j in bonds:
        deg[i] += 1
        deg[j] += 1
    B = len(bonds)
    D2 = sum(d * d for d in deg)
    return B, D2


def m_norm_sq(N, bonds, terms):
    bilinear = [(t[0], t[1], J) for t in terms]
    H = fw._build_bilinear(N, bonds, bilinear)
    L = fw.lindbladian_z_dephasing(H, [GAMMA] * N)
    M = fw.palindrome_residual(L, N * GAMMA, N)
    return float(np.linalg.norm(M)) ** 2


# Test Hamiltonians: one main-class, one single-body
TEST_HAMILTONIANS = [
    ('XY+YZ (main)',  [('X', 'Y'), ('Y', 'Z')], 'main'),
    ('IY+YI (single)', [('I', 'Y'), ('Y', 'I')], 'single_body'),
]


def main():
    print("EQ-031 step (b): topology dependence of ‖M(N)‖² scaling")
    print(f"  γ={GAMMA}, J={J}")
    print()
    print("Conjecture:")
    print(f"  main class         ‖M(N)‖² = c_H · B(G) · 4^(N−2)")
    print(f"  single-body class  ‖M(N)‖² = c_H · D2(G)/2 · 4^(N−2)")
    print()

    topologies = ['chain', 'ring', 'star', 'complete']
    Ns = [4, 5]

    # Anchor c_H from N=3 chain (simplest) for each Hamiltonian.
    # At N=3 chain: B=2, D2=6 → main predicts c_H·2·4=8c_H, single c_H·3·4=12c_H.
    print("Anchoring c_H from N=3 chain (B=2, D2=6, factor 4):")
    chain3 = topology_bonds('chain', 3)
    c_H_table = {}
    for label, terms, cls in TEST_HAMILTONIANS:
        m = m_norm_sq(3, chain3, terms)
        if cls == 'main':
            c_H = m / (2 * 4)  # divide by B · 4^(N-2)
        else:
            c_H = m / (3 * 4)  # D2/2 · 4^(N-2) = 3 · 4
        c_H_table[label] = c_H
        print(f"  {label}: ‖M(3, chain)‖² = {m:.4e}, c_H = {c_H:.4e}")
    print()

    print(f"  {'topology':>10s} {'N':>2s}  "
          f"{'B':>3s} {'D2':>4s}  "
          + "  ".join(f"{label:>32s}" for label, _, _ in TEST_HAMILTONIANS))
    print(f"  {'-' * 10} {'-' * 2}  {'-' * 3} {'-' * 4}  "
          + "  ".join("-" * 32 for _ in TEST_HAMILTONIANS))

    for topo in topologies:
        for N in Ns:
            bonds = topology_bonds(topo, N)
            B, D2 = graph_invariants(bonds, N)
            cells = []
            for label, terms, cls in TEST_HAMILTONIANS:
                m = m_norm_sq(N, bonds, terms)
                c_H = c_H_table[label]
                if cls == 'main':
                    pred = c_H * B * (4 ** (N - 2))
                else:
                    pred = c_H * (D2 / 2) * (4 ** (N - 2))
                ratio = m / pred if pred > 0 else float('nan')
                cells.append(f"  m={m:>10.3e}  m/p={ratio:.6f}")
            print(f"  {topo:>10s} {N:>2d}  {B:>3d} {D2:>4d}  " + "  ".join(cells))

    print()
    print("If meas/pred = 1.000000 across topologies, the conjecture holds.")
    print("If not, the deviation pattern tells us what additional topology")
    print("structure (cycles, hubs, triangles) the formula has missed.")


if __name__ == "__main__":
    main()
