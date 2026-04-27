#!/usr/bin/env python3
"""Demo of recommend_initial_state() — best ρ_0 for given H.

Tests the recommender against several Hamiltonians at N=3 with γ_T1 = 0.1·γ:
each shows the top-5 initial states ranked by Π-protected count.
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


def main():
    N = 3
    GAMMA, GAMMA_T1, J = 0.1, 0.01, 1.0
    bonds = [(i, i + 1) for i in range(N - 1)]

    cases = [
        ('truly XX+YY (Heisenberg-XY)', [('X', 'X', J), ('Y', 'Y', J)]),
        ('XY+YX (bond-flipped Z-free)', [('X', 'Y', J), ('Y', 'X', J)]),
        ('IY+YI (factorising)',         [('I', 'Y', J), ('Y', 'I', J)]),
        ('YZ+ZY (bond-flipped Z-cont)', [('Y', 'Z', J), ('Z', 'Y', J)]),
        ('XX+ZZ (Heisenberg-XZ)',       [('X', 'X', J), ('Z', 'Z', J)]),
        ('XX+YY+ZZ (full Heisenberg)',  [('X', 'X', J), ('Y', 'Y', J),
                                          ('Z', 'Z', J)]),
        ('XX+YZ (no-symmetry)',         [('X', 'X', J), ('Y', 'Z', J)]),
    ]

    print(f"recommend_initial_state demo, N={N}")
    print(f"  γ_deph={GAMMA}, γ_T1={GAMMA_T1}")
    print(f"  Catalog: 2^N + 6 = {2 ** N + 6} candidates")
    print(f"  Total Pauli observables (non-identity): {4 ** N - 1}")
    print()

    for label, terms in cases:
        H = fw._build_bilinear(N, bonds, terms)
        result = fw.recommend_initial_state(
            H, [GAMMA] * N, N, gamma_t1_l=[GAMMA_T1] * N, top_k=5,
        )
        best = result['best']
        max_protection_frac = best['n_protected'] / (4 ** N - 1)
        print(f"H = J · Σ_bonds [{label}]")
        print(f"  Best: {best['label']} → {best['n_protected']} protected "
              f"({max_protection_frac:.0%})")
        print(f"  Top 5:")
        for lbl, n, _ in result['top_k']:
            print(f"    {lbl:<20s}  {n:>3d} protected")
        print()


if __name__ == "__main__":
    main()
