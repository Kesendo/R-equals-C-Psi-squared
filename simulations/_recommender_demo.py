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
        pure_best = result['pure_protection_best']
        quantum_best = result.get('quantum_best')

        print(f"H = J · Σ_bonds [{label}]")
        print(f"  Theoretical max (pure state, N={N}): "
              f"{result['theoretical_max_pure']} ({result['saturation_status']})")
        print(f"  Combined best (with classical penalty): {best['label']} → "
              f"n_prot={best['n_protected']}, "
              f"bit_a_even={best['bit_a_even_fraction']:.3f}, "
              f"score={best['score']:.3f}")
        if quantum_best and quantum_best['label'] != best['label']:
            print(f"  Quantum best (no classical states): {quantum_best['label']} → "
                  f"n_prot={quantum_best['n_protected']}, "
                  f"bit_a_even={quantum_best['bit_a_even_fraction']:.3f}")
        for w in result['warnings']:
            print(f"  ⚠ {w}")
        print(f"  Top 5 by combined score:")
        for r in result['top_k']:
            classical_tag = " (classical)" if r['is_classical_diagonal'] else ""
            print(f"    {r['label']:<14s}  n_prot={r['n_protected']:>3d}  "
                  f"bit_a_even={r['bit_a_even_fraction']:.3f}  "
                  f"score={r['score']:>+5.3f}{classical_tag}")
        print()


if __name__ == "__main__":
    main()
