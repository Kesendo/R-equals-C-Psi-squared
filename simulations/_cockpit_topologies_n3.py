#!/usr/bin/env python3
"""Cockpit comparison across topologies at N=3.

The 4-panel cockpit was built and validated on the open chain. Most
primitives are topology-agnostic (Lebensader, cusp pattern, Y-parity
work directly on any L). Only the chiral K_full panel is chain-specific
(uses K = ⊗_{odd i} Z_i, the chain bipartition).

This script runs the cockpit on three topologies at N=3 with the same
6 'intact' bond bilinears + truly XX+YY:

  Chain:  bonds = [(0,1), (1,2)]            — bipartite, K = I⊗Z⊗I
  Ring:   bonds = [(0,1), (1,2), (2,0)]     — odd cycle, no bipartition
  Star:   bonds = [(0,1), (0,2)]            — bipartite, K_star = Z⊗I⊗I

For each (topology, H), report the cockpit's four panels and look for
patterns: which protections survive topology change?
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


def chiral_K_for_topology(topology, N):
    """Topology-specific chiral operator (or None if not bipartite).

    For chain N=3: K = I⊗Z⊗I (Z on the odd middle site)
    For star N=3 with center 0: K = Z⊗I⊗I (Z on center)
    For ring N=3: no bipartition (odd cycle), return None
    """
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    I2 = np.eye(2, dtype=complex)
    if topology == 'chain':
        return np.kron(I2, np.kron(Z, I2))
    elif topology == 'star':
        return np.kron(Z, np.kron(I2, I2))
    elif topology == 'ring':
        # Odd cycle: no bipartition. Could try K = Z⊗Z⊗Z (all sites)
        # but that doesn't anticommute cleanly with bond terms.
        return None
    else:
        raise ValueError(f"Unknown topology: {topology}")


def main():
    N = 3
    GAMMA_DEPH = 0.1
    GAMMA_T1 = 0.01
    J = 1.0

    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi = np.kron(plus, np.kron(minus, plus))
    rho_0 = np.outer(psi, psi.conj())

    topologies = {
        'chain': [(0, 1), (1, 2)],
        'ring':  [(0, 1), (1, 2), (2, 0)],
        'star':  [(0, 1), (0, 2)],
    }

    cases = [
        ('truly XX+YY', [('X', 'X', J), ('Y', 'Y', J)]),
        ('IY+YI',       [('I', 'Y', J), ('Y', 'I', J)]),
        ('XY+YX',       [('X', 'Y', J), ('Y', 'X', J)]),
        ('YZ+ZY',       [('Y', 'Z', J), ('Z', 'Y', J)]),
        ('XY+YZ',       [('X', 'Y', J), ('Y', 'Z', J)]),
    ]

    for topology, bonds in topologies.items():
        print(f"=" * 100)
        print(f"Topology: {topology}, bonds = {bonds}")
        print(f"=" * 100)
        print()

        K_top = chiral_K_for_topology(topology, N)

        print(f"  {'case':<14s}  {'Leb':<10s}  {'Π drop':>6s}  "
              f"{'Cusp pat':<11s}  {'Y-pres?':>7s}  {'Y-prot':>7s}  "
              f"{'K_stat':<10s}")
        print('-' * 90)

        for label, terms in cases:
            H = fw._build_bilinear(N, bonds, terms)
            panel = fw.cockpit_panel(
                H, [GAMMA_DEPH] * N, rho_0, N,
                gamma_t1_l=[GAMMA_T1] * N, t_max=8.0, dt=0.005,
            )
            leb = panel['lebensader']
            rating_short = leb['rating'].split(' ')[0]
            cusp_pat = panel['cusp']['pattern']
            yp = panel['y_parity']
            n_Y_prot = len(yp['Y_parity_protected'])
            y_pres = "yes" if yp['L_preserves_Y_parity'] else "no"

            # K-status using topology-specific K
            if K_top is not None:
                K_status = fw.k_classify_hamiltonian(H, N, K_full=K_top)
            else:
                K_status = "no-K"

            print(f"  {label:<14s}  {rating_short:<10s}  "
                  f"{leb['skeleton']['drop']:>6d}  "
                  f"{cusp_pat:<11s}  {y_pres:>7s}  {n_Y_prot:>7d}  "
                  f"{K_status:<10s}")
        print()

    print("=" * 100)
    print("Cross-topology reading:")
    print("=" * 100)
    print()
    print("Y-parity preservation is topology-agnostic — depends only on the")
    print("(H, dissipator) algebra, not on bond structure. If the 6 intact")
    print("Hamiltonians have '1 Y per term' on each bond, the Y-parity panel")
    print("should report yes regardless of how many bonds (chain has 2,")
    print("ring has 3, star has 2). Look for confirmation across topologies.")
    print()
    print("Chiral K depends on bipartition. Chain and star are bipartite")
    print("(distinct K_full). Ring at N=3 (odd cycle) has no clean K.")
    print()
    print("Lebensader rating may shift across topologies because the L")
    print("eigenvalue structure changes with bond count and connectivity,")
    print("which affects the Π-protected count under +T1.")


if __name__ == "__main__":
    main()
