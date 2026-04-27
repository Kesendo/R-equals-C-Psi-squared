#!/usr/bin/env python3
"""Demo of chiral_panel(): K_full = ⊗_{odd i} Z_i sublattice symmetry.

For the XY chain (no ZZ), K_full anti-commutes with H_xy, so K is a
symmetry of L = −i[H_xy, ·] + L_dephasing (Z-dephasing commutes with K
trivially). This gives a chiral skeleton: Pauli observables that are
strictly zero by K-mismatch with ρ_0.

Tested against the 6 Lebensader-intact cases at N=3, |+−+⟩, plus the
truly XX+YY (K-odd, well-behaved) and YZ+ZY (K-mixed, no clean symmetry).
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
    J = 1.0
    bonds = [(i, i + 1) for i in range(N - 1)]

    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi = np.kron(plus, np.kron(minus, plus))
    rho_0 = np.outer(psi, psi.conj())

    cases = [
        ('truly XX+YY', [('X', 'X', J), ('Y', 'Y', J)]),
        ('IY+YI',       [('I', 'Y', J), ('Y', 'I', J)]),
        ('XY+YX',       [('X', 'Y', J), ('Y', 'X', J)]),
        ('YZ+ZY',       [('Y', 'Z', J), ('Z', 'Y', J)]),
        ('XY+YZ',       [('X', 'Y', J), ('Y', 'Z', J)]),
        ('XY+ZY',       [('X', 'Y', J), ('Z', 'Y', J)]),
        ('YX+YZ',       [('Y', 'X', J), ('Y', 'Z', J)]),
        ('YX+ZY',       [('Y', 'X', J), ('Z', 'Y', J)]),
    ]

    print(f"Chiral panel demo, N={N}, |+−+⟩")
    print(f"K_full = I ⊗ Z ⊗ I (Z on the odd-indexed site only)")
    print()

    print(f"  {'case':<14s}  {'K_status':<10s}  {'L K-sym?':<9s}  "
          f"{'w_+':>7s}  {'w_-':>7s}  {'ρ_0 K-eig?':<11s}  "
          f"{'#static-0':>10s}  {'#dyn-prot':>10s}")
    print('-' * 100)

    for label, terms in cases:
        H = fw._build_bilinear(N, bonds, terms)
        panel = fw.chiral_panel(H, rho_0, N)

        K_status = panel['K_status']
        K_sym_L = panel['K_symmetric_L']
        w_p = panel['rho_decomposition']['w_plus']
        w_m = panel['rho_decomposition']['w_minus']
        rho_eig = panel['rho_is_K_eigenstate']
        n_static = len(panel['static_zero_at_t0'])
        n_dyn = len(panel['chiral_protected_dynamic'])
        sym_str = "yes" if K_sym_L else "no"
        eig_str = "yes" if rho_eig else "no"
        print(f"  {label:<14s}  {K_status:<10s}  {sym_str:<9s}  "
              f"{w_p:>7.4f}  {w_m:>7.4f}  {eig_str:<11s}  "
              f"{n_static:>10d}  {n_dyn:>10d}")

    print()
    print("Reading guide:")
    print("  K_status: 'K-even' if [K,H]=0; 'K-odd' if {K,H}=0; 'K-mixed' otherwise.")
    print("  L K-sym: True if K is a symmetry of L (block-diagonalises ρ).")
    print("  w_+, w_-: Tr(ρ_±²), K-symmetric and K-antisymmetric weights of ρ_0.")
    print("  ρ_0 K-eig?: True iff w_+ = 0 OR w_- = 0 (ρ_0 is purely one parity).")
    print("  #static-0: Pauli observables with ⟨P⟩(0) = 0 by K-structure of ρ_0.")
    print("             Does NOT imply dynamical protection.")
    print("  #dyn-prot: Pauli observables guaranteed to STAY zero, requiring")
    print("             BOTH L K-sym = yes AND ρ_0 a K-eigenstate.")
    print()

    print("With |+−+⟩ as initial state: not a K-eigenstate (K|+−+⟩ = |+++⟩),")
    print("so #dyn-prot = 0 for ALL cases — chiral protection inactive.")
    print("To activate chiral protection, the initial state must be K-symmetric")
    print("(e.g., (|+−+⟩+|+++⟩)/√2) or K-antisymmetric.")
    print()

    # Demonstrate with a K-eigenstate initial
    print("Demonstration: K-symmetric initial state (|+−+⟩+|+++⟩)/√2")
    psi_K_sym = (np.kron(plus, np.kron(minus, plus))
                 + np.kron(plus, np.kron(plus, plus))) / math.sqrt(2)
    psi_K_sym = psi_K_sym / np.linalg.norm(psi_K_sym)
    rho_K_sym = np.outer(psi_K_sym, psi_K_sym.conj())

    print(f"  {'case':<14s}  {'K_status':<10s}  {'w_+':>7s}  {'w_-':>7s}  "
          f"{'ρ_0 K-eig?':<11s}  {'#dyn-prot':>10s}")
    for label, terms in [('truly XX+YY', [('X', 'X', J), ('Y', 'Y', J)]),
                          ('XY+YX',       [('X', 'Y', J), ('Y', 'X', J)]),
                          ('XY+YZ',       [('X', 'Y', J), ('Y', 'Z', J)])]:
        H = fw._build_bilinear(N, bonds, terms)
        panel = fw.chiral_panel(H, rho_K_sym, N)
        w_p = panel['rho_decomposition']['w_plus']
        w_m = panel['rho_decomposition']['w_minus']
        eig = "yes" if panel['rho_is_K_eigenstate'] else "no"
        n_dyn = len(panel['chiral_protected_dynamic'])
        print(f"  {label:<14s}  {panel['K_status']:<10s}  "
              f"{w_p:>7.4f}  {w_m:>7.4f}  {eig:<11s}  {n_dyn:>10d}")
    print()
    print("On the K-symmetric initial state:")
    print("  - K-odd Hamiltonians (truly XX+YY, XY+YX) give chiral protection")
    print("    of all 32 K-odd Pauli observables (ρ_0 K-symmetric, P K-odd → 0)")
    print("  - K-mixed Hamiltonians (XY+YZ) do not (K not a symmetry of L)")


if __name__ == "__main__":
    main()
