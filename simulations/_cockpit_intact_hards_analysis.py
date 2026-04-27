#!/usr/bin/env python3
"""Algebraic characterisation of the 6 Lebensader-intact cases at N=3.

The cockpit 120-enum (commit cbb985c) found 6 cases with rating 'intact'
(drop ≤ 1 AND tail > 0.05) under |+−+⟩, γ_deph=0.1, γ_T1=0.01:

  Soft intacts (predicted):
    IY+YI  drop=0  tail=0.090  (factorising local-Y)
    XY+YX  drop=1  tail=0.080  (bond-flipped Z-free)

  Hard intacts (NEW class):
    XY+YZ  drop=0  tail=0.860
    XY+ZY  drop=0  tail=0.730
    YX+YZ  drop=0  tail=0.730
    YX+ZY  drop=0  tail=0.860

This script:
  1. Dumps the exact Π-protected Pauli set for each intact case
  2. Computes structural invariants (bit_a, bit_b totals; bond-flip orbits;
     Y-content)
  3. Compares the 6 against representative non-intact cases (truly XX+YY,
     a fragile soft YZ+ZY, a fragile hard XX+YZ)
  4. Looks for the algebraic rule that distinguishes intacts.
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


def bit_a_b(label):
    """Per-letter bit_a, bit_b for a 2-Pauli term label like 'XY'."""
    a_bits = [fw.bit_a(c) for c in label]
    b_bits = [fw.bit_b(c) for c in label]
    return a_bits, b_bits


def main():
    N = 3
    GAMMA_DEPH = 0.1
    GAMMA_T1 = 0.01
    J = 1.0
    bonds = [(i, i + 1) for i in range(N - 1)]

    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi = np.kron(plus, np.kron(minus, plus))
    rho_0 = np.outer(psi, psi.conj())

    intacts = [
        ('IY+YI', [('I', 'Y'), ('Y', 'I')], 'soft'),
        ('XY+YX', [('X', 'Y'), ('Y', 'X')], 'soft'),
        ('XY+YZ', [('X', 'Y'), ('Y', 'Z')], 'hard'),
        ('XY+ZY', [('X', 'Y'), ('Z', 'Y')], 'hard'),
        ('YX+YZ', [('Y', 'X'), ('Y', 'Z')], 'hard'),
        ('YX+ZY', [('Y', 'X'), ('Z', 'Y')], 'hard'),
    ]
    controls = [
        ('truly XX+YY', [('X', 'X'), ('Y', 'Y')], 'truly'),
        ('YZ+ZY', [('Y', 'Z'), ('Z', 'Y')], 'soft'),  # fragile soft
        ('XX+YZ', [('X', 'X'), ('Y', 'Z')], 'hard'),  # fragile hard
    ]

    print(f"Algebraic characterisation of intact cases at N={N}, |+−+⟩")
    print()
    print(f"Structural invariants per case:")
    print(f"  {'label':<14s}  {'cat':<6s}  {'term1':<6s}  {'term2':<6s}  "
          f"{'#Y':>3s}  {'#X':>3s}  {'#Z':>3s}  {'#I':>3s}  "
          f"{'bit_a sum':>9s}  {'bit_b sum':>9s}  {'bond-flip':>9s}")
    print('-' * 110)

    all_cases = intacts + controls
    for label, terms, cat in all_cases:
        t1 = ''.join(terms[0])
        t2 = ''.join(terms[1])
        full = t1 + t2
        ny = full.count('Y')
        nx = full.count('X')
        nz = full.count('Z')
        ni = full.count('I')

        a1, b1 = bit_a_b(t1)
        a2, b2 = bit_a_b(t2)
        bit_a_sum = sum(a1) + sum(a2)
        bit_b_sum = sum(b1) + sum(b2)
        bond_flipped = 'yes' if t1[::-1] == t2 else 'no'

        marker = "★" if (label, cat) in [(c[0], c[2]) for c in intacts] else " "
        print(f"  {marker}{label:<13s}  {cat:<6s}  {t1:<6s}  {t2:<6s}  "
              f"{ny:>3d}  {nx:>3d}  {nz:>3d}  {ni:>3d}  "
              f"{bit_a_sum:>9d}  {bit_b_sum:>9d}  {bond_flipped:>9s}")

    print()
    print("Protected-set dump (which Pauli observables stay strictly zero under +T1):")
    print()
    for label, terms, cat in intacts + controls[:1]:  # add truly XX+YY
        bilinear = [(t[0], t[1], J) for t in terms]
        H = fw._build_bilinear(N, bonds, bilinear)
        # Use lindbladian_z_plus_t1 directly via raw computation
        L_t1 = fw.lindbladian_z_plus_t1(H, [GAMMA_DEPH] * N, [GAMMA_T1] * N)
        M_basis = fw._vec_to_pauli_basis_transform(N)
        L_pauli = (M_basis.conj().T @ L_t1 @ M_basis) / (2 ** N)
        evals, V = np.linalg.eig(L_pauli)
        Vinv = np.linalg.inv(V)
        rho_pauli = fw.pauli_basis_vector(rho_0, N)
        c = Vinv @ rho_pauli

        n_eig = len(evals)
        used = np.zeros(n_eig, dtype=bool)
        clusters = []
        for i in range(n_eig):
            if used[i]:
                continue
            cl = [i]
            used[i] = True
            for j in range(i + 1, n_eig):
                if not used[j] and abs(evals[j] - evals[i]) < 1e-8:
                    cl.append(j)
                    used[j] = True
            clusters.append(cl)

        protected_labels = []
        for alpha in range(1, 4 ** N):
            max_S = 0.0
            for cl in clusters:
                S = sum(V[alpha, k] * c[k] for k in cl)
                max_S = max(max_S, abs(S))
            if max_S < 1e-9:
                lbl = ''.join(fw.PAULI_LABELS[idx] for idx in fw._k_to_indices(alpha, N))
                protected_labels.append(lbl)

        n_prot = len(protected_labels)
        print(f"  {label:<14s} ({cat}, n_protected={n_prot})")
        print(f"    {protected_labels}")
        print()

    # Look for shared substructure across intact protected sets
    print()
    print("Shared protected observables across the 6 intact cases:")
    intact_protected = []
    for label, terms, cat in intacts:
        bilinear = [(t[0], t[1], J) for t in terms]
        H = fw._build_bilinear(N, bonds, bilinear)
        L_t1 = fw.lindbladian_z_plus_t1(H, [GAMMA_DEPH] * N, [GAMMA_T1] * N)
        M_basis = fw._vec_to_pauli_basis_transform(N)
        L_pauli = (M_basis.conj().T @ L_t1 @ M_basis) / (2 ** N)
        evals, V = np.linalg.eig(L_pauli)
        Vinv = np.linalg.inv(V)
        rho_pauli = fw.pauli_basis_vector(rho_0, N)
        c = Vinv @ rho_pauli
        n_eig = len(evals)
        used = np.zeros(n_eig, dtype=bool)
        clusters = []
        for i in range(n_eig):
            if used[i]:
                continue
            cl = [i]
            used[i] = True
            for j in range(i + 1, n_eig):
                if not used[j] and abs(evals[j] - evals[i]) < 1e-8:
                    cl.append(j); used[j] = True
            clusters.append(cl)

        prot_set = set()
        for alpha in range(1, 4 ** N):
            max_S = max(
                abs(sum(V[alpha, k] * c[k] for k in cl)) for cl in clusters
            )
            if max_S < 1e-9:
                prot_set.add(''.join(
                    fw.PAULI_LABELS[idx] for idx in fw._k_to_indices(alpha, N)
                ))
        intact_protected.append((label, prot_set))

    common = set.intersection(*(s for _, s in intact_protected))
    print(f"  Shared by ALL 6 intacts: {len(common)} observables")
    print(f"    {sorted(common)}")
    print()

    # Pairwise: shared between just-soft-intacts vs just-hard-intacts
    soft_common = intact_protected[0][1] & intact_protected[1][1]
    hard_common = (intact_protected[2][1] & intact_protected[3][1]
                    & intact_protected[4][1] & intact_protected[5][1])
    print(f"  Soft intacts (IY+YI ∩ XY+YX): {len(soft_common)} observables")
    print(f"    {sorted(soft_common)}")
    print()
    print(f"  Hard intacts (4-way intersection): {len(hard_common)} observables")
    print(f"    {sorted(hard_common)}")


if __name__ == "__main__":
    main()
