#!/usr/bin/env python3
"""Combined protection analysis: how do multiple Z₂ symmetries layer?

For a chosen (H, ρ_0), pi_protected_observables counts ALL Pauli
observables that stay strictly zero forever (cluster cancellation in
L_pauli's eigenbasis). This count is the upper bound — it includes
contributions from Y-parity, bit_a, K, and any other symmetries that
happen to apply.

For the cockpit's view of which symmetries individually 'activate':

  ρ_0 = |+−+⟩  (Y-content=0, bit_a-mixed):
    Y-parity activates (28 Y-odd protected)
    bit_a doesn't activate (mixed)

  ρ_0 = |010⟩  (Y-content=0, bit_a-pure-even):
    Y-parity activates (28 Y-odd protected)
    bit_a activates (32 bit_a-odd protected)
    Combined union: 28 + 32 − overlap = 47 expected

  ρ_0 = |+++⟩  (Y-content=0, bit_a mixed via X):
    Y-parity activates
    bit_a doesn't (mixed)

This script:
  1. For each (H, ρ_0) pair, compute pi_protected count + protected set
  2. Compute Y-odd protected set + bit_a-odd protected set (independent)
  3. Verify: pi_protected ⊇ Y-odd ∪ bit_a-odd
  4. Quantify how much extra protection beyond symmetry analysis
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


def get_protected_set_full(H, gamma_l, gamma_t1_l, rho_0, N,
                            threshold=1e-9, cluster_tol=1e-8):
    """Cluster-cancellation protected Pauli labels under L_full."""
    if any(g != 0 for g in gamma_t1_l):
        L = fw.lindbladian_z_plus_t1(H, gamma_l, gamma_t1_l)
    else:
        L = fw.lindbladian_z_dephasing(H, gamma_l)
    M_basis = fw._vec_to_pauli_basis_transform(N)
    L_pauli = (M_basis.conj().T @ L @ M_basis) / (2 ** N)
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
            if not used[j] and abs(evals[j] - evals[i]) < cluster_tol:
                cl.append(j)
                used[j] = True
        clusters.append(cl)

    protected = set()
    for alpha in range(1, 4 ** N):
        max_S = max(
            abs(sum(V[alpha, k] * c[k] for k in cl)) for cl in clusters
        )
        if max_S < threshold:
            label = ''.join(fw.PAULI_LABELS[idx]
                             for idx in fw._k_to_indices(alpha, N))
            protected.add(label)
    return protected


def y_odd_set(N):
    """Y-parity-odd Pauli labels (#Y mod 2 = 1)."""
    out = set()
    for alpha in range(1, 4 ** N):
        idx = fw._k_to_indices(alpha, N)
        n_y = sum(1 for i in idx if i == (1, 1))
        if n_y % 2 == 1:
            label = ''.join(fw.PAULI_LABELS[i] for i in idx)
            out.add(label)
    return out


def bit_a_odd_set(N):
    """bit_a-parity-odd Pauli labels (#X+#Y mod 2 = 1)."""
    out = set()
    for alpha in range(1, 4 ** N):
        idx = fw._k_to_indices(alpha, N)
        n_xy = sum(1 for i in idx if i[0] == 1)
        if n_xy % 2 == 1:
            label = ''.join(fw.PAULI_LABELS[i] for i in idx)
            out.add(label)
    return out


def main():
    N = 3
    GAMMA, GAMMA_T1, J = 0.1, 0.01, 1.0
    bonds = [(i, i + 1) for i in range(N - 1)]

    # Initial states
    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    zero = np.array([1, 0], dtype=complex)
    one = np.array([0, 1], dtype=complex)

    # |+−+⟩
    psi_xneel = np.kron(plus, np.kron(minus, plus))
    rho_xneel = np.outer(psi_xneel, psi_xneel.conj())

    # |010⟩
    psi_010 = np.kron(zero, np.kron(one, zero))
    rho_010 = np.outer(psi_010, psi_010.conj())

    # |+++⟩
    psi_ppp = np.kron(plus, np.kron(plus, plus))
    rho_ppp = np.outer(psi_ppp, psi_ppp.conj())

    # |GHZ⟩ = (|000⟩+|111⟩)/√2
    psi_ghz = (np.kron(zero, np.kron(zero, zero))
                + np.kron(one, np.kron(one, one))) / math.sqrt(2)
    rho_ghz = np.outer(psi_ghz, psi_ghz.conj())

    initial_states = {
        '|+−+⟩': rho_xneel,
        '|010⟩': rho_010,
        '|+++⟩': rho_ppp,
        '|GHZ⟩': rho_ghz,
    }

    cases = [
        ('truly XX+YY', [('X', 'X', J), ('Y', 'Y', J)]),
        ('XY+YX',       [('X', 'Y', J), ('Y', 'X', J)]),
        ('IY+YI',       [('I', 'Y', J), ('Y', 'I', J)]),
        ('XX+ZZ',       [('X', 'X', J), ('Z', 'Z', J)]),  # Heisenberg-like
        ('XX+YZ',       [('X', 'X', J), ('Y', 'Z', J)]),
    ]

    Y_odd = y_odd_set(N)
    A_odd = bit_a_odd_set(N)
    print(f"At N={N}: |Y-odd| = {len(Y_odd)}, |bit_a-odd| = {len(A_odd)}, "
          f"|Y∪bit_a| = {len(Y_odd | A_odd)}, |Y∩bit_a| = {len(Y_odd & A_odd)}")
    print()

    for state_label, rho_0 in initial_states.items():
        print(f"=" * 100)
        print(f"Initial state: {state_label}")
        print(f"=" * 100)
        # Compute Y-content and bit_a-content of ρ_0
        rho_pauli = fw.pauli_basis_vector(rho_0, N)
        bit_a_arr, _ = fw._bit_a_b_classify_paulis(N)
        y_arr = np.zeros(4 ** N, dtype=int)
        for alpha in range(4 ** N):
            idx = fw._k_to_indices(alpha, N)
            n_y = sum(1 for i in idx if i == (1, 1))
            y_arr[alpha] = n_y % 2
        rho_Y_odd = float(np.linalg.norm(rho_pauli[y_arr == 1]))
        rho_A_odd = float(np.linalg.norm(rho_pauli[bit_a_arr == 1]))
        print(f"  Y-odd content of ρ_0: {rho_Y_odd:.4f}, "
              f"bit_a-odd content: {rho_A_odd:.4f}")
        print()

        print(f"  {'case':<14s}  {'pi_prot':>7s}  {'Y∩prot':>7s}  "
              f"{'A∩prot':>7s}  {'Y∪A':>4s}  {'Y∪A∩prot':>9s}  "
              f"{'extra':>5s}  {'rating':<25s}")
        print('-' * 105)

        for label, terms in cases:
            H = fw._build_bilinear(N, bonds, terms)
            prot = get_protected_set_full(
                H, [GAMMA] * N, [GAMMA_T1] * N, rho_0, N
            )
            n_prot = len(prot)

            # Intersection with Y-odd, bit_a-odd
            n_Y_in_prot = len(Y_odd & prot)
            n_A_in_prot = len(A_odd & prot)
            n_union_in_prot = len((Y_odd | A_odd) & prot)
            extra = n_prot - n_union_in_prot  # protected beyond Y∪bit_a

            # Rating: how much of the symmetry-predicted protection actually
            # activates?
            possible_union_size = len(Y_odd | A_odd) if (rho_Y_odd < 1e-10
                                                          and rho_A_odd < 1e-10) \
                                    else (len(Y_odd) if rho_Y_odd < 1e-10
                                            else (len(A_odd) if rho_A_odd < 1e-10
                                                    else 0))
            if n_prot >= possible_union_size and possible_union_size > 0:
                rating = "saturates symmetry"
            elif n_prot >= len(Y_odd) and rho_Y_odd < 1e-10:
                rating = "Y-parity active"
            elif n_prot > 0:
                rating = "partial"
            else:
                rating = "no protection"
            print(f"  {label:<14s}  {n_prot:>7d}  {n_Y_in_prot:>7d}  "
                  f"{n_A_in_prot:>7d}  {len(Y_odd | A_odd):>4d}  "
                  f"{n_union_in_prot:>9d}  {extra:>5d}  {rating:<25s}")

    print()
    print("Reading guide:")
    print("  pi_prot:    total protected (cluster cancellation, includes all symmetries)")
    print("  Y∩prot:     protected observables that are Y-odd")
    print("  A∩prot:     protected observables that are bit_a-odd")
    print("  Y∪A:        Pauli labels in Y-odd ∪ bit_a-odd (universal upper bound)")
    print("  Y∪A∩prot:   Y-odd OR bit_a-odd Paulis that are actually protected")
    print("  extra:      protected observables NOT in the Y∪A union (other symmetries)")


if __name__ == "__main__":
    main()
