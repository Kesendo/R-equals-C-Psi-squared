"""Verification of the F112 proof sketch.

Proof claim (to verify numerically):
  1. For c with bit_b-homogeneous Pauli support (all c's Pauli strings share
     the same bit_b = #Y + #Z mod 2), np.kron(c, c.conj()) lies entirely in
     the Pi^2-conjugation +1 eigenspace.
  2. Pi^2-conjugation +1 eigenspace = Pi-conjugation {+1, -1} eigenspaces
     together. So np.kron(c, c.conj()) has ZERO content in Pi +i or -i
     eigenspaces when c is bit_b-homogeneous.
  3. Hence dissipator contributes NOTHING to M_plus_half or M_minus_half;
     these come ONLY from the Hamiltonian part L_H.
  4. For Hermitian H (real Pauli coefficient expansion), L_H is real in the
     appropriate basis, hence M_+i and M_-i are complex conjugates,
     equal in Frobenius norm.

Tests:
  A) Verify claim 1 numerically: compute ||(I - P_{Pi^2=+1}) np.kron(c, c*)||
     for bit_b-homogeneous c, expect ~0; for mixed-bit_b c, expect > 0.
  B) Verify the M_plus_half / M_minus_half balance reduces to the L_H
     part alone when c is bit_b-homogeneous.
"""

import sys
sys.path.insert(0, 'simulations')

import numpy as np
from itertools import product

from framework.symmetry import build_pi_full
from framework.diagnostics.polarity_coordinates import polarity_coordinates_from_L
from framework.pauli import _vec_to_pauli_basis_transform

PAULI_2X2 = {
    'I': np.eye(2, dtype=complex),
    'X': np.array([[0, 1], [1, 0]], dtype=complex),
    'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
    'Z': np.array([[1, 0], [0, -1]], dtype=complex),
}


def pauli_string(letters):
    op = PAULI_2X2[letters[0]]
    for l in letters[1:]:
        op = np.kron(op, PAULI_2X2[l])
    return op


def bit_b(p):
    return (p.count('Y') + p.count('Z')) % 2


def project_Pi2(M, Pi, sign):
    """Project M onto Pi^2-conjugation eigenvalue 'sign' (= +1 or -1).

    P_sign(M) = (M + sign * Pi^2 M Pi^-2) / 2
    """
    Pi_inv = Pi.conj().T
    Pi2 = Pi @ Pi
    Pi2_inv = Pi_inv @ Pi_inv
    return (M + sign * Pi2 @ M @ Pi2_inv) / 2


def L_vec_to_pauli(L_vec, N):
    T = _vec_to_pauli_basis_transform(N)
    return (T.conj().T @ L_vec @ T) / (2 ** N)


def test_dissipator_Pi2_homogeneity(N, c_paulis, label, Pi):
    """Build c = sum random_complex * P over c_paulis; test np.kron(c, c.conj())
    Pi^2-conjugation content split."""
    d = 2 ** N
    rng = np.random.default_rng(seed=hash(label) & 0xffff_ffff)
    coeffs = rng.normal(size=len(c_paulis)) + 1j * rng.normal(size=len(c_paulis))
    c = sum(coef * pauli_string(p) for coef, p in zip(coeffs, c_paulis))
    diss_vec = np.kron(c, c.conj())
    diss_pauli = L_vec_to_pauli(diss_vec, N)

    diss_Pi2_plus = project_Pi2(diss_pauli, Pi, sign=+1)
    diss_Pi2_minus = project_Pi2(diss_pauli, Pi, sign=-1)

    norm_total = np.sum(np.abs(diss_pauli) ** 2)
    norm_plus = np.sum(np.abs(diss_Pi2_plus) ** 2)
    norm_minus = np.sum(np.abs(diss_Pi2_minus) ** 2)

    pct_plus = 100 * norm_plus / max(norm_total, 1e-15)
    pct_minus = 100 * norm_minus / max(norm_total, 1e-15)

    bit_bs = [bit_b(p) for p in c_paulis]
    bit_b_set = sorted(set(bit_bs))
    print(f"  {label:<50}  bit_b={bit_b_set}  ||diss||^2={norm_total:>10.4f}  "
          f"Pi^2=+1: {pct_plus:>6.2f}%   Pi^2=-1: {pct_minus:>6.2f}%")


def main():
    print("=" * 100)
    print("Verify claim 1: bit_b-homogeneous c => np.kron(c, c*) is 100% in Pi^2-conj +1 eigenspace")
    print("=" * 100)

    for N in [2, 3]:
        d = 2 ** N
        Pi = build_pi_full(N)
        print(f"\n--- N = {N} ---\n")

        # Build test sets: 5 bit_b=0 c's, 5 bit_b=1 c's, 5 mixed-bit_b c's
        all_paulis = list(product(['I', 'X', 'Y', 'Z'], repeat=N))
        bit_b_0 = [p for p in all_paulis if bit_b(p) == 0]
        bit_b_1 = [p for p in all_paulis if bit_b(p) == 1]

        rng = np.random.default_rng(seed=N * 100)

        # bit_b=0 homogeneous c's, varying k_pauli
        for k in [2, 3, 4]:
            for trial in range(3):
                indices = rng.choice(len(bit_b_0), size=min(k, len(bit_b_0)), replace=False)
                selected = [bit_b_0[i] for i in indices]
                names = ','.join(''.join(p) for p in selected)
                test_dissipator_Pi2_homogeneity(N, selected, f"bit_b=0 k={k}: {names}", Pi)

        print()
        # bit_b=1 homogeneous c's
        for k in [2, 3, 4]:
            for trial in range(3):
                indices = rng.choice(len(bit_b_1), size=min(k, len(bit_b_1)), replace=False)
                selected = [bit_b_1[i] for i in indices]
                names = ','.join(''.join(p) for p in selected)
                test_dissipator_Pi2_homogeneity(N, selected, f"bit_b=1 k={k}: {names}", Pi)

        print()
        # mixed bit_b c's
        for trial in range(5):
            # 1 from bit_b=0, 1 from bit_b=1
            i0 = rng.integers(0, len(bit_b_0))
            i1 = rng.integers(0, len(bit_b_1))
            selected = [bit_b_0[i0], bit_b_1[i1]]
            names = ','.join(''.join(p) for p in selected)
            test_dissipator_Pi2_homogeneity(N, selected, f"MIXED bit_b k=2: {names}", Pi)


if __name__ == '__main__':
    main()
