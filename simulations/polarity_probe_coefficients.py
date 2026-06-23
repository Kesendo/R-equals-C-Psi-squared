"""Polarity probe #11: which coefficient choices preserve balance.

Probe 10 found that c = P_alpha + i*P_beta preserves balance for ALL 136
Pauli-string pairs at N=2 (with coefficients exactly 1 and i). Probe 9
found that random complex coefficients SOMETIMES break balance.

This probe fixes a single Pauli pair and sweeps coefficient choices:
  - (1, 1)        real,real
  - (1, i)        real,pure-imag    (probe 10's choice)
  - (1, -1)       real,real-negative
  - (1, -i)       real,pure-imag-negative
  - (i, 1)        pure-imag,real
  - (1, 1+i)      real, mixed-complex
  - (1+i, 1-i)    mixed-complex
  - (random, random) fully random

For each coefficient choice and each of a few representative pairs.

Hypothesis: balance preserved iff each coefficient a_alpha is a "phase-locked"
multiple, i.e., a_alpha = e^(i*theta) * r_alpha with all r_alpha real (or
some other phase-consistency property).
"""

import sys
sys.path.insert(0, 'simulations')

import numpy as np
from itertools import product

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


def L_vec_to_pauli(L_vec, N):
    T = _vec_to_pauli_basis_transform(N)
    return (T.conj().T @ L_vec @ T) / (2 ** N)


def test_pair(N, p_alpha, p_beta, a, b, H, sigma):
    """Test c = a * P_alpha + b * P_beta."""
    d = 2 ** N
    c = a * pauli_string(p_alpha) + b * pauli_string(p_beta)
    Id = np.eye(d, dtype=complex)
    L_vec = -1j * (np.kron(H, Id) - np.kron(Id, H.T)) + 0.1 * np.kron(c, c.conj())
    L_pauli = L_vec_to_pauli(L_vec, N)
    result = polarity_coordinates_from_L(L_pauli, N, sigma)
    ns_M = result['norm_sq']['M']
    rel_asym = abs(result['asymmetry']) / max(ns_M, 1e-15)
    return ns_M, rel_asym


def main():
    N = 2
    d = 2 ** N
    sigma = 0.1
    H = np.zeros((d, d), dtype=complex)
    for letter in ['X', 'Y', 'Z']:
        H = H + pauli_string([letter, letter])

    # Representative pairs (Pauli strings as tuples)
    pairs = [
        (('X', 'I'), ('Y', 'I')),  # Klein same site
        (('X', 'I'), ('Z', 'I')),  # Klein different
        (('X', 'X'), ('Y', 'Y')),  # Klein (0,0) both
        (('X', 'Y'), ('Y', 'X')),  # Klein (1,1) both, y_par=1 both
        (('I', 'I'), ('X', 'X')),  # Identity + Pauli
    ]

    # Coefficient choices
    rng = np.random.default_rng(seed=2026)
    coefficient_choices = [
        ("a=1,    b=1   ", 1.0,    1.0),
        ("a=1,    b=-1  ", 1.0,   -1.0),
        ("a=1,    b=i   ", 1.0,    1.0j),
        ("a=1,    b=-i  ", 1.0,   -1.0j),
        ("a=i,    b=1   ", 1.0j,   1.0),
        ("a=1,    b=1+i ", 1.0,    1.0 + 1.0j),
        ("a=1+i,  b=1-i ", 1.0+1.0j, 1.0-1.0j),
        ("a=2,    b=3i  ", 2.0,    3.0j),
        ("a=2,    b=3   ", 2.0,    3.0),
        ("a=2+i,  b=1+3i", 2.0+1.0j, 1.0+3.0j),
    ]
    for trial in range(3):
        rval = complex(rng.normal(), rng.normal())
        cval = complex(rng.normal(), rng.normal())
        coefficient_choices.append((f"a=rand{trial}", rval, cval))

    print(f"N = {N}, H = XX + YY + ZZ Heisenberg, gamma = 0.1, sigma = {sigma}")
    print()
    print(f"For each pair (P_alpha, P_beta), sweep coefficient choices.")
    print()

    for p_a, p_b in pairs:
        label = f"c = a*{''.join(p_a)} + b*{''.join(p_b)}"
        print(f"=" * 80)
        print(label)
        print(f"=" * 80)
        print(f"{'coefficients':<22}  {'||M||^2':>12}  {'rel_asym':>12}  status")
        print("-" * 75)
        for coef_label, a, b in coefficient_choices:
            ns_M, rel_asym = test_pair(N, p_a, p_b, a, b, H, sigma)
            marker = "BALANCE" if rel_asym < 1e-10 else "BROKEN "
            print(f"{coef_label:<22}  {ns_M:>12.4f}  {rel_asym:>12.4e}  [{marker}]")
        print()


if __name__ == '__main__':
    main()
