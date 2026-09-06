#!/usr/bin/env python3
"""Demonstrate a sufficient operator identity for spectral palindromes.

For invertible Q, ``Q L Q^-1 + L + 2 c I = 0`` is sufficient for the
algebraic-multiplicity spectrum to be invariant under ``lambda -> -2c-lambda``.
The spectral symmetry alone is not a converse; it does not recover Q or the
operator identity.
"""

from itertools import product

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from scipy.spatial import cKDTree


def multiset_match(values, targets, tol=1e-7):
    """Return a global one-to-one tolerance matching and its worst residual."""
    values = np.asarray(values, dtype=complex)
    targets = np.asarray(targets, dtype=complex)
    if values.shape != targets.shape:
        return False, float("inf")
    value_points = np.column_stack((values.real, values.imag))
    target_points = np.column_stack((targets.real, targets.imag))
    neighbours = cKDTree(value_points).query_ball_point(target_points, tol)
    rows, cols = [], []
    for row, candidates in enumerate(neighbours):
        rows.extend([row] * len(candidates))
        cols.extend(candidates)
    graph = csr_matrix(
        (np.ones(len(rows), dtype=np.int8), (rows, cols)),
        shape=(len(targets), len(values)),
    )
    assignment = maximum_bipartite_matching(graph, perm_type="column")
    if np.count_nonzero(assignment >= 0) != len(values):
        return False, float("inf")
    # ``perm_type="column"`` returns the matched value-column for each
    # target-row.  This orientation matters for cycles longer than two.
    residual = float(np.max(np.abs(targets - values[assignment])))
    return residual <= tol, residual


def spectrum_palindromic(matrix, centre, tol=1e-7):
    values = np.linalg.eigvals(matrix)
    targets = -2 * centre - values
    return multiset_match(values, targets, tol)


def abstract_demo(n=12, seed=0):
    rng = np.random.RandomState(seed)
    q = np.diag([1.0] * (n // 2) + [-1.0] * (n - n // 2))
    centre = 1.0
    raw = rng.randn(n, n) + 1j * rng.randn(n, n)
    odd = (raw - q @ raw @ q) / 2
    generator = odd - centre * np.eye(n)
    residual = np.linalg.norm(q @ generator @ q + generator + 2 * centre * np.eye(n))
    palindrome, worst = spectrum_palindromic(generator, centre)
    return residual, palindrome, worst


def neural_demo(n=20, seed=42):
    """Construct an algebraically constrained E/I-labelled Jacobian.

    The construction enforces the Q identity directly. The E/I column signs
    satisfy a Dale-style sign constraint, but that constraint alone would not
    enforce the Q-related magnitudes or the time-constant pairing.
    """
    rng = np.random.RandomState(seed)
    half = n // 2
    q = np.block([[np.zeros((half, half)), np.eye(half)],
                  [np.eye(half), np.zeros((half, half))]])
    centre = 0.15
    raw = rng.uniform(0.0, 0.2, (n, n))
    raw[:, half:] *= -1
    np.fill_diagonal(raw, 0.0)
    odd = (raw - q @ raw @ q) / 2
    generator = odd - centre * np.eye(n)
    residual = np.linalg.norm(q @ generator @ q + generator + 2 * centre * np.eye(n))
    involution_residual = np.linalg.norm(q @ q - np.eye(n))
    palindrome, worst = spectrum_palindromic(generator, centre)
    return residual, involution_residual, palindrome, worst


def jordan_negative_control():
    """Same +/- eigenvalue multiset, incompatible Jordan block sizes."""
    matrix = np.array([[1.0, 1.0, 0.0, 0.0],
                       [0.0, 1.0, 0.0, 0.0],
                       [0.0, 0.0, -1.0, 0.0],
                       [0.0, 0.0, 0.0, -1.0]])
    palindrome, worst = spectrum_palindromic(matrix, 0.0)
    # Similarity to -L would preserve geometric multiplicities. Here +1 has
    # nullity 1 in L-I while -1 has nullity 2 in L+I.
    plus_nullity = 4 - np.linalg.matrix_rank(matrix - np.eye(4))
    minus_nullity = 4 - np.linalg.matrix_rank(matrix + np.eye(4))
    no_similarity = plus_nullity != minus_nullity
    return palindrome and no_similarity, worst, plus_nullity, minus_nullity


def quantum_demo(n=2, gamma=0.7, coupling=1.0):
    """Test the actual order-four Pauli-space Pi and its operator identity."""
    paulis = (
        np.eye(2, dtype=complex),
        np.array([[0, 1], [1, 0]], dtype=complex),
        np.array([[0, -1j], [1j, 0]], dtype=complex),
        np.array([[1, 0], [0, -1]], dtype=complex),
    )

    def tensor(items):
        result = np.array([[1.0]], dtype=complex)
        for item in items:
            result = np.kron(result, item)
        return result

    labels = list(product(range(4), repeat=n))
    basis = [tensor([paulis[index] for index in label]) for label in labels]
    dimension = 2**n
    hamiltonian = np.zeros((dimension, dimension), dtype=complex)
    for site in range(n - 1):
        for pauli_index in (1, 2, 3):
            factors = [paulis[0]] * n
            factors[site] = paulis[pauli_index]
            factors[site + 1] = paulis[pauli_index]
            hamiltonian += coupling * tensor(factors)

    generator = np.zeros((4**n, 4**n), dtype=complex)
    for column, operator in enumerate(basis):
        image = -1j * (hamiltonian @ operator - operator @ hamiltonian)
        for site in range(n):
            factors = [paulis[0]] * n
            factors[site] = paulis[3]
            z_site = tensor(factors)
            image += gamma * (z_site @ operator @ z_site - operator)
        for row, test_operator in enumerate(basis):
            generator[row, column] = np.trace(test_operator.conj().T @ image) / dimension

    pi_one = np.zeros((4, 4), dtype=complex)
    pi_one[1, 0] = 1       # I -> X
    pi_one[0, 1] = 1       # X -> I
    pi_one[3, 2] = 1j      # Y -> i Z
    pi_one[2, 3] = 1j      # Z -> i Y
    pi = tensor([pi_one] * n)
    identity = np.eye(4**n)
    centre = n * gamma
    operator_residual = np.linalg.norm(
        pi @ generator @ np.linalg.inv(pi) + generator + 2 * centre * identity
    )
    fourth_power_residual = np.linalg.norm(np.linalg.matrix_power(pi, 4) - identity)
    second_power_residual = np.linalg.norm(pi @ pi - identity)
    palindrome, worst = spectrum_palindromic(generator, centre)
    return operator_residual, fourth_power_residual, second_power_residual, palindrome, worst


def main():
    print("A SUFFICIENT OPERATOR IDENTITY FOR A SPECTRAL PALINDROME")
    print("This implication is not a converse: a paired eigenvalue multiset need not supply Q.\n")

    residual, palindrome, worst = abstract_demo()
    print(f"Abstract sufficient operator identity: residual={residual:.2e}, "
          f"palindrome={palindrome}, spectral residual={worst:.2e}")

    residual, involution, palindrome, worst = neural_demo()
    print(f"Constrained neural construction: identity={residual:.2e}, Q^2-I={involution:.2e}, "
          f"palindrome={palindrome}, spectral residual={worst:.2e}")
    print("Dale constraints alone are insufficient; the Q-related entries were imposed algebraically.\n")

    negative, worst, plus_nullity, minus_nullity = jordan_negative_control()
    print(f"Jordan negative control: spectral palindrome={negative}, residual={worst:.2e}, "
          f"geometric multiplicities (+1,-1)=({plus_nullity},{minus_nullity})")
    print("The mismatched Jordan structures rule out similarity to the reflected generator.\n")

    identity, fourth, second, palindrome, worst = quantum_demo()
    print(f"Quantum operator identity tested: residual={identity:.2e}, Pi^4-I={fourth:.2e}, "
          f"Pi^2-I={second:.2e}, palindrome={palindrome}, spectral residual={worst:.2e}")
    print("This Pi is tested as order four; it is not presented as an involution.")


if __name__ == "__main__":
    main()
