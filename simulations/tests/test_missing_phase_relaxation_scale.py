import math

import numpy as np
import pytest
import sympy

from simulations import missing_phase_relaxation_scale as mprs


def independent_action(h, gamma, x, seat=3):
    z = np.eye(7, dtype=complex)
    z[seat, seat] = -1.0
    return -1j * (h @ x - x @ h) + gamma * (z @ x @ z - x)


def row_stack(x):
    return np.asarray(x, dtype=complex).reshape(-1, order="C")


def test_two_sided_invariant_embedding_has_exact_double_cluster():
    epsilon = 0.125
    gamma = 0.3
    h = mprs.hopping_float(epsilon)
    k = mprs.k_reduced_float(epsilon, gamma)
    v = mprs.blind_vector_float(epsilon)

    eigenvalues, eigenvectors = np.linalg.eig(k)
    minus_index = np.argmin(abs(eigenvalues + 2j * math.sqrt(2)))
    plus_index = np.argmin(abs(eigenvalues - 2j * math.sqrt(2)))
    lambda_minus = eigenvalues[minus_index]
    lambda_plus = eigenvalues[plus_index]
    u_minus = eigenvectors[:, minus_index]
    u_plus = eigenvectors[:, plus_index]

    right_minus = np.outer(u_minus, np.conj(v))
    adjoint_minus = np.outer(v, np.conj(u_plus))
    right_plus = np.outer(u_plus, np.conj(v))
    adjoint_plus = np.outer(v, np.conj(u_minus))

    for x in (right_minus, adjoint_minus):
        assert np.linalg.norm(independent_action(h, gamma, x) - lambda_minus * x) < 2e-11
    for x in (right_plus, adjoint_plus):
        assert np.linalg.norm(independent_action(h, gamma, x) - lambda_plus * x) < 2e-11

    minus_cluster = np.column_stack((row_stack(right_minus), row_stack(adjoint_minus)))
    plus_cluster = np.column_stack((row_stack(right_plus), row_stack(adjoint_plus)))
    assert np.linalg.matrix_rank(minus_cluster, tol=1e-10) == 2
    assert np.linalg.matrix_rank(plus_cluster, tol=1e-10) == 2

    assert np.array_equal(mprs.row_stack_right_embedding(u_minus, v), row_stack(right_minus))
    assert np.array_equal(mprs.row_stack_adjoint_embedding(v, u_plus), row_stack(adjoint_minus))


@pytest.mark.parametrize(
    ("epsilon", "gamma"),
    [
        (sympy.Rational(1, 8), sympy.Rational(3, 10)),
        (sympy.Rational(-1, 10), sympy.Rational(3, 100)),
    ],
)
def test_exact_reduction_identity_for_basis_and_general_vectors(epsilon, gamma):
    h, z, v, k = mprs.exact_reduction_objects(epsilon, gamma)
    zero_vector = sympy.zeros(7, 1)
    zero_matrix = sympy.zeros(7, 7)

    assert h * v == zero_vector
    assert z * v == v

    vectors = [sympy.eye(7).col(index) for index in range(7)]
    vectors.append(sympy.Matrix([1, sympy.I, 2, 0, -1, 3 * sympy.I, 4]))
    for u in vectors:
        x = u * sympy.conjugate(v.T)
        target = k * u * sympy.conjugate(v.T)
        assert mprs.exact_a_action(h, z, gamma, x) - target == zero_matrix


def test_row_stack_mutation_distinguishes_the_invariant_side():
    epsilon = 0.125
    gamma = 0.3
    h = mprs.hopping_float(epsilon)
    z = mprs.watched_involution_float()
    independent_z = np.eye(7, dtype=complex)
    independent_z[3, 3] = -1.0
    assert np.array_equal(z, independent_z)

    k = mprs.k_reduced_float(epsilon, gamma)
    eigenvalues, eigenvectors = np.linalg.eig(k)
    index = np.argmin(abs(eigenvalues + 2j * math.sqrt(2)))
    eigenvalue = eigenvalues[index]
    u = eigenvectors[:, index]
    v = mprs.blind_vector_float(epsilon)
    generator = mprs.a_generator_float(epsilon, gamma)

    correct_row_stack = np.kron(u, np.conj(v))
    wrong_column_stack = np.kron(v, u)
    assert np.linalg.norm(generator @ correct_row_stack - eigenvalue * correct_row_stack) < 2e-11
    assert np.linalg.norm(generator @ wrong_column_stack - eigenvalue * wrong_column_stack) > 1e-4


def test_sign_seat_and_minus_identity_mutations_break_reduction():
    residuals = mprs.reduction_mutation_residuals(
        sympy.Rational(1, 8), sympy.Rational(3, 10)
    )
    assert residuals == {
        "commutator_sign": sympy.Rational(393267, 128),
        "operator_side": sympy.Rational(393267, 256),
        "watched_seat": sympy.Rational(7137, 400),
        "missing_minus_identity": sympy.Rational(2763, 200),
    }
