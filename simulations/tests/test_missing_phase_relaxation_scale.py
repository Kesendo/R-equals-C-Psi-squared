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


def independent_generator(h, gamma, seat=3):
    generator = np.empty((49, 49), dtype=complex)
    for coordinate in range(49):
        basis = np.zeros((7, 7), dtype=complex)
        basis[np.unravel_index(coordinate, (7, 7), order="C")] = 1.0
        generator[:, coordinate] = row_stack(
            independent_action(h, gamma, basis, seat=seat)
        )
    return generator


def test_float_construction_and_generator_use_literal_scale_and_normalization():
    epsilon = 0.125
    gamma = 0.3
    expected_h = np.zeros((7, 7), dtype=complex)
    for site, bond in enumerate([2.25, 2.0, 2.0, 2.0, 2.0, 2.0]):
        expected_h[site, site + 1] = expected_h[site + 1, site] = bond
    assert np.array_equal(mprs.hopping_float(epsilon), expected_h)

    expected_v = np.array(
        [1.0, 0.0, -1.125, 0.0, 1.125, 0.0, -1.125], dtype=complex
    )
    expected_v /= np.linalg.norm(expected_v)
    actual_v = mprs.blind_vector_float(epsilon)
    assert np.linalg.norm(actual_v) == pytest.approx(1.0, abs=1e-15)
    assert np.array_equal(actual_v, expected_v)

    identity = np.eye(7, dtype=complex)
    z = np.eye(7, dtype=complex)
    z[3, 3] = -1.0
    expected_generator = (
        -1j * (np.kron(expected_h, identity) - np.kron(identity, expected_h.T))
        + gamma * (np.kron(z, z.T) - np.eye(49, dtype=complex))
    )
    assert np.array_equal(
        mprs.a_generator_float(epsilon, gamma), expected_generator
    )


def test_full_generator_matches_columnwise_action_oracle():
    epsilon = 0.125
    gamma = 0.3
    expected_h = np.zeros((7, 7), dtype=complex)
    for site, bond in enumerate([2.25, 2.0, 2.0, 2.0, 2.0, 2.0]):
        expected_h[site, site + 1] = expected_h[site + 1, site] = bond

    oracle = independent_generator(expected_h, gamma)
    assert np.linalg.norm(mprs.a_generator_float(epsilon, gamma) - oracle) < 1e-12


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


@pytest.mark.parametrize(
    ("epsilon", "gamma"),
    [
        (0.125, sympy.Rational(3, 10)),
        (sympy.Rational(1, 8), 0.3),
    ],
)
def test_exact_reduction_rejects_inexact_scalars(epsilon, gamma):
    with pytest.raises(TypeError):
        mprs.exact_reduction_objects(epsilon, gamma)


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


def test_characteristic_polynomial_matches_displayed_q_exactly():
    lam, r, gamma = sympy.symbols("lambda r gamma")
    determinant = sympy.expand(mprs.characteristic_polynomial(lam, r, gamma))
    expected_q = (
        lam**6
        + 2 * gamma * lam**5
        + (4 * r**2 + 20) * lam**4
        + (8 * gamma * r**2 + 24 * gamma) * lam**3
        + (64 * r**2 + 96) * lam**2
        + (64 * gamma * r**2 + 64 * gamma) * lam
        + 192 * r**2
        + 64
    )

    assert sympy.expand(determinant - lam * expected_q) == 0
    assert sympy.rem(
        determinant,
        lam,
        domain=sympy.QQ.frac_field(r, gamma),
    ) == 0


def test_uniform_characteristic_factor_and_target_roots_are_exact_and_simple():
    lam, gamma = sympy.symbols("lambda gamma")
    determinant = sympy.expand(mprs.characteristic_polynomial(lam, 1, gamma))
    expected_q = (lam**2 + 8) * (
        lam**4 + 2 * gamma * lam**3 + 16 * lam**2 + 16 * gamma * lam + 32
    )
    assert sympy.expand(determinant - lam * expected_q) == 0

    for root in (-2 * sympy.sqrt(2) * sympy.I, 2 * sympy.sqrt(2) * sympy.I):
        assert sympy.simplify(expected_q.subs(lam, root)) == 0
        assert sympy.simplify(sympy.diff(expected_q, lam).subs(lam, root)) != 0


def test_minus_branch_series_solves_through_cubic_but_not_quartic():
    lam, r, epsilon, gamma = sympy.symbols("lambda r epsilon gamma")
    branch = mprs.minus_branch_series(epsilon, gamma)
    expected = (
        -2 * sympy.sqrt(2) * sympy.I
        - sympy.I / sympy.sqrt(2) * epsilon
        - (gamma / 2 + 3 * sympy.sqrt(2) * sympy.I / 16) * epsilon**2
        - (gamma / 2 - 19 * sympy.sqrt(2) * sympy.I / 64) * epsilon**3
    )
    assert sympy.expand(branch - expected) == 0

    residual = sympy.expand(
        mprs.characteristic_polynomial(lam, r, gamma).subs(
            {lam: branch, r: 1 + epsilon}
        )
    )
    for order in range(4):
        assert sympy.simplify(residual.coeff(epsilon, order)) == 0
    assert sympy.simplify(residual.coeff(epsilon, 4)) != 0


def test_gap_series_and_first_asymmetry_term_are_exact():
    epsilon, gamma = sympy.symbols("epsilon gamma")
    gap = mprs.gap_series(epsilon, gamma)
    assert sympy.expand(gap - gamma * epsilon**2 / 2 - gamma * epsilon**3 / 2) == 0
    assert sympy.expand(gap - gap.subs(epsilon, -epsilon) - gamma * epsilon**3) == 0


def test_gap_series_is_the_negative_real_part_of_the_branch_germ():
    epsilon, gamma = sympy.symbols("epsilon gamma", real=True)
    branch = mprs.minus_branch_series(epsilon, gamma)
    gap = mprs.gap_series(epsilon, gamma)
    assert sympy.simplify(gap + sympy.re(branch)) == 0

    for real_term in (-gamma * epsilon**2 / 2, -gamma * epsilon**3 / 2):
        sign_flipped_branch = branch - 2 * real_term
        assert sympy.simplify(gap + sympy.re(sign_flipped_branch)) != 0


def test_each_displayed_branch_component_is_required_through_cubic_order():
    lam, r, epsilon, gamma = sympy.symbols("lambda r epsilon gamma")
    branch = mprs.minus_branch_series(epsilon, gamma)
    components = {
        "constant_imaginary": -2 * sympy.sqrt(2) * sympy.I,
        "linear_imaginary": -sympy.I / sympy.sqrt(2) * epsilon,
        "quadratic_real": -gamma * epsilon**2 / 2,
        "quadratic_imaginary": -3 * sympy.sqrt(2) * sympy.I * epsilon**2 / 16,
        "cubic_real": -gamma * epsilon**3 / 2,
        "cubic_imaginary": 19 * sympy.sqrt(2) * sympy.I * epsilon**3 / 64,
    }

    def has_residual_through_cubic(mutated_branch):
        residual = sympy.expand(
            mprs.characteristic_polynomial(lam, r, gamma).subs(
                {lam: mutated_branch, r: 1 + epsilon}
            )
        )
        return any(
            sympy.simplify(residual.coeff(epsilon, order)) != 0
            for order in range(4)
        )

    for component in components.values():
        assert has_residual_through_cubic(branch - component)

    for real_term in (components["quadratic_real"], components["cubic_real"]):
        assert has_residual_through_cubic(branch - 2 * real_term)
