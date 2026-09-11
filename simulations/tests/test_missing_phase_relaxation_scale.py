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


def test_uniform_peripheral_census_is_complete_and_not_just_the_kernel():
    census = mprs.exact_uniform_peripheral_census(sympy.Rational(3, 10))
    assert census == {
        "blind_dimension": 3,
        "end_blind_algebra_dimension": 9,
        "identity_outer_dimension": 1,
        "kernel_dimension": 4,
        "peripheral_dimension": 10,
        "frequency_multiplicities": {
            "-4sqrt2i": 1,
            "-2sqrt2i": 2,
            "0": 4,
            "+2sqrt2i": 2,
            "+4sqrt2i": 1,
        },
    }
    assert census["peripheral_dimension"] != census["kernel_dimension"]


@pytest.mark.parametrize(
    "epsilon",
    [sympy.Rational(1, 8), sympy.Rational(-1, 8), sympy.Rational(1, 16)],
)
def test_punctured_certificate_has_only_the_two_stationary_zero_cost_modes(epsilon):
    certificate = mprs.exact_punctured_kernel_certificate(
        epsilon, sympy.Rational(3, 10)
    )
    assert certificate["kernel_dimension"] == 2
    assert certificate["stationary_rank"] == 2
    assert certificate["zero_cost_invariant_dimension"] == 2
    assert certificate["nonzero_imaginary_axis_dimension"] == 0
    polynomial = certificate["restricted_characteristic_polynomial"]
    assert polynomial.as_expr() == polynomial.gen**2


def _expected_effective_operators(gamma):
    root2 = sympy.sqrt(2)
    return {
        "first": {
            "-4sqrt2i": sympy.Matrix([[-sympy.I * root2]]),
            "-2sqrt2i": -sympy.I / root2 * sympy.eye(2),
            "0": sympy.zeros(4),
            "+2sqrt2i": sympy.I / root2 * sympy.eye(2),
            "+4sqrt2i": sympy.Matrix([[sympy.I * root2]]),
        },
        "second": {
            "-4sqrt2i": sympy.Matrix([[-gamma - 3 * root2 * sympy.I / 8]]),
            "-2sqrt2i": (-gamma / 2 - 3 * root2 * sympy.I / 16) * sympy.eye(2),
            "0": sympy.Matrix(
                [
                    [-gamma, 0, 0, gamma / 2],
                    [0, 0, 0, 0],
                    [0, 0, -gamma, gamma / 2],
                    [gamma / 2, 0, gamma / 2, -gamma / 2],
                ]
            ),
            "+2sqrt2i": (-gamma / 2 + 3 * root2 * sympy.I / 16) * sympy.eye(2),
            "+4sqrt2i": sympy.Matrix([[-gamma + 3 * root2 * sympy.I / 8]]),
        },
    }


def test_kato_feshbach_effective_operators_and_polynomials_are_exact():
    gamma = sympy.symbols("gamma", positive=True)
    actual = mprs.exact_effective_operators(gamma)
    expected = _expected_effective_operators(gamma)

    for order in ("first", "second"):
        for frequency, matrix in expected[order].items():
            assert actual[order][frequency] == matrix
            polynomial = actual["characteristic_polynomials"][order][frequency]
            assert polynomial == matrix.charpoly(polynomial.gen)

    zero_polynomial = actual["characteristic_polynomials"]["second"]["0"]
    lam = zero_polynomial.gen
    assert sympy.expand(
        zero_polynomial.as_expr()
        - lam**2 * (lam + gamma) * (2 * lam + 3 * gamma) / 2
    ) == 0


def test_kato_feshbach_projector_and_reduced_resolvent_identities_are_exact():
    gamma = sympy.symbols("gamma", positive=True)
    result = mprs.exact_effective_operators(gamma)
    for block in result["certificates"].values():
        l0 = block["generator"]
        projector = block["projector"]
        right = block["right_basis"]
        left = block["left_basis"]

        assert left.conjugate().T * right == sympy.eye(right.cols)
        assert projector * projector == projector
        assert l0 * projector == projector * l0
        assert block["basis_dimension"] == 49
        assert sum(len(part["indices"]) for part in block["resolvent_blocks"]) == 49
        assert len(
            {
                coordinate
                for part in block["resolvent_blocks"]
                for coordinate in part["indices"]
            }
        ) == 49
        for part in block["resolvent_blocks"]:
            size = len(part["indices"])
            # S=M+P has SP=PS=P.  Together with SS^-1=S^-1S=I and
            # P^2=P this is exactly the two-sided identity
            # M(S^-1-P)=(S^-1-P)M=I-P, without expanding the 16x16 inverse.
            assert part["shifted_is_invertible"]
            assert part["unshifted"] == part["shifted"] - part["projector"]
            assert part["shifted"] * part["projector"] == part["projector"]
            assert part["projector"] * part["shifted"] == part["projector"]


def test_effective_multiplicities_and_local_minimality_mutations_are_detected():
    gamma = sympy.symbols("gamma", positive=True)
    result = mprs.exact_effective_operators(gamma)
    expected = _expected_effective_operators(gamma)
    root2 = sympy.sqrt(2)
    expected_multiplicities = {
        "first": {
            "-4sqrt2i": {-root2 * sympy.I: 1},
            "-2sqrt2i": {-sympy.I / root2: 2},
            "0": {sympy.Integer(0): 4},
            "+2sqrt2i": {sympy.I / root2: 2},
            "+4sqrt2i": {root2 * sympy.I: 1},
        },
        "second": {
            "-4sqrt2i": {-gamma - 3 * root2 * sympy.I / 8: 1},
            "-2sqrt2i": {-gamma / 2 - 3 * root2 * sympy.I / 16: 2},
            "0": {sympy.Integer(0): 2, -gamma: 1, -3 * gamma / 2: 1},
            "+2sqrt2i": {-gamma / 2 + 3 * root2 * sympy.I / 16: 2},
            "+4sqrt2i": {-gamma + 3 * root2 * sympy.I / 8: 1},
        },
    }
    expected_geometric = {
        order: {
            frequency: dict(multiplicities)
            for frequency, multiplicities in by_frequency.items()
        }
        for order, by_frequency in expected_multiplicities.items()
    }

    for order in ("first", "second"):
        for frequency, matrix in result[order].items():
            assert matrix.eigenvals() == expected_multiplicities[order][frequency]
            for eigenvalue, geometric_expected in expected_geometric[order][frequency].items():
                geometric = len((matrix - eigenvalue * sympy.eye(matrix.rows)).nullspace())
                assert geometric == geometric_expected

    damping_mutant = result["first"]["-2sqrt2i"] - gamma * sympy.eye(2)
    assert damping_mutant.charpoly().as_expr() != result["characteristic_polynomials"]["first"]["-2sqrt2i"].as_expr()

    reactive_mutant = result["second"]["-4sqrt2i"].applyfunc(sympy.re)
    assert reactive_mutant.charpoly().as_expr() != result["characteristic_polynomials"]["second"]["-4sqrt2i"].as_expr()

    zero_minor = result["second"]["0"].extract([0, 1, 2], [0, 1, 2])
    assert zero_minor.charpoly().as_expr() != result["characteristic_polynomials"]["second"]["0"].as_expr()

    peripheral = mprs.exact_uniform_peripheral_census(sympy.Rational(3, 10))
    kernel_only = dict(peripheral)
    kernel_only["peripheral_dimension"] = kernel_only["kernel_dimension"]
    kernel_only["frequency_multiplicities"] = {"0": 4}
    assert kernel_only != peripheral
