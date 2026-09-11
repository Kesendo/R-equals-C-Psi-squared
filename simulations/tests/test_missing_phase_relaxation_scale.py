import math
import sys
from pathlib import Path

import numpy as np
import pytest
import sympy

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

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


def test_uniform_census_uses_ranks_and_rejects_a_duplicated_frequency_vector():
    gamma = sympy.Rational(3, 10)
    certificate = mprs.exact_uniform_peripheral_certificate(gamma)
    assert certificate["per_frequency_ranks"] == {
        "-4sqrt2i": 1,
        "-2sqrt2i": 2,
        "0": 4,
        "+2sqrt2i": 2,
        "+4sqrt2i": 1,
    }
    assert certificate["combined_operator_rank"] == 10
    assert certificate["zero_cost_invariant_dimension"] == 10
    assert certificate["combined_span_rank"] == 10

    mutated = {
        frequency: list(basis)
        for frequency, basis in certificate["operator_bases"].items()
    }
    mutated["-2sqrt2i"][1] = mutated["-2sqrt2i"][0]
    assert sum(len(basis) for basis in mutated.values()) == 10
    with pytest.raises(AssertionError, match="rank|span|complete"):
        mprs.certify_uniform_peripheral_bases(gamma, mutated)


@pytest.mark.parametrize(
    "gamma",
    [
        sympy.Integer(0),
        sympy.Rational(-3, 10),
        sympy.Symbol("gamma"),
        0.3,
    ],
)
def test_numeric_certificates_reject_nonpositive_or_inexact_gamma(gamma):
    with pytest.raises((TypeError, ValueError)):
        mprs.exact_uniform_peripheral_census(gamma)
    with pytest.raises((TypeError, ValueError)):
        mprs.exact_punctured_kernel_certificate(sympy.Rational(1, 8), gamma)


@pytest.mark.parametrize(
    "gamma",
    [
        sympy.Integer(0),
        sympy.Integer(-1),
        sympy.Rational(3, 10),
        sympy.Symbol("gamma"),
        0.3,
    ],
)
def test_symbolic_effective_certificate_requires_a_positive_real_indeterminate(gamma):
    with pytest.raises((TypeError, ValueError)):
        mprs.exact_effective_operators(gamma)


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
        assert block["denominator_nonzero_for_positive_gamma"]
        denominator_certificate = block["positive_denominator_certificate"]
        assert denominator_certificate["holds"]
        assert denominator_certificate["positive_real_root_count"] == 0
        if block["frequency"] == 0:
            assert block["resolvent_denominator"] == gamma**3
            assert denominator_certificate["zero_root_multiplicity"] == 3
            assert denominator_certificate["positive_axis_polynomial"].degree() == 0
        else:
            assert denominator_certificate[
                "common_real_zero_polynomial"
            ].degree() == 0
        assert sum(len(part["indices"]) for part in block["resolvent_blocks"]) == 49
        assert len(
            {
                coordinate
                for part in block["resolvent_blocks"]
                for coordinate in part["indices"]
            }
        ) == 49
        shifted_generator = l0 - block["frequency"] * sympy.eye(49)
        numerator = block["resolvent_numerator"]
        denominator = block["resolvent_denominator"]
        complement = sympy.eye(49) - projector
        assert block["shifted_off_block_residual"] == sympy.zeros(49)
        assert (
            shifted_generator * numerator - complement * denominator
        ).applyfunc(sympy.expand) == sympy.zeros(49)
        assert (
            numerator * shifted_generator - complement * denominator
        ).applyfunc(sympy.expand) == sympy.zeros(49)
        for part in block["resolvent_blocks"]:
            assert part["unshifted"] == part["shifted"] - part["projector"]
            assert part["shifted"] * part["projector"] == part["projector"]
            assert part["projector"] * part["shifted"] == part["projector"]


def test_off_block_mutation_escapes_the_old_local_check_but_fails_global_gate():
    gamma = sympy.symbols("gamma", positive=True)
    certificate = mprs.exact_effective_operators(gamma)["certificates"]["-2sqrt2i"]
    blocks = certificate["resolvent_blocks"]
    shifted_generator = (
        certificate["generator"] - certificate["frequency"] * sympy.eye(49)
    )
    mutant = shifted_generator.copy()
    source = blocks[0]["indices"][0]
    target = blocks[1]["indices"][0]
    mutant[source, target] += 1

    for part in blocks:
        indices = list(part["indices"])
        assert mutant.extract(indices, indices) == part["unshifted"]

    numerator = certificate["resolvent_numerator"]
    denominator = certificate["resolvent_denominator"]
    complement = sympy.eye(49) - certificate["projector"]
    assert (
        mutant * numerator - complement * denominator
    ).applyfunc(sympy.expand) != sympy.zeros(49)


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


def test_b_generator_is_priced_negative_adjoint_of_a():
    epsilon, gamma = 0.125, 0.3
    la = mprs.a_generator_float(epsilon, gamma)
    lb = mprs.b_generator_float(epsilon, gamma)
    assert np.linalg.norm(
        lb - (-2 * gamma * np.eye(49) - la.conj().T)
    ) < 1e-13


def test_b_boundary_has_no_peripheral_direction_and_control_can_create_one():
    boundary = mprs.exact_b_boundary_certificate()
    assert boundary["h_outer_times_u"] == sympy.Matrix([0, 4, 0, 0, 4, 0])
    assert boundary["has_peripheral_b_mode"] is False
    assert boundary["peripheral_dimension"] == 0

    control = mprs.exact_b_boundary_control()
    assert control["u_nonzero"] is True
    assert control["h_outer_times_u"] == sympy.zeros(6, 1)
    assert control["has_peripheral_b_mode"] is True
    assert control["peripheral_dimension"] == 1
    assert control["direct_generator_residual"] == sympy.zeros(49, 1)


def test_b_gap_uses_the_maximum_a_rate_not_the_a_gap():
    epsilon, gamma = 0.1, 0.3
    result = mprs.direct_float_gaps(epsilon, gamma)
    assert result["b_gap"] == pytest.approx(
        2 * gamma - result["max_a_rate"], rel=2e-10, abs=2e-12
    )
    assert abs(result["b_gap"] - (2 * gamma - result["a_gap"])) > 1e-3
    assert result["next_distinct_a_rate"] > result["a_gap"]


def test_each_b_rate_is_the_individually_paired_a_rate_with_the_price():
    result = mprs.direct_float_gaps(0.1, 0.3)
    paired_rates = result["paired_rate_residuals"]
    assert len(paired_rates) == 49
    assert max(paired_rates) < 2e-12


def test_a_rate_classifier_rejects_an_unstable_unclassified_eigenvalue():
    values = np.linalg.eigvals(mprs.a_generator_float(0.1, 0.3))
    tolerance = 128 * np.finfo(float).eps * max(
        1.0, np.linalg.norm(mprs.a_generator_float(0.1, 0.3), ord=2)
    )
    values[np.argmin(values.real)] = 1e-3 + 2j
    with pytest.raises(AssertionError, match="unstable|49|classif"):
        mprs.classify_a_rates(
            values,
            tolerance,
            expected_peripheral=2,
            expected_stationary=2,
        )


def test_direct_gap_route_cannot_bypass_complete_rate_classification(monkeypatch):
    original_eigvals = np.linalg.eigvals
    call_count = 0

    def eigenvalues_with_one_unstable_a_mode(matrix):
        nonlocal call_count
        values = original_eigvals(matrix)
        if call_count == 0:
            values[np.argmin(values.real)] = 1e-3 + 2j
        call_count += 1
        return values

    monkeypatch.setattr(
        mprs.np.linalg,
        "eigvals",
        eigenvalues_with_one_unstable_a_mode,
    )
    with pytest.raises(AssertionError, match="unstable|49|classif"):
        mprs.direct_float_gaps(0.1, 0.3)


@pytest.mark.parametrize("r", [0.9, 1.1])
def test_symmetric_double_end_detuning_preserves_ten_peripheral_directions(r):
    epsilon = r - 1.0
    exact_epsilon = sympy.Rational(-1, 10) if r < 1 else sympy.Rational(1, 10)
    certificate = mprs.exact_control_subspaces(
        exact_epsilon,
        sympy.Rational(3, 10),
        both_ends=True,
    )
    assert certificate["blind_dimension"] == 3
    assert mprs.float_peripheral_count(epsilon, 0.3, both_ends=True) == 10
    result = mprs.direct_float_gaps(epsilon, 0.3, both_ends=True)
    assert result["stationary_dimension"] == 4
    assert result["peripheral_dimension"] == 10
    assert result["nonstationary_peripheral_dimension"] == 6
    assert result["a_gap"] == pytest.approx(0.0, abs=2e-11)
    assert result["next_distinct_a_rate"] > 1e-3


def test_moving_watch_to_seat_two_opens_a_finite_gap():
    result = mprs.direct_float_gaps(0.0, 0.3, seat=2)
    assert result["stationary_dimension"] == 1
    assert result["peripheral_dimension"] == 1
    assert result["a_gap"] > 1e-3


def test_exact_control_certificates_reject_float_scalars_at_the_exact_door():
    with pytest.raises(TypeError):
        mprs.exact_control_subspaces(0.1, sympy.Rational(3, 10))
    with pytest.raises(TypeError):
        mprs.exact_control_subspaces(sympy.Rational(1, 10), 0.3)


def test_quadratic_defect_path_has_fourth_order_gap():
    values = [mprs.tracked_reduced_gap(e * e, 0.3) for e in (2**-5, 2**-6, 2**-7)]
    exponents = [math.log(values[i] / values[i + 1], 2) for i in range(2)]
    assert min(exponents) > 3.8


@pytest.mark.parametrize("epsilon", [0.125, -0.125, 0.3])
def test_site_zero_phase_gauge_is_epsilon_to_minus_two_minus_epsilon(epsilon):
    assert mprs.tracked_reduced_gap(epsilon, 0.3) == pytest.approx(
        mprs.tracked_reduced_gap(-2.0 - epsilon, 0.3),
        rel=2e-10,
        abs=2e-12,
    )
    assert abs(
        mprs.tracked_reduced_gap(epsilon, 0.3)
        - mprs.tracked_reduced_gap(-epsilon, 0.3)
    ) > 1e-8


@pytest.mark.parametrize("epsilon", [0.1, -0.1])
def test_new_a_builder_matches_historical_run_la_without_changing_it(epsilon):
    from simulations.missing_phase_long_time import Run

    old = Run(epsilon, 0.3).la
    new = mprs.a_generator_float(epsilon, 0.3)
    assert np.linalg.norm(old - new, ord=2) < 1e-14


def test_global_flip_carrier_copies_are_identical():
    blocks = mprs.direct_carrier_blocks(
        sympy.Rational(1, 8), sympy.Rational(3, 10)
    )
    assert blocks["A_11"] == blocks["A_66"]
    assert blocks["B_16"] == blocks["B_61"]
    assert blocks["A_11"] is not blocks["A_66"]
    assert blocks["B_16"] is not blocks["B_61"]


def test_legacy_real_part_sorted_index_is_not_a_branch_at_epsilon_zero():
    root2 = float(sympy.sqrt(2))
    values = np.array(
        [
            0j,
            0j,
            0j,
            0j,
            -2j * root2,
            -2j * root2,
            +2j * root2,
            +2j * root2,
            -4j * root2,
            +4j * root2,
        ]
    )
    pick = lambda xs: xs[np.argsort(-xs.real, kind="stable")[2]]
    seen = {
        pick(values[permutation])
        for permutation in (
            np.arange(10),
            np.arange(10)[::-1],
            np.roll(np.arange(10), 4),
        )
    }
    assert len(seen) > 1
    census = mprs.exact_uniform_peripheral_census(sympy.Rational(3, 10))
    assert census["frequency_multiplicities"]["+2sqrt2i"] == 2


def test_b_adjoint_mutation_breaks_the_entrywise_operator_map():
    result = mprs.b_operator_mutation_residuals(
        sympy.Rational(1, 8), sympy.Rational(3, 10)
    )
    assert result["correct"] == 0
    assert result["omitted_adjoint"] != 0


def test_b_missing_right_action_and_price_mutations_break_the_full_spectrum_map():
    result = mprs.b_spectrum_mutation_residuals(
        sympy.Rational(1, 8), sympy.Rational(3, 10)
    )
    assert result["correct"] == 0
    assert result["missing_right_action"] != 0
    assert result["wrong_price_sign"] != 0
