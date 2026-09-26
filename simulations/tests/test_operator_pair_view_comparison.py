"""Finite view comparison: exact signal ranks and Hilbert-parity motion."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
import sympy as sp
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from framework.chain_system import ChainSystem
from framework.lindblad import lindbladian_z_dephasing
from operator_pair_flow_atlas import exact_heisenberg_se_h, one_excitation_basis, single_excitation_flow
from operator_pair_view_comparison import (
    certified_rank_map,
    centre_return_signature,
    integer_hermitian_generator,
    odd_weight_curve,
    odd_weight_derivatives,
)


def _coordinates(rho: np.ndarray) -> np.ndarray:
    n = rho.shape[0]
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    return np.array(
        [rho[a, a].real for a in range(n)]
        + [rho[a, b].real for a, b in pairs]
        + [rho[a, b].imag for a, b in pairs],
        dtype=float,
    )


@pytest.mark.parametrize("n", [3, 4, 5])
def test_integer_hermitian_generator_matches_independent_complex_pair_action(n: int):
    gamma = tuple(1 if site in (0, n // 2) else 0 for site in range(n))
    h = ChainSystem(n, J=1.0, H_type="heisenberg").H
    pair = single_excitation_flow(h, gamma)
    integer_generator = integer_hermitian_generator(n, gamma)
    amplitudes = np.array([1, 1j, 0.5 + 0.25j] + [0] * (n - 3), complex)
    amplitudes /= np.linalg.norm(amplitudes)
    rho = np.outer(amplitudes, amplitudes.conj())
    pair_rhs = (pair.generator @ rho.flatten()).reshape(n, n)
    np.testing.assert_allclose(
        integer_generator @ _coordinates(rho), _coordinates(pair_rhs), atol=1e-12
    )


def test_n4_population_moments_match_independent_complex_pair_flow():
    n = 4
    gamma = (1, 2, 0, 1)
    h = ChainSystem(n, J=1.0, H_type="heisenberg").H
    pair_generator = single_excitation_flow(h, gamma).generator
    integer_generator = integer_hermitian_generator(n, gamma)

    for preparation in range(n):
        pair_state = np.zeros(n * n, dtype=complex)
        pair_state[preparation * n + preparation] = 1
        integer_state = np.zeros(n * n, dtype=np.int64)
        integer_state[preparation] = 1
        for _order in range(7):
            np.testing.assert_allclose(
                pair_state.reshape(n, n).diagonal().real,
                integer_state[:n],
                rtol=0,
                atol=1e-9,
            )
            pair_state = pair_generator @ pair_state
            integer_state = integer_generator @ integer_state


@pytest.mark.parametrize(
    "gamma, expected",
    [
        ((0, 1, 0), ((8, 4, 8), (4, 4, 4), (8, 4, 8))),
        ((1, 1, 1), ((8, 4, 8), (4, 4, 4), (8, 4, 8))),
        ((1, 1, 0), ((9, 9, 9), (9, 9, 9), (9, 9, 9))),
    ],
)
def test_n3_rank_maps_are_exact_and_profile_sensitive(gamma, expected):
    result = certified_rank_map(3, gamma)
    assert result.ranks == expected
    assert result.modular_ranks[0] == result.modular_ranks[1] == expected
    assert result.method == "exact-sympy-and-modular"


def test_equal_total_rate_control_keeps_the_centre_strong_parity():
    centre = certified_rank_map(3, (0, 3, 0))
    uniform = certified_rank_map(3, (1, 1, 1))
    assert centre.ranks == uniform.ranks
    assert sum(centre.gamma) == sum(uniform.gamma) == 3
    assert odd_weight_derivatives(centre.gamma, 4) == (0, 0, 0, 0)
    assert odd_weight_derivatives(uniform.gamma, 4) == (0, 0, 0, 32)


@pytest.mark.parametrize(
    "gamma, rank",
    [((0, 1, 1, 0), 15), ((1, 1, 1, 0), 16)],
)
def test_n4_no_fixed_seat_control(gamma, rank):
    result = certified_rank_map(4, gamma)
    expected = tuple((rank,) * 4 for _ in range(4))
    assert result.ranks == expected
    assert result.modular_ranks[0] == result.modular_ranks[1] == expected


def test_n5_centre_rank_map_reaches_exact_symmetry_and_kernel_upper_bounds():
    result = certified_rank_map(5, (0, 0, 1, 0, 0))
    expected = tuple(
        tuple(9 if a == 2 or b == 2 else 23 for b in range(5))
        for a in range(5)
    )
    assert result.ranks == expected
    assert result.modular_ranks[0] == result.modular_ranks[1] == expected
    assert result.upper_bounds == expected
    assert result.generator_rank == 22
    assert result.method == "modular-lower-plus-exact-upper"


def test_n5_broken_profile_has_full_exact_scalar_rank():
    result = certified_rank_map(5, (1, 0, 1, 0, 0))
    expected = tuple((25,) * 5 for _ in range(5))
    assert result.ranks == expected
    assert result.modular_ranks[0] == result.modular_ranks[1] == expected
    assert result.upper_bounds == expected


def test_n5_centre_jump_keeps_a_nonstationary_pure_odd_block():
    odd = sp.Matrix.hstack(
        sp.Matrix([1, 0, 0, 0, -1]) / sp.sqrt(2),
        sp.Matrix([0, 1, 0, -1, 0]) / sp.sqrt(2),
    )
    h = exact_heisenberg_se_h(5, 1)
    centre_jump = sp.diag(1, 1, -1, 1, 1)
    edge_jump = sp.diag(-1, 1, 1, 1, 1)
    h_odd = odd.T * h * odd
    assert h_odd == sp.Matrix([[2, 2], [2, 0]])
    assert (sp.eye(5) - odd * odd.T) * h * odd == sp.zeros(5, 2)
    assert odd.T * centre_jump * odd == sp.eye(2)
    assert (sp.eye(5) - odd * odd.T) * centre_jump * odd == sp.zeros(5, 2)
    assert (sp.eye(5) - odd * odd.T) * edge_jump * odd != sp.zeros(5, 2)
    initial = sp.diag(1, 0)
    assert h_odd * initial != initial * h_odd
    assert sp.trace(initial * initial) == 1
    assert sp.trace((sp.eye(2) / 2) ** 2) == sp.Rational(1, 2)


def test_same_n3_rank_map_can_hide_different_hilbert_parity_motion():
    assert odd_weight_derivatives((0, 1, 0), 4) == (0, 0, 0, 0)
    assert odd_weight_derivatives((1, 1, 1), 4) == (0, 0, 0, 32)
    assert odd_weight_derivatives((1, 1, 0), 4) == (0, 0, 0, 16)
    assert odd_weight_curve((0, 1, 0), [0.5])[0] == pytest.approx(0, abs=1e-12)
    assert odd_weight_curve((1, 1, 1), [0.5])[0] > 0.18


@pytest.mark.parametrize("gamma, expected_third", [((0, 1, 0), 32), ((1, 1, 1), 64)])
def test_equal_rank_maps_do_not_imply_equal_centre_return_curves(gamma, expected_third):
    signature = centre_return_signature(gamma)
    assert signature["third_derivative_at_zero"] == expected_third

    h = ChainSystem(3, J=1.0, H_type="heisenberg").H
    full = lindbladian_z_dephasing(h, gamma)
    d = 1 << 3
    rho = np.zeros((d, d), complex)
    centre_state = one_excitation_basis(3)[1]
    rho[centre_state, centre_state] = 1
    evolved = (expm(0.5 * full) @ rho.flatten()).reshape(d, d)
    assert signature["population_at_t_0_5"] == pytest.approx(
        evolved[centre_state, centre_state].real, abs=1e-12
    )


@pytest.mark.parametrize("gamma", [(0, 1, 0), (1, 1, 1), (1, 1, 0)])
def test_parity_curve_agrees_with_full_density_evolution(gamma):
    n = 3
    h = ChainSystem(n, J=1.0, H_type="heisenberg").H
    full = lindbladian_z_dephasing(h, gamma)
    d = 1 << n
    rho = np.zeros((d, d), complex)
    centre_state = one_excitation_basis(n)[1]
    rho[centre_state, centre_state] = 1
    evolved = (expm(0.5 * full) @ rho.flatten()).reshape(d, d)
    sector = evolved[np.ix_(one_excitation_basis(n), one_excitation_basis(n))]
    reversal = np.eye(n)[::-1]
    odd_projector = (np.eye(n) - reversal) / 2
    expected = np.trace(odd_projector @ sector).real
    assert odd_weight_curve(gamma, [0.5])[0] == pytest.approx(expected, abs=1e-12)


def test_unknown_case_is_not_misreported_as_exact():
    with pytest.raises(ValueError, match="predeclared"):
        certified_rank_map(5, (1, 1, 1, 1, 1))
