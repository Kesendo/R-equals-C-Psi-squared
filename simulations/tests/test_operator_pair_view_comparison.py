"""Signal ranks in the (1,1) block: exact maps, the two laws, and the parity motion."""

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
    blind_seat_reach,
    centre_return_recurrence,
    centre_return_signature,
    channel_order,
    channel_prediction,
    coupling_scan,
    cramer_return_rank,
    exact_hankel_rank,
    exact_rank_map,
    f157_blind,
    frozen_root_kernel_dims,
    general_bound,
    integer_hermitian_generator,
    middle_seat_rate_scan,
    node_modes,
    odd_weight_curve,
    odd_weight_derivatives,
    omega2_dimension,
    population_moments_commute,
    r90_locus_prediction,
    single_rate_omega2_dimension,
    single_seat_prediction,
    uniform_prediction,
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


def _full(n: int, value: int) -> tuple[tuple[int, ...], ...]:
    return tuple((value,) * n for _ in range(n))


def _cross(n: int, seats: set[int], on: int, off: int) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(on if a in seats or b in seats else off for b in range(n)) for a in range(n))


# ---------------------------------------------------------------- the generator

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
                pair_state.reshape(n, n).diagonal().real, integer_state[:n], rtol=0, atol=1e-9
            )
            pair_state = pair_generator @ pair_state
            integer_state = integer_generator @ integer_state


# ---------------------------------------------------------------- exact maps

@pytest.mark.parametrize(
    "gamma, expected",
    [
        ((0, 1, 0), _cross(3, {1}, 4, 8)),
        ((1, 1, 1), _cross(3, {1}, 4, 8)),
        ((1, 2, 1), _cross(3, {1}, 5, 9)),
        ((1, 1, 0), _full(3, 9)),
    ],
)
def test_n3_exact_rank_maps(gamma, expected):
    assert exact_rank_map(3, gamma) == expected


def test_the_remaining_page_maps_up_to_n7():
    assert exact_rank_map(3, (0, 3, 0)) == single_seat_prediction(3, 1) == _cross(3, {1}, 4, 8)
    assert exact_rank_map(5, (0, 1, 2, 3, 4)) == r90_locus_prediction(5, (0, 1, 2, 3, 4)) == _full(5, 23)
    assert exact_rank_map(7, (0, 0, 0, 1, 0, 0, 0)) == single_seat_prediction(7, 3) == _cross(7, {3}, 16, 46)
    assert exact_rank_map(7, (1,) * 7) == channel_prediction(7) == uniform_prediction(7) == _cross(7, {3}, 22, 46)


def test_the_n3_coincidence_is_two_laws_that_part_at_n5():
    # N = 3: centre-only (F157) and uniform (Omega_2) give the same map ...
    assert exact_rank_map(3, (0, 1, 0)) == exact_rank_map(3, (1, 1, 1))
    # ... and a symmetric profile with unequal rates keeps neither reduction
    assert exact_rank_map(3, (1, 2, 1)) == _cross(3, {1}, 5, 9)
    # N = 5: they agree off the centre (both N^2 - m, m = 2) and part on the cross by m(m-1) = 2
    centre, uniform = exact_rank_map(5, (0, 0, 1, 0, 0)), exact_rank_map(5, (1, 1, 1, 1, 1))
    assert centre == _cross(5, {2}, 9, 23)
    assert uniform == _cross(5, {2}, 11, 23)


# ---------------------------------------------------------------- law 1: one seat

@pytest.mark.parametrize("n, seat", [(3, 1), (5, 2), (5, 0), (5, 1), (6, 1), (6, 0)])
def test_single_seat_law_meets_the_exact_map(n: int, seat: int):
    gamma = tuple(1 if k == seat else 0 for k in range(n))
    assert exact_rank_map(n, gamma) == single_seat_prediction(n, seat)


def test_the_one_seat_bound_is_not_rate_free():
    # at 2J the curves through a middle seat of an even chain read N^2 - (N/2)^2,
    # every eigenvalue of L staying simple: residues vanish, no eigenvalues merge
    assert exact_rank_map(4, (0, 2, 0, 0)) == _cross(4, {1}, 12, 16)
    assert exact_rank_map(6, (0, 0, 2, 0, 0, 0)) == _cross(6, {2}, 27, 36)
    for n in (4, 6, 8):
        for row in middle_seat_rate_scan(n):
            assert row["one_seat_bound"] == n * n
            assert row["distinct_eigenvalues"] == row["dimension"] == n * n
            expected = n * n - (n // 2) ** 2 if row["rate"] == 2 else n * n
            assert row["reachable_dimension"] == row["return_rank"] == expected
            if n <= 6:
                assert row["map_meets_bound"] == (row["rate"] != 2)


def test_single_seat_law_needs_the_shared_nodes():
    # N = 6, light on seat 1: seat 4 shares the blind mode's node, so its row is capped too.
    exact = exact_rank_map(6, (0, 1, 0, 0, 0, 0))
    assert exact == _cross(6, {1, 4}, 25, 35)
    naive = _cross(6, {1}, 25, 35)            # the seat's own row only: must NOT match
    assert exact != naive


def test_single_seat_law_needs_the_kernel_inside_the_reachable_space():
    # N = 9, light on the centre: seats 1 and 7 share mode 3 with it.  Inside End(K_1) only
    # blind - s = 4 - 1 = 3 of the centre's blind projectors remain, so the row reads
    # (9 - 1)^2 - 3 = 61, not the 64 a form without that kernel term gives.
    n, seat = 9, 4
    generator = integer_hermitian_generator(n, tuple(1 if k == seat else 0 for k in range(n)))
    law = single_seat_prediction(n, seat)
    for a, b in ((1, 1), (1, 0), (0, 0), (4, 4), (0, 7)):
        assert exact_hankel_rank(generator, a, b) == law[a][b]
    assert law[1][1] == 61 and law[0][0] == 77 and law[4][4] == 25
    blind = f157_blind(n, seat)
    reach = [n - len(node_modes(n, a) & node_modes(n, seat)) for a in range(n)]
    without_kernel_term = min(reach[1] ** 2, n * n - blind)
    assert exact_hankel_rank(generator, 1, 1) != without_kernel_term    # 64: must NOT match


# ---------------------------------------------------------------- law 2: uniform light

@pytest.mark.parametrize("n", [3, 4, 5, 6, 7, 8])
def test_omega2_dimension_is_floor_half(n: int):
    assert omega2_dimension(n) == n // 2


@pytest.mark.parametrize("n", [3, 4, 5])
def test_uniform_law_meets_the_exact_map(n: int):
    exact = exact_rank_map(n, (1,) * n)
    assert exact == uniform_prediction(n)
    # the same reflection form without the floor(N/2) invisible modes must fail
    reflection_dim = ((n + 1) // 2) ** 2 + (n // 2) ** 2
    without_omega2 = _cross(n, {n // 2}, reflection_dim, n * n) if n % 2 else _full(n, n * n)
    assert exact != without_omega2


def test_uniform_n6_falls_below_the_bound_only_at_the_blind_seats():
    exact, bound = exact_rank_map(6, (1,) * 6), uniform_prediction(6)
    below = {(a, b) for a in range(6) for b in range(6) if exact[a][b] < bound[a][b]}
    above = {(a, b) for a in range(6) for b in range(6) if exact[a][b] > bound[a][b]}
    assert not above
    assert below == {(a, b) for a in range(6) for b in range(6) if a in (1, 4) or b in (1, 4)}
    assert exact == _cross(6, {1, 4}, 26, 33)


def test_uniform_reach_at_the_blind_seats():
    rows = {r["seat"]: r for r in blind_seat_reach(6)}
    assert (rows[0]["reachable_dimension"], rows[1]["reachable_dimension"]) == (33, 26)
    assert rows[1]["f157_blind"] == 1 and rows[1]["below_bound_by"] == 1 * (6 + 1)
    ten = {r["seat"]: r for r in blind_seat_reach(10)}   # a seat with blindness two
    assert ten[2]["f157_blind"] == 2 and ten[2]["below_bound_by"] == 2 * (10 + 1)
    control = blind_seat_reach(8)               # N = 8 has no blind seat: reach meets the bound
    assert all(r["below_bound_by"] == 0 for r in control)
    assert all(r["meets_channel_sum"] for n in (6, 8, 10) for r in blind_seat_reach(n))


# ---------------------------------------------------------------- uniform light, channel by channel

@pytest.mark.parametrize("n", [3, 4, 5, 6])
def test_population_transfer_is_diagonal_in_the_heisenberg_eigenbasis(n: int):
    assert population_moments_commute(n, 1) is None
    assert population_moments_commute(n, 2) is None
    # the XY chain's sine modes have no complementary pairing: the second moment already fails
    assert population_moments_commute(n, 1, "xy") == 2


@pytest.mark.parametrize("n", range(2, 14))
def test_channel_orders_sum_to_the_uniform_bound(n: int):
    assert sum(channel_order(n, q) for q in range(n)) == n * n - n // 2


@pytest.mark.parametrize("n", range(2, 31))
def test_every_blind_mode_is_in_the_n_plus_one_class(n: int):
    for seat in range(n):
        modes = node_modes(n, seat)
        assert len(modes) == f157_blind(n, seat)
        assert all((q - n) % 2 == 0 and channel_order(n, q) == n + 1 for q in modes)


@pytest.mark.parametrize("n, rate", [(3, 1), (4, 1), (4, 3), (5, 1), (5, 2), (6, 1), (6, 3)])
def test_channel_sum_is_the_uniform_map(n: int, rate: int):
    assert exact_rank_map(n, (rate,) * n) == channel_prediction(n)
    # the Omega_2 bound alone is too high at N = 6, where seats 1 and 4 are blind
    if n == 6:
        assert channel_prediction(n) != uniform_prediction(n)


# ---------------------------------------------------------------- the R90 locus and the upper bound

@pytest.mark.parametrize("n, gamma", [(3, (0, 1, 2)), (4, (0, 1, 2, 3))])
def test_frozen_divisor_hides_floor_half_modes_on_the_r90_locus(n, gamma):
    exact = exact_rank_map(n, gamma)
    assert exact == r90_locus_prediction(n, gamma) == _full(n, n * n - n // 2)
    # the bound without F140's frozen modes is too high there: must NOT match
    assert exact != general_bound(n, gamma)


def test_upper_bound_is_met_on_the_page_profiles():
    assert single_rate_omega2_dimension(4, (0, 1, 1, 0)) == 1
    assert exact_rank_map(4, (0, 1, 1, 0)) == general_bound(4, (0, 1, 1, 0)) == _full(4, 15)
    assert exact_rank_map(4, (1, 1, 1, 0)) == general_bound(4, (1, 1, 1, 0)) == _full(4, 16)
    assert exact_rank_map(3, (1, 2, 1)) == general_bound(3, (1, 2, 1))
    assert exact_rank_map(3, (1, 1, 0)) == general_bound(3, (1, 1, 0)) == _full(3, 9)
    assert exact_rank_map(5, (1, 0, 1, 0, 0)) == general_bound(5, (1, 0, 1, 0, 0)) == _full(5, 25)


def test_upper_bound_misses_a_mode_that_hides_from_one_site_only():
    # N = 5, gamma = (0,1,0,1,0): 13 on the cross, as bounded; 25, the bound, where
    # preparation and readout both sit on seats 1 or 3; 23 on every other curve.
    exact, bound = exact_rank_map(5, (0, 1, 0, 1, 0)), general_bound(5, (0, 1, 0, 1, 0))
    expected = tuple(tuple(13 if 2 in (a, b) else 25 if {a, b} <= {1, 3} else 23
                           for b in range(5)) for a in range(5))
    assert exact == expected
    assert bound == _cross(5, {2}, 13, 25)


def test_the_locus_misses_sit_at_exceptional_couplings():
    # at J = 1 the frozen root is defective: kernel dims of (L + 4 gamma-bar)^k
    assert frozen_root_kernel_dims(4, (0, 2, 0, 2)) == (2, 3, 4, 4)       # blocks 3 + 1
    assert frozen_root_kernel_dims(5, (1, 0, 1, 2, 1)) == (2, 3, 3, 3)    # blocks 2 + 1
    assert frozen_root_kernel_dims(4, (0, 1, 2, 3)) == (2, 2, 2, 2)       # semisimple, bound met
    # at J = 2 and 3 both profiles are semisimple and meet the R90 bound
    for coupling in (2, 3):
        assert frozen_root_kernel_dims(4, (0, 2, 0, 2), coupling) == (2, 2, 2, 2)
        assert frozen_root_kernel_dims(5, (1, 0, 1, 2, 1), coupling) == (2, 2, 2, 2)
        assert exact_rank_map(4, (0, 2, 0, 2), coupling=coupling) == _full(4, 14)
        assert exact_rank_map(5, (1, 0, 1, 2, 1), coupling=coupling) == _full(5, 23)


def test_the_end_sites_under_01010_meet_the_bound_off_j_equal_gamma():
    rows = {r["coupling"]: r for r in coupling_scan() if tuple(r["gamma"]) == (0, 1, 0, 1, 0)}
    bound = general_bound(5, (0, 1, 0, 1, 0))
    assert all(rows[j]["distinct_eigenvalues"] == 25 for j in (1, 2, 3))   # simple spectrum throughout
    assert tuple(tuple(r) for r in rows[1]["ranks"]) != bound
    assert tuple(tuple(r) for r in rows[2]["ranks"]) == bound
    assert tuple(tuple(r) for r in rows[3]["ranks"]) == bound


def test_the_r90_bound_is_not_tight_on_the_whole_locus():
    exact = exact_rank_map(4, (0, 2, 0, 2))
    assert r90_locus_prediction(4, (0, 2, 0, 2)) == _full(4, 14)
    assert exact == tuple(tuple(13 if {a, b} <= {0, 2} else 12 for b in range(4)) for a in range(4))
    assert exact_rank_map(5, (1, 0, 1, 2, 1)) == _full(5, 22)
    assert r90_locus_prediction(5, (1, 0, 1, 2, 1)) == _full(5, 23)


# ---------------------------------------------------------------- the parity motion

def test_n5_centre_jump_keeps_a_nonstationary_pure_odd_space():
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


def test_same_n3_rank_map_can_hide_different_parity_motion():
    assert odd_weight_derivatives((0, 1, 0), 4) == (0, 0, 0, 0)
    assert odd_weight_derivatives((0, 3, 0), 4) == (0, 0, 0, 0)
    assert odd_weight_derivatives((1, 1, 1), 4) == (0, 0, 0, 32)
    assert odd_weight_derivatives((1, 1, 0), 4) == (0, 0, 0, 16)
    assert odd_weight_curve((0, 1, 0), [0.5])[0] == pytest.approx(0, abs=1e-12)
    assert odd_weight_curve((1, 1, 1), [0.5])[0] > 0.18


@pytest.mark.parametrize("n, gamma", [(3, (0, 1, 0)), (3, (1, 1, 0)), (4, (1, 1, 1, 1)), (4, (0, 2, 0, 2)),
                                      (5, (0, 1, 0, 0, 0)), (5, (0, 1, 0, 1, 0))])
def test_the_return_rank_is_the_cramer_count(n, gamma):
    exact = exact_rank_map(n, gamma)
    assert all(cramer_return_rank(n, gamma, j) == exact[j][j] for j in range(n))


def test_equal_ranks_carry_different_recurrences():
    # central light: the Bloch cubic of the flow page; uniform light: The Reflection's p3 at gamma = J
    assert centre_return_recurrence((0, 1, 0)) == (1, 4, 40, 64, 0)
    assert centre_return_recurrence((1, 1, 1)) == (1, 8, 52, 96, 0)


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


def test_exact_rank_refuses_non_integer_rates():
    with pytest.raises(ValueError, match="integer"):
        exact_rank_map(3, (0, 0.5, 0))
