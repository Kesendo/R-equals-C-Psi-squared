"""Exact regression checks for the N=11 (1,1) density compression."""

import pytest
import sympy as sp

from simulations import n11_compressed_density_gate as gate


@pytest.fixture(scope="module")
def data():
    modes = gate.sine_modes(11)
    dyads = gate.frequency_dyads(11, (1, 6))
    return modes, dyads


def test_complete_frequency_space_and_parities(data):
    _, dyads = data
    assert dyads == ((1, 6), (3, 7), (5, 9), (6, 11))
    assert tuple(gate.reflection_parity(d) for d in dyads) == (-1, 1, 1, -1)


def test_physical_mirror_contrast_has_exact_signed_matrix(data):
    modes, dyads = data
    observed = gate.compressed_indicator(modes, dyads, 0) - gate.compressed_indicator(modes, dyads, 10)
    expected = -sp.sqrt(2) / 72 * sp.Matrix(
        [[0, 1, 1, 0], [1, 0, 0, 1], [1, 0, 0, 1], [0, 1, 1, 0]]
    )
    assert gate.matrix_equal(observed, expected)
    assert sorted(observed.eigenvals().items(), key=lambda item: item[0]) == [
        (-sp.sqrt(2) / 36, 1), (0, 2), (sp.sqrt(2) / 36, 1)
    ]
    left = gate.compressed_indicator(modes, dyads, 0)
    assert left[2, 0] == -sp.sqrt(2) / 144  # opposite parity: nonzero
    assert left[2, 1] == -sp.Rational(1, 72) - sp.sqrt(3) / 144  # same parity: nonzero
    assert observed[2, 1] == 0  # the same-parity contrast is exactly zero


def test_balanced_profile_breaks_f154_identity_but_obeys_diagonal_overlap_form(data):
    modes, dyads = data
    rates = (2,) + (1,) * 9 + (0,)
    diss = gate.compressed_dissipator(modes, dyads, rates)
    size = sum((gate.compressed_indicator(modes, dyads, site) for site in range(11)),
               sp.zeros(4))
    diagonal_overlap = gate.compressed_diagonal_overlap(modes, dyads, rates)
    assert diss[2, 0] == sp.sqrt(2) / 36
    assert size[2, 0] == 0
    assert gate.matrix_equal(diss, -4 * sp.eye(4) + 4 * diagonal_overlap)
    assert gate.matrix_equal(diss + 2 * size, -2 * (
        gate.compressed_indicator(modes, dyads, 0)
        - gate.compressed_indicator(modes, dyads, 10)
    ))


def test_off_locus_rayleigh_value_escapes_balanced_interval(data):
    modes, dyads = data
    rates = (1,) + (0,) * 10
    diss = gate.compressed_dissipator(modes, dyads, rates)
    assert diss[2, 2] == -(22 + 5 * sp.sqrt(3)) / 72
    assert sp.simplify(diss[2, 2] + sp.Rational(4, 11)) < 0


def test_zero_frequency_block_attains_both_centres_from_physical_cells():
    endpoints = gate.verify_zero_frequency_room()
    assert endpoints == {"zero_frequency_upper": 0, "zero_frequency_lower": -4}
    reading = gate.verify_n11()
    assert reading["zero_frequency_upper"] == 0
    assert reading["zero_frequency_lower"] == -4
    assert reading["nonzero_room_spectrum"] == {
        -4: 2,
        -sp.Rational(10, 3) - sp.sqrt(2) / 18: 1,
        -sp.Rational(10, 3) + sp.sqrt(2) / 18: 1,
    }
    result = gate.verify_n11()
    assert result["off_locus_gamma_legs"] == (4 + sp.sqrt(3)) / 24


def test_mutated_cell_rule_is_rejected_through_same_gate():
    # Wrong OR occupancy replaces XOR disagreement in the physical cell action.
    def wrong_or(ket_site, bra_site, site):
        return int(ket_site == site or bra_site == site)

    with pytest.raises(AssertionError, match="contrast matrix"):
        gate.verify_n11(wrong_or)


@pytest.mark.parametrize("wrong_rule", [
    lambda ket_site, bra_site, site: int(ket_site == site or bra_site == site),
    lambda ket_site, bra_site, site: int(ket_site == site),
])
def test_mutated_cell_rule_fails_the_endpoint_half_on_its_own(wrong_rule):
    # A diagonal cell carries rate under either wrong rule, so D I is no longer zero.
    with pytest.raises(AssertionError, match="not an eigenvector|upper endpoint"):
        gate.verify_zero_frequency_room(wrong_rule)
