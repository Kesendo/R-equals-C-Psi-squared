"""Independent controls for the coherent F124 -> F157 bond/seat producer."""

from __future__ import annotations

import importlib
import sys
from math import gcd
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm, expm_frechet


sys.path.insert(0, str(Path(__file__).parent))


def producer():
    try:
        return importlib.import_module("handshake_bond_seat_readout")
    except ModuleNotFoundError as exc:
        pytest.fail(f"coherent bond/seat producer is missing: {exc}")


def test_n3_transition_matrix_has_literal_bond_channels():
    m = producer().transition_matrix(3)
    expected = np.array(
        [[np.sqrt(0.5), 0.5, 0.0], [np.sqrt(0.5), -0.5, 0.0]]
    )
    np.testing.assert_allclose(m, expected, rtol=0.0, atol=16 * np.finfo(float).eps)


def test_n7_center_and_offcenter_have_different_exact_mode_inventories():
    read = producer().seat_inventory
    assert read(7, 3) == {
        "h": 4,
        "node_modes": (2, 4, 6),
        "visible_modes": (3, 5),
        "rank": 2,
    }
    assert read(7, 2) == {
        "h": 1,
        "node_modes": (),
        "visible_modes": (2, 3, 4, 5, 6),
        "rank": 5,
    }


def direct_population_slope(n: int, seat: int, bond: int, time: float,
                            prepared_mode: int = 1) -> float:
    """Independent N-by-N Hamiltonian derivative, with no modal dictionary."""
    h = np.diag(np.ones(n - 1), 1) + np.diag(np.ones(n - 1), -1)
    v = np.zeros((n, n))
    v[bond, bond + 1] = v[bond + 1, bond] = 1.0
    x = np.arange(1, n + 1)
    initial = np.sqrt(2.0 / (n + 1)) * np.sin(
        np.pi * prepared_mode * x / (n + 1)
    )
    u, du = expm_frechet(-1j * time * h, -1j * time * v,
                        compute_expm=True)
    z = (u @ initial)[seat]
    dz = (du @ initial)[seat]
    return float(2.0 * np.real(np.conj(z) * dz))


def test_modal_response_meets_independent_hamiltonian_derivative():
    predict = producer().predicted_slopes
    cases = ((3, 0, 0, 0.7), (3, 1, 0, 1.0),
             (7, 3, 0, 1.0), (7, 2, 0, 1.0), (7, 2, 5, 1.0))
    nonzero = []
    for n, seat, bond, time in cases:
        predicted = float(predict(n, seat, [time])[bond, 0])
        direct = direct_population_slope(n, seat, bond, time)
        budget = (128 * np.finfo(float).eps * n *
                  (1.0 + abs(time) + abs(predicted) + abs(direct)))
        assert abs(predicted - direct) <= budget
        nonzero.append(abs(direct))
    assert max(nonzero) > 0.01  # a zero-returning path cannot pass this gate


def test_short_time_quadratic_response_survives_cosine_cancellation():
    module = producer()
    tiny, tenfold = module.predicted_slopes(7, 2, [1e-9, 1e-8])[0]
    direct = module.direct_population_slope(7, 2, 0, 1e-9)
    assert abs(direct) > 1e-20
    np.testing.assert_allclose(tiny, direct, rtol=1e-11, atol=0.0)
    assert abs(tenfold / tiny - 100.0) < 1e-10


def test_short_time_centre_response_retains_fourth_order_sign():
    module = producer()
    time = 1e-8
    predicted = float(module.predicted_slopes(7, 3, [time])[0, 0])
    direct = module.direct_population_slope(7, 3, 0, time)
    leading = -np.sqrt(2.0) * time**4 / 48.0
    np.testing.assert_allclose(direct, leading, rtol=1e-10, atol=0.0)
    np.testing.assert_allclose(predicted, leading, rtol=1e-7, atol=0.0)


def test_centre_fourth_order_coefficient_from_exact_sine_modes():
    import sympy as sp

    hopping = sp.symbols("J", positive=True)
    n, seat, bond = 7, 3, 0

    def mode(k, site):
        return sp.sqrt(sp.Rational(2, n + 1)) * sp.sin(
            sp.pi * k * (site + 1) / (n + 1)
        )

    def element(k):
        return (mode(k, bond) * mode(1, bond + 1) +
                mode(k, bond + 1) * mode(1, bond))

    def gap(k):
        return 2 * hopping * (sp.cos(sp.pi / (n + 1)) -
                              sp.cos(sp.pi * k / (n + 1)))

    quadratic = mode(1, seat) * sum(
        mode(k, seat) * element(k) * gap(k) for k in range(2, n)
    )
    fourth = -mode(1, seat) / 12 * sum(
        mode(k, seat) * element(k) * gap(k)**3 for k in range(2, n)
    )
    assert sp.simplify(quadratic) == 0
    assert sp.simplify(fourth + sp.sqrt(2) * hopping**3 / 48) == 0


def test_producer_gate_includes_the_short_time_control():
    short = [row for row in producer().response_gate()["readings"]
             if row["time"] == 1e-9]
    assert len(short) == 1
    assert abs(short[0]["predicted"]) > 1e-20
    assert short[0]["error_ratio"] < 16.0


def test_producer_gate_checks_the_centre_fourth_order_germ():
    report = producer().response_gate()
    fourth_order = [row for row in report["readings"]
                    if (row["n"], row["seat"], row["bond"], row["time"])
                    == (7, 3, 0, 1e-8)]
    assert len(fourth_order) == 1
    assert report["fourth_order_relative_error"] < 1e-6
    assert report["passes"]


def test_modal_rank_report_matches_gcd_across_n3_to_n11():
    report = producer().modal_rank_report
    for n in range(3, 12):
        for seat in range(n):
            r = report(n, seat)
            h = gcd(seat + 1, n + 1)
            assert r["observed_full_rank"] == n - h
            assert r["observed_location_rank"] == n - h - 1
            assert r["passes"]


def test_rank_gate_rejects_a_wrong_bond_transition_matrix(monkeypatch):
    module = producer()
    monkeypatch.setattr(module, "transition_matrix",
                        lambda n: np.ones((n - 1, n)))
    assert not module.modal_rank_report(7, 3)["passes"]


def test_direct_derivative_has_quadratic_central_difference_convergence():
    direct = producer().direct_population_slope(7, 2, 0, 1.0)
    n = 7
    h = np.diag(np.ones(n - 1), 1) + np.diag(np.ones(n - 1), -1)
    v = np.zeros((n, n))
    v[0, 1] = v[1, 0] = 1.0
    psi1 = np.sqrt(2.0 / 8.0) * np.sin(np.pi * np.arange(1, 8) / 8.0)

    def population(epsilon):
        return abs((expm(-1j * (h + epsilon * v)) @ psi1)[2]) ** 2

    def difference(epsilon):
        return (population(epsilon) - population(-epsilon)) / (2 * epsilon)

    error_1 = abs(difference(0.02) - direct)
    error_2 = abs(difference(0.01) - direct)
    assert error_1 > 1e-7  # nonzero truncation error, so the ratio is informative
    assert 3.8 < error_1 / error_2 < 4.2


def test_preparing_psi2_breaks_the_psi1_response_prediction():
    module = producer()
    predicted = module.predicted_slopes(7, 2, [1.0])[0, 0]
    wrong_preparation = module.direct_population_slope(
        7, 2, 0, 1.0, prepared_mode=2
    )
    assert abs(wrong_preparation - predicted) > 0.02


def test_response_gate_rejects_a_corrupted_transition_matrix(monkeypatch):
    module = producer()
    good = module.response_gate()
    assert good["passes"]
    assert good["max_direct_slope"] > 0.01
    monkeypatch.setattr(module, "transition_matrix",
                        lambda n: np.zeros((n - 1, n)))
    bad = module.response_gate()
    assert not bad["passes"]
    assert bad["max_error_ratio"] > 1e6


def test_response_gate_covers_every_n7_bond_at_two_seats_and_times():
    rows = producer().response_gate()["readings"]
    measured = {(row["seat"], row["bond"], row["time"])
                for row in rows if row["n"] == 7 and row["time"] in (1.0, 1.3)}
    expected = {(seat, bond, time)
                for seat in (2, 3) for bond in range(6)
                for time in (1.0, 1.3)}
    assert measured == expected


def test_response_gate_rejects_interior_bond_mode_mutation(monkeypatch):
    module = producer()
    original = module.transition_matrix

    def corrupted(n):
        m = original(n)
        if n == 7:
            m[1, 2] += 0.1
            m[4, 2] -= 0.1
        return m

    monkeypatch.setattr(module, "transition_matrix", corrupted)
    assert not module.response_gate()["passes"]


@pytest.mark.parametrize("bad_value", [np.nan, np.inf, -np.inf])
def test_response_gate_rejects_nonfinite_later_acquisition(monkeypatch, bad_value):
    module = producer()
    original = module.direct_population_slope

    def corrupted(n, seat, bond, time, j_coupling=1.0, prepared_mode=1):
        if n == 9:
            return bad_value
        return original(n, seat, bond, time, j_coupling, prepared_mode)

    monkeypatch.setattr(module, "direct_population_slope", corrupted)
    assert not module.response_gate()["passes"]


def test_response_gate_checks_inverse_j_scaling_at_fixed_jt():
    rows = [row for row in producer().response_gate()["readings"]
            if row["J"] in (0.1, 10.0)]
    assert len(rows) == 2
    scaled = [row["J"] * row["direct"] for row in rows]
    assert abs(scaled[0] - scaled[1]) < 64 * np.finfo(float).eps
    assert all(row["error_ratio"] < 16.0 for row in rows)


def test_response_budget_keeps_the_zero_channel_at_weak_j():
    report = producer().response_gate()
    weak = [row for row in report["readings"]
            if row["J"] in (0.001, 0.0001)]
    assert len(weak) == 2
    assert all(row["n"] == 3 and row["seat"] == 1 for row in weak)
    assert report["passes"]
    assert all(row["error_ratio"] < 16.0 for row in weak)


def test_finite_defect_mirror_gate_distinguishes_a_wrong_pair():
    report = producer().reflection_report
    good = report()
    assert good["passes"]
    assert good["centre_difference"] < good["roundoff_budget"]
    assert good["offcentre_difference"] > 0.01
    assert not report(paired_bond=4)["passes"]


def test_producer_prints_its_scope_and_n7_findings(capsys):
    assert producer().main() == 0
    output = capsys.readouterr().out
    assert "gamma=0" in output
    assert "H_SE=J_hop*A_path" in output
    assert "24 N=7 all-bond comparisons" in output
    assert "J=1, t=1" in output
    assert "t=1e-8, J=1" in output
    assert "leading=-2.946278e-34" in output
    assert "N=7 location ranks: [5, 4, 5, 2, 5, 4, 5]" in output
    assert "centre modes: (3, 5)" in output
    assert "PTF alpha" in output
    assert "VERDICT: PASS" in output
