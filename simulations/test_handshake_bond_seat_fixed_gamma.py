"""Controls for a fixed positive dephasing rate in the bond-to-seat workflow."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp


sys.path.insert(0, str(Path(__file__).parent))


def producer():
    try:
        return importlib.import_module("handshake_bond_seat_fixed_gamma")
    except ModuleNotFoundError as exc:
        pytest.fail(f"fixed-gamma producer is missing: {exc}")


def test_fixed_positive_gamma_generator_matches_existing_uniform_block():
    module = producer()
    from coherence_horizon_se_block import L_se

    np.testing.assert_allclose(
        module.liouvillian(5, 2.0, 1.0, tuple(range(5))),
        L_se(5, 2.0, 1.0), rtol=0.0, atol=0.0
    )
    single = module.liouvillian(5, 2.0, 1.0, (2,))
    damping = np.real(np.diag(single)).reshape(5, 5)
    assert damping[2, 0] == -2.0
    assert damping[0, 1] == 0.0
    assert damping[2, 2] == 0.0
    with pytest.raises(ValueError):
        module.liouvillian(5, 2.0, 0.0, tuple(range(5)))


def test_fixed_gamma_rank_fixture_differs_from_coherent_node_count():
    module = producer()
    report = module.fixed_gamma_gate()
    assert report["passes"]
    assert report["uniform_ranks"] == [6, 6, 6, 3, 6, 6, 6]
    assert report["single_centre_ranks"] == [6, 6, 6, 3, 6, 6, 6]
    assert report["coherent_ranks"] == [5, 4, 5, 2, 5, 4, 5]
    assert report["coherent_sampled_ranks"] == report["coherent_ranks"]
    assert sorted(report["q_rank_profiles"]) == [1.0, 2.0, 10.0]
    assert report["q_time_grids"][1.0][0] == pytest.approx(0.2)
    assert report["q_time_grids"][1.0][-1] == pytest.approx(10.0)
    assert report["q_time_grids"][2.0][0] == pytest.approx(0.1)
    assert report["q_time_grids"][10.0][0] == pytest.approx(0.02)
    assert all(rows[profile] == [6, 6, 6, 3, 6, 6, 6]
               for rows in report["q_rank_profiles"].values()
               for profile in ("uniform", "single_centre"))
    assert report["rank_tolerance_margin"] > 1e8


def test_fixed_gamma_derivative_agrees_with_independent_tangent_ode():
    module = producer()
    report = module.fixed_gamma_gate()
    assert len(report["tangent_checks"]) >= 4
    assert report["max_tangent_error_ratio"] < 16.0
    assert any(row["profile"] == "uniform" for row in report["tangent_checks"])
    assert any(row["profile"] == "single_centre" for row in report["tangent_checks"])
    assert max(abs(row["direct"]) for row in report["tangent_checks"]) > 0.005


def test_fixed_gamma_breaks_uniform_bond_null_and_keeps_reflection():
    report = producer().fixed_gamma_gate()
    assert abs(report["uniform_bond_sum"]) > 1e-4
    assert abs(report["single_centre_bond_sum"]) > 1e-4
    assert report["centre_mirror_residual"] < report["reflection_budget"]
    assert report["offcentre_mirror_difference"] > 0.001
    assert report["exact_generator_reflection"]


def test_independent_tangent_gate_rejects_wrong_dephasing_rate(monkeypatch):
    module = producer()
    original = module.liouvillian

    def half_rate(n, hopping, gamma0, watched_sites):
        return original(n, hopping, gamma0 / 2.0, watched_sites)

    monkeypatch.setattr(module, "liouvillian", half_rate)
    bad = module.fixed_gamma_gate()
    assert not bad["passes"]
    assert bad["max_tangent_error_ratio"] > 1e6


@pytest.mark.parametrize("watched_sites", [(), (99,), (2, 2)])
def test_tangent_route_rejects_invalid_channel_profiles(watched_sites):
    with pytest.raises(ValueError):
        producer().independent_tangent_slope(
            5, 2, 0, 0.1, 2.0, 1.0, watched_sites
        )


def test_fixed_gamma_producer_prints_operating_regime(capsys):
    assert producer().main(("--unitized-control",)) == 0
    output = capsys.readouterr().out
    assert "gamma0=1 fixed" in output
    assert "J=10, Q=J/gamma0=10" in output
    assert "uniform ranks: [6, 6, 6, 3, 6, 6, 6]" in output
    assert "single-centre ranks: [6, 6, 6, 3, 6, 6, 6]" in output
    assert "Q=1 uniform ranks: [6, 6, 6, 3, 6, 6, 6]" in output
    assert "Q=2 single-centre ranks: [6, 6, 6, 3, 6, 6, 6]" in output
    assert "Q=1 gamma0*t:" in output
    assert "Q=2 gamma0*t:" in output
    assert "Q=10 gamma0*t:" in output
    assert "coherent-limit ranks: [5, 4, 5, 2, 5, 4, 5]" in output
    assert "VERDICT: PASS" in output


def test_repo_gamma005_j0075_fixture_uses_the_carrier_hopping_book():
    report = producer().repo_gamma005_fixture()
    assert report["gamma0"] == pytest.approx(0.05)
    assert report["hopping"] == pytest.approx(0.075)
    assert report["q"] == pytest.approx(1.5)
    assert report["times"][0] == pytest.approx(2.0)
    assert report["times"][-1] == pytest.approx(30.0)
    assert report["refined_sample_count"] == 37
    for profile in ("uniform", "single_centre"):
        assert report["ranks"][profile] == [6, 6, 6, 3, 6, 6, 6]
        assert report["refined_ranks"][profile] == [6, 6, 6, 3, 6, 6, 6]
    assert report["coherent_ranks"] == [5, 4, 5, 2, 5, 4, 5]
    assert report["max_tangent_error_ratio"] < 16.0
    assert report["centre_mirror_residual"] < report["reflection_budget"]
    assert report["passes"]
    assert report["uniform_bond_sum_t10"] == pytest.approx(-0.192335351, abs=1e-8)
    assert report["offcentre_end_bond_difference_t10"]["uniform"] == pytest.approx(
        0.184782141, abs=1e-8
    )


def test_repo_fixture_detects_wrong_dephasing_rate(monkeypatch):
    module = producer()
    original = module.liouvillian

    def half_rate(n, hopping, gamma0, watched_sites):
        return original(n, hopping, gamma0 / 2.0, watched_sites)

    monkeypatch.setattr(module, "liouvillian", half_rate)
    report = module.repo_gamma005_fixture()
    assert not report["passes"]
    assert report["max_tangent_error_ratio"] > 1e6


def test_default_fixed_gamma_run_names_repo_point(capsys):
    assert producer().main() == 0
    output = capsys.readouterr().out
    assert "gamma0=0.05 fixed; J_hop=0.075; Q=1.5" in output
    assert "H=(J_hop/2)*sum(XX+YY)" in output
    assert "t=[2, 4, 6, 8, 10, 12, 15, 18, 24, 30]" in output
    assert "end-bond slope separation 0.184782 (uniform illumination)" in output
    assert "VERDICT: PASS" in output


def test_repo_point_all_slopes_match_explicit_site_z_lindblad_ode():
    """Independent physical Z operators and matrix ODE, across both profiles."""
    n, hopping, gamma0 = 7, 0.075, 0.05
    times = np.array([2., 4., 6., 8., 10., 12., 15., 18., 24., 30.])
    h = hopping * (np.diag(np.ones(n - 1), 1) +
                   np.diag(np.ones(n - 1), -1))
    bonds = []
    for bond in range(n - 1):
        direction = np.zeros((n, n))
        direction[bond, bond + 1] = direction[bond + 1, bond] = 1.0
        bonds.append(direction)
    psi = np.sqrt(2.0 / (n + 1)) * np.sin(
        np.pi * np.arange(1, n + 1) / (n + 1)
    )
    initial = np.zeros((n, n, n), dtype=complex)
    initial[0] = np.outer(psi, psi)

    for watched_sites in (tuple(range(n)), (3,)):
        site_z = []
        for site in watched_sites:
            signs = np.ones(n)
            signs[site] = -1.0
            site_z.append(np.diag(signs))

        def rhs(_, flat):
            states = flat.reshape(n, n, n)
            deriv = np.empty_like(states)
            rho = states[0]
            for index, state in enumerate(states):
                out = -1j * (h @ state - state @ h)
                for z in site_z:
                    out += gamma0 * (z @ state @ z - state)
                if index:
                    v = bonds[index - 1]
                    out += -1j * (v @ rho - rho @ v)
                deriv[index] = out
            return deriv.ravel()

        rtol, atol = 2e-12, 2e-14
        solved = solve_ivp(rhs, (0.0, times[-1]), initial.ravel(),
                           t_eval=times, method="DOP853", rtol=rtol, atol=atol)
        assert solved.success
        states = solved.y.T.reshape(len(times), n, n, n)
        direct = np.array([
            np.diagonal(row[1:], axis1=1, axis2=2).real for row in states
        ])
        frechet = producer().response_tensor(n, times, hopping, gamma0,
                                             watched_sites)
        integration_budget = 8.0 * (atol + rtol * np.max(np.abs(frechet)))
        assert np.max(np.abs(direct - frechet)) <= integration_budget
        assert producer()._sample_ranks(direct)[0] == [6, 6, 6, 3, 6, 6, 6]


def test_repo_point_unit_conversion_keeps_q15_and_scales_additive_slopes():
    module = producer()
    times = np.array([2.0, 10.0, 30.0])
    for watched_sites in (tuple(range(7)), (3,)):
        physical = module.response_tensor(7, times, 0.075, 0.05,
                                          watched_sites)
        unitized = module.response_tensor(7, 0.05 * times, 1.5, 1.0,
                                          watched_sites)
        residual = np.max(np.abs(physical - 20.0 * unitized))
        roundoff_budget = (64.0 * np.finfo(float).eps *
                           max(1.0, float(np.max(np.abs(physical)))))
        assert residual <= roundoff_budget


def test_repo_pass_rejects_a_retained_direction_near_roundoff(monkeypatch):
    module = producer()
    original = module.response_tensor

    def weak_direction(n, times, hopping, gamma0, watched_sites):
        tensor = original(n, times, hopping, gamma0, watched_sites)
        if (n, hopping, gamma0, watched_sites, len(times)) == (
            7, 0.075, 0.05, tuple(range(7)), 10
        ):
            u, singular, vh = np.linalg.svd(tensor[:, :, 1],
                                             full_matrices=False)
            tolerance = (64.0 * np.finfo(float).eps *
                         max(tensor.shape[:2]) * singular[0])
            singular[-1] = 1000.0 * tolerance
            tensor[:, :, 1] = (u * singular) @ vh
        return tensor

    monkeypatch.setattr(module, "response_tensor", weak_direction)
    report = module.repo_gamma005_fixture()
    assert report["ranks"]["uniform"] == [6, 6, 6, 3, 6, 6, 6]
    assert report["rank_tolerance_margin"] < 2000.0
    assert not report["passes"]


def test_repo_pass_rejects_erased_end_bond_separation(monkeypatch):
    module = producer()
    original = module.response_tensor

    def erase_difference(n, times, hopping, gamma0, watched_sites):
        tensor = original(n, times, hopping, gamma0, watched_sites)
        if (n, hopping, gamma0, watched_sites, len(times)) == (
            7, 0.075, 0.05, tuple(range(7)), 10
        ):
            time_index = list(times).index(10.0)
            change = tensor[time_index, 0, 2] - tensor[time_index, 5, 2]
            tensor[time_index, 5, 2] += change
            tensor[time_index, 4, 2] -= change
        return tensor

    monkeypatch.setattr(module, "response_tensor", erase_difference)
    report = module.repo_gamma005_fixture()
    assert report["ranks"]["uniform"] == [6, 6, 6, 3, 6, 6, 6]
    assert report["offcentre_end_bond_difference_t10"]["uniform"] == 0.0
    assert report["uniform_bond_sum_t10"] == pytest.approx(-0.192335351,
                                                            abs=1e-8)
    assert not report["passes"]
