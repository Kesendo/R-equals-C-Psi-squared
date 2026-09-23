"""Controls for a fixed positive dephasing rate in the bond-to-seat workflow."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import pytest


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
    assert producer().main() == 0
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
