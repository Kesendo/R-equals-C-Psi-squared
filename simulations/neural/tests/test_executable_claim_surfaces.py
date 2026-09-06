"""Regression gates for executable claim surfaces repaired by the neural review."""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[3]


def _load(relative_path: str):
    path = ROOT / relative_path
    spec = importlib.util.spec_from_file_location(path.stem, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _run(relative_path: str) -> str:
    completed = subprocess.run(
        [sys.executable, str(ROOT / relative_path)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    return completed.stdout


def test_factor_two_census_distinguishes_the_two_involutions():
    module = _load("simulations/factor_two_standing_waves.py")
    values = module.load_eigenvalues(2)
    linear = module.orbit_census(values, 2 * module.GAMMA, "linear")
    composite = module.orbit_census(values, 2 * module.GAMMA, "composite")

    assert (linear.pairs, linear.fixed) == (6, 4)
    assert (composite.pairs, composite.fixed) == (3, 10)
    assert linear.assigned == composite.assigned == len(values)


def test_factor_two_report_pins_full_and_topology_censuses():
    module = _load("simulations/factor_two_standing_waves.py")
    output = module.render_report()
    assert "Committed CSV inputs: simulations/results/rmt_eigenvalues_N{2..7}.csv" in output
    assert "CSV metadata: columns Re, Im only" in output
    assert "dotnet run -c Release --project compute/RCPsiSquared.Compute -- rmt chain" in output
    assert "generator parameters in current source: chain, J=1.0, uniform gamma=0.05" in output
    assert "do not encode backend, source revision, command, or timestamp" in output
    assert "linear F1 total:       pairs=10903 fixed=34" in output
    assert "conjugate-composite:   pairs=9921 fixed=1998" in output
    assert "N=4 chain     linear: pairs=121 fixed=14" in output
    assert "N=4 star      linear: pairs=120 fixed=16" in output
    assert "N=4 ring      linear: pairs=116 fixed=24" in output
    assert "TOLERANCE STABILITY" in output
    assert "tol=1e-06" in output
    assert "tol=1e-08" in output
    assert "tol=1e-10" in output
    assert "mean decay=0.350000" in output
    for forbidden in ("finesse", "Beer-Lambert", "unpaired", "halve absorption"):
        assert forbidden.lower() not in output.lower()


def test_factor_two_committed_snapshot_matches_pure_render():
    module = _load("simulations/factor_two_standing_waves.py")
    expected = (ROOT / "simulations/results/factor_two_standing_waves.txt").read_text(encoding="utf-8")
    assert module.render_report() == expected


def test_general_demo_states_sufficiency_and_runs_negative_control():
    module = _load("simulations/palindrome_general.py")
    assert module.jordan_negative_control()[0]
    output = _run("simulations/palindrome_general.py")
    assert "sufficient operator identity" in output.lower()
    assert "not a converse" in output.lower()
    assert "dale constraints alone are insufficient" in output.lower()
    assert "jordan negative control" in output.lower()
    assert "quantum operator identity tested" in output.lower()
    assert "if and only if" not in output.lower()


def test_general_matcher_uses_global_multiplicity_assignment():
    module = _load("simulations/palindrome_general.py")
    # Greedy nearest-neighbour consumes 0.2 for the first target and then fails;
    # the valid assignment is 0.11 -> 0.0 and 0.2 -> 0.2.
    ok, worst = module.multiset_match(
        values=[0.0, 0.2], targets=[0.11, 0.2], tol=0.12
    )
    assert ok
    assert abs(worst - 0.11) < 1e-12


def test_general_matcher_reports_residual_in_row_to_column_orientation():
    module = _load("simulations/palindrome_general.py")
    # Each target has exactly one nearby value.  The resulting row->column
    # assignment is the non-self-inverse 3-cycle [1, 2, 0].
    ok, worst = module.multiset_match(
        values=[0.0, 10.0, 20.0],
        targets=[10.01, 20.02, 0.03],
        tol=0.1,
    )
    assert ok
    assert worst == pytest.approx(0.03)


def test_thermal_builder_includes_spontaneous_emission_at_zero():
    module = _load("simulations/thermal_blackbody.py")
    cold = module.build_thermal_liouvillian(1, [0.0], 0.0, [0.1])
    near = module.build_thermal_liouvillian(1, [0.0], 1e-9, [0.1])
    assert module.Sm[0, 1] == 1  # sigma- = |0><1|, emission
    assert module.Sp[1, 0] == 1  # sigma+ = |1><0|, absorption
    assert abs(cold).max() > 0.0
    assert abs(near - cold).max() < 1e-8

    excited = np.array([[0.0, 0.0], [0.0, 1.0]], dtype=complex).reshape(-1)
    flow = (cold @ excited).reshape(2, 2).real
    assert flow[0, 0] == pytest.approx(+0.1)
    assert flow[1, 1] == pytest.approx(-0.1)
    assert flow[0, 1] == flow[1, 0] == 0.0

    gate = getattr(module, "zero_temperature_emission_gate", lambda *_: {})
    signed_flow = gate(cold, 0.1)
    assert signed_flow["ground"] == pytest.approx(+0.1)
    assert signed_flow["excited"] == pytest.approx(-0.1)


@pytest.mark.parametrize("n_bar", [-1.0, np.nan, np.inf, -np.inf])
def test_thermal_builder_rejects_invalid_occupation(n_bar):
    module = _load("simulations/thermal_blackbody.py")
    with pytest.raises(ValueError):
        module.build_thermal_liouvillian(1, [0.0], n_bar, [0.1])


@pytest.mark.parametrize(
    "gammas,gamma_thermal",
    [
        ([-0.1], [0.1]),
        ([np.nan], [0.1]),
        ([np.inf], [0.1]),
        ([0.1], [-0.1]),
        ([0.1], [np.nan]),
        ([0.1], [np.inf]),
    ],
)
def test_thermal_builder_rejects_invalid_rates(gammas, gamma_thermal):
    module = _load("simulations/thermal_blackbody.py")
    with pytest.raises(ValueError):
        module.build_thermal_liouvillian(1, gammas, 0.0, gamma_thermal)


def test_repaired_standing_wave_report_uses_direct_observables():
    module = _load("simulations/standing_wave_analysis.py")
    output = module.render_report()
    assert "DIRECT PAULI-OBSERVABLE TIME TRACES" in output
    assert "half-range" in output
    for forbidden in ("osc%", "state weight in modes", "standing wave active"):
        assert forbidden.lower() not in output.lower()


def test_direct_pauli_trace_gate_has_hamiltonian_anchor_and_rejects_deleted_hamiltonian():
    module = _load("simulations/standing_wave_analysis.py")
    gate = getattr(module, "direct_trace_gate", lambda *_: {})
    anchors = gate()
    assert anchors
    assert anchors["w_iyy_t0"] == pytest.approx(2.0 / 3.0)
    assert anchors["w_iyy_dt0"] == pytest.approx(-2.0 / 15.0)
    assert anchors["w_state_dt0_norm"] > 0.1
    assert anchors["trace_dt0"] == pytest.approx(0.0, abs=1e-14)
    assert anchors["bell01_ixy_dt0"] == pytest.approx(-2.0)

    with pytest.raises(RuntimeError):
        gate(np.zeros_like(module.liouvillian()))

    original_j = module.J
    try:
        module.J = 0.0
        no_hamiltonian = module.liouvillian()
    finally:
        module.J = original_j
    with pytest.raises(RuntimeError):
        gate(no_hamiltonian)


def test_standing_wave_committed_snapshot_matches_pure_render():
    module = _load("simulations/standing_wave_analysis.py")
    expected = (ROOT / "simulations/results/standing_wave_analysis.txt").read_text(encoding="utf-8")
    assert module.render_report() == expected


def test_thermal_transition_surface_requires_defectiveness_gate():
    module = _load("simulations/thermal_ep_analysis.py")
    assert "not EP certificates" in (ROOT / "simulations/results/thermal_blackbody.txt").read_text(encoding="utf-8")
    assert module.oscillating_count([1 + 0j, 1 + 2j]) == 1


def test_thermal_committed_snapshot_matches_pure_render():
    module = _load("simulations/thermal_blackbody.py")
    expected = (ROOT / "simulations/results/thermal_blackbody.txt").read_text(encoding="utf-8")
    assert module.render_report() == expected


def test_veffect_control_runner_pins_refinement_and_transpose_outputs():
    module = _load("simulations/neural/veffect_controls.py")
    rows = module.coupling_controls()

    assert rows[(10, 0.05, 1e-6, False)] == (3, 12)
    assert rows[(10, 0.05, 1e-6, True)] == (3, 12)
    assert rows[(10, 0.05, 2.5e-7, False)] == (3, 12)
    assert rows[(10, 0.05, 6.25e-8, False)] == (3, 12)
    assert rows[(20, 0.05, 1e-6, False)] == (7, 62)
    assert rows[(20, 0.05, 2.5e-7, False)] == (8, 72)
    assert rows[(20, 0.05, 6.25e-8, False)] == (8, 73)
    assert rows[(20, 0.05, 1e-6, True)] == (7, 64)
    assert rows[(20, 0.0, 1e-6, False)] == (0, 0)
    assert rows[(20, 0.0, 1e-6, True)] == (2, 7)
    assert rows[(10, 0.0, 1e-6, True)] == (0, 0)

    output = _run("simulations/neural/veffect_controls.py")
    assert "N=20 c=0.05 eps=1e-06 direct: K_act=7 K_corr=62" in output
    assert "N=20 c=0.05 eps=1e-06 transpose: K_act=7 K_corr=64" in output
    drive_rows = module.drive_solver_controls()
    assert tuple(drive_rows) == (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0)
    assert all(row["baseline_counts"] == row["refined_counts"] for row in drive_rows.values())
    assert max(row["refined_residual"] for row in drive_rows.values()) < 1e-12
    assert drive_rows[4.0]["refined_counts"] == (19, 124)
    assert "Transposition preserves the exact spectrum" in output
    assert "solver tol 1e-14: all 13 drive-grid counts stable" in output


def test_withdrawn_lens_direct_run_emits_no_decomposition():
    output = _run("simulations/neural_framework_lens.py")
    assert "withdrawn" in output.lower()
    assert "neural_translation_gate.py" in output
    for forbidden in ("if and only if", "mag-only", "gap closed", "decomposition"):
        assert forbidden.lower() not in output.lower()


def test_repaired_cross_domain_surfaces_keep_scope_fences():
    forbidden_by_path = {
        "docs/THE_BRIDGE_WAS_ALWAYS_OPEN.md": (
            "The interaction is not chaotic",
            "The source is not chaotic",
            "Everything above this line follows from the incompleteness proof",
        ),
        "docs/GAMMA_TIME_DISTINCTION.md": (
            "gamma == experienced time",
            "gamma is the source of experienced time",
            "necessary and sufficient condition for experienced time",
        ),
        "experiments/PI_PAIR_FLUX_BALANCE.md": (
            "N ≡ 4 (mod 10)",
            "18 self-Π",
        ),
    }
    for relative_path, forbidden in forbidden_by_path.items():
        contents = (ROOT / relative_path).read_text(encoding="utf-8").lower()
        for phrase in forbidden:
            assert phrase.lower() not in contents


def test_gamma_time_producer_is_withdrawal_safe():
    output = _run("simulations/gamma_is_time_proof.py")
    assert "withdrawn" in output.lower()
    assert "docs/GAMMA_TIME_DISTINCTION.md" in output
    for forbidden in (
        "experienced time",
        "completeness",
        "exclusivity",
        "qed",
        "time exists",
        "time does not exist",
    ):
        assert forbidden.lower() not in output.lower()

    recorded = (ROOT / "simulations/results/gamma_is_time_proof.txt").read_text(
        encoding="utf-8"
    )
    assert recorded == output
