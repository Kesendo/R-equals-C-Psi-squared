"""Regression gates for executable claim surfaces repaired by the neural review."""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path


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
    output = _run("simulations/factor_two_standing_waves.py")
    assert "linear F1 total:       pairs=10903 fixed=34" in output
    assert "conjugate-composite:   pairs=9921 fixed=1998" in output
    assert "N=4 chain     linear: pairs=121 fixed=14" in output
    assert "N=4 star      linear: pairs=120 fixed=16" in output
    assert "N=4 ring      linear: pairs=116 fixed=24" in output
    for forbidden in ("finesse", "Beer-Lambert", "unpaired", "halve absorption"):
        assert forbidden.lower() not in output.lower()


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
