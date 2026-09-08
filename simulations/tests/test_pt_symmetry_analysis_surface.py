from pathlib import Path
import sys

import numpy as np


SCRIPT = Path(__file__).parents[1] / "pt_symmetry_analysis.py"
RESULT = Path(__file__).parents[1] / "results" / "pt_symmetry_analysis.txt"
PHASE3_SCRIPT = Path(__file__).parents[1] / "pt_palindrome_breaking.py"
PHASE3_RESULT = Path(__file__).parents[1] / "results" / "pt_palindrome_breaking.txt"
DOCUMENT = Path(__file__).parents[2] / "experiments" / "PT_SYMMETRY_ANALYSIS.md"
EXPERIMENTS_README = Path(__file__).parents[2] / "experiments" / "README.md"
SIMULATIONS = Path(__file__).parents[1]
if str(SIMULATIONS) not in sys.path:
    sys.path.insert(0, str(SIMULATIONS))

from pt_multiset_matching import multiset_reflection_error


def test_output_surface_reports_sectorwise_p_and_open_full_class():
    source = SCRIPT.read_text(encoding="utf-8")

    assert "sectorwise P" in source
    assert "full irreducible SRP class remains OPEN" in source
    assert "PALINDROME AXIS DEPARTURE" in source
    assert "single-eigenvector Petermann factors or angles" in source
    assert "operator relation remains exact" in source
    assert "within palindrome axis" in source
    assert "off palindrome axis; operator relation remains exact" in source
    assert "chiral phase" not in source
    assert "chiral breaking" not in source
    assert "CHIRAL SYMMETRY BREAKING" not in source
    assert "Class AIII" not in source
    assert "class AIII" not in source
    assert "Eigenvector coalescence" not in source
    assert "MIRROR-PARTNER ANGLE CONTROL" not in source
    assert "Petermann factor K and phase rigidity" not in source
    document = DOCUMENT.read_text(encoding="utf-8")
    for stale in ("CHIRAL (sublattice) symmetry", "chiral symmetry of the full operator structure",
                  "confirming the chiral symmetry"):
        assert stale not in document
    normalized = " ".join(document.split())
    assert "negative P-type relation" in normalized
    assert "does not assign a global irreducible symmetry class" in normalized


def test_generated_artifact_has_current_sectorwise_scope():
    output = RESULT.read_text(encoding="utf-8")

    assert "sectorwise P" in output
    assert "full irreducible SRP class remains OPEN" in output
    assert "palindrome axis" in output
    assert "operator relation remains exact" in output
    assert "AIII" not in output
    assert "chiral breaking" not in output.lower()
    assert "Hopf" not in output
    assert "No exceptional point" not in output
    assert "No classical EP" not in output
    assert "single-eigenvector Petermann factors or angles" in output
    assert "MIRROR-PARTNER ANGLE CONTROL" not in output
    assert "Petermann factor K and phase rigidity" not in output


def test_document_retires_basis_dependent_degenerate_eigenvector_readings():
    document = DOCUMENT.read_text(encoding="utf-8")
    normalized = " ".join(document.split())
    assert "does not report single-eigenvector Petermann factors or eigenvector angles" in normalized
    assert "basis-dependent inside a degenerate eigenspace" in normalized
    assert "invariant subspace-level conditioning" in normalized
    assert "| 1.463 | **402.7**" not in document
    assert "Above γ_crit: cos(theta) ~ 0.09" not in document


def test_document_does_not_infer_finite_time_norm_behavior_from_spectral_axis():
    document = DOCUMENT.read_text(encoding="utf-8")
    normalized = " ".join(document.split())
    assert "spectral-abscissa statement" in normalized
    assert "non-normal" in normalized
    assert "transient norm growth can occur" in normalized
    assert "not a physical Lindblad channel" in normalized
    assert "the system oscillates but does not grow" not in document
    assert "This is the explosion" not in document


def test_phase3_uses_trace_midpoint_and_names_only_axis_departure():
    source = PHASE3_SCRIPT.read_text(encoding="utf-8")
    output = PHASE3_RESULT.read_text(encoding="utf-8")
    document = DOCUMENT.read_text(encoding="utf-8")
    combined = source + output + document

    assert "midpoint = np.trace(L) / L.shape[0]" in source
    assert "from pt_multiset_matching import multiset_reflection_error" in source
    assert "trace(L)/dim(L) midpoint" in combined
    assert "imaginary-axis regime" in combined
    assert "chiral phase" not in combined.lower()
    assert "r = +0.987 correlation" in document
    assert "r = +0.988 correlation" not in document
    assert "| 0.300 | 7.9e-01 |" in document
    assert "at 1.01 γ_crit" in document
    assert "Finite-offset most-unstable-mode check" in source
    assert "No branch continuation is performed" in source
    assert "instability-branch check" not in source
    assert "unstable branch oscillates" not in source


def test_phase3_pairing_error_is_a_bijective_multiset_distance():
    asymmetric = np.array([1, 1, 1, -1, 2, -2, -2], dtype=complex)
    assert multiset_reflection_error(asymmetric, 0.0) == 2.0


def test_phase2_and_phase3_share_the_bijective_multiset_matcher():
    for script in (SCRIPT, PHASE3_SCRIPT):
        source = script.read_text(encoding="utf-8")
        assert "from pt_multiset_matching import multiset_reflection_error" in source
    phase2 = SCRIPT.read_text(encoding="utf-8")
    assert "max_pe = multiset_reflection_error(ev, midpoint=0.0)" in phase2
    assert "np.min(np.abs(ev - (-lam)))" not in phase2


def test_phase3_document_does_not_assign_causality_to_composite_epsilon_sweep():
    document = DOCUMENT.read_text(encoding="utf-8")

    for required in (
        "does not isolate the causal contribution",
        "2γ - 4epsilon/3",
        "eigenmode property of the full coupled generator",
        "scalar bottleneck distance",
        "does not identify a causal stability role",
        "Pi D Pi^-1 = -D - 2kappa I",
        "`D + kappa I` anti-commutes with `Pi`",
    ):
        assert required in document
    for stale in (
        "not the protection mechanism",
        "Breaking the palindrome makes the system MORE stable",
        "-(gamma - epsilon/3)",
        "Effect (2) dominates",
        "It makes it more stable",
        "damping wins",
        "Palindrome error = 2 * L_X",
        "protect the phase",
        "pairing pins every eigenvalue",
        "a pair meets itself in the mirror",
        "The Z and Y components anti-commute with Pi",
    ):
        assert stale not in document


def test_phase3_summary_tracks_current_correlation_without_causal_promotion():
    summary = EXPERIMENTS_README.read_text(encoding="utf-8")
    pt_row = next(line for line in summary.splitlines() if "[PT-Symmetry Analysis]" in line)

    assert "r = +0.987" in pt_row
    assert "causal contributions are not isolated" in pt_row
    assert "+0.988" not in pt_row
    assert "stabilizes the sampled system" not in pt_row


def test_phase3_artifact_is_machine_independent():
    source = PHASE3_SCRIPT.read_text(encoding="utf-8")
    output = PHASE3_RESULT.read_text(encoding="utf-8")

    assert "clock.time" not in source
    assert "Results: {OUT_PATH}" not in source
    assert str(Path(__file__).parents[2]) not in output
    assert "Source SHA256:" in output
    assert "Payload SHA256:" in output
