from pathlib import Path
import sys

import pytest


SIMULATIONS = Path(__file__).parents[1]
if str(SIMULATIONS) not in sys.path:
    sys.path.insert(0, str(SIMULATIONS))

import docs_verify


CASES = (
    (
        "pt_symmetry_analysis.py",
        "results/pt_symmetry_analysis.txt",
        '"results", "pt_symmetry_analysis.txt"',
    ),
    (
        "pt_palindrome_breaking.py",
        "results/pt_palindrome_breaking.txt",
        '"results", "pt_palindrome_breaking.txt"',
    ),
    (
        "entropy_production.py",
        "results/entropy_production.txt",
        '"results", "entropy_production.txt"',
    ),
    (
        "rmt_analysis.py",
        "results/rmt_analysis.txt",
        'RESULTS / "rmt_analysis.txt"',
    ),
    (
        "spectral_form_factor.py",
        "results/spectral_form_factor.txt",
        '"results", "spectral_form_factor.txt"',
    ),
    (
        "fragile_bridge_ep_signature.py",
        "results/fragile_bridge_ep_signature.txt",
        '"results", "fragile_bridge_ep_signature.txt"',
    ),
)


@pytest.mark.parametrize("producer_name,artifact_name,path_expression", CASES)
def test_producer_output_path_resolves_to_canonical_artifact(
    producer_name, artifact_name, path_expression
):
    producer_path = SIMULATIONS / producer_name
    artifact_path = SIMULATIONS / artifact_name
    source = producer_path.read_text(encoding="utf-8")

    assert docs_verify.producer_artifact_path_errors(
        source, producer_path, artifact_path
    ) == []


@pytest.mark.parametrize("producer_name,artifact_name,path_expression", CASES)
def test_output_path_mutation_is_rejected(
    producer_name, artifact_name, path_expression
):
    producer_path = SIMULATIONS / producer_name
    artifact_path = SIMULATIONS / artifact_name
    source = producer_path.read_text(encoding="utf-8")
    changed = source.replace(path_expression, path_expression.replace(".txt", "_moved.txt"))

    assert changed != source
    assert docs_verify.producer_artifact_path_errors(
        changed, producer_path, artifact_path
    )


def test_noncanonical_writer_cannot_be_masked_by_an_unused_canonical_sink():
    producer_path = SIMULATIONS / "entropy_production.py"
    artifact_path = SIMULATIONS / "results" / "entropy_production.txt"
    source = producer_path.read_text(encoding="utf-8")
    changed = source.replace(
        '_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)',
        'MOVED_PATH = OUT_PATH.replace(".txt", "_moved.txt")\n'
        '_outf = open(MOVED_PATH, "w", encoding="utf-8", buffering=1)\n'
        '_decoy = open(OUT_PATH, "a", encoding="utf-8")',
    )

    assert changed != source
    errors = docs_verify.producer_artifact_path_errors(
        changed, producer_path, artifact_path
    )
    assert errors
    assert any("MOVED_PATH" in error or "unused writable" in error for error in errors)


def test_phase3_integrity_rejects_source_and_numeric_artifact_mutations():
    producer_path = SIMULATIONS / "pt_palindrome_breaking.py"
    artifact_path = SIMULATIONS / "results" / "pt_palindrome_breaking.txt"
    source = producer_path.read_text(encoding="utf-8")
    matcher = (SIMULATIONS / "pt_multiset_matching.py").read_text(encoding="utf-8")
    artifact = artifact_path.read_text(encoding="utf-8")

    assert docs_verify.phase3_artifact_integrity_errors(source, artifact, matcher) == []
    assert docs_verify.phase3_artifact_integrity_errors(
        source + "\n# producer mutation", artifact, matcher
    )
    assert docs_verify.phase3_artifact_integrity_errors(
        source, artifact, matcher + "\n# matcher mutation"
    )
    assert docs_verify.phase3_artifact_integrity_errors(
        source, artifact.replace("0.1873101", "0.9999999", 1), matcher
    )
