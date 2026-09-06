"""Standalone safety and immutability gates for retired XOR event records."""

from pathlib import Path
import subprocess
import sys

import pytest


SIMULATIONS = Path(__file__).resolve().parents[3]
WARNING = "ARCHIVED EVENT RECORD - RETIRED INTERPRETATION"
ACTIVE_GATE = "simulations/f22_operator_charge.py"
ACTIVE_OUTPUT = "simulations/results/error_correction_palindrome_f22.txt"


@pytest.mark.parametrize(
    ("artifact", "first_retired_claim"),
    [
        ("error_correction_palindrome.txt", "carry most of the info lifetime"),
        ("xor_non_heisenberg_v2.txt", "GHZ showed 50% instead of 100%"),
    ],
)
def test_archival_warning_precedes_retired_claims(artifact, first_retired_claim):
    text = (SIMULATIONS / "results" / artifact).read_text(encoding="utf-8")

    assert text.startswith(WARNING)
    assert "They are not state weights." in text
    assert "No GHZ/W ranking" in text
    assert ACTIVE_GATE in text
    assert ACTIVE_OUTPUT in text
    assert 0 <= text.index(WARNING) < text.index(first_retired_claim)


@pytest.mark.parametrize(
    ("producer", "artifact"),
    [
        ("error_correction_palindrome.py", "error_correction_palindrome.txt"),
        ("xor_non_heisenberg_v2.py", "xor_non_heisenberg_v2.txt"),
    ],
)
def test_retired_producer_does_not_overwrite_annotated_record(producer, artifact):
    artifact_path = SIMULATIONS / "results" / artifact
    before = artifact_path.read_bytes()

    result = subprocess.run(
        [sys.executable, str(SIMULATIONS / producer)],
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )

    assert result.returncode == 0
    assert "RETIRED" in result.stdout
    assert artifact_path.read_bytes() == before
