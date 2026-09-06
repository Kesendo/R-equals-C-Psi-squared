"""Deterministic report regression for the N-to-infinity producer."""

from pathlib import Path
import subprocess
import sys


SIMULATIONS = Path(__file__).resolve().parents[3]
PRODUCER = SIMULATIONS / "n_infinity_analysis.py"
REPORT = SIMULATIONS / "results" / "n_infinity_analysis.txt"


def run_producer():
    result = subprocess.run(
        [sys.executable, str(PRODUCER)],
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    return REPORT.read_bytes()


def test_report_matches_snapshot_and_is_byte_stable_across_two_runs():
    expected = REPORT.read_bytes()
    first = run_producer()
    second = run_producer()

    assert first == expected
    assert second == first
    assert b"Date:" not in second
    assert b" eigenvalues," not in second
