"""Deterministic report regression for the N-to-infinity producer.

The producer writes to a temporary path (--out), so running the test never
rewrites the committed results file it compares against.
"""

from pathlib import Path
import subprocess
import sys


SIMULATIONS = Path(__file__).resolve().parents[3]
PRODUCER = SIMULATIONS / "n_infinity_analysis.py"
REPORT = SIMULATIONS / "results" / "n_infinity_analysis.txt"


def run_producer(out_path):
    result = subprocess.run(
        [sys.executable, str(PRODUCER), "--out", str(out_path)],
        capture_output=True,
        text=True,
        timeout=600,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    return out_path.read_bytes()


def test_report_matches_snapshot_and_is_byte_stable_across_two_runs(tmp_path):
    expected = REPORT.read_bytes()
    first = run_producer(tmp_path / "first.txt")
    second = run_producer(tmp_path / "second.txt")

    assert REPORT.read_bytes() == expected
    assert first == expected
    assert second == first
    assert b"Date:" not in second
    assert b" eigenvalues," not in second
    assert b"-0.0" not in second
