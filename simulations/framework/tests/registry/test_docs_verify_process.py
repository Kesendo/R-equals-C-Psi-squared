"""Process-contract regression for the scoped documentation smoke test."""

from pathlib import Path
import subprocess
import sys


def test_forced_failed_check_returns_nonzero():
    """Fails if docs_verify prints a FAIL but regresses to process exit zero."""
    script = Path(__file__).resolve().parents[3] / "docs_verify.py"
    result = subprocess.run(
        [sys.executable, str(script), "--force-failure"],
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )

    assert "FAIL: forced failure proves the process exits nonzero" in result.stdout
    assert result.returncode != 0
