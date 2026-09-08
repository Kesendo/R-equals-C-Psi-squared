from pathlib import Path

import numpy as np

from simulations import framework as fw


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "simulations" / "universal_carrier_demo.py"


def test_carrier_scaling_holds_q_fixed_by_scaling_j_with_gamma():
    source = SCRIPT.read_text(encoding="utf-8")
    assert "fixed Q=J/gamma=20" in source
    assert "J=2.0, gamma_0=gamma_10" in source
    assert "J=1.0, gamma_0=gamma_10" not in source

    base = np.sort(np.linalg.eigvals(fw.ChainSystem(N=3, J=1.0, gamma_0=0.05).L).real)
    scaled = np.sort(np.linalg.eigvals(fw.ChainSystem(N=3, J=2.0, gamma_0=0.10).L).real)
    np.testing.assert_allclose(scaled, 2 * base, atol=1e-10, rtol=1e-10)


def test_demo_does_not_hide_fixed_j_drift_behind_nearest_grid_tolerance():
    source = SCRIPT.read_text(encoding="utf-8")
    assert "J=1.0, gamma_0=gamma_10" not in source
    assert "abs(closest - pred) < 1e-4" not in source
