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

    # Exact route: 2*0.05 == 0.10 in binary, so the joint scaling doubles every generator entry
    # exactly; the generator is compared, not its eigensolver output.
    base = fw.ChainSystem(N=3, J=1.0, gamma_0=0.05).L
    scaled = fw.ChainSystem(N=3, J=2.0, gamma_0=0.10).L
    assert np.array_equal(scaled, 2 * base)
    assert "np.array_equal(chain.L, 2 * base_chain.L)" in source


def test_demo_does_not_hide_fixed_j_drift_behind_nearest_grid_tolerance():
    source = SCRIPT.read_text(encoding="utf-8")
    assert "J=1.0, gamma_0=gamma_10" not in source
    assert "abs(closest - pred) < 1e-4" not in source
