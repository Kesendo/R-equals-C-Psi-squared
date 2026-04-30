"""Tests for ChainSystem.gamma_probe_setup and estimate_gamma_from_cpsi."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_gamma_probe_setup_at_default_gamma():
    """Verify the optimal probe parameters at γ=0.05 match the manual analysis."""
    chain = fw.ChainSystem(N=2)
    setup = fw.gamma_probe_setup(chain,gamma_assumed=0.05, target_precision=0.01)
    # K_optimal ≈ 0.119 (4·γ·t ≈ 0.474 → K = γt ≈ 0.119)
    assert abs(setup['K_optimal'] - 0.119) < 0.01
    # cpsi_target ≈ 0.144
    assert abs(setup['cpsi_target'] - 0.144) < 0.005
    # K_cusp ≈ 0.0374 (4γt_cusp = -ln(f_cusp) ≈ 0.150 → K ≈ 0.0374)
    assert abs(setup['K_cusp'] - 0.0374) < 0.001
    # cpsi_cusp = 0.25 by construction
    assert setup['cpsi_cusp'] == 0.25


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_gamma_probe_K_invariance():
    """K_optimal = γ·t* should be γ-independent (γ-invariance in dimensionless form)."""
    chain = fw.ChainSystem(N=2)
    K_vals = []
    for gamma in [0.01, 0.05, 0.1, 0.5]:
        setup = fw.gamma_probe_setup(chain,gamma_assumed=gamma)
        K_vals.append(setup['K_optimal'])
    # All K_optimal values should agree (γ-invariance)
    assert all(abs(k - K_vals[0]) < 1e-6 for k in K_vals), \
        f"K_optimal varies with γ: {K_vals}"


def test_gamma_probe_default_gamma_uses_chain_gamma_0():
    chain = fw.ChainSystem(N=3, gamma_0=0.07)
    setup = fw.gamma_probe_setup(chain)
    assert setup['gamma_assumed'] == 0.07


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_estimate_gamma_from_cpsi_inverts_f25():
    """Round-trip: γ → CΨ(t) → estimate γ should recover the original."""
    chain = fw.ChainSystem(N=2)
    gamma_true = 0.07
    t = 2.0  # arbitrary probe time
    f = np.exp(-4 * gamma_true * t)
    cpsi = f * (1 + f**2) / 6.0
    gamma_est = fw.estimate_gamma_from_cpsi(chain, cpsi, t)
    assert abs(gamma_est - gamma_true) < 1e-9


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_estimate_gamma_rejects_out_of_range():
    chain = fw.ChainSystem(N=2)
    with pytest.raises(ValueError, match="≥ 1/3"):
        fw.estimate_gamma_from_cpsi(chain, 0.5, t=1.0)
    with pytest.raises(ValueError, match="≤ 0"):
        fw.estimate_gamma_from_cpsi(chain, -0.1, t=1.0)
    with pytest.raises(ValueError, match="t must be > 0"):
        fw.estimate_gamma_from_cpsi(chain, 0.1, t=0)


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_gamma_probe_setup_x_channel():
    """gamma_probe_setup with channel='X' should give K_cusp = 0.0867 vs Z's 0.0374."""
    chain = fw.ChainSystem(N=2)
    setup_z = fw.gamma_probe_setup(chain, gamma_assumed=0.05, channel='Z')
    setup_x = fw.gamma_probe_setup(chain, gamma_assumed=0.05, channel='X')
    setup_y = fw.gamma_probe_setup(chain, gamma_assumed=0.05, channel='Y')
    setup_d = fw.gamma_probe_setup(chain, gamma_assumed=0.05, channel='depolarizing')
    # K_cusp from F26 cusp condition (note: K_Y = K_Z, NOT K_X — doc has typo)
    assert abs(setup_z['K_cusp'] - 0.0374) < 0.001
    assert abs(setup_x['K_cusp'] - 0.0867) < 0.001
    assert abs(setup_y['K_cusp'] - 0.0374) < 0.001  # K_Y = K_Z, not K_X
    assert abs(setup_d['K_cusp'] - 0.0440) < 0.001


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_estimate_gamma_round_trip_x_channel():
    """Round-trip γ → CΨ_X(t) → estimate_γ for X-channel."""
    chain = fw.ChainSystem(N=2)
    gamma_true = 0.07
    t = 2.0
    cpsi = fw.cpsi_bell_plus(gamma_true, 0.0, 0.0, t)
    gamma_est = fw.estimate_gamma_from_cpsi(chain, cpsi, t, channel='X')
    assert abs(gamma_est - gamma_true) < 1e-9


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_gamma_probe_setup_kingston_data_consistency():
    """Kingston cusp-slowing F25 RMS residual was 0.0097; with 1% target,
    shots-needed should be reasonable (~10⁵-10⁶)."""
    chain = fw.ChainSystem(N=2)
    setup = fw.gamma_probe_setup(chain,gamma_assumed=0.05, target_precision=0.01)
    # Order of magnitude check: should be ~10⁶ shots for 1% precision
    assert 1e5 < setup['shots_needed'] < 1e7
