"""Tests for cockpit_panel (Lebensader analysis: skeleton + trace + cusp + chiral + Y-parity)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_cockpit_panel_yzzy_t1_drop_28_at_n3():
    """Reproduces the EQ-030 hardware-confirmed drop=28 for YZ+ZY +T1=0.005 at N=3."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(plus, minus), plus)
    r = fw.Receiver(psi, chain=chain)
    result = chain.cockpit_panel(r, terms=[('Y', 'Z'), ('Z', 'Y')],
                                  gamma_t1=0.005, t_max=5.0, dt=0.01)
    skeleton = result['lebensader']['skeleton']
    assert skeleton['drop'] == 28
    assert result['lebensader']['rating'].startswith('collapsed')


def test_cockpit_panel_truly_no_drop_pure_z():
    """XX+YY (truly) under pure-Z: no skeleton drop expected."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(plus, minus), plus)
    r = fw.Receiver(psi, chain=chain)
    result = chain.cockpit_panel(r, terms=[('X', 'X'), ('Y', 'Y')],
                                  t_max=5.0, dt=0.01)
    skeleton = result['lebensader']['skeleton']
    assert skeleton['drop'] == 0


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_cockpit_panel_rejects_non_hermitian_rho():
    chain = fw.ChainSystem(N=2, gamma_0=0.05)
    rho_nh = np.array([[0.5, 0.2, 0, 0],
                       [0.2, 0.3, 0, 0],
                       [0, 0, 0.1, 0.05],
                       [0, 0, 0.05j, 0.1]], dtype=complex)
    with pytest.raises(ValueError, match="Hermitian"):
        fw.cockpit_panel(chain.H, [0.05]*2, rho_nh, 2, t_max=0.5, dt=0.05)


def test_cockpit_panel_rejects_non_psd_rho():
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    rho_neg = np.zeros((8,8), dtype=complex)
    rho_neg[0,0] = 1.5
    rho_neg[1,1] = -0.5  # negative eigenvalue
    with pytest.raises(ValueError, match="positive semi-definite"):
        fw.cockpit_panel(chain.H, [0.05]*3, rho_neg, 3, t_max=0.5, dt=0.05)


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_cockpit_panel_rejects_trace_mismatch():
    chain = fw.ChainSystem(N=2, gamma_0=0.05)
    rho_no_trace = np.eye(4, dtype=complex) * 0.5  # trace = 2
    with pytest.raises(ValueError, match="trace"):
        fw.cockpit_panel(chain.H, [0.05]*2, rho_no_trace, 2, t_max=0.5, dt=0.05)


def test_cockpit_panel_terms_uses_chain_J():
    """terms-mode must scale by self.J, not hardcoded 1.0.

    A chain at J=2 with terms=YZ+ZY should give a different residual
    (scaled by J²=4) than the same chain at J=1.
    """
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(plus, minus), plus)

    chain1 = fw.ChainSystem(N=3, J=1.0, gamma_0=0.05)
    r1 = fw.Receiver(psi, chain=chain1)
    p1 = chain1.cockpit_panel(r1, terms=[('Y','Z'),('Z','Y')], gamma_t1=0.005,
                              t_max=2.0, dt=0.01)

    chain2 = fw.ChainSystem(N=3, J=2.0, gamma_0=0.05)
    r2 = fw.Receiver(psi, chain=chain2)
    p2 = chain2.cockpit_panel(r2, terms=[('Y','Z'),('Z','Y')], gamma_t1=0.005,
                              t_max=2.0, dt=0.01)

    # Different J → different θ-trajectory (Hamiltonian rescales coherent dynamics)
    theta1 = p1['_trajectory_for_inspection']['theta']
    theta2 = p2['_trajectory_for_inspection']['theta']
    assert not np.allclose(theta1, theta2), \
        "cockpit_panel(terms=...) must use self.J, not hardcoded 1.0"
