"""Tests for the bridge_panel workflow.

Verifies the six-angle aggregation structure: each sub-block returns
the expected reading on canonical chain configurations.

Reference: docs/THE_BRIDGE_WAS_ALWAYS_OPEN.md (canonical March-22 doc)
plus the 2026-05-01 six structural angles.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_bridge_panel_xy_chain_fully_visible():
    """Default XY chain (uniform J, γ₀-const, real H, bipartite NN-hopping)
    confirms all six angles. Bridge is fully visible.
    """
    chain = fw.ChainSystem(N=4, H_type='xy')
    panel = fw.bridge_panel(chain)
    assert panel['palindrome']['palindromic']
    assert panel['palindrome']['M_norm'] < 1e-10
    assert panel['taktgeber']['bidirectional']
    assert panel['k_partnership']['expected_to_hold']
    assert panel['bridge_visible']
    assert panel['occluded_at'] == []


def test_bridge_panel_heisenberg_chain_k_partnership_breaks():
    """Heisenberg open chain: K-partnership structurally breaks at boundary
    (per HANDSHAKE_ALGEBRA robustness ladder). Bridge is occluded at the
    K-partnership angle.
    """
    chain = fw.ChainSystem(N=4, H_type='heisenberg')
    panel = fw.bridge_panel(chain)
    # F1 palindrome still holds on chain.L (Heisenberg + Z-dephasing has
    # palindromic spectrum)
    assert panel['palindrome']['palindromic']
    # But K-partnership expectation is False
    assert not panel['k_partnership']['expected_to_hold']
    assert not panel['bridge_visible']
    assert 'K_partnership_structural_expected' in panel['occluded_at']


def test_bridge_panel_pi_decomposition_zero_for_default_h():
    """For default chain Hamiltonians (XX+YY or XX+YY+ZZ), all bilinears
    are Klein (0,0) Π²-even, so M_anti = 0 and anti_fraction = 0.
    """
    chain = fw.ChainSystem(N=4, H_type='xy')
    panel = fw.bridge_panel(chain)
    assert panel['pi_decomposition']['anti_fraction'] < 1e-10
    assert panel['pi_decomposition']['M_anti_norm'] < 1e-10


def test_bridge_panel_with_rho_0_adds_polarity():
    """Passing rho_0 adds the polarity sub-block."""
    chain = fw.ChainSystem(N=3, H_type='xy')
    psi = fw.polarity_state(3, +1)
    panel = fw.bridge_panel(chain, rho_0=psi)
    assert 'polarity' in panel
    assert abs(panel['polarity']['aggregate_polarity'] - 1.0) < 1e-12
    assert not panel['polarity']['on_boundary'].any()


def test_bridge_panel_bonding_mode_on_boundary():
    """Bonding-mode initial state sits on the boundary (X=0 per site)."""
    chain = fw.ChainSystem(N=4, H_type='xy')
    psi = fw.bonding_mode_state(4, 1)
    panel = fw.bridge_panel(chain, rho_0=psi)
    assert panel['polarity']['on_boundary'].all()
    assert abs(panel['polarity']['aggregate_polarity']) < 1e-12


def test_bridge_panel_with_terms_adds_klein_inheritance():
    """Passing terms adds the Klein-inheritance sub-block."""
    chain = fw.ChainSystem(N=4, H_type='xy')
    panel = fw.bridge_panel(chain, terms=[('I', 'I', 'Z'), ('Z', 'Z', 'Z')])
    assert 'klein_inheritance' in panel
    klein = panel['klein_inheritance']
    assert klein['klein_set'] == {(0, 1)}
    assert klein['dissipator_resonance_cell'] == (0, 1)
    assert klein['h_in_resonance']
    assert klein['klein_homogeneous']


def test_bridge_panel_klein_inhomogeneous_terms():
    """Klein-inhomogeneous terms (mixed Klein cells) show klein_homogeneous=False.

    XY: X(1,0)+Y(1,1) = (0,1). YZ: Y(1,1)+Z(0,1) = (1,0). Two different
    Klein cells → not homogeneous.
    """
    chain = fw.ChainSystem(N=4, H_type='xy')
    panel = fw.bridge_panel(chain, terms=[('X', 'Y'), ('Y', 'Z')])
    klein = panel['klein_inheritance']
    assert klein['klein_set'] == {(0, 1), (1, 0)}
    assert not klein['klein_homogeneous']


def test_bridge_panel_verify_k_runs_partnership_check():
    """verify_k_at parameter triggers numerical K-partnership verification."""
    chain = fw.ChainSystem(N=4, H_type='xy')
    panel = fw.bridge_panel(chain, verify_k_at=1)
    assert 'verification' in panel['k_partnership']
    assert panel['k_partnership']['verification']['partnership_holds']


def test_bridge_panel_angles_list_has_six_entries():
    """The angles list always has exactly six entries with correct kind labels."""
    chain = fw.ChainSystem(N=3, H_type='xy')
    panel = fw.bridge_panel(chain)
    assert len(panel['angles']) == 6
    names = {a['name'] for a in panel['angles']}
    assert names == {
        'F1_palindrome',
        'Handshake_idempotence',
        'Channel_not_memory_pi_decomp',
        'One_system_two_indices',
        'Algebra_is_inheritance',
        'Bidirectional_Taktgeber',
    }
