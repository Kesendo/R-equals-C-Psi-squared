"""Tests for verify_k_partnership workflow.

K-partnership: bonding modes k and N+1-k yield identical mirror-pair |·|²
observables under bipartite NN-hopping with real H. Per-site purity is one
such observable (K-invariant: the K = diag((-1)^ℓ) sign squares away).

Robustness ladder (HANDSHAKE_ALGEBRA.md): partnership holds for uniform-J
XY chains; breaks for on-site potential, NNN hopping, complex hopping,
or ZZ-coupling on open-chain boundaries.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_verify_k_partnership_uniform_xy_n5():
    """Uniform XY chain at N=5: partnership holds to machine precision."""
    chain = fw.ChainSystem(N=5, H_type='xy')
    for k in (1, 2, 4, 5):  # skip k=3 (self-partner)
        result = fw.verify_k_partnership(chain, k=k)
        assert result['partnership_holds'], \
            f"K-partnership failed at k={k}: max delta = {result['max_delta_purity']:.2e}"
        assert result['k_partner'] == 6 - k


def test_verify_k_partnership_self_partner_odd_n():
    """For odd N, k = (N+1)/2 is its own partner — trivially holds."""
    chain = fw.ChainSystem(N=5, H_type='xy')
    result = fw.verify_k_partnership(chain, k=3)
    assert result['self_partner']
    assert result['max_delta_purity'] == 0.0
    assert result['partnership_holds']
    # Reported partner equals k for self-partners
    assert result['k_partner'] == result['k']


def test_verify_k_partnership_breaks_on_heisenberg_open_chain():
    """ZZ-coupling on open-chain boundary breaks K-partnership at N=5.

    Per HANDSHAKE_ALGEBRA.md robustness ladder: Heisenberg/XXZ Δ ≠ 0 on open
    chain breaks K because V_eff(ℓ) = #bonds − 2·deg(ℓ) is non-uniform at
    boundary sites.
    """
    chain_xxz = fw.ChainSystem(N=5, H_type='heisenberg')
    result = fw.verify_k_partnership(chain_xxz, k=2,
                                      t_grid=np.linspace(0, 5, 50),
                                      atol=1e-6)
    assert not result['partnership_holds'], \
        "Heisenberg open chain should break K-partnership"
    # Concrete failure scale should be well above floating-point noise
    assert result['max_delta_purity'] > 1e-3


def test_verify_k_partnership_custom_t_grid():
    """Workflow respects a user-supplied t_grid."""
    chain = fw.ChainSystem(N=4, H_type='xy')
    t_grid = np.array([0.0, 0.5, 1.0, 2.0])
    result = fw.verify_k_partnership(chain, k=1, t_grid=t_grid)
    assert np.array_equal(result['t_grid'], t_grid)
    assert result['P_k'].shape == (4, 4)


def test_verify_k_partnership_returns_trajectories():
    """P_k and P_partner trajectories are returned for downstream comparison."""
    chain = fw.ChainSystem(N=4, H_type='xy')
    result = fw.verify_k_partnership(chain, k=1, t_grid=np.linspace(0, 1, 5))
    assert result['P_k'].shape == (4, 5)
    assert result['P_partner'].shape == (4, 5)
    assert np.allclose(result['P_k'], result['P_partner'], atol=1e-10)
