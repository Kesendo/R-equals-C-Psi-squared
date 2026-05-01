"""Handshake-algebra workflow: verify K-partnership numerically.

The Handshake Algebra (hypotheses/HANDSHAKE_ALGEBRA.md) names a structural
fact: two callers who independently specify the same handshake tuple
(N, k, t, basis) on a shared palindromic chain ARE the handshake — no
exchange step. The algebra is the agreement; there is no send primitive.

This workflow verifies the K-partnership branch of that algebra: bonding
modes k and N+1-k yield identical mirror-pair |·|²-observables under
bipartite NN-hopping with real H. Numerical check parallels
`simulations/_pi_partner_identity.py` at cockpit-N (≤5).

Public API:
  verify_k_partnership(chain, k, t_grid=None, atol=1e-12)
"""
from __future__ import annotations

import numpy as np

from ..pauli import bonding_mode_state
from ..symmetry import k_partner
from .ptf import _propagation_setup, _purity_trajectory, _site_paulis


def verify_k_partnership(chain, k, t_grid=None, atol=1e-12):
    """Verify K-partnership numerically: per-site purity for bonding:k and
    bonding:(N+1-k) coincide under chain.L.

    Per-site purity P_i = ½(1 + ⟨X_i⟩² + ⟨Y_i⟩² + ⟨Z_i⟩²) is K-invariant
    (the |·|² squares away the (-1)^ℓ Bloch-component sign that K induces).
    Identity should hold to machine precision for chain.H_type='xy' on
    uniform-J chains, and for any non-uniform real-J XX/XY profile.

    Args:
        chain: ChainSystem.
        k: bonding-mode index, 1 ≤ k ≤ chain.N.
        t_grid: time samples. Default linspace(0, 5, 50).
        atol: tolerance for 'partnership_holds' verdict.

    Returns:
        dict with:
          'k', 'k_partner': the two bonding-mode indices.
          'self_partner': True if k == k_partner (odd-N midpoint).
          'max_delta_purity': max |P_k(i, t) − P_partner(i, t)| over (i, t).
          'partnership_holds': bool (max_delta < atol).
          't_grid', 'P_k', 'P_partner': trajectories.
    """
    N = chain.N
    k_p = k_partner(N, k)
    self_partner = (k == k_p)

    if t_grid is None:
        t_grid = np.linspace(0.0, 5.0, 50)

    psi_k = bonding_mode_state(N, k)
    rho_k = np.outer(psi_k, psi_k.conj())
    site_paulis = _site_paulis(N)
    decomp_k = _propagation_setup(chain.L, rho_k)
    P_k = _purity_trajectory(decomp_k[0], decomp_k[1], decomp_k[3],
                              t_grid, site_paulis)

    if self_partner:
        P_partner = P_k.copy()
        max_delta = 0.0
    else:
        psi_p = bonding_mode_state(N, k_p)
        rho_p = np.outer(psi_p, psi_p.conj())
        decomp_p = _propagation_setup(chain.L, rho_p)
        P_partner = _purity_trajectory(decomp_p[0], decomp_p[1], decomp_p[3],
                                        t_grid, site_paulis)
        max_delta = float(np.abs(P_k - P_partner).max())

    return {
        'k': k,
        'k_partner': k_p,
        'self_partner': self_partner,
        'max_delta_purity': max_delta,
        'partnership_holds': max_delta < atol,
        't_grid': t_grid,
        'P_k': P_k,
        'P_partner': P_partner,
    }
