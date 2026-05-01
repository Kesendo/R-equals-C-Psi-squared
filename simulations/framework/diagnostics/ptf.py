"""PTF (Perspectival Time Field) structural readings.

For a slow-mode basis {M_s, W_s, λ_s} of the unperturbed Liouvillian L_A and
a bond-perturbation V_L = ∂L/∂J_b ("the dynamics of the dynamics"), this module
gives the first-order perturbation-theory matrix elements that drive the
per-site rate-of-painting α_i.

Background: hypotheses/PERSPECTIVAL_TIME_FIELD.md (Tier 2, downgraded from
Tier 1 by EQ-014). The α_i pattern is state-dependent; the closure law
Σ_i ln(α_i) ≈ 0 holds within ±0.05 in the perturbative window |δJ| ≤ 0.1
across five qualitatively distinct initial states tested at N=7. See
ON_THE_PAINTER_PRINCIPLE.md for the painter framing.

Public API:
  pt_matrix_elements(slow_modes_dict, V_L)
  pt_eigvec_shift(slow_modes_dict, V_L, mode_idx)
"""
from __future__ import annotations

import numpy as np


def pt_matrix_elements(slow_modes_dict, V_L):
    """First-order matrix elements ⟨W_s | V_L | M_{s'}⟩ on the slow-mode basis.

    Args:
        slow_modes_dict: output of `framework.slow_modes(chain, ...)` —
            must have 'right_eigvecs', 'left_covecs', 'eigenvalues'.
        V_L: bond-perturbation superoperator, e.g. from `bond_perturbation`.

    Returns:
        n_slow × n_slow complex matrix where entry [s, s'] = ⟨W_s | V_L | M_{s'}⟩.
        Diagonal entries are the first-order eigenvalue shifts δλ_s.
        Off-diagonal entries drive eigenvector mixing.
    """
    M = slow_modes_dict['right_eigvecs']     # d² × n_slow
    W = slow_modes_dict['left_covecs']        # n_slow × d²
    return W @ V_L @ M


def pt_eigvec_shift(slow_modes_dict, V_L, mode_idx, degenerate_tol=1e-12):
    """First-order eigenvector shift δM_s under perturbation V_L.

    δM_s = Σ_{s' ≠ s} [⟨W_{s'} | V_L | M_s⟩ / (λ_s − λ_{s'})] · M_{s'}

    Degenerate pairs (|λ_s − λ_{s'}| < tol) are excluded — first-order
    perturbation theory diverges there and a degenerate-subspace projection
    would be needed for a proper treatment.

    Args:
        slow_modes_dict: output of `framework.slow_modes(...)`.
        V_L: bond-perturbation superoperator.
        mode_idx: which slow mode (index into slow_modes_dict).
        degenerate_tol: skip terms with |λ_s − λ_{s'}| below this threshold.

    Returns:
        d² complex vector — the first-order shift to mode M_s under V_L.
    """
    M = slow_modes_dict['right_eigvecs']
    W = slow_modes_dict['left_covecs']
    evals = slow_modes_dict['eigenvalues']
    n_slow = len(evals)
    M_s = M[:, mode_idx]
    matrix_col = W @ V_L @ M_s            # ⟨W_{s'} | V_L | M_s⟩ for all s'
    denoms = evals[mode_idx] - evals
    mask = (np.arange(n_slow) != mode_idx) & (np.abs(denoms) > degenerate_tol)
    coeffs = np.zeros(n_slow, dtype=complex)
    coeffs[mask] = matrix_col[mask] / denoms[mask]
    return M @ coeffs
