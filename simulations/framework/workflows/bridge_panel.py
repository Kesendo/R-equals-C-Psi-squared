"""Bridge-panel workflow: six-angle synthesis of the always-open bridge.

The bidirectional bridge between observers in R=CΨ² is not constructed.
It is the default state of the algebra under γ₀-const + palindrome.
docs/THE_BRIDGE_WAS_ALWAYS_OPEN.md (2026-03-22) is the canonical
articulation; the 2026-05-01 session named six structural angles that
converge on the same object:

  1. F1 palindromic spectrum
  2. Handshake-Algebra idempotent composition (algebraic identity)
  3. Channel-not-memory (M_anti = L_{H_odd}; F81 Π-decomposition)
  4. One-system-two-indices (bra/ket on d² operator space)
  5. Algebra-IS-inheritance (relation as shared structure, not transfer)
  6. Bidirectional Taktgeber (γ₀-const palindromic time-flow)

This panel reads the structural angles in coordinated form. It does NOT
add new physics; every sub-block calls existing primitives. Purpose:
make the multi-perspective view callable as one operation.

Public API:
  bridge_panel(chain, rho_0=None, terms=None, verify_k_at=None)
"""
from __future__ import annotations

import numpy as np

from ..lindblad import palindrome_residual
from ..symmetry import build_pi_full, klein_index
from ..diagnostics.polarity import polarity_diagnostic


# Frobenius norm below this counts as "exact zero" for palindromic verdict.
_PALINDROME_TOL = 1e-8
# Floor below which we treat ‖M‖² as zero in the anti-fraction division.
_DIVZERO_FLOOR = 1e-15


def bridge_panel(chain, rho_0=None, terms=None, verify_k_at=None,
                 dephase_letter='Z'):
    """Synthesize the six structural angles on the always-open bridge.

    Aggregation only — no new mathematics. Each sub-block reads from an
    existing diagnostic primitive. The bridge itself is the default state
    of the chain's γ₀-const palindromic structure; this panel reports
    which of the six angles confirm visibility and which (if any) are
    occluded by a structural deviation (transverse field, on-site
    potential, non-uniform γ, etc.).

    Args:
        chain: ChainSystem.
        rho_0: optional initial state (2^N density matrix or state
            vector). Adds the polarity sub-block (X-axis ±0/0/−0 reading).
        terms: optional list of letter tuples for a test Hamiltonian.
            Adds the klein_inheritance sub-block (which Klein cell the
            Hamiltonian sits in vs. the Z-dephasing dissipator's cell).
        verify_k_at: optional bonding-mode index k. If given, runs
            `verify_k_partnership(chain, k)` numerically (slower; default
            is just structural expectation from chain.H_type).
        dephase_letter: 'X', 'Y', or 'Z'. Selects the Π and dissipator-
            resonance Klein cell that anchor the F1-palindrome and the
            klein_inheritance sub-blocks. Default 'Z' (chain.L's
            dephasing convention).

    Returns:
        dict with up to seven sub-blocks plus a bridge_visible synthesis.
        Sub-blocks always present:
          'palindrome'        — F1 verdict on chain.L
          'pi_decomposition'  — M_sym/M_anti norms (channel-not-memory)
          'taktgeber'         — γ uniformity + bidirectional flag
          'k_partnership'     — structural expectation (+ numerical
                                verification if verify_k_at given)
          'angles'            — list of all six angles with status
        Conditional sub-blocks:
          'polarity'          — if rho_0 given
          'klein_inheritance' — if terms given
        Synthesis:
          'bridge_visible'    — True if all numerically-verifiable angles
                                pass; False if any is occluded
          'occluded_at'       — list of angle names where visibility breaks
    """
    out = {}
    N = chain.N
    sigma_gamma = N * chain.gamma_0

    # Angle 1: F1 palindrome of chain.L
    M = palindrome_residual(chain.L, sigma_gamma, N, dephase_letter=dephase_letter)
    M_norm = float(np.linalg.norm(M))
    palindromic = M_norm < _PALINDROME_TOL
    out['palindrome'] = {
        'M_norm': M_norm,
        'sigma_gamma': sigma_gamma,
        'palindromic': palindromic,
    }

    # Angle 3 (numerically): Π-decomposition M_sym / M_anti
    Pi = build_pi_full(N, dephase_letter=dephase_letter)
    Pi_inv = Pi.conj().T
    M_pi_conj = Pi @ M @ Pi_inv
    M_sym = 0.5 * (M + M_pi_conj)
    M_anti = 0.5 * (M - M_pi_conj)
    M_sym_norm = float(np.linalg.norm(M_sym))
    M_anti_norm = float(np.linalg.norm(M_anti))
    if M_norm > _DIVZERO_FLOOR:
        anti_fraction = (M_anti_norm ** 2) / (M_norm ** 2)
    else:
        anti_fraction = 0.0
    out['pi_decomposition'] = {
        'M_sym_norm': M_sym_norm,
        'M_anti_norm': M_anti_norm,
        'anti_fraction': anti_fraction,
    }

    # Angle 6: Bidirectional Taktgeber
    out['taktgeber'] = {
        'gamma_0': chain.gamma_0,
        'gamma_uniform': True,
        'palindrome_center': sigma_gamma,
        'bidirectional': palindromic,
    }

    # Angle 2 (structural) + numerical (optional): K-partnership
    k_expected = (chain.H_type == 'xy')
    out['k_partnership'] = {
        'expected_to_hold': k_expected,
        'h_type': chain.H_type,
        'topology': chain.topology,
    }
    if verify_k_at is not None:
        from .handshake import verify_k_partnership
        out['k_partnership']['verification'] = verify_k_partnership(
            chain, k=verify_k_at
        )

    # Conditional: Polarity (Angle 6 connection — polarity-axis input gate)
    if rho_0 is not None:
        out['polarity'] = polarity_diagnostic(rho_0, N=N)

    # Conditional: Klein-inheritance (Angle 5 connection — F77-Klein layer)
    if terms is not None:
        klein_set = {klein_index(tuple(t)) for t in terms}
        dissipator_klein = klein_index((dephase_letter,))
        out['klein_inheritance'] = {
            'klein_set': klein_set,
            'dissipator_resonance_cell': dissipator_klein,
            'h_in_resonance': dissipator_klein in klein_set,
            'klein_homogeneous': len(klein_set) == 1,
        }

    # Six-angle status list. 'kind' classifies how visibility is established:
    #   'numerical'     — verifiable via a tolerance test (e.g. palindromic)
    #   'structural'    — observable computed but no occlusion verdict
    #                     (the norms describe the bridge; they don't break it)
    #   'algebraic'     — algebraic identity, true by construction
    #   'meta'          — meta-structural reading, not a numerical claim
    angles = [
        {'name': 'F1_palindrome',
         'kind': 'numerical',
         'visible': palindromic},
        {'name': 'Handshake_idempotence',
         'kind': 'algebraic',
         'visible': True},
        {'name': 'Channel_not_memory_pi_decomp',
         'kind': 'structural',
         'visible': True},
        {'name': 'One_system_two_indices',
         'kind': 'structural',
         'visible': True},
        {'name': 'Algebra_is_inheritance',
         'kind': 'meta',
         'visible': True},
        {'name': 'Bidirectional_Taktgeber',
         'kind': 'numerical',
         'visible': palindromic},
    ]
    out['angles'] = angles

    occluded = [a['name'] for a in angles
                if a['kind'] == 'numerical' and not a['visible']]
    if not k_expected:
        occluded.append('K_partnership_structural_expected')
    out['occluded_at'] = occluded
    out['bridge_visible'] = len(occluded) == 0

    return out
