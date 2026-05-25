"""Three-way polarity decomposition of M into {−1/2, 0, +1/2} coordinates.

Refinement of F81: F81 splits M = M_sym + M_anti by Π-conjugation parity
(eigenvalues ±1 of the linear map X ↦ Π·X·Π⁻¹). Π is order-4 on Liouville
space (Π⁴ = I) so the full Π-eigenvalue spectrum is {+1, −1, +i, −i}.
The +1 / −1 eigenspaces together form M_sym (Π²-even), and the +i / −i
eigenspaces together form M_anti (Π²-odd). This primitive refines the
±1 sub-split into the explicit +i / −i projections, giving the typed
polarity triple {−1/2, 0, +1/2} at d=2:

    M_zero       = (M + Π·M·Π⁻¹) / 2                      (0-axis, F81 M_sym)
    M_plus_half  = (M_anti − i · Π·M_anti·Π⁻¹) / 2        (Π eigenvalue +i, +1/2)
    M_minus_half = (M_anti + i · Π·M_anti·Π⁻¹) / 2        (Π eigenvalue −i, −1/2)

where M_anti = (M − Π·M·Π⁻¹) / 2 is the F81 antisymmetric part.

Frobenius-orthogonal invariant:

    ‖M‖² = ‖M_zero‖² + ‖M_plus_half‖² + ‖M_minus_half‖²

Connection to F81:
    F81 M_sym  = M_zero
    F81 M_anti = M_plus_half + M_minus_half (further split by Π ±i eigenvalue)

Working hypothesis (to be tested empirically by Task B):
    Hermitian H + pure Z-dephasing → ‖M_plus_half‖² = ‖M_minus_half‖²
    T1 cooling-only (γ_↓ ≠ γ_↑) → measurable asymmetry, F81 violation per F84

Outcome (Task B+C, 2026-05-25): Hermitian-H balance CONFIRMED across all six bilinear
H families; T1 asymmetry hypothesis REFUTED (asymmetry = 0.0 bit-exact across the 8 H
families x 3 dissipator settings tested). F84's F81-violation IS measurable but does
not split asymmetrically between +i and -i Π-eigenvalues. Structural reading: bra-ket
exchange symmetry of any Lindbladian (Hermitian-superoperator) M. Reflection doc at
reflections/POLARITY_COORDINATES.md.
"""
from __future__ import annotations

import numpy as np

from ..symmetry import build_pi_full
from .f81_pi_decomposition import pi_decompose_M


def polarity_coordinates(chain, terms, gamma_z=None, gamma_t1=None, gamma_pump=None, strict=None):
    """Three-way polarity decomposition of M = Π·L·Π⁻¹ + L + 2Σγ·I.

    Refines F81's binary sym/anti split into three orthogonal components:

        M_zero       = M_sym = (M + Π·M·Π⁻¹) / 2            (0-axis, Π²-symmetric)
        M_plus_half  = (M_anti − i·Π·M_anti·Π⁻¹) / 2        (+1/2 polarity, Π eigenvalue +i)
        M_minus_half = (M_anti + i·Π·M_anti·Π⁻¹) / 2        (−1/2 polarity, Π eigenvalue −i)

    Frobenius-orthogonal: ‖M‖² = ‖M_zero‖² + ‖M_plus_half‖² + ‖M_minus_half‖².

    The ±i projections are the standard Π-eigenvalue projectors restricted
    to the Π²-odd subspace (where Π acts with eigenvalues ±i). Π is unitary,
    so Π⁻¹ = Π†.

    Args:
        chain: ChainSystem providing N and the bond graph.
        terms: list of Pauli letter tuples; bilinear (a, b) or k-body
            (a, b, c, ...). Passed through to pi_decompose_M unchanged.
        gamma_z: per-site Z-dephasing rate (uniform if scalar; defaults to chain.gamma_0).
        gamma_t1: optional per-site T1 cooling (σ⁻ amplitude damping).
        gamma_pump: optional per-site T1 heating (σ⁺ amplitude damping).
        strict: forwarded to pi_decompose_M; if True, raises when F81
            violation > 1e-7. Defaults to True for pure Z-dephasing,
            False when any non-Z dissipator is given.

    Returns:
        dict with keys:
            'M':                   full 4^N × 4^N residual in Pauli basis.
            'M_zero':              0-axis component (Π²-symmetric, = F81 M_sym).
            'M_plus_half':         +1/2 polarity component (Π eigenvalue +i).
            'M_minus_half':        −1/2 polarity component (Π eigenvalue −i).
            'norm_sq':             dict of Frobenius norms² for M / M_zero / M_plus_half / M_minus_half.
            'asymmetry':           float ‖M_plus_half‖² − ‖M_minus_half‖² (zero for Hermitian H + pure Z-deph).
            'orthogonality_residual': float |‖M‖² − (‖M_zero‖² + ‖M_plus_half‖² + ‖M_minus_half‖²)|
                                      (machine precision when the invariant holds).
    """
    f81 = pi_decompose_M(
        chain, terms,
        gamma_z=gamma_z, gamma_t1=gamma_t1, gamma_pump=gamma_pump,
        strict=strict,
    )
    M = f81['M']
    M_sym = f81['M_sym']
    M_anti = f81['M_anti']

    # Reconstruct Π the same way pi_decompose_M does (build_pi_full from symmetry).
    # Π is unitary (signed permutation), so Π⁻¹ = Π†.
    Pi = build_pi_full(chain.N)
    Pi_inv = Pi.conj().T

    Pi_M_anti_Pi_inv = Pi @ M_anti @ Pi_inv

    M_plus_half = (M_anti - 1j * Pi_M_anti_Pi_inv) / 2
    M_minus_half = (M_anti + 1j * Pi_M_anti_Pi_inv) / 2

    M_zero = M_sym  # F81 M_sym is the 0-axis component by definition.

    norm_sq_M = float(np.sum(np.abs(M) ** 2))
    norm_sq_zero = float(np.sum(np.abs(M_zero) ** 2))
    norm_sq_plus = float(np.sum(np.abs(M_plus_half) ** 2))
    norm_sq_minus = float(np.sum(np.abs(M_minus_half) ** 2))

    orthogonality_residual = float(
        abs(norm_sq_M - (norm_sq_zero + norm_sq_plus + norm_sq_minus))
    )
    asymmetry = float(norm_sq_plus - norm_sq_minus)

    return {
        'M': M,
        'M_zero': M_zero,
        'M_plus_half': M_plus_half,
        'M_minus_half': M_minus_half,
        'norm_sq': {
            'M': norm_sq_M,
            'M_zero': norm_sq_zero,
            'M_plus_half': norm_sq_plus,
            'M_minus_half': norm_sq_minus,
        },
        'asymmetry': asymmetry,
        'orthogonality_residual': orthogonality_residual,
    }
