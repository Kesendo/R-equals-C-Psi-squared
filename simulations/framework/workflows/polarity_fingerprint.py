"""Joint F87 × F112 × F113 polarity-axis fingerprint workflow.

Combines three orthogonal diagnostics on the bit_b Z₂-axis into one fingerprint:
  - F87 trichotomy: truly / soft / hard from spec(L)-palindromy
  - F112 verdict: BALANCED / BROKEN from M_anti's Π +i / −i Frobenius split
  - F113 prediction: closed-form magnitude of the Z-drive × amplitude-damping break

For a Lindblad system on `chain` with Hamiltonian = sum of `terms` plus
optional dissipator (Z-dephasing γ_z + σ⁻ amplitude damping γ_t1 + σ⁺
pumping γ_pump), returns a dict capturing all three readings plus their
joint interpretation.

The F113 prediction requires identifying single-site Z-drive coefficients
in H; if `z_drive_omegas_per_site` is provided (one ω per site), F113's
closed form applies and we additionally extract γ_T1 via inversion.

Cockpit-discipline recurring question this answers: "for this Lindblad
system, where does it sit on the polarity axis?": useful for any new
Hamiltonian + dissipator combination being analyzed for hardware
relevance.

Returns a dict with keys:
  'f87_class':              str  ('truly' | 'soft' | 'hard')
  'f112_asymmetry':         float (signed ‖M_+1/2‖² − ‖M_−1/2‖²)
  'f112_rel_asymmetry':     float (|asym| / ‖M‖²)
  'f112_M_norm_sq':         float (‖M‖²)
  'f112_verdict':           str  ('BALANCED' | 'near-BALANCED' | 'BROKEN')
  'in_f112_typed_scope':    bool (Hermitian H + bit_b-homogeneous c)
  'h_bit_b_homogeneous':    bool (PauliHamiltonian check on terms)
  'c_bit_b_homogeneous':    bool (true iff dissipator has no σ⁻/σ⁺)
  'f113_applies':           bool (Z-drive omegas provided + Z-drive structure detected)
  'f113_predicted':         float | None (closed-form asymmetry prediction)
  'f113_extracted_gamma_t1': float | None (γ_T1 from F113 inversion of measured asym)
  'reading':                str (human-readable summary)
"""
from __future__ import annotations

from typing import Iterable, List, Optional, Sequence

import numpy as np

from ..diagnostics.f77_trichotomy import classify_pauli_pair
from ..diagnostics.polarity_coordinates import polarity_coordinates
from ..pauli_hamiltonian import PauliHamiltonian


def polarity_fingerprint(
    chain,
    terms: Sequence,
    gamma_z: Optional[float] = None,
    gamma_t1: Optional[float] = None,
    gamma_pump: Optional[float] = None,
    z_drive_omegas_per_site: Optional[Sequence[float]] = None,
    tol: float = 1e-10,
) -> dict:
    """Joint F87 × F112 × F113 polarity-axis fingerprint.

    Args:
        chain: ChainSystem providing N and the bond graph.
        terms: list of Pauli letter tuples for the Hamiltonian (bilinear or k-body).
        gamma_z: per-site Z-dephasing rate (defaults to chain.gamma_0).
        gamma_t1: optional σ⁻ amplitude-damping rate per site (standard physics
            convention; σ⁻ = [[0, 1], [0, 0]] lowering). If non-zero, c is
            bit_b-mixed and F112 typed scope is violated.
        gamma_pump: optional σ⁺ pumping rate per site.
        z_drive_omegas_per_site: optional per-site Z-drive amplitudes ω_l for
            H = Σ_l (ω_l/2)·Z_l. Only needed for F113 prediction; if None,
            F113 fields in the result are None.
        tol: relative-asymmetry threshold for BALANCED vs BROKEN verdict.

    Returns:
        dict with keys documented at module level.
    """
    if not terms:
        raise ValueError("terms must be non-empty for a meaningful fingerprint")

    # F87 trichotomy classification
    f87 = classify_pauli_pair(chain, terms)

    # F112 polarity coordinates on the framework's predicted L
    pol = polarity_coordinates(
        chain, terms,
        gamma_z=gamma_z,
        gamma_t1=gamma_t1,
        gamma_pump=gamma_pump,
    )
    asym = float(pol['asymmetry'])
    M_sq = float(pol['norm_sq']['M'])
    rel_asym = abs(asym) / max(M_sq, 1e-15)
    if rel_asym < tol:
        f112_verdict = 'BALANCED'
    elif rel_asym < 1e-6:
        f112_verdict = 'near-BALANCED'
    else:
        f112_verdict = 'BROKEN'

    # bit_b homogeneity checks
    H = PauliHamiltonian.from_letter_tuples(terms, chain_length=chain.N)
    h_bit_b_homog = H.is_bit_b_homogeneous
    # c is bit_b-homogeneous iff no σ⁻ / σ⁺ dissipators (Z-dephasing alone is
    # single-Pauli Z, trivially bit_b-homogeneous; σ⁻ = (X − iY)/2 is mixed)
    c_bit_b_homog = (gamma_t1 is None or gamma_t1 == 0) and \
                    (gamma_pump is None or gamma_pump == 0)
    in_typed_scope = h_bit_b_homog and c_bit_b_homog

    # F113 prediction (if applicable)
    f113_predicted = None
    f113_extracted_gt1 = None
    f113_applies = False
    if z_drive_omegas_per_site is not None:
        omegas = list(z_drive_omegas_per_site)
        if len(omegas) != chain.N:
            raise ValueError(
                f"z_drive_omegas_per_site length {len(omegas)} != chain.N {chain.N}"
            )
        gt1 = gamma_t1 if gamma_t1 is not None else 0.0
        gp = gamma_pump if gamma_pump is not None else 0.0
        # F113: asymmetry = (4^N / 2) · Σ_l ω_l · (γ_pump − γ_T1)
        # (uniform-rate dissipator assumed; per-site rates would extend trivially)
        prefactor = (4 ** chain.N) / 2.0
        net_rate_diff = float(gp - gt1)
        sum_omegas = float(sum(omegas))
        f113_predicted = prefactor * sum_omegas * net_rate_diff
        f113_applies = True
        # F113 inversion: γ_T1 = −asym / ((N/2)·4^N·ω_avg) for uniform-site case
        if abs(sum_omegas) > 1e-15 and gp == 0.0:
            f113_extracted_gt1 = -asym / (prefactor * sum_omegas)

    # Human-readable reading
    parts = [f"F87 = {f87}", f"F112 = {f112_verdict}"]
    parts.append(f"in_typed_scope = {in_typed_scope}")
    if f113_applies:
        parts.append(f"F113 predicted = {f113_predicted:+.4e}")
        if f113_extracted_gt1 is not None:
            parts.append(f"F113-extracted γ_T1 = {f113_extracted_gt1:+.5f}/μs")
    reading = "; ".join(parts)

    return {
        'f87_class': f87,
        'f112_asymmetry': asym,
        'f112_rel_asymmetry': rel_asym,
        'f112_M_norm_sq': M_sq,
        'f112_verdict': f112_verdict,
        'in_f112_typed_scope': in_typed_scope,
        'h_bit_b_homogeneous': h_bit_b_homog,
        'c_bit_b_homogeneous': c_bit_b_homog,
        'f113_applies': f113_applies,
        'f113_predicted': f113_predicted,
        'f113_extracted_gamma_t1': f113_extracted_gt1,
        'reading': reading,
    }
