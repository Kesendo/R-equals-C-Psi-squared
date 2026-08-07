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
closed form is evaluated at that declared drive, and γ_T1 is extracted by
inversion only when the declared drive is the one `terms` actually carries
(see 'f113_drive_matches_terms').

A term carrying an identity leg, such as ('Z', 'I'), is handled normally now:
`pi_decompose_M` used to skip such terms when building H_odd, leaving L_H_odd = 0
against a nonzero M_anti, so the strict F81 check raised and blamed a dissipator
that was not there. What ('Z', 'I') builds is a Z on the left site of every bond,
which is one site whenever the bonds share a left endpoint: a 2-site chain (one
bond) or a star at any N, where it is (N-1)*Z on the hub. A bare
1-body term ('Z',) is not expressible through `terms` at all and raises; for a drive
on one chosen site of a longer chain, build H and c explicitly and call
`polarity_coordinates_from_hc`.

Cockpit-discipline recurring question this answers: "for this Lindblad
system, where does it sit on the polarity axis?": useful for any new
Hamiltonian + dissipator combination being analyzed for hardware
relevance.

Returns a dict with keys:
  'f87_class':              str  ('truly' | 'soft' | 'hard')
  'f112_asymmetry':         float (signed ‖M_+1/2‖² − ‖M_−1/2‖²)
  'f112_rel_asymmetry':     float (|asym| / ‖M_anti‖², a contrast ratio in [0, 1]; 0.0 when
                            'f112_polarity_degenerate')
  'f112_M_norm_sq':         float (‖M‖²)
  'f112_M_anti_norm_sq':    float (‖M_anti‖² = ‖M_+1/2‖² + ‖M_−1/2‖², the polarity content
                            and the denominator above)
  'f112_polarity_degenerate': bool (there is no polarity content to balance: H is Π²-even
                            AND every c is bit_b-homogeneous, so both halves of M_anti
                            vanish as a theorem. Decided STRUCTURALLY from the Pauli
                            letters, with the float ‖M_anti‖² == 0.0 only as a fallback,
                            because the structural zero reaches exact 0.0 at N = 2 alone)
  'f112_verdict':           str  ('DEGENERATE' | 'BALANCED' | 'near-BALANCED' | 'BROKEN').
                            DEGENERATE means SILENT, not balanced: no verdict may be read
                            off such a row and counting it as an F112 confirmation counts
                            nothing. It outranks the ratio.
  'in_f112_typed_scope':    bool (Hermitian H + bit_b-homogeneous c)
  'h_bit_b_homogeneous':    bool (PauliHamiltonian check on terms)
  'c_bit_b_homogeneous':    bool (true iff dissipator has no σ⁻/σ⁺)
  'f113_applies':           bool (Z-drive omegas were provided; this does NOT check that
                            terms carries them, which is 'f113_drive_matches_terms')
  'f113_predicted':         float | None (closed-form asymmetry prediction at the DECLARED
                            drive, rate-weighted per site: (4^N/2)·Σ_l ω_l·(γ_pump,l − γ_T1,l))
  'f113_drive_matches_terms': bool | None (whether the declared drive equals the single-site
                            Z content of the H the measurement builds, Tr(Z_l H)/2^(N−1))
  'f113_true_z_moments':    list[float] | None (that content, so a mismatch is diagnosable)
  'f113_extracted_gamma_t1': float | None (Σ_l ω_l·γ_T1,l / Σ_l ω_l, obtained by inverting
                            the MEASURED asymmetry. None unless the declared drive matches
                            terms, Σ_l ω_l ≠ 0 and no pumping is on; `reading` names which
                            of the three failed. One asymmetry is one number: it cannot
                            resolve N per-site rates, though N independent drive profiles
                            can, one weighted combination each)
  'reading':                str (human-readable summary)
"""
from __future__ import annotations

import math
from typing import Iterable, List, Optional, Sequence

import numpy as np

from ..diagnostics.f77_trichotomy import classify_pauli_pair
from ..diagnostics.polarity_coordinates import (
    polarity_coordinates,
    is_structurally_degenerate,
)
from ..pauli_hamiltonian import PauliHamiltonian


def _true_z_moments(chain, terms) -> List[float]:
    """The ω_l the MEASUREMENT actually carries: ω_l = Tr(Z_l H) / 2^(N−1).

    Rebuilds H exactly as `pi_decompose_M` does (2-body terms over the bond graph,
    k-body over the sliding chain window, both at coefficient chain.J) so a caller's
    declared drive can be checked against it. Without that check the F113 inversion
    divides a measured asymmetry by a drive that is not in the Hamiltonian, which
    returns a plausible wrong rate rather than an error.
    """
    from ..pauli import _build_bilinear, _build_kbody_chain, site_op

    d = 2 ** chain.N
    H = np.zeros((d, d), dtype=complex)
    two_body = [t for t in terms if len(t) == 2]
    kbody = [t for t in terms if len(t) > 2]
    if two_body:
        H = H + _build_bilinear(
            chain.N, chain.bonds, [(t[0], t[1], chain.J) for t in two_body]
        )
    if kbody:
        H = H + _build_kbody_chain(
            chain.N, [tuple(t) + (chain.J,) for t in kbody]
        )
    scale = 2 ** (chain.N - 1)
    return [
        float(np.trace(site_op(chain.N, l, 'Z') @ H).real) / scale
        for l in range(chain.N)
    ]


def _materialise(value):
    """Turn an iterable input into a list once; leave None and scalars untouched.

    Everything downstream reads these inputs more than once (a shape check, the
    measurement, the F113 block), and a one-shot iterator silently empties after the
    first read.
    """
    if value is None or np.isscalar(value):
        return value
    return list(value)


def _reject_negative(values, name: str) -> None:
    """A negative rate reaches lindbladian_z_plus_t1 as sqrt(negative).

    Everything downstream then goes NaN, and the verdict branches compare NaN < tol, both
    False, so an all-NaN measurement falls through to a confident 'BROKEN'. Refuse it here
    instead, with the value named.
    """
    for l, v in enumerate(values):
        if not math.isfinite(v) or v < 0.0:
            raise ValueError(f"{name}[{l}] = {v} is not a valid rate (must be finite and >= 0)")


def _as_per_site(value, n: int, name: str) -> List[float]:
    """Broadcast a rate given as None / scalar / per-site sequence to a list of N.

    F113's closed form weights each site's moment by THAT site's net rate, so the
    per-site list is the shape the law needs; a scalar is the uniform special case.
    """
    if value is None:
        return [0.0] * n
    # np.isscalar, not isinstance(value, (int, float)): the measurement path
    # (pi_decompose_M) uses np.isscalar, and a bare isinstance rejects np.float32 /
    # np.int32 rates that the measurement accepts.
    if np.isscalar(value):
        out = [float(value)] * n
        _reject_negative(out, name)
        return out
    values = [float(v) for v in value]
    if len(values) != n:
        raise ValueError(f"{name} length {len(values)} != chain.N {n}")
    _reject_negative(values, name)
    return values


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
        gamma_z: per-site Z-dephasing rate (defaults to chain.gamma_0). It cannot move
            any returned value: the 2Σγ recentring cancels the dephasing content of M
            exactly (F1), so ‖M‖², the asymmetry and the extraction are identical for
            every gamma_z. It is also the one rate whose length is not validated.
        gamma_t1: optional σ⁻ amplitude-damping rate (standard physics convention;
            σ⁻ = [[0, 1], [0, 0]] lowering). Either a scalar, applied to every site,
            or a per-site sequence of length N. If non-zero, c is bit_b-mixed and
            F112 typed scope is violated. The per-site form is the one F113's law is
            written for; the scalar is its uniform special case.
        gamma_pump: optional σ⁺ pumping rate, scalar or per-site sequence of length N.
        z_drive_omegas_per_site: optional per-site Z-drive amplitudes ω_l for
            H = Σ_l (ω_l/2)·Z_l. Only needed for F113 prediction; if None,
            F113 fields in the result are None.
        tol: relative-asymmetry threshold for BALANCED vs BROKEN. The third verdict,
            near-BALANCED, occupies the fixed band [tol, 1e-6) and so is unreachable
            once tol ≥ 1e-6.

    Returns:
        dict with keys documented at module level.
    """
    # Materialise every iterable input ONCE, before anything reads it. A generator
    # consumed by a validation pass would reach the measurement empty, and an empty
    # rate list reads as "no T1 channel at all" while an empty drive makes the
    # drive-match check vacuously true. `terms` goes first and unconditionally: four
    # separate consumers read it below, and a generator is truthy even when empty, so
    # the non-empty guard cannot fire until it is a list.
    terms = list(terms)
    if not terms:
        raise ValueError("terms must be non-empty for a meaningful fingerprint")

    gamma_t1 = _materialise(gamma_t1)
    gamma_pump = _materialise(gamma_pump)
    omegas = (None if z_drive_omegas_per_site is None
              else [float(w) for w in z_drive_omegas_per_site])
    # Then validate the shapes, so a malformed input raises before the 4^N
    # decomposition below pays for it (tens of seconds at N = 6).
    gt1_per_site = _as_per_site(gamma_t1, chain.N, 'gamma_t1')
    gpump_per_site = _as_per_site(gamma_pump, chain.N, 'gamma_pump')
    if omegas is not None and len(omegas) != chain.N:
        raise ValueError(
            f"z_drive_omegas_per_site length {len(omegas)} != chain.N {chain.N}"
        )

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
    # Denominator is the POLARITY CONTENT ||M_anti||^2 = ||M_+||^2 + ||M_-||^2, the whole of
    # what asym is a difference of, so rel_asym is a contrast ratio in [0, 1] with no
    # tolerance constant in it. Replaces max(||M||^2, 1e-15), which was saturated wherever M
    # vanishes (F83: 'truly' H gives M = 0 identically), and which was also the wrong SCALE
    # even where it did not vanish: for Pi^2-even non-truly H such as YZ+ZY, ||M||^2 is large
    # while M_anti = L_{H_odd} is exactly zero. Mirrors
    # PolarityCoordinatesResult.RelativeAsymmetry in C#.
    M_anti_sq = float(pol['norm_sq']['M_plus_half'] + pol['norm_sq']['M_minus_half'])
    # STRUCTURAL first, float second, and the order is the whole point. The exact test reads
    # the Pauli letters and is N-independent; the float `== 0.0` below only fires at N = 2,
    # because the structurally-zero M_anti arrives as ~1e-33 at N=3 and ~1e-31 at N=5. With
    # only the float test a silent configuration reports a ratio of two noise quantities,
    # which can land anywhere: on a Pi^2-even H with per-bond-varying J at N=4 that reaches
    # rel ~ 0.29, i.e. a BROKEN verdict on the very hypothesis F112 proves to be balanced.
    structurally_silent = is_structurally_degenerate(terms, gamma_t1, gamma_pump)
    polarity_degenerate = structurally_silent or (M_anti_sq == 0.0)
    rel_asym = 0.0 if polarity_degenerate else abs(asym) / M_anti_sq
    if structurally_silent:
        # Not balanced: SILENT. There is no polarity content to be balanced or broken, so
        # neither verdict may be read off this row, and counting it as a confirmation of
        # F112 counts nothing.
        f112_verdict = 'DEGENERATE'
    elif rel_asym < tol:
        f112_verdict = 'BALANCED'
    elif rel_asym < 1e-6:
        f112_verdict = 'near-BALANCED'
    else:
        f112_verdict = 'BROKEN'

    # bit_b homogeneity checks
    H = PauliHamiltonian.from_letter_tuples(terms, chain_length=chain.N)
    h_bit_b_homog = H.is_bit_b_homogeneous
    # c is bit_b-homogeneous iff no σ⁻ / σ⁺ dissipators (Z-dephasing alone is
    # single-Pauli Z, trivially bit_b-homogeneous; σ⁻ = (X + iY)/2 is mixed)
    c_bit_b_homog = not any(gt1_per_site) and not any(gpump_per_site)
    in_typed_scope = h_bit_b_homog and c_bit_b_homog

    # F113 prediction (if applicable)
    f113_predicted = None
    f113_extracted_gt1 = None
    f113_applies = False
    f113_drive_matches_terms = None
    drive_mismatch = False
    f113_true_moments = None
    if omegas is not None:
        # F113: asymmetry = (4^N / 2) · Σ_l ω_l · (γ_pump,l − γ_T1,l), each site's
        # moment weighted by THAT site's net rate. The sum does not factorise into
        # (Σ_l ω_l) · (net rate), so a zero-total drive profile does not silence it.
        # The firing condition is Σ_l ω_l·(γ_pump,l − γ_T1,l) ≠ 0 and nothing weaker:
        # a non-uniform rate profile that is ORTHOGONAL to the drive profile gives
        # exactly zero (γ_T1 = (1,1,2,2) against ω = ±0.1 alternating, for instance).
        prefactor = (4 ** chain.N) / 2.0
        f113_predicted = prefactor * sum(
            w * (gpump_per_site[l] - gt1_per_site[l]) for l, w in enumerate(omegas)
        )
        f113_applies = True
        # The inversion divides the MEASURED asymmetry by the DECLARED drive, so it is
        # meaningless unless the two describe the same Hamiltonian. Check that first;
        # a declared drive absent from `terms` is what used to return a plausible
        # wrong rate (7.7x out) or a confident −0.0.
        f113_true_moments = _true_z_moments(chain, terms)
        # Relative to the scale of the whole drive, taken over BOTH vectors, so it
        # neither demands absolute precision at large J nor degenerates into an
        # absolute window at tiny J, where an absolute floor would accept a
        # sign-flipped drive and return a negative rate. The scale is the vector's and
        # not the component's, because a component can be a float-noise zero
        # (Tr(Z_l H) is a sum of 2^N terms and cancels to ~1e-17, not to 0.0, at
        # non-dyadic J) while its neighbours are O(1): per component the tolerance
        # there would be 1e-9·1e-17 and the correct declaration of 0.0 would be
        # refused. What is accepted is the ω_l the Hamiltonian carries, up to that
        # noise; a rounded declaration is a different drive and is refused, and
        # `f113_true_z_moments` in the result says what was expected.
        # The scale is the vector's, over both sides. One more floor is needed: when
        # EVERY true moment is float noise (a Hamiltonian with no single-site Z at all)
        # the max is itself ~1e-17, the window collapses again, and a correct all-zero
        # declaration gets refused with a message accusing the caller. Treat a moment
        # vector that is entirely noise as the zero drive it represents.
        raw_scale = max([abs(t) for t in f113_true_moments]
                        + [abs(w) for w in omegas])
        drive_scale = raw_scale if raw_scale > 1e-12 else 1.0
        f113_drive_matches_terms = all(
            abs(w - t) <= 1e-9 * drive_scale
            for w, t in zip(omegas, f113_true_moments)
        )
        # Even when it matches, one asymmetry is one number: what comes back is the
        # ratio Σ_l ω_l·γ_T1,l / Σ_l ω_l, never the N per-site rates. That ratio is a
        # weighted MEAN only when the ω_l share a sign; with mixed signs the weights
        # sum to 1 but are individually unbounded, so it can exceed every true rate or
        # turn negative. No guard for that here, because it cannot arise once the drive
        # matches: every term is built at the single coefficient chain.J, so the
        # single-site Z moments all carry sign(J), and the match above is relative, so
        # it cannot pass a drive whose sign differs from the moment's. The
        # hazard is real for anyone applying the closed form by hand, and is stated
        # where that happens (the F113 registry entry).
        sum_omegas = sum(omegas)
        if not f113_drive_matches_terms:
            drive_mismatch = True
        elif abs(sum_omegas) > 1e-15 and not any(gpump_per_site):
            f113_extracted_gt1 = -asym / (prefactor * sum_omegas)

    # Human-readable reading
    parts = [f"F87 = {f87}", f"F112 = {f112_verdict}"]
    parts.append(f"in_typed_scope = {in_typed_scope}")
    if f113_applies:
        parts.append(f"F113 predicted = {f113_predicted:+.4e}")
        if f113_extracted_gt1 is not None:
            parts.append(f"F113-extracted γ_T1 = {f113_extracted_gt1:+.5f}/μs")
        elif drive_mismatch:
            parts.append(
                "F113-extracted γ_T1 = unavailable (the declared z_drive_omegas_per_site "
                "is not the single-site Z content of terms, so the inversion would "
                "divide the measured asymmetry by a drive the Hamiltonian does not "
                "carry; pass the drive as a term, or build H and c explicitly and use "
                "polarity_coordinates_from_hc). The drive terms carries is "
                f"{[round(t, 12) for t in f113_true_moments]}"
            )
        elif any(gpump_per_site):
            parts.append(
                "F113-extracted γ_T1 = unavailable (pumping is on, so the measured "
                "asymmetry reads the net rate difference and cannot be resolved into "
                "a γ_T1)"
            )
        else:
            parts.append(
                "F113-extracted γ_T1 = unavailable (Σ_l ω_l = 0, so the inversion has "
                "no scale to divide by; the asymmetry itself is still readable)"
            )
    reading = "; ".join(parts)

    return {
        'f87_class': f87,
        'f112_asymmetry': asym,
        'f112_rel_asymmetry': rel_asym,
        'f112_M_norm_sq': M_sq,
        'f112_M_anti_norm_sq': M_anti_sq,
        'f112_polarity_degenerate': polarity_degenerate,
        'f112_verdict': f112_verdict,
        'in_f112_typed_scope': in_typed_scope,
        'h_bit_b_homogeneous': h_bit_b_homog,
        'c_bit_b_homogeneous': c_bit_b_homog,
        'f113_applies': f113_applies,
        'f113_predicted': f113_predicted,
        'f113_drive_matches_terms': f113_drive_matches_terms,
        'f113_true_z_moments': f113_true_moments,
        'f113_extracted_gamma_t1': f113_extracted_gt1,
        'reading': reading,
    }
