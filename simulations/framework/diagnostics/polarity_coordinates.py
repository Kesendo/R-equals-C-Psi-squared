"""Three-way polarity decomposition of M into {−1/2, 0, +1/2} coordinates.

Refinement of F81: F81 splits M = M_sym + M_anti by Π-conjugation parity
(eigenvalues ±1 of the linear map X ↦ Π·X·Π⁻¹). Π is order-4 on Liouville
space (Π⁴ = I), so the conjugation map A := Ad_Π has spectrum {+1, −1, +i, −i}.

CAREFUL (corrected 2026-07-29): because A has order 4, (1 ± A)/2 are NOT
eigenprojections, and the three pieces below are NOT the A-eigenspaces they were
once described as. Writing M = u + v + p + m for the true A-eigencomponents
(+1, −1, +i, −i), what this module builds is

    M_zero       = u + (1+i)p/2 + (1−i)m/2
    M_plus_half  =     (1+i)v/2 + (1−i)p/2
    M_minus_half =     (1−i)v/2 + (1+i)m/2

so M_zero is not the Π²-even part (for XY+YX at N=3 it is entirely Π²-ODD), and
the three are not mutually Frobenius-orthogonal. The Π²-even projector is
(1 + A²)/2. What survives untouched is the quantity everything downstream reads:
asymmetry = ½(‖P_{+i}M‖² − ‖P_{−i}M‖²), since the v contamination enters both
halves equally and cancels; F112's zero and F113's (4^N/2) are unaffected.

The typed polarity triple {−1/2, 0, +1/2} at d=2:

    M_zero       = (M + Π·M·Π⁻¹) / 2                      (0-axis, F81 M_sym)
    M_plus_half  = (M_anti − i · Π·M_anti·Π⁻¹) / 2        (Π eigenvalue +i, +1/2)
    M_minus_half = (M_anti + i · Π·M_anti·Π⁻¹) / 2        (Π eigenvalue −i, −1/2)

where M_anti = (M − Π·M·Π⁻¹) / 2 is the F81 antisymmetric part.

Norm-sum identity (NOT an orthogonality check):

    ‖M‖² = ‖M_zero‖² + ‖M_plus_half‖² + ‖M_minus_half‖²

This holds algebraically for ANY unitary Π and ANY M (the cross terms cancel
identically), so `orthogonality_residual` is machine zero by construction and
cannot fail. It is a wiring check on the arithmetic, not evidence about Π.

Connection to F81:
    F81 M_sym  = M_zero
    F81 M_anti = M_plus_half + M_minus_half (further split by Π ±i eigenvalue)

Working hypothesis (to be tested empirically by Task B):
    Hermitian H + pure Z-dephasing → ‖M_plus_half‖² = ‖M_minus_half‖²
    T1 cooling-only (γ_↓ ≠ γ_↑) → measurable asymmetry, F81 violation per F84

Outcome (Task B+C, 2026-05-25): Hermitian-H balance CONFIRMED across all six bilinear
H families; T1 asymmetry measured 0.0 bit-exact across the 8 H families x 3 dissipator
settings tested.

SCOPE CORRECTION (2026-07-29): that sweep covered BOND Hamiltonians only, for which
Tr(Z_l H) = 0, and the vanishing is a property of that family, NOT of Lindbladians in
general. The earlier "structural reading: bra-ket exchange symmetry of any Lindbladian"
is withdrawn: it is refuted by F113 below in this very module. Add single-site Z-drives
and the asymmetry is nonzero and exactly linear,
    asymmetry = (4^N / 2) · Σ_l ω_l · (γ_pump,l − γ_T1,l),
e.g. −2.08e-3 at N=2 and −1.248e-2 at N=3 for ω = 0.13 on EVERY site and γ_T1 = 0.001
(a drive on one site alone gives −1.04e-3 and −4.16e-3). Note the
chain-bound `polarity_coordinates(chain, terms, ...)` entry point cannot exhibit the
effect: pi_decompose_M keeps only len==2 and len>2 terms, so a single-site ('Z',) is
SILENTLY DROPPED (it returns a decomposition of the zero Hamiltonian, no error, unlike
the sister `classify_pauli_pair`, which raises for k<2). The L-bound entry point can.
Reflection doc at reflections/POLARITY_COORDINATES.md.
"""
from __future__ import annotations

import numpy as np

from ..symmetry import build_pi_full
from ..pauli import _vec_to_pauli_basis_transform
from .f81_pi_decomposition import pi_decompose_M


def polarity_coordinates(chain, terms, gamma_z=None, gamma_t1=None, gamma_pump=None, strict=None):
    """Three-way polarity decomposition of M = Π·L·Π⁻¹ + L + 2Σγ·I.

    Refines F81's binary sym/anti split into three components (NOT orthogonal, and
    NOT A-eigenspaces; see the module docstring):

        M_zero       = M_sym = (M + Π·M·Π⁻¹) / 2            (0-axis, F81 M_sym)
        M_plus_half  = (M_anti − i·Π·M_anti·Π⁻¹) / 2        (+1/2 polarity, Π eigenvalue +i)
        M_minus_half = (M_anti + i·Π·M_anti·Π⁻¹) / 2        (−1/2 polarity, Π eigenvalue −i)

    Norm-sum identity (holds for any unitary Π, so it cannot fail):
    ‖M‖² = ‖M_zero‖² + ‖M_plus_half‖² + ‖M_minus_half‖².

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
            'M_zero':              0-axis component (= F81 M_sym; NOT the Π²-even part).
            'M_plus_half':         +1/2 polarity component (Π eigenvalue +i).
            'M_minus_half':        −1/2 polarity component (Π eigenvalue −i).
            'norm_sq':             dict of Frobenius norms² for M / M_zero / M_plus_half / M_minus_half.
            'asymmetry':           float ‖M_plus_half‖² − ‖M_minus_half‖² (zero for Hermitian H + pure Z-deph).
            'orthogonality_residual': float |‖M‖² − (‖M_zero‖² + ‖M_plus_half‖² + ‖M_minus_half‖²)|
                                      (machine zero by construction; cannot fail).
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


def polarity_coordinates_from_L(L_pauli, N, sigma, Pi=None):
    """Three-way polarity decomposition starting from a Liouvillian L in Pauli basis.

    Sister primitive to <see cref="polarity_coordinates"/> that bypasses the
    framework's H + dissipator construction and accepts L directly. Useful for
    probing structurally exotic cases that pi_decompose_M doesn't construct:
      - Non-Hermitian H (complex bond coefficients; L = -i[H,·] is no longer skew-Hermitian)
      - Single-site terms (transverse fields h_l·σ_l that pi_decompose_M rejects)
      - Mixed dephase letters (γ_Z·D[Z] + γ_X·D[X] simultaneous; pi_decompose_M is single-letter)
      - Non-Lindblad dissipators (custom non-CP super-operators)

    Args:
        L_pauli: full 4^N × 4^N Liouvillian in Pauli basis (numpy complex array).
        N: chain length.
        sigma: total dephasing rate (Σγ); shifts the F1 palindrome around -σ.
        Pi: optional precomputed Π operator; defaults to build_pi_full(N).

    Returns:
        Same dict structure as polarity_coordinates, with 'M', 'M_zero',
        'M_plus_half', 'M_minus_half', 'norm_sq', 'asymmetry', 'orthogonality_residual'.

    No F81-violation check (the F81 identity is a 2-body statement about L_{H_odd};
    this entry point is for cases where that identity is expected to fail).
    """
    if Pi is None:
        Pi = build_pi_full(N)
    Pi_inv = Pi.conj().T
    d2 = 4 ** N

    Pi_L_Pi_inv = Pi @ L_pauli @ Pi_inv
    M = Pi_L_Pi_inv + L_pauli + (2.0 * sigma) * np.eye(d2, dtype=complex)

    Pi_M_Pi_inv = Pi @ M @ Pi_inv
    M_sym = (M + Pi_M_Pi_inv) / 2
    M_anti = (M - Pi_M_Pi_inv) / 2

    Pi_M_anti_Pi_inv = Pi @ M_anti @ Pi_inv
    M_plus_half = (M_anti - 1j * Pi_M_anti_Pi_inv) / 2
    M_minus_half = (M_anti + 1j * Pi_M_anti_Pi_inv) / 2

    M_zero = M_sym

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


def polarity_coordinates_from_hc(H, c_ops, gammas, N, sigma=None, Pi=None):
    """Polarity decomposition built from a standard Lindblad (H, c_ops, gammas) triple.

    Thin composition wrapper: constructs the full standard-Lindblad-channel
    Liouvillian

        L_vec = -i · (H ⊗ I − I ⊗ H^T)
              + Σ_k γ_k · [ kron(c_k, c_k^*) − (1/2)·( kron(c_k^† c_k, I)
                                                    + kron(I, (c_k^† c_k)^T) ) ]

    in vec(ρ) basis (the standard Lindblad / GKSL dissipator, trace-preserving
    for Hermitian H + arbitrary c), transforms to Pauli basis, and delegates
    to polarity_coordinates_from_L. Matches the chain-bound
    polarity_coordinates path bit-exactly for Hermitian Pauli-letter c with
    matching σ = Σ γ_k. Absorbs the build_L_standard_lindblad pattern that
    probe scripts 1, 5, 7, 9–14 hand-roll inline.

    F112 (Hermitian H + each c_k bit_b-homogeneous) predicts asymmetry = 0
    bit-exact. Asymmetry ≠ 0 here is the precise witness for non-Hermitian
    H, non-bit_b-homogeneous c, or both. To check bit_b-homogeneity of c
    when c is given as a PauliHamiltonian, use its is_bit_b_homogeneous
    property; this wrapper accepts c as raw matrices and does not perform
    the check.

    Args:
        H: Hilbert-space Hamiltonian as a 2^N × 2^N numpy complex array.
        c_ops: iterable of 2^N × 2^N collapse operators (numpy complex).
        gammas: iterable of rates matching c_ops (complex allowed for
            non-physical sweeps; standard Lindblad uses real ≥ 0).
        N: chain length.
        sigma: F1 palindrome center; defaults to sum(gammas). For uniform
            single-letter dephasing this is N · γ.
        Pi: optional precomputed Π operator.

    Returns:
        Same dict as polarity_coordinates_from_L.
    """
    c_list = list(c_ops)
    g_list = list(gammas)
    if len(c_list) != len(g_list):
        raise ValueError(f"len(c_ops)={len(c_list)} != len(gammas)={len(g_list)}")
    if sigma is None:
        sigma = float(np.real(sum(g_list)))
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L_vec = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for c, g in zip(c_list, g_list):
        c_dag_c = c.conj().T @ c
        anti = 0.5 * (np.kron(c_dag_c, Id) + np.kron(Id, c_dag_c.T))
        L_vec = L_vec + g * (np.kron(c, c.conj()) - anti)
    # The transform below is order='F' DELIBERATELY, against a row-stack L, and the
    # mismatch must not be "aligned" away.
    # The ±i projectors that split M_anti into M_plus_half / M_minus_half are exchanged
    # by the stacking twist, so this pairing decides which component is called +1/2 and
    # hence the SIGN of `asymmetry`. The pairing is mismatched on purpose, equivalent to
    # conjugating Π: either consistent pairing returns +2.08e-3 where this one returns
    # −2.08e-3 at ω=0.13, γ_T1=0.001, N=2. F113's stated direction (negative for cooling)
    # is written for THIS pairing, so changing it here means changing F113's sign too.
    # The asymmetry vanishes for bond Hamiltonians, so a check over Heisenberg/XXZ alone
    # will not catch a flip here.
    T = _vec_to_pauli_basis_transform(N, order='F')
    L_pauli = (T.conj().T @ L_vec @ T) / (2 ** N)
    return polarity_coordinates_from_L(L_pauli, N, sigma, Pi=Pi)
