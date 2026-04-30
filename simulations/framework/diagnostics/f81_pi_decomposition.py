"""F81 Π-decomposition of M into symmetric/antisymmetric parts."""
from __future__ import annotations

import numpy as np

from ..lindblad import (
    lindbladian_general,
    lindbladian_z_dephasing,
    lindbladian_z_plus_t1,
    palindrome_residual,
)
from ..pauli import (
    _build_bilinear,
    _build_kbody_chain,
    _vec_to_pauli_basis_transform,
    site_op,
)
from ..symmetry import (
    _pauli_tuple_is_truly,
    _pauli_tuple_pi2_class,
    build_pi_full,
)


def pi_decompose_M(chain, terms, gamma_z=None, gamma_t1=None, gamma_pump=None, strict=None):
    """F81: Decompose M under Π-conjugation into symmetric/antisymmetric parts.

    Recurring question: for a given 2-bilinear Hamiltonian H, what is the
    Π-conjugation behavior of M = Π·L·Π⁻¹ + L + 2Σγ·I? Is M its own Π-
    conjugate, or what is the explicit shift between M and Π·M·Π⁻¹?

    Theorem F81 (proved in PROOF_F81_PI_CONJUGATION_OF_M):

        Π · M · Π⁻¹ = M − 2 · L_{H_odd}    (Z-dephasing only)

    where H = H_even + H_odd is the Π²-parity decomposition of H (Π²-odd
    Pauli bilinears have bit_b(P)+bit_b(Q) ≡ 1 mod 2). For Π²-even-only H
    (truly + YZ-type non-truly), L_{H_odd} = 0 and Π·M·Π⁻¹ = M; M is its
    own Π-conjugate. For Π²-odd-containing H, M and Π·M·Π⁻¹ differ by
    exactly the explicit shift 2·L_{H_odd}.

    The Π-decomposition of M is:

        M_sym  = (M + Π·M·Π⁻¹) / 2 = Π·L·Π⁻¹ + L_diss + L_{H_even} + 2Σγ·I
        M_anti = (M − Π·M·Π⁻¹) / 2 = L_{H_odd}    (when F81 holds)

    with M_sym and M_anti Frobenius-orthogonal: ‖M‖² = ‖M_sym‖² + ‖M_anti‖².

    For non-Z dissipators (T1 amplitude damping in particular), F81 is no
    longer exact; the identity residual ‖M_anti − L_{H_odd}‖_F is positive
    and quantifies the non-Π²-symmetric content of the dissipator. This
    makes the F81 violation a quantitative diagnostic for non-Z noise on
    hardware. With T1 enabled, this method returns the violation as
    'f81_violation' instead of raising.

    F84 (proved in PROOF_F84_AMPLITUDE_DAMPING) extends F82 to thermal
    amplitude damping (cooling γ_↓ = gamma_t1 + heating γ_↑ = gamma_pump):

        f81_violation = √(Σ_l (γ_↓_l − γ_↑_l)²) · 2^(N−1)

    For γ_↓ = γ_↑ (detailed balance / thermal equilibrium): violation = 0
    even though both channels are active. The F81 violation isolates the
    net cooling rate, equivalent to the vacuum-fluctuation contribution.

    Args:
        terms: list of (a, b) Pauli letter tuples; H = J·Σ_l Σ_terms (a_l b_(l+1)).
        gamma_z: per-site Z-dephasing rate (uniform if scalar; defaults to chain.gamma_0).
        gamma_t1: optional per-site T1 cooling rates (σ⁻ amplitude damping). None = no cooling.
        gamma_pump: optional per-site T1 heating rates (σ⁺ amplitude damping). None = no heating.
            For thermal photon bath at temperature T: γ_↓ = γ_0·(n_th+1), γ_↑ = γ_0·n_th.
        strict: if True, raises when F81 violation > 1e-7. Defaults to True for
            pure Z-dephasing, False when any non-Z dissipator (gamma_t1 or
            gamma_pump) is given.

    Returns:
        dict with keys:
            'M':            full 4^N × 4^N residual in Pauli basis.
            'M_sym':        Π-symmetric component (= M − L_{H_odd} when F81 holds).
            'M_anti':       Π-antisymmetric component.
            'L_H_odd':      reference operator -i[H_odd, ·] of the Π²-odd part of H.
            'f81_violation': float ‖M_anti − L_{H_odd}‖_F (zero for Z-only).
            'norm_sq':      dict of 'M', 'M_sym', 'M_anti', 'L_H_odd' Frobenius norms².

    Raises (only when strict=True): RuntimeError if F81 violation > 1e-7.
    """
    gz = gamma_z if gamma_z is not None else chain.gamma_0
    gamma_l = [gz] * chain.N if np.isscalar(gz) else list(gz)
    Sigma_gamma = sum(gamma_l)
    if gamma_t1 is not None:
        gt1_l = [gamma_t1] * chain.N if np.isscalar(gamma_t1) else list(gamma_t1)
    else:
        gt1_l = [0.0] * chain.N
    if gamma_pump is not None:
        gpump_l = [gamma_pump] * chain.N if np.isscalar(gamma_pump) else list(gamma_pump)
    else:
        gpump_l = [0.0] * chain.N
    any_t1 = any(g != 0 for g in gt1_l)
    any_pump = any(g != 0 for g in gpump_l)
    any_amplitude_damping = any_t1 or any_pump

    if strict is None:
        strict = not any_amplitude_damping

    # Classify terms: F85 generalization to k-body. Use _pauli_tuple_pi2_class.
    all_terms_kbody = [tuple(t) + (chain.J,) for t in terms]
    bilinear_odd = []
    for term in terms:
        letters = tuple(term)
        if _pauli_tuple_is_truly(letters):
            continue
        if 'I' in letters:
            continue
        cls = _pauli_tuple_pi2_class(letters)
        if cls == 'pi2_odd':
            bilinear_odd.append(letters + (chain.J,))

    # Build H_full per term. 2-body terms use the bond graph
    # (chain/ring/star/K_N via F49); k-body terms use chain sliding-window
    # (F85 chain-only). If mixed body counts are present, build them
    # separately and add to preserve non-chain topology for the 2-body part.
    d_full = 2 ** chain.N
    H_full = np.zeros((d_full, d_full), dtype=complex)
    two_body_terms_raw = [t for t in terms if len(t) == 2]
    kbody_terms_raw = [t for t in terms if len(t) > 2]
    if two_body_terms_raw:
        two_body_all = [(t[0], t[1], chain.J) for t in two_body_terms_raw]
        H_full = H_full + _build_bilinear(chain.N, chain.bonds, two_body_all)
    if kbody_terms_raw:
        kbody_all_with_coeff = [tuple(t) + (chain.J,) for t in kbody_terms_raw]
        H_full = H_full + _build_kbody_chain(chain.N, kbody_all_with_coeff)
    # all_two_body kept for back-compat with later L_H_odd construction
    all_two_body = (len(kbody_terms_raw) == 0)
    if any_pump:
        # Build via lindbladian_general with explicit cooling + heating channels
        d = 2 ** chain.N
        c_ops = []
        # Z-dephasing as collapse operators
        for l, gz_l in enumerate(gamma_l):
            if gz_l != 0:
                c_ops.append(np.sqrt(gz_l) * site_op(chain.N, l, 'Z'))
        # Cooling (σ⁻)
        sm = np.array([[0, 1], [0, 0]], dtype=complex)
        for l, gd_l in enumerate(gt1_l):
            if gd_l != 0:
                ops = [np.eye(2, dtype=complex)] * chain.N
                ops[l] = sm
                op_full = ops[0]
                for op in ops[1:]:
                    op_full = np.kron(op_full, op)
                c_ops.append(np.sqrt(gd_l) * op_full)
        # Heating (σ⁺)
        sp_ = np.array([[0, 0], [1, 0]], dtype=complex)
        for l, gp_l in enumerate(gpump_l):
            if gp_l != 0:
                ops = [np.eye(2, dtype=complex)] * chain.N
                ops[l] = sp_
                op_full = ops[0]
                for op in ops[1:]:
                    op_full = np.kron(op_full, op)
                c_ops.append(np.sqrt(gp_l) * op_full)
        L_full = lindbladian_general(H_full, c_ops)
    elif any_t1:
        L_full = lindbladian_z_plus_t1(H_full, gamma_l, gt1_l)
    else:
        L_full = lindbladian_z_dephasing(H_full, gamma_l)
    M = palindrome_residual(L_full, Sigma_gamma, chain.N)

    d = 2 ** chain.N
    Id_d = np.eye(d, dtype=complex)
    T = _vec_to_pauli_basis_transform(chain.N)
    if bilinear_odd:
        # Same mixed-body splitting as H_full above
        two_body_odd = [t for t in bilinear_odd if len(t) == 3]
        kbody_odd = [t for t in bilinear_odd if len(t) > 3]
        H_odd = np.zeros((d_full, d_full), dtype=complex)
        if two_body_odd:
            tbo = [(t[0], t[1], t[-1]) for t in two_body_odd]
            H_odd = H_odd + _build_bilinear(chain.N, chain.bonds, tbo)
        if kbody_odd:
            H_odd = H_odd + _build_kbody_chain(chain.N, kbody_odd)
        L_H_odd_vec = -1j * (np.kron(H_odd, Id_d) - np.kron(Id_d, H_odd.T))
        L_H_odd = (T.conj().T @ L_H_odd_vec @ T) / d
    else:
        L_H_odd = np.zeros((4 ** chain.N, 4 ** chain.N), dtype=complex)

    Pi = build_pi_full(chain.N)
    Pi_inv = np.linalg.inv(Pi)
    PiMPi = Pi @ M @ Pi_inv
    M_sym = (M + PiMPi) / 2
    M_anti = (M - PiMPi) / 2

    f81_violation = float(np.linalg.norm(M_anti - L_H_odd))
    if strict and f81_violation > 1e-7:
        raise RuntimeError(
            f"F81 identity violated: ‖M_anti − L_H_odd‖ = {f81_violation:.4e}. "
            "Set strict=False to receive the violation as output, or check for "
            "non-Z dissipator presence."
        )

    return {
        'M': M,
        'M_sym': M_sym,
        'M_anti': M_anti,
        'L_H_odd': L_H_odd,
        'f81_violation': f81_violation,
        'norm_sq': {
            'M': float(np.linalg.norm(M) ** 2),
            'M_sym': float(np.linalg.norm(M_sym) ** 2),
            'M_anti': float(np.linalg.norm(M_anti) ** 2),
            'L_H_odd': float(np.linalg.norm(L_H_odd) ** 2),
        },
    }
