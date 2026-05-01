"""Lindbladian construction and palindrome-residual primitives.

Public API:
  lindbladian_general(H, c_ops)
  lindbladian_z_dephasing(H, gamma_l)
  lindbladian_z_plus_t1(H, gamma_l, gamma_t1_l)
  palindrome_residual(L, Sigma_gamma, N)
  palindrome_residual_norm_squared_factor(N, class)
  palindrome_residual_norm_squared_factor_graph(N, B, D2, class)
  palindrome_residual_norm_ratio_squared(N1, N2, class)

Note: palindrome_residual depends on build_pi_full from .symmetry. The import
is at the top of palindrome_residual to avoid circular-import risk.
"""
from __future__ import annotations

import math

import numpy as np

from .pauli import (
    site_op,
    _vec_to_pauli_basis_transform,
)


def lindbladian_general(H, c_ops):
    """General Lindbladian L = -i[H,·] + Σ_k (c_k(·)c_k† − ½{c_k†c_k, ·}) in vec form.

    Returns the 4^N × 4^N superoperator matrix using flatten('F') (column-stack)
    convention, compatible with palindrome_residual.
    """
    if not np.allclose(H, H.conj().T):
        raise ValueError("Hamiltonian H must be Hermitian.")
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for c in c_ops:
        c_dag_c = c.conj().T @ c
        L = L + (np.kron(c, c.conj())
                 - 0.5 * np.kron(c_dag_c, Id)
                 - 0.5 * np.kron(Id, c_dag_c.T))
    return L


def lindbladian_z_dephasing(H, gamma_l):
    """L = -i[H,·] + Σ_l γ_l (Z_l ρ Z_l - ρ).

    Pure Z-dephasing form for which the framework's palindrome
    Π·L·Π⁻¹ + L + 2Σγ·I = 0 holds (truly).
    """
    return lindbladian_pauli_dephasing(H, gamma_l, dephase_letter='Z')


def lindbladian_pauli_dephasing(H, gamma_l, dephase_letter='Z'):
    """L = -i[H,·] + Σ_l γ_l (P_l ρ P_l - ρ),  P ∈ {X, Y, Z}.

    Single-letter dephasing along the chosen axis. Each P_l is Hermitian
    with P² = I so D[√γ·P_l]ρ = γ·(P_l ρ P_l - ρ).

    The dissipator-resonance law (verified at N=4 k=3 over 294 Z₂³-homo-
    geneous pairs, 2026-05-01): F77-hardness lives exactly in the Klein
    cell that matches the dephase_letter's Klein index. Z's Klein index
    is (0, 1); X's is (1, 0); Y's is (1, 1).
    """
    if not np.allclose(H, H.conj().T):
        raise ValueError("Hamiltonian H must be Hermitian.")
    if dephase_letter not in ('X', 'Y', 'Z'):
        raise ValueError(
            f"dephase_letter must be 'X', 'Y', or 'Z'; got {dephase_letter!r}"
        )
    d = H.shape[0]
    N = int(round(math.log2(d)))
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for l, gamma in enumerate(gamma_l):
        if gamma == 0:
            continue
        Pl = site_op(N, l, dephase_letter)
        L = L + gamma * (np.kron(Pl, Pl.conj()) - np.kron(Id, Id))
    return L


def lindbladian_z_plus_t1(H, gamma_l, gamma_t1_l):
    """Z-dephasing + T1 amplitude damping.

    L(ρ) = -i[H, ρ] + Σ_l γ_l · (Z_l ρ Z_l − ρ)
                    + Σ_l γ^{T1}_l · (σ⁻_l ρ σ⁺_l − ½{σ⁺_l σ⁻_l, ρ})

    σ⁻ = (X + iY)/2 = [[0,1],[0,0]] = lowering operator (|1⟩→|0⟩). With γ^{T1}_l = 0 reduces to
    `lindbladian_z_dephasing`. T1 introduces palindrome-breaking; used by
    cockpit_panel to measure how Π-protected count shifts under amplitude damping.
    """
    if not np.allclose(H, H.conj().T):
        raise ValueError("Hamiltonian H must be Hermitian.")
    d = H.shape[0]
    N = int(round(math.log2(d)))
    c_ops = []
    for l, gamma in enumerate(gamma_l):
        if gamma == 0:
            continue
        c_ops.append(np.sqrt(gamma) * site_op(N, l, 'Z'))
    sigma_minus_2 = np.array([[0, 1], [0, 0]], dtype=complex)
    for l, gamma_t1 in enumerate(gamma_t1_l):
        if gamma_t1 == 0:
            continue
        ops = [np.eye(2, dtype=complex)] * N
        ops[l] = sigma_minus_2
        sigma_minus_l = ops[0]
        for op in ops[1:]:
            sigma_minus_l = np.kron(sigma_minus_l, op)
        c_ops.append(np.sqrt(gamma_t1) * sigma_minus_l)
    return lindbladian_general(H, c_ops)


def palindrome_residual(L, Sigma_gamma, N, dephase_letter='Z'):
    """Compute M = Π·L·Π⁻¹ + L + 2Σγ·I in Pauli-string basis.

    Returns the 4^N × 4^N residual matrix. For 'truly' Hamiltonians,
    ‖M‖ ≈ 0 to floating-point precision.

    Args:
        L: Liouvillian (4^N × 4^N) in vec form.
        Sigma_gamma: total dephasing rate Σ_l γ_l.
        N: chain length.
        dephase_letter: which dephasing letter the L was built for. Selects
            the right Π (Z-dephasing's Π flips bit_a, X's flips bit_b,
            Y's flips both). Default 'Z'.
    """
    from .symmetry import build_pi_full  # delayed import (Π lives in symmetry)
    Pi = build_pi_full(N, dephase_letter=dephase_letter)
    M = _vec_to_pauli_basis_transform(N)
    L_pauli = (M.conj().T @ L @ M) / (2 ** N)
    # Π is a unitary signed permutation: Π⁻¹ = Π†. Avoid np.linalg.inv.
    Pi_inv = Pi.conj().T
    return Pi @ L_pauli @ Pi_inv + L_pauli + 2 * Sigma_gamma * np.eye(4 ** N)


# ----------------------------------------------------------------------
# Palindrome-residual norm scaling — closed form (no compute)
# ----------------------------------------------------------------------

def palindrome_residual_norm_squared_factor_graph(N, B, D2, hamiltonian_class='main'):
    """F(N, G) such that ‖M(N, G)‖_F² = c_H · F(N, G) for any graph G.

      main class         ‖M(N, G)‖² = c_H · B(G) · 4^(N − 2)
      single-body class  ‖M(N, G)‖² = c_H · (D2(G) / 2) · 4^(N − 2)

    where B(G) = bond count, D2(G) = Σ_i deg_G(i)². Verified on chain, ring,
    star, K_N at N=4, 5 to machine precision. See
    experiments/OPERATOR_RIGIDITY_ACROSS_CUSP.md.
    """
    if N < 2:
        raise ValueError(f"N must be >= 2; got {N}")
    if B < 1:
        raise ValueError(f"B must be >= 1; got {B}")
    if D2 < 2 * B:
        raise ValueError(f"D2 ({D2}) inconsistent with B ({B}); D2 >= 2B = {2 * B}")
    if hamiltonian_class == 'main':
        return B * (4 ** (N - 2))
    if hamiltonian_class == 'single_body':
        return (D2 / 2) * (4 ** (N - 2))
    raise ValueError(f"hamiltonian_class must be 'main' or 'single_body'; got {hamiltonian_class!r}")


def palindrome_residual_norm_squared_factor(N, hamiltonian_class='main'):
    """Chain-specific F(N) such that ‖M(N)‖_F² = c_H · F(N).

      main class         (N − 1) · 4^(N − 2)
      single-body class  (2N − 3) · 4^(N − 2)
    """
    if N < 2:
        raise ValueError(f"N must be >= 2; got {N}")
    if hamiltonian_class == 'main':
        return (N - 1) * (4 ** (N - 2))
    if hamiltonian_class == 'single_body':
        return (2 * N - 3) * (4 ** (N - 2))
    raise ValueError(f"hamiltonian_class must be 'main' or 'single_body'; got {hamiltonian_class!r}")


def dissipator_c1_c2_from_pauli(alpha, beta, delta):
    """Closed-form (c1, c2) for a single-class dissipator from its Pauli decomposition.

    Given a Lindblad operator c = α·X + β·Y + δ·Z + ε·I (the ε·I component is
    Π-trivial and contributes nothing), the per-class palindrome-residual
    contribution at H=0 is

        ‖M(c)‖²_F = 4^(N-1) · [ c1 · Σ γ_l² + c2 · (Σ γ_l)² ]

    with closed forms (verified empirically on 25+ test classes including
    pure Paulis, σ⁻, σ⁺, real and complex superpositions, scaled operators):

        c1 = 16·|α|²·(|α|² + |β|² + |δ|²)
           + 32·|β|²·|δ|²
           + 16·Im(αβ*)² + 16·Im(αδ*)²

        c2 = 16·(|α|² + |β|² + |δ|²)²

    Structural interpretation (bit_b parity, see PROOF_BIT_B_PARITY_SYMMETRY):
        - α (X-component, bit_b-even): direct Π-conflict via |α|²·||c||² term
        - β, δ (Y, Z components, bit_b-odd): cross-conflict via |β|²·|δ|²
        - α-β and α-δ phases: extra bit_b-mixing terms via Im(αβ*)², Im(αδ*)²
        - β-δ phases: Π-trivial (no Re(βδ*) or Im(βδ*) terms)
        - Identity component: Π-trivial (drops out)

    Returns:
        (c1, c2) tuple of floats.
    """
    a2 = abs(alpha)**2
    b2 = abs(beta)**2
    d2 = abs(delta)**2
    norm_sq = a2 + b2 + d2
    im_ab = (alpha * np.conj(beta)).imag
    im_ad = (alpha * np.conj(delta)).imag
    c1 = 16 * a2 * norm_sq + 32 * b2 * d2 + 16 * im_ab**2 + 16 * im_ad**2
    c2 = 16 * norm_sq**2
    return float(c1), float(c2)


def cpsi_bell_plus(gamma_x, gamma_y, gamma_z, t):
    """F26: CΨ(t) closed form for Bell+ under arbitrary 3-axis Pauli noise.

    Bell+ = (|00⟩ + |11⟩)/√2 evolves under c_x = √γ_x · X_l, c_y = √γ_y · Y_l,
    c_z = √γ_z · Z_l dissipators on each of the 2 qubits. Then

        CΨ(t) = u · (1 + u² + v² + w²) / 12

        u = exp(−α·t),  α = 4·(γ_y + γ_z)
        v = exp(−β·t),  β = 4·(γ_x + γ_z)
        w = exp(−δ·t),  δ = 4·(γ_x + γ_y)

    Tier-1 proven (PROOF_MONOTONICITY_CPSI.md). Replaces Lindblad master
    equation solver for multi-axis Pauli noise on Bell+ states.

    F25 special case (γ_x = γ_y = 0): u = v = exp(−4γ_z·t), w = 1, recovers
    CΨ = u(1+u²)/6.

    Args:
        gamma_x, gamma_y, gamma_z: per-axis Pauli-channel rates (uniform on
            both qubits — Bell+ symmetry is preserved only under uniform noise).
        t: evolution time.

    Returns:
        CΨ(t) as float.
    """
    import math
    alpha = 4 * (gamma_y + gamma_z)
    beta = 4 * (gamma_x + gamma_z)
    delta = 4 * (gamma_x + gamma_y)
    u = math.exp(-alpha * t)
    v = math.exp(-beta * t)
    w = math.exp(-delta * t)
    return u * (1 + u * u + v * v + w * w) / 12.0


# K values (γ·t_cusp) at which CΨ first crosses 1/4, per noise channel.
# Computed directly from F26 cusp condition CΨ(t) = 1/4. K is γ-invariant.
#
# Note: F27 in ANALYTICAL_FORMULAS.md (and the table in
# PROOF_MONOTONICITY_CPSI.md) lists K_Y = 0.0867. That contradicts F26:
# pure Y gives α=4γ, β=0, δ=4γ → CΨ = u(1+u²)/6 (identical functional form
# to pure Z), so K_Y = K_Z = 0.0374. Doc typo. The true Bell+-symmetric
# pair is K_Y ↔ K_Z (both involve a non-zero α with one of {β, δ} = 0),
# while K_X is the asymmetric case (α=0, both β and δ non-zero).
CPSI_CUSP_K_PER_CHANNEL = {
    'Z':            0.0374,  # pure Z-dephasing
    'X':            0.0867,  # pure X-noise
    'Y':            0.0374,  # pure Y-noise (= K_Z by F26 functional symmetry)
    'depolarizing': 0.0440,  # γ/3 on each axis
}


def dissipator_d2_from_pauli(alpha1, beta1, delta1, alpha2, beta2, delta2):
    """Closed-form d2 cross-term constant for two single-class dissipators.

    Given two Lindblad operators c1 = α1·X+β1·Y+δ1·Z+ε1·I and c2 = α2·X+β2·Y+δ2·Z+ε2·I,
    the cross-term in ‖M(L_c1+L_c2)‖² has the bilinear form

        Cross = 4^(N-1) · [ d1 · Σ(γ1_l·γ2_l) + d2 · (Σγ1)·(Σγ2) ]

    The d2 part has a universal closed form (verified across 196 class-pair
    combinations, all topologies, σ_offset=0):

        d2 = 32 · ||c1_traceless||² · ||c2_traceless||²
           = 32 · (|α1|²+|β1|²+|δ1|²) · (|α2|²+|β2|²+|δ2|²)

    The d1 part is also bilinear in the two Pauli decompositions but the
    closed form is involved (mixes α-α phase products, β-δ Π-cross terms,
    and α-(β,δ) imaginary cross terms similar to single-class c1). Use
    `_dissipator_d1_numerical_fit` for d1 if needed.
    """
    a1 = abs(alpha1)**2 + abs(beta1)**2 + abs(delta1)**2
    a2 = abs(alpha2)**2 + abs(beta2)**2 + abs(delta2)**2
    return 32.0 * a1 * a2


# Hardware-relevant dissipator class table (verified at N=3, σ_offset=0).
# Each entry: name → ((alpha, beta, delta), c1, c2). For class K with rate
# γ_K_l per site l: per-class contribution = 4^(N-1)·[c1·Σγ²+c2·(Σγ)²].
HARDWARE_DISSIPATORS = {
    'T1':     {'pauli': (0.5, -0.5j, 0), 'c1': 3.0,  'c2': 4.0,
               'desc': "amplitude relaxation σ⁻ = (X−iY)/2"},
    'T1pump': {'pauli': (0.5, 0.5j, 0),  'c1': 3.0,  'c2': 4.0,
               'desc': "amplitude pumping σ⁺ = (X+iY)/2 (thermal excitation)"},
    'Tphi':   {'pauli': (0, 0, 1),       'c1': 0.0,  'c2': 16.0,
               'desc': "pure dephasing σ_z (Π-respecting with σ_offset=γ; "
                       "c1=c2=0 there, but c2=16 with σ_offset=0)"},
    'Xnoise': {'pauli': (1, 0, 0),       'c1': 16.0, 'c2': 16.0,
               'desc': "X-axis noise / cross-talk σ_x"},
    'Ynoise': {'pauli': (0, 1, 0),       'c1': 0.0,  'c2': 16.0,
               'desc': "Y-axis noise σ_y (Π-respecting with σ_offset=γ)"},
}

# Cross-class d1 table (numerically extracted, σ_offset=0).
# Symmetric: d1[(K1, K2)] = d1[(K2, K1)]. Self-pair gives 2·c1(K).
HARDWARE_DISSIPATOR_D1 = {
    ('T1',     'T1'):     6.0,
    ('T1',     'T1pump'): -2.0,    # cross compensates: heat bath less Π-breaking
    ('T1',     'Tphi'):   0.0,
    ('T1',     'Xnoise'): 8.0,
    ('T1',     'Ynoise'): 0.0,
    ('T1pump', 'T1pump'): 6.0,
    ('T1pump', 'Tphi'):   0.0,
    ('T1pump', 'Xnoise'): 8.0,
    ('T1pump', 'Ynoise'): 0.0,
    ('Tphi',   'Tphi'):   0.0,
    ('Tphi',   'Xnoise'): 0.0,
    ('Tphi',   'Ynoise'): 0.0,
    ('Xnoise', 'Xnoise'): 32.0,
    ('Xnoise', 'Ynoise'): 0.0,
    ('Ynoise', 'Ynoise'): 0.0,
}


def palindrome_residual_norm_ratio_squared(N1, N2, hamiltonian_class='main'):
    """Adjacent-N ratio ‖M(N2)‖²/‖M(N1)‖² (with N2 = N1+1).

      main class         4·k / (k − 1)
      single-body class  4·(2k − 1) / (2k − 3)
    """
    if N1 < 2:
        raise ValueError(f"N1 must be >= 2; got {N1}")
    if N2 != N1 + 1:
        raise NotImplementedError(f"only adjacent N supported; got N1={N1}, N2={N2}")
    k = N1
    if hamiltonian_class == 'main':
        return 4.0 * k / (k - 1)
    if hamiltonian_class == 'single_body':
        if 2 * k - 3 <= 0:
            raise ValueError(f"single-body formula requires N1 >= 2; got {N1}")
        return 4.0 * (2 * k - 1) / (2 * k - 3)
    raise ValueError(f"hamiltonian_class must be 'main' or 'single_body'; got {hamiltonian_class!r}")
