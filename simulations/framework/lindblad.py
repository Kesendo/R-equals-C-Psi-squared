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
    if not np.allclose(H, H.conj().T):
        raise ValueError("Hamiltonian H must be Hermitian.")
    d = H.shape[0]
    N = int(round(math.log2(d)))
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for l, gamma in enumerate(gamma_l):
        if gamma == 0:
            continue
        Zl = site_op(N, l, 'Z')
        L = L + gamma * (np.kron(Zl, Zl.conj()) - np.kron(Id, Id))
    return L


def lindbladian_z_plus_t1(H, gamma_l, gamma_t1_l):
    """Z-dephasing + T1 amplitude damping.

    L(ρ) = -i[H, ρ] + Σ_l γ_l · (Z_l ρ Z_l − ρ)
                    + Σ_l γ^{T1}_l · (σ⁻_l ρ σ⁺_l − ½{σ⁺_l σ⁻_l, ρ})

    σ⁻ = (X − iY)/2 = lowering operator (|1⟩→|0⟩). With γ^{T1}_l = 0 reduces to
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


def palindrome_residual(L, Sigma_gamma, N):
    """Compute M = Π·L·Π⁻¹ + L + 2Σγ·I in Pauli-string basis.

    Returns the 4^N × 4^N residual matrix. For 'truly' Hamiltonians,
    ‖M‖ ≈ 0 to floating-point precision.
    """
    from .symmetry import build_pi_full  # delayed import (Π lives in symmetry)
    Pi = build_pi_full(N)
    M = _vec_to_pauli_basis_transform(N)
    L_pauli = (M.conj().T @ L @ M) / (2 ** N)
    Pi_inv = np.linalg.inv(Pi)
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
