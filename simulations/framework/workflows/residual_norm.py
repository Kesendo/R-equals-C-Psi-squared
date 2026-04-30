"""Numerical ‖M‖² for a Pauli-pair Hamiltonian on a chain (with optional T1)."""
from __future__ import annotations

import numpy as np

from ..pauli import _build_bilinear


def residual_norm_squared(chain, terms, J_scale=None, gamma_t1=None):
    """Numerical ‖M‖² for a Pauli-pair Hamiltonian on this chain.

    Builds H, the Lindbladian (Z-dephasing or Z+T1 if `gamma_t1` is set),
    and the palindrome residual M; returns ‖M‖²_F. For 'truly' classes
    with gamma_t1=0, ≈ 0 to floating-point precision.

    Replaces the recurring boilerplate
    (_build_bilinear → lindbladian_* → palindrome_residual → norm).

    Args:
        terms: list of (a, b) letter tuples.
        J_scale: optional override for chain.J (default uses self.J).
        gamma_t1: optional T1 amplitude-damping rates. Scalar (uniform),
                  list of length N, or None / 0 (pure Z-dephasing).
    """
    from ..lindblad import (
        lindbladian_z_dephasing,
        lindbladian_z_plus_t1,
        palindrome_residual,
    )
    J = chain.J if J_scale is None else float(J_scale)
    bilinear = [(a, b, J) for (a, b) in terms]
    H = _build_bilinear(chain.N, chain.bonds, bilinear)
    gamma_z_l = [chain.gamma_0] * chain.N
    if gamma_t1 is None:
        L = lindbladian_z_dephasing(H, gamma_z_l)
    else:
        if np.isscalar(gamma_t1):
            gamma_t1_l = [float(gamma_t1)] * chain.N
        else:
            gamma_t1_l = [float(g) for g in gamma_t1]
            if len(gamma_t1_l) != chain.N:
                raise ValueError(
                    f"gamma_t1 list length {len(gamma_t1_l)} != N {chain.N}"
                )
        L = lindbladian_z_plus_t1(H, gamma_z_l, gamma_t1_l)
    Sigma_gamma = sum(gamma_z_l)
    M = palindrome_residual(L, Sigma_gamma, chain.N)
    return float(np.linalg.norm(M) ** 2)
