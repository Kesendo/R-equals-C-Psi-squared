"""d=0 substrate diagnostics: kernel of L, d=0/d=2 state decomposition,
and sector populations.

The framework's d²−2d=0 condition forces d=0 or d=2. d=2 is the qubit
dimension; d=0 is the substrate axis the qubit sits on. The polarity
layer (THE_POLARITY_LAYER.md) gives ±0/0 internal structure to that
axis at the d=2 projection level.

This module reads the d=0 substrate directly via the kernel of the
Liouvillian L. Eigenvectors of L with eigenvalue 0 are the stationary
modes: what does not decohere under L_D and does not rotate under L_H.
For uniform XY/Heisenberg chain + Z-dephasing they are the excitation-
sector projectors P_n (F4 in ANALYTICAL_FORMULAS.md), all living in the
{I, Z}^N Pauli sublattice (n_xy = 0, lens-immune per XOR_SPACE.md).
The natural observables on the d=0 axis are the sector populations
p_n = Tr(P_n · ρ): conserved under L, hardware-trivial to measure
(Z-basis tomography, group shots by Hamming weight).

Public API:
  stationary_modes(chain, L=None, tol=1e-9)
  d_zero_decomposition(rho, chain, tol=1e-9)
  sector_populations(rho_or_psi, N=None)
"""
from __future__ import annotations

import numpy as np

from ..pauli import _vec_to_pauli_basis_transform, to_density_matrix


def stationary_modes(chain, L=None, tol=1e-9):
    """Kernel of L — the d=0 substrate.

    Eigenvectors of L with |λ| < tol. For uniform XY/Heisenberg chain
    + Z-dephasing this matches F4's stationary-mode count: N+1 sector
    projectors P_n on {0, 1, ..., N}-excitation subspaces.

    Args:
        chain: ChainSystem.
        L: optional Liouvillian override. Default chain.L.
        tol: |λ| threshold for kernel-mode classification.

    Returns:
        dict with:
          'eigenvalues':       complex array, length n_kernel
          'right_eigvecs':     d² × n_kernel — right kernel basis
          'left_covecs':       n_kernel × d² — corresponding left covectors
          'indices':           original eig() indices
          'pauli_decomposition': n_kernel × 4^N — each mode's Pauli coefficients
                                 (vec[k] = (1/2^N)·Tr(σ_k · M_kernel_mode))
          'kernel_dimension':  n_kernel
    """
    if L is None:
        L = chain.L
    evals, R = np.linalg.eig(L)
    R_inv = np.linalg.inv(R)
    kernel_idx = np.where(np.abs(evals) < tol)[0]

    basis_transform = _vec_to_pauli_basis_transform(chain.N)
    pauli_decomp = np.zeros((len(kernel_idx), 4 ** chain.N), dtype=complex)
    for ii, k in enumerate(kernel_idx):
        pauli_decomp[ii] = (basis_transform.conj().T @ R[:, k]) / (2 ** chain.N)

    return {
        'eigenvalues': evals[kernel_idx],
        'right_eigvecs': R[:, kernel_idx],
        'left_covecs': R_inv[kernel_idx, :],
        'indices': kernel_idx,
        'pauli_decomposition': pauli_decomp,
        'kernel_dimension': int(len(kernel_idx)),
    }


def d_zero_decomposition(rho, chain, tol=1e-9):
    """Decompose ρ into d=0 (kernel of L) + d=2 (orthogonal) parts.

    ρ_d0 = lim_{t→∞} e^{Lt} ρ — the steady-state projection. What
    survives forever under L. Sits in the kernel subspace.

    ρ_d2 = ρ − ρ_d0 — the part that decoheres. Sits in the
    range-of-non-stationary-modes subspace.

    Args:
        rho: 2^N × 2^N density matrix.
        chain: ChainSystem.
        tol: kernel-mode threshold (forwarded to `stationary_modes`).

    Returns:
        dict with:
          'rho_d0':            the d=0 (steady-state) part of ρ
          'rho_d2':            the d=2 (decohering) part
          'd0_weight':         Tr(ρ_d0), what fraction lives in d=0
          'd2_norm':           ‖ρ_d2‖_F, magnitude of decohering content
          'kernel_dimension':  n_kernel
    """
    sm = stationary_modes(chain, tol=tol)
    rho_vec = rho.flatten('F').astype(complex)

    # Projector onto kernel via biorthogonal basis: P = R_kernel · W_kernel
    P_kernel = sm['right_eigvecs'] @ sm['left_covecs']

    rho_d0_vec = P_kernel @ rho_vec
    d_phys = 2 ** chain.N
    rho_d0 = rho_d0_vec.reshape(d_phys, d_phys, order='F')
    rho_d0 = 0.5 * (rho_d0 + rho_d0.conj().T)

    rho_d2 = rho - rho_d0

    return {
        'rho_d0': rho_d0,
        'rho_d2': rho_d2,
        'd0_weight': float(np.real(np.trace(rho_d0))),
        'd2_norm': float(np.linalg.norm(rho_d2)),
        'kernel_dimension': sm['kernel_dimension'],
    }


def sector_populations(rho_or_psi, N=None):
    """Sector populations p_n = Σ_{|k|=n} ρ_kk for n ∈ {0, 1, ..., N}.

    The natural observables on the d=0 axis. For uniform XY/Heisenberg
    + Z-dephasing, these are the projections of ρ onto the kernel modes
    P_n (excitation-sector projectors): p_n = Tr(P_n · ρ). Conserved by
    H = XY/Heisenberg and by Z-dephasing; drift in p_n(t) is a direct
    fingerprint of non-{Z, XY/Heisenberg} noise (T1 σ⁻ shifts mass to
    lower n; transverse fields h_x·X, h_y·Y mix sectors with Δn = ±1).

    Hardware-trivial measurement protocol: Z-basis tomography, group
    shots by Hamming weight, p_n = (#shots with n excitations) / total.
    No phase tomography needed.

    Args:
        rho_or_psi: 2^N × 2^N density matrix or 2^N pure-state vector.
        N: number of qubits (inferred from shape if None).

    Returns:
        dict with:
          'p':       (N+1)-array of real probabilities, p[n] = Tr(P_n · ρ)
          'mean_n':  ⟨n⟩ = Σ n·p_n (mean excitation count, conserved)
          'var_n':   Σ (n − ⟨n⟩)² · p_n (excitation-count variance)
          'entropy': −Σ p_n log p_n (Shannon, base e; 0·log0 := 0)
          'N':       number of qubits
    """
    rho, N = to_density_matrix(rho_or_psi, N=N)
    diag = np.real(np.diag(rho))
    d_phys = 2 ** N
    p = np.zeros(N + 1)
    for k in range(d_phys):
        n_k = bin(k).count('1')
        p[n_k] += diag[k]

    n_axis = np.arange(N + 1)
    mean_n = float(np.sum(n_axis * p))
    var_n = float(np.sum((n_axis - mean_n) ** 2 * p))
    nonzero = p > 0
    entropy = float(-np.sum(p[nonzero] * np.log(p[nonzero])))

    return {
        'p': p,
        'mean_n': mean_n,
        'var_n': var_n,
        'entropy': entropy,
        'N': N,
    }
