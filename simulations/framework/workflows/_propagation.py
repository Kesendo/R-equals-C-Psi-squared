"""Shared propagation helpers used by multiple workflows.

These primitives factor the spectral-propagation hot loop shared by:
  - ptf.py (per-site purity trajectories)
  - bridge_dynamics.py (per-site Bloch trajectories)

Vec convention: column-stack (matches `lindbladian_general`'s
`np.kron(H, Id) - np.kron(Id, H.T)`).
"""
from __future__ import annotations

import numpy as np


def propagation_setup(L, rho_0):
    """Eigendecompose L and project ρ_0 into the eigenbasis.

    Args:
        L: 4^N × 4^N Liouvillian (vec form, column-stack).
        rho_0: 2^N × 2^N density matrix.

    Returns:
        (evals, R, R_inv, c0) where c0 = R⁻¹ · vec(ρ_0). Subsequent
        evaluation: ρ(t) = (R · diag(exp(λ·t)) · c0).reshape(d, d, 'F').
    """
    rho_vec = rho_0.flatten('F').astype(complex)
    evals, R = np.linalg.eig(L)
    R_inv = np.linalg.inv(R)
    c0 = R_inv @ rho_vec
    return evals, R, R_inv, c0


def per_site_bloch_trajectory(evals, R, c0, t_grid, site_paulis_):
    """Per-site Bloch components ⟨X_i⟩, ⟨Y_i⟩, ⟨Z_i⟩ at each t.

    Spectral propagation: ρ(t) = R · diag(exp(λ·t)) · c0, reshaped to
    d×d (column-stack), Hermitized for numerical safety, then traced
    against per-site Pauli operators.

    Args:
        evals, R, c0: from `propagation_setup`.
        t_grid: 1D array of time samples.
        site_paulis_: tuple of (X_i, Y_i, Z_i) per site (use
            `pauli.site_paulis(N)`).

    Returns:
        (N, len(t_grid), 3) array. [i, ti, k]: site i, time ti, Bloch
        component k ∈ {0=X, 1=Y, 2=Z}.
    """
    N = len(site_paulis_)
    d = R.shape[0]
    d_phys = int(round(np.sqrt(d)))
    n_t = len(t_grid)
    trajectory = np.zeros((N, n_t, 3))
    for ti, t in enumerate(t_grid):
        rho_t = (R @ (np.exp(evals * t) * c0)).reshape(d_phys, d_phys, order='F')
        rho_t = 0.5 * (rho_t + rho_t.conj().T)
        for i, (Xi, Yi, Zi) in enumerate(site_paulis_):
            trajectory[i, ti, 0] = float(np.real(np.trace(Xi @ rho_t)))
            trajectory[i, ti, 1] = float(np.real(np.trace(Yi @ rho_t)))
            trajectory[i, ti, 2] = float(np.real(np.trace(Zi @ rho_t)))
    return trajectory
