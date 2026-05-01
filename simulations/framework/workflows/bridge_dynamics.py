"""Bridge-dynamics workflow: per-site Bloch trajectories + boundary crossings.

The always-open bridge (docs/THE_BRIDGE_WAS_ALWAYS_OPEN.md) is a static
structural fact of γ₀-const palindromic chains. This workflow exposes it
geometrically as the closed parametric curve on each site's Bloch ball,
indexed by t (the PTF Taktgeber). Under γ₀-const the trajectory is
bidirectional in t — the palindromic spectrum makes forward decay and
backward recurrence two readings of the same curve.

For polarity-anchored states (|+⟩^N at +0, |−⟩^N at −0), the Hamiltonian
dynamics traces a curve from one polarity end through the X=0 boundary
toward the other end, picking up YZ-content during transit (Tom 2026-05-01:
"the boundary is not a destination; it is what gets crossed"). The
boundary-crossing events expose the YZ-content "memory" the curve carries
across the boundary.

Public API:
  bloch_trajectory(chain, rho_0, t_grid, L=None)
  polarity_crossings(trajectory, t_grid, axis_index=0, tol=1e-6)
"""
from __future__ import annotations

import numpy as np

from ..pauli import site_paulis
from .ptf import _propagation_setup


def bloch_trajectory(chain, rho_0, t_grid, L=None):
    """Per-site Bloch trajectory (⟨X_i⟩, ⟨Y_i⟩, ⟨Z_i⟩) over t.

    Spectral propagation: eigendecompose L once, evaluate ρ(t) at each t
    via R · diag(exp(λ·t)) · c0. Per site, compute the three Bloch
    components by tracing against the cached site_paulis operators.

    Args:
        chain: ChainSystem (provides N, default L).
        rho_0: 2^N × 2^N density matrix or 2^N pure-state vector.
        t_grid: 1D array of time samples.
        L: optional Liouvillian override. Default chain.L.

    Returns:
        np.ndarray of shape (N, len(t_grid), 3). Slice [i, :, 0] is ⟨X_i⟩(t)
        — the polarity-axis trajectory of site i. [i, :, 1] is ⟨Y_i⟩(t),
        [i, :, 2] is ⟨Z_i⟩(t).
    """
    if L is None:
        L = chain.L
    N = chain.N
    d = 2 ** N
    arr = np.asarray(rho_0, dtype=complex)
    if arr.ndim == 1:
        rho_mat = np.outer(arr, arr.conj())
    elif arr.ndim == 2:
        rho_mat = arr
    else:
        raise ValueError(f"rho_0 must be 1D state or 2D density matrix; got {arr.ndim}D")

    evals, R, _R_inv, c0 = _propagation_setup(L, rho_mat)
    paulis = site_paulis(N)

    n_t = len(t_grid)
    trajectory = np.zeros((N, n_t, 3))
    for ti, t in enumerate(t_grid):
        rho_t = (R @ (np.exp(evals * t) * c0)).reshape(d, d, order='F')
        rho_t = 0.5 * (rho_t + rho_t.conj().T)
        for i, (Xi, Yi, Zi) in enumerate(paulis):
            trajectory[i, ti, 0] = float(np.real(np.trace(Xi @ rho_t)))
            trajectory[i, ti, 1] = float(np.real(np.trace(Yi @ rho_t)))
            trajectory[i, ti, 2] = float(np.real(np.trace(Zi @ rho_t)))
    return trajectory


def polarity_crossings(trajectory, t_grid, axis_index=0, tol=1e-6):
    """Find moments where a Bloch component crosses zero per site.

    For each site i, scan trajectory[i, :, axis_index] for sign changes.
    Each sign change is a "boundary crossing" — the trajectory crossing
    one of the three Bloch-ball equators. For axis_index=0 (default,
    polarity X), these are crossings of the +0 ↔ −0 boundary.

    The Bloch components at the crossing carry the "memory" the curve
    transports across the boundary (Y-content, Z-content; Tom 2026-05-01).

    Args:
        trajectory: N × len(t_grid) × 3 array from `bloch_trajectory`.
        t_grid: time samples corresponding to trajectory.
        axis_index: 0 = X (polarity), 1 = Y, 2 = Z.
        tol: minimum |Bloch_axis| at the bracketing samples for the
            crossing to count (filters numerical noise around zero).

    Returns:
        list of dicts, one per crossing event:
          'site': i (0..N-1)
          't_cross': interpolated zero-crossing time
          'bloch_at_cross': (X, Y, Z) tuple at the crossing
          'direction': '+→−' or '−→+' along the chosen axis
    """
    if axis_index not in (0, 1, 2):
        raise ValueError(f"axis_index must be 0, 1, or 2; got {axis_index}")
    N, n_t, _ = trajectory.shape
    events = []
    for i in range(N):
        sig = trajectory[i, :, axis_index]
        for k in range(n_t - 1):
            if abs(sig[k]) < tol or abs(sig[k + 1]) < tol:
                continue
            if sig[k] * sig[k + 1] < 0:
                frac = sig[k] / (sig[k] - sig[k + 1])
                t_cross = t_grid[k] + frac * (t_grid[k + 1] - t_grid[k])
                bloch_at = trajectory[i, k] + frac * (trajectory[i, k + 1] - trajectory[i, k])
                direction = '+→−' if sig[k] > 0 else '−→+'
                events.append({
                    'site': i,
                    't_cross': float(t_cross),
                    'bloch_at_cross': tuple(float(b) for b in bloch_at),
                    'direction': direction,
                })
    return events
