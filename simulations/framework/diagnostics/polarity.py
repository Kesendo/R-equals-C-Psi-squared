"""Polarity-axis diagnostic: read the +0 / 0 / −0 structure of a state.

The polarity layer (THE_POLARITY_LAYER.md) has three distinguishable
locations along the d=0 axis: +0 (X=+1, reflective end), 0 (X=0, boundary),
−0 (X=−1, reflective end). Reflection happens at the ends; 0 is the
boundary that gets crossed. The qubit is the natural window: ⟨X_i⟩ per
site is the coordinate along −0 → 0 → +0.

Public API:
  polarity_diagnostic(rho_or_psi, N=None)
"""
from __future__ import annotations

import numpy as np

from ..pauli import site_paulis, to_density_matrix


def polarity_diagnostic(rho_or_psi, N=None):
    """Read the polarity-axis content of a state.

    Args:
        rho_or_psi: 2^N × 2^N density matrix OR length-2^N state vector.
        N: chain length (inferred from input shape if not given).

    Returns:
        dict with:
          'polarity_axis': array length N, ⟨X_i⟩ ∈ [−1, +1] per site.
              +1 = at +0 end; −1 = at −0 end; 0 = on the boundary line.
          'off_axis': array length N, √(⟨Y_i⟩² + ⟨Z_i⟩²) per site.
              Distance from the polarity axis; 0 = on-axis, 1 = pure off-axis.
          'distance_to_plus_zero': array length N, 1 − ⟨X_i⟩.
              How far each site is from the +0 reflective end.
          'distance_to_minus_zero': array length N, 1 + ⟨X_i⟩.
              How far each site is from the −0 reflective end.
          'on_boundary': boolean array length N, True where ⟨X_i⟩ ≈ 0.
              Marks sites sitting on the d=0 boundary itself.
          'aggregate_polarity': float, mean(|polarity_axis|) — how strongly
              the state sits at the polarity ends overall (1 = pure ±0).
          'on_axis_fraction': float, 1 − mean(off_axis²) — what fraction of
              the per-site Bloch length lives on the X-axis.
          'site_blochs': array N×3, full (⟨X⟩, ⟨Y⟩, ⟨Z⟩) Bloch vectors per site.
    """
    rho, N = to_density_matrix(rho_or_psi, N)
    if 2 ** N != rho.shape[0]:
        raise ValueError(f"shape {rho.shape} inconsistent with N={N}")

    site_blochs = np.zeros((N, 3))
    for i, (Xi, Yi, Zi) in enumerate(site_paulis(N)):
        site_blochs[i, 0] = float(np.real(np.trace(Xi @ rho)))
        site_blochs[i, 1] = float(np.real(np.trace(Yi @ rho)))
        site_blochs[i, 2] = float(np.real(np.trace(Zi @ rho)))

    polarity_axis = site_blochs[:, 0]
    off_axis = np.sqrt(site_blochs[:, 1] ** 2 + site_blochs[:, 2] ** 2)
    distance_to_plus_zero = 1.0 - polarity_axis
    distance_to_minus_zero = 1.0 + polarity_axis
    on_boundary = np.abs(polarity_axis) < 1e-10
    aggregate_polarity = float(np.mean(np.abs(polarity_axis)))
    on_axis_fraction = 1.0 - float(np.mean(off_axis ** 2))

    return {
        'polarity_axis': polarity_axis,
        'off_axis': off_axis,
        'distance_to_plus_zero': distance_to_plus_zero,
        'distance_to_minus_zero': distance_to_minus_zero,
        'on_boundary': on_boundary,
        'aggregate_polarity': aggregate_polarity,
        'on_axis_fraction': on_axis_fraction,
        'site_blochs': site_blochs,
    }
