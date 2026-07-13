"""Sine-mode Slater primitives: the Python mirror of the C# canonicals
`compute/RCPsiSquared.Core/F86/JordanWigner/XyJordanWignerModes.cs` (the sine mode-matrix and
dispersion) and `compute/RCPsiSquared.Core/BlockSpectrum/JordanWigner/JwSlaterPairBasis.cs`
(the Slater determinant of sine modes).

The open-chain XY single-particle modes are u_k(j) = sin(pi*k*(j+1)/(N+1)) on interior sites
j = 0..N-1 (k = 1..N), with dispersion eps_k = 2*J*cos(pi*k/(N+1)); multi-particle sectors are
their Slater determinants. The F89 seed/twinning/cross-triple proofs use the rank-3 lift
D_tau(z1,z2,z3) = det[u_k(z_a)] (mode triple tau) with the L2 law ||D_tau||^2 = ((N+1)/2)^3 for
the RAW sine matrix. Three committed engines of record carry their own transcription
(`cross_triple_orthogonality.py` umat/slater/Ggrid, `y_zero_and_level_law.py` modes/slater/lift,
`resonant_n_twinning.py` a closure inside check_Y_is_zero); they stay frozen and are pinned by
`simulations/sine_slater_parity.py`.

Conventions (= the C# canonical; the axes the three transcriptions differ on):

- **Storage: interior only.** U[k-1, j], k = 1..N, j = 0..N-1 (site j, wall sites are NOT
  stored; the sine vanishes there identically). The legacy copies differ: cross_triple stores
  wall-padded columns x = 0..N+1 (walls zeroed), y_zero stores z = -1..N and shifts by +1 inside
  its slater, resonant stores interior like this module.
- **Normalization:** `normalized=True` (default) bakes sqrt(2/(N+1)) into each row, exactly the
  C# `XyJordanWignerModes` matrix; then the mode rows are orthonormal and every rank-r Slater
  vector over all interior r-subsets has unit L2 norm (in exact arithmetic). `normalized=False`
  gives the RAW sine the three legacy copies use, with the norm-squared law ((N+1)/2)^r carried
  separately (`slater_norm_sq_law`).
- **Float honesty:** unlike the integer-valued pencil (`weight_coherence_block`), these entries
  are transcendental; different but mathematically equal float programs (e.g. `k*x*pi/n` vs
  `pi*k*(j+1)/n`) differ in the last ULPs. Parity with the legacy copies is therefore pinned to
  tight tolerances, not entry-exact; this module's own program mirrors the C# operation order
  `norm * sin(pi * k * (j+1) / (N+1))`.
- The gauged ordering-sector lift (the (-1)^b factor on the -6 rung) is basis-dependent and
  stays with its consumers; this module provides the mode matrix and determinants only.

Tests: simulations/framework/tests/primitives/test_sine_slater.py.
Parity pins vs the frozen legacy copies: simulations/sine_slater_parity.py.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np

__all__ = [
    "sine_mode_matrix",
    "sine_dispersion",
    "slater_det",
    "slater_vector",
    "slater_norm_sq_law",
]


def sine_mode_matrix(n_sites: int, normalized: bool = True) -> np.ndarray:
    """U[k-1, j] = norm * sin(pi*k*(j+1)/(N+1)), k = 1..N, j = 0..N-1 (interior sites);
    norm = sqrt(2/(N+1)) if normalized (the C# XyJordanWignerModes matrix) else 1 (raw sine)."""
    n = n_sites + 1
    norm = np.sqrt(2.0 / n) if normalized else 1.0
    return np.array([[norm * np.sin(np.pi * k * (j + 1) / n) for j in range(n_sites)]
                     for k in range(1, n_sites + 1)])


def sine_dispersion(n_sites: int, j_coupling: float = 1.0) -> np.ndarray:
    """eps_k = 2*J*cos(pi*k/(N+1)), k = 1..N (the C# XyJordanWignerModes dispersion)."""
    n = n_sites + 1
    return np.array([2.0 * j_coupling * np.cos(np.pi * k / n) for k in range(1, n_sites + 1)])


def slater_det(u_matrix: np.ndarray, modes, sites) -> float:
    """det[ U[mode_i - 1, site_j] ] over the given 1-indexed mode set and 0-indexed interior
    sites (the C# JwSlaterPairBasis.SlaterDeterminant, general rank). Antisymmetric in sites and
    in modes; zero on coincidences."""
    modes = list(modes)
    sites = list(sites)
    if len(modes) != len(sites):
        raise ValueError("mode set and site set must have equal size")
    if not modes:
        return 1.0
    return float(np.linalg.det(np.array([[u_matrix[k - 1, z] for z in sites] for k in modes])))


def slater_vector(n_sites: int, modes, normalized: bool = True) -> np.ndarray:
    """The full rank-r Slater vector: slater_det over ALL ascending interior r-subsets of
    sites, r = len(modes), in combinations order. With normalized=True its L2 norm is 1 in
    exact arithmetic; with raw sine its norm-squared is slater_norm_sq_law(n_sites, r)."""
    u_matrix = sine_mode_matrix(n_sites, normalized)
    r = len(list(modes))
    return np.array([slater_det(u_matrix, modes, z)
                     for z in combinations(range(n_sites), r)])


def slater_norm_sq_law(n_sites: int, rank: int = 3) -> float:
    """||D_tau||^2 = ((N+1)/2)^rank for the RAW sine matrix (the (n/2)^3 law of the F89
    cross-triple / y_zero proofs at rank 3); equals 1 for the normalized matrix."""
    return ((n_sites + 1) / 2.0) ** rank
