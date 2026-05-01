"""Tests for polarity-axis primitives: polarity_state + polarity_diagnostic.

Structure (Tom 2026-05-01): three distinguishable locations along the
polarity axis — +0 end (X=+1), boundary (X=0), −0 end (X=−1) — with gap
spaces between each end and the boundary. Reflection happens at +0 or −0;
0 is the boundary, not the destination.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_polarity_state_uniform_plus():
    """polarity_state(N, +1) is |+⟩^N: at the +0 end, all sites X=+1."""
    psi = fw.polarity_state(3, +1)
    assert psi.shape == (8,)
    assert abs(np.linalg.norm(psi) - 1.0) < 1e-12
    # All amplitudes are +1/√8
    expected = np.ones(8, dtype=complex) / np.sqrt(8.0)
    assert np.allclose(psi, expected, atol=1e-12)


def test_polarity_state_uniform_minus():
    """polarity_state(N, -1) is |−⟩^N: at the −0 end."""
    psi = fw.polarity_state(2, -1)
    # |-,-> = (|0> - |1>) ⊗ (|0> - |1>) / 2 = (|00> - |01> - |10> + |11>) / 2
    expected = np.array([1.0, -1.0, -1.0, 1.0]) / 2.0
    assert np.allclose(psi, expected, atol=1e-12)


def test_polarity_state_x_neel():
    """polarity_state(N, [+1,-1,+1,...]) is the X-basis Néel."""
    psi = fw.polarity_state(3, [+1, -1, +1])
    # |+,-,+> via tensor product
    plus = np.array([1, 1]) / np.sqrt(2)
    minus = np.array([1, -1]) / np.sqrt(2)
    expected = np.kron(np.kron(plus, minus), plus)
    assert np.allclose(psi, expected, atol=1e-12)


def test_polarity_state_invalid_signs():
    with pytest.raises(ValueError):
        fw.polarity_state(3, [+1, -1])  # wrong length
    with pytest.raises(ValueError):
        fw.polarity_state(3, [+1, 0, -1])  # 0 not allowed


def test_polarity_diagnostic_at_plus_end():
    """A state at the +0 end has X=+1 per site, distance_to_+0 = 0."""
    psi = fw.polarity_state(3, +1)
    diag = fw.polarity_diagnostic(psi)
    assert np.allclose(diag['polarity_axis'], 1.0, atol=1e-12)
    assert np.allclose(diag['distance_to_plus_zero'], 0.0, atol=1e-12)
    assert np.allclose(diag['distance_to_minus_zero'], 2.0, atol=1e-12)
    assert not diag['on_boundary'].any()
    assert abs(diag['aggregate_polarity'] - 1.0) < 1e-12


def test_polarity_diagnostic_at_minus_end():
    """A state at the −0 end has X=−1 per site, distance_to_-0 = 0."""
    psi = fw.polarity_state(3, -1)
    diag = fw.polarity_diagnostic(psi)
    assert np.allclose(diag['polarity_axis'], -1.0, atol=1e-12)
    assert np.allclose(diag['distance_to_plus_zero'], 2.0, atol=1e-12)
    assert np.allclose(diag['distance_to_minus_zero'], 0.0, atol=1e-12)


def test_polarity_diagnostic_z_basis_state_on_boundary():
    """Z-basis state |000⟩ has ⟨X⟩=0 per site → ON the boundary."""
    psi = np.zeros(8, dtype=complex)
    psi[0] = 1.0
    diag = fw.polarity_diagnostic(psi, N=3)
    assert np.allclose(diag['polarity_axis'], 0.0, atol=1e-12)
    assert diag['on_boundary'].all()
    # Z-content is off-axis
    assert np.allclose(diag['site_blochs'][:, 2], 1.0, atol=1e-12)


def test_polarity_diagnostic_bonding_mode_on_boundary_z_only():
    """F65 bonding modes have ⟨X⟩=0 per site, ⟨Y⟩=0 per site, ⟨Z⟩≠0.

    Real-amplitude single-excitation states (the standard bonding modes) sit
    ON the X=0 boundary, AND have zero Y content because Y requires complex
    amplitudes. The off-axis content of a bonding mode is purely Z — the
    population direction. This is structural: real H + real initial state
    keeps ρ(t) real per site for all t, so ⟨Y_i⟩ ≡ 0.
    """
    chain_N = 4
    psi = fw.bonding_mode_state(chain_N, 1)
    diag = fw.polarity_diagnostic(psi, N=chain_N)
    # On the boundary along X
    assert np.allclose(diag['polarity_axis'], 0.0, atol=1e-12)
    assert diag['on_boundary'].all()
    # Y is zero per site (real amplitudes)
    assert np.allclose(diag['site_blochs'][:, 1], 0.0, atol=1e-12)
    # Z is non-trivial — the gap content lives entirely on Z
    assert (np.abs(diag['site_blochs'][:, 2]) > 0.1).any()


def test_polarity_diagnostic_accepts_density_matrix():
    """Diagnostic works on both |ψ⟩ and ρ inputs."""
    psi = fw.polarity_state(3, +1)
    rho = np.outer(psi, psi.conj())
    diag_psi = fw.polarity_diagnostic(psi)
    diag_rho = fw.polarity_diagnostic(rho)
    assert np.allclose(diag_psi['polarity_axis'], diag_rho['polarity_axis'])
    assert np.allclose(diag_psi['site_blochs'], diag_rho['site_blochs'])


def test_polarity_diagnostic_aggregate_polarity_for_neel():
    """X-Néel has |⟨X⟩|=1 per site → aggregate_polarity = 1."""
    psi = fw.polarity_state(3, [+1, -1, +1])
    diag = fw.polarity_diagnostic(psi)
    assert abs(diag['aggregate_polarity'] - 1.0) < 1e-12
