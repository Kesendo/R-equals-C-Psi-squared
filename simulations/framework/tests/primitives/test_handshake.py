"""Tests for handshake-algebra primitives: bonding_mode_state + k_partner.

These primitives encode the structural fact (HANDSHAKE_ALGEBRA.md) that a
handshake is the algebra acting on both sides — picking the same (N, k)
on both sides IS the handshake, no exchange step.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_bonding_mode_state_normalized_single_excitation():
    """ψ_k is normalized and lives entirely in the single-excitation sector."""
    N, k = 5, 2
    psi = fw.bonding_mode_state(N, k)
    assert psi.shape == (2 ** N,)
    assert abs(np.linalg.norm(psi) - 1.0) < 1e-12
    # Vacuum (idx 0) and multi-excitation amplitudes are zero.
    assert abs(psi[0]) < 1e-15
    se_indices = {2 ** (N - 1 - j) for j in range(N)}
    for idx in range(2 ** N):
        if idx not in se_indices:
            assert abs(psi[idx]) < 1e-15, f"non-SE index {idx} has amplitude"


def test_bonding_mode_state_n5_k1_amplitude_pattern():
    """N=5, k=1 amplitudes match the F65 closed form ψ_k(j) = √(2/6) sin(π(j+1)/6)."""
    N, k = 5, 1
    psi = fw.bonding_mode_state(N, k)
    expected = np.sqrt(2.0 / 6.0) * np.array(
        [np.sin(np.pi * (j + 1) / 6.0) for j in range(N)]
    )
    actual = np.array([psi[2 ** (N - 1 - j)].real for j in range(N)])
    assert np.allclose(actual, expected, atol=1e-12)


def test_bonding_mode_state_K_partner_identity():
    """ψ_{N+1-k}(j) = (-1)^j · ψ_k(j) at every site, for every k.

    This is the F65 amplitude identity that drives K-partnership.
    """
    N = 5
    for k in range(1, N + 1):
        psi_k = fw.bonding_mode_state(N, k)
        psi_p = fw.bonding_mode_state(N, N + 1 - k)
        for j in range(N):
            idx = 2 ** (N - 1 - j)
            sign = (-1) ** j
            assert np.allclose(psi_p[idx], sign * psi_k[idx], atol=1e-12), \
                f"K-identity broken at N={N}, k={k}, site={j}"


def test_bonding_mode_state_invalid_k():
    """k outside [1, N] raises."""
    with pytest.raises(ValueError):
        fw.bonding_mode_state(5, 0)
    with pytest.raises(ValueError):
        fw.bonding_mode_state(5, 6)


def test_bonding_mode_pair_state_normalized_includes_vacuum():
    """(|vac⟩ + |ψ_k⟩)/√2 is normalized and has both vacuum and SE content."""
    N, k = 5, 2
    psi = fw.bonding_mode_pair_state(N, k)
    assert abs(np.linalg.norm(psi) - 1.0) < 1e-12
    assert abs(psi[0]) > 0.5  # vacuum component dominates / matches
    # Has SE content
    se_norm_sq = sum(abs(psi[2 ** (N - 1 - j)]) ** 2 for j in range(N))
    assert se_norm_sq > 0.4  # meaningful SE component


def test_k_partner_basic():
    assert fw.k_partner(5, 1) == 5
    assert fw.k_partner(5, 2) == 4
    assert fw.k_partner(5, 3) == 3  # self-partner for odd N midpoint
    assert fw.k_partner(6, 3) == 4  # no self-partner for even N


def test_k_partner_involution():
    """k_partner is an involution: k_partner(N, k_partner(N, k)) == k."""
    for N in (4, 5, 6, 7):
        for k in range(1, N + 1):
            assert fw.k_partner(N, fw.k_partner(N, k)) == k


def test_k_partner_invalid_k():
    with pytest.raises(ValueError):
        fw.k_partner(5, 0)
    with pytest.raises(ValueError):
        fw.k_partner(5, 6)


def test_receiver_bonding_mode_classmethod():
    """Receiver.bonding_mode wraps bonding_mode_state with chain context."""
    chain = fw.ChainSystem(N=5, H_type='xy')
    r = fw.Receiver.bonding_mode(chain, k=2)
    assert r.N == 5
    assert r.chain is chain
    # F71 class for k=2 is antisymmetric: R·ψ_k = (-1)^(k+1) ψ_k = -ψ_2
    assert r.f71_class == -1


def test_receiver_bonding_mode_with_vacuum():
    """with_vacuum=True returns the pair state, has vacuum component."""
    chain = fw.ChainSystem(N=5, H_type='xy')
    r_pair = fw.Receiver.bonding_mode(chain, k=1, with_vacuum=True)
    r_pure = fw.Receiver.bonding_mode(chain, k=1, with_vacuum=False)
    assert abs(r_pair.psi[0]) > 0.5
    assert abs(r_pure.psi[0]) < 1e-15
