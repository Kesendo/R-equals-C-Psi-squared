"""Tests for the Receiver entity (F71 signature, validation, normalization)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_receiver_f71_zero_state():
    psi = np.zeros(2 ** 5, dtype=complex)
    psi[0] = 1.0
    r = fw.Receiver(psi)
    assert r.f71_class == +1
    sig = r.signature()
    assert sig['f71_eigenvalue'] == +1
    assert sig['bond_block_balanced'] is True  # N=5 → 2+2 balanced
    assert 'capacity-optimal' in sig['prediction']


def test_receiver_bonding_mode_k2_at_n5_is_f71_minus():
    """ψ_2 at N=5 is F71-antisymmetric per (-1)^(k+1) parity rule."""
    N, k = 5, 2
    psi = np.zeros(2 ** N, dtype=complex)
    norm = np.sqrt(2.0 / (N + 1))
    for i in range(N):
        psi[2 ** (N - 1 - i)] = norm * np.sin(np.pi * k * (i + 1) / (N + 1))
    psi /= np.linalg.norm(psi)
    r = fw.Receiver(psi)
    assert r.f71_class == -1


def test_receiver_at_n6_predicts_suboptimal():
    """At N=6 with unbalanced 3+2 split, F71-eigenstate is capacity-suboptimal."""
    psi = np.zeros(2 ** 6, dtype=complex)
    psi[0] = 1.0  # |0>^6, F71-symmetric
    r = fw.Receiver(psi)
    assert r.f71_class == +1
    sig = r.signature()
    assert sig['bond_block_dims'] == (3, 2)
    assert sig['bond_block_balanced'] is False
    assert 'capacity-suboptimal' in sig['prediction']


def test_receiver_rejects_2d_input():
    """Receiver must reject density matrices passed where psi is expected."""
    rho = np.eye(8, dtype=complex) / 8.0  # 8x8 max-mixed
    with pytest.raises(ValueError, match="1D state vector"):
        fw.Receiver(rho)


def test_receiver_rejects_unnormalized_psi():
    psi = np.array([1.0, 1.0, 0.0, 0.0], dtype=complex)  # norm sqrt(2)
    with pytest.raises(ValueError, match="normalized"):
        fw.Receiver(psi)


def test_receiver_from_psi_unnormalized():
    psi = np.array([1.0, 1.0, 0.0, 0.0], dtype=complex)
    r = fw.Receiver.from_psi_unnormalized(psi)
    assert np.isclose(np.linalg.norm(r.psi), 1.0)
    assert r.N == 2
