"""Smoke tests for the lean framework package.

Run with: pytest simulations/framework/tests/

Each test verifies one cockpit-level behavior end-to-end. They cover the
hot paths: trichotomy classification, F71 signature, Confirmations lookup,
and the full Lebensader cockpit_panel.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import framework as fw


# ----------------------------------------------------------------------
# ChainSystem
# ----------------------------------------------------------------------

def test_chainsystem_n5_chain_invariants():
    chain = fw.ChainSystem(N=5)
    assert chain.N == 5
    assert chain.B == 4
    assert chain.D2 == 14  # 1²+2²+2²+2²+1² = 14
    assert chain.degrees == [1, 2, 2, 2, 1]


def test_chainsystem_classify_known_pauli_pairs():
    chain = fw.ChainSystem(N=5)
    assert chain.classify_pauli_pair([('X', 'X'), ('Y', 'Y')]) == 'truly'
    assert chain.classify_pauli_pair([('Y', 'Z'), ('Z', 'Y')]) == 'soft'
    assert chain.classify_pauli_pair([('I', 'X'), ('I', 'Z')]) == 'hard'


def test_chainsystem_predict_residual_norm_squared():
    chain = fw.ChainSystem(N=5)
    # main: (N-1) · 4^(N-2) = 4 · 64 = 256
    assert chain.predict_residual_norm_squared(1.0, 'main') == 256.0
    # single_body at chain N=5: D2/2 · 4^(N-2) = 14/2 · 64 = 7·64 = 448
    assert chain.predict_residual_norm_squared(1.0, 'single_body') == 448.0


def test_chainsystem_topology_invariants():
    """Verify B and D2 for ring, star, K_N at N=5."""
    chain_chain = fw.ChainSystem(N=5, topology='chain')
    assert (chain_chain.B, chain_chain.D2) == (4, 14)

    chain_ring = fw.ChainSystem(N=5, topology='ring')
    assert chain_ring.B == 5
    assert chain_ring.D2 == 5 * 4  # all degree 2

    chain_star = fw.ChainSystem(N=5, topology='star')
    assert chain_star.B == 4
    # hub deg 4, 4 leaves deg 1: 16 + 4·1 = 20
    assert chain_star.D2 == 20

    chain_kn = fw.ChainSystem(N=5, topology='complete')
    assert chain_kn.B == 10  # C(5,2)
    assert chain_kn.D2 == 5 * 16  # all degree N-1=4


# ----------------------------------------------------------------------
# Receiver
# ----------------------------------------------------------------------

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


# ----------------------------------------------------------------------
# Confirmations
# ----------------------------------------------------------------------

def test_confirmations_has_seven_entries():
    names = fw.Confirmations.list_names()
    assert len(names) == 7
    assert 'palindrome_trichotomy' in names
    assert 'lebensader_skeleton_trace_decoupling' in names


def test_confirmations_lookup_palindrome_trichotomy():
    e = fw.Confirmations.lookup('palindrome_trichotomy')
    assert e['date'] == '2026-04-26'
    assert e['machine'] == 'ibm_marrakesh'
    assert e['job_id'] == 'd7mjnjjaq2pc73a1pk4g'
    assert e['measured_value']['delta_soft_minus_truly'] == -0.722


def test_confirmations_unknown_raises():
    with pytest.raises(KeyError):
        fw.Confirmations.lookup('does_not_exist')


def test_confirmations_by_machine():
    marrakesh = fw.Confirmations.by_machine('ibm_marrakesh')
    kingston = fw.Confirmations.by_machine('ibm_kingston')
    assert len(marrakesh) >= 4
    assert len(kingston) >= 3


# ----------------------------------------------------------------------
# Cockpit panel
# ----------------------------------------------------------------------

def test_cockpit_panel_yzzy_t1_drop_28_at_n3():
    """Reproduces the EQ-030 hardware-confirmed drop=28 for YZ+ZY +T1=0.005 at N=3."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(plus, minus), plus)
    r = fw.Receiver(psi, chain=chain)
    result = chain.cockpit_panel(r, terms=[('Y', 'Z'), ('Z', 'Y')],
                                  gamma_t1=0.005, t_max=5.0, dt=0.01)
    skeleton = result['lebensader']['skeleton']
    assert skeleton['drop'] == 28
    assert result['lebensader']['rating'].startswith('collapsed')


def test_cockpit_panel_truly_no_drop_pure_z():
    """XX+YY (truly) under pure-Z: no skeleton drop expected."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(plus, minus), plus)
    r = fw.Receiver(psi, chain=chain)
    result = chain.cockpit_panel(r, terms=[('X', 'X'), ('Y', 'Y')],
                                  t_max=5.0, dt=0.01)
    skeleton = result['lebensader']['skeleton']
    assert skeleton['drop'] == 0


# ----------------------------------------------------------------------
# Scaling formulas
# ----------------------------------------------------------------------

def test_palindrome_residual_norm_ratio_squared_n3_n4():
    # main: 4·k/(k-1) at k=3 → 6
    assert fw.palindrome_residual_norm_ratio_squared(3, 4, 'main') == 6.0
    # single-body: 4·(2k-1)/(2k-3) at k=3 → 20/3
    assert abs(fw.palindrome_residual_norm_ratio_squared(3, 4, 'single_body') - 20 / 3) < 1e-12


def test_palindrome_residual_norm_squared_factor_graph_topologies():
    # Chain N=5: B=4, D2=14
    assert fw.palindrome_residual_norm_squared_factor_graph(5, 4, 14, 'main') == 4 * 64
    # Ring N=5: B=5, D2=20
    assert fw.palindrome_residual_norm_squared_factor_graph(5, 5, 20, 'main') == 5 * 64
    assert fw.palindrome_residual_norm_squared_factor_graph(5, 5, 20, 'single_body') == 10 * 64


# ----------------------------------------------------------------------
# F71 / chain mirror
# ----------------------------------------------------------------------

def test_chain_mirror_state_involutory():
    for N in [3, 4, 5, 6]:
        R = fw.chain_mirror_state(N)
        np.testing.assert_allclose(R @ R, np.eye(2 ** N), atol=1e-12)


def test_bond_mirror_basis_dimensions():
    """Verify the parity-tied dimensional structure."""
    for N, expected_sym, expected_asym in [(3, 1, 1), (4, 2, 1), (5, 2, 2),
                                             (6, 3, 2), (7, 3, 3), (8, 4, 3)]:
        sym, asym = fw.bond_mirror_basis(N)
        assert (len(sym), len(asym)) == (expected_sym, expected_asym), \
            f"N={N}: got ({len(sym)}, {len(asym)}), expected ({expected_sym}, {expected_asym})"
