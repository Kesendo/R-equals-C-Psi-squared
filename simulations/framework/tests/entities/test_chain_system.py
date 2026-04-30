"""Tests for the ChainSystem entity (invariants, classification, residual_norm_squared)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_chainsystem_n5_chain_invariants():
    chain = fw.ChainSystem(N=5)
    assert chain.N == 5
    assert chain.B == 4
    assert chain.D2 == 14  # 1²+2²+2²+2²+1² = 14
    assert chain.degrees == [1, 2, 2, 2, 1]


def test_chainsystem_classify_known_pauli_pairs():
    chain = fw.ChainSystem(N=5)
    assert fw.classify_pauli_pair(chain, [('X', 'X'), ('Y', 'Y')]) == 'truly'
    assert fw.classify_pauli_pair(chain, [('Y', 'Z'), ('Z', 'Y')]) == 'soft'
    assert fw.classify_pauli_pair(chain, [('I', 'X'), ('I', 'Z')]) == 'hard'


def test_chainsystem_predict_residual_norm_squared():
    chain = fw.ChainSystem(N=5)
    # main: (N-1) · 4^(N-2) = 4 · 64 = 256
    assert fw.predict_residual_norm_squared(chain, 1.0, 'main') == 256.0
    # single_body at chain N=5: D2/2 · 4^(N-2) = 14/2 · 64 = 7·64 = 448
    assert fw.predict_residual_norm_squared(chain, 1.0, 'single_body') == 448.0


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


def test_chainsystem_n2_warns_about_structural_degeneracy():
    """N=2 is allowed but warns; fundamental ops still produce correct values."""
    import warnings
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        chain = fw.ChainSystem(N=2)
        user_warnings = [x for x in w if issubclass(x.category, UserWarning)]
        assert len(user_warnings) == 1
        assert 'structurally degenerate' in str(user_warnings[0].message)
    # Math still works at N=2 — fundamental vocabulary
    assert fw.classify_pauli_pair(chain, [('X','X'),('Y','Y'),('Z','Z')]) == 'truly'
    assert fw.classify_pauli_pair(chain, [('Y','Z'),('Z','Y')]) == 'soft'
    assert abs(fw.predict_residual_norm_squared_from_terms(
        chain, [('Y','Z'),('Z','Y')]) - 256.0) < 1e-6  # 2^4·2·8 = 256 at N=2


def test_chainsystem_n3plus_does_not_warn():
    import warnings
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        fw.ChainSystem(N=3)
        fw.ChainSystem(N=4)
        user_warnings = [x for x in w if issubclass(x.category, UserWarning)]
        assert len(user_warnings) == 0


def test_chainsystem_J_immutable_after_init():
    chain = fw.ChainSystem(N=3, J=1.0)
    with pytest.raises(AttributeError, match="immutable"):
        chain.J = 2.0


def test_chainsystem_gamma_0_immutable_after_init():
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    with pytest.raises(AttributeError, match="immutable"):
        chain.gamma_0 = 0.5


def test_chainsystem_topology_immutable_after_init():
    chain = fw.ChainSystem(N=3, topology='chain')
    with pytest.raises(AttributeError, match="immutable"):
        chain.topology = 'ring'


def test_residual_norm_squared_matches_classify_truly():
    """For Heisenberg (truly), ||M||^2 ≈ 0."""
    chain = fw.ChainSystem(N=4)
    norm_sq = fw.residual_norm_squared(chain, [('X','X'),('Y','Y'),('Z','Z')])
    assert norm_sq < 1e-10


def test_residual_norm_squared_matches_known_value():
    """At chain N=4, YZ+ZY soft pair has ||M||^2 = 12288 (= c_H_main(256) * F(48))."""
    chain = fw.ChainSystem(N=4)
    norm_sq = fw.residual_norm_squared(chain, [('Y','Z'),('Z','Y')])
    assert abs(norm_sq - 12288.0) < 1e-6


def test_residual_norm_squared_with_t1_matches_predict():
    """Numerical residual_norm_squared(gamma_t1=...) matches predict's closed form."""
    for N in [3, 4]:
        chain = fw.ChainSystem(N=N, gamma_0=0.05)
        for terms in [[('Y','Z'),('Z','Y')], [('I','Y'),('Y','I')],
                      [('X','X'),('Y','Y'),('Z','Z')]]:
            for gT1 in [0.001, 0.005, 0.01]:
                pred = fw.predict_residual_norm_squared_from_terms(
                    chain, terms, gamma_t1=gT1)
                num = fw.residual_norm_squared(chain, terms, gamma_t1=gT1)
                assert abs(pred - num) < 1e-6, \
                    f"N={N} terms={terms} gT1={gT1}: pred={pred} num={num}"
