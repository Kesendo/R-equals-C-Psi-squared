"""Tests for sine_slater (the Python mirror of XyJordanWignerModes + JwSlaterPairBasis).

Pins the mode-matrix formula by hand literals, the orthonormality/dispersion spec, the Slater
antisymmetry and coincidence-zero, and the ((N+1)/2)^r raw norm law the F89 cross-triple and
y_zero proofs lean on. Parity pins vs the frozen legacy copies live in the committed root gate
`simulations/sine_slater_parity.py`.
"""
from __future__ import annotations

import sys
from itertools import combinations
from math import comb, pi, sin, sqrt
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_mode_matrix_hand_literals():
    # N=3, n=4: raw entries sin(pi*k*(j+1)/4); k=1,j=0 -> sin(pi/4)=sqrt(2)/2, k=2,j=1 -> sin(pi)=0
    u = fw.sine_mode_matrix(3, normalized=False)
    assert u.shape == (3, 3)
    assert abs(u[0, 0] - sqrt(2) / 2) < 1e-15
    assert abs(u[1, 1]) < 1e-15                      # sin(pi) = 0
    assert abs(u[0, 1] - 1.0) < 1e-15                # sin(pi/2) = 1
    # normalized carries sqrt(2/(N+1)) per row (the C# XyJordanWignerModes program)
    un = fw.sine_mode_matrix(3, normalized=True)
    assert np.allclose(un, sqrt(2.0 / 4) * u, atol=1e-16)


def test_normalized_rows_are_orthonormal_raw_rows_carry_n_over_2():
    for n_sites in (4, 5, 8):
        un = fw.sine_mode_matrix(n_sites, normalized=True)
        assert np.allclose(un @ un.T, np.eye(n_sites), atol=1e-12)
        ur = fw.sine_mode_matrix(n_sites, normalized=False)
        assert np.allclose(ur @ ur.T, ((n_sites + 1) / 2.0) * np.eye(n_sites), atol=1e-12)


def test_dispersion_matches_the_hop_spectrum():
    # eps_k = 2*J*cos(pi*k/(N+1)) is exactly the weight_block_hop w=1 spectrum at J=1
    for n_sites in (4, 7):
        eps = np.sort(fw.sine_dispersion(n_sites))
        hop_spec = np.sort(np.linalg.eigvalsh(fw.weight_block_hop(n_sites, 1)))
        assert np.allclose(eps, hop_spec, atol=1e-12)


def test_slater_antisymmetry_and_coincidence_zero():
    n_sites = 6
    u = fw.sine_mode_matrix(n_sites, normalized=False)
    tau = (1, 3, 4)
    d = fw.slater_det(u, tau, (0, 2, 5))
    assert abs(fw.slater_det(u, tau, (2, 0, 5)) + d) < 1e-12   # site swap flips sign
    assert abs(fw.slater_det(u, (3, 1, 4), (0, 2, 5)) + d) < 1e-12  # mode swap flips sign
    assert abs(fw.slater_det(u, tau, (2, 2, 5))) < 1e-12       # coincidence
    assert fw.slater_det(u, (), ()) == 1.0                     # empty rank (C# convention)


def test_raw_norm_law_is_n_over_2_cubed():
    # the (n/2)^3 law of the F89 cross-triple / y_zero proofs, all mode triples
    for n_sites in (5, 7):
        want = fw.slater_norm_sq_law(n_sites, 3)
        assert want == ((n_sites + 1) / 2.0) ** 3
        for tau in combinations(range(1, n_sites + 1), 3):
            v = fw.slater_vector(n_sites, tau, normalized=False)
            assert abs(np.dot(v, v) - want) < 1e-9 * want
        # normalized vector is unit-norm
        v1 = fw.slater_vector(n_sites, (1, 2, 3), normalized=True)
        assert abs(np.dot(v1, v1) - 1.0) < 1e-12


def test_slater_vector_enumeration_order_and_length():
    n_sites, tau = 6, (2, 4, 5)
    v = fw.slater_vector(n_sites, tau, normalized=False)
    assert len(v) == comb(n_sites, 3)
    u = fw.sine_mode_matrix(n_sites, normalized=False)
    triples = list(combinations(range(n_sites), 3))
    for i in (0, 7, len(triples) - 1):
        assert v[i] == fw.slater_det(u, tau, triples[i])
