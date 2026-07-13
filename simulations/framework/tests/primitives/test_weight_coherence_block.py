"""Tests for weight_coherence_block (the Python mirror of the C# WeightCoherenceBlock canonical).

Each section pins one convention axis of the module docstring: spec (hop spectrum, AT rungs,
the one-term gamma program), the q-book conversion, the leg-adjoint relation, the delta term,
and the encoding adapters. The parity pins against the four FROZEN legacy builders live in the
committed root gate `simulations/weight_coherence_parity.py` (they import sibling root scripts,
which is out of pytest's import scope here).
"""
from __future__ import annotations

import sys
from itertools import combinations
from math import comb, cos, pi
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


# ---------------------------------------------------------------- G0: spec pins

def test_configs_ascending_and_counted():
    for n in (4, 5, 7):
        for w in (1, 2, 3):
            cfgs = fw.weight_block_configs(n, w)
            assert len(cfgs) == comb(n, w)
            assert cfgs == sorted(cfgs)
            assert all(bin(m).count("1") == w for m in cfgs)


def test_hop_w1_spectrum_is_open_chain_bloch():
    for n in (4, 5, 8):
        h = fw.weight_block_hop(n, 1)
        got = np.sort(np.linalg.eigvalsh(h))
        want = np.sort([2.0 * cos(k * pi / (n + 1)) for k in range(1, n + 1)])
        assert np.allclose(got, want, atol=1e-12)


def test_pencil_at_rungs_minus2_minus6():
    for n in (4, 5):
        for legs in ((1, 2), (2, 1)):
            a, _ = fw.weight_block_pencil(n, *legs)
            assert set(np.unique(a)) == {-2.0, -6.0}


def test_gamma_is_a_live_one_term_axis():
    n = 5
    a1, k1 = fw.weight_block_pencil(n, 1, 2, gamma=1.0)
    a3, k3 = fw.weight_block_pencil(n, 1, 2, gamma=0.3)
    a0, k0 = fw.weight_block_pencil(n, 1, 2, gamma=0.0)
    kets = fw.weight_block_configs(n, 1)
    bras = fw.weight_block_configs(n, 2)
    n_diff = np.array([bin(k ^ b).count("1") for k in kets for b in bras], dtype=float)
    # one-term product program: -2 * gamma * n_diff, entry-exact
    assert np.array_equal(a3, -2.0 * 0.3 * n_diff)
    assert np.array_equal(a0, np.zeros_like(a1))
    # gamma never touches the hop
    assert np.array_equal(k1, k3) and np.array_equal(k1, k0)


# ---------------------------------------------------------------- G1: the q-book pin

def test_book_is_only_a_knob_conversion_entry_exact():
    # q_carrier = 2 * q_octic is an exact IEEE power-of-two scaling, so the two books agree
    # ENTRY-EXACT for every (q, gamma, delta), including the delta*ZZ frequency.
    n = 4
    for q in (0.7, 1.5, 2.0 + 0.3j):
        for delta in (0.0, 0.35):
            lo = fw.weight_block_build(n, 1, 2, q, gamma=0.3, delta=delta, book="octic")
            lc = fw.weight_block_build(n, 1, 2, 2 * q, gamma=0.3, delta=delta, book="carrier")
            assert np.array_equal(lo, lc)


# ---------------------------------------------------------------- G2: the leg-adjoint pin

def _tau_indices(n, w_ket, w_bra):
    """Index map of the swap |a><b| -> |b><a|: position of (b, a) in the (w_bra, w_ket) block
    for each (a, b) of the (w_ket, w_bra) block, both in ascending-mask order."""
    d_bra = comb(n, w_bra)
    d_ket = comb(n, w_ket)
    return np.array([i_bra * d_ket + i_ket
                     for i_ket in range(d_ket) for i_bra in range(d_bra)], dtype=int)


def test_leg_adjoint_is_permuted_conjugate_with_q_conjugated():
    # tau L_(u,v)(q) tau^-1 = conj(L_(v,u)(conj(q))); at real q the conj(q) is invisible.
    n = 4
    tau = _tau_indices(n, 1, 2)
    for q in (1.5, 0.8 - 0.4j):
        for delta in (0.0, 0.3):
            l12 = fw.weight_block_build(n, 1, 2, q, delta=delta)
            l21 = fw.weight_block_build(n, 2, 1, np.conj(q), delta=delta)
            assert np.array_equal(l12, np.conj(l21)[np.ix_(tau, tau)])


def test_each_block_is_complex_symmetric():
    for legs in ((1, 2), (2, 1)):
        l_mat = fw.weight_block_build(5, *legs, 1.5, delta=0.2)
        assert np.array_equal(l_mat, l_mat.T)


# ---------------------------------------------------------------- G3/G4: delta pins

def test_delta_zero_reproduces_the_pure_pencil():
    n, q = 4, 1.5
    a, k = fw.weight_block_pencil(n, 1, 2)
    assembled = np.diag(a.astype(complex)) - 2j * q * k
    assert np.array_equal(fw.weight_block_build(n, 1, 2, q, delta=0.0), assembled)


def test_zz_hand_literals():
    # zz pinned by hand so the delta-value test below is not circular in weight_block_zz:
    # n=4 bonds (0,1),(1,2),(2,3); +1 equal bits, -1 differing.
    assert fw.weight_block_zz(4, 0b0000) == 3     # all equal
    assert fw.weight_block_zz(4, 0b1111) == 3     # all equal
    assert fw.weight_block_zz(4, 0b0101) == -3    # alternating: every bond differs
    assert fw.weight_block_zz(4, 0b0011) == 1     # +1 (bits 0,1 equal) -1 (1,2) +1 (2,3)
    assert fw.weight_block_zz(4, 0b0001) == 1     # -1 (0,1) +1 (1,2) +1 (2,3)
    # even under the global bit-flip (the C# Zz property)
    for n in (4, 5):
        full = (1 << n) - 1
        for c in range(1 << n):
            assert fw.weight_block_zz(n, c) == fw.weight_block_zz(n, c ^ full)


def test_delta_value_matches_the_explicit_zz_formula():
    # pins the VALUE of the delta term (not just its shape): the C# convention
    # -1j * q_octic * delta * (zz(ket) - zz(bra)) on the diagonal, entry-exact.
    n, q, delta = 5, 1.5, 0.7
    kets = fw.weight_block_configs(n, 1)
    bras = fw.weight_block_configs(n, 2)
    zz_diff = np.array([fw.weight_block_zz(n, k) - fw.weight_block_zz(n, b)
                        for k in kets for b in bras], dtype=float)
    want = np.diag(-1j * q * delta * zz_diff)
    got = fw.weight_block_build(n, 1, 2, q, delta=delta) - fw.weight_block_build(n, 1, 2, q, delta=0.0)
    assert np.array_equal(got, want)


def test_delta_term_is_diagonal_and_leaves_the_at_rate():
    n, q = 5, 1.5
    l0 = fw.weight_block_build(n, 1, 2, q, delta=0.0)
    ld = fw.weight_block_build(n, 1, 2, q, delta=0.7)
    diff = ld - l0
    assert np.array_equal(diff, np.diag(np.diag(diff)))          # diagonal only
    assert np.array_equal(np.real(np.diag(ld)), np.real(np.diag(l0)))  # Re(diag) untouched


# ---------------------------------------------------------------- adapters

def test_combinations_order_permutation_w2_example():
    # combinations(range(4), 2) -> masks 3, 5, 9, 6, 10, 12; ascending -> 3, 5, 6, 9, 10, 12
    perm = fw.combinations_order_permutation(4, 2)
    masks = [fw.mask_of_sites(c) for c in combinations(range(4), 2)]
    assert masks == [3, 5, 9, 6, 10, 12]
    ascending = fw.weight_block_configs(4, 2)
    assert [ascending[p] for p in perm] == masks
    # identity for w <= 1
    assert list(fw.combinations_order_permutation(6, 1)) == list(range(6))


def test_mask_site_roundtrip():
    for n in (4, 6):
        for w in (1, 2, 3):
            for m in fw.weight_block_configs(n, w):
                assert fw.mask_of_sites(fw.sites_of_mask(m, n)) == m
