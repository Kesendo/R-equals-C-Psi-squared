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
import pytest

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


# ------------------------------------------------- G5: the per-site gamma PROFILE
#
# Why the axis exists: at uniform gamma the AT rate depends only on n_diff, so every coherence of
# a sector shares it. That shared value is the hypothesis under the repo's "the sector sits at
# Re = -2*gamma" family and under the edge-normality lemma, and a builder that can only take a
# scalar cannot express its violation even deliberately. These pins are entry-EXACT: the profile
# values below are dyadic, so their sums carry no rounding and == is the right comparison.

DYADIC = [0.5, 0.25, 1.5, 0.125, 2.0, 0.75, 1.25]


def _profile(n):
    """Non-uniform, non-palindromic (so a site reversal changes it) and dyadic (so sums are exact)."""
    return np.array(DYADIC[:n])


def test_uniform_profile_equals_the_scalar_route_where_the_sum_is_EXACT():
    # The scalar and the constant-profile routes are the same real number by two arithmetics:
    # one multiply (gamma * n_diff) against a repeated sum. Where the sum is exact they are the
    # same double, and that is what is pinned here: dyadic rates at any n_diff, and gamma in
    # {0, 1} where every partial sum is an integer. The legs deliberately include a large n_diff
    # (up to 6) so the gate can actually fail.
    for n, legs in ((4, (0, 1)), (5, (1, 2)), (5, (2, 2)), (6, (1, 2)), (6, (1, 5)), (6, (3, 3))):
        for g in (1.0, 0.0, 0.5, 0.25, 0.125):
            a_scalar, k_scalar = fw.weight_block_pencil(n, *legs, gamma=g)
            a_prof, k_prof = fw.weight_block_pencil(n, *legs, gamma=np.full(n, g))
            assert np.array_equal(a_scalar, a_prof), f"n={n} legs={legs} gamma={g}"
            assert np.array_equal(k_scalar, k_prof)


def test_the_two_routes_to_a_uniform_rate_can_differ_in_the_last_bit():
    # Read the deviation rather than tolerate it. A NON-dyadic scalar reaches the same real
    # number by a route the value does not contain, so the two doubles may differ by one step;
    # this records the measured witness so nobody re-derives it as a bug. If the two routes ever
    # agree exactly here, that is a finding about the arithmetic and not a defect: read this note
    # before repairing anything.
    n, legs, g = 6, (1, 5), 0.3
    a_scalar, _ = fw.weight_block_pencil(n, *legs, gamma=g)
    a_prof, _ = fw.weight_block_pencil(n, *legs, gamma=np.full(n, g))
    assert not np.array_equal(a_scalar, a_prof)
    worst = int(np.argmax(np.abs(a_scalar - a_prof)))
    assert (a_scalar[worst], a_prof[worst]) == (-3.5999999999999996, -3.6)
    # and the whole disagreement is one representable step, nowhere else
    assert np.all(np.abs(a_scalar - a_prof) <= np.spacing(np.abs(a_scalar)))


def test_the_profile_agrees_with_the_full_liouvillian_only_when_reversed():
    # The cross-convention gate, the one this module had only in prose. `fw.lindbladian_z_dephasing`
    # builds the full 4^N L with site 0 as the LEFTMOST Kronecker factor (site l = bit n-1-l),
    # while this module is little-endian (site s = bit s). Reading a block off that L and comparing
    # it here is the obvious next move for a consumer, and unreversed it silently lands each rate
    # on the mirror site. Run at H = 0, so both objects ARE their dephasing diagonals: no
    # eigensolver, no Hamiltonian convention, no tolerance.
    n = 4
    g = _profile(n)
    d = 2 ** n
    l_full = fw.lindbladian_z_dephasing(np.zeros((d, d), dtype=complex), g[::-1])
    for w_ket, w_bra in ((0, 1), (1, 2)):
        kets = fw.weight_block_configs(n, w_ket)
        bras = fw.weight_block_configs(n, w_bra)
        flat = [k * d + b for k in kets for b in bras]
        a, _ = fw.weight_block_pencil(n, w_ket, w_bra, gamma=g)
        assert np.array_equal(a, np.real([l_full[f, f] for f in flat]))
        # negative control: the same array unreversed must DISAGREE, or this gate proves nothing
        # (_profile is non-palindromic exactly so that it can)
        l_wrong = fw.lindbladian_z_dephasing(np.zeros((d, d), dtype=complex), g)
        assert not np.array_equal(a, np.real([l_wrong[f, f] for f in flat]))
        # and why it stays hidden: at uniform gamma the reversal is invisible
        l_flat = fw.lindbladian_z_dephasing(np.zeros((d, d), dtype=complex), np.ones(n))
        a_flat, _ = fw.weight_block_pencil(n, w_ket, w_bra, gamma=1.0)
        assert np.array_equal(a_flat, np.real([l_flat[f, f] for f in flat]))


def test_the_profile_is_minus_twice_the_sum_over_disagreeing_sites():
    n, w_ket, w_bra = 5, 1, 2
    g = _profile(n)
    a, _ = fw.weight_block_pencil(n, w_ket, w_bra, gamma=g)
    kets = fw.weight_block_configs(n, w_ket)
    bras = fw.weight_block_configs(n, w_bra)
    expected = np.array([-2.0 * sum(g[s] for s in range(n) if ((k ^ b) >> s) & 1)
                         for k in kets for b in bras])
    assert np.array_equal(a, expected)          # site s is BIT s, this module's little-endian book


def test_the_profile_never_touches_the_hop():
    n = 5
    a_uniform, k_uniform = fw.weight_block_pencil(n, 1, 2, gamma=1.0)
    a_prof, k_prof = fw.weight_block_pencil(n, 1, 2, gamma=_profile(n))
    assert np.array_equal(k_uniform, k_prof)
    # a checker that reports nothing has proved nothing: the profile must have reached A at all
    assert not np.array_equal(a_uniform, a_prof)


def test_under_a_profile_the_rate_stops_being_a_function_of_n_diff_alone():
    # The physics the axis exists for, as a measurement rather than as prose.
    n = 5
    g = _profile(n)
    kets, bras = fw.weight_block_configs(n, 1), fw.weight_block_configs(n, 2)
    pairs = [(k, b) for k in kets for b in bras if bin(k ^ b).count("1") == 1]
    uniform = {fw.weight_block_disagreement_sum(n, 1.0, k, b) for k, b in pairs}
    profiled = {fw.weight_block_disagreement_sum(n, g, k, b) for k, b in pairs}
    assert len(uniform) == 1            # the collapse: one shared rate at fixed n_diff
    assert len(profiled) > 1            # the spread: what the fence is about


def test_the_profile_reaches_the_assembled_block_and_a_wrong_length_raises():
    n = 5
    g = _profile(n)
    l_uniform = fw.weight_block_build(n, 1, 2, 0.7, gamma=1.0, delta=0.85)
    l_prof = fw.weight_block_build(n, 1, 2, 0.7, gamma=g, delta=0.85)
    off = ~np.eye(l_prof.shape[0], dtype=bool)
    assert np.array_equal(l_uniform[off], l_prof[off])          # dephasing is diagonal
    assert not np.array_equal(np.diag(l_uniform), np.diag(l_prof))
    # the rate is the REAL part; the delta*ZZ frequency is imaginary and profile-free
    assert np.array_equal(np.diag(l_uniform).imag, np.diag(l_prof).imag)
    for bad in (np.ones(n - 1), np.ones((n, n))):
        with pytest.raises(ValueError):
            fw.weight_block_pencil(n, 1, 2, gamma=bad)
