"""Tests for framework.diagnostics.far_element_shape. The ranks and kind patterns are pinned against the C# route
(TrianglePalindromeLocusWitness, its kernels over Q and OperatorSchmidtRank, a separate implementation; its tests
pin the same ranks and patterns), the pencils and the homogeneity groups of C and E against an independent numpy
route (explicit 8 x 8 matrices, SVD nullspace), the primitives against operators written by hand. The sheet counts
(`extra`) come from this module alone; stage E's 0 is checked by hand in its test."""
from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from framework.diagnostics.far_element_shape import (  # noqa: E402
    far_element_shape, homogeneity_group, kind_patterns, operator_schmidt_rank, pencil, sheets)
from framework.diagnostics.f158_end_count import end_count  # noqa: E402


def _add(*ds):
    out = {}
    for d in ds:
        for k, v in d.items():
            out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v != 0}


def _bond(i, j, n, J, letters='XYZ'):
    t = {}
    for c in letters:
        s = ['I'] * n; s[i] = c; s[j] = c
        t[''.join(s)] = sp.nsimplify(J)
    return t


def _field(k, n, f):
    t = {}
    for c, x in zip('XY', f):
        if x:
            s = ['I'] * n; s[k] = c
            t[''.join(s)] = sp.nsimplify(x)
    return t


def triangle(J, fd, fu, fv):
    """Heisenberg triangle on (d, u, v) = (0, 1, 2), planar fields; the jump is Z on d."""
    return _add(_bond(0, 1, 3, J[0]), _bond(0, 2, 3, J[1]), _bond(1, 2, 3, J[2]),
                _field(0, 3, fd), _field(1, 3, fu), _field(2, 3, fv))


A = triangle((1, 2, 3), (1, 0), (2, 0), (3, 0))
B = triangle((2, 2, 3), (1, 1), (3, 4), (4, 3))
C = triangle((3, 6, -2), (0, 0), (3, 4), (5, 0))
D = triangle((1, 0, 3), (1, 0), (-1, -2), (0, 2))
E = triangle((1, 2, 6), (0, 0), (8, 0), (-5, 3))
S = triangle((2, 2, 3), (1, 2), (3, 1), (3, 1))


# ---- the primitives on operators written by hand ----------------------------------------------------------

def test_schmidt_rank_of_a_product_and_of_an_entangled_pair():
    product = {'IZY': 1, 'IZZ': 2, 'XZY': 1, 'XZZ': 2}            # (I + X) ⊗ Z ⊗ (Y + 2Z)
    assert [operator_schmidt_rank(product, 3, [l]) for l in range(3)] == [1, 1, 1]
    pair = {'XXY': 1, 'ZZY': 1}                                    # (XX + ZZ) on sites 0, 1, Y on 2: not symmetric
    assert [operator_schmidt_rank(pair, 3, [l]) for l in range(3)] == [2, 2, 1]
    assert operator_schmidt_rank({'II': 1, 'XX': 1, 'YY': 1, 'ZZ': 1}, 2, [0]) == 4   # 2·SWAP


def test_a_single_string_is_homogeneous_under_every_string():
    assert len(homogeneity_group({'XYZ': 1}, 3)) == 4 ** 3


def test_kind_patterns_read_against_the_axes():
    w = {'XIZ': 1, 'YXI': 1}
    assert kind_patterns(w, 3) == ['100', '110']
    assert kind_patterns(w, 3, axes='XZZ') == ['000', '110']     # against X on site 0: X quiet, Y lit


# ---- the triangle: ranks and patterns against the C# witness ----------------------------------------------

def test_the_triangle_ranks_and_patterns_match_the_witness():
    expected = {  # (ranks across d | uv, u | dv, v | du), patterns: TrianglePalindromeLocusWitness's readings
        'A': (A, (1, 1, 1), ['111']),
        'B': (B, (1, 4, 4), ['100', '111']),
        'C': (C, (2, 4, 4), ['100', '111']),
        'D': (D, (2, 3, 2), ['101', '111']),
        'E': (E, (2, 4, 4), ['100', '101', '110', '111']),
    }
    for name, (h, ranks, patterns) in expected.items():
        shape = far_element_shape(3, h, ['ZII'])
        assert shape is not None, name
        assert shape.single_site_ranks == ranks, name
        assert shape.patterns == patterns, name
        assert shape.product == (ranks == (1, 1, 1)), name


def test_no_shape_where_the_row_does_not_pair_through_one_element():
    assert far_element_shape(3, S, ['ZII']) is None                # the singlet: far 1, near 2
    assert far_element_shape(3, triangle((1, 2, 6), (1, 0), (8, 0), (-5, 3)), ['ZII']) is None   # broken


# ---- sheets and pencils -----------------------------------------------------------------------------------

def test_c_sits_on_the_pair_parity_sheet_and_e_on_none():
    # the groups against the numpy route: C {III, IZZ, ZII, ZZZ}, E {III, ZII}; the sheet counts are this module's
    c = far_element_shape(3, C, ['ZII'])
    assert sorted(homogeneity_group(c.element, 3)) == ['III', 'IZZ', 'ZII', 'ZZZ']
    assert c.sheets.extra == 1
    e = far_element_shape(3, E, ['ZII'])
    assert sorted(homogeneity_group(e.element, 3)) == ['III', 'ZII']   # the jump alone
    assert e.sheets.extra == 0


def test_the_pencils_against_a_numpy_route():
    # split by Z_u Z_v: H_− holds the hop (XX + YY of the d bonds) and the fields on u and v. The numpy route
    # (explicit matrices, SVD) reads C paired with the same element at λ = 2 and 1/3, E broken at both; at λ = 0,
    # H_+ alone, both ends grow, C to (12, 12) and E to (8, 8), which tells H_− from H_+.
    for lam in (2, sp.Rational(1, 3)):
        ec = pencil(3, C, ['ZII'], 'IZZ', lam)
        assert (ec.near, ec.far) == (1, 1)
        ee = pencil(3, E, ['ZII'], 'IZZ', lam)
        assert (ee.near, ee.far) == (1, 0)
    ec, ee = pencil(3, C, ['ZII'], 'IZZ', 0), pencil(3, E, ['ZII'], 'IZZ', 0)
    assert ((ec.near, ec.far), (ee.near, ee.far)) == ((12, 12), (8, 8))


def test_the_jump_split_is_given_to_every_row():
    # g = Z_d splits off the hop alone, and every far element anticommutes with Z_d: E keeps pairing along it
    ee = pencil(3, E, ['ZII'], 'ZII', 2)
    assert (ee.near, ee.far) == (1, 1)


def test_sheets_count_beyond_the_jumps():
    # stage E by hand: U = a YYZ + b ZXX, its group is every string commuting with YYZ·ZXX ∝ XZY, dimension 5, and
    # the jumps XII, IZI, IIY with XXX and YYY, which commute with every term of H, span all five: no sheet
    stage_e = _add(_bond(0, 1, 3, 1, 'XY'), _bond(1, 2, 3, 2, 'XY'))     # COLOURING stage E
    ec = end_count(3, stage_e, ['XII', 'IZI', 'IIY'])
    assert sheets(3, stage_e, ['XII', 'IZI', 'IIY'], ec.far_basis[0]).extra == 0
    d = far_element_shape(3, D, ['ZII'])
    assert 'IIZ' in homogeneity_group(d.element, 3)                 # v always lit: the element is Z_v-homogeneous
    assert d.sheets.extra == 3                                       # this module's count
