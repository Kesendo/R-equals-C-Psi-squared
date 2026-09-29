"""F158 §(f5) odd-word pre-filter: exact certificates, soundness against the two-end count."""
import itertools
from fractions import Fraction

import numpy as np
import pytest
import sympy as sp

import framework as fw

LET = {'I': sp.eye(2), 'X': sp.Matrix([[0, 1], [1, 0]]), 'Y': sp.Matrix([[0, -sp.I], [sp.I, 0]]),
       'Z': sp.diag(1, -1)}


def _dense(s):
    M = sp.Matrix([[1]])
    for ch in s:
        M = sp.kronecker_product(M, LET[ch])
    return M


def test_pauli_product_matches_matrices():
    for a, b in itertools.product('IXYZ', repeat=2):
        e, u = fw.pauli_product(a, b)
        assert LET[a] * LET[b] == sp.I ** e * LET[u]


def test_detuning_on_a_dephased_site_is_one_trace():
    # XY chain, Z dephasing everywhere, a field 3/10 on the last site: Tr(H Z_2) = 2^3 * 3/10
    out = fw.odd_word_obstruction({'XXI': 1, 'YYI': 1, 'IXX': 1, 'IYY': 1, 'IIZ': Fraction(3, 10)},
                                  ['ZII', 'IZI', 'IIZ'])
    assert out['word'] == 'H^1·A2' and out['trace'] == (Fraction(12, 5), 0)


def test_cascade_obstruction_is_the_squared_word():
    # the colouring page's cascade bond with a field eps X on the first site: Tr(H^2 IX) = 8 eps
    eps = Fraction(1, 3)
    h = {'XX': 1, 'YY': 1, 'YI': Fraction(3, 10), 'IZ': Fraction(11, 50), 'XI': eps}
    out = fw.odd_word_obstruction(h, ['IX'])
    assert out['word'] == 'H^2·A0' and out['trace'] == (8 * eps, 0)
    assert fw.odd_word_obstruction({k: v for k, v in h.items() if k != 'XI'}, ['IX']) is None


def test_silent_on_the_palindromic_heisenberg_chain():
    h = {s: 1 for s in ('XXI', 'YYI', 'ZZI', 'IXX', 'IYY', 'IZZ')}
    assert fw.odd_word_obstruction(h, ['ZII', 'IZI', 'IIZ']) is None


def _palindromic(h, jumps):
    """F158 exactly: dim ker ad_H on the strings commuting with every jump equals it on the lit ones"""
    n = len(jumps[0])
    H = sum((sp.nsimplify(c) * _dense(s) for s, c in h.items()), sp.zeros(2 ** n))

    def anti(s, t):
        return sum(x != 'I' and y != 'I' and x != y for x, y in zip(s, t)) % 2 == 1
    strings = [''.join(p) for p in itertools.product('IXYZ', repeat=n)]
    dims = []
    for want_anti in (False, True):
        cols = [_dense(s) for s in strings if all(anti(s, a) == want_anti for a in jumps)]
        M = sp.Matrix.hstack(*[(H * c - c * H).reshape(4 ** n, 1) for c in cols])
        dims.append(len(cols) - M.rank())
    return dims[0] == dims[1]


def test_sound_on_palindromic_rows_and_catches_broken_ones():
    # N = 3 path, ZZ or XX + YY bonds, one dephasing letter on each of two sites, fields (0.3, 0.22, 0.41):
    # the pre-filter must never fire where F158 says palindromic; on broken rows it may or may not
    rng = np.random.default_rng(3)
    mags = [Fraction(3, 10), Fraction(11, 50), Fraction(41, 100)]
    fired_on_pal, caught, broken = 0, 0, 0
    for _ in range(40):
        bset = [('Z',), ('X', 'Y')][rng.integers(2)]
        deph = [rng.choice(['X', 'Y', 'Z']) for _ in range(3)]
        fields = [rng.choice(['I', 'X', 'Y', 'Z']) for _ in range(3)]
        h = {}
        for a, b in ((0, 1), (1, 2)):
            for P in bset:
                s = ['I'] * 3
                s[a] = s[b] = P
                h[''.join(s)] = 1
        for l, P in enumerate(fields):
            if P != 'I':
                s = ['I'] * 3
                s[l] = P
                h[''.join(s)] = mags[l]
        jumps = [''.join(deph[k] if k == l else 'I' for k in range(3)) for l in range(3)]
        pal = _palindromic(h, jumps)
        out = fw.odd_word_obstruction(h, jumps, max_power=3)
        fired_on_pal += pal and out is not None
        broken += not pal
        caught += (not pal) and out is not None
    assert fired_on_pal == 0
    assert broken > 0 and caught > 0          # the filter is not silent everywhere


def test_needs_a_jump():
    with pytest.raises(ValueError):
        fw.odd_word_obstruction({'XX': 1}, [])
