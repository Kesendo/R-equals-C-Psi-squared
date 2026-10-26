"""Tests for framework.diagnostics.f158_end_count. The element is checked by an independent route (explicit 2^N
matrices: it commutes with H, anticommutes with every jump and is invertible), the counts against the C# twin's
known rows and against rows with an exact odd trace."""
from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from framework import cube_tables as ct  # noqa: E402
from framework.diagnostics.f158_end_count import end_count, end_count_of_setting, pauli_jump_letter  # noqa: E402


def op(s):
    return sp.kronecker_product(*[ct.PAULI[c] if c != 'I' else sp.eye(2) for c in s])


def matrix(terms):
    return sum((c * op(s) for s, c in terms.items()), sp.zeros(2 ** len(next(iter(terms)))))


def test_stage_e_weighted_element():
    # P3, a(XX+YY) on bond 0-1, b(XX+YY) on bond 1-2, jumps X, Z, Y: U = a YYZ + b ZXX (COLOURING stage E)
    a, b = 3, 4
    h = {'XXI': a, 'YYI': a, 'IXX': b, 'IYY': b}
    jumps = ['XII', 'IZI', 'IIY']
    e = end_count(3, h, jumps)
    assert e.pairs and e.near == e.far == 1
    k = e.element['YYZ']
    assert e.element == {'YYZ': k, 'ZXX': k * sp.Rational(b, a)}
    U, H = matrix(e.element), matrix(h)
    assert (H * U - U * H).is_zero_matrix
    assert all((op(j) * U + U * op(j)).is_zero_matrix for j in jumps)
    assert U * U == (k ** 2 * (1 + sp.Rational(b, a) ** 2)) * sp.eye(8)


def test_heisenberg_three_axes_breaks_and_two_axes_pairs():
    heis = {s: 1 for s in ('XXI', 'YYI', 'ZZI', 'IXX', 'IYY', 'IZZ')}
    assert not end_count(3, heis, ['XII', 'IZI', 'IIY']).pairs
    assert end_count(3, heis, ['XII', 'IZI', 'IIZ']).pairs


def test_longitudinal_field_on_a_dephased_site_breaks():
    # F158 (f5): Tr(H Z_0) != 0 rules the palindrome out
    assert not end_count(2, {'XX': 1, 'YY': 1, 'ZI': 1}, ['ZI', 'IZ']).pairs
    assert end_count(2, {'XX': 1, 'YY': 1, 'XI': 1}, ['ZI', 'IZ']).pairs


def test_setting_bridge():
    assert pauli_jump_letter(2 * ct.PAULI['Y']) == 'Y'
    assert pauli_jump_letter(ct.LOWER) is None
    s = ct.Setting(2, [('XX', (0, 1), 1)], [((0,), ct.LOWER, 1)])
    assert end_count_of_setting(s) is None
    s = ct.Setting(2, [('XX', (0, 1), 1), ('YY', (0, 1), 1)], [((l,), ct.PAULI['Z'], 1) for l in range(2)])
    assert end_count_of_setting(s).pairs
    # a zero rate removes its jump's condition: the count equals the count without that jump
    h = [('XX', (0, 1), 1), ('YY', (0, 1), 1), ('X', (1,), 1)]
    zero = end_count_of_setting(ct.Setting(2, h, [((0,), ct.PAULI['Z'], 1), ((1,), ct.PAULI['Z'], 0)]))
    absent = end_count_of_setting(ct.Setting(2, h, [((0,), ct.PAULI['Z'], 1)]))
    assert (zero.near, zero.far) == (absent.near, absent.far)


def _kernel_dim(M):
    from sympy.polys.matrices import DomainMatrix
    from sympy.polys.domains import QQ_I
    D = DomainMatrix.from_Matrix(M).convert_to(QQ_I)
    return M.rows - D.rank()


def test_a_row_with_unequal_nonzero_ends_against_the_generator():
    # Heisenberg P3, X dephasing on the middle site only, fields 2X on both ends: near 6, far 2, so the equality is
    # more than an emptiness test. Both counts checked against dim ker L and dim ker(L + 2 sigma) of the full
    # generator (row-major Lindblad formula, an independent route).
    h = {s: 1 for s in ('XXI', 'YYI', 'ZZI', 'IXX', 'IYY', 'IZZ')}
    h.update({'XII': 2, 'IIX': 2})
    e = end_count(3, h, ['IXI'])
    assert (e.near, e.far, e.pairs) == (6, 2, False)
    H, A, Id = matrix(h), op('IXI'), sp.eye(8)
    L = -sp.I * (sp.kronecker_product(H, Id) - sp.kronecker_product(Id, H.T))
    L += sp.kronecker_product(A, A.conjugate()) - sp.eye(64)
    assert _kernel_dim(L) == 6 and _kernel_dim(L + 2 * sp.eye(64)) == 2
