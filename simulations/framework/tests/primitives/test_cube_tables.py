"""Tests for framework.cube_tables (Device B's tables). Each number is checked against a route the module does
not use: the hand-computed Pauli transfer of sigma^-, the C# TransferTimesFour numbers, the trace argument for
the support sectors, and a direct commutator for one table entry."""
from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from framework import cube_tables as ct  # noqa: E402

X, Y, Z = ct.PAULI['X'], ct.PAULI['Y'], ct.PAULI['Z']


def test_transfer_of_sigma_minus_is_the_csharp_table():
    # 2 Tr(P_a D(P_b)) for D = D[|0><1|] at rate 1, rows a and columns b in I, X, Y, Z: the lights X, Y at 2,
    # Z at 4, and the arrow I -> Z with coefficient 4 (C# LetterCube.TransferTimesFour, "lights plus the arrow")
    expected = sp.Matrix([[0, 0, 0, 0], [0, -2, 0, 0], [0, 0, -2, 0], [4, 0, 0, -4]])
    assert ct.transfer_times_four(ct.LOWER, 1) == expected
    assert ct.transfer_times_four(ct.RAISE, 1) == sp.Matrix([[0, 0, 0, 0], [0, -2, 0, 0], [0, 0, -2, 0],
                                                             [-4, 0, 0, -4]])


def test_column_convention_is_pinned_by_the_arrow():
    # the arrow I -> Z sits in row Z, column I; the transposed convention would put it in row I, column Z
    T = ct.to_pauli(ct.dissipator(ct.LOWER, 1))
    assert T[3, 0] == 1 and T[0, 3] == 0


def test_z_dephasing_is_a_height_in_the_pauli_basis():
    g = sp.Rational(3, 7)
    assert ct.to_pauli(ct.dissipator(Z, g)) == sp.diag(0, -2 * g, -2 * g, 0)


def test_basis_change_round_trip():
    M = ct.dissipator(ct.LOWER, sp.Rational(1, 2)) + ct.dissipator(ct.RAISE, sp.Rational(1, 8))
    assert ct.to_cells(ct.to_pauli(M)) == M.applyfunc(sp.expand)


def test_the_closure_has_192_monomial_moves_and_8_cell_diagonal_ones():
    moves = ct.cube_moves()
    assert len(moves) == 192
    for _, M in moves:
        T = ct.to_pauli(M)
        assert all(sum(1 for j in range(4) if T[i, j] != 0) == 1 for i in range(4))
        assert all(sum(1 for i in range(4) if T[i, j] != 0) == 1 for j in range(4))
    assert sum(1 for _, M in moves if ct.is_cell_diagonal(M)) == 8
    right, left = ct.move_cosets('right'), ct.move_cosets('left')
    assert len(right) == 24 and len(left) == 24
    assert {k for k, _ in right} != {k for k, _ in left}      # the 8 cell-diagonal moves are not normal


def test_a_support_table_lives_in_the_traceless_sector():
    # every site factor of -i[P_S, .] is traceless, so the table has no component along I (x) B or A (x) I
    s = ct.Setting(2, [('XX', (0, 1), 1), ('YZ', (0, 1), sp.Rational(1, 3)), ('ZZ', (0, 1), 2)])
    L = ct.support_tables(s)[(0, 1)]
    Ia, Ib = sp.zeros(4, 4), sp.zeros(4, 4)
    for i in range(4):
        for j in range(4):
            for k in range(4):
                Ia[i, j] += L[k * 4 + i, k * 4 + j]
                Ib[i, j] += L[i * 4 + k, j * 4 + k]
    assert Ia.is_zero_matrix and Ib.is_zero_matrix


def test_one_support_table_entry_against_the_commutator():
    # rho = |01><10| (site 0 ket 0, site 1 ket 1; bra reversed); -i[XX + YY, rho] read off directly
    s = ct.Setting(2, [('XX', (0, 1), 1), ('YY', (0, 1), 1)])
    L = ct.support_tables(s)[(0, 1)]
    H = ct.hamiltonian_on(s, (0, 1))
    rho = sp.zeros(4, 4)
    rho[1, 2] = 1
    out = (-sp.I * (H * rho - rho * H)).applyfunc(sp.expand)
    # per-site order index of |i0 i1><j0 j1| is (2 i0 + j0) * 4 + (2 i1 + j1)
    def idx(i0, i1, j0, j1):
        return (2 * i0 + j0) * 4 + (2 * i1 + j1)
    col = idx(0, 1, 1, 0)
    for i0, i1, j0, j1 in [(a, b, c, d) for a in (0, 1) for b in (0, 1) for c in (0, 1) for d in (0, 1)]:
        assert L[idx(i0, i1, j0, j1), col] == out[2 * i0 + i1, 2 * j0 + j1]


def test_outside_reasons():
    assert ct.outside_reason(ct.Setting(4, [('XXXX', (0, 1, 2, 3), 1)])) is not None
    assert ct.outside_reason(ct.Setting(2, [('XY', (0, 1), sp.I)])) is not None
    assert ct.outside_reason(ct.Setting(2, [], [((0, 1), sp.eye(4), 1)])) is not None
    assert ct.outside_reason(ct.Setting(2, [('XX', (0, 1), 1)], [((0,), ct.LOWER, sp.Rational(1, 2))])) is None
