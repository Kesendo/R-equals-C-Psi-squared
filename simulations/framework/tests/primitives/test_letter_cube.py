"""Tests for the letter cube primitives (docs/THE_ONE_SQUARE.md §7, §8).

The move rule is checked against an independent route, the matrix commutator of the placed
term with the string, not against the letter arithmetic the primitive uses.
"""
from __future__ import annotations

import sys
from itertools import product
from math import comb
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw

LET = {'I': np.eye(2, dtype=complex), 'X': np.array([[0, 1], [1, 0]], dtype=complex),
       'Y': np.array([[0, -1j], [1j, 0]]), 'Z': np.array([[1, 0], [0, -1]], dtype=complex)}


def mat(label):
    M = np.array([[1]], dtype=complex)
    for c in label:
        M = np.kron(M, LET[c])
    return M


def labels(N):
    return [''.join(t) for t in product('IXYZ', repeat=N)]


def test_letter_coordinates_are_the_even_tetrahedron():
    assert {c: fw.cube_coords(c) for c in 'IXYZ'} == {
        'I': (0, 0, 0), 'X': (1, 0, 1), 'Y': (1, 1, 0), 'Z': (0, 1, 1)}
    for c in 'IXYZ':
        for i, P in enumerate('ZXY'):
            anti = np.array_equal(LET[c] @ LET[P], -(LET[P] @ LET[c]))
            assert fw.cube_coords(c)[i] == int(anti)


def test_points_fill_the_tetrahedron_with_multinomial_counts():
    for N in range(1, 6):
        pts = fw.cube_points(N)
        assert len(pts) == comb(N + 3, 3)
        assert sum(pts.values()) == 4 ** N
        seen = {}
        for s in labels(N):
            k = fw.cube_coords(s)
            seen[k] = seen.get(k, 0) + 1
            assert sum(k) % 2 == 0 and sum(k) == 2 * fw.cube_weight(s)
        assert seen == pts


def test_rate_matches_the_dephasing_dissipator():
    # the dissipator sum_P gamma_P sum_l (P_l s P_l - s) on a string is -rate * s
    N, g = 3, {'Z': 0.25, 'X': 0.5, 'Y': 1.0}
    for s in labels(N):
        S = mat(s)
        D = sum(g[P] * (mat(''.join(P if j == l else 'I' for j in range(N))) @ S
                        @ mat(''.join(P if j == l else 'I' for j in range(N))) - S)
                for P in 'ZXY' for l in range(N))
        assert np.array_equal(D, -fw.cube_rate(s, g['Z'], g['X'], g['Y']) * S)


@pytest.mark.parametrize("N", [2, 3])
def test_step_agrees_with_the_matrix_commutator(N):
    terms = [t for t in labels(N) if t != 'I' * N]
    for t in terms:
        T = mat(t)
        for s in labels(N):
            S = mat(s)
            C = T @ S - S @ T
            r = fw.cube_step(t, s)
            if r is None:
                assert not C.any()
            else:
                target, dk = r
                assert any(np.array_equal(C, 2 * ph * mat(target)) for ph in (1, -1, 1j, -1j))
                assert dk == tuple(b - a for a, b in zip(fw.cube_coords(s), fw.cube_coords(target)))


def test_two_site_terms_change_the_weight_by_one_and_xyz_bonds_step_on_an_axis():
    for P, Q in product('XYZ', repeat=2):
        for s in labels(3):
            r = fw.cube_step(P + Q + 'I', s)
            if r is None:
                continue
            assert abs(fw.cube_weight(r[0]) - fw.cube_weight(s)) == 1
            if P == Q:
                nz = [i for i in range(3) if r[1][i]]
                assert len(nz) == 1 and abs(r[1][nz[0]]) == 2 and 'ZXY'[nz[0]] != P


def test_fields_keep_the_weight_and_step_on_a_face_diagonal():
    for P in 'XYZ':
        steps = fw.cube_moves([P + 'II'], 3)
        i = 'ZXY'.index(P)
        assert steps and all(d[i] == 0 and sorted(map(abs, d)) == [0, 1, 1] and sum(d) == 0 for d in steps)


def test_chain_terms_and_heisenberg_moves():
    H = fw.PauliHamiltonian.from_letter_tuples([('X', 'X'), ('Y', 'Y'), ('Z', 'Z')], chain_length=4)
    placed = fw.chain_terms(H)
    assert sorted(placed) == sorted(['XXII', 'IXXI', 'IIXX', 'YYII', 'IYYI', 'IIYY', 'ZZII', 'IZZI', 'IIZZ'])
    assert fw.cube_moves(placed, 4) == {(2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, -2, 0), (0, 0, 2), (0, 0, -2)}
    ising = fw.chain_terms(fw.PauliHamiltonian.from_letter_tuples([('Z', 'Z')], chain_length=4))
    assert all(d[0] == 0 for d in fw.cube_moves(ising, 4))        # ZZ never changes k_Z


def _lindbladian_in_strings(H, N, g):
    # L in the (unnormalized) string basis, divided by d: exact for integer H and dyadic g
    d = 2 ** N
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for l in range(N):
        A = mat(''.join('Z' if j == l else 'I' for j in range(N)))
        L = L + g * (np.kron(A, A.conj()) - np.kron(Id, Id))
    # the basis in the framework's index order (I, X, Z, Y per site), the order cube_coords_table uses
    B = np.array([mat(fw.pauli._pauli_label(k, N)).reshape(-1) for k in range(4 ** N)]).T
    return (B.conj().T @ L @ B) / d


def test_the_budget_is_exact():
    # the Hermitian part of L in the string basis is -2 g diag(k_Z), with or without a field, so
    # d/dt ||c||^2 = 2 Re c^dag L c = -4 g sum_s k_Z(s) |c_s|^2 at every instant, i.e.
    # d/dt ln ||rho||^2 = -4 g <k_Z>: the budget, with no quadrature in it
    N, g = 3, 0.25
    kz = fw.cube_coords_table(N)[:, 0]
    heis = sum(mat(t) for t in fw.chain_terms(
        fw.PauliHamiltonian.from_letter_tuples([('X', 'X'), ('Y', 'Y'), ('Z', 'Z')], chain_length=N)))
    field = heis + sum(mat(''.join('X' if j == l else 'I' for j in range(N))) for l in range(N))
    for H in (heis, field):
        Ls = _lindbladian_in_strings(H, N, g)
        assert np.array_equal((Ls + Ls.conj().T) / 2, np.diag(-2 * g * kz).astype(complex))


def test_coords_table_and_centroid_follow_the_basis_order():
    N = 3
    K = fw.cube_coords_table(N)
    assert K.shape == (4 ** N, 3)
    for k in range(4 ** N):
        assert tuple(K[k]) == fw.cube_coords(fw.pauli._pauli_label(k, N))
    # equal weight on XIZ, at (1, 1, 2), and IYI, at (1, 1, 0): the centroid is their midpoint
    c = np.zeros(4 ** N, dtype=complex)
    s1, s2 = 'XIZ', 'IYI'
    i1 = next(k for k in range(4 ** N) if fw.pauli._pauli_label(k, N) == s1)
    i2 = next(k for k in range(4 ** N) if fw.pauli._pauli_label(k, N) == s2)
    c[i1], c[i2] = 1.0, 1.0j
    assert fw.cube_centroid(c, N) == (1.0, 1.0, 1.0)
    rho = mat('ZII')
    assert fw.cube_centroid(fw.pauli_basis_vector(rho, N), N) == (0.0, 1.0, 1.0)


@pytest.mark.parametrize("N", [2, 3])
def test_commutator_coefficient_rebuilds_the_commutator(N):
    for t in labels(N):
        T = mat(t)
        for s in labels(N):
            C = T @ mat(s) - mat(s) @ T
            r = fw.commutator_coefficient(t, s)
            if r is None:
                assert not C.any()
            else:
                assert r[1] in (2j, -2j)
                assert np.array_equal(C, r[1] * mat(r[0]))


def test_per_site_rate_matches_the_dissipator_and_the_uniform_case():
    N = 3
    rates = [(0.25, 0.5, 0.0), (1.0, 0.0, 0.25), (0.0, 0.5, 0.5)]
    for s in labels(N):
        S = mat(s)
        D = sum(g * (mat(''.join(P if j == l else 'I' for j in range(N))) @ S
                     @ mat(''.join(P if j == l else 'I' for j in range(N))) - S)
                for l in range(N) for P, g in zip('ZXY', rates[l]))
        assert np.array_equal(D, -fw.cube_rate_per_site(s, rates) * S)
        assert fw.cube_rate_per_site(s, [(0.25, 0.5, 1.0)] * N) == fw.cube_rate(s, 0.25, 0.5, 1.0)


def test_inputs_are_checked():
    with pytest.raises(ValueError):
        fw.cube_coords('XA')
    with pytest.raises(ValueError):
        fw.cube_step('XX', 'XXI')
    with pytest.raises(ValueError):
        fw.cube_centroid(np.zeros(16), 2)
