"""Tests for Lindblad primitives (palindrome residual norm scaling formulas)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


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
# Correlated Z dephasing (experiments/THE_MIRROR_UNDER_A_SHARED_FIELD.md).
# Integer couplings and integer Kossakowski entries keep every entry a
# multiple of 1/2, so the comparisons below are exact.
# ----------------------------------------------------------------------


def _xxz_chain(N, delta):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for k in range(N - 1):
        for P, c in (('X', 1), ('Y', 1), ('Z', delta)):
            H = H + c * fw.site_op(N, k, P) @ fw.site_op(N, k + 1, P)
    return H


_A3 = np.array([[4., 1., -2.], [1., 6., 3.], [-2., 3., 5.]])     # every covariance nonzero


def test_correlated_z_dephasing_diagonal_a_is_local_dephasing():
    H = _xxz_chain(3, 2)
    g = [3., 5., 2.]
    assert np.array_equal(fw.lindbladian_correlated_z_dephasing(H, np.diag(g)),
                          fw.lindbladian_z_dephasing(H, g))


def test_correlated_z_dephasing_rate_law():
    # the dissipator is diagonal on matrix units |i><j|, entry -1/2 u^T a u with u = z(i) - z(j)
    N, d = 3, 8
    L = fw.lindbladian_correlated_z_dephasing(np.zeros((d, d), dtype=complex), _A3)
    z = lambda i: np.array([1 - 2 * ((i >> (N - 1 - k)) & 1) for k in range(N)], dtype=float)
    want = [-0.5 * (z(i) - z(j)) @ _A3 @ (z(i) - z(j)) for i in range(d) for j in range(d)]
    assert np.array_equal(L, np.diag(want).astype(complex))


def test_correlated_z_dephasing_price_is_an_operator():
    # palindrome_residual(L, tr a, N) carries the norm of -{K - tr a, .}: ||M||^2 = 8 4^N sum_{k<l} a_kl^2, on XXZ chains at Delta = 1, 2
    N = 3
    for delta in (1, 2):
        L = fw.lindbladian_correlated_z_dephasing(_xxz_chain(N, delta), _A3)
        M = fw.palindrome_residual(L, np.trace(_A3), N)
        assert np.sum(np.abs(M) ** 2) == 8 * 4 ** N * (1 ** 2 + 2 ** 2 + 3 ** 2)
    # local noise: the price is the number 2 sigma and the residual vanishes
    L = fw.lindbladian_correlated_z_dephasing(_xxz_chain(N, 2), np.diag(np.diag(_A3)))
    assert not fw.palindrome_residual(L, np.trace(_A3), N).any()
    # the operator identity itself, Pi_Z L Pi_Z^-1 = -L - {K, .}, with Pi_Z(|i><j|) = |j><i_bar| a permutation
    # of matrix units; integer inputs, so the comparison is exact
    d = 2 ** N
    P = np.zeros((d * d, d * d))
    for i in range(d):
        for j in range(d):
            P[j * d + (d - 1 - i), i * d + j] = 1
    Zs = [np.kron(np.kron(np.eye(2 ** k), np.diag([1., -1.])), np.eye(2 ** (N - 1 - k))) for k in range(N)]
    K = sum(_A3[k, l] * Zs[l] @ Zs[k] for k in range(N) for l in range(N))
    anti_K = np.kron(K, np.eye(d)) + np.kron(np.eye(d), K.T)
    L = fw.lindbladian_correlated_z_dephasing(_xxz_chain(N, 2), _A3)
    assert not (P @ L @ P.T + L + anti_K).any()
    assert (P @ L @ P.T + L + 2 * np.trace(_A3) * np.eye(d * d)).any()     # and the number 2 sigma does not do it


def test_shared_field_block_law():
    # a_kk = g_k, a_kl = kappa: on block (p, q) the generator is local at g_k - kappa, minus 2 kappa (p - q)^2
    N, d = 3, 8
    g, kappa = [7., 12., 9.], 3.
    H = _xxz_chain(N, 2)
    a = fw.shared_field_kossakowski(g, kappa)
    assert np.array_equal(a, np.array([[7., 3., 3.], [3., 12., 3.], [3., 3., 9.]]))
    L = fw.lindbladian_correlated_z_dephasing(H, a)
    L_loc = fw.lindbladian_z_dephasing(H, [x - kappa for x in g])
    pop = [bin(i).count('1') for i in range(d)]
    shift = [-2 * kappa * (pop[i] - pop[j]) ** 2 for i in range(d) for j in range(d)]
    assert np.array_equal(L - L_loc, np.diag(shift).astype(complex))
