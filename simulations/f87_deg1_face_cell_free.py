#!/usr/bin/env python3
"""F87: the deg-1 face at m = 3 is CELL-FREE (Edge 1 of the PTF fresh-eyes chain, 2026-06-10).

THE CLAIM. For ANY Hermitian H (no diagonal-cell premise) under uniform Z-dephasing, the
m = 3 power-sum of the recentred generator M = A + gamma*Q is the exact monomial

    p_3(gamma) = P_{3,1} * gamma,    P_{3,1} = 6 * Sum_l Tr(H Z_l)^2 = 6 * 4^N * Sum_l c_l^2,

because the other three coefficients vanish identically for every Hermitian H:
  P_{3,0} = Tr(A^3)   = 0   (binomial antisymmetry: the j <-> 3-j terms of the A_L/A_R
                             expansion cancel pairwise since Tr((H^T)^j) = Tr(H^j));
  P_{3,2} = 3 Tr(Q^2 A) = 0 (Q^2 cross terms l != l' die on Tr(Z_l Z_l') = 0; the l = l'
                             terms die on Tr(H) - Tr(H^T) = 0);
  P_{3,3} = Tr(Q^3)   = 0   (any product of three single-site Z strings is a nonempty
                             Z string, traceless).
Consequence: ANY Hermitian H with a single-site-Z Pauli component (Sum c_l^2 > 0) breaks
the spectral palindrome at EVERY gamma > 0, with the closed-form leading coefficient
6*4^N*Sum c_l^2. This extends the girth ladder's ell = 1 face (PROOF_F87_WINDOWED_
MONOMIAL_CONVERSE.md Section 4, previously stated for diagonal-cell pairs) to all of
Hermitian operator space: the ell = 1 face is cell-free.

THE PTF TIE-IN. The Pi-break experiment (experiments/PTF_PALINDROME_BREAKING_
PERTURBATIONS.md, 2026-06-01) measured that a Z-field on the XY chain breaks the mirror
while the PTF closure survives. The mixed-cell H = XY chain + eps*Z_m sits OUTSIDE the
windowed diagonal-cell theorem's stated scope; this script closes that gap: the Z-field
row is now a theorem-grade all-gamma break with coefficient P_{3,1} = 6*4^N*eps^2.

Blocks (self-validating, every block raises on failure):
  1  XY chain + eps*Z_m at N = 4, 5: the three zero coefficients vanish (machine), and
     P_{3,1} == 6*4^N*eps^2 exactly, for eps in {0.3, 1.0} and two field sites.
  2  site-dependent field Sum_l eps_l Z_l: P_{3,1} == 6*4^N*Sum eps_l^2.
  3  RANDOM complex Hermitian H (seeded), arbitrary cell content: zeros vanish and
     P_{3,1} == 6*Sum_l Tr(H Z_l)^2 (the cell-free identity itself).
  4  control: pure XY chain (bipartite, no Z component): p_3 vanishes entirely.

Run: python simulations/f87_deg1_face_cell_free.py   (~10 s). 2026-06-10.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from framework.pauli import site_op  # noqa: E402


def build_super(H, N):
    d = 2 ** N
    A = -1j * (np.kron(H, np.eye(d)) - np.kron(np.eye(d), H.T))
    Q = sum(np.kron(site_op(N, l, 'Z'), site_op(N, l, 'Z').conj()) for l in range(N))
    return A, Q


def p3_coefficients(H, N):
    """(P30, P31, P32, P33) of p_3(gamma) = Tr((A + gamma*Q)^3), via cyclicity."""
    A, Q = build_super(H, N)
    A2 = A @ A
    P30 = complex(np.trace(A2 @ A))
    P31 = 3 * complex(np.trace(Q @ A2))
    P32 = 3 * complex(np.trace(Q @ Q @ A))
    P33 = complex(np.trace(Q @ Q @ Q))
    return P30, P31, P32, P33


def xy_chain(N, J=1.0):
    H = np.zeros((2 ** N, 2 ** N), dtype=complex)
    for l in range(N - 1):
        H += J * (site_op(N, l, 'X') @ site_op(N, l + 1, 'X')
                  + site_op(N, l, 'Y') @ site_op(N, l + 1, 'Y'))
    return H


def assert_monomial(H, N, want_P31, label, tol=1e-8):
    P30, P31, P32, P33 = p3_coefficients(H, N)
    scale = max(1.0, abs(want_P31))
    for name, val in [('P30', P30), ('P32', P32), ('P33', P33)]:
        assert abs(val) < tol * scale, f'{label}: {name} = {val} != 0'
    assert abs(P31.imag) < tol * scale and abs(P31.real - want_P31) < tol * scale, \
        f'{label}: P31 = {P31} != {want_P31}'
    print(f'  {label}: P30=P32=P33=0 (machine), P31 = {P31.real:.6f} == {want_P31:.6f}  OK',
          flush=True)


def block1():
    print('BLOCK 1 -- XY chain + eps*Z_m (the PTF Pi-break configuration)', flush=True)
    for N in (4, 5):
        for eps in (0.3, 1.0):
            for m in (0, N // 2):
                H = xy_chain(N) + eps * site_op(N, m, 'Z')
                assert_monomial(H, N, 6 * 4 ** N * eps ** 2,
                                f'N={N} eps={eps} site={m}')
    print('BLOCK 1 PASS', flush=True)


def block2():
    print('BLOCK 2 -- site-dependent Z field', flush=True)
    N = 5
    eps = [0.7, -0.2, 0.0, 1.3, 0.4]
    H = xy_chain(N) + sum(e * site_op(N, l, 'Z') for l, e in enumerate(eps))
    assert_monomial(H, N, 6 * 4 ** N * sum(e ** 2 for e in eps),
                    f'N={N} eps={eps}')
    print('BLOCK 2 PASS', flush=True)


def block3():
    print('BLOCK 3 -- random complex Hermitian H (the cell-free identity)', flush=True)
    rng = np.random.default_rng(87)
    for N in (3, 4):
        d = 2 ** N
        R = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
        H = R + R.conj().T
        want = 6 * sum(complex(np.trace(H @ site_op(N, l, 'Z'))).real ** 2
                       for l in range(N))
        assert_monomial(H, N, want, f'N={N} random Hermitian')
    print('BLOCK 3 PASS', flush=True)


def block4():
    print('BLOCK 4 -- control: pure XY chain (bipartite, no Z component)', flush=True)
    for N in (4, 5):
        assert_monomial(xy_chain(N), N, 0.0, f'N={N} pure XY')
    print('BLOCK 4 PASS', flush=True)


if __name__ == '__main__':
    block1()
    block2()
    block3()
    block4()
    print('\nALL BLOCKS PASS -- the ell = 1 face is cell-free: any Hermitian H with a')
    print('single-site-Z component breaks the palindrome at every gamma > 0, with')
    print('p_3(gamma) = 6*4^N*Sum_l c_l^2 * gamma exactly.')
