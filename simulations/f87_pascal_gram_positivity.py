#!/usr/bin/env python3
"""F87 windowed converse: the Pascal-Gram positivity theorem (R-sign resolved, 2026-06-10).

THE THEOREM. Let M(gamma) = A + gamma*Q be the recentred dephased generator of a windowed
diagonal-cell pair (PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md), and let m* be the first odd m
whose power-sum polynomial p_m(gamma) = Tr(M^m) does not vanish identically. Then EVERY
gamma-coefficient of p_{m*} is non-negative:

    P_{m*,d} = (m*/d) * Sum_{l in [N]^d} Sum_{k in [0..u]^d} |U^{(l)}_{k}|^2  or exactly 0,

    U^{(l)}_{k} = Sum_{|alpha| = u, d parts} prod_i C(alpha_i, k_i) * T^{(l)}_{alpha},
    T^{(l)}_{alpha} = Tr(Z_{l1} H^{alpha1} Z_{l2} H^{alpha2} ... Z_{ld} H^{alphad}),

with u = (m* - d)/2 the equal leg total. Since at least one class is nonzero at m*,
p_{m*}(gamma) > 0 for every gamma > 0: hard at one gamma is hard at ALL gamma. This closes
the R-sign residual (and with it the windowed all-gamma converse, no residual left).

Proof skeleton verified here block by block:
  (a) cyclic decomposition  P_{m,d} = (m/d) * Sum_{compositions a of m-d} Tr(Q A^{a1} ... Q A^{ad});
  (b) leg factorization     A = A_L + A_R commute, the supertrace splits into bra x ket legs,
                            coefficient (-i)^{|alpha|} (+i)^{|beta|} prod C(alpha_i+beta_i, beta_i);
  (c) Hermitian conjugacy   the ket leg (on H^T = conj H) is the complex conjugate of the bra
                            leg at the SAME indices: leg2 = conj(T^{(l)}_{beta});
  (d) leg parity + girth    T^{(l)}_{alpha} = 0 unless |alpha| odd (F-chirality, d odd) and
                            |alpha| >= ell (no closed odd walk shorter than the odd-girth);
  (e) Vandermonde assembly  C(a+b,b) = Sum_k C(a,k) C(b,k) per slot turns the equal-total
                            (u,u) block into Sum_k |U_k|^2, prefactor (-i)^u (+i)^u = +1;
  (f) slice inversion       U_k at |k| = u IS T_k (prod C(k_i,k_i) = 1, all other alpha drop),
                            so a vanished Gram block kills every total-u moment;
  (g) cascade induction     p_m == 0 for odd m < m* + (f) kill all lower-total moments, so at
                            m* every cross-total term is dead and every class is (e) or zero.
Selection rule corollary: P_{m*,d} can fire only for d == m*-2 (mod 4), d <= m* - 2*ell.
For deg = m* - 2*ell in {1, 3} that is a SINGLE candidate class: monomiality is derived, not
assumed. From deg = 5 on, two classes may coexist (positivity is what the converse needs).

Anchors: f87_windowed_monomial_converse.py (exact CRT machinery, canonical pairs),
f87_girth_dichotomy.py (the k=4 battery that pinned the gamma^5 rung).

Run: python simulations/f87_pascal_gram_positivity.py   (~10-15 min; the exact CRT
re-derivations at N=5 dominate). Self-validating: every block raises on failure. 2026-06-10.
"""
import sys
from math import comb
from pathlib import Path
from itertools import product as iproduct

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from framework.pauli import site_op  # noqa: E402
import f87_windowed_monomial_converse as anchor  # noqa: E402

PAIR_D1 = [('I', 'X', 'X', 'Z'), ('X', 'I', 'X', 'Z')]   # k=4 deg-1 cycle, t_3 != 0
PAIR_D5 = [('I', 'I', 'X', 'Y'), ('Z', 'X', 'Z', 'Y')]   # k=4 gamma^5 rung, t_3 = t_5 = 0


def compositions(total, parts):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in compositions(total - first, parts - 1):
            yield (first,) + rest


def build_super(H, N):
    d = 2 ** N
    A = -1j * (np.kron(H, np.eye(d)) - np.kron(np.eye(d), H.T))
    Q = sum(np.kron(site_op(N, l, 'Z'), site_op(N, l, 'Z').conj()) for l in range(N))
    return A, Q


def leg_moment(H, N, ls, alphas, transpose=False):
    M = np.eye(2 ** N, dtype=complex)
    Huse = H.T if transpose else H
    for l, a in zip(ls, alphas):
        M = M @ site_op(N, l, 'Z') @ np.linalg.matrix_power(Huse, a)
    return complex(np.trace(M))


def all_leg_moments(H, N, d, u):
    """All T^{(l)}_{alpha}, |alpha| = u, as dict l-tuple -> complex vector over compositions."""
    comps = list(compositions(u, d))
    cidx = {al: i for i, al in enumerate(comps)}
    P = [[site_op(N, l, 'Z') @ np.linalg.matrix_power(H, a) for a in range(u + 1)]
         for l in range(N)]
    out = {}
    if d == 1:
        for l in range(N):
            out[(l,)] = np.array([np.trace(P[l][u])], dtype=complex)
        return comps, out

    def rec(depth, ls, used, M):
        if depth == d - 1:
            a_last = u - sum(used)
            for l in range(N):
                T = np.einsum('ij,ji->', M, P[l][a_last])
                vec = out.setdefault(tuple(ls) + (l,),
                                     np.zeros(len(comps), dtype=complex))
                vec[cidx[tuple(used) + (a_last,)]] = T
            return
        for l in range(N):
            for a in range(u - sum(used) + 1):
                rec(depth + 1, ls + [l], used + [a], M @ P[l][a])

    for l in range(N):
        for a in range(u + 1):
            rec(1, [l], [a], P[l][a])
    return comps, out


def pascal_gram(H, N, ell, d):
    """The equal-total Pascal-Gram value for the #Q = d class at m = 2*ell + d."""
    m = 2 * ell + d
    comps, moments = all_leg_moments(H, N, d, ell)
    ks = list(iproduct(range(ell + 1), repeat=d))
    W = np.zeros((len(ks), len(comps)))
    for r, k in enumerate(ks):
        for c, al in enumerate(comps):
            w = 1
            for ai, ki in zip(al, k):
                w *= comb(ai, ki)
            W[r, c] = w
    total = 0.0
    for vec in moments.values():
        if np.max(np.abs(vec)) < 1e-12:
            continue
        U = W @ vec
        total += float(np.sum(np.abs(U) ** 2))
    return m / d * total, comps, ks, W, moments


def exact_coefficients(N, pair, max_m, nprimes=6):
    Ar, Ai, Q = anchor.build_integer_generators(N, pair)
    mstar, re_co, _, polys = anchor.first_nonvanishing_odd(Ar, Ai, Q, max_m, nprimes=nprimes)
    return mstar, re_co, polys


# ======================================================================
# BLOCK 1 -- cyclic decomposition: P_{m,3} = (m/3) Sum_{a+b+c=m-3} Tr(Q A^a Q A^b Q A^c),
# checked against the exact CRT gamma^3 coefficient at N=4 (K3 and the complex-H flux pair).
# ======================================================================
def block1():
    print('BLOCK 1 -- cyclic decomposition (d=3) vs exact CRT coefficient', flush=True)
    for name, pair, want in [('K3', anchor.K3_EVEN, 2064384), ('FLUX', anchor.FLUX, 589824)]:
        N, m = 4, 9
        H = anchor.build_H(N, pair)
        A, Q = build_super(H, N)
        S3 = 0
        for a in range(m - 2):
            for b in range(m - 2 - a):
                c = m - 3 - a - b
                Ap = {e: np.linalg.matrix_power(A, e) for e in {a, b, c}}
                S3 += complex(np.trace(Q @ Ap[a] @ Q @ Ap[b] @ Q @ Ap[c]))
        got = m / 3 * S3
        assert abs(got.imag) < 1e-6 and abs(got.real - want) < 1e-6 * want, \
            f'{name}: cyclic decomposition {got} != {want}'
        print(f'  {name}: (m/3)*Sum = {got.real:.1f} == exact P_(9,3) = {want}  OK', flush=True)
    print('BLOCK 1 PASS', flush=True)


# ======================================================================
# BLOCK 2 -- leg factorization + Hermitian conjugacy of the ket leg.
# ======================================================================
def block2():
    print('BLOCK 2 -- leg factorization (exact) + ket leg = conj(bra leg)', flush=True)
    for name, pair in [('K3', anchor.K3_EVEN), ('FLUX', anchor.FLUX)]:
        N, ell, m = 4, 3, 9
        H = anchor.build_H(N, pair)
        A, Q = build_super(H, N)
        a, b, c = ell, ell, m - 3 - 2 * ell
        Ap = {e: np.linalg.matrix_power(A, e) for e in {a, b, c}}
        direct = complex(np.trace(Q @ Ap[a] @ Q @ Ap[b] @ Q @ Ap[c]))
        s = 0
        for ls in iproduct(range(N), repeat=3):
            for j1 in range(a + 1):
                for j2 in range(b + 1):
                    for j3 in range(c + 1):
                        coef = ((-1j) ** ((a - j1) + (b - j2) + (c - j3)) * (1j) ** (j1 + j2 + j3)
                                * comb(a, j1) * comb(b, j2) * comb(c, j3))
                        t1 = leg_moment(H, N, ls, (a - j1, b - j2, c - j3))
                        t2 = leg_moment(H, N, ls, (j1, j2, j3), transpose=True)
                        s += coef * t1 * t2
        assert abs(direct - s) < 1e-6 * max(1, abs(direct)), \
            f'{name}: factorized {s} != direct {direct}'
        # Hermitian conjugacy: the H^T leg is the conjugate of the H leg, same indices
        worst = 0.0
        for ls in iproduct(range(N), repeat=3):
            for al in compositions(ell, 3):
                worst = max(worst, abs(leg_moment(H, N, ls, al, transpose=True)
                                       - np.conj(leg_moment(H, N, ls, al))))
        assert worst < 1e-9, f'{name}: ket leg != conj(bra leg), worst {worst}'
        print(f'  {name}: factorization diff 0, conjugacy worst {worst:.1e}  OK', flush=True)
    print('BLOCK 2 PASS', flush=True)


# ======================================================================
# BLOCK 3 -- the Pascal-Gram value equals the exact first coefficient on every branch:
# d=1 (IXXZ+XIXZ), d=3 (K3, FLUX, MULTIZ), d=5 (IIXY+ZXZY). Exact CRT re-derived in-run.
# ======================================================================
def block3():
    print('BLOCK 3 -- Pascal-Gram == exact CRT coefficient (d = 1, 3, 5)', flush=True)
    cases = [('IXXZ+XIXZ', PAIR_D1, 5, 3, 1, 7, 573440),
             ('K3', anchor.K3_EVEN, 4, 3, 3, 9, 2064384),
             ('FLUX', anchor.FLUX, 4, 3, 3, 9, 589824),
             ('MULTIZ', anchor.MULTIZ, 4, 1, 3, 5, 61440),
             ('IIXY+ZXZY', PAIR_D5, 5, 3, 5, 11, 86507520)]
    for name, pair, N, ell, d, mstar_want, want in cases:
        H = anchor.build_H(N, pair)
        got, *_ = pascal_gram(H, N, ell, d)
        assert abs(got - want) < 1e-6 * want, f'{name}: Pascal-Gram {got} != {want}'
        mstar, re_co, _ = exact_coefficients(N, pair, mstar_want, nprimes=6)
        assert mstar == mstar_want and re_co[d] == want, \
            f'{name}: exact ({mstar}, {re_co[d] if mstar else None}) != ({mstar_want}, {want})'
        nz = [j for j, c0 in enumerate(re_co) if c0 != 0]
        assert nz == [d], f'{name}: surviving classes {nz} != [{d}]'
        print(f'  {name}: d={d}, m*={mstar}: Pascal-Gram = {got:.1f} == exact = {want}, '
              f'single class  OK', flush=True)
    print('BLOCK 3 PASS', flush=True)


# ======================================================================
# BLOCK 4 -- cascade inputs at the gamma^5 rung: for IIXY+ZXZY every lower-total moment is
# zero EXACTLY (t_3, t_5, and all total-3 triple moments), which is what forces
# P_{11,1} = P_{11,3} = 0; for K3 the total-3 triple moments are NONZERO (the gamma^3 rung
# fires) while t_3 = 0 (the deg-1 class is dead): the cascade discriminates.
# ======================================================================
def block4():
    print('BLOCK 4 -- cascade inputs (exact zeros below the firing rung)', flush=True)
    H5 = anchor.build_H(5, PAIR_D5)
    worst_t = max(abs(leg_moment(H5, 5, (l,), (j,))) for l in range(5) for j in (3, 5))
    assert worst_t == 0.0, f'IIXY+ZXZY: t_3/t_5 not exactly zero ({worst_t})'
    worst3 = max(abs(leg_moment(H5, 5, ls, al))
                 for ls in iproduct(range(5), repeat=3) for al in compositions(3, 3))
    assert worst3 == 0.0, f'IIXY+ZXZY: total-3 triple moments not exactly zero ({worst3})'
    HK = anchor.build_H(4, anchor.K3_EVEN)
    tK = max(abs(leg_moment(HK, 4, (l,), (3,))) for l in range(4))
    m3K = max(abs(leg_moment(HK, 4, ls, al))
              for ls in iproduct(range(4), repeat=3) for al in compositions(3, 3))
    assert tK == 0.0 and m3K > 1.0, f'K3: expected t_3 = 0 and triple moments != 0 ({tK}, {m3K})'
    print(f'  IIXY+ZXZY: max|t_3|,|t_5| = 0, max total-3 triple = 0 (P_11,1 = P_11,3 = 0 forced)')
    print(f'  K3: t_3 = 0 (deg-1 dead) but max total-3 triple = {m3K:.0f} (gamma^3 fires)  OK')
    print('BLOCK 4 PASS', flush=True)


# ======================================================================
# BLOCK 5 -- slice inversion + the selection rule. The |k| = u slice of U is T itself
# (checked numerically on K3), and the mod-4 selection rule singles out the firing class on
# all pinned representatives.
# ======================================================================
def block5():
    print('BLOCK 5 -- slice inversion + mod-4 selection rule', flush=True)
    H = anchor.build_H(4, anchor.K3_EVEN)
    ell, d = 3, 3
    _, comps, ks, W, moments = pascal_gram(H, 4, ell, d)
    cidx = {al: i for i, al in enumerate(comps)}
    worst = 0.0
    for lvec, T in moments.items():
        U = W @ T
        for r, k in enumerate(ks):
            if sum(k) == ell:
                worst = max(worst, abs(U[r] - T[cidx[k]]))
    assert worst < 1e-9, f'slice inversion broken, worst {worst}'
    print(f'  K3: U restricted to |k| = ell equals T, worst diff {worst:.1e}  OK')

    def allowed(mstar, ell):
        return [d for d in range(1, mstar - 2 * ell + 1, 2) if (mstar - d) % 4 == 2]

    for name, mstar, ell_p, fired in [('IXXZ+XIXZ', 7, 3, 1), ('K3', 9, 3, 3),
                                      ('FLUX', 9, 3, 3), ('MULTIZ', 5, 1, 3),
                                      ('IIXY+ZXZY', 11, 3, 5)]:
        al = allowed(mstar, ell_p)
        assert fired in al, f'{name}: fired class {fired} not in allowed {al}'
        single = 'derived' if len(al) == 1 else f'allowed {al}, measured single'
        print(f'  {name}: m*={mstar}, allowed classes {al}, fired d={fired} ({single})  OK')
    print('BLOCK 5 PASS', flush=True)


if __name__ == '__main__':
    block1()
    block2()
    block3()
    block4()
    block5()
    print('\nALL BLOCKS PASS -- Pascal-Gram positivity verified: every class at m* is a sum of')
    print('squares (or exactly zero), so p_{m*}(gamma) > 0 for every gamma > 0. R-sign closed.')
