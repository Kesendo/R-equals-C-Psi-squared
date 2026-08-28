"""The pairing condition: one criterion for the dephasing palindrome, gated.

THE STATEMENT UNDER TEST. For L(rho) = -i[H, rho] + sum_l gamma_l (A_l rho A_l
- rho), every A_l a Pauli letter on one site, and sigma = sum_l gamma_l:

    det(x*I - L) == det((-x - 2*sigma)*I - L)     the char-poly palindrome
      <==>
    there is an INVERTIBLE U with  [U, H] = 0  and  U A_l U^-1 = -A_l  for all l.

The left side is the eigenvalue MULTISET statement, not the set statement, and
not F1's operator identity; it is what the committed GF(p) kernel decides.

WHAT U IS, and it is not only an algebraic device: the two conditions give
L(U) = -2*sigma*U outright, so U is an exact eigenmode of the Liouvillian at
-2*sigma, the mirror partner of the steady state at 0 and the far end of the
axis the palindrome folds about.

WHAT U BUYS, IN TWO STEPS, AND THE SECOND ONE IS CONDITIONAL.

    L_U L L_U^-1 = -L^dagger - 2*sigma      one-sided, needs only U (gate 3)
    Pi = L_U o T, Pi(rho) = U rho^T
    Pi L Pi^-1   = -L - 2*sigma             needs H = H^T as well (gates 3, 4)

A Y field makes H non-real and the second step FAILS there, which gate 4 pins
on a row inside the scored grid. The palindrome survives anyway, and the reason
is not modular and cannot be gated here: L preserves hermiticity, so its
characteristic polynomial is real and its spectrum is conjugation-closed, which
collapses the metric relation's reflection lambda -> -conj(lambda) - 2*sigma
onto the palindrome. That step is stated, not measured, and the page says so.

WHY THE CRITERION IS DECIDABLE AND NOT SCANNED. Both conditions are LINEAR in U,
so the admissible U form a subspace whose dimension is a rank over GF(p), and
the predicate is "that subspace contains an invertible element". No catalogue of
candidate operators, no group theory, no eigensolver.

VERDICT DISCIPLINE, and the two sides are not alike.
  The palindrome verdict is the exact GF(p) kernel of f138_exact_palindrome_test:
  a BREAK is proved, a HOLD is certified over three primes.
  The criterion verdict is EXACT when the subspace is empty, because a nullity
  read mod p can only come out too large. A NONEMPTY subspace can be an artifact
  of the modulus, so it is confirmed at all three primes before a row counts as
  admissible. Invertibility inside a nonempty subspace is then decided by
  sampling; with a degree-2^N determinant and p near 2^30 three independent
  singular draws on a subspace that holds an invertible element run below 1e-24
  at N=3 and below 5e-24 at N=4. Rows that come back nonempty-but-singular are
  counted and printed rather than folded in, and the residual doubt there points
  at a hidden FALSE POSITIVE, which is the load-bearing direction.
  Scoring runs in BOTH directions throughout.

Cost: about 30 minutes, dominated by gates 5 and 6.

Run:  python simulations/f138_pairing_condition.py
Out:  simulations/results/f138_pairing_condition.txt
"""
import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from f138_exact_palindrome_test import (PRIMES, det_mod_np_p,  # noqa: E402
                                        inv_p, sqrt_minus_one)
from f138_clause_two_sweep import (build_L, palindromic,  # noqa: E402
                                   predicate, sign_patterns, LETTER_NAME)

PRIME = PRIMES[0]
rng = np.random.default_rng(20260826)

P3 = [(0, 1), (1, 2)]
K3 = [(0, 1), (1, 2), (2, 0)]
BI = [(0, 1)]
P4 = [(0, 1), (1, 2), (2, 3)]
C4 = [(0, 1), (1, 2), (2, 3), (3, 0)]

# The row the open arc f138_converse_failures names, assembled once.
ROW = dict(n=3, edges=P3, deph=(0, 1, 0), fld=(1, 0, 1), signs=(1, 1, -1),
           mags=(30, 30, 30))
GEN3 = (30, 22, 41)

RESULTS = []


def gate(name, ok, detail=''):
    RESULTS.append((name, bool(ok)))
    print('  [%s] %-52s %s' % ('PASS' if ok else 'FAIL', name, detail))
    return ok


# --------------------------------------------------------------------------
# operators


def paulis(p):
    i = sqrt_minus_one(p)
    return [np.eye(2, dtype=np.int64),
            np.array([[0, 1], [1, 0]], dtype=np.int64),
            np.array([[0, (p - i) % p], [i, 0]], dtype=np.int64),
            np.array([[1, 0], [0, p - 1]], dtype=np.int64)]


def string_op(letters, p):
    PA = paulis(p)
    out = np.ones((1, 1), dtype=np.int64)
    for l in letters:
        out = np.kron(out, PA[l]) % p
    return out


def site_perm(sigma, n, p):
    d = 2 ** n
    P = np.zeros((d, d), dtype=np.int64)
    for b in range(d):
        bits = [(b >> (n - 1 - s)) & 1 for s in range(n)]
        nb = 0
        for s in range(n):
            nb |= bits[s] << (n - 1 - sigma[s])
        P[nb, b] = 1
    return P


def build_H(n, edges, fld, signs, mags, p, j_num=None, bond_terms=(1, 2, 3)):
    """The Hamiltonian build_L uses, re-derived independently here."""
    PA = paulis(p)
    d = 2 ** n
    H = np.zeros((d, d), dtype=np.int64)
    eye = np.eye(2, dtype=np.int64)
    for b, (a, c) in enumerate(edges):
        jj = 1 if j_num is None else (j_num[b] * inv_p(100, p)) % p
        for t in bond_terms:
            op = np.ones((1, 1), dtype=np.int64)
            for k in range(n):
                op = np.kron(op, PA[t] if k in (a, c) else eye) % p
            H = (H + jj * op) % p
    sg = signs if signs is not None else (1,) * n
    for s in range(n):
        if fld[s] == 0:
            continue
        h = (mags[s] * sg[s] % p) * inv_p(100, p) % p
        op = np.ones((1, 1), dtype=np.int64)
        for k in range(n):
            op = np.kron(op, PA[fld[s]] if k == s else eye) % p
        H = (H + h * op) % p
    return H


def jumps_of(n, deph, p):
    return [string_op(tuple(deph[s] if k == s else 0 for k in range(n)), p)
            for s in range(n) if deph[s] != 0]


def nullspace_mod(M, p):
    A = (M % p).astype(np.int64).copy()
    rows, cols = A.shape
    pivots, r = [], 0
    for c in range(cols):
        nz = np.nonzero(A[r:, c])[0]
        if nz.size == 0:
            continue
        piv = r + int(nz[0])
        if piv != r:
            A[[r, piv]] = A[[piv, r]]
        A[r] = (A[r] * inv_p(int(A[r, c]), p)) % p
        others = np.nonzero(A[:, c])[0]
        others = others[others != r]
        if others.size:
            f = A[others, c].copy()
            A[others] = (A[others] - f[:, None] * A[r]) % p
        pivots.append(c)
        r += 1
        if r == rows:
            break
    free = [c for c in range(cols) if c not in pivots]
    out = []
    for fc in free:
        v = np.zeros(cols, dtype=np.int64)
        v[fc] = 1
        for i, pc in enumerate(pivots):
            v[pc] = (-A[i, fc]) % p
        out.append(v)
    return (np.array(out, dtype=np.int64) if out
            else np.zeros((0, cols), dtype=np.int64))


def admissible_U(H, jumps, p, draws=3):
    """(subspace dimension at p, an invertible member or None)."""
    d = H.shape[0]
    eye = np.eye(d, dtype=np.int64)
    blocks = [(np.kron(eye, H.T) - np.kron(H, eye)) % p]
    for A in jumps:
        blocks.append((np.kron(eye, A.T) + np.kron(A, eye)) % p)
    B = nullspace_mod(np.concatenate(blocks, axis=0) % p, p)
    if B.shape[0] == 0:
        return 0, None
    for _ in range(draws):
        U = ((rng.integers(1, p, size=B.shape[0]) @ B) % p).reshape(d, d)
        if det_mod_np_p(U, p) % p:
            return B.shape[0], U
    return B.shape[0], None


def admissible_multi(n, edges, fld, signs, mags, deph, j_num=None,
                     bond_terms=(1, 2, 3)):
    """The criterion verdict, confirmed at every prime where it can be.

    A nullity read mod p can only come out too LARGE, so an EMPTY subspace at
    one prime already proves emptiness over Q(i) and needs no second opinion.
    A NONEMPTY one can be an artifact of the modulus, which is the escape
    channel f138_exact_palindrome_test documents, so it is confirmed at all
    three primes before the row is called admissible.
    """
    dims, U0 = [], None
    for p in PRIMES:
        H = build_H(n, edges, fld, signs, mags, p, j_num=j_num,
                    bond_terms=bond_terms)
        dim, U = admissible_U(H, jumps_of(n, deph, p), p)
        dims.append(dim)
        if dim == 0:
            return 0, None
        if U is None:
            return dim, None
        if p == PRIMES[0]:
            U0 = U
    return dims[0], U0


def vec_transpose(n, p):
    d = 2 ** n
    T = np.zeros((d * d, d * d), dtype=np.int64)
    for i in range(d):
        for j in range(d):
            T[i * d + j, j * d + i] = 1
    return T


def inverse_mod(M, p):
    n = M.shape[0]
    a = np.concatenate([M % p, np.eye(n, dtype=np.int64)], axis=1)
    for c in range(n):
        nz = np.nonzero(a[c:, c])[0]
        assert nz.size, 'singular'
        piv = c + int(nz[0])
        if piv != c:
            a[[c, piv]] = a[[piv, c]]
        a[c] = (a[c] * inv_p(int(a[c, c]), p)) % p
        rows = np.nonzero(a[:, c])[0]
        rows = rows[rows != c]
        if rows.size:
            f = a[rows, c].copy()
            a[rows] = (a[rows] - f[:, None] * a[c]) % p
    return a[:, n:]


def reflector(U, n, p):
    d = 2 ** n
    return (np.kron(U, np.eye(d, dtype=np.int64)) % p) @ vec_transpose(n, p) % p


# --------------------------------------------------------------------------
# gates


def gate1_row():
    print()
    print('## Gate 1: the arc row reproduces, and every control moves')
    print()
    r = ROW
    pred, why = predicate(r['n'], r['edges'], r['deph'], r['fld'])
    pairs = palindromic(r['n'], r['edges'], r['deph'], r['fld'],
                        signs=r['signs'], field_num=r['mags'])
    gate('the F138 predicate forbids the row', pred is False, '(%s)' % why)
    gate('the char-poly palindrome holds anyway', pairs is True)
    for lab, edges, mags in (('generic magnitudes', P3, GEN3),
                             ('equal pair = end+middle', P3, (30, 30, 41)),
                             ('bond+isolate', BI, r['mags'])):
        br = palindromic(3, edges, r['deph'], r['fld'], signs=r['signs'],
                         field_num=mags)
        gate('control breaks: %s' % lab, br is False)
    ok, _ = predicate(3, P3, (0, 1, 0), (0, 3, 0))
    suf = palindromic(3, P3, (0, 1, 0), (0, 3, 0), field_num=r['mags'])
    gate('the sufficient direction still holds', ok and suf)


def gate2_operator():
    """The named U pinned, and the physical reading of what it is."""
    print()
    print('## Gate 2: the named operator U = SWAP_02 . Z_0 Z_1 Z_2')
    print()
    print('  U is not only an algebraic device. [U,H]=0 and A U A = -U give')
    print('  L(U) = -2*sigma*U outright, so U is an exact eigenmode of the')
    print('  Liouvillian at -2*sigma: the mirror partner of the steady state')
    print('  at 0, sitting at the far end of the axis the palindrome folds')
    print('  about. That is what having a U MEANS.')
    print()
    r = ROW
    for p in PRIMES:
        U = site_perm((2, 1, 0), 3, p) @ string_op((3, 3, 3), p) % p
        H = build_H(3, r['edges'], r['fld'], r['signs'], r['mags'], p)
        A = jumps_of(3, r['deph'], p)[0]
        L, shift = build_L(3, r['edges'], r['deph'], r['fld'], p,
                           signs=r['signs'], field_num=r['mags'])
        inv = det_mod_np_p(U, p) % p != 0
        c1 = int(np.count_nonzero((U @ H % p - H @ U % p) % p))
        c2 = int(np.count_nonzero((U @ A % p + A @ U % p) % p))
        vecU = U.reshape(-1) % p
        eig = int(np.count_nonzero((L @ vecU % p + shift * vecU) % p))
        one = np.eye(2 ** 3, dtype=np.int64).reshape(-1) % p
        ker = int(np.count_nonzero(L @ one % p))
        gate('U invertible, [U,H]=0, U A + A U = 0 at p=%d' % p,
             inv and c1 == 0 and c2 == 0)
        gate('L(U) = -2*sigma*U and L(1) = 0 at p=%d' % p,
             eig == 0 and ker == 0)


def gate3_reflector():
    """What U buys, in two steps, both computed."""
    print()
    print('## Gate 3: the metric relation, then the similarity')
    print()
    r = ROW
    for p in PRIMES:
        H = build_H(3, r['edges'], r['fld'], r['signs'], r['mags'], p)
        dim, U = admissible_U(H, jumps_of(3, r['deph'], p), p)
        assert U is not None, 'no admissible U on the named row'
        L, shift = build_L(3, r['edges'], r['deph'], r['fld'], p,
                           signs=r['signs'], field_num=r['mags'])
        eye = np.eye(2 ** 3, dtype=np.int64)
        I2 = np.eye(L.shape[0], dtype=np.int64)
        i_img = sqrt_minus_one(p)
        LH = (((p - i_img) % p)
              * ((np.kron(H, eye) - np.kron(eye, H.T)) % p)) % p
        Ldag = (L - 2 * LH) % p
        LU = np.kron(U, eye) % p
        met = int(np.count_nonzero((((LU @ L % p) @ inverse_mod(LU, p) % p)
                                    - (-Ldag - shift * I2) % p) % p))
        gate('metric relation L_U L L_U^-1 = -L^dag - 2*sigma at p=%d' % p,
             met == 0)
        Pi = reflector(U, 3, p)
        sim = int(np.count_nonzero((((Pi @ L % p) @ inverse_mod(Pi, p) % p)
                                    - (-L - shift * I2) % p) % p))
        gate('similarity Pi L Pi^-1 = -L - 2*sigma at p=%d' % p, sim == 0,
             'subspace dim=%d' % dim)


def gate4_transpose_scope():
    """Where the transpose leg stops, and what carries the palindrome there.

    A Y field makes H non-real, so T no longer turns L into L^dagger and Pi is
    NOT a reflector. The palindrome survives anyway, and the reason cannot be
    seen in GF(p) at all: L preserves hermiticity, so its characteristic
    polynomial is real and its spectrum is conjugation-closed, which collapses
    the metric relation's reflection lambda -> -conj(lambda) - 2*sigma onto the
    palindrome. This gate pins the two halves that ARE modular and states the
    third rather than pretending to measure it.
    """
    print()
    print('## Gate 4: where the transpose leg stops')
    print()
    deph, fld, signs, mags = (0, 0, 1), (0, 0, 2), (1, 1, 1), (30, 30, 30)
    p = PRIME
    H = build_H(3, P3, fld, signs, mags, p)
    gate('the probe row has a NON-symmetric H (a Y field)',
         int(np.count_nonzero((H - H.T) % p)) != 0)
    dim, U = admissible_multi(3, P3, fld, signs, mags, deph)
    gate('an admissible U still exists there', U is not None,
         'subspace dim=%d' % dim)
    if U is None:
        return
    L, shift = build_L(3, P3, deph, fld, p, signs=signs, field_num=mags)
    I2 = np.eye(L.shape[0], dtype=np.int64)
    eye = np.eye(8, dtype=np.int64)
    i_img = sqrt_minus_one(p)
    LH = (((p - i_img) % p) * ((np.kron(H, eye) - np.kron(eye, H.T)) % p)) % p
    Ldag = (L - 2 * LH) % p
    LU = np.kron(U, eye) % p
    met = int(np.count_nonzero((((LU @ L % p) @ inverse_mod(LU, p) % p)
                                - (-Ldag - shift * I2) % p) % p))
    Pi = reflector(U, 3, p)
    sim = int(np.count_nonzero((((Pi @ L % p) @ inverse_mod(Pi, p) % p)
                               - (-L - shift * I2) % p) % p))
    gate('the metric relation still holds', met == 0)
    gate('the similarity Pi FAILS there, as it must', sim != 0,
         'residual entries=%d' % sim)
    gate('and the char-poly palindrome holds regardless',
         palindromic(3, P3, deph, fld, signs=signs, field_num=mags) is True)


def score(label, n, edges, cases, mags, j_num=None, gamma_num=None,
          bond_terms=(1, 2, 3)):
    fp = fn = singular = 0
    ex = []
    for deph, fld, signs in cases:
        dim, U = admissible_multi(n, edges, fld, signs, mags, deph,
                                  j_num=j_num, bond_terms=bond_terms)
        if dim and U is None:
            singular += 1
        pairs = palindromic(n, edges, deph, fld, signs=signs,
                            bond_terms=bond_terms,
                            field_num=mags, j_num=j_num, gamma_num=gamma_num)
        if (U is not None) and not pairs:
            fp += 1
            if len(ex) < 3:
                ex.append(('FALSE POSITIVE', deph, fld, signs))
        elif pairs and U is None:
            fn += 1
            if len(ex) < 3:
                ex.append(('FALSE NEGATIVE', deph, fld, signs))
    detail = 'rows=%d  nonempty-but-singular=%d' % (len(cases), singular)
    gate(label, fp == 0 and fn == 0, detail + '  FP=%d FN=%d' % (fp, fn))
    for tag, deph, fld, signs in ex:
        d = ''.join(LETTER_NAME[a] for a in deph)
        f = ''.join(LETTER_NAME[a] for a in fld)
        print('        %s deph=%s field=%s signs=%s' % (tag, d, f, str(signs)))
    return fp, fn


def grid3():
    return [(d, f, s)
            for d in itertools.product((0, 1, 2, 3), repeat=3)
            for f in itertools.product((0, 1, 2, 3), repeat=3)
            for s in sign_patterns(f)]


def gate5_full_grids():
    print()
    print('## Gate 5: scored on the full N=3 grids, both directions')
    print()
    g = grid3()
    for lab, edges, mags in (('P3, ends equal', P3, (30, 30, 30)),
                             ('K3, all equal', K3, (30, 30, 30)),
                             ('bond+isolate, all equal', BI, (30, 30, 30)),
                             ('P3, generic', P3, GEN3),
                             ('K3, generic', K3, GEN3)):
        score('criterion exact on %s' % lab, 3, edges, g, mags)


def gate6_bond_letters():
    print()
    print('## Gate 6: the bond-letter axis, where the F138 converse fails hardest')
    print()
    print('  F138 records 776 / 732 / 520 converse failures at ONE bond letter')
    print('  and 22 / 104 at two. The criterion is scored on the same axis, and')
    print('  it was never fitted to it.')
    print()
    cases = [(d, f, None)
             for d in itertools.product((0, 1, 2, 3), repeat=3)
             for f in itertools.product((0, 1, 2, 3), repeat=3)]
    for bt, name in (((1, 2, 3), 'XX+YY+ZZ'), ((1, 2), 'XX+YY'),
                     ((1, 3), 'XX+ZZ'), ((2, 3), 'YY+ZZ'),
                     ((1,), 'XX'), ((2,), 'YY'), ((3,), 'ZZ')):
        score('criterion exact at bond %s' % name, 3, P3, cases,
              GEN3, bond_terms=bt)


def gate7_scope():
    print()
    print('## Gate 7: gamma-blindness, the axis ceiling, and N=4')
    print()
    print('  Neither condition on U mentions gamma, so the criterion PREDICTS')
    print('  that the palindrome cannot depend on the gamma profile at all.')
    print('  The non-uniform row is that prediction under test. The two AXIS')
    print('  rows count distinct dephasing LETTERS, which is what the F138')
    print('  ceiling of two counts; the criterion carries no ceiling.')
    print()
    g = grid3()
    sample = [g[i] for i in rng.choice(len(g), 1200, replace=False)]
    score('non-uniform J per bond and gamma per site', 3, P3, sample,
          (30, 30, 30), j_num=[43, 55], gamma_num=[3, 7, 11])
    two_axes = [c for c in sample if len({x for x in c[0] if x}) == 2]
    three_axes = [c for c in sample if len({x for x in c[0] if x}) == 3]
    score('two distinct dephasing AXES', 3, P3, two_axes, (30, 30, 30))
    score('three distinct dephasing AXES, past the F138 ceiling', 3, P3,
          three_axes, (30, 30, 30))
    cases4 = []
    for _ in range(120):
        cases4.append((tuple(int(x) for x in rng.integers(0, 4, size=4)),
                       tuple(int(x) for x in rng.integers(0, 4, size=4)),
                       tuple(int(x) for x in rng.choice([1, -1], size=4))))
    score('N=4 path, equal magnitudes', 4, P4, cases4, (30, 30, 30, 30))
    score('N=4 path, generic magnitudes', 4, P4, cases4, (30, 22, 41, 17))
    score('N=4 ring, equal magnitudes', 4, C4, cases4, (30, 30, 30, 30))
    score('N=4 path, non-uniform J', 4, P4, cases4, (30, 30, 30, 30),
          j_num=[43, 55, 13])


def main():
    print('The pairing condition: one criterion for the dephasing palindrome')
    print('=' * 78)
    print(__doc__.split('\n\n')[1].strip())
    gate1_row()
    gate2_operator()
    gate3_reflector()
    gate4_transpose_scope()
    gate5_full_grids()
    gate6_bond_letters()
    gate7_scope()
    print()
    print('=' * 78)
    bad = [n for n, ok in RESULTS if not ok]
    print('  %d gates, %d passed, %d failed'
          % (len(RESULTS), len(RESULTS) - len(bad), len(bad)))
    for n in bad:
        print('    FAILED: %s' % n)
    print('  ALL GATES PASS' if not bad else '  SOMETHING MOVED, READ IT')


if __name__ == '__main__':
    main()
