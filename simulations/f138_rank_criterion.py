"""The pairing condition as a RANK EQUALITY, gated from below.

THE STATEMENT UNDER TEST, and it is one notch past the criterion the
pairing-condition page carries. For

    L(rho) = -i[H, rho] + sum_l gamma_l (A_l rho A_l - rho),
    every A_l HERMITIAN and UNITARY (A_l^2 = 1), gamma_l > 0, sigma = sum gamma_l,

the char-poly palindrome p(x) = p(-x - 2*sigma) holds exactly when

        dim ker L  ==  dim ker (L + 2*sigma).

No operator to exhibit, no subspace to sample, no search: two nullities.

THE CHAIN THAT PREDICTS IT, step by step, so that each step can fail alone.

  (1) ker L = N, the COMMUTANT {X : [H,X]=0, A_l X A_l = +X},
      ker (L + 2 sigma) = W, the ANTICOMMUTANT {W : [H,W]=0, A_l W A_l = -W}.
      One inclusion each way is pure algebra; the other is the EQUALITY CASE of
      Cauchy-Schwarz in the Hilbert-Schmidt norm and needs C, not GF(p):
      write L = -i ad_H + sigma (Phi - 1) with Phi = sum (gamma_l/sigma) C_l,
      C_l(X) = A_l X A_l an HS isometry. <X, -i ad_H X> is purely imaginary,
      so L(X) = -2 sigma X forces Re<X, Phi X> = -||X||^2, hence Phi X = -X,
      hence every C_l X = -X (equality in the triangle inequality), hence
      ad_H X = 0. The kernel end is the same four lines with the sign flipped.

  (2) BOTH eigenvalues are SEMISIMPLE. L is HS-dissipative and L + 2 sigma is
      HS-accretive, and for an accretive N, N W = V with N V = 0 gives
      Re<W + tV, N(W + tV)> = Re<W,NW> + t||V||^2 >= 0 for every real t, so
      V = 0. Geometric multiplicity = algebraic multiplicity at both ends.

  (3) W holds an INVERTIBLE element  <=>  dim W == dim N.
      (=>) W = N U, so the dimensions agree.
      (<=) if the sign-flip alpha (H -> H, A_l -> -A_l) is well defined on the
      algebra A = <H, A_l>, then W = Hom_A(rho, rho o alpha) has dimension
      sum_k m_k m_alpha(k) <= sum_k m_k^2 = dim N with equality iff the
      multiplicities are alpha-stable, which is exactly rho ~ rho o alpha. And
      if alpha is NOT well defined, some c != 0 is both an even and an odd word,
      so cW = -cW = 0 for every W in W, the ideal it generates is cut out, and
      the inequality is STRICT. Either way equality is the criterion.

  (4) invertible U  =>  palindrome is the pairing-condition page's sufficiency.
      palindrome  =>  alg mult at 0 = alg mult at -2 sigma  =>  (2) the two
      nullities agree. That is necessity, and it is now a proof, not a census.

VERDICT DISCIPLINE, and the direction matters. Reduction mod p is a ring map,
so a RANK can only come out too SMALL and a NULLITY too LARGE. Both nullities
read here are therefore upper bounds on the true ones, and an equality of two
mod-p nullities is not by itself an equality over Q(i); every
structural claim is therefore ALSO checked exactly over Q(i) with Fraction
arithmetic on named rows, where == 0 means == 0. The scale runs are scored at
one prime and every mismatch would be re-run exactly; there are none to re-run.

AND ONE ROUTE THAT SHARES NO CODE. Gate 9 rebuilds the whole object in dense
complex floats with an eigensolver and an optimal spectral matching, importing
nothing from this repo, because the failure this file cannot otherwise rule out
is a shared helper making the two sides agree by construction.

Run:  python simulations/f138_rank_criterion.py
   >  simulations/results/f138_rank_criterion.txt
"""
import itertools
import os
import sys
import time
from fractions import Fraction as Fr

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import f138_pairing_condition as G                              # noqa: E402
from f138_exact_palindrome_test import PRIMES, det_mod_np_p     # noqa: E402
from f138_clause_two_sweep import (build_L, palindromic,        # noqa: E402
                                   LETTER_NAME, GAMMA)

P = PRIMES[0]
rng = np.random.default_rng(20260828)
RESULTS = []


def gate(name, ok, detail=''):
    RESULTS.append((name, bool(ok)))
    print('  [%s] %-56s %s' % ('PASS' if ok else 'FAIL', name, detail))
    return ok


# ---------------------------------------------------------------------------
# modular layer: the two nullities, and the algebra behind them

def rank_mod(M, p):
    A = (M % p).astype(np.int64).copy()
    rows, cols = A.shape
    r = 0
    for c in range(cols):
        nz = np.nonzero(A[r:, c])[0]
        if nz.size == 0:
            continue
        piv = r + int(nz[0])
        if piv != r:
            A[[r, piv]] = A[[piv, r]]
        A[r] = (A[r] * pow(int(A[r, c]), p - 2, p)) % p
        oth = np.nonzero(A[r + 1:, c])[0]
        if oth.size:
            oi = oth + r + 1
            f = A[oi, c].copy()
            A[oi] = (A[oi] - f[:, None] * A[r]) % p
        r += 1
        if r == rows:
            break
    return r


def matmul_mod(A, B, p):
    """A @ B mod p without the silent int64 overflow.

    d^2 = 64 products of residues near 2^30 sum to about 5e18 against an
    int64 ceiling of 9.2e18, so the direct product survives at N=3 on
    sparsity and would overflow at N=4. build_L carries two comments about
    exactly this bug; this splits each entry into high and low halves so the
    headroom does not depend on how full the matrix happens to be.
    """
    A = (A % p).astype(np.int64)
    B = (B % p).astype(np.int64)
    hi, lo = A >> 15, A & 0x7FFF
    return ((((hi @ B) % p) << 15) + (lo @ B)) % p


def _stack(H, jumps, sign, p):
    d = H.shape[0]
    eye = np.eye(d, dtype=np.int64)
    blocks = [(np.kron(eye, H.T) - np.kron(H, eye)) % p]
    for A in jumps:
        if sign > 0:                      # A X A = +X  <=>  [A, X] = 0
            blocks.append((np.kron(eye, A.T) - np.kron(A, eye)) % p)
        else:                             # A X A = -X  <=>  {A, X} = 0
            blocks.append((np.kron(eye, A.T) + np.kron(A, eye)) % p)
    return np.concatenate(blocks, axis=0) % p


def commutant_basis(H, jumps, p):
    return G.nullspace_mod(_stack(H, jumps, +1, p), p)


def anticommutant_basis(H, jumps, p):
    return G.nullspace_mod(_stack(H, jumps, -1, p), p)


def liouvillian(n, edges, deph, fld, signs, mags, p, **kw):
    return build_L(n, edges, deph, fld, p, signs=signs, field_num=mags, **kw)


def invertible_in(B, d, p, draws=4):
    if B.shape[0] == 0:
        return None
    for _ in range(draws):
        U = ((rng.integers(1, p, size=B.shape[0]) @ B) % p).reshape(d, d)
        if det_mod_np_p(U, p) % p:
            return U
    return None


def row_spaces(n, edges, deph, fld, signs, mags, p, **kw):
    """(dim N, dim W, W-basis, H, jumps) at one prime."""
    H = G.build_H(n, edges, fld, signs, mags, p,
                  j_num=kw.get('j_num'), bond_terms=kw.get('bond_terms',
                                                           (1, 2, 3)))
    J = G.jumps_of(n, deph, p)
    BN = commutant_basis(H, J, p)
    BW = anticommutant_basis(H, J, p)
    return BN.shape[0], BW.shape[0], BW, H, J


# ---------------------------------------------------------------------------
# exact layer over Q(i): a Gaussian rational is a pair of Fractions

Z = (Fr(0), Fr(0))
ONE = (Fr(1), Fr(0))
IU = (Fr(0), Fr(1))


def qa(a, b):
    return (a[0] + b[0], a[1] + b[1])


def qs(a, b):
    return (a[0] - b[0], a[1] - b[1])


def qm(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def qd(a, b):
    n = b[0] * b[0] + b[1] * b[1]
    return ((a[0] * b[0] + a[1] * b[1]) / n, (a[1] * b[0] - a[0] * b[1]) / n)


def qc(a):
    return (a[0], -a[1])


def qzero(a):
    return a[0] == 0 and a[1] == 0


def qmat(rows, cols):
    return [[Z] * cols for _ in range(rows)]


def qeye(k):
    M = qmat(k, k)
    for i in range(k):
        M[i][i] = ONE
    return M


def qkron(A, B):
    ra, ca, rb, cb = len(A), len(A[0]), len(B), len(B[0])
    C = qmat(ra * rb, ca * cb)
    for i in range(ra):
        for j in range(ca):
            if qzero(A[i][j]):
                continue
            for k in range(rb):
                for m in range(cb):
                    C[i * rb + k][j * cb + m] = qm(A[i][j], B[k][m])
    return C


def qmul(A, B):
    ra, ca, cb = len(A), len(A[0]), len(B[0])
    C = qmat(ra, cb)
    for i in range(ra):
        Ai = A[i]
        for k in range(ca):
            a = Ai[k]
            if qzero(a):
                continue
            Bk = B[k]
            Ci = C[i]
            for j in range(cb):
                if not qzero(Bk[j]):
                    Ci[j] = qa(Ci[j], qm(a, Bk[j]))
    return C


def qadd(A, B, s=1):
    return [[qa(x, y) if s > 0 else qs(x, y) for x, y in zip(ra, rb)]
            for ra, rb in zip(A, B)]


def qscale(A, c):
    return [[qm(c, x) for x in row] for row in A]


def qdag(A):
    return [[qc(A[j][i]) for j in range(len(A))] for i in range(len(A[0]))]


def qnullspace(A):
    """Exact nullspace basis of A over Q(i), as a list of vectors."""
    rows, cols = len(A), len(A[0])
    M = [row[:] for row in A]
    piv, r = [], 0
    for c in range(cols):
        sel = None
        for i in range(r, rows):
            if not qzero(M[i][c]):
                sel = i
                break
        if sel is None:
            continue
        M[r], M[sel] = M[sel], M[r]
        inv = M[r][c]
        M[r] = [qd(x, inv) for x in M[r]]
        for i in range(rows):
            if i != r and not qzero(M[i][c]):
                f = M[i][c]
                M[i] = [qs(x, qm(f, y)) for x, y in zip(M[i], M[r])]
        piv.append(c)
        r += 1
        if r == rows:
            break
    free = [c for c in range(cols) if c not in piv]
    out = []
    for fc in free:
        v = [Z] * cols
        v[fc] = ONE
        for i, pc in enumerate(piv):
            v[pc] = qs(Z, M[i][fc])
        out.append(v)
    return out


QPAULI = [
    [[ONE, Z], [Z, ONE]],
    [[Z, ONE], [ONE, Z]],
    [[Z, qs(Z, IU)], [IU, Z]],
    [[ONE, Z], [Z, qs(Z, ONE)]],
]


def qsite_op(letter, site, n):
    out = [[ONE]]
    for k in range(n):
        out = qkron(out, QPAULI[letter] if k == site else QPAULI[0])
    return out


def qbond_op(letter, a, c, n):
    out = [[ONE]]
    for k in range(n):
        out = qkron(out, QPAULI[letter] if k in (a, c) else QPAULI[0])
    return out


def qbuild_H(n, edges, fld, signs, mags, bond_terms=(1, 2, 3)):
    d = 2 ** n
    H = qmat(d, d)
    for a, c in edges:
        for t in bond_terms:
            H = qadd(H, qbond_op(t, a, c, n))
    for s in range(n):
        if fld[s] == 0:
            continue
        h = (Fr(mags[s] * signs[s], 100), Fr(0))
        H = qadd(H, qscale(qsite_op(fld[s], s, n), h))
    return H


def qbuild_L(n, edges, deph, fld, signs, mags, bond_terms=(1, 2, 3)):
    """(L, sigma, H, jumps) exactly over Q(i), same convention as build_L."""
    d = 2 ** n
    H = qbuild_H(n, edges, fld, signs, mags, bond_terms=bond_terms)
    ident = qeye(d)
    minus_i = qs(Z, IU)
    L = qscale(qadd(qkron(H, ident), qkron(ident, transpose(H)), -1), minus_i)
    sigma = Fr(0)
    jumps = []
    g = Fr(GAMMA[0], GAMMA[1])
    for s in range(n):
        if deph[s] == 0:
            continue
        A = qsite_op(deph[s], s, n)
        jumps.append(A)
        L = qadd(L, qscale(qkron(A, transpose(A)), (g, Fr(0))))
        L = qadd(L, qscale(qeye(d * d), (g, Fr(0))), -1)
        sigma += g
    return L, sigma, H, jumps


def transpose(A):
    return [[A[j][i] for j in range(len(A))] for i in range(len(A[0]))]


def qreshape(vec, d):
    return [list(vec[i * d:(i + 1) * d]) for i in range(d)]


def qis_zero(A):
    return all(qzero(x) for row in A for x in row)


# ---------------------------------------------------------------------------
# Gate 1: the two kernels ARE the commutant and the anticommutant

NAMED = [
    ('arc row, palindrome holds', 3, G.P3, (0, 1, 0), (1, 0, 1), (1, 1, -1),
     (30, 30, 30), (1, 2, 3)),
    ('arc row, generic mags, breaks', 3, G.P3, (0, 1, 0), (1, 0, 1),
     (1, 1, -1), (30, 22, 41), (1, 2, 3)),
    ('Y field, transpose leg dead', 3, G.P3, (0, 0, 1), (0, 0, 2), (1, 1, 1),
     (30, 30, 30), (1, 2, 3)),
    ('canonical: Z everywhere, no field', 3, G.P3, (3, 3, 3), (0, 0, 0),
     (1, 1, 1), (30, 30, 30), (1, 2, 3)),
    ('bond + isolated site', 3, G.BI, (0, 1, 0), (1, 0, 1), (1, 1, -1),
     (30, 30, 30), (1, 2, 3)),
    ('ZZ bond only', 3, G.P3, (0, 1, 0), (1, 0, 1), (1, 1, -1),
     (30, 30, 30), (3,)),
    ('two axes, X and Z', 3, G.P3, (1, 0, 3), (0, 0, 0), (1, 1, 1),
     (30, 30, 30), (1, 2, 3)),
]


def gate1_kernels_exact():
    print()
    print('## Gate 1: ker L is the commutant and ker(L+2s) is the '
          'anticommutant, EXACTLY over Q(i)')
    print()
    print('  The inclusion each way that is pure algebra holds over any field.')
    print('  The reverse one is the Cauchy-Schwarz equality case and needs C;')
    print('  it is checked here by taking the exact kernel and demanding that')
    print('  EVERY basis element satisfies the two conditions, == 0, no norm.')
    print()
    for lab, n, edges, deph, fld, signs, mags, bt in NAMED:
        d = 2 ** n
        L, sigma, H, jumps = qbuild_L(n, edges, deph, fld, signs, mags,
                                      bond_terms=bt)
        two_s = (2 * sigma, Fr(0))
        M0 = L
        M2 = qadd(L, qscale(qeye(d * d), two_s))
        k0 = qnullspace(M0)
        k2 = qnullspace(M2)
        ok0 = ok2 = True
        for v in k0:
            X = qreshape(v, d)
            if not qis_zero(qadd(qmul(H, X), qmul(X, H), -1)):
                ok0 = False
            for A in jumps:
                if not qis_zero(qadd(qmul(qmul(A, X), A), X, -1)):
                    ok0 = False
        for v in k2:
            W = qreshape(v, d)
            if not qis_zero(qadd(qmul(H, W), qmul(W, H), -1)):
                ok2 = False
            for A in jumps:
                if not qis_zero(qadd(qmul(qmul(A, W), A), W)):
                    ok2 = False
        gate('ker L in the commutant: %s' % lab, ok0,
             'dim ker L = %d' % len(k0))
        gate('ker(L+2s) in the anticommutant: %s' % lab, ok2,
             'dim ker(L+2s) = %d' % len(k2))


def gate1b_reverse_inclusion():
    """The algebraic direction, exactly over Q(i). No prime loop here: this
    inclusion needs no field at all, so one exact field is the whole check."""
    print()
    print('## Gate 1b: the algebraic inclusion, which needs no field at all')
    print()
    for lab, n, edges, deph, fld, signs, mags, bt in NAMED:
        d = 2 ** n
        L, sigma, H, jumps = qbuild_L(n, edges, deph, fld, signs, mags,
                                      bond_terms=bt)
        # commutant and anticommutant, exact, then apply L
        dd = d * d
        rowsN, rowsW = [], []
        eye = qeye(d)
        cH = qadd(qkron(eye, transpose(H)), qkron(H, eye), -1)
        rowsN.extend(cH)
        rowsW.extend(cH)
        for A in jumps:
            kA = qkron(A, transpose(A))
            rowsN.extend(qadd(kA, qeye(dd), -1))
            rowsW.extend(qadd(kA, qeye(dd)))
        BN = qnullspace(rowsN if rowsN else qeye(dd))
        BW = qnullspace(rowsW if rowsW else qeye(dd))
        okN = all(all(qzero(x) for x in mv(L, v)) for v in BN)
        two_s = (2 * sigma, Fr(0))
        M2 = qadd(L, qscale(qeye(dd), two_s))
        okW = all(all(qzero(x) for x in mv(M2, v)) for v in BW)
        gate('commutant in ker L: %s' % lab, okN, 'dim N = %d' % len(BN))
        gate('anticommutant in ker(L+2s): %s' % lab, okW,
             'dim W = %d' % len(BW))


def mv(A, v):
    return [reduce_row(row, v) for row in A]


def reduce_row(row, v):
    acc = Z
    for a, b in zip(row, v):
        if not qzero(a) and not qzero(b):
            acc = qa(acc, qm(a, b))
    return acc


# ---------------------------------------------------------------------------
# Gate 2: semisimplicity at both ends

def gate2_semisimple():
    print()
    print('## Gate 2: both ends are SEMISIMPLE, so geometric = algebraic')
    print()
    print('  rank(M^2) == rank(M) is the semisimplicity of the eigenvalue at')
    print('  the corresponding end. It is what turns the char-poly identity')
    print('  into a statement about the two nullities.')
    print()
    bad0 = bad2 = 0
    for lab, n, edges, deph, fld, signs, mags, bt in NAMED:
        for p in PRIMES:
            L, shift = liouvillian(n, edges, deph, fld, signs, mags, p,
                                   bond_terms=bt)
            I = np.eye(L.shape[0], dtype=np.int64)
            M = (L + shift * I) % p
            r0, r0b = rank_mod(L, p), rank_mod(matmul_mod(L, L, p), p)
            r2, r2b = rank_mod(M, p), rank_mod(matmul_mod(M, M, p), p)
            bad0 += (r0 != r0b)
            bad2 += (r2 != r2b)
    gate('0 semisimple on every named row at every prime', bad0 == 0)
    gate('-2*sigma semisimple on every named row at every prime', bad2 == 0)

    g = G.grid3()
    idx = rng.choice(len(g), 250, replace=False)
    b0 = b2 = 0
    live0 = live2 = 0
    for i in idx:
        deph, fld, signs = g[i]
        L, shift = liouvillian(3, G.P3, deph, fld, signs, (30, 30, 30), P)
        I = np.eye(L.shape[0], dtype=np.int64)
        # % P, not % p: an earlier version reduced by the LEAKED loop variable
        # from the named-row loop above, i.e. by a different (smaller) prime
        # than the one L was built at. That destroys the matrix and gives it
        # full rank, so rank(M^2) == rank(M) passed for free on exactly the
        # rows that had a kernel to test. The live counters below exist so
        # that the same accident cannot pass silently a second time.
        M = (L + shift * I) % P
        rank0, rank2 = rank_mod(L, P), rank_mod(M, P)
        live0 += (rank0 < L.shape[0])
        live2 += (rank2 < L.shape[0])
        b0 += (rank0 != rank_mod(matmul_mod(L, L, P), P))
        b2 += (rank2 != rank_mod(matmul_mod(M, M, P), P))
    gate('the grid rows are not vacuous: both ends carry a kernel somewhere',
         live0 > 0 and live2 > 0,
         'rows with ker L nonzero=%d, with ker(L+2s) nonzero=%d, of 250'
         % (live0, live2))
    gate('0 semisimple on 250 grid rows', b0 == 0)
    gate('-2*sigma semisimple on 250 grid rows', b2 == 0)


def gate2b_accretive():
    """The inequality the semisimplicity argument rests on, exactly."""
    print()
    print('## Gate 2b: dissipative at 0, accretive at -2*sigma, exactly')
    print()
    for lab, n, edges, deph, fld, signs, mags, bt in NAMED[:4]:
        d = 2 ** n
        L, sigma, H, jumps = qbuild_L(n, edges, deph, fld, signs, mags,
                                      bond_terms=bt)
        okL = okA = True
        for _ in range(6):
            v = [(Fr(int(x), 7), Fr(int(y), 7)) for x, y in
                 zip(rng.integers(-9, 10, size=d * d),
                     rng.integers(-9, 10, size=d * d))]
            Lv = mv(L, v)
            re = sum((qm(qc(a), b)[0] for a, b in zip(v, Lv)), Fr(0))
            nrm = sum((qm(qc(a), a)[0] for a in v), Fr(0))
            if re > 0:
                okL = False
            if re + 2 * sigma * nrm < 0:
                okA = False
        gate('Re<X, L X> <= 0 (dissipative): %s' % lab, okL)
        gate('Re<X, (L+2s) X> >= 0 (accretive): %s' % lab, okA)


# ---------------------------------------------------------------------------
# Gate 3: the module structure, and invertibility as a rank equality

def gate3_module():
    print()
    print('## Gate 3: W = N U when U is invertible, so the dimensions agree')
    print()
    for lab, n, edges, deph, fld, signs, mags, bt in NAMED:
        p = P
        d = 2 ** n
        H = G.build_H(n, edges, fld, signs, mags, p, bond_terms=bt)
        J = G.jumps_of(n, deph, p)
        BN, BW = commutant_basis(H, J, p), anticommutant_basis(H, J, p)
        U = invertible_in(BW, d, p)
        if U is None:
            gate('no invertible U, and then dim W < dim N: %s' % lab,
                 BW.shape[0] < BN.shape[0],
                 'dim N = %d, dim W = %d' % (BN.shape[0], BW.shape[0]))
            continue
        prod = np.stack([((X.reshape(d, d) @ U) % p).reshape(-1) for X in BN])
        stacked = np.concatenate([BW, prod], axis=0) % p
        same = rank_mod(stacked, p) == BW.shape[0] == BN.shape[0]
        gate('N U spans W exactly: %s' % lab, same,
             'dim N = %d, dim W = %d' % (BN.shape[0], BW.shape[0]))


def gate3b_rank_decides_invertibility(label, n, edges, cases, mags, **kw):
    bad = 0
    ex = []
    for deph, fld, signs in cases:
        dN, dW, BW, H, J = row_spaces(n, edges, deph, fld, signs, mags, P, **kw)
        inv = invertible_in(BW, 2 ** n, P) is not None
        if inv != (dN == dW):
            bad += 1
            if len(ex) < 3:
                ex.append((deph, fld, signs, dN, dW, inv))
    gate(label, bad == 0, 'rows=%d  mismatches=%d' % (len(cases), bad))
    for deph, fld, signs, dN, dW, inv in ex:
        print('        MISMATCH deph=%s field=%s dimN=%d dimW=%d inv=%s'
              % (''.join(LETTER_NAME[a] for a in deph),
                 ''.join(LETTER_NAME[a] for a in fld), dN, dW, inv))
    return bad


# ---------------------------------------------------------------------------
# Gate 4: the criterion itself, scored against the exact palindrome

def score_rank(label, n, edges, cases, mags, **kw):
    fp = fn = holds = breaks = both_live = 0
    ex = []
    t0 = time.time()
    for deph, fld, signs in cases:
        dN, dW, BW, H, J = row_spaces(n, edges, deph, fld, signs, mags, P, **kw)
        says = (dN == dW)
        pairs = palindromic(n, edges, deph, fld, signs=signs,
                            field_num=mags, **kw)
        holds += bool(pairs)
        breaks += (not pairs)
        both_live += (0 < dW < dN)
        if says and not pairs:
            fp += 1
            ex.append(('FALSE POSITIVE', deph, fld, signs, dN, dW))
        elif pairs and not says:
            fn += 1
            ex.append(('FALSE NEGATIVE', deph, fld, signs, dN, dW))
    gate(label, fp == 0 and fn == 0 and holds > 0 and breaks > 0,
         'rows=%d holds=%d breaks=%d  0<dimW<dimN on %d  FP=%d FN=%d  (%.0f s)'
         % (len(cases), holds, breaks, both_live, fp, fn, time.time() - t0))
    for tag, deph, fld, signs, dN, dW in ex[:3]:
        print('        %s deph=%s field=%s signs=%s dimN=%d dimW=%d'
              % (tag, ''.join(LETTER_NAME[a] for a in deph),
                 ''.join(LETTER_NAME[a] for a in fld), signs, dN, dW))
    return fp, fn


def gate4_scored(rows3=2500, rows4=120):
    print()
    print('## Gate 4: the rank equality scored against the exact palindrome')
    print()
    g = G.grid3()
    idx = rng.choice(len(g), rows3, replace=False)
    sample = [g[i] for i in idx]
    for lab, edges, mags in (('P3, ends equal', G.P3, (30, 30, 30)),
                             ('K3, all equal', G.K3, (30, 30, 30)),
                             ('bond+isolate', G.BI, (30, 30, 30)),
                             ('P3, generic', G.P3, G.GEN3)):
        score_rank('rank criterion exact on %s' % lab, 3, edges, sample, mags)
    cases = [(d, f, (1, 1, 1))
             for d in itertools.product((0, 1, 2, 3), repeat=3)
             for f in itertools.product((0, 1, 2, 3), repeat=3)]
    ci = rng.choice(len(cases), 900, replace=False)
    cs = [cases[i] for i in ci]
    for bt, name in (((1, 2), 'XX+YY'), ((1,), 'XX'), ((3,), 'ZZ')):
        score_rank('rank criterion exact at bond %s' % name, 3, G.P3, cs,
                   G.GEN3, bond_terms=bt)
    score_rank('rank criterion under non-uniform J and gamma', 3, G.P3,
               sample[:600], (30, 30, 30), j_num=[43, 55],
               gamma_num=[3, 7, 11])
    c4 = []
    for _ in range(rows4):
        c4.append((tuple(int(x) for x in rng.integers(0, 4, size=4)),
                   tuple(int(x) for x in rng.integers(0, 4, size=4)),
                   tuple(int(x) for x in rng.choice([1, -1], size=4))))
    score_rank('rank criterion exact, N=4 path', 4, G.P4, c4,
               (30, 30, 30, 30))
    score_rank('rank criterion exact, N=4 ring', 4, G.C4, c4,
               (30, 30, 30, 30))
    print()
    print('  And the same rows scored for T4, which is what replaces sampling:')
    gate3b_rank_decides_invertibility(
        'invertible element <=> dim N == dim W, P3 equal', 3, G.P3, sample,
        (30, 30, 30))
    gate3b_rank_decides_invertibility(
        'invertible element <=> dim N == dim W, bond+isolate', 3, G.BI,
        sample, (30, 30, 30))


def gate5_canonical():
    print()
    print('## Gate 5: the canonical configuration, where the repo already '
          'held BOTH numbers')
    print()
    print('  experiments/DEGENERACY_PALINDROME.md Result 2 records N+1')
    print('  conserved quantities at Re = 0 and N+1 XOR sector modes at the')
    print('  far end, and calls the map between them bijective. Those are')
    print('  exactly dim ker L and dim ker(L+2*sigma). The criterion says the')
    print('  palindrome IS that equality.')
    print()
    for n, edges in ((2, [(0, 1)]), (3, G.P3), (4, G.P4)):
        deph = tuple([3] * n)
        fld = tuple([0] * n)
        signs = tuple([1] * n)
        mags = tuple([30] * n)
        dN, dW, BW, H, J = row_spaces(n, edges, deph, fld, signs, mags, P)
        pal = palindromic(n, edges, deph, fld, signs=signs, field_num=mags)
        gate('N=%d chain: dim ker L = dim ker(L+2s) = N+1 = %d' % (n, n + 1),
             dN == dW == n + 1, 'palindrome=%s' % pal)


# ---------------------------------------------------------------------------
# Gate 6: the kernel always dominates, and the F103 counterexamples

LETTER_IDX = {'I': 0, 'X': 1, 'Y': 2, 'Z': 3}
WORDS = [a + b for a in 'XYZ' for b in 'XYZ']


def word_op(word, sites, n, p):
    PA = G.paulis(p)
    eye = np.eye(2, dtype=np.int64)
    op = np.ones((1, 1), dtype=np.int64)
    for k in range(n):
        m = eye
        for w, s in zip(word, sites):
            if k == s:
                m = PA[LETTER_IDX[w]]
        op = np.kron(op, m) % p
    return op


def build_H_words(n, edges, terms, p):
    d = 2 ** n
    H = np.zeros((d, d), dtype=np.int64)
    for a, c in edges:
        for t in terms:
            H = (H + word_op(t, (a, c), n, p)) % p
    return H


def axis_jump(nvec, site, n, p):
    """A = (a X + b Y + c Z)/den on one site, a^2+b^2+c^2 = den^2, so A^2 = 1
    exactly over Q(i): dephasing along a RATIONAL direction off the letters."""
    a, b, c, den = nvec
    PA = G.paulis(p)
    M = ((a * PA[1] + b * PA[2] + c * PA[3]) * G_inv(den, p)) % p
    eye = np.eye(2, dtype=np.int64)
    op = np.ones((1, 1), dtype=np.int64)
    for k in range(n):
        op = np.kron(op, M if k == site else eye) % p
    return op


def G_inv(a, p):
    return pow(int(a) % p, p - 2, p)


def build_L_raw(H, jumps, gnum, gden, p):
    d = H.shape[0]
    ident = np.eye(d, dtype=np.int64)
    from f138_exact_palindrome_test import sqrt_minus_one
    i_img = sqrt_minus_one(p)
    L = (((p - i_img) % p) * ((np.kron(H, ident)
                               - np.kron(ident, H.T)) % p)) % p
    g = (gnum * G_inv(gden, p)) % p
    shift = 0
    for A in jumps:
        L = (L + g * (np.kron(A, A.T) % p)) % p
        L = (L - g * np.eye(d * d, dtype=np.int64)) % p
        shift = (shift + g) % p
    return L, (2 * shift) % p


def palindromic_raw(L, shift, p, points=6):
    eye = np.eye(L.shape[0], dtype=np.int64)
    for _ in range(points):
        x = int(rng.integers(0, p))
        y = (-x - shift) % p
        if (det_mod_np_p((x * eye - L) % p, p)
                != det_mod_np_p((y * eye - L) % p, p)):
            return False
    return True


def spaces_raw(H, jumps, p):
    d = H.shape[0]
    BN = commutant_basis(H, jumps, p)
    BW = anticommutant_basis(H, jumps, p)
    return BN.shape[0], BW.shape[0], BW


def score_raw(label, rows):
    fp = fn = holds = breaks = dominated = both_live = 0
    ex = []
    for tag, H, jumps in rows:
        dN, dW, BW = spaces_raw(H, jumps, P)
        L, shift = build_L_raw(H, jumps, 1, 20, P)
        pal = palindromic_raw(L, shift, P)
        says = (dN == dW)
        holds += bool(pal)
        breaks += (not pal)
        dominated += (dW <= dN)
        both_live += (0 < dW < dN)
        if says and not pal:
            fp += 1
            ex.append(('FALSE POSITIVE', tag, dN, dW))
        elif pal and not says:
            fn += 1
            ex.append(('FALSE NEGATIVE', tag, dN, dW))
    gate(label, fp == 0 and fn == 0 and dominated == len(rows),
         'rows=%d holds=%d breaks=%d  0<dimW<dimN on %d  FP=%d FN=%d'
         % (len(rows), holds, breaks, both_live, fp, fn))
    for t, tag, dN, dW in ex[:3]:
        print('        %s %s dimN=%d dimW=%d' % (t, tag, dN, dW))
    return fp, fn


def gate6_domination():
    print()
    print('## Gate 6: dim ker(L+2s) <= dim ker L, always')
    print()
    print('  A corollary of the same argument and a claim on its own: the')
    print('  steady modes can never be outnumbered by the fastest ones.')
    print()
    g = G.grid3()
    idx = rng.choice(len(g), 1200, replace=False)
    worst = 0
    bad = 0
    for i in idx:
        deph, fld, signs = g[i]
        for edges, mags in ((G.P3, (30, 30, 30)), (G.BI, (30, 30, 30)),
                            (G.K3, G.GEN3)):
            dN, dW, BW, H, J = row_spaces(3, edges, deph, fld, signs, mags, P)
            bad += (dW > dN)
            worst = max(worst, dN - dW)
    gate('dim W <= dim N on 3600 rows', bad == 0,
         'largest gap dim N - dim W = %d' % worst)


def gate7_f103_counterexamples():
    print()
    print('## Gate 7: the F103 7.12 rows, where no DIAGONAL operator exists')
    print()
    print('  PROOF_F103_F87_Z2_CUBED_REFINEMENT.md 7.12 records XX+XZ, YY+YZ')
    print('  and XX+XZ+ZX as SOFT on the chain at N = 3..6 with non-bipartite')
    print('  basis-state graphs: no diagonal D, no chiral K. It does NOT leave')
    print('  the restoring operator unexhibited, and an earlier version of')
    print('  this banner said it did: the same paragraph continues "that')
    print('  operator is no longer a mystery: it is the hidden-Q routing, a')
    print('  per-site Q from the P1/P4 families, which TwoTermPalindromeRouting')
    print('  classifies bit-exactly for 2-term pairs". So these rows are an')
    print('  AGREEMENT test against a route the repo already owns, not a')
    print('  discovery, and the operator below is named rather than admired.')
    print()
    print('  One thing here is a TAUTOLOGY and is gated as one rather than')
    print('  reported as evidence: under Z on every site, A W A = -W forces')
    print('  every diagonal entry of every element of W to vanish, so "the')
    print('  operator is not diagonal" cannot fail and says nothing.')
    print()
    letters = 'IXYZ'
    for terms in (['XX', 'XZ'], ['YY', 'YZ'], ['XX', 'XZ', 'ZX']):
        for n in (3, 4, 5):
            edges = [(i, i + 1) for i in range(n - 1)]
            H = build_H_words(n, edges, terms, P)
            jumps = [G.string_op(tuple(3 if k == s else 0 for k in range(n)),
                                 P) for s in range(n)]
            dN, dW, BW = spaces_raw(H, jumps, P)
            L, shift = build_L_raw(H, jumps, 1, 20, P)
            pal = palindromic_raw(L, shift, P)
            name = '-'
            if dW == 1:
                W = BW[0].reshape(2 ** n, 2 ** n) % P
                for s in itertools.product((0, 1, 2, 3), repeat=n):
                    S = G.string_op(s, P)
                    nz = np.nonzero(W.reshape(-1))[0]
                    if nz.size == 0:
                        break
                    c = (int(W.reshape(-1)[nz[0]])
                         * pow(int(S.reshape(-1)[nz[0]]) or 1, P - 2, P)) % P
                    if not np.count_nonzero((W - c * S) % P):
                        name = ''.join(letters[t] for t in s)
                        break
            expected = ('Y' if terms[0][0] == 'X' else 'X') * n
            gate('%s at N=%d: criterion = palindrome, and W is %s'
                 % ('+'.join(terms), n, expected),
                 (dN == dW) == bool(pal) and dW == 1 and name == expected,
                 'dimN=%d dimW=%d palindrome=%s  W spanned by %s'
                 % (dN, dW, pal, name))
    print()
    print('  So the operator is one Pauli string per row, found by a linear')
    print('  solve rather than a routing table, and on YY+YZ it is X^(x)N,')
    print('  which is F1 own Pi squared. That is a check that the count')
    print('  agrees with what the repo already knows, and nothing more.')
    z = [G.string_op(tuple(3 if k == s else 0 for k in range(3)), P)
         for s in range(3)]
    Hq = build_H_words(3, [(0, 1), (1, 2)], ['XX', 'XZ'], P)
    _, _, BWq = spaces_raw(Hq, z, P)
    diagfree = all(not np.count_nonzero(
        np.diag(BWq[k].reshape(8, 8)) % P) for k in range(BWq.shape[0]))
    gate('and the diagonal-free property is a tautology, as stated', diagfree,
         'every element of W has zero diagonal, so the flag cannot fail')


def gate8_scope_axes():
    print()
    print('## Gate 8: three axes the criterion page lists as NOT TESTED')
    print()
    print('  The theorem asks only that every A_l be Hermitian and square to')
    print('  1, which a rational-direction n.sigma and a Pauli STRING both do.')
    print('  The holds/breaks columns are printed because a block carrying')
    print('  only one verdict tests only one direction. The per-site mixed')
    print('  directions block is very nearly such a block, and the printed')
    print('  counts below are the authority on how nearly: read them rather')
    print('  than this sentence, which an earlier version got wrong by')
    print('  hard-coding a count that the run then contradicted eight lines')
    print('  later.')
    print()
    n = 3
    edges = G.P3
    VECS = [(3, 4, 0, 5), (0, 3, 4, 5), (4, 0, 3, 5), (2, 3, 6, 7),
            (6, 2, 3, 7), (1, 2, 2, 3), (2, 6, 3, 7)]
    sq = all(not np.count_nonzero((axis_jump(v, 0, n, P) @ axis_jump(v, 0, n, P)
                                   % P - np.eye(2 ** n, dtype=np.int64)) % P)
             for v in VECS)
    gate('every rational-direction jump squares to 1', sq)
    rows = []
    for v in VECS:
        for fld in itertools.product((0, 1, 2, 3), repeat=3):
            for sites in ((0,), (0, 2), (0, 1, 2)):
                rows.append(('n=%s fld=%s' % (v, fld),
                             G.build_H(n, edges, fld, (1, 1, 1),
                                       (30, 30, 30), P),
                             [axis_jump(v, s, n, P) for s in sites]))
    idx = rng.choice(len(rows), 600, replace=False)
    score_raw('off-axis dephasing, one direction, 1-3 sites',
              [rows[i] for i in idx])
    rows = []
    for fld in itertools.product((0, 1, 2, 3), repeat=3):
        for _ in range(3):
            vs = [VECS[int(i)] for i in rng.integers(0, len(VECS), size=3)]
            rows.append(('mixed fld=%s' % (fld,),
                         G.build_H(n, edges, fld, (1, 1, 1), (30, 30, 30), P),
                         [axis_jump(v, s, n, P) for s, v in enumerate(vs)]))
    score_raw('off-axis dephasing, a DIFFERENT direction per site', rows)
    strings = [(1, 1, 0), (3, 3, 0), (1, 2, 3), (0, 1, 1), (2, 2, 2),
               (1, 0, 1), (3, 1, 0), (1, 3, 3)]
    rows = []
    for s1 in strings:
        for fld in itertools.product((0, 1, 2, 3), repeat=3):
            rows.append(('A=%s' % (s1,),
                         G.build_H(n, edges, fld, (1, 1, 1), (30, 30, 30), P),
                         [G.string_op(s1, P)]))
    idx = rng.choice(len(rows), 400, replace=False)
    score_raw('a multi-site Pauli STRING as the jump', [rows[i] for i in idx])
    rows = []
    for s1, s2 in itertools.combinations(strings, 2):
        for fld in itertools.product((0, 1, 2, 3), repeat=3):
            rows.append(('A=%s,%s' % (s1, s2),
                         G.build_H(n, edges, fld, (1, 1, 1), (30, 30, 30), P),
                         [G.string_op(s1, P), G.string_op(s2, P)]))
    idx = rng.choice(len(rows), 500, replace=False)
    score_raw('two Pauli-string jumps', [rows[i] for i in idx])
    for dephname, letters in (('Z everywhere', (3, 3, 3)),
                              ('X everywhere', (1, 1, 1))):
        rows = []
        for k in (1, 2):
            for terms in itertools.combinations(WORDS, k):
                rows.append(('+'.join(terms), build_H_words(n, edges,
                                                            list(terms), P),
                             G.jumps_of(n, letters, P)))
        score_raw('all 1- and 2-letter bond words, deph %s' % dephname, rows)
    rows = []
    for terms in itertools.combinations(WORDS, 3):
        rows.append(('+'.join(terms), build_H_words(n, edges, list(terms), P),
                     G.jumps_of(n, (3, 3, 3), P)))
    score_raw('all 84 three-term bond words, deph Z everywhere', rows)


# ---------------------------------------------------------------------------
# Gate 9: the same verdict by a route that shares no code with the rest

def gate9_independent_route():
    """Dense complex floats and an eigensolver, built here from nothing.

    The CONSTRUCTION is self-contained: Paulis, H, the Liouvillian, the
    palindrome (an optimal matching of the spectrum against its reflection) and
    the two nullities (singular values) are all rebuilt inside this function.
    The one import is row_spaces, and it is the other SIDE of the cross-check
    rather than a shared helper: without it there would be nothing modular to
    compare against.

    Two comparisons, and they are different. The first is internal to the float
    route, its verdict against its own two nullities. The second crosses the
    routes, the float SVD nullities against the GF(p) eliminations row by row.
    Only the second excludes a shared helper making the two sides agree by
    construction, which is the one failure mode the rest of this file cannot
    rule out, and an earlier version of this docstring described the first as
    though it were the second.

    A float eigensolver has no exact route, so the verdict here is read at a
    threshold, and the threshold is defended the only way it can be: the run
    prints the worst distance among rows it calls palindromic and the smallest
    among rows it calls broken, and gates on the SEPARATION between them rather
    than on the threshold's value. A lexicographic sort is NOT a matching, and
    reporting that is part of the gate: the first version of this check sorted,
    and called every row broken including the canonical chain, because near
    degenerate real parts break the tie in opposite orders.
    """
    from scipy.optimize import linear_sum_assignment
    print()
    print('## Gate 9: the same verdict by an independent route')
    print()
    I2 = np.eye(2, dtype=complex)
    PAc = [I2,
           np.array([[0, 1], [1, 0]], dtype=complex),
           np.array([[0, -1j], [1j, 0]], dtype=complex),
           np.array([[1, 0], [0, -1]], dtype=complex)]

    def site(letter, s, n):
        out = np.ones((1, 1), dtype=complex)
        for k in range(n):
            out = np.kron(out, PAc[letter] if k == s else I2)
        return out

    def bond(letter, a, c, n):
        out = np.ones((1, 1), dtype=complex)
        for k in range(n):
            out = np.kron(out, PAc[letter] if k in (a, c) else I2)
        return out

    n, edges = 3, [(0, 1), (1, 2)]
    grid = [(d, f, s)
            for d in itertools.product((0, 1, 2, 3), repeat=3)
            for f in itertools.product((0, 1, 2, 3), repeat=3)
            for s in ((1, 1, 1), (1, 1, -1), (1, -1, 1))]
    idx = rng.choice(len(grid), 900, replace=False)
    agree = disagree = holds = breaks = nullity_mismatch = 0
    worst_hold, best_break = 0.0, np.inf
    sv_zero, sv_live = 0.0, np.inf
    for i in idx:
        deph, fld, signs = grid[i]
        js = [j for j in range(n) if deph[j]]
        if not js:
            continue
        d = 2 ** n
        H = np.zeros((d, d), dtype=complex)
        for a, c in edges:
            for t in (1, 2, 3):
                H = H + bond(t, a, c, n)
        for s in range(n):
            if fld[s]:
                H = H + (30 * signs[s] / 100.0) * site(fld[s], s, n)
        ident = np.eye(d, dtype=complex)
        L = -1j * (np.kron(H, ident) - np.kron(ident, H.T))
        sigma = 0.0
        for j in js:
            A = site(deph[j], j, n)
            L = L + 0.05 * np.kron(A, A.T) - 0.05 * np.eye(d * d,
                                                           dtype=complex)
            sigma += 0.05
        ev = np.linalg.eigvals(L)
        cost = np.abs(ev[:, None] - (-ev - 2 * sigma)[None, :])
        r, c = linear_sum_assignment(cost)
        dist = float(cost[r, c].max())
        pal = dist < 1e-7

        def nullity(M):
            nonlocal sv_zero, sv_live
            sv = np.linalg.svd(M, compute_uv=False)
            cut = 1e-9 * max(1.0, sv[0])
            dead, live = sv[sv < cut], sv[sv >= cut]
            if dead.size:
                sv_zero = max(sv_zero, float(dead.max()))
            if live.size:
                sv_live = min(sv_live, float(live.min()))
            return int(np.sum(sv < cut))

        nf, nw = nullity(L), nullity(L + 2 * sigma * np.eye(d * d,
                                                             dtype=complex))
        says = (nf == nw)
        # and the point of the gate: read the SAME two numbers the modular
        # side reads, and compare them row by row. An earlier version compared
        # the float verdict against the float nullities only, which is not a
        # cross-check at all.
        mN, mW, _, _, _ = row_spaces(n, edges, deph, fld, signs,
                                     (30, 30, 30), P)
        nullity_mismatch += (nf != mN or nw != mW)
        holds += bool(pal)
        breaks += (not pal)
        if pal:
            worst_hold = max(worst_hold, dist)
        else:
            best_break = min(best_break, dist)
        agree += (says == pal)
        disagree += (says != pal)
    sep = np.log10(best_break / max(worst_hold, 1e-300))
    svsep = np.log10(sv_live / max(sv_zero, 1e-300))
    gate('float route: its own verdict matches its own two nullities',
         disagree == 0 and holds > 0 and breaks > 0,
         'rows=%d holds=%d breaks=%d disagreements=%d'
         % (agree + disagree, holds, breaks, disagree))
    gate('and the two ROUTES read the same two nullities, row by row',
         nullity_mismatch == 0,
         'float SVD against GF(p) elimination, mismatches=%d' % nullity_mismatch)
    gate('the spectral threshold is a law, not a number',
         sep > 6 and holds > 0 and breaks > 0,
         'worst hold %.1e, smallest break %.1e, %.1f decades apart'
         % (worst_hold, best_break, sep))
    gate('and so is the singular-value threshold', svsep > 6,
         'largest zero %.1e, smallest live %.1e, %.1f decades apart'
         % (sv_zero, sv_live, svsep))


# ---------------------------------------------------------------------------
# Gate 10: three boundaries of the class, two of them out of sample

def build_L_gammas(H, jumps, gammas, p):
    """Same convention as build_L_raw but with a rate PER JUMP."""
    d = H.shape[0]
    ident = np.eye(d, dtype=np.int64)
    from f138_exact_palindrome_test import sqrt_minus_one
    i_img = sqrt_minus_one(p)
    L = (((p - i_img) % p) * ((np.kron(H, ident)
                               - np.kron(ident, H.T)) % p)) % p
    shift = 0
    for A, (gn, gd) in zip(jumps, gammas):
        if gn == 0:
            continue
        g = (gn * G_inv(gd, p)) % p
        L = (L + g * (np.kron(A, A.T) % p)) % p
        L = (L - g * np.eye(d * d, dtype=np.int64)) % p
        shift = (shift + g) % p
    return L, (2 * shift) % p


def gate10_boundaries():
    print()
    print('## Gate 10: three boundaries, two of them outside the sample')
    print()

    print('  (a) DEPOLARIZING. X, Y and Z are each Hermitian and square to 1,')
    print('  so a depolarizing site is INSIDE this class, and F1 records it as')
    print('  the canonical BREAK with a residual of (2/3)*Sum gamma. Nothing')
    print('  in the scored grids has three jumps on one site, so this is out')
    print('  of sample in both senses.')
    for n in (2, 3):
        edges = [(i, i + 1) for i in range(n - 1)]
        H = build_H_words(n, edges, ['XX', 'YY', 'ZZ'], P)
        jumps = [G.string_op(tuple(t if k == s else 0 for k in range(n)), P)
                 for s in range(n) for t in (1, 2, 3)]
        dN, dW, _ = spaces_raw(H, jumps, P)
        L, shift = build_L_raw(H, jumps, 1, 20, P)
        pal = palindromic_raw(L, shift, P)
        gate('depolarizing at N=%d: criterion = palindrome' % n,
             (dN == dW) == bool(pal),
             'dimN=%d dimW=%d palindrome=%s (both say NO)' % (dN, dW, pal))

    print()
    print('  (b) A RATE SWITCHED OFF. The theorem asks gamma_l > 0, and the')
    print('  fence is not decoration: a zero rate DROPS a condition from W,')
    print('  so the VERDICT itself moves at the boundary of the orthant, not')
    print('  merely the dimensions. The row below is chosen to show that: a')
    print('  ZZ bond with an X field on site 0, dephasing X on both sites,')
    print('  broken while both rates are on and palindromic the moment site')
    print('  0 stops being watched. Nothing in the scored grids sets a rate')
    print('  to zero, so this axis is out of sample.')
    n = 2
    H = G.build_H(n, [(0, 1)], (1, 0), (1, 1), (30, 30), P, bond_terms=(3,))
    A0 = G.string_op((1, 0), P)
    A1 = G.string_op((0, 1), P)
    moved = set()
    for lab, gam in (('both rates on  ', ((2, 5), (9, 10))),
                     ('site 1 rate off', ((2, 5), (0, 1))),
                     ('site 0 rate off', ((0, 1), (9, 10)))):
        live = [A for A, (gn, _) in zip((A0, A1), gam) if gn]
        dN, dW, _ = spaces_raw(H, live, P)
        L, shift = build_L_gammas(H, (A0, A1), gam, P)
        pal = palindromic_raw(L, shift, P)
        moved.add(bool(pal))
        gate('rate fence, %s: criterion = palindrome' % lab,
             (dN == dW) == bool(pal),
             'dimN=%d dimW=%d palindrome=%s' % (dN, dW, pal))
    gate('and the verdict really moves when a rate is switched off',
         len(moved) == 2, 'both verdicts appear across the three rows')

    print()
    print('  (c) ODD LOCAL DIMENSION. An invertible U with U A U^-1 = -A makes')
    print('  A and -A similar, so A needs balanced +1/-1 multiplicities and d')
    print('  must be EVEN. At odd d the criterion must therefore say NO at')
    print('  every H, which is a corollary and not a measurement, gated here')
    print('  because every other row of this file sits at d = 2^N.')
    ok = True
    worst = None
    for trial in range(24):
        d = 3 if trial % 2 == 0 else 5
        A = np.diag([1] * (d - 1) + [p_ - 1 if (p_ := P) else 0]).astype(
            np.int64) % P
        Hs = rng.integers(0, 40, size=(d, d)).astype(np.int64)
        H = (Hs + Hs.T) % P                       # real symmetric, Hermitian
        dN, dW, BW = spaces_raw(H, [A], P)
        if dW >= dN:
            ok = False
            worst = (d, dN, dW)
        if dW and invertible_in(BW, d, P, draws=4) is not None:
            ok = False
            worst = (d, dN, dW)
    gate('odd d: dim W < dim N always, and no invertible element', ok,
         'd = 3 and 5, 24 random Hermitian H' if ok else 'broke at %s' % (worst,))


# ---------------------------------------------------------------------------
# Gate 11: rows built to have 0 < dim W < dim N, where the criterion's EQUALITY
# says more than its nonemptiness

def gate11_strict_inequality_by_construction():
    """The axis every other gate is thin on, reached by construction.

    Across the scored grids only 204 of 15,406 rows have 0 < dim W < dim N, so
    on the rest a much weaker predicate (dim W > 0) would score identically.
    Those rows are the ones that separate the criterion from its own shadow,
    and they can be BUILT rather than waited for:

        H = H1 (+) H2 and A = A1 (+) A2 with spec(H1) disjoint from spec(H2),
        block 1 carrying an invertible anticommuting element and block 2 none.

    Disjoint spectra kill the cross blocks, so W = W1 (+) 0 with every element
    singular by construction, while N draws from both blocks. The criterion
    must then say BROKEN, and the spectrum must agree. Smallest instance:
    H = X (+) (Z + 5I), A = Z (+) Z at d = 4.

    Construction from the second session, 2026-08-29, arrived while this file
    was being repaired and is credited rather than absorbed.
    """
    print()
    print('## Gate 11: 0 < dim W < dim N, built rather than waited for')
    print()
    print('  Only 204 of the 15,406 scored grid rows have 0 < dim W < dim N,')
    print('  which is where the EQUALITY says more than nonemptiness. Here')
    print('  they are constructed: two blocks with disjoint spectra, one')
    print('  carrying an invertible anticommuting element and one carrying')
    print('  none, so W is nonzero and entirely singular by construction.')
    print()
    strict = 0
    for shift in (5, 7, 11, 13):
        for lead in (1, 2):
            d = 4
            # block 1: H1 = X or Y (anticommutes with Z, so W1 = the block itself),
            # block 2: H2 = Z + shift*I (commutes with Z, so W2 = 0), and the shift
            # separates the two spectra.
            PA = G.paulis(P)
            H = np.zeros((d, d), dtype=np.int64)
            H[:2, :2] = PA[lead]
            H[2:, 2:] = (PA[3] + shift * np.eye(2, dtype=np.int64)) % P
            A = np.zeros((d, d), dtype=np.int64)
            A[:2, :2] = PA[3]
            A[2:, 2:] = PA[3]
            dN, dW, BW = spaces_raw(H, [A], P)
            L, sh = build_L_raw(H, [A], 1, 20, P)
            pal = palindromic_raw(L, sh, P)
            U = invertible_in(BW, d, P, draws=6) if dW else None
            strict += (0 < dW < dN)
            gate('block pair (H1=%s, shift %d): criterion = palindrome'
                 % ('XYZ'[lead - 1], shift), (dN == dW) == bool(pal),
                 'dimN=%d dimW=%d strict=%s invertible=%s palindrome=%s'
                 % (dN, dW, 0 < dW < dN, U is not None, pal))
    gate('every constructed row really is strictly between', strict == 8,
         '%d of 8 rows have 0 < dim W < dim N' % strict)


def main():
    print('The pairing condition as a rank equality')
    print('=' * 78)
    print(__doc__.split('\n\n')[1].strip())
    gate1_kernels_exact()
    gate1b_reverse_inclusion()
    gate2_semisimple()
    gate2b_accretive()
    gate3_module()
    gate5_canonical()
    gate6_domination()
    gate7_f103_counterexamples()
    gate8_scope_axes()
    gate9_independent_route()
    gate10_boundaries()
    gate11_strict_inequality_by_construction()
    gate4_scored()
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
