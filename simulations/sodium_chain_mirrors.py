"""Chains of sodium nuclei: the two mirrors of one spin 3/2 under a coupling between nuclei.

Verifier for experiments/THE_SODIUM_CHAIN.md. Nucleus k is the two qubits 2k, 2k+1 in Souza et al.'s labelling
(|3/2>, |1/2>, |-1/2>, |-3/2> = |00>, |01>, |10>, |11>), its relaxation Pauli-string dephasing with the jumps ZZ at
3J0/2, ZX and ZY at 3J1/2, XI and YI at 3J2/2 (experiments/THE_SODIUM_NUCLEUS_TWO_MIRRORS.md), an independent
gradient on every nucleus. S(rho) = (⊗ZZ_k)·rho^T.

Exactness. H is taken times 4 and the rates as even integers, so every entry of L is a small Gaussian integer and
float64 arithmetic on it is exact; such residuals are compared with == 0. The sqrt(3) part of the dipolar coupling is
carried as its own integer matrix. Over GF(p), p = 1 mod 12 (i and sqrt(3) exist, so the reduction of Z[i, sqrt3] is
a ring map), statements of three strengths are read: a characteristic polynomial of a matrix defined over the ring
that is not palindromic modulo p is not palindromic before the reduction (exact); a dimension of an invariant subspace
over GF(p) is an upper bound on the complex one; a palindromic characteristic polynomial modulo p, and any reading on
a subspace computed modulo p, is evidence, not a proof (the failure of K_0 is proved separately, G9). Two primes must
agree. Where the page proves a statement from the lemma (the ladder, the pairing on K_m), the GF(p) reading follows
from the exactly gated lemma (G4) and checks the code, not the mathematics.

The model: identical nuclei (one splitting, one offset) for the equal-splitting truncation and the full coupling, the
full coupling at zero splitting, where it is the Hamiltonian; two splittings and one offset for the unequal
truncation; two splittings and two offsets for the Ising coupling, which is what a coupling reduces to between nuclei
of different resonance frequencies. The common rates J0 >= J1 >= J2, the order a motional spectral density gives;
the triples per nucleus (G4b, G7b) are arbitrary.

G10 reads the full coupling at two nuclei: the invariant part of the two-sided kernel of its sqrt3 part, and what
K_0 holds beyond it; the float spectra of that part are printed as a reading.

G11 puts two nuclei under a shared gradient (one Kossakowski matrix g·[[1, c], [c, 1]] per channel, everything times
the denominator of c): the price identity with an operator price K, the conserved sign pattern t without a
Hamiltonian, the one sector (-,-,-,-) where the price is a number, and over GF(p) the breaks elsewhere.

Run: python simulations/sodium_chain_mirrors.py [--n3]   (n = 2 and 3 nuclei; --n3 adds the sharp subspaces at n = 3, slow)
"""
from __future__ import annotations

import sys
from itertools import product

import numpy as np
import sympy as sp

PRIMES = (10009, 10069)          # both 1 mod 12
RESULTS = []


def check(name, ok):
    RESULTS.append((name, bool(ok)))
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)


# ---------------------------------------------------------------- strings

PA = {'I': np.eye(2), 'X': np.array([[0, 1], [1, 0]]), 'Y': np.array([[0, -1j], [1j, 0]]), 'Z': np.diag([1, -1])}


def string(word):
    m = np.eye(1)
    for c in word:
        m = np.kron(m, PA[c])
    return m.astype(complex)


def anticommutes(a, b):
    return sum(1 for x, y in zip(a, b) if x != 'I' and y != 'I' and x != y) % 2 == 1


def place(nq, at, letters):
    w = ['I'] * nq
    for i, c in enumerate(letters):
        w[at + i] = c
    return ''.join(w)


# ---------------------------------------------------------------- G1: the dictionary, exact

def gate_dictionary():
    SP = {'I': sp.eye(2), 'X': sp.Matrix([[0, 1], [1, 0]]), 'Y': sp.Matrix([[0, -sp.I], [sp.I, 0]]),
          'Z': sp.Matrix([[1, 0], [0, -1]])}

    def s(w):
        m = sp.Matrix([[1]])
        for c in w:
            m = sp.kronecker_product(m, SP[c])
        return m

    ms = [sp.Rational(3, 2), sp.Rational(1, 2), -sp.Rational(1, 2), -sp.Rational(3, 2)]
    Iz = sp.diag(*ms)
    Ip = sp.zeros(4, 4)
    for k in range(1, 4):
        Ip[k - 1, k] = sp.sqrt(sp.Rational(15, 4) - ms[k] * (ms[k] + 1))
    Ix, Iy = (Ip + Ip.T) / 2, (Ip - Ip.T) / (2 * sp.I)
    r3 = sp.sqrt(3)
    zero = sp.zeros(4, 4)
    check('G1 I_x = (sqrt3/2)·IX + (XX + YY)/2', sp.simplify(Ix - (r3 / 2 * s('IX') + (s('XX') + s('YY')) / 2)) == zero)
    check('G1 I_y = (sqrt3/2)·IY + (YX - XY)/2', sp.simplify(Iy - (r3 / 2 * s('IY') + (s('YX') - s('XY')) / 2)) == zero)
    check('G1 I_z = ZI + IZ/2', Iz - (s('ZI') + s('IZ') / 2) == zero)
    check('G1 T20 = (3/sqrt6)·ZZ', sp.simplify((3 * Iz ** 2 - sp.Rational(15, 4) * sp.eye(4)) / sp.sqrt(6)
                                               - 3 / sp.sqrt(6) * s('ZZ')) == zero)
    check('G1 control: I_x without the sqrt3 part is wrong', sp.simplify(Ix - (s('XX') + s('YY')) / 2) != zero)
    # the secular dipolar coupling of two nuclei, 3 Iz Iz - I.I, in strings (times 4)
    k = sp.kronecker_product
    dip = 4 * (2 * k(Iz, Iz) - k(Ix, Ix) - k(Iy, Iy))
    terms = {}
    for w in product('IXYZ', repeat=4):
        w = ''.join(w)
        c = sp.nsimplify(sp.simplify((s(w).H * dip).trace() / 16))
        if c != 0:
            terms[w] = c
    rational = {w: int(c) for w, c in terms.items() if c.is_rational}
    root = {w: int(sp.nsimplify(c / r3)) for w, c in terms.items() if not c.is_rational}
    check('G1 the coupling splits into a rational part and sqrt3 times an integer part',
          all(sp.simplify(terms[w] - (rational.get(w, 0) + r3 * root.get(w, 0))) == 0 for w in terms))
    check('G1 22 strings, 8 of them with sqrt3', len(terms) == 22 and len(root) == 8)
    return rational, root


# ---------------------------------------------------------------- the generator

def lindbladian(n, coupling, J, q=(3,), w=(5,), root=None, r3=0.0):
    """L in the row-major vec basis; H times 4: q_k·4·ZZ_k + w_k·(4·ZI + 2·IZ) per nucleus (q_k the coefficient of ZZ,
    w_k that of I_z, cycled over the nuclei), the coupling (a dict of 4-qubit strings, times 4) between neighbours, the
    sqrt3 part (root) at the factor r3 (0: left out). J is one triple (J0, J1, J2) for every nucleus, or a tuple of
    triples, one per nucleus. The rates must be even, so that 3J/2 is an integer."""
    Js = list(J) if isinstance(J[0], tuple) else [tuple(J)] * n
    assert len(Js) == n and all(j % 2 == 0 for t in Js for j in t), J
    nq, dim = 2 * n, 4 ** n
    H = np.zeros((dim, dim), complex)
    for k in range(n):
        H += 4 * q[k % len(q)] * string(place(nq, 2 * k, 'ZZ'))
        H += w[k % len(w)] * (4 * string(place(nq, 2 * k, 'ZI')) + 2 * string(place(nq, 2 * k, 'IZ')))
    for k in range(n - 1):
        for s_, c in coupling.items():
            H += c * string(place(nq, 2 * k, s_))
        if root:
            for s_, c in root.items():
                H += r3 * c * string(place(nq, 2 * k, s_))
    Id = np.eye(dim)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    LJ0 = np.zeros_like(L)
    for k in range(n):
        J0, J1, J2 = Js[k]
        for s_, g in (('ZZ', 3 * J0 // 2), ('ZX', 3 * J1 // 2), ('ZY', 3 * J1 // 2), ('XI', 3 * J2 // 2), ('YI', 3 * J2 // 2)):
            P = string(place(nq, 2 * k, s_))
            D = g * (np.kron(P, P.conj()) - np.eye(dim * dim))
            L += D
            if s_ == 'ZZ':
                LJ0 += D
    return L, LJ0


def s_perm(n):
    """S(rho) = Q·rho^T as a signed permutation of the vec basis: (S v)[a] = sign[a]·v[perm[a]]"""
    dim = 4 ** n
    Qd = np.real(np.diag(string('ZZ' * n)))
    perm = np.array([(a % dim) * dim + a // dim for a in range(dim * dim)])      # (i, j) <- (j, i)
    sign = np.array([Qd[a // dim] for a in range(dim * dim)])
    return perm, sign


def conj_S(L, perm, sign):
    """S L S^-1 for the signed permutation S"""
    return (sign[:, None] * sign[None, :]) * L[np.ix_(perm, perm)]


def bright(n):
    """per nucleus the diagonal of P_k on the cells: 1 where the cell is charged by ZZ_k"""
    nq, dim = 2 * n, 4 ** n
    out = []
    for k in range(n):
        z = np.real(np.diag(string(place(nq, 2 * k, 'ZZ'))))
        out.append(((1 - np.kron(z, z)) / 2).astype(int))
    return out


# ---------------------------------------------------------------- GF(p)

def to_mod(M, p):
    re, im = np.rint(M.real).astype(np.int64), np.rint(M.imag).astype(np.int64)
    assert np.abs(M.real - re).max() == 0 and np.abs(M.imag - im).max() == 0
    i_p = next(x for x in range(p) if x * x % p == p - 1)
    return (re + im * i_p) % p


def rref_mod(A, p):
    A = A.copy() % p
    rows, cols = A.shape
    r, pivots = 0, []
    for c in range(cols):
        if r == rows:
            break
        nz = np.nonzero(A[r:, c])[0]
        if nz.size == 0:
            continue
        pr = r + nz[0]
        A[[r, pr]] = A[[pr, r]]
        A[r] = A[r] * pow(int(A[r, c]), p - 2, p) % p
        others = np.nonzero(A[:, c])[0]
        others = others[others != r]
        if others.size:
            A[others] = (A[others] - np.outer(A[others, c], A[r])) % p
        pivots.append(c)
        r += 1
    return A[:r], pivots


def null_mod(A, p, cols):
    R, piv = rref_mod(A, p) if A.shape[0] else (np.zeros((0, cols), np.int64), [])
    free = [c for c in range(cols) if c not in piv]
    B = np.zeros((cols, len(free)), np.int64)
    for j, f in enumerate(free):
        B[f, j] = 1
        for i, c in enumerate(piv):
            B[c, j] = (-R[i, f]) % p
    return B


def det_mod(A, p):
    A = A.copy() % p
    n = A.shape[0]
    d = 1
    for c in range(n):
        nz = np.nonzero(A[c:, c])[0]
        if nz.size == 0:
            return 0
        pr = c + nz[0]
        if pr != c:
            A[[c, pr]] = A[[pr, c]]
            d = -d
        d = d * int(A[c, c]) % p
        inv = pow(int(A[c, c]), p - 2, p)
        rows = np.arange(c + 1, n)
        f = A[rows, c] * inv % p
        A[rows] = (A[rows] - np.outer(f, A[c])) % p
    return d % p


def palindromic_mod(M, c, p):
    """det(x - M) == (-1)^k det(-x - c - M) at k + 1 points over GF(p): the spectrum is symmetric about -c/2"""
    k = M.shape[0]
    I = np.eye(k, dtype=np.int64)
    for x in range(1, k + 2):
        a = det_mod(x * I - M, p)
        b = det_mod(((-x - c) % p) * I - M, p) * (1 if k % 2 == 0 else -1) % p
        if a != b:
            return False
    return True


def sharp_subspace(Lp, E, O, p):
    """basis (over GF(p)) of the largest L-invariant subspace inside the coordinate subspace E"""
    C = Lp[np.ix_(O, E)]
    LEE = Lp[np.ix_(E, E)]
    C, _ = rref_mod(C, p)
    while True:
        Cn, _ = rref_mod(np.vstack([C, C @ LEE % p]), p)
        if Cn.shape[0] == C.shape[0]:
            break
        C = Cn
    return null_mod(C, p, len(E))


def restricted(Lp, E, B, p):
    """M with L B = B M over GF(p) (B spans an L-invariant subspace)"""
    LB = Lp[np.ix_(E, E)] @ B % p
    chosen = []                                 # k independent rows of B, greedily
    for i in range(B.shape[0]):
        trial = chosen + [i]
        if rref_mod(B[trial], p)[0].shape[0] == len(trial):
            chosen = trial
        if len(chosen) == B.shape[1]:
            break
    Bs = B[chosen]
    k = Bs.shape[0]
    aug = np.hstack([Bs, LB[chosen]]) % p
    R, _ = rref_mod(aug, p)
    return R[:, k:]


# ---------------------------------------------------------------- main

def main():
    sys.stdout.reconfigure(encoding='utf-8')
    rational, root = gate_dictionary()
    check('G2 the sqrt3 terms are exactly the terms anticommuting with ZZ_1·ZZ_2 (satellite x central)',
          set(root) == {w for w in set(rational) | set(root) if anticommutes(w, 'ZZZZ')})
    each = {w: c for w, c in rational.items() if not (anticommutes(w, 'ZZII') or anticommutes(w, 'IIZZ'))}
    check('G2 the rational terms anticommuting with each ZZ_k are IXIX and IYIY',
          set(rational) - set(each) == {'IXIX', 'IYIY'})

    # G3: the secular truncations, against the projection onto the energy shells of q1·ZZ_1 + q2·ZZ_2
    full_rat = sum(c * string(w) for w, c in rational.items())
    full_root = sum(c * string(w) for w, c in root.items())
    extras = {'IXIX': -3, 'IYIY': -3, 'ZXZX': -3, 'ZYZY': -3}

    def projection(M, q1, q2):
        h = np.real(np.diag(q1 * string('ZZII') + q2 * string('IIZZ')))
        return sum(np.diag((h == v).astype(float)) @ M @ np.diag((h == v).astype(float)) for v in np.unique(h))

    equal2 = {w: 2 * each.get(w, 0) + extras.get(w, 0) for w in set(each) | set(extras)}
    check('G3 equal splitting: twice the secular part = twice the terms commuting with each ZZ_k, '
          'plus IXIX, IYIY, ZXZX, ZYZY at -3',
          np.abs(2 * projection(full_rat, 3, 3) - sum(c * string(w) for w, c in equal2.items())).max() == 0)
    check('G3 equal splitting: the sqrt3 part has no secular part', np.abs(projection(full_root, 3, 3)).max() == 0)
    uneq_form = sum(c * string(w) for w, c in each.items())
    check('G3 unequal splitting: the secular part is the rational terms commuting with each ZZ_k',
          np.abs(projection(full_rat, 3, 5) - uneq_form).max() == 0 and np.abs(projection(full_root, 3, 5)).max() == 0)
    check('G3 control: the two truncations differ', np.abs(projection(full_rat, 3, 3) - uneq_form).max() > 0)

    # (coupling, the splittings q_k, the offsets w_k)
    couplings = {
        'Ising': ({'ZIZI': 4, 'ZIIZ': 2, 'IZZI': 2, 'IZIZ': 1}, (3, 5), (5, 7)),
        'dipolar, unequal splittings': (each, (3, 5), (5,)),
        'dipolar, equal splittings (times 2)': (equal2, (3,), (5,)),
    }
    J = (6, 4, 2)
    for n in (2, 3):
        perm, sign = s_perm(n)
        P = bright(n)
        c = 6 * n * (J[1] + J[2])
        for name, (cp, q, w) in couplings.items():
            L, LJ0 = lindbladian(n, cp, J, q=q, w=w)
            Lneg, _ = lindbladian(n, cp, (-J[0], J[1], J[2]), q=q, w=w)
            res = conj_S(L, perm, sign) + Lneg + c * np.eye(L.shape[0])
            check(f'G4 n={n} {name}: S·L(J0)·S^-1 = -L(-J0) - 6n(J1+J2) exactly', np.abs(res).max() == 0)
            keeps = all(np.abs(L[np.ix_(Pk == 1, Pk == 0)]).max() == 0 and np.abs(L[np.ix_(Pk == 0, Pk == 1)]).max() == 0
                        for Pk in P)
            expect = not name.startswith('dipolar, equal')
            check(f'G5 n={n} {name}: every nucleus keeps its brightness ([L, P_k] = 0) is {expect}', keeps == expect)
        # the full dipolar coupling: its sqrt3 part alone already breaks the identity (value-free: any nonzero factor)
        Lr, _ = lindbladian(n, {}, (0, 0, 0), root=root, r3=1.0)
        L0, _ = lindbladian(n, {}, (0, 0, 0))
        Lroot = Lr - L0
        check(f'G4 n={n} control: the sqrt3 part of the dipolar coupling breaks S (S·L_r·S^-1 + L_r != 0)',
              np.abs(conj_S(Lroot, perm, sign) + Lroot).max() > 0)
        Q = 'ZZ' * n
        jumps = [place(2 * n, 2 * k, j) for k in range(n) for j in ('ZX', 'ZY', 'XI', 'YI')]
        lit = [''.join(w) for w in product('IXYZ', repeat=2 * n) if all(anticommutes(''.join(w), j) for j in jumps)]
        check(f'G6 n={n} at J0 = 0 the strings anticommuting with every jump are Q = ⊗ZZ_k alone', lit == [Q])
        jumps0 = jumps + [place(2 * n, 2 * k, 'ZZ') for k in range(n)]
        lit0 = [''.join(w) for w in product('IXYZ', repeat=2 * n) if all(anticommutes(''.join(w), j) for j in jumps0)]
        check(f'G6 n={n} control: with the J0 jumps present no string anticommutes with every jump', lit0 == [])
        Hr = sum(c * string(place(2 * n, 2 * k, w)) for k in range(n - 1) for w, c in root.items())
        Qm = string(Q)
        check(f'G6 n={n} F158 at J0 = 0, value-free: [H, ⊗ZZ_k] != 0 for the full coupling, from its sqrt3 part alone '
              f'(the rest commutes), so the far end is empty',
              np.abs(Hr @ Qm - Qm @ Hr).max() > 0)
        # the exact K_0 route: on the cells with no bright nucleus the dissipator's most negative entry, -2 s_0, sits on
        # Q alone; the identity (eigenvalue 0) would need a partner there, and Q is no eigenvector since [H, Q] != 0
        Ddiag = {''.join(w): -2 * sum(g for j, g in zip(jumps, [3 * J[1] // 2, 3 * J[1] // 2, 3 * J[2] // 2, 3 * J[2] // 2] * n)
                                      if anticommutes(''.join(w), j))
                 for w in product('IXYZ', repeat=2 * n)
                 if not any(anticommutes(''.join(w), place(2 * n, 2 * k, 'ZZ')) for k in range(n))}
        low = min(Ddiag.values())
        check(f'G9 n={n} on the cells with no bright nucleus the dissipator reaches -2 s_0 = {low} on Q alone',
              low == -6 * n * (J[1] + J[2]) and [w for w, v in Ddiag.items() if v == low] == [Q])

    # G6b, a second route at n = 2: at J0 = 0 the full coupling's characteristic polynomial is not palindromic about
    # -sigma (a failure modulo p is a failure, exact), the equal-splitting truncation's is (evidence)
    n = 2
    J00 = (0, J[1], J[2])
    p = PRIMES[0]
    r3p = next(x for x in range(p) if x * x % p == 3)
    Lfull = (to_mod(lindbladian(n, rational, J00, q=(0,))[0], p)
             + r3p * to_mod(lindbladian(n, {}, J00, q=(0,), root=root, r3=1.0)[0] - lindbladian(n, {}, J00, q=(0,))[0], p)) % p
    sig2 = 6 * n * (J[1] + J[2])
    check('G6b n=2 the full coupling at J0 = 0 is not palindromic about -sigma (exact)', not palindromic_mod(Lfull, sig2, p))
    check('G6b n=2 control: the equal-splitting truncation at J0 = 0 is (proved by the lemma; a code check)',
          palindromic_mod(to_mod(lindbladian(n, equal2, J00)[0], p), sig2, p))

    # G4b/G7b: rates per nucleus (nuclei of different sites or isotopes): the lemma with 6·sum_k(J1_k + J2_k), the
    # ladder indexed by the SET B of bright nuclei, 2 s_B = 6·sum_k(J1_k + J2_k) + 6·sum_{k in B} J0_k
    Jk = ((6, 4, 2), (10, 2, 4), (4, 8, 6))
    for nn in (2, 3):
        perm, sign = s_perm(nn)
        Js = Jk[:nn]
        cp, q, w = couplings['Ising']
        L, _ = lindbladian(nn, cp, Js, q=q, w=w)
        Lneg, _ = lindbladian(nn, cp, tuple((-a, b, c_) for a, b, c_ in Js), q=q, w=w)
        cc = 6 * sum(b + c_ for _, b, c_ in Js)
        check(f'G4b n={nn} Ising, rates per nucleus: S·L(J0)·S^-1 = -L(-J0) - 6·sum(J1_k + J2_k) exactly',
              np.abs(conj_S(L, perm, sign) + Lneg + cc * np.eye(L.shape[0])).max() == 0)
    Js = Jk[:2]
    P = bright(2)
    L, _ = lindbladian(2, couplings['Ising'][0], Js, q=couplings['Ising'][1], w=couplings['Ising'][2])
    cc = 6 * sum(b + c_ for _, b, c_ in Js)
    for p in PRIMES:
        ok = True
        for pattern in product((0, 1), repeat=2):
            idx = np.where(np.all([P[k] == b for k, b in enumerate(pattern)], axis=0))[0]
            ok &= palindromic_mod(to_mod(L, p)[np.ix_(idx, idx)], cc + 6 * sum(Js[k][0] for k in range(2) if pattern[k]), p)
        check(f'G7b n=2 Ising, rates per nucleus: every pattern block B palindromic about -s_B over GF({p})', ok)
    idx = np.where(np.all([P[k] == b for k, b in enumerate((1, 0))], axis=0))[0]
    check('G7b control: the block B = {0} is not palindromic about -s_{1} (the other bright nucleus)',
          not palindromic_mod(to_mod(L, PRIMES[0])[np.ix_(idx, idx)], cc + 6 * Js[1][0], PRIMES[0]))

    # G7: the ladder at n = 2, a code check (it follows from G4/G5): the characteristic polynomial of every bright pattern block
    P = bright(n)
    for name in ('Ising', 'dipolar, unequal splittings'):
        cp, q, w = couplings[name]
        L, _ = lindbladian(n, cp, J, q=q, w=w)
        for p in PRIMES:
            Lp = to_mod(L, p)
            ok = True
            for pattern in product((0, 1), repeat=n):
                idx = np.where(np.all([P[k] == b for k, b in enumerate(pattern)], axis=0))[0]
                m = sum(pattern)
                ok &= palindromic_mod(Lp[np.ix_(idx, idx)], 6 * m * J[0] + 6 * n * (J[1] + J[2]), p)
            check(f'G7 n=2 {name}: every bright pattern block palindromic about -s_m over GF({p})', ok)
        idx = np.where(np.all([P[k] == b for k, b in enumerate((1, 0))], axis=0))[0]
        check(f'G7 n=2 {name} control: the m = 1 block is not palindromic about -s_0',
              not palindromic_mod(to_mod(L, PRIMES[0])[np.ix_(idx, idx)], 6 * n * (J[1] + J[2]), PRIMES[0]))

    # G8: the sharp subspaces K_m (equal splittings; the full coupling at zero splitting), n = 2 at two rate points,
    # n = 3 with --n3
    for Jp in (J, (8, 6, 2)):
        sharp_gate(2, rational, root, couplings, Jp,
                   {'dipolar, equal splittings (times 2)': ((24, 128, 24), (True, True, True)),
                    'full dipolar': ((10, 12, 1), (False, True, True))})
    kernel_gate(rational, root, J)
    shared_gate()
    if '--n3' in sys.argv:
        sharp_gate(3, rational, root, couplings, J,
                   {'dipolar, equal splittings (times 2)': ((64, 140, 84, 16), (True, True, True, True)),
                    'full dipolar': ((12, 16, 0, 0), (False, True, None, None))})

    fails = [nm for nm, ok in RESULTS if not ok]
    print(f'\n{len(RESULTS) - len(fails)}/{len(RESULTS)} PASS')
    return 1 if fails else 0


def kernel_gate(rational, root, J):
    """G10, the full coupling at two nuclei: what carries the readings. The premises of the page's proposition exactly
    (H3, the sqrt3 part, is real and anticommutes with Q; rank H3 = 4 of 16), then U_m, the largest L-invariant subspace
    inside the operators x with H3·x = 0 = x·H3 and m bright nuclei, over GF(p): its dimensions and its palindrome
    (a code check of the proposition, which holds without the reduction)."""
    n, nq, dim = 2, 4, 16
    H3 = sum(c * string(place(nq, 0, w)) for w, c in root.items())
    Qm = string('ZZ' * n)
    check('G10 H3 is symmetric and anticommutes with Q (realness: next check)',
          np.abs(H3 - H3.T).max() == 0 and np.abs(Qm @ H3 @ Qm + H3).max() == 0)
    check('G10 rank H3 = 4 of 16 (exact)', sp.Matrix(np.rint(H3.real).astype(int)).rank() == 4
          and np.abs(H3.imag).max() == 0)
    W = sum(bright(n))
    dims, pal = {}, {}
    for p in PRIMES:
        r3p = next(x for x in range(p) if x * x % p == 3)
        Lp = (to_mod(lindbladian(n, rational, J, q=(0,))[0], p)
              + r3p * to_mod(lindbladian(n, {}, J, q=(0,), root=root, r3=1.0)[0] - lindbladian(n, {}, J, q=(0,))[0], p)) % p
        H3p = to_mod(H3, p)
        left = np.kron(H3p, np.eye(dim, dtype=np.int64)) % p           # x -> H3·x on the row-major vec
        right = np.kron(np.eye(dim, dtype=np.int64), H3p.T) % p         # x -> x·H3
        dims[p], pal[p] = [], []
        for m in range(n + 1):
            E, O = np.where(W == m)[0], np.where(W != m)[0]
            C, _ = rref_mod(np.vstack([np.eye(dim * dim, dtype=np.int64)[O], left, right]) % p, p)
            while True:
                Cn, _ = rref_mod(np.vstack([C, C @ Lp % p]) % p, p)
                if Cn.shape[0] == C.shape[0]:
                    break
                C = Cn
            B = null_mod(C, p, dim * dim)
            dims[p].append(B.shape[1])
            if m == 0:
                MU = restricted(Lp, E, B[E], p)
                if p == PRIMES[0]:
                    check('G10 control: the part with no bright nucleus is not palindromic about -s_1 (a code check)',
                          not palindromic_mod(MU, 6 * J[0] + 6 * n * (J[1] + J[2]), p))
                # what K_0 holds beyond U_0: char(K_0) = char(U_0)·x(x + 48)^2((x + 36)^2 + 1600) at J = (6, 4, 2)
                Bk = sharp_subspace(Lp, E, O, p)
                MK = restricted(Lp, E, Bk, p)
                ok = J == (6, 4, 2)
                for x in range(1, MK.shape[0] + 2):
                    lhs = det_mod(x * np.eye(MK.shape[0], dtype=np.int64) - MK, p)
                    rhs = det_mod(x * np.eye(MU.shape[0], dtype=np.int64) - MU, p) * x * (x + 48) ** 2 * ((x + 36) ** 2 + 1600) % p
                    ok &= lhs == rhs
                check(f'G10 K_0 beyond U_0 holds exactly 0, a double -48 and the pair -36 ± 40i over GF({p}) (evidence)', ok)
            pal[p].append(palindromic_mod(restricted(Lp, E, B[E], p), 6 * m * J[0] + 6 * n * (J[1] + J[2]), p)
                          if B.shape[1] else None)
    check('G10 n=2 the invariant part inside the two-sided kernel of H3 has dimensions (5, 12, 1) over both primes '
          '(upper bounds)', all(tuple(dims[p]) == (5, 12, 1) for p in PRIMES))
    check('G10 n=2 it is palindromic about -s_m on every m over both primes (proved; a code check)',
          all(tuple(pal[p]) == (True, True, True) for p in PRIMES))
    # a reading, printed, not gated: the float spectra of U_m (the dimensions from GF(p)), which show where the
    # mirror does more than conjugation: an eigenvalue off Re = -s_m needs a partner of another real part
    L = lindbladian(n, rational, J, q=(0,), root=root, r3=np.sqrt(3))[0]
    for m, k in zip(range(n + 1), dims[PRIMES[0]]):
        Pm = np.diag((W == m).astype(float))
        _, sv, vh = np.linalg.svd(np.vstack([np.kron(H3, np.eye(dim)), np.kron(np.eye(dim), H3.T), np.eye(dim * dim) - Pm]))
        B = vh.conj().T[:, -int(np.sum(sv < 1e-9) + (dim * dim - len(sv))):] if len(sv) < dim * dim else vh.conj().T[:, sv < 1e-9]
        for _ in range(100):
            if B.shape[1] <= k:
                break
            R = L @ B - B @ (B.conj().T @ (L @ B))
            _, s2, v2 = np.linalg.svd(R)
            s2 = np.concatenate([s2, np.zeros(max(0, B.shape[1] - len(s2)))])
            B = B @ v2.conj().T[:, np.argsort(s2)[:max(k, int(np.sum(s2 < 1e-8 * max(1, s2.max()))))]]
        ev = np.sort_complex(np.linalg.eigvals(B.conj().T @ L @ B))
        s_m = 3 * m * J[0] + 3 * n * (J[1] + J[2])
        off = [complex(round(e.real, 3), round(e.imag, 3)) for e in ev if abs(e.real + s_m) > 1e-6]
        inv = np.abs(L @ B - B @ (B.conj().T @ L @ B)).max()
        print(f'      reading: U_{m} (dim {B.shape[1]}, invariance residual {inv:.0e}), centre -{s_m}: '
              f'eigenvalues off Re = -s_m: {off or "none"}')


def sharp_gate(n, rational, root, couplings, J, cases):
    W = sum(bright(n))
    for name, (dims_expected, pairs_expected) in cases.items():
        dims, pairs = {}, {}
        for p in PRIMES:
            if name == 'full dipolar':
                Lrat, _ = lindbladian(n, rational, J, q=(0,))
                Lr, _ = lindbladian(n, {}, J, q=(0,), root=root, r3=1.0)
                Lz, _ = lindbladian(n, {}, J, q=(0,))
                r3p = next(x for x in range(p) if x * x % p == 3)
                Lp = (to_mod(Lrat, p) + r3p * to_mod(Lr - Lz, p)) % p
            else:
                cp, q, w = couplings[name]
                Lp = to_mod(lindbladian(n, cp, J, q=q, w=w)[0], p)
            dims[p], pairs[p] = [], []
            for m in range(n + 1):
                E = np.where(W == m)[0]
                O = np.where(W != m)[0]
                B = sharp_subspace(Lp, E, O, p)
                dims[p].append(B.shape[1])
                if B.shape[1]:
                    M = restricted(Lp, E, B, p)
                    pairs[p].append(palindromic_mod(M, 6 * m * J[0] + 6 * n * (J[1] + J[2]), p))
                else:
                    pairs[p].append(None)
        check(f'G8 n={n} J={J} {name}: dim K_m = {dims_expected} over both primes (upper bounds)',
              all(tuple(dims[p]) == dims_expected for p in PRIMES))
        check(f'G8 n={n} J={J} {name}: K_m palindromic about -s_m {pairs_expected} over both primes '
              f'(evidence; K_0 proved by G9)',
              all(tuple(pairs[p]) == pairs_expected for p in PRIMES))


# ---------------------------------------------------------------- G11: a shared gradient

CHANNELS = ('ZZ', 'ZX', 'ZY', 'XI', 'YI')


def kossakowski(ops, a):
    """sum_ij a_ij (A_i rho A_j^dag - 1/2 {A_j^dag A_i, rho}) and the price part sum_ij a_ij A_j^dag A_i"""
    dim = ops[0].shape[0]
    Id = np.eye(dim)
    D = np.zeros((dim * dim, dim * dim), complex)
    Ks = np.zeros((dim, dim), complex)
    for i, Ai in enumerate(ops):
        for j, Aj in enumerate(ops):
            if a[i][j] == 0:
                continue
            ji = Aj.conj().T @ Ai
            D += a[i][j] * (np.kron(Ai, Aj.conj()) - 0.5 * np.kron(ji, Id) - 0.5 * np.kron(Id, ji.T))
            Ks += a[i][j] * ji
    return D, Ks


def shared_lindbladian(J, cn, cd, H=None, imag=0):
    """Two nuclei under a shared gradient, everything times cd: per channel the Kossakowski matrix g·[[1, c], [c, 1]]
    on the two nuclei's copies of the string, c = cn/cd (cn one integer, or one per channel in the order ZZ, ZX, ZY,
    XI, YI). imag = e adds the time-odd cross-correlation of the J1 channel, the Hermitian entries ±i·e/cd·g between
    ZX on one nucleus and ZY on the other. Returns cd·L, cd·L_J0 and cd·K, K the price operator of the channels
    other than ZZ."""
    assert all(j % 2 == 0 for j in J), J
    cns = tuple(cn) if isinstance(cn, (tuple, list)) else (cn,) * 5
    dim = 16
    Id = np.eye(dim)
    L = np.zeros((dim * dim, dim * dim), complex)
    if H is not None:
        L += -1j * cd * (np.kron(H, Id) - np.kron(Id, H.T))
    LJ0 = np.zeros_like(L)
    K = np.zeros((dim, dim), complex)
    rates = dict(zip(CHANNELS, (3 * J[0] // 2, 3 * J[1] // 2, 3 * J[1] // 2, 3 * J[2] // 2, 3 * J[2] // 2)))
    for s_, c_ in zip(CHANNELS, cns):
        g = rates[s_]
        ops = (string(place(4, 0, s_)), string(place(4, 2, s_)))
        D, Ks = kossakowski(ops, [[g * cd, g * c_], [g * c_, g * cd]])
        L += D
        if s_ == 'ZZ':
            LJ0 += D
        else:
            K += Ks
    if imag:
        g = rates['ZX']
        ops = (string(place(4, 0, 'ZX')), string(place(4, 0, 'ZY')), string(place(4, 2, 'ZX')), string(place(4, 2, 'ZY')))
        a = [[0, 0, 0, 1j * imag * g], [0, 0, -1j * imag * g, 0], [0, 1j * imag * g, 0, 0], [-1j * imag * g, 0, 0, 0]]
        D, Ks = kossakowski(ops, a)
        L += D
        K += Ks
    return L, LJ0, K


def anticommutator_super(K):
    return np.kron(K, np.eye(K.shape[0])) + np.kron(np.eye(K.shape[0]), K.T)


def shared_gate():
    """G11: the price becomes an operator, a sign pattern is conserved, and one sector keeps the mirror"""
    sys.stdout.reconfigure(encoding='utf-8')
    perm, sign = s_perm(2)
    G = [string(s_ + s_) for s_ in CHANNELS[1:]]                       # Γ_a ⊗ Γ_a, a = 1..4 (real, symmetric)
    R = [np.kron(g, g) for g in G]                                       # X -> G X G
    q_w = 4 * 3 * (string('ZZII') + string('IIZZ')) + 5 * (4 * (string('ZIII') + string('IIZI'))
                                                          + 2 * (string('IZII') + string('IIIZ')))
    Iz1, Iz2 = string('ZIII') + 0.5 * string('IZII'), string('IIZI') + 0.5 * string('IIIZ')
    ising = 4 * Iz1 @ Iz2
    flipflop = (string('XXXX') + string('YYYY') + string('XXYY') + string('YYXX')
                + string('YXYX') - string('YXXY') - string('XYYX') + string('XYXY'))
    points = (((0, 4, 6), 1, 2), ((2, 4, 6), 1, 1), ((2, 4, 6), -1, 3), ((4, 2, 2), 1, 2), ((0, 4, 6), -1, 1))
    ok = True
    for J, cn, cd in points:
        for H in (None, q_w, q_w + 2 * ising + 3 * flipflop,
                  q_w + string('ZXZX') + string('ZYZY') + string('IXIX') + string('IYIY')):
            L, LJ0, K = shared_lindbladian(J, cn, cd, H)
            ok &= np.abs(conj_S(L, perm, sign) - (-L + 2 * LJ0 - anticommutator_super(K))).max() == 0
    Lb, LJ0b, Kb = shared_lindbladian((2, 4, 6), 1, 2, H=q_w + 2 * string('IXXX'))
    check('G11 shared gradient: S L S^-1 = -L + 2 L_J0 - {K, .} exactly, K = sum_(a>=1) g_a(2 + 2c Γ_a⊗Γ_a), at five '
          '(J, c), without a Hamiltonian, with offset and splitting, with an Ising coupling and the central-transition '
          'flip-flop added, and with the strings of the equal-splitting truncation (real, commuting with Q); a sqrt3 '
          'string (IXXX, anticommuting with Q) breaks it', ok
          and np.abs(conj_S(Lb, perm, sign) - (-Lb + 2 * LJ0b - anticommutator_super(Kb))).max() > 0)
    Lc0, _, K0 = shared_lindbladian((2, 4, 6), 0, 1)
    Lc1, _, K1 = shared_lindbladian((2, 4, 6), 1, 2)
    # times cd = 2: sum over the four channels a >= 1 of g_a·(2·cd + 2·cn·Γ_a⊗Γ_a), g = 6 (ZX, ZY) and 9 (XI, YI)
    closed = 4 * (6 + 6 + 9 + 9) * np.eye(16) + 2 * (6 * (string('ZXZX') + string('ZYZY'))
                                                      + 9 * (string('XIXI') + string('YIYI')))
    check('G11 the price is a number exactly without correlation: K = 2·sum g_a at c = 0; at c = 1/2 it is the closed form '
          'sum_(a>=1) g_a(2 + 2c Γ_a⊗Γ_a), not a multiple of the identity',
          np.abs(K0 - 2 * (3 * 4 + 3 * 6) * np.eye(16)).max() == 0
          and np.abs(K1 - closed).max() == 0 and np.abs(K1 - K1[0, 0] * np.eye(16)).max() > 0)
    check('G11 the four Γ_a⊗Γ_a (a >= 1) commute pairwise, and ZZ⊗ZZ is their product',
          all(np.abs(x @ y - y @ x).max() == 0 for x in G for y in G)
          and np.abs(G[0] @ G[1] @ G[2] @ G[3] - string('ZZZZ')).max() == 0)
    table = {(ch, k): tuple(np.abs(string(place(4, 2 * k, ch)) @ g + g @ string(place(4, 2 * k, ch))).max() == 0 for g in G)
             for ch in CHANNELS for k in (0, 1)}
    check('G11 the jumps against the four Γ_a⊗Γ_a: a jump of channel a >= 1 commutes with its own and anticommutes with '
          'the other three, on either nucleus; the J0 jump ZZ anticommutes with all four; so every jump changes the '
          'signs of ket and bra alike',
          all(table[(ch, k)] == tuple(b != a for b in range(4)) for a, ch in enumerate(CHANNELS[1:]) for k in (0, 1))
          and all(table[('ZZ', k)] == (True,) * 4 for k in (0, 1)))
    Pk = [np.diag(b.astype(float)) for b in bright(2)]
    Lm = shared_lindbladian((2, 4, 6), 1, 2)[0]
    L0m = shared_lindbladian((2, 4, 6), 0, 1)[0]
    Lc0only = shared_lindbladian((2, 4, 6), (1, 0, 0, 0, 0), 2, H=q_w + 2 * ising)[0]
    Lc1only = shared_lindbladian((2, 4, 6), (0, 1, 1, 0, 0), 2)[0]
    Lc2only = shared_lindbladian((2, 4, 6), (0, 0, 0, 1, 1), 2)[0]
    check('G11 the ladder mixes: with c1 alone, with c2 alone and with all on, the generator commutes with neither '
          'nucleus\'s bright projector (with both at c = 0, and with c0 alone, here with an Ising coupling)',
          all(np.abs(Lx @ P - P @ Lx).max() > 0 for Lx in (Lm, Lc1only, Lc2only) for P in Pk)
          and all(np.abs(L0m @ P - P @ L0m).max() == 0 for P in Pk)
          and all(np.abs(Lc0only @ P - P @ Lc0only).max() == 0 for P in Pk))
    W = Pk[0].diagonal() + Pk[1].diagonal()
    one, rest = np.where(W == 1)[0], np.where(W != 1)[0]
    check('G11 the correlated sandwiches change both nuclei\'s brightness at once, so the cells with one bright nucleus '
          'stay an invariant block at every correlation (they leak nowhere, nothing leaks in), while those with none and '
          'with two mix, with c1 alone, with c2 alone and with all on',
          all(np.abs(Lx[np.ix_(rest, one)]).max() == 0 and np.abs(Lx[np.ix_(one, rest)]).max() == 0
              for Lx in (Lm, Lc1only, Lc2only))
          and all(np.abs(Lx[np.ix_(np.where(W == 0)[0], np.where(W == 2)[0])]).max() > 0 for Lx in (Lm, Lc1only, Lc2only)))
    Iz2 = string('IIZI') + 0.5 * string('IIIZ')
    rot = True
    for (a, b), q in ((('ZX', 'ZY'), 1), (('XI', 'YI'), 2)):
        A, B = string(place(4, 2, a)), string(place(4, 2, b))
        rot &= (np.abs((Iz2 @ A - A @ Iz2) - 1j * q * B).max() == 0 and np.abs((Iz2 @ B - B @ Iz2) + 1j * q * A).max() == 0) \
            or (np.abs((Iz2 @ A - A @ Iz2) + 1j * q * B).max() == 0 and np.abs((Iz2 @ B - B @ Iz2) - 1j * q * A).max() == 0)
    zz2 = string(place(4, 2, 'ZZ'))
    ff_full = (string('XXXX') + string('YYYY') + string('XXYY') + string('YYXX')
               + string('YXYX') - string('YXXY') - string('XYYX') + string('XYXY'))
    check('G11 a rotation of nucleus 2 about the field turns its jumps of order q into each other by the angle q·phi '
          '([I_z2, ZX] and [I_z2, ZY] carry one into the other with factor 1, XI and YI with factor 2, ZZ untouched); the '
          'splitting and the Ising coupling commute with I_z2, the central-transition flip-flop does not',
          rot and np.abs(Iz2 @ zz2 - zz2 @ Iz2).max() == 0
          and np.abs(Iz2 @ q_w - q_w @ Iz2).max() == 0 and np.abs(Iz2 @ ising - ising @ Iz2).max() == 0
          and np.abs(Iz2 @ ff_full - ff_full @ Iz2).max() > 0)

    Lnone = shared_lindbladian((2, 4, 6), 0, 2, H=q_w + 2 * ising)[0]          # the same point at c = 0, times 2
    diff = Lc0only - Lnone
    blocks = {}
    for pattern in ((0, 0), (1, 0), (0, 1), (1, 1)):
        idx = np.where(np.all([Pk[k].diagonal() == b for k, b in enumerate(pattern)], axis=0))[0]
        blocks[pattern] = np.abs(diff[np.ix_(idx, idx)]).max()
    check('G11 with c0 alone the generator differs from the one at c = 0 only on the cells where both nuclei are bright: '
          'zero on the blocks with no bright nucleus and with one, so there the lemma\'s ladder holds exactly; non-zero on '
          'the block with two, there with both signs',
          blocks[(0, 0)] == 0 and blocks[(1, 0)] == 0 and blocks[(0, 1)] == 0 and blocks[(1, 1)] > 0
          and np.abs(diff - np.diag(np.diag(diff))).max() == 0
          and np.diag(diff).real.min() < 0 < np.diag(diff).real.max())
    cons = all(np.abs(shared_lindbladian(J, cn, cd)[0] @ r - r @ shared_lindbladian(J, cn, cd)[0]).max() == 0
               for J, cn, cd in points for r in R)
    check('G11 without a Hamiltonian the pattern t = s_ket ⊙ s_bra is conserved: L commutes with X -> Γ_a⊗Γ_a X '
          'Γ_a⊗Γ_a for a = 1..4, at every point, local and correlated parts alike', cons)
    Iz1, Iz2 = string('ZIII') + 0.5 * string('IZII'), string('IIZI') + 0.5 * string('IIIZ')
    terms = {'offset': Iz1 + Iz2, 'splitting': string('ZZII') + string('IIZZ'), 'Ising Iz1·Iz2': Iz1 @ Iz2,
             'central-transition flip-flop': string('XXXX') + string('YYYY') + string('XXYY') + string('YYXX')
             + string('YXYX') - string('YXXY') - string('XYYX') + string('XYXY')}
    broken = {}
    for name, Hx in terms.items():
        LH = -1j * (np.kron(Hx, np.eye(16)) - np.kron(np.eye(16), Hx.T))
        broken[name] = any(np.abs(LH @ r - r @ LH).max() > 0 for r in R)
    check(f'G11 every Hamiltonian term breaks the conservation: {", ".join(broken)}', all(broken.values()))
    sat_ff = (string('IXIX') + string('IYIY')) @ (np.eye(16) - string('ZZZZ')) / 2
    LH = -1j * (np.kron(sat_ff, np.eye(16)) - np.kron(np.eye(16), sat_ff.T))
    check('G11 the satellite flip-flop of the equal-splitting truncation, (IXIX + IYIY)(1 - Q)/2, lies in the algebra of '
          'the Γ_a⊗Γ_a (IYIY = -ZXZX·Q) and keeps the conservation',
          np.abs(string('IYIY') + G[0] @ string('ZZZZ')).max() == 0 and all(np.abs(LH @ r - r @ LH).max() == 0 for r in R))

    # a separate correlation per channel (anisotropic motion) keeps all of it; a time-odd cross-correlation does not
    cv = (2, 3, 3, -2, -2)
    Lv, LJ0v, Kv = shared_lindbladian((0, 4, 6), cv, 6)
    check('G11 one correlation per channel, c = (2, 3, 3, -2, -2)/6: the identity holds with K = sum g_a(2 + 2c_a Γ_a⊗Γ_a) '
          'and the pattern t is conserved',
          np.abs(conj_S(Lv, perm, sign) - (-Lv + 2 * LJ0v - anticommutator_super(Kv))).max() == 0
          and all(np.abs(Lv @ r - r @ Lv).max() == 0 for r in R))
    Li, _, _ = shared_lindbladian((0, 4, 6), 1, 2, imag=1)
    check('G11 a time-odd cross-correlation (±i between ZX on one nucleus and ZY on the other) breaks the conservation of t',
          any(np.abs(Li @ r - r @ Li).max() > 0 for r in R))
    Hk = 2 * G[0] + 3 * G[2] - G[3] + 5 * string('ZZZZ')
    Lk, LJ0k, Kk = shared_lindbladian((0, 4, 6), 1, 2, H=Hk)
    check('G11 a real Hamiltonian built from the Γ_a⊗Γ_a and Q keeps the identity and the conservation of t',
          np.abs(conj_S(Lk, perm, sign) - (-Lk + 2 * LJ0k - anticommutator_super(Kk))).max() == 0
          and all(np.abs(Lk @ r - r @ Lk).max() == 0 for r in R))

    Iz_tot = Iz1 + Iz2
    adIz = np.kron(Iz_tot, np.eye(16)) - np.kron(np.eye(16), Iz_tot.T)
    comm_off = lambda L_: np.abs(L_ @ adIz - adIz @ L_).max() == 0
    check('G11 the common offset commutes with the dissipator when the correlation is one per coherence order q '
          '(c_ZX = c_ZY, c_XI = c_YI; the secular case), at J0 = 2 and J0 = 0, at a common c and at c = (2, 3, 3, -2, '
          '-2)/6, and not when c_ZX != c_ZY: it turns the modes and moves no rate',
          all(comm_off(shared_lindbladian(J, cn, cd)[0]) for J in ((2, 4, 6), (0, 4, 6)) for cn, cd in ((1, 2), (cv, 6)))
          and not comm_off(shared_lindbladian((2, 4, 6), (2, 3, 1, -2, -2), 6)[0]))

    def projector(t):
        P = np.eye(256)
        for ta, r in zip(t, R):
            P = P @ (np.eye(256) + ta * r) / 2
        return P
    Pm = projector((-1, -1, -1, -1))
    check('G11 S keeps each of the sixteen sectors t',
          all(np.abs(conj_S(projector(t), perm, sign) - projector(t)).max() == 0 for t in product((1, -1), repeat=4)))
    sand = True
    for ch in CHANNELS[1:]:
        A1, A2 = string(place(4, 0, ch)), string(place(4, 2, ch))
        local = np.kron(A1, A1.conj()) + np.kron(A2, A2.conj())
        cross = np.kron(A1, A2.conj()) + np.kron(A2, A1.conj())
        sand &= np.abs(local @ Pm).max() == 0 and np.abs(cross @ Pm).max() == 0
    check('G11 on the sector (-,-,-,-) every sandwich of the channels a >= 1 vanishes, local (A_k rho A_k) and '
          'correlated (A_1 rho A_2 + A_2 rho A_1) alike: the jumps carry these modes nowhere', sand)
    scal = all(np.abs(shared_lindbladian(J, cn, cd)[0] @ Pm + 6 * (J[1] + J[2]) * cd * Pm).max() == 0
               for J, cn, cd in (((0, 4, 6), 1, 2), ((0, 4, 6), 1, 1), ((0, 4, 6), -1, 2), ((0, 2, 8), 1, 3),
                                 ((0, 4, 6), cv, 6), ((0, 4, 6), 0, 1)))
    Lsh = shared_lindbladian((2, 4, 6), 1, 2)[0]
    check('G11 at J0 = 0 without a Hamiltonian the sector (-,-,-,-) is an eigenspace of L at -6(J1 + J2), at every '
          'correlation tried (c = 0, 1/2, 1, -1/2, 1/3, one per q) and two rate points; not at J0 = 2, c = 1/2',
          scal and np.abs(Lsh @ Pm + 6 * 10 * 2 * Pm).max() > 0)
    LHk = -1j * (np.kron(Hk, np.eye(16)) - np.kron(np.eye(16), Hk.T))
    Lwith = shared_lindbladian((0, 4, 6), 1, 2, H=Hk)[0]
    check('G11 with the real Hamiltonian built from the Γ_a⊗Γ_a and Q the sector (-,-,-,-) stays invariant and L there is '
          '-6(J1 + J2) plus the anti-Hermitian Hamiltonian part: the modes gain frequencies and keep the rate',
          np.abs(Lwith @ Pm - Pm @ Lwith @ Pm).max() == 0
          and np.abs(Lwith @ Pm - (-6 * 10 * 2 * Pm + 2 * LHk @ Pm)).max() == 0)
    okp = np.abs(conj_S(Pm, perm, sign) - Pm).max() == 0
    for J, cn, cd in points:
        _, _, K = shared_lindbladian(J, cn, cd)
        price = 2 * cd * 2 * (3 * J[1] + 3 * J[2])          # K(s) + K(-s) = 4·sum_(a>=1) g_a, times cd
        okp &= np.abs((anticommutator_super(K) - price * np.eye(256)) @ Pm).max() == 0
    check('G11 on the sector t = (-,-,-,-), ket and bra opposite in every sign, the price is the number 4·sum_(a>=1) '
          'g_a at every c, and S keeps the sector; so at J0 = 0 without a Hamiltonian it pairs about -6(J1 + J2), '
          'the lemma\'s global centre at J0 = 0, where every rung of the ladder sits', okp)

    for p in PRIMES:
        i16 = pow(16, p - 2, p)
        res, dims = {}, set()
        for J, cn, cd, Hx in (((0, 4, 6), 1, 2, None), ((0, 4, 6), 1, 1, None), ((2, 4, 6), 1, 2, None),
                              ((0, 4, 6), (2, 3, 3, -2, -2), 6, None), ((0, 4, 6), 1, 2, 'Hk'),
                              ((2, 4, 6), 0, 1, None)):
            L, _, _ = shared_lindbladian(J, cn, cd, H=Hk if Hx else None)
            key = (J, cn, cd, Hx)
            Lp = to_mod(L, p)
            row = []
            for t in product((1, -1), repeat=4):
                P16f = 16 * projector(t)
                assert np.abs(P16f.imag).max() == 0 and np.abs(P16f.real - np.rint(P16f.real)).max() == 0
                P16 = np.rint(P16f.real).astype(np.int64)
                B = null_mod((16 * np.eye(256, dtype=np.int64) - P16) % p, p, 256)
                dims.add(B.shape[1])
                M = restricted(Lp, np.arange(256), B, p)
                c2 = (-2 * int(np.trace(M) % p) * i16) % p                 # symmetric about the sector's mean
                row.append(palindromic_mod(M, c2, p))
            res[key] = row
        all_minus = 15                                                          # index of (-1, -1, -1, -1)
        check(f'G11 over GF({p}): each of the sixteen sectors t has dimension 16', dims == {16})
        check(f'G11 over GF({p}) at J0 = 0, c = 1/2 and c = 1: the sector (-,-,-,-) palindromic about its mean '
              f'(evidence; the algebra above proves it), each of the other fifteen not about its mean (exact)',
              all(res[k][all_minus] and not any(res[k][:all_minus])
                  for k in (((0, 4, 6), 1, 2, None), ((0, 4, 6), 1, 1, None))))
        check(f'G11 over GF({p}) at J0 = 0: the sector (-,-,-,-) palindromic about its mean also with one correlation '
              f'per channel and with the real Hamiltonian built from the Γ_a⊗Γ_a and Q (evidence)',
              res[((0, 4, 6), (2, 3, 3, -2, -2), 6, None)][all_minus] and res[((0, 4, 6), 1, 2, 'Hk')][all_minus])
        check(f'G11 over GF({p}) at J0 = 2, c = 1/2: the sector (-,-,-,-) is not palindromic about its mean (exact): '
              f'the correlated J0 jump keeps no number on it; at c = 0 it is (the control, evidence)',
              not res[((2, 4, 6), 1, 2, None)][all_minus] and res[((2, 4, 6), 0, 1, None)][all_minus])
        Lone = to_mod(shared_lindbladian((2, 4, 6), 1, 2)[0], p)
        Lone0 = to_mod(shared_lindbladian((0, 4, 6), 1, 2)[0], p)
        mone = []
        trunc1 = q_w + string('ZXZX') + string('ZYZY') + string('IXIX') + string('IYIY') + 2 * ising
        for Lx in (Lone, Lone0, to_mod(shared_lindbladian((2, 4, 6), (0, 1, 1, 0, 0), 2)[0], p),
                   to_mod(shared_lindbladian((2, 4, 6), (0, 0, 0, 1, 1), 2)[0], p),
                   to_mod(shared_lindbladian((2, 4, 6), 1, 2, H=trunc1)[0], p)):
            Mo = Lx[np.ix_(one, one)]
            co = (-2 * int(np.trace(Mo) % p) * pow(len(one), p - 2, p)) % p
            mone.append(not palindromic_mod(Mo, co, p))
        check(f'G11 over GF({p}): the block with one bright nucleus is not palindromic about its mean at c = 1/2 in every '
              f'channel (J0 = 2 and J0 = 0), with c1 = 1/2 alone and c2 = 1/2 alone (J0 = 2), and with the equal-splitting '
              f'truncation, an Ising coupling, offset and splitting (exact): it stays separate and loses its centre', all(mone))
        trunc = q_w + string('ZXZX') + string('ZYZY') + string('IXIX') + string('IYIY') + 2 * ising
        nocentre = []
        for cn, Hx in (((0, 1, 1, 0, 0), None), ((0, 0, 0, 1, 1), None), (1, trunc)):
            Lx = to_mod(shared_lindbladian((0, 4, 6), cn, 2, H=Hx)[0], p)
            cx = (-2 * int(np.trace(Lx) % p) * pow(256, p - 2, p)) % p
            nocentre.append(not palindromic_mod(Lx, cx, p))
        check(f'G11 over GF({p}) at J0 = 0: no centre with c1 = 1/2 alone, with c2 = 1/2 alone, and at c = 1/2 with the '
              f'equal-splitting truncation, an Ising coupling, offset and splitting (exact)', all(nocentre))
        Lg, _, _ = shared_lindbladian((0, 4, 6), 1, 2)
        Lg0, _, _ = shared_lindbladian((0, 4, 6), 0, 1)
        cg = (-2 * int(np.trace(to_mod(Lg, p)) % p) * pow(256, p - 2, p)) % p
        cg0 = (-2 * int(np.trace(to_mod(Lg0, p)) % p) * pow(256, p - 2, p)) % p
        Lc0 = to_mod(shared_lindbladian((2, 4, 6), (1, 0, 0, 0, 0), 2, H=q_w + 2 * ising)[0], p)
        Pb = bright(2)
        rungs = []
        for pattern in ((0, 0), (1, 0), (0, 1), (1, 1)):
            idx = np.where(np.all([Pb[k] == b for k, b in enumerate(pattern)], axis=0))[0]
            Mb = Lc0[np.ix_(idx, idx)]
            c2b = (-2 * int(np.trace(Mb) % p) * pow(len(idx), p - 2, p)) % p
            rungs.append(palindromic_mod(Mb, c2b, p))
        check(f'G11 over GF({p}) with c0 alone (c0 = 1/2, J = (2, 4, 6), offset, splitting and an Ising coupling): the '
              f'blocks with no bright nucleus and with one are palindromic about their means (evidence), the block with two '
              f'is not (exact): the correlated ZZ term charges only the two-bright cells, with both signs',
              rungs == [True, True, True, False])
        check(f'G11 over GF({p}) at J0 = 0: the whole spectrum is not palindromic about its mean at c = 1/2 (exact), '
              f'and is at c = 0 (the control, proved by the lemma)',
              not palindromic_mod(to_mod(Lg, p), cg, p) and palindromic_mod(to_mod(Lg0, p), cg0, p))


if __name__ == '__main__':
    sys.exit(main())
