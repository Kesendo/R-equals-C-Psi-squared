"""Where the two ends differ, a far element splits off what pairs: exact gate over Q(i).

The question it answers. F158 (docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md) writes the two ends of
L(rho) = -i[H, rho] + sum_l g_l (A_l rho A_l - rho) as N = {X : [H, X] = 0, A_l X A_l = X} and
W = {W : [H, W] = 0, A_l W A_l = -W for every l}. Where W holds an invertible element the spectrum pairs about
-sigma. What does a far element do on a row with 0 < dim W < dim N, where every far element is singular?

  (1) The range R and the corange C (the range of W^dagger) of a far element are invariant under H and under
      every jump.
  (2) For every operator rho, L(rho W) = [(-L^dagger - 2 sigma)(rho)] W: right multiplication by W carries
      -L^dagger - 2 sigma on the operators whose bra lies in R onto L on the operators whose bra lies in C,
  (3) so the spectra of L on those two pieces are mirror images under the LINE reflection
      lambda -> -conj(lambda) - 2 sigma, with multiplicity: p_C(x) = (-1)^n conj(p_R)(-x - 2 sigma).
  (4) The split. For a Hermitian far element of maximal rank, H = R + R^perp splits the system: the system on R
      pairs about -sigma (point reflection), the coherences between R and R^perp in both orders pair about -sigma,
      and the system on R^perp has no far element and carries the whole deficit dim N - dim W. The three parts
      make up the spectrum of L.

Everything is exact: Pauli matrices over the Gaussian rationals, rational rates, far and near bases from
framework.diagnostics.f158_end_count (exact kernels over Q(i)), characteristic polynomials by sympy DomainMatrix,
nullities by exact ranks. No tolerance appears anywhere.

What would make it fail, each fed through the same door as the real input: W + 1, which is not a far element,
breaks (2); the two pieces of a singular NEAR element do not reflect; the piece on the complement of the range
does not reflect onto itself; with several jumps at unequal rates, the reflection about the largest single rate
instead of sigma fails; the point reflection fails on pieces where only the line reflection is implied; and a
Hermitian far element of less than maximal rank leaves far elements in the system on its complement.

Rows: the Heisenberg triangle d, u, v with the jump Z on d alone (the cat cell |000><111| at fields along z
summing to zero, a resonant two-level block, the singlet class S with a planar field on d, Z letters -1 and +1 on
u and v at couplings (-1, 1, -1)), an XY chain of three sites with every site dephased at unequal rates and z
fields summing to zero (the cat of F22), and two qubits built as two blocks (F158's gate-11 construction), with
one jump and with two at unequal rates.
"""
from __future__ import annotations

import itertools
import sys
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ_I
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from framework.diagnostics.f158_end_count import end_count, _kernel, _span  # noqa: E402

OUT = ROOT / "results" / "far_element_mirror_gate.txt"
LINES = []
FAILS = []


def say(msg=""):
    print(msg)
    LINES.append(msg)


def check(cond, what):
    say(f"   {'PASS' if cond else 'FAIL'}  {what}")
    if not cond:
        FAILS.append(what)


I = sp.I
PAULI = {'I': sp.Matrix([[1, 0], [0, 1]]), 'X': sp.Matrix([[0, 1], [1, 0]]),
         'Y': sp.Matrix([[0, -I], [I, 0]]), 'Z': sp.Matrix([[1, 0], [0, -1]])}
ZERO = None


def mat(s):
    M = sp.Matrix([[1]])
    for L in s:
        M = sp.kronecker_product(M, PAULI[L])
    return M


def op(terms, n):
    d = 2 ** n
    M = sp.zeros(d, d)
    for s, c in terms.items():
        M += sp.Rational(c) * mat(s)
    return M.applyfunc(sp.expand)


def heis_triangle(J, fd, fu, fv):
    """Sites (d, u, v) = (0, 1, 2); Heisenberg bonds J = (J_du, J_dv, J_uv); fields as vectors."""
    h = {}
    for (i, j, c) in ((0, 1, J[0]), (0, 2, J[1]), (1, 2, J[2])):
        if c:
            for L in 'XYZ':
                s = ['I'] * 3; s[i] = s[j] = L
                h[''.join(s)] = h.get(''.join(s), 0) + c
    for site, f in ((0, fd), (1, fu), (2, fv)):
        for L, a in zip('XYZ', f):
            if a:
                s = ['I'] * 3; s[site] = L
                h[''.join(s)] = h.get(''.join(s), 0) + a
    return h


def pauli_terms(M, n):
    """Exact Pauli-string coefficients of a Hermitian matrix with rational entries: c_s = Tr(P_s M) / 2^n."""
    out = {}
    for letters in itertools.product('IXYZ', repeat=n):
        s = ''.join(letters)
        c = sp.expand((mat(s) * M).trace() / 2 ** n)
        if c != 0:
            out[s] = sp.Rational(c)
    return out


def four_blocks():
    """Three qubits, blocks by the bits of qubits 1 and 2: X_0, X_0 + 3, Z_0 + 7, Z_0 + 11 on qubit 0."""
    P = {b: sp.diag(*[1 if k == b else 0 for k in range(4)]) for b in range(4)}
    h = [PAULI['X'], PAULI['X'] + 3 * PAULI['I'], PAULI['Z'] + 7 * PAULI['I'], PAULI['Z'] + 11 * PAULI['I']]
    M = sum((sp.kronecker_product(h[b], P[b]) for b in range(4)), sp.zeros(8, 8))
    return pauli_terms(M, 3)


def to_matrix(v, n):
    """The operator of a basis vector {string: exact Gaussian rational}, scaled so its first coefficient is 1."""
    d = 2 ** n
    lead = sp.sympify(next(iter(v.values())))
    M = sp.zeros(d, d)
    for s, c in v.items():
        M += sp.sympify(c) / lead * mat(s)
    return M.applyfunc(sp.expand)


def superop(H, jumps, rates):
    """Row-major vec: vec(A rho B) = (A kron B^T) vec(rho)."""
    d = H.shape[0]
    Id = sp.eye(d)
    L = -I * (sp.kronecker_product(H, Id) - sp.kronecker_product(Id, H.T))
    for A, g in zip(jumps, rates):
        L += g * (sp.kronecker_product(A, A.T) - sp.kronecker_product(Id, Id))
    return L.applyfunc(sp.expand)


def apply_L(H, jumps, rates, rho):
    out = -I * (H * rho - rho * H)
    for A, g in zip(jumps, rates):
        out += g * (A * rho * A - rho)
    return out


def apply_reflected(H, jumps, rates, rho):
    """(-L^dagger - 2 sigma)(rho) = -i[H, rho] - sum g (A rho A - rho) - 2 sigma rho."""
    sigma = sum(rates)
    out = -I * (H * rho - rho * H)
    for A, g in zip(jumps, rates):
        out -= g * (A * rho * A - rho)
    return out - 2 * sigma * rho


def is_zero(M):
    return all(sp.expand(e) == 0 for e in M)


def intertwines(H, jumps, rates, W):
    d = H.shape[0]
    for i, j in itertools.product(range(d), repeat=2):
        E = sp.zeros(d, d); E[i, j] = 1
        if not is_zero(apply_L(H, jumps, rates, E * W) - apply_reflected(H, jumps, rates, E) * W):
            return False
    return True


def column_basis(M):
    cols = M.columnspace()
    return sp.Matrix.hstack(*cols) if cols else sp.zeros(M.shape[0], 0)


def invariant(Bm, ops):
    r = Bm.rank()
    return all(sp.Matrix.hstack(Bm, (O * Bm).applyfunc(sp.expand)).rank() == r for O in ops)


def block_basis(K, B):
    """vec basis of span{|k><b| : k a column of K, b a column of B}: kets in span(K), bras in span(B)."""
    d = K.shape[0]
    return sp.Matrix.hstack(*[(K[:, a] * B[:, b].H).applyfunc(sp.expand).reshape(d * d, 1)
                              for a in range(K.shape[1]) for b in range(B.shape[1])])


def restricted(Lsup, P):
    """The matrix M with Lsup*P == P*M exactly (raises if the span of P is not invariant)."""
    G = P.H * P
    M = (G.inv() * (P.H * (Lsup * P))).applyfunc(sp.expand)
    assert is_zero(Lsup * P - P * M), "piece not invariant"
    return M


def nullity(M):
    return M.shape[0] - DomainMatrix.from_Matrix(M.applyfunc(sp.expand)).convert_to(QQ_I).rank()


x = sp.Symbol('x')


def charpoly(M):
    n = M.shape[0]
    dm = DomainMatrix.from_Matrix(M.applyfunc(sp.expand)).convert_to(QQ_I)
    coeffs = [QQ_I.to_sympy(c) for c in dm.charpoly()]
    return sp.Poly(sum(c * x ** (n - k) for k, c in enumerate(coeffs)), x)


def line_reflect(p, sigma):
    """(-1)^n conj(p)(-x - 2 sigma): the char poly of -M^dagger - 2 sigma when p is that of M."""
    n = p.degree()
    pc = sp.Poly([sp.conjugate(c) for c in p.all_coeffs()], x)
    return sp.Poly(sp.expand((-1) ** n * pc.as_expr().subs(x, -x - 2 * sigma)), x)


def point_reflect(p, sigma):
    n = p.degree()
    return sp.Poly(sp.expand((-1) ** n * p.as_expr().subs(x, -x - 2 * sigma)), x)


def same(p, q):
    return sp.expand(p.as_expr() - q.as_expr()) == 0


def product(*ps):
    return sp.Poly(sp.expand(sp.Mul(*[p.as_expr() for p in ps])), x)


# ---------------------------------------------------------------------------------------------------------------
ROWS = [
    ("triangle, every field along z summing to zero (the cat cell)", 3,
     heis_triangle((2, 2, 1), (0, 0, -1), (0, 0, -1), (0, 0, 2)), ['ZII'], [sp.Rational(3, 10)]),
    ("triangle, swap-symmetric, fields along z not summing to zero (a resonant two-level block)", 3,
     heis_triangle((2, 2, 1), (0, 0, 1), (0, 0, -1), (0, 0, -1)), ['ZII'], [sp.Rational(3, 10)]),
    ("triangle, the singlet class S with a planar field on d", 3,
     heis_triangle((2, 2, 3), (1, 0, 0), (1, 2, 1), (1, 2, 1)), ['ZII'], [sp.Rational(3, 10)]),
    ("triangle, Z letters -1 and +1 on u and v, couplings (-1, 1, -1), no field on d", 3,
     heis_triangle((-1, 1, -1), (0, 0, 0), (0, 0, -1), (0, 0, 1)), ['ZII'], [sp.Rational(3, 10)]),
    ("XY chain of three, every site dephased at unequal rates, z fields summing to zero (F22's cat)", 3,
     {'XXI': 1, 'YYI': 1, 'IXX': 1, 'IYY': 1, 'ZII': 1, 'IZI': 2, 'IIZ': -3}, ['ZII', 'IZI', 'IIZ'],
     [sp.Rational(1, 10), sp.Rational(2, 10), sp.Rational(3, 10)]),
    ("three qubits built as four blocks, two of them carrying an invertible far element (gate 11 widened)", 3,
     four_blocks(), ['ZII'], [sp.Rational(2, 5)]),
    ("two qubits built as two blocks, H = (1 + Z_1) X_0 + (1 - Z_1)(Z_0 + 5), jump Z_0 (F158 gate 11)", 2,
     {'XI': 1, 'XZ': 1, 'ZI': 1, 'ZZ': -1, 'II': 5, 'IZ': -5}, ['ZI'], [sp.Rational(1, 2)]),
    ("the same blocks with a second jump Z_0 Z_1 at another rate", 2,
     {'XI': 1, 'XZ': 1, 'ZI': 1, 'ZZ': -1, 'II': 5, 'IZ': -5}, ['ZI', 'ZZ'], [sp.Rational(1, 2), sp.Rational(1, 5)]),
]

say("Where the two ends differ, a far element splits off what pairs (exact over Q(i))")
say()
tally = dict(read=0, distinct=0, point_fails=0, book_fails=0, mutation=0, complement=0, near=0, nonmax=0)
for (label, n, hterms, jstr, rates) in ROWS:
    say(label)
    d = 2 ** n
    sigma = sum(rates)
    H = op(hterms, n)
    jumps = [mat(s) for s in jstr]
    hq = {s: sp.Rational(c) for s, c in hterms.items()}
    ec = end_count(n, hq, jstr)
    say(f"   jumps {jstr} at rates {[str(g) for g in rates]}, sigma = {sigma}; (near, far) = ({ec.near}, {ec.far})")
    check(0 < ec.far < ec.near, "a strict row: 0 < dim W < dim N, so no far element is invertible (F158 Lemma 3)")
    Lsup = superop(H, jumps, rates)
    pL = charpoly(Lsup)
    check(not same(pL, point_reflect(pL, sigma)), "the whole spectrum does not pair about -sigma (exact char poly)")
    Ws = [to_matrix(v, n) for v in ec.far_basis]
    if len(Ws) > 1:
        Ws.append(sum(((k + 2) * W for k, W in enumerate(Ws)), sp.zeros(d, d)).applyfunc(sp.expand))
        Ws.append((Ws[0] + I * Ws[1]).applyfunc(sp.expand))
    nb = len(ec.far_basis)
    for k, W in enumerate(Ws):
        tag = (f"far element {k + 1}" if k < nb else
               "a combination of the far basis" if k == nb else "a complex combination of the far basis")
        tally['read'] += 1
        r = W.rank()
        R, C = column_basis(W), column_basis(W.H)
        same_rc = sp.Matrix.hstack(R, C).rank() == r
        tally['distinct'] += (not same_rc)
        say(f"   {tag}: rank {r} of {d}, range {'=' if same_rc else '!='} corange")
        check(is_zero(H * W - W * H) and all(is_zero(A * W * A + W) for A in jumps) and r < d,
              f"{tag} lies in W and is singular")
        check(invariant(R, [H] + jumps) and invariant(C, [H] + jumps),
              f"{tag}: range and corange invariant under H and every jump (1)")
        check(intertwines(H, jumps, rates, W),
              f"{tag}: L(rho W) = [(-L^dag - 2 sigma) rho] W on all {d * d} matrix units (2)")
        PR, PC = block_basis(sp.eye(d), R), block_basis(sp.eye(d), C)
        pR, pC = charpoly(restricted(Lsup, PR)), charpoly(restricted(Lsup, PC))
        check(same(pC, line_reflect(pR, sigma)),
              f"{tag}: pieces of dimension {PR.shape[1]}: p_C(x) = (-1)^n conj(p_R)(-x - 2 sigma) exactly (3)")
        if len(rates) > 1:
            wrong = not same(pC, line_reflect(pR, max(rates)))
            tally['book_fails'] += wrong
            say(f"   {tag}: about the largest single rate instead of sigma the reflection "
                f"{'fails' if wrong else 'holds'} (read)")
        pt = same(pC, point_reflect(pR, sigma))
        tally['point_fails'] += (not pt)
        say(f"   {tag}: the point reflection p_C(x) = (-1)^n p_R(-x - 2 sigma) {'holds' if pt else 'fails'} (read)")
        if k == 0:
            fired = not intertwines(H, jumps, rates, W + sp.eye(d))
            tally['mutation'] += fired
            say(f"   mutation (W + 1, not in W): the identity (2) {'fails' if fired else 'still holds'} (read)")
            perp = R.H.nullspace()
            if perp:
                Rp = sp.Matrix.hstack(*perp)
                pP = charpoly(restricted(Lsup, block_basis(sp.eye(d), Rp)))
                cf = not same(pP, line_reflect(pP, sigma))
                tally['complement'] += cf
                say(f"   control: the piece on the complement of the range reflects onto itself: "
                    f"{'no' if cf else 'yes'} (read)")
    near = [to_matrix(v, n) for v in _kernel(hq, _span(n, jstr, lit=False))]
    sing = [Nm for Nm in near if Nm.rank() < d]
    if sing:
        Nm = sing[0]
        pRN = charpoly(restricted(Lsup, block_basis(sp.eye(d), column_basis(Nm))))
        pCN = charpoly(restricted(Lsup, block_basis(sp.eye(d), column_basis(Nm.H))))
        nf = not same(pCN, line_reflect(pRN, sigma))
        tally['near'] += nf
        say(f"   control: a singular near element of rank {Nm.rank()}: its two pieces are line reflections: "
            f"{'no' if nf else 'yes'} (read)")
    # (4) the split by a Hermitian far element of maximal rank among those read
    # candidates: Hermitian parts of the far elements read, and of near times far (N W lies in W: F158's grading)
    projs = []
    for Nm in near:                       # spectral projectors of Hermitian near elements lie in N (polynomials in it)
        Nh = (Nm + Nm.H).applyfunc(sp.expand)
        for lam, mult, vecs in Nh.eigenvects():
            if lam.is_rational and 0 < len(vecs) < d:
                V = sp.Matrix.hstack(*vecs)
                projs.append((V * (V.H * V).inv() * V.H).applyfunc(sp.expand))
    pool = Ws + [(Pm * W).applyfunc(sp.expand) for Pm in projs for W in Ws[:nb]]
    herm = [(W + W.H).applyfunc(sp.expand) for W in pool] + [(I * (W - W.H)).applyfunc(sp.expand) for W in pool]
    herm = [h for h in herm if not is_zero(h)]
    herm.append(sum(((k + 2) * h for k, h in enumerate(herm)), sp.zeros(d, d)).applyfunc(sp.expand))
    ranks = [h.rank() for h in herm]
    Wh = herm[ranks.index(max(ranks))]
    r = Wh.rank()
    R = column_basis(Wh)
    Rp = sp.Matrix.hstack(*R.H.nullspace())
    say(f"   the split by a Hermitian far element of rank {r}: dim R = {r}, dim R^perp = {d - r}")
    check(is_zero(Wh - Wh.H) and is_zero(H * Wh - Wh * H) and all(is_zero(A * Wh * A + Wh) for A in jumps),
          "the splitting element is Hermitian and lies in W")
    check(invariant(R, [H] + jumps) and invariant(Rp, [H] + jumps), "R and R^perp invariant under H and every jump")
    pRR = charpoly(restricted(Lsup, block_basis(R, R)))
    check(same(pRR, point_reflect(pRR, sigma)), f"the system on R ({r * r} modes) pairs about -sigma")
    pX1 = charpoly(restricted(Lsup, block_basis(Rp, R)))
    pX2 = charpoly(restricted(Lsup, block_basis(R, Rp)))
    check(same(pX1, line_reflect(pX1, sigma)), "the coherences with ket in R^perp and bra in R: closed under the line "
          "reflection")
    pXX = product(pX1, pX2)
    check(same(pXX, point_reflect(pXX, sigma)),
          f"the coherences between R and R^perp in both orders ({2 * r * (d - r)} modes) pair about -sigma")
    Lrem = restricted(Lsup, block_basis(Rp, Rp))
    m = Lrem.shape[0]
    dn_rem, dw_rem = nullity(Lrem), nullity(Lrem + 2 * sigma * sp.eye(m))
    check(dn_rem - dw_rem == ec.near - ec.far and dw_rem == 0,
          f"the system on R^perp: (near, far) = ({dn_rem}, {dw_rem}), the whole deficit "
          f"{ec.near - ec.far} and no far element")
    check(same(product(pRR, pXX, charpoly(Lrem)), pL),
          f"the three parts make up the spectrum of L: {d * d - (d - r) ** 2} of its {d * d} eigenvalues in the parts "
          "that pair")
    smaller = [h for h, rk in zip(herm, ranks) if 0 < rk < r]
    if smaller:
        Rs = column_basis(smaller[0])
        Rsp = sp.Matrix.hstack(*Rs.H.nullspace())
        Ls = restricted(Lsup, block_basis(Rsp, Rsp))
        fs = nullity(Ls + 2 * sigma * sp.eye(Ls.shape[0]))
        tally['nonmax'] += (fs > 0)
        say(f"   control: a Hermitian far element of rank {smaller[0].rank()} leaves {fs} far elements in the system "
            "on its complement (read)")
    say()

say("controls across the rows")
check(tally['read'] >= 10, f"far elements read: {tally['read']}")
check(tally['distinct'] >= 1, f"far elements whose range and corange differ: {tally['distinct']}")
check(tally['mutation'] == len(ROWS), f"the mutation W + 1 breaks (2) on {tally['mutation']} of {len(ROWS)} rows")
check(tally['near'] >= 1, f"the pieces of a singular near element fail to reflect on {tally['near']} rows: "
      "(3) needs a FAR element")
check(tally['complement'] >= 1, f"the piece on the complement of the range fails to reflect on {tally['complement']} "
      "rows: (3) is about the piece, not the whole")
check(tally['book_fails'] >= 1, f"with several jumps at unequal rates the reflection about the largest rate fails on "
      f"{tally['book_fails']} far elements: the centre is sigma, the sum of the rates")
check(tally['point_fails'] >= 1, f"the point reflection fails on {tally['point_fails']} pieces: (3) is the line "
      "reflection")
check(tally['nonmax'] >= 1, f"a Hermitian far element of less than maximal rank leaves far elements in its complement "
      f"on {tally['nonmax']} rows: (4) needs maximal rank")
say()
say("ALL CHECKS PASS" if not FAILS else f"{len(FAILS)} CHECK(S) FAILED")
OUT.write_text("\n".join(LINES) + "\n", encoding="utf-8")
sys.exit(1 if FAILS else 0)
