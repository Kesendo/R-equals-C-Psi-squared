# r90_frozen_divisor_gate.py
#
# Gate for docs/proofs/PROOF_R90_FROZEN_DIVISOR.md: the R90-locus frozen
# divisor (J-independent corner-block eigenvalues at lambda = -4*gbar with
# multiplicity floor(N/2), proved by the tauQ pencil argument).
#
# Checks (all must pass; prints "R90 frozen divisor gate: ALL GREEN"):
#   G0  block builder cross-check against the framework Liouvillian (N=4)
#   G1  the mirror identity tauQ (M + 4gbar) tauQ = -(M + 4gbar) + 8gbar P_diag
#       on the locus, machine zero (N = 4, 5, 6); antilinear variant must FAIL
#   G2  corner-block frozen census: multiplicity floor(N/2) at lambda = -4gbar
#       for J in {0.6, 1.0, 2.3} (N = 3..6); plus the FULL (p,q) block census at
#       N = 4, 5, 6 against both roots: exactly the four corner blocks carry
#       floor(N/2), each at the root its fold parity selects, all others none
#   G3  the pencil kernel dim(ker of (E', K_{D-})) on O+ equals floor(N/2)
#   G4  eigenvector predictions: v_D = 0 and tauQ-even O-part (machine zero);
#       plus dim(ker K cap anti-diagonal of O) = 1 at even N, 0 at odd (N = 3..8)
#   G5  partial balance yields nothing (N = 5 with the middle off the common
#       mean; N = 6 with two of three pair sums equal)
#   G6  N=3 closed form (sympy, exact): det(-2g2 I - M_(1,2)) =
#       512 J^4 (g1+g3)^2 (g1+g3-2g2) (4J^2 + (g1+g3)(g1+g3-2g2))
#   G7  N=6 exact arithmetic (Gaussian rationals): on-locus 36x36 corner det
#       is exactly zero; transverse vanishing order is 3
#   G8  the cofactor theorem (proof doc Section 6):
#       (a) exact Gaussian-rational: coeff of eps^m of det(eps I - Mtilde)
#           equals (-1)^N (4 gbar)^ceil(N/2) det((X P_{O+} X)|_{V-}) at N=4,5
#           (Heisenberg, non-AP locus profile, rational J); float check XY N=4
#       (b) symbolic N=3 corner cofactor = 2^12 gbar^2 J^4 (3J^2 - d1^2)
#       (c) leading coefficient det((K P_{O+} K)|_{V-}) nonzero, N = 3..10,
#           Heisenberg and XY
#   G9  the two boundary clocks (proof doc Section 7):
#       (a) SE eigenbasis identification: Heis = DCT-II (lam_k = 4cos(k pi/N)
#           + N - 5, u_k(a) ~ cos((2a-1)k pi/(2N))), XY = DST-I
#           (lam_k = 4cos(k pi/(N+1)), u_k(a) ~ sin(a k pi/(N+1))), N = 3..10
#       (b) BB^T = (1-1/M) 11^T/N + (I+R)/(2M), M = N (Heis) / N+1 (XY);
#           G = B^T B spectrum {1, 1/M x (ceil(N/2)-1), 0};
#           pdet(G) = M^{-floor((N-1)/2)}
#       (c) D_N = (-1)^{N(N-1)/2} prod (lam_i-lam_j)^2 M^{-floor((N-1)/2)}
#           against the direct determinant, N = 3..10 both chains
#       (d) exact rational assembly at the uniform point (N = 4, 5 Heisenberg,
#           sympy discriminant); uniform corner multiplicity = floor(N/2) at
#           every sampled J (N = 3..6)
#   G10 the J-valuation ladder (proof doc Section 6, last consequence):
#       (a) exact discriminators: ord_J of the cofactor = 18 at N=6
#           (= 2 floor(N/2) ceil(N/2), not 16 = 4(N-2)) and = 24 at N=7
#           (not 20), Gaussian-rational at J = 1/1000 vs 1/2000
#       (b) the exact second-order rung: rank(P_anti K P_{D-}) = 1 at even N,
#           0 at odd N (only the distance-1 middle pair is one hop from the
#           anti-diagonal), N = 4..7
#   G11 the valuation law (proof doc Section 8), integer-exact:
#       (a) Lemma 7 (the level grading ell = distance to the anti-diagonal):
#           no Y entry across two levels, J^0 part level-diagonal, N = 3..10
#           both chains; Lemma 8 (the level census) against the closed form
#       (b) Lemma 9 (the transport bound): exhaustively over ALL maximal row
#           sets, ord_J det Y_I >= its transport bound, the bound is >=
#           floor(N^2/4) for every row set, and the measured minimum over row
#           sets is exactly floor(N^2/4) (N = 3..6)
#       (c) the theorem: ord_J det((X P_{O+} X)|_{V-}) = 2 floor(N^2/4),
#           N = 3..10, Heisenberg and XY
#       (d) the per-pair law: ord_J S_(c,c') = d_c + d_c' on the D- Schur
#           complement, all pairs, N = 3..7
#       (e) the corollary: at J = 0 the frozen root carries 2 floor(N/2), twice
#           the generic multiplicity, so exactly floor(N/2) modes depart
#           (float census; (a)-(d) and (f)-(g) are integer/Fraction exact)
#       (f) Lemma 6 entry by entry, (X P_{O+} X)|_{V-} = Y^T D Y against an
#           independent cell-basis build; and the generality of the hypothesis:
#           a third h (non-uniform R-symmetric bonds, arbitrary R-symmetric
#           diagonal) obeys the same grading and the same valuation, N = 4..6
#       (g) the reduction of the open half: the leading Schur matrix is a
#           positive-weight real Gram matrix of the ghat vectors, hence
#           positive semidefinite, and it is positive definite at N = 4..6
#   G12 the exceptional couplings (proof doc Section 9), exact:
#       (a) N=3 over Q(sqrt(3)): the kernel dimensions of Mtilde, Mtilde^2,
#           Mtilde^3 at J* = d1/sqrt(3) are (1, 2, 2) for three profiles, i.e.
#           geometric 1, algebraic 2, one block of size exactly 2; a
#           non-exceptional control gives (1, 1, 1)
#       (b) N=4 at a profile with a RATIONAL exceptional coupling: (2, 3, 3)
#           against (2, 2, 2) at both neighbouring couplings, with the
#           algebraic count confirmed a second time from the characteristic
#           polynomial (two independent routes)
#   G13 the counting question (proof doc Section 12), exact Sturm counts of the
#       real nonzero roots: 1 pair at N=3 (both profiles), 2 at N=4, 2 at N=5,
#       and 2 versus 4 for two generic N=6 profiles, so the count is NOT a
#       function of N
#   G14 the pointed grading (proof doc Section 8.1): chi_x admissible and
#       constant on every parity basis vector (N = 3..9); chi_x(c,c) =
#       max(d_c, d_x) = the L1 walk length; the census collapse
#       F_j = [chi_x(c,c) >= j]; the order law ord = max(d_c, d_x) for
#       every (cell, pair); each support = the interval between that
#       pair's own sites; and the strict nesting giving triangularity
#
# Runtime: ~8-12 min. Standalone except G0 (imports framework once).
import sys
import math
import random
from fractions import Fraction
from itertools import combinations

import numpy as np

random.seed(2026)
FAILURES = []


def check(name, ok, detail=""):
    status = "ok" if ok else "FAIL"
    print(f"  [{status}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILURES.append(name)


# ---------- shared builders (site 0 = leftmost bit) ----------

def build_H(N):
    """Heisenberg XXX open chain, J = 1, Pauli convention (XX+YY+ZZ per bond)."""
    d = 2 ** N
    H = np.zeros((d, d))
    for a in range(N - 1):
        pa, pb = N - 1 - a, N - 2 - a
        for i in range(d):
            ba, bb = (i >> pa) & 1, (i >> pb) & 1
            H[i, i] += (1.0 if ba == bb else -1.0)
            if ba != bb:
                H[i ^ (1 << pa) ^ (1 << pb), i] += 2.0
    return H


def corner_pieces(N, gl, zz=True):
    """(1,1) block in the cell basis (a,b): K (J=1 Hamiltonian part) and
    Gamma (rate matrix, diagonal; Gamma[(a,b)] = g_a + g_b off-diagonal, 0 on
    the diagonal cells). L_block(J) = J*K - 2*Gamma. zz=False drops the ZZ
    diagonal of h (the XY chain): the theorem needs only h real, symmetric,
    R-invariant."""
    Hfull = build_H(N)
    basis = [1 << (N - 1 - a) for a in range(N)]
    h = np.array([[Hfull[basis[a], basis[b]] for b in range(N)] for a in range(N)])
    if not zz:
        np.fill_diagonal(h, 0.0)
    n2 = N * N
    K = np.zeros((n2, n2), dtype=complex)
    G = np.zeros((n2, n2))
    for a in range(N):
        for b in range(N):
            r = a * N + b
            for c in range(N):
                K[r, c * N + b] += -1j * h[a, c]
                K[r, a * N + c] += 1j * h[c, b]
            if a != b:
                G[r, r] = gl[a] + gl[b]
    return K, G


def general_block(N, gl, J, p, q):
    """Any joint-popcount block (p,q), for the census controls."""
    d = 2 ** N
    H = build_H(N) * J
    cells = [(i, j) for i in range(d) for j in range(d)
             if bin(i).count('1') == p and bin(j).count('1') == q]
    idx = {c: k for k, c in enumerate(cells)}
    n = len(cells)
    M = np.zeros((n, n), dtype=complex)
    for a2, (i, j) in enumerate(cells):
        for k in range(d):
            if H[i, k] != 0 and (k, j) in idx:
                M[a2, idx[(k, j)]] += -1j * H[i, k]
            if H[k, j] != 0 and (i, k) in idx:
                M[a2, idx[(i, k)]] += 1j * H[k, j]
        t = i ^ j
        M[a2, a2] += -2 * sum(gl[m] for m in range(N) if (t >> (N - 1 - m)) & 1)
    return M


def frozen_count_at(Mfun, target, J_list=(0.6, 1.0, 2.3), tol=1e-9):
    """Number of eigenvalues equal to `target` at every J in J_list."""
    counts = [int(np.sum(np.abs(np.linalg.eigvals(Mfun(J)) - target) < tol))
              for J in J_list]
    return min(counts), counts


def tauQ_perm(N):
    n2 = N * N
    TQ = np.zeros((n2, n2))
    for a in range(N):
        for b in range(N):
            TQ[(N - 1 - b) * N + (N - 1 - a), a * N + b] = 1.0
    return TQ


def r90_profile(N):
    """A non-AP R90-fixed profile: pairs sum to 2*gbar, gbar = 0.09."""
    tg = 0.18
    half = [0.04, 0.06, 0.07][:N // 2]
    return half + ([tg / 2] if N % 2 else []) + [tg - g for g in reversed(half)]


# ---------- G0: cross-check the block builder against the framework ----------

print("G0  block builder vs framework Liouvillian (N=4, generic gammas)")
sys.path.insert(0, 'simulations')
import framework as fw  # noqa: E402
gl0 = [0.04, 0.09, 0.05, 0.13]
cs = fw.ChainSystem(N=4, gamma_0=0.05, J=1.0, H_type='heisenberg')
Lfw = fw.lindbladian_z_dephasing(cs.H * 0.83, np.array(gl0))
d = 16
cells = [i * d + j for i in range(d) for j in range(d)
         if bin(i).count('1') == 1 and bin(j).count('1') == 1]
sub = Lfw[np.ix_(cells, cells)]
K0, G0m = corner_pieces(4, gl0)
own = 0.83 * K0 - 2 * G0m
# cell orders coincide: both enumerate (i,j) with ascending i then j, and the
# corner basis (a,b) maps to i = e_a, j = e_b with a ascending = i descending;
# align via the basis map before comparing.
basis = [1 << (4 - 1 - a) for a in range(4)]
order = sorted(range(16), key=lambda r: (basis[r // 4], basis[r % 4]))
own_sorted = own[np.ix_(order, order)]
check("corner block equals framework sub-block", np.max(np.abs(own_sorted - sub)) < 1e-12,
      f"max dev {np.max(np.abs(own_sorted - sub)):.1e}")

# ---------- G1: the mirror identity ----------

print("G1  mirror identity on the locus (linear tauQ; antilinear must fail)")
for N in (4, 5, 6):
    gl = r90_profile(N)
    gbar = sum(gl) / N
    K, G = corner_pieces(N, gl)
    M = 0.83 * K - 2 * G
    TQ = tauQ_perm(N)
    Pd = np.zeros((N * N, N * N))
    for a in range(N):
        Pd[a * N + a, a * N + a] = 1.0
    Mc = M + 4 * gbar * np.eye(N * N)
    lin = np.max(np.abs(TQ @ Mc @ TQ - (-Mc + 8 * gbar * Pd)))
    anti = np.max(np.abs(TQ @ Mc.conj() @ TQ - (-Mc + 8 * gbar * Pd)))
    check(f"N={N} linear identity", lin < 1e-12, f"residual {lin:.1e}")
    check(f"N={N} antilinear variant fails", anti > 1e-2, f"residual {anti:.1e}")

# ---------- G2: frozen census ----------

print("G2  corner-block multiplicity floor(N/2) at -4*gbar; full (p,q) census")
for N in (3, 4, 5, 6):
    gl = r90_profile(N)
    gbar = sum(gl) / N
    K, G = corner_pieces(N, gl)
    mmin, counts = frozen_count_at(lambda J: J * K - 2 * G, -4 * gbar)
    check(f"N={N} corner multiplicity = {N // 2}", mmin == N // 2, f"counts {counts}")
# G2a2: the FULL block census, not one control block. Over every joint-popcount
# block (p,q), count eigenvalues at BOTH frozen roots (-4gbar, and its one-sided
# fold image 4gbar-2sigma). Expected: exactly the four corner blocks p,q in
# {1, N-1} carry floor(N/2), with the root selected by the parity of how many
# one-sided folds separate the block from (1,1) (each fold sends r -> -r-2sigma);
# every other block carries nothing at either root.
for N in (4, 5, 6):
    gl = r90_profile(N)
    gbar = sum(gl) / N
    sig = sum(gl)
    m = N // 2
    roots = (-4 * gbar, 4 * gbar - 2 * sig)
    bad = []
    for p in range(N + 1):
        for q in range(N + 1):
            got = [N * N, N * N]
            for J in (0.6, 1.0, 2.3):
                ev = np.linalg.eigvals(general_block(N, gl, J, p, q))
                for t, r in enumerate(roots):
                    got[t] = min(got[t], int(np.sum(np.abs(ev - r) < 1e-9)))
            corner = p in (1, N - 1) and q in (1, N - 1)
            # number of one-sided folds from (1,1): one per index sitting at N-1
            folds = (p == N - 1) + (q == N - 1)
            want = [0, 0]
            if corner:
                want[folds % 2] = m
            if abs(roots[0] - roots[1]) < 1e-12:      # N = 4: the roots coincide
                if got != ([m, m] if corner else [0, 0]):
                    bad.append((p, q, got))
            elif got != want:
                bad.append((p, q, got, want))
    check(f"N={N} full (p,q) census: only the four corners, root by fold parity",
          not bad, f"deviations {bad}")
# G2b: XY variant (h without the ZZ diagonal): same multiplicity
for N in (4, 5):
    gl = r90_profile(N)
    gbar = sum(gl) / N
    Kxy, Gxy = corner_pieces(N, gl, zz=False)
    mmin, counts = frozen_count_at(lambda J: J * Kxy - 2 * Gxy, -4 * gbar)
    check(f"N={N} XY corner multiplicity = {N // 2}", mmin == N // 2, f"counts {counts}")
# ---------- G3 + G4: pencil kernel dims and eigenvector form ----------

print("G3/G4  pencil kernel dims; v_D = 0 and tauQ-even O-part")
for N in (4, 5, 6):
    gl = r90_profile(N)
    gbar = sum(gl) / N
    n2 = N * N
    K, G = corner_pieces(N, gl)
    O = [a * N + b for a in range(N) for b in range(N) if a != b]
    D = [a * N + a for a in range(N)]
    TQ = tauQ_perm(N)

    def parity_basis(idx, sign):
        sub = TQ[np.ix_(idx, idx)]
        w, V = np.linalg.eigh(sub)
        return V[:, np.abs(w - sign) < 1e-9]

    Op = parity_basis(O, +1)
    Dm = parity_basis(D, -1)
    Ep = (K - 2 * G + 4 * gbar * np.eye(n2))[np.ix_(O, O)]
    A = np.vstack([Ep @ Op, Dm.conj().T @ K[np.ix_(D, O)] @ Op])
    s = np.linalg.svd(A, compute_uv=False)
    null = int(np.sum(s < 1e-9 * max(1, s[0])))
    check(f"N={N} pencil kernel dim = {N // 2}", null == N // 2, f"dim {null}")

    M = 1.7 * K - 2 * G
    w, V = np.linalg.eig(M)
    idx = [i for i in range(n2) if abs(w[i] - (-4 * gbar)) < 1e-8]
    ok_D, ok_par = True, True
    for i in idx:
        v = V[:, i]
        ok_D &= np.linalg.norm(v[D]) < 1e-8
        vO = np.zeros(n2, dtype=complex)
        vO[O] = v[O]
        ok_par &= np.linalg.norm(vO - TQ @ vO) < 1e-8 * np.linalg.norm(vO)
    check(f"N={N} frozen eigenvectors: v_D = 0, tauQ-even", ok_D and ok_par)

# G4b: the placement side claim (Section 9, first bullet): the frozen modes are
# not built out of ker(ad_H). dim(ker K cap span of the anti-diagonal cells OF O,
# i.e. (a, Ra) with a != Ra) is 1 at even N and 0 at odd N, far below floor(N/2).
# (Note the reading matters: including the odd-N centre cell, which is diagonal
# and not in O, would give 1 at every N.)
for N in range(3, 9):
    K, _ = corner_pieces(N, [0.0] * N)
    idx = [a * N + (N - 1 - a) for a in range(N) if a != N - 1 - a]
    E = np.zeros((N * N, len(idx)), dtype=complex)
    for j, i in enumerate(idx):
        E[i, j] = 1.0
    s = np.linalg.svd(K @ E, compute_uv=False)
    nul = int(np.sum(s < 1e-10 * max(1.0, s[0])))
    want = 1 if N % 2 == 0 else 0
    check(f"N={N} dim(ker K cap anti-diagonal of O) = {want}", nul == want,
          f"dim {nul}, floor(N/2) = {N // 2}")

# ---------- G5: partial balance yields nothing ----------

print("G5  partial balance yields nothing")
gl = [0.04, 0.06, 0.05, 0.12, 0.14]        # s1 = s2 = 0.18, middle off the mean
gbar_c = 0.09                               # the would-be common mean rate scale
K, G = corner_pieces(5, gl)
mmin, counts = frozen_count_at(lambda J: J * K - 2 * G, -2 * 0.18)
check("N=5 s1=s2 but middle off: none", mmin == 0, f"counts {counts}")
gl = [0.04, 0.06, 0.09, 0.11, 0.12, 0.14]  # s1 = s2 = 0.18, s3 = 0.20
K, G = corner_pieces(6, gl)
mmin, counts = frozen_count_at(lambda J: J * K - 2 * G, -2 * 0.18)
check("N=6 two of three pairs: none", mmin == 0, f"counts {counts}")

# ---------- G6: N=3 closed form, exact ----------

print("G6  N=3 closed form (sympy exact, block (1,2))")
import sympy as sp  # noqa: E402
g1, g2, g3, Js, lam = sp.symbols('g1 g2 g3 J lam')
glist = [g1, g2, g3]
X = sp.Matrix([[0, 1], [1, 0]])
Y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
Z = sp.Matrix([[1, 0], [0, -1]])
I2 = sp.eye(2)


def kron3(a, b, c):
    return sp.Matrix(np.kron(np.kron(np.array(a.tolist(), dtype=object),
                                     np.array(b.tolist(), dtype=object)),
                             np.array(c.tolist(), dtype=object)))


Hs = sp.zeros(8, 8)
for (a, b) in ((0, 1), (1, 2)):
    for P in (X, Y, Z):
        m = [I2, I2, I2]
        m[a] = P
        m[b] = P
        Hs += Js * kron3(*m)
cells = [(i, j) for i in range(8) for j in range(8)
         if bin(i).count('1') == 1 and bin(j).count('1') == 2]
n = len(cells)
Ms = sp.zeros(n, n)
for a2, (i, j) in enumerate(cells):
    for b2, (k, l) in enumerate(cells):
        v = -sp.I * (Hs[i, k] * (1 if j == l else 0) - (1 if i == k else 0) * Hs[l, j])
        if a2 == b2:
            t = i ^ j
            v += -2 * sum(glist[m2] for m2 in range(3) if (t >> (2 - m2)) & 1)
        Ms[a2, b2] = v
detv = sp.expand((-2 * g2 * sp.eye(n) - Ms).det(method='berkowitz'))
closed = 512 * Js**4 * (g1 + g3)**2 * (g1 + g3 - 2 * g2) * \
    (4 * Js**2 + (g1 + g3) * (g1 + g3 - 2 * g2))
check("det(-2 g2 I - M_(1,2)) equals the closed form",
      sp.simplify(detv - sp.expand(closed)) == 0)

# ---------- G7: N=6 exact arithmetic ----------

print("G7  N=6 exact Gaussian-rational: on-locus zero + transverse order 3")


class GQ:
    __slots__ = ('a', 'b')

    def __init__(self, a=0, b=0):
        self.a = Fraction(a)
        self.b = Fraction(b)

    def __add__(s, o):
        return GQ(s.a + o.a, s.b + o.b)

    def __sub__(s, o):
        return GQ(s.a - o.a, s.b - o.b)

    def __mul__(s, o):
        return GQ(s.a * o.a - s.b * o.b, s.a * o.b + s.b * o.a)

    def __truediv__(s, o):
        dd = o.a * o.a + o.b * o.b
        return GQ((s.a * o.a + s.b * o.b) / dd, (s.b * o.a - s.a * o.b) / dd)

    def is_zero(s):
        return s.a == 0 and s.b == 0


def det_exact(M):
    n = len(M)
    M = [row[:] for row in M]
    det = GQ(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if not M[r][c].is_zero()), None)
        if piv is None:
            return GQ(0)
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            det = det * GQ(-1)
        det = det * M[c][c]
        inv = GQ(1) / M[c][c]
        for r in range(c + 1, n):
            if M[r][c].is_zero():
                continue
            f = M[r][c] * inv
            for c2 in range(c, n):
                M[r][c2] = M[r][c2] - f * M[c][c2]
    return det


def corner_exact(N, gv, Jv, root):
    """root*I - L_block on the (1,1) corner, exact (rates and hops rational)."""
    Hfull = build_H(N)
    basis = [1 << (N - 1 - a) for a in range(N)]
    h = [[Fraction(int(Hfull[basis[a], basis[b]])) for b in range(N)] for a in range(N)]
    n2 = N * N
    M = [[GQ(0) for _ in range(n2)] for _ in range(n2)]
    for a in range(N):
        for b in range(N):
            r = a * N + b
            for c in range(N):
                M[r][c * N + b] = M[r][c * N + b] + GQ(0, -Jv * h[a][c])
                M[r][a * N + c] = M[r][a * N + c] + GQ(0, Jv * h[c][b])
            if a != b:
                M[r][r] = M[r][r] + GQ(-2 * (gv[a] + gv[b]), 0)
    A = [[GQ(0) for _ in range(n2)] for _ in range(n2)]
    for r in range(n2):
        for c in range(n2):
            A[r][c] = GQ(-M[r][c].a, -M[r][c].b)
        A[r][r] = A[r][r] + GQ(root, 0)
    return A


tg = Fraction(9, 50)
q1, q2, q3 = Fraction(2, 50), Fraction(3, 50), Fraction(7, 100)
Jq = Fraction(5, 7)


def prof6(eps):
    return [q1, q2, q3, tg - q3, tg - q2, tg - q1 + eps]


def det6(eps):
    gv = prof6(eps)
    return det_exact(corner_exact(6, gv, Jq, -2 * tg))


v0 = det6(Fraction(0))
check("on-locus corner det exactly zero", v0.is_zero())
e1, e2 = Fraction(1, 10**6), Fraction(1, 2 * 10**6)
d1, d2 = det6(e1), det6(e2)


def mag(gq):
    return math.sqrt(float(gq.a)**2 + float(gq.b)**2)


order = math.log(mag(d1) / mag(d2)) / math.log(2.0)
check("transverse vanishing order = 3", abs(order - 3.0) < 0.01, f"order {order:.4f}")

# ---------- G8: the cofactor theorem ----------

print("G8  cofactor closed form: exact N=4,5; symbolic N=3; leading coeff N=3..10")


def tq_index(N):
    return [(N - 1 - (r % N)) * N + (N - 1 - (r // N)) for r in range(N * N)]


def vm_columns(N):
    """Unnormalized V- basis: e_r - e_{tq(r)} over tq-2-cycles (D- and O-)."""
    tq = tq_index(N)
    cols, seen = [], set()
    for r in range(N * N):
        rp = tq[r]
        if rp == r or r in seen:
            continue
        seen.add(r)
        seen.add(rp)
        cols.append((r, rp))
    return cols


def se_h_frac(N):
    h = [[Fraction(0)] * N for _ in range(N)]
    for a in range(N - 1):
        h[a][a + 1] = h[a + 1][a] = Fraction(2)
    for a in range(N):
        h[a][a] = Fraction(sum(-1 if (a == b or a == b + 1) else 1
                               for b in range(N - 1)))
    return h


def gq_matmul(A, B):
    n, k, m2 = len(A), len(B), len(B[0])
    C = [[GQ(0) for _ in range(m2)] for _ in range(n)]
    for i in range(n):
        for t in range(k):
            a = A[i][t]
            if a.is_zero():
                continue
            for j in range(m2):
                if not B[t][j].is_zero():
                    C[i][j] += a * B[t][j]
    return C


def cofactor_exact_direct(N, gv, Jv):
    """coeff of eps^m of det(eps I - Mtilde): dets at n2+1 nodes + Newton
    interpolation, exact; asserts coeffs eps^0..eps^(m-1) vanish."""
    n2, m = N * N, N // 2
    gbar = sum(gv) / N
    # corner_exact(root) returns root*I - L_block, so eps*I - Mtilde =
    # (eps - 4*gbar)*I - L_block = corner_exact(root = eps - 4*gbar)
    nodes = [Fraction(k, 7) for k in range(1, n2 + 2)]
    dets = []
    for e in nodes:
        A = corner_exact(N, gv, Jv, e - 4 * gbar)
        dets.append(det_exact(A))

    def interp(vals):
        n = len(nodes)
        dd = [list(vals)]
        for k in range(1, n):
            prev = dd[-1]
            dd.append([(prev[i + 1] - prev[i]) / (nodes[i + k] - nodes[i])
                       for i in range(n - k)])
        coeffs = [Fraction(0)] * n
        basis = [Fraction(1)] + [Fraction(0)] * (n - 1)
        for k in range(n):
            c = dd[k][0]
            for i in range(k + 1):
                coeffs[i] += c * basis[i]
            if k == n - 1:
                break
            new = [Fraction(0)] * n
            for i in range(k + 1):
                new[i + 1] += basis[i]
                new[i] += -nodes[k] * basis[i]
            basis = new
        return coeffs

    ca = interp([dv.a for dv in dets])
    cb = interp([dv.b for dv in dets])
    for k in range(m):
        assert ca[k] == 0 and cb[k] == 0, f"eps^{k} coeff nonzero"
    assert cb[m] == 0, "cofactor not real"
    return ca[m]


def cofactor_exact_formula(N, gv, Jv):
    """(-1)^N (4 gbar)^ceil(N/2) det((X P+ X)|_{V-}), exact."""
    n2 = N * N
    h = se_h_frac(N)
    gbar = sum(gv) / N
    X = [[GQ(0) for _ in range(n2)] for _ in range(n2)]
    for a in range(N):
        for b in range(N):
            r = a * N + b
            for c in range(N):
                X[r][c * N + b] += GQ(0, -Jv * h[a][c])
                X[r][a * N + c] += GQ(0, Jv * h[c][b])
            if a != b:
                X[r][r] += GQ(-2 * (gv[a] + gv[b]) + 4 * gbar)
    tq = tq_index(N)
    P = [[GQ(0) for _ in range(n2)] for _ in range(n2)]
    for r in range(n2):
        if r // N == r % N:
            continue
        P[r][r] += GQ(Fraction(1, 2))
        P[r][tq[r]] += GQ(Fraction(1, 2))
    Y = gq_matmul(gq_matmul(X, P), X)
    cols = vm_columns(N)
    vm = len(cols)
    Red = [[Y[r1][c1] - Y[r1][c2] - Y[r2][c1] + Y[r2][c2]
            for (c1, c2) in cols] for (r1, r2) in cols]
    dv = det_exact(Red)
    assert dv.b == 0, "reduced det not real"
    sign = 1 if N % 2 == 0 else -1
    return sign * (4 * gbar) ** (N - N // 2) * dv.a / Fraction(2) ** vm


for N, deltas, Jv in [
    (4, [Fraction(1, 50), Fraction(-3, 100)], Fraction(9, 5)),
    (5, [Fraction(1, 50), Fraction(3, 100)], Fraction(4, 3)),
]:
    gbq = Fraction(9, 100)
    gv = [gbq + d for d in deltas] + ([gbq] if N % 2 else []) + \
         [gbq - d for d in reversed(deltas)]
    c_dir = cofactor_exact_direct(N, gv, Jv)
    c_form = cofactor_exact_formula(N, gv, Jv)
    check(f"N={N} exact cofactor = closed form", c_dir == c_form,
          f"value {c_dir}")

# G8 float check, XY N=4: formula vs product of nonzero eigenvalues route
glxy = r90_profile(4)
gbxy = sum(glxy) / 4
Kxy, Gxy = corner_pieces(4, glxy, zz=False)
Mtxy = 1.7 * Kxy - 2 * Gxy + 4 * gbxy * np.eye(16)
pxy = np.poly(np.linalg.eigvals(Mtxy))
c_dir_xy = pxy[16 - 2].real
TQ4 = tauQ_perm(4)
Pd4 = np.diag([1.0 if r // 4 == r % 4 else 0.0 for r in range(16)])
Xxy = Mtxy - 4 * gbxy * Pd4
Ppl = 0.5 * ((np.eye(16) - Pd4) + TQ4 @ (np.eye(16) - Pd4))
# V- basis (orthonormal, float)
wv, Vv = np.linalg.eigh(TQ4)
Vm4 = Vv[:, np.abs(wv + 1) < 1e-9]
c_form_xy = (4 * gbxy) ** 2 * np.linalg.det(
    Vm4.conj().T @ (Xxy @ Ppl @ Xxy) @ Vm4).real
check("N=4 XY float cofactor = closed form",
      abs(c_dir_xy - c_form_xy) < 1e-8 * abs(c_dir_xy),
      f"rel dev {abs(c_dir_xy - c_form_xy) / abs(c_dir_xy):.1e}")

# G8b: symbolic N=3 corner cofactor
g1s, g3s = sp.symbols('g1c g3c', positive=True)
gbar_s = (g1s + g3s) / 2
gl_s = [g1s, gbar_s, g3s]
h3 = sp.Matrix([[0, 2, 0], [2, -2, 2], [0, 2, 0]])
M3 = sp.zeros(9, 9)
for a in range(3):
    for b in range(3):
        r = a * 3 + b
        for c in range(3):
            M3[r, c * 3 + b] += -sp.I * Js * h3[a, c]
            M3[r, a * 3 + c] += sp.I * Js * h3[c, b]
        if a != b:
            M3[r, r] += -2 * (gl_s[a] + gl_s[b])
Mt3 = M3 + 4 * gbar_s * sp.eye(9)
epss = sp.symbols('epss')
p3 = sp.Poly(sp.expand(Mt3.charpoly(epss).as_expr()), epss)
co = p3.all_coeffs()[::-1]
d1s = (g1s - g3s) / 2
target = 2**12 * gbar_s**2 * Js**4 * (3 * Js**2 - d1s**2)
check("N=3 symbolic: eps^0 coeff zero", sp.simplify(co[0]) == 0)
check("N=3 symbolic corner cofactor = 2^12 gbar^2 J^4 (3J^2 - d1^2)",
      sp.simplify(sp.expand(co[1]) - sp.expand(target)) == 0)

# G8c: leading coefficient nonzero N=3..10, Heisenberg and XY (float)
for N in range(3, 11):
    n2 = N * N
    TQn = tauQ_perm(N)
    Pdn = np.diag([1.0 if r // N == r % N else 0.0 for r in range(n2)])
    Ppn = 0.5 * ((np.eye(n2) - Pdn) + TQn @ (np.eye(n2) - Pdn))
    wv, Vv = np.linalg.eigh(TQn)
    Vmn = Vv[:, np.abs(wv + 1) < 1e-9]
    for zz, lbl in ((True, "Heis"), (False, "XY")):
        K, _ = corner_pieces(N, [0.0] * N, zz=zz)
        Redn = Vmn.conj().T @ (K @ Ppn @ K) @ Vmn
        s = np.linalg.svd(Redn, compute_uv=False)
        check(f"N={N} {lbl} leading coeff det((K P+ K)|V-) nonzero",
              s[-1] > 1e-8 * s[0], f"smin/smax {s[-1] / s[0]:.1e}")

# ---------- G9: the two clocks ----------

print("G9  the two boundary clocks: DCT-II/DST-I, BB^T law, pdet(G), D_N, "
      "uniform assembly")


def se_h_np(N, zz=True):
    """Single-excitation h (hopping 2 + optional ZZ diagonal), float."""
    h = np.zeros((N, N))
    for a in range(N - 1):
        h[a, a + 1] = h[a + 1, a] = 2.0
    if zz:
        for a in range(N):
            h[a, a] = sum(-1 if (a == b or a == b + 1) else 1
                          for b in range(N - 1))
    return h


for zz, lbl in ((True, "Heis"), (False, "XY")):
    dev_clock = dev_bbt = dev_spec = dev_dn = 0.0
    for N in range(3, 11):
        M = N if zz else N + 1
        h = se_h_np(N, zz)
        lam, U = np.linalg.eigh(h)
        # (a) clock identification
        if zz:
            ks = np.arange(N)
            lpred = 4 * np.cos(ks * np.pi / N) + N - 5
            vecs = [np.cos((2 * np.arange(1, N + 1) - 1) * k * np.pi / (2 * N))
                    for k in ks]
        else:
            ks = np.arange(1, N + 1)
            lpred = 4 * np.cos(ks * np.pi / (N + 1))
            vecs = [np.sin(np.arange(1, N + 1) * k * np.pi / (N + 1)) for k in ks]
        for lk, v in zip(lpred, vecs):
            idx = int(np.argmin(np.abs(lam - lk)))
            v = v / np.linalg.norm(v)
            u = U[:, idx]
            dev_clock = max(dev_clock, abs(lam[idx] - lk),
                            min(np.linalg.norm(u - v), np.linalg.norm(u + v)))
        # (b) BB^T law + G spectrum + pdet
        B = U ** 2
        Rm = np.eye(N)[::-1]
        cand = (1 - 1 / M) * np.ones((N, N)) / N + (np.eye(N) + Rm) / (2 * M)
        dev_bbt = max(dev_bbt, np.max(np.abs(B @ B.T - cand)))
        G = B.T @ B
        w = np.sort(np.linalg.eigvalsh(G))[::-1]
        p = (N + 1) // 2
        dev_spec = max(dev_spec, abs(w[0] - 1),
                       float(np.max(np.abs(w[1:p] - 1 / M))),
                       float(np.max(np.abs(w[p:]))),
                       abs(float(np.prod(w[:p])) - M ** (-((N - 1) // 2))))
        # (c) D_N closed form vs direct
        K, _ = corner_pieces(N, [0.0] * N, zz=zz)
        n2 = N * N
        TQn = tauQ_perm(N)
        wv, Vv = np.linalg.eigh(TQn)
        Vmn = Vv[:, np.abs(wv + 1) < 1e-9]
        Pdn = np.diag([1.0 if r // N == r % N else 0.0 for r in range(n2)])
        Ppn = 0.5 * ((np.eye(n2) - Pdn) + TQn @ (np.eye(n2) - Pdn))
        sgn, ld = np.linalg.slogdet(Vmn.conj().T @ (K @ Ppn @ K) @ Vmn)
        disc = np.prod([(lam[i] - lam[j]) ** 2 for i in range(N)
                        for j in range(i + 1, N)])
        rhs = (-1) ** (N * (N - 1) // 2) * disc * M ** (-((N - 1) // 2))
        dev_dn = max(dev_dn, abs(ld - np.log(abs(rhs))),
                     0.0 if sgn.real * np.sign(rhs) > 0 else 1.0)
    check(f"{lbl} clock identification N=3..10", dev_clock < 1e-10,
          f"max dev {dev_clock:.1e}")
    check(f"{lbl} BB^T law N=3..10", dev_bbt < 1e-12, f"max dev {dev_bbt:.1e}")
    check(f"{lbl} G spectrum + pdet law N=3..10", dev_spec < 1e-10,
          f"max dev {dev_spec:.1e}")
    check(f"{lbl} D_N closed form N=3..10", dev_dn < 1e-8, f"max dev {dev_dn:.1e}")

# (d) exact rational assembly at the uniform point, N = 4, 5 Heisenberg
for N, Jv in ((4, Fraction(9, 5)), (5, Fraction(4, 3))):
    gbq = Fraction(9, 100)
    c_dir = cofactor_exact_direct(N, [gbq] * N, Jv)
    hs = sp.Matrix(N, N, lambda a, b: sp.Integer(int(se_h_np(N)[a, b])))
    xs = sp.symbols('xdisc')
    disc_i = sp.discriminant(hs.charpoly(xs).as_expr(), xs)
    p = (N + 1) // 2
    rhs_q = (Fraction(-1) ** N) * (4 * gbq) ** p * Jv ** (N * (N - 1)) \
        * (Fraction(-1) ** (N * (N - 1) // 2)) * Fraction(int(disc_i)) \
        / Fraction(N) ** ((N - 1) // 2)
    check(f"N={N} exact uniform cofactor = clock closed form", c_dir == rhs_q,
          f"value {c_dir}")

# uniform tightness at every sampled J
for N in (3, 4, 5, 6):
    K, G9G = corner_pieces(N, [0.07] * N)
    gbar9 = 0.07
    mmin, counts = frozen_count_at(lambda J: J * K - 2 * G9G, -4 * gbar9,
                                   J_list=(0.3, 1.0, 2.7))
    check(f"N={N} uniform corner multiplicity = {N // 2}",
          all(c == N // 2 for c in counts), f"counts {counts}")

# ---------- G10: the J-valuation ladder ----------

print("G10 J-valuation ladder: exact N=6/N=7 discriminators; second-order rung")

# (a) exact discriminator at N=6
gbq10 = Fraction(9, 100)
d10 = [Fraction(1, 50), Fraction(-3, 100), Fraction(1, 25)]
gv10 = [gbq10 + d for d in d10] + [gbq10 - d for d in reversed(d10)]
cJ1 = cofactor_exact_formula(6, gv10, Fraction(1, 1000))
cJ2 = cofactor_exact_formula(6, gv10, Fraction(1, 2000))
ordJ6 = math.log2(abs(cJ1 / cJ2))
check("N=6 exact J-valuation = 18 (2 floor(N/2) ceil(N/2)), not 16",
      abs(ordJ6 - 18) < 0.05, f"measured {ordJ6:.4f}")

# (a2) second discriminator at N=7: 24 (= 2 m (N-m)), not 20 (= 4(N-2));
# needs J deep below every delta scale (1/100 still overshoots), so 1/1000
d10b = [Fraction(1, 50), Fraction(3, 100), Fraction(-1, 25)]
gv10b = [gbq10 + d for d in d10b] + [gbq10] + [gbq10 - d for d in reversed(d10b)]
cJ1b = cofactor_exact_formula(7, gv10b, Fraction(1, 1000))
cJ2b = cofactor_exact_formula(7, gv10b, Fraction(1, 2000))
ordJ7 = math.log2(abs(cJ1b / cJ2b))
check("N=7 exact J-valuation = 24 (not 20)", abs(ordJ7 - 24) < 0.1,
      f"measured {ordJ7:.4f}")

# (b) the exact second-order rung: the D- Schur complement at order J^2 is
# -C(dagger)C with C = P_anti K P_{D-}; rank C = 1 at even N (only the
# distance-1 middle pair reaches the anti-diagonal in one hop), 0 at odd N
for N in (4, 5, 6, 7):
    K, _ = corner_pieces(N, [0.0] * N)
    n2 = N * N
    anti = [a * N + (N - 1 - a) for a in range(N) if a != N - 1 - a]
    TQ = tauQ_perm(N)
    D = [r for r in range(n2) if r // N == r % N]

    def pb(idx, sign):
        sub = TQ[np.ix_(idx, idx)]
        w2, V2 = np.linalg.eigh(sub)
        cols = V2[:, np.abs(w2 - sign) < 1e-9]
        E = np.zeros((n2, cols.shape[1]))
        E[idx, :] = cols
        return E

    Dm_b = pb(D, -1)
    Cmat = (K @ Dm_b)[anti, :]
    s = np.linalg.svd(Cmat, compute_uv=False)
    rk = int(np.sum(s > 1e-10 * max(1.0, s[0])))
    check(f"N={N} rank(P_anti K P_D-) = {1 if N % 2 == 0 else 0}",
          rk == (1 if N % 2 == 0 else 0), f"rank {rk}")

# ---------- G11: the valuation law (the walk to the anti-diagonal) ----------

print("G11 valuation law: the Gram form, the level grading, the transport count")

# Everything here is integer-exact. X is gamma-bar free, so only the
# antisymmetric part delta of the profile enters, and scaling (z, delta)
# jointly only scales the determinant: integer deltas suffice. With z := -i*J
# the matrix X(z) = z*k - 2*Delta is REAL, k = i*K the one-hop part, and
# ord_z = ord_J.

WC_DELTA = {3: [2], 4: [2, -3], 5: [2, 3], 6: [2, -3, 5], 7: [2, -3, 5],
            8: [2, -3, 5, -7], 9: [2, -3, 5, -7], 10: [2, -3, 5, -7, 11]}


def wc_delta(N):
    """Antisymmetric integer delta on the R90 locus."""
    d = [0] * N
    for i, v in enumerate(WC_DELTA[N][:N // 2]):
        d[i], d[N - 1 - i] = v, -v
    return d


def wc_se_h(N, zz=True):
    h = [[0] * N for _ in range(N)]
    for a in range(N - 1):
        h[a][a + 1] = h[a + 1][a] = 2
    if zz:
        for a in range(N):
            h[a][a] = sum(-1 if (a == b or a == b + 1) else 1
                          for b in range(N - 1))
    return h


def wc_build(N, zz=True, h=None, delta=None):
    """Y_{y,x} = u_y^T X w_x as [c0, c1] in z, with row/column levels.

    Rows: the O+ basis (anti-diagonal cells alone, off-diagonal tauQ pairs
    summed). Columns: the V- basis (tauQ pairs differenced). Levels are
    ell(a,b) = |a + b - (N+1)|, the distance of the cell to the anti-diagonal.
    h and delta override the defaults (used for the generality check).
    """
    h = wc_se_h(N, zz) if h is None else h
    delta = wc_delta(N) if delta is None else delta
    n2 = N * N
    k = [[0] * n2 for _ in range(n2)]
    dl = [0] * n2
    for a in range(N):
        for b in range(N):
            r = a * N + b
            for c in range(N):
                k[c * N + b][r] += h[a][c]
                k[a * N + c][r] -= h[c][b]
            if a != b:
                dl[r] = delta[a] + delta[b]
    tq = [(N - 1 - r % N) * N + (N - 1 - r // N) for r in range(n2)]
    ell = [abs(r // N + r % N - (N - 1)) for r in range(n2)]
    vm, op = [], []
    for r in range(n2):
        a, b, t = r // N, r % N, tq[r]
        if t != r and r < t:
            vm.append(((r, t), (1, -1), ell[r]))
        if a == b:
            continue
        if t == r:
            op.append(((r,), (1,), ell[r]))
        elif r < t:
            op.append(((r, t), (1, 1), ell[r]))
    Y = [[[0, 0] for _ in vm] for _ in op]
    for iy, (sy, cy, _) in enumerate(op):
        for ix, (sx, cx, _) in enumerate(vm):
            c0 = c1 = 0
            for ry, ay in zip(sy, cy):
                for rx, ax in zip(sx, cx):
                    c1 += ay * ax * k[ry][rx]
                    if ry == rx:
                        c0 -= 2 * ay * ax * dl[rx]
            Y[iy][ix] = [c0, c1]
    diag_cols = [i for i, (s, _, _) in enumerate(vm) if s[0] // N == s[0] % N]
    return (Y, [lv for _, _, lv in op], [lv for _, _, lv in vm],
            [len(s) for s, _, _ in op], diag_cols)


def wc_det_int(M):
    """Bareiss fraction-free determinant of an integer matrix."""
    n = len(M)
    if n == 0:
        return 1
    M = [row[:] for row in M]
    sign, prev = 1, 1
    for c in range(n - 1):
        if M[c][c] == 0:
            piv = next((r for r in range(c + 1, n) if M[r][c] != 0), None)
            if piv is None:
                return 0
            M[c], M[piv] = M[piv], M[c]
            sign = -sign
        for r in range(c + 1, n):
            for c2 in range(c + 1, n):
                M[r][c2] = (M[r][c2] * M[c][c] - M[r][c] * M[c][c2]) // prev
            M[r][c] = 0
        prev = M[c][c]
    return sign * M[n - 1][n - 1]


def wc_det_frac(M):
    n = len(M)
    M = [row[:] for row in M]
    det = Fraction(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if M[r][c] != 0), None)
        if piv is None:
            return Fraction(0)
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            det = -det
        det *= M[c][c]
        for r in range(c + 1, n):
            if M[r][c] == 0:
                continue
            f = M[r][c] / M[c][c]
            for c2 in range(c, n):
                M[r][c2] -= f * M[c][c2]
    return det


def wc_ord_from_values(vals, deg):
    """ord_z of the degree<=deg polynomial through (i, vals[i]), exactly."""
    if all(v == 0 for v in vals):
        return None
    m = deg + 1
    A = [[Fraction(i ** j) for j in range(m)] + [Fraction(vals[i])]
         for i in range(m)]
    for c in range(m):
        piv = next(r for r in range(c, m) if A[r][c] != 0)
        A[c], A[piv] = A[piv], A[c]
        inv = 1 / A[c][c]
        A[c] = [v * inv for v in A[c]]
        for r in range(m):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [v - f * w for v, w in zip(A[r], A[c])]
    for i in range(m):
        if A[i][m] != 0:
            return i
    return None


def wc_minor_ord(G, rowsel, colsel):
    n = len(colsel)
    vals = [wc_det_int([[G[y][x][0] + z * G[y][x][1] for x in colsel]
                        for y in rowsel]) for z in range(n + 1)]
    return wc_ord_from_values(vals, n)


def wc_gram_ord(G, nrm, ci, cj):
    """ord_z det of the (ci, cj) submatrix of G^T D G, D = 1/||u||^2."""
    n, nr = len(ci), len(nrm)
    vals = []
    for z in range(2 * n + 1):
        cI = [[G[y][x][0] + z * G[y][x][1] for y in range(nr)] for x in ci]
        cJ = [[G[y][x][0] + z * G[y][x][1] for y in range(nr)] for x in cj]
        vals.append(wc_det_frac(
            [[sum(Fraction(cI[i][y] * cJ[j][y], nrm[y]) for y in range(nr))
              for j in range(n)] for i in range(n)]))
    return wc_ord_from_values(vals, 2 * n)


def wc_transport(rows_I, cols):
    """sum_j max(0, F_j): the flow lower bound of Lemma 9."""
    top = max(list(cols) + list(rows_I))
    return sum(max(0, sum(1 for c in cols if c >= j)
                   - sum(1 for r in rows_I if r >= j))
               for j in range(1, top + 1))


# (a) Lemma 7: the level grading (no entry across two levels; the J^0 part is
#     level-diagonal), and Lemma 8: the level census against the closed form
for zz, lbl in ((True, "Heis"), (False, "XY  ")):
    far = zero = cens = 0
    for N in range(3, 11):
        G, rows, cols, nrm, _ = wc_build(N, zz)
        for y in range(len(rows)):
            for x in range(len(cols)):
                c0, c1 = G[y][x]
                d = abs(rows[y] - cols[x])
                far += 1 if (d >= 2 and (c0 or c1)) else 0
                zero += 1 if (d >= 1 and c0) else 0
        m = N // 2
        for j in range(N):
            eps = 1 if (j >= 1 and (N - 1 - j) % 2 == 0) else 0
            n_pred = 0 if j == 0 else N - j
            M_pred = 2 * m if j == 0 else N - j - eps
            cens += 0 if (sum(1 for c in cols if c == j) == n_pred and
                          sum(1 for r in rows if r == j) == M_pred) else 1
    check(f"{lbl} grading: no G entry across two levels, N=3..10", far == 0,
          f"violations {far}")
    check(f"{lbl} grading: the J^0 part is level-diagonal, N=3..10", zero == 0,
          f"violations {zero}")
    check(f"{lbl} census n_j / M_j = closed form, N=3..10", cens == 0,
          f"deviations {cens}")

# the census identity behind the bound: sum over diagonal levels = floor(N^2/4)
lvl_sum = all(sum(j for j in range(1, N) if (N - 1 - j) % 2 == 0)
              == (N * N) // 4 for N in range(3, 21))
check("sum of the diagonal levels d_c = floor(N^2/4), N=3..20", lvl_sum)

# (b) Lemma 9: every maximal minor of G obeys its transport bound, and the
#     minimum over ALL row sets is exactly floor(N^2/4)
for N in range(3, 7):
    G, rows, cols, nrm, _ = wc_build(N, True)
    tgt = (N * N) // 4
    viol, best, wb = 0, None, None
    for I in combinations(range(len(rows)), len(cols)):
        fb = wc_transport([rows[y] for y in I], cols)
        wb = fb if wb is None else min(wb, fb)
        o = wc_minor_ord(G, list(I), list(range(len(cols))))
        if o is None:
            continue
        viol += 1 if o < fb else 0
        best = o if best is None else min(best, o)
    check(f"N={N} every maximal minor of G >= its transport bound", viol == 0,
          f"violations {viol} over {math.comb(len(rows), len(cols))} minors")
    check(f"N={N} transport bound over all row sets >= floor(N^2/4) = {tgt}",
          wb >= tgt, f"min bound {wb}")
    check(f"N={N} measured min ord_J det G_I = {tgt}", best == tgt,
          f"measured {best}")

# (c) the theorem: ord_J det((X P_{O+} X)|_{V-}) = 2 floor(N^2/4)
for zz, lbl in ((True, "Heis"), (False, "XY  ")):
    devs = []
    for N in range(3, 11):
        G, rows, cols, nrm, _ = wc_build(N, zz)
        o = wc_gram_ord(G, nrm, list(range(len(cols))), list(range(len(cols))))
        if o != 2 * ((N * N) // 4):
            devs.append((N, o))
    check(f"{lbl} ord_J det((X P_O+ X)|_V-) = 2 floor(N^2/4), N=3..10",
          not devs, f"deviations {devs}")

# (d) the per-pair law: ord_J S_(c,c') = d_c + d_c' on the D- Schur complement
for N in range(3, 8):
    G, rows, cols, nrm, dcols = wc_build(N, True)
    ocols = [x for x in range(len(cols)) if x not in dcols]
    base = wc_gram_ord(G, nrm, ocols, ocols)
    bad = [(cols[c], cols[c2],
            wc_gram_ord(G, nrm, ocols + [c], ocols + [c2]) - base)
           for c in dcols for c2 in dcols
           if wc_gram_ord(G, nrm, ocols + [c], ocols + [c2]) - base
           != cols[c] + cols[c2]]
    check(f"N={N} T_OO is a unit at J=0 (ord 0)", base == 0, f"ord {base}")
    check(f"N={N} ord_J S_(c,c') = d_c + d_c' for all pairs", not bad,
          f"deviations {bad}")

# (f) Lemma 6 itself, entry by entry: (X P_{O+} X)|_{V-} = Y^T D Y against an
#     independent build of the left side straight in the cell basis; plus the
#     generality of the hypothesis (h real symmetric R-invariant tridiagonal):
#     a third family with NON-UNIFORM R-symmetric bonds and an arbitrary
#     R-symmetric diagonal must obey the same grading and the same valuation
def wc_direct_T(N, delta, z, h):
    """P_{V-} X P_{O+} X P_{V-} built directly in the cell basis (unnormalized
    V- columns), the independent left-hand side of Lemma 6."""
    k, dl = [[0] * (N * N) for _ in range(N * N)], [0] * (N * N)
    for a in range(N):
        for b in range(N):
            r = a * N + b
            for c in range(N):
                k[c * N + b][r] += h[a][c]
                k[a * N + c][r] -= h[c][b]
            if a != b:
                dl[r] = delta[a] + delta[b]
    n2 = N * N
    X = [[Fraction(z) * k[i][j] for j in range(n2)] for i in range(n2)]
    for r in range(n2):
        X[r][r] -= 2 * dl[r]
    tq = [(N - 1 - r % N) * N + (N - 1 - r // N) for r in range(n2)]
    vmc, opc = [], []
    for r in range(n2):
        a, b, t = r // N, r % N, tq[r]
        if t != r and r < t:
            vmc.append((r, t))
        if a == b:
            continue
        opc.append((r,) if t == r else ((r, t) if r < t else None))
    opc = [o for o in opc if o]

    def proj(v):
        out = [Fraction(0)] * n2
        for sup in opc:
            c = sum(v[i] for i in sup)
            for i in sup:
                out[i] += Fraction(c, len(sup))
        return out

    T = []
    for r1, t1 in vmc:
        w = [Fraction(0)] * n2
        w[r1], w[t1] = Fraction(1), Fraction(-1)
        Xw = [sum(X[r][c] * w[c] for c in range(n2)) for r in range(n2)]
        PXw = proj(Xw)
        XPXw = [sum(X[r][c] * PXw[c] for c in range(n2)) for r in range(n2)]
        T.append([XPXw[r0] - XPXw[t0] for r0, t0 in vmc])
    return T


for N, hb, hd in ((4, [3, 5, 3], [2, -6, -6, 2]),
                  (5, [3, 5, 5, 3], [1, 4, 9, 4, 1]),
                  (6, [2, 7, 3, 7, 2], [5, -2, 8, 8, -2, 5])):
    hcus = [[0] * N for _ in range(N)]
    for a in range(N - 1):
        hcus[a][a + 1] = hcus[a + 1][a] = hb[a]
    for a in range(N):
        hcus[a][a] = hd[a]
    check(f"N={N} the third h is R-invariant (bonds and diagonal palindromic)",
          hb == hb[::-1] and hd == hd[::-1])
    for lbl, hh in (("Heisenberg", wc_se_h(N)), ("non-uniform bonds", hcus)):
        G_, rows_, cols_, nrm_, _ = wc_build(N, h=hh)
        zq = Fraction(3, 7)
        Td = wc_direct_T(N, wc_delta(N), zq, hh)
        n_ = len(cols_)
        rhs = [[sum(Fraction((G_[y][i][0] + zq * G_[y][i][1])
                             * (G_[y][j][0] + zq * G_[y][j][1]), nrm_[y])
                    for y in range(len(rows_))) for j in range(n_)]
               for i in range(n_)]
        dev = max(abs(Td[i][j] - rhs[i][j]) for i in range(n_) for j in range(n_))
        far_ = sum(1 for y in range(len(rows_)) for x in range(n_)
                   if abs(rows_[y] - cols_[x]) >= 2
                   and (G_[y][x][0] or G_[y][x][1]))
        o_ = wc_gram_ord(G_, nrm_, list(range(n_)), list(range(n_)))
        check(f"N={N} {lbl}: Lemma 6 identity exact, grading holds, "
              f"ord = 2 floor(N^2/4)",
              dev == 0 and far_ == 0 and o_ == 2 * ((N * N) // 4),
              f"dev {dev}, grading violations {far_}, ord {o_}")

# (g) the open half, reduced: the leading Schur matrix is a POSITIVE-weight
#     real Gram matrix of the vectors ghat_I(c) = [z^{d_c}] det Y_{I, O- + c},
#     so it is positive semidefinite and nonsingular exactly when those
#     floor(N/2) vectors are linearly independent
def wc_coeffs(vals, deg):
    m = deg + 1
    A = [[Fraction(i ** j) for j in range(m)] + [Fraction(vals[i])]
         for i in range(m)]
    for c in range(m):
        piv = next(r for r in range(c, m) if A[r][c] != 0)
        A[c], A[piv] = A[piv], A[c]
        inv = 1 / A[c][c]
        A[c] = [v * inv for v in A[c]]
        for r in range(m):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [v - f * w for v, w in zip(A[r], A[c])]
    return [A[i][m] for i in range(m)]


def wc_lead(G, nrm, ci, cj, want):
    n, nr = len(ci), len(nrm)
    vals = []
    for z in range(2 * n + 1):
        cI = [[G[y][x][0] + z * G[y][x][1] for y in range(nr)] for x in ci]
        cJ = [[G[y][x][0] + z * G[y][x][1] for y in range(nr)] for x in cj]
        vals.append(wc_det_frac(
            [[sum(Fraction(cI[i][y] * cJ[j][y], nrm[y]) for y in range(nr))
              for j in range(n)] for i in range(n)]))
    return wc_coeffs(vals, 2 * n)[want]


for N in (4, 5, 6):
    G, rows, cols, nrm, dcols = wc_build(N, True)
    ocols = [x for x in range(len(cols)) if x not in dcols]
    nR, nC = len(rows), len(ocols) + 1
    gh, wt = {c: {} for c in dcols}, {}
    for I in combinations(range(nR), nC):
        wt[I] = Fraction(1)
        for y in I:
            wt[I] /= nrm[y]
        for c in dcols:
            cs = ocols + [c]
            vals = [wc_det_int([[G[y][x][0] + z * G[y][x][1] for x in cs]
                                for y in I]) for z in range(nC + 1)]
            co = wc_coeffs(vals, nC)
            gh[c][I] = co[cols[c]] if cols[c] < len(co) else Fraction(0)
    Gram = [[sum(wt[I] * gh[c][I] * gh[c2][I] for I in wt) for c2 in dcols]
            for c in dcols]
    Sd = [[wc_lead(G, nrm, ocols + [c], ocols + [c2], cols[c] + cols[c2])
           for c2 in dcols] for c in dcols]
    base = wc_lead(G, nrm, ocols, ocols, 0)
    m = len(dcols)
    check(f"N={N} leading Schur matrix = positive-weight real Gram",
          all(Sd[i][j] == Gram[i][j] for i in range(m) for j in range(m))
          and base > 0, f"det T_OO(0) = {base}")
    S0 = [[Fraction(Sd[i][j], base) for j in range(m)] for i in range(m)]
    minors = [wc_det_frac([row[:kk + 1] for row in S0[:kk + 1]])
              for kk in range(m)]
    check(f"N={N} it is positive definite (all leading minors > 0)",
          all(mm > 0 for mm in minors), f"minors {[str(mm) for mm in minors]}")

# (e) what the order counts: at J = 0 the frozen root carries TWICE the
#     generic multiplicity (the anti-diagonal cells off the centre), so exactly
#     floor(N/2) modes depart as the coupling turns on
for N in range(3, 7):
    gl11 = r90_profile(N)
    K11, G11m = corner_pieces(N, gl11)
    gbar11 = sum(gl11) / N

    def blk(J, K=K11, G=G11m):
        return J * K - 2 * G

    m0, c0 = frozen_count_at(blk, -4 * gbar11, J_list=(0.0,))
    mJ, cJ = frozen_count_at(blk, -4 * gbar11)
    check(f"N={N} frozen root multiplicity at J=0 is 2 floor(N/2)",
          m0 == 2 * (N // 2), f"counts {c0} vs {cJ} at J != 0")
    check(f"N={N} exactly floor(N/2) modes depart", m0 - mJ == N // 2,
          f"{m0} -> {mJ}")

# ---------- G12: the exceptional couplings are defective ----------

print("G12 exceptional couplings: real, and a Jordan block sits on each")

# Exact throughout. At a root of the cofactor q(0)(J) the algebraic multiplicity
# is >= floor(N/2)+1 (q(0) is the coefficient of eps^floor(N/2), Section 6),
# while the pencil argument still gives only floor(N/2) geometric dimensions.
# A float rank test at a tuned degeneracy is worthless, so: N=3 runs over
# Q(sqrt(3)) (where J* = d1/sqrt(3) is irrational), and N=4 runs at a profile
# whose exceptional coupling is RATIONAL, so both multiplicities are computed
# in plain Gaussian rationals and the algebraic one comes from the
# characteristic polynomial itself, independently of the cofactor theorem.

class QSq:
    """a + b*sqrt(s), with a and b Gaussian rational (pairs of Fractions)."""
    __slots__ = ('a', 'b', 's')

    def __init__(self, a=(0, 0), b=(0, 0), s=1):
        self.a = (Fraction(a[0]), Fraction(a[1]))
        self.b = (Fraction(b[0]), Fraction(b[1]))
        self.s = s

    @staticmethod
    def _m(x, y):
        return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])

    def __add__(u, v):
        return QSq((u.a[0] + v.a[0], u.a[1] + v.a[1]),
                   (u.b[0] + v.b[0], u.b[1] + v.b[1]), u.s)

    def __sub__(u, v):
        return QSq((u.a[0] - v.a[0], u.a[1] - v.a[1]),
                   (u.b[0] - v.b[0], u.b[1] - v.b[1]), u.s)

    def __mul__(u, v):
        bb = QSq._m(u.b, v.b)
        aa = QSq._m(u.a, v.a)
        return QSq((aa[0] + u.s * bb[0], aa[1] + u.s * bb[1]),
                   tuple(x + y for x, y in zip(QSq._m(u.a, v.b),
                                               QSq._m(u.b, v.a))), u.s)

    def inv(u):
        nb = QSq._m(u.b, u.b)
        den = (u.a[0] * u.a[0] - u.a[1] * u.a[1] - u.s * nb[0],
               2 * u.a[0] * u.a[1] - u.s * nb[1])
        d = den[0] * den[0] + den[1] * den[1]
        ic = (den[0] / d, -den[1] / d)
        return QSq(QSq._m(u.a, ic), QSq._m((-u.b[0], -u.b[1]), ic), u.s)

    def is_zero(u):
        return u.a == (0, 0) and u.b == (0, 0)


def xc_build(N, gl, Jcoef, s, sqrt_J=True):
    """Mtilde = L_block + 4 gbar over Q(sqrt(s), i), s a NON-SQUARE.

    The coupling is J = Jcoef*sqrt(s) when sqrt_J is True, and J = Jcoef (a plain
    rational, carried in the a-slot) when it is False. s must not be a perfect
    square: in a + b*sqrt(1) the element a = b has zero norm, the ring is not a
    field, and the elimination below would divide by zero.
    """
    assert s > 0 and round(s ** 0.5) ** 2 != s, "s must be a positive non-square"
    h = se_h_frac(N)
    gbar = sum(gl) / N
    n2 = N * N
    M = [[QSq(s=s) for _ in range(n2)] for _ in range(n2)]

    def add_i(r, c, coeff, in_a):                 # add J*(-i) to the entry
        cur = M[r][c]
        if in_a:
            M[r][c] = QSq((cur.a[0], cur.a[1] - coeff), cur.b, s)
        else:
            M[r][c] = QSq(cur.a, (cur.b[0], cur.b[1] - coeff), s)

    in_a = not sqrt_J
    for a in range(N):
        for b in range(N):
            r = a * N + b
            for c in range(N):
                if h[a][c]:
                    add_i(r, c * N + b, Jcoef * h[a][c], in_a)
                if h[c][b]:
                    add_i(r, a * N + c, -Jcoef * h[c][b], in_a)
            diag = 4 * gbar - (2 * (gl[a] + gl[b]) if a != b else 0)
            cur = M[r][r]
            M[r][r] = QSq((cur.a[0] + diag, cur.a[1]), cur.b, s)
    return M


def xc_matmul(A, B, s):
    n = len(A)
    C = [[QSq(s=s) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for k in range(n):
            if A[i][k].is_zero():
                continue
            a = A[i][k]
            for j in range(n):
                C[i][j] = C[i][j] + a * B[k][j]
    return C


def xc_null(M):
    n = len(M)
    M = [row[:] for row in M]
    rank = 0
    for c in range(n):
        piv = next((i for i in range(rank, n) if not M[i][c].is_zero()), None)
        if piv is None:
            continue
        M[rank], M[piv] = M[piv], M[rank]
        inv = M[rank][c].inv()
        for i in range(rank + 1, n):
            if M[i][c].is_zero():
                continue
            f = M[i][c] * inv
            for j in range(c, n):
                M[i][j] = M[i][j] - f * M[rank][j]
        rank += 1
        if rank == n:
            break
    return n - rank


def xc_jordan(N, gl, Jcoef, s, sqrt_J=True):
    """(nullity Mt, nullity Mt^2, nullity Mt^3): geometric, then the algebraic
    multiplicity once it stabilises, with the increments giving the block sizes.

    The plateau is what licenses reading the last entry as the algebraic
    multiplicity: once two consecutive kernel dimensions agree they agree
    forever, so no deeper block can hide. This is asserted here rather than
    left to the caller, so the routine cannot silently under-report a block of
    size 3 or more if it is reused at a larger N."""
    M = xc_build(N, gl, Jcoef, s, sqrt_J)
    M2 = xc_matmul(M, M, s)
    M3 = xc_matmul(M2, M, s)
    n1, n2, n3 = xc_null(M), xc_null(M2), xc_null(M3)
    assert n2 == n3, f"kernel chain has not stabilised by power 3: {n1, n2, n3}"
    return n1, n2, n3


def xc_nullity(N, gl, Jcoef, s, sqrt_J=True):
    return xc_null(xc_build(N, gl, Jcoef, s, sqrt_J))


def xc_alg_mult(N, gl, J):
    """Order of vanishing of det(eps I - Mtilde) at eps = 0, rational J."""
    h = se_h_frac(N)
    gbar = sum(gl) / N
    n2 = N * N
    base = [[GQ(0) for _ in range(n2)] for _ in range(n2)]
    for a in range(N):
        for b in range(N):
            r = a * N + b
            for c in range(N):
                base[r][c * N + b] += GQ(0, -J * h[a][c])
                base[r][a * N + c] += GQ(0, J * h[c][b])
            if a != b:
                base[r][r] += GQ(-2 * (gl[a] + gl[b]))
            base[r][r] += GQ(4 * gbar)
    nodes = [Fraction(k + 1, 7) for k in range(n2 + 1)]
    vals = []
    for e in nodes:
        A = [[(GQ(e) - base[i][j]) if i == j else GQ(0) - base[i][j]
              for j in range(n2)] for i in range(n2)]
        vals.append(det_exact(A))
    m = len(nodes)

    def solve(comp):
        A = [[Fraction(nodes[i] ** j) for j in range(m)] + [comp[i]]
             for i in range(m)]
        for c in range(m):
            piv = next(r for r in range(c, m) if A[r][c] != 0)
            A[c], A[piv] = A[piv], A[c]
            inv = 1 / A[c][c]
            A[c] = [v * inv for v in A[c]]
            for r in range(m):
                if r != c and A[r][c] != 0:
                    f = A[r][c]
                    A[r] = [v - f * w for v, w in zip(A[r], A[c])]
        return [A[i][m] for i in range(m)]

    re, im = solve([v.a for v in vals]), solve([v.b for v in vals])
    return next(i for i in range(m) if re[i] != 0 or im[i] != 0)


# (a) N=3 over Q(sqrt(3)): at the closed-form exceptional coupling J* = d1/sqrt(3)
#     the kernel dimensions of the powers give the whole Jordan signature, with
#     no appeal to the cofactor theorem: 1, 2, 2 means geometric 1, algebraic 2,
#     one block of size exactly 2.
gb3 = Fraction(9, 100)
for d1 in (Fraction(1, 50), Fraction(1, 20), Fraction(7, 100)):
    gl3 = [gb3 + d1, gb3, gb3 - d1]
    sig = xc_jordan(3, gl3, d1 / 3, 3)          # J* = (d1/3)*sqrt(3)
    check(f"N=3 d1={d1}: at J* the nullities are (1, 2, 2), one 2x2 block",
          sig == (1, 2, 2), f"nullities {sig}")
sig_ctrl = xc_jordan(3, [gb3 + Fraction(1, 50), gb3, gb3 - Fraction(1, 50)],
                     Fraction(1, 7), 3)
check("N=3 control at a non-exceptional coupling: (1, 1, 1), semisimple",
      sig_ctrl == (1, 1, 1), f"nullities {sig_ctrl}")

# (b) N=4 at a profile whose exceptional coupling is RATIONAL: the same signature
#     by the same route, and the algebraic multiplicity confirmed a second time
#     from the characteristic polynomial itself; two neighbouring couplings as
#     controls
sc = Fraction(1, 1000)
gl4 = [Fraction(9, 100) + 1 * sc, Fraction(9, 100) + 3 * sc,
       Fraction(9, 100) - 3 * sc, Fraction(9, 100) - 1 * sc]
sig4 = xc_jordan(4, gl4, Fraction(1, 1000), 2, sqrt_J=False)
alg = xc_alg_mult(4, gl4, Fraction(1, 1000))
check("N=4 at the rational J* = 1/1000: nullities (2, 3, 3), one 2x2 block "
      "beside one 1x1", sig4 == (2, 3, 3), f"nullities {sig4}")
check("N=4 the algebraic multiplicity 3 again from the characteristic "
      "polynomial (two independent routes agree)", alg == sig4[-1],
      f"charpoly {alg}, powers {sig4[-1]}")
for Jc in (Fraction(1, 999), Fraction(1, 1001)):
    s4, a2 = xc_jordan(4, gl4, Jc, 2, sqrt_J=False), xc_alg_mult(4, gl4, Jc)
    check(f"N=4 control J = {Jc}: (2, 2, 2) and charpoly 2, semisimple",
          s4 == (2, 2, 2) and a2 == 2, f"nullities {s4}, charpoly {a2}")

# (c) the counting question of Section 12, exactly: the cofactor is even in
#     z = -i J, so put w = z^2 and strip the J-monomial; a real coupling means
#     w < 0, and Sturm's theorem counts those exactly. The point of the check is
#     the NEGATIVE result: the count is not a function of N, since two generic
#     profiles on the same N = 6 locus disagree.

def st_trim(p):
    while len(p) > 1 and p[-1] == 0:
        p = p[:-1]
    return p


def st_rem(a, b):
    a = a[:]
    while len(a) >= len(b) and any(x != 0 for x in a):
        d, c = len(a) - len(b), a[-1] / b[-1]
        for i in range(len(b)):
            a[d + i] -= c * b[i]
        a = st_trim(a)
        if len(a) < len(b):
            break
    return st_trim(a)


def st_chain(p):
    p = st_trim(p)
    chain = [p, st_trim([p[i] * i for i in range(1, len(p))])]
    while len(chain[-1]) > 1:
        r = st_rem(chain[-2], chain[-1])
        if all(x == 0 for x in r):
            break
        chain.append([-x for x in r])
    return chain


def st_var(chain, at_minus_inf, x=None):
    sg = []
    for p in chain:
        v = (p[-1] * (-1) ** (len(p) - 1) if at_minus_inf
             else sum(c * (x ** i) for i, c in enumerate(p)))
        if v != 0:
            sg.append(1 if v > 0 else -1)
    return sum(1 for i in range(len(sg) - 1) if sg[i] * sg[i + 1] < 0)


def st_real_pairs(N, delta):
    """Distinct real +- pairs of nonzero exceptional couplings, exact."""
    Y, rows, cols, nrm, _ = wc_build(N, True, delta=delta)
    n = len(cols)
    vals = []
    for z in range(2 * n + 1):
        col = [[Y[y][x][0] + z * Y[y][x][1] for y in range(len(rows))]
               for x in range(n)]
        vals.append(wc_det_frac(
            [[sum(Fraction(col[i][y] * col[j][y], nrm[y])
                  for y in range(len(rows))) for j in range(n)]
             for i in range(n)]))
    c = wc_coeffs(vals, 2 * n)
    assert all(c[i] == 0 for i in range(1, len(c), 2)), "cofactor not even in z"
    w = st_trim([c[2 * k] for k in range(len(c) // 2 + 1)])
    m = min(i for i, x in enumerate(w) if x != 0)          # strip the monomial
    ch = st_chain(w[m:])
    return st_var(ch, True) - st_var(ch, False, Fraction(0)), len(w[m:]) - 1


print("G13 the counting question: real exceptional couplings, exact Sturm")
for N, half, want in ((3, [2], 1), (3, [5], 1), (4, [2, -3], 2),
                      (5, [2, 3], 2), (6, [2, -3, 5], 2), (6, [1, 2, 3], 4)):
    delta = [0] * N
    for i, v in enumerate(half[:N // 2]):
        delta[i], delta[N - 1 - i] = v, -v
    pairs, deg = st_real_pairs(N, delta)
    check(f"N={N} delta-half={half}: {want} real pair(s) of {deg} candidates",
          pairs == want, f"counted {pairs}, w-degree {deg}")
check("the count is NOT a function of N: two generic N=6 profiles disagree",
      st_real_pairs(6, [2, -3, 5, -5, 3, -2])[0]
      != st_real_pairs(6, [1, 2, 3, -3, -2, -1])[0])

# ---------- G14: the pointed grading and sharpness (proof doc Section 8.1) ----

print("G14 pointed grading chi_x: admissibility, the census, the staircase")


def pg_chi(N, x, a, b):
    return abs(a - x) + abs(b - (N + 1 - x))


def pg_sets(N):
    """(rows, columns, anti-diagonal rows, diagonal columns) as index lists."""
    _, rows, cols, _, dcols = wc_build(N, True)
    vm, op = [], []
    n2 = N * N
    tq = [(N - 1 - r % N) * N + (N - 1 - r // N) for r in range(n2)]
    for r in range(n2):
        a, b, t = r // N, r % N, tq[r]
        if t != r and r < t:
            vm.append((r, t))
        if a == b:
            continue
        if t == r:
            op.append((r,))
        elif r < t:
            op.append((r, t))
    anti = [y for y in range(len(op)) if len(op[y]) == 1]
    return vm, op, anti, dcols, rows, cols


for N in range(3, 10):
    vm, op, anti, dcols, rows, cols = pg_sets(N)
    ocols = [i for i in range(len(vm)) if i not in dcols]
    bad_adm, bad_F, bad_h = 0, 0, 0
    for y in anti:
        x = op[y][0] // N + 1
        keep = [yy for yy in range(len(op)) if yy not in anti or yy == y]
        # admissibility has TWO halves (Lemma 7); test both, not just the first
        # (a) tauQ-invariance: constant on every parity basis vector
        for cells in op + [v for v in vm]:
            vals = {pg_chi(N, x, r // N + 1, r % N + 1) for r in cells}
            bad_adm += 0 if len(vals) == 1 else 1
        # (b) one hop moves chi_x by at most one
        for a in range(1, N + 1):
            for b in range(1, N + 1):
                base = pg_chi(N, x, a, b)
                for a2, b2 in ((a + 1, b), (a - 1, b), (a, b + 1), (a, b - 1)):
                    if 1 <= a2 <= N and 1 <= b2 <= N:
                        bad_adm += (0 if abs(pg_chi(N, x, a2, b2) - base) <= 1
                                    else 1)
        for c in dcols:
            site = vm[c][0] // N + 1
            h = pg_chi(N, x, site, site)
            bad_h += 0 if h == max(abs(N + 1 - 2 * site),
                                   abs(N + 1 - 2 * x)) else 1
            A = ocols + [c]
            pc = [pg_chi(N, x, vm[i][0] // N + 1, vm[i][0] % N + 1) for i in A]
            pr = [pg_chi(N, x, op[i][0] // N + 1, op[i][0] % N + 1)
                  for i in keep]
            for j in range(1, h + 3):
                F = sum(1 for v in pc if v >= j) - sum(1 for v in pr if v >= j)
                bad_F += 0 if F == (1 if h >= j else 0) else 1
    check(f"N={N} chi_x is admissible: tauQ-invariant AND one-hop Lipschitz",
          bad_adm == 0, f"violations {bad_adm}")
    check(f"N={N} chi_x(c,c) = max(d_c, d_x), the L1 walk length", bad_h == 0,
          f"violations {bad_h}")
    check(f"N={N} the census collapses: F_j = [chi_x(c,c) >= j]", bad_F == 0,
          f"violations {bad_F}")

# the case split of Lemma 10's proof, exhaustively. It stood inverted in the
# document once (a review catch), so it is pinned here: chi_x(c,c) is d_c when
# x lies between c and R(c), and d_x otherwise, and both readings agree with
# max(d_c, d_x).
split_bad = []
for N in range(3, 13):
    for c in range(1, N + 1):
        for x in range(1, N + 1):
            if x == N + 1 - x:
                continue
            chi = abs(c - x) + abs(c - (N + 1 - x))
            dc, dx = abs(N + 1 - 2 * c), abs(N + 1 - 2 * x)
            between = min(c, N + 1 - c) <= x <= max(c, N + 1 - c)
            if chi != (dc if between else dx) or chi != max(dc, dx):
                split_bad.append((N, c, x))
check("Lemma 10 case split: chi_x(c,c) = d_c iff x lies between c and R(c), "
      "N = 3..12", not split_bad, f"deviations {split_bad[:5]}")

# the order law itself, and the staircase it produces, exactly (N = 4..7)
for N in range(4, 8):
    delta = wc_delta(N)
    Y, rows, cols, nrm, dcols = wc_build(N, True)
    vm, op, anti, _, _, _ = pg_sets(N)
    ocols = [i for i in range(len(cols)) if i not in dcols]
    bad_ord, staircase = [], {}
    for y in anti:
        x = op[y][0] // N + 1
        dx = abs(N + 1 - 2 * x)
        keep = [yy for yy in range(len(op)) if yy not in anti or yy == y]
        for c in dcols:
            site = vm[c][0] // N + 1
            dc = abs(N + 1 - 2 * site)
            o = wc_minor_ord(Y, keep, ocols + [c])
            if o != max(dc, dx):
                bad_ord.append((x, site, o, max(dc, dx)))
            staircase.setdefault(dc, set())
            if o == dc:
                staircase[dc].add(x)
    check(f"N={N} ord det Y_(I_x, A_c) = max(d_c, d_x) for every pair and cell",
          not bad_ord, f"deviations {bad_ord}")
    ds = sorted(staircase)
    nested = all(staircase[ds[i]] < staircase[ds[i + 1]]
                 for i in range(len(ds) - 1))
    intervals = all(staircase[d] == {x for x in
                                     [op[y][0] // N + 1 for y in anti]
                                     if abs(N + 1 - 2 * x) <= d}
                    for d in ds)
    check(f"N={N} the supports are the intervals between each pair's own sites",
          intervals, f"supports {[sorted(staircase[d]) for d in ds]}")
    check(f"N={N} and they are STRICTLY nested, so the matrix is triangular",
          nested, f"sizes {[len(staircase[d]) for d in ds]}")

# ---------- verdict ----------

print()
if FAILURES:
    print(f"R90 frozen divisor gate: {len(FAILURES)} FAILURE(S): {FAILURES}")
    sys.exit(1)
print("R90 frozen divisor gate: ALL GREEN")
