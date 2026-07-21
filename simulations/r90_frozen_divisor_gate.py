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
#       for J in {0.6, 1.0, 2.3} (N = 3..6); control block (2,2) carries none
#   G3  the pencil kernel dim(ker of (E', K_{D-})) on O+ equals floor(N/2)
#   G4  eigenvector predictions: v_D = 0 and tauQ-even O-part (machine zero)
#   G5  partial balance yields nothing (N = 5 with the middle off the common
#       mean; N = 6 with two of three pair sums equal)
#   G6  N=3 closed form (sympy, exact): det(-2g2 I - M_(1,2)) =
#       512 J^4 (g1+g3)^2 (g1+g3-2g2) (4J^2 + (g1+g3)(g1+g3-2g2))
#   G7  N=6 exact arithmetic (Gaussian rationals): on-locus 36x36 corner det
#       is exactly zero; transverse vanishing order is 3
#
# Runtime: ~1-2 min. Standalone except G0 (imports framework once).
import sys
import math
import random
from fractions import Fraction

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

print("G2  corner-block multiplicity floor(N/2) at -4*gbar; (2,2) control empty")
for N in (3, 4, 5, 6):
    gl = r90_profile(N)
    gbar = sum(gl) / N
    K, G = corner_pieces(N, gl)
    mmin, counts = frozen_count_at(lambda J: J * K - 2 * G, -4 * gbar)
    check(f"N={N} corner multiplicity = {N // 2}", mmin == N // 2, f"counts {counts}")
for N in (4, 5, 6):
    gl = r90_profile(N)
    gbar = sum(gl) / N
    mmin, counts = frozen_count_at(lambda J: general_block(N, gl, J, 2, 2), -4 * gbar)
    check(f"N={N} (2,2) control carries none", mmin == 0, f"counts {counts}")
# G2b: XY variant (h without the ZZ diagonal): same multiplicity
for N in (4, 5):
    gl = r90_profile(N)
    gbar = sum(gl) / N
    Kxy, Gxy = corner_pieces(N, gl, zz=False)
    mmin, counts = frozen_count_at(lambda J: J * Kxy - 2 * Gxy, -4 * gbar)
    check(f"N={N} XY corner multiplicity = {N // 2}", mmin == N // 2, f"counts {counts}")
# G2c: the GammaFold corollary: the antidiagonal corner (1, N-1) freezes at
# 4*gbar - 2*sigma with the same multiplicity
for N in (4, 5):
    gl = r90_profile(N)
    gbar = sum(gl) / N
    sigma = sum(gl)
    mmin, counts = frozen_count_at(
        lambda J: general_block(N, gl, J, 1, N - 1), 4 * gbar - 2 * sigma)
    check(f"N={N} (1,{N-1}) corner at 4*gbar-2*sigma, mult {N // 2}",
          mmin == N // 2, f"counts {counts}")

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

# ---------- verdict ----------

print()
if FAILURES:
    print(f"R90 frozen divisor gate: {len(FAILURES)} FAILURE(S): {FAILURES}")
    sys.exit(1)
print("R90 frozen divisor gate: ALL GREEN")
