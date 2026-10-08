"""F50's exceptional couplings as plane crossings: the Schur bound, the Krein bound, the p = 1 count.

Gate for the section "The count as a plane crossing" of docs/proofs/PROOF_WEIGHT1_DEGENERACY.md.

Objects, in the diagonal joint-popcount block (p, p) of a number-conserving real H under uniform Z-dephasing
(Pauli J, units gamma = 1, x = J/gamma): the block is B(x) = -i x A + D with A = H (x) I - I (x) H^T real symmetric,
D = -2*Hamming. tau is the transpose (|x><y| -> |y><x|). Populations P (Hamming 0) and coherences Q (Hamming >= 2).
Blocks are run for p <= N/2; the X^N mirror carries (p,p) to (N-p,N-p).

Rows (an exact row compares to 0 or to an integer; a float row states its error model and prints the ratio c
to it; a measured row prints what it measured):
  R1  tau B^dagger tau = B exactly, and A_PP = 0 exactly.  Control: a DM term XY - YX breaks the identity in
      this basis (on a tree it is gauge-equivalent to a real H; on the ring it is a flux, see R5).
  R2  The Schur complement of B + 2 onto P is real symmetric and equals 2 I - 2 R(x) with
      R = (x^2/2) A_PQ (Delta + x^2 A_QQ Delta^-1 A_QQ)^-1 A_QP, Delta = 2 Hamming - 2, rank R = C(N,p) - 1.
      Exact over Q (N = 3, 4).  Control: a dissipator that is not tau-even gives Im S != 0.
  R3  R(x2) - R(x1) is positive semidefinite for x1 < x2 with kernel of dimension 1 (exact over Q, N = 3, 4;
      the control with one negated Delta cell runs through the same exact PSD test and must fail).
      At N = 5, 6 in float: the most negative eigenvalue against eps*||R||*cond(T), and the second-smallest
      eigenvalue (the kernel is one-dimensional) above that noise.
  R4  #{eigenvalues of R(x) below 1} = #{eigenvalues of B(x) with Re > -2 by more than the noise} at every grid
      point (eigenvalues within noise of the line are counted and printed: a complex pair running along the
      line at XY N = 3 and on the N = 4 ring; no r_j within noise of 1), the net drop of the first count from the Zeno end to x = 20 equals C(N,p) - 1
      on the chain (Heisenberg N = 3..6, XY N = 3..5), is at most C(N,p) - 1 on XXZ Delta = 1/2, the star, K4
      and the ring (N = 4), and bisected the chain's crossings are E to the six pinned decimals (N = 3, 4).
  R5  kappa(v) = ||v_even||^2 - ||v_odd||^2: every complex eigenvector has kappa = 0 (against eps*||B||/|Im|),
      every eigenvector obeys height >= 1 - kappa, every mode below the plane has kappa >= 1 - height > 0, and the
      smallest complex height (grid plus a local minimisation over x) is >= 1.  Control: the N = 5 ring with a
      DM flux has a complex mode far below the plane; the chain's DM term is a gauge and does not.
  R6  p = 1 Hamiltonian end, the Heisenberg chain: single-excitation H = (N-1) I - 2 L_path (cosine modes),
      W_ab = sum_x E_a^2 E_b^2 = J/N + (I' - F')/(2N) entrywise, spectrum {1, 1/N x floor((N-1)/2), 0 x ceil((N-1)/2)};
      the (1,1) block's real modes end at heights 0, 2(N-1)/N, 2 with the deviation falling as (gamma/J)^2 (ratio 4
      between x = 100 and 200).  The XY chain beside it: W = (J + (I + F)/2)/(N+1), heights 2N/(N+1) and 2.
  R7  Ring (1,1), momentum sector K (rho_xy = e^{iKx} phi(y-x)): the sector matrices' spectra are the block's
      spectrum as a multiset (N = 5); for N = 60 the slow mode is the bound state k = 2 - 2 sqrt(1 - 4 x^2 sin^2(K/2))
      (against eps*||B_K||); at N = 4 the K = pi/2 sector is the N = 2 pair lambda^2 + 4 lambda + 16 x^2 (both
      roots), and the K = pi sector has |det(B_K + 2)|^2 = 256 at thirteen rational x, so exactly, and never meets
      the plane.
  R8  Zeno end: (A^2)_PP = 8 Lap(token graph) exactly; the slow eigenvalues of B are -2 x^2 ell_j + O(x^4)
      (ratio 16 between x = 0.02 and 0.01); the p-token graph's Laplacian carries the path's spectrum and the same
      gap (measured, N = 5..7).

Run: python simulations/f50_plane_crossing_count.py   (about two minutes; prints "ALL GATES PASS" or a FAIL line)
"""
import sys
from fractions import Fraction
from itertools import combinations
from math import comb, sqrt, sin, pi
import numpy as np
from scipy.optimize import minimize_scalar, linear_sum_assignment

sys.stdout.reconfigure(encoding="utf-8")
EPS = np.finfo(float).eps
FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}{(': ' + detail) if detail else ''}")
    if not ok:
        FAILS.append(name)


# ----------------------------------------------------------------------------------------------------------
# builders
# ----------------------------------------------------------------------------------------------------------
def bonds_chain(N):
    return [(i, i + 1) for i in range(N - 1)]


def bonds_ring(N):
    return [(i, (i + 1) % N) for i in range(N)]


def bonds_star(N):
    return [(0, i) for i in range(1, N)]


def bonds_complete(N):
    return [(i, j) for i in range(N) for j in range(i + 1, N)]


def ham_matrix(N, bonds, dm=0.0, delta=1.0):
    """H = sum_bonds XX + YY + delta*ZZ + dm*(XY - YX); integer entries when dm = 0 and delta integer."""
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex if dm else float)
    bit = lambda x, l: (x >> (N - 1 - l)) & 1
    for a, b in bonds:
        for x in range(d):
            if bit(x, a) != bit(x, b):
                y = x ^ (1 << (N - 1 - a)) ^ (1 << (N - 1 - b))
                H[y, x] += 2
                if dm:
                    H[y, x] += 2j * dm * (1 if bit(x, a) else -1)   # XY - YX: Hermitian, imaginary
                H[x, x] -= delta
            else:
                H[x, x] += delta
    return H


def block(N, H, p):
    """A = ad_H in |x><y| coordinates (real symmetric for real H), Hamming vector, index list."""
    d = 2 ** N
    pc = lambda x: bin(x).count('1')
    idx = [(x, y) for x in range(d) for y in range(d) if pc(x) == p and pc(y) == p]
    pos = {e: i for i, e in enumerate(idx)}
    n = len(idx)
    A = np.zeros((n, n), dtype=H.dtype)
    for c, (x, y) in enumerate(idx):
        for xp in range(d):
            if H[xp, x]:
                A[pos[(xp, y)], c] += H[xp, x]
        for yp in range(d):
            if H[y, yp]:
                A[pos[(x, yp)], c] -= H[y, yp]
    ham = np.array([bin(x ^ y).count('1') for x, y in idx])
    return A, ham, idx


def tau_perm(idx):
    pos = {e: i for i, e in enumerate(idx)}
    return np.array([pos[(y, x)] for x, y in idx])


def B_of(A, ham, x):
    return -1j * x * A + np.diag(-2.0 * ham)


def token_laplacian(N, p, bonds):
    confs = [sum(1 << i for i in c) for c in combinations(range(N), p)]
    pos = {c: i for i, c in enumerate(confs)}
    Lap = np.zeros((len(confs), len(confs)))
    for c in confs:
        for a, b in bonds:
            if ((c >> a) & 1) != ((c >> b) & 1):
                d = c ^ (1 << a) ^ (1 << b)
                Lap[pos[c], pos[c]] += 1
                Lap[pos[c], pos[d]] -= 1
    return Lap, confs


# ----------------------------------------------------------------------------------------------------------
# exact rational linear algebra (small blocks)
# ----------------------------------------------------------------------------------------------------------
def frac_matrix(M):
    return [[Fraction(int(round(v))) for v in row] for row in M]


def frac_solve(M, Bm):
    """Solve M X = Bm over Fractions by Gauss-Jordan with pivoting; raises if singular."""
    n = len(M)
    m = len(Bm[0])
    aug = [list(M[i]) + list(Bm[i]) for i in range(n)]
    for col in range(n):
        piv = next((r for r in range(col, n) if aug[r][col] != 0), None)
        if piv is None:
            raise ZeroDivisionError("singular")
        aug[col], aug[piv] = aug[piv], aug[col]
        pv = aug[col][col]
        aug[col] = [v / pv for v in aug[col]]
        for r in range(n):
            if r != col and aug[r][col] != 0:
                f = aug[r][col]
                aug[r] = [a - f * b for a, b in zip(aug[r], aug[col])]
    return [row[n:n + m] for row in aug]


def frac_mul(X, Y):
    return [[sum(X[i][k] * Y[k][j] for k in range(len(Y))) for j in range(len(Y[0]))] for i in range(len(X))]


def frac_det(Mf):
    n = len(Mf)
    M2 = [list(r) for r in Mf]
    d = Fraction(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if M2[r][c] != 0), None)
        if piv is None:
            return Fraction(0)
        if piv != c:
            M2[c], M2[piv] = M2[piv], M2[c]
            d = -d
        d *= M2[c][c]
        for r in range(c + 1, n):
            f = M2[r][c] / M2[c][c]
            if f:
                M2[r] = [a - f * b for a, b in zip(M2[r], M2[c])]
    return d


def frac_psd(S):
    """Exact LDL^T test of a symmetric rational matrix: PSD iff no negative pivot and no zero pivot with a
    nonzero row.  Returns (is_psd, rank)."""
    n = len(S)
    assert all(S[i][j] == S[j][i] for i in range(n) for j in range(n)), "frac_psd needs a symmetric input"
    M = [list(r) for r in S]
    rank = 0
    for k in range(n):
        if M[k][k] < 0:
            return False, rank
        if M[k][k] == 0:
            if any(M[k][j] != 0 for j in range(k + 1, n)):
                return False, rank
            continue
        rank += 1
        pv = M[k][k]
        for i in range(k + 1, n):
            if M[i][k] != 0:
                f = M[i][k] / pv
                for j in range(k, n):
                    M[i][j] -= f * M[k][j]
    return True, rank


def exact_schur_and_R(A, ham, x, delta_override=None):
    """Direct Schur complement of B + 2 onto P (via the realification of Delta + i x A_QQ) and the closed-form R,
    as Fraction matrices on the populations: (S_re, S_im, R).  delta_override replaces the diagonal Delta on Q."""
    P = [i for i, h in enumerate(ham) if h == 0]
    Q = [i for i, h in enumerate(ham) if h > 0]
    Af = frac_matrix(A)
    sub = lambda rows, cols: [[Af[i][j] for j in cols] for i in rows]
    APQ, AQP, AQQ = sub(P, Q), sub(Q, P), sub(Q, Q)
    nq = len(Q)
    dvals = [Fraction(2 * int(ham[q]) - 2) for q in Q] if delta_override is None else delta_override
    Delta = [[dvals[i] if i == j else Fraction(0) for j in range(nq)] for i in range(nq)]
    xf = Fraction(x)
    big = [[Delta[i][j] for j in range(nq)] + [-xf * AQQ[i][j] for j in range(nq)] for i in range(nq)] + \
          [[xf * AQQ[i][j] for j in range(nq)] + [Delta[i][j] for j in range(nq)] for i in range(nq)]
    rhs = [[AQP[i][j] for j in range(len(P))] for i in range(nq)] + [[Fraction(0)] * len(P) for _ in range(nq)]
    sol = frac_solve(big, rhs)                       # (Delta + i x A_QQ)^-1 A_QP as (Re; Im)
    G_re, G_im = sol[:nq], sol[nq:]
    S_re = frac_mul(APQ, G_re)
    S_im = frac_mul(APQ, G_im)
    npop = len(P)
    S_re = [[(Fraction(2) if i == j else Fraction(0)) - xf * xf * S_re[i][j] for j in range(npop)] for i in range(npop)]
    S_im = [[-xf * xf * S_im[i][j] for j in range(npop)] for i in range(npop)]
    Dinv_A = [[AQQ[i][j] / Delta[i][i] for j in range(nq)] for i in range(nq)]
    T = frac_mul(AQQ, Dinv_A)
    T = [[Delta[i][j] + xf * xf * T[i][j] for j in range(nq)] for i in range(nq)]
    R = frac_mul(APQ, frac_solve(T, AQP))
    R = [[xf * xf / 2 * v for v in row] for row in R]
    return S_re, S_im, R


def float_R(A, ham, x, with_cond=False):
    P = ham == 0
    Q = ~P
    APQ, AQP, AQQ = A[np.ix_(P, Q)], A[np.ix_(Q, P)], A[np.ix_(Q, Q)]
    Delta = 2.0 * ham[Q] - 2.0
    T = np.diag(Delta) + x * x * AQQ @ (AQQ / Delta[:, None])
    R = 0.5 * x * x * APQ @ np.linalg.solve(T, AQP)
    R = (R + R.T) / 2
    return (R, np.linalg.cond(T)) if with_cond else R


def n_below(A, ham, x, noise):
    """Eigenvalues above the line by more than the noise, and the number lying within the noise of the line
    (a complex pair running along the line, the trace-midpoint law at N = 2, 3 and on the N = 4 ring)."""
    w = np.linalg.eigvals(B_of(A, ham, x))
    return int(np.sum(w.real > -2.0 + noise)), int(np.sum(np.abs(w.real + 2.0) <= noise))


def n_pos(A, ham, x):
    r = np.linalg.eigvalsh(float_R(A, ham, x))
    return int(np.sum(r < 1.0)), float(np.abs(r - 1.0).min())


# ----------------------------------------------------------------------------------------------------------
# R1: tau-self-adjointness (exact), A_PP = 0 (exact)
# ----------------------------------------------------------------------------------------------------------
print("R1  tau B^dagger tau = B and A_PP = 0, exactly")
for N in (3, 4):
    for p in range(1, N // 2 + 1):
        A, ham, idx = block(N, ham_matrix(N, bonds_chain(N)), p)
        t = tau_perm(idx)
        B = B_of(A, ham, 0.5)
        P = ham == 0
        check(f"N={N} ({p},{p}) residual", np.array_equal(B.conj().T[np.ix_(t, t)], B) and not A[np.ix_(P, P)].any(), "exactly 0")
A, ham, idx = block(3, ham_matrix(3, bonds_chain(3), dm=0.5), 1)
t = tau_perm(idx)
B = B_of(A, ham, 0.5)
check("control: DM term XY - YX (dm=0.5) breaks tau-self-adjointness in this basis", not np.array_equal(B.conj().T[np.ix_(t, t)], B),
      f"max residual {np.abs(B.conj().T[np.ix_(t, t)] - B).max():.3f}")

# ----------------------------------------------------------------------------------------------------------
# R2 / R3: the Schur complement exactly, Loewner monotonicity exactly (N <= 4), float at N = 5, 6
# ----------------------------------------------------------------------------------------------------------
print("R2  Schur complement onto populations = 2I - 2R, R real symmetric, rank C(N,p) - 1 (exact over Q)")
print("R3  R(x2) - R(x1) >= 0 with one-dimensional kernel (exact over Q)")
for N in (3, 4):
    H = ham_matrix(N, bonds_chain(N))
    for p in range(1, N // 2 + 1):
        A, ham, idx = block(N, H, p)
        Rs = {}
        for x in (Fraction(1, 3), Fraction(1, 2), Fraction(2), Fraction(5)):
            S_re, S_im, R = exact_schur_and_R(A, ham, x)
            n = len(R)
            im_zero = all(v == 0 for row in S_im for v in row)
            sym = all(S_re[i][j] == S_re[j][i] for i in range(n) for j in range(n))
            eq = all(S_re[i][j] == (2 if i == j else 0) - 2 * R[i][j] for i in range(n) for j in range(n))
            psd, rank = frac_psd(R)
            check(f"N={N} ({p},{p}) x={x}: Im S = 0, S symmetric, S = 2I - 2R, R >= 0, rank R = C(N,p) - 1",
                  im_zero and sym and eq and psd and rank == comb(N, p) - 1, f"rank {rank}")
            Rs[x] = R
        xs = sorted(Rs)
        for a, b in zip(xs, xs[1:]):
            diff = [[Rs[b][i][j] - Rs[a][i][j] for j in range(n)] for i in range(n)]
            psd, rank = frac_psd(diff)
            check(f"N={N} ({p},{p}) R({b}) - R({a}) >= 0, kernel dim 1", psd and rank == n - 1, f"rank {rank} of {n}")
# controls, through the same exact routines
A, ham, idx = block(3, ham_matrix(3, bonds_chain(3)), 1)
Q = [i for i, h in enumerate(ham) if h > 0]
dv = [Fraction(2 * int(ham[q]) - 2) for q in Q]
dv_odd = list(dv)
dv_odd[0] += 1                                        # one cell of a tau-pair dephased differently: not tau-even
S_re, S_im, R = exact_schur_and_R(A, ham, Fraction(1, 2), delta_override=dv_odd)
check("control: a dissipator that is not tau-even gives Im S != 0", any(v != 0 for row in S_im for v in row))
dv_neg = list(dv)
dv_neg[0] = -dv_neg[0]
xs_c = [Fraction(1, 4), Fraction(1, 2), Fraction(1), Fraction(2), Fraction(4)]
Rm = {x: exact_schur_and_R(A, ham, x, delta_override=dv_neg)[2] for x in xs_c}
broken = [not frac_psd([[Rm[b][i][j] - Rm[a][i][j] for j in range(3)] for i in range(3)])[0] for a, b in zip(xs_c, xs_c[1:])]
check("control: one negated Delta cell breaks Loewner monotonicity (the exact PSD test fails on some step)", any(broken),
      f"failing steps {sum(broken)} of {len(broken)}")

print("R3f Loewner at N = 5, 6 in float: min eig(R(x2) - R(x1)) against eps*||R||*cond(T); second-smallest above it")
for N, ps in ((5, (1, 2)), (6, (1, 2, 3))):
    H = ham_matrix(N, bonds_chain(N))
    for p in ps:
        A, ham, idx = block(N, H, p)
        xs = [0.1, 0.3, 0.6, 1.0, 1.5, 2.5, 4.0, 8.0]
        Rs = [float_R(A, ham, x, with_cond=True) for x in xs]
        worst_neg, worst_second = 0.0, np.inf
        for (Ra, ca), (Rb, cb) in zip(Rs, Rs[1:]):
            ev = np.sort(np.linalg.eigvalsh(Rb - Ra))
            noise = EPS * np.abs(Rb).max() * max(ca, cb)
            worst_neg = max(worst_neg, -ev[0] / noise)
            worst_second = min(worst_second, ev[1] / noise)
        check(f"N={N} ({p},{p}) min eig >= -c*noise with c = {worst_neg:.3f}; second eigenvalue / noise >= {worst_second:.1e}",
              worst_neg < 10 and worst_second > 1e3)

# ----------------------------------------------------------------------------------------------------------
# R4: the count, n_pos(S) = n_below(B) on a grid; net drops; crossings = E
# ----------------------------------------------------------------------------------------------------------
print("R4  #{r_j < 1} = #{Re lambda > -2} on the grid; net drop from the Zeno end to x = 20; crossings = E")
E_known = {3: [1.249621, 1.732051],
           4: [0.745022, 0.745439, 0.910180, 1.545936, 1.572303, 1.867978, 2.088800, 2.197368]}


def count_rows(label, N, H, ps, expect_equal, points=400):
    found = []
    for p in ps:
        A, ham, idx = block(N, H, p)
        grid = np.geomspace(0.05, 20.0, points)
        normB = 2 * ham.max() + 20.0 * np.abs(A).sum(axis=0).max()
        mism, min_margin, on_line = 0, np.inf, 0
        noise = 1e3 * EPS * normB
        for x in grid:
            a, ma = n_pos(A, ham, x)
            b, nb = n_below(A, ham, x, noise)
            mism += (a != b)
            on_line += nb
            min_margin = min(min_margin, ma / EPS)
        check(f"{label} N={N} ({p},{p}) n_pos(S) = n_below(B) at every grid point; nearest approach of an r_j to 1: {min_margin:.1e} eps; "
              f"eigenvalues within noise of the line: {on_line}", mism == 0 and min_margin > 1e3, f"{mism} mismatches of {len(grid)}")
        drops = n_pos(A, ham, grid[0])[0] - n_pos(A, ham, grid[-1])[0]
        ok = drops == comb(N, p) - 1 if expect_equal else drops <= comb(N, p) - 1
        check(f"{label} N={N} ({p},{p}) net drop {'=' if expect_equal else '<='} C(N,p) - 1 = {comb(N, p) - 1}", ok,
              f"{drops}; n_pos {n_pos(A, ham, grid[0])[0]} at the Zeno end, {n_pos(A, ham, grid[-1])[0]} at x = 20")
        if label == "Heisenberg chain" and N <= 4:
            def locate(a, b, ca, cb):
                if ca == cb:
                    return []
                if ca - cb == 1 or b / a < 1 + 1e-13:
                    lo, hi = a, b
                    for _ in range(60):
                        m = sqrt(lo * hi)
                        if n_pos(A, ham, m)[0] == ca:
                            lo = m
                        else:
                            hi = m
                    return [1.0 / lo] * (ca - cb)
                m = sqrt(a * b)
                cm = n_pos(A, ham, m)[0]
                return locate(a, m, ca, cm) + locate(m, b, cm, cb)
            counts = [(x, n_pos(A, ham, x)[0]) for x in grid]
            for (xa, ca), (xb, cb) in zip(counts, counts[1:]):
                found += locate(xa, xb, ca, cb)
    if found:
        found.sort()
        ok = len(found) == len(E_known[N]) and all(abs(f - e) <= 5e-7 for f, e in zip(found, E_known[N]))
        check(f"{label} N={N} the crossings are E(N, chain) to the six pinned decimals", ok, ", ".join(f"{f:.6f}" for f in found))


for N in (3, 4, 5, 6):
    count_rows("Heisenberg chain", N, ham_matrix(N, bonds_chain(N)), range(1, N // 2 + 1), True, 400 if N < 6 else 150)
for N in (3, 4, 5):
    count_rows("XY chain", N, ham_matrix(N, bonds_chain(N), delta=0.0), range(1, N // 2 + 1), True, 200)
count_rows("XXZ chain Delta=1/2", 4, ham_matrix(4, bonds_chain(4), delta=0.5), (1, 2), False, 200)
count_rows("star", 4, ham_matrix(4, bonds_star(4)), (1, 2), False, 200)
count_rows("K4", 4, ham_matrix(4, bonds_complete(4)), (1, 2), False, 200)
count_rows("ring", 4, ham_matrix(4, bonds_ring(4)), (1, 2), False, 200)

# ----------------------------------------------------------------------------------------------------------
# R5: the Krein signature
# ----------------------------------------------------------------------------------------------------------
print("R5  kappa = ||v_even||^2 - ||v_odd||^2: complex modes neutral and at height >= 1; height >= 1 - kappa")


def krein_rows(N, bonds, dm=0.0):
    H = ham_matrix(N, bonds, dm=dm)
    out = dict(min_complex=np.inf, kappa_c=0.0, bound=np.inf, below=np.inf, noise=0.0)
    for p in range(1, N // 2 + 1):
        A, ham, idx = block(N, H, p)
        t = tau_perm(idx)
        normB = 2 * ham.max() + 20.0 * np.abs(A).sum(axis=0).max()
        split = 10 * sqrt(EPS * normB)                       # an EP's float splitting is ~ sqrt(eps ||B||)
        out['noise'] = max(out['noise'], EPS * normB)

        def complex_height(x):
            w = np.linalg.eigvals(B_of(A, ham, x))
            c = w[np.abs(w.imag) > split]
            return (-c.real / 2).min() if len(c) else 9.0

        grid = np.geomspace(0.05, 20.0, 120)
        for x in grid:
            w, V = np.linalg.eig(B_of(A, ham, x))
            for j in range(len(w)):
                v = V[:, j] / np.linalg.norm(V[:, j])
                kappa = float(np.real(np.vdot(v, v[t])))
                height = -w[j].real / 2
                out['bound'] = min(out['bound'], height - (1 - kappa))
                if abs(w[j].imag) > split:
                    out['kappa_c'] = max(out['kappa_c'], abs(kappa) * abs(w[j].imag) / (EPS * normB))
                elif height < 1:
                    out['below'] = min(out['below'], kappa - (1 - height))
        vals = [complex_height(x) for x in grid]
        i = int(np.argmin(vals))
        lo, hi = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
        r = minimize_scalar(complex_height, bounds=(lo, hi), method='bounded', options=dict(xatol=1e-10))
        out['min_complex'] = min(out['min_complex'], vals[i], r.fun)
    return out


for N in (3, 4, 5):
    o = krein_rows(N, bonds_chain(N))
    check(f"N={N} chain: height >= 1 - kappa for every eigenvector and kappa >= 1 - height below the plane, slack against eps*||B||: c = {max(-o['bound'], -o['below']) / o['noise']:.2f}",
          o['bound'] > -100 * o['noise'] and o['below'] > -100 * o['noise'])
    check(f"N={N} chain: complex modes neutral, |kappa|*|Im|/(eps*||B||) <= c, c = {o['kappa_c']:.2f}", o['kappa_c'] < 100)
    check(f"N={N} chain: smallest complex height (grid + local minimisation) >= 1 (against eps*||B||)", o['min_complex'] >= 1 - 100 * o['noise'], f"{o['min_complex']:.6f}")
for N in (4, 5):
    o = krein_rows(N, bonds_ring(N))
    check(f"N={N} ring: complex modes neutral (c = {o['kappa_c']:.2f}) and smallest complex height >= 1 (against eps*||B||)",
          o['kappa_c'] < 100 and o['min_complex'] >= 1 - 100 * o['noise'], f"{o['min_complex']:.9f}")
o = krein_rows(5, bonds_ring(5), dm=0.3)
check("control: N=5 ring with a DM flux (dm=0.3) has a complex mode below the plane", o['min_complex'] < 1,
      f"smallest complex height {o['min_complex']:.4f}")
o = krein_rows(4, bonds_chain(4), dm=0.7)
check("the chain's DM term is a gauge: complex modes keep height >= 1 although kappa != 0 in this basis",
      o['min_complex'] >= 1 and o['kappa_c'] > 1e6, f"smallest complex height {o['min_complex']:.4f}")

# ----------------------------------------------------------------------------------------------------------
# R6: p = 1 Hamiltonian end, the chain's W and the limiting heights
# ----------------------------------------------------------------------------------------------------------
print("R6  p = 1 Hamiltonian end: W = sum_x E_a^2 E_b^2 in closed form, so m_1 = 0 at every N")
print("    Heisenberg chain: single-excitation H = (N-1) I - 2 L_path, cosine modes, W = J/N + (I' - F')/(2N)")
cs = {}
for N in range(2, 61):
    adj = np.zeros((N, N))
    for i in range(N - 1):
        adj[i, i + 1] = adj[i + 1, i] = 1.0
    ends = np.zeros((N, N))
    ends[0, 0] = ends[-1, -1] = 1.0
    Hse = 2 * adj + 2 * ends + (N - 5) * np.eye(N)
    Lpath = np.diag(adj.sum(axis=1)) - adj
    assert np.array_equal(Hse, (N - 1) * np.eye(N) - 2 * Lpath)
    E_true = np.array([[(1 / sqrt(N)) if m == 0 else sqrt(2 / N) * np.cos(pi * m * (x + 0.5) / N) for m in range(N)] for x in range(N)])
    res_eig = np.abs(Hse @ E_true - E_true @ np.diag((N - 1) - 2 * (2 - 2 * np.cos(pi * np.arange(N) / N)))).max()
    W = (E_true ** 2).T @ (E_true ** 2)
    closed = np.full((N, N), 1.0 / N)
    for a in range(1, N):
        closed[a, a] += 1 / (2 * N)
        closed[a, N - a] -= 1 / (2 * N)
    res_W = np.abs(W - closed).max()
    mu = np.sort(np.linalg.eigvalsh(closed))[::-1]
    target = [1.0] + [1.0 / N] * ((N - 1) // 2) + [0.0] * (N // 2)
    cs[N] = (res_eig / (EPS * N), res_W / EPS, np.abs(mu - target).max() / EPS)
worst = max(max(v) for v in cs.values())
check("N = 2..60 Heisenberg chain: cosine modes are eigenvectors (c against eps*N), W equals the closed form entrywise and has spectrum {1, 1/N, 0} (c against eps)",
      worst < 50, "c at N = 10, 20, 40, 60: " + "; ".join(", ".join(f"{v:.1f}" for v in cs[n]) for n in (10, 20, 40, 60)))
worst_xy = 0.0
for N in range(2, 61):
    adj = np.zeros((N, N))
    for i in range(N - 1):
        adj[i, i + 1] = adj[i + 1, i] = 1.0
    E_true = np.array([[sqrt(2 / (N + 1)) * np.sin(pi * (m + 1) * (x + 1) / (N + 1)) for m in range(N)] for x in range(N)])
    worst_xy = max(worst_xy, np.abs(adj @ E_true - E_true @ np.diag(2 * np.cos(pi * np.arange(1, N + 1) / (N + 1)))).max() / (EPS * N))
    W = (E_true ** 2).T @ (E_true ** 2)
    F = np.eye(N)[::-1]
    worst_xy = max(worst_xy, np.abs(W - (np.ones((N, N)) + 0.5 * (np.eye(N) + F)) / (N + 1)).max() / EPS)
check("N = 2..60 XY chain: sine modes, W = (J + (I + F)/2)/(N + 1) entrywise (c against eps)", worst_xy < 50, f"c = {worst_xy:.1f}")
for N in (3, 4, 5, 6):
    A, ham, idx = block(N, ham_matrix(N, bonds_chain(N)), 1)
    heights = {}
    for x in (100.0, 200.0):
        w = np.linalg.eigvals(B_of(A, ham, x))
        heights[x] = np.sort(-w.real[np.abs(w.imag) < 1e-6] / 2)
    target = np.sort([0.0] + [2.0 * (N - 1) / N] * ((N - 1) // 2) + [2.0] * (N // 2))
    d100 = np.abs(heights[100.0] - target).max()
    d200 = np.abs(heights[200.0] - target).max()
    check(f"N={N} (1,1) real modes at the Hamiltonian end -> heights 0, 2(N-1)/N, 2; deviation ratio x=100/x=200 = 4 (1/x^2 law)",
          len(heights[200.0]) == N and abs(d100 / d200 - 4) < 0.04, f"deviation {d100:.2e}, {d200:.2e}, ratio {d100 / d200:.3f}")

# ----------------------------------------------------------------------------------------------------------
# R7: ring (1,1), momentum sectors
# ----------------------------------------------------------------------------------------------------------
print("R7  ring (1,1) sector K: bound state k = 2 - 2 sqrt(1 - 4 x^2 sin^2(K/2)) (N -> infinity); N = 4 sectors")


def sector_matrix(N, K, x):
    """B_K on phi(r), r in Z_N, from rho_xy = e^{iKx} phi(y - x):
    -2 i x [(e^{iK} - 1) phi(r-1) + (e^{-iK} - 1) phi(r+1)] - 4 (1 - delta_r0) phi(r); hopping modulus 4 x sin(K/2)."""
    M = np.zeros((N, N), dtype=complex)
    for r in range(N):
        M[r, (r - 1) % N] += -2j * x * (np.exp(1j * K) - 1)
        M[r, (r + 1) % N] += -2j * x * (np.exp(-1j * K) - 1)
        if r:
            M[r, r] = -4.0
    return M


A, ham, idx = block(5, ham_matrix(5, bonds_ring(5)), 1)
for x in (0.3, 1.0):
    full = np.linalg.eigvals(B_of(A, ham, x))
    sect = np.concatenate([np.linalg.eigvals(sector_matrix(5, 2 * pi * m / 5, x)) for m in range(5)])
    r, c = linear_sum_assignment(np.abs(full[:, None] - sect[None, :]))
    dist = np.abs(full[r] - sect[c]).max() / (EPS * (4 + 4 * x))
    check(f"N=5 ring (1,1) at x={x}: the sector matrices' spectra are the block spectrum as a multiset (c = {dist:.1f} against eps*||B||)", dist < 100)
N = 60
worst = 0.0
for m in (1, 2, 5, 15, 30):
    K = 2 * pi * m / N
    for x in (0.2, 0.5, 0.8):
        if 4 * x * x * sin(K / 2) ** 2 >= 0.99:
            continue
        w = np.linalg.eigvals(sector_matrix(N, K, x))
        pred = -2 * (2 - 2 * sqrt(1 - 4 * x * x * sin(K / 2) ** 2))
        worst = max(worst, abs(max(w.real) - pred) / (EPS * (4 + 8 * x * sin(K / 2))))
check(f"N=60 ring: slow mode of every tested (K, x) on the bound-state curve (c = {worst:.1f} against eps*||B_K||)", worst < 100)
worst = 0.0
for x in (0.1, 0.3, 0.45, 0.6, 1.0):
    w = np.linalg.eigvals(sector_matrix(4, pi / 2, x))
    pair = np.roots([1, 4, 16 * x * x])
    worst = max(worst, max(np.abs(w - pr).min() for pr in pair) / (EPS * (4 + 8 * x)))
check(f"N=4 ring, K=pi/2: both roots of lambda^2 + 4 lambda + 16 x^2 are in the sector (c = {worst:.1f}; the N = 2 pair as a factor, EP on the plane at gamma/J = 2)", worst < 100)


def det_Bpi_plus_2_abs2(x):
    """|det(B_pi + 2)|^2 at N = 4 via the realification [[Re, -Im],[Im, Re]], whose determinant is |det|^2;
    B_pi + 2 = diag(2, -2, -2, -2) + 4 i x (S + S^-1) on Z_4."""
    Re = [[2, 0, 0, 0], [0, -2, 0, 0], [0, 0, -2, 0], [0, 0, 0, -2]]
    Im = [[0, 4 * x, 0, 4 * x], [4 * x, 0, 4 * x, 0], [0, 4 * x, 0, 4 * x], [4 * x, 0, 4 * x, 0]]
    big = [[Fraction(Re[i][j]) for j in range(4)] + [Fraction(-Im[i][j]) for j in range(4)] for i in range(4)] + \
          [[Fraction(Im[i][j]) for j in range(4)] + [Fraction(Re[i][j]) for j in range(4)] for i in range(4)]
    return frac_det(big)


vals = [det_Bpi_plus_2_abs2(Fraction(k, 7)) for k in range(0, 39, 3)]
check("N=4 ring, K=pi: |det(B_K + 2)|^2 = 256 at thirteen rational x (a polynomial of degree <= 8 in x: the constant 16^2, no root; the sector never meets the plane)",
      len(vals) == 13 and all(v == 256 for v in vals))
hs = [-max(np.linalg.eigvals(sector_matrix(4, pi, x)).real) / 2 for x in (20.0, 200.0)]
check(f"N=4 ring, K=pi: the slow mode approaches the plane from below as 1/x^2 ((1-h(20))/(1-h(200)) = {(1 - hs[0]) / (1 - hs[1]):.1f})",
      hs[0] < hs[1] < 1 and abs((1 - hs[0]) / (1 - hs[1]) - 100) < 1)

# ----------------------------------------------------------------------------------------------------------
# R8: Zeno end, the token-graph Laplacian
# ----------------------------------------------------------------------------------------------------------
print("R8  Zeno end: (A^2)_PP = 8 Lap exactly; slow eigenvalues of B = -2 x^2 ell_j + O(x^4); token spectrum carries the path's")
for N, p in ((4, 1), (4, 2), (5, 2)):
    A, ham, idx = block(N, ham_matrix(N, bonds_chain(N)), p)
    P = [i for i, h in enumerate(ham) if h == 0]
    Lap, confs = token_laplacian(N, p, bonds_chain(N))
    # align the token order with the block's population order (the block numbers site l as bit N-1-l)
    pops = [idx[i][0] for i in P]
    conv = {c: sum(1 << (N - 1 - l) for l in range(N) if (c >> l) & 1) for c in confs}
    order = [pops.index(conv[c]) for c in confs]
    A2 = (A @ A)[np.ix_(P, P)][np.ix_(order, order)]
    check(f"N={N} ({p},{p}) (A^2)_PP = 8 Lap(token graph) exactly", np.array_equal(A2, 8 * Lap))
    ell = np.sort(np.linalg.eigvalsh(Lap))
    devs = {}
    for x in (0.01, 0.02):
        w = np.sort(np.linalg.eigvals(B_of(A, ham, x)).real)[::-1][:len(P)]
        devs[x] = np.abs(w - (-2 * x * x * ell)).max()
    ratio = devs[0.02] / devs[0.01]
    check(f"N={N} ({p},{p}) slow eigenvalues = -2 x^2 ell_j with deviation ratio x=0.02/x=0.01 = 16 (x^4 law)", abs(ratio - 16) < 0.16,
          f"deviation {devs[0.01]:.2e}, {devs[0.02]:.2e}, ratio {ratio:.3f}")
for N, p in ((5, 2), (6, 2), (6, 3), (7, 3)):
    Lap, _ = token_laplacian(N, p, bonds_chain(N))
    ev = np.sort(np.linalg.eigvalsh(Lap))
    path = 2 - 2 * np.cos(pi * np.arange(N) / N)
    c = max(np.abs(ev - e).min() for e in path) / (EPS * 4)
    check(f"N={N} p={p}: the token graph's Laplacian carries the path Laplacian's spectrum (c = {c:.1f}) and the same gap {ev[1]:.6f} (measured)",
          c < 100 and abs(ev[1] - path[1]) < 100 * EPS * 4)

print()
if FAILS:
    print("FAILED:", *FAILS, sep="\n  ")
    sys.exit(1)
print("ALL GATES PASS")
