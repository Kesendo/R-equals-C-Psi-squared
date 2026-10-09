"""F50's Zeno end: the population generator is the ferromagnet, and on every tree its gap is lambda_1(T).

Gate for the Zeno-end paragraphs of docs/proofs/PROOF_WEIGHT1_DEGENERACY.md, section "The count as a plane
crossing": "The Zeno end and the dispersion", Theorem E with its proof ("By the lift the token gap is at most"), and the Zeno half of "The Zeno end and the Hamiltonian end side by side" (its overlaps are read by
f50_isotropy_and_energy_mode.py, rows E3 and E4).

Objects (Pauli book, H = sum_b J_b (X_a X_b + Y_a Y_b + Delta Z_a Z_b), uniform Z-dephasing, x = J/gamma):
in the diagonal joint-popcount block (p, p) the populations P see, as gamma -> infinity, the generator
-(A^2)_PP / (4 gamma) with A = H (x) I - I (x) H^T; the exclusion graph of block p has the popcount-p
configurations as vertices and an edge of weight c_b for every bond b whose two sites disagree.
Site l is bit N-1-l throughout, as in f50_plane_crossing_count.py (whose helpers are copied, not imported:
importing a gate runs it).

Rows (an exact row compares integers, rationals or Gaussian rationals with ==; a float row states its error model and prints
the ratio to it; a reading prints and gates nothing):
  Z1  The exclusion Laplacian is the ferromagnet: 2 Lap(block p) == block p of sum_b c_b (I - X X - Y Y - Z Z)
      built from Kronecker products, integer weights, on the chain, ring, star, complete graph and random graphs,
      every p.  Control: the same operator built without its Y Y term goes red.
  Z2  (A^2)_PP is 8 times that Laplacian with the bond weights |t_b|^2/4, t_b the bond's hop amplitude, so the Zeno
      generator -(A^2)_PP/(4 gamma) is -(2/gamma) times it:
      (A^2)_PP == 8 Lap(J_b^2) exactly at Delta = 0, 1/2, 1, 2, -1 on chains with integer bond profiles and
      Delta = 3/2 on a random graph; unchanged by terms diagonal in the configurations (the longitudinal
      fields (2, -1, 3, 0, -2) and next-nearest ZZ couplings); with a Dzyaloshinskii-Moriya term D_b the weight is J_b^2 + D_b^2.
      The rows in Delta, fields and next-nearest ZZ are structural (the diagonal of H cancels in A_QP, so they
      check the construction); the Liouvillian rows of Z6 test the blindness.  The assembled generator commutes
      with S^- exactly on three of these rows.  Controls: the weights J_b in
      place of J_b^2 miss the first identity; a correlated hop, its amplitude doubled when a third site is
      occupied, breaks the commutation (SU(2) gone).
  Z2b Under a rational profile of rates gamma_l (real H: Delta = 1/2 with fields and a next-nearest ZZ on a chain;
      Delta = 3/2 with fields on a graph with a triangle, since on a graph without triangles the next term vanishes for every
      H: its diagonal parts cancel pair by pair and the rest needs three configurations pairwise one hop apart) minus the leading term of the Schur complement is exact: A_PQ Gamma_Q^-1 A_QP == Lap with the
      weights 4 J_b^2/(gamma_i + gamma_j), Gamma the cells' rates 2 sum gamma over the disagreeing sites, and the
      next term A_PQ Gamma^-1 A_QQ Gamma^-1 A_QP == 0 exactly (the transpose makes every odd order vanish for a
      real H).  Controls: the uniform-rate weights 2 J_b^2/gamma_mean miss the first identity; a
      Dzyaloshinskii-Moriya term on a bond of the triangle (a flux, H no longer real) makes the next term nonzero,
      computed exactly in Gaussian rationals; that term is purely imaginary and antisymmetric there, so i M3 moves
      no real part (row Z6c reads what it does to a degenerate level); and on a chain and a 4-ring that carry a
      Dzyaloshinskii-Moriya term and no triangle the next term vanishes, exactly, H not real.
  Z3  For the uniform Heisenberg graph Lap(block p) == ((#bonds) I - H_p)/2 exactly: the Hamiltonian read from
      the top of its spectrum.  Control: Delta = 2 goes red.
  Z4  The lift is SU(2)'s lowering: (S^-)^(p-1) sum_i f(i)|i> == (p-1)! sum_S F(S)|S>, F(S) = sum_{i in S} f(i),
      and [sum_b c_b (1 - SWAP_b), S^-] == 0, and Lap F == lift of (L_site f), all exact on integer f.
      Controls: S^- with one site missing breaks the lift; the anisotropic sum_b c_b (I - XX - YY - 2 ZZ) breaks
      the commutation with S^-; a site Laplacian with one diagonal entry raised breaks Lap F == lift(L_site f).
  Z5  The theorem, without an eigensolver: for every tree (paths N = 2..9, stars N = 3..9, random integer-weighted
      trees N = 3..9, two chains with random integer bond profiles) and every block 1 <= p <= N-1, the number of
      eigenvalues below a rational q is read exactly (Jacobi's rule on Bareiss minors of b*Lap - a*I over the
      integers; an eigensolver only places the rationals, which the counts then verify).  Below q_lo < lambda_1(T)
      only the zero mode; below q_hi in (lambda_1, lambda_1 + 1e-6 max(1, lambda_1)) exactly as
      many as the tree's Laplacian has; below q* placed 1 to 2 ppm under max_leaf lambda_1(T - v), again exactly the tree's count
      (the new multiplets lie at or above lambda_1(T - v)).  Each bracket is verified exactly on the site Laplacians.
      On the stars with N >= 4 max_leaf lambda_1(T - v) = lambda_1(T), so q* sits below lambda_1 there and adds nothing to the
      q_lo bracket; the row's label says so.
      Controls, through the same row, each counting only with exact brackets and a definite count that differs:
      the blocks of the N = 7 chain with its bond (3, 4) doubled, and with it halved, counted against the
      unchanged chain's brackets, miss (through q_hi and through q_lo); and, on what the theorem adds, the term
      -c (1 - SWAP_01)(1 - SWAP_23), which vanishes on the ferromagnetic multiplet and the one-magnon band and
      moves only multiplets of lower spin, added to every block of the N = 6 chain: at c = 1 a new multiplet
      falls below lambda_1 (every bracket misses), at bond weights 100 and c = 95 one lands between lambda_1 and
      lambda_1(T - v) (q* alone misses).  A cross-check compares the exact counter with a float count far from
      the eigenvalues.
  Z6  The Zeno limit of the full Liouvillian is blind to Delta: on the chains N = 4, 5 the C(N,p) slowest
      eigenvalues of B(x) are -2 x^2 ell_j + O(x^4) with ell_j the Laplacian's spectrum at Delta = 0, 1/2, 1, 2,
      and at Delta = 1 with the fields (2, -1, 3, 0) and a next-nearest ZZ on N = 4.  The error model is the
      expansion in even powers (for a real H, as on every row here, every odd order vanishes by the transpose of
      Lemma 1): deviation/x^4 = a + b x^2 + ...
      A correction in x^k makes the successive differences over x = 0.005, 0.01, 0.02 stand in the ratio 2^k;
      the law is k = 2, and the gate decides it against k = 1 and k = 3 at the half-integer powers,
      2^(3/2) < ratio < 2^(5/2) (measured within 3 % of 4, the float route's rounding included).  A leading generator off by a relative
      delta adds 2 delta ell/x^2, k = -2.  Control on every row: ell (1 + 1e-7), fed through the same
      comparison, is not read as the law.  Under a random profile of rates gamma_l and couplings J_b at
      Delta = 7/10 (chain N = 5, every block) the slow eigenvalues are minus the spectrum of the Laplacian with
      weights 4 J_b^2/(gamma_i + gamma_j); the deviation times the cube of the overall rate scale s, over
      s = 50, 100, 200, follows s^k with the law k = -2, decided the same way, 2^(-5/2) < ratio < 2^(-3/2).
      Control: the weights 2 J_b^2/gamma_mean of a uniform profile, fed through the same comparison, give
      k = 2, a wrong leading term, in every block.
  Z6c What the next term does where it does not vanish: on the triangle with the same Dzyaloshinskii-Moriya D
      on every bond along it (a flux 3 arctan(D/J) through it), block p = 1, the leading pair is degenerate and
      the gamma^-2 term splits it into a complex pair, |Im lambda| falling as gamma^-2 (ratio 2^2 per doubling of
      gamma, decided at the half-integer powers), the real parts of the pair printed; on the ring N = 5 with the
      same D, no triangle, |Im lambda| falls as gamma^-4 (ratio 2^4); each row is read against the other's law,
      which it must not meet, as its control.  The exact counterpart without a flux, the next
      term zero, is row Z2b's.
  Z7  Reading, not gated: rings N = 4..9, the complete graph and random graphs (where no leaf removal is available)
      under the same exact counter, each bracket's lower end checked to have only the zero mode below it: whether
      any block has an eigenvalue below lambda_1, and the multiplicity at lambda_1 per block (on the N = 4 ring a
      singlet joins the two lifted modes at lambda_1 = 2 in block p = 2).  The equality holds there by the route
      the proof recognizes afterwards (the exclusion process a quotient of the interchange process, whose gap
      Caputo, Liggett and Richthammer proved equal to lambda_1); this row reads it exactly.
  Z8  Reading, not gated: the chain prefactors that experiments/CHAIN_GAP_SECTOR_DIAGNOSTIC.md and
      hypotheses/F1_DISSIPATION_GAP_PATTERN.md quote, in their book H = (J/4) sum sigma.sigma, gamma (Z rho Z - rho):
      the slowest rate over the diagonal blocks (every gap printed lies below 2 gamma, so by Theorem A it is a real
      mode of a diagonal block), as <n_XY> N^2/Q^2 and gap N^2/(gamma Q^2), for the XXX and the XX chain at Q = 2,
      N = 4, 5, 6, the N = 5 sweep Q = 0.5 to 2.5 and the Q = 0.5 mean over N = 4, 5, 6, beside the Zeno values
      N^2 (1 - cos(pi/N))/8 and /4.

Run: python simulations/f50_zeno_end_ferromagnet.py   (about half a minute; prints "ALL GATES PASS" or a FAIL line)
"""
import sys
from fractions import Fraction
from itertools import combinations
from math import comb, factorial
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
FAILS = []
RNG = np.random.default_rng(20261008)


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}{(': ' + detail) if detail else ''}")
    if not ok:
        FAILS.append(name)


# ----------------------------------------------------------------------------------------------------------
# builders (site l is bit N-1-l)
# ----------------------------------------------------------------------------------------------------------
def bit(N, x, l):
    return (x >> (N - 1 - l)) & 1


def confs_of(N, p):
    return [x for x in range(2 ** N) if bin(x).count("1") == p]


def chain(N, w=None):
    return [(i, i + 1, 1 if w is None else w[i]) for i in range(N - 1)]


def ring(N):
    return [(i, (i + 1) % N, 1) for i in range(N)]


def star(N):
    return [(0, i, 1) for i in range(1, N)]


def complete(N):
    return [(i, j, 1) for i in range(N) for j in range(i + 1, N)]


def random_tree(N):
    return [(int(RNG.integers(0, v)), v, int(RNG.integers(1, 6))) for v in range(1, N)]


def random_graph(N, extra):
    """a random spanning tree plus `extra` random chords, integer weights"""
    E = {(min(a, b), max(a, b)): c for a, b, c in random_tree(N)}
    while len(E) < N - 1 + extra:
        a, b = sorted(RNG.choice(N, 2, replace=False).tolist())
        E.setdefault((a, b), int(RNG.integers(1, 6)))
    return [(a, b, c) for (a, b), c in E.items()]


def exclusion_laplacian(N, p, edges):
    """Integer Laplacian of the exclusion graph of block p: an edge of weight c per bond whose sites disagree."""
    cf = confs_of(N, p); pos = {x: k for k, x in enumerate(cf)}
    L = [[0] * len(cf) for _ in cf]
    for a, b, c in edges:
        for x in cf:
            if bit(N, x, a) != bit(N, x, b):
                y = x ^ (1 << (N - 1 - a)) ^ (1 << (N - 1 - b))
                L[pos[x]][pos[x]] += c
                L[pos[x]][pos[y]] -= c
    return L


def site_laplacian(N, edges):
    L = [[0] * N for _ in range(N)]
    for a, b, c in edges:
        L[a][a] += c; L[b][b] += c; L[a][b] -= c; L[b][a] -= c
    return L


I2 = np.eye(2, dtype=complex)
PX = np.array([[0, 1], [1, 0]], dtype=complex)
PY = np.array([[0, -1j], [1j, 0]], dtype=complex)
PZ = np.array([[1, 0], [0, -1]], dtype=complex)
SIGMA_MINUS = np.array([[0, 0], [1, 0]], dtype=complex)   # |0> -> |1>: adds one excitation


def kron_op(N, site_ops):
    out = np.ones((1, 1), dtype=complex)
    for l in range(N):
        out = np.kron(out, site_ops.get(l, I2))
    return out


def ferromagnet_times_two(N, edges, with_yy=True):
    """sum_b c_b (I - XX - YY - ZZ) = 2 sum_b c_b (1 - SWAP_b), from Kronecker products (integer entries)."""
    d = 2 ** N
    F = np.zeros((d, d), dtype=complex)
    for a, b, c in edges:
        F += c * (np.eye(d) - kron_op(N, {a: PX, b: PX}) - kron_op(N, {a: PZ, b: PZ})
                  - (kron_op(N, {a: PY, b: PY}) if with_yy else 0))
    return F


def ham_matrix(N, edges, delta=Fraction(1), dm=0):
    """H = sum_b c_b (XX + YY + delta ZZ) + dm (XY - YX) on the first bond; entries exact (Fraction) in a dict-of-dicts."""
    d = 2 ** N
    H = [dict() for _ in range(d)]

    def add(r, s, v):
        H[r][s] = H[r].get(s, 0) + v
    for k, (a, b, c) in enumerate(edges):
        for x in range(d):
            if bit(N, x, a) != bit(N, x, b):
                y = x ^ (1 << (N - 1 - a)) ^ (1 << (N - 1 - b))
                add(y, x, 2 * c)
                if dm and k == 0:
                    add(y, x, 2j * dm * (1 if bit(N, x, a) else -1))
                add(x, x, -delta * c)
            else:
                add(x, x, delta * c)
    return H


def add_diagonal(N, H, fields=(), zz=()):
    """Add sum_l h_l Z_l + sum K Z_a Z_b (terms diagonal in the configurations) to H in place."""
    for x in range(2 ** N):
        z = [1 - 2 * bit(N, x, l) for l in range(N)]
        v = sum(Fraction(h) * z[l] for l, h in enumerate(fields)) + sum(Fraction(K) * z[a] * z[b] for a, b, K in zz)
        if v:
            H[x][x] = H[x].get(x, 0) + v
    return H


def ham_correlated_hop(N, edges, k):
    """XX + YY + ZZ, the hop on the first bond doubled when site k is occupied (number conserving, not SU(2))."""
    d = 2 ** N
    H = [dict() for _ in range(d)]
    for idx_b, (a, b, c) in enumerate(edges):
        for x in range(d):
            if bit(N, x, a) != bit(N, x, b):
                y = x ^ (1 << (N - 1 - a)) ^ (1 << (N - 1 - b))
                amp = 2 * c * (2 if (idx_b == 0 and bit(N, x, k)) else 1)
                H[y][x] = H[y].get(x, 0) + amp
                H[x][x] = H[x].get(x, 0) - c
            else:
                H[x][x] = H[x].get(x, 0) + c
    return H


class GQ:
    """An exact Gaussian rational a + b i with Fraction parts."""
    __slots__ = ("re", "im")

    def __init__(self, re, im=0):
        if isinstance(re, complex):
            re, im = Fraction(re.real), Fraction(re.imag)   # integer-valued floats only (asserted below)
        self.re, self.im = Fraction(re), Fraction(im)

    @staticmethod
    def of(v):
        if isinstance(v, GQ):
            return v
        if isinstance(v, complex):
            assert v.real == int(v.real) and v.imag == int(v.imag), v
        return GQ(v)

    def __add__(self, o):
        o = GQ.of(o); return GQ(self.re + o.re, self.im + o.im)
    __radd__ = __add__

    def __sub__(self, o):
        o = GQ.of(o); return GQ(self.re - o.re, self.im - o.im)

    def __rsub__(self, o):
        return GQ.of(o) - self

    def __neg__(self):
        return GQ(-self.re, -self.im)

    def __mul__(self, o):
        o = GQ.of(o); return GQ(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = GQ.of(o)
        n = o.re * o.re + o.im * o.im
        return GQ((self.re * o.re + self.im * o.im) / n, (self.im * o.re - self.re * o.im) / n)

    def __eq__(self, o):
        o = GQ.of(o); return self.re == o.re and self.im == o.im

    def __ne__(self, o):
        return not self.__eq__(o)

    __hash__ = None


def to_gq(H):
    return [{k: GQ.of(v) for k, v in row.items()} for row in H]


def zeno_terms(N, H, p, gam):
    """Exact A_PQ Gamma^-1 A_QP and A_PQ Gamma^-1 A_QQ Gamma^-1 A_QP on block (p, p); Gamma_xy = 2 sum gamma over
    the sites where x and y disagree (gam rational)."""
    cf = confs_of(N, p); pos = {x: k for k, x in enumerate(cf)}
    H_cols = [dict() for _ in range(2 ** N)]
    for r in range(2 ** N):
        for c_, v in H[r].items():
            H_cols[c_][r] = v

    def A_vec(vec):
        out = {}
        for (x, y), a in vec.items():
            for xp, v in H_cols[x].items():
                if xp in pos:
                    out[(xp, y)] = out.get((xp, y), 0) + v * a
            for yp, v in H[y].items():
                if yp in pos:
                    out[(x, yp)] = out.get((x, yp), 0) - v * a
        return {k_: v for k_, v in out.items() if v != 0}

    def rate(x, y):
        return 2 * sum(gam[l] for l in range(N) if bit(N, x, l) != bit(N, y, l))

    n = len(cf)
    M1 = [[0] * n for _ in range(n)]
    M3 = [[0] * n for _ in range(n)]
    for j, x in enumerate(cf):
        first = {k_: v / rate(*k_) for k_, v in A_vec({(x, x): 1}).items()}           # Gamma^-1 A_QP, all in Q
        second = A_vec(first)
        for (u, w), v in second.items():
            if u == w:
                M1[pos[u]][j] += v
        mid = {k_: v / rate(*k_) for k_, v in second.items() if k_[0] != k_[1]}     # Gamma^-1 A_QQ Gamma^-1 A_QP
        third = A_vec(mid)
        for (u, w), v in third.items():
            if u == w:
                M3[pos[u]][j] += v
    return M1, M3


def ad_block_squared_pp(N, H, p):
    """(A^2)_PP for A = H (x) I - I (x) H^T on block (p, p), exact (Fraction / complex Fraction pairs)."""
    cf = confs_of(N, p); pos = {x: k for k, x in enumerate(cf)}

    def A_col(x, y):
        out = {}
        for xp, v in H_cols[x].items():
            if xp in pos:
                out[(xp, y)] = out.get((xp, y), 0) + v
        for yp, v in H[y].items():            # (rho H)_{x, yp} = sum_y rho_{x y} H_{y yp}
            if yp in pos:
                out[(x, yp)] = out.get((x, yp), 0) - v
        return out
    H_cols = [dict() for _ in range(2 ** N)]  # H_cols[x][xp] = H[xp][x]
    for r in range(2 ** N):
        for s, v in H[r].items():
            H_cols[s][r] = v
    n = len(cf)
    M = [[0] * n for _ in range(n)]
    for j, x in enumerate(cf):
        first = A_col(x, x)
        second = {}
        for (u, v), a1 in first.items():
            for cell, a2 in A_col(u, v).items():
                second[cell] = second.get(cell, 0) + a2 * a1
        for i, z in enumerate(cf):
            M[i][j] = second.get((z, z), 0)
    return M


def block_of(Mfull, N, p):
    cf = confs_of(N, p)
    return Mfull[np.ix_(cf, cf)]


# ----------------------------------------------------------------------------------------------------------
# exact eigenvalue counting: Jacobi's rule on the leading principal minors, Bareiss over the integers
# ----------------------------------------------------------------------------------------------------------
def leading_minors(M):
    A = [row[:] for row in M]; n = len(A); minors = []; prev = 1
    for k in range(n):
        if A[k][k] == 0:
            return None
        minors.append(A[k][k])
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev   # exact division (Bareiss)
        prev = A[k][k]
    return minors


def count_below(L, q):
    """# eigenvalues of the integer symmetric matrix L strictly below the rational q (exact; None on a zero minor)."""
    a, b = q.numerator, q.denominator
    n = len(L)
    m = leading_minors([[b * L[i][j] - (a if i == j else 0) for j in range(n)] for i in range(n)])
    if m is None:
        return None
    seq = [1] + m
    return sum(1 for k in range(n) if (seq[k] > 0) != (seq[k + 1] > 0))


def rational_in(lo, hi):
    """a rational in the middle half of (lo, hi); floats only place it, the counts verify it exactly"""
    for den in (10 ** 4, 10 ** 6, 10 ** 8, 10 ** 10, 10 ** 12):
        q = Fraction(round((lo + hi) / 2 * den), den)
        if lo + (hi - lo) / 4 < float(q) < hi - (hi - lo) / 4:
            return q
    raise RuntimeError("window too narrow")


def remove_vertex(N, edges, v):
    keep = [k for k in range(N) if k != v]; ren = {k: n for n, k in enumerate(keep)}
    return N - 1, [(ren[a], ren[b], c) for a, b, c in edges if v not in (a, b)]


def leaves(N, edges):
    deg = [0] * N
    for a, b, c in edges:
        deg[a] += 1; deg[b] += 1
    return [v for v in range(N) if deg[v] == 1]


def site_spectrum(N, edges):
    return np.linalg.eigvalsh(np.array(site_laplacian(N, edges), float))


# ----------------------------------------------------------------------------------------------------------
# Z1: the exclusion Laplacian is the ferromagnet
# ----------------------------------------------------------------------------------------------------------
print("Z1  2 Lap(exclusion graph, block p) == block p of sum_b c_b (I - XX - YY - ZZ), exact")
graphs_z1 = [(f"chain N={N}", N, chain(N)) for N in (3, 4, 5, 6)] + [(f"ring N={N}", N, ring(N)) for N in (4, 5)] + \
            [(f"star N={N}", N, star(N)) for N in (4, 5)] + [("K5", 5, complete(5))] + \
            [(f"random graph N=6 #{t}", 6, random_graph(6, 3)) for t in range(3)]
for label, N, E in graphs_z1:
    F2 = ferromagnet_times_two(N, E)
    ok = np.array_equal(F2.imag, np.zeros_like(F2.imag))
    for p in range(N + 1):
        ok &= np.array_equal(block_of(F2.real, N, p), 2 * np.array(exclusion_laplacian(N, p, E), float))
    check(f"{label}: every block p = 0..{N}", ok)
N, E = 5, random_graph(5, 2)
F2bad = ferromagnet_times_two(N, E, with_yy=False)
bad = any(not np.array_equal(block_of(F2bad.real, N, p), 2 * np.array(exclusion_laplacian(N, p, E), float)) for p in range(1, N))
check("control: the operator built without its YY term differs from 2 Lap in some block 1 <= p <= N-1", bad)

# ----------------------------------------------------------------------------------------------------------
# Z2: the Zeno generator of every XXZ is the Laplacian with weights J_b^2
# ----------------------------------------------------------------------------------------------------------
print("Z2  (A^2)_PP == 8 Lap(weights J_b^2) exactly, every Delta")
for N, p, prof in ((4, 1, [1, 1, 1]), (4, 2, [1, 1, 1]), (5, 2, [1, 2, 3, 1]), (6, 3, [2, 1, 1, 3, 1])):
    E = chain(N, prof)
    Lsq = exclusion_laplacian(N, p, [(a, b, c * c) for a, b, c in E])
    for delta in (Fraction(0), Fraction(1, 2), Fraction(1), Fraction(2), Fraction(-1)):
        M = ad_block_squared_pp(N, ham_matrix(N, E, delta=delta), p)
        same = all(M[i][j] == 8 * Lsq[i][j] for i in range(len(M)) for j in range(len(M)))
        check(f"chain N={N} profile {prof} ({p},{p}) Delta = {delta}", same)
for N, p in ((5, 2),):
    E = random_graph(N, 2)
    Lsq = exclusion_laplacian(N, p, [(a, b, c * c) for a, b, c in E])
    M = ad_block_squared_pp(N, ham_matrix(N, E, delta=Fraction(3, 2)), p)
    check(f"random graph N={N} ({p},{p}) Delta = 3/2", all(M[i][j] == 8 * Lsq[i][j] for i in range(len(M)) for j in range(len(M))))
N = 5
E = chain(N, [1, 2, 3, 1])
H_diag = add_diagonal(N, ham_matrix(N, E, delta=Fraction(1, 2)), fields=(2, -1, 3, 0, -2), zz=((0, 2, 2), (1, 3, -1), (2, 4, 3)))
ok = True
for p in (1, 2):
    Lsq = exclusion_laplacian(N, p, [(a, b, c * c) for a, b, c in E])
    M = ad_block_squared_pp(N, H_diag, p)
    ok &= all(M[i][j] == 8 * Lsq[i][j] for i in range(len(M)) for j in range(len(M)))
check("chain N=5, Delta = 1/2 with the fields (2, -1, 3, 0, -2) and next-nearest ZZ: (A^2)_PP == 8 Lap(J_b^2), blocks 1, 2", ok)
N = 4
E = chain(N)
ok = True
for p in (1, 2):
    Lw = exclusion_laplacian(N, p, [(a, b, c * c + (1 if a == 0 else 0)) for a, b, c in E])   # D_0 = 1 on bond 0
    M = ad_block_squared_pp(N, ham_matrix(N, E, dm=1), p)
    ok &= all(M[i][j] == 8 * Lw[i][j] for i in range(len(M)) for j in range(len(M)))
check("chain N=4 with a DM term D = 1 on bond 0: (A^2)_PP == 8 Lap(J_b^2 + D_b^2), blocks 1, 2", ok)
E = chain(5, [1, 2, 3, 1])
Lw = exclusion_laplacian(5, 2, E)
M = ad_block_squared_pp(5, ham_matrix(5, E), 2)
check("control: the weights J_b in place of J_b^2 miss (A^2)_PP on the chain N=5 with profile [1, 2, 3, 1]",
      any(M[i][j] != 8 * Lw[i][j] for i in range(len(M)) for j in range(len(M))))


def zeno_full(N, H):
    """The Zeno generator (A^2)_PP of every block assembled on the 2^N configurations (exact lists)."""
    Z = [[0] * 2 ** N for _ in range(2 ** N)]
    for p in range(N + 1):
        cf = confs_of(N, p)
        M = ad_block_squared_pp(N, H, p)
        for i, x in enumerate(cf):
            for j, y in enumerate(cf):
                Z[x][y] = M[i][j]
    return Z


def commutes_with_lowering(N, Z):
    Sm = np.real(sum(kron_op(N, {l: SIGMA_MINUS}) for l in range(N))).astype(int).tolist()
    d = 2 ** N
    ZS = [[sum(Z[i][k] * Sm[k][j] for k in range(d) if Sm[k][j]) for j in range(d)] for i in range(d)]
    SZ = [[sum(Sm[i][k] * Z[k][j] for k in range(d) if Sm[i][k]) for j in range(d)] for i in range(d)]
    return all(ZS[i][j] == SZ[i][j] for i in range(d) for j in range(d))


rows = [("chain N=5, Delta = 1/2, fields and next-nearest ZZ", 5, H_diag),
        ("chain N=4, DM D = 1 on bond 0", 4, ham_matrix(4, chain(4), dm=1)),
        ("random graph N=5, Delta = 2", 5, ham_matrix(5, random_graph(5, 2), delta=Fraction(2)))]
for label, N, H in rows:
    check(f"{label}: the assembled Zeno generator commutes with S^- exactly", commutes_with_lowering(N, zeno_full(N, H)))
check("control: a correlated hop (doubled when site 2 is occupied) breaks the commutation with S^-",
      not commutes_with_lowering(5, zeno_full(5, ham_correlated_hop(5, chain(5), 2))))

# ----------------------------------------------------------------------------------------------------------
# Z2b: a rational profile of rates, the leading Schur term exact and the next one zero
# ----------------------------------------------------------------------------------------------------------
print("Z2b A_PQ Gamma^-1 A_QP == Lap(4 J_b^2/(gamma_i + gamma_j)) and A_PQ Gamma^-1 A_QQ Gamma^-1 A_QP == 0, exact")
gam_q = [Fraction(int(v), 7) for v in RNG.integers(3, 15, size=5)]
TRI = [(0, 1, 1), (1, 2, 2), (0, 2, 1), (2, 3, 1), (3, 4, 2)]
for label, N, E, H in (
        ("chain N=5, Delta = 1/2, fields, next-nearest ZZ", 5, chain(5, [1, 2, 1, 3]),
         add_diagonal(5, ham_matrix(5, chain(5, [1, 2, 1, 3]), delta=Fraction(1, 2)), fields=(1, 0, -2, 1, 3), zz=((0, 2, 1),))),
        ("triangle graph N=5 (bonds 01, 12, 02, 23, 34), Delta = 3/2, fields", 5, TRI,
         add_diagonal(5, ham_matrix(5, TRI, delta=Fraction(3, 2)), fields=(1, -2, 0, 3, 1)))):
    ok1, ok3 = True, True
    for p in (1, 2):
        M1, M3 = zeno_terms(N, H, p, gam_q)
        Lw = exclusion_laplacian(N, p, [(a, b, Fraction(4) * Fraction(c) ** 2 / (gam_q[a] + gam_q[b])) for a, b, c in E])
        ok1 &= all(M1[i][j] == Lw[i][j] for i in range(len(M1)) for j in range(len(M1)))
        ok3 &= all(v == 0 for row in M3 for v in row)
    check(f"{label}, gamma_l = {[str(g) for g in gam_q]}: leading term == Lap(4 J_b^2/(gamma_i + gamma_j)), blocks 1, 2", ok1)
    check(f"{label}: the next term vanishes, blocks 1, 2", ok3)
E = chain(5, [1, 2, 1, 3])
M1, _ = zeno_terms(5, ham_matrix(5, E), 2, gam_q)
gbar = sum(gam_q) / 5
Lwrong = exclusion_laplacian(5, 2, [(a, b, Fraction(2) * Fraction(c) ** 2 / gbar) for a, b, c in E])
check("control: the uniform-rate weights 2 J_b^2/gamma_mean miss the leading term",
      any(M1[i][j] != Lwrong[i][j] for i in range(len(M1)) for j in range(len(M1))))
_, M3f = zeno_terms(5, to_gq(ham_matrix(5, TRI, delta=Fraction(3, 2), dm=1)), 2, gam_q)
check("control: a DM flux through the triangle (H not real) makes the next term nonzero (exact, Gaussian rationals)",
      any(v != 0 for row in M3f for v in row))
# that term is i M3 with M3 purely imaginary and antisymmetric: real and antisymmetric on the populations, so it moves
# no real part (L commutes with X -> X^dagger, which conjugates the populations)
check("the flux term M3 is purely imaginary and antisymmetric, so i M3 is real and antisymmetric (exact)",
      all(GQ.of(v).re == 0 for row in M3f for v in row)
      and all(GQ.of(M3f[r][c]) == -GQ.of(M3f[c][r]) for r in range(len(M3f)) for c in range(len(M3f))))
# a complex H on a graph without a triangle: the next term vanishes all the same
for label, Nn, Hn in (("chain N=5 with a DM term", 5, ham_matrix(5, chain(5, [1, 2, 1, 3]), delta=Fraction(1, 2), dm=1)),
                      ("ring N=4 with a DM term (a flux)", 4, ham_matrix(4, ring(4), delta=Fraction(1), dm=1))):
    zero = all(v == 0 for p in (1, 2) for row in zeno_terms(Nn, to_gq(Hn), p, gam_q[:Nn])[1] for v in row)
    check(f"{label}, H not real, no triangle: the next term vanishes, blocks 1, 2 (exact, Gaussian rationals)", zero)

# ----------------------------------------------------------------------------------------------------------
# Z3: the uniform Heisenberg graph: Lap = (#bonds - H)/2, the Hamiltonian read from the top
# ----------------------------------------------------------------------------------------------------------
print("Z3  uniform Heisenberg: 2 Lap(block p) == #bonds I - H_p exactly")


def dense_block(H, N, p):
    cf = confs_of(N, p); pos = {x: k for k, x in enumerate(cf)}
    n = len(cf); M = [[0] * n for _ in range(n)]
    for x in cf:
        for y, v in H[x].items():
            if y in pos:
                M[pos[x]][pos[y]] = v
    return M


for label, N, E in ((f"chain N=6", 6, chain(6)), ("ring N=5", 5, ring(5)), ("star N=5", 5, star(5)), ("K4", 4, complete(4))):
    ok = True
    for p in range(N + 1):
        Hp = dense_block(ham_matrix(N, E, delta=Fraction(1)), N, p)
        L = exclusion_laplacian(N, p, E)
        n = len(L)
        ok &= all(2 * L[i][j] == (len(E) if i == j else 0) - Hp[i][j] for i in range(n) for j in range(n))
    check(f"{label}: every block", ok)
N = 6; E = chain(N); p = 3
Hp = dense_block(ham_matrix(N, E, delta=Fraction(2)), N, p); L = exclusion_laplacian(N, p, E); n = len(L)
check("control: Delta = 2 breaks the identity", any(2 * L[i][j] != (len(E) if i == j else 0) - Hp[i][j] for i in range(n) for j in range(n)))

# ----------------------------------------------------------------------------------------------------------
# Z4: the lift is SU(2)'s lowering
# ----------------------------------------------------------------------------------------------------------
print("Z4  (S^-)^(p-1) sum_i f(i)|i> == (p-1)! sum_S F(S)|S>;  [ferromagnet, S^-] == 0;  Lap F == lift(L_site f)")


def lowering(N, missing=None):
    return sum(kron_op(N, {l: SIGMA_MINUS}) for l in range(N) if l != missing)


for label, N, E in (("chain N=6", 6, chain(6, [1, 2, 1, 3, 1])), ("random tree N=7", 7, random_tree(7)), ("ring N=5", 5, ring(5))):
    f = RNG.integers(-5, 6, size=N)
    v = np.zeros(2 ** N, dtype=complex)
    for i in range(N):
        v[1 << (N - 1 - i)] = f[i]
    Sm = lowering(N)
    F2 = ferromagnet_times_two(N, E)
    ok_comm = np.array_equal(F2 @ Sm - Sm @ F2, np.zeros_like(F2))
    ok_lift, ok_eig = True, True
    Lsite = np.array(site_laplacian(N, E))
    w = v.copy()
    for p in range(1, N):
        cf = confs_of(N, p)
        F = np.array([sum(f[l] for l in range(N) if bit(N, x, l)) for x in cf])
        ok_lift &= np.array_equal(w[cf], factorial(p - 1) * F.astype(complex)) and \
            np.array_equal(np.delete(w, cf), np.zeros(2 ** N - len(cf)))
        Lf = Lsite @ f
        liftLf = np.array([sum(Lf[l] for l in range(N) if bit(N, x, l)) for x in cf])
        ok_eig &= np.array_equal(np.array(exclusion_laplacian(N, p, E)) @ F, liftLf)
        w = Sm @ w
    check(f"{label}: [2 sum c_b (1 - SWAP_b), S^-] == 0", ok_comm)
    check(f"{label}: (S^-)^(p-1) f == (p-1)! F in every block 1..N-1", ok_lift)
    check(f"{label}: Lap F == lift(L_site f) in every block 1..N-1", ok_eig)
N = 5; f = np.array([3, -1, 4, 1, -5])
v = np.zeros(2 ** N, dtype=complex)
for i in range(N):
    v[1 << (N - 1 - i)] = f[i]
w = lowering(N, missing=2) @ v
cf = confs_of(N, 2)
F = np.array([sum(f[l] for l in range(N) if bit(N, x, l)) for x in cf])
check("control: S^- with site 2 missing does not give the lift at p = 2", not np.array_equal(w[cf], F.astype(complex)))
# controls on the commutation and on the lift identity, through the same doors
Ec = chain(6, [1, 2, 1, 3, 1])
Sm6 = lowering(6)
F2a = ferromagnet_times_two(6, Ec) - sum(c * kron_op(6, {a: PZ, b: PZ}) for a, b, c in Ec)
check("control: the anisotropic sum_b c_b (I - XX - YY - 2 ZZ) does not commute with S^-",
      not np.array_equal(F2a @ Sm6 - Sm6 @ F2a, np.zeros_like(F2a)))
f6 = np.array([2, -3, 1, 4, -1, -3])
Lbad = np.array(site_laplacian(6, Ec)); Lbad[2, 2] += 1
Lf6 = Lbad @ f6
F6 = np.array([sum(f6[l] for l in range(6) if bit(6, x, l)) for x in confs_of(6, 2)])
liftbad = np.array([sum(Lf6[l] for l in range(6) if bit(6, x, l)) for x in confs_of(6, 2)])
check("control: a site Laplacian with one diagonal entry raised breaks Lap F == lift(L_site f) at p = 2",
      not np.array_equal(np.array(exclusion_laplacian(6, 2, Ec)) @ F6, liftbad))

# ----------------------------------------------------------------------------------------------------------
# Z5: the theorem on trees, by exact counts
# ----------------------------------------------------------------------------------------------------------
print("Z5  every tree, every block: only the zero mode below q_lo < lambda_1; the tree's count below q_hi and below q*")


def tree_rows(label, N, E, token_edges=None, extra=None, expect=True, miss_only_at=None):
    ev = site_spectrum(N, E)
    l1 = ev[1]
    above = ev[ev > l1 * (1 + 1e-9) + 1e-12]
    l2 = above[0] if len(above) else l1 + 1.0
    w = 1e-6 * max(1.0, l1)
    q_lo = rational_in(l1 - w, l1)
    q_hi = rational_in(l1, min(l1 + w, (l1 + l2) / 2))
    LT = site_laplacian(N, E)
    m_lo, m_hi = count_below(LT, q_lo), count_below(LT, q_hi)
    bracket_ok = (m_lo == 1 and m_hi is not None and m_hi >= 2)
    # q*: just under max over leaves of lambda_1(T - v), verified exactly on that leaf's Laplacian (N >= 3)
    qs, m_s = None, None
    if N >= 3:
        best = max(leaves(N, E), key=lambda v: site_spectrum(*remove_vertex(N, E, v))[1])
        Nm, Em = remove_vertex(N, E, best)
        lv = site_spectrum(Nm, Em)[1]
        qs = rational_in(lv * (1 - 2e-6), lv * (1 - 1e-6))
        bracket_ok &= count_below(site_laplacian(Nm, Em), qs) == 1
        m_s = count_below(LT, qs)
        bracket_ok &= m_s is not None
    ok = bracket_ok
    misses = []
    for p in range(1, N):
        Lt = exclusion_laplacian(N, p, E if token_edges is None else token_edges)
        if extra is not None:
            Xp = extra(N, p)
            Lt = [[Lt[i][j] + Xp[i][j] for j in range(len(Lt))] for i in range(len(Lt))]
        got = (count_below(Lt, q_lo), count_below(Lt, q_hi), count_below(Lt, qs) if qs is not None else None)
        if None in got[:2] or (qs is not None and got[2] is None) or got != (m_lo, m_hi, m_s):
            misses.append((p, got))
    ok &= not misses
    qs_note = ""
    if qs is not None:
        qs_note = f"; below q* = {float(qs):.6f}, under lambda_1(T - leaf): {m_s}"
        if float(qs) < l1:
            qs_note += " (q* below lambda_1: the leaf bound adds nothing here)"
    if expect:
        check(f"{label}: brackets exact, every block counts as the tree (zero mode; {m_hi - 1} at lambda_1 = {l1:.6f}"
              + qs_note + ")", ok, f"misses {misses}" if misses else "")
    else:
        # a control counts only with exact brackets and a definite count that differs (a zero minor is no miss);
        # miss_only_at names the brackets that must differ, the others agreeing
        want = (m_lo, m_hi, m_s)
        definite = [(p, g) for p, g in misses if None not in g[:2] and (qs is None or g[2] is not None)]
        if miss_only_at is not None:
            definite = [(p, g) for p, g in definite if {k for k in range(3) if g[k] != want[k]} == set(miss_only_at)]
        check(f"{label}: brackets exact, some block counts differently from the tree", bracket_ok and bool(definite),
              f"misses {misses}" if misses else "")
    return ok


for N in range(2, 10):
    tree_rows(f"path N={N}", N, chain(N))
for N in range(3, 10):
    tree_rows(f"star N={N}", N, star(N))
for N in range(3, 10):
    for t in range(2):
        tree_rows(f"random tree N={N} #{t}", N, random_tree(N))
for N in (6, 8):
    prof = [int(v) for v in RNG.integers(1, 4, size=N - 1)]
    tree_rows(f"chain N={N}, Zeno weights J_b^2 for J_b = {prof}", N, chain(N, [j * j for j in prof]))
# control: the same door, the blocks of a reweighted chain counted against the original chain's brackets
E_ctrl = chain(7)
tree_rows("control: chain N=7 brackets, blocks built with the bond (3, 4) doubled", 7, E_ctrl,
          token_edges=[(a, b, 2 * c if a == 3 else c) for a, b, c in E_ctrl], expect=False)
E_ctrl2 = chain(7, [2] * 6)
tree_rows("control: chain N=7 (bonds 2) brackets, blocks built with the bond (3, 4) halved", 7, E_ctrl2,
          token_edges=[(a, b, 1 if a == 3 else c) for a, b, c in E_ctrl2], expect=False)


# controls on what the theorem adds: -c (1 - SWAP_01)(1 - SWAP_23) is SU(2)-invariant and nonzero only on states
# with a singlet on both pairs, so it leaves the ferromagnetic multiplet and the one-magnon band in place and moves
# only the multiplets of lower spin, the ones the leaf bound puts at or above lambda_1(T - v)
def four_site(c):
    def extra(N, p):
        A, B = exclusion_laplacian(N, p, [(0, 1, 1)]), exclusion_laplacian(N, p, [(2, 3, 1)])
        n = len(A)
        return [[-c * sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    return extra


tree_rows("control: chain N=6 with -(1 - SWAP_01)(1 - SWAP_23) in every block, a new multiplet below lambda_1", 6,
          chain(6), extra=four_site(1), expect=False, miss_only_at=(0, 1, 2))
tree_rows("control: chain N=6, bonds 100, with -95 (1 - SWAP_01)(1 - SWAP_23), a new multiplet between lambda_1 and "
          "lambda_1(T - v): the q* count alone differs", 6, chain(6, [100] * 5), extra=four_site(95), expect=False,
          miss_only_at=(2,))
# the exact counter against a float count on a matrix with a known spectrum, far from the eigenvalues
Lp = exclusion_laplacian(6, 3, chain(6))
fl = np.linalg.eigvalsh(np.array(Lp, float))
qq = Fraction(7, 3)
check("cross-check: the exact counter agrees with a float count at q = 7/3 on the N=6 half-filling block",
      count_below(Lp, qq) == int(np.sum(fl < 7 / 3)) and np.abs(fl - 7 / 3).min() > 1e-3)

# ----------------------------------------------------------------------------------------------------------
# Z6: the full Liouvillian's Zeno limit is blind to Delta
# ----------------------------------------------------------------------------------------------------------
print("Z6  slow eigenvalues of B(x) = -2 x^2 ell_j + O(x^4) at Delta = 0, 1/2, 1, 2, the next order even (differences of "
      "deviation/x^4 in the ratio 4); a rate profile (deviation*s^3, differences in the ratio 1/4)")


def float_block(N, H, p):
    cf = confs_of(N, p)
    idx = [(x, y) for x in cf for y in cf]
    pos = {e: i for i, e in enumerate(idx)}
    n = len(idx)
    Hd = np.zeros((2 ** N, 2 ** N), dtype=complex)
    for r in range(2 ** N):
        for s, v in H[r].items():
            Hd[r, s] = complex(v)
    A = np.zeros((n, n), dtype=complex)
    for c, (x, y) in enumerate(idx):
        for xp in cf:
            if Hd[xp, x]:
                A[pos[(xp, y)], c] += Hd[xp, x]
        for yp in cf:
            if Hd[y, yp]:
                A[pos[(x, yp)], c] -= Hd[y, yp]
    ham = np.array([bin(x ^ y).count("1") for x, y in idx])
    return A, ham


XS = (0.005, 0.01, 0.02)


def x4_law(slow, ell_used):
    """deviation/x^4 over XS and the ratio of its successive differences, which names the power of the leading
    correction: 4 for x^2, 16 for x^4, 1/2 for x^-1 (an odd order), 1/4 for the x^-2 term 2 delta ell/x^2 that a
    leading generator off by a relative delta adds"""
    rho = [np.abs(slow[x] - np.sort(-2 * x * x * ell_used)).max() / x ** 4 for x in XS]
    return rho, (rho[2] - rho[1]) / (rho[1] - rho[0])


def names_x2(ratio):
    """the x^2 law against its neighbours x^1 and x^3, decided at the half-integer powers"""
    return 2 ** 1.5 < ratio < 2 ** 2.5


for N, p in ((4, 1), (4, 2), (5, 2)):
    E = chain(N)
    ell = np.sort(np.linalg.eigvalsh(np.array(exclusion_laplacian(N, p, E), float)))
    cases = [(str(d), ham_matrix(N, E, delta=d)) for d in (Fraction(0), Fraction(1, 2), Fraction(1), Fraction(2))]
    if N == 4:
        cases.append(("1 with fields (2, -1, 3, 0) and ZZ on (0, 2)",
                      add_diagonal(N, ham_matrix(N, E, delta=Fraction(1)), fields=(2, -1, 3, 0), zz=((0, 2, 1),))))
    for delta, Hc in cases:
        A, ham = float_block(N, Hc, p)
        slow = {x: np.sort(np.sort(np.linalg.eigvals(-1j * x * A + np.diag(-2.0 * ham)).real)[::-1][:comb(N, p)])
                for x in XS}
        rho, ratio = x4_law(slow, ell)
        _, ratio_c = x4_law(slow, ell * (1 + 1e-7))
        check(f"chain N={N} ({p},{p}) Delta = {delta}: deviation/x^4 = {rho[0]:.4f}, {rho[1]:.4f}, {rho[2]:.4f} at "
              f"x = 0.005, 0.01, 0.02, differences in the ratio {ratio:.4f} (x^2: 4, decided at the half-integer powers); "
              f"control ell (1 + 1e-7): {ratio_c:.4f}", names_x2(ratio) and not names_x2(ratio_c))

# Z6b: a profile of rates and couplings, the weights 4 J_b^2/(gamma_i + gamma_j)
N = 5
J_prof = [Fraction(int(v)) for v in RNG.integers(1, 4, size=N - 1)]
gam_prof = RNG.uniform(0.5, 2.0, size=N)
E = [(i, i + 1, J_prof[i]) for i in range(N - 1)]
H7 = ham_matrix(N, E, delta=Fraction(7, 10))


def profile_rows(weights_of, label, expect_law):
    ratios = []
    for p in range(1, N):
        A, _ = float_block(N, H7, p)
        cf = confs_of(N, p)
        rate = np.array([sum(gam_prof[l] for l in range(N) if bit(N, x, l) != bit(N, y, l)) for x in cf for y in cf])
        devs = {}
        for scale in (50.0, 100.0, 200.0):
            w = weights_of(scale)
            Lw = np.zeros((len(cf), len(cf)))
            pos = {x: k for k, x in enumerate(cf)}
            for (a, b, c), wb in zip(E, w):
                for x in cf:
                    if bit(N, x, a) != bit(N, x, b):
                        y = x ^ (1 << (N - 1 - a)) ^ (1 << (N - 1 - b))
                        Lw[pos[x], pos[x]] += wb
                        Lw[pos[x], pos[y]] -= wb
            ev = np.linalg.eigvals(-1j * A + np.diag(-2.0 * scale * rate))
            slow = np.sort(np.sort(ev.real)[::-1][:len(cf)])
            devs[scale] = np.abs(slow - np.sort(-np.linalg.eigvalsh(Lw))).max()
        rho = [devs[sc] * sc ** 3 for sc in (50.0, 100.0, 200.0)]
        ratios.append((rho[2] - rho[1]) / (rho[1] - rho[0]))
    holds = [2 ** -2.5 < r < 2 ** -1.5 for r in ratios]    # s^-2 against s^-3 and s^-1, at the half-integer powers
    ok = all(holds) if expect_law else not any(holds)
    check(f"{label}: deviation*s^3 over s = 50, 100, 200, differences in the ratio (s^-2: 1/4, decided at the "
          f"half-integer powers; a wrong leading term: 4) per block p = 1..{N - 1}: {', '.join(f'{r:.4f}' for r in ratios)}", ok)


profile_rows(lambda scale: [4 * float(J_prof[i]) ** 2 / (scale * (gam_prof[i] + gam_prof[i + 1])) for i in range(N - 1)],
             f"chain N={N}, Delta = 7/10, J_b = {[int(j) for j in J_prof]}, gamma_l = {np.round(gam_prof, 3).tolist()}: "
             "weights 4 J_b^2/(gamma_i + gamma_j), the s^-3 law", True)
profile_rows(lambda scale: [2 * float(J_prof[i]) ** 2 / (scale * gam_prof.mean()) for i in range(N - 1)],
             "control: the uniform-profile weights 2 J_b^2/gamma_mean miss the s^-3 law", False)

# Z6c: what the next term does where it does not vanish, on a triangle with a flux
print("Z6c the triangle with a flux: the gamma^-2 term splits the degenerate pair (|Im| ~ gamma^-2); no triangle, no such term")


def ham_flux_ring(N, J, D):
    """J (XX + YY + ZZ) around the ring 0 -> 1 -> ... -> N-1 -> 0 with the same Dzyaloshinskii-Moriya D on every
    bond, oriented along the ring: a flux N arctan(D/J) through it"""
    edges = [(l, (l + 1) % N, J) for l in range(N)]
    H = ham_matrix(N, edges, delta=Fraction(1))
    for a, b, c in edges:
        for x in range(2 ** N):
            if bit(N, x, a) != bit(N, x, b):
                y = x ^ (1 << (N - 1 - a)) ^ (1 << (N - 1 - b))
                H[y][x] = H[y].get(x, 0) + 2j * D * (1 if bit(N, x, a) else -1)
    return H


for N, k_law, k_other in ((3, 2, 4), (5, 4, 2)):
    A_f, ham_f = float_block(N, ham_flux_ring(N, 1, 1), 1)
    im, re = [], []
    for g in (40.0, 80.0, 160.0):
        ev = np.linalg.eigvals(-1j * A_f + np.diag(-2.0 * g * ham_f))
        pair = sorted(ev, key=lambda z: -z.real)[1:3]
        im.append(max(abs(z.imag) for z in pair)); re.append(tuple(round(float(z.real * g), 4) for z in pair))
    ratios = [im[0] / im[1], im[1] / im[2]]
    check(f"{'triangle' if N == 3 else 'ring N=5, no triangle'} with D = J on every bond, block p = 1, the leading "
          f"pair: Re*gamma = {re} at gamma = 40, 80, 160; |Im| per doubling of gamma falls by {ratios[0]:.3f}, "
          f"{ratios[1]:.3f} (gamma^-{k_law}: {2 ** k_law}, decided at the half-integer powers; control: not read as "
          f"gamma^-{k_other})",
          all(2 ** (k_law - 0.5) < r < 2 ** (k_law + 0.5) for r in ratios)
          and not all(2 ** (k_other - 0.5) < r < 2 ** (k_other + 0.5) for r in ratios))

# ----------------------------------------------------------------------------------------------------------
# Z7: readings beyond trees
# ----------------------------------------------------------------------------------------------------------
print("Z7  reading: graphs without the leaf argument, the same exact counter")
for label, N, E in [(f"ring N={N}", N, ring(N)) for N in range(4, 10)] + [("K6", 6, complete(6))] + \
                   [(f"random graph N=7 #{t}", 7, random_graph(7, 4)) for t in range(2)]:
    ev = site_spectrum(N, E); l1 = ev[1]
    q_lo = rational_in(l1 * (1 - 1e-6), l1); q_hi = rational_in(l1, l1 * (1 + 1e-6))
    LT = site_laplacian(N, E); m_lo, m_hi = count_below(LT, q_lo), count_below(LT, q_hi)
    assert m_lo == 1, (label, m_lo)   # the bracket's lower end has only the zero mode below it
    rows = [(count_below(exclusion_laplacian(N, p, E), q_lo), count_below(exclusion_laplacian(N, p, E), q_hi)) for p in range(1, N)]
    gap_same = all(r[0] == m_lo for r in rows)
    mult = [r[1] - r[0] for r in rows]
    print(f"    {label}: lambda_1 = {l1:.6f} (multiplicity {m_hi - m_lo}); nothing below it in any block: {gap_same}; "
          f"multiplicity at lambda_1 per block p = 1..{N - 1}: {mult}")

# ----------------------------------------------------------------------------------------------------------
# Z8: the chain prefactors of CHAIN_GAP_SECTOR_DIAGNOSTIC and F1_DISSIPATION_GAP_PATTERN, in their spin book
# ----------------------------------------------------------------------------------------------------------
print("Z8  reading: the chain prefactors in the book H = (J/4) sum sigma.sigma, gamma (Z rho Z - rho), gamma = 1/2")


def spin_book_gap(N, Q, delta, g=0.5):
    """Slowest nonzero rate over the diagonal blocks of the chain, H = (J/4)(XX + YY + delta ZZ), J = Q g."""
    best = np.inf
    H = ham_matrix(N, chain(N), delta=Fraction(delta))
    for p in range(1, N):
        A, ham = float_block(N, H, p)
        ev = np.linalg.eigvals(-1j * (Q * g / 4) * A + np.diag(-2.0 * g * ham)).real
        r = -ev[-ev > 1e-9]
        best = min(best, r.min())
    return best


g = 0.5
for N in (4, 5, 6):
    z8 = N * N * (1 - np.cos(np.pi / N)) / 8
    gx, g0 = spin_book_gap(N, 2.0, 1), spin_book_gap(N, 2.0, 0)
    print(f"    Q = 2, N = {N}: <n_XY> N^2/Q^2 = {gx * N * N / (2 * g * 4):.4f} (XXX), {g0 * N * N / (2 * g * 4):.4f} (XX), "
          f"Zeno {z8:.4f};  gap N^2/(gamma Q^2) = {gx * N * N / (g * 4):.4f}, {g0 * N * N / (g * 4):.4f}, Zeno {2 * z8:.4f};  "
          f"gaps below 2 gamma: {max(gx, g0) < 2 * g}")
print("    N = 5 XXX, <n_XY> N^2/Q^2 at Q = 0.5, 1, 1.5, 2, 2.5: "
      + ", ".join(f"{spin_book_gap(5, Q, 1) * 25 / (2 * g * Q * Q):.4f}" for Q in (0.5, 1.0, 1.5, 2.0, 2.5))
      + f"  (Zeno {25 * (1 - np.cos(np.pi / 5)) / 8:.4f})")
m = np.mean([spin_book_gap(N, 0.5, 1) * N * N / (g * 0.25) for N in (4, 5, 6)])
print(f"    Q = 0.5, XXX, gap N^2/(gamma Q^2) mean over N = 4, 5, 6: {m:.4f}  "
      f"(Zeno mean {np.mean([N * N * (1 - np.cos(np.pi / N)) / 4 for N in (4, 5, 6)]):.4f})")

print()
if FAILS:
    print("FAILED:", *FAILS, sep="\n  ")
    sys.exit(1)
print("ALL GATES PASS")
