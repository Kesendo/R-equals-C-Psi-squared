"""
The mirror under a shared field: the palindrome's price becomes an operator
===========================================================================
experiments/THE_MIRROR_UNDER_A_SHARED_FIELD.md. N qubits, Z dephasing whose site fields may be
correlated: the Kossakowski form

    D_a(rho) = sum_kl a_kl (Z_k rho Z_l - 1/2 {Z_l Z_k, rho}),   a real symmetric,

with a diagonal for local dephasing (F1's setting) and a = gamma [(1 - eta) 1 + eta 11^T] for a field
shared by every pair with correlation coefficient eta (the N-site form of
experiments/SINGLET_UNDER_A_SHARED_FIELD.md). Pi_Z(rho) = rho^T X^N is the canonical palindromizer
(docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md); K = sum_kl a_kl Z_l Z_k is sum_mu c_mu^dag c_mu for any
choice of jumps realising a.

The script checks, exactly:
  A. the rate law: D_a is diagonal on the matrix units |i><j|, with entry -1/2 u^T a u, u = z(i) - z(j),
     z_k = +-1 the Z eigenvalue of site k; for the shared field -2 gamma [(1 - eta) |D| + eta (p - q)^2],
     |D| the number of sites where i and j differ, p, q the popcounts of ket and bra; on the aligned coherences and their Pi_Z
     images, the price-pair page's pair sum 2 sigma + 4 sum_{k<l same side} a_kl;
  B. the block law: on every joint-popcount block (p, q) of a popcount-conserving H,
     L_eta = L_local((1 - eta) gamma) - 2 gamma eta (p - q)^2 Id, entry by entry, and nothing leaves a block;
     so every (p, p) block runs at the rate (1 - eta) gamma, and at eta = 1 it is unitary; and the same at
     site-dependent rates with one common covariance, a = diag(g) + c (11^T - 1): local at g_k - c,
     shifted by -2c (p - q)^2, pairing about (sigma - N c) + c (d^2 + m^2), d = p - q, m = p + q - N,
     the pairing also with DM bonds on the open chain (N = 3, 4, characteristic polynomials);
     B'''. the four edge blocks (0, 1), (1, 0), (N - 1, N), (N, N - 1) run exactly as under local noise for
     every covariance pattern, their images not;
     and for thirteen covariance patterns that are not one number (N = 3 and 4) the blocks that pair with their
     images are exactly those whose two popcount sectors each see one value of K (sufficient by the price
     identity; popcounts 1 and N - 1 qualify exactly when the covariance row sums are equal, 16 of 25 for such
     patterns at N = 4, and popcount 2 alone for one sign pattern, 9 of 25);
  C. the price identity: Pi_Z L Pi_Z^-1 = -L - {K, .} for every real symmetric a (positive definite and
     indefinite) and every real H with [H, X^N] = 0, and K = sigma 1 for the diagonal part of a and not for a itself (F1; the iff is the page's argument);
     controls: a field on one site, an
     imaginary X^N-even term, and a Hermitian non-real a (Absorption section 2's pair (Z_0 +- i Z_1)),
     whose K is a number and whose residual is exactly twice its Pi_Z-fixed imaginary part, each break
     it through the same construction; for the imaginary bond, F108's Pi_5bilinear (Pi_Z after conjugation
     by Y^N) holds it exactly, and the residuals add in squares, ||M||^2 = ||M_H||^2 + ||M_D||^2;
  D. the residual: M_D = Pi_Z L_D Pi_Z^-1 + L_D + 2 tr(a) = -{K - tr(a), .}, with
     ||M_D||_F^2 = 8 4^N sum_{k<l} a_kl^2 (shared field: 4^(N+1) N(N-1) gamma^2 eta^2), the norm
     convention fixed by reproducing the registry's T1 residual 4^(N-1) [3 sum g^2 + 4 (sum g)^2]; M_D is
     Pi_Z-invariant, so F81's relation gains no dissipator term;
  E. the blockwise pairing: Pi_Z restricted to block (p, q) carries L there onto -L - 2 c(p, q) on block
     (q, N - p), c(p, q) = N (1 - eta) gamma + gamma eta (d^2 + m^2); the bra fold |i><j| -> |i><j_bar|
     taken antiunitarily (F89d's form) carries the same centre onto block (p, N - q) at Delta = 0, 2/5,
     1, 3, and taken linearly with the staggered gauge (the fold lattice's linear leg) it is exact at
     Delta = 0 and broken at Delta = 2/5; by a second route that never uses Pi_Z, exact characteristic
     polynomials over Q(i), p_B(y) = (-1)^n p_A(-2c - y) for every block and its image at N = 3 and 4;
     about the one centre F137 allows the whole spectrum, trace/dim = -N gamma, exactly the blocks with
     d^2 + m^2 = N pair (none at N = 3, four at N = 4), and the whole characteristic polynomial at N = 2
     and 3 pairs there only at eta = 0; for eta > 0 the GHZ coherence's partner would sit at
     2 gamma eta N (N - 1) > 0, which breaks the global palindrome at every N (read at N = 2..5), and for
     any covariance pattern at +4 sum_{k<l} a_kl (read on random a);
  F. a uniform longitudinal field h sum_k Z_k on top: the pair sum gains 2 i h (2p - N), at the shared
     field and at site-dependent local rates, the same operator identity and the same characteristic-
     polynomial route (without the imaginary part no block pairs); a real non-uniform field sum_k h_k Z_k
     is not a constant per block pair, and of the blocks 8 of 16 (N = 3) and 12 of 25 (N = 4) still pair
     with their image at the shift the traces allow for a generic field, 12 and 16 for an anti-palindromic
     one (h_k + h_{N-1-k} constant): for the anti-palindromic field exactly the blocks with an empty or
     full ket or bra (Pi_Z, and the site reversal composed with it, leave the full residual a multiple of 1
     there), for the generic field the empty- or full-ket blocks plus, at even N, (N/2, 0) and (N/2, N),
     whose images are their own conjugate blocks, and these two drop out at site-dependent rates;
  G. the two ladders of the arc sideways_spin_ladder on the open XY chain: Pi_Z W Pi_Z^-1 = S+ exactly
     (the transpose leaves W fixed, the bra flip rho -> rho X^N carries it onto S+), with
     W(rho) = sum_l c_l^dag rho c_l and S+(rho) = sum_l (-1)^l c_l^dag rho c_l^dag (sites l counted from
     0; counted from 1, the image is -S+); W commutes with L for every real symmetric a; [S+, L] =
     -[S+, {K, .}], so S+ commutes for every diagonal a and fails at every correlated one (the page's
     argument; checked here on random a, the shared field and a single small covariance), the whole
     failure sitting in the off-diagonal part of a; under the shared field the failure is a shift only,
     [L, S+] = -8 kappa (d + 1) S+ on each source block, and for random a it is not; control: the
     Heisenberg chain breaks W through the Hamiltonian.
  H. the dark cells: for a positive semidefinite a, a matrix unit |i><j| has dissipator
     rate zero exactly when u = z(i) - z(j) lies in the kernel of a, read on random positive definite a,
     the shared field at eta = 3/10 and 1, the anti-correlated end eta = -1/(N - 1) (a_kk = N - 1,
     a_kl = -1, every row summing to zero), a random weighted graph Laplacian (a >= 0, every row summing
     to zero, its weights not all equal from N = 3 on), and a random rank-deficient a whose kernel is
     one random {0, +-1} direction that can single out particular sites from N = 3 on: 2^N cells for a
     positive definite a, C(2N, N) at eta = 1, and 2^N + 2 at that end, the populations plus
     |0...0><1...1| and its transpose; and on all six the GHZ coherence dark exactly when every row of
     a sums to zero (N = 2..5). The dark classes, the configurations grouped by a z: on all six the dark
     cells number the sum of the squared class sizes, at eta = 1 the classes are the popcount sectors and
     at the anti-correlated end the one class beyond single configurations is the GHZ pair; anti-correlated
     pairs a = [[1, -1], [-1, 1]] multiply the dark cells, (2^2 + 2)^(N/2) at N = 2 and 4, the largest
     class taking 00 or 11 on every pair (N = 2 is the anti-correlated end itself); a ring of covariances
     at gamma = 2|kappa| (its Laplacian D - A, its signless Laplacian D + A) leaves the GHZ pair at
     N = 3, 4, 5 (N = 3 the anti-correlated end again), the Neel pair at N = 4 and nothing beyond the
     populations at N = 3 and 5, and at N = 4 the flip on every other site carries the one set of dark
     cells onto the other; the open XY chain at N = 4 conserves each popcount, keeps the span of the
     popcount sectors (the commutator of every dark cell at eta = 1 stays inside them), commutes with the
     GHZ cell and carries the Neel cell out of the dark cells.
  I. which popcounts see one value of K: K(A) = sigma + 2T - 4 cut(A), and with a_kl = kbar + x_k + x_l
     + Omega_kl (sum x = 0, every row of Omega summing to zero) cut(A) = kbar p(N - p) + (N - 2p) x(A)
     - 2 Omega(A); K is one number on popcount p (1 <= p <= N - 1) exactly when x = 0 at p = 1, N - 1,
     Omega = 0 at p = N/2, and x = Omega = 0 elsewhere, checked against K = z^T a z itself on every integer pattern with
     entries in -2..2 at N = 3 and 4 and in -1..1 at N = 5, each at unequal rates, and on additive,
     equal-row-sum, mixed and random patterns at N = 6, 7, 8; the cut form of K, the decomposition
     (sum x = 0, the rows of Omega summing to zero, the parts adding back to the covariances) and the cut formula read set by
     set at N = 4..6; the shapes seen (two at N = 3, four at N = 4, three at N = 5); on every qualifying popcount K_p = sigma - N kbar + kbar (N - 2p)^2, the shared field at
     the mean covariance; for an additive C, K = sigma - 2 sum y + 2 (N - 2p) sum_k y_k z_k exactly; the
     page's sign pattern s = (1, 1, 1, -1) additive; on the open XY chain at N = 4 and 5, K one number on
     popcount p exactly when it commutes with the chain there; and the ZZ faces, ZZ couplings on the bonds
     of a ring, a star and an open chain one number on popcount p exactly where the law says (ring
     {0, 1, N - 1, N}, star {0, N/2, N} at even N, chain {0, N} from N = 4).
  J. the third moment about the one candidate centre -sigma: for every real symmetric a (either sign) and
     every traceless Hermitian H, Tr((L + sigma)^3) = -4^N [4 tr(C(2a + Lambda)C) - 3 tr(C Theta) - 6 xi^T a xi],
     C the covariances, Lambda the diagonal, Theta_kl = Tr(H^2 Z_k Z_l)/2^N, xi_k = Tr(H Z_k)/2^N, read against
     the exact trace of (2L + 2 sigma)^3 for random complex H at N = 2..4, for X^N-even H (xi = 0), and for
     local noise with fields, where it is F117's third-moment face 6 * 4^N sum gamma_k xi_k^2 carried to
     site-dependent rates; at H = 0 for N = 2..5; the noise part as 12 tr(Lambda C^2) + 8 tr(C^3) and in triangle
     form at N = 2..8, and as 8 tr(CaC) + 4 tr(C Lambda C) > 0 on positive definite draws (an indefinite a turns it
     negative, -48 at N = 3); Theta on the open XY and Heisenberg chains at N = 3..6; the shared field's
     closed forms on both chains at N = 2..4; a negative covariance between the end sites of the XY chain
     at N = 3, where the GHZ argument says nothing: the third moment sees the break, and the whole
     characteristic polynomial does not pair; and the shared field at the anti-correlated end with J = gamma
     at N = 3, where the third moment is exactly 0, the fifth is -11520 gamma^5 (the value of the companion
     script shared_field_odd_moments.py, met here by matrix powers), the far-end mode at p = 1 sits below
     -2 sigma, and the whole characteristic polynomial does not pair.
  K. the far-end modes: for a popcount-conserving X^N-even H and a pattern on which K takes one value K_p on
     popcount p, P_p X^N (the bra-flipped populations; F22's X^N P_p is its adjoint) is an eigenvector at -2 K_p, read on the XXZ chain under one common covariance at
     N = 2..5 for eta = -1/5 and 3/10 (K_p = sigma - N kappa + kappa (N - 2p)^2, above sigma for kappa < 0 at
     p = floor(N/2), so a rate above 2 sigma whose partner about -sigma would grow), on an additive pattern
     with a negative sum at N = 4 (K = sigma - N kbar on popcount 2), the law K_p > sigma exactly when
     T ((N - 2p)^2 - N) > 0 on every popcount where K is one number, on every integer pattern in -1..1 at
     N = 3..5 (T the summed covariance), and the end-site covariance at N = 3..5, where K is one number only
     on popcounts 0 and N and there lies below sigma, so neither the GHZ coherence nor a far-end mode decides, while the closed-form third
     moment on the XY chain is nonzero.

Every input is an integer (gamma = 10, eta in {3/10, -2/10, 9/10, 1}, so a has integer entries; the
anti-correlated end of section H at gamma = N - 1; the Hamiltonian integer), and the script works
with 2L, whose entries are Gaussian integers. Every comparison is an integer comparison; no tolerance
appears in any gate (a float eigenvalue only selects which random indefinite matrix to test). In
section H the kernel test a u = 0 is the commutation with every jump, since ker a = ker B for a = B^T B.

Import-inert; prints only, and exits with status 1 if any check fails.
"""

import itertools
import sys
from fractions import Fraction
from math import comb

import numpy as np
from scipy.sparse import csr_matrix
from sympy import I as SYM_I, Matrix as SymMatrix, Poly, Rational, symbols
from sympy.polys.domains import QQ_I
from sympy.polys.matrices import DomainMatrix

GAMMA = 10
ETAS = (Fraction(3, 10), Fraction(-2, 10), Fraction(9, 10), Fraction(1))
X = np.array([[0, 1], [1, 0]], dtype=np.int64)
Z = np.array([[1, 0], [0, -1]], dtype=np.int64)
YR = np.array([[0, -1], [1, 0]], dtype=np.int64)        # Y = i * YR
I2 = np.eye(2, dtype=np.int64)
CHECKS, CLAIMS = [], []


def claim(label, ok, detail=""):
    CLAIMS.append(bool(ok))
    print(f"  {label} {detail}  {'holds' if ok else 'FAILS'}")


def check(label, page, value):
    ok = value == page
    CHECKS.append(ok)
    print(f"  {label:<70s} page {page}  computed {value}  {'match' if ok else 'DIFFERS'}")


def op(P, k, N):
    out = np.ones((1, 1), dtype=P.dtype)
    for s in range(N):
        out = np.kron(out, P if s == k else np.eye(2, dtype=P.dtype))
    return out


def zval(i, k, N):
    """Z eigenvalue of site k (site 0 the leading tensor factor) in basis state i."""
    return 1 - 2 * ((i >> (N - 1 - k)) & 1)


def popcount(i):
    return bin(i).count("1")


class GaussH:
    """A Hermitian H = Hr + i Hi with integer Hr (symmetric) and Hi (antisymmetric)."""

    def __init__(self, Hr, Hi=None):
        self.r = Hr
        self.i = np.zeros_like(Hr) if Hi is None else Hi

    def __add__(self, other):
        return GaussH(self.r + other.r, self.i + other.i)


def xxz_chain(N, J, JZ, bonds=None):
    """sum over bonds of J (XX + YY) + JZ ZZ, integer. YY = -(YR x YR)."""
    d = 2 ** N
    H = np.zeros((d, d), dtype=np.int64)
    for k, l in bonds if bonds is not None else [(k, k + 1) for k in range(N - 1)]:
        H += J * (op(X, k, N) @ op(X, l, N) - op(YR, k, N) @ op(YR, l, N)) + JZ * op(Z, k, N) @ op(Z, l, N)
    return GaussH(H)


def two_l(H, a, N, a_im=None):
    """2L as (real part, imaginary part), integer, row-major vec: vec(A rho B) = kron(A, B^T) vec(rho).
    a is the real symmetric part of the Kossakowski matrix; a_im (antisymmetric, optional) its imaginary
    part, whose anticommutator terms cancel, leaving i sum a_im_kl Z_k rho Z_l."""
    d = 2 ** N
    Id = np.eye(d, dtype=np.int64)
    re = 2 * (np.kron(H.i, Id) - np.kron(Id, H.i.T))      # -i [i Hi, rho] = [Hi, rho]
    im = -2 * (np.kron(H.r, Id) - np.kron(Id, H.r.T))     # -i [Hr, rho]
    Zs = [op(Z, k, N) for k in range(N)]
    for k in range(N):
        for l in range(N):
            if a[k][l] == 0:
                continue
            ZlZk = Zs[l] @ Zs[k]
            re += int(a[k][l]) * (2 * np.kron(Zs[k], Zs[l].T) - np.kron(ZlZk, Id) - np.kron(Id, ZlZk.T))
    if a_im is not None:
        for k in range(N):
            for l in range(N):
                if a_im[k][l]:
                    im += 2 * int(a_im[k][l]) * np.kron(Zs[k], Zs[l].T)
    return re, im


def kossakowski_shared(N, eta):
    m = GAMMA * ((1 - eta) * np.eye(N, dtype=object) + eta * np.ones((N, N), dtype=object))
    assert all(Fraction(x).denominator == 1 for x in m.flat)
    return np.array([[int(x) for x in row] for row in m], dtype=np.int64)


def k_operator(a, N):
    Zs = [op(Z, k, N) for k in range(N)]
    return sum(int(a[k][l]) * Zs[l] @ Zs[k] for k in range(N) for l in range(N))


def dark_cells(a, N):
    """The dark cells: the matrix units |i><j| whose dissipator rate is zero, read off the diagonal of 2L
    with H = 0."""
    d = 2 ** N
    re, _ = two_l(GaussH(np.zeros((d, d), dtype=np.int64)), a, N)
    diag = np.diag(re)
    return {(i, j) for i in range(d) for j in range(d) if diag[i * d + j] == 0}


def dark_classes(a, N):
    """The configurations grouped by a z(i), each class sorted, in order of first appearance."""
    groups = {}
    for i in range(2 ** N):
        key = tuple(int(x) for x in a @ np.array([zval(i, k, N) for k in range(N)], dtype=np.int64))
        groups.setdefault(key, []).append(i)
    return list(groups.values())


_ZM = {}


def k_on_configs(c, diag, N):
    """K(i) = z(i)^T a z(i) with a = c off the diagonal and diag on it, on every configuration i, integers."""
    Zm = _ZM.get(N)
    if Zm is None:
        Zm = _ZM[N] = np.array([[zval(i, k, N) for k in range(N)] for i in range(2 ** N)], dtype=np.int64)
    a = np.array(c, dtype=np.int64).copy()
    for k in range(N):
        a[k, k] = diag[k]
    return np.einsum("ik,kl,il->i", Zm, a, Zm)


def qualifying_popcounts(c, diag, N):
    kv = k_on_configs(c, diag, N)
    out = []
    for p in range(N + 1):
        vals = {int(kv[i]) for i in range(2 ** N) if popcount(i) == p}
        if len(vals) == 1:
            out.append(p)
    return out


def three_parts(c, N):
    """C_kl = kbar + x_k + x_l + R_kl, Fractions; N >= 3."""
    T = sum(Fraction(int(c[k][l])) for k in range(N) for l in range(k + 1, N))
    kbar = T / comb(N, 2)
    r = [sum(Fraction(int(c[k][l])) for l in range(N) if l != k) for k in range(N)]
    x = [(r[k] - (N - 1) * kbar) / (N - 2) for k in range(N)]
    Rm = [[(Fraction(int(c[k][l])) - kbar - x[k] - x[l]) if k != l else Fraction(0) for l in range(N)] for k in range(N)]
    return T, kbar, x, Rm


def predicted_popcounts(c, N):
    if N == 2:
        return [0, 1, 2]
    _, _, x, Rm = three_parts(c, N)
    x0 = all(v == 0 for v in x)
    r0 = all(v == 0 for row in Rm for v in row)
    out = []
    for p in range(N + 1):
        if p in (0, N) or (x0 and r0) or (p in (1, N - 1) and x0) or (2 * p == N and 2 <= p <= N - 2 and r0):
            out.append(p)
    return out


def pi_perm(N):
    """Pi_Z(|i><j|) = |i><j|^T X^N = |j><i_bar|: the index map, a pure permutation of matrix units."""
    d = 2 ** N
    perm = np.empty(d * d, dtype=np.int64)
    for i, j in itertools.product(range(d), range(d)):
        perm[i * d + j] = j * d + (d - 1 - i)
    return perm


def conj_by_pi(M, perm):
    out = np.empty_like(M)
    out[np.ix_(perm, perm)] = M
    return out


def conj_by_signed(M, perm, sign):
    """(P S) M (P S)^-1 for the signed permutation |u> -> sign[u] |perm[u]> of matrix units."""
    out = np.empty_like(M)
    out[np.ix_(perm, perm)] = M * np.outer(sign, sign)
    return out


def blocks(N):
    d = 2 ** N
    idx = {}
    for i, j in itertools.product(range(d), range(d)):
        idx.setdefault((popcount(i), popcount(j)), []).append(i * d + j)
    return {key: np.array(v) for key, v in idx.items()}


def centre(N, eta, p, q):
    """c(p, q) = N (1 - eta) gamma + 2 gamma eta [(p - N/2)^2 + (q - N/2)^2], exact."""
    h = Fraction(N, 2)
    return N * (1 - eta) * GAMMA + 2 * GAMMA * eta * ((p - h) ** 2 + (q - h) ** 2)


def random_psd(N, rng):
    """A positive definite integer Kossakowski matrix with EVERY off-diagonal entry nonzero, so that no
    draw degenerates into local dephasing (where every identity below collapses onto F1)."""
    while True:
        B = rng.integers(-2, 3, size=(N, N))
        a = (B.T @ B + np.eye(N, dtype=np.int64)).astype(np.int64)
        if all(a[k, l] != 0 for k in range(N) for l in range(N) if k != l):
            return a


def random_psd_with_kernel(N, rng):
    """a = B^T B with every row of B orthogonal to one random v in {0, +-1}^N of support >= 2 that is not
    all-equal, so ker a = span(v) (checked by an exact rank) and the cells with u = +-2 v are dark:
    unlike {0}, the all-ones direction and its complement, such a kernel can single out particular sites
    from N = 3 on (at N = 2 the only admissible v is +-(1, -1), the sum-zero line), so a site-order
    mismatch between the rate and the kernel would show where it does."""
    while True:
        v = rng.integers(-1, 2, size=N).astype(np.int64)
        supp = [k for k in range(N) if v[k] != 0]
        if len(supp) < 2 or (len(supp) == N and len(set(v.tolist())) == 1):
            continue
        s0 = supp[0]
        gens = [np.eye(N, dtype=np.int64)[k] for k in range(N) if v[k] == 0]
        for l in supp[1:]:
            g = np.zeros(N, dtype=np.int64)
            g[l], g[s0] = v[s0], -v[l]
            gens.append(g)
        G = np.array(gens, dtype=np.int64)
        B = rng.integers(-2, 3, size=(N - 1, N - 1)).astype(np.int64) @ G
        a = (B.T @ B).astype(np.int64)
        if not (a @ v).any() and SymMatrix(a.tolist()).rank() == N - 1:
            return a, v


def random_xn_even(N, rng):
    d = 2 ** N
    M = rng.integers(-3, 4, size=(d, d))
    M = M + M.T
    XN = np.ones((1, 1), dtype=np.int64)
    for _ in range(N):
        XN = np.kron(XN, X)
    return GaussH((M + XN @ M @ XN).astype(np.int64))


Y_SYM = symbols("y")


def random_symmetric(N, rng):
    """A real symmetric integer matrix with entries of either sign, not necessarily positive semidefinite."""
    B = rng.integers(-4, 5, size=(N, N))
    return (B + B.T).astype(np.int64)


def random_hermitian(N, rng):
    """A traceless Hermitian H with integer real and imaginary parts and no symmetry at all."""
    d = 2 ** N
    M = rng.integers(-3, 4, size=(d, d))
    A = rng.integers(-3, 4, size=(d, d))
    Hr = (M + M.T).astype(np.int64)
    Hr[0, 0] -= int(np.trace(Hr))
    return GaussH(Hr, (A - A.T).astype(np.int64))


def power_trace(re, im, shift, n):
    """Tr((2L + shift)^n) for n = 3 or 5 as (real, imaginary) Python integers, exact (object arithmetic)."""
    assert n in (3, 5)
    size = re.shape[0]
    R = re.astype(object) + shift * np.eye(size, dtype=np.int64).astype(object)
    Im = im.astype(object)
    R2, I2 = R.dot(R) - Im.dot(Im), R.dot(Im) + Im.dot(R)
    if n == 3:
        return int(np.sum(R2 * R.T) - np.sum(I2 * Im.T)), int(np.sum(R2 * Im.T) + np.sum(I2 * R.T))
    R4, I4 = R2.dot(R2) - I2.dot(I2), R2.dot(I2) + I2.dot(R2)
    return int(np.sum(R4 * R.T) - np.sum(I4 * Im.T)), int(np.sum(R4 * Im.T) + np.sum(I4 * R.T))


def noise_part(a):
    """4 tr(C (2a + Lambda) C), C the covariances, Lambda the diagonal: the noise term of the third moment's bracket."""
    N = len(a)
    A = [[int(a[k][l]) for l in range(N)] for k in range(N)]
    C = [[A[k][l] if k != l else 0 for l in range(N)] for k in range(N)]
    W = [[2 * A[l][m] + (A[l][l] if l == m else 0) for m in range(N)] for l in range(N)]
    return 4 * sum(C[k][l] * W[l][m] * C[m][k] for k in range(N) for l in range(N) for m in range(N))


def third_moment_closed(H, a, N):
    """Tr((2L + 2 sigma)^3) = 8 Tr((L + sigma)^3) from the closed form -4^N [4 tr(C(2a + Lambda)C) - 3 tr(C Theta)
    - 6 xi^T a xi], with Theta_kl = Tr(H^2 Z_k Z_l)/2^N and xi_k = Tr(H Z_k)/2^N, H traceless; in integers."""
    d = 2 ** N
    assert int(np.trace(H.r)) == 0, "the closed form is written for a traceless H"
    Zs = [op(Z, k, N) for k in range(N)]
    Hsq = H.r @ H.r - H.i @ H.i                          # the real part of H^2; its imaginary part has a zero diagonal
    T2 = [[int(np.trace(Hsq @ Zs[k] @ Zs[l])) for l in range(N)] for k in range(N)]
    T1 = [int(np.trace(H.r @ Zs[k])) for k in range(N)]
    h_two = sum(int(a[k][l]) * T2[k][l] for k in range(N) for l in range(N) if k != l)
    h_one = sum(int(a[k][l]) * T1[k] * T1[l] for k in range(N) for l in range(N))
    return -(8 * d * d * noise_part(a) - 24 * d * h_two - 48 * h_one)


def charpoly(re, im):
    n = re.shape[0]
    rows = [[QQ_I(int(re[r, c]), int(im[r, c])) for c in range(n)] for r in range(n)]
    return Poly(DomainMatrix(rows, (n, n), QQ_I).charpoly(), Y_SYM, domain=QQ_I)


def pairs_hold(re, im, idx, shift_of, N):
    """The blocks A = (p, q) whose partner B = (q, N - p) has p_B(y) == (-1)^n p_A(s - y), s = shift_of(p, q)
    the pair sum in 2L units (a sympy number); a list of the A that pair."""
    held = []
    cache = {}
    for (p, q), ix in idx.items():
        B = (q, N - p)
        for key in ((p, q), B):
            if key not in cache:
                jx = idx[key]
                cache[key] = charpoly(re[np.ix_(jx, jx)], im[np.ix_(jx, jx)])
        pA, pB = cache[(p, q)], cache[B]
        n = len(ix)
        s = shift_of(p, q)
        rhs = (-1) ** n * pA.compose(Poly(s - Y_SYM, Y_SYM, domain=QQ_I))
        if pB == rhs:
            held.append((p, q))
    return held


def trace_shift(re, im, idx, N):
    """The only pair sum two blocks can have, (tr A + tr B)/dim, in 2L units, as a function of A = (p, q)."""
    def s(p, q):
        ix, jx = idx[(p, q)], idx[(q, N - p)]
        n = len(ix)
        tre = int(np.trace(re[np.ix_(ix, ix)])) + int(np.trace(re[np.ix_(jx, jx)]))
        tim = int(np.trace(im[np.ix_(ix, ix)])) + int(np.trace(im[np.ix_(jx, jx)]))
        return Rational(tre, n) + Rational(tim, n) * SYM_I
    return s


def comm_zero(A, B):
    """[A, B] == 0 exactly, sparse integer products."""
    A, B = csr_matrix(A), csr_matrix(B)
    return (A @ B - B @ A).count_nonzero() == 0


def comm(A, B):
    A, B = csr_matrix(A), csr_matrix(B)
    return (A @ B - B @ A).toarray()


def main():
    rng = np.random.default_rng(20261010)

    print("A. The rate law: D_a is diagonal on |i><j| with entry -1/2 d^T a d (2L units: -d^T a d)")
    for N in (2, 3, 4, 5):
        d = 2 ** N
        for name, a in (("random PSD a", random_psd(N, rng)), ("shared field, eta 3/10", kossakowski_shared(N, ETAS[0]))):
            re, _ = two_l(GaussH(np.zeros((d, d), dtype=np.int64)), a, N)
            off = re - np.diag(np.diag(re))
            want = np.empty(d * d, dtype=np.int64)
            for i, j in itertools.product(range(d), range(d)):
                dv = np.array([zval(i, k, N) - zval(j, k, N) for k in range(N)], dtype=np.int64)
                want[i * d + j] = -int(dv @ a @ dv)
            claim(f"N={N} {name}: off-diagonal entries 0, diagonal = -u^T a u", not off.any() and np.array_equal(np.diag(re), want))
        for eta in ETAS:
            a = kossakowski_shared(N, eta)
            re, _ = two_l(GaussH(np.zeros((d, d), dtype=np.int64)), a, N)
            ok = True
            for i, j in itertools.product(range(d), range(d)):
                h, p, q = popcount(i ^ j), popcount(i), popcount(j)
                ok &= Fraction(int(re[i * d + j, i * d + j]), 2) == -2 * GAMMA * ((1 - eta) * h + eta * (p - q) ** 2)
            claim(f"N={N} eta {eta}: rate = 2 gamma [(1 - eta) |D| + eta (p - q)^2] on every |i><j|", ok)
    # the price-pair page's law: the coherence aligned on D and the one aligned on D_bar (Pi_Z partners) pay
    # 2 sigma + 4 sum_{k<l same side} a_kl together, its P + 2 sum Cov_ij [same side] with Cov_ij = 2 a_kl
    for N in (3, 4):
        d = 2 ** N
        a = random_psd(N, rng)
        re, _ = two_l(GaussH(np.zeros((d, d), dtype=np.int64)), a, N)
        ok = True
        for D in range(1, d - 1):
            i, j = 0, D                                    # |0...0><1_D 0_Dbar|, aligned on D
            jb, ib = D, d - 1                              # its Pi_Z image |j><i_bar| = |1_D 0_Dbar><1...1|
            pair = Fraction(-int(re[i * d + j, i * d + j]) - int(re[jb * d + ib, jb * d + ib]), 2)
            side = [(D >> (N - 1 - k)) & 1 for k in range(N)]
            same = sum(int(a[k, l]) for k in range(N) for l in range(k + 1, N) if side[k] == side[l])
            ok &= pair == 2 * int(np.trace(a)) + 4 * same
        claim(f"N={N} random PSD a: aligned D and aligned D_bar pay 2 sigma + 4 sum_(same side) a_kl (price-pair law)", ok)
    a = kossakowski_shared(3, ETAS[0])
    re, _ = two_l(GaussH(np.zeros((8, 8), dtype=np.int64)), a, 3)
    loc = all(Fraction(int(re[r, r]), 2) == -2 * GAMMA * popcount((r // 8) ^ (r % 8)) for r in range(64))
    claim("control, N=3 eta 3/10: the local rate 2 gamma |D| does NOT describe the shared field", not loc)

    print("\nB. The block law: L_eta = L_local((1 - eta) gamma) - 2 gamma eta (p - q)^2 on block (p, q)")
    for N in (2, 3, 4, 5):
        d = 2 ** N
        H = xxz_chain(N, 5, 2)                      # J = 5, Delta = 2/5
        idx = blocks(N)
        allb = np.arange(d * d)
        for eta in ETAS:
            a = kossakowski_shared(N, eta)
            aloc = np.diag([int(GAMMA * (1 - eta))] * N).astype(np.int64)
            re, im = two_l(H, a, N)
            rl, il = two_l(H, aloc, N)
            ok_law = ok_stay = ok_dose = ok_unit = True
            for (p, q), ix in idx.items():
                shift = int(-4 * GAMMA * eta * (p - q) ** 2)            # 2L units
                ok_law &= np.array_equal(re[np.ix_(ix, ix)] - rl[np.ix_(ix, ix)], shift * np.eye(len(ix), dtype=np.int64))
                ok_law &= np.array_equal(im[np.ix_(ix, ix)], il[np.ix_(ix, ix)])
                rest = np.setdiff1d(allb, ix)
                ok_stay &= not re[np.ix_(rest, ix)].any() and not im[np.ix_(rest, ix)].any()
                if p == q:
                    ok_dose &= np.array_equal(re[np.ix_(ix, ix)], rl[np.ix_(ix, ix)])
                    if eta == 1:
                        ok_unit &= not re[np.ix_(ix, ix)].any()
            claim(f"N={N} eta {eta}: block law entry by entry, nothing leaves a block", ok_law and ok_stay)
            claim(f"N={N} eta {eta}: every (p, p) block is the local one at (1 - eta) gamma", ok_dose)
            if eta == 1:
                claim(f"N={N} eta 1: every (p, p) block is unitary (dissipative part 0)", ok_unit)
    N = 4
    a = random_psd(N, rng)
    re, _ = two_l(xxz_chain(N, 5, 2), a, N)
    rl, _ = two_l(xxz_chain(N, 5, 2), np.diag(np.diag(a)), N)
    scalar = all(len(set(np.diag(re[np.ix_(ix, ix)] - rl[np.ix_(ix, ix)]))) == 1 for ix in blocks(N).values())
    claim("control, N=4 random PSD a: the correlated part is NOT one constant per block", not scalar)

    print("\nB'. One common covariance at site-dependent rates: a = diag(g) + c (11^T - 1)")
    for N in (3, 4, 5):
        g = [7, 12, 9, 15, 10][:N]
        cc = 3
        a = (np.diag(g) + cc * (np.ones((N, N), dtype=np.int64) - np.eye(N, dtype=np.int64))).astype(np.int64)
        aloc = np.diag([x - cc for x in g]).astype(np.int64)
        H = xxz_chain(N, 5, 2)
        re, im = two_l(H, a, N)
        rl, il = two_l(H, aloc, N)
        idx = blocks(N)
        perm = pi_perm(N)
        cre, cim = conj_by_pi(re, perm), conj_by_pi(im, perm)
        ok_law = ok_pair = True
        for (p, q), ix in idx.items():
            ok_law &= np.array_equal(re[np.ix_(ix, ix)] - rl[np.ix_(ix, ix)], -4 * cc * (p - q) ** 2 * np.eye(len(ix), dtype=np.int64))
            ok_law &= np.array_equal(im[np.ix_(ix, ix)], il[np.ix_(ix, ix)])
            jx = idx[(q, N - p)]
            cent = (sum(g) - N * cc) + cc * ((p - q) ** 2 + (p + q - N) ** 2)
            ok_pair &= np.array_equal(cre[np.ix_(jx, jx)], -re[np.ix_(jx, jx)] - 4 * cent * np.eye(len(jx), dtype=np.int64))
            ok_pair &= np.array_equal(cim[np.ix_(jx, jx)], -im[np.ix_(jx, jx)])
        claim(f"N={N} rates {g}, covariance {cc}: local at g - c, shifted by -2c(p - q)^2 on every block", ok_law)
        claim(f"N={N} rates {g}, covariance {cc}: Pi_Z pairs the blocks about (sigma - N c) + c (d^2 + m^2)", ok_pair)
        if N == 3:
            held = pairs_hold(re, im, idx, lambda p, q: Rational(-4 * ((sum(g) - N * cc) + cc * ((p - q) ** 2 + (p + q - N) ** 2))), N)
            claim("N=3 rates [7, 12, 9], covariance 3: characteristic polynomials pair about the same centres", len(held) == len(idx), f"({len(held)}/{len(idx)})")
    # DM bonds on the open chain: X^N-odd, so Pi_Z is not the route; a rotation of each site about z, a diagonal
    # unitary, turns J (XX + YY) + D (XY - YX) into hopping and commutes with every Z_k, so the blocks still pair
    for N in (3, 4):
        g = [7, 12, 9, 15][:N]
        cc = 3
        a = (np.diag(g) + cc * (np.ones((N, N), dtype=np.int64) - np.eye(N, dtype=np.int64))).astype(np.int64)
        d = 2 ** N
        dmi = sum(3 * (op(X, k, N) @ op(YR, k + 1, N) - op(YR, k, N) @ op(X, k + 1, N)) for k in range(N - 1))   # XY - YX = i * dmi
        Hd = xxz_chain(N, 5, 2) + GaussH(np.zeros((d, d), dtype=np.int64), dmi)
        XN = np.ones((1, 1), dtype=np.int64)
        for _ in range(N):
            XN = np.kron(XN, X)
        odd = not np.array_equal(XN @ Hd.i @ XN, Hd.i)
        re, im = two_l(Hd, a, N)
        idx = blocks(N)
        held = pairs_hold(re, im, idx, lambda p, q: Rational(-4 * ((sum(g) - N * cc) + cc * ((p - q) ** 2 + (p + q - N) ** 2))), N)
        claim(f"N={N} rates {g}, covariance {cc}, open chain with DM bonds 3 (XY - YX) (odd under X^N): every block pairs "
              f"with its image about the same centres", odd and len(held) == len(idx), f"({len(held)}/{len(idx)})")
        if N == 3:
            # the control through the same door: on the triangle the hopping phases multiply to a flux that no
            # diagonal unitary removes, and blocks stop pairing
            ring = [(0, 1), (1, 2), (2, 0)]
            dmr = sum(3 * (op(X, k, N) @ op(YR, l, N) - op(YR, k, N) @ op(X, l, N)) for k, l in ring)
            re, im = two_l(xxz_chain(N, 5, 2, bonds=ring) + GaussH(np.zeros((d, d), dtype=np.int64), dmr), a, N)
            held = pairs_hold(re, im, idx, lambda p, q: Rational(-4 * ((sum(g) - N * cc) + cc * ((p - q) ** 2 + (p + q - N) ** 2))), N)
            claim(f"N=3 the same on the triangle with DM bonds of one orientation (the control): 12 of 16 blocks pair",
                  len(held) == 12, f"({len(held)}/{len(idx)})")

    print("\nB'''. The edge blocks: a coherence that disagrees at one site l pays 2 a_ll whatever the covariances")
    for N in (3, 4, 5):
        a = random_psd(N, rng)
        aloc = np.diag(np.diag(a)).astype(np.int64)
        H = xxz_chain(N, 5, 2)
        re, im = two_l(H, a, N)
        rl, il = two_l(H, aloc, N)
        idx = blocks(N)
        edges = [(0, 1), (1, 0), (N - 1, N), (N, N - 1)]
        same = all(np.array_equal(re[np.ix_(idx[e], idx[e])], rl[np.ix_(idx[e], idx[e])])
                   and np.array_equal(im[np.ix_(idx[e], idx[e])], il[np.ix_(idx[e], idx[e])]) for e in edges)
        img = idx[(1, N)]
        differs = not np.array_equal(re[np.ix_(img, img)], rl[np.ix_(img, img)])
        claim(f"N={N} random PSD a: the four edge blocks run exactly as under local noise, and the image (1, N) of (0, 1) does not",
              same and differs)

    print("\nB''. A general covariance pattern: a block pairs with its image when K is one number on both its popcount sectors")
    # sufficient by the price identity: Pi_Z's residual {K, .} is then one number on the image; on popcount 1 and
    # N - 1 that means equal covariance row sums, which at N = 3 forces one common covariance
    patterns = {3: [("one covariance on the bond (1, 2)", [[10, 0, 0], [0, 10, 3], [0, 3, 10]]),
                    ("one covariance end to end", [[10, 0, 3], [0, 10, 0], [3, 0, 10]]),
                    ("two different covariances", [[10, 3, 0], [3, 10, -3], [0, -3, 10]]),
                    ("three different covariances", [[10, 3, 1], [3, 12, 2], [1, 2, 9]])],
                4: [("a decaying pattern", [[10, 3, 1, 0], [3, 10, 3, 1], [1, 3, 10, 3], [0, 1, 3, 10]]),
                    ("one covariance on the bond (0, 1)", [[10, 3, 0, 0], [3, 10, 0, 0], [0, 0, 10, 0], [0, 0, 0, 10]]),
                    ("one covariance on the middle bond (1, 2)", [[10, 0, 0, 0], [0, 10, 3, 0], [0, 3, 10, 0], [0, 0, 0, 10]]),
                    ("two disjoint shared pairs (0, 1), (2, 3)", [[10, 3, 0, 0], [3, 10, 0, 0], [0, 0, 10, 3], [0, 0, 3, 10]]),
                    ("a ring of nearest-neighbour covariances", [[10, 3, 0, 3], [3, 10, 3, 0], [0, 3, 10, 3], [3, 0, 3, 10]]),
                    ("a perfect-matching pattern (3, 1, -2)", [[10, 3, 1, -2], [3, 10, -2, 1], [1, -2, 10, 3], [-2, 1, 3, 10]]),
                    ("the sign pattern 3 s_k s_l, s = (1, 1, 1, -1)", [[10, 3, 3, -3], [3, 10, 3, -3], [3, 3, 10, -3], [-3, -3, -3, 10]]),
                    ("the sign pattern 3 s_k s_l, s = (1, -1, 1, -1)", [[10, -3, 3, -3], [-3, 10, -3, 3], [3, -3, 10, -3], [-3, 3, -3, 10]])]}
    patterns[3].append(("the sign pattern 3 s_k s_l, s = (1, -1, 1)", [[10, -3, 3], [-3, 10, -3], [3, -3, 10]]))
    for N, cases in patterns.items():
        idx = blocks(N)
        d = 2 ** N
        for name, a in cases:
            a = np.array(a, dtype=np.int64)
            Kd = np.diag(k_operator(a, N))
            const = {c for c in range(N + 1) if len({int(Kd[i]) for i in range(d) if popcount(i) == c}) == 1}
            want = {(p, q) for (p, q) in idx if p in const and q in const}
            re, im = two_l(xxz_chain(N, 5, 2), a, N)
            held = pairs_hold(re, im, idx, trace_shift(re, im, idx, N), N)
            claim(f"N={N} {name}: the blocks that pair are exactly those whose two popcount sectors see one value of K "
                  f"(K constant on popcounts {sorted(const)})", set(held) == want, f"({len(held)}/{len(idx)})")
            # and they pair about the stated centre -(K_p + K_q)/2 (2L units: pair sum -2 (K_p + K_q))
            kval = {c: int(Kd[next(i for i in range(d) if popcount(i) == c)]) for c in const}
            sub = {key: ix for key, ix in idx.items() if key in want}
            held_k = pairs_hold(re, im, sub, lambda p, q: Rational(-2 * (kval[p] + kval[q])), N)
            claim(f"N={N} {name}: every qualifying block pairs with its image about -(K_p + K_q)/2", set(held_k) == want)

    print("\nC. The price identity: Pi_Z L Pi_Z^-1 + L + {K, .} = 0 (2L units: + 2{K, .})")
    for N in (2, 3, 4, 5):
        d = 2 ** N
        Id = np.eye(d, dtype=np.int64)
        perm = pi_perm(N)
        a = random_psd(N, rng)
        K = k_operator(a, N)
        antiK = np.kron(K, Id) + np.kron(Id, K.T)
        for name, H in (("XXZ chain", xxz_chain(N, 5, 2)), ("random real X^N-even H", random_xn_even(N, rng))):
            re, im = two_l(H, a, N)
            r_re = conj_by_pi(re, perm) + re + 2 * antiK
            r_im = conj_by_pi(im, perm) + im
            claim(f"N={N} random PSD a, {name}: residual exactly 0", not r_re.any() and not r_im.any())
        while True:                                                   # an indefinite symmetric a, every entry nonzero
            B = rng.integers(-4, 5, size=(N, N))
            ai = (B + B.T).astype(np.int64)
            if (ai != 0).all() and np.linalg.eigvalsh(ai.astype(float)).min() < 0:
                break
        Ki = k_operator(ai, N)
        re, im = two_l(xxz_chain(N, 5, 2), ai, N)
        r_re = conj_by_pi(re, perm) + re + 2 * (np.kron(Ki, Id) + np.kron(Id, Ki.T))
        claim(f"N={N} indefinite symmetric a (negative entries on the diagonal allowed): residual exactly 0",
              not r_re.any() and not (conj_by_pi(im, perm) + im).any())
        Hf = xxz_chain(N, 5, 2) + GaussH(3 * op(Z, 0, N))
        re, im = two_l(Hf, a, N)
        bad = (conj_by_pi(re, perm) + re + 2 * antiK).any() or (conj_by_pi(im, perm) + im).any()
        claim(f"control N={N}: a Z field on site 0 ([H, X^N] != 0) breaks it", bad)
        if N >= 2:
            Hy = xxz_chain(N, 5, 2) + GaussH(np.zeros((d, d), dtype=np.int64), 3 * op(YR, 0, N) @ op(Z, 1, N))
            re, im = two_l(Hy, a, N)
            bad = (conj_by_pi(re, perm) + re + 2 * antiK).any() or (conj_by_pi(im, perm) + im).any()
            claim(f"control N={N}: the imaginary X^N-even term 3 Y_0 Z_1 (H not real) breaks it", bad)
            # F108's Pi_5bilinear is Pi_Z after conjugation by Y^N; that conjugation flips every Z_k, so it leaves the
            # correlated dissipator alone, and the composite is the signed permutation |i><j| -> s_i s_j |j_bar><i|,
            # Y^N |i> = (i^N) s_i |i_bar>; it carries L with the imaginary bond to -L - {K, .} exactly
            YRN = np.ones((1, 1), dtype=np.int64)
            for _ in range(N):
                YRN = np.kron(YRN, YR)
            sgn_site = [int(YRN[d - 1 - i, i]) for i in range(d)]
            perm5 = np.empty(d * d, dtype=np.int64)
            sign5 = np.empty(d * d, dtype=np.int64)
            for i, j in itertools.product(range(d), range(d)):
                perm5[i * d + j] = (d - 1 - j) * d + i
                sign5[i * d + j] = sgn_site[i] * sgn_site[j]
            ok5 = not (conj_by_signed(re, perm5, sign5) + re + 2 * antiK).any() and not (conj_by_signed(im, perm5, sign5) + im).any()
            claim(f"N={N}: F108's Pi_5bilinear (Pi_Z after conjugation by Y^N) carries L with that term to -L - {{K, .}} exactly", ok5)
            # the residuals add in squares: M = M_H + M_D, M_D real and diagonal, M_H's diagonal imaginary
            rh, ih = two_l(Hy, np.zeros((N, N), dtype=np.int64), N)
            mh_re, mh_im = conj_by_pi(rh, perm) + rh, conj_by_pi(ih, perm) + ih
            md = -2 * (antiK - 2 * int(np.trace(a)) * np.eye(d * d, dtype=np.int64))
            m_re, m_im = conj_by_pi(re, perm) + re + 4 * int(np.trace(a)) * np.eye(d * d, dtype=np.int64), conj_by_pi(im, perm) + im
            sq = lambda A: int((A.astype(object) ** 2).sum())
            claim(f"N={N}: with that term, ||M||^2 = ||M_H||^2 + ||M_D||^2 exactly (M_H != 0 with no real diagonal, "
                  f"M_D = -{{K - sigma, .}} real and diagonal)",
                  np.array_equal(m_re, mh_re + md) and np.array_equal(m_im, mh_im) and sq(mh_re) + sq(mh_im) > 0
                  and not np.diag(mh_re).any() and sq(m_re) + sq(m_im) == sq(mh_re) + sq(mh_im) + sq(md))
        aloc = np.diag(np.diag(a))
        Kloc = k_operator(aloc, N)
        claim(f"N={N}: K = tr(a) 1 for the diagonal part of a (F1's 2 sigma), K != scalar for a itself",
              np.array_equal(Kloc, int(np.trace(a)) * Id) and len(set(np.diag(K))) > 1)
        if N >= 2:
            # Absorption section 2's pair sqrt(3)(Z_0 + i Z_1), (Z_0 - i Z_1), with local Z at rate 4 on the other sites:
            # Kossakowski diag 4, off-diagonal -2i on (0, 1); K = 4N 1 is a number, the identity still fails
            ar = np.diag([4] * N).astype(np.int64)
            ai = np.zeros((N, N), dtype=np.int64)
            ai[0, 1], ai[1, 0] = -2, 2
            re, im = two_l(xxz_chain(N, 5, 2), ar, N, ai)
            re0, im0 = two_l(xxz_chain(N, 5, 2), ar, N)
            T = im - im0                                                   # the imaginary Kossakowski part, i T
            Kc = k_operator(ar, N)
            antiKc = np.kron(Kc, Id) + np.kron(Id, Kc.T)
            r_re, r_im = conj_by_pi(re, perm) + re + 2 * antiKc, conj_by_pi(im, perm) + im
            claim(f"control N={N}: a Hermitian non-real a (the pair (Z_0 +- i Z_1)) has K = {4 * N} 1 and still breaks it: "
                  "Pi_Z fixes its imaginary part T, so the residual is exactly 2T",
                  T.any() and np.array_equal(conj_by_pi(T, perm), T) and not r_re.any() and np.array_equal(r_im, 2 * T))

    print("\nD. The residual: ||M_D||_F^2 = 8 4^N sum_{k<l} a_kl^2")
    for N in (2, 3, 4, 5):
        d = 2 ** N
        perm = pi_perm(N)
        for name, a in (("random PSD a", random_psd(N, rng)), ("shared field, eta 3/10", kossakowski_shared(N, ETAS[0]))):
            re, _ = two_l(GaussH(np.zeros((d, d), dtype=np.int64)), a, N)
            M2 = conj_by_pi(re, perm) + re + 4 * int(np.trace(a)) * np.eye(d * d, dtype=np.int64)   # 2 M_D
            norm2 = Fraction(int((M2.astype(object) ** 2).sum()), 4)
            pred = 8 * 4 ** N * sum(int(a[k, l]) ** 2 for k in range(N) for l in range(k + 1, N))
            check(f"N={N} {name}: ||M_D||_F^2", pred, norm2)
            claim(f"N={N} {name}: M_D is Pi_Z-invariant (Pi_Z M_D Pi_Z^-1 = M_D) and nonzero",
                  np.array_equal(conj_by_pi(M2, perm), M2) and M2.any())
        eta = ETAS[0]
        a = kossakowski_shared(N, eta)
        check(f"N={N} shared field: 4^(N+1) N(N-1) gamma^2 eta^2", 4 ** (N + 1) * N * (N - 1) * GAMMA ** 2 * eta ** 2,
              Fraction(8 * 4 ** N * sum(int(a[k, l]) ** 2 for k in range(N) for l in range(k + 1, N))))
    # the norm convention against the registry's T1 closed form ||M(T1)||_F^2 = 4^(N-1) [3 sum g^2 + 4 (sum g)^2]
    # (ANALYTICAL_FORMULAS, F1's T1 residual, at sigma = 0 for pure T1), same Pi_Z, same vec, same Frobenius norm
    sm = np.array([[0, 1], [0, 0]], dtype=np.int64)                  # sigma^- = |0><1|
    for N in (2, 3, 4):
        d = 2 ** N
        Id = np.eye(d, dtype=np.int64)
        g = [3, 5, 2, 7][:N]
        re = np.zeros((d * d, d * d), dtype=np.int64)
        for k in range(N):
            c = op(sm, k, N)
            n = c.T @ c
            re += g[k] * (2 * np.kron(c, c) - np.kron(n, Id) - np.kron(Id, n.T))   # 2 L_D, c rho c^dag = kron(c, conj c)
        M2 = conj_by_pi(re, pi_perm(N)) + re
        check(f"N={N} T1 rates {g}: ||M(T1)||_F^2 against the registry's closed form",
              Fraction(4 ** (N - 1) * (3 * sum(x * x for x in g) + 4 * sum(g) ** 2)), Fraction(int((M2.astype(object) ** 2).sum()), 4))

    print("\nE. The blockwise mirror")
    for N in (2, 3, 4, 5):
        d = 2 ** N
        perm = pi_perm(N)
        idx = blocks(N)
        H = xxz_chain(N, 5, 2)
        for eta in ETAS:
            a = kossakowski_shared(N, eta)
            re, im = two_l(H, a, N)
            cre, cim = conj_by_pi(re, perm), conj_by_pi(im, perm)
            ok = True
            for (p, q), ix in idx.items():
                jx = idx[(q, N - p)]
                c2 = 4 * centre(N, eta, p, q)                              # 2 * (2c) in 2L units
                assert c2.denominator == 1
                eye = np.eye(len(jx), dtype=np.int64)
                ok &= np.array_equal(cre[np.ix_(jx, jx)], -re[np.ix_(jx, jx)] - int(c2) * eye)
                ok &= np.array_equal(cim[np.ix_(jx, jx)], -im[np.ix_(jx, jx)])
            claim(f"N={N} eta {eta}: Pi_Z maps block (p, q) onto -L - 2c(p, q) on block (q, N - p)", ok)
    for N in (2, 3, 4, 5):
        d = 2 ** N
        fold = np.empty(d * d, dtype=np.int64)
        for i, j in itertools.product(range(d), range(d)):
            fold[i * d + j] = i * d + (d - 1 - j)                      # |i><j| -> |i><j_bar|, F89d's P
        odd = lambda x: sum((x >> (N - 1 - l)) & 1 for l in range(1, N, 2))
        gauge = np.array([(-1) ** (odd(i) + odd(j)) for i, j in itertools.product(range(d), range(d))], dtype=np.int64)
        idx = blocks(N)

        def fold_holds(re_t, im_t, re, im, eta):
            ok = True
            for (p, q) in idx:
                jx = idx[(p, N - q)]
                eye = np.eye(len(jx), dtype=np.int64)
                ok &= np.array_equal(re_t[np.ix_(jx, jx)], -re[np.ix_(jx, jx)] - int(4 * centre(N, eta, p, q)) * eye)
                ok &= np.array_equal(im_t[np.ix_(jx, jx)], -im[np.ix_(jx, jx)])
            return ok

        for eta in (ETAS[0], ETAS[2]):
            for jz in (0, 2, 5, 15):                                    # Delta = 0, 2/5, 1, 3 at J = 5
                re, im = two_l(xxz_chain(N, 5, jz), kossakowski_shared(N, eta), N)
                ok = fold_holds(conj_by_pi(re, fold), conj_by_pi(-im, fold), re, im, eta)     # P conj(L) P^T
                claim(f"N={N} eta {eta} Delta {Fraction(jz, 5)}: the antiunitary bra fold (F89d) maps (p, q) onto -conj - 2c(p, q) on (p, N - q)", ok)
            # the linear leg of the fold lattice is P with the staggered gauge D rho D, D = prod over odd sites of Z
            # (PROOF_CODIM1_BY_ADDITIVITY section 7): exact at Delta = 0, and it must break at Delta != 0
            for jz, want in ((0, True), (2, False)):
                re, im = two_l(xxz_chain(N, 5, jz), kossakowski_shared(N, eta), N)
                ok = fold_holds(conj_by_signed(re, fold, gauge), conj_by_signed(im, fold, gauge), re, im, eta)
                claim(f"N={N} eta {eta}: the linear gauged bra fold is {'exact' if want else 'broken'} at Delta = {Fraction(jz, 5)}", ok == want)
    for N in (3, 4):
        H = xxz_chain(N, 5, 2)
        idx = blocks(N)
        for eta in (ETAS[0], ETAS[1]):
            re, im = two_l(H, kossakowski_shared(N, eta), N)
            held = pairs_hold(re, im, idx, lambda p, q: Rational(-4 * centre(N, eta, p, q)), N)
            claim(f"N={N} eta {eta}: characteristic polynomials pair about -c(p, q) for every block and its image", len(held) == len(idx), f"({len(held)}/{len(idx)})")
            # the one centre F137 allows the whole spectrum is trace/dim = -N gamma; a block and its image pair
            # there exactly when d^2 + m^2 = N (none at N = 3, which is not a sum of two squares)
            held_s = pairs_hold(re, im, idx, lambda p, q: Rational(-4 * N * GAMMA), N)
            want = sorted(k for k in idx if (k[0] - k[1]) ** 2 + (k[0] + k[1] - N) ** 2 == N)
            claim(f"N={N} eta {eta}: about -N gamma only the blocks with d^2 + m^2 = N pair", sorted(held_s) == want, f"({len(held_s)}/{len(idx)}: {sorted(held_s)})")
    # the GHZ coherence |0..0><1..1| is the one-dimensional block (0, N); paired about F137's -N gamma its
    # partner would sit at -2 N gamma - lambda = 2 gamma eta N (N - 1), above zero for eta > 0
    for N in (2, 3, 4, 5):
        d = 2 ** N
        for eta in ETAS:
            re, im = two_l(xxz_chain(N, 5, 2), kossakowski_shared(N, eta), N)
            u = 0 * d + (d - 1)
            lam = Fraction(int(re[u, u]), 2)
            ok = lam == -2 * GAMMA * ((1 - eta) * N + eta * N * N) and im[u, u] == 0
            partner = -2 * N * GAMMA - lam
            claim(f"N={N} eta {eta}: GHZ coherence at {lam}; a partner about -N gamma would sit at {partner}, above zero (forbidden) exactly when eta > 0",
                  ok and (partner > 0) == (eta > 0) and partner == 2 * GAMMA * eta * N * (N - 1))
        # any covariance pattern: real part -2 sigma - 4 S, partner about -sigma at +4 S, S = sum_{k<l} a_kl
        a = random_psd(N, rng)
        S = sum(int(a[k, l]) for k in range(N) for l in range(k + 1, N))
        re, im = two_l(xxz_chain(N, 5, 2), a, N)
        u = d - 1
        lam = Fraction(int(re[u, u]), 2)
        claim(f"N={N} random a (summed covariance {S}): GHZ coherence's real part -2 sigma - 4S = {lam}, partner about -sigma at +4S = {4 * S}",
              lam == -2 * int(np.trace(a)) - 4 * S)
    for N in (2, 3):
        d = 2 ** N
        for eta in (Fraction(0), ETAS[0], ETAS[1], ETAS[2]):
            re, im = two_l(xxz_chain(N, 5, 2), kossakowski_shared(N, eta), N)
            check(f"N={N} eta {eta}: trace(L)/dim in units of gamma (F137's candidate centre)", Fraction(-N), Fraction(int(np.trace(re)), 2 * d * d * GAMMA))
            p_all = charpoly(re, im)
            pal = p_all == p_all.compose(Poly(Rational(-4 * N * GAMMA) - Y_SYM, Y_SYM, domain=QQ_I))
            claim(f"N={N} eta {eta}: the whole characteristic polynomial pairs about -N gamma", pal == (eta == 0), f"({'pairs' if pal else 'does not pair'})")
    # the spectral pairing does not need a real H: ρ -> ρ^dag X^N carries L to -L^dag - {K, .} for any H
    # commuting with X^N, and L preserves hermiticity; here a complex X^N-even H, for which the Pi_Z
    # similarity fails, still pairs every block with its image about sigma(p, q). On an open chain every
    # popcount sector's hopping graph is a tree and a diagonal unitary removes the phases, so the H sits on
    # the triangle, where the single-excitation hopping carries a loop product that is not real: a flux
    # no diagonal unitary removes (each popcount-1 state is one configuration, so any diagonal unitary
    # multiplies the three amplitudes by phases whose product around the loop is 1)
    N = 3
    tri = [(0, 1), (1, 2), (2, 0)]
    third = {(0, 1): 2, (1, 2): 0, (2, 0): 1}
    dm = sum(3 * op(Z, third[(k, l)], N) @ (op(X, k, N) @ op(YR, l, N) - op(YR, k, N) @ op(X, l, N)) for k, l in tri)  # = -i * (3 Z_m (X_k Y_l - Y_k X_l))
    Hc = xxz_chain(N, 5, 2, bonds=tri) + GaussH(np.zeros((8, 8), dtype=np.int64), dm)
    XN3 = np.kron(np.kron(X, X), X)
    even = np.array_equal(XN3 @ Hc.r @ XN3, Hc.r) and np.array_equal(XN3 @ Hc.i @ XN3, Hc.i)
    one = [1 << (N - 1 - k) for k in range(N)]
    loop = (1, 0)
    for k, l in tri:
        a_r, a_i = int(Hc.r[one[l], one[k]]), int(Hc.i[one[l], one[k]])
        loop = (loop[0] * a_r - loop[1] * a_i, loop[0] * a_i + loop[1] * a_r)
    claim(f"N=3 the complex H on the triangle, XXZ + 3 sum Z_m (X_k Y_l - Y_k X_l) (m the third site): commutes with X^N, "
          f"and its single-excitation loop product {loop[0]} + {loop[1]} i is not real", even and loop[1] != 0)
    for eta in (ETAS[0], ETAS[1]):
        re, im = two_l(Hc, kossakowski_shared(N, eta), N)
        cre, cim = conj_by_pi(re, pi_perm(N)), conj_by_pi(im, pi_perm(N))
        idx = blocks(N)
        simil = all(np.array_equal(cre[np.ix_(idx[(q, N - p)], idx[(q, N - p)])],
                                   -re[np.ix_(idx[(q, N - p)], idx[(q, N - p)])] - int(4 * centre(N, eta, p, q)) * np.eye(len(idx[(q, N - p)]), dtype=np.int64))
                    and np.array_equal(cim[np.ix_(idx[(q, N - p)], idx[(q, N - p)])], -im[np.ix_(idx[(q, N - p)], idx[(q, N - p)])])
                    for (p, q) in idx)
        held = pairs_hold(re, im, idx, lambda p, q: Rational(-4 * centre(N, eta, p, q)), N)
        claim(f"N=3 eta {eta}, that complex X^N-even H on the triangle: the Pi_Z similarity fails, "
              f"yet every block pairs with its image about -sigma(p, q)", (not simil) and len(held) == 16, f"({len(held)}/16)")
    check("N=4 eta 3/10: sigma(2, 2), the pair centre being its negative", Fraction(28), centre(4, ETAS[0], 2, 2))
    check("N=4 eta 3/10: sigma(0, 4), the pair centre being its negative", Fraction(76), centre(4, ETAS[0], 0, 4))

    print("\nF. A uniform field h sum Z on top: the pair sum gains 2 i h (2p - N)")
    hf = 3
    for N in (2, 3, 4, 5):
        d = 2 ** N
        perm = pi_perm(N)
        idx = blocks(N)
        field = GaussH(hf * sum(op(Z, k, N) for k in range(N)))
        for eta in (Fraction(0), ETAS[0], ETAS[2]):
            re, im = two_l(xxz_chain(N, 5, 2) + field, kossakowski_shared(N, eta), N)
            cre, cim = conj_by_pi(re, perm), conj_by_pi(im, perm)
            ok = True
            for (p, q), ix in idx.items():
                jx = idx[(q, N - p)]
                eye = np.eye(len(jx), dtype=np.int64)
                ok &= np.array_equal(cre[np.ix_(jx, jx)], -re[np.ix_(jx, jx)] - int(4 * centre(N, eta, p, q)) * eye)
                ok &= np.array_equal(cim[np.ix_(jx, jx)], -im[np.ix_(jx, jx)] + 4 * hf * (2 * p - N) * eye)
            claim(f"N={N} eta {eta} h {hf}: Pi_Z maps block (p, q) onto -L - 2c + 2ih(2p - N) on (q, N - p)", ok)
    for N in (2, 3, 4, 5):                                   # local noise at site-dependent rates: centre sigma
        d = 2 ** N
        perm = pi_perm(N)
        idx = blocks(N)
        aloc = np.diag([3, 5, 2, 7, 4][:N]).astype(np.int64)
        sigma = int(np.trace(aloc))
        re, im = two_l(xxz_chain(N, 5, 2) + GaussH(hf * sum(op(Z, k, N) for k in range(N))), aloc, N)
        cre, cim = conj_by_pi(re, perm), conj_by_pi(im, perm)
        ok = True
        for (p, q), ix in idx.items():
            jx = idx[(q, N - p)]
            eye = np.eye(len(jx), dtype=np.int64)
            ok &= np.array_equal(cre[np.ix_(jx, jx)], -re[np.ix_(jx, jx)] - 4 * sigma * eye)
            ok &= np.array_equal(cim[np.ix_(jx, jx)], -im[np.ix_(jx, jx)] + 4 * hf * (2 * p - N) * eye)
        claim(f"N={N} local rates {[int(x) for x in np.diag(aloc)]}, h {hf}: block pairs sum to -2 sigma + 2ih(2p - N)", ok)
    N = 3
    for eta in (Fraction(0), ETAS[0]):
        re, im = two_l(xxz_chain(N, 5, 2) + GaussH(hf * sum(op(Z, k, N) for k in range(N))), kossakowski_shared(N, eta), N)
        held = pairs_hold(re, im, blocks(N), lambda p, q: Rational(-4 * centre(N, eta, p, q)) + 4 * hf * (2 * p - N) * SYM_I, N)
        claim(f"N=3 eta {eta} h 3: characteristic polynomials pair about -c(p, q) + i h (2p - N)", len(held) == 16, f"({len(held)}/16)")
        held_r = pairs_hold(re, im, blocks(N), lambda p, q: Rational(-4 * centre(N, eta, p, q)), N)
        claim(f"N=3 eta {eta} h 3: without the imaginary part of the centre none pairs", len(held_r) == 0, f"({len(held_r)}/16)")
    # a real non-uniform field sum_k h_k Z_k: not a constant on a block, so the block identity fails; the
    # spectral pairing (at the only shift the traces allow) survives on some block pairs and not on others
    for N, h, n_held in ((3, [3, 1, 4], 8), (4, [3, 1, 4, 5], 12), (3, [1, 2, 3], 12), (4, [3, 1, 4, 2], 16)):
        kind = "anti-palindromic" if len({h[k] + h[N - 1 - k] for k in range(N)}) == 1 else "generic"
        re, im = two_l(xxz_chain(N, 5, 2) + GaussH(sum(hk * op(Z, k, N) for k, hk in enumerate(h))), kossakowski_shared(N, ETAS[0]), N)
        perm = pi_perm(N)
        idx = blocks(N)
        cre, cim = conj_by_pi(re, perm), conj_by_pi(im, perm)
        scalar = all(len(set(np.diag(cim[np.ix_(idx[(q, N - p)], idx[(q, N - p)])] + im[np.ix_(idx[(q, N - p)], idx[(q, N - p)])]))) == 1
                     for (p, q) in idx)
        claim(f"control N={N} field {h}: the field is not one constant per block pair, the identity fails", not scalar)
        held = pairs_hold(re, im, idx, trace_shift(re, im, idx, N), N)
        ends = all((p, q) in held for p in (0, N) for q in range(N + 1))
        claim(f"N={N} {kind} field {h}: {n_held} of {len(idx)} blocks still pair with their image, every block with an empty or full ket among them",
              len(held) == n_held and ends, f"({len(held)}/{len(idx)})")
        # the operator reason: Pi_Z leaves the field term -2i h.z(i) on the image of |i><j|, one number where the
        # ket is a single configuration; (site reversal) o Pi_Z (the chain reflection, a symmetry of the uniform chain and of
        # the shared field) leaves -i[(h + h_R).z(i) + (h_R - h).z(j)], one number where the bra is single, if
        # h + h_R is constant, that is for an anti-palindromic field
        d = 2 ** N
        rev = lambda x: int(format(x, f"0{N}b")[::-1], 2)
        rperm = np.empty(d * d, dtype=np.int64)
        for i, j in itertools.product(range(d), range(d)):
            rperm[i * d + j] = rev(i) * d + rev(j)
        comp = rperm[pi_perm(N)]                                       # |i><j| -> R|j><i_bar|R
        cim_r = conj_by_pi(im, comp)
        cre_r = conj_by_pi(re, comp)

        def scalar_on(crex, cimx, p, q):
            """The full residual conj(L) + L on the image block, real and imaginary parts, is a multiple of 1."""
            jx = idx[(q, N - p)]
            ok = True
            for M in (crex[np.ix_(jx, jx)] + re[np.ix_(jx, jx)], cimx[np.ix_(jx, jx)] + im[np.ix_(jx, jx)]):
                ok &= not (M - np.diag(np.diag(M))).any() and len(set(np.diag(M))) == 1
            return ok
        ket_single = all(scalar_on(cre, cim, p, q) for p in (0, N) for q in range(N + 1))
        bra_single = all(scalar_on(cre_r, cim_r, p, q) for q in (0, N) for p in range(N + 1))
        claim(f"N={N} {kind} field {h}: the field term is one number on every empty- or full-ket block under Pi_Z"
              + (", and on every empty- or full-bra block under (site reversal) o Pi_Z" if kind == "anti-palindromic" else ""),
              ket_single and (bra_single if kind == "anti-palindromic" else True))
        if kind == "generic":
            # (N/2, q) maps onto (q, N/2), its own Hermitian-conjugate block, so it pairs with its image exactly
            # when its spectrum is symmetric about the vertical line Re = -sigma(p, q); with q = 0 or N the bra is
            # one configuration, the dissipator one number on the block at uniform rates, the spectrum lies on
            # that line, and it does
            expect = {(p, q) for (p, q) in idx if p in (0, N)} | ({(N // 2, 0), (N // 2, N)} if N % 2 == 0 else set())
            claim(f"N={N} generic field {h}: the survivors are the empty- or full-ket blocks"
                  + (f" and (N/2, 0), (N/2, N), whose image is their own conjugate block" if N % 2 == 0 else ""),
                  set(held) == expect, f"({len(expect)} blocks)")
            if N % 2 == 0:
                g = [7, 12, 9, 15][:N]
                anu = (np.diag(g) + 3 * (np.ones((N, N), dtype=np.int64) - np.eye(N, dtype=np.int64))).astype(np.int64)
                re2, im2 = two_l(xxz_chain(N, 5, 2) + GaussH(sum(hk * op(Z, k, N) for k, hk in enumerate(h))), anu, N)
                held2 = pairs_hold(re2, im2, idx, trace_shift(re2, im2, idx, N), N)
                claim(f"N={N} generic field {h}, site-dependent rates {g} with covariance 3: (N/2, 0), (N/2, N) drop out, "
                      "the empty- or full-ket blocks stay", set(held2) == {(p, q) for (p, q) in idx if p in (0, N)},
                      f"({len(held2)}/{len(idx)})")
        if kind == "anti-palindromic":
            union = {(p, q) for (p, q) in idx if p in (0, N) or q in (0, N)}
            claim(f"N={N} anti-palindromic field {h}: the surviving blocks are exactly those with an empty or full ket or bra",
                  set(held) == union, f"({len(union)} blocks)")

    print("\nG. The two ladders on the XY chain (Jordan-Wigner c_l)")
    for N in (3, 4, 5):
        d = 2 ** N
        Id = np.eye(d, dtype=np.int64)
        sm = np.array([[0, 1], [0, 0]], dtype=np.int64)            # |0><1|, |1> occupied
        cs = []
        for l in range(N):
            c = np.ones((1, 1), dtype=np.int64)
            for s in range(N):
                c = np.kron(c, Z if s < l else (sm if s == l else I2))
            cs.append(c)
        W = sum(np.kron(c.T, c.T) for c in cs)                      # c^dag rho c: kron(c^dag, c^T), c real
        Sp = sum((-1) ** l * np.kron(c.T, c) for l, c in enumerate(cs))   # c^dag rho c^dag: kron(c^dag, (c^dag)^T)
        H = xxz_chain(N, 1, 0)
        claim(f"N={N}: Pi_Z W Pi_Z^-1 = S+ exactly", np.array_equal(conj_by_pi(W, pi_perm(N)), Sp))
        bra_flip = np.array([i * d + (d - 1 - j) for i in range(d) for j in range(d)], dtype=np.int64)   # R: rho -> rho X^N
        transpose = np.array([j * d + i for i in range(d) for j in range(d)], dtype=np.int64)            # D: rho -> rho^T
        claim(f"N={N}: the transpose leaves W fixed and the bra flip R(rho) = rho X^N alone carries W onto S+",
              np.array_equal(conj_by_pi(W, transpose), W) and np.array_equal(conj_by_pi(W, bra_flip), Sp))
        while True:
            Bi = rng.integers(-4, 5, size=(N, N))
            a_ind = (Bi + Bi.T).astype(np.int64)
            if (a_ind != 0).all() and np.linalg.eigvalsh(a_ind.astype(float)).min() < 0:
                break
        rei, imi = two_l(H, a_ind, N)
        claim(f"N={N} indefinite symmetric a: [L, W] = 0 still", comm_zero(rei, W) and comm_zero(imi, W))
        a_sd = (np.diag([7, 12, 9, 15, 10][:N]) + 3 * (np.ones((N, N), dtype=np.int64) - np.eye(N, dtype=np.int64))).astype(np.int64)
        for name, a in (("random PSD a", random_psd(N, rng)), ("shared field eta 3/10", kossakowski_shared(N, ETAS[0])),
                        ("shared field at site-dependent rates (covariance 3)", a_sd)):
            re, im = two_l(H, a, N)
            ok_w = comm_zero(re, W) and comm_zero(im, W)
            claim(f"N={N} {name}: [L, W] = 0", ok_w)
            K = k_operator(a, N)
            antiK = np.kron(K, Id) + np.kron(Id, K.T)
            claim(f"N={N} {name}: [S+, L] = -[S+, {{K, .}}] (the price identity carried through Pi_Z)",
                  np.array_equal(comm(Sp, re), -2 * comm(Sp, antiK)) and comm_zero(Sp, im))
            aloc = np.diag(np.diag(a))
            rl, il = two_l(H, aloc, N)
            ok_loc = comm_zero(rl, Sp) and comm_zero(il, Sp)
            claim(f"N={N} {name}, diagonal part only: [L, S+] = 0", ok_loc)
            roff, _ = two_l(GaussH(np.zeros((d, d), dtype=np.int64)), a - aloc, N)
            br = comm(re, Sp)
            claim(f"N={N} {name}: [L, S+] != 0 and equals [D_offdiagonal, S+]",
                  br.any() and np.array_equal(br, comm(roff, Sp)) and comm_zero(im, Sp))
            # is the failure only a constant shift on each source block (S+ still an intertwiner of blocks)?
            dd = np.array([popcount(i) - popcount(j) for i in range(d) for j in range(d)], dtype=np.int64)
            Spd = Sp.astype(np.int64)
            if name.startswith("shared"):
                kap = int(a[0, 1])
                claim(f"N={N} {name}: [L, S+] = -8 kappa (d + 1) S+ on every source block (p, q), d = p - q: S+ still "
                      "carries eigenvectors to eigenvectors, the eigenvalue shifted",
                      np.array_equal(br, -16 * kap * Spd * (dd + 1)[np.newaxis, :]))
            else:
                shift_only = True
                for (p, q), ix in blocks(N).items():
                    sub_br, sub_sp = br[:, ix], Spd[:, ix]
                    if (sub_br[sub_sp == 0] != 0).any():
                        shift_only = False
                        break
                    ratios = {Fraction(int(x), int(y)) for x, y in zip(sub_br[sub_sp != 0], sub_sp[sub_sp != 0])}
                    if len(ratios) > 1:
                        shift_only = False
                        break
                claim(f"N={N} {name}: [L, S+] is not S+ times one constant per source block (no block shift)", not shift_only)
        a1 = np.diag([GAMMA] * N).astype(np.int64)
        a1[0, N - 1] = a1[N - 1, 0] = 1                             # one small covariance, between the two ends
        re, im = two_l(H, a1, N)
        claim(f"N={N} a single covariance a_0,{N - 1} = 1 beside local rates {GAMMA}: W still commutes, S+ does not",
              comm_zero(re, W) and comm_zero(im, W) and not comm_zero(re, Sp))
        re, im = two_l(xxz_chain(N, 1, 1), np.diag([GAMMA] * N).astype(np.int64), N)
        claim(f"control N={N}: the Heisenberg chain breaks W through the H door (local noise)",
              not (comm_zero(re, W) and comm_zero(im, W)))

    print("\nH. The dark cells: for a positive semidefinite a, |i><j| has rate 0 exactly when u = z(i) - z(j) lies in the kernel of a")
    for N in (2, 3, 4, 5):
        d = 2 ** N
        anti = (N * np.eye(N, dtype=np.int64) - np.ones((N, N), dtype=np.int64)).astype(np.int64)
        # the anti-correlated end eta = -1/(N - 1) at gamma = N - 1: a_kk = N - 1, a_kl = -1, every row sum 0
        ak, v = random_psd_with_kernel(N, rng)
        while True:                                         # a weighted graph Laplacian: a >= 0, rows sum to 0
            lap = np.zeros((N, N), dtype=np.int64)
            ws = []
            for k in range(N):
                for l in range(k + 1, N):
                    w = int(rng.integers(1, 4))
                    ws.append(w)
                    lap[k, k] += w
                    lap[l, l] += w
                    lap[k, l] -= w
                    lap[l, k] -= w
            if N == 2 or len(set(ws)) > 1:                  # not uniform from N = 3 on (N = 2 has one weight)
                break
        cases = (("random positive definite a", random_psd(N, rng), 2 ** N),
                 ("shared field, eta 3/10", kossakowski_shared(N, ETAS[0]), 2 ** N),
                 ("shared field, eta 1", kossakowski_shared(N, Fraction(1)), comb(2 * N, N)),
                 ("anti-correlated end, eta -1/(N - 1)", anti, 2 ** N + 2),
                 (f"rank-deficient random a, ker a = span{tuple(int(x) for x in v)}", ak, None),
                 ("random weighted graph Laplacian", lap, None))
        ghz_ok = True
        for name, a, closed in cases:
            dark = dark_cells(a, N)
            kernel = {(i, j) for i in range(d) for j in range(d)
                      if not (a @ np.array([zval(i, k, N) - zval(j, k, N) for k in range(N)], dtype=np.int64)).any()}
            if closed is None:
                claim(f"N={N} {name}: the cells of rate zero are exactly those with a u = 0, coherences among them",
                      dark == kernel and len(dark) > d)
            else:
                claim(f"N={N} {name}: the cells of rate zero are exactly those with a u = 0", dark == kernel)
                check(f"N={N} {name}: dark cells", closed, len(dark))
            ghz_ok &= ((0, d - 1) in dark) == (not (a @ np.ones(N, dtype=np.int64)).any())
        claim(f"N={N}: on all six, the GHZ coherence is dark exactly when every row of a sums to zero (uniform at the anti-correlated end, not uniform on the Laplacian from N = 3)", ghz_ok)
        coh = {cell for cell in dark_cells(anti, N) if cell[0] != cell[1]}
        claim(f"N={N} anti-correlated end: the dark coherences are exactly |0...0><1...1| and its transpose",
              coh == {(0, d - 1), (d - 1, 0)})
        claim(f"N={N} on all six: the dark cells are exactly the cells inside the classes of a z",
              all(dark_cells(a, N) == {(i, j) for c in dark_classes(a, N) for i in c for j in c} for _, a, _ in cases))
        claim(f"N={N} eta 1: the dark classes are the popcount sectors",
              sorted(dark_classes(kossakowski_shared(N, Fraction(1)), N))
              == sorted([i for i in range(d) if popcount(i) == p] for p in range(N + 1)))
        claim(f"N={N} anti-correlated end: the one class beyond single configurations is the GHZ pair",
              [c for c in dark_classes(anti, N) if len(c) > 1] == [[0, d - 1]])

    print("\nH, continued. Classes that multiply, and the ring at the edge of positive semidefiniteness")
    for N in (2, 4):
        d = 2 ** N
        apairs = np.zeros((N, N), dtype=np.int64)
        for k in range(0, N, 2):
            apairs[k, k] = apairs[k + 1, k + 1] = 1
            apairs[k, k + 1] = apairs[k + 1, k] = -1
        check(f"N={N} anti-correlated pairs: dark cells, (2^2 + 2)^(N/2)", 6 ** (N // 2), len(dark_cells(apairs, N)))
        alike = [i for i in range(d) if all(((i >> k) & 1) == ((i >> (k + 1)) & 1) for k in range(0, N, 2))]
        claim(f"N={N} anti-correlated pairs: the largest class takes 00 or 11 on every pair",
              max(dark_classes(apairs, N), key=len) == alike)
    for N in (3, 4, 5):
        d = 2 ** N
        adj = np.zeros((N, N), dtype=np.int64)
        for k in range(N):
            adj[k, (k + 1) % N] = adj[(k + 1) % N, k] = 1
        lap, signless = 2 * np.eye(N, dtype=np.int64) - adj, 2 * np.eye(N, dtype=np.int64) + adj
        coh = {c for c in dark_cells(lap, N) if c[0] != c[1]}
        claim(f"N={N} ring, kappa = -1 at gamma = 2 (D - A): the dark coherences are the GHZ pair's", coh == {(0, d - 1), (d - 1, 0)})
        coh = {c for c in dark_cells(signless, N) if c[0] != c[1]}
        neel = sum(1 << k for k in range(0, N, 2))
        if N % 2 == 0:
            claim(f"N={N} ring, kappa = +1 at gamma = 2 (D + A): the dark coherences are the Neel pair's",
                  coh == {(neel, neel ^ (d - 1)), (neel ^ (d - 1), neel)})
        else:
            claim(f"N={N} ring, kappa = +1 at gamma = 2 (D + A): no dark coherence, no alternating pattern closes", coh == set())
        if N == 4:
            m = 0b0101                                       # X on every other site
            claim("N=4: the flip on every other site carries the dark cells of D + A onto those of D - A",
                  {(i ^ m, j ^ m) for i, j in dark_cells(signless, N)} == dark_cells(lap, N))
            H = xxz_chain(N, 1, 0).r                            # the open XY chain
            unit = lambda i, j: csr_matrix(([1], ([i], [j])), shape=(d, d)).toarray()
            claim("N=4 open XY chain: it conserves each popcount",
                  all(comm_zero(H, np.diag([1 if popcount(i) == p else 0 for i in range(d)])) for p in range(N + 1)))
            sectors = dark_cells(kossakowski_shared(N, Fraction(1)), N)
            claim("N=4 open XY chain: the commutator of every dark cell at eta = 1 stays inside the popcount sectors",
                  all({(r, c) for r, c in zip(*np.nonzero(comm(H, unit(i, j))))} <= sectors for i, j in sectors))
            out = {(r, c) for r, c in zip(*np.nonzero(comm(H, unit(neel, neel ^ (d - 1)))))}
            claim("N=4 open XY chain: it commutes with the GHZ cell |0000><1111| and carries the Neel cell out of the dark cells",
                  comm_zero(H, unit(0, d - 1)) and bool(out) and not out <= dark_cells(signless, N))

    print("\nI. Which popcounts see one value of K: x = 0 at p = 1, N - 1; Omega = 0 at p = N/2; x = Omega = 0 elsewhere")
    lrng = np.random.default_rng(20261011)

    def sym(N, entries):
        c = np.zeros((N, N), dtype=np.int64)
        t = 0
        for k in range(N):
            for l in range(k + 1, N):
                c[k, l] = c[l, k] = entries[t]
                t += 1
        return c

    for N, lo, hi in ((3, -2, 2), (4, -2, 2), (5, -1, 1)):
        diag = [3 + 2 * k for k in range(N)]                           # unequal rates: the diagonal enters only through sigma
        shapes, ok, d = {}, True, 2 ** N
        for entries in itertools.product(range(lo, hi + 1), repeat=N * (N - 1) // 2):
            c = sym(N, entries)
            got = qualifying_popcounts(c, diag, N)
            ok &= got == predicted_popcounts(c, N)
            shapes[tuple(got)] = shapes.get(tuple(got), 0) + 1
        claim(f"N={N} every pattern with entries in {lo}..{hi} at unequal rates: the qualifying popcounts are the law's", ok)
        expect = {3: {(0, 3), (0, 1, 2, 3)}, 4: {(0, 4), (0, 1, 3, 4), (0, 2, 4), (0, 1, 2, 3, 4)},
                  5: {(0, 5), (0, 1, 4, 5), (0, 1, 2, 3, 4, 5)}}[N]
        claim(f"N={N}: the shapes seen are exactly {sorted(len(s) ** 2 for s in expect)} blocks", set(shapes) == expect)

    def additive(N):
        y = lrng.integers(-3, 4, N)
        return np.array([[0 if k == l else y[k] + y[l] for l in range(N)] for k in range(N)], dtype=np.int64), y

    def circulant(N):
        c = np.zeros((N, N), dtype=np.int64)
        for s_ in range(1, N // 2 + 1):
            w = int(lrng.integers(-3, 4))
            for k in range(N):
                l = (k + s_) % N
                if 2 * s_ == N and k >= s_:
                    continue
                c[k, l] += w
                c[l, k] += w
        return c

    def random_pattern(N):
        e = lrng.integers(-3, 4, N * (N - 1) // 2)
        return sym(N, e)

    for N in (6, 7, 8):
        diag = [1 + k for k in range(N)]
        ok, mean_ok, add_ok, d = True, True, True, 2 ** N
        for _ in range(60):
            for c in (additive(N)[0], circulant(N), additive(N)[0] + circulant(N), random_pattern(N)):
                got = qualifying_popcounts(c, diag, N)
                ok &= got == predicted_popcounts(c, N)
                T, kbar, _, _ = three_parts(c, N)
                kv = k_on_configs(c, diag, N)
                sigma = sum(diag)
                for p_ in got:
                    kp = int(next(kv[i] for i in range(d) if popcount(i) == p_))
                    mean_ok &= comb(N, 2) * kp == comb(N, 2) * sigma - N * T + T * (N - 2 * p_) ** 2
            c, y = additive(N)
            kv = k_on_configs(c, diag, N)
            add_ok &= all(int(kv[i]) == sum(diag) - 2 * int(y.sum())
                          + 2 * (N - 2 * popcount(i)) * sum(int(y[k]) * zval(i, k, N) for k in range(N)) for i in range(d))
        claim(f"N={N} additive, circulant (equal row sums), their sums and random patterns: the qualifying popcounts are the law's", ok)
        claim(f"N={N}: on every qualifying popcount K_p = sigma - N kbar + kbar (N - 2p)^2, the shared field at the mean covariance", mean_ok)
        claim(f"N={N} additive C: K = sigma - 2 sum y + 2 (N - 2p) sum_k y_k z_k on every configuration", add_ok)

    for N in (4, 5, 6):
        diag = [2 + 3 * k for k in range(N)]
        ident_ok = True
        for _ in range(40):
            c = random_pattern(N)
            T, kbar, x, Rm = three_parts(c, N)
            ident_ok &= sum(x) == 0 and all(sum(Rm[k]) == 0 for k in range(N))
            ident_ok &= all(Fraction(int(c[k][l])) == kbar + x[k] + x[l] + Rm[k][l]
                            for k in range(N) for l in range(N) if k != l)
            kv = k_on_configs(c, diag, N)
            for i in range(2 ** N):
                exc = [k for k in range(N) if zval(i, k, N) == -1]
                cut = sum(int(c[k][l]) for k in exc for l in range(N) if l not in exc)
                pp = len(exc)
                ident_ok &= int(kv[i]) == sum(diag) + 2 * int(T) - 4 * cut
                ident_ok &= cut == kbar * pp * (N - pp) + (N - 2 * pp) * sum(x[k] for k in exc) \
                    - 2 * sum(Rm[k][l] for k in exc for l in exc if k < l)
        claim(f"N={N} random patterns at unequal rates: K = sigma + 2T - 4 cut(A), the parts of the covariances add back "
              f"with sum x = 0 and zero row sums of Omega, and cut(A) = kbar p(N - p) + (N - 2p) x(A) - 2 Omega(A), set by set", ident_ok)

    s4 = (1, 1, 1, -1)
    sign = np.array([[0 if k == l else 3 * s4[k] * s4[l] for l in range(4)] for k in range(4)], dtype=np.int64)
    _, kbar4, x4, Rm = three_parts(sign, 4)
    y4 = [x4[k] + kbar4 / 2 for k in range(4)]
    claim("the page's sign pattern s = (1, 1, 1, -1) is additive, y proportional to (1, 1, 1, -3)",
          all(v == 0 for row in Rm for v in row) and predicted_popcounts(sign, 4) == [0, 2, 4]
          and y4[0] == y4[1] == y4[2] and y4[3] == -3 * y4[0] != 0
          and qualifying_popcounts(sign, [10] * 4, 4) == [0, 2, 4])

    for N in (4, 5):
        d = 2 ** N
        H = xxz_chain(N, 1, 0).r
        diag = [GAMMA] * N
        agree = True
        for c in (random_pattern(N), circulant(N), additive(N)[0], sign if N == 4 else random_pattern(N)):
            kv = k_on_configs(c, diag, N)
            qual = qualifying_popcounts(c, diag, N)
            for p_ in range(N + 1):
                idx = [i for i in range(d) if popcount(i) == p_]
                Hs = H[np.ix_(idx, idx)]
                Ks = np.diag([int(kv[i]) for i in idx])
                agree &= comm_zero(Hs, Ks) == (p_ in qual)
        claim(f"N={N} open XY chain: K is one number on popcount p exactly when it commutes with the chain there", agree)

    for N in (4, 5, 6, 7):
        bonds = {"ring": [(k, (k + 1) % N) for k in range(N)], "star": [(0, l) for l in range(1, N)],
                 "chain": [(k, k + 1) for k in range(N - 1)]}
        want = {"ring": [0, 1, N - 1, N], "star": [0, N // 2, N] if N % 2 == 0 else [0, N], "chain": [0, N]}
        for name, bl in bonds.items():
            w = np.zeros((N, N), dtype=np.int64)
            for a_, b_ in bl:
                w[a_, b_] = w[b_, a_] = 1
            got = qualifying_popcounts(w, [0] * N, N)
            claim(f"N={N} ZZ on the bonds of the {name}: one number on popcounts {got}, the law's {want[name]}",
                  got == want[name] == predicted_popcounts(w, N))

    print("\nJ. the third moment about the one candidate centre: Tr((L + sigma)^3) = -4^N [4 tr(C(2a + Lambda)C) - 3 tr(C Theta) - 6 xi^T a xi]")
    jrng = np.random.default_rng(20261012)
    for N in (2, 3, 4):
        agree = True
        for _ in range(3):
            a = random_symmetric(N, jrng)
            H = random_hermitian(N, jrng)
            re, im = two_l(H, a, N)
            got = power_trace(re, im, 2 * int(np.trace(a)), 3)
            agree &= got == (third_moment_closed(H, a, N), 0)
        claim(f"N={N} three random real symmetric a of either sign, three random complex traceless H: the closed form is the trace, the imaginary part zero", agree)
    for N in (2, 3, 4):
        Zs = [op(Z, k, N) for k in range(N)]
        XN = np.ones((1, 1), dtype=np.int64)
        for _ in range(N):
            XN = np.kron(XN, X)
        agree, no_c = True, True
        for _ in range(2):
            H = random_xn_even(N, jrng)
            t = int(np.trace(H.r))                       # a multiple of 4; the shift keeps H commuting with X^N
            H.r[0, 0] -= t // 2
            H.r[-1, -1] -= t // 2
            no_c &= (int(np.trace(H.r)) == 0 and comm_zero(H.r, XN)
                     and all(int(np.trace(H.r @ Zs[k])) == 0 for k in range(N)))
            a = random_symmetric(N, jrng)
            re, im = two_l(H, a, N)
            agree &= power_trace(re, im, 2 * int(np.trace(a)), 3) == (third_moment_closed(H, a, N), 0)
        claim(f"N={N} two traceless H commuting with X^N: Tr(H Z_k) = 0 at every site (xi = 0), and the closed form is the trace", agree and no_c)
    for N in (2, 3, 4):
        diag = [int(v) for v in jrng.integers(1, 9, size=N)]
        h = [int(v) for v in jrng.integers(-3, 4, size=N)]
        Zs = [op(Z, k, N) for k in range(N)]
        H = xxz_chain(N, 2, 1) + GaussH(sum(h[k] * Zs[k] for k in range(N)))
        a = np.diag(diag).astype(np.int64)
        re, im = two_l(H, a, N)
        got = power_trace(re, im, 2 * sum(diag), 3)
        want = 8 * 6 * 4 ** N * sum(diag[k] * h[k] ** 2 for k in range(N))
        claim(f"N={N} local noise at rates {diag}, the XXZ chain plus fields {h}: the trace is 8 * 6 * 4^N sum gamma_k h_k^2, F117's third-moment face carried to site-dependent rates",
              got == (want, 0) and third_moment_closed(H, a, N) == want)
    for N in (2, 3, 4, 5):
        agree = True
        for _ in range(4):
            a = random_symmetric(N, jrng)
            re, _ = two_l(GaussH(np.zeros((2 ** N, 2 ** N), dtype=np.int64)), a, N)
            dg = [int(v) + 2 * int(np.trace(a)) for v in np.diag(re)]
            agree &= sum(v ** 3 for v in dg) == -8 * 4 ** N * noise_part(a)
        claim(f"N={N} H = 0, four random real symmetric a of either sign: the cube sum of the cell entries is -8 * 4^N * 4 tr(C(2a + Lambda)C)", agree)
    agree = True
    for N in range(2, 9):
        for _ in range(4):
            a = random_symmetric(N, jrng)
            A = [[int(a[k][l]) for l in range(N)] for k in range(N)]
            edges = sum(A[k][l] ** 2 * (A[k][k] + A[l][l]) for k in range(N) for l in range(k + 1, N))
            tri = sum(A[k][l] * A[l][m] * A[k][m] for k in range(N) for l in range(k + 1, N) for m in range(l + 1, N))
            C = np.array([[A[k][l] if k != l else 0 for l in range(N)] for k in range(N)], dtype=object)
            Dl = np.diag([A[k][k] for k in range(N)]).astype(object)
            agree &= noise_part(a) == 12 * edges + 48 * tri == 12 * int(np.trace(Dl.dot(C).dot(C))) + 8 * int(np.trace(C.dot(C).dot(C)))
    claim("N=2..8 the noise term is 12 sum_{k<l} a_kl^2 (a_kk + a_ll) + 48 sum_triangles a_kl a_lm a_km = 12 tr(Lambda C^2) + 8 tr(C^3)", agree)
    agree, low = True, None
    for N in range(2, 9):
        for _ in range(4):
            a = random_psd(N, jrng)
            A = np.array([[int(a[k][l]) for l in range(N)] for k in range(N)], dtype=object)
            C = A.copy()
            for k in range(N):
                C[k, k] = 0
            Dl = np.diag([A[k][k] for k in range(N)]).astype(object)
            caC, cdC = int(np.trace(C.dot(A).dot(C))), int(np.trace(C.dot(Dl).dot(C)))
            agree &= noise_part(a) == 8 * caC + 4 * cdC and caC >= 0 and cdC > 0
            low = noise_part(a) if low is None else min(low, noise_part(a))
    claim("N=2..8 four positive definite a each: the noise term is 8 tr(CaC) + 4 tr(C Lambda C), the first >= 0 and the second > 0", agree, f"(smallest {low})")
    a_ind = np.array([[0, -1, -1], [-1, 0, -1], [-1, -1, 0]], dtype=np.int64)
    check("N=3 the indefinite a with zero diagonal and every covariance -1: noise term 8 tr(C^3)", -48, noise_part(a_ind))
    for N in (3, 4, 5, 6):
        Zs = [op(Z, k, N) for k in range(N)]
        for jz, name in ((0, "XY"), (1, "Heisenberg")):
            H = xxz_chain(N, 1, jz)
            Hsq = H.r @ H.r - H.i @ H.i
            h2 = {(k, l): Fraction(int(np.trace(Hsq @ Zs[k] @ Zs[l])), 2 ** N) for k in range(N) for l in range(k + 1, N)}
            want = {(k, l): (-2 if l == k + 1 else (2 * jz if l == k + 2 else 0)) for k in range(N) for l in range(k + 1, N)}
            claim(f"N={N} open {name} chain at J = 1: Tr(H^2 Z_k Z_l)/2^N is -2 on the bonds" + (", +2 on next-nearest pairs" if jz else "") + ", 0 elsewhere", h2 == want)
    agree = True
    for N in (2, 3, 4):
        for J in (1, 2):
            for eta in ETAS:
                kap = int(GAMMA * eta)
                g = GAMMA
                for jz in (0, 1):
                    re, im = two_l(xxz_chain(N, J, jz * J), kossakowski_shared(N, eta), N)
                    got = power_trace(re, im, 2 * N * g, 3)
                    if jz == 0:
                        m3 = -4 ** N * 4 * (N - 1) * kap * (2 * N * (N - 2) * kap ** 2 + 3 * N * g * kap + 3 * J ** 2)
                    else:
                        m3 = -4 ** N * 4 * kap * (3 * N * (N - 1) * g * kap + 2 * N * (N - 1) * (N - 2) * kap ** 2 + 3 * J ** 2)
                    agree &= got == (8 * m3, 0)
    claim("N=2..4, J in {1, 2}, every eta: the shared field's third moment on the open XY chain is -4^N 4 (N - 1) kappa [2N(N - 2) kappa^2 + 3N gamma kappa + 3J^2], "
          "on the open Heisenberg chain -4^N 4 kappa [3N(N - 1) gamma kappa + 2N(N - 1)(N - 2) kappa^2 + 3J^2]", agree)
    # where the GHZ argument says nothing: a negative covariance between two sites that share no bond of the XY chain
    a_nb = np.array([[10, 0, -3], [0, 10, 0], [-3, 0, 10]], dtype=np.int64)
    re, im = two_l(xxz_chain(3, 5, 0), a_nb, 3)
    check("N=3 open XY chain J = 5, covariance -3 between the end sites only: Tr((2L + 2 sigma)^3)", (-8 * 4 ** 3 * 2160, 0), power_trace(re, im, 60, 3))
    p_all = charpoly(re, im)
    claim("N=3 the same: the summed covariance is negative, and the whole characteristic polynomial does not pair about -sigma",
          sum(int(a_nb[k, l]) for k in range(3) for l in range(k + 1, 3)) < 0
          and p_all != p_all.compose(Poly(Rational(-120) - Y_SYM, Y_SYM, domain=QQ_I)))
    # where the third moment is blind: the shared field at the anti-correlated end, gamma = J = 10 at N = 3
    re, im = two_l(xxz_chain(3, 10, 0), kossakowski_shared(3, Fraction(-1, 2)), 3)
    check("N=3 open XY chain J = 10, shared field at the anti-correlated end (kappa = -5): Tr((2L + 2 sigma)^3)", (0, 0), power_trace(re, im, 60, 3))
    m5 = power_trace(re, im, 60, 5)
    check("N=3 the same: Tr((2L + 2 sigma)^5), 32 * 10^5 times the odd-moments script's -11520 at gamma = 1", (-36864000000, 0), m5)
    p_all = charpoly(re, im)
    claim("N=3 the same: the whole characteristic polynomial does not pair about -sigma (the fifth moment sees what the third does not)",
          p_all != p_all.compose(Poly(Rational(-120) - Y_SYM, Y_SYM, domain=QQ_I)) and m5[0] != 0)
    v1 = np.zeros(64, dtype=np.int64)
    for i in range(8):
        if popcount(i) == 1:
            v1[i * 8 + (7 - i)] = 1
    claim("N=3 the same: the far-end mode P_1 X^N is an eigenvector at -2 K_1 = -80 = -8 gamma, below -2 sigma = -60",
          np.array_equal(re @ v1, -160 * v1) and not np.any(im @ v1))

    print("\nK. the far-end modes: P_p X^N is an eigenvector at -2 K_p wherever K takes one value K_p on popcount p")

    def far_end_vec(N, p_):
        """vec(P_p X^N) in the row-major vec of two_l: the cells |i><i-bar| with popcount(i) = p_"""
        d = 2 ** N
        v = np.zeros(d * d, dtype=np.int64)
        for i in range(d):
            if popcount(i) == p_:
                v[i * d + (d - 1 - i)] = 1
        return v

    for N in (2, 3, 4, 5):
        d = 2 ** N
        for eta in (Fraction(-2, 10), Fraction(3, 10)):
            re, im = two_l(xxz_chain(N, 5, 2), kossakowski_shared(N, eta), N)
            sigma = N * GAMMA
            agree, below = True, []
            for p_ in range(N + 1):
                v = far_end_vec(N, p_)
                kp = sigma - N * GAMMA * eta + GAMMA * eta * (N - 2 * p_) ** 2
                agree &= (np.array_equal(re @ v, -4 * kp * v) and not np.any(im @ v)) if kp.denominator == 1 else False
                if kp > sigma:
                    below.append(p_)
            want = [p_ for p_ in range(N + 1) if ((N - 2 * p_) ** 2 > N if eta > 0 else (N - 2 * p_) ** 2 < N)]
            claim(f"N={N} eta {eta} XXZ chain: P_p X^N is an eigenvector at -2 K_p, K_p = sigma - N kappa + kappa (N - 2p)^2, for every p; "
                  f"K_p > sigma (a rate above 2 sigma, whose partner about -sigma would grow) at p = {below}",
                  agree and below == want and len(below) > 0 and (eta < 0) == (N // 2 in below))
    # the law: on a popcount where K is one number, K_p > sigma exactly when T ((N - 2p)^2 - N) > 0
    for N in (3, 4, 5):
        agree, seen = True, set()
        idx_pairs = [(k, l) for k in range(N) for l in range(k + 1, N)]
        for vals in itertools.product((-1, 0, 1), repeat=len(idx_pairs)):
            c = np.zeros((N, N), dtype=np.int64)
            for (k, l), v in zip(idx_pairs, vals):
                c[k, l] = c[l, k] = v
            kv = k_on_configs(c, [GAMMA] * N, N)
            T = sum(vals)
            for p_ in qualifying_popcounts(c, [GAMMA] * N, N):
                kp = int(kv[next(i for i in range(2 ** N) if popcount(i) == p_)])
                above = kp > N * GAMMA
                agree &= above == (T * ((N - 2 * p_) ** 2 - N) > 0)
                seen.add(above)
        claim(f"N={N} every pattern in -1..1: on every popcount where K is one number, K_p > sigma exactly when T ((N - 2p)^2 - N) > 0",
              agree and seen == {True, False})
    # an additive pattern at even N with a negative sum: K takes one value on popcount N/2, sigma - N kbar > sigma
    y = (-1, -1, 0, 0)
    a_add = np.array([[GAMMA if k == l else y[k] + y[l] for l in range(4)] for k in range(4)], dtype=np.int64)
    re, im = two_l(xxz_chain(4, 5, 2), a_add, 4)
    k2 = Fraction(int(k_on_configs(a_add - np.diag(np.diag(a_add)), [GAMMA] * 4, 4)[0b0011]))
    v = far_end_vec(4, 2)
    T = sum(int(a_add[k, l]) for k in range(4) for l in range(k + 1, 4))
    claim(f"N=4 additive pattern y = {y} (summed covariance {T}, so the GHZ argument is silent): K is one number {k2} on popcount 2, "
          f"above sigma = {4 * GAMMA}, and P_2 X^N is an eigenvector at -2 K_2",
          T < 0 and qualifying_popcounts(a_add - np.diag(np.diag(a_add)), [GAMMA] * 4, 4) == [0, 2, 4]
          and k2 == 4 * GAMMA - 4 * Fraction(T, 6) and k2 > 4 * GAMMA
          and np.array_equal(re @ v, -4 * int(k2) * v) and not np.any(im @ v))
    # the end-site covariance, where neither eigenvector argument reaches: K takes one value only on popcounts 0 and N,
    # and there it is sigma + 2T < sigma
    for N in (3, 4, 5):
        a_end = np.diag([GAMMA] * N).astype(np.int64)
        a_end[0, N - 1] = a_end[N - 1, 0] = -3
        c_end = a_end - np.diag(np.diag(a_end))
        kv = k_on_configs(c_end, [GAMMA] * N, N)
        claim(f"N={N} covariance -3 between the end sites only: K takes one value only on popcounts {qualifying_popcounts(c_end, [GAMMA] * N, N)}, "
              f"and there it is sigma - 6 < sigma, so no far-end mode and no GHZ coherence decides; the closed-form third moment "
              f"on the XY chain at J = 5 is {third_moment_closed(xxz_chain(N, 5, 0), a_end, N) // 8}, nonzero",
              qualifying_popcounts(c_end, [GAMMA] * N, N) == [0, N] and int(kv[0]) == N * GAMMA - 6 == int(kv[2 ** N - 1])
              and third_moment_closed(xxz_chain(N, 5, 0), a_end, N) != 0)

    held, matched = sum(CLAIMS), sum(CHECKS)
    print(f"\n{matched} of {len(CHECKS)} printed values match exactly; {held} of {len(CLAIMS)} claims hold.")
    if held < len(CLAIMS) or matched < len(CHECKS):
        sys.exit(1)


if __name__ == "__main__":
    main()
