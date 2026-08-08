"""Gate for the exceptional couplings of the frozen band.

The frozen depth on a band block is floor(N/2) for all but finitely many couplings, and
this gate is about the "finitely many": that set is NOT empty, and at a coupling in it one
more mode freezes.

Self-contained: every object is built from scratch here, including one construction from the
full 4^N Liouvillian, so agreement with the other gates is evidence and not a shared bug.

Run:
    python simulations/exceptional_couplings.py            # the default set
    python simulations/exceptional_couplings.py --deep     # + N = 7, the N = 6 whole-block
                                                           #   enumeration and the N = 7 rung
    python simulations/exceptional_couplings.py --slow     # + the two exact counts that cost
                                                           #   about half an hour: N=7 rung 3, N=8
    python simulations/exceptional_couplings.py --rungs    # + the numeric per-rung counts at
                                                           #   N = 8 up to the middle block; it
                                                           #   turns on --deep and --slow too,
                                                           #   about an hour and a half in all
The run prints how many checks it made, so no count is carried in prose anywhere.

What it checks, block by block:

  E1  the builders, from below.  The block generator built sector by sector agrees entry for
      entry with the same block cut out of the dense 4^N Liouvillian, and the split
      L + 4*gbar = C + i*J*A0 has integer C and A0 and is linear in J.
  E2  the rational reformulation.  Writing v = x + i*J*w turns (C + i*J*A0)v = 0 into
      C x - J^2 A0 w = 0 and A0 x + C w = 0, in which only J^2 appears, so at a rational J^2
      the kernel is a rank over QQ.  Checked as an algebraic identity and against the
      complex nullity.
  E3  N = 5, exactly.  Over the whole band, the exact rational frozen count at J^2 = 1,
      J^2 = 3/2 and J^2 = 2.  At J^2 = 3/2 the two diagonal blocks (2,2) and (3,3) carry 3
      where floor(5/2) = 2, and every other block and every other coupling carries 2.
  E4  the singlet subspace.  V = ker(Psi) n ker(S+) n ker(S-) has integer defining matrices,
      both halves of L restrict to it exactly, and no nonzero combination of the floor's
      seeds lies in it, since the seeds at rung p are Phi(corner seeds) and Psi Phi is
      injective there.  So a frozen vector inside V is an EXTRA one.
  E5  the exact determinant polynomial.  q(z) = det(C_V + z A_V) in QQ[z] is even in z, and
      the real roots of q(iJ), counted by Sturm sequences, are the exceptional couplings:
      3 of them at N = 5, 6 at N = 6, 11 at N = 7 (deep).  At N = 5 the smallest is the root
      of 2J^2 - 3, which is E3's J^2 = 3/2 reached by a completely different route.
  E6  the shape of the failure, and the two ladders that give it that shape.  The extra mode
      sits on the diagonal rungs l .. N-l with l = 2 and on no side line, which is an eta
      multiplet that is a SPIN SINGLET.  Then the operator facts the singlet reading rests on,
      exactly and in integers: S+- commutes with L on EVERY block, [S+, S-] = 2*S_z so the two
      generate su(2), and a break-input (one next-nearest bond, which breaks Sigma h Sigma =
      -h) that fails the first of those and leaves the dissipator half untouched.  Together
      they give mult(m=0) >= mult(m=1) on the frozen space, which is what closes open item 4
      of the note: a coupling visible only on a side line does not exist.
  E7  the one coupling the bridge needs.  Nullity of a pencil is constant off a finite set
      and minimal there, and the floor holds everywhere, so a single coupling with nullity
      <= floor(N/2) fixes the generic value.  An exact GF(p) rank at J = 1 supplies one at
      every N here: a rank over GF(p) can only be smaller than the rank over QQ, so a
      measured nullity floor(N/2) mod p bounds the true count from above.
  E8  the WHOLE block, with no subspace assumed.  ker C n ker A0 is the rational, coupling-
      independent part of every kernel, and a congruence removing it leaves an exact
      polynomial for the whole block.  At N = 5 that enumerates the entire band at every
      coupling at once: three exceptional couplings, on (2,2) and (3,3) only, side lines
      carrying none at ANY J; at N = 6 (deep) the block (2,2) carries six.  Both counts are
      the singlet count, and since a singlet root is a root of the whole block's polynomial
      (V meets ker C n ker A0 only in zero), equal counts make the two sets EQUAL.  Every
      positive root is simple wherever this block enumerates, so the count at an exceptional
      coupling is exactly floor(N/2) + 1 rather than at least that, at N = 5 and N = 6.

  E9  a higher rung carries couplings of its own.  A block (p, p) sees only the rungs
      l <= min(p, N-p), so the block (2,2) can carry the rung-2 couplings and no others.  At
      N = 6 the eta-lowest-weight singlets of the block (3,3) carry four couplings of their
      own, exactly, sharing no positive root with the rung-2 polynomial, and each of the four
      raises (3,3) alone, which is the multiplet of a rung-3 lowest weight at eta-spin 0.
      Deep adds the N = 7 counterpart, where that eta-spin is 1/2 and the multiplet is two
      blocks wide.  So the counts of E5 and E8 are rung-2 counts and not the count.
      Then OUTSIDE the band, which no check of this file reached before 2026-08-08, every one
      of them having looped over band(N): at none of the thirteen couplings of N = 5 and N = 6
      does a block with p - q ODD or |p - q| > 2 gain a mode.  The N = 6 half of that measures
      nothing on its own, since no off-band block there crosses the frozen value at ANY J; the
      N = 5 half is the one the positive control certifies, and that control is built from
      those blocks' OWN crossings rather than from the side lines.  And the trap beside it,
      since the exact route DOES return positive roots on an off-band block: those are ordinary
      level crossings, and what says so is exact and has no threshold.  A band block's frozen
      space contains a COUPLING-INDEPENDENT part, ker C n ker A0, of dimension floor(N/2);
      off the band that intersection is zero.
  E10 the two exact counts the other modes are too short for (slow).  The same Sturm route as
      E5, at the rung-3 singlets of the block (3,3) at N = 7 and at the rung-2 singlets at
      N = 8, which is where the exact route runs out.  Kept behind its own flag so that the
      numbers it produces are in the record without lengthening --deep.
  E11 past the exact route: the middle block at N = 8, numerically (rungs).  Detection is the
      finite spectrum of the same pencil rather than a polynomial, and every candidate is
      verified on its own by a SIGN BRACKET of the determinant, so an accepted root is
      certified real.  Nothing certifies that no root was MISSED, so every BLOCK count here is
      a verified LOWER bound and the block says so.  Three things make that word earn its keep
      rather than be asserted.  The expected values in the VALIDATION rows are read out of the
      exact blocks of the same run and none of them is typed here, which is why --rungs turns
      on --deep and --slow; the recorded 40 and 52 are of course literals, being the
      measurement itself and having no exact source to be read from.  The subtraction that
      turns block counts into rung counts is itself validated wherever an exact rung count
      exists, at N = 6 and N = 7, and no such check can exist for rung 4 at N = 8.  And the
      one premise a lower bound silently needs, that no two bracket windows OVERLAP, since
      then one root could certify two candidates, is measured as a margin rather than assumed.
Reading: PROOF_FROZEN_BAND_SO4 Section 5 (the pencil and the bridge) and Section 8
(what is open), ETA_CEILING_REDUCTION (the per-rung certificate this replaces the limit of).
"""
import itertools
import sys
import time

import numpy as np
import scipy.linalg as sla
from sympy import Matrix, Poly, QQ, Rational, factor_list, gcd, interpolate, oo, symbols
from sympy.polys.matrices import DomainMatrix

RUNGS = "--rungs" in sys.argv
# --rungs IMPLIES the two exact tiers, and not for convenience.  E11 validates a numeric
# detector against counts this file proves exactly, and a validation compared against integers
# typed into E11 could not catch a stale one: the run that is the record has to COMPUTE both
# sides.  So the exact counts are published into EXACT_SINGLET / EXACT_BLOCK by the blocks that
# prove them, and E11 reads them from there and fails if one it needs is absent.
DEEP = "--deep" in sys.argv or RUNGS
SLOW = "--slow" in sys.argv or RUNGS

EXACT_SINGLET = {}       # (N, p) -> exact Sturm count on the rung-p singlet space of (p,p)
EXACT_BLOCK = {}         # (N, p, q) -> exact count of the WHOLE block, over QQ
EXACT_POLY = {}          # (N, p) -> the exact polynomial itself, for gcd questions
PRIMES = (998244353, 1004535809)
FAILURES = []
T0 = time.time()


CHECKS = []


def check(label, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {label}" + (f"  ({detail})" if detail else ""))
    CHECKS.append(label)
    if not ok:
        FAILURES.append(label)


# ---------------------------------------------------------------- builders

def subsets(N, p):
    return [sum(1 << i for i in c) for c in itertools.combinations(range(N), p)]


def bits(x, N):
    return [i for i in range(N) if (x >> i) & 1]


def h_sector(N, p, J=1):
    """The open XY chain restricted to the p-excitation sector, occupation-priced."""
    basis = subsets(N, p)
    idx = {b: i for i, b in enumerate(basis)}
    H = np.zeros((len(basis), len(basis)), dtype=np.int64)
    for A in basis:
        for l in range(N - 1):
            for (x, y) in ((l, l + 1), (l + 1, l)):
                if (A >> x) & 1 and not (A >> y) & 1:
                    H[idx[A ^ (1 << x) ^ (1 << y)], idx[A]] += J
    return H, basis


def block_parts(N, p, q, gam=1):
    """C and A0 with L_(p,q) + 4*gbar = C + i*J*A0, both integer at gam = 1.

    C is the rate diagonal shifted by the frozen root, A0 the Hamiltonian part at J = 1."""
    Hp, bp = h_sector(N, p)
    Hq, bq = h_sector(N, q)
    A0 = -(np.kron(Hp, np.eye(len(bq), dtype=np.int64))
           - np.kron(np.eye(len(bp), dtype=np.int64), Hq.T))
    rate = np.array([[2 * gam * bin(A ^ B).count("1") for B in bq] for A in bp]).reshape(-1)
    C = np.diag(-rate + 4 * gam).astype(np.int64)
    return C, A0


def dense_block(N, p, q, J, gam=1):
    """The same block cut out of the FULL 4^N Liouvillian, built from scratch."""
    dim = 1 << N
    H = np.zeros((dim, dim))
    for A in range(dim):
        for l in range(N - 1):
            for (x, y) in ((l, l + 1), (l + 1, l)):
                if (A >> x) & 1 and not (A >> y) & 1:
                    H[A ^ (1 << x) ^ (1 << y), A] += J
    L = -1j * (np.kron(H, np.eye(dim)) - np.kron(np.eye(dim), H.T))
    for l in range(N):
        z = np.array([1.0 if not (A >> l) & 1 else -1.0 for A in range(dim)])
        L += gam * (np.kron(np.diag(z), np.diag(z)) - np.eye(dim * dim))
    rows = [A * dim + B for A in subsets(N, p) for B in subsets(N, q)]
    return L[np.ix_(rows, rows)]


def jw_sign(A, l):
    return -1 if bin(A & ((1 << l) - 1)).count("1") % 2 else 1


def phi_matrix(N, p, q):
    """Phi(rho) = sum_l d+_l rho d_l : (p,q) -> (p+1,q+1)."""
    bp, bq = subsets(N, p), subsets(N, q)
    iP = {b: i for i, b in enumerate(subsets(N, p + 1))}
    iQ = {b: i for i, b in enumerate(subsets(N, q + 1))}
    M = np.zeros((len(iP) * len(iQ), len(bp) * len(bq)), dtype=np.int64)
    for i, A in enumerate(bp):
        for j, B in enumerate(bq):
            for l in range(N):
                if (A >> l) & 1 or (B >> l) & 1:
                    continue
                M[iP[A | (1 << l)] * len(iQ) + iQ[B | (1 << l)], i * len(bq) + j] \
                    += jw_sign(A, l) * jw_sign(B, l)
    return M


def splus_matrix(N, p, q):
    """S+ = sum_l (-1)^l d+_l rho d+_l : (p,q) -> (p+1,q-1)."""
    bp, bq = subsets(N, p), subsets(N, q)
    iP = {b: i for i, b in enumerate(subsets(N, p + 1))}
    iQ = {b: i for i, b in enumerate(subsets(N, q - 1))}
    M = np.zeros((len(iP) * len(iQ), len(bp) * len(bq)), dtype=np.int64)
    for i, A in enumerate(bp):
        for j, B in enumerate(bq):
            for l in range(N):
                if (A >> l) & 1 or not (B >> l) & 1:
                    continue
                sg = jw_sign(A, l) * jw_sign(B ^ (1 << l), l) * (-1) ** l
                M[iP[A | (1 << l)] * len(iQ) + iQ[B ^ (1 << l)], i * len(bq) + j] += sg
    return M


def ladder_plus(N, p, q):
    """S+ : (p,q) -> (p+1,q-1), and the ZERO map where the target block does not exist."""
    if p + 1 > N or q - 1 < 0:
        return np.zeros((0, len(subsets(N, p)) * len(subsets(N, q))), dtype=np.int64)
    return splus_matrix(N, p, q)


def ladder_minus(N, p, q):
    """S- : (p,q) -> (p-1,q+1), as the transpose of the S+ that lands in this block."""
    if p - 1 < 0 or q + 1 > N:
        return np.zeros((0, len(subsets(N, p)) * len(subsets(N, q))), dtype=np.int64)
    return ladder_plus(N, p - 1, q + 1).T


def general_parts(N, p, q, bonds, gammas):
    """(D, A) with L_(p,q) = D + i*A, for an arbitrary bond list [(x,y,J)] and rate profile.

    The gate's block_parts is the uniform case with the frozen root already subtracted; this is
    the unshifted generator at an arbitrary profile, and it exists so that the ladder checks can
    ask Lemma 2.1's stronger question (any gamma_l) and Lemma 2.3's break-input (a bond that
    joins two sites of the SAME parity)."""
    def hop(r):
        # The JORDAN-WIGNER STRING is the whole reason this is not h_sector with a bond list.
        # On an open chain with nearest-neighbour bonds the strings cancel and h_sector may omit
        # them; across a bond that skips sites they do not, and omitting them makes EVERY added
        # bond break the ladder by 2*J whatever its parity, so a break-input built on the naive
        # hop tests "a bond was added" and not "Sigma h Sigma = -h failed".  With the string the
        # measured verdicts are the ones PROOF_FROZEN_BAND_SO4 Lemma 2.3 predicts and the ones
        # the `sideways_spin_ladder` arc recorded, the ring among them.
        basis = subsets(N, r)
        idx = {b: i for i, b in enumerate(basis)}
        H = np.zeros((len(basis), len(basis)))
        for A in basis:
            for (u, v, J) in bonds:
                for (x, y) in ((u, v), (v, u)):
                    if x == y:
                        # an ON-SITE term, (l, l, h): diagonal, and it breaks Sigma h Sigma = -h
                        # because Sigma-oddness forces a zero diagonal.  It is one of the five
                        # geometries the arc names, and the first version of this builder could
                        # not express it: the hop condition below is never true for x == y, so a
                        # self-bond was silently a no-op.
                        if (A >> x) & 1:
                            H[idx[A], idx[A]] += J
                        break
                    if (A >> x) & 1 and not (A >> y) & 1:
                        lo, hi = min(x, y), max(x, y)
                        between = bin(A & (((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1))).count("1")
                        H[idx[A ^ (1 << x) ^ (1 << y)], idx[A]] += J * (-1) ** between
        return H

    Hp, Hq = hop(p), hop(q)
    A = -(np.kron(Hp, np.eye(len(subsets(N, q)))) - np.kron(np.eye(len(subsets(N, p))), Hq.T))
    D = np.diag([-2.0 * sum(gammas[l] for l in bits(X ^ Y, N))
                 for X in subsets(N, p) for Y in subsets(N, q)])
    return D, A


def band(N):
    return sorted({(p, q) for p in range(N + 1) for q in range(N + 1)
                   if abs(p - q) in (0, 2)} - {(0, 0), (N, N)})


# ------------------------------------------------- exact ranks and nullities

def rational_nullity(C, A0, J2):
    """The exact frozen count at the coupling J^2 = J2, as a rank over QQ.

    v = x + i*J*w turns (C + i*J*A0)v = 0 into the rational pair below, so the real form has
    twice the complex nullity and only J^2 enters."""
    d = C.shape[0]
    R = Matrix.zeros(2 * d, 2 * d)
    for i in range(d):
        R[i, i] = Rational(int(C[i, i]))
        R[d + i, d + i] = Rational(int(C[i, i]))
    for i in range(d):
        for j in range(d):
            if A0[i, j]:
                R[i, d + j] = -J2 * Rational(int(A0[i, j]))
                R[d + i, j] = Rational(int(A0[i, j]))
    dm = DomainMatrix.from_Matrix(R).convert_to(QQ)
    nul = 2 * d - len(dm.rref()[1])
    assert nul % 2 == 0, "the real form of a complex kernel has even dimension"
    return nul // 2


def modpow(b, e, p):
    acc, b = 1, b % p
    while e:
        if e & 1:
            acc = acc * b % p
        b = b * b % p
        e >>= 1
    return acc


def sqrt_minus_one(p):
    for t in range(2, 200):
        r = modpow(t, (p - 1) // 4, p)
        if r * r % p == p - 1:
            return r
    raise RuntimeError(f"no square root of -1 mod {p}")


def rank_mod_p(re, im, p):
    """Exact rank over GF(p) with i sent to a square root of -1.  A rank over GF(p) can only
    be SMALLER than the rank over QQ, so a nullity read here bounds the true nullity from
    above, which is the direction a ceiling needs."""
    r = sqrt_minus_one(p)
    a = (((re % p) + r * (im % p)) % p).astype(np.int64)
    rows, cols = a.shape
    rank = 0
    for col in range(cols):
        if rank >= rows:
            break
        nz = np.nonzero(a[rank:, col])[0]
        if nz.size == 0:
            continue
        piv = rank + int(nz[0])
        if piv != rank:
            a[[rank, piv]] = a[[piv, rank]]
        inv = modpow(int(a[rank, col]), p - 2, p)
        a[rank, col:] = (a[rank, col:] * inv) % p
        below = a[rank + 1:, col]
        hit = np.nonzero(below)[0]
        if hit.size:
            a[rank + 1 + hit, col:] = (a[rank + 1 + hit, col:]
                                       - below[hit, None] * a[rank, col:][None, :]) % p
        rank += 1
    return rank


def modp_nullity(C, A0, J):
    re, im = C, (J * A0)
    return C.shape[0] - max(rank_mod_p(re, im, p) for p in PRIMES)


# ------------------------------------------------------- the singlet subspace

def singlet_basis(N, p):
    """An integer basis of V = ker(Psi) n ker(S+) n ker(S-) in the block (p,p), as columns."""
    stack = np.vstack([phi_matrix(N, p - 1, p - 1).T,
                       splus_matrix(N, p, p),
                       splus_matrix(N, p - 1, p + 1).T])
    ns = Matrix(stack.tolist()).nullspace()
    return Matrix.hstack(*ns) if ns else Matrix.zeros(stack.shape[1], 0)


def restrict(M, B):
    """The matrix of an operator that preserves span(B), in the basis B, exactly."""
    A = Matrix(np.asarray(M).astype(int).tolist())
    AB = A * B
    _, piv = B.T.rref()                      # pivot columns of B^T are independent rows of B
    X = B[list(piv), :].solve(AB[list(piv), :])
    ok = (B * X - AB).is_zero_matrix
    return X, ok


def det_poly(CV, AV, z):
    """q(z) = det(C_V + z A_V) in QQ[z], by interpolation from exact rational determinants."""
    deg = CV.shape[0]
    pts = []
    for k in range(deg + 1):
        zk = Rational(k - deg // 2)
        pts.append((zk, DomainMatrix.from_Matrix(CV + zk * AV).convert_to(QQ).det()))
    return Poly(interpolate(pts, z), z)


def exceptional_polynomial(N, p, J):
    """The exact polynomial in J whose real roots are the exceptional couplings seen by the
    eta-lowest-weight singlets of the block (p,p), together with dim V."""
    z = symbols('z')
    C, A0 = block_parts(N, p, p)
    B = singlet_basis(N, p)
    CV, okc = restrict(C, B)
    AV, oka = restrict(A0, B)
    q = det_poly(CV, AV, z)
    co = list(reversed(q.all_coeffs()))
    even = all(c == 0 for k, c in enumerate(co) if k % 2)
    re = sum(c * (-1) ** (k // 2) * J ** k for k, c in enumerate(co) if k % 2 == 0)
    im = sum(c * (-1) ** (k // 2) * J ** k for k, c in enumerate(co) if k % 2 == 1)
    g = Poly(re, J) if even else Poly(gcd(Poly(re, J), Poly(im, J)), J)
    return g, B.cols, q, (okc and oka), even


# ------------------------------------------------------------------- E1

def e1_builders():
    print("E1  the builders, from below")
    worst = 0.0
    # (5,2,3) and (5,1,2) are OFF the band, p - q odd.  They are here because E9(g) makes a
    # statement about the roots those blocks carry, and a builder check that only ever saw band
    # blocks cannot support it.
    for (N, p, q) in [(4, 2, 2), (4, 3, 1), (5, 2, 2), (5, 3, 3), (5, 2, 0),
                      (5, 2, 3), (5, 1, 2)]:
        C, A0 = block_parts(N, p, q)
        for J in (1, 2):
            mine = C - 4 * np.eye(C.shape[0], dtype=np.int64) + 1j * J * A0
            worst = max(worst, float(np.max(np.abs(mine - dense_block(N, p, q, J)))))
    check("(a) block generator = the dense 4^N Liouvillian, cut to the block",
          worst == 0.0, f"max |difference| = {worst}")
    C, A0 = block_parts(5, 2, 2)
    lin = all(np.array_equal((C + 1j * J * A0).imag, J * A0) for J in (1, 2, 3))
    integral = C.dtype == np.int64 and A0.dtype == np.int64
    check("(b) L + 4*gbar = C + i*J*A0 with C, A0 integer and the J dependence linear",
          lin and integral)
    diag_only = np.array_equal(C, np.diag(np.diag(C)))
    vals = sorted(set(np.diag(C).tolist()))
    check("(c) C is diagonal, entries 4*gamma*(1 - k) on the disagreement k",
          diag_only and vals == [-4, 0, 4], f"entries {vals}")
    # the chain is bipartite, so each sector's hopping is chiral-ODD: with the site parity
    # G_p = diag(prod over the occupied sites of (-1)^site), P = G_p (x) G_q fixes C and flips
    # A0.  det P = +-1, so det(C + z A0) = det(P (C + z A0) P) = det(C - z A0) and the
    # determinant is EVEN in z, on every block and with no interpolation.  E5 measures the
    # evenness case by case; this is the reason.
    worst_p, all_inv = 0, True
    for (N, p, q) in [(5, 2, 2), (5, 3, 3), (6, 2, 2), (6, 3, 3), (6, 2, 4), (7, 3, 3)]:
        C, A0 = block_parts(N, p, q)
        g = [np.diag([(-1) ** sum(l for l in bits(A, N)) for A in subsets(N, r)])
             for r in (p, q)]
        P = np.kron(g[0], g[1])
        worst_p = max(worst_p, int(np.max(np.abs(P @ C @ P - C))),
                      int(np.max(np.abs(P @ A0 @ P + A0))))
        all_inv = all_inv and np.array_equal(P @ P, np.eye(P.shape[0], dtype=np.int64))
    check("(d) the chiral involution P = G_p (x) G_q fixes C and flips A0 entry for entry, so "
          "det(C + z A0) is EVEN in z on every block, side lines included",
          worst_p == 0 and all_inv, f"worst entry {worst_p}, P^2 = I on all: {all_inv}")


# ------------------------------------------------------------------- E2

def e2_reformulation():
    print("E2  the rational reformulation, only J^2 enters")
    rng = np.random.default_rng(20260728)
    C, A0 = block_parts(5, 2, 2)
    d = C.shape[0]
    x, w = rng.normal(size=d), rng.normal(size=d)
    J = 1.3
    lhs = (C + 1j * J * A0) @ (x + 1j * J * w)
    rhs = (C @ x - J * J * (A0 @ w)) + 1j * J * (A0 @ x + C @ w)
    check("(a) (C + iJ A0)(x + iJ w) = (C x - J^2 A0 w) + iJ (A0 x + C w)",
          float(np.max(np.abs(lhs - rhs))) < 1e-10)
    for (N, p, J2) in [(4, 2, Rational(1)), (5, 2, Rational(1)), (5, 2, Rational(3, 2))]:
        C, A0 = block_parts(N, p, p)
        exact = rational_nullity(C, A0, J2)
        s = np.linalg.svd(C + 1j * float(J2) ** 0.5 * A0, compute_uv=False)
        num = int(np.sum(s < 1e-9 * max(s[0], 1.0)))
        check(f"(b) N={N} ({p},{p}) J^2={J2}: exact rational count = complex nullity",
              exact == num, f"{exact} = {num}")
    # C is gamma times an integer matrix and A0 carries no gamma, so C + iJ A0 =
    # gamma*(C_1 + i(J/gamma) A0): the exceptional set is a set of RATIOS J/gamma and every
    # number tabulated anywhere is at gamma = 1.  Checked exactly, not by scaling the output.
    ratio_ok, moved = True, {}
    for gam in (2, 3):
        C, A0 = block_parts(5, 2, 2, gam=gam)
        C1, A01 = block_parts(5, 2, 2, gam=1)
        structural = np.array_equal(C, gam * C1) and np.array_equal(A0, A01)
        here = rational_nullity(C, A0, Rational(3, 2) * gam ** 2)
        there = rational_nullity(C, A0, Rational(3, 2))
        moved[gam] = (here, there)
        ratio_ok = ratio_ok and structural and here == 3 and there == 2
    check("(c) only the ratio J/gamma enters: at gamma = 2 and 3 the N=5 point sits at "
          "J^2 = (3/2)*gamma^2 and no longer at 3/2, exactly", ratio_ok,
          f"(nullity at (3/2)g^2, at 3/2) = {moved}")


# ------------------------------------------------------------------- E3

def e3_n5_exact():
    print("E3  N = 5, the whole band, exactly")
    N, m = 5, 2
    couplings = [Rational(1), Rational(3, 2), Rational(2)]
    extra = []
    table_ok = True
    for (p, q) in band(N):
        C, A0 = block_parts(N, p, q)
        counts = [rational_nullity(C, A0, t) for t in couplings]
        if (p, q) in {(2, 2), (3, 3)}:
            ok = counts == [m, m + 1, m]
            if counts[1] == m + 1:
                extra.append((p, q))
        else:
            ok = counts == [m, m, m]
        table_ok = table_ok and ok
        if not ok:
            check(f"    block ({p},{q}) unexpected", False, f"{counts}")
    check("(a) the exact table over the band at J^2 = 1, 3/2, 2", table_ok)
    check("(b) at J^2 = 3/2 the count exceeds floor(N/2) on exactly two blocks",
          sorted(extra) == [(2, 2), (3, 3)], f"{sorted(extra)}")
    check("(c) so the exceptional set is NOT empty, and the ceiling is a generic-J statement",
          len(extra) > 0)


# ------------------------------------------------------------------- E4

def e4_singlets():
    print("E4  the singlet subspace, and why a frozen vector in it is an extra one")
    dims = {}
    for N in (5, 6):
        C, A0 = block_parts(N, 2, 2)
        B = singlet_basis(N, 2)
        dims[N] = B.cols
        _, okc = restrict(C, B)
        _, oka = restrict(A0, B)
        check(f"(a) N={N}: both halves of L restrict to V exactly", okc and oka,
              f"dim V = {B.cols} of {C.shape[0]}")
    check("(b) V is a proper subspace at both N", all(v > 0 for v in dims.values())
          and dims[5] == 35 and dims[6] == 84, f"{dims}")
    # the seeds are Phi(corner seeds) and Psi Phi is injective there, so none lies in ker Psi
    for N in (5, 6):
        C1, A01 = block_parts(N, 1, 1)
        s = np.linalg.svd(C1 + 1j * A01, compute_uv=False)
        k = int(np.sum(s < 1e-9 * max(s[0], 1.0)))
        _, _, vh = np.linalg.svd(C1 + 1j * A01)
        F = vh[len(s) - k:].conj().T                       # the corner frozen space
        climbed = phi_matrix(N, 1, 1) @ F
        back = phi_matrix(N, 1, 1).T @ climbed             # Psi Phi on the seeds
        sv = np.linalg.svd(back, compute_uv=False)
        check(f"(c) N={N}: Psi Phi is injective on the corner frozen space, so no seed "
              f"combination is eta-lowest-weight at rung 2",
              k == N // 2 and float(sv[-1]) > 1e-6, f"{k} seeds, smallest sv {sv[-1]:.3f}")


# ------------------------------------------------------------------- E5

def e5_polynomial():
    print("E5  the exact determinant polynomial and its real roots")
    J = symbols('J')
    expect = {5: 3, 6: 6, 7: 11}
    for N in ([5, 6] + ([7] if DEEP else [])):
        g, dimV, q, restrict_ok, even = exceptional_polynomial(N, 2, J)
        check(f"(a) N={N}: q(z) = det(C_V + z A_V) is even in z, so q(iJ) is real",
              even and restrict_ok, f"dim V = {dimV}, deg q = {q.degree()}")
        npos = g.count_roots(0, oo) - (1 if g.eval(0) == 0 else 0)
        EXACT_SINGLET[(N, 2)] = npos
        EXACT_POLY[(N, 2)] = g
        check(f"(b) N={N}: exact count of exceptional couplings J > 0, by Sturm sequences",
              npos == expect[N], f"{npos}, expected {expect[N]}")
        roots = sorted(float(r) for r in g.real_roots() if r.is_real and float(r) > 1e-9)
        # every exact root must make the FULL block singular, which the subspace did not assume
        C, A0 = block_parts(N, 2, 2)
        worst = 0.0
        for r in roots:
            s = np.linalg.svd(C + 1j * r * A0, compute_uv=False)      # descending
            worst = max(worst, float(s[-(N // 2 + 1)] / s[0]))        # the (m+1)-th SMALLEST
        check(f"(c) N={N}: at every exact root the full block is singular beyond the floor, "
              f"so it carries AT LEAST floor(N/2)+1 frozen modes (a numeric read of a lower "
              f"bound; the exact count is E8's)", worst < 1e-9,
              f"largest surplus singular value {worst:.2e}")
        if N == 5:
            fac = [f for f, _ in factor_list(g.as_expr())[1]
                   if Poly(f, J).degree() == 2 and Poly(f, J).count_roots(0, oo)]
            match = bool(fac) and Poly(fac[0], J).as_expr().equals(2 * J ** 2 - 3)
            check("(d) N=5: the smallest exceptional coupling is the root of 2J^2 - 3, which "
                  "is E3's J^2 = 3/2 by a different route", match,
                  f"{fac[0] if fac else None}")
        # (a) to (c) repeat per N; (d) to (g) are per-N SPECIALS and so must be distinct.  They
        # were duplicated and out of order until 2026-08-08: N = 5 and N = 7 both printed a
        # (d), and N = 6 printed (f) before (e).
        if N == 7:
            check("(g) N=7: the smallest exceptional coupling is the value the note tabulates",
                  abs(roots[0] - 0.952056678) < 1e-9, f"{roots[0]:.12f}")
        if N == 6:
            check("(e) N=6: the smallest exceptional coupling is the value the note tabulates",
                  abs(roots[0] - 0.749042443688) < 1e-11, f"{roots[0]:.12f}")
            degs = sorted(Poly(f, J).degree() for f, _ in factor_list(g.as_expr())[1]
                          if any(float(r) > 1e-9 for r in Poly(f, J).real_roots() if r.is_real))
            check("(f) N=6: BOTH factors carrying a coupling have degree 12, so every one of "
                  "the six is algebraic of that degree and no small closed form appears",
                  degs == [12, 12], f"degrees {degs}")


# ------------------------------------------------------------------- E6

def e6_shape():
    print("E6  the shape of the failure: an l = 2 multiplet that is a spin singlet")
    J = symbols('J')
    g, _, _, _, _ = exceptional_polynomial(6, 2, J)
    roots = sorted(float(r) for r in g.real_roots() if r.is_real and float(r) > 1e-9)
    J0 = roots[-1]
    counts = {}
    for (p, q) in band(6):
        C, A0 = block_parts(6, p, q)
        s = np.linalg.svd(C + 1j * J0 * A0, compute_uv=False)
        counts[(p, q)] = int(np.sum(s < 1e-9 * max(s[0], 1.0)))
    raised = sorted(k for k, v in counts.items() if v > 3)
    check("(a) N=6: at this coupling the extra mode occupies the diagonal rungs 2 .. N-2 and "
          "nothing else (one coupling only; the all-J statement is E8, at N = 5)",
          raised == [(2, 2), (3, 3), (4, 4)], f"raised on {raised}")
    check("(b) so, a lowest weight at rung l spanning N - 2l + 1 blocks, it is an eta "
          "multiplet seeded at l = 2 with eta-spin N/2 - 2: the arithmetic of (a)",
          len(raised) == 6 - 2 * 2 + 1)
    check("(c) no side line is raised: every band block with p != q keeps the generic "
          "floor(N/2). The SINGLET reading of this is (d) to (h) below, which compute the "
          "premise this check used to carry in its own label and never verified",
          all(counts[(p, q)] == 3 for (p, q) in band(6) if p != q))

    # (d) to (f): the premise, computed.  Until 2026-08-08 the label of (c) read "since S+
    # commutes with L", and no commutator was computed anywhere in this file: S+ entered it
    # only through singlet_basis.  The three checks below are what that sentence asserts, and
    # (f) is the input built to make (d) fail, because a guard nobody can break is worth
    # nothing.  They also supply what open item 4 of the note needs; see the note.
    worst, npairs = 0, 0
    for Nc in (4, 5, 6):
        for p in range(Nc + 1):
            for q in range(Nc + 1):
                C, A0 = block_parts(Nc, p, q)
                for (M, tp, tq) in ((ladder_plus(Nc, p, q), p + 1, q - 1),
                                    (ladder_minus(Nc, p, q), p - 1, q + 1)):
                    if M.shape[0] == 0:
                        continue
                    Ct, At = block_parts(Nc, tp, tq)
                    worst = max(worst, int(np.max(np.abs(M @ C - Ct @ M))),
                                int(np.max(np.abs(M @ A0 - At @ M))))
                    npairs += 1
    check("(d) N = 4, 5, 6, EVERY block and both ladders: S+-.L = L.S+- entry for entry. The "
          "shift 4*gbar is a scalar, so C and A0 carry the whole question, and both are "
          "integer, so this is an exact 0 and not a tolerance",
          worst == 0, f"worst |entry| = {worst} over {npairs} block pairs")

    worst2, nblk = 0, 0
    for Nc in (4, 5, 6):
        for p in range(Nc + 1):
            for q in range(Nc + 1):
                d = len(subsets(Nc, p)) * len(subsets(Nc, q))
                lhs = np.zeros((d, d), dtype=np.int64)
                if p + 1 <= Nc and q - 1 >= 0:
                    lhs -= ladder_minus(Nc, p + 1, q - 1) @ ladder_plus(Nc, p, q)
                if p - 1 >= 0 and q + 1 <= Nc:
                    lhs += ladder_plus(Nc, p - 1, q + 1) @ ladder_minus(Nc, p, q)
                worst2 = max(worst2, int(np.max(np.abs(lhs - (p - q) * np.eye(d, dtype=np.int64)))))
                nblk += 1
    # This is what makes "integer spin" mean anything.  "S+ shifts p - q by 2" alone does not
    # exclude a module living at one weight and nothing else; the algebra does, since such a
    # module would make [S+, S-] = 2*S_z read 0 = 2.  What FOLLOWS from it, that the weight
    # multiplicities of a finite-dimensional su(2) module obey mult(0) >= mult(1), is
    # representation theory and is deliberately NOT claimed in this label: it is measured
    # separately in (h).  A label that carries its own consequence is the defect this block
    # was written to repair.
    # Count the trivial slots rather than deriving them from a formula: only the blocks where
    # BOTH ladder products are absent are 0 == 0, and that is the corners, not every p = q.  A
    # p = q block with both products PRESENT is the Cartan relation cancelling, which is the
    # most substantive cell in the check.  The first version of this annotation said 18 by a
    # formula that never looked at a block; it is 6.
    trivial = sum(1 for Nc in (4, 5, 6) for p in range(Nc + 1) for q in range(Nc + 1)
                  if not (p + 1 <= Nc and q - 1 >= 0) and not (p - 1 >= 0 and q + 1 <= Nc))
    check("(e) N = 4, 5, 6, every block: [S+, S-] = (p - q)*I = 2*S_z, exactly, so the two "
          "ladders generate su(2) and the weights are the block labels",
          worst2 == 0, f"worst |entry| = {worst2} over {nblk} blocks, of which {trivial} are "
          f"the corners where both ladder products are absent")

    def ladder_break(Nb, bonds, gam):
        """(dissipator half, turning half) of the S+- commutator, at an arbitrary profile."""
        wd = wh = 0.0
        for p in range(Nb + 1):
            for q in range(Nb + 1):
                D, A = general_parts(Nb, p, q, bonds, gam)
                for (M, tp, tq) in ((ladder_plus(Nb, p, q), p + 1, q - 1),
                                    (ladder_minus(Nb, p, q), p - 1, q + 1)):
                    if M.shape[0] == 0:
                        continue
                    Dt, At = general_parts(Nb, tp, tq, bonds, gam)
                    wd = max(wd, float(np.max(np.abs(M @ D - Dt @ M))))
                    wh = max(wh, float(np.max(np.abs(M @ A - At @ M))))
        return wd, wh

    Nb, gam = 5, [0.25, 1.5, 0.75, 2.0, 0.5]
    nn = [(l, l + 1, 0.5 + 0.25 * l) for l in range(Nb - 1)]
    # The break-input is a PARITY test, not an "a bond was added" test, and the difference is
    # the whole content of Lemma 2.3.  Both kinds of extra bond are run: same-parity ones must
    # break Sigma h Sigma = -h and opposite-parity ones must NOT.  The ring is run at odd and
    # even N for the same reason, and reproduces the verdicts the `sideways_spin_ladder` arc
    # recorded on 2026-08-07 (breaks at odd N, holds at even).
    # Keyed on the BONDS, never on the printed label: the first version selected the conjuncts
    # with k.startswith("same") / ("opposite"), so renaming a print string emptied both lists
    # and made two all(...) vacuously true.  A guard that a rename can silence is not a guard.
    nn_break = ladder_break(Nb, nn, gam)
    # sites are 0 .. Nb-1; same parity: (0,2), (1,3), (2,4), (0,4); opposite: (0,3), (1,4)
    extras = [(0, 2, 0.9), (1, 3, 0.9), (2, 4, 0.9), (0, 4, 0.9), (0, 3, 0.9), (1, 4, 0.9)]
    same, opp, broke = [], [], {"nearest neighbour": nn_break}
    for (x, y, Jb) in extras:
        v = ladder_break(Nb, nn + [(x, y, Jb)], gam)
        broke[f"bond ({x},{y})"] = v
        ((same if (x - y) % 2 == 0 else opp)).append(v[1])
    onsite = ladder_break(Nb, nn + [(2, 2, 0.7)], gam)
    broke["on-site (2,2)"] = onsite
    rings = {}
    for Nr in (4, 5, 6):
        rings[Nr] = ladder_break(Nr, [(l, (l + 1) % Nr, 1.0) for l in range(Nr)], [1.0] * Nr)[1]
    # The dissipator column is bond-INDEPENDENT by construction: general_parts builds D from
    # the rate profile and the disagreement set alone, and S+- leaves that set fixed, which is
    # Lemma 2.1 entry-wise.  So the two rows agreeing there is not a contrast between them, and
    # the label must not draw one.  What the break-input tests is the TURNING half only.  The
    # threshold is exact (> 0): a small next-nearest J would break the physics by a small
    # number, and a fixed cutoff would then fail on a correct object.
    # The arc `sideways_spin_ladder` (compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs,
    # opened 2026-08-07) ran this break the day before, on five geometries: a next-nearest
    # bond, an on-site potential, the ring at odd N (it holds at even N) and the star.  This
    # is one of those five, kept here because this gate must be able to fail on its own.
    check("(f) the break-input, so (d) is a guard that CAN fail, and it tests the PARITY rather "
          "than the presence of a bond: on the open chain with nearest-neighbour bonds the "
          "turning half commutes exactly; a SAME-parity bond breaks Sigma h Sigma = -h and must "
          "break it; an OPPOSITE-parity one must not, and is exactly 0; an ON-SITE term breaks "
          "it, Sigma-oddness forcing a zero diagonal; and the ring breaks at odd N and holds at "
          "even (PROOF_FROZEN_BAND_SO4 Lemma 2.3, and four of the five verdicts the "
          "sideways_spin_ladder arc recorded). The gamma profile is non-uniform throughout, but "
          "that is the DISSIPATOR's business: its column is bond-independent by construction "
          "and is printed, not contrasted",
          nn_break[1] == 0.0 and all(v > 0 for v in same) and all(v == 0.0 for v in opp)
          and onsite[1] > 0 and rings[5] > 0 and rings[4] == 0.0 and rings[6] == 0.0,
          f"turning half: NN {nn_break[1]}, same-parity {same}, opposite-parity {opp}, on-site "
          f"{onsite[1]}, rings {rings}; dissipator column (bond-free) "
          f"{sorted({v[0] for v in broke.values()})}")

    # (g) the OTHER ladder, which the note's item-4 argument needs and which had the same gap
    # (h) closes on the spin side: the step from a diagonal block to the MIDDLE one is Phi, not
    # S+, and for ker(L - lambda) to be an sl(2)_eta module BOTH Phi and its adjoint Psi must
    # commute with L.  PROOF_FROZEN_BAND_SO4 proves [K, Phi] = 0 (Lemma 2.2) and [Kcount, Psi]
    # = 0 (Theorem 4.1, the integer count, not L).  Neither is [L, Psi] = 0.
    worst3, npair3 = 0, 0
    for Nc in (4, 5, 6):
        for p in range(Nc + 1):
            for q in range(Nc + 1):
                C, A0 = block_parts(Nc, p, q)
                legs = []
                if p + 1 <= Nc and q + 1 <= Nc:
                    legs.append((phi_matrix(Nc, p, q), p + 1, q + 1))
                if p - 1 >= 0 and q - 1 >= 0:
                    legs.append((phi_matrix(Nc, p - 1, q - 1).T, p - 1, q - 1))
                for (M, tp, tq) in legs:
                    Ct, At = block_parts(Nc, tp, tq)
                    worst3 = max(worst3, int(np.max(np.abs(M @ C - Ct @ M))),
                                 int(np.max(np.abs(M @ A0 - At @ M))))
                    npair3 += 1
    check("(g) N = 4, 5, 6, every block: the eta ladder Phi AND its adjoint Psi commute with L "
          "entry for entry, so ker(L - lambda) is an sl(2)_eta module too. Every eta multiplet "
          "spans the rungs N/2 - j to N/2 + j and so contains the MIDDLE rung at every j, which "
          "is the step from a diagonal block to (floor(N/2), floor(N/2))",
          worst3 == 0, f"worst |entry| = {worst3} over {npair3} block pairs")

    # (h) the consequence, MEASURED rather than carried in a label.  mult(0) >= mult(1) is
    # representation theory once (d), (e) hold; what a gate can do is read it off the object.
    # Both the generic couplings and the exceptional ones, since the exceptional ones are the
    # whole point: at each, the diagonal block at a given p + q must carry at least as many
    # frozen modes as the side line beside it.
    viol, seen_pairs = [], 0
    for Jx in [0.6, 1.3172583, 2.4] + roots:
        for p in range(1, 6):
            nul = {}
            for b in ((p, p), (p + 1, p - 1), (p - 1, p + 1)):
                C, A0 = block_parts(6, *b)
                s = np.linalg.svd(C + 1j * Jx * A0, compute_uv=False)
                nul[b] = int(np.sum(s < 1e-9 * max(s[0], 1.0)))
            seen_pairs += 1
            if nul[(p, p)] < max(nul[(p + 1, p - 1)], nul[(p - 1, p + 1)]):
                viol.append((round(Jx, 9), p, nul))
    # TWO controls, and the first is the one that matters, because two earlier versions of this
    # comment asserted there was none.  There IS an object that breaks the inequality: add the
    # bond (1,3) or (2,4) to the chain and at N = 6 the side lines carry 2 where the diagonal
    # carries 1, at p = 1 and p = 5.  That is mult(m=1) > mult(m=0), on an object where the
    # theorem has nothing to say because [L, S+-] is no longer 0.
    #
    # Say WHY carefully, because the first version of this comment said it wrong and said it in
    # the same shape the file corrects elsewhere.  This control is built in the pricing
    # block_parts uses, by OCCUPATION, and there an added bond breaks the ladder whatever its
    # parity.  The check COMPUTES that for the two bonds it uses and requires it (the ladder
    # residual is in the detail line).  Read off-gate on 2026-08-08 and recorded as an
    # observation, not gated: (0,2), (1,3), (2,4), (0,4), (0,3), (1,4), (2,5) all give exactly
    # 2.0 at N = 6, J = 1.  Parity does no work in this pricing, and that is exactly what
    # general_parts exists to fix for E6(f).  So the honest statement is: any added bond makes
    # this a non-su(2) object, and on two of them the inequality visibly fails.  Which two is
    # not predicted here; they were found.
    #
    # Not every mutation does it.  The bond (0,2) collapses the sides to 0 and breaks nothing
    # visible; the STAR breaks Sigma-oddness and satisfies the inequality anyway, priced either
    # way: by occupation (no Jordan-Wigner string) N = 6 reads diag 12/11/12/11/12 against sides
    # 5/0/1/0/5, and by order (the fermionic hop general_parts builds) diag 12/13/17/13/12
    # against sides 7/7/7/7/7, the same split PROOF_FROZEN_BAND_SO4 section 2 names as 2/5/9
    # against 4/7/11.  So a mutation that fails to break this proves nothing, and the specific
    # one that does is the control.
    #
    # The second control is the PAIRING, which is what could be silently wrong: the whole
    # content is that (p,p) and (p+1,p-1) sit at the same p + q and so are two weights of ONE
    # module.  Pair the diagonal with the next DIAGONAL block instead, the eta ladder's step and
    # a different p + q, and the same inequality must hold generically and fail at the
    # exceptional couplings.
    def wrong_pairing(Js):
        """Violations of the same inequality across the ETA step, a DIFFERENT p + q."""
        out = []
        for Jx in Js:
            for p in range(1, 5):
                n = []
                for tgt in ((p, p), (p + 1, p + 1)):
                    C, A0 = block_parts(6, *tgt)
                    s = np.linalg.svd(C + 1j * Jx * A0, compute_uv=False)
                    n.append(int(np.sum(s < 1e-9 * max(s[0], 1.0))))
                if n[0] < n[1]:
                    out.append((round(Jx, 9), p, n[0], n[1]))
        return out

    wrong = wrong_pairing(roots)
    wrong_generic = wrong_pairing([0.6, 1.3172583, 2.4])

    def broken_object_violates(extra):
        """Does THIS extra bond break mult(0) >= mult(1)?  For (1,3) and (2,4) at N = 6 it does.

        NOT a parity statement, and an earlier version of this line said it was: (0,2) has the
        same parity as (1,3) and produces no violation at all.  What is true is that any added
        bond makes [L, S+-] nonzero in this pricing, so the su(2) theorem stops protecting the
        inequality; whether it then breaks is a fact about the particular bond, and this
        function is how the two that do were found rather than predicted.

        Built in the SAME pricing as block_parts, by occupation, because that is the object
        (h) checks.  The distinction is load-bearing here and not pedantry: priced by ORDER,
        with the Jordan-Wigner string general_parts carries, the same bond does NOT produce a
        violation, so a break-input taken from the other builder would have proved nothing."""
        N6, bonds = 6, [(l, l + 1, 1.0) for l in range(5)] + [extra]

        def hop(r):
            basis = subsets(N6, r)
            idx = {b: i for i, b in enumerate(basis)}
            H = np.zeros((len(basis), len(basis)))
            for A in basis:
                for (u, v, Jb) in bonds:
                    for (x, y) in ((u, v), (v, u)):
                        if (A >> x) & 1 and not (A >> y) & 1:
                            H[idx[A ^ (1 << x) ^ (1 << y)], idx[A]] += Jb
            return H

        def nullity(p, q):
            Hp, Hq = hop(p), hop(q)
            A = -(np.kron(Hp, np.eye(len(subsets(N6, q))))
                  - np.kron(np.eye(len(subsets(N6, p))), Hq.T))
            rate = np.array([[2 * bin(X ^ Y).count("1") for Y in subsets(N6, q)]
                             for X in subsets(N6, p)]).reshape(-1)
            s = np.linalg.svd(np.diag(-rate + 4.0) + 1j * A, compute_uv=False)
            return int(np.sum(s < 1e-9 * max(s[0], 1.0)))

        # The PREMISE of the control, computed and not asserted: on this mutated object the
        # ladder no longer commutes with L, so ker(L - lambda) is not an su(2) module and the
        # theorem is silent.  Without this the control would be a mutation with no stated
        # reason to expect anything, and the label would be carrying its cause again.
        ladder_resid = 0.0
        for p in range(N6 + 1):
            for q in range(N6 + 1):
                A = -(np.kron(hop(p), np.eye(len(subsets(N6, q))))
                      - np.kron(np.eye(len(subsets(N6, p))), hop(q).T))
                for (M, tp, tq) in ((ladder_plus(N6, p, q), p + 1, q - 1),
                                    (ladder_minus(N6, p, q), p - 1, q + 1)):
                    if M.shape[0] == 0:
                        continue
                    At = -(np.kron(hop(tp), np.eye(len(subsets(N6, tq))))
                           - np.kron(np.eye(len(subsets(N6, tp))), hop(tq).T))
                    ladder_resid = max(ladder_resid, float(np.max(np.abs(M @ A - At @ M))))

        out = []
        for p in range(1, 6):
            d, up, dn = nullity(p, p), nullity(p + 1, p - 1), nullity(p - 1, p + 1)
            if d < max(up, dn):                 # report BOTH sides: the condition maxes them
                out.append((p, d, up, dn))
        return out, ladder_resid

    broken = {e[:2]: broken_object_violates(e) for e in ((1, 3, 1.0), (2, 4, 1.0))}
    check("(h) N=6, at three generic couplings and at the six RUNG-2 couplings (the rung-3 ones "
          "of E9 are a different check's and are not in reach here): the diagonal "
          "block (p,p) carries at least as many frozen modes as BOTH side lines at the same "
          "p + q, which is mult(m=0) >= mult(m=1) read off the object, the inequality open "
          "item 4 of the note rests on and the one (e) deliberately does not assert. TWO "
          "controls, both required. The OBJECT control: add the bond (1,3) or (2,4). The "
          "ladder commutator on that object is COMPUTED here and is nonzero, so ker(L - lambda) "
          "is not an su(2) module and the theorem is silent; the inequality then does break, at "
          "p = 1 and p = 5. It is not that any broken object must break it, and the comment "
          "above names four added bonds that leave it standing: these two were FOUND, not "
          "predicted. The PAIRING control: pair the diagonal with the next DIAGONAL block, a "
          "different p + q, and the same inequality must hold generically and FAIL at the "
          "exceptional couplings",
          not viol and not wrong_generic and len(wrong) > 0
          and all(len(v) > 0 and r > 0 for v, r in broken.values()),
          f"{seen_pairs} weight triples, violations {viol}; object control "
          f"{'; '.join(f'{k}: {v} at ladder residual {r:.3g}' for k, (v, r) in broken.items())}; "
          f"wrong pairing: {len(wrong_generic)} violations generic (must be 0), {len(wrong)} "
          f"exceptional (must be > 0), e.g. {wrong[:2]}")


# ------------------------------------------------------------------- E7

def e7_bridge():
    print("E7  the one coupling the bridge needs")
    for N in ([4, 5, 6] + ([7] if DEEP else [])):
        m, worst = N // 2, None
        ok = True
        for (p, q) in band(N):
            C, A0 = block_parts(N, p, q)
            if C.shape[0] > (5000 if DEEP else 1300):
                continue
            nul = modp_nullity(C, A0, 1)
            if nul != m:
                ok, worst = False, (p, q, nul)
        check(f"(a) N={N}: an exact GF(p) rank at J = 1 gives nullity floor(N/2) on every "
              f"band block it reaches", ok, f"{worst if worst else 'all ' + str(m)}")
    print("  [--] a GF(p) rank bounds the count from ABOVE over QQ, so one coupling fixes "
          "the generic value, and the floor closes it to an equality. That step is the "
          "pencil argument of Proposition 5.3 and is not measured here.")


# ------------------------------------------------------------------- E8

def fixed_kernel(C, A0):
    """ker C n ker A0, exactly.  The frozen modes the floor supplies are annihilated by BOTH
    halves separately, so this rational subspace is the J-independent part of every kernel."""
    ns = Matrix(np.vstack([C, A0]).astype(int).tolist()).nullspace()
    return Matrix.hstack(*ns) if ns else Matrix.zeros(C.shape[0], 0)


def deflate(C, A0, K):
    """A rational congruence removing the fixed kernel.  M(z) = C + z*A0 is complex SYMMETRIC,
    so putting a basis of K first makes the first dim K rows AND columns of T^T M T vanish,
    and rank M = rank of what is left, at every z."""
    d, m = C.shape[0], K.cols
    _, piv = K.T.rref()
    keep = [i for i in range(d) if i not in set(piv)]
    T = Matrix.hstack(K, Matrix.eye(d)[:, keep])
    CT = T.T * Matrix(C.astype(int).tolist()) * T
    AT = T.T * Matrix(A0.astype(int).tolist()) * T
    clean = all(CT[i, j] == 0 and AT[i, j] == 0 and CT[j, i] == 0 and AT[j, i] == 0
                for i in range(m) for j in range(d))
    return CT[m:, m:], AT[m:, m:], clean, T.det() != 0


def crossings_of_module(N, p, q, loose=1e-6):
    """The couplings at which the FROZEN value is an eigenvalue of this block, numerically.

    On a band block this is empty by construction: the coupling-independent part of the frozen
    space is ker C n ker A0, which the deflation removes.  Off the band it is the set of isolated
    crossings, and E9(g) is what says those are not freezings.

    It UNDER-counts, deliberately and by a named amount.  A multiple root leaves ordqz with an
    imaginary ratio around 1e-5, above the candidacy cutoff, so the triple root at 2/sqrt3 on
    (2,3) at N = 5 is dropped and this returns five where E8(g) proves six.  Nothing here rests
    on the count: E9(e) needs the list to be NON-EMPTY, which is a verdict a cutoff cannot move.
    The j > 0.05 filter is the z = 0 smear of the second trap in the note, not a tuning knob."""
    C, A0 = block_parts(N, p, q)
    Cg, Ag, _ = deflate_orthogonal(C.astype(float), A0.astype(float))
    if Cg.shape[0] == 0:
        return []
    _, _, al, be, _, _ = sla.ordqz(Cg, -Ag, output="complex")
    fin = np.abs(be) > 1e-10 * max(np.abs(be).max(), 1e-300)
    z2 = (al[fin] / be[fin]) ** 2
    ok = (np.abs(z2.imag) / np.maximum(np.abs(z2), 1e-300) < loose) & (z2.real < -1e-12)
    out = []
    for j in np.sort(np.sqrt(-z2.real[ok])):
        if j > 0.05 and not (out and j - out[-1] <= 1e-7 * max(j, 1.0)):
            out.append(float(j))
    return out


def full_block_roots(N, p, q, J):
    """The exceptional couplings of the WHOLE block, exactly: the real roots of the deflated
    determinant, with no subspace assumed anywhere."""
    z = symbols('z')
    C, A0 = block_parts(N, p, q)
    K = fixed_kernel(C, A0)
    CG, AG, clean, invertible = deflate(C, A0, K)
    n = CG.rows
    pts = [(Rational(k - n // 2),
            DomainMatrix.from_Matrix(CG + Rational(k - n // 2) * AG).convert_to(QQ).det())
           for k in range(n + 1)]
    qz = Poly(interpolate(pts, z), z)
    co = list(reversed(qz.all_coeffs()))
    even = all(c == 0 for k, c in enumerate(co) if k % 2)
    R = Poly(sum(c * (-1) ** (k // 2) * J ** k for k, c in enumerate(co) if k % 2 == 0), J)

    def pos(f):
        pf = Poly(f, J)
        return pf.count_roots(0, oo) - (1 if pf.eval(0) == 0 else 0)

    simple = all(mult == 1 or pos(f) == 0 for f, mult in factor_list(R.as_expr())[1])
    roots = sorted(float(r) for r in R.real_roots() if r.is_real and float(r) > 1e-9)
    return dict(m=K.cols, clean=clean and invertible, even=even, count=pos(R.as_expr()),
                simple=simple, roots=roots, poly=R)


def e8_full_block():
    print("E8  the WHOLE block, exactly, at every coupling at once")
    N, m = 5, 2
    J = symbols('J')
    res = {b: full_block_roots(N, *b, J) for b in band(N)}
    for b, r in res.items():
        EXACT_BLOCK[(N, *b)] = r["count"]
    ok = all(r["clean"] and r["m"] == m and r["even"] for r in res.values())
    check("(a) N=5: ker C n ker A0 is rational of dimension floor(N/2) on every band block, "
          "and the congruence that removes it is clean", ok)
    raised = sorted(b for b, r in res.items() if r["count"])
    counts = {b: res[b]["count"] for b in raised}
    check("(b) N=5: over ALL couplings, only the two diagonal blocks carry exceptional "
          "couplings, and exactly three each",
          raised == [(2, 2), (3, 3)] and set(counts.values()) == {3}, f"{counts}")
    # the comparison against E5's singlet polynomial is by DIVISIBILITY, so it is exact:
    # every irreducible factor of the singlet polynomial carrying a positive root must divide
    # the whole block's polynomial, and the counts must agree.
    g, _, _, _, _ = exceptional_polynomial(N, 2, J)

    def pos_factors(poly):
        out = []
        for f, _mult in factor_list(Poly(poly, J).as_expr())[1]:
            pf = Poly(f, J)
            if pf.count_roots(0, oo) - (1 if pf.eval(0) == 0 else 0) > 0:
                out.append(pf)
        return out

    full = Poly(res[(2, 2)]["poly"], J)
    divides = all(full.rem(f).is_zero for f in pos_factors(g))
    same = divides and res[(2, 2)]["count"] == 3 and res[(3, 3)]["count"] == 3
    check("(c) N=5: every factor of E5's SINGLET polynomial that carries a coupling divides "
          "the whole block's polynomial, and the counts agree, so the singlets are not merely "
          "among the exceptional modes but all of them", same,
          f"{[round(x, 9) for x in res[(2, 2)]['roots']]}")
    check("(d) N=5: every positive root is simple, so at an exceptional coupling the count is "
          "EXACTLY floor(N/2)+1 and not merely at least that",
          all(r["simple"] for r in res.values()))
    check("(e) N=5: no side line carries an exceptional coupling at any J, which the sampled "
          "couplings of E3 could not have decided",
          all(res[b]["count"] == 0 for b in band(N) if b[0] != b[1]))
    if DEEP:
        r = full_block_roots(6, 2, 2, J)
        EXACT_BLOCK[(6, 2, 2)] = r["count"]
        check("(f) N=6: the whole block (2,2) carries exactly six exceptional couplings, all "
              "simple, which is again the singlet count",
              r["count"] == 6 and r["simple"] and r["clean"], f"{[round(x, 6) for x in r['roots']]}")
        # (g) the OFF-BAND count, exactly, because E9(g) makes a statement about these roots and
        # a number carried in a label is not a measurement.  It also pins the gap between the
        # exact route and the numeric one: the pencil detector of E9(e) finds FEWER, because the
        # triple root at 2/sqrt3 leaves ordqz with an imaginary ratio around 1e-5 and the
        # candidacy cutoff there is 1e-6.  That gap is why E9(e) needs only non-emptiness.
        off = {b: full_block_roots(5, *b, J) for b in ((2, 3), (1, 2))}
        found = crossings_of_module(5, 2, 3)
        numeric = len(found)
        # The root the numeric route drops is pinned EXACTLY, not to a tolerance: it is the
        # factor 3J^2 - 4 of the block's own polynomial, and its multiplicity is read off the
        # factorisation rather than asserted in prose.  A float comparison here would have been
        # a gate on a residual, where the comparison installed below has residual exactly 0.0
        # because both floats come from the same CRootOf.  (A float compare against the
        # LITERAL 2/sqrt(3) would have been off by one ulp, 2.2e-16, and would have needed a
        # threshold; the literal is not what is used.)
        # Watch the two senses of "six": full_block_roots COUNTS distinct roots by Sturm and
        # LISTS them with multiplicity, so `count` is 6 while `roots` has 8 entries.  The first
        # version of this check asserted one dropped entry and got three, all equal, which is
        # the better statement anyway: the number dropped IS the multiplicity of the factor.
        fl = {Poly(f, J).as_expr(): m for f, m in factor_list(off[(2, 3)]["poly"].as_expr())[1]}
        triple = fl.get(3 * J ** 2 - 4)
        missing = [r for r in off[(2, 3)]["roots"] if all(abs(r - f) > 1e-6 for f in found)]
        # TIE the dropped value to the factor rather than printing the identification beside it:
        # both floats descend from the same exact algebraic number, 2*sqrt(3)/3, so they are
        # bit-identical and == is the right comparison.  Without this the check would pass on
        # any multiplicity-3 factor that happened to sit there.
        # The first two terms of `tied` are CRASH GUARDS and not measurements: a quadratic
        # literal always has one positive root, and `missing` is already pinned non-empty by
        # len(missing) == triple below.  They are here so that max()/[0] cannot raise, because
        # a check that raises aborts the run instead of printing a FAIL, which is worse than a
        # check that cannot fail.  The measurement is the third term.
        tri_root = [float(r) for r in Poly(3 * J ** 2 - 4, J).real_roots() if float(r) > 0]
        tied = len(tri_root) == 1 and bool(missing) and max(missing) == tri_root[0] == min(missing)
        check("(g) N=5, OFF the band: the whole blocks (2,3) and (1,2) carry six and four "
              "DISTINCT positive roots, over all couplings at once. They are real and the "
              "eigenvalue there is -4*gbar; E9(g) is what says they are crossings and not "
              "freezings. The numeric pencil route finds fewer, and the difference is a "
              "multiple root and not a disagreement: what it drops is the whole multiplicity-3 "
              "factor 3J^2 - 4, the number of dropped entries IS that multiplicity, and the "
              "value dropped IS that factor's positive root, compared exactly because both "
              "floats descend from the same algebraic number. The one part NOT exact is the "
              "numeric count itself, and it is pinned as an equality on purpose: if the pencil "
              "ever resolves the triple root this must fail loudly and be read, not pass",
              off[(2, 3)]["count"] == 6 and off[(1, 2)]["count"] == 4
              and off[(2, 3)]["clean"] and off[(1, 2)]["clean"] and numeric == 5
              and triple == 3 and len(missing) == triple and tied,
              f"(2,3) exact {off[(2, 3)]['count']} distinct against numeric {numeric}; "
              f"{len(missing)} dropped entries, all equal: {tied}, value "
              f"{missing[0] if missing else 'NONE'}, factor multiplicity {triple}; (1,2) exact "
              f"{off[(1, 2)]['count']}")


# ------------------------------------------------------------------- E9

def e9_higher_rungs():
    print("E9  a higher rung carries couplings of its own: rung 3 at N = 6")
    J = symbols('J')
    g3, dimV3, q3, ok3, even3 = exceptional_polynomial(6, 3, J)
    n3 = g3.count_roots(0, oo) - (1 if g3.eval(0) == 0 else 0)
    EXACT_SINGLET[(6, 3)] = n3
    roots3 = sorted(float(r) for r in g3.real_roots() if r.is_real and float(r) > 1e-9)
    check("(a) N=6: the eta-lowest-weight singlets of the block (3,3) carry exactly four "
          "exceptional couplings, by Sturm on an exact polynomial", n3 == 4 and ok3 and even3,
          f"dim V = {dimV3}, deg q = {q3.degree()}, roots {[round(r, 9) for r in roots3]}")
    g2, _, _, _, _ = exceptional_polynomial(6, 2, J)
    # the two polynomials necessarily share the root J = 0, which is not a coupling, so the
    # exact statement is that their gcd carries no POSITIVE root.
    common = Poly(gcd(g2, g3), J)
    shared = common.count_roots(0, oo) - (1 if common.eval(0) == 0 else 0)
    check("(b) and they are NEW: the gcd of the rung-3 and rung-2 polynomials carries no "
          "positive root, so none of the four is one of E5's six", shared == 0,
          f"gcd of degree {common.degree()}, {shared} positive roots")
    # a rung-3 lowest weight at N = 6 has eta-spin N/2 - 3 = 0, so its multiplet is a single
    # block: it must raise (3,3) and nothing else, where E6's rung-2 mode raises 2 .. N-2.
    shape_ok, seen = True, {}
    for J0 in roots3:
        counts = {}
        for (p, q) in band(6):
            C, A0 = block_parts(6, p, q)
            s = np.linalg.svd(C + 1j * J0 * A0, compute_uv=False)
            counts[(p, q)] = int(np.sum(s < 1e-9 * max(s[0], 1.0)))
        raised = sorted(k for k, v in counts.items() if v > 3)
        seen[round(J0, 9)] = raised
        shape_ok = shape_ok and raised == [(3, 3)] and counts[(3, 3)] == 4
    check("(c) at each of the four the nullity rises on (3,3) and on NO other block, which is "
          "the multiplet of a rung-3 lowest weight at N = 6 (eta-spin 0) and cannot be a "
          "rung-2 one, those occupying 2 .. N-2 (a numeric read of the shape)", shape_ok,
          f"{seen}")
    g2roots = sorted(float(r) for r in g2.real_roots() if r.is_real and float(r) > 1e-9)
    also33 = True
    for J0 in g2roots:
        C, A0 = block_parts(6, 3, 3)
        sv = np.linalg.svd(C + 1j * J0 * A0, compute_uv=False)
        also33 = also33 and int(np.sum(sv < 1e-9 * max(sv[0], 1.0))) > 3
    check("(d) the six rung-2 couplings raise the block (3,3) as well, so that block carries "
          "the six AND the four and the exceptional set at N = 6 has at least ten members: "
          "the count on (2,2) is the rung-2 count and not the count",
          also33 and n3 == 4 and shared == 0, f"six rung-2 roots raise (3,3): {also33}")

    # (f) and (g): OUTSIDE the band.  Every side-line statement this file makes, here and in
    # E6 and E8(e), loops over band(N), so until 2026-08-08 no block with |p - q| > 2 and no
    # block with p - q ODD had been looked at, at any N.  That is exactly the case the m = 0
    # argument of the note's open item 4 cannot reach: a half-integer weight lives there and
    # its multiplet has no zero-weight member, so an exceptional coupling hiding there would be
    # invisible to every diagonal block including the middle one.
    def offband_scan(N, couplings, generic):
        """Which off-band blocks gain a mode at these couplings, over the generic baseline."""
        base, hits = {}, []
        for p in range(N + 1):
            for q in range(N + 1):
                if abs(p - q) in (0, 2):
                    continue
                C, A0 = block_parts(N, p, q)
                s = np.linalg.svd(C + 1j * generic * A0, compute_uv=False)
                base[(p, q)] = int(np.sum(s < 1e-9 * max(s[0], 1.0)))
        for J0 in couplings:
            for b in base:
                C, A0 = block_parts(N, *b)
                s = np.linalg.svd(C + 1j * J0 * A0, compute_uv=False)
                if int(np.sum(s < 1e-9 * max(s[0], 1.0))) > base[b]:
                    hits.append((round(J0, 9), b))
        return hits, base

    # 1.3172583 is a control coupling and nothing more: it is not near any known root at N = 5
    # or N = 6, band or OFF-band: the nearest is 1.372601, a crossing of (2,3), 0.0553 away. The
    # off-band roots count here because this baseline is taken over off-band blocks; the first
    # version of this comment quoted 0.0885, the distance to the nearest BAND root, before
    # E8(g) made the off-band ones known. If a future N puts a root beside it, the baseline
    # inflates silently, so any extension of this check re-picks it.
    generic = 1.3172583
    ten = sorted(roots3 + g2roots)
    n5roots = sorted(float(r) for r in EXACT_POLY[(5, 2)].real_roots()
                     if r.is_real and float(r) > 1e-9)
    offband_crossings = crossings_of_module(5, 2, 3)
    hits6, base6 = offband_scan(6, ten, generic)
    hits5, base5 = offband_scan(5, n5roots, generic)
    # The POSITIVE CONTROL, and it has to be built from the same detector on the same blocks,
    # not from the side lines: those are band blocks and say nothing about whether the detector
    # fires off-band.  At N = 5 the off-band blocks DO carry isolated crossings of the frozen
    # value, at couplings of their own (see (g)); feeding those to the same scan must produce
    # hits.  Without this the whole check passes on random couplings, which is what the first
    # version at N = 6 did: no off-band block there crosses at ANY J, so there was nothing to
    # see and ten random numbers scored the same as the ten exceptional ones.
    ctrl_hits, _ = offband_scan(5, offband_crossings, generic)
    check("(e) POSITIVE CONTROL for (f), on the SAME blocks with the SAME detector: at N = 5 "
          "the off-band blocks do cross the frozen value, at couplings of their own, and the "
          "scan finds them. Without this (f) passes on any list of numbers, which is what its "
          "first version did at N = 6, where no off-band block crosses at ANY J",
          len(ctrl_hits) > 0,
          f"{len(ctrl_hits)} off-band hits at the {len(offband_crossings)} crossing couplings "
          f"of (2,3), e.g. {ctrl_hits[:3]}")
    check("(f) N=5 and N=6: at NONE of the thirteen exceptional couplings does a block with "
          "p - q ODD or with |p - q| > 2 gain a mode, so the extra mode stays inside the band "
          "and not merely inside the blocks that were scanned. The N=6 half measures nothing "
          "on its own, since no off-band block there crosses at any J; the N=5 half is the one "
          "(e) certifies",
          not hits5 and not hits6,
          f"N=5 {hits5} over {len(base5)} blocks at {len(n5roots)} couplings, N=6 {hits6} over "
          f"{len(base6)} blocks at {len(ten)} couplings")

    # (g) the trap, recorded because walking into it costs an afternoon and a false alarm.  Run
    # the EXACT route of E8 on an OFF-BAND block and it returns positive roots: at N = 5 the
    # block (2,3) yields six, (1,2) and (3,4) four each, against a complete count of three for
    # the whole chain.  They are real, and the eigenvalue there really is -4*gbar.  They are
    # also content-free, and the discriminator is exact and needs no threshold at all: on a
    # band block the frozen space is ker C n ker A0, which does not depend on the coupling and
    # is there at EVERY J, while off the band that intersection is ZERO and the roots are
    # isolated crossings of a value that is nothing in particular for that block to pass
    # through.  A CROSSING IS NOT A FREEZING.  Integer matrices, exact GF(p) rank, and the rank
    # over GF(p) can only be smaller than over QQ, so a zero read here IS zero over QQ.
    def common_kernel_dim(N, p, q):
        """dim(ker C n ker A0), exactly over QQ.

        Called on 143 blocks, the 6 corners (0,0) and (N,N) falling through both branches
        below; about thirteen seconds in all.

        The first version read this over GF(p), where rank_p <= rank_QQ and so nullity_p >=
        nullity_QQ.  That direction is sound for the OFF-BAND half, where a measured 0 forces 0
        over QQ, and unsound for the BAND half, where it would give only "at most floor(N/2)"
        while the label asserts equality.  Rather than state the two halves differently, the
        whole thing is exact."""
        C, A0 = block_parts(N, p, q)
        S = np.vstack([C, A0]).astype(int)
        dm = DomainMatrix.from_Matrix(Matrix(S.tolist())).convert_to(QQ)
        return S.shape[1] - len(dm.rref()[1])

    inband, offb, skipped = {}, {}, []
    for Nk in (5, 6, 7):
        for p in range(Nk + 1):
            for q in range(Nk + 1):
                if len(subsets(Nk, p)) * len(subsets(Nk, q)) > 1300:
                    skipped.append((Nk, p, q))       # never fires below N = 8; announced, not silent
                    continue
                if abs(p - q) in (0, 2) and (p, q) not in ((0, 0), (Nk, Nk)):
                    inband[(Nk, p, q)] = common_kernel_dim(Nk, p, q)
                elif abs(p - q) not in (0, 2):
                    offb[(Nk, p, q)] = common_kernel_dim(Nk, p, q)
    band_ok = all(v == k[0] // 2 for k, v in inband.items())
    off_ok = all(v == 0 for v in offb.values())
    check("(g) the trap, and its exact resolution: the exact route of E8 run on an OFF-BAND "
          "block returns positive roots, real, and the eigenvalue there really is -4*gbar (the "
          "count is E8(g) under --deep), yet they are ordinary level crossings and not "
          "freezings. The discriminator carries no threshold and is exact over QQ: a band "
          "block's frozen space contains a COUPLING-INDEPENDENT part, ker C n ker A0, of "
          "dimension floor(N/2), and off the band that intersection is ZERO",
          band_ok and off_ok and not skipped,
          f"N=5,6,7: {len(inband)} band blocks all at floor(N/2), {len(offb)} off-band blocks "
          f"all at 0; worst off-band {max(offb.values())}; blocks skipped by the size cap "
          f"{skipped}")
    if DEEP:
        # at N = 7 a rung-3 lowest weight has eta-spin N/2 - 3 = 1/2, so its multiplet is TWO
        # blocks, 3 .. N-3.  This coupling is the smallest the (3,3) enumeration offers and is
        # absent from the rung-2 list, so it must raise (3,3) and (4,4) and neither (2,2) nor
        # (5,5).  A dip ratio, so nine digits of the root are enough and no tolerance is fixed.
        J0, width = 0.473255094, {}
        for (p, q) in [(2, 2), (3, 3), (4, 4), (5, 5)]:
            C, A0 = block_parts(7, p, q)

            def sm(x, C=C, A0=A0):
                sv = np.sort(np.linalg.svd(C + 1j * x * A0, compute_uv=False))
                return sv[3] / sv[-1]

            width[(p, q)] = sm(J0) / max(min(sm(J0 * 0.999), sm(J0 * 1.001)), 1e-300)
        check("(h) N=7: the smallest coupling of the block (3,3) collapses the fourth "
              "singular value on (3,3) and (4,4) and on neither (2,2) nor (5,5), the two-block "
              "multiplet of a rung-3 lowest weight at eta-spin 1/2 (a numeric read)",
              width[(3, 3)] < 1e-3 and width[(4, 4)] < 1e-3
              and width[(2, 2)] > 0.1 and width[(5, 5)] > 0.1,
              f"dip ratios {[(b, f'{w:.1e}') for b, w in width.items()]}")


def e10_slow_exact():
    print("E10  the two exact counts the other modes are too short for")
    J = symbols('J')
    for (tag, N, p, expect, small) in [("a", 7, 3, 13, 0.473255094), ("b", 8, 2, 15, 0.790268421)]:
        g, dimV, q, ok, even = exceptional_polynomial(N, p, J)
        npos = g.count_roots(0, oo) - (1 if g.eval(0) == 0 else 0)
        EXACT_SINGLET[(N, p)] = npos
        EXACT_POLY[(N, p)] = g
        roots = sorted(float(r) for r in g.real_roots() if r.is_real and float(r) > 1e-9)
        # the smallest is TABULATED, so it is asserted and not merely printed
        pinned = bool(roots) and abs(roots[0] - small) < 1e-9
        check(f"({tag}) N={N}, the rung-{p} singlets of the block ({p},{p}): the exact Sturm "
              f"count of the exceptional couplings, and the smallest of them is the value the "
              f"note tabulates", npos == expect and ok and even and pinned,
              f"dim V = {dimV}, deg q = {q.degree()}, {npos} roots, smallest "
              f"{roots[0]:.9f}" if roots else f"{npos} roots")


# ------------------------------------------------------------------- E11

def deflate_orthogonal(C, A0):
    """An orthonormal restriction to (ker C n ker A0)^perp.

    C + z A0 is real SYMMETRIC, so the common kernel is annihilated on both sides and the rank
    is the rank of the restriction at every z.  Orthonormal, so the deflation itself spends no
    conditioning.  This is E8's `deflate` with the exact congruence replaced by a numerical
    one, which is what buys the blocks the interpolation cannot reach."""
    S = np.vstack([C, A0]).astype(float)
    _, s, vh = sla.svd(S, full_matrices=True)
    tol = max(S.shape) * np.finfo(float).eps * (s[0] if s.size else 1.0)
    s_full = np.zeros(S.shape[1])
    s_full[: s.size] = s
    Q = vh[s_full > tol].conj().T
    return Q.T @ C @ Q, Q.T @ A0 @ Q, C.shape[0] - Q.shape[1]


def det_sign(Cg, Ag, x):
    """The sign of det(C_G + i x A_G), real because the polynomial is even in z (E1(d) for the
    block; the orthogonal deflation below preserves it, since the chiral involution commutes
    with C and anticommutes with A0 and so restricts to the complement of their common kernel).
    Carried as a phase, so the modulus may overflow freely.

    The second value is |sin(phase)|, the LINEAR distance from the real axis, and the sign
    returned is meaningful exactly while it is small.  |cos| - 1 would be the natural-looking
    reading and is the wrong one: cos is flat at its extrema, so it reports a phase error d as
    d^2/2 and returns exactly 0.0 for every d below sqrt(eps), about 1.5e-8 rad.  That is cos's
    own blindness and not an error estimate; the drift actually accumulated by summing n angles
    is about sqrt(n)*pi*eps, five decades smaller, so the old metric returned 0.0 on every run
    it was ever asked.  A quantity that reads zero because it cannot see is not a
    measurement."""
    lu, piv = sla.lu_factor(Cg + 1j * x * Ag)
    d = np.diag(lu)
    phase = float(np.sum(np.angle(d))) + np.pi * (int(np.sum(piv != np.arange(len(piv)))) % 2)
    return (1.0 if np.cos(phase) > 0 else -1.0), abs(float(np.sin(phase)))


def numeric_block_roots(N, p, q, loose=1e-3, smallj=0.05, eps=1e-5):
    """The exceptional couplings of the WHOLE block (p,q), numerically but each one certified.

    Detection: the finite spectrum of the deflated pencil.  z = iJ is a root of an even real
    polynomial, so a purely imaginary z is a real NEGATIVE z^2, a well-conditioned question
    where "is the real part small" is not.  Candidacy is deliberately LOOSE, since a cluster of
    close roots at dimension 4900 pushes a genuine z^2 off the real axis; a loose net costs two
    LU factorisations per false candidate and a tight one loses roots silently.

    Certification: a SIGN BRACKET of the determinant across J(1 +- eps).  A sign change over a
    window narrower than the gap to the next root proves an odd-order root inside it, so every
    ACCEPTED root is real.  Nothing here certifies that no root was missed, so the returned
    count is a LOWER bound -- but only if nothing is counted TWICE, and exactly ONE premise
    stands between this routine and that word:

      * DISJOINTNESS, and it is SUFFICIENT, which is all the count needs.  Two windows that
        both flip can only be one root counted twice if they both contain it; and if they both
        contain rho then j2(1-eps) <= rho <= j1(1+eps), so j2 - j1 <= eps*(j1 + j2).  Disjoint
        windows therefore contain distinct roots and the count is a lower bound.  It is NOT
        equivalent: two overlapping windows each holding a root of its own outside the overlap
        also count correctly, so overlap would only mean the premise is no longer measured, not
        that the bound has failed.  `window_ratio` -- the smallest gap between consecutive
        accepted roots over the distance eps*(j_i + j_i+1) the two windows reach toward each
        other -- is that premise as a number: above 1 it holds, and the value says by how much.

    `mults`, how many pencil eigenvalues each accepted root merged, is kept beside it as a
    detector diagnostic and NOT as a second premise.  A root arrives as the pair +-z, so the
    healthy value is 2; a 1 means one copy was lost to the candidacy filter, which is
    UNDER-detection and cannot threaten a lower bound, and a split pair that survives as two
    candidates is already caught by disjointness above.  It is worth printing because the split
    happens: the z = 0 smear at N = 8 (4,4) arrives as two unmerged copies, and the run shows
    them as merged 1 where every accepted root shows merged 2.

    Rejected on purpose: J below `smallj`.  C is diagonal and vanishes wherever the
    disagreement is 2, so z = 0 is a root of high multiplicity whose numerical image smears to
    about 1e-4 (Jordan blocks of size four, the fourth root of machine epsilon).  Counting the
    smear is how a law that is not there gets restored.  The threshold sits two decades above
    the smear it removes, so `smear_ratio` reports the true separation actually measured: the
    smallest accepted root over the largest rejected one, which is the number that says whether
    the cut had room, where 0.05 on its own says nothing."""
    C, A0 = block_parts(N, p, q)
    Cg, Ag, ker = deflate_orthogonal(C.astype(float), A0.astype(float))
    _, _, alpha, beta, _, _ = sla.ordqz(Cg, -Ag, output="complex")
    fin = np.abs(beta) > 1e-10 * max(np.abs(beta).max(), 1e-300)
    z2 = (alpha[fin] / beta[fin]) ** 2
    rel = np.abs(z2.imag) / np.maximum(np.abs(z2), 1e-300)
    ok = (rel < loose) & (z2.real < -1e-12)
    merged = []                                  # every root gives +-z, so z^2 appears twice
    for j, r in sorted(zip(np.sqrt(-z2.real[ok]), rel[ok])):
        if merged and j - merged[-1][0] <= 1e-7 * max(j, 1.0):
            merged[-1][1] += 1
            merged[-1][2] = min(merged[-1][2], float(r))
        else:
            merged.append([float(j), 1, float(r)])
    kept, mults, rejected, stray = [], [], [], 0.0
    for j, mult, r in merged:
        if j < smallj:
            rejected.append((j, mult, "inside the z = 0 smear"))
            continue
        a, sa = det_sign(Cg, Ag, j * (1 - eps))
        b, sb = det_sign(Cg, Ag, j * (1 + eps))
        stray = max(stray, sa, sb)
        if a != b:
            kept.append(j)
            mults.append(mult)
        else:
            rejected.append((j, mult, "the sign bracket did not flip"))
    # the premise, as a measured margin.  Each window reaches eps*J toward its neighbour, so
    # what a gap has to clear is eps*(J_i + J_i+1), the sum of the two REACHES and half the sum
    # of the two widths.  The wording matters: 2.8 here is 1.4 combined window widths.
    ratios = [(kept[i + 1] - kept[i]) / (eps * (kept[i] + kept[i + 1]))
              for i in range(len(kept) - 1)]
    # only the SMEAR rejections belong in the smallj margin.  A candidate the bracket refused
    # at an ordinary coupling says nothing about whether the cut near zero had room, and mixing
    # the two would make the margin fail for a reason it does not name.
    smeared = [j for j, _, why in rejected if why == "inside the z = 0 smear"]
    # the phase drift has an error model: a sum of n angles, each rounded, drifts by about
    # sqrt(n)*pi*eps.  What is worth gating is the RATIO to that, held across the decades of n
    # this block spans, not a fixed threshold that eleven decades of headroom would pass.
    model = (Cg.shape[0] ** 0.5) * np.pi * np.finfo(float).eps
    return dict(roots=kept, mults=mults, rejected=rejected, ker=ker, dim=C.shape[0], stray=stray,
                parity=sorted(set(mults)), n=Cg.shape[0], stray_over_model=stray / model,
                window_ratio=min(ratios) if len(kept) > 1 else None,
                smear_ratio=(min(kept) / max(smeared)) if kept and smeared else None)


def nested(small, large, tol=2e-5):
    """Every root of the smaller block matched to one of the larger, ONE TO ONE.

    The matching has to be injective or the nesting proves less than it looks: two roots of the
    small list both matching one root of the large list would pass a per-root test while the
    set difference is one larger than the count difference.  (Greedy first-fit, so it can still
    report a miss where a valid injection exists; that error is a spurious FAILURE, which is
    the safe direction.)

    The tolerance is TWICE the bracket's reach, and RELATIVE throughout.  Two blocks place the
    same root independently, each within eps = 1e-5 relative of it, so their two values can
    differ by 2e-5 relative and a tolerance of 1e-5 would report a false miss on a genuinely
    shared root.  It is still far from merging distinct roots: the smallest spacing measured at
    N = 8 is 5.6e-5 relative, a factor 2.8.  Relative and not max(j, 1) on purpose, since the
    roots crowd at the SMALL-J end, exactly where a floor of 1 would quietly turn the tolerance
    absolute and eat most of that factor."""
    taken, miss = set(), []
    for j in small:
        hit = next((i for i, k in enumerate(large)
                    if i not in taken and abs(j - k) <= tol * j), None)
        if hit is None:
            miss.append(j)
        else:
            taken.add(hit)
    return not miss, miss


def e11_numeric_rungs():
    print("E11  past the exact route: the middle block at N = 8, numerically")
    J = symbols('J')
    # (a) the detector against this file's own exact roots, digit for digit and not by count.
    #     N = 6, block (3,3) carries the six of rung 2 (E5) and the four of rung 3 (E9), and
    #     both lists are exact here, so a numeric count of ten is not enough: the ten VALUES
    #     must be the ten exact ones.
    exact6 = sorted([float(r) for g in (exceptional_polynomial(6, 2, J)[0],
                                        exceptional_polynomial(6, 3, J)[0])
                     for r in g.real_roots() if r.is_real and float(r) > 1e-9])
    r6 = numeric_block_roots(6, 3, 3)
    want6 = EXACT_SINGLET[(6, 2)] + EXACT_SINGLET[(6, 3)]      # the six and the four, from E5/E9
    same = (len(r6["roots"]) == len(exact6)
            and all(abs(a - b) < 1e-9 * max(a, 1.0) for a, b in zip(r6["roots"], exact6)))
    check("(a) N=6 (3,3): the pencil detector returns exactly the roots the exact polynomials "
          "of E5 and E9 give, value for value and not merely as many of them",
          same and len(exact6) == want6,
          f"{len(r6['roots'])} numeric vs {len(exact6)} exact ({want6} expected from the two "
          f"Sturm counts), worst relative difference "
          f"{max((abs(a - b) / max(a, 1.0) for a, b in zip(r6['roots'], exact6)), default=0):.1e}"
          f" against the 1e-9 the check requires, {len(r6['rejected'])} rejected")
    # (b) the counts the exact route reaches, by the numeric route.  Every expected value is
    #     READ OUT of the exact blocks of this same run and none is typed here, so a stale
    #     literal cannot hide in the validation.  The N = 8 entry is the one that matters: the
    #     detector is validated AT the N it is then trusted at, on the block below the ones
    #     being measured, against E10(b)'s Sturm count.
    #     Where the WHOLE block was enumerated over QQ (N = 5, and N = 6 under deep) the
    #     comparison is like for like.  From N = 7 on only the singlet space is exact, so the
    #     comparison additionally assumes the singlets are all of it, which is the note's own
    #     open identification and is said here rather than buried.
    table = [(5, 2, 2, lambda: EXACT_BLOCK[(5, 2, 2)], "E8, the whole block over QQ"),
             (6, 2, 2, lambda: EXACT_BLOCK[(6, 2, 2)], "E8(f), the whole block over QQ"),
             (7, 2, 2, lambda: EXACT_SINGLET[(7, 2)], "E5, the singlet space"),
             (7, 3, 3, lambda: EXACT_SINGLET[(7, 2)] + EXACT_SINGLET[(7, 3)],
              "E5 + E10(a), rung 2 plus rung 3 on the singlet spaces"),
             (8, 2, 2, lambda: EXACT_SINGLET[(8, 2)], "E10(b), the singlet space")]
    valblocks = []
    for (N, p, q, want_fn, where) in table:
        try:
            want = want_fn()
        except KeyError as miss:
            check(f"(b) N={N} ({p},{q}): the exact count this validation needs was computed in "
                  f"this run", False, f"missing {miss}")
            continue
        r = numeric_block_roots(N, p, q)
        valblocks.append(r)
        # equality, and both ways of failing it are worth having.  FEWER means the detector
        # missed a root it should have found.  MORE, from N = 7 on where the exact number is
        # the SINGLET count, would mean an exceptional coupling outside the singlet spaces,
        # which is this note's open identification resolving in the interesting direction and
        # is a finding to write up rather than a bug to fix.
        check(f"(b) N={N} ({p},{q}): the numeric count is the exact one, {want} from {where}",
              len(r["roots"]) == want, f"{len(r['roots'])}, ker {r['ker']} of {r['dim']}, "
              f"worst |sin(phase)| {r['stray']:.1e}, {len(r['rejected'])} rejected, "
              f"window margin {r['window_ratio']:.1f}x")
    # (b2) the METHOD, not the detector: wherever an exact rung count exists, the difference of
    #      two numeric block counts must reproduce it.  This is the only place the subtraction
    #      itself is tested, and it tests both premises it rests on at once, completeness of the
    #      smaller list and disjointness of the two rungs.  It exists at N = 6 and N = 7 and can
    #      exist for no rung 4 anywhere: at N = 8 the rung-4 multiplet has eta-spin 0 and lives
    #      in (4,4) alone, so no larger block can read the twelve back.  A Sturm count on the
    #      rung-4 singlet space of (4,4) would be a second route, and has not been run.
    for (N, want_fn, where) in [(6, lambda: EXACT_SINGLET[(6, 3)], "E9(a)"),
                                (7, lambda: EXACT_SINGLET[(7, 3)], "E10(a)")]:
        got = len(numeric_block_roots(N, 3, 3)["roots"]) - len(numeric_block_roots(N, 2, 2)["roots"])
        want = want_fn()
        check(f"(b2) N={N}: the DIFFERENCE of the numeric (3,3) and (2,2) counts is the exact "
              f"rung-3 count, {want} from {where}, so the subtraction that produces the N = 8 "
              f"rung numbers is itself validated where an exact rung count exists",
              got == want, f"{got} vs {want}")
    # (b3) DISJOINTNESS of the two rungs, exactly, at N = 7.  A difference is a rung count only
    #      if the two rungs share no coupling; a shared one sits in the smaller list, is
    #      subtracted away, and leaves the true rung count above the difference.  E9(b) settles
    #      this at N = 6 by a gcd of the two exact polynomials.  Both polynomials exist at N = 7
    #      as well, in this same run, so the gcd is available at the LAST N where it is, and
    #      taking it turns (b2) at N = 7 from a test of completeness and disjointness together
    #      into two separate results.  At N = 8 there is no exact polynomial and nothing checks
    #      it, which is why the note carries disjointness as an assumption there.
    t0 = time.time()
    common = Poly(gcd(EXACT_POLY[(7, 2)], EXACT_POLY[(7, 3)]), J)
    shared = common.count_roots(0, oo) - (1 if common.eval(0) == 0 else 0)
    check("(b3) N=7: the gcd of the exact rung-2 and rung-3 polynomials carries no positive "
          "root, so the two rungs share no coupling and the difference at N = 7 is a rung count "
          "for that reason and not by assumption", shared == 0,
          f"gcd of degree {common.degree()}, {shared} positive roots ({time.time() - t0:.0f} s)")
    # (c) the two blocks past the exact route.  A rung-l multiplet spans the blocks l .. N-l,
    #     so a block (p,p) sees only the rungs l <= min(p, N-p) and a bigger block sees strictly
    #     more.  The per-rung count is therefore a DIFFERENCE of nested block counts, and the
    #     nesting has to be checked rather than assumed.
    counts, roots, res8 = {}, {}, {}
    for (p, q) in [(2, 2), (3, 3), (4, 4)]:
        t0 = time.time()
        r = numeric_block_roots(8, p, q)
        counts[(p, q)], roots[(p, q)], res8[(p, q)] = len(r["roots"]), r["roots"], r
        print(f"      N=8 ({p},{q}) dim {r['dim']} ker {r['ker']}: {len(r['roots'])} accepted, "
              f"{len(r['rejected'])} rejected, |sin(phase)| {r['stray']:.1e}, smallest accepted "
              f"{min(r['roots']):.9f}, largest {max(r['roots']):.9f} ({time.time() - t0:.0f} s)")
        for j, mult, why in r["rejected"]:
            print(f"        REJECTED J = {j:.6f} (merged {mult}): {why}")
    # the deflation is numerical here where E8 does it exactly, so the dimension it removes is
    # asserted and not merely printed: ker C n ker A0 is the floor's own frozen space and must
    # come out at floor(N/2) on every block.
    check("(c0) N=8: the numerical deflation removes exactly floor(N/2) = 4 dimensions on each "
          "of the three blocks, which is the rational common kernel E8 computes over QQ at the "
          "N it can reach", all(res8[b]["ker"] == 4 for b in res8),
          f"{ {b: res8[b]['ker'] for b in res8} }")
    # the sign returned by det_sign is meaningful only while the phase stays near a multiple of
    # pi, and |sin| is that distance linearly.  The gate is the RATIO to the error model
    # sqrt(n)*pi*eps, held across every block this run brackets, n running from 96 to 4896: a
    # fixed threshold with eleven decades of headroom would pass whatever happened, while a
    # ratio that stays O(1) across two decades of n is the model holding.
    seen = list(res8.values()) + valblocks
    ratio = max(r["stray_over_model"] for r in seen)
    check("(c1) the determinant's phase stays on the real axis to the size the rounding model "
          "predicts, on every block bracketed in this run, so the sign each certification rests "
          "on is unambiguous", ratio < 10,
          f"worst |sin(phase)| / (sqrt(n)*pi*eps) = {ratio:.1f} over n from "
          f"{min(r['n'] for r in seen)} to {max(r['n'] for r in seen)}; worst absolute "
          f"{max(r['stray'] for r in seen):.1e}")
    ok23, miss23 = nested(roots[(2, 2)], roots[(3, 3)])
    ok34, miss34 = nested(roots[(3, 3)], roots[(4, 4)])
    check("(c) N=8: the roots of each block sit inside the next one's, value for value and one "
          "to one, so the three lists are nested and a difference of counts is a count of what "
          "the bigger block adds", ok23 and ok34,
          f"missing from (3,3): {miss23}, from (4,4): {miss34}")
    # (d) the ONE premise the word "bound" rests on, measured on each of the three blocks.  Two
    #     windows that both flip can only be one root counted twice if they both contain it,
    #     and then their gap is at most eps*(J_i + J_i+1); conversely two flipping windows that
    #     do not overlap contain two distinct roots.  So disjointness is not a proxy for the
    #     bound, it is equivalent to it.  A margin, so the number is printed and not the verdict
    #     alone; and 2.8 here means 1.4 combined window widths, since each window REACHES eps*J.
    #     The guard on the accepted count is not decoration: with one root there is no gap to
    #     measure and the ratio would be vacuously fine.
    wrs = [res8[b]["window_ratio"] for b in res8]
    wr = min([w for w in wrs if w is not None], default=None)
    check("(d) N=8: no two bracket windows overlap on any of the three blocks, which is "
          "SUFFICIENT for each count to be a lower bound and is the premise this run measures "
          "rather than assumes", wr is not None and None not in wrs and wr > 1,
          f"smallest gap between consecutive accepted roots is {wr:.1f}x the distance eps*J the "
          f"two windows reach toward each other, so {wr / 2:.1f}x their combined width"
          if wr is not None else "a block offered fewer than two roots, so there is no gap")
    # (e) the smallj cut, READ rather than trusted: the threshold 0.05 says nothing on its own,
    #     the separation it sits in says everything.  Only smear rejections enter, since a
    #     candidate the bracket refused at an ordinary coupling is not evidence about the cut.
    srs = [res8[b]["smear_ratio"] for b in res8]
    sr = min([s for s in srs if s is not None], default=None)
    par = sorted(set(m for b in res8 for m in res8[b]["parity"]))
    check("(e) N=8: the smallj cut has room, so nothing near it was decided by the threshold: "
          "the smallest ACCEPTED root sits decades above the largest candidate rejected as the "
          "z = 0 smear", sr is not None and sr > 100,
          f"ratio {sr:.0f}x, cut at 0.05, smallest accepted "
          f"{min(min(roots[b]) for b in roots):.9f}" if sr is not None
          else "no block produced a smear candidate, so there is no separation to read")
    print(f"  [--] detector health, not a premise: every accepted root merged {par} pencil "
          f"eigenvalues, and 2 is the healthy value since a root arrives as the pair +-z. A 1 "
          f"would mean a copy was lost, which is under-detection and cannot inflate a count; a "
          f"SPLIT pair is caught by (d). The z = 0 smear at (4,4) is the split case and shows "
          f"as merged 1 above.")
    # (f) one-sided on purpose.  Each count is a LOWER bound, so a later run finding one more
    #     root is new information and not a regression; a run finding one FEWER is.  Pinning
    #     equality would invert that and make the gate fail on the good news.  The recorded 15
    #     is read from E10(b) rather than typed, since an exact source for it exists; 40 and 52
    #     are typed, being the measurement itself.
    rec = {(2, 2): EXACT_SINGLET[(8, 2)], (3, 3): 40, (4, 4): 52}
    check("(f) N=8: every block count is at least the recorded one, 15 on (2,2), 40 on (3,3) "
          "and 52 on the middle block (4,4); the test is one-sided because the counts are",
          all(counts[b] >= rec[b] for b in rec), f"{counts}, recorded {rec}")
    rung3 = counts[(3, 3)] - counts[(2, 2)]
    rung4 = counts[(4, 4)] - counts[(3, 3)]
    # (g) is a different question from (f) and needs its own check, not a branch inside it.
    #     (f) asks whether the detector REGRESSED; (g) asks whether the note is CURRENT.  A run
    #     finding a 41st root on (3,3) passes (f), which is right, and must fail (g), because
    #     the table in the note is then wrong.  Making (g) conditional on counts == rec instead
    #     produced a check that could not fail: under that guard 52 - 40 = 12 is arithmetic, and
    #     40 - 15 = 25 follows from E10(b) having already passed.  The differences are printed
    #     rather than asserted for the same reason.
    check("(g) N=8: the counts are EXACTLY the ones the note records, so its table is current. "
          "This is not (f) again: (f) fails on a regression, this fails on new information, and "
          "a run that finds one more root should fail here and be written up",
          counts == rec, f"{counts} against {rec}; the differences are "
          f"(3,3) - (2,2) = {rung3} and (4,4) - (3,3) = {rung4}")
    print("  [--] a sign bracket certifies each ACCEPTED root and certifies nothing about a "
          "MISSED one, so every BLOCK count above is a verified lower bound. A DIFFERENCE of "
          "two lower bounds is bounded in neither direction, and the two above are differences: "
          "what the run supports with nothing further assumed is that at least 25 certified "
          "roots of (3,3) lie outside the accepted (2,2) list and at least 12 of (4,4) outside "
          "the accepted (3,3) list. Assume the SMALLER list complete and the difference becomes "
          "a lower bound on the rung; equality needs BOTH lists complete AND the two rungs to "
          "share no coupling. The second condition is settled exactly by a gcd at N = 6 in "
          "E9(b) and at N = 7 in (b3) above, which is the last N where both exact polynomials "
          "exist, and nothing settles it at N = 8. A growth law fitted to these numbers is "
          "fitted to counts that may be truncated.")


def main():
    print(f"exceptional couplings gate{'  (deep)' if DEEP else ''}"
          f"{'  (slow)' if SLOW else ''}{'  (rungs)' if RUNGS else ''}\n")
    for fn in (e1_builders, e2_reformulation, e3_n5_exact, e4_singlets, e5_polynomial,
               e6_shape, e7_bridge, e8_full_block, e9_higher_rungs,
               *([e10_slow_exact] if SLOW else ()),
               *([e11_numeric_rungs] if RUNGS else ())):
        fn()
        print()
    total = len(FAILURES)
    print(f"{len(CHECKS)} checks in {time.time() - T0:.1f} s")
    if total:
        print(f"exceptional couplings gate: {total} of {len(CHECKS)} FAILED")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    print(f"exceptional couplings gate: ALL GREEN, {len(CHECKS)} checks")


if __name__ == "__main__":
    main()
