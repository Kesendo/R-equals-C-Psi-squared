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
  E6  the shape of the failure.  The extra mode sits on the diagonal rungs l .. N-l with
      l = 2 and on no side line, which is an eta multiplet that is a SPIN SINGLET.
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
  E10 the two exact counts the other modes are too short for (slow).  The same Sturm route as
      E5, at the rung-3 singlets of the block (3,3) at N = 7 and at the rung-2 singlets at
      N = 8, which is where the exact route runs out.  Kept behind its own flag so that the
      numbers it produces are in the record without lengthening --deep.
  E11 past the exact route: the middle block at N = 8, numerically (rungs).  Detection is the
      finite spectrum of the same pencil rather than a polynomial, and every candidate is
      verified on its own by a SIGN BRACKET of the determinant, so an accepted root is
      certified real.  Nothing certifies that no root was MISSED, so every count here is a
      verified LOWER bound and the block says so.  Three things make that word earn its keep
      rather than be asserted.  The expected values it validates against are READ OUT of the
      exact blocks of the same run and none is typed here, which is why --rungs turns on --deep
      and --slow.  The subtraction that turns block counts into rung counts is itself validated
      wherever an exact rung count exists, at N = 6 and N = 7, and no such check can exist for
      rung 4.  And the two premises a lower bound silently needs are measured as margins: that
      each root's +-z pair MERGED, since a split pair is counted twice and inflates the count,
      and that no two bracket windows OVERLAP, since then one root certifies two candidates.
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
    for (N, p, q) in [(4, 2, 2), (4, 3, 1), (5, 2, 2), (5, 3, 3), (5, 2, 0)]:
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
        if N == 7:
            check("(d) N=7: the smallest exceptional coupling is the value the note tabulates",
                  abs(roots[0] - 0.952056678) < 1e-9, f"{roots[0]:.12f}")
        if N == 6:
            check("(f) N=6: the smallest exceptional coupling is the value the note tabulates",
                  abs(roots[0] - 0.749042443688) < 1e-11, f"{roots[0]:.12f}")
            degs = sorted(Poly(f, J).degree() for f, _ in factor_list(g.as_expr())[1]
                          if any(float(r) > 1e-9 for r in Poly(f, J).real_roots() if r.is_real))
            check("(e) N=6: BOTH factors carrying a coupling have degree 12, so every one of "
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
    check("(c) and a SPIN SINGLET, since S+ commutes with L and no side line is raised",
          all(counts[(p, q)] == 3 for (p, q) in band(6) if p != q))


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
            counts[(p, q)] = int(np.sum(s < 1e-8 * max(s[0], 1.0)))
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
        also33 = also33 and int(np.sum(sv < 1e-8 * sv[0])) > 3
    check("(d) the six rung-2 couplings raise the block (3,3) as well, so that block carries "
          "the six AND the four and the exceptional set at N = 6 has at least ten members: "
          "the count on (2,2) is the rung-2 count and not the count",
          also33 and n3 == 4 and shared == 0, f"six rung-2 roots raise (3,3): {also33}")
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
        check("(e) N=7: the smallest coupling of the block (3,3) collapses the fourth "
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
    """The sign of det(C_G + i x A_G), real because E1(d) makes the polynomial even in z.
    Carried as a phase, so the modulus may overflow freely; the second value reports how far
    the phase strays from the real axis, which is a measured check and not an assumption."""
    lu, piv = sla.lu_factor(Cg + 1j * x * Ag)
    d = np.diag(lu)
    phase = float(np.sum(np.angle(d))) + np.pi * (int(np.sum(piv != np.arange(len(piv)))) % 2)
    c = np.cos(phase)
    return (1.0 if c > 0 else -1.0), abs(abs(c) - 1.0)


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
    count is a LOWER bound -- but only if nothing is counted TWICE, and that is not free.  Two
    premises stand between this routine and the word "bound", and both are returned as measured
    margins rather than assumed:

      * PARITY.  The pencil spectrum carries every root as the pair +-z, so |z| arrives twice
        and the two copies are merged.  If a pair is resolved to worse than the merge threshold
        it splits into two candidates; if the split is narrower than the bracket window, BOTH
        windows straddle the same true root, both flip, and one root is counted twice.  So
        `mults` reports how many pencil eigenvalues each accepted root merged, and anything but
        an even number at least 2 breaks the bound.  This is not hypothetical: the z = 0 smear
        at N = 8 (4,4) arrives as two unmerged copies, and only the smallj rule catches them.
      * DISJOINTNESS.  Two accepted windows that overlap can certify each other from one root.
        `window_ratio` is the smallest gap between consecutive accepted roots divided by the
        width of the two windows it must separate, so a value above 1 is the premise holding
        and the number says by how much.

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
    # the two premises, as measured margins.  A gap must separate two windows of half-width
    # eps*J each, so the width to clear is eps*(J_i + J_i+1).
    ratios = [(kept[i + 1] - kept[i]) / (eps * (kept[i] + kept[i + 1]))
              for i in range(len(kept) - 1)]
    below = [j for j, _, _ in rejected]
    return dict(roots=kept, mults=mults, rejected=rejected, ker=ker, dim=C.shape[0], stray=stray,
                parity_ok=all(m >= 2 and m % 2 == 0 for m in mults),
                window_ratio=min(ratios) if ratios else float("inf"),
                smear_ratio=(min(kept) / max(below)) if kept and below else float("inf"))


def nested(small, large, tol=1e-5):
    """Every root of the smaller block matched to one of the larger, ONE TO ONE.

    The matching has to be injective or the nesting proves less than it looks: two roots of the
    small list both matching one root of the large list would pass a per-root test while the
    set difference is one larger than the count difference.  The tolerance is the bracket's own
    reach, 1e-5 relative, since that is how well an accepted value is placed; anything tighter
    would trust the values better than they are certified."""
    taken, miss = set(), []
    for j in small:
        hit = next((i for i, k in enumerate(large)
                    if i not in taken and abs(j - k) <= tol * max(j, 1.0)), None)
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
    same = (len(r6["roots"]) == len(exact6)
            and all(abs(a - b) < 1e-9 * max(a, 1.0) for a, b in zip(r6["roots"], exact6)))
    check("(a) N=6 (3,3): the pencil detector returns exactly the ten roots the exact "
          "polynomials of E5 and E9 give, value for value and not merely ten of them",
          same and len(exact6) == 10,
          f"{len(r6['roots'])} numeric vs {len(exact6)} exact, worst relative difference "
          f"{max((abs(a - b) / max(a, 1.0) for a, b in zip(r6['roots'], exact6)), default=0):.1e}"
          f", {len(r6['rejected'])} rejected")
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
    for (N, p, q, want_fn, where) in table:
        try:
            want = want_fn()
        except KeyError as miss:
            check(f"(b) N={N} ({p},{q}): the exact count this validation needs was computed in "
                  f"this run", False, f"missing {miss}")
            continue
        r = numeric_block_roots(N, p, q)
        check(f"(b) N={N} ({p},{q}): the numeric count is the exact one, {want} from {where}",
              len(r["roots"]) == want, f"{len(r['roots'])}, ker {r['ker']} of {r['dim']}, "
              f"worst phase stray {r['stray']:.1e}, {len(r['rejected'])} rejected, "
              f"window margin {r['window_ratio']:.1f}x")
    # (b2) the METHOD, not the detector: wherever an exact rung count exists, the difference of
    #      two numeric block counts must reproduce it.  This is the only place the subtraction
    #      itself is tested, and it tests both premises it rests on at once, completeness of the
    #      smaller list and disjointness of the two rungs.  It exists at N = 6 and N = 7 and can
    #      exist for no rung 4 anywhere: at N = 8 the rung-4 multiplet has eta-spin 0 and lives
    #      in (4,4) alone, so the twelve has no second route of any kind.
    for (N, want_fn, where) in [(6, lambda: EXACT_SINGLET[(6, 3)], "E9(a)"),
                                (7, lambda: EXACT_SINGLET[(7, 3)], "E10(a)")]:
        got = len(numeric_block_roots(N, 3, 3)["roots"]) - len(numeric_block_roots(N, 2, 2)["roots"])
        want = want_fn()
        check(f"(b2) N={N}: the DIFFERENCE of the numeric (3,3) and (2,2) counts is the exact "
              f"rung-3 count, {want} from {where}, so the subtraction that produces the N = 8 "
              f"rung numbers is itself validated where an exact rung count exists",
              got == want, f"{got} vs {want}")
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
              f"{len(r['rejected'])} rejected, phase stray {r['stray']:.1e}, smallest accepted "
              f"{min(r['roots']):.9f}, largest {max(r['roots']):.9f} ({time.time() - t0:.0f} s)")
        for j, mult, why in r["rejected"]:
            print(f"        REJECTED J = {j:.6f} (merged {mult}): {why}")
    ok23, miss23 = nested(roots[(2, 2)], roots[(3, 3)])
    ok34, miss34 = nested(roots[(3, 3)], roots[(4, 4)])
    check("(c) N=8: the roots of each block sit inside the next one's, value for value and one "
          "to one, so the three lists are nested and a difference of counts is a count of what "
          "the bigger block adds", ok23 and ok34,
          f"missing from (3,3): {miss23}, from (4,4): {miss34}")
    # (d) the two premises the word "bound" rests on, measured on each of the three blocks.
    #     PARITY first: a root arrives as the pair +-z, and a pair that fails to merge is
    #     counted twice, which breaks the bound in the direction that inflates it.  Then
    #     DISJOINTNESS: overlapping windows let one root certify two candidates.  Both are
    #     margins, so the numbers are printed and not just the verdict.
    par = all(res8[b]["parity_ok"] for b in res8)
    wr = min(res8[b]["window_ratio"] for b in res8)
    sr = min(res8[b]["smear_ratio"] for b in res8)
    check("(d) N=8: every accepted root merged an EVEN number of pencil eigenvalues, at least "
          "two, so no root is counted twice, and no two bracket windows overlap", par and wr > 1,
          f"multiplicities {sorted(set(m for b in res8 for m in res8[b]['mults']))}, smallest "
          f"gap between consecutive accepted roots is {wr:.1f}x the width of the two windows "
          f"it separates")
    check("(e) N=8: the smallj cut has room, so nothing near it was decided by the threshold: "
          "the smallest ACCEPTED root and the largest REJECTED one are decades apart",
          sr > 100, f"ratio {sr:.0f}x, cut at 0.05, smallest accepted "
          f"{min(min(roots[b]) for b in roots):.9f}")
    # (f) one-sided on purpose.  Each count is a LOWER bound, so a later run finding one more
    #     root is new information and not a regression; a run finding one FEWER is.  Pinning
    #     equality would invert that and make the gate fail on the good news.
    rec = {(2, 2): 15, (3, 3): 40, (4, 4): 52}
    check("(f) N=8: every block count is at least the recorded one, 15 on (2,2), 40 on (3,3) "
          "and 52 on the middle block (4,4); the test is one-sided because the counts are",
          all(counts[b] >= rec[b] for b in rec), f"{counts}, recorded {rec}")
    if any(counts[b] > rec[b] for b in rec):
        print(f"  [!!] a block count EXCEEDS the recorded one: {counts} against {rec}. That is "
              f"new information, not a failure, and the note's numbers need updating.")
    rung3 = counts[(3, 3)] - counts[(2, 2)]
    rung4 = counts[(4, 4)] - counts[(3, 3)]
    check("(g) N=8: the differences are 25 and 12, the first reproducing the rung-3 count "
          "measured at this N on the singlet space of (3,3) rather than on the whole block, "
          "and the second a rung no block below (4,4) could see",
          rung3 == 25 and rung4 == 12, f"(3,3) - (2,2) = {rung3}, (4,4) - (3,3) = {rung4}")
    print("  [--] a sign bracket certifies each ACCEPTED root and certifies nothing about a "
          "MISSED one, so every BLOCK count above is a verified lower bound. A DIFFERENCE of "
          "two lower bounds is bounded in neither direction, and the two above are differences: "
          "what the run supports with nothing further assumed is that at least 25 certified "
          "roots of (3,3) lie outside the accepted (2,2) list and at least 12 of (4,4) outside "
          "the accepted (3,3) list. Assume the SMALLER list complete and the difference becomes "
          "a lower bound on the rung; equality needs BOTH lists complete AND the two rungs to "
          "share no coupling, which E9(b) checks exactly at N = 6 by a gcd and which nothing "
          "checks here. A growth law fitted to these numbers is fitted to counts that may be "
          "truncated.")


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
