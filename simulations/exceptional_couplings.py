"""Gate for the exceptional couplings of the frozen band.

The frozen depth on a band block is floor(N/2) for all but finitely many couplings, and
this gate is about the "finitely many": that set is NOT empty, and at a coupling in it one
more mode freezes.

Self-contained: every object is built from scratch here, including one construction from the
full 4^N Liouvillian, so agreement with the other gates is evidence and not a shared bug.

Run:
    python simulations/exceptional_couplings.py            # 35 checks, about a minute
    python simulations/exceptional_couplings.py --deep     # + N = 7 and the N = 6 whole-block
                                                           #   enumeration, 41 checks, ~9 min

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

Reading: PROOF_FROZEN_BAND_SO4 Section 5 (the pencil and the bridge) and Section 8
(what is open), ETA_CEILING_REDUCTION (the per-rung certificate this replaces the limit of).
"""
import itertools
import sys
import time

import numpy as np
from sympy import Matrix, Poly, QQ, Rational, factor_list, gcd, interpolate, oo, symbols
from sympy.polys.matrices import DomainMatrix

DEEP = "--deep" in sys.argv
PRIMES = (998244353, 1004535809)
FAILURES = []
T0 = time.time()


def check(label, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {label}" + (f"  ({detail})" if detail else ""))
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
    check("(b) so it is an eta multiplet seeded at l = 2, of eta-spin N/2 - 2",
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
        check("(f) N=6: the whole block (2,2) carries exactly six exceptional couplings, all "
              "simple, which is again the singlet count",
              r["count"] == 6 and r["simple"] and r["clean"], f"{[round(x, 6) for x in r['roots']]}")


def main():
    print(f"exceptional couplings gate{'  (deep)' if DEEP else ''}\n")
    for fn in (e1_builders, e2_reformulation, e3_n5_exact, e4_singlets, e5_polynomial,
               e6_shape, e7_bridge, e8_full_block):
        fn()
        print()
    total = len(FAILURES)
    print(f"{time.time() - T0:.1f} s")
    if total:
        print(f"exceptional couplings gate: {total} FAILED")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    print("exceptional couplings gate: ALL GREEN")


if __name__ == "__main__":
    main()
