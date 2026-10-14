"""The second centre-rung mode of the N = 8 chain as a polynomial in gamma.

Block (N/2, N/2) of the XXZ chain at Delta = 1, T(g) = -iA + g(D + N), A = [H, .], D + N = N - 2h on a cell
at Hamming distance h. The pencil is singular (a kernel at every g), so its kernel over Q(g) has a basis of
polynomial vectors (Kronecker); the degree of one is the length of a chain v_0, ..., v_eps with A v_0 = 0,
iA v_k = (D + N) v_{k-1}, (D + N) v_eps = 0. A and D + N are real, so v_k = i^k w_k with w_k real, and the
real pencil is A + mu (D + N), mu = i g. Transposition followed by R = right multiplication by X^N anticommutes
with A and with D + N and squares to 1 on the X-even half (EE + OO), splitting it: the mode whose chain starts in
H's commutant on the X-odd states has w_even symmetric on OO and w_odd antisymmetric on EE (the -1 eigenspace,
630 + 595); the orbit sum of the Hamming-N/2 cells is the length-0 chain of the +1 eigenspace.

Rows:
  D1  the Wong sequence F_0 = ker A (35-dimensional at N = 8), F_k = A^{-1}((D+N) F_{k-1}) in that class
      first meets ker(D+N), with dim F = 134 there, at
      eps = 1 for N = 4 (C + gC') and at eps = 267 for N = 8, modulo two primes; a meet at k is the nonzero top
      coefficient of a degree-k kernel vector, so the first meet is the minimal index.
  D2  (--border) the bordered determinant det [[A + mu B, u], [l^T, 0]] of the class's real pencil, at 640
      integer mu modulo 1000003, has degree 565 and the root mu = 0 with multiplicity 34 (gated). Reading, not gated: the adjugate of a
      pencil of rank one less than full is c(mu) times the mode times the left null vector, here the constant
      orbit sum, so 267 roots are the zeros of l on the mode and 565 - 267 - 34 = 264 lie at mu != 0, where the
      pencil holds at least one further mode.
  D3  (--pointwise) the class's kernel computed directly at mu = 1..600 modulo 1000003 is one-dimensional at
      every point, and each of four coordinates divided by a random functional of the kernel vector is a
      rational function of degree 267 (the larger of numerator and denominator degree), not less: the degree by a route without the chain.

Exact modular arithmetic throughout. D1 takes under two minutes; D2 about twenty and D3 about forty minutes on eight processes.
Usage: python simulations/f50_centre_mode_degree.py [--border] [--pointwise]
"""
import sys
import numpy as np
from concurrent.futures import ProcessPoolExecutor

FAILS = []


def check(name, ok):
    print(("ok    " if ok else "FAIL  ") + name, flush=True)
    if not ok:
        FAILS.append(name)


def pencil(N, p):
    """Class-O blocks mod p: A_e: Sym(OO)->Asym(OO), B_e: Sym(OO)->Sym(EE), A_o: Asym(EE)->Sym(EE),
    B_o: Asym(EE)->Asym(OO), in the bases u_i u_j^T (+/- transpose) of the pair vectors x +/- xbar."""
    st = [x for x in range(2 ** N) if bin(x).count("1") == N // 2]
    n = len(st); idx = {x: i for i, x in enumerate(st)}
    H = np.zeros((n, n), dtype=np.int64)
    for c, x in enumerate(st):
        for i in range(N - 1):
            bi, bj = (x >> i) & 1, (x >> (i + 1)) & 1
            H[c, c] += 1 if bi == bj else -1
            if bi != bj:
                H[idx[x ^ (1 << i) ^ (1 << (i + 1))], c] += 2
    full = 2 ** N - 1
    reps = [x for x in st if x < (x ^ full)]
    m = len(reps)
    Ev = np.zeros((n, m), dtype=np.int64); Ov = np.zeros((n, m), dtype=np.int64)
    for k, x in enumerate(reps):
        Ev[idx[x], k] = 1; Ev[idx[x ^ full], k] = 1
        Ov[idx[x], k] = 1; Ov[idx[x ^ full], k] = -1
    ham = np.array([[bin(x ^ y).count("1") for y in st] for x in st], dtype=np.int64)
    bcell = N - 2 * ham
    inv4 = pow(4, p - 2, p)
    sym = [(i, j) for i in range(m) for j in range(i, m)]
    asy = [(i, j) for i in range(m) for j in range(i + 1, m)]

    def ops(U, pairs, s):
        out = []
        for i, j in pairs:
            W = np.outer(U[:, i], U[:, j])
            out.append(W + W.T if (s and i != j) else (W if s else W - W.T))
        return out

    def mat(src, f, U, pairs):
        cols = []
        for W in src:
            M = (U.T @ f(W) @ U) % p * inv4 % p
            cols.append([M[i, j] for i, j in pairs])
        return np.array(cols, dtype=np.int64).T % p

    Aop = lambda W: H @ W - W @ H
    Bop = lambda W: bcell * W
    es, os_ = ops(Ov, sym, True), ops(Ev, asy, False)
    return mat(es, Aop, Ov, asy), mat(es, Bop, Ev, sym), mat(os_, Aop, Ev, sym), mat(os_, Bop, Ov, asy)


def rref(M, p):
    M = M.copy() % p; r = 0; piv = []
    rows, cols = M.shape
    for c in range(cols):
        if r == rows:
            break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0:
            continue
        k = r + nz[0]
        if k != r:
            M[[r, k]] = M[[k, r]]
        M[r] = M[r] * pow(int(M[r, c]), p - 2, p) % p
        others = np.nonzero(M[:, c])[0]
        others = others[others != r]
        if len(others):
            M[others] = (M[others] - np.outer(M[others, c], M[r])) % p
        piv.append(c); r += 1
    return M[:r], piv


def nullspace(M, p):
    R, piv = rref(M, p); cols = M.shape[1]
    ps = set(piv)
    free = [c for c in range(cols) if c not in ps]
    Z = np.zeros((cols, len(free)), dtype=np.int64)
    for t, f in enumerate(free):
        Z[f, t] = 1
        for r, c in enumerate(piv):
            Z[c, t] = (-R[r, f]) % p
    return Z


def rank(M, p):
    return 0 if M.shape[1] == 0 else len(rref(M, p)[1])


def inv_mod(M, p):
    k = M.shape[0]
    R, piv = rref(np.hstack([M % p, np.eye(k, dtype=np.int64)]), p)
    assert piv[:k] == list(range(k))
    return R[:, k:]


def first_meet(N, p, dmax):
    A_e, B_e, A_o, B_o = pencil(N, p)
    _, rows_o = rref(A_o.T.copy(), p)
    Lo = inv_mod(A_o[rows_o], p)
    Y = nullspace(A_o.T.copy(), p).T
    _, cols_e = rref(A_e.copy(), p)
    Ie = inv_mod(A_e[:, cols_e], p)
    K0 = nullspace(A_e, p)
    F = K0
    for d in range(dmax + 1):
        BF = (B_e if d % 2 == 0 else B_o) @ F % p
        if F.shape[1] - rank(BF, p):
            return d, F.shape[1], K0.shape[1]
        if d % 2 == 0:
            S = BF @ nullspace(Y @ BF % p, p) % p
            F = Lo @ S[rows_o] % p
        else:
            W = np.zeros((A_e.shape[1], BF.shape[1]), dtype=np.int64)
            W[cols_e] = Ie @ BF % p
            F = rref(np.hstack([W, K0]).T.copy(), p)[0].T
    return None, F.shape[1], K0.shape[1]


P_B = 1000003


def _border_det(g):
    N = 8
    A_e, B_e, A_o, B_o = _border_det.blocks
    s0, r0 = A_e.shape[1], A_e.shape[0]
    n = s0 + A_o.shape[1]
    M = np.zeros((n + 1, n + 1), dtype=np.int64)
    M[:r0, :s0] = A_e; M[r0:n, s0:n] = A_o
    M[:r0, s0:n] = (M[:r0, s0:n] + g * B_o) % P_B
    M[r0:n, :s0] = (M[r0:n, :s0] + g * B_e) % P_B
    rng = np.random.default_rng(11)
    M[:n, n] = rng.integers(1, P_B, n); M[n, :n] = rng.integers(1, P_B, n)
    d = 1
    for c in range(n + 1):
        nz = np.nonzero(M[c:, c])[0]
        if len(nz) == 0:
            return g, 0
        r = c + nz[0]
        if r != c:
            M[[c, r]] = M[[r, c]]; d = -d
        piv = int(M[c, c]); d = d * piv % P_B
        rows = c + 1 + np.nonzero(M[c + 1:, c])[0]
        if len(rows):
            f = M[rows, c] * pow(piv, P_B - 2, P_B) % P_B
            M[rows, c:] = (M[rows, c:] - np.outer(f, M[c, c:])) % P_B
    return g, d % P_B


def _class_pencil(p):
    A_e, B_e, A_o, B_o = pencil(8, p)
    s0, r0 = A_e.shape[1], A_e.shape[0]
    n = s0 + A_o.shape[1]
    A = np.zeros((n, n), dtype=np.int64); B = np.zeros_like(A)
    A[:r0, :s0] = A_e; A[r0:, s0:] = A_o; B[:r0, s0:] = B_o; B[r0:, :s0] = B_e
    return A, B


def _kernel_at(g):
    A, B = _kernel_at.AB
    Z = nullspace((A + g * B) % P_B, P_B)
    if Z.shape[1] != 1:
        return g, None
    v = Z[:, 0]
    ell = np.random.default_rng(7).integers(1, P_B, v.shape[0])
    inv = pow(int(ell @ v % P_B), P_B - 2, P_B)
    return g, [int(v[c]) * inv % P_B for c in (0, 315, 633, 1224)]


def _init_kernel():
    _kernel_at.AB = _class_pencil(P_B)


def _init():
    _border_det.blocks = pencil(8, P_B)


if __name__ == "__main__":
    print("D1  the chain's first meet with ker(D + N)", flush=True)
    d4, _, _ = first_meet(4, 1000003, 6)
    check(f"N = 4: first meet at eps = {d4}, expected 1 (C + gC')", d4 == 1)
    for p in (1000003, 1000033):
        d8, dimF, c0 = first_meet(8, p, 400)
        check(f"N = 8 mod {p}: commutant 35 ({c0}), first meet at eps = {d8} with dim F = {dimF}, expected 267 and 134",
              (d8, dimF, c0) == (267, 134, 35))
    if "--border" in sys.argv:
        print("D2  the bordered determinant of the class's pencil", flush=True)
        with ProcessPoolExecutor(max_workers=8, initializer=_init) as ex:
            pts = sorted(ex.map(_border_det, range(1, 641)))
        gs = [g for g, _ in pts]; ys = [y for _, y in pts]

        def fits(m):
            V = np.array([[pow(g, k, P_B) for k in range(m + 1)] + [y] for g, y in zip(gs, ys)], dtype=np.int64)
            return rank(V, P_B) == m + 1
        lo, hi = 0, len(gs) - 2
        while lo < hi:
            mid = (lo + hi) // 2
            lo, hi = (lo, mid) if fits(mid) else (mid + 1, hi)
        V = np.array([[pow(g, k, P_B) for k in range(lo + 1)] + [y] for g, y in zip(gs[: lo + 1], ys[: lo + 1])],
                     dtype=np.int64)
        coef = [int(r[-1]) for r in rref(V, P_B)[0]]
        mult0 = next(k for k, c in enumerate(coef) if c)
        check(f"degree {lo} over 640 points, root mu = 0 of multiplicity {mult0}; expected 565 and 34", (lo, mult0) == (565, 34))
    if "--pointwise" in sys.argv:
        print("D3  the degree from kernels at single points mu = i g", flush=True)
        with ProcessPoolExecutor(max_workers=8, initializer=_init_kernel) as ex:
            pts = sorted(ex.map(_kernel_at, range(1, 601)))
        check("the kernel is one-dimensional at mu = 1..600", all(v is not None for _, v in pts))
        pts = [(g, v) for g, v in pts if v is not None]
        gs = [g for g, _ in pts]
        for j, c in enumerate((0, 315, 633, 1224)):
            r = np.array([v[j] for _, v in pts], dtype=np.int64)

            def fits_ratio(m):
                V = np.array([[pow(g, k, P_B) for k in range(m + 1)] for g in gs], dtype=np.int64)
                M = np.hstack([V, (-(r[:, None] * V)) % P_B])
                return M.shape[1] - rank(M, P_B) > 0
            check(f"coordinate {c}: a ratio of degree 267 fits, 266 does not", fits_ratio(267) and not fits_ratio(266))
    print("\nALL GATES PASS" if not FAILS else f"\n{len(FAILS)} FAILED")
    sys.exit(1 if FAILS else 0)
