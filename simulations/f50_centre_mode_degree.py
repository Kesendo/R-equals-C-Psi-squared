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
  D4  (--real, floating point) which roots lie at real g. For real g the class's pencil is equivalent to the real
      pencil [[A_e, -g B_o], [g B_e, A_o]]; bordered with two random borders, QZ gives 298 eigenvalues common to
      both, their matching distances six decades and more below the unshared ones' (34 at 0, 264 elsewhere,
      the count D2's reading gives), and 16 of the 264 come out exactly real, +-g at eight couplings, the others at least a tenth of
      their modulus off the axis. (The X-odd half is D5's.) On the full block (4,4), 4900 dimensions, the four smallest singular values of
      T(g) at each lie six decades and more below the fifth (two modes more than the generic two, one per
      class, the X-odd half adding none). Off each coupling by g t, t = 1e-3 and 1e-4, the two eigenvalues of T beside
      the generic two scale linearly in t (ratio 10 within 5 percent), so the pair grows linearly off g rather
      than branching, and at g exactly four eigenvalues lie below 1e-6 (a Jordan block would split at rounding by about
      sqrt(eps |T|) ~ 1e-7), so the four modes are the whole eigenvalue; the pair classifies the coupling: an oscillating pair on the
      line -8g whose frequency passes through 0, or a real palindromic pair -8g +- r whose r passes through 0.

  D5  (--odd) the X-odd half, U on E><O and V on O><E: A keeps each side, R (-1 on E><O, +1 on O><E) reflects
      D + N, so D + N maps one side to the other, and transposition carries one side onto the other. Where h_E
      and h_O share no eigenvalue the half is singular exactly at g = 1/|omega|, omega a real nonzero eigenvalue
      of G = t S^-1 (D + N) on O><E, S(U) = h_E U - U h_O. N = 4: g = 1 twice (the Proposition), the kernel
      two-dimensional there and the two pairs changing sign in opposite directions. N = 8, modulo two
      primes: D + 8 has rank 595 between the sides and G has 566 nonzero eigenvalues (the stable rank of its
      powers). Floating point: the same 566 (a magnitude gap of at least six decades), 48 real (imaginary part
      exactly 0), 48 distinct couplings, split 24 + 24 by the chain reflection, which commutes with G (exactly, in
      integer bases), none among
      D4's eight; the half rebuilt from the cells of the block (4,4), [H, .] and D + 8 projected onto it, equals
      the construction and stays invariant, exactly in integer bases (a second route to the matrices read below);
      the non-real eigenvalue nearest the axis is genuinely non-real (along real g the half's smallest singular
      value stays far above eps |T| near it). Each coupling is an exceptional point: nullity one, which with the
      even multiplicity sigma(rho) = rho^T X^N forces makes a 2 x 2 Jordan block (two eigenvalues at the
      sqrt(eps |T|) rounding scale, far above eps |T|), the pair +-mu growing as the square root of the offset
      (exponent within 0.05 of 1/2, a crossing giving 1), real on one side and oscillating on the other (offsets taken away from the nearest other
      coupling: 0.295030 and 0.295382 lie 1.2e-3 apart, relatively).

  D6  (--krein) the real couplings of the X-odd half from below. S' = t S t and M = t (D + N) are both symmetric,
      so G = S'^-1 M is self-adjoint for both (S' G = M, M G = M S'^-1 M). M is the pair form
      V -> sum beta_ij V_ij V_ji with beta_ij = N - 2 h(x_i, x_j), exactly (integer bases), and M's spectrum has
      the inertia that form predicts: signature C(N, N/2)/2 whatever H is. The construction gives D5's 48 couplings,
      and so does the half's pencil without the elimination. The M-signs of the real eigenvalues (all simple) and
      M on the zero root space add up to that signature, complex pairs carrying none: chain N = 4, 6, 8 with real
      counts 2, 14, 48, zero-root parts 1, 0, 1, so 2, 10, 34 forced; at N = 8 per reflection sector 24, 24 real
      against 14, 20 forced. At N = 6 there is no rung cell and M is invertible, so every H with S invertible has
      at least 10: forty random connected XXZ graphs give 10 to 20, never fewer. Where S is singular the bound has
      no hold, and the half's pencil reads 15 on the chain at Delta = 1/2 and 8 on a five-site chain beside an
      isolated site. Real means exactly on the axis, the rest lying at least 1e-5 of the modulus off it (a split
      double root would sit between).

D1 to D3 and D5's counts in exact modular arithmetic, D4 and D5's couplings in floating point, D6's inertia of M
by eigvalsh against the pair form's exact count, its sign characteristics in floating point. D1 takes under two minutes; D2 about twenty and D3 about
forty minutes on eight processes; D4 about fifteen minutes; D5 about ten minutes; D6 about two minutes.
Usage: python simulations/f50_centre_mode_degree.py [--border] [--pointwise] [--real] [--odd] [--krein]
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


def _float_blocks():
    """Class-O blocks over the reals (the pencil() construction without the modulus), N = 8."""
    N = 8
    st = [x for x in range(2 ** N) if bin(x).count("1") == N // 2]
    n = len(st); idx = {x: i for i, x in enumerate(st)}
    H = np.zeros((n, n))
    for c, x in enumerate(st):
        for i in range(N - 1):
            bi, bj = (x >> i) & 1, (x >> (i + 1)) & 1
            H[c, c] += 1 if bi == bj else -1
            if bi != bj:
                H[idx[x ^ (1 << i) ^ (1 << (i + 1))], c] += 2
    full = 2 ** N - 1
    reps = [x for x in st if x < (x ^ full)]
    m = len(reps)
    Ev = np.zeros((n, m)); Ov = np.zeros((n, m))
    for k, x in enumerate(reps):
        Ev[idx[x], k] = 1; Ev[idx[x ^ full], k] = 1
        Ov[idx[x], k] = 1; Ov[idx[x ^ full], k] = -1
    ham = np.array([[bin(x ^ y).count("1") for y in st] for x in st])
    bcell = (N - 2 * ham).astype(float)
    sym = [(i, j) for i in range(m) for j in range(i, m)]
    asy = [(i, j) for i in range(m) for j in range(i + 1, m)]

    def ops(U, pairs, s):
        out = []
        for i, j in pairs:
            W = np.outer(U[:, i], U[:, j])
            out.append(W + W.T if (s and i != j) else (W if s else W - W.T))
        return out

    def mat(src, f, U, pairs):
        I = np.array([q[0] for q in pairs]); J = np.array([q[1] for q in pairs])
        return np.array([((U.T @ f(W) @ U) / 4)[I, J] for W in src]).T

    Aop = lambda W: H @ W - W @ H
    Bop = lambda W: bcell * W
    es, os_ = ops(Ov, sym, True), ops(Ev, asy, False)
    blocks = mat(es, Aop, Ov, asy), mat(es, Bop, Ev, sym), mat(os_, Aop, Ev, sym), mat(os_, Bop, Ov, asy)
    return blocks, H, ham


def real_couplings():
    import scipy.linalg as sl
    (A_e, B_e, A_o, B_o), H, ham = _float_blocks()
    r0, s0 = A_e.shape
    M0 = np.zeros((r0 + A_o.shape[0], s0 + A_o.shape[1])); M1 = np.zeros_like(M0)
    M0[:r0, :s0] = A_e; M0[r0:, s0:] = A_o
    M1[:r0, s0:] = -B_o; M1[r0:, :s0] = B_e
    nn = M0.shape[0]

    def roots(seed):
        rng = np.random.default_rng(seed)
        P0 = np.zeros((nn + 1, nn + 1)); P1 = np.zeros_like(P0)
        P0[:nn, :nn] = M0; P1[:nn, :nn] = M1
        P0[:nn, nn] = rng.standard_normal(nn); P0[nn, :nn] = rng.standard_normal(nn)
        w = sl.eigvals(P0, -P1)
        return w[np.isfinite(w)]
    ra, rb = roots(1), roots(2)
    # a root of c does not depend on the border: sort the distances to the other border's spectrum and cut at the gap
    d = np.array([np.min(np.abs(rb - z)) / max(1.0, abs(z)) for z in ra])
    o = np.argsort(d)
    eps = np.finfo(float).eps                       # distances below the rounding unit count as the rounding unit
    k = int(np.argmax(np.log(np.maximum(d[o][1:], eps)) - np.log(np.maximum(d[o][:-1], eps)))) + 1
    match_gap = (d[o][k - 1], d[o][k])
    common = ra[o[:k]]
    a = np.sort(np.abs(common))
    k0 = int(np.argmax(np.log(np.maximum(a[1:], eps)) - np.log(np.maximum(a[:-1], eps)))) + 1
    zero_gap = (a[k0 - 1], a[k0])
    nz = common[np.abs(common) >= a[k0]]
    real = np.sort(nz[nz.imag == 0.0].real)          # QZ on a real pencil returns a simple real root with imag exactly 0
    others = nz[nz.imag != 0.0]
    margin = np.min(np.abs(others.imag) / np.abs(others))
    return (len(common), k0, len(nz), real, match_gap, zero_gap, margin), H, ham


def _init():
    _border_det.blocks = pencil(8, P_B)


def odd_half(N, p=None):
    """The X-odd half of the block (N/2, N/2): operators U on |E><O| and V on |O><E| (35 x 35 at N = 8, in
    the bases e_k = x + xbar, o_k = x - xbar of the popcount-N/2 states, normalized when p is None).
    A acts within each side, A1(U) = h_E U - U h_O; D + N maps O><E to E><O by B (R anticommutes with it),
    and transposition t carries the E><O side to the O><E side. Returns (A1, B, t, H, P_E, P_O) with the
    chain reflection restricted to E and O; mod p (integer, unnormalized) or in floating point."""
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
    reps = [x for x in st if x < (x ^ full)]; m = len(reps)
    e = np.zeros((n, m), dtype=np.int64); o = np.zeros((n, m), dtype=np.int64)
    for k, x in enumerate(reps):
        e[idx[x], k] = e[idx[x ^ full], k] = 1
        o[idx[x], k] = 1; o[idx[x ^ full], k] = -1
    ham = np.array([[bin(x ^ y).count("1") for y in st] for x in st])
    bc = (N - 2 * ham).astype(np.int64)
    rev = lambda x: int(format(x, f"0{N}b")[::-1], 2)
    P = np.zeros((n, n), dtype=np.int64)
    for x in st:
        P[idx[rev(x)], idx[x]] = 1
    if p is None:
        e, o, H, bc, P = e / np.sqrt(2), o / np.sqrt(2), H.astype(float), bc.astype(float), P.astype(float)
    hE, hO = e.T @ H @ e, o.T @ H @ o
    Im = np.eye(m, dtype=e.dtype)
    A1 = np.kron(Im, hE) - np.kron(hO.T, Im)            # column-major vec: U -> hE U - U hO
    B = np.zeros((m * m, m * m), dtype=e.dtype)
    for k in range(m * m):
        V = np.zeros(m * m, dtype=e.dtype); V[k] = 1; V = V.reshape((m, m), order="F")
        B[:, k] = (e.T @ (bc * (o @ V @ e.T)) @ o).reshape(-1, order="F")
    t = np.zeros((m * m, m * m), dtype=e.dtype)
    for k in range(m * m):
        t[(k // m) + (k % m) * m, k] = 1
    if p is not None:
        A1, B = A1 % p, B % p
    return A1, B, t, (H, e, o, bc), e.T @ P @ e, o.T @ P @ o


def _odd_svd(args):
    """At a coupling g of the X-odd half T(g) = -iA + g(D + N): its two smallest singular values and its norm, its
    three smallest eigenvalue moduli, and the two eigenvalues nearest 0 at g(1 + s d), d = 1e-5, 1e-6, -1e-5, with
    s = +-1 pointing away from the nearest other coupling (two of them lie 1.2e-3 apart, relatively)."""
    # Only the bottom of the spectrum is needed, so everything goes through one LU factorization per matrix:
    # the two smallest singular values by Lanczos on (T^H T)^-1, the eigenvalues nearest 0 by shift-invert at 0.
    import scipy.linalg as sl
    from scipy.sparse.linalg import LinearOperator, eigsh, eigs
    g, sgn = args
    Aodd, Dodd = _odd_svd.mats
    T = -1j * Aodd + g * Dodd
    n = T.shape[0]
    lu = sl.lu_factor(T)
    inv_gram = LinearOperator((n, n), dtype=complex,
                              matvec=lambda x: sl.lu_solve(lu, sl.lu_solve(lu, x, trans=2)))
    s = np.sort(1 / np.sqrt(np.abs(eigsh(inv_gram, k=2, which="LM", return_eigenvectors=False))))
    norm = np.sqrt(eigsh(LinearOperator((n, n), dtype=complex, matvec=lambda x: T.conj().T @ (T @ x)),
                         k=1, which="LM", return_eigenvectors=False)[0].real)
    near0 = lambda M, k: eigs(M, k=k, sigma=0, which="LM", return_eigenvectors=False)
    lam = np.sort(np.abs(near0(T, 3)))
    pairs = []
    for d in (1e-5, 1e-6, -1e-5):
        w = near0(-1j * Aodd + g * (1 + sgn * d) * Dodd, 4)
        pairs.append(w[np.argsort(np.abs(w))[:2]])
    return g, (s[0], s[1], norm, lam, pairs)


def _init_odd():
    A1, B, t, _, _, _ = odd_half(8)
    mm = A1.shape[0]; Z = np.zeros((mm, mm))
    _odd_svd.mats = (np.block([[A1, Z], [Z, -(t @ A1 @ t)]]), np.block([[Z, B], [t @ B @ t, Z]]))


def krein_half(N, bonds):
    """The X-odd half on |O><E| for XXZ bonds (i, j, J, Delta), normalized bases, as the symmetric pair (S', M):
    S' = t S t from the Sylvester operator, M[(i, j), (j, i)] = w_ij = N - 2 h(x_i, x_j) from the dephasing,
    G = S'^-1 M. Column-major V, as odd_half. Also the pair weights w."""
    st = [x for x in range(2 ** N) if bin(x).count("1") == N // 2]
    n = len(st); idx = {x: i for i, x in enumerate(st)}
    H = np.zeros((n, n))
    for c, x in enumerate(st):
        for i, j, J, Dl in bonds:
            bi, bj = (x >> i) & 1, (x >> j) & 1
            H[c, c] += J * Dl * (1 if bi == bj else -1)
            if bi != bj:
                H[idx[x ^ (1 << i) ^ (1 << j)], c] += 2 * J
    full = 2 ** N - 1
    reps = [x for x in st if x < (x ^ full)]; m = len(reps)
    e = np.zeros((n, m)); o = np.zeros((n, m))
    for k, x in enumerate(reps):
        e[idx[x], k] = e[idx[x ^ full], k] = 1 / np.sqrt(2)
        o[idx[x], k] = 1 / np.sqrt(2); o[idx[x ^ full], k] = -1 / np.sqrt(2)
    hE, hO = e.T @ H @ e, o.T @ H @ o
    Im = np.eye(m)
    S = np.kron(Im, hE) - np.kron(hO.T, Im)
    t = np.zeros((m * m, m * m))
    for k in range(m * m):
        t[(k // m) + (k % m) * m, k] = 1
    w = np.array([[N - 2 * bin(x ^ y).count("1") for y in reps] for x in reps])
    M = np.zeros((m * m, m * m))
    for i in range(m):
        for j in range(m):
            M[i + j * m, j + i * m] = w[i, j]
    return t @ S @ t, M, w


def krein_count(Sp, M, P=None):
    """Real nonzero eigenvalues of G = S'^-1 M (restricted to the columns of P if given): their number, the sum of
    their sign characteristics sign(v^T M v), and the inertia of M on G's zero root space and on the whole."""
    import scipy.linalg as sl
    if P is not None:
        Sp, M = P.T @ Sp @ P, P.T @ M @ P
    G = np.linalg.solve(Sp, M)
    w, V = np.linalg.eig(G)
    tiny = 1e-7 * np.abs(w).max()
    real = (np.abs(w) > tiny) & (w.imag == 0.0)          # a real eigenvalue of a real matrix comes back with imag 0
    # the sign of v^T M v is the block's sign only at a simple real eigenvalue (at a defective one v is neutral),
    # so the real eigenvalues must be separated and their vectors far from neutral
    wr = np.sort(w[real].real)
    simple = bool(len(wr) < 2 or np.min(np.diff(wr) / np.maximum(np.abs(wr[1:]), np.abs(wr[:-1]))) > 1e-6)
    q = np.array([V[:, k].real @ M @ V[:, k].real / (V[:, k].real @ V[:, k].real) for k in np.where(real)[0]])
    simple = simple and bool(np.all(np.abs(q) > 1e-6))
    # a real double root can come back as a complex pair split by rounding, which imag == 0 would miss: no nonzero
    # eigenvalue may sit closer to the axis than 1e-5 of its modulus without being exactly on it
    rel = np.abs(w.imag) / np.maximum(np.abs(w), tiny)
    simple = simple and not bool(np.any((np.abs(w) > tiny) & (w.imag != 0.0) & (rel < 1e-5)))
    eps_sum = int(np.sum(np.sign(q)))
    side = (int(np.sum((w[real].real > 0) & (q > 0))), int(np.sum((w[real].real > 0) & (q < 0))),
            int(np.sum((w[real].real < 0) & (q > 0))), int(np.sum((w[real].real < 0) & (q < 0))))
    _, Z, sdim = sl.schur(G, output="real", sort=lambda x, y: abs(x + 1j * y) < tiny)

    def sig(X):
        if X.size == 0:
            return 0
        ev = np.linalg.eigvalsh((X + X.T) / 2); c = 1e-9 * max(1.0, np.abs(ev).max())
        return int((ev > c).sum() - (ev < -c).sum())
    Z0 = Z[:, :sdim]
    return int(real.sum()), eps_sum, sig(Z0.T @ M @ Z0), sig(M), simple, side, wr


def half_couplings(Sp, M):
    """The couplings of the X-odd half without the elimination, so also where S is singular: with V scaled by i the
    half's T(g) is the real pencil [[S, 0], [0, -t S t]] - g [[0, B], [-t B t, 0]], B = t M from |O><E| to |E><O|
    and S = t S' t. Returns the number of finite real g > 0, multiplicity counted, and the gap that separates them
    from the rest: QZ returns a simple real root with imaginary part exactly 0, while a real double root may come
    back split by rounding, so roots within 1e-6 of the axis (relative) count as real, and the gap is the largest
    relative distance among those against the smallest among the others."""
    import scipy.linalg as sl
    k = Sp.shape[0]
    t = np.zeros((k, k)); m = int(round(np.sqrt(k)))
    for c in range(k):
        t[(c // m) + (c % m) * m, c] = 1
    S, B = t @ Sp @ t, t @ M
    Z = np.zeros((k, k))
    Ar = np.block([[S, Z], [Z, -Sp]]); Br = np.block([[Z, B], [-(t @ B @ t), Z]])
    g = sl.eigvals(Ar, Br)
    g = g[np.isfinite(g)]
    g = g[np.abs(g) > 1e-7 * max(1.0, np.abs(g).max())]
    rel = np.abs(g.imag) / np.abs(g)
    on = rel < 1e-6
    gap = (float(rel[on].max()) if on.any() else 0.0, float(rel[~on].min()) if (~on).any() else 1.0)
    return int(np.sum(on & (g.real > 0))), gap


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
    if "--real" in sys.argv:
        print("D4  the roots at real coupling (floating point)", flush=True)
        (nc, n0, nnz, real, mg, zg, margin), H, ham = real_couplings()
        check(f"two borders share {nc} roots (largest shared distance {mg[0]:.1e}, smallest unshared {mg[1]:.1e}), "
              f"{n0} at 0 (largest {zg[0]:.1e}, next {zg[1]:.1e}) and {nnz} elsewhere; expected 298, 34, 264",
              (nc, n0, nnz) == (298, 34, 264) and mg[1] > 1e6 * mg[0] and zg[1] > 1e6 * zg[0])
        pos = real[real > 0]
        check(f"{len(real)} of them real (imaginary part exactly 0), the others at least {margin:.2f} of their modulus "
              "off the axis; symmetric under g -> -g: " + ", ".join(f"{x:.6f}" for x in pos),
              len(real) == 16 and np.allclose(np.sort(-real), real, rtol=1e-9) and margin > 0.1)
        n = H.shape[0]; Id = np.eye(n)
        A = np.kron(H, Id) - np.kron(Id, H.T)
        bd = np.diag((8 - 2 * ham).ravel().astype(float))
        for g in pos:
            s0 = np.linalg.svd(-1j * A + g * bd, compute_uv=False)[-5:][::-1]
            # a fifth eigenvalue joining in a Jordan block would split at rounding by about sqrt(eps |T|) ~ 1e-7,
            # so exactly four eigenvalues below 1e-6 at g makes the four kernel vectors the whole eigenvalue 0
            at0 = int(np.sum(np.abs(np.linalg.eigvals(-1j * A + g * bd)) < 1e-6))
            pairs = []
            for t in (1e-3, 1e-4):
                w = np.linalg.eigvals(-1j * A + (1 + t) * g * bd)
                o = np.argsort(np.abs(w))
                pairs.append(w[o[2:4]])
            ratio = np.mean(np.abs(pairs[0])) / np.mean(np.abs(pairs[1]))
            pr = pairs[1]
            kind = ("oscillating pair on the line" if np.all(np.abs(pr.real) < 1e-6 * np.abs(pr)) else
                    "real pair crossing the line" if np.all(np.abs(pr.imag) < 1e-6 * np.abs(pr)) else "other")
            check(f"g = {g:.6f}: singular values 4th {s0[3]:.1e}, 5th {s0[4]:.1e}; {at0} eigenvalues below 1e-6; "
                  f"the pair beside the generic two grows linearly off g (ratio {ratio:.3f} for offsets 1e-3 : 1e-4); {kind}",
                  s0[4] > 1e6 * s0[3] and at0 == 4 and abs(ratio / 10 - 1) < 0.05 and kind != "other")
    if "--odd" in sys.argv:
        print("D5  the X-odd half: a Green's matrix at every even N, read at N = 8", flush=True)
        # N = 4: the Proposition's single point gamma = J, twice
        A1, B, t, _, _, _ = odd_half(4)
        w4 = np.linalg.eigvals(t @ np.linalg.solve(A1, B))
        r4 = w4[(w4.imag == 0) & (np.abs(w4) > 1e-9)].real
        check(f"N = 4: real nonzero eigenvalues of G give gamma = 1/|omega| = {np.round(1 / np.abs(r4), 12)}, expected [1, 1]",
              len(r4) == 2 and np.allclose(1 / np.abs(r4), 1.0, rtol=1e-12))
        # at gamma = J two pairs meet at once: nullity two, and off J one pair is real on one side, the other on the other
        m4 = A1.shape[0]; Z4 = np.zeros((m4, m4))
        Ao4 = np.block([[A1, Z4], [Z4, -(t @ A1 @ t)]]); Do4 = np.block([[Z4, B], [t @ B @ t, Z4]])
        s4 = np.linalg.svd(-1j * Ao4 + Do4, compute_uv=False)
        mu2 = []
        for d in (-1e-6, 1e-6):
            w = np.linalg.eigvals(-1j * Ao4 + (1 + d) * Do4); w = w[np.argsort(np.abs(w))[:4]]
            v = (w ** 2).real
            mu2.append(v[np.argsort(np.abs(v))][[0, 2]])         # one mu^2 per pair, the pairs told apart by size
        check(f"N = 4 at gamma = J: kernel two-dimensional (singular values {s4[-1]:.0e}, {s4[-2]:.0e}, next {s4[-3]:.2f}); "
              f"the two pairs' mu^2 are {mu2[0]} below J and {mu2[1]} above: each changes sign, in opposite directions",
              s4[-3] > 1e12 * s4[-2] and np.all(mu2[0] * mu2[1] < 0) and mu2[0][0] * mu2[0][1] < 0)
        for p in (1000003, 1000033):
            A1p, Bp, tp, _, _, _ = odd_half(8, p)
            Gp = tp @ (inv_mod(A1p, p) @ Bp % p) % p
            G2 = Gp @ Gp % p
            rB, r2, r3 = rank(Bp, p), rank(G2, p), rank(G2 @ Gp % p, p)
            check(f"N = 8 mod {p}: rank of D + 8 from O><E to E><O {rB}, G has {r2} = {r3} nonzero eigenvalues "
                  "(the stable rank of its powers, multiplicity counted); expected 595 and 566", (rB, r2, r3) == (595, 566, 566))
        A1, B, t, _, PE, PO = odd_half(8)
        G = t @ np.linalg.solve(A1, B)
        w, vec = np.linalg.eig(G)
        a = np.sort(np.abs(w)); eps = np.finfo(float).eps
        k0 = int(np.argmax(np.log(np.maximum(a[1:], eps)) - np.log(np.maximum(a[:-1], eps)))) + 1
        nz = np.abs(w) >= a[k0]
        real = nz & (w.imag == 0.0)                      # a real eigenvalue of a real matrix comes back with imag exactly 0
        gam = np.sort(1 / np.abs(w[real].real))
        Rop = np.kron(PE, PO)                            # the chain reflection on O><E coordinates
        par = np.array([np.real(np.vdot(vec[:, k], Rop @ vec[:, k]) / np.vdot(vec[:, k], vec[:, k]))
                        for k in np.where(real)[0]])
        offr = np.abs(w[nz & ~real].imag) / np.abs(w[nz & ~real])
        check(f"N = 8 floating point: {int(nz.sum())} nonzero eigenvalues (gap {a[k0 - 1]:.1e} to {a[k0]:.1e}), "
              f"{int(real.sum())} real, giving {len(np.unique(np.round(gam, 9)))} couplings; the reflection "
              f"(exact commutation below) splits them {int((par > 0).sum())} even, {int((par < 0).sum())} odd; "
              "expected 566, 48, 48, 24 + 24",
              int(nz.sum()) == 566 and int(real.sum()) == 48 and len(np.unique(np.round(gam, 9))) == 48
              and a[k0] > 1e6 * a[k0 - 1] and np.all(np.abs(np.abs(par) - 1) < 1e-6)
              and int((par > 0).sum()) == 24)
        print("    couplings: " + ", ".join(f"{g:.6f}" for g in gam), flush=True)
        _init_odd()
        Aodd, Dodd = _odd_svd.mats
        # a second route to the half itself: the cell-basis generator of the block (4,4) in the unnormalized
        # integer bases; every entry and partial sum is an integer far below 2^53, so the float products are exact
        A1i, Bi, ti, (Hf, ef, of, bcf), PEi, POi = odd_half(8, 1 << 40)
        Hf, ef, of, bcf, A1i, Bi, ti = (x.astype(float) for x in (Hf, ef, of, bcf, A1i, Bi, ti))
        A1i[A1i > (1 << 39)] -= 1 << 40; Bi[Bi > (1 << 39)] -= 1 << 40           # undo the modulus's sign shift
        mq = ef.shape[1]; Zq = np.zeros_like(A1i)
        Q = np.array([np.outer(*((ef[:, i], of[:, j]) if s == 0 else (of[:, i], ef[:, j]))).ravel()
                      for s in (0, 1) for j in range(mq) for i in range(mq)]).T     # same column-major order as U, V
        cellA = np.array([(Hf @ q.reshape(Hf.shape) - q.reshape(Hf.shape) @ Hf).ravel() for q in Q.T]).T
        cellD = bcf.ravel()[:, None] * Q
        # Q^T Q = 4, the integer h's are twice the normalized ones and the integer B four times, so exactly:
        okA = np.array_equal(Q.T @ cellA, 2 * np.block([[A1i, Zq], [Zq, -(ti @ A1i @ ti)]]))
        okD = np.array_equal(Q.T @ cellD, np.block([[Zq, Bi], [ti @ Bi @ ti, Zq]]))
        inv = np.array_equal(4 * cellA, Q @ (Q.T @ cellA)) and np.array_equal(4 * cellD, Q @ (Q.T @ cellD))
        check("the half built from the cells: [H, .] and D + 8 projected onto it equal the construction and leave "
              "it invariant, exactly in the integer bases", okA and okD and inv)
        PEi, POi = PEi.astype(float), POi.astype(float)
        R_eo, R_oe = np.kron(POi, PEi), np.kron(PEi, POi)                  # the reflection on E><O and O><E
        refl = (np.array_equal(A1i @ R_eo, R_eo @ A1i) and np.array_equal(Bi @ R_oe, R_eo @ Bi)
                and np.array_equal(ti @ R_eo, R_oe @ ti))
        check("the chain reflection commutes with S, D + 8 and t, hence with G, exactly in the integer bases", refl)
        x_even = (0.099186, 0.282490, 0.624282, 0.784448, 0.873035, 1.334206, 1.518087, 6.551270)   # D4
        sep = min(abs(g / x - 1) for g in gam for x in x_even)
        check(f"none of the 48 among the X-even half's eight (nearest relative distance {sep:.1e})", sep > 1e-4)
        near = w[nz & ~real & (w.imag > 0)]; near = near[np.abs(near.imag) / np.abs(near) < 1e-2]
        print(f"    {len(near)} non-real eigenvalue(s) within 1e-2 of the axis to examine", flush=True)
        from scipy.optimize import minimize_scalar
        for z in near:
            g0, rel = 1 / abs(z.real), abs(z.imag) / abs(z)
            f = lambda g: np.linalg.svd(-1j * Aodd + g * Dodd, compute_uv=False)[-1]
            res = minimize_scalar(f, bounds=(g0 * (1 - 5 * rel), g0 * (1 + 5 * rel)), method="bounded",
                                  options={"xatol": 1e-12})
            nT = np.linalg.norm(-1j * Aodd + res.x * Dodd, 2)
            check(f"the non-real eigenvalue nearest the axis ({rel:.1e} of its modulus off it, at 1/|Re| = {g0:.6f}): "
                  f"the smallest singular value found along real g is {res.fun:.1e}, {res.fun / (eps * nT):.0e} "
                  "times eps |T|: the pair is non-real, not a real double split by rounding",
                  res.fun > 1e6 * eps * nT)
        away = [1.0 if (k == len(gam) - 1 or (k > 0 and gam[k] - gam[k - 1] < gam[k + 1] - gam[k])) else -1.0
                for k in range(len(gam))]
        rows = [_odd_svd(a) for a in zip(gam, away)]
        bad = []; nreal = 0; worst_expo = 0.0
        for g, (s0, s1, nT, lam, (pa, pb, pc)) in rows:
            # nullity one: a backward-stable SVD reads an exactly singular matrix at a few eps |T|, the next six
            # decades above. Algebraic multiplicity two: a 2 x 2 Jordan block at 0 splits at rounding by about
            # sqrt(eps |T|), so two eigenvalues below 10 sqrt(eps |T|) yet far above the eps |T| at which a semisimple
            # double zero would split, and the third three decades above.
            # Off g the pair is +-mu with mu^2 linear in the offset: |mu| at 1e-5 : 1e-6 is sqrt 10,
            # and mu^2 = -(product of the pair) changes sign across g (a real pair on one side, oscillating on the other).
            r = np.sqrt(eps * nT)
            m2a, m2c = -(pa[0] * pa[1]).real, -(pc[0] * pc[1]).real
            grow = np.mean(np.abs(pa)) / np.mean(np.abs(pb))
            # The Jordan block itself follows from nullity one and the even multiplicity sigma forces; the growth law
            # is read as a classification: the exponent log10(grow) is 1/2 at an exceptional point and 1 at a crossing,
            # so the band (0.45, 0.55) separates them (second order in the offset and rounding move it by at most
            # 0.007 over the 48)
            expo = np.log10(grow)
            ok = (s0 < 10 * eps * nT and s1 > 1e6 * s0 and 1e3 * eps * nT < lam[0] and lam[1] < 10 * r
                  and lam[2] > 1e3 * lam[1] and abs(expo - 0.5) < 0.05 and m2a * m2c < 0)
            worst_expo = max(worst_expo, abs(expo - 0.5))
            nreal += m2a > 0
            if not ok:
                bad.append((g, s0 / (eps * nT), lam / r, grow, m2a, m2c))
        check("every coupling is an exceptional point of the X-odd half: nullity one (smallest singular value at most "
              "10 eps |T|, the next six decades above), algebraic multiplicity two (two eigenvalues between "
              "1e3 eps |T| and 10 sqrt(eps |T|), the third three decades above), the pair +-mu growing as the square root of the "
              f"offset (exponent within {worst_expo:.3f} of 1/2; a crossing gives 1) and mu^2 changing sign across g "
              f"({nreal} of 48 with the real pair on the side away from "
              "the neighbour)" + (f"; failing: {bad}" if bad else ""), not bad)
    if "--krein" in sys.argv:
        print("D6  the X-odd half's real couplings from below: the transpose's two forms", flush=True)
        from math import comb
        chain = lambda N: [(i, i + 1, 1.0, 1.0) for i in range(N - 1)]

        def inertia(X):
            ev = np.linalg.eigvalsh(X); c = 1e-9 * max(1.0, np.abs(ev).max())
            return int((ev > c).sum()), int((ev < -c).sum()), int((np.abs(ev) <= c).sum())
        for N in (4, 6, 8):
            m = comb(N, N // 2) // 2
            # the integer half of odd_half (unnormalized bases, so B is four times the normalized one) against the
            # pair form: t B sends V_ij to beta_ij at (j, i), nothing else, exactly
            _, Bi, ti, _, _, _ = odd_half(N, 1 << 40)
            Bi = Bi.astype(np.int64); Bi[Bi > (1 << 39)] -= 1 << 40
            Sp, M, w = krein_half(N, chain(N))
            Mi = np.zeros_like(Bi)
            for i in range(m):
                for j in range(m):
                    Mi[i + j * m, j + i * m] = 4 * w[i, j]
            exact = np.array_equal(ti.astype(np.int64) @ Bi, Mi)
            # the inertia the pair form predicts (m diagonal squares beta_ii = N > 0, one positive and one negative
            # square per pair i < j with beta_ij != 0, two zeros per pair with beta_ij = 0) against M's own spectrum
            off = [w[i, j] for i in range(m) for j in range(i + 1, m)]
            nz = sum(1 for x in off if x != 0)
            pred = (m + nz, nz, 2 * (len(off) - nz))
            check(f"N = {N}: t(D + N) on the half is the pair form exactly; its spectrum has inertia {inertia(M)}, "
                  f"the pair form's count {pred}: signature {m} = C({N}, {N // 2})/2",
                  exact and inertia(M) == pred)
        # D5's couplings by this construction, and by the pencil without the elimination
        A1, B, t, _, _, _ = odd_half(8)
        w5 = np.linalg.eigvals(t @ np.linalg.solve(A1, B))
        g5 = np.sort(1 / np.abs(w5[(w5.imag == 0.0) & (np.abs(w5) > 1e-6)].real))
        Sp, M, _ = krein_half(8, chain(8))
        nr, es, s0, sM, simple, side, wr = krein_count(Sp, M)
        g6 = np.sort(1 / np.abs(wr))
        hc, hgap = half_couplings(Sp, M)
        check(f"N = 8 chain: the same {len(g6)} couplings as D5's construction (largest relative difference "
              f"{np.max(np.abs(g6 / g5 - 1)):.1e}); the pencil without the elimination finds {hc} (on the axis within "
              f"{hgap[0]:.0e}, the rest at least {hgap[1]:.0e} off it)",
              len(g5) == len(g6) == 48 and np.allclose(g5, g6, rtol=1e-9) and hc == 48 and hgap[0] == 0.0 and hgap[1] > 1e-5)
        # sign characteristics: on the real eigenvalues of G = S'^-1 M they and M on the zero root space add up to
        # sig(M), complex pairs carrying none, so #real >= sig(M) - s0; read per reflection sector at N = 8
        for N, expect in ((4, (2, 1)), (6, (14, 0)), (8, (48, 1))):
            m = comb(N, N // 2) // 2
            Sp, M, _ = krein_half(N, chain(N))
            nr, es, s0, sM, simple, side, _ = krein_count(Sp, M)
            check(f"N = {N} chain: {nr} real nonzero omega, simple, their M-signs summing to {es} (omega > 0: {side[0]} "
                  f"positive, {side[1]} negative; omega < 0: {side[2]} positive, {side[3]} negative), M on the zero root "
                  f"space {s0}, total {es + s0} = sig M = {m}; forced {m - s0}, expected (real, s0) = {expect}",
                  simple and sM == m and es + s0 == m and (nr, s0) == expect and nr >= m - s0
                  and side == {4: (1, 0, 1, 0), 6: (6, 2, 6, 0), 8: (20, 4, 21, 3)}[N])
        _, _, _, _, PE, PO = odd_half(8)
        rv, rV = np.linalg.eigh(np.kron(PE, PO))
        Sp, M, _ = krein_half(8, chain(8))
        rows = [krein_count(Sp, M, rV[:, np.isclose(rv, s)]) for s in (1, -1)]
        check("N = 8 chain per reflection sector (+, -): real " + ", ".join(str(r[0]) for r in rows) + "; signature of M "
              + ", ".join(str(r[3]) for r in rows) + "; M on the zero root space " + ", ".join(str(r[2]) for r in rows)
              + "; forced " + ", ".join(str(r[3] - r[2]) for r in rows) + "; expected 24, 24; 15, 20; 1, 0; 14, 20",
              all(r[4] and r[1] + r[2] == r[3] for r in rows) and [r[0] for r in rows] == [24, 24]
              and [r[3] for r in rows] == [15, 20] and [r[2] for r in rows] == [1, 0])
        # N = 6: no rung cell (h is even, N - 2h = 2 mod 4), M invertible, s0 = 0: at least 10 wherever S is invertible
        rng = np.random.default_rng(7)
        counts = []

        def connected(E):
            seen, todo = {0}, [0]
            while todo:
                a = todo.pop()
                for i, j in E:
                    for x, y in ((i, j), (j, i)):
                        if x == a and y not in seen:
                            seen.add(y); todo.append(y)
            return len(seen) == 6
        while len(counts) < 40:
            E = [(i, j) for i in range(6) for j in range(i + 1, 6) if rng.random() < 0.5]
            if not connected(E):
                continue
            Sp, M, _ = krein_half(6, [(i, j, rng.uniform(0.3, 2), rng.uniform(-2, 2)) for i, j in E])
            if np.linalg.svd(Sp, compute_uv=False)[-1] < 1e-6:
                continue
            nr, es, s0, sM, simple, _, _ = krein_count(Sp, M)
            hc, hgap = half_couplings(Sp, M)
            counts.append((nr, simple and es + s0 == sM == 10 and s0 == 0, hc, hgap))
        check(f"N = 6, 40 random connected XXZ graphs with S invertible: real counts from {min(c[0] for c in counts)} to "
              f"{max(c[0] for c in counts)}, never below 10, the sign identity holding with s0 = 0 and the pencil without "
              f"the elimination agreeing in each (on the axis within {max(c[3][0] for c in counts):.0e}, the rest at least "
              f"{min(c[3][1] for c in counts):.0e} off it)",
              all(c[0] >= 10 and c[1] and c[2] == c[0] and c[3][0] == 0.0 and c[3][1] > 1e-5 for c in counts))
        # where S is singular the bound has no hold, and the pencil reads what is there: the chain at Delta = 1/2
        # (h_E and h_O share an eigenvalue) stays above 10, five chained sites beside an isolated sixth fall below
        rows = []
        for name, b in (("the chain at Delta = 1/2", [(i, i + 1, 1.0, 0.5) for i in range(5)]),
                        ("a five-site chain beside an isolated site", [(i, i + 1, 1.0, 1.0) for i in range(4)])):
            Sp, M, _ = krein_half(6, b)
            hc, hgap = half_couplings(Sp, M)
            rows.append((name, np.linalg.svd(Sp, compute_uv=False)[-1], hc, hgap))
        check("N = 6 with S singular, read by the pencil: " + "; ".join(
              f"{r[0]}: smallest singular value of S {r[1]:.0e}, {r[2]} couplings with multiplicity" for r in rows)
              + "; expected 15 and 8",
              all(r[1] < 1e-12 and r[3][0] == 0.0 and r[3][1] > 1e-5 for r in rows) and [r[2] for r in rows] == [15, 8])
    print("\nALL GATES PASS" if not FAILS else f"\n{len(FAILS)} FAILED")
    sys.exit(1 if FAILS else 0)
