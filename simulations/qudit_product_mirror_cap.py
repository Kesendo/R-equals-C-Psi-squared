#!/usr/bin/env python3
"""The qudit product-mirror cap: the operator realization of F121 (2026-06-11).

F121 (PROOF_QUDIT_PARTIAL_PALINDROME) counts the partial palindrome of the
full-Cartan dephasing dissipator at d > 2: paired ceiling
Sigma_k d^N C(N,k) (d-1)^min(k, N-k), full iff d = 2. This verifier banks the
OPERATOR side of that count:

(1) THE PRODUCT CAP (theorem): any per-site mirror W = tensor_l q_l that
    intertwines W L_D = (-L_D - 2N gamma) W carries ONE lit grade
    c_l = (lit out) + (lit in) in {0, 1, 2} per site (rate additivity asks
    sum c_l = N on every block product; varying one site pins its grade),
    with rank <= d, 2d, d^2 - d. sum c_l = N forces #(c=0) = #(c=2) = m, so
        rank W <= P(d, N) = max_m (2d)^(N-2m) * (d^3 - d^2)^m.
    A (0,2) pair beats a (1,1) pair iff d^3 - d^2 > 4d^2 iff d >= 6 (tie at
    d = 5): P = (2d)^N for every N when d <= 5 (the strict dark <-> lit swap),
    larger from d = 6 (P_dark tensor P_lit at d=6, N=2: rank 180 > 144).
    Full iff P = d^(2N) iff d^2 - 2d = 0 iff d = 2: the trunk's third
    appearance.
(2) THE OPERATOR: the qubit palindromizer's formula generalizes VERBATIM,
    Pi_d(rho) = rho^T * Shift^(tensor N) (F118: Pi_Z = rho^T * X^(tensor N)).
    It attains the cap for d <= 5: on the shift-aligned subspace (per-site
    letters {(x,x)} U {(a, a-1)}, dimension (2d)^N, Pi_d-closed) the
    intertwining residual is EXACTLY ZERO; on the complement this realization
    fails at O(gamma). Two chiralities Pi_d^+/- (the two
    shift directions); at d = 2 the two off-diagonals coincide, the
    chiralities merge, and the mirror is full: that degeneracy IS the qubit
    magic.
(3) THE MIRROR GROUP: <Pi_d, D> (D = transpose) has order 2d^2 with
    ord(Pi_d) = 2d, and D conjugation EXCHANGES the two shift factors:
    <Pi_d, D> ~ Z_d wr Z_2 (wreath product). At d = 2 this is D_4: the F118
    mirror group is the d = 2 column of a d-indexed family. D does NOT
    preserve the aligned subspace for d > 2 (it swaps the chiralities).
(4) GLOBAL REACHER: the combinatorial ceiling IS reachable by a global
    (non-product) partial isometry (explicit greedy rung matching, exact
    intertwining on its support); since no product reaches past P, the gap
    ceiling - P (= 18 at d=3, N=2; first nonzero at N = 2) is exactly the
    non-product part of the partial palindrome.

Blocks: A cap arithmetic (d = 2..5); B Pi_d P_aligned exact on the aligned
subspace (d = 3, N = 1..3; d = 4, N = 1..2); C random per-site class swaps
attain (2d)^N at (3,2), and the d = 6 projector P_dark tensor P_lit; D the
global ceiling-reacher; E the group law |<Pi_d, D>| = 2d^2 and
ord(Pi_d) = 2d (d = 2..5); F d = 2 degeneracy (aligned subspace =
everything, Pi_2 = the F118 palindromizer); G the product lemma from below
(integer grade factors reach rank d, 2d, d^2 - d with residual exactly 0, a
mixed-grade factor breaks the identity, and P(d, N) equals the enumeration
of every grade pattern at d = 2..8, N = 1..5).

All assertions exact (permutation/integer arithmetic; float only in eig-free
residual norms and in SVD ranks of 0/1-scale matrices whose rank the block
structure bounds). Companion: PROOF_QUDIT_PARTIAL_PALINDROME.md section 6.
"""
import sys
from itertools import product as iproduct
from math import comb
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

GAMMA = 1.0


def basis(d, N):
    states = list(iproduct(range(d), repeat=N))
    pairs = [(i, j) for i in states for j in states]
    return states, pairs


def L_diss_diag(d, N):
    _, pairs = basis(d, N)
    return np.array([-2.0 * GAMMA * sum(1 for a, b in zip(i, j) if a != b)
                     for (i, j) in pairs])


def pi_d_matrix(d, N, chirality=+1):
    """Pi_d: E_ij -> E_(j, i - chirality)  (rho -> rho^T * Shift^chirality)."""
    _, pairs = basis(d, N)
    idx = {p: n for n, p in enumerate(pairs)}
    n = len(pairs)
    W = np.zeros((n, n))
    for col, (i, j) in enumerate(pairs):
        tgt = (j, tuple((a - chirality) % d for a in i))
        W[idx[tgt], col] = 1.0
    return W


def transpose_matrix(d, N):
    _, pairs = basis(d, N)
    idx = {p: n for n, p in enumerate(pairs)}
    n = len(pairs)
    Dm = np.zeros((n, n))
    for col, (i, j) in enumerate(pairs):
        Dm[idx[(j, i)], col] = 1.0
    return Dm


def aligned_mask(d, N, chirality=+1):
    _, pairs = basis(d, N)
    def good(a, b):
        return a == b or b == (a - chirality) % d
    return np.array([all(good(a, b) for a, b in zip(i, j)) for (i, j) in pairs])


def ceiling(d, N):
    return sum(d ** N * comb(N, k) * (d - 1) ** min(k, N - k) for k in range(N + 1))


def product_cap(d, N):
    """P(d, N) = max_m (2d)^(N-2m) * (d^3 - d^2)^m, the product-mirror cap."""
    return max((2 * d) ** (N - 2 * m) * (d ** 3 - d ** 2) ** m for m in range(N // 2 + 1))


def block_a():
    print("BLOCK A  cap arithmetic: (2d)^N = P <= ceiling <= d^(2N) for d <= 5")
    for d in (2, 3, 4, 5):
        for N in (1, 2, 3):
            cap = product_cap(d, N)
            ceil = ceiling(d, N)
            total = d ** (2 * N)
            assert cap == (2 * d) ** N
            assert cap <= ceil <= total
            assert (cap == total) == (d == 2), (d, N)
            assert (cap == ceil) == (d == 2 or N == 1), (d, N)
        print(f"  d={d}: cap (2d)^N {'= full at d=2 OK' if d == 2 else f'< ceiling for N>=2'}"
              f"  (N=2: cap {(2*d)**2}, ceiling {ceiling(d,2)}, total {d**4})")
    print("BLOCK A PASS")


def block_b():
    print("BLOCK B  Pi_d P_aligned exact on the shift-aligned subspace")
    for d, Ns in ((3, (1, 2, 3)), (4, (1, 2))):
        for N in Ns:
            diag = L_diss_diag(d, N)
            W = pi_d_matrix(d, N)
            L = np.diag(diag)
            res = W @ L + L @ W + 2 * N * GAMMA * W
            mask = aligned_mask(d, N)
            sub = np.where(mask)[0]
            sub_set = set(sub)
            assert len(sub) == (2 * d) ** N
            closed = all(int(np.argmax(np.abs(W[:, c]))) in sub_set for c in sub)
            assert closed, f"aligned subspace not Pi_d-closed at d={d}, N={N}"
            r_sub = np.max(np.abs(res[:, sub]))
            assert r_sub == 0.0, f"nonzero residual on aligned: {r_sub}"
            comp = np.where(~mask)[0]
            r_comp = np.max(np.abs(res[:, comp]))
            assert r_comp > 0.5, "complement unexpectedly mirrored"
            print(f"  d={d} N={N}: aligned dim {(2*d)**N}, residual 0.0 EXACT on aligned, "
                  f"{r_comp:.1f} on complement")
    print("BLOCK B PASS")


def block_c():
    print("BLOCK C  strict class swaps attain (2d)^N; from d = 6 the (0,2) pairing beats them")
    rng = np.random.default_rng(7)
    d, N = 3, 2
    lit_count = d * d - d
    for trial in range(5):
        sel = rng.permutation(lit_count)[:d]
        U1 = np.linalg.qr(rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d)))[0]
        U2 = np.linalg.qr(rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d)))[0]
        q = np.zeros((d * d, d * d), dtype=complex)
        for r_ in range(d):
            for c_ in range(d):
                q[d + sel[r_], c_] = U1[r_, c_]
                q[c_, d + sel[r_]] = U2[c_, r_]
        Wq = np.kron(q, q)
        rank = np.linalg.matrix_rank(Wq)
        assert rank <= (2 * d) ** N
        assert rank == (2 * d) ** N   # these constructions attain it
    print(f"  5 random strict class-swap products at (3,2): rank = 36 = cap")
    d6, N6 = 6, 2
    diag6 = L_diss_diag(d6, N6)
    _, pairs6 = basis(d6, N6)
    support6 = np.array([i[0] == j[0] and i[1] != j[1] for i, j in pairs6])
    d6_rank = int(np.sum(support6))
    # W=P_dark tensor P_lit is diagonal on this support, so the residual diagonal is
    # W L_D + L_D W + 2N gamma W = (2*diag+2N gamma)*support.
    residual6 = (2 * diag6 + 2 * N6 * GAMMA) * support6
    assert d6_rank == 180 == product_cap(d6, N6) and d6_rank > (2 * d6) ** N6
    assert np.max(np.abs(residual6)) == 0.0
    print("  d=6,N=2: P_dark tensor P_lit stays on h=1, residual 0, rank 180 = P > 144 = (2d)^N")
    print("  -> the swap is the optimum only up to d = 5; the cap is P(d, N)")
    print("BLOCK C PASS")


def block_d():
    print("BLOCK D  a global ceiling-reacher")
    d, N = 3, 2
    diag = L_diss_diag(d, N)
    n = len(diag)
    rungs = {}
    for m, r_ in enumerate(diag):
        rungs.setdefault(round(-r_ / (2 * GAMMA)), []).append(m)
    Wg = np.zeros((n, n))
    paired = 0
    for k in sorted(rungs):
        kk = N - k
        if k > kk:
            continue
        if k == kk:
            for x in rungs[k]:
                Wg[x, x] = 1.0
            paired += len(rungs[k])
        else:
            m = min(len(rungs[k]), len(rungs[kk]))
            for x, y in zip(rungs[k][:m], rungs[kk][:m]):
                Wg[y, x] = 1.0
                Wg[x, y] = 1.0
            paired += 2 * m
    assert paired == ceiling(d, N) == 54
    L = np.diag(diag)
    sup = np.abs(Wg).sum(axis=0) > 0
    res = np.max(np.abs((Wg @ L + L @ Wg + 2 * N * GAMMA * Wg)[:, sup]))
    assert res == 0.0
    print(f"  (3,2): paired = 54 = ceiling, exact intertwining on support; "
          f"the non-product part is 54 - 36 = 18")
    print("BLOCK D PASS")


def block_e():
    print("BLOCK E  the mirror group law: |<Pi_d, D>| = 2d^2, ord(Pi_d) = 2d")
    for d in (2, 3, 4, 5):
        N = 1
        W = pi_d_matrix(d, N)
        Dm = transpose_matrix(d, N)
        n = W.shape[0]
        def order_of(M, maxo=4 * d + 2):
            P = np.eye(n)
            for o in range(1, maxo + 1):
                P = P @ M
                if np.array_equal(P, np.eye(n)):
                    return o
            return None
        assert order_of(W) == 2 * d
        assert order_of(Dm) == 2
        elems = [np.eye(n)]
        keys = {elems[0].tobytes()}
        frontier = [np.eye(n)]
        while frontier:
            new = []
            for E in frontier:
                for g in (W, Dm):
                    C = g @ E
                    kbytes = C.tobytes()
                    if kbytes not in keys:
                        keys.add(kbytes)
                        elems.append(C)
                        new.append(C)
            frontier = new
        assert len(elems) == 2 * d * d, (d, len(elems))
        # D exchanges the two chiralities' aligned subspaces (d > 2)
        if d > 2:
            mp = aligned_mask(d, N, +1)
            mm = aligned_mask(d, N, -1)
            perm = np.argmax(Dm, axis=0)
            img_plus = set(perm[np.where(mp)[0]])
            assert img_plus == set(np.where(mm)[0]), "D does not swap chiralities"
        print(f"  d={d}: |<Pi_d, D>| = {2*d*d} = 2d^2, ord(Pi_d) = {2*d} = 2d"
              + ("  (D swaps the +/- chiralities)" if d > 2 else "  (= D4, F118)"))
    print("  => <Pi_d, D> ~ Z_d wr Z_2; the F118 mirror group D4 is the d = 2 column")
    print("BLOCK E PASS")


def block_f():
    print("BLOCK F  d = 2 degeneracy: aligned = everything, the mirror is full")
    d = 2
    for N in (1, 2, 3):
        mask = aligned_mask(d, N)
        assert mask.all(), "d=2 aligned subspace is not the full space"
        diag = L_diss_diag(d, N)
        W = pi_d_matrix(d, N)
        L = np.diag(diag)
        res = np.max(np.abs(W @ L + L @ W + 2 * N * GAMMA * W))
        assert res == 0.0
    print("  d=2, N=1..3: the shift-aligned subspace IS the whole space and the")
    print("  intertwining is globally exact: Pi_2 = rho^T X^(tensor N), the F118 palindromizer")
    print("BLOCK F PASS")


def grade_factor(d, grades, rng):
    """Integer per-site factor on the d^2 letters (index a*d + b for |a><b|), nonzero exactly
    where the block's lit grade (lit out + lit in) lies in `grades`."""
    lit = np.array([(l // d) != (l % d) for l in range(d * d)], dtype=int)
    grade = lit[:, None] + lit[None, :]
    mask = np.isin(grade, list(grades))
    return mask * rng.integers(1, 10, size=(d * d, d * d)), grade


def rank_mod_p(M, p=1_000_000_007):
    A = [[int(x) % p for x in row] for row in M]
    n, m, r = len(A), len(A[0]), 0
    for c in range(m):
        piv = next((i for i in range(r, n) if A[i][c]), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = pow(A[r][c], p - 2, p)
        A[r] = [x * inv % p for x in A[r]]
        for i in range(n):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [(x - f * y) % p for x, y in zip(A[i], A[r])]
        r += 1
    return r


def block_g():
    print("BLOCK G  the product lemma from below")
    rng = np.random.default_rng(11)
    for d in (3, 4, 5, 6):
        want = {0: d, 1: 2 * d, 2: d * d - d}
        for c, r in want.items():
            q, _ = grade_factor(d, [c], rng)
            # rank mod p <= rank over Q <= the block bound, so equality at the bound is exact
            assert rank_mod_p(q) == r, (d, c)
        for pattern in ((1, 1), (0, 2)):
            q0, g = grade_factor(d, [pattern[0]], rng)
            q1, _ = grade_factor(d, [pattern[1]], rng)
            # residual of W = q0 (x) q1 at N = 2: entry W * 2 gamma (2 - Ham(out) - Ham(in))
            h = g[:, :, None, None] + g[None, None, :, :]      # Ham(out) + Ham(in) per product entry
            W = np.einsum("ab,cd->abcd", q0, q1)
            assert np.max(np.abs(W * (2 - h))) == 0, (d, pattern)
        qm, g = grade_factor(d, [0, 1], rng)
        q1, _ = grade_factor(d, [1], rng)
        hm = g[:, :, None, None] + g[None, None, :, :]
        Wm = np.einsum("ab,cd->abcd", qm, q1)
        assert np.max(np.abs(Wm * (2 - hm))) > 0, "a mixed-grade factor must break the identity"
        print(f"  d={d}: grade ranks {want[0]}, {want[1]}, {want[2]} (mod p, exact); "
              f"(1,1) and (0,2) residual 0; mixed factor breaks")
    for d in range(2, 9):
        for N in range(1, 6):
            best = 0
            for pattern in iproduct((0, 1, 2), repeat=N):
                if sum(pattern) == N:
                    r = 1
                    for c in pattern:
                        r *= (d, 2 * d, d * d - d)[c]
                    best = max(best, r)
            assert best == product_cap(d, N), (d, N)
            assert (product_cap(d, N) == (2 * d) ** N) == (d <= 5 or N == 1), (d, N)
            assert (product_cap(d, N) == d ** (2 * N)) == (d == 2), (d, N)
    for d in range(2, 31):
        assert (d ** 3 - d ** 2 > 4 * d * d) == (d >= 6)
    assert 5 ** 3 - 5 ** 2 == 4 * 25
    print("  P(d, N) = pattern maximum at d = 2..8, N = 1..5; swap optimal iff d <= 5 (tie at 5);"
          " full iff d = 2")
    print(f"  P(6,2) = {product_cap(6, 2)}, P(7,2) = {product_cap(7, 2)}, P(8,2) = {product_cap(8, 2)}")
    print("BLOCK G PASS")


def main():
    print("=" * 78)
    print("THE QUDIT PRODUCT-MIRROR CAP + THE SHIFT MIRROR (operator side of F121)")
    print("=" * 78)
    block_a()
    block_b()
    block_c()
    block_d()
    block_e()
    block_f()
    block_g()
    print("=" * 78)
    print("ALL BLOCKS PASS")
    print("=" * 78)


if __name__ == "__main__":
    main()
