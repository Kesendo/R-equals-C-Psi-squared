#!/usr/bin/env python3
"""The qudit product-mirror cap: the operator realization of F121 (2026-06-11).

F121 (PROOF_QUDIT_PARTIAL_PALINDROME) counts the partial palindrome of the
full-Cartan dephasing dissipator at d > 2: paired ceiling
Sigma_k d^N C(N,k) (d-1)^min(k, N-k), full iff d = 2. This verifier banks the
OPERATOR side of that count:

(1) PRODUCT CAP (theorem): any per-site mirror W = tensor_l q_l (one-sided,
    two-sided, antilinear, site-dependent) that intertwines the dissipator
    palindrome W L_D = (-L_D - 2N*gamma) W on its support pairs at most
    (2d)^N of the d^(2N) coherences. Proof: rate additivity forces each q_l
    to be a strict per-site class swap (dark {i=j, d dims} <-> lit {i != j,
    d^2-d dims}); rank(q_l) <= min(d, d^2-d) + min(d, d^2-d) = 2d for d >= 2.
    The cap is FULL iff (2d)^N = d^(2N) iff d^2 - 2d = 0 iff d = 2: the
    QUBIT_NECESSITY trunk polynomial, third appearance (per-site split,
    ceiling column, operator cap).
(2) THE OPERATOR: the qubit palindromizer's formula generalizes VERBATIM,
    Pi_d(rho) = rho^T * Shift^(tensor N) (F118: Pi_Z = rho^T * X^(tensor N)).
    It attains the cap exactly: on the shift-aligned subspace (per-site
    letters {(x,x)} U {(a, a-1)}, dimension (2d)^N, Pi_d-closed) the
    intertwining residual is EXACTLY ZERO; on the complement it fails at
    O(gamma) (the provably unpaired part). Two chiralities Pi_d^+/- (the two
    shift directions); at d = 2 the two off-diagonals coincide, the
    chiralities merge, and the mirror is full: that degeneracy IS the qubit
    magic.
(3) THE MIRROR GROUP: <Pi_d, D> (D = transpose) has order 2d^2 with
    ord(Pi_d) = 2d, and D conjugation EXCHANGES the two shift factors:
    <Pi_d, D> ~ Z_d wr Z_2 (wreath product). At d = 2 this is D_4: the F118
    mirror group is the d = 2 column of a d-indexed family. D does NOT
    preserve the aligned subspace for d > 2 (it swaps the chiralities).
(4) HONESTY: the combinatorial ceiling IS reachable by a global
    (non-product) partial isometry (explicit greedy rung matching, exact
    intertwining on its support); the gap ceiling - (2d)^N (= 18 at d=3,
    N=2; first nonzero at N = 2) is therefore exactly the NON-PRODUCT part
    of the partial palindrome.

Blocks: A cap + trunk arithmetic (d = 2..5); B Pi_d exact on the aligned
subspace (d = 3, N = 1..3; d = 4, N = 1..2); C random per-site class swaps
never exceed the cap (rank check); D the global ceiling-reacher; E the
group law |<Pi_d, D>| = 2d^2 and ord(Pi_d) = 2d (d = 2..5); F d = 2
degeneracy (aligned subspace = everything, Pi_2 = the F118 palindromizer).

All assertions exact (permutation/integer arithmetic; float only in eig-free
residual norms). Companion: PROOF_QUDIT_PARTIAL_PALINDROME.md section 6.
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


def block_a():
    print("BLOCK A  cap + trunk arithmetic")
    for d in (2, 3, 4, 5):
        for N in (1, 2, 3):
            cap = (2 * d) ** N
            ceil = ceiling(d, N)
            total = d ** (2 * N)
            assert cap <= ceil <= total
            assert (cap == total) == (d == 2), (d, N)
            assert (cap == ceil) == (d == 2 or N == 1), (d, N)
        print(f"  d={d}: cap (2d)^N {'= full iff d=2 OK' if d == 2 else f'< ceiling for N>=2 OK'}"
              f"  (N=2: cap {(2*d)**2}, ceiling {ceiling(d,2)}, total {d**4})")
    print("BLOCK A PASS")


def block_b():
    print("BLOCK B  Pi_d exact on the shift-aligned subspace")
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
    print("BLOCK C  random per-site class swaps never exceed the cap (rank)")
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
    print(f"  5 random class-swap products at (3,2): rank = 36 = cap, never above")
    print("BLOCK C PASS")


def block_d():
    print("BLOCK D  the global (non-product) ceiling-reacher")
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
          f"the gap 54 - 36 = 18 is exactly the NON-PRODUCT part")
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


def main():
    print("=" * 78)
    print("THE QUDIT PRODUCT-MIRROR CAP (operator side of F121)")
    print("=" * 78)
    block_a()
    block_b()
    block_c()
    block_d()
    block_e()
    block_f()
    print("=" * 78)
    print("ALL BLOCKS PASS")
    print("=" * 78)


if __name__ == "__main__":
    main()
