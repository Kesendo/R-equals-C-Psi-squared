"""Gate for the bi-chromatic blocks of docs/proofs/PROOF_DEPHASING_FRONT_RENEWAL.md (section "The same split on
every bi-chromatic block").

Uniform Z dephasing at rate gamma per site, H number-conserving, L_H = -i[H, .]. The Liouvillian keeps the
joint-popcount blocks
(p, q), the cells |i><j| with popcount(i) = p and popcount(j) = q, and the dephasing charges a cell the rate
2 gamma Hamming(i, j); the Hamming distance is the cell's height, and the number of distinct heights in a block is
its chromaticity, the number of distinct dephasing rates in it (F74, PROOF_CHROMATICITY, PROOF_UNIFORM_LAW). The
settings: the Heisenberg chain, the XXZ chain at Delta = 2 with Z fields, a ring topology with unequal couplings,
the XY star, a seeded random graph with integer couplings, anisotropies and fields, a single hopping bond (XX + YY
on sites 0 and 1), a DM chain (X_i Y_j - Y_i X_j bonds, purely imaginary hopping, so a purely imaginary H) and the
same DM chain with a Z field (a complex H).

  T1  the chromaticity, a re-check of the owned count. Block (p, q) holds the Hamming distances p + q - 2k,
      max(0, p + q - N) <= k <= min(p, q), so min(p, q, N - p, N - q) + 1 of them: one on the 4N blocks of the
      grid's border (F153's pinned ones), two on the 4N - 8 blocks one step inside it (N >= 3; the block (1, 1)
      alone at N = 2). Brute force over every cell, N = 1..10; control: the count min(p, q) + 1 is refuted by
      the same brute force at every N.
  T2  the split. On every bi-chromatic block the generator, built on all 4^N cells from H and the dephasing's
      action on each cell (the signs of Z_l on ket and bra), closes on the block, and its dissipator equals
      -2 gamma h_+ + 4 gamma P cell by cell, h_+ the upper height and P the projector onto the cells at the lower
      one (integer entries, compared with ==). Both are construction checks: the closure is number
      conservation, built into every setting, and the dissipator's form is the Hamming distance's definition.
      The controls show that the comparison and the closure test can fail: on the three-height block (2, 2) at
      N = 4 the same form, written out by hand, leaves 4 gamma on the 30 cells below the top; site rates
      (1, 2, 1, 1) give the dissipator on (1, 1) three values where the uniform mean rate gives two, so no split
      with one refill exists there; an X field leaks out of the block.
  T3  the renewal as assembled, exactly, in the Laplace domain over Q(i), at two points (gamma, z), x0 a fixed
      seed vector: (z - L)^-1 x0 from the block's generator equals R0 x0 + 4 gamma R0 P Y, with R0 = (z + 2 gamma h_+ - L_H)^-1 and Y the
      solution of the floor's closed equation (1 - 4 gamma P R0 P) Y = P R0 x0. The identity is algebraic once
      the split holds; the controls carry its weight: a refill of 2 gamma instead of 4 gamma; the same
      construction on the three-height block (2, 2) at N = 4.
  T4  the kernel K(tau) = P e^{tau L_H} P. Its first-order term P L_H P is anti-Hermitian in every setting. On the
      diagonal blocks (1, 1) and (N - 1, N - 1) it is zero, and the Laplace transform P (z - L_H)^-1 P at a real
      z is real, entrywise non-negative, with every column summing to 1/z (what the transform of the stochastic
      |G|^2 must satisfy); the same test fails on every other bi-chromatic block, its negative leg. There the
      off-diagonal part of P L_H P, the hops that stay on the floor, is non-zero in every hopping setting and
      zero for an Ising-only H (ZZ and fields, the control that this part needs a hop); the whole term is purely
      imaginary with a transform that is not real for each real H, and real with a negative entry and a real
      transform for the bare DM chain, the control that "not real" fails for a purely imaginary H while "not stochastic" (a
      non-zero anti-Hermitian first-order term) holds throughout. Construction checks among these, forced by the
      setup rather than read: the anti-Hermitian form (H = H^dagger), the zero on the diagonal blocks (L_H has no
      population-to-population entry), the purely imaginary term for a real H, the DM leg's realness and negative
      entry (H's reality class), and the Ising control's missing hops (a diagonal H). The reads are the
      probability test with its negative leg, the non-zero hops in every hopping setting, the single bond
      included, and the transform that is not real in each real-H setting.
  T5  the frozen modes, a re-check of PROOF_R90_FROZEN_DIVISOR §10 and §7. At uniform gamma on (1, 1),
      L + 4 gamma = L_H + 4 gamma P, the generator of the renewal's S(t), and its kernel equals the intersection
      of the kernels of L_H and P as subspaces (the nullity of the three stacked equals both nullities), of
      dimension N // 2 on the Heisenberg and the XY chain at N = 3..6, exact ranks over Q(i). Control: on the
      R90 locus off the uniform point, rates (1, 2, 0, 1) at N = 4 and (1, 2, 1, 0, 1) at N = 5, the kernel of
      L + 4 gamma-bar and the intersection differ. The equality at the uniform point is forced by
      Re<v, (L_H + 4 gamma P) v> = 4 gamma |P v|^2; the N // 2 count and the control are the reads.

Usage: python simulations/f126_bichromatic_blocks.py
"""
import sys

import numpy as np
from sympy import QQ, QQ_I
from sympy.polys.matrices import DomainMatrix

FAILS = []


def check(name, ok):
    print(("PASS  " if ok else "FAIL  ") + name, flush=True)
    if not ok:
        FAILS.append(name)


def chromaticity(N, p, q):
    return min(p, q, N - p, N - q) + 1


def heights_closed(N, p, q):
    return {p + q - 2 * k for k in range(max(0, p + q - N), min(p, q) + 1)}


def hamiltonian(N, bonds=(), fields=(), dm=(), xfield=()):
    """H on all 2^N states: bonds (i, j, J, Delta) as J (XX + YY + Delta ZZ), fields h_l Z_l, DM bonds (i, j, D) as
    D (X_i Y_j - Y_i X_j), whose hop |0_i 1_j> -> |1_i 0_j> carries -2iD, and an optional transverse field
    hx_l X_l (the control that breaks the excitation number)."""
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for x in range(d):
        for i, j, J, Dl in bonds:
            bi, bj = (x >> i) & 1, (x >> j) & 1
            H[x, x] += J * Dl * (1 if bi == bj else -1)
            if bi != bj:
                H[x ^ (1 << i) ^ (1 << j), x] += 2 * J
        for i, j, D in dm:
            bi, bj = (x >> i) & 1, (x >> j) & 1
            if (bi, bj) == (0, 1):
                H[x ^ (1 << i) ^ (1 << j), x] += -2j * D
            elif (bi, bj) == (1, 0):
                H[x ^ (1 << i) ^ (1 << j), x] += 2j * D
        for l, h in enumerate(fields):
            H[x, x] += h * (1 - 2 * ((x >> l) & 1))
        for l, hx in enumerate(xfield):
            H[x ^ (1 << l), x] += hx
    assert np.array_equal(H, H.conj().T)
    return H


def full_liouvillian(N, H, gammas):
    """L on every cell |i><j| (index i * 2^N + j): -i[H, .] + sum_l gamma_l (Z_l . Z_l - .), the dissipator
    read off the signs of Z_l on ket and bra, not off the Hamming distance."""
    d = 2 ** N
    I = np.eye(d)
    LH = -1j * (np.kron(H, I) - np.kron(I, H.T))
    z = np.array([[1 - 2 * ((x >> l) & 1) for l in range(N)] for x in range(d)])
    D = np.zeros(d * d)
    for i in range(d):
        for j in range(d):
            D[i * d + j] = sum(g * (z[i, l] * z[j, l] - 1) for l, g in enumerate(gammas))
    return LH, LH + np.diag(D)


def block_cells(N, p, q):
    d = 2 ** N
    kets = [x for x in range(d) if bin(x).count("1") == p]
    bras = [x for x in range(d) if bin(x).count("1") == q]
    cells = [(i, j) for i in kets for j in bras]
    return cells, np.array([i * d + j for i, j in cells]), np.array([bin(i ^ j).count("1") for i, j in cells])


def to_dm(A):
    """A numpy matrix with Gaussian-integer entries as an exact matrix over Q(i); refuses anything else."""
    rows, cols = A.shape
    out = []
    for r in range(rows):
        row = []
        for c in range(cols):
            re, im = A[r, c].real, A[r, c].imag
            assert re == int(re) and im == int(im), (r, c, A[r, c])
            row.append(QQ_I(int(re), int(im)))
        out.append(row)
    return DomainMatrix(out, (rows, cols), QQ_I)


def diag_dm(vals):
    n = len(vals)
    return DomainMatrix([[vals[r] if r == c else QQ_I(0) for c in range(n)] for r in range(n)], (n, n), QQ_I)


def chain(N, J=1, Dl=1):
    return [(i, i + 1, J, Dl) for i in range(N - 1)]


def settings(N):
    """(name, bonds, fields, dm, kind): kind "real" for a real H, "imaginary" for a purely imaginary one (the bare DM
    chain), "complex" otherwise (the DM chain with a Z field)."""
    rng = np.random.default_rng(126 + N)
    E = [(i, j) for i in range(N) for j in range(i + 1, N) if rng.random() < 0.6]
    E = sorted(set(E) | {(i, i + 1) for i in range(N - 1)})
    dm = [(i, i + 1, 1) for i in range(N - 1)]
    return [("Heisenberg chain", chain(N), (), (), "real"),
            ("XXZ chain, Delta = 2, Z fields", chain(N, 1, 2), (1, 0, -2, 3, 1, -1)[:N], (), "real"),
            ("ring topology, couplings 1, 2, 1, 3, ...", [(i, (i + 1) % N, (1, 2, 1, 3, 1, 2)[i], 1) for i in range(N)],
             (), (), "real"),
            ("XY star", [(0, k, 1, 0) for k in range(1, N)], (), (), "real"),
            ("random graph, integer J, Delta, fields",
             [(i, j, int(rng.integers(1, 4)), int(rng.integers(0, 3))) for i, j in E],
             tuple(int(v) for v in rng.integers(-2, 3, N)), (), "real"),
            ("single hopping bond (0, 1)", [(0, 1, 1, 0)], (), (), "real"),
            ("DM chain", (), (), dm, "imaginary"),
            ("DM chain with a Z field", (), (2, 0, -1, 1, 0, 1)[:N], dm, "complex")]


def bichromatic(N):
    return [(p, q) for p in range(N + 1) for q in range(N + 1) if chromaticity(N, p, q) == 2]


def floor_index(h):
    lo = int(min(h))
    return [k for k in range(len(h)) if h[k] == lo]


if __name__ == "__main__":
    # ---------------------------------------------------------------- T1
    print("T1  the chromaticity of the joint-popcount blocks (a re-check of the owned count)", flush=True)
    ok_sets, ok_counts, refuted = True, True, []
    for N in range(1, 11):
        d = 2 ** N
        table = np.array([bin(x).count("1") for x in range(d)])
        ham = table[np.arange(d)[:, None] ^ np.arange(d)[None, :]]
        sizes = {}
        wrong_hit = False
        for p in range(N + 1):
            for q in range(N + 1):
                mask = (table[:, None] == p) & (table[None, :] == q)
                seen = set(int(h) for h in np.unique(ham[mask]))
                if seen != heights_closed(N, p, q) or len(seen) != chromaticity(N, p, q):
                    ok_sets = False
                sizes[len(seen)] = sizes.get(len(seen), 0) + 1
                if len(seen) != min(p, q) + 1:
                    wrong_hit = True
        two_expected = 4 * N - 8 if N >= 3 else (1 if N == 2 else 0)
        if sizes.get(1, 0) != 4 * N or sizes.get(2, 0) != two_expected:
            ok_counts = False
        refuted.append(wrong_hit)
    check("N = 1..10, every block (p, q) by brute force: its heights are {p + q - 2k}, min(p, q, N - p, N - q) + 1 "
          "of them", ok_sets)
    check("N = 1..10: chromaticity 1 on 4N blocks, 2 on 4N - 8 for N >= 3 (the block (1, 1) alone at N = 2)",
          ok_counts)
    check("control: the count min(p, q) + 1 is refuted by the brute force at every N = 1..10 "
          f"(refuted at N = {[N for N, r in zip(range(1, 11), refuted) if r]})", refuted == [True] * 10)

    # ---------------------------------------------------------------- T2
    print("\nT2  the split: dissipator = -2 gamma h_+ + 4 gamma P on every bi-chromatic block", flush=True)
    for N in (3, 4, 5):
        blocks = bichromatic(N)
        for name, bonds, fields, dm, _ in settings(N):
            H = hamiltonian(N, bonds, fields, dm)
            bad = []
            for gamma in (1, 3):
                LH, L = full_liouvillian(N, H, [gamma] * N)
                for p, q in blocks:
                    cells, idx, h = block_cells(N, p, q)
                    others = np.setdiff1d(np.arange(L.shape[0]), idx)
                    leak = np.abs(L[np.ix_(others, idx)]).max()
                    lo, hi = sorted(set(int(v) for v in h))
                    diss = L[np.ix_(idx, idx)] - LH[np.ix_(idx, idx)]
                    split = np.diag(-2.0 * gamma * hi + 4.0 * gamma * (h == lo))
                    if not (leak == 0.0 and hi == lo + 2 and np.array_equal(diss, split)):
                        bad.append((gamma, p, q))
            check(f"N = {N}, {name}: all {len(blocks)} bi-chromatic blocks closed under L, dissipator equal to the "
                  "split cell by cell, gamma = 1 and 3" + (f"; failing {bad}" if bad else ""), not bad)
    N = 4
    H = hamiltonian(N, chain(N))
    LH, L = full_liouvillian(N, H, [1] * N)
    cells, idx, h = block_cells(N, 2, 2)
    res = (L[np.ix_(idx, idx)] - LH[np.ix_(idx, idx)]) - np.diag(-2.0 * 4 + 4.0 * (h == 0))
    off = np.abs(res - np.diag(np.diag(res))).max()
    check("control: the block (2, 2) at N = 4 (heights 0, 2, 4), the same form with h_+ = 4 and P on height 0: the "
          f"residual is 4 gamma exactly on the {int((h < 4).sum())} cells below the top height ({int((h == 2).sum())} at "
          f"height 2, {int((h == 0).sum())} at height 0) and 0 on the {int((h == 4).sum())} at the top",
          off == 0.0 and np.array_equal(np.diag(res), 4.0 * (h < 4)) and int((h == 2).sum()) == 24
          and int((h == 0).sum()) == 6)
    cells, idx, h = block_cells(N, 1, 1)
    LHp, Lp = full_liouvillian(N, H, [1, 2, 1, 1])
    vp = sorted(set(np.diag(Lp[np.ix_(idx, idx)] - LHp[np.ix_(idx, idx)]).real.tolist()))
    LHm, Lm = full_liouvillian(N, H, [1.25] * 4)
    vm = sorted(set(np.diag(Lm[np.ix_(idx, idx)] - LHm[np.ix_(idx, idx)]).real.tolist()))
    check(f"control: site rates (1, 2, 1, 1) give the dissipator on (1, 1) at N = 4 the values {vp}, three, so no "
          f"scalar plus one refill exists; the uniform mean rate 5/4 gives {vm}, two",
          vp == [-6.0, -4.0, 0.0] and vm == [-5.0, 0.0])
    Hx = hamiltonian(N, chain(N), xfield=(1, 0, 0, 0))
    _, Lx = full_liouvillian(N, Hx, [1] * N)
    cells, idx, h = block_cells(N, 1, 2)
    others = np.setdiff1d(np.arange(Lx.shape[0]), idx)
    leak = np.abs(Lx[np.ix_(others, idx)]).max()
    check(f"control: an X field on site 0 leaks out of the block (1, 2) at N = 4 (largest entry {leak}), not zero",
          leak > 0.0)

    # ---------------------------------------------------------------- T3
    print("\nT3  the renewal, exact over Q(i) in the Laplace domain", flush=True)
    points = [(QQ_I(QQ(3, 2), 0), QQ_I(QQ(7, 3), QQ(1, 5))), (QQ_I(QQ(1, 3), 0), QQ_I(QQ(1, 2), -2))]

    def pieces(N, H, p, q, gamma):
        LH, _ = full_liouvillian(N, H, [0] * N)
        cells, idx, h = block_cells(N, p, q)
        n = len(idx)
        d = 2 ** N
        zs = [[1 - 2 * ((x >> l) & 1) for l in range(N)] for x in range(d)]
        Dexact = [gamma * QQ_I(sum(zs[i][l] * zs[j][l] - 1 for l in range(N)), 0) for i, j in cells]
        LHb = to_dm(LH[np.ix_(idx, idx)])
        x0 = DomainMatrix([[QQ_I(k % 3 - 1, (k * 7) % 5 - 2)] for k in range(n)], (n, 1), QQ_I)
        return LHb, LHb + diag_dm(Dexact), h, x0, n

    def renewal(LHb, h, x0, n, hi, lo, refill, gamma, z):
        A0 = diag_dm([z + QQ_I(2, 0) * gamma * QQ_I(int(hi), 0)] * n) - LHb
        floor = [k for k in range(n) if h[k] == lo]
        f = len(floor)
        E = DomainMatrix([[QQ_I(1) if r == floor[c] else QQ_I(0) for c in range(f)] for r in range(n)], (n, f), QQ_I)
        C = A0.lu_solve(E)
        r0x = A0.lu_solve(x0)
        Cff = DomainMatrix([[C[floor[r], c].element for c in range(f)] for r in range(f)], (f, f), QQ_I)
        rhs = DomainMatrix([[r0x[floor[r], 0].element] for r in range(f)], (f, 1), QQ_I)
        Y = (diag_dm([QQ_I(1)] * f) - Cff * refill).lu_solve(rhs)
        return r0x + C * Y * refill

    for N in (3, 4, 5):
        blocks = bichromatic(N)
        for name, bonds, fields, dm, _ in settings(N):
            H = hamiltonian(N, bonds, fields, dm)
            bad = []
            for gamma, z in points:
                for p, q in blocks:
                    LHb, Lb, h, x0, n = pieces(N, H, p, q, gamma)
                    lo, hi = sorted(set(int(v) for v in h))
                    left = (diag_dm([z] * n) - Lb).lu_solve(x0)
                    if left != renewal(LHb, h, x0, n, hi, lo, QQ_I(4, 0) * gamma, gamma, z):
                        bad.append((p, q))
            check(f"N = {N}, {name}: the renewal equals the resolvent exactly on all {len(blocks)} bi-chromatic "
                  "blocks, at (gamma, z) = (3/2, 7/3 + i/5) and (1/3, 1/2 - 2i)" + (f"; failing {bad}" if bad else ""),
                  not bad)
    N = 4
    H = hamiltonian(N, chain(N))
    gamma, z = points[0]
    LHb, Lb, h, x0, n = pieces(N, H, 1, 2, gamma)
    left = (diag_dm([z] * n) - Lb).lu_solve(x0)
    check("control: the block (1, 2) at N = 4 with the refill 2 gamma instead of 4 gamma differs from the resolvent",
          left != renewal(LHb, h, x0, n, 3, 1, QQ_I(2, 0) * gamma, gamma, z))
    LHb, Lb, h, x0, n = pieces(N, H, 2, 2, gamma)
    left = (diag_dm([z] * n) - Lb).lu_solve(x0)
    check("control: the three-height block (2, 2) at N = 4, envelope at height 4 and refill 4 gamma onto height 0, "
          "differs from the resolvent", left != renewal(LHb, h, x0, n, 4, 0, QQ_I(4, 0) * gamma, gamma, z))

    # ---------------------------------------------------------------- T4
    print("\nT4  the kernel on the floor: a probability ladder only where the floor is populations", flush=True)
    zr = QQ_I(QQ(5, 2), 0)

    def floor_reads(N, H, p, q):
        """The first-order term F = P L_H P on the floor, and the transform P (z - L_H)^-1 P at z = 5/2 read as
        (real, the probability test: real, entrywise non-negative, columns summing to 1/z)."""
        LH, _ = full_liouvillian(N, H, [0] * N)
        cells, idx, h = block_cells(N, p, q)
        floor = floor_index(h)
        LHblk = LH[np.ix_(idx, idx)]
        F = LHblk[np.ix_(floor, floor)]
        n, f = len(idx), len(floor)
        E = DomainMatrix([[QQ_I(1) if r == floor[c] else QQ_I(0) for c in range(f)] for r in range(n)], (n, f), QQ_I)
        C = (diag_dm([zr] * n) - to_dm(LHblk)).lu_solve(E)
        K = [[C[floor[r], c].element for c in range(f)] for r in range(f)]
        real = all(K[r][c].y == 0 for r in range(f) for c in range(f))
        nonneg = all(K[r][c].x >= 0 for r in range(f) for c in range(f))
        sums = all(sum((K[r][c] for r in range(f)), QQ_I(0)) == QQ_I(1, 0) / zr for c in range(f))
        return F, real, real and nonneg and sums

    for N in (3, 4, 5):
        blocks = bichromatic(N)
        for name, bonds, fields, dm, kind in settings(N):
            H = hamiltonian(N, bonds, fields, dm)
            bad, reads = [], []
            for p, q in blocks:
                F, real, prob = floor_reads(N, H, p, q)
                antiherm = np.array_equal(F, -F.conj().T)
                hops = F - np.diag(np.diag(F))
                if p == q:
                    ok = antiherm and not F.any() and prob
                else:
                    ok = antiherm and hops.any() and not prob
                    if kind == "real":
                        ok = ok and not F.real.any() and not real
                    elif kind == "imaginary":
                        # L_H is real: the first-order term is real and antisymmetric, so it has a negative entry,
                        # and the kernel is real; "not real" is not the general statement, "not stochastic" is
                        ok = ok and not F.imag.any() and bool((F.real < 0).any()) and real
                    else:
                        reads.append(real)
                if not ok:
                    bad.append((p, q))
            diag_blocks = [b for b in blocks if b[0] == b[1]]
            head = (f"N = {N}, {name}: P L_H P anti-Hermitian on all {len(blocks)}; on {diag_blocks} it is zero and the "
                    "transform passes the probability test (real, non-negative, columns summing to 1/z); on the other "
                    f"{len(blocks) - len(diag_blocks)} the test fails and the hops on the floor are non-zero")
            tail = {"real": ", the term purely imaginary and the transform not real",
                    "imaginary": ", the term real with a negative entry and the transform real (the control: not "
                                 "complex, still not stochastic)",
                    "complex": f" (the transform real on {sum(reads)} of {len(reads)}, read)"}[kind]
            check(head + tail + (f"; failing {bad}" if bad else ""), not bad)
    # control: an Ising-only H (ZZ bonds and Z fields, no hop) has no hop on the floor
    for N in (3, 4, 5):
        H = np.diag([sum((1 if ((x >> i) & 1) == ((x >> (i + 1)) & 1) else -1) for i in range(N - 1))
                     + sum(hz * (1 - 2 * ((x >> l) & 1)) for l, hz in enumerate((1, 0, -2, 3, 1)[:N]))
                     for x in range(2 ** N)]).astype(complex)
        off = [(p, q) for p, q in bichromatic(N) if p != q]
        firsts = [floor_reads(N, H, p, q)[0] for p, q in off]
        hop_free = all(not (F - np.diag(np.diag(F))).any() for F in firsts)
        diag_on = sum(1 for F in firsts if np.diag(F).any())
        check(f"control: N = {N}, an Ising chain with Z fields (no hop): no hop on the floor of any of the {len(off)} "
              f"off-diagonal bi-chromatic blocks (the diagonal, the energy differences, non-zero on {diag_on})",
              hop_free and diag_on > 0)

    # ---------------------------------------------------------------- T5
    print("\nT5  the frozen modes at the uniform point: ker(L + 4 gamma) = ker L_H and ker P on (1, 1)", flush=True)

    def frozen_pieces(N, H, gammas, shift):
        """On (1, 1): M = L + shift, L_H, P as exact matrices; the nullities of M, of [L_H; P] and of [M; L_H; P]."""
        LH, L = full_liouvillian(N, H, gammas)
        cells, idx, h = block_cells(N, 1, 1)
        n = len(idx)
        LHb = to_dm(LH[np.ix_(idx, idx)])
        Pd = diag_dm([QQ_I(1) if h[k] == 0 else QQ_I(0) for k in range(n)])
        M = to_dm(L[np.ix_(idx, idx)]) + diag_dm([shift] * n)
        return n - M.rank(), n - LHb.vstack(Pd).rank(), n - M.vstack(LHb).vstack(Pd).rank()

    for name, mk in (("Heisenberg chain", lambda N: chain(N)), ("XY chain", lambda N: chain(N, 1, 0))):
        rows, ok = [], True
        for N in (3, 4, 5, 6):
            H = hamiltonian(N, mk(N))
            k_m, k_i, k_t = frozen_pieces(N, H, [1] * N, QQ_I(4, 0))
            rows.append((N, k_m, k_i, k_t))
            ok = ok and k_m == k_i == k_t == N // 2
        check(f"{name}, gamma = 1: on (1, 1) the nullities of L + 4 gamma, of [L_H; P] and of all three stacked are "
              "equal and N // 2, so the kernel is the intersection: "
              + ", ".join(f"N = {N}: {a}, {b}, {c}" for N, a, b, c in rows), ok)
    rows, ok = [], True
    for N, prof in ((4, (1, 2, 0, 1)), (5, (1, 2, 1, 0, 1))):
        H = hamiltonian(N, chain(N))
        k_m, k_i, k_t = frozen_pieces(N, H, list(prof), QQ_I(4, 0))
        rows.append((N, prof, k_m, k_i, k_t))
        ok = ok and k_m >= N // 2 and k_t < k_m
    check("control: on the R90 locus off the uniform point (mean rate 1, the Heisenberg chain), the kernel of "
          "L + 4 gamma-bar on (1, 1) is not the intersection: "
          + ", ".join(f"N = {N}, rates {pr}: nullities {a}, {b}, stacked {c}" for N, pr, a, b, c in rows), ok)

    print("\nALL GATES PASS" if not FAILS else f"\n{len(FAILS)} FAILED")
    sys.exit(1 if FAILS else 0)
