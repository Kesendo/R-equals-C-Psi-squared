"""XOR sector verification: the invariant, exact checks.

The XOR sector is the span of the N+1 operators v_k = X^N P_k, k = 0..N, where
P_k is the projector onto popcount-k basis states. Every check below is exact:
H and the dephasing rates are integers, so every entry of L is a Gaussian
integer and every residual in float64 is an exact integer computation.

1. Eigenvectors. Each v_k is a right AND a left eigenvector of L at -2*Sigma
   (Sigma = sum of the rates) whenever [H, X^N] = 0 and [H, sum_l Z_l] = 0:
   Heisenberg, XY and XXZ couplings on chain, ring, star and complete graphs,
   non-uniform rates, N = 2..5. Residuals ||L v + 2 Sigma v|| and
   ||L^dag v + 2 Sigma v|| must be exactly 0.0.
2. Multiplicity. The eigenvalues -2*Sigma and 0 each have algebraic multiplicity
   exactly N+1 (semisimple), by ranks over GF(p) with i -> sqrt(-1) mod p:
   rank mod p can only drop, so nullity(M) <= N+1 and nullity(M^2) <= N+1 mod p,
   together with the N+1 explicit eigenvectors, give nullity exactly N+1 over Q
   and no Jordan chains. Heisenberg and XY chains and XXZ on the complete graph
   with random integer couplings, N = 2..4.
3. Shares. Since v_k are left and right eigenvectors, their span reduces L and
   the orthogonal projection onto it is an invariant share. GHZ_N minus its
   stationary populations lies 100% in the sector at every N; W_N lies 0% in it
   for N >= 3 (100% at N = 2, where Hamming distance 2 is N). Exact in Fractions.
4. Result 4: every v_k is purely off-diagonal (support |a><a-bar|).
5. Controls. A transverse field, or XYZ with Jx != Jy, leaves no v_k an
   eigenvector; a non-uniform longitudinal field or a DM term keeps exactly
   k = 0 and k = N (the GHZ coherences); a uniform field h*sum Z keeps all N+1,
   shifted to -2*Sigma + 2ih(N - 2k).

Run: python simulations/xor_verify.py   (prints PASS/FAIL, exits 1 on any FAIL)
"""

import sys
from fractions import Fraction
from itertools import combinations

import numpy as np

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)

P_MOD = 1000037            # prime, 1 mod 4
I_MOD = None               # a square root of -1 mod P_MOD, set below


def site_op(op, s, n):
    m = np.array([[1]], dtype=complex)
    for q in range(n):
        m = np.kron(m, op if q == s else I2)
    return m


def bonds(n, topo):
    if topo == "chain":
        return [(i, i + 1) for i in range(n - 1)]
    if topo == "ring":
        return [(i, (i + 1) % n) for i in range(n)]
    if topo == "star":
        return [(0, i) for i in range(1, n)]
    if topo == "complete":
        return list(combinations(range(n), 2))
    raise ValueError(topo)


def build_h(n, topo, jxy, jz, jx=None):
    """sum over bonds of jx*XX + jxy*YY + jz*ZZ (jx defaults to jxy). Integer couplings."""
    d = 2 ** n
    h = np.zeros((d, d), dtype=complex)
    for b, (i, j) in enumerate(bonds(n, topo)):
        cxy = jxy[b] if isinstance(jxy, (list, tuple)) else jxy
        cz = jz[b] if isinstance(jz, (list, tuple)) else jz
        cx = cxy if jx is None else jx
        h += cx * site_op(SX, i, n) @ site_op(SX, j, n)
        h += cxy * site_op(SY, i, n) @ site_op(SY, j, n)
        h += cz * site_op(SZ, i, n) @ site_op(SZ, j, n)
    return h


def build_l(h, gammas, n):
    """Row-major vec(rho): L = -i(H x I - I x H^T) + sum_l g_l (Z_l x Z_l - I)."""
    d = 2 ** n
    lv = -1j * (np.kron(h, np.eye(d)) - np.kron(np.eye(d), h.T))
    for l in range(n):
        z = site_op(SZ, l, n)
        lv += gammas[l] * (np.kron(z, z) - np.eye(d * d))
    return lv


def popcount_projector(n, k):
    d = 2 ** n
    p = np.zeros((d, d), dtype=complex)
    for a in range(d):
        if bin(a).count("1") == k:
            p[a, a] = 1
    return p


def xn(n):
    m = np.array([[1]], dtype=complex)
    for _ in range(n):
        m = np.kron(m, SX)
    return m


def xor_vectors(n):
    x = xn(n)
    return [(x @ popcount_projector(n, k)).reshape(-1) for k in range(n + 1)]


PASS = 0
FAIL = 0


def check(name, ok, detail=""):
    global PASS, FAIL
    if ok:
        PASS += 1
        print(f"  PASS: {name}")
    else:
        FAIL += 1
        print(f"  FAIL: {name} {detail}")


def residuals(lv, v, target):
    right = np.linalg.norm(lv @ v - target * v)
    left = np.linalg.norm(lv.conj().T @ v - np.conj(target) * v)
    return right, left


def rank_mod_p(m_re, m_im):
    """Rank over GF(p) of the Gaussian-integer matrix m_re + i*m_im, i -> I_MOD."""
    a = (np.mod(m_re, P_MOD) + I_MOD * np.mod(m_im, P_MOD)) % P_MOD
    a = a.astype(np.int64)
    rows, cols = a.shape
    rank = 0
    for c in range(cols):
        piv = None
        for r in range(rank, rows):
            if a[r, c] != 0:
                piv = r
                break
        if piv is None:
            continue
        a[[rank, piv]] = a[[piv, rank]]
        inv = pow(int(a[rank, c]), P_MOD - 2, P_MOD)
        a[rank] = (a[rank] * inv) % P_MOD
        col = a[:, c].copy()
        col[rank] = 0
        nz = np.nonzero(col)[0]
        if nz.size:
            a[nz] = (a[nz] - (col[nz, None] * a[rank][None, :]) % P_MOD) % P_MOD
        rank += 1
        if rank == rows:
            break
    return rank


def gaussian_parts(m):
    re = np.rint(m.real).astype(np.int64)
    im = np.rint(m.imag).astype(np.int64)
    assert np.array_equal(re, m.real) and np.array_equal(im, m.imag), "matrix is not Gaussian-integer"
    return re, im


def exact_share(rho, n):
    """Fractions: share of (rho minus its projection onto span{P_k}) inside span{v_k}."""
    d = 2 ** n
    rho = [[Fraction(x) for x in row] for row in rho]
    # stationary part: projection onto span{P_k} (orthogonal, P_k diagonal)
    ns = [row[:] for row in rho]
    for k in range(n + 1):
        idx = [a for a in range(d) if bin(a).count("1") == k]
        mean = sum(rho[a][a] for a in idx) / len(idx)
        for a in idx:
            ns[a][a] -= mean
    total = sum(x * x for row in ns for x in row)
    # projection onto span{v_k}: v_k has ones at (a-bar, a) for popcount(a) = k
    inside = Fraction(0)
    full = d - 1
    for k in range(n + 1):
        idx = [a for a in range(d) if bin(a).count("1") == k]
        overlap = sum(ns[full ^ a][a] for a in idx)
        inside += overlap * overlap / len(idx)
    return inside, total


def main():
    global I_MOD
    I_MOD = pow(2, (P_MOD - 1) // 4, P_MOD)
    while (I_MOD * I_MOD) % P_MOD != P_MOD - 1:
        I_MOD = pow(I_MOD + 1, (P_MOD - 1) // 4, P_MOD)

    rng = np.random.default_rng(20260926)

    print("=" * 70)
    print("CHECK 1: X^N P_k are right and left eigenvectors at -2*Sigma (exact)")
    print("=" * 70)
    for n in [2, 3, 4, 5]:
        gammas = list(range(1, n + 1))                     # non-uniform, integer
        sigma = sum(gammas)
        nb = {t: len(bonds(n, t)) for t in ["chain", "ring", "star", "complete"]}
        systems = [
            ("Heisenberg chain", build_h(n, "chain", 1, 1)),
            ("XY chain", build_h(n, "chain", 1, 0)),
            ("XXZ(Delta=3) chain", build_h(n, "chain", 2, 6)),
            ("Heisenberg ring", build_h(n, "ring", 1, 1)),
            ("Heisenberg star", build_h(n, "star", 1, 1)),
            ("XXZ complete, random J and Delta",
             build_h(n, "complete", [int(x) for x in rng.integers(1, 5, nb["complete"])],
                     [int(x) for x in rng.integers(-4, 5, nb["complete"])])),
        ]
        for name, h in systems:
            lv = build_l(h, gammas, n)
            worst = max(max(residuals(lv, v, -2 * sigma)) for v in xor_vectors(n))
            check(f"N={n} {name}, gamma={gammas}: max residual = {worst}", worst == 0.0)

    print("\n" + "=" * 70)
    print("CHECK 2: multiplicity of -2*Sigma and of 0 is exactly N+1, semisimple (GF(p) ranks)")
    print("=" * 70)
    for n in [2, 3, 4]:
        gammas = list(range(1, n + 1))
        sigma = sum(gammas)
        nbc = len(bonds(n, "complete"))
        for name, h in [("Heisenberg chain", build_h(n, "chain", 1, 1)),
                        ("XY chain", build_h(n, "chain", 1, 0)),
                        ("XXZ complete, random J and Delta",
                         build_h(n, "complete", [int(x) for x in rng.integers(1, 5, nbc)],
                                 [int(x) for x in rng.integers(-4, 5, nbc)]))]:
            lv = build_l(h, gammas, n)
            dim = lv.shape[0]
            # the lower bound N+1 on this very instance: v_k at -2*Sigma, P_k at 0
            pk = [popcount_projector(n, k).reshape(-1) for k in range(n + 1)]
            low = max(max(np.linalg.norm(lv @ v + 2 * sigma * v) for v in xor_vectors(n)),
                      max(np.linalg.norm(lv @ p) for p in pk))
            check(f"N={n} {name}: the N+1 explicit eigenvectors at -2*Sigma and at 0, max residual = {low}",
                  low == 0.0)
            for label, shift in [("-2*Sigma", 2 * sigma), ("0", 0)]:
                m = lv + shift * np.eye(dim)
                re, im = gaussian_parts(m)
                m2 = m @ m
                re2, im2 = gaussian_parts(m2)
                null1 = dim - rank_mod_p(re, im)
                null2 = dim - rank_mod_p(re2, im2)
                check(f"N={n} {name}: nullity at {label} mod p = {null1}, of the square = {null2} "
                      f"(expect N+1 = {n + 1} for both)", null1 == n + 1 and null2 == n + 1)

    print("\n" + "=" * 70)
    print("CHECK 3: invariant shares in the XOR sector (exact Fractions)")
    print("=" * 70)
    for n in [2, 3, 4, 5]:
        d = 2 ** n
        w = [[Fraction(0)] * d for _ in range(d)]
        for i in range(n):
            for j in range(n):
                w[1 << i][1 << j] = Fraction(1, n)
        gi, gt = exact_share([[Fraction(1, 2) if (a in (0, d - 1) and b in (0, d - 1)) else Fraction(0)
                               for b in range(d)] for a in range(d)], n)
        wi, wt = exact_share(w, n)
        check(f"N={n}: GHZ non-stationary share = {gi / gt} (expect 1)", gt > 0 and gi == gt)
        expect_w = Fraction(1) if n == 2 else Fraction(0)
        check(f"N={n}: W non-stationary share = {wi / wt} (expect {expect_w})", wt > 0 and wi / wt == expect_w)

    print("\n" + "=" * 70)
    print("CHECK 4: every X^N P_k is purely off-diagonal (Result 4)")
    print("=" * 70)
    for n in [2, 3, 4, 5]:
        d = 2 ** n
        diag_weight = max(np.abs(np.diag(v.reshape(d, d))).max() for v in xor_vectors(n))
        check(f"N={n}: largest diagonal entry of any v_k = {diag_weight}", diag_weight == 0.0)

    print("\n" + "=" * 70)
    print("CHECK 5: controls (which H keep which v_k)")
    print("=" * 70)
    for n in [3, 4]:
        gammas = list(range(1, n + 1))
        sigma = sum(gammas)
        base = build_h(n, "chain", 1, 1)
        field_nonuni = sum((l + 1) * site_op(SZ, l, n) for l in range(n))
        dm = np.zeros_like(base)
        for (i, j) in bonds(n, "chain"):
            dm += site_op(SX, i, n) @ site_op(SY, j, n) - site_op(SY, i, n) @ site_op(SX, j, n)
        controls = [
            ("transverse field X_0", base + site_op(SX, 0, n), set()),
            ("XYZ Jx=2, Jy=1, Jz=1", build_h(n, "chain", 1, 1, jx=2), set()),
            ("random Hermitian (integer)", None, set()),
            ("non-uniform longitudinal field", base + field_nonuni, {0, n}),
            ("DM term", base + dm, {0, n}),
        ]
        for name, h, keep in controls:
            if h is None:
                a = rng.integers(-3, 4, (2 ** n, 2 ** n)) + 1j * rng.integers(-3, 4, (2 ** n, 2 ** n))
                h = a + a.conj().T
            lv = build_l(h, gammas, n)
            kept = set()
            for k, v in enumerate(xor_vectors(n)):
                lvv = lv @ v
                # eigenvector at any eigenvalue: L v parallel to v on v's support, zero elsewhere
                support = v != 0
                if np.all(lvv[~support] == 0):
                    ratios = lvv[support] / v[support]
                    if np.all(ratios == ratios[0]):
                        kept.add(k)
            check(f"N={n} {name}: v_k eigenvectors for k in {sorted(kept)} (expect {sorted(keep)})",
                  kept == keep)
        h_uni = base + 3 * sum(site_op(SZ, l, n) for l in range(n))
        lv = build_l(h_uni, gammas, n)
        worst = 0.0
        for k, v in enumerate(xor_vectors(n)):
            target = -2 * sigma + 2j * 3 * (n - 2 * k)
            worst = max(worst, *residuals(lv, v, target))
        check(f"N={n} uniform field h=3: every v_k at -2*Sigma + 2ih(N-2k), max residual = {worst}",
              worst == 0.0)

    print(f"\n{PASS} passed, {FAIL} failed")
    return 0 if FAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
