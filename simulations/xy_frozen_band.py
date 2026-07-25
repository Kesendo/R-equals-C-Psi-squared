# -*- coding: utf-8 -*-
"""xy_frozen_band.py

Verifier for experiments/XY_FROZEN_BAND.md: where the F140 frozen root -4*gbar lives once the
chain is XY rather than Heisenberg, and what decides it.

Prints "XY frozen band: ALL GREEN" when every check passes. Runtime ~4-8 min.

  V1  the band law, by full (p,q) census at N = 4..7: the blocks carrying -4*gbar are exactly
      |p - q| in {0, 2} minus the two 1x1 corners (0,0) and (N,N), each carrying floor(N/2)
      and never a part of it; count 3*(N-1); the fold root's carriers are the image under the
      one-sided fold p -> N-p
  V2  the band EDGE, at the N that can see it: |p-q| = 4 carries nothing at N = 8, 9, 10, on
      two locus profiles and two couplings. Needed because at N <= 7 the rival "width
      floor(N/2)" is the SAME predicate (and so is "width 3", since p+q even forces |p-q| even)
  V3  why the edge sits there: the 4*gbar recentering makes the rate part tauQ-odd EXACTLY on
      the cells with |A^B| = 2, and |A^B| >= |p-q|, so a block with |p-q| >= 4 holds no such
      cell at all. Chain-independent, pure arithmetic plus combinatorics
  V4  the GATE: the off-diagonal band exists iff the single-excitation matrix h is bipartite
      (its spectrum symmetric about 0). Four independent ways: two R-invariant diagonals added
      to the chain, three N each, and the ring, where bipartiteness is the parity of N. The
      CORNER survives all five rows, which is the control that separates the two mechanisms
  V6  the OTHER falsified candidate, as a count: tauQ fixes a cell only at B = R(A), which
      forces |A| = |B|, so the off-diagonal band pairs have no fixed cell at all and the corner's
      room shortage predicts zero where floor(N/2) is measured; on the diagonal it misses in
      both directions
  V5  a FALSIFIED candidate, kept because it is the obvious one: T(rho) = d_a rho d_b on a
      chiral pair commutes with the Hamiltonian part exactly (checked to machine zero, and
      against the predicted commutator off the chiral locus) yet does NOT carry the corner's
      frozen vectors into the (0,2) frozen subspace. The dephasing is what it does not respect

Exactness: J and every gamma_l sit on the dyadic grid 1/1024, so 4*N*1024*(block entries) are
integers and the multiplicity is an exact GF(p) rank at two primes = 1 (mod 4) wherever the
block is small enough; above that the read is an SVD nullity, and V1 reports how many of its
reads were exact. Never an eigenvalue count: the departing modes crowd the root at spacing
J^(2d), which is exactly where a float spectrum lies.
"""
from itertools import combinations
from math import comb

import numpy as np

GRID = 1024
PRIMES = (998244353, 1004535809)
EXACT_MAX = 260
FLOAT_TOL = 1e-7

FAILURES = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILURES.append(name)


# ---------- profiles and single-excitation matrices ----------

def locus(n, gbar_num=92, half=(29, -13, 7, 19)):
    """An R90 locus profile as integer numerators over GRID: every reflection pair sums to
    2*gbar exactly, and the |delta| are distinct so no unintended pair balances."""
    g = [gbar_num] * n
    for i in range(n // 2):
        d = half[i % len(half)]
        g[i] = gbar_num + d
        g[n - 1 - i] = gbar_num - d
    return g


def h_chain(n, diag=None):
    h = np.zeros((n, n), dtype=np.int64)
    for a in range(n - 1):
        h[a, a + 1] = h[a + 1, a] = 1
    if diag is not None:
        for a in range(n):
            h[a, a] = diag[a]
    return h


def h_ring(n):
    h = h_chain(n)
    h[0, n - 1] = h[n - 1, 0] = 1
    return h


def is_bipartite_spectrum(h):
    lam = np.sort(np.linalg.eigvalsh(h.astype(float)))
    return bool(np.allclose(lam, -lam[::-1], atol=1e-9))


# ---------- blocks ----------

def sector(n, k, h):
    """H restricted to the k-excitation sector: hopping from h's off-diagonal plus h's
    on-site energies. The on-site term is what an R-invariant diagonal probe rides on;
    dropping it would make V4's diagonal cases silent no-ops."""
    st = list(combinations(range(n), k))
    idx = {s: i for i, s in enumerate(st)}
    H = np.zeros((len(st), len(st)), dtype=np.int64)
    for s, cfg in enumerate(st):
        occ = set(cfg)
        for site in cfg:
            H[s, s] += h[site, site]
            for tgt in range(n):
                if h[site, tgt] and tgt not in occ:
                    H[idx[tuple(sorted((occ - {site}) | {tgt}))], s] += h[site, tgt]
    return st, idx, H


def scaled_block(n, h, gnum, jnum, p, q):
    """4*N*GRID * (L_(p,q)) as Gaussian integers, in the cell basis (ket, bra)."""
    stk, ik, Hk = sector(n, p, h)
    stb, ib, Hb = sector(n, q, h)
    f = 4 * n
    d = len(stk) * len(stb)
    re = np.zeros((d, d), dtype=np.int64)
    im = np.zeros((d, d), dtype=np.int64)
    for a, ka in enumerate(stk):
        for b, kb in enumerate(stb):
            col = a * len(stb) + b
            for c in range(len(stk)):
                if Hk[c, a]:
                    im[c * len(stb) + b, col] += -f * jnum * Hk[c, a]
            for c in range(len(stb)):
                if Hb[b, c]:
                    im[a * len(stb) + c, col] += f * jnum * Hb[b, c]
            re[col, col] += -2 * f * sum(gnum[l] for l in set(ka) ^ set(kb))
    return re, im


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
    r = sqrt_minus_one(p)
    a = (((re % p) + r * (im % p)) % p).astype(object)
    d = a.shape[0]
    rank = 0
    for col in range(d):
        piv = next((row for row in range(rank, d) if a[row, col] % p), None)
        if piv is None:
            continue
        if piv != rank:
            a[[rank, piv]] = a[[piv, rank]]
        inv = modpow(int(a[rank, col]), p - 2, p)
        a[rank, col:] = (a[rank, col:] * inv) % p
        for row in range(rank + 1, d):
            f = int(a[row, col]) % p
            if f:
                a[row, col:] = (a[row, col:] - f * a[rank, col:]) % p
        rank += 1
        if rank == d:
            break
    return rank


def frozen(n, h, gnum, jnum, p, q, root_num):
    """dim ker(L_(p,q) - root), the root given as its 4*N*GRID-scaled integer."""
    re, im = scaled_block(n, h, gnum, jnum, p, q)
    d = re.shape[0]
    re = re - root_num * np.eye(d, dtype=np.int64)
    if d <= EXACT_MAX:
        return d - max(rank_mod_p(re, im, pp) for pp in PRIMES), True
    s = np.linalg.svd((re + 1j * im).astype(complex), compute_uv=False)
    return int(np.sum(s < FLOAT_TOL * max(1.0, s[0]))), False


def roots(n, gnum):
    sig = sum(gnum)
    return -16 * sig, (16 - 8 * n) * sig       # scaled -4*gbar and 4*gbar - 2*sigma


def band(n):
    return {(p, q) for p in range(n + 1) for q in range(n + 1)
            if abs(p - q) in (0, 2)} - {(0, 0), (n, n)}


# ---------- V1 ----------

print("V1  the band law, by full (p,q) census")
for n in (4, 5, 6, 7):
    gnum, m = locus(n), n // 2
    r0, rf = roots(n, gnum)
    at0, atf, partial, exact_reads, total_reads = set(), set(), [], 0, 0
    for p in range(n + 1):
        for q in range(n + 1):
            k0, e0 = frozen(n, h_chain(n), gnum, 768, p, q, r0)
            kf, ef = frozen(n, h_chain(n), gnum, 768, p, q, rf)
            total_reads += 2
            exact_reads += int(e0) + int(ef)
            if k0:
                at0.add((p, q))
            if kf:
                atf.add((p, q))
            if (k0 and k0 != m) or (kf and kf != m):
                partial.append((p, q, k0, kf))
    check(f"N={n} the -4gbar carriers are exactly the band |p-q| in {{0,2}} minus the 1x1 corners",
          at0 == band(n), f"{len(at0)} carriers, symmetric difference {sorted(at0 ^ band(n))}")
    check(f"N={n} every carrier carries the WHOLE floor(N/2) = {m}, never a part",
          not partial, f"partial {partial}")
    check(f"N={n} the count is 3*(N-1) = {3 * (n - 1)}", len(at0) == 3 * (n - 1), f"{len(at0)}")
    check(f"N={n} the fold root's carriers are the image under p -> N-p",
          {(n - p, q) for (p, q) in at0} == atf, f"{len(atf)} fold carriers")
    print(f"       ({exact_reads}/{total_reads} reads exact over GF(p), rest SVD nullity)")

# ---------- V2 ----------

print("\nV2  the band edge, at the N that can separate width 2 from width floor(N/2)")
for n in (8, 9, 10):
    for half, jnum, lbl in (((29, -13, 7, 19), 768, "profile A, J=3/4"),
                            ((-5, 23, -31, 3), 768, "profile B, J=3/4"),
                            ((29, -13, 7, 19), 1536, "profile A, J=3/2")):
        gnum, m = locus(n, half=half), n // 2
        r0, _ = roots(n, gnum)
        k2, _ = frozen(n, h_chain(n), gnum, jnum, 0, 2, r0)
        k4, _ = frozen(n, h_chain(n), gnum, jnum, 0, 4, r0)
        k6, _ = frozen(n, h_chain(n), gnum, jnum, 0, 6, r0)
        check(f"N={n} {lbl}: |p-q|=2 carries {m}, |p-q|=4 and 6 carry nothing "
              f"(width floor(N/2) = {m} would want 4 to carry)",
              k2 == m and k4 == 0 and k6 == 0, f"kernels {k2}/{k4}/{k6}")

# ---------- V3 ----------

print("\nV3  why the edge sits there: the recentering is tauQ-odd exactly on |A^B| = 2")
for n in (4, 5, 6, 7):
    gl = np.array(locus(n), dtype=float) / GRID
    gbar = gl.mean()
    R = lambda s: tuple(sorted(n - 1 - x for x in s))
    worst = {}
    for p in range(n + 1):
        for q in range(n + 1):
            for a in combinations(range(n), p):
                for b in combinations(range(n), q):
                    size = len(set(a) ^ set(b))
                    rate = lambda c: -2.0 * sum(gl[l] for l in set(c[0]) ^ set(c[1])) + 4.0 * gbar
                    res = abs(rate((R(b), R(a))) + rate((a, b)))
                    worst[size] = max(worst.get(size, 0.0), res)
    check(f"N={n} odd at |A^B| = 2, and NOT odd at any other size",
          worst.get(2, 1.0) < 1e-12 and all(v > 1e-9 for s, v in worst.items() if s != 2),
          ", ".join(f"|A^B|={s}: {v:.1e}" for s, v in sorted(worst.items())))
    reach = {(p, q) for p in range(n + 1) for q in range(n + 1)
             if any(len(set(a) ^ set(b)) == 2
                    for a in combinations(range(n), p) for b in combinations(range(n), q))}
    check(f"N={n} the blocks holding such a cell are exactly the band", reach == band(n),
          f"|reach| = {len(reach)}")

# ---------- V4 ----------

print("\nV4  the gate: the off-diagonal band exists iff h is bipartite; the corner does not care")
cases = []
for n in (5, 6, 7):
    cases.append((n, h_chain(n), "chain, bipartite"))
    cases.append((n, h_chain(n, [3 if a in (0, n - 1) else 0 for a in range(n)]),
                  "chain + R-invariant end diagonal"))
    cases.append((n, h_chain(n, [2 if a in (1, n - 2) else 0 for a in range(n)]),
                  "chain + R-invariant inner diagonal"))
for n in (5, 6, 7, 8):
    cases.append((n, h_ring(n), f"ring, N {'even' if n % 2 == 0 else 'odd'}"))
for n, h, lbl in cases:
    gnum, m = locus(n), n // 2
    r0, _ = roots(n, gnum)
    bip = is_bipartite_spectrum(h)
    corner, _ = frozen(n, h, gnum, 768, 1, 1, r0)
    off = [frozen(n, h, gnum, 768, p, q, r0)[0]
           for (p, q) in ((0, 2), (2, 0)) if comb(n, p) * comb(n, q) <= 700]
    want = m if bip else 0
    check(f"N={n} {lbl}: bipartite={bip}, corner {corner} = floor(N/2), off-diagonal {off} = {want}",
          corner == m and all(v == want for v in off))

# ---------- V5 ----------

print("\nV5  the falsified candidate: the chiral intertwiner does NOT transport the corner")


def c_op(j, s):
    if not (s >> j) & 1:
        return None
    return (-1.0 if bin(s & ((1 << j) - 1)).count("1") % 2 else 1.0), s ^ (1 << j)


def cdag_op(j, s):
    if (s >> j) & 1:
        return None
    return (-1.0 if bin(s & ((1 << j) - 1)).count("1") % 2 else 1.0), s | (1 << j)


def liou(n, p, q, h, gl, gbar, J=0.75):
    st = {k: [sum(1 << i for i in c) for c in combinations(range(n), k)] for k in (p, q)}
    idx = {k: {s: i for i, s in enumerate(v)} for k, v in st.items()}
    H = {}
    for k in (p, q):
        M = np.zeros((len(st[k]), len(st[k])))
        for col, s in enumerate(st[k]):
            for kk in range(n):
                a1 = c_op(kk, s)
                if a1 is None:
                    continue
                for jj in range(n):
                    if h[jj, kk] == 0:
                        continue
                    a2 = cdag_op(jj, a1[1])
                    if a2 is not None:
                        M[idx[k][a2[1]], col] += h[jj, kk] * a1[0] * a2[0]
        H[k] = M
    cells = [(a, b) for a in st[p] for b in st[q]]
    pos = {c: i for i, c in enumerate(cells)}
    d = len(cells)
    L = np.zeros((d, d), dtype=complex)
    for col, (a, b) in enumerate(cells):
        for r in range(len(st[p])):
            if H[p][r, idx[p][a]]:
                L[pos[(st[p][r], b)], col] += -1j * J * H[p][r, idx[p][a]]
        for r in range(len(st[q])):
            if H[q][idx[q][b], r]:
                L[pos[(a, st[q][r])], col] += 1j * J * H[q][idx[q][b], r]
        diff = a ^ b
        L[col, col] += -2.0 * sum(gl[l] for l in range(n) if (diff >> l) & 1) + 4.0 * gbar
    return cells, pos, L


def transport_matrix(n, ua, ub, cc, c02, p02):
    T = np.zeros((len(c02), len(cc)), dtype=complex)
    for j, (ket, bra) in enumerate(cc):
        x = ket.bit_length() - 1
        if ua[x] == 0.0:
            continue
        for m in range(n):
            if ub[m] == 0.0:
                continue
            r = cdag_op(m, bra)
            if r is not None and (0, r[1]) in p02:
                T[p02[(0, r[1])], j] += ua[x] * ub[m] * r[0]
    return T


for n in (5, 6, 7):
    h = h_chain(n).astype(float)
    gl = np.array(locus(n), dtype=float) / GRID
    gbar = gl.mean()
    lam, U = np.linalg.eigh(h)
    cc, _, L11 = liou(n, 1, 1, h, gl, gbar)
    c02, p02, L02 = liou(n, 0, 2, h, gl, gbar)
    _, _, H11 = liou(n, 1, 1, h, np.zeros(n), 0.0)
    _, _, H02 = liou(n, 0, 2, h, np.zeros(n), 0.0)
    u, s, vh = np.linalg.svd(L11)
    K11 = vh[len(s) - int(np.sum(s < 1e-9 * s[0])):].conj().T
    chiral = [(a, b) for a in range(n) for b in range(n) if abs(lam[a] + lam[b]) < 1e-9]

    # the self-check first: without it a coding slip would masquerade as the falsification
    worst_int, worst_pred = 0.0, 0.0
    for (a, b) in chiral[:4]:
        T = transport_matrix(n, U[:, a], U[:, b], cc, c02, p02)
        worst_int = max(worst_int, np.linalg.norm(H02 @ T - T @ H11))
    for (a, b) in [(0, 1), (0, 0)]:
        T = transport_matrix(n, U[:, a], U[:, b], cc, c02, p02)
        worst_pred = max(worst_pred, np.linalg.norm(
            (H02 @ T - T @ H11) - 1j * 0.75 * (lam[a] + lam[b]) * T))
    check(f"N={n} the map really is the intertwiner: [L_H, T] = 0 on a chiral pair, and equals "
          f"i*J*(eps_a+eps_b)*T off it", worst_int < 1e-10 and worst_pred < 1e-10,
          f"{worst_int:.1e} / {worst_pred:.1e}")

    worst_res, biggest_img = 0.0, 0.0
    for (a, b) in chiral[:4]:
        T = transport_matrix(n, U[:, a], U[:, b], cc, c02, p02)
        for c in range(K11.shape[1]):
            img = T @ K11[:, c]
            nrm = np.linalg.norm(img)
            if nrm > 1e-10:
                biggest_img = max(biggest_img, nrm)
                worst_res = max(worst_res, np.linalg.norm(L02 @ img) / nrm)
    check(f"N={n} and yet it does NOT land in the frozen subspace: the image is nonzero but "
          f"(L_(0,2) + 4 gbar) does not kill it", biggest_img > 1e-6 and worst_res > 1e-3,
          f"worst relative residual {worst_res:.2e}, largest image norm {biggest_img:.2e}")

print()
print("V6  the other falsified candidate: the corner's room shortage does not extend")
for n in (5, 6, 7):
    m = n // 2
    R = lambda t: tuple(sorted(n - 1 - x for x in t))
    rows = []
    for (p, q) in ((1, 1), (2, 2), (3, 3), (0, 2), (1, 3)):
        if p > n or q > n or abs(p - q) not in (0, 2) or (p, q) in {(0, 0), (n, n)}:
            continue
        parts = [(p, q)] if p == q else [(p, q), (q, p)]
        cl = [(a, b) for (x, y) in parts
              for a in combinations(range(n), x) for b in combinations(range(n), y)]
        # the shortage argument counts tauQ-fixed cells that also carry the odd recentering
        fixed2 = sum(1 for c in cl
                     if (R(c[1]), R(c[0])) == c and len(set(c[0]) ^ set(c[1])) == 2)
        predicted = fixed2 - fixed2 // 2          # surplus minus the tax that halves it
        rows.append(((p, q), len(cl), fixed2, predicted))
    offdiag_all_zero = all(pr == 0 for (pq, _, _, pr) in rows if pq[0] != pq[1])
    misses = [(pq, pr) for (pq, _, _, pr) in rows if pr != m]
    check(f"N={n} the shortage predicts ZERO in every off-diagonal band pair, where floor(N/2) "
          f"= {m} is measured", offdiag_all_zero,
          "; ".join(f"{pq} dim {d} fixed2 {f} -> predicts {pr}" for (pq, d, f, pr) in rows))
    check(f"N={n} and on the diagonal it misses too, so it is not merely incomplete",
          bool(misses), f"misses {misses}")

print()
if FAILURES:
    print(f"XY frozen band: {len(FAILURES)} FAILURE(S)")
    for f in FAILURES:
        print(f"  - {f}")
    raise SystemExit(1)
print("XY frozen band: ALL GREEN")
