# -*- coding: utf-8 -*-
"""xy_frozen_band.py

Verifier for experiments/XY_FROZEN_BAND.md: where the F140 frozen root -4*gbar lives once the
chain is XY rather than Heisenberg, and what decides it.

Prints "XY frozen band: ALL GREEN" when every check passes. Runtime about two minutes; the
opt-in "--deep" run adds minutes more. This file is a RUNNER, in the house style of the gates:
importing it executes every check.

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
  V4  the SEED gate: the OFF-DIAGONAL band exists iff the single-excitation matrix h is bipartite
      (its spectrum symmetric about 0). Four independent ways: two R-invariant diagonals added
      to the chain, three N each, and the ring, where bipartiteness is the parity of N. Every row
      is read in all three columns, because the CORNER and the DIAGONAL band survive all five,
      which is the control separating this gate from the ladder's (V8)
  V7  where the number comes from, as far as it is understood: at J = 0 every band block is
      diagonal and its cells sitting exactly at the root are precisely the |A^B| = 2 cells whose
      disagreement set is a balanced pair. Their NUMBER differs wildly from block to block, and
      every block collapses to the same floor(N/2) once the coupling turns on. In the (0,2)
      block alone the J = 0 count is already floor(N/2), nothing departs, and the frozen
      subspace converges onto the balanced-pair cells as J -> 0 (though at J = 3/4 it has
      rotated well off them, which the check pins from both sides)
  V6  the OTHER falsified candidate, as a count: tauQ fixes a cell only at B = R(A), which
      forces |A| = |B|, so the off-diagonal band pairs have no fixed cell at all and the corner's
      room shortage predicts zero where floor(N/2) is measured; on the diagonal it misses in
      both directions
  V8  WHAT FILLS THE BAND: Phi(rho) = sum_l d+_l rho d_l and Psi(rho) = sum_l d_l rho d+_l
      commute with the WHOLE Liouvillian (the dephasing cannot see a particle added at the same
      site on both sides, and the ket's +eps_a pays for the bra's -eps_a), close into an sl(2)
      with [Psi, Phi] = N - p - q, and carry the frozen subspace of a rung ONTO the one above.
      Three ladders of N-1 rungs, seeded at (1,1), (0,2) and (2,0), is the band. Includes the
      row that discriminates (a non-bipartite h keeps the DIAGONAL band, so the gate is on the
      seed, not on the ladder), the Heisenberg control, the ring, where the climber fails
      and the band is full anyway, and the PRICING rows: the ring's wrap bond and a range-2
      bond, each read twice on one and the same h, once priced by occupation (sigma+ sigma-)
      and once by order (the Jordan-Wigner sign of what is jumped over). Order-honest pricing
      makes the climber commute EXACTLY and fills every rung; occupation pricing tears it, and
      the two occupation rows then part, the ring keeping its band and the hop-over losing
      everything but its two end rungs. The two N = 4 exceptions the missing ladder leaves are
      read at five couplings on two profiles, beside the coincidence of the frozen root with
      its fold partner that leaves room for them, and beside their absence at N = 5, 6, 7
  V9  the CEILING, one N at a time: an exact GF(p) rank is an UPPER bound on the rational
      nullity (rank can only drop under the reduction), which is the direction the proved floor
      cannot supply and a float nullity never could. So an exact read of floor(N/2) certifies
      the ceiling at that point, and, the entries being polynomial in (J, gbar, delta), off a
      proper algebraic subset of them. Every rung of all three ladders at N = 5, 6, 7, the two
      sl(2) shapes as stacked-matrix nullities at every step, the room-1 free ceiling read
      two-sided (it needs offsets distinct in ABSOLUTE VALUE, and the shipped profile breaks
      that from N = 10),
      the uniform point with and without an R-invariant diagonal, and the U*(D - (p-1))
      identity entry for entry. "--deep" adds N = 8, "--deep 9" the N = 9 middle rungs
  V5  a FALSIFIED candidate, kept because it is the obvious one: Xi(rho) = d_a rho d_b on a
      chiral pair commutes with the Hamiltonian part exactly (checked to machine zero, and
      against the predicted commutator off the chiral locus) yet does NOT carry the corner's
      frozen vectors into the (0,2) frozen subspace. The dephasing is what it does not respect

Exactness: J and every gamma_l sit on the dyadic grid 1/1024, so 4*N*1024*(block entries) are
integers and the multiplicity is an exact GF(p) rank at two primes = 1 (mod 4) wherever the
block is small enough; above that the read is an SVD nullity, and V1 reports how many of its
reads were exact. Never an eigenvalue count: the departing modes crowd the root at spacing
J^(2d), which is exactly where a float spectrum lies.

A float nullity is only as good as the GAP around it, not the threshold: the N=8 ring puts a
departing mode at 7.8e-8, under the 1e-7 cut this file used to apply, and the block read 5 where
the exact rank is 4. See frozen(): the read now has to land in a moat, escalates to the exact
rank when it does not, and refuses rather than guesses when it cannot afford to.
"""
import sys
import time
from itertools import combinations
from math import comb

import numpy as np

GRID = 1024
PRIMES = (998244353, 1004535809)
EXACT_MAX = 1300
FLOAT_TOL = 1e-7
# The moat a float nullity read has to land in, and how far it is worth escalating to the exact
# rank when it does not. See frozen() for what a read inside the moat is worth and what a read
# in it is not.
ZERO_TOL = 1e-12
NONZERO_TOL = 1e-6
ESCALATE_MAX = 5000

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


def h_hopover(n):
    """The chain plus every range-2 bond, R-invariant. Read through the SPIN sector builder it
    is bare next-nearest-neighbour XX; through the fermionic one it is an order-honest range-2
    hop. Same matrix h either way: with one excitation nothing is ever jumped over."""
    h = h_chain(n)
    for a in range(n - 2):
        h[a, a + 2] = h[a + 2, a] = 1
    return h


def is_bipartite_spectrum(h):
    lam = np.sort(np.linalg.eigvalsh(h.astype(float)))
    return bool(np.allclose(lam, -lam[::-1], atol=1e-9))


# ---------- blocks ----------

def sector(n, k, h, fermi=False):
    """H restricted to the k-excitation sector: hopping from h's off-diagonal plus h's
    on-site energies. The on-site term is what an R-invariant diagonal probe rides on;
    dropping it would make V4's diagonal cases silent no-ops.

    fermi=False PRICES A HOP BY OCCUPATION: the amplitude is h[site, tgt] whatever stands in
    between, which is what sigma+ sigma- does. fermi=True prices it by ORDER instead, applying
    the Slater sign of c+_tgt c_site in the same ordered-tuple convention _ladder uses. The two
    agree on every open nearest-neighbour bond, since nothing is ever jumped over there, and
    part on the wrap bond and on any range >= 2 bond. That difference is a row of the table
    (V8 f4), not a convention: it decides whether the climber commutes."""
    st = list(combinations(range(n), k))
    idx = {s: i for i, s in enumerate(st)}
    H = np.zeros((len(st), len(st)), dtype=np.int64)
    for s, cfg in enumerate(st):
        occ = set(cfg)
        for site in cfg:
            H[s, s] += h[site, site]
            for tgt in range(n):
                if h[site, tgt] and tgt not in occ:
                    rest = occ - {site}
                    sg = (-1) ** (sum(1 for x in occ if x < site)
                                  + sum(1 for x in rest if x < tgt)) if fermi else 1
                    H[idx[tuple(sorted(rest | {tgt}))], s] += sg * h[site, tgt]
    return st, idx, H


def scaled_block(n, h, gnum, jnum, p, q, fermi=False):
    """4*N*GRID * (L_(p,q)) as Gaussian integers, in the cell basis (ket, bra)."""
    stk, ik, Hk = sector(n, p, h, fermi)
    stb, ib, Hb = sector(n, q, h, fermi)
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
    """Exact rank over GF(p), with i sent to a square root of -1. Vectorized on int64: the
    trailing submatrix is updated in one numpy operation per pivot rather than row by row, which
    is what lets the caps above reach the blocks the census used to hand to the SVD. Products
    stay below p^2 < 2^60, so nothing overflows. Rectangular input is allowed: V9 reads the
    nullity of a STACKED map that way."""
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


def frozen(n, h, gnum, jnum, p, q, root_num, fermi=False):
    """dim ker(L_(p,q) - root), the root given as its 4*N*GRID-scaled integer.

    Above EXACT_MAX the read is an SVD nullity, and a bare threshold is not enough: at N = 8
    the ring puts a departing mode at 7.8e-8, under the old 1e-7 cut, and the block read 5
    where the exact rank is 4. So the float read must land inside a MOAT (nothing between
    ZERO_TOL and NONZERO_TOL); if it does not, escalate to the exact rank while that is
    affordable, and refuse the read otherwise rather than report a number it cannot support.
    """
    re, im = scaled_block(n, h, gnum, jnum, p, q, fermi)
    d = re.shape[0]
    re = re - root_num * np.eye(d, dtype=np.int64)
    if d <= EXACT_MAX:
        return d - max(rank_mod_p(re, im, pp) for pp in PRIMES), True
    s = np.linalg.svd((re + 1j * im).astype(complex), compute_uv=False)
    s = s / max(s[0], 1.0)
    k = int(np.sum(s < ZERO_TOL))
    in_moat = not np.any((s >= ZERO_TOL) & (s <= NONZERO_TOL))
    if in_moat:
        return k, False
    if d <= ESCALATE_MAX:
        return d - max(rank_mod_p(re, im, pp) for pp in PRIMES), True
    raise RuntimeError(f"undecidable nullity read at N={n} ({p},{q}) d={d}: singular values in "
                       f"the moat {s[(s >= ZERO_TOL) & (s <= NONZERO_TOL)]}")


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
        k6, _ = frozen(n, h_chain(n), gnum, jnum, 0, 6, r0)
        # EVERY |p-q| = 4 block that fits, not just the p = 0 one: at N = 8 that is all of them
        # (5 of 5), which is the N the discriminator is actually run for.
        wide = [(p, p + 4) for p in range(n - 3) if comb(n, p) * comb(n, p + 4) <= 800]
        k4 = [frozen(n, h_chain(n), gnum, jnum, p, q, r0)[0] for (p, q) in wide]
        check(f"N={n} {lbl}: |p-q|=2 carries {m}, |p-q|=4 carries nothing in {len(wide)} of the "
              f"{n - 3} blocks, |p-q|=6 nothing (width floor(N/2) = {m} would want 4 to carry)",
              k2 == m and all(v == 0 for v in k4) and k6 == 0 and len(wide) >= 2,
              f"|p-q|=2: {k2}; |p-q|=4: {dict(zip(wide, k4))}; |p-q|=6: {k6}")

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

# The reach is not only measured. |A^B| = p + q - 2|A n B|, so a cell with |A^B| = 2 exists in
# (p,q) iff |A n B| = (p+q)/2 - 1 is an integer (p+q even) lying in [0, min(p,q)] with a union
# |A n B| + 2 that fits in N. The three bounds read |p-q| <= 2, p+q >= 2 and p+q <= 2N-2, i.e.
# the band minus (0,0) and (N,N). Enumeration cannot go past N = 7; the closed form can.
for n in (8, 12, 20, 40):
    pred = {(p, q) for p in range(n + 1) for q in range(n + 1)
            if (p + q) % 2 == 0 and abs(p - q) <= 2 and 2 <= p + q <= 2 * n - 2}
    check(f"N={n} the closed-form reach criterion is exactly the band, and gives 3(N-1) = "
          f"{3 * (n - 1)} blocks", pred == band(n) and len(pred) == 3 * (n - 1),
          f"|reach| = {len(pred)}")

# ---------- V4 ----------

print("\nV4  the SEED gate: the off-diagonal band exists iff h is bipartite, while the corner and "
      "the\n    DIAGONAL band do not care (that column is the ladder's, gated in V8)")
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
    # The DIAGONAL band is the column that separates the two gates: it must stay full in every
    # row, bipartite or not. Without it the table reads as though one gate ran the whole band.
    dia = [frozen(n, h, gnum, 768, p, p, r0)[0]
           for p in (2, 3) if comb(n, p) ** 2 <= 800]      # 800 so N=8 still reads (2,2)
    want = m if bip else 0
    check(f"N={n} {lbl}: bipartite={bip}, corner {corner} = floor(N/2), diagonal band {dia} = "
          f"floor(N/2), off-diagonal {off} = {want}",
          corner == m and all(v == want for v in off)
          and bool(dia) and all(v == m for v in dia))

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
    """The map the note calls Xi; not the proof's T, which is its bordered cofactor."""
    Xi = np.zeros((len(c02), len(cc)), dtype=complex)
    for j, (ket, bra) in enumerate(cc):
        x = ket.bit_length() - 1
        if ua[x] == 0.0:
            continue
        for m in range(n):
            if ub[m] == 0.0:
                continue
            r = cdag_op(m, bra)
            if r is not None and (0, r[1]) in p02:
                Xi[p02[(0, r[1])], j] += ua[x] * ub[m] * r[0]
    return Xi


for n in (5, 6, 7):
    h = h_chain(n).astype(float)
    gl = np.array(locus(n), dtype=float) / GRID
    gbar = gl.mean()
    lam, U = np.linalg.eigh(h)
    cc, _, L11 = liou(n, 1, 1, h, gl, gbar)
    c02, p02, L02 = liou(n, 0, 2, h, gl, gbar)
    _, _, H11 = liou(n, 1, 1, h, np.zeros(n), 0.0)
    _, _, H02 = liou(n, 0, 2, h, np.zeros(n), 0.0)
    # the corner's frozen basis: DIMENSION from the exact read, vectors from the SVD, which is
    # the discipline frozen_basis() follows and the one this block used to bypass.
    u, s, vh = np.linalg.svd(L11)
    k11 = frozen(n, h_chain(n), locus(n), 768, 1, 1, roots(n, locus(n))[0])[0]
    K11 = vh[len(s) - k11:].conj().T
    chiral = [(a, b) for a in range(n) for b in range(n) if abs(lam[a] + lam[b]) < 1e-9]

    # the self-check first: without it a coding slip would masquerade as the falsification
    worst_int, worst_pred = 0.0, 0.0
    for (a, b) in chiral[:4]:
        Xi = transport_matrix(n, U[:, a], U[:, b], cc, c02, p02)
        worst_int = max(worst_int, np.linalg.norm(H02 @ Xi - Xi @ H11))
    for (a, b) in [(0, 1), (0, 0)]:
        Xi = transport_matrix(n, U[:, a], U[:, b], cc, c02, p02)
        worst_pred = max(worst_pred, np.linalg.norm(
            (H02 @ Xi - Xi @ H11) - 1j * 0.75 * (lam[a] + lam[b]) * Xi))
    check(f"N={n} the map really is the intertwiner: [L_H, Xi] = 0 on a chiral pair, and equals "
          f"i*J*(eps_a+eps_b)*Xi off it", worst_int < 1e-10 and worst_pred < 1e-10,
          f"{worst_int:.1e} / {worst_pred:.1e}")

    worst_res, biggest_img = 0.0, 0.0
    span = []
    for (a, b) in chiral[:4]:
        Xi = transport_matrix(n, U[:, a], U[:, b], cc, c02, p02)
        for c in range(K11.shape[1]):
            img = Xi @ K11[:, c]
            nrm = np.linalg.norm(img)
            if nrm > 1e-10:
                biggest_img = max(biggest_img, nrm)
                worst_res = max(worst_res, np.linalg.norm(L02 @ img) / nrm)
                span.append(img / nrm)
    # how big is the image span, against the frozen subspace of dimension floor(N/2) it was
    # supposed to land in? A hand-read number in the note until it was gated.
    sv = np.linalg.svd(np.array(span).T, compute_uv=False)
    img_rank = int(np.sum(sv > 1e-9 * sv[0]))
    check(f"N={n} the images span far more than the floor(N/2) = {n // 2} dimensions of the "
          f"target frozen subspace", img_rank > n // 2,
          f"image span rank {img_rank} in the {comb(n, 2)}-dimensional (0,2) block")
    check(f"N={n} and yet it does NOT land in the frozen subspace: the image is nonzero but "
          f"(L_(0,2) + 4 gbar) does not kill it", biggest_img > 1e-6 and worst_res > 1e-3,
          f"worst relative residual {worst_res:.2e}, largest image norm {biggest_img:.2e}")

print()
print("V7  where the number comes from, as far as it is understood")
for n in (5, 6, 7):
    gnum, m = locus(n), n // 2
    sig = sum(gnum)
    r0, _ = roots(n, gnum)
    # (a) at J = 0 the block is diagonal; which cells sit exactly at the root, and how many?
    counts, sizes_seen = {}, set()
    for (p, q) in ((0, 2), (1, 1), (2, 2), (1, 3), (3, 3)):
        if p > n or q > n or abs(p - q) not in (0, 2) or (p, q) in {(0, 0), (n, n)}:
            continue
        zero = [(a, b) for a in combinations(range(n), p) for b in combinations(range(n), q)
                if 2 * n * sum(gnum[l] for l in set(a) ^ set(b)) == 4 * sig]
        counts[(p, q)] = len(zero)
        sizes_seen |= {len(set(a) ^ set(b)) for (a, b) in zero}
        # every such cell's disagreement set must be a balanced pair {l, R(l)}
        bad = [(a, b) for (a, b) in zero
               if sorted(set(a) ^ set(b)) != sorted({min(set(a) ^ set(b)),
                                                     n - 1 - min(set(a) ^ set(b))})]
        if bad:
            counts[(p, q)] = -1
    check(f"N={n} at J = 0 the cells on the root are exactly the balanced-pair cells, and only "
          f"|A^B| = 2 ones", all(v >= 0 for v in counts.values()) and sizes_seen == {2},
          f"sizes {sorted(sizes_seen)}")
    # Both operands here are J = 0 quantities, so this check is a J = 0 statement and its name
    # says so: the collapse ITSELF, at J > 0, is V1's census and is not re-read here.
    check(f"N={n} at J = 0 their number differs block to block, while the (0,2) seed already "
          f"sits at floor(N/2) = {m} (the collapse of the others is V1's census, not this read)",
          len(set(counts.values())) > 1 and counts[(0, 2)] == m,
          "; ".join(f"{pq}: {v}" for pq, v in counts.items()))

    # (b) the (0,2) block: how does its frozen subspace sit against the pair cells, as J varies?
    # "Continuation" can only be a J -> 0 statement, and the two-sided form matters: at the
    # coupling this note reads everywhere else the subspace has rotated almost entirely OFF the
    # pair cells, so a bare "overlap is nonsingular" would pass on nearly orthogonal subspaces.
    def ker(jnum):
        re, im = scaled_block(n, h_chain(n), gnum, jnum, 0, 2)
        d = re.shape[0]
        M = (re + 1j * im).astype(complex) - r0 * np.eye(d)
        _, sv, vh = np.linalg.svd(M)
        return vh[d - int(np.sum(sv < 1e-9 * max(1.0, sv[0]))):].conj().T
    cells02 = list(combinations(range(n), 2))
    Ppair = np.zeros((len(cells02), m))
    for j, i in enumerate([i for i, c in enumerate(cells02) if c[1] == n - 1 - c[0]]):
        Ppair[i, j] = 1.0
    Kbig, Ksmall = ker(768), ker(4)          # J = 3/4 and J = 1/256
    cos_big = np.linalg.svd(Ppair.T @ Kbig, compute_uv=False)
    cos_small = np.linalg.svd(Ppair.T @ Ksmall, compute_uv=False)
    check(f"N={n} in (0,2) the frozen dimension is floor(N/2) = {m} at every coupling, and the "
          f"subspace CONVERGES onto the pair cells as J -> 0 while sitting well off them at "
          f"J = 3/4",
          Kbig.shape[1] == m and Ksmall.shape[1] == m
          and cos_small.min() > 0.9 and cos_big.min() < 0.9,
          f"principal cosines with the pair-cell span: J=3/4 "
          f"{np.array2string(cos_big, precision=4)}, J=1/256 "
          f"{np.array2string(cos_small, precision=4)}")

    # (c) The J = 0 census above is a statement about the PROFILE, not about the locus. A locus
    # profile whose pair offsets sum to 2*gbar puts a FOUR-site set on the root: the note used to
    # say that cannot happen on the locus, and it can. Both halves are asserted, since the point
    # is that the extra root cells exist AND that the band law survives them.
# (c2) the one diagonal: A^B is not a new object, it is the rung k = popcount(i^j) of the single
# dephasing diagonal Q (rate -2*gamma*k). Three consequences, all exact: the frozen root IS the
# k = 2 rung and the locus is what puts a BALANCED pair back on it while an unbalanced pair falls
# off; k and |p-q| always share a parity; and Phi leaves k exactly fixed, which is the whole
# reason the dephasing commutes with it.
for n in (4, 5, 6):
    gnum, gbar = locus(n), sum(locus(n)) / n
    pr = [(l, n - 1 - l) for l in range(n // 2)]
    on_rung = all(abs(-2 * sum(gnum[l] for l in p) + 4 * gbar) < 1e-9 for p in pr)
    off = [s for s in combinations(range(n), 2) if frozenset(s) not in {frozenset(p) for p in pr}]
    off_rung = all(abs(-2 * sum(gnum[l] for l in s) + 4 * gbar) > 1e-9 for s in off)
    par = sum(1 for p in range(n + 1) for q in range(n + 1)
              for a in combinations(range(n), p) for b in combinations(range(n), q)
              if (len(set(a) ^ set(b)) - abs(p - q)) % 2)
    kmoved = sum(1 for p in range(n) for q in range(n)
                 for a in combinations(range(n), p) for b in combinations(range(n), q)
                 for l in range(n) if l not in a and l not in b
                 and len(set(a + (l,)) ^ set(b + (l,))) != len(set(a) ^ set(b)))
    # the builder's own rate diagonal against the Absorption Theorem's CARRIER-VECTOR form,
    # -Re(lambda) = 2*sum_x gamma_x*<Delta_x>: the band note uses that law, it does not define
    # a rate of its own. Checked against the block this file actually builds.
    absorb = 0
    for (p, q) in ((1, 1), (0, 2), (2, 2)):
        re, im = scaled_block(n, h_chain(n), gnum, 768, p, q)
        stk = list(combinations(range(n), p))
        stb = list(combinations(range(n), q))
        for a, ka in enumerate(stk):
            for b, kb in enumerate(stb):
                col = a * len(stb) + b
                want = -2 * 4 * n * sum(gnum[x] * (1 if x in (set(ka) ^ set(kb)) else 0)
                                        for x in range(n))
                absorb = max(absorb, abs(int(re[col, col]) - want))
    # and the EVEN rail is not the number-conserving rail: depth 2 also carries |p-q| = 2, which
    # is where the off-diagonal ladders live. (Corrects reflections/THE_VIEW_ONTO_THE_MEMORY.md.)
    even_dn = set()
    for p in range(n + 1):
        for q in range(n + 1):
            for a in combinations(range(n), p):
                for b in combinations(range(n), q):
                    if len(set(a) ^ set(b)) % 2 == 0:
                        even_dn.add(abs(p - q))
    # Three of the six operands below (par, kmoved, even_dn) are SET identities, true for every
    # N and every profile, and the name says which they are. They are here because they are the
    # combinatorial half of the ladder's mechanism, not because they witness the operator: the
    # operator-level statement, that Phi commutes with L, is V8 (a)'s exact zero.
    check(f"N={n} the root is the k=2 rung of the one diagonal (balanced pairs ON it, unbalanced "
          f"two-sets OFF it), the rate diagonal IS the gamma-vector Absorption law, and, by set "
          f"algebra alone, k and |p-q| share a parity, adding a site to both indices leaves k "
          f"fixed, and the EVEN rail carries |p-q| = 2 too",
          on_rung and off_rung and par == 0 and kmoved == 0 and absorb == 0
          and {0, 2} <= even_dn,
          f"rung -4*gbar = {-4*gbar/GRID:.4f}; a balanced pair pays "
          f"{-2*sum(gnum[l] for l in pr[0])/GRID:.4f}, the unbalanced two-set {off[0]} pays "
          f"{-2*sum(gnum[l] for l in off[0])/GRID:.4f}; parity misses {par}, "
          f"Phi images moving k {kmoved}, Absorption mismatch {absorb}, "
          f"|p-q| seen on the even rail {sorted(even_dn)}")

# (d) the third structure in the (0,2) block, stated in the note and until now only asserted: the
# involution that makes M odd there is ANTILINEAR. The block's ket sector is the vacuum, so its
# Hamiltonian part is +i*J*H_2 alone; theta = (site reversal) o (complex conjugation) turns that
# sign, and on the locus reversal turns the rate diagonal's sign too, at |A^B| = 2. Two-sided: it
# must ALSO survive the diagonal perturbation that empties the band, which is why it is not the gate.
for n in (5, 6, 7):
    gnum, m = locus(n), n // 2
    r0, _ = roots(n, gnum)
    cells = list(combinations(range(n), 2))
    perm = [cells.index(tuple(sorted(n - 1 - x for x in c))) for c in cells]
    rows = []
    for lbl, h in (("chain", h_chain(n)),
                   ("chain + end diagonal", h_chain(n, [3 if a in (0, n - 1) else 0
                                                        for a in range(n)]))):
        re02, im02 = scaled_block(n, h, gnum, 768, 0, 2)
        M = (re02 + 1j * im02).astype(complex) - r0 * np.eye(len(cells))
        th = np.conj(M)[np.ix_(perm, perm)]          # theta M theta, theta antilinear
        k = frozen(n, h, gnum, 768, 0, 2, r0)[0]
        rows.append((lbl, np.abs(th + M).max() / np.abs(M).max(), k))
    check(f"N={n} in (0,2) the odd involution is ANTILINEAR, and it survives the perturbation that "
          f"empties the band, so it is neither the counter nor the gate",
          all(r < 1e-12 for _, r, _ in rows) and rows[0][2] == m and rows[1][2] == 0,
          "; ".join(f"{lbl}: |theta M theta + M| = {r:.1e}, frozen {k}" for lbl, r, k in rows))

nsp = 8
gs = [0] * nsp
for i, d in enumerate((60, 55, 45, 24)):              # distinct, sum = 184 = 2*gbar
    gs[i], gs[nsp - 1 - i] = 92 + d, 92 - d
sigs = sum(gs)
by_size = {}
for (p, q) in ((0, 4), (0, 2)):
    for a in combinations(range(nsp), p):
        for b in combinations(range(nsp), q):
            s = set(a) ^ set(b)
            if 2 * nsp * sum(gs[l] for l in s) == 4 * sigs:
                by_size[len(s)] = by_size.get(len(s), 0) + 1
rs, _ = roots(nsp, gs)
band_k = [frozen(nsp, h_chain(nsp), gs, 768, p, q, rs)[0] for (p, q) in ((1, 1), (0, 2))]
wide_k = [frozen(nsp, h_chain(nsp), gs, 768, p, q, rs)[0] for (p, q) in ((0, 4), (1, 5))]
check(f"N={nsp} on a locus profile whose pair offsets sum to 2*gbar a FOUR-site set sits on the "
      f"root at J = 0, and the band law survives it anyway",
      by_size.get(4, 0) > 0 and all(v == nsp // 2 for v in band_k)
      and all(v == 0 for v in wide_k),
      f"J=0 root cells by |A^B| {dict(sorted(by_size.items()))}; at J=3/4 the band carries "
      f"{band_k} and |p-q|=4 carries {wide_k}")

print()
print("V6  the other falsified candidate: the corner's room shortage does not extend")
for n in (5, 6, 7):
    m = n // 2
    gnum6 = locus(n)
    r6, _ = roots(n, gnum6)
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
    # The prediction side is structurally zero off the diagonal (tauQ fixes a cell only at
    # |A| = |B|), so on its own it could not fail. The contrast the check is named for is
    # therefore READ here rather than referred to: the same blocks' multiplicities, measured.
    offdiag_all_zero = all(pr == 0 for (pq, _, _, pr) in rows if pq[0] != pq[1])
    measured = [((p, q), frozen(n, h_chain(n), gnum6, 768, p, q, r6)[0])
                for (p, q) in ((0, 2), (1, 3))]
    check(f"N={n} the shortage predicts ZERO in every off-diagonal band pair, and floor(N/2) "
          f"= {m} is what those same blocks measure", offdiag_all_zero
          and all(v == m for _, v in measured),
          "; ".join(f"{pq} dim {d} fixed2 {f} -> predicts {pr}" for (pq, d, f, pr) in rows)
          + f"; measured {measured}")
    # The diagonal half, scoped to the diagonal ONLY. Including the off-diagonal pairs here
    # would make the check vacuous, since their prediction is 0 != m by construction.
    diag = [(pq, pr) for (pq, _, _, pr) in rows if pq[0] == pq[1]]
    diag_misses = [(pq, pr) for (pq, pr) in diag if pr != m]
    want_miss = {5: False, 6: True, 7: True}[n]
    check(f"N={n} on the diagonal the shortage "
          f"{'misses, and in BOTH directions' if want_miss else 'happens to MATCH everywhere, so it is not always wrong'}",
          bool(diag_misses) == want_miss
          and (n != 6 or sorted(diag_misses) == [((2, 2), 0), ((3, 3), 6)]),
          f"diagonal predictions {diag}, misses {diag_misses}")

print()
print("V8  the collapse: an exact ladder carries the seed's count along the band")

# Phi(rho) = sum_l d+_l rho d_l and Psi(rho) = sum_l d_l rho d+_l, in the cell basis. Adding
# (removing) one particle at the SAME site on both sides leaves the disagreement set A^B alone,
# so the dephasing diagonal is untouched cell by cell; and sum over SITES = sum over free-fermion
# EIGENMODES, each of which intertwines the Hamiltonian part because the ket's +eps_a and the
# bra's -eps_a cancel. Hence both commute with the whole Liouvillian, for ANY gamma profile.


def _ladder(n, p, q, up):
    """Phi: (p,q) -> (p+1,q+1) when up, else Psi: (p,q) -> (p-1,q-1). Jordan-Wigner signs."""
    s = 1 if up else -1
    stk = list(combinations(range(n), p))
    stb = list(combinations(range(n), q))
    stk2 = list(combinations(range(n), p + s)); ik2 = {t: i for i, t in enumerate(stk2)}
    stb2 = list(combinations(range(n), q + s)); ib2 = {t: i for i, t in enumerate(stb2)}
    M = np.zeros((len(stk2) * len(stb2), len(stk) * len(stb)))
    for a, ka in enumerate(stk):
        for b, kb in enumerate(stb):
            for l in range(n):
                if up and (l in ka or l in kb):
                    continue
                if not up and (l not in ka or l not in kb):
                    continue
                A2 = tuple(sorted(ka + (l,))) if up else tuple(x for x in ka if x != l)
                B2 = tuple(sorted(kb + (l,))) if up else tuple(x for x in kb if x != l)
                sg = (-1) ** (sum(1 for x in A2 if x < l) + sum(1 for x in B2 if x < l))
                M[ik2[A2] * len(stb2) + ib2[B2], a * len(stb) + b] += sg
    return M


def _dense(n, h, gnum, jnum, p, q, fermi=False):
    re, im = scaled_block(n, h, gnum, jnum, p, q, fermi)
    return (re + 1j * im).astype(complex)


def _heis(n, gnum, jnum, p, q):
    """The XXZ block (XX+YY+ZZ), the control: its H is quartic, so no ladder can intertwine."""
    def sec(k):
        st = list(combinations(range(n), k))
        idx = {t: i for i, t in enumerate(st)}
        H = np.zeros((len(st), len(st)), dtype=np.int64)
        for s, cfg in enumerate(st):
            occ = set(cfg)
            for b in range(n - 1):
                H[s, s] += (1 if b in occ else -1) * (1 if b + 1 in occ else -1)
                if (b in occ) != (b + 1 in occ):
                    H[idx[tuple(sorted(occ ^ {b, b + 1}))], s] += 2
        return st, H
    stk, Hk = sec(p)
    stb, Hb = sec(q)
    f, nb = 4 * n, len(stb)
    d = len(stk) * nb
    re = np.zeros((d, d), dtype=np.int64)
    im = np.zeros((d, d), dtype=np.int64)
    for a, ka in enumerate(stk):
        for b, kb in enumerate(stb):
            col = a * nb + b
            for c in range(len(stk)):
                if Hk[c, a]:
                    im[c * nb + b, col] += -f * jnum * Hk[c, a]
            for c in range(nb):
                if Hb[b, c]:
                    im[a * nb + c, col] += f * jnum * Hb[b, c]
            re[col, col] += -2 * f * sum(gnum[l] for l in set(ka) ^ set(kb))
    return re, im


def frozen_basis(n, h, gnum, jnum, p, q, root_num, fermi=False):
    """A basis of ker(L_(p,q) - root), its DIMENSION fixed by the exact/SVD read of frozen()."""
    k = frozen(n, h, gnum, jnum, p, q, root_num, fermi)[0]
    d = comb(n, p) * comb(n, q)
    if k == 0:
        return np.zeros((d, 0), dtype=complex)
    M = _dense(n, h, gnum, jnum, p, q, fermi) - root_num * np.eye(d)
    return np.linalg.svd(M)[2][d - k:].conj().T


def cosines(A, B):
    if A.shape[1] == 0 or B.shape[1] == 0:
        return np.array([0.0])
    return np.linalg.svd(np.linalg.qr(A)[0].conj().T @ np.linalg.qr(B)[0],
                         compute_uv=False)


# (a) the intertwining itself, exactly (integer matrices, so an exact zero is available), on the
# band and OFF it; against the Heisenberg control, where it must visibly fail.
for n in (5, 6):
    gnum, jn = locus(n), 768
    worst_up = worst_dn = 0.0
    # Both an h with a symmetric spectrum and one without: the climber asks that H be QUADRATIC,
    # not that h be bipartite, so a non-bipartite row has to intertwine just as exactly.
    for h in (h_chain(n), h_chain(n, [3 if a in (0, n - 1) else 0 for a in range(n)])):
        for (p, q) in ((1, 1), (0, 2), (2, 2), (1, 3), (1, 2), (0, 3)):
            if p + 1 > n or q + 1 > n:
                continue
            U = _ladder(n, p, q, True)
            worst_up = max(worst_up, np.abs(_dense(n, h, gnum, jn, p + 1, q + 1) @ U
                                            - U @ _dense(n, h, gnum, jn, p, q)).max())
            D = _ladder(n, p + 1, q + 1, False)
            worst_dn = max(worst_dn, np.abs(_dense(n, h, gnum, jn, p, q) @ D
                                            - D @ _dense(n, h, gnum, jn, p + 1, q + 1)).max())
    def _hc(p, q):
        re, im = _heis(n, gnum, jn, p, q)
        return (re + 1j * im).astype(complex)
    hz = [np.abs(_hc(p + 1, q + 1) @ _ladder(n, p, q, True)
                 - _ladder(n, p, q, True) @ _hc(p, q)).max()
          / np.abs(_hc(p, q)).max() for (p, q) in ((1, 1), (0, 2))]
    check(f"N={n} Phi and Psi commute with L EXACTLY on XY (band and off-band alike, bipartite "
          f"h or not), and visibly do NOT on Heisenberg",
          worst_up == 0.0 and worst_dn == 0.0 and min(hz) > 1e-2,
          f"XY max |[L, Phi]| = {worst_up:.1e}, |[L, Psi]| = {worst_dn:.1e}; "
          f"Heisenberg relative {min(hz):.2e}")
        # (a2) how far "quadratic" may be weakened, answered with the one term the BLOCK picture
    # cannot hold: a pairing term c+_a c+_{a+1} is quadratic and does NOT conserve number, so it
    # has no (p,q) block at all. Read in the full 2^N x 2^N Fock space instead, no blocks. The
    # condition on the Hamiltonian is number-conserving quadratic, not quadratic, and the
    # dephasing half of the argument is read here too, where it must hold with no condition.
    if n == 5:
        NF, gf, jf = 4, [0.13, 0.29, 0.07, 0.19], 0.7

        def _cdag(l):
            M = np.zeros((2 ** NF, 2 ** NF))
            for s in range(2 ** NF):
                if not (s >> l) & 1:
                    M[s | (1 << l), s] = (-1) ** sum(1 for x in range(l) if (s >> x) & 1)
            return M
        CD = [_cdag(l) for l in range(NF)]
        CC = [m.T for m in CD]
        mask = np.array([[-2 * sum(gf[l] for l in range(NF)
                                   if ((a >> l) & 1) != ((b >> l) & 1))
                          for b in range(2 ** NF)] for a in range(2 ** NF)])

        def _phi(r):
            return sum(CD[l] @ r @ CC[l] for l in range(NF))
        rng = np.random.default_rng(7)

        def _res(HF):
            w = 0.0
            for _ in range(6):
                r = (rng.normal(size=(2 ** NF,) * 2)
                     + 1j * rng.normal(size=(2 ** NF,) * 2))
                lhs = -1j * (HF @ _phi(r) - _phi(r) @ HF) + mask * _phi(r)
                rhs = _phi(-1j * (HF @ r - r @ HF) + mask * r)
                w = max(w, np.abs(lhs - rhs).max() / max(np.abs(rhs).max(), 1.0))
            return w

        def _quad(hm):
            return sum(hm[a, b] * (CD[a] @ CC[b]) for a in range(NF) for b in range(NF))
        # The reason is an index relabelling and asks NOTHING of h, so h is varied all the way
        # to non-Hermitian: the eigenmode picture the note also gives would need normality here.
        hs = rng.normal(size=(NF, NF))
        hh = rng.normal(size=(NF, NF)) + 1j * rng.normal(size=(NF, NF))
        variants = {"real symmetric": _quad(hs + hs.T),
                    "complex Hermitian": _quad(hh + hh.conj().T),
                    "non-Hermitian": _quad(hh)}
        pair = sum(jf * (CD[a] @ CC[a + 1] + CD[a + 1] @ CC[a]) for a in range(NF - 1)) \
            + sum(0.4 * (CD[a] @ CD[a + 1] + CC[a + 1] @ CC[a]) for a in range(NF - 1))
        good = {k: _res(v) for k, v in variants.items()}
        bad = _res(pair)
        nonnormal = np.abs(variants["non-Hermitian"] @ variants["non-Hermitian"].conj().T
                           - variants["non-Hermitian"].conj().T
                           @ variants["non-Hermitian"]).max()
        dep = max(np.abs(mask * _phi(r) - _phi(mask * r)).max() for r in
                  [rng.normal(size=(2 ** NF,) * 2) + 1j * rng.normal(size=(2 ** NF,) * 2)
                   for _ in range(6)])
        check("N=4 in the FULL Fock space: the condition is number-conserving quadratic and NOT "
              "quadratic (a pairing term is quadratic and breaks the climber), it asks nothing "
              "of h (a non-Hermitian, non-normal h commutes just as exactly), and the dephasing "
              "half holds with no condition at all",
              max(good.values()) < 1e-12 < bad and dep < 1e-12 and nonnormal > 1e-2,
              f"{ {k: f'{v:.1e}' for k, v in good.items()} }, with pairing {bad:.1e}, "
              f"dephasing alone {dep:.1e}, non-normality of the third h {nonnormal:.1e}")

    # (b) the sl(2) relation that makes the ladder finite
    bad = 0.0
    for (p, q) in ((1, 1), (2, 2), (1, 3)):
        if p + 1 > n or q + 1 > n:
            continue
        C = _ladder(n, p + 1, q + 1, False) @ _ladder(n, p, q, True) \
            - _ladder(n, p - 1, q - 1, True) @ _ladder(n, p, q, False)
        bad = max(bad, np.abs(C - (n - p - q) * np.eye(C.shape[0])).max())
    check(f"N={n} [Psi, Phi] = (N - p - q) on the (p,q) block, so the ladder is an sl(2) and "
          f"must terminate", bad == 0.0, f"max deviation {bad:.1e}")

# (c) the transport: Phi carries the frozen subspace ONTO the one above, not merely into it.
# The control at the end is what stops this passing on a subspace that merely fails to be
# orthogonal: the same climber applied to a NON-frozen subspace must not land on F1.
for n, jn, jl in ((5, 768, "3/4"), (6, 768, "3/4"), (6, 1536, "3/2"), (7, 768, "3/4")):
    h, gnum, m = h_chain(n), locus(n), n // 2
    r0, _ = roots(n, gnum)
    rows, ok = [], True
    for (p, q) in ((1, 1), (0, 2)):
        if comb(n, p + 1) * comb(n, q + 1) > 520:
            continue
        F0 = frozen_basis(n, h, gnum, jn, p, q, r0)
        F1 = frozen_basis(n, h, gnum, jn, p + 1, q + 1, r0)
        img = _ladder(n, p, q, True) @ F0
        c = cosines(F1, img)
        d = comb(n, p) * comb(n, q)
        junk = np.linalg.svd(_dense(n, h, gnum, jn, p, q)
                             - r0 * np.eye(d))[2][:F0.shape[1]].conj().T
        cj = cosines(F1, _ladder(n, p, q, True) @ junk)
        rows.append(f"({p},{q})->({p+1},{q+1}) {F0.shape[1]}->{F1.shape[1]}, "
                    f"min cos {c.min():.9f}, control {cj.max():.3f}")
        ok = ok and F0.shape[1] == m and F1.shape[1] == m and c.min() > 1 - 1e-7 \
            and cj.max() < 0.9
    # len(rows) is asserted so that a size filter which happened to skip every case would FAIL
    # here rather than pass on an empty loop, the way (d) below asserts its list lengths.
    check(f"N={n} J={jl} Phi maps the frozen subspace ONTO the one above (and the same map "
          f"applied off it does not), on both ladder seeds",
          ok and len(rows) == 2, "; ".join(rows))

# (d) the ladder's ends and its full length: this is where 3(N-1) comes from.
for n in (5, 6):
    h, gnum, m = h_chain(n), locus(n), n // 2
    r0, _ = roots(n, gnum)
    diag = [frozen(n, h, gnum, 768, p, p, r0)[0] for p in range(n + 1)
            if comb(n, p) ** 2 <= 520]
    off = [frozen(n, h, gnum, 768, p, p + 2, r0)[0] for p in range(n - 1)
           if comb(n, p) * comb(n, p + 2) <= 520]
    top = frozen_basis(n, h, gnum, 768, n - 1, n - 1, r0)
    seed = frozen_basis(n, h, gnum, 768, 1, 1, r0)
    kill_top = np.abs(_ladder(n, n - 1, n - 1, True) @ top).max()
    kill_seed = np.abs(_ladder(n, 1, 1, False) @ seed).max()
    check(f"N={n} the diagonal ladder runs p = 1..N-1 and the off-diagonal p = 0..N-2, every "
          f"rung floor(N/2) = {m}, both ends empty: that is the 3(N-1) count",
          diag[0] == 0 and diag[-1] == 0 and all(v == m for v in diag[1:-1])
          and all(v == m for v in off) and len(diag) == n + 1 and len(off) == n - 1,
          f"diagonal {diag}, off-diagonal {off}")
    check(f"N={n} Phi annihilates the TOP rung and Psi the seed, so each multiplet has exactly "
          f"N-1 rungs", kill_top < 1e-9 and kill_seed < 1e-9 and top.shape[1] == m,
          f"|Phi F(N-1,N-1)| = {kill_top:.1e}, |Psi F(1,1)| = {kill_seed:.1e}")

# (e) THE DISCRIMINATOR. Phi needs h quadratic, not h bipartite. So on a chain carrying an
# R-invariant diagonal, where V4 measures the off-diagonal band EMPTY, the DIAGONAL ladder must
# still be FULL. Both halves are asserted, since either alone would pass on a wrong picture.
for n in (5, 6, 7):
    gnum, m = locus(n), n // 2
    r0, _ = roots(n, gnum)
    for lbl, dg in (("end", [3 if a in (0, n - 1) else 0 for a in range(n)]),
                    ("inner", [2 if a in (1, n - 2) else 0 for a in range(n)])):
        h = h_chain(n, dg)
        dd = [(p, frozen(n, h, gnum, 768, p, p, r0)[0]) for p in (1, 2, 3)
              if comb(n, p) ** 2 <= 520]
        oo = [(p, frozen(n, h, gnum, 768, p, p + 2, r0)[0]) for p in (0, 1)
              if comb(n, p) * comb(n, p + 2) <= 520]
        check(f"N={n} chain + R-invariant {lbl} diagonal (h NOT bipartite): the diagonal ladder "
              f"stays FULL at floor(N/2) = {m} while the off-diagonal one is empty",
              not is_bipartite_spectrum(h) and all(v == m for _, v in dd)
              and all(v == 0 for _, v in oo) and len(dd) >= 2 and len(oo) >= 1,
              f"diagonal {dd}, off-diagonal {oo}")

# (f) where the J = 0 numbers come from: floor(N/2) pairs times the room beside the pair. The
# room is what collapses; the pair count is what survives.
for n in range(4, 10):
    m, bad = n // 2, []
    for (p, q) in sorted(band(n)):
        if p > q:
            continue
        per = {}
        for a in combinations(range(n), p):
            for b in combinations(range(n), q):
                s = set(a) ^ set(b)
                if len(s) == 2 and sorted(s) == sorted({min(s), n - 1 - min(s)}):
                    per[min(s)] = per.get(min(s), 0) + 1
        want = 2 * comb(n - 2, p - 1) if p == q else comb(n - 2, min(p, q))
        if len(per) != m or any(v != want for v in per.values()):
            bad.append(((p, q), per, want))
    check(f"N={n} the J = 0 root cells factor as floor(N/2) pairs x room(p,q), room = "
          f"2*C(N-2,p-1) on the diagonal and C(N-2,p) off it", not bad, f"misses {bad[:3]}")

# (f2) the Heisenberg row of the table, from an independently built block: the census F140
# already knows, which is what makes the failed intertwining above a control and not a bug.
for n in (5, 6):
    gnum, m = locus(n), n // 2
    r0, _ = roots(n, gnum)
    rows = []
    # (n-1, n-1) is in the list because it is F140's SECOND corner: leaving it out would let
    # "everything but the corner is empty" pass while the top rung quietly carried floor(N/2),
    # and the note's table reads that rung in its diagonal-band column.
    for (p, q) in ((1, 1), (0, 2), (2, 2), (1, 3), (n - 1, n - 1)):
        if comb(n, p) * comb(n, q) > 520:
            continue
        # exact, like every other multiplicity here: the Heisenberg blocks are integer too, so
        # there is no reason for the control to be the one read taken on a float threshold.
        re, im = _heis(n, gnum, 768, p, q)
        d = comb(n, p) * comb(n, q)
        re = re - r0 * np.eye(d, dtype=np.int64)
        rows.append(((p, q), d - max(rank_mod_p(re, im, pp) for pp in PRIMES)))
    check(f"N={n} Heisenberg: BOTH corners carry floor(N/2) = {m} and every band block between "
          f"them carries nothing, so the row that kills the ladder is the row that empties the "
          f"band down to F140's two corners",
          dict(rows)[(1, 1)] == m and dict(rows)[(n - 1, n - 1)] == m
          and all(v == 0 for pq, v in rows if pq not in ((1, 1), (n - 1, n - 1))),
          f"{rows}")

# (f3) the honest limit: on the RING the climber does not commute (the Jordan-Wigner string
# wraps, so Phi moves between the periodic and antiperiodic fermion sectors) and does not carry
# the corner up, yet the diagonal band is full. Sufficient, not known to be necessary.
for n in (5, 6, 7):
    gnum, m = locus(n), n // 2
    r0, _ = roots(n, gnum)
    h = h_ring(n)
    U = _ladder(n, 1, 1, True)
    Unos = np.abs(_ladder(n, 1, 1, True))        # the same climber without the JW string signs
    res = min(np.abs(_dense(n, h, gnum, 768, 2, 2) @ V - V @ _dense(n, h, gnum, 768, 1, 1)).max()
              / max(np.abs(_dense(n, h, gnum, 768, 2, 2)).max(), 1.0) for V in (U, Unos))
    F1 = frozen_basis(n, h, gnum, 768, 1, 1, r0)
    F2 = frozen_basis(n, h, gnum, 768, 2, 2, r0)
    c = cosines(F2, _ladder(n, 1, 1, True) @ F1)
    # every diagonal rung the ring can afford, not just the first two: "the band is full there"
    # is the whole point of the exception, so it has to be read past rung 2.
    rungs = [(p, frozen(n, h, gnum, 768, p, p, r0)[0]) for p in range(1, n)
             if comb(n, p) ** 2 <= 800]
    check(f"N={n} ring: the climber does NOT commute and does NOT carry the corner ONTO the rung "
          f"above, yet every diagonal rung read is full at floor(N/2) = {m}: the ladder is "
          f"sufficient, not necessary",
          res > 1e-2 and c.min() < 0.9 and all(v == m for _, v in rungs) and len(rungs) >= 3,
          f"residual {res:.2e}, cos {np.array2string(np.sort(c)[::-1], precision=4)}, "
          f"rungs {rungs}")

# (f4) THE PRICING ROWS, and they are the cleanest discriminator in this file: hold the bond
# set fixed and change only what a hop OVER an occupied site costs. The single-excitation
# matrix h is literally the same either way, since with one excitation nothing is ever jumped
# over, so the corner is the same F140 corner and the seed gate reads the same spectrum of h.
# Order-honest pricing (the Slater sign, sector(fermi=True)) leaves H quadratic in the fermions
# and the climber commutes EXACTLY; occupation pricing (sigma+ sigma-, the JW string dropped)
# does not. What the band then does is NOT the same in the two cases, which is why both rows
# are here: the ring keeps its band anyway, the hop-over loses everything but the two ends.
_SPIN_EXPECT = {
    # name -> (predicate on the spin-priced rung list, what the row says)
    "ring": (lambda g, m: all(v == m for v in g),
             "keeps every rung full anyway"),
    "hop-over": (lambda g, m: g[0] == m and g[-1] == m and all(v == 0 for v in g[1:-1]),
                 "loses every rung but the two ends"),
}
for n in (5, 6, 7):
    gnum, m = locus(n), n // 2
    r0, _ = roots(n, gnum)
    for lbl, h in (("ring", h_ring(n)), ("hop-over", h_hopover(n))):
        out = {}
        for fermi in (True, False):
            res = 0.0
            for (p, q) in ((1, 1), (0, 2), (2, 2), (1, 2)):
                if p + 1 > n or q + 1 > n or comb(n, p + 1) * comb(n, q + 1) > 900:
                    continue
                U = _ladder(n, p, q, True)
                up = _dense(n, h, gnum, 768, p + 1, q + 1, fermi)
                res = max(res, np.abs(up @ U - U @ _dense(n, h, gnum, 768, p, q, fermi)).max()
                          / max(np.abs(up).max(), 1.0))
            rungs = [frozen(n, h, gnum, 768, p, p, r0, fermi)[0] for p in range(1, n)
                     if comb(n, p) ** 2 <= EXACT_MAX]
            out[fermi] = (res, rungs)
        (rf, gf), (rs, gs) = out[True], out[False]
        pred, says = _SPIN_EXPECT[lbl]
        check(f"N={n} the {lbl} priced by ORDER: the climber commutes EXACTLY and every "
              f"diagonal rung is floor(N/2) = {m}; priced by OCCUPATION, on the same h, it "
              f"tears and the band {says}",
              rf == 0.0 and all(v == m for v in gf) and rs > 1e-2 and pred(gs, m)
              and len(gf) >= 4,
              f"order-honest |[L, Phi]| = {rf:.1e}, rungs {gf}; occupation-priced "
              f"{rs:.2e}, rungs {gs}")

# (f4a) the off-diagonal column of the two order-honest rows, which the rung census above does
# NOT imply: the (0,2) seed's two-excitation sector is where the two pricings part, so V4's
# occupation-priced reading of it carries nothing over. The seed gate is read here in its own
# right, and it holds: the order-honest ring follows the parity of N, exactly as its
# occupation-priced twin does, while the hop-over's h has no symmetric spectrum at any N here.
for n in (5, 6, 7, 8):
    gnum, m = locus(n), n // 2
    r0, _ = roots(n, gnum)
    want = m if n % 2 == 0 else 0
    got = {lbl: frozen(n, h, gnum, 768, 0, 2, r0, True)[0]
           for lbl, h in (("ring", h_ring(n)), ("hop-over", h_hopover(n)))}
    check(f"N={n} the off-diagonal seed of the order-honest rows: the ring carries "
          f"{want} on the parity of N and the hop-over, whose h is never symmetric, carries 0",
          got["ring"] == want and got["hop-over"] == 0
          and is_bipartite_spectrum(h_ring(n)) == (n % 2 == 0)
          and not is_bipartite_spectrum(h_hopover(n)), f"{got}")

# (f4'') what the missing ladder costs, read where it shows: at N = 4 the occupation-priced
# ring carries FOUR on its middle rung against floor(N/2) = 2, and the occupation-priced
# hop-over carries ONE in the off-diagonal seed where its h-spectrum is not symmetric and the
# seed gate says nothing should carry. Both are exact, both hold at five couplings and on two
# locus profiles, and both are N = 4 only: at N = 5, 6, 7 the same reads are floor(N/2) and 0
# on both profiles. Order-honest pricing has neither exception, which is the point of reading
# them side by side: a rung count is held equal by the ladder, and where the ladder is gone
# nothing holds it. N = 4 is also the one N where the frozen root and its fold partner
# coincide (4*gbar - 2*sigma = -4*gbar exactly at N = 4), which is what leaves room for a count
# above floor(N/2) there at all; the coincidence is checked beside the exceptions, since it
# says the overshoot is not forbidden at N = 4, not that it has to happen.
COUPLINGS = (256, 512, 768, 1536, 3072)
PROFILES = ((29, -13, 7, 19), (41, -23, 11, 5))

# (f4a2) N = 4 read WHOLE, because there it is not a handful of exceptions but a different
# reading: the frozen root and its fold partner are the same number, so every block's kernel is
# fed by both families and the table's N >= 5 cells do not apply. The census below is what the
# note carries in place of per-cell footnotes. Note in particular that the occupation-priced
# hop-over, which at N >= 5 keeps only F140's two corners, has its WHOLE diagonal here and half
# the count off it; and that Heisenberg keeps four corner blocks, its own two and the two fold
# corners it only meets at this N.
_N4_CENSUS = {
    "chain": {(0, 2): 2, (1, 1): 2, (1, 3): 2, (2, 0): 2, (2, 2): 2,
              (2, 4): 2, (3, 1): 2, (3, 3): 2, (4, 2): 2},
    "ring by occupation": {(0, 2): 2, (1, 1): 2, (1, 3): 2, (2, 0): 2, (2, 2): 4,
                           (2, 4): 2, (3, 1): 2, (3, 3): 2, (4, 2): 2},
    "ring by order": {(0, 2): 2, (1, 1): 2, (1, 3): 2, (2, 0): 2, (2, 2): 2,
                      (2, 4): 2, (3, 1): 2, (3, 3): 2, (4, 2): 2},
    "hop-over by occupation": {(0, 2): 1, (1, 1): 2, (1, 3): 2, (2, 0): 1, (2, 2): 2,
                               (2, 4): 1, (3, 1): 2, (3, 3): 2, (4, 2): 1},
    "hop-over by order": {(0, 2): 0, (1, 1): 2, (1, 3): 0, (2, 0): 0, (2, 2): 2,
                          (2, 4): 0, (3, 1): 0, (3, 3): 2, (4, 2): 0},
    "Heisenberg": {(0, 2): 0, (1, 1): 2, (1, 3): 2, (2, 0): 0, (2, 2): 0,
                   (2, 4): 0, (3, 1): 2, (3, 3): 2, (4, 2): 0},
}
for half in PROFILES:
    g4 = locus(4, half=half)
    r4, rfold4 = roots(4, g4)

    def _hk(p, q, _g=g4, _r=r4):
        d = comb(4, p) * comb(4, q)
        re, im = _heis(4, _g, 768, p, q)
        re = re - _r * np.eye(d, dtype=np.int64)
        return d - max(rank_mod_p(re, im, pp) for pp in PRIMES)
    readers = {
        "chain": lambda p, q: frozen(4, h_chain(4), g4, 768, p, q, r4, False)[0],
        "ring by occupation": lambda p, q: frozen(4, h_ring(4), g4, 768, p, q, r4, False)[0],
        "ring by order": lambda p, q: frozen(4, h_ring(4), g4, 768, p, q, r4, True)[0],
        "hop-over by occupation":
            lambda p, q: frozen(4, h_hopover(4), g4, 768, p, q, r4, False)[0],
        "hop-over by order":
            lambda p, q: frozen(4, h_hopover(4), g4, 768, p, q, r4, True)[0],
        "Heisenberg": _hk,
    }
    got = {k: {pq: f(*pq) for pq in sorted(band(4))} for k, f in readers.items()}
    check(f"N=4 offsets {half[:2]} the WHOLE band, model by model, where the frozen root and "
          f"its fold partner coincide and the table's N >= 5 reading does not apply",
          got == _N4_CENSUS and r4 == rfold4,
          "; ".join(f"{k} {sorted(v.items())}" for k, v in got.items())
          if got != _N4_CENSUS else f"six models x {len(band(4))} blocks, roots both {r4}")
for half in PROFILES:
    g4 = locus(4, half=half)
    r4, rfold4 = roots(4, g4)
    reads = [(J, frozen(4, h_ring(4), g4, J, 2, 2, r4, False)[0],
              frozen(4, h_ring(4), g4, J, 2, 2, r4, True)[0],
              frozen(4, h_hopover(4), g4, J, 0, 2, r4, False)[0],
              frozen(4, h_hopover(4), g4, J, 0, 2, r4, True)[0]) for J in COUPLINGS]
    check(f"N=4 offsets {half[:2]} the occupation-priced exceptions, which the order-honest "
          f"pricing does not have, at all {len(COUPLINGS)} couplings: the ring's middle rung "
          f"overshoots floor(N/2) = 2 and the hop-over's seed carries where its h-spectrum is "
          f"not symmetric",
          all(rs == 4 and ro == 2 and hs == 1 and ho == 0 for _, rs, ro, hs, ho in reads)
          and r4 == rfold4,
          f"(J, ring occ, ring ord, hop occ, hop ord) {reads}; the two roots coincide at "
          f"N = 4, both {r4}")
for n in (5, 6, 7):
    rows, ok = [], True
    for half in PROFILES:
        g, (r0, rfold) = locus(n, half=half), roots(n, locus(n, half=half))
        a = frozen(n, h_ring(n), g, 768, 2, 2, r0, False)[0]
        b = frozen(n, h_hopover(n), g, 768, 0, 2, r0, False)[0]
        rows.append(f"offsets {half[:2]}: ring (2,2) {a}, hop-over (0,2) {b}")
        ok = ok and a == n // 2 and b == 0 and r0 != rfold
    check(f"N={n} and both exceptions are N = 4 only, on both profiles: the same two reads are "
          f"floor(N/2) = {n // 2} and 0 here, where the fold root no longer coincides with the "
          f"frozen one", ok, "; ".join(rows))

# (f4') and the transport that goes with it: on the order-honest ring the climber carries the
# corner ONTO the rung above, which is exactly what (f3) measures it failing to do on the
# occupation-priced one. Same bonds, same corner, and the ladder present or absent.
for n in (5, 6):
    gnum, m = locus(n), n // 2
    r0, _ = roots(n, gnum)
    rows, ok = [], True
    for lbl, h in (("ring", h_ring(n)), ("hop-over", h_hopover(n))):
        F1 = frozen_basis(n, h, gnum, 768, 1, 1, r0, True)
        F2 = frozen_basis(n, h, gnum, 768, 2, 2, r0, True)
        c = cosines(F2, _ladder(n, 1, 1, True) @ F1)
        rows.append(f"{lbl} {F1.shape[1]}->{F2.shape[1]}, min cos {c.min():.9f}")
        ok = ok and F1.shape[1] == m and F2.shape[1] == m and c.min() > 1 - 1e-7
    check(f"N={n} order-honest pricing restores the transport on BOTH deformations: Phi maps "
          f"the corner's frozen subspace ONTO the rung above", ok, "; ".join(rows))

# (g) the ladder is wider than the theorem: it knows nothing about the R90 locus.
for n in (5, 6):
    h, gnum = h_chain(n), [37, 91, 5, 143, 61, 17][:n]
    r0, _ = roots(n, gnum)
    U = _ladder(n, 1, 1, True)
    res = np.abs(_dense(n, h, gnum, 768, 2, 2) @ U - U @ _dense(n, h, gnum, 768, 1, 1)).max()
    frz = [frozen(n, h, gnum, 768, p, q, r0)[0] for (p, q) in ((1, 1), (0, 2), (2, 2))]
    check(f"N={n} off the R90 locus the ladder still commutes exactly, and freezes nothing: "
          f"the symmetry is wider than the theorem it serves",
          res == 0.0 and all(v == 0 for v in frz), f"residual {res:.1e}, dims {frz}")

print()
print("V9  the CEILING: what one exact rank buys, and where it is free")

# The floor is an argument and holds at every coupling; the ceiling was a census. The gap is
# shorter than it looks. With the dyadic scaling cleared the block is a Gaussian-integer matrix,
# and reduction mod p is a ring homomorphism into a field, so rank can only DROP:
#
#     nullity over GF(p)  >=  nullity over Q(i).
#
# An exact read of floor(N/2) therefore CERTIFIES the rational nullity is at most floor(N/2),
# which is the direction the ceiling was missing and the one a float nullity can never supply.
# The floor supplies the other direction, so the two meet and one prime is enough: an unlucky
# prime can only make the read too LARGE, and a too-large read certifies nothing instead of
# certifying something false.


def _rank_reference(re, im, p):
    """The object-dtype elimination this file used before the vectorized one landed. Kept for
    exactly one check: the fast routine carries every read above, so it needs a witness."""
    r = sqrt_minus_one(p)
    a = (((re % p) + r * (im % p)) % p).astype(object)
    rows, cols = a.shape
    rank = 0
    for col in range(cols):
        piv = next((row for row in range(rank, rows) if a[row, col] % p), None)
        if piv is None:
            continue
        if piv != rank:
            a[[rank, piv]] = a[[piv, rank]]
        inv = modpow(int(a[rank, col]), p - 2, p)
        a[rank, col:] = (a[rank, col:] * inv) % p
        for row in range(rank + 1, rows):
            f = int(a[row, col]) % p
            if f:
                a[row, col:] = (a[row, col:] - f * a[rank, col:]) % p
        rank += 1
        if rank == rows:
            break
    return rank


def _M(n, h, gnum, jnum, p, q, root_num, fermi=False):
    re, im = scaled_block(n, h, gnum, jnum, p, q, fermi)
    d = re.shape[0]
    return re - root_num * np.eye(d, dtype=np.int64), im


def certify(n, h, gnum, jnum, p, q, root_num, extra=None, fermi=False):
    """The certified nullity: an exact GF(p) read, of the block alone or of the block STACKED
    with a second integer map, which reads dim(ker M cap ker extra) in one rank.

    The astype is not decoration. _ladder builds its matrix without a dtype, so it is float64,
    and stacking it under an int64 block silently upcasts the whole thing; the reduction at
    r*(im mod p) then lands above 2^53 and the residues come out wrong, biased towards full
    rank. That is invisible in a check that asserts a nullity of ZERO, which is exactly what
    the two sl(2) shapes assert, so the cast is what keeps them from passing vacuously. The
    entries of _ladder are 0 and +-1, so the cast itself is exact."""
    re, im = _M(n, h, gnum, jnum, p, q, root_num, fermi)
    if extra is not None:
        re = np.vstack([re, extra]).astype(np.int64)
        im = np.vstack([im, np.zeros_like(extra)]).astype(np.int64)
    return re.shape[1] - max(rank_mod_p(re, im, pp) for pp in PRIMES)


# (a) the witness for the routine that carries every read in this file
for n, p, q in ((5, 2, 2), (4, 1, 3)):
    gnum = locus(n)
    r0, _ = roots(n, gnum)
    re, im = _M(n, h_chain(n), gnum, 768, p, q, r0)
    check(f"N={n} ({p},{q}) the vectorized elimination and the reference one agree exactly",
          rank_mod_p(re, im, PRIMES[0]) == _rank_reference(re, im, PRIMES[0]),
          f"rank {rank_mod_p(re, im, PRIMES[0])} of {re.shape[0]}")

# (b) the certificate itself, every rung of all three ladders. What this settles is the ceiling
# at THIS point of the (J, gbar, delta) space, and, since the entries are polynomial there,
# everywhere off a proper algebraic subset of it. That is JOINT genericity in coupling AND
# profile: it is weaker than the corner's own tightness theorem (PROOF_R90_FROZEN_DIVISOR
# Section 7), which holds at every locus profile because a closed form does the work there.
for n in (5, 6, 7):
    gnum, m = locus(n), n // 2
    r0, _ = roots(n, gnum)
    diag = {p: certify(n, h_chain(n), gnum, 768, p, p, r0) for p in range(1, n)}
    off = {p: certify(n, h_chain(n), gnum, 768, p, p + 2, r0) for p in range(0, n - 1)
           if comb(n, p) * comb(n, p + 2) <= 1300}
    check(f"N={n} EXACT nullity floor(N/2) = {m} on every diagonal rung, so the ceiling there is "
          f"certified and not measured",
          all(v == m for v in diag.values()) and len(diag) == n - 1, f"{sorted(diag.items())}")
    # The rung count is asserted for the same reason it is in V8 (d): the size filter above must
    # not be able to empty this loop and leave all() true on nothing.
    check(f"N={n} the same certificate on the off-diagonal ladder, where it settles the ceiling "
          f"and leaves the floor untouched (the transposed ladder follows by the fold)",
          all(v == m for v in off.values()) and len(off) >= n - 2,
          f"{sorted(off.items())}")

# (c) the two sl(2) shapes of the same fact, exactly rather than by principal cosines, and at
# EVERY step of the ladder rather than the first. ker(M) cap ker(Psi) = 0 says no irreducible
# component has its lowest weight above the seed; Phi injective on a frozen subspace of the same
# dimension as the one above says Phi lands ONTO it.
for n in (5, 6, 7):
    gnum, m = locus(n), n // 2
    r0, _ = roots(n, gnum)
    lw = {p: certify(n, h_chain(n), gnum, 768, p, p, r0, _ladder(n, p, p, False))
          for p in range(2, n)}
    inj = {p: certify(n, h_chain(n), gnum, 768, p, p, r0, _ladder(n, p, p, True))
           for p in range(1, n - 1)}
    check(f"N={n} ker(M) cap ker(Psi) = 0 on every rung above the seed", all(v == 0 for v in lw.values()),
          f"{sorted(lw.items())}")
    check(f"N={n} Phi is injective on every rung's frozen subspace, hence onto the next one, at "
          f"every step", all(v == 0 for v in inj.values()), f"{sorted(inj.items())}")
    # The two checks above assert a nullity of ZERO, and a stacked read that has silently gone
    # wrong returns zero too, so on their own they cannot tell a certificate from a corruption.
    # These two run the SAME path and must come back NONZERO: at the two ends of the ladder the
    # climber and the descender annihilate the frozen subspace whole, which V8 reads as a float
    # norm and this reads exactly. They are the sl(2) statement and the witness at once.
    ends = {"Phi at the top rung": certify(n, h_chain(n), gnum, 768, n - 1, n - 1, r0,
                                           _ladder(n, n - 1, n - 1, True)),
            "Psi at the seed": certify(n, h_chain(n), gnum, 768, 1, 1, r0,
                                       _ladder(n, 1, 1, False))}
    check(f"N={n} both ends annihilated, exactly: the stacked read returns the whole "
          f"floor(N/2) = {m} there, which is also what keeps the two zero-checks honest",
          all(v == m for v in ends.values()), f"{ends}")

# (d) where the ceiling is free, and the condition that makes it so. At room 1 the J = 0 block is
# DIAGONAL, and if no non-mirror pair of offsets cancels, its zeros are exactly the floor(N/2)
# balanced-pair cells; rank can only rise, so floor(N/2) caps the nullity at every N with no
# computation at all. The condition is not decoration: the profile this file ships repeats its
# offsets from N = 10 on, and there the count breaks. Both sides are read.
for n in (5, 6, 7, 8, 9, 10, 12):
    m = n // 2
    for half, lbl in (((29, -13, 7, 19, 41, 53), "distinct offsets"),
                      ((29, -13, 7, 19), "the shipped offsets, which repeat from N = 10")):
        gnum = locus(n, half=half)
        r0, _ = roots(n, gnum)
        re0, im0 = _M(n, h_chain(n), gnum, 0, 0, 2, r0)
        zeros = int(np.sum((np.diag(re0) == 0) & (np.diag(im0) == 0)))
        offd = int(np.abs(re0 - np.diag(np.diag(re0))).max() + np.abs(im0).max())
        want = m if (lbl.startswith("distinct") or n < 10) else None
        if want is not None:
            check(f"N={n} at J=0 the room-1 seed (0,2) is diagonal with exactly floor(N/2) = {m} "
                  f"zeros ({lbl})", zeros == m and offd == 0,
                  f"{zeros} zeros, off-diagonal weight {offd}")
        else:
            check(f"N={n} with repeating offsets a NON-MIRROR pair balances too and the free "
                  f"ceiling loses its count", zeros > m, f"{zeros} zeros against floor(N/2) = {m}")

# and the condition is on the offsets' ABSOLUTE VALUES, not on the offsets: two offsets that are
# negatives of one another balance a non-mirror pair just as a repeated one does, so three
# distinct numbers are not enough. This is the case the N = 10 row above cannot see.
for half, lbl, want_free in (((5, -5, 7), "(5, -5, 7): distinct, but two of them cancel", False),
                             ((29, -13, 7), "(29, -13, 7): distinct in absolute value", True)):
    n, m = 6, 3
    gnum = locus(n, half=half)
    r0, _ = roots(n, gnum)
    re0, im0 = _M(n, h_chain(n), gnum, 0, 0, 2, r0)
    zeros = int(np.sum((np.diag(re0) == 0) & (np.diag(im0) == 0)))
    check(f"N={n} the free ceiling's condition is on |offsets|, {lbl}",
          (zeros == m) == want_free, f"{zeros} zeros against floor(N/2) = {m}")

# (e) the uniform point delta = 0 lies ON the locus, and the count survives there. Since a
# special point of a polynomial family can only carry MORE kernel, an all-N argument at the
# uniform point alone would cap every rung generically. That move is the proof document's own at
# the corner (Section 7, where it is sharper: a closed form, every J != 0); what is open is the
# rungs above it. Note what does NOT survive the trip: at delta = 0 EVERY pair balances, so the
# balanced-pair count that indexes the modes on a generic profile is not what counts here.
for n in (5, 6, 7):
    gu, m = [92] * n, n // 2
    r0, _ = roots(n, gu)
    plain = {p: certify(n, h_chain(n), gu, 768, p, p, r0) for p in range(1, n)
             if comb(n, p) ** 2 <= 1300}
    dg = [3 if a in (0, n - 1) else 0 for a in range(n)]
    diag_h = {p: certify(n, h_chain(n, dg), gu, 768, p, p, r0) for p in range(1, n)
              if comb(n, p) ** 2 <= 1300}
    check(f"N={n} the uniform profile carries floor(N/2) = {m} on every diagonal rung, so the "
          f"all-N question may be asked there", all(v == m for v in plain.values()),
          f"{sorted(plain.items())}")
    check(f"N={n} and it is not h's chiral pairs doing the counting: an R-invariant diagonal "
          f"leaves none and the count stands", all(v == m for v in diag_h.values()),
          f"{sorted(diag_h.items())}")

# (f) what the uniform point IS, entry for entry: |A^B| = p + q - 2*D with D the cell's double
# occupancy, so the recentred rate diagonal is exactly U*(D - (p-1)) with U = 4*gamma. Hence
# M = i*J*K + U*D - U*(p-1): a frozen mode is an eigenvector of a Hubbard operator (hopping
# inverted on the bra copy, interaction U) at the eigenvalue U*(p-1), and Phi is the eta-pairing
# operator, which raises THAT energy by exactly U while leaving the Liouvillian's eigenvalue
# alone. The two statements are about two different operators, E = M + U*(p-1) and M.
for n in (4, 5, 6):
    gu = [92] * n
    r0, _ = roots(n, gu)
    worst = 0
    for p in range(1, n):
        if comb(n, p) ** 2 > 1300:
            continue
        re, im = _M(n, h_chain(n), gu, 768, p, p, r0)
        st = list(combinations(range(n), p))
        f, u = 4 * n, 4 * 92
        want = np.array([f * u * (len(set(a) & set(b)) - (p - 1)) for a in st for b in st])
        worst = max(worst, int(np.abs(np.diag(re) - want).max()),
                    int(np.abs(re - np.diag(np.diag(re))).max()))
    check(f"N={n} at the uniform point the recentred rate diagonal is exactly U*(D - (p-1)), "
          f"U = 4*gamma, and nothing else", worst == 0, f"max deviation {worst}")

# (g) the reach, opt-in because it costs minutes rather than seconds: the same certificate at the
# N the default run cannot afford. "python xy_frozen_band.py --deep" adds N = 8 (every diagonal
# rung, the largest block 4900 x 4900, a few minutes); "--deep 9" adds N = 9, whose two middle
# rungs are 15876 x 15876 and take about an hour and forty minutes each, the whole run about
# three and a half hours. At N = 9 one of the two
# settles both, since rungs 4 and 5 sit at weights -1 and +1 of the same sl(2) module and weight
# multiplicities are symmetric.
if "--deep" in sys.argv:
    i = sys.argv.index("--deep")
    rest = sys.argv[i + 1:]
    if any(not a.isdigit() for a in rest):
        raise SystemExit(f"--deep takes N values and nothing else, got {rest}")
    deep_n = [int(a) for a in rest] or [8]
    if any(n < 4 or n > 9 for n in deep_n):
        raise SystemExit("--deep is meant for N = 8 and 9; N = 10 needs 290 GB for one block")
    print(f"\nV9-deep  the certificate at N = {deep_n}")
    for n in deep_n:
        gnum, m = locus(n), n // 2
        r0, _ = roots(n, gnum)
        got = {}
        for p in range(1, n):
            t = time.time()
            got[p] = certify(n, h_chain(n), gnum, 768, p, p, r0)
            print(f"       N={n} ({p},{p}) d={comb(n, p) ** 2} nullity {got[p]} "
                  f"[{time.time() - t:.0f}s]", flush=True)
        check(f"N={n} EXACT nullity floor(N/2) = {m} on every diagonal rung",
              all(v == m for v in got.values()) and len(got) == n - 1,
              f"{sorted(got.items())}")
        # The rival reading of the N = 4 ring exception is "even N overshoots", and the middle
        # rung of an even N is where it would show. At N = 8 it does not: both pricings read
        # floor(N/2) there, so the overshoot stays with the one N whose two roots coincide.
        if n % 2 == 0:
            t = time.time()
            ring = {lbl: certify(n, h_ring(n), gnum, 768, n // 2, n // 2, r0, fermi=fermi)
                    for lbl, fermi in (("occupation", False), ("order", True))}
            print(f"       N={n} ring ({n//2},{n//2}) d={comb(n, n // 2) ** 2} {ring} "
                  f"[{time.time() - t:.0f}s]", flush=True)
            check(f"N={n} the ring's middle rung reads floor(N/2) = {m} under BOTH pricings, so "
                  f"the N = 4 overshoot is not an even-N effect",
                  all(v == m for v in ring.values()), f"{ring}")

print()
if FAILURES:
    print(f"XY frozen band: {len(FAILURES)} FAILURE(S)")
    for f in FAILURES:
        print(f"  - {f}")
    raise SystemExit(1)
print("XY frozen band: ALL GREEN")
