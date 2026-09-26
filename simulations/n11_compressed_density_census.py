"""Measured census around the N=11 compressed-density proof (docs/proofs/PROOF_N11_COMPRESSED_DENSITY.md).

The exact N=11 statements live in simulations/n11_compressed_density_gate.py (SymPy) and in the C#
witness CompressedDensityN11Witness (exact Q(sqrt2, sqrt3), all 55 rooms). This producer holds the
float readings the proof quotes beyond that, each with its error model:

  (A) per-room lower endpoint: mult(-4 gbar) against the number of non-fixed chiral-transpose orbits
      (a,b) -> (M-b, M-a), uniform open XY (1,1) block, N = 2..14, three random locus profiles each;
  (B) parity-mixed (1,1) rooms and their firing site contrast, N = 2..30, frequencies grouped at 40 digits;
  (C) the zero-frequency room for random reflection-symmetric simple h, real and complex, N = 2..11:
      D_Omega0 = -4 gbar (I - W^T W) and the multiplicities of -4 gbar and 0;
  (D) the finite-J (1,1) Liouvillian block at N = 11: the small-J escape below -4 gbar and the kernel of
      L + 4 gbar (F140);
  (E) the simple-spectrum hypothesis has teeth: an N = 6 reflection-symmetric h with degenerate levels.

Error model for every eigenvalue count: eigvalsh is backward stable, so a computed eigenvalue sits within
c * eps * ||D|| of the true one with c of order the dimension. A multiplicity is counted inside
TOL = 1e3 * eps * ||D|| * dim, and the next eigenvalue's distance is reported in units of
eps * ||D|| * dim ("rounding units"); a miscount would show as a small ratio. Singular values of L + 4 gbar
are printed against eps * ||L|| * dim(L) as readings only. In (C) the eigenvectors of h carry the
conditioning 1/(min level gap), so the residual of the identity is printed both in rounding units and
multiplied by that gap; both are readings, and no law is claimed for either.

The kernel of L + 4 gbar in (D) is not read off singular values: F140's J^(2d) ladder puts a genuine sixth
singular value near 1e-10 at J = 0.3, inside any rounding tolerance. It is certified exactly instead. J is
rational and the rates are integers, so L + 4 gbar has entries in Q(i); reducing modulo a prime p = 1 mod 4,
with i sent to a square root of -1, is a ring map, so rank over F_p <= rank over Q(i) and dim ker over F_p is
an UPPER bound on the true kernel. F140 gives the LOWER bound floor(N/2) = 5. Where the two meet the kernel is
exactly 5; the certificate is run at two primes.

Writes simulations/results/n11_compressed_density_census.txt.
Run: python simulations/n11_compressed_density_census.py
"""
from pathlib import Path

from fractions import Fraction

import mpmath as mp
import numpy as np

EPS = np.finfo(float).eps
OUT = Path(__file__).resolve().parent / "results" / "n11_compressed_density_census.txt"
LINES = []


def say(text=""):
    print(text)
    LINES.append(text)


def modes_xy(n):
    m = n + 1
    z = np.arange(n)
    psi = {k: np.sqrt(2 / m) * np.sin(np.pi * k * (z + 1) / m) for k in range(1, n + 1)}
    energy = {k: 4 * np.cos(np.pi * k / m) for k in range(1, n + 1)}
    return psi, energy


def rooms_exactish(n):
    """(1,1) frequency rooms of the uniform XY chain, grouped at 40 digits; returns rooms and margins."""
    mp.mp.dps = 40
    m = n + 1
    e = {k: 4 * mp.cos(mp.pi * k / m) for k in range(1, n + 1)}
    fr = sorted(((e[a] - e[b], (a, b)) for a in range(1, n + 1) for b in range(1, n + 1)),
                key=lambda t: t[0])
    rooms, cur = [], [fr[0][1]]
    same_max, distinct_min = mp.mpf(0), mp.mpf(100)
    for i in range(1, len(fr)):
        d = fr[i][0] - fr[i - 1][0]
        if d < mp.mpf(10) ** -30:
            cur.append(fr[i][1])
            same_max = max(same_max, d)
        else:
            rooms.append(cur)
            cur = [fr[i][1]]
            distinct_min = min(distinct_min, d)
    rooms.append(cur)
    return rooms, float(same_max), float(distinct_min)


def cell_rates(g):
    n = len(g)
    return np.array([[-2 * (g[a] + g[b]) if a != b else 0.0 for b in range(n)] for a in range(n)])


def compress(vectors, g):
    """<v_i, D v_j> over the physical (1,1) cells, D diagonal with the cell rates."""
    r = cell_rates(g).ravel()
    v = np.array([x.ravel() for x in vectors])
    d = (v.conj() * r) @ v.T
    return (d + d.conj().T) / 2


def count_at(eigs, target, scale):
    dist = np.abs(eigs - target)
    tol = 1e3 * scale
    mult = int(np.sum(dist < tol))
    rest = np.sort(dist[dist >= tol])
    return mult, (rest[0] / scale if len(rest) else np.inf)


def locus_profile(rng, n, gbar=1.0, signed=False):
    g = np.zeros(n)
    for l in range(n // 2):
        x = rng.uniform(-3, 5) if signed else rng.uniform(0, 2 * gbar)
        g[l], g[n - 1 - l] = x, 2 * gbar - x
    if n % 2:
        g[n // 2] = gbar
    return g


def section_a():
    say("(A) per-room lower endpoint, uniform open XY (1,1): mult(-4 gbar) vs non-fixed chiral-transpose orbits")
    rng = np.random.default_rng(20260926)
    for n in range(2, 15):
        m = n + 1
        psi, _ = modes_xy(n)
        rooms, _, margin = rooms_exactish(n)
        mismatches, zero_rooms, low_rooms, worst = 0, 0, 0, np.inf
        for _ in range(3):
            g = locus_profile(rng, n)
            for room in rooms:
                vecs = [np.outer(psi[a], psi[b]) for a, b in room]
                d = compress(vecs, g)
                eigs = np.linalg.eigvalsh(d)
                scale = EPS * max(1.0, np.abs(d).max()) * len(room)
                mult, gap = count_at(eigs, -4.0, scale)
                zero_mult, _ = count_at(eigs, 0.0, scale)
                orbits = len({frozenset({(a, b), (m - b, m - a)}) for a, b in room if (a, b) != (m - b, m - a)})
                mismatches += mult != orbits
                zero_rooms += zero_mult > 0
                low_rooms += mult > 0
                worst = min(worst, gap)
        say(f"  N={n:2d}: rooms={len(rooms):3d} x 3 profiles; mult(-4gbar) != #orbits in {mismatches} rooms; "
            f"rooms reaching -4gbar per profile={low_rooms // 3}; rooms reaching 0 per profile={zero_rooms // 3}; "
            f"min next-eigenvalue gap = {worst:.2e} rounding units; min distinct frequency spacing = {margin:.2e}")


def section_b():
    say("(B) parity-mixed (1,1) rooms on the uniform open XY chain, N = 2..30 (40-digit grouping)")
    hits = []
    for n in range(2, 31):
        psi, _ = modes_xy(n)
        rooms, same_max, margin = rooms_exactish(n)
        mixed = [r for r in rooms if len({(a + b) % 2 for a, b in r}) == 2]
        firing = 0
        for room in mixed:
            v = np.array([psi[a] * psi[b] for a, b in room]).T
            c0 = -2 * (np.outer(v[0], v[0]) - np.outer(v[n - 1], v[n - 1]))
            firing += np.abs(c0).max() > 1e-6
        if mixed:
            hits.append(n)
        say(f"  N={n:2d} M={n + 1:2d}: mixed rooms={len(mixed):3d}, with C_0 != 0: {firing:3d}; "
            f"min distinct frequency spacing {margin:.2e}")
    say(f"  N with parity-mixed (1,1) rooms: {hits}")


def section_c():
    say("(C) zero-frequency room for random reflection-symmetric simple h: D = -4 gbar (I - W^T W)")
    rng = np.random.default_rng(11)
    for n in range(2, 12):
        refl = np.eye(n)[::-1]
        for cplx in (False, True):
            a = rng.normal(size=(n, n)) + (1j * rng.normal(size=(n, n)) if cplx else 0)
            h = a + a.conj().T
            h = (h + refl @ h @ refl) / 2
            e, u = np.linalg.eigh(h)
            g = locus_profile(rng, n, gbar=1.0, signed=not cplx)
            vecs = [np.outer(u[:, k], u[:, k].conj()) for k in range(n)]
            d = compress(vecs, g)
            w = np.abs(u) ** 2
            scale = EPS * max(1.0, np.abs(d).max()) * n
            residual = np.abs(d - (-4.0) * (np.eye(n) - w.T @ w)).max() / scale
            eigs = np.linalg.eigvalsh(d)
            low, low_gap = count_at(eigs, -4.0, scale)
            top, top_gap = count_at(eigs, 0.0, scale)
            sv = np.linalg.svd(w, compute_uv=False)
            rank = int(np.sum(sv > 1e-10 * sv[0]))
            gap_h = np.diff(e).min()
            say(f"  N={n:2d} {'complex' if cplx else 'real   '} ({'nonneg' if cplx else 'signed'} profile): "
                f"residual {residual:.2f} rounding units, times min level gap {residual * gap_h:.2f}; "
                f"min level gap {gap_h:.3f}; "
                f"mult(-4gbar)={low} (floor(N/2)={n // 2}), next gap {low_gap:.2e}; mult(0)={top}, next gap "
                f"{top_gap:.2e}; rank W={rank} (ceil(N/2)={(n + 1) // 2}), kept sv >= {sv[rank - 1]:.2e}, "
                f"dropped <= {sv[rank] if rank < n else 0.0:.2e}")


def block11(n, j, g):
    hop = np.zeros((n, n))
    for i in range(n - 1):
        hop[i, i + 1] = hop[i + 1, i] = 2 * j
    ident = np.eye(n)
    lh = -1j * (np.kron(hop, ident) - np.kron(ident, hop.T))
    return lh + np.diag(cell_rates(g).ravel())


PRIMES = (1000000009, 998244353)   # both 1 mod 4


def sqrt_minus_one(prime):
    for base in range(2, prime):
        root = pow(base, (prime - 1) // 4, prime)
        if root * root % prime == prime - 1:
            return root
    raise ValueError(prime)


def block11_mod_p(n, j, rates, prime):
    """L + 4 on the (1,1) cells over F_p: -i J (h x I - I x h^T) + diag(cell rates) + 4, h = 2 x hopping."""
    i_p = sqrt_minus_one(prime)
    j_p = j.numerator * pow(j.denominator, -1, prime) % prime
    hop = [[2 if abs(a - b) == 1 else 0 for b in range(n)] for a in range(n)]
    size = n * n
    mat = [[0] * size for _ in range(size)]
    for a in range(n):
        for b in range(n):
            row = a * n + b
            rate = -2 * (rates[a] + rates[b]) if a != b else 0
            mat[row][row] = (rate + 4) % prime
            for c in range(n):
                if hop[a][c]:
                    mat[row][c * n + b] = (mat[row][c * n + b] - i_p * j_p * hop[a][c]) % prime
                if hop[c][b]:
                    mat[row][a * n + c] = (mat[row][a * n + c] + i_p * j_p * hop[c][b]) % prime
    return mat


def rank_mod_p(mat, prime):
    m = [row[:] for row in mat]
    rank, cols = 0, len(m[0])
    for col in range(cols):
        pivot = next((r for r in range(rank, len(m)) if m[r][col]), None)
        if pivot is None:
            continue
        m[rank], m[pivot] = m[pivot], m[rank]
        inv = pow(m[rank][col], -1, prime)
        m[rank] = [x * inv % prime for x in m[rank]]
        for r in range(len(m)):
            if r != rank and m[r][col]:
                f = m[r][col]
                m[r] = [(x - f * y) % prime for x, y in zip(m[r], m[rank])]
        rank += 1
    return rank


def section_d():
    say("(D) finite-J (1,1) Liouvillian block, uniform open XY, N = 11")
    n = 11
    g = np.array([2, 2] + [1] * 7 + [0, 0], float)
    for j in (0.0, 0.05, 0.2, 1.0):
        say(f"  rates (2,2,1,...,1,0,0), J={j}: min Re lambda = {np.linalg.eigvals(block11(n, j, g)).real.min():.4f}")
    g = np.array([2] + [1] * 9 + [0], float)
    rates = [2] + [1] * 9 + [0]
    for j in (Fraction(3, 10), Fraction(1, 2), Fraction(1), Fraction(37, 10), Fraction(25)):
        kernels = [n * n - rank_mod_p(block11_mod_p(n, j, rates, prime), prime) for prime in PRIMES]
        lmat = block11(n, float(j), g) + 4 * np.eye(n * n)
        sv = np.sort(np.linalg.svd(lmat, compute_uv=False))
        scale = EPS * np.linalg.norm(lmat, 2) * n * n
        say(f"  rates (2,1,...,1,0), J={j}: dim ker(L+4) over F_p = {kernels} at p = {PRIMES}, F140 lower bound "
            f"{n // 2}, so dim ker = {n // 2 if max(kernels) == n // 2 else 'not certified'}; six smallest sv = "
            f"{', '.join(f'{x:.2e}' for x in sv[:6])}; next = {sv[6]:.2e}; rounding unit = {scale:.2e}; "
            f"J^20 = {float(j) ** 20:.2e}")


def section_e():
    say("(E) simple-spectrum tooth: N = 6, dimers on (0,1) and (4,5) at coupling 1, middle dimer (2,3) at 0.3")
    n = 6
    h = np.zeros((n, n))
    for a, b in ((0, 1), (4, 5)):
        h[a, b] = h[b, a] = 1.0
    h[2, 3] = h[3, 2] = 0.3
    g = np.array([2, 2, 1, 1, 0, 0], float)
    ad = np.kron(h, np.eye(n)) - np.kron(np.eye(n), h.T)
    w, u = np.linalg.eigh(ad)
    say(f"  one-body levels {np.round(np.linalg.eigvalsh(h), 6).tolist()}; gbar = {g.mean():.1f}")
    dmat = np.diag(cell_rates(g).ravel())
    for value in np.unique(np.round(w, 9)):
        p = u[:, np.round(w, 9) == value]
        lo = np.linalg.eigvalsh(p.T @ dmat @ p).min()
        if lo < -4 - 1e-9:
            say(f"  complete ad_H room at omega={value:+.1f}, dim {p.shape[1]}: min eigenvalue {lo:.6f} < -4gbar = -4")


def main():
    section_a()
    section_b()
    section_c()
    section_d()
    section_e()
    OUT.write_bytes(("\n".join(LINES) + "\n").encode("utf-8"))
    print(f"written {OUT}")


if __name__ == "__main__":
    main()
