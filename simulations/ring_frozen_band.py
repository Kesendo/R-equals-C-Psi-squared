"""Gate for docs/proofs/PROOF_RING_FROZEN_BAND.md: the XY ring's frozen band at odd N, and at even N
above the corner.

Self-contained: builds every object from scratch and imports nothing from the framework or
from the other frozen-band gates, so agreement with them is evidence and not a shared bug.

Run:
    python simulations/ring_frozen_band.py
    python simulations/ring_frozen_band.py --deep    # + the exact rungs at N = 7, several minutes

What it checks:

  R1  the reduction: the spin ring's (p,p) sector Hamiltonian, priced by occupation, IS the
      fermion ring's with wrap sign (-1)^(p-1), priced by order (integer matrices, entry for entry).
  R2  the resonance lemma: two allowed momenta of the periodic or the pi-flux ring sum to pi
      (mod 2 pi) at no odd N (exact integer enumeration, N = 3..41), and at every even N some do.
  R3  the identity N*Dhat + Dmom = l^2 on LW_l cap V0, both rings, every rung, odd N = 5, 7, an
      operator statement, matrix against matrix; the basis of LW_l cap V0 comes from an SVD, so
      the residual is a rounding residual rather than an exact zero, gated against its error
      model, eps times the dimension, below 1e3 of it; the even-N control N = 4, 6 misses by
      more than 0.1.
  R4  the floor it yields: min K on LW_l cap V0 is l(N-l)/N, attained, and K = 1 does not occur
      at l >= 2 for odd N >= 5.
  R5  the corner at large J: span{P_k - P_-k} lies in the kernel of the compressed rate
      operator for every locus profile, and is the whole kernel on random locus profiles at odd
      N = 5..15, read inside a moat (no singular value between 1e-11 and 1e-6 relative: the
      zeros sit at rounding, the lifted values are second order in the offsets and sit near
      1e-5 and above at these profiles); at even
      N the kernel is larger than floor(N/2) (order one is loose there).
  R6  exact counts, rank modulo two primes (an upper bound that meets F140's floor): the ring
      corner carries floor(N/2) on a locus profile with distinct offsets, both rings, N = 3..21;
      every diagonal rung carries floor(N/2) on the locus and N - 1 at the uniform point, all of
      it lowest weight at rung 1, at N = 3 and 5 (N = 7 under --deep).
  R7  the even-N reads the proof leaves open: N = 4 locus [2, 4, 2] with the extra two lowest
      weight at rung 2 of the pi-flux ring; the uniform point at N = 6 carries lowest weight at
      rung 2 of the periodic ring, so it is no anchor there. The mod-p nullities are upper
      bounds; the lowest-weight reads are crossed by a float SVD inside a moat, an independent
      second route (not a bound from below).
  R8  the exceptional profiles: an offset on the single pair (l, N-1-l) gives the corner
      floor(N/2) + (gcd(2l+1, N) - 1)/2, both rings, two couplings, N = 9, 15, 21, every l
      (nullities modulo primes, upper bounds; above the floor they are reads); and, exact in
      integers, F157's node law at the seat l, (2l+1)m = N (mod 2N), picks exactly one member of
      each missed momentum pair {m, N-m}, m(2l+1) = 0 (mod N), odd N to 199; the same check run with
      the XY book's nodes, (l+1)m = 0 (mod N+1), must fail. Which member is the node is read in
      floats, to N = 31 inside a moat: the zeros of the offset's transform at 2k = 2 pi m / N are the
      missed set, the nodes of the Heisenberg chain's single-excitation eigenvectors at the seat l are
      the node law's set (both read through the same functions the integer check uses), on a
      range holding 14 seats with gcd(2l+1, N) > 1, the seats with pairs for the read to tell apart.
  R9  even N: k -> pi - k maps the allowed momenta onto themselves; the identity
      N*Dhat = l(l+1) - |S v|^2 - (1/2) sum_a |g_a + g_abar|^2 - sum_{a != abar} |Y_{a abar}|^2 on
      LW cap V0 (S the staggered spin ladder, g_a = Y_aa, Y_ab : v -> c_a v c_b^dagger), N = 4, 6, 8,
      both rings, with the unstaggered ladder as the control on every rung l >= 2; its floor l(N-1-l)/N attained; the
      corner at the uniform point exactly 3N/2 - 2 (periodic) and 3N/2 - 1 (pi flux), N = 4..20,
      and the same number as the nullity of the commutant with no diagonal; the spin ring's
      uniform rungs 1 and 2 at N = 8, carrying 10 and 11; and, read in floats inside a moat, the large-J kernel's
      excess at even N = 4..14 as one mode per orbit of k -> -k, k -> k + pi; and the corner's split
      under U v(a,b) = (-1)^(a+b) v(b,a) (PROOF_R90_FROZEN_DIVISOR 9.1), exact at even N = 4..14 on
      integer locus profiles, with U = +1, the populations, empty and U = -1 carrying N/2; U = +1
      reaching the root at the uniform point, N = 4..12; and, exact in sympy, the determinant of the
      populations' sector on one N = 6 periodic locus profile, c*J^10*(13J^2 - 35816)^2.
"""
import sys, random
from itertools import combinations
import numpy as np

DEEP = "--deep" in sys.argv
FAILS = []


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(label)


# ---------------- models ----------------

def h_ring(n, wrap=1):
    h = np.zeros((n, n), dtype=np.int64)
    for a in range(n - 1):
        h[a, a + 1] = h[a + 1, a] = 1
    h[0, n - 1] = h[n - 1, 0] = wrap
    return h


def sector(n, k, h, fermi):
    """h second-quantized on the k-particle sector; fermi prices a hop by order (Slater sign),
    otherwise by occupation (sigma+ sigma-)."""
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
    return st, H


def jw(cfg, l):
    return -1 if sum(1 for x in cfg if x < l) % 2 else 1


def psi(n, p):
    """Psi = Phi^dagger : (p,p) -> (p-1,p-1), rho -> sum_l c_l rho c_l^dagger."""
    st = list(combinations(range(n), p))
    lo = list(combinations(range(n), p - 1))
    ilo = {s: i for i, s in enumerate(lo)}
    M = np.zeros((len(lo) ** 2, len(st) ** 2), dtype=np.int64)
    for i, A in enumerate(st):
        for j, B in enumerate(st):
            for l in set(A) & set(B):
                A2 = tuple(x for x in A if x != l)
                B2 = tuple(x for x in B if x != l)
                M[ilo[A2] * len(lo) + ilo[B2], i * len(st) + j] += jw(A2, l) * jw(B2, l)
    return M


def null(M, tol=1e-9):
    u, s, vh = np.linalg.svd(M)
    sf = np.zeros(M.shape[1]); sf[:len(s)] = s
    k = int(np.sum(sf < tol * max(s[0] if len(s) else 1.0, 1.0)))
    return vh[M.shape[1] - k:].conj().T


def moat_nullity(M, lo=1e-10, hi=1e-4):
    """float nullity read only inside a moat: None if a singular value sits between lo and hi."""
    s = np.linalg.svd(M, compute_uv=False)
    sf = np.zeros(M.shape[1]); sf[:len(s)] = s
    sf = sf / max(sf.max(), 1.0)
    if np.any((sf > lo) & (sf < hi)):
        return None
    return int(np.sum(sf <= lo))


def momenta(n, wrap):
    return [2 * np.pi * m / n for m in range(n)] if wrap == 1 else \
           [(2 * m + 1) * np.pi / n for m in range(n)]


def dgamma(n, l, P):
    st = list(combinations(range(n), l))
    idx = {s: i for i, s in enumerate(st)}
    M = np.zeros((len(st), len(st)), complex)
    for s, cfg in enumerate(st):
        occ = set(cfg)
        for site in cfg:
            M[s, s] += P[site, site]
            for tgt in range(n):
                if tgt == site or tgt in occ or P[tgt, site] == 0:
                    continue
                rest = occ - {site}
                sg = (-1) ** (sum(1 for x in occ if x < site) + sum(1 for x in rest if x < tgt))
                M[idx[tuple(sorted(rest | {tgt}))], s] += sg * P[tgt, site]
    return M


def lw_v0(n, l, wrap):
    st, H = sector(n, l, h_ring(n, wrap), True)
    d = len(st)
    A = (np.kron(H, np.eye(d)) - np.kron(np.eye(d), H.T)).astype(float)
    V0 = null(A)
    Z = null(psi(n, l).astype(float) @ V0)
    Q, _ = np.linalg.qr(V0 @ Z)
    return st, Q


# ---------------- exact ranks ----------------

PRIMES = [2147483629, 2147483549]


def sqrtm1(p):
    for g in range(2, 500):
        x = pow(g, (p - 1) // 4, p)
        if x * x % p == p - 1:
            return x
    raise ValueError(p)


def rank_mod(M, p):
    A = M.copy() % p
    r, (rows, cols) = 0, A.shape
    for c in range(cols):
        piv = np.nonzero(A[r:, c])[0]
        if len(piv) == 0:
            continue
        k = r + piv[0]
        if k != r:
            A[[r, k]] = A[[k, r]]
        A[r] = (A[r] * pow(int(A[r, c]), p - 2, p)) % p
        nz = np.nonzero(A[:, c])[0]
        nz = nz[nz != r]
        if len(nz):
            A[nz] = (A[nz] - (A[nz, c][:, None] * A[r][None, :]) % p) % p
        r += 1
        if r == rows:
            break
    return r


def nullity(re, im, extra=None):
    best = None
    for p in PRIMES:
        s = sqrtm1(p)
        M = (re % p + (s * (im % p)) % p) % p
        if extra is not None:
            M = np.vstack([M, extra % p])
        k = M.shape[1] - rank_mod(M.astype(np.int64), p)
        best = k if best is None else min(best, k)
    return best


def locus(n, rnd, base=500):
    g = [base] * n
    for i, d in enumerate(rnd.sample(range(5, 400), n // 2)):
        s = rnd.choice((1, -1))
        g[i], g[n - 1 - i] = base + s * d, base - s * d
    return g


def scaled(n, p, wrap, fermi, jnum, g):
    """the block (p,p) of L + 4*gbar, scaled by the rate grid, as (re, im) integers."""
    st, H = sector(n, p, h_ring(n, wrap), fermi)
    d = len(st)
    I = np.eye(d, dtype=np.int64)
    im = -jnum * (np.kron(H, I) - np.kron(I, H.T))
    gb4 = 4 * sum(g) // n
    re = np.diag(np.array([[-2 * sum(g[l] for l in set(A) ^ set(B)) + gb4 for B in st]
                           for A in st]).reshape(-1))
    return re.astype(np.int64), im.astype(np.int64)


def rung_counts(n, wrap, fermi, jnum, g):
    tot, lw = [], []
    for p in range(1, n):
        re, im = scaled(n, p, wrap, fermi, jnum, g)
        tot.append(nullity(re, im))
        lw.append(nullity(re, im, psi(n, p)))
    return tot, lw


# ---------------- checks ----------------

print("R1  the spin ring's diagonal sectors are the fermion ring's with wrap sign (-1)^(p-1)")
ok = True
for n in range(3, 9):
    for p in range(1, n):
        _, Hs = sector(n, p, h_ring(n, 1), False)
        _, Hf = sector(n, p, h_ring(n, (-1) ** (p - 1)), True)
        ok &= bool(np.array_equal(Hs, Hf))
check("N = 3..8, every p: integer sector matrices equal entry for entry", ok)
_, a = sector(5, 2, h_ring(5, 1), False)
_, b = sector(5, 2, h_ring(5, 1), True)
check("control: at p = 2 the periodic fermion ring differs from the spin ring", not np.array_equal(a, b))

print("R2  the resonance lemma: no two allowed momenta sum to pi at odd N")
odd_ok, even_ok = True, True
for n in range(3, 42):
    # momenta in units of pi/n: periodic 2m, pi-flux 2m+1; sum = pi (mod 2 pi) <=> sum = n (mod 2n)
    for off in (0, 1):
        ks = [2 * m + off for m in range(n)]
        hit = any((x + y - n) % (2 * n) == 0 for x in ks for y in ks)
        if n % 2:
            odd_ok &= not hit
        else:
            even_ok &= hit
check("odd N = 3..41, periodic and pi-flux: never", odd_ok)
check("even N = 4..40, periodic and pi-flux: always some pair", even_ok)

print("R3  the identity N*Dhat + Dmom = l^2 on LW_l cap V0")
worst, ctrl = 0.0, np.inf
for n in (5, 7, 4, 6):
    for wrap in (1, -1):
        V = [np.exp(1j * k * np.arange(n)) / np.sqrt(n) for k in momenta(n, wrap)]
        for l in range(1, n // 2 + 1):
            st, Q = lw_v0(n, l, wrap)
            Dh = np.array([[len(set(x) & set(y)) for y in st] for x in st], float).reshape(-1)
            ND = n * Q.conj().T @ (Dh[:, None] * Q)
            Dm = sum(np.kron(Nk, Nk.T) for Nk in (dgamma(n, l, np.outer(v, v.conj())) for v in V))
            res = float(np.max(np.abs(ND + Q.conj().T @ Dm @ Q - l * l * np.eye(Q.shape[1]))))
            if n % 2:
                worst = max(worst, res / (np.finfo(float).eps * Q.shape[0]))
            else:
                ctrl = min(ctrl, res)
check("odd N = 5, 7, both rings, every rung: residual / (eps * dim) below 1e3",
      worst < 1e3, f"worst {worst:.1f}")
check("control: even N = 4, 6 miss by more than 0.1", ctrl > 0.1, f"smallest miss {ctrl:.2f}")

print("R4  the floor min K = l(N-l)/N, attained, and no K = 1 above the seed at odd N >= 5")
ok_min, ok_one = True, True
for n in (5, 7):
    for wrap in (1, -1):
        for l in range(1, n // 2 + 1):
            st, Q = lw_v0(n, l, wrap)
            K = np.array([[len(set(x) - set(y)) for y in st] for x in st], float).reshape(-1)
            ev = np.linalg.eigvalsh(Q.conj().T @ (K[:, None] * Q))
            ok_min &= abs(ev[0] - l * (n - l) / n) < 1e-10
            if l >= 2:
                ok_one &= not np.any(np.abs(ev - 1) < 1e-8)
check("N = 5, 7, both rings: min K equals l(N-l)/N", ok_min)
check("N = 5, 7, both rings, l >= 2: K = 1 absent", ok_one)

print("R5  the corner at large J: span{P_k - P_-k} is the kernel of the compressed rate")
rng = np.random.default_rng(11)
ok_in, ok_dim, ok_even = True, True, True
for n in (5, 7, 9, 11, 13, 15, 4, 6, 8):
    for wrap in (1, -1):
        h = h_ring(n, wrap).astype(float)
        I = np.eye(n)
        V0 = null(np.kron(h, I) - np.kron(I, h.T))
        for trial in range(2):
            d = rng.uniform(-0.3, 0.3, n // 2)
            g = np.ones(n)
            for i in range(n // 2):
                g[i], g[n - 1 - i] = 1 + d[i], 1 - d[i]
            R = np.array([[(2 * (g[a] + g[b]) if a != b else 0) - 4 for b in range(n)]
                          for a in range(n)]).reshape(-1)
            B1 = V0.conj().T @ (R[:, None] * V0)
            ker = null(B1)
            kmoat = moat_nullity(B1, lo=1e-11, hi=1e-6)
            if n % 2:
                ks = momenta(n, wrap)
                Vm = [np.exp(1j * k * np.arange(n)) / np.sqrt(n) for k in ks]
                for i, k in enumerate(ks):
                    j = min(range(n), key=lambda j: abs(np.exp(1j * ks[j]) - np.exp(-1j * k)))
                    if i < j:
                        X = (np.outer(Vm[i], Vm[i].conj()) - np.outer(Vm[j], Vm[j].conj())).reshape(-1)
                        Xc = V0.conj().T @ X
                        ok_in &= np.linalg.norm(V0.conj().T @ (R * (V0 @ Xc))) < 1e-10 \
                            and np.linalg.norm(V0 @ Xc - X) < 1e-10
                ok_dim &= ker.shape[1] == n // 2 and kmoat == n // 2
            else:
                ok_even &= ker.shape[1] > n // 2
check("odd N = 5..15, both rings, random locus profiles: P_k - P_-k in the kernel", ok_in)
check("odd N = 5..15: the kernel has dimension floor(N/2) exactly, read inside the moat", ok_dim)
check("control: even N = 4, 6, 8: the kernel exceeds floor(N/2)", ok_even)

print("R6  exact counts, rank modulo two primes against F140's floor")
rnd = random.Random(1)
ok = True
for n in range(3, 22):
    g = locus(n, rnd)
    jnum = rnd.randint(300, 900)
    for wrap in (1, -1):
        re, im = scaled(n, 1, wrap, True, jnum, g)
        ok &= nullity(re, im) == n // 2
check("ring corner, both rings, N = 3..21, a locus profile with distinct offsets: floor(N/2)", ok)
for n in ([3, 5, 7] if DEEP else [3, 5]):
    g = locus(n, random.Random(2))
    for label, wrap, fermi in (("spin ring", 1, False), ("periodic", 1, True), ("pi flux", -1, True)):
        tot, lw = rung_counts(n, wrap, fermi, 768, g)
        check(f"N={n} {label}, locus: every rung floor(N/2), all of it lowest weight at rung 1",
              tot == [n // 2] * (n - 1) and lw == [n // 2] + [0] * (n - 2), f"{tot} LW {lw}")
    tot, lw = rung_counts(n, 1, False, 768, [500] * n)
    check(f"N={n} spin ring, uniform point: every rung N - 1, all of it lowest weight at rung 1",
          tot == [n - 1] * (n - 1) and lw == [n - 1] + [0] * (n - 2), f"{tot} LW {lw}")

print("R7  the even-N reads the proof leaves open")
g4 = locus(4, random.Random(3))
tot, lw = rung_counts(4, 1, False, 768, g4)
check("N=4 spin ring, locus: [2, 4, 2]", tot == [2, 4, 2], f"{tot}")
tot, lw = rung_counts(4, -1, True, 768, g4)
check("N=4 pi-flux ring, locus: two lowest-weight modes at rung 2", lw == [2, 2, 0], f"LW {lw}")
tot, lw = rung_counts(6, 1, True, 768, [500] * 6)
check("N=6 periodic ring, uniform point: lowest weight at rung 2 (no anchor there)", lw[1] > 0, f"LW {lw}")


def float_lw(n, p, wrap, fermi, J, g):
    re, im = scaled(n, p, wrap, fermi, J, g)
    M = (re + 1j * im).astype(complex)
    return moat_nullity(np.vstack([M, psi(n, p).astype(complex)]))


a = float_lw(6, 2, 1, True, 768, [500] * 6)
b = float_lw(4, 2, -1, True, 768, g4)
check("from below: the same two lowest-weight reads by float SVD inside the moat",
      a == lw[1] and b == 2, f"N=6 uniform periodic rung 2: {a}, N=4 locus pi flux rung 2: {b}")

print("R8  the exceptional profiles: one offset pair at composite N")
from math import gcd
ok = True
for n in (9, 15, 21):
    for l in range(n // 2):
        g = [500] * n
        g[l], g[n - 1 - l] = 537, 463
        want = n // 2 + (gcd(2 * l + 1, n) - 1) // 2
        for wrap in (1, -1):
            for jnum in (311, 977):
                re, im = scaled(n, 1, wrap, True, jnum, g)
                ok &= nullity(re, im) == want
check("corner on a single-pair profile: floor(N/2) + (gcd(2l+1, N) - 1)/2, N = 9, 15, 21, every l", ok)


def bijects(pairs, nodes, want):
    return (len(pairs) == len(nodes) == want and all(len(p & nodes) == 1 for p in pairs)
            and all(any(m in p for p in pairs) for m in nodes))


def missed(n, l):
    """The indices m of the doubled momenta 2k = 2 pi m / N at which the offset's transform vanishes."""
    return {m for m in range(1, n) if m * (2 * l + 1) % n == 0}


def heis_nodes(n, l):
    """F157's node law at the seat l, (2l+1)m = N (mod 2N); the float read below tests it."""
    return {m for m in range(1, n) if (m * (2 * l + 1) - n) % (2 * n) == 0}


ok, xy_ok = True, True
for n in range(3, 200, 2):
    for l in range(n // 2):
        g = gcd(2 * l + 1, n)
        pairs = {frozenset((m, n - m)) for m in missed(n, l)}
        xy_nodes = {m for m in range(1, n + 1) if m * (l + 1) % (n + 1) == 0}
        ok &= bijects(pairs, heis_nodes(n, l), (g - 1) // 2)
        xy_ok &= bijects(pairs, xy_nodes, (g - 1) // 2)
check("F157's node law picks exactly one member of each missed pair, odd N to 199", ok)
check("control: the same bijection with the XY book's nodes fails", not xy_ok)
ok, cases = True, 0
for n in range(3, 32, 2):
    x = np.arange(n)
    lap = np.diag([1.0] + [2.0] * (n - 2) + [1.0]) - np.eye(n, k=1) - np.eye(n, k=-1)
    _, vecs = np.linalg.eigh(lap)                     # ascending: the mode m = 0..N-1 in order
    for l in range(n // 2):
        d = np.zeros(n); d[l], d[n - 1 - l] = 1.0, -1.0
        dhat = np.array([abs((d * np.exp(-2j * np.pi * m * x / n)).sum()) / n for m in range(n)])
        amp = np.abs(vecs[l, :])
        if np.any((dhat > 1e-10) & (dhat < 1e-4)) or np.any((amp > 1e-10) & (amp < 1e-4)):
            ok = False                                # a value inside the moat: no read
            continue
        zeros = {m for m in range(1, n) if dhat[m] <= 1e-10}
        blind = {m for m in range(1, n) if amp[m] <= 1e-10}
        ok &= zeros == missed(n, l)
        ok &= blind == heis_nodes(n, l)
        cases += gcd(2 * l + 1, n) > 1                # seats where a pair's two members are told apart
check(f"floats, N to 31, in a moat: the transform's zeros are the missed set, the chain's nodes at the seat "
      f"are F157's node law, {cases} seats with gcd(2l+1, N) > 1", ok and cases > 0)

print("R9  even N: the resonance as the staggered spin ladder, and what it leaves to the corner")


def bar(n, wrap, m):
    """the index of pi - k_m among the allowed momenta, at even N."""
    return (n // 2 - m) % n if wrap == 1 else (n // 2 - 1 - m) % n


def mode_ann(n, p, v):
    """sum_l conj(v_l) c_l from the p- to the (p-1)-particle sector, Slater signs as in psi."""
    hi, lo = list(combinations(range(n), p)), list(combinations(range(n), p - 1))
    ilo = {s: i for i, s in enumerate(lo)}
    C = np.zeros((len(lo), len(hi)), complex)
    for j, A in enumerate(hi):
        for l in A:
            A2 = tuple(x for x in A if x != l)
            C[ilo[A2], j] += np.conj(v[l]) * jw(A2, l)
    return C


ok_bar = all(abs(np.exp(1j * momenta(n, w)[bar(n, w, m)]) - np.exp(1j * (np.pi - momenta(n, w)[m]))) < 1e-12
             for n in range(4, 42, 2) for w in (1, -1) for m in range(n))
check("even N = 4..40, both rings: k -> pi - k maps the allowed momenta onto themselves", ok_bar)

worst, ctrl, ok_min, cases = 0.0, np.inf, True, []
for n, ls in ((4, (1, 2)), (6, (1, 2, 3)), (8, (1, 2, 3) if DEEP else (1, 2))):
    for wrap in (1, -1):
        V = [np.exp(1j * k * np.arange(n)) / np.sqrt(n) for k in momenta(n, wrap)]
        for l in ls:
            st, Q = lw_v0(n, l, wrap)
            d = Q.shape[1]
            C = [mode_ann(n, l, v) for v in V]                 # l -> l-1
            Cup = [mode_ann(n, l + 1, v) for v in V] if l < n else None
            Dh = np.array([[len(set(x) & set(y)) for y in st] for x in st], float).reshape(-1)
            lhs = n * Q.conj().T @ (Dh[:, None] * Q)
            for stag in (True, False):
                b_ = (lambda m: bar(n, wrap, m)) if stag else (lambda m: m)
                Y = lambda a, b: np.kron(C[a], C[b].conj())
                rhs = l * (l + 1) * np.eye(d, dtype=complex)
                S = sum(np.kron(C[b_(b)], Cup[b].T) for b in range(n)) @ Q
                rhs -= S.conj().T @ S
                for a in range(n):
                    G = (Y(a, a) + Y(b_(a), b_(a))) @ Q
                    rhs -= 0.5 * G.conj().T @ G
                    if b_(a) != a:
                        T = Y(a, b_(a)) @ Q
                        rhs -= T.conj().T @ T
                res = float(np.max(np.abs(lhs - rhs)))
                if stag:
                    worst = max(worst, res / (np.finfo(float).eps * Q.shape[0]))
                elif l >= 2:
                    ctrl = min(ctrl, res)
            K = np.array([[len(set(x) - set(y)) for y in st] for x in st], float).reshape(-1)
            ev = np.linalg.eigvalsh(Q.conj().T @ (K[:, None] * Q))
            ok_min &= abs(ev[0] - l * (n - 1 - l) / n) < 1e-10
            if l >= 2:
                cases.append((n, ev[0]))
check("N = 4, 6 every rung, N = 8 rungs 1, 2 (3 under --deep), both rings: N*Dhat = l(l+1) - |S v|^2 "
      "- (1/2)sum|g_a + g_abar|^2 - sum|Y_a,abar|^2 on LW cap V0, residual / (eps * dim) below 1e3", worst < 1e3, f"worst {worst:.1f}")
check("control: with the ladder unstaggered (abar -> a) it misses by more than 0.1 on every rung l >= 2 "
      "(at l = 1 it holds by accident on the periodic ring at N = 4)", ctrl > 0.1,
      f"smallest miss {ctrl:.2f}")
check("min K = l(N-1-l)/N on LW cap V0, attained, every rung read; above the seed that is 1 at N = 6 "
      "and above 1 at N = 8", ok_min and any(n == 8 for n, _ in cases),
      ", ".join(f"N={n}: {e:.4f}" for n, e in cases))

ok_u, ok_foot = True, True
for n in range(4, 22, 2):
    for wrap in (1, -1):
        re, im = scaled(n, 1, wrap, True, 577, [500] * n)
        want = 3 * n // 2 - 2 if wrap == 1 else 3 * n // 2 - 1
        ok_u &= nullity(re, im) == want
        hw = h_ring(n, wrap)
        comm = (np.kron(hw, np.eye(n, dtype=np.int64)) - np.kron(np.eye(n, dtype=np.int64), hw.T))
        diag = np.zeros((n, n * n), dtype=np.int64)
        for a in range(n):
            diag[a, a * n + a] = 1
        ok_foot &= nullity(comm, np.zeros_like(comm), diag) == want
check("uniform point, corner, N = 4..20: 3N/2 - 2 on the periodic ring, 3N/2 - 1 on the pi-flux ring",
      ok_u)
check("the same numbers as the commutant elements with no diagonal, frozen at every J", ok_foot)
cnt = []
for p in (1, 2):
    re, im = scaled(8, p, 1, False, 577, [500] * 8)
    cnt.append(nullity(re, im))
check("spin ring, N = 8, uniform point: rung 1 carries 10, rung 2 carries 11 (the two corners alternate)",
      cnt == [10, 11], f"{cnt}")

ok_orb = True
for n in range(4, 16, 2):
    for wrap in (1, -1):
        ks = momenta(n, wrap)
        neg = lambda m: min(range(n), key=lambda j: abs(np.exp(1j * ks[j]) - np.exp(-1j * ks[m])))
        shift = lambda m: min(range(n), key=lambda j: abs(np.exp(1j * ks[j]) + np.exp(1j * ks[m])))
        free = {m for m in range(n) if neg(m) != m}
        orbits, seen = 0, set()
        for m in free:
            if m not in seen:
                orbits += 1
                seen |= {m, neg(m), shift(m), neg(shift(m))}
        h = h_ring(n, wrap).astype(float)
        V0 = null(np.kron(h, np.eye(n)) - np.kron(np.eye(n), h.T))
        d = rng.uniform(-0.3, 0.3, n // 2)
        g = np.ones(n)
        for i in range(n // 2):
            g[i], g[n - 1 - i] = 1 + d[i], 1 - d[i]
        R = np.array([[(2 * (g[a] + g[b]) if a != b else 0) - 4 for b in range(n)]
                      for a in range(n)]).reshape(-1)
        k = moat_nullity(V0.conj().T @ (R[:, None] * V0), lo=1e-11, hi=1e-6)
        ok_orb &= k == n // 2 + orbits
check("read: at even N = 4..14 the large-J kernel exceeds N/2 by one per orbit of k -> -k, k -> k + pi",
      ok_orb)

ok_c, ok_split = True, True
rnd9 = random.Random(3)
for n in range(4, 16, 2):
    U = np.zeros((n * n, n * n), dtype=np.int64)
    for a in range(n):
        for b in range(n):
            U[b * n + a, a * n + b] = (-1) ** (a + b)
    g = locus(n, rnd9)
    for wrap in (1, -1):
        re, im = scaled(n, 1, wrap, True, 311, g)
        ok_c &= not np.any(U @ re - re @ U) and not np.any(U @ im - im @ U)
        I = np.eye(n * n, dtype=np.int64)
        ok_split &= nullity(re, im, U - I) == 0 and nullity(re, im, U + I) == n // 2
check("even N = 4..14, both rings, integer locus profiles: U v(a,b) = (-1)^(a+b) v(b,a) commutes exactly "
      "with the corner", ok_c)
check("there U = +1 (the populations) carries no frozen mode and U = -1 carries N/2", ok_split)

ok_ur = True
for n in range(4, 14, 2):
    U = np.zeros((n * n, n * n), dtype=np.int64)
    for a in range(n):
        for b in range(n):
            U[b * n + a, a * n + b] = (-1) ** (a + b)
    for wrap in (1, -1):
        re, im = scaled(n, 1, wrap, True, 311, [500] * n)
        ok_ur &= nullity(re, im, U - np.eye(n * n, dtype=np.int64)) > 0
check("uniform point, even N = 4..12, both rings: U = +1 reaches the root, so it is no anchor for the populations", ok_ur)

import sympy as sp
from sympy.polys.matrices import DomainMatrix
Js = sp.symbols("J")
g6, n = [60, 31, 64, 16, 49, 20], 6
gb = sp.Rational(sum(g6), n)
h6 = h_ring(n, 1)
M6 = sp.zeros(n * n, n * n)
for a in range(n):
    for b in range(n):
        r = a * n + b
        M6[r, r] += 4 * gb if a == b else 4 * gb - 2 * (g6[a] + g6[b])
        for c in range(n):
            if h6[a, c]:
                M6[r, c * n + b] += -sp.I * Js * int(h6[a, c])
            if h6[c, b]:
                M6[r, a * n + c] += sp.I * Js * int(h6[c, b])
cols = []
for a in range(n):
    v = sp.zeros(n * n, 1); v[a * n + a] = 1; cols.append(v)
for a in range(n):
    for b in range(a + 1, n):
        v = sp.zeros(n * n, 1); v[a * n + b] = 1; v[b * n + a] = (-1) ** (a + b); cols.append(v)
Bp = sp.Matrix.hstack(*cols)
d6 = sp.factor(DomainMatrix.from_Matrix(Bp.T * M6 * Bp).det().as_expr())
ratio = sp.simplify(d6 / (Js ** 10 * (13 * Js ** 2 - 35816) ** 2))
check("N = 6, periodic, rates (60, 31, 64, 16, 49, 20): det of the populations' sector is c*J^10*(13J^2 - 35816)^2, "
      "a root at J^2 = 35816/13", ratio.is_number and ratio != 0, f"{d6}")

print()
print("ring frozen band gate:", "ALL GREEN" if not FAILS else f"{len(FAILS)} FAILED: {FAILS}")
sys.exit(1 if FAILS else 0)
