"""Gate for docs/proofs/PROOF_RING_FROZEN_BAND.md: the XY ring's frozen band at odd N, and at even N
above the corner and in the corners of both fixed-flux rings.

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
  R10 even N, the populations' sector U = +1 at second order from the uniform point (Proposition 7):
      its uniform kernel is the orbits {k, pi - k} minus one, modulo two primes at N = 4..20, an upper
      bound that the orbit sums with zero diagonal meet;
      the second-order matrix W2 by the closed form (Cauchy rows from the nonresonant offsets, one
      class sum per resonant one, X_r in its displayed form, no matrix inverse) equals W2 through the pseudo-inverse of the
      uniform corner at N = 4..22, both rings, random couplings and profiles, the deviation gated
      against the size of a rounding through that route, eps |K1|^2 |K0^+|; it is real, symmetric and
      (gbar > 0) negative semidefinite through the pseudo-inverse, of nullity 0 on the periodic ring and 1 on the
      pi-flux ring, read against eps |K1|^2 |K0^+| inside a moat; a control with the self-paired orbit mis-weighted must miss; in integers to N = 200 the
      resonant offsets cover the classes and the nonresonant ones match them in number on the periodic
      ring and fall one short on the pi-flux ring; the one null direction there is sum_k P_k / sin k,
      v(a,b) = sgn(b - a) for a + b odd; the reflection l -> N-1-l maps K0 to itself, K1 to -K1 and that
      direction to minus itself, so <w, W3 w> = 0; read only, the fourth is real, positive
      at gbar > 0 and proportional to 1/(J^2 gbar) at N = 6..12, meeting the smallest eigenvalue at offset
      1e-2; and, read off the locus on zero-mean profiles, the same second-order nullities.
  R11 even N, the pi-flux ring's populations' sector at fourth order (Proposition 8): the rotation
      identity {D, j^-1} = [h, [Lap^-1 D, Lam]], Lam = sum_k cot k P_k, at N = 4..24 on random zero-mean
      offsets, to the size of a rounding, with the inverse Laplacian dropped as the control that must miss;
      in the corner, K1 w = K0 y with y = -(2i/J)[phi, Lam] free of populations and of the uniform kernel;
      on the locus the second-order form couples y to no t-direction, and off the locus it does (the
      control); c4 by its definition equals the population cost (4 gbar / |w|^2) min_z |diag(R0^dag K1 (y + z))|^2
      over Z+ and over Z' at N = 4..20 on random locus profiles and couplings; read, it equals the closed form
      through tau, and y pays nothing in the nonresonant sectors; the populations of y and of z in Z+ at the
      resonant q meet their formulas, and mu_C is real on the locus; on the single-pair profile e_0 - e_(N-1)
      c4 by its definition and by the closed form meet the formula by hand at N = 4..24 with the sign of gbar,
      and the smallest eigenvalue of the U = +1 block is c4 e^4 (1 + a e^2) at N = 4..12, a route with no
      pseudo-inverse; and, read, c4 drops to rounding where tau changes sign on the locus, away from its poles,
      at N = 6, 8, 10.
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
check("uniform point, even N = 4..12, both rings, one coupling each: U = +1 reaches the root, so the uniform "
      "point anchors it only at second order (R10)", ok_ur)

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

print("R10 even N: the populations' sector at second order from the uniform point")


def corner_f(n, wrap, J, gbar, D):
    """the corner (1,1) of L + 4*gbar in floats, cell (a, b) at a*n + b, rates gbar + D."""
    h = h_ring(n, wrap).astype(float)
    I = np.eye(n)
    g = gbar + np.asarray(D, float)
    rate = np.array([[4 * gbar if a == b else 4 * gbar - 2 * (g[a] + g[b]) for b in range(n)]
                     for a in range(n)]).reshape(-1)
    return np.diag(rate).astype(complex) - 1j * J * (np.kron(h, I) - np.kron(I, h))


def sym_basis(n):
    """orthonormal basis of U = +1: v(a,b) = (-1)^(a+b) v(b,a)."""
    cols = []
    for a in range(n):
        v = np.zeros(n * n); v[a * n + a] = 1; cols.append(v)
    for a in range(n):
        for b in range(a + 1, n):
            v = np.zeros(n * n); v[a * n + b] = 1; v[b * n + a] = (-1) ** (a + b)
            cols.append(v / np.sqrt(2))
    return np.array(cols).T


def orbits_pi(n, wrap):
    """A-orbits O = {k, pi - k} of total momentum pi, momenta k = pi*kap/n, kap mod 2n."""
    kap = [2 * m + (0 if wrap == 1 else 1) for m in range(n)]
    out = []
    for c in kap:
        O = tuple(sorted({c % (2 * n), (n - c) % (2 * n)}))
        if O not in out:
            out.append(O)
    return out


def w2_routes(n, wrap, J, gbar, amp, fixed_weight=1):
    """W2 on the uniform kernel of U = +1, two routes: the site route through a pseudo-inverse,
    and the closed form of Proposition 7 (no inverse). Coordinates: z = sum_O s_O T_O,
    T_O = sum_{k in O} phi_k phi_{pi-k}^T, kernel sum_O |O| s_O = 0."""
    B = sym_basis(n)
    Sig = np.diag([(-1.0) ** a for a in range(n)])
    phi = lambda c: np.exp(1j * np.pi * c / n * np.arange(n)) / np.sqrt(n)
    orbs = orbits_pi(n, wrap)
    T = []
    for O in orbs:
        S = sum(np.outer(phi(c), phi((n - c) % (2 * n))) for c in O)
        T.append(B.T @ (S @ Sig).reshape(-1))
    T = np.array(T).T
    wts = np.array([len(O) for O in orbs], float)
    Ck = np.linalg.svd(wts[None, :])[2][1:].T
    Zs = T @ Ck
    a = np.arange(n) + 0.5
    D = sum(amp[j - 1] * np.sin(2 * np.pi * j / n * a) for j in range(1, n // 2 + 1))
    K0 = B.T @ corner_f(n, wrap, J, gbar, np.zeros(n)) @ B
    K1 = B.T @ (corner_f(n, wrap, J, gbar, D) - corner_f(n, wrap, J, gbar, np.zeros(n))) @ B
    W_site = -Zs.conj().T @ K1 @ np.linalg.pinv(K0, rcond=1e-12) @ K1 @ Zs
    c4 = 4 * gbar
    eps = lambda c: 2 * np.cos(np.pi * c / n)
    kap = [2 * m + (0 if wrap == 1 else 1) for m in range(n)]
    W = np.zeros((len(orbs), len(orbs)), complex)
    for j in range(1, n // 2 + 1):
        rr = 2 * j                                         # r = pi*rr/n
        dh2 = (amp[j - 1] / 2) ** 2 if j < n // 2 else amp[j - 1] ** 2   # |dhat(r)|^2
        mult = 2 if j < n // 2 else 1                       # r and -r, or r = pi once
        sr = np.sin(np.pi * j / n)
        res = any(abs(np.sin(np.pi * c / n) ** 2 - sr ** 2) < 1e-12 for c in kap)
        if not res:
            if j == n // 2:
                continue
            X = np.array([len(O) * np.cos(np.pi * j / n) * np.sin(np.pi * O[0] / n)
                          / (2j * J * np.sqrt(n) * sr * (np.sin(np.pi * O[0] / n) ** 2 - sr ** 2)) for O in orbs])
            Xsum = np.array([len(O) / np.sqrt(n) * (1 / (-1j * J * (eps(O[0]) - eps(O[0] - rr)))
                                                    + 1 / (-1j * J * (eps(O[0] + rr) - eps(O[0]))))
                             for O in orbs])
            X_OK.append(np.abs(X - Xsum).max() <= 1e-12 * np.abs(Xsum).max())
            W += mult * 4 * c4 * dh2 * np.outer(X, X)
        else:
            half = j
            Cl = {half % (2 * n), (n - half) % (2 * n), (-half) % (2 * n), (n + half) % (2 * n)}
            y = np.array([np.sqrt(2) * (fixed_weight if len(O) == 1 else 1) if set(O) & Cl else 0.0
                          for O in orbs])
            W += -mult * 4 * dh2 / (c4 * 2 / n) * np.outer(y, y)
    W_th = Ck.T @ W @ Ck
    LOW_OK.append(np.abs(K0 @ Zs).max() < 1e-12 and np.linalg.matrix_rank(Zs) == len(orbs) - 1)
    return W_site, W_th, Zs, K0, K1, B


X_OK, LOW_OK = [], []
rnd10 = np.random.default_rng(10)
worst, ok_neg, ok_null, ok_ctl, gap = 0.0, True, True, True, 1.0
nulls = {}
for n in range(4, 23, 2):
    for wrap in (1, -1):
        J, gbar = rnd10.uniform(0.5, 2), rnd10.uniform(0.5, 2)
        amp = rnd10.uniform(0.5, 1.5, n // 2) * rnd10.choice((-1, 1), n // 2)
        Ws, Wt, Zs, K0, K1, B = w2_routes(n, wrap, J, gbar, amp)
        size = np.linalg.norm(K1, 2) ** 2 * np.linalg.norm(np.linalg.pinv(K0, rcond=1e-12), 2)
        worst = max(worst, np.abs(Ws - Wt).max() / (np.finfo(float).eps * size))
        tol = 1e3 * np.finfo(float).eps * size
        ev = np.linalg.eigvalsh((Ws + Ws.conj().T) / 2) if Ws.size else np.zeros(0)
        ok_neg &= np.abs(Ws.imag).max() < tol and np.abs(Ws - Ws.T).max() < tol and (ev.size == 0 or ev.max() < tol)
        k0 = int(np.sum(np.abs(ev) < tol)); mid = np.any((np.abs(ev) >= tol) & (np.abs(ev) < 1e-6 * size))
        gap = min(gap, np.abs(ev)[np.abs(ev) >= tol].min() / size) if np.any(np.abs(ev) >= tol) else gap
        nulls[(n, wrap)] = k0
        ok_null &= (not mid) and k0 == (0 if wrap == 1 else 1)
        if n % 4 == 0 and wrap == 1:
            _, Wbad, *_ = w2_routes(n, wrap, J, gbar, amp, fixed_weight=2)
            ok_ctl &= np.abs(Ws - Wbad).max() / (np.finfo(float).eps * size) > 1e8
check("even N = 4..22, both rings, random J, gbar and locus offsets with every Fourier component: W2 by the "
      "closed form equals W2 through the pseudo-inverse, the deviation below 1e2 eps times |K1|^2 |K0^+|, "
      "the size of a rounding through that route; the displayed form of X_r equals its sum 1/a_k + 1/a_(k+r)",
      worst < 1e2 and all(X_OK), f"worst {worst:.2f}")
check("W2 through the pseudo-inverse is real, symmetric and, at gbar > 0, negative semidefinite, to 1e3 eps "
      "|K1|^2 |K0^+|", ok_neg)
check("its nullity, read against that scale inside a moat (nothing between 1e3 eps and 1e-6 of it): 0 on the "
      "periodic ring, 1 on the pi-flux ring", ok_null, f"smallest nonzero |eigenvalue| / scale {gap:.1e}; {nulls}")
check("control: weighting the self-paired orbit {pi/2} like a pair breaks the agreement at 4 | N, periodic", ok_ctl)

ok_k = True
for n in range(4, 21, 2):
    U = np.zeros((n * n, n * n), dtype=np.int64)
    for a in range(n):
        for b in range(n):
            U[b * n + a, a * n + b] = (-1) ** (a + b)
    for wrap in (1, -1):
        want = len(orbits_pi(n, wrap)) - 1
        for jn in (311, 7):
            re, im = scaled(n, 1, wrap, True, jn, [500] * n)
            ok_k &= nullity(re, im, U - np.eye(n * n, dtype=np.int64)) == want
check("uniform point, N = 4..20, both rings, two couplings: U = +1 carries exactly the orbits {k, pi - k} "
      "minus one (nullity modulo two primes, an upper bound, met by the orbit sums with zero diagonal, "
      "which K0 annihilates in floats, N = 4..22)", ok_k and all(LOW_OK))

ok_cnt = True
for n in range(4, 202, 2):
    for wrap in (1, -1):
        kap = [2 * m + (0 if wrap == 1 else 1) for m in range(n)]
        cls = {min(c % n, (-c) % n) for c in kap if c % n}          # classes with sin k != 0, by |k| mod pi
        nonres = [j for j in range(1, n // 2) if j % 2 == (1 if wrap == 1 else 0)]
        res = [j for j in range(1, n // 2 + 1) if j % 2 == (0 if wrap == 1 else 1)]
        ok_cnt &= len(nonres) == len(cls) - (0 if wrap == 1 else 1)
        ok_cnt &= sorted(j % n for j in res) == sorted(cls)
check("integers, even N = 4..200: the classes with sin k != 0 are exactly the resonant r/2, and the "
      "nonresonant r with cos(r/2) != 0 match them in number on the periodic ring and fall one short on the "
      "pi-flux ring", ok_cnt)

ok_t = True
for n in range(4, 23, 2):
    J, gbar = rnd10.uniform(0.5, 2), rnd10.uniform(0.5, 2)
    amp = rnd10.uniform(0.5, 1.5, n // 2)
    Ws, Wt, Zs, K0, K1, B = w2_routes(n, -1, J, gbar, amp)
    v = np.zeros((n, n))
    for a in range(n):
        for b in range(n):
            if (a + b) % 2:
                v[a, b] = np.sign(b - a)
    ks = momenta(n, -1)
    vk = sum(np.outer(np.exp(1j * k * np.arange(n)), np.exp(-1j * k * np.arange(n))) / n / np.sin(k) for k in ks)
    z = B.T @ v.reshape(-1); z = z / np.linalg.norm(z)
    ok_t &= np.abs(vk / vk[0, 1] - v).max() < 1e-12
    ok_t &= np.abs(K0 @ z).max() < 1e-12
    W2z = -K1 @ np.linalg.pinv(K0, rcond=1e-12) @ K1 @ z
    ok_t &= np.abs(Zs.conj().T @ W2z).max() < 1e-12 * max(np.abs(Ws).max(), 1.0)
check("pi-flux ring, N = 4..22: the one null direction is v(a,b) = sgn(b - a) for a + b odd, which is "
      "sum_k P_k / sin k, the inverse of the particle current", ok_t)

ok_l4, ok_sc, ok_3 = True, True, True
for n in (6, 8, 10, 12):
    D = np.zeros(n); d = rnd10.standard_normal(n // 2)
    for l in range(n // 2):
        D[l], D[n - 1 - l] = d[l], -d[l]
    vals = {}
    for (J, gbar) in ((1.0, 1.0), (2.0, 1.0), (1.0, 2.0), (3.0, 1.0), (1.0, 4.0)):
        B = sym_basis(n)
        K0 = B.T @ corner_f(n, -1, J, gbar, np.zeros(n)) @ B
        K1 = B.T @ (corner_f(n, -1, J, gbar, D) - corner_f(n, -1, J, gbar, np.zeros(n))) @ B
        R = np.linalg.pinv(K0, rcond=1e-12)
        Z = null(K0)
        P = Z @ Z.conj().T
        v = np.zeros((n, n))
        for a in range(n):
            for b in range(n):
                if (a + b) % 2:
                    v[a, b] = np.sign(b - a)
        t = B.T @ v.reshape(-1); t = t / np.linalg.norm(t)
        W2 = -P @ K1 @ R @ K1 @ P
        W3 = P @ K1 @ R @ K1 @ R @ K1 @ P
        W4 = -P @ K1 @ R @ K1 @ R @ K1 @ R @ K1 @ P
        l3 = t.conj() @ W3 @ t
        l4 = t.conj() @ W4 @ t - t.conj() @ W3 @ np.linalg.pinv(W2, rcond=1e-10) @ W3 @ t
        ev = np.linalg.eigvals(K0 + 1e-2 * K1)
        e0 = ev[np.argmin(abs(ev))] / 1e-8
        Rf = np.zeros((n * n, n * n))
        for x in range(n):
            for y in range(n):
                Rf[(n - 1 - x) * n + (n - 1 - y), x * n + y] = 1
        Rb = B.T @ Rf @ B
        ok_3 &= np.abs(Rf @ B - B @ Rb).max() < 1e-13            # the reflection keeps U = +1
        ok_3 &= np.abs(Rb @ K0 - K0 @ Rb).max() < 1e-13 and np.abs(Rb @ K1 + K1 @ Rb).max() < 1e-12 \
            and np.abs(Rb @ t + t).max() < 1e-13
        ok_3 &= abs(l3) < 1e-10 * np.abs(W3).max()
        ok_l4 &= abs(l4.imag) < 1e-10 * abs(l4) and l4.real > 0 and abs(e0 / l4 - 1) < 1e-2
        vals[(J, gbar)] = l4.real
    ok_sc &= all(abs(vals[(J, g)] * J * J * g / vals[(1.0, 1.0)] - 1) < 1e-8 for (J, g) in vals)
check("pi-flux ring, N = 6..12, one locus profile each: the reflection maps K1 to -K1, K0 and U to themselves and "
      "that direction to minus itself, and <w, W3 w> = 0", ok_3)
check("read: the fourth-order coefficient there is real and positive and meets the smallest eigenvalue of the "
      "corner's U = +1 block at offset 1e-2 over its fourth power within 1 per cent", ok_l4)
check("read: it scales as 1/(J^2 gbar), at J/gbar = 1/4, 1/2, 1, 2, 3", ok_sc)

ok_off = True
for n in range(4, 17, 2):
    for wrap in (1, -1):
        B = sym_basis(n)
        D = rnd10.standard_normal(n); D -= D.mean()
        K0 = B.T @ corner_f(n, wrap, 1.1, 1.0, np.zeros(n)) @ B
        K1 = B.T @ corner_f(n, wrap, 1.1, 1.0, D) @ B - K0
        Z = null(K0)
        W2 = -Z.conj().T @ K1 @ np.linalg.pinv(K0, rcond=1e-12) @ K1 @ Z
        k = moat_nullity(W2, lo=1e-11, hi=1e-6)
        ok_off &= k == (0 if wrap == 1 else 1)
check("read, off the locus, N = 4..16, random zero-mean profiles: the same second-order nullities, 0 periodic and "
      "1 pi flux, inside a moat (the proof does not use the locus)", ok_off)



# ---------------- R11: the pi-flux ring's populations' sector at fourth order ----------------

print()
print("R11 even N, pi-flux ring: the inverse current is only turned at first order, and the trace keeps it from the root")


def pflux(n):
    """plane waves of the pi-flux ring, the inverse current j^-1, Lam = (1/2) h j^-1 = sum cot k P_k,
    and the periodic lattice Laplacian on the sites."""
    ks = np.array(momenta(n, -1))
    Ph = np.exp(1j * np.outer(np.arange(n), ks)) / np.sqrt(n)
    jinv = Ph @ np.diag(1 / np.sin(ks)) @ Ph.conj().T
    Lam = Ph @ np.diag(1 / np.tan(ks)) @ Ph.conj().T
    Lap = -2 * np.eye(n) + np.roll(np.eye(n), 1, 0) + np.roll(np.eye(n), -1, 0)
    return ks, Ph, jinv, Lam, Lap


def turn(D, Lam, Lap, J):
    """the first-order dressing y = -(2i/J) [phi, Lam], phi = Lap^-1 D (zero mean): i[phi, Lam] is the first-order
    change of Lam under the site phases e^(i phi), a turn about Z."""
    phi = np.linalg.lstsq(Lap, D, rcond=None)[0]
    return -2j / J * (np.diag(phi) @ Lam - Lam @ np.diag(phi))


def locus_profile(n, rnd):
    D = np.zeros(n); d = rnd.standard_normal(n // 2)
    for l in range(n // 2):
        D[l], D[n - 1 - l] = d[l], -d[l]
    return D


def fourth(n, J, gbar, D):
    """c4 by its definition (Proposition 8), and the population cost minimized over Z' (the kernel without w)
    and over Z+ (the traceless combinations of class sums)."""
    ks, Ph, jinv, Lam, Lap = pflux(n)
    B = sym_basis(n)
    K0 = B.T @ corner_f(n, -1, J, gbar, np.zeros(n)) @ B
    K1 = B.T @ corner_f(n, -1, J, gbar, D) @ B - K0
    R = np.linalg.pinv(K0, rcond=1e-12)
    Z = null(K0)
    P = Z @ Z.conj().T
    w = B.T @ jinv.reshape(-1)
    nw = np.linalg.norm(w); t = w / nw
    W2 = -P @ K1 @ R @ K1 @ P
    W3 = P @ K1 @ R @ K1 @ R @ K1 @ P
    W4 = -P @ K1 @ R @ K1 @ R @ K1 @ R @ K1 @ P
    c4 = t.conj() @ W4 @ t - t.conj() @ W3 @ np.linalg.pinv(W2, rcond=1e-10) @ W3 @ t
    y = B.T @ turn(D, Lam, Lap, J).reshape(-1)
    G = (R.conj().T @ K1)[:n, :]            # the first n coordinates of sym_basis are the populations
    def cost(Zs):
        a0 = G @ y
        if Zs.shape[1]:
            Am = G @ Zs
            a0 = a0 + Am @ np.linalg.lstsq(Am, -a0, rcond=None)[0]
        return 4 * gbar * np.linalg.norm(a0) ** 2 / nw ** 2
    Zp = Z @ null((t.conj() @ Z)[None, :])  # the uniform kernel without w, whose response is exactly zero
    cls = classes(n)
    S = np.array([B.T @ (Ph[:, c] @ Ph[:, c].conj().T).reshape(-1) for c in cls]).T
    Zplus = S @ null(np.array([[len(c) for c in cls]], float)) if len(cls) > 1 else np.zeros((len(y), 0))
    return dict(c4=c4, cost=cost(Zp), cost_plus=cost(Zplus), K0=K0, K1=K1, R=R, Z=Z, w=w, y=y, B=B, G=G,
                ks=ks, Ph=Ph, Zplus=Zplus, cls=cls, S=S)


def classes(n):
    """the classes {k, -k, pi + k, pi - k} of the pi-flux momenta k = pi (2m + 1) / n, as index lists."""
    out, seen = [], set()
    for m in range(n):
        if m in seen:
            continue
        kap = 2 * m + 1
        C = sorted({((c % (2 * n)) - 1) // 2 for c in (kap, -kap, n + kap, n - kap)})
        seen |= set(C)
        out.append(C)
    return out


def theta(dh, n, k):
    """Theta(2k) = sum_p g(k - p) g(k + p), g(u) = dhat(u) / sin(u/2), p over the pi-flux momenta."""
    G = 0.0
    for p in momenta(n, -1):
        u, v = k - p, k + p
        mu, mv = int(round(u * n / (2 * np.pi))) % n, int(round(v * n / (2 * np.pi))) % n
        if mu and mv:
            G += dh[mu] * dh[mv] / (np.sin(u / 2) * np.sin(v / 2))
    return G


def mu(dh, n, j, J):
    """mu_C = -i Theta(q) / (J sin(q/2) dhat(q)) at the resonant q = 2 pi j / n."""
    q = 2 * np.pi * j / n
    return -1j * theta(dh, n, q / 2) / (J * np.sin(q / 2) * dh[j])


def c4_closed(D, J, gbar):
    """c4 = (2 / (N J^2 gbar)) |tau|^2 / sum_C |O_C|^2 / omega_C, tau = sum_C |O_C| Theta / (sin(q_C/2) dhat(q_C)),
    one class per resonant pair {j, N - j}, |O_C| = 2, and the class at q = pi, |O_C| = 1."""
    n = len(D); dh = np.fft.fft(D) / n
    seen, tau, s2 = set(), 0.0, 0.0
    for j in range(1, n, 2):
        if j in seen:
            continue
        pair = {j, n - j}; seen |= pair
        osz = 1.0 if len(pair) == 1 else 2.0
        q = 2 * np.pi * j / n
        tau += osz * theta(dh, n, q / 2) / (np.sin(q / 2) * dh[j])
        s2 += osz ** 2 / sum(abs(dh[i]) ** 2 for i in pair)
    return 2 * abs(tau) ** 2 / s2 / (n * J * J * gbar), tau


def c4_pair(n, J, gbar):
    """the single-pair profile e_0 - e_(N-1), by hand: tau = -(2i(N-2)/N) sum_C |O_C| / sin^2(q_C/2),
    omega_C = (4/N^2) (number of resonant q in C) sin^2(q_C/2)."""
    s1 = s2 = 0.0; seen = set()
    for j in range(1, n, 2):
        if j in seen:
            continue
        pair = {j, n - j}; seen |= pair
        osz = 1.0 if len(pair) == 1 else 2.0
        sn2 = np.sin(np.pi * j / n) ** 2
        s1 += osz / sn2
        s2 += osz ** 2 / (len(pair) * 4 / n ** 2 * sn2)
    return 2 * (2 * (n - 2) / n) ** 2 * s1 ** 2 / s2 / (n * J * J * gbar)


rnd11 = np.random.default_rng(11)
eps = np.finfo(float).eps

# the rotation identity, in single-particle space
worst, ctrl = 0.0, np.inf
for n in range(4, 25, 2):
    ks, Ph, jinv, Lam, Lap = pflux(n)
    h = h_ring(n, -1).astype(float)
    for _ in range(3):
        D = rnd11.standard_normal(n); D -= D.mean()
        phi = np.linalg.lstsq(Lap, D, rcond=None)[0]
        lhs = np.diag(D) @ jinv + jinv @ np.diag(D)
        X = np.diag(phi) @ Lam - Lam @ np.diag(phi)
        scale = eps * n * n * np.abs(D).max() * np.abs(jinv).max() * np.abs(Lam).max()
        worst = max(worst, np.abs(lhs - (h @ X - X @ h)).max() / scale)
        Xc = np.diag(D) @ Lam - Lam @ np.diag(D)              # control: without the inverse Laplacian
        ctrl = min(ctrl, np.abs(lhs - (h @ Xc - Xc @ h)).max() / np.abs(lhs).max())
check("pi-flux ring, N = 4..24, random zero-mean offsets: {D, j^-1} = [h, [Lap^-1 D, Lam]], Lam = sum_k cot k P_k, "
      "to the size of a rounding (worst / eps-scale < 1e2); without the inverse Laplacian it misses by > 0.1",
      worst < 1e2 and ctrl > 0.1, f"worst {worst:.2g}, control {ctrl:.2g}")

# in the corner: K1 w = K0 y, y population-free and outside the uniform kernel
worst2, pop2 = 0.0, 0.0
for n in range(4, 17, 2):
    for J, gbar in ((1.0, 1.0), (1.7, 0.6)):
        D = locus_profile(n, rnd11)
        f = fourth(n, J, gbar, D)
        K0, K1, w, y, Z = f["K0"], f["K1"], f["w"], f["y"], f["Z"]
        sc = eps * n * n * np.abs(K0).max() * np.linalg.norm(y)
        worst2 = max(worst2, np.abs(K1 @ w - K0 @ y).max() / sc,
                     np.abs(Z.conj().T @ y).max() / (eps * n * n * np.linalg.norm(y)))
        pop2 = max(pop2, np.abs(y[:n]).max())
check("the corner's U = +1, N = 4..16, locus profiles, two couplings: K1 w = K0 y with w = j^-1 and "
      "y = -(2i/J)[phi, Lam], to the size of a rounding (worst / eps-scale < 1e2), y with no part in the uniform "
      "kernel to that size and with populations exactly 0", worst2 < 1e2 and pop2 == 0.0,
      f"worst {worst2:.2g}, populations {pop2}")

# the reflection decouples y from the t-directions on the locus; off the locus it does not
ok3, ctl3 = True, np.inf
for n in range(6, 17, 2):
    for lab in ("locus", "off"):
        if lab == "locus":
            D = locus_profile(n, rnd11)
        else:
            D = rnd11.standard_normal(n); D -= D.mean()
        f = fourth(n, 1.3, 0.7, D)
        Xf = -f["K1"] @ f["R"] @ f["K1"]
        tds = [f["B"].T @ (f["Ph"] @ np.diag(np.sin(m * f["ks"])) @ f["Ph"].conj().T).reshape(-1)
               for m in range(1, n, 2)]
        y = f["y"]
        cr = max(abs(y.conj() @ Xf @ td) / (np.linalg.norm(y) * np.linalg.norm(td) * np.abs(Xf).max()) for td in tds)
        if lab == "locus":
            ok3 &= cr < 1e3 * eps * n * n
        else:
            ctl3 = min(ctl3, cr)
check("N = 6..16: on the locus the second-order form couples y to no t-direction (sum_k f(k) P_k, f odd), "
      "below 1e3 eps N^2 of its scale; off the locus the coupling is of order one (> 1e-3), the control",
      ok3 and ctl3 > 1e-3, f"off-locus min {ctl3:.2g}")

# c4 is the population cost, over Z+ as displayed and over Z'; the closed form and the nonresonant silence (read)
dev4 = dev5 = nonres = 0.0
for n in range(4, 21, 2):
    for _ in range(3):
        J, gbar = rnd11.uniform(0.5, 2.0), rnd11.uniform(0.3, 2.0)
        D = locus_profile(n, rnd11)
        f = fourth(n, J, gbar, D)
        c4 = f["c4"]
        dev4 = max(dev4, abs(c4 / f["cost_plus"] - 1), abs(c4 / f["cost"] - 1))
        dev5 = max(dev5, abs(c4 / c4_closed(D, J, gbar)[0] - 1))
        py = np.fft.fft(f["G"] @ f["y"]) / np.sqrt(n)
        nonres = max(nonres, np.abs(py[0::2]).max() / np.abs(py).max())
check("N = 4..20, random locus profiles and couplings: c4 by its definition equals the population cost "
      "(4 gbar / |w|^2) min_z |diag(R0^dag K1 (y + z))|^2 over Z+ and over Z' within 1e-8 (all three through "
      "one pseudo-inverse of K0)", dev4 < 1e-8, f"worst relative {dev4:.2g}")
check("read: c4 equals the closed form (2 / (N J^2 gbar)) |tau|^2 / sum_C |O_C|^2 / omega_C within 1e-8, and "
      "y's population response vanishes in the nonresonant sectors (below 1e-12 of the resonant)",
      dev5 < 1e-8 and nonres < 1e-12, f"worst relative {dev5:.2g}, nonresonant {nonres:.2g}")

# the populations at the resonant q: y, and z in Z+ through its class sums; mu_C real on the locus
dev6, imu = 0.0, 0.0
for n in range(6, 17, 2):
    J, gbar = 1.3, 0.7
    D = locus_profile(n, rnd11)
    dh = np.fft.fft(D) / n
    f = fourth(n, J, gbar, D)
    z = f["Zplus"] @ rnd11.standard_normal(f["Zplus"].shape[1])
    Bz = (f["B"] @ z).reshape(n, n)
    s = np.real(np.einsum("lk,lm,mk->k", f["Ph"].conj(), Bz, f["Ph"]))      # z = sum_k s_k P_k
    py = np.fft.fft(f["G"] @ f["y"]) / np.sqrt(n)
    pz = np.fft.fft(f["G"] @ z) / np.sqrt(n)
    for j in range(1, n, 2):
        q = 2 * np.pi * j / n
        cl = next(c for c in f["cls"] if (j - 1) // 2 in c or ((2 * n - j) % (2 * n) - 1) // 2 in c)
        Y = s[cl].sum() / (len(cl) // 2)
        m = mu(dh, n, j, J)
        want_y = np.sqrt(n) / (2 * gbar) * dh[j] * m
        want_z = -np.sqrt(n) / (2 * gbar) * dh[j] * Y
        dev6 = max(dev6, abs(py[j] - want_y) / abs(py).max(), abs(pz[j] - want_z) / abs(pz).max())
        imu = max(imu, abs(m.imag) / abs(m))
check("N = 6..16, locus profiles: at every resonant q the population that R0^dag K1 y holds is (N^(1/2) / (2 gbar)) "
      "dhat(q) mu_C and that of R0^dag K1 z, z random in Z+, is -(N^(1/2) / (2 gbar)) dhat(q) Y_C, within 1e-10 of the largest, and mu_C is real "
      "(below 1e-10)", dev6 < 1e-10 and imu < 1e-10, f"worst {dev6:.2g}, Im mu {imu:.2g}")

# the single-pair profile e_0 - e_(N-1): by hand against both routes, both signs of gbar
dev7, sgn7 = 0.0, True
for n in range(4, 25, 2):
    D = np.zeros(n); D[0], D[n - 1] = 1.0, -1.0
    for J, gbar in ((1.3, 0.7), (-0.9, -0.7)):
        f = fourth(n, J, gbar, D)
        hand = c4_pair(n, J, gbar)
        sgn7 &= np.sign(f["c4"].real) == np.sign(gbar)
        dev7 = max(dev7, abs(f["c4"] / hand - 1), abs(c4_closed(D, J, gbar)[0] / hand - 1))
check("the single-pair profile e_0 - e_(N-1), N = 4..24, (J, gbar) = (1.3, 0.7) and (-0.9, -0.7): c4 by its definition and by the "
      "closed form equals the formula by hand within 1e-8, with the sign of gbar", dev7 < 1e-8 and sgn7,
      f"worst relative {dev7:.2g}")

# the single-pair profile through the spectrum itself: no pseudo-inverse. The scalar on w is even in the
# offset e, so lambda / (e^4 c4) = 1 + a e^2 + b e^4 + ...: gated as that law, (lambda / (e^4 c4) - 1) / e^2
# agreeing across e = 4e-2, 3e-2, 2e-2. The error model: the spread is b e^2, measured below 0.3 per cent
# of a at N = 4..12, and the eigensolver's rounding, eps |M| / (c4 e^6) at e = 2e-2, below 4e-4 of a; the
# threshold, 1 per cent, sits above both and fails on an e or e^3 next order
dev9 = 0.0
for n in range(4, 13, 2):
    D = np.zeros(n); D[0], D[n - 1] = 1.0, -1.0
    J, gbar = 1.3, 0.7
    B = sym_basis(n)
    a2 = []
    for e in (4e-2, 3e-2, 2e-2):
        ev = np.linalg.eigvals(B.T @ corner_f(n, -1, J, gbar, e * D) @ B)
        lam = ev[np.argmin(abs(ev))]
        a2.append((lam / e ** 4 / c4_pair(n, J, gbar) - 1).real / e ** 2)
    dev9 = max(dev9, (max(a2) - min(a2)) / max(abs(np.mean(a2)), 1e-2))
check("the single-pair profile, N = 4..12: the smallest eigenvalue of the corner's U = +1 block is c4 e^4 "
      "(1 + a e^2 + ...): (lambda / (e^4 c4) - 1) / e^2 agrees across e = 4e-2, 3e-2, 2e-2 within 1 per cent",
      dev9 < 1e-2, f"worst spread {dev9:.2g}")

# c4 vanishes on the hypersurface tau = 0 of the locus


def min_dhat(D):
    n = len(D); dh = np.fft.fft(D) / n
    return min(abs(dh[j]) for j in range(1, n, 2)) / np.abs(dh).max()


ok8, rat8 = True, []
for n in (6, 8, 10):
    found = False
    for tries in range(400):
        Da, Db = locus_profile(n, rnd11), locus_profile(n, rnd11)

        def T(s):
            with np.errstate(divide="ignore", invalid="ignore"):
                return (c4_closed((1 - s) * Da + s * Db, 1.0, 1.0)[1] / 1j).real
        if not T(0.0) * T(1.0) < 0:
            continue
        lo, hi = 0.0, 1.0
        for _ in range(60):
            mid = (lo + hi) / 2
            tm = T(mid)
            if not np.isfinite(tm):
                break
            if T(lo) * tm <= 0:
                hi = mid
            else:
                lo = mid
        Dz = (1 - lo) * Da + lo * Db
        if min_dhat(Dz) < 1e-2:           # a pole of tau, where some dhat(q_C) vanishes, not a zero
            continue
        ref = max(fourth(n, 1.0, 1.0, Da)["c4"].real, fourth(n, 1.0, 1.0, Db)["c4"].real)
        rat8.append(abs(fourth(n, 1.0, 1.0, Dz)["c4"].real) / ref)
        found = True
        break
    ok8 &= found
check("read: N = 6, 8, 10, on a segment of locus profiles where tau / i changes sign with every resonant dhat "
      "bounded away from zero, c4 by its definition drops below 1e-6 of its value at the ends",
      ok8 and max(rat8, default=1.0) < 1e-6, f"ratios {[f'{r:.1g}' for r in rat8]}")

print()
print("ring frozen band gate:", "ALL GREEN" if not FAILS else f"{len(FAILS)} FAILED: {FAILS}")
sys.exit(1 if FAILS else 0)
