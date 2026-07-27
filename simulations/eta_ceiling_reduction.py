"""Gate for F141, F142, F143 and the eta-ceiling reduction.

Self-contained: builds every object from scratch and imports nothing from
the framework, so agreement with the other gates is evidence and not a
shared bug.

Run:
    python simulations/eta_ceiling_reduction.py            # about a minute
    python simulations/eta_ceiling_reduction.py --deep     # + the larger rungs

What it checks, block by block:

  V1  the two SU(2)s.  Phi (eta-pairing, unstaggered) commutes with the
      Liouvillian for any number-conserving quadratic H and ANY rate
      profile; S+ (spin, staggered) commutes iff the single-excitation
      matrix h is odd under Sigma = diag((-1)^l).  Sigma-oddness is
      SUFFICIENT for the off-diagonal floor through that ladder and NOT
      necessary for the off-diagonal band to be occupied: block (c) shows
      the star filling it without being Sigma-odd.
  V2  the band as ONE SO(4) irrep: S+ carries the frozen space at (1,1)
      ONTO the frozen space at (2,0), principal cosine 1.
  V3  F141: the disagreement count K is eta-invariant, [K, Phi] = 0
      exactly, so K acts only on the sl(2) multiplicity spaces and the
      spectrum on rung p is the union of the lowest-weight spectra at
      rungs l <= min(p, N-p).  The upper limit is NOT p: a multiplet
      seeded at l spans rungs l .. N-l, so taking the union to p
      overcounts on the upper half of the ladder.
  V4  the overlap closed form 2M*T_abcd, as an integer, against the
      analytic sine modes.
  V5  F143: the reduced seed operator G = (1/M)(J + (I+R)/2), its
      spectrum, and the gap 1/(N+1).
  V6  the large-J reduction: order 1 bounds the multiplicity from above,
      every EVEN order vanishes by transpose parity, and order 3 closes
      the one N where order 1 is loose.
  V7  the ceiling as an exact mod-p certificate on the reduced object.
  V8  the measured law min K = l(N-l)/(N+1) on LW_l intersected with
      V0 = ker(ad_h).  The intersection is load-bearing: on the FULL
      lowest-weight space the minimum is 0, since K is diagonal there
      with integer entries and preserves ker Psi.
  V9  the maximizer in closed form: one vector per pair of chiral pairs,
      plus the exact rational nullity, which needs Q and not a prime
      because it asserts that vectors EXIST.

Conventions.  Open XY chain of N sites, H = (J/2) sum_b (X X + Y Y), so the
single-excitation matrix h has hopping J on neighbours.  Site-resolved
Z-dephasing with rates gamma_l, gbar their mean.  Cells |A><B| with A, B
site SETS; blocks (p, q) by joint popcount.  The frozen root is -4*gbar.
M := N + 1 throughout, and p is a popcount, never a prime; primes are
written q_p.
"""

import argparse
import itertools
import sys
from fractions import Fraction as Fr
from math import comb, gcd

import numpy as np

TOL = 1e-10
# No square root of -1 is needed here (the matrices are integer, not
# Gaussian-integer), so these are just two large primes.
PRIMES = (998244353, 1004535809)
FAILURES = []


def check(label, ok, detail=""):
    if not ok:
        FAILURES.append(label)
    print(f"  [{'OK  ' if ok else 'FAIL'}] {label}{('  ' + detail) if detail else ''}")
    return ok


# --------------------------------------------------------------------------
# basis, Hamiltonians, blocks
# --------------------------------------------------------------------------

def subsets(N, p):
    return [sum(1 << i for i in c) for c in itertools.combinations(range(N), p)]


def bits(x, N):
    return [i for i in range(N) if (x >> i) & 1]


def h_chain(N, J=1.0, diag=None, extra=()):
    h = np.zeros((N, N))
    for l in range(N - 1):
        h[l, l + 1] = h[l + 1, l] = J
    for (a, b, v) in extra:
        h[a, b] = h[b, a] = v
    if diag is not None:
        h += np.diag(np.asarray(diag, dtype=float))
    return h


def h_sector(N, p, h):
    """H restricted to the p-excitation sector, in the cell basis."""
    basis = subsets(N, p)
    idx = {b: i for i, b in enumerate(basis)}
    H = np.zeros((len(basis), len(basis)))
    for A in basis:
        for x in range(N):
            for y in range(N):
                if x != y and h[x, y] and (A >> x) & 1 and not (A >> y) & 1:
                    H[idx[A ^ (1 << x) ^ (1 << y)], idx[A]] += h[x, y]
        H[idx[A], idx[A]] += sum(h[l, l] for l in bits(A, N))
    return H, basis


def block_L(N, p, q, h, gammas):
    Hp, bp = h_sector(N, p, h)
    Hq, bq = h_sector(N, q, h)
    rate = np.array([[2.0 * sum(gammas[l] for l in bits(A ^ B, N)) for B in bq]
                     for A in bp])
    L = (-1j * (np.kron(Hp, np.eye(len(bq))) - np.kron(np.eye(len(bp)), Hq.T))
         - np.diag(rate.reshape(-1)).astype(complex))
    return L, bp, bq


def frozen_space(N, p, q, h, gammas):
    L, bp, bq = block_L(N, p, q, h, gammas)
    lam = -4.0 * np.mean(gammas)
    u, s, vh = np.linalg.svd(L - lam * np.eye(L.shape[0]))
    k = int(np.sum(s < TOL * max(s[0], 1.0)))
    return (vh[len(s) - k:].conj().T if k else np.zeros((L.shape[0], 0))), k


def jw_sign(A, l, N):
    return -1 if bin(A & ((1 << l) - 1)).count("1") % 2 else 1


def phi_matrix(N, p, q):
    """Phi(rho) = sum_l d+_l rho d_l : (p,q) -> (p+1,q+1)."""
    bp, bq = subsets(N, p), subsets(N, q)
    bP, bQ = subsets(N, p + 1), subsets(N, q + 1)
    iP = {b: i for i, b in enumerate(bP)}
    iQ = {b: i for i, b in enumerate(bQ)}
    M = np.zeros((len(bP) * len(bQ), len(bp) * len(bq)))
    for i, A in enumerate(bp):
        for j, B in enumerate(bq):
            for l in range(N):
                if (A >> l) & 1 or (B >> l) & 1:
                    continue
                M[iP[A | (1 << l)] * len(bQ) + iQ[B | (1 << l)], i * len(bq) + j] \
                    += jw_sign(A, l, N) * jw_sign(B, l, N)
    return M


def splus_matrix(N, p, q):
    """S+ = sum_l (-1)^l d+_l rho d+_l : (p,q) -> (p+1,q-1)."""
    bp, bq = subsets(N, p), subsets(N, q)
    bP, bQ = subsets(N, p + 1), subsets(N, q - 1)
    iP = {b: i for i, b in enumerate(bP)}
    iQ = {b: i for i, b in enumerate(bQ)}
    M = np.zeros((len(bP) * len(bQ), len(bp) * len(bq)))
    for i, A in enumerate(bp):
        for j, B in enumerate(bq):
            for l in range(N):
                if (A >> l) & 1 or not (B >> l) & 1:
                    continue
                sg = jw_sign(A, l, N) * jw_sign(B ^ (1 << l), l, N) * (-1) ** l
                M[iP[A | (1 << l)] * len(bQ) + iQ[B ^ (1 << l)], i * len(bq) + j] += sg
    return M


def rel(x, ref):
    return float(np.max(np.abs(x)) / max(np.max(np.abs(ref)), 1.0))


# --------------------------------------------------------------------------
# the mode side: exact energy classes, the integer overlap, the reduction
# --------------------------------------------------------------------------

def modes(N, J=1.0):
    """Analytic sine modes: eps_a = 2J cos(a pi / M), v_a(l) ~ sin(a l pi / M)."""
    M = N + 1
    a = np.arange(1, N + 1)
    l = np.arange(1, N + 1)
    V = np.sqrt(2.0 / M) * np.sin(np.pi * np.outer(a, l) / M)   # V[a-1, l-1]
    eps = 2.0 * J * np.cos(np.pi * a / M)
    return eps, V


def two_m_T(a, b, c, d, twoM):
    """2M * sum_l v_a v_b v_c v_d, exactly, as an integer.

    The condition is DIVISIBILITY by 2M, not equality to zero: a+b+c-d runs
    up to 3M-4 and can hit 2M.  Zero when a+b+c+d is odd, because the
    summand is odd under l -> M-l.
    """
    if (a + b + c + d) % 2:
        return 0
    s = 0
    for sb in (1, -1):
        for sc in (1, -1):
            for sd in (1, -1):
                if (a + sb * b + sc * c + sd * d) % twoM == 0:
                    s += sb * sc * sd
    return s


def cyclotomic(n):
    """Integer coefficients of Phi_n, highest degree first."""
    num = [1] + [0] * (n - 1) + [-1]              # x^n - 1
    for d in range(1, n):
        if n % d == 0:
            num = poly_div(num, cyclotomic(d))
    return num


def poly_div(num, den):
    num = list(num)
    out = []
    while len(num) >= len(den):
        c = num[0] // den[0]
        out.append(c)
        for j in range(len(den)):
            num[j] -= c * den[j]
        num.pop(0)
    return out


def poly_rem(num, phi):
    num = list(num)
    dp = len(phi) - 1
    for i in range(len(num) - dp):
        c = num[i]
        if c:
            for j in range(len(phi)):
                num[i + j] -= c * phi[j]
    return tuple(num[len(num) - dp:])


def energy_classes(N, l):
    """l-subsets of modes grouped by EXACT Slater energy (cyclotomic, no tolerance)."""
    twoM = 2 * (N + 1)
    phi = cyclotomic(twoM)
    buckets = {}
    for S in itertools.combinations(range(1, N + 1), l):
        v = [0] * twoM
        for k in S:
            v[k] += 1
            v[twoM - k] += 1
        buckets.setdefault(poly_rem(v[::-1], phi), []).append(S)
    return list(buckets.values())


def v0_pairs(N, l):
    pairs = [(A, B) for g in energy_classes(N, l) for A in g for B in g]
    return pairs, {p: i for i, p in enumerate(pairs)}


def sign_remove(S, m):
    return -1 if sum(1 for j in S if j < m) % 2 else 1


def integer_pieces(N, l):
    """Z = 2M*D and Psi on the V0 basis, both exact integer matrices."""
    twoM = 2 * (N + 1)
    pairs, ipair = v0_pairs(N, l)
    n = len(pairs)
    mv = {}
    for S in set(x for pr in pairs for x in pr):
        out = []
        for b in S:
            rest = tuple(x for x in S if x != b)
            sb = sign_remove(S, b)
            for a in range(1, N + 1):
                if a in rest:
                    continue
                newS = tuple(sorted(rest + (a,)))
                out.append((a, b, newS, sb * sign_remove(newS, a)))
        mv[S] = out
    Z = np.zeros((n, n), dtype=np.int64)
    for col, (A, B) in enumerate(pairs):
        for (a, b, A2, sA) in mv[A]:
            for (c, d, B2, sB) in mv[B]:
                r = ipair.get((A2, B2))
                if r is None:
                    continue
                t = two_m_T(a, b, c, d, twoM)
                if t:
                    Z[r, col] += sA * sB * t
    lo_pairs, ilo = v0_pairs(N, l - 1)
    P = np.zeros((len(lo_pairs), n), dtype=np.int64)
    for col, (A, B) in enumerate(pairs):
        for m in set(A) & set(B):
            r = ilo.get((tuple(x for x in A if x != m), tuple(x for x in B if x != m)))
            if r is not None:
                P[r, col] += sign_remove(A, m) * sign_remove(B, m)
    return Z, P, pairs


def bareiss_rank(rows):
    """Exact fraction-free rank over Q of an integer matrix. Needed where mod q
    is the WRONG direction: a mod-q rank can only drop, so it can certify that
    no vector exists but never that one does."""
    m = [list(r) for r in rows]
    R, C = len(m), len(m[0])
    rank, prev = 0, 1
    for col in range(C):
        piv = next((r for r in range(rank, R) if m[r][col]), None)
        if piv is None:
            continue
        m[rank], m[piv] = m[piv], m[rank]
        for r in range(rank + 1, R):
            for cc in range(col + 1, C):
                m[r][cc] = (m[rank][col] * m[r][cc] - m[r][col] * m[rank][cc]) // prev
            m[r][col] = 0
        prev = m[rank][col]
        rank += 1
        if rank == R:
            break
    return rank


def rank_mod_p(A, q_p):
    a = np.mod(A, q_p).astype(np.int64)
    rows, cols = a.shape
    r = 0
    for c in range(cols):
        if r >= rows:
            break
        nz = np.nonzero(a[r:, c])[0]
        if nz.size == 0:
            continue
        i = r + nz[0]
        if i != r:
            a[[r, i]] = a[[i, r]]
        a[r] = (a[r] * pow(int(a[r, c]), q_p - 2, q_p)) % q_p
        below = np.nonzero(a[r + 1:, c])[0]
        if below.size:
            idx = r + 1 + below
            a[idx] = (a[idx] - a[idx, c][:, None] * a[r][None, :]) % q_p
        r += 1
    return r


def reduced_orders(N, p, gam4=1.0):
    """(dim V0, order-1 count, order-2 count, order-3 count) of the frozen branch."""
    basis = subsets(N, p)
    H, _ = h_sector(N, p, h_chain(N, 1.0))
    E, U = np.linalg.eigh(H)
    d = len(basis)
    Dm = np.array([[gam4 * (1.0 - p + bin(A & B).count("1")) for B in basis]
                   for A in basis])
    deg = [(a, b) for a in range(d) for b in range(d) if abs(E[a] - E[b]) < 1e-9]
    Bm = np.empty((d * d, d * d))
    for a in range(d):
        for b in range(d):
            Bm[a * d + b] = (U.T @ (U[:, a][:, None] * Dm * U[:, b][None, :]) @ U).reshape(-1)
    rows = np.array([Bm[a * d + b] for (a, b) in deg])
    B1 = rows[:, [c * d + e for (c, e) in deg]]
    w1, Vv = np.linalg.eigh((B1 + B1.T) / 2)
    ker1 = Vv[:, np.abs(w1) < 1e-8]
    m1 = ker1.shape[1]
    with np.errstate(divide="ignore", invalid="ignore"):
        Sm = 1.0j / (E[:, None] - E[None, :])
    Sm[~np.isfinite(Sm)] = 0.0
    for (a, b) in deg:
        Sm[a, b] = 0.0
    flatS = Sm.reshape(-1)
    if m1 == 0:
        return len(deg), 0, 0, 0, 0.0
    M2 = (rows * flatS[None, :]) @ rows.T
    B2 = -(ker1.T @ M2 @ ker1)
    m2 = int(np.sum(np.abs(np.linalg.eigvals(B2)) < 1e-8))
    # counting zero EIGENVALUES of B2 is weaker than B2 being the zero matrix
    # (a nonzero nilpotent would pass), and Prop 5.2 claims the latter, so
    # measure the norm and let the check assert what the proof says.
    b2norm = float(np.max(np.abs(B2))) if B2.size else 0.0
    left = rows * flatS[None, :]
    m3 = int(np.sum(np.abs(np.linalg.eigvals(ker1.T @ (left @ Bm @ left.T) @ ker1)) < 1e-8))
    return len(deg), m1, m2, m3, b2norm


# --------------------------------------------------------------------------
# the checks
# --------------------------------------------------------------------------

def _star_free_fermion_nullity(N):
    """Zero modes of the two-excitation hopping on the star, priced by ORDER
    (Jordan-Wigner signs) instead of by occupation."""
    A_ = [tuple(c) for c in itertools.combinations(range(N), 2)]
    ix = {x: i for i, x in enumerate(A_)}
    M = np.zeros((len(A_), len(A_)))
    for A in A_:
        st = set(A)
        for x in A:
            for y in range(N):
                if y in st or not (x == 0 or y == 0):
                    continue
                lo, hi = min(x, y), max(x, y)
                sg = -1.0 if sum(1 for z in st if lo < z < hi and z != x) % 2 else 1.0
                M[ix[tuple(sorted((st - {x}) | {y}))], ix[A]] += sg
    return int(np.sum(np.abs(np.linalg.eigvalsh(M)) < 1e-9))


def v1_two_su2():
    print("\nV1  the two SU(2)s, and the seed gate as Sigma-oddness")
    J, gam = 0.7, 0.13
    profiles = [("uniform", lambda N: np.full(N, gam)),
                ("R90 locus", lambda N: np.array([gam + 0.4 * gam * (1 if 2 * i < N - 1
                                                  else (-1 if 2 * i > N - 1 else 0))
                                                  for i in range(N)]))]
    for (pname, prof) in profiles:
        for N in (4, 5, 6):
            g = prof(N)
            h = h_chain(N, J)
            for (p, q) in [(1, 1), (2, 2)]:
                Lo, _, _ = block_L(N, p, q, h, g)
                Lu, _, _ = block_L(N, p + 1, q + 1, h, g)
                P = phi_matrix(N, p, q)
                check(f"(a) Phi commutes, {pname}, N={N}, ({p},{q})",
                      rel(Lu @ P - P @ Lo, P) < 1e-10)
    def star(N):
        """Centre 0 joined to every other site.  Bipartite, so its spectrum is
        symmetric about zero, but its two colours are centre-versus-leaves and
        NOT even-versus-odd, so it is not Sigma-odd.  This is the case that
        SEPARATES the two conditions; without it the block would only show
        them agreeing."""
        h = np.zeros((N, N))
        for i in range(1, N):
            h[0, i] = h[i, 0] = J
        return h

    variants = [
        ("chain", lambda N: h_chain(N, J), True),
        ("chain + alternating bonds", lambda N: h_chain(N, J, extra=[(l, l + 1, 0.55 * J)
                                                                     for l in range(1, N - 1, 2)]), True),
        ("chain + on-site defect", lambda N: h_chain(N, J, diag=[0.9] + [0.0] * (N - 1)), False),
        ("chain + range-2 bonds", lambda N: h_chain(N, J, extra=[(l, l + 2, 0.6 * J)
                                                                 for l in range(N - 2)]), False),
        ("star", star, False),
    ]
    seen = []
    for (name, hf, want) in variants:
        for N in (4, 5, 6):
            h = hf(N)
            sg = np.diag([(-1.0) ** l for l in range(N)])
            stag = bool(np.max(np.abs(sg @ h @ sg + h)) < 1e-12)
            e = np.sort(np.linalg.eigvalsh(h))
            sym = bool(np.max(np.abs(e + e[::-1])) < 1e-10)
            g = np.full(N, gam)
            Lo, _, _ = block_L(N, 1, 1, h, g)
            Lu, _, _ = block_L(N, 2, 0, h, g)
            S = splus_matrix(N, 1, 1)
            comm = rel(Lu @ S - S @ Lo, S) < 1e-10
            seen.append((sym, stag))
            check(f"(b) {name}, N={N}: sym-spectrum={sym}, Sigma-odd={stag}, "
                  f"S+ commutes={comm}", stag == comm == want)
    check("(b) the two conditions are SEPARATED, not merely correlated: some "
          "Hamiltonian has a symmetric spectrum and is not Sigma-odd",
          any(sym and not stag for (sym, stag) in seen),
          "the star is that case, so Sigma-oddness is strictly stronger")

    # (c) Sigma-oddness is SUFFICIENT for the off-diagonal floor via the S+ ladder,
    # and NOT NECESSARY for the band to exist at all.  The star fills the (2,0)
    # block while failing Sigma-oddness, so it is a counterexample to necessity
    # and NOT a witness for replacing the measured band-existence gate.
    for N in (5, 6, 7):
        h = star(N)
        g = np.full(N, gam)
        _, k20 = frozen_space(N, 2, 0, h, g)
        sg = np.diag([(-1.0) ** l for l in range(N)])
        stag = bool(np.max(np.abs(sg @ h @ sg + h)) < 1e-12)
        # The count is the SPIN sector's, which is what L is built from.  The
        # star's bonds are not adjacent in any ordering, so the JW strings do
        # not cancel and it is not a free-fermion model; pricing by order gives
        # the Slater count 1 + C(N-2,2) instead, which is what a quick
        # recomputation produces.  Pin both so the difference is on the record.
        ferm = _star_free_fermion_nullity(N)
        check(f"(c) star, N={N}: spin-sector count {k20}, free-fermion count "
              f"{ferm} = 1+C(N-2,2) = {1 + comb(N - 2, 2)}",
              ferm == 1 + comb(N - 2, 2) and k20 != ferm,
              "the two pricings differ here, so the model is not free-fermion")
        check(f"(c) star, N={N}: Sigma-odd={stag} yet the (2,0) block carries "
              f"{k20} frozen modes", (not stag) and k20 > 0,
              "so Sigma-oddness is sufficient for the LADDER, not necessary for "
              "the BAND; it does not replace the measured band-existence gate")

    # (d) Lemma 2.1 claims BOTH ladders are blind to the watching at ANY profile.
    # Phi is exercised on a non-uniform profile above; do the same for S+.
    for N in (4, 5, 6):
        gg = np.array([gam * (1 + 0.3 * i) for i in range(N)])
        h = h_chain(N, J)
        Lo, _, _ = block_L(N, 1, 1, h, gg)
        Lu, _, _ = block_L(N, 2, 0, h, gg)
        S = splus_matrix(N, 1, 1)
        check(f"(d) S+ commutes at a GRADED profile too, N={N}",
              rel(Lu @ S - S @ Lo, S) < 1e-10)


def v2_so4_triplet():
    print("\nV2  the band as one SO(4) irrep: S+ carries (1,1) onto (2,0)")
    J, gam = 0.7, 0.13
    for N in (4, 5, 6):
        g = np.full(N, gam)
        h = h_chain(N, J)
        F11, k11 = frozen_space(N, 1, 1, h, g)
        F20, k20 = frozen_space(N, 2, 0, h, g)
        worst = 0.0
        if k11 and k20:
            Q1, _ = np.linalg.qr(splus_matrix(N, 1, 1) @ F11)
            Q2, _ = np.linalg.qr(F20)
            worst = float(np.linalg.svd(Q1.conj().T @ Q2, compute_uv=False).min())
        check(f"N={N}: dim(1,1)={k11} dim(2,0)={k20}=floor(N/2)={N // 2}",
              k11 == k20 == N // 2 and worst > 1 - 1e-8,
              f"worst principal cosine {worst:.9f}")


def v2b_splus_injective():
    """The off-diagonal FLOOR needs S+ to be injective on the corner seeds, and
    commutation alone does not give that.  Three exact facts do, and they hold
    at every N rather than at the N a census reaches."""
    print("\nV2(b) S+ is INJECTIVE on the seed space, derived rather than measured")
    J, gam = 0.7, 0.13
    for N in range(4, 11):
        h = h_chain(N, J)
        eps, V = np.linalg.eigh(h)
        Sg = np.diag([(-1.0) ** l for l in range(N)])
        S = splus_matrix(N, 1, 1)
        pairs = list(itertools.combinations(range(N), 2))
        # (i) for SYMMETRIC X the map is exactly the commutator, read in Lambda^2
        rng = np.random.default_rng(11)
        worst = 0.0
        for _ in range(8):
            Mr = rng.normal(size=(N, N))
            X = (Mr + Mr.T) / 2
            C = Sg @ X - X @ Sg
            worst = max(worst, float(np.max(np.abs(
                S @ X.reshape(-1) - np.array([C[a, b] for (a, b) in pairs])))))
        # (ii) the seeds are Sigma-ODD, so the commutator is 2*Sigma*X
        seeds = [np.outer(V[:, a], V[:, a]) - np.outer(V[:, b], V[:, b])
                 for a in range(N) for b in range(a + 1, N)
                 if abs(eps[a] + eps[b]) < 1e-10]
        odd = max(float(np.max(np.abs(Sg @ X @ Sg + X))) for X in seeds)
        # (iii) hence injective on the span
        img = np.array([S @ (X / np.linalg.norm(X)).reshape(-1) for X in seeds]).T
        smin = float(np.linalg.svd(img, compute_uv=False).min())
        check(f"N={N}: S+ = [Sigma, .] on symmetric X ({worst:.1e}), seeds "
              f"Sigma-odd ({odd:.1e}), smallest singular value of the image "
              f"{smin:.6f}", worst < 1e-12 and odd < 1e-12 and smin > 1e-6,
              f"{len(seeds)} seeds = floor(N/2) = {N // 2}"
              if len(seeds) == N // 2 else "SEED COUNT WRONG")


def v3_eta_invariance():
    print("\nV3  F141: the disagreement count is eta-invariant, [K, Phi] = 0")
    J, gam = 0.7, 0.13
    for N in (4, 5, 6):
        for prof in ([gam] * N, [gam * (1 + 0.3 * i) for i in range(N)]):
            g = np.array(prof)
            h = h_chain(N, J)
            for p in (1, 2):
                if p + 1 >= N:
                    continue
                bp, bu = subsets(N, p), subsets(N, p + 1)
                P = phi_matrix(N, p, p)
                ko = np.array([2.0 * sum(g[l] for l in bits(A ^ B, N))
                               for A in bp for B in bp])
                ku = np.array([2.0 * sum(g[l] for l in bits(A ^ B, N))
                               for A in bu for B in bu])
                Hp, _ = h_sector(N, p, h)
                Hu, _ = h_sector(N, p + 1, h)
                ado = np.kron(Hp, np.eye(len(bp))) - np.kron(np.eye(len(bp)), Hp.T)
                adu = np.kron(Hu, np.eye(len(bu))) - np.kron(np.eye(len(bu)), Hu.T)
                nm = "uniform" if len(set(prof)) == 1 else "graded"
                check(f"(a) N={N} p={p} {nm}: [dissipator, Phi] = 0 (the rate "
                      f"2*sum gamma over A xor B, of which the integer count K "
                      f"is the uniform case) and [ad_h, Phi] = 0",
                      rel(ku[:, None] * P - P * ko[None, :], P) < 1e-12
                      and rel(adu @ P - P @ ado, P) < 1e-12)

    # (b) the consequence: on V0 the spectrum of K at rung p is the union of the
    # lowest-weight spectra at rungs l <= min(p, N-p).  The upper limit is NOT p:
    # a multiplet seeded at l spans rungs l .. N-l, so it does not reach rung p
    # once p > N-l, and the union over l <= p would overcount on the upper half.
    for N in (4, 5, 6, 7):
        for p in range(1, N):
            pairs_p, _ = v0_pairs(N, p)
            if len(pairs_p) > 400:
                continue
            Zp, _, _ = integer_pieces(N, p)
            Kp = p * np.eye(Zp.shape[0]) - Zp / (2.0 * (N + 1))
            full = np.sort(np.round(np.linalg.eigvalsh((Kp + Kp.T) / 2), 8))
            acc = []
            for l in range(0, min(p, N - p) + 1):
                Zl, Pl, _ = integer_pieces(N, l) if l >= 1 else (np.zeros((1, 1)), None, None)
                if l == 0:
                    acc.append(0.0)
                    continue
                nl = Zl.shape[0]
                Kl = l * np.eye(nl) - Zl / (2.0 * (N + 1))
                _, s, vh = np.linalg.svd(Pl.astype(float), full_matrices=True)
                nk = nl - int(np.sum(s > 1e-9))
                if nk == 0:
                    continue
                Lb = vh[nl - nk:].T
                Rm = Lb.T @ ((Kl + Kl.T) / 2) @ Lb
                acc += list(np.round(np.linalg.eigvalsh((Rm + Rm.T) / 2), 8))
            acc = np.sort(np.array(acc))
            ok = len(acc) == len(full) and np.max(np.abs(acc - full)) < 1e-6
            check(f"(b) N={N} p={p}: spec on the rung = union over l <= "
                  f"min(p, N-p) = {min(p, N - p)}", ok,
                  f"dims {len(acc)} vs {len(full)}")


def v3b_full_space_is_vacuous():
    """The docs justify intersecting with V0 by saying the criterion is vacuous
    on the FULL lowest-weight space.  That is a load-bearing claim and was not
    computed anywhere, so compute it."""
    print("\nV3(e) why the intersection with V0 is load-bearing: on the FULL "
          "lowest-weight space min K = 0 and the eigenvalue 1 is present")
    for N in (4, 6, 7, 8):
        p = 2
        b = subsets(N, p)
        d = len(b)
        kv = np.array([float(bin(A ^ B).count("1")) / 2 for A in b for B in b])
        P = phi_matrix(N, p - 1, p - 1).T
        _, s, vh = np.linalg.svd(P, full_matrices=True)
        nk = d * d - int(np.sum(s > 1e-9))
        L = vh[d * d - nk:].T
        K = L.T @ (kv[:, None] * L)
        w = np.linalg.eigvalsh((K + K.T) / 2)
        has1 = bool(np.any(np.abs(w - 1) < 1e-9))
        check(f"N={N} l={p}: dim ker(Psi)={nk}, min K = {w.min():.3f}, "
              f"eigenvalue 1 present = {has1}",
              abs(w.min()) < 1e-9 and has1)


def v4_overlap_closed_form(deep):
    print("\nV4  the integer overlap 2M*T_abcd against the analytic modes")
    for N in range(3, 11 if not deep else 14):
        M = N + 1
        _, V = modes(N)
        e = 0.0
        for a in range(1, N + 1):
            for b in range(1, N + 1):
                for c in range(1, N + 1):
                    for d in range(1, N + 1):
                        num = 2 * M * float(np.sum(V[a - 1] * V[b - 1] * V[c - 1] * V[d - 1]))
                        e = max(e, abs(num - two_m_T(a, b, c, d, 2 * M)))
        check(f"N={N}: all {N ** 4} quadruples", e < 1e-9, f"max error {e:.2e}")


def v5_seed_spectrum(deep):
    print("\nV5  F143: G = (1/M)(J + (I+R)/2), its spectrum, and the gap 1/(N+1)")
    for N in list(range(2, 13)) + ([20, 40] if deep else []):
        M = N + 1
        _, V = modes(N)
        W = (V ** 2).T                      # W[l-1, a-1] = v_a(l)^2
        G = W.T @ W
        Jm = np.ones((N, N))
        R = np.zeros((N, N))
        for a in range(1, N + 1):
            R[a - 1, M - a - 1] = 1.0
        closed = (Jm + (np.eye(N) + R) / 2) / M
        w = np.sort(np.linalg.eigvalsh(G))
        pred = np.sort(np.array([0.0] * (N // 2) + [1.0 / M] * ((N + 1) // 2 - 1) + [1.0]))
        check(f"N={N}: G matches the closed form, spectrum as predicted",
              np.max(np.abs(G - closed)) < 1e-12 and np.max(np.abs(w - pred)) < 1e-12,
              f"|G-closed| {np.max(np.abs(G - closed)):.1e}")


def v6_reduction(deep):
    print("\nV6  the large-J reduction: order 1 tight, even orders zero, order 3 closes N=5")
    for N in (4, 5, 6, 7):
        for p in range(1, N):
            if comb(N, p) > 40:
                continue
            dim, m1, m2, m3, b2 = reduced_orders(N, p)
            # V6 builds V0 by a float degeneracy tolerance while V7/V8 use the
            # exact cyclotomic classes.  That is the very trap this arc records,
            # so assert the two agree rather than assume it.
            exact_dim = len(v0_pairs(N, p)[0])
            fl = N // 2
            ok = (m3 == fl and m2 == m1 and m1 >= fl and b2 < 1e-9
                  and dim == exact_dim)
            note = "order 1 tight" if m1 == fl else f"order 1 loose by {m1 - fl}, order 3 closes"
            check(f"N={N} p={p}: dim V0={dim}, orders {m1}/{m2}/{m3}, floor={fl}, "
                  f"|order 2| = {b2:.1e} (the ZERO MATRIX, not merely zero "
                  f"eigenvalues); float dim V0 == exact dim V0 = {exact_dim}",
                  ok, note)


def v7_ceiling_certificate(deep):
    print("\nV7  the ceiling as an exact mod-p certificate on the reduced object")
    cap = 700 if not deep else 2200
    for N in range(4, 17 if deep else 14):
        for l in range(2, N // 2 + 1):
            pairs, _ = v0_pairs(N, l)
            if len(pairs) > cap:
                continue
            Z, P, _ = integer_pieces(N, l)
            n = Z.shape[0]
            A = np.vstack([Z - 2 * (N + 1) * (l - 1) * np.eye(n, dtype=np.int64), P])
            r = max(rank_mod_p(A, q) for q in PRIMES)
            expect_full = not (N == 5 and l == 2)
            check(f"N={N} l={l}: dim V0={n}, rank {r}, full={r == n}",
                  (r == n) == expect_full,
                  "" if expect_full else "N=5 is the known resonance, a lowest weight sits at K=1")


def v8_min_k_law(deep):
    print("\nV8  the measured law min K = l(N-l)/(N+1) on the lowest-weight spaces")
    cap = 700 if not deep else 2200
    for N in range(4, 17 if deep else 14):
        for l in range(2, N // 2 + 1):
            pairs, _ = v0_pairs(N, l)
            if len(pairs) > cap:
                continue
            Z, P, _ = integer_pieces(N, l)
            n = Z.shape[0]
            K = l * np.eye(n) - Z / (2.0 * (N + 1))
            _, s, vh = np.linalg.svd(P.astype(float), full_matrices=True)
            nk = n - int(np.sum(s > 1e-9))
            if nk == 0:
                continue
            Lb = vh[n - nk:].T
            Rm = Lb.T @ ((K + K.T) / 2) @ Lb
            got = float(np.linalg.eigvalsh((Rm + Rm.T) / 2).min())
            pred = Fr(l * (N - l), N + 1)
            check(f"N={N} l={l}: min K = {Fr(got).limit_denominator(999)}, "
                  f"predicted {pred}", abs(got - float(pred)) < 1e-9)


def v9_maximizer(deep):
    print("\nV9  the maximizer in closed form: one vector per pair of chiral pairs")
    for N in range(5, 13 if not deep else 17):
        M = N + 1
        Z, P, pairs = integer_pieces(N, 2)
        ipair = {pr: i for i, pr in enumerate(pairs)}
        n = Z.shape[0]
        chir = [(a, M - a) for a in range(1, M // 2 + 1) if a < M - a and M - a <= N]
        good = 0
        total = 0
        for (p1, s1), (q1, r1) in itertools.combinations(chir, 2):
            pp, qq, rr, ss = sorted([p1, s1, q1, r1])
            if pp + ss != M or qq + rr != M:
                continue
            v = np.zeros(n, dtype=np.int64)
            for (A, B, c) in [((pp, qq), (pp, qq), 1), ((pp, rr), (pp, rr), -1),
                              ((qq, ss), (qq, ss), -1), ((rr, ss), (rr, ss), 1),
                              ((pp, ss), (qq, rr), -2), ((qq, rr), (pp, ss), -2)]:
                v[ipair[(tuple(sorted(A)), tuple(sorted(B)))]] += c
            total += 1
            if np.array_equal(Z @ v, 12 * v) and not np.any(P @ v):
                good += 1
        check(f"N={N}: {total} chiral-pair pairs, Z v = 12 v and Psi v = 0 for {good}",
              total == comb(N // 2, 2) and good == total)

    # The closed-form family is a LOWER bound on the nullity.  The surplus at the
    # resonant N is a separate number and needs the rank over Q, not mod q: a
    # mod-q rank can only drop, so it cannot certify that extra vectors exist.
    print("      and the full rational nullity, where the closed form is not all of it")
    for (N, want) in ((10, 10), (11, 14), (12, 15), (14, 29)):
        if not deep and N > 11:
            continue
        Z, P, _ = integer_pieces(N, 2)
        n = Z.shape[0]
        A = np.vstack([Z - 12 * np.eye(n, dtype=np.int64), P]).tolist()
        nul = n - bareiss_rank(A)
        closed = comb(N // 2, 2)
        note = ("matches the closed form" if nul == closed
                else f"SURPLUS {nul - closed} over the closed form's {closed}, M={N + 1}")
        check(f"N={N}: exact rational nullity {nul}", nul == want, note)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--deep", action="store_true", help="larger N and larger rungs")
    args = ap.parse_args()
    print("eta-ceiling reduction gate (F141, F142, F143)")
    v1_two_su2()
    v2_so4_triplet()
    v2b_splus_injective()
    v3_eta_invariance()
    v3b_full_space_is_vacuous()
    v4_overlap_closed_form(args.deep)
    v5_seed_spectrum(args.deep)
    v6_reduction(args.deep)
    v7_ceiling_certificate(args.deep)
    v8_min_k_law(args.deep)
    v9_maximizer(args.deep)
    print()
    if FAILURES:
        print(f"eta-ceiling reduction gate: {len(FAILURES)} FAILURES")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    print("eta-ceiling reduction gate: ALL GREEN")


if __name__ == "__main__":
    main()
