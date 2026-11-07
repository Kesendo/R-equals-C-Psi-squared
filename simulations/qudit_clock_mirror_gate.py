"""Gate for docs/proofs/PROOF_QUDIT_CLOCK_MIRROR.md: the qudit chain under clock dephasing.

Jumps: the clock Z = diag(1, w, w^2, ...), w = exp(2 pi i / d), and its adjoint Z^dagger on every
site, both at rate gamma. Mirror: Pi(rho) = rho^T . (X^{d/2})^{(x)N}, X the cyclic shift.

What is exact, and why. The identities live in the ring Z[zeta_m], m = lcm(d, 4), and are
reduced modulo primes p = 1 mod m, where zeta_m is an integer of order m (a ring homomorphism
Z[zeta_m] -> GF(p); complex conjugation is the Galois map zeta -> zeta^-1 and is written into the
construction, conj(w^j) = w^-j, never applied to residues). Every mirror used here is a bare
permutation of the matrix units, so a residual entry is a reindexed entry of L plus an entry of
L or L^dagger plus the centre: an element alpha of Z[zeta_m] whose every Galois conjugate is at
most B = 2 max|L| + c in modulus, max|L| taken entry-wise on a majorant built with every root of
unity replaced by its modulus 1. If alpha is zero modulo degree-one primes above p1 and p2, then
p1 p2 divides its norm, which is at most B^phi(m); so where B^phi(m) < p1 p2 (checked beside every
identity, G2, G9, G10, G12) zero modulo the two primes is zero. The odd centred power sums (G5 and the
moment lines of G9, G10) have heights far above p1 p2: their vanishing modulo two primes is
consistency, not exactness (and G2, G9, G10 imply it). The obstruction at odd d (G4) runs the
other way and is exact as it stands: if the spectrum paired about any centre c, the trace would
force c = -2 Tr L / dim and every odd centred power sum Tr (L + c/2)^k would vanish; one that is
nonzero modulo a prime is nonzero, a no for every element and every centre.

  G1  one site, every d = 2..8: the rates 2(1 - cos(2 pi a / d)) of the letters X^a Z^b,
      exact in sympy: a -> a + d/2 sends r to 4 - r exactly when d is even; at odd d the
      rate 0 has d letters and the largest rate 2d letters, so no centre exists
  G2  even d = 2, 4, 6, 8: Pi L Pi^-1 = -L - 4 N gamma exactly (mod two primes), for the
      SU(d) Heisenberg (SWAP) chain, the clock-shift chain J (X X^dag + h.c.) + Delta (Z Z^dag + h.c.),
      and that chain with a shift field h (X + X^dag), at N = 1, 2 (N = 3 at d = 2, 4)
  G3  control: a clock field h (Z + Z^dag) breaks the identity (X^{d/2} turns Z into -Z)
  G4  odd d = 3, 5: an odd centred power sum is nonzero mod p, with H = 0, the SWAP chain and
      (at d = 3) the clock-shift chain: no palindrome about any centre
  G5  even d, the settings of G2 with d^(2N) <= 1296: the odd centred power sums k = 3, 5 vanish mod p (implied by
      G2; consistency through the spectrum's moments, not exactness)
  G6  continuous phase noise, the jump n = diag(0, 1, ..., d-1): a cell |i><j| decays at
      (i - j)^2 / 2 per unit rate (Fractions); the multiset is symmetric about some centre only at
      d = 2 (exact integer check, d = 2..8)
  G7  the hidden qubit, even d = 2..8: in the labels j = (d/2) b + j', Z_d = Z_b (x) diag(w^j')
      and X^{d/2} = X_b (x) 1, so X^{d/2} Z_d X^{-d/2} = -Z_d (exact mod p)
  G8  d = 4: the clock pair's dissipator is the two-qubit Pauli-string dephasing with jumps
      Z (x) 1 and Z (x) Z at the same rate (exact mod p, superoperators compared entry-wise)
  G9  the one-sided element, H only commuting with U = (X^{d/2})^{(x)N}: with a field
      h i (X - X^dag), Hermitian and not real, U rho reflects L onto -L^dag - 4N gamma exactly,
      the transpose mirror of G2 does not, and the odd centred power sums still vanish
  G10 the single jump Z without Z^dag: with Q = P^{(x)N}, P|j> = |-j>, and [Q, H] = [U, H] = 0,
      U . Ad_Q reflects L onto -L^dag - 2 N gamma exactly; a field i(X - X^dag) breaks P and the
      identity, while the odd centred power sums still vanish (G12 explains them)
  G11 the conventions the reduction depends on: at H = 0 the builders act on a matrix unit
      |j><k| as gamma sum_l (w^(j_l - k_l) + w^(k_l - j_l) - 2) (the pair) and
      gamma sum_l (w^(j_l - k_l) - 1) (the single jump), exactly; a wrong conjugation fails here
  G12 the lone jump through the transpose: with [U, H] = 0 and P H-bar P = H (a condition
      other than [Q, H] = 0, neither implying the other; met by i(X - X^dag)), M(rho) = U Q rho^T Q^-1 gives M L M^-1 = -L - 2N gamma
      exactly; control: a bond i(Z Z - h.c.) commutes with U, has P H-bar P = -H, and breaks it
  G13 the numerical-range no-go, its exact inequalities: the self-adjoint part of the dissipator
      bounds every eigenvalue, Re lambda >= -N r_max, while a palindrome about the trace-forced
      centre needs 0's partner at -2 N (mean rate); r_max < 2 mean at odd d for the clock pair
      (and its half, the lone jump), r_max = 2 mean at even d, and r < 2 r (d - 1)/d for uniform
      full-Cartan dephasing at d >= 3

Runtime: one to three minutes.
"""
import sys
from math import gcd
from fractions import Fraction
from itertools import product

import numpy as np
import sympy as sp

FAILS = []
COUNT = 0


def check(name, ok):
    global COUNT
    COUNT += 1
    print(("  ok    " if ok else "  FAIL  ") + name)
    if not ok:
        FAILS.append(name)


def lcm(a, b):
    return a * b // gcd(a, b)


def primes_1_mod(m, count=2, start=20000):
    out, p = [], start - (start % m) + 1
    while len(out) < count:
        if p > 2 and sp.isprime(p):
            out.append(p)
        p += m
    return out


def root_of_order(m, p):
    for g in range(2, p):
        z = pow(g, (p - 1) // m, p)
        if all(pow(z, m // q, p) != 1 for q in sp.primefactors(m)):
            return z
    raise ValueError("no root")


class Ring:
    """Z[zeta_m] reduced mod p: w = zeta^(m/d) of order d, i = zeta^(m/4)."""

    def __init__(self, d, p):
        m = lcm(d, 4)
        z = root_of_order(m, p)
        self.p, self.d = p, d
        self.w = pow(z, m // d, p)
        self.i = pow(z, m // 4, p)


def kron_list(mats, p):
    out = np.array([[1]], dtype=np.int64)
    for M in mats:
        out = np.kron(out, M) % p
    return out


def site(op, l, N, d, p):
    return kron_list([op if m == l else np.eye(d, dtype=np.int64) for m in range(N)], p)


def clock_ops(R):
    d, p = R.d, R.p
    X = np.roll(np.eye(d, dtype=np.int64), 1, axis=0)                   # X|j> = |j+1>
    Z = np.diag([pow(R.w, j, p) for j in range(d)]).astype(np.int64)
    Zd = np.diag([pow(R.w, (-j) % d, p) for j in range(d)]).astype(np.int64)   # Z^dagger = Z^-1
    return X, Z, Zd


def hamiltonian(kind, N, R, J=3, Delta=2, h=5):
    """Integer couplings; every H here is real symmetric as a complex matrix."""
    d, p = R.d, R.p
    X, Z, Zd = clock_ops(R)
    Xd = X.T.copy()
    D = d ** N
    H = np.zeros((D, D), dtype=np.int64)
    if kind == "zero":
        return H
    if kind == "swap":
        for l in range(N - 1):
            for idx in range(D):
                dig = list(np.unravel_index(idx, (d,) * N))
                dig[l], dig[l + 1] = dig[l + 1], dig[l]
                H[np.ravel_multi_index(dig, (d,) * N), idx] += J
        return H % p
    for l in range(N - 1):
        A = site(X, l, N, d, p) @ site(Xd, l + 1, N, d, p) % p        # real; A^dagger = A^T
        B = site(Z, l, N, d, p) @ site(Zd, l + 1, N, d, p) % p        # diagonal
        Bd = site(Zd, l, N, d, p) @ site(Z, l + 1, N, d, p) % p       # B^dagger, its conjugate
        H = (H + J * (A + A.T) + Delta * (B + Bd)) % p
    if kind in ("shift-field", "clock-field"):
        for l in range(N):
            if kind == "shift-field":
                H = (H + h * (site(X, l, N, d, p) + site(Xd, l, N, d, p))) % p
            else:
                H = (H + h * (site(Z, l, N, d, p) + site(Zd, l, N, d, p))) % p
    return H % p


def liouvillian(H, N, R, gamma):
    """Column-stacking vec: vec(A rho B) = (B^T (x) A) vec(rho). Jumps Z_l, Z_l^dag at rate gamma;
    they are unitary, so -1/2{J^dag J, rho} = -rho per jump. conj(Z) = Z^dag entrywise."""
    d, p = R.d, R.p
    X, Z, Zd = clock_ops(R)
    D = d ** N
    I = np.eye(D, dtype=np.int64)
    mi = (-R.i) % p
    L = (mi * (np.kron(I, H) - np.kron(H.T, I))) % p
    for l in range(N):
        Zl, Zdl = site(Z, l, N, d, p), site(Zd, l, N, d, p)
        # jump Z_l: conj(Z_l) (x) Z_l ; jump Z_l^dag: conj(Z_l^dag) (x) Z_l^dag = Z_l (x) Z_l^dag
        L = (L + gamma * (np.kron(Zdl, Zl) + np.kron(Zl, Zdl))) % p
        L = (L - 2 * gamma * np.eye(D * D, dtype=np.int64)) % p
    return L % p


def mirror_perm(N, d):
    """Pi(rho) = rho^T U, U = (X^{d/2})^{(x)N} a permutation |j> -> |u(j)>. On matrix units:
    |r><c| -> |c><r| U = |c><u^-1(r)|. Column-stacking index of (row, col) = col * D + row."""
    D = d ** N
    s = d // 2
    def u_inv(r):
        dig = np.unravel_index(r, (d,) * N)
        return int(np.ravel_multi_index([(x - s) % d for x in dig], (d,) * N))
    P = np.empty(D * D, dtype=np.int64)
    for r, c in product(range(D), range(D)):
        P[c * D + r] = u_inv(r) * D + c
    return P


def identity_residual_zero(L, N, R, gamma):
    p = R.p
    P = mirror_perm(N, R.d)
    M = np.empty_like(L)
    M[np.ix_(P, P)] = L                       # M = Pi L Pi^-1
    c = 4 * N * gamma
    T = (M + L + c * np.eye(len(L), dtype=np.int64)) % p
    return not T.any()


def majorant_max(kind, N, d, gamma, n_jumps, imag_field=False, J=3, Delta=2, h=5, hi=4):
    """Entry-wise majorant of |L| in every Galois embedding: each root of unity replaced by 1,
    each coupling by its modulus (the same couplings hamiltonian() and with_imag_shift_field() use).
    The anticommutator and every jump term are diagonal for the clock (|conj(Z) (x) Z| = 1)."""
    D = d ** N
    X = np.roll(np.eye(d), 1, axis=0)
    def site_f(op, l):
        out = np.array([[1.0]])
        for m in range(N):
            out = np.kron(out, op if m == l else np.eye(d))
        return out
    Ha = np.zeros((D, D))
    if kind == "swap":
        for l in range(N - 1):
            for idx in range(D):
                dig = list(np.unravel_index(idx, (d,) * N))
                dig[l], dig[l + 1] = dig[l + 1], dig[l]
                Ha[np.ravel_multi_index(dig, (d,) * N), idx] += J
    elif kind != "zero":
        for l in range(N - 1):
            A = site_f(X, l) @ site_f(X.T, l + 1)
            Ha += J * (A + A.T) + 2 * Delta * np.eye(D)
        for l in range(N):
            if kind == "shift-field":
                Ha += h * (site_f(X, l) + site_f(X.T, l))
            if kind == "clock-field":
                Ha += 2 * h * np.eye(D)
    if imag_field:
        for l in range(N):
            Ha += hi * (site_f(X, l) + site_f(X.T, l))
    I = np.eye(D)
    La = np.kron(I, Ha) + np.kron(Ha.T, I) + 2 * n_jumps * gamma * np.eye(D * D)
    return La.max()


def zero_is_zero(Lmax, c, d, primes):
    """B^phi(m) < p1 p2 with B = 2 max|L| + c: a residual zero mod both primes is zero."""
    B = 2 * Lmax + c
    phi = int(sp.totient(lcm(d, 4)))
    return B ** phi < primes[0] * primes[1]


def odd_power_sums(L, N, R, gamma, ks=(3, 5)):
    p = R.p
    n = len(L)
    trL = int(np.trace(L) % p)
    # Tr L = -dim * 2 N gamma exactly (the mean of 1 - cos over Z_d is 1); the forced centre is
    # c = 4 N gamma, and c/2 = 2 N gamma
    assert trL == (-n * 2 * N * gamma) % p, "trace is not -dim * 2 N gamma"
    A = (L + 2 * N * gamma * np.eye(n, dtype=np.int64)) % p
    out, Pw, k = {}, A.copy(), 1
    A2 = (A @ A) % p
    while k < max(ks):
        Pw = (Pw @ A2) % p
        k += 2
        if k in ks:
            out[k] = int(np.trace(Pw) % p)
    return out


gamma = 1
print("G1  one site: the clock rates and the half-shift")
for d in range(2, 9):
    a = sp.symbols('a')
    r = [sp.nsimplify(2 * (1 - sp.cos(2 * sp.pi * k / d))) for k in range(d)]
    if d % 2 == 0:
        ok = all(sp.simplify(r[(k + d // 2) % d] - (4 - r[k])) == 0 for k in range(d))
        check(f"d={d}: a -> a + d/2 sends the rate r(a) to 4 - r(a) for every column", ok)
    else:
        rmax = max(r, key=lambda x: float(x))
        n_max = sum(1 for x in r if sp.simplify(x - rmax) == 0)
        n_min = sum(1 for x in r if sp.simplify(x) == 0)
        check(f"d={d}: rate 0 has {n_min} column, the largest rate {n_max}: no centre "
              "(each column holds d letters)", n_min == 1 and n_max == 2)

print("G2/G3/G5  even d: the identity, its control, the odd moments")
settings = []
for d in (2, 4, 6, 8):
    Ns = (1, 2, 3) if d in (2, 4) else (1, 2)
    for N in Ns:
        for kind in ("swap", "clock-shift", "shift-field"):
            if kind == "swap" and N == 1:
                continue
            settings.append((d, N, kind))
def label(kind, N):
    if N == 1 and kind == "clock-shift":
        return "H = 0 (one site, the dissipator alone)"
    if N == 1 and kind == "shift-field":
        return "shift field alone (one site)"
    return kind


for d, N, kind in settings:
    ok_id, ok_mom = True, True
    pr = primes_1_mod(lcm(d, 4))
    ok_b = zero_is_zero(majorant_max(kind, N, d, gamma, 2 * N), 4 * N * gamma, d, pr)
    for p in pr:
        R = Ring(d, p)
        L = liouvillian(hamiltonian(kind, N, R), N, R, gamma)
        ok_id &= identity_residual_zero(L, N, R, gamma)
        if d ** (2 * N) <= 1296:
            ok_mom &= all(v == 0 for v in odd_power_sums(L, N, R, gamma).values())
    check(f"d={d} N={N} {label(kind, N)}: Pi L Pi^-1 = -L - 4N gamma mod two primes, and the "
          "norm bound makes that exact", ok_id and ok_b)
    if d ** (2 * N) <= 1296:
        check(f"d={d} N={N} {label(kind, N)}: Tr (L + 2N gamma)^k = 0 mod p for k = 3, 5 "
              "(consistency)", ok_mom)
for d, N in ((4, 2), (6, 2)):
    p = primes_1_mod(lcm(d, 4))[0]
    R = Ring(d, p)
    L = liouvillian(hamiltonian("clock-field", N, R), N, R, gamma)
    check(f"control d={d} N={N}: a clock field h (Z + Z^dag) breaks the identity",
          not identity_residual_zero(L, N, R, gamma))

print("G4  odd d: the moment obstruction")
for d, N, kind in ((3, 1, "zero"), (3, 2, "zero"), (3, 2, "swap"), (3, 2, "clock-shift"),
                   (5, 1, "zero"), (5, 2, "swap")):
    sums = []
    for p in primes_1_mod(lcm(d, 4)):
        R = Ring(d, p)
        L = liouvillian(hamiltonian(kind, N, R), N, R, gamma)
        sums.append(odd_power_sums(L, N, R, gamma, ks=(3, 5, 7)))
    nz = any(v != 0 for s in sums for v in s.values())
    check(f"d={d} N={N} H={kind}: an odd centred power sum (k = 3, 5, 7) is nonzero mod p: "
          "no palindrome about any centre", nz)

print("G6  continuous phase noise n = diag(0..d-1): cell rates (i - j)^2 / 2")
for d in range(2, 9):
    rates = sorted(Fraction((i - j) ** 2, 2) for i in range(d) for j in range(d))
    sym = all(rates[k] + rates[-1 - k] == rates[0] + rates[-1] for k in range(len(rates)))
    check(f"d={d}: the cell rates are {'symmetric' if sym else 'not symmetric'} about a centre",
          sym == (d == 2))

print("G7  the hidden qubit in every even qudit")
for d, p in [(d, p) for d in (2, 4, 6, 8) for p in primes_1_mod(lcm(d, 4))]:
    R = Ring(d, p)
    X, Z, Zd = clock_ops(R)
    h = d // 2
    Zb = np.diag([1, p - 1]).astype(np.int64)
    Xb = np.array([[0, 1], [1, 0]], dtype=np.int64)
    Zp = np.diag([pow(R.w, j, p) for j in range(h)]).astype(np.int64)
    ok_z = np.array_equal(np.kron(Zb, Zp) % p, Z)
    ok_x = np.array_equal(np.kron(Xb, np.eye(h, dtype=np.int64)), np.linalg.matrix_power(X, h))
    Xh = np.linalg.matrix_power(X, h)
    ok_anti = np.array_equal((Xh @ Z @ Xh.T) % p, (-Z) % p)
    check(f"d={d} p={p}: Z_d = Z_b (x) diag(w^j'), X^(d/2) = X_b (x) 1, and X^(d/2) anticommutes with Z_d",
          ok_z and ok_x and ok_anti)

print("G8  d = 4 is two qubits: the clock pair is Pauli-string dephasing")
for p in primes_1_mod(4):
    R = Ring(4, p)
    X, Z, Zd = clock_ops(R)
    clock_D = (np.kron(Zd, Z) + np.kron(Z, Zd) - 2 * np.eye(16, dtype=np.int64)) % p
    Z1 = np.kron(np.diag([1, p - 1]), np.eye(2, dtype=np.int64)) % p       # Z on the high qubit b
    Z1Z0 = np.kron(np.diag([1, p - 1]), np.diag([1, p - 1])) % p
    pauli_D = (np.kron(Z1, Z1) + np.kron(Z1Z0, Z1Z0) - 2 * np.eye(16, dtype=np.int64)) % p
    check(f"d=4 p={p}: D[Z_4] + D[Z_4^dag] = D[Z (x) 1] + D[Z (x) Z], entry by entry",
          np.array_equal(clock_D, pauli_D))

print("G9  the one-sided element: U rho with U = (X^{d/2})^{(x)N}, H only commuting with U")
def with_imag_shift_field(H, N, R, h=4):
    """Adds h i (X - X^dag) on every site: Hermitian, commutes with X^{d/2}, not real."""
    d, p = R.d, R.p
    X, _, _ = clock_ops(R)
    for l in range(N):
        Xl = site(X, l, N, d, p)
        H = (H + h * R.i * (Xl - Xl.T)) % p
    return H % p


def left_mult_identity_zero(H, N, R, gamma):
    """(I (x) U) L (I (x) U^-1) = -L^dagger - 2 sigma, sigma = 2 N gamma; L^dagger is built from its
    formula: -H in the commutator, the adjoint jumps (the pair {Z, Z^dag} is its own adjoint set)."""
    d, p = R.d, R.p
    X, _, _ = clock_ops(R)
    U = kron_list([np.linalg.matrix_power(X, d // 2)] * N, p)
    D = d ** N
    S = np.kron(np.eye(D, dtype=np.int64), U)
    Si = np.kron(np.eye(D, dtype=np.int64), U.T)
    L = liouvillian(H, N, R, gamma)
    Ldag = liouvillian((-H) % p, N, R, gamma)
    T = (S @ L % p @ Si + Ldag + 4 * N * gamma * np.eye(D * D, dtype=np.int64)) % p
    return not T.any(), L


for d, N, kind in ((4, 2, "swap"), (4, 2, "clock-shift"), (6, 2, "swap"), (4, 1, "clock-shift")):
    ok_l, ok_t, ok_m = True, True, True
    pr = primes_1_mod(lcm(d, 4))
    ok_l &= zero_is_zero(majorant_max(kind, N, d, gamma, 2 * N, imag_field=True), 4 * N * gamma, d, pr)
    for p in pr:
        R = Ring(d, p)
        H = with_imag_shift_field(hamiltonian(kind, N, R), N, R)
        z, L = left_mult_identity_zero(H, N, R, gamma)
        ok_l &= z
        ok_t &= not identity_residual_zero(L, N, R, gamma)       # the transpose mirror needs real H
        ok_m &= all(v == 0 for v in odd_power_sums(L, N, R, gamma).values())
    lab = "one site, i(X - X^dag) field alone" if N == 1 else f"{kind} + i(X - X^dag) field"
    check(f"d={d} N={N} {lab}: U rho reflects L onto "
          "-L^dag - 4N gamma, exact (norm bound)", ok_l)
    check(f"d={d} N={N} {lab}: the transpose mirror rho^T U "
          "does not (H not real)", ok_t)
    ok_c = True
    for p in pr:
        R = Ring(d, p)
        Hc = hamiltonian("clock-field", N, R) if N > 1 else hamiltonian("clock-field", 1, R)
        ok_c &= not left_mult_identity_zero(Hc, N, R, gamma)[0]
    check(f"control d={d} N={N}: with a clock field instead, the one-sided identity fails", ok_c)
    check(f"d={d} N={N} {lab}: the odd centred power sums "
          "vanish mod p (consistency)", ok_m)

print("G10 the single clock jump Z, without Z^dag: U composed with the parity P|j> = |-j>")
def liouvillian_single(H, N, R, gamma):
    d, p = R.d, R.p
    _, Z, Zd = clock_ops(R)
    D = d ** N
    I = np.eye(D, dtype=np.int64)
    L = (((-R.i) % p) * (np.kron(I, H) - np.kron(H.T, I))) % p
    for l in range(N):
        L = (L + gamma * np.kron(site(Zd, l, N, d, p), site(Z, l, N, d, p))) % p
    return (L - N * gamma * np.eye(D * D, dtype=np.int64)) % p


def liouvillian_single_dagger(H, N, R, gamma):
    """L^dagger from its formula: -H in the commutator, the jump Z^dag in place of Z."""
    d, p = R.d, R.p
    _, Z, Zd = clock_ops(R)
    D = d ** N
    I = np.eye(D, dtype=np.int64)
    Hm = (-H) % p
    L = (((-R.i) % p) * (np.kron(I, Hm) - np.kron(Hm.T, I))) % p
    for l in range(N):
        L = (L + gamma * np.kron(site(Z, l, N, d, p), site(Zd, l, N, d, p))) % p
    return (L - N * gamma * np.eye(D * D, dtype=np.int64)) % p


def single_identity_zero(H, N, R, gamma):
    """(L_U . Ad_Q) L (L_U . Ad_Q)^-1 = -L^dagger - 2 sigma, sigma = N gamma, Q = P^{(x)N}:
    Ad_Q turns the jump Z into Z^dag and keeps H when [Q, H] = 0; L_U then reverses it."""
    d, p = R.d, R.p
    X, _, _ = clock_ops(R)
    D = d ** N
    P1 = np.zeros((d, d), dtype=np.int64)
    for j in range(d):
        P1[(-j) % d, j] = 1
    Q = kron_list([P1] * N, p)
    U = kron_list([np.linalg.matrix_power(X, d // 2)] * N, p)
    M = np.kron(np.eye(D, dtype=np.int64), U) @ np.kron(Q, Q) % p      # Q real: conj(Q) = Q
    Mi = np.kron(Q.T, Q.T) @ np.kron(np.eye(D, dtype=np.int64), U.T) % p
    L = liouvillian_single(H, N, R, gamma)
    T = (M @ L % p @ Mi + liouvillian_single_dagger(H, N, R, gamma)
         + 2 * N * gamma * np.eye(D * D, dtype=np.int64)) % p
    return not T.any(), L


def single_odd_power_sums(L, N, R, gamma, ks=(3, 5)):
    p, n = R.p, len(L)
    assert int(np.trace(L) % p) == (-n * N * gamma) % p, "trace is not -dim * N gamma"
    A = (L + N * gamma * np.eye(n, dtype=np.int64)) % p
    A2, Pw, k, out = (A @ A) % p, A.copy(), 1, {}
    while k < max(ks):
        Pw, k = (Pw @ A2) % p, k + 2
        if k in ks:
            out[k] = int(np.trace(Pw) % p)
    return out


for d, N, kind in ((4, 2, "swap"), (4, 2, "clock-shift"), (4, 2, "shift-field"), (6, 2, "swap"),
                   (6, 2, "clock-shift")):
    ok_i, ok_m = True, True
    pr = primes_1_mod(lcm(d, 4))
    ok_i &= zero_is_zero(majorant_max(kind, N, d, gamma, N), 2 * N * gamma, d, pr)
    for p in pr:
        R = Ring(d, p)
        z, L = single_identity_zero(hamiltonian(kind, N, R), N, R, gamma)
        ok_i &= z
        ok_m &= all(v == 0 for v in single_odd_power_sums(L, N, R, gamma).values())
    check(f"d={d} N={N} {kind}, jump Z alone: U . Ad_P reflects L onto -L^dag - 2N gamma, "
          "exact (norm bound)", ok_i)
    check(f"d={d} N={N} {kind}, jump Z alone: Tr (L + N gamma)^k = 0 mod p for k = 3, 5 "
          "(consistency)", ok_m)
for d in (4, 6):
    ok_break, ok_m = True, True
    for p in primes_1_mod(lcm(d, 4)):
        R = Ring(d, p)
        H = with_imag_shift_field(hamiltonian("swap", 2, R), 2, R)
        z, L = single_identity_zero(H, 2, R, gamma)
        ok_break &= not z
        ok_m &= all(v == 0 for v in single_odd_power_sums(L, 2, R, gamma).values())
    check(f"control d={d} N=2 swap + i(X - X^dag), jump Z alone: the field breaks the parity route's "
          "identity (the identity is what is checked)", ok_break)
    check(f"reading d={d} N=2 swap + i(X - X^dag), jump Z alone: the odd centred power sums still "
          "vanish mod p (consistency; G12 gives the element)", ok_m)

print("G12 the lone jump through the transpose: M(rho) = U Q rho^T Q^-1, PH-bar P = H")
def transpose_parity_perm(N, d):
    """|r><c| -> U Q |c><r| Q^-1 = |u(q(c))><q(r)|, q(j) = -j per site, u the shift by d/2."""
    D = d ** N
    s = d // 2
    def q(j):
        dig = np.unravel_index(j, (d,) * N)
        return int(np.ravel_multi_index([(-x) % d for x in dig], (d,) * N))
    def u(j):
        dig = np.unravel_index(j, (d,) * N)
        return int(np.ravel_multi_index([(x + s) % d for x in dig], (d,) * N))
    P = np.empty(D * D, dtype=np.int64)
    for r, c in product(range(D), range(D)):
        P[c * D + r] = q(r) * D + u(q(c))
    return P


def with_imag_clock_bond(H, N, R, h=3):
    """Adds h i (Z_1 Z_2 - h.c.): Hermitian, commutes with U, and P H-bar P = -(this term)."""
    d, p = R.d, R.p
    _, Z, Zd = clock_ops(R)
    B = site(Z, 0, N, d, p) @ site(Z, 1, N, d, p) % p
    Bd = site(Zd, 0, N, d, p) @ site(Zd, 1, N, d, p) % p
    return (H + h * R.i * (B - Bd)) % p


def transpose_route_zero(H, N, R, gamma):
    P = transpose_parity_perm(N, R.d)
    L = liouvillian_single(H, N, R, gamma)
    M = np.empty_like(L)
    M[np.ix_(P, P)] = L
    return not ((M + L + 2 * N * gamma * np.eye(len(L), dtype=np.int64)) % R.p).any()


for d, N, kind in ((4, 2, "swap"), (4, 2, "clock-shift"), (6, 2, "swap")):
    ok_i, ok_v = True, True
    pr = primes_1_mod(lcm(d, 4))
    ok_i &= zero_is_zero(majorant_max(kind, N, d, gamma, N, imag_field=True), 2 * N * gamma, d, pr)
    for p in pr:
        R = Ring(d, p)
        H = with_imag_shift_field(hamiltonian(kind, N, R), N, R)
        ok_i &= transpose_route_zero(H, N, R, gamma)
        ok_v &= not single_identity_zero(H, N, R, gamma)[0]
    check(f"d={d} N={N} {kind} + i(X - X^dag), jump Z alone: U Q rho^T Q^-1 reflects L onto "
          "-L - 2N gamma, exact (norm bound), where the parity route of G10 fails", ok_i and ok_v)
for d in (4, 6):
    p = primes_1_mod(lcm(d, 4))[0]
    R = Ring(d, p)
    H = with_imag_clock_bond(hamiltonian("swap", 2, R), 2, R)
    check(f"control d={d} N=2 swap + i(Z Z - h.c.), jump Z alone: P H-bar P = -H on that term, "
          "the transpose route breaks", not transpose_route_zero(H, 2, R, gamma))

print("G13 the numerical-range no-go: exact inequalities")
for d in range(2, 10):
    r = [2 * (1 - sp.cos(2 * sp.pi * k / d)) for k in range(d)]
    rmax = max(r, key=lambda x: float(x))
    mean = sp.nsimplify(sp.simplify(sum(r) / d))
    rel = sp.simplify(rmax - 2 * mean)
    if d % 2 == 0:
        check(f"d={d}: clock pair, r_max = 2 * mean = 4: the bound meets 0's partner exactly",
              rel == 0 and mean == 2)
    else:
        check(f"d={d}: clock pair, r_max = 2(1 + cos(pi/d)) < 4 = 2 * mean (the lone jump's "
              "self-adjoint part is half of it): no palindrome for any H", bool(rel < 0) and mean == 2)
    if d >= 3:
        check(f"d={d}: uniform full-Cartan, the rate r is below 2 r (d - 1)/d: no palindrome for any H",
              Fraction(1) < Fraction(2 * (d - 1), d))

print("G11 the conventions: the builders on matrix units at H = 0")
for d, N, p in [(d, N, p) for d, N in ((3, 1), (4, 1), (4, 2), (6, 2), (5, 1))
                for p in primes_1_mod(lcm(d, 4))]:
    R = Ring(d, p)
    D = d ** N
    H0 = np.zeros((D, D), dtype=np.int64)
    Lp = liouvillian(H0, N, R, gamma)
    Ls = liouvillian_single(H0, N, R, gamma)
    ok_p, ok_s = True, True
    for j, k in product(range(D), range(D)):
        dj, dk = np.unravel_index(j, (d,) * N), np.unravel_index(k, (d,) * N)
        idx = k * D + j                                   # column stacking: |j><k|
        ep = sum(pow(R.w, (a - b) % d, p) + pow(R.w, (b - a) % d, p) - 2 for a, b in zip(map(int, dj), map(int, dk))) % p
        es = sum(pow(R.w, (a - b) % d, p) - 1 for a, b in zip(map(int, dj), map(int, dk))) % p
        col_p, col_s = Lp[:, idx].copy(), Ls[:, idx].copy()
        ok_p &= col_p[idx] == (gamma * ep) % p
        ok_s &= col_s[idx] == (gamma * es) % p
        col_p[idx] = col_s[idx] = 0
        ok_p &= not col_p.any()
        ok_s &= not col_s.any()
    check(f"d={d} N={N} p={p}: the pair builder multiplies |j><k| by gamma sum (w^(j-k) + w^(k-j) - 2), the "
          "single builder by gamma sum (w^(j-k) - 1), and nothing else", ok_p and ok_s)

print()
if FAILS:
    print(f"qudit clock mirror gate: {len(FAILS)} of {COUNT} FAILED")
    sys.exit(1)
print(f"qudit clock mirror gate: ALL GREEN, {COUNT} checks")
