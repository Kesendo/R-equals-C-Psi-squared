"""The transport formula for the odd centred traces under loss and gain, and the seventh trace's part of
degree two in the rates (F137; docs/proofs/PROOF_TRANSPORT_FORMULA_LOSS_AND_GAIN.md).

Setting and conventions of simulations/fifth_moment_flux_area.py, whose generator this gate imports:
H = sum over bonds of J(XX + YY) + D_z(X_a Y_b - Y_a X_b) + Delta Z_a Z_b on any graph, t_ab = J + i D_z the
amplitude of the hop a -> b, loss gamma_down and gain gamma_up per site, no fields, no dephasing,
M = 2L + c, s = gamma_down + gamma_up, a = gamma_down - gamma_up.

  (T) the transport formula, order by order: for l != k the coefficient of s_l a_k in Tr(M^m) is
          m * sum over p + q + p' + q' = m - 2 of (-2i)^(p+q) (2i)^(p'+q') (p+p')! (q+q')! / (p! q! p'! q'!)
              * (A+_pq conj(A+_p'q') - A-_pq conj(A-_p'q')),
      A+_pq = Tr(s+_l H^p s-_k H^q), A-_pq = Tr(s-_l H^p s+_k H^q): the Taylor coefficients of
      t * integral_0^t (|f+|^2 - |f-|^2) dx. Checked at m = 5 and m = 7 for every ordered pair, the left side
      read off the generator with s = 2 e_l, a = 2 e_k.
  (D) the degrees: Tr(M^7) has no part of degree 0 or 6 in the rates
  (L) the part of degree two of Tr(M^7) = sum over ordered pairs l != k of s_l a_k B_lk,
          B_lk = 10752 * 4^N * [ sum_j Im(t_lk t_kj t_jl) * (W_lkj + Q_lkj)
                                 - 15 sum_{i,j} |t_li|^2 Im(t_ik t_kj t_ji) - 15 sum_{i,j} |t_kj|^2 Im(t_lj t_ji t_il)
                                 + 15 sum Im(t_lk t_kd t_db t_ba t_al) - 15 sum Im(t_lm t_mk t_kj t_ji t_il)
                                 + 5 sum_{i,j} Im(t_lk t_kj t_ji t_il) (D_lk - D_li - D_ij - D_jk)
                                 + 5 sum_{m,j} Im(t_lm t_mk t_kj t_jl) (D_lj + D_jk) ],
          W_lkj = 30 S - 10 (d_l + d_k) - 5 d_j - 2 |t_lk|^2,   S = sum over bonds |t|^2,  d_v = sum_m |t_vm|^2,
          Q_lkj = 15 (ZZ bonds off the triangle) + 5 (one end at l or k, the other off it) + 10 (one end at j,
                  the other off it) + 5 (D_lj^2 + D_jk^2) + 10/3 D_lk^2 + 2 D_lj D_jk - 8/3 D_lk (D_lj + D_jk),
      the first three counting D^2 of each such bond, D_ab the coefficient of Z_a Z_b, and the inner sites
      distinct and off {l, k}. Against the generator on graphs at N = 3, 4, 5.
  (R) recorded values: the unit triangle at flux pi/2 under a temperature gradient at one coupling, where
      Tr((L + c/2)^7) = 430080 is all of degree four, the part of degree two being zero on that line of rate
      points; a ring of five, where the part of degree two is 2520 * 4^N * Im(loop) * (pentagram - pentagon);
      the two settings at N = 4 whose fifth trace the fifth trace's gate finds zero (a triangle at one ratio
      with a fourth site at another, a ring of four with flux, ZZ and four ratios): their seventh trace is
      not zero, all of it of degree two, and (L) gives it; a triangle with three unequal weights and three
      rate points on a line off the origin where the part of degree two is zero all the same; gate M11's own
      ring of four (cube_moves_gate.py), Tr((L + c/2)^7) = -20380416/25, all of degree two, zero without its ZZ
  (V) the closed form is held as 21 integer coefficients (the ZZ ones at three times their value, each term
      of the two linear forms of the rings of four on its own); each of them moved by one is noticed by a row
  (C) controls: with a longitudinal field on one site Tr(M^7) is no longer even in the rates; an XX term on
      one bond moves the part of degree two off (L)

numpy alone, integer inputs. Every comparison is exact: the traces are taken modulo ten primes (float64
holding integers below 2^53) and joined by the Chinese remainder theorem, the degrees in the rates are
separated by scaling the rates, the transport sum is a sum of Gaussian rationals in Fractions.
"""
import itertools
import math
import random
import sys
from fractions import Fraction

import numpy as np

import fifth_moment_flux_area as g


def primes(count, top=1 << 20):
    out, n = [], top - 1
    while len(out) < count:
        if n % 4 == 1 and all(n % q for q in range(3, int(n ** 0.5) + 1, 2)):
            out.append(n)
        n -= 2
    return out


PRIMES = primes(10)
MODULUS = math.prod(PRIMES)


def sqrt_minus_one(p):
    for x in range(2, p):
        r = pow(x, (p - 1) // 4, p)
        if r * r % p == p - 1:
            return r


def trace_power(M, m):
    """Tr(M^m) for a matrix of Gaussian integers, m = 5 or 7, exact while |Tr| < MODULUS / 2."""
    assert np.array_equal(M.real, np.round(M.real)) and np.array_equal(M.imag, np.round(M.imag))
    residues = []
    for p in PRIMES:
        i = sqrt_minus_one(p)
        A = ((M.real.astype(np.int64) + i * M.imag.astype(np.int64)) % p).astype(np.float64)
        assert A.shape[0] * (p - 1) ** 2 < 2 ** 53
        P2 = np.mod(A @ A, p)
        P3 = np.mod(A @ P2, p)
        right = np.mod(P2 @ P2, p) if m == 7 else P2
        tr = 0
        for r in range(A.shape[0]):
            tr = (tr + int(np.mod(P3[r] @ right[:, r], p))) % p
        residues.append(tr)
    x, mod = 0, 1
    for r, p in zip(residues, PRIMES):
        x += mod * ((r - x) * pow(mod, -1, p) % p)
        mod *= p
    return x - mod if x > mod // 2 else x


def solve(rows, rhs):
    A = [[Fraction(x) for x in row] for row in rows]
    b = [Fraction(x) for x in rhs]
    for c in range(len(b)):
        piv = A[c][c]
        A[c] = [x / piv for x in A[c]]
        b[c] /= piv
        for r in range(len(b)):
            if r != c:
                f = A[r][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
                b[r] -= f * b[c]
    return b


def rate_degrees(n, H, s, a, m=7):
    """Tr(M^m) split by its degree in the rates: {0: ..., 2: ..., ..., m - 1: ...}, m odd."""
    degrees = list(range(0, m, 2))
    lams = list(range(1, len(degrees) + 1))
    vals = [trace_power(g.generator(n, H, [lam * x for x in s], [lam * x for x in a]), m) for lam in lams]
    parts = solve([[lam ** e for e in degrees] for lam in lams], vals)
    assert all(x.denominator == 1 for x in parts)
    return dict(zip(degrees, (int(x) for x in parts)))


# ---- (T) the transport sum -------------------------------------------------------------------------------

def gauss(z):
    assert z.real == round(z.real) and z.imag == round(z.imag)
    return (Fraction(int(round(z.real))), Fraction(int(round(z.imag))))


def gmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def gconj(x):
    return (x[0], -x[1])


def transport(n, H, l, k, m):
    """The coefficient of s_l a_k in Tr(M^m) by the transport formula, as a Gaussian rational."""
    powers = [np.eye(1 << n, dtype=complex)]
    for _ in range(m - 2):
        powers.append(powers[-1] @ H)
    amps = {}
    for sign, first, second in ((+1, g.SP, g.SM), (-1, g.SM, g.SP)):
        left, right = g.site(first, l, n), g.site(second, k, n)
        for p in range(m - 1):
            for q in range(m - 1 - p):
                amps[sign, p, q] = gauss(np.trace(left @ powers[p] @ right @ powers[q]))
    unit = [(Fraction(1), Fraction(0)), (Fraction(0), Fraction(1)), (Fraction(-1), Fraction(0)), (Fraction(0), Fraction(-1))]
    total = (Fraction(0), Fraction(0))
    for p, q, pp in itertools.product(range(m - 1), repeat=3):
        qq = m - 2 - p - q - pp
        if qq < 0:
            continue
        weight = Fraction(math.factorial(p + pp) * math.factorial(q + qq),
                          math.factorial(p) * math.factorial(q) * math.factorial(pp) * math.factorial(qq))
        # (-2i)^(p+q) (2i)^(pp+qq) = 2^(m-2) i^(3(p+q) + (pp+qq))
        phase = unit[(3 * (p + q) + pp + qq) % 4]
        diff = tuple(x - y for x, y in zip(gmul(amps[+1, p, q], gconj(amps[+1, pp, qq])),
                                           gmul(amps[-1, p, q], gconj(amps[-1, pp, qq]))))
        term = gmul(phase, diff)
        total = (total[0] + weight * term[0], total[1] + weight * term[1])
    scale = m * 2 ** (m - 2)
    return (scale * total[0], scale * total[1])


# ---- (L) the closed form ---------------------------------------------------------------------------------

CONSTANTS = dict(bond=30, touch=10, third=5, own=2, borrowed_l=15, borrowed_k=15, pent_near=15, pent_far=15,
                 zz_off=45, zz_lk_out=15, zz_j_out=30, zz_side=15, zz_own=10, zz_sides=6, zz_own_side=8,
                 sq_own=15, sq_li=15, sq_ij=15, sq_jk=15, sq_far_lj=15, sq_far_jk=15)


def closed_form(n, hop, zz, s, a, C=CONSTANTS):
    """The part of degree two of Tr(M^7), by (L). The ZZ constants are held at three times their value."""
    t = [[0j] * n for _ in range(n)]
    for (x, y), v in hop.items():
        t[x][y] += v
        t[y][x] += v.conjugate()
    D = [[0] * n for _ in range(n)]
    for (x, y), v in zz.items():
        D[x][y] += v
        D[y][x] += v
    im = lambda z: int(round(z.imag))
    sq = lambda z: int(round(abs(z) ** 2))
    S = sum(sq(t[x][y]) for x, y in itertools.combinations(range(n), 2))
    deg = [sum(sq(t[v][w]) for w in range(n)) for v in range(n)]
    total = 0           # three times the bracket
    for l, k in itertools.permutations(range(n), 2):
        rest = [v for v in range(n) if v not in (l, k)]
        B = 0
        for j in rest:
            loop = im(t[l][k] * t[k][j] * t[j][l])
            W = C["bond"] * S - C["touch"] * (deg[l] + deg[k]) - C["third"] * deg[j] - C["own"] * sq(t[l][k])
            Q = 0
            for x, y in itertools.combinations(range(n), 2):
                e = {x, y}
                if e == {l, k}:
                    w = C["zz_own"]
                elif e in ({l, j}, {j, k}):
                    w = C["zz_side"]
                elif not e & {l, k, j}:
                    w = C["zz_off"]
                elif j in e:
                    w = C["zz_j_out"]
                else:
                    w = C["zz_lk_out"]
                Q += w * D[x][y] ** 2
            Q += C["zz_sides"] * D[l][j] * D[j][k] - C["zz_own_side"] * D[l][k] * (D[l][j] + D[j][k])
            B += loop * (3 * W + Q)
        for i, j in itertools.permutations(rest, 2):
            B -= 3 * C["borrowed_l"] * sq(t[l][i]) * im(t[i][k] * t[k][j] * t[j][i])
            B -= 3 * C["borrowed_k"] * sq(t[k][j]) * im(t[l][j] * t[j][i] * t[i][l])
            near = im(t[l][k] * t[k][j] * t[j][i] * t[i][l])
            B += near * (C["sq_own"] * D[l][k] - C["sq_li"] * D[l][i] - C["sq_ij"] * D[i][j] - C["sq_jk"] * D[j][k])
            far = im(t[l][i] * t[i][k] * t[k][j] * t[j][l])
            B += far * (C["sq_far_lj"] * D[l][j] + C["sq_far_jk"] * D[j][k])
        for x, y, z in itertools.permutations(rest, 3):
            B += 3 * C["pent_near"] * im(t[l][k] * t[k][x] * t[x][y] * t[y][z] * t[z][l])
            B -= 3 * C["pent_far"] * im(t[l][x] * t[x][k] * t[k][y] * t[y][z] * t[z][l])
        total += s[l] * a[k] * B
    assert (10752 * total) % 3 == 0
    return 10752 * 4 ** n * total // 3


FAILED = []


def check(label, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {label}" + (f"  ({detail})" if detail else ""), flush=True)
    if not ok:
        FAILED.append(label)


ROWS = [
    ("triangle, |t|^2 = 5, 10, 5, ZZ 5, -2, 1", 3,
     {(0, 1): 1 + 2j, (1, 2): 3 - 1j, (0, 2): 2 + 1j}, {(0, 1): 5, (1, 2): -2, (0, 2): 1}, [2, 1, 3], [1, 3, 1]),
    ("triangle with a hopping pendant", 4,
     {(0, 1): 1 + 0j, (1, 2): 2 + 1j, (2, 0): 1j, (3, 1): 2 - 3j}, {}, [3, 1, 2, 5], [1, 2, 2, 1]),
    ("triangle, ZZ only on three bonds to a fourth site", 4,
     {(0, 1): 1 + 0j, (1, 2): 2 + 1j, (2, 0): 1j}, {(0, 3): 2, (1, 3): -3, (2, 3): 1}, [3, 1, 2, 5], [1, 2, 2, 1]),
    ("ring of four with flux, ZZ on its edges and both diagonals", 4,
     {(0, 1): 1 + 0j, (1, 2): 2 + 1j, (2, 3): 1 - 1j, (3, 0): 2j},
     {(0, 1): 3, (1, 2): -1, (2, 3): 2, (0, 3): 5, (0, 2): 7, (1, 3): -4}, [2, 1, 3, 1], [1, 3, 1, 2]),
    ("K4, complex hoppings, six different ZZ", 4,
     {(0, 1): 1 + 2j, (0, 2): 2 - 1j, (0, 3): 1 + 1j, (1, 2): 3 + 1j, (3, 1): -1 + 2j, (2, 3): 2 + 3j},
     {(0, 1): 1, (0, 2): -2, (0, 3): 3, (1, 2): 2, (1, 3): -1, (2, 3): 4}, [2, 1, 3, 1], [1, 3, 1, 2]),
    ("triangle with a ZZ bond off it, pendants and ZZ on its edges (N = 5)", 5,
     {(0, 1): 1 + 1j, (1, 2): 2 + 0j, (2, 0): 1j},
     {(3, 4): 3, (0, 3): 2, (2, 4): -1, (0, 1): 1, (1, 2): -2, (0, 2): 3}, [2, 1, 3, 1, 4], [1, 3, 1, 2, 1]),
    ("ring of five with a chord, ZZ on every hopping bond and on two others (N = 5)", 5,
     {(0, 1): 1 + 1j, (1, 2): 2 + 0j, (2, 3): 1 - 1j, (3, 4): 1j, (4, 0): 1 + 0j, (0, 2): 1 + 2j},
     {(0, 1): 2, (1, 2): -1, (2, 3): 3, (3, 4): 1, (0, 4): -2, (0, 2): 1, (1, 3): 2, (1, 4): -3},
     [2, 1, 3, 1, 4], [1, 3, 1, 2, 2]),
]


def main():
    rnd = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 7)

    # (T) the transport formula at m = 5 and 7, every ordered pair
    settings = [(n, hop, zz) for _, n, hop, zz, _, _ in ROWS[:1] + ROWS[3:5]]
    settings += [(n,) + g.random_setting(rnd, n)[:2] for n in (3, 3, 4)]
    seen = 0
    for n, hop, zz in settings:
        H = g.hamiltonian(n, hop, zz)
        ok, nonzero = True, 0
        for m in (5, 7):
            for l, k in itertools.permutations(range(n), 2):
                s = [2 if v == l else 0 for v in range(n)]
                a = [2 if v == k else 0 for v in range(n)]
                read = rate_degrees(n, H, s, a, m)[2]
                want = transport(n, H, l, k, m)
                ok = ok and want[1] == 0 and 4 * want[0] == read
                nonzero += read != 0
        check(f"(T) N={n}, {len(hop)} hopping and {len(zz)} ZZ bonds: the transport sum is the coefficient of "
              f"s_l a_k in Tr M^5 and Tr M^7, every ordered pair", ok, f"{nonzero} of them not zero")
        seen += nonzero
    check("(T) is not read on zeros alone (at least forty coefficients not zero)", seen >= 40, f"{seen}")

    # (D), (L) the closed form on the rows and on random graphs
    rows = list(ROWS)
    for n, trials in ((3, 2), (4, 3), (5, 2)):
        for trial in range(trials):
            pairs = list(itertools.combinations(range(n), 2))
            hop = {p: complex(rnd.randint(-2, 2), rnd.randint(-2, 2)) for p in pairs if rnd.random() < 0.9}
            hop = {p: v for p, v in hop.items() if v != 0}
            zz = {p: rnd.randint(-3, 3) for p in pairs if rnd.random() < 0.6}
            rows.append((f"random graph, N = {n}, {len(hop)} hopping and {len(zz)} ZZ bonds", n, hop, zz,
                         [2 * rnd.randint(1, 3) for _ in range(n)], [2 * rnd.randint(-2, 2) for _ in range(n)]))
    read_rows, silent = [], 0
    for name, n, hop, zz, s, a in rows:
        parts = rate_degrees(n, g.hamiltonian(n, hop, zz), s, a)
        want = closed_form(n, hop, zz, s, a)
        read_rows.append((n, hop, zz, s, a, parts[2]))
        silent += parts[2] == 0
        check(f"(D), (L) {name}: degree two {parts[2]}", parts[0] == 0 and parts[6] == 0 and parts[2] == want,
              f"closed form {want}; degree four {parts[4]}")
    check("the rows are not silent (at least ten of them with a part of degree two that is not zero)",
          len(rows) - silent >= 10, f"{len(rows) - silent} of {len(rows)}")

    # (V) every constant moved by one is noticed
    unnoticed = []
    for name in CONSTANTS:
        moved = dict(CONSTANTS)
        moved[name] += 1
        if all(closed_form(n, hop, zz, s, a, moved) == read for n, hop, zz, s, a, read in read_rows):
            unnoticed.append(name)
    check(f"(V) each of the {len(CONSTANTS)} integer coefficients of the closed form, moved by one, misses the generator on a row",
          not unnoticed, f"unnoticed: {unnoticed}" if unnoticed else "")

    # (R) recorded values
    unit = {(0, 1): 1 + 0j, (1, 2): 1 + 0j, (2, 0): 1j}
    dg, ug = [4, 6, 8], [2, 4, 6]                # gamma = 2, nbar = (1, 2, 3)
    sg, ag = [x + y for x, y in zip(dg, ug)], [x - y for x, y in zip(dg, ug)]
    parts = rate_degrees(3, g.hamiltonian(3, unit, {}), sg, ag)
    check("(R) a temperature gradient on the unit triangle at flux pi/2: Tr((L + c/2)^7) = 430080 is all of degree "
          "four, the part of degree two zero by the closed form and in the generator",
          parts[2] == 0 and closed_form(3, unit, {}, sg, ag) == 0 and parts[4] == 430080 * 128 and parts[0] == parts[6] == 0)
    ring = {(0, 1): 1 + 0j, (1, 2): 1 + 1j, (2, 3): 1 + 0j, (3, 4): 2 + 0j, (4, 0): 1 + 2j}
    d5, u5 = [3, 1, 4, 1, 5], [2, 7, 1, 8, 2]
    s5, a5 = [x + y for x, y in zip(d5, u5)], [x - y for x, y in zip(d5, u5)]
    loop = 1
    for x in range(5):
        loop *= g.t_of(ring, x, (x + 1) % 5)
    cross = lambda x, y: d5[x] * u5[y] - d5[y] * u5[x]      # twice the area of (origin, point x, point y)
    pentagon = sum(cross(x, (x + 1) % 5) for x in range(5))
    pentagram = sum(cross(x, (x + 2) % 5) for x in range(5))
    parts = rate_degrees(5, g.hamiltonian(5, ring, {}), s5, a5)
    want = 128 * 2520 * 4 ** 5 * int(loop.imag) * (pentagram - pentagon)
    check("(R) a ring of five: the part of degree two of Tr((L + c/2)^7) is 2520 * 4^N * Im(loop) * (pentagram - "
          "pentagon)", parts[2] == want and want != 0 and closed_form(5, ring, {}, s5, a5) == want, f"{parts[2] // 128}")

    skew = {(0, 1): 1 + 0j, (1, 2): 1 + 1j, (2, 0): 2j}       # three different weights: 98, 91, 77
    dl, ul = [3, 5, 2], [1, 3, 0]                             # gamma_up = gamma_down - 2
    sl, al = [x + y for x, y in zip(dl, ul)], [x - y for x, y in zip(dl, ul)]
    parts = rate_degrees(3, g.hamiltonian(3, skew, {}), sl, al)
    check("(R) three unequal weights and three rate points on a line that misses the origin: the part of degree two "
          "can still be zero, Tr((L + c/2)^7) = 322560 all of degree four",
          parts[2] == 0 and closed_form(3, skew, {}, sl, al) == 0 and parts[4] == 322560 * 128 and parts[0] == parts[6] == 0
          and len({Fraction(x, y) for x, y in zip(ul, dl)}) == 3)

    # the two settings of fifth_moment_flux_area.py with a fifth trace of zero, rates times 72
    d4, u4 = [36, 72, 108, 72], [9, 18, 27, 8]          # ratio 1/4 on the triangle, 1/9 on the fourth site
    s4, a4 = [x + y for x, y in zip(d4, u4)], [x - y for x, y in zip(d4, u4)]
    hop4, zz4 = {(0, 1): 1 + 0j, (1, 2): 1 + 0j, (2, 0): 1j, (3, 0): 1 + 0j}, {(0, 3): 1}
    parts = rate_degrees(4, g.hamiltonian(4, hop4, zz4), s4, a4)
    check("(R) a triangle at one ratio with a fourth site at another, joined by a hopping bond: the seventh trace is "
          "-44800, all of degree two, and the closed form gives it",
          parts[2] == closed_form(4, hop4, zz4, s4, a4) == -44800 * 128 * 72 ** 2 and parts[0] == parts[4] == parts[6] == 0,
          f"Tr((L + c/2)^7) = {Fraction(parts[2], 128 * 72 ** 2)}")
    ring4 = {(0, 1): 1 + 0j, (1, 2): 2 + 0j, (2, 3): 3 + 0j, (3, 0): 17j}
    zz4 = {(0, 1): 7, (1, 2): 7, (2, 3): 7, (0, 3): 7}
    d4, u4 = [36, 72, 108, 72], [9, 24, 20, 8]          # four ratios
    s4, a4 = [x + y for x, y in zip(d4, u4)], [x - y for x, y in zip(d4, u4)]
    parts = rate_degrees(4, g.hamiltonian(4, ring4, zz4), s4, a4)
    check("(R) a ring of four with flux pi/2, ZZ on every bond and four ratios: the seventh trace is -68239360, all "
          "of degree two, and the closed form gives it; without the ZZ it is zero",
          parts[2] == closed_form(4, ring4, zz4, s4, a4) == -68239360 * 128 * 72 ** 2 and parts[0] == parts[4] == parts[6] == 0
          and trace_power(g.generator(4, g.hamiltonian(4, ring4, {}), s4, a4), 7) == 0,
          f"Tr((L + c/2)^7) = {Fraction(parts[2], 128 * 72 ** 2)}")

    # gate M11's own ring of four (cube_moves_gate.py): hoppings 1, 2, 3, 17i/10, Delta = 7/10, four ratios;
    # H times 10 and the rates times 840, each part of the trace scaled back by its degrees
    ring11 = {(0, 1): 10 + 0j, (1, 2): 20 + 0j, (2, 3): 30 + 0j, (3, 0): 17j}
    zz11 = {(0, 1): 7, (1, 2): 7, (2, 3): 7, (0, 3): 7}
    d11, u11 = [420, 840, 1260, 1680], [105, 280, 168, 120]     # 840 * (1/2, 1, 3/2, 2) and 840 * (1/8, 1/3, 1/5, 1/7)
    s11, a11 = [x + y for x, y in zip(d11, u11)], [x - y for x, y in zip(d11, u11)]
    parts = rate_degrees(4, g.hamiltonian(4, ring11, zz11), s11, a11)
    seventh = Fraction(parts[2], 128 * 10 ** 5 * 840 ** 2)
    check("(R) gate M11's ring of four at Delta = 7/10, flux pi/2 and four ratios: Tr((L + c/2)^7) = -20380416/25, "
          "all of degree two, by the closed form and in the generator; without the ZZ it is zero",
          seventh == Fraction(-20380416, 25) and parts[2] == closed_form(4, ring11, zz11, s11, a11)
          and parts[0] == parts[4] == parts[6] == 0
          and trace_power(g.generator(4, g.hamiltonian(4, ring11, {}), s11, a11), 7) == 0, f"{seventh}")

    # (C) control
    name, n, hop, zz, s, a = ROWS[0]
    Hf = g.hamiltonian(n, hop, zz) + 3 * g.site(g.Z, 0, n)
    neg_s, neg_a = [-x for x in s], [-x for x in a]
    check("(C) control, a longitudinal field on one site: Tr M^7 is no longer even in the rates",
          trace_power(g.generator(n, Hf, s, a), 7) != trace_power(g.generator(n, Hf, neg_s, neg_a), 7))
    # 3 X_0 X_1 and not 2: the part of 2 X_0 X_1 that conserves the number moves t_01 to 2 + 2i, the loop becomes
    # real and every odd trace is zero by a pairing, which would be a miss for another reason
    Hx = g.hamiltonian(n, hop, zz) + 3 * g.site(g.X, 0, n) @ g.site(g.X, 1, n)
    xx = rate_degrees(n, Hx, s, a)[2]
    check("(C) control, an XX term on one bond (the excitation number not conserved): the part of degree two "
          "is not zero and leaves the closed form", xx != 0 and xx != closed_form(n, hop, zz, s, a), f"{xx}")

    print("ALL PASS" if not FAILED else f"FAILED: {FAILED}")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
