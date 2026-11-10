"""Gate: the seventh centred trace at fourth order in the rates on the triangle, the line form of the seventh trace,
and the ninth on equal heights (docs/proofs/PROOF_SEVENTH_TRACE_DEGREE_FOUR_TRIANGLE.md). Integers, every
comparison exact; traces modulo primes p = 1 mod 4 joined by the Chinese remainder theorem.

Setting as in F137 and the transport formula's proof: hopping bonds with ZZ, loss gd_l and gain gu_l per site, no
fields, no dephasing, s_l = (gd_l + gu_l)/2, omega_ab = gd_a*gu_b - gd_b*gu_a, Im = Im(t01*t12*t20), the loop of the
triangle 0 -> 1 -> 2 -> 0 with t_ab the hop a -> b.

  (F) the part of Tr((L + c/2)^7) of fourth degree in the rates, fitted exactly over all 126 monomials of the gate's
      rates on the unit triangle, has exactly the coefficients of 13440*Im*sum_cyc omega_{l,l+1}*(s_l^2 + s_{l+1}^2
      + 2*s_{l+2}^2); on three other triangles (unequal |t|, other flux, ZZ) it is that form at random rates, and
      the seventh trace split into all eight degrees has only degrees two and four
  (W) the same form by a second road: closed walks on the Pauli letters, one arrow 2|Z><I|, three hops -i*ad_H,
      the heights diag(1, 0, 0, -1) summed as h_3, times 7: an identity on the unit triangle, the 144 walks counted,
      at random rates on the others
  (K) past the triangle: with spectator sites (bonds to the triangle, ZZ to it) and with two triangles, the part of
      fourth degree is the sum over the triangles of Theorem 1's form 4^N*Im*sum_cyc omega_{l,l+1}*(630*|s|^2 - 420*(s_l^2 +
      s_{l+1}^2) - 210*s_{l+2}^2), |s|^2 over every site
  (G) on a line of rate points gu = e + m*gd the whole seventh trace is -e*[3360*Im*(1+m)^2*V(gd) +
      10752*<gd, b>], V the Vandermonde (x0-x1)(x1-x2)(x2-x0), b_l = sum_k beta_lk the flow of Theorem B's beta; on a
      vertical line gd = e it is e*[3360*Im*V(gu) + 10752*<gu, b>]; t01 = t12 = 2+i, t02 = 3-4i with ZZ 2, 2, 3 as an
      example, b = Im*(473/3, 0, -473/3), the constant 7568/15, a point of the cubic curve; F137's 430080 of the
      temperature gradient and the transport proof's 322560 as Theorem 1 and Theorem 2; a line at m = -1, b not zero;
      the shear d -> d + s keeps the fifth and seventh traces and moves the ninth
  (S) a triangle with equal weights, without the full symmetry of the sites (t01 = t12 = 11+3i, t20 = 5, ZZ 21 on (0, 2)): b = 0, and on
      equal heights at S = 4, arrows (3, -1, 2), the fifth and seventh traces zero and the ninth -53828247306240
      and on equal heights the whole ninth trace -6082560*S*(d0 - d2)*Q(S, d), Q an indefinite quadratic
  (E) on the triangle with equal |t| = r and one Delta, on equal heights gd_l + gu_l = 2S with arrows
      d_l = (gd_l - gu_l)/2: Tr((L + c/2)^9) = 36864*Im*(r^2 + 7*Delta^2)*S*V(d), and with two equal rate points the
      seventh and ninth traces are zero (the site-permutation carrier)
  (V) each integer constant above, moved by one, is noticed by a row

Run: python simulations/seventh_moment_degree_four.py [seed]
"""
import itertools
import random
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import seventh_moment_transport as T  # noqa: E402
g = T.g

FAILS = []


def check(label, ok, detail=""):
    print(("[PASS] " if ok else "[FAIL] ") + label + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(label)


# ---- exact traces of any power ----------------------------------------------------------------------------

def _primes(count, top=1 << 20):
    out, p = [], top
    while len(out) < count:
        p -= 1
        if p % 4 == 1 and all(p % q for q in range(2, int(p ** 0.5) + 1)):
            out.append(p)
    return out


PR = _primes(40)


def _trace_mod(M, m, p):
    i = T.sqrt_minus_one(p)
    A = ((M.real.astype(np.int64) + i * M.imag.astype(np.int64)) % p).astype(np.float64)
    assert A.shape[0] * (p - 1) ** 2 < 2 ** 53
    R, base, e = None, A, m
    while e:
        if e & 1:
            R = base if R is None else np.mod(R @ base, p)
        e >>= 1
        if e:
            base = np.mod(base @ base, p)
    return int(round(float(np.trace(R)))) % p


def trace_power(M, m):
    """Tr(M^m) of a Gaussian-integer matrix; the CRT over 24 and over 30 primes must agree."""
    assert np.array_equal(M.real, np.round(M.real)) and np.array_equal(M.imag, np.round(M.imag))

    def crt(k):
        x, mod = 0, 1
        for p in PR[:k]:
            r = _trace_mod(M, m, p)
            x += mod * ((r - x) * pow(mod, -1, p) % p)
            mod *= p
        return x - mod if x > mod // 2 else x
    a = crt(24)
    assert a == crt(30)
    return a


def centred_trace(n, hop, zz, gd, gu, m, scale=1):
    """Tr((L + c/2)^m) at integer gd, gu (and integer hops, ZZ); scale multiplies H and the rates (Tr scales by
    scale^m), so a rational point is reached through integers."""
    H = g.hamiltonian(n, {k: v * scale for k, v in hop.items()}, {k: v * scale for k, v in zz.items()})
    s = [scale * (x + y) for x, y in zip(gd, gu)]
    a = [scale * (x - y) for x, y in zip(gd, gu)]
    return Fraction(trace_power(g.generator(n, H, s, a), m), 2 ** m * scale ** m)


def rate_parts(n, hop, zz, gd, gu, m):
    """Tr((L + c/2)^m) split by degree in the rates (even degrees 0 .. m-1)."""
    degs = list(range(0, m, 2))
    vals = [centred_trace(n, hop, zz, [l * x for x in gd], [l * x for x in gu], m) for l in range(1, len(degs) + 1)]
    sol = T.solve([[l ** d for d in degs] for l in range(1, len(degs) + 1)], vals)
    return dict(zip(degs, sol))


# ---- the closed forms, constants held by name --------------------------------------------------------------

C = dict(r4=13440, own=1, opp=2, prod=630, kown=420, kopp=210, line4=3360, line2=10752, e9=36864, e9zz=7)


def gauss_t(hop, a, b):
    """t_ab as an exact pair (re, im) of Fractions; the hops are Gaussian integers, so the floats are exact."""
    t = g.t_of(hop, a, b)
    t = complex(t)
    assert t.real == int(t.real) and t.imag == int(t.imag)
    return Fraction(int(t.real)), Fraction(int(t.imag))


def gmul(x, y):
    return x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0]


def loop_im(hop):
    """Im(t01*t12*t20), exactly."""
    return gmul(gmul(gauss_t(hop, 0, 1), gauss_t(hop, 1, 2)), gauss_t(hop, 2, 0))[1]


def t2(hop, a, b):
    re, im = gauss_t(hop, a, b)
    return re * re + im * im


def omega(gd, gu, a, b):
    return gd[a] * gu[b] - gd[b] * gu[a]


def r4_form(hop, gd, gu, c=C):
    s = [Fraction(x + y, 2) for x, y in zip(gd, gu)]
    return c["r4"] * loop_im(hop) * sum(omega(gd, gu, l, (l + 1) % 3) * (c["own"] * (s[l] ** 2 + s[(l + 1) % 3] ** 2)
                                                                        + c["opp"] * s[(l + 2) % 3] ** 2) for l in range(3))


def triangle_sum_form(n, hop, gd, gu, c=C):
    """The part of fourth degree on any graph: a sum over its triangles (a, b, d), each with the loop a -> b -> d -> a,
    |s|^2 over every site of the graph."""
    s = [Fraction(x + y, 2) for x, y in zip(gd, gu)]
    all2 = sum(x * x for x in s)
    total = Fraction(0)
    for tri in itertools.combinations(range(n), 3):
        if any(t2(hop, tri[i], tri[(i + 1) % 3]) == 0 for i in range(3)):
            continue
        im = gmul(gmul(gauss_t(hop, tri[0], tri[1]), gauss_t(hop, tri[1], tri[2])), gauss_t(hop, tri[2], tri[0]))[1]
        total += im * sum(omega(gd, gu, tri[l], tri[(l + 1) % 3])
                          * (c["prod"] * all2 - c["kown"] * (s[tri[l]] ** 2 + s[tri[(l + 1) % 3]] ** 2)
                             - c["kopp"] * s[tri[(l + 2) % 3]] ** 2) for l in range(3))
    return 4 ** n * total


def flow_per_im(hop, zz):
    """b_l/Im, b_l = sum_k beta_lk the flow of Theorem B's beta on a triangle alone (W + Q per pair, nothing else at
    N = 3), exactly."""
    D = lambda a, b: Fraction(zz.get((a, b), zz.get((b, a), 0)))
    S = t2(hop, 0, 1) + t2(hop, 1, 2) + t2(hop, 0, 2)
    deg = lambda v: sum(t2(hop, v, w) for w in range(3) if w != v)
    b = [Fraction(0)] * 3
    for l in range(3):
        k, j = (l + 1) % 3, (l + 2) % 3
        W = 30 * S - 10 * (deg(l) + deg(k)) - 5 * deg(j) - 2 * t2(hop, l, k)
        Q = 5 * (D(l, j) ** 2 + D(j, k) ** 2) + Fraction(10, 3) * D(l, k) ** 2 + 2 * D(l, j) * D(j, k) - Fraction(8, 3) * D(l, k) * (D(l, j) + D(j, k))
        b[l] += W + Q
        b[k] -= W + Q
    return b


def vandermonde(x):
    return (x[0] - x[1]) * (x[1] - x[2]) * (x[2] - x[0])


def line_form(hop, zz, gd, e, m, c=C):
    bh = flow_per_im(hop, zz)
    return -e * loop_im(hop) * (c["line4"] * (1 + m) ** 2 * vandermonde(gd) + c["line2"] * sum(x * y for x, y in zip(gd, bh)))


def vertical_form(hop, zz, gu, e, c=C):
    bh = flow_per_im(hop, zz)
    return e * loop_im(hop) * (c["line4"] * vandermonde(gu) + c["line2"] * sum(x * y for x, y in zip(gu, bh)))


def ninth_equal_heights(hop, r2, delta, S, arrows, c=C):
    return c["e9"] * loop_im(hop) * (r2 + c["e9zz"] * delta ** 2) * S * vandermonde(arrows)


# ---- (W) the walk sum on the letters ------------------------------------------------------------------------

LET = {"I": np.eye(2), "X": np.array([[0, 1], [1, 0]], complex), "Y": np.array([[0, -1j], [1j, 0]]), "Z": np.diag([1., -1])}
HEIGHT = {"I": 1, "X": 0, "Y": 0, "Z": -1}
STRINGS = ["".join(p) for p in itertools.product("IXYZ", repeat=3)]


def walks(hop, zz):
    """Every closed walk P -> (arrow at k) -> P' -> Q1 -> Q2 -> P with three steps of -i*ad_H, as (k, strings, amp);
    the matrix entries of -i*ad_H on normalised Pauli strings are Gaussian integers here, read exactly."""
    mats = {p: np.kron(np.kron(LET[p[0]], LET[p[1]]), LET[p[2]]) for p in STRINGS}
    H = g.hamiltonian(3, hop, zz)
    G = {}
    for q in STRINGS:
        img = -1j * (H @ mats[q] - mats[q] @ H)
        for p in STRINGS:
            v = complex(np.trace(mats[p] @ img) / 8)
            assert v.real == round(v.real) and v.imag == round(v.imag)
            if v != 0:
                G[(p, q)] = (int(round(v.real)), int(round(v.imag)))
    succ = {}
    for (p, q), v in G.items():
        succ.setdefault(q, []).append((p, v))
    out = []
    for k in range(3):
        for p in STRINGS:
            if p[k] != "I":
                continue
            p1 = p[:k] + "Z" + p[k + 1:]
            for q1, a1 in succ.get(p1, []):
                for q2, a2 in succ.get(q1, []):
                    a3 = G.get((p, q2))
                    if a3 is not None:
                        out.append((k, (p1, q1, q2, p), gmul(gmul(a3, a2), a1)))
    return out


def h3(vals):
    tot = 0
    for combo in itertools.combinations_with_replacement(range(len(vals)), 3):
        tot += vals[combo[0]] * vals[combo[1]] * vals[combo[2]]
    return tot


def walk_sum(walk_list, s, dl):
    """7 * sum over walks of amp * 2*delta_k * h_3(heights), as (real, imaginary) parts, for any ring of s, dl."""
    eta = lambda p: sum(s[l] * HEIGHT[p[l]] for l in range(3))
    re = im = 0
    for k, path, amp in walk_list:
        w = 7 * 2 * dl[k] * h3([eta(x) for x in path])
        re += amp[0] * w
        im += amp[1] * w
    return re, im


# ---- fitting ------------------------------------------------------------------------------------------------

def solve_exact(A, b):
    """Exact elimination with pivoting on an overdetermined system; every surplus row must be consistent."""
    n = len(A[0])
    M = [[Fraction(v) for v in row] + [Fraction(r)] for row, r in zip(A, b)]
    piv, r = [], 0
    for c in range(n):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        inv = 1 / M[r][c]
        M[r] = [v * inv for v in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [vi - f * vr for vi, vr in zip(M[i], M[r])]
        piv.append(c)
        r += 1
    assert r == n, f"rank {r} < {n}"
    assert all(M[i][n] == 0 for i in range(r, len(M))), "inconsistent rows"
    sol = [Fraction(0)] * n
    for i, c in enumerate(piv):
        sol[c] = M[i][n]
    return sol


MONOS = list(itertools.combinations_with_replacement(range(6), 4))


def fit_degree_four(hop, zz, rnd, npts=140):
    rows, rhs = [], []
    for _ in range(npts):
        gd = [rnd.randint(1, 9) for _ in range(3)]
        gu = [rnd.randint(0, 9) for _ in range(3)]
        s2 = [x + y for x, y in zip(gd, gu)]
        a2 = [x - y for x, y in zip(gd, gu)]
        part = T.rate_degrees(3, g.hamiltonian(3, hop, zz), s2, a2)[4]      # Tr(M^7), M = 2L + c, gate rates
        x = s2 + a2
        rows.append([x[mono[0]] * x[mono[1]] * x[mono[2]] * x[mono[3]] for mono in MONOS])
        rhs.append(part)
    return dict(zip(MONOS, solve_exact(rows, rhs)))


def form_coefficients(hop):
    """The closed form in the gate's rates s' = gd + gu, a' = gd - gu, times 128 (Tr M^7 against Tr((L + c/2)^7)),
    as a dict over MONOS, by exact interpolation of the form itself on the same monomials."""
    import sympy as sp
    sv = sp.symbols("s0:3")
    av = sp.symbols("a0:3")
    gdv = [(sv[i] + av[i]) / 2 for i in range(3)]
    guv = [(sv[i] - av[i]) / 2 for i in range(3)]
    sh = [sv[i] / 2 for i in range(3)]
    li = loop_im(hop)
    im = sp.Rational(li.numerator, li.denominator)
    om = lambda a, b: gdv[a] * guv[b] - gdv[b] * guv[a]
    form = 128 * C["r4"] * im * sum(om(l, (l + 1) % 3) * (C["own"] * (sh[l] ** 2 + sh[(l + 1) % 3] ** 2)
                                                          + C["opp"] * sh[(l + 2) % 3] ** 2) for l in range(3))
    poly = sp.Poly(sp.expand(form), *sv, *av)
    out = {mono: Fraction(0) for mono in MONOS}
    for exps, coef in poly.terms():
        mono = tuple(sorted(sum(([i] * e for i, e in enumerate(exps)), [])))
        out[mono] = Fraction(int(coef.p), int(coef.q))
    return out


# ---- main ---------------------------------------------------------------------------------------------------

def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    rnd = random.Random(seed)
    unit = {(0, 1): 1 + 0j, (1, 2): 1 + 0j, (2, 0): 1j}
    tris = [("unequal |t|, ZZ", {(0, 1): 2 + 1j, (1, 2): 1 - 1j, (2, 0): 3 + 0j}, {(0, 1): 2, (1, 2): -1, (0, 2): 3}),
            ("equal |t|^2 = 5, ZZ 2", {(0, 1): 2 + 1j, (1, 2): 2 + 1j, (2, 0): 2 + 1j}, {(0, 1): 2, (1, 2): 2, (0, 2): 2}),
            ("a pure DM bond", {(0, 1): 1 + 0j, (1, 2): 2 + 0j, (2, 0): 3j}, {(1, 2): 1})]

    # (F) the fit on the unit triangle, coefficient for coefficient, then the form on other triangles
    fitted = fit_degree_four(unit, {}, rnd)
    want = form_coefficients(unit)
    check("(F) the part of fourth degree on the unit triangle, fitted over all 126 monomials, has exactly the "
          "coefficients of 13440*Im*sum_cyc omega_{l,l+1}*(s_l^2 + s_{l+1}^2 + 2*s_{l+2}^2)",
          fitted == want and any(v != 0 for v in want.values()),
          f"{sum(1 for v in want.values() if v != 0)} nonzero coefficients")
    rows_f = []
    for name, hop, zz in tris:
        ok, seen = True, False
        for _ in range(3):
            gd = [rnd.randint(1, 9) for _ in range(3)]
            gu = [rnd.randint(0, 9) for _ in range(3)]
            meas = rate_parts(3, hop, zz, gd, gu, 7)[4]
            rows_f.append((hop, gd, gu, meas))
            ok &= meas == r4_form(hop, gd, gu)
            seen |= meas != 0
        check(f"(F) {name}: the part of fourth degree is the same form at three random rate points", ok and seen)

    # the degrees the proof takes from F137's maps: all eight degrees 0 .. 7 split by eight scalings, only 2 and 4
    ok = True
    for name, hop, zz in tris:
        gd = [rnd.randint(1, 9) for _ in range(3)]
        gu = [rnd.randint(0, 9) for _ in range(3)]
        vals = [centred_trace(3, hop, zz, [l * x for x in gd], [l * x for x in gu], 7) for l in range(1, 9)]
        parts = T.solve([[l ** d for d in range(8)] for l in range(1, 9)], vals)
        ok &= all(parts[d] == 0 for d in (0, 1, 3, 5, 6, 7)) and parts[2] != 0 and parts[4] != 0
    check("(F) on the three other triangles the seventh trace split into all eight degrees 0 to 7: only degrees two "
          "and four are not zero", ok)

    # (W) the walk sum: as a polynomial on the unit triangle, at random rates on the other three
    import sympy as sp
    sv, dv = sp.symbols("s0:3"), sp.symbols("d0:3")
    wl = walks(unit, {})
    re, im = walk_sum(wl, sv, dv)
    form = 26880 * sum(dv[l] * (sv[(l + 1) % 3] - sv[(l + 2) % 3]) * (sv[0] ** 2 + sv[1] ** 2 + sv[2] ** 2
                                                                      - sv[(l + 1) % 3] * sv[(l + 2) % 3]) for l in range(3))
    site0 = [w for w in wl if w[0] == 0]
    re0, _ = walk_sum(site0, sv, dv)
    check("(W) on the unit triangle the walk sum (one arrow, three hops, h_3 of the heights, times 7) is Theorem 1 as a "
          "polynomial; 144 walks, 48 per arrow site, each of amplitude +-8; the arrow at site 0 alone gives "
          "7*2*8*240*d0*(s1 - s2)(s0^2 + s1^2 + s2^2 - s1*s2)",
          sp.expand(re - form) == 0 and im == 0 and len(wl) == 144 and len(site0) == 48
          and all(abs(a[0]) == 8 and a[1] == 0 for _, _, a in wl)
          and sp.expand(re0 - 7 * 2 * 8 * 240 * dv[0] * (sv[1] - sv[2]) * (sv[0] ** 2 + sv[1] ** 2 + sv[2] ** 2 - sv[1] * sv[2])) == 0)
    ok = True
    for name, hop, zz in tris:
        wl = walks(hop, zz)
        for _ in range(2):
            gd = [rnd.randint(1, 9) for _ in range(3)]
            gu = [rnd.randint(0, 9) for _ in range(3)]
            re, im = walk_sum(wl, [Fraction(x + y, 2) for x, y in zip(gd, gu)], [Fraction(x - y, 2) for x, y in zip(gd, gu)])
            ok &= im == 0 and re == r4_form(hop, gd, gu)
    check("(W) on the three other triangles the walk sum is the form at random rates, the flux-even walks cancelling", ok)

    # (K) the part of fourth degree on graphs past the triangle: spectators, two triangles
    specs = [("a site hanging on corner 0 by a hopping bond, ZZ to corner 2", 4,
              {(0, 1): 1 + 0j, (1, 2): 2 + 1j, (2, 0): 1j, (3, 0): 1 + 1j}, {(0, 1): 1, (2, 3): 2}),
             ("two triangles sharing the edge 0-1, ZZ across", 4,
              {(0, 1): 1 + 1j, (1, 2): 1 + 0j, (2, 0): 2j, (1, 3): 2 + 0j, (3, 0): 1 - 1j}, {(2, 3): 1, (0, 1): 2}),
             ("two spectators, one by hopping, one by ZZ alone", 5,
              {(0, 1): 1 + 1j, (1, 2): 1 + 0j, (2, 0): 2j, (3, 1): 1 + 0j}, {(2, 4): 1, (0, 2): 1})]
    rows_k = []
    for name, n, hop, zz in specs:
        ok, seen = True, False
        for _ in range(2 if n == 5 else 3):
            gd = [rnd.randint(1, 6) for _ in range(n)]
            gu = [rnd.randint(0, 6) for _ in range(n)]
            meas = rate_parts(n, hop, zz, gd, gu, 7)[4]
            rows_k.append((n, hop, gd, gu, meas))
            ok &= meas == triangle_sum_form(n, hop, gd, gu)
            seen |= meas != 0
        ok &= seen
        check(f"(K) {name} (N = {n}): the part of fourth degree is the sum over the triangles of Theorem 1's form, the "
              "spectators' heights in the product alone", ok)

    # (G) the line form and the vertical form, unequal weights included
    rows_g = []
    ok = True
    for name, hop, zz in [("unit", unit, {})] + tris:
        for _ in range(2):
            e, m = rnd.randint(1, 4), rnd.randint(0, 3)
            gd = [rnd.randint(0, 6) for _ in range(3)]
            gu = [e + m * x for x in gd]
            meas = centred_trace(3, hop, zz, gd, gu, 7)
            rows_g.append(("line", hop, zz, gd, e, m, meas))
            ok &= meas == line_form(hop, zz, gd, e, m)
            gu2 = [rnd.randint(0, 6) for _ in range(3)]
            e2 = rnd.randint(1, 4)
            meas2 = centred_trace(3, hop, zz, [e2] * 3, gu2, 7)
            rows_g.append(("vertical", hop, zz, gu2, e2, None, meas2))
            ok &= meas2 == vertical_form(hop, zz, gu2, e2)
    check("(G) on lines gu = e + m*gd and on vertical lines gd = e the seventh trace is the line form, on all four "
          "triangles", ok and any(r[-1] != 0 for r in rows_g))
    old = ({(0, 1): 2 + 1j, (1, 2): 2 + 1j, (0, 2): 3 - 4j}, {(0, 1): 2, (1, 2): 2, (0, 2): 3})
    b = flow_per_im(*old)
    ratio = Fraction(C["line2"]) * b[0] / C["line4"]
    check("(G) t01 = t12 = 2+i, t02 = 3-4i, ZZ 2, 2, 3: b = Im*(473/3, 0, -473/3), and the seventh trace on a line is "
          "a multiple of e*(g0 - g2)*[(1+m)^2*(g0 - g1)(g1 - g2) - 7568/15]",
          b == [Fraction(473, 3), 0, Fraction(-473, 3)] and ratio == Fraction(7568, 15), f"{b}, {ratio}")
    # a point of that curve: m = 0, g1 - g2 = 1, g0 - g1 = 7568/15; scaled by 15 to integers
    gd = [Fraction(7568, 15) + 1, 1, 0]
    gu = [Fraction(1)] * 3
    sc = 15                                  # H and the rates together: the trace scales by 15^7, the zero stays
    hop15 = {k: sc * v for k, v in old[0].items()}
    zz15 = {k: sc * v for k, v in old[1].items()}
    seventh = centred_trace(3, hop15, zz15, [int(sc * x) for x in gd], [int(sc * x) for x in gu], 7)
    off = [int(sc * x) for x in gd]
    off[0] += 1                              # one site moved (moving all three is a translation, along the cylinder)
    near = centred_trace(3, hop15, zz15, off, [int(sc * x) for x in gu], 7)
    along = centred_trace(3, hop15, zz15, [int(sc * x) + 1 for x in gd], [int(sc * x) for x in gu], 7)
    check("(G) on a point of that curve the seventh trace is zero, also after moving every gd by one (the surface is a "
          "cylinder along (1, 1, 1)); moving one site by one, it is not", seventh == 0 and along == 0 and near != 0)

    # two recorded values by Theorem 1 and Theorem 2: F137's temperature gradient, the transport proof's 322560
    tg_gd, tg_gu = [4, 6, 8], [2, 4, 6]
    tg = centred_trace(3, unit, {}, tg_gd, tg_gu, 7)
    check("(G) F137's temperature gradient on the unit triangle (gd = 2(n+1), gu = 2n, n = 1, 2, 3): the seventh trace "
          "is 430080, Theorem 1 at one arrow, and the line form gives it", tg == 430080 == r4_form(unit, tg_gd, tg_gu)
          and tg == line_form(unit, {}, tg_gd, -2, 1), f"{tg}")
    hop322 = {(0, 1): 1 + 0j, (1, 2): 1 + 1j, (2, 0): 2j}
    t322 = centred_trace(3, hop322, {}, [3, 5, 2], [1, 3, 0], 7)
    check("(G) the transport proof's triangle t = 1, 1+i, 2i at gd = (3, 5, 2), gu = gd - 2: the seventh trace is 322560, "
          "the line form with <gd, b> = 0", t322 == 322560 == line_form(hop322, {}, [3, 5, 2], -2, 1), f"{t322}")

    # a line at m = -1 on the unequal triangle, where b is not zero
    hopM, zzM = tris[0][1], tris[0][2]
    gdM = [2, 5, 1]
    lineM = centred_trace(3, hopM, zzM, gdM, [3 - x for x in gdM], 7)
    check("(G) a line gu = 3 - gd (m = -1) on the unequal triangle, b not zero: the seventh trace is the line form",
          lineM == line_form(hopM, zzM, gdM, 3, -1) and lineM != 0 and any(v != 0 for v in flow_per_im(hopM, zzM)),
          f"{lineM}")

    # the shear d -> d + s at every site keeps omega and s: the fifth and seventh traces stay, the ninth moves
    hopH, zzH = tris[0][1], tris[0][2]
    gdH, guH = [5, 3, 6], [1, 3, 2]
    sH = [(x + y) // 2 for x, y in zip(gdH, guH)]
    gdH2, guH2 = [x + v for x, v in zip(gdH, sH)], [y - v for y, v in zip(guH, sH)]
    before = [centred_trace(3, hopH, zzH, gdH, guH, m) for m in (5, 7, 9)]
    after = [centred_trace(3, hopH, zzH, gdH2, guH2, m) for m in (5, 7, 9)]
    check("(G) the shear d_l -> d_l + s_l on the unequal triangle: the fifth and seventh traces unchanged, the ninth "
          "changed", before[:2] == after[:2] and before[1] != 0 and before[2] != after[2], f"{before} -> {after}")

    # (S) a triangle with equal weights, without the full symmetry of the sites: fifth and seventh zero on equal heights, ninth not
    hopS = {(0, 1): 11 + 3j, (1, 2): 11 + 3j, (2, 0): 5 + 0j}
    zzS = {(0, 2): 21}
    bS = flow_per_im(hopS, zzS)
    gdS, guS = [4 + a for a in (3, -1, 2)], [4 - a for a in (3, -1, 2)]
    tS = [centred_trace(3, hopS, zzS, gdS, guS, m) for m in (5, 7, 9)]
    check("(S) t01 = t12 = 11+3i, t20 = 5, ZZ 21 on (0, 2): b = 0, Im = 330; at S = 4, arrows (3, -1, 2) the fifth and "
          "seventh traces are zero and the ninth is -53828247306240",
          bS == [0, 0, 0] and loop_im(hopS) == 330 and tS == [0, 0, -53828247306240], f"{bS}, {tS}")
    # the whole ninth trace of that triangle on equal heights, against its factored form at random integer points
    def ninth_S(S, d):
        q = (2436 * S * S + 4851 * d[0] ** 2 + 5846 * d[0] * d[1] + 2890 * d[0] * d[2] - 260 * d[1] ** 2
             + 5846 * d[1] * d[2] + 4851 * d[2] ** 2 + 2122512)
        return -6082560 * S * (d[0] - d[2]) * q
    ok, seen = True, False
    for _ in range(3):
        S = rnd.randint(2, 9)
        d = [rnd.randint(-S, S) for _ in range(3)]
        meas = centred_trace(3, hopS, zzS, [S + x for x in d], [S - x for x in d], 9)
        ok &= meas == ninth_S(S, d)
        seen |= meas != 0
    check("(S) on equal heights that triangle's ninth trace is -6082560*S*(d0 - d2)*Q(S, d), Q the indefinite quadratic "
          "of the proof's section 6, at three random points", ok and seen and tS[2] == ninth_S(4, [3, -1, 2]))

    # (E) equal heights on equal triangles, and two equal rate points
    eqs = [("unit, Delta = 0", unit, {}, 1, 0), ("unit, Delta = 1", unit, {(0, 1): 1, (1, 2): 1, (0, 2): 1}, 1, 1),
           ("unit, Delta = 3", unit, {(0, 1): 3, (1, 2): 3, (0, 2): 3}, 1, 3),
           ("t = 2+i, Delta = 2", tris[1][1], tris[1][2], 5, 2)]
    rows_e = []
    ok = True
    for name, hop, zz, r2, delta in eqs:
        S = rnd.randint(3, 8)
        arrows = [rnd.randint(-S + 1, S - 1) for _ in range(3)]
        while len(set(arrows)) < 3:
            arrows = [rnd.randint(-S + 1, S - 1) for _ in range(3)]
        gd = [S + a for a in arrows]
        gu = [S - a for a in arrows]
        fifth = centred_trace(3, hop, zz, gd, gu, 5)
        seventh = centred_trace(3, hop, zz, gd, gu, 7)
        ninth = centred_trace(3, hop, zz, gd, gu, 9)
        rows_e.append((hop, r2, delta, S, arrows, ninth))
        ok &= fifth == 0 and seventh == 0 and ninth == ninth_equal_heights(hop, r2, delta, S, arrows) and ninth != 0
    check("(E) equal triangles on equal heights: the fifth and seventh traces are zero and the ninth is "
          "36864*Im*(r^2 + 7*Delta^2)*S*V(d), not zero", ok)
    # the whole ninth trace on equal heights of the unit triangle, every degree, fitted exactly in (S, a0, a1, a2)
    monos9 = {d: list(itertools.combinations_with_replacement(range(4), d)) for d in (2, 4, 6, 8)}
    rows9 = {d: [] for d in monos9}
    rhs9 = {d: [] for d in monos9}
    for _ in range(100):
        S = rnd.randint(3, 12)
        arrows = [rnd.randint(-S + 1, S - 1) for _ in range(3)]
        parts = rate_parts(3, unit, {}, [S + a for a in arrows], [S - a for a in arrows], 9)
        x = [S] + arrows
        for d in monos9:
            rows9[d].append([np.prod([x[i] for i in mono], dtype=object) for mono in monos9[d]])
            rhs9[d].append(parts[d])
    fit9 = {d: dict(zip(monos9[d], solve_exact(rows9[d], rhs9[d]))) for d in (2, 4, 6)}
    want4 = {mono: Fraction(0) for mono in monos9[4]}
    # 36864*S*V(d) = -36864*S*(d0-d1)(d0-d2)(d1-d2) = -36864*S*(d0^2 d1 - d0^2 d2 - d0 d1^2 + d0 d2^2 + d1^2 d2 - d1 d2^2)
    for mono, sign in (((0, 1, 1, 2), 1), ((0, 1, 1, 3), -1), ((0, 1, 2, 2), -1), ((0, 1, 3, 3), 1), ((0, 2, 2, 3), 1), ((0, 2, 3, 3), -1)):
        want4[mono] = Fraction(-36864 * sign)
    check("(E) the ninth trace on equal heights of the unit triangle, fitted in degrees two, four and six: nothing in degree two or "
          "six, none in degree eight at any of the points, and degree four is exactly 36864*S*V(d)",
          all(v == 0 for v in fit9[2].values()) and all(v == 0 for v in fit9[6].values()) and fit9[4] == want4
          and all(v == 0 for v in rhs9[8]))
    ok = True
    for name, hop, zz, r2, delta in eqs:
        gd = [rnd.randint(1, 6) for _ in range(2)]
        gu = [rnd.randint(0, 6) for _ in range(2)]
        gd = [gd[0], gd[1], gd[0]]
        gu = [gu[0], gu[1], gu[0]]
        ok &= centred_trace(3, hop, zz, gd, gu, 7) == 0 and centred_trace(3, hop, zz, gd, gu, 9) == 0
    check("(E) control: equal triangles with two equal rate points (the site-permutation carrier) have seventh and ninth "
          "traces zero", ok)

    # (V) every constant moved by one is noticed
    unnoticed = []
    for key in C:
        moved = dict(C)
        moved[key] += 1
        hit = any(meas != r4_form(hop, gd, gu, moved) for hop, gd, gu, meas in rows_f) or \
            any(meas != triangle_sum_form(n, hop, gd, gu, moved) for n, hop, gd, gu, meas in rows_k) or \
            any(meas != (line_form(hop, zz, x, e, m, moved) if kind == "line" else vertical_form(hop, zz, x, e, moved))
                for kind, hop, zz, x, e, m, meas in rows_g) or \
            any(ninth != ninth_equal_heights(hop, r2, delta, S, arrows, moved) for hop, r2, delta, S, arrows, ninth in rows_e)
        if not hit:
            unnoticed.append(key)
    check(f"(V) each of the {len(C)} constants, moved by one, misses a row", not unnoticed,
          f"unnoticed: {unnoticed}" if unnoticed else "")

    print("ALL PASS" if not FAILS else f"{len(FAILS)} FAILED")
    return 0 if not FAILS else 1


if __name__ == "__main__":
    sys.exit(main())
