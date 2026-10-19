"""
The sodium nucleus: two mirrors in the quadrupolar relaxation of a spin 3/2
==========================================================================
experiments/THE_SODIUM_NUCLEUS_TWO_MIRRORS.md. One spin 3/2 (23Na), Lindblad generator with the five
rank-2 jumps T_{2,m}, m = -2..2, at rates J0, J1, J1, J2, J2 (the secular Redfield form of
quadrupolar relaxation at infinite temperature), a static quadrupolar splitting wQ·T20 and an offset
w·Iz. All exact, in sympy, with J0, J1, J2, wQ, w as symbols: no eigensolver, no tolerance.

Checks:
  - the generator is block diagonal in the coherence order p = m - m';
  - the closed rates at wQ = w = 0 are those of Souza et al., arXiv:1105.2539, Eqs. (15)-(27),
    with their constant set to C = 3 by this normalization of T_{2,m}: populations 0, 2CJ1, 2CJ2, 2C(J1+J2);
    central transition and triple quantum C(J1+J2); satellites C(J0+J1), C(J0+J1)+2CJ2;
    double quantum C(J0+J2), C(J0+J2)+2CJ1;
  - the dark mirror: the populations pair about -3(J1+J2), and the central transition and the triple
    quantum sit exactly at Re = -3(J1+J2), for every wQ and w;
  - the bright mirror: after removing the block frequency p·w, the satellite and double-quantum
    polynomials have real coefficients and root sum -6(J0+J1+J2), so their roots pair about
    -3(J0+J1+J2) under λ -> -conj(λ) - 2c, for every wQ and w;
  - T20 is diagonal, (3, -3, -3, 3)/sqrt(6) on m = 3/2, 1/2, -1/2, -3/2: one Z on the outer/inner bit.
  - why spin 3/2: the five Hermitian components of T_{2,m} are Dirac matrices (ZZ, ZX, ZY, XI, YI in the
    two-qubit labelling); their sixteen products diagonalize the generator at wQ = w = 0, class by class; S(X) = Gamma0·X^T gives S L S^-1 = -L + 2·L_J0 - 6(J1+J2)
    at every wQ and w, both mirrors in one map;
  - beyond spin 3/2 (I = 1 .. 7/2): T20² is scalar only at 3/2; J0 alone counts more zero-price cells
    than top-price cells;
  - no global centre beyond spin 3/2 at wQ = w = 0, for every real rate triple other than all zero, at each I = 1, 2, 5/2, ..., 6:
    a centre would sit at -s (the trace) and make every odd centred moment m_k = Tr (L + s)^k vanish;
    m3, m5, m7 are exact homogeneous polynomials (a rational similarity, block by block), and m7 is
    non-zero at every real common zero of m3 and m5 (resultant, subresultant, interval arithmetic on
    exact isolating intervals), the counts of common zeros, of positive ones and of ordered ones (J0 > J1 > J2) pinned;
    at 3/2 the same route's closed forms vanish exactly on the theorem's zero set (one way only: three vanishing moments would not prove a centre); the jump identities behind
    the offset argument (T_{2,-m} = ±T_{2,m}^†, [Iz, T_{2,m}] = m·T_{2,m}) checked at every spin;
  - off resonance (w != 0, wQ = 0) at every spin I >= 2: the diagonal of -L on |m><m'| as a closed polynomial in m
    and I(I+1), checked against the generator in every cell at I = 1 .. 6, the mean rate likewise; the traces of
    the blocks of order 2I, 2I-1, 2I-2 must be 1, 2, 3 times -s, three linear forms in J0, J1, J2 whose determinant
    I(I-1)(2I-3)²(2I-1)³(2I+5)(4I+3)/60 vanishes for I >= 1 only at 1 and 3/2; there the kernel is the ray
    (0, 2, 1) and the plane J0 = 0;
  - every offset and splitting, every spin I >= 2: the populations block W carries neither J0, wQ nor w,
    and the orders p, -p have conjugate polynomials, so a centre leaves the odd-multiplicity rates of W symmetric;
    nu3, nu5 of W exact in I and t = J2/J1 by interpolation (degree <= 4k + 1); nu3 at J1 = 0 and the
    resultant in t vanish at no spin I >= 2, so W simple forbids a centre; the route reaches 3/2 (identically 0)
    and the spin-1 rays (t = 1/8, 5/4); the repeated-rate ratios are finitely many (integer I: J2 = 0 is an
    irreducible chain; half-integer I to 61/2: a squarefree point mod 1000003), six real ones at I = 2;
  - a static splitting at w = 0, I = 1, 2, 5/2, ..., 6 (not 3/2): mu3, mu5, mu7, mu9 as exact polynomials in
    J0, J1, J2 and q = wQ/sqrt(6) (over Q(i), real traces), their common real zeros by exact elimination: none
    at I >= 2 other than all rates zero, six rays at I = 1, five of them with wQ real and a centre by the
    characteristic polynomial (nine eigenvalues: the four moments are the whole condition), two with
    non-negative rates, both J0 = 0; a detuned splitting and an offset w = 1 break both as complex maps while
    the decay rates stay palindromic; the same two rays again by hand from the block polynomials; at both,
    simple block eigenvalues, semisimple, the pairing crossing the coherence orders, just below the
    single-quantum exceptional point, 3 of 9 paired without the splitting;
  - the global centre: the trace fixes 2s = 3J0 + 6J1 + 6J2; the coefficients of P(-x - 2s) - P(x) have the
    lex Groebner basis {J0·wQ², J0·(9J1² + w²)(9J2² + 4w²)}, whose real zeros are the theorem; the three
    reflections at a vanishing density; and sixteen exact points, both outcomes, each by the identity at the
    point and by F158's end count dim ker L = dim ker(L + 2s) as a second route.

Import-inert; prints only, and exits with status 1 if any check fails.
"""

import sys
from collections import Counter

import sympy as sp

CLAIMS = []


def claim(label, ok, detail=""):
    CLAIMS.append(bool(ok))
    print(f"  [{'ok' if ok else 'FAIL'}] {label} {detail}")


def spin_ops(I):
    """Levels m = I .. -I, Iz and the rank-2 tensors T_{2,m} of a spin I."""
    ms = [I - k for k in range(int(2 * I) + 1)]
    d = len(ms)
    Iz = sp.diag(*ms)
    Ip = sp.zeros(d)
    for k in range(1, d):
        Ip[k - 1, k] = sp.sqrt(I * (I + 1) - ms[k] * (ms[k] + 1))
    Im = Ip.T
    T = {0: (3 * Iz * Iz - I * (I + 1) * sp.eye(d)) / sp.sqrt(6),
         1: -(Iz * Ip + Ip * Iz) / 2, -1: (Iz * Im + Im * Iz) / 2,
         2: Ip * Ip / 2, -2: Im * Im / 2}
    return ms, Iz, T


def generator(I, J0, J1, J2, wQ=0, w=0):
    ms, Iz, T = spin_ops(I)
    d = len(ms)
    Id = sp.eye(d)
    kron = sp.kronecker_product
    H = wQ * T[0] + w * Iz
    L = -sp.I * (kron(H, Id) - kron(Id, H.T))
    for g, A in [(J0, T[0]), (J1, T[1]), (J1, T[-1]), (J2, T[2]), (J2, T[-2])]:
        AdA = A.H * A
        L += g * (kron(A, A.conjugate()) - kron(AdA, Id) / 2 - kron(Id, AdA.T) / 2)
    return ms, T, L


def rational_blocks(I):
    """The three rate parts of the generator at wQ = w = 0, one block per coherence order, in rationals.

    The diagonal similarity D with D·I+·D^-1 = (c_k²) above the diagonal and D·I-·D^-1 = 1 below it makes
    every jump rational; rho -> D·rho·D^-1 is a similarity of the superoperator, so the spectrum is kept.
    """
    ms = [I - k for k in range(int(2 * I) + 1)]
    d = len(ms)
    Iz = sp.diag(*ms)
    Ip, Im = sp.zeros(d), sp.zeros(d)
    for k in range(1, d):
        Ip[k - 1, k] = I * (I + 1) - ms[k] * (ms[k] + 1)
        Im[k, k - 1] = 1
    T0 = 3 * Iz * Iz - I * (I + 1) * sp.eye(d)             # sqrt(6)·T20, so its rate carries 1/6
    T1, Tm1 = -(Iz * Ip + Ip * Iz) / 2, (Iz * Im + Im * Iz) / 2
    T2, Tm2 = Ip * Ip / 2, Im * Im / 2
    jumps = [(0, T0, T0, sp.Rational(1, 6)), (1, T1, -Tm1, 1), (1, Tm1, -T1, 1), (2, T2, Tm2, 1), (2, Tm2, T2, 1)]
    blocks = []
    for p in range(-(d - 1), d):
        cells = [(i, i + p) for i in range(d) if 0 <= i + p < d]
        index = {c: t for t, c in enumerate(cells)}
        parts = [sp.zeros(len(cells)) for _ in range(3)]
        for part, A, Adag, scale in jumps:
            AdA = Adag * A
            for (i, j), col in index.items():
                for (k, l), row in index.items():      # <k|A E_ij Adag - (AdA E_ij + E_ij AdA)/2|l>
                    v = A[k, i] * Adag[j, l] - (AdA[k, i] * int(j == l) + int(k == i) * AdA[j, l]) / 2
                    if v:
                        parts[part][row, col] += scale * v
        blocks.append(parts)
    return d * d, blocks


def odd_moments(I, syms):
    """m_k = Tr (L + s)^k, k = 3, 5, 7, at wQ = w = 0, as exact homogeneous polynomials in J0, J1, J2.

    Fitted on the principal lattice (i, j, 1), i + j <= k, unisolvent for degree <= k in two variables,
    homogenized, and checked exactly at points off the lattice, J2 = 0 among them.
    """
    from sympy.polys.matrices import DomainMatrix
    n, blocks = rational_blocks(I)

    def moments(J):
        mats = [sum((sp.Rational(J[t]) * parts[t] for t in range(3)), sp.zeros(parts[0].rows)) for parts in blocks]
        s = -sum(M.trace() for M in mats) / n
        out = {}
        for M in mats:
            M = M + s * sp.eye(M.rows)
            Q = M
            for k in range(2, 8):
                Q = Q * M
                if k % 2:
                    out[k] = out.get(k, 0) + Q.trace()
        return out

    J0, J1, J2 = syms
    polys = {}
    for k in (3, 5, 7):
        monos = [(i, j) for i in range(k + 1) for j in range(k + 1 - i)]
        A = DomainMatrix([[sp.QQ(p[0] ** a * p[1] ** b) for a, b in monos] for p in monos], (len(monos),) * 2, sp.QQ)
        b = DomainMatrix([[sp.QQ.from_sympy(sp.Rational(moments((p[0], p[1], 1))[k]))] for p in monos], (len(monos), 1), sp.QQ)
        c = A.lu_solve(b).to_Matrix()
        polys[k] = sp.expand(sum(c[t] * J0 ** a * J1 ** b_ * J2 ** (k - a - b_) for t, (a, b_) in enumerate(monos)))
        for test in [(3, -2, 5), (-7, 4, 2), (1, 1, 0), (2, 0, 0), (0, 3, 0)]:
            assert moments(test)[k] == polys[k].subs({J0: test[0], J1: test[1], J2: test[2]}), (I, k, test)
    return polys


def common_real_zeros(m, syms):
    """Every real projective common zero of m3, m5: the signs of J0, J1 (at J2 = 1), whether m7 != 0 there,
    and whether J0 > J1 > J2 = 1 holds at it.

    J2 = 0 must carry none (resultant of the two binary forms non-zero, m3(1, 0, 0) non-zero). On J2 = 1 the
    common zeros sit at the real roots of the resultant in J1 (the leading coefficients in J1 are constants, so
    it specializes), and since the degree-1 subresultant's leading coefficient is checked non-zero at each root,
    there J1 = -s10/s11. Signs and m7 are read by interval arithmetic on exact isolating intervals.
    """
    import mpmath
    iv = mpmath.iv
    iv.dps = 60
    J0, J1, J2 = syms
    a0, b0 = m[3].subs({J2: 0, J1: 1}), m[5].subs({J2: 0, J1: 1})
    assert sp.resultant(a0, b0, J0) != 0 and m[3].subs({J0: 1, J1: 0, J2: 0}) != 0, "a common zero on J2 = 0"
    a, b, c = (sp.expand(m[k].subs(J2, 1)) for k in (3, 5, 7))
    assert sp.Poly(a, J1).LC().is_Rational and sp.Poly(b, J1).LC().is_Rational
    res = sp.Poly(sp.resultant(a, b, J1), J0)
    assert res.degree() > 0, "m3 and m5 share a factor: the resultant vanishes"
    S1 = sp.Poly(sp.subresultants(a, b, J1)[-2], J1)
    assert S1.degree() == 1
    s11, s10 = S1.coeff_monomial(J1), S1.coeff_monomial(1)

    def ev(e, *args):
        tot = iv.mpf(0)
        for mono, co in sp.Poly(e, *(J0, J1)[:len(args)]).terms():
            p, q = sp.fraction(co)
            t = iv.mpf(int(p)) / iv.mpf(int(q))
            for v, e_ in zip(args, mono):
                t = t * v ** e_
            tot = tot + t
        return tot

    def sign(v):
        return "+" if v.a > 0 else "-" if v.b < 0 else "?"

    out = []
    for (lo, hi), _ in res.intervals(eps=sp.Rational(1, 10**40)):
        p_lo, q_lo = sp.fraction(lo)
        p_hi, q_hi = sp.fraction(hi)
        X = iv.mpf([(iv.mpf(int(p_lo)) / iv.mpf(int(q_lo))).a, (iv.mpf(int(p_hi)) / iv.mpf(int(q_hi))).b])
        den = ev(s11, X)
        assert not den.a <= 0 <= den.b, "degree-1 subresultant degenerate at a root"
        x1 = ev(-s10, X) / den
        m7 = ev(c, X, x1)
        assert "?" not in (sign(X), sign(x1)), "a sign not decided"
        out.append((sign(X), sign(x1), not m7.a <= 0 <= m7.b, X.a > x1.b and x1.a > 1))
    return out


def split_moments(I, syms, ks=(3, 5, 7, 9)):
    """mu_k = Tr (L + s)^k with the static splitting wQ·T20 at offset w = 0, exact over Q(i), as polynomials in
    J0, J1, J2 and q = wQ/sqrt(6), homogeneous of degree k.

    The splitting is diagonal and commutes with the rational similarity, so it adds -i·q·(t_a - t_b) to the
    cell (a, b), t_m = 3m² - I(I+1); real and imaginary parts are carried as two rational matrices. Fitted on
    the principal lattice (i, j, l, J2 = 1), i + j + l <= k, checked off the lattice, and every trace is real.
    """
    from sympy.polys.matrices import DomainMatrix
    n, blocks = rational_blocks(I)
    ms = [I - k for k in range(int(2 * I) + 1)]
    d = len(ms)
    t = [3 * m * m - I * (I + 1) for m in ms]
    QQ = sp.QQ
    dm = []
    for p, parts in zip(range(-(d - 1), d), blocks):
        cells = [(i, i + p) for i in range(d) if 0 <= i + p < d]
        size = len(cells)
        conv = [DomainMatrix([[QQ.from_sympy(P[r, c]) for c in range(size)] for r in range(size)], (size, size), QQ)
                for P in parts]
        h = DomainMatrix.diag([QQ.from_sympy(-(t[i] - t[j])) for i, j in cells], QQ, (size, size))
        dm.append((conv, h, size))

    def tr(M):
        return sum((M[i, i].element for i in range(M.shape[0])), QQ(0))

    def moments(J, q):
        A_s = [conv[0] * QQ.from_sympy(sp.Rational(J[0])) + conv[1] * QQ.from_sympy(sp.Rational(J[1]))
               + conv[2] * QQ.from_sympy(sp.Rational(J[2])) for conv, _, _ in dm]
        s = -sum((tr(A) for A in A_s), QQ(0)) / n
        out = {k: [QQ(0), QQ(0)] for k in ks}
        for A, (_, h, size) in zip(A_s, dm):
            A = A + DomainMatrix.eye(size, QQ) * s
            B = h * QQ.from_sympy(sp.Rational(q))
            X, Y = A, B
            for k in range(2, max(ks) + 1):
                X, Y = X * A - Y * B, X * B + Y * A
                if k in out:
                    out[k][0] += tr(X)
                    out[k][1] += tr(Y)
        return out

    J0, J1, J2, q = syms
    polys = {}
    for k in ks:
        monos = [(a, b, c) for a in range(k + 1) for b in range(k + 1 - a) for c in range(k + 1 - a - b)]
        A = DomainMatrix([[QQ(p[0] ** a * p[1] ** b * p[2] ** c) for a, b, c in monos] for p in monos], (len(monos),) * 2, QQ)
        vals = [moments((p[0], p[1], 1), p[2])[k] for p in monos]
        assert all(v[1] == 0 for v in vals), "an imaginary trace"
        sol = A.lu_solve(DomainMatrix([[v[0]] for v in vals], (len(monos), 1), QQ)).to_Matrix()
        polys[k] = sp.expand(sum(sol[r] * J0 ** a * J1 ** b * q ** c * J2 ** (k - a - b - c) for r, (a, b, c) in enumerate(monos)))
        for test in [(3, -2, 5, 7), (-7, 4, 2, -3), (1, 1, 0, 2), (2, 0, 0, 1), (0, 3, 0, 0)]:
            got = moments(test[:3], test[3])[k]
            assert got[1] == 0 and sp.Rational(got[0]) == polys[k].subs(dict(zip(syms, test))), (I, k, test)
    return polys


def split_centre_points(P, syms):
    """Exact elimination of mu3 = mu5 = mu7 = mu9 = 0 at w = 0 with a splitting.

    In u = q², mu3 = c + l·u. Chart J2 = 1, branch l != 0: u = -c/l, g_k = l^deg·mu_k(u = -c/l) with the spurious
    powers of l divided out; candidate J0 are the roots of gcd(Res_J1(g5, g7), Res_J1(g5, g9)). A linear factor is
    solved exactly; a non-linear one must have its common zeros with g5, g7, g9 on l = 0 (the ideal with 1 - z·l
    added is the unit ideal). Branch l = 0 (then c = 0 too): mu5, mu7, mu9 must share a root in u, here a gcd of resultants.
    The line J2 = 0 by Groebner bases. Returns the rational points (J0, J1, u) of the chart and whether every
    other case is empty.
    """
    J0, J1, J2, q = syms
    u = sp.Symbol("u")

    def in_u(e):
        return sp.expand(sp.expand(e).replace(lambda w_: w_.is_Pow and w_.base == q, lambda w_: u ** (w_.exp // 2)))

    F = {k: in_u(P[k].subs(J2, 1)) for k in P}
    assert all(not f.has(q) for f in F.values()), "an odd power of q"
    p3 = sp.Poly(F[3], u)
    assert p3.degree() == 1
    c, l = p3.coeff_monomial(1), p3.coeff_monomial(u)
    Lp = sp.Poly(l, J0, J1)
    assert Lp.degree(J1) == 1 and Lp.coeff_monomial(J1).is_Rational, "l is not linear in J1 with a constant coefficient"

    def sub(f):
        pf = sp.Poly(f, u)
        g = sp.Poly(sum(co * (-c) ** e[0] * l ** (pf.degree() - e[0]) for e, co in pf.terms()), J0, J1)
        while g.rem(Lp).is_zero:
            g = g.quo(Lp)
        return g.as_expr()

    g5, g7, g9 = sub(F[5]), sub(F[7]), sub(F[9])
    r57, r59 = sp.Poly(sp.resultant(g5, g7, J1), J0), sp.Poly(sp.resultant(g5, g9, J1), J0)
    assert r57.degree() > 0 and r59.degree() > 0
    points, rest_empty = [], True
    for f, _ in sp.factor_list(sp.gcd(r57, r59).as_expr())[1]:
        if not f.has(J0):
            continue
        if sp.Poly(f, J0).degree() == 1:
            x0 = sp.solve(f, J0)[0]
            h = sp.gcd(sp.gcd(sp.Poly(g5.subs(J0, x0), J1), sp.Poly(g7.subs(J0, x0), J1)), sp.Poly(g9.subs(J0, x0), J1))
            for x1 in sp.Poly(h, J1).all_roots():
                if x1.is_real and l.subs({J0: x0, J1: x1}) != 0:
                    uv = -c.subs({J0: x0, J1: x1}) / l.subs({J0: x0, J1: x1})
                    assert all(sp.simplify(F[k].subs({J0: x0, J1: x1, u: uv})) == 0 for k in (3, 5, 7, 9))
                    points.append((x0, x1, uv))
        else:
            z = sp.Symbol("z")              # Rabinowitsch: no common zero with l != 0
            rest_empty &= sp.groebner([f, g5, g7, g9, 1 - z * l], z, J1, J0, order="grevlex").exprs == [1]
    x1_of = sp.solve(l, J1)[0]
    M = [sp.expand(F[k].subs(J1, x1_of)) for k in (5, 7, 9)]
    rest_empty &= sp.gcd(sp.gcd(sp.Poly(sp.expand(c.subs(J1, x1_of)), J0), sp.Poly(sp.resultant(M[0], M[1], u), J0)),
                         sp.Poly(sp.resultant(M[0], M[2], u), J0)).degree() == 0
    line = [in_u(P[k].subs({J2: 0, J1: 1})) for k in (3, 5, 7, 9)]
    rest_empty &= sp.groebner(line, u, J0, order="lex").exprs == [1]
    corner = [sp.Poly(in_u(P[k].subs({J0: 1, J1: 0, J2: 0})), u) for k in (3, 5, 7, 9)]
    g = corner[0]
    for f in corner[1:]:
        g = sp.gcd(g, f)
    rest_empty &= g.degree() == 0
    return points, rest_empty


def build():
    J0, J1, J2, wQ, w = sp.symbols("J0 J1 J2 wQ w", real=True)
    ms, T, L = generator(sp.Rational(3, 2), J0, J1, J2, wQ, w)
    return (J0, J1, J2, wQ, w), ms, T, sp.simplify(L)


def main():
    R = sp.Rational
    (J0, J1, J2, wQ, w), ms, T, L = build()
    x, y = sp.symbols("x y")
    d = 4
    print("the jump T20 on the four levels m = 3/2, 1/2, -1/2, -3/2:", [sp.nsimplify(T[0][k, k] * sp.sqrt(6)) for k in range(d)], "/ sqrt(6)")
    claim("T20 is one Z on the outer/inner bit (equal on |m| = 3/2, opposite on |m| = 1/2)",
          T[0].is_diagonal() and T[0][0, 0] == T[0][3, 3] == -T[0][1, 1] == -T[0][2, 2])

    blocks = {}
    for p in range(-3, 4):
        idx = [a * d + b for a in range(d) for b in range(d) if ms[a] - ms[b] == p]
        rest = [k for k in range(d * d) if k not in idx]
        leak = L.extract(idx, rest)
        claim(f"block p = {p:+d} is invariant (no entry to other orders)", leak.is_zero_matrix)
        blocks[p] = sp.expand(L.extract(idx, idx).charpoly(x).as_expr())

    C = 3
    c1, c2 = C * (J1 + J2), C * (J0 + J1 + J2)
    print("\nclosed rates at wQ = w = 0, against Souza et al. Eqs. (15)-(27) with C = 3")
    expect = {
        0: [0, 2 * C * J1, 2 * C * J2, 2 * C * (J1 + J2)],
        1: [C * (J1 + J2), C * (J0 + J1), C * (J0 + J1) + 2 * C * J2],
        2: [C * (J0 + J2), C * (J0 + J2) + 2 * C * J1],
        3: [C * (J1 + J2)],
    }
    for p, rates in expect.items():
        target = sp.expand(sp.prod([x + r for r in rates]))
        got = sp.expand(blocks[p].subs({wQ: 0, w: 0}))
        claim(f"p = {p}: decay rates {rates}", sp.expand(got - target) == 0)

    print("\nthe dark mirror about -3(J1 + J2), for every wQ and w")
    pops = sp.Poly(blocks[0], x)
    claim("populations: roots 0, -6J1, -6J2, -6(J1+J2) for every wQ, w (they pair about -c1)",
          sp.expand(blocks[0] - sp.expand(x * (x + 6 * J1) * (x + 6 * J2) * (x + 6 * J1 + 6 * J2))) == 0)
    claim("central transition: (x + c1 + i·w) divides the p = 1 polynomial",
          sp.rem(sp.Poly(blocks[1], x), sp.Poly(x + c1 + sp.I * w, x)).is_zero)
    claim("triple quantum: the p = 3 polynomial is x + c1 + 3i·w", sp.expand(blocks[3] - (x + c1 + 3 * sp.I * w)) == 0)

    print("\nthe bright mirror about -3(J0 + J1 + J2), for every wQ and w")
    sat = sp.quo(sp.Poly(blocks[1], x), sp.Poly(x + c1 + sp.I * w, x)).as_expr()
    for p, poly in [(1, sat), (2, blocks[2])]:
        shifted = sp.Poly(sp.expand(poly.subs(x, y - p * sp.I * w)), y)
        coeffs = shifted.all_coeffs()
        real = all(sp.simplify(sp.im(sp.expand(c))) == 0 for c in coeffs)
        claim(f"p = {p}: after the shift x = y - {p}i·w the quadratic has real coefficients", real,
              f"({[sp.factor(c) for c in coeffs]})")
        claim(f"p = {p}: root sum is -2·3(J0 + J1 + J2)", sp.expand(-coeffs[1] / coeffs[0] + 2 * c2) == 0)

    claim("the static splitting enters both quadratics only as +6·wQ² in the constant term",
          all(sp.expand(sp.Poly(sp.expand(poly.subs(x, y - p * sp.I * w)), y).all_coeffs()[2]
                        - sp.Poly(sp.expand(poly.subs(x, y - p * sp.I * w)), y).all_coeffs()[2].subs(wQ, 0) - 6 * wQ ** 2) == 0
              for p, poly in [(1, sat), (2, blocks[2])]))
    claim("the negative orders are the conjugates: P_-p(x) = conj(P_p)(x) coefficient by coefficient",
          all(sp.expand(sp.Poly(blocks[-p], x).as_expr()
                        - sum(sp.conjugate(c) * x ** k for k, c in enumerate(reversed(sp.Poly(blocks[p], x).all_coeffs())))) == 0
              for p in (1, 2, 3)))
    # the centres read off the generator's own polynomials, not from the closed forms above
    c1_read = sp.expand(blocks[3] - x).subs(w, 0)                               # p = 3: x + c1 + 3iw
    c2_read = sp.Poly(sp.expand(blocks[2].subs(x, y - 2 * sp.I * w)), y).all_coeffs()[1] / 2
    claim("the centres read off the generator differ by exactly 3·J0, so they coincide iff J0 = 0",
          sp.expand(c2_read - c1_read - 3 * J0) == 0, f"(c1 = {sp.factor(c1_read)}, c2 = {sp.factor(c2_read)})")
    D0 = lambda r, c: L.subs({J1: 0, J2: 0, wQ: 0, w: 0})[r * d + c, r * d + c]
    claim("the J0 dissipator alone charges a bright cell 3·J0 and a dark cell 0 (satellite, double quantum; "
          "central, triple quantum)", [D0(0, 1), D0(0, 2), D0(1, 2), D0(0, 3)] == [-3 * J0, -3 * J0, 0, 0])
    for p, poly, Jp in [(1, sat, J2), (2, blocks[2], J1)]:
        co = sp.Poly(sp.expand(poly.subs(x, y - p * sp.I * w)), y).all_coeffs()
        disc = sp.expand(co[1] ** 2 / 4 - co[2])
        claim(f"p = {p}: the pair is real while 6wQ² < 9·{Jp}², a conjugate pair beyond (discriminant/4 = 9{Jp}² - 6wQ²)",
              sp.expand(disc - (9 * Jp ** 2 - 6 * wQ ** 2)) == 0)
    pop_rates = sorted((sp.expand(-r) for r in sp.roots(sp.Poly(blocks[0], x))), key=str)
    claim("the triple-quantum rate, read off the p = 3 block, is the mean of the two longitudinal rates "
          "6J1, 6J2 read off the population block",
          {6 * J1, 6 * J2} <= set(pop_rates) and sp.expand(sp.expand(blocks[3] - x).subs(w, 0) - (6 * J1 + 6 * J2) / 2) == 0
          and sp.diag(*ms).trace() == 0 and (sp.diag(*ms) * T[0]).trace() == 0)   # Iz has no weight on 1 or on T20

    print("\nthe global centre: one centre for the whole spectrum iff J0 = 0, or wQ = w = 0 and J1·J2 = 0")
    claim("Tr L = -24·J0 - 48·J1 - 48·J2 at every wQ and w, so 2s = 3J0 + 6J1 + 6J2",
          sp.expand(L.trace() + 24 * J0 + 48 * J1 + 48 * J2) == 0)
    s2 = 3 * J0 + 6 * J1 + 6 * J2
    P = sp.expand(sp.prod([blocks[p] for p in range(-3, 4)]))
    coeffs = [c for c in sp.Poly(sp.expand(P.subs(x, -x - s2) - P), x).all_coeffs() if sp.expand(c) != 0]
    G = sp.groebner(coeffs, J0, J1, J2, wQ, w, order="lex")
    target = sp.groebner([J0 * wQ ** 2, J0 * (9 * J1 ** 2 + w ** 2) * (9 * J2 ** 2 + 4 * w ** 2)], J0, J1, J2, wQ, w, order="lex")
    claim("P(-x - 2s) - P(x): its coefficients generate the ideal with lex Groebner basis "
          "{J0·wQ², J0·(9J1² + w²)(9J2² + 4w²)} (real zeros: J0 = 0, or wQ = w = 0 and J1·J2 = 0)",
          len(coeffs) == 14 and list(G.exprs) == list(target.exprs), f"({len(coeffs)} coefficients)")
    Dpoly = sp.Poly(sp.expand(P.subs(x, -x - s2) - P), x)
    claim("the difference has degree 13 and leading coefficient 144·J0·wQ²",
          Dpoly.degree() == 13 and sp.expand(Dpoly.LC() - 144 * J0 * wQ ** 2) == 0)
    # the closed rates at wQ = w = 0, taken from the gated list above, sixteen with multiplicity
    rates = (expect[0] + [expect[1][0]] * 2 + expect[1][1:] * 2 + expect[2] * 2 + expect[3] * 2)
    claim("the closed rates count sixteen and sum to 16·s", len(rates) == 16 and sp.expand(sum(rates) - 8 * s2) == 0)
    for z in (J0, J1, J2):
        R0, S0 = [sp.expand(sp.sympify(r).subs(z, 0)) for r in rates], sp.expand(s2.subs(z, 0))
        claim(f"{z} = 0, wQ = w = 0: the sixteen rates are invariant under r -> {S0} - r",
              Counter(R0) == Counter(sp.expand(S0 - r) for r in R0))

    cases = [("J0 alone", {J0: 1, J1: 0, J2: 0, wQ: 0, w: 0}, True),
             ("J1 alone", {J0: 0, J1: 1, J2: 0, wQ: 0, w: 0}, True),
             ("J1 and J2, no J0", {J0: 0, J1: R(2, 5), J2: R(7, 3), wQ: 0, w: 0}, True),
             ("J1 and J2, no J0, with wQ and w", {J0: 0, J1: R(2, 5), J2: R(7, 3), wQ: 1, w: 3}, True),
             ("extreme narrowing J0 = J1 = J2", {J0: 1, J1: 1, J2: 1, wQ: 0, w: 0}, False),
             ("J0 with J1, no J2", {J0: 1, J1: R(1, 3), J2: 0, wQ: 0, w: 0}, True),
             ("J0 with J2, no J1", {J0: R(3, 4), J1: 0, J2: R(5, 2), wQ: 0, w: 0}, True),
             ("J0 with J2, no J1, with wQ", {J0: R(3, 4), J1: 0, J2: R(5, 2), wQ: 1, w: 0}, False),
             ("J0 with J2, no J1, with w", {J0: R(3, 4), J1: 0, J2: R(5, 2), wQ: 0, w: 1}, False),
             ("J0 with J1, no J2, with wQ", {J0: 1, J1: R(1, 3), J2: 0, wQ: R(1, 2), w: 0}, False),
             ("J0 with J1, no J2, with w", {J0: 1, J1: R(1, 3), J2: 0, wQ: 0, w: 3}, False),
             ("J0 alone with a static splitting", {J0: 1, J1: 0, J2: 0, wQ: 1, w: 0}, False),
             ("J0 alone, w = sqrt(6)·wQ: a real eigenvalue at -2s, twice, against four zeros",
              {J0: 1, J1: 0, J2: 0, wQ: sp.sqrt(6), w: 6}, False),
             ("generic, all three positive", {J0: 1, J1: R(37, 100), J2: R(81, 100), wQ: 0, w: 0}, False),
             ("all three positive, other ratios", {J0: R(15, 2), J1: 21, J2: 2, wQ: 0, w: 0}, False),
             ("all three positive, small J0", {J0: R(1, 5), J1: R(16, 5), J2: R(8, 3), wQ: 0, w: 0}, False)]
    for name, vals, expected in cases:
        zero = all(sp.expand(sp.radsimp(c.subs(vals))) == 0 for c in coeffs)
        Lp = L.subs(vals).applyfunc(lambda z: sp.nsimplify(sp.expand(z)))
        two_s = s2.subs(vals)                                                   # s2 holds 2s
        ends = (Lp.rank(simplify=True) == (Lp + two_s * sp.eye(d * d)).rank(simplify=True))  # F158: the two ends
        claim(f"{name}: global centre {'exists' if expected else 'does not exist'} "
              "(the polynomial identity at the point, and F158's end count)",
              zero == ends == expected)

    print("\nwhy spin 3/2: the five quadrupole components are Dirac matrices")
    G = [T[0], T[1] + T[1].H, sp.I * (T[1] - T[1].H), T[2] + T[2].H, sp.I * (T[2] - T[2].H)]
    sq = [sp.simplify(g * g) for g in G]
    claim("the five Hermitian components of T_{2,m} anticommute pairwise and square to 3/2, 3, 3, 3, 3 times 1",
          all((G[a] * G[b] + G[b] * G[a]).applyfunc(sp.simplify).is_zero_matrix for a in range(5) for b in range(a + 1, 5))
          and [sq[k] - v * sp.eye(d) == sp.zeros(d) for k, v in enumerate([R(3, 2), 3, 3, 3, 3])] == [True] * 5)
    # the sixteen products of at most two distinct Gammas; a jump Gamma_a charges every product it anticommutes
    # with the price 2·(rate)·Gamma_a², which is 3J0 for T20 and 3J1, 3J1, 3J2, 3J2 for the four others
    from itertools import combinations
    price = [3 * J0, 3 * J1, 3 * J1, 3 * J2, 3 * J2]
    L0 = L.subs({wQ: 0, w: 0})
    ok = True
    for k in range(0, 3):
        for S_ in combinations(range(5), k):
            B = sp.prod([G[a] for a in S_], start=sp.eye(d)) if S_ else sp.eye(d)
            vecB = B.reshape(d * d, 1)
            rate = sum(price[c] for c in range(5) if (G[c] * B + B * G[c]).applyfunc(sp.simplify).is_zero_matrix)
            ok &= (L0 * vecB + rate * vecB).applyfunc(sp.expand).is_zero_matrix
    claim("at wQ = w = 0 the 16 Clifford products (1, Gamma_a, Gamma_a·Gamma_b) diagonalize the generator; "
          "each decays at the summed price of the jumps it anticommutes with", ok)
    Xp, Yp, Zp, Ip_ = sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.diag(1, -1), sp.eye(2)
    strings = [sp.kronecker_product(*pq) for pq in [(Zp, Zp), (Zp, Xp), (Zp, Yp), (Xp, Ip_), (Yp, Ip_)]]

    def proportional(A, B):
        k = next((r, c) for r in range(d) for c in range(d) if B[r, c] != 0)
        return (A - sp.simplify(A[k] / B[k]) * B).applyfunc(sp.simplify).is_zero_matrix and A[k] != 0

    claim("in the labelling 3/2, 1/2, -1/2, -3/2 = 00, 01, 10, 11 the five Gammas are ZZ, ZX, ZY, XI, YI up to scalars",
          all(proportional(G[a], strings[a]) for a in range(5)))
    claim("and T20 = (3/sqrt6)·ZZ, Iz = ZI + IZ/2 there (the Hamiltonian as strings)",
          (T[0] - 3 / sp.sqrt(6) * strings[0]).applyfunc(sp.simplify).is_zero_matrix
          and sp.diag(*ms) == sp.kronecker_product(Zp, Ip_) + sp.kronecker_product(Ip_, Zp) / 2)
    blindcell = lambda r, c: abs(ms[r]) == abs(ms[c])
    table = {(): "population", (0,): "population", (1, 2): "population", (3, 4): "population",
             (1, 3): "dark", (1, 4): "dark", (2, 3): "dark", (2, 4): "dark",
             (1,): "satellite", (2,): "satellite", (0, 1): "satellite", (0, 2): "satellite",
             (3,): "double", (4,): "double", (0, 3): "double", (0, 4): "double"}
    ok = True
    for S_, kind in table.items():
        B = sp.prod([G[a] for a in S_], start=sp.eye(d)) if S_ else sp.eye(d)
        cells = [(r, c) for r in range(d) for c in range(d) if sp.simplify(B[r, c]) != 0]
        orders = {abs(ms[r] - ms[c]) for r, c in cells}
        blind = {blindcell(r, c) for r, c in cells}
        want = {"population": ({0}, {True}), "dark": ({1, 3}, {True}),
                "satellite": ({1}, {False}), "double": ({2}, {False})}[kind]
        ok &= (orders <= want[0]) and blind == want[1]
    claim("the sixteen products sit in their classes: four populations, four central-transition/triple-quantum, "
          "four satellite and four double-quantum operators", ok)
    G0 = T[0] * sp.sqrt(6) / 3
    kron = sp.kronecker_product
    perm = sp.zeros(d * d)
    for a_ in range(d):
        for b_ in range(d):
            perm[a_ * d + b_, b_ * d + a_] = 1
    S = kron(G0, sp.eye(d)) * perm                       # X -> Gamma0 · X^T
    LJ0 = sp.expand(L.diff(J0) * J0)
    claim("one map for both mirrors: S(X) = Gamma0·X^T satisfies S L S^-1 = -L + 2·L_J0 - 6(J1 + J2) at every wQ, w",
          (S * L * S.inv() + L - 2 * LJ0 + 6 * (J1 + J2) * sp.eye(d * d)).applyfunc(sp.expand).is_zero_matrix)
    S2 = kron(sp.eye(d), G0.T) * perm                    # X -> X^T · Gamma0, the other order
    claim("so does X -> X^T·Gamma0 (F118's two orders, Pi_Y = F·rho^T and Pi_Z = rho^T·F, with Gamma0 for F)",
          (S2 * L * S2.inv() + L - 2 * sp.expand(L.diff(J0) * J0) + 6 * (J1 + J2) * sp.eye(d * d)).applyfunc(sp.expand).is_zero_matrix)
    Pb = -LJ0 / (3 * J0)
    bright = sp.diag(*[0 if blindcell(r, c) else 1 for r in range(d) for c in range(d)])
    claim("the J0 dissipator is -3J0 times the projector on the bright cells, which commutes with L at every wQ, w "
          "and with S, so S maps each sector to itself",
          (Pb - bright).applyfunc(sp.expand).is_zero_matrix and (Pb * L - L * Pb).applyfunc(sp.expand).is_zero_matrix
          and (S * Pb - Pb * S).applyfunc(sp.expand).is_zero_matrix)

    print("\nbeyond spin 3/2 (wQ = w = 0, exact)")
    spins = [R(1), R(3, 2), R(2), R(5, 2), R(3), R(7, 2)]
    claim("T20^2 is a multiple of 1 only at I = 3/2 among I = 1 .. 7/2",
          [I_ for I_ in spins if (lambda t: t * t == (t[0, 0] ** 2) * sp.eye(t.shape[0]))(spin_ops(I_)[2][0])] == [R(3, 2)])
    for I_ in spins:
        ms_, _, T_ = spin_ops(I_)
        t = [3 * m * m for m in ms_]
        prices = Counter(sp.Rational(1, 2) * (ta - tb) ** 2 for ta in t for tb in t)
        top = max(prices)
        claim(f"I = {I_}, J0 alone: {prices[0]} cells at price 0 against {prices[top]} at the top price, "
              f"so {'a palindrome' if prices[0] == prices[top] else 'no palindrome'}",
              (prices[0] == prices[top]) == (I_ == R(3, 2))
              and (I_ != R(3, 2) or sorted(prices.values()) == [8, 8]))

    ms5, _, L5 = generator(R(5, 2), 0, R(2, 5), R(7, 3))
    blind = sp.diag(*[1 if abs(a_) == abs(b_) else 0 for a_ in ms5 for b_ in ms5])
    claim("at I = 5/2 the cells blind to J0 (|m| = |m'|) are not an invariant set: J1, J2 mix them with the rest",
          not ((sp.eye(len(ms5) ** 2) - blind) * L5 * blind).applyfunc(lambda z: sp.nsimplify(sp.expand(z))).is_zero_matrix)

    print("\nno global centre beyond spin 3/2 (wQ = w = 0, every real rate triple other than all zero, exact)")
    K0, K1, K2 = sp.symbols("J0 J1 J2", real=True)
    m = {}
    for I_ in [R(3, 2)] + [R(k, 2) for k in range(2, 13) if k != 3]:
        m[I_] = odd_moments(I_, (K0, K1, K2))
    claim("I = 3/2: m3 = 0, m5 = 14580·J0·J1²·J2², m7 = 229635·J0·J1²·J2²·(J0² + 2J1² + 2J2²), zero exactly "
          "where the theorem puts a centre (the route does not say no where the theorem says yes)",
          m[R(3, 2)][3] == 0 and sp.expand(m[R(3, 2)][5] - 14580 * K0 * K1**2 * K2**2) == 0
          and sp.expand(m[R(3, 2)][7] - 229635 * K0 * K1**2 * K2**2 * (K0**2 + 2 * K1**2 + 2 * K2**2)) == 0)
    ordered_want = {R(4): 1, R(9, 2): 1, R(5): 1, R(11, 2): 2, R(6): 2}
    want = {R(1): (5, 0), R(4): (11, 6), R(9, 2): (11, 6), R(5): (11, 6), R(11, 2): (11, 5), R(6): (11, 5)}
    for I_ in m:
        _, Iz_, T_ = spin_ops(I_)
        claim(f"I = {I_}: every jump comes with its adjoint, T_(2,-m) = (-1)^m T_(2,m)^dagger, and [Iz, T_(2,m)] = m T_(2,m)",
              all((T_[-q] - (-1) ** q * T_[q].H).applyfunc(sp.simplify).is_zero_matrix
                  and (Iz_ * T_[q] - T_[q] * Iz_ - q * T_[q]).applyfunc(sp.simplify).is_zero_matrix for q in (-2, -1, 0, 1, 2)))
    for I_ in [I_ for I_ in m if I_ != R(3, 2)]:
        zeros = common_real_zeros(m[I_], (K0, K1, K2))
        positive = sum(1 for z in zeros if z[0] == "+" and z[1] == "+")
        ordered = sum(1 for z in zeros if z[3])
        claim(f"I = {I_}: m3, m5 share {len(zeros)} real zeros on J2 = 1 and none on J2 = 0, {positive} of them "
              f"with J0, J1 > 0 ({ordered} with J0 > J1 > J2); m7 is non-zero at every one, so no rates other than "
              f"all zero give a centre",
              all(z[2] for z in zeros) and (len(zeros), positive) == want.get(I_, (3, 0))
              and ordered == ordered_want.get(I_, 0))

    print("\nwith an offset w != 0, every spin I >= 2 (wQ = 0, every real rate triple other than all zero, exact)")
    Is, ms_ = sp.symbols("I m")
    X_ = Is * (Is + 1)
    # diagonal of the sum of T^dag·T over each channel's jumps at level m, and the diagonal of T20
    A0 = (3 * ms_**2 - X_) ** 2 / 6
    A1 = ((2 * ms_ + 1) ** 2 * (X_ - ms_ * (ms_ + 1)) + (2 * ms_ - 1) ** 2 * (X_ - ms_ * (ms_ - 1))) / 4
    A2 = ((X_ - ms_ * (ms_ + 1)) * (X_ - (ms_ + 1) * (ms_ + 2)) + (X_ - ms_ * (ms_ - 1)) * (X_ - (ms_ - 1) * (ms_ - 2))) / 4
    t20 = (3 * ms_**2 - X_) / sp.sqrt(6)
    gm = K0 * A0 + K1 * A1 + K2 * A2

    def cell_rate(a, b):
        """The diagonal of -L on |a><b| at wQ = w = 0, closed in I: (g_a + g_b)/2 - J0·t_a·t_b."""
        return sp.expand(gm.subs(ms_, a) / 2 + gm.subs(ms_, b) / 2 - K0 * t20.subs(ms_, a) * t20.subs(ms_, b))
    s_closed = (K0 + 2 * K1 + 2 * K2) * X_ * (4 * X_ - 3) / 30
    spins = [R(k, 2) for k in range(2, 13)]
    cells_ok, mean_ok = True, True
    for I_ in spins:
        ms_I, _, L_ = generator(I_, K0, K1, K2)
        d_ = len(ms_I)
        for a_ in range(d_):
            for b_ in range(d_):
                if sp.expand(L_[a_ * d_ + b_, a_ * d_ + b_] + cell_rate(ms_I[a_], ms_I[b_]).subs(Is, I_)) != 0:
                    cells_ok = False
        mean_ok &= sp.expand(-L_.trace() / d_**2 - s_closed.subs(Is, I_)) == 0
    claim("the cell diagonals of -L, (g_m + g_m')/2 - J0·t_m·t_m', closed polynomials in m and I(I+1), equal the "
          "generator's diagonal in every cell at I = 1, 3/2, ..., 6", cells_ok)
    claim("the mean rate s = -Tr L / d² = (J0 + 2J1 + 2J2)·X(4X - 3)/30, X = I(I + 1), at I = 1, 3/2, ..., 6", mean_ok)
    E_top = [cell_rate(Is, -Is) - s_closed,
             cell_rate(Is, 1 - Is) + cell_rate(Is - 1, -Is) - 2 * s_closed,
             cell_rate(Is, 2 - Is) + cell_rate(Is - 1, 1 - Is) + cell_rate(Is - 2, -Is) - 3 * s_closed]
    traces_ok = True
    for I_ in (R(2), R(5, 2), R(3)):
        ms_I, _, L_ = generator(I_, K0, K1, K2)
        d_ = len(ms_I)
        s_I = -L_.trace() / d_**2
        for k_, E_ in enumerate(E_top):
            p_ = d_ - 1 - k_
            block = [a_ * d_ + a_ + p_ for a_ in range(d_ - p_)]    # cells |m_a><m_b|, b = a + p, so m - m' = p
            traces_ok &= sp.expand(-sum(L_[c_, c_] for c_ in block) - (k_ + 1) * s_I - E_.subs(Is, I_)) == 0
    claim("the three trace conditions, as written in I, equal the generator's own traces of its blocks of order 2I, "
          "2I - 1, 2I - 2 minus 1, 2, 3 times s at I = 2, 5/2, 3", traces_ok)
    M_top = sp.Matrix([[sp.diff(e, k) for k in (K0, K1, K2)] for e in E_top])
    det_top = sp.factor(M_top.det())
    claim("the traces of the blocks of order 2I, 2I - 1, 2I - 2 equal 1, 2, 3 times s only where the matrix of the three "
          "linear forms is singular, det = I(I - 1)(2I - 3)²(2I - 1)³(2I + 5)(4I + 3)/60, non-zero at every I >= 2",
          sp.expand(det_top - Is * (Is - 1) * (2 * Is - 3) ** 2 * (2 * Is - 1) ** 3 * (2 * Is + 5) * (4 * Is + 3) / 60) == 0,
          f"[{det_top}]")
    ker1 = M_top.subs(Is, 1).nullspace()
    ker32 = M_top.subs(Is, R(3, 2)).nullspace()
    claim("at the two roots the forms leave what the moments and the theorem decide: at I = 1 the ray (0, 2, 1), at "
          "I = 3/2 the plane J0 = 0, which is the theorem's centre",
          len(ker1) == 1 and list(ker1[0] / ker1[0][2]) == [0, 2, 1]
          and len(ker32) == 2 and all(v[0] == 0 for v in ker32))

    print("\nevery offset and splitting, every spin I >= 2: a centre needs a repeated population rate (exact)")
    tt = sp.Symbol("t")
    wQs = sp.Symbol("wQ", real=True)
    ws = sp.Symbol("w", real=True)

    def populations(I_):
        """The populations block (order 0) of the J1 and J2 parts, rational similarity of the generator's."""
        _, blocks_ = rational_blocks(I_)
        d_ = int(2 * I_) + 1
        return blocks_[d_ - 1][1], blocks_[d_ - 1][2]
    blind_ok, conj_ok = True, True
    for I_ in (R(2), R(5, 2)):
        ms_I, _, L_ = generator(I_, K0, K1, K2, wQs, ws)
        d_ = len(ms_I)
        pop = [a_ * d_ + a_ for a_ in range(d_)]
        blind_ok &= not (L_.extract(pop, pop).free_symbols & {K0, wQs, ws})
        for p_ in range(1, d_):
            up_ = [a_ * d_ + a_ + p_ for a_ in range(d_ - p_)]
            dn_ = [(a_ + p_) * d_ + a_ for a_ in range(d_ - p_)]
            cu = sp.Poly(L_.extract(up_, up_).charpoly(x).as_expr(), x).all_coeffs()
            cd = sp.Poly(L_.extract(dn_, dn_).charpoly(x).as_expr(), x).all_coeffs()
            conj_ok &= all(sp.expand(sp.conjugate(a_) - b_) == 0 for a_, b_ in zip(cu, cd))
    claim("the populations block depends on J1 and J2 only (not on J0, the splitting or the offset), and the "
          "blocks of order p and -p have complex-conjugate characteristic polynomials, I = 2, 5/2, wQ and w symbols",
          blind_ok and conj_ok)

    def nu_closed(k, extra=3):
        """nu_k = Tr (W + c)^k at J1 = 1, J2 = t, c = -Tr W / d, as a polynomial in I and t. Each entry of W is a
        polynomial of degree 4 in m and I and the trace sums over m = -I .. I, so nu_k has degree <= 4k + 1 in I:
        interpolation at 4k + 2 spins fixes it, and three more spins check it."""
        deg = 4 * k + 1
        spins_ = [R(n, 2) for n in range(2, 2 + deg + 1 + extra)]
        forms = {}
        for I_ in spins_:
            W1_, W2_ = populations(I_)
            vals = []
            for j_ in range(k + 1):
                W_ = W1_ + j_ * W2_
                A_ = W_ - (W_.trace() / W_.rows) * sp.eye(W_.rows)
                vals.append((R(j_), (A_**k).trace()))
            forms[I_] = sp.Poly(sp.interpolate(vals, tt), tt)
        out, agree = 0, True
        for e in range(k + 1):
            q = sp.expand(sp.interpolate([(I_, forms[I_].coeff_monomial(tt**e)) for I_ in spins_[:deg + 1]], Is))
            agree &= all(q.subs(Is, I_) == forms[I_].coeff_monomial(tt**e) for I_ in spins_[deg + 1:])
            out += q * tt**e
        return sp.expand(out), agree
    nu3, ok3 = nu_closed(3)
    nu5, ok5 = nu_closed(5)
    claim("nu3, nu5 of the populations block as exact polynomials in I and t = J2/J1 (degrees <= 13 and 21 in I), "
          "fixed by 14 and 22 spins and checked at three more", ok3 and ok5)
    lead3 = sp.factor(nu3.coeff(tt, 3))                         # nu3 at J1 = 0, J2 = 1
    res_ = sp.factor(sp.resultant(nu3, nu5, tt))

    def no_spin_root(poly_expr):
        """True if the polynomial in I vanishes at no spin I >= 2 (exact: values at spins up to a bound past every
        real root, and no root beyond the bound)."""
        out = sp.expand(poly_expr) != 0                         # an identically zero polynomial decides nothing
        for f_, _ in sp.factor_list(poly_expr, Is)[1]:
            P_ = sp.Poly(f_, Is)
            if P_.degree() == 0:
                continue
            B = int(sp.ceiling(max([sp.N(r_) for r_ in P_.real_roots()] + [2]))) + 1
            out &= P_.count_roots(B, None) == 0 and all(P_.eval(R(n, 2)) != 0 for n in range(4, 2 * B + 1))
        return out
    claim("the populations block is palindromic at no spin I >= 2 for real rates other than J1 = J2 = 0: nu3 at "
          "J1 = 0 and the resultant of nu3 and nu5 in t vanish at no spin I >= 2",
          no_spin_root(lead3) and no_spin_root(res_), f"[resultant: {len(sp.factor_list(res_, Is)[1])} factors in I]")
    z1 = sp.Poly(sp.gcd(sp.Poly(nu3.subs(Is, 1), tt), sp.Poly(nu5.subs(Is, 1), tt)), tt)
    claim("the route reaches the known cases: nu3 and nu5 vanish identically at I = 3/2, and at I = 1 they share "
          "exactly the ratios J2/J1 = 1/8, 5/4 (the populations of the spin-1 rays) and -1",
          nu3.subs(Is, R(3, 2)) == 0 and nu5.subs(Is, R(3, 2)) == 0
          and set(sp.roots(z1, multiple=True)) == {R(1, 8), R(5, 4), -1})
    jacobi_ok = True
    for I_ in [R(n) for n in range(2, 11)]:
        W1_, _ = populations(I_)
        jacobi_ok &= all(W1_[i_, i_ + 1] * W1_[i_ + 1, i_] > 0 for i_ in range(W1_.rows - 1)) \
            and all(W1_[i_, j_] == 0 for i_ in range(W1_.rows) for j_ in range(W1_.rows) if abs(i_ - j_) > 1)
    from sympy.polys.matrices import DomainMatrix
    from sympy.polys.domains import GF
    prime = 1000003
    simple_half = []
    for n in range(5, 62, 2):
        W1_, W2_ = populations(R(n, 2))
        for j_ in (1, 2, 3, 5, 7):
            M_ = W1_ + j_ * W2_
            den = sp.ilcm(*[sp.fraction(v_)[1] for v_ in M_])
            cp = DomainMatrix.from_Matrix((den * M_).applyfunc(sp.Integer)).convert_to(GF(prime)).charpoly()
            f_ = sp.Poly([int(c_) % prime for c_ in cp], x, modulus=prime)
            if sp.degree(sp.gcd(f_, f_.diff()), x) == 0:
                simple_half.append(n)
                break
    claim("the repeated rates lie on finitely many ratios J2/J1 (plus J1 = 0): at integer I, J2 = 0 makes the "
          "populations block an irreducible tridiagonal chain, so simple (I = 2 .. 10 shown, the argument holds at "
          "every integer I); at half-integer I = 5/2 .. 61/2 a rational point with a characteristic polynomial "
          "squarefree mod 1000003", jacobi_ok and simple_half == list(range(5, 62, 2)))
    def flip_sectors(I_):
        """Characteristic polynomials (in x) of the populations block at J1 = 1, J2 = t on the m -> -m even and
        odd populations; the block is taken from the generator itself, where it is real symmetric."""
        ms_I, _, L_ = generator(I_, 0, 1, tt)
        d_ = len(ms_I)
        pop = [a_ * d_ + a_ for a_ in range(d_)]
        W_ = L_.extract(pop, pop).applyfunc(lambda z_: sp.nsimplify(sp.expand(z_)))
        sym_ok = (W_ - W_.T).applyfunc(sp.expand).is_zero_matrix
        cols = {+1: [], -1: []}
        for i_ in range(d_):
            j_ = d_ - 1 - i_
            for sgn in (+1, -1):
                if i_ < j_ or (i_ == j_ and sgn == 1):
                    v_ = sp.zeros(d_, 1)
                    v_[i_] += 1
                    v_[j_] += sgn
                    cols[sgn].append(v_)
        out = []
        for sgn in (+1, -1):
            B_ = sp.Matrix.hstack(*cols[sgn])
            out.append(sp.expand(((B_.T * B_).inv() * B_.T * W_ * B_).charpoly(x).as_expr()))
        return sym_ok, out
    sym2, (pe2, po2) = flip_sectors(R(2))
    cross2 = set(sp.Poly(sp.resultant(pe2, po2, x), tt).real_roots())
    within2 = set(sp.Poly(sp.discriminant(pe2, x), tt).real_roots())
    claim("the exceptional set is more than the two lines: at I = 2 the populations block (real symmetric) has a "
          "repeated rate at six real ratios J2/J1, five where an even rate under m -> -m crosses an odd one (two of "
          "them positive) and one, J2/J1 = -3/14, where two even rates meet",
          sym2 and len(cross2) == 5 and sum(1 for r_ in cross2 if r_ > 0) == 2 and within2 == {R(-3, 14)}
          and not sp.Poly(sp.discriminant(po2, x), tt).real_roots(), str(sorted(sp.N(r_, 4) for r_ in cross2)))
    W1_, W2_ = populations(R(7, 2))
    claim("at I = 7/2 and J1 = J2 (the populations of extreme narrowing, rotation invariant) two ranks share a "
          "rate: the populations decay at 0, 30, 84 twice, 150, 210 twice, 240 (times J1)",
          sp.expand((W1_ + W2_).charpoly(x).as_expr()
                    - x * (x + 30) * (x + 84)**2 * (x + 150) * (x + 210)**2 * (x + 240)) == 0)

    print("\na static splitting beyond spin 3/2 (w = 0, exact): centres appear at spin 1 only, on five rays")
    Q4 = sp.symbols("J0 J1 J2 q", real=True)
    mu32 = split_moments(R(3, 2), Q4)
    claim("I = 3/2 with the splitting: mu3 = 1296·J0·q², zero exactly at J0 = 0 (the theorem's centre) or q = 0 "
          "(the rule without a splitting)",
          sp.expand(mu32[3] - 1296 * Q4[0] * Q4[3] ** 2) == 0)
    spin1_rays = {(4, 0, R(-1, 36)), (-6, 2, R(1, 12)), (-12, -4, R(5, 12)), (0, -1, R(1, 48)), (0, R(4, 5), R(1, 60)),
                  (0, 8, R(7, 4))}
    for I_ in [R(1)] + [R(k, 2) for k in range(4, 13)]:
        pts, rest_empty = split_centre_points(split_moments(I_, Q4), Q4)
        want_pts = spin1_rays if I_ == 1 else set()
        claim(f"I = {I_}: mu3, mu5, mu7, mu9 vanish together on {len(pts)} rays (J0, J1, J2 = 1, q²) and nowhere else "
              f"but at all rates zero (outside the charts, L = -i[H, .], trivially centred)",
              rest_empty and set(pts) == want_pts, str(sorted(pts)) if pts else "")

    def centred_at(J0_, J1_, J2_, wQ2, w_=0):
        # the offset shifts the block of order p by -i·p·w: P_w(x) = prod_p P_p(x + i·p·w), P_p the block at w = 0
        _, _, L_ = generator(R(1), J0_, J1_, J2_, sp.sqrt(wQ2), 0)
        P_ = sp.Integer(1)
        for p in range(-2, 3):
            idx = [i * 3 + (i + p) for i in range(3) if 0 <= i + p < 3]
            P_ *= L_.extract(idx, idx).charpoly(x).as_expr().subs(x, x + sp.I * p * w_)
        P_ = sp.Poly(sp.expand(P_), x)
        s_ = P_.all_coeffs()[1] / P_.degree()
        return sp.expand(sp.simplify((-1) ** P_.degree() * P_.as_expr().subs(x, -x - 2 * s_) - P_.as_expr())) == 0, P_
    real_rays = [r for r in spin1_rays if r[2] > 0]
    claim("I = 1: at five of the six rays q² > 0, and at each of them the characteristic polynomial satisfies "
          "(-1)^n P(-x - 2s) = P(x), n = 9, exactly (at nine eigenvalues the four odd moments are the whole condition)",
          len(real_rays) == 5 and all(centred_at(a, b, 1, 6 * u_)[0] for a, b, u_ in real_rays))
    phys = sorted(r for r in real_rays if r[0] >= 0 and r[1] >= 0)
    ok_a, P_a = centred_at(0, R(4, 5), 1, R(1, 10))
    ok_b, P_b = centred_at(0, 8, 1, R(21, 2))
    claim("I = 1, non-negative rates: exactly two rays, J0 = 0 on both: J1 = 4J2/5 at wQ² = J2²/10, spectrum "
          "{0, -1², -(6/5)³, -(7/5)², -12/5}·J2 about -6/5·J2; J1 = 8J2 at wQ² = 21J2²/2, {0, -5², -6³, -7², -12}·J2 "
          "about -6·J2",
          phys == [(0, R(4, 5), R(1, 60)), (0, 8, R(7, 4))] and ok_a and ok_b
          and sp.expand(P_a.as_expr() - x * (x + 1) ** 2 * (x + R(6, 5)) ** 3 * (x + R(7, 5)) ** 2 * (x + R(12, 5))) == 0
          and sp.expand(P_b.as_expr() - x * (x + 5) ** 2 * (x + 6) ** 3 * (x + 7) ** 2 * (x + 12)) == 0)
    def rates_palindromic(J0_, J1_, J2_, wQ2, w_):
        _, _, L_ = generator(R(1), J0_, J1_, J2_, sp.sqrt(wQ2), 0)
        re_ = []
        for p in range(-2, 3):
            idx = [i * 3 + (i + p) for i in range(3) if 0 <= i + p < 3]
            for ev, am in sp.roots(sp.Poly(L_.extract(idx, idx).charpoly(x).as_expr(), x)).items():
                re_ += [sp.re(sp.nsimplify(ev) - sp.I * p * w_)] * am
        s_ = -sum(re_) / len(re_)
        return sorted(re_) == sorted(-r - 2 * s_ for r in re_)
    claim("I = 1: off the tuned splitting the centre is gone at both rays (wQ² detuned by 6/5 and by 21/20), and an "
          "offset w = 1 breaks the complex map at both, while the nine decay rates (real parts) stay palindromic in "
          "every frame tried (w = 1, 5)",
          not centred_at(0, R(4, 5), 1, R(6, 50))[0] and not centred_at(0, 8, 1, R(21, 2) * R(21, 20))[0]
          and not centred_at(0, R(4, 5), 1, R(1, 10), 1)[0] and not centred_at(0, 8, 1, R(21, 2), 1)[0]
          and all(rates_palindromic(0, R(4, 5), 1, R(1, 10), w_) and rates_palindromic(0, 8, 1, R(21, 2), w_) for w_ in (1, 5)))

    # the second route, by hand: the block polynomials at I = 1 in the rates r = -lambda
    K0_, K1_, K2_ = sp.symbols("J0 J1 J2", real=True)
    W = sp.Symbol("wQ", positive=True)
    _, _, L1 = generator(R(1), K0_, K1_, K2_, W, 0)
    blk = {}
    for p in (0, 1, 2):
        idx = [i * 3 + (i + p) for i in range(3) if 0 <= i + p < 3]
        blk[p] = sp.factor(sp.Poly(sp.expand(L1.extract(idx, idx).charpoly(x).as_expr().subs(x, -x)), x).monic().as_expr())
    two_s = sp.Rational(4, 3) * (K0_ / 2 + K1_ + K2_)
    sq = sp.Poly(sp.expand(16 * blk[1]), x)
    sigma1, prod1 = -sq.coeff_monomial(x) / 16, sq.coeff_monomial(1) / 16
    cq = (K1_ + 2 * K2_) / 2
    claim("I = 1 by hand: populations 0, 3J1/2, (J1 + 4J2)/2; double quantum (J1 + 2J2)/2; the single-quantum "
          "pair has sum (3J0 + 3J1 + 2J2)/2 and carries 3wQ²/2 in its product, the only place wQ enters; "
          "2s = 4(J0/2 + J1 + J2)/3",
          sp.expand(blk[0] - x * (x - R(3, 2) * K1_) * (x - (K1_ + 4 * K2_) / 2)) == 0
          and sp.expand(blk[2] - (x - cq)) == 0 and sp.expand(sigma1 - (3 * K0_ + 3 * K1_ + 2 * K2_) / 2) == 0
          and sp.expand(prod1 - prod1.subs(W, 0) - R(3, 2) * W ** 2) == 0 and sp.expand(-L1.trace() * 2 / 9 - two_s) == 0)
    s_ = two_s / 2
    j0_forced = sp.solve(sp.Eq(sigma1, 3 * s_ - cq), K0_)
    hand = []
    for single in (sp.Eq(R(3, 2) * K1_, s_), sp.Eq((K1_ + 4 * K2_) / 2, s_)):
        k1 = sp.solve(single.subs(K0_, 0), K1_)[0]
        w2 = sp.solve(sp.Eq(prod1.subs({K0_: 0, K1_: k1}), (s_ * (2 * s_ - cq)).subs({K0_: 0, K1_: k1})), W)
        hand.append((k1, sorted({sp.factor(v ** 2) for v in w2}, key=str)))
    claim("I = 1 by hand: one single-quantum rate at s and the other mirroring the double quantum forces J0 = 0; the "
          "populations then fix J1 = 4J2/5 or 8J2, the single-quantum product wQ² = J2²/10 or 21J2²/2: the same two rays",
          j0_forced == [0] and hand == [(4 * K2_ / 5, [K2_ ** 2 / 10]), (8 * K2_, [21 * K2_ ** 2 / 2])])

    def block_structure(J0_, J1_, J2_, wQ2):
        _, _, L_ = generator(R(1), J0_, J1_, J2_, sp.sqrt(wQ2), 0)
        L_ = L_.applyfunc(lambda z: sp.nsimplify(sp.expand(z)))
        out = {}
        for p in range(-2, 3):
            idx = [i * 3 + (i + p) for i in range(3) if 0 <= i + p < 3]
            B = L_.extract(idx, idx)
            for ev, am in sp.roots(sp.Poly(B.charpoly(x).as_expr(), x)).items():
                out.setdefault(ev, []).append((p, am, len((B - ev * sp.eye(B.rows)).nullspace())))
        return out
    bs = block_structure(0, R(4, 5), 1, R(1, 10))
    bs2 = block_structure(0, 8, 1, R(21, 2))
    claim("I = 1 at both rays: every eigenvalue is simple in its block and semisimple (no exceptional point); at "
          "J1 = 4J2/5 the mirror pairs 0 with -12/5 among the populations and -1 (orders ±1) with -7/5 (orders ±2), "
          "-6/5 on the centre; at J1 = 8J2, -7 (orders ±1) with -5 (orders ±2), -6 on the centre",
          all(am == 1 and gm == 1 for b_ in (bs, bs2) for v in b_.values() for _, am, gm in v)
          and sorted(p for p, _, _ in bs[-1]) == [-1, 1] and sorted(p for p, _, _ in bs[R(-7, 5)]) == [-2, 2]
          and sorted(p for p, _, _ in bs[R(-6, 5)]) == [-1, 0, 1] and bs[0] == [(0, 1, 1)] and bs[R(-12, 5)] == [(0, 1, 1)]
          and sorted(p for p, _, _ in bs2[-7]) == [-1, 1] and sorted(p for p, _, _ in bs2[-5]) == [-2, 2]
          and sorted(p for p, _, _ in bs2[-6]) == [-1, 0, 1])
    disc = sp.factor(sp.Poly(sp.expand(16 * blk[1]), x).discriminant())
    paired = []
    for J1_, w2 in ((R(4, 5), R(1, 10)), (8, R(21, 2))):
        for wq2 in (0, w2):
            ev = sp.roots(sp.Poly(sp.expand(generator(R(1), 0, J1_, 1, sp.sqrt(wq2), 0)[2].charpoly(x).as_expr()), x), multiple=True)
            sh = -sum(ev) / len(ev)
            paired.append(sum(1 for e in ev if any(sp.simplify(e + f + 2 * sh) == 0 for f in ev)))
    claim("I = 1: the single-quantum pair's exceptional point (zero discriminant) is 6wQ² = J1² at J0 = 0; the rays sit "
          "at 6wQ²/J1² = 15/16 and 63/64, just below it; without the splitting 3 of 9 eigenvalues pair about -s, with it 9; "
          "the pair is drawn from {7/10, 3/2}·J2 at the first ray",
          sp.expand(disc.subs(K0_, 0) - 256 * (K1_ ** 2 - 6 * W ** 2)) == 0
          and [6 * R(1, 10) / R(4, 5) ** 2, 6 * R(21, 2) / 64] == [R(15, 16), R(63, 64)] and paired == [3, 9, 3, 9]
          and set(sp.roots(sp.Poly(blk[1].subs({K0_: 0, K1_: R(4, 5), K2_: 1, W: 0}), x))) == {R(7, 10), R(3, 2)})

    held = sum(CLAIMS)
    print(f"\n{held} of {len(CLAIMS)} claims hold.")
    if held < len(CLAIMS):
        sys.exit(1)


if __name__ == "__main__":
    main()
