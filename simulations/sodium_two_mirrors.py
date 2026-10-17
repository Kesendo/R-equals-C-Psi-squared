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
    than top-price cells; no global centre at eight sample rate triples;
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

    def centred(I_, vals):
        _, _, L_ = generator(I_, *vals)
        P_ = sp.Poly(sp.expand(L_.applyfunc(sp.nsimplify).charpoly(x).as_expr()), x)
        if not all(c.is_Rational for c in P_.all_coeffs()):          # else the identity is not decided exactly
            return None
        n_ = P_.degree()
        s_ = P_.all_coeffs()[1] / n_
        return sp.expand((-1) ** n_ * P_.as_expr().subs(x, -x - 2 * s_) - P_.as_expr()) == 0

    points = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (0, R(2, 5), R(7, 3)), (1, R(1, 3), 0), (R(3, 4), 0, R(5, 2)),
              (1, 1, 1), (1, R(37, 100), R(81, 100))]
    for I_ in spins:
        got = [centred(I_, v) for v in points]
        want = [v[0] * v[1] * v[2] == 0 if I_ == R(3, 2) else False for v in points]
        claim(f"I = {I_}: one centre for the whole spectrum at {sum(got)} of the 8 sample rate triples "
              f"({'exactly where one density vanishes' if I_ == R(3, 2) else 'none'})", got == want)

    held = sum(CLAIMS)
    print(f"\n{held} of {len(CLAIMS)} claims hold.")
    if held < len(CLAIMS):
        sys.exit(1)


if __name__ == "__main__":
    main()
