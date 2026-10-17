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
  - the global centre: the trace fixes 2s = 3J0 + 6J1 + 6J2; the coefficients of P(-x - 2s) - P(x) have the
    lex Groebner basis {J0·wQ², J0·(9J1² + w²)(9J2² + 4w²)}, whose real zeros are the theorem; the three
    reflections at a vanishing density; and sixteen exact points, both maps, both outcomes.

Import-inert; prints only, and exits with status 1 if any check fails.
"""

import sys
from collections import Counter

import sympy as sp

CLAIMS = []


def claim(label, ok, detail=""):
    CLAIMS.append(bool(ok))
    print(f"  [{'ok' if ok else 'FAIL'}] {label} {detail}")


def build():
    J0, J1, J2, wQ, w = sp.symbols("J0 J1 J2 wQ w", real=True)
    I = sp.Rational(3, 2)
    ms = [sp.Rational(3, 2), sp.Rational(1, 2), -sp.Rational(1, 2), -sp.Rational(3, 2)]
    d = 4
    Iz = sp.diag(*ms)
    Ip = sp.zeros(d)
    for k in range(1, d):
        Ip[k - 1, k] = sp.sqrt(I * (I + 1) - ms[k] * (ms[k] + 1))
    Im = Ip.T
    T = {0: (3 * Iz * Iz - I * (I + 1) * sp.eye(d)) / sp.sqrt(6),
         1: -(Iz * Ip + Ip * Iz) / 2, -1: (Iz * Im + Im * Iz) / 2,
         2: Ip * Ip / 2, -2: Im * Im / 2}
    Id = sp.eye(d)
    kron = sp.kronecker_product
    H = wQ * T[0] + w * Iz
    L = -sp.I * (kron(H, Id) - kron(Id, H.T))
    for g, A in [(J0, T[0]), (J1, T[1]), (J1, T[-1]), (J2, T[2]), (J2, T[-2])]:
        AdA = A.H * A
        L += g * (kron(A, A.conjugate()) - kron(AdA, Id) / 2 - kron(Id, AdA.T) / 2)
    return (J0, J1, J2, wQ, w), ms, T, sp.simplify(L)


def main():
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
          {6 * J1, 6 * J2} <= set(pop_rates) and sp.expand(sp.expand(blocks[3] - x).subs(w, 0) - (6 * J1 + 6 * J2) / 2) == 0)

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
    # the closed rates at wQ = w = 0, taken from the gated list above, sixteen with multiplicity
    rates = (expect[0] + [expect[1][0]] * 2 + expect[1][1:] * 2 + expect[2] * 2 + expect[3] * 2)
    claim("the closed rates count sixteen and sum to 16·s", len(rates) == 16 and sp.expand(sum(rates) - 8 * s2) == 0)
    for z in (J0, J1, J2):
        R0, S0 = [sp.expand(sp.sympify(r).subs(z, 0)) for r in rates], sp.expand(s2.subs(z, 0))
        claim(f"{z} = 0, wQ = w = 0: the sixteen rates are invariant under r -> {S0} - r",
              Counter(R0) == Counter(sp.expand(S0 - r) for r in R0))

    def spectrum(vals):
        roots = []
        for p in range(-3, 4):
            roots += [r for r, k in sp.roots(sp.Poly(blocks[p].subs(vals), x)).items() for _ in range(k)]
        return roots

    def has_global_centre(vals):
        """Exact, under both maps; they must agree, the spectrum being closed under conjugation."""
        lam = spectrum(vals)
        s = -sum(sp.re(r) for r in lam) / len(lam)   # a pairing fixes its centre at minus the mean
        canon = lambda z: sp.nsimplify(sp.expand(sp.radsimp(z)))
        own = Counter(canon(r) for r in lam)
        conj = own == Counter(canon(-sp.conjugate(r) - 2 * s) for r in lam)
        flat = own == Counter(canon(-r - 2 * s) for r in lam)
        return conj, flat

    R = sp.Rational
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
        conj, flat = has_global_centre(vals)
        zero = all(sp.expand(sp.radsimp(c.subs(vals))) == 0 for c in coeffs)
        claim(f"{name}: global centre {'exists' if expected else 'does not exist'} "
              "(roots under both maps, and the polynomial identity at the point)",
              conj == flat == zero == expected)

    held = sum(CLAIMS)
    print(f"\n{held} of {len(CLAIMS)} claims hold.")
    if held < len(CLAIMS):
        sys.exit(1)


if __name__ == "__main__":
    main()
