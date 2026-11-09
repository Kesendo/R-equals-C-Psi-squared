"""The third centred moment under cooling and heating, and a conditioned field that pairs (F137), exact.

L = -i[H, .] + sum_l (cool_l D[sigma^-_l] + heat_l D[sigma^+_l]), sigma^- = |0><1| with Z|0> = |0>, H a sum of
Pauli strings other than the identity with real coefficients (a constant does not enter L), N qubits, d = 2^N, the centre -c/2 = Tr L / d^2. The law:

  Tr((L + c/2)^3) = 12*4^(N-1) * sum_l s_l (hX_l^2 + hY_l^2)  -  3*2^(N-1) * sum_m (cool_m - heat_m) Tr(Z_m H^2),

s_l = (cool_l + heat_l)/2, hX_l and hY_l the coefficients of X_l and Y_l in H. Two routes meet: the left side as the
trace of the cubed generator (matrices in exact rationals, computational basis), the right side from Pauli-string
algebra alone (the field coefficients, and Tr(Z_m H^2) as 2^N times the signed products of term pairs whose product
is +-Z_m).

The conditioned field H = X_0 + X_0 Z_1 = 2 X_0 |0><0|_1, site 0 cooled and heated, site 1 under one kind of jump
only: every odd centred power sum up to the dimension 16 vanishes exactly at the rates where the field term and
the arrow term cancel, so the characteristic polynomial of L + c/2 is even (Newton's identities) and the spectrum
pairs. Symbolically, with H = h (X_0 + X_0 Z_1) + w Z_1 and the rates a, b (site 0) and g (site 1) as symbols,
every odd coefficient of the characteristic polynomial of L + c/2 vanishes on a + b = 2g, and the coefficient of
x^13 is -8 h^2 (a + b - 2g): the pairing holds exactly when s_0 = g, at every field and every frequency of site 1.
H conserves Z_1 and site 1 moves only one way, so L is block triangular in site 1's ket-bra sectors; that is
checked entry by entry. The controls must be nonzero: the third moment at rates off the cancellation, and, where
the third moment cancels and the pairing still fails, the fifth.

Run: python simulations/thermal_third_moment.py   (prints "thermal_third_moment: ALL GATES PASS", exit 0)
"""
from fractions import Fraction as F
import random
import sys

# ---------- Pauli-string algebra (the right side) ----------
PM = {('X', 'Y'): (1j, 'Z'), ('Y', 'X'): (-1j, 'Z'), ('Y', 'Z'): (1j, 'X'),
      ('Z', 'Y'): (-1j, 'X'), ('Z', 'X'): (1j, 'Y'), ('X', 'Z'): (-1j, 'Y')}


def pmul(a, b):
    phase, out = 1, []
    for x, y in zip(a, b):
        if x == 'I':
            out.append(y)
        elif y == 'I' or x == y:
            out.append('I' if x == y else x)
        else:
            p, r = PM[(x, y)]
            phase *= p
            out.append(r)
    return phase, ''.join(out)


def tr_z_h2(n, terms, m):
    """Tr(Z_m H^2) / 2^N, from the term pairs whose product is a real multiple of Z_m."""
    target = 'I' * m + 'Z' + 'I' * (n - m - 1)
    total = F(0)
    for a, ca in terms:
        for b, cb in terms:
            ph, prod = pmul(a, b)
            if prod == target:
                assert ph in (1, -1, 1j, -1j)
                if ph in (1, -1):
                    total += int(ph.real) * F(ca) * F(cb)
                # an imaginary phase comes from an anticommuting pair, cancelled by its reversed product
    return total


def law(n, terms, cool, heat):
    terms = [(s, c) for s, c in terms if s != 'I' * n]   # a constant does not enter L
    field = F(0)
    for s, c in terms:
        supp = [k for k, ch in enumerate(s) if ch != 'I']
        if len(supp) == 1 and s[supp[0]] in 'XY':
            field += F(cool[supp[0]] + heat[supp[0]], 2) * F(c) ** 2
    arrow = sum(F(cool[m] - heat[m]) * tr_z_h2(n, terms, m) for m in range(n)) * 2 ** n
    return 12 * 4 ** (n - 1) * field - 3 * 2 ** (n - 1) * arrow, field, arrow


# ---------- exact matrices (the left side) ----------
# entries are Gaussian rationals held as (re, im) pairs of Fractions
Z0 = (F(0), F(0))
LET = {'I': [[(F(1), F(0)), Z0], [Z0, (F(1), F(0))]],
       'X': [[Z0, (F(1), F(0))], [(F(1), F(0)), Z0]],
       'Y': [[Z0, (F(0), F(-1))], [(F(0), F(1)), Z0]],
       'Z': [[(F(1), F(0)), Z0], [Z0, (F(-1), F(0))]]}
SM = [[Z0, (F(1), F(0))], [Z0, Z0]]          # sigma^- = |0><1|
SP = [[Z0, Z0], [(F(1), F(0)), Z0]]          # sigma^+ = |1><0|


def gmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def gadd(x, y):
    return (x[0] + y[0], x[1] + y[1])


def gscale(x, s):
    return (x[0] * s, x[1] * s)


def kron(a, b):
    n, m = len(a), len(b)
    return [[gmul(a[i // m][j // m], b[i % m][j % m]) for j in range(n * m)] for i in range(n * m)]


def kron_all(ms):
    r = [[(F(1), F(0))]]
    for m in ms:
        r = kron(r, m)
    return r


def mat_add(a, b, s=1):
    return [[gadd(x, gscale(y, s)) for x, y in zip(r, q)] for r, q in zip(a, b)]


def mat_mul(a, b):
    n, k = len(a), len(b[0])
    out = []
    for i in range(n):
        row = []
        for j in range(k):
            acc = Z0
            for t in range(len(b)):
                x = a[i][t]
                if x == Z0:
                    continue
                y = b[t][j]
                if y == Z0:
                    continue
                acc = gadd(acc, gmul(x, y))
            row.append(acc)
        out.append(row)
    return out


def transpose(a):
    return [list(r) for r in zip(*a)]


def conj(a):
    return [[(x[0], -x[1]) for x in r] for r in a]


def eye(dim):
    return [[(F(1), F(0)) if i == j else Z0 for j in range(dim)] for i in range(dim)]


def generator(n, terms, cool, heat):
    """L acting on vec(rho), row-major: vec(A rho B) = (A kron B^T) vec(rho)."""
    d = 2 ** n
    Id = eye(d)
    H = [[Z0] * d for _ in range(d)]
    for s, c in terms:
        H = mat_add(H, kron_all([LET[ch] for ch in s]), F(c))
    comm = mat_add(kron(H, Id), kron(Id, transpose(H)), -1)
    L = [[gmul((F(0), F(-1)), x) for x in r] for r in comm]
    for l in range(n):
        for A0, r in ((SM, cool[l]), (SP, heat[l])):
            if r == 0:
                continue
            A = kron_all([A0 if k == l else LET['I'] for k in range(n)])
            AdA = mat_mul(transpose(conj(A)), A)
            D = mat_add(mat_add(kron(A, conj(A)), kron(AdA, Id), F(-1, 2)), kron(Id, transpose(AdA)), F(-1, 2))
            L = mat_add(L, D, F(r))
    return L


def centred(L):
    dim = len(L)
    tr = Z0
    for i in range(dim):
        tr = gadd(tr, L[i][i])
    assert tr[1] == 0
    shift = (-tr[0] / dim, F(0))
    return [[gadd(x, shift) if i == j else x for j, x in enumerate(r)] for i, r in enumerate(L)]


def trace(a):
    t = Z0
    for i in range(len(a)):
        t = gadd(t, a[i][i])
    return t


def odd_power_sums(M, top):
    out, P = {}, M
    for m in range(1, top + 1):
        if m > 1:
            P = mat_mul(P, M)
        if m % 2:
            out[m] = trace(P)
    return out


def third(n, terms, cool, heat):
    M = centred(generator(n, terms, cool, heat))
    t = trace(mat_mul(mat_mul(M, M), M))
    assert t[1] == 0, t
    return t[0]


def chain(n, templates):
    out = []
    for tpl, c in templates:
        for k in range(n - len(tpl) + 1):
            out.append(('I' * k + tpl + 'I' * (n - k - len(tpl)), c))
    return out


FAIL = []


def check(ok, msg):
    print(("  ok   " if ok else "  FAIL ") + msg)
    if not ok:
        FAIL.append(msg)


def gate_law():
    print("the law: Tr((L + c/2)^3) from the matrices against the Pauli-algebra closed form")
    heis3 = chain(3, [('XX', 1), ('YY', 1), ('ZZ', 1)])
    named = [
        # (label, N, terms, cool, heat, value recorded in the repo or None)
        ("Heisenberg N=3, uniform X field 1, cooling 1 (F137: 288)", 3, heis3 + chain(3, [('X', 1)]), [1] * 3, [0] * 3, F(288)),
        ("Heisenberg N=3, uniform Z field 1, cooling 1 (F137: -768 as written)", 3, heis3 + chain(3, [('Z', 1)]), [1] * 3, [0] * 3, F(-768)),
        ("Heisenberg N=3, Z fields (1/2, -1/2, 0), cooling 1 (CAUGHT_ERRORS: 96)", 3,
         heis3 + [('ZII', F(1, 2)), ('IZI', F(-1, 2))], [1] * 3, [0] * 3, F(96)),
        ("Heisenberg N=3, Z fields (3/10, -1/5, 1/2), cooling 1 (F137, gate M11: -384/5 as written)", 3,
         heis3 + [('ZII', F(3, 10)), ('IZI', F(-1, 5)), ('IIZ', F(1, 2))], [1] * 3, [0] * 3, F(-384, 5)),
        ("transverse DM YZ - ZY, N=3, cooling 1 (a two-site term: no third moment)", 3,
         chain(3, [('YZ', 1), ('ZY', -1)]), [1] * 3, [0] * 3, F(0)),
        ("X and Y fields on one site with two rates (the arrow pair X, Y anticommutes: no arrow)", 2,
         [('XI', 1), ('YI', 2)], [3, 2], [1, 1], None),
        ("X_0 + X_0 Z_1 off the cancellation", 2, [('XI', 1), ('XZ', 1)], [1, 2], [1, 0], None),
        ("a constant beside a Z field and XX (the law reads the traceless part)", 2,
         [('II', 1), ('ZI', 1), ('XX', F(1, 2))], [1, 0], [0, 0], F(0)),
    ]
    rng = random.Random(20261109)
    randoms = []
    coef = lambda: rng.choice([-3, -2, -1, 1, 2, 3])
    for k in range(16):
        n = 2 if k < 8 else 3
        terms, seen = [], set()

        def put(s):
            if s != 'I' * n and s not in seen:
                seen.add(s)
                terms.append((s, coef()))
        # a transverse field on a drawn site, a drawn string T beside T*Z_m where the two commute, then drawn strings
        site = rng.randrange(n)
        put('I' * site + rng.choice('XY') + 'I' * (n - site - 1))
        while True:
            t = ''.join(rng.choice('IXYZ') for _ in range(n))
            m = rng.randrange(n)
            ph, partner = pmul(t, 'I' * m + 'Z' + 'I' * (n - m - 1))
            if t != 'I' * n and ph in (1, -1):
                put(t)
                put(partner)
                break
        for _ in range(rng.randint(1, 4)):
            put(''.join(rng.choice('IXYZ') for _ in range(n)))
        randoms.append((f"drawn #{k}", n, terms, [rng.randint(0, 3) for _ in range(n)], [rng.randint(0, 2) for _ in range(n)], None))
    both_nonzero = 0
    for label, n, terms, cool, heat, recorded in named + randoms:
        t = third(n, terms, cool, heat)
        value, field, arrow = law(n, terms, cool, heat)
        if field != 0 and arrow != 0:
            both_nonzero += 1
        rec = "" if recorded is None else f", recorded {recorded}"
        check(t == value and (recorded is None or t == recorded),
              f"{label}: trace {t}, law {value} (field part {field}, arrow part {arrow}){rec}")
    check(both_nonzero >= 3, f"{both_nonzero} settings carry both a nonzero field part and a nonzero arrow part")


def gate_conditioned_field():
    print("the conditioned field X_0 (1 + Z_1): every odd centred power sum up to the dimension 16, exactly")
    cases = [
        ("X_0 + X_0 Z_1, site 0 at 1/1, site 1 cooled at 1", [('XI', 1), ('XZ', 1)], [1, 1], [1, 0]),
        ("X_0 + X_0 Z_1, site 0 cooled at 2, site 1 cooled at 1", [('XI', 1), ('XZ', 1)], [2, 1], [0, 0]),
        ("X_0 - X_0 Z_1, site 0 at 1/1, site 1 heated at 1", [('XI', 1), ('XZ', -1)], [1, 0], [1, 1]),
        ("Y_0 + Y_0 Z_1, site 0 at 1/1, site 1 cooled at 1", [('YI', 1), ('YZ', 1)], [1, 1], [1, 0]),
    ]
    for label, terms, cool, heat in cases:
        value, _, _ = law(2, terms, cool, heat)
        sums = odd_power_sums(centred(generator(2, terms, cool, heat)), 16)
        nonzero = [m for m, t in sums.items() if t != Z0]
        check(value == 0 and not nonzero, f"{label}: law {value}, odd power sums nonzero at {nonzero or 'no m <= 16'}")
    print("  block triangularity in site 1's ket-bra sectors (a, b), X_0 + X_0 Z_1 with site 1 cooled only")
    L = generator(2, [('XI', 1), ('XZ', 1)], [1, 1], [1, 0])
    # vec index (ket, bra) with ket = 2*q0 + q1 (site 0 the high bit), likewise bra; row-major: idx = ket*4 + bra
    sector = lambda idx: ((idx // 4) & 1, (idx % 4) & 1)
    leaks = [(i, j) for i in range(16) for j in range(16) if L[i][j] != Z0 and sector(i) != sector(j)
             and not (sector(j) == (1, 1) and sector(i) == (0, 0))]
    moves = [(i, j) for i in range(16) for j in range(16) if L[i][j] != Z0 and sector(j) == (1, 1) and sector(i) == (0, 0)]
    check(not leaks and moves, f"entries leaving a sector: {len(leaks)} besides the arrow (1,1) -> (0,0), which has {len(moves)}")


def gate_symbolic():
    print("the conditioned field with symbols: H = h (X_0 + X_0 Z_1) + w Z_1, site 0 at a/b, site 1 cooled at g")
    import sympy as sp
    from sympy.polys.matrices import DomainMatrix
    a, b, g, h, w = sp.symbols('a b g h w')

    def piece(terms, cool, heat):
        L = generator(2, terms, cool, heat)
        return sp.Matrix(16, 16, lambda i, j: sp.Rational(L[i][j][0].numerator, L[i][j][0].denominator)
                         + sp.I * sp.Rational(L[i][j][1].numerator, L[i][j][1].denominator))
    L = (h * piece([('XI', 1), ('XZ', 1)], [0, 0], [0, 0]) + w * piece([('IZ', 1)], [0, 0], [0, 0])
         + a * piece([], [1, 0], [0, 0]) + b * piece([], [0, 0], [1, 0]) + g * piece([], [0, 1], [0, 0]))
    M = L - L.trace() / 16 * sp.eye(16)
    cp = [sp.expand(k.as_expr() if hasattr(k, 'as_expr') else k) for k in DomainMatrix.from_Matrix(M).charpoly()]
    odd = {16 - k: cp[k] for k in range(17) if (16 - k) % 2 == 1}
    on_line = all(sp.expand(co.subs(a, 2 * g - b)) == 0 for co in odd.values())
    first = sp.expand(odd[13] + 8 * h ** 2 * (a + b - 2 * g)) == 0
    check(on_line, "every odd coefficient of the characteristic polynomial vanishes on a + b = 2g")
    check(first, f"the coefficient of x^13 is -8 h^2 (a + b - 2g): {sp.factor(odd[13])}")

    # the four blocks of site 1's ket-bra sectors (row-major vec index ket*4 + bra, site 1 the low bit of each)
    x = sp.Symbol('x')
    s0 = (a + b) / 2

    def block_poly(mat, kb):
        idx = [i for i in range(16) if ((i // 4) & 1, (i % 4) & 1) == kb]
        co = DomainMatrix.from_Matrix(mat.extract(idx, idx)).charpoly()
        co = [k.as_expr() if hasattr(k, 'as_expr') else k for k in co]
        return sp.expand(sum(c * x ** (len(co) - 1 - k) for k, c in enumerate(co)))
    full = sp.expand(sum((k.as_expr() if hasattr(k, 'as_expr') else k) * x ** (16 - n)
                         for n, k in enumerate(DomainMatrix.from_Matrix(L).charpoly())))
    p = {kb: block_poly(L, kb) for kb in ((0, 0), (0, 1), (1, 0), (1, 1))}
    check(sp.expand(p[(0, 0)] * p[(0, 1)] * p[(1, 0)] * p[(1, 1)] - full) == 0,
          "the four blocks' characteristic polynomials multiply to L's")
    pc = {kb: block_poly(M, kb) for kb in ((0, 1), (1, 0))}
    check(sp.expand(pc[(0, 1)] - pc[(1, 0)].subs(x, -x)) == 0,
          "the one-sided blocks pair with each other about the centre, at every rate, h and w")
    check(sp.expand(p[(1, 1)] - (x + g) * (x + s0 + g) ** 2 * (x + 2 * s0 + g)) == 0,
          "the field-free block holds -g, -s0 - g twice and -2 s0 - g")
    q, r = sp.div(sp.Poly(p[(0, 0)], x), sp.Poly(x * (x + s0), x))
    check(r.is_zero and q.degree() == 2 and sp.expand(q.all_coeffs()[1] - 3 * s0) == 0,
          f"the field's block holds 0, -s0 and a pair of mean -(3/2) s0: the pair's polynomial {sp.factor(q.as_expr())}")


def gate_controls():
    print("the controls, each a nonzero exact value")
    t = third(2, [('XI', 1), ('XZ', 1)], [1, 2], [1, 0])
    check(t == -48, f"X_0 + X_0 Z_1 at site 1 cooled at 2: Tr((L + c/2)^3) = {t}, the law's 48*(1 - 2) = -48")
    for label, terms, cool, heat in [
        ("ZZ added", [('XI', 1), ('XZ', 1), ('ZZ', 1)], [1, 1], [1, 0]),
        ("site 1 cooled at 2 and heated at 1", [('XI', 1), ('XZ', 1)], [1, 2], [1, 1]),
        ("field 2 X_0 + X_0 Z_1, site 1 cooled at 2 (the field on in both sectors)", [("XI", 2), ("XZ", 1)], [1, 2], [1, 0]),
    ]:
        fifth_breaks(label, 2, terms, cool, heat)
    fifth_breaks("two neighbours controlled by site 1, X_0 (1 + Z_1) + X_2 (1 + Z_1), site 1 cooled at 1, the ends at 1/1", 3,
                 [('XII', 1), ('XZI', 1), ('IIX', 1), ('IZX', 1)], [1, 1, 1], [1, 0, 1])


def fifth_breaks(label, n, terms, cool, heat):
    value, _, _ = law(n, terms, cool, heat)
    sums = odd_power_sums(centred(generator(n, terms, cool, heat)), 5)
    check(value == 0 and sums[3] == Z0 and sums[5] != Z0,
          f"{label}: the third moment cancels (law {value}, trace {sums[3][0]}), the fifth does not ({sums[5][0]})")


if __name__ == "__main__":
    gate_law()
    gate_conditioned_field()
    gate_symbolic()
    gate_controls()
    if FAIL:
        print(f"thermal_third_moment: {len(FAIL)} GATE(S) FAILED")
        sys.exit(1)
    print("thermal_third_moment: ALL GATES PASS")
