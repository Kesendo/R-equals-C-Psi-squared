"""The far element's locality and sheets on the carriers of PROOF_PALINDROME_COMPLEMENT_CONNECTION, and blindness beside.

Reads, with framework's far_element_shape (simulations/framework/diagnostics/far_element_shape.py), every carrier the
section "Blind dephased seats, and palindromes beyond colourings and site symmetries" and its Open items hold, and the
two grids of TrianglePalindromeLocusWitness (its Grid and PlanarGrid, copied here as data): where a row pairs through
one far element W, W's operator Schmidt ranks, and its sheets (the Pauli strings g with g·W·g = ±W beyond the jumps
and the g that do not split H). Blindness is F157's as the section states it: the Krylov space of the dephased seat
under the J-weighted Laplacian of the Heisenberg bonds is smaller than N. Exact throughout: the ends over Q(i), the
Krylov ranks over Q, and the N = 7 carrier W = V·W0·V^-1 built from its closed form over Q(i), no threshold.

Asserted, each a sentence of the section: the carriers' ranks across every cut (at N = 4 and N = 7 the multi-site
cuts too, N = 7's element checked in the far end exactly); which of them a site permutation fixing the watched seat
turns into a product (B's alone; at N = 7 all 720 read exactly, none); the pencil of C's sheet (C pairs at lam = 2 and 1/3, E breaks); and on the two
grids, counted as distinct Hamiltonians since the grids repeat rows, that every blind row that pairs through one
element has a sheet, which of those carriers are a permutation times a product (B's swap), and that the rows
without one are E's eight. A reading, not a proof: the grids are finite, the
triangle is N = 3, and every class is read at the row named.
Run: python simulations/far_element_sheets_census.py   (about two minutes; writes results/far_element_sheets_census.txt)
"""
from __future__ import annotations

import itertools
import sys
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ_I

sys.path.insert(0, str(Path(__file__).parent))
from framework.diagnostics.far_element_shape import far_element_shape, operator_schmidt_rank, pencil, sheets  # noqa: E402
from framework.pauli import pauli_product  # noqa: E402

OUT = []


def say(*a):
    line = " ".join(str(x) for x in a)
    print(line, flush=True)
    OUT.append(line)


def heis(n, bonds, fields, letters='XYZ'):
    h = {}
    for (i, j), w in bonds.items():
        for c in letters:
            s = ['I'] * n; s[i] = c; s[j] = c; k = ''.join(s); h[k] = h.get(k, 0) + sp.nsimplify(w)
    for site, f in fields.items():
        for c, x in zip('XYZ', f):
            if x:
                s = ['I'] * n; s[site] = c; k = ''.join(s); h[k] = h.get(k, 0) + sp.nsimplify(x)
    return {k: v for k, v in h.items() if v != 0}


def blind(n, bonds, seat=0):
    L = sp.zeros(n)
    for (i, j), w in bonds.items():
        w = sp.nsimplify(w); L[i, i] += w; L[j, j] += w; L[i, j] -= w; L[j, i] -= w
    e = sp.zeros(n, 1); e[seat] = 1
    return sp.Matrix.hstack(*[L ** k * e for k in range(n)]).rank() < n


def mul(A, B):
    out = {}
    for s_, a in A.items():
        for t_, b in B.items():
            e, u = pauli_product(s_, t_)
            out[u] = out.get(u, 0) + a * b * sp.I ** e
    return {k: v for k, v in ((k, sp.nsimplify(sp.expand(v))) for k, v in out.items()) if v != 0}


def perm_op(n, perm):
    """The permutation operator of site l -> perm[l], as Pauli strings, a product of SWAPs = (II+XX+YY+ZZ)/2."""
    ops, a = {'I' * n: sp.Integer(1)}, list(range(n))
    for i in range(n):
        j = a.index(perm[i])
        if j != i:
            a[i], a[j] = a[j], a[i]
            sw = {}
            for c in 'IXYZ':
                st = ['I'] * n; st[i] = c; st[j] = c; sw[''.join(st)] = sp.Rational(1, 2)
            ops = mul(ops, sw)
    return ops


def all_cuts(n):
    return [list(c) for k in range(1, n // 2 + 1) for c in itertools.combinations(range(n), k)
            if not (2 * k == n and 0 not in c)]


Q = sp.Rational

# ---- 1. the carriers ---------------------------------------------------------------------------------------
say("1. carriers: single-site Schmidt ranks, sheets beyond the jumps, blindness of the dephased seat")
tri = lambda J: {(0, 1): J[0], (0, 2): J[1], (1, 2): J[2]}
carriers = [
    ("A, colouring", 3, tri((1, 2, 3)), {0: (1, 0, 0), 1: (2, 0, 0), 2: (3, 0, 0)}),
    ("B, SWAP times letters", 3, tri((2, 2, 3)), {0: (1, 1, 0), 1: (3, 4, 0), 2: (4, 3, 0)}),
    ("C, e2 = 0", 3, tri((3, 6, -2)), {1: (3, 4, 0), 2: (5, 0, 0)}),
    ("D, d at an end, sum zero", 3, tri((1, 0, 3)), {0: (1, 0, 0), 1: (-1, -2, 0), 2: (0, 2, 0)}),
    ("E", 3, tri((1, 2, 6)), {1: (8, 0, 0), 2: (-5, 3, 0)}),
    ("E's corner", 3, tri((2, -2, 0)), {1: (1, 4, 0), 2: (4, 1, 0)}),
    ("N = 4 family a=1 e=2 h=1 g=3", 4, {(1, 2): 2, (0, 1): 1, (0, 2): -1, (1, 3): -1, (2, 3): 1},
     {1: (-1, 0, 0), 2: (1, 0, 0), 3: (0, 3, 0)}),
]
rows_1 = {}
for name, n, bonds, fields in carriers:
    sh = far_element_shape(n, heis(n, bonds, fields), ['Z' + 'I' * (n - 1)])
    assert sh is not None, name
    b = blind(n, bonds)
    cuts = {tuple(c): operator_schmidt_rank(sh.element, n, c) for c in all_cuts(n)}
    perm_products = [pm for pm in ((0,) + r for r in itertools.permutations(range(1, n)))
                     if all(operator_schmidt_rank(mul(perm_op(n, pm), sh.element), n, [l]) == 1 for l in range(n))]
    rows_1[name] = (sh.single_site_ranks, sh.sheets.extra, b, min(cuts.values()), perm_products)
    say(f"   {name:30s} ranks {sh.single_site_ranks}  every cut {cuts}  sheets {sh.sheets.extra}  blind {b}  "
        f"site permutations making it a product: {perm_products or 'none'}")
assert rows_1["A, colouring"][0] == (1, 1, 1)
assert rows_1["B, SWAP times letters"][0] == (1, 4, 4) and rows_1["B, SWAP times letters"][4] == [(0, 2, 1)]
nowhere = ("C, e2 = 0", "D, d at an end, sum zero", "E", "E's corner", "N = 4 family a=1 e=2 h=1 g=3")
assert all(rows_1[k][3] > 1 and rows_1[k][4] == [] for k in nowhere)    # no cut, no permutation makes a product
assert rows_1["D, d at an end, sum zero"][:3] == ((2, 3, 2), 3, False)
assert rows_1["E"][1:3] == (0, False) and rows_1["E's corner"][1:3] == (0, False)
assert rows_1["N = 4 family a=1 e=2 h=1 g=3"][1:3] == (2, True)
assert rows_1["C, e2 = 0"][1:3] == (1, True)
say("   the pencil of C's sheet Z_uZ_v (H_- the hop and the fields on u and v), read at lam = 2, 1/3 and 0:")
for name in ("C, e2 = 0", "E"):
    n, bonds, fields = [(c[1], c[2], c[3]) for c in carriers if c[0] == name][0]
    got = {str(lam): (lambda e: (e.near, e.far))(pencil(3, heis(3, bonds, fields), ['ZII'], 'IZZ', lam))
           for lam in (2, Q(1, 3), 0)}
    say(f"      {name}: {got}")
    assert (got['2'], got['1/3']) == (((1, 1), (1, 1)) if name.startswith("C") else ((1, 0), (1, 0)))

# ---- 2. the N = 7 carrier, exact from its closed form -------------------------------------------------------
N = 7
say("2. N = 7: W = V·W0·V^-1, W0 = SWAP(123<->456)·(X+Y)^7, V = T on (4,5,6), exact over Q(i)")
P1 = {'I': [[1, 0], [0, 1]], 'X': [[0, 1], [1, 0]], 'Y': [[0, -1j], [1j, 0]], 'Z': [[1, 0], [0, -1]]}
def g(z):
    return QQ_I(int(z.real), int(z.imag)) if isinstance(z, complex) else QQ_I(z, 0)
def kron(a, b):
    return [[a[i][j] * b[k][l] for j in range(len(a)) for l in range(len(b))] for i in range(len(a)) for k in range(len(b))]
def op(s):
    m = [[QQ_I(1, 0)]]
    for c in s:
        m = kron(m, [[g(complex(x)) if isinstance(x, complex) else QQ_I(x, 0) for x in r] for r in P1[c]])
    return m
def mat(a, b):
    n_ = len(a); bt = list(zip(*b))
    return [[sum((a[i][k] * bt[j][k] for k in range(n_) if a[i][k] and bt[j][k]), QQ_I(0, 0)) for j in range(n_)] for i in range(n_)]
def add(a, b, cb=1):
    return [[x + cb * y for x, y in zip(ra, rb)] for ra, rb in zip(a, b)]
def scal(c, a):
    return [[c * x for x in r] for r in a]
I8 = [[QQ_I(int(i == j), 0) for j in range(8)] for i in range(8)]
def q3(i, j):
    m = [[QQ_I(0, 0)] * 8 for _ in range(8)]
    for c in 'XYZ':
        s = ['I'] * 3; s[i] = c; s[j] = c; m = add(m, op(''.join(s)))
    return m
Pp = scal(QQ_I(Q(1, 6), 0), add(add(add(q3(0, 1), q3(0, 2)), q3(1, 2)), scal(QQ_I(3, 0), I8)))
Dd = scal(QQ_I(Q(1, 7), 0), add(add(scal(QQ_I(-10, 0), q3(0, 1)), scal(QQ_I(8, 0), q3(0, 2))), scal(QQ_I(2, 0), q3(1, 2))))
t = add(Pp, Dd); tinv = add(Pp, scal(QQ_I(Q(7, 144), 0), Dd))
assert mat(t, tinv) == I8
dim = 2 ** N
def expand_on_last3(m8):   # I on sites 0..3 times m8 on 4,5,6, as a sparse dict
    out = {}
    for hi in range(16):
        for a in range(8):
            for b in range(8):
                if m8[a][b]:
                    out[(hi * 8 + a, hi * 8 + b)] = m8[a][b]
    return out
Tm, Tinv = expand_on_last3(t), expand_on_last3(tinv)
# W0 = SWAP · (X+Y)^7: (X+Y) on one site is [[0, 1-i], [1+i, 0]]; its tensor power maps |b> to |b xor 1..1> with a phase
xy = {(0, 1): QQ_I(1, -1), (1, 0): QQ_I(1, 1)}
def w0_entry(col):
    row, ph = 0, QQ_I(1, 0)
    for k in range(N):
        bit = (col >> (N - 1 - k)) & 1
        ph = ph * xy[(1 - bit, bit)]
        row |= (1 - bit) << (N - 1 - k)
    bits = [(row >> (N - 1 - k)) & 1 for k in range(N)]
    for a, b in ((1, 4), (2, 5), (3, 6)):
        bits[a], bits[b] = bits[b], bits[a]
    return sum(bt << (N - 1 - k) for k, bt in enumerate(bits)), ph
W0 = {}
for col in range(dim):
    r, ph = w0_entry(col)
    W0[(r, col)] = ph
def spmul(A, B):
    byrow = {}
    for (i, k), v in B.items():
        byrow.setdefault(i, []).append((k, v))
    out = {}
    for (i, k), v in A.items():
        for (j, w) in byrow.get(k, []):
            out[(i, j)] = out.get((i, j), QQ_I(0, 0)) + v * w
    return {k: v for k, v in out.items() if v}
W = spmul(spmul(Tm, W0), Tinv)
# Pauli coefficients: c_s = 2^-N · sum_i <i| s |j> W[j, i] with s a monomial matrix (x-mask flips, z-mask signs)
letters = 'IXYZ'
coef = {}
for s in (''.join(w) for w in itertools.product(letters, repeat=N)):
    xm = sum((1 << (N - 1 - k)) for k, c in enumerate(s) if c in 'XY')
    acc = QQ_I(0, 0)
    for i in range(dim):
        j = i ^ xm
        w = W.get((j, i))
        if not w:
            continue
        ph = QQ_I(1, 0)            # <i| s |j> for s = tensor of letters, the row bit of each site is i's
        for k, c in enumerate(s):
            bi = (i >> (N - 1 - k)) & 1
            if c == 'Z' and bi: ph = -ph
            elif c == 'Y': ph = ph * (QQ_I(0, -1) if bi == 0 else QQ_I(0, 1))
        acc += ph * w
    if acc:
        coef[s] = sp.Rational(acc.x.numerator, acc.x.denominator) / dim + sp.I * sp.Rational(acc.y.numerator, acc.y.denominator) / dim
JA = {(1, 2): 1, (1, 3): 2, (2, 3): 3}; JB = {(4, 5): Q(11, 7), (4, 6): Q(22, 7), (5, 6): Q(9, 7)}
bonds7 = {(0, k): 1 for k in range(1, 7)} | JA | JB
h7 = heis(N, bonds7, {1: (1, 0, 0), 2: (1, 0, 0), 3: (1, 0, 0), 4: (0, 1, 0), 5: (0, 1, 0), 6: (0, 1, 0)})
sh7 = sheets(N, h7, ['Z' + 'I' * 6], coef)
comm = mul(h7, coef); anti = mul(coef, h7)
assert all(sp.expand(comm.get(k, 0) - anti.get(k, 0)) == 0 for k in set(comm) | set(anti))     # [H, W] = 0
zw = mul({'Z' + 'I' * 6: 1}, coef); wz = mul(coef, {'Z' + 'I' * 6: 1})
assert all(sp.expand(zw.get(k, 0) + wz.get(k, 0)) == 0 for k in set(zw) | set(wz))             # {Z0, W} = 0
r7 = tuple(operator_schmidt_rank(coef, N, [l]) for l in range(N))
cuts7 = {tuple(c): operator_schmidt_rank(coef, N, c) for c in all_cuts(N)}
say(f"   W lies in the far end: [H, W] = 0 and {{Z0, W}} = 0 exactly; rank one across {[c for c, r in cuts7.items() if r == 1]} alone "
    f"of the {len(cuts7)} cuts, the least rank across the others {min(r for r in cuts7.values() if r > 1)}")
assert [c for c, r in cuts7.items() if r == 1] == [(0,)]
def realigned_rank(Wd, site):
    """Operator Schmidt rank of the matrix Wd ({(row, col): value}) across {site} | rest, by realignment."""
    sh = N - 1 - site
    rows = {}
    for (i, j), v in Wd.items():
        a, b = (i >> sh) & 1, (j >> sh) & 1
        ri = i & ~(1 << sh); rj = j & ~(1 << sh)
        rows.setdefault((a, b), {})[(ri, rj)] = v
    keys = sorted({k for r in rows.values() for k in r})
    if not rows:
        return 0
    from sympy.polys.matrices import DomainMatrix
    M = DomainMatrix([[rows.get(rk, {}).get(k, QQ_I(0, 0)) for k in keys] for rk in sorted(rows)], (len(rows), len(keys)), QQ_I)
    return M.rank()
assert all(realigned_rank(W, l) == r for l, r in zip(range(N), r7))      # the realignment agrees with the string route
perm_hits = []
for rest in itertools.permutations(range(1, N)):
    pm = (0,) + rest                                          # site l -> pm[l]; P permutes the bits of the row index
    def prow(i):
        bits = [(i >> (N - 1 - k)) & 1 for k in range(N)]
        nb = [0] * N
        for k in range(N):
            nb[pm[k]] = bits[k]
        return sum(b << (N - 1 - k) for k, b in enumerate(nb))
    PW = {(prow(i), j): v for (i, j), v in W.items()}
    if all(realigned_rank(PW, l) == 1 for l in range(1, N)):
        perm_hits.append(pm)
say(f"   site permutations fixing the centre that make P·W a product: {perm_hits or 'none'} (all 720 read exactly)")
assert perm_hits == []
say(f"   W has {len(coef)} strings; single-site ranks {r7}; sheets beyond the jump {sh7.extra}; blind {blind(N, bonds7)}")
say("   sheet splits:", ", ".join(gs for gs, _ in sh7.splits))
assert sh7.extra == 2 and r7 == (1, 4, 4, 4, 4, 4, 4) and blind(N, bonds7)   # (X+Y) on the centre times the rest

# ---- 3. the witness's two grids ----------------------------------------------------------------------------
say("3. the triangle witness's grids: every row that pairs through one far element, by blindness, product, sheets")
letter_couplings = [(1, 2, 3), (2, 2, 3), (2, -2, 3), (3, 6, -2), (1, -2, -2), (2, 2, -1), (0, 2, 3), (2, 5, Q(-10, 7))]
letter_fields = [(0, 1, 1), (0, 1, -1), (1, 1, 1), (0, 2, 3), (2, 0, 0), (1, 2, 2)]
planar_couplings = [(1, 2, 3), (2, 2, 3), (3, 6, -2), (2, 5, Q(-10, 7)), (1, 2, 6), (2, 1, 6), (1, 3, 3), (2, 4, 12),
                    (2, -2, 0), (4, -4, 0), (0, 2, 3), (1, 2, 7), (-4, -1, -1)]
planar_pairs = [((1, 2), (2, 4)), ((0, 0), (3, 1)), ((3, 4), (5, 0)), ((3, 4), (4, 3)), ((1, 2), (3, 1)),
                ((8, 0), (-5, 3)), ((8, 0), (-5, -3)), ((Q(24, 5), Q(32, 5)), (Q(-27, 5), Q(-11, 5))), ((-5, 3), (8, 0)),
                ((16, 0), (-10, 6)), ((8, 0), (-5, 4)), ((3, 1), (3, -1)), ((1, 4), (4, 1)), ((9, 7), (9, -7)),
                ((3, 1), (3, -2)), ((6, 0), (-3, 3))]
grid = []
for w in (''.join(x) for x in itertools.product('XYZ', repeat=3)):
    for J in letter_couplings:
        for hh in letter_fields:
            grid.append((f"letter {w} J{J} h{hh}", J, {l: tuple(hh[l] if c == w[l] else 0 for c in 'XYZ') for l in range(3)}))
for J in planar_couplings:
    for fu, fv in planar_pairs:
        grid.append((f"planar J{J} fu{fu} fv{fv}", J, {1: (fu[0], fu[1], 0), 2: (fv[0], fv[1], 0)}))
tally, sheetless, seen = {}, [], set()
for name, J, f in grid:
    if sum(1 for x in J if x == 0) > 1:
        continue
    ham = heis(3, tri(J), f)
    key_h = frozenset(ham.items())
    if key_h in seen:
        continue                                                # the grids repeat rows; count each Hamiltonian once
    seen.add(key_h)
    sh = far_element_shape(3, ham, ['ZII'])
    if sh is None:
        continue
    swapped = (not sh.product) and all(operator_schmidt_rank(mul(perm_op(3, (0, 2, 1)), sh.element), 3, [l]) == 1
                                       for l in range(3))
    kind = "product" if sh.product else ("swap times a product" if swapped else "no product")
    key = (blind(3, tri(J)), kind, sh.sheets.extra > 0)
    tally[key] = tally.get(key, 0) + 1
    if sh.sheets.extra == 0:
        sheetless.append((name, J, f))
say(f"   {len(seen)} distinct connected Hamiltonians read")
for (b, p, s), c in sorted(tally.items()):
    say(f"   blind {b!s:5}  {p:20}  sheet {s!s:5}: {c}")
say("   rows without a sheet:", "; ".join(n for n, _, _ in sheetless))
assert not any(b and not s for (b, p, s) in tally)                     # every blind row has a sheet
def on_E(J, f):
    s_, t_, r_ = [sp.nsimplify(x) for x in J]
    fu = [sp.nsimplify(x) for x in f.get(1, (0, 0, 0))[:2]]; fv = [sp.nsimplify(x) for x in f.get(2, (0, 0, 0))[:2]]
    dot = fu[0] * fv[0] + fu[1] * fv[1]; e2 = s_ * t_ + s_ * r_ + t_ * r_
    r1 = r_ * (s_ - t_) ** 2 - s_ * t_ * (s_ + t_)
    r3 = (s_ + r_) * (fu[0] ** 2 + fu[1] ** 2) + (t_ + r_) * (fv[0] ** 2 + fv[1] ** 2) + 2 * (s_ + t_ + r_) * dot
    planar = all(x == 0 for x in f.get(0, (0, 0, 0))) and all(len(v) == 3 and v[2] == 0 for v in f.values())
    return planar and r1 == 0 and dot + 2 * e2 == 0 and r3 == 0 and (fu[0] * fv[1] - fu[1] * fv[0]) != 0
assert len(sheetless) == 8 and all(on_E(J, f) for _, J, f in sheetless)   # E's rows alone: R1 = R2 = R3 = 0, outside A
say("   each of them on R1 = R2 = R3 = 0, outside A")
say("ALL ASSERTIONS PASS")
Path(__file__).parent.joinpath("results", "far_element_sheets_census.txt").write_bytes(("\n".join(OUT) + "\n").encode("utf-8"))
