"""A palindrome at N = 4 with neither a colouring nor a site symmetry, exact, and the blind seat behind it.

Family: Heisenberg bonds on sites 0..3 with the jump Z on site 0 alone,
    H = e·s1·s2 + a·(s0 - s3)·(s1 - s2) - h·X1 + h·X2 + g·Y3
(J01 = a, J02 = -a, J13 = -a, J23 = a, J12 = e, J03 = 0; s_l the vector of Pauli matrices on site l).
Part 1, exact elimination in (a, e, h, g): the two ends of F158 (the real kernels of W -> [H, W] on the strings
whose site-0 letter is in {I, Z}, near, or {X, Y}, far) by Gaussian elimination with case splits; the generic
leaf has near = far = 1 on every sub-leaf, and so has the leaf h = 0 (where g·Y3 is the only field and Y^4 colours
it); with a != 0 and g != 0 the only other leaf, 4a² + h² = 0, has no real point. So every real row with a != 0,
g != 0 pairs, for every e and h, and with h != 0 it is no colouring.
Part 2, independent gate: both ends from dense 16 x 16 Pauli matrices in exact Fraction arithmetic at rational
points of the family (they pair) and with J03 != 0 or J23 detuned (they do not), sharing no code with part 1.
Part 3, the blind seat: on the triangle the Krylov determinant of the dephased seat d under the signed one-magnon
Laplacian is -(s - t)·e2 (the factor of the common-carrier proposition), and on the family seat 0 is blind.
Run: python simulations/n4_pair_channel_family.py"""
import itertools, sys
from fractions import Fraction as Fr
import sympy as sp
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
from sympy.polys.groebnertools import groebner as _groebner
from sympy.polys.orderings import grevlex

a, e, h, g = SYMS = sp.symbols('a e h g')
RING, *_ = ring('a e h g', QQ, grevlex)
COUPLINGS = {(0, 1): a, (0, 2): -a, (1, 2): e, (1, 3): -a, (2, 3): a}
FIELDS = ((1, 'X', -h), (2, 'X', h), (3, 'Y', g))

# ---- part 1: exact elimination ----------------------------------------------------------------------------
def _mul1(p, q):
    if p == 0: return 0, q
    if q == 0: return 0, p
    if p == q: return 0, 0
    return (1 if (p, q) in ((1, 2), (2, 3), (3, 1)) else 3), 6 - p - q

def _mul(s, t):
    ph, out = 0, []
    for p, q in zip(s, t):
        k, c = _mul1(p, q); ph += k; out.append(c)
    return ph % 4, tuple(out)

def _anti(s, t): return sum(1 for p, q in zip(s, t) if p and q and p != q) % 2 == 1

def system(far):
    cols = [s for s in itertools.product(range(4), repeat=4) if (s[0] in (1, 2)) == far]
    terms = []
    for (p, q), w in COUPLINGS.items():
        for c in (1, 2, 3):
            s = [0] * 4; s[p] = c; s[q] = c; terms.append((w, tuple(s)))
    for l, letter, w in FIELDS:
        s = [0] * 4; s[l] = 'XYZ'.index(letter) + 1; terms.append((w, tuple(s)))
    rows = {}
    for w, t in terms:
        for j, s in enumerate(cols):
            if _anti(t, s):
                k, u = _mul(t, s)
                rows.setdefault(u, {}); rows[u][j] = sp.expand(rows[u].get(j, 0) + (1 if k == 1 else -1) * w)
    return [{j: RING.from_expr(v) for j, v in r.items() if v != 0} for r in rows.values()], len(cols)

_GB, _FAC = {}, {}
def _gb(eqs):
    key = tuple(sorted(str(x) for x in eqs))
    if key not in _GB: _GB[key] = _groebner(list(eqs), RING) if eqs else []
    return _GB[key]
def _factors(v):
    if v not in _FAC: _FAC[v] = v.factor_list()
    return _FAC[v]
def _nf(p, G): return p.rem(G) if G else p

def eliminate(rows, ncols, eqs, neqs, leaves):
    """Gaussian elimination with case splits: a pivot is a constant or has every factor recorded nonzero;
    otherwise each factor's vanishing and the all-nonzero case are separate branches; rows are divided only by
    recorded-nonzero factors; a branch is dropped only when 1 or a recorded factor lies in its ideal."""
    G = _gb(eqs)
    if any(x.is_ground and x != 0 for x in G): return
    if any(_nf(n, G) == 0 for n in neqs): return
    rows = [{j: _nf(v, G) for j, v in r.items()} for r in rows]
    rows = [r for r in ({j: v for j, v in r.items() if v != 0} for r in rows) if r]
    keys = {n.monic() for n in neqs}
    while rows:
        piv = next(((i, j) for i, r in enumerate(rows) for j, v in r.items() if v.is_ground), None)
        if piv is None:
            piv = next(((i, j) for i, r in enumerate(rows) for j, v in r.items()
                        if all(f.monic() in keys for f, _ in _factors(v)[1])), None)
        if piv is None:
            v = min(((max(sum(m) for m in v.monoms()), len(v), v) for r in rows for v in r.values()),
                    key=lambda c: (c[0], c[1]))[2]
            facs = [f for f, _ in _factors(v)[1] if not f.is_ground]
            for f in facs: eliminate([dict(r) for r in rows], ncols, eqs + [f], neqs, leaves)
            eliminate([dict(r) for r in rows], ncols, eqs, neqs + facs, leaves)
            return
        i, j = piv
        pr = rows.pop(i); p = pr[j]; new = []
        for r in rows:
            if j not in r: new.append(r); continue
            x = r[j]; nr = {}
            for k in set(r) | set(pr):
                val = _nf(p * r.get(k, RING.zero) - x * pr.get(k, RING.zero), G)
                if val != 0: nr[k] = val
            if nr:
                vals = list(nr.values()); d = vals[0]
                for v2 in vals[1:]:
                    d = d.gcd(v2)
                    if d.is_ground: break
                keep = RING.one
                if not d.is_ground:
                    for fac, mult in _factors(d)[1]:
                        if fac.monic() in keys: keep *= fac ** mult
                if keep != 1: nr = {k: v2.exquo(keep) for k, v2 in nr.items()}
                new.append(nr)
        rows = new; ncols -= 1
    leaves.append((list(eqs), list(neqs), ncols))

# ---- part 2: independent dense gate ------------------------------------------------------------------------
_P = {'I': ((1, 0), (0, 0), (0, 0), (1, 0)), 'X': ((0, 0), (1, 0), (1, 0), (0, 0)),
      'Y': ((0, 0), (0, -1), (0, 1), (0, 0)), 'Z': ((1, 0), (0, 0), (0, 0), (-1, 0))}

def _string(word):
    m = [[(Fr(1), Fr(0))]]
    for c in word:
        q = _P[c]; n = len(m); out = [[None] * (2 * n) for _ in range(2 * n)]
        for i in range(n):
            for j in range(n):
                x = m[i][j]
                for bi in range(2):
                    for bj in range(2):
                        y = q[2 * bi + bj]
                        out[2 * i + bi][2 * j + bj] = (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])
        m = out
    return m

STR = {w: _string(w) for w in (''.join(t) for t in itertools.product('IXYZ', repeat=4))}

def _rank(M, ncols):
    M = [r[:] for r in M if any(r)]; r = 0
    for c in range(ncols):
        piv = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c] / M[r][c]; M[i] = [u - f * v for u, v in zip(M[i], M[r])]
        r += 1
    return r

def dense_ends(J, fields, far_rows_only=False):
    """J: {(p, q): weight}, fields: [(site, letter, weight)]; exact near and far over Q(i) (or, with
    far_rows_only, the real rows of the far system, for stacking several rows' conditions)."""
    H = [[(Fr(0), Fr(0)) for _ in range(16)] for _ in range(16)]
    def add(w, c):
        S = STR[w]
        for i in range(16):
            for j in range(16):
                if S[i][j] != (0, 0): H[i][j] = (H[i][j][0] + c * S[i][j][0], H[i][j][1] + c * S[i][j][1])
    for (p, q), c in J.items():
        for L in 'XYZ':
            w = ['I'] * 4; w[p] = L; w[q] = L; add(''.join(w), Fr(c))
    for l, L, c in fields:
        w = ['I'] * 4; w[l] = L; add(''.join(w), Fr(c))
    def mm(A, B):
        return [[(sum(A[i][k][0] * B[k][j][0] - A[i][k][1] * B[k][j][1] for k in range(16)),
                  sum(A[i][k][0] * B[k][j][1] + A[i][k][1] * B[k][j][0] for k in range(16))) for j in range(16)] for i in range(16)]
    out = []
    for far in ((True,) if far_rows_only else (False, True)):
        basis = [w for w in STR if (w[0] in 'XY') == far]
        comms = []
        for w in basis:
            A = mm(H, STR[w]); B = mm(STR[w], H)
            comms.append([(A[i][j][0] - B[i][j][0], A[i][j][1] - B[i][j][1]) for i in range(16) for j in range(16)])
        n = len(basis); rows = []
        for x in range(256):
            rows.append([comms[k][x][0] for k in range(n)] + [-comms[k][x][1] for k in range(n)])
            rows.append([comms[k][x][1] for k in range(n)] + [comms[k][x][0] for k in range(n)])
        if far_rows_only: return rows, 2 * n
        out.append((2 * n - _rank(rows, 2 * n)) // 2)
    return tuple(out)

# ---- part 3: the blind seat --------------------------------------------------------------------------------
def laplacian(n, J):
    L = sp.zeros(n, n)
    for (p, q), w in J.items():
        L[p, p] += w; L[q, q] += w; L[p, q] -= w; L[q, p] -= w
    return L

def krylov_det(n, J, seat):
    L = laplacian(n, J); v = sp.zeros(n, 1); v[seat] = 1
    return sp.factor(sp.Matrix.hstack(*[L ** k * v for k in range(n)]).det())

def main():
    fails = 0
    far_rows, nf = system(True); near_rows, nn = system(False)
    leaves = []; eliminate(far_rows, nf, [], [], leaves)
    generic = [lf for lf in leaves if not lf[0]]
    print('1. far leaves:', len(leaves))
    for E, N, k in leaves:
        nl = []; eliminate(near_rows, nn, E, N, nl)
        print('   far %3d | =0: %s | !=0: %s | near on it: %s' % (
            k, [str(x.as_expr()) for x in E], sorted({str(x.as_expr()) for x in N}), sorted({k2 for _, _, k2 in nl})))
    if len(generic) != 1 or generic[0][2] != 1:
        fails += 1; print('FAIL the generic far leaf is not one-dimensional')
    else:
        nl = []; eliminate(near_rows, nn, [], generic[0][1], nl)
        nz = sorted({str(x.as_expr()) for x in generic[0][1]})
        if {k for _, _, k in nl} != {1}: fails += 1; print('FAIL near on the generic leaf is not 1 everywhere')
        else: print('   generic leaf (only %s nonzero): far 1, near 1, so it pairs' % ', '.join(nz))
    # every leaf on which neither a nor g is forced to vanish is one of exactly three: generic, h = 0, 4a^2 + h^2 = 0
    live = [[str(x.as_expr()) for x in E] for E, N, k in leaves
            if not any(str(x.as_expr()) in ('a', 'g') for x in E)]
    if sorted(map(tuple, live)) != sorted([(), ('h',), ('4*a**2 + h**2',)]):
        fails += 1; print('FAIL the leaves with a and g free are', live)
    else:
        sq = [lf for lf in leaves if [str(x.as_expr()) for x in lf[0]] == ['4*a**2 + h**2']][0]
        if 'a' not in {str(x.as_expr()) for x in sq[1]}: fails += 1; print('FAIL a is not recorded nonzero on 4a^2 + h^2 = 0')
        else: print('   leaves with a and g free: exactly the generic one, h = 0, and 4a² + h² = 0 (a recorded nonzero: no real point)')
    h0 = [lf for lf in leaves if [str(x.as_expr()) for x in lf[0]] == ['h']]
    if len(h0) != 1 or h0[0][2] != 1: fails += 1; print('FAIL the leaf h = 0 is not one leaf of far count 1')
    else:
        nl = []; eliminate(near_rows, nn, h0[0][0], h0[0][1], nl)
        if {k for _, _, k in nl} != {1}: fails += 1; print('FAIL near on the leaf h = 0 is not 1 everywhere')
        else: print('   leaf h = 0 (a, g nonzero): far 1, near 1; there the only field is g·Y3 and Y^4 colours it')
    # 2. the gate
    def J_of(A, Ee, c=0, d23=0):
        return {(0, 1): A, (0, 2): -A, (1, 2): Ee, (1, 3): -A, (2, 3): A + d23, (0, 3): c}
    pts = [(1, 2, 3, 5), (-2, 0, 1, -1), (3, -1, Fr(1, 2), 2), (Fr(-1, 3), 4, 2, Fr(7, 2)), (5, 1, Fr(-2, 5), 1)]
    for A, Ee, H_, G_ in pts:
        nr, fr = dense_ends(J_of(A, Ee), [(1, 'X', -H_), (2, 'X', H_), (3, 'Y', G_)])
        ok = nr == fr >= 1
        print('2. a=%s e=%s h=%s g=%s: near %d far %d %s' % (A, Ee, H_, G_, nr, fr, 'pairs' if ok else 'BREAKS'))
        if not ok: fails += 1
    for A, Ee, H_, G_, c, d in [(1, 2, 3, 5, 1, 0), (1, 2, 3, 5, 0, 1), (-2, 0, 1, -1, 2, 0)]:
        nr, fr = dense_ends(J_of(A, Ee, c, d), [(1, 'X', -H_), (2, 'X', H_), (3, 'Y', G_)])
        print('   control J03=%s, J23 detuned by %s: near %d far %d %s' % (c, d, nr, fr, 'breaks' if nr != fr else 'PAIRS'))
        if nr == fr: fails += 1
    # the pair's fields must be opposite: unequal magnitudes on sites 1 and 2 break it
    for A, Ee, H_, G_, H2 in [(1, 2, 3, 5, 2), (-2, 0, 1, -1, -1)]:
        nr, fr = dense_ends(J_of(A, Ee), [(1, 'X', -H_), (2, 'X', H2), (3, 'Y', G_)])
        print('   control field on site 2 = %s instead of %s: near %d far %d %s' % (H2, H_, nr, fr, 'breaks' if nr != fr else 'PAIRS'))
        if nr == fr: fails += 1
    # the carrier depends on h: on the line g = h at a = -1, e = 2, far is 1 at h = 1 and at h = 2, and 0 for both at once
    r1, n1 = dense_ends(J_of(-1, 2), [(1, 'X', -1), (2, 'X', 1), (3, 'Y', 1)], far_rows_only=True)
    r2, _ = dense_ends(J_of(-1, 2), [(1, 'X', -2), (2, 'X', 2), (3, 'Y', 2)], far_rows_only=True)
    f1, f2, f12 = [(n1 - _rank(r, n1)) // 2 for r in (r1, r2, r1 + r2)]
    print('   carrier on the line g = h at a = -1, e = 2: far %d at h = 1, %d at h = 2, %d for both at once' % (f1, f2, f12))
    if (f1, f2, f12) != (1, 1, 0): fails += 1
    # the crossing holds for every J03: an S = 1 triplet at the S = 2 quintet's level, multiplicity 8 at e + J03
    xs = sp.Symbol('x')
    for A, Ee, c in [(1, 2, 0), (1, 2, 1), (-1, 2, Fr(7, 10)), (2, -1, 3)]:
        Hb = sp.zeros(16, 16)
        for (p, q), w in J_of(A, Ee, c).items():
            for L in 'XYZ':
                wd = ['I'] * 4; wd[p] = L; wd[q] = L; S = STR[''.join(wd)]
                Hb += sp.Rational(w) * sp.Matrix(16, 16, lambda i, j: sp.Rational(S[i][j][0]) + sp.I * sp.Rational(S[i][j][1]))
        lev = sp.Rational(Ee) + sp.Rational(c)
        mult = 16 - (Hb - lev * sp.eye(16)).rank()
        print('   J03 = %s: level e + J03 = %s has multiplicity %d (quintet 5 + triplet 3)' % (c, lev, mult))
        if mult != 8: fails += 1
    # 3. the blind seat
    s_, t_, r_ = sp.symbols('s t r')
    tri = krylov_det(3, {(0, 1): s_, (0, 2): t_, (1, 2): r_}, 0)
    print('3. triangle, Krylov determinant of the dephased seat:', tri)
    if sp.expand(tri + (s_ - t_) * (s_ * t_ + s_ * r_ + t_ * r_)) != 0: fails += 1; print('FAIL triangle identity')
    fam = krylov_det(4, {(0, 1): a, (0, 2): -a, (1, 2): e, (1, 3): -a, (2, 3): a}, 0)
    print('   family, Krylov determinant of seat 0:', fam, '(identically zero: seat 0 blind)')
    if fam != 0: fails += 1
    print('ALL PASS' if fails == 0 else f'FAIL ({fails})')
    sys.exit(1 if fails else 0)

if __name__ == '__main__':
    main()
