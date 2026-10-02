"""Independent gate for n3_palindrome_covering.py. It shares no code with the covering: the two ends are built
here from dense 8x8 Pauli matrices with complex coefficients over Q (real and imaginary parts as separate
unknowns), as the kernels of W -> [H, W] restricted to d-letter {I, Z} (near) and {X, Y} (far), ranked exactly
with Fraction elimination; the classes A, B, C are evaluated directly on the point (B on the same four axes
x, y, x + y, x - y the covering uses, which suffice for letter fields).

Four checks:
 1. on-leaf: for every far leaf of the certificate, exact rational points on its region (free variables drawn,
    the rest solved from a lex Groebner basis); the leaf's stored far nullity must equal the dense one, and
    every connected palindromic point found must lie in A, B or C;
 2. special loci over all 27 letter words (no symmetry reduction): equal, opposite and reciprocal couplings,
    one coupling zero, equal, opposite and zero fields; every connected palindromic point in A, B or C, and
    every connected point of A, B or C palindromic (the converse); and each connected point, moved by the
    symmetries to its word's orbit representative and scaled into its chart, lies in a far leaf of the
    certificate whose stored far count equals the dense one (the leaves cover);
 3. targeted points of B and of C (all letters, signs, several rational e2 = 0 couplings): every one pairs, and
    with C left out of the classes the C points are reported outside, so the class check can fail;
 4. a planted mutation: dropping the J3 ZZ term from the dense H must break check 1 somewhere.
Check 1 also requires at least one leaf inside each of A, B and C to be reached at a rational point.
Run:  python simulations/n3_palindrome_covering_gate.py [workers=20] [budget_seconds=1800]"""
import itertools, json, os, random, sys
from fractions import Fraction as Fr
import sympy as sp

CERT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'n3_palindrome_covering.json')
J1, J2, J3, h0, h1, h2 = PARAMS = sp.symbols('J1 J2 J3 h0 h1 h2')

# Gaussian rationals as (re, im) Fraction pairs
P1 = {'I': [[(1, 0), (0, 0)], [(0, 0), (1, 0)]], 'X': [[(0, 0), (1, 0)], [(1, 0), (0, 0)]],
      'Y': [[(0, 0), (0, -1)], [(0, 1), (0, 0)]], 'Z': [[(1, 0), (0, 0)], [(0, 0), (-1, 0)]]}

def cmul(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])

def kron(A, B):
    n, m = len(A), len(B)
    return [[cmul(A[i // m][j // m], B[i % m][j % m]) for j in range(n * m)] for i in range(n * m)]

def string(w):
    M = P1[w[0]]
    for c in w[1:]: M = kron(M, P1[c])
    return M

STR = {w: string(w) for w in (''.join(t) for t in itertools.product('IXYZ', repeat=3))}

def hamiltonian(pat, v, drop_zz3=False):
    H = [[(Fr(0), Fr(0))] * 8 for _ in range(8)]
    def add(w, c):
        S = STR[w]
        for i in range(8):
            for j in range(8):
                if S[i][j] != (0, 0):
                    H[i][j] = (H[i][j][0] + c * S[i][j][0], H[i][j][1] + c * S[i][j][1])
    for (a, b, J) in ((0, 1, v['J1']), (0, 2, v['J2']), (1, 2, v['J3'])):
        for L in 'XYZ':
            if drop_zz3 and (a, b) == (1, 2) and L == 'Z': continue
            w = ['I'] * 3; w[a] = L; w[b] = L; add(''.join(w), J)
    for l, hh in enumerate((v['h0'], v['h1'], v['h2'])):
        w = ['I'] * 3; w[l] = pat[l]; add(''.join(w), hh)
    return H

def matmul(A, B):
    return [[(sum(A[i][k][0] * B[k][j][0] - A[i][k][1] * B[k][j][1] for k in range(8)),
              sum(A[i][k][0] * B[k][j][1] + A[i][k][1] * B[k][j][0] for k in range(8))) for j in range(8)] for i in range(8)]

def nullity(rows, n):
    M = [r[:] for r in rows]; r = 0
    for c in range(n):
        piv = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c] / M[r][c]; M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        r += 1
    return n - r

def end(H, far):
    """dim over C of {W in span of strings with d-letter in {X,Y} (far) or {I,Z} (near) : HW = WH}."""
    basis = [w for w in STR if (w[0] in 'XY') == far]
    comms = []
    for w in basis:
        S = STR[w]; A = matmul(H, S); B = matmul(S, H)
        comms.append([(A[i][j][0] - B[i][j][0], A[i][j][1] - B[i][j][1]) for i in range(8) for j in range(8)])
    n = len(basis); rows = []
    for e in range(64):                       # unknowns x_k + i y_k: real part and imaginary part rows
        rows.append([comms[k][e][0] for k in range(n)] + [-comms[k][e][1] for k in range(n)])
        rows.append([comms[k][e][1] for k in range(n)] + [comms[k][e][0] for k in range(n)])
    return nullity(rows, 2 * n) // 2

def in_class(pat, v, without=None):
    """A, B or C evaluated at the point (independent of the covering's generator lists); `without` names a class
    to leave out, for the planted mutation."""
    f = [tuple(v[h] * (1 if pat[l] == L else 0) for L in 'XYZ') for l, h in enumerate(('h0', 'h1', 'h2'))]
    for p in (0, 1):                                                       # A: every field along X, or along Y
        if without != 'A' and all(fl[q] == 0 for fl in f for q in range(3) if q != p): return 'A'
    rots = [lambda x: (x[0], -x[1], -x[2]), lambda x: (-x[0], x[1], -x[2]),
            lambda x: (x[1], x[0], -x[2]), lambda x: (-x[1], -x[0], -x[2])]
    if without != 'B' and v['J1'] == v['J2'] and any(R(f[0]) == f[0] and R(f[1]) == f[2] for R in rots): return 'B'
    e2 = v['J1'] * v['J2'] + v['J1'] * v['J3'] + v['J2'] * v['J3']
    if without != 'C' and e2 == 0 and f[0] == (0, 0, 0) and f[1][2] == 0 and f[2][2] == 0 and \
            sum(c * c for c in f[1]) == sum(c * c for c in f[2]): return 'C'
    return None

def connected(v):
    return sum(1 for k in ('J1', 'J2', 'J3') if v[k] == 0) <= 1

def check_point(pat, v, drop=False):
    H = hamiltonian(pat, v, drop)
    nr, fr = end(H, False), end(H, True)
    pal = fr >= 1 and nr == fr
    bad = pal and connected(v) and not in_class(pat, v)
    return nr, fr, pal, bad

# the symmetries of step 2: the quarter turn about z (X -> Y, Y -> -X on every site) and the swap u <-> v
def _quarter(word, v):
    w, out = '', dict(v)
    for l, (c, h) in enumerate(zip(word, ('h0', 'h1', 'h2'))):
        if c == 'X': w += 'Y'
        elif c == 'Y': w += 'X'; out[h] = -v[h]
        else: w += c
    return w, out

def _swap(word, v):
    out = dict(v)
    out['J1'], out['J2'], out['h1'], out['h2'] = v['J2'], v['J1'], v['h2'], v['h1']
    return word[0] + word[2] + word[1], out

def to_chart(word, v, charts, reps):
    """The point moved to its orbit representative and scaled into the first chart of that word containing it."""
    cands = [(word, v)]
    for _ in range(3):
        cands += [_quarter(w, x) for w, x in cands] + [_swap(w, x) for w, x in cands]
    for w, x in cands:
        if w not in reps: continue
        a, b = ('J2', 'J1') if charts[w][0].startswith('J2') else ('J1', 'J2')
        if x[a] != 0:
            return w, charts[w][0], {k: val / x[a] for k, val in x.items()}
        if x[b] != 0:
            return w, charts[w][1], {k: val / x[b] for k, val in x.items()}
    return None

def _exact(expr):
    """A polynomial as a function of the six parameters evaluated in Fraction arithmetic, coefficients included."""
    terms = [(tuple(mon), Fr(int(c.p), int(c.q))) for mon, c in sp.Poly(sp.sympify(expr), *PARAMS).terms()]
    def f(*x):
        total = Fr(0)
        for mon, c in terms:
            t = c
            for xi, d in zip(x, mon):
                if d: t *= xi ** d
            total += t
        return total
    return f

def compile_leaf(lf):
    """Equations and recorded-nonzero factors as exact functions of the six parameters."""
    return [_exact(e) for e in lf['eqs']], [_exact(n) for n in lf['neqs']]

def leaf_points(eqs, neqs, rng, tries=30):
    """Exact rational points on {eqs = 0, neqs != 0}: draw the free variables of a lex basis, solve the rest."""
    G = sp.groebner([sp.sympify(e) for e in eqs], *PARAMS, order='lex')
    if list(G.exprs) == [1]: return []
    pts = []
    for _ in range(tries):
        v = {}
        ok = True
        for x in reversed(PARAMS):                     # lex: the last variable is eliminated last
            polys = [sp.expand(g.subs(v)) for g in G.exprs]
            polys = [q for q in polys if q != 0 and q.free_symbols <= {x}]
            if polys:
                g = polys[0]
                for q in polys[1:]: g = sp.gcd(g, q)
                rts = [r for r in sp.Poly(g, x).ground_roots() if r.is_rational] if sp.Poly(g, x).degree() > 0 else []
                if not rts: ok = False; break
                v[x] = rng.choice(rts)
            else:
                v[x] = sp.Rational(rng.randint(-7, 7), rng.randint(1, 3))
        if not ok: continue
        if any(sp.sympify(e).subs(v) != 0 for e in eqs): continue
        if any(sp.sympify(n).subs(v) == 0 for n in neqs): continue
        pts.append({str(k): Fr(int(val.p), int(val.q)) for k, val in v.items()})
        if len(pts) >= 3: break
    return pts

def leaf_job(job):
    """One far leaf: its rational points, the dense far nullity there, palindromic points against the classes,
    and the J3-ZZ mutation."""
    pat, chart, idx, lf, seed = job
    rng = random.Random(seed)
    pts = leaf_points(lf['eqs'], lf['neqs'], rng)
    why = lf.get('decided_on_far_region', '')
    kind = why[len('inside ')] if why.startswith('inside ') else None
    out = {'tested': bool(pts), 'fails': [], 'pal': 0, 'mutant': False, 'class': kind}
    for v in pts:
        nr, fr, pal, bad = check_point(pat, v)
        out['pal'] += pal
        if fr != lf['far'] or bad:
            out['fails'].append((pat, chart, idx, {k: str(x) for k, x in v.items()}, lf['far'], nr, fr))
        if lf['far'] > 0 and check_point(pat, v, drop=True)[1] != lf['far']:
            out['mutant'] = True
    return out

def main(workers=20, budget=1800):
    import time
    from multiprocessing import Pool, TimeoutError
    rng = random.Random(5)
    fails = 0
    jobs = []
    cert = json.load(open(CERT))
    if cert['summary']['open'] != 0: fails += 1; print('FAIL the certificate reports open regions', flush=True)
    jobs = [(lf['pattern'], lf['chart'], lf['leaf'], lf, rng.randrange(10 ** 9)) for lf in cert['far_leaves']]
    # 1. on-leaf, every far leaf, a shared wall-clock budget; a leaf without a result in it counts as untested
    tested = untested = pal_pts = 0
    mutant = False
    reached = {'A': 0, 'B': 0, 'C': 0}
    inside = {k: sum(1 for lf in cert['far_leaves'] if lf.get('decided_on_far_region', '').startswith('inside ' + k))
              for k in reached}
    pool = Pool(workers)
    handles = [pool.apply_async(leaf_job, (j,)) for j in jobs]
    t_end = time.time() + budget
    for j, h in zip(jobs, handles):
        try:
            left = t_end - time.time()
            if left <= 0 and not h.ready(): raise TimeoutError
            r = h.get(timeout=max(left, 0.01))
        except TimeoutError:
            untested += 1; continue
        if not r['tested']: untested += 1; continue
        tested += 1; pal_pts += r['pal']; mutant |= r['mutant']
        if r['class']: reached[r['class']] += 1
        for f in r['fails']:
            fails += 1; print('FAIL on-leaf', f, flush=True)
    pool.terminate()
    print(f'1. on-leaf: {len(jobs)} far leaves, {tested} tested at exact rational points, {untested} without a point '
          f'(none found, or not done within {budget} s), {pal_pts} palindromic points', flush=True)
    print('   leaves inside a class reached at a rational point: ' +
          ', '.join(f'{k} {reached[k]} of {inside[k]}' for k in reached), flush=True)
    for k in reached:                                  # every class must be met at least once on its own leaves
        if reached[k] == 0: fails += 1; print(f'FAIL no leaf inside {k} reached', flush=True)
    # 2. special loci over all 27 words; each point is also located in the certificate: moved to its orbit
    # representative and scaled into its chart, some far leaf must contain it and carry its dense far count
    vals = [Fr(-2), Fr(-1), Fr(1), Fr(2), Fr(3), Fr(1, 2)]
    n2 = pal2 = located = unscaled_miss = 0
    charts = cert['charts']; reps = set(charts)
    by_chart = {}
    for lf in cert['far_leaves']:
        by_chart.setdefault((lf['pattern'], lf['chart']), []).append((lf, None))
    for pat in (''.join(t) for t in itertools.product('XYZ', repeat=3)):
        for _ in range(60):
            a, b, c = (rng.choice(vals) for _ in range(3))
            kind = rng.randrange(6)
            J = [a, b, c]
            if kind == 1: J[1] = J[0]
            elif kind == 2 and a + b != 0: J[2] = -a * b / (a + b)
            elif kind == 3: J[rng.randrange(3)] = Fr(0)
            elif kind == 4: J[1] = -J[0]
            s = rng.choice(vals)
            hs = [rng.choice([Fr(0), s, -s, rng.choice(vals)]) for _ in range(3)]
            if kind == 5: hs[0] = Fr(0); hs[2] = rng.choice([s, -s]); hs[1] = s
            v = dict(zip(('J1', 'J2', 'J3', 'h0', 'h1', 'h2'), J + hs))
            nr, fr, pal, bad = check_point(pat, v)
            n2 += 1; pal2 += pal
            if bad: fails += 1; print('FAIL special', pat, v, nr, fr, flush=True)
            if connected(v):
                w, ch, x = to_chart(pat, v, charts, reps)
                args = [x[k] for k in ('J1', 'J2', 'J3', 'h0', 'h1', 'h2')]
                hit = []
                lst = by_chart[(w, ch)]
                for i, (lf, comp) in enumerate(lst):
                    if comp is None: comp = compile_leaf(lf); lst[i] = (lf, comp)
                    if all(f(*args) == 0 for f in comp[0]) and all(f(*args) != 0 for f in comp[1]): hit.append(lf['far'])
                if not hit or any(k != fr for k in hit):
                    fails += 1; print('FAIL cover', pat, v, w, ch, 'leaves', hit, 'dense far', fr, flush=True)
                else: located += 1
                # mutation: the same point without the scaling into the chart must miss the leaves somewhere
                raw = [v[k] for k in ('J1', 'J2', 'J3', 'h0', 'h1', 'h2')]
                if raw != args and not any(all(f(*raw) == 0 for f in comp[0]) and all(f(*raw) != 0 for f in comp[1])
                                           for _, comp in lst if comp is not None):
                    unscaled_miss += 1
            if connected(v) and in_class(pat, v) and not pal:      # the converse: every class point pairs
                fails += 1; print('FAIL converse', pat, v, in_class(pat, v), nr, fr, flush=True)
    print(f'2. special loci: {n2} points over the 27 words, {pal2} palindromic; every connected palindromic point '
          f'in A, B or C and every connected point of A, B or C palindromic, unless FAIL lines above; {located} '
          f'connected points located in a far leaf of the certificate carrying their far count; unscaled, '
          f'{unscaled_miss} of them miss every leaf, so the location can fail', flush=True)
    if unscaled_miss == 0: fails += 1; print('FAIL the cover check cannot fail', flush=True)
    # 3. targeted B and C points: every one pairs; with C left out of the classes, the C points must be flagged
    nb = nc = flagged = 0
    for (a, b) in [(3, 6), (1, -2), (-1, 2), (2, 5), (Fr(1, 2), 3), (-3, 7)]:
        a, b = Fr(a), Fr(b)
        r = -a * b / (a + b)
        for d in 'XYZ':
            for pu, pv in (('X', 'Y'), ('Y', 'X')):
                for sg in (1, -1):
                    for h in (Fr(1), Fr(5, 2)):
                        v = dict(J1=a, J2=b, J3=r, h0=Fr(0), h1=h, h2=sg * h)
                        pat = d + pu + pv
                        nr, fr, pal, bad = check_point(pat, v)
                        nc += 1
                        if in_class(pat, v) != 'C' and a != b: fails += 1; print('FAIL not C', pat, v, flush=True)
                        if not pal or bad: fails += 1; print('FAIL C point', pat, v, nr, fr, flush=True)
                        if pal and not in_class(pat, v, without='C'): flagged += 1
    for _ in range(3000):
        a = rng.choice(vals); c = rng.choice(vals + [Fr(0)])
        pat = ''.join(rng.choice('XYZ') for _ in range(3))
        hs = [rng.choice([Fr(0), Fr(1), Fr(-1), Fr(2)]) for _ in range(3)]
        v = dict(zip(('J1', 'J2', 'J3', 'h0', 'h1', 'h2'), [a, a, c] + hs))
        if in_class(pat, v) != 'B': continue
        nr, fr, pal, bad = check_point(pat, v)
        nb += 1
        if not pal: fails += 1; print('FAIL B point', pat, v, nr, fr, flush=True)
    print(f'3. targeted: {nc} C points and {nb} B points, all palindromic unless FAIL lines above; with C left out '
          f'of the classes {flagged} of the C points would be reported outside A, B and C', flush=True)
    if flagged == 0 or nb == 0: fails += 1; print('FAIL the class check cannot fail on C, or no B point', flush=True)
    # 4. mutation of the physics
    print(f'4. mutation (the J3 ZZ term dropped) changes a stored far nullity: {"yes" if mutant else "NO"}', flush=True)
    if not mutant: fails += 1
    print('GATE', 'PASS' if fails == 0 else f'FAIL ({fails})', flush=True)
    sys.exit(1 if fails else 0)

if __name__ == '__main__':
    main(*(int(a) for a in sys.argv[1:3]))
