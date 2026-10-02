"""The palindromic locus of the N = 3 Heisenberg triangle with one dephased site, covered exactly.

Sites (d, u, v) = (0, 1, 2), bonds J1 = J_du, J2 = J_dv, J3 = J_uv real (connected: at most one zero), the
jump Z on d alone, u and v undephased, a single-letter field of free real magnitude h0, h1, h2 on each site.
The two ends of F158 are the real solution spaces of [H, W] = 0 over the Pauli strings whose d-letter is in
{I, Z} (near) or {X, Y} (far); the row is palindromic iff near = far >= 1.

Every entry of the linear system is linear and homogeneous in (J1, J2, J3, h0, h1, h2), so nullities are
invariant under a common nonzero scaling. A connected point has J1 != 0 or (J1 = 0 and J2 != 0), so the two
charts J1 = 1 and (J1 = 0, J2 = 1) cover every connected point up to scaling; so do J2 = 1 and (J2 = 0, J1 = 1).

Per field pattern (one orbit representative of the 27 letter words under the two symmetries that fix the
jump: the swap u <-> v with J1 <-> J2, and the global quarter turn about z, X <-> Y; a site without a field is
h = 0 inside any word) and chart: Gaussian elimination of the far system with exact case splits gives far
leaves (equalities, inequations, exact nullity k). Each leaf with k >= 1 is decided either on the far region
itself (inside a class, disconnected, or without a real point) or by eliminating the near system on it; every
near region with nullity k is palindromic and is classified against
  A  a colouring: some letter p in {X, Y} carries every nonzero field,
  B  J1 = J2 and a pi rotation R_n about an axis n perpendicular to z with R_n f_d = f_d, R_n f_u = f_v,
  C  e2 = J1J2 + J1J3 + J2J3 = 0, no field on d, f_u and f_v perpendicular to z, |f_u| = |f_v|,
by Rabinowitsch (every generator of the class vanishes on the region). A region in none is reported OPEN.

The certificate (every far leaf, its decision, every palindromic region and the class it fell into) is
written to simulations/results/n3_palindrome_covering.json; per-leaf work files go to the local-only
simulations/_n3_palindrome_covering_work/, so an interrupted run resumes. Run:
  python simulations/n3_palindrome_covering.py [workers]
The independent gate is n3_palindrome_covering_gate.py."""
import itertools, json, os, sys, time, hashlib
import sympy as sp
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
from sympy.polys.groebnertools import groebner as _groebner
from sympy.polys.orderings import grevlex

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '_n3_palindrome_covering_work')             # restartable work files, local only
CERT = os.path.join(HERE, 'results', 'n3_palindrome_covering.json')  # the certificate
J1, J2, J3, h0, h1, h2 = PARAMS = sp.symbols('J1 J2 J3 h0 h1 h2')
PATTERNS = ['XXX', 'XXY', 'XXZ', 'XYY', 'XYZ', 'XZZ', 'ZXX', 'ZXY', 'ZXZ', 'ZZZ']
CHARTS = {'J1=1': ['J1 - 1'], 'J1=0,J2=1': ['J1', 'J2 - 1'], 'J2=1': ['J2 - 1'], 'J2=0,J1=1': ['J2', 'J1 - 1']}
# two charts per pattern, either pair covering every connected point; XXY takes the J2 pair, on which its far
# elimination finishes in seconds (on the J1 pair the case splits run for hours)
CHART_PAIR = {pat: ('J2=1', 'J2=0,J1=1') if pat == 'XXY' else ('J1=1', 'J1=0,J2=1') for pat in PATTERNS}

# ---- the linear system ---------------------------------------------------------------------------------
def _mul1(a, b):
    if a == 0: return 0, b
    if b == 0: return 0, a
    if a == b: return 0, 0
    return (1 if (a, b) in ((1, 2), (2, 3), (3, 1)) else 3), 6 - a - b

def _mul(s, t):
    ph, out = 0, []
    for a, b in zip(s, t):
        q, c = _mul1(a, b); ph += q; out.append(c)
    return ph % 4, tuple(out)

def _anticommute(s, t):
    return sum(1 for a, b in zip(s, t) if a and b and a != b) % 2 == 1

def terms(pat):
    out = []
    for (a, b, J) in ((0, 1, J1), (0, 2, J2), (1, 2, J3)):
        for c in (1, 2, 3):
            s = [0, 0, 0]; s[a] = c; s[b] = c; out.append((J, tuple(s)))
    for site, hh in enumerate((h0, h1, h2)):
        s = [0, 0, 0]; s[site] = 'XYZ'.index(pat[site]) + 1; out.append((hh, tuple(s)))
    return out

def allowed(far):
    d = (1, 2) if far else (0, 3)
    return [s for s in itertools.product(range(4), repeat=3) if s[0] in d]

def system(pat, far):
    """[t, s] = 2ts = 2 i^k u for anticommuting strings, k odd; with real coefficients c_s (each end is closed
    under the adjoint, so Hermitian elements span it) dividing by 2i leaves real entries +-h_t."""
    cols = allowed(far)
    rows = {}
    for (h, t) in terms(pat):
        for j, s in enumerate(cols):
            if _anticommute(t, s):
                k, u = _mul(t, s)
                rows.setdefault(u, {})
                rows[u][j] = sp.expand(rows[u].get(j, 0) + (1 if k == 1 else -1) * h)
    return [{j: v for j, v in r.items() if v != 0} for r in rows.values()], len(cols)

# ---- exact elimination with case splits ----------------------------------------------------------------
RING, *_ = ring('J1 J2 J3 h0 h1 h2', QQ, grevlex)
_GB, _FAC = {}, {}

def _to_ring(e):
    return RING.from_expr(sp.sympify(e)) if not hasattr(e, 'ring') else e

def _gb(eqs):
    key = tuple(sorted(str(e) for e in eqs))
    if key not in _GB:
        _GB[key] = _groebner(list(eqs), RING) if eqs else []
    return _GB[key]

def _factors(v):
    if v not in _FAC:
        _FAC[v] = v.factor_list()
    return _FAC[v]

def _nf(p, G):
    return p.rem(G) if G else p

class Leaf:
    def __init__(self, eqs, neqs, nullity):
        self.eqs = [e.as_expr() for e in eqs]
        self.neqs = [n.as_expr() for n in neqs]
        self.nullity = nullity

def _eliminate(rows, ncols, eqs, neqs, leaves):
    G = _gb(eqs)
    if any(g.is_ground and g != 0 for g in G):
        return                                   # 1 in the ideal: empty branch
    if any(_nf(n, G) == 0 for n in neqs):
        return                                   # a factor recorded nonzero vanishes on the branch
    rows = [{j: _nf(v, G) for j, v in r.items()} for r in rows]
    rows = [r for r in ({j: v for j, v in r.items() if v != 0} for r in rows) if r]
    neqkeys = {n.monic() for n in neqs}
    while rows:
        piv = next(((i, j) for i, r in enumerate(rows) for j, v in r.items() if v.is_ground), None)
        if piv is None:
            piv = next(((i, j) for i, r in enumerate(rows) for j, v in r.items()
                        if all(f.monic() in neqkeys for f, _ in _factors(v)[1])), None)
        if piv is None:                          # split: some factor of the pivot vanishes, or none does
            cands = [(max(sum(mm) for mm in v.monoms()), len(v), i, j, v)
                     for i, r in enumerate(rows) for j, v in r.items()]
            v = min(cands, key=lambda c: (c[0], c[1]))[4]
            facs = [f for f, _ in _factors(v)[1] if not f.is_ground]
            for f in facs:
                _eliminate([dict(r) for r in rows], ncols, eqs + [f], neqs, leaves)
            _eliminate([dict(r) for r in rows], ncols, eqs, neqs + facs, leaves)
            return
        i, j = piv
        pr = rows.pop(i); p = pr[j]; new = []
        for r in rows:
            if j not in r:
                new.append(r); continue
            e = r[j]; nr = {}
            for k in set(r) | set(pr):
                val = _nf(p * r.get(k, RING.zero) - e * pr.get(k, RING.zero), G)
                if val != 0: nr[k] = val
            if nr:
                # divide out only factors recorded nonzero (a factor that may vanish stays: dividing by it
                # would invent a rank where it vanishes)
                vals = list(nr.values()); g = vals[0]
                for v2 in vals[1:]:
                    g = g.gcd(v2)
                    if g.is_ground: break
                keep = RING.one
                if not g.is_ground:
                    for fac, mult in _factors(g)[1]:
                        if fac.monic() in neqkeys: keep *= fac ** mult
                if keep != 1:
                    nr = {k: v2.exquo(keep) for k, v2 in nr.items()}
                new.append(nr)
        rows = new
        ncols -= 1
    leaves.append(Leaf(eqs, neqs, ncols))

def run(pat, far, eqs=(), neqs=()):
    rows, n = system(pat, far)
    rows = [{j: _to_ring(v) for j, v in r.items()} for r in rows]
    leaves = []
    _eliminate(rows, n, [_to_ring(e) for e in eqs], [_to_ring(x) for x in neqs], leaves)
    return leaves

# ---- the three classes and the region tests ------------------------------------------------------------
t_ = sp.Symbol('t')
AX = {'X': (1, 0, 0), 'Y': (0, 1, 0), 'Z': (0, 0, 1)}
ROT = {'x': lambda v: (v[0], -v[1], -v[2]), 'y': lambda v: (-v[0], v[1], -v[2]),
       'x+y': lambda v: (v[1], v[0], -v[2]), 'x-y': lambda v: (-v[1], -v[0], -v[2])}

def field(pat, l):
    return tuple((h0, h1, h2)[l] * c for c in AX[pat[l]])

def components(pat):
    """Each class as a list of generators that must all vanish. For B the four axes x, y, x+y, x-y are the
    only ones a pi rotation about which maps single-letter fields to single-letter fields."""
    comps = [('A' + p, [(h0, h1, h2)[l] for l in range(3) if pat[l] != p]) for p in 'XY']
    for n, R in ROT.items():
        fd, fu, fv = field(pat, 0), field(pat, 1), field(pat, 2)
        gens = [J1 - J2] + [a - b for a, b in zip(R(fd), fd)] + [a - b for a, b in zip(R(fu), fv)]
        comps.append(('B' + n, [g for g in map(sp.expand, gens) if g != 0]))
    gens = [J1*J2 + J1*J3 + J2*J3] + list(field(pat, 0)) + [field(pat, 1)[2], field(pat, 2)[2]]
    gens += [sp.expand(sum(c * c for c in field(pat, 1)) - sum(c * c for c in field(pat, 2)))]
    comps.append(('C', [g for g in gens if g != 0]))
    return sorted(comps, key=lambda c: len(c[1]))

_GBS = {}
def vanishes(g, E, N):
    """g vanishes on {E = 0, every N != 0} over C: g in <E>, or 1 in <E, 1 - t*prod(N)*g> (Rabinowitsch)."""
    key = tuple(map(str, E))
    if key not in _GBS:
        _GBS[key] = sp.groebner(list(E), *PARAMS, order='grevlex', domain='QQ') if E else None
    if _GBS[key] is not None and _GBS[key].reduce(sp.expand(g))[1] == 0:
        return True
    if not E:
        return sp.expand(g) == 0
    prodN = sp.Mul(*N) if N else sp.Integer(1)
    return list(sp.groebner(list(E) + [1 - t_ * prodN * g], t_, *PARAMS, order='grevlex', domain='QQ').exprs) == [1]

def inside(pat, E, N):
    return next((name for name, comp in components(pat) if all(vanishes(g, E, N) for g in comp)), None)

def disconnected(E, N):
    return sum(1 for j in (J1, J2, J3) if vanishes(j, E, N)) >= 2

def no_real_point(E, N):
    """An equation whose terms all have one sign and even exponents, one term built only of variables recorded
    nonzero, has no real zero there."""
    nonzero = {sp.Symbol(str(n)) for n in N if sp.sympify(n).is_Symbol}
    for e in E:
        tm = sp.Poly(e, *PARAMS).terms()
        if all(all(x % 2 == 0 for x in mon) for mon, _ in tm) and (all(c > 0 for _, c in tm) or all(c < 0 for _, c in tm)):
            if any(all(PARAMS[i] in nonzero for i, x in enumerate(mon) if x) for mon, _ in tm):
                return True
    return False

# Regions the sign test above cannot see as real-empty, each with a sum of squares q = sum w_i p_i^2 (w_i > 0)
# that vanishes on the region and whose strict term p_k has every factor recorded nonzero or incompatible with
# the region; then q > 0 at any real point of it, so it has none. Keyed by the region's digest, checked on use.
_R = sp.Rational
SOS_CERTIFICATES = {
    '7c30fb386dd80034': ([(_R(1, 8), J1**3*h1), (_R(1, 4), J1*h1**2), (_R(1, 8), J1*h1*h2**2), (_R(1, 4), J1*h2**2),
                          (_R(1, 8), h1**2*h2), (_R(1, 8), h2)], 5),                       # XXY, chart J2 = 1
    '40fd2c968eafb49e': ([(1, h1*h2*(J1 - 1)**3), (16, J1*(J1 + 1)*(J1**2 + 1))], 0),        # XXY, chart J2 = 1
    'fddf347ffef5baa5': ([(_R(1, 5), h1**3*h2), (_R(3, 5), h1**2), (_R(1, 5), h1*h2**3)], 1),  # ZXY, chart J1 = 1
    '3b6972a79868ef6c': ([(_R(1, 4), J2**2*h1*h2), (_R(1, 8), J2**2*h1), (_R(1, 8), J2*h2**3), (_R(1, 8), J2*h2),
                          (_R(1, 8), h1**3), (_R(1, 4), h1*h2)], 4),                       # ZXY, chart J1 = 1
}

def sos_certified(E, N):
    cert = SOS_CERTIFICATES.get(_digest(E, N))
    if cert is None:
        return False
    terms, k = cert
    if not all(w > 0 for w, _ in terms) or not vanishes(sp.expand(sum(w * p**2 for w, p in terms)), E, N):
        return False
    for fac, _ in sp.factor_list(terms[k][1])[1]:
        recorded = any(sp.expand(fac - n) == 0 or sp.expand(fac + n) == 0 for n in N)
        if not recorded and list(sp.groebner(list(E) + [fac], *PARAMS, order='grevlex', domain='QQ').exprs) != [1]:
            return False
    return True

def decide_region(pat, E, N):
    """Why a palindromic region is no counterexample: 'disconnected', 'no real point' or 'inside <class>'."""
    cls = inside(pat, E, N)
    if cls: return 'inside ' + cls
    if disconnected(E, N): return 'disconnected'
    if no_real_point(E, N): return 'no real point'
    if sos_certified(E, N): return 'no real point (sum-of-squares certificate)'
    return None

# ---- driver ------------------------------------------------------------------------------------------
def _digest(E, N):
    return hashlib.sha256(json.dumps([sorted(map(str, E)), sorted(map(str, N))]).encode()).hexdigest()[:16]

def far_task(job):
    pat, chart = job
    path = os.path.join(OUT, f'far_{pat}_{chart}.json')
    if not os.path.exists(path):
        t = time.time()
        leaves = run(pat, True, CHARTS[chart])
        json.dump([{'eqs': list(map(str, l.eqs)), 'neqs': list(map(str, l.neqs)), 'far': l.nullity} for l in leaves],
                  open(path, 'w'), indent=0)
        print(f'far {pat} {chart}: {len(leaves)} leaves, {time.time() - t:.0f} s', flush=True)
    return pat, chart, json.load(open(path))

def zzz_near_floor():
    """For the pattern ZZZ the near end holds I, M, M^2, M^3 (M = Z0 + Z1 + Z2, eigenvalues +-1, +-3) at every
    parameter point: H conserves M and every power is diagonal, so commutes with the jump. Checked here on the
    near system's own rows, and their independence by rank. Returns the floor, 4."""
    cols = allowed(False); rows, n = system('ZZZ', False)
    Ms = {(3, 0, 0): 1, (0, 3, 0): 1, (0, 0, 3): 1}
    powers, P = [], {(0, 0, 0): 1}
    for _ in range(4):
        powers.append(P)
        nxt = {}
        for s1, c1 in P.items():
            for s2, c2 in Ms.items():
                ph, u = _mul(s1, s2)
                nxt[u] = nxt.get(u, 0) + c1 * c2 * (1, 1j, -1, -1j)[ph]
        P = {k: v for k, v in nxt.items() if v != 0}
    vecs = []
    for P in powers:
        assert all(abs(v.imag) == 0 if isinstance(v, complex) else True for v in P.values())
        vec = [int(P.get(c, 0).real if isinstance(P.get(c, 0), complex) else P.get(c, 0)) for c in cols]
        assert all(sp.expand(sum(r.get(j, 0) * vec[j] for j in range(n))) == 0 for r in rows)
        vecs.append(vec)
    assert sp.Matrix(vecs).rank() == 4
    return 4

def leaf_task(job):
    pat, chart, idx, lf = job
    path = os.path.join(OUT, f'leaf_{pat}_{chart}_{idx:03d}.json')
    if os.path.exists(path):
        return json.load(open(path))
    t = time.time()
    E = [sp.sympify(e) for e in lf['eqs']]; N = [sp.sympify(n) for n in lf['neqs']]
    k = lf['far']
    res = {'pattern': pat, 'chart': chart, 'leaf': idx, 'far': k, 'digest': _digest(E, N)}
    if pat == 'ZZZ' and k < zzz_near_floor():
        why = f'not palindromic: near >= 4 > far = {k} (I, M, M^2, M^3)'
    else:
        why = decide_region(pat, E, N)
    if why:
        res['decided_on_far_region'] = why
        res['palindromic_regions'] = []
    else:
        regs = []
        for nl in run(pat, False, E, N):
            if nl.nullity != k:
                continue                          # near != far: not palindromic
            regs.append({'eqs': list(map(str, nl.eqs)), 'neqs': list(map(str, nl.neqs)),
                         'verdict': decide_region(pat, nl.eqs, nl.neqs) or 'OPEN'})
        res['palindromic_regions'] = regs
    res['open'] = sum(1 for r in res['palindromic_regions'] if r['verdict'] == 'OPEN')
    res['seconds'] = round(time.time() - t, 1)
    json.dump(res, open(path, 'w'), indent=0)
    print(f"leaf {pat} {chart} {idx}: far {k}, {res.get('decided_on_far_region') or str(len(res['palindromic_regions'])) + ' palindromic regions'}"
          f"{', OPEN ' + str(res['open']) if res['open'] else ''}, {res['seconds']} s", flush=True)
    return res

def main(workers):
    from multiprocessing import Pool
    os.makedirs(OUT, exist_ok=True)
    jobs = [(p, c) for p in PATTERNS for c in CHART_PAIR[p]]
    with Pool(workers) as pool:
        fars = pool.map(far_task, jobs)
        tasks = [(p, c, i, lf) for p, c, far in fars for i, lf in enumerate(far) if lf['far'] > 0]
        tasks.sort(key=lambda t: sum(len(x) for x in t[3]['eqs']))
        print(f'{len(tasks)} far leaves with far >= 1', flush=True)
        results = list(pool.imap_unordered(leaf_task, tasks))
    summary = {'far_leaves': sum(len(f) for _, _, f in fars), 'far_leaves_with_far_ge_1': len(results),
               'decided_on_far_region': {}, 'palindromic_regions': {}, 'open': sum(r['open'] for r in results)}
    for r in results:
        if 'decided_on_far_region' in r:
            w = r['decided_on_far_region']; summary['decided_on_far_region'][w] = summary['decided_on_far_region'].get(w, 0) + 1
        for g in r['palindromic_regions']:
            summary['palindromic_regions'][g['verdict']] = summary['palindromic_regions'].get(g['verdict'], 0) + 1
    by_leaf = {(r['pattern'], r['chart'], r['leaf']): r for r in results}
    leaves = []
    for p, c, far in fars:
        for i, lf in enumerate(far):
            entry = {'pattern': p, 'chart': c, 'leaf': i, 'eqs': lf['eqs'], 'neqs': lf['neqs'], 'far': lf['far']}
            if lf['far'] > 0:
                r = by_leaf[(p, c, i)]
                entry.update({k: r[k] for k in ('digest', 'palindromic_regions', 'open')})
                if 'decided_on_far_region' in r: entry['decided_on_far_region'] = r['decided_on_far_region']
            leaves.append(entry)
    json.dump({'summary': summary, 'charts': {p: list(CHART_PAIR[p]) for p in PATTERNS}, 'far_leaves': leaves},
              open(CERT, 'w', newline='\n'), indent=0)
    print(json.dumps(summary, indent=1))
    print('OPEN regions:', summary['open'])

if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 20)
