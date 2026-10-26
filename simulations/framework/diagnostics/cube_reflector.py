"""The Theodolite: does a product of dressed cube moves reflect the generator?

See experiments/THE_CUBE_SITE_BY_SITE.md. Per site one of the cube's moves P,
dressed with a cell weight w (a diagonal map on the Z-basis cells) after or before it, E_l = w·P or P·w: 40
distinct families of the 192 moves. The product E = ⊗E_l is
asked to carry L onto
  -L - c      a Pi-reflection ('pi'), the point reflection about -c/2, or
  -L^† - c    a reflection ('refl'), the line reflection across Re = -c/2,
with c = -Σ_l tr L_l / 2. On a given graph the question splits exactly into one linear equation per site,
E_l·L_l° = 𝒯(L_l°)·E_l (L_l° the traceless part, 𝒯 the target map), and one per support set,
(⊗_{l∈S} E_l)·L_S = 𝒯(L_S)·(⊗ E_l). Conventions are those of framework.cube_tables.

The letters answer first (framework.diagnostics.cube_letters: monomial sandwiches, every equation binomial, solved
in closed form without a Groebner basis). The pre-check, the cube search below and the general machine (any
invertible site maps, sympy Groebner bases over Q(i)) only shake hands: did the letters stay silent because
nothing site-local exists, or because the element is not a monomial sandwich?

Public API:
  site_family(L_site, relation)                 basis of {E: E·L° = 𝒯(L°)·E}
  precheck(setting, relation)                   ('NONE', why) | ('PASS', None) | ('SKIPPED', why)
  cube_search(setting, relation, max_calls)     ('FOUND', hit) | ('NOT FOUND', None) | ('BUDGET', None)
  machine(setting, relation, timeout)           'NONE' | 'FOUND' | 'OPEN'
  answer(setting, relation, ...)                the driver -> Verdict
  Cell(jumps, fields, windows), every_n(cell, relation, rate_free)
                                                'UNIFORM' | 'TURNED' | ('PERIODIC', p) | None
"""
from __future__ import annotations

import multiprocessing as mp
from dataclasses import dataclass, field
from functools import lru_cache
from itertools import product

import sympy as sp

from .. import cube_tables as ct

RELATIONS = ('pi', 'refl')


def target(M, relation):
    """𝒯 of a relation. 'copy' (E·L = L·E, a symmetry, c = 0) is read by the letters alone (cube_letters): the
    identity is always a copy, so the pre-check, the cube search and the machine below have nothing to decide."""
    if relation == 'pi':
        return -M
    if relation == 'refl':
        return -M.H
    if relation == 'copy':
        return M
    raise ValueError(f"relation must be 'pi', 'refl' or 'copy', not {relation!r}")


def _reflections_only(relation):
    if relation not in RELATIONS:
        raise ValueError(f"the handshake asks 'pi' or 'refl', not {relation!r}; a copy is read by "
                         "cube_letters.read / cube_letters.every_n")


def traceless(M):
    return M - (M.trace() / M.shape[0]) * sp.eye(M.shape[0])


def constant_c(setting):
    return sp.nsimplify(-sum(ct.site_table(setting, l).trace() for l in range(setting.n)) / 2)


def site_family(L_site, relation):
    """Basis of the site maps E with E·L° = 𝒯(L°)·E (a nullspace over Q(i)); E as 4x4 cell maps."""
    L0 = traceless(L_site)
    return _nullspace_of([(L0, target(L0, relation))])


def _nullspace_of(pairs):
    """Basis of {E: E·A = B·E for every (A, B)}."""
    I4 = sp.eye(4)
    M = sp.Matrix.vstack(*[sp.kronecker_product(I4, A.T) - sp.kronecker_product(B, I4) for A, B in pairs])
    return [sp.Matrix(4, 4, list(v)) for v in M.nullspace()]


def _as_operator(v):
    """A cell vector (row-major) as the 2x2 operator it stacks."""
    return sp.Matrix([[v[0], v[1]], [v[2], v[3]]])


# ---------------------------------------------------------------- pre-check

def precheck(setting, relation, max_charts=4096):
    """Exact necessary condition. Applying the relation to the identity fixes one operator per site:
        unital L (L(I) = 0):  W = E(I)     ('pi': L W = -c W;    'refl': L^† W = -c W)
        every Lindbladian:    W = E^†(I)   ('pi': L^† W = -c W)  or  W = E^-1(I)  ('refl': L W = -c W)
    The site equation puts each factor W_l in its site's eigenspace at tr L_l / 2, and every support set needs
    [H_S, ⊗_{l∈S} W_l] = 0. Decided by a chart cover (per site the first nonzero coordinate set to 1)."""
    n = setting.n
    sites = [ct.site_table(setting, l) for l in range(n)]
    vec_i = sp.Matrix([1, 0, 0, 1])
    unital = all((L * vec_i).applyfunc(sp.expand).is_zero_matrix for L in sites)
    supports = sorted({tuple(sorted(s)) for _, s, _ in ct.normalized_terms(setting.terms) if len(s) >= 2})
    hams = {S: ct.hamiltonian_on(setting, S) for S in supports}
    rows = (['column'] if unital else []) + ['row']
    for row in rows:
        use_l = (relation, row) in (('pi', 'column'), ('refl', 'row'))
        spaces = []
        for l, L in enumerate(sites):
            op = L if use_l else L.H
            ns = (op - (L.trace() / 2) * sp.eye(4)).applyfunc(sp.expand).nullspace()
            if not ns:
                return 'NONE', f"site {l} has no operator in its eigenspace at tr L/2 ({row} of the identity)"
            spaces.append(ns)
        n_charts = 1
        for ns in spaces:
            n_charts *= len(ns)
        if n_charts > max_charts:
            return 'SKIPPED', f"{n_charts} charts exceed the cap {max_charts}"
        if not supports:
            continue
        found_any = False
        for chart in product(*[range(len(ns)) for ns in spaces]):
            W, gens = [], []
            for l, (ns, j) in enumerate(zip(spaces, chart)):
                a = sp.symbols(f'a{l}_0:{len(ns)}')
                coeffs = [0] * j + [1] + list(a[j + 1:])
                gens += list(a[j + 1:])
                W.append(_as_operator(sum((cf * v for cf, v in zip(coeffs, ns)), sp.zeros(4, 1))))
            eqs = []
            for S in supports:
                WS = sp.kronecker_product(*[W[l] for l in S])
                eqs += [e for e in (hams[S] * WS - WS * hams[S]).applyfunc(sp.expand) if e != 0]
            if not eqs:
                found_any = True
                break
            if not gens:
                continue
            if sp.groebner(eqs, *gens, order='grevlex', domain=sp.QQ_I).exprs != [1]:
                found_any = True
                break
        if not found_any:
            return 'NONE', f"no product of site operators commutes with every support set ({row} of the identity)"
    return 'PASS', None


# ---------------------------------------------------------------- cube search

def _dressed(P, D, side):
    return D * P if side == 'after' else P * D


@lru_cache(maxsize=None)
def _families():
    """The distinct families {w·P} (weight after the move, 24 right cosets) and {P·w} (weight before it, 24 left
    cosets): 40, since 8 coincide. [(name, P, side)]."""
    units = [sp.diag(*[1 if i == j else 0 for i in range(4)]) for j in range(4)]
    out, spans = [], []
    for side, cos in (('after', ct.move_cosets('right')), ('before', ct.move_cosets('left'))):
        for name, P in cos:
            span = sp.Matrix.hstack(*[sp.Matrix(list(_dressed(P, U, side))) for U in units])
            if any(span.rank() == G.rank() == sp.Matrix.hstack(span, G).rank() for G in spans):
                continue
            spans.append(span)
            out.append((name, P, side))
    return tuple(out)


def cube_candidates(L_site, relation):
    """Per family of dressed cube moves, the cell weights w with w·P (or P·w) in the site family; kept when an
    invertible w is possible. [(name, P, weight basis, side)]."""
    L0 = traceless(L_site)
    TL = target(L0, relation)
    d = sp.symbols('d0:4')
    out = []
    for name, P, side in _families():
        E = _dressed(P, sp.diag(*d), side)
        eqs = (E * L0 - TL * E).applyfunc(sp.expand)
        A = sp.Matrix([[sp.diff(e, x) for x in d] for e in eqs])
        ns = A.nullspace()
        if not ns or any(all(v[j] == 0 for v in ns) for j in range(4)):
            continue
        out.append((name, P, [sp.diag(*list(v)) for v in ns], side))
    return out


class _Budget(Exception):
    pass


def _system(assigned, cands, choice, tables, relation, tie=None):
    """The cube-shape system for the assigned sites (tie: site -> variable group): the support-set equations and
    t·det − 1 per group (core), and the normalisation w_00 = 1 per group (norm). An invertible cell weight has every
    diagonal entry nonzero, so the normalisation loses no solution."""
    tie = tie or {l: l for l in assigned}
    groups = sorted({tie[l] for l in assigned})
    E, gens, core, cut = {}, [], [], []
    for g in groups:
        name, P, basis, side = cands[g][choice[g]]
        y = sp.symbols(f'y{g}_0:{len(basis)}')
        D = sum((yy * B for yy, B in zip(y, basis)), sp.zeros(4, 4))
        E[g] = _dressed(P, D, side)
        t = sp.Symbol(f't{g}')
        gens += list(y) + [t]
        core.append(t * sp.prod([D[i, i] for i in range(4)]) - 1)
        cut.append(D[0, 0] - 1)
    for S, LS in tables:
        ES = sp.kronecker_product(*[E[tie[l]] for l in S])
        core += [e for e in (ES * LS - target(LS, relation) * ES).applyfunc(sp.expand) if e != 0]
    return core, cut, gens, E


def _solvable(eqs, gens):
    return sp.groebner(eqs, *gens, order='grevlex', domain=sp.QQ_I).exprs != [1]


def cube_search(setting, relation, max_calls=2000, tie=None, groups=None):
    """Cube moves per site (or per tied group), dressed with cell weights. Backtracking over the sites; at each
    completed support set an exact Groebner check of every completed support set, each weight normalised to
    w_00 = 1. A hit proves a cube-move answer exists; a miss after every combination proves that none exists with
    a weight along Z, before or after the move ('NOT FOUND'; weights along X or Y are not searched); 'BUDGET' when
    max_calls ran out."""
    n = setting.n
    tie = tie or {l: l for l in range(n)}
    groups = groups if groups is not None else sorted(set(tie.values()))
    rep = {g: min(l for l in range(n) if tie[l] == g) for g in groups}
    cands = {g: cube_candidates(ct.site_table(setting, rep[g]), relation) for g in groups}
    if any(not cands[g] for g in groups):
        return 'NOT FOUND', None
    tables = ct.support_tables(setting)
    order = list(groups)
    pos = {g: i for i, g in enumerate(order)}
    done_at = {g: [] for g in order}
    for S, LS in tables.items():
        done_at[order[max(pos[tie[l]] for l in S)]].append((S, LS))
    calls = 0
    choice = {}

    def rec(i):
        nonlocal calls
        if i == len(order):
            return True
        g = order[i]
        for ci in range(len(cands[g])):
            choice[g] = ci
            if done_at[g]:
                calls += 1
                if calls > max_calls:
                    raise _Budget
                completed = [st for h in order[:i + 1] for st in done_at[h]]
                assigned = [l for l in range(n) if pos[tie[l]] <= i]
                core, cut, gens, _ = _system(assigned, cands, choice, completed, relation, tie)
                if not _solvable(core + cut, gens):
                    continue
            if rec(i + 1):
                return True
        del choice[g]
        return False

    try:
        if not rec(0):
            return 'NOT FOUND', None
    except _Budget:
        return 'BUDGET', None
    core, norm, gens, E = _system(list(range(n)), cands, choice, list(tables.items()), relation, tie)
    names = {g: (cands[g][choice[g]][0], cands[g][choice[g]][3]) for g in groups}
    return 'FOUND', {'core': core + norm, 'gens': gens, 'E': E, 'names': names, 'tie': tie}


_VALUES = (1, 2, -1, 3, sp.Rational(1, 2), -2, 5, sp.Rational(1, 3))


def _consistent(core, gens, subs):
    eqs = [sp.expand(e.subs(subs)) for e in core]
    eqs = [e for e in eqs if e != 0]
    if any(e.is_number for e in eqs):
        return False
    free = [g for g in gens if g not in subs]
    return not eqs or sp.groebner(eqs, *free, order='grevlex', domain=sp.QQ_I).exprs != [1]


def extract(hit):
    """A point of the found family (support-set equations, t·det − 1, w_00 = 1 per site). The weight
    coordinates are fixed one by one to small rationals, each choice kept when the system stays consistent; the
    ones no small rational fits are determined algebraically and solved last. The point is verified by exact
    substitution. Returns {group: 4x4 matrix} or None."""
    core, gens = list(hit['core']), list(hit['gens'])
    subs = {}
    for y in (v for v in gens if not v.name.startswith('t')):
        for val in _VALUES:
            subs[y] = val
            if _consistent(core, gens, subs):
                break
            del subs[y]
    rest = [v for v in gens if v not in subs]
    eqs = [e for e in (sp.expand(e.subs(subs)) for e in core) if e != 0]
    sols = sp.solve(eqs, rest, dict=True) if eqs else [{}]
    for sol in sols:
        full = {**subs, **sol}
        point = {g: M.subs(full).applyfunc(sp.nsimplify) for g, M in hit['E'].items()}
        if any(M.free_symbols for M in point.values()):
            continue
        if all(sp.simplify(M.det()) != 0 for M in point.values()) and \
                all(sp.simplify(e.subs(full)) == 0 for e in core):
            return point
    return None


def reading(hit, point=None):
    """Per site the cube move and, when a point is known, its cell weight in the cube's words."""
    out = []
    for l in sorted(hit['tie']):
        g = hit['tie'][l]
        name, side = hit['names'][g]
        entry = {'site': l, 'move': name, 'weight_side': side}
        if point is not None:
            E = point[g]
            P = next(M for nm, M, sd in _families() if nm == name and sd == side)
            D = (E * P.inv() if side == 'after' else P.inv() * E).applyfunc(sp.nsimplify)
            w = [D[i, i] for i in range(4)]
            entry['weight'] = w
            if w[0] != 0 and sp.simplify(w[0] * w[3] - w[1] * w[2]) == 0:
                u, v = sp.nsimplify(w[2] / w[0]), sp.nsimplify(w[1] / w[0])
                entry['turn'] = sp.simplify(sp.arg(u / v) / 2)
                entry['thermal'] = sp.simplify(sp.sqrt(sp.Abs(u * v)))
                entry['imaginary_turn'] = sp.simplify(sp.sqrt(sp.Abs(u / v)))
                entry['eta_phase'] = sp.simplify(sp.arg(u * v) / 2)
            else:
                entry['cell_only_ratio'] = sp.simplify(w[0] * w[3] / (w[1] * w[2])) if w[1] * w[2] != 0 else None
        out.append(entry)
    return out


# ---------------------------------------------------------------- the machine (handshake)

def _machine_core(setting, relation):
    n = setting.n
    fams = [site_family(ct.site_table(setting, l), relation) for l in range(n)]
    if any(not f for f in fams):
        return 'NONE'
    E, gens, eqs = [], [], []
    for l, f in enumerate(fams):
        x = sp.symbols(f'x{l}_0:{len(f)}')
        El = sum((xx * B for xx, B in zip(x, f)), sp.zeros(4, 4))
        t = sp.Symbol(f't{l}')
        E.append(El)
        gens += list(x) + [t]
        eqs.append(t * El.det() - 1)
    for S, LS in ct.support_tables(setting).items():
        ES = sp.kronecker_product(*[E[l] for l in S])
        eqs += [e for e in (ES * LS - target(LS, relation) * ES).applyfunc(sp.expand) if e != 0]
    G = sp.groebner(eqs, *gens, order='grevlex', domain=sp.QQ_I)
    return 'NONE' if G.exprs == [1] else 'FOUND'


def _machine_worker(setting, relation, conn):
    conn.send(_machine_core(setting, relation))
    conn.close()


def machine(setting, relation, timeout=600):
    """General invertible site maps, the unsliced joint ideal: 'NONE' (basis [1], exact), 'FOUND', or 'OPEN'
    when the timeout (seconds, a separate process) runs out. timeout=None runs in this process."""
    if timeout is None:
        return _machine_core(setting, relation)
    ctx = mp.get_context('spawn')
    parent, child = ctx.Pipe(duplex=False)
    proc = ctx.Process(target=_machine_worker, args=(setting, relation, child))
    proc.start()
    if parent.poll(timeout):
        result = parent.recv()
        proc.join()
        return result
    proc.terminate()
    proc.join()
    return 'OPEN'


# ---------------------------------------------------------------- the driver

@dataclass
class Verdict:
    """kind is about SITE-LOCAL elements only: 'NONE' means no product of site maps carries the relation, not that
    the spectrum fails to pair. pairing is the spectrum's answer where one is known: 'PAIRS' whenever a site-local
    element is found, and after a silence F158's exact two-end count for Pauli-letter jumps ('PAIRS' | 'BREAKS',
    with the carrier in element); for other jumps (sigma^+-, thermal) 'PAIRS' when the letters read an element for
    the other relation, since either relation pairs the spectrum, and None otherwise. element is F158's U, and it
    carries the line reflection by one-sided multiplication, rho -> rho·U maps L to -L^† - 2 sigma ('refl'); it is
    attached only to a 'refl' verdict (the same spectrum pairs under 'pi', through no U of this form)."""
    kind: str                    # 'CUBE MOVES' | 'NO CUBE MOVE' | 'NONE' | 'OPEN' | 'OUTSIDE'
    relation: str
    c: object = None
    moves: list = None
    point: dict = None
    note: str = ''
    pairing: str = None          # 'PAIRS' | 'BREAKS' | None
    element: dict = None         # F158's U {Pauli string: coefficient}: rho -> rho·U gives -L^† - 2 sigma


def _hand_over(verdict, setting):
    """After a silence: does the spectrum pair at all? F158's exact two-end count, for Pauli-letter jumps."""
    from .f158_end_count import end_count_of_setting
    if verdict.kind == 'NONE':
        verdict.note += '; so no site-local element (NONE is about site maps, not the spectrum)'
    ec = end_count_of_setting(setting)
    if ec is None:
        # either relation pairs the spectrum about -c/2 (L's spectrum is closed under conjugation), so the other
        # relation's letters decide it when they read an element
        from . import cube_letters as cl
        other = 'refl' if verdict.relation == 'pi' else 'pi'
        rd = cl.read(setting, other)
        if rd is not None and cl.verify(setting, other, rd.maps):
            verdict.pairing = 'PAIRS'
            verdict.note += (f"; the spectrum PAIRS: the letters read a site-local element for '{other}' "
                             f"(a jump is not a Pauli letter, F158 does not apply)")
            return verdict
        verdict.note += ("; pairing not decided here (a jump is not a Pauli letter, F158 does not apply, and the "
                         f"letters read no element for '{other}')")
        return verdict
    verdict.pairing = 'PAIRS' if ec.pairs else 'BREAKS'
    verdict.note += f'; F158 (exact): near {ec.near}, far {ec.far}, the spectrum '
    if not ec.pairs:
        verdict.note += 'does not pair'
        return verdict
    verdict.note += 'PAIRS'
    if ec.element is None:
        verdict.note += ' (no basis vector of the far kernel squares to a scalar; none named)'
    elif verdict.relation == 'refl':
        verdict.element = ec.element
        verdict.note += ', carried by rho -> rho·U (element), -L^† - 2 sigma'
    else:
        verdict.note += "; F158's U carries the 'refl' relation, ask that one for it"
    return verdict


def answer(setting, relation, max_calls=2000, timeout=600):
    """The letters answer first (framework.diagnostics.cube_letters: a monomial sandwich, every equation binomial,
    solved in closed form); when they are silent the handshake decides: the pre-check (exact), then the cube search
    and the machine (Groebner bases). Where no site-local element is found, the hand-over asks F158's exact
    two-end count whether the spectrum pairs anyway (framework.diagnostics.f158_end_count)."""
    from . import cube_letters as cl
    _reflections_only(relation)
    reason = ct.outside_reason(setting)
    if reason:
        return Verdict('OUTSIDE', relation, note=reason)
    c = constant_c(setting)
    rd = cl.read(setting, relation)
    if rd is not None and cl.verify(setting, relation, rd.maps):
        return Verdict('CUBE MOVES', relation, c, moves=cl.describe(rd), point=rd.maps,
                       note='read from the letters: a monomial sandwich, every equation binomial', pairing='PAIRS')
    pre, why = precheck(setting, relation)
    if pre == 'NONE':
        return _hand_over(Verdict('NONE', relation, c, note=f'pre-check (exact): {why}'),
                          setting)
    status, hit = cube_search(setting, relation, max_calls)
    if status == 'FOUND':
        point = extract(hit)
        return Verdict('CUBE MOVES', relation, c, moves=reading(hit, point), point=point,
                       note='the letters are silent; found by the cube search' + ('' if point else
                                                                                   ', weights not extracted'),
                       pairing='PAIRS')
    result = machine(setting, relation, timeout)
    if result == 'FOUND':
        return Verdict('NO CUBE MOVE', relation, c,
                       note=f'cube search {status} (weights along Z); a general site map exists', pairing='PAIRS')
    if result == 'NONE':
        return _hand_over(Verdict('NONE', relation, c,
                                  note=f'machine (exact); cube search {status}'), setting)
    return _hand_over(Verdict('OPEN', relation, c, note=f'cube search {status}; machine over {timeout}s'), setting)


# ---------------------------------------------------------------- every N

@dataclass
class Cell:
    """One site of a chain repeated: jumps [(A, rate)], fields [(letter, coefficient)], windows
    [(letters, offsets, coefficient)] with offsets from 0, e.g. ('XZX', (0, 1, 2), 3)."""
    jumps: list = field(default_factory=list)
    fields: list = field(default_factory=list)
    windows: list = field(default_factory=list)

    def normalized(self):
        """The same cell with identity letters removed from its windows, each window's offsets shifted to start at 0
        (the chain is translation invariant), and one-letter windows moved to the fields."""
        fields, windows = list(self.fields), []
        for letters, off, coef in ct.normalized_terms(self.windows):
            if len(off) == 1:
                fields.append((letters, coef))
            else:
                m = min(off)
                windows.append((letters, tuple(o - m for o in off), coef))
        return Cell(list(self.jumps), fields, windows)

    def span(self):
        return max(max(off) for _, off, _ in self.windows) + 1 if self.windows else 1

    def chain(self, length):
        terms = [(letter, (l,), coef) for l in range(length) for letter, coef in self.fields]
        for letters, off, coef in self.windows:
            for s in range(length - max(off)):
                terms.append((letters, tuple(s + o for o in off), coef))
        jumps = [((l,), A, rate) for l in range(length) for A, rate in self.jumps]
        return ct.Setting(length, terms, jumps)


def _cell_family(cell, relation, rate_free):
    one = cell.chain(1)
    if not rate_free:
        return site_family(ct.site_table(one, 0), relation)
    pairs = []
    for A, _ in cell.jumps:
        D0 = traceless(ct.dissipator(sp.Matrix(A), 1))
        pairs.append((D0, target(D0, relation)))
    if cell.fields:
        F0 = traceless(ct.site_table(ct.Setting(1, [(lt, (0,), c) for lt, c in cell.fields], []), 0))
        pairs.append((F0, target(F0, relation)))
    return _nullspace_of(pairs) if pairs else [sp.Matrix(4, 4, lambda i, j: 1 if 4 * i + j == k else 0)
                                              for k in range(16)]


def _uniform_or_turned(cell, relation, rate_free, turned):
    fam = _cell_family(cell, relation, rate_free)
    if not fam:
        return False
    x = sp.symbols(f'x0:{len(fam)}')
    E = sum((xx * B for xx, B in zip(x, fam)), sp.zeros(4, 4))
    t = sp.Symbol('t')
    gens, eqs = list(x) + [t], [t * E.det() - 1]
    span = cell.span()
    tables = ct.support_tables(cell.chain(span))
    if turned:
        z, zi = sp.symbols('z zi')
        gens += [z, zi]
        eqs.append(z * zi - 1)
        g = lambda m: sp.diag(1, zi ** m, z ** m, 1) if m >= 0 else sp.diag(1, z ** -m, zi ** -m, 1)
        site = ct.site_table(cell.chain(1), 0)
        L0 = traceless(site)
        for M in (L0, target(L0, relation)):
            eqs += [e for e in (g(1) * M - M * g(1)).applyfunc(sp.expand) if e != 0]
        maps = [g(m) * E * g(-m) for m in range(span)]
    else:
        maps = [E] * span
    for S, LS in tables.items():
        if turned:
            gS = sp.kronecker_product(*[g(1)] * len(S))
            for M in (LS, target(LS, relation)):
                eqs += [e for e in (gS * M - M * gS).applyfunc(sp.expand) if e != 0]
        ES = sp.kronecker_product(*[maps[l] for l in S])
        eqs += [e for e in (ES * LS - target(LS, relation) * ES).applyfunc(sp.expand) if e != 0]
    return _solvable(eqs, gens)


def every_n(cell, relation, rate_free=False, max_period=4, max_calls=2000):
    """For a chain built from one repeated cell: 'UNIFORM' (one dressed map everywhere, exact), 'TURNED' (the map
    turned about Z by a fixed angle from site to site, exact), ('PERIODIC', p) for p = 2..max_period, 'STEPPED'
    (the letters' turns step from site to site without a turn of the map), or None: then an answer holds at the
    computed N only. The letters answer first (cube_letters.every_n, exact on every open chain); when they are
    silent the handshake below decides (site families or cube moves, Groebner bases; PERIODIC found only)."""
    from . import cube_letters as cl
    _reflections_only(relation)
    read_out = cl.every_n(cell, relation, rate_free=rate_free, max_period=max_period)
    if read_out is not None:
        return ('PERIODIC', read_out[1]) if read_out[0] == 'PERIODIC' else read_out[0]
    cell = cell.normalized()
    if _uniform_or_turned(cell, relation, rate_free, turned=False):
        return 'UNIFORM'
    if _uniform_or_turned(cell, relation, rate_free, turned=True):
        return 'TURNED'
    if rate_free:
        return None
    span = cell.span()
    for p in range(2, max_period + 1):
        length = p + span - 1
        setting = cell.chain(length)
        if precheck(setting, relation)[0] == 'NONE':
            continue
        tie = {l: l % p for l in range(length)}
        status, _ = cube_search(setting, relation, max_calls, tie=tie, groups=list(range(p)))
        if status == 'FOUND':
            return ('PERIODIC', p)
    return None
