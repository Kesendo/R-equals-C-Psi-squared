"""The letter layer of Device B (a working name): site-local mirrors read off the letters, in closed form.

See experiments/THE_CUBE_SITE_BY_SITE.md. The element is a sandwich E(rho) = A·rho^(T)·B with A = ⊗A_l and
B = ⊗B_l, every site factor monomial: diagonal 'd' = diag(1, x) or anti-diagonal 'a' = [[0, x], [1, 0]] = X·diag(1, x).
On the letters such a factor is a pure turn: it multiplies sigma^+ by a power of x or swaps sigma^+ and sigma^- with
phases, and keeps or flips Z. Every condition is therefore binomial in the turns, a monomial equal to a constant:
  site l      E_l·L_l° = 𝒯(L_l°)·E_l (heights, arrows and field of the site);
  support S   A_S·h_S^(T)·A_S⁻¹ = s·h_S and B_S⁻¹·h_S^(T)·B_S = s·h_S, the left and right factors apart,
              s = -1 for ('pi', plain) and ('refl', transposed), +1 otherwise (the commutator part of the sector
              split, exactly: the bond terms are traceless on every site, so no constant can enter).
A binomial system is solved by integer elimination on the exponents (Euclid on the columns) and radicals of the
constants: exact, in closed form, without a Groebner basis. Every N follows from the same equations with a step per
site, x_l = x_0·r^l and y_l = y_0·s^l: two positions of a window fix the step, so an answer holds on every open chain.

Public API:
  factor(kind, x), site_map(kA, x, kB, y, transposed)
  solve_binomials(eqs, unknowns)            one solution (free turns set to 1) or None
  read(setting, relation)                   Reading | None   (None: the letters are silent)
  verify(setting, relation, maps)           exact check of every site and support equation
  describe(reading)                         the equations in words
  every_n(cell, relation, rate_free, max_period)
        ('UNIFORM', info) | ('TURNED', info) | ('PERIODIC', p, info) | ('STEPPED', info) | None
  turn_of(kinds, transposed, r, s)          the turn about Z carrying one site's map onto the next, or None
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import product

import sympy as sp

from .. import cube_tables as ct
from . import cube_reflector as cr

KINDS = ('d', 'a')
_T = sp.zeros(4, 4)
for _r, _c in ((0, 0), (1, 2), (2, 1), (3, 3)):
    _T[_r, _c] = 1


def factor(kind, x):
    return sp.diag(1, x) if kind == 'd' else sp.Matrix([[0, x], [1, 0]])


def site_map(kA, x, kB, y, transposed):
    """rho -> A·rho^(T)·B on one site, in cells."""
    E = ct.lr(factor(kA, x), factor(kB, y))
    return E * _T if transposed else E


def h_sign(relation, transposed):
    return -1 if (relation, transposed) in (('pi', False), ('refl', True)) else 1


class Silent(Exception):
    """An equation with more than two monomials: the letters do not read this assignment."""


def _binomials(exprs, unknowns, origin):
    """Each expression = 0 with the unknowns nonzero, as binomials ({unknown: exponent}, constant, origin) meaning
    prod u^e = constant; None when an expression is a single monomial (no nonzero solution)."""
    out = []
    for e in exprs:
        num = sp.expand(sp.fraction(sp.together(sp.expand(e)))[0])
        if num == 0:
            continue
        terms = sp.Poly(num, *unknowns).terms()
        if len(terms) == 1:
            return None
        if len(terms) > 2:
            raise Silent(f'{origin}: {len(terms)} monomials')
        (m1, c1), (m2, c2) = terms
        exps = {u: a - b for u, a, b in zip(unknowns, m1, m2) if a != b}
        out.append((exps, sp.nsimplify(-c2 / c1), origin))
    return out


def site_binomials(pieces, relation, transposed, kA, kB, x, y, origin):
    """The site equation for every traceless piece (one piece: the whole site table; rate-free: one per channel)."""
    E = site_map(kA, x, kB, y, transposed)
    out = []
    for L0 in pieces:
        b = _binomials(list(E * L0 - cr.target(L0, relation) * E), [x, y], origin)
        if b is None:
            return None
        out += b
    return out


def support_binomials(h, relation, transposed, kinds, xs, side, origin):
    """side 'A': A_S·h^(T)·A_S⁻¹ = s·h; side 'B': B_S⁻¹·h^(T)·B_S = s·h."""
    M = sp.kronecker_product(*[factor(k, x) for k, x in zip(kinds, xs)])
    hT = h.T if transposed else h
    lhs = M * hT * M.inv() if side == 'A' else M.inv() * hT * M
    return _binomials(list(lhs - h_sign(relation, transposed) * h), list(xs), origin)


# ---------------------------------------------------------------- binomial systems

def _simp(c):
    return sp.nsimplify(sp.simplify(sp.radsimp(c)))


def solve_binomials(eqs, unknowns):
    """prod u^e = c for every ({u: e}, c, origin): one exact solution {u: value}, every free turn set to 1 and every
    root principal, or None. Euclid on the exponent columns keeps the row operations unimodular; the constants are
    carried multiplicatively; the solution is checked against every equation before it is returned."""
    n = len(unknowns)
    eqs = [(e, sp.nsimplify(sp.sympify(c)), origin) for e, c, origin in eqs]
    rows = [[[e.get(u, 0) for u in unknowns], c] for e, c, _ in eqs]
    r = 0
    pivots = []
    for col in range(n):
        while True:
            nz = [i for i in range(r, len(rows)) if rows[i][0][col] != 0]
            if not nz:
                break
            i0 = min(nz, key=lambda i: abs(rows[i][0][col]))
            rows[r], rows[i0] = rows[i0], rows[r]
            clean = True
            for i in range(r + 1, len(rows)):
                a = rows[i][0][col]
                if a:
                    q = a // rows[r][0][col]
                    rows[i] = [[ai - q * ar for ai, ar in zip(rows[i][0], rows[r][0])],
                               _simp(rows[i][1] / rows[r][1] ** q)]
                    if rows[i][0][col]:
                        clean = False
            if clean:
                break
        if any(rows[i][0][col] for i in range(r, len(rows))):
            pivots.append((r, col))
            r += 1
    for i in range(r, len(rows)):            # rows emptied by the elimination: the constant must be 1
        if _simp(rows[i][1] - 1) != 0:
            return None
    val = {j: sp.Integer(1) for j in range(n)}
    for ri, col in reversed(pivots):
        e, c = rows[ri]
        rest = sp.Mul(*[val[j] ** e[j] for j in range(col + 1, n) if e[j]])
        val[col] = _simp((c / rest) ** sp.Rational(1, e[col]))
    sol = {unknowns[j]: val[j] for j in range(n)}
    for e, c, _ in eqs:
        if _simp(sp.Mul(*[sol[u] ** k for u, k in e.items()]) - c) != 0:
            return None
    return sol


# ---------------------------------------------------------------- reading a setting

@dataclass
class Reading:
    transposed: bool
    kinds: list                 # per site (kind of A_l, kind of B_l)
    x: list                     # per site the turn of A_l
    y: list                     # per site the turn of B_l
    maps: list                  # per site the 4x4 cell map
    equations: list = field(default_factory=list)


def read(setting, relation):
    """A monomial sandwich carrying L onto 𝒯(L) - c, read from the letters, or None when the letters are silent
    (no such sandwich exists, or an equation had more than two monomials)."""
    n = setting.n
    pieces = [[cr.traceless(ct.site_table(setting, l))] for l in range(n)]
    supports = list(ct.support_tables(setting))
    hams = {S: ct.hamiltonian_on(setting, S) for S in supports}
    xs, ys = sp.symbols(f'x0:{n}'), sp.symbols(f'y0:{n}')
    done_at = {l: [S for S in supports if max(S) == l] for l in range(n)}
    for transposed in (False, True):
        allowed = []
        for l in range(n):
            opts = []
            for kA, kB in product(KINDS, KINDS):
                try:
                    b = site_binomials(pieces[l], relation, transposed, kA, kB, xs[l], ys[l], ('site', l))
                except Silent:
                    continue
                if b is not None:
                    opts.append(((kA, kB), b))
            allowed.append(opts)
        if any(not o for o in allowed):
            continue
        choice = {}

        def system(upto):
            eqs = []
            for l in range(upto + 1):
                eqs += allowed[l][choice[l]][1]
                for S in done_at[l]:
                    for side, j, us in (('A', 0, [xs[k] for k in S]), ('B', 1, [ys[k] for k in S])):
                        kinds = [allowed[k][choice[k]][0][j] for k in S]
                        b = support_binomials(hams[S], relation, transposed, kinds, us, side, ('support', S, side))
                        if b is None:
                            return None
                        eqs += b
            return eqs

        def rec(l):
            if l == n:
                return True
            for ci in range(len(allowed[l])):
                choice[l] = ci
                try:
                    eqs = system(l)
                except Silent:
                    continue
                if eqs is None:
                    continue
                if solve_binomials(eqs, [u for k in range(l + 1) for u in (xs[k], ys[k])]) is None:
                    continue
                if rec(l + 1):
                    return True
            choice.pop(l, None)
            return False

        if rec(0):
            eqs = system(n - 1)
            sol = solve_binomials(eqs, [u for k in range(n) for u in (xs[k], ys[k])])
            kinds = [allowed[l][choice[l]][0] for l in range(n)]
            x, y = [sol[v] for v in xs], [sol[v] for v in ys]
            maps = [site_map(kinds[l][0], x[l], kinds[l][1], y[l], transposed) for l in range(n)]
            return Reading(transposed, kinds, x, y, maps, eqs)
    return None


def verify(setting, relation, maps):
    """Exact: every site equation and every support-set equation of the sector split, and every map invertible."""
    for l in range(setting.n):
        L0 = cr.traceless(ct.site_table(setting, l))
        if not (maps[l] * L0 - cr.target(L0, relation) * maps[l]).applyfunc(sp.simplify).is_zero_matrix:
            return False
    for S, LS in ct.support_tables(setting).items():
        ES = sp.kronecker_product(*[maps[l] for l in S])
        if not (ES * LS - cr.target(LS, relation) * ES).applyfunc(sp.simplify).is_zero_matrix:
            return False
    return all(sp.simplify(M.det()) != 0 for M in maps)


def _monomial(e):
    return '·'.join(f'{u}' if k == 1 else f'{u}^{k}' for u, k in sorted(e.items(), key=lambda t: str(t[0])))


def describe(reading):
    """The element and the equations it was read from, in words."""
    form = 'A·ρᵀ·B' if reading.transposed else 'A·ρ·B'
    out = [f'ρ ↦ {form} with A = ⊗A_l, B = ⊗B_l']
    for l, ((kA, kB), x, y) in enumerate(zip(reading.kinds, reading.x, reading.y)):
        fa = 'diag(1, x)' if kA == 'd' else '[[0, x], [1, 0]]'
        fb = 'diag(1, y)' if kB == 'd' else '[[0, y], [1, 0]]'
        out.append(f'site {l}: A = {fa} with x = {x}; B = {fb} with y = {y}')
    for e, c, origin in reading.equations:
        where = f'site {origin[1]}' if origin[0] == 'site' else \
            f"support {origin[1]}, {'left' if origin[2] == 'A' else 'right'} factor"
        out.append(f'  {where}: {_monomial(e)} = {c}')
    return out


# ---------------------------------------------------------------- every N

def _cell_pieces(cell, relation, rate_free):
    """The traceless site pieces: the whole one-site table, or with rate_free every channel and the fields apart."""
    one = cell.chain(1)
    if not rate_free:
        return [cr.traceless(ct.site_table(one, 0))]
    pieces = [cr.traceless(ct.dissipator(sp.Matrix(A), 1)) for A, _ in cell.jumps]
    if cell.fields:
        pieces.append(cr.traceless(ct.site_table(ct.Setting(1, [(lt, (0,), c) for lt, c in cell.fields], []), 0)))
    return pieces


def _chain_binomials(cell, relation, transposed, kinds, length, xs, ys, pieces):
    s = cell.chain(length)
    eqs = []
    for l in range(length):
        b = site_binomials(pieces, relation, transposed, kinds[l][0], kinds[l][1], xs[l], ys[l], ('site', l))
        if b is None:
            return None
        eqs += b
    for S in ct.support_tables(s):
        h = ct.hamiltonian_on(s, S)
        for side, j, us in (('A', 0, [xs[k] for k in S]), ('B', 1, [ys[k] for k in S])):
            b = support_binomials(h, relation, transposed, [kinds[k][j] for k in S], us, side, ('support', S, side))
            if b is None:
                return None
            eqs += b
    return eqs


def _substitute(eqs, subst):
    out = []
    for e, c, origin in eqs:
        ne = {}
        for u, k in e.items():
            for v, m in subst[u].items():
                ne[v] = ne.get(v, 0) + k * m
        out.append(({v: k for v, k in ne.items() if k}, c, origin))
    return out


def turn_of(kinds, transposed, r, s):
    """z with u = diag(1, z) and Ad_u∘E_l∘Ad_u⁻¹ = E_{l+1} up to a factor, the sandwich's turns stepping by r (left)
    and s (right); candidates are the square roots of the steps and their inverses, each checked exactly. None if no
    turn about Z carries one site's map onto the next."""
    kA, kB = kinds
    x, y = sp.symbols('x y', nonzero=True)
    E0 = site_map(kA, x, kB, y, transposed)
    E1 = site_map(kA, x * r, kB, y * s, transposed)
    cands = {_simp(root) for step in (r, s, 1 / r, 1 / s) for root in (sp.sqrt(step), -sp.sqrt(step))}
    for z in sorted(cands, key=str):
        u, ui = sp.diag(1, z), sp.diag(1, 1 / z)
        M = (ct.lr(u, u.H) * E0 * ct.lr(ui, ui.H)).applyfunc(sp.simplify)
        ratio, ok = None, True
        for a, b in zip(M, E1):
            a, b = sp.simplify(a), sp.simplify(b)
            if a == 0 and b == 0:
                continue
            if a == 0 or b == 0:
                ok = False
                break
            q = sp.simplify(a / b)
            if ratio is None:
                ratio = q
            elif sp.simplify(q - ratio) != 0:
                ok = False
                break
        if ok:
            return z
    return None


def every_n(cell, relation, rate_free=False, max_period=4):
    """A chain built from one repeated cell, read from the letters for every length (open chains):
      ('UNIFORM', info)        one sandwich on every site;
      ('TURNED', info)         the map turned about Z by info['turn'] from site to site (exact);
      ('PERIODIC', p, info)    kinds and turns repeating with period p (2..max_period);
      ('STEPPED', info)        the turns step by info['r'], info['s'], which no turn about Z produces;
      None                     the letters are silent.
    rate_free: the element must hold at every rate of every channel (each channel apart)."""
    cell = cell.normalized()
    span = cell.span()
    pieces = _cell_pieces(cell, relation, rate_free)
    u_, v_ = sp.symbols('u_ v_')

    def site_ok(transposed, kA, kB):
        try:
            return site_binomials(pieces, relation, transposed, kA, kB, u_, v_, ('site', 0)) is not None
        except Silent:
            return False

    allowed = {tr: [k for k in product(KINDS, KINDS) if site_ok(tr, *k)] for tr in (False, True)}
    length = span + 1
    xs, ys = sp.symbols(f'x0:{length}'), sp.symbols(f'y0:{length}')
    x0, y0, r, s = sp.symbols('X0 Y0 r s')
    geo = {}
    for l in range(length):
        geo[xs[l]] = {x0: 1, r: l} if l else {x0: 1}
        geo[ys[l]] = {y0: 1, s: l} if l else {y0: 1}
    systems = []
    for transposed in (False, True):
        for kinds in allowed[transposed]:
            try:
                eqs = _chain_binomials(cell, relation, transposed, [kinds] * length, length, xs, ys, pieces)
            except Silent:
                continue
            if eqs is not None:
                systems.append((transposed, kinds, _substitute(eqs, geo)))
    for transposed, kinds, g in systems:
        sol = solve_binomials(g + [({r: 1}, 1, 'uniform'), ({s: 1}, 1, 'uniform')], [x0, y0, r, s])
        if sol is not None:
            return ('UNIFORM', {'transposed': transposed, 'kinds': kinds, 'x': sol[x0], 'y': sol[y0]})
    stepped = None
    for transposed, kinds, g in systems:
        sol = solve_binomials(g, [x0, y0, r, s])
        if sol is None:
            continue
        info = {'transposed': transposed, 'kinds': kinds, 'x': sol[x0], 'y': sol[y0], 'r': sol[r], 's': sol[s]}
        z = turn_of(kinds, transposed, sol[r], sol[s])
        if z is not None:
            info['turn'] = z
            return ('TURNED', info)
        stepped = stepped or info
    for p in range(2, max_period + 1):
        length = p + span - 1
        xs, ys = sp.symbols(f'x0:{length}'), sp.symbols(f'y0:{length}')
        cx, cy = sp.symbols(f'P0:{p}'), sp.symbols(f'Q0:{p}')
        per = {}
        for l in range(length):
            per[xs[l]] = {cx[l % p]: 1}
            per[ys[l]] = {cy[l % p]: 1}
        for transposed in (False, True):
            for pattern in product(allowed[transposed], repeat=p):
                try:
                    eqs = _chain_binomials(cell, relation, transposed, [pattern[l % p] for l in range(length)],
                                           length, xs, ys, pieces)
                except Silent:
                    continue
                if eqs is None:
                    continue
                sol = solve_binomials(_substitute(eqs, per), list(cx) + list(cy))
                if sol is not None:
                    return ('PERIODIC', p, {'transposed': transposed, 'kinds': list(pattern),
                                            'x': [sol[v] for v in cx], 'y': [sol[v] for v in cy]})
    return ('STEPPED', stepped) if stepped is not None else None
