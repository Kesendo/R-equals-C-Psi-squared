"""The letter cube's tables, site by site: one-site and per-window generator tables, and the 192 one-site moves.

Device B (a working name), experiments/THE_CUBE_SITE_BY_SITE.md. One site's
operators are written on the four cells |i><j| in the order (0,0), (0,1), (1,0), (1,1), stacked row-major, so
the superoperator rho -> A rho B is kron(A, B^T) and a product of site maps is a Kronecker product. The same 4x4
map in the Pauli basis {I, X, Y, Z} has column P = the image of the letter P, normalised as
T[a, b] = Tr(P_a M(P_b)) / 2; 4*T is C# LetterCube.TransferTimesFour. sigma^- = |0><1|.

Exact arithmetic throughout: sympy, entries in Q(i). Not exported from `framework` (it would load sympy on every
import); import it as `from framework import cube_tables`.

Public API:
  Setting(n, terms, jumps)          terms: [(letters, sites, coefficient)], jumps: [(sites, A, rate)]
  normalized_terms(terms)           each term reduced to its support (identity letters removed)
  outside_reason(setting)           why B cannot read a setting, or None
  site_table(setting, l)            4x4 one-site generator in cells (fields and jumps on site l)
  support_tables(setting)           {sorted support tuple: 4^k x 4^k table in cells, per-site order}
  hamiltonian_on(setting, support)  2^k x 2^k Hamiltonian of the terms with exactly that support
  to_pauli(M), to_cells(T)          basis change of a one-site 4x4 map
  transfer_times_four(jump, rate)   4x the Pauli transfer matrix of one jump (C# TransferTimesFour)
  cube_moves()                      the 192 one-site moves, [(name, 4x4 cell map)]
  move_cosets(side)                 24 representatives modulo the 8 cell-diagonal moves, [(name, map)]
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from itertools import product

import sympy as sp

LETTERS = 'IXYZ'
PAULI = {
    'I': sp.eye(2),
    'X': sp.Matrix([[0, 1], [1, 0]]),
    'Y': sp.Matrix([[0, -sp.I], [sp.I, 0]]),
    'Z': sp.Matrix([[1, 0], [0, -1]]),
}
LOWER = sp.Matrix([[0, 1], [0, 0]])      # sigma^- = |0><1|
RAISE = sp.Matrix([[0, 0], [1, 0]])      # sigma^+ = |1><0|


@dataclass
class Setting:
    """A finite graph of n sites. terms: (letters, sites, coefficient), e.g. ('XZX', (0, 1, 2), 3) puts X on
    site 0, Z on 1, X on 2; a one-letter term is a field. jumps: (sites, A, rate) with A a 2x2 matrix on one site."""
    n: int
    terms: list = field(default_factory=list)
    jumps: list = field(default_factory=list)


def _in_q_i(x) -> bool:
    x = sp.nsimplify(x)
    re, im = x.as_real_imag()
    return bool(sp.sympify(re).is_rational) and bool(sp.sympify(im).is_rational)


def normalized_terms(terms):
    """Each term reduced to the sites where its letter is not I. A Pauli term is traceless exactly on those sites,
    so they are its support; an all-identity term is a constant and drops out of the generator."""
    out = []
    for letters, sites, coef in terms:
        kept = [(c, s) for c, s in zip(letters, sites) if c != 'I']
        if kept:
            out.append((''.join(c for c, _ in kept), tuple(s for _, s in kept), coef))
    return out


def outside_reason(setting: Setting):
    """Why B cannot read this setting (scope of the spec), or None."""
    for letters, sites, coef in setting.terms:
        if len(sites) != len(letters) or len(set(sites)) != len(sites) or set(letters) - set(LETTERS):
            return f"term {letters} on {sites}: letters and sites do not match"
    for letters, sites, coef in normalized_terms(setting.terms):
        if len(sites) > 3:
            return f"term {letters} on {sites}: more than three sites"
        if not _in_q_i(coef) or sp.im(sp.nsimplify(coef)) != 0:
            return f"term {letters}: coefficient {coef} is not a real rational (H must be Hermitian, in Q(i))"
    for sites, A, rate in setting.jumps:
        if len(sites) != 1:
            return f"jump on {sites}: only one-site jumps are read"
        if not (_in_q_i(rate) and sp.im(sp.nsimplify(rate)) == 0 and sp.nsimplify(rate) >= 0):
            return f"jump on {sites}: rate {rate} is not a nonnegative rational"
    for l in range(setting.n):
        if not all(_in_q_i(v) for v in site_table(setting, l)):
            return f"site {l}: its dissipator has entries outside Q(i)"
    return None


def lr(A, B):
    """rho -> A rho B on one site, in cells."""
    return sp.kronecker_product(A, B.T)


def dissipator(A, rate):
    AdA = A.H * A
    return rate * (lr(A, A.H) - sp.Rational(1, 2) * (lr(AdA, sp.eye(2)) + lr(sp.eye(2), AdA)))


def site_table(setting: Setting, l: int):
    """The one-site generator of site l: -i[h_l, .] for its one-letter terms, plus its jumps."""
    h = sp.zeros(2, 2)
    for letters, sites, coef in normalized_terms(setting.terms):
        if len(sites) == 1 and sites[0] == l:
            h += sp.nsimplify(coef) * PAULI[letters]
    L = -sp.I * (lr(h, sp.eye(2)) - lr(sp.eye(2), h))
    for sites, A, rate in setting.jumps:
        if tuple(sites) == (l,):
            L += dissipator(sp.Matrix(A), sp.nsimplify(rate))
    return L.applyfunc(sp.expand)


def _per_site_perm(k: int):
    """Row-major vec of a 2^k x 2^k matrix (i_1..i_k, j_1..j_k) -> per-site order ((i_1 j_1), ..., (i_k j_k))."""
    d = 2 ** k
    P = sp.zeros(d * d, d * d)
    for bits in product(range(2), repeat=2 * k):
        i, j = bits[:k], bits[k:]
        row_major = int(''.join(map(str, i)), 2) * d + int(''.join(map(str, j)), 2)
        per_site = 0
        for a, b in zip(i, j):
            per_site = per_site * 4 + 2 * a + b
        P[per_site, row_major] = 1
    return P


def hamiltonian_on(setting: Setting, support: tuple):
    """The Hamiltonian of all terms whose support is exactly `support` (sorted), on those sites in that order."""
    support = tuple(sorted(support))
    k = len(support)
    H = sp.zeros(2 ** k, 2 ** k)
    for letters, sites, coef in normalized_terms(setting.terms):
        if tuple(sorted(sites)) != support:
            continue
        on = dict(zip(sites, letters))
        H += sp.nsimplify(coef) * sp.kronecker_product(*[PAULI[on[s]] for s in support])
    return H


def support_tables(setting: Setting) -> dict:
    """{support: table} for every support set of the Hamiltonian with at least two sites, terms summed."""
    supports = sorted({tuple(sorted(sites)) for _, sites, _ in normalized_terms(setting.terms) if len(sites) >= 2})
    out = {}
    for S in supports:
        k = len(S)
        H = hamiltonian_on(setting, S)
        Id = sp.eye(2 ** k)
        Lrm = -sp.I * (sp.kronecker_product(H, Id) - sp.kronecker_product(Id, H.T))
        P = _per_site_perm(k)
        out[S] = (P * Lrm * P.T).applyfunc(sp.expand)
    return out


# Pauli basis vectors in cells: column b is vec(P_b), row-major
_B = sp.Matrix.hstack(*[sp.Matrix([PAULI[c][0, 0], PAULI[c][0, 1], PAULI[c][1, 0], PAULI[c][1, 1]])
                        for c in LETTERS])


def to_pauli(M):
    """A one-site 4x4 map in cells -> its matrix T in the Pauli basis, T[a, b] = Tr(P_a M(P_b)) / 2."""
    return (_B.H * M * _B / 2).applyfunc(sp.expand)


def to_cells(T):
    return (_B * T * _B.H / 2).applyfunc(sp.expand)


def transfer_times_four(jump, rate):
    """4x the Pauli transfer matrix of rate·D[jump]: the numbers of C# LetterCube.TransferTimesFour."""
    return (4 * to_pauli(dissipator(sp.Matrix(jump), sp.nsimplify(rate)))).applyfunc(sp.expand)


def _key(M):
    for v in M:
        if v != 0:
            return tuple(sp.nsimplify(x / v) for x in M)
    return tuple(M)


@lru_cache(maxsize=None)
def cube_moves():
    """The 192 one-site moves (C# LetterCube.Closure): generated by X·rho, Z·rho, rho·X, rho·Z, the Hadamard on
    both sides, the quarter-turn S·rho·S^dagger and the transpose, modulo a global phase. [(name, cell map)],
    each name the shortest word of generators."""
    S = sp.diag(1, sp.I)
    T = sp.zeros(4, 4)
    for r, c in ((0, 0), (1, 2), (2, 1), (3, 3)):
        T[r, c] = 1
    gens = [('X·', lr(PAULI['X'], sp.eye(2))), ('Z·', lr(PAULI['Z'], sp.eye(2))),
            ('·X', lr(sp.eye(2), PAULI['X'])), ('·Z', lr(sp.eye(2), PAULI['Z'])),
            ('H', lr(PAULI['X'] + PAULI['Z'], PAULI['X'] + PAULI['Z']) / 2), ('S', lr(S, S.H)), ('θ', T)]
    seen = {_key(sp.eye(4)): ('1', sp.eye(4))}
    frontier = [('1', sp.eye(4))]
    while frontier:
        nxt = []
        for name, M in frontier:
            for gname, g in gens:
                P = (g * M).applyfunc(sp.nsimplify)
                k = _key(P)
                if k not in seen:
                    word = gname if name == '1' else f'{gname} {name}'
                    seen[k] = (word, P)
                    nxt.append((word, P))
        frontier = nxt
    return tuple(seen.values())


def is_cell_diagonal(M) -> bool:
    return all(M[i, j] == 0 for i in range(4) for j in range(4) if i != j)


@lru_cache(maxsize=None)
def move_cosets(side='right'):
    """24 representatives of the 192 moves modulo the 8 cell-diagonal ones: right cosets {w·P} (the weight after
    the move) or left cosets {P·w} (the weight before it). The 8 cell-diagonal moves are not a normal subgroup,
    so the two lists differ."""
    reps = []
    for name, P in cube_moves():
        same = (lambda Q: (P * Q.inv())) if side == 'right' else (lambda Q: (Q.inv() * P))
        if not any(is_cell_diagonal(same(Q).applyfunc(sp.nsimplify)) for _, Q in reps):
            reps.append((name, P))
    return tuple(reps)
