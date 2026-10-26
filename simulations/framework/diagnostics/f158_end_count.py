"""F158's two-end count for Pauli-string jumps: does the spectrum pair about -sigma, and what carries it?

The question it answers: H a sum of Pauli strings with rational coefficients, every jump a Pauli string at a
positive rate. Is the spectrum of L palindromic about -sigma (sigma the sum of the rates), and if so, which
operator carries it? This is the hand-over after the Theodolite's silence: the letters and the pre-check decide
whether a SITE-LOCAL element exists; this decides whether the spectrum pairs at all.

F158 (docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md): the palindrome holds exactly when dim ker L =
dim ker(L + 2 sigma). For Pauli-string jumps each end is the commutator with H on a span of Pauli strings:
  near end   the DARK strings, commuting with every jump;
  far end    the LIT strings, anticommuting with every jump.
So both counts are kernels of [H, .] restricted to a span, decided here exactly by a rank over Q(i) (sympy
DomainMatrix), no prime, no float, no upper bound. The C# twin is compute/MirrorWorld/EndCount.cs (ranks modulo
primes, upper bounds, elements lifted by rational reconstruction); this module is its small-N exact sibling.

An element U of the far kernel that is invertible carries the palindrome by ONE-SIDED multiplication: U commutes
with H and anticommutes with every jump, so S(rho) = rho·U keeps the Hamiltonian part and turns the dissipator D into
-D - 2 sigma, S^-1·L·S = -L^† - 2 sigma (the line reflection, relation 'refl'). Conjugation by U is a symmetry of L,
not a reflection. The far kernel's basis is returned, and when one basis vector squares to a nonzero multiple of the
identity it is named as the element.

Public API:
  pauli_jump_letter(A)                 'X' | 'Y' | 'Z' if the 2x2 matrix A is a Pauli letter up to a scalar, else None
  end_count(n, h_terms, jumps)         EndCount(near, far, pairs, far_basis, element)
      h_terms: {pauli string: rational coefficient}; jumps: list of Pauli strings (rates only need be positive)
  end_count_of_setting(setting)        the same for a framework.cube_tables.Setting (zero-rate jumps dropped), or
                                       None when a jump is not a Pauli letter (sigma^+-, thermal: F158's count does
                                       not apply)
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import sympy as sp
from sympy.polys.domains import QQ_I
from sympy.polys.matrices import DomainMatrix

from ..pauli import pauli_product

__all__ = ["EndCount", "pauli_jump_letter", "end_count", "end_count_of_setting"]

_LETTERS = 'IXYZ'
_I_POW = (QQ_I(1, 0), QQ_I(0, 1), QQ_I(-1, 0), QQ_I(0, -1))


@dataclass
class EndCount:
    near: int                 # dim ker L       (dark strings commuting with H)
    far: int                  # dim ker(L+2s)   (lit strings commuting with H)
    pairs: bool               # near == far: the palindrome about -sigma (F158)
    far_basis: list           # basis of the far kernel, each {string: coefficient}
    element: dict | None      # a basis vector squaring to a nonzero scalar, if one exists


def pauli_jump_letter(A):
    """The Pauli letter A is proportional to, or None."""
    A = sp.Matrix(A)
    for letter, P in (('X', sp.Matrix([[0, 1], [1, 0]])), ('Y', sp.Matrix([[0, -sp.I], [sp.I, 0]])),
                      ('Z', sp.Matrix([[1, 0], [0, -1]]))):
        k = next((A[i] / P[i] for i in range(4) if P[i] != 0), None)
        if k is not None and k != 0 and (A - k * P).applyfunc(sp.simplify).is_zero_matrix:
            return letter
    return None


def _anticommutes(s, t):
    return sum(1 for a, b in zip(s, t) if a != 'I' and b != 'I' and a != b) % 2 == 1


def _span(n, jumps, lit):
    """Every Pauli string commuting (lit False) or anticommuting (lit True) with every jump."""
    return [s for s in (''.join(w) for w in product(_LETTERS, repeat=n))
            if all(_anticommutes(s, j) == lit for j in jumps)]


def _gauss(c):
    c = sp.nsimplify(c)
    re, im = c.as_real_imag()
    return QQ_I(sp.Rational(re), sp.Rational(im))


def _kernel(h_terms, span):
    """Basis of {sum c_s s : s in span, [H, sum c_s s] = 0}, exactly."""
    if not span:
        return []
    h = [(P, _gauss(c)) for P, c in h_terms.items() if c != 0 and set(P) != {'I'}]
    rows, entries = {}, {}
    for col, s in enumerate(span):
        for P, c in h:
            if not _anticommutes(P, s):
                continue
            e, u = pauli_product(P, s)           # [P, s] = 2 P s when they anticommute
            r = rows.setdefault(u, len(rows))
            entries[(r, col)] = entries.get((r, col), QQ_I(0, 0)) + 2 * c * _I_POW[e]
    if not rows:
        return [{s: sp.Integer(1)} for s in span]
    M = DomainMatrix([[entries.get((r, c), QQ_I(0, 0)) for c in range(len(span))] for r in range(len(rows))],
                     (len(rows), len(span)), QQ_I)
    basis = M.nullspace().to_Matrix()
    out = []
    for i in range(basis.rows):
        v = {s: sp.nsimplify(basis[i, k]) for k, s in enumerate(span) if basis[i, k] != 0}
        out.append(v)
    return out


def _square_scalar(v, n):
    """The scalar v^2 is, or None when v^2 is not a multiple of the identity."""
    sq = {}
    for s, a in v.items():
        for t, b in v.items():
            e, u = pauli_product(s, t)
            sq[u] = sq.get(u, 0) + a * b * sp.I ** e
    sq = {u: sp.expand(c) for u, c in sq.items() if sp.expand(c) != 0}
    ident = 'I' * n
    if set(sq) == {ident}:
        return sq[ident]
    return None


def end_count(n, h_terms, jumps):
    """F158's two ends for H = sum h_terms and Pauli-string jumps, exactly."""
    far_basis = _kernel(h_terms, _span(n, jumps, lit=True))
    near = len(_kernel(h_terms, _span(n, jumps, lit=False)))
    element = next((v for v in far_basis if _square_scalar(v, n) not in (None, 0)), None)
    return EndCount(near, len(far_basis), near == len(far_basis), far_basis, element)


def end_count_of_setting(setting):
    """end_count for a cube_tables.Setting, or None when a jump is not a Pauli letter or a rate is not positive."""
    from .. import cube_tables as ct
    n = setting.n
    jumps = []
    for sites, A, rate in setting.jumps:
        if sp.nsimplify(rate) == 0:
            continue                             # a rate set to zero removes that jump's condition (F158)
        letter = pauli_jump_letter(A)
        if letter is None or len(sites) != 1 or not sp.nsimplify(rate) > 0:
            return None
        w = ['I'] * n
        w[sites[0]] = letter
        jumps.append(''.join(w))
    h = {}
    for letters, sites, coef in ct.normalized_terms(setting.terms):
        w = ['I'] * n
        for c, s in zip(letters, sites):
            w[s] = c
        key = ''.join(w)
        h[key] = h.get(key, 0) + sp.nsimplify(coef)
    return end_count(n, h, jumps)
