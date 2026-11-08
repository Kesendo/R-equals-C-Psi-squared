"""The shape of F158's far element: is it a product, on which sheets does it sit, and does the pairing hold there?

The question it answers: a row pairs through a one-dimensional far end (F158, f158_end_count), so every mirror
carries the identity onto that end's one element W. What shape is W?

  operator Schmidt rank   the rank of W's coefficient matrix across a split of the sites. A mirror of rank r sends a
                          product to rank at most r, so W's rank bounds every mirror's from below, and where W is
                          no product no mirror that factorises over the sites (linear or antilinear, after a frame
                          change per site, a transposition or a relabelling by conjugation) carries the row: the
                          lemma of F116 §4 (docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md) and, in rank form, of
                          finding 1 of experiments/THE_CUBE_SITE_BY_SITE.md. A mirror conditioned on the bits of a
                          site is not covered.
  kind patterns           per string the bits, site by site, 1 where the letter anticommutes with that site's axis
                          letter (Z by default: the per-site share of the cube's k_Z), 0 where it commutes.
  homogeneity group       the Pauli strings g under which every string of W has one sign, g·W·g = ±W. Every jump is
                          in it (W anticommutes with each), and, the far end being one-dimensional, so is every g
                          that commutes or anticommutes with all of H (then g·W·g is again a far element, so ±W).
                          What lies beyond those is a SHEET: g splits H into H_+ (terms commuting with g) and H_−
                          (anticommuting), and W commutes with each part alone, since g·H·g = H_+ − H_− and
                          g·W·g = ±W turn [H, W] = 0 into [H_+ − H_−, W] = 0 as well. So W stays in the far end of
                          every member of the pencil H_+ + λ·H_−; whether the row still PAIRS there (the two ends
                          staying equal) is read by `pencil`, a count and not a consequence.

Public API:
  operator_schmidt_rank(element, n, left)   exact rank over Q(i) across left | rest
  kind_patterns(element, n, axes=None)      sorted distinct patterns, axes a string of one letter per site
  homogeneity_group(element, n)             every Pauli string g with g·W·g = ±W (2^dim of them)
  sheets(n, h_terms, jumps, element)        Sheets(dim_group, dim_given, extra, splits)
  pencil(n, h_terms, jumps, g, lam)         end_count of H_+ + λ·H_− for the split by g
  far_element_shape(n, h_terms, jumps)      FarElementShape, or None unless near = far = 1
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import sympy as sp
from sympy.polys.domains import QQ_I
from sympy.polys.matrices import DomainMatrix

from .f158_end_count import end_count, _gauss

__all__ = ["Sheets", "FarElementShape", "operator_schmidt_rank", "kind_patterns", "homogeneity_group", "sheets",
           "pencil", "far_element_shape"]


def _anticommutes(s, t):
    return sum(1 for a, b in zip(s, t) if a != 'I' and b != 'I' and a != b) % 2 == 1


def operator_schmidt_rank(element, n, left):
    """Rank over Q(i) of c[left letters, right letters] for W = sum c_s s."""
    left = sorted(set(left))
    if any(l < 0 or l >= n for l in left):
        raise ValueError("a site of the split is outside the chain")
    right = [l for l in range(n) if l not in left]
    rows, cols, entries = {}, {}, {}
    for s, c in element.items():
        if c == 0:
            continue
        r = rows.setdefault(''.join(s[l] for l in left), len(rows))
        k = cols.setdefault(''.join(s[l] for l in right), len(cols))
        entries[(r, k)] = _gauss(c)
    if not entries:
        return 0
    M = DomainMatrix([[entries.get((r, k), QQ_I(0, 0)) for k in range(len(cols))] for r in range(len(rows))],
                     (len(rows), len(cols)), QQ_I)
    return M.rank()


def kind_patterns(element, n, axes=None):
    """Distinct per-site kind bits of W's strings: 1 where the letter anticommutes with the site's axis letter."""
    axes = axes or 'Z' * n
    if len(axes) != n:
        raise ValueError("one axis letter per site")
    def bit(letter, axis):
        return '1' if letter != 'I' and letter != axis else '0'
    return sorted({''.join(bit(s[l], axes[l]) for l in range(n)) for s, c in element.items() if c != 0})


def _bits(s):
    return [1 if c in 'XY' else 0 for c in s] + [1 if c in 'ZY' else 0 for c in s]


def _letters(v, n):
    return ''.join('IXZY'[v[i] + 2 * v[n + i]] for i in range(n))


def homogeneity_group(element, n):
    """Every Pauli string g such that all strings of W share one commutation sign with g: g anticommutes with
    s + s_0 never, for every string s of W, a linear condition over GF(2) on g's 2N bits (the symplectic form),
    so the group is a nullspace, enumerated from its basis."""
    strings = [s for s, c in element.items() if c != 0]
    if not strings:
        raise ValueError("an empty element has no homogeneity group")
    s0 = _bits(strings[0])
    rows = []
    for s in strings[1:]:
        d = [(a + b) % 2 for a, b in zip(_bits(s), s0)]
        rows.append(d[n:] + d[:n])          # <g, d> = g_x·d_z + g_z·d_x, so the row is d with its halves swapped
    # reduced row echelon form over GF(2), then a nullspace basis
    m, pivots, r = [list(x) for x in rows], [], 0
    for c in range(2 * n):
        p = next((i for i in range(r, len(m)) if m[i][c]), None)
        if p is None:
            continue
        m[r], m[p] = m[p], m[r]
        for i in range(len(m)):
            if i != r and m[i][c]:
                m[i] = [(a + b) % 2 for a, b in zip(m[i], m[r])]
        pivots.append(c)
        r += 1
    basis = []
    for free in (c for c in range(2 * n) if c not in pivots):
        v = [0] * (2 * n)
        v[free] = 1
        for k, pc in enumerate(pivots):
            v[pc] = m[k][free]
        basis.append(v)
    out = []
    for coeffs in product((0, 1), repeat=len(basis)):
        v = [0] * (2 * n)
        for a, b in zip(coeffs, basis):
            if a:
                v = [(x + y) % 2 for x, y in zip(v, b)]
        out.append(_letters(v, n))
    return out


def _rank2(vectors):
    rows = [list(v) for v in vectors]
    r = 0
    for c in range(len(rows[0]) if rows else 0):
        p = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                rows[i] = [(a + b) % 2 for a, b in zip(rows[i], rows[r])]
        r += 1
    return r


@dataclass
class Sheets:
    dim_group: int    # dim of the homogeneity group (a subspace of the 2N-bit string space)
    dim_given: int    # dim of the part every far element has: the jumps and the g that do not split H
    extra: int        # dim_group - dim_given: the sheets
    splits: list      # one (g, H_- terms) per distinct split of H by a g in the group, the given ones included


def sheets(n, h_terms, jumps, element):
    """The homogeneity group of W measured against what every far element of this row is given."""
    G = homogeneity_group(element, n)
    terms = [P for P, c in h_terms.items() if c != 0 and set(P) != {'I'}]
    given = [j for j in jumps] + [g for g in G if len({_anticommutes(g, t) for t in terms}) <= 1]
    dim_group = len(G).bit_length() - 1
    dim_given = _rank2([_bits(g) for g in given]) if given else 0
    splits, seen = [], set()
    for g in G:
        signs = tuple(_anticommutes(g, t) for t in terms)
        if 0 < sum(signs) < len(signs) and signs not in seen:
            seen.add(signs)
            splits.append((g, sorted(t for t, s in zip(terms, signs) if s)))
    return Sheets(dim_group, dim_given, dim_group - dim_given, splits)


def pencil(n, h_terms, jumps, g, lam):
    """end_count of H_+ + λ·H_−, the terms of H split by whether they anticommute with the Pauli string g."""
    lam = sp.nsimplify(lam)
    h = {P: (lam * c if _anticommutes(g, P) else c) for P, c in h_terms.items()}
    return end_count(n, h, jumps)


@dataclass
class FarElementShape:
    element: dict             # the far end's one element, {string: coefficient}
    single_site_ranks: tuple  # operator Schmidt rank across {l} | rest, l = 0 .. N-1
    product: bool             # rank one across every single-site split, so a product operator
    patterns: list            # kind_patterns against Z on every site
    sheets: Sheets


def far_element_shape(n, h_terms, jumps):
    """The shape of the far element where the row pairs through one (near = far = 1), else None."""
    ec = end_count(n, h_terms, jumps)
    if not (ec.near == 1 and ec.far == 1):
        return None
    w = {s: c for s, c in ec.far_basis[0].items() if c != 0}
    ranks = tuple(operator_schmidt_rank(w, n, [l]) for l in range(n))
    return FarElementShape(w, ranks, all(r == 1 for r in ranks), kind_patterns(w, n), sheets(n, h_terms, jumps, w))
