"""The letter cube: where a Pauli string sits, what it costs, and how a Hamiltonian moves it.

docs/THE_ONE_SQUARE.md §7 and §8. A Pauli string σ on N sites sits at the lattice point

    (k_Z, k_X, k_Y) = (n_X + n_Y, n_Y + n_Z, n_X + n_Z),

k_P = number of sites whose letter anticommutes with P. The points fill the even lattice
points of the tetrahedron with vertices I^N, X^N, Y^N, Z^N, and w = (k_Z + k_X + k_Y)/2 is the
number of non-identity letters. Under dephasing along P at the same rate γ_P on every site a
string decays at 2·Σ_P γ_P·k_P (the absorption theorem's cost law on this figure).

The Hamiltonian moves strings (§8). A placed Pauli term T acting on σ gives 0 when they commute
and 2·(phase)·(T·σ) when they anticommute, one string; for a two-site term the weight w changes
by exactly one, and a bond P⊗P steps ±2 along one of the two axes other than P. A one-site field
keeps w and steps on a face diagonal.

Public API:
  cube_coords(label)                  -> (k_Z, k_X, k_Y)
  cube_coords_table(N)                -> (4^N, 3) array, row k = the point of basis index k
  cube_weight(label)                  -> w
  cube_rate(label, gamma_z, gamma_x=0, gamma_y=0) -> 2 Σ γ_P k_P
  cube_rate_per_site(label, rates)    -> the same cost with a rate triple per site
  cube_points(N)                      -> {(k_Z, k_X, k_Y): number of strings}
  cube_step(term, label)              -> (target_label, (dk_Z, dk_X, dk_Y)) or None
  commutator_coefficient(term, label) -> (target_label, c) with [term, string] = c·target, or None
  chain_terms(H)                      -> placed term labels of a PauliHamiltonian on an open chain
  cube_moves(terms, N)                -> set of steps (dk_Z, dk_X, dk_Y) the terms make on N sites
  cube_edges(terms, N)                -> set of pairs of lattice points the terms join
  cube_centroid(coeffs, N)            -> centroid of a Pauli-coefficient vector, weights |c|²

Labels are strings over 'IXYZ', site 0 leftmost, as fw._pauli_label writes them; coefficient
vectors are in fw.pauli_basis_vector's order. k_Z is the XY-weight, the total bit_a of
fw.total_bit_a; fw.PauliTerm.n_x etc. count letters of a Hamiltonian template, not of a string.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import product as iproduct
from math import factorial

import numpy as np

from .pauli import _pauli_label

AXES = ('Z', 'X', 'Y')                       # coordinate order (k_Z, k_X, k_Y)

# letter product without phase, and whether two letters anticommute
_OTHER = {frozenset('XY'): 'Z', frozenset('YZ'): 'X', frozenset('XZ'): 'Y'}


def _mul(a, b):
    if a == 'I':
        return b
    if b == 'I':
        return a
    if a == b:
        return 'I'
    return _OTHER[frozenset((a, b))]


def _anti(a, b):
    return a != 'I' and b != 'I' and a != b


# phase of a letter product: a·b = _PHASE[(a, b)] · _mul(a, b)
_PHASE = {(a, b): 1 for a in 'IXYZ' for b in 'IXYZ'}
_PHASE.update({('X', 'Y'): 1j, ('Y', 'Z'): 1j, ('Z', 'X'): 1j,
               ('Y', 'X'): -1j, ('Z', 'Y'): -1j, ('X', 'Z'): -1j})


def _check(label):
    if not label or any(c not in 'IXYZ' for c in label):
        raise ValueError(f"a Pauli label is a nonempty string over 'IXYZ'; got {label!r}")


def cube_coords(label):
    """(k_Z, k_X, k_Y): how many letters anticommute with Z, with X, with Y."""
    _check(label)
    nx, ny, nz = label.count('X'), label.count('Y'), label.count('Z')
    return (nx + ny, ny + nz, nx + nz)


def cube_weight(label):
    """w = n_X + n_Y + n_Z = (k_Z + k_X + k_Y)/2, the number of non-identity letters."""
    _check(label)
    return len(label) - label.count('I')


def cube_rate(label, gamma_z, gamma_x=0.0, gamma_y=0.0):
    """Decay rate of the string under dephasing along Z, X, Y at these rates on every site:
    2 (γ_Z k_Z + γ_X k_X + γ_Y k_Y). With site-dependent rates the cost is a sum over sites
    and not a function of the lattice point; this function does not cover that case."""
    kz, kx, ky = cube_coords(label)
    return 2 * (gamma_z * kz + gamma_x * kx + gamma_y * ky)


def cube_rate_per_site(label, rates):
    """Decay rate of the string when site l dephases along Z, X, Y at rates[l] = (γ_Z, γ_X, γ_Y):
    2 Σ_l Σ_P γ_P^l [letter l anticommutes with P]. This is the cost when the rates differ between
    sites; it is then not a function of the lattice point alone."""
    _check(label)
    if len(rates) != len(label):
        raise ValueError(f"one rate triple per site; got {len(rates)} for {len(label)} sites")
    return 2 * sum(g * _anti(c, P) for c, trip in zip(label, rates) for P, g in zip(AXES, trip))


@lru_cache(maxsize=None)
def _coords_table(N):
    K = np.array([cube_coords(_pauli_label(k, N)) for k in range(4 ** N)], dtype=int)
    K.flags.writeable = False
    return K


def cube_coords_table(N):
    """(4^N, 3) integer array: row k is the lattice point of the basis string with flat index k
    (fw._pauli_label order). A read-only cached array; copy it before modifying."""
    return _coords_table(N)


def cube_points(N):
    """The occupied lattice points and how many strings sit at each:
    N! / (n_I! n_X! n_Y! n_Z!) for the composition the point determines."""
    out = {}
    for nx, ny, nz in iproduct(range(N + 1), repeat=3):
        ni = N - nx - ny - nz
        if ni < 0:
            continue
        cnt = factorial(N) // (factorial(ni) * factorial(nx) * factorial(ny) * factorial(nz))
        out[(nx + ny, ny + nz, nx + nz)] = cnt
    return out


def cube_step(term, label):
    """The move a placed term (a label of the same length, e.g. 'XXI') makes on a string.

    Returns None when the two commute (the commutator is zero), otherwise the string the
    commutator is proportional to and the step (dk_Z, dk_X, dk_Y). The commutator is then
    2·(phase)·(term·string)."""
    _check(term)
    _check(label)
    if len(term) != len(label):
        raise ValueError(f"term and string must have the same length; got {term!r}, {label!r}")
    n_anti = sum(_anti(t, s) for t, s in zip(term, label))
    if n_anti % 2 == 0:
        return None
    target = ''.join(_mul(t, s) for t, s in zip(term, label))
    k0, k1 = cube_coords(label), cube_coords(target)
    return target, tuple(b - a for a, b in zip(k0, k1))


def commutator_coefficient(term, label):
    """[term, string] as c·target with the phase: None when they commute, otherwise the target
    string and c = 2·(phase of term·string), which is ±2i for two Hermitian strings that
    anticommute. Enough to rebuild the entries of ad_H in the string basis."""
    r = cube_step(term, label)
    if r is None:
        return None
    ph = 1
    for t, s in zip(term, label):
        ph *= _PHASE[(t, s)]
    return r[0], 2 * ph


def chain_terms(H):
    """Place the templates of a framework PauliHamiltonian on every window of an open chain,
    as the framework's k-body builders do: a k-letter template sits on sites i..i+k-1. Unlike those
    builders it also accepts one-letter templates (fields, placed on every site), and it ignores the
    coefficients, so a term of coefficient zero still contributes its moves."""
    N = H.chain_length
    out = []
    for t in H.terms:
        k = len(t.letters)
        for i in range(N - k + 1):
            out.append('I' * i + ''.join(t.letters) + 'I' * (N - k - i))
    return out


def _all_labels(N):
    return [''.join(t) for t in iproduct('IXYZ', repeat=N)]


def cube_moves(terms, N):
    """Every step (dk_Z, dk_X, dk_Y) the placed terms make on any string of N sites."""
    steps = set()
    for s in _all_labels(N):
        for t in terms:
            r = cube_step(t, s)
            if r is not None:
                steps.add(r[1])
    return steps


def cube_edges(terms, N):
    """The pairs of distinct lattice points the placed terms join (a move that lands on the
    point it started from, possible for terms on three or more sites, gives no edge)."""
    edges = set()
    for s in _all_labels(N):
        for t in terms:
            r = cube_step(t, s)
            if r is not None:
                a, b = cube_coords(s), cube_coords(r[0])
                if a != b:
                    edges.add(tuple(sorted((a, b))))
    return edges


def cube_centroid(coeffs, N):
    """Centroid ⟨(k_Z, k_X, k_Y)⟩ of a 4^N Pauli-coefficient vector (fw.pauli_basis_vector order),
    weighted by |c|². Along a trajectory under Z-dephasing at rate γ on every site,
    d/dt ln‖ρ‖² = −4γ·⟨k_Z⟩ (THE_ONE_SQUARE §8)."""
    c = np.asarray(coeffs)
    if c.shape != (4 ** N,):
        raise ValueError(f"expected a vector of length 4^N = {4 ** N}; got shape {c.shape}")
    p = np.abs(c) ** 2
    tot = p.sum()
    if tot == 0:
        raise ValueError("the zero vector has no centroid")
    return tuple(float(x) for x in (p @ _coords_table(N)) / tot)
