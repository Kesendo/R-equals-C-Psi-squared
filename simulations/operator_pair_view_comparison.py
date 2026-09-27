"""Signal ranks in the single-excitation block, and the owned objects that predict them.

The object is the (1,1) joint-popcount block of the open Heisenberg chain, the
N x N Haken-Strobl density block, with J = 1 and integer local Z rates.  A
preparation |a><a| and a population readout |b><b| give one scalar curve
y(t) = tr(|b><b| e^(tL) |a><a|); the rank of its moment-Hankel matrix is the
minimal order of the linear recurrence the curve obeys: the number of poles of
its Laplace transform, each counted with its order.

Every rank here is computed EXACTLY over Q (sympy DomainMatrix over ZZ on the
integer Hankel matrix), so no modular bound is involved.  A mode is hidden from a
curve in three ways, and the predictions below are upper bounds built from them:
the preparation cannot reach it, it has no population content on either side,
or it shares its eigenvalue with a mode the curve already sees (the whole kernel
reads as one constant).

  one seat j (F157):  rank(a -> b) <= min(rho_a, rho_b),
      rho_a = (N - s_a)^2 - (blind(j) - s_a),  blind(j) = (gcd(2j+1, N) - 1)/2,
      s_a = the modes with a node at both a and j.  A preparation at a lives in
      End(K_a), K_a the H-invariant span of e_a and e_j, of dimension N - s_a;
      inside it the kernel is the identity on K_a plus the projectors on the
      j-blind modes that K_a contains, and a single curve sees all of it as one
      constant.  On the seat itself s = blind(j), giving (N - blind)^2.
  uniform light (Omega_2, PROOF_UNIFORM_LAW B1):  the operators that commute
      with H and have no diagonal, dim floor(N/2), are the whole -4 gamma root
      space (F140's frozen modes at its uniform point), shared by L and L-dagger
      and orthogonal to every population preparation and readout.  Bound: the
      symmetry-allowed dimension minus floor(N/2).
  the R90 locus (F140):  on anti-palindromic profiles the frozen divisor holds
      floor(N/2) modes at -4 gamma-bar with zero weight on the diagonal cells, on
      both sides since L-dagger = L(-J) lies on the locus too.  Bound:
      N^2 - (dim ker L - 1) - floor(N/2).
  any other profile:  N^2 - (dim ker L - 1) - (operators of Omega_2 whose cells
      all sit at one nonzero price), an upper bound only.

Under uniform light the channels give more than a bound.  In F126's renewal the
kernel is |G(t)|^2, and on this chain it is diagonal in B0's Neumann modes
psi_q: the pairs (k, l) and (N-l, N-k) have the same gap, and the cross terms of
their products cancel.  So y(a -> b) = sum_q psi_q(a) psi_q(b) Y_q(t), channel q
has N + 1 poles when q = N mod 2, N otherwise, and one for q = 0, and the rank
is at most the sum over the modes with no node at a or at b.  F157's node
condition puts every blind mode in the N + 1 class.  Checked exactly: every
population moment matrix commutes with H_SE (on the XY chain it does not), and
the channel sum meets every uniform map and reachable dimension computed.

Also: the N = 3 parity motion (the reflection-odd weight after a centre
preparation).

Usage: python simulations/operator_pair_view_comparison.py [--output-dir DIR]
"""

from __future__ import annotations

import argparse
import json
from math import gcd
from numbers import Integral
from pathlib import Path
from typing import Sequence

import numpy as np
import sympy as sp
from scipy.linalg import expm
from sympy import ZZ
from sympy.polys.matrices import DomainMatrix


# Profiles shown on the page, with the label the figures use.
CASES = {
    (3, (0, 1, 0)): "centre only",
    (3, (1, 1, 1)): "uniform",
    (3, (1, 2, 1)): "symmetric, unequal",
    (3, (1, 1, 0)): "one-sided",
    (3, (0, 3, 0)): "centre only, same total as uniform",
    (4, (0, 1, 1, 0)): "symmetric, no fixed seat",
    (4, (1, 1, 1, 0)): "one-sided",
    (4, (1, 1, 1, 1)): "uniform",
    (5, (0, 0, 1, 0, 0)): "centre only",
    (5, (1, 0, 0, 0, 0)): "one end seat",
    (5, (0, 1, 0, 0, 0)): "one off-centre seat",
    (5, (1, 1, 1, 1, 1)): "uniform",
    (5, (1, 0, 1, 0, 0)): "one-sided",
    (5, (0, 1, 0, 1, 0)): "symmetric, two seats",
    (4, (0, 2, 0, 0)): "one middle seat at 2J",
    (6, (1, 0, 0, 0, 0, 0)): "one end seat",
    (6, (0, 1, 0, 0, 0, 0)): "one off-centre seat",
    (6, (0, 0, 2, 0, 0, 0)): "one middle seat at 2J",
    (6, (1, 1, 1, 1, 1, 1)): "uniform",
    (7, (0, 0, 0, 1, 0, 0, 0)): "centre only",
    (7, (1, 1, 1, 1, 1, 1, 1)): "uniform",
    (9, (0, 0, 0, 0, 1, 0, 0, 0, 0)): "centre only",
    (9, (1, 1, 1, 1, 1, 1, 1, 1, 1)): "uniform",
    (3, (0, 1, 2)): "R90 locus",
    (4, (0, 1, 2, 3)): "R90 locus",
    (5, (0, 1, 2, 3, 4)): "R90 locus",
    (4, (0, 2, 0, 2)): "R90 locus",
    (5, (1, 0, 1, 2, 1)): "R90 locus",
}


def _integer_profile(n: int, gamma_per_site: Sequence[int]) -> tuple[int, ...]:
    if n < 2 or len(gamma_per_site) != n:
        raise ValueError("n >= 2 and one site rate per site required")
    if any(not isinstance(rate, Integral) or rate < 0 for rate in gamma_per_site):
        raise ValueError("nonnegative integer site rates required")
    return tuple(int(rate) for rate in gamma_per_site)


def _single_excitation_h(n: int, book: str = "heisenberg", coupling: int = 1) -> np.ndarray:
    """Integer H_SE on an open chain: J sum(XX+YY+ZZ) gives J((N-1) Id - 2 Laplacian),
    J sum(XX+YY) gives 2J x adjacency, without the ZZ term's diagonal."""
    if book not in ("heisenberg", "xy"):
        raise ValueError("book must be 'heisenberg' or 'xy'")
    if not isinstance(coupling, Integral) or coupling < 1:
        raise ValueError("a positive integer coupling J is required")
    h = np.zeros((n, n), dtype=np.int64)
    if book == "heisenberg":
        for site in range(n):
            degree = int(site > 0) + int(site < n - 1)
            h[site, site] = (n - 1) - 2 * degree
    for site in range(n - 1):
        h[site, site + 1] = h[site + 1, site] = 2
    return h * int(coupling)


def integer_hermitian_generator(n: int, gamma_per_site: Sequence[int],
                                book: str = "heisenberg", coupling: int = 1) -> np.ndarray:
    """Exact integer L in coordinates (populations, Re upper, Im upper)."""
    gamma = _integer_profile(n, gamma_per_site)
    h = _single_excitation_h(n, book, coupling)
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    labels = [("population", a, a) for a in range(n)]
    labels += [("real", a, b) for a, b in pairs]
    labels += [("imag", a, b) for a, b in pairs]
    d = n * n
    generator = np.zeros((d, d), dtype=np.int64)

    for source, (kind, a, b) in enumerate(labels):
        real = np.zeros((n, n), dtype=np.int64)
        imag = np.zeros((n, n), dtype=np.int64)
        if kind == "population":
            real[a, a] = 1
        elif kind == "real":
            real[a, b] = real[b, a] = 1
        else:
            imag[a, b], imag[b, a] = 1, -1

        # -i[H, real + i*imag] = [H, imag] - i[H, real].
        image_real = h @ imag - imag @ h
        image_imag = -(h @ real - real @ h)
        for left in range(n):
            for right in range(n):
                if left != right:
                    cost = 2 * (gamma[left] + gamma[right])
                    image_real[left, right] -= cost * real[left, right]
                    image_imag[left, right] -= cost * imag[left, right]
        for target, (target_kind, left, right) in enumerate(labels):
            generator[target, source] = (
                image_real[left, right] if target_kind != "imag" else image_imag[left, right]
            )
    return generator


def _sparse_rows(generator: np.ndarray) -> list[list[tuple[int, int]]]:
    return [[(j, int(x)) for j, x in enumerate(row) if x] for row in generator]


def _apply(rows: list[list[tuple[int, int]]], vector: list[int]) -> list[int]:
    return [sum(x * vector[j] for j, x in row) for row in rows]


def exact_hankel_rank(generator: np.ndarray, preparation: int, readout: int) -> int:
    """Rank over Q of the d x d moment-Hankel matrix of y(t) = <readout| e^(tL) |preparation>.

    `preparation` and `readout` index the population coordinates.  The moments
    m_k = tr(O L^k rho_0) are exact integers; the rank is taken by DomainMatrix
    over ZZ, fraction-free, so nothing is rounded and no prime is involved."""
    d = generator.shape[0]
    rows = _sparse_rows(generator)
    state = [0] * d
    state[preparation] = 1
    moments = []
    for _ in range(2 * d - 1):
        moments.append(state[readout])
        state = _apply(rows, state)
    hankel = DomainMatrix([[ZZ(moments[i + j]) for j in range(d)] for i in range(d)], (d, d), ZZ)
    return int(hankel.rank())


def exact_krylov_dimension(generator: np.ndarray, preparation: int) -> int:
    """dim span{L^k rho_0} over Q for rho_0 = |a><a|: every mode the preparation reaches."""
    d = generator.shape[0]
    rows = _sparse_rows(generator)
    state = [0] * d
    state[preparation] = 1
    columns = []
    for _ in range(d):
        columns.append(state)
        state = _apply(rows, state)
    krylov = DomainMatrix([[ZZ(columns[j][i]) for j in range(d)] for i in range(d)], (d, d), ZZ)
    return int(krylov.rank())


def exact_rank_map(n: int, gamma_per_site: Sequence[int],
                   coupling: int = 1) -> tuple[tuple[int, ...], ...]:
    """Exact rank of every (preparation a, readout b) population signal."""
    generator = integer_hermitian_generator(n, gamma_per_site, coupling=coupling)
    return tuple(tuple(exact_hankel_rank(generator, a, b) for b in range(n)) for a in range(n))


# ----------------------------------------------------------------------
# the two laws
# ----------------------------------------------------------------------

def f157_blind(n: int, seat: int) -> int:
    """Uniform Heisenberg chain, ZZ on: blind(j) = (gcd(2j+1, N) - 1)/2 (F157)."""
    return (gcd(2 * seat + 1, n) - 1) // 2


def node_modes(n: int, site: int) -> set[int]:
    """Neumann modes of the uniform Heisenberg chain with an exact node at `site`:
    psi_m(j) ~ cos(pi m (2j+1)/(2N)) vanishes iff (2j+1) m = N mod 2N."""
    return {m for m in range(n) if ((2 * site + 1) * m - n) % (2 * n) == 0}


def single_seat_prediction(n: int, seat: int) -> tuple[tuple[int, ...], ...]:
    """Light on one seat j.  A preparation at site a can only reach End(K_a), with
    K_a the H-invariant span of e_a and e_j (the jump at j needs e_j), whose
    dimension is N - s_a, s_a the modes with a node at BOTH a and j (The Blind
    Site's intersection law).  Inside End(K_a) the kernel is the identity on K_a
    plus the projectors on the blind(j) - s_a blind modes of j that K_a contains,
    and one curve sees that whole kernel as a single constant.  A readout at b is
    bounded the same way by duality."""
    blind = f157_blind(n, seat)
    shared = [len(node_modes(n, a) & node_modes(n, seat)) for a in range(n)]
    rho = [(n - s) ** 2 - (blind - s) for s in shared]
    return tuple(tuple(min(rho[a], rho[b]) for b in range(n)) for a in range(n))


def on_r90_locus(gamma_per_site: Sequence[int]) -> bool:
    """Anti-palindromic profile: gamma_l + gamma_(N-1-l) the same for every l."""
    g = tuple(gamma_per_site)
    return all(g[k] + g[-1 - k] == g[0] + g[-1] for k in range(len(g)))


def _kernel_dimension(n: int, gamma: Sequence[int]) -> int:
    generator = integer_hermitian_generator(n, gamma)
    return n * n - int(DomainMatrix([[ZZ(int(x)) for x in row] for row in generator],
                                    (n * n, n * n), ZZ).rank())


def r90_locus_prediction(n: int, gamma_per_site: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    """On the R90 locus F140's frozen divisor holds floor(N/2) modes at -4 gamma-bar
    with zero weight on the diagonal cells (PROOF_R90_FROZEN_DIVISOR section 3), on
    both sides since L-dagger = L(-J) is on the locus as well: no population curve
    sees them.  For a profile that is not reflection-symmetric."""
    gamma = _integer_profile(n, gamma_per_site)
    if not on_r90_locus(gamma) or gamma == tuple(reversed(gamma)):
        raise ValueError("an anti-palindromic profile that is not reflection-symmetric is required")
    bound = n * n - (_kernel_dimension(n, gamma) - 1) - n // 2
    return tuple((bound,) * n for _ in range(n))


def omega2_dimension(n: int) -> int:
    """Exact dim over Q of {X Hermitian : [H_SE, X] = 0, diag X = 0}."""
    h = sp.Matrix(_single_excitation_h(n).tolist())
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    total = 0
    # A Hermitian X splits into a real symmetric and a real antisymmetric part,
    # each commuting with the real symmetric H on its own.
    for sign in (+1, -1):
        unknowns = sp.symbols(f"x0:{len(pairs)}")
        x = sp.zeros(n)
        for k, (a, b) in enumerate(pairs):
            x[a, b] = unknowns[k]
            x[b, a] = sign * unknowns[k]
        commutator = h * x - x * h
        equations = [commutator[a, b] for a in range(n) for b in range(n)]
        system = sp.Matrix([[sp.diff(eq, u) for u in unknowns] for eq in equations])
        total += len(pairs) - system.rank()
    return total


def single_rate_omega2_dimension(n: int, gamma_per_site: Sequence[int]) -> int:
    """dim of the operators that commute with H_SE, have no diagonal, and sit on
    cells of ONE nonzero dephasing cost 2(gamma_a + gamma_b).  Each is a shared
    L / L-dagger mode (D multiplies it by -cost, [H, .] kills it) with no
    population content, invisible to every population signal.  Cost 0 is left
    out: those operators are stationary and belong to the kernel term."""
    gamma = _integer_profile(n, gamma_per_site)
    h = sp.Matrix(_single_excitation_h(n).tolist())
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    total = 0
    for cost in sorted({2 * (gamma[a] + gamma[b]) for a, b in pairs} - {0}):
        cells = [(a, b) for a, b in pairs if 2 * (gamma[a] + gamma[b]) == cost]
        for sign in (+1, -1):
            unknowns = sp.symbols(f"y0:{len(cells)}")
            x = sp.zeros(n)
            for k, (a, b) in enumerate(cells):
                x[a, b] = unknowns[k]
                x[b, a] = sign * unknowns[k]
            commutator = h * x - x * h
            equations = [commutator[a, b] for a in range(n) for b in range(n)]
            system = sp.Matrix([[sp.diff(eq, u) for u in unknowns] for eq in equations])
            total += len(cells) - system.rank()
    return total


def general_bound(n: int, gamma_per_site: Sequence[int]) -> tuple[tuple[int | None, ...], ...]:
    """An UPPER bound for any profile, from the three hiding mechanisms: N^2 minus
    the conserved directions beyond the one every curve keeps (a single curve sees
    the whole kernel as one constant) minus the operators of Omega_2 whose cells
    all sit at one nonzero price (each a shared L / L-dagger mode with no
    population content); and on the centre's row and column of an odd chain with
    a reflection-symmetric profile, where a reflection-invariant preparation or
    readout cannot reach the reflection-odd operators, the reflection-invariant
    dimension ceil(N/2)^2 + floor(N/2)^2 in place of N^2.  That branch is used only
    where the kernel is one-dimensional, and entries it does not cover are None.
    It is not tight in general: on the R90 locus F140's frozen modes hide more, and
    some modes hide from one site and not another."""
    gamma = _integer_profile(n, gamma_per_site)
    kernel = _kernel_dimension(n, gamma)
    invisible = single_rate_omega2_dimension(n, gamma)
    symmetric = n % 2 == 1 and gamma == tuple(reversed(gamma))
    reflection_dim = ((n + 1) // 2) ** 2 + (n // 2) ** 2
    centre = n // 2

    def entry(a: int, b: int) -> int | None:
        if symmetric and centre in (a, b):
            return reflection_dim - invisible if kernel == 1 else None
        return n * n - (kernel - 1) - invisible

    return tuple(tuple(entry(a, b) for b in range(n)) for a in range(n))


def uniform_prediction(n: int) -> tuple[tuple[int, ...], ...]:
    """N^2 - dim Omega_2 off the centre; at the centre of odd N, where preparation
    or readout is reflection-invariant, (even^2 + odd^2) - dim Omega_2."""
    o2 = omega2_dimension(n)
    centre = n // 2 if n % 2 else None
    reflection_dim = ((n + 1) // 2) ** 2 + (n // 2) ** 2
    return tuple(tuple(reflection_dim - o2 if centre in (a, b) else n * n - o2
                       for b in range(n)) for a in range(n))


# ----------------------------------------------------------------------
# uniform light, channel by channel
# ----------------------------------------------------------------------

def channel_order(n: int, q: int) -> int:
    """Poles of channel q of the uniformly lit Heisenberg chain.  Its renewal kernel
    psi_q^T |G(t)|^2 psi_q carries the gaps E_k - E_l of the pairs whose product
    psi_k psi_l has a component on psi_q, all of the form
    8 sin(pi q/2N) sin(pi x/2N) with x = q mod 2 and |x| <= N, each with positive
    weight; the renewal keeps the pole count.  So: 1 for q = 0, N + 1 when
    q = N mod 2, N otherwise."""
    if not 0 <= q < n:
        raise ValueError("0 <= q < N required")
    return 1 if q == 0 else n + int((q - n) % 2 == 0)


def channel_prediction(n: int) -> tuple[tuple[int, ...], ...]:
    """Uniform light: y(a -> b) = sum_q psi_q(a) psi_q(b) Y_q(t), so its rank is at
    most the sum of the channel orders over the modes with no node at a or at b,
    and equal to it wherever no two of those channels share a pole."""
    nodes = [node_modes(n, site) for site in range(n)]
    return tuple(tuple(sum(channel_order(n, q) for q in range(n) if q not in nodes[a] | nodes[b])
                       for b in range(n)) for a in range(n))


def population_moments_commute(n: int, gamma: int, book: str = "heisenberg") -> int | None:
    """Uniform integer rate gamma.  None when every population moment matrix
    M_k[b][a] = (L^k |a><a|)_bb for k < N^2 commutes with H_SE, which by
    Cayley-Hamilton is every k, so the population transfer is diagonal in the
    eigenbasis of H_SE; otherwise the first k at which it does not."""
    generator = integer_hermitian_generator(n, (gamma,) * n, book)
    h = [[int(x) for x in row] for row in _single_excitation_h(n, book)]
    rows = _sparse_rows(generator)
    states = [[int(i == a) for i in range(n * n)] for a in range(n)]
    for k in range(n * n):
        m = [[states[a][b] for a in range(n)] for b in range(n)]
        mh = [[sum(m[i][t] * h[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
        hm = [[sum(h[i][t] * m[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
        if mh != hm:
            return k
        states = [_apply(rows, state) for state in states]
    return None


def frozen_root_kernel_dims(n: int, gamma_per_site: Sequence[int], coupling: int = 1,
                            powers: int = 4) -> tuple[int, ...]:
    """dim ker (L + 4 gamma-bar)^k for k = 1..powers, exactly: the integer matrix
    N L + 4 (sum gamma) Id has the same kernels.  The sequence stops growing at the
    algebraic multiplicity of the frozen root, its first entry is the geometric one,
    and its increments give the Jordan block sizes."""
    gamma = _integer_profile(n, gamma_per_site)
    generator = integer_hermitian_generator(n, gamma, coupling=coupling)
    d = n * n
    shift = 4 * sum(gamma)
    rows = [[ZZ(n * int(generator[i, j]) + (shift if i == j else 0)) for j in range(d)]
            for i in range(d)]
    m = DomainMatrix(rows, (d, d), ZZ)
    power, dims = m, []
    for _ in range(powers):
        dims.append(d - int(power.rank()))
        power = power * m
    return tuple(dims)


def coupling_scan(couplings: Sequence[int] = (1, 2, 3)) -> list[dict]:
    """The maps that fall below their bound at J = 1, and one that meets it, read
    at several integer couplings: the frozen root's kernel sequence on the R90
    locus, the number of distinct eigenvalues of L, and the values of the map."""
    s = sp.symbols("s")
    out = []
    for n, gamma in ((4, (0, 2, 0, 2)), (5, (1, 0, 1, 2, 1)), (4, (0, 1, 2, 3)),
                     (5, (0, 1, 0, 1, 0))):
        for coupling in couplings:
            generator = integer_hermitian_generator(n, gamma, coupling=coupling)
            d = n * n
            chi = DomainMatrix([[ZZ(int(x)) for x in row] for row in generator], (d, d), ZZ).charpoly()
            distinct = sp.Poly(sp.sqf_part(sp.Poly(chi, s).as_expr()), s).degree()
            ranks = exact_rank_map(n, gamma, coupling=coupling)
            out.append({"n": n, "gamma": list(gamma), "coupling": int(coupling),
                        "frozen_root_kernel_dims": (list(frozen_root_kernel_dims(n, gamma, coupling))
                                                    if on_r90_locus(gamma) else None),
                        "distinct_eigenvalues": int(distinct),
                        "ranks": [list(r) for r in ranks]})
    return out


def middle_seat_rate_scan(n: int, rates: Sequence[int] = (1, 2, 3, 4)) -> list[dict]:
    """Light on the middle seat j = N/2 - 1 of an even chain at several integer
    rates: the reachable dimension from |j><j|, the rank of the return, the
    one-seat bound on that return, and the number of distinct eigenvalues of L
    (the square-free part of its exact characteristic polynomial).  A return below
    the bound with every eigenvalue simple means residues vanish, not that two
    eigenvalues merge."""
    if n % 2 or n < 4:
        raise ValueError("an even N >= 4 is required")
    seat = n // 2 - 1
    s = sp.symbols("s")
    out = []
    for rate in rates:
        gamma = tuple(int(rate) if k == seat else 0 for k in range(n))
        generator = integer_hermitian_generator(n, gamma)
        d = n * n
        chi = DomainMatrix([[ZZ(int(x)) for x in row] for row in generator], (d, d), ZZ).charpoly()
        distinct = sp.Poly(sp.sqf_part(sp.Poly(chi, s).as_expr()), s).degree()
        bound = single_seat_prediction(n, seat)
        out.append({"n": n, "seat": seat, "rate": int(rate),
                    "reachable_dimension": exact_krylov_dimension(generator, seat),
                    "return_rank": exact_hankel_rank(generator, seat, seat),
                    "one_seat_bound": bound[seat][seat],
                    "map_meets_bound": (exact_rank_map(n, gamma) == bound) if n <= 6 else None,
                    "distinct_eigenvalues": int(distinct), "dimension": d})
    return out


def blind_seat_reach(n: int) -> list[dict]:
    """Uniform light: the reachable dimension from |j><j| at every seat up to the
    centre (the rest by reflection), beside the uniform bound, F157's blindness at
    that seat, and the channel sum over the modes with no node at j.  Every blind
    mode q satisfies q = N mod 2 and costs N + 1, so that sum is
    N^2 - floor(N/2) - blind(j)(N + 1); at the centre of an odd chain it is the
    reflection bound, since N^2 - ceil(N/2)^2 - floor(N/2)^2 = m(N + 1)."""
    generator = integer_hermitian_generator(n, (1,) * n)
    bound = uniform_prediction(n)
    channels = channel_prediction(n)
    out = []
    for seat in range((n + 1) // 2):
        reach = exact_krylov_dimension(generator, seat)
        closed_form = n * n - n // 2 - f157_blind(n, seat) * (n + 1)
        if channels[seat][seat] != closed_form:
            raise AssertionError("the channel sum and its closed form disagree")
        out.append({"seat": seat, "f157_blind": f157_blind(n, seat),
                    "uniform_bound": bound[seat][seat], "reachable_dimension": reach,
                    "below_bound_by": bound[seat][seat] - reach,
                    "channel_sum": channels[seat][seat],
                    "meets_channel_sum": reach == channels[seat][seat]})
    return out


# ----------------------------------------------------------------------
# N = 3 parity motion (the second reading of the same flow)
# ----------------------------------------------------------------------

def _twice_odd_readout() -> np.ndarray:
    """2*Tr(P_odd rho) for N=3 in the integer Hermitian coordinates."""
    readout = np.zeros(9, dtype=np.int64)
    readout[0] = readout[2] = 1
    readout[4] = -2  # Re rho_02, after the populations and the (0,1) coordinate.
    return readout


def odd_weight_derivatives(gamma_per_site: Sequence[int], count: int) -> tuple[sp.Rational, ...]:
    """Exact derivatives at t=0 of Tr(P_odd rho(t)) from a centre preparation, N=3."""
    gamma = _integer_profile(3, gamma_per_site)
    if count < 1:
        raise ValueError("positive derivative count required")
    generator = sp.Matrix(integer_hermitian_generator(3, gamma).tolist())
    state = sp.zeros(9, 1)
    state[1] = 1
    readout = sp.Matrix([int(x) for x in _twice_odd_readout()]).T
    values = []
    for _ in range(count):
        values.append(sp.Rational((readout * state)[0], 2))
        state = generator * state
    return tuple(values)


def odd_weight_curve(gamma_per_site: Sequence[int], times: Sequence[float]) -> tuple[float, ...]:
    """Hilbert reflection-odd population after preparing site 1 at N=3."""
    generator = integer_hermitian_generator(3, gamma_per_site).astype(float)
    initial = np.zeros(9)
    initial[1] = 1
    readout = _twice_odd_readout().astype(float) / 2
    return tuple(float(readout @ (expm(float(t) * generator) @ initial)) for t in times)


def centre_return_signature(gamma_per_site: Sequence[int]) -> dict[str, float | int]:
    """N=3 centre population: an exact third derivative and one finite-time read."""
    generator = integer_hermitian_generator(3, gamma_per_site)
    state = np.zeros(9, dtype=np.int64)
    state[1] = 1
    initial = state.copy()
    for _ in range(3):
        state = generator @ state
    return {
        "third_derivative_at_zero": int(state[1]),
        "population_at_t_0_5": float((expm(0.5 * generator.astype(float)) @ initial)[1]),
    }


def cramer_return_rank(n: int, gamma_per_site: Sequence[int], seat: int) -> int:
    """The return's rank by the Cramer identity, F157's count one level up: the
    (j, j) entry of (s - L)^(-1) is the principal cofactor over det(s - L) for any
    matrix, so the rank is N^2 - deg gcd(chi(L), chi(L with that population row and
    column struck)).  Integer characteristic polynomials, exact."""
    generator = integer_hermitian_generator(n, gamma_per_site)
    d = n * n
    s = sp.symbols("s")
    full = [[int(x) for x in row] for row in generator]
    keep = [k for k in range(d) if k != seat]
    cut = [[full[r][c] for c in keep] for r in keep]
    chi = sp.Poly(DomainMatrix([[ZZ(x) for x in row] for row in full], (d, d), ZZ).charpoly(), s)
    chi_cut = sp.Poly(DomainMatrix([[ZZ(x) for x in row] for row in cut], (d - 1, d - 1), ZZ).charpoly(), s)
    return d - sp.gcd(chi, chi_cut).degree()


def centre_return_recurrence(gamma_per_site: Sequence[int]) -> tuple[int, ...]:
    """Monic minimal recurrence of the centre return y(t) = <c| rho(t) |c> from
    rho(0) = |c><c|, read from the exact moments: coefficients from the highest
    power of s down.
    The recurrence is solved on the leading Hankel block of the curve's rank and
    checked on every moment the generator's dimension provides."""
    gamma = _integer_profile(len(gamma_per_site), gamma_per_site)
    n = len(gamma)
    centre = n // 2
    generator = integer_hermitian_generator(n, gamma)
    rank = exact_hankel_rank(generator, centre, centre)
    rows = _sparse_rows(generator)
    d = n * n
    state = [int(i == centre) for i in range(d)]
    moments = []
    for _ in range(2 * d):
        moments.append(state[centre])
        state = _apply(rows, state)
    block = sp.Matrix(rank, rank, lambda i, j: moments[i + j])
    rhs = sp.Matrix(rank, 1, lambda i, _: -moments[i + rank])
    low = block.LUsolve(rhs)
    for i in range(len(moments) - rank):
        if moments[i + rank] + sum(low[k] * moments[i + k] for k in range(rank)) != 0:
            raise AssertionError("the recurrence fails on a later moment")
    if any(not c.is_integer for c in low):
        raise AssertionError("a non-integer recurrence coefficient")
    return (1,) + tuple(int(low[k]) for k in reversed(range(rank)))


# ----------------------------------------------------------------------
# figures
# ----------------------------------------------------------------------

def _draw_rank_map(ax, ranks, *, vmin: int, vmax: int, title: str):
    image = ax.imshow(ranks, cmap="viridis", vmin=vmin, vmax=vmax)
    n = len(ranks)
    ax.set(title=title, xlabel="readout site b", ylabel="preparation site a",
           xticks=range(n), yticks=range(n))
    threshold = (vmin + vmax) / 2
    for a in range(n):
        for b in range(n):
            value = ranks[a][b]
            ax.text(b, a, str(value), ha="center", va="center",
                    color="black" if value > threshold else "white",
                    fontsize=11, fontweight="bold")
    return image


def _render_n3(maps: dict, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    shown = [((0, 1, 0), "centre only (0,1,0)\nF157: (3−1)² = 4, 9−1 = 8"),
             ((1, 1, 1), "uniform (1,1,1)\nΩ₂: 5−1 = 4, 9−1 = 8"),
             ((1, 2, 1), "symmetric, unequal (1,2,1)\nthe upper bound, met: 5 and 9"),
             ((1, 1, 0), "one-sided (1,1,0)\nfull: 9")]
    palette = {(0, 1, 0): "#2864a6", (1, 1, 1): "#6b3c91", (1, 1, 0): "#c45c3c",
               (1, 2, 1): "#3f8f5a"}
    fig = plt.figure(figsize=(15.5, 8.4))
    grid = fig.add_gridspec(2, 5, width_ratios=(1, 1, 1, 1, 0.055),
                            height_ratios=(1, 1.0), hspace=0.45, wspace=0.34)
    fig.suptitle("N = 3: one rank, two laws, and the parity motion they do not show", fontsize=14)
    image = None
    for col, (gamma, title) in enumerate(shown):
        image = _draw_rank_map(fig.add_subplot(grid[0, col]), maps[(3, gamma)],
                               vmin=4, vmax=9, title=title)
    fig.colorbar(image, cax=fig.add_subplot(grid[0, 4]), label="exact signal rank")

    ax = fig.add_subplot(grid[1, :4])
    times = np.linspace(0, 1.4, 141)
    for gamma, _ in shown:
        ax.plot(times, odd_weight_curve(gamma, times), color=palette[gamma], lw=2.4,
                label=f"γ = {gamma}")
    ax.set(title="prepared at the centre: weight in the reflection-odd space",
           xlabel="time t", ylabel="Tr(P₋ ρ(t))", xlim=(0, 1.4), ylim=(-0.015, 0.52))
    ax.legend(loc="upper left", frameon=False)
    ax.grid(alpha=0.22)
    fig.subplots_adjust(left=0.06, right=0.95, top=0.88, bottom=0.08)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=170)
    plt.close(fig)


def _render_n5(maps: dict, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    shown = [((0, 0, 1, 0, 0), "centre only\nF157: 9 on the cross, 23 elsewhere"),
             ((1, 1, 1, 1, 1), "uniform\nΩ₂: 11 on the cross, 23 elsewhere"),
             ((1, 0, 1, 0, 0), "one-sided\nfull: 25")]
    fig = plt.figure(figsize=(15.0, 5.6))
    grid = fig.add_gridspec(1, 4, width_ratios=(1, 1, 1, 0.055), wspace=0.33)
    fig.suptitle("N = 5: the two laws that agreed at N = 3 part at the centre", fontsize=14)
    image = None
    for col, (gamma, title) in enumerate(shown):
        image = _draw_rank_map(fig.add_subplot(grid[0, col]), maps[(5, gamma)],
                               vmin=9, vmax=25, title=title)
    fig.colorbar(image, cax=fig.add_subplot(grid[0, 3]), label="exact signal rank")
    fig.subplots_adjust(left=0.06, right=0.94, top=0.80, bottom=0.12)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=170)
    plt.close(fig)


# ----------------------------------------------------------------------
# the run
# ----------------------------------------------------------------------

def run_view_comparison(output_dir: Path) -> dict:
    maps = {}
    for (n, gamma) in CASES:
        maps[(n, gamma)] = exact_rank_map(n, gamma)

    law_checks = []
    for (n, gamma), ranks in maps.items():
        support = [k for k, g in enumerate(gamma) if g]
        if len(support) == 1:
            law, predicted = f"F157, one seat {support[0]}", single_seat_prediction(n, support[0])
        elif len(set(gamma)) == 1:
            law, predicted = "Omega_2, PROOF_UNIFORM_LAW B1", uniform_prediction(n)
        elif on_r90_locus(gamma) and gamma != tuple(reversed(gamma)):
            law, predicted = "F140 frozen divisor, R90 locus", r90_locus_prediction(n, gamma)
        else:
            law, predicted = "upper bound only", general_bound(n, gamma)
        misses = sorted({(a, b) for a in range(n) for b in range(n)
                         if predicted[a][b] is not None and predicted[a][b] != ranks[a][b]})
        law_checks.append({"n": n, "gamma": list(gamma), "law": law,
                           "predicted": [list(r) for r in predicted],
                           "exact": [list(r) for r in ranks],
                           "entries_below_the_bound": [list(m) for m in misses
                                                       if ranks[m[0]][m[1]] < predicted[m[0]][m[1]]],
                           "entries_above_the_bound": [list(m) for m in misses
                                                       if ranks[m[0]][m[1]] > predicted[m[0]][m[1]]],
                           "met": not misses})

    omega2 = {n: omega2_dimension(n) for n in range(3, 10)}
    reach = {n: blind_seat_reach(n) for n in (6, 8, 9, 10, 12)}

    channel_checks = []
    for (n, gamma), ranks in maps.items():
        if len(set(gamma)) == 1:
            channel_checks.append({"n": n, "rate": gamma[0], "meets_channel_sum":
                                   ranks == channel_prediction(n)})
    for n in (3, 4, 5, 6):
        for rate in (2, 3):
            channel_checks.append({"n": n, "rate": rate, "meets_channel_sum":
                                   exact_rank_map(n, (rate,) * n) == channel_prediction(n)})
    commutation = {book: {n: population_moments_commute(n, 1, book) for n in range(3, 8)}
                   for book in ("heisenberg", "xy")}
    middle_seat = {n: middle_seat_rate_scan(n) for n in (4, 6, 8)}
    couplings = coupling_scan()
    recurrences = {str(gamma): list(centre_return_recurrence(gamma))
                   for gamma in ((0, 1, 0), (1, 1, 1))}

    output_dir.mkdir(parents=True, exist_ok=True)
    n3_figure = output_dir / "n3_rank_and_parity.png"
    n5_figure = output_dir / "n5_rank_maps.png"
    _render_n3(maps, n3_figure)
    _render_n5(maps, n5_figure)

    parity = {}
    for gamma in ((0, 1, 0), (1, 1, 1), (1, 2, 1), (1, 1, 0), (0, 3, 0)):
        parity[str(gamma)] = {
            "derivatives_at_zero": [int(v) for v in odd_weight_derivatives(gamma, 4)],
            "weight_at_t_0_5": odd_weight_curve(gamma, (0.5,))[0],
            **centre_return_signature(gamma),
        }
    summary = {
        "scope": "open Heisenberg chain, J = 1, integer local Z rates, the (1,1) block; "
                 "every rank exact over Q",
        "rank_maps": [{"n": n, "gamma": list(g), "label": CASES[(n, g)],
                       "ranks": [list(r) for r in maps[(n, g)]]} for (n, g) in CASES],
        "law_checks": law_checks,
        "omega2_dimension": omega2,
        "uniform_reach_by_seat": reach,
        "uniform_channel_checks": channel_checks,
        "first_noncommuting_population_moment": commutation,
        "middle_seat_rate_scan": middle_seat,
        "coupling_scan": couplings,
        "n3_centre_return_recurrence": recurrences,
        "n3_centre_preparation": parity,
        "figures": [n3_figure.name, n5_figure.name],
    }
    (output_dir / "operator_pair_view_comparison.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path,
                        default=Path(__file__).parent / "results" / "operator_pair_view_comparison")
    args = parser.parse_args()
    summary = run_view_comparison(args.output_dir)
    for check in summary["law_checks"]:
        if check["met"]:
            verdict = "met at every entry it predicts"
        else:
            below = sorted(a for a, b in check["entries_below_the_bound"] if a == b)
            verdict = (f"exact rank BELOW the bound at {len(check['entries_below_the_bound'])} entries"
                       f" (the returns at seats {below}); above it at "
                       f"{len(check['entries_above_the_bound'])}")
        print(f"N = {check['n']} gamma = {tuple(check['gamma'])}: {check['law']}: {verdict}")
    print(f"dim Omega_2, N = 3..9: {summary['omega2_dimension']}")
    for n, rows in summary["uniform_reach_by_seat"].items():
        print(f"uniform N = {n}: " + "; ".join(
            f"seat {r['seat']} (blind {r['f157_blind']}): reach {r['reachable_dimension']}"
            f" of bound {r['uniform_bound']}" for r in rows))
        print(f"    reach = channel sum = N^2 - floor(N/2) - blind(j)(N+1) at every seat: "
              f"{all(r['meets_channel_sum'] for r in rows)}")
    for check in summary["uniform_channel_checks"]:
        print(f"uniform N = {check['n']}, rate {check['rate']}: the map is the channel sum: "
              f"{check['meets_channel_sum']}")
    for book, firsts in summary["first_noncommuting_population_moment"].items():
        print(f"{book}: first population moment not commuting with H_SE, N = 3..7: "
              + ", ".join("none" if k is None else str(k) for k in firsts.values()))
    for n, rows in summary["middle_seat_rate_scan"].items():
        print(f"N = {n}, light on the middle seat {rows[0]['seat']}: " + "; ".join(
            f"rate {r['rate']}: reach {r['reachable_dimension']}, return {r['return_rank']}"
            f" of {r['one_seat_bound']}, {r['distinct_eigenvalues']} of {r['dimension']}"
            f" eigenvalues distinct" for r in rows))
    for row in summary["coupling_scan"]:
        values = sorted({x for r in row["ranks"] for x in r})
        print(f"N = {row['n']} gamma = {tuple(row['gamma'])} J = {row['coupling']}: map values {values}, "
              f"{row['distinct_eigenvalues']} distinct eigenvalues, frozen-root kernel dims "
              f"{row['frozen_root_kernel_dims']}")
    for gamma, coefficients in summary["n3_centre_return_recurrence"].items():
        print(f"N = 3 gamma = {gamma}: centre-return recurrence, from s^4 down: {coefficients}")
    print(f"wrote {', '.join(summary['figures'])} and the JSON to {args.output_dir}")


if __name__ == "__main__":
    main()
