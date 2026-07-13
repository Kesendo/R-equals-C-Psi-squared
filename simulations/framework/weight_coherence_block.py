"""The (w_ket, w_bra) joint-popcount coherence pencil: the Python mirror of the C# canonical
`compute/RCPsiSquared.Core/F89PathK/WeightCoherenceBlock.cs`.

The Z-dephasing Liouvillian L = -i[H, rho] + D[rho] of an open XY chain preserves the joint
popcount (popcount(ket), popcount(bra)), so it acts on each coherence sector |a><b| as a finite
block of size C(N, w_ket) * C(N, w_bra). This primitive builds that block in the C# canonical
convention, so exploration scripts stop re-implementing it (four independent transcriptions
existed by 2026-07; they stay untouched as frozen engines of record and are pinned against this
module by `simulations/weight_coherence_parity.py`).

Conventions (= the C# canonical; every axis stated because the four legacy transcriptions differ):

- **Bit convention: little-endian.** Bit s of a config mask = site s (`1 << s`), matching the C#
  builder and all four F89 legacy scripts. NOTE the sibling primitive `coherence_block.py` in this
  package is the OPPOSITE (big-endian, site 0 = MSB) and uses the J-book with hop amplitude J;
  the two primitives serve different objects (that one: the (n, n+1) chromaticity block split by
  bond for F86/Q-scale work; this one: the general (w_ket, w_bra) pencil mirror of the C# core).
  Do not compose their bases without an explicit bit-reversal.
- **Basis order:** kets = configs(n, w_ket) ascending by mask, outer; bras = configs(n, w_bra)
  ascending by mask, inner; index = i_ket * C(n, w_bra) + i_bra. The legacy scripts enumerate
  `itertools.combinations` site tuples instead; the two orders coincide for w <= 1 but not in
  general (at n = 4, w = 2, combinations order gives masks 3, 5, 9, 6, ... which is not
  ascending; small n can coincide accidentally, e.g. n = 3, w = 2), so
  `combinations_order_permutation` is provided as the exact adapter and is always applied.
- **Pencil:** L(q) = diag(A) - 2j * q_octic * K + (delta term), with A = -2 * gamma * n_diff
  (the Absorption-Theorem diagonal, n_diff = popcount(ket ^ bra)) and
  K = kron(H_ket, I) - kron(I, H_bra) the real unit-hop difference. Ket excitations hop
  -2iq, bra excitations +2iq (nearest-neighbour, Pauli-excluded), exactly the C# entries.
- **q-book:** `book='octic'` (default) reads q with the C# / octic normalization (hop -2iq,
  q = J of H = J*Sigma(XX+YY)); `book='carrier'` reads q as the carrier-clock knob
  q_carrier = 2 * q_octic (the unit-hop convention of the seed-existence census builders).
  The book converts the KNOB at the API boundary (q_octic = q / 2, an exact IEEE scaling) and
  nothing else, so it applies uniformly to the hop and the delta*ZZ frequency. Canonical
  statement of the factor 2: docs/GLOSSARY.md, "The coupling ratio q and Q".
- **gamma is a live axis** (deliberate divergence from the C# signature, which fixes gamma = 1):
  the diagonal is the one-term product -2 * gamma * n_diff, and gamma = 0 gives the pure
  Hamiltonian pencil.
- **delta:** the XXZ ZZ-anisotropy of H = J*Sigma(XX+YY) + J*delta*Sigma(ZZ); a DIAGONAL
  Hermitian frequency -1j * q_octic * delta * (zz(ket) - zz(bra)) that leaves
  Re(diag) = -2 * gamma * n_diff untouched (the C# (q, delta) overload). The C# per-site
  longitudinal-field overload is NOT mirrored here (no Python consumer asks for it).
- **Leg-adjoint relation** (each block is complex-symmetric, L^T = L): with tau the swap
  |a><b| -> |b><a| mapping block (u, v) to block (v, u),
      tau L_(u,v)(q) tau^-1 = conj(L_(v,u)(conj(q))),
  which for REAL q reads tau L_(u,v) tau^-1 = conj(L_(v,u)) = L_(v,u)^dagger. For complex q the
  conjugation also sends q -> conj(q); do not use the real-q form off the real axis.

The three legacy sign conventions this module does NOT adopt (see the parity gate's adapter
table): `resonant_n_twinning.blocks` and `cross_triple_orthogonality.block_space` return the
NEGATED real hop -K in combinations order; `seed_existence_nullity_check.build` returns
C = -1j * K in combinations order at the carrier book.

Scope: the block BUILDERS only (the C# `Build`/`Configs`/`Zz`). The C# class's permutation
carriers (`BraComplementPermutation`, `KetComplementPermutation`, `ReflectionPermutation`, the
sector-CSR assembler) and the longitudinal-field overload stay C#-side; port them only when a
Python consumer actually asks. Decision rule for new work: reach for THIS module in Python
exploration (gamma sweeps, quick pencils, adapters to the legacy conventions); reach for the C#
core for witnesses, the fold/reflection carriers, large-N sector work, and anything that lands
as claim evidence (C# witness first).

Tests: simulations/framework/tests/primitives/test_weight_coherence_block.py.
Parity pins vs the four frozen legacy builders: simulations/weight_coherence_parity.py.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np

__all__ = [
    "weight_block_configs",
    "weight_block_hop",
    "weight_block_pencil",
    "weight_block_build",
    "weight_block_zz",
    "mask_of_sites",
    "sites_of_mask",
    "combinations_order_permutation",
]


def weight_block_configs(n: int, w: int) -> list[int]:
    """All n-bit masks with exactly w set bits, ascending (the C# `Configs` order)."""
    return [m for m in range(1 << n) if bin(m).count("1") == w]


def mask_of_sites(sites) -> int:
    """Occupied-site tuple -> little-endian bitmask (site s -> bit s)."""
    m = 0
    for s in sites:
        m |= 1 << s
    return m


def sites_of_mask(mask: int, n: int) -> tuple:
    """Little-endian bitmask -> ascending occupied-site tuple."""
    return tuple(s for s in range(n) if (mask >> s) & 1)


def combinations_order_permutation(n: int, w: int) -> np.ndarray:
    """perm[i] = the ascending-mask index of the i-th `itertools.combinations(range(n), w)` tuple.

    The exact adapter between the legacy scripts' combinations order and this module's
    ascending-mask order. Identity for w <= 1; a genuine permutation from w = 2 on.
    """
    ascending = {m: i for i, m in enumerate(weight_block_configs(n, w))}
    return np.array([ascending[mask_of_sites(c)] for c in combinations(range(n), w)], dtype=int)


def weight_block_hop(n: int, w: int) -> np.ndarray:
    """Nearest-neighbour XY unit hop on the w-excitation sector (real symmetric), ascending-mask
    basis. Spec pin: at w = 1 this is the open-chain path adjacency, spectrum 2*cos(k*pi/(n+1))."""
    states = weight_block_configs(n, w)
    idx = {m: i for i, m in enumerate(states)}
    h = np.zeros((len(states), len(states)))
    for m in states:
        for s in range(n - 1):
            pair = (1 << s) | (1 << (s + 1))
            if bin(m & pair).count("1") == 1:  # exactly one of the two sites occupied
                h[idx[m ^ pair], idx[m]] += 1.0
    return h


def weight_block_pencil(n: int, w_ket: int, w_bra: int, gamma: float = 1.0):
    """The (w_ket, w_bra) pencil pieces (A, K): A the real AT diagonal -2*gamma*n_diff as a
    1-D vector over the (ket outer, bra inner) ascending-mask basis, K the real unit-hop
    difference kron(H_ket, I) - kron(I, H_bra). L(q) = diag(A) - 2j*q_octic*K (see module doc)."""
    kets = weight_block_configs(n, w_ket)
    bras = weight_block_configs(n, w_bra)
    a = np.array([-2.0 * gamma * bin(k ^ b).count("1") for k in kets for b in bras])
    k_mat = np.kron(weight_block_hop(n, w_ket), np.eye(len(bras))) - np.kron(
        np.eye(len(kets)), weight_block_hop(n, w_bra))
    return a, k_mat


def weight_block_zz(n: int, c: int) -> int:
    """zz(c) = sum over open-chain bonds of <c|Z_b Z_{b+1}|c> (+1 equal bits, -1 differing);
    even under the global bit-flip (the C# `Zz`)."""
    return sum(1 if ((c >> b) & 1) == ((c >> (b + 1)) & 1) else -1 for b in range(n - 1))


def weight_block_build(n: int, w_ket: int, w_bra: int, q, gamma: float = 1.0,
                       delta: float = 0.0, book: str = "octic") -> np.ndarray:
    """The assembled complex (w_ket, w_bra) block L, mirroring the C# `Build(n, wKet, wBra, q, delta)`
    with gamma as a live axis. `book` selects how the knob q is read ('octic' = C# hop -2iq;
    'carrier' = unit-hop knob, q_octic = q/2); the conversion is the only thing the book changes."""
    if book == "octic":
        q_octic = q
    elif book == "carrier":
        q_octic = q / 2
    else:
        raise ValueError(f"unknown q-book {book!r}: use 'octic' or 'carrier'")
    a, k_mat = weight_block_pencil(n, w_ket, w_bra, gamma)
    l_mat = np.diag(a.astype(complex)) - 2j * q_octic * k_mat
    if delta != 0.0:
        kets = weight_block_configs(n, w_ket)
        bras = weight_block_configs(n, w_bra)
        zz_diff = np.array([weight_block_zz(n, k) - weight_block_zz(n, b)
                            for k in kets for b in bras], dtype=float)
        l_mat += np.diag(-1j * q_octic * delta * zz_diff)
    return l_mat
