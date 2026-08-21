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
- **Pencil:** L(q) = diag(A) - 2j * q_octic * K + (delta term), with
  A = -2 * sum_s gamma_s * [ket_s != bra_s] (the Absorption-Theorem diagonal, which at UNIFORM
  gamma collapses to the familiar -2 * gamma * n_diff, n_diff = popcount(ket ^ bra)) and
  K = kron(H_ket, I) - kron(I, H_bra) the real unit-hop difference. Ket excitations hop
  -2iq, bra excitations +2iq (nearest-neighbour, Pauli-excluded), exactly the C# entries.
- **q-book:** `book='octic'` (default) reads q with the C# / octic normalization (hop -2iq,
  q = J of H = J*Sigma(XX+YY)); `book='carrier'` reads q as the carrier-clock knob
  q_carrier = 2 * q_octic (the unit-hop convention of the seed-existence census builders).
  The book converts the KNOB at the API boundary (q_octic = q / 2, an exact IEEE scaling) and
  nothing else, so it applies uniformly to the hop and the delta*ZZ frequency. Canonical
  statement of the factor 2: docs/GLOSSARY.md, "The coupling ratio q and Q".
- **gamma is a live axis, scalar OR a per-site PROFILE of length n.** A scalar is the uniform
  case, where the diagonal is the one-term product -2 * gamma * n_diff and gamma = 0 gives the pure
  Hamiltonian pencil. A profile is `gamma[s]` = the rate of site s = BIT s, the little-endian
  convention above. Why it exists rather than the scalar being enough: at uniform gamma the diagonal
  entry depends only on n_diff, so the coherences of a constant-n_diff block all share it, and that
  shared value is the hypothesis under the repo's "the sector sits at Re = -2*gamma" family and under
  the edge-normality lemma. Under a profile it is a subset sum, and the sharing generically stops (a profile whose subset
  sums collide still shares; the C# gate asserts only that more than one value survives). What that costs
  SPECTRALLY is not stated here: the diagonal is not the spectrum, and
  docs/proofs/PROOF_EDGE_BLOCK_DEFECTIVE_UNDER_PROFILE.md owns the answer. The law itself is the Absorption
  Theorem's vector form, held in the repo since 2026-05-29, not anything added here.
  The scalar is a Python-only convenience and a different arithmetic ROUTE to the same real number
  (one multiply against a repeated sum); the two agree exactly where the sum is exact and may differ
  in the last bit otherwise, which the G5 tests read rather than tolerate.
- **REVERSE the profile for the big-endian siblings.** `fw.lindbladian_z_dephasing` in this package
  builds the full 4^N L with site 0 as the leftmost Kronecker factor (site l = bit n-1-l); so does
  `simulations/d10_block_closure_verify.py`, and the C# `PerBlockLiouvillianBuilder` shares the index
  convention without ever materialising the full matrix. (`simulations/edge_block_defective_ep_gate.py`
  is NOT in this list: it builds an n x n single-excitation block directly in the SITE basis, so no
  bitmask and no reversal question arise there. The big-endian sibling `coherence_block.py`, named in
  the bit-convention bullet above, is also absent for its own reason: it takes a scalar rate only, so
  it has no profile to reverse.)
  Reading a block off a full L built by a big-endian sibling and comparing it here is the
  obvious next move, and unreversed each rate lands on the mirror site: invisible at uniform gamma
  (a constant array is its own reverse) and invisible in the SPECTRUM whenever H is
  reflection-symmetric. Gated by `test_the_profile_agrees_with_the_full_liouvillian_only_when_reversed`.
- **delta:** the XXZ ZZ-anisotropy of H = J*Sigma(XX+YY) + J*delta*Sigma(ZZ); a DIAGONAL
  Hermitian frequency -1j * q_octic * delta * (zz(ket) - zz(bra)) that leaves
  Re(diag) = -2 * sum_s gamma_s * [ket_s != bra_s] untouched (the C# (q, delta) overload), under
  a profile exactly as at uniform gamma: at REAL q the delta term is imaginary and rate-free (off the
  real axis -1j * q_octic * delta * zz_diff acquires a real part, as the leg-adjoint bullet below also warns). The C# per-site
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
    "weight_block_disagreement_sum",
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


def weight_block_disagreement_sum(n: int, gamma, a: int, b: int) -> float:
    """sum_s gamma_s * [a_s != b_s], the gamma-weighted count of bra-ket disagreements. NOT the
    Absorption-Theorem rate: the pencil's diagonal ENTRY is -2 times it, and that is not an eigenvalue.
    The AT reads -Re(lambda) = 2 * sum_l gamma_l * <Delta_l> over eigenmodes, with the expectation
    <Delta_l> in [0,1]; the hopping is what turns the sharp bit into that expectation. A SCALAR gamma is the uniform case and collapses
    to gamma * n_diff, computed that way (one multiply, so the uniform route keeps the arithmetic it
    always had). A profile is summed in ascending-site order. Site s is bit s (module docstring)."""
    x = a ^ b
    if np.ndim(gamma) == 0:
        return float(gamma) * bin(x).count("1")
    g = np.asarray(gamma, dtype=float)
    if g.shape != (n,):
        raise ValueError(f"gamma profile shape {g.shape} != (n,) = ({n},)")
    return float(sum(g[s] for s in range(n) if (x >> s) & 1))


def weight_block_pencil(n: int, w_ket: int, w_bra: int, gamma=1.0):
    """The (w_ket, w_bra) pencil pieces (A, K): A the real AT diagonal as a 1-D vector over the
    (ket outer, bra inner) ascending-mask basis, K the real unit-hop difference
    kron(H_ket, I) - kron(I, H_bra). L(q) = diag(A) - 2j*q_octic*K (see module doc).

    `gamma` is either a scalar (uniform, A = -2*gamma*n_diff) or a per-site PROFILE of length n
    (A = -2 * sum_s gamma_s * [ket_s != bra_s], `weight_block_disagreement_sum`). The profile
    reaches A and nothing else: dephasing is diagonal in the computational basis, so K is
    gamma-free either way and the q-linear split survives with the profile confined to A."""
    kets = weight_block_configs(n, w_ket)
    bras = weight_block_configs(n, w_bra)
    a = np.array([-2.0 * weight_block_disagreement_sum(n, gamma, k, b) for k in kets for b in bras])
    k_mat = np.kron(weight_block_hop(n, w_ket), np.eye(len(bras))) - np.kron(
        np.eye(len(kets)), weight_block_hop(n, w_bra))
    return a, k_mat


def weight_block_zz(n: int, c: int) -> int:
    """zz(c) = sum over open-chain bonds of <c|Z_b Z_{b+1}|c> (+1 equal bits, -1 differing);
    even under the global bit-flip (the C# `Zz`)."""
    return sum(1 if ((c >> b) & 1) == ((c >> (b + 1)) & 1) else -1 for b in range(n - 1))


def weight_block_build(n: int, w_ket: int, w_bra: int, q, gamma=1.0,
                       delta: float = 0.0, book: str = "octic") -> np.ndarray:
    """The assembled complex (w_ket, w_bra) block L, mirroring the C# `Build` at field = null
    (the C# longitudinal-field knob is not mirrored here), with gamma as a live axis, scalar
    (uniform) or a per-site profile of length n. `book` selects how the knob q is read
    ('octic' = C# hop -2iq; 'carrier' = unit-hop knob, q_octic = q/2); the conversion is the only
    thing the book changes, and it does not touch gamma."""
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
