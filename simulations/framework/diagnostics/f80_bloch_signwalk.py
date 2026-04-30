"""F80 structural-identity prediction of M's spectrum (chain Π²-odd 2-body)."""
from __future__ import annotations

from collections import Counter

import numpy as np

from ..pauli import _build_bilinear
from ..symmetry import _pauli_pair_is_truly


def predict_M_spectrum_pi2_odd(chain, terms, c=1.0):
    """F80 structural-identity prediction of M's spectrum (chain Π²-odd 2-body).

    Theorem F80 (verified bit-exact at N=3..7, all four pure Π²-odd Pauli
    pairs and their sums, see PROOF_F80_BLOCH_SIGNWALK):

        Spec(M)_{nontrivial} = { ±2i · λ : λ ∈ Spec(H_non-truly) }
        mult_M(±2i·λ) = mult_H_non-truly(λ) · 2^N

    H_non-truly is the input Hamiltonian with truly-class bilinears
    (XX, YY, ZZ, or any (a, b) ⊆ {I, X}) dropped, since they contribute
    M_term = 0 and are invisible to the residual.

    The remaining bilinears must all be Π²-odd 2-body: P, Q ∈ {X, Y, Z}
    with bit_b(P) + bit_b(Q) ≡ 1 (mod 2), i.e., one of (X,Y), (X,Z),
    (Y,X), (Z,X). Π²-even non-truly bilinears (Y,Z) and (Z,Y) are not in
    F80's verified scope and are rejected.

    Args:
        terms: list of (a, b) Pauli letter tuples. Truly terms are dropped;
            Π²-odd 2-body terms are summed into H_non-truly.
        c: bond coupling magnitude (default 1.0); each term is scaled by c.

    Returns:
        dict {eigenvalue (complex, purely imaginary): multiplicity (int)}
        for M's spectrum. If H_non-truly = 0 (all terms truly), returns
        {0+0j: 4^N}; the trivially-paired truly case.

    Raises:
        ValueError: if any non-truly term is not Π²-odd 2-body, or
            contains an identity letter (single-body falls under F78).
    """
    # Filter truly, validate remaining as Π²-odd 2-body
    non_truly_terms = []
    for (a, b) in terms:
        if _pauli_pair_is_truly(a, b):
            continue
        if 'I' in (a, b):
            raise ValueError(
                f"term ({a},{b}) contains identity (single-body); "
                "F80 covers chain Π²-odd 2-body only; see F78 for single-body"
            )
        from ..pauli import bit_b, _resolve
        ab_idx = _resolve(a)
        bb_idx = _resolve(b)
        parity = (bit_b(ab_idx) + bit_b(bb_idx)) % 2
        if parity != 1:
            raise ValueError(
                f"term ({a},{b}) is Π²-even non-truly (parity 0); "
                "F80 covers Π²-odd only; Π²-even non-truly clusters "
                "are richer and not in F80's verified scope"
            )
        non_truly_terms.append((a, b, c))

    if not non_truly_terms:
        return {0 + 0j: 4 ** chain.N}

    # Build H_non-truly and compute many-body spectrum
    H_nt = _build_bilinear(chain.N, chain.bonds, non_truly_terms)
    H_evs = np.linalg.eigvalsh(H_nt)
    h_counts = Counter(round(float(e), 10) for e in H_evs)

    # Translate to M-spectrum: Spec(M) = ±2i · Spec(H), mult ×2^N
    # H is Hermitian so eigenvalues are real; ±2i·λ for each H eigenvalue.
    # Note: H already has ±λ pairs (particle-hole symmetry of Majorana
    # bilinear), so 2i·λ already carries the ± structure.
    bra_factor = 2 ** chain.N
    m_spectrum = {}
    for h_ev, h_mult in h_counts.items():
        m_ev = 2j * h_ev
        m_spectrum[m_ev] = h_mult * bra_factor
    return m_spectrum
