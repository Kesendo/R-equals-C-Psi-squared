"""F83/F85 closed-form Π-decomposition of M (anti-fraction, ‖M‖²,‖M_anti‖²,‖M_sym‖²)."""
from __future__ import annotations

import numpy as np

from ..pauli import _build_bilinear, _build_kbody_chain
from ..symmetry import _pauli_tuple_is_truly, _pauli_tuple_pi2_class


def predict_pi_decomposition(chain, terms):
    """F83 (k-body via F85): predict ‖M‖², ‖M_anti‖², ‖M_sym‖², and
    anti-fraction from H alone, without computing M.

    Theorem F83/F85 (proven in PROOF_F83_PI_DECOMPOSITION_RATIO and
    PROOF_F85_KBODY_GENERALIZATION):

        ‖M‖²_F        = 4·‖H_odd‖²·2^N + 8·‖H_even_nontruly‖²·2^N
        ‖M_anti‖²_F  = 2·‖H_odd‖²·2^N
        ‖M_sym‖²_F   = 2·‖H_odd‖²·2^N + 8·‖H_even_nontruly‖²·2^N
        anti-fraction = ‖M_anti‖² / ‖M‖² = 1 / (2 + 4·r)

    with r = ‖H_even_nontruly‖²/‖H_odd‖². Special cases:
        r = 0  (pure Π²-odd):              anti = 1/2  (F81 50/50)
        r = ∞  (pure Π²-even non-truly):  anti = 0    (F81 100/0)
        r = 1  (equal-Frobenius mix):     anti = 1/6  (5/6+1/6 finding)

    F85 generalization to k-body: terms can be tuples of any length k ≥ 2.
    Truly classification (M = 0 for term alone): #Y even AND #Z even.
    Π²-class determines factor c ∈ {0, 1, 2}, replacing F49's 2-body-
    coincidental n_YZ formula with the structurally correct c(k) formula.

    For 2-body, c(k) = n_YZ(k) coincidentally (Π²-odd ↔ n_YZ=1,
    Π²-even non-truly ↔ n_YZ=2). For k ≥ 3, the n_YZ formula breaks
    and only the Π²-class matters.

    γ_z-independent. Topology-independent (verified chain/ring/star/K_N
    at 2-body; chain at k = 3, 4 for higher-body).

    Args:
        terms: list of Pauli-letter tuples. Each tuple has length k ≥ 1
            and contains letters from {'I', 'X', 'Y', 'Z'}. Mixed body
            counts in the same call are supported.

    Returns:
        dict with keys:
            'M_sq', 'M_anti_sq', 'M_sym_sq', 'anti_fraction',
            'h_odd_sq', 'h_even_nontruly_sq', 'r'.
    """
    # Group non-truly terms by Π²-class (truly drops, Π²-odd, Π²-even non-truly)
    odd_terms = []
    even_nontruly_terms = []
    for term in terms:
        letters = tuple(term)
        if _pauli_tuple_is_truly(letters):
            continue
        if 'I' in letters:
            continue  # single-body or partial-identity outside F83/F85 scope
        cls = _pauli_tuple_pi2_class(letters)
        term_with_coeff = letters + (chain.J,)
        if cls == 'pi2_odd':
            odd_terms.append(term_with_coeff)
        elif cls == 'pi2_even_nontruly':
            even_nontruly_terms.append(term_with_coeff)

    # Build sub-Hamiltonians per Π²-class. Topology handling:
    #   - 2-body terms in the group → _build_bilinear with chain.bonds
    #     (chain/ring/star/K_N, F49 scope).
    #   - k-body (k ≥ 3) terms → _build_kbody_chain (chain sliding-window,
    #     F85 chain-only scope).
    # If a group contains both 2-body and k-body terms, build them
    # separately and add to avoid silently degrading non-chain topology
    # for the 2-body part.
    def _build_group(group_terms):
        if not group_terms:
            return None
        two_body = [t for t in group_terms if len(t) == 3]  # 2 letters + 1 coeff
        kbody = [t for t in group_terms if len(t) > 3]
        d = 2 ** chain.N
        H = np.zeros((d, d), dtype=complex)
        if two_body:
            two_body_clean = [(t[0], t[1], t[2]) for t in two_body]
            H = H + _build_bilinear(chain.N, chain.bonds, two_body_clean)
        if kbody:
            H = H + _build_kbody_chain(chain.N, kbody)
        return H

    H_odd = _build_group(odd_terms)
    h_odd_sq = float(np.real(np.trace(H_odd.conj().T @ H_odd))) if H_odd is not None else 0.0
    H_even = _build_group(even_nontruly_terms)
    h_even_nontruly_sq = float(np.real(np.trace(H_even.conj().T @ H_even))) if H_even is not None else 0.0

    d_pow = 2 ** chain.N
    m_sq = 4 * h_odd_sq * d_pow + 8 * h_even_nontruly_sq * d_pow
    m_anti_sq = 2 * h_odd_sq * d_pow
    m_sym_sq = m_sq - m_anti_sq

    if h_odd_sq < 1e-15:
        r = float('inf')
        anti_fraction = 0.0
    else:
        r = h_even_nontruly_sq / h_odd_sq
        anti_fraction = m_anti_sq / m_sq if m_sq > 0 else 0.0

    return {
        'M_sq': float(m_sq),
        'M_anti_sq': float(m_anti_sq),
        'M_sym_sq': float(m_sym_sq),
        'anti_fraction': float(anti_fraction),
        'h_odd_sq': float(h_odd_sq),
        'h_even_nontruly_sq': float(h_even_nontruly_sq),
        'r': float(r),
    }


def predict_pi_decomposition_anti_fraction(chain, terms):
    """F83 anti-fraction convenience wrapper around `predict_pi_decomposition`.

    Returns just the float anti-fraction = ‖M_anti‖²/‖M‖² = 1/(2+4r) for
    the given Hamiltonian terms. See `predict_pi_decomposition` for the
    full F83 closed-form prediction (M, M_anti, M_sym norms plus inputs).
    """
    return predict_pi_decomposition(chain, terms)['anti_fraction']
