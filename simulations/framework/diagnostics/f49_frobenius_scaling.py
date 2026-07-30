"""F49/F85 Π²-class Frobenius identity for ‖M‖² closed forms."""
from __future__ import annotations

import numpy as np

from ..lindblad import palindrome_residual_norm_squared_factor_graph
from ..pauli import _build_bilinear, _require_placeable_body_counts
from ..symmetry import _pauli_tuple_is_truly, _pauli_tuple_pi2_class


def predict_residual_norm_squared(chain, c_H, hamiltonian_class='main'):
    """Closed-form ‖M(N, G)‖² = c_H · F(N, G) without computing M.

    Uses palindrome_residual_norm_squared_factor_graph with this chain's
    topology invariants (B, D2). For chains, B = N-1 and D2 = 4N - 6.
    """
    factor = palindrome_residual_norm_squared_factor_graph(
        chain.N, chain.B, chain.D2, hamiltonian_class
    )
    return c_H * factor


def predict_residual_norm_squared_from_terms(chain, terms, gamma_t1=None):
    """Closed-form ‖M‖² from terms via F85 Π²-class Frobenius identity.

    Theorem F85 (proved in PROOF_F85_KBODY_GENERALIZATION) generalizes
    the F49 formula to k-body terms via Π²-class:

        ‖M(L_Z)‖²    = Σ_k 4·c(k)·‖H_k‖²_F·2^N
        ‖M(L_Z+T1)‖² = ‖M(L_Z)‖² + 4^(N−1) · [3·Σ γT1² + 4·(Σ γT1)²]

    where c(k) is the F85 Π²-class factor:

        c(truly term)              = 0     (M = 0 by Master Lemma)
        c(Π²-odd non-truly)        = 1     (factor 4·2^N)
        c(Π²-even non-truly)       = 2     (factor 8·2^N)

    Truly criterion (k-body): #Y even AND #Z even. At 2-body this
    reduces to "a == b OR {a,b} ⊆ {I,X}" (matches `_pauli_pair_is_truly`).

    Each non-truly term contributes independently because distinct
    Pauli strings within a Π²-class are orthogonal in Frobenius, AND
    Π²-classes are mutually orthogonal in operator space.

    For 2-body, F85's c(k) coincides with F49's n_YZ_k (Π²-odd ↔
    n_YZ=1, Π²-even non-truly ↔ n_YZ=2). For k ≥ 3, n_YZ is no
    longer the determining quantity (e.g., YYY has n_YZ=3 but c=1).

    T1 contribution is Hamiltonian-independent, γ_Z-independent, and
    orthogonal to the Hamiltonian palindrome residual (no cross-term).
    Verified at N=3..6 for arbitrary {γ_T1_l} distributions.

    Topology: 2-body uses the chain's bond graph (chain/ring/star/K_N
    all supported via F49); k ≥ 3 uses chain sliding-window semantics
    (F85 chain-only). Mixed-body term lists may degrade non-chain
    topology silently for the k-body group.

    Args:
        terms: list of Pauli-letter tuples. Each tuple has length
            2 ≤ k ≤ chain.N; a body count outside that range raises, since
            neither builder can place it. Mixed body counts in the same call
            are supported, and `terms` is materialised once so a generator
            gives the same answer as the list.
        gamma_t1: optional T1 amplitude-damping rates. Scalar (uniform
                  across sites), list of length N, or None / 0 (no T1).
    """
    # F85 Hamiltonian palindrome part. Decompose into:
    #   1. Drop per-term truly-class terms via _pauli_tuple_is_truly
    #      (M = 0 by Master Lemma; truly = #Y even AND #Z even).
    #   2. Group remaining terms by Π²-class (pi2_odd or pi2_even_nontruly).
    #      For 2-body, this matches F49's n_YZ grouping (n_YZ=1 ↔ Π²-odd,
    #      n_YZ=2 ↔ Π²-even non-truly). For k-body (k ≥ 3), the n_YZ-based
    #      F49 formula breaks; only Π²-class determines the factor c ∈ {1, 2}.
    #   3. Within each class, Pauli strings can overlap so per-class Frobenius
    #      norm is the norm of the combined sub-Hamiltonian.
    #   4. Sum per-class contributions: factor 4·2^N for Π²-odd (c=1) and
    #      8·2^N for Π²-even non-truly (c=2).
    # Materialise once: the loop below walks `terms` a second time for the
    # other Π²-class, and a generator was empty by then, silently dropping a
    # whole class (XY+YZ returned 512 instead of 1536).
    terms = [tuple(t) for t in terms]
    _require_placeable_body_counts(chain, terms)
    z_part = 0.0
    for cls in ('pi2_odd', 'pi2_even_nontruly'):
        group_terms = []
        for term in terms:
            letters = tuple(term)
            if _pauli_tuple_is_truly(letters):
                continue
            if _pauli_tuple_pi2_class(letters) != cls:
                continue
            group_terms.append(letters + (chain.J,))
        if not group_terms:
            continue
        # Build sub-Hamiltonian per Π²-class. 2-body terms use the
        # bond graph (chain/ring/star/K_N via F49); k-body terms use
        # chain sliding-window (F85 chain-only). If both present in
        # the same group, build them separately and add (preserves
        # non-chain topology for the 2-body part).
        two_body = [t for t in group_terms if len(t) == 3]
        kbody = [t for t in group_terms if len(t) > 3]
        d = 2 ** chain.N
        H_group = np.zeros((d, d), dtype=complex)
        if two_body:
            two_body_clean = [(t[0], t[1], t[2]) for t in two_body]
            H_group = H_group + _build_bilinear(chain.N, chain.bonds, two_body_clean)
        if kbody:
            from ..pauli import _build_kbody_chain
            H_group = H_group + _build_kbody_chain(chain.N, kbody)
        H_group_frob_sq = float(np.real(np.trace(H_group.conj().T @ H_group)))
        c_factor = 1 if cls == 'pi2_odd' else 2
        z_part += 4 * (2 ** chain.N) * c_factor * H_group_frob_sq

    # T1 dissipator contribution (Hamiltonian-independent, additive)
    if gamma_t1 is None:
        t1_part = 0.0
    else:
        if np.isscalar(gamma_t1):
            gamma_t1_l = [float(gamma_t1)] * chain.N
        else:
            gamma_t1_l = [float(g) for g in gamma_t1]
            if len(gamma_t1_l) != chain.N:
                raise ValueError(
                    f"gamma_t1 list length {len(gamma_t1_l)} != N {chain.N}"
                )
        sum_g2 = sum(g * g for g in gamma_t1_l)
        sum_g = sum(gamma_t1_l)
        t1_part = (4 ** (chain.N - 1)) * (3 * sum_g2 + 4 * sum_g * sum_g)

    return z_part + t1_part
