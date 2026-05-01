"""F77 Π²-class trichotomy classifier (truly / soft / hard) for Hamiltonian terms.

Generalizes to arbitrary k-body chain Hamiltonians via term-length dispatch:
length-2 terms use the bond-bilinear builder, length-k terms (k ≥ 2) use
the chain sliding-window k-body builder. Mixed body counts are summed.
"""
from __future__ import annotations

import numpy as np

from ..lindblad import lindbladian_z_dephasing, palindrome_residual
from ..pauli import _build_bilinear, _build_kbody_chain


def classify_pauli_pair(chain, terms, J_scale=1.0, op_tol=1e-10, spec_tol=1e-6):
    """Trichotomy classification (truly / soft / hard) for a Pauli-pair or k-body H.

    Builds H from the given terms on this chain (length-2 terms via the
    bond-bilinear builder, length-k≥2 terms via the chain sliding-window
    k-body builder), computes the palindrome residual M and tests the
    spectrum-pairing error of L.

    Args:
        terms: list of letter tuples, e.g. [('X','X'), ('Y','Y')] for k=2,
            or [('X','Y','Z'), ('Y','Z','X')] for k=3, or mixed.
        J_scale: bond coupling for the test Hamiltonian (default 1.0).
        op_tol: tolerance for ‖M‖ < op_tol → truly.
        spec_tol: tolerance for spectrum-pairing closure → soft (else hard).

    Returns:
        'truly' | 'soft' | 'hard'

    Raises:
        ValueError: if any term length exceeds chain.N or is < 2.
    """
    if not terms:
        return 'truly'
    body_counts = {len(t) for t in terms}
    for k in body_counts:
        if k < 2 or k > chain.N:
            raise ValueError(
                f"term body count {k} outside [2, chain.N={chain.N}]"
            )

    # Build H summing over body-count groups
    d = 2 ** chain.N
    H_test = np.zeros((d, d), dtype=complex)
    bilinear_terms = [t for t in terms if len(t) == 2]
    kbody_terms = [t for t in terms if len(t) > 2]
    if bilinear_terms:
        bilinear = [(t[0], t[1], J_scale) for t in bilinear_terms]
        H_test = H_test + _build_bilinear(chain.N, chain.bonds, bilinear)
    if kbody_terms:
        kbody_with_coeff = [tuple(t) + (J_scale,) for t in kbody_terms]
        H_test = H_test + _build_kbody_chain(chain.N, kbody_with_coeff)

    L_test = lindbladian_z_dephasing(H_test, [chain.gamma_0] * chain.N)
    Sigma_gamma = chain.N * chain.gamma_0
    M = palindrome_residual(L_test, Sigma_gamma, chain.N)
    op_norm = float(np.linalg.norm(M))
    if op_norm < op_tol:
        return 'truly'

    evals = np.linalg.eigvals(L_test)
    used = np.zeros(len(evals), dtype=bool)
    max_err = 0.0
    for i in range(len(evals)):
        if used[i]:
            continue
        target = -evals[i] - 2 * Sigma_gamma
        dists = np.abs(evals - target)
        for j in range(len(evals)):
            if used[j]:
                dists[j] = np.inf
        best_j = int(np.argmin(dists))
        if best_j != i:
            used[i] = True
            used[best_j] = True
        else:
            used[i] = True
        max_err = max(max_err, float(dists[best_j]))

    if max_err < spec_tol:
        return 'soft'
    return 'hard'
