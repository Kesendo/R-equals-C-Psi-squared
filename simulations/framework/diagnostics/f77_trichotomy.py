"""F77 Π²-class trichotomy classifier (truly / soft / hard) for Hamiltonian terms."""
from __future__ import annotations

import numpy as np

from ..lindblad import lindbladian_z_dephasing, palindrome_residual
from ..pauli import _build_bilinear


def classify_pauli_pair(chain, terms, J_scale=1.0, op_tol=1e-10, spec_tol=1e-6):
    """Trichotomy classification (truly / soft / hard) for a Pauli-pair H.

    Builds H from the given terms on this chain's bonds, computes the
    palindrome residual M and the spectrum-pairing error.

    Args:
        terms: list of (a, b) letter tuples, e.g. [('X','X'), ('Y','Y')].
        J_scale: bond coupling for the test Hamiltonian (default 1.0).

    Returns:
        'truly' | 'soft' | 'hard'
    """
    bilinear = [(t[0], t[1], J_scale) for t in terms]
    H_test = _build_bilinear(chain.N, chain.bonds, bilinear)
    L_test = lindbladian_z_dephasing(H_test, [chain.gamma_0] * chain.N)
    Sigma_gamma = chain.N * chain.gamma_0
    M = palindrome_residual(L_test, Sigma_gamma, chain.N)
    op_norm = float(np.linalg.norm(M))
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
    if op_norm < op_tol:
        return 'truly'
    if max_err < spec_tol:
        return 'soft'
    return 'hard'
