"""F77 Π²-class trichotomy classifier (truly / soft / hard) for Hamiltonian terms.

Generalizes to arbitrary k-body chain Hamiltonians via term-length dispatch:
length-2 terms use the bond-bilinear builder, length-k terms (k ≥ 2) use
the chain sliding-window k-body builder. Mixed body counts are summed.
"""
from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

from ..lindblad import lindbladian_pauli_dephasing, palindrome_residual
from ..pauli import _build_bilinear, _build_kbody_chain


def _greedy_pairing_error(evals, sigma_gamma):
    """Nearest-unused-partner pass over the spectrum. Upper bound only.

    Kept because it is O(n²) with a tiny constant and bounds the exact answer
    from above, which seeds the search in `spectrum_pairing_error`.
    """
    used = np.zeros(len(evals), dtype=bool)
    max_err = 0.0
    for i in range(len(evals)):
        if used[i]:
            continue
        # written as a SUM, bit-identical to the D matrix in `spectrum_pairing_error`:
        # the subtracted form differs from it by ~1 ulp, which can prune the very edge
        # this pass used and leave the binary search with no feasible candidate.
        dists = np.abs(evals + evals[i] + 2 * sigma_gamma)
        dists[used] = np.inf
        j = int(np.argmin(dists))
        used[i] = True
        if j != i:
            used[j] = True
        max_err = max(max_err, float(dists[j]))
    return max_err


def spectrum_pairing_error(evals, sigma_gamma):
    """How far the spectrum is from being closed under λ ↦ −λ − 2Σγ.

    The palindrome asks whether the eigenvalue MULTISET is invariant under the
    mirror, so the honest error is the bottleneck distance between the spectrum
    and its mirror image: the smallest ε admitting a bijection that moves no
    eigenvalue further than ε. Computed exactly, by binary search on the
    distinct distances with a bipartite-matching feasibility test.

    A greedy nearest-partner pass answers a strictly harder question, since it
    also forces the bijection to be an involution (it pairs i and j off
    together), and greedy pairing is not optimal for the bottleneck objective
    anyway. Both effects can only INFLATE the error, so a greedy classifier is
    conservative: it can call a soft system hard, never a hard system soft. The
    gap is real on ordinary inputs (0.1 for the single-term systems IZ and ZI at
    N = 2 and 3, and for the pair IIZ+IZI at N = 4, where exact gives 0.3 against
    greedy's 0.4), but in every case swept the EXACT error itself already sits
    far above `spec_tol`, so the gap never straddles the threshold and no verdict
    is known to have turned on it. Exact costs 0.027-0.042 s at N = 5 on an idle machine, against
    1.7 s for the eigendecomposition that produced `evals` on a generic pair, so
    roughly 2.5% there. Near-diagonal systems invert that ratio, since their
    eigendecomposition is itself fast (IZ+ZI: 0.015 s eigen, 0.027 s pairing).
    Either way the absolute cost is small, which is why it is simply done
    properly here.
    """
    n = len(evals)
    if n == 0:
        return 0.0
    # |λ_i − (−λ_j − 2σ)| written as a sum: exactly symmetric (IEEE addition commutes),
    # so D[i,j] == D[j,i] to the bit and the greedy bound below is a true D entry.
    D = np.abs(evals[:, None] + evals[None, :] + 2 * sigma_gamma)
    # optimal ≤ greedy, so everything above the greedy value can be discarded
    # before the search: that also strips the dense end of the threshold range.
    cand = np.unique(D[D <= _greedy_pairing_error(evals, sigma_gamma)])
    if len(cand) == 0:
        return float(D.min())

    def feasible(idx):
        matching = maximum_bipartite_matching(csr_matrix(D <= cand[idx]),
                                              perm_type='column')
        return bool(np.all(matching >= 0))

    lo, hi = 0, len(cand) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        if feasible(mid):
            hi = mid
        else:
            lo = mid + 1
    return float(cand[lo])


def classify_pauli_pair(chain, terms, J_scale=1.0, op_tol=1e-10, spec_tol=1e-6,
                        dephase_letter='Z'):
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
        dephase_letter: dissipator letter for the Lindbladian and the
            corresponding Π palindrome operator. Default 'Z'.

    The dissipator-resonance law (verified at N=4 k=3 over 294 Z₂³-homo-
    geneous pairs, 2026-05-01): F77-hardness lives exactly in the Klein
    cell that matches the dephase_letter's Klein index. Z's Klein index
    is (0, 1); X's is (1, 0); Y's is (1, 1). SU(2)-rotation-equivalent.

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

    L_test = lindbladian_pauli_dephasing(
        H_test, [chain.gamma_0] * chain.N, dephase_letter=dephase_letter
    )
    Sigma_gamma = chain.N * chain.gamma_0
    M = palindrome_residual(L_test, Sigma_gamma, chain.N,
                            dephase_letter=dephase_letter)
    op_norm = float(np.linalg.norm(M))
    if op_norm < op_tol:
        return 'truly'

    evals = np.linalg.eigvals(L_test)
    max_err = spectrum_pairing_error(evals, Sigma_gamma)

    if max_err < spec_tol:
        return 'soft'
    return 'hard'
