"""Import-safe multiset diagnostics for palindromic spectra."""

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching


def multiset_reflection_error(eigenvalues, midpoint):
    """Bottleneck distance to reflection about a fixed midpoint.

    The perfect-matching condition preserves algebraic multiplicity.  A
    nearest-neighbour check does not: several source eigenvalues may otherwise
    reuse the same reflected target and hide a multiplicity mismatch.
    """
    eigenvalues = np.asarray(eigenvalues, dtype=complex)
    if len(eigenvalues) == 0:
        return 0.0
    reflected = 2 * midpoint - eigenvalues
    distances = np.abs(eigenvalues[:, None] - reflected[None, :])
    candidates = np.unique(distances)
    lo, hi = 0, len(candidates) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        matching = maximum_bipartite_matching(
            csr_matrix(distances <= candidates[mid]), perm_type="column")
        if (matching >= 0).all():
            hi = mid
        else:
            lo = mid + 1
    return float(candidates[lo])
