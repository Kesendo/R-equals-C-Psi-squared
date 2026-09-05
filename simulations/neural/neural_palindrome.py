"""Matrix diagnostics for the conditional identity Q J Q.T + J + 2 s I = 0."""

import numpy as np
from scipy.optimize import linear_sum_assignment


def _validated_permutation(perm):
    perm = np.asarray(perm)
    if (
        perm.ndim != 1
        or not np.issubdtype(perm.dtype, np.integer)
        or not np.array_equal(np.sort(perm), np.arange(perm.size))
    ):
        raise ValueError("perm must contain each integer index exactly once")
    return perm


def permutation_matrix(perm):
    """Return the permutation matrix with Q[i, perm[i]] = 1."""
    perm = _validated_permutation(perm)
    Q = np.zeros((perm.size, perm.size))
    Q[np.arange(perm.size), perm] = 1.0
    return Q


def is_involution(perm):
    """Return whether perm is a valid permutation whose square is identity."""
    try:
        perm = _validated_permutation(perm)
    except ValueError:
        return False
    return bool(np.array_equal(perm[perm], np.arange(perm.size)))


def _matrix_and_permutation(J, perm):
    J = np.asarray(J)
    Q = permutation_matrix(perm)
    if J.ndim != 2 or J.shape != Q.shape:
        raise ValueError("J must be square with one row per permutation index")
    return J, Q


def _normalized_residual(residual, J):
    norm_J = np.linalg.norm(J, "fro")
    # At the zero matrix, report the absolute residual instead of dividing by zero.
    return float(np.linalg.norm(residual, "fro") / (norm_J if norm_J else 1.0))


def scalar_center_residual(J, perm, s):
    """Return ||Q J Q.T + J + 2 s I||_F / ||J||_F (absolute if J is zero).

    The scalar s is supplied by the caller; no diagonal entries are discarded.
    Floating arithmetic residuals are returned without tolerance or clipping.
    """
    J, Q = _matrix_and_permutation(J, perm)
    if np.ndim(s) != 0:
        raise ValueError("s must be a scalar")
    residual = Q @ J @ Q.T + J + 2 * s * np.eye(J.shape[0])
    return _normalized_residual(residual, J)


def fitted_offdiagonal_residual(J, perm):
    """Return the legacy fitted off-diagonal residual, normalized by ||J||_F.

    Fit a separate S for each diagonal seat, then discard the diagonal. This
    instrument does not test the identity at any one scalar center s.
    At J = 0, return the absolute residual.
    """
    J, Q = _matrix_and_permutation(J, perm)
    QJQ = Q @ J @ Q.T
    S = -(np.diag(QJQ) + np.diag(J)) / 2
    residual = QJQ + J + 2 * np.diag(S)
    np.fill_diagonal(residual, 0)
    return _normalized_residual(residual, J)


def spectral_pairing_error(values, s):
    """Return the largest complex distance in a minimum-total-cost assignment.

    Match the eigenvalue multiset to -values - 2*s, preserving multiplicity.
    This reports the maximum assigned cost, not a bottleneck-optimal distance.
    """
    values = np.asarray(values, dtype=complex)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("values must be a nonempty one-dimensional array")
    if np.ndim(s) != 0:
        raise ValueError("s must be a scalar")
    targets = -values - 2 * s
    costs = np.abs(values[:, None] - targets[None, :])
    rows, columns = linear_sum_assignment(costs)
    return float(np.max(costs[rows, columns]))
