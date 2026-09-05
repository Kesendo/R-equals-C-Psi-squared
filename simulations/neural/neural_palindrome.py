"""Matrix diagnostics for the conditional identity Q J Q.T + J + 2 s I = 0."""

import numpy as np
from scipy.linalg import schur, subspace_angles
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


def build_exact_weights(n, n_exc, tau_e, tau_i, density=0.3, seed=42):
    """Return (W, signs, perm) using veffect_exact's original RandomState recipe.

    Equal E/I populations give a fixed-point-free involution and the weight
    identity W[Q(i), Q(j)] = -(tau[Q(i)] / tau[i]) * W[i, j]. The original
    upper-triangle iteration, including partner overwrites, is preserved.
    Unbalanced populations are accepted for legacy callers, but need not give
    the scalar-center identity, since unpaired seats retain their own leak.
    """
    rng = np.random.RandomState(seed)
    signs = np.ones(n)
    inh_idx = rng.choice(n, n - n_exc, replace=False)
    signs[inh_idx] = -1
    e_idx = list(np.where(signs > 0)[0])
    i_idx = list(np.where(signs < 0)[0])
    perm = np.arange(n)
    for k in range(min(len(e_idx), len(i_idx))):
        perm[e_idx[k]] = i_idx[k]
        perm[i_idx[k]] = e_idx[k]

    W = np.zeros((n, n))
    mask = rng.random((n, n)) < density
    np.fill_diagonal(mask, False)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            qi, qj = perm[i], perm[j]
            if i < j:
                if mask[i, j]:
                    base = rng.exponential(0.3)
                    W[i, j] = signs[j] * base
                    row_tau = tau_e if signs[i] > 0 else tau_i
                    partner_tau = tau_e if signs[qi] > 0 else tau_i
                    W[qi, qj] = -(partner_tau / row_tau) * W[i, j]
    mx = np.max(np.abs(W))
    if mx > 0:
        W /= mx
    return W, signs, perm


def build_linear_jacobian(W, signs, tau_e, tau_i, alpha):
    """Build J[i,i] = -1/tau[i], J[i,j] = alpha*W[i,j]/tau[i]."""
    n = len(signs)
    J = np.zeros((n, n))
    for i in range(n):
        row_tau = tau_e if signs[i] > 0 else tau_i
        J[i, i] = -1.0 / row_tau
        for j in range(n):
            if i != j:
                J[i, j] = alpha * W[i, j] / row_tau
    return J


def make_exact_network(n=10, tau_e=5.0, tau_i=10.0, alpha=0.5,
                       seed=42, density=0.3):
    """Return (J, perm, s) for the balanced E/I construction, centered at -s."""
    if n <= 0 or n % 2:
        raise ValueError("n must be a positive even number for balanced E/I pairs")
    if tau_e <= 0 or tau_i <= 0:
        raise ValueError("time constants must be positive")
    W, signs, perm = build_exact_weights(n, n // 2, tau_e, tau_i, density, seed)
    J = build_linear_jacobian(W, signs, tau_e, tau_i, alpha)
    return J, perm, 0.5 * (1 / tau_e + 1 / tau_i)


def exact_ensemble_census():
    """Return (alpha, scalar-pass, oscillatory, unstable) for seeds 0..199.

    Fixed inputs: N=10, tau_E=5, tau_I=10, density=0.3. Scalar-pass means
    relative Frobenius residual < 1e-13; oscillatory means any |Im(mu)| >
    1e-8. Instability uses max Re(mu) > 0, with no stability tolerance.
    These are finite-ensemble numerical counts, not distributional claims.
    """
    rows = []
    for alpha in (0.5, 1.5, 3.0, 5.0, 10.0):
        exact = oscillatory = unstable = 0
        for seed in range(200):
            J, perm, s = make_exact_network(alpha=alpha, seed=seed)
            values = np.linalg.eigvals(J)
            exact += int(scalar_center_residual(J, perm, s) < 1e-13)
            oscillatory += int(np.any(np.abs(values.imag) > 1e-8))
            unstable += int(np.max(values.real) > 0)
        rows.append((alpha, exact, oscillatory, unstable))
    return rows


def _eigenvalue_clusters(values, tol):
    """Connected components under the absolute complex-distance tolerance."""
    remaining = set(range(len(values)))
    clusters = []
    while remaining:
        indices = [min(remaining)]
        remaining.remove(indices[0])
        for i in indices:
            neighbors = sorted(j for j in remaining if abs(values[i] - values[j]) <= tol)
            indices.extend(neighbors)
            remaining.difference_update(neighbors)
        clusters.append(values[indices])
    return clusters


def partner_subspace_error(J, perm, s, cluster_tol=1e-7):
    """Largest sine of principal angles between Q U_C and U_{-C-2s}.

    Complex eigenvalues are clustered by connected components at absolute
    distance cluster_tol (default 1e-7 in eigenvalue units). Ordered complex
    Schur vectors span each cluster's invariant subspace, including algebraic
    multiplicity and generalized eigenvectors at defective eigenvalues.
    The comparison is independent of the basis within each selected subspace.

    A tolerance cluster is a numerical resolution choice, not proof of exact
    degeneracy. A missing partner or wrong selected algebraic dimension raises
    ValueError rather than comparing unequal-dimensional subspaces. The empty
    matrix has no angles and returns 0.0.
    """
    J, Q = _matrix_and_permutation(J, perm)
    if np.ndim(s) != 0 or not np.isfinite(s):
        raise ValueError("s must be a finite scalar")
    if not np.isfinite(cluster_tol) or cluster_tol <= 0:
        raise ValueError("cluster_tol must be finite and positive")
    if J.shape[0] == 0:
        return 0.0
    T, _ = schur(J, output="complex")
    clusters = _eigenvalue_clusters(np.diag(T), cluster_tol)

    def invariant_subspace(centers, multiplicity):
        _, vectors, selected = schur(
            J, output="complex",
            sort=lambda value: bool(np.any(np.abs(value - centers) <= cluster_tol)),
        )
        if selected != multiplicity:
            raise ValueError(
                f"Schur selected dimension {selected} differs from "
                f"algebraic cluster multiplicity {multiplicity}"
            )
        return vectors[:, :selected]

    largest = 0.0
    for cluster in clusters:
        source = invariant_subspace(cluster, len(cluster))
        partner = invariant_subspace(-cluster - 2 * s, len(cluster))
        angles = subspace_angles(Q @ source, partner)
        if angles.size:
            largest = max(largest, float(np.max(np.sin(angles))))
    return largest
