"""Matrix diagnostics for the conditional identity Q J Q.T + J + 2 s I = 0."""

from numbers import Integral, Real

import numpy as np
from scipy.linalg import schur, subspace_angles
from scipy.optimize import linear_sum_assignment


def make_balanced_dale_network(n=50, n_exc=25, density=0.3, seed=42):
    """Return balanced Dale weights with the original drive-sweep RNG recipe.

    Balance refers only to the E/I population; no palindrome is imposed.
    """
    if not isinstance(n, Integral) or isinstance(n, bool) or n <= 0 or n % 2:
        raise ValueError("n must be a positive even integer")
    if (not isinstance(n_exc, Integral) or isinstance(n_exc, bool)
            or n_exc != n // 2):
        raise ValueError("n_exc must equal n // 2 for balanced E/I populations")
    if not isinstance(density, Real) or not np.isfinite(density) or not 0 <= density <= 1:
        raise ValueError("density must be finite and in [0, 1]")
    if not isinstance(seed, Integral) or isinstance(seed, bool):
        raise ValueError("seed must be an integer")
    rng = np.random.RandomState(seed)
    signs = np.ones(n)
    signs[rng.choice(n, n - n_exc, replace=False)] = -1
    mask = rng.random((n, n)) < density
    np.fill_diagonal(mask, False)
    weights = rng.exponential(0.3, (n, n))
    W = np.where(mask, weights * signs[None, :], 0.0)
    maximum = np.max(np.abs(W))
    if maximum > 0:
        W /= maximum
    return W, signs


def typed_sigmoid(inputs: np.ndarray, signs: np.ndarray,
                  a_e: float = 1.3, theta_e: float = 4.0,
                  a_i: float = 2.0, theta_i: float = 3.7) -> np.ndarray:
    """Vectorized E/I sigmoid, using the drive model's per-type parameters."""
    slopes = np.where(signs > 0, a_e, a_i)
    thresholds = np.where(signs > 0, theta_e, theta_i)
    return 1.0 / (1.0 + np.exp(np.clip(-slopes * (inputs - thresholds), -500, 500)))


def solve_fixed_point(W, signs, alpha=0.3, drive=4.0, tol=1e-12,
                      max_iter=5000, *, a_e=1.3, theta_e=4.0,
                      a_i=2.0, theta_i=3.7):
    """Return (x, max|x-F(x)|) after synchronous fixed-point convergence.

    Start at x=0.3. Stop only when max|x_next-x| < tol and the fresh equation
    residual at x_next is <= tol. No convergence is assumed for arbitrary
    coupling: exhausted iterations raise with the last fresh equation residual.
    """
    W = np.asarray(W)
    signs = np.asarray(signs)
    if (W.ndim != 2 or W.shape[0] == 0 or W.shape[0] != W.shape[1]
            or not np.isrealobj(W) or not np.all(np.isfinite(W))):
        raise ValueError("W must be a nonempty finite real square matrix")
    if signs.shape != (W.shape[0],) or not np.all(np.isin(signs, [-1, 1])):
        raise ValueError("signs must have one +1 or -1 per row of W")
    for name, value in (("alpha", alpha), ("drive", drive), ("tol", tol),
                        ("a_e", a_e), ("theta_e", theta_e),
                        ("a_i", a_i), ("theta_i", theta_i)):
        if not isinstance(value, Real) or not np.isfinite(value):
            raise ValueError(f"{name} must be a finite real number")
    if tol <= 0:
        raise ValueError("tol must be positive")
    if (not isinstance(max_iter, Integral) or isinstance(max_iter, bool)
            or max_iter <= 0):
        raise ValueError("max_iter must be a positive integer")

    def evaluate(x):
        return typed_sigmoid(alpha * W @ x + drive, signs, a_e, theta_e, a_i, theta_i)

    x = np.full(W.shape[0], 0.3)
    for _ in range(max_iter):
        candidate = evaluate(x)
        update_residual = float(np.max(np.abs(candidate - x)))
        x = candidate
        last_residual = float(np.max(np.abs(x - evaluate(x))))
        if update_residual < tol and last_residual <= tol:
            return x, last_residual
    raise RuntimeError(f"fixed point did not converge in {max_iter} iterations; "
                       f"last residual={last_residual:.16g}")


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
    """Return (J, perm, s) for the balanced E/I construction, centered at -s.

    Require positive even integer n, positive finite time constants, finite
    real alpha, density in [0, 1], and an integer RandomState seed.
    """
    if not isinstance(n, Integral) or n <= 0 or n % 2:
        raise ValueError("n must be a positive even integer for balanced E/I pairs")
    for name, value in (("tau_e", tau_e), ("tau_i", tau_i),
                        ("alpha", alpha), ("density", density)):
        if not isinstance(value, Real) or not np.isfinite(value):
            raise ValueError(f"{name} must be a finite real number")
    for name, value in (("tau_e", tau_e), ("tau_i", tau_i)):
        if value <= 0:
            raise ValueError(f"{name} must be positive")
    if not 0 <= density <= 1:
        raise ValueError("density must be in [0, 1]")
    if not isinstance(seed, Integral) or isinstance(seed, bool):
        raise ValueError("seed must be an integer")
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
