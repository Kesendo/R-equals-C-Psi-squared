"""Test 1a: trace the N=6 16+16 dip cluster back to eps=0 via eigenvector continuity.

At eps=-0.83 we know the 32 cluster modes (smallest |Im| family, Re=[-0.3906,
-0.2094]). We want to identify what they are at eps=0.

Method: starting at eps=-0.83, eigendecompose L, identify cluster eigenvectors.
Step eps -> eps + delta (delta=+0.05 toward 0), eigendecompose new L, match new
modes to old cluster via max-overlap |<V_new_i | V_old_j>|. Trace each cluster
mode through eps in [-0.83, 0].

Output: for each eps, the cluster's mean |Im|, mean Re, Re-spread, plus whether
the cluster has merged with a larger family (e.g., the 30+30 cluster).

If at eps=0 the 16+16 collapses into a single 32-fold-degenerate cluster, that
is the uniform-gamma parent. If it merges into a larger cluster (e.g., 45+45),
it is a sub-family of that. If it stays distinct, it's an independent family.

Investigation only.
"""
import sys
from collections import Counter

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
N = 6
GAMMA0 = 0.05
J = 0.075
CLUSTER_SIZE = 32


def site_op(op, k, N):
    m = np.array([[1.0 + 0j]])
    for i in range(N):
        m = np.kron(m, op if i == k else I2)
    return m


def chain_H(N):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        for P in (X, Y):
            H += site_op(P, b, N) @ site_op(P, b + 1, N)
    return H


def build_L(H, gamma, N):
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(Z, k, N)
        L += gamma[k] * (np.kron(Zk, Zk.conj()) - np.eye(d * d, dtype=complex))
    return L


def sym_shape(N):
    i = np.arange(N, dtype=float)
    u = (i - (N - 1) / 2.0) ** 2
    u = u - u.mean()
    return u / np.max(np.abs(u))


def get_cluster_at_eps(eps, prev_V_cluster=None, H=None, u=None):
    """Build L at given eps, find cluster modes either by smallest |Im|
    (if prev_V_cluster is None) or by max overlap with previous cluster."""
    gamma = GAMMA0 * (1.0 + eps * u)
    L = build_L(H, gamma, N)
    sigma = float(np.sum(gamma))
    ev, V = np.linalg.eig(L)
    if prev_V_cluster is None:
        # Smallest |Im| cluster (initial step at dip)
        aim = np.abs(ev.imag)
        order = np.argsort(aim)
        # Filter to non-real
        nonreal = [i for i in order if aim[i] >= 1e-13]
        cluster_idx = np.array(nonreal[:CLUSTER_SIZE])
    else:
        # Match by max overlap |<V_new_j | V_old_k>|
        # Compute overlap matrix between all V_new columns and prev cluster columns
        # Take top-CLUSTER_SIZE V_new columns with largest projection norm onto prev_V_cluster
        # Use pseudo-inverse for non-orthogonal V
        # Compute projections via SVD-based pseudo-inverse if needed
        # Simpler: for each V_new column, compute sum_j |<V_new_i | prev_V_cluster_j>|^2
        # That tells us how much the new column projects onto the prev cluster subspace
        # But V_new is not unitary, so naive dot product gives misleading values.
        # Use Gram matrix: ||V_new_i^* @ prev_V_cluster||
        # For non-normal eigenvectors, the right way: take the projector P = prev_V @ prev_Vleft
        # and see which V_new columns survive P.
        # We need prev_V_left (left eigenvectors) for proper projector.
        # Stored in prev_V_cluster_full from previous call.
        raise NotImplementedError("Use get_cluster_with_projector instead")
    return L, ev, V, cluster_idx, sigma


def get_cluster_with_projector(eps, projector, H, u):
    """Build L at eps, find CLUSTER_SIZE V columns with largest residual after
    projecting onto previous cluster subspace.

    projector: (D, D) matrix representing the previous-step cluster subspace
              (orthonormalized via SVD of prev cluster vectors)
    """
    gamma = GAMMA0 * (1.0 + eps * u)
    L = build_L(H, gamma, N)
    sigma = float(np.sum(gamma))
    ev, V = np.linalg.eig(L)

    # Project each V column onto the previous-cluster subspace and measure overlap
    # projector is the orthonormal basis (D, k) where k=CLUSTER_SIZE
    overlaps = np.linalg.norm(projector.conj().T @ V, axis=0)
    # Pick CLUSTER_SIZE columns with largest overlap
    cluster_idx = np.argsort(-overlaps)[:CLUSTER_SIZE]
    return L, ev, V, cluster_idx, sigma


def main():
    H = J * chain_H(N)
    u = sym_shape(N)

    eps_list = np.arange(-0.83, 0.01, 0.05)  # -0.83 to 0 in steps of 0.05
    eps_list = np.round(eps_list, 4)

    print(f"N={N}, tracking 32-mode cluster backward from eps=-0.83 to eps=0")
    print(f"step Δeps = 0.05, total {len(eps_list)} points\n")
    print(f"{'eps':>8} {'sigma':>8} {'<Re>':>10} {'Re-spread':>10} "
          f"{'<|Im|>':>12} {'|Im|-spread':>12} {'overlap':>10}")
    print("-" * 90)

    # Initial step at eps=-0.83: smallest |Im| cluster
    eps0 = eps_list[0]
    gamma0 = GAMMA0 * (1.0 + eps0 * u)
    L0 = build_L(H, gamma0, N)
    sigma0 = float(np.sum(gamma0))
    ev0, V0 = np.linalg.eig(L0)
    aim0 = np.abs(ev0.imag)
    order0 = np.argsort(aim0)
    nonreal0 = [i for i in order0 if aim0[i] >= 1e-13]
    cluster_idx0 = np.array(nonreal0[:CLUSTER_SIZE])

    eigvals_cluster = ev0[cluster_idx0]
    Re_mean = float(np.mean(eigvals_cluster.real))
    Re_spread = float(np.std(eigvals_cluster.real))
    Im_mean = float(np.mean(np.abs(eigvals_cluster.imag)))
    Im_spread = float(np.std(np.abs(eigvals_cluster.imag)))
    # Orthonormalize cluster vectors for projector
    V_cluster_raw = V0[:, cluster_idx0]
    # SVD to get orthonormal basis for cluster subspace
    Q0, _, _ = np.linalg.svd(V_cluster_raw, full_matrices=False)
    print(f"{eps0:>+8.4f} {sigma0:>8.4f} {Re_mean:>+10.4f} {Re_spread:>10.4e} "
          f"{Im_mean:>12.4e} {Im_spread:>12.4e} {'(start)':>10}")
    sys.stdout.flush()

    Q_prev = Q0
    history = [(float(eps0), Re_mean, Re_spread, Im_mean, Im_spread, 1.0)]

    for eps in eps_list[1:]:
        gamma = GAMMA0 * (1.0 + eps * u)
        L = build_L(H, gamma, N)
        sigma = float(np.sum(gamma))
        ev, V = np.linalg.eig(L)
        # Overlap of each V column onto Q_prev subspace: ||Q_prev^H V_i||
        overlaps = np.linalg.norm(Q_prev.conj().T @ V, axis=0)
        # Normalize V columns (V may have non-unit columns) so overlap is meaningful
        V_norms = np.linalg.norm(V, axis=0)
        overlaps_normalized = overlaps / np.maximum(V_norms, 1e-30)
        # Pick CLUSTER_SIZE columns with largest overlap
        cluster_idx = np.argsort(-overlaps_normalized)[:CLUSTER_SIZE]
        eigvals_cluster = ev[cluster_idx]
        Re_mean = float(np.mean(eigvals_cluster.real))
        Re_spread = float(np.std(eigvals_cluster.real))
        Im_mean = float(np.mean(np.abs(eigvals_cluster.imag)))
        Im_spread = float(np.std(np.abs(eigvals_cluster.imag)))
        # Measure subspace overlap: how well does Q_prev still span the new cluster?
        V_cluster_new = V[:, cluster_idx]
        Q_new_basis, _, _ = np.linalg.svd(V_cluster_new, full_matrices=False)
        # Principal angle: smallest singular value of Q_prev^H Q_new_basis
        # If subspaces coincide, all singular values are 1; if disjoint, 0.
        Sv = np.linalg.svd(Q_prev.conj().T @ Q_new_basis, compute_uv=False)
        subspace_overlap = float(np.mean(Sv))
        print(f"{eps:>+8.4f} {sigma:>8.4f} {Re_mean:>+10.4f} {Re_spread:>10.4e} "
              f"{Im_mean:>12.4e} {Im_spread:>12.4e} {subspace_overlap:>10.4f}")
        sys.stdout.flush()
        history.append((float(eps), Re_mean, Re_spread, Im_mean, Im_spread,
                        subspace_overlap))
        Q_prev = Q_new_basis

    # Summary
    print("\n--- summary ---")
    print(f"  eps=-0.83 (start): <Re>={history[0][1]:+.4f}, Re-spread={history[0][2]:.3e}, "
          f"<|Im|>={history[0][3]:.3e}")
    eps_zero_row = history[-1]
    print(f"  eps=0 (end):       <Re>={eps_zero_row[1]:+.4f}, Re-spread={eps_zero_row[2]:.3e}, "
          f"<|Im|>={eps_zero_row[3]:.3e}, subspace_overlap={eps_zero_row[5]:.4f}")
    if eps_zero_row[2] < 1e-3:
        print(f"  -> at eps=0 the cluster is Re-degenerate (single Re value)")
    else:
        print(f"  -> at eps=0 the cluster has Re-spread {eps_zero_row[2]:.3e} "
              f"(multiple Re values)")


if __name__ == "__main__":
    main()
