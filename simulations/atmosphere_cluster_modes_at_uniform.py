"""Test B-redux: mode-ID of the 16+16 cluster at uniform gamma (eps=0).

From the continuity scan we know: the 32-mode cluster at eps=-0.83 (|Im|=7e-5,
Re=[-0.343, -0.257]) continues to eps=0 as a 32-mode invariant subspace at
|Im|=0.0223, same Re-split (~+/-0.043 around -sigma), 99.91% subspace overlap.

This script identifies the cluster at eps=0 by overlap with the eps=-0.83 cluster
(rather than by smallest |Im|, since at eps=0 the cluster is NOT among the
smallest-|Im| modes), then projects each cluster eigenvector onto the Pauli basis
and reports dominant strings, XY-weight distribution, Pi2-parity distribution.

If at eps=0 the dominant Pauli structure matches eps=-0.83 (XY-wt 1 vs 5 split
in two sub-clusters), the cluster identity is intrinsic. If different, the
identity is eps-dependent.

Investigation only.
"""
import sys
import itertools
from collections import Counter

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}
LABELS = ('I', 'X', 'Y', 'Z')
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


def pauli_basis(N):
    D = 4 ** N
    P = np.empty((D, D), dtype=complex)
    labels = []
    for idx, letters in enumerate(itertools.product(LABELS, repeat=N)):
        op = np.array([[1.0 + 0j]])
        for ch in letters:
            op = np.kron(op, PAULIS[ch])
        P[:, idx] = op.flatten() / np.sqrt(2 ** N)
        labels.append(''.join(letters))
    return P, labels


def xy_weight(lbl):
    return sum(1 for c in lbl if c in 'XY')


def pi2_parity(lbl):
    return sum(1 for c in lbl if c in 'YZ') % 2


def spatial_palindrome(lbl):
    return int(lbl == lbl[::-1])


def main():
    H = J * chain_H(N)
    u = sym_shape(N)
    print(f"N={N}, J={J}, gamma0={GAMMA0}")

    # Step 1: identify cluster at eps=-0.83 by smallest |Im|
    print("\n--- Step 1: cluster at eps=-0.83 (dip) ---")
    eps_dip = -0.83
    gamma_dip = GAMMA0 * (1.0 + eps_dip * u)
    L_dip = build_L(H, gamma_dip, N)
    print(f"  building L_dip, eigendecomposing...")
    sys.stdout.flush()
    ev_dip, V_dip = np.linalg.eig(L_dip)
    aim_dip = np.abs(ev_dip.imag)
    order_dip = np.argsort(aim_dip)
    nonreal_dip = [i for i in order_dip if aim_dip[i] >= 1e-13]
    cluster_idx_dip = np.array(nonreal_dip[:CLUSTER_SIZE])
    print(f"  cluster at dip: <Re>={np.mean(ev_dip[cluster_idx_dip].real):+.4f}, "
          f"<|Im|>={np.mean(np.abs(ev_dip[cluster_idx_dip].imag)):.4e}")

    # Step 2: build cluster subspace projector
    Vc_dip = V_dip[:, cluster_idx_dip]
    Q_dip, _, _ = np.linalg.svd(Vc_dip, full_matrices=False)
    print(f"  cluster subspace orthonormalized (dim {Q_dip.shape[1]})")

    # Free up memory
    del L_dip
    del V_dip

    # Step 3: identify cluster at eps=0 by overlap with Q_dip
    print("\n--- Step 2: cluster at eps=0 (uniform gamma) via overlap ---")
    gamma_uniform = np.full(N, GAMMA0)
    L_uniform = build_L(H, gamma_uniform, N)
    print(f"  building L_uniform, eigendecomposing...")
    sys.stdout.flush()
    ev_uniform, V_uniform = np.linalg.eig(L_uniform)
    del L_uniform

    # Find columns with largest overlap with Q_dip
    overlaps = np.linalg.norm(Q_dip.conj().T @ V_uniform, axis=0)
    V_norms = np.linalg.norm(V_uniform, axis=0)
    overlaps_norm = overlaps / np.maximum(V_norms, 1e-30)
    cluster_idx_u = np.argsort(-overlaps_norm)[:CLUSTER_SIZE]
    eigvals_u = ev_uniform[cluster_idx_u]
    print(f"  cluster at uniform: <Re>={np.mean(eigvals_u.real):+.4f}, "
          f"<|Im|>={np.mean(np.abs(eigvals_u.imag)):.4e}")
    print(f"  Re values found:")
    re_round = np.round(eigvals_u.real, 4)
    re_count = Counter(re_round.tolist())
    for re, c in sorted(re_count.items()):
        print(f"    Re={re:+.4f}: {c} modes")
    aim_u = np.abs(eigvals_u.imag)
    print(f"  |Im| range in cluster: [{aim_u.min():.5e}, {aim_u.max():.5e}]")

    # Step 4: project cluster eigenvectors onto Pauli basis
    print("\n--- Step 3: Pauli decomposition of uniform-gamma cluster ---")
    Vc_u = V_uniform[:, cluster_idx_u]
    P_basis, labels = pauli_basis(N)
    # Pauli basis vectors are orthonormal (already normalized by 1/sqrt(2^N))
    C = P_basis.conj().T @ Vc_u  # (4^N, CLUSTER_SIZE)
    norms = np.sum(np.abs(C) ** 2, axis=0)
    print(f"  decomposition norm-check (should be ~1 each): "
          f"mean={np.mean(norms):.4f}, min={np.min(norms):.4f}, max={np.max(norms):.4f}")

    # Aggregate weight per Pauli string
    total_weight = np.sum(np.abs(C) ** 2, axis=1)
    top_alpha = np.argsort(-total_weight)[:15]
    print(f"  top 15 Pauli strings by total cluster weight (sum over {CLUSTER_SIZE} modes):")
    for a in top_alpha:
        w = float(total_weight[a])
        print(f"    {labels[a]}: weight {w:.3f}  "
              f"(XY-wt={xy_weight(labels[a])}, Pi2={pi2_parity(labels[a])}, "
              f"palindromic={spatial_palindrome(labels[a])})")

    # Per-mode top-Pauli classification
    xy_c = Counter()
    pi2_c = Counter()
    pal_c = Counter()
    for i in range(CLUSTER_SIZE):
        ci = np.abs(C[:, i]) ** 2
        top = int(np.argmax(ci))
        xy_c[xy_weight(labels[top])] += 1
        pi2_c[pi2_parity(labels[top])] += 1
        pal_c[spatial_palindrome(labels[top])] += 1
    print(f"  cluster sector distribution by each mode's top Pauli:")
    print(f"    XY-weight: {dict(sorted(xy_c.items()))}")
    print(f"    Pi2-parity: {dict(sorted(pi2_c.items()))}")
    print(f"    palindromic-vs-not: {dict(sorted(pal_c.items()))}")

    # Per-mode XY-weight distribution (split by Re sub-cluster)
    print(f"\n  per-mode XY-weight distribution, split by Re sub-cluster:")
    re_groups = {}
    for i in range(CLUSTER_SIZE):
        r = round(float(eigvals_u[i].real), 3)
        re_groups.setdefault(r, []).append(i)
    for r in sorted(re_groups.keys()):
        idxs = re_groups[r]
        xy_dist = Counter()
        for i in idxs:
            ci = np.abs(C[:, i]) ** 2
            top = int(np.argmax(ci))
            xy_dist[xy_weight(labels[top])] += 1
        print(f"    Re={r:+.3f} ({len(idxs)} modes): XY-wt dist {dict(sorted(xy_dist.items()))}")

    # Detail: 3 modes per Re sub-cluster, top 5 Paulis each
    print(f"\n  detail: top 5 Pauli strings for 3 modes from each Re sub-cluster:")
    for r in sorted(re_groups.keys()):
        idxs = re_groups[r]
        print(f"\n  --- Re={r:+.3f} sub-cluster ---")
        for i in idxs[:3]:
            ci = np.abs(C[:, i]) ** 2
            top5 = np.argsort(-ci)[:5]
            top_str = "; ".join(
                f"{labels[t]}({100 * ci[t]:.1f}%)" for t in top5)
            print(f"    mode {i}: Re={eigvals_u[i].real:+.4f} |Im|={abs(eigvals_u[i].imag):.4e}  top: {top_str}")


if __name__ == "__main__":
    main()
