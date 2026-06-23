"""Mode identification: which Pauli sector do the cluster modes live in?

The EP confirmation (atmosphere_cluster_ep.py) showed the cluster is NOT an
exceptional point: the 32 cluster eigenvectors stay rank-32/32 at every eps,
including at the |Im| dip. It is a symmetry-protected degenerate near-
coalescence (N=5 reaches full coalescence; N=6 stops at |Im|~7e-5 in our grid).
The next question: which symmetry protects the degeneracy?

For N=5 at eps=-0.3357 (just before coalescence) and N=6 at eps=-0.83 (the dip):
full eigendecomposition, take the cluster eigenvectors, decompose each into the
Pauli basis. Identify the dominant Pauli strings; characterize XY-weight (bit_a,
Absorption) and Pi^2-parity (bit_b, w_YZ mod 2). The 16 (N=6) or 12 (N=5)
degenerate pairs should live in a specific sector with a recognizable structure.
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
GAMMA0 = 0.05
J = 0.075


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
    """Returns (P shape (4^N, 4^N) with vec(sigma_alpha) columns, labels list)."""
    D = 4 ** N
    P = np.empty((D, D), dtype=complex)
    labels = []
    for idx, letters in enumerate(itertools.product(LABELS, repeat=N)):
        op = np.array([[1.0 + 0j]])
        for ch in letters:
            op = np.kron(op, PAULIS[ch])
        P[:, idx] = op.flatten()
        labels.append(''.join(letters))
    return P, labels


def xy_weight(lbl):
    return sum(1 for c in lbl if c in 'XY')


def pi2_parity(lbl):
    return sum(1 for c in lbl if c in 'YZ') % 2


def spatial_palindrome(lbl):
    """Is the label palindromic (left-right symmetric)?"""
    return int(lbl == lbl[::-1])


for N, dip_eps, cluster_size in ((5, -0.3357, 24), (6, -0.83, 32)):
    print(f"=== N={N}, eps={dip_eps}, cluster of {cluster_size} modes ===")
    H = J * chain_H(N)
    u = sym_shape(N)
    gamma = GAMMA0 * (1.0 + dip_eps * u)
    ev, V = np.linalg.eig(build_L(H, gamma, N))
    aim = np.abs(ev.imag)
    order = np.argsort(aim)
    cluster_idx = [i for i in order if aim[i] >= 1e-13][:cluster_size]
    Vc = V[:, cluster_idx]
    P, labels = pauli_basis(N)
    scale = 1.0 / np.sqrt(2 ** N)
    C = scale * (P.conj().T @ Vc)
    norms = np.sum(np.abs(C) ** 2, axis=0)
    print(f"  decomposition norm-check (should be ~1 each): "
          f"mean={float(np.mean(norms)):.4f}, "
          f"min={float(np.min(norms)):.4f}, max={float(np.max(norms)):.4f}")

    # aggregate weight per Pauli string across cluster modes
    total_weight = np.sum(np.abs(C) ** 2, axis=1)
    top_alpha = np.argsort(-total_weight)[:15]
    total_sum = float(total_weight.sum())
    print(f"  top 15 Pauli strings by total cluster weight "
          f"(of {cluster_size}.0 total):")
    for a in top_alpha:
        w = float(total_weight[a])
        print(f"    {labels[a]}: weight {w:.3f}  "
              f"(XY-wt={xy_weight(labels[a])}, Pi2={pi2_parity(labels[a])}, "
              f"palindromic={spatial_palindrome(labels[a])})")

    # per-mode top Paulis
    print(f"  per cluster mode, top 3 Pauli strings:")
    for i in range(cluster_size):
        ci = np.abs(C[:, i]) ** 2
        top3 = np.argsort(-ci)[:3]
        top_str = "; ".join(
            f"{labels[t]}({100 * ci[t]:.0f}%)" for t in top3)
        print(f"    mode {i}: |Im|={aim[cluster_idx[i]]:.3e}  "
              f"Re={ev[cluster_idx[i]].real:+.4f}  top: {top_str}")

    # sector distribution across cluster modes' dominant Paulis
    xy_c = Counter()
    pi2_c = Counter()
    pal_c = Counter()
    for i in range(cluster_size):
        ci = np.abs(C[:, i]) ** 2
        top = int(np.argmax(ci))
        xy_c[xy_weight(labels[top])] += 1
        pi2_c[pi2_parity(labels[top])] += 1
        pal_c[spatial_palindrome(labels[top])] += 1
    print(f"  cluster sector distribution (by each mode's top Pauli):")
    print(f"    XY-weight: {dict(sorted(xy_c.items()))}")
    print(f"    Pi2-parity: {dict(sorted(pi2_c.items()))}")
    print(f"    palindromic-vs-not: {dict(sorted(pal_c.items()))}")
    print()
    sys.stdout.flush()
