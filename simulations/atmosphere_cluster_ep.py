"""EP confirmation: is the cluster's near-coalescence a true exceptional point?

The fine eps-scan (atmosphere_cluster_fine.py) showed the degenerate cluster's
|Im| swept to the real axis: N=5 coalesces fully (|Im| < 1e-13 near eps=-0.33),
N=6 dips to 7e-5 near eps=-0.83. A conjugate pair meeting the real axis is
generically an exceptional point. This confirms it directly.

At a true EP the eigenvectors coalesce too, so the right-eigenvector matrix V
becomes singular (cond(V) diverges) and the cluster's own 32 eigenvectors lose
rank. Full eigendecomposition of L at N=6 for eps approaching the dip and away
from it; if cond(V) and the cluster-eigenvector rank-deficiency grow toward
eps=-0.83, the near-coalescence is an EP. Investigation only.
"""
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
N = 6
GAMMA0 = 0.05
J = 0.075
EPS = [-0.60, -0.75, -0.83, -0.85]


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


H = J * chain_H(N)
u = sym_shape(N)
print(f"N={N}, EP confirmation via eigenvector coalescence")
print("a true EP -> cond(V) diverges and the 32 cluster eigenvectors "
      "lose rank\n")
for eps in EPS:
    gamma = GAMMA0 * (1.0 + eps * u)
    ev, V = np.linalg.eig(build_L(H, gamma, N))
    aim = np.abs(ev.imag)
    order = np.argsort(aim)
    nonreal = [i for i in order if aim[i] >= 1e-13]
    cluster_idx = nonreal[:32]
    cluster_im = aim[cluster_idx[0]]
    condV = float(np.linalg.cond(V))
    Vcl = V[:, cluster_idx].copy()
    Vcl /= np.linalg.norm(Vcl, axis=0, keepdims=True)
    s = np.linalg.svd(Vcl, compute_uv=False)
    rank = int(np.sum(s > 1e-8 * s[0]))
    print(f"eps={eps:+.2f}:  cluster |Im| = {cluster_im:.3e}   "
          f"cond(V) = {condV:.3e}")
    print(f"            cluster-eigenvector singular values "
          f"{s[0]:.2e} .. {s[-1]:.2e},  rank = {rank}/32")
    print()
    sys.stdout.flush()
