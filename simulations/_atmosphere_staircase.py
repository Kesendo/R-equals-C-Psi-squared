"""Mapping the n_osc staircase over the middle-peaked gamma window.

Following the cluster thread (#17). The original even/odd scan showed N=6
n_osc=[3876,3894,3892,3860,3916,...] for eps=[-0.8,-0.6,-0.4,-0.2,0,...]:
multiple dips, not just one. The cluster we have characterized lives at the
eps=-0.83 dip. This scans eps in [-1, 0] at Delta-eps=0.05 for N=5 and N=6,
tracking n_osc, n_real, and the smallest non-real |Im|, to see whether the
staircase has multiple cluster-dip events at different eps - candidates for
additional clusters in other XY-weight sectors.
Investigation only.
"""
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
GAMMA0 = 0.05
J = 0.075
EPS = np.linspace(-1.0, 0.0, 21)
N_LIST = (5, 6)


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


results = {}
for N in N_LIST:
    H = J * chain_H(N)
    u = sym_shape(N)
    rows = []
    print(f"=== N={N} ({4 ** N} modes) ===")
    print(f"  {'eps':>7}  {'n_real':>8}  {'n_osc':>8}  "
          f"{'smallest non-real |Im|':>22}")
    for eps in EPS:
        gamma = GAMMA0 * (1.0 + eps * u)
        ev = np.linalg.eigvals(build_L(H, gamma, N))
        aim = np.abs(ev.imag)
        n_real_10 = int(np.sum(aim < 1e-10))
        n_osc = int(np.sum(aim > 1e-7))
        nonreal = aim[aim >= 1e-13]
        smin = float(nonreal[0]) if nonreal.size else float('nan')
        rows.append((float(eps), n_real_10, n_osc, smin))
        print(f"  {eps:>+7.3f}  {n_real_10:>8}  {n_osc:>8}  {smin:>22.3e}")
        sys.stdout.flush()
    results[N] = rows
    print()

OUTDIR = os.path.join("simulations", "results", "atmosphere_staircase")
os.makedirs(OUTDIR, exist_ok=True)
fig, axes = plt.subplots(2, 2, figsize=(13, 9))
for col, N in enumerate(N_LIST):
    rows = results[N]
    eps = [r[0] for r in rows]
    n_osc = [r[2] for r in rows]
    smin = [r[3] for r in rows]
    axes[0][col].plot(eps, n_osc, "-o", ms=4, color="tab:purple")
    axes[0][col].set_xlabel("eps")
    axes[0][col].set_ylabel("n_osc")
    axes[0][col].set_title(f"N={N}: n_osc vs eps (staircase)")
    axes[0][col].grid(alpha=0.3)
    axes[1][col].semilogy(eps, smin, "-o", ms=4, color="tab:red")
    axes[1][col].set_xlabel("eps")
    axes[1][col].set_ylabel("smallest non-real |Im|")
    axes[1][col].set_title(f"N={N}: cluster |Im| vs eps")
    axes[1][col].grid(alpha=0.3, which="both")
fig.suptitle("Staircase mapping: multiple cluster dips along the "
             "middle-peaked gamma direction")
fig.tight_layout()
path = os.path.join(OUTDIR, "atmosphere_staircase.png")
fig.savefig(path, dpi=130)
print(f"saved {path}")
