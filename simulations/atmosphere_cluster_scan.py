"""Following the cluster thread: the degenerate near-EP cluster.

atmosphere_n6_gap.py found, at the extreme middle-peaked palindromic gamma
(eps=-0.8, N=6), 32 modes (16 conjugate pairs) at |Im lambda| ~ 5.6e-4, a sharp
degenerate cluster, neither genuinely real nor a genuine oscillation. It was
absent at eps=-0.6. This traces the cluster: a fine eps-scan for N=5 and N=6,
binning |Im lambda| per eps. It shows where the cluster is born and whether its
|Im| heads toward zero (an exceptional point) as eps grows more negative.

eps is bounded below by -1: the chain-end gamma reaches 0 there, below it the
ends would have negative (gain) dephasing. Investigation only.
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
N_LIST = (5, 6)
EPS = np.linspace(-1.0, -0.3, 15)


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
    print("  bins of |Im|: [ <1e-10 | 1e-10..1e-4 | 1e-4..1e-3 | "
          "1e-3..1e-2 | >=1e-2 ]")
    for eps in EPS:
        gamma = GAMMA0 * (1.0 + eps * u)
        ev = np.linalg.eigvals(build_L(H, gamma, N))
        aim = np.abs(ev.imag)
        b0 = int(np.sum(aim < 1e-10))
        b1 = int(np.sum((aim >= 1e-10) & (aim < 1e-4)))
        b2 = int(np.sum((aim >= 1e-4) & (aim < 1e-3)))
        b3 = int(np.sum((aim >= 1e-3) & (aim < 1e-2)))
        b4 = int(np.sum(aim >= 1e-2))
        n_low = b1 + b2                          # the flicker zone [1e-10,1e-3)
        nonreal = np.sort(aim[aim >= 1e-10])
        smin = float(nonreal[0]) if nonreal.size else float('nan')
        smallest3 = nonreal[:3]
        rows.append((eps, b0, n_low, smin))
        print(f"  eps={eps:+.2f}: bins=[{b0},{b1},{b2},{b3},{b4}]  "
              f"flicker-zone(n_low)={n_low}  "
              f"smallest non-real |Im| = "
              f"{', '.join(f'{x:.2e}' for x in smallest3)}")
        sys.stdout.flush()
    results[N] = rows
    print()

OUTDIR = os.path.join("simulations", "results", "atmosphere_cluster_scan")
os.makedirs(OUTDIR, exist_ok=True)
fig, axes = plt.subplots(2, 2, figsize=(13, 9))
for col, N in enumerate(N_LIST):
    rows = results[N]
    eps = [r[0] for r in rows]
    n_low = [r[2] for r in rows]
    smin = [r[3] for r in rows]
    axes[0][col].plot(eps, n_low, "-o", ms=4, color="tab:purple")
    axes[0][col].set_xlabel("eps")
    axes[0][col].set_ylabel("flicker-zone size  (modes in [1e-10, 1e-3))")
    axes[0][col].set_title(f"N={N}: cluster size vs eps")
    axes[0][col].grid(alpha=0.3)
    axes[1][col].semilogy(eps, smin, "-o", ms=4, color="tab:red")
    axes[1][col].set_xlabel("eps")
    axes[1][col].set_ylabel("smallest non-real |Im|")
    axes[1][col].set_title(f"N={N}: lowest non-real mode (toward 0 = an EP)")
    axes[1][col].grid(alpha=0.3)
fig.suptitle("Following the cluster thread: the near-EP flicker cluster "
             "under a middle-peaked gamma-profile")
fig.tight_layout()
path = os.path.join(OUTDIR, "atmosphere_cluster_scan.png")
fig.savefig(path, dpi=130)
print(f"saved {path}")
