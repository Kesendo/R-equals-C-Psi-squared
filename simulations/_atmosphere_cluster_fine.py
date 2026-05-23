"""Resolving the cluster |Im| minimum: a true EP or an avoided coalescence?

Follows _atmosphere_cluster_scan.py. That found a persistent, sharply degenerate
cluster (N=6: 16 mode-pairs; N=5: 12 pairs) whose |Im| the middle-peaked
palindromic gamma-profile sweeps, dipping to a sharp minimum: N=6 near eps=-0.85
(|Im| ~ 2.35e-4), N=5 near eps=-0.35 (|Im| ~ 1e-3). At Delta-eps=0.05 it did not
reach 0. This refines each dip with Delta-eps ~ 0.01: does the cluster |Im| reach
zero (eigenvalue coalescence, an exceptional point) or bottom at a nonzero value
(an avoided coalescence)? Investigation only.
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
WINDOWS = {5: np.linspace(-0.45, -0.25, 15),
           6: np.linspace(-0.92, -0.78, 15)}


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
for N in (5, 6):
    H = J * chain_H(N)
    u = sym_shape(N)
    rows = []
    print(f"=== N={N} ===")
    for eps in WINDOWS[N]:
        gamma = GAMMA0 * (1.0 + eps * u)
        ev = np.linalg.eigvals(build_L(H, gamma, N))
        aim = np.sort(np.abs(ev.imag))
        n_real = int(np.sum(aim < 1e-10))
        nonreal = aim[aim >= 1e-13]
        cluster_im = float(nonreal[0]) if nonreal.size else float('nan')
        s6 = nonreal[:6]
        rows.append((float(eps), n_real, cluster_im))
        print(f"  eps={eps:+.4f}: n_real={n_real}  cluster |Im| = "
              f"{cluster_im:.3e}   smallest6: "
              f"{', '.join(f'{x:.2e}' for x in s6)}")
        sys.stdout.flush()
    ims = [r[2] for r in rows]
    imin = int(np.argmin(ims))
    print(f"  --> minimum cluster |Im| = {ims[imin]:.3e} "
          f"at eps={rows[imin][0]:+.4f}")
    print()
    results[N] = rows

OUTDIR = os.path.join("simulations", "results", "atmosphere_cluster_fine")
os.makedirs(OUTDIR, exist_ok=True)
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
for ax, N in zip(axes, (5, 6)):
    rows = results[N]
    eps = [r[0] for r in rows]
    ims = [r[2] for r in rows]
    ax.semilogy(eps, ims, "-o", ms=4, color="tab:red")
    ax.set_xlabel("eps")
    ax.set_ylabel("cluster |Im|")
    ax.set_title(f"N={N}: cluster |Im| near its dip")
    ax.grid(alpha=0.3, which="both")
fig.suptitle("Resolving the cluster |Im| minimum: "
             "EP (-> 0) or avoided coalescence")
fig.tight_layout()
path = os.path.join(OUTDIR, "atmosphere_cluster_fine.png")
fig.savefig(path, dpi=130)
print(f"saved {path}")
