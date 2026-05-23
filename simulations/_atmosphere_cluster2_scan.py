"""Second-cluster fine scan: N=6 around the eps ~ -0.22 dip (28 pairs).

The staircase (_atmosphere_staircase.py) showed N=6 has multiple n_osc dip
events along eps in [-1, 0]. The deepest at eps <= -0.85 is the 16-pair
cluster we characterized. The second-deepest is the 28-pair dip at
eps in [-0.25, -0.20]. This is its fine-scan companion to _atmosphere_
cluster_fine.py: track the smallest non-real |Im| (with proper sort!)
across eps in [-0.30, -0.15], to locate the second cluster's coalescence
ep* and its dip depth. Sets up the second mode-ID step.
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
N = 6
GAMMA0 = 0.05
J = 0.075
EPS = np.linspace(-0.30, -0.15, 16)


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
print(f"N={N}, eps in [-0.30, -0.15], Delta-eps ~ 0.01")
print("locating the second-cluster (~28-pair) dip\n")
rows = []
for eps in EPS:
    gamma = GAMMA0 * (1.0 + eps * u)
    ev = np.linalg.eigvals(build_L(H, gamma, N))
    aim = np.sort(np.abs(ev.imag))
    n_real = int(np.sum(aim < 1e-10))
    nonreal = aim[aim >= 1e-13]
    smin = float(nonreal[0]) if nonreal.size else float('nan')
    s6 = nonreal[:6]
    rows.append((float(eps), n_real, smin))
    print(f"  eps={eps:+.4f}: n_real={n_real}  smallest |Im|={smin:.3e}   "
          f"6 smallest: {', '.join(f'{x:.2e}' for x in s6)}")
    sys.stdout.flush()

ims = [r[2] for r in rows]
imin = int(np.argmin(ims))
print(f"\n--> minimum smallest-|Im| = {ims[imin]:.3e} at eps={rows[imin][0]:+.4f}")

OUTDIR = os.path.join("simulations", "results", "atmosphere_cluster2_scan")
os.makedirs(OUTDIR, exist_ok=True)
fig, ax = plt.subplots(figsize=(8, 5))
eps = [r[0] for r in rows]
ax.semilogy(eps, ims, "-o", ms=4, color="tab:red")
ax.set_xlabel("eps")
ax.set_ylabel("smallest non-real |Im|")
ax.set_title("N=6: second-cluster trace, eps in [-0.30, -0.15]")
ax.grid(alpha=0.3, which="both")
fig.tight_layout()
path = os.path.join(OUTDIR, "atmosphere_cluster2_scan.png")
fig.savefig(path, dpi=130)
print(f"saved {path}")
