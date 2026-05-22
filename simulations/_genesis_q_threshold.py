"""The genesis N=1 -> N=2, gated by Q = J/gamma_0: the first oscillation.

docs/carbon master-question thread (2026-05-22, Tom + Claude).

The heat axis is closed: _real_axis_liftoff.py showed the bath moves only the
real parts; Im (the oscillations) is rigid against it. Bath and Hamiltonian are
orthogonal. What remains is gamma_0 and J, that is Q = J/gamma_0, the framework
axis.

The genesis "how does a qubit get a baby" is then the step N=1 -> N=2 (Tom's
early reframing). N=1: one qubit, no bond, H=0, only Z-dephasing -> every
eigenvalue real, no oscillation, no heartbeat. N=2: the first J-bond couples
two qubits, a mode CAN oscillate. Whether it does is overdamped vs underdamped:
gamma_0 too strong -> overdamped, the eigenvalue stays real, no heartbeat; J
strong enough -> underdamped, the eigenvalue lifts off the real axis, the first
heartbeat.

The earlier lift-off hunt (_real_axis_liftoff.py) was on the wrong axis, the
bath temperature. This is the right one: Q. This scans Q for N=2,3,4 and finds
the lift-off threshold Q*, the genesis. Pure F1 system: H + Z-dephasing, no
amplitude channel. Investigation only.
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


def site_op(op, k, N):
    m = np.array([[1.0 + 0j]])
    for i in range(N):
        m = np.kron(m, op if i == k else I2)
    return m


def chain_H(N):
    """XX + YY on every bond, unit coupling."""
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        for P in (X, Y):
            H += site_op(P, b, N) @ site_op(P, b + 1, N)
    return H


def build_L(H, gamma0, N):
    """-i[H,.] + Z-dephasing at rate gamma0 per site (pure F1 system)."""
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(Z, k, N)
        L += gamma0 * (np.kron(Zk, Zk.conj()) - np.eye(d * d, dtype=complex))
    return L


def spectrum(N, Q, gamma0=1.0):
    H = (Q * gamma0) * chain_H(N)
    return np.linalg.eigvals(build_L(H, gamma0, N))


IM_TOL = 1e-6
Q_values = np.linspace(0.0, 3.0, 301)

n_complex = {2: [], 3: [], 4: []}
specs = {2: [], 3: []}
for N in (2, 3, 4):
    for Q in Q_values:
        ev = spectrum(N, Q)
        n_complex[N].append(int(np.sum(np.abs(ev.imag) > IM_TOL)))
        if N in specs:
            specs[N].append(ev)

qstars = {}
for N in (2, 3, 4):
    nc = n_complex[N]
    idx = next((i for i, c in enumerate(nc) if c > 0), None)
    qstars[N] = Q_values[idx] if idx is not None else None
    print(f"N={N}: Q* (first lift-off) = {qstars[N]},  "
          f"#complex at Q=3 -> {nc[-1]} of {4 ** N}")

for N in (2, 3, 4):
    nc = n_complex[N]
    steps = [(Q_values[i], nc[i - 1], nc[i])
             for i in range(1, len(nc)) if nc[i] != nc[i - 1]]
    print(f"N={N} staircase, {len(steps)} steps:")
    for q, a, b in steps:
        print(f"  Q={q:.2f}:  {a} -> {b}  (+{b - a})")

OUTDIR = os.path.join("simulations", "results", "genesis_q_threshold")
os.makedirs(OUTDIR, exist_ok=True)

fig, axes = plt.subplots(1, 3, figsize=(19, 5.6))

axA = axes[0]
for N, col in [(2, "tab:red"), (3, "tab:green"), (4, "tab:blue")]:
    axA.plot(Q_values, n_complex[N], "-", lw=1.8, color=col, label=f"N={N}")
    if qstars[N] is not None:
        axA.axvline(qstars[N], color=col, ls=":", lw=1.1)
axA.axvline(1.0, color="black", ls="--", lw=0.9, label="Q=1 (balance)")
axA.set_xlabel("Q = J / gamma_0")
axA.set_ylabel("# complex eigenvalues  (|Im| > 1e-6)")
axA.set_title("Oscillations born vs Q   (lift-off = first heartbeat)")
axA.legend()
axA.grid(alpha=0.3)

for ax, N in [(axes[1], 2), (axes[2], 3)]:
    for i, Q in enumerate(Q_values):
        ev = specs[N][i]
        ax.scatter(np.full(len(ev), Q), ev.imag, s=5, c="tab:purple",
                   alpha=0.35, linewidths=0)
    ax.axhline(0.0, color="crimson", ls=":", lw=1.2)
    if qstars[N] is not None:
        ax.axvline(qstars[N], color="gray", ls="--", lw=1.1)
    ax.set_xlabel("Q = J / gamma_0")
    ax.set_ylabel("Im lambda")
    ax.set_title(f"N={N}: Im lambda vs Q   (fork from Im=0 = oscillation born)")
    ax.grid(alpha=0.3)

fig.suptitle("The genesis on the Q axis: where does the first oscillation "
             "lift off the real axis?")
fig.tight_layout()
path = os.path.join(OUTDIR, "genesis_q_threshold.png")
fig.savefig(path, dpi=130)
print(f"saved {path}")
