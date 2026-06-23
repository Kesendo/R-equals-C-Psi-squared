"""The genesis at the framework's normal values: gamma_0 = 0.05, J around 0.075.

docs/carbon master-question thread (2026-05-22, Tom + Claude).

genesis_q_threshold.py found (gamma_0 = 1): the first oscillations are born
unconditionally at Q = 0+ (the instant J leaves zero), then a small staircase
of discrete later lift-offs (N=2 at Q~0.51, N=3 at Q~0.71, N=4 at Q~0.45,
0.48, 0.94).

Tom: gamma_0 is the fixed constant, J the only knob, Q = J/gamma_0 the scale.
The normal values: gamma_0 = 0.05, J = 0.075, that is Q = 1.5, the framework's
peak anchor.

This re-runs at gamma_0 = 0.05. The Liouvillian factors exactly as
L = gamma_0 * Ltilde(Q), so every eigenvalue scales by gamma_0 and the
Q-structure (the staircase, the birth) is gamma_0-independent. The run shows
the birth at the real absolute scale and confirms Q is the scale. The operating
point Q = 1.5 (J = 0.075) is marked. Investigation only.
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
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        for P in (X, Y):
            H += site_op(P, b, N) @ site_op(P, b + 1, N)
    return H


def build_L(H, gamma0, N):
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(Z, k, N)
        L += gamma0 * (np.kron(Zk, Zk.conj()) - np.eye(d * d, dtype=complex))
    return L


GAMMA0 = 0.05
J_OP = 0.075
Q_OP = J_OP / GAMMA0          # 1.5
IM_TOL = 1e-6 * GAMMA0        # complex-detection tol, scaled with gamma_0


def spectrum(N, Q):
    H = (Q * GAMMA0) * chain_H(N)
    return np.linalg.eigvals(build_L(H, GAMMA0, N))


Q_values = np.linspace(0.0, 4.0, 401)

n_complex = {2: [], 3: [], 4: []}
specs = {2: [], 3: []}
for N in (2, 3, 4):
    for Q in Q_values:
        ev = spectrum(N, Q)
        n_complex[N].append(int(np.sum(np.abs(ev.imag) > IM_TOL)))
        if N in specs:
            specs[N].append(ev)

for N in (2, 3, 4):
    nc = n_complex[N]
    steps = [(Q_values[i], nc[i - 1], nc[i])
             for i in range(1, len(nc)) if nc[i] != nc[i - 1]]
    print(f"N={N} staircase, {len(steps)} steps:")
    for q, a, b in steps:
        print(f"  Q={q:.2f} (J={q * GAMMA0:.4f}):  {a} -> {b}  (+{b - a})")

print(f"\nOperating point: gamma_0 = {GAMMA0}, J = {J_OP}  (Q = {Q_OP})")
for N in (2, 3, 4):
    ev = spectrum(N, Q_OP)
    osc = np.abs(ev.imag)[np.abs(ev.imag) > IM_TOL]
    print(f"  N={N}: {len(osc)} of {4 ** N} modes oscillate;  "
          f"frequency |Im| = {osc.min():.4f} .. {osc.max():.4f}")

OUTDIR = os.path.join("simulations", "results", "genesis_real_gamma")
os.makedirs(OUTDIR, exist_ok=True)

fig, axes = plt.subplots(1, 3, figsize=(19, 5.6))

axA = axes[0]
for N, col in [(2, "tab:red"), (3, "tab:green"), (4, "tab:blue")]:
    axA.plot(Q_values, n_complex[N], "-", lw=1.8, color=col, label=f"N={N}")
axA.axvline(Q_OP, color="black", ls="--", lw=1.3, label=f"J={J_OP}  (Q={Q_OP})")
axA.set_xlabel("Q = J / gamma_0")
axA.set_ylabel("# oscillating modes")
axA.set_title(f"Oscillations born vs Q   (gamma_0 = {GAMMA0})")
axA.legend()
axA.grid(alpha=0.3)

for ax, N in [(axes[1], 2), (axes[2], 3)]:
    for i, Q in enumerate(Q_values):
        ev = specs[N][i]
        ax.scatter(np.full(len(ev), Q), ev.imag, s=5, c="tab:purple",
                   alpha=0.35, linewidths=0)
    ax.axhline(0.0, color="crimson", ls=":", lw=1.2)
    ax.axvline(Q_OP, color="black", ls="--", lw=1.3)
    ax.set_xlabel("Q = J / gamma_0")
    ax.set_ylabel("Im lambda")
    ax.set_title(f"N={N}: Im lambda vs Q   (gamma_0 = {GAMMA0}, the birth)")
    ax.grid(alpha=0.3)

fig.suptitle(f"The genesis at the normal values: gamma_0 = {GAMMA0}, "
             f"operating point J = {J_OP}  (Q = {Q_OP})")
fig.tight_layout()
path = os.path.join(OUTDIR, "genesis_real_gamma.png")
fig.savefig(path, dpi=130)
print(f"saved {path}")
