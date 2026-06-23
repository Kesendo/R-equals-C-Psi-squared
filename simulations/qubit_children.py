"""Can a qubit have children? On-site clocks (the source) plus the exchange.

docs/carbon master-question thread (2026-05-22, Tom + Claude).

Tom's question: is the genesis a birth, or are qubits just sources? The genesis
model had H = bonds only, so a lone qubit was inert by construction and any
oscillation looked born. A real qubit has a level splitting, an on-site clock,
it is a source, always rotating.

This puts the on-site term back: H = sum_l h_l Z_l + J*sum_bonds(XX+YY), with
distinct, incommensurate on-site fields h_l = 0.5*sqrt({1,2,3,5}). At J=0 the
qubits oscillate on their own clocks (the source); the all-I/Z Pauli strings
carry no X/Y content, the Z-fields do not rotate them, so they stay real, the
silent modes (2^N of them, no accidental zeros thanks to incommensurate h_l).

The question, "can a qubit have a child": as J turns on, does a silent real
mode lift off into oscillation? A real mode that starts oscillating is a child,
an oscillation that exists only because of the coupling, that no single qubit
had. #real vs J: a drop below the J=0 value of 2^N means children.

Pure F1 Liouvillian (-i[H,.] + Z-dephasing). Investigation only.
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


def onsite_fields(N):
    """Distinct, incommensurate per-qubit on-site fields: the qubits' clocks."""
    return 0.5 * np.sqrt(np.array([1.0, 2.0, 3.0, 5.0])[:N])


def build_H(N, J):
    """H = sum_l h_l Z_l  +  J * sum_bonds (XX + YY)."""
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    h = onsite_fields(N)
    for l in range(N):
        H += h[l] * site_op(Z, l, N)
    for b in range(N - 1):
        for P in (X, Y):
            H += J * site_op(P, b, N) @ site_op(P, b + 1, N)
    return H


def build_L(H, gamma0, N):
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(Z, k, N)
        L += gamma0 * (np.kron(Zk, Zk.conj()) - np.eye(d * d, dtype=complex))
    return L


GAMMA0 = 1.0
IM_TOL = 1e-6
J_values = np.linspace(0.0, 3.0, 301)

n_real = {2: [], 3: [], 4: []}
specs = {2: [], 3: []}
for N in (2, 3, 4):
    for J in J_values:
        ev = np.linalg.eigvals(build_L(build_H(N, J), GAMMA0, N))
        n_real[N].append(int(np.sum(np.abs(ev.imag) < IM_TOL)))
        if N in specs:
            specs[N].append(ev)

for N in (2, 3, 4):
    nr = n_real[N]
    print(f"N={N}: #real at J=0 -> {nr[0]} (2^N = {2 ** N}),  "
          f"at J=3 -> {nr[-1]},  min over scan -> {min(nr)}  "
          f"(children = {nr[0] - min(nr)})")

OUTDIR = os.path.join("simulations", "results", "qubit_children")
os.makedirs(OUTDIR, exist_ok=True)

fig, axes = plt.subplots(1, 3, figsize=(19, 5.6))

axA = axes[0]
for N, col in [(2, "tab:red"), (3, "tab:green"), (4, "tab:blue")]:
    axA.plot(J_values, n_real[N], "-", lw=1.8, color=col, label=f"N={N}")
    axA.axhline(2 ** N, color=col, ls=":", lw=1.0)
axA.set_xlabel("J   (gamma_0 = 1)")
axA.set_ylabel("# real (silent) modes")
axA.set_title("Silent modes vs J   (drop below 2^N = children)")
axA.legend()
axA.grid(alpha=0.3)

for ax, N in [(axes[1], 2), (axes[2], 3)]:
    for i, J in enumerate(J_values):
        ev = specs[N][i]
        ax.scatter(np.full(len(ev), J), ev.imag, s=5, c="tab:purple",
                   alpha=0.35, linewidths=0)
    ax.axhline(0.0, color="crimson", ls=":", lw=1.2)
    ax.set_xlabel("J")
    ax.set_ylabel("Im lambda")
    ax.set_title(f"N={N}: Im lambda vs J   (fork from Im=0 = a child)")
    ax.grid(alpha=0.3)

fig.suptitle("Can a qubit have children?  On-site clocks present, "
             "scanning the coupling J")
fig.tight_layout()
path = os.path.join(OUTDIR, "qubit_children.png")
fig.savefig(path, dpi=130)
print(f"saved {path}")
