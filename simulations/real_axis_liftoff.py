"""Track the real-axis modes across the balance: do oscillations get born?

docs/carbon master-question thread (2026-05-22, Tom + Claude).

Tom's correction to the previous artifact: the spectrum's mirror is the
vertical one, conjugation lambda <-> conj(lambda), reflected across the real
axis Im=0. There is NO horizontal palindrome: the gain bath breaks F1 (proven
for Z-dephasing only; the residuals 0.05-0.13 in break_heat_minus_scan.py
already said so). The earlier "best centre" was a fit to a broken palindrome.

In the correct, vertical frame: a pair is a conjugate pair, lambda with
conj(lambda), which together are one oscillation. The real axis Im=0 holds the
non-oscillating modes. "How does a qubit get a baby" becomes: a mode lifts off
the real axis into a conjugate pair, an oscillation born, a heartbeat.

This tracks, as the balance tips (s from 0 to 1; the temperature axis of
break_heat_minus_scan.py is s<->-s symmetric, so one side suffices):
  - the count of real eigenvalues, |Im| < 1e-6, vs s. A drop = modes lifting
    off, oscillations born.
  - Im(lambda) vs s for every eigenvalue (a bifurcation diagram): a fork
    opening from the Im=0 line is a lift-off.

build_L_gain etc. reused from break_heat_minus_scan.py. Investigation only.
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
DOWN = np.array([[0, 1], [0, 0]], dtype=complex)
UP = np.array([[0, 0], [1, 0]], dtype=complex)
PM = {"X": X, "Y": Y, "Z": Z}


def site_op(op, k, N):
    m = np.array([[1.0 + 0j]])
    for i in range(N):
        m = np.kron(m, op if i == k else I2)
    return m


def chain_H(N, combo="XX+YY", J=1.0):
    t1, t2 = combo.split("+")
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        for t in (t1, t2):
            H += J * site_op(PM[t[0]], b, N) @ site_op(PM[t[1]], b + 1, N)
    return H


def build_L_gain(H, gamma_z, gamma_down, gamma_up):
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    n = int(np.log2(d))
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for k in range(n):
        jumps = []
        if gamma_z > 0:
            jumps.append(np.sqrt(gamma_z) * site_op(Z, k, n))
        if gamma_down > 0:
            jumps.append(np.sqrt(gamma_down) * site_op(DOWN, k, n))
        if gamma_up > 0:
            jumps.append(np.sqrt(gamma_up) * site_op(UP, k, n))
        for Lop in jumps:
            LdL = Lop.conj().T @ Lop
            L += np.kron(Lop.conj(), Lop)
            L -= 0.5 * np.kron(Id, LdL)
            L -= 0.5 * np.kron(LdL.T, Id)
    return L


def spectrum(N, s, g_tot=0.1, gz=0.05):
    g_down = g_tot * (1.0 + s) / 2.0
    g_up = g_tot * (1.0 - s) / 2.0
    return np.linalg.eigvals(build_L_gain(chain_H(N), gz, g_down, g_up))


REAL_TOL = 1e-6
s_values = np.linspace(0.0, 1.0, 101)

specs = {3: [], 4: []}
real_count = {3: [], 4: [], 5: []}
for N in (3, 4, 5):
    for s in s_values:
        ev = spectrum(N, s)
        real_count[N].append(int(np.sum(np.abs(ev.imag) < REAL_TOL)))
        if N in specs:
            specs[N].append(ev)

for N in (3, 4, 5):
    rc = real_count[N]
    print(f"N={N} ({4 ** N} eigenvalues): #real  s=0 -> {rc[0]},  "
          f"s=1 -> {rc[-1]},  min {min(rc)},  max {max(rc)}")

OUTDIR = os.path.join("simulations", "results", "real_axis_liftoff")
os.makedirs(OUTDIR, exist_ok=True)

fig, axes = plt.subplots(1, 3, figsize=(19, 5.6))

axA = axes[0]
for N, col in [(3, "tab:green"), (4, "tab:orange"), (5, "tab:purple")]:
    axA.plot(s_values, real_count[N], "-", lw=1.7, color=col, label=f"N={N}")
axA.set_xlabel("s   (0 = balance, 1 = pure decay)")
axA.set_ylabel("# real eigenvalues  (|Im| < 1e-6)")
axA.set_title("Real-mode count across the balance")
axA.legend()
axA.grid(alpha=0.3)

for ax, N in [(axes[1], 3), (axes[2], 4)]:
    for i, s in enumerate(s_values):
        ev = specs[N][i]
        ax.scatter(np.full(len(ev), s), ev.imag, s=5, c="tab:blue",
                   alpha=0.35, linewidths=0)
    ax.axhline(0.0, color="crimson", ls=":", lw=1.2)
    ax.set_ylim(-3.0, 3.0)
    ax.set_xlabel("s")
    ax.set_ylabel("Im lambda")
    ax.set_title(f"N={N}: Im lambda vs s   (fork from Im=0 = lift-off)")
    ax.grid(alpha=0.3)

fig.suptitle("Do modes lift off the real axis when the balance tips?   "
             "s=0 balance -> s=1 pure decay")
fig.tight_layout()
path = os.path.join(OUTDIR, "real_axis_liftoff.png")
fig.savefig(path, dpi=130)
print(f"saved {path}")
