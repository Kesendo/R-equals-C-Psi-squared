"""Inspect the real spectrum artifact at the balance s=0.

docs/carbon master-question thread (2026-05-22, Tom + Claude).

Tom's correction: the genesis is not "hidden", it is "inconspicuous", so plain
that nobody looks. The paired-count of _break_heat_minus_scan.py was waved away
as "erratic, threshold-sensitive". That is the trap: a scalar proxy looks like
noise while the real artifact carries the structure. A "baby" is a new mirror
pair, and the pairing got reduced to a wobbly number and skipped.

So: stop counting, inspect the artifact. This plots the Liouvillian eigenvalue
spectrum in the complex plane at s=0 (the balance, where modes are degenerate /
merged) and just off it (s=0.05, 0.20), where the merged modes split. The
palindrome centre (best-fit) is marked; a mirror pair is two eigenvalues
point-symmetric through (centre, 0). The question, visible in the picture: when
a degenerate s=0 mode splits, are the two children a mirror pair? One mode -> a
pair, one becomes two: a qubit getting a baby.

build_L_gain etc. reused from _break_heat_minus_scan.py. Investigation only.
"""
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

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


def residuals(ev, centre):
    reflected = 2.0 * complex(centre) - ev
    return np.abs(ev[None, :] - reflected[:, None]).min(axis=1)


def spectrum(N, s, g_tot=0.1, gz=0.05):
    g_down = g_tot * (1.0 + s) / 2.0
    g_up = g_tot * (1.0 - s) / 2.0
    return np.linalg.eigvals(build_L_gain(chain_H(N), gz, g_down, g_up))


def best_centre(ev):
    lo, hi = ev.real.min() - 1.0, ev.real.max() + 1.0
    opt = minimize_scalar(lambda c: residuals(ev, c).max(),
                          bounds=(lo, hi), method="bounded")
    return opt.x


def n_distinct(ev, dec=4):
    return len({(round(e.real, dec), round(e.imag, dec)) for e in ev})


OUTDIR = os.path.join("simulations", "results", "s0_spectrum_artifact")
os.makedirs(OUTDIR, exist_ok=True)

S_SHOW = [(0.00, "s=0  balance", "black", "o", 28),
          (0.05, "s=0.05", "tab:red", "x", 34),
          (0.20, "s=0.20", "tab:blue", "+", 40)]

fig, axes = plt.subplots(1, 2, figsize=(14, 8))
for ax, N in zip(axes, (3, 4)):
    specs = {s: spectrum(N, s) for s, *_ in S_SHOW}
    for s, lbl, col, mk, sz in S_SHOW:
        ev = specs[s]
        ax.scatter(ev.real, ev.imag, s=sz, c=col, marker=mk, alpha=0.55,
                   linewidths=1.0, label=f"{lbl}  ({n_distinct(ev)} distinct)")
    c0 = best_centre(specs[0.00])
    ax.axvline(c0, color="gray", ls=":", lw=1.2)
    ax.axhline(0.0, color="gray", ls=":", lw=0.8)
    ax.scatter([c0], [0.0], c="gray", marker="*", s=170, zorder=6,
               label=f"centre {c0:.3f}")
    ax.set_title(f"N = {N}   Liouvillian spectrum: balance and just off it")
    ax.set_xlabel("Re lambda")
    ax.set_ylabel("Im lambda")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")
fig.suptitle("Spectrum artifact at the balance s=0: "
             "do merged modes split into mirror pairs?")
fig.tight_layout()
path = os.path.join(OUTDIR, "s0_spectrum.png")
fig.savefig(path, dpi=130)
print(f"saved {path}")

for N in (3, 4, 5):
    cells = "  ".join(f"s={s:.2f}: {n_distinct(spectrum(N, s)):>4} distinct"
                      for s, *_ in S_SHOW)
    print(f"N={N}  ({4 ** N} eigenvalues)   {cells}")
