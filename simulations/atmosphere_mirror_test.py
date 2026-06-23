"""The atmosphere test: the F71-twin for gamma.

docs/ atmosphere thread (2026-05-22, Tom + Claude).

THE_ATMOSPHERE_AND_THE_CANCELLED_FORMULAS names "atmosphere" = gamma_0, the
uniform Z-dephasing. The message: only the right atmosphere works. F71 found,
for non-uniform J: the framework's spatial mirror survives a palindromic
J-profile and deviates only under the anti-palindromic part. This applies the
same decomposition to the atmosphere gamma itself.

A gamma-profile gamma_i = gamma_0 * (1 + eps * u_i) is split by shape:
  u_sym  palindromic       u(i) =  u(N-1-i)
  u_anti anti-palindromic  u(i) = -u(N-1-i)
Both shapes are mean-zero, so Sum gamma_i = N*gamma_0 stays fixed; only the
profile shape changes. For each shape, eps is scanned and three quantities
are measured on the full Liouvillian L:

  d_spatial : ||[L, R]|| / ||L||, R the spatial-reflection superoperator.
              The F71 spatial mirror. Hypothesis: 0 for palindromic gamma,
              growing for anti-palindromic.
  d_F1      : the spectral palindrome residual (lambda -> -lambda - 2 sigma
              matching error). The F1 mirror. Control: expect ~0 for any
              per-site gamma (Master Lemma).
  n_osc     : count of oscillating modes (|Im lambda| > tol). Genesis
              control: gamma touches Re, so expect it gamma-shape-independent.

Pure F1 system: H = uniform XY chain, per-site Z-dephasing. Investigation only.
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
J = 0.075                       # uniform coupling; Q = J/gamma_0 = 1.5
IM_TOL = 1e-7
EPS_VALUES = np.linspace(-0.8, 0.8, 33)
N_LIST = (3, 4, 5)


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


def build_L(H, gamma, N):
    """-i[H,.] + per-site Z-dephasing at rates gamma[k]  (row-stacking vec)."""
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(Z, k, N)
        L += gamma[k] * (np.kron(Zk, Zk.conj()) - np.eye(d * d, dtype=complex))
    return L


def reflection_super(N):
    """Spatial reflection: reverse qubit order, lifted to Liouville space."""
    d = 2 ** N
    P = np.zeros((d, d), dtype=complex)
    for x in range(d):
        r = 0
        for i in range(N):
            r |= ((x >> i) & 1) << (N - 1 - i)
        P[r, x] = 1.0
    return np.kron(P, P)


def shapes(N):
    """Palindromic and anti-palindromic mean-zero unit shapes."""
    i = np.arange(N, dtype=float)
    c = (N - 1) / 2.0
    u_anti = i - c
    u_anti = u_anti / np.max(np.abs(u_anti))
    u_sym = (i - c) ** 2
    u_sym = u_sym - u_sym.mean()
    u_sym = u_sym / np.max(np.abs(u_sym))
    return u_sym, u_anti


def palindrome_residual(ev, sigma):
    """Max matching error of the spectrum under lambda -> -lambda - 2 sigma."""
    partners = -ev - 2.0 * sigma
    diff = np.abs(ev[None, :] - partners[:, None])
    return float(np.max(np.min(diff, axis=1)))


results = {}
for N in N_LIST:
    H = J * chain_H(N)
    Rhat = reflection_super(N)
    u_sym, u_anti = shapes(N)
    d_spatial = {"sym": [], "anti": []}
    d_f1 = {"sym": [], "anti": []}
    n_osc = {"sym": [], "anti": []}
    for label, u in (("sym", u_sym), ("anti", u_anti)):
        for eps in EPS_VALUES:
            gamma = GAMMA0 * (1.0 + eps * u)
            L = build_L(H, gamma, N)
            comm = L @ Rhat - Rhat @ L
            d_spatial[label].append(
                float(np.linalg.norm(comm) / np.linalg.norm(L)))
            ev = np.linalg.eigvals(L)
            d_f1[label].append(palindrome_residual(ev, float(gamma.sum())))
            n_osc[label].append(int(np.sum(np.abs(ev.imag) > IM_TOL)))
    results[N] = dict(d_spatial=d_spatial, d_f1=d_f1, n_osc=n_osc)

    f1max = max(max(d_f1["sym"]), max(d_f1["anti"]))
    osc_all = n_osc["sym"] + n_osc["anti"]
    print(f"N={N}:")
    print(f"  F1 palindrome residual, max over scan = {f1max:.2e}   "
          f"(control, expect ~0)")
    print(f"  n_osc, range over scan = [{min(osc_all)}, {max(osc_all)}]  "
          f"of {4 ** N}   (control, expect flat)")
    print(f"  d_spatial, palindromic gamma:      "
          f"[{min(d_spatial['sym']):.2e}, {max(d_spatial['sym']):.2e}]")
    print(f"  d_spatial, anti-palindromic gamma: "
          f"[{min(d_spatial['anti']):.2e}, {max(d_spatial['anti']):.2e}]")
    print()

OUTDIR = os.path.join("simulations", "results", "atmosphere_mirror_test")
os.makedirs(OUTDIR, exist_ok=True)

fig, axes = plt.subplots(1, 3, figsize=(18, 5.4))
for ax, N in zip(axes, N_LIST):
    ds = results[N]["d_spatial"]
    ax.plot(EPS_VALUES, ds["sym"], "-o", ms=3, color="tab:green",
            label="palindromic gamma-profile")
    ax.plot(EPS_VALUES, ds["anti"], "-o", ms=3, color="tab:red",
            label="anti-palindromic gamma-profile")
    ax.axhline(0.0, color="gray", ls=":", lw=1)
    ax.set_xlabel("eps   (gamma-profile amplitude)")
    ax.set_ylabel("spatial-mirror deviation  ||[L,R]|| / ||L||")
    ax.set_title(f"N={N}")
    ax.legend()
    ax.grid(alpha=0.3)
fig.suptitle("The atmosphere test: the F71 spatial mirror under a "
             "non-uniform gamma-profile   (uniform J, Q=1.5)")
fig.tight_layout()
path = os.path.join(OUTDIR, "atmosphere_mirror_test.png")
fig.savefig(path, dpi=130)
print(f"saved {path}")
