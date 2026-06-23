"""N=5 follow-up for atmosphere_mirror_test.py: is the oscillating-mode-count
movement at N=5 real, or tolerance-boundary flicker?

The main test found n_osc flat at N=3,4 but ranging [880,904] at N=5. This
isolates it: n_osc per shape vs eps at four tolerances, plus the |Im lambda|
distribution near zero. A clean gap means n_osc is well-defined and any
movement is real; a smear means the count is tolerance-sensitive flicker.
Investigation only.
"""
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
N = 5
GAMMA0 = 0.05
J = 0.075
TOLS = [1e-5, 1e-6, 1e-7, 1e-8]


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


def shapes(N):
    i = np.arange(N, dtype=float)
    c = (N - 1) / 2.0
    ua = i - c
    ua = ua / np.max(np.abs(ua))
    us = (i - c) ** 2
    us = us - us.mean()
    us = us / np.max(np.abs(us))
    return us, ua


H = J * chain_H(N)
u_sym, u_anti = shapes(N)
eps_values = np.linspace(-0.8, 0.8, 9)

print(f"N={N}, J={J}, gamma_0={GAMMA0}")
print("bins of |Im lambda|:  [ <1e-10 | 1e-10..1e-7 | 1e-7..1e-4 | >1e-4 ]\n")

tol_sensitive = False
for label, u in (("palindromic", u_sym), ("anti-palindromic", u_anti)):
    print(f"--- {label} gamma-profile ---")
    for eps in eps_values:
        gamma = GAMMA0 * (1.0 + eps * u)
        ev = np.linalg.eigvals(build_L(H, gamma, N))
        aim = np.abs(ev.imag)
        counts = [int(np.sum(aim > t)) for t in TOLS]
        if len(set(counts)) > 1:
            tol_sensitive = True
        b = [int(np.sum(aim < 1e-10)),
             int(np.sum((aim >= 1e-10) & (aim < 1e-7))),
             int(np.sum((aim >= 1e-7) & (aim < 1e-4))),
             int(np.sum(aim >= 1e-4))]
        print(f"  eps={eps:+.2f}:  n_osc @tol(1e-5,1e-6,1e-7,1e-8)={counts}"
              f"   |Im| bins={b}")
    print()

print("VERDICT:",
      "tol-sensitive (count depends on the cutoff -> boundary flicker)"
      if tol_sensitive
      else "tol-robust (count independent of the cutoff -> the movement is real)")
