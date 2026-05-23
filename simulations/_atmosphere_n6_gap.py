"""N=6 gap check: when n_osc drops, do the modes jump cleanly to Im=0, or hide
just below the oscillation threshold?

Tom's question on _atmosphere_evenodd.py. For the palindromic gamma-profile at
N=6, several eps: sort |Im lambda|, count modes in the suspicious window
[1e-11, 1e-3], and show the 8 |Im| values bracketing the real/oscillating
boundary. An empty window with a wide bracket gap means the modes genuinely
jump oscillating <-> real; values inside the window mean they hide below the
counting tolerance. Investigation only.
"""
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
N = 6
GAMMA0 = 0.05
J = 0.075


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
d2 = 4 ** N
print(f"N={N}, 4^N={d2}, J={J}, gamma_0={GAMMA0}, palindromic gamma-profile\n")

for eps in (0.0, -0.2, -0.6, -0.8):
    gamma = GAMMA0 * (1.0 + eps * u)
    ev = np.linalg.eigvals(build_L(H, gamma, N))
    aim = np.sort(np.abs(ev.imag))
    n_real = int(np.sum(aim < 1e-7))
    n_osc = d2 - n_real
    window = int(np.sum((aim >= 1e-11) & (aim < 1e-3)))
    lo = max(0, n_real - 4)
    bracket = aim[lo:n_real + 4]
    print(f"eps={eps:+.1f}:  n_osc={n_osc}, n_real={n_real}, "
          f"sum={n_osc + n_real}")
    print(f"   modes with |Im| in [1e-11, 1e-3]: {window}   "
          f"(0 = clean gap, the modes jump)")
    print(f"   |Im| bracketing the real/osc boundary: "
          f"{', '.join(f'{x:.2e}' for x in bracket)}")
    print()
    sys.stdout.flush()
