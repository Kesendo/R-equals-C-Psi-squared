"""Why does the n_osc movement appear at N=5? Even/odd, or spectrum density?

Follow-up to _atmosphere_mirror_test.py / _atmosphere_n5_check.py: a middle-
peaked palindromic gamma-profile left the oscillating-mode count n_osc flat at
N=3 (odd) and N=4 (even), and moved it at N=5 (odd). Tom's hypothesis: even/odd,
with N=3 the framework's special minimal case, so N=5 is the first "normal" odd
chain. Competing explanation: spectrum density (N=5 has 1024 modes, N=3 only 64).

N=6 decides:
  even/odd  -> N=6 even  -> n_osc flat
  density   -> N=6 (4096 modes, denser than N=5) -> n_osc moves

This scans the palindromic gamma-profile gamma_i = gamma_0 (1 + eps * u_sym(i))
for N=3,4,5,6 and reports n_osc vs eps. Investigation only.
"""
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
GAMMA0 = 0.05
J = 0.075
IM_TOL = 1e-7
EPS = np.linspace(-0.8, 0.8, 9)
N_LIST = (3, 4, 5, 6)


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
    """Palindromic mean-zero unit shape; minimum at the chain centre."""
    i = np.arange(N, dtype=float)
    u = (i - (N - 1) / 2.0) ** 2
    u = u - u.mean()
    return u / np.max(np.abs(u))


print(f"J={J}, gamma_0={GAMMA0}, palindromic gamma-profile, eps in [-0.8, 0.8]")
print(f"(eps < 0 peaks gamma at the chain centre)\n")
for N in N_LIST:
    H = J * chain_H(N)
    u = sym_shape(N)
    counts = []
    for eps in EPS:
        gamma = GAMMA0 * (1.0 + eps * u)
        ev = np.linalg.eigvals(build_L(H, gamma, N))
        counts.append(int(np.sum(np.abs(ev.imag) > IM_TOL)))
    parity = "odd" if N % 2 else "even"
    moved = max(counts) - min(counts)
    print(f"N={N} ({parity}, {4 ** N} modes):  n_osc = {counts}")
    print(f"      uniform (eps=0) = {counts[len(EPS) // 2]},   "
          f"range = [{min(counts)}, {max(counts)}],   moved by {moved}")
    print()
    sys.stdout.flush()

print("Tom's even/odd hypothesis predicts N=6 flat (moved by 0);")
print("the density explanation predicts N=6 moves.")
