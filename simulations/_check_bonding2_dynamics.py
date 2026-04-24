#!/usr/bin/env python3
"""F75 verification: MM(t) dynamics for bonding:2 at N=5 under full Heisenberg + dephasing.

Confirms that MM(0) = 2 h(1/4) - h(1/2) + 2 h(1/4) - h(1/2) = 1.2451 is reached
as the Lindblad evolution max: MM oscillates and its revival at t approx pi/Delta
returns to just above MM(0), with dephasing then damping it.

At N=5 bonding:2, Delta = 2 sqrt(5) (odd-sector Heisenberg gap), so revival
sits at t approx 0.7; this explains the C# brecher observed PeakT_MM = 0.60
(coarse measurement grid rounds to nearest 0.1).
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from scipy.linalg import expm

N = 5
d = 2 ** N
gamma = 0.05

I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)


def site_op(op, k):
    ops = [I2] * N
    ops[k] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


H = np.zeros((d, d), dtype=complex)
for b in range(N - 1):
    H += (site_op(sx, b) @ site_op(sx, b + 1)
          + site_op(sy, b) @ site_op(sy, b + 1)
          + site_op(sz, b) @ site_op(sz, b + 1))

Id = np.eye(d, dtype=complex)
L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
for k in range(N):
    Zk = site_op(sz, k)
    L += gamma * (np.kron(Zk, Zk.conj()) - np.eye(d * d, dtype=complex))

norm = np.sqrt(2.0 / (N + 1))
psi = np.zeros(d, dtype=complex)
for j in range(N):
    amp = norm * np.sin(np.pi * 2 * (j + 1) / (N + 1))
    idx = 1 << (N - 1 - j)
    psi[idx] = amp
rho0 = np.outer(psi, psi.conj())
rho0 = (rho0 + rho0.conj().T) / 2
rho0 /= np.trace(rho0).real


def evolve_vec(rho_vec, t):
    return expm(L * t) @ rho_vec


def ptrace_keep(rho, keep):
    keep = list(keep)
    trace_out = [q for q in range(N) if q not in keep]
    dims = [2] * N
    r = rho.reshape(dims + dims)
    cn = N
    for q in sorted(trace_out, reverse=True):
        r = np.trace(r, axis1=q, axis2=q + cn)
        cn -= 1
    dk = 2 ** len(keep)
    return r.reshape((dk, dk))


def vn(rho):
    e = np.real(np.linalg.eigvalsh(rho))
    e = e[e > 1e-12]
    return float(-np.sum(e * np.log2(e))) if len(e) else 0.0


def mi(rho_full, A, B):
    return vn(ptrace_keep(rho_full, A)) + vn(ptrace_keep(rho_full, B)) - vn(ptrace_keep(rho_full, list(A) + list(B)))


rho_vec = rho0.flatten()
print("t       MM       MI(0,4)  MI(1,3)  pop_0   pop_1   pop_3   pop_4")
for t in (0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.6, 0.7, 0.8, 1.0, 1.4):
    rho_t = evolve_vec(rho_vec, t).reshape(d, d)
    rho_t = (rho_t + rho_t.conj().T) / 2
    mi04 = mi(rho_t, [0], [4])
    mi13 = mi(rho_t, [1], [3])
    pops = [ptrace_keep(rho_t, [k])[1, 1].real for k in range(N)]
    print(f"{t:.3f}  {mi04 + mi13:.4f}  {mi04:.4f}  {mi13:.4f}  {pops[0]:.4f}  {pops[1]:.4f}  {pops[3]:.4f}  {pops[4]:.4f}")

mm_max, t_max = 0.0, 0.0
for t in np.linspace(0, 2.0, 401):
    rho_t = evolve_vec(rho_vec, t).reshape(d, d)
    rho_t = (rho_t + rho_t.conj().T) / 2
    m = mi(rho_t, [0], [4]) + mi(rho_t, [1], [3])
    if m > mm_max:
        mm_max, t_max = m, t

print(f"\nMax MM over t in [0, 2]: {mm_max:.4f} at t = {t_max:.3f}")
print(f"MM(0) analytic (F75):     1.2451")
print(f"Ratio numerical / F75:    {mm_max / 1.2451:.4f}")
