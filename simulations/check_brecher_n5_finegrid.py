#!/usr/bin/env python3
"""Fine-grid trajectory check for |+-+-+> at N=5 under uniform J and uniform gamma.

Compares against C# brecher output to determine if Python's original coarse
grid (np.linspace(0.1, 15.0, 40), step ~0.38) missed the true Peak SumMI.
"""
import sys
import numpy as np
from scipy.linalg import expm

# Minimal setup copied from shadow_lens_broken.py
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)
plus = (up + dn) / np.sqrt(2)
minus = (up - dn) / np.sqrt(2)

N = 5
d = 32
d2 = 1024


def site_op(op, k):
    ops = [I2] * N
    ops[k] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def build_H(J_vec):
    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        block = (site_op(sx, b) @ site_op(sx, b + 1)
                 + site_op(sy, b) @ site_op(sy, b + 1)
                 + site_op(sz, b) @ site_op(sz, b + 1))
        H += J_vec[b] * block
    return H


def build_L(H, gammas):
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(sz, k)
        L += gammas[k] * (np.kron(Zk, Zk.conj()) - np.eye(d2))
    return L


def evolve(L, rho0, t):
    v = expm(L * t) @ rho0.flatten()
    rho = v.reshape(d, d)
    return (rho + rho.conj().T) / 2


def ptrace_keep(rho, keep):
    keep = list(keep)
    trace_out = [q for q in range(N) if q not in keep]
    dims = [2] * N
    reshaped = rho.reshape(dims + dims)
    current_n = N
    for q in sorted(trace_out, reverse=True):
        reshaped = np.trace(reshaped, axis1=q, axis2=q + current_n)
        current_n -= 1
    d_k = 2 ** len(keep)
    return reshaped.reshape((d_k, d_k))


def von_neumann(rho):
    evals = np.linalg.eigvalsh(rho)
    evals = np.real(evals)
    evals = evals[evals > 1e-12]
    if len(evals) == 0:
        return 0.0
    return float(-np.sum(evals * np.log2(evals)))


def mutual_information(rho_full, sites_A, sites_B):
    rho_A = ptrace_keep(rho_full, sites_A)
    rho_B = ptrace_keep(rho_full, sites_B)
    rho_AB = ptrace_keep(rho_full, list(sites_A) + list(sites_B))
    return von_neumann(rho_A) + von_neumann(rho_B) - von_neumann(rho_AB)


def sum_mi_adjacent(rho_full):
    total = 0.0
    for i in range(N - 1):
        total += mutual_information(rho_full, [i], [i + 1])
    return total


# |+-+-+> initial state
psi = plus
for ket in [minus, plus, minus, plus]:
    psi = np.kron(psi, ket)

rho0 = np.outer(psi, psi.conj())
rho0 = (rho0 + rho0.conj().T) / 2
rho0 /= np.trace(rho0).real

# Fine grid
J_vec = [1.0] * 4
gammas = [0.05] * 5

H = build_H(J_vec)
L = build_L(H, gammas)

print("t, SumMI")
peak_t = 0.0
peak_smi = 0.0
t_vals = np.concatenate([
    np.linspace(0.0, 1.0, 101),   # fine early: t=0.00 to 1.0 in 0.01 steps
    np.linspace(1.1, 3.0, 20),    # medium: t=1.1 to 3.0 in 0.1 steps
    np.linspace(3.5, 15.0, 24),   # coarse late: every 0.5
])

for t in t_vals:
    rho_t = evolve(L, rho0, t)
    smi = sum_mi_adjacent(rho_t)
    if smi > peak_smi:
        peak_smi = smi
        peak_t = t
    if t <= 0.5 or abs(t % 0.1) < 0.02:
        print(f"  t={t:.3f}  SumMI={smi:.6f}")

print()
print(f"FINE-GRID PEAK: SumMI={peak_smi:.4f} at t={peak_t:.3f}")
print(f"Python original coarse-grid (40 pts in [0.1, 15]): 1.3206 at t=0.10")
print(f"C# brecher with fine grid:                         2.5694 at t=0.20")
