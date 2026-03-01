"""
Bridge Closure: Pre-Shared Entanglement = Shared Randomness
============================================================
Proves that CΨ fingerprints cannot carry information beyond what
a classical pre-shared key provides, when A and B are separated
with zero coupling.

Result: Bridge hypothesis is dead.
See experiments/BRIDGE_CLOSURE.md for the full analysis.

Date: 2026-03-01
Authors: Thomas Wicht, Claude (Anthropic)
"""

import numpy as np
from qutip import (
    basis, tensor, ket2dm, sigmax, sigmay, sigmaz, qeye,
    mesolve, expect, entropy_vn
)

zero = basis(2, 0)
one = basis(2, 1)
plus = (zero + one).unit()
minus = (zero - one).unit()

gamma = 0.05
H = 0 * tensor(sigmax(), sigmax())  # J = 0
c_ops = [
    np.sqrt(gamma) * tensor(sigmaz(), qeye(2)),
    np.sqrt(gamma) * tensor(qeye(2), sigmaz()),
]
times = np.linspace(0, 10, 500)

# === Test 1: B's action is invisible to A ===
bell = (tensor(zero, zero) + tensor(one, one)).unit()
result = mesolve(H, ket2dm(bell), times, c_ops, [])
rho_t2 = result.states[100]  # t ≈ 2.0

# Branch: B does nothing vs B measures Z vs B measures X
P0_B = tensor(qeye(2), zero * zero.dag())
P1_B = tensor(qeye(2), one * one.dag())
P_plus_B = tensor(qeye(2), plus * plus.dag())
P_minus_B = tensor(qeye(2), minus * minus.dag())

rho_Bz = P0_B * rho_t2 * P0_B.dag() + P1_B * rho_t2 * P1_B.dag()
rho_Bx = P_plus_B * rho_t2 * P_plus_B.dag() + P_minus_B * rho_t2 * P_minus_B.dag()

r_nothing = mesolve(H, rho_t2, times, c_ops, [])
r_Bz = mesolve(H, rho_Bz, times, c_ops, [])
r_Bx = mesolve(H, rho_Bx, times, c_ops, [])

print("=== B's action invisible to A ===")
for i in [0, 50, 100, 200, 400]:
    d_z = np.linalg.norm((r_nothing.states[i].ptrace(0) - r_Bz.states[i].ptrace(0)).full())
    d_x = np.linalg.norm((r_nothing.states[i].ptrace(0) - r_Bx.states[i].ptrace(0)).full())
    print(f"  t={times[i]:.1f}: ||Δ(nothing,Bz)||={d_z:.2e}, ||Δ(nothing,Bx)||={d_x:.2e}")
    assert d_z < 1e-12 and d_x < 1e-12

# === Test 2: Product states indistinguishable to A ===
states = {
    "|++>": tensor(plus, plus),
    "|+0>": tensor(plus, zero),
    "|+->": tensor(plus, minus),
}

print("\n=== A cannot distinguish |++>, |+0>, |+-> ===")
rho_As = {}
for name, psi in states.items():
    rho_A = ket2dm(psi).ptrace(0)
    rho_As[name] = rho_A
    print(f"  {name}: rho_A = {rho_A.full().diagonal().real}")

for n1, n2 in [("|++>", "|+0>"), ("|++>", "|+->"), ("|+0>", "|+->")]:
    d = np.linalg.norm((rho_As[n1] - rho_As[n2]).full())
    print(f"  ||rho_A({n1}) - rho_A({n2})|| = {d:.2e}")
    assert d < 1e-12

# === Test 3: Bell+ gives A zero information forever ===
print("\n=== Bell+ qubit = noise forever ===")
for t_idx in [0, 100, 200, 400, 499]:
    rho_A = result.states[t_idx].ptrace(0)
    S = entropy_vn(rho_A, 2)
    info = 1.0 - S
    print(f"  t={times[t_idx]:.1f}: S(rho_A)={S:.6f}, info={info:.6f} bits")
    assert abs(info) < 1e-10

print("\nAll tests passed. Bridge is dead.")
print("Pre-shared entanglement without a channel = shared randomness.")
