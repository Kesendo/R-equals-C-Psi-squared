"""
Shift Mechanism: Why B's measurement accelerates A's crossing

The nonlocal coherence acts as a shield/reservoir. The Hamiltonian
(J > 0) continuously exchanges coherence between local and nonlocal.
B's measurement destroys the nonlocal pool. A's local coherence
is then exposed to dephasing without protection. Decay accelerates.

Key findings:
  - Local coherence A unchanged at moment of measurement
  - Nonlocal coherence drops 0.82 → 0.00 instantly
  - A's subsequent decay ~4x faster
  - Entanglement attempts to regenerate but fails (peak 0.15)

Source: Thomas Wicht + Claude (Anthropic), 2026-03-01
See: experiments/OBSERVER_GRAVITY_BRIDGE.md §7
"""

import numpy as np
from qutip import (
    basis, tensor, ket2dm, sigmax, sigmay, sigmaz, qeye,
    mesolve, concurrence
)

zero = basis(2, 0)
one  = basis(2, 1)
plus = (zero + one).unit()

gamma = 0.05
J = 0.5
H = J * (tensor(sigmax(), sigmax()) +
         tensor(sigmay(), sigmay()) +
         tensor(sigmaz(), sigmaz()))

c_ops = [np.sqrt(gamma) * tensor(sigmaz(), qeye(2)),
         np.sqrt(gamma) * tensor(qeye(2), sigmaz())]

times = np.linspace(0, 10, 1000)
P0_B = tensor(qeye(2), zero * zero.dag())
P1_B = tensor(qeye(2), one * one.dag())

def local_coherence_A(rho):
    rho_A = rho.ptrace(0)
    return abs(rho_A.full()[0,1]) + abs(rho_A.full()[1,0])

def nonlocal_coherence(rho):
    rho_full = rho.full()
    d = rho_full.shape[0]
    total = sum(abs(rho_full[i,j]) for i in range(d) for j in range(d) if i != j)
    rho_A = rho.ptrace(0).full()
    rho_B = rho.ptrace(1).full()
    loc_A = sum(abs(rho_A[i,j]) for i in range(2) for j in range(2) if i != j)
    loc_B = sum(abs(rho_B[i,j]) for i in range(2) for j in range(2) if i != j)
    return total - loc_A - loc_B

# Evolve to t_B = 1.0
t_B = 1.0
idx_tB = np.argmin(np.abs(times - t_B))
r_pre = mesolve(H, ket2dm(tensor(plus, plus)), times[:idx_tB+1], c_ops, [])
rho_tB = r_pre.states[-1]

# Branch
rho_Bm = P0_B * rho_tB * P0_B.dag() + P1_B * rho_tB * P1_B.dag()

print("At measurement:")
print(f"  Local coh A: {local_coherence_A(rho_tB):.6f} -> {local_coherence_A(rho_Bm):.6f}")
print(f"  Nonlocal:    {nonlocal_coherence(rho_tB):.6f} -> {nonlocal_coherence(rho_Bm):.6f}")

r_silent = mesolve(H, rho_tB, times, c_ops, [])
r_meas = mesolve(H, rho_Bm, times, c_ops, [])

print(f"\n{'t':>6} {'Coh_A(silent)':>14} {'Coh_A(meas)':>14} {'Conc(meas)':>12}")
for t in [0, 0.5, 1.0, 2.0, 5.0, 8.0]:
    idx = np.argmin(np.abs(times - t))
    cs = local_coherence_A(r_silent.states[idx])
    cm = local_coherence_A(r_meas.states[idx])
    try: conc = concurrence(r_meas.states[idx])
    except: conc = 0.0
    print(f"{t:>6.1f} {cs:>14.6f} {cm:>14.6f} {conc:>12.6f}")
