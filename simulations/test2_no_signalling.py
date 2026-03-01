"""
Test #2: No-Signalling Test for Bridge Protocol
================================================
Test #1 (2026-02-24) showed all 15 local observables distinguish Bell+ from
product states WITH Heisenberg coupling (J=0.5).

Test #2: Same observables, J=0. No physical coupling between qubits.
Each qubit experiences independent local dephasing only.

Question: Can qubit A's local observables distinguish whether the JOINT state
is Bell+ vs |++> vs |+0> vs |00>, when there is NO interaction between A and B?

Result: No-signalling holds exactly. rho_A unchanged. But CΨ drops from
0.500 to 0.250 because global purity C drops while Ψ stays constant.
See experiments/NO_SIGNALLING_BOUNDARY.md for the full analysis.

Date: 2026-03-01
Authors: Thomas Wicht, Claude (Anthropic)
"""

import numpy as np
from qutip import (
    basis, tensor, ket2dm, sigmax, sigmay, sigmaz, qeye,
    mesolve, expect, entropy_vn
)

# Parameters
gamma = 0.05
J_bridge = 0.0  # NO coupling — this is the test
t_max = 5.0
n_steps = 200
times = np.linspace(0, t_max, n_steps)

# Basis
zero = basis(2, 0)
one = basis(2, 1)
plus = (zero + one).unit()
minus = (zero - one).unit()

states = {
    "Bell+":  (tensor(zero, zero) + tensor(one, one)).unit(),
    "|++>":   tensor(plus, plus),
    "|+0>":   tensor(plus, zero),
    "|00>":   tensor(zero, zero),
    "|+->":   tensor(plus, minus),
}

H = 0 * tensor(sigmax(), sigmax())  # zero Hamiltonian

c_ops = [
    np.sqrt(gamma) * tensor(sigmaz(), qeye(2)),
    np.sqrt(gamma) * tensor(qeye(2), sigmaz()),
]

# === Part 1: Time evolution with local observables ===

results = {}
for name, psi in states.items():
    rho0 = ket2dm(psi)
    result = mesolve(H, rho0, times, c_ops, [])
    sx_vals, sy_vals, sz_vals = [], [], []
    purity_A, entropy_A, cpsi_vals = [], [], []

    for i, t in enumerate(times):
        rho_t = result.states[i]
        rho_A = rho_t.ptrace(0)
        sx_vals.append(expect(sigmax(), rho_A))
        sy_vals.append(expect(sigmay(), rho_A))
        sz_vals.append(expect(sigmaz(), rho_A))
        purity_A.append((rho_A * rho_A).tr().real)
        entropy_A.append(entropy_vn(rho_A, 2))
        C = (rho_t * rho_t).tr().real
        evals = rho_A.eigenenergies()
        Psi = max(evals).real
        cpsi_vals.append(C * Psi)

    results[name] = {
        'sx': np.array(sx_vals), 'sy': np.array(sy_vals),
        'sz': np.array(sz_vals), 'purity': np.array(purity_A),
        'entropy': np.array(entropy_A), 'cpsi': np.array(cpsi_vals),
    }

# === Part 2: The real no-signalling test ===

rho_bell = ket2dm(states["Bell+"])
rho_A_before = rho_bell.ptrace(0)

P0_B = tensor(qeye(2), zero * zero.dag())
P1_B = tensor(qeye(2), one * one.dag())
rho_after_meas = P0_B * rho_bell * P0_B.dag() + P1_B * rho_bell * P1_B.dag()
rho_A_after = rho_after_meas.ptrace(0)

dist = np.linalg.norm((rho_A_before - rho_A_after).full())
print(f"||Δρ_A|| = {dist:.10f}")  # → 0.0

# === Part 3: CΨ analysis ===

C_before = (rho_bell * rho_bell).tr().real        # 1.0
C_after = (rho_after_meas * rho_after_meas).tr().real  # 0.5
Psi = max(rho_A_before.eigenenergies()).real       # 0.5

print(f"CΨ before: {C_before * Psi:.3f}")  # 0.500
print(f"CΨ after:  {C_after * Psi:.3f}")   # 0.250
print(f"ΔCΨ:       {(C_after - C_before) * Psi:.3f}")  # -0.250
