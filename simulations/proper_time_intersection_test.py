"""
PROPER TIME INTERSECTION TEST
==============================
Tests whether observers with different noise rates (gamma) see different
CΨ windows on the same AB pair. The connection (CoA) is always there,
but its visible expression (CΨ) depends on who is looking.

Date: 2026-03-08
See: experiments/OBSERVER_DEPENDENT_VISIBILITY.md
"""
import numpy as np
import sys
sys.path.insert(0, '.')
import gpt_code as gpt

configs = [
    {"gamma_A": 0.03, "gamma_B": 0.05, "gamma_S": 0.05, "label": "slow A"},
    {"gamma_A": 0.05, "gamma_B": 0.05, "gamma_S": 0.05, "label": "equal"},
    {"gamma_A": 0.10, "gamma_B": 0.05, "gamma_S": 0.05, "label": "fast A"},
    {"gamma_A": 0.20, "gamma_B": 0.05, "gamma_S": 0.05, "label": "very fast A"},
]

results = {}
for cfg in configs:
    H = gpt.star_hamiltonian_n(n_observers=2, J_SA=1.0, J_SB=2.0)
    L_ops = gpt.dephasing_ops_n([cfg["gamma_S"], cfg["gamma_A"], cfg["gamma_B"]])
    psi = gpt.bell_phi_plus()
    psi = np.kron(psi, gpt.plus_state())
    rho = gpt.density_from_statevector(psi)
    dt = 0.005
    data = {"t": [], "tau_A": [], "tau_B": [], "cpsi_AB": []}
    for step in range(2001):
        t = step * dt
        if step % 10 == 0:
            rho_AB = gpt.partial_trace_keep(rho, keep=[1, 2], n_qubits=3)
            c_ab = gpt.concurrence_two_qubit(rho_AB)
            psi_ab = gpt.psi_norm(rho_AB)
            data["t"].append(t)
            data["tau_A"].append(cfg["gamma_A"] * t)
            data["tau_B"].append(cfg["gamma_B"] * t)
            data["cpsi_AB"].append(c_ab * psi_ab)
        if step < 2000:
            rho = gpt.rk4_step(rho, H, L_ops, dt)
    results[cfg["label"]] = data
