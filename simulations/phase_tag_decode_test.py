"""
PHASE-TAG-AND-DECODE TEST
==========================
Tests whether local phase operations on mediator S produce readable
changes in the AB reduced state. Result: yes, but this is standard
Heisenberg exchange dynamics, not a new transport mechanism.

Date: 2026-03-08
See: experiments/WHATS_INSIDE_THE_WINDOWS.md
"""
import numpy as np
import sys
sys.path.insert(0, '.')
import gpt_code as gpt


def apply_on_S(rho, operator):
    """Apply a single-qubit operator on qubit S (index 0) in 3-qubit system."""
    U = np.kron(np.kron(operator, np.eye(2)), np.eye(2))
    return U @ rho @ U.conj().T


def trace_distance(rho1, rho2):
    """Trace distance between two density matrices."""
    diff = rho1 - rho2
    evals = np.linalg.eigvalsh(diff)
    return 0.5 * np.sum(np.abs(evals))


# Interventions on S
I2 = np.eye(2, dtype=complex)
phi = np.pi / 4

interventions = {
    "I (nothing)":   I2,
    "Rz(+pi/4)":    np.array([[np.exp(-1j*phi/2), 0],
                               [0, np.exp(1j*phi/2)]], dtype=complex),
    "Rz(-pi/4)":    np.array([[np.exp(1j*phi/2), 0],
                               [0, np.exp(-1j*phi/2)]], dtype=complex),
    "Rz(+pi/2)":    np.array([[np.exp(-1j*np.pi/4), 0],
                               [0, np.exp(1j*np.pi/4)]], dtype=complex),
    "X (bit-flip)":  np.array([[0, 1], [1, 0]], dtype=complex),
    "Rx(pi/2)":      np.array([[np.cos(np.pi/4), -1j*np.sin(np.pi/4)],
                                [-1j*np.sin(np.pi/4), np.cos(np.pi/4)]],
                               dtype=complex),
}


def run_tagged_evolution(rho_init, H, L_ops, write_t, readout_t,
                          intervention_op, dt=0.005):
    """Evolve to write_t, apply intervention on S, evolve to readout_t."""
    rho = rho_init.copy()
    steps_to_write = int(write_t / dt)
    steps_to_read = int(readout_t / dt)
    for step in range(steps_to_read + 1):
        if step == steps_to_write:
            rho = apply_on_S(rho, intervention_op)
        if step < steps_to_read:
            rho = gpt.rk4_step(rho, H, L_ops, dt)
    return gpt.partial_trace_keep(rho, keep=[1, 2], n_qubits=3)
