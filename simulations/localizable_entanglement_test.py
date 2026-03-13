"""
LOCALIZABLE ENTANGLEMENT TEST
==============================
Benchmark CΨ against localizable entanglement (LE) and concurrence of
assistance (CoA) on the star topology. Tests whether CΨ is a witness
for coherent, localizable pairwise entanglement.

Date: 2026-03-08
See: experiments/LOCALIZABLE_ENTANGLEMENT_BENCHMARK.md
"""
import numpy as np
import sys
sys.path.insert(0, '.')
# NOTE: gpt_code module was renamed to star_topology_v3
import gpt_code as gpt


def partial_trace_S(rho_3q, outcome_proj):
    I4 = np.eye(4, dtype=complex)
    P = np.kron(outcome_proj, I4)
    rho_post = P @ rho_3q @ P
    rho_post_r = rho_post.reshape(2, 4, 2, 4)
    rho_AB = np.trace(rho_post_r, axis1=0, axis2=2)
    return rho_AB


def measurement_basis_projectors(theta, phi):
    c = np.cos(theta / 2)
    s = np.sin(theta / 2)
    ep = np.exp(1j * phi)
    v0 = np.array([c, ep * s], dtype=complex)
    v1 = np.array([s, -ep * c], dtype=complex)
    P0 = np.outer(v0, v0.conj())
    P1 = np.outer(v1, v1.conj())
    return P0, P1


def localizable_entanglement_AB(rho_3q, n_angles=20):
    best_avg_C = 0.0
    best_angles = (0, 0)
    thetas = np.linspace(0, np.pi, n_angles)
    phis = np.linspace(0, 2 * np.pi, n_angles)
    for theta in thetas:
        for phi in phis:
            P0, P1 = measurement_basis_projectors(theta, phi)
            rho_AB_0 = partial_trace_S(rho_3q, P0)
            rho_AB_1 = partial_trace_S(rho_3q, P1)
            p0 = np.real(np.trace(rho_AB_0))
            p1 = np.real(np.trace(rho_AB_1))
            avg_C = 0.0
            if p0 > 1e-10:
                avg_C += p0 * gpt.concurrence_two_qubit(rho_AB_0 / p0)
            if p1 > 1e-10:
                avg_C += p1 * gpt.concurrence_two_qubit(rho_AB_1 / p1)
            if avg_C > best_avg_C:
                best_avg_C = avg_C
                best_angles = (theta, phi)
    return best_avg_C, best_angles


def concurrence_of_assistance(rho_AB):
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    rho_tilde = np.kron(sy, sy) @ rho_AB.conj() @ np.kron(sy, sy)
    R = rho_AB @ rho_tilde
    eigenvalues = np.sort(np.real(np.sqrt(np.maximum(
        np.linalg.eigvals(R), 0))))[::-1]
    return float(np.sum(eigenvalues))
