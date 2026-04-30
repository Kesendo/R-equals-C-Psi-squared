"""Trotter + Z-dephasing prediction of multi-Pauli expectation tables per Hamiltonian category."""
from __future__ import annotations

import numpy as np
from scipy.linalg import expm

from ..pauli import ur_pauli, site_op


def _xneel_state(N):
    """X-Néel |+, -, +, -, ...⟩ on N qubits as 2^N density matrix."""
    ket_p = np.array([1, 1], dtype=complex) / np.sqrt(2)
    ket_m = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = ket_p if N >= 1 else None
    for q in range(1, N):
        psi = np.kron(psi, ket_m if q % 2 == 1 else ket_p)
    return np.outer(psi, psi.conj())


def _trotter_step_unitary(N, bonds, terms, J, delta_t):
    """First-order Trotter step unitary for ΣJ·terms over bonds."""
    U_step = np.eye(2 ** N, dtype=complex)
    for (P, Q) in terms:
        for (l, m) in bonds:
            ops = [ur_pauli('I')] * N
            ops[l] = ur_pauli(P)
            ops[m] = ur_pauli(Q)
            op_full = ops[0]
            for op in ops[1:]:
                op_full = np.kron(op_full, op)
            U_step = expm(-1j * J * delta_t * op_full) @ U_step
    return U_step


def _vec_F(M):
    return M.flatten('F')


def _unvec_F(v, d):
    return v.reshape((d, d), order='F')


def predict_signature_table(
    chain,
    terms_per_category,
    t=0.8,
    n_trotter=3,
    gamma_z=None,
    rho_0=None,
    q0=0,
    q2=None,
    paulis=('I', 'X', 'Y', 'Z'),
):
    """Predict 2-qubit Pauli expectation table for multiple Hamiltonian categories.

    For each category's Hamiltonian (sum of bilinears at chain bonds), evolve
    the initial state under a first-order Trotter circuit with n_trotter steps,
    inserting per-step Z-dephasing dissipation at rate γ_Z. Return expectation
    values of all 16 (or specified subset) 2-qubit Paulis on (q0, q2).

    This matches the IBM hardware Trotter circuit (PauliEvolutionGate × n_trotter
    + thermal channels) and is the model that fit the April 26 Marrakesh
    soft_break run to 0.0014 on ⟨X₀Z₂⟩.

    Args:
        chain: ChainSystem providing N, J, bonds.
        terms_per_category: dict {name: [(P, Q), ...]} — Hamiltonian terms per
            category (each term is a 2-Pauli-letter bilinear at every chain bond).
        t: total evolution time.
        n_trotter: number of Trotter steps.
        gamma_z: per-step Z-dephasing rate (defaults to chain.gamma_0).
        rho_0: initial density matrix (defaults to X-Néel |+,-,+,...⟩).
        q0, q2: qubit positions for 2-qubit observables (default endpoints).
        paulis: tuple of Pauli labels to enumerate. Default ('I','X','Y','Z').

    Returns:
        dict {category_name: {(p0, p2): expectation, ...}}
    """
    N = chain.N
    bonds = chain.bonds
    J = chain.J
    if gamma_z is None:
        gamma_z = chain.gamma_0
    if rho_0 is None:
        rho_0 = _xneel_state(N)
    if q2 is None:
        q2 = N - 1

    delta_t = t / n_trotter
    d = 2 ** N
    Id_d = np.eye(d, dtype=complex)

    # Z-dephasing channel matrix in vec basis (per Trotter step)
    L_deph = np.zeros((4 ** N, 4 ** N), dtype=complex)
    for site in range(N):
        Zl = site_op(N, site, 'Z')
        L_deph += gamma_z * (np.kron(Zl, Zl.conj()) - np.kron(Id_d, Id_d))
    M_deph = expm(L_deph * delta_t)

    # 2-qubit Pauli observables on (q0, q2)
    def two_qubit_pauli(p0, p2):
        ops = [ur_pauli('I')] * N
        ops[q0] = ur_pauli(p0)
        ops[q2] = ur_pauli(p2)
        out = ops[0]
        for op in ops[1:]:
            out = np.kron(out, op)
        return out

    obs_table = {(p0, p2): two_qubit_pauli(p0, p2) for p0 in paulis for p2 in paulis}

    predictions = {}
    for category, terms in terms_per_category.items():
        U_step = _trotter_step_unitary(N, bonds, terms, J, delta_t)
        rho = rho_0.copy()
        for _ in range(n_trotter):
            rho = U_step @ rho @ U_step.conj().T
            rho = _unvec_F(M_deph @ _vec_F(rho), d)
        cat_predictions = {}
        for key, obs in obs_table.items():
            cat_predictions[key] = float(np.real(np.trace(rho @ obs)))
        predictions[category] = cat_predictions

    return predictions
