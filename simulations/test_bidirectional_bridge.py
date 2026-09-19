#!/usr/bin/env python3
"""Finite transformed-state comparison for one N=3 calculation.

The script evolves |+>^3 under the named generator, applies the named Pi map
to each sampled operator, and evaluates the chosen scalar readout on the
original density matrix and on that transformed operator.  Pi(rho) here is a
linear operator image, not a second density matrix: this computation does not
establish trace, Hermiticity, or positivity of the image.  It is not two
ontological regimes and is not evidence for the Spectral Midpoint construction.
"""

import numpy as np
from scipy.linalg import expm
from itertools import product as iprod

# === Pauli matrices ===
I2 = np.eye(2, dtype=complex)
X = np.array([[0,1],[1,0]], dtype=complex)
Y = np.array([[0,-1j],[1j,0]], dtype=complex)
Z = np.array([[1,0],[0,-1]], dtype=complex)
PAULIS = [I2, X, Y, Z]

def kron_list(ops):
    """Tensor product of list of operators."""
    result = ops[0]
    for op in ops[1:]:
        result = np.kron(result, op)
    return result

def build_heisenberg_ham(N, J=1.0):
    """Heisenberg XXX Hamiltonian on N-qubit chain."""
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for i in range(N-1):
        for P in [X, Y, Z]:
            ops = [I2]*N
            ops[i] = P
            ops[i+1] = P
            H += J * kron_list(ops)
    return H

def lindblad_rhs(rho, H, gammas, N):
    """drho/dt = -i[H,rho] + sum_k gamma_k (Z_k rho Z_k - rho)."""
    d = 2**N
    drho = -1j * (H @ rho - rho @ H)
    for k in range(N):
        ops = [I2]*N
        ops[k] = Z
        Zk = kron_list(ops)
        drho += gammas[k] * (Zk @ rho @ Zk - rho)
    return drho

def rk4_step(rho, H, gammas, N, dt):
    k1 = lindblad_rhs(rho, H, gammas, N)
    k2 = lindblad_rhs(rho + 0.5*dt*k1, H, gammas, N)
    k3 = lindblad_rhs(rho + 0.5*dt*k2, H, gammas, N)
    k4 = lindblad_rhs(rho + dt*k3, H, gammas, N)
    return rho + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)

def partial_trace(rho, N, keep):
    """Partial trace: keep specified qubits, trace out the rest."""
    d = 2**N
    n_keep = len(keep)
    d_keep = 2**n_keep
    d_trace = d // d_keep
    
    # Build the partial trace by explicit summation
    rho_reduced = np.zeros((d_keep, d_keep), dtype=complex)
    
    trace_qubits = [q for q in range(N) if q not in keep]
    
    for i in range(d_keep):
        for j in range(d_keep):
            # Convert i,j to binary for kept qubits
            i_bits = [(i >> (n_keep-1-k)) & 1 for k in range(n_keep)]
            j_bits = [(j >> (n_keep-1-k)) & 1 for k in range(n_keep)]
            
            # Sum over traced qubits
            for t_val in iprod([0,1], repeat=len(trace_qubits)):
                # Build full index
                row_bits = [0]*N
                col_bits = [0]*N
                for k, q in enumerate(keep):
                    row_bits[q] = i_bits[k]
                    col_bits[q] = j_bits[k]
                for k, q in enumerate(trace_qubits):
                    row_bits[q] = t_val[k]
                    col_bits[q] = t_val[k]
                
                row_idx = sum(b << (N-1-q) for q, b in enumerate(row_bits))
                col_idx = sum(b << (N-1-q) for q, b in enumerate(col_bits))
                rho_reduced[i, j] += rho[row_idx, col_idx]
    
    return rho_reduced

def purity(rho):
    return np.real(np.trace(rho @ rho))

def l1_coherence(rho):
    """L1 norm of off-diagonal elements."""
    d = rho.shape[0]
    return np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho)))

def selected_scalar_readout(operator_2qubit):
    """Chosen algebraic Tr(A^2)*l1(A)/3 readout; no state claim for Pi(rho)."""
    p = purity(operator_2qubit)
    l1 = l1_coherence(operator_2qubit)
    return p * l1 / 3.0  # d=4, so d-1=3

def mutual_information(rho, N, qA, qB):
    """MI between qubit qA and qubit qB."""
    rhoA = partial_trace(rho, N, [qA])
    rhoB = partial_trace(rho, N, [qB])
    rhoAB = partial_trace(rho, N, [qA, qB])
    
    def von_neumann(r):
        evals = np.real(np.linalg.eigvalsh(r))
        evals = evals[evals > 1e-15]
        return -np.sum(evals * np.log2(evals))
    
    return von_neumann(rhoA) + von_neumann(rhoB) - von_neumann(rhoAB)


# === The Pi operator: constructs the algebraic image ===
def pi_map_index(a):
    """Map single-site Pauli index under Pi: I->X, X->I, Y->iZ, Z->iY."""
    if a == 0: return 1, 1.0        # I -> X, factor +1
    if a == 1: return 0, 1.0        # X -> I, factor +1
    if a == 2: return 3, 1j         # Y -> iZ, factor +i
    if a == 3: return 2, 1j         # Z -> iY, factor +i

def pauli_decompose(rho, N):
    """Decompose density matrix into Pauli coefficients."""
    d = 2**N
    coeffs = {}
    for indices in iprod(range(4), repeat=N):
        ops = [PAULIS[i] for i in indices]
        P = kron_list(ops)
        c = np.trace(rho @ P) / d  # c = Tr(rho * P) / 2^N
        if abs(c) > 1e-15:
            coeffs[indices] = c
    return coeffs

def pauli_recompose(coeffs, N):
    """Recompose density matrix from Pauli coefficients."""
    d = 2**N
    rho = np.zeros((d, d), dtype=complex)
    for indices, c in coeffs.items():
        ops = [PAULIS[i] for i in indices]
        P = kron_list(ops)
        rho += c * P
    return rho

def apply_pi(rho, N):
    """Apply the canonical linear Pi map to an operator."""
    coeffs = pauli_decompose(rho, N)
    pi_coeffs = {}
    for indices, c in coeffs.items():
        new_indices = []
        total_factor = 1.0
        for a in indices:
            new_a, factor = pi_map_index(a)
            new_indices.append(new_a)
            total_factor *= factor
        pi_coeffs[tuple(new_indices)] = c * total_factor
    return pauli_recompose(pi_coeffs, N)


# === Main test ===
def run_test():
    N = 3
    d = 2**N
    gamma_base = 0.05
    eps = 0.001
    
    # Sacrifice-zone formula: all noise on qubit 0
    gammas = [N * gamma_base - (N-1) * eps] + [eps] * (N-1)
    assert abs(sum(gammas) - 0.15) < 1e-15

    print(f"=== FINITE PI-TRANSFORMED-STATE COMPARISON (N={N}) ===")
    print("Named Pi map; chosen scalar readout on the original density matrix "
          "and its linear operator image, not a density matrix.")
    print("This is not two regimes and does not establish a physical bridge.")
    print(f"Gammas: {[f'{g:.3f}' for g in gammas]}")
    print(f"Sum_gamma = {sum(gammas):.4f}; finite centre-rate coordinate only")
    print()
    
    H = build_heisenberg_ham(N)
    
    # Initial state: |+>^N
    psi = np.ones(d, dtype=complex) / np.sqrt(d)
    rho = np.outer(psi, psi.conj())
    
    dt = 0.02
    t_max = 15.0
    t_meas = 0.25
    
    rows = []
    t = 0
    next_meas = 0
    
    while t <= t_max + 1e-9:
        if t >= next_meas - 1e-9:
            # Mutual information belongs to the physical density matrix.
            mi01 = mutual_information(rho, N, 0, 1)
            mi02 = mutual_information(rho, N, 0, 2)
            mi12 = mutual_information(rho, N, 1, 2)
            adjacent_sum_mi = mi01 + mi12
            
            rho01_original = partial_trace(rho, N, [0, 1])
            rho02_original = partial_trace(rho, N, [0, 2])
            rho12_original = partial_trace(rho, N, [1, 2])
            scalar01_original = selected_scalar_readout(rho01_original)
            scalar02_original = selected_scalar_readout(rho02_original)
            scalar12_original = selected_scalar_readout(rho12_original)

            # Pi(rho) is a linear operator image, not a physical state.
            operator_pi = apply_pi(rho, N)
            operator01_pi = partial_trace(operator_pi, N, [0, 1])
            operator02_pi = partial_trace(operator_pi, N, [0, 2])
            operator12_pi = partial_trace(operator_pi, N, [1, 2])
            scalar01_pi = selected_scalar_readout(operator01_pi)
            scalar02_pi = selected_scalar_readout(operator02_pi)
            scalar12_pi = selected_scalar_readout(operator12_pi)

            row = {
                "t": t,
                "mi": (mi01, mi02, mi12),
                "adjacent_sum_mi": adjacent_sum_mi,
                "original": (scalar01_original, scalar02_original,
                             scalar12_original),
                "pi_operator": (scalar01_pi, scalar02_pi, scalar12_pi),
            }
            rows.append(row)

            next_meas += t_meas
        
        rho = rk4_step(rho, H, gammas, N, dt)
        t += dt
    
    best_index = max(range(len(rows)), key=lambda index: rows[index]["adjacent_sum_mi"])
    print(f"{'T':>5}  {'MI01':>7} {'MI02':>7} {'MI12':>7}  |  "
          f"{'R01_orig':>9} {'R02_orig':>9} {'R12_orig':>9}  |  "
          f"{'R01_PiOp':>9} {'R02_PiOp':>9} {'R12_PiOp':>9}  |  "
          f"{'orig<.25':>8} {'PiOp<.25':>8}")
    print("-" * 132)
    for index, row in enumerate(rows):
        mi01, mi02, mi12 = row["mi"]
        original = row["original"]
        pi_operator = row["pi_operator"]
        original_below = sum(value < 0.25 for value in original)
        pi_below = sum(value < 0.25 for value in pi_operator)
        marker = " <-- GLOBAL ARGMAX OF ADJACENT-PAIR SUM" if index == best_index else ""
        print(f"{row['t']:5.2f}  {mi01:7.4f} {mi02:7.4f} {mi12:7.4f}  |  "
              f"{original[0]:9.4f} {original[1]:9.4f} {original[2]:9.4f}  |  "
              f"{pi_operator[0]:9.4f} {pi_operator[1]:9.4f} {pi_operator[2]:9.4f}  |  "
              f"{original_below:>8} {pi_below:>8}{marker}")

    best = rows[best_index]
    print()
    print(f"Global argmax of adjacent-pair SumMI = {best['adjacent_sum_mi']:.6f} "
          f"at t = {best['t']:.2f}")
    print("Finite transformed-state comparison only: the Pi image is not "
          "assumed trace-one, Hermitian, or positive. The scalar comparison "
          "is not evidence for the Spectral Midpoint construction.")

if __name__ == "__main__":
    run_test()
