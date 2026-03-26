#!/usr/bin/env python3
"""
Test: Bidirectional Bridge Between Two Sides
=============================================
Theory (Tom): Side A and Side B each have their own CΨ = 1/4.
The qubit (d=2) is the mediator. MI peaks when BOTH sides cross.

We compute CΨ from OUR perspective and from the Π perspective
simultaneously, and check if both cross 1/4 at the PeakMI time.

N=3, sacrifice-zone formula, |+>^3 initial state.
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

def cpsi(rho_2qubit):
    """CΨ = Tr(rho^2) * L1/(d-1) for a 2-qubit subsystem."""
    p = purity(rho_2qubit)
    l1 = l1_coherence(rho_2qubit)
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


# === The Pi operator: maps to the "other side" ===
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
    """Apply Pi operator to density matrix: transform to the other side."""
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
    # gammas = [0.148, 0.001, 0.001]
    
    print(f"=== BIDIRECTIONAL BRIDGE TEST (N={N}) ===")
    print(f"Gammas: {[f'{g:.3f}' for g in gammas]}")
    print(f"Sum_gamma = {sum(gammas):.4f}, Midpoint = {sum(gammas):.4f}")
    print()
    
    H = build_heisenberg_ham(N)
    
    # Initial state: |+>^N
    psi = np.ones(d, dtype=complex) / np.sqrt(d)
    rho = np.outer(psi, psi.conj())
    
    dt = 0.02
    t_max = 15.0
    t_meas = 0.25
    
    print(f"{'T':>5}  {'MI01':>7} {'MI02':>7} {'MI12':>7}  |  "
          f"{'CP01_A':>7} {'CP02_A':>7} {'CP12_A':>7}  |  "
          f"{'CP01_B':>7} {'CP02_B':>7} {'CP12_B':>7}  |  "
          f"{'A<.25':>5} {'B<.25':>5}")
    print("-" * 120)
    
    best_mi = 0
    best_t = 0
    t = 0
    next_meas = 0
    
    while t <= t_max + 1e-9:
        if t >= next_meas - 1e-9:
            # Compute MI from our side
            mi01 = mutual_information(rho, N, 0, 1)
            mi02 = mutual_information(rho, N, 0, 2)
            mi12 = mutual_information(rho, N, 1, 2)
            sum_mi = mi01 + mi12  # adjacent pairs only
            
            # CΨ from OUR side (standard)
            rho01_A = partial_trace(rho, N, [0, 1])
            rho02_A = partial_trace(rho, N, [0, 2])
            rho12_A = partial_trace(rho, N, [1, 2])
            cpsi01_A = cpsi(rho01_A)
            cpsi02_A = cpsi(rho02_A)
            cpsi12_A = cpsi(rho12_A)
            
            # Apply Pi: transform to the OTHER side
            rho_pi = apply_pi(rho, N)
            
            # CΨ from the PI side
            rho01_B = partial_trace(rho_pi, N, [0, 1])
            rho02_B = partial_trace(rho_pi, N, [0, 2])
            rho12_B = partial_trace(rho_pi, N, [1, 2])
            cpsi01_B = cpsi(rho01_B)
            cpsi02_B = cpsi(rho02_B)
            cpsi12_B = cpsi(rho12_B)
            
            # Count crossings
            a_below = sum(1 for c in [cpsi01_A, cpsi02_A, cpsi12_A] if c < 0.25)
            b_below = sum(1 for c in [cpsi01_B, cpsi02_B, cpsi12_B] if c < 0.25)
            
            marker = ""
            if sum_mi > best_mi:
                best_mi = sum_mi
                best_t = t
            if abs(t - best_t) < 0.01 and t > 0.5:
                marker = " <-- PEAK MI"
            
            print(f"{t:5.2f}  {mi01:7.4f} {mi02:7.4f} {mi12:7.4f}  |  "
                  f"{cpsi01_A:7.4f} {cpsi02_A:7.4f} {cpsi12_A:7.4f}  |  "
                  f"{cpsi01_B:7.4f} {cpsi02_B:7.4f} {cpsi12_B:7.4f}  |  "
                  f"{a_below:>5} {b_below:>5}{marker}")
            
            next_meas += t_meas
        
        rho = rk4_step(rho, H, gammas, N, dt)
        t += dt
    
    print()
    print(f"Peak SumMI = {best_mi:.6f} at t = {best_t:.2f}")
    print()
    
    # Summary: geometric mean at peak
    print("=== AT PEAK MI: GEOMETRIC MEAN OF BOTH PERSPECTIVES ===")

if __name__ == "__main__":
    run_test()

