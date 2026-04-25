#!/usr/bin/env python3
"""V-Effect weight-sector immunity: which XY-weight blocks of L break the
palindromic relation Π·L·Π⁻¹ = −L − 2Σγ·I, and why?

View from Level 1: Heisenberg saw exchange-coupled atoms, palindrome is the
spectral signature he could not directly read but every mode of his system
inherits. We verify here, with the palindrome as the anchor, where it
holds and where it breaks under V-Effect (two bonds at N=3).

Empirical claim from V_EFFECT_PALINDROME (March 2026):
- w=0 (all I, Z) and w=3 (all X, Y) are immune
- w=1 and w=2 (boundary, mixed) break for certain Hamiltonian combinations

Logical question: derive WHY in terms of the Π-action on each weight sector.

Π per site: I↔X (sign 1), Y↔Z (sign i). Maps weight w → N−w on each site.
For total: maps w → N−w. So Π pairs (w=0, w=3) and (w=1, w=2) at N=3.

Test cases:
1. H = XX + YY + ZZ on each bond (Heisenberg). Both Z₂ parities preserved.
   Expected: palindrome holds in all sectors.
2. H = XX + YY on each bond (XY model, Δ=0 XXZ). Both Z₂ parities preserved.
   Expected: palindrome holds in all sectors.
3. H = XX on bond (0,1) + XY on bond (1,2). Breaks bit_b (XY has 1 Y).
   Expected: w=0, w=3 still palindromic. w=1, w=2 break.
4. H = XX + XZ. Breaks bit_a (XZ has 1 X+Y across both — actually XZ has bit_a=1).
   Expected: also breaks the boundary sectors.
"""
import math
import sys

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


# Pauli matrices, indexed 0=I, 1=X, 2=Y, 3=Z
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, X, Y, Z]
NAMES = ['I', 'X', 'Y', 'Z']

# Π per site: I↔X (sign 1), Y↔Z (sign i)
PI_PERM = {0: 1, 1: 0, 2: 3, 3: 2}
PI_SIGN = {0: 1, 1: 1, 2: 1j, 3: 1j}

# bit_a = 1 if Pauli is X or Y (decaying under Z-dephasing); 0 for I, Z
BIT_A = {0: 0, 1: 1, 2: 1, 3: 0}
# bit_b = 1 if Pauli is Y or Z (Π² eigenvalue −1); 0 for I, X
BIT_B = {0: 0, 1: 0, 2: 1, 3: 1}


def pauli_string(indices):
    out = PAULIS[indices[0]]
    for i in indices[1:]:
        out = np.kron(out, PAULIS[i])
    return out


def k_to_indices(k, N):
    indices = []
    kk = k
    for _ in range(N):
        indices.append(kk % 4)
        kk //= 4
    return tuple(reversed(indices))


def indices_to_k(indices):
    k = 0
    for i in indices:
        k = k * 4 + i
    return k


def total_weight_a(indices):
    return sum(BIT_A[i] for i in indices)


def build_H_bonds(N, bonds_with_terms):
    """bonds_with_terms = list of (bond_indices, [pauli_pair_label]).
    bond_indices = (i, j); pauli_pair_label like 'XX', 'XY', etc.
    """
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    name_to_idx = {'I': 0, 'X': 1, 'Y': 2, 'Z': 3}
    for (i, j), terms in bonds_with_terms:
        for term in terms:
            ai, aj = name_to_idx[term[0]], name_to_idx[term[1]]
            ops = [I2] * N
            ops[i] = PAULIS[ai]
            ops[j] = PAULIS[aj]
            mat = ops[0]
            for o in ops[1:]:
                mat = np.kron(mat, o)
            H += mat
    return H


def build_L_vec(H, gamma_l, N):
    """Liouvillian in column-stack vec basis: |L(ρ)⟩ = L_vec |ρ⟩."""
    d = 2**N
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for l in range(N):
        if gamma_l[l] == 0:
            continue
        Zl = pauli_string([0]*l + [3] + [0]*(N-l-1))
        L = L + gamma_l[l] * (np.kron(Zl, Zl.conj()) - np.kron(Id, Id))
    return L


def build_pi_pauli(N):
    """Π in Pauli-string basis: 4^N × 4^N matrix."""
    d2 = 4**N
    Pi = np.zeros((d2, d2), dtype=complex)
    for k in range(d2):
        indices = k_to_indices(k, N)
        new_indices = tuple(PI_PERM[i] for i in indices)
        sign = 1
        for i in indices:
            sign *= PI_SIGN[i]
        new_k = indices_to_k(new_indices)
        Pi[new_k, k] = sign
    return Pi


def build_M_pauli_to_vec(N):
    """Matrix M such that vec(σ_string_α) = M[:, α] (column-stack convention)."""
    d2 = 4**N
    d = 2**N
    M = np.zeros((d*d, d2), dtype=complex)
    for k in range(d2):
        indices = k_to_indices(k, N)
        sigma = pauli_string(indices)
        M[:, k] = sigma.flatten('F')  # column-stack
    return M


def vec_to_pauli_basis(L_vec, M):
    """Transform L from vec basis to Pauli basis."""
    M_inv = np.linalg.solve(M.conj().T @ M, M.conj().T)
    return M_inv @ L_vec @ M


def palindrome_residual(L_pauli, Pi, Sigma_gamma, d2):
    """Compute Π·L·Π⁻¹ + L + 2Σγ·I. Should be zero if palindrome holds."""
    Pi_inv = np.linalg.inv(Pi)
    return Pi @ L_pauli @ Pi_inv + L_pauli + 2 * Sigma_gamma * np.eye(d2)


def block_norms_by_weight(matrix, N):
    """Compute block-norm of matrix per (XY-weight, XY-weight) pair."""
    d2 = 4**N
    weights = [total_weight_a(k_to_indices(k, N)) for k in range(d2)]
    blocks = {}
    for w_row in range(N + 1):
        rows = [k for k, w in enumerate(weights) if w == w_row]
        for w_col in range(N + 1):
            cols = [k for k, w in enumerate(weights) if w == w_col]
            block = matrix[np.ix_(rows, cols)]
            blocks[(w_row, w_col)] = np.linalg.norm(block)
    return blocks


def report(label, L_pauli, Pi, Sigma_gamma, N):
    d2 = 4**N
    R = palindrome_residual(L_pauli, Pi, Sigma_gamma, d2)
    total = np.linalg.norm(R)
    print(f"\n{label}")
    print(f"  Total residual ‖Π·L·Π⁻¹ + L + 2Σγ·I‖ = {total:.4e}")
    blocks = block_norms_by_weight(R, N)

    # Diagonal blocks (w_row = w_col)
    print(f"  Per-weight diagonal blocks:")
    for w in range(N + 1):
        size = sum(1 for k in range(d2) if total_weight_a(k_to_indices(k, N)) == w)
        print(f"    w = {w} ({size}×{size}): residual norm = {blocks[(w, w)]:.4e}")

    # Cross-weight blocks (w + w' = N, the Π-paired sectors)
    print(f"  Π-paired off-diagonal blocks (w_row + w_col = N):")
    for w in range(N + 1):
        if (w, N - w) in blocks and w <= N - w:
            print(f"    (w={w}, w={N-w}): {blocks[(w, N-w)]:.4e}")


def run_at_N(N):
    gamma = 0.1
    gamma_l = [gamma] * N
    Sigma_gamma = sum(gamma_l)

    print("=" * 78)
    print(f"V-Effect weight-sector immunity (N={N}, γ_l={gamma_l}, Σγ={Sigma_gamma:.3f})")
    print(f"Π pairs (w, N−w): N={N} gives pairs ", end="")
    pair_strs = []
    for w in range(N+1):
        if w <= N-w:
            if w == N-w:
                pair_strs.append(f"({w}↔{w})self")
            else:
                pair_strs.append(f"({w}↔{N-w})")
    print(", ".join(pair_strs))
    print("=" * 78)

    Pi = build_pi_pauli(N)
    M = build_M_pauli_to_vec(N)

    bonds = [(i, i+1) for i in range(N-1)]

    # Test 1: Heisenberg (XX + YY + ZZ on each bond)
    bonds_with_terms = [(b, ['XX', 'YY', 'ZZ']) for b in bonds]
    H = build_H_bonds(N, bonds_with_terms)
    L_vec = build_L_vec(H, gamma_l, N)
    L_pauli = vec_to_pauli_basis(L_vec, M)
    report("Test 1: H = XX+YY+ZZ on each bond (Heisenberg, both parities preserved)",
           L_pauli, Pi, Sigma_gamma, N)

    # Test 2: XY-model (XX + YY on each bond)
    bonds_with_terms = [(b, ['XX', 'YY']) for b in bonds]
    H = build_H_bonds(N, bonds_with_terms)
    L_vec = build_L_vec(H, gamma_l, N)
    L_pauli = vec_to_pauli_basis(L_vec, M)
    report("Test 2: H = XX+YY on each bond (XY model, both parities preserved)",
           L_pauli, Pi, Sigma_gamma, N)

    # Test 3: XX on (0,1) + XY on (1,2). Breaks bit_b parity.
    bonds_with_terms = [((0, 1), ['XX']), ((1, 2), ['XY'])]
    H = build_H_bonds(N, bonds_with_terms)
    L_vec = build_L_vec(H, gamma_l, N)
    L_pauli = vec_to_pauli_basis(L_vec, M)
    report("Test 3: H = XX on (0,1) + XY on (1,2) (asymmetric, breaks bit_b)",
           L_pauli, Pi, Sigma_gamma, N)

    # Test 4: XX + ZZ on each bond (XXZ Δ=large, both parities preserved)
    bonds_with_terms = [(b, ['XX', 'ZZ']) for b in bonds]
    H = build_H_bonds(N, bonds_with_terms)
    L_vec = build_L_vec(H, gamma_l, N)
    L_pauli = vec_to_pauli_basis(L_vec, M)
    report("Test 4: H = XX+ZZ on each bond (both parities preserved, no YY)",
           L_pauli, Pi, Sigma_gamma, N)

    # Test 5: ZZ on (0,1) + XY on (1,2). Single XY mixed in.
    bonds_with_terms = [((0, 1), ['ZZ']), ((1, 2), ['XY'])]
    H = build_H_bonds(N, bonds_with_terms)
    L_vec = build_L_vec(H, gamma_l, N)
    L_pauli = vec_to_pauli_basis(L_vec, M)
    report("Test 5: H = ZZ on (0,1) + XY on (1,2) (only mixed bond breaks parity)",
           L_pauli, Pi, Sigma_gamma, N)

    print("\n" + "=" * 78)
    print("Reading: where the residual is exactly 0, the palindrome holds in")
    print("that block. Where it is non-zero, that sector breaks. Compare across")
    print("Hamiltonian forms to see WHICH ones break which sectors.")
    print("=" * 78)


def main():
    for N in [3, 4]:
        run_at_N(N)
        print()


if __name__ == "__main__":
    main()
