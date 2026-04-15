"""
Primordial Qubit Follow-up: Probe 1, [L, Pi^2] = 0 at N >= 3?
===============================================================

At N=2, [L, Pi^2] = 0 exactly (primordial_bit_a_bit_b.py). Does this
hold at N=3 and beyond?

Pi^2 in the Pauli basis: diagonal with entries (-1)^{w_YZ_total},
where w_YZ = number of Y or Z factors in the Pauli string.

If [L, Pi^2] = 0 at all N: the C2xC2 sector decomposition is universal.
If it breaks at N >= 3: structural confirmation is N=2-specific.

Setup: chain S-...-B with XX+YY coupling, gamma only on B (outermost).

Date: 2026-04-15
Authors: Tom and Claude (Code)
"""

import numpy as np
from itertools import product as iproduct
from pathlib import Path
import sys

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
paulis = [I2, X, Y, Z]
pauli_names = ['I', 'X', 'Y', 'Z']

n_xy_site = [0, 1, 1, 0]  # bit a per site
w_yz_site = [0, 0, 1, 1]  # bit b per site


def kron(*args):
    out = args[0]
    for a in args[1:]:
        out = np.kron(out, a)
    return out


def liouvillian_comp(H, jumps):
    """Build Liouvillian in computational (column-vec) basis."""
    d = H.shape[0]
    Idd = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Idd, H) - np.kron(H.T, Idd))
    for Lk in jumps:
        LdL = Lk.conj().T @ Lk
        L += (np.kron(Lk.conj(), Lk)
              - 0.5 * (np.kron(Idd, LdL) + np.kron(LdL.T, Idd)))
    return L


def build_pi2_comp(N):
    """Build Pi^2 superoperator in computational (column-vec) basis.

    Pi^2 acts on a Pauli string by multiplying by (-1)^{w_YZ_total}.
    In the computational basis, Pi^2 acts on rho (density matrix) by
    conjugation with the Hilbert-space operator Pi^2_H.

    Per site: Pi^2 maps I->I, X->X, Y->-Y, Z->-Z.
    In the Hilbert space: Pi^2_H = X on each qubit (since X*I*X=I,
    X*X*X=X, X*Y*X=-Y, X*Z*X=-Z).

    Pi^2_super(rho) = Pi^2_H rho Pi^2_H^dagger.
    In column-vec: Pi^2_super = conj(Pi^2_H) (x) Pi^2_H.
    Since Pi^2_H = X^{(x)N} is real: Pi^2_super = Pi^2_H (x) Pi^2_H.
    """
    # Pi^2 in Hilbert space = X on each qubit
    Pi2_H = np.eye(1, dtype=complex)
    for _ in range(N):
        Pi2_H = np.kron(Pi2_H, X)

    # Superoperator: conj(Pi2_H) (x) Pi2_H
    Pi2_super = np.kron(Pi2_H.conj(), Pi2_H)
    return Pi2_super


def build_chain_hamiltonian(N, J):
    """Build XX+YY chain Hamiltonian for N qubits."""
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for bond in range(N - 1):
        for P in [X, Y]:
            ops = [I2] * N
            ops[bond] = P
            ops[bond + 1] = P
            term = ops[0]
            for op in ops[1:]:
                term = np.kron(term, op)
            H += J * 0.5 * term
    return H


results_dir = Path("simulations/results/primordial_bit_a_bit_b")

print("=" * 72)
print("Probe 1: [L, Pi^2] = 0 at N >= 3?")
print("=" * 72)

log_lines = []

for N in [2, 3, 4, 5]:
    dim_H = 2**N      # Hilbert space dim
    dim_L = 4**N      # Liouville space dim

    if dim_L > 100000:  # skip if > ~100k (memory)
        print(f"\nN={N}: dim={dim_L}, skipping (memory)")
        continue

    print(f"\nN={N}: dim_H={dim_H}, dim_L={dim_L}")
    print("-" * 40)

    # Hamiltonian: XX+YY chain
    H = build_chain_hamiltonian(N, J=1.0)

    # Jump: gamma only on last qubit (B)
    gamma_B = 0.1
    ops_jump = [I2] * N
    ops_jump[N - 1] = Z
    L_jump_op = np.sqrt(gamma_B) * ops_jump[0]
    for op in ops_jump[1:]:
        L_jump_op = np.kron(L_jump_op, op)

    # Build L
    L = liouvillian_comp(H, [L_jump_op])

    # Build Pi^2
    Pi2 = build_pi2_comp(N)

    # Test commutator
    comm = L @ Pi2 - Pi2 @ L
    comm_norm = np.linalg.norm(comm)
    L_norm = np.linalg.norm(L)
    rel_comm = comm_norm / L_norm if L_norm > 0 else 0

    commutes = comm_norm < 1e-10
    result = "COMMUTES" if commutes else f"BREAKS (rel={rel_comm:.2e})"

    line = f"  ||[L, Pi^2]|| = {comm_norm:.6e}  (relative: {rel_comm:.6e})  {result}"
    print(line)
    log_lines.append(f"N={N}: {line}")

    if not commutes:
        print(f"  [L, Pi^2] != 0 at N={N}.")
        # Characterize: which part of L causes the breaking?
        # L = L_H + L_D. Test each.
        L_H = liouvillian_comp(H, [])
        L_D = L - L_H

        comm_H = np.linalg.norm(L_H @ Pi2 - Pi2 @ L_H)
        comm_D = np.linalg.norm(L_D @ Pi2 - Pi2 @ L_D)
        print(f"  ||[L_H, Pi^2]|| = {comm_H:.6e}")
        print(f"  ||[L_D, Pi^2]|| = {comm_D:.6e}")
        log_lines.append(f"  L_H: {comm_H:.6e}, L_D: {comm_D:.6e}")
        continue

    # If commuting: sector decomposition
    # Pi^2 eigenvalues: +1 (even) and -1 (odd)
    Pi2_eigvals, Pi2_eigvecs = np.linalg.eigh(Pi2.real)  # Pi2 is real symmetric

    even_mask = Pi2_eigvals > 0.5
    odd_mask = Pi2_eigvals < -0.5
    n_even = np.sum(even_mask)
    n_odd = np.sum(odd_mask)

    print(f"  Even sector: {n_even} modes, Odd sector: {n_odd} modes")
    log_lines.append(f"  Even: {n_even}, Odd: {n_odd}")

    # Project L into sectors
    V_even = Pi2_eigvecs[:, even_mask]
    V_odd = Pi2_eigvecs[:, odd_mask]

    L_even = V_even.T.conj() @ L @ V_even
    L_odd = V_odd.T.conj() @ L @ V_odd

    # Cross-sector coupling (should be zero)
    cross = np.linalg.norm(V_even.T.conj() @ L @ V_odd)
    print(f"  Cross-sector coupling: {cross:.6e}")

    # Eigenvalues per sector
    ev_even = np.linalg.eigvals(L_even)
    ev_odd = np.linalg.eigvals(L_odd)

    # Count per Re(lambda) class
    def classify_modes(eigvals, gamma):
        classes = {'conserved': 0, 'correlation': 0, 'mirror': 0}
        for e in eigvals:
            if abs(e.real) < 1e-6:
                classes['conserved'] += 1
            elif abs(e.real + 2 * gamma) < 1e-5:
                classes['correlation'] += 1
            else:
                classes['mirror'] += 1
        return classes

    cls_even = classify_modes(ev_even, gamma_B)
    cls_odd = classify_modes(ev_odd, gamma_B)

    print(f"  Even: conserved={cls_even['conserved']}, mirror={cls_even['mirror']}, "
          f"correlation={cls_even['correlation']}")
    print(f"  Odd:  conserved={cls_odd['conserved']}, mirror={cls_odd['mirror']}, "
          f"correlation={cls_odd['correlation']}")
    log_lines.append(f"  Even: {cls_even}")
    log_lines.append(f"  Odd:  {cls_odd}")

    # Total check
    total = sum(cls_even.values()) + sum(cls_odd.values())
    print(f"  Total: {total} (expected {dim_L})")

# Also test with Heisenberg XXX (adding ZZ) and dephasing on ALL sites
print("\n" + "=" * 72)
print("Additional: N=3 with dephasing on ALL sites (uniform gamma)")
print("=" * 72)

N = 3
H_xxx = build_chain_hamiltonian(N, J=1.0)
# Add ZZ coupling
for bond in range(N - 1):
    ops = [I2] * N
    ops[bond] = Z
    ops[bond + 1] = Z
    term = ops[0]
    for op in ops[1:]:
        term = np.kron(term, op)
    H_xxx += 0.5 * term

# Dephasing on all sites
jumps_all = []
for site in range(N):
    ops = [I2] * N
    ops[site] = Z
    jump_op = np.sqrt(gamma_B) * ops[0]
    for op in ops[1:]:
        jump_op = np.kron(jump_op, op)
    jumps_all.append(jump_op)

L_all = liouvillian_comp(H_xxx, jumps_all)
Pi2 = build_pi2_comp(N)

comm_norm = np.linalg.norm(L_all @ Pi2 - Pi2 @ L_all)
L_norm = np.linalg.norm(L_all)
print(f"  N=3 XXX, uniform gamma: ||[L, Pi^2]|| = {comm_norm:.6e} "
      f"(relative: {comm_norm/L_norm:.6e})")
log_lines.append(f"\nN=3 XXX uniform: ||[L, Pi^2]|| = {comm_norm:.6e}")

# Save
with open(results_dir / 'pi2_commutation_scaling.txt', 'w', encoding='utf-8') as f:
    f.write("Probe 1: [L, Pi^2] = 0 scaling test\n")
    f.write("=" * 72 + "\n\n")
    f.write('\n'.join(log_lines))

print(f"\nResults saved to {results_dir / 'pi2_commutation_scaling.txt'}")
