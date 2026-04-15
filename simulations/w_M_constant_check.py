"""
Primordial Qubit Follow-up: Probe 2, w_M = 0.75 explanation
=============================================================

Observation from Inside-Outside Probe 4: at N=3 chain S-M-B, the
middle qubit M has constant Pauli weight w_M = 0.75 across all
J/gamma values, while w_S and w_B vary.

Test whether this is:
  (a) Chain-N=3 specific (changes with topology)
  (b) Universal (holds across topologies and N)
  (c) A numerical artifact

Date: 2026-04-15
Authors: Tom and Claude (Code)
"""

import numpy as np
from pathlib import Path
from itertools import product as iproduct
import sys

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
paulis = [I2, X, Y, Z]


def kron(*args):
    out = args[0]
    for a in args[1:]:
        out = np.kron(out, a)
    return out


def liouvillian(H, jumps):
    d = H.shape[0]
    Idd = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Idd, H) - np.kron(H.T, Idd))
    for Lk in jumps:
        LdL = Lk.conj().T @ Lk
        L += (np.kron(Lk.conj(), Lk)
              - 0.5 * (np.kron(Idd, LdL) + np.kron(LdL.T, Idd)))
    return L


def build_pauli_basis(n_qubits):
    ops = []
    labels_per_site = []
    for idx in iproduct(range(4), repeat=n_qubits):
        op = paulis[idx[0]]
        for k in idx[1:]:
            op = np.kron(op, paulis[k])
        ops.append(op)
        labels_per_site.append(idx)
    return ops, labels_per_site


def site_weight(M, site, basis_ops, labels, d, n_qubits):
    """Fraction of Pauli expansion weight that has non-I at given site."""
    coeffs = np.zeros(len(basis_ops), dtype=complex)
    for k, P in enumerate(basis_ops):
        coeffs[k] = np.trace(P.conj().T @ M) / d
    weights = np.abs(coeffs) ** 2
    total = weights.sum()
    if total < 1e-14:
        return 0.0
    w = sum(weights[k] for k in range(len(basis_ops)) if labels[k][site] != 0)
    return w / total


def analyze_topology(name, H, jumps, n_qubits, gamma_B, J_values):
    """Analyze slow-mode site weights across J values."""
    d = 2**n_qubits
    dim_L = d**2
    basis, labels = build_pauli_basis(n_qubits)

    print(f"\n  {name}:")
    print(f"  {'J':>8}", end="")
    for s in range(n_qubits):
        print(f"  {'w_'+str(s):>8}", end="")
    print()

    for J in J_values:
        H_scaled = J * H  # H is normalized to J=1
        jump_ops = [np.sqrt(gamma_B) * j for j in jumps]
        L = liouvillian(H_scaled, jump_ops)
        eigvals, eigvecs = np.linalg.eig(L)

        # Find slowest nonzero modes
        mode_data = []
        for k in range(dim_L):
            re = eigvals[k].real
            if abs(re) > 1e-8:
                M = eigvecs[:, k].reshape((d, d), order='F')
                ws = [site_weight(M, s, basis, labels, d, n_qubits)
                      for s in range(n_qubits)]
                mode_data.append((abs(re), ws))

        mode_data.sort(key=lambda x: x[0])
        # Average over the 4 slowest modes
        n_avg = min(4, len(mode_data))
        if n_avg == 0:
            continue
        avg_ws = [np.mean([mode_data[i][1][s] for i in range(n_avg)])
                  for s in range(n_qubits)]

        print(f"  {J:8.2f}", end="")
        for w in avg_ws:
            print(f"  {w:8.4f}", end="")
        print()

    return


gamma_B = 0.1
J_values = [0.1, 1.0, 10.0]

print("=" * 72)
print("Probe 2: w_M = 0.75 constant check across topologies")
print("=" * 72)

# Topology 1: N=3 linear chain S-M-B, gamma on B
print("\n--- N=3 Linear Chain S-M-B, gamma on B ---")
H_chain3 = 0.5 * (kron(X, X, I2) + kron(Y, Y, I2)
                  + kron(I2, X, X) + kron(I2, Y, Y))
jumps_chain3 = [kron(I2, I2, Z)]
analyze_topology("chain S-M-B (gamma on B)", H_chain3, jumps_chain3, 3, gamma_B, J_values)

# Topology 2: N=3 linear chain, gamma on S (other end)
print("\n--- N=3 Linear Chain S-M-B, gamma on S ---")
jumps_chain3_S = [kron(Z, I2, I2)]
analyze_topology("chain S-M-B (gamma on S)", H_chain3, jumps_chain3_S, 3, gamma_B, J_values)

# Topology 3: N=3 linear chain, gamma on M (middle)
print("\n--- N=3 Linear Chain S-M-B, gamma on M ---")
jumps_chain3_M = [kron(I2, Z, I2)]
analyze_topology("chain S-M-B (gamma on M)", H_chain3, jumps_chain3_M, 3, gamma_B, J_values)

# Topology 4: N=3 star (M coupled to both S and B)
print("\n--- N=3 Star (M central), gamma on B ---")
H_star3 = 0.5 * (kron(X, X, I2) + kron(Y, Y, I2)   # S-M bond
                + kron(I2, X, X) + kron(I2, Y, Y))   # M-B bond (same as chain)
# Star = chain for 3 qubits
analyze_topology("star (= chain at N=3, gamma on B)", H_star3, jumps_chain3, 3, gamma_B, J_values)

# Topology 5: N=4 linear chain S-M1-M2-B, gamma on B
print("\n--- N=4 Linear Chain S-M1-M2-B, gamma on B ---")
H_chain4 = 0.5 * (kron(X, X, I2, I2) + kron(Y, Y, I2, I2)
                  + kron(I2, X, X, I2) + kron(I2, Y, Y, I2)
                  + kron(I2, I2, X, X) + kron(I2, I2, Y, Y))
jumps_chain4 = [kron(I2, I2, I2, Z)]
analyze_topology("chain S-M1-M2-B (gamma on B)", H_chain4, jumps_chain4, 4, gamma_B, J_values)

# Topology 6: N=3 ring S-M-B-S, gamma on B
print("\n--- N=3 Ring S-M-B-S, gamma on B ---")
H_ring3 = H_chain3 + 0.5 * (kron(X, I2, X) + kron(Y, I2, Y))  # add S-B bond
analyze_topology("ring S-M-B-S (gamma on B)", H_ring3, jumps_chain3, 3, gamma_B, J_values)

# =========================================================================
# VERDICT
# =========================================================================
print("\n" + "=" * 72)
print("VERDICT: w_M = 0.75")
print("=" * 72)

results_dir = Path("simulations/results/primordial_bit_a_bit_b")
with open(results_dir / 'w_M_check.txt', 'w', encoding='utf-8') as f:
    f.write("Probe 2: w_M = 0.75 constant check\n")
    f.write("See console output for full topology comparison.\n")

print("Results saved.")
