"""
Primordial Qubit TFD: Probe 1, KMS / Detailed Balance Test
============================================================

Question: is there a state rho_* in the Lindblad steady-state manifold
such that the centered Liouvillian L_c satisfies quantum detailed balance
(QDB) with respect to rho_*?

QDB condition (sigma-DBC, Alicki 1976):
    Tr(sigma L(A)^dagger B) = Tr(sigma A^dagger L(B))   for all A, B

In the Pauli basis with metric G_ab = Tr(sigma P_a^dagger P_b) / d:
    L^T G = G L

For sigma = I/d: G = I, so QDB reduces to L = L^T (L symmetric).
This fails because L_H is antisymmetric.

For non-trivial sigma in the steady-state manifold {II, ZZ, ZI+IZ}:
G is non-trivial, and QDB might be satisfiable.

If QDB holds for some sigma: L_c is self-adjoint w.r.t. the GNS inner
product defined by sigma, enabling a TFD-like construction.

If QDB fails for ALL sigma in the manifold: the TFD route via standard
modular theory is blocked.

Date: 2026-04-15
Authors: Tom and Claude (Code)
"""

import numpy as np
from itertools import product as iproduct
from pathlib import Path
import sys, os

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, sx, sy, sz]


def kron(*args):
    out = args[0]
    for a in args[1:]:
        out = np.kron(out, a)
    return out


def build_pauli_matrices_2q():
    """Build 16 two-qubit Pauli string matrices."""
    pmats = []
    labels = []
    for i in range(4):
        for j in range(4):
            pmats.append(np.kron(PAULIS[i], PAULIS[j]))
            labels.append(['I','X','Y','Z'][i] + ['I','X','Y','Z'][j])
    return pmats, labels


def liouvillian_pauli_basis(H, jumps, pmats, d):
    """Build Liouvillian in the Pauli basis.
    L_ab = (1/d) Tr(P_a^dag L(P_b))
    """
    num = len(pmats)
    L = np.zeros((num, num), dtype=complex)

    for b in range(num):
        # Hamiltonian part: -i[H, P_b]
        comm = -1j * (H @ pmats[b] - pmats[b] @ H)

        # Dissipator part
        diss = np.zeros_like(pmats[0])
        for Lk in jumps:
            LdL = Lk.conj().T @ Lk
            diss += (Lk @ pmats[b] @ Lk.conj().T
                     - 0.5 * (LdL @ pmats[b] + pmats[b] @ LdL))

        total = comm + diss
        for a in range(num):
            L[a, b] = np.trace(pmats[a].conj().T @ total) / d

    return L


# =========================================================================
# Setup: N=2 Heisenberg chain with Z-dephasing
# =========================================================================
N = 2
d = 2**N  # = 4
num = 4**N  # = 16

J_coup = 1.0
gamma = 0.1
Sigma_gamma = 2 * gamma  # two sites, uniform

pmats, labels = build_pauli_matrices_2q()

# Hamiltonian: Heisenberg XXX
H = J_coup * (kron(sx, sx) + kron(sy, sy) + kron(sz, sz))

# Jump operators: sqrt(gamma) * Z on each site
jumps = [np.sqrt(gamma) * kron(sz, I2),
         np.sqrt(gamma) * kron(I2, sz)]

# Build L in Pauli basis
L = liouvillian_pauli_basis(H, jumps, pmats, d)

# Centered Liouvillian
L_c = L + Sigma_gamma * np.eye(num)

print("=" * 72)
print("Probe 1: KMS / Quantum Detailed Balance Test")
print("=" * 72)
print(f"N=2, J={J_coup}, gamma={gamma}, Sigma_gamma={Sigma_gamma}")
print(f"Liouvillian dim: {num}x{num}")
print()

# =========================================================================
# Step 1: Identify steady-state manifold
# =========================================================================
print("Step 1: Steady-state manifold of L")
print("-" * 72)

eigvals_L, eigvecs_L = np.linalg.eig(L)
steady_indices = [i for i in range(num) if abs(eigvals_L[i]) < 1e-10]
print(f"  Steady-state dimension: {len(steady_indices)}")

# The steady states in operator form
steady_ops = []
for idx in steady_indices:
    vec = eigvecs_L[:, idx].real  # steady states are real in Pauli basis
    op = sum(c * P for c, P in zip(vec, pmats))
    steady_ops.append((vec, op, labels))
    # Identify which Pauli strings contribute
    nonzero = [(labels[k], vec[k]) for k in range(num) if abs(vec[k]) > 1e-10]
    print(f"  Mode {idx}: {nonzero}")

# =========================================================================
# Step 2: Build the metric G(sigma) and test QDB for parametric sigma
# =========================================================================
print()
print("Step 2: Quantum Detailed Balance test over steady-state manifold")
print("-" * 72)

# Parametrize sigma = (1/4)(II + alpha*(ZI+IZ) + beta*ZZ)
# Positivity constraint: sigma must have non-negative eigenvalues
# sigma eigenvalues for a state (1/4)(II + alpha*(ZI+IZ) + beta*ZZ):
# Diagonal in computational basis: (1+alpha+beta, 1-alpha+beta, etc.)/4
# Need: 1+2*alpha+beta >= 0, 1-2*alpha+beta >= 0, 1-beta >= 0 (doubly degenerate)
# So: beta <= 1, |alpha| <= (1+beta)/2

def make_sigma(alpha, beta, pmats):
    """Build sigma = (1/4)(II + alpha*(ZI+IZ) + beta*ZZ)."""
    II = pmats[0]   # index 0: II
    ZI = pmats[12]  # index 12: ZI
    IZ = pmats[3]   # index 3: IZ
    ZZ = pmats[15]  # index 15: ZZ
    return (II + alpha * (ZI + IZ) + beta * ZZ) / 4


def qdb_violation(L_mat, sigma, pmats, d):
    """Compute QDB violation: ||L^T G - G L|| / ||G L||.
    G_ab = (1/d) Tr(sigma P_a^dagger P_b)
    """
    num = len(pmats)
    G = np.zeros((num, num), dtype=complex)
    for a in range(num):
        for b in range(num):
            G[a, b] = np.trace(sigma @ pmats[a].conj().T @ pmats[b]) / d

    LtG = L_mat.T @ G
    GL = G @ L_mat
    diff = LtG - GL
    norm_GL = np.linalg.norm(GL)
    if norm_GL < 1e-15:
        return 1.0, G
    return np.linalg.norm(diff) / norm_GL, G


# Sweep over (alpha, beta) grid
print(f"\n  Sweeping (alpha, beta) for sigma = (1/4)(II + alpha*(ZI+IZ) + beta*ZZ)")
print(f"  Testing QDB for L and L_c separately\n")

best_L = {'violation': 1.0, 'alpha': 0, 'beta': 0}
best_Lc = {'violation': 1.0, 'alpha': 0, 'beta': 0}

alpha_range = np.linspace(-0.9, 0.9, 37)
beta_range = np.linspace(-0.9, 0.9, 37)

for alpha in alpha_range:
    for beta in beta_range:
        # Check positivity
        if beta > 1 - 1e-8:
            continue
        if abs(alpha) > (1 + beta) / 2 - 1e-8:
            continue

        sigma = make_sigma(alpha, beta, pmats)

        # Check positive semidefinite
        eigs = np.linalg.eigvalsh(sigma)
        if np.min(eigs) < -1e-10:
            continue

        viol_L, _ = qdb_violation(L, sigma, pmats, d)
        viol_Lc, _ = qdb_violation(L_c, sigma, pmats, d)

        if viol_L < best_L['violation']:
            best_L = {'violation': viol_L, 'alpha': alpha, 'beta': beta}
        if viol_Lc < best_Lc['violation']:
            best_Lc = {'violation': viol_Lc, 'alpha': alpha, 'beta': beta}

print(f"  Best QDB violation for L:")
print(f"    alpha={best_L['alpha']:.4f}, beta={best_L['beta']:.4f}, "
      f"violation={best_L['violation']:.6e}")

print(f"\n  Best QDB violation for L_c:")
print(f"    alpha={best_Lc['alpha']:.4f}, beta={best_Lc['beta']:.4f}, "
      f"violation={best_Lc['violation']:.6e}")

# =========================================================================
# Step 3: Fine-grid search around the best point
# =========================================================================
print()
print("Step 3: Fine-grid search around best L_c point")
print("-" * 72)

a0, b0 = best_Lc['alpha'], best_Lc['beta']
delta = 0.1
alpha_fine = np.linspace(max(-0.99, a0 - delta), min(0.99, a0 + delta), 51)
beta_fine = np.linspace(max(-0.99, b0 - delta), min(0.99, b0 + delta), 51)

for alpha in alpha_fine:
    for beta in beta_fine:
        if beta > 1 - 1e-8 or abs(alpha) > (1 + beta) / 2 - 1e-8:
            continue
        sigma = make_sigma(alpha, beta, pmats)
        eigs = np.linalg.eigvalsh(sigma)
        if np.min(eigs) < -1e-10:
            continue

        viol_Lc, _ = qdb_violation(L_c, sigma, pmats, d)
        if viol_Lc < best_Lc['violation']:
            best_Lc = {'violation': viol_Lc, 'alpha': alpha, 'beta': beta}

print(f"  Refined best QDB violation for L_c:")
print(f"    alpha={best_Lc['alpha']:.6f}, beta={best_Lc['beta']:.6f}, "
      f"violation={best_Lc['violation']:.6e}")

# =========================================================================
# Step 4: Test at Sigma_gamma = 0 (the mirror)
# =========================================================================
print()
print("Step 4: QDB test at Sigma_gamma = 0 (L_c = L_H)")
print("-" * 72)

# At gamma = 0, L = L_H (pure Hamiltonian). L_H is antisymmetric.
# QDB for L_H w.r.t. I/d: L_H^T = -L_H, so L_H^T G = -G L_H
# This means L_H is ANTI-self-adjoint, not self-adjoint.
# But at beta=0 (infinite T), KMS is trivially satisfied.

# Build L_H only (no dissipation)
L_H = liouvillian_pauli_basis(H, [], pmats, d)

# For sigma = I/d: G = I
viol_LH_trivial, _ = qdb_violation(L_H, np.eye(d) / d, pmats, d)
print(f"  L_H QDB violation at sigma=I/d: {viol_LH_trivial:.6e}")
print(f"  (Expected: nonzero, because L_H is antisymmetric)")

# For L_H, test if L_H^T G + G L_H = 0 (anti-self-adjoint)
G_trivial = np.eye(num)
anti_sa = np.linalg.norm(L_H.T @ G_trivial + G_trivial @ L_H)
print(f"  ||L_H^T + L_H|| = {anti_sa:.6e} (anti-self-adjoint test)")

# Check: is L_H antisymmetric in Pauli basis?
print(f"  ||L_H + L_H^T|| = {np.linalg.norm(L_H + L_H.T):.6e}")
print(f"  ||L_H - L_H^T|| = {np.linalg.norm(L_H - L_H.T):.6e} (should be ~ 2*||L_H||)")

# =========================================================================
# Step 5: Modular Hamiltonian test
# =========================================================================
print()
print("Step 5: Test if L_c could be a modular Hamiltonian")
print("-" * 72)

# For a thermal state sigma_beta = exp(-beta*H) / Z:
# The modular Hamiltonian is K = beta * L_H (acting on operator space)
# K is anti-self-adjoint w.r.t. sigma_beta inner product (generates unitaries)

# Test: at what beta does sigma_beta sit in the steady-state manifold?
# sigma_beta = exp(-beta*H) / Z is diagonal in H-eigenbasis
# For Heisenberg H on 2 qubits: eigenvalues are {J, J, J, -3J} (triplet+singlet)

H_eigvals = np.linalg.eigvalsh(H)
print(f"  H eigenvalues: {H_eigvals}")

# sigma_beta in Pauli basis at various beta
for beta in [0.0, 1.0, 5.0, 10.0, 1/(2*gamma)]:
    exp_bH = np.diag(np.exp(-beta * H_eigvals))
    # Transform back to computational basis
    _, H_eigvecs = np.linalg.eigh(H)
    sigma_b = H_eigvecs @ exp_bH @ H_eigvecs.conj().T
    sigma_b /= np.trace(sigma_b)

    # Check if sigma_b is in the steady-state manifold (commutes with Z_i)
    comm_Z1 = np.linalg.norm(sigma_b @ kron(sz, I2) - kron(sz, I2) @ sigma_b)
    comm_Z2 = np.linalg.norm(sigma_b @ kron(I2, sz) - kron(I2, sz) @ sigma_b)

    # Is it a steady state of L?
    vec_sigma = np.zeros(num, dtype=complex)
    for k in range(num):
        vec_sigma[k] = np.trace(pmats[k].conj().T @ sigma_b) / d
    L_sigma = L @ vec_sigma
    is_steady = np.linalg.norm(L_sigma) < 1e-8

    beta_label = f"1/Sg={beta:.1f}" if abs(beta - 1/(2*gamma)) < 0.01 else f"{beta:.1f}"
    print(f"  beta={beta_label:>10}: [Z1,sigma]={comm_Z1:.6e}, "
          f"[Z2,sigma]={comm_Z2:.6e}, steady={is_steady}")

# =========================================================================
# VERDICT
# =========================================================================
print()
print("=" * 72)
print("VERDICT: KMS / Quantum Detailed Balance")
print("=" * 72)

verdict = []
if best_Lc['violation'] < 1e-6:
    verdict.append(f"QDB SATISFIED for L_c at alpha={best_Lc['alpha']:.4f}, "
                   f"beta={best_Lc['beta']:.4f}")
    verdict.append("TFD construction has a starting point.")
elif best_Lc['violation'] < 0.1:
    verdict.append(f"QDB approximately satisfied (violation={best_Lc['violation']:.4e})")
    verdict.append("Weak TFD structure may exist.")
else:
    verdict.append(f"QDB NOT satisfied for any sigma in the steady-state manifold")
    verdict.append(f"Best violation: {best_Lc['violation']:.4e} at "
                   f"alpha={best_Lc['alpha']:.4f}, beta={best_Lc['beta']:.4f}")
    verdict.append("The naive TFD route via standard modular theory is blocked.")
    verdict.append("")
    verdict.append("L_c contains both Hamiltonian (antisymmetric) and dissipative")
    verdict.append("(symmetric, diagonal) parts. No single state makes both")
    verdict.append("self-adjoint simultaneously. This is the algebraic obstruction.")

for line in verdict:
    print(line)

# Save
results_dir = Path("simulations/results/primordial_qubit_kms")
with open(results_dir / 'kms_test_results.txt', 'w', encoding='utf-8') as f:
    f.write("Probe 1: KMS / Quantum Detailed Balance Test\n")
    f.write("=" * 72 + "\n\n")
    f.write(f"Best QDB violation for L: {best_L['violation']:.6e}\n")
    f.write(f"  at alpha={best_L['alpha']:.4f}, beta={best_L['beta']:.4f}\n\n")
    f.write(f"Best QDB violation for L_c: {best_Lc['violation']:.6e}\n")
    f.write(f"  at alpha={best_Lc['alpha']:.6f}, beta={best_Lc['beta']:.6f}\n\n")
    f.write("Verdict:\n")
    for line in verdict:
        f.write(line + "\n")

print(f"\nResults saved to {results_dir / 'kms_test_results.txt'}")
