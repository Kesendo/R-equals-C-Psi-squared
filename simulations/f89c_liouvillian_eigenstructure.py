"""F89c structural lemma verification: Liouvillian eigenstructure for path-k blocks.

For a (k+1)-qubit block under H_B = J·Σ_b (X_b X_{b+1} + Y_b Y_{b+1}) + uniform
Z-dephasing γ per site, the full Lindbladian L_super has dimension (2^(k+1))².

Claim (F89c): the eigenvalues organize as λ_mode = -Γ_mode·γ/2 - i·ω_mode·J where
  Γ_mode ∈ {2, 6} per coherence (corresponding to overlap-count of the basis-state
    pair: 1 differing site → 2γ, 3 differing sites → 6γ)
  ω_mode = (E_μ - E_ν)/J for H_B eigenstates μ, ν

S(t) = Σ_l 2|(ρ_l)_{0,1}|² has rates 2·Re(λ) on |.|² level → effective Γ ∈ {4, 8, 12}γ
in the |amplitude|² envelope.

Test: compute L_super for path-1 (2 qubits) and path-2 (3 qubits), eigendecompose,
list distinct decay rates and frequencies. Confirm structural prediction.
"""
from __future__ import annotations

import sys

import numpy as np
from scipy.linalg import eig

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

J = 1.0  # for cleaner ratios
GAMMA = 1.0  # for cleaner ratios


def pauli_at(P: np.ndarray, site: int, n_qubits: int) -> np.ndarray:
    op = np.array([[1.0]], dtype=complex)
    I2 = np.eye(2, dtype=complex)
    for q in range(n_qubits):
        op = np.kron(op, P if q == site else I2)
    return op


X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def build_path_k_H_B(n_qubits: int, j: float = J) -> np.ndarray:
    """H_B for path-(n_qubits-1) on n_qubits."""
    n = n_qubits
    H = np.zeros((2**n, 2**n), dtype=complex)
    for b in range(n - 1):
        H += j * (
            pauli_at(X, b, n) @ pauli_at(X, b + 1, n)
            + pauli_at(Y, b, n) @ pauli_at(Y, b + 1, n)
        )
    return H


def build_liouvillian(H: np.ndarray, gamma: float, n_qubits: int) -> np.ndarray:
    """L_super on vec(ρ) for H + uniform Z-dephasing γ on each site.

    L[ρ] = -i[H, ρ] + γ·Σ_l (Z_l ρ Z_l - ρ)
    Vec convention: vec(ABC) = (C^T ⊗ A) vec(B).
    """
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for l in range(n_qubits):
        Zl = pauli_at(Z, l, n_qubits)
        L += gamma * (np.kron(Zl.T, Zl) - np.kron(Id, Id))
    return L


def analyze_block(n_qubits: int, label: str) -> None:
    print(f"\n## {label} ({n_qubits} qubits, path-{n_qubits-1} topology)\n")
    H = build_path_k_H_B(n_qubits)
    L = build_liouvillian(H, GAMMA, n_qubits)

    # H_B eigenvalues (SE sector only, single-excitation)
    # SE basis indices: those with popcount=1
    se_indices = [i for i in range(2**n_qubits) if bin(i).count("1") == 1]
    H_se = H[np.ix_(se_indices, se_indices)]
    H_se_eigs = np.linalg.eigvalsh(H_se)
    print(f"  H_B SE eigenvalues / J: {sorted(H_se_eigs.real)}")

    # DE eigenvalues
    de_indices = [i for i in range(2**n_qubits) if bin(i).count("1") == 2]
    if len(de_indices) > 0:
        H_de = H[np.ix_(de_indices, de_indices)]
        H_de_eigs = np.linalg.eigvalsh(H_de)
        print(f"  H_B DE eigenvalues / J: {sorted(H_de_eigs.real)}")

    # L_super full eigenvalues
    L_eigs, _ = eig(L)
    real_parts = -L_eigs.real  # decay rate (positive)
    imag_parts = L_eigs.imag  # frequency

    # Round and group
    distinct_rates = sorted(set(np.round(real_parts, 6)))
    print(f"  L_super dim: {L.shape}")
    print(f"  Distinct decay rates Γ/γ (-Re(λ)/γ, all should be 0, 2, 4, 6, 8, ...): {distinct_rates}")

    # Group eigenvalues by rate
    print(f"\n  Eigenvalue table (rate, frequency / J):")
    for rate in distinct_rates:
        mask = np.abs(real_parts - rate) < 1e-6
        freqs = sorted(set(np.round(imag_parts[mask], 6)))
        print(f"    Γ/γ = {rate:6.3f}: {len(freqs)} distinct freqs ω/J: {freqs}")


# Run for path-1 (2 qubits) - should match F89 all-isolated structure
analyze_block(2, "Path-1 / 2-qubit block (F89 all-isolated case)")

# Run for path-2 (3 qubits) - should show richer structure with rate 6γ appearing
analyze_block(3, "Path-2 / 3-qubit block (F89c subclass)")

# Run for path-3 (4 qubits)
analyze_block(4, "Path-3 / 4-qubit block (F89c subclass)")
