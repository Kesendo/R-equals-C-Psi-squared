"""Shared infrastructure for F89 path-k closed-form derivation scripts.

Generic n_block-parametric building blocks for the per-block Lindbladian
super-operator + ρ_cc partial-trace machinery used by:
  - _f89_path2_multi_exp_derive.py
  - _f89_path3_multi_exp_derive.py
  - _f89_pathk_survey.py
  - _f89_mixed_topology_additive.py
  - _f89_vac_se_parseval_closed.py
  - _f89_path6_full_chain.py

See `experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md` for derivation context.
"""

from __future__ import annotations

import numpy as np

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_at(P: np.ndarray, site: int, n_qubits: int) -> np.ndarray:
    """Single-site Kronecker placement: I⊗...⊗P_site⊗...⊗I on n_qubits."""
    op = np.array([[1.0]], dtype=complex)
    for q in range(n_qubits):
        op = np.kron(op, P if q == site else I2)
    return op


def build_block_H(J: float, n_block: int) -> np.ndarray:
    """XY Hamiltonian on n_block sites: H_B = J·Σ_b (X_b X_{b+1} + Y_b Y_{b+1})."""
    d = 2**n_block
    H = np.zeros((d, d), dtype=complex)
    for b in range(n_block - 1):
        H += J * (
            kron_at(X, b, n_block) @ kron_at(X, b + 1, n_block)
            + kron_at(Y, b, n_block) @ kron_at(Y, b + 1, n_block)
        )
    return H


def build_block_L(J: float, gamma: float, n_block: int) -> np.ndarray:
    """L_block = -i[H_B, ·] + Σ_l γ(Z_l ρ Z_l - ρ) on n_block-qubit block.

    Returns (d² × d²) super-operator on column-major vec(ρ):
    vec(M)[b·d + a] = M[a, b].
    """
    d = 2**n_block
    H = build_block_H(J, n_block)
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for l in range(n_block):
        Zl = kron_at(Z, l, n_block)
        L += gamma * (np.kron(Zl.T, Zl) - np.kron(Id, Id))
    return L


def block_bit_pos(n_block: int) -> list[int]:
    """State-index bit-position weights: |b_0 b_1 ... b_{n-1}⟩ → Σ b_i·2^(n-1-i)."""
    return [2 ** (n_block - 1 - i) for i in range(n_block)]


def state_idx(bits: list[int], bit_pos: list[int] | None = None) -> int:
    """Compute basis-state index from bit list. Pass bit_pos to avoid recomputing."""
    if bit_pos is None:
        bit_pos = block_bit_pos(len(bits))
    return sum(bit_pos[i] * bits[i] for i in range(len(bits)))


def compute_rho_block_0(n_block: int, N: int) -> np.ndarray:
    """ρ_block(0) = Tr_{N-n_block bare}(ρ_cc) for any path-k block at N qubits.

    ρ_cc = (|S_1⟩⟨S_2| + h.c.)/2; partial trace yields:
      Term 1 (popcount(c)=0): Σ_{i ∈ block} Σ_{j<k both in block} |SE_i^B⟩⟨DE_{jk}^B|
      Term 2 (popcount(c)=1): N_E·Σ_{j ∈ block} |0^B⟩⟨SE_j^B|
    All scaled by pre = 1/√(N·C(N, 2)), then h.c. + /2.
    """
    if N < n_block:
        raise ValueError(f"N={N} must be ≥ n_block={n_block}")
    d = 2**n_block
    bit_pos = block_bit_pos(n_block)
    N_E = N - n_block
    pre = 1.0 / np.sqrt(N * N * (N - 1) / 2)

    rho = np.zeros((d, d), dtype=complex)

    for i in range(n_block):
        bits = [0] * n_block
        bits[i] = 1
        idx_se = state_idx(bits, bit_pos)
        for j in range(n_block):
            for k in range(j + 1, n_block):
                bits_de = [0] * n_block
                bits_de[j] = 1
                bits_de[k] = 1
                rho[idx_se, state_idx(bits_de, bit_pos)] += pre

    for j in range(n_block):
        bits = [0] * n_block
        bits[j] = 1
        rho[0, state_idx(bits, bit_pos)] += pre * N_E

    rho = (rho + rho.conj().T) / 2.0
    return rho


def reduce_block_to_site_01(rho_block: np.ndarray, l: int, n_block: int) -> complex:
    """⟨0_l|Tr_{block\\{l}}(ρ_block)|1_l⟩ for block site l."""
    other = [s for s in range(n_block) if s != l]
    bit_pos = block_bit_pos(n_block)
    val = 0.0 + 0.0j
    for c in range(2 ** (n_block - 1)):
        bits_other = [(c >> (n_block - 2 - i)) & 1 for i in range(n_block - 1)]
        idx_0 = sum(bit_pos[other[i]] * bits_other[i] for i in range(n_block - 1))
        idx_1 = idx_0 + bit_pos[l]
        val += rho_block[idx_0, idx_1]
    return val


def bare_site_initial_01(N: int) -> float:
    """(ρ_l)_{0,1}(0) for any bare site l: (N-1) / (2·√(N·C(N,2)))."""
    return (N - 1) / (2.0 * np.sqrt(N * N * (N - 1) / 2))


def per_site_reduction_matrix(n_block: int) -> np.ndarray:
    """Build w[l, idx] for vectorised per-site reduction.

    w[l] is a column vector of length d² such that Σ_idx w[l, idx]·vec(rho)[idx]
    equals reduce_block_to_site_01(rho, l, n_block). Used in eigenmode
    projection workflows.
    """
    d = 2**n_block
    bit_pos = block_bit_pos(n_block)
    w = np.zeros((n_block, d * d), dtype=complex)
    for l in range(n_block):
        other = [s for s in range(n_block) if s != l]
        for c in range(2 ** (n_block - 1)):
            bits_other = [(c >> (n_block - 2 - i)) & 1 for i in range(n_block - 1)]
            idx_0 = sum(bit_pos[other[i]] * bits_other[i] for i in range(n_block - 1))
            idx_1 = idx_0 + bit_pos[l]
            w[l, idx_1 * d + idx_0] = 1.0
    return w
