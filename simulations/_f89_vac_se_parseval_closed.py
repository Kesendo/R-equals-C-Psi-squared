"""F89 (vac, SE) self-contribution closed form via Parseval orthogonality.

For any pure-path-k block of N_block = k+1 sites embedded in N qubits, the
(vac, SE)_B sector self-contribution to S(t) (i.e. the part of the per-site
spatial-sum coherence that comes ONLY from the (vac, SE) block of L_super,
ignoring cross-products with the (SE, DE) sector) has the closed form:

    S^{(vac,SE)}_block(t; k, N) = (k+1) · (N-k-1)² / (N² · (N-1)) · exp(-4γ₀ t)

Pure exponential decay at rate 4γ₀, no oscillation — because Parseval
orthogonality of H_B^SE eigenstates ψ_k eliminates all k ≠ k' cross-terms
when summed over the (k+1) block sites:

    Σ_l ψ_k(l)·ψ_{k'}(l) = δ_{k,k'}

This is the "smooth backbone" of S_(k)(t). The remaining structure (oscillations
at frequencies E_k − E_{k'} from the (SE, DE) sector and cross-products between
(vac, SE) and (SE, DE)) lives in the H_B-mixed (SE, DE) sub-block where no
clean orthogonality argument applies.

Derivation (in 4 lines):
  1. ρ_block(0)|_{(vac,SE)} = N_E·pre·√(N_block)/2 · |0⟩⟨α_{N_block}|
     where |α_{N_block}⟩ = (1/√N_block) Σ_j |SE_j⟩
  2. Decompose |α⟩ in H_B^SE eigenbasis: |α⟩ = Σ_k ⟨ψ_k|α⟩ |ψ_k⟩
  3. Time evolve each |0⟩⟨ψ_k| → exp(+iE_k t − 2γt) |0⟩⟨ψ_k|
  4. Σ_l 2|(ρ_l(t))^{(vac,SE)}_{0,1}|² = 2·Σ_k|coef_k|²·exp(-4γt)
                                        = N_block·(N_E·pre)²/2·exp(-4γt) (Parseval)
                                        = (k+1)·(N-k-1)²/(N²(N-1))·exp(-4γt)

Verification: compare against numerical (vac, SE)-only contribution from the
path-k script (with (SE, DE) sector zeroed out in ρ_block(0)).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO = Path(__file__).resolve().parents[1]

# ----------- Closed form ----------------------------------------------------


def S_vac_se_closed_form(k: int, N: int, gamma: float, t_array: np.ndarray) -> np.ndarray:
    """Closed-form (vac, SE) block self-contribution.

    S^{(vac,SE)}_block(t; k, N) = (k+1) · (N-k-1)² / (N² · (N-1)) · exp(-4γ₀ t)
    """
    n_block = k + 1
    n_E = N - n_block
    prefactor = n_block * n_E * n_E / (N * N * (N - 1))
    return prefactor * np.exp(-4.0 * gamma * t_array)


# ----------- Numerical verification: zero out (SE, DE) in ρ_block(0) -------


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_at(P: np.ndarray, site: int, n: int) -> np.ndarray:
    op = np.array([[1.0]], dtype=complex)
    for q in range(n):
        op = np.kron(op, P if q == site else I2)
    return op


def numerical_vac_se_block(k: int, N: int, J: float, gamma: float, t_array: np.ndarray) -> np.ndarray:
    """Compute (vac, SE)-only block contribution numerically by zeroing (SE, DE) in ρ_block(0)."""
    n_block = k + 1
    D = 2**n_block

    # Build block L_super
    H = np.zeros((D, D), dtype=complex)
    for b in range(n_block - 1):
        H += J * (
            kron_at(X, b, n_block) @ kron_at(X, b + 1, n_block)
            + kron_at(Y, b, n_block) @ kron_at(Y, b + 1, n_block)
        )
    Id = np.eye(D, dtype=complex)
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for l in range(n_block):
        Zl = kron_at(Z, l, n_block)
        L += gamma * (np.kron(Zl.T, Zl) - np.kron(Id, Id))

    bit_pos = [2 ** (n_block - 1 - i) for i in range(n_block)]

    def state_idx(bits):
        return sum(bit_pos[i] * bits[i] for i in range(n_block))

    N_E = N - n_block
    pre = 1.0 / np.sqrt(N * N * (N - 1) / 2)

    # Build ρ_block(0) with ONLY (vac, SE) sector populated (zero out (SE, DE) term)
    rho = np.zeros((D, D), dtype=complex)
    for j in range(n_block):
        bits = [0] * n_block
        bits[j] = 1
        rho[0, state_idx(bits)] += pre * N_E
    rho = (rho + rho.conj().T) / 2.0

    vec = rho.flatten(order="F")
    eigvals, R = np.linalg.eig(L)
    c = np.linalg.solve(R, vec)

    def reduce_to_site(rho_block, l):
        other = [s for s in range(n_block) if s != l]
        val = 0.0 + 0.0j
        for cc in range(2 ** (n_block - 1)):
            bits_other = [(cc >> (n_block - 2 - i)) & 1 for i in range(n_block - 1)]
            idx_0 = sum(bit_pos[other[i]] * bits_other[i] for i in range(n_block - 1))
            idx_1 = idx_0 + bit_pos[l]
            val += rho_block[idx_0, idx_1]
        return val

    S = np.zeros_like(t_array, dtype=float)
    for ti, t in enumerate(t_array):
        vec_t = R @ (np.exp(eigvals * t) * c)
        rho_t = vec_t.reshape((D, D), order="F")
        S[ti] = sum(2.0 * abs(reduce_to_site(rho_t, l)) ** 2 for l in range(n_block))

    return S


def main() -> None:
    J, gamma = 0.075, 0.05
    t_array = np.linspace(0, 30, 301)

    print(f"# F89 (vac, SE) closed form via Parseval, J={J}, γ={gamma}\n")
    print("# Closed form: S^(vac,SE)_block(t; k, N) = (k+1)(N-k-1)²/(N²(N-1)) · exp(-4γ₀ t)\n")

    print("## Verification: closed form vs numerical (with (SE, DE) zeroed out)")
    print("| k | N | (k+1)(N-k-1)²/(N²(N-1)) | max |closed - numerical| |")
    print("|---|---|---|---|")
    for k in [1, 2, 3, 4, 5]:
        for N in [k + 2, k + 4, k + 6]:  # spans different N_E
            if k + 1 > N:
                continue
            S_closed = S_vac_se_closed_form(k, N, gamma, t_array)
            S_numerical = numerical_vac_se_block(k, N, J, gamma, t_array)
            max_diff = np.max(np.abs(S_closed - S_numerical))
            prefactor = (k + 1) * (N - k - 1) ** 2 / (N * N * (N - 1))
            print(f"| {k} | {N} | {prefactor:.6f} | {max_diff:.3e} |")

    print()
    print("## All differences should be at machine precision (~1e-14).")
    print("# Parseval orthogonality is exact, so the closed form is bit-exact.")


if __name__ == "__main__":
    main()
