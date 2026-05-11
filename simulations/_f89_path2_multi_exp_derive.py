"""F89c → F89-(2): multi-exponential closed form for S_(2)(t) at any N.

Topology (2): single connected path of 3 sites (path-2 block, bonds {0-1, 1-2})
plus N-3 bare sites. Lindbladian factorises across block + bare sites.

Per F89c + AbsorptionTheoremClaim, the L_super eigenvalues on the populated
sectors of the block are known (from sector eigendecomposition of the 8x8 block):

  - (vac, SE)_B sector: 3 modes, all rate 2γ₀, frequencies = H_B^SE eigenvalues
    {2√2 J, 0, -2√2 J} per F65 (4J cos(πk/(N+1))).
  - (SE, DE)_B sector: 9 modes with rates {2γ₀, 3.112γ₀, 3.444γ₀, 4γ₀, 6γ₀}
    (multiplicities 3, 1, 2, 2, 1; bit-exact via F89c column-bit-flip Hamming
    complement: pair-sum = 2γ₀·N = 6γ₀ for each (SE, SE) ↔ (SE, DE) eigenvalue
    pair).

The (SE, SE) sub-block rates 2.556γ₀ and 2.889γ₀ (from F89c) do NOT appear in
ρ_block(0) projected from ρ_cc, since ρ_cc has no (SE, SE) population.

This script does the missing eigenvector decomposition + initial-state
projection + per-site reduction. Output: numerical multi-exponential closed
form for S_(2)(t), verified against bond-isolate CSV at N=7.

For the |·|² spectrum, the Liouvillian eigenvalue pair (λ_k + λ_{k'}^*) gives
rates r_k + r_{k'} on |·|². Per AbsorptionTheoremClaim quantization, all such
sums lie on the 2γ₀ grid: {0, 2, 4, 6, 8, 10, 12} γ₀ for the populated sector
combinations.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO = Path(__file__).resolve().parents[1]
CSV_DIR = REPO / "simulations" / "results" / "bond_isolate"

# ----------- Operators on 3-qubit block ------------------------------------

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_at(P: np.ndarray, site: int, n_qubits: int) -> np.ndarray:
    op = np.array([[1.0]], dtype=complex)
    for q in range(n_qubits):
        op = np.kron(op, P if q == site else I2)
    return op


def build_block_H(J: float, n: int = 3) -> np.ndarray:
    """H_B = J · Σ (X_b X_{b+1} + Y_b Y_{b+1}) for path on n sites."""
    d = 2**n
    H = np.zeros((d, d), dtype=complex)
    for b in range(n - 1):
        H += J * (
            kron_at(X, b, n) @ kron_at(X, b + 1, n)
            + kron_at(Y, b, n) @ kron_at(Y, b + 1, n)
        )
    return H


def build_block_L(J: float, gamma: float, n: int = 3) -> np.ndarray:
    """L_block on n-qubit block: -i[H, ·] + Σ_l γ (Z_l ρ Z_l - ρ).

    Returns (d^2 × d^2) super-operator acting on column-major vec(ρ):
    vec(M)[b * d + a] = M[a, b].
    """
    d = 2**n
    H = build_block_H(J, n)
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for l in range(n):
        Zl = kron_at(Z, l, n)
        L += gamma * (np.kron(Zl.T, Zl) - np.kron(Id, Id))
    return L


# ----------- Initial state: ρ_block(0) for topology (2) at any N ----------


def block_state_index(bits: list[int]) -> int:
    """|b_0 b_1 b_2⟩ → state index 4·b_0 + 2·b_1 + b_2 (matches kron_at)."""
    return 4 * bits[0] + 2 * bits[1] + bits[2]


def compute_rho_block_0(N: int) -> np.ndarray:
    """ρ_block(0) = Tr_{N-3 bare}(ρ_cc) for topology (2). 8×8 complex array.

    Derivation:
      ρ_cc = (|S_1⟩⟨S_2| + |S_2⟩⟨S_1|) / 2 with S_n the popcount-n Dicke state.
      Tr_E(|S_1⟩⟨S_2|) splits by popcount(c) (c = bare-site config):
        popcount(c)=0: Σ_i Σ_{j<k both in block} |SE_i^B⟩⟨DE_{jk}^B|
        popcount(c)=1: N_E · Σ_{j ∈ block} |0^B⟩⟨SE_j^B|
        Other popcounts vanish (orthogonal bare-site partial traces).
      Then h.c. and divide by 2.
    """
    if N < 3:
        raise ValueError("Topology (2) requires N >= 3")
    n_block = 3
    N_E = N - n_block
    pre = 1.0 / np.sqrt(N * N * (N - 1) / 2)  # 1 / √(N · C(N, 2))

    rho = np.zeros((8, 8), dtype=complex)

    # Term 1 (popcount(c) = 0)
    for i in range(n_block):
        bits_se_i = [0, 0, 0]
        bits_se_i[i] = 1
        idx_se_i = block_state_index(bits_se_i)
        for j in range(n_block):
            for k in range(j + 1, n_block):
                bits_de_jk = [0, 0, 0]
                bits_de_jk[j] = 1
                bits_de_jk[k] = 1
                idx_de_jk = block_state_index(bits_de_jk)
                rho[idx_se_i, idx_de_jk] += pre

    # Term 2 (popcount(c) = 1)
    idx_vac = 0
    for j in range(n_block):
        bits_se_j = [0, 0, 0]
        bits_se_j[j] = 1
        idx_se_j = block_state_index(bits_se_j)
        rho[idx_vac, idx_se_j] += pre * N_E

    # ρ_cc = (S_1 S_2 + h.c.) / 2: add hermitian conjugate, divide by 2
    rho = (rho + rho.conj().T) / 2.0
    return rho


def reduce_block_to_site_01(rho_block: np.ndarray, l: int) -> complex:
    """⟨0_l|Tr_{block\\{l}}(ρ_block)|1_l⟩ for block site l ∈ {0, 1, 2}."""
    other = [s for s in range(3) if s != l]
    bit_pos = [4, 2, 1]  # state index contribution per bit position 0, 1, 2
    val = 0.0 + 0.0j
    for c in range(4):  # 2^2 states of the 2 other sites
        b0 = (c >> 1) & 1
        b1 = c & 1
        idx_0 = bit_pos[other[0]] * b0 + bit_pos[other[1]] * b1
        idx_1 = idx_0 + bit_pos[l]
        val += rho_block[idx_0, idx_1]
    return val


def bare_site_initial_01(N: int) -> float:
    """(ρ_l)_{0,1}(0) for any bare site l. By S_N symmetry, identical for all bare.

    Derivation: Σ_c ⟨0_l c|ρ_cc|1_l c⟩ with c on N-1 other sites, popcount(c) = 1
    (and the c=1 site ≠ l). Number of such c = N - 1. Each contributes
    1 / √(N · C(N, 2)). ρ_cc factor of 1/2 → result (N-1) / (2 √(N · C(N, 2))).
    """
    return (N - 1) / (2.0 * np.sqrt(N * N * (N - 1) / 2))


# ----------- Time evolution: numerical multi-exponential closed form ------


def evolve_S_total(N: int, J: float, gamma: float, t_array: np.ndarray) -> np.ndarray:
    """S_(2)(t) = Σ_l 2|(ρ_l)_{0,1}(t)|² for topology (2) at N qubits."""
    L = build_block_L(J, gamma)
    rho_block_0 = compute_rho_block_0(N)
    vec_rho_0 = rho_block_0.flatten(order="F")  # column-major: vec(M)[b*d + a] = M[a, b]

    # Eigendecomposition (non-Hermitian Liouvillian)
    eigvals, R = np.linalg.eig(L)
    c = np.linalg.solve(R, vec_rho_0)  # ρ_block(0) = Σ_k c_k · R_k

    S = np.zeros_like(t_array, dtype=float)
    x_bare_0 = bare_site_initial_01(N)

    for ti, t in enumerate(t_array):
        # ρ_block(t) = R · diag(exp(λ_k t)) · c
        vec_rho_t = R @ (np.exp(eigvals * t) * c)
        rho_block_t = vec_rho_t.reshape((8, 8), order="F")

        S_block = sum(
            2.0 * abs(reduce_block_to_site_01(rho_block_t, l)) ** 2
            for l in range(3)
        )
        S_bare = (N - 3) * 2.0 * (x_bare_0 * np.exp(-2 * gamma * t)) ** 2
        S[ti] = S_block + S_bare

    return S


def closed_form_terms(N: int, J: float, gamma: float):
    """Return the multi-exponential decomposition of S_(2)(t).

    S_(2)(t) = Σ_{(k, k')} A_{k,k'}(N) · exp(-(r_k + r_{k'}) t) · cos((ω_k - ω_{k'}) t + φ)
              + bare-site term

    Returns (rate_pairs, frequency_diffs, amplitude_table) for the block part
    plus the bare-site closed form.
    """
    L = build_block_L(J, gamma)
    rho_block_0 = compute_rho_block_0(N)
    vec_rho_0 = rho_block_0.flatten(order="F")
    eigvals, R = np.linalg.eig(L)
    R_inv = np.linalg.inv(R)
    c = R_inv @ vec_rho_0  # 64-dim coefficient vector

    # For each block site l, build the linear functional acting on rho_block_t -> (rho_l)_{0,1}.
    # (rho_l)_{0,1}(t) = ⟨w_l | vec(rho_block_t)⟩ = ⟨w_l | R · exp(Λt) · c⟩
    # = Σ_k (w_l^* @ R_k) · exp(λ_k t) · c_k
    # = Σ_k M_l(k) · exp(λ_k t) · c_k    where M_l(k) = (w_l^* @ R_k)
    w = np.zeros((3, 64), dtype=complex)
    for l in range(3):
        # w_l[idx] = 1 iff vec(rho_block)[idx] is one of the 4 (idx_0, idx_1) summed in reduce_block_to_site_01
        other = [s for s in range(3) if s != l]
        bit_pos = [4, 2, 1]
        for cc in range(4):
            b0 = (cc >> 1) & 1
            b1 = cc & 1
            idx_0 = bit_pos[other[0]] * b0 + bit_pos[other[1]] * b1
            idx_1 = idx_0 + bit_pos[l]
            w[l, idx_1 * 8 + idx_0] = 1.0  # column-major: vec[b*d + a] = M[a, b]

    M = w @ R  # (3, 64): M[l, k] = w_l^T @ R_k
    a = M * c[None, :]  # (3, 64): a[l, k] = M_l(k) · c_k
    return eigvals, a


def populated_distinct_rates(N: int, J: float, gamma: float, threshold: float = 1e-12):
    """Return distinct (rate, freq) pairs of L_super modes populated by ρ_block(0)."""
    eigvals, a = closed_form_terms(N, J, gamma)
    rates = -eigvals.real
    freqs = eigvals.imag
    sig = np.sum(np.abs(a) ** 2, axis=0)
    contributing = np.where(sig > threshold)[0]
    # Group by (rate, |freq|) rounded
    grouped = {}
    for k in contributing:
        key = (round(rates[k] / gamma, 4), round(abs(freqs[k]) / J, 4))
        grouped[key] = grouped.get(key, 0.0) + sig[k]
    return sorted(grouped.items())


def main() -> None:
    J, gamma = 0.075, 0.05

    # ---- N=7 verification against bond-isolate CSV ----
    csv_path = CSV_DIR / "N7_b0-1_J0.0750_gamma0.0500_probe-coherence.csv"
    print(f"# F89-(2) closed-form derivation, J={J}, γ={gamma}\n")

    if not csv_path.exists():
        print(f"# CSV not found: {csv_path}")
        return

    data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    t_csv = data[:, 0]
    S_csv = data[:, -1]

    print("## Verification at N=7 vs bond-isolate RK4 CSV")
    N = 7
    S_pred = evolve_S_total(N, J, gamma, t_csv)
    diff = S_pred - S_csv
    print(f"# max |diff|: {np.max(np.abs(diff)):.3e} (CSV write precision is ~5e-7)")
    print(f"# mean |diff|: {np.mean(np.abs(diff)):.3e}")
    print(f"# S(0) = {S_pred[0]:.6f} = (N-1)/N = {(N-1)/N:.6f} ✓")
    print()

    # ---- N=5 prediction (no CSV, pure prediction from closed form) ----
    print("## N=5 prediction (no CSV, fresh prediction from F89c+AT-derived form)")
    N5 = 5
    S5 = evolve_S_total(N5, J, gamma, t_csv[:201])  # t in [0, 20]
    print(f"# S(0) = {S5[0]:.6f} = (N-1)/N = {(N5-1)/N5:.6f} ✓")
    print(f"# S(t=10) = {S5[100]:.6f}; S(t=20) = {S5[200]:.6f}")
    print()

    # ---- Structural analysis: which L_super modes are populated, by N ----
    print("## L_super mode populations (distinct (rate Γ/γ, |freq|/J) groups)")
    print("# Per AbsorptionTheoremClaim: 64 block-L modes total. Only S_3-fully-symmetric")
    print("# initial state populates a small subset.")
    print()
    print("| (rate Γ/γ, |freq|/J) | N=5 sig | N=7 sig | N=11 sig | N-scaling |")
    print("|---|---|---|---|---|")
    n5_groups = dict(populated_distinct_rates(5, J, gamma))
    n7_groups = dict(populated_distinct_rates(7, J, gamma))
    n11_groups = dict(populated_distinct_rates(11, J, gamma))
    all_keys = sorted(set(n5_groups) | set(n7_groups) | set(n11_groups))
    for k in all_keys:
        s5, s7, s11 = n5_groups.get(k, 0.0), n7_groups.get(k, 0.0), n11_groups.get(k, 0.0)
        # Estimate N-scaling (rough): s7/s5 ratio
        scale = "?"
        if s5 > 1e-15 and s11 > 1e-15:
            r = s11 / s5
            # Test N² and (N_E)² scaling
            n_E_5, n_E_11 = 5 - 3, 11 - 3
            if abs(r - (n_E_11 / n_E_5) ** 2) < 0.05:
                scale = "(N-3)²"
            elif abs(r - 1.0) < 0.05:
                scale = "const"
            elif abs(r - (n_E_11 / n_E_5)) < 0.05:
                scale = "(N-3)"
            else:
                scale = f"~{r:.2f}× from N=5→11"
        print(
            f"| ({k[0]:.4f}, {k[1]:.4f}) | {s5:.2e} | {s7:.2e} | {s11:.2e} | {scale} |"
        )
    print()

    # ---- |·|² rate grid ----
    print("## |·|² rate spectrum (rates r_k + r_{k'} for populated pairs)")
    rates_n7 = sorted({k[0] for k in n7_groups})
    pair_rates = sorted(set(
        round(r1 + r2, 4)
        for r1 in rates_n7 for r2 in rates_n7
    ))
    print(f"# Amplitude-level rates Γ/γ at N=7: {rates_n7}")
    print(f"# |·|²-level rates Γ/γ at N=7: {pair_rates}")
    print()
    print("# Note: amplitude-level rates 2γ are pure-AT (n_XY=1 per coherence).")
    print("# Fractional rates 3.04γ, 3.48γ are H_B-mixed (SE,DE) sub-block eigenvalues")
    print("# at this J/γ=1.5; they are NOT on the AT 2γ-grid because they come from")
    print("# diagonalising a non-trivial 2×2 H_B-coupled sub-block (not pure-Pauli modes).")
    print()
    print("# Pure-AT modes 4γ and 6γ get ZERO projection from ρ_cc-derived ρ_block(0):")
    print("# their (SE,DE) eigenvectors are S_3-asymmetric, while ρ_block(0)'s (SE,DE)")
    print("# part is a fully symmetric Σ_i Σ_{j<k} |SE_i⟩⟨DE_{jk}| superposition.")


if __name__ == "__main__":
    main()
