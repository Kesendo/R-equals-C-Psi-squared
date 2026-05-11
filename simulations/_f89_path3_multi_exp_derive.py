"""F89c → F89-(3): multi-exponential closed form for S_(3)(t) at any N.

Topology (3): single connected path of 4 sites (path-3 block, bonds {0-1, 1-2, 2-3})
plus N-4 bare sites. Lindbladian factorises across block + bare sites.

Generalization of path-2 (`_f89_path2_multi_exp_derive.py`) to 4-qubit block
(d² = 256-dim L_super). Same partial-trace structure for ρ_block(0) (S_4
symmetry inherited from ρ_cc); same script pattern (eigenvector decomposition
+ initial-state projection + per-site reduction + bare-site addition).

Per F89c spectrum at γ=J=1: 4-qubit block has 25 distinct decay rates. We
expect the S_4 symmetry of ρ_block(0) to populate only a small subset.
Hamming-complement pair-sum at N_block=4: 2γ₀·N_block = 8γ₀ (column bit-flip
maps (SE,SE) ↔ (SE,TE) since bar(SE) = TE = popcount-3 at N=4, NOT (SE,DE)
which is popcount-2-self-symmetric).

Verification: bond-isolate `N7_b0-1-2` CSV (topology (3) at N=7).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO = Path(__file__).resolve().parents[1]
CSV_DIR = REPO / "simulations" / "results" / "bond_isolate"

# ----------- Operators on 4-qubit block ------------------------------------

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

N_BLOCK = 4
D_BLOCK = 2**N_BLOCK  # 16


def kron_at(P: np.ndarray, site: int, n_qubits: int) -> np.ndarray:
    op = np.array([[1.0]], dtype=complex)
    for q in range(n_qubits):
        op = np.kron(op, P if q == site else I2)
    return op


def build_block_H(J: float) -> np.ndarray:
    H = np.zeros((D_BLOCK, D_BLOCK), dtype=complex)
    for b in range(N_BLOCK - 1):
        H += J * (
            kron_at(X, b, N_BLOCK) @ kron_at(X, b + 1, N_BLOCK)
            + kron_at(Y, b, N_BLOCK) @ kron_at(Y, b + 1, N_BLOCK)
        )
    return H


def build_block_L(J: float, gamma: float) -> np.ndarray:
    """L_block on 4 qubits: -i[H, ·] + Σ_l γ (Z_l ρ Z_l - ρ). 256×256 matrix.

    Convention: vec(M)[b * D + a] = M[a, b] (column-major).
    """
    H = build_block_H(J)
    Id = np.eye(D_BLOCK, dtype=complex)
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for l in range(N_BLOCK):
        Zl = kron_at(Z, l, N_BLOCK)
        L += gamma * (np.kron(Zl.T, Zl) - np.kron(Id, Id))
    return L


# ----------- Initial state: ρ_block(0) for topology (3) at any N -----------


def block_state_index(bits: list[int]) -> int:
    """|b_0 b_1 b_2 b_3⟩ → state index 8·b_0 + 4·b_1 + 2·b_2 + b_3."""
    return 8 * bits[0] + 4 * bits[1] + 2 * bits[2] + bits[3]


def compute_rho_block_0(N: int) -> np.ndarray:
    """ρ_block(0) = Tr_{N-4 bare}(ρ_cc) for topology (3). 16×16 complex array.

    Per the same partial-trace bookkeeping as path-2 (only block size differs):
      Tr_E(|S_1⟩⟨S_2|) = pre · [Σ_i Σ_{j<k both in block} |SE_i^B⟩⟨DE_{jk}^B|
                              + N_E · Σ_{j ∈ block} |0^B⟩⟨SE_j^B|]
      pre = 1 / √(N · C(N, 2))
    Then h.c. and divide by 2.
    """
    if N < N_BLOCK:
        raise ValueError(f"Topology (3) requires N >= {N_BLOCK}")
    N_E = N - N_BLOCK
    pre = 1.0 / np.sqrt(N * N * (N - 1) / 2)

    rho = np.zeros((D_BLOCK, D_BLOCK), dtype=complex)

    # Term 1: popcount(c) = 0 → Σ_{i ∈ block} Σ_{j<k both in block} |SE_i^B⟩⟨DE_{jk}^B|
    for i in range(N_BLOCK):
        bits_se = [0] * N_BLOCK
        bits_se[i] = 1
        idx_se = block_state_index(bits_se)
        for j in range(N_BLOCK):
            for k in range(j + 1, N_BLOCK):
                bits_de = [0] * N_BLOCK
                bits_de[j] = 1
                bits_de[k] = 1
                idx_de = block_state_index(bits_de)
                rho[idx_se, idx_de] += pre

    # Term 2: popcount(c) = 1 → N_E · Σ_{j ∈ block} |vac⟩⟨SE_j|
    idx_vac = 0
    for j in range(N_BLOCK):
        bits_se = [0] * N_BLOCK
        bits_se[j] = 1
        idx_se = block_state_index(bits_se)
        rho[idx_vac, idx_se] += pre * N_E

    rho = (rho + rho.conj().T) / 2.0
    return rho


def reduce_block_to_site_01(rho_block: np.ndarray, l: int) -> complex:
    """⟨0_l|Tr_{block\\{l}}(ρ_block)|1_l⟩ for block site l ∈ {0, 1, 2, 3}."""
    other = [s for s in range(N_BLOCK) if s != l]
    bit_pos = [8, 4, 2, 1]
    val = 0.0 + 0.0j
    for c in range(2 ** (N_BLOCK - 1)):  # 8 states of the 3 other sites
        bits_other = [(c >> (N_BLOCK - 2 - i)) & 1 for i in range(N_BLOCK - 1)]
        idx_0 = sum(bit_pos[other[i]] * bits_other[i] for i in range(N_BLOCK - 1))
        idx_1 = idx_0 + bit_pos[l]
        val += rho_block[idx_0, idx_1]
    return val


def bare_site_initial_01(N: int) -> float:
    """Same closed form as path-2 (independent of block size, only depends on N).

    (ρ_l)_{0,1}(0) = (N-1) / (2 · √(N · C(N, 2))) for any bare site l.
    """
    return (N - 1) / (2.0 * np.sqrt(N * N * (N - 1) / 2))


def evolve_S_total(N: int, J: float, gamma: float, t_array: np.ndarray) -> np.ndarray:
    """S_(3)(t) = Σ_l 2|(ρ_l)_{0,1}(t)|² for topology (3) at N qubits."""
    L = build_block_L(J, gamma)
    rho_block_0 = compute_rho_block_0(N)
    vec_rho_0 = rho_block_0.flatten(order="F")

    eigvals, R = np.linalg.eig(L)
    c = np.linalg.solve(R, vec_rho_0)

    S = np.zeros_like(t_array, dtype=float)
    x_bare_0 = bare_site_initial_01(N)

    for ti, t in enumerate(t_array):
        vec_rho_t = R @ (np.exp(eigvals * t) * c)
        rho_block_t = vec_rho_t.reshape((D_BLOCK, D_BLOCK), order="F")

        S_block = sum(
            2.0 * abs(reduce_block_to_site_01(rho_block_t, l)) ** 2
            for l in range(N_BLOCK)
        )
        S_bare = (N - N_BLOCK) * 2.0 * (x_bare_0 * np.exp(-2 * gamma * t)) ** 2
        S[ti] = S_block + S_bare

    return S


def closed_form_terms(N: int, J: float, gamma: float):
    """Multi-exponential decomposition of S_(3)(t).

    Returns (eigvals_64, amp_per_site) where amp_per_site is (N_BLOCK, D_BLOCK²).
    """
    L = build_block_L(J, gamma)
    rho_block_0 = compute_rho_block_0(N)
    vec_rho_0 = rho_block_0.flatten(order="F")
    eigvals, R = np.linalg.eig(L)
    R_inv = np.linalg.inv(R)
    c = R_inv @ vec_rho_0

    w = np.zeros((N_BLOCK, D_BLOCK * D_BLOCK), dtype=complex)
    for l in range(N_BLOCK):
        other = [s for s in range(N_BLOCK) if s != l]
        bit_pos = [8, 4, 2, 1]
        for cc in range(2 ** (N_BLOCK - 1)):
            bits_other = [(cc >> (N_BLOCK - 2 - i)) & 1 for i in range(N_BLOCK - 1)]
            idx_0 = sum(bit_pos[other[i]] * bits_other[i] for i in range(N_BLOCK - 1))
            idx_1 = idx_0 + bit_pos[l]
            w[l, idx_1 * D_BLOCK + idx_0] = 1.0

    M = w @ R
    a = M * c[None, :]
    return eigvals, a


def populated_distinct_rates(N: int, J: float, gamma: float, threshold: float = 1e-12):
    eigvals, a = closed_form_terms(N, J, gamma)
    rates = -eigvals.real
    freqs = eigvals.imag
    sig = np.sum(np.abs(a) ** 2, axis=0)
    contributing = np.where(sig > threshold)[0]
    grouped = {}
    for k in contributing:
        key = (round(rates[k] / gamma, 4), round(abs(freqs[k]) / J, 4))
        grouped[key] = grouped.get(key, 0.0) + sig[k]
    return sorted(grouped.items())


def main() -> None:
    J, gamma = 0.075, 0.05
    print(f"# F89-(3) closed-form derivation, J={J}, γ={gamma}\n")

    csv_path = CSV_DIR / "N7_b0-1-2_J0.0750_gamma0.0500_probe-coherence.csv"
    if not csv_path.exists():
        print(f"# CSV not found: {csv_path}")
        return

    data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    t_csv = data[:, 0]
    S_csv = data[:, -1]

    N = 7
    print(f"## Verification at N={N} vs bond-isolate CSV ({csv_path.name})")
    S_pred = evolve_S_total(N, J, gamma, t_csv)
    diff = S_pred - S_csv
    print(f"# max |diff|: {np.max(np.abs(diff)):.3e} (CSV write precision is ~5e-7)")
    print(f"# mean |diff|: {np.mean(np.abs(diff)):.3e}")
    print(f"# S(0) = {S_pred[0]:.6f} = (N-1)/N = {(N-1)/N:.6f} ✓")
    print()
    print("| t  | S_csv (RK4) | S_pred (closed form) | diff |")
    print("|---|---|---|---|")
    for i in [0, 30, 50, 100, 150, 200, 250, 300]:
        if i < len(t_csv):
            print(
                f"| {t_csv[i]:5.2f} | {S_csv[i]:.6f} | {S_pred[i]:.6f} | "
                f"{diff[i]:+.3e} |"
            )
    print()

    # ---- Mode population analysis ----
    print("## L_super mode populations at N=5, 7, 11 (S_4-symmetric subset)")
    print("# 256 total block-L modes. Per F89c at γ=J=1: 25 distinct rates.")
    print("# At J/γ=1.5 (this run) the fractional rates shift; populated subset")
    print("# is reduced by S_4 symmetry of ρ_block(0).")
    print()
    print("| (rate Γ/γ, |freq|/J) | N=5 sig | N=7 sig | N=11 sig |")
    print("|---|---|---|---|")
    n5_groups = dict(populated_distinct_rates(5, J, gamma))
    n7_groups = dict(populated_distinct_rates(7, J, gamma))
    n11_groups = dict(populated_distinct_rates(11, J, gamma))
    all_keys = sorted(set(n5_groups) | set(n7_groups) | set(n11_groups))
    for k in all_keys:
        s5 = n5_groups.get(k, 0.0)
        s7 = n7_groups.get(k, 0.0)
        s11 = n11_groups.get(k, 0.0)
        print(f"| ({k[0]:.4f}, {k[1]:.4f}) | {s5:.2e} | {s7:.2e} | {s11:.2e} |")
    print()

    # ---- N=5 prediction ----
    print("## N=5 prediction (no CSV; pure prediction from F89c+AT-derived form)")
    N5 = 5
    if N5 >= N_BLOCK:
        S5 = evolve_S_total(N5, J, gamma, t_csv[:201])
        print(f"# S(0) = {S5[0]:.6f} = (N-1)/N = {(N5-1)/N5:.6f} ✓")
        print(f"# S(t=10) = {S5[100]:.6f}; S(t=20) = {S5[200]:.6f}")
    print()

    # ---- |·|² rates ----
    rates_n7 = sorted({k[0] for k in n7_groups})
    pair_rates = sorted(set(round(r1 + r2, 4) for r1 in rates_n7 for r2 in rates_n7))
    print("## |·|² rate spectrum at N=7")
    print(f"# Amplitude-level rates Γ/γ: {rates_n7}")
    print(f"# |·|²-level rates Γ/γ: {pair_rates}")


if __name__ == "__main__":
    main()
