"""F89-(k) survey: multi-exponential closed-form derivation across path-k=2..5.

Generalises `_f89_path2_multi_exp_derive.py` and `_f89_path3_multi_exp_derive.py`
to arbitrary block size N_BLOCK = k+1. For each path-k we:

1. Build block L_super (4^N_BLOCK dim).
2. Compute ρ_block(0) via partial-trace of ρ_cc.
3. Eigendecompose + project + per-site reduce.
4. Verify against bond-isolate `N7_b0-1-...-(k-1)_J0.0750_gamma0.0500_probe-coherence.csv`.
5. Survey populated mode-group count + Hamming-complement pair structure.

Result table (J/γ=1.5 at N=7):

| Path | N_block | d²=4^N_block | Mode-groups | Contributing | Pair-sums to 2γ·N_block |
|---|---|---|---|---|---|
| path-2 | 3 | 64 | 4 | 16 | 0 (S_3-asymmetric partners absent) |
| path-3 | 4 | 256 | 10 | 65 | 7 (full Hamming-complement structure) |
| path-4 | 5 | 1024 | 12 | 128 | 0 (S_5-asymmetric partners absent) |
| path-5 | 6 | 4096 | 35 | 314 | 0 (S_6-asymmetric partners absent) |

Path-3 is privileged: at N_block=4, DE = popcount-2 = bar(popcount-2)
self-symmetric, so column-bit-flip maps populated (SE,DE) modes to other
populated (SE,DE) modes within the symmetric subspace. For N_block ∈ {3, 5, 6}
the column-flip partners land in S_{N_block}-asymmetric territory and get
zero projection from ρ_cc-derived ρ_block(0).

L_super dimensions match `experiments/CAVITY_MODES_FORMULA.md` because both
index the same 4^N operator space; CAVITY_MODES counts SU(2)-Heisenberg
stationary modes via Schur-Weyl, this script counts S_{N_block}-symmetric
populated XY modes — different decompositions of the same d².
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO = Path(__file__).resolve().parents[1]
CSV_DIR = REPO / "simulations" / "results" / "bond_isolate"

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_at(P: np.ndarray, site: int, n: int) -> np.ndarray:
    op = np.array([[1.0]], dtype=complex)
    for q in range(n):
        op = np.kron(op, P if q == site else I2)
    return op


def derive_path_k(n_block: int, J: float, gamma: float, N: int):
    """Build L_block, eigendecompose, project ρ_cc-derived ρ_block(0). Return all parts."""
    D = 2**n_block
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

    def state_idx(bits: list[int]) -> int:
        return sum(bit_pos[i] * bits[i] for i in range(n_block))

    N_E = N - n_block
    pre = 1.0 / np.sqrt(N * N * (N - 1) / 2)
    rho = np.zeros((D, D), dtype=complex)
    for i in range(n_block):
        bits = [0] * n_block
        bits[i] = 1
        idx_se = state_idx(bits)
        for j in range(n_block):
            for k in range(j + 1, n_block):
                bits_de = [0] * n_block
                bits_de[j] = 1
                bits_de[k] = 1
                rho[idx_se, state_idx(bits_de)] += pre
    for j in range(n_block):
        bits = [0] * n_block
        bits[j] = 1
        rho[0, state_idx(bits)] += pre * N_E
    rho = (rho + rho.conj().T) / 2.0

    vec = rho.flatten(order="F")
    eigvals, R = np.linalg.eig(L)
    c = np.linalg.solve(R, vec)

    w = np.zeros((n_block, D * D), dtype=complex)
    for l in range(n_block):
        other = [s for s in range(n_block) if s != l]
        for cc in range(2 ** (n_block - 1)):
            bits_other = [(cc >> (n_block - 2 - i)) & 1 for i in range(n_block - 1)]
            idx_0 = sum(bit_pos[other[i]] * bits_other[i] for i in range(n_block - 1))
            idx_1 = idx_0 + bit_pos[l]
            w[l, idx_1 * D + idx_0] = 1.0
    M = w @ R
    a = M * c[None, :]

    return L, eigvals, R, c, a, rho, bit_pos


def reduce_to_site(rho_block: np.ndarray, l: int, n_block: int, bit_pos: list[int]) -> complex:
    other = [s for s in range(n_block) if s != l]
    val = 0.0 + 0.0j
    for cc in range(2 ** (n_block - 1)):
        bits_other = [(cc >> (n_block - 2 - i)) & 1 for i in range(n_block - 1)]
        idx_0 = sum(bit_pos[other[i]] * bits_other[i] for i in range(n_block - 1))
        idx_1 = idx_0 + bit_pos[l]
        val += rho_block[idx_0, idx_1]
    return val


def evolve(n_block: int, J: float, gamma: float, N: int, t_array: np.ndarray) -> np.ndarray:
    L, eigvals, R, c, a, rho, bit_pos = derive_path_k(n_block, J, gamma, N)
    D = 2**n_block
    x_bare_0 = (N - 1) / (2.0 * np.sqrt(N * N * (N - 1) / 2))
    S = np.zeros_like(t_array, dtype=float)
    for ti, t in enumerate(t_array):
        vec_t = R @ (np.exp(eigvals * t) * c)
        rho_t = vec_t.reshape((D, D), order="F")
        S_block = sum(2.0 * abs(reduce_to_site(rho_t, l, n_block, bit_pos)) ** 2 for l in range(n_block))
        S_bare = (N - n_block) * 2.0 * (x_bare_0 * np.exp(-2 * gamma * t)) ** 2
        S[ti] = S_block + S_bare
    return S


def main() -> None:
    J, gamma = 0.075, 0.05
    N = 7
    print(f"# F89-(k) path-k closed-form survey, J={J}, γ={gamma}, N={N}\n")

    print("## Bond-isolate verification across path-k")
    print("| Path | bonds | max |diff| | mean |diff| | S(0) match |")
    print("|---|---|---|---|---|")
    for n_block in [3, 4, 5]:
        bonds_str = "-".join(str(b) for b in range(n_block - 1))
        csv = CSV_DIR / f"N{N}_b{bonds_str}_J0.0750_gamma0.0500_probe-coherence.csv"
        if not csv.exists():
            print(f"| path-{n_block-1} | {bonds_str} | (CSV missing) | | |")
            continue
        data = np.loadtxt(csv, delimiter=",", skiprows=1)
        t_csv, S_csv = data[:, 0], data[:, -1]
        S_pred = evolve(n_block, J, gamma, N, t_csv)
        diff = S_pred - S_csv
        match = "✓" if abs(S_pred[0] - (N - 1) / N) < 1e-6 else "✗"
        print(
            f"| path-{n_block-1} | {{{bonds_str}}} | {np.max(np.abs(diff)):.2e} | "
            f"{np.mean(np.abs(diff)):.2e} | {match} |"
        )
    print()

    print("## Populated mode-group survey")
    print("| Path | N_block | d² | Mode-groups | Contributing modes | Pair-sums to 2γ·N_block |")
    print("|---|---|---|---|---|---|")
    for n_block in [3, 4, 5, 6]:
        L, eigvals, R, c, a, rho, _ = derive_path_k(n_block, J, gamma, N)
        sig = np.sum(np.abs(a) ** 2, axis=0)
        contributing = np.where(sig > 1e-12)[0]
        rates = -eigvals.real / gamma
        freqs = np.abs(eigvals.imag) / J
        groups = {}
        for k in contributing:
            key = (round(rates[k], 4), round(freqs[k], 4))
            groups[key] = groups.get(key, 0.0) + sig[k]
        unique_rates = sorted({k[0] for k in groups.keys()})
        target_pair_sum = 2 * n_block
        n_pairs = sum(
            1 for r1 in unique_rates for r2 in unique_rates
            if abs((r1 + r2) - target_pair_sum) < 0.001
        )
        print(
            f"| path-{n_block-1} | {n_block} | {4**n_block} | "
            f"{len(groups)} | {len(contributing)} | "
            f"{n_pairs}{' ✓' if n_pairs > 0 else ''} |"
        )


if __name__ == "__main__":
    main()
