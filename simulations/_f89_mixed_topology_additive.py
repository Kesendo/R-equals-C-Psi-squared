"""F89 mixed topology: additive identity from Lindbladian factorisation.

For any topology T = (k_1, k_2, ..., k_m) at N qubits, the spatial-sum
coherence decomposes as:

    S_T(t) = Σ_i S_(k_i)(t)  −  (m − 1) · N · S_bare(t; N)

where:
- S_(k_i)(t) = closed-form per-pure-path-k_i contribution (k=1 analytic via
  all-isolated formula; k=2..6 numerical via path-k script)
- S_bare(t; N) = (N − 1) / N² · exp(−4γ₀ t) per bare site (closed form)
- m = number of disjoint blocks
- Subtraction term cancels overcounting of bare-site contributions: each
  per-block S_(k_i) includes (N − k_i − 1) bare sites worth of contribution,
  but the actual mixed topology has only (N − Σ_i (k_i + 1)) bare sites.
  Total double-count: (m − 1) · N.

Derivation: Lindbladian L = Σ_blocks L_block + Σ_bare L_l factorises across
disjoint blocks. ρ(t) = ⊗ exp(L_block · t)[Tr_else(ρ_cc)]. Per-site reduction
ρ_l(t) = Tr_else(ρ(t)) = exp(L_block_l · t)[ρ_block_l(0)] depends only on
the block containing l. Hence S_T(t) = Σ_l 2|(ρ_l)_{0,1}|² is sum of
per-block contributions.

The per-block ρ_block(0) = Tr_E(ρ_cc) depends on N (via the 1/√(N·C(N,2))
prefactor and N_E factor) but NOT on which OTHER blocks are present. So
the per-block S_(k_i) function is the same regardless of mixed-topology
context.

Verification: 14 topology classes at N=7 against bond-isolate CSVs.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO = Path(__file__).resolve().parents[1]
CSV_DIR = REPO / "simulations" / "results" / "bond_isolate"

# ----------- Pauli operators ------------------------------------------------

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_at(P: np.ndarray, site: int, n: int) -> np.ndarray:
    op = np.array([[1.0]], dtype=complex)
    for q in range(n):
        op = np.kron(op, P if q == site else I2)
    return op


# ----------- Per-block path-k closed form (numerical, k >= 2) --------------


def path_k_block_S(n_block: int, N: int, J: float, gamma: float, t_array: np.ndarray) -> np.ndarray:
    """Per-block contribution to S(t) from a (k+1)-site path block embedded in N qubits.

    Returns the contribution from the n_block block-sites only (excludes bare).
    For (k+1)-site path block, this is the k-th member of the path-k family.
    """
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

    def reduce_to_site(rho_block, l):
        other = [s for s in range(n_block) if s != l]
        val = 0.0 + 0.0j
        for cc in range(2 ** (n_block - 1)):
            bits_other = [(cc >> (n_block - 2 - i)) & 1 for i in range(n_block - 1)]
            idx_0 = sum(bit_pos[other[i]] * bits_other[i] for i in range(n_block - 1))
            idx_1 = idx_0 + bit_pos[l]
            val += rho_block[idx_0, idx_1]
        return val

    S_block_contrib = np.zeros_like(t_array, dtype=float)
    for ti, t in enumerate(t_array):
        vec_t = R @ (np.exp(eigvals * t) * c)
        rho_t = vec_t.reshape((D, D), order="F")
        S_block_contrib[ti] = sum(2.0 * abs(reduce_to_site(rho_t, l)) ** 2 for l in range(n_block))

    return S_block_contrib


# ----------- Closed-form components ----------------------------------------


def S_bare_per_site(N: int, gamma: float, t_array: np.ndarray) -> np.ndarray:
    """S_bare(t; N) per bare site (closed form): (N-1)/N² · exp(-4γt)."""
    return (N - 1) / (N * N) * np.exp(-4.0 * gamma * t_array)


_PATH_K_CACHE: dict[tuple, np.ndarray] = {}


def S_path_k(k: int, N: int, J: float, gamma: float, t_array: np.ndarray) -> np.ndarray:
    """Pure-topology S_(k)(t) at N qubits.

    For k=1: analytic all-isolated closed form S_(1)^1, N(t).
    For k>=2: block contribution (numerical, cached) + (N − k − 1) bare contributions.
    """
    if k == 1:
        # All-isolated formula at m = 1: S_(1)(t) = [(N-1)/N + 4(N-2)(cos(4Jt)-1)/(N²(N-1))] · exp(-4γt)
        baseline = (N - 1) / N
        correction = 4.0 * (N - 2) * (np.cos(4.0 * J * t_array) - 1.0) / (N * N * (N - 1))
        return (baseline + correction) * np.exp(-4.0 * gamma * t_array)

    n_block = k + 1
    cache_key = (n_block, N, J, gamma, len(t_array), float(t_array[0]), float(t_array[-1]))
    if cache_key not in _PATH_K_CACHE:
        print(f"  [computing path-{k} block (d²={4**n_block})...]", flush=True)
        _PATH_K_CACHE[cache_key] = path_k_block_S(n_block, N, J, gamma, t_array)
    S_block = _PATH_K_CACHE[cache_key]
    n_bare = N - n_block
    return S_block + n_bare * S_bare_per_site(N, gamma, t_array)


# ----------- Topology classifier --------------------------------------------


def classify_topology(bonds: list[int]) -> tuple[int, ...]:
    """Bond set → sorted tuple of connected-component-lengths."""
    if not bonds:
        return ()
    sorted_bonds = sorted(bonds)
    components = []
    cur_len = 1
    for i in range(1, len(sorted_bonds)):
        if sorted_bonds[i] == sorted_bonds[i - 1] + 1:
            cur_len += 1
        else:
            components.append(cur_len)
            cur_len = 1
    components.append(cur_len)
    return tuple(sorted(components))


def S_topology_additive(topology: tuple[int, ...], N: int, J: float, gamma: float, t_array: np.ndarray) -> np.ndarray:
    """Apply the additive identity: S_T = Σ_i S_(k_i) − (m−1)·N·S_bare."""
    m = len(topology)
    S = np.zeros_like(t_array, dtype=float)
    for k_i in topology:
        S += S_path_k(k_i, N, J, gamma, t_array)
    S -= (m - 1) * N * S_bare_per_site(N, gamma, t_array)
    return S


# ----------- Main: verify against 14 N=7 CSVs -------------------------------


def parse_bonds_from_csv(name: str) -> list[int]:
    """Parse CSV name like 'N7_b0-1-3_J...csv' → bonds [0, 1, 3]."""
    bond_str = name.split("_b")[1].split("_")[0]
    return [int(b) for b in bond_str.split("-")]


def main() -> None:
    J, gamma = 0.075, 0.05
    N = 7
    skip_path_6 = "--with-path6" not in sys.argv  # path-6 = 16384-dim eigendecomp (~13 GB)
    print(f"# F89 mixed-topology additive identity verification, N={N}, J={J}, γ={gamma}\n")
    print("# Identity: S_T(t) = Σ_i S_(k_i)(t) − (m − 1)·N·S_bare(t; N)")
    print("# where m = number of blocks, k_i = path-length of block i.\n")
    if skip_path_6:
        print("# Skipping topology (6) (path-6 d²=16384 eigendecomp). Pass --with-path6 to include.\n")

    csvs = sorted(CSV_DIR.glob(f"N{N}_b*_J0.0750_gamma0.0500_probe-coherence.csv"))
    if skip_path_6:
        csvs = [c for c in csvs if classify_topology(parse_bonds_from_csv(c.name)) != (6,)]
    print(f"## Verifying {len(csvs)} CSVs across topology classes at N={N}\n")
    print("| CSV (bonds) | Topology | m blocks | max |diff| | mean |diff| |")
    print("|---|---|---|---|---|")

    by_topology = {}
    for csv in csvs:
        bonds = parse_bonds_from_csv(csv.name)
        topology = classify_topology(bonds)
        if topology not in by_topology:
            by_topology[topology] = []
        by_topology[topology].append(csv)

    overall_max_diff = 0.0
    for topology in sorted(by_topology, key=lambda t: (sum(t), len(t), t)):
        for csv in by_topology[topology]:
            data = np.loadtxt(csv, delimiter=",", skiprows=1)
            t_csv, S_csv = data[:, 0], data[:, -1]
            S_pred = S_topology_additive(topology, N, J, gamma, t_csv)
            diff = S_pred - S_csv
            max_d = np.max(np.abs(diff))
            mean_d = np.mean(np.abs(diff))
            overall_max_diff = max(overall_max_diff, max_d)
            bonds = parse_bonds_from_csv(csv.name)
            print(
                f"| b{'-'.join(str(b) for b in bonds)} | {topology} | "
                f"{len(topology)} | {max_d:.2e} | {mean_d:.2e} |"
            )

    print()
    print(f"## Overall max |diff| across all {len(csvs)} CSVs: {overall_max_diff:.3e}")
    print(f"# CSV write precision is ~5e-7. If all matches are at this floor,")
    print(f"# the additive identity is verified at machine precision.")


if __name__ == "__main__":
    main()
