"""F89 mixed topology: additive identity verification across N=7 CSVs.

S_T(t) = Σ_i S_(k_i)(t) − (m − 1) · N · S_bare(t; N) for any T = (k_1, ..., k_m).

See `experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md` § "Mixed-topology additive
identity" for the derivation from Lindbladian factorisation.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from _f89_pathk_lib import (
    build_block_L,
    compute_rho_block_0,
    reduce_block_to_site_01,
)

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO = Path(__file__).resolve().parents[1]
CSV_DIR = REPO / "simulations" / "results" / "bond_isolate"


def path_k_block_S(n_block: int, N: int, J: float, gamma: float, t_array: np.ndarray) -> np.ndarray:
    """Per-block S(t) contribution from a n_block-site path block in N qubits."""
    D = 2**n_block
    L = build_block_L(J, gamma, n_block)
    rho = compute_rho_block_0(n_block, N)
    vec = rho.flatten(order="F")
    eigvals, R = np.linalg.eig(L)
    c = np.linalg.solve(R, vec)

    S_block_contrib = np.zeros_like(t_array, dtype=float)
    for ti, t in enumerate(t_array):
        vec_t = R @ (np.exp(eigvals * t) * c)
        rho_t = vec_t.reshape((D, D), order="F")
        S_block_contrib[ti] = sum(
            2.0 * abs(reduce_block_to_site_01(rho_t, l, n_block)) ** 2
            for l in range(n_block)
        )

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
