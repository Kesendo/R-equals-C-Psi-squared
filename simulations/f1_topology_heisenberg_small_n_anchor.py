"""Python small-N anchor for chain/ring/star Heisenberg F1 spectrum metrics.

Extends `f1_chain_heisenberg_small_n_anchor.py` to 3 topologies at N=3..6,
giving us cross-topology + cross-N data for:
  - The `gap × N² ≈ 2.20` chain-Heisenberg scaling discovered 2026-05-19
  - The "Im/σ = 1 for single-hub geometries" pattern (N=3 chain, all-N star)
  - MinReal = −2σ universal (5+ datapoints per topology)
  - Kernel dim = N+1 for connected, Π_c(|c|+1) for disconnected (Tier1Candidate)

Hamiltonian: H = (J/4) Σ_b (X_i X_j + Y_i Y_j + Z_i Z_j), J=1, γ=0.5 per site.
Matches the C# F1GeneralTopologyN{7,8,9}BlockSpectrumChainTests convention.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, "simulations")

from framework.lindblad import lindbladian_z_dephasing  # noqa: E402

from f1_chain_heisenberg_small_n_anchor import (  # noqa: E402
    chain_bonds,
    compute_palindromic_pairing_distances,
    compute_spectrum_structure,
    heisenberg_graph_h,
)


def ring_bonds(N: int) -> list[tuple[int, int]]:
    return chain_bonds(N) + [(N - 1, 0)]


def star_bonds(N: int) -> list[tuple[int, int]]:
    return [(0, i) for i in range(1, N)]


def run(topology_name: str, bonds: list[tuple[int, int]], N: int,
        J: float = 1.0, gamma: float = 0.5) -> dict:
    sigma = N * gamma
    H = heisenberg_graph_h(N, bonds, J=J)
    gammas = [gamma] * N
    t0 = time.time()
    L = lindbladian_z_dephasing(H, gammas)
    t_build = time.time() - t0
    t0 = time.time()
    spectrum = np.linalg.eigvals(L)
    t_eig = time.time() - t0

    palindromic = compute_palindromic_pairing_distances(spectrum, sigma)
    structure = compute_spectrum_structure(spectrum)

    return {
        "N": N,
        "TopologyName": f"{topology_name} ({len(bonds)} bonds)",
        "TotalWallSeconds": t_build + t_eig,
        "ComputeSpectrumWallSeconds": t_eig,
        **palindromic,
        **structure,
        "JValue": J,
        "GammaValue": gamma,
        "SigmaShift": -2.0 * sigma,
        "ComputedBy": "python_dense_eig_anchor",
    }


def main() -> None:
    print("F1 chain/ring/star Heisenberg small-N Python anchor")
    print("=" * 100)

    topologies = [
        ("chain", chain_bonds),
        ("ring", ring_bonds),
        ("star", star_bonds),
    ]
    Ns = [3, 4, 5, 6]

    all_results: dict[tuple[str, int], dict] = {}
    for topology, bonds_fn in topologies:
        for N in Ns:
            if topology == "ring" and N < 3:
                continue
            if topology == "star" and N < 2:
                continue
            bonds = bonds_fn(N)
            print(f"\n  {topology} N={N} ({len(bonds)} bonds)...", end=" ", flush=True)
            t0 = time.time()
            metrics = run(topology, bonds, N)
            print(f"{time.time()-t0:.1f}s wall, eig {metrics['ComputeSpectrumWallSeconds']:.1f}s")
            all_results[(topology, N)] = metrics

    out_dir = Path("simulations/results/f1_n8_n9_metrics")
    out_dir.mkdir(parents=True, exist_ok=True)
    for (topology, N), metrics in all_results.items():
        out_path = out_dir / f"{topology}_N{N}_python.json"
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(metrics, f, indent=2)

    print("\n" + "=" * 100)
    print("Summary table (J=1, γ=0.5):")
    print(f"{'Top':>6} {'N':>3} {'MinReal':>10} {'|Im|max':>9} {'Im/σ':>6} {'KernDim':>8} "
          f"{'N+1':>4} {'DissGap':>10} {'gap·N²':>8} {'DistBin':>8} {'MaxPair':>10}")
    for topology, _ in topologies:
        for N in Ns:
            if (topology, N) not in all_results:
                continue
            m = all_results[(topology, N)]
            sigma = N * 0.5
            print(f"{topology:>6} {N:3d} {m['MinReal']:10.5f} {m['MaxImag']:9.4f} "
                  f"{m['MaxImag']/sigma:6.4f} {m['KernelDimension']:8d} {N+1:4d} "
                  f"{m['DissipationGap']:10.5f} {m['DissipationGap']*N*N:8.4f} "
                  f"{m['DistinctBinnedEigenvalueCount']:8d} {m['MaxPairingDistance']:10.3e}")


if __name__ == "__main__":
    main()
