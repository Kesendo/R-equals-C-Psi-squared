"""Python small-N anchor for chain Heisenberg F1 spectrum metrics (N=3..6).

Computes the same metric structure as the C# F1SpectrumStatistics for direct
comparison against the C# block-spectrum infrastructure (N=7,8,9). Validates
that small-N spectrum properties scale predictably and gives us an independent
implementation cross-check on the C# MklDirect bridge path.

Hamiltonian: same as the C# N=7,8,9 tests:
  H = (J/4) Σ_b (X_i X_j + Y_i Y_j + Z_i Z_j)   J = 1
  uniform Z-dephasing γ = 0.5 per site
  σ = N·γ

The N=3..6 chain spectra are computed via direct dense numpy eig on the full
Liouvillian L = -i[H, ·] + Σ_l γ_l (Z_l ρ Z_l - ρ) with framework's
lindbladian_z_dephasing primitive. Dense at N=6 is 4096² ≈ 256 MB, ~30s total.

Output: simulations/results/f1_n8_n9_metrics/chain_N{3,4,5,6}_python.json,
matching F1SpectrumStatistics.TopologyMetrics record structure as closely as
possible (the C# pairing tolerance fields and the block-structure fields are
omitted since this is dense eig, no block decomposition).
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

import framework as fw  # noqa: E402
from framework.lindblad import lindbladian_z_dephasing  # noqa: E402


def chain_bonds(N: int) -> list[tuple[int, int]]:
    """Bond list (i, i+1) for the open chain on N sites."""
    return [(i, i + 1) for i in range(N - 1)]


def heisenberg_graph_h(N: int, bonds: list[tuple[int, int]], J: float = 1.0) -> np.ndarray:
    """H = (J/4) Σ_b (X_i X_j + Y_i Y_j + Z_i Z_j) on the given bond list.

    Matches the C# F1GeneralTopologyN{7,8,9}BlockSpectrumChainTests convention.
    Shared by the chain anchor script and the cross-topology extension
    `_f1_topology_heisenberg_small_n_anchor.py`.
    """
    terms = [("X", "X", J / 4.0), ("Y", "Y", J / 4.0), ("Z", "Z", J / 4.0)]
    return fw._build_bilinear(N, bonds, terms)


def heisenberg_chain_h(N: int, J: float = 1.0) -> np.ndarray:
    """Convenience wrapper: Heisenberg H on the open chain at N."""
    return heisenberg_graph_h(N, chain_bonds(N), J=J)


def compute_palindromic_pairing_distances(
    spectrum: np.ndarray, sigma: float
) -> dict[str, float | int]:
    """Match each eigenvalue λ to its nearest neighbour in {−2σ − λ'}.

    Mirror of F1SpectrumStatistics.ComputePalindromicMetrics: sort both lists
    by real part, walk in O(n) pairing each λ with the closest mirror.
    """
    n = len(spectrum)
    target = -2.0 * sigma - spectrum
    idx_spec = np.argsort(spectrum.real * 1e6 + spectrum.imag)
    idx_targ = np.argsort(target.real * 1e6 + target.imag)
    sorted_spec = spectrum[idx_spec]
    sorted_targ = target[idx_targ]
    distances = np.abs(sorted_spec - sorted_targ)
    finite = distances[np.isfinite(distances)]
    median = float(np.median(finite))
    outlier_count = int(np.sum(finite > 100.0 * max(median, 1e-300)))
    return {
        "MaxPairingDistance": float(np.max(finite)),
        "MeanPairingDistance": float(np.mean(finite)),
        "MedianPairingDistance": median,
        "P99PairingDistance": float(np.percentile(finite, 99)),
        "MinPairingDistance": float(np.min(finite)),
        "OutlierPairCount": outlier_count,
    }


def compute_spectrum_structure(
    spectrum: np.ndarray, kernel_tol: float = 1e-9
) -> dict[str, float | int]:
    re = spectrum.real
    im = spectrum.imag
    abs_re = np.abs(re)
    abs_im = np.abs(im)
    abs_lam = np.abs(spectrum)

    kernel_mask = abs_lam < kernel_tol
    real_mask = abs_im < kernel_tol
    pure_imag_mask = (abs_re < kernel_tol) & (abs_im >= kernel_tol)
    non_kernel_decay = abs_re[(abs_re > kernel_tol)]
    diss_gap = float(np.min(non_kernel_decay)) if non_kernel_decay.size > 0 else float("nan")

    # Binned distinct count to 1e-9 precision
    binned = np.round(re * 1e9).astype(np.int64) * 10_000_000_000 + np.round(im * 1e9).astype(np.int64)
    distinct = int(len(np.unique(binned)))

    return {
        "SpectrumSize": int(len(spectrum)),
        "MinReal": float(np.min(re)),
        "MaxReal": float(np.max(re)),
        "MinImag": float(np.min(im)),
        "MaxImag": float(np.max(im)),
        "DissipationGap": diss_gap,
        "KernelDimension": int(np.sum(kernel_mask)),
        "PureImaginaryCount": int(np.sum(pure_imag_mask)),
        "RealEigenvalueCount": int(np.sum(real_mask)),
        "DistinctBinnedEigenvalueCount": distinct,
    }


def run_chain_N(N: int, J: float = 1.0, gamma: float = 0.5) -> dict:
    """Compute the full F1 spectrum metrics for chain Heisenberg at given N."""
    d = 2 ** N
    print(f"\n--- N={N} chain Heisenberg (4^{N} = {d**2} Liouvillian dim) ---")
    sigma = N * gamma
    print(f"σ = N·γ = {N}·{gamma} = {sigma}")
    print(f"σ_shift convention: −2σ = {-2*sigma}")

    t0 = time.time()
    H = heisenberg_chain_h(N, J=J)
    gammas = [gamma] * N
    L = lindbladian_z_dephasing(H, gammas)
    t_build = time.time() - t0

    print(f"H built and L assembled: {t_build:.2f}s, L shape {L.shape}")
    t0 = time.time()
    spectrum = np.linalg.eigvals(L)
    t_eig = time.time() - t0
    print(f"Eig: {t_eig:.2f}s, {len(spectrum)} eigenvalues")

    palindromic = compute_palindromic_pairing_distances(spectrum, sigma)
    structure = compute_spectrum_structure(spectrum)
    print(f"  MinReal = {structure['MinReal']:+.6e} (predicted −2σ = {-2*sigma})")
    print(f"  MaxReal = {structure['MaxReal']:+.6e}")
    print(f"  Im range = ±{structure['MaxImag']:.6f}")
    print(f"  Im / σ = {structure['MaxImag']/sigma:.6f}")
    print(f"  Kernel dim = {structure['KernelDimension']} (predicted N+1 = {N+1})")
    print(f"  Dissipation gap = {structure['DissipationGap']:.6e}")
    print(f"  Distinct binned λ = {structure['DistinctBinnedEigenvalueCount']}")
    print(f"  Max palindromic pair dist = {palindromic['MaxPairingDistance']:.3e}")

    return {
        "N": N,
        "TopologyName": f"chain ({N-1} bonds)",
        "TotalWallSeconds": t_build + t_eig,
        "ComputeSpectrumWallSeconds": t_eig,
        "EffectiveSpeedupOverDense": 1.0,  # dense path, no speedup
        **palindromic,
        **structure,
        # block-structure fields N/A for dense eig
        "SectorCount": -1,
        "PrimarySectorCount": -1,
        "MaxBlockSize": -1,
        "MaxBlockSectorPCol": -1,
        "MaxBlockSectorPRow": -1,
        "Top3BlockSizes": [],
        "TotalBlockCubicCost": -1,
        "JValue": J,
        "GammaValue": gamma,
        "SigmaShift": -2.0 * sigma,
        "HamiltonianClass": f"Heisenberg XXX (XX+YY+ZZ) at J={J}, uniform Z-dephasing γ={gamma}",
        "Bonds": [{"Site1": i, "Site2": i + 1, "J": 1.0} for i in range(N - 1)],
        "ComputedBy": "python_dense_eig_anchor",
    }


def main() -> None:
    print("F1 chain Heisenberg small-N Python anchor")
    print("=" * 78)
    print("Cross-check anchor for the C# block-spectrum N=7,8,9 work.")
    print()

    results = {}
    for N in (3, 4, 5, 6):
        metrics = run_chain_N(N)
        results[N] = metrics

    out_dir = Path("simulations/results/f1_n8_n9_metrics")
    out_dir.mkdir(parents=True, exist_ok=True)
    for N, metrics in results.items():
        out_path = out_dir / f"chain_N{N}_python.json"
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(metrics, f, indent=2)
        print(f"Wrote {out_path}")

    print("\n" + "=" * 78)
    print("Summary table (chain Heisenberg J=1, γ=0.5):")
    print(f"{'N':>3}  {'MinReal':>12}  {'-2σ pred':>12}  {'|Im max|':>9}  {'Im/σ':>6}  "
          f"{'KernDim':>7}  {'N+1':>4}  {'DissGap':>10}  {'MaxPair':>10}")
    for N, m in results.items():
        sigma = N * 0.5
        print(f"{N:3d}  {m['MinReal']:12.6f}  {-2*sigma:12.6f}  "
              f"{m['MaxImag']:9.4f}  {m['MaxImag']/sigma:6.4f}  "
              f"{m['KernelDimension']:7d}  {N+1:4d}  "
              f"{m['DissipationGap']:10.6f}  {m['MaxPairingDistance']:10.3e}")


if __name__ == "__main__":
    main()
