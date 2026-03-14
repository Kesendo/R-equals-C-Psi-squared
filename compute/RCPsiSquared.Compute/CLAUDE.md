# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

C# compute engine for Liouvillian spectral analysis of small open quantum networks (2-7+ qubits) under Lindblad dynamics with dephasing noise. Part of the larger R=CΨ² research project investigating mirror symmetry of decay rate spectra.

The engine constructs the Lindblad superoperator as a dense d²×d² matrix (d=2^N), diagonalizes it, and analyzes the resulting eigenvalue spectrum for mirror symmetry and band structure.

## Build and Run

```bash
dotnet restore
dotnet run -c Release
```

No solution file — single project. No test framework (xUnit/NUnit); validation is embedded in `Program.cs` which runs benchmarks, topology surveys, and stress tests. Results are written to `simulations/results_csharp_compute.txt` (hardcoded absolute path in Program.cs).

## Dependencies

- **.NET 9.0** (target framework)
- **MathNet.Numerics 5.0.0** — dense linear algebra, eigendecomposition
- **MathNet.Numerics.Providers.MKL 5.0.0** — optional Intel MKL native LAPACK (5-10x speedup). Program.cs tries `Control.UseNativeMKL()` at startup and falls back gracefully.
- `AllowUnsafeBlocks` is enabled in the csproj.

## Architecture

Four domain modules + one runner, all in namespace `RCPsiSquared.Compute`:

- **PauliOps.cs** — Single-qubit Pauli matrices (I, X, Y, Z) as `Matrix<Complex>`. `At(op, qubit, n)` places an operator on a specific qubit via tensor products. `PauliBasis(n)` generates all 4^N Pauli strings.
- **Topology.cs** — Defines `Bond` record struct and topology generators: `Star`, `Chain`, `Ring`, `Complete`, `BinaryTree`. Each returns `Bond[]` with isotropic Heisenberg (X+Y+Z) coupling. `BuildHamiltonian()` constructs H from bonds.
- **Liouvillian.cs** — `Build(nQubits, bonds, gammaPerQubit)` constructs the full Lindblad superoperator: -i[H,ρ] + dephasing dissipators (√γ_k Z_k per qubit). `GetOscillatoryRates()` extracts sorted decay rates from eigenvalues with Im > threshold.
- **MirrorAnalysis.cs** — `CheckSymmetry(rates, center)` scores mirror quality around a given center. Returns (Score 0-1, matched pair count).
- **Program.cs** — Top-level runner (not a class): benchmark N=2-7, topology survey at N=4-6, non-uniform coupling/noise stress tests, N=8 OOM attempt.

## Key Physics Concepts

- The Liouvillian is a d²×d² superoperator (vectorized Lindblad equation). For N=7, this is 16384×16384 (~8GB RAM).
- Oscillatory rates = -Re(λ) for eigenvalues with nonzero Im part. These form the decay rate spectrum.
- Mirror symmetry center = Σγ_k (sum of all dephasing rates). The spectrum is symmetric around this point under pure dephasing.
- Boundary formula: rates range from 2γ to 2(N-1)γ for uniform γ.

## Memory Constraints

| N | Matrix | RAM needed |
|---|--------|------------|
| 6 | 4096² | ~1 GB |
| 7 | 16384² | ~8 GB |
| 8 | 65536² | ~64 GB (OOM on most machines) |

## Relationship to Python Simulations

The parent `simulations/` directory contains ~37 Python scripts that implement the same physics with numpy/scipy. The C# engine was built for performance at larger N. Both implementations should produce identical eigenvalue spectra for the same parameters. The Streamlit app in `simulations/app/` provides an interactive UI.
