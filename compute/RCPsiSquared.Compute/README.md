# RCPsiSquared.Compute

C# compute engine for Liouvillian spectral analysis.
Uses MathNet.Numerics with Intel MKL for native LAPACK performance.

## Requirements

- .NET 10.0 SDK
- ~8 GB RAM for N=7, ~128 GB for N=8

## Build and Run

```bash
cd compute/RCPsiSquared.Compute
dotnet run -c Release
```

Results are written to `simulations/results/csharp_compute.txt`.

## Architecture

| File | Purpose |
|------|---------|
| PauliOps.cs | Pauli matrices, tensor products, N-qubit Pauli basis |
| Topology.cs | Star, Chain, Ring, Complete, Tree bond generators |
| Liouvillian.cs | Lindblad superoperator: KroneckerProduct build (N<=6) and direct raw build (N>=7) |
| MklDirect.cs | Direct MKL P/Invoke with pinned arrays, bypasses .NET 2 GB marshaling limit |
| MirrorAnalysis.cs | Mirror symmetry scoring and spectral statistics |
| Program.cs | Benchmark N=2-7, topology survey, stress tests |

## Two build paths

- **N <= 6**: `Liouvillian.Build()` uses MathNet KroneckerProduct + standard MKL Evd. Fast and simple.
- **N >= 7**: `Liouvillian.BuildDirectRaw()` fills the superoperator element-wise into a raw `Complex[]`, then `MklDirect.EigenvaluesRaw()` calls MKL's `z_eigen` directly with pinned pointers. This bypasses the .NET P/Invoke 2 GB array marshaling limit and uses all CPU cores.

## Performance vs Python

With MKL active, expect 5-10x speedup for eigendecomposition
and 10-50x for matrix construction (no Python interpreter overhead).
N=7 that took 18 minutes in Python should complete in 2-5 minutes.
