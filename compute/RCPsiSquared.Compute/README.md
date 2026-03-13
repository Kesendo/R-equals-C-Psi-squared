# RCPsiSquared.Compute

C# compute engine for Liouvillian spectral analysis.
Uses MathNet.Numerics with Intel MKL for native LAPACK performance.

## Requirements

- .NET 9.0 SDK
- ~8GB RAM for N=7, ~128GB for N=8

## Build and Run

```bash
cd compute/RCPsiSquared.Compute
dotnet restore
dotnet run -c Release
```

Results are written to `simulations/results_csharp_compute.txt`.

## Architecture

| File | Purpose |
|------|---------|
| PauliOps.cs | Pauli matrices, tensor products, N-qubit Pauli basis |
| Topology.cs | Star, Chain, Ring, Complete, Tree bond generators |
| Liouvillian.cs | Lindblad superoperator construction and eigendecomposition |
| MirrorAnalysis.cs | Mirror symmetry scoring and spectral statistics |
| Program.cs | Test runner: benchmark, topology survey, N=8 attempt |

## Performance vs Python

With MKL active, expect 5-10x speedup for eigendecomposition
and 10-50x for matrix construction (no Python interpreter overhead).
N=7 that took 18 minutes in Python should complete in 2-5 minutes.
