# RCPsiSquared.Compute

C# compute engine for Liouvillian spectral analysis of open quantum networks (2-8 qubits). Constructs the Lindblad superoperator, diagonalizes it via native LAPACK, and verifies palindromic mirror symmetry of the decay rate spectrum.

## Requirements

- .NET 10.0 SDK
- Intel MKL (via MathNet NuGet, auto-restored)
- OpenBLAS 0.3.31 (native DLLs, manual setup — see below)
- ~8 GB RAM for N=7, ~73 GB for N=8

## Setup

### 1. Restore NuGet packages

```bash
cd compute/RCPsiSquared.Compute
dotnet restore
```

### 2. Download OpenBLAS native DLLs

The eigenvalue-only path (required for N=8) calls LAPACK `zgeev_` directly via OpenBLAS. Two builds are needed: LP64 for validation and ILP64 for N=8.

**Linux/macOS/Git Bash:**
```bash
cd compute/RCPsiSquared.Compute
mkdir -p native

# LP64 (32-bit integer LAPACK, for validation at N≤7)
curl -LO https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.31/OpenBLAS-0.3.31-x64.zip
unzip OpenBLAS-0.3.31-x64.zip win64/bin/libopenblas.dll -d .
cp win64/bin/libopenblas.dll native/libopenblas.dll

# ILP64 (64-bit integer LAPACK, required for N=8)
curl -LO https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.31/OpenBLAS-0.3.31-x64-64.zip
unzip OpenBLAS-0.3.31-x64-64.zip win64-64/bin/libopenblas.dll -d .
cp win64-64/bin/libopenblas.dll native/libopenblas64.dll

# Cleanup
rm -rf *.zip win64/ win64-64/
```

The `.csproj` copies both DLLs to the output directory on build. They are gitignored.

### 3. Build and validate

```bash
dotnet build -c Release

# Quick validation (~6 seconds): compares eigenvalues from three backends at N=5
dotnet run -c Release -- validate
```

Expected output:
```
Match: YES - eigenvalue-only path validated!
Match: YES - ILP64 path validated! N=8 is safe.
```

## Running

```bash
# Full suite: N=2-7 benchmark + topology survey + stress tests + N=8
# Warning: N=7 alone takes ~92 minutes
dotnet run -c Release

# N=8 only (skip N=2-7): needs ~73 GB free RAM, takes ~10 hours
# Close all other applications first
dotnet run -c Release -- n8

# Validation only (~6 seconds)
dotnet run -c Release -- validate
```

Results are written to `simulations/results/csharp_compute.txt`.

## Architecture

| File | Purpose |
|------|---------|
| PauliOps.cs | Pauli matrices, tensor products, N-qubit Pauli basis |
| Topology.cs | Star, Chain, Ring, Complete, Tree bond generators |
| Liouvillian.cs | Lindblad superoperator: three build paths by N |
| MklDirect.cs | Direct LAPACK P/Invoke: LP64 + ILP64 with backend auto-detection |
| MirrorAnalysis.cs | Mirror symmetry scoring and spectral statistics |
| Program.cs | Benchmark, topology survey, stress tests, validation, N=8 |

## Three compute paths

| N | Matrix | Build path | Eigen path | Memory |
|---|--------|-----------|-----------|--------|
| 2-6 | ≤4096² | `Build()` MathNet Kronecker | MathNet `Evd()` via MKL | ≤1 GB |
| 7 | 16384² | `BuildDirectRaw()` element-wise | MKL `z_eigen` (with eigenvectors) | ~8 GB |
| 8 | 65536² | `BuildDirectNative()` parallel, native memory | OpenBLAS `zgeev_` ILP64 (eigenvalues only) | ~73 GB |

### Why N=8 needs special handling

1. **Array limit**: 65536² = 4.29B elements > .NET `int.MaxValue`. Solved with `NativeMemory.AllocZeroed`.
2. **LAPACK int overflow**: internal offset `lda×n` = 4.29B > `int32`. Solved with ILP64 OpenBLAS (`USE64BITINT`).
3. **Eigenvector memory**: MathNet's `z_eigen` always computes eigenvectors (3×64 GB = 192 GB). Solved with direct `zgeev_` JOBVL='N',JOBVR='N' (64 GB + 1 GB workspace).

## Performance vs Python

With MKL active, expect 5-10x speedup for eigendecomposition and 10-50x for matrix construction. N=7 that took 18 minutes in Python completes in ~92 minutes in C# (larger matrix via direct build). N=8 took 10.6 hours via OpenBLAS ILP64 on 24 cores.

All timings measured on Intel Core Ultra 9 285k (24 cores), 128 GB RAM, Windows 11.
