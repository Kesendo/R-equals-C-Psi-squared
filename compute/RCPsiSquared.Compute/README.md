# RCPsiSquared.Compute

C# compute engine for Liouvillian spectral analysis of open quantum networks (2-8 qubits). Constructs the Lindblad superoperator, diagonalizes it via native LAPACK, and verifies palindromic mirror symmetry of the decay rate spectrum.

## What it does

This engine answered the central question of the R=CPsi^2 project: is the Liouvillian decay spectrum exactly palindromic? It diagonalizes the full superoperator for N=2 through N=8 (matrices up to 65536x65536, 87,376 eigenvalues total) and confirms 100% palindromic symmetry across all tested topologies, coupling models, and noise configurations. Zero exceptions.

The results feed directly into the [Mirror Symmetry Proof](../../docs/proofs/MIRROR_SYMMETRY_PROOF.md), which provides the analytical explanation via the conjugation operator Pi. This engine supplied the numerical verification that motivated and validated that proof.

For time-domain dynamics (information flow, relay protocols, sacrifice-zone formula scaling), see [RCPsiSquared.Propagate](../RCPsiSquared.Propagate/README.md).

## Requirements

- .NET 10.0 SDK
- Intel MKL (via MathNet NuGet, auto-restored)
- OpenBLAS 0.3.31 (native DLLs, manual setup; see below)
- ~8 GB RAM for N=7, ~73 GB for N=8

## Setup

### 1. Restore NuGet packages

```
cd compute/RCPsiSquared.Compute
dotnet restore
```

### 2. Download OpenBLAS native DLLs

The eigenvalue-only path (required for N=8) calls LAPACK `zgeev_` directly via OpenBLAS. Two builds are needed: LP64 for validation and ILP64 for N=8.

**Windows (PowerShell):**
```powershell
cd compute\RCPsiSquared.Compute
New-Item -ItemType Directory -Force native

# LP64 (32-bit integer LAPACK, for validation at N<=7)
Invoke-WebRequest -Uri "https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.31/OpenBLAS-0.3.31-x64.zip" -OutFile OpenBLAS-LP64.zip
Expand-Archive OpenBLAS-LP64.zip -DestinationPath temp-lp64
Copy-Item temp-lp64\win64\bin\libopenblas.dll native\libopenblas.dll

# ILP64 (64-bit integer LAPACK, required for N=8)
Invoke-WebRequest -Uri "https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.31/OpenBLAS-0.3.31-x64-64.zip" -OutFile OpenBLAS-ILP64.zip
Expand-Archive OpenBLAS-ILP64.zip -DestinationPath temp-ilp64
Copy-Item temp-ilp64\win64-64\bin\libopenblas.dll native\libopenblas64.dll

# Cleanup
Remove-Item OpenBLAS-LP64.zip, OpenBLAS-ILP64.zip, temp-lp64, temp-ilp64 -Recurse -Force
```

**Linux/macOS/Git Bash (Windows DLLs only - see note):**
```
cd compute/RCPsiSquared.Compute
mkdir -p native

# LP64 (32-bit integer LAPACK, for validation at N<=7)
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

**Note:** The native DLLs are Windows binaries. The Linux/macOS instructions above are for Git Bash on Windows. Running on native Linux/macOS would require building OpenBLAS from source with the corresponding shared library names (libopenblas.so / libopenblas.dylib). This has not been tested.

The `.csproj` copies both DLLs to the output directory on build. They are gitignored.

### 3. Build and validate

```
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

```
# Full suite: N=2-7 benchmark + topology survey + stress tests + N=8
# Warning: N=7 alone takes ~92 minutes
dotnet run -c Release

# N=8 only (skip N=2-7): needs ~73 GB free RAM, takes ~10 hours
# Close all other applications first
dotnet run -c Release -- n8

# Validation only (~6 seconds)
dotnet run -c Release -- validate

# RMT eigenvalue export: all complex eigenvalues as CSV (N=2-7)
# N=2-6 in ~1 min, N=7 in ~95 min
dotnet run -c Release -- rmt

# Cavity modes at zero noise (N=2-7, ~40 min for N=7)
dotnet run -c Release -- cavity

# Cavity topology/coupling tests (Ring, Complete, non-uniform J)
dotnet run -c Release -- cavity tests
```

Results are written to:
- `simulations/results/csharp_compute.txt` (main suite)
- `simulations/results/rmt_eigenvalues_N{2..7}.csv` (RMT export)
- `simulations/results/cavity_modes_zero_noise.txt` (cavity modes)
- `simulations/results/cavity_modes_tests.txt` (cavity tests)

## Architecture

| File | Purpose |
|------|---------|
| PauliOps.cs | Pauli matrices, tensor products, N-qubit Pauli basis |
| Topology.cs | Star, Chain, Ring, Complete, Tree bond generators |
| Liouvillian.cs | Lindblad superoperator: three build paths by N. `GetAllEigenvalues()` for RMT export. `GetCavityModes()` for zero-noise eigenfrequency analysis |
| MklDirect.cs | Direct LAPACK P/Invoke: LP64 + ILP64 with backend auto-detection |
| MirrorAnalysis.cs | Mirror symmetry scoring and spectral statistics |
| Program.cs | Benchmark, topology survey, stress tests, validation, N=8, RMT export |

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

## Performance

The C# engine's main advantage is matrix construction: `BuildDirectRaw` (element-wise, no Kronecker products) is orders of magnitude faster than Python's numpy Kronecker approach. Eigendecomposition speed is comparable since both call LAPACK under the hood.

| N | Build (C#) | Eigen (C#) | Engine |
|---|------------|------------|--------|
| 6 | 8.7s | 56s | MathNet + MKL |
| 7 | 0.1s | 92min | Direct build + MKL z_eigen |
| 8 | 5.6s | 10.6h | Native memory + OpenBLAS ILP64 zgeev |

N=7 and N=8 are only feasible in C# due to the direct build path (no Kronecker, no 2GB .NET limit). Python's numpy/scipy cannot build the 16384x16384 or 65536x65536 Liouvillian in reasonable time or memory.

All timings measured on Intel Core Ultra 9 285k (24 cores), 128 GB RAM, Windows 11.

## Relationship to experiments

| Experiment document | Uses Compute for |
|--------------------|-----------------|
| [Non-Heisenberg Palindrome](../../experiments/NON_HEISENBERG_PALINDROME.md) | Eigendecomposition across XY, Ising, XXZ, DM models. Two Pi families (P1, P4) |
| [XOR Space](../../experiments/XOR_SPACE.md) | Eigenvalue-to-mode mapping. GHZ vs W initial state spectral decomposition |
| [Standing Wave Analysis](../../experiments/STANDING_WAVE_ANALYSIS.md) | Spectral structure behind XX/YY oscillation and ZZZ static behavior |
| [Structural Cartography](../../experiments/STRUCTURAL_CARTOGRAPHY.md) | 3D spectral manifold, rate count statistics, topology survey |
| [Crossing Taxonomy](../../experiments/CROSSING_TAXONOMY.md) | K-invariance from Lindblad eigenvalue scaling |
| [Cavity Modes Formula](../../experiments/CAVITY_MODES_FORMULA.md) | Zero-noise eigenfrequencies N=2-7. Closed-form via Clebsch-Gordan. Topology comparison |
| [IBM Cavity Spectral](../../experiments/IBM_CAVITY_SPECTRAL_ANALYSIS.md) | Sacrifice zone protects cavity modes at 2.81x. Real IBM T2* data |
| [Cavity Mode Localization](../../experiments/CAVITY_MODE_LOCALIZATION.md) | Eigenvector Pauli decomposition: protected modes are center-localized (r = 0.994) |
| [Random Matrix Theory](../../experiments/RANDOM_MATRIX_THEORY.md) | RMT analysis of 21,832 eigenvalues (N=2-7): Poisson level statistics, chiral symmetry class AIII. Uses `--rmt` CSV export |

## See also

- [RCPsiSquared.Propagate](../RCPsiSquared.Propagate/README.md): Time-domain engine for dynamics, sacrifice-zone formula, relay protocol (N=5 through N=15)
- [Mirror Symmetry Proof](../../docs/proofs/MIRROR_SYMMETRY_PROOF.md): The analytical proof that this engine's numerical results motivated
- [Complete Mathematical Documentation](../../docs/proofs/COMPLETE_MATHEMATICAL_DOCUMENTATION.md): The Tafelwerk
- [Engineering Blueprint](../../publications/ENGINEERING_BLUEPRINT.md): Seven design rules derived from spectral analysis
