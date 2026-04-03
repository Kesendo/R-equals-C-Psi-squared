# RCPsiSquared.Propagate

C# time-domain propagation engine for Lindblad master equation dynamics on N-qubit systems (tested through N=13, targeting N=15). Integrates drho/dt via RK4 with zero-allocation hot loop, computing mutual information, CΨ, concurrence, and Pauli correlators at specified measurement times.

## What it does

Where `RCPsiSquared.Compute` diagonalizes the Liouvillian to get the *spectrum*, this engine *propagates* density matrices in time to observe the *dynamics*. Three modes serve different questions:

- **Profile mode** (March 2026): Evaluate any spatial dephasing profile on a Heisenberg chain. This is the engine behind the sacrifice-zone formula scaling results (N=5 through N=13, targeting N=15).
- **Default mode**: Mediator bridge topology tests - cross-bridge information flow, coupling/noise sweeps, standing wave correlators (N=5 and N=11).
- **Pull mode**: Scaling curves, coupling optimization, and the relay protocol (+83% MI improvement, N=3 through N=11).

## Requirements

- .NET 10.0 SDK
- Intel MKL (via MathNet NuGet, auto-restored)
- RAM scales with N (see table below)

| N | Density matrix | RAM (approx) | Runtime (t=20, dt=0.05) |
|---|---------------|-------------|-------------------------|
| 5 | 32x32 | <1 MB | <1s |
| 7 | 128x128 | ~10 MB | ~5s |
| 9 | 512x512 | ~200 MB | ~2 min |
| 11 | 2048x2048 | ~4 GB | ~10 min |
| 13 | 8192x8192 | ~16 GB | 1-6 hours (profile-dependent) |
| 15 | 32768x32768 | ~64 GB | ~10+ hours (estimated) |

## Setup

```
cd compute/RCPsiSquared.Propagate
dotnet restore
dotnet build -c Release
```

No native DLL setup needed (unlike RCPsiSquared.Compute). MKL handles the matrix multiplications via NuGet.

## Running

```
# Default: Mediator bridge tests (Test 0-3)
dotnet run -c Release

# Pull principle tests (scaling curve, 2:1 optimization, relay protocol)
dotnet run -c Release -- pull

# Profile mode: evaluate a single gamma profile (sacrifice-zone formula, optimizer output, etc.)
dotnet run -c Release -- profile <N> <g1,g2,...,gN> [--tmax 20] [--dt 0.05]
```

Default and pull results are written to `simulations/results/mediator_bridge_scale.txt` or `simulations/results/pull_principle.txt`. Profile mode writes to stdout only (single machine-parseable RESULT line).

### Profile mode examples

```
# N=7 sacrifice-zone formula: gamma_edge = 7*0.05 - 6*0.001 = 0.344
dotnet run -c Release -- profile 7 0.001,0.001,0.001,0.001,0.001,0.001,0.344

# N=7 V-shape baseline
dotnet run -c Release -- profile 7 0.080,0.070,0.060,0.050,0.060,0.070,0.080

# N=13 sacrifice-zone formula: gamma_edge = 13*0.05 - 12*0.001 = 0.638
dotnet run -c Release -- profile 13 0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.638
```

Output format:
```
RESULT SumMI=0.408000 PeakMI=0.012345 PeakT=3.50 CPsi01=0.187654 Purity=0.001234 SumMI5=0.350000
```

Fields: SumMI (sum of MI for all adjacent pairs at best time), PeakMI (peak MI between endpoints 0 and N-1), PeakT (time of peak SumMI), CPsi01 (CΨ of qubit pair 0-1 at peak), Purity (global purity at peak), SumMI5 (SumMI at t=5.0).

Initial state is always |+>^N. Product states are the optimal choice because each qubit must respond independently to its local dephasing rate (GHZ states are completely blind to spatial dephasing profiles, see [gamma as Signal](../../experiments/GAMMA_AS_SIGNAL.md)). Hamiltonian is a Heisenberg chain with J=1.0 on all bonds.

**Note:** Only one profile evaluation can run at a time (the executable locks while running). Plan batch runs sequentially.

## Test suite

### Default mode (mediator bridge)

| Test | What it measures | N |
|------|-----------------|---|
| Test 0 | Validation: C# vs known values (Level 2, 5-qubit, no Python needed) | 5 |
| Test 1 | Cross-bridge information flow: MI(BridgeA:BridgeB) at Level 3 | 11 |
| Test 2a | Meta-mediator coupling sweep (J_meta = 0.25 to 3.0) | 11 |
| Test 2b | Meta-mediator noise sweep (γ_M = 0.0 to 0.5) | 11 |
| Test 3 | Standing wave correlators (XX, ZZ) across the full bridge | 11 |

### Pull mode (`-- pull`)

| Test | What it measures | N |
|------|-----------------|---|
| Part 1 | Scaling curve: MI vs N for uniform chain vs hierarchical mediator | 3-11 |
| Part 2 | Pull optimization: 5 coupling configurations (symmetric, 2:1 bridges, 2:1 meta, 2:1 all, reversed) | 11 |
| Part 3 | Relay protocol: passive vs staged γ vs staged + 2:1 pull | 11 |

## Architecture

| File | Purpose |
|------|---------|
| LindbladPropagator.cs | RK4 integrator with zero-allocation hot loop. Precomputed dephasing mask. 2 BLAS calls per RHS evaluation. |
| DensityMatrixTools.cs | Partial trace, mutual information, CΨ, concurrence, purity, expectation values |
| PauliOps.cs | Pauli matrices, tensor products, N-qubit operator construction |
| Topology.cs | Chain, MediatorBridge(level), bond generators with configurable J |
| Program.cs | Test dispatch: profile mode (single evaluation), mediator bridge suite, pull principle suite |

## Key design decisions

**Zero-allocation RK4.** All six work matrices (k1-k4, tmp, drho) are pre-allocated in the constructor. The hot loop runs `Parallel.For` over flat arrays with no GC pressure. This matters at N=11 where each density matrix is 2048×2048 = 4M complex entries.

**Precomputed dephasing mask.** The Z-dephasing dissipator reduces to element-wise multiplication: `drho[i,j] += mask[i,j] * rho[i,j]` where mask depends only on the XOR of basis indices. Computed once, stored as flat double array in column-major order.

**Staged propagation for relay protocol.** The relay protocol changes γ profiles between stages. Each stage builds a new `LindbladPropagator` with different gammas, then propagates for one stage duration. The density matrix carries over between stages.

## Topology: The hierarchical mediator bridge

`Topology.MediatorBridge(level)` builds a self-similar bridge topology:

```
Level 2 (N=5):   (0-1) - M(2) - (3-4)
Level 3 (N=11):  (0-1) - m1(2) - (3-4) - M(5) - (6-7) - m2(8) - (9-10)
```

Each level wraps the previous in a pair-mediator-pair structure. Level 3 connects two Level 2 bridges through a meta-mediator M(5).

## Relationship to experiments

| Experiment document | Uses Propagate for |
|--------------------|--------------------|
| [Resonant Return](../../experiments/RESONANT_RETURN.md) | Sacrifice-zone formula validation (profile mode), N=5 through N=13 |
| [Signal Analysis: Scaling](../../experiments/SIGNAL_ANALYSIS_SCALING.md) | Formula scaling curve, quadratic growth analysis |
| [IBM Sacrifice Zone](../../experiments/IBM_SACRIFICE_ZONE.md) | Simulation baselines for IBM hardware comparison |
| [Relay Protocol](../../experiments/RELAY_PROTOCOL.md) | Staged γ relay, +83% MI improvement |
| [Scaling Curve](../../experiments/SCALING_CURVE.md) | MI vs N, hierarchical vs uniform |
| [Star Topology Observers](../../experiments/STAR_TOPOLOGY_OBSERVERS.md) | Entanglement echo, Bohr frequencies |

## See also

- [RCPsiSquared.Compute](../RCPsiSquared.Compute/README.md): Spectral analysis engine (diagonalization, palindrome verification, N=2 through N=8)
- [Engineering Blueprint](../../publications/ENGINEERING_BLUEPRINT.md): Seven design rules derived from Propagate results
- [Complete Mathematical Documentation](../../docs/proofs/COMPLETE_MATHEMATICAL_DOCUMENTATION.md): The Tafelwerk
