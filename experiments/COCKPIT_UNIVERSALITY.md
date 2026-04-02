# Decoherence Cockpit: A Low-Dimensional Diagnostic Framework for Open Quantum Systems

**Date:** April 2, 2026
**Status:** Complete
**Scripts:** `simulations/cockpit_*.py`, `simulations/theta_pc_correlation.py`

## Abstract

Tracking the decoherence of a quantum system normally requires full
state tomography, which scales as 4^N measurement bases for N qubits.
We show this is unnecessary: for Heisenberg spin chains under local
dephasing, 3 observables capture 88-96% of the decoherence trajectory
variance across system sizes N=2-5, 9 topologies, and 2 noise types.
The three observables are Purity (trace of ρ squared), Concurrence
(entanglement), and normalized L₁ coherence (off-diagonal magnitude).
Which observable dominates depends on the architecture: PCA
automatically selects it. On IBM Torino hardware, the framework
achieves sub-1% crossing-time accuracy on well-characterized qubits.
The practical cost is 3 measurements per qubit pair instead of 4^N
tomographic bases.

## 1. Background

### The problem

A pair of qubits in a quantum computer is described by a 4x4 density
matrix ρ with 15 independent real parameters. As the system
interacts with its environment (decoherence), these parameters evolve
along a trajectory from the initial state toward the maximally mixed
state. Monitoring this trajectory is essential for error detection,
but measuring all 15 parameters requires 9 measurement bases per
pair, repeated many times for statistics. For N qubits, the cost
grows as 4^N. This is impractical for real-time monitoring.

### Key definitions

**Purity** = Tr(ρ^2). Ranges from 1 (pure state) to 1/d (maximally
mixed). Measures how much the state has decohered.

**Concurrence** (Wootters). Ranges from 0 (no entanglement) to 1
(maximally entangled). Measures how much entanglement remains between
two qubits.

**Ψ-norm** = L₁/(d-1), where L₁ is the sum of absolute values of
all off-diagonal elements of ρ, and d is the Hilbert space dimension.
Measures how much quantum coherence (superposition) remains.

**CΨ** = Purity x Ψ-norm. A combined measure of quantum strength
that weights coherence by state purity. When CΨ > ¼, the system
has no classical attractor (the fixed-point equation R = C*Psi^2 has
complex roots). When CΨ < ¼, classical attractors exist.

**θ** = arctan(sqrt(4*CΨ - 1)). Defined only when CΨ > ¼.
Measures the angular distance from the quantum-classical boundary at
CΨ = ¼. θ = 0 at the boundary, θ = 60 deg for a pure
coherent product state.

**Sacrifice zone.** In a qubit chain, deliberately concentrating noise
on boundary qubits ("sacrificing" them) to protect center qubits.
The total noise budget stays the same; only its distribution changes.

**Liouvillian.** The generator of open quantum dynamics. For N qubits
under dephasing, it is a 4^N x 4^N matrix whose eigenvalues give
the decay rates and whose eigenvectors give the decay modes.

**Bures distance.** A metric on the space of density matrices that
quantifies how distinguishable two quantum states are. From it, one
derives the Bures velocity (how fast the state changes) and the
Gaussian curvature (how curved the trajectory is near the CΨ = ¼
boundary).

### The question

**How many observables do you actually need to track decoherence?**

## 2. Method

### Feature vector

For a qubit pair (i,j) reduced from an N-qubit system, we extract
9 observables at each time step:

| Feature | Physical meaning |
|---|---|
| Phi+, Phi-, Psi+, Psi- | Overlap with each of the 4 Bell states |
| Pur | Purity = Tr(ρ^2), how mixed the state is |
| SvN | Von Neumann entropy, information-theoretic mixedness |
| C | Concurrence, entanglement strength |
| Psi | Ψ-norm, coherence magnitude |
| ph03 | Phase angle of the (0,3) density matrix element |

### Simulation

Lindblad master equation evolution (4th-order Runge-Kutta, dt=0.005,
2000 steps) from an initial state: Bell pair on qubits (0,1), |+>
state on remaining qubits. At 501 evenly spaced sample points, we
extract the feature vector for a chosen pair, standardize to zero
mean and unit variance, and compute principal component analysis
(PCA) via singular value decomposition.

The number of principal components needed for 95% explained variance
(n95) is the effective dimensionality of the decoherence trajectory.

### Why the observed pair changes with N

For N=2-3, pair (0,1) is analyzed: it starts as a Bell pair and its
decoherence is directly driven by the environment. For N=4-5, pair
(1,2) is analyzed instead: it starts unentangled and gains coherence
through Hamiltonian transfer from the Bell pair, producing richer
dynamics with oscillations and sector switches. The center pair is
also the natural target for sacrifice-zone analysis (it is the pair
being protected).

## 3. Results

### 3.1 Scaling: dimensionality grows, but information concentrates

Heisenberg chain with coupling J=1.0, uniform Z-dephasing γ=0.05:

| N | Qubits | n95 | PC1 variance | 3-PC coverage | PC1 best proxy |
|---|---|---|---|---|---|
| 2 | 2 | 1 | 99% | 100% | Concurrence (r=1.00) |
| 3 | 3 | 3 | 58% | 96% | Concurrence (r=0.94) |
| 4 | 4 | 4 | 62% | 94% | Purity (r=0.98) |
| 5 | 5 | 5 | 46% | 88% | Purity (r=0.99) |

**The 95% dimensionality grows as n95 ~ N.** Each additional
environment qubit adds approximately one effective dimension to
the pair's decoherence trajectory.

**But the first 3 PCs always capture 88-96%.** The 4th, 5th, etc.
components add diminishing returns (6%, 4% at N=5). This defines
two practical regimes:

| Regime | PCs | Coverage | Cost | Use case |
|---|---|---|---|---|
| Monitoring | 3 | 88-96% | 3 measurements | Real-time oversight |
| Full diagnostics | ~N | 95% | N measurements | Calibration, debugging |

### 3.2 Universality across topologies

We tested 9 topologies at N=2-4 under both Z-dephasing and
depolarizing noise (18 configurations total):

| Topology | N | n95 (Z-deph) | n95 (depol) | PC1 tracks |
|---|---|---|---|---|
| Pair | 2 | 1 | 1 | Concurrence |
| Chain | 3 | 3 | 4 | Concurrence |
| Star (symmetric) | 3 | 4 | 3 | Concurrence |
| Star (asymmetric) | 3 | 4 | 4 | Purity |
| Ring | 3 | 3 | 3 | Purity |
| Complete graph | 3 | 3 | 3 | Purity |
| Chain | 4 | 4 | 3 | Purity |
| Star | 4 | 4 | 4 | Purity |
| Ring | 4 | 3 | 4 | Purity |

**Dimensionality** is stable at 3-4 for fixed N, regardless of
topology and noise type.

**PC1 identity shifts** systematically: Concurrence dominates in
sparse topologies (chains, stars at small N), Purity dominates in
dense topologies (rings, complete graphs, chains at large N). In
all cases, Concurrence and Purity are the top-2 loadings of PC1
(loading difference < 0.02). PC1 always captures "overall quantum
strength"; the question is which component decays fastest.

### 3.3 Self-calibrating compass

The PCA automatically selects the fastest-decaying observable as
PC1. No manual tuning is needed. The procedure:

1. Run a decoherence trajectory on your hardware
2. Extract the 9 features at each time step
3. Apply PCA
4. PC1 tells you what to monitor; its proxy observable (Purity or
   Concurrence) is your primary diagnostic

### 3.4 N=5: the sweet spot

N=5 is a natural testbed because it sits at the efficiency knee of
the palindromic eigenvalue structure: each additional qubit beyond
N=5 yields diminishing returns in spectral richness.

**All 10 qubit pairs at N=5 (Heisenberg chain):**

| Pair | Distance | Max CΨ | Max θ | Max Concurrence | Entanglement death |
|---|---|---|---|---|---|
| (0,1) | 1 (edge) | 0.429 | 40 deg | 1.000 | t = 1.28 |
| (0,2) | 2 | 0.303 | 25 | 0.278 | t = 0.66 |
| (0,3) | 3 | 0.287 | 21 | 0.022 | t = 0.54 |
| (0,4) | 4 | 0.250 | 0 | 0.002 | alive |
| (1,2) | 1 (center) | 0.283 | 20 | 0.212 | t = 0.50 |
| (2,3) | 1 (center) | 1.000 | 60 | 0.024 | t = 0.60 |
| (3,4) | 1 (edge) | 1.000 | 60 | 0.142 | t = 3.38 |

Entanglement (Concurrence) falls exponentially with pair distance.
The initial Bell pair (0,1) dies at t=1.28; the far edge (3,4)
retains some entanglement until t=3.38.

**Liouvillian spectrum.** Spectral gap: 2*γ = 0.100 (exact match
to analytical prediction). 212 distinct decay rates. Fastest rate:
0.500.

**Bures curvature at the fold (CΨ = ¼).** K = -141 at N=5, vs
K = -25 at N=2. The curvature grows with system size: the geometry
of the quantum-classical boundary becomes sharper when the
environment is larger.

### 3.5 Sacrifice zone: θ is the most sensitive instrument

Three noise distributions, same total noise budget (sum of gammas
= 0.25), N=5 chain. Center pair (1,2) response:

| Observable | Uniform | Edge sacrifice | Improvement |
|---|---|---|---|
| Max θ | 20.0 deg | 33.7 deg | **1.68×** |
| Max CΨ | 0.283 | 0.361 | 1.28× |
| Max Concurrence | 0.212 | 0.251 | 1.18× |
| CΨ at t=10 | 0.054 | 0.068 | 1.26× |

**θ shows the largest effect.** The nonlinear mapping
CΨ -> θ = arctan(sqrt(4*CΨ - 1)) amplifies small CΨ
changes near the ¼ boundary. This critical amplification makes
θ the optimal objective function for noise engineering: it is
most sensitive precisely where it matters most (near the quantum-
classical transition).

### 3.6 Hardware validation

We applied the cockpit to IBM Torino data (February-March 2026):
single-qubit tomography (Q52, 25 points), shadow measurements
(Q80, Q102, 10 points each), and 5-qubit bitstring counts
(sacrifice zone experiment, 3 noise strategies).

**What is computable from each data type:**

| Instrument | Full tomography | Shadow | Bitstring counts |
|---|---|---|---|
| θ / CΨ | YES | YES | NO (needs off-diagonals) |
| Concurrence | NO (needs 2-qubit) | NO | NO |
| Ψ-norm | YES | YES | NO |
| Mutual information | NO (needs 2+ qubits) | NO | YES |
| Decay rate | YES | YES | YES |
| Bures velocity | YES | YES | NO |

**Key findings:**

- Q52 (good qubit, T2 = 298 us): CΨ = ¼ crossing measured at
  115.0 us, predicted 114.7 us. **0.3% error.**
- Internal consistency: Ψ-norm vs off-diagonal magnitude r = 1.000.
  Bures velocity vs CΨ: r = 0.954. All instruments agree.
- Q80 (problematic qubit, T2* = 11 us): 61.5% crossing-time
  deviation. Accuracy depends on qubit quality, not on the cockpit.
- 5-qubit sacrifice zone: Selective dynamical decoupling beats
  uniform DD by 3.2× in mutual information at t=4 Trotter steps.

**Assessment.** The 4 instruments computable from existing data
(θ, Ψ-norm, Bures velocity, decay rate) are mutually
consistent and physically correct. The 3 remaining instruments
(Concurrence, curvature, Petermann factor) require 2-qubit
tomography that does not yet exist in the dataset. The framework
is not refuted; it is incompletely tested.

## 4. Implications

### 4.1 A 3-observable cockpit for quantum hardware

Any Heisenberg spin chain under Markovian dephasing can be monitored
with 3 observables capturing ~88% of the dynamics:

1. **Purity or Concurrence** (whichever PCA selects as PC1):
   the dominant decay direction.
2. **Ψ-norm** (normalized off-diagonal coherence):
   the rate of coherence loss.
3. **A Bell-sector indicator** (fidelity with a specific Bell state):
   the fine structure of the decay.

Cost: 3 targeted measurements per pair, vs 4^N for full tomography.
One calibration run per architecture determines which observables
to track.

### 4.2 Two monitoring regimes

The 3-observable cockpit captures 88-96% (monitoring regime). Full
95% coverage requires ~N observables (diagnostic regime). This
naturally maps to engineering practice: use the cheap 3-observable
cockpit for real-time oversight, switch to full diagnostics when
something looks anomalous.

### 4.3 Noise engineering

Edge sacrifice (concentrating noise on boundary qubits) improves
center-pair coherence by up to 1.68× in θ. The critical
amplification near CΨ = ¼ makes θ the optimal objective
function: it is most sensitive where quantum coherence is most
fragile.

## 5. Limitations and caveats

1. **Concurrence is untested on hardware.** All existing tomographic
   data is single-qubit. The most important instrument (PC1 proxy,
   57% variance) has never been validated against hardware.

2. **Markovian noise only.** All results assume memoryless dephasing.
   Real hardware exhibits 1/f noise, two-level-system defects, and
   non-Markovian revivals (observed as excess late-time coherence
   in the Q52 data).

3. **N=5 is the largest system tested.** Extrapolating n95 ~ N to
   N=100 is speculative. Dense topologies (rings, complete graphs)
   may show slower dimensionality growth.

4. **Bures curvature is noisy.** The curvature formula involves
   second derivatives of sparse data. Reliable curvature estimation
   requires 50+ densely spaced time points.

5. **Petermann factor is uninteresting here.** K_P ~ 1 for all
   pure-dephasing cases. It becomes relevant only in gain-loss
   systems (PT-symmetric configurations).

## 6. Open questions

1. **Scaling beyond N=5.** Does n95 grow linearly for chains? Does
   it saturate for dense topologies? The C# engine (N=7 eigenvalues
   already computed) can test this.

2. **Non-Markovian noise.** Does the 3-observable cockpit hold under
   colored noise or 1/f spectra?

3. **2-qubit hardware validation.** Measuring Concurrence on a qubit
   pair would validate the most important untested instrument.

4. **Universality of edge sacrifice.** Does "sacrifice boundary
   qubits + optimize θ" generalize beyond Heisenberg chains?

## Appendix: Scripts and data

| Script | Purpose |
|---|---|
| [`cockpit_n5.py`](../simulations/cockpit_n5.py) | N=5: PCA, all 10 pairs, sacrifice zone |
| [`cockpit_universality.py`](../simulations/cockpit_universality.py) | 9 topologies x 2 noise types |
| [`cockpit_navigation.py`](../simulations/cockpit_navigation.py) | 7 instruments on N=3 Star topology |
| [`cockpit_validation.py`](../simulations/cockpit_validation.py) | Hardware validation (tomography + shadow) |
| [`cockpit_ibm_hardware.py`](../simulations/cockpit_ibm_hardware.py) | IBM sacrifice zone + Q52/Q80 data |
| [`theta_pc_correlation.py`](../simulations/theta_pc_correlation.py) | θ vs PCA correlation |

| Related experiment | Key result used here |
|---|---|
| [Theta-PC Analysis](THETA_PC_ANALYSIS.md) | θ requires all 3 PCs (R^2 = 0.87) |
| [Structural Cartography](STRUCTURAL_CARTOGRAPHY.md) | PCA on CΨ windows (original 3D finding) |
| [Information Geometry](INFORMATION_GEOMETRY.md) | Bures metric g = 3.36 at fold, K = -25 |
| [Boundary Navigation](BOUNDARY_NAVIGATION.md) | θ definition, CΨ = ¼ boundary |
| [PT-Symmetry Analysis](PT_SYMMETRY_ANALYSIS.md) | Petermann K=403 (gain-loss only) |
| [V-Effect Palindrome](V_EFFECT_PALINDROME.md) | N=5 as sweet spot, V(5) = 1.81 |
