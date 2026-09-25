# Decoherence Cockpit: A Low-Dimensional Diagnostic Framework for Open Quantum Systems

**Date:** April 2, 2026
**Status:** Complete
**Scripts:** `simulations/cockpit_*.py`, [`simulations/theta_pc_correlation.py`](../simulations/theta_pc_correlation.py)

## What this document is about

This document records a PCA of selected simulated feature dashboards.
The N=3-5 rows span 88-96% variance in the first three principal components;
the N=2 row is 100%. The broader table contains nine topology-size
configurations under two noise types. Purity, Concurrence, and coherence
magnitude are computed features in those simulations, not three established
hardware settings; this PCA does not establish measurement cost or replace
tomography. The PC1 proxy (Concurrence) is not yet hardware-validated. The
separate IBM Torino crossing-time comparison is reported in its own scope.

---

## Abstract

For the tested Heisenberg spin-chain simulations under local dephasing, the
N=3-5 rows span 88-96% of the variance in the first three principal components
of a selected feature dashboard; the N=2 row is 100%. The broader comparison
contains nine topology-size configurations under two noise types. The dashboard
includes Purity (trace of ρ squared), Concurrence (entanglement), and normalized
L₁ coherence (off-diagonal magnitude). Which computed feature best correlates
with PC1 depends on the architecture. This is a dimensionality description of
simulated features, not a hardware measurement protocol; measurement cost
remains unestablished.

## 1. Background

### The problem

A pair of qubits in a quantum computer is described by a 4x4 density
matrix ρ with 15 independent real parameters. Under the number-conserving
Heisenberg Hamiltonian and local Z dephasing used here,
the full state does not generally approach the globally maximally mixed state: the conserved
excitation-sector weights for the N=5 initial state are
`[1/16, 3/16, 1/4, 1/4, 3/16, 1/16]`. Complete mixing within those fixed
sectors would give the displayed pair limit
`diag(11/40, 9/40, 9/40, 11/40)`, with purity `101/400`, rather than `I/4`.
The simulated reduced state therefore follows a model- and initial-state-
dependent trajectory. Reconstructing all 15 parameters requires multiple measurement
settings per pair, repeated many times for statistics. For N qubits, full
reconstruction has exponential cost, with the exact setting count depending
on the tomography protocol. The PCA below does not convert its computed
features into a smaller certified setting count.

### Key definitions

**Purity** = Tr(ρ^2). Ranges from 1 (pure state) to 1/d (maximally
mixed). Measures how much the state has decohered.

**Concurrence** (Wootters). Ranges from 0 (no entanglement) to 1
(maximally entangled). Measures how much entanglement remains between
two qubits.

**Ψ-norm** = L₁/(d-1), where L₁ is the sum of absolute values of
all off-diagonal elements of ρ, and d is the Hilbert space dimension.
Measures how much quantum coherence (superposition) remains.

**CΨ** = Purity x Ψ-norm. It is the selected product coordinate used here.
For the separate scalar recursion `R = C(Ψ+R)^2`, `CΨ = ¼` is its algebraic
discriminant-zero point. That coordinate is not a physical quantum/classical boundary
for the simulated density-matrix dynamics, and the recursion does
not supply attractors of the Lindblad evolution.

**θ** = arctan(sqrt(4*CΨ - 1)). Defined only when CΨ > ¼.
This is a nonlinear remapping of the algebraic discriminant coordinate,
not an angular distance to a physical quantum/classical boundary. It is
zero at `CΨ = ¼` and 60 deg for a pure coherent product state.

**Sacrifice zone.** In a qubit chain, deliberately concentrating noise
on boundary qubits ("sacrificing" them) to protect center qubits.
The total noise budget stays the same; only its distribution changes.

**Liouvillian.** The generator of open quantum dynamics. For N qubits
under dephasing, it is a 4^N x 4^N matrix whose eigenvalues give
the decay rates and whose eigenvectors give the decay modes.

**Bures distance.** A metric on the space of density matrices that
quantifies how distinguishable two quantum states are. The N=5 producer
uses adjacent Bures distances to estimate a one-dimensional Bures
path-metric coefficient `g_path(CΨ)`. It then prints a coordinate-shape
second derivative `S_CΨ = -(2g_path)^-1 d²log(g_path)/dCΨ²`. This finite
path-coordinate proxy is branch-, coordinate-, and stencil-dependent; no
two-dimensional metric or intrinsic/Gaussian curvature is computed.

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
(PCA, the standard technique that finds the directions of maximum variance in high-dimensional data) via singular value decomposition.

The number of principal components needed for 95% explained variance
(n95) is the effective dimensionality of the decoherence trajectory.

### Why the observed pair changes with N

For N=2-4, pair (0,1) is analyzed: it starts as a Bell pair. The finite scaling
table uses pair (1,2) only at N=5; that pair starts unentangled and
receives coherence through Hamiltonian evolution from the Bell pair. That
changed focus pair confounds a direct size-only comparison. Pair (1,2) is also
the target in the separate finite sacrifice-profile comparison.

## 3. Results

### 3.1 Scaling: four finite rows with a changed N=5 focus

Heisenberg chain with coupling J=1.0, uniform Z-dephasing γ=0.05:

| N | Qubits | n95 | PC1 variance | 3-PC coverage | PC1 best proxy |
|---|---|---|---|---|---|
| 2 | 2 | 1 | 99% | 100% | Concurrence (r=1.00) |
| 3 | 3 | 3 | 58% | 96% | Concurrence (r=0.94) |
| 4 | 4 | 4 | 62% | 94% | Purity (r=0.98) |
| 5 | 5 | 5 | 46% | 88% | Purity (r=0.99) |

These four displayed rows do not establish an n95 scaling law: the focus is
pair (0,1) at N=2-4 but a different focus pair at N=5. The N=3-5 rows span
88-96% in the first three PCs; the N=2 row is 100%. At N=5 the fourth and
fifth components add 6% and 4% of this selected-dashboard variance.

| Regime | PCs | Coverage | Cost | Use case |
|---|---|---|---|---|
| Three-PC summary | 3 | 88-96% | Not established | Simulated dashboard |
| 95% PCA summary | 1, 3, 4, 5 in the four rows | 95% | Not established | Simulated dashboard |

### 3.2 Finite topology-size configurations

We tested nine topology-size configurations at N=2-4 under both Z-dephasing
and depolarizing noise (18 noise/configuration rows total):

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

The best-correlated feature label changes across these rows. This correlation
label is not a stable physical identity and PCA is not a decay-rate ranking.
The finite table records covariance in standardized simulated features.

### 3.3 PCA directions and candidate feature proxies

PCA finds covariance directions in the nine standardized simulated features.
The subsequent correlation table names a candidate feature proxy for each PC;
it does not select the fastest-decaying observable and does not identify a
hardware monitor. Several input features themselves require reduced-state
reconstruction, so the table does not define a low-cost measurement protocol.

### 3.4 N=5: the largest row in this local pair scan

N=5 is the largest system in the all-pairs table below. It is a useful
testbed, but this finite cockpit calculation does not choose one preferred
system size or license a trend outside the sampled rows.

**All 10 qubit pairs at N=5 (Heisenberg chain):**

| Pair | Distance | Max CΨ | Max θ | Max Concurrence | C > 0.01 threshold status |
|---|---|---|---|---|---|
| (0,1) | 1 (edge) | 0.429 | 40 deg | 1.000 | t = 1.28 |
| (0,2) | 2 | 0.303 | 25 | 0.278 | t = 0.66 |
| (0,3) | 3 | 0.287 | 21 | 0.022 | t = 0.54 |
| (0,4) | 4 | 0.250 | 0 | 0.002 | never above 0.01 |
| (1,2) | 1 (center) | 0.283 | 20 | 0.212 | t = 0.50 |
| (1,3) | 2 | 0.260 | 11 | 0.000 | never above 0.01 |
| (1,4) | 3 | 0.277 | 18 | 0.000 | never above 0.01 |
| (2,3) | 1 (center) | 1.000 | 60 | 0.024 | t = 0.60 |
| (2,4) | 2 | 1.000 | 60 | 0.040 | t = 1.66 |
| (3,4) | 1 (edge) | 1.000 | 60 | 0.142 | t = 3.38 |

The ten displayed maxima are position-dependent and nonmonotone; no
exponential distance fit was performed. “Never above 0.01” is distinct from
remaining alive above the threshold at the end of the sampled trajectory.

**Liouvillian spectrum.** Spectral gap: 2*γ = 0.100 (exact match
to analytical prediction). 212 distinct decay rates. Fastest rate:
0.500, matching the independent uniform-dephasing control
`2*N*gamma = 0.500` for `X^⊗N`.

**Bures path-metric coefficient and coordinate-shape second derivative near
the sampled CΨ = ¼ point.** The two displayed `S_CΨ` readings are finite
path-coordinate proxies: -141 at N=5 and -25 at N=2, from sparse,
nonmonotone sampled trajectories, not Gaussian or intrinsic curvature. They
do not establish growth with system size, sharper boundary geometry, or
state-space divergence.

### 3.5 Finite N=5 sensitivity comparison

Three noise distributions, same total noise budget (sum of gammas
= 0.25), N=5 chain. Center pair (1,2) response:

| Observable | Uniform | Edge sacrifice | Improvement |
|---|---|---|---|
| Max θ | 20.0 deg | 33.7 deg | **1.68×** |
| Max CΨ | 0.283 | 0.361 | 1.28× |
| Max Concurrence | 0.212 | 0.251 | 1.18× |
| CΨ at t=10 | 0.054 | 0.068 | 1.26× |

**θ has the largest displayed relative change among these finite N=5 rows.**
The nonlinear mapping
CΨ -> θ = arctan(sqrt(4*CΨ - 1)) amplifies small CΨ
changes near the ¼ boundary. That remapping explains the displayed
sensitivity; this finite comparison does not establish θ as an optimal
hardware objective or identify where noise engineering matters most.

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

- Q52 (good qubit, calibration-era run): Q52 is a qualitative crossing record, not a precision match: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the generalized prediction 0.936. The legacy 115.0-versus-114.7 comparison
  recomputed the same hardware record; it was not an independent prediction.
- Internal consistency: Ψ-norm vs off-diagonal magnitude r = 1.000, which
  is the identity Ψ = 2|ρ01| at d = 2 rather than a check; Bures velocity
  vs CΨ r = 0.954, which could have come out otherwise, though an
  exponential decay makes a high correlation close to forced.
- Q80: 61.5% crossing-time deviation against a T2* = 11 us measured six
  days earlier, and 1.9% against a same-day Ramsey T2* = 17.36 us
  ([IBM Run 3](IBM_RUN3_PALINDROME.md)). Accuracy depends on how fresh
  the input calibration is, not on the cockpit.
- 5-qubit sacrifice zone: Selective dynamical decoupling beats
  uniform DD by 3.2× in mutual information at t=4 Trotter steps.

**Assessment.** The 4 instruments computable from existing data
(θ, Ψ-norm, Bures velocity, decay rate) give physically correct values.
None of the agreements is a strong check. θ = arctan(√(4CΨ−1)) and
Ψ = 2|ρ01| are functions of quantities already read, so their consistency
is arithmetic; and under dephasing the Bures velocity goes as |ρ01|, which
makes its high correlation with CΨ close to forced. What that correlation
tests is the exponential decay model, not the instruments'
independence. The 3 remaining instruments
 (Concurrence, the path-coordinate proxy, Petermann factor) require 2-qubit
tomography that does not yet exist in the dataset. The framework
is not refuted; it is incompletely tested.

## 4. Implications

### 4.1 A three-PC summary of selected simulated features

Across N=3-5 under Markovian dephasing, the first three PCs capture 88-96%
of the selected dashboard variance; the N=2 row is 100%. The PCs are linear
combinations of standardized features. The named Purity, Concurrence, and
Ψ-norm links below are candidate correlations, not selected hardware
observables:

1. **Purity or Concurrence** (whichever correlates most with PC1):
   a candidate feature proxy for that covariance direction.
2. **Ψ-norm** (normalized off-diagonal coherence):
   a coherence-magnitude feature proxy, not a decay-rate estimate.
3. **A Bell-sector indicator** (fidelity with a specific Bell state):
   the fine structure of the decay.

These are computed features, and several require reduced-state reconstruction.
The simulations do not establish a three-setting hardware protocol or its
measurement cost.

### 4.2 Dashboard summaries, not monitoring regimes

The three-PC and 95%-variance summaries are not two established monitoring
regimes. Whether either PCA summary yields a cheaper hardware-monitoring
protocol remains open.

### 4.3 Noise engineering

In the tested N=5 profiles, the edge-sacrifice row has a 1.68× larger computed
maximum θ than the uniform row and a larger maximum center-pair CΨ. The profiles
change several site rates together, so this finite comparison neither isolates
a cause nor establishes θ as an optimal hardware objective.

## 5. Limitations and caveats

1. **Concurrence is untested on hardware.** All existing tomographic
   data is single-qubit. Its correlation with a simulated PC score has not
   been validated as a hardware diagnostic.

2. **Markovian simulation only.** The simulations here use memoryless
   dephasing. Detuning is the preferred explanation for the phase component,
   but the Q52 late-time excess mechanism remains unresolved absent a
   Q52-specific fit/control. Only the universal-boundary/non-Markovian-witness
   interpretation is closed; colored-noise and memory-kernel tests remain open.

3. **N=5 is the largest system tested.** The four rows do not establish
   n95 ~ N because the focus changes at N=5. Scaling beyond this finite,
   confounded comparison remains open.

4. **The coordinate-shape second derivative is numerically fragile.** It
   applies two derivatives to a sparse, nonmonotone path coordinate. Dense
   sampling plus explicit monotone-branch handling would be needed even to
   stabilize this coordinate-dependent proxy; it is not intrinsic curvature.

5. **The old blanket Petermann null is refuted.** The early cockpit sample found K_P near 1,
   but pure Z-dephasing Liouvillians can be strongly non-normal: a simple N=5 mode has
   `||P||=sqrt(375)`, and real-axis defective seeds are certified at N=5,7,9. Petermann
   readings are therefore relevant beyond gain-loss/PT systems. A single-eigenvector value
   is not basis-invariant inside a degenerate eigenspace; use it only for a simple isolated mode,
   and use subspace/Jordan diagnostics at degeneracy.

## 6. Open questions

1. **Scaling beyond N=5.** Does n95 grow linearly for chains? Does
   it saturate for dense topologies? The C# engine (N=7 eigenvalues
   already computed) can test this.

2. **Non-Markovian noise.** Does the three-PC variance summary persist under
   colored noise or 1/f spectra in the same simulated feature dashboard?

3. **2-qubit hardware validation.** Measuring Concurrence on a qubit
   pair would validate the most important untested instrument.

4. **Edge-profile scope.** Does the finite edge-profile difference persist
   beyond Heisenberg chains, and which objective would be appropriate?

## Appendix: Scripts and data

| Script | Purpose |
|---|---|
| [`cockpit_n5.py`](../simulations/cockpit_n5.py) | N=5: PCA, all-pair threshold-status inventory, sacrifice profiles |
| [`cockpit_universality.py`](../simulations/cockpit_universality.py) | Nine topology-size configurations x 2 noise types |
| [`cockpit_navigation.py`](../simulations/cockpit_navigation.py) | 7 instruments on N=3 Star topology |
| [`cockpit_validation.py`](../simulations/cockpit_validation.py) | Hardware validation (tomography + shadow) |
| [`cockpit_ibm_hardware.py`](../simulations/cockpit_ibm_hardware.py) | IBM sacrifice zone + Q52/Q80 data |
| [`theta_pc_correlation.py`](../simulations/theta_pc_correlation.py) | θ vs PCA correlation |

| Related experiment | Key result used here |
|---|---|
| [Cockpit Scaling](COCKPIT_SCALING.md) | Finite selected-feature PCA to N=7-11; no size law established |
| [Theta-PC Analysis](THETA_PC_ANALYSIS.md) | θ requires all 3 PCs (R^2 = 0.87) |
| [Structural Cartography](STRUCTURAL_CARTOGRAPHY.md) | PCA on CΨ windows (original 3D finding) |
| [Information Geometry](INFORMATION_GEOMETRY.md) | Bures path-metric coefficient and its coordinate-shape second derivative; the Bures geodesic test in the full state space |
| [Boundary Navigation](BOUNDARY_NAVIGATION.md) | θ definition, CΨ = ¼ boundary |
| [PT-Symmetry Analysis](PT_SYMMETRY_ANALYSIS.md) | Gain-loss example; not the only setting with strong non-normality |
| [V-Effect Palindrome](V_EFFECT_PALINDROME.md) | Finite V-Effect census; F6 Q-edge gain is a separate ratio |
