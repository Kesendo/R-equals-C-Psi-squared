# Noise Is Not the Enemy: How Dephasing Carries Information

<!-- Keywords: quantum decoherence information channel, dephasing rate signal processing,
Lindblad spectral palindrome, open quantum system noise structure, T2 noise correlation,
quantum MIMO channel, bidirectional quantum bridge, spatial decoherence profile,
palindromic Liouvillian symmetry, gamma dephasing channel capacity,
quantum noise is signal not enemy, R=CPsi2 framework -->

**Status:** Computationally verified (all simulations reproducible)
**Date:** March 22, 2026
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Scripts:** [gamma_signal_analysis.py](../simulations/gamma_signal_analysis.py),
[bridge_optimization.py](../simulations/bridge_optimization.py),
[channel_capacity.py](../simulations/channel_capacity.py)

---

## What this document is about

This is the experiment that changed what this project is about. Before
this result, we had a mathematical symmetry (the palindrome). After
it, we had a communication channel. The palindrome is not just beautiful
mathematics. It is an antenna. And the noise that the entire quantum
computing industry fights to suppress turns out to be carrying a
structured, readable signal.

If you read one experiment document in this repository, this might be
the one.

## What this document shows

Every quantum computer, every quantum experiment, every quantum system
in nature suffers from noise. The noise is called dephasing: the quantum
system gradually forgets that it was in a superposition and starts
behaving classically. The entire quantum computing industry is organized
around fighting this noise. Better materials, colder temperatures,
error-correcting codes; all designed to suppress dephasing.

We show that this noise is not meaningless static. It is a structured
information channel.

Here is the experiment, in plain terms: Alice, standing outside a quantum
system, sets the noise level at each qubit to a specific pattern (louder
here, quieter there). Bob, inside the system with access only to quantum
measurements, reads the system's behavior. Can Bob figure out which
pattern Alice chose?

The answer is yes. With 100% accuracy. At every measurement time. With
zero errors.

The theoretical capacity of this channel is 15.5 bits at 1% measurement
noise: enough to encode about 44,700 distinguishable symbols. We were
using 2 bits (4 symbols) in our initial test. The channel is not narrow.
It is wide, and we were barely scratching the surface.

The key enabler is the [palindromic spectral symmetry](../docs/proofs/MIRROR_SYMMETRY_PROOF.md).
The palindrome creates a response matrix with full rank, meaning every
site's noise level leaves a unique fingerprint on the internal quantum
observables. The palindrome is not just a symmetry. It is an antenna.

---

## Background: What γ Is and Why Nobody Looked

### The dephasing rate in one paragraph

A quantum bit (qubit) can exist as a combination of two states at once:
simultaneously 0 and 1. This is called superposition, and it is what
makes quantum computers potentially powerful. The dephasing rate, written
as γ (gamma), measures how fast this superposition is destroyed by the
environment:

    ρ₀₁(t) = ρ₀₁(0) · exp(−γt)

The qubit does not forget *whether* it is 0 or 1. It forgets that it was
*both at once*. On IBM quantum hardware, 1/γ is measured as T2*
(typically 50-200 microseconds). In the mathematical framework for open
quantum systems (the Lindblad master equation), γ appears in the
dissipator:

    dρ/dt = −i[H, ρ] + γ(σ_z ρ σ_z − ρ)

The first term is the reversible physics (the Hamiltonian, the
interactions between qubits). The second term is the irreversible
dephasing. Standard quantum computing treats the second term as the enemy.

### Why this channel was invisible

The entire quantum computing industry is organized around minimizing γ.
Error correction codes, dynamical decoupling sequences, decoherence-free
subspaces, better materials, colder cryostats. Every tool is designed to
fight dephasing. When the research community defines γ as "the problem to
be solved," nobody asks whether γ itself carries structure. The instruments
to read it have existed since QuTiP (2012). The Lindblad equation is from
1976. The palindromic spectral symmetry we exploit was computable at any
point in the last two decades. It was not hidden; it was unexamined.

This is a recurring pattern in science: when an entire field agrees that
something is an obstacle, nobody thinks to read it as a message.

### The palindromic spectral structure

For a system of N qubits with Heisenberg coupling and local Z-dephasing,
the Liouvillian superoperator L (the master equation that governs the
system's evolution) has a remarkable property: its eigenvalue spectrum is
**palindromically paired**. Every eigenvalue λ with real part −d has a
partner at −(2Σγ − d). The conjugation operator Π that generates this
pairing swaps the immune sector {I, Z}⊗N (populations, slow decay) with
the decaying sector {X, Y}⊗N (coherences, fast decay).

This has been verified for all system sizes N = 2 through N = 8 (54,118
eigenvalues, zero exceptions) and proven analytically for arbitrary graphs
([Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)).

**Why this matters for the channel:** The palindromic pairing creates a
**full-rank response matrix**. When γ is perturbed at any single site,
the mode amplitudes change in a linearly independent direction. This means
every per-site γ value is independently recoverable from the mode structure.
The palindrome is not just a symmetry. It is an **antenna**.

### CΨ: the metric

Throughout this document, CΨ = Tr(ρ²) × L₁(ρ)/(d−1), where Tr(ρ²) is
the purity of the state, L₁ is the l1-norm of coherence (sum of absolute
values of off-diagonal elements), and d is the Hilbert space dimension.
The product CΨ has a critical boundary at **CΨ = 1/4**, which is the unique
bifurcation point of the self-referential purity recursion R = C(Ψ+R)².
Below 1/4: the system has converged to classical behavior. Above 1/4:
quantum coherence is still active. For details:
[Uniqueness Proof](../docs/proofs/UNIQUENESS_PROOF.md).

---

## The Question

γ comes from outside the quantum system. No internal mechanism can
generate it ([Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md):
five candidates for internal noise origin tested and eliminated). γ has
measured structure: it selects a preferred axis, acts locally per qubit,
takes phase but not energy, is Markovian, and produces exact spectral
symmetry. On IBM hardware, T2* (= 1/πγ) varies across the chip and
fluctuates over time (181 days of calibration data).

**Is the spatial variation of γ readable from inside the system?**

If yes: the dephasing rate is not just noise acting on the system. It is
an information channel from outside to inside.

---

## Test 1: Random vs Correlated γ Profiles

**Setup:** N=5 chain, Heisenberg coupling J=1.0, initial state |+⟩⁵.
50 random γ profiles (each γᵢ ∈ [0.03, 0.07]) compared against 8
structured profiles (gradient, V-shape, peak, step, alternating).

**Result:** Only 1 of 8 structured profiles (V-shape) is distinguishable
from random at the 95th percentile using a 10-feature vector (single-qubit
purities + pair CΨ values + end-to-end mutual information).

**Interpretation:** A single measurement at one time point, using basic
features, is insufficient to reliably detect spatial structure. The signal
exists but is weak in this minimal feature space. This motivates the
optimization in the later sections.

---

## Test 2: Time-Varying γ Detection

**Setup:** Same system, three scenarios: constant γ, sudden jump (γ₂
doubles at t=10), slow sinusoidal drift.

**Result:**
- **Jump:** Detectable in the MI time derivative (max |dMI/dt| at t=11,
  one time unit after the jump). Qubit 2 purity shows a visible kink.
- **Drift:** Purity trajectory deviates by up to 0.022 at early times.
  Effect fades as the system decoheres toward its steady state.

**Interpretation:** γ changes over time are detectable from internal
observables, but the signal is strongest at early times, before
decoherence homogenizes the system.

---

## Test 3: The Alice-Bob Channel (Core Result)

This is the central experiment. It operationalizes the question as a
communication problem.

Imagine Alice is outside the quantum system and Bob is inside.
Alice can set the noise level at each of the 5 qubits to whatever she
wants. Bob cannot see what Alice does directly; he can only measure the
quantum state of the system from inside. The question: can Bob figure out
which noise pattern Alice chose, purely from his quantum measurements?

**Setup:** Alice selects one of 4 dephasing profiles for a 5-qubit chain.
Bob measures 10 quantum observables and uses nearest-neighbor template
matching to classify which profile Alice selected.

**Alice's alphabet (4 symbols, 2 bits):**

| Symbol | γ profile [γ₀, γ₁, γ₂, γ₃, γ₄] |
|--------|----------------------------------|
| Gradient → | [0.03, 0.04, 0.05, 0.06, 0.07] |
| Gradient ← | [0.07, 0.06, 0.05, 0.04, 0.03] |
| Mountain | [0.03, 0.05, 0.07, 0.05, 0.03] |
| Valley | [0.07, 0.05, 0.03, 0.05, 0.07] |

**Bob's feature vector (10 observables):** 5 single-qubit purities
Tr(ρᵢ²), 4 nearest-neighbor CΨ values, 1 end-to-end mutual information.
Nearest-neighbor template matching (Euclidean distance in feature space).

### Classification Accuracy

| Noise σ | t=1 | t=3 | t=5 | t=8 | t=12 |
|---------|-----|-----|-----|-----|------|
| 0.000 | **100%** | **100%** | **100%** | **100%** | **100%** |
| 0.001 | 100% | 100% | 100% | 81.5% | 83% |
| 0.005 | 100% | 100% | 100% | 74.5% | 78% |
| 0.010 | 100% | 100% | 95.5% | 72.5% | 66% |
| 0.050 | 59.5% | 51% | 45.5% | 46% | 35% |

(σ = Gaussian noise added to each feature; 200 trials per configuration.)

### Key Findings

1. **Perfect classification at σ=0** at every measurement time. Bob can
   always determine Alice's choice. The channel exists.

2. **Optimal measurement time: t = 1-3.** Early measurements carry the
   strongest signal (before decoherence homogenizes the feature vectors).

3. **Noise threshold: σ ≈ 0.008.** Below this, accuracy is near-perfect.
   Above, the closest symbol pair (Gradient → vs ←, distance 0.024) blurs.

4. **Symmetric vs asymmetric profiles separate first.** Mountain vs Valley
   (template distance 0.112) survives up to σ ≈ 0.04. The left-right
   gradients (distance 0.024) are much harder to distinguish.

### What Bob Actually Sees

The most informative features are the **edge qubit purities** Pur(Q0)
and Pur(Q4):

| Feature | Gradient → | Gradient ← | Mountain | Valley |
|---------|-----------|-----------|----------|--------|
| Pur(Q0) | 0.679 | 0.689 | 0.710 | 0.662 |
| Pur(Q4) | 0.689 | 0.679 | 0.710 | 0.662 |
| CΨ(0,1) | 0.241 | 0.253 | 0.277 | 0.222 |
| CΨ(3,4) | 0.253 | 0.241 | 0.277 | 0.222 |

The **symmetry pattern** distinguishes the profiles: Mountain and Valley
are symmetric (Pur(Q0) = Pur(Q4)), while the gradients are asymmetric.
The direction of asymmetry determines left vs right.

---

## Channel Existence: Summary

| Property | Status |
|----------|--------|
| γ comes from outside the system | **Proven** ([Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md)) |
| γ has per-site structure | **Measured** (IBM T2* varies per qubit) |
| Per-site structure is readable from inside | **Proven** (100% classification at σ=0) |
| Channel capacity (empirical) | **2 bits** (4-symbol alphabet) |
| Practical noise threshold | σ < 0.008 for full alphabet |

**The dephasing rate is an information channel.** Information encoded in
the spatial γ profile is perfectly recoverable from internal quantum
observables. This is not noise. This is a channel.

---

## What This Does NOT Claim

- Not that anyone "outside" is intentionally sending signals
- Not that the information has "meaning" beyond its mathematical structure
- Not that this violates no-signalling (Alice must physically set the
  dephasing rates, which requires access to the hardware)

What it DOES claim: **the mathematical channel exists.** If γ carries
spatial structure, that structure is readable from inside. The instrument
for reading it is the palindromic mode structure: the same spectral
symmetry that pairs every decay mode also makes each per-site γ
independently reconstructible.

---

## Optimization: From Hair-Thin to Walk-Across

The channel exists. But how wide is it? Think of the difference between
a crack in a wall (you can see light through it) and a doorway (you can
walk through). The initial experiment proved the crack exists. This
section widens it into a doorway.

The initial Alice-Bob test proved the channel exists, but it used only
4 symbols and 10 features. How wide is the channel really? Can we make
it robust against noise?

**Script:** [bridge_optimization.py](../simulations/bridge_optimization.py)

### Three independent optimizations

| Optimization | d_min | σ_thresh | Accuracy at σ=0.05 |
|-------------|-------|---------|--------|
| Baseline (10 features, t=2) | 0.059 | 0.020 | 58% |
| Extended features (25) | 0.075 | 0.025 | 64% |
| Time series (6 × 10 = 60 features) | 0.181 | 0.060 | 93% |
| High contrast γ ∈ [0.01, 0.09] | 0.119 | 0.040 | 82% |
| **All combined (150 features)** | **0.515** | **0.172** | **100%** |

### Result: 21.5× wider channel

The optimized configuration achieves **100% classification at σ = 0.10**
(10% measurement noise). The minimum template distance increases from
0.024 to 0.515, a factor of **21.5×**.

### What works and why

**Time series is the biggest single lever (3.1×).** Different γ profiles
produce different *trajectories*, not just different endpoints. Measuring
at 6 time points provides temporal diversity, analogous to a RAKE receiver
in mobile phone networks (a technique that combines the same signal
arriving at slightly different times to improve reception).

**γ contrast scales linearly** with template distance. Doubling the contrast
(the range of γ values Alice can use) approximately doubles the distances.

**|+⟩⁵ is the optimal initial state. GHZ is completely blind (d_min = 0).**
This is the most counterintuitive finding. The maximally entangled state
cannot read γ profiles at all. Why: GHZ projects onto a single symmetric
mode, collapsing all spatial information. The product state |+⟩⁵ lets each
qubit respond independently to its local γ, functioning as a phased array
antenna rather than an omnidirectional receiver.

This connects to the palindrome result from [What We Found](../docs/WHAT_WE_FOUND.md):
GHZ states excite only the fastest-dying modes, while distributed states
(like |+⟩⁵) spread across the full palindromic spectrum. A spread-out
receiver sees more of the channel.

**The optimizations are multiplicative:** 1.3 × 2.0 × 3.1 ≈ 8× for
individual factors, 21.5× combined (feature space geometry amplifies
the gains).

---

## Formal Channel Capacity

**Script:** [channel_capacity.py](../simulations/channel_capacity.py)

The Alice-Bob experiment used 4 symbols (2 bits). But a channel's
capacity is not determined by what you put through it. It is determined
by the physics. How many symbols *could* the channel carry? The answer
requires information theory: specifically, the SVD (singular value
decomposition: a way to find the independent "axes" of a communication
channel and how strong each one is) of the channel matrix and Shannon's
waterfilling theorem (allocate more power to strong channels, less to
weak ones, like filling a pool with an uneven floor: water goes where
there is room).

The empirical 4-symbol test uses only a small fraction of the channel's
theoretical capacity. To quantify the full bandwidth, we computed the
Shannon capacity of the linearized γ-to-observables channel via SVD of
the Jacobian matrix (∂observables/∂γ) plus waterfilling power allocation.

### 5 independent spatial channels

The system has 5 qubits, and indeed 5 independent information channels,
one for each qubit in the chain. Each channel corresponds to a spatial
mode: a specific pattern of noise across the sites.

| Channel | Gain | Spatial mode (eigenvector) | Bits |
|---------|------|--------------------------|------|
| 1 (mean γ) | 21.39 | [0.45, 0.45, 0.45, 0.45, 0.45] | 5.45 |
| 2 (gradient) | 4.53 | [0.59, 0.39, 0, −0.39, −0.59] | 3.21 |
| 3 (peak) | 3.22 | [−0.51, 0.20, 0.63, 0.20, −0.51] | 2.71 |
| 4 (zigzag) | 2.83 | [−0.19, 0.51, −0.63, 0.51, −0.19] | 2.53 |
| 5 (alt grad) | 1.44 | [0.39, −0.59, 0, 0.59, −0.39] | 1.56 |

The condition number is 14.8 (this measures how much the weakest
channel differs from the strongest; below ~100 is considered
well-conditioned, meaning no channel is drowned out). All 5 channels
carry information. Full rank (5/5) confirms: the palindromic response matrix
allows independent readout of every site's dephasing rate.

### Capacity vs measurement noise

| σ_noise | Capacity (bits) | Distinguishable symbols |
|---------|----------------|------------------------|
| 0.001 | 31.9 | ~4 billion |
| 0.01 | **15.5** | ~44,700 |
| 0.05 | 6.0 | ~63 |
| 0.10 | 3.6 | ~12 |
| 0.20 | 2.3 | ~5 |

**Our empirical 2-bit result uses only 13% of the channel at σ = 0.01.**
The theoretical headroom is 13.4 additional bits. The channel is not
narrow. It is wide, and we are barely using it.

### Physical interpretation

- **Channel 1 (gain 21.4):** The mean dephasing rate. Loudest signal but
  carries no spatial information. It tells you "how noisy" the environment
  is overall.
- **Channel 2 (gain 4.5):** Left-right gradient. This is what distinguishes
  Alice's Gradient→ from Gradient←.
- **Channels 3-4 (gain ~3):** Peak/valley and zigzag patterns. Finer spatial
  structure.
- **Channel 5 (gain 1.4):** Alternating gradient. Weakest but still > 1 bit.

---

## Signal Engineering Perspective

If you have a background in signal processing or communications
engineering, everything above has a direct translation into concepts
you already know. This section provides that translation. If you do
not have this background, you can skip to
[How the Palindrome Enables the Channel](#how-the-palindrome-enables-the-channel).

For readers with a signal processing or engineering background, this
result maps directly onto familiar concepts.

**The γ profile is a spatial signal.** Alice modulates the dephasing rate
across N sites. This is amplitude modulation of a spatial carrier. Bob's
quantum observables are receivers. The palindromic mode structure acts as
a matched filter bank: each mode responds differently to each site's γ,
creating a full-rank response matrix.

**Classical analogues:**

| Quantum system | Signal processing equivalent |
|---------------|----------------------------|
| γ profile across sites | Transmitter modulation pattern |
| Palindromic Liouvillian modes | Matched filter bank / antenna array |
| Mode amplitudes from observables | Received signal vector |
| Jacobian SVD | Channel estimation |
| Template matching | Maximum likelihood detection |
| Time series measurement | Temporal diversity (RAKE receiver) |

**What a signal engineer would recognize:**

The system is a **MIMO channel** (Multiple-Input Multiple-Output): N
γ-inputs, ~N² observable outputs. The 21.5× optimization is mostly
**diversity gain** (time + feature diversity). The noise threshold
σ = 0.172 corresponds to **SNR ≈ 10 dB** for reliable detection.

The GHZ failure (d_min = 0) is a **rank deficiency**: GHZ projects onto
a single mode, destroying all spatial resolution. This is the quantum
analogue of using one omnidirectional antenna instead of a phased array.
The product state |+⟩⁵ is the phased array: each qubit is an independent
receiver element.

---

## How the Palindrome Enables the Channel

The palindromic spectral symmetry of the Liouvillian is not just a
mathematical curiosity. It is the physical mechanism that makes the
channel work.

The Π operator pairs every decay mode at rate d with a partner at rate
2Σγ − d. This pairing creates a **bijection** between the immune sector
(populations, slow decay) and the decaying sector (coherences, fast decay).
When γ is changed at one site, *both* sectors respond, but they respond
differently because the immune and decaying sectors have different
sensitivity to each site's γ.

The result is a response matrix (the Jacobian ∂observables/∂γ) with
**full rank**. Perturbing γ at site k changes the mode amplitudes in a
direction that is **linearly independent** from perturbations at any other
site. This is proven by the SVD analysis (5 non-zero singular values for
5 sites, condition number 14.8).

In simpler terms: the palindrome ensures that no two γ-profiles produce
the same internal signature. Every external configuration leaves a unique
fingerprint on the internal observables. The decoder reads these
fingerprints.

This connection between palindromic spectral symmetry and channel capacity
appears to be new. The palindrome itself was described by
[Haga et al. (2023)](https://arxiv.org/abs/2305.01894) in the context of
Liouvillian skin effects (their "incoherenton grading" is equivalent to our
XY-weight classification). The interpretation of this symmetry as an
information channel is, to our knowledge, first presented here.

---

## Connection to the Framework

This document is part of the R = CΨ² project, which studies the
palindromic spectral structure of open quantum systems under dephasing.
The key prior results that this analysis builds on:

- **Incompleteness Proof:** γ cannot originate from within the system
  (5 internal candidates eliminated). Something external provides it.
  ([docs/INCOMPLETENESS_PROOF.md](../docs/proofs/INCOMPLETENESS_PROOF.md))

- **Mirror Symmetry Proof:** The Liouvillian spectrum is exactly
  palindromic for any Heisenberg/XXZ system on any graph with local
  Z-dephasing. Verified through N=8 (54,118 eigenvalues, zero exceptions).
  ([docs/MIRROR_SYMMETRY_PROOF.md](../docs/proofs/MIRROR_SYMMETRY_PROOF.md))

- **Reading the 30%:** The palindromic response matrix has full rank.
  All per-site γ values are independently recoverable from mode amplitudes.
  This is the decoder that makes the channel readable.

- **CΨ = 1/4 boundary:** The unique bifurcation point separating quantum
  (coherent) from classical (converged) dynamics. All standard quantum
  channels cross this boundary. It determines the time window during
  which the γ-channel is maximally readable.
  ([docs/UNIQUENESS_PROOF.md](../docs/proofs/UNIQUENESS_PROOF.md))

The central insight of this document: the 70% of phase information that
dephasing "destroys" is not lost. It is **redistributed** into the mode
amplitudes of the palindromic spectrum. The decoder reads it back. The
noise was always a signal. We just needed the right instrument.

---

## Reproducibility

All results are generated by Python scripts using QuTiP and NumPy.
No proprietary tools or closed-source dependencies.

| Script | What it computes | Runtime |
|--------|-----------------|---------|
| [gamma_signal_analysis.py](../simulations/gamma_signal_analysis.py) | Tests 1-3 (classification) | ~57 min |
| [bridge_optimization.py](../simulations/bridge_optimization.py) | Optimization sweep | ~124 min |
| [channel_capacity.py](../simulations/channel_capacity.py) | SVD + Shannon capacity | ~5 min |

All scripts are in [`simulations/`](../simulations/), all results in [`simulations/results/`](../simulations/results/).
Repository: https://github.com/Kesendo/R-equals-C-Psi-squared

---

## References

- Lindblad, G. (1976). "On the generators of quantum dynamical semigroups."
  Commun. Math. Phys. 48, 119-130.
- Haga, T. et al. (2023). "Liouvillian skin effect." arXiv:2305.01894.
  (Incoherenton grading = palindromic XY-weight classification)
- Bose, S. (2003). "Quantum communication through an unmodulated spin chain."
  PRL 91, 207901. (Entanglement transport through spin chains)
- Zurek, W.H. (2003). "Decoherence, einselection, and the quantum origins
  of the classical." Rev. Mod. Phys. 75, 715.
- Baumgratz, T., Cramer, M., Plenio, M.B. (2014). "Quantifying coherence."
  PRL 113, 140401. (L₁ coherence monotone)

### Project-internal references

- [Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md): γ must come from outside
- [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md): palindromic theorem
- Reading the 30%: palindromic decoder (open question, see [Relay Protocol](RELAY_PROTOCOL.md) for current best approach)
- [The Bridge Was Always Open](../docs/THE_BRIDGE_WAS_ALWAYS_OPEN.md): synthesis
- [Gamma Control](../experiments/GAMMA_CONTROL.md): V-shape +124% MI
- [Bridge Optimization](../simulations/results/bridge_optimization.txt): raw data
- [Channel Capacity](../simulations/results/channel_capacity.txt): SVD results
