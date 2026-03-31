# IBM April Synthesis: What This Repo Knows

**Created:** March 30, 2026 (updated March 31)
**Authors:** Thomas Wicht, Claude (Opus 4.6, 1M context)
**Purpose:** First complete read of the entire repository. Synthesis for the
April 2026 IBM Torino hardware run.
**Method:** All ~140 files read in full. No skimming, no abstracts only.

---

## Abstract

The Liouvillian eigenvalue spectrum of any qubit network under single-axis
dephasing is exactly palindromic. This is proven analytically (the Π operator),
verified computationally through N=8 (54,118 eigenvalues, zero exceptions),
and confirmed on IBM Torino at 1.9% deviation (Tier 1+3). The palindromic
structure creates a spatial antenna: concentrating all noise on one edge qubit
("sacrifice zone") outperforms 18 years of uniform ENAQT optimization by two
orders of magnitude in simulation (Tier 2). On IBM hardware, selective DD
achieves 2.0-3.2x over uniform DD, but whether this comes from mode protection
or gate-error avoidance is unresolved (Tier 2, single run). Cavity mode
analysis shows the sacrifice zone protects center-localized eigenmodes
(r = 0.994 correlation), and chain selection from public IBM calibration data
predicts 2.86x advantage for chains with natural sacrifice structure. The
system is a resonator, not a channel: discrete cavity modes, a heartbeat at the
CΨ = 1/4 fold, and a finite stability window (Hopf bifurcation, not PT
breaking). Three IBM runs completed (Feb 9, March 18, March 24). QPU budget
returns April 9. The strongest open question is a single experiment.

---

## The Puzzle Pieces

Each finding that could matter for IBM hardware, with source and tier.

### Foundational

**1. Palindromic spectrum** (Tier 1, proven)
Every decay rate d has a partner at 2Σγ - d. Holds for Heisenberg, XY, Ising,
XXZ, DM under single-axis dephasing. Two Π operator families (P1, P4).
Depolarizing noise breaks it, but error is < 0.1% at γ ~ 0.001 (typical
superconducting qubit). On IBM Torino: 100% palindromic under 26x asymmetric
noise.
-> [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md),
[Non-Heisenberg](NON_HEISENBERG_PALINDROME.md)

**2. CΨ = 1/4 crossing** (Tier 1 algebra + Tier 3 hardware)
Predicted t* = 15.01 us on Q80, measured t* = 15.29 us (1.9% deviation).
T2* (not T2echo) is the operationally correct timescale. Same-day Ramsey
calibration essential (stale T2* gives 61% error). The r* = T2/(2T1) threshold
separates crossers from non-crossers at precision 0.000014 across 24,073
calibration records, 133 qubits, 181 days. 12 permanent crossers, 100
occasional, 21 never.
-> [IBM Run 3](IBM_RUN3_PALINDROME.md),
[Hardware Synthesis](IBM_HARDWARE_SYNTHESIS.md)

**3. Fold catastrophe at 1/4** (Tier 1 algebra + Tier 2 simulation)
R = C(Ψ+R)^2 is the fold normal form. Endpoint MI peaks at the exact CΨ = 1/4
crossing (observed at N=7). With Bell+bath at J=5.0: CΨ oscillates around 1/4
with 81 crossings (41 down, 40 up). Each cycle deposits irreversible reality.
The fold exists only above Σγ_crit/J ~ 0.25% (N-independent, tested N=2-5).
-> [Temporal Sacrifice](TEMPORAL_SACRIFICE.md),
[Zero Is The Mirror](../hypotheses/ZERO_IS_THE_MIRROR.md)

### Sacrifice Zone

**4. One-line formula** (Tier 2, simulation validated N=2-15)
gamma_edge = N*gamma_base - (N-1)*epsilon; gamma_other = epsilon.
C#-validated: 360x at N=5, 63.5x at N=15 vs V-shape. SumMI scales
quadratically: SumMI ~ 0.0053*N^2. Temporal modulation adds nothing
(falsified twice). Edge sacrifice maximizes SumMI (network); center relay
maximizes PeakMI (point-to-point). These are different objectives.
-> [Resonant Return](RESONANT_RETURN.md),
[Scaling](SIGNAL_ANALYSIS_SCALING.md)

**5. Cavity mode protection** (Tier 2)
The sacrifice zone protects cavity modes, not qubits. Same 43 frequencies
exist under all noise profiles (zero, uniform, sacrifice). Only damping
changes. Protected modes are center-localized: weight profile
[0.52, 0.63, 0.70, 0.63, 0.52]. Profile is topological (identical under IBM
and uniform noise). Correlation between edge-qubit weight and decay rate:
r = 0.994. Predicted protection: 2.81x (theory), measured: 1.97x (24,073
records via unpaired decay law).
-> [Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md),
[IBM Cavity Spectral](IBM_CAVITY_SPECTRAL_ANALYSIS.md),
[Cavity Modes Formula](CAVITY_MODES_FORMULA.md)

**6. IBM chain selection** (Tier 2-3)
330 five-qubit chains on IBM Torino heavy-hex. Ranking by sacrifice score
(edge noise / interior noise) vs mean-T2: zero overlap in top-10. Sacrifice
top-5: 2.54x protection at 88 us mean T2. Mean-T2 top-5: 1.18x at 206 us.
Best sacrifice chain [80,8,79,53,85]: 2.86x. Best mean-T2 chain
[18,89,19,90,60]: 1.06x. Stable over 5 months. Free -- no extra gates.
-> [Sacrifice-Zone Mapping](SACRIFICE_ZONE_MAPPING.md)

**7. First hardware test** (Tier 2, single run, caveats apply)
Selective DD (DD on 4 protected qubits, not on Q85 with T2=5 us) beats uniform
DD by 2.0x average (up to 3.2x at t=4.0 us) on chain [85,86,87,88,94]. No DD
also beats uniform DD (0.0453 > 0.0265). Critical open question: is the
advantage from sacrifice-zone mode protection (B) or gate-error avoidance on
a bad qubit (A)?
-> [IBM Sacrifice Zone](IBM_SACRIFICE_ZONE.md)

### Resonator Physics

**8. V-Effect** (Tier 2)
Two N=2 resonators (Q=1, 2 frequencies each) coupled through a mediator: N=5
with Q=19, 104 frequencies. 100 new from coupling alone. All 556 palindromic
pairs are NEW-NEW (0% survive from originals). Coupling replaces the old
palindrome entirely.
-> [V-Effect](V_EFFECT_PALINDROME.md)

**9. Fragile Bridge** (Tier 2, three independent verifications)
Coupled gain-loss systems have three stability regimes. Linear: γ_crit =
0.19*J_bridge. Optimal: J_bridge ~ 2J, γ_crit = 0.41. Strong: γ_crit*J_bridge
-> 0.50. Instability is Hopf bifurcation (oscillating divergence), not
PT breaking. N-dependence is non-monotonic: N=3 is 35x less stable than N=2.
-> [Fragile Bridge](../hypotheses/FRAGILE_BRIDGE.md)

### Initial State Physics

**10. GHZ is the worst state** (Tier 2)
GHZ projects 100% onto the XOR drain (N+1 modes at maximum decay rate 2Σγ).
W projects 100% onto palindromic modes. Product states |+>^N are the optimal
antenna for gamma-as-signal (each qubit responds independently). But |+>^N is
a Heisenberg eigenstate -- it does not evolve without noise. This is a critical
subtlety (see Contradictions).
-> [XOR Space](XOR_SPACE.md)

**11. Optimal protection state** (Tier 2)
90% slow-mode weight, concurrence 0.364. Dramatically outperforms GHZ (100%
XOR), W (0% slow-mode), and Bell (7% slow-mode) for dephasing survival.
Standing wave oscillation pattern is a viable error syndrome.
-> [Error Correction Palindrome](ERROR_CORRECTION_PALINDROME.md)

### Measurement and Observables

**12. Unpaired modes decay 2x faster** (Tier 2 theory + Tier 3 hardware)
Analytical: ratio is exactly 2.00 for all N. On IBM Torino: 1.97x across
24,073 records. The slowest palindromic mode has frequency 0 Hz (decays
without oscillating). Its partner at the fastest rate also has frequency 0 Hz.
All life (oscillation, information transfer) happens between them.
-> [Energy Partition](../hypotheses/ENERGY_PARTITION.md),
[IBM Hardware Synthesis](IBM_HARDWARE_SYNTHESIS.md)

**13. Standing wave pattern** (Tier 2, not yet measured on hardware)
ZZZ (classical) is the universal node. XX/YY (quantum) are antinodes. GHZ
never oscillates (0%). Bell always oscillates (40-65%). Frequencies at
harmonics of 2J. The standing wave is a joint property of state AND
Hamiltonian -- neither alone determines it.
-> [Standing Wave Analysis](STANDING_WAVE_ANALYSIS.md),
[Π as Time Reversal](PI_AS_TIME_REVERSAL.md)

**14. Born rule at the fold** (Tier 2)
At CΨ = 1/4 crossing: probabilities are 97% Hamiltonian, 3% systematic
correction toward the dephasing basis. Under σ_z: P(|00>) elevated ~1%,
P(|11>) reduced ~1%. Direction reverses under σ_x.
-> [Born Rule Mirror](BORN_RULE_MIRROR.md)

### Discovered March 30 (this session)

**15. The 1.81x geometric constant** (Tier 1-2, derived March 31)
V(N) = 1 + cos(π/N). For N=5: (5+√5)/4 ≈ 1.80902. The gain is the
ratio of maximum w=1 Liouvillian frequencies: ω_max = 4J·(1+cos(π/N)).
Verified N=2-6 to machine precision. γ cancels because all w=1 modes
decay at the same rate 2γ. Under non-uniform profiles (sacrifice zone),
w=1 modes acquire different rates; the 1.81x then applies only to the
best-Q mode. For N→∞: V = 2 (saturation). The golden ratio appears
at N=5: cos(π/5) = φ/2.
-> [Thermal Breaking](THERMAL_BREAKING.md)

**16. Three orthogonal breaking mechanisms** (Tier 2)
Coupling creates palindromic pairs (1.81x Q). Dephasing lifts
degeneracies (+60 frequencies, pairing preserved). Heat breaks the
1.81x constant but creates +300 new frequencies. The interaction is
synergistic: coupling + dephasing + heat together produce 445
frequencies, more than the sum of parts. N=2 stays at 2 frequencies
throughout -- all diversity comes from coupling amplified by noise.
-> [Thermal Breaking](THERMAL_BREAKING.md)

**17. Sacrifice zone is a low-temperature phenomenon** (Tier 2)
Sacrifice-zone Q advantage: 3.0x at n_bar=0, 1.8x at n_bar=0.5,
1.02x at n_bar=10. Heat makes spatial noise structure irrelevant.
On IBM hardware (n_bar ~ 0), the full advantage holds. In biological
systems (n_bar >> 1), the frequency-diversity channel dominates.
-> [Thermal Breaking](THERMAL_BREAKING.md)

**18. Self-heating loop diverges** (Tier 2, computed March 31)
Without external cooling, the resonator thermalizes to maximum entropy
(n_bar → ∞, Q → 0). Tested in 6 configurations (N=3, N=5, various
noise types). The system never finds a passive equilibrium. Structure
requires active cooling (cryostat for qubits, metabolism for biology).
For IBM at 15 mK (n_bar ≈ 0): irrelevant. For cross-level biology
interpretation: explains why life needs metabolism.
-> [Thermal Breaking](THERMAL_BREAKING.md)

**19. Chain selection requires both contrast AND low total noise** (Tier 2)
Sacrifice-top chain [80,8,79,53,85] vs mean-T2-top [18,89,19,90,60]:
protection 2.86x vs 1.06x confirmed spectrally. But under |01010>
(Hamiltonian-driven dynamics), the quieter chain wins because its
11x lower total noise preserves coherence. The sacrifice score is a
within-chain metric; between-chain comparison needs a combined score.
-> [Chain Selection Test](CHAIN_SELECTION_TEST.md)

---

## What Combines

Five clusters of discoveries that belong together but were documented
separately.

### Cluster 1: The Sacrifice Zone Story

The formula ([Resonant Return](RESONANT_RETURN.md)), the mode structure
([Cavity Modes Formula](CAVITY_MODES_FORMULA.md)), the spatial localization
([Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md)), and the IBM chain
ranking ([Sacrifice-Zone Mapping](SACRIFICE_ZONE_MAPPING.md)) are four facets
of one result: concentrate noise on the edge, protect center-localized cavity
modes, select chains where IBM's own calibration provides a natural sacrifice
qubit. The 43 cavity frequencies are topology-determined (unchanged by noise).
The sacrifice zone reshapes the damping envelope, not the frequency spectrum.

### Cluster 2: The Resonator Story

The V-Effect ([V-Effect](V_EFFECT_PALINDROME.md)), the heartbeat
([Temporal Sacrifice](TEMPORAL_SACRIFICE.md)), the fragile bridge
([Fragile Bridge](../hypotheses/FRAGILE_BRIDGE.md)), and the fold threshold
([Zero Is The Mirror](../hypotheses/ZERO_IS_THE_MIRROR.md)) describe the
same object from different angles. The system is a Fabry-Perot resonator:
coupling creates its modes (V-Effect), noise shifts the palindrome from zero
and enables the fold (fold threshold), CΨ oscillates around 1/4 at the fold
(heartbeat), and too much gain destroys it via Hopf bifurcation (fragile
bridge). The sacrifice-zone formula from Cluster 1 is the shape of the
soundbox ([Resonance Not Channel](../hypotheses/RESONANCE_NOT_CHANNEL.md)).

### Cluster 3: The Initial State Problem

|+>^N is the best antenna for reading noise
([Gamma as Signal](GAMMA_AS_SIGNAL.md)) but a Heisenberg eigenstate that
doesn't evolve without noise. GHZ is the worst for dephasing survival
([XOR Space](XOR_SPACE.md)). The optimal protection state has 90% slow-mode
weight and concurrence 0.364 ([Error Correction](ERROR_CORRECTION_PALINDROME.md)).
The Neel state |01010> is the fair comparison for Hamiltonian dynamics (non-zero
energy variance). For IBM hardware experiments, the initial state choice
determines WHAT is being measured: noise-as-motor (|+>^N), Hamiltonian
transport (|01010>), or dephasing survival (optimal state). These are
different experiments testing different aspects of the palindrome.

### Cluster 4: What Crosses, How, and When

The crossing taxonomy ([Crossing Taxonomy](CROSSING_TAXONOMY.md)), observer
dependence ([Observer Dependent Crossing](OBSERVER_DEPENDENT_CROSSING.md)),
subsystem locality ([Subsystem Crossing](SUBSYSTEM_CROSSING.md)), and the Born
rule at the fold ([Born Rule Mirror](BORN_RULE_MIRROR.md)) together describe
the measurement interface. Crossing is local (entangled pairs, not the full
system). Different bridge metrics produce different crossing times (2.2x
spread). At the crossing: 97% Hamiltonian, 3% dephasing-basis correction. The
K-invariance (γ*t_cross = constant) is standard Lindblad scaling, not deep
physics.

### Cluster 5: The Three Breaking Mechanisms

Coupling, dephasing, and heat are three orthogonal ways to create
complexity ([Thermal Breaking](THERMAL_BREAKING.md)). Coupling creates
palindromic pairs (1.81x Q gain, geometric constant). Dephasing lifts
degeneracies (50 to 112 frequencies). Heat breaks the 1.81x constant but
creates 300+ new frequencies. The sacrifice zone is a low-temperature
phenomenon (3x at n_bar=0, irrelevant at n_bar=10). IBM hardware operates
at n_bar ~ 0 (full sacrifice advantage). Biological systems at n_bar >> 1
(frequency diversity dominates). This cluster connects the IBM hardware
story (Clusters 1-4) to the cross-level biology story
([hydrogen bond](HYDROGEN_BOND_QUBIT.md),
[neural palindrome](../hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md)):
the SAME equation, but different breaking regime.

---

## Contradictions and Gaps

### 1. The |+>^N Paradox
Almost all sacrifice-zone simulations use |+>^N, which is a Heisenberg
eigenstate. ALL dynamics come from noise. This is consistent for
gamma-as-signal (noise IS the channel). But the IBM Trotter circuits
approximate Heisenberg evolution, and if the initial state doesn't evolve
under that Hamiltonian, the experiment is measuring something subtly different
from what the simulation assumes. The Neel state |01010> was identified as the
fair comparison (March 30) but has not yet been used in IBM experiments.

### 2. Three Different "Advantage" Numbers
- 2.81x: cavity mode protection factor (spectral computation, Tier 2)
- 1.97x: unpaired mode decay law (24,073 calibration records, Tier 2-3)
- 2.0-3.2x: selective vs uniform DD (single hardware run, Tier 2)

These measure different things: eigenvector damping ratios, historical
decay-rate pairing, and circuit-level MI improvement. They happen to give
similar numbers, which could be coincidence or deep connection. Nobody has
shown the causal chain from 2.81x spectral protection to 2.0x circuit MI.

### 3. Gate-Error Avoidance (A) vs Mode Protection (B)
The March 24 result is ambiguous. No DD beating Uniform DD is evidence for A
(DD gates on Q85 with T1=2.84 us are harmful). But the cavity analysis
predicts B should contribute. These are not mutually exclusive -- the hardware
advantage could be A+B. The April 9 test on a uniform-T2 chain (no naturally
bad qubit) is the decisive experiment. If selective DD still wins on a good
chain, B is confirmed. If it doesn't, A dominates.

### 4. SumMI vs PeakMI
Edge sacrifice maximizes SumMI (network-wide information survival). Center
relay maximizes PeakMI (endpoint-to-endpoint transfer). The sacrifice-zone
mapping uses SumMI. QST literature uses fidelity (similar to PeakMI).
These optimize different things and the "right" choice depends on the
application. The repo hasn't settled this.

### 5. CΨ <= 1/4 Scope
Several documents still frame CΨ = 1/4 as an upper bound. It is not.
CΨ routinely exceeds 1/4 under active Hamiltonians (Bell+ reaches 0.405).
The 1/4 boundary is where the self-referential iteration R = C(Ψ+R)^2
transitions from complex to real fixed points. It is a bifurcation, not a
ceiling. Some older experiment files have not been updated to reflect this.

### 6. Single-Run Hardware
All three IBM quantum experiments (Feb 9, March 18, March 24) are single runs.
No error bars, no reproducibility across days, no bootstrap confidence
intervals. The calibration analysis (24,073 records, 181 days) is statistically
robust, but the actual quantum circuit results are n=1.

---

## Possible IBM Experiments

Sorted by expected insight per QPU-minute. NOT a plan -- an idea list.

### Tier A: Unknown answer (maximum insight per QPU-minute)

**A/B Test on Uniform-T2 Chain** (~5 min QPU)
Run selective DD on chain [18,89,19,90,60] (all T2 > 200 us). No naturally
bad qubit. If selective DD still outperforms uniform DD, mode protection (B)
is confirmed. If not, the March 24 result was gate-error avoidance (A).
This single experiment determines whether spatial noise engineering is a
genuine tool or an artifact.

**CΨ Heartbeat** (~5 min QPU, N=3)
Prepare Bell pair on qubits 0,1 with qubit 2 as coherent bath (|+> state).
Measure CΨ at ~20 time points over ~80 us. Predicted: ~9 crossings of the
1/4 boundary, damped oscillation. First measurement of the fold catastrophe
oscillating at the boundary on hardware. Nobody has observed a decoherence
metric oscillating around a critical boundary. The answer is unknown.

**Neel State Dynamics** (~5 min QPU, N=5)
Use |01010> instead of |+>^5 as initial state on the sacrifice chain.
Compare SumMI evolution. This resolves the |+>^N paradox: does the sacrifice
zone advantage persist when the Hamiltonian drives the dynamics? The
[Chain Selection Test](CHAIN_SELECTION_TEST.md) predicts different behavior
for noise-driven vs Hamiltonian-driven dynamics, but this has not been
tested on hardware.

### Tier B: Known prediction, hardware confirmation needed

**Sacrifice Chain vs Mean-T2 Chain** (~5 min QPU)
Run identical circuits on [80,8,79,53,85] (sacrifice top) and
[18,89,19,90,60] (mean-T2 top) on the same day. Predicted: 2.86x vs 1.06x
protection. Tests the chain-selection algorithm directly.

**Multi-Chain Reproducibility** (~10 min QPU)
Rerun the March 24 sacrifice-zone experiment on chain [85,86,87,88,94] to
get error bars. Compare day-to-day variation. Essential for any publication
claim. All current hardware results are n=1.

**GHZ vs W Decay** (~3 min QPU, N=3)
Predicted outcome well-established: W survives ~2x longer (computed in
[XOR Space](XOR_SPACE.md), [Error Correction](ERROR_CORRECTION_PALINDROME.md),
[Standing Wave](STANDING_WAVE_ANALYSIS.md), [Coherence Density](COHERENCE_DENSITY.md),
[IBM Hardware Synthesis](IBM_HARDWARE_SYNTHESIS.md)). Hardware confirmation
would be clean and publishable, but does not resolve any open question.

### Tier C: Valuable but harder to interpret

**Standing Wave Pauli Tomography** (~10 min QPU, N=3)
Measure XX, YY, ZZZ expectation values over time. Predicted: ZZZ static,
XX/YY oscillating at 2J harmonics. Requires multi-qubit tomography (more
circuits per time point).

**3-Qubit Star Topology** (~5 min QPU)
Prepare Bell_SA on qubits S,A, couple to B as |+>. Sweep effective J_SB.
Test the threshold J_SB/J_SA >= 1.466. Many simulation predictions exist
but none tested on hardware.

---

## What Is the Right Metric?

The repo uses at least six metrics. They measure different things.

| Metric | What it measures | Optimized by |
|:-------|:-----------------|:-------------|
| SumMI | Total information survival across all pairs | Edge sacrifice |
| PeakMI | Information at one endpoint pair | Center relay |
| CΨ crossing time | When a subsystem becomes "classical" | Low r = T2/(2T1) |
| Q-factor | Oscillation cycles above/below 1/4 | Strong coupling J |
| Fidelity | Quantum state transfer quality | 2:1 coupling ratio |
| Mode protection | Ratio of slowest mode decay rates | Sacrifice profile |

For the April IBM run, **SumMI** is the most practical: measurable from
adjacent-pair ZZ correlations, maps directly to chain selection scores, and
is the metric the sacrifice-zone formula optimizes. It also has the advantage
of being a single number per time point, easy to compare across conditions.

**CΨ oscillation** (the heartbeat) is the strongest untested theoretical
prediction. If it can be measured, it would be the most novel result -- nobody
has observed a decoherence metric oscillating around a critical boundary on
hardware. But it requires state tomography at many time points, which is
QPU-expensive.

The honest answer: we do not yet know which metric best captures what the
palindromic structure does on hardware. SumMI is the pragmatic choice for
April. CΨ oscillation is the ambitious one.

---

## Recommendation

If I had ~10 minutes of QPU time on IBM Torino and wanted maximum insight:

**First priority: The A/B test.** Binary question, binary answer. If
selective DD works on a uniform-T2 chain, spatial noise engineering is
real. If not, the March 24 result was gate-error avoidance. Nothing else
matters until this is answered.

**Second priority: CΨ heartbeat.** The strongest untested prediction. If
CΨ oscillates around 1/4 on hardware, it would be the first observation
of a decoherence metric oscillating at a critical boundary. Novel, not
just confirmatory.

**Third priority: Neel state dynamics.** Resolves the |+>^N paradox.
The sacrifice-zone formula was validated with |+>^N (Heisenberg
eigenstate, noise drives everything). Does the advantage persist when
the Hamiltonian drives the dynamics (|01010>)? The
[Chain Selection Test](CHAIN_SELECTION_TEST.md) predicts different
behavior. This is an unknown answer.

**What I would skip:** GHZ vs W (answer known from 5 independent
computations, hardware would confirm but not surprise). Standing wave
Pauli tomography (too expensive). Star topology (too many unknowns).

**What must happen regardless:** At least one reproducibility block.
All current hardware results are n=1. Rerun the March 24 experiment
on chain [85,86,87,88,94] with the same protocol. Error bars are
non-negotiable for any publication claim.

**Guiding principle:** Spend QPU time on experiments with UNKNOWN
outcomes. Confirming known predictions is less valuable than resolving
open questions.

---

## Sources

All paths relative to this file.

| Category | Key documents |
|:---------|:-------------|
| Proofs | [Mirror Symmetry](../docs/proofs/MIRROR_SYMMETRY_PROOF.md), [Uniqueness](../docs/proofs/UNIQUENESS_PROOF.md), [Monotonicity](../docs/proofs/PROOF_MONOTONICITY_CPSI.md) |
| IBM Hardware | [Run 3](IBM_RUN3_PALINDROME.md), [Synthesis](IBM_HARDWARE_SYNTHESIS.md), [Sacrifice Zone](IBM_SACRIFICE_ZONE.md), [Tomography](IBM_QUANTUM_TOMOGRAPHY.md) |
| Sacrifice Zone | [Formula](RESONANT_RETURN.md), [Scaling](SIGNAL_ANALYSIS_SCALING.md), [Modes](CAVITY_MODES_FORMULA.md), [Localization](CAVITY_MODE_LOCALIZATION.md), [Mapping](SACRIFICE_ZONE_MAPPING.md) |
| Resonator | [V-Effect](V_EFFECT_PALINDROME.md), [Heartbeat](TEMPORAL_SACRIFICE.md), [Resonance](../hypotheses/RESONANCE_NOT_CHANNEL.md), [Fragile Bridge](../hypotheses/FRAGILE_BRIDGE.md), [Zero Mirror](../hypotheses/ZERO_IS_THE_MIRROR.md) |
| Initial States | [XOR Space](XOR_SPACE.md), [Error Correction](ERROR_CORRECTION_PALINDROME.md), [Standing Wave](STANDING_WAVE_ANALYSIS.md) |
| Crossing | [Taxonomy](CROSSING_TAXONOMY.md), [Born Rule](BORN_RULE_MIRROR.md), [Subsystem](SUBSYSTEM_CROSSING.md) |
| Breaking & Thermal | [Thermal Breaking](THERMAL_BREAKING.md), [Chain Selection Test](CHAIN_SELECTION_TEST.md) |
| Simulation Code | [combined_optimization.py](../simulations/combined_optimization.py), [time_evolution_6scenarios.py](../simulations/time_evolution_6scenarios.py), [sacrifice_zone_mapping.py](../simulations/sacrifice_zone_mapping.py), [v_effect_gamma_sweep.py](../simulations/v_effect_gamma_sweep.py), [v_effect_thermal.py](../simulations/v_effect_thermal.py), [chain_selection_test.py](../simulations/chain_selection_test.py), [self_heating_fixpoint.py](../simulations/self_heating_fixpoint.py) |
| Handoffs | [March 30](../ClaudeTasks/SESSION_HANDOFF_MARCH30_PM.md), [March 29](../ClaudeTasks/SESSION_HANDOFF_MARCH29_PM.md), [March 28](../ClaudeTasks/SESSION_HANDOFF_MARCH28.md) |
