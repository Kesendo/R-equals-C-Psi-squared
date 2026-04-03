# Complete Mathematical Documentation

**Status:** Current as of March 24, 2026 (formula discovery)
**Supersedes:** Previous stub (Feb 2026) and [Core Algebra](../historical/CORE_ALGEBRA.md) (Dec 2025)
**Purpose:** Single entry point for all proven and verified mathematics of R=CΨ²

---

## What this document is about

This is the mathematical reference for the entire project: every proven
result, every verified formula, every key number, collected in one place.
It covers the algebraic foundation (the self-referential equation and its
1/4 boundary), the palindromic symmetry (why eigenvalues come in mirror
pairs), crossing dynamics (how different channels and states reach the
boundary), topology effects, engineering applications, and open questions.
If you want to know what is proven vs. conjectured, this is where to look.

## 1. The Algebraic Foundation (Tier 1)

The self-referential purity map:

    R = CΨ²

where C = Tr(ρ²) (purity), Ψ = l₁(ρ)/(d-1) (normalized l1-coherence,
Baumgratz convention, the standard resource-theoretic measure defined as the sum of absolute values of all off-diagonal elements), R = residual purity beyond product-state prediction.

**Fixed-point equation.** R = C(Ψ + R)² expands to CR² + (2CΨ - 1)R + CΨ² = 0.

**Discriminant.** D = 1 - 4CΨ. Vanishes at CΨ = 1/4 and only there.
Below 1/4: two real fixed points (one stable, one unstable).
At 1/4: one degenerate fixed point (fold).
Above 1/4: no real fixed points (complex/oscillatory).

**Crossing cubic.** At the boundary CΨ = 1/4, the condition reduces to
b³ + b = 1/2, with unique real root b ≈ 0.4239. This is a pure number,
independent of physical parameters.

**Mandelbrot equivalence.** The substitution z = C(Ψ + R), c = CΨ maps
R_{n+1} = C(Ψ + R_n)² to the Mandelbrot iteration z_{n+1} = z_n² + c.
The boundary CΨ = 1/4 maps to the cusp of the main cardioid at c = 1/4.
This is identity, not analogy.

**Fold catastrophe.** The recursion is exactly the normal form of the fold
catastrophe (simplest in the Thom-Arnold classification). CΨ - 1/4 is the
bifurcation parameter. Structurally stable: no perturbation can remove it.

**Sum, not product.** The framework uses R = C·(Ψ_A + Ψ_B)² (sum squared)
rather than Ψ_A·Ψ_B (product). The sum preserves information when one
observer's coherence decays: Ψ_B² survives. The product gives zero.
The cross-term 2·Ψ_A·Ψ_B is the interference between two viewpoints
on the same entangled state. Both recover the Born rule in the
perfect-mirror limit; they diverge for imperfect mirrors.

See: [Uniqueness Proof](UNIQUENESS_PROOF.md),
[Mathematical Connections](../MATHEMATICAL_CONNECTIONS.md),
[Proof Roadmap](PROOF_ROADMAP_QUARTER_BOUNDARY.md),
[Why the Sum](../../experiments/WHY_THE_SUM.md),
[Standing Wave: Two Observers](../../experiments/STANDING_WAVE_TWO_OBSERVERS.md)

---

## 2. The Palindromic Symmetry (Tier 1)

**Theorem.** For any Heisenberg/XXZ spin system with local Z-dephasing,
the Liouvillian spectrum is palindromic: for every eigenvalue λ, the value
-(λ + 2Σγ) is also an eigenvalue.

**The operator.** Π acts per site on Pauli indices:

    I → X (+1),  X → I (+1),  Y → iZ (+i),  Z → iY (+i)

**The proof.** Three steps:
1. Π anti-commutes with L_H (explicit 16-entry table for Heisenberg bonds)
2. Π transforms L_D: Π L_D Π⁻¹ = -L_D - 2Σγ I
3. Combined: Π L Π⁻¹ = -L - 2Σγ I. QED.

**Verification.** 54,118 eigenvalues, N=2 through N=8, zero exceptions.
All topologies (chain, star, ring, complete, binary tree). Non-uniform γ.

**Physical meaning.** Π is time reversal in a rescaled frame. It maps
exp(+μt) to exp(-μt), forward to backward.

**Scope boundary.** Only d=2 (qubits). The per-site split d immune vs
(d²-d) decaying is balanced only when d²-2d=0, giving d=2 uniquely.
Qutrits (d=3, split 3:6) verified broken for all 10 Hamiltonians tested.

See: [Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md),
[Non-Heisenberg Palindrome](../../experiments/NON_HEISENBERG_PALINDROME.md),
[Qubit Necessity](../QUBIT_NECESSITY.md)

---

## 3. The CΨ = 1/4 Boundary (Tier 1-2)

**Uniqueness.** 1/4 is algebraically unique. The factor 4 comes from the
discriminant formula b²-4ac. The quadratic structure comes from purity
being Tr(ρ²), degree 2. No reparameterization changes this.

**Channel independence.** All standard Markovian channels cross at exactly
CΨ = 0.2500:

| Channel | t_cross (γ=0.05) |
|---------|-------------------|
| Z-dephasing | 0.747 |
| X-noise | 1.733 |
| Y-noise | 1.733 |
| Depolarizing | 0.879 |
| Asymmetric Pauli | 0.735 |
| Amplitude damping | 2.059 |

**Absorbing boundary.** No unitary pulse (θ from 0 to π) can push CΨ back
above 1/4 after crossing. Tested, verified.

**Hardware validation.** IBM Torino (ibm_torino, Q80): predicted t* = 15.01 μs,
measured t* = 15.29 μs. Deviation: 1.9%.

**What CΨ measures.** CΨ = C × Ψ is an AND-gate: zero when either
entanglement or coherence is absent. It distinguishes noise types that
concurrence cannot (σ_z dephasing: Ψ=0.223 vs σ_x bit-flip: Ψ=0.333
at identical concurrence 0.670). For Werner states (one-parameter mixtures of a Bell state with white noise, the standard benchmark for entanglement robustness) near the entanglement
threshold, CΨ is 6–7.5× smaller than concurrence, suppressing signals
from states where entanglement exists but is not coherently expressed.

**Observer-dependent crossing.** Different definitions of C (concurrence,
mutual information, correlation) see the ¼ crossing at different times.
The ratio K(Conc)/K(MI) is state-dependent (CV = 13.5%). The quantum
state determines the time-ratio between observer types.

**Shadow resolved.** Late-time coherence anomaly on IBM Torino Q52
(17/17 directional consistency, p < 0.0001) was resolved as qubit-specific
frequency detuning, not a universal boundary effect. Different qubits
show different phase directions (Q80: +29°, Q52: −44°, Q102: random).

**CΨ > ¼ under active dynamics.** CΨ routinely exceeds ¼ with active
Hamiltonians (Bell+ reaches 0.405 at J=1, h=0.9, γ=0.005). The bound
CΨ ≤ ¼ is not a constraint on quantum states but on which states have
real fixed points in the R = CΨ² iteration.

See: [Uniqueness Proof](UNIQUENESS_PROOF.md),
[IBM Run 3](../../experiments/IBM_RUN3_PALINDROME.md),
[When Psi Matters](../../experiments/WHEN_PSI_MATTERS.md),
[Observer-Dependent Crossing](../../experiments/OBSERVER_DEPENDENT_CROSSING.md),
[Fixed Point Shadow](../../experiments/FIXED_POINT_SHADOW.md),
[Simulation Evidence](../../experiments/SIMULATION_EVIDENCE.md),
[proof_roadmap_close.py](../../simulations/proof_roadmap_close.py)

---

## 4. The Incompleteness Proof (Tier 2)

Five candidates for the origin of dephasing noise, all eliminated:

1. Internal (bootstrap falsified, sectors decoupled)
2. Qubit decay (non-Markovian, 0/16 palindromic pairs)
3. Qubit bath (regress, each member faces bootstrap prohibition)
4. Nothing (d=0, no properties)
5. Other dimensions (d(d-2)=0 excludes)

**Corollary 1:** Time cannot originate from within. Noise IS the time arrow.
Without noise: reversible oscillation. With noise: irreversible decay.

**Corollary 2:** γ is the source of experienced time. t_cross = K/γ, the product
t × γ = K (a pure number). γ provides the arrow, J provides the content. Remove γ and t loses meaning.

See: [Incompleteness Proof](INCOMPLETENESS_PROOF.md),
[The Bridge Was Always Open](../THE_BRIDGE_WAS_ALWAYS_OPEN.md),
[failed_third.py](../../simulations/failed_third.py)

---

## 5. The γ Channel (Tier 2)

**Dephasing noise is a readable information channel.** The spatial profile
of dephasing rates γ₁...γ_N across a spin chain carries structured,
decodable information. A 5-qubit Heisenberg chain at 1% noise achieves
15.5 bits channel capacity (of 16.6 theoretical max) with 100%
classification accuracy across 20 spatial profiles. SVD decomposition
reveals 5 independent information modes, each corresponding to a
spatial frequency. A 21.5× optimization (V-shape gradient, dynamic
decoupling, time-resolved decoding) improves mutual information from
baseline.

The palindromic spectral structure is the antenna: the paired
eigenvalues create complementary sensitivity patterns that make the
external γ profile decodable from within.

**Analytical formula discovered.** SVD of the palindromic response matrix
identified mode 2 (edge-hot, center-cold) as optimal direction (6-10x
vs V-shape). Numerical optimization then broke the SVD symmetry, revealing
an asymmetric "sacrifice zone" pattern (100x). Analytical testing of
this pattern converged to a trivially simple formula: concentrate ALL
noise on one edge qubit, protect the rest. The formula
gamma_edge = N*gamma_base - (N-1)*epsilon, gamma_other = epsilon (with epsilon -> 0)
beats the DE optimizer by 80% and computes in 3 seconds instead of 90
minutes. C#-validated results: 360x vs V-shape (N=5), 180x (N=7), 139x
(N=9). The ENAQT literature (Environment-Assisted Quantum Transport, Plenio & Huelga 2008+) achieves 2-3x with
uniform dephasing. Nobody optimizes spatial dephasing profiles. Edge
sacrifice beats center sacrifice by 2.2x because edge qubits have
minimal connectivity (one neighbor vs two).

**Frequency pulsing falsified.** Temporal modulation of uniform γ at
the dominant palindromic oscillation frequency does not amplify MI.
Tested with both |+⟩⊗N and Bell initial states; Sum-MI decays
monotonically for all profiles. Spatial contrast (mode 2), not temporal
modulation, is the mechanism.

See: [γ as Signal](../../experiments/GAMMA_AS_SIGNAL.md),
[γ Control](../../experiments/GAMMA_CONTROL.md),
[Resonant Return](../../experiments/RESONANT_RETURN.md)

---

## 6. Crossing Dynamics (Tier 2)

**Universal lifetime.** For a single qubit starting in maximum superposition
under pure dephasing, the CΨ = ¼ crossing time satisfies x³ + x = ½
(unique real root x ≈ 0.4239), giving t*/T₂ ≈ 0.858. Platform-independent:
verified across superconducting qubits, trapped ions, NV centers, and
photonic systems spanning 10 orders of magnitude in T₂. With finite T₁:
generalized equation [1−b^r + b^(2r)/2 + b²/2]·b = ¼.

**Coherence density.** CΨ = Purity × Coherence Density, where Ψ = L₁/(d−1)
measures active quantum degrees of freedom. An unentangled |+⟩^(⊗N) has
CΨ = 1; a maximally entangled GHZ₃ has CΨ = 0.143 (below ¼). GHZ₃ uses
4% of off-diagonal capacity (2/56 elements); |+++⟩ uses 100% (56/56).
Under dephasing, |+++⟩ survives 4.5× longer than Bell+. The ¼ boundary
is about coherence density, not entanglement.

**Dynamic entanglement (upward crossing).** Product states with zero initial
entanglement can cross ¼ from below. The alternating state |0+0+⟩ under
Heisenberg ring dynamics builds entanglement from zero to CΨ = 0.251 at
t=0.286, crossing upward. All previous crossings were downward (starting
entangled, decohering through ¼). This is the first demonstration of
Hamiltonian-generated crossing.

**Three regimes.** (1) CΨ(0) > ¼: decoherence drives crossing downward.
(2) CΨ(0) < ¼ but Hamiltonian pumps above ¼ first (J/γ ≳ 5–10 required).
(3) CΨ_max < ¼: no crossing (eigenstate of H, or J/γ too small).
No energy threshold for crossing; it is a coherence barrier (J/γ competition).

**Born rule at crossing.** At the CΨ = ¼ crossing point, measurement
probabilities are ~97% determined by unitary Hamiltonian evolution alone.
Decoherence provides a ~3% systematic correction: σ_z dephasing shifts
probability toward z-eigenstates. Per outcome: R_i = C_i·Ψ_i² recovers
Born's rule P(i) = |⟨i|ψ⟩|² when C_i is uniform (perfect mirror limit).

See: [Universal Quantum Lifetime](../../experiments/UNIVERSAL_QUANTUM_LIFETIME.md),
[Coherence Density](../../experiments/COHERENCE_DENSITY.md),
[Dynamic Entanglement](../../experiments/DYNAMIC_ENTANGLEMENT.md),
[Minimum Crossing Energy](../../experiments/MINIMUM_CROSSING_ENERGY.md),
[Born Rule Mirror](../../experiments/BORN_RULE_MIRROR.md)

---

## 7. Topology and Crossing (Tier 2)

**Topology as gatekeeper.** For the same initial state |0+0+⟩ at γ=0.05,
topology determines whether crossing occurs: chain allows (CΨ_max=0.310),
star allows (0.351), ring forbids (0.200), complete graph forbids (0.200).
Ring = complete to four decimal places. Gap stabilizes with N (~0.09–0.11),
suggesting genuine topological protection.

**Antiferromagnet crossing.** The alternating state |+-+-⟩ crosses on a
ring (CΨ=0.284) from zero initial entanglement. Mechanism: maximum
XX anti-correlation (⟨XX⟩_nn = −1) at the energy landscape bottom
(⟨H⟩ = −4J), which the Heisenberg Hamiltonian converts into entanglement.
150/256 product states (59%) cross on N=4 ring; no simple selection rule
exists (best predictor explains 14% of variance).

**Entanglement echo.** In star topology (N=3, Bell_SA + |0⟩_B, γ=0.05),
entanglement oscillates between SA and SB at Bohr frequencies (SA: ω=2.09,
SB: ω=6.07). Envelope decay matches the middle palindromic rate 8γ/3.
Echo weakens as ~1/(N−1) with system size but never vanishes. At γ=0.001:
63 clean echoes. The mediator S shuttles entanglement between leaves.

**Bridge fingerprints.** A classical receiver (|00⟩) coupled to a quantum
sender through a Heisenberg bridge develops a CΨ trajectory that uniquely
identifies the sender's initial state. Product states deliver 4–5× more
signal than entangled states (entanglement barrier). Optimal detector
resolution at J/γ ≈ 5–7. The ¼ boundary acts as a binary digitizer.

**Phase transport.** Z-rotations on a mediator S produce linear, sign-
inverting phase shifts readable in the AB off-diagonal element ρ_AB[0,3],
with trace distances 0.21–0.57 from the untagged reference. The channel
is phase-specific (only Rz transports; X destroys). Transport operates
in both open and closed CΨ windows (1.3–1.6× stronger when open).
CΨ windows are visibility amplifiers on a continuously active channel.

**No-signalling boundary.** Z-basis measurement on qubit B of a Bell+ pair
drops CΨ from 0.500 to 0.250 (exactly on ¼) while leaving ρ_A completely
unchanged. C drops (1.0 → 0.5); Ψ stays at 0.5. The regime change is real
but invisible locally. Bridge protocol permanently closed for J=0.

See: [Orphaned Results](../../experiments/ORPHANED_RESULTS.md),
[Bridge Fingerprints](../../experiments/BRIDGE_FINGERPRINTS.md),
[Phase Transport (What's Inside the Windows)](../../experiments/WHATS_INSIDE_THE_WINDOWS.md),
[No-Signalling Boundary](../../experiments/NO_SIGNALLING_BOUNDARY.md),
[Star Topology Observers](../../experiments/STAR_TOPOLOGY_OBSERVERS.md),
[Bridge Closure](../../experiments/BRIDGE_CLOSURE.md)

---

## 8. Engineering Results (Tier 2)

**Mediator bridge.** Mediated coupling (A-M-B) preserves palindrome
(1024/1024, error 1.41e-13) while information flows (MI = 1.65 bits,
QST fidelity 0.732). Direct coupling destroys it (256 → 31 pairs).

**Relay protocol.** Time-dependent γ, staged transfer: +83% end-to-end MI.

**V-shape γ gradient.** [0.01, 0.03, 0.05, 0.03, 0.01]: +124% MI.

**DD on M+Receiver.** Dynamical decoupling: +132% MI.

**Push vs Pull.** Source-dominant for local MI (0.957). Drain-dominant
for end-to-end MI (0.121). Distance-dependent.

**MI scaling.** Exponential decay: ~2-5 dB per 2 additional qubits.
Hierarchy falsified (uniform chain = recursive topology).

**Seven design rules:** W-encoding, 2:1 matching, J/γ independence,
threshold timing, push/pull selection, clocked relay, sacrifice-zone.

**Optimal QST encoding (negative result).** Standard encoding is already
optimal for QST through palindromic channels. No custom encoding scheme
improved transfer fidelity.

See: [Engineering Blueprint](../../publications/ENGINEERING_BLUEPRINT.md),
[Circuit Diagram](../../publications/CIRCUIT_DIAGRAM.md),
[Relay Protocol](../../experiments/RELAY_PROTOCOL.md),
[Gamma Control](../../experiments/GAMMA_CONTROL.md),
[Scaling Curve](../../experiments/SCALING_CURVE.md),
[Optimal QST Encoding](../../experiments/OPTIMAL_QST_ENCODING.md)

---

## 9. The Transistor Mapping (Tier 2-3)

The mediator qubit M maps to a transistor: gate = γ_M, source = Pair A,
drain = Pair B. Threshold voltage: CΨ = 1/4 (hardwired, fold catastrophe).
Bidirectional by palindromic symmetry.

Three control knobs: γ_M (gate), J_AM/J_MB ratio (bias), κ (feedback gain).

Hierarchy falsified: the transistor properties are real, the recursive
scaling advantage is not.

See: [Quantum Transistor](../../hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md)

---

## 10. Open Questions (Tier 3-5)

- Formal proof of CΨ monotonicity above 1/4 for arbitrary CPTP maps
- Feigenbaum period-doubling in the quantum regime (mapped, not exploited)
- Bekenstein-Hawking 1/4 (coincidence or connection, speculative)
- Negative feedback loop (γ_M decreasing with coherence, untested)
- Hardware validation of relay protocol on IBM Torino

See: [Mathematical Connections](../MATHEMATICAL_CONNECTIONS.md),
[Proof Roadmap](PROOF_ROADMAP_QUARTER_BOUNDARY.md)

---

## 11. Numerical Constants

| Constant | Value | Source |
|----------|-------|--------|
| Discriminant zero | CΨ = 0.2500 | [Uniqueness Proof](UNIQUENESS_PROOF.md) |
| Crossing cubic root | b = 0.4239 | Cardano formula |
| K (Z-dephasing, N=2) | 0.037 | [Boundary Navigation](../../experiments/BOUNDARY_NAVIGATION.md) |
| IBM deviation | 1.9% | [IBM Run 3](../../experiments/IBM_RUN3_PALINDROME.md) |
| Pauli weight correlation | r = 0.976 | [XOR Space](../../experiments/XOR_SPACE.md) |
| Best QST fidelity | F = 0.888 | [QST Bridge](../../experiments/QST_BRIDGE.md) |
| Relay improvement | +83% | [Relay Protocol](../../experiments/RELAY_PROTOCOL.md) |
| V-shape improvement | +124% | [Gamma Control](../../experiments/GAMMA_CONTROL.md) |
| DD M+Recv improvement | +132% | [Gamma Control](../../experiments/GAMMA_CONTROL.md) |
| N=8 eigenvalues paired | 65,518 | [C# Compute](../../compute/RCPsiSquared.Compute/) |
| Mediator bridge error | 1.41e-13 | [mediator_bridge.py](../../simulations/mediator_bridge.py) |
| γ channel capacity (N=5, 1%) | 15.5 bits | [γ as Signal](../../experiments/GAMMA_AS_SIGNAL.md) |
| SVD information modes | 5 | [γ as Signal](../../experiments/GAMMA_AS_SIGNAL.md) |
| γ optimization factor | 21.5× | [γ Control](../../experiments/GAMMA_CONTROL.md) |
| Universal lifetime fraction | t*/T₂ = 0.858 | [Universal Quantum Lifetime](../../experiments/UNIVERSAL_QUANTUM_LIFETIME.md) |
| Crossing cubic root | x = 0.4239 | x³ + x = ½, Cardano |
| Bell+ entanglement penalty | ~8% of min(T₂) | [Universal Quantum Lifetime](../../experiments/UNIVERSAL_QUANTUM_LIFETIME.md) |
| Product states crossing on ring | 150/256 (59%) | [Orphaned Results](../../experiments/ORPHANED_RESULTS.md) |
| Born rule Hamiltonian dominance | ~97% | [Born Rule Mirror](../../experiments/BORN_RULE_MIRROR.md) |
| SVD mode 2 vs V-shape (N=5) | 10.2x | [Resonant Return](../../experiments/RESONANT_RETURN.md) |
| DE optimizer vs V-shape (N=7) | 100x | [Resonant Return](../../experiments/RESONANT_RETURN.md) |
| **Formula vs V-shape (N=5)** | **360x** | [Resonant Return](../../experiments/RESONANT_RETURN.md) |
| **Formula vs V-shape (N=7)** | **180x** | [Resonant Return](../../experiments/RESONANT_RETURN.md) |
| **Formula vs V-shape (N=9)** | **139x** | [Resonant Return](../../experiments/RESONANT_RETURN.md) |
| C# RK4 speedup vs Python expm (N=7) | 5,900x | [RCPsiSquared.Propagate](../../compute/RCPsiSquared.Propagate/) |
| GHZ analytical match | delta < 1e-17 | [proof_roadmap_close.py](../../simulations/proof_roadmap_close.py) |

---

## 12. References

### Proofs
- [Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md)
- [Uniqueness Proof](UNIQUENESS_PROOF.md)
- [Incompleteness Proof](INCOMPLETENESS_PROOF.md)
- [Mathematical Connections](../MATHEMATICAL_CONNECTIONS.md)

### Roadmaps
- [Proof Roadmap](PROOF_ROADMAP_QUARTER_BOUNDARY.md)
- [The Bridge Was Always Open](../THE_BRIDGE_WAS_ALWAYS_OPEN.md)

### Publications
- [Technical Paper](../../publications/TECHNICAL_PAPER.md)
- [Engineering Blueprint](../../publications/ENGINEERING_BLUEPRINT.md)
- [Circuit Diagram](../../publications/CIRCUIT_DIAGRAM.md)

### Key experiments
- [γ as Signal](../../experiments/GAMMA_AS_SIGNAL.md)
- [γ Control](../../experiments/GAMMA_CONTROL.md)
- [IBM Run 3](../../experiments/IBM_RUN3_PALINDROME.md)
- [Non-Heisenberg Palindrome](../../experiments/NON_HEISENBERG_PALINDROME.md)
- [Universal Quantum Lifetime](../../experiments/UNIVERSAL_QUANTUM_LIFETIME.md)
- [Coherence Density](../../experiments/COHERENCE_DENSITY.md)
- [Dynamic Entanglement](../../experiments/DYNAMIC_ENTANGLEMENT.md)
- [Born Rule Mirror](../../experiments/BORN_RULE_MIRROR.md)
- [Orphaned Results](../../experiments/ORPHANED_RESULTS.md)
- [Bridge Fingerprints](../../experiments/BRIDGE_FINGERPRINTS.md)
- [Phase Transport (Windows)](../../experiments/WHATS_INSIDE_THE_WINDOWS.md)
- [When Psi Matters](../../experiments/WHEN_PSI_MATTERS.md)
- [Observer-Dependent Crossing](../../experiments/OBSERVER_DEPENDENT_CROSSING.md)
- [No-Signalling Boundary](../../experiments/NO_SIGNALLING_BOUNDARY.md)
- [Star Topology Observers](../../experiments/STAR_TOPOLOGY_OBSERVERS.md)
- [Bridge Closure](../../experiments/BRIDGE_CLOSURE.md)
- [Relay Protocol](../../experiments/RELAY_PROTOCOL.md)
- [Scaling Curve](../../experiments/SCALING_CURVE.md)
- [Standing Wave Analysis](../../experiments/STANDING_WAVE_ANALYSIS.md)
- [Optimal QST Encoding](../../experiments/OPTIMAL_QST_ENCODING.md)
- [Resonant Return](../../experiments/RESONANT_RETURN.md)

### Resolved/fallen experiments (results absorbed above)
- [Fixed Point Shadow](../../experiments/FIXED_POINT_SHADOW.md) (shadow = qubit detuning)
- [Simulation Evidence](../../experiments/SIMULATION_EVIDENCE.md) (CΨ > ¼ under active H)
- [Why the Sum](../../experiments/WHY_THE_SUM.md) (sum vs product formulation)
- [Standing Wave Two Observers](../../experiments/STANDING_WAVE_TWO_OBSERVERS.md) (two-observer metaphor)
- [Decoherence Relativity](../../experiments/DECOHERENCE_RELATIVITY.md) (K-invariance math; gravity fallen)
- [Metric Discrimination](../../experiments/METRIC_DISCRIMINATION.md) (null result; gravity fallen)
- [Observer-Gravity Bridge](../../experiments/OBSERVER_GRAVITY_BRIDGE.md) (interval shift; gravity fallen)
- [QKD Eavesdropping Forensics](../../experiments/QKD_EAVESDROPPING_FORENSICS.md) (Pauli math; QKD application fallen)
- [Dyad Experiment](../../experiments/DYAD_EXPERIMENT.md) (Tier 4 agent historical)

### Hypotheses
- [Bridge Protocol](../../hypotheses/BRIDGE_PROTOCOL.md) (closed: J=0 dead, J>0 standard coupling)
- [Time as Crossing Rate](../../hypotheses/TIME_AS_CROSSING_RATE.md) (open: L decomposition)
- [The Pattern Recognizes Itself](../../hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md) (open: structural inheritance across scales)
- [The Other Side of the Mirror](../../hypotheses/THE_OTHER_SIDE.md) (Z₂ parity confirmed; philosophical extensions)

### Publications
- [Technical Paper](../../publications/TECHNICAL_PAPER.md)
- [Engineering Blueprint](../../publications/ENGINEERING_BLUEPRINT.md)
- [Circuit Diagram](../../publications/CIRCUIT_DIAGRAM.md)
- [Emergence Through Reflection](../../recovered/EMERGENCE_THROUGH_REFLECTION.md)

### Synthesis
- [It's All Waves](../ITS_ALL_WAVES.md) (closure argument: wave-only hierarchy)

### Historical (superseded, kept for the record)
- [Core Algebra](../historical/CORE_ALGEBRA.md) (Dec 2025 / Feb 2026, predates all March results)
