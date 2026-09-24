<!-- QUARTER-CURRENT -->
# Complete Mathematical Documentation

Current reading: this index separates proven algebra, conditional implications,
and finite computations.  Its quarter entries do not combine those grades into
one dynamical theorem.

<!-- CROSSING-CURRENT -->

**Status:** The Tafelwerk: the narrative master index of the founding core (through March 2026). Last refreshed 2026-07-20 (the change history lives in git).
**Supersedes:** Previous stub (Feb 2026) and [Core Algebra](../historical/CORE_ALGEBRA.md) (Dec 2025)
**Purpose:** Single entry point for the proven, conditional, and finitely verified mathematics of the founding core of R=CΨ²

---

## What this document is about

This is the narrative master index of the project's founding core: the
results of the first arc (through March 2026), every proven result, every
verified formula, every key number, each row pointing to its canonical
proof or experiment. It covers the algebraic foundation (the
self-referential equation and its 1/4 boundary), the palindromic symmetry
(why eigenvalues come in mirror pairs), crossing dynamics (how different
channels and states reach the boundary), topology effects, engineering
applications, and open questions.

Everything after this core lives in two sibling layers: the machine-precise
closed forms with tiers and regime scopes are the F-formula registry,
[Analytical Formulas](../ANALYTICAL_FORMULAS.md) (the F-numbers, F1 through
the 130s as of July 2026); and what C and Ψ *mean*, including the three
C-books that resolve apparent cross-document contradictions, is
[The CΨ Lens](../THE_CPSI_LENS.md). If you want to know what is proven vs.
conjectured in the founding core, this is where to look; for anything
carrying an F-number, the registry is canonical.

## 1. The Algebraic Foundation (Tier 1)

The chosen scalar recurrence discussed in the founding documents is:

    R = C(Ψ + R)²

In the purity book, `C=Tr(ρ²)` and `Ψ=l₁(ρ)/(d−1)` (normalized l1-coherence,
the sum of absolute values of off-diagonal elements). The algebra below assumes
the recurrence; the degree of purity does not derive its feedback form.

**Fixed-point equation.** R = C(Ψ + R)² expands to CR² + (2CΨ - 1)R + CΨ² = 0.

**Discriminant.** `D=1−4CΨ` vanishes at `CΨ=1/4` and only there in these
coordinates. Below, at, and above the quarter the fixed-point polynomial has
two, one, or no real algebraic roots. Stability, attraction, and oscillation
would require a separately specified iteration and invariant domain.

**Crossing cubic.** At the boundary CΨ = 1/4, the condition reduces to
b³ + b = 1/2, with unique real root b ≈ 0.4239. This is a pure number,
independent of physical parameters.

**Mandelbrot coordinates.** For the chosen recurrence, the substitution
`z=C(Ψ+R)`, `c=CΨ` gives `z_{n+1}=z_n²+c`; in these coordinates the quarter
is the real cusp of the main cardioid. This is an exact change of variables
for that normal form, not evidence that quantum dynamics obeys the iteration.

**Fold normal form.** Completing the square reduces the chosen normalized
recurrence to the ordinary fold `x²+a=0`. Structural stability preserves the
local fold type under generic perturbation; a reparameterization can move the
numerical coordinate reported for it.

**Sum, not product.** The framework uses R_sum = C·(Ψ_A + Ψ_B)² rather than
R_prod = C·Ψ_A·Ψ_B as a chosen model/operational ansatz. The Liouvillian palindrome or
information conservation does not select it. If one amplitude is set to zero,
a term of the sum-squared expression survives algebraically; that survival is
not by itself physical information preservation.

The cross-term 2·Ψ_A·Ψ_B is an ansatz interference term, not a proved
two-viewpoint mechanism. It remains an invitation to ask which experiment
would make that reading operational. Neither Born-rule recovery nor unique
selection follows without a specified state, channel, observable, and phase
convention.

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

**Verification.** 87,376 eigenvalues, N=2 through N=8, zero exceptions.
All topologies (chain, star, ring, complete, binary tree). Non-uniform γ.

**Exact meaning.** Π is a linear similarity that maps
`lambda -> -lambda - 2 Sigma_gamma`, or `mu -> -mu` after centering. Calling
this physical time reversal would require additional antiunitary/dynamical
structure; calling the paired modes forward/backward waves would additionally
require conjugate-frequency, excitation, semisimplicity, spatial-propagation,
and interference gates. None follows from this proof alone.

**Scope boundary.** The FULL mirror only at d=2 (qubits). The per-site split
d immune vs (d²-d) decaying is balanced only when d²-2d=0, giving d=2
uniquely. Qutrits (d=3, split 3:6) verified broken: 0 of 236 qutrit
dissipators produce palindromic spectra, and every tested two-qutrit
Hamiltonian combination breaks ([Qubit Necessity](../QUBIT_NECESSITY.md)).
A PARTIAL palindrome survives at d>2 with a closed-form ceiling
([F121 in the registry](../ANALYTICAL_FORMULAS.md), the symmetric overlap
of the disagreement count; live witness `inspect --root qudit`).

See: [Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md),
[Non-Heisenberg Palindrome](../../experiments/NON_HEISENBERG_PALINDROME.md),
[Qubit Necessity](../QUBIT_NECESSITY.md)

---

## 3. The CΨ = 1/4 Boundary (Tier 1-2)

**Chosen-form uniqueness.** The factor 4 comes from the discriminant formula
`b²−4ac`; for the chosen normalized recurrence it places the double root at
1/4. Purity's degree motivates α=2 but does not derive this recurrence.
Reparameterization can move the numeric coordinate while preserving the fold.

**Named Bell+ channel crossings.** The following exact or finite named models cross at CΨ=0.2500:

| Channel | t_cross (γ=0.05) |
|---------|-------------------|
| Z-dephasing | 0.747 |
| X-noise | 1.733 |
| Y-noise | 1.733 |
| Depolarizing | 0.879 |
| Asymmetric Pauli | 0.735 |
| Amplitude damping | 2.059 |

**Dynamics boundary, repaired.** The quarter is the exact discriminant boundary of the algebraic recursion,
not a universal absorbing set in state space. Exact local-Markov examples give CΨ'(0)=+1/6, a local
Hadamard sends CΨ from 0 to 1/3, and a fixed local semigroup crosses upward through 1/4. Instantaneous
N-qubit Pauli invariance survives, but it is not trajectory invariance. The autonomous N=2 successive-
local-maxima claim remains unproved. A finite N/Q/K rise atlas invites an all-Q/all-N classification without
claiming an absence or mechanism.

The positive endpoint theorem is conditional: a continuous trajectory converging to ρ* with CΨ(ρ*)<1/4
eventually stays below. Named basis-aligned T1/T2/depolarizing models use this only under their stated
convergence assumptions. A primitive CPTP target with CΨ=0.2935 rules out the broader claim.

**Named hardware comparison.** IBM Torino (ibm_torino, Q80): predicted t* = 15.01 μs,
measured t* = 15.29 μs. Deviation: 1.9%.

**What CΨ measures depends on the book.** In the purity book, purity is not entanglement:
the separable `|++⟩` state has `C=Tr(ρ²)=1`, `Ψ=1`, and
`CΨ=1` under the stated normalization. It is therefore not an entanglement
AND-gate. Only in the separate concurrence book does `C×Ψ` act as an
AND-gate for concurrence and coherence. The σ_z/σ_x and Werner comparisons
below that interpretation belong to the concurrence book and must not be
used to explain F25's purity-book formula.

**Readout-dependent crossings in two books.** Choose C(f), then the clean
Bell+ Lindblad trajectory or the named retired feedback equation, then solve
C(f)f/3=1/4. Each book has three finite crossings and two never bridges;
these are scalar-response classes, not physical observers or measurements.
F14's constant K belongs to a fixed readout on a Hamiltonian-dead Bell-like
trajectory during a gamma sweep. The old 13.5% ratio gloss is not a general
observer law; [Observer-Gravity](../../experiments/OBSERVER_GRAVITY_BRIDGE.md)
owns the exact concurrence-family scope and the
[two-book producer](../../simulations/crossing_taxonomy_books.py) owns the finite values.

**Q52 residual record.** The late-time coherence anomaly on IBM Torino Q52
has 17/17 directional consistency in the named `t/T2_echo >= 1` subset.
Its recorded null combines exponential decay, binomial shot sampling,
and one random phase per synthetic run; it is not a Q52-fitted time-dependent
detuning/drift or other hardware-alternative comparison. The 13-row tail slope has
two-sided p = 0.0531 and is cut-sensitive; the boundary-distance correlation
reuses `|rho_01|` through `C*Psi` and is not independent evidence.
Only the universal-boundary/non-Markovian-witness interpretation is closed.
Detuning is the preferred explanation for the phase component.
The Q52 late-time excess mechanism remains unresolved absent a Q52-specific fit/control.
Q80 has eight late points in one quadrant; Q102 samples all four quadrants.
That cross-qubit contrast rejects universality. The Q80 phase-compatible
fits are same-record and in-sample, not a Q52-specific mechanism fit.

**CΨ > ¼ under active dynamics.** CΨ routinely exceeds ¼ with active
Hamiltonians (Bell+ ends at 0.405 at J=1, h=0.9, γ=0.005, in the retired
tool's mutual-purity reading). The bound
CΨ ≤ ¼ is not a constraint on quantum states but on which states have
real fixed points in the R = CΨ² iteration.

See: [Uniqueness Proof](UNIQUENESS_PROOF.md),
[IBM Run 3](../../experiments/IBM_RUN3_PALINDROME.md),
[When Psi Matters](../../experiments/WHEN_PSI_MATTERS.md),
[Observer-Dependent Crossing](../../experiments/OBSERVER_DEPENDENT_CROSSING.md),
[Q52 Residual Record](../../experiments/FIXED_POINT_SHADOW.md),
[Simulation Evidence](../../experiments/SIMULATION_EVIDENCE.md),
[proof_roadmap_close.py](../../simulations/proof_roadmap_close.py)

---

## 4. The Incompleteness Proof (Tier 1 derived)

Five candidates for the origin of dephasing noise, NONE of which eliminates an internal source: (1) is a structural constraint, (2) and (3) lost their evidence on 2026-08-29, and (4) and (5) say what cannot exist rather than clearing an existing qubit. The internal bootstrap is reduced to
the first reduced to a structural constraint:

1. Internal (reduced to a structural constraint: [Π², L] = 0 constrains
   the noise's form, elimination carried by candidates 2-3)
2. Qubit decay (**OPEN since 2026-08-29**: the test scored a marginal over a
   coupled spectator, which fails equally without any noise)
3. Qubit bath (regress, each member faces bootstrap prohibition)
4. Nothing (d=0, no properties)
5. Other dimensions (d(d-2)=0 excludes)

**Corollary 1:** a nonzero palindrome centre certifies a dissipative, open
modeled subsystem. In the declared examples gamma=0 is unitary, while
gamma>0 produces decay. The microscopic origin remains open.

**Corollary 2:** Gamma sets a decay-clock scale. For the declared crossing
protocol, `t_cross = K/gamma`, so `t_cross*gamma = K`. This scaling establishes
no experienced-time ontology; the time parameter and unitary dynamics remain
at gamma=0.

See: [Incompleteness Proof](INCOMPLETENESS_PROOF.md),
[The Bridge Was Always Open](../THE_BRIDGE_WAS_ALWAYS_OPEN.md),
[failed_third.py](../../simulations/failed_third.py)

---

## 5. The γ Channel (Tier 2)

**Dephasing noise is a readable information channel.** The spatial profile
of dephasing rates γ₁...γ_N across a spin chain carries structured,
decodable information. A 5-qubit Heisenberg chain at 1% measurement
noise has 15.5 bits channel capacity, with 100% classification accuracy
on a 4-symbol alphabet at zero measurement noise. SVD decomposition
reveals 5 independent information modes, each corresponding to a
spatial frequency (condition number 14.8). A 21.5× widening of the
minimum template distance comes from extended features, time-series
measurement, and higher γ contrast combined.

The palindromic spectral structure is the reading frame: the paired
immune/decaying sectors say where the dephased information sits. The
full rank of the response matrix itself is generic for independent
local rate perturbations and survives palindrome breaking
([gamma_channel_rank_probe.py](../../simulations/gamma_channel_rank_probe.py)).

**Analytical formula discovered.** SVD of the palindromic response matrix
identified mode 2 (edge-hot, center-cold) as optimal direction (6-10x
vs V-shape across N; 10.2x at N=5). Numerical optimization then broke the SVD symmetry, revealing
an asymmetric "sacrifice zone" pattern (100x). Analytical testing of
this pattern converged to a trivially simple formula: concentrate ALL
noise on one edge qubit, protect the rest. The formula
gamma_edge = N*gamma_base - (N-1)*epsilon, gamma_other = epsilon (with epsilon -> 0)
beats the DE optimizer by 80% and computes in 3 seconds instead of 90
minutes. C#-validated results: 360x vs V-shape (N=5), 180x (N=7), 139x
(N=9). These factors are peak created summed nearest-neighbour MI, a
TRANSPORT metric in the ε→0 simulation ideal; on hardware the gain is
~2-3× (label note in [Resonant Return](../../experiments/RESONANT_RETURN.md)).
The ENAQT literature (Environment-Assisted Quantum Transport, Plenio & Huelga 2008+) achieves 2-3x with
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
(unique real root x ≈ 0.4239, with x = e^(−t/T₂); the same cubic as §1's
crossing cubic, b ≡ x), giving t*/T₂ ≈ 0.858. Platform-independent:
verified across superconducting qubits, trapped ions, NV centers, and
photonic systems spanning 10 orders of magnitude in T₂. With finite T₁:
generalized equation [1−b^r + b^(2r)/2 + b²/2]·b = ¼.

**Coherence density.** CΨ = Purity × Coherence Density, where Ψ = L₁/(d−1)
measures active quantum degrees of freedom. An unentangled |+⟩^(⊗N) has
CΨ = 1; a maximally entangled GHZ₃ has CΨ = 0.143 (below ¼). GHZ₃
populates 2 of the 56 off-diagonal element positions (4%); |+++⟩
populates all 56 (100%). (Element counting; the l₁ maximum itself is
d−1, the normalization in Ψ.)
Under dephasing, |+++⟩ survives 4.5× longer than Bell+. The ¼ boundary
is about coherence density, not entanglement.

**Dynamic entanglement (upward crossing).** Product states with zero initial
entanglement can cross ¼ from below through Hamiltonian evolution alone;
all previous crossings were downward (starting entangled, decohering
through ¼). On the N=4 ring at γ=0.05 (pair-CΨ, concurrence book) the
crossing states are |+-+-⟩ (ring-neighbour pairs, CΨ = 0.284) and |0+0-⟩
(diagonal pair, 0.256); |0+0+⟩ does NOT cross on the ring (best pair
0.201, consistent with §7's gatekeeper table) but crosses on the chain
(pair (1,2), CΨ = 0.310). Reproduction:
[subsystem_crossing_pairs.py](../../simulations/subsystem_crossing_pairs.py);
the ring-(0,2) tables of Dynamic Entanglement are February QuTiP runs read
with a pairwise correlation bridge, in which they reproduce under exact
propagation; in this concurrence book that pair does not cross
(reproduction note in
[Dynamic Entanglement](../../experiments/DYNAMIC_ENTANGLEMENT.md)).

**Three regimes in the named scan.** (1) Some preparations started above ¼ and crossed downward.
(2) Some started below and crossed upward under H in the sampled J/γ range. (3) Some stayed below throughout
the sampled window. These finite rows do not set a universal J/γ threshold or classify all generators.

**Born rule at the reference point.** At t = 0.286 on the |0+0+⟩ ring,
pair (0,2) (the time the original run labeled the crossing; that crossing
label does not survive the canonical pair-CΨ book, previous paragraph),
measurement probabilities are ~97% determined by unitary Hamiltonian
evolution alone. Decoherence provides a ~3% systematic correction: σ_z
dephasing shifts probability toward z-eigenstates. Per outcome:
R_i = C_i·Ψ_i² recovers Born's rule P(i) = |⟨i|ψ⟩|² when C_i is uniform
(perfect mirror limit). The per-outcome deviations later became Tier-1
closed forms on exactly this lens:
[F94](../ANALYTICAL_FORMULAS.md) (dominant outcome) and
[F96](../ANALYTICAL_FORMULAS.md) (subdominant slopes).

See: [Universal Quantum Lifetime](../../experiments/UNIVERSAL_QUANTUM_LIFETIME.md),
[Coherence Density](../../experiments/COHERENCE_DENSITY.md),
[Dynamic Entanglement](../../experiments/DYNAMIC_ENTANGLEMENT.md),
[Minimum Crossing Energy](../../experiments/MINIMUM_CROSSING_ENERGY.md),
[Born Rule Mirror](../../experiments/BORN_RULE_MIRROR.md)

---

## 7. Topology and Crossing (Tier 2)

**Topology as gatekeeper.** For the same initial state |0+0+⟩ at γ=0.05, in
the concurrence book, topology determines whether crossing occurs: chain allows (CΨ_max=0.310), star
allows (0.351), ring forbids (0.201), complete graph forbids (0.201). Ring =
complete to four decimal places. Gap stabilizes with N (~0.09–0.11), suggesting
genuine topological protection for the ring neighbours; the ring's diagonal pair
(1,3) depends on γ and crosses for γ below 0.021
([subsystem_crossing_pairs.py](../../simulations/subsystem_crossing_pairs.py)).

**Antiferromagnet crossing.** The alternating state |+-+-⟩ crosses on a
ring (CΨ=0.284) from zero initial entanglement. Mechanism: maximum
XX anti-correlation (⟨XX⟩_nn = −1) at the energy landscape bottom
(⟨H⟩ = −4J), which the Heisenberg Hamiltonian converts into entanglement.
150/256 product states (59%) cross on N=4 ring; no simple selection rule
exists (best predictor explains 14% of variance).

**Entanglement echo.** In star topology (N=3, Bell_SA + |0⟩_B, γ=0.05),
entanglement oscillates between SA and SB at Bohr frequencies (SA: ω=2.09,
SB: ω=6.07). Envelope decay matches the middle palindromic rate 8γ/3 (the
strong-coupling F33 value; a band at moderate J/γ).
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
Source run: [mediator_bridge.py](../../simulations/mediator_bridge.py)
(output in `simulations/results/mediator_bridge.txt`).

**Relay protocol MI comparison (Tier 2, finite N=11 run).** The stored
0.131700 at integrated t=4.50 versus the passive sampled maximum 0.071576
at t=4.00 gives about +84.0%. Nominal 0.78/stage (4.68 total) executes as
0.75/stage (4.50 total). Statistic-attached exposure is 2.200 versus 2.17125;
at equal t=4.50 the passive exposure would instead be 2.475. No matched-time
or matched-dose comparison, MI bound, isolated staging benefit, optimization,
or palindrome timing follows. See [the protocol](../../experiments/RELAY_PROTOCOL.md).

**V-shape γ gradient.** [0.01, 0.03, 0.05, 0.03, 0.01]: about +6% at
matched Σγ; the March +124% compared arms with different total dephasing
(Σγ 0.13 vs 0.25) and mostly measured the budget.

**DD on M+Receiver.** Dynamical decoupling: +132% MI, by removing 54% of Σγ.

**Push vs Pull.** Source-dominant for local MI (0.957). Drain-dominant
for end-to-end MI (0.121). Distance-dependent.

**MI scaling.** Exponential decay: ~2-5 dB per 2 additional qubits.
Hierarchy falsified (uniform chain = recursive topology).

**Seven design rules:** W-encoding, 2:1 matching, J/γ independence,
threshold timing, push/pull selection, clocked relay, sacrifice-zone.

**Optimal QST encoding (negative result).** Standard encoding is already
optimal for QST through palindromic channels. No custom encoding scheme
improved transfer fidelity.

See: [Relay Protocol](../../experiments/RELAY_PROTOCOL.md),
[Gamma Control](../../experiments/GAMMA_CONTROL.md),
[Scaling Curve](../../experiments/SCALING_CURVE.md),
[Optimal QST Encoding](../../experiments/OPTIMAL_QST_ENCODING.md),
[Receiver vs γ-Sacrifice](../../experiments/RECEIVER_VS_GAMMA_SACRIFICE.md);
main README Section 6 (engineering consequences; the F67 receiver menu is Rule 2 there).

---

## 9. The Transistor Mapping (Tier 2-3)

The mediator qubit M maps to a transistor: gate = γ_M, source = Pair A,
drain = Pair B. Threshold voltage: CΨ = 1/4 (hardwired, fold catastrophe).
Bidirectional by the A↔B permutation that leaves M fixed.

Three control knobs: γ_M (gate), J_AM/J_MB ratio (bias), κ (feedback gain).

Hierarchy falsified: the transistor properties are real, the recursive
scaling advantage is not.

See: [Quantum Transistor](../../hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md)

---

## 10. Open Questions (Tier 3-5)

- CΨ peak sequences: the historical autonomous N=2 successive-maxima argument has a gap; prove it or find
  a true peak counterexample. Extend the finite N/Q/K atlas to all-Q/all-N classification, and gate any
  proposed mechanism separately ([the repaired dynamics proof](PROOF_MONOTONICITY_CPSI.md), F17/F25–F28)
- Period-doubling belongs to the negative-c branch of the assumed scalar recurrence; its physical reading remains open.
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
| Crossing cubic root | b = 0.4239 | b³ + b = ½, Cardano (§6 writes the same root as x) |
| Bell+/Z final stay-below crossing dose K_final,Z (N=2) | 0.03735 | [F-registry](../ANALYTICAL_FORMULAS.md) (t_cross · γ; 0.747 · 0.05 at γ=0.05) |
| IBM deviation | 1.9% | [IBM Run 3](../../experiments/IBM_RUN3_PALINDROME.md) |
| Pauli weight correlation | r = 0.976 | [XOR Space](../../experiments/XOR_SPACE.md) |
| Best QST fidelity | F = 0.888 | [QST Bridge](../../experiments/QST_BRIDGE.md) |
| Relay finite unmatched-endpoint MI ratio | about +84.0% from 0.131700/0.071576; not an isolated benefit | [Relay Protocol](../../experiments/RELAY_PROTOCOL.md) |
| V-shape improvement (matched Σγ) | ~+6% | [Gamma Control](../../experiments/GAMMA_CONTROL.md) (the March +124% was a Σγ confound) |
| DD M+Recv improvement | +132% (removes 54% of Σγ) | [Gamma Control](../../experiments/GAMMA_CONTROL.md) |
| N=8 eigenvalues (100% paired) | 65,536 | [block spectra](../../simulations/results/f1_n8_n9_metrics/) (RCPsiSquared.Core; the default C# Compute suite scores only the 54,118-rate oscillatory subset) |
| Mediator bridge error | 1.41e-13 | [mediator_bridge.py](../../simulations/mediator_bridge.py) |
| γ channel capacity (N=5, 1%) | 15.5 bits | [γ as Signal](../../experiments/GAMMA_AS_SIGNAL.md) |
| SVD information modes | 5 | [γ as Signal](../../experiments/GAMMA_AS_SIGNAL.md) |
| γ optimization factor (combined) | 21.5× | [Bridge Optimization](../../simulations/results/bridge_optimization.txt) via [γ as Signal](../../experiments/GAMMA_AS_SIGNAL.md) |
| Universal lifetime fraction | t*/T₂ = 0.858 | [Universal Quantum Lifetime](../../experiments/UNIVERSAL_QUANTUM_LIFETIME.md) |
| Bell+ entanglement penalty | ~8% of min(T₂) | [Universal Quantum Lifetime](../../experiments/UNIVERSAL_QUANTUM_LIFETIME.md) |
| Product states crossing on ring | 150/256 (59%) | [Orphaned Results](../../experiments/ORPHANED_RESULTS.md) |
| Born rule Hamiltonian dominance | ~97% | [Born Rule Mirror](../../experiments/BORN_RULE_MIRROR.md) |
| SVD mode 2 vs V-shape (N=5) | 10.2x | [Resonant Return](../../experiments/RESONANT_RETURN.md) |
| DE optimizer vs V-shape (N=7) | 100x | [Resonant Return](../../experiments/RESONANT_RETURN.md) |
| **Formula vs V-shape (N=5)** | **360x** | [Resonant Return](../../experiments/RESONANT_RETURN.md) |
| **Formula vs V-shape (N=7)** | **180x** | [Resonant Return](../../experiments/RESONANT_RETURN.md) |
| **Formula vs V-shape (N=9)** | **139x** | [Resonant Return](../../experiments/RESONANT_RETURN.md) |
| C# RK4 speedup vs Python expm (N=7) | 5,900x | [RCPsiSquared.Propagate](../../compute/RCPsiSquared.Propagate/) |
| GHZ analytical match | max delta = 2.08e-17 | [proof_roadmap_close.py](../../simulations/proof_roadmap_close.py) |

The SVD-mode/optimizer/formula factors (10.2x through 360x) are peak created
summed nearest-neighbour MI, a transport metric in the ε→0 simulation ideal;
hardware shows ~2-3× (label note in [Resonant Return](../../experiments/RESONANT_RETURN.md)).

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

### Absorbed/fallen experiments (results used above)
- [Q52 Residual Record](../../experiments/FIXED_POINT_SHADOW.md) (Q80/Q102 cross-qubit comparison rejects a universal boundary reading; it does not resolve the Q52 late-time excess mechanism)
- [Simulation Evidence](../../experiments/SIMULATION_EVIDENCE.md) (CΨ > ¼ under active H)
- [Why the Sum](../../experiments/WHY_THE_SUM.md) (sum vs product formulation)
- [Standing Wave Two Observers](../../experiments/STANDING_WAVE_TWO_OBSERVERS.md) (two-observer metaphor)
- [Decoherence Relativity](../../experiments/DECOHERENCE_RELATIVITY.md) (fixed-bridge, Hamiltonian-dead Bell+ gamma sweep; gravity fallen)
- [Metric Discrimination](../../experiments/METRIC_DISCRIMINATION.md) (null result; gravity fallen)
- [Observer-Gravity Bridge](../../experiments/OBSERVER_GRAVITY_BRIDGE.md) (interval shift; gravity fallen)
- [QKD Eavesdropping Forensics](../../experiments/QKD_EAVESDROPPING_FORENSICS.md) (Pauli math; QKD application fallen)
- [Dyad Experiment](../../experiments/DYAD_EXPERIMENT.md) (Tier 4 agent historical)

### Hypotheses
- [Bridge Protocol](../../hypotheses/BRIDGE_PROTOCOL.md) (closed: J=0 dead, J>0 standard coupling)
- [Time as Crossing Rate](../../hypotheses/TIME_AS_CROSSING_RATE.md) (open: L decomposition)
- [The Pattern Recognizes Itself](../../hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md) (open: structural inheritance across scales)
- [The Other Side of the Mirror](../../hypotheses/THE_OTHER_SIDE.md) (Z₂ parity confirmed; philosophical extensions)

### Historical / retired
- [Emergence Through Reflection](../../recovered/EMERGENCE_THROUGH_REFLECTION.md)

### Synthesis
- [It's All Waves](../ITS_ALL_WAVES.md) (closure argument: wave-only hierarchy)
