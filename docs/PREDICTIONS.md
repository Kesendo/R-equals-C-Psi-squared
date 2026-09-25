# Predictions and Interpretive Questions of R = CΨ²

<!-- CROSSING-CURRENT -->

> **Status:** Living reference document. Predictions are individually
> labeled with tier and verification status. The
> [README's hardware table](../README.md) gives the one-page view.

**Date:** 2026-02-09 (created), last refreshed 2026-07-16 (the change history lives in git)
**Depends on:** All experiment documents

**Tier:** Mixed (Tier 1-4, labeled per prediction)
**Status:** Collection with per-item tier labels
**Scope:** All predictions with epistemic tiers and falsification criteria
**Does not establish:** That all predictions are equally established

---

## What this document is about

This is the master catalog of every testable prediction the R=CΨ²
framework has produced, organized by scope and evidence: the
earliest Torino q52 hardware record (§1, with five scoped entries); the
[Confirmations registry](../compute/RCPsiSquared.Core/Confirmations/ConfirmationsRegistry.cs)
is authoritative for registered hardware entries; computationally verified (§2),
the Q52 residual record (§3; interpretation closed, late-time excess mechanism open),
testable with current hardware (§4),
testable in principle (§5), speculative (§6), unverified agent claims
(§7), null results (§8), falsified (§8b), and closed/reopened
hypotheses (§9). Each entry has a tier label, a falsification
criterion, and a link to its source experiment; the counts live in the
Summary by Tier table at the end.

---

## 1. What the Hardware Saw: One Torino Record, Five Scoped Entries

**Scope of this section:** Mixed Tier 1-2: five mixed-scope entries are organized around one 2026-02-09 Torino q52 dataset.
Four entries read or fit the record; the algebraic r → 0 row is context, not a hardware test.
Two registry entries share this raw record. These entries are not five independent hardware tests.

| Prediction | Predicted | Measured | Status | Source |
|------------|-----------|----------|--------|--------|
| C·Ψ = ¼ crossing during free decoherence | Crossing exists | t*/T₂* = 1.04 | **CONFIRMED** | [IBM Quantum Tomography](../experiments/IBM_QUANTUM_TOMOGRAPHY.md) |
| Generalized crossing equation | t*/T₂* = 0.94 (at r = 0.46) | t*/T₂* = 1.04 (11% deviation) | **Finite same-record fitted comparison; not independent confirmation** | [IBM Quantum Tomography](../experiments/IBM_QUANTUM_TOMOGRAPHY.md) |
| T₂* ≠ T₂ for free induction decay | T₂* < T₂ | T₂*/T₂ = 0.37 (factor 2.7×) | **CONFIRMED** | [IBM Quantum Tomography](../experiments/IBM_QUANTUM_TOMOGRAPHY.md) |
| x³ + x = ½ is the r → 0 limit of crossing fraction | 0.858 (pure dephasing) | N/A (algebraic context) | **Algebraic r → 0 limit; not a hardware test** | [Universal Quantum Lifetime](../experiments/UNIVERSAL_QUANTUM_LIFETIME.md) |
| Absorption Theorem ratio Re(λ)/(−2γ⟨n_XY⟩) | = 1 | 1.03 (3%, Q52; two fits to the same N=1 decay record) | **Same-record N=1 fit consistency; N ≥ 2 ladder spacing unmeasured** | [Q52 hardware record](../experiments/IBM_ABSORPTION_THEOREM.md), [Absorption Theorem](ANALYTICAL_FORMULAS.md#at-absorption-theorem-tier-1-proven), [proof](proofs/PROOF_ABSORPTION_THEOREM.md) |

**Hardware:** ibm_torino, T₁ = 221 μs, T₂(echo) = 298 μs, T₂*(FID) = 110 μs.

**Calibration-model context:** The 24,073 historical calibration rows (181 days,
133 qubits) feed the same free-|+⟩ model that classifies them. In this
self-classifying proxy, 2,417 / 24,073 rows are proxy-below; 112 / 133 (84%)
qubit histories contain at least one proxy-below row, with 12 frequent
proxy-below histories. No tomography was performed for these rows. They are
not observed crossing frequencies and are not independent validation.

**The Absorption Theorem** Re(λ) = −2γ⟨n_XY⟩ (the last row above) is Tier-1
**proven** ([the Absorption Theorem proof](proofs/PROOF_ABSORPTION_THEOREM.md)).
It was verified on 1,342 modes (CV = 0). The registered Q52 ratio 1.03 compares
two fits to the same N=1 coherence-decay record.
It is same-record fit consistency, not a ladder measurement.
The N ≥ 2 ladder spacing remains unmeasured. The theorem also falsifies the old
"E = mγ²" guess: the decay law
is **linear** in γ, not quadratic (see [Falsified Predictions](#8b-falsified-predictions)
below).

**The Confirmations registry is the authoritative list of 24 registered hardware entries** (`fw.Confirmations` / [ConfirmationsRegistry.cs](../compute/RCPsiSquared.Core/Confirmations/ConfirmationsRegistry.cs)): ibm_torino + ibm_marrakesh + ibm_kingston, spanning the 2026-02 Torino calibration runs through the 2026-07 Kingston flights: palindrome trichotomy, F25 cusp trajectory, F57 K_dwell γ-invariance, F83/F95, block-CΨ saturation, the F120 moment tower, the F84 heating-leg attribution, the concentrator site contrast, and the F129 standing fringe. Registration does not imply an independent dataset or independent test. The quarter-crossing and absorption entries share the 2026-02-09 Torino q52 raw record, and the five entries in §1 are not five independent confirmations. The registry remains the single live record for registered entries.

---

## 2. Computationally Verified

**Tier: 2, Reproducible via delta_calc MCP tools or standalone Lindblad simulation**

| Prediction | Value | Falsified if | Source |
|------------|-------|-------------|--------|
| Fixed-book Bell+ gamma sweep | Retained 0.039 ± 0.001 feedback fit over 50×; clean Hamiltonian-dead Wootters book K=ln(4/3)/8 | K varies within the same fixed preparation/readout/book | [Gravitational Invariance](../experiments/GRAVITATIONAL_INVARIANCE.md) |
| θ decreases continuously to 0 at C·Ψ = ¼ | Smooth trajectory observed | Discontinuity at boundary | [Boundary Navigation](../experiments/BOUNDARY_NAVIGATION.md) |
| Two real fixed points emerge below ¼ | Topology change confirmed | No bifurcation at ¼ | [Dynamic Fixed Points](../experiments/DYNAMIC_FIXED_POINTS.md) |
| Operator feedback: γ_eff = γ₀(1 − κ⟨O_int⟩) | Modulates γ_eff (~10% at tested params); preservation is parameter-dependent, not a clean separation | Mechanism produces unphysical results | [Operator Feedback](../experiments/OPERATOR_FEEDBACK.md) |
| Ψ_interaction does not shift ¼ boundary | Δδ ≈ −8 × 10⁻⁴ | Boundary shifts under bidirectional coupling | [Core Algebra](historical/CORE_ALGEBRA.md) §8 |
| Readout-dependent quarter equalities | C(f), then evolution book, then C(f)f/3=1/4; clean t=0.5931/0.7192/1.4384 and feedback t=0.6529/0.7735/1.4384 at γ=0.05 | Shared producer fails the named scalar roots | [Observer-Dependent Crossing](../experiments/OBSERVER_DEPENDENT_CROSSING.md) |
| Two never bridges in the finite taxonomy | mutual_purity C=0.5 and overlap C=0.25 stay below the selected quarter level | These fixed bridges cross in either named book | [Observer-Dependent Crossing](../experiments/OBSERVER_DEPENDENT_CROSSING.md) |
| Noise coverage of the taxonomy | One local σ_z sweep, printed by the retired tool under the σ_x and σ_y names too; real σ_x and σ_y keep every bridge's C curve but hold Ψ at 1/3 and move every crossing; depolarizing turns the correlation bridge Type B | A real σ_x run keeps the σ_z crossing times | [Noise Robustness](../experiments/NOISE_ROBUSTNESS.md) |
| N-scaling barrier | Ψ(0) = l₁/(2ᴺ − 1) blocks crossing for GHZ N≥3 and W N≥4 | GHZ N≥3 or W N≥4 cross in this named readout/book | [N-Scaling Barrier](../experiments/N_SCALING_BARRIER.md) |
| W N=3 crosses, GHZ N=3 does not | W: Ψ(0)=0.286 > ¼, GHZ: Ψ(0)=0.143 < ¼ | Both cross or both fail | [N-Scaling Barrier](../experiments/N_SCALING_BARRIER.md) |
| Type A survives at N=3,4 | Correlation C=1.0 for W N=3 (until t≈2.3) and W N=4 (until t≈1.5) | Correlation C drops below 1.0 at larger N | [N-Scaling Barrier](../experiments/N_SCALING_BARRIER.md) |
| Subsystem pairs cross when full system cannot | Bell+⊗Bell+ N=4: pairs (0,1) and (2,3) cross at t=0.080 despite full-system Ψ=0.200 | Pairs fail to cross | [Subsystem Crossing](../experiments/SUBSYSTEM_CROSSING.md) |
| Bell+⊗Bell+ cross-pairs never cross | N=4 ring, γ=0.05, concurrence book: the cross-pairs start at C = 0, l₁ = 0; the Hamiltonian entangles them (concurrence up to 0.56 in t ≤ 5), but CΨ peaks at 0.147 | A cross-pair reaches CΨ ≥ 1/4 (concurrence book) | [Subsystem Crossing](../experiments/SUBSYSTEM_CROSSING.md) + [subsystem_crossing_pairs.py](../simulations/subsystem_crossing_pairs.py) |
| Product state: Ψ=1 but C=0 means no crossing | \|+⟩⊗⁴: every pair has Ψ(0)=1.0 and C=0.000 at all times | \|+⟩⊗⁴ pairs develop nonzero C | [Subsystem Crossing](../experiments/SUBSYSTEM_CROSSING.md) |
| GHZ pair-level coherence is zero | GHZ N=4 traced to any pair: l₁=0.000 at all times | GHZ pairs carry nonzero off-diagonal coherence | [Subsystem Crossing](../experiments/SUBSYSTEM_CROSSING.md) |
| Initial \|+⟩⊗ᴺ is a Heisenberg eigenstate | Zero initial Hamiltonian variance does not imply zero Lindblad dynamics; a dissipator may act | Initial variance nonzero in the named model | [Dynamic Entanglement](../experiments/DYNAMIC_ENTANGLEMENT.md) |
| Product states generate upward crossings | Reproduced (canonical pair-CΨ book, γ=0.05): chain \|0+0+⟩ pair (1,2) 0.310; ring \|+-+-⟩ 0.284, \|0+0-⟩ 0.256 | No product state reaches CΨ ≥ ¼ | [Dynamic Entanglement](../experiments/DYNAMIC_ENTANGLEMENT.md) (reproduction note) + [subsystem_crossing_pairs.py](../simulations/subsystem_crossing_pairs.py) |
| Dephasing kills most dynamic crossings | \|0+0+⟩ ring at γ=0.05: in the concurrence book no pair crosses (best ≈0.20); in the pairwise bridge (P_AB − P_A·P_B)/(1 − P_A·P_B) under exact propagation only the diagonal (0,2) crosses (t ≈ 0.285, max 0.320) | All pairs cross equally under dephasing | [Dynamic Entanglement](../experiments/DYNAMIC_ENTANGLEMENT.md) (reproduction note) |
| Four finite log-coherence traces | Historical slope variations <0.01% on the declared grid, not arbitrary-Hamiltonian log-linearity | Stored finite traces fail reproduction | [Algebraic Exploration](../experiments/ALGEBRAIC_EXPLORATION.md) |
| One memory-feedback log-coherence trace | Historical 24.5% slope variation at κ=0.5, τ=1.0; curvature alone is not a memory certificate | The named finite record fails reproduction | [Core Algebra](historical/CORE_ALGEBRA.md) §11 |
| Coherence-purity bound holds throughout Lindblad trajectory | 0 violations for Bell+ (d=4) and W (d=8) | Trajectory violates C ≥ Ψ²(d-1)/d + 1/d | [Algebraic Exploration](../experiments/ALGEBRAIC_EXPLORATION.md) |
| CΨ after Eve intercept-resend depends on Eve's measurement basis | R(θ_Eve) = [sin²θ + \|sin2θ\|]²/18, closed form | R independent of θ_Eve | [QKD Eavesdropping Forensics](../experiments/QKD_EAVESDROPPING_FORENSICS.md) |
| Concurrence is basis-blind under intercept-resend | Conc = 1−f for all θ_Eve (exact) | Concurrence varies with θ_Eve at fixed f | [QKD Eavesdropping Forensics](../experiments/QKD_EAVESDROPPING_FORENSICS.md) |
| CΨ distinguishes Eve from channel noise at equal Concurrence | 4 causes at Conc=0.80 give CΨ ∈ {0.058, 0.068, 0.115} | CΨ identical for Eve and noise at same Conc | [QKD Eavesdropping Forensics](../experiments/QKD_EAVESDROPPING_FORENSICS.md) |
| Off-diagonal ratio breaks θ_Eve degeneracy | \|ρ₀₁\|/\|ρ₀₃\| = cot(θ_Eve), monotonic | Ratio non-monotonic or θ_Eve-independent | [QKD Eavesdropping Forensics](../experiments/QKD_EAVESDROPPING_FORENSICS.md) |
| All four Bell states give identical CΨ under Eve attack | CΨ(Φ+)=CΨ(Φ−)=CΨ(Ψ+)=CΨ(Ψ−) for all θ_Eve | Any Bell state differs | [QKD Eavesdropping Forensics](../experiments/QKD_EAVESDROPPING_FORENSICS.md) |
| R(θ_Eve) has azimuthal symmetry | R depends on θ_Eve only, not φ | R varies with φ | [QKD Eavesdropping Forensics](../experiments/QKD_EAVESDROPPING_FORENSICS.md) |
| R(θ_Eve) maximum at ~60°, not σ_x | R_max = 0.145 at θ_Eve ≈ 60° | Maximum at 90° (σ_x) | [QKD Eavesdropping Forensics](../experiments/QKD_EAVESDROPPING_FORENSICS.md) |
| MI and Correlation are θ_Eve-independent | ΔMI, ΔCorr identical at all θ_Eve for fixed f | MI or Corr varies with θ_Eve | [QKD Eavesdropping Forensics](../experiments/QKD_EAVESDROPPING_FORENSICS.md) |
| Bridge framework detects Eve at stealth angle | Conc/MI/Corr show 17-24% delta at CΨ stealth | Eve invisible to all metrics at stealth | [QKD Eavesdropping Forensics](../experiments/QKD_EAVESDROPPING_FORENSICS.md) |
| ξ-curvature identical for Eve and noise | d²ξ/dp² difference = 0.0% under added depolarization | Curvature discriminates Eve from noise | [QKD Eavesdropping Forensics](../experiments/QKD_EAVESDROPPING_FORENSICS.md) |
| Off-diagonal ratio invertible under noise | \|ρ₀₁\|/\|ρ₀₃\| is a unique function of θ_Eve at fixed (p,f), not cot(θ) but still invertible | Ratio degenerate or constant | [QKD Eavesdropping Forensics](../experiments/QKD_EAVESDROPPING_FORENSICS.md) |
| θ_stealth deterministic function of (p,f) | θ_stealth(0.10, 0.20) = 72.2°, computable from calibration | θ_stealth unpredictable | [QKD Eavesdropping Forensics](../experiments/QKD_EAVESDROPPING_FORENSICS.md) |
| No-signalling: rho_A unchanged by remote measurement | \|\|Δρ_A\|\| = 0 for Bell+ under B Z-measurement (averaged) | rho_A changes (no-signalling violated) | [No-Signalling Boundary](../experiments/NO_SIGNALLING_BOUNDARY.md) |
| CΨ drops to ¼ under remote measurement | CΨ: 0.500 → 0.250 (C drops 1.0→0.5, Ψ unchanged at 0.5) | CΨ unchanged, or drops to value ≠ ¼ | [No-Signalling Boundary](../experiments/NO_SIGNALLING_BOUNDARY.md) |
| Joint scalar change invisible to the separated local subsystem | At J=0 the uncommunicated trace-preserving B operation leaves A's marginal unchanged | A's marginal changes with B's uncommunicated choice | [No-Signalling Boundary](../experiments/NO_SIGNALLING_BOUNDARY.md) |
| Critical slowing at the cusp: asymptotic K(ε, tol) | (1/2)·ln(4ε/tol) + [−4 + (1/2)·ln(16·tol)]·√ε, zero fit parameters; finite-ε residuals remain, while the correction coefficient agrees at 0.5-2% across the tested tolerance grid | Expansion fails to approach the direct iteration count in its stated scale-separated regime | [Boundary Navigation](../experiments/BOUNDARY_NAVIGATION.md) |
| Cusp dwell time is γ-invariant | K_dwell = γ·t_dwell = 1.080088·δ for Bell+, std < 2×10⁻¹⁷ across γ ∈ [0.1, 10] | K_dwell varies with γ in rescaled units | [Boundary Navigation](../experiments/BOUNDARY_NAVIGATION.md) |

The separate universal-noise claim does not hold: the retired tool's three
channel columns are one σ_z run, and under a real σ_x the correlation
bridge crosses only after its plateau ends (C = 0.75). See
[Noise Robustness](../experiments/NOISE_ROBUSTNESS.md).

---

## 3. The Q52 Residual Record (interpretation closed; Q52 mechanism open)

**Tier: 2. Three finite Q52/Q80/Q102 records are retained.
Only the universal-boundary/non-Markovian-witness interpretation is closed.
Detuning is the preferred explanation for the phase component.
The Q52 late-time excess mechanism remains unresolved absent a Q52-specific fit/control.**

The Q52 record contains several descriptive statistics, each with a narrower
scope than the former shadow interpretation:

| Observation | Current reading | Source |
|-------------|-----------------|--------|
| Late-time excess | Zero exceedances in 10,000 draws of the recorded null described below | [Residual Analysis](../experiments/RESIDUAL_ANALYSIS.md) |
| Directional consistency | 17/17 rows at `t/T2_echo >= 1` in Re+/Im-; the nominal sign probability assumes iid uniform phases and a pre-specified quadrant | [Residual Analysis](../experiments/RESIDUAL_ANALYSIS.md) |
| Finite tail fit | 13-row `>= 1.5` slope `+0.00819/T2_echo`, two-sided p = 0.0531; non-monotone and cut-sensitive | [Q52 Residual Record](../experiments/FIXED_POINT_SHADOW.md) |
| Boundary-distance correlation | `r = -0.9955`, with `\|rho_01\|` reused through `C*Psi`; descriptive, not independent evidence | [Q52 Residual Record](../experiments/FIXED_POINT_SHADOW.md) |
| Algebraic-root phase comparison | The phase of `R-` and the phase of `rho_01` are phases of different objects; no dynamics map connects them here | [Q52 Residual Record](../experiments/FIXED_POINT_SHADOW.md) |

The recorded 10,000-draw null combines exponential decay, binomial shot sampling, and one random phase per synthetic run; it is not a Q52-fitted time-dependent detuning/drift or other hardware-alternative comparison.
The 13-row ≥ 1.5 slope is +0.00819/T2_echo with two-sided p = 0.0531; the amplitudes are non-monotone and the result is cut-sensitive.
Because |ρ₀₁| enters CΨ, r = −0.9955 is a same-record algebraic coupling, not independent boundary evidence.
The algebraic R₋ phase comparison supplies no dynamical mapping to ρ₀₁.

**What the comparison establishes.** [Q52 Residual Record](../experiments/FIXED_POINT_SHADOW.md)
reports exploratory Q80-only fits. The fixed-T2 phase-line/intercept-only errors are 0.0356/0.0508 (1.4x);
the free-complex/envelope-only errors are 0.0138/0.0487 (3.5x). Both Q80 fit
comparisons are in-sample on the same record; neither is a held-out prediction
or a Q52 mechanism fit. Q80 has 8/8 late points in Quadrant 1, while Q102
samples all four quadrants with no common direction. That comparison rejects a
universal boundary reading and motivates the phase-detuning hypothesis; it does
not identify what produces Q52's late-time magnitude excess.

---

## 4. Testable with Current Hardware

**Tier: 3, Concrete protocols, testable on existing quantum hardware or NMR**

> Note: the QKD-forensics rows below verify the *mathematics* (the closed forms reproduce bit-exact), but the eavesdropping-*detection* application they were framed for is retired: the source [QKD Eavesdropping Forensics](../experiments/QKD_EAVESDROPPING_FORENSICS.md) is marked Fallen for the application ("cannot replace or improve standard QKD detection"). The closed forms stand; the forensic protocol claim does not.

| Prediction | Specific value | Test protocol | Falsified if | Source |
|------------|---------------|---------------|-------------|--------|
| Critical slowing at CΨ = ¼ | Diverging convergence period | Tune system toward ¼, measure convergence time | No critical slowing | [Mandelbrot Connection](../experiments/MANDELBROT_CONNECTION.md) |
| Bridge fingerprints: initial state determines crossing trajectory | State-specific C(t), Ψ(t) curves | Prepare different initial states, run tomography through ¼ | All states show identical crossing | [Bridge Fingerprints](../experiments/BRIDGE_FINGERPRINTS.md) |
| Discriminate a specified memory model from Markovian alternatives | The Markovian H=Y control gives ξ″(0)=−7; no iff test follows | Specify alternative generators and compare full trajectories | Diagnostic cannot distinguish its declared alternatives | [Algebraic Exploration](../experiments/ALGEBRAIC_EXPLORATION.md), [Core Algebra](historical/CORE_ALGEBRA.md) §11 |
| QKD eavesdropping forensics | R(θ_E) = [sin²θ+\|sin2θ\|]²/18, ~500 pairs for 3.8σ (noiseless, naive Eve) | Prepare Bell+, intercept-resend on Bob, joint tomography | CΨ carries no θ_E information | [QKD Eavesdropping Forensics](../experiments/QKD_EAVESDROPPING_FORENSICS.md) |
| Multi-metric forensics detects Eve at stealth | MI/Conc/Corr > 17% delta even at CΨ stealth angle | Simultaneous tomography with multiple bridge metrics | All metrics fooled at stealth angle | [QKD Eavesdropping Forensics](../experiments/QKD_EAVESDROPPING_FORENSICS.md) |
| CΨ identifies Eve's measurement basis | R(θ_E) = [sin²θ_E + \|sin 2θ_E\|]²/18 | Joint tomography of Bell+ after intercept-resend | CΨ constant across all θ_E | [QKD Eavesdropping Forensics](../experiments/QKD_EAVESDROPPING_FORENSICS.md) |
| Concurrence = 1−f (basis-independent) | Linear, θ_E-independent | Concurrence measurement on partially intercepted pairs | Concurrence depends on θ_E | [QKD Eavesdropping Forensics](../experiments/QKD_EAVESDROPPING_FORENSICS.md) |
| CΨ azimuthally symmetric under Eve attack | R depends only on θ_E, not φ | Vary φ at fixed θ_E, compute CΨ | CΨ varies with φ | [QKD Eavesdropping Forensics](../experiments/QKD_EAVESDROPPING_FORENSICS.md) |
| Off-diagonal ratio (noiseless = cot(θ_E)) | \|ρ₀₁\|/\|ρ₀₃\| monotonic, breaks degeneracy; noiseless it equals cot(θ_E), under channel noise it deviates but stays invertible (see §2) | Tomographic element extraction | Ratio not monotonic in θ_E | [QKD Eavesdropping Forensics](../experiments/QKD_EAVESDROPPING_FORENSICS.md) |
| Eve σ_z ≈ dephasing in CΨ | CΨ = 0.058 for both at Conc = 0.80 (noiseless) | Compare Eve σ_z with pure dephasing | Distinguishable CΨ values | [QKD Eavesdropping Forensics](../experiments/QKD_EAVESDROPPING_FORENSICS.md) |
| All four Bell states identical under Eve | R(θ_E) same for Φ+, Φ−, Ψ+, Ψ− | Repeat with all Bell states | CΨ depends on initial Bell state | [QKD Eavesdropping Forensics](../experiments/QKD_EAVESDROPPING_FORENSICS.md) |
| Strategic Eve has stealth angle | θ_stealth ≈ 42-74° makes CΨ(Eve+noise) ≈ CΨ(noise) | Optimize θ_E against CΨ detection | No stealth angle exists | [QKD Eavesdropping Forensics](../experiments/QKD_EAVESDROPPING_FORENSICS.md) |
| Moment-tower pump channel (F120, added 2026-06-11; **RUN same day on ibm_kingston**) | slope of ⟨H_p^j⟩ from I/d = (1/d)·Σ_l Δγ_l·t_j(l); first firing rung = girth ℓ ⟹ palindrome hard at m\* = 2ℓ+1; girth-2 witness has slope⟨H⟩ = 0 exactly next to firing slope⟨H²⟩ | Basis-state-averaged preparation, free evolution under the chip's own damping, Pauli-polynomial readout of ⟨H_p^j⟩ early-window slopes | Slopes deviate from the closed form at the device's true pump vector, or the null rung fires | **Structure CONFIRMED** (double null z = 1.47/0.04, row-exact ⟨H²⟩ identity, girth 2 read, site tracking, 0.3-5.7% reproducibility); the in-situ pump ≤ Γ bound holds on all qubits (prep-conditioned, same-day correction of an initial cross-epoch misread); finding: minute-scale T1 telegraphing on q13 and q9; the protocol is self-arbitrating (pump, Γ, γ↑ from one circuit set). [Experiment](../experiments/F120_MOMENT_TOWER_KINGSTON.md), [proof + protocol](proofs/PROOF_MOMENT_TOWER_PUMP_CHANNEL.md), registry entry `f120_moment_tower_kingston_june2026` |

---

## 5. Testable in Principle

**Tier: 3-5, Require hardware or conditions not currently available**

| Prediction | Direction | Would require | Falsified if | Source |
|------------|-----------|---------------|-------------|--------|
| θ measures proximity to ¼ boundary | θ = arctan(√(4CΨ−1)) | Continuous CΨ monitoring near boundary | θ uncorrelated with distance to ¼ | [Boundary Navigation](../experiments/BOUNDARY_NAVIGATION.md) |
| Fractal structure in coherence decay near ¼ | Self-similar patterns | High-resolution time series near boundary | Smooth exponential decay | [Mandelbrot Connection](../experiments/MANDELBROT_CONNECTION.md) |

---

## 6. Speculative

**Tier: 5, No current path to testing; included for completeness**

| Prediction | Implication | Would require | Falsified if | Source |
|------------|------------|---------------|--------------|--------|
| Experienced time = rate of ¼ crossings | High C → more crossings/sec → denser time | Subjective time measurement against coupling strength | Time perception independent of coupling | [Observer-Dependent Crossing](../experiments/OBSERVER_DEPENDENT_CROSSING.md) |
| Anesthesia = C → 0 for environmental coupling | Zero crossings → zero experienced time | Neural coupling measurement during anesthesia | Time perception persists with C = 0 | [Observer-Dependent Crossing](../experiments/OBSERVER_DEPENDENT_CROSSING.md) |

---

## 6a. Interpretive Questions, Not Predictions from F14

<!-- CROSSING-INTERPRETIVE -->

> Interpretive invitation, not a result: these six questions require a new
> physical model. F14's fixed-readout Hamiltonian-dead Bell+ gamma sweep
> supplies neither a spatial metric nor a universal cosmological curve.

| Question | Missing physical link | Source |
|----------|-----------------------|--------|
| Could a horizon correlate with high coherence? | A spatial open-system/gravity model | [Self-Consistency](../recovered/SELF_CONSISTENCY_SCHWARZSCHILD.md) |
| Could voids have a distinctive coherence signature? | A calibrated environment-to-channel map | [Self-Consistency](../recovered/SELF_CONSISTENCY_SCHWARZSCHILD.md) |
| Could a CMB feature relate to a quarter readout? | A quantitative cosmological observable | [Black/White Holes](../recovered/BLACK_WHITE_HOLES_BIGBANG.md) |
| Could an evaporation endpoint show a coherent burst? | A dynamical evaporation model | [Black/White Holes](../recovered/BLACK_WHITE_HOLES_BIGBANG.md) |
| Could a Page curve connect to a scalar recrossing? | An information-recovery model, not a threshold analogy | [Black/White Holes](../recovered/BLACK_WHITE_HOLES_BIGBANG.md) |
| Could black/white-hole imagery describe two directions of a model? | A justified physical map; no universal curve is established | [Black/White Holes](../recovered/BLACK_WHITE_HOLES_BIGBANG.md) |

<!-- CROSSING-CURRENT -->

## 7. Unverified Agent Claims

**Tier: 4, Generated by 120B-parameter local LLM agents using tools no longer available. Could not be independently reproduced (2026-02-08). These are hypotheses, not results.**

See [Mathematical Findings](../experiments/MATHEMATICAL_FINDINGS.md), Sections 8-9, for full context.

| Claim | Claimed value | Why unverified | Falsified if | Source |
|-------|--------------|----------------|-------------|--------|
| C_int ≫ C_ext (33:1 ratio) | 0.427 vs 0.013 | Not reproduced by MCP tools | Ratio ≈ 1 | [Mathematical Findings](../experiments/MATHEMATICAL_FINDINGS.md) §9 |
| t_coh ~ N (linear scaling) | Linear, N = 2 to 6 | Not reproduced by MCP tools | Exponential decay with N | [Mathematical Findings](../experiments/MATHEMATICAL_FINDINGS.md) §9 |
| δ requires dynamics (H ≠ 0) | δ = 0 when H = 0 | Not reproduced by MCP tools | δ > 0 with H = 0 | [Mathematical Findings](../experiments/MATHEMATICAL_FINDINGS.md) §8 |
| C = 0.5 is optimal observer | Maximum R at C = 0.5 | Not independently tested | Peak at C ≠ 0.5 | [Mathematical Findings](../experiments/MATHEMATICAL_FINDINGS.md) |

These claims may be correct. They may also be artifacts of the agent's training data or tool usage. Until independently verified, they carry no epistemic weight.

**Current status** (two of the four are settled):
- **C_int ≫ C_ext (33:1): REFUTED.** [Mathematical Findings](../experiments/MATHEMATICAL_FINDINGS.md) §9: "The 33:1 ratio claimed by the agents does not exist" (proper Lindblad simulation; 21 noise distributions on Bell+ gave identical dynamics).
- **"δ requires dynamics (H≠0)": REFUTED.** [The Genesis of an Oscillation](THE_GENESIS_OF_AN_OSCILLATION.md): the oscillation is J-driven with no threshold (born at Q=0+); §9 records a finite positive-J interval-shift scan; it does not prove an all-J detectability theorem.
- **"C=0.5 is optimal observer": still unverified.** The literal "max R at C=0.5" claim has not been re-tested. Equal normalized E/I fractions give one half arithmetically, but [Complexity Threshold](../hypotheses/COMPLEXITY_THRESHOLD.md) supplies no persistence optimum or neural boundary from that equality.
- **"t_coh ~ N linear": still unverified** (no later test located).

---

## 8. Null Results

**Tier: 2, Computationally verified null result**

| Prediction | Result | Implication | Source |
|------------|--------|-------------|--------|
| Single-system sims discriminate metric forms | **Null:** Cannot distinguish | A γ-form discrimination null; not a framework failure (the earlier "equivalence principle" reading rested on the now-retired gravity interpretation) | [Metric Discrimination](../experiments/METRIC_DISCRIMINATION.md) |

---

## 8b. Falsified Predictions

**Tier: 2, Predictions the framework made and then refuted by its own mathematics, simulations or hardware.**

| Prediction | Why falsified | Correct result | Source |
|------------|---------------|----------------|--------|
| E = mγ² (decay energy quadratic in γ) | The decay law is **linear** in γ, not quadratic | Absorption Theorem: Re(λ) = −2γ⟨n_XY⟩ (linear; verified on 1,342 modes, CV = 0; IBM ratio 1.03) | [the Absorption Theorem proof](proofs/PROOF_ABSORPTION_THEOREM.md) |
| Dephasing survival is basis-dependent: σ_x dephasing moves the surviving \|0+0+⟩ pair from (0,2) to (1,3) | Under σ_x the same single pair (0,2) crosses and (1,3) stays below | N=4 ring, γ=0.05, pairwise bridge under exact propagation: (0,2) max 0.320 under σ_z and 0.335 under σ_x; (1,3) max 0.224 and 0.240 | [Dynamic Entanglement](../experiments/DYNAMIC_ENTANGLEMENT.md) §5.3 + [delta_calc_pairwise_bridge.py](../simulations/delta_calc_pairwise_bridge.py) |
| Cross-pairs stay incoherent: Bell+⊗Bell+ cross-pairs keep C = 0, l₁ = 0 at all times | The Hamiltonian entangles them | N=4 ring, γ=0.05, t ≤ 5: concurrence up to 0.56 and l₁ up to 0.84; they still never cross (concurrence-book CΨ ≤ 0.147) | [Subsystem Crossing](../experiments/SUBSYSTEM_CROSSING.md) §3.3 + [subsystem_crossing_pairs.py](../simulations/subsystem_crossing_pairs.py) |
| Strong dynamics needed for CΨ > ¼ (threshold at h ≈ 0.9) | The sweep's C·Ψ column sat below ¼ because Ψ was fixed at 0.27, not for lack of dynamics | The same runs, read with the density matrix's own Ψ (concurrence × l₁/3), stay above ¼ at h = 0.7 as at h = 1.0 | [Simulation Evidence](../experiments/SIMULATION_EVIDENCE.md) §2 note + [Operator Feedback](../experiments/OPERATOR_FEEDBACK.md) §4 + [delta_calc_feedback_runs.py](../simulations/delta_calc_feedback_runs.py) |

---

## 9. Closed Channel-Free Hypotheses and Separate Coupled Readouts

**Tier: 2, Computationally verified closure for J=0**

| Hypothesis | Result | Why closed | What survives | Source |
|------------|--------|-----------|---------------|--------|
| Bridge protocol (dynamic: B signals A via CΨ crossing) | **Dead for J=0.** No-signalling holds exactly. ρ_A unchanged. CΨ regime change invisible to A. | C is global (ρ_AB), not local. No single-qubit measurement accesses it. | QKD forensics with a channel | [No-Signalling Boundary](../experiments/NO_SIGNALLING_BOUNDARY.md) |
| Bridge protocol (pre-encoded crossing schedule) | No demonstrated advantage at J=0 | Joint fingerprints are not locally available; no-signalling does not equate all entangled correlations with shared randomness | Coupled finite readouts, not a channel-free reopening | [Bridge Closure](../experiments/BRIDGE_CLOSURE.md) |

**Tier: 2, finite inter-qubit J-coupling readings (J > 0), not a channel-free reopening.** The
"environments" in the cited experiment are γ values (0.01-0.50), not
gravitational fields; the experiment's own gravity reading is fallen
(recorded in [Observer × Gravity Bridge](../experiments/OBSERVER_GRAVITY_BRIDGE.md),
which keeps its historical filename). The J-coupling results below stand.

| Prediction | Result | Source |
|-----------|--------|--------|
| Fixed-book Bell+ concurrence crossing | Hamiltonian-dead Bell+ with Wootters concurrence and equal local Z-dephasing: K_conc=ln(4/3)/8=0.0359602590564726. The fixed book includes generator, state, readout, target and crossing convention; no arbitrary-H/readout/channel factorization follows. A fixed-J gamma sweep generally moves Q=J/γ. See the [F14 scope](ANALYTICAL_FORMULAS.md#f14-k-invariance-tier-2-lindblad-scaling) | [Observer × Gravity Bridge](../experiments/OBSERVER_GRAVITY_BRIDGE.md) |
| K is state-dependent (varies with initial entanglement) | Confirmed, in closed form: K_conc(α) = ln(4 sin²(2α)/3)/8 on cosα\|00⟩ + sinα\|11⟩, initial equality at α = 30°, not a positive-time downward crossing | [Observer × Gravity Bridge](../experiments/OBSERVER_GRAVITY_BRIDGE.md) |
| Interval shift in the sampled positive-J range | No threshold seen in this finite scan; J = 0.001 gives about 0.03% shift. Detectability needs an acquisition/error model | [Observer × Gravity Bridge](../experiments/OBSERVER_GRAVITY_BRIDGE.md) |
| Different local readouts in the named run | \|++⟩ reaches the selected local quarter level; Bell+ does not in that run; neither is an experienced clock | [Observer × Gravity Bridge](../experiments/OBSERVER_GRAVITY_BRIDGE.md) |
| Concurrence-book α < 30° has no positive downward crossing | The Hamiltonian-dead cosα\|00⟩+sinα\|11⟩ family starts below the selected level; no claim about time | [Observer × Gravity Bridge](../experiments/OBSERVER_GRAVITY_BRIDGE.md) |
| Finite coupled-model crossing-time shift | Δt = -0.218 at J=0.01, a 2.54% shift against 8.5837; this is not a channel-capacity measurement | [Observer × Gravity Bridge](../experiments/OBSERVER_GRAVITY_BRIDGE.md) |
| Toy repeated-trajectory jitter estimate | N_min=(σ/Δt)² gives about 21 at σ=1.0 (11.7% of t₀), or ~1550 at full-t₀ jitter; no acquisition/error criterion or one-bit channel is established | [Observer × Gravity Bridge](../experiments/OBSERVER_GRAVITY_BRIDGE.md) |
| Coupling accelerates local crossing vs single qubit (0.07x at J=1) | Confirmed. \|+,0⟩ with J=1: t=0.64 vs single \|+⟩: t=8.58 | [Observer × Gravity Bridge](../experiments/OBSERVER_GRAVITY_BRIDGE.md) |
| No energy threshold for crossing: same ⟨H⟩, different crossing | Confirmed. ⟨H⟩ = J for all α, crossing depends on CΨ(0) | [Minimum Crossing Energy](../experiments/MINIMUM_CROSSING_ENERGY.md) |
| α_critical = 30° exactly, CΨ(0) = ¼ exactly | Confirmed. Binary search to 10⁻¹⁵ precision | [Minimum Crossing Energy](../experiments/MINIMUM_CROSSING_ENERGY.md) |
| Product-state pumping in the finite run | \|01> and \|10> reach 0.308/YES in the current producer at J/γ=20 (historical table: 0.309). The separate \|+,0>, \|0,+>, \|+,1> rows are unresolved: current producer 0.077/NO versus historical 0.295/YES | [Minimum Crossing Energy](../experiments/MINIMUM_CROSSING_ENERGY.md) |
| Three named initial eigenstates: finite CΨ_max=0 | \|++⟩, \|00⟩, \|11⟩ in the recorded concurrence run; \|++⟩ dephases and its purity changes, while \|00⟩ and \|11⟩ are stationary | [Minimum Crossing Energy](../experiments/MINIMUM_CROSSING_ENERGY.md) |
| Finite \|01> product J/γ sweep | The stated scan crosses between J/γ=5 (0.248/NO) and 10 (0.286/YES), not a state-universal threshold | [Minimum Crossing Energy](../experiments/MINIMUM_CROSSING_ENERGY.md) |

---

## Summary by Tier

| Tier | Count | Examples |
|------|-------|---------|
| **Registered hardware entries** | **24** in the [Confirmations registry](../compute/RCPsiSquared.Core/Confirmations/ConfirmationsRegistry.cs) (ibm_torino + Marrakesh/Kingston, 2026-02 to 2026-07); §1 has five mixed-scope entries around one 2026-02-09 Torino q52 record, not five records; the quarter-crossing and absorption registrations share that dataset | §1 distinguishes qualitative/direct observations (¼ crossing and T₂*/T₂), the generalized crossing equation (same-record fitted comparison; no independent hardware test), algebraic context that is not a hardware test (x³+x=½), and absorption same-record N=1 fit consistency; N ≥ 2 ladder spacing remains unmeasured. Registry: the 3 Torino runs, palindrome trichotomy, F25 cusp, F57 K_dwell, F83/F95, block-CΨ saturation, F120 moment tower, F84 heating leg, concentrator site contrast, F129 standing fringe |
| **Proven (analytical)** | 1 | Absorption Theorem: Re(λ) = −2γ⟨n_XY⟩, verified on 1,342 modes, CV=0 |
| **Computational records** | 51 (38 rows in §2 + 13 finite J>0 rows in §9, with local scope and unresolved labels) | From γ·t_cross invariance, the θ trajectory, and the N-scaling barrier through the QKD closed-form family and no-signalling to the J>0 clock results; the rows themselves are the list |
| **Q52 residual record** | 1 finite record | Interpretation closed only at the universal-boundary/non-Markovian-witness level; detuning is preferred for the phase component; the Q52 late-time excess mechanism remains open pending a Q52-specific fit/control (§3) |
| **Testable now** | 13 rows in §4 | Critical slowing, fingerprints, specified memory-model discrimination, the QKD forensics family (math verified; application retired), stealth angle existence, F120 moment tower (already flown and registered) |
| **Testable in principle** | 2 | θ compass, fractal decay |
| **Speculative** | 2 | Experienced time as crossing rate, anesthesia as C → 0 |
| **Interpretive questions** | 6 | Horizon coherence, voids, CMB, burst, Page curve, black/white-hole imagery; no F14 predictions |
| **Unverified agent claims** | 4 | 33:1 ratio, linear scaling, H≠0 requirement, optimal C |
| **Null result** | 1 | Metric discrimination |
| **Closed hypothesis (J=0)** | 2 | Bridge dynamic (no-signalling), Bridge pre-encoded (no demonstrated advantage; not an equivalence of all entangled and classical correlations). Separate J>0 inter-qubit readings, not a channel-free reopening. |
| **Unestablished noise generalization** | 1 | Universal five-bridge noise-independence remains unsupported by the retained coverage |
| **Falsified predictions** | 4 | E=mγ² (wrong: α=2γ⟨n_XY⟩, linear not quadratic); σ_x dephasing moving the surviving \|0+0+⟩ pair to (1,3) (σ_x keeps (0,2)); Bell+⊗Bell+ cross-pairs staying incoherent (the Hamiltonian entangles them); the h ≈ 0.9 field threshold (the sweep's column sat below ¼ through a fixed Ψ) |

---

*This document consolidates predictions from across the R = CΨ² framework.*
*For the proven algebra, see [Core Algebra](historical/CORE_ALGEBRA.md). For the interpretive framework, see [Interpretive Framework](../hypotheses/archive/INTERPRETIVE_FRAMEWORK.md).*
*For the phase boundary analysis, see [Dynamic Fixed Points](../experiments/DYNAMIC_FIXED_POINTS.md).*
*For weaknesses and honest self-assessment, see [Weaknesses](WEAKNESSES_OPEN_QUESTIONS.md).*
