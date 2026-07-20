# Testable Predictions of R = CΨ²

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
framework has produced, organized by verification status: confirmed on
hardware (§1, with the [Confirmations registry](../compute/RCPsiSquared.Core/Confirmations/ConfirmationsRegistry.cs)
as the authoritative live record), computationally verified (§2), the
resolved shadow anomaly (§3), testable with current hardware (§4),
testable in principle (§5), speculative (§6), unverified agent claims
(§7), null results (§8), falsified (§8b), and closed/reopened
hypotheses (§9). Each entry has a tier label, a falsification
criterion, and a link to its source experiment; the counts live in the
Summary by Tier table at the end.

---

## 1. Empirically Verified on Real Quantum Hardware

**Tier: 2, Measured on IBM Torino (Heron r1), qubit 52, 2026-02-09**

| Prediction | Predicted | Measured | Status | Source |
|------------|-----------|----------|--------|--------|
| C·Ψ = ¼ crossing during free decoherence | Crossing exists | t*/T₂* = 1.04 | **CONFIRMED** | [IBM Quantum Tomography](../experiments/IBM_QUANTUM_TOMOGRAPHY.md) |
| Generalized crossing equation | t*/T₂* = 0.94 (at r = 0.46) | t*/T₂* = 1.04 (11% deviation) | **Partially confirmed** | [IBM Quantum Tomography](../experiments/IBM_QUANTUM_TOMOGRAPHY.md) |
| T₂* ≠ T₂ for free induction decay | T₂* < T₂ | T₂*/T₂ = 0.37 (factor 2.7×) | **CONFIRMED** | [IBM Quantum Tomography](../experiments/IBM_QUANTUM_TOMOGRAPHY.md) |
| x³ + x = ½ is the r → 0 limit of crossing fraction | 0.858 (pure dephasing) | Consistent | **Algebraically confirmed** | [Universal Quantum Lifetime](../experiments/UNIVERSAL_QUANTUM_LIFETIME.md) |
| Absorption Theorem ratio Re(λ)/(−2γ⟨n_XY⟩) | = 1 | 1.03 (3%, Q52; detuning oscillations, not cavity fringes) | **CONFIRMED** | [Absorption Theorem](ANALYTICAL_FORMULAS.md#at-absorption-theorem-tier-1-proven), [proof](proofs/PROOF_ABSORPTION_THEOREM.md) |

**Hardware:** ibm_torino, T₁ = 221 μs, T₂(echo) = 298 μs, T₂*(FID) = 110 μs.

**Supporting evidence:** 24,073 historical calibration records (181 days, 133 qubits) validate the theory curve C_min(r). 10.1% of snapshots below crossing threshold, 84% of qubits cross at least once, 12 qubits cross almost every day.

**The Absorption Theorem** Re(λ) = −2γ⟨n_XY⟩ (the last row above) is Tier-1 **proven** ([the Absorption Theorem proof](proofs/PROOF_ABSORPTION_THEOREM.md)) and verified on 1,342 modes (CV = 0); the 1.03 ratio is its IBM confirmation. It also falsifies the old "E = mγ²" guess: the decay law is **linear** in γ, not quadratic (see [Falsified Predictions](#8b-falsified-predictions) below).

**The authoritative list of hardware-confirmed predictions is the Confirmations registry** (`fw.Confirmations` / [ConfirmationsRegistry.cs](../compute/RCPsiSquared.Core/Confirmations/ConfirmationsRegistry.cs)): 24 entries (ibm_torino + ibm_marrakesh + ibm_kingston, spanning the 2026-02 Torino calibration runs through the 2026-07 Kingston flights: palindrome trichotomy, F25 cusp trajectory, F57 K_dwell γ-invariance, F83/F95, block-CΨ saturation, the F120 moment tower, the F84 heating-leg attribution, the concentrator site contrast, the F129 standing fringe). Section 1 here details the earliest Torino set, itself registered; the registry is the single live record.

---

## 2. Computationally Verified

**Tier: 2, Reproducible via delta_calc MCP tools or standalone Lindblad simulation**

| Prediction | Value | Falsified if | Source |
|------------|-------|-------------|--------|
| γ · t_cross = constant across decoherence rates | 0.039 ± 0.001 over 50× range | Product varies with γ | [Gravitational Invariance](../experiments/GRAVITATIONAL_INVARIANCE.md) |
| θ decreases continuously to 0 at C·Ψ = ¼ | Smooth trajectory observed | Discontinuity at boundary | [Boundary Navigation](../experiments/BOUNDARY_NAVIGATION.md) |
| Two real fixed points emerge below ¼ | Topology change confirmed | No bifurcation at ¼ | [Dynamic Fixed Points](../experiments/DYNAMIC_FIXED_POINTS.md) |
| Operator feedback: γ_eff = γ₀(1 − κ⟨O_int⟩) | Modulates γ_eff (~10% at tested params); preservation is parameter-dependent, not a clean separation | Mechanism produces unphysical results | [Operator Feedback](../experiments/OPERATOR_FEEDBACK.md) |
| Ψ_interaction does not shift ¼ boundary | Δδ ≈ −8 × 10⁻⁴ | Boundary shifts under bidirectional coupling | [Core Algebra](historical/CORE_ALGEBRA.md) §8 |
| Observer-dependent crossing time | t_cross = 0.652 / 0.773 / 1.437 for mutual_info / concurrence / correlation | All bridge types give same t_cross | [Observer-Dependent Crossing](../experiments/OBSERVER_DEPENDENT_CROSSING.md) |
| Two observers never see crossing | mutual_purity (C=0.5), overlap (C=0.25): C·Ψ < ¼ always | These observers eventually cross | [Observer-Dependent Crossing](../experiments/OBSERVER_DEPENDENT_CROSSING.md) |
| Crossing taxonomy noise-robust | Type A/B/C identical under σ_x, σ_y, σ_z | Taxonomy changes with jump operator | [Noise Robustness](../experiments/NOISE_ROBUSTNESS.md) |
| N-scaling barrier | Ψ(0) = l1/(2^N−1) blocks crossing for GHZ N≥3 and W N≥4 | Standard states cross at large N | [N-Scaling Barrier](../experiments/N_SCALING_BARRIER.md) |
| W N=3 crosses, GHZ N=3 does not | W: Ψ(0)=0.286 > ¼, GHZ: Ψ(0)=0.143 < ¼ | Both cross or both fail | [N-Scaling Barrier](../experiments/N_SCALING_BARRIER.md) |
| Type A survives at N=3,4 | Correlation C=1.0 for W N=3 (until t≈2.3) and W N=4 (until t≈1.5) | Correlation C drops below 1.0 at larger N | [N-Scaling Barrier](../experiments/N_SCALING_BARRIER.md) |
| Subsystem pairs cross when full system cannot | Bell+xBell+ N=4: pairs (0,1) and (2,3) cross at t=0.080 despite full-system Psi=0.200 | Pairs fail to cross | [Subsystem Crossing](../experiments/SUBSYSTEM_CROSSING.md) |
| Non-entangled pairs never cross | Bell+xBell+ cross-pairs (0,2) etc: C=0, l1=0 at all times | Cross-pairs eventually develop coherence | [Subsystem Crossing](../experiments/SUBSYSTEM_CROSSING.md) |
| Product state: Psi=1 but C=0 means no crossing | \|+⟩^4: every pair has Psi=1.0 and C=0.000 permanently | Product state pairs develop nonzero C | [Subsystem Crossing](../experiments/SUBSYSTEM_CROSSING.md) |
| GHZ pair-level coherence is zero | GHZ N=4 traced to any pair: l1=0.000 at all times | GHZ pairs carry nonzero off-diagonal coherence | [Subsystem Crossing](../experiments/SUBSYSTEM_CROSSING.md) |
| \|+⟩^N is Heisenberg eigenstate | Energy variance = 0, no dynamics, C = 0 forever | \|+⟩^N evolves nontrivially | [Dynamic Entanglement](../experiments/DYNAMIC_ENTANGLEMENT.md) |
| Product state \|0+0+⟩ generates crossings | All 6 pairs cross under unitary Heisenberg evolution | No pair reaches CΨ >= 1/4 | [Dynamic Entanglement](../experiments/DYNAMIC_ENTANGLEMENT.md) |
| Dephasing kills most dynamic crossings | Only pair (0,2) crosses with gamma=0.05; others reach max 0.247 | All pairs cross equally under dephasing | [Dynamic Entanglement](../experiments/DYNAMIC_ENTANGLEMENT.md) |
| Dephasing survival is basis-dependent | σ_z dephasing spares \|0⟩-qubits (0,2), kills \|+⟩-qubits (1,3) | Crossing pattern independent of noise basis | [Dynamic Entanglement](../experiments/DYNAMIC_ENTANGLEMENT.md) |
| ξ = ln(Ψ) linear under Markovian dephasing | Slope variation < 0.01% across all tested configs | Variation exceeds 1% for any Markovian channel | [Algebraic Exploration](../experiments/ALGEBRAIC_EXPLORATION.md) |
| ξ linearity breaks under non-Markovian noise | Slope variation 24.5% under memory-kernel feedback (κ=0.5, τ=1.0), vs < 0.01% Markovian | ξ stays linear under memory kernel feedback | [Core Algebra](historical/CORE_ALGEBRA.md) §11 |
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
| CΨ regime change invisible to local subsystem | A cannot detect the quantum→boundary transition | A detects the regime change (no-signalling violated) | [No-Signalling Boundary](../experiments/NO_SIGNALLING_BOUNDARY.md) |
| Critical slowing at the cusp: closed-form K(ε, tol) | (1/2)·ln(4ε/tol) + [−4 + (1/2)·ln(16·tol)]·√ε, zero fit parameters, matches 0.5-2% across 10 ε-decades × 5 tol-decades | Formula fails at any tested scale | [Boundary Navigation](../experiments/BOUNDARY_NAVIGATION.md) |
| Cusp dwell time is γ-invariant | K_dwell = γ·t_dwell = 1.080088·δ for Bell+, std < 2×10⁻¹⁷ across γ ∈ [0.1, 10] | K_dwell varies with γ in rescaled units | [Boundary Navigation](../experiments/BOUNDARY_NAVIGATION.md) |

---

## 3. The Shadow Anomaly (resolved: qubit-specific detuning)

**Tier: 2. The anomaly is real; its cause is settled: qubit-specific
frequency detuning (an uncompensated Z-type offset), not a boundary
property.**

Detected in IBM Torino tomography data (2026-02-09); statistically
significant on every signature measured:

| Observation | Value | Source |
|-------------|-------|--------|
| Late-time excess coherence | p < 0.0001 (10,000-run Monte Carlo null) | [Residual Analysis](../experiments/RESIDUAL_ANALYSIS.md) |
| Directional consistency Re⁺/Im⁻ | 17/17 measurements (P = 6 × 10⁻¹¹) | [Residual Analysis](../experiments/RESIDUAL_ANALYSIS.md) |
| Rising coherence trend in classical regime | +0.008/T₂ slope | [Residual Analysis](../experiments/RESIDUAL_ANALYSIS.md) |
| Boundary correlation | r = −0.9955 (\|ρ₀₁\| vs distance from ¼) | [Fixed Point Shadow](../experiments/FIXED_POINT_SHADOW.md) |
| Shadow direction on Q52 | FP⁻ phase = −12°, residual phase = −48° (same quadrant) | [Fixed Point Shadow](../experiments/FIXED_POINT_SHADOW.md) |

**The cause** ([Fixed Point Shadow](../experiments/FIXED_POINT_SHADOW.md),
retrodiction): an explicit detuning model beats standard Lindblad ~4× on
Q80. The discriminating data are the shadow directions, which are
qubit-specific, not universal: Q52 sits in Quadrant 4 at −44°, Q80 shows
8/8 points in Quadrant 1 at +29° ± 10°, Q102 shows no consistent
direction (0/8). A universal boundary structure would put every qubit in
the same quadrant, so that reading is not supported; a fixed SPAM offset
and a TLS reading are both superseded by the explicit detuning model
(which also carries the Occam weight). The excess coherence is a real
device effect with a known, per-qubit cause.

---

## 4. Testable with Current Hardware

**Tier: 3, Concrete protocols, testable on existing quantum hardware or NMR**

> Note: the QKD-forensics rows below verify the *mathematics* (the closed forms reproduce bit-exact), but the eavesdropping-*detection* application they were framed for is retired: the source [QKD Eavesdropping Forensics](../experiments/QKD_EAVESDROPPING_FORENSICS.md) is marked Fallen for the application ("cannot replace or improve standard QKD detection"). The closed forms stand; the forensic protocol claim does not.

| Prediction | Specific value | Test protocol | Falsified if | Source |
|------------|---------------|---------------|-------------|--------|
| Strong dynamics needed for CΨ > ¼ | Threshold at h ≈ 0.9 | Parameter sweep of transverse field strength | CΨ > ¼ at low h | [Simulation Evidence](../experiments/SIMULATION_EVIDENCE.md) |
| Critical slowing at CΨ = ¼ | Diverging convergence period | Tune system toward ¼, measure convergence time | No critical slowing | [Mandelbrot Connection](../experiments/MANDELBROT_CONNECTION.md) |
| Bridge fingerprints: initial state determines crossing trajectory | State-specific C(t), Ψ(t) curves | Prepare different initial states, run tomography through ¼ | All states show identical crossing | [Bridge Fingerprints](../experiments/BRIDGE_FINGERPRINTS.md) |
| ξ curvature detects non-Markovian noise | ξ = ln(Ψ) linear iff Markovian; curves iff memory effects present | Measure Ψ(t) via tomography, compute ξ(t), fit linearity | ξ is linear on hardware known to have non-Markovian signatures | [Algebraic Exploration](../experiments/ALGEBRAIC_EXPLORATION.md) (linearity), [Core Algebra](historical/CORE_ALGEBRA.md) §11 (curvature) |
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
| Event horizon = maximum coherence (τ = 0) | Coherence peaks near horizon | Analog black hole experiments | Coherence minimum at horizon | [Self-Consistency Schwarzschild](../recovered/SELF_CONSISTENCY_SCHWARZSCHILD.md) |
| Intergalactic voids = most quantum regions | Low γ → high CΨ | Space-based quantum correlation measurements | No gravitational environment dependence | [Self-Consistency Schwarzschild](../recovered/SELF_CONSISTENCY_SCHWARZSCHILD.md) |
| Fractal structure in coherence decay near ¼ | Self-similar patterns | High-resolution time series near boundary | Smooth exponential decay | [Mandelbrot Connection](../experiments/MANDELBROT_CONNECTION.md) |

---

## 6. Speculative

**Tier: 5, No current path to testing; included for completeness**

| Prediction | Implication | Would require | Falsified if | Source |
|------------|------------|---------------|--------------|--------|
| CMB = universal CΨ = ¼ crossing | Big Bang as phase transition from complex to classical | Cosmological extension of framework |  | [Black/White Holes](../recovered/BLACK_WHITE_HOLES_BIGBANG.md) |
| Black hole evaporation ends with coherent burst | Final state not thermal but shows coherence | Observation of BH end-states |  | [Black/White Holes](../recovered/BLACK_WHITE_HOLES_BIGBANG.md) |
| Page time = re-crossing ¼ from below | Information recovery begins at re-entry to complex regime | Quantitative Page curve model |  | [Black/White Holes](../recovered/BLACK_WHITE_HOLES_BIGBANG.md) |
| Black/white holes = opposite directions on universal curve | τ → 0 from both sides | Resolution of information paradox |  | [Black/White Holes](../recovered/BLACK_WHITE_HOLES_BIGBANG.md) |
| Experienced time = rate of ¼ crossings | High C → more crossings/sec → denser time | Subjective time measurement against coupling strength | Time perception independent of coupling | [Observer-Dependent Crossing](../experiments/OBSERVER_DEPENDENT_CROSSING.md) |
| Anesthesia = C → 0 for environmental coupling | Zero crossings → zero experienced time | Neural coupling measurement during anesthesia | Time perception persists with C = 0 | [Observer-Dependent Crossing](../experiments/OBSERVER_DEPENDENT_CROSSING.md) |

---

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

**Current status** (three of the four are settled):
- **C_int ≫ C_ext (33:1): REFUTED.** [Mathematical Findings](../experiments/MATHEMATICAL_FINDINGS.md) §9: "The 33:1 ratio claimed by the agents does not exist" (proper Lindblad simulation; 21 noise distributions on Bell+ gave identical dynamics).
- **"δ requires dynamics (H≠0)": REFUTED.** [The Genesis of an Oscillation](THE_GENESIS_OF_AN_OSCILLATION.md): the oscillation is J-driven with no threshold (born at Q=0+); §9 of this doc already carries the "interval shift continuous in J, no threshold" result that undercuts it.
- **"C=0.5 is optimal observer": now structurally confirmed in spirit.** C=0.5 (half-occupation) is the framework's universal building-block ratio (V-Effect; [Complexity Threshold](../hypotheses/COMPLEXITY_THRESHOLD.md): "population balance C=0.5 with broken magnitude symmetry"). The literal "max R at C=0.5" calculus claim was not separately re-tested.
- **"t_coh ~ N linear": still unverified** (no later test located).

---

## 8. Null Results

**Tier: 2, Computationally verified null result**

| Prediction | Result | Implication | Source |
|------------|--------|-------------|--------|
| Single-system sims discriminate metric forms | **Null:** Cannot distinguish | A γ-form discrimination null; not a framework failure (the earlier "equivalence principle" reading rested on the now-retired gravity interpretation) | [Metric Discrimination](../experiments/METRIC_DISCRIMINATION.md) |

---

## 8b. Falsified Predictions

**Tier: 2, Predictions the framework made and then refuted by its own mathematics or hardware.**

| Prediction | Why falsified | Correct result | Source |
|------------|---------------|----------------|--------|
| Crossing taxonomy is noise-dependent | It is noise-**independent**: Type A/B/C identical under σ_x, σ_y, σ_z | See §2 "crossing taxonomy noise-robust" | [Noise Robustness](../experiments/NOISE_ROBUSTNESS.md) |
| E = mγ² (decay energy quadratic in γ) | The decay law is **linear** in γ, not quadratic | Absorption Theorem: Re(λ) = −2γ⟨n_XY⟩ (linear; verified on 1,342 modes, CV = 0; IBM ratio 1.03) | [the Absorption Theorem proof](proofs/PROOF_ABSORPTION_THEOREM.md) |

---

## 9. Closed Hypotheses (J=0) and Reopened (J>0)

**Tier: 2, Computationally verified closure for J=0**

| Hypothesis | Result | Why closed | What survives | Source |
|------------|--------|-----------|---------------|--------|
| Bridge protocol (dynamic: B signals A via CΨ crossing) | **Dead for J=0.** No-signalling holds exactly. ρ_A unchanged. CΨ regime change invisible to A. | C is global (ρ_AB), not local. No single-qubit measurement accesses it. | QKD forensics with a channel | [No-Signalling Boundary](../experiments/NO_SIGNALLING_BOUNDARY.md) |
| Bridge protocol (pre-encoded: CΨ fingerprints > classical keys) | **Dead for J=0.** A's info ⊆ {ρ_A(0), E_A}. Entanglement without a channel = shared randomness. | Fingerprints require ρ_AB. Qubit carries less info than schedule. | Interval shift for J>0 | [Bridge Closure](../experiments/BRIDGE_CLOSURE.md) |

**Tier: 2, Reopened via inter-qubit J-coupling (J > 0).** The
"environments" in the cited experiment are γ values (0.01-0.50), not
gravitational fields; the experiment's own gravity reading is fallen
(recorded in [Observer × Gravity Bridge](../experiments/OBSERVER_GRAVITY_BRIDGE.md),
which keeps its historical filename). The J-coupling results below stand.

| Prediction | Result | Source |
|-----------|--------|--------|
| t_cross = K(Observer, State) / γ, K is γ-invariant | Confirmed. CV = 0.00% across 6 γ values | [Observer × Gravity Bridge](../experiments/OBSERVER_GRAVITY_BRIDGE.md) |
| K(Conc)/K(MI) = 1.2125 across all γ values | Confirmed. CV = 0.00% | [Observer × Gravity Bridge](../experiments/OBSERVER_GRAVITY_BRIDGE.md) |
| K ratio is state-dependent (varies with initial entanglement) | Confirmed. CV = 13.5% across α = 35°-45° | [Observer × Gravity Bridge](../experiments/OBSERVER_GRAVITY_BRIDGE.md) |
| Interval shift Δt continuous in J, no threshold | Confirmed. J = 0.001 gives 0.03% shift | [Observer × Gravity Bridge](../experiments/OBSERVER_GRAVITY_BRIDGE.md) |
| Product states have local clock (Bell+ does not) | Confirmed. \|++⟩ crosses, Bell+ never | [Observer × Gravity Bridge](../experiments/OBSERVER_GRAVITY_BRIDGE.md) |
| States with α < 30° never cross, no observer time | Confirmed. CΨ stays below ¼ | [Observer × Gravity Bridge](../experiments/OBSERVER_GRAVITY_BRIDGE.md) |
| Crossing time shift encodes 1 bit (B measures vs not) | Confirmed. Δt = -0.218 at J=0.01 (2.87%) | [Observer × Gravity Bridge](../experiments/OBSERVER_GRAVITY_BRIDGE.md) |
| 21 pairs sufficient for 1 bit at 100% timing jitter | Confirmed. N_min = (σ/Δt)² | [Observer × Gravity Bridge](../experiments/OBSERVER_GRAVITY_BRIDGE.md) |
| Reservoir mechanism: nonlocal coherence destroyed by B's measurement | Confirmed. 0.82 → 0.00, A decays 4x faster. CORRECTED: coupling accelerates local decay vs isolation (not protects). Measurement cuts return flow of oscillation. | [Observer × Gravity Bridge](../experiments/OBSERVER_GRAVITY_BRIDGE.md) |
| Coupling accelerates local crossing vs single qubit (0.07x at J=1) | Confirmed. \|+,0⟩ with J=1: t=0.64 vs single \|+⟩: t=8.58 | [Observer × Gravity Bridge](../experiments/OBSERVER_GRAVITY_BRIDGE.md) |
| Measurement damage is timing-dependent (oscillation phase) | Confirmed. Max damage at t_B≈1.0 (61%), min at t_B≈0.01 (0%) | [Observer × Gravity Bridge](../experiments/OBSERVER_GRAVITY_BRIDGE.md) |
| No energy threshold for crossing: same ⟨H⟩, different crossing | Confirmed. ⟨H⟩ = J for all α, crossing depends on CΨ(0) | [Minimum Crossing Energy](../experiments/MINIMUM_CROSSING_ENERGY.md) |
| α_critical = 30° exactly, CΨ(0) = ¼ exactly | Confirmed. Binary search to 10⁻¹⁵ precision | [Minimum Crossing Energy](../experiments/MINIMUM_CROSSING_ENERGY.md) |
| Product states cross via Hamiltonian pumping (CΨ: 0 → 0.31) | Confirmed. \|0,1⟩ crosses, \|0,0⟩ does not | [Minimum Crossing Energy](../experiments/MINIMUM_CROSSING_ENERGY.md) |
| Eigenstates never cross (CΨ_max = 0, no dynamics) | Confirmed. \|+,+⟩, \|0,0⟩, \|1,1⟩ all CΨ_max = 0 | [Minimum Crossing Energy](../experiments/MINIMUM_CROSSING_ENERGY.md) |
| Critical J/γ ≈ 5-10 for \|0,1⟩ product state | Confirmed. CΨ_max crosses ¼ between J/γ = 5 and 10 | [Minimum Crossing Energy](../experiments/MINIMUM_CROSSING_ENERGY.md) |

---

## Summary by Tier

| Tier | Count | Examples |
|------|-------|---------|
| **Confirmed on hardware** | **24** in the [Confirmations registry](../compute/RCPsiSquared.Core/Confirmations/ConfirmationsRegistry.cs) (ibm_torino + Marrakesh/Kingston, 2026-02 to 2026-07); §1 details the 5 earliest Torino rows (the ¼ crossing and the Absorption ratio carry their own registry entries) | §1: ¼ crossing, T₂*/T₂, crossing eqn, x³+x=½, Absorption ratio 1.03. Registry: the 3 Torino runs, palindrome trichotomy, F25 cusp, F57 K_dwell, F83/F95, block-CΨ saturation, F120 moment tower, F84 heating leg, concentrator site contrast, F129 standing fringe |
| **Proven (analytical)** | 1 | Absorption Theorem: Re(λ) = −2γ⟨n_XY⟩, verified on 1,342 modes, CV=0 |
| **Computationally verified** | 55 (39 rows in §2 + 16 reopened J>0 rows in §9) | From γ·t_cross invariance, the θ trajectory, and the N-scaling barrier through the QKD closed-form family and no-signalling to the J>0 clock results; the rows themselves are the list |
| **Resolved anomaly (detuning)** | 1 (five measured signatures) | The late-time shadow: excess coherence, directionality, rising trend, boundary correlation, per-qubit shadow direction; cause settled as qubit-specific detuning (§3) |
| **Testable now** | 14 rows in §4 | Critical slowing, fingerprints, field threshold, ξ Markovianity diagnostic, the QKD forensics family (math verified; application retired), stealth angle existence, F120 moment tower (already flown and registered) |
| **Testable in principle** | 4 | θ compass, analog BH, voids, fractal decay |
| **Speculative** | 6 | CMB, BH burst, Page time, BH/WH unification, experienced time as crossing rate, anesthesia as C → 0 |
| **Unverified agent claims** | 4 | 33:1 ratio, linear scaling, H≠0 requirement, optimal C |
| **Null result** | 1 | Metric discrimination |
| **Closed hypothesis (J=0)** | 2 | Bridge dynamic (no-signalling), Bridge pre-encoded (= shared randomness). Reopened for J>0 via inter-qubit J-coupling. |
| **Falsified prediction** | 2 | Taxonomy noise-dependent (wrong: noise-independent), E=mγ² (wrong: α=2γ⟨n_XY⟩, linear not quadratic) |

---

*This document consolidates predictions from across the R = CΨ² framework.*
*For the proven algebra, see [Core Algebra](historical/CORE_ALGEBRA.md). For the interpretive framework, see [Interpretive Framework](../hypotheses/archive/INTERPRETIVE_FRAMEWORK.md).*
*For the phase boundary analysis, see [Dynamic Fixed Points](../experiments/DYNAMIC_FIXED_POINTS.md).*
*For weaknesses and honest self-assessment, see [Weaknesses](WEAKNESSES_OPEN_QUESTIONS.md).*
