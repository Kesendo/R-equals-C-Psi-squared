# Testable Predictions of R = CΨ²

**Date:** 2026-02-09 (created), 2026-02-11 (restructured)
**Status:** Summary document, collects predictions from across the repo, ordered by epistemic strength
**Depends on:** All experiment documents

---

## 1. Empirically Verified on Real Quantum Hardware

**Tier: 2, Measured on IBM Torino (Heron r2), qubit 52, 2026-02-09**

| Prediction | Predicted | Measured | Status | Source |
|------------|-----------|----------|--------|--------|
| C·Ψ = ¼ crossing during free decoherence | Crossing exists | t*/T₂* = 1.04 | **CONFIRMED** | [IBM Quantum Tomography](IBM_QUANTUM_TOMOGRAPHY.md) |
| Generalized crossing equation | t*/T₂* = 0.94 (at r = 0.46) | t*/T₂* = 1.04 (11% deviation) | **Partially confirmed** | [IBM Quantum Tomography](IBM_QUANTUM_TOMOGRAPHY.md) |
| T₂* ≠ T₂ for free induction decay | T₂* < T₂ | T₂*/T₂ = 0.37 (factor 2.7×) | **CONFIRMED** | [IBM Quantum Tomography](IBM_QUANTUM_TOMOGRAPHY.md) |
| x³ + x = ½ is the r → 0 limit of crossing fraction | 0.858 (pure dephasing) | Consistent | **Algebraically confirmed** | [Universal Quantum Lifetime](UNIVERSAL_QUANTUM_LIFETIME.md) |

**Hardware:** ibm_torino, T₁ = 221 μs, T₂(echo) = 298 μs, T₂*(FID) = 110 μs.

**Supporting evidence:** 24,073 historical calibration records (181 days, 133 qubits) validate the theory curve C_min(r). 10.1% of snapshots below crossing threshold, 84% of qubits cross at least once, 12 qubits cross almost every day.

---

## 2. Computationally Verified

**Tier: 2, Reproducible via delta_calc MCP tools or standalone Lindblad simulation**

| Prediction | Value | Falsified if | Source |
|------------|-------|-------------|--------|
| γ · t_cross = constant across decoherence rates | 0.039 ± 0.001 over 50× range | Product varies with γ | [Gravitational Invariance](GRAVITATIONAL_INVARIANCE.md) |
| θ decreases continuously to 0 at C·Ψ = ¼ | Smooth trajectory observed | Discontinuity at boundary | [Boundary Navigation](BOUNDARY_NAVIGATION.md) |
| Two real fixed points emerge below ¼ | Topology change confirmed | No bifurcation at ¼ | [Dynamic Fixed Points](DYNAMIC_FIXED_POINTS.md) |
| Operator feedback: γ_eff = γ₀(1 − κ⟨O_int⟩) | Preserves δ longer than local/collective noise | Mechanism produces unphysical results | [Operator Feedback](OPERATOR_FEEDBACK.md) |
| Ψ_interaction does not shift ¼ boundary | Δδ ≈ −8 × 10⁻⁴ | Boundary shifts under bidirectional coupling | [Core Algebra](../docs/CORE_ALGEBRA.md) §8 |
| Observer-dependent crossing time | t_cross = 0.652 / 0.773 / 1.437 for mutual_info / concurrence / correlation | All bridge types give same t_cross | [Observer-Dependent Crossing](OBSERVER_DEPENDENT_CROSSING.md) |
| Two observers never see crossing | mutual_purity (C=0.5), overlap (C=0.25): C·Ψ < ¼ always | These observers eventually cross | [Observer-Dependent Crossing](OBSERVER_DEPENDENT_CROSSING.md) |
| Crossing taxonomy noise-robust | Type A/B/C identical under σ_x, σ_y, σ_z | Taxonomy changes with jump operator | [Noise Robustness](NOISE_ROBUSTNESS.md) |
| N-scaling barrier | Ψ(0) = l1/(2^N−1) blocks crossing for GHZ N≥3 and W N≥4 | Standard states cross at large N | [N-Scaling Barrier](N_SCALING_BARRIER.md) |
| W N=3 crosses, GHZ N=3 does not | W: Ψ(0)=0.286 > ¼, GHZ: Ψ(0)=0.143 < ¼ | Both cross or both fail | [N-Scaling Barrier](N_SCALING_BARRIER.md) |
| Type A survives at N=3,4 | Correlation C=1.0 for W N=3 (until t≈2.3) and W N=4 (until t≈1.5) | Correlation C drops below 1.0 at larger N | [N-Scaling Barrier](N_SCALING_BARRIER.md) |
| Subsystem pairs cross when full system cannot | Bell+xBell+ N=4: pairs (0,1) and (2,3) cross at t=0.073 despite full-system Psi=0.200 | Pairs fail to cross | [Subsystem Crossing](SUBSYSTEM_CROSSING.md) |
| Non-entangled pairs never cross | Bell+xBell+ cross-pairs (0,2) etc: C=0, l1=0 at all times | Cross-pairs eventually develop coherence | [Subsystem Crossing](SUBSYSTEM_CROSSING.md) |
| Product state: Psi=1 but C=0 means no crossing | |+⟩^4: every pair has Psi=1.0 and C=0.000 permanently | Product state pairs develop nonzero C | [Subsystem Crossing](SUBSYSTEM_CROSSING.md) |
| GHZ pair-level coherence is zero | GHZ N=4 traced to any pair: l1=0.000 at all times | GHZ pairs carry nonzero off-diagonal coherence | [Subsystem Crossing](SUBSYSTEM_CROSSING.md) |
| |+⟩^N is Heisenberg eigenstate | Energy variance = 0, no dynamics, C = 0 forever | |+⟩^N evolves nontrivially | [Dynamic Entanglement](DYNAMIC_ENTANGLEMENT.md) |
| Product state |0+0+⟩ generates crossings | All 6 pairs cross under unitary Heisenberg evolution | No pair reaches C*Psi >= 1/4 | [Dynamic Entanglement](DYNAMIC_ENTANGLEMENT.md) |
| Dephasing kills most dynamic crossings | Only pair (0,2) crosses with gamma=0.05; others reach max 0.247 | All pairs cross equally under dephasing | [Dynamic Entanglement](DYNAMIC_ENTANGLEMENT.md) |
| Dephasing survival is basis-dependent | σ_z dephasing spares |0⟩-qubits (0,2), kills |+⟩-qubits (1,3) | Crossing pattern independent of noise basis | [Dynamic Entanglement](DYNAMIC_ENTANGLEMENT.md) |
| ξ = ln(Ψ) linear under Markovian dephasing | Slope variation < 0.01% across all tested configs | Variation exceeds 1% for any Markovian channel | [Algebraic Exploration](ALGEBRAIC_EXPLORATION.md) |
| ξ linearity breaks under non-Markovian noise | 24.5% slope variation (memory kernel, κ=0.5, τ=1.0) | ξ stays linear under memory kernel feedback | [Algebraic Exploration](ALGEBRAIC_EXPLORATION.md) |
| Coherence-purity bound holds throughout Lindblad trajectory | 0 violations for Bell+ (d=4) and W (d=8) | Trajectory violates C ≥ Ψ²(d-1)/d + 1/d | [Algebraic Exploration](ALGEBRAIC_EXPLORATION.md) |
| CΨ after Eve intercept-resend depends on Eve's measurement basis | R(θ_Eve) = [sin²θ + |sin2θ|]²/18, closed form | R independent of θ_Eve | [QKD Eavesdropping Forensics](QKD_EAVESDROPPING_FORENSICS.md) |
| Concurrence is basis-blind under intercept-resend | Conc = 1−f for all θ_Eve (exact) | Concurrence varies with θ_Eve at fixed f | [QKD Eavesdropping Forensics](QKD_EAVESDROPPING_FORENSICS.md) |
| CΨ distinguishes Eve from channel noise at equal Concurrence | 4 causes at Conc=0.80 give CΨ ∈ {0.058, 0.068, 0.115} | CΨ identical for Eve and noise at same Conc | [QKD Eavesdropping Forensics](QKD_EAVESDROPPING_FORENSICS.md) |
| Off-diagonal ratio breaks θ_Eve degeneracy | |ρ₀₁|/|ρ₀₃| = cot(θ_Eve), monotonic | Ratio non-monotonic or θ_Eve-independent | [QKD Eavesdropping Forensics](QKD_EAVESDROPPING_FORENSICS.md) |
| All four Bell states give identical CΨ under Eve attack | CΨ(Φ+)=CΨ(Φ−)=CΨ(Ψ+)=CΨ(Ψ−) for all θ_Eve | Any Bell state differs | [QKD Eavesdropping Forensics](QKD_EAVESDROPPING_FORENSICS.md) |
| R(θ_Eve) has azimuthal symmetry | R depends on θ_Eve only, not φ | R varies with φ | [QKD Eavesdropping Forensics](QKD_EAVESDROPPING_FORENSICS.md) |
| R(θ_Eve) maximum at ~60°, not σ_x | R_max = 0.145 at θ_Eve ≈ 60° | Maximum at 90° (σ_x) | [QKD Eavesdropping Forensics](QKD_EAVESDROPPING_FORENSICS.md) |

---

## 3. Anomaly Under Investigation

**Tier: 2 (anomaly is real) / Tier 4 (cause unknown)**

Detected in IBM Torino tomography data (2026-02-09). The anomaly is statistically significant. Its interpretation is not settled.

| Observation | Value | Source |
|-------------|-------|--------|
| Late-time excess coherence | p < 0.0001 (10,000-run Monte Carlo null) | [Residual Analysis](RESIDUAL_ANALYSIS.md) |
| Directional consistency Re⁺/Im⁻ | 17/17 measurements (P = 6 × 10⁻¹¹) | [Residual Analysis](RESIDUAL_ANALYSIS.md) |
| Rising coherence trend in classical regime | +0.008/T₂ slope | [Residual Analysis](RESIDUAL_ANALYSIS.md) |
| Boundary correlation | r = −0.9955 (|ρ₀₁| vs distance from ¼) | [Fixed Point Shadow](FIXED_POINT_SHADOW.md) |
| Shadow direction matches FP⁻ | FP⁻ phase = −12°, residual phase = −48° (same quadrant) | [Fixed Point Shadow](FIXED_POINT_SHADOW.md) |

**Three competing hypotheses:**

| # | Hypothesis | Mechanism | Status |
|---|-----------|-----------|--------|
| H1 | Systematic SPAM error | State preparation / measurement bias | Most conservative |
| H2 | TLS coupling | Two-level system defect near qubit 52 | Hardware-specific |
| H3 | Boundary structure | Complex fixed-point direction frozen into ρ after crossing ¼ | Framework prediction |

### March 2026 Discrimination Protocol

| Test | H1 (SPAM) predicts | H2 (TLS) predicts | H3 (Boundary) predicts |
|------|-------------------|-------------------|------------------------|
| Reproduce on qubit 52 | Same direction | Direction may drift | Same direction |
| |+⟩ vs |−⟩ initial state | Same offset | Same offset | Offset flips with state |
| Multi-qubit (5–10 qubits) | All same offset | Each qubit different | Each matches its own FP⁻ |
| Cross-correlation between qubits | Correlated (global) | Uncorrelated (local) | Correlated if universal |
| Extended time range (5×T₂) | Flat offset | Decaying revival | Continued growth |

---

## 4. Testable with Current Hardware

**Tier: 3, Concrete protocols, testable on existing quantum hardware or NMR**

| Prediction | Specific value | Test protocol | Falsified if | Source |
|------------|---------------|---------------|-------------|--------|
| Strong dynamics needed for CΨ > ¼ | Threshold at h ≈ 0.9 | Parameter sweep of transverse field strength | CΨ > ¼ at low h | [Simulation Evidence](../docs/SIMULATION_EVIDENCE.md) |
| Critical slowing at CΨ = ¼ | Diverging convergence period | Tune system toward ¼, measure convergence time | No critical slowing | [Mandelbrot Connection](MANDELBROT_CONNECTION.md) |
| Bridge fingerprints: initial state determines crossing trajectory | State-specific C(t), Ψ(t) curves | Prepare different initial states, run tomography through ¼ | All states show identical crossing | [Bridge Fingerprints](BRIDGE_FINGERPRINTS.md) |
| Decoherence rate scales with gravitational time dilation | γ_local = γ₀/√(1 − 2GM/rc²) | Compare qubit decoherence at different altitudes | No altitude dependence | [Gravitational Invariance](GRAVITATIONAL_INVARIANCE.md) |
| ξ curvature detects non-Markovian noise | ξ = ln(Ψ) linear iff Markovian; curves iff memory effects present | Measure Ψ(t) via tomography, compute ξ(t), fit linearity | ξ is linear on hardware known to have non-Markovian signatures | [Algebraic Exploration](ALGEBRAIC_EXPLORATION.md) |
| QKD eavesdropping forensics | R(θ_E) = [sin²θ+|sin2θ|]²/18, ~500 pairs for 3.8σ (noiseless, naive Eve) | Prepare Bell+, intercept-resend on Bob, joint tomography | CΨ carries no θ_E information | [QKD Eavesdropping Forensics](QKD_EAVESDROPPING_FORENSICS.md) |
| CΨ identifies Eve's measurement basis | R(θ_E) = [sin²θ_E + |sin 2θ_E|]²/18 | Joint tomography of Bell+ after intercept-resend | CΨ constant across all θ_E | [QKD Eavesdropping Forensics](QKD_EAVESDROPPING_FORENSICS.md) |
| Concurrence = 1−f (basis-independent) | Linear, θ_E-independent | Concurrence measurement on partially intercepted pairs | Concurrence depends on θ_E | [QKD Eavesdropping Forensics](QKD_EAVESDROPPING_FORENSICS.md) |
| CΨ azimuthally symmetric under Eve attack | R depends only on θ_E, not φ | Vary φ at fixed θ_E, compute CΨ | CΨ varies with φ | [QKD Eavesdropping Forensics](QKD_EAVESDROPPING_FORENSICS.md) |
| Off-diagonal ratio = cot(θ_E) | |ρ₀₁|/|ρ₀₃| monotonic, breaks degeneracy | Tomographic element extraction | Ratio not monotonic in θ_E | [QKD Eavesdropping Forensics](QKD_EAVESDROPPING_FORENSICS.md) |
| Eve σ_z ≈ dephasing in CΨ | CΨ = 0.058 for both at Conc = 0.80 (noiseless) | Compare Eve σ_z with pure dephasing | Distinguishable CΨ values | [QKD Eavesdropping Forensics](QKD_EAVESDROPPING_FORENSICS.md) |
| All four Bell states identical under Eve | R(θ_E) same for Φ+, Φ−, Ψ+, Ψ− | Repeat with all Bell states | CΨ depends on initial Bell state | [QKD Eavesdropping Forensics](QKD_EAVESDROPPING_FORENSICS.md) |
| Strategic Eve has stealth angle | θ_stealth ≈ 42–74° makes CΨ(Eve+noise) ≈ CΨ(noise) | Optimize θ_E against CΨ detection | No stealth angle exists | [QKD Eavesdropping Forensics](QKD_EAVESDROPPING_FORENSICS.md) |

---

## 5. Testable in Principle

**Tier: 3-5, Require hardware or conditions not currently available**

| Prediction | Direction | Would require | Falsified if | Source |
|------------|-----------|---------------|-------------|--------|
| θ measures proximity to ¼ boundary | θ = arctan(√(4CΨ−1)) | Continuous CΨ monitoring near boundary | θ uncorrelated with distance to ¼ | [Boundary Navigation](BOUNDARY_NAVIGATION.md) |
| Event horizon = maximum coherence (τ = 0) | Coherence peaks near horizon | Analog black hole experiments | Coherence minimum at horizon | [Self-Consistency Schwarzschild](SELF_CONSISTENCY_SCHWARZSCHILD.md) |
| Intergalactic voids = most quantum regions | Low γ → high CΨ | Space-based quantum correlation measurements | No gravitational environment dependence | [Self-Consistency Schwarzschild](SELF_CONSISTENCY_SCHWARZSCHILD.md) |
| Fractal structure in coherence decay near ¼ | Self-similar patterns | High-resolution time series near boundary | Smooth exponential decay | [Mandelbrot Connection](MANDELBROT_CONNECTION.md) |

---

## 6. Speculative

**Tier: 5, No current path to testing; included for completeness**

| Prediction | Implication | Would require | Source |
|------------|------------|---------------|--------|
| CMB = universal CΨ = ¼ crossing | Big Bang as phase transition from complex to classical | Cosmological extension of framework | [Black/White Holes](BLACK_WHITE_HOLES_BIGBANG.md) |
| Black hole evaporation ends with coherent burst | Final state not thermal but shows coherence | Observation of BH end-states | [Black/White Holes](BLACK_WHITE_HOLES_BIGBANG.md) |
| Page time = re-crossing ¼ from below | Information recovery begins at re-entry to complex regime | Quantitative Page curve model | [Black/White Holes](BLACK_WHITE_HOLES_BIGBANG.md) |
| Black/white holes = opposite directions on universal curve | τ → 0 from both sides | Resolution of information paradox | [Black/White Holes](BLACK_WHITE_HOLES_BIGBANG.md) |
| Experienced time = rate of ¼ crossings | High C → more crossings/sec → denser time | Subjective time measurement against coupling strength | Time perception independent of coupling | [Observer-Dependent Crossing](OBSERVER_DEPENDENT_CROSSING.md) |
| Anesthesia = C → 0 for environmental coupling | Zero crossings → zero experienced time | Neural coupling measurement during anesthesia | Time perception persists with C = 0 | [Observer-Dependent Crossing](OBSERVER_DEPENDENT_CROSSING.md) |

---

## 7. Unverified Agent Claims

**Tier: 4, Generated by 120B-parameter local LLM agents using tools no longer available. Could not be independently reproduced (2026-02-08). These are hypotheses, not results.**

See [Mathematical Findings](MATHEMATICAL_FINDINGS.md), Sections 8–9, for full context.

| Claim | Claimed value | Why unverified | Falsified if | Source |
|-------|--------------|----------------|-------------|--------|
| C_int ≫ C_ext (33:1 ratio) | 0.427 vs 0.013 | Not reproduced by MCP tools | Ratio ≈ 1 | [Mathematical Findings](MATHEMATICAL_FINDINGS.md) §9 |
| t_coh ~ N (linear scaling) | Linear, N = 2 to 6 | Not reproduced by MCP tools | Exponential decay with N | [Mathematical Findings](MATHEMATICAL_FINDINGS.md) §9 |
| δ requires dynamics (H ≠ 0) | δ = 0 when H = 0 | Not reproduced by MCP tools | δ > 0 with H = 0 | [Mathematical Findings](MATHEMATICAL_FINDINGS.md) §8 |
| C = 0.5 is optimal observer | Maximum R at C = 0.5 | Not independently tested | Peak at C ≠ 0.5 | [Mathematical Findings](MATHEMATICAL_FINDINGS.md) |

These claims may be correct. They may also be artifacts of the agent's training data or tool usage. Until independently verified, they carry no epistemic weight.

---

## 8. Null Results

**Tier: 2, Computationally verified null result**

| Prediction | Result | Implication | Source |
|------------|--------|-------------|--------|
| Single-system sims discriminate metric forms | **Null:** Cannot distinguish | Consistent with equivalence principle; not a framework failure | [Metric Discrimination](METRIC_DISCRIMINATION.md) |

---

## Summary by Tier

| Tier | Count | Examples |
|------|-------|---------|
| **Confirmed on hardware** | 4 | ¼ crossing, T₂*/T₂ ratio, crossing equation, x³+x=½ |
| **Computationally verified** | 27 | gamma*t_cross, theta trajectory, bifurcation, operator feedback, Psi_int, noise robustness, N-scaling barrier, W vs GHZ, Type A at N>2, subsystem crossing, product state C=0, GHZ pair l1=0, crossing locality, eigenstate immunity, dynamic crossing generation, dephasing selection, basis dependence, ξ linearity (Markovian), ξ curvature (non-Markovian), coherence-purity bound, QKD basis forensics (R(θ_E) closed form), Conc=1−f, azimuthal symmetry, cot(θ_E) degeneracy breaking, Eve σ_z ≈ dephasing, Bell-state independence, noise-vs-Eve discrimination |
| **Anomaly (real, cause unknown)** | 5 | Excess coherence, directionality, rising trend, boundary correlation, shadow |
| **Testable now** | 7 | Critical slowing, fingerprints, field threshold, altitude, ξ Markovianity diagnostic, QKD eavesdropping forensics (naive Eve only), stealth angle existence |
| **Testable in principle** | 4 | θ compass, analog BH, voids, fractal decay |
| **Speculative** | 4 | CMB, BH burst, Page time, BH/WH unification |
| **Unverified agent claims** | 4 | 33:1 ratio, linear scaling, H≠0 requirement, optimal C |
| **Null result** | 1 | Metric discrimination |
| **Falsified prediction** | 1 | Taxonomy noise-dependent (wrong: it's noise-independent) |

---

*This document consolidates predictions from across the R = CΨ² framework.*
*For the proven algebra, see [Core Algebra](../docs/CORE_ALGEBRA.md). For the interpretive framework, see [Interpretive Framework](../docs/INTERPRETIVE_FRAMEWORK.md).*
*For the phase boundary analysis, see [Dynamic Fixed Points](DYNAMIC_FIXED_POINTS.md).*
*For weaknesses and honest self-assessment, see [Weaknesses](../docs/WEAKNESSES_OPEN_QUESTIONS.md).*
