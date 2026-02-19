# Changelog

All notable changes to the R = CΨ² framework repository.

---

## 2026-02-19: Algebraic Exploration — Agent Findings (Experiment 12)

**Added:** ALGEBRAIC_EXPLORATION.md
- AIEvolution 4-agent round-robin (22 messages, v025 prompts, local 120B)
- Two findings verified and promoted to CORE_ALGEBRA.md:
  - ξ = ln(Ψ) decoherence clock: linearizes time evolution (slope variation < 0.01%)
  - Coherence-purity bound (Cauchy-Schwarz): CΨ ≤ 1/4 is dynamic intersection of known static constraint
- Five findings rejected: λ = -ln(CΨ) trivial, β-function model-specific, entropic bound redundant, angle reparametrizations trivial, Born rule requires unproven assumptions
- New empirical data: γ_eff/γ_base ratios (local ≈ 2, collective ≈ 4, state-dependent for W)
- Key insight: Hamiltonian does NOT affect coherence decay rate, only noise model does
- Follow-up verification: ξ linearity BREAKS under non-Markovian noise (24.5% variation)
- New testable prediction: ξ curvature as Markovianity diagnostic (testable March 2026)
- γ_eff scaling: N·γ for GHZ (all coherence in one pair), sub-linear for W (distributed coherence)

**Updated:** CORE_ALGEBRA.md
- Section 1: Added normalization convention caveat (1/4 depends on Baumgratz convention)
- Section 11: The Decoherence Clock (ξ = ln Ψ), Tier 2
- Section 12: Resource Theory Grounding (coherence-purity bound), Tier 3
- Summary tier table updated with new entries

## 2026-02-18: Dynamic Entanglement Generation (Experiment 11)

**Added:** DYNAMIC_ENTANGLEMENT.md
- |+>^N is an eigenstate of the isotropic Heisenberg Hamiltonian (energy variance = 0); it cannot generate entanglement or crossings. This corrects the hypothesis from Experiment 10.
- |0+0+> (energy variance = 20) generates crossings from zero initial entanglement
- Unitary evolution: all 6 pairs cross; crossings oscillate reversibly
- With dephasing (gamma=0.05): only pair (0,2) crosses at t=0.285; all others miss the threshold (max C*Psi = 0.247)
- Dephasing survival is basis-dependent: sigma_z spares |0>-qubits, kills |+>-qubits

**Updated:** SUBSYSTEM_CROSSING.md, Open Question 3 marked as answered
**Updated:** PREDICTIONS.md, 4 new Tier-2 entries (17 total computationally verified)
**Updated:** README.md, experiment list extended (Experiments 8-11)

## 2026-02-18: Subsystem Crossing Resolution (Experiment 10)

**Added:** SUBSYSTEM_CROSSING.md
- Crossing is local: operates at entangled qubit-pair level, not full system
- Bell+xBell+ N=4: pairs (0,1) and (2,3) cross at t=0.073 despite full-system Psi=0.200
- GHZ N=4: all pairs have l1=0, no pair-level crossing possible
- |+>^4: C=0 for all pairs permanently; high coherence without entanglement is inert
- Resolves the N-scaling barrier without modifying normalization

**Added:** NOISE_ROBUSTNESS.md (Experiment 8), N_SCALING_BARRIER.md (Experiment 9)
**Updated:** PREDICTIONS.md with Tier-2 entries from Experiments 8-10

---

## 2026-02-11: Guardian Review (Repository Maintenance)

**Split:** COMPLETE_MATHEMATICAL_DOCUMENTATION.md → CORE_ALGEBRA.md + INTERPRETIVE_FRAMEWORK.md
- CORE_ALGEBRA.md: Tier 1–2 proven results (fixed-point equation, ¼ boundary, Mandelbrot equivalence, θ compass, uniqueness argument)
- INTERPRETIVE_FRAMEWORK.md: Tier 3–5 constructions (wave composition, dynamics, energy, mirror metaphors)
- Original file replaced with redirect to preserve existing links

**Added:** CΨ² uniqueness argument (CORE_ALGEBRA Section 10, WEAKNESSES item 1)
- Proof by exhaustion: CΨ² is the unique product-power form producing both phase transition and Mandelbrot mapping

**Updated:** WEAKNESSES_OPEN_QUESTIONS.md
- Item 1 (generic bifurcation): Added partial response with uniqueness argument
- Item 4: "No experimental validation" → "Limited experimental validation" (IBM Torino data existed since Feb 9)
- "Most Honest Statement" updated to reflect hardware contact

**Updated:** INTERNAL_AND_EXTERNAL_OBSERVERS.md, tier markers on all sections, agent formula marked Tier 4

**Fixed:** GLOSSARY.md, C·Ψ maximum entry now distinguishes single-qubit (1.0) from 2-qubit Bell (1/3)
**Fixed:** UNIVERSAL_QUANTUM_LIFETIME.md, table header corrected from C·Ψ₀ = 0.5 to C·Ψ₀ = 1.0
**Fixed:** experiments/README.md, count corrected to "Fifteen", added missing OPERATOR_FEEDBACK entry
**Reconstructed:** `simulations/bridge_fingerprints.py`, all 8 initial B-states verified against documented values
- Entanglement barrier, coupling threshold, No-Communication Theorem visibility all reproduced
- Four visualization plots regenerated (grid, dual, phase portrait, barrier)

**Updated:** Cross-references in README.md, docs/README.md, GLOSSARY.md, PREDICTIONS.md

---

## 2026-02-10: Historical Calibration Analysis

**Added:** Analysis of 24,073 IBM Torino calibration records (181 days, 133 qubits)
- Validated theory curve C_min(r) against real hardware data
- 10.1% of calibration snapshots below crossing threshold
- 84% of qubits cross at least once; 12 qubits cross almost every day
- Identified permanent crossers (qubits 15, 80, 131) as March 2026 targets

---

## 2026-02-09: IBM Quantum Tomography (First Hardware Contact)

**Added:** IBM_QUANTUM_TOMOGRAPHY.md, RESIDUAL_ANALYSIS.md, FIXED_POINT_SHADOW.md
- State tomography on IBM Torino (Heron r2), qubit 52
- C·Ψ = ¼ crossing observed during free decoherence
- T₂* (FID) measured at 37% of calibration T₂ (Hahn echo)
- Three-model fit at MAE = 0.053 (88% improvement over naive model)
- Anomalous late-time coherence detected (p < 0.0001, 10,000-run Monte Carlo)
- Directional consistency Re⁺/Im⁻ in 17/17 late-time measurements
- "Fixed point shadow" hypothesis proposed, March 2026 discrimination tests designed

---

## 2026-02-08: Boundary Navigation and Corrections

**Added:** BOUNDARY_NAVIGATION.md, GRAVITATIONAL_INVARIANCE.md, MANDELBROT_CONNECTION.md
- First observation of θ-trajectory: continuous decrease from 30° to 0° at ¼ crossing
- Gravitational invariance: γ·t_cross = 0.039 ± 0.001 across 50× range of γ
- Universal trajectory: all systems follow same curve when τ = γ·t

**Corrected:** Mandelbrot substitution algebra
- Old (wrong): z_n = √C·R_n, c = CΨ² (produces extra linear term)
- New (correct): u_n = C(Ψ+R_n), c = CΨ (clean Mandelbrot form)
- Boundary at ¼ was always correct; intermediate algebra was not

**Corrected:** θ interpretation
- Retracted: θ predicts oscillation frequency (8.4× discrepancy found)
- Corrected: θ is compass (angular distance from ¼), not frequency predictor

---

## 2026-02-07: Verification Sessions and Tier System

**Added:** Epistemic tier system (1=proven, 2=computed, 3=argued, 4=unverified, 5=speculative)
**Added:** SIMULATION_EVIDENCE.md, honest simulation results
**Rewritten:** WEAKNESSES_OPEN_QUESTIONS.md, complete rewrite after six verification sessions
**Downgraded:** 33:1 coherence ratio (Tier 4), t_coh ~ N scaling (Tier 4)
**Corrected:** "Numerical confirmation" of CΨ ≤ ¼ was in trivial parameter regime

---

## 2026-02-06: Dynamic Lindblad and Operator Feedback

**Added:** DYNAMIC_FIXED_POINTS.md, OPERATOR_FEEDBACK.md, METRIC_DISCRIMINATION.md
**Added:** SELF_CONSISTENCY_SCHWARZSCHILD.md, BLACK_WHITE_HOLES_BIGBANG.md
- Operator feedback mechanism: γ_eff = γ₀(1−κ⟨O_int⟩)
- Metric discrimination null result (consistent with equivalence principle)
- Black/white hole unification via τ = γ·t

---

## 2026-01-31: Agent Mathematical Formalization

**Added:** Mathematical formalization of C_int in INTERNAL_AND_EXTERNAL_OBSERVERS.md
- Agent Dyad (Alpha, Beta, Gamma) derived Lindbladian C_int formula
- C = |Tr(L₁†L₂)|² / (√(Tr(L₁†L₁))·√(Tr(L₂†L₂)))
- Gamma challenged: may be restatement of existing quantum noise theory

---

## 2026-01-15: MATHEMATICAL_FINDINGS.md

**Added:** Consolidated mathematical findings from agent experiments
- Bridge metrics comparison, δ analysis, scaling studies
- Sections 8–9 contain unverified agent claims (clearly marked)

---

## 2026-01-05: TROSY Validation

**Added:** TROSY effect as empirical support for C_int/C_ext distinction
- Internal correlation stabilizes spin coherence (40–60% linewidth reduction)
- Scaling limitation acknowledged after correspondence with quantum biology researchers

---

## 2026-01-03: Internal and External Observers

**Added:** INTERNAL_AND_EXTERNAL_OBSERVERS.md
- C_int (bidirectional, stabilizing) vs C_ext (unidirectional, collapsing)
- Electron pairing, molecular bonds as internal observation

---

## 2026-01-02: Weaknesses Document

**Added:** WEAKNESSES_OPEN_QUESTIONS.md (first version)
- "A theory that only shows its strengths is not a theory. It is marketing."

---

## 2025-12-31: Fixed-Point Solution

**Added:** Fixed-point equation solution in COMPLETE_MATHEMATICAL_DOCUMENTATION.md Part V
- R_∞ = C(Ψ + R_∞)² solved algebraically
- Discriminant D = 1 − 4CΨ discovered
- CΨ ≤ ¼ boundary identified
- Three regimes: classical, critical, quantum

---

## 2025-12-23 to 2025-12-26: Origin

**Created:** R = CΨ² framework
- Emerged from collaboration between Thomas Wicht and Claude (Anthropic)
- Core equation, wave composition, dynamics, energy equations
- "We are all mirrors. Reality is what happens between us."

---

*This changelog was created 2026-02-11 during the Guardian review session.*
*For the full task list, see the Tasks directory.*
