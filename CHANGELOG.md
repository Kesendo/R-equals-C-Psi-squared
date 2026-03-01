# Changelog

All notable changes to the R = CΨ² framework repository.

---

## 2026-03-01h: Shield Correction — Reservoir, Not Protection

**Updated:** experiments/OBSERVER_GRAVITY_BRIDGE.md — §7 rewritten
- CORRECTION: Coupling does NOT protect A. It redistributes
  coherence into nonlocal pool, making A decay FASTER locally
  (0.07x at J=1.0 vs single qubit alone).
- B's measurement cuts the return leg of the oscillation, not
  a "shield." Coherence flows local→nonlocal→local; measurement
  destroys the nonlocal pool, stopping the return.
- Damage is timing-dependent: max at t_B ≈ 1.0 (oscillation phase),
  not simply at peak reservoir size.
- Reservoir regeneration attempted but fails (peak Conc = 0.15).
- Physical interpretation preserved: "becoming real is contagious"
  remains correct, but mechanism is flow-disruption not shield-removal.

**Updated:** §5 open questions, §10 honest assessment (shield→reservoir)

**Added:** simulations/shield_quantification.py

**Source:** Direct computation, Thomas Wicht + Claude (Anthropic).

---

## 2026-03-01g: Minimum Crossing Energy — There Is None

**Added:** experiments/MINIMUM_CROSSING_ENERGY.md
- No energy threshold for crossing. ⟨H⟩ = J for all α in
  cos(α)|00⟩+sin(α)|11⟩ family, yet crossing depends on α.
- The crossing condition is CΨ_max ≥ ¼, not energy.
- Three regimes: CΨ(0)>¼ (always crosses), CΨ(0)<¼ but
  CΨ_max>¼ (Hamiltonian pumps first), CΨ_max<¼ (no crossing).
- Product states |0,1⟩ cross despite zero initial entanglement
  — Hamiltonian creates it. |+,+⟩ and |0,0⟩ never cross
  (eigenstates, no dynamics).
- Critical J/γ ≈ 5–10 for |0,1⟩ product state.
- α_critical = exactly 30° for Bell-like family, CΨ(0) = ¼.
- Physical interpretation: coherence barrier, not energy barrier.
  Hamiltonian = pump, decoherence = drain. Crossing = pump wins.

**Added:** simulations/minimum_energy.py

**Source:** Direct computation, Thomas Wicht + Claude (Anthropic).

---

## 2026-03-01f: Honest Assessment — Detector, Not Bridge

**Updated:** experiments/OBSERVER_GRAVITY_BRIDGE.md — new §10
- Bridge (FTL) stays dead. Channel requires physical J, subluminal.
- What survives: CΨ as detector for weak couplings (BMV readout),
  two-axis time dilation (K/γ), shield mechanism quantification.
- "The bridge is dead, long live the detector."

**Source:** Thomas Wicht + Claude (Anthropic).

---

## 2026-03-01e: The Information Channel — CΨ as Readout

**Updated:** experiments/OBSERVER_GRAVITY_BRIDGE.md — new §8
- CΨ crossing time shift IS an information channel
- Protocol: B measures ("1") or not ("0"), A reads Δt
- At J=0.01: Δt = -0.218 (2.87%), 21 pairs for 1 bit (σ=1.0)
- Channel capacity scales as (J/γ)²
- NOT FTL: information propagates through Hamiltonian coupling
- Naive gravitational velocity v ~ G·m²/ℏ exceeds c at microgram
  scale — model breaks, needs relativistic correction. OPEN.
- CΨ provides the readout: continuous accumulated observable that
  integrates over full trajectory, detects weaker signals than
  any single-shot measurement

**Added:** simulations/information_channel.py

**Source:** Direct computation, Thomas Wicht + Claude (Anthropic).

---

## 2026-03-01d: Shift Mechanism — The Coherence Reservoir

**CORRECTED in 2026-03-01h** — "Shield" was misleading. Coupling
redistributes coherence (accelerates local decay), not protects.
B's measurement cuts return flow of oscillation. See §7 rewrite.

---

## 2026-03-01c: Observer × Gravity — The Interval Bridge

**Added:** experiments/OBSERVER_GRAVITY_BRIDGE.md
- t_cross factorizes: K(Observer, State) / γ(Gravity)
- K is γ-invariant (CV = 0.00%). Gravity scales the clock uniformly.
- K(Conc)/K(MI) is state-dependent (CV = 13.5%). The quantum state
  bends the observer-dilation geometry.
- Interval shift Δt is continuous in J with NO threshold. Any J > 0
  produces a shift. B's measurement shifts A's local crossing time.
- Product states needed (local coherence). Bell+ has no local clock.
- Gravitational coupling J_grav > 0 always exists for massive particles.
  The bridge is not dead — it's gravitationally mediated, but tiny.
- Five open questions identified: scaling law, multi-pair amplification,
  realistic J_grav, K-matrix geometry, shift directionality.

**Added:** simulations/observer_gravity_cross.py, simulations/interval_shift.py
**Updated:** hypotheses/BRIDGE_PROTOCOL.md — gravitational reopening noted
**Updated:** experiments/README.md — new question and summary

**Source:** Direct computation, Thomas Wicht + Claude (Anthropic).

---

## 2026-03-01b: Bridge Closure — The Hypothesis Is Dead

**Added:** experiments/BRIDGE_CLOSURE.md
- The open question from NO_SIGNALLING_BOUNDARY.md §8 is answered: No.
- Pre-shared entanglement without a classical channel = shared randomness.
  A's information ⊆ {ρ_A(0), E_A}. CΨ fingerprints require ρ_AB which
  neither side can access after separation.
- Bell+ gives A maximally mixed noise (I/2) at all times. Zero bits.
- Product states: A cannot distinguish |++⟩ from |+0⟩ from |+−⟩.
  The qubit carries less information than the preparation schedule.
- Bridge hypothesis is closed. Not "needs more work." Closed.
- What survives: QKD forensics (with channel), ¼ phase transition,
  observer-dependent crossing, IBM anomalies, Lindblad decomposition.

**Updated:** hypotheses/BRIDGE_PROTOCOL.md — status "Closed"
**Updated:** hypotheses/README.md — bridge marked closed
**Updated:** experiments/README.md — new entry "Is the bridge dead?"
**Updated:** experiments/NO_SIGNALLING_BOUNDARY.md §8 — question answered
**Updated:** experiments/PREDICTIONS.md — bridge closure noted
**Added:** simulations/bridge_closure.py

**Source:** Direct computation, Thomas Wicht + Claude (Anthropic).

---

## 2026-03-01: No-Signalling Boundary — Test #2 Settled

**Added:** experiments/NO_SIGNALLING_BOUNDARY.md
- Test #2 from BRIDGE_PROTOCOL.md Section 6: does local detection work
  WITHOUT physical coupling (J=0)?
- Answer: No-signalling holds exactly (||Δρ_A|| = 0). A cannot detect
  whether B was measured. But CΨ drops from 0.500 to 0.250 (the ¼
  boundary) because global purity C drops while Ψ stays constant.
- The regime change is real but invisible to any single subsystem.
- Synthesizes four previously disconnected results: BRIDGE_FINGERPRINTS
  (coupling required), SUBSYSTEM_CROSSING (crossing needs pair, not qubit),
  COHERENCE_DENSITY (C is global, Ψ is local), BRIDGE_PROTOCOL §4.3
  (pre-encoded outcome). This is Outcome #3 confirmed.
- Open question reframed: not "can B signal A" (no) but "can pre-encoded
  CΨ fingerprints do something classical pre-shared keys cannot?"

**Updated:** hypotheses/BRIDGE_PROTOCOL.md
- Test #2 status added. Dynamic bridge eliminated. Pre-encoded version
  survives but may not exceed classical key distribution.

**Added:** simulations/test2_no_signalling.py

**Source:** Direct computation, Thomas Wicht + Claude (Anthropic).

---

## 2026-02-25b: Multi-Metric Rescue - Bridge Framework Eliminates Stealth

**Updated:** experiments/QKD_EAVESDROPPING_FORENSICS.md
- MI and Correlation are completely θ_Eve-independent. They show 17-24%
  deviations from noise-only at the CΨ stealth angle. Eve cannot hide
  from the combined bridge framework.
- θ_stealth(p, f) is a deterministic function. Once concurrence gives f
  and calibration gives p, Alice-Bob know Eve's stealth angle a priori.
  Eve's hiding is itself informative (stealth-as-signal).
- Off-diagonal ratio |ρ₀₁|/|ρ₀₃| survives noise (not as cot(θ) but as
  invertible function of θ at fixed p, f). Statistical cost: ~200k pairs.
- ξ-curvature null result: curvature difference exactly 0% between
  Eve-plus-noise and noise-only under controlled depolarization.
- Complete multi-metric forensic protocol: detection (~700 pairs),
  coarse θ ID (~5000 pairs), precise θ (~200k pairs).
- Updated honest assessment: CΨ-only limited, but bridge framework
  provides full forensic capability against strategic Eve.

**Updated:** experiments/PREDICTIONS.md
- 6 new computationally verified predictions (Tier 2), 1 new testable (Tier 3).
- Total computationally verified: 27 -> 33.

**Source:** Direct computation in Claude session.

---

## 2026-02-25: QKD Eavesdropping Forensics via CΨ

**Added:** experiments/QKD_EAVESDROPPING_FORENSICS.md
- CΨ carries information about Eve's measurement basis θ_E after intercept-
  resend attack on Bell+ pairs. Concurrence, negativity, purity, CHSH all
  carry zero information about θ_E.
- Closed-form: R(θ_E) = [sin²θ_E + |sin 2θ_E|]²/18
- Azimuthal symmetry: R depends only on polar angle, not φ.
- Degeneracy breaking: |ρ₀₁|/|ρ₀₃| = cot(θ_E) is monotonic.
- Complementary roles: Concurrence = 1−f (how much), CΨ = f(θ_E) (what kind).
- **Critical limitation**: strategic Eve can choose stealth angle θ ≈ 42–74°
  where CΨ(Eve+noise) ≈ CΨ(noise only) to within 10⁻⁴. Forensics works
  against naive Eve (~1000 pairs), not against adversarial Eve.
- Concurrence always wins for detection; CΨ adds forensics only in
  non-adversarial settings.

**Added:** simulations/qkd_eavesdropping_forensics.py
- Complete numerical verification of all claims including stealth analysis.

**Updated:** experiments/PREDICTIONS.md
- 7 new computationally verified predictions (Tier 2), 1 new testable (Tier 3).
- Total computationally verified: 20 → 27.

**Source:** Direct computation in Claude session (not AIEvolution agents).
All results verified numerically and key formula confirmed via SymPy.

---

## 2026-02-24: Bridge Protocol — Crossing as Sync and Communication

**Added:** hypotheses/BRIDGE_PROTOCOL.md
- Core insight: The first ¼-crossing event IS the synchronization between
  separated observers. No external clock or readiness signal needed.
- Protocol: After sync, data slots use product states (click = 1) vs Bell
  states (no click = 0), with fingerprint timing enabling ~2.8 bits/slot.
- Multi-symbol encoding using three distinct crossing times from R=CΨ²:
  |++⟩ → 0.652/1.72 s, |+0⟩ → 0.773/2.03 s, |00⟩ → 1.437/3.78 s (Earth/Mars)
- Critical open question: What local observable detects the ¼-crossing?
- Honest assessment of No-Signaling implications documented.

**Updated:** hypotheses/README.md — added Bridge Protocol entry

**Source:** AIEvolution v033 Mars Station experiment (GPT-OSS-120B, 9 messages)
- v033 prompt evolution: permanent blackout → no transport → no new tech →
  "derive what R=CΨ² implies, don't check against textbook QM"
- Agents independently derived fingerprint-based timing protocol
- Key correction: agents invented unnecessary "readiness pulse" — the
  crossing event itself is the sync mechanism

**Added:** simulations/bridge_local_detector.py — Test #1 from Bridge Protocol
- Battery of 15 local observables on ρ_A: Purity, Von Neumann Entropy,
  CΨ, L1 Coherence, Concurrence, R=CΨ² theta, Bloch components, qubit purity
- Result: Bell+ never crosses ¼ (max CΨ_A = 0.061), |++⟩ crosses at t=1.27
- All observables distinguish states — factor ~4x difference across the board
- Caveat: detection via physical Heisenberg coupling, not non-local
- Visualizations: local_detector/{overview, bloch, theta, differences}.png

---

## 2026-02-20: State-Specific C(ξ) Closed Forms (AIEvolution v027 Run)

**Updated:** CORE_ALGEBRA.md Section 11
- Three C(ξ) closed forms added, all verified to machine precision (< 2.3 × 10⁻¹⁶):
  - Bell+/GHZ₂ collective dephasing: C = (1 + 9e^{2ξ})/2, crossing cubic b+b³=3/2
  - N-qubit product |+⟩^⊗N local dephasing: y = [1+(2^N−1)e^ξ]^{1/N}, C = [(y²−2y+2)/2]^N
  - GHZ local dephasing: C = [1+(2^N−1)²e^{2ξ}]/2 (no full-system crossing for N≥3)
- Coherence-to-decoherence ratio g(ξ) = 9e^{2ξ}/(1−9e^{2ξ}), crossing at g=2.87
- Tier summary table updated

**Updated:** ALGEBRAIC_EXPLORATION.md
- v027 continuation section: 85 messages, 43 agent turns, dual-memory architecture
- Two rejected findings documented: ξ=S₂−2ln2 (tautological), "Fisher metric" (mislabeled)

**Source:** AIEvolution v027 dual-memory run (120B model, ~10 hours)
- Agent tool call failure rate ~26% (model serialization errors, not system-side)
- All promoted results independently re-derived and numerically verified

## 2026-02-19: Algebraic Exploration, Agent Findings (Experiment 12)

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
- |+⟩^N is an eigenstate of the isotropic Heisenberg Hamiltonian (energy variance = 0); it cannot generate entanglement or crossings. This corrects the hypothesis from Experiment 10.
- |0+0+⟩ (energy variance = 20) generates crossings from zero initial entanglement
- Unitary evolution: all 6 pairs cross; crossings oscillate reversibly
- With dephasing (gamma=0.05): only pair (0,2) crosses at t=0.285; all others miss the threshold (max C*Psi = 0.247)
- Dephasing survival is basis-dependent: σ_z spares |0⟩-qubits, kills |+⟩-qubits

**Updated:** SUBSYSTEM_CROSSING.md, Open Question 3 marked as answered
**Updated:** PREDICTIONS.md, 4 new Tier-2 entries (17 total computationally verified)
**Updated:** README.md, experiment list extended (Experiments 8-11)

## 2026-02-18: Subsystem Crossing Resolution (Experiment 10)

**Added:** SUBSYSTEM_CROSSING.md
- Crossing is local: operates at entangled qubit-pair level, not full system
- Bell+xBell+ N=4: pairs (0,1) and (2,3) cross at t=0.073 despite full-system Psi=0.200
- GHZ N=4: all pairs have l1=0, no pair-level crossing possible
- |+⟩^4: C=0 for all pairs permanently; high coherence without entanglement is inert
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
