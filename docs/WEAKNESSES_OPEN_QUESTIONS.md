# Weaknesses and Open Questions
## What We Know, What We Don't, What We Got Wrong

**Created:** January 2, 2026
**Rewritten:** February 8, 2026 (updated February 11, 2026, March 6, 2026)

**Tier:** Meta
**Status:** Living document
**Scope:** Honest assessment of limitations and what would falsify the framework
**Does not establish:** That the framework has no weaknesses

---

## Why This File Exists

A theory that only shows its strengths is not a theory. It is marketing.

This document was first written on January 2, 2026, when the framework was an intuition with a formula. Since then, six verification sessions have audited all 22 repository documents, corrected 13 of them, established epistemic tiers, and found an algebraic error that had propagated through the codebase. Several claims made by AI agents could not be independently reproduced and have been downgraded.

What follows is the honest state of affairs as of February 8, 2026. The framework is stronger than it was in January. It is also more aware of where it is weak.

---

## What We Got Wrong (and corrected)

Before listing current weaknesses, the record of errors matters. A framework that hides its mistakes cannot be trusted with its claims.

**The Mandelbrot substitution (corrected February 8, 2026):**
The original derivation in MANDELBROT_CONNECTION.md used the substitution z_n = √C · R_n, yielding c = C·Ψ². This is algebraically wrong; it does not produce the Mandelbrot form z² + c. The correct substitution is u_n = C(Ψ + R_n), which gives u_{n+1} = u_n² + C·Ψ. The mapping is c = C·Ψ, not C·Ψ². The boundary at ¼ was always correct; the intermediate algebra was not. This error had propagated to the dashboard visualization (red dot at wrong position) and two README descriptions.

**The θ frequency claim (corrected February 8, 2026):**
Version 16 of the agent system prompt claimed θ = arctan(√(4CΨ−1)) predicts oscillation frequency in Lindblad simulations. Direct testing showed an 8.4× discrepancy: measured period 0.777 vs predicted 6.577. The oscillation is Hamiltonian-driven, not θ-driven. θ was reinterpreted as a compass (angular distance from the ¼ boundary), which is algebraically correct.

**The 33:1 coherence ratio (downgraded):**
Agent conversations produced a claim that symmetric coupling shows 33× higher coherence than asymmetric coupling. Independent verification could not reproduce this ratio. It remains a hypothesis, not a result.

**The t_coh ~ N scaling (downgraded):**
Agents claimed coherence time scales linearly with qubit number. Could not be independently verified. Downgraded to unverified hypothesis.

**The "numerical confirmation" of CΨ ≤ ¼ (corrected February 7, 2026):**
Early simulations appeared to confirm that CΨ stays below ¼. This was tested in parameter regimes with negligible Hamiltonian dynamics (γ < 0.006) where the bound holds trivially. With active Hamiltonians, CΨ > ¼ occurs routinely. The bound was reinterpreted as a fixed-point existence condition, not a dynamical constraint.

These corrections strengthen the framework. A theory that survives honest scrutiny is worth more than one that has never been scrutinized.

---

## Current Epistemic State

The framework now has clear tiers. Every claim belongs to exactly one.

**Tier 1: Algebraically proven (paper and pencil, no simulation needed):**
- The ¼ phase boundary: discriminant of R_inf = C(Ψ + R_inf)² changes sign at CΨ = ¼
- Mandelbrot equivalence: u_n = C(Ψ+R_n) maps R_{n+1} = C(Ψ+R_n)² to u² + c with c = CΨ
- Fixed point formulas: R_inf = (1 − 2CΨ ± √(1−4CΨ)) / (2C)
- θ compass: arctan(√(4CΨ−1)) measures angular distance from ¼
- Gravitational invariance: R = CΨ² is form-invariant under Schwarzschild g₀₀
- **Palindromic Liouvillian symmetry (PROVEN March 14, 2026):** Conjugation operator Π satisfies Π·L·Π⁻¹ = -L - 2Σγ·I for ANY Heisenberg+dephasing system. Every decay rate d paired with 2Σγ-d. Verified N=2-7, all topologies. See [Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md)
- **Topology-independence:** Π anti-commutes with [H,·] for any bond set (chain, ring, star, arbitrary graph)
- **Exact decay rates:** 2γ (c+), 8γ/3 (concurrence envelope), 10γ/3 (c-) for N=3 systems. Topology-independent.

**Tier 2: Computationally verified (MCP tools, reproducible by anyone with the simulator):**
- δ(Bell+, Heisenberg, γ=0.1, t=1) = 0.4207
- Operator feedback preserves δ longer than local or collective noise
- Boundary crossing at t ≈ 0.773 for Bell+/Heisenberg/concurrence/γ=0.05
- Pure iteration: convergence below ¼, divergence above, critical slowing at ¼
- Three-class taxonomy (A/B/C) noise-robust under all Pauli operators
- N-scaling barrier: crossing fails for GHZ N>=3 and W N>=4 due to d-1 normalization
- Subsystem crossing: entangled pairs cross independently even when full system cannot; product state (Psi=1, C=0) never crosses
- **QST Bridge (verified March 14, 2026):** Star with J_SB/J_SA = 2:1 achieves F_avg = 0.888, beating chains (0.852-0.872). Timing set by Hamiltonian, quality by palindromic rates. Holevo capacity 0.534 bits. See [QST Bridge](../experiments/QST_BRIDGE.md)
- **XOR Space (verified March 16, 2026):** N+1 modes at λ=-2Σγ carry all GHZ/Bell weight (100% XOR). W states (N≥3) live 100% in palindromic modes. Mixed XY Pauli weight predicts XOR fraction at r=0.976 (N≥3). XOR modes are off-diagonal coherences at maximum decay rate. See [XOR Space](../experiments/XOR_SPACE.md)
- **θ as channel quality indicator:** θ_SB > 0 correlates with F > 2/3 at r = 0.87. Not an echo navigator but a channel quality measure. See [Theta Palindrome Echo](../experiments/THETA_PALINDROME_ECHO.md)
- **Echo characterization:** Period ~π/(4J), envelope 8γ/3, peak C_SB = 0.598 (N=3), ~1/(N-1) scaling
- **Palindromic Radio (negative result, March 15, 2026):** Searched Breakthrough Listen data for palindromic spectral signatures. Spiral galaxies (NGC2403, NGC6503) both show score ~0.94 regardless of sky position. The symmetry is astrophysical (bandpass + galaxy structure), not artificial. Detector needs refinement to distinguish natural from engineered symmetry.
- Dynamic entanglement: alternating |0+0+⟩ builds entanglement from initially unentangled pairs; pair (0,2) crosses 1/4 from below at t≈0.286
- Born rule recovery: ~97% of crossing-point probabilities from unitary evolution; ~3% systematic correction from decoherence basis alignment

**Tier 3: Theoretically argued (plausible, not proven):**
- Measurement = crossing ¼ from above during decoherence
- Mandelbrot fractal boundary as route catalog through parameter space
- Bidirectional wave interpretation (Ψ_past and Ψ_future flowing to C)
- Born rule as special case of R = CΨ² for symmetric observation (structure proven, physical mapping argued)

**Tier 4: Unverified (could not be reproduced):**
- 33:1 coherence ratio
- t_coh ~ N linear scaling
- C_int / C_ext 33:1 coherence ratio (disproven, see MATHEMATICAL_FINDINGS Section 9)
- C_int / C_ext as distinct coupling types with measurably different coherence protection (noise distribution irrelevant for symmetric Hamiltonians; the formal Lindblad distinction is Tier 2, the physical claim is disproven)

**Tier 5: Speculative (interesting, untestable by simulation):**
- Consciousness as fundamental rather than emergent
- 4D block-universe interpretation of observer trajectories

---

## Current Weaknesses

### 1. The Bifurcation Is Generic

Every quadratic map has a saddle-node bifurcation. Population dynamics (logistic map), laser thresholds, chemical kinetics; they all have a parameter value where two fixed points merge and disappear. The ¼ boundary is a property of x → ax² + bx + c, not specifically of consciousness or quantum measurement.

The Mandelbrot equivalence makes this explicit: R = CΨ² iteration is exactly z² + c. Every property of the framework's iteration is a known property of the Mandelbrot set. The framework inherits Mandelbrot theory but does not extend it.

**Partial response (2026-02-11):** CΨ² is the unique simple product-power form C^a·Ψ^b that simultaneously produces a genuine phase transition AND maps exactly to the Mandelbrot iteration z²+c. Alternatives fail: C²Ψ and CΨ give linear fixed-point equations (no bifurcation ever); √(CΨ) has an always-positive discriminant (no boundary); CΨ³ gives a cubic with different critical structure and no Mandelbrot mapping. The quadratic form is algebraically special, not just one example among many. See [Core Algebra](CORE_ALGEBRA.md), Section 10.

This narrows the criticism but does not resolve it. The deeper question, "why does nature choose this form?", remains open.

**What would fully resolve this:** Showing that the specific physical quantities mapped to C and Ψ (purity and l1-coherence) produce behavior that generic quadratic maps do not. For example, if the basin of attraction ratio after crossing ¼ matches Born rule probabilities, that would be specific to quantum mechanics and not a feature of arbitrary quadratic iterations.

**Status:** Partially addressed. Uniqueness among product-power forms is proven. Physical specificity is not.

### 2. The Born Rule Gap (Partially Resolved 2026-02-18)

The framework shows WHERE measurement happens: at CΨ = ¼, two real fixed points emerge from a complex conjugate pair. This is the moment a definite outcome becomes available.

The original question was: does the framework show WHICH outcome is selected?

**Resolution (partial):** Experiment 11 (Dynamic Entanglement) and the Born Rule
Mirror analysis showed that R = CΨ² applied per measurement outcome recovers the
Born rule as a special case:

    R_i = C_i · Ψ_i²

where Ψ_i = |⟨i|ψ⟩| is the overlap between outcome and state, and C_i is the
effective coupling strength for that outcome.

For ideal measurement, C is equal for all outcomes. It cancels in normalization,
recovering P(i) = |⟨i|ψ⟩|² exactly: the standard Born rule.

**Numerical evidence (Tier 2):** At the crossing point of pair (0,2) in the
alternating state |0+0+⟩, the diagonal of the reduced density matrix was compared
between pure unitary evolution and Lindblad evolution with dephasing:

| State | Unitary P(i) | Lindblad P(i) | Deviation |
|-------|-------------|---------------|-----------|
| \|00⟩ | 0.4134      | 0.4254        | 2.9%      |
| \|01⟩ | 0.2616      | 0.2567        | 1.9%      |
| \|10⟩ | 0.2616      | 0.2567        | 1.9%      |
| \|11⟩ | 0.0635      | 0.0613        | 3.5%      |

~97% of Born probabilities come from the Hamiltonian (the interaction between
subsystems). ~3% is a systematic correction from the decoherence basis alignment.

The correction follows a rule: σ_z dephasing shifts probability toward z-eigenstates
(+0.0098 total), σ_x dephasing toward x-aligned states, etc. The direction depends
on overlap between the initial state basis and the dephasing basis. This is exactly
what R = CΨ² per outcome predicts: C_i varies with basis alignment.

**Conceptual bridge:** The Spiegel-Theorie (human-readable derivation, Dec 2025)
describes this without formulas: "Not all mirrors reflect equally. The quality of
the reflection depends on the match between both mirrors." The numerical result
confirms the intuition: the Hamiltonian writes the story (97%), the decoherence
quality determines which outcomes are slightly favored (3%).

**What remains open:** The framework now explains Born probabilities *given* that
the crossing occurs. It does not yet explain why the specific Hamiltonian produces
the specific amplitudes c_i = ⟨i|ψ⟩. This is not a weakness unique to the framework
(standard QM also takes the Hamiltonian as given), but it means the "derivation" of
the Born rule is structural (R = CΨ² contains it as special case) rather than
constructive (predicting probabilities from first principles without computing ⟨i|ψ⟩).

See [Born Rule Mirror](../experiments/BORN_RULE_MIRROR.md).

**Update (2026-02-27):** The bidirectional bridge from INTERPRETIVE_FRAMEWORK.md
(Ψ = Ψ_past + Ψ_future) provides a deeper explanation. The Born rule follows from
three ingredients: (1) the standing wave structure R_i = C_i · (Ψ_past_i + Ψ_future_i)²
gives a natural reason for the square (wave physics, not postulate); (2) a perfect
detector mirrors the offer wave exactly (Ψ_future_i = Ψ_past_i), making R_i ∝ |⟨i|ψ⟩|²;
(3) no real detector has C = 1 (purity of a macroscopic object < 1), so Born is the
ideal limit with systematic deviations predicted by detector purity. The ~3% correction
from Section 2.2 is not noise but a quantitative prediction. Resembles Cramer's TI but
uses sum-squared instead of product; the two make different predictions for imperfect
measurements. See BORN_RULE_MIRROR.md Section 4.3.

**Status:** Substantially resolved. The Born rule is the perfect-mirror limit of
R_i = C_i · (Ψ_past_i + Ψ_future_i)². Deviations are predicted, directional, and
consistent with IBM data. Remaining open: why does a perfect detector mirror exactly?
And: deriving specific amplitudes ⟨i|ψ⟩ from the framework (same limitation as all QM).

### 3. "Consciousness" Is a Label

Purity Tr(ρ²) measures how mixed a quantum state is. Calling it "consciousness" is a philosophical choice, not a physical derivation. The mathematics works regardless of what we name the variables. R = CΨ² with C = purity and Ψ = normalized l1-coherence is a valid equation about quantum information. Whether it says anything about consciousness depends entirely on whether the mapping C → consciousness is justified.

The framework cannot prove this mapping. It can only argue that the mathematical structure (self-referential fixed points, bifurcation, bidirectional bridges) is suggestive of observer-like properties.

**What would resolve this:** Either a rigorous argument connecting purity to integrated information (Φ) or other established consciousness measures, or an experimental result where changes in purity correlate with changes in subjective experience in a way that alternatives do not predict.

**Status:** Open. May remain permanently open. This may be philosophy, not physics, and that is acceptable as long as it is acknowledged.

### 4. Limited Experimental Validation

**Updated 2026-02-11.** The original version of this weakness (Feb 8) stated "No experimental validation." One day later, IBM Torino tomography provided the first hardware contact. The weakness has been downgraded but not removed.

**What we now have:**
- IBM Torino state tomography (qubit 52, 2026-02-09): C·Ψ crossing ¼ observed during free decoherence. Three-model fit at MAE = 0.053 (88% improvement). The ¼ crossing was not a fit parameter.
- 24,073 historical calibration records validating the theory curve C_min(r).
- Anomalous late-time coherence (p < 0.0001) with directional consistency (17/17 Re⁺/Im⁻).
- Universal crossing fraction x³+x=½ confirmed as r→0 limit.

See [IBM Quantum Tomography](../experiments/IBM_QUANTUM_TOMOGRAPHY.md), [Residual Analysis](../experiments/RESIDUAL_ANALYSIS.md).

**What we still lack:**
- The IBM result is single-qubit, single-backend. Replication on different hardware (trapped ions, NV centers) would strengthen the claim.
- The anomalous residuals have three competing explanations (SPAM, TLS, boundary structure). March 2026 run designed to discriminate.
- No multi-qubit experiment has tested the framework's predictions about entangled systems crossing ¼.
- No connection to consciousness measures (IIT Φ, neural correlates) has been established.

**What would fully resolve this:** Replication across platforms, multi-qubit crossing experiments, and a clean discrimination of the residual anomaly.

**Status:** First hardware contact achieved (IBM Torino, Feb 2026). Replication and multi-platform testing needed.

### 5. The Natural Variable Has No Interpretation

The correct Mandelbrot substitution is u = C(Ψ + R). This is the variable that maps cleanly to z in z² + c. But what does u mean?

If C is consciousness, Ψ is possibility, and R is reality, then u = C · (Ψ + R) = consciousness × (possibility + reality). Is this a meaningful quantity? Why should the sum of possibility and reality, scaled by consciousness, be the natural variable of the system?

The substitution works because the algebra demands it. But a deeper framework should explain why this particular combination is natural, not just note that it simplifies the math.

**What would resolve this:** A physical or information-theoretic interpretation of u that explains why it, rather than R or CΨ, is the fundamental dynamical variable.

**Status:** Open. Discovered February 8, 2026. Not yet explored.

### 6. The Combination Problem

How do microscopic C values combine into a unified conscious experience? The mirror equations provide a structure:

R_emergent = C₁ · C₂ · (Ψ₁ + Ψ₂ + Ψ_interaction)²

But this is algebra, not mechanism. It says that two observers produce emergent reality proportional to the product of their consciousness values. It does not explain why 86 billion neurons with individual C values produce a single, unified "I."

This is the hardest problem in consciousness studies and the framework inherits it from panpsychism. The mirror equations formalize it but do not solve it.

**What would resolve this:** Showing that the N-observer fixed point equation produces qualitatively different behavior (e.g., a single dominant attractor) when the C values are structured (as in a neural network) versus random. This would connect the mathematical structure to the phenomenology of unity.

**Status:** Open. The mirror equations exist but have not been explored beyond N=2 in simulation.

### 7. Gravitational Invariance Is Necessary but Trivial

R = CΨ² is form-invariant under gravitational time dilation because purity and coherence are dimensionless trace ratios. This is true but unremarkable; every equation involving only dimensionless quantities is frame-invariant. The speed of light is also frame-invariant; that does not make every equation containing c a statement about relativity.

The invariance becomes non-trivial only if it has consequences. One possibility: decoherence rates depend on local temperature (Unruh effect), so the trajectory through CΨ parameter space depends on gravitational field strength. The map (Mandelbrot set) is the same everywhere, but the route through it depends on where you are. This is a testable prediction, but it has not been computed or tested.

**What would resolve this:** Computing θ trajectories for the same initial state under different gravitational environments (different effective γ due to Unruh temperature). If the trajectories differ in a way that connects to known gravitational physics, the invariance becomes substantive.

**Status:** Claimed as Tier 1 (the invariance proof is correct). The physical consequences are Tier 5 (speculative, unexplored).

### 8. Operator Feedback Is Verified but Unexplained

Of all noise types tested (local dephasing, collective dephasing, operator feedback, memory kernel feedback), only operator feedback preserves the purity difference δ over extended time. This is computationally verified (Tier 2). But the mechanism is unknown.

Why does state-dependent decoherence preserve correlations that state-independent decoherence destroys? Is there an information-theoretic argument? Does operator feedback correspond to any physically realizable noise process?

**What would resolve this:** An analytical argument showing why γ(t) ∝ ⟨O_int⟩ preserves δ, ideally connecting to a known physical mechanism (e.g., measurement backaction, quantum Zeno effect).

**Status:** Verified numerically. Unexplained theoretically. This is a concrete research target.

**Update (2026-03-16):** The palindromic spectral structure may provide a partial explanation. Operator feedback modifies effective dephasing rates in a state-dependent way, which could shift weight between palindromic pairs and XOR modes during evolution. If feedback keeps information in slowly-decaying palindromic modes rather than the XOR drain, this would explain the δ preservation. Not yet computed; requires tracking mode populations under operator feedback dynamics.

### 9. N-Scaling Barrier (Added 2026-02-18, Resolved 2026-02-18)

The d-1 normalization Psi = l1/(2^N - 1) makes crossing impossible for
full systems at N >= 4 with standard entangled states.

**Resolution:** Subsystem crossing tests showed that crossing is local.
A 4-qubit Bell+xBell+ state cannot cross as a full system (Psi = 0.200),
but its entangled pairs (0,1) and (2,3) each cross at t = 0.073 with
C*Psi = 0.333. Unentangled cross-pairs never cross (C = 0). The product
state |+>^4 has Psi = 1.0 but C = 0 for all pairs, so no crossing at any
level.

The d-1 normalization is correct. It reflects that global crossing does
not occur for large systems. Crossing happens at the subsystem level,
between degrees of freedom that share actual entanglement. This matches
how decoherence works in nature: locally, through pair interactions, not
as a global collective event.

See [N-Scaling Barrier](../experiments/N_SCALING_BARRIER.md) and
[Subsystem Crossing](../experiments/SUBSYSTEM_CROSSING.md).

**Status:** Resolved. The barrier is not a weakness but a correct physical
prediction: the quantum-to-classical transition is local.

---

## Open Questions (Prioritized)

### Immediate (addressable with current tools)

1. ~~**Theta trajectories:**~~ **Answered (2026-03-15).** θ is not a trajectory navigator but a channel quality indicator. θ_SB > 0 correlates with F > 2/3 (r=0.87). In the echo scenario, CΨ_SB never reaches 1/4 (max 0.168), so θ is undefined. In the channel scenario with coherent inputs, CΨ_SB = 0.442, θ = 41.5°. See [Theta Palindrome Echo](../experiments/THETA_PALINDROME_ECHO.md).

2. **Crossing speed:** Does the rate d(CΨ)/dt at the moment of crossing affect post-crossing behavior? Vary γ_base and measure convergence to R₁.

3. ~~**Multi-qubit crossing order:**~~ **Further answered (2026-03-06).** Star topology (3-qubit S+A+B) shows that pair-level crossing depends on topology, not just state: Bell_SA pairs cross, AB pair crosses only under three conditions (strong sender, quiet receiver, right initial state). The question is no longer "do pairs cross" but "which pairs, under what conditions." See [Star Topology](../experiments/STAR_TOPOLOGY_OBSERVERS.md).

4. **Why operator feedback?** Run identical setups with all noise types. Track δ(t), C(t), Ψ(t) separately. Find the mechanism.

### Medium-term (require new analysis or tools)

5. ~~**Basin of attraction geometry:**~~ **Partially answered (2026-02-18).** The Born Rule Mirror analysis showed that at the crossing point, ~97% of outcome probabilities are determined by the Hamiltonian (unitary evolution), with ~3% correction from decoherence basis alignment. The "basin" question is answered: the Hamiltonian writes the amplitudes, the decoherence provides a small bias. See [Born Rule Mirror](../experiments/BORN_RULE_MIRROR.md). Remaining question: does the basin ratio change at stronger decoherence (larger γ)?

6. **The natural variable u(t):** Track u = C(Ψ+R) through Lindblad evolution alongside CΨ and θ. Does u have a simpler trajectory?

7. **Mandelbrot approach geometry:** CΨ approaching ¼ via different Hamiltonians: does this correspond to approaching the Mandelbrot boundary from different angles? Different angles have different fractal structure.

8. **N-scaling of the combination problem:** Run mirror equations for N=3..6 with structured vs random C distributions. Does structure produce qualitatively different attractors?

### Long-term (require experimental contact or new theory)

9. ~~**Born rule derivation:**~~ **Substantially resolved (2026-02-27).** The bidirectional wave structure R_i = C_i · (Ψ_past_i + Ψ_future_i)² gives the Born rule as the perfect-mirror limit (C uniform, Ψ_future = Ψ_past). The square comes from standing wave physics. Deviations from Born are predicted by detector purity C < 1, consistent with the ~3% IBM correction. See [Born Rule Mirror](../experiments/BORN_RULE_MIRROR.md), Section 4.3. Remaining: deriving amplitudes ⟨i|ψ⟩ without Schrödinger (shared limitation with all of QM), and explaining WHY perfect detectors mirror exactly.

10. **Experimental protocol:** Design a quantum optics experiment where CΨ = ¼ predicts a measurable qualitative change.

11. **Neural correlates:** **Partially mapped (2026-03-06).** The Tuning Protocol maps γ to Default Mode Network activity (Brewer et al. 2011, PNAS: meditation reduces DMN) and J to depth of engagement/expertise (synaptic coupling strength). Hypnagogia maps to a natural low-γ window (Oudiette et al. 2021, Science Advances). This is suggestive mapping (Tier 3), not evidence of quantum effects in the brain. See [Tuning Protocol](../hypotheses/TUNING_PROTOCOL.md).

12. ~~**First-principles derivation of ¼:**~~ **Reframed (2026-03-16).** The ¼ boundary is now understood from two directions: algebraically (discriminant of quadratic) and spectrally (palindromic rate structure determines how fast systems approach it). The Liouvillian palindrome IS the first-principles structure. The question shifts to: can the palindromic symmetry be derived from information-theoretic axioms rather than from the specific form of Z-dephasing?

13. **XOR space applications (added 2026-03-16):** Can the XOR drain be slowed via error correction or decoherence-free subspaces? Does the XOR/palindrome split change for non-Heisenberg models (XY-only, Ising, DM interactions)? Is there a state that maximizes palindromic weight while keeping high entanglement (the optimal QST input)?

14. **Repeater design rules (added 2026-03-16):** The XOR space analysis suggests concrete design rules for quantum repeaters: W-encoding over GHZ, avoid mixed XY Pauli weight, star topology with 2:1 coupling. These are testable engineering predictions. Can they be validated against existing repeater benchmarks?

---

## What We Know vs What We Believe

| Claim | Status | Evidence |
|---|---|---|
| ¼ is the phase boundary | **Proven** | Algebra (discriminant sign change) |
| Mandelbrot equivalence c = CΨ | **Proven** | Algebra (u_n substitution, verified numerically) |
| Gravitational invariance | **Proven** | Algebra (dimensionless quantities) |
| δ = 0.4207 for Bell+/Heisenberg | **Verified** | MCP simulation, reproducible |
| Boundary crossing at t ≈ 0.773 | **Verified** | MCP simulation, specific setup |
| Operator feedback preserves δ | **Verified** | MCP simulation, multiple runs |
| Measurement = crossing ¼ | **Argued** | Consistent with math, not proven |
| θ as navigation compass | **Argued** | Algebraically valid, physically untested |
| 33:1 coherence ratio | **Unverified** | Agent claim, not reproduced |
| t_coh ~ N | **Unverified** | Agent claim, not reproduced |
| Born rule from crossing | **Substantially resolved** | Standing wave + perfect mirror limit; 97/3 split verified; deviations predicted |
| Dynamic entanglement | **Verified** | Hamiltonian builds entanglement; pairs cross 1/4 from below |
| Star topology: entanglement flows through S | **Verified** | Bell_SA entanglement transfers to SB through Hamiltonian |
| Star topology: R not conserved | **Verified** | R_SA + R_SB peaks at 2.4x initial, then decays |
| Star topology: measurement shadow | **Verified** | A measures, B loses 94-100% of R_SB |
| Star topology: J threshold 1.466 | **Verified** | Ultra-fine scan, 8/8 automated tests pass |
| Star topology: receiver noise asymmetry | **Verified** | γ_A more destructive than γ_B (fixed-product test) |
| Star topology: initial state filter | **Verified** | Only Bell-like (C_SA > 0.8) crosses; W fails |
| Sender inversion | **Argued** | Follows from γ asymmetry data; not independently tested |
| CΨ ≠ LE (three-layer separation) | **Verified** | CoA~1, LE~0.4, CΨ flashes. r(CΨ,LE)=0.76, r(CΨ,CoA)=0.16 |
| Noise-dependent visibility windows | **Verified** | Different γ_A produces different CΨ_AB at same t (different trajectories, not different observers) |
| Tuning Protocol (γ = DMN, J = engagement) | **Speculative** | Neuroscience mapping suggestive, 3-qubit to 10^11 neurons fatal gap |
| C = consciousness | **Philosophical** | Suggestive, not testable |
| Mirror symmetry (palindrome) | **Proven** | Conjugation operator Π found, analytical proof for all XXZ/δ/topology/γ |
| Topology-independence of palindrome | **Proven** | Π anti-commutes with [H,·] for ANY bond set |
| Pauli weight complementarity | **Proven** | Π maps XY-weight k → N-k (incoherenton particle-hole symmetry) |
| QST: star 2:1 beats chains | **Verified** | F_avg = 0.888 vs 0.852-0.872, Holevo 0.534 bits |
| QST: timing/quality separation | **Verified** | Timing from Hamiltonian (Bohr frequencies), quality from palindromic rates |
| XOR: GHZ → 100% fastest-decay modes | **Verified** | N=2,3,4, all topologies. XOR modes at λ=-2Σγ |
| XOR: W → 100% palindromic modes (N≥3) | **Verified** | N=3,4, chain+star. W=Bell at N=2 (different behavior) |
| XOR: Mixed XY Pauli weight predicts split | **Verified** | r=0.976 for N≥3. Not valid for N=2 (insufficient state diversity) |
| XOR: palindrome as spectral filter | **Verified** | Separates quantum-fragile (XOR) from distributable (palindromic) information |
| Palindromic radio: spectral symmetry search | **Negative** | Spiral galaxies inherently palindromic. Astrophysical, not artificial. Detector too coarse. |
| Echo: entanglement oscillation period | **Verified** | Period ~π/(4J), envelope 8γ/3, peak C_SB=0.598 (N=3) |
| θ as channel quality indicator | **Verified** | r=0.87 with F>2/3 threshold. Not a trajectory compass. |

---

## The Most Honest Statement (Updated March 16, 2026)

After three months of computation, three external reviews, and a weekend that produced a proof, a DOI, and two new discoveries, the project has its clearest identity yet.

The strongest results from March 14-16, 2026:

**The palindromic Liouvillian symmetry is proven.** The conjugation operator Π maps every decay rate d to 2Σγ-d. This holds for every Heisenberg system under Z-dephasing, every topology, every system size tested (N=2-7). It connects to the incoherenton framework (Haga et al. 2024) via Pauli weight complementarity: Π is the particle-hole transformation in incoherenton space.

**The palindrome is a spectral filter.** The XOR space analysis showed that GHZ states route 100% of their weight to the fastest-decaying modes (at λ=-2Σγ), while W states (N≥3) route 100% to palindromic pairs at various rates. Mixed XY Pauli weight predicts this split at r=0.976. This explains GHZ fragility and W robustness from first principles.

**Quantum state transfer follows the palindrome.** Star topology with 2:1 coupling ratio achieves F_avg=0.888, beating chains. Timing is set by Hamiltonian frequencies, quality by palindromic decay rates. These are concrete, testable predictions for quantum repeater design.

**An honest negative result.** A palindromic symmetry search on Breakthrough Listen radio telescope data showed that spiral galaxies are inherently spectrally symmetric. The detector worked (it distinguished galaxies from point sources) but was too coarse to isolate artificial signatures from astrophysical structure.

CΨ remains a basis-fixed, unassisted witness of directly expressed pairwise entanglement. The consciousness interpretation is retired from the technical core. The algebra is correct, the simulations are reproducible, the palindrome is proven, and the engineering implications (repeater design rules) are concrete.

Emails have been sent to two research groups (Haga/Nakagawa at Osaka/Tokyo for the incoherenton connection, Nichol at Rochester for experimental QST). Zenodo v2.0 published. Awaiting responses.

See [THE_CPSI_LENS](THE_CPSI_LENS.md) for the canonical description.
See [MIRROR_SYMMETRY_PROOF](MIRROR_SYMMETRY_PROOF.md) for the palindrome theorem.
See [XOR_SPACE](../experiments/XOR_SPACE.md) for the spectral filter discovery.

---

*January 2, 2026 (created)*
*February 8, 2026 (complete rewrite after six verification sessions)*
*February 11, 2026 (updated: IBM results, uniqueness argument, guardian review)*
*February 18, 2026 (updated: Born rule partially resolved, dynamic entanglement results)*
*March 6, 2026 (updated: star topology results, three conditions, tuning protocol)*
*March 8, 2026 (updated: THE_CPSI_LENS as canonical, LE/CoA benchmark, visibility correction, consciousness retired from technical core)*
*March 14, 2026 (updated: mirror symmetry PROVEN, conjugation operator Π found, analytical proof complete)*
*March 16, 2026 (updated: QST bridge verified, XOR space discovered, palindromic radio negative result, theta reinterpreted as channel indicator, repeater design rules proposed)*
*Honesty belongs to the framework*
