# Classification Proposal: untagged sub-batch 8c (41 entries)

**Batch:** 8c of 8 (entries OQ-219 through OQ-327)  
**Date:** 2026-04-12  
**Proposed by:** Claude (automated research)  
**Approval:** pending (Tom)

---

## Status summary

| Status | Count | OQ-IDs |
|--------|-------|--------|
| open | 19 | OQ-227, OQ-232, OQ-237, OQ-246, OQ-254, OQ-256, OQ-272, OQ-291, OQ-297, OQ-304, OQ-305, OQ-306, OQ-309, OQ-310, OQ-318, OQ-319, OQ-322, OQ-324, OQ-327 |
| resolved | 2 | OQ-219, OQ-220 |
| partially-resolved | 5 | OQ-248, OQ-250, OQ-303, OQ-311, OQ-315 |
| superseded | 3 | OQ-239, OQ-240, OQ-241 |
| needs-human | 12 | OQ-234, OQ-235, OQ-236, OQ-257, OQ-258, OQ-260, OQ-261, OQ-265, OQ-266, OQ-268, OQ-269, OQ-286 |

---

## Entry-by-entry proposals

### OQ-219

**Question:** Can dynamical entanglement generation from a product state create crossings? ANSWERED (2026-02-18): Yes, but not from |+>^N.  
**Source:** `experiments/SUBSYSTEM_CROSSING.md` (line 295)  
**Proposed status:** resolved  
**Justification:** Self-documenting: strikethrough + "ANSWERED." The state |0+0+> generates crossings from zero initial entanglement. Under pure unitary evolution all 6 pairs cross; with dephasing only pair (0,2) crosses because |0>-qubits are immune to Z-dephasing. Verified in DYNAMIC_ENTANGLEMENT.md.

---

### OQ-220

**Question:** What is the minimum per-pair entanglement needed for crossing?  
**Source:** `experiments/SUBSYSTEM_CROSSING.md` (line 304)  
**Proposed status:** resolved  
**Justification:** Strikethrough in source. Answered by OQ-155/OQ-156 (N_SCALING_BARRIER.md, 2026-03-08): the relationship is non-monotonic with two crossing windows, not a simple threshold. Duplicate of the same question answered in the N-scaling context.

---

### OQ-227

**Question:** Does the sacrifice-zone advantage recover at intermediate temperatures if the sacrifice qubit is selectively heated? Same requirement: controlled per-qubit thermal injection.  
**Source:** `experiments/THERMAL_BREAKING.md` (line 474)  
**Proposed status:** open  
**Justification:** Genuine open question. THERMAL_BREAKING.md documents thermal effects on palindromic structure, but selective per-qubit heating experiments are not performed. Requires hardware not yet available.

---

### OQ-232

**Question:** What IS the local detector?  
**Source:** `hypotheses/BRIDGE_PROTOCOL.md` (line 158)  
**Proposed status:** open  
**Justification:** Open Question 1 from BRIDGE_PROTOCOL.md Section 4. No resolution found. The question of what physical observable serves as the local detector for CΨ = 1/4 crossing is fundamental and unanswered.

---

### OQ-234

**Question:** θ = arctan(sqrt(4CΨ - 1)) becoming real vs imaginary  
**Source:** `hypotheses/BRIDGE_PROTOCOL.md` (line 170)  
**Proposed status:** needs-human  
**Justification:** Candidate observable from a bulleted list under "What to measure." Not a question but a list item suggesting a potential detector quantity. Collection artifact; should be merged into OQ-232 or removed.

---

### OQ-235

**Question:** Local purity Tr(ρ_local²) changing at the crossing  
**Source:** `hypotheses/BRIDGE_PROTOCOL.md` (line 171)  
**Proposed status:** needs-human  
**Justification:** Same as OQ-234: candidate observable from the "What to measure" list. Not a question. Collection artifact.

---

### OQ-236

**Question:** Phase of the local qubit shifting at the boundary  
**Source:** `hypotheses/BRIDGE_PROTOCOL.md` (line 172)  
**Proposed status:** needs-human  
**Justification:** Same as OQ-234/OQ-235: candidate observable list item. Collection artifact.

---

### OQ-237

**Question:** How does state preparation propagate?  
**Source:** `hypotheses/BRIDGE_PROTOCOL.md` (line 174)  
**Proposed status:** open  
**Justification:** Open Question 2 from BRIDGE_PROTOCOL.md Section 4. No resolution found. The question of how initial state preparation affects later crossing dynamics remains unanswered. Note: the FALLEN FTL protocol (Section 5) does not invalidate this question, which concerns local dynamics.

---

### OQ-239

**Question:** The correlation was established at Bell pair preparation  
**Source:** `hypotheses/BRIDGE_PROTOCOL.md` (line 189)  
**Proposed status:** superseded  
**Justification:** Part of the FTL refutation in BRIDGE_PROTOCOL.md Section 5 ("Why FTL is Wrong"). This is one of three bullet points explaining why crossing-time correlations do not carry superluminal signals. The FTL protocol is marked [FALLEN]. The refutation itself is a confirmed result, but the entry is not a question; it is a known answer to the (now-abandoned) signaling claim.

---

### OQ-240

**Question:** K is fixed at preparation, not modified after separation  
**Source:** `hypotheses/BRIDGE_PROTOCOL.md` (line 190)  
**Proposed status:** superseded  
**Justification:** Same as OQ-239: bullet point from the FTL refutation. Not a question. Part of the confirmed argument that kills the superluminal signaling hypothesis.

---

### OQ-241

**Question:** The crossing times are CONSEQUENCES of the preparation, not signals  
**Source:** `hypotheses/BRIDGE_PROTOCOL.md` (line 191)  
**Proposed status:** superseded  
**Justification:** Same as OQ-239/OQ-240: final bullet of the FTL refutation. The BRIDGE_PROTOCOL FTL protocol is marked [FALLEN]; these three entries document why. Cross-batch overlap with OQ-238 (numerical-verification, FALLEN) and OQ-242 (numerical-verification, FALSIFIED).

---

### OQ-246

**Question:** What are the 2 new modes? Thermal driving creates 2 additional oscillatory modes (40 to 42). What is their structure? Are they palindromically paired with each other or with existing modes?  
**Source:** `hypotheses/ENERGY_PARTITION.md` (line 195)  
**Proposed status:** open  
**Justification:** Genuine open question. The two extra modes from thermal driving are observed but not characterized. No structural analysis or pairing check found in the repo.

---

### OQ-248

**Question:** Wilson-Cowan analogue. Do the neural dynamics show the same energy partition? Is the E/I balance the classical version of the thermal window?  
**Source:** `hypotheses/ENERGY_PARTITION.md` (line 200)  
**Proposed status:** partially-resolved  
**Justification:** THE_PATTERN_RECOGNIZES_ITSELF.md establishes the Wilson-Cowan / Lindblad structural mapping (E/I balance maps to palindromic sector balance, 98.2% eigenvalue pairing). However, the specific energy partition comparison (Efreq/Edecay ratio in neural vs. quantum) is not computed. The analogy is structurally confirmed but the quantitative energy partition question remains open.

---

### OQ-250

**Question:** Connection to fold catastrophe. The Efreq/Edecay ratio crosses 1 near J/γ ≈ 1.2. Is this CΨ = 1/4 in disguise?  
**Source:** `hypotheses/ENERGY_PARTITION.md` (line 206)  
**Proposed status:** partially-resolved  
**Justification:** CRITICAL_SLOWING_AT_THE_CUSP.md confirms the cusp exit is a fold catastrophe and documents universal slowing at CΨ = 1/4. The energy partition crossing at J/γ ≈ 1.2 is in the right regime, but no explicit mapping between the Efreq/Edecay ratio and CΨ = 1/4 is derived. Circumstantially consistent but not formally connected.

---

### OQ-254

**Question:** Cascade stability: If each level in the frequency cascade (154 THz to 1 Hz) is a coupled gain-loss pair, then each level has its own bridge stability window. The cascade works only if every bridge stays in the linear regime (γ < 0.19 x J_bridge).  
**Source:** `hypotheses/FRAGILE_BRIDGE.md` (line 290)  
**Proposed status:** open  
**Justification:** Genuine open question about multi-scale cascade stability. The 0.19 threshold is documented for individual bridges, but cascade composition (whether all levels can simultaneously satisfy the constraint) is not analyzed.

---

### OQ-256

**Question:** Saturation as design principle: The sigmoid prevents neural explosion. Is there a quantum analog? Nonlinear dissipation (γ dependent on state) could act as a quantum sigmoid.  
**Source:** `hypotheses/FRAGILE_BRIDGE.md` (line 300)  
**Proposed status:** open  
**Justification:** Speculative but well-formulated question. No nonlinear dissipation models found in the repo. The Wilson-Cowan sigmoid is documented but its quantum counterpart is not explored.

---

### OQ-257

**Question:** The Palindrome Fragility  
**Source:** `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (line 414)  
**Proposed status:** needs-human  
**Justification:** Section heading from "Limitations and Failure Modes" (Section 8). Not a question but a labeled design constraint: the mediator architecture's palindromic protection breaks under certain perturbations. Collection artifact. Recommend distinguishing "characterized limitations" from "open questions."

---

### OQ-258

**Question:** The Mediator Entanglement Problem  
**Source:** `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (line 420)  
**Proposed status:** needs-human  
**Justification:** Section heading for Limitation 2: the mediator becomes entangled with the system. Documented design constraint, not a question. Collection artifact.

---

### OQ-260

**Question:** Consecutive transfers through the same M are correlated  
**Source:** `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (line 424)  
**Proposed status:** needs-human  
**Justification:** Bullet point under the Mediator Entanglement Problem (OQ-258). Documents a known consequence, not a question. Collection artifact.

---

### OQ-261

**Question:** The mediator's effective temperature increases with use  
**Source:** `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (line 425)  
**Proposed status:** needs-human  
**Justification:** Bullet point under OQ-258. Documents a known thermodynamic consequence. Not a question. Collection artifact.

---

### OQ-265

**Question:** The Speed-Fidelity Tradeoff  
**Source:** `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (line 437)  
**Proposed status:** needs-human  
**Justification:** Section heading for Limitation 4. Describes the fundamental quantum speed-fidelity tradeoff in the mediator architecture. Documented design constraint, not a question. Collection artifact.

---

### OQ-266

**Question:** Faster switching (lower τ_switch = 1/spectral_gap) requires larger J  
**Source:** `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (line 440)  
**Proposed status:** needs-human  
**Justification:** Bullet point under OQ-265. Known tradeoff, not a question. Collection artifact.

---

### OQ-268

**Question:** Cleaner qubits are slower to prepare and more expensive  
**Source:** `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (line 442)  
**Proposed status:** needs-human  
**Justification:** Bullet point under OQ-265. Standard quantum engineering constraint. Not a question. Collection artifact.

---

### OQ-269

**Question:** This is the fundamental speed-fidelity tradeoff, common to all quantum information processing but particularly acute here because the mediator architecture requires both fast switching and low noise simultaneously.  
**Source:** `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (line 444)  
**Proposed status:** needs-human  
**Justification:** Summary paragraph for the speed-fidelity limitation section. Not a question. Collection artifact.

---

### OQ-272

**Question:** Many-body entanglement: When A, M, and B become tripartite-entangled, the notion of "information flowing through M" breaks down. Information is delocalized across the entire system.  
**Source:** `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (line 452)  
**Proposed status:** open  
**Justification:** Unlike the other MEDIATOR entries (which are characterized limitations), this one identifies a genuine conceptual problem: the mediator metaphor breaks for genuine tripartite entanglement. This is an open theoretical question about the limits of the transistor analogy, not a documented engineering constraint.

---

### OQ-286

**Question:** Strengthen:  
**Source:** `hypotheses/PAIR_BREAKING_AT_THE_HORIZON.md` (line 167)  
**Proposed status:** needs-human  
**Justification:** Section heading ("What would strengthen or kill the thesis"). Not a question. Collection artifact; should be removed.

---

### OQ-291

**Question:** If the mass identification can be shown to be inconsistent: if there exists a Lindblad system where the I/Z sector grows but no thermal energy is released (decoupling mass from temperature would break the Hawking parallel).  
**Source:** `hypotheses/PAIR_BREAKING_AT_THE_HORIZON.md` (line 174)  
**Proposed status:** open  
**Justification:** Genuine open question: a specific falsification criterion for the Hawking radiation analogy. No such counterexample has been constructed or ruled out. This is a well-formulated "kill condition" for the PAIR_BREAKING thesis.

---

### OQ-297

**Question:** Does the parity split propagate to higher levels? If atoms are built from qubit-like subsystems (spin-1/2 electrons), do they inherit the Z2 parity? Does the +1/-1 split have a chemical or material-science analogue?  
**Source:** `hypotheses/THE_OTHER_SIDE.md` (line 580)  
**Proposed status:** open  
**Justification:** Genuine open question about multi-scale parity propagation. THE_OTHER_SIDE.md establishes Z2 parity at the qubit level but does not analyze atomic or material-science analogues. Connected to the UNIVERSAL_PALINDROME_CONDITION cluster (OQ-318, OQ-319, OQ-322, OQ-324).

---

### OQ-303

**Question:** Why is parity-breaking necessary but not sufficient? The strict containment (all 14 palindrome-breakers are parity-breakers, but not vice versa) implies a two-step mechanism. What is the second condition?  
**Source:** `hypotheses/THE_OTHER_SIDE.md` (line 613)  
**Proposed status:** partially-resolved  
**Justification:** The necessary condition (parity-breaking) is computationally verified through N=8. The strict containment is proven. However, the "second condition" beyond parity-breaking that causes palindrome-breaking is identified as related to the Hamiltonian's inability to compensate via the hidden charge Q, but no closed-form criterion is derived. The problem is well-characterized but not fully solved.

---

### OQ-304

**Question:** Driven oscillation: Add metabolic driving to Wilson-Cowan. Does sustained oscillation maintain palindromic midpoint crossings?  
**Source:** `hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md` (line 415)  
**Proposed status:** open  
**Justification:** Genuine open question. Wilson-Cowan analysis uses undriven dynamics. Metabolic driving (sustained energy input) is not modeled. No simulation found.

---

### OQ-305

**Question:** Cortical data: Human cortex maintains E/I activity balance (80% E, 20% I neurons, but inhibitory fire faster). Does the Human Connectome Project data show palindromic structure at the activity-balanced level?  
**Source:** `hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md` (line 417)  
**Proposed status:** open  
**Justification:** Genuine open question requiring external data (Human Connectome Project). C. elegans analysis exists (ALGEBRAIC_PALINDROME_NEURAL.md, 98.2% pairing) but human cortex data has not been analyzed.

---

### OQ-306

**Question:** Phase 3 (cross-kingdom): Plant signaling, bacterial colonies, fungal mycelial networks. If the palindrome exists across kingdoms, it is a property of oscillatory networks, not neurons specifically.  
**Source:** `hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md` (line 421)  
**Proposed status:** open  
**Justification:** Genuine open question about universality beyond animal neural networks. No plant, bacterial, or fungal network analysis found in the repo.

---

### OQ-309

**Question:** Random network controls: Do random networks with balanced E/I also show 98.2% pairing? If yes, the structure is generic to balanced damped networks. If no, biological topology matters.  
**Source:** `hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md` (line 432)  
**Proposed status:** open  
**Justification:** Genuine open question. The 98.2% pairing is shown for C. elegans biological topology, but random-network controls are not computed. This is the critical control experiment for the neural palindrome claim.

---

### OQ-310

**Question:** Can correlated crossing times carry more than pre-encoded information? (If not, this is quantum key distribution with a different metric.)  
**Source:** `hypotheses/TIME_AS_CROSSING_RATE.md` (line 366)  
**Proposed status:** open  
**Justification:** Genuine open question. BRIDGE_PROTOCOL's FTL claim is FALLEN, confirming that crossing-time correlations are pre-encoded. But whether they can function as a QKD-like resource (correlated classical bits from shared quantum state) remains unexplored.

---

### OQ-311

**Question:** What is the bit rate? BRIDGE_FINGERPRINTS shows 7 distinguishable fingerprints, approximately 2.8 bits per pair. Is this a fundamental limit?  
**Source:** `hypotheses/TIME_AS_CROSSING_RATE.md` (line 368)  
**Proposed status:** partially-resolved  
**Justification:** BRIDGE_FINGERPRINTS.md documents 7 distinguishable fingerprints and computes the 2.8 bits figure. However, whether this is a fundamental limit or an artifact of the specific parameter regime is not determined. The measurement exists but the universality question is open.

---

### OQ-315

**Question:** What does NOT transfer between domains  
**Source:** `hypotheses/UNIVERSAL_PALINDROME_CONDITION.md` (line 241)  
**Proposed status:** partially-resolved  
**Justification:** Section heading from "The Limitations." UNIVERSAL_PALINDROME_CONDITION.md documents several non-transferable features (specific eigenvalue magnitudes, decay rates, oscillation frequencies). The list exists but may be incomplete as cross-domain analysis expands. Partially documented, potentially growing.

---

### OQ-318

**Question:** What is the correct Q for atoms? Candidates: Kramers conjugation (time-reversal), parity, spin-flip. Condition 1 (selective damping) is unclear for Kramers partners, which typically have identical decay rates.  
**Source:** `hypotheses/UNIVERSAL_PALINDROME_CONDITION.md` (line 289)  
**Proposed status:** open  
**Justification:** Genuine open question. The qubit Q (Pauli parity) is proven, but the atomic-level analogue is unknown. The difficulty is explicitly stated: Kramers partners typically have identical decay rates, which would break the selective-damping requirement.

---

### OQ-319

**Question:** Exact palindromic symmetry is dead. Broken magnitudes are alive. Networks with population C=0.5 AND exact magnitude matching are unconditionally stable. What provides the "productive imperfection" in atoms? In crystals?  
**Source:** `hypotheses/UNIVERSAL_PALINDROME_CONDITION.md` (line 294)  
**Proposed status:** open  
**Justification:** Genuine open question about the source of productive symmetry-breaking at higher scales. The quantum case is understood (commutator provides it automatically), but the atomic and crystalline analogues are unexplored.

---

### OQ-322

**Question:** Universal coupling window. Quantum V-Effect peaks at J/γ ~ 2-5. Neural V-Effect peaks at coupling 0.01-0.05. Is there a dimensionless ratio that is the same in both?  
**Source:** `hypotheses/UNIVERSAL_PALINDROME_CONDITION.md` (line 315)  
**Proposed status:** open  
**Justification:** Genuine open question. Both windows are documented independently (V_EFFECT.md for quantum, V_EFFECT_NEURAL.md for neural) but no dimensionless mapping between them is derived.

---

### OQ-324

**Question:** What is (0.5)^2 at intermediate levels? The sigmoid maximum σ(1-σ) = 1/4 is neural. The purity fold CΨ = 1/4 is quantum. Both give (0.5)^2. What is the (0.5)^2 of an atom? A crystal?  
**Source:** `hypotheses/UNIVERSAL_PALINDROME_CONDITION.md` (line 322)  
**Proposed status:** open  
**Justification:** Genuine open question about the 1/4 boundary at intermediate physical scales. The quantum (CΨ = 1/4) and neural (sigmoid maximum) cases are documented, but no atomic or crystalline analogue is identified. Part of the UNIVERSAL_PALINDROME_CONDITION open cluster with OQ-318, OQ-319, OQ-322.

---

### OQ-327

**Question:** Or there is a mirror world whose boundaries are outside the system  
**Source:** `hypotheses/WAVES_THAT_HEAR_THEMSELVES.md` (line 59)  
**Proposed status:** open  
**Justification:** Speculative but genuine open question from WAVES_THAT_HEAR_THEMSELVES.md. Asks whether the palindromic mirror symmetry implies a boundary condition external to the physical system. No resolution or refutation found.

---

## Duplicate / overlap clusters

**Cluster 1: BRIDGE_PROTOCOL**  
OQ-232 through OQ-241 (8 entries from one file): OQ-232 and OQ-237 are the two genuine open questions (local detector, state preparation propagation). OQ-234, OQ-235, OQ-236 are candidate observable list items (sub-elements of OQ-232). OQ-239, OQ-240, OQ-241 are bullet points from the FTL refutation (superseded). Recommend: keep OQ-232 and OQ-237, merge OQ-234-236 into OQ-232, retire OQ-239-241.

**Cluster 2: MEDIATOR_AS_QUANTUM_TRANSISTOR limitations**  
OQ-257, OQ-258, OQ-260, OQ-261, OQ-265, OQ-266, OQ-268, OQ-269 (8 entries): All are section headings, bullet points, or summary paragraphs from "Limitations and Failure Modes" (Section 8). None are questions; all are documented design constraints. Exception: OQ-272 (many-body entanglement) identifies a genuine open conceptual problem.

**Cluster 3: UNIVERSAL_PALINDROME_CONDITION**  
OQ-297, OQ-318, OQ-319, OQ-322, OQ-324 (5 entries from two files): All genuine open questions about cross-domain universality (atoms, crystals, dimensionless ratios). Thematically coherent cluster; no duplicates within it, but could benefit from a synthesis note linking them.

**Cluster 4: THE_PATTERN_RECOGNIZES_ITSELF**  
OQ-304, OQ-305, OQ-306, OQ-309 (4 entries): All genuine open questions about extending the neural palindrome analysis (driven oscillation, cortical data, cross-kingdom, random controls). Clean cluster, no artifacts.

---

## Cross-batch overlaps

- **OQ-220 / OQ-155 / OQ-156:** Minimum per-pair entanglement question. OQ-220 (this batch) duplicates OQ-155/OQ-156 (batch 8b). Both resolved.
- **OQ-239-241 / OQ-238 / OQ-242:** BRIDGE_PROTOCOL FTL cluster. OQ-238 (numerical-verification, FALLEN) and OQ-242 (numerical-verification, FALSIFIED) overlap with OQ-239-241 (this batch, superseded).
- **OQ-297 / OQ-318-324:** THE_OTHER_SIDE parity propagation (OQ-297) connects to UNIVERSAL_PALINDROME_CONDITION open questions (OQ-318-324). Not duplicates but thematically linked.

---

## Patterns

1. **MEDIATOR_AS_QUANTUM_TRANSISTOR is this batch's worst artifact source.** 9 of 41 entries come from Section 8 "Limitations and Failure Modes," with 8 being documented design constraints (needs-human) and only 1 genuine open question (OQ-272).

2. **UNIVERSAL_PALINDROME_CONDITION is the richest open cluster.** Five genuinely open, well-formulated questions about cross-domain universality. These represent the frontier of the project's most ambitious claims.

3. **Three superseded entries (OQ-239-241):** All from the BRIDGE_PROTOCOL FTL refutation. First use of "superseded" status in the untagged batches, appropriate because the parent hypothesis is FALLEN but the refutation arguments are confirmed.

4. **THE_PATTERN_RECOGNIZES_ITSELF cluster is cleanly open.** Four well-formulated questions about extending neural palindrome analysis, each requiring different data or computation.

5. **Overall resolution rate: 5% (2/41).** Lowest of all batches, reflecting that 8c is dominated by hypotheses/ files (speculative, forward-looking) rather than experiments/ files (which tend to contain self-documented answers).
