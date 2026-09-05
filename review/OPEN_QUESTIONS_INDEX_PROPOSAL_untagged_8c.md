# Classification Proposal: untagged sub-batch 8c (41 entries)

**Batch:** 8c of 8 (entries OQ-219 through OQ-327)  
**Date:** 2026-04-12  
**Proposed by:** Claude (automated research)  
**Approval:** pending (Tom)

---

## Status summary

| Status | Count | OQ-IDs |
|--------|-------|--------|
| open | 19 | OQ-227, OQ-232, OQ-237, OQ-246, OQ-248, OQ-250, OQ-254, OQ-256, OQ-272, OQ-291, OQ-297, OQ-304, OQ-305, OQ-306, OQ-310, OQ-318, OQ-322, OQ-324, OQ-327 |
| resolved | 3 | OQ-219, OQ-220, OQ-309 |
| partially-resolved | 2 | OQ-303, OQ-311 |
| superseded | 4 | OQ-239, OQ-240, OQ-241, OQ-319 |
| needs-human | 13 | OQ-234, OQ-235, OQ-236, OQ-257, OQ-258, OQ-260, OQ-261, OQ-265, OQ-266, OQ-268, OQ-269, OQ-286, OQ-315 |

---

## Entry-by-entry proposals

### OQ-219

**Question:** Can dynamical entanglement generation from a product state create crossings? ANSWERED (2026-02-18): Yes, but not from |+>^N.  
**Source:** `experiments/SUBSYSTEM_CROSSING.md` (line 295)  
**Proposed status:** resolved  
**Justification:** Self-documenting: strikethrough + "ANSWERED." The state |0+0+> generates crossings from zero initial entanglement. Under pure unitary evolution all 6 pairs cross; with dephasing only pair (0,2) crosses because |0>-qubits are immune to Z-dephasing. Verified in DYNAMIC_ENTANGLEMENT.md. *(2026-07-20: the (0,2)-crossing detail did not reproduce under the canonical pair-CΨ book; the YES survives via the chain and |+-+->/|0+0-> on the ring, see the reproduction note in DYNAMIC_ENTANGLEMENT.md.)*

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

**Question:** What are the two extra resolved quantum roots in the N=3 combined-bath sweep (40 to 42), and how do their invariant subspaces contribute to a specified response?\
**Source:** [hypotheses/ENERGY_PARTITION.md](../hypotheses/ENERGY_PARTITION.md), section "4. Open questions"\
**Proposed status:** open\
**Justification:** The thermal table is a finite eigenvalue census at J=1, dephasing γ=0.1 and thermal Γ=0.1. Tracking roots and readout contributions remains a separate task; spectral sums are not energies.

---

### OQ-248

**Question:** Can a specified neural generator support a comparable spectral or response diagnostic under controlled external input?\
**Source:** [hypotheses/ENERGY_PARTITION.md](../hypotheses/ENERGY_PARTITION.md), section "3. What a neural comparison would require"\
**Proposed status:** open\
**Justification:** E/I balance supplies no thermal calibration or biological F36 result. The [current neural report](../docs/neural/V_EFFECT_NEURAL.md) has synthetic coupling/drive censuses, not a neural energy-partition or 2× law. Specify units, a converged operating point, readout and matched controls before comparing.

---

### OQ-250

**Question:** Does a spectral frequency/decay-sum crossover coincide with state-dependent CΨ=1/4 on the same specified quantum trajectory?\
**Source:** [hypotheses/ENERGY_PARTITION.md](../hypotheses/ENERGY_PARTITION.md), section "4. Open questions"\
**Proposed status:** open\
**Justification:** The scalar recursion's algebraic fold supplies no identification with the spectral ratio. Define generator, state, time and normalization and evaluate the two independently; neither a nearby J/γ value nor the neural quarter arithmetic resolves this gate.

---

### OQ-254

**Question:** Can a specified sequence of quantum gain-loss bridges remain stable under joint coupling?\
**Source:** [hypotheses/FRAGILE_BRIDGE.md](../hypotheses/FRAGILE_BRIDGE.md), section "6. Open questions"\
**Proposed status:** open\
**Justification:** The N=2 weak-bridge fit is not a compositional criterion or a biological frequency-cascade model. Build the combined generator and compare its spectrum with isolated-bridge predictions.

---

### OQ-256

**Question:** Can an explicitly specified quantum gain model with saturation bound its dynamics?\
**Source:** [hypotheses/FRAGILE_BRIDGE.md](../hypotheses/FRAGILE_BRIDGE.md), section "6. Open questions"\
**Proposed status:** open\
**Justification:** The Wilson-Cowan invariant activity square follows from its bounded sigmoid and positive time constants. It establishes neither biological safety nor quantum stability. A candidate quantum model needs a physical state domain, generator and independent stability gate.

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

**Question:** For a specified driven neural model, can the F36 transport relation predict an observed response?\
**Source:** [hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md](../hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md), section "A research program that can reject the pattern"\
**Proposed status:** open\
**Justification:** External-drive censuses already exist in [the neural report](../docs/neural/V_EFFECT_NEURAL.md). P has no metabolic calibration. The remaining gate needs a converged operating point, full scalar identity, specified perturbation/readout and agreement between linear transport and the measured response.

---

### OQ-305

**Question:** Can a specified biological cortical circuit and operating point satisfy both F36 conditions?\
**Source:** [hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md](../hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md), section "The biological result is a support obstruction"\
**Proposed status:** open\
**Justification:** No biological neural network in the repository is known to pass F36. The full committed C. elegans chemical model fails its necessary support condition. Cortical activity balance or a connectome alone does not supply an effective Jacobian; test its support, diagonal pair sums and scaled magnitudes.

---

### OQ-306

**Question:** Can a specified non-neural biological generator realize the same scalar conjugation identity?\
**Source:** [hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md](../hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md), section "A research program that can reject the pattern"\
**Proposed status:** open\
**Justification:** This is a candidate-substrate search, not an established cross-kingdom law. Identify a generator, map Q and scalar shift, then reject failures entrywise. Oscillation by itself supplies none of those hypotheses.

---

### OQ-309

**Question:** Do matched random controls and the full scalar gate support the proposed biological pairing interpretation?\
**Source:** [hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md](../hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md), section "The biological result is a support obstruction"\
**Proposed status:** resolved\
**Justification:** The [current matched audit](../docs/neural/ALGEBRAIC_PALINDROME_NEURAL.md#3-empirical-c-elegans-null) supplies the support null and normalization/control results. A matching percentage does not establish biological F36, and degree-preserving controls can leave a coupling-norm diagnostic unchanged. The original comparison is not a missing-control confirmation task; a different specified biological candidate remains open.

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

**Question:** What does not transfer between domains?\
**Source:** [hypotheses/UNIVERSAL_PALINDROME_CONDITION.md](../hypotheses/UNIVERSAL_PALINDROME_CONDITION.md), section "The exact identity and its reach"\
**Proposed status:** needs-human\
**Justification:** This is a comparison heading, not a distinct open question. The current source separates conditional complex pairing and generalized-eigenspace transport from biological realization, dynamics, rate ratios and mechanism claims. Merge this collection artifact into the source's concrete gates.

---

### OQ-318

**Question:** What generator and linear conjugation could realize the scalar identity in a specified atomic model?\
**Source:** [hypotheses/UNIVERSAL_PALINDROME_CONDITION.md](../hypotheses/UNIVERSAL_PALINDROME_CONDITION.md), section "The hypothesis and the next gates"\
**Proposed status:** open\
**Justification:** A named candidate such as parity, spin-flip or time reversal is insufficient. Specify whether the map is linear, its domain and the full identity. In the involutive permutation setting equal diagonal rates are not an obstruction: all paired diagonal sums must equal one scalar −2s, and effective coupling must separately reverse.

---

### OQ-319

**Question:** Does exact palindromic symmetry force silence or stability, making broken magnitudes necessary for activity?\
**Source:** [hypotheses/UNIVERSAL_PALINDROME_CONDITION.md](../hypotheses/UNIVERSAL_PALINDROME_CONDITION.md), section "The exact identity and its reach"\
**Proposed status:** superseded\
**Justification:** The [canonical gate](../simulations/neural/neural_translation_gate.py) contains exact palindromes with nonreal spectra and constructed unstable instances. The premise of necessary productive imperfection is false. A new atomic or crystalline mechanism would require its own generator and test.

---

### OQ-322

**Question:** Can specified quantum and neural coupling responses be compared through a physically justified dimensionless parameter?\
**Source:** [hypotheses/UNIVERSAL_PALINDROME_CONDITION.md](../hypotheses/UNIVERSAL_PALINDROME_CONDITION.md), section "Coupling and drive: finite censuses, open mechanism"\
**Proposed status:** open\
**Justification:** No universal window is established. The [neural report](../docs/neural/V_EFFECT_NEURAL.md) gives finite model/seed/resolution-dependent censuses; coupling and drive are different protocols. Define observables and units, refine the census and test a proposed mapping against controls before claiming a common optimum.

---

### OQ-324

**Question:** Does any independently defined observable at another scale have a physical threshold corresponding to quantum CΨ=1/4?\
**Source:** [hypotheses/UNIVERSAL_PALINDROME_CONDITION.md](../hypotheses/UNIVERSAL_PALINDROME_CONDITION.md), section "What the quarter does and does not transfer"\
**Proposed status:** open\
**Justification:** Normalized E/I fractions multiply to 1/4 at equality by arithmetic. A logistic sigmoid's maximum slope is a/4 for steepness a. Neither identity locates a neural transition or provides an atomic threshold. Supply a physical observable, units and an independent transition gate.

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
OQ-297, OQ-318, OQ-319, OQ-322, OQ-324 (5 entries from two files): OQ-319's silence/stability premise is superseded. The remaining entries require specified substrates, observables and independent gates; they are not evidence for cross-domain universality.

**Cluster 4: THE_PATTERN_RECOGNIZES_ITSELF**  
OQ-304, OQ-305, OQ-306, OQ-309 (4 entries): OQ-304/305/306 retain response and candidate-substrate gates. OQ-309 resolves the original biological comparison through the current support and instrument audit.

---

## Cross-batch overlaps

- **OQ-220 / OQ-155 / OQ-156:** Minimum per-pair entanglement question. OQ-220 (this batch) duplicates OQ-155/OQ-156 (batch 8b). Both resolved.
- **OQ-239-241 / OQ-238 / OQ-242:** BRIDGE_PROTOCOL FTL cluster. OQ-238 (numerical-verification, FALLEN) and OQ-242 (numerical-verification, FALSIFIED) overlap with OQ-239-241 (this batch, superseded).
- **OQ-297 / OQ-318-324:** THE_OTHER_SIDE parity propagation (OQ-297) connects to UNIVERSAL_PALINDROME_CONDITION open questions (OQ-318-324). Not duplicates but thematically linked.

---

## Patterns

1. **MEDIATOR_AS_QUANTUM_TRANSISTOR is this batch's worst artifact source.** 9 of 41 entries come from Section 8 "Limitations and Failure Modes," with 8 being documented design constraints (needs-human) and only 1 genuine open question (OQ-272).

2. **UNIVERSAL_PALINDROME_CONDITION mixes gates and rejected premises.** Candidate-substrate and comparison gates remain open; exact pairing does not require silence or stability (OQ-319).

3. **Superseded premises:** OQ-239-241 belong to the BRIDGE_PROTOCOL FTL refutation; OQ-319's exact-pairing silence/stability premise fails the constructed neural controls.

4. **THE_PATTERN_RECOGNIZES_ITSELF separates algebra from biology.** Three candidate/response gates remain open; the original control comparison is resolved without a biological pairing confirmation.

5. **Overall resolution rate: 7% (3/41).** Resolved status includes a negative disposition; it does not imply confirmation of the premise.
