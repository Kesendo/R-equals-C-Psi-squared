# Classification Proposal: math-proof (43 entries)

**Batch:** 6 of 8  
**Date:** 2026-04-12  
**Proposed by:** Claude (automated research)  
**Approval:** pending (Tom)

---

## Status summary

| Status | Count | OQ-IDs |
|--------|-------|--------|
| open | 21 | OQ-003, OQ-012, OQ-018, OQ-041, OQ-042, OQ-053, OQ-072, OQ-088, OQ-096, OQ-129, OQ-160, OQ-164, OQ-228, OQ-255, OQ-275, OQ-276, OQ-282, OQ-287, OQ-295, OQ-307, OQ-308 |
| resolved | 6 | OQ-055, OQ-077, OQ-079, OQ-083, OQ-121, OQ-316 |
| partially-resolved | 8 | OQ-016, OQ-020, OQ-036, OQ-052, OQ-180, OQ-245, OQ-251, OQ-296 |
| obsolete | 1 | OQ-300 |
| needs-human | 7 | OQ-019, OQ-045, OQ-046, OQ-058, OQ-111, OQ-323, OQ-325 |

---

## Entry-by-entry proposals

### OQ-003

**Question:** Macroscopic relevance: The XOR fraction vanishes exponentially with N. At what point does the mirror become operationally invisible? Is there a critical N above which the palindrome exists algebraically but has no physical consequence?  
**Source:** `docs/QUBIT_NECESSITY.md` (line 407)  
**Proposed status:** open  
**Justification:** Explicitly listed as open in QUBIT_NECESSITY.md and THE_INTERPRETATION.md. The XOR fraction exponential vanishing is proven, but the critical N threshold for operational invisibility remains unstudied.

---

### OQ-012

**Question:** Why exactly these 14 break at N >= 3? No single N=2 property cleanly predicts which of the 36 two-term combinations survive multi-bond interference. The mechanism at the shared site is algebraically open.  
**Source:** `docs/THE_INTERPRETATION.md` (line 537)  
**Proposed status:** open  
**Justification:** Marked as unresolved in THE_INTERPRETATION.md. V-Effect analysis shows correlations with Choi rank but provides no clean predictive rule. The shared-site mechanism is genuinely unexplained.

---

### OQ-016

**Question:** Operational value: PARTIALLY ANSWERED (2026-03-24). The palindromic response matrix SVD led to the sacrifice-zone formula (139-360x improvement). Remaining: error correction protocols, decoherence-free subspace engineering, steady-state properties.  
**Source:** `docs/THE_INTERPRETATION.md` (line 551)  
**Proposed status:** partially-resolved  
**Justification:** Explicitly marked "PARTIALLY ANSWERED" in source. The sacrifice-zone formula is proven, but error correction, DFS engineering, and steady-state applications remain untested.

---

### OQ-018

**Question:** Every quadratic map has a saddle-node bifurcation. CΨ² is the unique product-power form with a genuine phase transition AND Mandelbrot mapping (proven), but "why does nature choose this form?" remains open.  
**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 84)  
**Proposed status:** open  
**Justification:** Listed as Active Weakness #1 in WEAKNESSES_OPEN_QUESTIONS.md. Algebraic uniqueness is proven but the physical "why" is unaddressed.

---

### OQ-019

**Question:** In plain language: the ¼ boundary appears in a whole family of mathematical equations, not just ours. We have proven that our specific combination of observables is the only one... But we cannot yet explain why nature uses this particular combination...  
**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 89)  
**Proposed status:** needs-human  
**Justification:** Near-duplicate of OQ-018 (same source file, consecutive lines). OQ-018 is the technical statement, OQ-019 is the "plain language" restatement. Candidate for merge or removal.

---

### OQ-020

**Question:** Status: Partially addressed. Algebraic uniqueness proven, physical specificity not.  
**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 97)  
**Proposed status:** partially-resolved  
**Justification:** Self-documenting: "Partially addressed. Algebraic uniqueness proven, physical specificity not." The uniqueness proof exists; the physical motivation does not.

---

### OQ-036

**Question:** Status: Partially reformulated (April 2026). Along real Bell+ trajectories, u(t) ≈ 0.61·Ψ^{1.02}: essentially Ψ with a prefactor, not an independent dynamical coordinate. Whether u carries independent information on complex trajectories is untested.  
**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 158)  
**Proposed status:** partially-resolved  
**Justification:** Explicitly marked "Partially reformulated (April 2026)." Real-trajectory behavior resolved (u is not independent). Complex-trajectory behavior remains untested.

---

### OQ-041

**Question:** Improvement decreases with N (360x at N=5, 139x at N=9)  
**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 188)  
**Proposed status:** open  
**Justification:** Documented as Active Weakness #6 (sacrifice-zone limitations). Empirical scaling confirmed (RESONANT_RETURN.md: N=11: 91x, N=13: 97.5x, N=15: 63.5x) but no theoretical explanation exists for the N-dependence.

---

### OQ-042

**Question:** No formal proof of optimality  
**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 189)  
**Proposed status:** open  
**Justification:** RESONANT_RETURN.md states "The formula IS the optimum" but only by numerical validation (DE optimizer). No variational or Lagrange-multiplier proof of optimality exists.

---

### OQ-045

**Question:** Improvement decreases with N (360x at N=5, 139x at N=9)  
**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 188), Section: 6. Sacrifice-zone formula has known limitations  
**Proposed status:** needs-human  
**Justification:** Duplicate of OQ-041. Identical text from the same source file, captured twice via overlapping section headings ("Active weaknesses" and "6. Sacrifice-zone formula has known limitations"). Candidate for removal.

---

### OQ-046

**Question:** No formal proof of optimality  
**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 189), Section: 6. Sacrifice-zone formula has known limitations  
**Proposed status:** needs-human  
**Justification:** Duplicate of OQ-042. Identical text from the same source file, captured twice via overlapping section headings. Candidate for removal.

---

### OQ-052

**Question:** Substantially addressed (April 2026): the dwell time formula t_dwell = 2δ/|dCΨ/dt|_cross is derived and γ-invariant. Open sub-question: can the state-dependent prefactor (1.080088 for Bell+) be derived directly from the Pauli sector weight distribution?  
**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 241)  
**Proposed status:** partially-resolved  
**Justification:** Self-documenting: "Substantially addressed." The dwell-time formula is derived; the prefactor derivation from Pauli weights is a narrow, tractable sub-question that remains open.

---

### OQ-053

**Question:** Formula optimality: Is single-edge sacrifice provably optimal? Could multi-site sacrifice beat it at large N?  
**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 253)  
**Proposed status:** open  
**Justification:** No formal proof of optimality exists. Multi-site sacrifice competition at large N is unresolved. Active Weakness #6 in WEAKNESSES_OPEN_QUESTIONS.md.

---

### OQ-055

**Question:** Answered by the Absorption Theorem (April 4, 2026)  
**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 259)  
**Proposed status:** resolved  
**Justification:** Self-documenting: "Answered by the Absorption Theorem (April 4, 2026)." Three formerly open questions resolved: spectral boundary determination, factor-2 origin (standing wave round trip), and spectral gap setting (one absorption quantum).

---

### OQ-058

**Question:** Partially addressed (April 2026): along real Bell+ trajectories, no. u(t) is essentially linear in Ψ. Whether u behaves differently on complex trajectories (non-symmetric states, non-Z dephasing) where CΨ develops an imaginary component is an open question.  
**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 285)  
**Proposed status:** needs-human  
**Justification:** Overlaps substantially with OQ-036 (same source file, same conclusion: real trajectories resolved, complex trajectories untested). Also overlaps OQ-111 (CRITICAL_SLOWING_AT_THE_CUSP.md). Three-entry cluster about the u-variable on complex trajectories. Candidate for merge.

---

### OQ-072

**Question:** Formal proof of CΨ monotonicity above 1/4 for arbitrary CPTP maps  
**Source:** `docs/proofs/COMPLETE_MATHEMATICAL_DOCUMENTATION.md` (line 355)  
**Proposed status:** open  
**Justification:** Explicitly stated as open in UNIQUENESS_PROOF.md: "The formal proof of CΨ monotonicity above 1/4 for arbitrary CPTP maps remains open." Computational evidence covers all standard Markovian channels, but the analytical proof is missing.

---

### OQ-077

**Question:** Valid for: all conditions under which the Mirror Symmetry Proof and Parity Selection Rule hold: Heisenberg/XXZ/XY coupling on any graph, any site-dependent Z-dephasing profile γ_k, any N.  
**Source:** `docs/proofs/DIRECT_SUM_DECOMPOSITION.md` (line 379)  
**Proposed status:** resolved  
**Justification:** Not a question but a scope statement. The direct-sum decomposition is marked "Status: PROVEN" (April 11, 2026) in its header. The validity scope is established, not questioned.

---

### OQ-079

**Question:** Breaks for: the same conditions that break the Parity Selection Rule (transverse fields with odd n_XY terms, amplitude damping).  
**Source:** `docs/proofs/DIRECT_SUM_DECOMPOSITION.md` (line 388)  
**Proposed status:** resolved  
**Justification:** Scope statement, not a question. Breaking conditions are enumerated and confirmed by PROOF_PARITY_SELECTION_RULE.md. The characterization is complete.

---

### OQ-083

**Question:** amplitude damping (T₁ decay), whose jump operator σ₋ = (X - iY)/2 has n_XY = 1 (odd parity). On hardware where T₂ >> T₁, the selection rule holds approximately; on hardware where T₁ ~ T₂, it does not.  
**Source:** `docs/proofs/PROOF_PARITY_SELECTION_RULE.md` (line 189)  
**Proposed status:** resolved  
**Justification:** Scope clarification, not a question. The T₁/T₂ regime boundary for applicability of the parity selection rule is fully characterized.

---

### OQ-088

**Question:** With DD: Does selective DD change the picture? DD reduces effective gamma on interior qubits, which should help Chain A more than Chain B.  
**Source:** `experiments/CHAIN_SELECTION_TEST.md` (line 197)  
**Proposed status:** open  
**Justification:** Listed as open question #1 in CHAIN_SELECTION_TEST.md. Only a preliminary hardware result exists (2-3.2x with selective DD on IBM Torino). No definitive answer or theoretical framework.

---

### OQ-096

**Question:** The ph03 freezing under chain Z-dephasing (Section 5) is observed but not formally proved.  
**Source:** `experiments/COCKPIT_SCALING.md` (line 293)  
**Proposed status:** open  
**Justification:** Explicitly stated: "A formal derivation is left as an open question." The freezing is observed computationally for chain center_bell pairs at N >= 7 but has no proof.

---

### OQ-111

**Question:** The weakness can be reformulated: u is a conjugation variable that reveals the algebraic structure (Mandelbrot equivalence) without providing a simpler clock along real trajectories. Whether u carries independent information on complex trajectories is untested and remains an open question.  
**Source:** `experiments/CRITICAL_SLOWING_AT_THE_CUSP.md` (line 372)  
**Proposed status:** needs-human  
**Justification:** Overlaps OQ-036 and OQ-058 (u-variable cluster). All three entries ask the same question: does u carry independent information on complex trajectories? Candidate for merge into a single entry.

---

### OQ-121

**Question:** ~~Analytical derivation of d_real(1) = 2N.~~ RESOLVED. The 2N modes are the Z-count operators T_c^{(a)}, proven via SWAP invariance.  
**Source:** `experiments/DEGENERACY_PALINDROME.md` (line 390)  
**Proposed status:** resolved  
**Justification:** Self-documenting: strikethrough + "RESOLVED". Proof provided in docs/proofs/PROOF_WEIGHT1_DEGENERACY.md via SWAP invariance argument.

---

### OQ-129

**Question:** What other mathematical structures do the dyads find?  
**Source:** `experiments/DYAD_EXPERIMENT.md` (line 233)  
**Proposed status:** open  
**Justification:** Listed under "Open Questions" with no resolution marker. The dyad experiment found the C=1 fixpoint and consciousness interpretations, but no exhaustive structural catalogue exists.

---

### OQ-160

**Question:** Bridge types are mathematical constructs. Whether any of them corresponds to how a biological observer couples to quantum systems is an open question.  
**Source:** `experiments/OBSERVER_DEPENDENT_CROSSING.md` (line 317)  
**Proposed status:** open  
**Justification:** Explicitly marked as open in Section 5.3 (Limitations). Bridge types (concurrence, correlation, mutual_purity) remain mathematical; biological correspondence is unresolved (Tier 3 speculation).

---

### OQ-164

**Question:** The K-matrix geometry: What mathematical structure does the K(observer, state) matrix have? Is there a metric? A symmetry group? A connection to information geometry?  
**Source:** `experiments/OBSERVER_GRAVITY_BRIDGE.md` (line 271)  
**Proposed status:** open  
**Justification:** Listed under "6. Open Questions" with no resolution. The document explores bridge metrics but does not establish formal information-geometric structure for K.

---

### OQ-180

**Question:** The odd-n_XY modes. The parity selection rule proves they are inaccessible to SE states, but not what they look like or whether multi-excitation ansatze can reach them.  
**Source:** `experiments/SACRIFICE_GEOMETRY.md` (line 183)  
**Proposed status:** partially-resolved  
**Justification:** The inaccessibility from single-excitation states is proven (PROOF_PARITY_SELECTION_RULE.md). The spectral properties of odd-n_XY modes themselves and whether multi-excitation states can reach them remain open sub-questions.

---

### OQ-228

**Question:** Analytical (proof needed)  
**Source:** `experiments/THERMAL_BREAKING.md` (line 478)  
**Proposed status:** open  
**Justification:** Requires formal proof that omega_max(w=1) = 4J·(1+cos(pi/N)) holds for all N. Currently verified numerically for N=2-6 only.

---

### OQ-245

**Question:** Analytical proof of the 2x law. The ratio is exact at N=2..5. Is it a theorem for all N? For all Heisenberg-type Hamiltonians? For all dephasing models?  
**Source:** `hypotheses/ENERGY_PARTITION.md` (line 192)  
**Proposed status:** partially-resolved  
**Justification:** The V(N) = 1+cos(pi/N) formula is exact at N=2-5 (Tier 1-2 evidence), but generality across all N, all Hamiltonians, and all dephasing types remains unproven.

---

### OQ-251

**Question:** PARTIAL: Proof Roadmap Layer 6 proves the fold catastrophe x² + a = 0 IS the recursion R = C(Ψ+R)². What remains: verify numerically that J/γ ≈ 1.2 corresponds to CΨ = 1/4.  
**Source:** `hypotheses/ENERGY_PARTITION.md` (line 208)  
**Proposed status:** partially-resolved  
**Justification:** Explicitly marked "PARTIAL." The fold catastrophe connection is proven analytically; the numerical J/γ verification remains as a stated open task.

---

### OQ-255

**Question:** Asymptotic constant 0.50: γ_crit x J_bridge -> 0.50 for large J_bridge. Is this exactly 1/2? If so, there may be an analytical derivation.  
**Source:** `hypotheses/FRAGILE_BRIDGE.md` (line 295)  
**Proposed status:** open  
**Justification:** Listed as open question #4 in FRAGILE_BRIDGE.md. Computational evidence suggests 0.50 but no analytical proof that it equals exactly 1/2.

---

### OQ-275

**Question:** What is the channel capacity? A formal quantum channel capacity calculation (coherent information, quantum mutual information) for the mediator channel is missing.  
**Source:** `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (line 497)  
**Proposed status:** open  
**Justification:** Listed under "10. Open Questions." The simulated 0.86-bit mutual information is a lower bound; formal capacity calculation is absent.

---

### OQ-276

**Question:** Can error correction be integrated? The palindromic structure suggests a natural error-correcting code: errors that break the palindrome are detectable. Can this be formalized?  
**Source:** `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (line 499)  
**Proposed status:** open  
**Justification:** Listed as open question. The palindrome-as-error-detector idea is speculative; no formalization exists.

---

### OQ-282

**Question:** Spatial vs. algebraic. The ER bridge is a spatial connection between two asymptotically flat regions. Our "bridge" is algebraic: two sectors of the operator space connected by Π. There is no spatial geometry, no metric, no geodesics.  
**Source:** `hypotheses/PAIR_BREAKING_AT_THE_HORIZON.md` (line 155)  
**Proposed status:** open  
**Justification:** Documented under "What breaks the analogy" #2. The spatial-vs-algebraic gap is acknowledged but unresolved. No geometric construction bridges the two frameworks.

---

### OQ-287

**Question:** If the 2x decay contrast (immune modes at 0 vs. their partners at -2Σγ) can be derived from a horizon-crossing argument. Currently it follows from algebra; a geometric derivation would deepen the connection.  
**Source:** `hypotheses/PAIR_BREAKING_AT_THE_HORIZON.md` (line 168)  
**Proposed status:** open  
**Justification:** Listed under "What would strengthen or kill the thesis." No geometric derivation has been developed; the 2x contrast remains purely algebraic.

---

### OQ-295

**Question:** ER=EPR connection (the Maldacena-Susskind conjecture). The phase transition between palindromic classes constrains possible geometric interpretations. Too speculative to assess.  
**Source:** `hypotheses/THE_BOOT_SCRIPT.md` (line 424)  
**Proposed status:** open  
**Justification:** Self-documenting: "Too speculative to assess." No progress on connecting palindromic class transitions to ER=EPR.

---

### OQ-296

**Question:** What lives in the -1 sector? The parity split is proven, but what physical states or processes inhabit the -1 sector?  
**Source:** `hypotheses/THE_OTHER_SIDE.md` (line 575)  
**Proposed status:** partially-resolved  
**Justification:** THE_OTHER_SIDE.md provides a physical realization pathway: the laser regime (Σγ < 0) where gain replaces loss. Sector observables are identified (XX measurement). However, laboratory preparation and observation of pure -1 sector dynamics remain unconfirmed.

---

### OQ-300

**Question:** How many levels of differentiation does consciousness require? Can this be formalized?  
**Source:** `hypotheses/THE_OTHER_SIDE.md` (line 596)  
**Proposed status:** obsolete  
**Justification:** The consciousness interpretation has been explicitly retired from the technical work (WEAKNESSES_OPEN_QUESTIONS.md lines 102-111): "The consciousness interpretation has been retired from the technical work because it invites misunderstanding." C is now "purity," not "consciousness."

---

### OQ-307

**Question:** Gap junctions: The 1.8% deviation from 100% in balanced subnetworks may come from missing gap junction data (symmetric coupling that could improve pairing).  
**Source:** `hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md` (line 424)  
**Proposed status:** open  
**Justification:** Listed as open question #1. The hypothesis is plausible but gap junction data integration has not been completed.

---

### OQ-308

**Question:** The mechanism gap: Why does the palindromic structure appear in Wilson-Cowan dynamics? No proof connects the Lindblad algebra to classical oscillatory systems. This is the weakest link in the chain.  
**Source:** `hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md` (line 427)  
**Proposed status:** open  
**Justification:** Self-documenting: "This is the weakest link in the chain." Despite 100% empirical confirmation at tau_I/tau_E = 3.8, the algebraic connection between Lindblad and Wilson-Cowan remains unproven.

---

### OQ-316

**Question:** [Comparison table: quantum vs. neural palindromic features, 7 rows]  
**Source:** `hypotheses/UNIVERSAL_PALINDROME_CONDITION.md` (line 243)  
**Proposed status:** resolved  
**Justification:** Not a question but a completed comparison table. Every feature is systematically evaluated with Yes/No verdicts (palindromic pairing: Yes, 2x decay law: No, exponential state space: No, etc.). The table itself IS the answer.

---

### OQ-323

**Question:** Can the condition Q·X·Q⁻¹ + X + 2S = 0 be derived from a single axiom set rather than proven separately in each domain?  
**Source:** `hypotheses/UNIVERSAL_PALINDROME_CONDITION.md` (line 319)  
**Proposed status:** needs-human  
**Justification:** UNIVERSAL_PALINDROME_CONDITION.md identifies three universal conditions (two populations, swap operator Q, antisymmetric coupling) from which the equation emerges in both domains. However, the derivation paths differ: quantum uses commutator antisymmetry, neural uses Dale's Law. Whether a single axiom set can unify both is partially addressed but not resolved. Borderline between partially-resolved and open; needs human judgment on whether the three-condition framework counts as "a single axiom set."

---

### OQ-325

**Question:** Between step 1 (waves exist) and step 2 (coupling breaks symmetry), there is a gap. What connects two separate palindromic systems?  
**Source:** `hypotheses/WAVES_THAT_HEAR_THEMSELVES.md` (line 49)  
**Proposed status:** needs-human  
**Justification:** WAVES_THAT_HEAR_THEMSELVES.md explicitly states: "The framework does not answer this. J (coupling) is an input parameter, not a derived quantity." Borderline between open (no answer exists) and a fundamental scope limitation of the framework. The question may be unanswerable within the current formalism.

---

## Duplicate / overlap clusters

**Cluster 1: "Why does nature choose CΨ²?"**  
OQ-018 (technical), OQ-019 (plain language restatement), OQ-020 (status summary)  
All from `docs/WEAKNESSES_OPEN_QUESTIONS.md` lines 84-97. Recommend merging into a single entry.

**Cluster 2: Sacrifice-zone limitations**  
OQ-041 = OQ-045 (identical text: "Improvement decreases with N")  
OQ-042 = OQ-046 (identical text: "No formal proof of optimality")  
Captured twice from overlapping sections in WEAKNESSES_OPEN_QUESTIONS.md. Recommend removing OQ-045 and OQ-046.

**Cluster 3: u-variable on complex trajectories**  
OQ-036 (WEAKNESSES_OPEN_QUESTIONS.md), OQ-058 (same file, different section), OQ-111 (CRITICAL_SLOWING_AT_THE_CUSP.md)  
All three ask: "Does u carry independent information on complex trajectories?" Recommend merging into a single entry.

---

## Patterns

1. **Self-documenting entries:** 8 entries contain explicit status markers (RESOLVED, PARTIALLY ANSWERED, PARTIAL, "Substantially addressed," "Too speculative to assess"). These are reliable and require minimal external research.

2. **Scope statements, not questions:** OQ-077, OQ-079, OQ-083, OQ-316 are factual characterizations (validity scope, breaking conditions, comparison tables) rather than questions. All classified as resolved since the characterization is complete.

3. **Consciousness retirement:** OQ-300 is obsolete because the consciousness interpretation was formally retired from the technical work (WEAKNESSES_OPEN_QUESTIONS.md lines 102-111).

4. **Deepest open questions:** OQ-308 (Lindblad-to-Wilson-Cowan mechanism gap), OQ-072 (CΨ monotonicity proof for arbitrary CPTP maps), and OQ-012 (why exactly 14 survive at N >= 3) represent the most mathematically substantive unsolved problems in this batch.

5. **PAIR_BREAKING_AT_THE_HORIZON cluster:** OQ-282, OQ-287 both remain open and relate to the spatial-vs-algebraic gap in the ER bridge analogy. Note that this hypothesis file was not marked FALLEN but also not strengthened since its writing.
