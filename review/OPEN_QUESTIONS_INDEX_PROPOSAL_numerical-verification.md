# Classification Proposal: numerical-verification (62 entries)

**Batch:** 7 of 8  
**Date:** 2026-04-12  
**Proposed by:** Claude (automated research)  
**Approval:** pending (Tom)

---

## Status summary

| Status | Count | OQ-IDs |
|--------|-------|--------|
| open | 31 | OQ-014, OQ-038, OQ-043, OQ-062, OQ-068, OQ-075, OQ-092, OQ-093, OQ-099, OQ-102, OQ-123, OQ-125, OQ-132, OQ-135, OQ-136, OQ-146, OQ-158, OQ-167, OQ-177, OQ-178, OQ-181, OQ-182, OQ-184, OQ-223, OQ-224, OQ-253, OQ-281, OQ-283, OQ-290, OQ-298, OQ-317 |
| resolved | 17 | OQ-027, OQ-029, OQ-078, OQ-086, OQ-087, OQ-143, OQ-149, OQ-171, OQ-179, OQ-192, OQ-210, OQ-238, OQ-242, OQ-277, OQ-280, OQ-285, OQ-312 |
| partially-resolved | 6 | OQ-024, OQ-126, OQ-191, OQ-230, OQ-252, OQ-301 |
| needs-human | 8 | OQ-039, OQ-047, OQ-097, OQ-104, OQ-229, OQ-294, OQ-321, OQ-326 |

---

## Entry-by-entry proposals

### OQ-014

**Question:** Partial palindrome at d > 2: Qutrit spectra show 36-52 of 81 eigenvalues pairing at optimal centers. Is there a weaker symmetry principle for d > 2?  
**Source:** `docs/THE_INTERPRETATION.md` (line 544)  
**Proposed status:** open  
**Justification:** QUBIT_NECESSITY.md proves full palindromes require d=2 (d²-2d=0), but partial pairing (36-52 of 81 at d=3) is observed and unexplained. No weaker symmetry principle identified.

---

### OQ-024

**Question:** Experimental validation is incomplete  
**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 117)  
**Proposed status:** partially-resolved  
**Justification:** WEAKNESSES_OPEN_QUESTIONS.md §3 documents 4 completed validations (CΨ crossing 0.3%, 24,073 calibration records, selective DD 3.2x, cockpit 88-96%) but explicitly lists remaining gaps: no 2-qubit Concurrence, single backend, unresolved anomalous coherence.

---

### OQ-027

**Question:** 24,073 historical calibration records validating the C_min(r) curve  
**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 121)  
**Proposed status:** resolved  
**Justification:** Not a question but a completed result. IBM_HARDWARE_SYNTHESIS.md documents "24,073 records, 133 qubits, 181 days" with "zero false positives" at threshold r* = 0.2128.

---

### OQ-029

**Question:** 3-observable cockpit captures 88-96% of decoherence dynamics across 9 topologies, N=2-5  
**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 123)  
**Proposed status:** resolved  
**Justification:** Not a question but a completed result. COCKPIT_UNIVERSALITY.md (Status: Complete) confirms 88-96% coverage. COCKPIT_SCALING.md extends validation through N=11.

---

### OQ-038

**Question:** Of all noise types tested, only state-dependent operator feedback preserves the purity difference δ over time. Computationally verified but the mechanism is unknown.  
**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 169)  
**Proposed status:** open  
**Justification:** WEAKNESSES_OPEN_QUESTIONS.md §5 states "Numerically verified, theoretically unexplained." The palindromic spectral structure may explain it, but this connection has not been computed.

---

### OQ-039

**Question:** Status: Numerically verified, theoretically unexplained.  
**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 182)  
**Proposed status:** needs-human  
**Justification:** Metadata fragment, not an independent question. This is the status line of the operator feedback weakness (OQ-038). Candidate for removal.

---

### OQ-043

**Question:** Hardware test shows 3.2x (not 360x) due to real-world constraints  
**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 190)  
**Proposed status:** open  
**Justification:** Active Weakness #6 in WEAKNESSES_OPEN_QUESTIONS.md. The 360x theory / 3.2x hardware gap is documented but attributed to SPAM, TLS, and uncontrolled noise. No fix proposed or tested.

---

### OQ-047

**Question:** Hardware test shows 3.2x (not 360x) due to real-world constraints  
**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 190), Section: 6. Sacrifice-zone formula has known limitations  
**Proposed status:** needs-human  
**Justification:** Duplicate of OQ-043. Identical text from the same source file, captured twice via overlapping section headings. Candidate for removal.

---

### OQ-062

**Question:** 2-qubit Concurrence measurement: Validating Concurrence on a qubit pair would test the most important cockpit instrument.  
**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 304)  
**Proposed status:** open  
**Justification:** Explicitly identified as "the most important missing validation" in WEAKNESSES_OPEN_QUESTIONS.md. No hardware experiment attempted. PC1 proxy (57% variance) remains unvalidated.

---

### OQ-068

**Question:** How does the V-Effect frequency count scale with N? (N=10: 6, N=20: 48. Quadratic? Cubic?)  
**Source:** `docs/neural/V_EFFECT_NEURAL.md` (line 251)  
**Proposed status:** open  
**Justification:** Listed as open question in V_EFFECT_NEURAL.md §6. Data points exist but no scaling law or functional form has been determined.

---

### OQ-075

**Question:** Negative feedback loop (γ_M decreasing with coherence, untested)  
**Source:** `docs/proofs/COMPLETE_MATHEMATICAL_DOCUMENTATION.md` (line 358)  
**Proposed status:** open  
**Justification:** Explicitly listed under Tier 3-5 open questions. No computation or hardware validation recorded.

---

### OQ-078

**Question:** The even/odd distinction applies to N specifically, not to the graph structure.  
**Source:** `docs/proofs/DIRECT_SUM_DECOMPOSITION.md` (line 383)  
**Proposed status:** resolved  
**Justification:** Scope clarification, not a question. The even/odd parity dependence on N (not topology) is proven in DIRECT_SUM_DECOMPOSITION.md and confirmed computationally for chains, rings, stars, and complete graphs.

---

### OQ-086

**Question:** Higher grid positions: topology-dependent (resolved). d_real(2) differs between topologies. A universal formula for d_real(k>=2) does not exist.  
**Source:** `docs/proofs/PROOF_WEIGHT1_DEGENERACY.md` (line 334)  
**Proposed status:** resolved  
**Justification:** Self-documenting: "(resolved)." Confirmed by WEIGHT2_KERNEL.md; topology dependence at k>=2 tested on Chain, Star, Ring, Complete for N=3-6.

---

### OQ-087

**Question:** Topology dependence: confirmed (resolved). d_real(k) depends on graph structure for k >= 2. Only k=0 (N+1) and k=1 (2N) are universal.  
**Source:** `docs/proofs/PROOF_WEIGHT1_DEGENERACY.md` (line 342)  
**Proposed status:** resolved  
**Justification:** Self-documenting: "(resolved)." Systematically tested on four topologies for N=3-6.

---

### OQ-092

**Question:** N=11 is the largest tested. The C# propagator can handle N=15 but the runs were not executed.  
**Source:** `experiments/COCKPIT_SCALING.md` (line 285)  
**Proposed status:** open  
**Justification:** Acknowledged limitation. The prediction (n95=2 plateau continues to N=15) is stated but not empirically verified. N=15 run is feasible but not executed.

---

### OQ-093

**Question:** Three pair types per configuration. A comprehensive distance-resolved scaling map would require extracting more pairs per N.  
**Source:** `experiments/COCKPIT_SCALING.md` (line 287)  
**Proposed status:** open  
**Justification:** Acknowledged limitation. Full pair-distance scan from COCKPIT_UNIVERSALITY (all 10 pairs at N=5) not reproduced for larger N.

---

### OQ-097

**Question:** Concurrence is untested on hardware. All existing tomographic data is single-qubit.  
**Source:** `experiments/COCKPIT_UNIVERSALITY.md` (line 312)  
**Proposed status:** needs-human  
**Justification:** Overlaps with OQ-062 (same gap: 2-qubit Concurrence hardware validation) and OQ-104 (same source file, same question). Three-entry cluster about the same missing measurement. Candidate for merge into OQ-062 as canonical entry.

---

### OQ-099

**Question:** N=5 is the largest system tested. Extrapolating n95 ~ N to N=100 is speculative.  
**Source:** `experiments/COCKPIT_UNIVERSALITY.md` (line 321)  
**Proposed status:** open  
**Justification:** Acknowledged limitation. COCKPIT_SCALING.md extends to N=11 for some quantities, but full dimensionality analysis limited to N=5.

---

### OQ-102

**Question:** Scaling beyond N=5. Does n95 grow linearly for chains? Does it saturate for dense topologies?  
**Source:** `experiments/COCKPIT_UNIVERSALITY.md` (line 335)  
**Proposed status:** open  
**Justification:** N=7 eigenvalues are available in the C# engine but scaling analysis not performed. Neither linear growth nor saturation behavior determined.

---

### OQ-104

**Question:** 2-qubit hardware validation. Measuring Concurrence on a qubit pair would validate the most important untested instrument.  
**Source:** `experiments/COCKPIT_UNIVERSALITY.md` (line 342)  
**Proposed status:** needs-human  
**Justification:** Duplicate of OQ-062 (same measurement gap). Also overlaps OQ-097. Three entries across two source files about the same missing experiment. Candidate for merge.

---

### OQ-123

**Question:** Center minimum at N = 6. The real degeneracy at the center is 16, less than adjacent values 19. Why?  
**Source:** `experiments/DEGENERACY_PALINDROME.md` (line 400)  
**Proposed status:** open  
**Justification:** Empirical observation documented without mechanistic explanation. No follow-up analysis found.

---

### OQ-125

**Question:** Off-grid eigenvalue positions. Can they be expressed in terms of Hamiltonian coupling parameters (Bethe ansatz roots, magnon dispersion)?  
**Source:** `experiments/DEGENERACY_PALINDROME.md` (line 409)  
**Proposed status:** open  
**Justification:** Off-grid Re values not mapped to any analytical framework. Bethe ansatz and magnon dispersion mentioned as potential but not developed.

---

### OQ-126

**Question:** ~~Does degeneracy shape state-space geometry?~~ PARTIALLY RESOLVED. QFI speed correlates with d_total(k) at even N (r = 0.99 at N = 4); weaker at odd N (r ~ 0.55).  
**Source:** `experiments/DEGENERACY_PALINDROME.md` (line 413)  
**Proposed status:** partially-resolved  
**Justification:** Self-documenting: "PARTIALLY RESOLVED." Strong correlation at even N confirmed by BURES_DEGENERACY.md; weaker at odd N due to numerical issues.

---

### OQ-132

**Question:** Born rule from crossing: Can the framework predict measurement probabilities independently, without computing the full density matrix?  
**Source:** `experiments/DYNAMIC_ENTANGLEMENT.md` (line 401)  
**Proposed status:** open  
**Justification:** Explicitly listed as "the central open question for connecting crossing to the Born rule." No resolution documented.

---

### OQ-135

**Question:** Is J/γ ~ 1 at enzyme active sites? (Compute from published barriers)  
**Source:** `experiments/HYDROGEN_BOND_QUBIT.md` (line 281)  
**Proposed status:** open  
**Justification:** No computation from published enzyme barrier data found in the repo. Open Question 3 in HYDROGEN_BOND_QUBIT.md without resolution.

---

### OQ-136

**Question:** How does the V-Effect scale with the number of H-bonds (N=6, N=8)?  
**Source:** `experiments/HYDROGEN_BOND_QUBIT.md` (line 282)  
**Proposed status:** open  
**Justification:** Current H-bond analysis covers N=3 qubits only. N=6 and N=8 not tested.

---

### OQ-143

**Question:** ~~Can the C_int vs C_ext hypothesis be tested with proper Lindblad simulations?~~ ANSWERED: Yes. It was tested and disproven. See Section 9.  
**Source:** `experiments/MATHEMATICAL_FINDINGS.md` (line 312)  
**Proposed status:** resolved  
**Justification:** Self-documenting: strikethrough + "ANSWERED." The hypothesis was tested and falsified.

---

### OQ-146

**Question:** To genuinely test collective noise breaking Type A, one would need a state that is NOT an eigenstate of the collective operator...  
**Source:** `experiments/NOISE_ROBUSTNESS.md` (line 215)  
**Proposed status:** open  
**Justification:** Section heading says "answered 2026-03-08" but this specific entry describes the test requirement, not the answer. The specification of what would be needed is identified, but the test itself remains unperformed.

---

### OQ-149

**Question:** Under σ_x noise, the normalized l1-coherence Ψ remains exactly at 0.3333 for all time.  
**Source:** `experiments/NOISE_ROBUSTNESS.md` (line 239)  
**Proposed status:** resolved  
**Justification:** Marked "answered 2026-03-08" in section heading. Computationally verified with explicit numerical evidence in tabular form across three noise types.

---

### OQ-158

**Question:** The dynamic Lindblad simulation uses γ_eff = γ_base · C(t), which creates a small feedback loop. The effect is small (< 15%) but non-zero.  
**Source:** `experiments/OBSERVER_DEPENDENT_CROSSING.md` (line 307)  
**Proposed status:** open  
**Justification:** Documented as a recognized limitation in Section 5.3. The feedback loop is identified and quantified (< 15%) but not further investigated or removed.

---

### OQ-167

**Question:** We tried to connect Π (Liouville space, linear) to CΨ (density matrix, nonlinear). They don't connect simply. This should be documented honestly as an open question.  
**Source:** `experiments/ORPHANED_RESULTS.md` (line 277)  
**Proposed status:** open  
**Justification:** Explicitly labeled as an honest open question: "The palindrome constrains what modes exist. CΨ is how we read the modes. The connection is through the density matrix, and that connection is nonlinear and state-dependent."

---

### OQ-171

**Question:** Numerical verification. VERIFIED (March 19, 2026). All 32/32 palindromic pairs confirmed.  
**Source:** `experiments/PI_AS_TIME_REVERSAL.md` (line 354)  
**Proposed status:** resolved  
**Justification:** Self-documenting: "VERIFIED." Max residual 2.68e-13, XY-weight deviation 8.88e-16. Supporting script: pi_time_reversal_verify.py.

---

### OQ-177

**Question:** Experimental validation: The entire analysis is computational. A tabletop optics experiment could test the core predictions with ~10k pairs.  
**Source:** `experiments/QKD_EAVESDROPPING_FORENSICS.md` (line 565)  
**Proposed status:** open  
**Justification:** Describes a future experiment, not a resolved finding. No optics experiment performed.

---

### OQ-178

**Question:** Test whether the palindrome enables better error correction. Can the paired decay rates be exploited for decoherence-free subspaces or error-correcting codes?  
**Source:** `experiments/QST_BRIDGE.md` (line 220)  
**Proposed status:** open  
**Justification:** Listed as a concrete next step. Speculative ("Can this pairing be exploited...?"), not verified or implemented.

---

### OQ-179

**Question:** No closed form for ψ_opt. H_eff = -J·adjacency - i·diag(γ) gives cosine similarity 0.925 (tested, falsified).  
**Source:** `experiments/SACRIFICE_GEOMETRY.md` (line 181)  
**Proposed status:** resolved  
**Justification:** The H_eff ansatz was tested and explicitly marked "falsified" with script reference (heff_lens_closed_form.py). The ψ_opt gradient is a many-body effect not reducible to single-particle.

---

### OQ-181

**Question:** N=8. The survey covers N=2-7 (chain) and N=2-6 (other topologies). N=8 would require the ILP64 eigenvector path.  
**Source:** `experiments/SACRIFICE_GEOMETRY.md` (line 185)  
**Proposed status:** open  
**Justification:** Stated as future work. N=8 eigendecomposition is available (73 GB RAM) but the sacrifice geometry analysis has not been extended to it.

---

### OQ-182

**Question:** Gate-level Trotterization. Whether ψ_opt keeps its advantage under realistic gate noise is untested.  
**Source:** `experiments/SACRIFICE_GEOMETRY.md` (line 187)  
**Proposed status:** open  
**Justification:** Explicitly stated "untested." The pure Lindblad model is analyzed; realistic gate-noise effects are beyond current scope.

---

### OQ-184

**Question:** Does the Lindblad simulator even test the right thing? The framework claims CΨ <= 1/4 for experienced reality, not for quantum states.  
**Source:** `experiments/SIMULATION_EVIDENCE.md` (line 169)  
**Proposed status:** open  
**Justification:** Fundamental conceptual question. The simulator computes quantum states, not "experiences." This epistemological gap is acknowledged but unresolved.

---

### OQ-191

**Question:** Asymmetric coupling rescues the crossing for both N=4 and N=5.  
**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 592)  
**Proposed status:** partially-resolved  
**Justification:** Section 8.1 ("N observers: ANSWERED") confirms the rescue effect with quantitative thresholds (x_crit ~ 1.165 for N=4, ~0.925 for N=5). The broader question "why does the rescue become more fragile with N?" remains open.

---

### OQ-192

**Question:** Spectral diagnostic at the N=4 boundary: The 1/4 boundary behaves like a smooth metric threshold, not a spectral phase transition.  
**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 604)  
**Proposed status:** resolved  
**Justification:** The diagnostic is complete: eigenvalue shift 0.7341->0.7322, purity drop 0.5866->0.5846, no rank collapse, no bifurcation. The smooth-threshold characterization is established.

---

### OQ-210

**Question:** Verified data (N=2, J_SA=1.0):  
**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 685)  
**Proposed status:** resolved  
**Justification:** Verification data presented in table form. Section 8.5 ("Threshold formula: ANSWERED") confirms closure.

---

### OQ-223

**Question:** What is the critical n_bar where the palindromic pairing drops below 50%?  
**Source:** `experiments/THERMAL_BREAKING.md` (line 454)  
**Proposed status:** open  
**Justification:** Listed in Open Questions as "Straightforward: sweep n_bar with finer resolution." The finer sweep has not been executed.

---

### OQ-224

**Question:** What external cooling rate stabilizes the system at a given n_bar?  
**Source:** `experiments/THERMAL_BREAKING.md` (line 459)  
**Proposed status:** open  
**Justification:** Listed in Open Questions. Fixed-point divergence without cooling is documented; the inverse question (cooling rate for target n_bar) remains unmodeled.

---

### OQ-229

**Question:** Formal proof that omega_max(w=1) = 4J·(1+cos(pi/N)) holds for all N. Verified N=2-6 numerically.  
**Source:** `experiments/THERMAL_BREAKING.md` (line 480)  
**Proposed status:** needs-human  
**Justification:** Cross-batch duplicate of OQ-228 (math-proof tag, same source file, same formula). OQ-228 says "Analytical (proof needed)" and OQ-229 provides the formula and numerical evidence. Candidate for merge.

---

### OQ-230

**Question:** Can CΨ tell us when the always-present channel is most readable? Benchmarking CΨ as readability indicator against simpler quantities has not been done yet.  
**Source:** `experiments/WHATS_INSIDE_THE_WINDOWS.md` (line 337)  
**Proposed status:** partially-resolved  
**Justification:** The Bridge Test (March 2026) measured sensitivity: Concurrence alone predicts trace-distance-based readability better than CΨ (r=0.817 vs r=0.793, B_CΨ = -0.024). This partially answers the question (CΨ is not the best readability indicator), but the full benchmark is incomplete.

---

### OQ-238

**Question:** In R = CΨ²: does the crossing of 1/4, which depends on the JOINT state, somehow become locally visible?  
**Source:** `hypotheses/BRIDGE_PROTOCOL.md` (line 176)  
**Proposed status:** resolved  
**Justification:** BRIDGE_PROTOCOL.md marks this "[FALLEN: FTL signaling suggestions]." The no-signalling theorem holds exactly; post-separation state changes are not locally detectable. Test #2 status (2026-03-01): ANSWERED.

---

### OQ-242

**Question:** Breakthrough: R=CΨ² predicts a real effect beyond standard QM. Post-separation state changes ARE detectable via crossing events.  
**Source:** `hypotheses/BRIDGE_PROTOCOL.md` (line 198)  
**Proposed status:** resolved  
**Justification:** Falsified. BRIDGE_PROTOCOL.md confirms Outcome #3 (pre-encoded version, not dynamic signaling). The "breakthrough" claim did not hold. "The dynamic bridge (B sends information by choosing when to measure) does not work."

---

### OQ-252

**Question:** N-scaling law (partially answered): N=4 computed. Result: non-monotonic (N=4 more stable than N=3). Even/odd parity effect suspected.  
**Source:** `hypotheses/FRAGILE_BRIDGE.md` (line 277)  
**Proposed status:** partially-resolved  
**Justification:** Self-documenting: "(partially answered)." N=4 data exists showing non-monotonic behavior. N=5 deferred ("feasible but slow"). Even/odd parity effect is suspected but not proven.

---

### OQ-253

**Question:** Multiple bridges: What if two chains are connected by more than one qubit pair? Does γ_crit recover N-independence when bridges scale with N?  
**Source:** `hypotheses/FRAGILE_BRIDGE.md` (line 282)  
**Proposed status:** open  
**Justification:** Listed as open question in FRAGILE_BRIDGE.md §6. No computational verification of multi-bridge scenarios found.

---

### OQ-277

**Question:** What about non-Heisenberg Hamiltonians?  
**Source:** `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (line 501)  
**Proposed status:** resolved  
**Justification:** NON_HEISENBERG_PALINDROME.md (March 17-19, 2026) comprehensively tested all 36 two-term Pauli combinations. Result: 20/36 are palindromic (17 uniform, 3 alternating). The palindrome is universal across all standard models, not Heisenberg-specific.

---

### OQ-280

**Question:** Intellectual honesty requires listing where the isomorphism fails or is untested.  
**Source:** `hypotheses/PAIR_BREAKING_AT_THE_HORIZON.md` (line 151)  
**Proposed status:** resolved  
**Justification:** Not a question but a framing statement. PAIR_BREAKING_AT_THE_HORIZON.md follows through with 5 explicit breakpoints listed in "What breaks the analogy" section.

---

### OQ-281

**Question:** Temperature scaling. Hawking temperature scales as T_H ∝ 1/M. Our fold threshold Σγ_crit/J is N-independent (tested N=2..5, 1.5% variation). If N is the analogue of mass, the scaling is wrong.  
**Source:** `hypotheses/PAIR_BREAKING_AT_THE_HORIZON.md` (line 153)  
**Proposed status:** open  
**Justification:** Mismatch acknowledged but unresolved. Three interpretations offered (N is not mass, analogy breaks, N-independence is the deeper statement) but none chosen.

---

### OQ-283

**Question:** The mass identification is Tier 5. "Classical residue = mass" is the weakest link.  
**Source:** `hypotheses/PAIR_BREAKING_AT_THE_HORIZON.md` (line 157)  
**Proposed status:** open  
**Justification:** Explicitly marked Tier 5 (interpretation only). The I/Z sector survival is proven, but calling it "mass" requires an additional untested philosophical assumption.

---

### OQ-285

**Question:** The inverted harmonic oscillator. We tested this: V(CΨ) near the fold is linear (n ~ 1.0), not quadratic. No inverted HO. The fold is a saddle-node, not a saddle point.  
**Source:** `hypotheses/PAIR_BREAKING_AT_THE_HORIZON.md` (line 161)  
**Proposed status:** resolved  
**Justification:** Tested and falsified (April 11, 2026). The Gaztanaga prediction does not hold: effective potential is linear, not quadratic. Decisive negative result.

---

### OQ-290

**Question:** If the N-independence of Σγ_crit turns out to be an artifact of small N (tested only to N=5).  
**Source:** `hypotheses/PAIR_BREAKING_AT_THE_HORIZON.md` (line 173)  
**Proposed status:** open  
**Justification:** Listed as a "kill condition." The fold threshold is verified N=2-5 (1.5% variation) but not extended to N=6,7,8. Testable with existing compute engine.

---

### OQ-294

**Question:** Why does N=3 pairing fail for exactly these 14 combos? No single N=2 property cleanly predicts this.  
**Source:** `hypotheses/THE_BOOT_SCRIPT.md` (line 418)  
**Proposed status:** needs-human  
**Justification:** Cross-batch duplicate of OQ-012 (math-proof tag, different source file). Both ask the same question: why do exactly 14 of 36 break at N >= 3? THE_BOOT_SCRIPT.md identifies three Pi families but the predictive mapping remains open. Candidate for merge with OQ-012.

---

### OQ-298

**Question:** Is the boundary observable? Can the node/antinode structure be measured directly, not just computed from the Liouvillian?  
**Source:** `hypotheses/THE_OTHER_SIDE.md` (line 585)  
**Proposed status:** open  
**Justification:** Standing waves are computationally proven, and the laser regime (Σγ < 0) provides a physical realization, but direct experimental measurement has not been attempted.

---

### OQ-301

**Question:** The hidden symmetry Q: 12 of 26 parity-breaking Hamiltonians still preserve the palindrome. What is Q?  
**Source:** `hypotheses/THE_OTHER_SIDE.md` (line 601)  
**Proposed status:** partially-resolved  
**Justification:** NON_HEISENBERG_PALINDROME.md identifies three conjugation operator families (uniform P1/P4, alternating M1xM2xM1, non-local entangled Pi). These are Q candidates. However, the explicit mapping from the 12 parity-breaking preservers onto these families has not been constructed.

---

### OQ-312

**Question:** How robust is the crossing-time correlation to local noise at A and B independently?  
**Source:** `hypotheses/TIME_AS_CROSSING_RATE.md` (line 370)  
**Proposed status:** resolved  
**Justification:** TIME_AS_CROSSING_RATE.md reports "Answered (2026-03-06)": STAR_TOPOLOGY_OBSERVERS.md gamma_A vs gamma_B scan shows receiver noise is far more destructive than sender noise. Robust to sender noise, fragile to receiver noise.

---

### OQ-317

**Question:** Antiferromagnet test: conditions 1 and 3 are problematic. Heisenberg exchange is SYMMETRIC under sublattice swap (not antisymmetric).  
**Source:** `hypotheses/UNIVERSAL_PALINDROME_CONDITION.md` (line 281)  
**Proposed status:** open  
**Justification:** The three universal conditions may not be met for antiferromagnets. The correct Q and "antisymmetry" for magnetic systems are unidentified.

---

### OQ-321

**Question:** V-Effect scaling with N. Neural: 0+0=6 (N=10), 0+0=48 (N=20). How does the number of V-Effect frequencies scale with N?  
**Source:** `hypotheses/UNIVERSAL_PALINDROME_CONDITION.md` (line 311)  
**Proposed status:** needs-human  
**Justification:** Overlaps with OQ-068 (same question about V-Effect frequency scaling with N, different source files: neural doc vs UNIVERSAL_PALINDROME_CONDITION). OQ-321 adds neural-specific data points. Candidate for merge with OQ-068.

---

### OQ-326

**Question:** Either we live in a simulation (the external is a simulator)  
**Source:** `hypotheses/WAVES_THAT_HEAR_THEMSELVES.md` (line 58)  
**Proposed status:** needs-human  
**Justification:** Fragment of a philosophical speculation about why coupling exists. Borderline between "open" (it's a genuine philosophical question the repo poses) and "obsolete" (it's unfalsifiable within the current framework). The framework explicitly acknowledges: "J is an input parameter, not a derived quantity." Needs human judgment on whether this belongs in an inventory of technical open questions.

---

## Duplicate / overlap clusters

**Cluster 1: 2-qubit Concurrence hardware validation**  
OQ-062 (WEAKNESSES_OPEN_QUESTIONS.md), OQ-097 (COCKPIT_UNIVERSALITY.md §5), OQ-104 (COCKPIT_UNIVERSALITY.md §6)  
Three entries across two source files, all asking for the same missing experiment. Recommend merging into OQ-062 as the canonical entry (it contains the most context: "most important missing validation").

**Cluster 2: Sacrifice-zone hardware gap**  
OQ-043 = OQ-047 (identical text from overlapping sections in WEAKNESSES_OPEN_QUESTIONS.md)  
Same pattern as OQ-041/OQ-045 and OQ-042/OQ-046 from Batch 6. Recommend removing OQ-047.

**Cluster 3: V-Effect frequency scaling**  
OQ-068 (V_EFFECT_NEURAL.md), OQ-321 (UNIVERSAL_PALINDROME_CONDITION.md)  
Same question with partially overlapping data. Recommend merging into OQ-068 with OQ-321's neural data points added.

**Cluster 4: Cross-batch duplicates**  
OQ-229 duplicates OQ-228 (math-proof, same formula from same source file)  
OQ-294 duplicates OQ-012 (math-proof, same "14 combos" question from different source file)

---

## Patterns

1. **High resolution rate.** 17 of 62 entries (27%) are resolved, the highest rate of any batch so far. Many numerical-verification entries are self-documenting results (VERIFIED, ANSWERED, falsified) rather than open questions, which made classification straightforward.

2. **Falsified claims properly documented.** OQ-238, OQ-242 (BRIDGE_PROTOCOL FTL claims), OQ-285 (inverted HO) are all correctly marked FALLEN/falsified in source. The repo's intellectual honesty convention works well for this tag.

3. **Concurrence hardware gap is the dominant open cluster.** OQ-062, OQ-097, OQ-104 (plus OQ-031 from Batch 5) all point to the same missing experiment: 2-qubit state tomography on hardware. This is the single most repeated open question across the entire inventory.

4. **Many entries are limitations, not questions.** OQ-092 (N=11 largest tested), OQ-093 (three pair types), OQ-099 (N=5 largest), OQ-181 (N=8 not tested), OQ-158 (feedback loop < 15%) describe scope boundaries rather than asking something. They're correctly classified as "open" since the limitation is acknowledged but not resolved.

5. **PAIR_BREAKING_AT_THE_HORIZON cluster.** OQ-280, OQ-281, OQ-283, OQ-285, OQ-290 form a coherent cluster from a single hypothesis file. OQ-280 and OQ-285 are resolved (listing completed, HO falsified); OQ-281, OQ-283, OQ-290 remain open (temperature scaling, mass identification, N-independence artifact).
