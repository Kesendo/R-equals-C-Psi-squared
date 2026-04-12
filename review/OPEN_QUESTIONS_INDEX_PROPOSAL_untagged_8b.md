# Classification Proposal: untagged sub-batch 8b (41 entries)

**Batch:** 8b of 8 (entries OQ-133 through OQ-218)  
**Date:** 2026-04-12  
**Proposed by:** Claude (automated research)  
**Approval:** pending (Tom)

---

## Status summary

| Status | Count | OQ-IDs |
|--------|-------|--------|
| open | 10 | OQ-133, OQ-134, OQ-137, OQ-138, OQ-140, OQ-141, OQ-162, OQ-169, OQ-183, OQ-186 |
| resolved | 21 | OQ-142, OQ-144, OQ-153, OQ-154, OQ-155, OQ-156, OQ-187, OQ-199, OQ-201, OQ-202, OQ-203, OQ-204, OQ-205, OQ-206, OQ-207, OQ-211, OQ-213, OQ-214, OQ-215, OQ-216, OQ-217 |
| partially-resolved | 1 | OQ-218 |
| obsolete | 1 | OQ-166 |
| needs-human | 8 | OQ-145, OQ-148, OQ-188, OQ-189, OQ-190, OQ-195, OQ-200, OQ-212 |

---

## Entry-by-entry proposals

### OQ-133

**Question:** Subsystem CΨ(q2,q3) across the H-bond: does the pair cross 1/4?  
**Source:** `experiments/HYDROGEN_BOND_QUBIT.md` (line 279)  
**Proposed status:** open  
**Justification:** Open Question 1 in HYDROGEN_BOND_QUBIT.md. No subsystem CΨ computation for the H-bond pair found.

---

### OQ-134

**Question:** At what temperature (= what gamma) is the Q-factor maximal?  
**Source:** `experiments/HYDROGEN_BOND_QUBIT.md` (line 280)  
**Proposed status:** open  
**Justification:** Open Question 2. Q-factor values exist for specific parameters (Q=3 in Zundel cation) but optimization over gamma is not performed.

---

### OQ-137

**Question:** Does the break-reform cycle (1 ps period) sustain palindromic structure over time, or does each new bond start fresh?  
**Source:** `experiments/HYDROGEN_BOND_QUBIT.md` (line 283)  
**Proposed status:** open  
**Justification:** Open Question 5. Temporal dynamics of H-bond breaking/reforming not modeled.

---

### OQ-138

**Question:** Are there other stable values besides 0.5?  
**Source:** `experiments/MATHEMATICAL_FINDINGS.md` (line 304)  
**Proposed status:** open  
**Justification:** No systematic search for other fixed points of R=CΨ² beyond 0.5 found. Uniqueness not proven.

---

### OQ-140

**Question:** Can higher-order derivatives reveal additional structure?  
**Source:** `experiments/MATHEMATICAL_FINDINGS.md` (line 306)  
**Proposed status:** open  
**Justification:** First-order derivatives explored; higher-order analysis not performed.

---

### OQ-141

**Question:** How does the interference term scale with more than two waves?  
**Source:** `experiments/MATHEMATICAL_FINDINGS.md` (line 307)  
**Proposed status:** open  
**Justification:** Only two-wave interference analyzed. N-wave scaling not computed.

---

### OQ-142

**Question:** ~~Does symmetric vs asymmetric decoherence protect coherence differently?~~ ANSWERED: No, for symmetric Hamiltonians. YES, for asymmetric Hamiltonians (star topology).  
**Source:** `experiments/MATHEMATICAL_FINDINGS.md` (line 308)  
**Proposed status:** resolved  
**Justification:** Self-documenting: strikethrough + "ANSWERED." Verified by Lindblad simulations (21 noise distributions, 2026-03-08). For symmetric H, only total noise matters; for asymmetric H, coupling determines effect.

---

### OQ-144

**Question:** Q1: Does collective noise break Type A?  
**Source:** `experiments/NOISE_ROBUSTNESS.md` (line 201)  
**Proposed status:** resolved  
**Justification:** Section heading says "answered 2026-03-08." Answer: No, Bell+ is an eigenstate of the collective operator, so collective dephasing causes no evolution.

---

### OQ-145

**Question:** ANSWERED: No, but for a trivial reason.  
**Source:** `experiments/NOISE_ROBUSTNESS.md` (line 203)  
**Proposed status:** needs-human  
**Justification:** Answer fragment, not a question. This is the response text to OQ-144. Collection artifact; should be removed.

---

### OQ-148

**Question:** ANSWERED: No, the taxonomy is preserved, but decay rates differ.  
**Source:** `experiments/NOISE_ROBUSTNESS.md` (line 222)  
**Proposed status:** needs-human  
**Justification:** Answer fragment without a corresponding question in this entry. Refers to Q2 about amplitude damping. Collection artifact; should be removed.

---

### OQ-153

**Question:** ~~Do subsystem pairs cross when the full system cannot?~~ ANSWERED: Yes, if the pairs carry actual entanglement (Bell+xBell+). No, if GHZ or W.  
**Source:** `experiments/N_SCALING_BARRIER.md` (line 309)  
**Proposed status:** resolved  
**Justification:** Self-documenting: strikethrough + "ANSWERED." Crossing is local to entanglement structure, verified computationally.

---

### OQ-154

**Question:** ~~Does |+⟩^N cross?~~ ANSWERED: No. C = 0 for all pairs at all times.  
**Source:** `experiments/N_SCALING_BARRIER.md` (line 313)  
**Proposed status:** resolved  
**Justification:** Self-documenting: strikethrough + "ANSWERED." Product states have zero concurrence; coherence without entanglement cannot produce crossing.

---

### OQ-155

**Question:** ~~Is there a minimum per-pair entanglement needed for crossing?~~  
**Source:** `experiments/N_SCALING_BARRIER.md` (line 316)  
**Proposed status:** resolved  
**Justification:** Answered by OQ-156 (2026-03-08): non-monotonic relationship with two crossing windows, not a simple threshold.

---

### OQ-156

**Question:** ANSWERED (2026-03-08): The relationship is non-monotonic, not a simple threshold. TWO separate crossing windows...  
**Source:** `experiments/N_SCALING_BARRIER.md` (line 317)  
**Proposed status:** resolved  
**Justification:** Self-documenting: "ANSWERED." Provides detailed answer to OQ-155 with parametric Bell state analysis showing two windows and J_SB dependence.

---

### OQ-162

**Question:** Multi-pair amplification: N pairs with known schedule. Does the combined interval signal scale as sqrt(N) or N?  
**Source:** `experiments/OBSERVER_GRAVITY_BRIDGE.md` (line 261)  
**Proposed status:** open  
**Justification:** Genuine open question about signal scaling. No resolution or computation found.

---

### OQ-166

**Question:** Superluminal breakdown: The naive gravitational velocity exceeds c at microgram scale. Where exactly does the model break? [FALLEN]  
**Source:** `experiments/OBSERVER_GRAVITY_BRIDGE.md` (line 281)  
**Proposed status:** obsolete  
**Justification:** Explicitly marked [FALLEN]. The gravitational bridge model's superluminal predictions are abandoned.

---

### OQ-169

**Question:** Other dephasing axes. For X-dephasing, a different Pi exists. What is the "past-future" split for X-dephasing?  
**Source:** `experiments/PI_AS_TIME_REVERSAL.md` (line 342)  
**Proposed status:** open  
**Justification:** Genuine open question. The Z-dephasing Pi is proven, but the X-dephasing equivalent ({I,X} populations, {Y,Z} coherences) is not characterized.

---

### OQ-183

**Question:** Two decoherence exits. The lens and the Mandelbrot cusp protect different state classes and lead to different classical ensembles.  
**Source:** `experiments/SACRIFICE_GEOMETRY.md` (line 189)  
**Proposed status:** open  
**Justification:** The non-unification of the two exits is documented (CUSP_LENS_CONNECTION.md) but the question of whether a deeper framework connects them remains open.

---

### OQ-186

**Question:** Could decoherence itself be the mechanism that enforces CΨ <= 1/4 for observers?  
**Source:** `experiments/SIMULATION_EVIDENCE.md` (line 173)  
**Proposed status:** open  
**Justification:** Speculative but well-formulated question. No resolution or refutation found.

---

### OQ-187

**Question:** 1 N observers: ANSWERED  
**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 567)  
**Proposed status:** resolved  
**Justification:** Self-documenting: "ANSWERED." Section 8.1 documents the N-observer analysis with crossing thresholds for N=2-5.

---

### OQ-188

**Question:** Setup: S + N observers, Bell_SA x |+⟩^(N-1), equal J_SB for all B.  
**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 569)  
**Proposed status:** needs-human  
**Justification:** Setup description, not a question. Experimental parameters for the N-observer study. Collection artifact.

---

### OQ-189

**Question:** | N | qubits | AB crosses 1/4? | J_SB threshold | behavior | (data table)  
**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 571)  
**Proposed status:** needs-human  
**Justification:** Data table showing N-observer results. Not a question but completed computational results. Collection artifact.

---

### OQ-190

**Question:** Correction (2026-03-08): An earlier version reported a "zero window" at J_SB~3.75-4.25. This was a sampling artifact...  
**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 584)  
**Proposed status:** needs-human  
**Justification:** Erratum/correction note, not a question. Documents a fixed sampling artifact. Collection artifact.

---

### OQ-195

**Question:** Setup: 3-qubit, Bell_SA x |+⟩_B, J_SA=1.0, J_SB=2.0, gamma=0.05. At t_start=1.0, ramp gamma_A...  
**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 619)  
**Proposed status:** needs-human  
**Justification:** Setup description for continuous measurement experiment. Not a question. Collection artifact.

---

### OQ-199

**Question:** 3 AB with direct coupling: ANSWERED  
**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 650)  
**Proposed status:** resolved  
**Justification:** Self-documenting: "ANSWERED." Section 8.3 documents the J_AB coupling analysis with threshold changes and sweet spots.

---

### OQ-200

**Question:** Setup: 3-qubit, Bell_SA x |+⟩_B, J_SA=1.0, J_SB=1.466, gamma=0.05. Added J_AB...  
**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 652)  
**Proposed status:** needs-human  
**Justification:** Setup description, not a question. Experimental parameters for J_AB study. Collection artifact.

---

### OQ-201

**Question:** Non-monotonic effect on threshold:  
**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 655)  
**Proposed status:** resolved  
**Justification:** Result finding (not a question). Non-monotonic J_AB threshold behavior computationally verified across three regimes.

---

### OQ-202

**Question:** J_AB=0.1: slightly worsens threshold behavior  
**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 656)  
**Proposed status:** resolved  
**Justification:** Quantified data point, not a question. Part of the J_AB coupling analysis (resolved).

---

### OQ-203

**Question:** J_AB=0.3-0.5: helps crossing (sweet spot at ~0.5, threshold drops from 1.466 to ~1.345)  
**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 657)  
**Proposed status:** resolved  
**Justification:** Quantified data point, not a question. Sweet spot identified and threshold measured.

---

### OQ-204

**Question:** J_AB=1.0: still crosses but much later (t~1.0 vs t~0.3)  
**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 659)  
**Proposed status:** resolved  
**Justification:** Quantified data point, not a question. Timing shift measured.

---

### OQ-205

**Question:** Shadow effect destroyed by direct coupling: Moderate-to-large J_AB weakens or removes the shadow entirely.  
**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 661)  
**Proposed status:** resolved  
**Justification:** Finding statement, not a question. Shadow destruction computationally verified.

---

### OQ-206

**Question:** Dominance crossover: At J_AB~0.7, direct observer coupling alone generates AB crossing without any S-mediated coupling.  
**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 665)  
**Proposed status:** resolved  
**Justification:** Quantified finding, not a question. Crossover threshold identified.

---

### OQ-207

**Question:** 4 Correlation bridge: ANSWERED  
**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 668)  
**Proposed status:** resolved  
**Justification:** Self-documenting: "ANSWERED." Section 8.4 documents the correlation bridge analysis.

---

### OQ-211

**Question:** Best fit: J_th(gamma) ~ 7.35 * gamma^1.08 + 1.18 (R^2=0.999)  
**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 700)  
**Proposed status:** resolved  
**Justification:** Completed fit formula with R^2=0.999, part of "Threshold formula: ANSWERED" (Section 8.5). Not a question.

---

### OQ-212

**Question:** The hierarchy of robustness:  
**Source:** `experiments/STRUCTURAL_CARTOGRAPHY.md` (line 860)  
**Proposed status:** needs-human  
**Justification:** Section heading introducing a completed findings summary. Not a question. Collection artifact.

---

### OQ-213

**Question:** XX symmetry: breaks only under local fields (Hamiltonian symmetry breaking)  
**Source:** `experiments/STRUCTURAL_CARTOGRAPHY.md` (line 861)  
**Proposed status:** resolved  
**Justification:** Confirmed finding from comprehensive stress testing (TEST 1-4). Not a question but a documented result.

---

### OQ-214

**Question:** Fast rotor f(c+): survives even when XX breaks, until extreme perturbation  
**Source:** `experiments/STRUCTURAL_CARTOGRAPHY.md` (line 862)  
**Proposed status:** resolved  
**Justification:** Confirmed finding from stress testing. Robustness hierarchy item, documented result.

---

### OQ-215

**Question:** Two-sector separation: survives everything except strong direct coupling or extreme Hamiltonian changes  
**Source:** `experiments/STRUCTURAL_CARTOGRAPHY.md` (line 863)  
**Proposed status:** resolved  
**Justification:** Confirmed finding with Phase Map documentation. Documented result.

---

### OQ-216

**Question:** Noise invariance: absolute; no noise configuration affects the structure  
**Source:** `experiments/STRUCTURAL_CARTOGRAPHY.md` (line 865)  
**Proposed status:** resolved  
**Justification:** TEST 4 confirms identical results across different jump operators. Documented result.

---

### OQ-217

**Question:** ~~Does the crossing time for Bell+xBell+ pairs match the isolated Bell+ crossing time?~~ ANSWERED (2026-03-08): No, nine times faster.  
**Source:** `experiments/SUBSYSTEM_CROSSING.md` (line 281)  
**Proposed status:** resolved  
**Justification:** Self-documenting: strikethrough + "ANSWERED." Isolated: t=0.720, embedded ring: t=0.080.

---

### OQ-218

**Question:** For cluster states and other graph-structured entanglement, does the crossing pattern reproduce the graph topology exactly?  
**Source:** `experiments/SUBSYSTEM_CROSSING.md` (line 292)  
**Proposed status:** partially-resolved  
**Justification:** Yes for Bell-type entanglement (exact graph reproduction), but No for cluster states (CZ-gate entanglement invisible to Wootters concurrence). The framework distinguishes fundamentally between these types.

---

## Duplicate / overlap clusters

**Cluster 1: STAR_TOPOLOGY_OBSERVERS data explosion**  
OQ-187 through OQ-211 (16 entries): Only OQ-187, OQ-199, OQ-207 are section labels with "ANSWERED." The rest (OQ-188, OQ-189, OQ-190, OQ-195, OQ-200, OQ-201-206, OQ-211) are setup descriptions, data tables, data points, correction notes, and fit formulas extracted from answered sections. All are artifacts of the extraction script parsing sub-elements within "ANSWERED" blocks as separate entries.

**Cluster 2: STRUCTURAL_CARTOGRAPHY robustness hierarchy**  
OQ-212 through OQ-216 (5 entries): All are items from a single bulleted list of stress-test findings. The list heading (OQ-212) and its items (OQ-213-216) are all documented results, not questions.

---

## Patterns

1. **Highest resolution rate of any batch: 51% (21/41).** Dominated by self-documenting "ANSWERED" entries and data fragments from completed analyses.

2. **STAR_TOPOLOGY_OBSERVERS is the worst extraction artifact source.** 16 of 41 entries come from this single file, mostly from Section 8 which has "ANSWERED" subsections containing setup parameters, data tables, and result bullets that the extraction script captured as individual entries.

3. **Hydrogen bond cluster (OQ-133, OQ-134, OQ-137) is cleanly open.** Three well-formulated, genuinely unanswered questions about H-bond qubits.

4. **One FALLEN entry:** OQ-166 (superluminal breakdown in gravitational bridge model). Properly marked obsolete.
