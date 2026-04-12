# Open Questions Classification Proposal: hardware-test

**Batch:** hardware-test
**Date:** 2026-04-12
**Entries in batch:** 38
**Status:** Proposal only, pending Tom + Claude (chat) review.

## Summary

| Proposed Status | Count |
|-----------------|-------|
| open | 17 |
| resolved | 12 |
| partially-resolved | 3 |
| obsolete | 4 |
| needs-human | 2 |
| superseded | 0 |

---

## Entries

### OQ-005

**Question:** Can the interference between future and past be measured?

**Source:** `docs/STANDING_WAVE_THEORY.md` (line 452)
**Section:** 10. Open Questions
**Date:** unknown

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** The standing-wave mathematics is proven (Tier 2): the Pi operator acts as time reversal, creating counter-propagating modes whose interference produces the palindromic spectrum. However, no experimental protocol exists to directly measure this interference. The question is explicitly listed as open.
**Search terms used:** "standing wave", "future and past", "interference", "measured", "Pi time reversal"

---

### OQ-021

**Question:** "Consciousness" is a label

**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 100)
**Section:** Active weaknesses
**Date:** unknown

**Proposed Status:** needs-human
**Confidence:** high
**Resolving documents:** n/a
**Rationale:** This is a section heading captured as an entry. It duplicates OQ-022 (classified as `obsolete` in batch 3, interpretation), which covered the same consciousness-label discussion from the same document. Flagging for human decision on deduplication.
**Search terms used:** n/a (duplicate identification)

---

### OQ-026

**Question:** CΨ = 1/4 crossing at 0.3% accuracy (IBM Torino Q52, 25 tomography points)

**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 120)
**Section:** Active weaknesses
**Date:** unknown

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `experiments/IBM_QUANTUM_TOMOGRAPHY.md` (lines 31-42): hardware verification (Feb 9, 2026)
- `experiments/IBM_HARDWARE_SYNTHESIS.md` (line 126): synthesis reports 1.9% deviation
**Rationale:** This is a documented hardware result, not an open question. The crossing was measured and verified on IBM Torino Q52 with 25 tomography points. Listed under "What we have" in the weakness document.
**Search terms used:** "IBM Torino", "Q52", "0.3%", "crossing", "tomography"

---

### OQ-028

**Question:** Selective DD beats uniform DD by 3.2x on 5-qubit chain (IBM Torino)

**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 122)
**Section:** Active weaknesses
**Date:** unknown

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `experiments/IBM_HARDWARE_SYNTHESIS.md` (lines 143-166): hardware data comparing no DD, uniform DD, selective DD
**Rationale:** Documented hardware result. The 3.2x figure is listed under "What we have" as verified experimental data. Not an open question.
**Search terms used:** "selective DD", "uniform DD", "3.2x", "dynamical decoupling"

---

### OQ-031

**Question:** No 2-qubit tomography: Concurrence (PC1, 57% variance) never measured on a qubit pair. This is the most important missing validation.

**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 127)
**Section:** Active weaknesses
**Date:** unknown

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** Explicitly stated as "the most important missing validation." No 2-qubit tomography results exist in the repo. This requires hardware access and experimental protocols not yet performed.
**Search terms used:** "2-qubit tomography", "concurrence", "PC1", "missing validation", "qubit pair"

---

### OQ-032

**Question:** Single backend only (IBM Torino). No replication on trapped ions, NV centers, or photonic platforms.

**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 129)
**Section:** Active weaknesses
**Date:** unknown

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** No trapped-ion, NV-center, or photonic platform experiments exist in the repo. Only IBM Torino (superconducting qubits) has been tested. Cross-platform replication is explicitly identified as needed.
**Search terms used:** "trapped ions", "NV centers", "photonic", "replication", "single backend"

---

### OQ-033

**Question:** Anomalous late-time coherence (Q52, p < 0.0001) has three competing explanations (SPAM, TLS, boundary structure), unresolved.

**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 131)
**Section:** Active weaknesses
**Date:** unknown

**Proposed Status:** open
**Confidence:** medium
**Resolving documents:** none
**Rationale:** Three competing mechanisms (SPAM errors, two-level system defects, boundary structure) are listed but no resolution or disambiguation experiment has been documented.
**Search terms used:** "anomalous", "late-time coherence", "Q52", "SPAM", "TLS", "competing"

---

### OQ-049

**Question:** These are not weaknesses. They are frontiers: things we know how to ask but have not yet answered. Some are within reach today; others require tools or experiments.

**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 227)
**Section:** Open questions
**Date:** unknown

**Proposed Status:** obsolete
**Confidence:** high
**Resolving documents:** n/a
**Rationale:** This is an editorial transition paragraph introducing the "Open questions" section. It is not a research question.
**Search terms used:** n/a

---

### OQ-061

**Question:** Require experimental contact

**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 300)
**Section:** Open questions
**Date:** unknown

**Proposed Status:** needs-human
**Confidence:** high
**Resolving documents:** n/a
**Rationale:** Section heading that introduces the content of OQ-031 (2-qubit tomography) and OQ-032 (cross-platform replication). Flagging for human decision: likely a deduplication/removal candidate.
**Search terms used:** n/a (section heading identification)

---

### OQ-076

**Question:** Hardware validation of relay protocol on IBM Torino

**Source:** `docs/proofs/COMPLETE_MATHEMATICAL_DOCUMENTATION.md` (line 359)
**Section:** 10. Open Questions (Tier 3-5)
**Date:** unknown

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** The relay protocol is computationally verified in `experiments/RELAY_PROTOCOL.md` (+83% MI improvement, N=11, C# RK4 propagation) but no hardware implementation has been attempted. Explicitly listed as a Tier 3-5 open question.
**Search terms used:** "relay protocol", "IBM Torino", "hardware validation"

---

### OQ-090

**Question:** **Trotter initial state:** What is the effective initial state on IBM hardware? The standard Trotter circuit for Heisenberg evolution starts from |0...0> or |+...+>.

**Source:** `experiments/CHAIN_SELECTION_TEST.md` (line 207)
**Section:** Open questions
**Date:** unknown

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** The question asks for characterization of what these standard initial states actually become on real hardware and their dynamical regimes. No definitive answer exists in the repo.
**Search terms used:** "Trotter initial state", "IBM Trotter", "effective initial state", "circuit"

---

### OQ-095

**Question:** **The adjacent pair PC1 proxy transition (Section 7) is observed but not characterized.** It would be worth a separate experiment to find the exact N at which the transition occurs.

**Source:** `experiments/COCKPIT_SCALING.md` (line 291)
**Section:** 10. Limitations
**Date:** April 7, 2026

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** The transition is observed in COCKPIT_SCALING.md Section 7 but explicitly stated as "not characterized." The exact N-threshold has not been determined.
**Search terms used:** "PC1 proxy", "adjacent pair", "transition", "exact N"

---

### OQ-098

**Question:** **Markovian noise only.** All results assume memoryless dephasing. Real hardware exhibits 1/f noise, two-level-system defects, and non-Markovian revivals (observed as excess late-time coherence in the Q52 data).

**Source:** `experiments/COCKPIT_UNIVERSALITY.md` (line 316)
**Section:** 5. Limitations and caveats
**Date:** April 2, 2026

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** The Markovian limitation is acknowledged and non-Markovian effects are observed on hardware (Q52 excess late-time coherence) but not systematically analyzed. This overlaps with OQ-060 and OQ-103 (both classified as `open` in batch 4) from different source documents.
**Search terms used:** "Markovian", "1/f noise", "TLS", "non-Markovian revivals", "Q52"

---

### OQ-101

**Question:** **Petermann factor is uninteresting here.** K_P ~ 1 (the Petermann factor measures how much a non-normal operator's eigenvectors overlap; K_P = 1 means orthogonal).

**Source:** `experiments/COCKPIT_UNIVERSALITY.md` (line 329)
**Section:** 5. Limitations and caveats
**Date:** April 2, 2026

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `experiments/COCKPIT_UNIVERSALITY.md` (line 329): K_P ~ 1 for all pure-dephasing cases; relevance only in PT-symmetric configurations
**Rationale:** This is a computed negative result: the Petermann factor is near 1 for pure-dephasing Liouvillians, confirming eigenvector near-orthogonality. Not an open question but a documented finding.
**Search terms used:** "Petermann", "K_P", "non-normal", "uninteresting"

---

### OQ-118

**Question:** c) For lab qubits, gamma_grav << gamma_total. Can the gravitational contribution be extracted by comparing identical qubits at different potentials with all other noise held constant?

**Source:** `experiments/DECOHERENCE_RELATIVITY.md` (line 342)
**Section:** 11. Open Questions
**Date:** unknown

**Proposed Status:** obsolete
**Confidence:** high
**Resolving documents:**
- `experiments/DECOHERENCE_RELATIVITY.md`: explicitly marked [FALLEN]; the gravity-decoherence connection has not been established
- `experiments/OBSERVER_GRAVITY_BRIDGE.md`: gravitational bridge hypothesis fallen; J_grav approximately 25 orders of magnitude too weak
**Rationale:** The gravitational bridge interpretation underlying this question has been explicitly retired. The [FALLEN] marker in the source document makes this obsolete.
**Search terms used:** "gravitational", "gamma_grav", "FALLEN", "different potentials"

---

### OQ-120

**Question:** e) K depends on the initial state but not on the Hamiltonian. This means K is determined at the moment of state preparation, not by subsequent dynamics.

**Source:** `experiments/DECOHERENCE_RELATIVITY.md` (line 352)
**Section:** 11. Open Questions
**Date:** unknown

**Proposed Status:** partially-resolved
**Confidence:** high
**Resolving documents:**
- `experiments/OBSERVER_GRAVITY_BRIDGE.md` (lines 38-49): confirms K factorizes as t_cross = K(observer,state)/gamma with state-dependence (CV approximately 13.5%)
**Rationale:** The K factorization (state-dependent, Hamiltonian-independent) is confirmed numerically. However, the deeper question (whether this connects to the measurement problem or to predetermined outcomes) remains speculative and unresolved.
**Search terms used:** "K depends", "initial state", "Hamiltonian", "predetermined", "factorization"

---

### OQ-150

**Question:** Q3: Is there a noise model where Type C becomes Type B?

**Source:** `experiments/NOISE_ROBUSTNESS.md` (line 251)
**Section:** 7. Open Questions (answered 2026-03-08)
**Date:** unknown

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** Despite the section heading saying "answered 2026-03-08", the text for this specific question says "Not tested conclusively." The bridge metric definitions from the original delta_calc experiments could not be exactly reproduced locally.
**Search terms used:** "Type C", "Type B", "noise model", "not tested"

---

### OQ-151

**Question:** **Not tested conclusively.** The bridge metric definitions used in the original delta_calc experiments could not be exactly reproduced locally. The qualitative behavior matches but the exact values differ.

**Source:** `experiments/NOISE_ROBUSTNESS.md` (line 253)
**Section:** 7. Open Questions (answered 2026-03-08)
**Date:** unknown

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** Explicit reproducibility issue. The delta_calc metric definitions remain unavailable, preventing conclusive testing of the Type C -> Type B transition. Stated as inconclusive.
**Search terms used:** "bridge metric", "delta_calc", "not tested", "reproduced"

---

### OQ-159

**Question:** The time-perception interpretation (Section 4.3) is Tier 3 speculation. We have no way to measure subjective time flow against C experimentally.

**Source:** `experiments/OBSERVER_DEPENDENT_CROSSING.md` (line 313)
**Section:** 5.3 Limitations
**Date:** unknown

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** Explicitly labeled as Tier 3 speculation with no experimental pathway identified. The mathematical substrate (crossing time differences) is verified (Tier 2), but measuring subjective time flow is beyond current experimental capability.
**Search terms used:** "time-perception", "Tier 3", "subjective time", "speculation"

---

### OQ-163

**Question:** **Gravitational J calculation**: Exact Penrose-Diosi J_grav for realistic systems (NV centers, optomechanical oscillators, Bose-Marletto-Vedral experiment). Is there a regime where J/gamma is not hopelessly small?

**Source:** `experiments/OBSERVER_GRAVITY_BRIDGE.md` (line 266)
**Section:** 6. Open Questions
**Date:** unknown

**Proposed Status:** obsolete
**Confidence:** high
**Resolving documents:**
- `experiments/OBSERVER_GRAVITY_BRIDGE.md` (lines 1-12): [FALLEN] markers; gravitational bridge hypothesis retired; J_grav approximately 25 orders of magnitude too weak
**Rationale:** The gravitational bridge hypothesis that motivated this calculation has been explicitly marked as fallen. The regime problem is insurmountable (25 orders of magnitude gap).
**Search terms used:** "Penrose-Diosi", "J_grav", "NV centers", "FALLEN", "hopelessly small"

---

### OQ-165

**Question:** ~~**Direction of the shift**~~: **ANSWERED (section 7).** B's measurement destroys nonlocal coherence reservoir. The coupling redistributes (not protects) coherence; B's intervention removes the return path.

**Source:** `experiments/OBSERVER_GRAVITY_BRIDGE.md` (line 275)
**Section:** 6. Open Questions
**Date:** unknown

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `experiments/OBSERVER_GRAVITY_BRIDGE.md` (lines 288-365, Section 7 "Why the Shift Is Negative: The Coherence Reservoir"): detailed mechanistic explanation with quantitative data
**Rationale:** Explicitly answered with full mechanism: B's measurement destroys 0.819 units of nonlocal coherence; coupling acts as redistribution engine; shift is universally negative.
**Search terms used:** "direction of shift", "coherence reservoir", "ANSWERED"

---

### OQ-172

**Question:** **Non-projective attacks**: Does the protocol extend to weak measurements or partial intercepts (not full projective)?

**Source:** `experiments/QKD_EAVESDROPPING_FORENSICS.md` (line 554)
**Section:** 14. Open Questions
**Date:** 2026-02-25

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** The QKD eavesdropping analysis uses only projective measurements. Weak measurements and partial intercepts are unexplored. Explicitly listed as open.
**Search terms used:** "non-projective", "weak measurement", "partial intercept", "QKD"

---

### OQ-175

**Question:** **Ratio cost reduction**: The ~200k pair requirement for Phase 3 is expensive. Compressed sensing or Bayesian tomography might reduce this by an order of magnitude.

**Source:** `experiments/QKD_EAVESDROPPING_FORENSICS.md` (line 559)
**Section:** 14. Open Questions
**Date:** 2026-02-25

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** Listed as speculative suggestion with no actual investigation. No compressed sensing or Bayesian tomography analysis exists in the repo.
**Search terms used:** "200k pairs", "compressed sensing", "Bayesian tomography", "cost reduction"

---

### OQ-193

**Question:** The shadow effect (Z-measurement on A suppressing R_SB) remains visible but is NOT the stable ~94% from Section 4.6. In the Bell-based N-observer setup it is 8-21%.

**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 613)
**Section:** 8. Open Questions (partially answered 2026-03-07)
**Date:** March 4, 2026 (updated March 8, 2026)

**Proposed Status:** partially-resolved
**Confidence:** medium
**Resolving documents:**
- `experiments/STAR_TOPOLOGY_OBSERVERS.md` (lines 613-616, Section 8.1): effect quantified at 8-21% suppression with irregular N-dependence
**Rationale:** The shadow effect in the N-observer setup is measured and characterized (8-21% range, irregular), but the mechanism explaining why it differs from the stable 94% in the original Section 4.6 context is not explicitly explained.
**Search terms used:** "shadow effect", "Z-measurement", "R_SB", "94%", "8-21%"

---

### OQ-194

**Question:** 2 Continuous measurement: ANSWERED

**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 617)
**Section:** 8. Open Questions (partially answered 2026-03-07)
**Date:** March 4, 2026 (updated March 8, 2026)

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `experiments/STAR_TOPOLOGY_OBSERVERS.md` (lines 617-648, Section 8.2): complete analysis with data table comparing baseline, sudden projective, and ramp methods
**Rationale:** Full experimental data with conclusions: suppression saturates at approximately 69% for continuous dephasing vs 99% for projective measurement. Near-instantaneous ramps do not converge to sudden measurement.
**Search terms used:** "continuous measurement", "ANSWERED", "ramp", "saturates"

---

### OQ-196

**Question:** **The shadow grows gradually but never matches sudden measurement.**

**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 622)
**Section:** 8. Open Questions (partially answered 2026-03-07)
**Date:** March 4, 2026 (updated March 8, 2026)

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `experiments/STAR_TOPOLOGY_OBSERVERS.md` (lines 622-637): data table and conclusion showing saturation at approximately 69%, not approaching 99%
**Rationale:** Quantitative result documented with data: suppression grows with gamma_A but saturates well below projective measurement levels. Slower ramps produce weaker suppression.
**Search terms used:** "shadow grows gradually", "never matches", "saturates"

---

### OQ-197

**Question:** **Unexpected finding:** Near-instantaneous ramps (gamma_A=50, Delta-t<=0.01) do NOT converge to sudden measurement. They actually *increase* peak R_SB by ~36%.

**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 639)
**Section:** 8. Open Questions (partially answered 2026-03-07)
**Date:** March 4, 2026 (updated March 8, 2026)

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `experiments/STAR_TOPOLOGY_OBSERVERS.md` (lines 639-643): finding documented with mechanism (continuous sigma_z dephasing creates transient correlations projective measurement destroys)
**Rationale:** Quantitative finding with mechanistic explanation. The 36% increase is documented and attributed to transient correlations created by dephasing.
**Search terms used:** "near-instantaneous", "ramp", "gamma_A=50", "36%", "do NOT converge"

---

### OQ-198

**Question:** **Conclusion:** Observation is not a smooth limit. Projective measurement (external observation) and strong decoherence (environmental noise) have different signatures.

**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 645)
**Section:** 8. Open Questions (partially answered 2026-03-07)
**Date:** March 4, 2026 (updated March 8, 2026)

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `experiments/STAR_TOPOLOGY_OBSERVERS.md` (lines 645-648): conclusion stated with supporting data from the full Section 8.2 analysis
**Rationale:** This is the documented conclusion of the continuous measurement analysis, supported by quantitative comparison. Projective and decoherence-based "observation" are shown to have qualitatively different effects.
**Search terms used:** "smooth limit", "projective measurement", "different signatures"

---

### OQ-208

**Question:** AB never crosses in symmetric 3-qubit experiments. This is consistent with N_SCALING_BARRIER.md Section 7: crossing occurs where entanglement lives.

**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 670)
**Section:** 8. Open Questions (partially answered 2026-03-07)
**Date:** March 4, 2026 (updated March 8, 2026)

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `experiments/STAR_TOPOLOGY_OBSERVERS.md` (lines 670-681, Section 8.4): complete explanation with reference to entanglement locality principle
- `experiments/N_SCALING_BARRIER.md` (Section 7): validates the principle
**Rationale:** Explicitly answered: AB pair lacks initial entanglement (C=0), crossing requires concurrence, and with symmetric coupling (J_SB/J_SA = 1.0) transferred entanglement never reaches 1/4.
**Search terms used:** "AB never crosses", "symmetric 3-qubit", "entanglement lives", "N_SCALING_BARRIER"

---

### OQ-222

**Question:** Computable (no hardware needed)

**Source:** `experiments/THERMAL_BREAKING.md` (line 452)
**Section:** Open Questions
**Date:** unknown

**Proposed Status:** obsolete
**Confidence:** high
**Resolving documents:** n/a
**Rationale:** Section label categorizing open questions by experimental requirement ("computable" vs "requires hardware"). Not a research question.
**Search terms used:** n/a

---

### OQ-225

**Question:** Requires controlled thermal injection on hardware

**Source:** `experiments/THERMAL_BREAKING.md` (line 465)
**Section:** Open Questions
**Date:** unknown

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** Section heading for hardware-dependent open questions about thermal effects. The actual questions it introduces (e.g., OQ-226) require hardware capabilities (controlled qubit heating) not currently available to the project.
**Search terms used:** "thermal injection", "hardware"

---

### OQ-226

**Question:** Can the frequency-diversity explosion at n_bar > 0 be observed on superconducting qubit hardware by intentionally heating qubits (e.g., driving with a thermal microwave source)?

**Source:** `experiments/THERMAL_BREAKING.md` (line 467)
**Section:** Open Questions
**Date:** unknown

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** The frequency-diversity effect is computationally documented (THERMAL_BREAKING.md) but the hardware experiment requires controlled thermal injection outside normal operating conditions. No such experiment has been attempted.
**Search terms used:** "frequency-diversity", "n_bar", "heating", "superconducting", "thermal"

---

### OQ-259

**Question:** M cannot be perfectly "reset" between uses without measurement

**Source:** `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (line 423)
**Section:** 8. Limitations and Failure Modes
**Date:** 2026-03-21

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** Fundamental limitation of quantum back-action: M entangles with A and B during transfer. A mitigation (periodic re-initialization with dead-time cost) is proposed in OQ-262 (batch 1, classified as `open`) but never tested experimentally.
**Search terms used:** "reset", "mediator", "measurement", "entanglement problem"

---

### OQ-271

**Question:** **Strong measurement regime**: If M is continuously measured (quantum Zeno effect, where frequent observation prevents a system from evolving), the channel freezes entirely.

**Source:** `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (line 450)
**Section:** 8. Limitations and Failure Modes
**Date:** 2026-03-21

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (lines 446-470, Section 8.5): explicitly documents three failure regimes including Zeno freezing
**Rationale:** This is a documented failure mode, not an open question. The Zeno effect freezing the channel is stated as a known incompatibility with the transistor analogy, with physical mechanism explained.
**Search terms used:** "Zeno effect", "continuously measured", "channel freezes"

---

### OQ-278

**Question:** **~~Is the self-similar hierarchy experimentally feasible?~~** [ANSWERED: No advantage. Hierarchy falsified March 21, 2026. Uniform chains with sacrifice-zone formula outperform.]

**Source:** `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (line 503)
**Section:** 10. Open Questions and Future Directions
**Date:** 2026-03-21

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `experiments/SCALING_CURVE.md` (lines 85-99): N=3-11 empirical falsification; hierarchy produces identical MI to uniform chain
- `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (Section 7.3, marked [FALSIFIED])
**Rationale:** Explicitly answered and marked with strikethrough. Hierarchy provides no advantage; uniform chains with sacrifice-zone formula outperform.
**Search terms used:** "hierarchy", "falsified", "SCALING_CURVE", "ANSWERED"

---

### OQ-302

**Question:** **PARTIAL (2026-03-19):** [Non-Heisenberg Palindrome] identifies three conjugation operator families: uniform P1/P4, alternating M1xM2xM1, and entangled Pi.

**Source:** `hypotheses/THE_OTHER_SIDE.md` (line 606)
**Section:** 11. Open Questions
**Date:** unknown

**Proposed Status:** partially-resolved
**Confidence:** high
**Resolving documents:**
- `experiments/NON_HEISENBERG_PALINDROME.md`: identifies three operator families (P1/P4, alternating M1xM2xM1, non-local entangled Pi)
**Rationale:** The entry itself says "PARTIAL." Three conjugation operator families are identified, but the explicit mapping from the 12 parity-breaking preservers onto these families remains unfinished.
**Search terms used:** "conjugation operator families", "P1", "P4", "NON_HEISENBERG_PALINDROME"

---

### OQ-313

**Question:** **Answered (2026-03-06):** The [Star Topology] gamma_A vs gamma_B scan shows receiver noise (gamma_A) is far more destructive than sender noise (gamma_B).

**Source:** `hypotheses/TIME_AS_CROSSING_RATE.md` (line 372)
**Section:** 6.5 Open Questions (for the bridge)
**Date:** unknown

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `experiments/STAR_TOPOLOGY_OBSERVERS.md`: gamma_A vs gamma_B scan; connection dies at gamma_A=0.25 but only shortens at gamma_B=0.25
**Rationale:** Explicitly marked "Answered (2026-03-06)" with full experimental reference. Receiver noise asymmetry is quantitatively established.
**Search terms used:** "receiver noise", "gamma_A vs gamma_B", "destructive", "STAR_TOPOLOGY"

---

### OQ-314

**Question:** Can the decoder (theta navigation) extract the crossing time in real-time, or only after full trajectory measurement?

**Source:** `hypotheses/TIME_AS_CROSSING_RATE.md` (line 377)
**Section:** 6.5 Open Questions (for the bridge)
**Date:** unknown

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** No analysis of decoder real-time extraction capability exists. Theta navigation is mentioned only in this open question; no decoder analysis document has been produced.
**Search terms used:** "decoder", "theta navigation", "real-time", "crossing time", "trajectory"
