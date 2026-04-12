# Open Questions Classification Proposal: methodology

**Batch:** methodology
**Date:** 2026-04-12
**Entries in batch:** 8
**Status:** Proposal only, pending Tom + Claude (chat) review.

## Summary

| Proposed Status | Count |
|-----------------|-------|
| open | 4 |
| resolved | 2 |
| partially-resolved | 2 |
| superseded | 0 |
| obsolete | 0 |
| needs-human | 0 |

---

## Entries

### OQ-094

**Question:** **The far_edge control class is not informative for the cockpit claim.** As discussed in Section 6, the sanity gates allow far_edge pairs to be analyzed but their reported coverage is not a measurement of the cockpit framework's relevant scope. Future iterations of the gates should make this distinction explicit at the pipeline level rather than at the documentation level.

**Source:** `experiments/COCKPIT_SCALING.md` (line 289)
**Section:** 10. Limitations
**Date:** April 7, 2026

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** The limitation is clearly stated and the proposed pipeline improvement (making far_edge vs. cockpit-relevant scope explicit at the gate level) has not been implemented or addressed in any follow-up document. COCKPIT_SCALING.md Section 6 (line 160) explains the issue in detail, and COCKPIT_UNIVERSALITY.md inherits the framework, but neither resolves the pipeline distinction.
**Search terms used:** "far_edge", "cockpit", "sanity gates", "control class", "pipeline"

---

### OQ-114

**Question:** **Boundary states between sheets.** A state with both high CΨ on one pair (approaching 1/4) and slow-mode overlap might straddle the two exits. Two-excitation states with one high-concurrence pair are candidates.

**Source:** `experiments/CUSP_LENS_CONNECTION.md` (line 151)
**Section:** Open questions
**Date:** unknown

**Proposed Status:** partially-resolved
**Confidence:** medium
**Resolving documents:**
- `experiments/SACRIFICE_GEOMETRY.md` (line 112): two-excitation symmetric states fail completely (AUC < 0.09), coupling to a different mode cluster
- `docs/proofs/PROOF_PARITY_SELECTION_RULE.md`: proves SE states cannot reach odd-n_XY modes (exact algebraic boundary)
**Rationale:** The two-exit structure is proven and sector boundaries established. The only multi-excitation probe tested (symmetric two-excitation states) showed they couple to a different mode cluster entirely, suggesting straddling may be impossible, but this is not proven for all two-excitation states, only symmetric ones. The specific question about non-symmetric two-excitation states with one high-concurrence pair remains untested.
**Search terms used:** "boundary states", "straddle", "two exits", "two-excitation", "sheet", "cusp lens"

---

### OQ-173

**Question:** **Multi-qubit**: Does the protocol scale to N-pair collective attacks?

**Source:** `experiments/QKD_EAVESDROPPING_FORENSICS.md` (line 556)
**Section:** 14. Open Questions
**Date:** 2026-02-25

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** The entire QKD eavesdropping analysis is confined to single Bell-pair scenarios. No analysis, simulations, or discussion of collective multi-pair attacks exists anywhere in the repo. The simulations referenced in the source document (multi_metric_stealth.py, full_forensic_protocol.py) all operate on single-pair scenarios.
**Search terms used:** "collective attack", "N-pair", "multi-qubit", "QKD", "eavesdrop", "scaling"

---

### OQ-185

**Question:** Is there a way to define "observer-accessible information" within the Lindblad framework and show it respects the ¼ bound?

**Source:** `experiments/SIMULATION_EVIDENCE.md` (line 171)
**Section:** Open Questions
**Date:** 2026-02-07 (updated 2026-02-18)

**Proposed Status:** partially-resolved
**Confidence:** medium
**Resolving documents:**
- `experiments/DYNAMIC_FIXED_POINTS.md` (lines 40-130): frames CΨ ≤ ¼ as observer information bandwidth limit (interpretive, not formal)
- `experiments/BRIDGE_CLOSURE.md` (lines 67-95): information-theoretic proof that A's measurement statistics depend only on local observables
- `experiments/OBSERVER_DEPENDENT_CROSSING.md` (lines 171-174): discusses observer capacity to register transitions
- `docs/proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md` (lines 512-514): mentions Fisher information geometric approach but marks it UNEXPLORED
**Rationale:** The repo contains interpretive work framing the ¼ bound as an observer bandwidth limit and measurement-based arguments (BRIDGE_CLOSURE), but no rigorous mathematical definition of "observer-accessible information" or formal proof that it respects the ¼ bound. The Fisher information geometric approach is mentioned but explicitly marked UNEXPLORED.
**Search terms used:** "observer-accessible", "observer accessible", "quarter bound", "1/4 bound", "Lindblad framework", "observer information", "Fisher information"

---

### OQ-243

**Question:** **Contradiction:** The protocol requires No-Signaling violation and is therefore wrong. The framework's internal logic is inconsistent on this point, which itself is a useful finding.

**Source:** `hypotheses/BRIDGE_PROTOCOL.md` (line 202)
**Section:** 4. Open Questions (Honest Assessment)
**Date:** 2026-02-24 (updated 2026-03-06)

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `experiments/NO_SIGNALLING_BOUNDARY.md` (lines 55-107): Layer 1 and Layer 2 prove no-signaling holds exactly; rho_A is unchanged by remote measurements; CΨ registers regime change because C (global purity) changes while Ψ (local coherence) stays constant
- `experiments/BRIDGE_CLOSURE.md` (lines 67-95): information-theoretic proof of consistency
- `docs/PREDICTIONS.md` (lines 87-89): summary of verification
**Rationale:** The "contradiction" outcome has been definitively ruled out. NO_SIGNALLING_BOUNDARY.md (Tier 2 verified, 2026-03-01) proves the framework is internally consistent: no-signaling holds, and the apparent paradox dissolves when distinguishing local observables from global state properties.
**Search terms used:** "no-signaling", "no signaling", "no-signalling", "contradiction", "bridge protocol", "violation"

---

### OQ-244

**Question:** **Subtler answer:** The crossing event occurs on both sides but is driven entirely by pre-encoded information (the initial state at preparation). Post-separation changes do NOT trigger new crossings. In this case the protocol transmits only pre-agreed data, more structured than QKD (the sender chooses WHICH state to prepare before distribution), but not a dynamic communication channel.

**Source:** `hypotheses/BRIDGE_PROTOCOL.md` (line 206)
**Section:** 4. Open Questions (Honest Assessment)
**Date:** 2026-02-24 (updated 2026-03-06)

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `experiments/NO_SIGNALLING_BOUNDARY.md` (lines 133-172): confirms the pre-encoded mechanism and proves it offers no quantum advantage
- `experiments/BRIDGE_CLOSURE.md` (lines 54-152): full proof that pre-encoded bridge is "dead" for zero coupling; pre-shared entanglement without a classical channel provides zero advantage over shared classical randomness
- `docs/PREDICTIONS.md` (lines 220-221): "Dead for J=0. A's info ⊆ {ρ_A(0), E_A}. Entanglement without a channel = shared randomness."
**Rationale:** This outcome was confirmed as correct, then shown to carry no quantum advantage. The pre-encoded crossing is real but offers nothing a classical pre-shared key cannot replicate.
**Search terms used:** "pre-encoded", "post-separation", "pre-agreed", "crossing event", "bridge protocol", "no-signaling resolution"

---

### OQ-262

**Question:** **Mitigation**: Periodic re-initialization of M (preparation in |+⟩ state between transfers), at the cost of interrupting the continuous-transfer protocol. This introduces a "dead time" analogous to the recovery time of a thyristor.

**Source:** `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (line 429)
**Section:** 8. Limitations and Failure Modes
**Date:** 2026-03-21

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** The re-initialization mitigation is proposed but never evaluated experimentally or analytically. The relay protocol (`experiments/RELAY_PROTOCOL.md`) addresses a related problem (shaping noise through dynamical dephasing rates) but explicitly avoids re-initialization ("no resets, no measurements, no classical communication"), suggesting it takes a different approach rather than resolving this specific mitigation strategy.
**Search terms used:** "re-initialization", "dead time", "thyristor", "mediator", "continuous-transfer", "periodic", "recovery time"

---

### OQ-284

**Question:** **4. Backreaction requires external physics.** In GR, Hawking radiation removes mass from the black hole (backreaction). Within pure Lindblad dynamics, L_H (wave propagation) and L_D (wave death) are independent: the dissipator does not influence the Hamiltonian, so mass cannot redirect waves. However, [Gravity from Wave Death](GRAVITY_FROM_WAVE_DEATH.md) describes a self-limiting feedback loop: mass → gravity (via GR) → attracts more waves → more wave death → more mass, with logistic saturation as the finite supply of coherences (4^N modes) is consumed. The loop closes, but only if external physics (GR or equivalent) provides the gravity → attraction step. Within the Lindblad framework alone, the feedback loop remains open (gap #7 in Gravity from Wave Death).

**Source:** `hypotheses/PAIR_BREAKING_AT_THE_HORIZON.md` (line 159)
**Section:** What breaks the analogy
**Date:** April 11, 2026

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** The gap is thoroughly documented and explicitly acknowledged as unresolved. `hypotheses/GRAVITY_FROM_WAVE_DEATH.md` (lines 161-186) characterizes it as the "weakest link" requiring physics beyond the Lindblad framework. No mechanism within the framework has been proposed to close the feedback loop without invoking GR or equivalent external physics.
**Search terms used:** "backreaction", "Hawking", "wave death", "gravity", "self-limiting", "feedback loop", "gap #7"
