# Open Questions Classification Proposal: closed-form

**Batch:** closed-form
**Date:** 2026-04-12
**Entries in batch:** 9
**Status:** Proposal only, pending Tom + Claude (chat) review.

## Summary

| Proposed Status | Count |
|-----------------|-------|
| open | 3 |
| resolved | 3 |
| partially-resolved | 1 |
| superseded | 1 |
| needs-human | 1 |

---

## Entries

### OQ-040

**Question:** Sacrifice-zone formula has known limitations

**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 184)
**Section:** Active weaknesses
**Date:** unknown

**Proposed Status:** partially-resolved
**Confidence:** high
**Resolving documents:**
- `experiments/RESONANT_RETURN.md`: confirms N-scaling limitation (360x at N=5, 139x at N=9) and notes formal proof of optimality is open
- `experiments/COCKPIT_UNIVERSALITY.md` (Section 3.5-3.6): documents hardware gap (3.2x vs 360x theory) and theta sensitivity
- `experiments/COCKPIT_SCALING.md`: scaling analysis to N=11
**Rationale:** This is a section heading captured as an entry, not a single question. The four limitations it introduces (N-scaling decay, no optimality proof, hardware gap, theta sensitivity) are all individually documented and characterized across the repo, but two remain unresolved: formal optimality proof and theta re-optimization. Classified as partially-resolved because the limitations are understood but not all overcome.
**Search terms used:** "sacrifice-zone formula", "known limitations", "optimality", "improvement decreases"

---

### OQ-044

**Question:** θ is 1.68x more sensitive than CΨ as objective function (April 2026 cockpit finding), but the formula was not optimized for θ

**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 191)
**Section:** Active weaknesses
**Date:** unknown

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none (the finding is documented but the optimization gap is not addressed)
**Rationale:** The 1.68x sensitivity is empirically confirmed in `experiments/COCKPIT_UNIVERSALITY.md` (Section 3.5, lines 228-230), and the sacrifice-zone formula's optimization target (Sum-MI, not theta) is explicit in `experiments/RESONANT_RETURN.md` (Tests 7-8). However, no document in the repo attempts to re-optimize the formula with theta as the objective function. The weakness is documented, not resolved.
**Search terms used:** "theta", "1.68", "objective function", "optimized for theta", "Sum-MI", "cockpit"

---

### OQ-048

**Question:** θ is 1.68x more sensitive than CΨ as objective function (April 2026 cockpit finding), but the formula was not optimized for θ

**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 191)
**Section:** 6. Sacrifice-zone formula has known limitations
**Date:** unknown

**Proposed Status:** needs-human
**Confidence:** high
**Resolving documents:** n/a
**Rationale:** This is a duplicate of OQ-044: identical verbatim text from the same source file and line number (191). The inventory collected it twice because the extraction matched two overlapping section headings ("Active weaknesses" and "6. Sacrifice-zone formula has known limitations") that both span this line. Flagging for human decision on whether to merge with OQ-044 or remove.
**Search terms used:** n/a (duplicate identification)

---

### OQ-100

**Question:** **Bures curvature is noisy.** The curvature formula involves second derivatives of sparse data. Reliable curvature estimation requires 50+ densely spaced time points.

**Source:** `experiments/COCKPIT_UNIVERSALITY.md` (line 325)
**Section:** 5. Limitations and caveats
**Date:** April 2, 2026

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** The Bures curvature is computed and shown to be finite (K approximately -25 at CΨ approximately 1/4) in `experiments/INFORMATION_GEOMETRY.md`, but the specific data-sparsity limitation remains unaddressed. No follow-up document increases the time-point density or proposes a noise-robust curvature estimator.
**Search terms used:** "Bures curvature", "curvature", "second derivative", "sparse data", "time points", "INFORMATION_GEOMETRY"

---

### OQ-109

**Question:** Active Weakness #4 (the natural variable u): REFORMULATED

**Source:** `experiments/CRITICAL_SLOWING_AT_THE_CUSP.md` (line 366)
**Section:** 8. Consequences for Open Questions
**Date:** unknown

**Proposed Status:** superseded
**Confidence:** high
**Resolving documents:**
- `experiments/CRITICAL_SLOWING_AT_THE_CUSP.md` (lines 365-372): provides the reformulation
- `docs/WEAKNESSES_OPEN_QUESTIONS.md` (lines 275-295, April 2026 update): u(t) approximately 0.61 * Psi^1.02 on Bell+ trajectories; u is an algebraic conjugation variable, not an independent dynamical coordinate
**Rationale:** The original question ("does u(t) provide a simpler trajectory than CΨ?") was answered negatively for real Bell+ trajectories and reformulated to "does u carry independent information on complex trajectories?" The text itself says "REFORMULATED" and the newer formulation lives in WEAKNESSES_OPEN_QUESTIONS.md lines 285-295.
**Search terms used:** "natural variable u", "u(t)", "Active Weakness #4", "reformulated", "conjugation variable"

---

### OQ-122

**Question:** **Closed form for inner positions: topology-dependent (partially resolved).** d_real(2) differs between topologies (Chain=14, Star=16, Complete=36 at N=4). No universal formula exists for k >= 2. The weight-2 kernel vectors transform under mixed S_N representations, not the trivial representation as at k=1. See [Weight-2 Kernel](WEIGHT2_KERNEL.md).

**Source:** `experiments/DEGENERACY_PALINDROME.md` (line 394)
**Section:** Open questions
**Date:** April 3, 2026

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `experiments/WEIGHT2_KERNEL.md`: provides deep analysis showing d_real(2) is topology-dependent, with data for Chain, Star, Complete at N=4
- `docs/proofs/PROOF_PARITY_SELECTION_RULE.md`: proves the structural reason (weight-2 kernel vectors transform under mixed S_N representations)
**Rationale:** Negative resolution: the question "is there a universal closed form for k >= 2?" is answered definitively no. The entry itself already contains the resolution ("No universal formula exists for k >= 2") and cites the structural reason (representation theory). The palindrome property d_real(k) = d_real(N-k) is universal, but individual values require graph structure.
**Search terms used:** "d_real", "weight-2", "inner positions", "closed form", "topology-dependent", "WEIGHT2_KERNEL"

---

### OQ-139

**Question:** What happens with complex numbers in the formula?

**Source:** `experiments/MATHEMATICAL_FINDINGS.md` (line 305)
**Section:** Open Questions
**Date:** 2026-01-30

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** The R = CΨ^2 formula has been explored exclusively with real-valued inputs throughout the repo. No document discusses complex-valued extensions, imaginary components, or analytic continuation. The mathematical analysis in MATHEMATICAL_FINDINGS.md (Sections 1-7) covers limits, derivatives, integrals, and optimization, all for real numbers.
**Search terms used:** "complex numbers", "complex", "imaginary", "analytic continuation", "MATHEMATICAL_FINDINGS"

---

### OQ-209

**Question:** 5 Threshold formula: ANSWERED

**Source:** `experiments/STAR_TOPOLOGY_OBSERVERS.md` (line 683)
**Section:** 8. Open Questions (partially answered 2026-03-07)
**Date:** March 4, 2026 (updated March 8, 2026)

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `experiments/STAR_TOPOLOGY_OBSERVERS.md` (lines 683-708): contains the answer in-place with verification table and fitted formulas
- `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (line 439): cross-references the threshold formula in speed-fidelity tradeoff analysis
**Rationale:** The entry itself says "ANSWERED" and contains the resolution: J_th(gamma) approximately 7.35 * gamma^1.08 + 1.18 (R^2 = 0.999), verified at N=2 with gamma in [0.001, 0.2]. The answer was produced within the same document and cross-referenced in downstream work.
**Search terms used:** "threshold formula", "J_th", "star topology", "ANSWERED"

---

### OQ-264

**Question:** **Implication**: The transistor architecture works best for **small, localized** subsystems (2-4 qubits per side). For larger systems, longer uniform chains with the sacrifice-zone formula outperform hierarchical topologies (Section 7.3, falsified).

**Source:** `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (line 435)
**Section:** 8. Limitations and Failure Modes
**Date:** 2026-03-21

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (lines 374-400): Section 7.3 marked [FALSIFIED], hierarchical architecture produces identical MI to uniform chain at N=3-11
- `experiments/SCALING_CURVE.md` (lines 85-99): empirical verification at N=3, 5, 7, 9, 11 showing identical MI
**Rationale:** The entry is itself the documented conclusion of a falsified hypothesis. Section 7.3 is explicitly marked [FALSIFIED] and the empirical evidence is in SCALING_CURVE.md: "The hierarchy exists in the naming convention, not in the physics." The limitation is fully characterized; no open question remains.
**Search terms used:** "transistor architecture", "hierarchical", "uniform chain", "falsified", "Section 7.3", "SCALING_CURVE"
