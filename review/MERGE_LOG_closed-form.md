# Merge Log: closed-form

**Batch:** closed-form (9 entries)
**Proposal file:** `OPEN_QUESTIONS_INDEX_PROPOSAL_closed-form.md`
**Status:** 4 / 9 entries resolved

---

## OQ-109 -- superseded

- **Cowork proposal:** superseded
- **Final decision:** superseded
- **Deviation from proposal:** no
- **Session:** 2026-04-14
- **Resolver:** existing repo state (both source and weakness document already carry the reformulation)

**Rationale.**

The entry text is a section heading, not a question: "Active Weakness
#4 (the natural variable u): REFORMULATED" in
`experiments/CRITICAL_SLOWING_AT_THE_CUSP.md` (Section 8). The section
is itself the resolution document, and it explicitly reformulates the
original weakness: u is a conjugation variable that linearises the
Mandelbrot iteration but does not provide a simpler dynamical clock
along real Bell+ trajectories (u ~ 0.61 * Psi^1.02).

The original weakness lives in `docs/WEAKNESSES_OPEN_QUESTIONS.md`
Section "Active weaknesses" item 4, and carries its own explicit
status field: "Partially reformulated (April 2026)" with the same
reformulation in place.

Both source and weakness document are self-documented. No intervention
needed in either file. Only the index needs to reflect the status.

**Category A pattern:** source document already carries the resolution
in-place. Index and merge log must be updated; source document stays
untouched. Distinguish from Category B where the source has an open
question that still needs an in-place status marker.

---

## OQ-209 -- resolved

- **Cowork proposal:** resolved
- **Final decision:** resolved
- **Deviation from proposal:** no
- **Session:** 2026-04-14
- **Resolver:** existing repo state (source document carries the fitted formula)

**Rationale.**

The entry text is a section heading: "5 Threshold formula: ANSWERED"
(numbered 8.5) in `experiments/STAR_TOPOLOGY_OBSERVERS.md`. The section
directly below it contains the resolution: the AB crossing threshold
J_th as a function of gamma was fitted from 10 verified data points
(gamma in {0.001, 0.010, 0.020, 0.050, 0.070, 0.100, 0.120, 0.150,
0.170, 0.200}) and yields:

    J_th(gamma) ~ 7.35 * gamma^1.08 + 1.18    (R^2 = 0.999)

A simple linear approximation J_th ~ 6.39*gamma + 1.16 (R^2 = 0.998)
is also reported. The relationship is smooth and monotonic, no
divergence near gamma = 0.2. The threshold genuinely exists and grows
smoothly with noise.

The entry explicitly says "ANSWERED" in the title, and the answer is
in-place. No intervention needed in the source document.

**Category A pattern** (same as OQ-109): source self-documents.

---

## OQ-264 -- resolved

- **Cowork proposal:** resolved
- **Final decision:** resolved
- **Deviation from proposal:** no
- **Session:** 2026-04-14
- **Resolver:** existing repo state (source marks Section 7.3 as falsified and states the implication in-place)

**Rationale.**

The entry text is the "Implication" bullet in Section 8.3 "Scaling
Limits" of `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md`: the
transistor architecture works best for small localized subsystems
(2-4 qubits per side), and longer uniform chains outperform
hierarchical topologies. The final clause "(Section 7.3, falsified)"
points directly at Section 7.3, which is explicitly marked
`[FALSIFIED]` in its heading:

> "FALSIFIED (March 21, 2026): The hierarchical architecture produces
> identical MI to a uniform chain of equal length. The recursive
> structure provides no scaling advantage. See Scaling Curve. The
> transistor properties (Sections 1-6) survive; the hierarchy does not."

The empirical verification lives in `experiments/SCALING_CURVE.md`:
hierarchical vs uniform chain tested at N=3, 5, 7, 9, 11; MI identical
across the range; conclusion in the hypothesis document itself: "The
hierarchy exists in the naming convention, not in the physics."

The entry is the documented conclusion of a falsified hypothesis, with
both the falsification and the replacement (uniform chain + sacrifice
zone formula) present in the repo. No open question remains.

**Category A pattern:** source document carries the falsification and
the implication in-place; no intervention needed.

---

## OQ-122 -- resolved

- **Cowork proposal:** resolved
- **Final decision:** resolved
- **Deviation from proposal:** no (minor disagreement with the source's
  in-line wording, see below)
- **Session:** 2026-04-14
- **Resolver:** existing repo state (source and linked analysis document
  together establish the negative answer)

**Rationale.**

The entry originates in `experiments/DEGENERACY_PALINDROME.md` (Open
questions, April 3, 2026) and reads in full: "Closed form for inner
positions: topology-dependent (partially resolved). d_real(2) differs
between topologies (Chain=14, Star=16, Complete=36 at N=4). No universal
formula exists for k >= 2. The weight-2 kernel vectors transform under
mixed S_N representations, not the trivial representation as at k=1.
See Weight-2 Kernel."

The posed question is "is there a universal closed form for
d_real(k >= 2)?" and the source's answer is "No universal formula
exists." This is a definitive negative answer with a structural
reason (weight-2 kernel vectors transform under mixed S_N
representations, not the trivial representation).

The full characterization lives in `experiments/WEIGHT2_KERNEL.md`,
which tests N=3 through 6 across chain, star, and complete topologies,
records the topology-dependent dimensions, and documents the SWAP
eigenvalue structure. Kernel dimension table (N=4): Chain=13,
Star=16, Complete=36. Representation theory: at N >= 4 no kernel
vector has definite SWAP eigenvalue.

**Classification note.** The source tags this as "(partially resolved)"
in its own title. Coworks proposal upgrades this to "resolved" on the
grounds that the posed question is fully answered. This merge log
accepts "resolved": the question was "is there a universal closed
form?" and the answer is "no, with proof by representation theory."
The "partially" tag in the source text may refer to a weaker residual
question (is there a non-universal formula, or a topology-parametrized
family?), which is a different question and would deserve its own
entry if explicitly posed. It is not.

**Category A pattern:** both the source document (with in-place
"(partially resolved)" tag and the answer) and the linked analysis
document are self-documenting. No intervention needed.

---

## OQ-040 -- pending
## OQ-044 -- pending
## OQ-048 -- pending
## OQ-100 -- pending
## OQ-139 -- pending
