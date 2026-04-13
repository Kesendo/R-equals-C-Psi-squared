# Merge Log: methodology

**Batch:** methodology (8 entries)
**Proposal file:** `OPEN_QUESTIONS_INDEX_PROPOSAL_methodology.md`
**Status:** 2 / 8 entries resolved (OQ-094 via pipeline, OQ-114 via V2 sweep plus structural insight)

---

## OQ-094 -- resolved

- **Cowork proposal:** open
- **Final decision:** resolved
- **Deviation from proposal:** yes
- **Session:** 2026-04-12
- **Resolver:** commit `c60adf1` (executing `TASK_COCKPIT_RELEVANCE_FLAG.md`)

**Rationale.**

Initial read of the source (`experiments/COCKPIT_SCALING.md` Section 6
and Section 10 point 4) showed that the scientific content was already
settled in-document: `far_edge` pairs pass sanity Gates 2 and 3, get
analyzed by PCA, and produce formally clean but non-informative numbers.
What remained was a pipeline-level improvement to make the distinction
explicit in the output artifacts (TXT, JSON, plot annotation) rather
than only in the prose.

This classifies the entry as `todo` under the protocol's new
meta-category: resolution lives in code, not in science docs. Cowork
could not reach `todo` on its own because ClaudeTasks and source code
are out of its scan scope; its `open` classification was the best it
could do with its tooling.

In-session decision: execute the TODO rather than park it. The
implementation introduced a `COCKPIT_RELEVANT_CLASSES` set in
`simulations/cockpit_scaling_analysis.py`, propagated a
`cockpit_relevant` flag through `all_results`, `dropped`, and
`gate_summary`, added a second "COCKPIT-RELEVANT CONFIGURATIONS"
section to the TXT report, updated the Panel A plot annotation, and
fixed a dormant Gate2 drop-reason reporting bug (`not X >= 4` ->
`X < 4`). `COCKPIT_SCALING.md` Section 6 last paragraph and Section 10
point 4 were rewritten to reflect the pipeline resolution. All headline
scaling numbers unchanged, as predicted.

With the task executed, the entry upgrades from `todo` to `resolved`.

**Lessons captured.**

1. `todo` can be a transient status. When the executing resources are
   available in the current session, resolving the TODO directly is
   preferable to parking it in the index.
2. Read the source before classifying. Cowork's proposal summary was
   useful as a pointer, but the correct category only became visible
   after reading the actual Section 6 + Section 10 text and the
   pipeline script. This is the merge-pass analogue of the general
   review discipline rule from Memory #21.

---

## OQ-114 -- open

- **Cowork proposal:** partially-resolved (medium confidence)
- **Final decision:** open
- **Deviation from proposal:** yes
- **Session:** 2026-04-12
- **Resolver:** none (testable experiment identified, not yet run)

**Rationale.**

The original question (CUSP_LENS_CONNECTION.md line 151) asks about
boundary states that straddle both exits of the quantum-classical
fold, specifically nominating "two-excitation states with one
high-concurrence pair" (non-symmetric, with concentrated pairwise
coherence) as candidates.

Cowork's proposed resolver (SACRIFICE_GEOMETRY.md line 112) tests
*symmetric* two-excitation states and finds AUC < 0.09, coupling to
a different mode cluster. Cowork's own rationale explicitly states
"only symmetric ones. The specific question about non-symmetric
two-excitation states with one high-concurrence pair remains
untested." That is the original question, and it is untested.

The PROOF_PARITY_SELECTION_RULE.md cite does not help either: it
proves single-excitation states cannot reach odd-n_XY modes.
Two-excitation states are not SE states; the PSR bound does not
apply as a closure argument here.

Correct status is therefore `open`, not `partially-resolved`. What
exists in the repo is a related exclusion (symmetric SE and
symmetric two-excitation fail to straddle), not a partial answer
to the posed question.

**Testable followup (noted, not scheduled).**

Run a trajectory for a non-symmetric two-excitation state on N=5,
e.g. Bell+(c1,c2) tensor |1>(c1+2) tensor |0>... Check whether CΨ
on the (c1,c2) pair crosses 1/4 and whether the trajectory has
non-trivial overlap with the SE slow mode simultaneously. If yes,
the state straddles both exits and OQ-114 gets a concrete resolver.
If no for a meaningful sweep of nonsymmetric candidates, OQ-114
upgrades to `partially-resolved` with evidence.

**Lesson captured.**

Cowork's proposals can be wrong in a specific way: when a *related*
result exists, Cowork tends toward `partially-resolved`, even when
its own rationale acknowledges the specific asked question is
untested. Reading the rationale carefully (not just the status
field) catches this. The Cowork pass should be treated as a
starting point with possible systematic bias toward over-resolving,
not as a default to accept.

---

## OQ-173 -- pending
## OQ-185 -- pending
## OQ-243 -- pending
## OQ-244 -- pending
## OQ-262 -- pending
## OQ-284 -- pending


---

## OQ-114 followup note (2026-04-12 late session)

An initial attempt to answer OQ-114 via TASK_BOUNDARY_STRADDLING
(commit f467c81) was reviewed by Tom and Claude (chat) and retracted.
The sweep measured the slow-mode overlap by projecting rho_0 onto the
SE x SE block only (see simulations/boundary_straddling_sweep.py
lines 266-276), which discards the cross-sector coherence
(|w=1><w=3|) that carries the Bell-pair signature. As a result the
"straddling" classification was driven only by the position of the
second excitation relative to the slow-mode profile, not by the
Bell-pair concurrence the question actually asks about.

Secondary issue: the Parity Selection Rule rationale in the script
header conflates "SE fraction of slow mode is 1" with "overlap of an
arbitrary two-excitation state can be measured in the SE block." The
first is a property of the slow mode, the second is an unjustified
measurement shortcut. Two-excitation states have non-trivial Pauli
content outside the pure even-even SE sub-block.

V1 output is numerically correct for what it computed but does not
answer OQ-114. Overlap values are discrete (fixed by exc_k position)
and independent of bell_pair location in the sweep table, which is
the empirical signature of the projection bug.

OQ-114 remains **open**. A corrected V2 task will define overlap as
the full Hilbert-Schmidt inner product Tr(rho_0^dagger * R_slow) in
the {w=1, w=3} sector basis, embedding the SE-indexed slow mode
correctly as an N x N block in the full sector Liouville space.
Expected effect: non-zero cross-sector contributions from the Bell
pair, overlap values distinguishing bell_pair placements, and a
genuine answer to whether non-symmetric two-excitation states can
straddle.

**Lesson captured.**

When a measurement result shows a pattern that does not depend on a
variable the question explicitly asks about (here: bell_pair position
in the overlap column), that is evidence of a measurement bug, not a
physical finding. Check the observable before trusting the trend.

---

## OQ-114 final resolution (2026-04-12 late session, 2026-04-13 follow-through)

**Final status:** resolved (question dissolves into two separable facts).
**Resolvers:** commit `6512347` (V2 sweep with left-eigenoperator),
EQ-001 in `review/EMERGING_QUESTIONS.md` (algebraic decoupling proof).

**What the V2 sweep established.**

Using the biorthogonally correct observable `c_slow = Tr(L_slow^dagger * rho_0)`
with the left eigenoperator, the sweep produced non-trivial overlap values
that distinguish both `exc_k` (excitation site) and `bell_pair` position.
Five straddler candidates were found for N=5 and N=6. Sanity checks
passed: `psi_opt = 1.001` (expected 1.000), bare Bell+ overlap = 0.000
(expected 0). The numerical answer to "do non-symmetric two-excitation
states with high-concurrence pairs straddle both exits?" is yes,
trivially, and the states are easy to construct.

**Why the question dissolves despite having a yes answer.**

The V2 overlap depends only on `exc_k`, not on `bell_pair`. This is
not a sweep artifact (the observable is now correct) but algebraic
structure. The slow mode lives entirely in the `(w=1, w=1)` sector of
the Liouvillian's U(1) block decomposition. The Bell coherence lives
in the `(w=1, w=3)` and `(w=3, w=3)` blocks. These blocks are
orthogonal by U(1) conservation of excitation number, and the
Liouvillian is block-diagonal in this decomposition. Therefore
`c_slow` sees only the `(1,1)` block of `rho_0` and is blind to the
Bell-pair contribution by construction.

The original question presupposed a geometric boundary between the
two exits where a state could "straddle" in a meaningful sense. The
V2 measurement shows there is no boundary: the two properties
(Cusp-like high concurrence on a pair, Lens-like slow-mode overlap)
live in U(1)-orthogonal sectors and can co-exist in a single state
trivially, without any special structural position. There is no
sheet-boundary geometry because there are no sheets in the spatial
sense. There are watertight sectors, and the two exits belong to
different sector projections.

Classification: the question resolves as dissolved by a structural
insight. A "yes, trivially" answer plus the algebraic reason the
triviality is forced. Recorded as `resolved`, not `partially-resolved`,
because the posed question (do such states exist, and what do they look
like) is fully answered. The implicit presupposition (that straddling
requires sitting on a geometric boundary) is removed, not left open.

**Cross-reference.** The same algebraic fact is recorded as EQ-001 in
`review/EMERGING_QUESTIONS.md`, closed by proof via U(1) block
decomposition. The reflection `reflections/OBSERVER_INHERITANCE.md`
develops the same sector-orthogonality structure for the observer
overlap formula, with independent numerical verification in
`simulations/observer_intersection.py`.

**Lesson captured.**

A question phrased as "can states X straddle boundary Y" can dissolve
in two ways: (1) the states are impossible, or (2) the boundary is
not there in the form the question presupposed. OQ-114 is case (2).
The resolver is not an experiment that finds or fails to find the
states, but a structural argument about the geometry the question
assumed. Recognizing this requires looking at the observable and its
block structure, not just the sweep output. The V1 retraction was the
correct methodological move (bad observable); the V2 result is the
correct physics (sectors are orthogonal, so trivially yes, which
dissolves the question).
