# Merge Log: methodology

**Batch:** methodology (8 entries)
**Proposal file:** `OPEN_QUESTIONS_INDEX_PROPOSAL_methodology.md`
**Status:** 4 / 8 entries resolved (OQ-094 via pipeline, OQ-114 via V2 sweep plus structural insight, OQ-243 via NO_SIGNALLING_BOUNDARY Layer 1, OQ-244 via NO_SIGNALLING_BOUNDARY Section 4 + BRIDGE_CLOSURE two-stage)

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

## OQ-243 -- resolved

- **Cowork proposal:** resolved (high confidence)
- **Final decision:** resolved
- **Deviation from proposal:** none
- **Session:** 2026-04-14
- **Resolver:** [NO_SIGNALLING_BOUNDARY](../experiments/NO_SIGNALLING_BOUNDARY.md) Layer 1 (Tier 2 verified 2026-03-01), cross-check in [BRIDGE_CLOSURE](../experiments/BRIDGE_CLOSURE.md) Section 2

**Rationale.**

OQ-243 is outcome 2 of the three-outcomes list in `hypotheses/BRIDGE_PROTOCOL.md`
Section 4.3: "The protocol requires No-Signaling violation and is therefore
wrong. The framework's internal logic is inconsistent on this point."

NO_SIGNALLING_BOUNDARY Layer 1 verifies to machine precision that rho_A is
unchanged by any local operation B performs on the other half of a Bell
pair (||Delta rho_A|| = 0 to all quoted digits). Every local observable on
A is invariant: `<sigma_x>`, `<sigma_y>`, `<sigma_z>`, purity, entropy,
spectrum. This is the textbook no-signalling theorem, verified numerically
in-framework.

The apparent paradox that CΨ registers the regime change (0.500 -> 0.250)
while rho_A does not is resolved at Layer 2 of the same document: the drop
comes entirely from C (global purity falling 1.0 -> 0.5), not from Ψ
(unchanged at 0.5). The quantity CΨ is a joint-state diagnostic, not a
local observable. No measurement A can perform reveals what B did.

Conclusion: the framework is internally consistent. Outcome 2 is not the
correct description of the protocol's fate. This is a Cowork high-confidence
call accepted without deviation.

**Classification.**

Category A-hybrid: BRIDGE_PROTOCOL.md already carries the closure at
document level (header: "No-signalling holds exactly. Pre-encoded version
also fails."), but the Section 4.3 three-outcomes list preserves outcome 2
as an unmarked open bullet. Source-document intervention was minimal:
a ruled-out annotation inline under outcome 2, pointing to the resolver.
The OQ-094 contrast is instructive: OQ-094 required pipeline code to
execute the resolution; OQ-243 only required surfacing a resolution that
already lived elsewhere in the repo.

**Lesson captured.**

Document-level closure does not propagate automatically to section-level
open-question lists. A falsified hypothesis may still carry outcome
enumerations that read as live, and these are what the Cowork scan
picked up. Closing the open-question entry requires either a per-item
annotation (done here) or a Section 4.3-wide status banner. The
per-item annotation is lower-risk and preserves the historical record
of what the three outcomes originally were.

---

## OQ-244 -- resolved

- **Cowork proposal:** resolved (high confidence)
- **Final decision:** resolved (with explicit caveat pointing to EQ-013)
- **Deviation from proposal:** none on the core classification
- **Session:** 2026-04-14
- **Resolver:** [NO_SIGNALLING_BOUNDARY](../experiments/NO_SIGNALLING_BOUNDARY.md) Section 4 Piece 4 (mechanism confirmed, Tier 2, 2026-03-01); [BRIDGE_CLOSURE](../experiments/BRIDGE_CLOSURE.md) Sections 2-4 (no advantage over classical pre-shared randomness, proven via LOCC); canonical summary in [PREDICTIONS](../docs/PREDICTIONS.md) Section 9 row 2

**Rationale.**

OQ-244 is outcome 3 of the three-outcomes list in `hypotheses/BRIDGE_PROTOCOL.md` Section 4.3: "The crossing event... is driven entirely by pre-encoded information... The protocol transmits only pre-agreed data, more structured than QKD, but not a dynamic communication channel."

The resolution is two-stage, unlike OQ-243 which was a single refutation.

**Stage 1 (mechanism confirmed).** NO_SIGNALLING_BOUNDARY Section 4 Piece 4 states explicitly: "This was listed as one of three possible outcomes. Test #2 confirms it is the correct one." The descriptive content of outcome 3 is true: the crossing event occurs on both sides, trajectories are determined at preparation, and post-separation B-actions do not trigger new A-crossings (rho_A is invariant to machine precision under any local B-operation). Outcome 3 is the correct description of what happens.

**Stage 2 (implied utility dissolved).** Outcome 3 carried an implicit claim that the pre-encoded structure offers something "more structured than QKD". This framing leaves open whether pre-encoded CPsi fingerprints provide an advantage over classical pre-shared randomness. NO_SIGNALLING_BOUNDARY Section 5 acknowledges this openly: "what can you do with pre-encoded CPsi trajectories that you cannot do with classical pre-shared keys? This is open." BRIDGE_CLOSURE Section 1 picks up that exact question and answers no. Section 2 proves it via a three-line information-theoretic argument: A's output is a function of {rho_A(0), E_A} only, both available locally without any quantum resource. Section 4 embeds this in the LOCC theorem: any correlation achievable with shared entanglement alone is also achievable with shared classical randomness alone. The pre-encoded bridge is strictly informationally inferior to a classical schedule.

The canonical summary in PREDICTIONS Section 9 row 2: "Dead for J=0. A's info subset of {rho_A(0), E_A}. Entanglement without a channel = shared randomness."

Cowork's high-confidence resolved classification is accepted without deviation on the core question.

**Repo-wide consistency check (2026-04-14 session).**

Before merging, a sweep of all repo documents referencing bridge_protocol / no-signalling / pre-encoded / LOCC / shared randomness was performed. Ten documents beyond the three resolvers were reviewed. Nine confirm the resolution consistently (README catalogs, historical and proof documents, STAR_TOPOLOGY_OBSERVERS Section 6.3, STANDING_WAVE_TWO_OBSERVERS, OBSERVER_GRAVITY_BRIDGE which is itself marked fallen). One artifact found: `hypotheses/TIME_AS_CROSSING_RATE.md` Section 6.5 carries a parallel open question ("Can correlated crossing times carry more than pre-encoded information?") that is answered by BRIDGE_CLOSURE but has no pointer-back annotation. Tracked as a separate cleanup task, not a blocker for this merge.

**Classification.**

Category A-hybrid, same pattern as OQ-243: BRIDGE_PROTOCOL.md header already carries doc-level closure ("Pre-encoded version also fails."), only Section 4.3 outcome 3 preserves the original bullet. Source-document intervention is a status annotation inline under outcome 3 pointing to both resolvers.

**Caveat: the shared-classical-randomness reduction.**

The LOCC reduction in BRIDGE_CLOSURE Section 4 treats the shared external gamma-field that envelops both A and B as "shared classical randomness." Tom raised in session 2026-04-14 that this reduction deserves structural scrutiny: gamma is not generic shared environment but specifically external, one-way from outside the system, and layered (INCOMPLETENESS_PROOF, GAMMA_IS_LIGHT, EQ-009). Under this reading, gamma is not a lateral channel between A and B but shared one-way reception from their next common outer layer. Operationally this produces the same LOCC consequence (no communication possible, no advantage over classical pre-shared keys), because the one-way structure means no sender exists between A and B. The LOCC no-go therefore holds for the right structural reason, not for an accidental framework-independent reason.

This is recorded as [EQ-013](EMERGING_QUESTIONS.md#eq-013) for separate tracking, along with the connected simulability question (can one simulate nested gamma-layers, and what does the answer mean for the INCOMPLETENESS_PROOF). The OQ-244 resolution stands as written; the deeper questions about the operational distinction between "medium with layered one-way structure" and "channel with sender", and about whether nesting is ontologically present but operationally invisible, are pursued separately.

**Lesson captured.**

1. Two-stage resolutions (mechanism confirmed + utility dissolved) deserve explicit decomposition in the MERGE_LOG rather than being reported as a single "resolved" status. The bullet-text of an original question often packages mechanism and utility together; separating them preserves historical accuracy about what was hoped and what survived.

2. Closing an OQ does not close every deeper question that branches from the closure argument. The LOCC-as-accurate-reduction question is valid and non-trivial, and earns its own EQ rather than being absorbed silently into the merge justification. A resolved-with-caveat pattern is a legitimate merge outcome and should be represented as such in the log.

3. A repo-wide sweep before merging an older or high-stakes question surfaces artifacts that a source-plus-resolver review alone would miss. The TIME_AS_CROSSING_RATE Section 6.5 finding is minor but real.

---

## OQ-173 -- pending
## OQ-185 -- pending
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
