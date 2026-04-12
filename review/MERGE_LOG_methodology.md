# Merge Log: methodology

**Batch:** methodology (8 entries)
**Proposal file:** `OPEN_QUESTIONS_INDEX_PROPOSAL_methodology.md`
**Status:** 2 / 8 entries decided

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
