# Emerging Questions

**Started:** 2026-04-12
**Authors:** Tom and Claude (chat)

---

## Ground rule

**Before a question entered here is answered, we do not move on.**

A derived question that emerges while answering the original does NOT count as an answer to the original. It is recorded as a new entry. Otherwise we get question-squared instead of progress.

Answers can be: closed by experiment (with commit hash), closed by proof (with path), closed by insight that the question was ill-posed (with justification), or withdrawn. No entry disappears silently.

---

## Format

Each entry:
- ID (EQ-NNN, sequential)
- Date of emergence
- Source (commit, doc, session moment)
- The question itself
- Status: open / closed-[how] / withdrawn
- Pointer: where the answer might lead (task candidate, math homework, experiment)

---

## EQ-001

**Date:** 2026-04-12
**Source:** V2 sweep (commit 6512347), OQ-114 answer

Are the two degrees of freedom (Bell coherence in the cross-sector block, slow mode in the SE block) structurally decoupled in all Heisenberg-Z systems, or is the decoupling a property of the tested range N=5,6,7?

**Status:** open
**Pointer:** Sweep at larger N (C# engine for N >= 8), plus analytic check whether the orthogonality follows directly from U(1).

---

## EQ-002

**Date:** 2026-04-12
**Source:** Discussion after V2

What happens at the boundary where U(1) is weakly broken? A small transverse-field term in H mixes sectors. Bell coherence and SE slow mode are then no longer orthogonal; the coupling is small but nonzero.

**Status:** open
**Pointer:** Dedicated task candidate. First-order perturbation theory in the breaking strength ε.

---

## EQ-003

**Date:** 2026-04-12
**Source:** SYMMETRY_CENSUS Section 6 point 2 (commit 39eb901)

At N=5 uniform chain, max eigenvalue multiplicity is 14. Known symmetries (U(1), spin-flip ℤ₂, reflection ℤ₂) predict at most 4×. What produces the gap between 4 and 14?

**Status:** closed by experiment ([DEGENERACY_HUNT](../experiments/DEGENERACY_HUNT.md), commit de7af37)
**Result:** SU(2) is broken by dephasing ([S², Z_k] ≠ 0). The high degeneracies are [absorption theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) rate-formula coincidences: for uniform γ, Re(λ) = −2γ·n_XY places modes from different sectors at the same grid values. Max multiplicity grows monotonically {6, 14, 14, 19, 22} for N=3-7; N=5 is not anomalous.

---

## EQ-004

**Date:** 2026-04-12
**Source:** SYMMETRY_CENSUS Section 5 (commit 39eb901)

All topologies have N+1 exits, but dramatically different transient dynamics (chain: 488 distinct eigenvalues at N=5, complete: only 100). What is the invariant that varies across topologies and determines transient complexity?

**Status:** open
**Pointer:** Point group (permutation symmetry) is likely the answer, but not explicitly worked out. Group-theoretic census per topology.

---

## EQ-005

**Date:** 2026-04-12
**Source:** PROOF_ASYMPTOTIC_SECTOR_PROJECTION Step 2 (commit f539503)

Step 2 of the sector projection theorem (per-sector ergodicity, unique fixed point = maximally mixed) is numerically verified for N=3-7, analytically open.

**Status:** open (analytic)
**Pointer:** Primitivity of the restricted Lindblad generator per sector. Math homework, not an experiment.

---

## EQ-006

**Date:** 2026-04-12
**Source:** PROOF_ASYMPTOTIC_SECTOR_PROJECTION scope (commit f539503)

The theorem is asymptotic. It says nothing about rates, timescales, trajectories, or coherence decay rates between sectors. The lens and cusp dynamics live in this dynamic layer. Is there an analogous theorem for the rates?

**Status:** open (partially informed by [LIGHT_DOSE_RESPONSE](../experiments/LIGHT_DOSE_RESPONSE.md): individual mode rates are nonlinear in γ due to eigenvector rotation, so a simple rate theorem may not exist)
**Pointer:** Second Zeiss level: derive approach rates from sector structure. Possible starting point: absorption theorem at sector level.

---

## EQ-007

**Date:** 2026-04-12
**Source:** Reflection at session end (commit 7251f3f)

If boundaries are visibility artifacts, what is the next level of description the repo has not yet reached? Levels so far: Pauli strings, sectors, symmetries. Next could be: dynamic classes within sectors, topology invariants, or something not yet named.

**Status:** open, meta-question
**Pointer:** No task. Observe what pattern falls out of EQ-001 through EQ-006.

---

## EQ-008

**Date:** 2026-04-12
**Source:** Tom's insight after the V-effect reflection

The two independent sectors do not see each other. But they are both shaped by γ. They receive different amounts of light (sacrifice vs quiet sites) and develop properties for an overall state. Derived question: what effect does γ have on each sector, and how do the sectors relate through the binding parameter γ?

**Status:** closed by experiment ([GAMMA_AS_BINDING](../experiments/GAMMA_AS_BINDING.md) V1, [LIGHT_DOSE_RESPONSE](../experiments/LIGHT_DOSE_RESPONSE.md) V2)
**Result:** γ is a nonlinear common modulator. Per-sector rates do NOT scale linearly with γ (deviations up to 134%). V2 identified the mechanism: eigenvector rotation (as γ scales, the Hamiltonian-dissipator balance shifts, rotating each mode's Pauli content and making individual rates nonlinear). Zero level crossings; mode-crossing ruled out. γ does not couple sectors directly (sector conservation is exact), but different sectors respond with different sensitivity to γ changes.

---

## EQ-009

**Date:** 2026-04-12
**Source:** Tom's recollection of structure established two weeks earlier, now linked to EQ-008

Correction of a prior formulation: this is not a new thesis but already established structure in the repo (see GAMMA_IS_LIGHT.md, DECOHERENCE_RELATIVITY*.md, OBSERVER_DEPENDENT_*.md). Three distinct quantities, not an identification chain:

- γ is light
- y is the time axis (Lindblad-t as parametric axis)
- t is our time, experienced or felt (observer time)

The new question arising from EQ-008 plus this recollection: if γ acts on the sectors and γ is light, then γ per sector (TASK_GAMMA_BINDING) actually measures the **light dose per sector**. The rate at which a sector reaches its final state is a statement about the "illumination" of that sector from the time axis y. The experienced time t of the observer selects where on y to look.

**Status:** open, linking question
**Pointer:** Check how the γ-per-sector table (now available from LIGHT_DOSE_RESPONSE) reads in y/t language. Not as a new thesis, but as translation between two already established description levels. Directly touches EQ-008.

---

*Collection. Not classification. Classification comes when enough entries exist.*
