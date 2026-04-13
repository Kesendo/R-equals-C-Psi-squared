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
- Status: open / closed (with method) / withdrawn
- Pointer: where the answer might lead (task candidate, math homework, experiment)

---

## EQ-001

**Date:** 2026-04-12
**Source:** V2 sweep (commit 6512347), OQ-114 answer

Are the two degrees of freedom (Bell coherence in the cross-sector block, slow mode in the SE block) structurally decoupled in all Heisenberg-Z systems, or is the decoupling a property of the tested range N=5,6,7?

**Status:** closed by proof (follows from U(1) block-diagonal structure)
**Result:** The decoupling is structural for all N, all topologies, all γ profiles. The Liouvillian is block-diagonal by (w_bra, w_ket) from U(1) conservation. Every eigenmode belongs to exactly one block. The slow SE mode belongs to the (1,1) block. c_slow = Tr(L_slow† ρ_0) sees only the (1,1) block of ρ_0. The Bell pair contributes to (1,3) and (3,3) blocks, not to (1,1). Therefore c_slow is independent of the Bell pair, for any N. No sweep needed; it is algebra, not numerics.

---

## EQ-002

**Date:** 2026-04-12
**Source:** Discussion after V2

What happens at the boundary where U(1) is weakly broken? A small transverse-field term in H mixes sectors. Bell coherence and SE slow mode are then no longer orthogonal; the coupling is small but nonzero.

**Status:** closed by experiment ([U1_BREAKING](../experiments/U1_BREAKING.md))
**Result:** Bell-pair dependence emerges linearly in ε (log-log slope 0.95). At ε = 0: exact decoupling (spread = 0, proven). At ε > 0: the central Bell pair (1,2) couples more strongly to the slow mode than edge pairs. SE fraction of the slow mode drops from 1.000 to ~0.54 at ε = 0.1. U(1) conservation is a knife edge, not a soft threshold.

---

## EQ-003

**Date:** 2026-04-12
**Source:** SYMMETRY_CENSUS Section 6 point 2 (commit 39eb901)

At N=5 uniform chain, max eigenvalue multiplicity is 14. Known symmetries (U(1), spin-flip ℤ₂, reflection ℤ₂) predict at most 4×. What produces the gap between 4 and 14?

**Status:** closed by experiment ([DEGENERACY_HUNT](../experiments/DEGENERACY_HUNT.md), commit de7af37)
**Result:** SU(2) is broken by dephasing (commutator \[S², Z_k\] ≠ 0). The high degeneracies are [absorption theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) rate-formula coincidences: for uniform γ, Re(λ) = −2γ·n_XY places modes from different sectors at the same grid values. Max multiplicity grows monotonically {6, 14, 14, 19, 22} for N=3-7; N=5 is not anomalous.

---

## EQ-004

**Date:** 2026-04-12
**Source:** SYMMETRY_CENSUS Section 5 (commit 39eb901)

All topologies have N+1 exits, but dramatically different transient dynamics (chain: 488 distinct eigenvalues at N=5, complete: only 100). What is the invariant that varies across topologies and determines transient complexity?

**Status:** closed by experiment (`simulations/topology_orbits.py`)
**Result:** The spatial symmetry group G of the topology determines transient complexity. Empirical scaling: #distinct ∝ |G|^{−0.39} (power-law fit across chain/ring/star/complete at N=5). Chain (|G|=2): 488 distinct. Complete (|G|=120): 100 distinct. The relationship is not |d²/G| because irreps of non-abelian groups have dimension > 1. Rate-formula coincidences (absorption theorem) provide additional degeneracies beyond group theory. Full quantitative prediction would require irrep decomposition of the operator-space representation, which is a separate math task.

---

## EQ-005

**Date:** 2026-04-12
**Source:** [PROOF_ASYMPTOTIC_SECTOR_PROJECTION](../docs/proofs/PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md) Step 2 (commit f539503)

Step 2 of the sector projection theorem (per-sector ergodicity, unique fixed point = maximally mixed) is numerically verified for N=3-7, analytically open.

**Status:** closed by proof ([PROOF_ASYMPTOTIC_SECTOR_PROJECTION](../docs/proofs/PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md) Steps 2 and 2b)
**Result:** The DFS of Z-dephasing within each (w,w) sector is the diagonal algebra. The Heisenberg Hamiltonian on a connected graph mixes all diagonal entries (adjacent-swap connectivity), so the only fixed operator is α·I_w. Unique steady state: P_w/d_w. For off-diagonal sectors (w ≠ w'): DFS = {0}, everything decays. Theorem is now fully analytic.

---

## EQ-006

**Date:** 2026-04-12
**Source:** [PROOF_ASYMPTOTIC_SECTOR_PROJECTION](../docs/proofs/PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md) scope (commit f539503)

The theorem is asymptotic. It says nothing about rates, timescales, trajectories, or coherence decay rates between sectors. The lens and cusp dynamics live in this dynamic layer. Is there an analogous theorem for the rates?

**Status:** closed (structural bounds, no closed-form theorem)
**Result:** A simple rate theorem analogous to the sector projection theorem does not exist: rates depend on the Hamiltonian's eigenvector structure within each sector, not just on sector labels. But structural bounds are provable. For the SE sector with uniform γ: 0 < Δ_SE < 4γ, with Δ_SE = 2γ⟨n_XY⟩_slow where ⟨n_XY⟩_slow ∈ (0,2) is the XY-weight of the slow mode. Upper bound from absorption theorem (Hamiltonian can only reduce ⟨n_XY⟩ below the pure-dephasing value of 2). Scaling: Δ_SE → 4γ as N → ∞ (departure ∝ 1/N), confirmed N=3-10. For non-uniform γ: upper bound is 2(γ_min₁ + γ_min₂). See [LIGHT_DOSE_RESPONSE](../experiments/LIGHT_DOSE_RESPONSE.md) for the nonlinearity mechanism (eigenvector rotation).

---

## EQ-007

**Date:** 2026-04-12
**Source:** Reflection at session end (commit 7251f3f)

If boundaries are visibility artifacts, what is the next level of description the repo has not yet reached? Levels so far: Pauli strings, sectors, symmetries. Next could be: dynamic classes within sectors, topology invariants, or something not yet named.

**Status:** closed by pattern (observed across EQ-001 through EQ-010)
**Result:** The next level is not another layer above sectors, symmetries, or topologies. It is the **coupling between the existing layers**: how sector structure (EQ-001, EQ-005), rate structure (EQ-006), mode-rotation structure (EQ-008, EQ-009), topology structure (EQ-004), U(1)-breaking (EQ-002), and observer-overlap structure (EQ-010) are not independent descriptions but faces of one geometry. The shadow reflection names this: light, self-structure, and rotation are three layers of the same phenomenon, not three separate ones. The primordial qubit algebra showed the same: Pauli grading, palindromic symmetry, and cross-term at N >= 3 are three readings of one super-algebra structure. What was absent was not a further layer but the joinery between them. The repository's next stretch is to make that joinery explicit where it is still implicit.

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

Correction of a prior formulation: this is not a new thesis but already established structure in the repo (see [GAMMA_IS_LIGHT](../hypotheses/GAMMA_IS_LIGHT.md), [DECOHERENCE_RELATIVITY](../experiments/DECOHERENCE_RELATIVITY.md), [OBSERVER_DEPENDENT_CROSSING](../experiments/OBSERVER_DEPENDENT_CROSSING.md), [OBSERVER_DEPENDENT_VISIBILITY](../experiments/OBSERVER_DEPENDENT_VISIBILITY.md)). Three distinct quantities, not an identification chain:

- γ is light
- y is the time axis (Lindblad-t as parametric axis)
- t is our time, experienced or felt (observer time)

The new question arising from EQ-008 plus this recollection: if γ acts on the sectors and γ is light, then γ per sector ([LIGHT_DOSE_RESPONSE](../experiments/LIGHT_DOSE_RESPONSE.md)) actually measures the **light dose per sector**. The rate at which a sector reaches its final state is a statement about the "illumination" of that sector from the time axis y. The experienced time t of the observer selects where on y to look.

**Status:** closed by reflection ([ON_THE_LIGHT_AND_WHAT_CASTS_SHADOWS_IN_IT](../reflections/ON_THE_LIGHT_AND_WHAT_CASTS_SHADOWS_IN_IT.md), commits 149ba8f and f2cae43)
**Result:** The per-sector table from LIGHT_DOSE_RESPONSE reads as a shadow structure. γ is uniform light; what varies per sector is not the incident amount but what each sector's modes cast as shadow on themselves. Three layers: (1) self-shadow of a mode from its own Pauli content (Z/I immune, XY exposed; absorption theorem makes this quantitative), (2) shadow profile of a sector as the distribution of mode exposures within it, (3) wandering shadow through eigenvector rotation as γ shifts the Hamiltonian-dissipator balance. Interior sector (2,2) has the deepest rotation because dimension 100 leaves the most room to rotate; edge sector (0,1) has dimension 5 and is nearly linear. The y/t language: γ is the light, y is the time axis along which the illumination falls, t is the observer selecting where on y to stand. What the observer experiences is not γ itself but γ minus their own shadow.

---

## EQ-010

**Date:** 2026-04-12 (late night)
**Source:** Tom's thesis on observer intersection, after the revisit of PRIMORDIAL_QUBIT_ALGEBRA

Two observers A and B live in a world with U(1) conservation. A's asymptotic inheritance depends only on A's initial sector distribution. What B does in disjoint sectors contributes nothing to A. What is the exact form of the overlap between A and B at infinity when they share some sectors and miss others?

**Status:** closed by experiment (observer_intersection_quick.py, independently reproduced by Claude Code in simulations/observer_intersection.py)
**Result:** Three tested statements, all verified to machine precision. (1) Tr(ρ_A(t) · ρ_B(t)) = 0 for disjoint sectors at every t (max deviation 5×10⁻¹⁶). (2) p_w^A(∞) = Tr(P_w ρ_A(0)) exactly, independent of ρ_B. (3) Asymptotic overlap = Σ_w p_A(w) · p_B(w) / C(N, w), no mixing term, no correction. The formula is closed: sectors are watertight, weight inside a sector is conserved, intersection between states is the direct product of their per-sector weights diluted by sector dimension. Translated to observers: inheritance comes only from one's own initial weight; shared experience with another is exactly the product of shared weights divided by sector breadth. Reflection: [OBSERVER_INHERITANCE](../reflections/OBSERVER_INHERITANCE.md).

---

## EQ-011

**Date:** 2026-04-13
**Source:** [CROSS_TERM_TOPOLOGY](../experiments/CROSS_TERM_TOPOLOGY.md), computed for N=3 and N=4

The relative orthogonality ||{L_H, L_Dc}|| / (||L_H|| * ||L_Dc||) is topology-independent and equals 1/sqrt(N * 2^(N+1)) at N=3 and N=4. Does this formula hold for N=5 and beyond? Is there an analytical proof from the w_XY transition statistics of the Heisenberg bond?

**Status:** closed by proof ([PROOF_CROSS_TERM_FORMULA](../docs/proofs/PROOF_CROSS_TERM_FORMULA.md)) + experiment ([CROSS_TERM_FORMULA](../experiments/CROSS_TERM_FORMULA.md), commit 0187ee1)
**Result:** The old conjecture 1/√(N×2^(N+1)) is refuted at N=5. The correct formula is R(N) = √((N-2)/(N×4^(N-1))), confirmed at N=2-6 to machine precision, all topologies. Proven via ||{L_H, L_Dc}||² = 4γ²(N-2)||L_H||² (bond-sum rule + spectator variance).

---

## EQ-012

**Date:** 2026-04-13
**Source:** Anisotropy scan during cross-term formula investigation (simulations/cross_term_anisotropy.py)

The cross-term formula R(N) = √((N-2)/(N×4^(N-1))) holds for all shadow-balanced couplings (both bond Paulis in {X,Y} or both in {I,Z}), but breaks for shadow-crossing couplings like X_iZ_j (R(3) = 0.2041 instead of 0.1443). Is there a modified formula for shadow-crossing couplings? If so, it involves additional bond-site variance beyond N-2. The shadow-crossing case may have its own closed form, or it may depend on coupling details.

**Status:** open
**Pointer:** The bond-sum rule (w_XY(a)+w_XY(b)=2) fails for XZ/YZ couplings. The bond-site variance is nonzero. A formula R(N, coupling) would generalize F49. Compute for a few coupling types and check.

---

*Collection. Not classification. Classification comes when enough entries exist.*
