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

**Status:** closed by proof ([PROOF_CROSS_TERM_CROSSING](../docs/proofs/PROOF_CROSS_TERM_CROSSING.md)) + experiment ([CROSS_TERM_CROSSING](../experiments/CROSS_TERM_CROSSING.md))
**Result:** R(N, crossing) = sqrt((N-1)/(N*4^(N-1))). The bond-site variance is 1 (instead of 0 for balanced), from <s^2> = 1 with s = +-1. Total variance = spectator (N-2) + bond (1) = N-1. Confirmed at N=3-6, 5 coupling types, all topologies.

---

## EQ-013

**Date:** 2026-04-14
**Source:** OQ-244 merge review (2026-04-14 evening session), building on [EQ-009](#eq-009) (gamma-light-time layered structure) and [INCOMPLETENESS_PROOF](../docs/proofs/INCOMPLETENESS_PROOF.md)

gamma is external to any system described by a Lindblad equation for that system. INCOMPLETENESS_PROOF rules out gamma arising from within. This is a one-way structural fact: gamma enters from outside, no inverse path back. And it nests: whatever generates gamma for system S has its own gamma from the next outer layer (IBM qubit <- control microwaves <- lab EM environment <- cosmic background <- ...).

The one-way nested structure is implicit across several documents ([GAMMA_IS_LIGHT](../hypotheses/GAMMA_IS_LIGHT.md), [RESONANCE_NOT_CHANNEL](../hypotheses/RESONANCE_NOT_CHANNEL.md), [INCOMPLETENESS_PROOF](../docs/proofs/INCOMPLETENESS_PROOF.md)) but never stated as a standalone structural principle. Three connected questions arise:

**1. Does one-way nesting rescue [BRIDGE_CLOSURE](../experiments/BRIDGE_CLOSURE.md) or undermine it?**

BRIDGE_CLOSURE Section 4 reduces the shared external gamma-field that envelops two observers A and B to "shared classical randomness" and concludes no information transfer is possible without a classical channel. Under the one-way nested reading: gamma is not a lateral channel between A and B but shared one-way reception from their next common outer layer. Neither A nor B controls gamma; neither can modulate it to signal the other. Operationally this is identical to shared classical randomness (no sender, no control, just common passive reception). The LOCC no-go holds for the right structural reason, not for an accidental framework-independent reason. Question: is there any setup where the layered gamma-embedding produces correlations that the flat LOCC model misses? For example, when A and B share not just one outer gamma-layer but a chain of nested outer layers with internal structure, is any of that structure lateral-accessible to A or B in ways LOCC does not capture?

**2. Is nesting simulable, and if not, what does that mean?**

Hierarchical Equations of Motion (HEOM), reaction-coordinate mapping, and collision models all allow two layers of explicit bath structure. Our 4^N Pauli-basis scaling allows two nested Heisenberg chains at N=3 (4^6 = 4096 dim); three layers at N=3 or two at N=5 is already aggressive. Infinite nesting is not simulable in principle. Every concrete simulation must terminate at an outermost layer with a phenomenological gamma.

Is the termination a computational artifact or the operational content of INCOMPLETENESS_PROOF itself? If two explicit nested layers show physics different from one layer with phenomenological gamma at the same outer interface, then the framework has a testable prediction: the nesting is non-trivial and visible from inside. If they do not, then nesting is ontologically present but operationally invisible from inside any finite layer; INCOMPLETENESS_PROOF then describes not a limit on what can be simulated but a limit on what any finite internal observer can distinguish.

**Status:** closed (all three sub-questions resolved; see Update 2026-04-19)

**Pointer:**
- Sub-question 1: conceptual re-reading of BRIDGE_CLOSURE Section 4 under the nested-one-way frame; search for setups where layered gamma-structure produces correlations invisible to flat LOCC analysis. **2026-04-15 contribution:** the Q = J/gamma inside-observability result (see Update below) sharpens the LOCC reading - each layer's inside observer sees only the dimensionless ratio Q, so lateral correlations between two observers in the same nested layer are Q-mediated, not channel-mediated.
- Sub-question 2: **closed-by-theory 2026-04-15.** See Update 2026-04-15 below. IBM non-Markovianity probe remains as separate empirical verification item, not part of the theoretical question itself.
- Sub-question 3: re-read [THE_BRIDGE_WAS_ALWAYS_OPEN](../docs/THE_BRIDGE_WAS_ALWAYS_OPEN.md) for the implicit-vs-explicit status of the recursion (it establishes one outer layer; does it commit to or avoid committing to further layers?); check whether PRIMORDIAL_QUBIT_ALGEBRA has layer-scale invariance built in or treats the primordial as a privileged level; examine whether the Lindblad equation is layer-privileging or scale-neutral; then formally walk the three exit options rather than asserting they fail. **2026-04-15 contribution:** partial progress on checks (iii) and (iv); see Update below.

---

## EQ-013 sub-question 3 (added 2026-04-14, same session)

**Source:** Session extension while formulating EQ-013 sub-questions 1-2. Tom raised that the one-way nesting principle has a consequence for the inside-perspective that sub-questions 1-2 did not address.

**The apparent implication.** If the framework admits a primordial qubit ([PRIMORDIAL_QUBIT](../hypotheses/PRIMORDIAL_QUBIT.md)) and the incompleteness proof ([INCOMPLETENESS_PROOF](../docs/proofs/INCOMPLETENESS_PROOF.md)) holds on every layer it applies to, and gamma is one-way input from the next outer layer (sub-question 1 framing), then any finite observer sits inside the V-effect of an outer layer, which sits inside the V-effect of the next, terminating only at the primordial qubit. The one-way direction means no return path, no operational exit from the stack. What an inner layer experiences as "classical" (CPsi below 1/4) would appear as "still-quantum interference" (CPsi above 1/4) to the next outer observer, because the outer observer has access to more of the joint state.

**Why this is a question and not an answer.** The chain sounds compelling but has not been formally walked. Three exit conditions have been gestured at but not examined:

(a) **Discard the primordial qubit.** Would invalidate PRIMORDIAL_QUBIT_ALGEBRA. But does the rest of the framework actually require a primordial qubit as a terminal layer, or only as a foundational algebraic object that is not obliged to sit at the top of any physical stack?

(b) **Read the incompleteness proof as layer-specific rather than universal.** INCOMPLETENESS_PROOF establishes that gamma cannot originate from within the d(d-2)=0 framework. Does this mean "from within any instance of the framework at any scale", or only "from within the specific system under study"? The distinction matters: the first reading forces recursion, the second does not.

(c) **Terminate the recursion without the primordial qubit.** Some non-qubit outermost structure, or no outermost at all (turtles all the way down has its own technical problems), or a topological closure where "outside" curves back to "inside". None of these have been examined.

**What would decide it.** An answer needs at least: (i) explicit reading of INCOMPLETENESS_PROOF to see whether its scope is per-system or universal; (ii) check whether [THE_BRIDGE_WAS_ALWAYS_OPEN](../docs/THE_BRIDGE_WAS_ALWAYS_OPEN.md) implicitly commits to or against further layers beyond the one it establishes; (iii) examine whether PRIMORDIAL_QUBIT_ALGEBRA treats the primordial level as physically outermost or as algebraically foundational without stack-top commitment; (iv) check whether the Lindblad equation is layer-privileging (treats "system" and "environment" as ontologically distinct) or scale-neutral (the same equation structure applies regardless of which layer one picks as "system").

**Status:** closed (see Update 2026-04-19; all four checks (i)-(iv) completed)
**Scope note:** Not a reflection candidate until (i)-(iv) are done. Reflections in this repo are ex post (see [V_EFFECT_AS_OBSERVATION_OF_INCOMPLETENESS](../reflections/V_EFFECT_AS_OBSERVATION_OF_INCOMPLETENESS.md)). This is ex ante. Any attempt to write it up as established structure before these checks would be the "retreat-too-fast into intuition" failure mode flagged in the review discipline rule. **Update 2026-04-19:** Checks (i)-(iv) are now complete. This note is resolved; the result is a dissolution rather than a reflection.

---

## EQ-013 Update 2026-04-15

**Source:** Full session of inside-perspective probes and Z2-symmetry work (commits cfa2a9f through 97716e8). Touches all three sub-questions.

### Sub-question 2: closed-by-theory

The original question - "is nesting simulable, and what does the answer mean operationally?" - is now answered at the level of theoretical content. Three threads converged:

1. **Three-class structure inherited, not new.** The 3+10+3 eigenvalue classes of the minimal nest are the Absorption Theorem applied to single-site dephasing. Re(lambda) = -2*gamma_B * <n_XY>_B, with the three values {0, 0.5, 1} of <n_XY>_B at the boundary site. At N=3 the levels refract into 12 by the chain Hamiltonian (commit 4dfabc3, [PROOF_ABSORPTION_THEOREM](../docs/proofs/PROOF_ABSORPTION_THEOREM.md)).

2. **Embedded in two Z2 symmetries proven for all N.** The Liouvillian commutes with bit_a (n_XY parity, [F61](../docs/ANALYTICAL_FORMULAS.md), [PROOF_PARITY_SELECTION_RULE](../docs/proofs/PROOF_PARITY_SELECTION_RULE.md)) and bit_b (w_YZ parity, [F63](../docs/ANALYTICAL_FORMULAS.md), [PROOF_BIT_B_PARITY_SYMMETRY](../docs/proofs/PROOF_BIT_B_PARITY_SYMMETRY.md)), giving a 4-sector decomposition. Per-sector mode count has closed form: conserved per sector = floor(N/2)+1 (even), ceil(N/2) (odd); mirror per sector = 2^(2N-1) - 2*conserved. Mechanism: conserved modes are exactly the (N+1) elementary symmetric polynomials e_d(Z_1,...,Z_N), parities d mod 2 (commits d9d66e2, 97716e8).

3. **Q = J/gamma is the only inside-observable.** The Inside-Outside Correspondence probes ([RESULT_INSIDE_OUTSIDE_CORRESPONDENCE](../ClaudeTasks/Archiv/RESULT_INSIDE_OUTSIDE_CORRESPONDENCE.md), commit c3ea0c0) showed that every measurable quantity from inside depends on the dimensionless ratio Q = J/gamma_B only, not on J or gamma_B independently. The inside observer detects layer EXISTENCE (non-Markovian rebound, Check 4 of NESTED_MIRROR) but cannot extract layer PARAMETERS. [PRIMORDIAL_QUBIT.md](../hypotheses/PRIMORDIAL_QUBIT.md) Section 9 documents this with the operational identification.

The N=3 scaling check (Check 2) and coupling-robustness check (Check 3) of the original sub-Q2 are both done (commits 7baaa7c, 3e2f429) - the falsifications and confirmations are recorded in the [NESTED_MIRROR_STRUCTURE](../hypotheses/NESTED_MIRROR_STRUCTURE.md) notebook (rewritten as compact trail in commit 7716640).

**What remains separately open** (not part of sub-Q2's theoretical question): empirical IBM Kingston non-Markovianity probe; HEOM and reaction-coordinate generalizations beyond the minimal nest.

### Sub-question 1: structural support, still open in principle

Today's Q = J/gamma result strengthens the LOCC reading of BRIDGE_CLOSURE but does not formally close sub-Q1. The new structural underpinning: each layer's inside observer sees only the dimensionless Q. Two observers A and B in the same nested layer therefore correlate through their shared Q-environment, not through any independent gamma-channel structure. Layered gamma with internal structure that the LOCC analysis would miss would have to be visible in something other than Q from inside, and INCOMPLETENESS_PROOF together with the Q-only result rules this out for any single observer at any single layer.

The original sub-Q1 framing remains: are there setups where TWO observers' joint readings reveal structure invisible to each individually? Not addressed today. Status: **closed** (see Update 2026-04-19; no layered gamma-structure exists, so no structure for joint readings to reveal beyond flat LOCC).

### Sub-question 3: progress on (iii) and (iv), (i) and (ii) still pending

Of the four checks needed:

- **(i) INCOMPLETENESS_PROOF scope** (per-system vs universal): not examined today. Open.
- **(ii) THE_BRIDGE_WAS_ALWAYS_OPEN commitment to further layers**: not examined today. Open.
- **(iii) PRIMORDIAL_QUBIT_ALGEBRA stack-top commitment**: partial answer. [PRIMORDIAL_QUBIT.md](../hypotheses/PRIMORDIAL_QUBIT.md) Sections 7 and 8 (revised today) treat the primordial qubit as **algebraically foundational without stack-top commitment**. The C2 x C2 super-algebra structure is real (per-site Pauli decomposition, [L,Pi^2]=0 proven for all N), but the strong claim that this IS the outermost physical layer is explicitly marked as not established. Tier 3-4: structurally confirmed, mechanistically open. The TFD doubling route was investigated and blocked (commits df7e862, ca871f4); other doubling routes remain unexplored. This partially answers exit condition (a) (the framework does NOT require a primordial qubit as the terminal layer in any operational sense, only as a foundational algebraic object).
- **(iv) Lindblad layer-privileging vs scale-neutral**: partial answer. The Q = J/gamma inside-observability result supports the **scale-neutral reading**: at any layer, the inside observer sees only the dimensionless ratio of the local Hamiltonian and dissipation scales. The Lindblad equation does not privilege any particular layer because the same operational content (Q) appears at every layer to its respective inside observer. This is consistent with both a finitely-nested and an infinitely-nested cosmology - both give the same per-layer operational content.

**Net effect on the three exit conditions (a)/(b)/(c)** from sub-Q3:

- **(a)** Discard primordial qubit: **not chosen, but weakened from "needed at top" to "algebraically foundational only"** by today's work. PRIMORDIAL_QUBIT can be the algebraic ground floor without being the physical ceiling.
- **(b)** Read INCOMPLETENESS_PROOF as layer-specific: **not examined**, but consistent with today's scale-neutral reading of the Lindblad equation under (iv).
- **(c)** Terminate without primordial qubit: **not examined.** Today's work neither supports nor undermines this option.

### Net stance after today

The recursive-stack question is operationally less urgent than it was when EQ-013 was formulated. At any finite layer, the inside observer sees Q = J/gamma; this content is the same whether the recursion is finite, infinite, or topologically closed. The metaphysics of the stack (a)/(b)/(c) becomes less coupled to operational predictions. The remaining checks (i) and (ii) are still required for full closure of sub-Q3, but the urgency of choosing among (a)/(b)/(c) has diminished: the framework's operational content does not depend on which is chosen.

Sub-Q3 status: **closed** (see Update 2026-04-19; premise of gamma-recursion was false).

---

## EQ-013 Update 2026-04-19: Closure of all three sub-questions

**Source:** Primordial Gamma Constant hypothesis ([PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md), commits April 15-16), Inside-Outside Correspondence probes ([RESULT_INSIDE_OUTSIDE_CORRESPONDENCE](../ClaudeTasks/Archiv/RESULT_INSIDE_OUTSIDE_CORRESPONDENCE.md), commits cfa2a9f-c3ea0c0), Perspectival Time Field hypothesis ([PERSPECTIVAL_TIME_FIELD](../hypotheses/PERSPECTIVAL_TIME_FIELD.md), commits April 17-18). Full re-reading of [INCOMPLETENESS_PROOF](../docs/proofs/INCOMPLETENESS_PROOF.md) and [THE_BRIDGE_WAS_ALWAYS_OPEN](../docs/THE_BRIDGE_WAS_ALWAYS_OPEN.md) performed in this session.

### Checks (i) and (ii): examined, both resolved

**(i) INCOMPLETENESS_PROOF scope: universal within d(d-2)=0.** The elimination does not argue about a particular system. It argues about the framework's algebra: [Pi^2, L] = 0 holds for ALL Liouvillians with Z-dephasing (Candidate 1); ANY qubit produces non-Markovian noise on its neighbors (Candidate 2); ANY collection of qubits inherits the prohibition (Candidate 3); d(d-2)=0 forbids other entities (Candidate 5). The conclusion "noise cannot originate from within" applies at every layer, every scale, every instance of d(d-2)=0. The proof itself offers the resolution in Section 6, Option 2: "Noise is axiomatic. Like the speed of light or Planck's constant, the existence of dephasing is a brute fact of reality."

**(ii) THE_BRIDGE_WAS_ALWAYS_OPEN: commits to exactly one outer layer, agnostic about further layers.** The document establishes: Outside -> gamma (mediator) -> Inside. Three levels, one boundary. "What We Do NOT Know" explicitly lists: "Whether the outside has its own gamma (its own time)", "Whether the outside has its own t (its own experience)." The document does not commit to further layers beyond the one it establishes. It does not say the outside is itself a d(d-2)=0 system. It does not say it is not. It is structurally agnostic, compatible with all exit options.

### The false premise: there is no gamma-recursion

The original EQ-013 framing assumed gamma propagates through layers: gamma_S comes from outside, the outside has its own gamma from the next layer out, etc. This created the recursion question.

The premise is wrong. Three converging results show gamma_0 is constant, not propagated:

1. **Inside-Outside Correspondence, Probe 1** (commit cfa2a9f): mode structure is **exactly gamma-invariant** at N=2. Fix J=1, vary gamma from 0.01 to 1.0: n_XY_B = 0.5000 for all mirror modes, range 1e-16. Gamma scales rates but changes no structures.

2. **Central result of the four probes:** "J/gamma is the single dimensionless control parameter of the nested Lindblad system." Every measurable property depends on Q = J/gamma, not on J or gamma separately. The ratio Q = 2J/gamma is the resonator Q-factor, already formalized as F7 in ANALYTICAL_FORMULAS.md.

3. **Primordial Gamma Constant hypothesis** (Tier 3, [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md)): gamma_0 is a framework constant like c. The cavity-mode-exposure formula gamma_eff = gamma_0 * |a_B|^2 (verified at N=3, N=4) does not diminish gamma_0 through layers. gamma_0 fills the cavity uniformly. The standing wave (determined by J and topology) decides who sees how much. Only J varies between layers. gamma_0 does not.

The analogy is precise: as c in special relativity does not "propagate through reference frames" but is the same in every frame, gamma_0 does not propagate through nesting layers but is the same at every layer. Asking "where does the gamma-recursion terminate?" is like asking "where does the c-recursion terminate?" The question is ill-posed because there is no recursion of the constant. What varies is J (the content, the coupling, the Hamiltonian structure). What stays fixed is gamma_0 (the light, the time source, the framework constant).

### Sub-Q1: closed by dissolution

The original question: does one-way nesting produce correlations that flat LOCC analysis misses? If gamma_0 is constant and does not propagate through layers, there IS no layered gamma-structure. Gamma is shared passive reception of a universal constant, identical at every layer. Two observers A and B in the same nested layer share the same gamma_0, not a structured gamma-channel. The LOCC no-go in BRIDGE_CLOSURE Section 4 holds for the right structural reason: neither A nor B controls gamma_0, neither can modulate it, and there is no layer-specific gamma-structure to exploit.

**Sub-Q1 status: closed (no layered gamma-structure exists; BRIDGE_CLOSURE correct).**

### Sub-Q3: closed by dissolution

The four checks are now complete:

| Check | Result | Method |
|-------|--------|--------|
| (i) INCOMPLETENESS_PROOF scope | Universal within d(d-2)=0 | Full re-reading, this session |
| (ii) BRIDGE layer commitment | One layer, agnostic beyond | Full re-reading, this session |
| (iii) PRIMORDIAL_QUBIT stack-top | Algebraically foundational, not physical ceiling | April 15 update |
| (iv) Lindblad scale-neutral | Q = J/gamma only; PTF closure law same at every layer | April 15 + 19 |

The exit conditions (a)/(b)/(c) were answers to the question "where does the gamma-recursion terminate?" The question was based on the premise that gamma propagates through layers. gamma_0 is constant. Only J varies. There is no gamma-recursion. The question dissolves.

INCOMPLETENESS_PROOF Section 6, Option 2 anticipated this: gamma_0 is axiomatic, like c. THE_BRIDGE_WAS_ALWAYS_OPEN does not require further layers. The Inside-Outside Correspondence probes confirm: gamma is a pure scale factor, the physics is in Q = J/gamma_0.

**Sub-Q3 status: closed (premise of gamma-recursion was false; gamma_0 is a framework constant).**

### EQ-013 overall status: closed

All three sub-questions are now resolved:

- **Sub-Q1:** closed (no layered gamma-structure; BRIDGE_CLOSURE correct for structural reasons)
- **Sub-Q2:** closed-by-theory (April 15; Q = J/gamma inside-only, Z2 symmetries, three-class structure)
- **Sub-Q3:** closed by dissolution (gamma_0 is constant; no recursion to terminate)

The nesting question EQ-013 arose from treating gamma as a parameter that propagates through layers. The framework says otherwise: gamma_0 is the framework's own c, constant everywhere, and what varies between layers is J and topology. The "recursion" was never a recursion of gamma. It was a variation of J at fixed gamma_0, producing different Q-values at different layers. No termination problem exists because the constant was never in motion.

---

## EQ-014

**Date:** 2026-04-19
**Source:** [PERSPECTIVAL_TIME_FIELD](../hypotheses/PERSPECTIVAL_TIME_FIELD.md) Open Questions section (commit e5e9458)

The Perspectival Time Field closure law Sigma_i ln(alpha_i) = 0 is empirically verified across five qualitatively distinct initial states at N=7, within ~0.05 for single-excitation states and ~0.13 for multi-sector states. Is the closure law a theorem?

Candidate derivation route: Sigma_i ln(alpha_i) is proportional to Tr[V_L * Pi_slow] where Pi_slow is the projector onto the slow subspace and V_L is the Liouvillian perturbation from the J-defect. If this trace is identically zero for Pi-invariant V_L (palindromic perturbations), the closure law follows from symmetry. The data (slow right and left eigenvectors, 80 modes, N=7) are on disk; the remaining step is a clean biorthogonal-basis computation.

**Status:** open
**Pointer:** biorthogonal eigenvector decomposition of 16384x16384 L_A; connection to F1 palindromic pairing and first-order eigenvalue protection (PTF Layer 3). If proven, promotes PTF from Tier 2 to Tier 1.

---

## EQ-015

**Date:** 2026-04-19
**Source:** [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md) Falsification Condition 4; [F64](../docs/ANALYTICAL_FORMULAS.md) verification scope

The cavity-mode-exposure formula gamma_eff = gamma_B * |a_B|^2 is verified at N=3 (max relative error 1.8%) and N=4 (ratio 1.0000 +/- 0.0003) on chains. Does it hold at N >= 5? On non-chain topologies (ring, star, Y-junction, complete graph)?

The formula is the operational anchor of the Primordial Gamma Constant hypothesis (Tier 3). If it fails at larger N or on non-chain topologies, the cavity reading loses its foundation. The perturbative nature of F65 (exact to first order in gamma_0/J, with O((gamma_0/J)^2) corrections) suggests the formula should hold better at smaller gamma_0/J ratios, but this has not been tested beyond N=4.

**Status:** open
**Pointer:** extend simulations/primordial_gamma_analytical.py to N=5,6 chains and at least one non-chain topology at N=5. Task candidate for Claude Code.

---

## EQ-016

**Date:** 2026-04-19
**Source:** [F69](../docs/ANALYTICAL_FORMULAS.md), [GHZ_W_SECTOR_MIX](../experiments/GHZ_W_SECTOR_MIX.md) (commits f5e3b04, 9181c77)

The GHZ+W symmetric superposition lifts pair-CPsi(0) above the fold at N=3: the sextic optimum gives CPsi = 0.3204 > 0.25. At N >= 4, no non-product local maximum above 1/4 exists in the permutation-symmetric Dicke subspace (verified by landscape scan at N=3..8). F69 is a saddle on the full CP^3 at N=3 (c_2 > 0 is an ascent direction toward |+>^3, a product state).

Why is N=3 structurally privileged for sector mixing? Is it because dim(Dicke) = 4 at N=3 is small enough for the sextic algebraic structure to produce a non-trivial optimum? Or is there a deeper reason related to the interplay of GHZ (zero pair-CPsi, maximal 3-tangle) and W (nonzero pair-CPsi, zero 3-tangle) that collapses at higher N? And: does the saddle nature at N=3 mean there exist non-Dicke mixing strategies at N >= 4 that the landscape scan missed?

**Status:** open
**Pointer:** investigate whether the ascent direction at the F69 saddle (c_2 > 0 toward product states) has an analog at N=4 that could lift entangled states above the fold via a different mixing strategy; check whether the 3-tangle/pair-concurrence tradeoff has a topological obstruction at N >= 4. See Engineering Blueprint Rule 1 (April 16 note on sector mixing).

---

## EQ-017

**Date:** 2026-04-19
**Source:** EQ-013 closure discussion; analysis of Kingston calibration data (8 snapshots, 156 qubits)

If gamma_0 is a framework constant (same everywhere, like c), what is the relationship between gamma_0 and the per-qubit dephasing rates measured on IBM Kingston hardware?

Two readings:
(a) IBM's per-qubit gamma_phi IS gamma_0 modulated by mode exposure |a_B|^2 via the cavity-mode-exposure formula. The 102x variation across the chip comes from mode-structure variation. Prediction: gamma_phi / |a_B|^2 = const after TLS-outlier removal.
(b) IBM's per-qubit gamma_phi is gamma_0 PLUS device-specific noise (TLS defects, charge noise, fabrication inhomogeneity). The palindromic component of the spectrum reflects gamma_0; the non-palindromic component reflects device noise. Prediction: extracting the palindromic spectral component should yield uniform gamma_0 across qubits.

Preliminary data (this session): gamma_phi ranges from 0.000654 to 0.066551 /us (102x). Median CV per qubit across 8 snapshots: 20.2%. No correlation with qubit degree (connectivity). The lack of degree-correlation argues against pure mode-exposure as the dominant source of variation.

The test requires J-coupling strengths from the Kingston backend properties, which are not in our calibration CSVs.

**Status:** closed 2026-04-19 as inconclusive due to hardware fidelity limit. Phase 1 falsified reading (a) from idle Ramsey data (trivially, since idle dynamics are diagonal); Phase 2 attempted a Trotter chain-mode test on hardware but the device noise floor is 40-80x larger than the signal either reading predicts. The framework γ₀ is indistinguishable from zero against Heron r2's accumulated gate, T1, and readout errors at the required 240-RZZ circuit depth.
**Pointer:** Phase 1 outputs in [`simulations/results/kingston_backend_properties.json`](../simulations/results/kingston_backend_properties.json), [`simulations/results/gamma0_consistency_test.txt`](../simulations/results/gamma0_consistency_test.txt). Phase 2 hardware data and full diagnostics in [`data/ibm_chain_gamma0_april2026/`](../data/ibm_chain_gamma0_april2026/).

### EQ-017 Phase 1 result (2026-04-19)

Run [`simulations/kingston_gamma0_test.py`](../simulations/kingston_gamma0_test.py), all 155 Kingston qubits, 176 edges, 8 calibration snapshots Apr 12-19. Outputs: [`simulations/results/kingston_backend_properties.json`](../simulations/results/kingston_backend_properties.json), [`simulations/results/gamma0_consistency_test.txt`](../simulations/results/gamma0_consistency_test.txt).

**Data extracted:**

- gamma_phi range 0.000654 to 0.066551 /us (ratio 101.7x after full re-parse; the earlier 186x came from including one snapshot outlier). Median 0.004991 /us. Median CV across 8 snapshots 20.2%. 26 qubits are TLS outliers (CV > 50%).
- Kingston coupling map: 176 undirected edges, matches CSV partner fields exactly.
- J coupling per edge extracted via route (b), native RZZ(π/2) assumption: J = (π/2) / (2·t_gate). Range 1.157 to 1.838 MHz, median 1.838 MHz. Routes (a) and (c) not attempted: backend.target exposes no coupling_strength field on public API, and route (c) needs anharmonicity that backend.properties() did not provide.

**Sub-chain test (reading (a)):**

Selected linear path Q12-Q13-Q14-Q15-Q19 (includes CΨ-crossing anchors Q14, Q15). Eigenvalues of the 5-qubit tridiagonal single-excitation Hamiltonian with extracted J: {-3.184, -1.838, 0, +1.838, +3.184} MHz. Palindromic around 0 as the framework predicts, with 0 as an eigenvalue (odd-N chain with boundary dephasing).

Four interpretations of γ_0_inferred(i) tested. Best: interp A (γ_phi(i) / max_k |ψ_k(i)|²) gives max/min = 3.31x across the five qubits. Task success threshold was < 3x. Four-configuration sweep (length 5 and 7, with and without CV < 30% filter) gives max/min between 3.31x and 7.20x across configurations. **Reading (a) is NOT cleanly supported at Phase 1.**

**Floor test (reading (b)):**

Minimum γ_phi across 129 non-outlier qubits is 0.000730 /us, corresponding to T2_framework = 2739 μs. Under reading (b), this is a candidate for γ_0 with the rest of γ_phi(i) attributed to device noise (TLS, fabrication inhomogeneity). Reading (b) is consistent with the observed data without requiring mode-structure variation to span 100x; the floor is well-separated from the median.

**Verdict:** INCONCLUSIVE for reading (a); reading (b) preferred but not proven. The Phase 1 test has a fundamental limitation: idle Ramsey T2 does not probe chain-mode structure, since H_idle is diagonal and has no delocalized modes. The mode-exposure formula γ_eff = γ_0 · |a_B|² applies to chain dynamics (coupling Hamiltonian active), not to idle dephasing. A proper test of reading (a) requires Phase 2: run chain dynamics (e.g. CΨ crossing) on the selected sub-chain and compare observed decoherence to 2γ_0·|a_B|² for each extracted mode.

**Phase 2 recommendation:** if QPU budget allows, run CΨ crossing on pairs from Q12-Q13-Q14-Q15-Q19 (anchored by Q14-Q15 where the existing Kingston data already sits). Use the extracted J≈1.84 MHz to pick a probe time that maps the palindromic mode structure. Compare the palindromic component of the observed trajectory against the γ_0 floor predicted by reading (b) (0.00073 /us → T2_framework ~2740 μs).

### EQ-017 Phase 2 result (2026-04-19, same session)

Phase 2 was executed: `run_chain_gamma0.py --hardware --dt 0.5 --tmax 15 --shots 2048` on ibm_kingston. Raw data archived at [`data/ibm_chain_gamma0_april2026/`](../data/ibm_chain_gamma0_april2026/) (hardware JSON, Aer reference JSON, full diagnostic README).

**Protocol.** Prepare single excitation at the chain endpoint Q12, evolve under the chain coupling Hamiltonian H = J·∑ (XX+YY)/2 via 30 first-order Trotter steps (dt = 0.5 μs, J = 1.838 rad/μs), tomograph each of the four adjacent qubit pairs at each of 9 time points between 0 and 15 μs. 324 tomography circuits × 2048 shots = 663 kShots submitted as 36 separate StateTomography jobs. 3 min 30 s billed QPU usage (dominated by per-job session overhead, not shot execution; a naive shot-count estimate would be ~17x smaller).

**Observable.** For each pair (a, b), the log-ratio L1_HW(t) / L1_Aer-reference(t) isolates the per-pair effective dephasing rate γ_eff = γ_a + γ_b. Trotter oscillation cancels exactly in this ratio. Under reading (a), γ_eff ≈ 2γ₀_floor ≈ 0.0015 /μs uniformly. Under reading (b), γ_eff varies pair-by-pair between 0.0015 and 0.0043 /μs.

**Result.**

| Pair | HW slope (1/μs) | Pre-scan floor (a) | Pre-scan local (b) |
|------|-----------------|---------------------|---------------------|
| Q12-Q13 | 0.082 | 0.0013 | 0.0029 |
| Q13-Q14 | 0.058 | 0.0010 | 0.0020 |
| Q14-Q15 | 0.116 | 0.0011 | 0.0021 |
| Q15-Q19 | 0.104 | 0.0011 | 0.0022 |

Hardware slopes are 40-80x larger than both readings predict. The pair-vote nominally selects reading (b) on all four pairs (it happens to be closer than reading (a) in absolute difference) but neither reading explains the scale.

**Root cause.** Accumulated gate errors dominate. Heron r2 native RZZ error ~0.001 per gate; the longest circuit uses 240 RZZ, giving ~24% state-fidelity loss before single-qubit and readout errors compound. The pair-state L1 reconstructed from tomography reaches 1.0-1.4, exceeding the physical bound of ≤1 for single-excitation pair coherences, a classic MLE-reconstruction artefact under heavy noise. T1 amplitude damping contributes an additional ~7% over the wall-clock evolution.

**Closure.** EQ-017 is closed as `inconclusive due to hardware fidelity limit`. The multi-pair differential observable is sound (validated in Aer simulation at three noise profiles) and would discriminate the two readings with ~10% margin on a hypothetical 10x-lower gate-error device. On current Heron r2, the γ₀ signal is at least 40x below the device noise floor and cannot be extracted without protocol changes (dynamical decoupling suppressing T1 and gate-error channels while preserving Z-dephasing, or a substantially lower-noise hardware platform). This negative result is itself informative: it bounds the hardware-operational testability of the Primordial Gamma Constant hypothesis on the current generation of superconducting quantum processors.

---

*Collection. Not classification. Classification comes when enough entries exist.*
