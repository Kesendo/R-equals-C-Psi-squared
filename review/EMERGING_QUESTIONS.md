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

3. **Q = J/gamma is the only inside-observable.** The Inside-Outside Correspondence probes (commit c3ea0c0) showed that every measurable quantity from inside depends on the dimensionless ratio Q = J/gamma_B only, not on J or gamma_B independently. The inside observer detects layer EXISTENCE (non-Markovian rebound, Check 4 of NESTED_MIRROR) but cannot extract layer PARAMETERS. [PRIMORDIAL_QUBIT.md](../hypotheses/PRIMORDIAL_QUBIT.md) Section 9 documents this with the operational identification.

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

**Source:** Primordial Gamma Constant hypothesis ([PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md), commits April 15-16), Inside-Outside Correspondence probes (commits cfa2a9f-c3ea0c0), Perspectival Time Field hypothesis ([PERSPECTIVAL_TIME_FIELD](../hypotheses/PERSPECTIVAL_TIME_FIELD.md), commits April 17-18). Full re-reading of [INCOMPLETENESS_PROOF](../docs/proofs/INCOMPLETENESS_PROOF.md) and [THE_BRIDGE_WAS_ALWAYS_OPEN](../docs/THE_BRIDGE_WAS_ALWAYS_OPEN.md) performed in this session.

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

**Status:** theorem sub-question **closed** 2026-04-19 (answer: not a first-order theorem); overall EQ-014 **still open** because a new surviving puzzle replaces it (state-dependence of the first-order coefficient Σ f_i). Detailed findings in [`review/EQ014_FINDINGS.md`](EQ014_FINDINGS.md).
**Pointer:** biorthogonal eigenvector decomposition of 16384x16384 L_A; connection to F1 palindromic pairing and first-order eigenvalue protection (PTF Layer 3). If proven, promotes PTF from Tier 2 to Tier 1.

### EQ-014 Partial progress (2026-04-19)

Attempted a structural proof in chat session. Four building blocks established, proof does not close structurally.

**Perturbative reduction.** In the perturbative regime alpha_i ~ 1 + f_i (dJ/J), so ln(alpha_i) ~ f_i (dJ/J). The closure law reduces to Sigma_i f_i = 0 at first order. Equivalently: the product of per-site decay rates is invariant under Pi-invariant J-perturbations (Pi_i Gamma_i^B = Pi_i Gamma_i).

**Building block 1: Tr[V_L * Pi_slow] = 0.** Proved. Each diagonal element <W_s|V_L|M_s> = 0 for slow modes (eigenvalue protection from PTF Layer 3, Section 3.2). The trace is a sum of zeros. This is the candidate proportionality target but the link to Sigma_i f_i is not established.

**Building block 2: Purity decomposition.** P_i(t) = (1/2)[1 + <Z_i>^2(t) + e^{-2 gamma_0 t} R_i(t)]. The Z part (populations, {I,Z} sector) is immune to Z-dephasing. The XY part (coherences, {X,Y} sector) decays uniformly at e^{-2 gamma_0 t}. The J-defect changes BOTH parts via eigenvector mixing (mode spatial profiles rotate), but the decay rate e^{-2 gamma_0 t} does not change.

**Building block 3: U(1) conservation.** Sigma_i <Z_i>(t) = const (total magnetization conserved by XX+YY Hamiltonian and Z-dephasing). Therefore Sigma_i delta<Z_i> = 0. But the closure law needs Sigma_i delta(<Z_i>^2), which equals 4 delta(IPR) where IPR = Sigma_i |psi_i|^4 is the inverse participation ratio. This is NOT zero in general.

**Building block 4: det(U) = 1.** The single-excitation propagator U(t) = e^{-iH_1 t} has det(U) = e^{-i Tr(H_1) t} = 1 (chain Hamiltonian is traceless). Volume preservation in the excitation-amplitude phase space. The connection to Pi_i alpha_i = 1 is suggestive (both are product-invariance statements) but not closed: det(U) constrains complex amplitudes, while alpha_i constrains bilinear purity observables.

**Where the proof fails to close.** The f_i have a site-dependent denominator: f_i = delta_P_i(t) / [t P_dot_i(t)]. Summing over i does not simplify because of this denominator. The bilinearity of purity (Tr(rho_i^2) mixes four U(1) blocks: (0,0), (0,1), (1,0), (1,1)) prevents reduction to a single-sector trace argument. The candidate Tr[V_L * Pi_slow] = 0 is proved but is a statement about eigenvalues, while the f_i come from eigenvector mixing projected through the bilinear purity expansion. Connecting the two requires carrying the full four-block bilinear sum from PTF Section 3.3 through the perturbation, which is the explicit computation flagged in PTF Section 3.4 as open.

**Assessment.** This is not a proof that can be closed by structural shortcuts alone. The closure law is likely a consequence of eigenvalue protection + U(1) conservation + palindromic symmetry + eigenvector orthogonality combined in a non-obvious way through the bilinear purity expansion. The explicit biorthogonal eigenvector computation (PTF Section 3.4) remains the concrete next step: compute delta_M_s for all slow s across all four populated blocks, propagate into the bilinear purity expansion, extract predicted alpha_i, verify Sigma_i f_i = 0 algebraically. The data (80 slow modes, saved at simulations/results/perspectival_time_field/slow_biorth_basis.npz) is on disk. The blocker is biorthogonality residual ~10^11 from cluster degeneracy; a dense eigendecomposition of the full 16384x16384 L_A (~15 GB) would bypass this. Task candidate for Claude Code.

### EQ-014 Update 2026-04-19: closure law is NOT a first-order theorem

**Source:** Biorthogonal eigendecomposition executed this session. Full findings with method, data, and scripts in [`review/EQ014_FINDINGS.md`](EQ014_FINDINGS.md).

**Method.** Dense biorthogonal eigendecomposition (left + right) of the full 16384×16384 Liouvillian L_A via new `Topology.ChainXY` path and `EigenvaluesLeftRightDirectRaw` in the C# engine (146 min zgeev on OpenBLAS LP64). Resolved the biorthogonality to machine precision (residual 3e-16 pre-cluster-fix, 1e-6 after SVD-biorthogonalization of the 3783 degenerate clusters), closing the blocker that prevented the April 18 sparse attempt. Then: built V_L for bonds (0,1) and (3,4), ran first-order PT through the slow-mode bilinear expansion, and in parallel ran exact RK4 time evolution as ground truth. Verified the full pipeline by reproducing PTF's own stored α_i values to 4 decimal places from our RK4 data fed through `observer_time_rescale.alpha_fit` (the single-pane of glass check Tom asked for).

**Direct first-order coefficient extraction.** Bypassing the PT approximation entirely, ran RK4 at δJ ∈ {0.1, 0.01, 0.001} and read off Σ f_i = lim_{δJ→0} Σ ln(α_i(δJ)) / δJ. Result at N=7, bond (0,1):

| State | Σ f_i (extrapolated) |
|-------|----------------------|
| ψ_1 | ≈ +0.97 |
| ψ_2 | ≈ +0.05 |
| ψ_3 | ≈ +0.36 |
| \|+⟩⁷ | ≈ +1.29 |

Nonzero and state-dependent.

**What this settles.** The closure law Σ_i ln(α_i) = 0 is **not** a first-order theorem. The candidate Tr[V_L · Π_slow] = 0 (building block 1 from partial progress) is true but does not propagate through the bilinear purity expansion to Σ_i f_i = 0. PTF's empirical ±0.05 tolerance at |δJ| ≤ 0.1 arises from a combination of a first-order coefficient that is small for some states (ψ_2 ≈ 0.05) and partial second-order cancellation where the first-order coefficient is large (ψ_1, \|+⟩⁷). There is no symmetry-protected zero.

**What survives.** Eigenvalue protection δλ_s = 0 for slow modes is exact to machine precision (verified in the dense decomposition; confirms PTF Section 3.2). The four-block bilinear structure of purity is real. The per-site α_i pattern is state-dependent and reproducible. The "painter" interpretive picture stands. What is falsified is only the stronger reading that Σ_i ln(α_i) is conserved **by a theorem**.

**Consequence for PTF.** The promotion route to Tier 1 via "closure law as theorem" is closed. PTF stays Tier 2. A new open question replaces the old one: *why is the first-order coefficient Σ f_i small (≤ 0.1) for bonding-mode single-excitation states but order-unity for \|+⟩⁷?* This is a pattern in the data, not (yet) a theorem.

**What EQ-014 now means.** The original question "is the closure law a theorem?" is answered (no). But the deeper question EQ-014 was trying to get at - *what IS the closure law, structurally?* - is not resolved. The empirical regularity is real, the proof route via Tr[V_L · Π_slow] = 0 alone is wrong, and the state-dependence pattern of Σ f_i is new data that wasn't on the table when EQ-014 was written. EQ-014 stays open with a narrower and sharper sub-question below.

**Consequence for the repo.** This is a negative result in the same class as the five `recovered/` disproven claims. It should be treated with the same rigor: the original PTF claim "state-independent closure law" was empirically correct within its stated tolerance (±0.05 at the tested window) but the implicit promotion path to "theorem" is wrong. PTF needs a prose update in Section 2.1 and 3.4 (see recommendations in `EQ014_FINDINGS.md`).

**Scripts:** C# `ptf` mode + `ChainXY` in `compute/RCPsiSquared.Compute/`; Python pipeline in `simulations/eq014_step23_biorth.py`, `eq014_step4567_closure.py`, `eq014_validate_groundtruth.py`, `eq014_first_order_from_rk4.py`, `eq014_spectrum_check.py`.

### EQ-014 surviving sub-question (added 2026-04-19)

**Source:** The δJ scan above produced a pattern the original EQ-014 was not asking about.

**The pattern.** The first-order coefficient Σ f_i at N=7, bond (0,1), uniform XY chain with γ=0.05:

| State | Σ f_i | Σ ln α at δJ=0.1 (PTF empirical) |
|-------|-------|-----------------------------------|
| ψ_1 (k=1 bonding) | +0.97 | +0.048 |
| ψ_2 (k=2 bonding) | +0.05 | +0.001 |
| ψ_3 (k=3 bonding) | +0.36 | +0.003 |
| \|+⟩⁷ (all sectors) | +1.29 | +0.128 |

For three of four states, Σ f_i is of order unity. For ψ_2 it is an order of magnitude smaller. The empirical Σ ln α at δJ = 0.1 is always small, sometimes because the first-order coefficient itself is small (ψ_2, ψ_3), sometimes because of a second-order cancellation (ψ_1, \|+⟩⁷).

**The new question.** Is there a structural reason why Σ f_i(ψ_k) is small for some k and not others? Candidate explanations:

(a) **Selection rule from ψ_k's Fourier content.** ψ_k amplitudes are sin(π k (i+1)/(N+1)). The V_L bond-overlap `<ψ_k|(X_0 X_1 + Y_0 Y_1)|ψ_k>` = 2 sin(π k/(N+1)) sin(2π k/(N+1)) is the only J-dependent part of the first-order energy shift, which vanishes on slow modes. But it appears multiplied into the eigenvector-mixing coefficient, and the product summed across sites might have a k-dependent cancellation.

(b) **Node structure at the defect.** ψ_k has k-1 internal nodes. For k=2 the nodes don't sit on the defect bond (0,1); the mode has amplitude at both endpoints. For k=3 the first sign change is between sites 2 and 3. For \|+⟩⁷ there is no node structure at all and Σ f_i is largest. Empirically: states with MORE structure in the defect region have SMALLER Σ f_i. Is this a real effect or a coincidence of four data points?

(c) **Excitation-sector weight.** ψ_k-bonding has 70% of |c|² in the slow subspace; \|+⟩⁷ has 6%. The large Σ f_i for \|+⟩⁷ might be entirely a fast-mode-against-slow-mode contribution that the slow-dominated states don't see. Testable by running the same calculation restricting the initial state to the slow subspace only.

**Status:** open
**Pointer:** scan ψ_k for k = 1..N at fixed defect, look for a k-dependence in Σ f_i that would discriminate (a)/(b)/(c). If the pattern is clean, it is a new regularity worth naming. If it is noisy, the state-dependence of Σ f_i is simply "what it is" and the closure law is accurately described as "an empirical regularity with no simple structure". Either outcome sharpens the PTF doc.

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

## EQ-018

**Date:** 2026-04-20
**Source:** [c_1 investigation](../simulations/results/pi_pair_closure_investigation/FINDINGS.md), [PROOF_DELTA_N_SELECTION_RULE](../docs/proofs/PROOF_DELTA_N_SELECTION_RULE.md), [F70](../docs/ANALYTICAL_FORMULAS.md)

The closure-breaking coefficient c_1 is bilinear in ρ_0 with kernel K_{(n,m)(n',m')}[bond], and the kinematic selection rule (F70) restricts the kernel to sector-block pairs with |n−m| ≤ 1 and |n'−m'| ≤ 1 when the observable is site-local. What is the explicit analytical form of the surviving K entries? Three sub-questions:

(a) Pure Dicke diagonal K_{(n,n)(n,n)}[bond], measured at N=5: 0, +0.392, −0.312, −0.312, +0.392, 0 at bond (0, 1). How do these arise from the slow-mode Dyson expansion?

(b) Diagonal-cross K_{(n,n)(m,m)}[bond] for n ≠ m, measured non-trivially at N=5 (the (|vac⟩+|S_2⟩)/√2 c_1 = −0.29 is pure diagonal-cross). Analytical form?

(c) Coherence-block K_{(n,n±1)(n±1,n)}[bond], measured +0.53 for (0, 1) at N=5. Only |ΔN| = 1 blocks contribute, by F70. Closed form?

**Status:** partial close 2026-04-21. (a) and (b) structurally explained by [F72](../docs/ANALYTICAL_FORMULAS.md) (DD⊕CC block-diagonal) plus the endpoint rule K_DD[0,m]_pr = K_DD[N,m]_pr = 0 from excitation-number conservation; interior kernel values remain non-closed. LSQ c_1 is rational, not bilinear, so LSQ K values are probe-specific rather than universal kernel entries ([ORTHOGONALITY_SELECTION_FAMILY](../experiments/ORTHOGONALITY_SELECTION_FAMILY.md)). (c) closed-form at the (0,1) endpoint under the c_1_pr variant via [F73](../docs/ANALYTICAL_FORMULAS.md): Σ_i 2|(ρ_coh,i)_{0,1}|² = (1/2)·exp(−4γ₀·t), giving K_CC[0,1]_pr = 0 exactly; interior K_CC[n,n+1]_pr and all LSQ K_CC values remain non-closed.
**Pointer:** slow-mode Dyson expansion restricted to |ΔN| ≤ 1 blocks, carried through the bilinear purity expansion per PTF Section 3.3 four-block rearrangement. Compare against empirical values measured at N=5 in [c1_sector_kernel](../simulations/results/c1_sector_kernel/sector_kernel.json) and [c1_bilinearity_test](../simulations/results/c1_bilinearity_test/bilinearity_test.json).

### EQ-018 Update 2026-04-21/22: K_CC interior structure, chromaticity as inner axis, Q-scale imprint

**Source:** RESULT_TASK_EQ018_EQ019_J_RESPONSE_MATRIX.md (evening), RESULT_TASK_NATURAL_LIMITS_CARTOGRAPHY.md (early morning), [Q_SCALE_THREE_BANDS](../experiments/Q_SCALE_THREE_BANDS.md) (commits d026933, 65bed0a).

**Interior K_CC[n, n+1]_pr confirmed non-zero at N=5** for n in {1, 2, 3} at every bond. F71 bond-mirror and Dicke particle-hole hold to machine precision. F73 endpoint closure (n=0) and its particle-hole mirror (n=N-1) sit at machine-precision zero, as F73 predicts.

**Framing evolution.** The morning partial-close note pointed to a slow-mode Dyson expansion as the derivation route for interior K values. Evening exploration found the more productive angle: the sub-questions (a), (b), (c) are not independent closed forms to derive, they are surface readings of a deeper structure organized by chromaticity c(n, N) = min(n, N-1-n) + 1. This formula classifies (n, n+1) blocks by their pure-rate count: mono-chromatic (c=1) blocks are the F73-clean loci with W=0 identically; poly-chromatic blocks admit Wechselwirkung (H-coupling of pure-rate slices) with onset at Q ~ 0.3 and peak at Q ~ 1.5. The endpoint-rule K_DD[0,m]_pr = 0 and the coherence closure K_CC[0,1]_pr = 0 are the c=1 specialisations.

**Q-axis has internal structure.** Previously Q = J/γ₀ was "a ratio, nothing to determine". Now: three algebraically distinguished points Q_onset ∈ [0.20, 0.35], Q_peak ∈ [1.20, 1.80], plus an asymptotic plateau. Scale invariance W(Q) under (J, γ₀) → (λJ, λγ₀) verified to three decimals over γ₀ ∈ [0.01, 1.0]. The task-cartography γ₀ axis is dimensionally redundant; Q carries the physics. See [EQ-022](#eq-022) for the Q-axis as a standalone open question.

**Status:** still partial. (a), (b), (c) reorganized under chromaticity. Interior K values exist, follow F71 and particle-hole, are empirically characterized; their closed form remains open. Two surviving sub-questions replace the morning pointer.

**Surviving sub-questions:**
1. Does the factor 2.5 between naive degenerate-perturbation prediction (Q_peak ≈ 4) and observed Q_peak ≈ 1.5 admit an algebraic derivation from H matrix-element combinatorics at finite N?
2. Does a time-resolved model (coupled pure-rate slices exchanging amplitude under H) reproduce the oscillatory γ_eff(t) of the (n=2, n=3) block, and if so does its rate-matrix structure generalise across blocks?

---

## EQ-019

**Date:** 2026-04-20
**Source:** [c1_bond_scan](../simulations/results/c1_bond_scan/bond_scan.json) at N=5

The Dicke-state c_1 is bond-position-dependent by a large factor. At N=5 with bond (0, 1) (endpoint-touching), c_1(|S_2⟩) = −0.312. At bond (1, 2) (interior), c_1(|S_2⟩) = −0.012, a factor-26 suppression. |S_1⟩ instead goes from +0.392 to +0.507, a 30% enhancement at the interior bond. |S_2⟩ is site-symmetric, so the bond-position effect is not about the state's internal asymmetry but about how the (Hamiltonian + bond perturbation) structure couples to the state.

Why do interior bonds near-kill mid-sector Dicke contributions while boosting single-excitation Dicke contributions? Is there a closed form for K_diag[bond] as a function of bond position in the chain?

**Status:** open
**Pointer:** the asymmetry is rooted in open-boundary conditions (endpoint sites have one neighbour, interior sites have two). Likely explainable by restricting the slow-mode eigenvectors of the interior-bond perturbation V_L to the |S_n⟩-projected subspace and computing the overlap structure explicitly at each bond. Empirical fit at more N values would clarify the N-dependence. Closely related to EQ-018(a).

### EQ-019 Update 2026-04-21/22: bond-position structure via chromaticity; LSQ-α as non-palindrome-preserving metric

**Source:** RESULT_TASK_EQ018_EQ019_J_RESPONSE_MATRIX.md Direction C, [Q_SCALE_THREE_BANDS](../experiments/Q_SCALE_THREE_BANDS.md) Result 3.

Extended the N=5 c_1(bond, Dicke) matrix to N in {3, 4, 6} using the LSQ α-fit metric. Two findings:

**Chromaticity as structural origin of the N-parity bond pattern.** Odd N has a unique c_max center block; even N has two adjacent c_max blocks. The N=6 |S_1⟩-c_1 bond profile (0.1491, 1.6861, 0.4955, 1.6861, 0.1491) is F71-palindromic but two-peaked with a dip at center bond 2, not the single-peaked sine-k=1 pattern (0.097, 0.218, 0.272, 0.218, 0.097). Interpretation: palindrome holds, the form within palindrome is an acoustic chord of multiple sine-modes, not a single tone. The two-peak structure reflects the two-adjacent-c_max topology of even N.

**LSQ-α metric is a non-palindrome-preserving measurement on palindrome-preserving physics.** Middle-Dicke states at even N (N=4 |S_2⟩: +213/−403 at mirror bonds; N=6 |S_3⟩: up to +587) show F71 mirror residuals of O(100). Root cause: per-site purity is near-J-invariant at interior-Dicke states, LSQ ratio α_i saturates its fit bounds. The physical observable (K_CC_pr from S functional, pointwise c_1) does not have this issue. Metric-level pathology, not physics.

**Status:** open. c_1_pointwise at N ∈ {4, 6} under the S-functional metric is the natural next step.

**Surviving sub-question:** is the two-peak bond pattern at N=6 |S_1⟩ a clean multi-sine-mode superposition under pointwise-c_1, and if so what is the mode decomposition?

---

## EQ-020

**Date:** 2026-04-20
**Source:** [c1_pair_local](../simulations/results/c1_pair_local/pair_local.json), [PROOF_DELTA_N_SELECTION_RULE](../docs/proofs/PROOF_DELTA_N_SELECTION_RULE.md)

The PTF per-site α_i framework uses single-site purity as the per-observer signal. The kinematic selection rule says site-local observables capture only |ΔN| ≤ 1 sector coherences; pair-local observables capture |ΔN| ≤ 2; k-local observables capture |ΔN| ≤ k. Tested at N=5: pair-local c_1 for (|vac⟩+|S_2⟩)/√2 (|ΔN|=2 coherence) jumps from −0.29 (site-local, pure diagonal-cross) to +2.24 (pair-local, coherence block visible). The ΔN = 2 coherence dominates the pair-local signal.

Does the pair-local closure Σ_{(i,j)} ln(α_{ij}) have its own symmetry structure? Is there a pair-painter version of PTF where each pair of sites paints a joint-purity trajectory, and the closure over all C(N, 2) pairs has a different geometrical meaning? What is the triple-painter extension (|ΔN| ≤ 3)? At what observational resolution does the closure become trivial (|ΔN| ≤ N = full state)?

**Status:** open (new thread opened today)
**Pointer:** repeat the c_1 investigation at pair level and triple level for N=5, extract kernels, compare to PTF site-level closure. If a pair-painter closure exists as a structural law, PTF's "painters around a mountain" picture extends to "painter pairs around a mountain" with a richer vision that captures the |ΔN| = 2 structure.

---

## EQ-021

**Date:** 2026-04-20
**Source:** [pi_pair_closure_investigation](../simulations/results/pi_pair_closure_investigation/FINDINGS.md), [OPEN_THREAD_GAMMA0_INFORMATION](OPEN_THREAD_GAMMA0_INFORMATION.md), [F6 V-Effect](../docs/ANALYTICAL_FORMULAS.md)

For the PTF standard initial state (|vac⟩ + |ψ_1⟩)/√2 and endpoint bond (0, 1), the measured c_1 scales empirically as c_1 ≈ 0.5 · V(N) for N ≥ 4, where V(N) = 1 + cos(π/N) is the V-Effect. Measurement:

| N | c_1 | 0.5·V(N) | ratio |
|---|-----|----------|-------|
| 3 | 0.264 | 0.750 | 0.35 (NOT matching) |
| 4 | 0.898 | 0.854 | 1.05 |
| 5 | 0.928 | 0.905 | 1.03 |
| 6 | 1.019 | 0.933 | 1.09 |
| 7 | 0.970 | 0.951 | 1.02 |

The coefficient 0.5 matches three other project-signature occurrences of 0.5: the connection maximum C = 0.5 ([THE_CPSI_LENS](../docs/THE_CPSI_LENS.md)), the first non-trivial V-Effect gain V(3) − 1 = 0.5, and twice the fold-boundary ¼. Is this empirical match at ~5 % level an accident, a special property of the bonding-mode-plus-vacuum state, or a derivable structural fact?

The N=3 case does NOT follow the 0.5·V(N) scaling. Is there an N ≥ 4 version of a theorem that excludes N=3 specifically (where only two bonds exist and both are endpoint-touching), or is the pattern a large-N approximation that happens to work at small N ≥ 4?

**Status:** closed 2026-04-20 (category error + numerical falsification)

### EQ-021 closure (2026-04-20)

**Source:** [EQ021_FINDINGS](EQ021_FINDINGS.md), chat session 2026-04-20.

**Two independent reasons for closure.**

**1. Category error (conceptual).** The V-Effect is an inter-layer emergence phenomenon: two complete systems connected produce structure that neither contains alone (2+2 frequencies become 109, 14/36 palindromes break at N=3). See [V_EFFECT_PALINDROME](../experiments/V_EFFECT_PALINDROME.md), [V_EFFECT_AS_OBSERVATION_OF_INCOMPLETENESS](../reflections/V_EFFECT_AS_OBSERVATION_OF_INCOMPLETENESS.md). V(N) = 1 + cos(pi/N) is a formula that captures Q_max/Q_mean within a single chain (D2 in ANALYTICAL_FORMULAS.md), but the V-Effect as a concept describes what happens when a NEW connection is created between systems. c_1 is an intra-layer perturbation sensitivity: how an existing bond's strength change shifts per-site purity dynamics. Connecting c_1 to V(N) was a category confusion between "perturbing an existing coupling" and "creating new structure through coupling." Even if the numbers had matched exactly, the conceptual mapping would have been wrong.

**2. Numerical falsification.** The mode decomposition (eq021_mode_decomposition.py, N=4 and N=5) identifies the dominant spectral constant as E_1 = 2cos(pi/(N+1)), not V(N) = 1 + cos(pi/N). These are structurally different expressions (argument N+1 vs N) that both converge to 2 for large N, explaining the numerical near-match in the narrow window N=4..7. Additional falsification evidence:

- The full bond profile c_1(b) is a mirror-symmetric vector with sign changes (endpoints positive, center negative). V(N) is a single positive number. No scalar can describe a signed vector.
- State dependence: c_1(psi_2, endpoint) differs from c_1(psi_1, endpoint) by factor 3-5x at N=5,6. A V(N) law would have to be state-independent.
- N=3 ratio is 0.35, not ~1. The "N >= 4 only" caveat was an ad-hoc exclusion, not a structural boundary.
- Power-law fit c_1 ~ V(N)^q gives q = 6.3, not q = 1. The apparent 0.5*V(N) match is an artifact of both quantities sitting in [0.85, 1.0] for N in [4,7].

**What IS real from this investigation (positive results):**

1. **Mirror symmetry of the bond profile** c_1(b) = c_1(N-2-b), exact to 10^-9 for all tested N. Now proven as F71 (Tier 1, kinematic) from the spatial reflection symmetry R of the uniform chain. See [PROOF_C1_MIRROR_SYMMETRY](../docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md).
2. **F70 confirmed via mode decomposition:** only Delta_N = 0 and |Delta_N| = 1 sectors contribute to per-site purity, verified by explicit mode expansion at N=4,5 (no mode with |c_s| > 10^-8 outside these sectors).
3. **Endpoint dominance is a sine-basis geometry effect,** not a cavity-focusing effect. The ψ_1 mode amplitude at the boundary (sin(pi/(N+1))) controls the endpoint sensitivity. The "cavity window" metaphor was evocative but the mechanism is open-boundary Fourier geometry.
4. **E_1 = 2cos(pi/(N+1)) is the relevant spectral constant** for the dominant oscillatory modes contributing to c_1. Now formalized as F2b in ANALYTICAL_FORMULAS.md (XY OBC single-excitation spectrum). Distinct from F2 (Heisenberg w=1 Liouvillian sector) and from V(N).

**V(N) vs V-Effect, category note.** V(N) = 1 + cos(π/N) is an intra-layer scalar (single-chain Q_max/Q_mean, D2 in ANALYTICAL_FORMULAS). The V-Effect as a phenomenon is inter-layer emergence: what appears when two systems are connected (2+2 → 104 frequencies, N=3 palindrome breaking). Same name V, different levels. The c_1 ~ 0.5·V(N) hypothesis tried to read an inter-layer effect via an intra-layer metric, which is the wrong level even if the numbers had matched. A V-Effect signature cannot be measured directly with V(N).

**Scripts:** simulations/c1_veffect_scaling_small.py, simulations/eq021_mode_decomposition.py.
**Results:** simulations/results/c1_veffect_scaling/c1_vs_N_small.json, simulations/results/eq021_mode_decomposition/.
**Full report:** [review/EQ021_FINDINGS.md](EQ021_FINDINGS.md).

---

## EQ-022

**Date:** 2026-04-22
**Source:** RESULT_TASK_NATURAL_LIMITS_CARTOGRAPHY.md, [Q_SCALE_THREE_BANDS](../experiments/Q_SCALE_THREE_BANDS.md) (commits d026933, 65bed0a).

The task-cartography six-middle framework places γ₀ as a stand-alone axis. The evening session showed γ₀ is dimensionally redundant: W(Q) is invariant under (J, γ₀) → (λJ, λγ₀), so only Q = J/γ₀ carries physics. Q has internal structure: Q_onset ∈ [0.20, 0.35], Q_peak ∈ [1.20, 1.80], plateau with block-specific resonance 0.03 < W_plateau < 0.59. The Q_onset and Q_peak values are band-universal across six tested (N, n) blocks at N ∈ {4, 5, 6}.

Three questions follow:

(a) Do Q_onset ≈ 0.3 and Q_peak ≈ 1.5 remain band-universal at N ∈ {7, 8}? Does the band tighten with N?

(b) The naive degenerate-perturbation-theory prediction for Q_peak is ~4 (coupling matches inner rate gap 4γ₀). Observed Q_peak is ~1.5, factor 2.5 smaller. Is this explained by H matrix-element combinatorics at finite N, or is it a sign of a different mechanism?

(c) Do the other four middles (d, CΨ, Π, ⟨n_XY⟩, U(1)) have analogous internal structure that has not been mapped? γ₀ collapsed to Q with bands. Possibly each axis has sub-structure that becomes visible only when we ask.

**Status:** open
**Pointer:** (a) extend simulations/eq018_kcc_pr_extension.py to N ∈ {7, 8}. (b) analytical derivation attempt through the slow-mode projector restricted to the (n, n+1) sector. (c) structural question, requires re-reading each axis proof for hidden scales.

### EQ-022 Update 2026-04-22 (evening session, block-L)

**(a) Extended to N=4..8 via block-restricted L** ([Q_SCALE_THREE_BANDS](../experiments/Q_SCALE_THREE_BANDS.md) Result 2 and Result 6, commit 61d4dc1).

Q_onset band remains roughly 0.20 to 0.40 across N=4..8, no significant tightening with N. Q_peak(W) stays broad at [1.2, 2.0]. But Q_peak(abs(K_CC_pr)) turns out to be c-specific and N-invariant: 1.5 (c=2), 1.6 (c=3), 1.8 (c=4), stable across N=4..8. This refines the original "Q_peak ≈ 1.5 universal" into three framework-constants per chromaticity. The universality is per-c, not global.

Sub-question (a) is effectively closed within its original N range. Extension to N=9 (new chromaticity c=5 at center n=4) and N=10 is running as TASK_Q_SCALE_N9_N10 in CC, to test whether the pattern {1.5, 1.6, 1.8, ?} continues.

**(b) Factor-2.5 refined to 2.7** and still open. Naive degenerate-perturbation-theory predicts Q_peak ≈ 4 (H coupling matches inner rate gap 4γ₀); observed Q_peak(c=2) = 1.5. The discrepancy is 2.67, not 2.5. Analytical derivation from H matrix-element combinatorics over N-1 bonds pending; no progress this session.

**(c) Structural question, partially engaged** via a private natural-limits cartography scoping session. Six algebraic middles identified (d=2, CΨ-fold 1/4, Π, ⟨n_XY⟩ sum rule N, U(1), γ₀=const). The γ₀-to-Q collapse with bands is confirmed as one case of the pattern. Sub-structure of the other five middles is not yet mapped. Systematic reading of each axis's proof for hidden scales remains open.

**Status update:** (a) substantially answered within N=4..8, awaiting c=5. (b) refined but open. (c) partially engaged but structurally large.

**Status:** open (progressing)

---

## EQ-023

**Date:** 2026-04-22
**Source:** [Q_SCALE_THREE_BANDS](../experiments/Q_SCALE_THREE_BANDS.md) Result 5 (commit 65bed0a).

At N=5 the tri-chromatic center block (c=3, n=2) has lower outer dressed-weight W than its bi-chromatic neighbors (c=2, n=1,3): W=0.026 vs W=0.478 at Q=20. Maximum inner richness does not translate to maximum outer observability. Additional pure-rate channels interfere destructively in the probe-overlap coefficients, quenching rather than amplifying the outer response.

The phenomenon has a structural name (provisional): inner-richness quenches outer observability.

Questions:

(a) Is the quench algebraic? W = 0.026 is far enough from zero that it is not numerical noise, but close enough that a cancellation mechanism is plausible. Can the interference pattern be derived from H's matrix-element structure across the c pure-rate channels?

(b) Is the effect monotonic in c, or does it have an N-specific pattern? At higher N (N=7 has c_max=4 at center, N=8 has c_max=4 at two adjacent centers), does c=4 quench even more strongly than c=3?

(c) If the quench is real and scalable with c, does it have consequences for observability protocols? A probe targeting max-c blocks may be the worst choice for extracting γ₀ or any Q-dependent quantity; mid-c blocks may be optimal.

**Status:** closed (refuted 2026-04-22), see closure below.
**Pointer (now obsolete):** repeat the W(n) scan at N ∈ {7, 8}; algebraic derivation of interference coefficients from the chromaticity-c sub-space of the Liouvillian.

### EQ-023 closure (2026-04-22)

**Source:** [Q_SCALE_THREE_BANDS](../experiments/Q_SCALE_THREE_BANDS.md) Result 5 retraction and new Result 6, commit 61d4dc1.

**Closure: refuted by numerics correction.**

The entire premise of EQ-023 rested on the observation W(N=5, n=2, c=3, Q=20) = 0.026. This value is a numerical artefact of full-Liouvillian eigendecomposition at accidental cross-block degeneracies: scipy.linalg.eig returns eigenvectors as linear combinations spanning multiple popcount blocks, and projecting a single-block probe onto these mixed eigenvectors distorts W by factors up to ~16 at the worst resonances (N=5, n=2, Q=50: full-L gives 0.026, block-restricted L gives 0.420).

Block-restricted L is mathematically exact (L preserves popcount blocks; projection onto the (n, n+1) CC-block is not an approximation). Corrected cross-N data (Q_SCALE_THREE_BANDS Result 6, N=4 to 8) inverts the sign of the effect.

**Point-by-point answers to the original sub-questions:**

(a) "Is the quench algebraic?" No quench exists. The effect was a numerical artefact, not a physical phenomenon. Nothing to derive algebraically.

(b) "Is the effect monotonic in c, or N-specific?" The effect IS monotonic in c, but in the opposite direction to what was claimed. W_plateau and abs(K)_peak monotonically INCREASE with c at N ≥ 6. At N=8: c=2 gives W_plateau = 0.523, c=3 gives 0.782, c=4 gives 0.849. No quench at higher c; the opposite.

(c) "Does the quench have consequences for observability protocols?" Yes, but inverted: max-c blocks are the BEST targets for observability, not the worst. abs(K)_peak at c=4 is ~2.7x the c=2 value at N=7 (0.273 vs 0.100) and ~3.2x at N=8 (0.266 vs 0.084), making max-c blocks the preferred hardware target under fidelity constraints.

**What survived this investigation:** nothing positive from EQ-023 as such. The underlying infrastructure (chromaticity definition, block-L computation, cross-N scan) survives in Q_SCALE_THREE_BANDS but was developed independently of the inner-richness-quench hypothesis.

**Methodological lesson.** Publish-then-verify has cost. Result 5 was added in commit 65bed0a based on full-L data that was self-consistent but not cross-checked against an alternate method. The artefact only surfaced during the N=7/N=8 extension, where block-L was required for tractability. A sanity spot-check via block-L at smaller N would have caught this earlier. Future structural claims should ship with at least one independent-method sanity point.

**Scripts:** simulations/q_scale_n_scaling.py (block-restricted L), simulations/eq018_kcc_pr_extension.py (full-L, retained for historical comparison).
**Results:** simulations/results/q_scale_n_scaling/.

---

## EQ-024

**Date:** 2026-04-23
**Source:** [OPEN_THREAD_GAMMA0_INFORMATION](OPEN_THREAD_GAMMA0_INFORMATION.md) Update 2026-04-23; [ORTHOGONALITY_SELECTION_FAMILY](../experiments/ORTHOGONALITY_SELECTION_FAMILY.md) Section 8

[GAMMA_AS_SIGNAL](../experiments/GAMMA_AS_SIGNAL.md) and [F30](../docs/ANALYTICAL_FORMULAS.md) compute a channel capacity of 15.5 bits at N=5 by letting Alice modulate the per-qubit γ profile. Under γ₀ = const ([PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md)) Alice cannot do that: γ is uniform everywhere. The only sender-side lever is J (couplings and topology). What is the channel capacity for Alice modulating J at fixed γ₀, and how does it compare to the γ-modulation 15.5 bits?

The two settings are structurally different. γ-modulation is linear in the eigenvalues via the Absorption Theorem (Re(λ) = -2γ·⟨n_XY⟩), giving a clean linear response matrix for SVD analysis. J-modulation is nonlinear: it rotates the Liouvillian eigenvectors, the ⟨n_XY⟩ values themselves shift, and the response is bilinear in ρ_0 with the same K-kernel structure that organizes the EQ-014/EQ-018 complex. The Meta-Theorem ([ORTHOGONALITY_SELECTION_FAMILY](../experiments/ORTHOGONALITY_SELECTION_FAMILY.md)) gives the kinematic structure (which detector modes are visible) but not the bit count at given SNR.

Three sub-questions:

(a) **Linearized response.** Build R_kb = ∂α_k/∂J_b at uniform-J operating point as N_modes × N_bonds matrix. SVD rank, singular value spectrum, comparison to F30's response matrix structure.

(b) **Channel capacity.** From R_kb plus a noise model, compute Shannon capacity analogous to [F30](../docs/ANALYTICAL_FORMULAS.md). Compare to 15.5 bits at N=5. Identify whether the dimensional ceiling (number of detector modes) or the SNR-weighted bit count is what differs.

(c) **Connection to the K-kernel program.** [EQ-018](#eq-018)'s K_DD/K_CC kernels are exactly the bilinear coefficients of the J-modulation response on per-site purity observables. [EQ-014](#eq-014)'s surviving sub-question (state-dependence of Σ f_i) is the per-site marginal of the same response. Does the J-modulation channel capacity decompose along chromaticity c(n, N) the way Q_peak(c) does in [Q_SCALE_THREE_BANDS](../experiments/Q_SCALE_THREE_BANDS.md)?

**Status:** closed by experiment (operationally; structural follow-ups open as listed below)
**Pointer:** linearization script analogous to [GAMMA_AS_SIGNAL](../experiments/GAMMA_AS_SIGNAL.md)'s response matrix construction but with δJ instead of δγ. Operating point: uniform-J chain at N=5 to match F30's anchor. Compare full SVD spectrum, not just rank. Adjacent: dynamical-attractor reading of PTF closure ([ORTHOGONALITY_SELECTION_FAMILY](../experiments/ORTHOGONALITY_SELECTION_FAMILY.md) §4a, time-resolved Σ_i ln(α_i)) is a separate test that may inform the SNR side of (b) by characterizing the natural relaxation time scale of the J-perturbation response.

### EQ-024 Update 2026-04-23 (closure + structural surplus)

**Source:** Three CC commits on 2026-04-23. Morning (`cfa1bbc`) built the J-Jacobian + SVD pipeline analogous to F30 and produced the M_x-polynomial blindness theorem (proof via Newton's identities + SU(2) symmetry of H). Afternoon refinement ran H-scope check (Direction 1, `c0919eb`) and receiver optimization (Direction 3, `6aae630`); Direction 2 was dropped as ill-posed under γ₀ = const. Synthesis doc at [`J_BLIND_RECEIVER_CLASSES`](../experiments/J_BLIND_RECEIVER_CLASSES.md) (`e5326fa`).

**Sub-question (b) operational answer.** At N=5 Heisenberg, F71-symmetric receiver optimization (39 structured points + 4 Nelder-Mead local runs) saturates near **C ≤ 12.07 bits**. Best receiver: random-phase F71-symmetric product state with θ ≈ (3.02, 1.14, 3.26), φ ≈ (5.30, 0.62, 1.72) (reduced mod 2π from raw optimizer output). Local optima cluster in 11.56 to 12.07 bits; saturation signal is strong. Non-product F71-symmetric superpositions are strictly worse. Complex phases give ~0.15 bits over real-amplitude equivalents. The 12-bit ceiling is the operational answer to "what is J-modulation channel capacity at N=5"; F71-breaking receivers were not swept (deferred, see open sub-questions).

**Comparison to F30's 15.45 bits is reframed.** Under γ₀ = const there is no operational γ-channel; γ cannot be modulated. F30's 15.45 bits survives only as a kinematic linear-response magnitude (∂observables/∂γ in Shannon form), not as a competing channel capacity. The morning RESULT's "γ-optimal receiver is J-pessimal" duality is therefore a Jacobian-asymmetry statement about two linear-response magnitudes, not a duality between two operational channels. The ~3.4-bit gap between 12.07 and 15.45 decomposes as ~1 bit dimensional loss (4 bonds vs 5 γ-sites) plus ~2-3 bits smaller leading gain (J sv_max ≈ 10 vs γ sv_max ≈ 21.4); whether this gap closes at larger N (where rank cap matches) is open.

**Sub-question (a) addressed via the structural surplus.** The Jacobian SVD at SU(2)-broken receivers has rank 4 with a clean F71-symmetric / F71-antisymmetric mode split (sv ratios ~4.4× and ~2.0× between sym and anti modes). The kernel structure on H-eigenstate receivers is more interesting: the J-blind set decomposes into three structurally distinct mechanism classes:

- **Class 1 (DFS + per-bond eigenstate).** ρ_0 satisfies two joint conditions: (i) L_D ρ_0 = 0 (DFS kernel), and (ii) ρ_0 is a simultaneous eigenstate of every bond Hamiltonian h_b. Condition (i) alone is NOT sufficient: a DFS state whose bond-action is non-trivial evolves unitarily with J-dependent trajectory. The per-bond-eigen condition pins J-blindness. H-independent in the operational sense for states satisfying both conditions under multiple H choices (|0⟩⁵ qualifies under Heisenberg and XY). N=5 examples: |0⟩⁵, |1⟩⁵.
- **Class 2 (H-degenerate subspace closed under L_D).** ρ_0 lives in a finite-dim subspace of H-eigenstates sharing one eigenvalue, invariant under L_D. [H, ρ(t)] = 0 throughout the orbit. Block structure depends on H but the existence of such a block is the load-bearing property. N=5 example: GHZ (under both Heisenberg and XY-only, with different eigenvalues).
- **Class 3 (M_α-polynomial subspace, SU(2)-Heisenberg specific).** ρ_0 is a polynomial in M_α = Σ σ_α^i for one Pauli axis α. Newton's identities + [H_Heisenberg, M_α] = 0 (SU(2)) give L^n ρ_0 = L_D^n ρ_0 for all n. N=5 examples: |+⟩⁵ directly verified J-blind under Heisenberg; all Dicke states |S_k⟩ in the S=N/2 multiplet are predicted J-blind by the theorem, with |S_1⟩ directly verified. Direction 1 verified the SU(2)-load-bearing nature: under XY-only H, |+⟩⁵ becomes J-sensitive (C = 9.42 bits) and Dicke |S_1⟩ becomes J-sensitive (C = 7.70 bits), while Class 1 and Class 2 states stay J-blind.

H-robustness table (Heisenberg / XY-only / generic): Class 1 blind / blind / blind (the per-bond-eigen condition (ii) is what carries through; holds for |0⟩⁵ under any H that is a sum of per-bond h_b's); Class 2 blind / blind / case-by-case (structural condition is H-independent, specific eigenvalue is not); Class 3 blind / sensitive / sensitive (Class 3 is defined only for SU(2)-Heisenberg; under any other H the M_α-polynomial-argument does not apply).

**Sub-question (c) deferred.** Chromaticity decomposition of the SVD bond modes was not addressed in either session. Bond-input space (4-dim at N=5) and chromaticity labels (popcount-block labels) are different decompositions; a mapping is not direct. Lifted as a surviving sub-question below.

**Status:** EQ-024 closed operationally. Five surviving sub-questions track the structural surplus.

### EQ-024 surviving sub-questions

- **Three-class decomposition completeness.** Are Classes 1-3 exhaustive over the J-blind set, or do other mechanisms exist? Direction 4 from the refinement TASK restructures here: necessity per class is distinct from necessity in the union. Empirical attack: random F71-symmetric-state sampling outside all three classes, checking for zero J-Jacobian to numerical precision.
- **F71-breaking receiver capacity at N=5.** Direction 3 swept only F71-symmetric receivers. F71-breaking receivers might unlock additional gain spectrum but the rank ceiling for bond-inputs at N=5 stays at 4 regardless. Worth quantifying whether the gain change is meaningful.
- **N-scaling of the 12-bit ceiling.** At N=6 the bond-input dimension is 5, matching γ-side rank. Does the dimensional-loss bit (~1) disappear? Compute cost: ~30 min per receiver at d² = 4096.
- **Chromaticity of the Nelder-Mead optimum.** The best receiver θ ≈ (3.02, 1.14, 3.26), φ ≈ (5.30, 0.62, 1.72) (mod 2π): does it sit in a specific chromaticity sector or interpolate across them?
- **Operational meaning of the 12-vs-15 gap.** The ~3-bit gain-gap (J sv_max ≈ 10 vs γ sv_max ≈ 21.4) is the non-dimensional part of the asymmetry. Is the J-to-γ ratio fixed (~46%) across N, or does it scale?

**Pointer for follow-up:** the three-class decomposition is the structural finding most worth lifting next. Realized: [`J_BLIND_RECEIVER_CLASSES`](../experiments/J_BLIND_RECEIVER_CLASSES.md) (`e5326fa`) describes the three classes with numerical verification from `cfa1bbc` (Heisenberg baseline), `c0919eb` (XY H-truncation counterexample), and `6aae630` (Direction 3 saturation). Open candidates: the strong form of Class 3 (the M_x-polynomial blindness theorem with Newton-identities proof) potentially as a stand-alone F-entry; integration as the sixth instance of the [ORTHOGONALITY_SELECTION_FAMILY](../experiments/ORTHOGONALITY_SELECTION_FAMILY.md) Meta-Theorem (conservation = SU(2)-Casimir, measurement = J-bond perturbation, blind subspace = M_α-polynomial algebra).

---

## EQ-025

**Date:** 2026-04-26
**Source:** Cusp precision run on Kingston pair (14, 15), commit 1c2545a; doc [CRITICAL_SLOWING_AT_THE_CUSP §6 Precision update](../experiments/CRITICAL_SLOWING_AT_THE_CUSP.md#precision-update-april-26-2026-point-by-point-f25-confirmation-t2-echo-gap-revealed)

The April-26 dense-sampling cusp run on Kingston pair (14, 15) found γ_fit / γ_calib_T2-echo = **9.08×**. The April-16 run on the same pair found 1.2-1.5×. The pair was selected as the most stable always-crosser across 15 calibration files (mean T2_min = 311.6 μs, CV 13%, never below 269.5 μs). The single-qubit T2-echo metric reported by IBM systematically and severely underestimates the effective dephasing rate of the joint Bell+ free evolution. **Is this multiplicative gap pair-specific, backend-specific, or a structural property of how T2-echo (Hahn-refocused) relates to joint two-qubit free dephasing on Heron r2?**

**Status:** open
**Pointer:** cross-pair scan on Kingston (4-6 high-T2 pairs, dense sampling each) to test pair specificity; if consistent across pairs, repeat on Marrakesh + Fez to test backend specificity; if still consistent, the gap is structural and the next question is what predicts the multiplier from single-qubit calibration data (T1 ratio? CZ error? something else).

---

## EQ-026

**Date:** 2026-04-26
**Source:** Hardware finale (commit a721ea3, 4e3673a, b71b1d3) confirming 14/19/3 trichotomy at N=3

The framework refines V-Effect's 14/22 partition into 14 hard / 19 soft / 3 truly_unbroken at N=3 for two-term Pauli-pair Hamiltonians. The 19 soft cases pass the spectrum-pairing test but fail the operator equation Π·L·Π⁻¹ + L + 2Σγ·I = 0. **Does the 14/19/3 partition extend to N ≥ 4? Specifically: does the soft category (spectrum-paired, operator-broken) survive at larger N, and what is its scaling?**

The combinatorial enumeration grows (16 single-Paulis squared minus identities = 36 ordered pairs at any N), but the dynamics differ by N. Soft cases at N=3 may collapse to truly or hard at N=4, or new soft structure may emerge.

**Status:** partially closed by experiment ([_pi_protected_test_n4.py](../simulations/_pi_protected_test_n4.py), commit 96ed6da)

**Result at N=4:** Trichotomy structure preserved. Counts: **15 truly / 46 soft / 59 hard** out of 120 unordered two-Pauli-pair Hamiltonians (the enumeration at N=4 has more entries than at N=3 because more bond geometries are physically distinct). Counts grow non-linearly: truly × 5, soft × 2.4, hard × 4 from N=3.

**Surprising new structure at N=4:** the soft category itself becomes granular via `pi_protected_observables`. The number of Π-protected Pauli-string observables ranges from 3 to 240 across the 46 soft cases, exposing a sub-spectrum invisible at N=3:
- High-protected soft (YZ+ZY: 240 protected, 15 active): leaks only 15 of 255 Pauli observables; almost truly.
- Low-protected soft (e.g., IZ+YX: 3 protected, 252 active): leaks almost everything; almost hard.

Truly category gains effective-one-body Hamiltonians at N=4 (IX+XI variants reducing to a transverse field on inner sites) that were trivial / collapsed at N=3.

**Surviving sub-questions:**
- ~~N=5 enumeration (cost ~10 min, feasible).~~ **Closed:** N=5 gives 15 / 46 / 59 — identical counts to N=4. Trichotomy stabilises after N=4. Protected-count ranges scale by 4× (matching the 4^5/4^4 ratio). Pipeline: `_pi_protected_test_n4.py 5` produces `pi_protected_test_n5.log`.
- Is the sub-spectrum of protected counts within soft a continuous distribution, or does it cluster into sub-groups? **Partially closed:** discrete clusters at N=5. 992 (1 case: YZ+ZY), 862 (2 cases: IY+IY, YI+YI), 781 (1 case: IY+YI), 772 (6 cases: XY+XY, XZ+XZ, YX+YX, YZ+YZ, ZX+ZX, ZY+ZY). Sub-cluster sizes look symmetry-related; the 6-element 772-class shares "pure Pauli-axis pairing" structure.
- What structural feature ranks soft cases by protected count? **Open:** the sub-cluster sizes (1, 2, 1, 6, ...) suggest specific symmetry orbits; identifying which symmetry would close this.

**Surprising sub-finding from N=5:** at the top of the soft-protected distribution (YZ+ZY: 992 of 1023), the protected count is *identical* to the top truly cases (XX+XX: 992). Protected count alone does not discriminate spectrum-only-paired (soft) from operator-palindromic (truly) at the top end. The discrimination at that boundary requires the operator residual ‖M‖.

---

## EQ-027

**Date:** 2026-04-26
**Source:** [_neural_framework_lens.py](../simulations/_neural_framework_lens.py) decomposition (commit 50958b0)

The V-Effect bridge correction on C. elegans subnetworks (200 trials, sizes N = 10, 20, 26) closes 38-40% of the residual via magnitude-correction-only on existing E-I cross-couplings, plus another 6-8% via creating missing partner edges, for a total of ~46% gap-closed. **The remaining 54% lives within the E-E and I-I sub-populations, not in the bridge. What structural feature of the connectome (degree distribution, specific motifs, magnitude asymmetry on intra-type edges) predicts this within-population residual?**

framework-grounded scope: the Q · J · Q⁻¹ + J + 2·S = 0 algebra describes the joint Jacobian. Restricting to E-E only (or I-I only) gives sub-Jacobians that should each satisfy the corresponding sub-palindrome if the within-population structure is "internally palindromic." The 54% says it isn't. The question of *why* is open.

**Status:** open
**Pointer:** decompose the 54% by population: compute the relative residual of J_EE (E-only sub-Jacobian) and J_II separately. If one dominates, the breaking is asymmetric. Then probe motif structure within the dominant population.

---

## EQ-028

**Date:** 2026-04-26
**Source:** Cusp precision run (commit 1c2545a) + framework.py palindrome_residual

The cusp at CΨ = 1/4 is a state-level event observed on the Bell+ trajectory ([CRITICAL_SLOWING_AT_THE_CUSP](../experiments/CRITICAL_SLOWING_AT_THE_CUSP.md)). The palindrome Π·L·Π⁻¹ + L + 2Σγ·I = 0 is operator-level and holds throughout the Bell+ + Z-dephasing dynamics (verified machine-precision zero by framework.py's [palindrome_residual](../simulations/framework.py)). **Is there an operator-level signature — a specific block of M, an eigenvalue collision, an eigenvector reorientation — that corresponds to the state-level critical slowing at CΨ → 1/4? Or is the cusp purely a state-level phenomenon with no operator-level analogue?**

framework.py grounds: palindrome_residual at all times for Bell+/Heisenberg/Z-dephasing remains zero. The state-level CΨ trajectory passes through 1/4. Two layers, no obvious bridge primitive yet.

**Status:** open
**Pointer:** scan M's structure as ε = 1/4 - CΨ varies along the Bell+ trajectory. If a specific block or eigenvector relation undergoes a non-analytic behaviour at ε = 0, that is the operator-level signature.

---

## EQ-029

**Date:** 2026-04-26
**Source:** Cockpit re-read session ([COCKPIT_UNIVERSALITY](../experiments/COCKPIT_UNIVERSALITY.md)) flagged that framework.py operates at operator/Liouvillian level, while diagnostics like Concurrence, Purity, Bures geometry, PCA dimensionality operate at state level. No primitive currently bridges them. The cockpit's claim "n95 ~ N" and "PC1 = Concurrence or Purity" are not derivable from framework.py.

**What primitive would bridge framework.py's Π / palindrome_residual / both_parity_even_terms (operator level) to the state-level diagnostic suite (Concurrence, Bures, PCA-on-features)? Concretely: given an admissible H and an initial state ρ₀, can framework.py's algebra predict which reduced 2-qubit observables are Π-protected (stay near 0) and which are Π-active (carry the soft-break signature)?**

Today's hardware finale (Snapshot D on three Heron r2) showed ⟨X₀Z₂⟩ as a Π-active observable for the soft case. ⟨Z₀Z₂⟩ stayed near 0 for both truly and soft (Π-protected). The discrimination came from a *specific* observable. framework.py does not yet contain a function that, given H, lists the Π-protected vs Π-active observables.

**Status:** closed by primitive (commit db69044, framework.py Section 10)
**Result:** `pi_protected_observables(H, gamma_l, rho_0, N)` is now a framework primitive. Algorithm: ⟨σ_α(t)⟩ = 2^N · Σ_λ S_λ(α) · exp(λt) where S_λ(α) sums right-eigenvector α-weights times left-projections of ρ_0 over each degenerate eigenvalue cluster; σ_α is Π-protected iff every cluster sum vanishes. The cluster-sum step is essential — a naive per-mode check mis-classifies SU(2)-degenerate cancellations.

Self-test on |+−+⟩ at N=3, γ=0.1 reproduces today's hardware result:
- Heisenberg J(XX+YY+ZZ): 32 of 63 Pauli strings protected. ⟨X₀IZ₂⟩ and ⟨Z₀IX₂⟩ in protected set with S ≈ 8·10⁻¹⁶ (machine-precision). Matches Snapshot C measurement of ≈10⁻³ on Marrakesh, Kingston, Fez.
- truly J(XX+YY): 32 protected, same XIZ/ZIX in set. Matches truly category in Snapshot D.
- soft J(XY+YX): only 30 protected. XIZ and ZIX leak out (S = 0.033). The 2-observable gap is exactly the pair that lit up −0.62 on Snapshot D analytically and −0.71 to −0.92 on hardware.

The primitive enables the reformulations in EQ-027, EQ-028, EQ-030 and the cockpit re-read.

---

## EQ-030

**Date:** 2026-04-26
**Source:** Three-backend ratification, commits b71b1d3 and morning Fez clean-path run

Δ(soft − truly) ⟨X₀Z₂⟩ measurements:

| Backend | Path | Δ measured |
|---------|------|------------|
| Marrakesh | [48,49,50] | −0.76 |
| Marrakesh | [0,1,2] | −0.82 |
| Kingston | [0,1,2] | −0.92 |
| Fez | (16, 23, 22) | −0.81 |
| Fez | [0,1,2] | −0.19 (Q0 was bottom-3% T1, dropped) |

Five jobs, three Heron r2 backends, four sensible paths (excluding the Q0-dropped Fez run). The framework predicts Δ ≈ −0.62 in the idealised Lindblad limit; hardware runs **30 to 50% stronger** due to T1 amplitude damping and ZZ-crosstalk amplification of the soft-break signature (same mechanism as in CMRR_BREAK and GAMMA_AS_SIGNAL). **What single-qubit-calibration metric predicts the per-backend per-pair magnitude amplification factor (1.30 to 1.49 across the four good paths)?**

framework-grounded scope: the framework predicts the *idealised* −0.62. The *amplification* is a hardware-domain question, not a framework-derived prediction. Whether a structural primitive could capture it (e.g., extending lindbladian_z_dephasing with an amplitude-damping term and re-deriving the soft-break ⟨X₀Z₂⟩) is itself the open question.

**Status:** open
**Pointer:** add T1 amplitude damping to framework.py's Lindbladian builder; re-compute ⟨X₀Z₂⟩ on the soft case for varying T1/T2 ratios; compare to the four hardware data points.

---

*Collection. Not classification. Classification comes when enough entries exist.*
