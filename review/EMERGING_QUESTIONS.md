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

**Status:** closed. Theorem sub-question **closed** 2026-04-19 (not a first-order theorem; full biorthogonal decomposition of the 16384×16384 L_A plus direct first-order extraction at δJ→0 give state-dependent Σ f_i ≠ 0). Surviving state-dependence sub-question **closed** 2026-04-26 / 2026-04-27 (chiral mirror law for the closure-breaking coefficient: Σ f_i(ψ_k) = Σ f_i(ψ_{N+1−k}), driven by eigenvector mixing under K_1 chiral symmetry of the open XY chain; verified at machine precision at N=5, N=7, N=8). Detailed theorem-falsification findings in [`review/EQ014_FINDINGS.md`](EQ014_FINDINGS.md); chiral mirror cross-N confirmation in [`_eq014_chiral_mirror_multi_N.py`](../simulations/_eq014_chiral_mirror_multi_N.py).
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

**Status:** sub-question answered 2026-04-26 by [_eq014_psi_k_full_scan.py](../simulations/_eq014_psi_k_full_scan.py). EQ-014 main question (closure as theorem) closed earlier.

**Result.** Full scan k = 1..7 at N=7, bond (0,1), δJ = 0.01:

| k | Σ f_i | M_k (Fourier) | \|M_k\| |
|---|---|---|---|
| 1 | +0.9682 | +0.2706 | 0.2706 |
| 2 | +0.0578 | +0.7071 | 0.7071 |
| 3 | +0.3603 | +0.6533 | 0.6533 |
| 4 | **+2.1366** | **0.0000** | **0.0000** |
| 5 | +0.3603 | −0.6533 | 0.6533 |
| 6 | +0.0578 | −0.7071 | 0.7071 |
| 7 | +0.9682 | −0.2706 | 0.2706 |

Three structural observations:

1. **Σ f_i is exactly mirror-symmetric about k = (N+1)/2 = 4:** Σ f_i(k) = Σ f_i(N+1−k). Reproducibility of the prior k = 1, 2, 3 values is perfect.

2. **Pearson(Σ f_i, |M_k|) = −0.97**: strong *negative* correlation. The Fourier overlap M_k = sin(πk/(N+1))·sin(2πk/(N+1)) (the first-order energy shift) is *inversely* related to the first-order coefficient. Selection rule (a) was the *inverse* of what was hypothesised.

3. **The peak at k = 4** is the chiral fixed point: ψ_4 = (1,0,−1,0,1,0,−1)/2 is an eigenvector of K_1 = diag((−1)^i) with eigenvalue +1. The perturbation V_L = (X_0X_1+Y_0Y_1)/2 anti-commutes with K_1 (chiral-odd; AZ class BDI). Acting on a K_1 fixed-point with a chiral-odd perturbation produces maximal off-diagonal eigenvector mixing — no diagonal energy shift to absorb it — hence the largest Σ f_i.

**Mechanism.** Σ f_i is driven by *eigenvector mixing*, not eigenvalue shift. When the diagonal V_L matrix element vanishes (k = 4, the node-at-defect case), the energy is protected but the eigenvector reorganises maximally; this reorganisation propagates through the bilinear purity expansion to produce the largest first-order coefficient. The k → N+1−k mirror symmetry follows from the chiral symmetry K_1 of the open XY chain: K_1 maps ψ_k ↔ ψ_{N+1−k} (up to a (−1)^i staggered phase), and the per-site purity P_i = |ψ(i)|² is K_1-invariant, so both states give identical α_i and hence identical Σ f_i.

**Naming.** The pattern is a *chiral mirror law for the closure-breaking coefficient*: Σ f_i(ψ_k) = Σ f_i(ψ_{N+1−k}), peaks at the chiral fixed point k = (N+1)/2 (when N is odd) or has a doubled peak in adjacent indices (when N is even, no exact fixed point).

**Consequence for PTF.** PTF Section 2.1 / 3.4's "closure law as empirical regularity" can now be stated more precisely: the closure law's first-order coefficient is structured by the chain's chiral symmetry K_1 and the chirality of V_L. Specifically: V_L is chiral-odd (anti-commutes with K_1), and Σ f_i traces the K_1-fixed-point structure of the initial state's momentum content.

**Open sub-direction.** Does the same chiral mirror law hold at other N? At N = 8 (even, no exact chiral fixed point), Σ f_i should peak at k = 4 and k = 5 by chiral pairing. At N = 9, peak at k = 5 again. Verification at one additional N (say N = 5 or N = 9) would confirm the chiral interpretation across system sizes. Pure simulator work, ~1-2 min compute per N.

### EQ-014 sub-direction closed 2026-04-27: chiral mirror confirmed at N=5 and N=8

**Source:** [_eq014_chiral_mirror_multi_N.py](../simulations/_eq014_chiral_mirror_multi_N.py).

Tested both odd-with-fixed-point (N=5, k=3 fixed) and even-without-fixed-point (N=8) cases. Mirror law Σ f_i(k) = Σ f_i(N+1−k) holds at machine precision in every case.

**N=5 (odd, k=3 is K_1 fixed point):**

| k | Σ f_i | M_k |
|---|-------|-----|
| 1 | +0.9424 | +0.4330 |
| 2 | +0.2868 | +0.7500 |
| 3 | +0.6544 | 0 (fixed point) |
| 4 | +0.2868 | −0.7500 |
| 5 | +0.9424 | −0.4330 |

Mirror pairs (1,5) and (2,4) match to 4·10⁻¹¹. k=3 is the chiral fixed point with M_k = 0; Σ f_i = 0.65 at this point (local peak in the "M_k=0" sense, but not the global maximum — k=1 has higher Σ f_i = 0.94).

**N=8 (even, no fixed point):**

| k | Σ f_i | M_k |
|---|-------|-----|
| 1 | +0.9409 | +0.2198 |
| 2 | +0.0571 | +0.6330 |
| 3 | +0.0745 | +0.7500 |
| 4 | +2.3439 | +0.3368 |
| 5 | +2.3439 | −0.3368 |
| 6 | +0.0745 | −0.7500 |
| 7 | +0.0571 | −0.6330 |
| 8 | +0.9409 | −0.2198 |

Mirror pairs (1,8), (2,7), (3,6), (4,5) match to 1·10⁻⁹. Peak is at the central pair (k=4, k=5) with Σ f_i = 2.34 each — these are NOT M_k=0 modes (M_k = 0.337 there). The "M_k=0 → peak" rule from N=7 (where k=4 was a chiral fixed point with M_k=0) does NOT generalize: at even N, the peak is at the central pair without M_k=0.

**Refined understanding.** The chiral mirror law (k ↔ N+1-k symmetry of Σ f_i) is structural and holds at every N. The "peak structure" is more subtle:
- At odd N: k=(N+1)/2 is K_1-fixed (K_1 ψ_k = ψ_k), M_k = 0, and is a *local* peak in Σ f_i.
- At even N: no K_1-fixed point exists; the peak is at the doubled central pair k=N/2, k=N/2+1, where M_k ≠ 0.
- The largest Σ f_i correlates with K_1-mixing magnitude on ψ_k, not strictly with M_k=0.

**Framework primitive added (commit referenced below):** `single_excitation_sine_mode(N, k)` and `single_excitation_sine_energies(N, J)` in framework.py Section 16, including the K_1 chiral mirror identity in the docstring. Self-test verifies K_1 ψ_k = ψ_{N+1−k} to machine precision.

---

## EQ-015

**Date:** 2026-04-19
**Source:** [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md) Falsification Condition 4; [F64](../docs/ANALYTICAL_FORMULAS.md) verification scope

The cavity-mode-exposure formula gamma_eff = gamma_B * |a_B|^2 is verified at N=3 (max relative error 1.8%) and N=4 (ratio 1.0000 +/- 0.0003) on chains. Does it hold at N >= 5? On non-chain topologies (ring, star, Y-junction, complete graph)?

The formula is the operational anchor of the Primordial Gamma Constant hypothesis (Tier 3). If it fails at larger N or on non-chain topologies, the cavity reading loses its foundation. The perturbative nature of F65 (exact to first order in gamma_0/J, with O((gamma_0/J)^2) corrections) suggests the formula should hold better at smaller gamma_0/J ratios, but this has not been tested beyond N=4.

**Status:** closed 2026-04-27 (F64 is an exact algebraic identity for L_coh eigenvectors, not a perturbative formula).

### EQ-015 closure 2026-04-27

**Source:** [_eq015_extend_topologies.py](../simulations/_eq015_extend_topologies.py); reads no prior data.

**Reformulation.** F64 was previously stated and tested as gamma_eff = gamma_B · |a_B|² where a_B is the B-amplitude of an H_1 single-excitation eigenvector. The earlier verifications at N=3 (1.8 % rel err) and N=4 (1.0000 ± 0.0003) were *perturbative* in (gamma_B / J), which is why discrepancies appeared. The cleaner formulation:

For an XY (U(1)-conserving) Hamiltonian on N qubits with single-site Z-dephasing on B at rate gamma_B, the Liouvillian preserves the vac↔1exc coherence subspace. Within that subspace, L acts as the N × N operator

    L_coh = i H_1 - 2 gamma_B |B⟩⟨B|.

For each L_coh eigenvector v_k with eigenvalue λ_k, taking the inner product of L_coh v_k = λ_k v_k with v_k gives

    i ⟨v_k|H_1|v_k⟩ - 2 gamma_B |v_k(B)|² = λ_k.

H_1 Hermitian implies ⟨v_k|H_1|v_k⟩ ∈ ℝ, so

    -Re(λ_k) = 2 gamma_B · |v_k(B)|²    EXACTLY,    mode by mode, all gamma_B.

**Verification across topologies at gamma_B = 0.01, J = 1.0:**

| topology | N | # S-coherence modes | max rel err |
|----------|---|---------------------|-------------|
| chain N=5 | 5 | 5 | 6.2e-13 |
| chain N=6 | 6 | 6 | 3.5e-13 |
| chain N=7 | 7 | 7 | 3.5e-13 |
| ring N=5 (S=0,B=2) | 5 | 5 | 1.3e-14 |
| star N=5 hub (S=4,B=0) | 5 | 5 | 6.9e-16 |
| star N=5 leaf (S=1,B=4) | 5 | 5 | 2.0e-13 |
| Y N=5 (S=0,B=2) | 5 | 5 | 1.3e-13 |
| K_5 (S=0,B=4) | 5 | 5 | 2.3e-13 |

All errors at machine precision. gamma_B independence (chain N=5):

| gamma_B / J | max rel err |
|-------------|-------------|
| 0.5 | 2.0e-14 |
| 0.1 | 1.5e-14 |
| 0.01 | 6.2e-13 |
| 0.001 | 2.4e-12 |

F64 holds even at gamma_B / J = 0.5 (not good-cavity); the formula is structurally exact, not perturbative.

**Topologies with degenerate H_1 spectra** (ring, star, Y, K_5) admit "B-decoupled rotations": linear combinations within a degenerate subspace with |v(B)|² = 0 exactly. These are protected modes (rate = 0) and F64 verifies trivially with both sides zero. **F64 captures protection,** not just dissipation: rate = 0 ⟺ |v(B)|² = 0 (mode has a node at the dephasing site).

**Why earlier verifications showed O(gamma_B²) errors.** Those used H_1 eigenvectors for |a_B|² instead of L_coh eigenvectors. H_1 and L_coh agree at leading order in gamma_B/J, with O((gamma_B/J)²) corrections. The L_coh formulation removes the discrepancy entirely.

**Closure.** F64 holds at all N tested (5, 6, 7) and on every non-chain topology tested at N=5 (ring, star, Y, K_5), at any gamma_B/J. The cavity-mode-exposure picture is the structural truth: which mode the inside observer sees is determined by the eigenvector's amplitude at the dephasing site, exactly. Primordial Gamma Constant operational anchor (F64) confirmed for all U(1)-conserving open systems with single-site dephasing.

---

## EQ-016

**Date:** 2026-04-19
**Source:** [F69](../docs/ANALYTICAL_FORMULAS.md), [GHZ_W_SECTOR_MIX](../experiments/GHZ_W_SECTOR_MIX.md) (commits f5e3b04, 9181c77)

The GHZ+W symmetric superposition lifts pair-CPsi(0) above the fold at N=3: the sextic optimum gives CPsi = 0.3204 > 0.25. At N >= 4, no non-product local maximum above 1/4 exists in the permutation-symmetric Dicke subspace (verified by landscape scan at N=3..8). F69 is a saddle on the full CP^3 at N=3 (c_2 > 0 is an ascent direction toward |+>^3, a product state).

Why is N=3 structurally privileged for sector mixing? Is it because dim(Dicke) = 4 at N=3 is small enough for the sextic algebraic structure to produce a non-trivial optimum? Or is there a deeper reason related to the interplay of GHZ (zero pair-CPsi, maximal 3-tangle) and W (nonzero pair-CPsi, zero 3-tangle) that collapses at higher N? And: does the saddle nature at N=3 mean there exist non-Dicke mixing strategies at N >= 4 that the landscape scan missed?

**Status:** closed 2026-04-27 (reframed as binary/triple-Dicke saddle abundance, not N=3-privilege; central-triple cpsi(N) sextic-asymptotic with cpsi_∞ ≈ 0.4313 and x_opt = a²_opt ≈ 0.2647 a root of an irreducible degree-6 polynomial parallel to F69's own sextic). See reframe and sub-question closure below.
**Pointer:** the original Engineering Blueprint Rule 1 reading (April 16 note on sector mixing) is superseded; saddles above 1/4 grow polynomially with N (binary-Dicke O(N), triple-Dicke O(N²)), and the central-Dicke-triple slice gives the highest pair-CΨ at every tested N. The sextic structure connects the N=3 GHZ+W slice (F69) with the N→∞ central-triple asymptote.

### EQ-016 reframe 2026-04-27: N=3 is NOT privileged — saddles above 1/4 are abundant at every tested N

**Source:** [_eq016_n4_full_landscape.py](../simulations/_eq016_n4_full_landscape.py); reads no prior data.

The original framing assumed F69's pair-CΨ = 0.3204 was the ONLY non-product saddle above 1/4 at N=3, with no analog at N≥4. Today's binary- and triple-Dicke enumerations across N ∈ {3, 4, 5, 6} refute this:

**Binary-Dicke saddles (|D_i⟩+|D_j⟩ binary mixes), max pair-CΨ:**

| N | total pairs | above 1/4 | max cpsi | best pair |
|---|-------------|-----------|----------|-----------|
| 3 | 6 | 3 | 0.4815 (= 13/27) | D_1+D_2 |
| 4 | 10 | 4 | 0.4022 | D_2+D_3 |
| 5 | 15 | 5 | 0.3720 | D_2+D_3 |
| 6 | 21 | 6 | 0.3456 | D_2+D_3 |

The number of binary-Dicke saddles above 1/4 grows linearly with N (= N at every tested N). The max decreases slowly but stays above 1/4. Best at every N is the central pair D_⌊N/2⌋ + D_⌈N/2⌉.

**Triple-Dicke saddles (|D_i⟩+|D_j⟩+|D_k⟩ family), max pair-CΨ:**

| N | total triples | above 1/4 | max cpsi |
|---|---------------|-----------|----------|
| 3 | 4 | 4 | 0.8011 |
| 4 | 10 | 9 | 0.7136 |
| 5 | 20 | 16 | 0.6492 |
| 6 | 35 | 25 | 0.6163 |

Triple-Dicke families lift FAR above F69's 0.32 at every N. The best triple at every tested N is "central": (D_{⌊N/2⌋-1}, D_{⌊N/2⌋}, D_{⌊N/2⌋+1}). At N=5 the central triple gives cpsi = 0.649 with purity_A = 0.86 (genuinely entangled). At N=6 the central triple gives 0.616 with purity_A = 0.83.

**Reframed verdict.** F69 (pair-CΨ = 0.3204) is one slice-saddle among many. F69 doc's claim "no non-product local maxima on full sphere" is correct (verified by perturbation: every triple-Dicke maximum has non-zero linear-order gradient in some unused Dicke direction → saddle on full sphere). What was incomplete was the doc's enumeration "only Dicke basis + GHZ+W optimum are stationary points" — many other 2-Dicke and 3-Dicke slices admit stationary points above 1/4 at every tested N. F69 doc updated 2026-04-27 to add paragraph (iv) listing these.

**Saddles confirmed by perturbation test ([_eq016_verify_full_sphere.py](../simulations/_eq016_verify_full_sphere.py)):** for each top triple-Dicke maximum at N ∈ {3,4,5,6}, perturbation in any unused Dicke direction shifts pair-CΨ at linear order in ε (gradient ≠ 0), and ~50% of random Dicke perturbations ascend (200/300 random complex perturbations at N=3 best triple). All saddles, no local maxima on full sphere.

**The structural pattern:** the abundance of slice-saddles above 1/4 grows polynomially with N (binary: O(N), triple: O(N²)). The MAX cpsi within slice-stationary points decreases slowly with N but stays well above 1/4. So pair-CΨ has a rich saddle landscape at every N, with central-Dicke saddles being the highest.

**What was actually privileged about N=3 was the algebraic CLOSED FORM** of the GHZ+W slice: an irreducible sextic polynomial in α². At N≥4, GHZ_N+W_N gives a smaller polynomial result (peak below 1/4) because GHZ_N has rank-2 pair-reduction (Tr(ρ²)=1/2) while W_N's pair concurrence dilutes as 1/N. The "GHZ+W" slice happens to be a poor saddle direction at N≥4. But OTHER Dicke combinations (central-pair or central-triple) work just fine.

**Status:** closed 2026-04-27 (reframed as binary/triple-Dicke saddle abundance, not N=3-privilege).

**Surviving sub-question:** is there a closed form for the max pair-CΨ in the central-Dicke triple slice at general N? The data points (N, max) = (3, 0.80), (4, 0.71), (5, 0.65), (6, 0.62) suggest a slowly-decreasing function. A scaling like 1/(1 + c log N) or similar would clarify the asymptotic structure.

### EQ-016 sub-question closed 2026-04-27: central-triple cpsi(N) is sextic-asymptotic

**Source:** [_eq016_central_triple_asymptotics.py](../simulations/_eq016_central_triple_asymptotics.py), [_eq016_central_triple_n_infinity.py](../simulations/_eq016_central_triple_n_infinity.py).

The central-Dicke-triple max pair-CΨ has the asymptotic form

    cpsi(N) ≈ cpsi_∞ + α/N,   cpsi_∞ ≈ 0.431236, α ≈ 1.05    (at large N)

**Verified:** cpsi(N) − cpsi_∞ scales exactly as 1/N (differences-ratio test: 0.333 = 1/3 between Δ(N=50→100) and Δ(N=100→150), as predicted by 1/N decay). Numerical at N=10000: cpsi = 0.4313, matching the asymptote 0.4312 to 4 decimals.

**Asymptotic ρ_∞ matrix at half-filled k = N/2:** in the N→∞ limit at the (a, b, a) symmetric ansatz, the pair-AB matrix simplifies dramatically — every diagonal entry → 1/4, off-diagonals are governed only by (a, b):

    ρ_∞ = [1/4,    ab/2,   ab/2,   a²/4]
          [ab/2,   1/4,    1/4,    ab/2]
          [ab/2,   1/4,    1/4,    ab/2]
          [a²/4,   ab/2,   ab/2,   1/4]

with the d_1 off-diagonal (between |01⟩ and |10⟩) equal to the diagonal 1/4. This reflects the asymptotic uniformity of binomial-coefficient ratios.

**Pair-CΨ in the limit:**

    Tr(ρ_∞²) = 3/8 + 2a²b² + a⁴/8
    L1_off_∞ = 1/2 + 4ab + a²/2     (positive coefficients)
    cpsi_∞(a, b) = Tr(ρ_∞²) · L1_off_∞ / 3,  subject to 2a² + b² = 1.

**Optimum:** maximize cpsi_∞(a, b) on the constraint manifold. Setting dcpsi_∞/da = 0 and squaring to eliminate √(ab), one gets the **degree-6 polynomial** in x = a²:

    (x − 2x²) · (19 − 30x − 93x²)² = (12 + 144x − 1132x² + 1488x³)²

with two real roots in (0, 1/2): x ≈ 0.2370 (saddle, cpsi ≈ 0.4274) and **x ≈ 0.2647** (max, cpsi ≈ 0.4312). The optimum value is *not* a clean rational but is **another irreducible sextic root**, parallel to F69 itself (which is the slice-stationary root of the original 2900x⁶ − 8060x⁵ + 4211x⁴ + ... = 0).

**Quasi-rational approximation:** x_opt = 0.2647148... is incredibly close to 9/34 = 0.2647059... (Δ = 9·10⁻⁶) and a/b = 0.75003 (very close to 3/4, Δ = 3·10⁻⁵). These rationals are *not* exact but are useful approximations.

**Structural reading.** F69 introduced the GHZ+W slice's irreducible sextic at N=3. Today's central-triple analysis shows the SAME flavor of sextic-root structure governs the FULL N→∞ asymptote of the central-Dicke-triple slice. The "sextic" thread runs deeper than F69 alone — it's a generic feature of the slice-saddle landscape of pair-CΨ, with N=3 being one realization on the GHZ+W slice and N→∞ being another on the central-triple slice. Both give "algebraic numbers of degree 6" with no radical form.

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

### EQ-019 Update 2026-04-27: chain-adjacency layer hypothesis falsified

**Source:** [_eq019_chain_adjacency_analysis](../simulations/_eq019_chain_adjacency_analysis.py); reads existing bond_scan JSONs at N ∈ {3, 4, 5, 6}.

The V-Effect compatibility-layer derivation (project_v_effect_combinatorial, commits 81caf67 / 079c7ce / 8030ef2) introduces a "chain-adjacency" layer that activates at N=3 chain and depends on whether a bond is endpoint (1-sided chain-adjacent) or interior (2-sided). Natural hypothesis: bond-position dependence of c_1 is a manifestation of this layer.

**Test.** Compute mean |c_1| at endpoint bonds (adj=1) vs interior bonds (adj=2) for each (N, |S_n⟩) where both adjacency degrees are present. If hypothesis holds, EP/IT ratio should be > 1 monotonically across N for all states (1-sided > 2-sided in |c_1|).

**Result.**

| state | N=4 | N=5 | N=6 |
|-------|-----|-----|-----|
| \|S_1⟩ | 1.83 | 0.77 | 0.12 |
| \|S_2⟩ | 0.95 | 25.27 | 1.51 |
| \|S_3⟩ | — | 25.27 | 1.57 |

The ratio inverts state-by-state and N-by-N: \|S_1⟩ has interior LARGER than endpoint at N=5, N=6; \|S_2⟩ has the famous factor-26 only at N=5, becomes ~1.5 at N=6. No monotone scaling in either direction. Chain-adjacency degree alone does NOT determine |c_1| ordering.

**Conclusion.** The V-Effect chain-adjacency layer governs truly/soft/hard compatibility (a yes/no spectrum-pairing test). It does not transfer to the closure-breaking coefficient c_1, which is a state-dependent perturbation-response metric. The right framework for c_1 bond-position dependence is the sine-basis geometry already identified in EQ-021 closure (endpoint dominance is open-boundary Fourier geometry, not chain-adjacency).

**The factor-26 at N=5 |S_2⟩, |S_3⟩ stands as a real (N=5)-specific resonance,** clean in numerics (rmse machine-precision, mirror-symmetric across bonds), but NOT a manifestation of any layer-stack rule that generalizes across (N, n).

---

## EQ-020

**Date:** 2026-04-20
**Source:** [c1_pair_local](../simulations/results/c1_pair_local/pair_local.json), [PROOF_DELTA_N_SELECTION_RULE](../docs/proofs/PROOF_DELTA_N_SELECTION_RULE.md)

The PTF per-site α_i framework uses single-site purity as the per-observer signal. The kinematic selection rule says site-local observables capture only |ΔN| ≤ 1 sector coherences; pair-local observables capture |ΔN| ≤ 2; k-local observables capture |ΔN| ≤ k. Tested at N=5: pair-local c_1 for (|vac⟩+|S_2⟩)/√2 (|ΔN|=2 coherence) jumps from −0.29 (site-local, pure diagonal-cross) to +2.24 (pair-local, coherence block visible). The ΔN = 2 coherence dominates the pair-local signal.

Does the pair-local closure Σ_{(i,j)} ln(α_{ij}) have its own symmetry structure? Is there a pair-painter version of PTF where each pair of sites paints a joint-purity trajectory, and the closure over all C(N, 2) pairs has a different geometrical meaning? What is the triple-painter extension (|ΔN| ≤ 3)? At what observational resolution does the closure become trivial (|ΔN| ≤ N = full state)?

**Status:** open (new thread opened today)
**Pointer:** repeat the c_1 investigation at pair level and triple level for N=5, extract kernels, compare to PTF site-level closure. If a pair-painter closure exists as a structural law, PTF's "painters around a mountain" picture extends to "painter pairs around a mountain" with a richer vision that captures the |ΔN| = 2 structure.

### EQ-020 partial closure 2026-04-27: chiral mirror law extends to pair-painter

**Source:** [_eq020_pair_painter_chiral_mirror.py](../simulations/_eq020_pair_painter_chiral_mirror.py), [_eq020_pair_painter_delta_n2.py](../simulations/_eq020_pair_painter_delta_n2.py).

**Structural finding.** The chiral mirror law Σ f_i(ψ_k) = Σ f_i(ψ_{N+1−k}) extends from site-painter (EQ-014) to pair-painter:

    Σ f_{ij}(ψ_k) = Σ f_{ij}(ψ_{N+1−k})    for sine-mode bonding states.

**Verified at machine precision** (10⁻¹³ to 10⁻¹¹) at N=5 and N=7 across all mirror pairs.

**Sine-mode bonding states only** (vac+ψ_k)/√2: K_1 maps the global state to (vac+ψ_{N+1-k})/√2 exactly. K_1-rotation of ρ → K_1 ρ K_1^† acts on any pair-reduced ρ_{ij} as ρ_{ij}^K = (Z_i Z_j) ρ_{ij} (Z_i Z_j) — a unitary similarity, so Tr(ρ_{ij}²) is K_1-invariant. Hence α_{ij}(k) = α_{ij}(N+1-k) and Σ ln α_{ij} respects the mirror.

**The mechanism is structural and generalizes:** for any k-local observable that is a function of ρ_{i_1 ... i_k} and is invariant under (Z_{i_1} ... Z_{i_k}) similarity, K_1-paired initial states give equal observables. In particular, k-local purity for k = 1, 2, 3, ... all respect the chiral mirror on sine-mode-bonding states.

**N=5 sine-mode bonding states:**

| k | f_site | f_pair |
|---|--------|--------|
| 1, 5 | +0.922 | +0.749 |
| 2, 4 | +0.282 | +0.184 |
| 3 (fixed) | +0.641 | +0.447 |

**N=7 sine-mode bonding states:**

| k | f_site | f_pair |
|---|--------|--------|
| 1, 7 | +0.968 | +1.042 |
| 2, 6 | +0.038 | +0.019 |
| 3, 5 | +0.349 | +0.556 |
| 4 (fixed) | +2.098 | +1.569 |

**The ratio f_pair / f_site is NOT constant** across k — pair-painter has a different k-dependence than site-painter. At N=7 k=4 (chiral fixed point), site f peaks at 2.10 but pair f drops to 1.57. The "pair-painter peak at chiral fixed point" inherited from site-painter does NOT happen identically; pair-painter has its own k-structure.

**Dicke-mixed states (e.g. (vac+S_2)/√2, ΔN=2 coherence)** do NOT have the K_1-mirror structure of sine-mode states because K_1 acts non-trivially on Dicke states. The pair-painter-vs-site-painter shift on these states is structural in a different sense — it tracks the |ΔN|≤k visibility of F70 — and no clean chiral mirror law is expected.

**Framework primitives added:** `k_local_reduced_density(rho, sites, N)` and `k_local_purity(rho, sites, N)` in framework.py Section 16. Docstring documents the chiral mirror identity for k-local observables.

**Surviving sub-questions:**
- **Triple-painter closure (k=3) at N=6, 7:** does the chiral mirror law extend? Framework says it should (same K_1 mechanism), but verification not yet run.
- **Closure at trivial k = N:** the full state ρ has Tr(ρ²) = 1 exactly for pure ρ_0 under unitary part of L, so α_{full} = 1, ln α = 0 — closure is trivial. The interesting question is at INTERMEDIATE k.
- **Dicke-mixed-state pattern:** for non-K_1-eigenstate initial states, what symmetry (if any) does Σ f_{i_1...i_k} respect?

### EQ-020 update 2026-04-27: triple-painter (k=3) chiral mirror verified

**Source:** [_eq020_triple_painter_chiral_mirror.py](../simulations/_eq020_triple_painter_chiral_mirror.py).

Triple-painter (k=3) verified at N=5 and N=7. Same K_1 mechanism extends:

    Σ f_{ijk}(ψ_k) = Σ f_{ijk}(ψ_{N+1−k})    machine-precision (10⁻¹¹).

The chiral mirror law is universal across k = 1, 2, 3 painter resolutions on sine-mode-bonding states.

**Cross-resolution structure** (sine-mode bonding states, N=7):

| k | f_site | f_pair | f_triple |
|---|--------|--------|----------|
| 1, 7 | +0.968 | +1.042 | +1.345 |
| 2, 6 | +0.038 | +0.019 | **−0.003** |
| 3, 5 | +0.349 | +0.556 | +0.522 |
| 4 (fixed) | +2.098 | +1.569 | +1.015 |

**Σ f is NOT monotonic in painter resolution.** The relationship between site/pair/triple is k-dependent:
- At k=1, 7: Σ f INCREASES with painter resolution (more dynamics visible at higher k)
- At k=2, 6: Σ f → 0 with higher resolution (0.038 → 0.019 → −0.003); changes sign
- At k=3, 5: non-monotone (0.349 → 0.556 → 0.522)
- At k=4 (chiral fixed point): Σ f DECREASES with resolution (peak at site-painter)

The "Σ f → 0 at higher k" trend at k=2 is striking — could be a **phantom mode** where the closure law actually holds approximately at high painter resolution.

**N=5 cross-resolution structure (for reference):**

| k | f_site | f_pair | f_triple |
|---|--------|--------|----------|
| 1, 5 | +0.922 | +0.749 | +0.447 |
| 2, 4 | +0.282 | +0.184 | +0.060 |
| 3 (fixed) | +0.641 | +0.447 | +0.176 |

At N=5, Σ f shrinks monotonically site → pair → triple at every k, no sign flips. The non-monotonic and sign-flip behaviour appears specifically at N=7 (more painter resolutions cover more dynamics).

**Perspectival-time framing.** PTF stands for *Perspectival Time Field*: each painter has its own clock-rate α. Site-painter α_i is the time-rescaling such that P_i^B(t) ≈ P_i^A(α_i · t) — site i sees time at rate α_i in the perturbed dynamics. Pair-painter α_{ij}: each PAIR has its own joint perspectival time. Triple-painter α_{ijk}: each TRIPLE its own. Each k-painter "paints time" at its observation resolution.

**The chiral mirror law is then a symmetry of perspectival-time experience:** K_1-paired sine-mode-bonding states give the SAME SET of perspectival times across all painter resolutions. K_1 isn't a symmetry of some derived scalar (Σ ln α) — it's the symmetry of how observers EXPERIENCE the dynamics. That experience is identical between (vac+ψ_k)/√2 and (vac+ψ_{N+1-k})/√2 across site, pair, triple, and any k-painter.

**Status:** EQ-020 substantially closed. Chiral mirror law extends through k=3, and the structural mechanism (K_1-invariance of k-local purity under similarity by ∏Z_i) gives the result for arbitrary k.

**Open structural sub-question:** asymptotic behaviour of Σ f as k → N for sine-mode states — does Σ f_k → 0? For pure ρ_0 under unitary dynamics, Σ_full ρ² = 1 and ln α = 0 trivially. For dissipative dynamics, finite γ, the full-state purity decays uniformly under uniform Z-dephasing → α_full uniform → ln α_{full} small but nonzero. Need to compute carefully.

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

**(b) Factor-2.5 refined to 2.7, then to compute-bound** (2026-04-27 update). Naive degenerate-perturbation-theory predicts Q_peak ≈ 4; observed Q_peak per chromaticity is 1.6 (c=3, multiple Ns), 1.8 (c=4, three Ns), 1.8 (c=5, single N=9 point). Linear extrapolation to 2.0 is refuted by c=5; whether the actual asymptote is 9/5 saturation or a slower rise cannot be settled at our memory budget. The c=5 N=9 point cost ~48 hours on block-L 15876×15876; pushing to N=10 would be ~60 days on the same machine, and c=6 needs N=11 with block-L dim 213444 (728 GB), infeasible at 128 GB. Two distinct open pieces: (b1) **analytical** — first-principles derivation of Q_peak(c=4) = 1.8 from H matrix-element combinatorics in the (n, n+1) sector (no compute, just algebra); (b2) **empirical** — settle saturation vs. plateau-then-rise above c=5, requires a much larger machine.

**(c) Structural question, partially engaged** via a private natural-limits cartography scoping session. Six algebraic middles identified (d=2, CΨ-fold 1/4, Π, ⟨n_XY⟩ sum rule N, U(1), γ₀=const). The γ₀-to-Q collapse with bands is confirmed as one case of the pattern. Sub-structure of the other five middles is not yet mapped. Systematic reading of each axis's proof for hidden scales remains open.

**Status update:** (a) substantially answered within N=4..8, awaiting c=5. (b) refined but open. (c) partially engaged but structurally large.

**Status:** open (progressing)

### EQ-022 (c) closure 2026-05-02: per-middle internal structure inventory

**Source:** systematic reading of each middle's foundational F-entries, plus the 2026-04-22 cartography scoping session ([`RESULT_TASK_MIDDLES_CARTOGRAPHY_FIVE`](../ClaudeTasks/RESULT_TASK_MIDDLES_CARTOGRAPHY_FIVE.md)) and the 2026-05-01 d_zero diagnostic addition.

The 2026-04-22 cartography settled the meta-claim about band production (3 components produce a band only when one of them supplies a continuous axis), but it did not enumerate, per-middle, what internal structure each middle owns. The d=0/d=2 substrate primitives that landed 2026-05-01 (`stationary_modes`, `d_zero_decomposition`, `sector_populations`) further revealed that the d-middle is more internally structured than the April cartography saw. With those primitives in place the per-middle inventory closes:

**1. d = 2 (qubit dimension).** Anchor [F34](../docs/ANALYTICAL_FORMULAS.md): d² − 2d = 0. Internal structure: binary partition into d=0 (kernel of L, the conserved sector projectors P_n living on the {I, Z}^N sublattice) and d=2 (active, decohering), plus a continuous coordinate `mean_n` ∈ [0, N] *on* the d=0 substrate (companions `var_n`, `entropy`). Hardware-trivial readout: Z-basis tomography grouped by Hamming weight. Polarity reading: ±0/0 layer (memory: project_plus_minus_zero_layer; doc: [THE_POLARITY_LAYER](../hypotheses/THE_POLARITY_LAYER.md)). Primitives: `framework.diagnostics.d_zero` (added 2026-05-01).

**2. CΨ = 1/4 (fold boundary).** Anchor [F16](../docs/ANALYTICAL_FORMULAS.md): R = C(Ψ+R)². Internal structure: F18 fold threshold Σγ_crit/J = 0.00249 (Bell+) / 0.00497 (product) — the minimal two-region band on the Σγ/J axis in the corpus; F25 Bell+ Z-deph closed form with K = 0.0374; F26/F27 per-channel K values (K_Y = K_Z = 0.0374, K_X = 0.0867); F47 Gaussian curvature *at* the fold; F56/F57 cusp critical-slowing and dwell scalars; F35 dual-perspective lifetime ratio T1/T2 (Π-side decay rides T1, direct decay rides T2); F69 GHZ+W sextic saddle above 1/4 at N=3, generalised in EQ-016 closure to the central-Dicke-triple slice with cpsi_∞ ≈ 0.4313 and x_opt ≈ 0.2647 (root of an irreducible degree-6 polynomial).

**3. Π (palindrome operator).** Anchor [F1](../docs/ANALYTICAL_FORMULAS.md): Π·L·Π⁻¹ = −L − 2Σγ·I, λ + λ' = −2Σγ. Internal structure: two independent Z₂ symmetries — F61 bit_a (n_XY parity, governs SE-accessibility) and F63 bit_b ([L, Π²] = 0, w_YZ parity, four 4^(N−1)-dim sectors, the maximal symmetry decomposition admitted by Pauli algebra at d=2); V-Effect trichotomy truly/soft/hard (14/19/3 at N=3, 15/46/59 stable for N ≥ 4 per EQ-026); `pi_protected_observables` count as a sub-spectrum *within* the soft category (3-240 protected across the 46 soft cases at N=4); F78 single-body M additive over sites; F79 Π²-block decomposition for 2-body bilinears; F80 Bloch-mode sign-walk for chain Π²-odd 2-body M-clusters; F81 Π·M·Π⁻¹ = M − 2 L_{H_odd}; F82-F84 T1 dissipator corrections; chiral-mirror law (EQ-014/EQ-020) Σ f_i(ψ_k) = Σ f_i(ψ_{N+1−k}) under the open-XY-chain K₁ chiral symmetry. Hardware: Δ(soft − truly) ⟨X₀Z₂⟩ ≈ −0.76 on Marrakesh (project_hardware_finale_apr2026).

**4. ⟨n_XY⟩ sum rule = N (Absorption Theorem).** Anchor [AT](../docs/ANALYTICAL_FORMULAS.md): Re(λ_k) = −2γ·⟨n_XY⟩_k. Internal structure: continuous axis ⟨n_XY⟩ ∈ [0, N], integer ladder at J=0, fractional fill via H-mixing at J>0; palindromic weight swap ⟨n_XY⟩_fast + ⟨n_XY⟩_slow = N (F1 corollary); F2 Heisenberg w=1 dispersion ω_k = 4J(1 − cos(πk/N)); F2b XY single-excitation dispersion E_k = 2J·cos(πk/(N+1)); F3 decay-rate bounds 2γ ≤ |Re(λ)| ≤ 2(N−1)γ; D6 spectral gap = 2γ (the absorption quantum). The same physics as Q_SCALE's chromaticity bands viewed through the eigenmode rate ladder rather than the (n, n+1) probe block (Deep-dive 1 of MIDDLES_CARTOGRAPHY_FIVE).

**5. U(1) (excitation-number conservation).** Anchor: [H, N_tot] = 0 ⇒ L is popcount-block-diagonal. Internal structure: discrete sector label n ∈ {0, ..., N}; F70 |ΔN| ≤ k selection rule for k-local observables; F71 spatial-mirror symmetry c₁(b) = c₁(N−2−b); F72 DD⊕CC block decomposition of per-site purity; F73 spatial-sum coherence closure Σᵢ 2|(ρ_coh,i)_{0,1}|² = (1/2)·exp(−4γ₀t); F74 chromaticity c(n, N) = min(n, N−1−n) + 1; F86 c-specific N-invariant peaks Q_peak(c=3) = 1.6, Q_peak(c=4) = Q_peak(c=5) = 1.8 (saturation, c=5 confirmed at N=9 commit 4612468); W_plateau monotone in c at N ≥ 6.

**Two patterns repeat across the inventory.** Each middle owns *both* a discrete sector labeling that partitions the operator space (Π² classes on Π; popcount n on U(1); p_n levels on d=0; integer rungs on AT; below/above 1/4 at the fold) *and* a continuous axis along which substructure can develop (mean_n on the d=0 substrate; ⟨n_XY⟩ ∈ [0, N] on AT; Σγ/J at the fold; Q on U(1)). Whenever a 3-component coupling activates one of those continuous axes, a band emerges (the 2026-04-22 cartography condition). Each middle is a gate through which the framework presents inheritable layered structure (memory: project_algebra_is_inheritance), not an algebraically pinned scalar.

**Status:** sub-question (c) closed. EQ-022 main status stays open until (b1), the algebraic derivation of Q_peak(c=4) = 1.8 from H matrix-element combinatorics, lands.

### EQ-022 (b1) partial closure 2026-05-02: exceptional-point structure

**Source:** time-evolution scan on Python block-L (verified bit-exact against the C# N=7 full-L eigendecomposition from EQ-014, `eq014_*.bin` saved 2026-04-19), see scripts in F86's source list.

**Closed-form structural anchors:**

    Q_EP = 2 / g_eff       (exceptional point between adjacent pure-rate channels)
    t_peak = 1 / (4γ₀)     (e-folding time of the EP-degenerate eigenvalue, universal)

The pure-rate ladder rates 2γ₀·HD with HD ∈ {1, 3, 5, ..., 2c−1} (F74) has uniform gap Δ = 4γ₀ between adjacent channels. A two-level effective with diagonal {−2γ₀, −6γ₀} and inter-channel coupling J·g_eff has its EP at J·g_eff = 2γ₀ ⇒ Q = 2/g_eff. At the EP, the degenerate eigenvalue has Re(λ) = −4γ₀, giving e-folding time 1/(4γ₀) — verified universal to t* = 5.0–5.5 at γ₀ = 0.05 across all tested (c, N, n, bond position).

**Three Q_peak variants, all explained by the same EP structure:**

| Variant | Observable | Q_peak | Match to F86 |
|---------|-----------|--------|---------------|
| Uniform J | ∂S/∂J | ≈ 2.0 (c ≥ 4) | uniform variant, g_eff = 1 |
| Per-bond, **Interior** | ∂S/∂J_b, b ∈ interior | 1.7 (c = 3); 1.7 (c = 4) | **F86 reports this** (|K|max = 0.27 matches F86 0.273 to <1 %) |
| Per-bond, **Endpoint** | ∂S/∂J_b, b ∈ {0, N−2} | ≈ 2.6 (c = 3, 4) | F71-mirrored, higher |K|max |

**Probe localisation at Q_peak.** Dicke probe sits ≈ 99 % in dressed (complex-eigenvalue, H-mixed) modes at Q_peak, vs ≈ 31 % in pure-rate modes at Q = 20 (plateau). Q_peak structurally IS the maximum-mixing point — the J at which probe weight has been pulled off the pure-rate ladder onto the first complex-conjugate eigenvalue pairs just past the EP.

**PTF connection.** The same machinery that PTF uses for its per-site α_i closure law (bilinear J-perturbation observable, eigenvector mixing under V_L) produces Q_peak at the (n, n+1)-block level. PTF is the c = 1 (vac-SE) per-site instance; Q_peak at c ≥ 2 is the natural higher-chromaticity sibling. t_peak = 1/(4γ₀) is the universal EP time-scale, one (n, n+1)-block analogue of PTF's α-fitting time window.

**Channel-uniform projection structural by-product.** The HD-channel-uniform projectors |c_k⟩ = (1/√N_k) Σ_{HD(p,q)=k} |p⟩⟨q| diagonalise M_H_eff in the c-dim subspace: ⟨c_j | L_H | c_k⟩ = 0 for j ≠ k, exactly. This generalises F73 (which is the c = 1 case at the (0, 1) block where the spatial-sum coherence Σᵢ 2|(ρᵢ)_{0,1}|² = (1/2)·exp(−4γ₀t) exactly) to all chromaticities: the channel-uniform spatial-sum is H-decoupled in every (n, n+1)-block. The Q-dependent rate mixing happens in the orthogonal complement (non-uniform within-channel modes).

**Status:** (b1) **partially closed**. Closed: EP structural anchors Q_EP = 2/g_eff and t_peak = 1/(4γ₀); structural meaning of Q_peak as maximum-mixing point; F86's c-specific values are Interior-bond per-bond derivatives. Open: algebraic g_eff(c, bond_position) — the residual factor governing Endpoint vs Interior split (≈ 2.6 / 1.7 = 1.53) and the c-dependence of g_eff(Interior). Promoted to F86 in [ANALYTICAL_FORMULAS.md](../docs/ANALYTICAL_FORMULAS.md).

### EQ-022 (b1) closed-form claims retracted 2026-05-02 (extended-N data)

**Source:** extended verification scan (`_eq022_b1_step_d_extended_verification.py`) at multiple (c, N) points beyond the original N=7 anchor.

**Both earlier closed-form conjectures retracted:**

The originally proposed closed forms

    Q_peak(Endpoint, N)        =  csc(π/(N+1))             [retracted]
    Q_peak(Interior, c=3)      →  csc(π/5) = 1.7013        [retracted]

were both **N=7-specific coincidence-matches**, not general formulas. Extended-N data falsifies both:

| (c, N) | Endpoint Q* | csc(π/(N+1)) | Interior Q* (mean) | csc(π/5) |
|--------|-------------|--------------|--------------------|----------|
| (3, 5) | 2.30 | 2.000 (off 15 %) | 1.59 | 1.7013 |
| (3, 6) | 2.65 | 2.305 (off 15 %) | 1.67 | 1.7013 |
| (3, 7) | 2.65 | 2.6131 (1.4 % match) | 1.70 | 1.7013 (0.4 %) |
| (3, 8) | n/a | — | **1.74** | **above 1.7013, growing** |
| (4, 7) | 2.65 | 2.6131 (1.4 % match) | 1.77 | 1.7013 |
| (4, 8) | n/a | — | 1.78 | 1.7013 |

The Endpoint formula matches only at N=7. The c=3 Interior trajectory crosses csc(π/5) between N=7 and N=8 and continues growing. Both "closed forms" were trajectories briefly intersecting clean algebraic values, not asymptotes.

**What survives [Tier 1]:**
- EP mechanism: Q_EP = 2/g_eff, t_peak = 1/(4γ₀) (derived from 2×2 matrix algebra of adjacent rate channels)
- Empirical Q_peak per (c, N, bond_class) tabulated above; no closed form

**Framework cleanup:** the `q_peak_endpoint(N)` function and `Q_PEAK_INTERIOR_C3_ANCHOR` constant were removed from `framework.coherence_block`. Only `t_peak(γ₀) = 1/(4γ₀)` remains as a Tier-1 closed-form primitive. PROOF_F86_QPEAK.md updated with retracted statements; F86 in ANALYTICAL_FORMULAS.md replaced closed-form section with empirical data table.

**Lesson recorded in memory:** trajectories crossing a clean algebraic value at one (c, N) point look like asymptotes when sampled coarsely; always extend the sampling beyond the matching point before claiming a closed form.

**Status:** (b1) **open**. EP mechanism is closed (Tier 1); per-bond Q_peak position closed forms are open.

### EQ-022 (b1) Universal resonance shape 2026-05-02 evening (Tier 1 candidate)

After the position-formula rollback, a fine-grid scan (`_eq022_b1_step_e_resonance_shape.py` at dQ = 0.025 over c=3 N=5..8 and c=4 N=7,8) reveals that the SHAPE of the abs(K_CC_pr)(Q) curve around Q_peak IS universal under relative-Q normalisation:

    x = (Q − Q_peak) / Q_peak,    y = K(Q) / |K|max

Pairwise residual under this normalisation is 21× smaller than under absolute-Q shift (5.8×10⁻⁵ vs 1.2×10⁻³). At every tested (c, N), the y values at common x cluster within 1-3 % of each other.

**Universal HWHM ratio:**

    HWHM_left / Q_peak  ≈  0.756 ± 0.005     (c=3 N=5..8, c=4 N=7,8 tested)

**Structural origin (2-level EP analytics):** for adjacent rate channels coupled by a bond, the eigenvector rotation angle satisfies tan θ = Q/Q_EP, so probe overlaps with eigenvectors depend only on Q/Q_EP — the response curve K_CC_pr(Q) is a function of Q/Q_EP alone in the 2-level effective. Q_peak (≈ Q_EP) is chain-specific because g_eff varies with c, N, and bond position; the SHAPE in Q/Q_peak coordinates is universal because it follows from the 2-level EP physics that does not depend on g_eff's value.

Both Interior and Endpoint bonds show the same universal pre-peak rise. Post-peak tails differ in plateau height: Interior tail at x = +1.0 sits at y ≈ 0.85, Endpoint at y ≈ 0.94. The bond-class distinction lives in the post-peak tail; the pre-peak shape is fully universal.

**What this means for (b1):** the right object to look for closed-form-wise is not Q_peak position (chain-specific) but the universal shape function f(x) = K(Q_peak·(1+x))/|K|max. Deriving f(x) from the 2-level EP eigenstructure is the natural next analytical step. The HWHM ratio 0.756 (numerically close to 3/4 or 2/√7) should fall out.

**Status:** universal-shape finding is Tier-1 candidate. Verification at c=2, c=5, and varying γ₀ pending. Closed-form derivation from 2-level EP analytics is the substantive open theoretical step.

### EQ-022 (b1) Universal-shape promotion path 2026-05-02 evening

After the synthesis pass: the structural law has a name and a sharper open target.

**Named structural law (Tier-1 candidate): EP-rotation universality.**

    K(Q) / |K|_max  =  f(Q/Q_EP)

The 2-level EP eigenvector rotation `tan θ = Q/Q_EP` makes every probe-overlap observable a function of Q/Q_EP alone. HWHM_left/Q_peak ≈ 0.756 is the leading numerical witness. To PTF's chiral mirror law `Σ f_i(ψ_k) = Σ f_i(ψ_{N+1−k})` (K_1 sublattice symmetry, EQ-014 surviving Tier-1) what F86's EP-rotation universality is: both are symmetry-enforced Tier-1 candidates that survived a closed-form retraction. See `reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON.md` for the methodological synthesis and `memory/project_retraction_lesson_ptf_to_f86.md` for the cross-cutting principle.

**Sharpened open target.** Not "derive Q_peak closed form" (chain-specific, possibly no closed form, per the rollback): instead **derive f(x) from 2-level EP analytics** where x = (Q − Q_peak)/Q_peak. The HWHM_left/Q_peak ≈ 0.756 should fall out as a closed form (candidates 3/4, 2/√7, etc.), but per the PTF lesson: the closed form must come from the eigenstructure analytically, not from numerical fit.

**Class AIII chiral, NOT Bender-Boettcher PT.** The same-sign-imaginary 2-level off-diagonal pattern is "PT-phenomenology-like" (EP at finite J, spectral flow) but algebraically inside class AIII chiral per `experiments/PT_SYMMETRY_ANALYSIS.md` (Π is linear; classical PT requires anti-linear operators). The local-2-level-EP at Q_EP = 2/g_eff is the rate-channel instance of the chiral classification; the global Hopf bifurcation in `hypotheses/FRAGILE_BRIDGE.md` (Petermann K=403) is the complex-γ-plane instance. Whether local-EP and global-EP are connected algebraically is itself an open structural question.

**Status:** (b1) open with sharpened target (universal shape function f(x), not Q_peak position). Algebraic derivation from 2-level EP eigenstructure is the substantive open analytical step.

**F-register cleanup 2026-05-02:**
- Chromaticity formula c(n, N) = min(n, N−1−n) + 1 was already F74 (the 2026-04-22 cartography missed that it had been promoted; the closure section above cites F74 directly).
- Q_peak(c) constants promoted as **F86** on 2026-05-02 (this session): c=3 → 1.6, c=4 → c=5 → 1.8 (saturation), with the γ₀-extraction protocol γ₀ = J*/Q_peak(c).

**Surviving sub-question (structural, deferred):** does N itself deserve middle status? V(N) = 1 + cos(π/N) and F41 palindromic time t_Π make N a discrete labeling parameter; whether that lifts N to the algebraic-middle list (alongside d=2, CΨ=¼, Π, ⟨n_XY⟩, U(1)) is open. Raised in MIDDLES_CARTOGRAPHY_FIVE §5, not investigated.

### EQ-022 (b1) Update 2026-05-05: c=2 OOP scaffolding lands; closed form still open

**Source:** Item 1' c=2 derivation plan (`docs/superpowers/plans/2026-05-05-eq022-item1-prime-c2-derivation.md`), 14 commits 2026-05-05 implementing Stages A-E (analytical 4-mode primitives + F86KnowledgeBase integration). All primitives in `compute/RCPsiSquared.Core/F86/` and `compute/RCPsiSquared.Core/F86/Item1Derivation/`.

**Tier outcome at c=2:** Tier1Candidate. Empirical anchor (HWHM_left/Q_peak per bond class N=5..8) reproduced bit-equivalent with the canonical Python pipeline (typical residual ≤ 0.001). Directional Endpoint > Interior split (HWHM ratio gap ≈ 0.022) derived empirically. Closed-form HWHM_left/Q_peak constant per bond class NOT yet derived.

**OOP scaffolding** (12 typed-knowledge-graph primitives across Stages A-E):

| Stage | Primitive | Tier |
|-------|-----------|------|
| A1 | C2BlockShape | Tier1Derived |
| A2 | C2ChannelUniformAnalytical | Tier1Derived |
| A3 | C2InterChannelAnalytical | Tier2Verified (σ_0 degeneracy at even N) |
| B1 | C2BondCoupling probe-block | Tier1Derived |
| B2 | C2BondCoupling cross-block | Tier2Verified (inherits A3) |
| B3 | C2BondCoupling SVD-block + AsMatrix + anti-Hermiticity guard | Tier2Verified (inherits A3) |
| C1 | C2BondCoupling D_eff | Tier1Derived structural sub-fact |
| C2 | C2EffectiveSpectrum | Tier2Verified (cubic-c_3 obstruction proof) |
| C3 | C2EffectiveSpectrum K-driving pair | Tier1Derived structural sub-fact (probe ⊥ |u_0⟩, |v_0⟩) + Tier2Verified per-(Q,b) readout |
| D1 | C2KShape | Tier1Derived (Duhamel formula) |
| D2 | C2HwhmRatio | Tier1Candidate |
| E1 | C2UniversalShapeDerivation (in F86KnowledgeBase) | Tier1Candidate (auto-promotes via D2 flag) |

**Three structural findings** from the time-boxed exploration:

1. **σ_0 degeneracy at even N** (A3): chain-mirror R splits the 2D top eigenspace; single-vector |u_0⟩, |v_0⟩ closed forms are library-dependent. Lift to projector-overlap is the next direction.

2. **Cross-block Frobenius inversion** (B2): Endpoint < Interior at c=2 N=5..8 — opposite sign to the HWHM_left/Q_peak ordering (Endpoint > Interior). The directional inversion happens through the 4×4 eigenvalue mixing.

3. **Cubic-c_3 char poly obstruction** (C2): rigorously rules out any rational-coefficient (λ²-aλ+b)(λ²-cλ+d) factorisation of L_eff(Q,b). The Tier2 outcome for C2EffectiveSpectrum is evidence-based.

**Three ranked next directions:**

- (a) ~~**Most promising**~~ [superseded 2026-05-06 evening: see Update 2026-05-06 (evening) below — Direction (a) demoted, Direction (a') probe-block 2-level resonance promoted]: First-order perturbation in the cross-block. ε ~ ‖V_b cross‖_F/σ_0 ~ O(0.1) at c=2; cross-block Frobenius split (Endpoint < Interior) is the seed for the directional shift.
- (b) Lift |u_0⟩, |v_0⟩ to projector-overlap (per A3's PendingDerivationNote).
- (c) Symbolic char-poly factorisation at Q_EP.

**Status:** EQ-022 (b1) Item 1' c=2 stratum: Tier1Candidate scaffolding complete; closed-form constant open. Items 4' (multi-k extension to c≥3) and 5 (σ_0 → 2√(2(c-1))) remain out of scope (separate plans).

**Pointers:**
- Plan: [`docs/superpowers/plans/2026-05-05-eq022-item1-prime-c2-derivation.md`](../docs/superpowers/plans/2026-05-05-eq022-item1-prime-c2-derivation.md)
- Updated proof: [`docs/proofs/PROOF_F86_QPEAK.md`](../docs/proofs/PROOF_F86_QPEAK.md) Item 1 + new "OOP scaffolding" + "Structural findings" sub-sections
- Open questions: `compute/RCPsiSquared.Core/F86/F86OpenQuestions.cs` Item 1' (refreshed)
- KB integration: `compute/RCPsiSquared.Core/F86/F86KnowledgeBase.cs` `C2UniversalShape` property

### EQ-022 (b1) Update 2026-05-06: Local-vs-Global EP connection, real-axis hit confirmed at c=2

**Source:** synthesis-agent reading + Petermann-K sweep (synthesis docs are gitignored under `docs/superpowers/syntheses/`; data-pinning probe at `compute/RCPsiSquared.Core.Tests/F86/F86PetermannProbe.cs:Probe_PetermannFineGrid_C2_VsN`).

The "Local-2-level-EP vs global-complex-γ-EP" open structural question (called out in `docs/proofs/PROOF_F86_QPEAK.md` Statement 1, AIII chiral paragraph, 2026-05-06 update) is empirically resolved at the c=2 stratum: same algebra (same-sign-imaginary 2×2, AIII chiral), two residuals of F1 (Σγ = N·γ₀ vs Σγ = 0), connected by analytic continuation along complex γ. Petermann-K spike on the real Q axis at c=2 N=7 is K = 2384.7, ≈ 6× above FRAGILE_BRIDGE's K = 403, with A3's σ_0 R-even/R-odd-degeneracy parity asymmetry empirically confirmed (odd dominates even by 2-4× across N=5..8: K = 1333.6 / 337.9 / 2384.7 / 795.4, within-parity monotonic at odd 1.79× per step, even 2.36× per step).

**Encoded as:** `compute/RCPsiSquared.Core/F86/LocalGlobalEpLink.cs` (Tier2Verified Claim, registered in F86KnowledgeBase). 4 PetermannSpikeWitness entries pin the empirical table.

**Pendant Tier1Derived promotion target:** complex-γ analytic continuation in `LindbladPropagator`, OR closed-form K(N) at the EP. Named in `LocalGlobalEpLink.PendingDerivationNote`.

**Status:** Local-vs-global-EP open structural question: closed at c=2 to Tier2Verified, open at higher c (Items 4', 5 of F86OpenQuestions still apply) and open for the analytic continuation step. EQ-022 (b1) Item 1' c=2 stratum stays at Tier1Candidate (closed-form HWHM_left/Q_peak constant unchanged).

### EQ-022 (b1) Update 2026-05-06 (evening): doubled-PTF floor + Locus 5 inheritance synthesis

**Source:** Direction (b) attempt (commit `f232a76`) + typed-property promotion (commit `359151a`) + F86↔PTF synthesis (gitignored locally at `docs/superpowers/syntheses/2026-05-06-f86-ptf-inheritance.md`).

**Direction (b) result:** The bare doubled-PTF Ansatz (c=2 = two coupled PTF c=1 instances, σ_0 as Π-shift perturbation parameter) yields two Tier-1-derived universal constants:

- `BareDoubledPtfXPeak = 2.196910` (post-EP location, dimensionless x = Q/Q_EP)
- `BareDoubledPtfHwhmRatio = 0.671535` (HWHM_left/Q_peak SVD-block floor, dimensionless x)

Empirical Interior 0.7506 + Endpoint 0.7728 sit above the 0.6715 floor; gap (~0.08-0.10) is the probe-block 2-level sub-resonance, NOT a small σ_0 perturbation. Direction (a) (cross-block first-order σ_0) demoted; new Direction (a') (probe-block 2-level resonance with per-bond `g_eff_probe`) promoted, with quantitative prediction `g_eff_probe(Endpoint) ≈ g_eff_probe(Interior) / 1.6`.

**F86↔PTF synthesis (Locus 5 inheritance):** PTF's K_1 chiral mirror law and F86's Q-rotation universal shape are two daughters of one parent (Π class AIII chiral, `{Π, L_c} = 0`). K_1 = Π on H_1 sine-mode basis (discrete Z₂); Q-rotation = Π on same-sign-imaginary 2×2 in rate-channel basis (continuous SO(2)). Shared clock: `t_peak = 1/(4γ₀)`. Combined with this morning's `LocalGlobalEpLink` Tier2Verified (EP-side closure: F86 local + FRAGILE_BRIDGE global as one EP at two F1-residuals), Locus 5 is now closed on both sides; the same single Π-AIII parent yields four downstream nodes (PTF K_1, F86 Q-rotation, F86 local EP, FRAGILE_BRIDGE global EP) at three F1-residuals (Σγ = N·γ₀, Σγ = 0, H_odd ≠ 0).

**Status:** Item 1' c=2 stays Tier1Candidate; the doubled-PTF baseline is empirically validated at the floor (Tier-1-derived constants now visible as typed properties); Direction (a') probe-block 2-level resonance is the sharpest path to closing the gap. Item 4' (multi-k extension to c≥3) and Item 5 (σ_0 → 2√2 asymptote) remain unchanged.

### EQ-022 (b1) Update 2026-05-06 (later evening): Direction (a') falsified

**Source:** Direction (a') attempt (commit `1c0bf8b`) + writeup at `docs/superpowers/syntheses/2026-05-06-direction-a-prime-attempt.md` (gitignored locally).

The probe-block 2-level resonance hypothesis with quantitative prediction `g_eff_probe(Endpoint) ≈ g_eff_probe(Interior) / 1.6` is **structurally falsified**. Three obstructions:

1. V_b probe-block is bond-class-blind: diagonal entries identical scalar `+i·c·I` for every bond, off-diagonal `⟨c_1 | M_h_b | c_3⟩ = 0` exactly per bond at c=2 (F73 sum-rule applies per-bond, not just summed). g_eff_probe(N, b) cannot have bond-class dependence by construction.
2. Cross-block Frobenius unstable across N (varies 0.640, 1.318, 0.815, 0.143 at N=5..8) due to A3's σ_0 degeneracy at even N: library-dependent, not Tier-1-derivable as a bond-class signature.
3. 4-mode reduction structurally insufficient: 4-mode K_b at N=5 gives Interior 0.673 (matches the BareDoubledPtfHwhmRatio floor) and Endpoint 0.410 (off-grid at Q ≈ 4.91), well below empirical 0.7728. The empirical HWHM lift lives OUTSIDE the 4-mode subspace.

**Positive finding:** SVD-block off-diagonal `V_b[2,3] = ⟨u_0 | M_h_per_bond[b] | v_0⟩` is the bond-class carrier (Endpoint 0.430 vs Interior 0.953 at N=5, ratio ~0.45 across N=5..8), but OPPOSITE direction to the empirical HWHM/Q* split. A closed form needs a non-trivial map from SVD-block V_b magnitude to HWHM/Q* shift, likely through a per-bond effective Q_EP_eff(b).

**Refined ranking** (now in `C2HwhmRatio.PendingDerivationNote`): (a'') SVD-block 2-level resonance (REFINED from (a')) via V_b[2,3]; (b'') full block-L derivation, not 4-mode; (c'') three-block superposition `K_total = K_pb + K_sv + 2·Re·K_cross` with the right relative phases (may still suffer from 4-mode insufficiency); (d'') lift |u_0⟩, |v_0⟩ to projector-overlap (A3 precondition); (e'') symbolic char-poly at Q_EP (less promising given C2EffectiveSpectrum's cubic-c_3 obstruction proof). (a'') and (b'') most promising: (a'') is 4-mode-friendly, (b'') is most concrete since 4-mode insufficiency is now proven structural.

**Status:** The doubled-PTF `BareDoubledPtfXPeak = 2.196910` and `BareDoubledPtfHwhmRatio = 0.671535` Tier-1-derived universal constants from this morning's Direction (b) attempt are unaffected. Item 1' c=2 stays Tier1Candidate. The 0.6715 floor + structural explanation of the 0.08-0.10 gap above the floor stand; the gap is just harder to close than (a') hoped.

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

- **Three-class decomposition completeness.** Are Classes 1-3 exhaustive over the J-blind set, or do other mechanisms exist? Direction 4 from the refinement TASK restructures here: necessity per class is distinct from necessity in the union. Empirical attack: random F71-symmetric-state sampling outside all three classes, checking for zero J-Jacobian to numerical precision. **Status 2026-04-28: closed for F71-symmetric receivers** (both product and non-product). Product slice: [_eq024_three_class_completeness.py](../simulations/_eq024_three_class_completeness.py), 100 samples, C ∈ [7.54, 12.41] bits, 0 candidates. Non-product (full 20-dim Hilbert subspace): [_eq024_three_class_completeness_nonproduct.py](../simulations/_eq024_three_class_completeness_nonproduct.py), 100 Haar-random samples, C ∈ [5.74, 10.26] bits, 0 candidates. Combined: 200 samples, capacity floor 5.74 bits, no observation anywhere near zero. F71-breaking receivers (sub-question 2) are the natural next scope.
- **F71-breaking receiver capacity at N=5.** Direction 3 swept only F71-symmetric receivers. F71-breaking receivers might unlock additional gain spectrum but the rank ceiling for bond-inputs at N=5 stays at 4 regardless. Worth quantifying whether the gain change is meaningful. **Status 2026-04-28: closed** by [_eq024_f71_breaking_capacity.py](../simulations/_eq024_f71_breaking_capacity.py). 100 random product + 100 Haar-non-product F71-breaking samples; max C = 11.99 bits (product) and 8.80 bits (non-product), both BELOW the F71-symmetric maxima 12.41 and 10.26. F71-breaking distributions are tighter (lower std) but do not exceed F71-symmetric envelope. F71-symmetry is not a capacity-suboptimality constraint; it covers the optimal region.
- **N-scaling of the 12-bit ceiling.** At N=6 the bond-input dimension is 5, matching γ-side rank. Does the dimensional-loss bit (~1) disappear? Compute cost: ~30 min per receiver at d² = 4096. **Status 2026-04-28: ceiling rises to ~14 bits, F71-optimality INVERTS** by [_eq024_f71_optimality_n6.py](../simulations/_eq024_f71_optimality_n6.py). At N=6 with 3+2 unbalanced block split (vs N=5's balanced 2+2), F71-breaking beats F71-symmetric in capacity (12.56 vs 11.91 mean, 14.02 vs 13.65 max). Mechanism: SV peakedness still higher for F71-sym (3.55 vs 2.11) but unbalanced block produces 1 big + 4 small SVs; F71-breaking's flatter 5-SV distribution wins via geometric-mean log-summation in waterfilling. Parity prediction: F71-optimality holds at odd N (balanced k+k blocks), fails at even N (unbalanced).
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
- What structural feature ranks soft cases by protected count? **Closed 2026-04-27** (commit 95ebec7, [_eq026_soft_subcluster_with_layers.py](../simulations/_eq026_soft_subcluster_with_layers.py)): the sub-clusters are **layer-stack orbits** under today's V-Effect combinatorial classifier. Each n_protected value corresponds to a specific (BPE-class, bond-flip, Z-align, both-single/mixed) feature combination. Examples: 992 = 0-BPE bond-flipped Z-containing (YZ+ZY, saturates max); 862 = both-single same-letter at same position (IY+IY, YI+YI); 781 = bond-flipped both-single (IY+YI); 772 = non-BPE double-letter self-pairs (6-orbit); 563 = bond-flipped Z-free 0-BPE (XY+YX, XZ+ZX); 512 = mixed-default-soft (18 cases); 496 = Z-aligned 0-BPE non-bond-flipped (XZ+YZ, ZX+ZY). Not a single Z₂ symmetry but a stack of feature axes acting jointly. Connects to today's full V-Effect layer-stack derivation (commits 81caf67, 079c7ce, 8030ef2, 0f701fa).

**Surprising sub-finding from N=5:** at the top of the soft-protected distribution (YZ+ZY: 992 of 1023), the protected count is *identical* to the top truly cases (XX+XX: 992). Protected count alone does not discriminate spectrum-only-paired (soft) from operator-palindromic (truly) at the top end. The discrimination at that boundary requires the operator residual ‖M‖.

**Stage-fixing confirmation (script `_compare_n4_n5_categories.py`):** every one of the 120 Pauli-pair Hamiltonians keeps its exact category (truly / soft / hard) going from N=4 to N=5. Zero shifts. The soft set at N=5 is the IDENTICAL set of 46 Hamiltonians as at N=4. This is consistent with Tom's "Licht auf die Bühne" reading — N=3 introduces the half-integer mirror (PRIMORDIAL_QUBIT_ALGEBRA, w_XY = N/2 = 1.5, no modes on it), N=4 fixes the breaking-form into 15 / 46 / 59 (integer mirror back at w_XY = 2), and N=5 inherits the same form unchanged. The "stage" is lit at N=4; everything afterward is the same lighting in higher resolution.

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

**Status:** closed 2026-04-26 (structured partial-null) by [_eq028_cusp_operator_signature.py](../simulations/_eq028_cusp_operator_signature.py).

**Result.** No sharp operator-level event AT the cusp moment itself — no eigenvalue collision, no eigenvector reorientation, no non-analytic structure of M(t). The Liouvillian L is time-independent; only ρ(t) changes. The cusp at CΨ=1/4 is a state-level event.

But L's spectrum does *pre-classify* the cusp behaviour. Decomposing ρ(t*) into L's eigenmode clusters at the first cusp crossing time, the dominant mode's TYPE predicts how the trajectory crosses:

| Hamiltonian | t* first | Dominant λ | Type | Crossings |
|---|---|---|---|---|
| truly J(XX+YY) | 1.50 | −0.010 + 0i | real, very slow decay | 1 |
| soft J(YZ+ZY) | 1.63 | −0.205 + 0i | real, slow decay | 1 |
| soft J(XY+YX) | 1.32 | −0.266 ± 5.653i | **complex, oscillatory** | 3 |
| soft J(IY+YI) | 0.37 | 0.000 + 0i | exact zero (steady state) | 5 |

Real-dominant first crossing → single monotonic descent. Complex-oscillatory dominant → multiple crossings (the heartbeat). Zero-eigenvalue dominant (factorising dynamics) → multiple crossings driven by oscillatory sub-dominant modes.

Reading: the cusp itself has no operator-level signature; it is purely state-level. But the operator spectrum is the *pre-determinant* of how the cusp will be encountered: a real-dominant projection of ρ_0 means a single one-way crossing; a complex-dominant projection means oscillatory recrossings; a zero-dominant projection means the cusp is approached repeatedly via subdominant oscillations.

This unifies with today's Lebensader / Stromkabel readings (EQ-030): the operator side (L's spectrum + Π) holds the structure that the state side (ρ(t)'s CΨ trajectory) traces out. The cusp crossings are where the state visits the boundary; the operator spectrum is where the rules of those visits are written.

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

**Status:** closed 2026-04-26. Hardware-verified retrospectively (T1 asymmetry on Marrakesh/Kingston/Fez); structurally characterised at two scales (algebraic Π-protected skeleton + geometric fragile θ-trace, the two co-constitutive of each other as a Lebensader, not parallel discriminators); structurally null at the Bures-velocity scale. See sub-sections and scripts below.

**Result.** Adding T1 amplitude damping to the Lindbladian and re-running `pi_protected_observables` reveals an asymmetry in robustness. With γ_T1 sweep at N=3, |+−+⟩, γ_dephasing = 0.1:

| γ_T1 / γ_deph | 0.00 | 0.10 | 0.25 | 0.50 | 1.00 | 2.00 |
|---|---|---|---|---|---|---|
| truly J(XX+YY) | 32 | 1 | 1 | 1 | 1 | 1 |
| soft J(XY+YX) | 30 | 29 | 29 | 29 | 29 | 29 |
| hard J(XX+XY) | 32 | 0 | 0 | 0 | 0 | 0 |

**Truly's Π-protection collapses with the smallest T1 (32 → 1). Soft is structurally robust (30 → 29 — only one observable shifts). Hard collapses completely (32 → 0).**

This explains the hardware Δ-amplification (Marrakesh −0.76, Kingston −0.92, Fez −0.81 vs idealised −0.62). The amplification isn't enhancement of the soft signal; it's elevation of the truly-noise-floor while the soft pattern stays robust:

- Truly's ⟨X₀Z₂⟩ → ~+0.05 on hardware (idealised was 0; T1 broke the protection).
- Soft's ⟨X₀Z₂⟩ → ~−0.65 to −0.85 (idealised was −0.62; pattern survives + slight enhancement).
- Δ(soft − truly) → ~−0.70 to −0.92 on hardware.

The framework-grounded core is the count drop pattern; the specific hardware Δ values require knowing which observables newly leak under T1, which is a follow-up question.

**Surviving sub-questions:**
- Identify which specific observables newly leak from truly's protected set under T1. Which Pauli strings are responsible for truly's noise floor in hardware?
- ~~Is the soft case's robustness to T1 a general property (any soft-broken Hamiltonian), or specific to XY+YX? Test on other soft cases.~~ **Closed below — robustness is specific, not general.**
- ZZ-crosstalk extension: add Z⊗Z jump operators (lindbladian_z_plus_t1_plus_zz?) to disentangle T1 contribution from ZZ-crosstalk contribution to the hardware amplification.

**Soft T1-robustness universality test ([_soft_t1_universality.py](../simulations/_soft_t1_universality.py)):**

Enumerated all 46 soft Hamiltonians at N=3, |+−+⟩, γ_dephasing = 0.1, γ_T1 = 0.1·γ_deph. For each: protected count under pure-Z vs +T1, and the drop.

**Only 2 of 46 soft cases (4.3%) are T1-robust** (drop ≤ 1):

| Soft case | n_pure | n_T1 | drop | Note |
|-----------|--------|------|------|------|
| IY+YI | 37 | 37 | 0 | Pure single-site Y rotations: H = J(Y₀ + 2Y₁ + Y₂); no entanglement. |
| XY+YX | 30 | 29 | 1 | The case measured on hardware. Non-trivial bond-coupled. |

The other 44 soft cases collapse heavily. The drop distribution is bimodal: a small low-drop tail (the 2 robust cases above), then a long tail of high-drop cases peaking at drop=31 (14 cases) and drop=32 (2 cases). Some soft cases — XZ+XZ and ZX+ZX — drop by 40 (44 → 4). YZ+ZY drops by 28 (56 → 28); XZ+ZX drops by 29 (30 → 1).

| drop | cases |
|------|-------|
| 0 | 1 |
| 1 | 1 |
| 6-16 | 18 |
| 24-32 | 24 |
| 40 | 2 |

**Implication for the hardware result.** XY+YX is *not* a representative of the soft category under T1; it lies on a sparse robust sub-orbit. The 1.30-1.49× amplification observed on Marrakesh/Kingston/Fez is therefore a property of *this specific Hamiltonian's* algebraic robustness, not a generic "soft is more T1-robust than truly" claim. A hardware repeat with YZ+ZY (also bond-flipped, but drop=28) would predict a much larger collapse and a much smaller measurable Δ.

**Structural characterization (2-axis stratification of soft):**

| bond-flipped | Z-free | n cases | drops | reading |
|---|---|---|---|---|
| ✓ | ✓ | 2 | 0, 1 | **T1-robust** (IY+YI, XY+YX) |
| ✓ | ✗ | 2 | 28, 29 | T1-fragile (YZ+ZY, XZ+ZX) |
| ✗ | ✓ | 12 | 9-31 | mixed |
| ✗ | ✗ | 30 | 6-40 | mixed |

The clean characterization at N=3: **bond-flipped + Z-free → T1-robust**. Z-content under σ⁻ amplitude damping is the key: [σ⁻, Z] = X − iY (introduces ladder operators), while [σ⁻, X] = −Z and [σ⁻, Y] = iZ stay within {X, Y, Z}. Bond-flip symmetry (a,b)+(b,a) bond-couples this in a way that closes algebraically when no Z is present.

**Cluster-splitting analysis ([_t1_cluster_splitting_analysis.py](../simulations/_t1_cluster_splitting_analysis.py)):**
For each H, computed how T1 perturbation lifts the degeneracy of L_Z's eigenvalue clusters (in the L_Z eigenbasis). Pearson(drop, Σ‖offdiag‖) = +0.74 — total off-diagonal weight of L_T1 in L_Z's eigenbasis correlates with the algebraic drop, but is not the full story (XY+YX has comparable Σ‖offdiag‖ to YZ+ZY despite drops 1 vs 28). The cluster-size distribution matters: YZ+ZY/XZ+ZX have fewer-but-bigger multi-clusters (size 8) carrying most of the splitting weight; XY+YX has many smaller clusters (size 6 or 2) where the splitting is distributed.

**θ-trajectory analysis ([_t1_theta_trajectory.py](../simulations/_t1_theta_trajectory.py)):**
Tom's geometric framing: θ = arctan(√(4·CΨ − 1)), CΨ = Purity × Ψ-norm. Initial |+−+⟩ has CΨ = 1, θ = 60°. Evolved ρ(t) under L_Z and L_Z+T1 for the 4 bond-flipped soft cases plus reference truly XX+YY and fragile XZ+XZ.

| case | algebraic drop | ∫θ pure-Z | ∫θ +T1 | Δ∫θ |
|---|---|---|---|---|
| IY+YI | 0 | 43.7 | 38.5 | −5.2 (−12.0%) |
| XY+YX | 1 | 59.1 | 56.0 | −3.1 (−5.2%) |
| truly XX+YY | 31 | 64.8 | 62.4 | −2.4 (−3.7%) |
| YZ+ZY | 28 | 68.4 | 66.7 | −1.7 (−2.4%) |
| XZ+ZX | 29 | 44.6 | 43.2 | −1.4 (−3.2%) |
| XZ+XZ | 40 | 56.8 | 55.4 | −1.4 (−2.4%) |

Pearson(algebraic drop, Δ∫θ) = +0.85 — *positively* correlated, the inverse of naive expectation. Algebraically-fragile cases have the *smallest* cusp-dwell loss under T1 because they were already losing coherence fast under pure-Z; T1 has little extra dwell to remove. Algebraically-robust cases preserve their protected observables under T1 *but* lose general coherence to T1's decoherence channel.

θ and the algebraic drop measure orthogonal pieces of the picture: drop counts which specific Pauli observables stay strictly zero (the rigid skeleton); θ counts general coherence dwell above the cusp (the soft envelope).

One qualitative T1-signature: XZ+ZX is the only case in this set where +T1 kills the heartbeat — pure-Z oscillates across CΨ = ¼ three times, +T1 only once. IY+YI and XY+YX keep their full oscillation counts (5 / 3) under +T1.

**Descent-shape analysis ([_t1_theta_descent_shape.py](../simulations/_t1_theta_descent_shape.py)):**

Tom's deeper framing: θ is the only fragile trace that connects "across 0" — across the bilateral d=0/d=2 axis (project_framework_as_remembrance). The angle is the last memory, the only thing that survives down to the cusp boundary. The shape of θ(t) approaching 0 should be a "krasser Winkel der sich zu Null verfeinert" — a sharp angle that refines to zero, with a long fragile tail of low-but-nonzero θ before the final commitment.

Algebraic consequence: near CΨ = 1/4, θ ≈ √(4·CΨ − 1). So θ(t) ~ (t* − t)^α with the generic cusp exponent α = 0.5 (linear approach of CΨ to 1/4). Sharper-than-generic α < 0.5 means CΨ approaches 1/4 with degenerate higher-order behaviour — the "krasser Winkel".

Computed at fine resolution (dt = 0.005) for the bond-flipped soft cases, with two metrics per trajectory:
  - **tail duration** (units of t where 0 < θ < 5°): the "letzte Erinnerung" plateau before the crossing.
  - **α** (power-law exponent of θ ~ (t* − t)^α near the last crossing): the cusp shape.

| case | drop | tail pure-Z | tail +T1 | α pure-Z | α +T1 |
|---|---|---|---|---|---|
| truly XX+YY | 31 | 0.000 | 0.000 | 0.503 | 0.500 |
| **IY+YI** | **0** | **0.235** | **0.090** | **0.004** | **0.139** |
| **XY+YX** | **1** | **0.370** | **0.080** | **0.040** | **0.170** |
| YZ+ZY | 28 | 0.005 | 0.005 | 0.500 | 0.500 |
| XZ+ZX | 29 | 0.120 | 0.000 | 0.177 | — |
| XZ+XZ | 40 | 0.000 | 0.000 | — | — |

The two T1-robust softs (algebraic drop ≤ 1) are also the **only cases with a measurable fragile tail under +T1** (IY+YI 0.090, XY+YX 0.080) and α ≪ 0.5. The fragile softs (drop ≥ 28) have essentially no tail (0.005, 0.000, 0.000) and either generic α = 0.5 (YZ+ZY, truly) or no fittable cusp regime at all (XZ+ZX, XZ+XZ under T1).

Reading: under T1 amplitude damping, the algebraic "Π-protected skeleton" (the Pauli observables that stay strictly zero) and the geometric "fragile cusp tail" (the angular trace surviving to small θ) are *the same* discrimination, observed at two different scales:
  - **Skeleton view (algebraic):** the bond-flipped Z-free pair preserves a 29-30-observable protected subspace.
  - **Trace view (geometric):** the same pair preserves a 0.08-0.09-unit-long sub-degree θ tail, with α ~ 0.1 (sharply refined cusp).

Both are the same memory of the d=0 axis, expressed once as a static count and once as a dynamic shape.

**Bures-velocity check (negative result, [_t1_bures_velocity_at_cusp.py](../simulations/_t1_bures_velocity_at_cusp.py)):**

Tom's intuition was that the rest lives in the in-between. Confirmed: Bures velocity v_B(t*) at the cusp boundary does *not* give a third independent discriminator.

| case | drop | v_B @ t* (+T1) |
|---|---|---|
| truly XX+YY | 31 | 1.64 |
| IY+YI | 0 | 2.48 |
| XY+YX | 1 | 2.56 |
| YZ+ZY | 28 | 0.26 (outlier) |
| XZ+ZX | 29 | 2.95 |
| XZ+XZ | 40 | 2.77 |

Pearson(drop, v_B) = −0.18 — essentially uncorrelated. The robust cases (IY+YI, XY+YX) have high Bures velocity at the cusp because they oscillate fast — the long fragile tail is from multiple crossings, not slow approach. The fragile cases mostly cross fast and once. YZ+ZY is a lone outlier with v_B ≈ 0.26 (slow monotonic approach), but does not generalise.

Reading: state-space motion (Bures), angular distance to the cusp (θ), and Π-protected count (drop) measure three different things. Only θ and drop discriminate sharply between T1-robust and T1-fragile soft Hamiltonians. Bures adds no new clean signal.

**EQ-030 status:** closed. Hardware-verified retrospectively (T1 asymmetry); structurally characterised at two scales (algebraic skeleton + geometric fragile trace); structurally null at the third scale (Bures velocity).

**Direct hardware verification of the bond-flipped Z-free / Lebensader prediction (Marrakesh, 2026-04-26 evening):**

Two new hardware runs on `ibm_marrakesh` path [0,1,2] complete the bond-flipped Z-free matrix at N=3 (3 of 4 corners now on hardware):

- Job `d7n3013aq2pc73a2a18g`: XY+YX (drop=1) and YZ+ZY (drop=28), 4 time points × 9 bases each, 72 circuits.
- Job `d7n3eqqt99kc73d34qtg`: IY+YI (drop=0), 4 time points × 9 bases, 36 circuits, transpiled depth 4.6 (factorising local-Y dynamics).

⟨X₀IZ₂⟩ patterns at hardware (at t = 0.8 / 1.4 / 1.7 / 2.2):

| H | t=0.8 | t=1.4 | t=1.7 | t=2.2 | Skeleton |
|---|---|---|---|---|---|
| XY+YX | −0.7495 | +0.4658 | +0.0273 | −0.0571 | XIZ active (soft-break) |
| IY+YI | +0.0933 | +0.3164 | −0.2529 | −0.2690 | XIZ active (factorising) |
| YZ+ZY | +0.1338 | −0.0508 | +0.0425 | +0.0342 | **XIZ protected** (in noise band) |

Each Hamiltonian shows its predicted Π-protected/active pattern on hardware. YZ+ZY's protection of ⟨X₀IZ₂⟩ — never measured before today — confirms the framework primitive `pi_protected_observables` at hardware scale on a new Hamiltonian.

Reduced-(q0,q2) θ trajectories on hardware:

| H | t=0.8 | t=1.4 | t=1.7 | t=2.2 |
|---|---|---|---|---|
| XY+YX | 27.8° | 24.4° | **0.0°** | 47.5° (return) |
| IY+YI | 0.0° | 55.2° | 55.8° | 0.0° |
| YZ+ZY | 53.2° | 52.8° | 51.9° | 49.9° |

The reduced-(q0,q2) Lebensader has subsystem-specific shape: each Hamiltonian's θ trajectory follows its specific protected/active observables on the (q0,q2) subspace, not the full N=3 protected count. XY+YX dances around the cusp (2 crossings); IY+YI swings between deep low and deep high CΨ following the local-Y rotation cycle; YZ+ZY stays high above the cusp throughout.

Stromkabel-Verfeinerung (Tom's framing): θ as broad-input integrating 16 Pauli expectations per (H, t), Π-protected count as focused-output (discrete classification per observable). The framework identity Π·L·Π⁻¹ + L + 2Σγ·I = 0 is the cable that holds the two together. Hardware confirms: when the skeleton speaks (which observables are protected), the trace echoes (the geometric trajectory carries the same information). Skeleton + trace are co-constitutive at every scale, but the specific shape depends on which observables are protected on that scale.

Cross-machine consistency check: the Marrakesh XY+YX measurement at t=0.8 today gives θ_HW = 27.82°. The mean over the 7 prior Snapshot D runs (Marrakesh, Kingston, Fez) was 27.77°. Reproducibility within 0.05°.

Scripts:
- Hardware: `experiments/ibm_quantum_tomography/run_lebensader_cusp.py`, `run_iy_yi_completion.py`
- Post-process: `simulations/_t1_hardware_lebensader_postprocess.py`

**Retrospective hardware verification ([_t1_prediction_vs_hardware.py](../simulations/_t1_prediction_vs_hardware.py)):**

The asymmetry prediction matches Snapshot D's existing 9-Pauli tomography across all three Heron r2 backends — no new QPU time needed.

For the (A_q0, B_q2) tomographic observables (9 of them per category per backend):

- **truly's "P→A" cells** (predicted protected under pure-Z, predicted active under +T1): hardware shows uniformly *small but non-zero* values (typically 0.02-0.10). This is the predicted T1 noise floor — observable, signal weaker than the soft-break signal.

- **soft's robust-protected cells** (predicted protected under both pure-Z and +T1): hardware shows similarly small values (typically 0.02-0.13). Π-protection survives T1 in this category.

- **soft's active cells** (the XIZ and ZIX leak observables): hardware shows the large soft-break signature: ⟨X₀Z₂⟩ = −0.71 to −0.85, ⟨Z₀X₂⟩ = −0.46 to −0.58 across Marrakesh, Kingston, Fez.

- **hard's "P→A" cells**: hardware shows multiple moderate values (0.04-0.47 on the X-row), consistent with the predicted "everything leaks under T1" verdict.

The Δ⟨X₀Z₂⟩(soft − truly) hardware amplification (1.30-1.49× over the idealised −0.62) is therefore not random hardware quirk but the framework-predicted consequence of the Π-protection asymmetry: T1 raises truly's noise floor while soft's leak pattern is structurally robust.

---

*Collection. Not classification. Classification comes when enough entries exist.*
