# PROOF: Heisenberg Bonds off the Letters, a Colour as a Half-Turn

**Status:** Tier 1 derived: the half-turn law, for any coupling matrix; the split of the far equation into bond equations and field equations, for Heisenberg bonds under any integer dephasing axes with every site watched; the reduction that follows (wherever the bond equations admit only the colours, the palindrome holds exactly when an admitted colour exists, at any fields); and the bond equations admitting only the colours at N = 2, on the chain at N = 3, and on every connected graph when the axes are all parallel or pairwise parallel or perpendicular (there by Theorems 1 to 3 of the complement connection), so that with all axes parallel the rule is exact on every connected graph at fields of any direction. Exact row by row and not proven beyond: on 270 rows (the chain at N = 2 to 5, the ring at N = 3 to 5, the star and K₄ at N = 4) the bond equations admit exactly the colours (`TiltedAxesTests`). The holonomies of the last section read the same structure as angles, consistent with the proven part and proving nothing of their own: the triangle's exact at four triples, the others in floating point ([`simulations/heisenberg_colour_holonomies.py`](../../simulations/heisenberg_colour_holonomies.py)).
**Date:** 2026-10-27
**Authors:** Thomas Wicht, Claude (Opus 5.5)
**Builds on:** [F158](PROOF_PALINDROME_TWO_END_COUNT.md) (the palindrome holds exactly when F158's far kernel 𝒲 holds an invertible element; its §(h) covers jumps along any axis n̂·σ); [the palindrome as a colouring](../../experiments/THE_PALINDROME_AS_A_COLOURING.md) (colours on a circle); [F138](../ANALYTICAL_FORMULAS.md#f138), its paragraph "Off the letters"; [PROOF_PALINDROME_COMPLEMENT_CONNECTION](PROOF_PALINDROME_COMPLEMENT_CONNECTION.md) (Theorem 2, the complement connection; Theorems 1 and 3, the same rule at letter axes).

Dephase every site of a Heisenberg graph along an integer axis that need not be a Pauli letter. F138's paragraph "Off the letters" reads the palindrome there with colours: one direction perpendicular to every axis, with every field along it, so that the axes lie in one plane and the fields along its normal. A colour carries the palindrome; this page asks whether the palindrome needs one. It proves why a bond accepts a colour, splits the condition for the far kernel into a part the bonds own and a part the fields own, and shows that the fields' part is settled wherever the bonds admit nothing but the colours. Whether Heisenberg bonds ever admit more is the open question; XX + YY bonds do, through a sum of strings (F138's turned stage E row).

<!-- STAGE 0 SWEEP RECORD (for us, not for the reader: hidden on GitHub, read by every session, agent and hook)

A scout swept the stores by name for five topics: Heisenberg bonds under axes off the letters and any necessity there; the commutator identity of Lemma 2 or any grading by identity factors; conjugation by m·σ as a half-turn on a coupling matrix; solid angles, Bargmann invariants, Pancharatnam or Berry phases; spinor overlaps as a basis for tilted dephasing. A review round added two finds.

- docs/ANALYTICAL_FORMULAS.md, F138's paragraph "Off the letters": the rule, the cross form (W·a)×b = 0 = (W·b)×a for diagonal W, the coplanar row and the turned stage E row, the 4600 draws; F158 §(h) covers n̂·σ jumps. Nothing derives the cross form; nothing on the grading.
- docs/proofs/: PROOF_PALINDROME_COMPLEMENT_CONNECTION proves the rule for letter axes only (Theorems 1, 3, 5); its Theorem 5 and class B use a π rotation about an axis fixing every Heisenberg bond, Lemma 1 at W = 1; PROOF_PALINDROME_TWO_END_COUNT's element built from the half-turn about the bisector (X + Y)/√2 is Lemma 1 on one site; MIRROR_SYMMETRY_PROOF records that off the letters the clauses are not the law. "Holonomy" there is a plain cycle product, given as a sign (−1 on the triangle with three axes), never as an angle.
- experiments/: THE_PALINDROME_AS_A_COLOURING, "Colours on the circle", circle colours without the off-letters rule; TOPOLOGICAL_EDGE_MODES holds a Berry phase and the Pancharatnam closing rule for a biorthogonal path in parameter space, another object.
- compute/: TiltedAxes (the rule and the draws), Sextant.CircleColour, EndCount's HalfTurnAxis and EndCountTests' R(v) = 2(n·v)n − |n|²v (the same half-turn); the typed home of Theorem 2, PalindromeComplementConnectionClaim with ComplementConnectionGraph (Core) and ComplementConnectionWitness (Diagnostics), all at letter axes; SeedHolonomyClaim, an EP holonomy, another object.
- The OpenArcs registry: f138_clause_two_sweep names this question (where the colour is necessary off the letters, Heisenberg first).
- hypotheses/, reflections/, docs/GLOSSARY.md, docs/CAUGHT_ERRORS.md, fw.Confirmations, simulations/framework: nothing on these topics. No solid-angle or Bargmann reading of a holonomy anywhere.

Checked for adjacency: Theorems 1 to 3 of the complement connection give the bond equations' solutions wherever one global turn makes every axis a letter. The complement-connection proof and F138's paragraph are about the same necessity and did not cite each other; this page cites both, and both now point here. Theorem 5's π rotation, the bisector element and EndCount's HalfTurnAxis are Lemma 1 at W = 1 or on one site. The triangle's holonomy −1 of Theorem 3 (d) is the solid angle π/2 of the octant (last section).
-->

## Setup

N sites on a connected graph, H = Σ_{(i,j)} J_ij·σ_i·σ_j + Σ_l h_l·σ_l with every J_ij ≠ 0 and real and every field h_l real, and on every site l one jump n_l·σ_l along a nonzero integer axis n_l at a positive rate (the unit direction at rate γ|n_l|²). A **colour** is a direction m perpendicular to every axis; it is **admitted** when every field lies along it. By F158 the palindrome holds exactly when 𝒲 = {W : [H, W] = 0, (n_l·σ_l) W (n_l·σ_l) = −|n_l|²·W for every l} holds an invertible element, and an admitted colour gives one, ⊗(m·σ), which squares to |m|^{2N}. An operator anticommuting with n·σ on a site is there a combination of the w·σ with w ⊥ n, so 𝒲 lies in the span of the products ⊗_l (w_l·σ), w_l ⊥ n_l: the **far span**, one 2-dimensional slot per site. Its condition [H, W] = 0 is the **far equation**.

## Lemma 1 (the half-turn law)

For a two-site coupling B = Σ_{k,l} W_kl σ^k ⊗ σ^l with any real 3 × 3 matrix W and nonzero real vectors a, b: (a·σ)⊗(b·σ) commutes with B exactly when R_a·W·R_b = W, where R_v = (2vvᵀ − |v|²)/|v|² is the half-turn about v (symmetric, so R_v = R_vᵀ). Equivalently, (Wᵀ·a)×b = 0 = (W·b)×a. For a field h·σ on one site, a·σ commutes with it exactly when R_a·h = h, that is when h ∥ a or h = 0.

*Proof.* (a·σ) σ^k (a·σ) = Σ_l (2a_k a_l − |a|² δ_kl) σ^l, so conjugating by a·σ turns the Pauli vector by the half-turn about a, times |a|². Conjugating B by (a·σ)⊗(b·σ), which squares to |a|²|b|², replaces W by R_a W R_b; the product commutes with B exactly when that leaves B unchanged. For unit a, b, R_a W R_b = W reads 2(aᵀWb)·abᵀ = aaᵀW + Wbbᵀ; multiplying by aᵀ on the left gives aᵀW = (aᵀWb)·bᵀ, by b on the right W·b = (aᵀWb)·a, and conversely those two give the identity back. The field is the one-site case. ∎

F138 writes the cross form for diagonal W, where the transpose does not matter. Without it the form needs a symmetric W: for XX + XY under a = b = x it accepts, and the bond does not commute. For Heisenberg bonds, W = J·1, the law asks R_a = R_b: one colour on both ends of every bond, a single colour per component, and every field along it.

## Lemma 2 (the far equation splits)

For any vectors a, b (complex allowed),

  [σ_i·σ_j, (a·σ_i)(b·σ_j)] = 2i·(a×b)·(σ_j − σ_i),  [h·σ_l, w·σ_l] = 2i·(h×w)·σ_l .

So the commutator of a bond with a product in the far span has exactly one identity factor, at i or at j, and the commutator of a field with it has none. Sort [H, W] by the set of identity factors: W in the far span lies in 𝒲 exactly when

- **(B)** for every site i, Σ_{j~i} J_ij·(w_i × w_j) placed at slot j, the other slots as in W, vanishes (the terms with the identity at i), and
- **(F)** Σ_l (h_l × w_l) placed at slot l vanishes (the terms with no identity),

read linearly on W. (B) does not see the fields, and it is the far equation of the field-free H: its solutions are the far kernel of the bonds alone.

*Proof.* The first identity is a direct expansion, (σ^k)(a·σ) = a_k + i(e_k × a)·σ on each site, summed over k; the second is the one-site product rule. Operators with different sets of identity factors are linearly independent, so each set's part of [H, W] vanishes on its own; the bonds fill the sets of one identity, the fields the set of none. ∎

## The reduction

**Theorem.** If on the graph the bond equations (B) admit only the colours, that is their solutions are spanned by the products ⊗(m·σ) over the colours m, then the palindrome holds exactly when an admitted colour exists.

*Proof.* An admitted colour is in 𝒲 and invertible. Conversely, suppose no colour is admitted. If the axes span space there is no colour, (B) has no solution, and 𝒲 = 0. If some field has a component along its own axis, the palindrome fails by F158 §(f5): the word of H and that one jump has trace 2^N·h_l·n_l ≠ 0. Otherwise every field is perpendicular to its axis, and the colours span either one line m (the axes span a plane) or, when all axes are parallel to one n, the plane S ⊥ n, where the solutions of (B) are the symmetric tensors over S.

- *One line.* W = c·m^⊗N, and (F) reads Σ_l (h_l × m) at slot l, m elsewhere. The part of slot l perpendicular to m comes from the term l alone, so h_l × m = 0 for every l: every field along m, which is then admitted.
- *A plane.* A field h_l ⊥ n lies in S, and h_l × w for w ∈ S points along n; the terms with n at slot l come from the field at l alone, so W contracted at slot l with the in-plane turn n × h_l vanishes. For a symmetric tensor over S, read as a homogeneous polynomial P of degree N in two variables, that is a directional derivative ∂_q P = 0, so P is a power of the one linear form that vanishes on q, and W ∝ ĥ_l^⊗N with ĥ_l the direction of h_l. Two fields of different directions leave W = 0; fields of one direction m make m an admitted colour.

In every case without an admitted colour the palindrome fails. ∎

## Where the bond equations admit only the colours

**Per bond with two different axes.** Put m = n_i × n_j and p_i = n_i × m, p_j = n_j × m, so that {m, p_i} spans n_i^⊥ and {m, p_j} spans n_j^⊥. The output of (B) at i in slot j along n_j comes only from the bond (i, j) (every other term leaves slot j in n_j^⊥), and likewise the output of (B) at j in slot i along n_i. So the bond's 2 × 2 slot matrix of W (all other slots fixed) pairs to zero with det(·, ·, n_j) and with det(·, ·, n_i). On the pairs (m, m) and (p_i, p_j) both determinants vanish, since p_i, p_j, n_i, n_j all lie in the plane ⊥ m; on the mixed pairs (m, p_j) and (p_i, m) the two determinants form a 2 × 2 system of determinant |m|⁴·(cos²θ_ij − 1) ≠ 0 for unit axes, θ_ij the angle between them. So **the slot pair lies in span{m⊗m, p_i⊗p_j}**.

**A leaf.** If i has the single neighbour j, (B) at i is one term, and its whole output vanishes, not only the part along n_j: Σ_{a,b} M_ab·(e_a × e_b) = 0 for the bond's slot matrix M read as a 3 × 3 matrix, so M is symmetric. Its rows lie in n_i^⊥ and its columns in n_j^⊥, so with two different axes M is a multiple of m·mᵀ and **both slots carry m**; with one axis on both, M is a symmetric tensor over the common plane.

**N = 2.** Two leaves: different axes give m^⊗2, one axis gives the symmetric tensors over its plane, the colours either way.

**The chain at N = 3.** Both ends are leaves. If n₀ ∦ n₁ and n₁ ∦ n₂, slot 1 must carry n₀ × n₁ and n₁ × n₂, which are parallel exactly when the three axes are coplanar; otherwise W = 0, and when they are coplanar W = c·m^⊗3. If n₀ ∥ n₁ and n₁ ∦ n₂, the leaf 2 fixes slots 1 and 2 to m = n₁ × n₂, and the leaf 0's symmetric slot matrix with slot 1 along m forces slot 0 along m, which is perpendicular to n₀; the case n₀ ∦ n₁, n₁ ∥ n₂ is its mirror image. If all three are parallel, the two leaves make W symmetric in slots 0, 1 and in slots 1, 2, hence fully symmetric over the plane. In every case the solutions are the colours, which do solve (B) by Lemma 1.

So **on the chain at N = 3 the rule is exact**: Heisenberg bonds of any nonzero weights, axes anywhere, fields of any direction; the palindrome holds exactly when the axes have a common perpendicular m with every field along it.

**Axes the letters reach by one turn.** A global turn changes neither Heisenberg bonds nor the count, so where the axes are all parallel, or pairwise parallel or perpendicular, the field-free far end is read at letter axes by [PROOF_PALINDROME_COMPLEMENT_CONNECTION](PROOF_PALINDROME_COMPLEMENT_CONNECTION.md), on any connected graph: with all axes along one n, Theorems 1 and 2 give N + 1 good components, the popcount shells, whose flat sections span the symmetric tensors over the plane ⊥ n; with two letters Theorem 3 gives one, the third letter's colouring, and with three none. Each is the colours' dimension, so (B) admits only the colours there, and the reduction applies. In particular, **with all axes parallel the rule is exact on every connected graph at fields of any direction**: the palindrome holds exactly when every field is perpendicular to the common axis and all of them are parallel to each other. That reaches the fields Theorem 1 leaves out under Heisenberg bonds, a field along X + Y on every site among them; under XX + YY bonds it is not claimed here.

**Beyond.** On cycles and at interior sites with several neighbours, the terms of (B) from different neighbours can cancel in principle (each p_i⊗p_j part leaves a multiple of m at slot j), and that is how XX + YY bonds carry sums such as stage E's. For Heisenberg bonds the count says they do not: `TiltedAxesTests.Heisenberg_Bond_Equations_Admit_Exactly_The_Colours` builds (B) in integers on 270 rows (30 each on the chain at N = 2 to 5, the ring at N = 3 to 5, and the star and K₄ at N = 4), with random signed weights and axes generic, coplanar, parallel or in clusters of two and three, and finds the nullity over GF(p), an upper bound, equal on every row to the colours' dimension (N + 1, 1 or 0), an exact lower bound. That is exact per row and a measurement as a law.

## The same reading as holonomies

Theorem 2 of the complement connection holds for any Hermitian H with Z on every site, so it applies after turning each site's axis to z. Its far end is the flat sections of a connection on the bitstrings, d_y/d_x = H_xy/H_x̄ȳ, and a component carries one exactly when (i) |H_xy| = |H_x̄ȳ| on its edges, (ii) H_xx = H_x̄x̄ on its vertices and (iii) every cycle has trivial holonomy. In the turned frame a basis state is a product of spinors along ±n_l, and σ_i·σ_j = 2·SWAP − 1 makes every entry of a bond a product of two spinor overlaps; flipping every bit sends each spinor to its antipode. The holonomies then read as angles, each the doubled phase of a loop of overlaps. With α_l = s_l·n̂_l and s_l = (−1)^{x_l} at the starting vertex x:

| loop from x | holonomy | reading | status |
|---|---|---|---|
| triangle of pair flips, x → x⊕e₀⊕e₁ → x⊕e₀⊕e₂ → x | ((d + i·det)/(d − i·det))², d = 1 + a·b + b·c + c·a, (a, b, c) = (s₀n̂₀, −s₁n̂₁, −s₂n̂₂) | twice the solid angle of the three axes | derived; exact at four triples and all eight vertices |
| chain, the two pair flips in turn | e^{∓4iφ}, φ the angle at n₁ between the planes (n₀, n₁) and (n₁, n₂), the sign set by the bit of site 1 | the lune between n₁ and −n₁ | floating point, random weights |
| chain, single flips at the sites 2 and 1 in turn | e^{2i·arg(u_⊥·v_⊥ + i·a·(u×v))}, a = α₁, u = J₀₁α₀ − J₁₂α₂, v = J₀₁α₀ + J₁₂α₂ | the triple product det(n₀, n₁, n₂) | derived; floating point, random weights |
| all axes along z, fields at the angles θ₀ (site 0) and θ₂ (site 2) in the plane, 000 → 100 → 010 → 001 → 000 | e^{2i(θ₂ − θ₀)} | twice the angle between the fields | derived; floating point, weights 0.5 and −2 |

The reverse orientation of each loop gives the inverse. The overlaps around the triangle group into two Bargmann loops ⟨a|b⟩⟨b|c⟩⟨c|a⟩, whose phase is half the solid angle (Pancharatnam), and the complement reverses every loop's orientation, which doubles the phase in the ratio; the bond weights cancel in it (checked in floating point at random weights). Along z a field flip from 0 to 1 at site l has the entry h_x − i·h_y and its complement h_x + i·h_y, the ratio e^{−2iθ_l}, while the hops have ratio 1; that is the fourth row. On the triangle the values of d for the four sign classes of the vertices sum to 4, so they cannot all vanish, and with the axes not coplanar (det ≠ 0, no two axes parallel, every pair-flip entry nonzero) some vertex carries a nontrivial loop. At the letters X, Y, Z the solid angle is the octant's π/2 and the holonomy is the −1 of Theorem 3 (d).

The triangle and the single-flip square are trivial at every vertex exactly when the axes are coplanar (at a single vertex the triangle can also be trivial where d = 0, as for the axes (1, 1, 1), (−1, 1, 1), (1, −1, 1) at the vertices 000 and 111); the lune is trivial also at a right angle φ = π/2, where the letters X, Z, Y have no colour, which is why Theorem 3 needs its parity step there, and the single-flip square, the triple product, is that step off the letters. A single flip at a site is the neighbours' field turning it: its entry is ⟨α_l|m_l·σ|−α_l⟩ with m_l = Σ_{k~l} J_lk·α_k, the complement turns m_l with the spinors, so the ratio is −w/w̄ and condition (i) holds without fields. A field does not turn with the complement: condition (ii) asks it to be perpendicular to its own axis, and then condition (i) asks it to be perpendicular to every neighbour's axis, the local half of the reduction's field step.

## Scope, and what would falsify it

Lemma 1 holds for any real coupling matrix and any nonzero real vectors. Lemma 2 and the reduction hold for Heisenberg bonds of nonzero real weights on a connected graph with every site watched along an integer axis, fields of any direction. The rule is proven at N = 2, on the chain at N = 3, and on every connected graph with the axes all parallel or pairwise parallel or perpendicular. A counterexample to the measured part would be a Heisenberg row, every site watched, whose bond equations admit a solution that is no colour, or, through the reduction, a palindrome with no admitted colour; `TiltedAxesTests` would report it as a nullity above the colours' dimension. Unwatched sites are outside: there the SWAP rows of F138 (b) pair without a colour at letter axes already. Other bonds are outside Lemma 2's split as stated: their commutators with a product can leave terms without an identity factor.

## Open

1. The bond equations on a general connected graph: that the cancellations between neighbours at an interior site never happen under Heisenberg bonds. Settled where the axes are all parallel or pairwise parallel or perpendicular (Theorems 1 to 3). Otherwise the per-bond lemma and the leaf step pin a leaf on an axis different from its neighbour's, and what remains is the cancellation at interior sites with several neighbours that are not leaves, on trees and cycles alike, and leaves sharing their neighbour's axis.
2. The three-site ring as a theorem: the solid-angle holonomy shows a nontrivial loop at some vertex, and closing the argument needs that vertex's component to carry every far element (a connectivity statement for the turned hopping graph).
