# PROOF: the signed balance of the flip masks (F173)

**Status:** Tier 1 derived. Under Z dephasing on every site (any strictly positive rates), for a Hamiltonian that is a real combination of Pauli strings in which every flip set (the X/Y mask of a string) carries exactly one string and the odd-weight diagonal part vanishes, the palindrome holds exactly when the signed Cayley graph of the masks is balanced. Three closed forms follow: two templates on the open chain, one template on the ring, and any number of templates on the ring, which reads N only through whether 2^(v₂(N)) exceeds the least valuation of the masks. Typed as `SignedMaskBalanceClaim`, with a live battery (`inspect --root signedbalance`) against Theorem 2 walked on the bitstrings.
**Date:** 2026-11-11
**Authors:** Thomas Wicht, Claude (Opus 5.5)
**Depends on:** Theorem 2 of [the complement connection](PROOF_PALINDROME_COMPLEMENT_CONNECTION.md) (the palindrome as the goodness of every component of the basis-state graph), and through it [F158](../ANALYTICAL_FORMULAS.md#f158).

<!-- STAGE 0 SWEEP RECORD
The question: which hopping cells pair under Z dephasing when the ratio H_xy / H_x̄ȳ of every edge is a sign, and where along N does that change.
docs/ANALYTICAL_FORMULAS.md: F115 (statement 3 reads Theorem 2 for the diagonal cell, every term of odd #Y + #Z, H_x̄ȳ = −H_xy, the Cayley graph of the window masks, bipartite exactly when no odd subset XORs to zero; statement 2 names the obstruction as the minimal odd relation among the shifts); F103 to F111 (F108 Part 1: every Π²_Z-even string commutes with X^⊗N, the colouring, so no such set is hard; F110 Aspect A: a Klein-homogeneous set outside the diagonal cell is coloured, its (1,1) cell by Y^⊗N; nothing on sets mixing the two parities); F116, F117, F138, F158, F169 (other mechanisms, no signed relation; F116's ring router, PROOF_CEILING_GOLDEN_ROUTER §8, exists exactly at 4 | N on the window mask 1 + x², whose flip sets carry strings of mixed signs, outside the assumption, the same factor (1 + x)² as Theorem C's XIY row).
docs/proofs: PROOF_PALINDROME_COMPLEMENT_CONNECTION (Theorem 2, the base; Theorem 3's proof already notes that a single string's ratio is one sign at every x); PROOF_F103_F87_Z2_CUBED_REFINEMENT §7.4 (bipartite iff some linear φ has φ(m) = 1 on every mask), §7.5 (FHF = −H in the diagonal cell), §7.7 (gcd, quotients, every relation a multiple s·(b, a), parity s(1)·(a(1) + b(1))), §7.7 to §7.8 (the generator has minimal degree, not minimal weight; cancellation from k = 6); PROOF_F111_HARD_CELL_PURE_D_TEMPLATE (pure-D templates; pairs of two mixed templates answered by Theorem 2, hard above the window floor N ≥ 2k − 2 exactly by F115's valuation criterion); PROOF_F87_WINDOWED_MONOMIAL_CONVERSE (the girth ladder of the trace route, a different object); MIRROR_SYMMETRY_PROOF (DM on odd cycles: the uniform map fails where tested, a statement about the map).
experiments: SOFTNESS_IS_N_DEPENDENT (one witness by hand, three strings on one flip set, the crossing as a square that first fits); THE_PALINDROME_CLASSIFIER (IXZX + XIZX soft below N = 6); F115_OBSTRUCTION_DISTRIBUTION (nothing on signs); V_EFFECT_FINE_STRUCTURE and TWO_TERM_PALINDROME_KLEIN_ROUTING (XY + ZZ carried by the alternating map X⊗Y⊗X…, hence on even rings; odd rings unstated); NON_HEISENBERG_PALINDROME (the alternating map does not exist on an odd cycle); THE_PALINDROME_AS_A_COLOURING (its frame theorem under uniform Z dephasing: the palindrome holds exactly when some diagonal phase gauge makes H commute with X^⊗N; it notes that F103's 2-colouring and its colouring are siblings that do not cite each other).
OpenArcs registry: f115_valuation_gate_width, gap (1) settled through Theorem 2, gap (2) the window floor, which Theorem B settles (the arc retired with this proof).
Typed layer: PalindromeComplementConnectionClaim (Theorem 2); ChiralKClaim (the sublattice sign within one excitation) and BipartiteGaugeCriterionClaim (F151, a gauge on coherence blocks) checked for adjacency, other objects; Diagnostics PalindromeMaskClassifier (the bipartite test, found to call bit_b = 0 sets hard and fixed the same day, CAUGHT_ERRORS 2026-11-11), WindowedObstructionScan and PalindromeSoftCertifier (diagonal cell only); MirrorWorld Hardness (the valuation), Stitch.LocalHolonomy (reads +XIY+XXX broken from five sites, no closed form); the MirrorWorld README's sieve row (on the window family's ring it counts 24 cells that pair by a law in N, read by no every-N route of the sieve; of these, 16 carry one string per flip set and fall under Theorem A, read for every N by Theorem D, and the other 8, X·(X ± Y)·X and relatives, carry two); MirrorWorld Stitch.FreeEveryN (a flip set carrying two strings on the ring, read through gcd(f, x^N + 1), the cyclotomic factors of x^N + 1 other than 1 + x that Theorem D's verdict does not see), checked for adjacency.
Glossary, CAUGHT_ERRORS, the whole repo: no "signed graph", "Harary" or "balance" in this sense.
Verdict: the all-negative case (bipartite) and the all-positive case (coloured) are held, in two places that do not cite each other; the mixed-sign criterion, the closed form for two templates with any signs, the chain threshold as the generator's footprint, and the ring laws for one template and for any number are not written anywhere. Adjacent: the colouring page's diagonal gauge, F103 §7.5's FK = X^⊗N·Z^φ and F115 statement 3 are one object; under the assumption the gauge can be taken to be a ±1 character.
-->

## The setting

N sites, each dephased along Z at a strictly positive rate. H is a real combination of Pauli strings. A string P has an X/Y mask (the sites carrying X or Y) and a sign

  s(P) = (#Y + #Z) mod 2,

so that X^⊗N·P·X^⊗N = (−1)^s(P)·P; s is the repo's bit_b, the Π²_Z parity. The **flip sets** of H are the masks of its strings; a string with mask 0 is diagonal.

**Assumption (one string per flip set).** No two strings of H have the same nonzero mask, and the odd-weight part of the diagonal of H vanishes (every diagonal string of odd weight has coefficient zero once equal strings are collected).

A **template** is a string of w letters, placed at every window: on the open chain at the starts 0, …, N − w, on the ring at every start, wrapped. Its mask is read as a polynomial m(x) over GF(2), bit j ↦ x^j, with the identity letters at its ends kept inside the window. The placed masks are x^p·m on the chain and x^p·m mod (x^N + 1) on the ring. Two templates whose masks are not translates of each other never place the same mask on the chain, nor on the ring of N ≥ 2·max(w₁, w₂) − 1 sites (a mask in an arc of at most (N + 1)/2 sites has one largest gap, which fixes how it sits on the ring); one template places N distinct masks on the ring of N ≥ 2w − 1 sites (a mask invariant under a rotation by k is a union of cosets of the subgroup generated by gcd(k, N), each spanning at least N − gcd(k, N) + 1 ≥ N/2 + 1 consecutive sites, which a window of w ≤ N/2 + 1/2 sites cannot hold).

## Theorem A (the balance)

Under the assumption, the palindrome holds exactly when **no GF(2) relation among the placed masks uses an odd number of negative strings**: for every subset S of the strings with ⊕_{P∈S} mask(P) = 0, Σ_{P∈S} s(P) is even. Equivalently, the Cayley graph of the masks, each edge signed by its string, is balanced in Harary's sense. Where it fails, no component is good and the far end is zero.

**Proof.** Theorem 2 of the complement connection: the palindrome holds exactly when every component of Γ is good. Take an edge x – x ⊕ m, m the mask of the string P = c·X^m Z^z (up to a phase). Its entries are H_{x,x⊕m} = c·u·(−1)^((x⊕m)·z) and H_{x̄,x̄⊕m} = c·u·(−1)^((x̄⊕m)·z), u the phase of P, so both are nonzero, of modulus |c| (condition (i)), and their ratio is (−1)^|z| = (−1)^s(P), the same at every x (|z| = #Y + #Z). No other string contributes to that edge, by the assumption. So Γ_H is the Cayley graph of the masks on each coset of their span, every edge present at every vertex. Condition (ii) is the vanishing of the odd-weight diagonal part, since X^⊗N·Z^v·X^⊗N = (−1)^|v|·Z^v.

A closed walk from x uses the strings P with multiplicities n_P and closes exactly when ⊕ (n_P mod 2)·mask(P) = 0, that is when ε_P = n_P mod 2 is a relation. Its holonomy is Π (−1)^(n_P·s(P)) = (−1)^(Σ ε_P·s(P)). Conversely every relation ε is walked from every vertex (the strings commute as flips, and each step stays inside Γ_H). So if every relation has even negative count, every cycle has trivial holonomy and every component is good; if one relation has odd negative count, every component carries a cycle of holonomy −1, none is good, and the far end is zero while the near end, the number of components, is not. ∎

**Remark (the same proof, wider).** The proof uses only that the ratio on each flip set is one sign at every x and that its entry never vanishes. A flip set may carry several strings of one sign whose entries never cancel; F115 statement 3 reads its y-parity-mismatched pairs that way. Where strings of different signs share a flip set, the ratio reads bits of x (the bits where their Z letters differ), and even cycles, squares among them, can fail: that is the mechanism of [SOFTNESS_IS_N_DEPENDENT](../../experiments/SOFTNESS_IS_N_DEPENDENT.md), outside this theorem.

**The two pure cases.**

- **Every sign negative** (every string of odd #Y + #Z): balanced means bipartite, no relation of odd size. This is [F103 §7.4](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) and F115 statement 3, the diagonal cell's hardness.
- **Every sign positive** (every string Π²_Z-even): every relation is balanced, and X^⊗N commutes with H, a colouring: [F108 Part 1](PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md). F110 Aspect A's Klein-homogeneous cells fall on both sides: the cells of even #Y + #Z are all positive, and its (1, 1) cell, odd popcount and odd #Y + #Z, is all negative with odd masks, where every relation uses an even number of strings (an odd number of odd-popcount masks cannot XOR to zero), coloured by Y^⊗N.

A balanced signing in general is a linear φ over GF(2) with φ(mask(P)) = s(P) for every string, and X^⊗N·Z^φ is then a lit string commuting with H. [The colouring page](../../experiments/THE_PALINDROME_AS_A_COLOURING.md) reads the palindrome under uniform Z dephasing as some diagonal phase gauge making H commute with X^⊗N; under the assumption Theorem A says that gauge can be taken to be a ±1 character, the lit string X^⊗N·Z^φ (F103 §7.5's FK in the all-negative case).

## Theorem B (two templates on the open chain)

Let H place two templates t₁, t₂ of widths w₁, w₂, with nonzero masks m₁, m₂ that are not translates of each other and signs s₁, s₂, plus any diagonal part with vanishing odd-weight part, every coefficient nonzero. Let g = gcd(m₁, m₂), q_i = m_i/g and

  c = s₁·q₂(1) + s₂·q₁(1)  (mod 2).

If c = 0, the chain pairs at every N. If c = 1, it pairs exactly for N < N* and breaks for every N ≥ N*, where

  N* = max(w₁ + deg q₂, w₂ + deg q₁).

**Proof.** A relation among the placed masks is a pair of polynomials u₁, u₂ (u_i's coefficient at x^p says whether the placement of t_i at p is used) with u₁·m₁ = u₂·m₂, deg u_i ≤ N − w_i. Dividing by g, u₁·q₁ = u₂·q₂ with gcd(q₁, q₂) = 1, so u₁ = r·q₂ and u₂ = r·q₁ for some polynomial r. The number of negative strings it uses is s₁·|u₁| + s₂·|u₂| ≡ s₁·u₁(1) + s₂·u₂(1) = r(1)·c (mod 2), the weight of a polynomial being its value at 1 modulo 2. If c = 0, every relation is balanced at every N, and Theorem A gives the palindrome. If c = 1, an unbalanced relation needs r(1) = 1. The choice r = 1 fits exactly when deg q₂ ≤ N − w₁ and deg q₁ ≤ N − w₂, that is N ≥ N*. Any r ≠ 0 has deg(r·q_j) = deg r + deg q_j ≥ deg q_j, so it needs N ≥ N* + deg r ≥ N*. Below N* every relation is balanced, from N* on r = 1 is unbalanced. ∎

**With both signs negative** (F115's diagonal cell), c = q₁(1) + q₂(1). The two quotients are coprime, so at most one is divisible by 1 + x; c = 1 exactly when one is and the other is not, that is when v_{1+x}(m₁) ≠ v_{1+x}(m₂): F115 statement 1, the valuation criterion. N* is its window floor, exactly; N ≥ 2k − 2, the floor F115 statement 3 and PROOF_F111 use, bounds it for every pair of width k (both masks are divisible by 1 + x, so deg q_i ≤ k − 2). At k = 4 six of the fourteen valuation-different mask pairs have N* = 5 and eight have N* = 6, the classifier's "six of the fourteen valuation-different mask pairs are already hard at N=5". For IXZX + XIZX, m₁ = x + x³, m₂ = 1 + x³, g = 1 + x, q₁ = x + x², q₂ = 1 + x + x², both of degree 2, both widths 4: N* = 6, soft below and hard from N = 6, as [THE_PALINDROME_CLASSIFIER](../../experiments/THE_PALINDROME_CLASSIFIER.md) records.

**With mixed signs**, outside F115:

- +XIY + XXX: m₁ = 1 + x² (s₁ = 1), m₂ = 1 + x + x² (s₂ = 0), g = 1, c = q₂(1) = 1, N* = 5. The relation r = 1 uses XIY at 0, 1, 2 and XXX at 0, 2: five strings, three of them negative, the only unbalanced relation at N = 5 (no square and no triangle carries it).
- +XIY + XIZ: m₂ = 1 (s₂ = 1), g = 1, c = 1·1 + 1·0 = 1, N* = 5: a triangle, {0, 2} = {0} ⊕ {2}.
- +XIY + XXZ: m₂ = 1 + x, g = 1 + x, q₁ = 1 + x, q₂ = 1, c = 1, N* = 4.
- +XIX + XXX: both signs positive, c = 0, every chain pairs.

On the window family of [the sieve](../../compute/MirrorWorld/Sieve.cs) (three-site templates with a letter at each end, two of them, signs ±) 864 cells are two-template cells with nonzero masks of different shapes: 432 have c = 0 and pair at every N, 432 have c = 1 and are, among the family's 528 mixed chain cells, the ones whose flip sets carry one string each, each pairing below its N* (4 or 5) and broken from it.

## Theorem C (one template on the ring)

Let H place one template t of width w with nonzero mask m and sign s on the ring of N ≥ 2w − 1 sites, plus any diagonal part with vanishing odd-weight part, every coefficient nonzero. Write v = v_{1+x} for the (1 + x)-adic valuation and v₂ for the 2-adic valuation of N. The ring breaks exactly when

  s = 1 and v(m) ≥ 2^(v₂(N)).

**Proof.** A relation is a polynomial u of degree below N with u·m ≡ 0 modulo x^N + 1. With h = gcd(m, x^N + 1), u is a multiple of (x^N + 1)/h, u = r·(x^N + 1)/h with deg r < deg h. Its negative count is s·u(1) = s·r(1)·((x^N + 1)/h)(1) (mod 2). Write N = 2^a·N' with N' odd: x^N + 1 = (x^N' + 1)^(2^a), and x^N' + 1 has 1 as a simple root, so v(x^N + 1) = 2^a. Then ((x^N + 1)/h)(1) = 1 exactly when v(h) = 2^a, and v(h) = min(v(m), 2^a). So an unbalanced relation exists exactly when s = 1 and v(m) ≥ 2^a, and then h ≠ 1 and r = 1 is one (deg h ≥ 1). ∎

The examples, each read by Theorem 2 on the bitstrings at N = 2w − 1 (from 3) to 10:

- XY (and XY + ZZ): m = 1 + x, v = 1, s = 1: breaks at every odd N, pairs at every even N. The even side is the alternating map X⊗Y⊗X… of V_EFFECT_FINE_STRUCTURE.
- XIY: m = (1 + x)², s = 1: breaks unless 4 | N. XIIIY: m = (1 + x)⁴: breaks unless 8 | N (its ring needs N ≥ 9; read at 9 and 10).
- XIIY: m = 1 + x³ = (1 + x)(1 + x + x²), v = 1: breaks at odd N.
- XYX, XXIY: v = 0, never breaks. YY, XIX, YIIIY: s = 0, never breaks.

## Theorem D (any number of templates on the ring)

Let H place templates t₁, …, t_k of widths w_i, with nonzero masks m_i no two of which are translates and signs s_i, on the ring of N ≥ 2·max w_i − 1 sites, plus any diagonal part with vanishing odd-weight part, every coefficient nonzero. Write T = 2^(v₂(N)), v_i = v(m_i) and v₀ = min v_i. The ring breaks exactly when

  T ≤ v₀ and some s_i = 1,  or  T > v₀ and two templates of valuation v₀ differ in sign or a template of greater valuation has s = 1.

The second condition does not involve N, so the ring reads N only through whether 2^(v₂(N)) exceeds v₀, and each mask only through its valuation. With one template it is Theorem C.

**Proof.** A relation is a vector (u₁, …, u_k) of polynomials modulo x^N + 1 (u_i's coefficient at x^p says whether the placement of t_i at p is used; the placements are distinct masks, by the setting) with Σ u_i·m_i ≡ 0 modulo x^N + 1, and its negative count is Σ s_i·u_i(1) mod 2. As in Theorem C, x^N + 1 = y^T·F with y = 1 + x and F(1) = 1. A relation reduced modulo y^T is a relation in L = GF(2)[x]/(y^T) with the same negative count, since u(1) is read modulo y. Conversely, if ū is a relation in L, then u_i = F·ū_i is one modulo x^N + 1 (Σ F·ū_i·m_i is F times a multiple of y^T), and u_i(1) = ū_i(1). So an unbalanced relation exists on the ring exactly when one exists in L.

Write m_i = y^(v_i)·e_i with e_i(1) = 1; in L, e_i is a unit, and m_i = 0 when v_i ≥ T. If v₀ ≥ T, every m_i is 0 in L, every vector is a relation, and the unit vector at i has negative count s_i. If v₀ < T, choose i₀ with v_{i₀} = v₀; for each j ≠ i₀ the vector ρ_j with u_j = 1, u_{i₀} = y^(v_j − v₀)·e_j·e_{i₀}^(−1) and the other entries 0 is a relation (u_{i₀}·m_{i₀} = y^(v_j)·e_j, which is m_j, also where v_j ≥ T and both are 0), and its negative count is s_j + s_{i₀}·[v_j = v₀], since y^(v_j − v₀) is 1 at x = 1 exactly when v_j = v₀ and a unit is 1 there. For any relation u, u − Σ_{j≠i₀} u_j·ρ_j is a relation whose only entry a sits at i₀, with a·y^(v₀)·e_{i₀} = 0 in L, so y^(T − v₀) divides a and a(1) = 0, T − v₀ being at least 1. Evaluation at 1 is a ring map, so the negative count of u is Σ_{j≠i₀} u_j(1)·(the negative count of ρ_j), and an unbalanced relation exists exactly when some ρ_j is unbalanced: s_j + s_{i₀} odd for a j of valuation v₀, which is two templates of valuation v₀ differing in sign (one of the two then differs from s_{i₀}), or s_j = 1 for a j of greater valuation. ∎

**At odd N** (T = 1) a mask has v ≥ 1 exactly when its weight is even, and the law reads: the ring breaks exactly when an even-weight mask carries an odd sign or two odd-weight masks carry different signs. **Where T > v₀**, the verdict for two templates is the long chain's of Theorem B: q_i(1) = 1 exactly when v_i is the least of the two, so c = s₂ when v₁ < v₂ and c = s₁ + s₂ when v₁ = v₂, which is what the law gives. Where T ≤ v₀, every mask is divisible by (1 + x)^T and the ring breaks on any odd sign.

The examples:

- +XXZ ± XZY, +XXZ ± YZX and their relatives, the 16 cells of the sieve's 24 law-in-N ring cells that carry one string per flip set: XXZ has v = 1 and s = 1, XZY v = 2 and s = 0. At odd N, T = 1 ≤ v₀ = 1 and XXZ's sign is odd, so the ring breaks; at even N, T ≥ 2 > v₀, XXZ alone has the least valuation and XZY, of greater valuation, has s = 0, so it pairs. Exactly at even N, at every N.
- XX + XIY: v = 1, 2 and s = 0, 1, broken on every ring: at T = 1 ≤ v₀ XIY's sign is odd, and from T = 2 on XIY has the greater valuation and s = 1.
- XX + XIIY (v = 1, 1, s = 0, 1) breaks on every ring; XY + XIIY (v = 1, 1, both signs odd) breaks at odd N and pairs at every even N.

**What the ring does not see.** The negative count is read at x = 1, so of the factors of x^N + 1 only 1 + x enters the verdict; the other cyclotomic factors drop out. Where one flip set carries several strings they do not: the window family's X·(X ± Y)·X, whose flip set 1 + x + x² carries two strings, pairs on the ring wherever gcd(1 + x + x², x^N + 1) = 1, read for every such N from N = 8 on by [the stitch's](../../compute/MirrorWorld/Stitch.cs) free chain, and breaks at every 3 | N read, N = 6, 9 and 12 (the end count at fixed N), exactly where 1 + x + x² divides x^N + 1.

## What is gated, and where

Typed as `SignedMaskBalanceClaim` (`compute/RCPsiSquared.Core/Symmetry/SignedMaskBalanceClaim.cs`, parent `PalindromeComplementConnectionClaim`), which holds the readings as GF(2) computations on masks and signs (`Balanced` at one N, `Chain` for Theorem B, `RingBreaks` for Theorems C and D). Its live battery, `SignedMaskBalanceWitness` (`inspect --root signedbalance`, `compute/RCPsiSquared.Diagnostics/Foundation/SignedMaskBalanceWitness.cs`, pinned by `compute/RCPsiSquared.Diagnostics.Tests/Foundation/SignedMaskBalanceWitnessTests.cs`), holds them against Theorem 2 walked on the bitstrings (`ComplementConnectionGraph`, every entry a Gaussian integer, every edge cross-multiplied), which never sees a mask relation or a sign, the templates placed with their identity letters:

- Theorem B at width three: every two-template cell of the window family with nonzero masks of different shapes, both signs of the second coefficient, at N* − 1, N*, N* + 1 (from N = 3) where c = 1 and at N = 3 to 6 where c = 0: 864 cells (432 breaking, 432 pairing), 3024 readings, the balance at one N agreeing at each.
- Theorem B at unequal widths, where N* pairs each width with the other template's quotient: 368 pairs of a two-site or three-site template with a wider one, 1036 readings from N* − 1 to N* + 1 where both templates are placed (where c = 0, N* is still the first N that holds a relation); XZ + XIXIX breaks from 6 and XY + XXIX from 5, where the widths swapped would say 9 and 7; IXZX + XIZX, both of width 4, pairs at N = 5 and breaks at 6 and 7.
- The F115 face: on the diagonal-cell masks at k = 3 and 4 (19 pairs of different shape), Theorem B breaks exactly where F115's own valuation (`WindowedObstructionScan.IsHardPair`) says the masks differ, and Theorem 2 pairs at N* − 1 and breaks at N*; at k = 4 the floors are six at N = 5 and eight at N = 6.
- Theorem C on 49 templates (the two-site and three-site ones with a nonzero mask and seven longer ones), N = 2w − 1 (from 3) to 10: 290 rings; and XY + ZZ, its diagonal of even weight, paired exactly at even N from 3 to 10.
- Theorem D against Theorem 2 on every pair of different shapes from the two-site and window templates with a nonzero mask (648 pairs, both signs of the second coefficient), and on the pairs with at least one of nine templates four or five sites wide (one sign), N = 2·max w − 1 (from 3) to 10: 1694 cells, 9140 readings. Against the balance, N = 2·max w − 1 (from 3) to 32, so that T runs from 1 to 32: on 1046 pairs up to five sites wide, 28232 readings, and on 416 triples of different shapes from sixteen templates with valuations 0 to 4, 10418 readings.
- Theorem A at one N on the 20 cells of three templates of different shapes in a pool of seven, chain N = 4 to 6 and ring N = 5, 6, the ring readings also by Theorem D.
- The worked examples, the all-positive face among them: XXI + XIX pairs on the chain at N = 4 to 6, where its odd relation fits, and on the ring of five, and the same masks with one sign negative, XXI + XIY, break from N* = 4; by Theorem D, the 16 ring cells above pair exactly at even N and XX + XIY breaks on every ring, N = 5 to 64.

Five controls inside the witness must disagree, and do: against Theorem 2, N* with each width beside its own quotient, c with the two quotients swapped, and the ring law without 2^(v₂(N)); against the balance, Theorem D with the parity of N in place of 2^(v₂(N)) and Theorem D read from a template of the greatest valuation in place of the least.
