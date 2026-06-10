# PROOF: the F87 windowed converse closes to a positive monomial (two-reflection theorem)

**Status:** Tier1Derived spine + Tier1Candidate full theorem (proven modulo **one** residual, R-sign, and only on the t_ℓ = 0 branch). General N. Sharpened 2026-06-09 (review wave): the deg-1 positivity gained the Pauli-coefficient closed form P_{3,1} = 6·4^N·Σ_l c_l². Sharpened again 2026-06-10 (the **girth dichotomy**, §4): the deg-1 class is in closed form at every m via the supertrace factorization, P_{2ℓ+1,1} = (2ℓ+1)·C(2ℓ,ℓ)·Σ_l t_ℓ² with t_ℓ = Tr(Z_l H^ℓ); the t_ℓ ≠ 0 branch is hard-at-all-γ **outright**, the t_ℓ = 0 branch climbs a girth ladder (m\* = 2ℓ + deg, deg odd; the γ³ rung when it fires, higher rungs when it does not, first γ⁵ witness IIXY+ZXZY) with R-sign open in ladder form; R-deg is retired (it was a k = 3 truth, refuted at k = 4 by IXXZ+XIXZ and replaced by something stronger).
**Date:** 2026-06-09
**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Builds on:**
- [PROOF_F103_F87_Z2_CUBED_REFINEMENT.md](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md): §7.5 (+N-Perron / −N-reflection mode), §7.7 ((1+x)-valuation), §7.10–§7.11 (what reaches the spectrum); the genericity result this strengthens.
- the F115 / WindowedHardnessClaim (1+x)-valuation criterion.
- ChiralKClaim ([compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs](../../compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs)): the chiral-K proof of bipartite ⟹ soft, re-proven here by a second (two-reflection) route.

## Abstract

A windowed diagonal-cell Pauli pair builds a chain Hamiltonian H, and the question F87 asks of it is whether the dephased Liouvillian's spectrum still pairs about the palindrome centre, the property [PROOF_F103](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) calls *soft*, or whether it genuinely fails to, which is *hard*. Recenter the Liouvillian at that centre and call the recentred superoperator **M(γ) = L + σ = A + γQ**, with A = −i[H, ·] the Hamiltonian commutator and Q = Σ_l Z_l ⊗ Z_l the (diagonal) dephasing generator, so that the spectrum of M is symmetric about 0 exactly when the pair is soft. The theorem is then a statement about M's odd power-sums.

**Theorem (windowed converse, monomial form).** For a non-bipartite windowed diagonal-cell pair, the first nonvanishing odd power-sum of M is a strictly positive monomial in γ. A positive monomial has no positive real root, so that power-sum is nonzero for every γ > 0, hence spec(M) ≠ spec(−M) for all γ > 0: the pair is **hard at every operating point**, not merely at all-but-finitely-many γ.

The proof splits cleanly into a rigorous spine and one sharp residual. **RIGOROUS-GENERAL** (general N, no premise): the all-odd word parity, the threshold #A ≥ 2ℓ (#A the total commutator-factor count, ℓ the unsigned odd-girth), the soft re-proof bipartite ⟹ soft, and the entire deg-1 class in closed form (the supertrace factorization through the moments t_j = Tr(Z_l H^j) and the sum-of-squares identity P_{2ℓ+1,1} = (2ℓ+1)·C(2ℓ,ℓ)·Σ_l t_ℓ², whose ℓ = 1 face is P_{3,1} = 6·4^N·Σ_l c_l²). These give the **girth ladder**: a hard pair with t_ℓ ≠ 0 is hard at every γ > 0 outright (m\* = 2ℓ+1, deg = 1, positive by the sum of squares); a hard pair with t_ℓ ≡ 0 fires at a higher odd rung, **proven modulo one residual** R-sign in ladder form (the first surviving class is single and positive; its k = 3 face is P_{2ℓ+3,3} > 0, and the γ³ rung can itself be silent, first γ⁵ witness IIXY+ZXZY), verified bit-exact over the entire N=4 k=3 Z diagonal cell, the full k=4 cell census, and the N=5 / N=6 representatives. The first nonvanishing moment sits at **m\* = 2ℓ + deg with deg odd, deg = 1 ⟺ t_ℓ ≠ 0** (the k = 3 cell realizes deg ∈ {1, 3} with the deg = 1 channel exactly the single-site-Z lift; k = 4 adds higher deg = 1 channels, IXXZ+XIXZ first, and the γ⁵ rung). The tempting m\* = 3ℓ guess is a low-ℓ coincidence, false already at ℓ = 5 where m\* = 13, not 15.

The computational anchors are the self-validating scripts [`simulations/f87_windowed_monomial_converse.py`](../../simulations/f87_windowed_monomial_converse.py) (the spine and the k=3 cell) and [`simulations/f87_girth_dichotomy.py`](../../simulations/f87_girth_dichotomy.py) (the factorization, the closed form, and the k=4 dichotomy battery).

## §1 The recentered object

Work on the d²-dimensional coherence space, d = 2^N, in the operator basis |i⟩⟨j| ↦ e_i ⊗ e_j (i-major, the row-stacking convention: left multiplication kron(·, I) acts on the bra index i and right multiplication kron(I, ·) acts on the ket index j). Recenter the Liouvillian at the palindrome centre σ = Nγ and write

  **M = A + γQ.**

The commutator part is

  A = −i[H, ·] = A_L + A_R,   A_L = −i(H ⊗ I)  (the left / bra hop),   A_R = +i(I ⊗ Hᵀ)  (the right / ket hop),

and the dephasing part is

  Q = Σ_l Z_l ⊗ Z_l,   diagonal,   Q_{ij} = N − 2·popcount(i ⊕ j) ∈ [−N, N].

(Each Ad_{Z_l} has eigenvalues ±1, so Q's diagonal entry reads off the Hamming distance between the bra and ket labels; this is the Absorption-Theorem form of the dephaser, [PROOF_ABSORPTION_THEOREM.md](PROOF_ABSORPTION_THEOREM.md).) By construction M(γ) equals the framework Lindbladian `lindbladian_pauli_dephasing(H, [γ]*N, 'Z')` recentred by Nγ, bit-for-bit, so the generator is the same one [PROOF_F103](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) classifies, just shifted so that softness becomes a statement about symmetry around 0.

**The standard power-sum lemma.** A finite multiset of complex numbers is symmetric about 0 if and only if all of its odd power-sums vanish. One direction is immediate, a symmetric multiset pairs λ with −λ and λ^m + (−λ)^m = 0 for odd m; the converse is the classical fact that the power-sums determine the multiset (Newton's identities recover the elementary symmetric functions, hence the characteristic polynomial, hence the roots-with-multiplicity). Applied to spec(M),

  **spec(M) symmetric about 0 ⟺ p_m(γ) := Tr(M^m) = 0 for every odd m.**

Expand each power-sum in the dephasing strength,

  p_m(γ) = Σ_j P_{m,j} γ^j,   P_{m,j} = Σ (traces of length-m words in {A_L, A_R, Q} with exactly j factors of Q),

since M = (A_L + A_R) + γQ and a length-m word picks up one factor of γ per Q it contains. The whole theorem is read off the structure of the coefficients P_{m,j}: which (m, j) first fail to vanish, how many j survive there, and the sign of the survivor.

## §2 The two reflections (general N, including complex H)

Two involutions on coherence space do all the parity bookkeeping. Let F = X^⊗N (so F² = I), and set

  𝓕 = F ⊗ F,   R = I ⊗ F.

Both square to the identity, so each is its own inverse. Their action on the three building blocks is a fixed sign table, and the table holds for **complex H** (a flux pair, odd #Y, has an X…Y bond that makes H Gaussian-integer, not real), which is the case the older chiral-K route did not have to face.

**The driving lemma.** Every diagonal-cell term has #Y + #Z odd (that is the (0,1) Klein cell's defining bit_b parity, [PROOF_F103 §Notation](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md)), and conjugation by X flips the sign of each Y and each Z while fixing each X and I. So **F H F = −H**, term by term. For complex H this is all we need: transposing gives **F Hᵀ F = (F H F)ᵀ = −Hᵀ**, with no appeal to H = Hᵀ. Likewise F Z_l F = −Z_l on every site.

**The sign table.** From F H F = −H and F Hᵀ F = −Hᵀ,

  𝓕 A_L 𝓕 = −A_L,   𝓕 A_R 𝓕 = −A_R,   𝓕 Q 𝓕 = +Q,

and from R = I ⊗ F (which touches only the ket factor) together with F Z_l F = −Z_l,

  R A_L R = +A_L,   R A_R R = −A_R,   R Q R = −Q.

(𝓕 sees both tensor legs, so it flips A_L and A_R alike and leaves Q untouched, two Z-sign flips cancelling; R sees only the ket leg, so it flips A_R and Q and fixes A_L.) These are RIGOROUS-GENERAL operator identities; the verification script checks all six bit-exactly, including the complex-H flux pair IXY+XIY and at N=5.

**The all-odd consequence (RIGOROUS-GENERAL).** Conjugation by an involution preserves each word's trace, Tr(GWG⁻¹) = Tr(W). For a length-m word W with counts (#A_L, #A_R, #Q), the sign table gives

  Tr(W) = (−1)^{#A} Tr(W)   (from 𝓕, with #A = #A_L + #A_R),
  Tr(W) = (−1)^{#A_R + #Q} Tr(W)   (from R).

A trace that equals its own negative is 0, so a word survives only if #A is even and #A_R + #Q is even. Now take m odd. #A even forces #Q = m − #A odd; #A_R + #Q even with #Q odd forces #A_R odd; and #A even with #A_R odd forces #A_L odd. Hence

  **every word surviving an odd power-sum has #A_L, #A_R, #Q all odd.**

This is proven per word from the sign table; the script confirms it two ways, by checking the count-arithmetic implication over every triple summing to small odd m, and by an exhaustive word census at a cheap scale (the multi-Z lift at m\* = 5, all 3^5 words: every nonzero word has all-odd counts and the census sums to the exact p_5).

## §3 The threshold and the soft re-proof

Read a surviving trace as a closed double-walk. Since Tr(W) = Σ_{ij} ⟨ij|W|ij⟩, and A_L moves only the bra index, A_R only the ket index, Q neither, a diagonal element ⟨ij|W|ij⟩ is nonzero only if the bra index returns to i after its #A_L left-hops (a closed walk of length #A_L on H's hopping graph G_H) and the ket index returns to j after its #A_R right-hops (a closed walk of length #A_R). Here G_H is the graph on the 2^N basis states with an edge for every nonzero off-diagonal entry of H **and a self-loop at i for every nonzero diagonal entry H_{ii}**; a self-loop is an odd closed walk of length 1, so a diagonal lift has ℓ = 1 and a graph with a loop is never bipartite. With that convention, an odd-length closed walk exists in a graph if and only if the graph is **non-bipartite**, and the minimal odd length is the **unsigned odd-girth** ℓ.

**Path-existence, not signed cancellation.** The point worth stressing is that Q is diagonal, so a trace's support is governed by *index-trajectory existence*, the unsigned question (|H|^k)_{ii} > 0, and is immune to the signed XX+YY cancellation that the [PROOF_F103 §7.10](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) retraction warned about. The flux pair is the clean witness: its signed third-moment amplitude (H³)_{ii} cancels to 0, yet the unsigned three-walk count (|H|³)_{ii} = 6 > 0, so a closed odd walk genuinely exists and the threshold below still binds. (The script asserts exactly this, sgn (H³)_{ii} = 0 while unsigned > 0, for the flux pair.) Two conclusions follow.

- **Bipartite ⟹ soft (RIGOROUS-GENERAL).** If G_H is bipartite there is no odd closed walk, so #A_L odd is impossible, so by §2 no word survives any odd power-sum, so every odd p_m ≡ 0, so spec(M) is symmetric about 0: **soft**. This is a second, independent proof of bipartite ⟹ soft, complementary to the chiral-K route of [ChiralKClaim](../../compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs) and [PROOF_F103 §7.1](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md). Where the chiral-K argument exhibits a similarity W with W L W⁻¹ = −L − 2σ, this route never builds an operator: it simply observes there are no surviving words. The script's soft control (the bipartite pair XXZ+ZXX) confirms all odd p_m vanish exactly.

- **Non-bipartite ⟹ the threshold.** If G_H is non-bipartite then both the bra walk and the ket walk must reach the odd-girth, #A_L ≥ ℓ and #A_R ≥ ℓ, so every word surviving an odd moment has **#A ≥ 2ℓ**. Combined with #Q ≥ 1 odd from §2,

  **m\* ≥ 2ℓ + 1.**

  (That the minimal class #A = 2ℓ is also the one that actually fires is part of the law m\* = 2ℓ + deg, settled by the girth dichotomy of §4; the spine needs only the inequality.)

The threshold #A ≥ 2ℓ is RIGOROUS-GENERAL (it is the unsigned odd-girth path-existence statement, verified including the flux pair).

## §4 Monomial structure and degree-1 positivity

At the first nonvanishing odd moment m\*, the surviving (#A, #Q) classes have #A even and ≥ 2ℓ and #Q = m\* − #A odd, so p_{m\*} is a sum of **odd** γ-powers (one per surviving #Q). It is a **monomial** exactly when a single value of #Q survives. The single-Q coefficient has a clean cyclic form: collecting the words with #Q = 1 and using cyclicity of the trace to slide the lone Q to the front,

  **P_{m,1} = m · Tr(Q · A^{m−1}).**

**Degree-1 positivity (RIGOROUS-GENERAL, closed form).** At the lowest possible degree, m = 3 and #A = 2, the single-Q coefficient is P_{3,1} = 3·Tr(A²Q). Anti-Hermiticity of A (it is −i times a Hermitian commutator) gives a first identity: the diagonal of A² is manifestly non-positive, (A²)_{xx} = Σ_y A_{xy} A_{yx} = −Σ_y |A_{xy}|² ≤ 0, so with deg_A(x) = Σ_y |A_{xy}|² the coherence-space out-weight at x and Q_x the diagonal of Q,

  P_{3,1} = 3·Tr(A²Q) = −3·Σ_x deg_A(x)·Q_x = 6·Σ_x deg_A(x)·(w(x) − N/2),

where w(x) = (N − Q_x)/2 = popcount of the bra-ket difference at coherence index x. This form is exact but does not yet show the sign: the weights w(x) − N/2 are mixed-sign, and the sum lands positive only after cancellation. The sign comes from evaluating the same trace tensorially. With A = −i(H ⊗ I − I ⊗ Hᵀ),

  A² = −H² ⊗ I + 2·H ⊗ Hᵀ − I ⊗ (Hᵀ)²,

and tracing against Q = Σ_l Z_l ⊗ Z_l kills the two single-leg terms on Tr(Z_l) = 0, leaving Tr(A²Q) = 2·Σ_l Tr(HZ_l)·Tr(HᵀZ_l) = 2·Σ_l Tr(HZ_l)², real even for a complex flux H since H is Hermitian. Reading off Pauli coefficients (Tr(HZ_l) = 2^N·c_l, with c_l the coefficient of the single-site string Z_l in H),

  **P_{3,1} = 6·4^N·Σ_l c_l²,**

manifestly non-negative, and strictly positive exactly when H carries a single-site-Z component. For the single-site-Z diagonal lift (the ℓ = 1, deg = 1 case) this gives 6·4⁴·(1² + 2² + 1²) = 9216 (N=4) and 6·4⁵·(1² + 2² + 2² + 1²) = 61440 (N=5), each matching the exact P_{3,1} bit-for-bit (the script checks both forms against the exact γ¹ coefficient of p₃).

**The factorization at every m (2026-06-10).** The tensor evaluation is not an m = 3 accident. Expanding A^{2k} = (−1)^k(H ⊗ I − I ⊗ Hᵀ)^{2k} binomially and tracing each leg against Q = Σ_l Z_l ⊗ Z_l gives the supertrace factorization

  **Tr(Q · A^{2k}) = (−1)^k · Σ_l Σ_{j=0}^{2k} (−1)^j C(2k, j) · t_j^{(l)} · t_{2k−j}^{(l)},   t_j^{(l)} = Tr(Z_l H^j):**

the entire deg-1 class is a bilinear form in the **Z_l-weighted moments** of H. Two RIGOROUS-GENERAL kills then collapse it. The F-chirality (𝓕's driving lemma F H F = −H together with F Z_l F = −Z_l) gives t_j = (−1)^{j+1} t_j, so every even-j moment vanishes. And the unsigned odd-girth kills t_j for j < ℓ outright on a pure off-diagonal H (the diagonal of H^j is a closed-walk count, and no closed walk is shorter than the girth). At the first candidate moment m = 2ℓ+1 the binomial sum therefore holds a single surviving term, and it is a **sum of squares**:

  **P_{2ℓ+1, 1} = (2ℓ+1) · C(2ℓ, ℓ) · Σ_l (t_ℓ^{(l)})²  ≥ 0.**

**The girth dichotomy.** Everything about the deg-1 class now reads off one computable quantity, the girth moment t_ℓ:

- **t_ℓ ≠ 0 at some site ⟹ m\* = 2ℓ + 1 with deg = 1, proven outright.** The threshold kills every odd moment below 2ℓ+1, and at 2ℓ+1 the #Q = 3 class needs #A = 2ℓ − 2 < 2ℓ, so p_{2ℓ+1} = P_{2ℓ+1,1}·γ exactly: a monomial whose coefficient is the sum of squares above, positive with no further input. **Hard at every γ > 0, no residual.**
- **t_ℓ ≡ 0 at every site ⟹ the deg-1 class is dead at 2ℓ+1 and at 2ℓ+3** (the binomial terms of P_{2ℓ+3,1} pair t_ℓ with t_{ℓ+2}, and t_ℓ = 0 kills them all). If the γ³ class fires there, p_{2ℓ+3} = P_{2ℓ+3,3}·γ³ exactly: the **monomial property is proven at that rung**. And if it does not fire, the ladder continues to higher odd rungs: the k = 4 census produced the first γ⁵ witness, IIXY+ZXZY, whose γ³ class also dies and whose first nonvanishing moment is p₁₁ = 86507520·γ⁵ at m\* = 11 = 2ℓ+5, again a single positive power. The open residual (§5) is the ladder's coherence: that the first surviving class is always single and positive.

The single-site-Z lift is the ℓ = 1 face of the dichotomy: t_1 = 2^N·c_l, and the sum-of-squares formula reproduces P_{3,1} = 6·4^N·Σ_l c_l² verbatim. The taxonomy sentence this proof first carried, "deg = 1 only for a single-site-Z diagonal lift," is a **k = 3 fact, not a general one**: the k = 3 cell happens to have t₃ ≡ 0 for all 16 of its pure cycles, but at k = 4 exactly 20 of the 192 hard pure-cycle pairs carry t₃ ≠ 0 and fire deg = 1 at m\* = 7. The first representative is IXXZ+XIXZ: t₃ = ±64 at N = 5, p₇ = 573440·γ bit-exactly equal to 7·C(6,3)·Σ t₃², single power, positive. The same windowed GF(2) machinery that classifies hardness also classifies the branches: production of a Z_l word requires the syzygy polynomial B = (p₂/g)q₁ + (p₁/g)q₂ of the templates' XY- and YZ-mask polynomials to be a monomial, and when it is, the word signs follow the Cartier-Foata graph-inversion character, cancelling for the k = 3 flux channels (and for every same-y-parity-1 pair, where word reversal is a sign-reversing involution) and surviving for the k ≥ 4 deg-1 cycles.

## §5 The residual (open: R-sign, on the t_ℓ = 0 branch only)

This section originally carried two residuals. The first, **R-deg**, asserted that for every pure hopping cycle the deg-1 class dies and the monomial lifts to degree 3; the 2026-06-10 girth dichotomy (§4) **dissolved it**. As formulated, R-deg was a k = 3 truth: the k = 4 census produces twenty pure-cycle pairs (IXXZ+XIXZ first among them) where the deg-1 class fires at m\* = 2ℓ+1 with a sum-of-squares-positive coefficient, so "the deg-1 class always cancels for cycles" is simply false in general. What replaces it is stronger than what it asked for: the deg-1 class is now in closed form at every m, the monomial property is **proven in both branches** of the dichotomy, and the t_ℓ ≠ 0 branch needs no residual at all. The old k = 3 evidence stands unchanged and is re-derived by the dichotomy (all 16 pure cycles of the N=4 Z cell have t₃ = 0 exactly, which is why the k = 3 world only ever showed deg = 3 cycles): the K3 pair has p_9 = 2064384·γ³ at N=4 and 16515072·γ³ at N=5, the complex-H flux pair p_9 = 589824·γ³, and the ℓ=5 / N=6 representative p_1…p_11 = 0 with p_13 = 50381979648·γ³ (the `--heavy` certificate).

What remains is one residual, in ladder form:

- **R-sign (the first surviving class is single and positive).** For a hard pair with t_ℓ ≡ 0, the first nonvanishing odd moment sits at some m\* = 2ℓ + deg with odd deg ≥ 3, and the open statement is that the surviving class there is **single** (the lower #Q classes vanish together, which is proven for #Q = 1 but not yet for the intermediate rungs once deg ≥ 5) **and positive**. The k = 3 face of this is "P_{2ℓ+3,3} > 0 for genuine cycles," which is [PROOF_F103 §7.5](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md)'s +N-population-Perron top-skew: the gain channel's +N Perron mode is unpaired because its −N reflection partner is absent on a non-bipartite graph, and the skew lands with a definite, positive sign. But the k = 4 census shows the γ³ rung can itself be silent (IIXY+ZXZY has P_{2ℓ+3,3} = 0 and fires at γ⁵ instead, positively), so the residual is genuinely the ladder statement, not a single coefficient. Verified: all 16 pure cycles of the N=4 Z cell (16/16 Perron), the stratified k=4 battery (every sampled pair a single positive power on the ladder), and the N=5/N=6 representatives. What is missing is a closed-form identity writing the first surviving P_{m\*,deg} as a manifestly-non-negative functional, the higher-#Q analogue of §4's sum-of-squares form for #Q = 1.

**What the residual buys.** Once R-sign holds, the first nonvanishing odd moment is a single positive γ-power on both branches, and the structural law is complete:

  **m\* = 2ℓ + deg,   deg odd,   deg = 1 ⟺ the girth moment fires (t_ℓ ≠ 0).**

At k = 3 the cell realizes only deg ∈ {1, 3}, with the deg = 1 channel exactly the single-site-Z lift (t₁ = 2^N c_l); at k = 4 the higher girth moments join the deg = 1 channel (the B-monomial production of §4) and the γ⁵ rung appears. The earlier **m\* = 3ℓ** guess remains a low-ℓ coincidence, false at ℓ = 5 where m\* = 13, not 15. The verification is carried by two scripts: the anchor [`simulations/f87_windowed_monomial_converse.py`](../../simulations/f87_windowed_monomial_converse.py) (Blocks 2, 3, 4, 6 RIGOROUS-GENERAL; Block 7 the k=3 cell-wide certificate) and the dichotomy battery [`simulations/f87_girth_dichotomy.py`](../../simulations/f87_girth_dichotomy.py) (the factorization checked exactly, the closed form at all four branch representatives, the full k=4 census with all 20 deg-1 cycles verified coefficient-exact, and the k=3 regression).

## §6 Consequence

Put the pieces together. For a non-bipartite windowed diagonal-cell pair the first nonvanishing odd power-sum is a positive monomial c·γ^deg with c > 0 and deg odd: **outright** on the t_ℓ ≠ 0 branch (the sum-of-squares coefficient), and **modulo R-sign in ladder form** on the t_ℓ = 0 branch (the existence of m\* is the proven genericity result; what is open is that the first surviving class there is single and positive, as every verified case is). A positive monomial has no positive real root, so

  **p_{m\*}(γ) > 0 for all γ > 0   ⟹   spec(M) ≠ spec(−M) for all γ > 0   ⟹   non-bipartite is hard at every operating point.**

This is the upgrade. [PROOF_F103 §7.6](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) closed the windowed converse to "hard for all but finitely many γ" by a degenerate-perturbation-theory plus analyticity argument: a first-order break (c ≠ 0) forces some recentred characteristic-polynomial coefficient Δ_j(γ) to be a nonzero polynomial, whose zero-set is finite. The residual lemma there, the part flagged Tier1Candidate, was exactly "no positive γ is one of those finitely-many accidental soft points," a resultant / Sturm question on Δ_j. The monomial theorem answers it structurally: the first nonvanishing odd moment is a positive monomial, hence nonzero at every γ > 0, so the spectrum is asymmetric at every operating point and there are no accidental soft points left to rule out. "All but finitely many γ" becomes "all γ > 0," **unconditionally on the t_ℓ ≠ 0 branch and modulo R-sign on the t_ℓ = 0 branch** (until 2026-06-10 this read "modulo R-deg + R-sign"; the girth dichotomy retired R-deg).

The result is carried by two typed claims and two scripts. The Tier1Derived spine, the threshold and the all-odd word parity and the soft re-proof and the deg-1 closed forms (the c_l² identity and its general-m sum-of-squares form, §4), is the node `WindowedConverseThresholdClaim`. The Tier1Candidate full lemma, now the monomial-and-positive statement gated on R-sign alone, is the node `WindowedConverseAllGammaClaim` (the same residual [PROOF_F103](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) isolated, sharpened twice: from a 700-point numerical spot-check to a polynomial-root statement, and from two residuals to one). The reproducible anchors are [`simulations/f87_windowed_monomial_converse.py`](../../simulations/f87_windowed_monomial_converse.py) and [`simulations/f87_girth_dichotomy.py`](../../simulations/f87_girth_dichotomy.py), self-validating: every block raises on failure, prints a single PASS line on success, and the process exits 0 only if the whole spine / residual ledger holds.
