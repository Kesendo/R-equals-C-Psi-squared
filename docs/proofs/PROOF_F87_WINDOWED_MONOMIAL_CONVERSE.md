# PROOF: the F87 windowed converse closes to a positive monomial (two-reflection theorem)

**Status:** Tier1Derived, **no residual** (closed 2026-06-10, §5). General N. Sharpened 2026-06-09 (review wave): the deg-1 positivity gained the Pauli-coefficient closed form P_{3,1} = 6·4^N·Σ_l c_l². Sharpened 2026-06-10 (the **girth dichotomy**, §4): the deg-1 class is in closed form at every m via the supertrace factorization, P_{2ℓ+1,1} = (2ℓ+1)·C(2ℓ,ℓ)·Σ_l t_ℓ² with t_ℓ = Tr(Z_l H^ℓ); the t_ℓ ≠ 0 branch is hard-at-all-γ **outright**, the t_ℓ = 0 branch climbs a girth ladder (m\* = 2ℓ + deg, deg odd); R-deg retired. **Closed later the same day (the Pascal-Gram positivity theorem, §5): the supertrace factorization generalizes to every #Q class, every γ-coefficient of the first nonvanishing odd power-sum is an equal-total sum of squares or exactly zero, so p_{m\*}(γ) > 0 for every γ > 0. R-sign is resolved, and the all-γ converse holds unconditionally.**
**Date:** 2026-06-09
**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Builds on:**
- [F103 refinement](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md): §7.5 (+N-Perron / −N-reflection mode), §7.7 ((1+x)-valuation), §7.10–§7.11 (what reaches the spectrum); the genericity result this strengthens.
- the F115 / WindowedHardnessClaim (1+x)-valuation criterion.
- ChiralKClaim ([compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs](../../compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs)): the chiral-K proof of bipartite ⟹ soft, re-proven here by a second (two-reflection) route.

## Abstract

A windowed diagonal-cell Pauli pair builds a chain Hamiltonian H, and the question F87 asks of it is whether the dephased Liouvillian's spectrum still pairs about the palindrome centre, the property [F103](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) calls *soft*, or whether it genuinely fails to, which is *hard*. Recenter the Liouvillian at that centre and call the recentred superoperator **M(γ) = L + σ = A + γQ**, with A = −i[H, ·] the Hamiltonian commutator and Q = Σ_l Z_l ⊗ Z_l the (diagonal) dephasing generator, so that the spectrum of M is symmetric about 0 exactly when the pair is soft. The theorem is then a statement about M's odd power-sums.

**Theorem (windowed converse, positivity form).** For a non-bipartite windowed diagonal-cell pair, every γ-coefficient of the first nonvanishing odd power-sum of M is non-negative: each surviving #Q class is an equal-total **Pascal-Gram sum of squares**, and every other class vanishes exactly. At least one class is positive, so p_{m\*}(γ) > 0 for every γ > 0, hence spec(M) ≠ spec(−M) for all γ > 0: the pair is **hard at every operating point**, not merely at all-but-finitely-many γ.

The proof is rigorous end to end (general N, no premise, no residual): the all-odd word parity, the threshold #A ≥ 2ℓ (#A the total commutator-factor count, ℓ the unsigned odd-girth), the soft re-proof bipartite ⟹ soft, the deg-1 class in closed form (the supertrace factorization through the moments t_j = Tr(Z_l H^j) and the sum-of-squares identity P_{2ℓ+1,1} = (2ℓ+1)·C(2ℓ,ℓ)·Σ_l t_ℓ², whose ℓ = 1 face is P_{3,1} = 6·4^N·Σ_l c_l²), and then (§5) the same factorization for **every** #Q class: each coefficient P_{m,d} factorizes through d-leg moments T^{(l⃗)}_{α⃗} = Tr(Z_{l₁}H^{α₁}···Z_{l_d}H^{α_d}), and at the first nonvanishing odd moment a parity + girth + cascade argument leaves only the equal-leg-total block, which the Vandermonde identity assembles into Σ|U|² with prefactor +1. The structural law is **m\* = 2ℓ + deg with deg odd, deg = 1 ⟺ t_ℓ ≠ 0** (the **girth dichotomy**): the k = 3 cell realizes deg ∈ {1, 3} with the deg = 1 channel exactly the single-site-Z lift; k = 4 adds higher deg = 1 channels (IXXZ+XIXZ first) and the γ⁵ rung (IIXY+ZXZY). A selection rule (§5) even pins the surviving class uniquely for deg ≤ 3, so the older monomial form of the theorem is *derived* there; from deg = 5 on two classes may coexist and positivity, which is all the converse needs, carries alone. The tempting m\* = 3ℓ guess is a low-ℓ coincidence, false already at ℓ = 5 where m\* = 13, not 15.

The computational anchors are the self-validating scripts [`simulations/f87_windowed_monomial_converse.py`](../../simulations/f87_windowed_monomial_converse.py) (the spine and the k=3 cell), [`simulations/f87_girth_dichotomy.py`](../../simulations/f87_girth_dichotomy.py) (the factorization, the closed form, and the k=4 dichotomy battery), and [`simulations/f87_pascal_gram_positivity.py`](../../simulations/f87_pascal_gram_positivity.py) (the Pascal-Gram chain of §5, block by block).

## §1 The recentered object

Work on the d²-dimensional coherence space, d = 2^N, in the operator basis |i⟩⟨j| ↦ e_i ⊗ e_j (i-major, the row-stacking convention: left multiplication kron(·, I) acts on the bra index i and right multiplication kron(I, ·) acts on the ket index j). Recenter the Liouvillian at the palindrome centre σ = Nγ and write

  **M = A + γQ.**

The commutator part is

  A = −i[H, ·] = A_L + A_R,   A_L = −i(H ⊗ I)  (the left / bra hop),   A_R = +i(I ⊗ Hᵀ)  (the right / ket hop),

and the dephasing part is

  Q = Σ_l Z_l ⊗ Z_l,   diagonal,   Q_{ij} = N − 2·popcount(i ⊕ j) ∈ [−N, N].

(Each Ad_{Z_l} has eigenvalues ±1, so Q's diagonal entry reads off the Hamming distance between the bra and ket labels; this is the Absorption-Theorem form of the dephaser, [Absorption Theorem](PROOF_ABSORPTION_THEOREM.md).) By construction M(γ) equals the framework Lindbladian `lindbladian_pauli_dephasing(H, [γ]*N, 'Z')` recentred by Nγ, bit-for-bit, so the generator is the same one [F103](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) classifies, just shifted so that softness becomes a statement about symmetry around 0.

**The standard power-sum lemma.** A finite multiset of complex numbers is symmetric about 0 if and only if all of its odd power-sums vanish. One direction is immediate, a symmetric multiset pairs λ with −λ and λ^m + (−λ)^m = 0 for odd m; the converse is the classical fact that the power-sums determine the multiset (Newton's identities recover the elementary symmetric functions, hence the characteristic polynomial, hence the roots-with-multiplicity). Applied to spec(M),

  **spec(M) symmetric about 0 ⟺ p_m(γ) := Tr(M^m) = 0 for every odd m.**

Expand each power-sum in the dephasing strength,

  p_m(γ) = Σ_j P_{m,j} γ^j,   P_{m,j} = Σ (traces of length-m words in {A_L, A_R, Q} with exactly j factors of Q),

since M = (A_L + A_R) + γQ and a length-m word picks up one factor of γ per Q it contains. The whole theorem is read off the structure of the coefficients P_{m,j}: which (m, j) first fail to vanish, how many j survive there, and the sign of the survivor.

## §2 The two reflections (general N, including complex H)

Two involutions on coherence space do all the parity bookkeeping. Let F = X^⊗N (so F² = I), and set

  𝓕 = F ⊗ F,   R = I ⊗ F.

Both square to the identity, so each is its own inverse. Their action on the three building blocks is a fixed sign table, and the table holds for **complex H** (a flux pair, odd #Y, has an X…Y bond that makes H Gaussian-integer, not real), which is the case the older chiral-K route did not have to face.

**The driving lemma.** Every diagonal-cell term has #Y + #Z odd (that is the (0,1) Klein cell's defining bit_b parity, [F103 §Notation](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md)), and conjugation by X flips the sign of each Y and each Z while fixing each X and I. So **F H F = −H**, term by term. For complex H this is all we need: transposing gives **F Hᵀ F = (F H F)ᵀ = −Hᵀ**, with no appeal to H = Hᵀ. Likewise F Z_l F = −Z_l on every site.

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

**Path-existence, not signed cancellation.** The point worth stressing is that Q is diagonal, so a trace's support is governed by *index-trajectory existence*, the unsigned question (|H|^k)_{ii} > 0, and is immune to the signed XX+YY cancellation that the [F103 §7.10](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) retraction warned about. The flux pair is the clean witness: its signed third-moment amplitude (H³)_{ii} cancels to 0, yet the unsigned three-walk count (|H|³)_{ii} = 6 > 0, so a closed odd walk genuinely exists and the threshold below still binds. (The script asserts exactly this, sgn (H³)_{ii} = 0 while unsigned > 0, for the flux pair.) Two conclusions follow.

- **Bipartite ⟹ soft (RIGOROUS-GENERAL).** If G_H is bipartite there is no odd closed walk, so #A_L odd is impossible, so by §2 no word survives any odd power-sum, so every odd p_m ≡ 0, so spec(M) is symmetric about 0: **soft**. This is a second, independent proof of bipartite ⟹ soft, complementary to the chiral-K route of [ChiralKClaim](../../compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs) and [F103 §7.1](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md). Where the chiral-K argument exhibits a similarity W with W L W⁻¹ = −L − 2σ, this route never builds an operator: it simply observes there are no surviving words. The script's soft control (the bipartite pair XXZ+ZXX) confirms all odd p_m vanish exactly.

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

**The m = 3 face is cell-free (2026-06-10, the PTF tie-in).** Nothing in the m = 3 computation used the diagonal cell. The three companion coefficients vanish for *every* Hermitian H: Tr(A³) dies by the binomial antisymmetry of the A_L/A_R expansion (the j ↔ 3−j terms cancel on Tr((Hᵀ)^j) = Tr(H^j)), Tr(Q²A) dies on Tr(Z_l Z_{l'}) = 0 for l ≠ l' and on Tr(H) − Tr(Hᵀ) = 0 for l = l', and Tr(Q³) is a sum of traces of nonempty Z-strings. So for any Hermitian H, p₃(γ) = 6·4^N·Σ_l c_l²·γ exactly: a single-site-Z Pauli component breaks the palindrome at every γ > 0, whatever cell the rest of H lives in. The verifier is [`simulations/f87_deg1_face_cell_free.py`](../../simulations/f87_deg1_face_cell_free.py) (the PTF Π-break configuration XY chain + εZ_m, site-dependent fields, and a random complex Hermitian H, all exact, with the bipartite control vanishing). This is what makes the Z-field row of the PTF Π-break experiment ([PTF under palindrome-breaking perturbations](../../experiments/PTF_PALINDROME_BREAKING_PERTURBATIONS.md)) theorem-grade: its measured mirror break is all-γ with leading coefficient 6·4^N·ε²·γ. (The higher rungs of the ladder do use the cell parity; cell-freeness is established here for the m = 3 face only.)

**The factorization at every m (2026-06-10).** The tensor evaluation is not an m = 3 accident. Expanding A^{2k} = (−1)^k(H ⊗ I − I ⊗ Hᵀ)^{2k} binomially and tracing each leg against Q = Σ_l Z_l ⊗ Z_l gives the supertrace factorization

  **Tr(Q · A^{2k}) = (−1)^k · Σ_l Σ_{j=0}^{2k} (−1)^j C(2k, j) · t_j^{(l)} · t_{2k−j}^{(l)},   t_j^{(l)} = Tr(Z_l H^j):**

the entire deg-1 class is a bilinear form in the **Z_l-weighted moments** of H. Two RIGOROUS-GENERAL kills then collapse it. The F-chirality (𝓕's driving lemma F H F = −H together with F Z_l F = −Z_l) gives t_j = (−1)^{j+1} t_j, so every even-j moment vanishes. And the unsigned odd-girth kills t_j for j < ℓ outright on a pure off-diagonal H (the diagonal of H^j is a closed-walk count, and no closed walk is shorter than the girth). At the first candidate moment m = 2ℓ+1 the binomial sum therefore holds a single surviving term, and it is a **sum of squares**:

  **P_{2ℓ+1, 1} = (2ℓ+1) · C(2ℓ, ℓ) · Σ_l (t_ℓ^{(l)})²  ≥ 0.**

**The girth dichotomy.** Everything about the deg-1 class now reads off one computable quantity, the girth moment t_ℓ:

- **t_ℓ ≠ 0 at some site ⟹ m\* = 2ℓ + 1 with deg = 1, proven outright.** The threshold kills every odd moment below 2ℓ+1, and at 2ℓ+1 the #Q = 3 class needs #A = 2ℓ − 2 < 2ℓ, so p_{2ℓ+1} = P_{2ℓ+1,1}·γ exactly: a monomial whose coefficient is the sum of squares above, positive with no further input. **Hard at every γ > 0, no residual.**
- **t_ℓ ≡ 0 at every site ⟹ the deg-1 class is dead at 2ℓ+1 and at 2ℓ+3** (the binomial terms of P_{2ℓ+3,1} pair t_ℓ with t_{ℓ+2}, and t_ℓ = 0 kills them all). If the γ³ class fires there, p_{2ℓ+3} = P_{2ℓ+3,3}·γ³ exactly: the **monomial property is proven at that rung**. And if it does not fire, the ladder continues to higher odd rungs: the k = 4 census produced the first γ⁵ witness, IIXY+ZXZY, whose γ³ class also dies and whose first nonvanishing moment is p₁₁ = 86507520·γ⁵ at m\* = 11 = 2ℓ+5, again a single positive power. That the ladder always lands on non-negative coefficients, whichever rung fires, is the Pascal-Gram positivity theorem of §5 (it was the day's open residual, R-sign in ladder form, for the hours between the two waves).

The single-site-Z lift is the ℓ = 1 face of the dichotomy: t_1 = 2^N·c_l, and the sum-of-squares formula reproduces P_{3,1} = 6·4^N·Σ_l c_l² verbatim. The taxonomy sentence this proof first carried, "deg = 1 only for a single-site-Z diagonal lift," is a **k = 3 fact, not a general one**: the k = 3 cell happens to have t₃ ≡ 0 for all 16 of its pure cycles, but at k = 4 exactly 20 of the 192 hard pure-cycle pairs carry t₃ ≠ 0 and fire deg = 1 at m\* = 7. The first representative is IXXZ+XIXZ: t₃ = ±64 at N = 5, p₇ = 573440·γ bit-exactly equal to 7·C(6,3)·Σ t₃², single power, positive. The same windowed GF(2) machinery that classifies hardness also classifies the branches: production of a Z_l word requires the syzygy polynomial B = (p₂/g)q₁ + (p₁/g)q₂ of the templates' XY- and YZ-mask polynomials to be a monomial, and when it is, the word signs follow the Cartier-Foata graph-inversion character, cancelling for the k = 3 flux channels (and for every same-y-parity-1 pair, where word reversal is a sign-reversing involution) and surviving for the k ≥ 4 deg-1 cycles. The involution is the Pauli transpose antiautomorphism θ(ρ) = ρᵀ (word reversal = transpose × (−1)^{n_Y(word)}, sign-reversing exactly at y-parity 1): the same θ whose 2-letter commutator face is [F114](../ANALYTICAL_FORMULAS.md)'s D·L_σ·D = (−1)^{n_Y(σ)+1}·L_σ and whose Hermitian-conjugacy face is F113's Lemma C (already cited as step 3 of §5).

## §5 The residual resolved: the Pascal-Gram positivity theorem (2026-06-10)

This section originally carried two residuals, and both dissolved on the same day. The first, **R-deg**, asserted that for every pure hopping cycle the deg-1 class dies and the monomial lifts to degree 3; the girth dichotomy (§4) **dissolved it** in the morning wave. As formulated, R-deg was a k = 3 truth: the k = 4 census produces twenty pure-cycle pairs (IXXZ+XIXZ first among them) where the deg-1 class fires at m\* = 2ℓ+1 with a sum-of-squares-positive coefficient, so "the deg-1 class always cancels for cycles" is simply false in general. The old k = 3 evidence stands unchanged and is re-derived by the dichotomy (all 16 pure cycles of the N=4 Z cell have t₃ = 0 exactly): the K3 pair has p_9 = 2064384·γ³ at N=4 and 16515072·γ³ at N=5, the complex-H flux pair p_9 = 589824·γ³, and the ℓ=5 / N=6 representative p_1…p_11 = 0 with p_13 = 50381979648·γ³ (the `--heavy` certificate).

The second, **R-sign in ladder form** ("at the first nonvanishing odd moment the surviving class is single and positive"), is resolved by generalizing §4's supertrace factorization from the #Q = 1 class to **every** #Q class. The resolution proves more than the residual asked: positivity holds for *all* classes at m\* simultaneously, and singleness, where it holds, is *derived* rather than assumed.

**Theorem (Pascal-Gram positivity).** Let m\* be the first odd m with p_m(γ) = Tr(M^m) not identically zero. Then for every d, either P_{m\*,d} = 0 exactly, or

  **P_{m\*,d} = (m\*/d) · Σ_{l⃗ ∈ [N]^d} Σ_{k⃗} |U^{(l⃗)}_{k⃗}|² ≥ 0,   U^{(l⃗)}_{k⃗} = Σ_{|α⃗| = u} ∏_i C(α_i, k_i) · T^{(l⃗)}_{α⃗},**

with u = (m\* − d)/2 the common leg total and T^{(l⃗)}_{α⃗} = Tr(Z_{l₁}H^{α₁}Z_{l₂}H^{α₂}···Z_{l_d}H^{α_d}) the **d-leg moments** of H. Since p_{m\*} ≢ 0, at least one class is positive, and p_{m\*}(γ) > 0 for every γ > 0.

**Proof.** Seven steps, each verified exactly in [`simulations/f87_pascal_gram_positivity.py`](../../simulations/f87_pascal_gram_positivity.py).

1. **Cyclic decomposition.** Collecting the words with d factors of Q by trace-cyclicity (rotate each so a Q leads; each composition term is hit m\*/d times, symmetric necklaces included),

   P_{m,d} = (m/d) · Σ_{a⃗} Tr(Q A^{a₁} Q A^{a₂} ··· Q A^{a_d}),  a⃗ running over the compositions of m − d into d parts.

2. **Leg factorization.** A_L and A_R commute, so A^a = Σ_j C(a,j) A_L^{a−j} A_R^j; Q = Σ_l Z_l ⊗ Z_l; and the supertrace of a tensor word splits as Tr(bra)·Tr(ket). Writing α_i = a_i − j_i (bra exponents) and β_i = j_i (ket exponents), the bijection (a⃗, j⃗) ↔ (α⃗, β⃗) turns the composition sum into a free double sum:

   Σ_{a⃗} Tr(Q A^{a₁} ··· Q A^{a_d}) = Σ_{l⃗} Σ_{α⃗,β⃗: |α⃗|+|β⃗| = m−d} (−i)^{|α⃗|} (+i)^{|β⃗|} ∏_i C(α_i + β_i, β_i) · T^{(l⃗)}_{α⃗} · T̃^{(l⃗)}_{β⃗},

   where T̃ is the moment built on Hᵀ.

3. **Hermitian conjugacy.** H is Hermitian, so Hᵀ = H̄, and the Z_l are real: the ket-leg moment is the entrywise conjugate of the bra-leg moment at the *same* indices, **T̃^{(l⃗)}_{β⃗} = conj(T^{(l⃗)}_{β⃗})**. (No index reversal needed; this is the transpose-trick sibling of F113's Lemma C.)

4. **Leg parity and leg girth.** Conjugating the single trace by F = X^⊗N (F H F = −H, F Z_l F = −Z_l) gives T^{(l⃗)}_{α⃗} = (−1)^{d + |α⃗|} T^{(l⃗)}_{α⃗}: a d-leg moment vanishes unless **|α⃗| ≡ d (mod 2)**, so for odd d only odd totals survive. And the diagonal entries of the leg product are weighted closed-walk counts of total length |α⃗| (the Z insertions do not move the index), so on a pure off-diagonal H no leg total below the unsigned odd-girth survives: **T^{(l⃗)}_{α⃗} = 0 for |α⃗| < ℓ**. (Diagonal lifts have ℓ = 1 and the bound is vacuous.)

5. **Vandermonde assembly of the equal-total block.** On the block |α⃗| = |β⃗| = u the phase is (−i)^u(+i)^u = +1, and the slot-wise Vandermonde identity C(α+β, β) = Σ_k C(α,k)·C(β,k) factorizes the double sum:

   Σ_{|α⃗|=|β⃗|=u} ∏_i C(α_i+β_i, β_i) · T_{α⃗} · conj(T_{β⃗}) = Σ_{k⃗} [Σ_{α⃗} ∏C(α_i,k_i) T_{α⃗}] · conj[Σ_{β⃗} ∏C(β_i,k_i) T_{β⃗}] = Σ_{k⃗} |U^{(l⃗)}_{k⃗}|².

6. **Slice inversion.** For |k⃗| = u, the only composition with ∏C(α_i,k_i) ≠ 0 and |α⃗| = u is α⃗ = k⃗ itself, so that slice of U *is* T: **U^{(l⃗)}_{k⃗} = T^{(l⃗)}_{k⃗} when |k⃗| = u**. Hence a vanished Gram block kills every total-u moment outright; no transform inversion is needed.

7. **Cascade induction.** Claim: all d-leg moments at odd total u ≥ ℓ vanish whenever 2u + d < m\*. Induct on u. At the rung m = 2u + d (odd, < m\*), the surviving (|α⃗|, |β⃗|) splits are pairs of odd totals ≥ ℓ summing to 2u; every unequal split has its smaller total v < u with 2v + d < m\*, dead by inductive hypothesis; so P_{m,d} is the equal-total Pascal-Gram form alone. But p_m ≡ 0 below m\*, so P_{m,d} = 0, so every U vanishes, so by step 6 every total-u moment vanishes. ∎(claim)

   At m\* itself, for each d: if u = (m\* − d)/2 is odd and ≥ ℓ, the cross splits are dead by the cascade and P_{m\*,d} is the equal-total form, ≥ 0. Otherwise (u even, or u < ℓ) *every* split has a total killed by parity, girth, or the cascade, and P_{m\*,d} = 0 exactly. ∎

**The selection rule (monomiality derived where it holds).** The surviving classes at m\* are exactly those d with u = (m\* − d)/2 odd and ≥ ℓ, i.e. **d ≡ m\* − 2 (mod 4) and d ≤ m\* − 2ℓ**. For deg := m\* − 2ℓ ∈ {1, 3} this set is the single class d = deg: at the first two rungs of the ladder the moment is a positive **monomial**, now derived rather than residual (every k = 3 case and every k = 4 deg ≤ 3 case in the censuses is single for this reason). From deg = 5 on, two classes may coexist (at m\* = 2ℓ+5 the rule allows d ∈ {1, 5}); the γ⁵ witness IIXY+ZXZY is in fact single because its t₅ also vanishes, but nothing forces that, and nothing needs to: positivity of all classes is what the converse uses. The old k = 3 face of R-sign, "P_{2ℓ+3,3} > 0 for genuine cycles" with its [F103 §7.5](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) Perron-skew reading, is now the d = 3 instance of the theorem.

**Verification.** [`simulations/f87_pascal_gram_positivity.py`](../../simulations/f87_pascal_gram_positivity.py) checks every step against the exact CRT machinery: the cyclic decomposition and the leg factorization term by term (diff exactly 0), the Hermitian conjugacy (exactly 0), the Pascal-Gram value against the exact first coefficient on all five branch representatives, d = 1 (IXXZ+XIXZ, 573440), d = 3 (K3 2064384, flux 589824, multi-Z 61440), d = 5 (IIXY+ZXZY, 86507520), the cascade's forced zeros (t₃ = t₅ = 0 and all total-3 triple moments exactly 0 for the γ⁵ witness, while K3's total-3 triple moments are nonzero exactly where its γ³ rung fires), the slice inversion, and the selection rule on every representative.

**The structural law, complete.** The first nonvanishing odd moment sits at

  **m\* = 2ℓ + deg,   deg odd,   deg = 1 ⟺ the girth moment fires (t_ℓ ≠ 0),   all surviving classes ≥ 0.**

At k = 3 the cell realizes only deg ∈ {1, 3}, with the deg = 1 channel exactly the single-site-Z lift (t₁ = 2^N c_l); at k = 4 the higher girth moments join the deg = 1 channel (the B-monomial production of §4) and the γ⁵ rung appears. The earlier **m\* = 3ℓ** guess remains a low-ℓ coincidence, false at ℓ = 5 where m\* = 13, not 15.

## §6 Consequence

Put the pieces together. For a non-bipartite windowed diagonal-cell pair the existence of m\* is the proven genericity result, and the Pascal-Gram positivity theorem (§5) makes every γ-coefficient of p_{m\*} a sum of squares or exactly zero, with at least one positive. A polynomial with non-negative coefficients, not all zero, has no positive real root, so

  **p_{m\*}(γ) > 0 for all γ > 0   ⟹   spec(M) ≠ spec(−M) for all γ > 0   ⟹   non-bipartite is hard at every operating point.**

This is the upgrade. [F103 §7.6](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) closed the windowed converse to "hard for all but finitely many γ" by a degenerate-perturbation-theory plus analyticity argument: a first-order break (c ≠ 0) forces some recentred characteristic-polynomial coefficient Δ_j(γ) to be a nonzero polynomial, whose zero-set is finite. The residual lemma there, the part flagged Tier1Candidate, was exactly "no positive γ is one of those finitely-many accidental soft points," a resultant / Sturm question on Δ_j. The positivity theorem answers it structurally: the first nonvanishing odd moment has non-negative coefficients throughout, hence is nonzero at every γ > 0, so the spectrum is asymmetric at every operating point and there are no accidental soft points left to rule out. "All but finitely many γ" becomes "all γ > 0," **unconditionally on both branches** (until 2026-06-10 this read "modulo R-deg + R-sign"; the girth dichotomy retired R-deg in the morning wave and the Pascal-Gram positivity theorem resolved R-sign the same day).

The result is carried by two typed claims and three scripts. The Tier1Derived spine, the threshold and the all-odd word parity and the soft re-proof and the deg-1 closed forms (the c_l² identity and its general-m sum-of-squares form, §4), is the node `WindowedConverseThresholdClaim`. The full lemma, now Tier1Derived with no residual, is the node `WindowedConverseAllGammaClaim` (the residual [F103](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) isolated was sharpened three times: from a 700-point numerical spot-check to a polynomial-root statement, from two residuals to one, and finally to none). The reproducible anchors are [`simulations/f87_windowed_monomial_converse.py`](../../simulations/f87_windowed_monomial_converse.py), [`simulations/f87_girth_dichotomy.py`](../../simulations/f87_girth_dichotomy.py), and [`simulations/f87_pascal_gram_positivity.py`](../../simulations/f87_pascal_gram_positivity.py), self-validating: every block raises on failure, prints a single PASS line on success, and the process exits 0 only if the whole ledger holds.
