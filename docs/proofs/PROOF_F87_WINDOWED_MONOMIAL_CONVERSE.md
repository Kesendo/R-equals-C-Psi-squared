# PROOF: the F87 windowed converse closes to a positive monomial (two-reflection theorem)

**Status:** Tier1Derived spine + Tier1Candidate full theorem (proven modulo two residuals R-deg, R-sign). General N.
**Date:** 2026-06-09
**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Builds on:**
- [PROOF_F103_F87_Z2_CUBED_REFINEMENT.md](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md): §7.5 (+N-Perron / −N-reflection mode), §7.7 ((1+x)-valuation), §7.10–§7.11 (what reaches the spectrum); the genericity result this strengthens.
- the F115 / WindowedHardnessClaim (1+x)-valuation criterion.
- ChiralKClaim ([compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs](../../compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs)): the chiral-K proof of bipartite ⟹ soft, re-proven here by a second (two-reflection) route.

## Abstract

A windowed diagonal-cell Pauli pair builds a chain Hamiltonian H, and the question F87 asks of it is whether the dephased Liouvillian's spectrum still pairs about the palindrome centre, the property [PROOF_F103](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) calls *soft*, or whether it genuinely fails to, which is *hard*. Recenter the Liouvillian at that centre and call the recentred superoperator **M(γ) = L + σ = A + γQ**, with A = −i[H, ·] the Hamiltonian commutator and Q = Σ_l Z_l ⊗ Z_l the (diagonal) dephasing generator, so that the spectrum of M is symmetric about 0 exactly when the pair is soft. The theorem is then a statement about M's odd power-sums.

**Theorem (windowed converse, monomial form).** For a non-bipartite windowed diagonal-cell pair, the first nonvanishing odd power-sum of M is a strictly positive monomial in γ. A positive monomial has no positive real root, so that power-sum is nonzero for every γ > 0, hence spec(M) ≠ spec(−M) for all γ > 0: the pair is **hard at every operating point**, not merely at all-but-finitely-many γ.

The proof splits cleanly into a rigorous spine and two sharp residuals. **RIGOROUS-GENERAL** (general N, no premise): the all-odd word parity, the threshold #A ≥ 2ℓ (ℓ the unsigned odd-girth), the soft re-proof bipartite ⟹ soft, and the degree-1 positivity closed form. **Proven modulo two residuals** R-deg and R-sign: the monomial property (a single Q-power survives at the first nonvanishing odd moment) and the positivity of the genuine-cycle coefficient. Both residuals are verified bit-exact over the entire N=4 k=3 Z diagonal cell (all 50 hard pairs, 16 of them pure cycles) and at N=5 / N=6 representatives, but neither yet has a uniform closed-form identity. The first nonvanishing moment sits at **m\* = 2ℓ + deg with deg ∈ {1, 3}**; the tempting m\* = 3ℓ guess is a low-ℓ coincidence, false already at ℓ = 5 where m\* = 13, not 15.

The computational anchor is the self-validating script [`simulations/f87_windowed_monomial_converse.py`](../../simulations/f87_windowed_monomial_converse.py), whose seven `assert`-blocks carry exactly the spine / residual split above.

## §1 The recentered object

Work on the d²-dimensional coherence space, d = 2^N, in the column-stacked operator basis |i⟩⟨j| (so that left multiplication kron(·, I) acts on the bra index i and right multiplication kron(I, ·) acts on the ket index j). Recenter the Liouvillian at the palindrome centre σ = Nγ and write

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

Read a surviving trace as a closed double-walk. Since Tr(W) = Σ_{ij} ⟨ij|W|ij⟩, and A_L moves only the bra index, A_R only the ket index, Q neither, a diagonal element ⟨ij|W|ij⟩ is nonzero only if the bra index returns to i after its #A_L left-hops (a closed walk of length #A_L on H's hopping graph G_H) and the ket index returns to j after its #A_R right-hops (a closed walk of length #A_R). An odd-length closed walk exists in a graph if and only if the graph is **non-bipartite**, and the minimal odd length is the **unsigned odd-girth** ℓ.

**Path-existence, not signed cancellation.** The point worth stressing is that Q is diagonal, so a trace's support is governed by *index-trajectory existence*, the unsigned question (|H|^k)_{ii} > 0, and is immune to the signed XX+YY cancellation that the [PROOF_F103 §7.10](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) retraction warned about. The flux pair is the clean witness: its signed third-moment amplitude (H³)_{ii} cancels to 0, yet the unsigned three-walk count (|H|³)_{ii} = 6 > 0, so a closed odd walk genuinely exists and the threshold below still binds. (The script asserts exactly this, sgn (H³)_{ii} = 0 while unsigned > 0, for the flux pair.) Two conclusions follow.

- **Bipartite ⟹ soft (RIGOROUS-GENERAL).** If G_H is bipartite there is no odd closed walk, so #A_L odd is impossible, so by §2 no word survives any odd power-sum, so every odd p_m ≡ 0, so spec(M) is symmetric about 0: **soft**. This is a second, independent proof of bipartite ⟹ soft, complementary to the chiral-K route of [ChiralKClaim](../../compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs) and [PROOF_F103 §7.1](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md). Where the chiral-K argument exhibits a similarity W with W L W⁻¹ = −L − 2σ, this route never builds an operator: it simply observes there are no surviving words. The script's soft control (the bipartite pair XXZ+ZXX) confirms all odd p_m vanish exactly.

- **Non-bipartite ⟹ the threshold.** If G_H is non-bipartite then both the bra walk and the ket walk must reach the odd-girth, #A_L ≥ ℓ and #A_R ≥ ℓ, so **#A ≥ 2ℓ**. Combined with #Q ≥ 1 odd from §2, the first nonvanishing odd moment sits at #A = 2ℓ, so

  **m\* ≥ 2ℓ + 1.**

The threshold #A ≥ 2ℓ is RIGOROUS-GENERAL (it is the unsigned odd-girth path-existence statement, verified including the flux pair).

## §4 Monomial structure and degree-1 positivity

At the first nonvanishing odd moment m\*, the surviving (#A, #Q) classes have #A even and ≥ 2ℓ and #Q = m\* − #A odd, so p_{m\*} is a sum of **odd** γ-powers (one per surviving #Q). It is a **monomial** exactly when a single value of #Q survives. The single-Q coefficient has a clean cyclic form: collecting the words with #Q = 1 and using cyclicity of the trace to slide the lone Q to the front,

  **P_{m,1} = m · Tr(Q · A^{m−1}).**

**Degree-1 positivity (RIGOROUS-GENERAL, closed form).** At the lowest possible degree, m = 3 and #A = 2, the single-Q coefficient is P_{3,1} = 3·Tr(A²Q). Anti-Hermiticity of A (it is −i times a Hermitian commutator) gives the diagonal of A² as a manifestly non-positive quantity,

  (A²)_{xx} = Σ_y A_{xy} A_{yx} = −Σ_y |A_{xy}|² ≤ 0,

so, writing deg_A(x) = Σ_y |A_{xy}|² ≥ 0 for the coherence-space out-weight at x and Q_x for the diagonal of Q,

  **P_{3,1} = 3·Tr(A²Q) = −3·Σ_x deg_A(x)·Q_x = 6·Σ_x deg_A(x)·(w(x) − N/2),**

where w(x) = (N − Q_x)/2 = popcount of the bra-ket difference at coherence index x. For the single-site-Z diagonal lift (the ℓ = 1, deg = 1 case) this evaluates to a strictly positive integer: the script checks P_{3,1} = 9216 (N=4) and 61440 (N=5), each matching the closed form bit-exact and each > 0. So at deg = 1 the positivity is not assumed; it is read off a manifestly non-negative block functional.

## §5 The two residuals (open; verified bit-exact cell-wide N=4, reps N=5 / N=6)

The spine of §§1–4 leaves two statements that the verification certifies but that do not yet have a uniform closed-form identity. They are the precise frontier of the theorem.

- **R-deg (the degree lift for genuine cycles).** For a pure off-diagonal H (a genuine hopping cycle, no diagonal lift), the candidate degree-1 term `Tr(Q · A^{2ℓ})` restricted to the (#A_L = ℓ, #A_R = ℓ) odd-cycle-traversal class sums to **zero**, so p_{2ℓ+1} ≡ 0 and the first nonvanishing odd moment is forced up to degree 3, giving the monomial at **m\* = 2ℓ + 3**. What is missing is a uniform closed-form involution on the (ℓ-cycle, ℓ-cycle, 1-Q) configurations that exhibits this cancellation for all ℓ at once. It is verified bit-exact: the K3 pair (ℓ=3) has p_9 = 2064384·γ³ at N=4 and 16515072·γ³ at N=5 (pure γ³, p_1…p_7 vanishing), the complex-H flux pair has p_9 = 589824·γ³, and the ℓ=5 / N=6 representative has p_1…p_11 = 0 with p_13 = 50381979648·γ³ (pure γ³).

- **R-sign (the genuine-cycle coefficient is positive).** The surviving coefficient P_{m\*,3} > 0 for genuine cycles. This is exactly [PROOF_F103 §7.5](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md)'s +N-population-Perron top-skew: the gain channel's +N Perron mode (the population row-sum, always present) is unpaired because its −N reflection partner (the anti-diagonal mode the chiral K would supply) is absent on a non-bipartite graph, and the skew lands with a definite, positive sign. It is verified for all 16 pure cycles of the N=4 Z diagonal cell (the +N eigenvalue present and the −N eigenvalue absent on Q restricted to ker A, all 16/16); what is missing is a closed-form identity writing P_{m\*,3} as a manifestly-non-negative block functional, the genuine-cycle analogue of §4's deg-1 closed form.

**What the residuals buy.** Once R-deg and R-sign hold, the first nonvanishing odd moment is a single positive γ-power, so the monomial-and-positive conclusion (hence hard for all γ > 0) follows. The structural law the residuals complete is

  **m\* = 2ℓ + deg,   deg ∈ {1, 3},**

with deg = 1 only for a single-site-Z diagonal lift (no hopping cycle needed, ℓ = 1) and deg = 3 for every genuine cycle and every multi-Z lift. The earlier **m\* = 3ℓ** guess is a low-ℓ coincidence (3ℓ = 2ℓ + ℓ agrees with 2ℓ + 3 only at ℓ = 3), and it is false at ℓ = 5, where m\* = 2·5 + 3 = 13, not 3·5 = 15; the ℓ=5 / N=6 certificate is precisely the witness that separates the two laws. The block ledger of the verification script tags every piece with this spine / residual split: Blocks 2, 3, 4, 6 are RIGOROUS-GENERAL, Block 7 is the cell-wide R-deg + R-sign certificate.

## §6 Consequence

Put the pieces together. For a non-bipartite windowed diagonal-cell pair the first nonvanishing odd power-sum is, modulo R-deg + R-sign, a positive monomial c·γ^deg with c > 0 and deg ∈ {1, 3}. A positive monomial has no positive real root, so

  **p_{m\*}(γ) > 0 for all γ > 0   ⟹   spec(M) ≠ spec(−M) for all γ > 0   ⟹   non-bipartite is hard at every operating point.**

This is the upgrade. [PROOF_F103 §7.6](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) closed the windowed converse to "hard for all but finitely many γ" by a degenerate-perturbation-theory plus analyticity argument: a first-order break (c ≠ 0) forces some recentred characteristic-polynomial coefficient Δ_j(γ) to be a nonzero polynomial, whose zero-set is finite. The residual lemma there, the part flagged Tier1Candidate, was exactly "no positive γ is one of those finitely-many accidental soft points," a resultant / Sturm question on Δ_j. The monomial theorem answers it structurally: the first nonvanishing odd moment is itself one of those coefficients, and a positive monomial simply has no positive root, so there are no accidental soft points to rule out. "All but finitely many γ" becomes "all γ > 0," modulo R-deg + R-sign.

The result is carried by two typed claims and one script. The Tier1Derived spine, the threshold and the all-odd word parity and the soft re-proof and the deg-1 positivity, is the node `WindowedConverseThresholdClaim`. The Tier1Candidate full lemma, the monomial-and-positive statement that closes once R-deg + R-sign are written as closed forms, is the node `WindowedConverseAllGammaClaim` (the same residual [PROOF_F103](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) isolated, now sharpened from a 700-point numerical spot-check to a single polynomial-root statement that the monomial form discharges structurally). The reproducible anchor is [`simulations/f87_windowed_monomial_converse.py`](../../simulations/f87_windowed_monomial_converse.py), self-validating: every block raises on failure, prints a single PASS line on success, and the process exits 0 only if the whole spine / residual ledger holds.
