# PROOF F155: The Polarity Break of the No-Jump Generator Is a Diagonal Bilinear Form

**Status:** Tier 1 derived, universal N. Every step is an identity or an exact algebraic rule, and every step is gated.
**Date:** 2026-08-19
**Authors:** Thomas Wicht, Claude (Opus 5)
**Scripts:** [`simulations/f155_polarity_break_bilinear.py`](../../simulations/f155_polarity_break_bilinear.py) (20 gates, ~5 s) for §(a) to §(f); [`simulations/f155_dephase_siblings.py`](../../simulations/f155_dephase_siblings.py) (10 gates, ~2 s) for §(g)
**Builds on:**
- [F112 non-Hermitian extension](PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md), whose section (b) parity reduction is this theorem's bit_b-even case
- [F113 coefficient derivation](PROOF_F113_COEFFICIENT_DERIVATION.md), the special case this generalises
- F38 (Π² = (−1)^bit_b on Pauli strings) and [the bit_b parity proof](PROOF_BIT_B_PARITY_SYMMETRY.md) (bit_b is the character of conjugation by X^⊗N, the global spin flip)
- **F118 and [the Π = R·D factorisation](PROOF_PI_FACTORS_AS_R_TIMES_D.md)**, F119, and F114 (typed as `compute/RCPsiSquared.Core/Symmetry/CommutatorDConjugationSign.cs`): the parents of the swap rule below, see §(d)

## What the repo already held, and what this adds

The sweep, run 2026-08-19 by agents over `docs/proofs/`, `docs/ANALYTICAL_FORMULAS.md`, `experiments/`, `hypotheses/` including `archive/`, `recovered/`, `reflections/`, `compute/` including the OpenArcs and both Confirmations registries, `docs/CAUGHT_ERRORS.md` and `simulations/framework/`:

- **Owned and proven, and this theorem contains it:** "L_σ has zero Π-conjugation ±i content whenever σ is bit_b-even", section (b) of [the F112 non-Hermitian proof](PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md), universal N, from F38.
- **The swap rule's parents are owned, and two of them are typed.** F118 factors Π_Z = R ∘ D with D the transpose superoperator, which on the Pauli basis is diag((−1)^{n_Y}); an antiautomorphism is precisely what exchanges left multiplication with right multiplication, and [the antilinear-triangle proof](PROOF_ANTILINEAR_TRIANGLE.md) says the same for the F118 edge mirrors R and 𝓕R, where conjugation by 𝒦 swaps ρ·F with F·ρ (a statement about those mirrors, not about Ad_Π on a one-sided multiplication). F114 owns the commutator combination D·L_σ·D = (−1)^{n_Y+1}·L_σ. [The factorisation proof](PROOF_PI_FACTORS_AS_R_TIMES_D.md) §(f) owns the trichotomy "Π_Z·L_σ·Π_Z⁻¹ = −L_σ, or +L_σ, or an anticommutator superoperator, neither ±L_σ". **What §(d) below adds is the per-side form with both signs, which pins the anticommutator sign that trichotomy leaves open**, and what follows from it. This parentage was missed by the session's first sweep, which asked about the consequence (the ±i content) instead of the mechanism (how Π conjugates a one-sided multiplication); a second sweep found it.
- **Owned and free:** the hypothesis-free half of the bilinear expansion of a squared-norm difference, in [the F113 proof](PROOF_F113_COEFFICIENT_DERIVATION.md). Its collapse to `4·Re⟨·,·⟩` is specific to F113's configuration and is not used here.
- **Measured only, and derived here:** the anticommutator lift's own balance, recorded 2026-08-19 in [`hypotheses/THE_DRAIN_HAS_NO_CHIRALITY.md`](../../hypotheses/THE_DRAIN_HAS_NO_CHIRALITY.md) and gate G9b of [`f112_gain_loss_carrier_check.py`](../../simulations/f112_gain_loss_carrier_check.py); and the bit_b sorting of A, recorded the same day as a one-directional measurement at N = 2, 3.
- **A near neighbour that was REFUTED and is a different object:** [`reflections/POLARITY_COORDINATES.md`](../../reflections/POLARITY_COORDINATES.md) conjectured that "any operator of the form A ⊗ B* + B ⊗ A* is automatically balanced" and killed it at 240 of 240 configurations. That is the CONJUGATE pairing; the lifts here are the TRANSPOSE pairing. The conjugate arm is the positive control of gate G9b in the sibling checker, where it is macroscopic.
- **Nothing found:** `fw.Confirmations` and the C# `ConfirmationsRegistry` (no hardware confirmation touches this quantity), `hypotheses/archive/`, `recovered/`.

## Abstract

[F112](PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md) says when the polarity asymmetry of a Lindbladian vanishes; [F113](PROOF_F113_COEFFICIENT_DERIVATION.md) gives its magnitude when a single-site Z-drive meets amplitude damping. Both leave open the NO-JUMP generator of a general non-Hermitian Hamiltonian, ρ ↦ −i(Hρ − ρH†), which is the conditional generator of post-selected dynamics and the operator PT and gain-loss models are usually written with. This proof closes that case for arbitrary H, universal N.

Write H = A + iB with A and B Hermitian and expand both over Pauli strings, A = Σ_σ a_σ σ and B = Σ_σ b_σ σ with real coefficients. Then the asymmetry is a **bilinear form, diagonal in the Pauli basis, supported on the bit_b-odd strings only**:

    asymmetry = 4^(N+1) · Σ_{σ : bit_b(σ) odd} (−1)^(#Z(σ)) · a_σ · b_σ

No pair of distinct strings contributes. F113 is the Z-diagonal special case of this, coefficient and sign included; the balance of the commutator of an arbitrary matrix falls out in one line, as the aggregate face of F112's per-pair theorem rather than that theorem; and the drain's freedom from chirality, measured that morning, is one of the two self-terms that drop because they are real.

The engine is one rule about the mirror, sharpening what F118's Π = R·D factorisation already implies: **Π-conjugation exchanges the two multiplications**, paying a per-letter sign. The two signs agree exactly on bit_b-even strings, which is why those are silent, and are opposite on the odd ones, where their difference is the coefficient above.

## (a) Statement

Let N ≥ 1, let Π be the F1 palindrome operator in the Z-dephasing convention (`framework.symmetry.build_pi_full`, F118's Π_Z: ρ ↦ ρᵀ·F), and let Ad_Π(M) = Π M Π⁻¹ on superoperators. For a superoperator M write

    asymmetry(M) := ‖P_{+i} M‖² − ‖P_{−i} M‖²

with P_{±i} the eigenprojections of Ad_Π, which has order 4. This is what F112 and F113 read and what `polarity_coordinates_from_L` returns.

**Theorem (F155).** For H = A + iB with A, B Hermitian, and G_H : ρ ↦ −i(Hρ − ρH†) = −i[A, ρ] + {B, ρ},

    asymmetry(G_H) = 4^(N+1) · Σ_{σ : bit_b(σ) odd} (−1)^(#Z(σ)) · a_σ · b_σ ,

where bit_b(σ) = #Y(σ) + #Z(σ) mod 2 and a_σ, b_σ are the Pauli coefficients of A and B.

**Which pairing, and it is not decoration.** The framework reads a row-stack generator against the order='F' Pauli transform. That pairing conjugates the bare multiplication superoperators by D = diag((−1)^{#Y}), F118's transpose superoperator, which is the same as conjugating Π. F113's registered numbers live in this pairing and so does the statement above. In the untwisted pairing the swap rule's two signs exchange and the law carries (−1)^{#Y}, so the whole value flips sign on the bit_b-odd support. Both spellings are gated (G9). **Convention-free:** the magnitude, the zero-versus-nonzero verdict, and the ratio of signs between any two odd strings. **Not convention-free:** the attachment to #Z rather than to #Y. [F113's entry](../ANALYTICAL_FORMULAS.md) pins the same freedom for its own number, and neither Π nor conj(Π) is privileged.

**What bit_b and #Z are, in the repo's own coordinates.** bit_b is the character of conjugation by X^⊗N, the global spin flip ([the bit_b parity proof](PROOF_BIT_B_PARITY_SYMMETRY.md), `docs/GLOSSARY.md` under the polarity cube). On the bit_b-odd support (−1)^{#Z} = −(−1)^{#Y}, and (−1)^{#Y} is the transpose parity, the cube's third axis. So the law reads: **only the spin-flip-odd content of A and B pairs, each term weighted by its transpose parity**, with a minus in the pairing above and a plus in the untwisted one, which is the same single global sign the paragraph before this one is about.

## (b) Step 1: the asymmetry is one imaginary part

**Lemma 1.** For every superoperator M, asymmetry(M) = Im⟨M, Ad_Π M⟩, with ⟨X, Y⟩ = Σ conj(X_ij) Y_ij.

*Proof.* Ad_Π has order 4 and is unitary for the Frobenius inner product (Π is a unitary signed permutation, so ⟨ΠXΠ⁻¹, ΠYΠ⁻¹⟩ = ⟨X, Y⟩). Its eigenprojections are P_λ = ¼ Σ_{k=0}^{3} λ^(−k) Ad_Π^k, so

    P_{+i} − P_{−i} = ½(−i·Ad_Π + i·Ad_Π³) = (i/2)(Ad_Π³ − Ad_Π) .

The P are orthogonal projections, so ‖P_{+i}M‖² − ‖P_{−i}M‖² = ⟨M, (P_{+i} − P_{−i})M⟩; unitarity gives ⟨M, Ad_Π³M⟩ = conj⟨M, Ad_Π M⟩; and with z := ⟨M, Ad_Π M⟩ the bracket is (i/2)(conj(z) − z) = Im z.  ∎

No hypothesis on M. Gate G1 checks it against the framework's own returned value on random superoperators at N = 1, 2, 3.

## (c) Step 2: everything factorises over sites

Π_N = π^{⊗N} with π the single-site map I ↦ X, X ↦ I, Z ↦ i·Y, Y ↦ i·Z. This holds by construction rather than by measurement: `build_pi_full` walks the sites, applies the per-letter action to each and multiplies the phases, which is the tensor power written out. Gate G2 confirms it bit-exactly at N = 2, 3, 4 as a wiring check. Ad_Π therefore factorises too. Writing 𝕃_σ and ℝ_σ for the two multiplication superoperators, both factorise as well, and

    −i[σ, ·] = −i(𝕃_σ − ℝ_σ) ,   {σ, ·} = 𝕃_σ + ℝ_σ ,

so G_H = Σ_σ a_σ·(−i)(𝕃_σ − ℝ_σ) + Σ_σ b_σ·(𝕃_σ + ℝ_σ) is a real-coefficient combination of the 2·4^N lifts.

## (d) Step 3: the mirror exchanges the two sides

**Lemma 2 (the swap rule).** In the pairing named in §(a), per single-qubit letter and bit-exact,

    Ad_π(𝕃_p) = ε_L(p)·ℝ_p ,   Ad_π(ℝ_p) = ε_R(p)·𝕃_p ,

with (ε_L, ε_R) = (+1, +1) for I and X, (−1, +1) for Y, (+1, −1) for Z; hence for a string σ,

    Ad_Π(𝕃_σ) = (−1)^(#Y(σ))·ℝ_σ ,   Ad_Π(ℝ_σ) = (−1)^(#Z(σ))·𝕃_σ .

In the untwisted pairing the two exponents exchange.

*Proof.* The per-site statement is a finite check on four letters and two sides of 4×4 signed-permutation matrices, bit-exact (gate G3a); the string statement is its tensor product by Step 2 (gate G3b, entry-wise deviation exactly 0.0 at N = 2 and 3); the untwisted spelling is gate G9.  ∎

**Where this comes from, and what is new.** That Π-conjugation should exchange the two sides is not a surprise and not new: F118 factors Π_Z = R ∘ D with D the transpose, and an antiautomorphism exchanges left with right multiplication by construction; [the antilinear-triangle proof](PROOF_ANTILINEAR_TRIANGLE.md) records the sibling swap ρ·F ↔ F·ρ under 𝒦 for F118's edge mirrors; F114 owns the commutator combination D·L_σ·D = (−1)^{n_Y+1}·L_σ; and [the factorisation proof](PROOF_PI_FACTORS_AS_R_TIMES_D.md) §(f) already records that in the mixed-parity cells Π_Z·L_σ·Π_Z⁻¹ lands on "an anticommutator superoperator, neither ±L_σ". **What Lemma 2 adds is the per-side form carrying both signs**, which pins the sign that trichotomy leaves open and makes the anticommutator lift computable rather than merely named. Applying it twice gives ε_L·ε_R = (−1)^{bit_b} per site, the sign pattern of F38 (recovered as a pattern; F38's statement about Π² itself is the parent, not the child).

## (e) Step 4: orthogonality, and the assembly

**Lemma 3 (lift orthogonality).** ⟨𝕃_σ, 𝕃_τ⟩ = ⟨ℝ_σ, ℝ_τ⟩ = 4^N·δ_{στ}, and ⟨𝕃_σ, ℝ_τ⟩ = Tr(σ)·Tr(τ), which vanishes unless both strings are the identity.

*Proof.* As superoperator traces, Tr(𝕃_A† 𝕃_B) = 2^N·Tr(A†B) and Tr(𝕃_A† ℝ_B) = Tr(A†)·Tr(B); for Pauli strings Tr(σ†τ) = 2^N δ_{στ} and Tr(σ) = 2^N δ_{σ,I}. Gate G3c checks all 256 string pairs at N = 2 against the integers, exactly.  ∎

Expand z = ⟨G_H, Ad_Π G_H⟩ with Step 2, apply Lemma 2, and keep what Lemma 3 leaves. Writing E_Y = (−1)^{#Y} and E_Z = (−1)^{#Z}, the same-string terms are

- **commutator with commutator:** conj(−i a_σ)(−i a_σ)·[⟨𝕃_σ, −E_Z 𝕃_σ⟩ + ⟨−ℝ_σ, E_Y ℝ_σ⟩] = −a_σ²·4^N·(E_Z + E_Y), **real**;
- **anticommutator with anticommutator:** b_σ²·4^N·(E_Z + E_Y), **real**;
- **the two cross terms, summed:** 2i·a_σ b_σ·4^N·(E_Z − E_Y).

Two bookkeeping notes that change nothing in the result. The identity string is bit_b-even and contributes nothing, but it is the one place where Lemma 3's L-R overlap does not vanish (⟨𝕃_I, ℝ_I⟩ = 4^N) and where 𝕃_I = ℝ_I makes the commutator lift identically zero; both facts move real terms only. And distinct strings drop entirely, which is the diagonality (gate G4: 210 pairs at N = 2 and 3906 at N = 3, all zero).

Taking imaginary parts kills both self-terms. That is not a small remark: the first IS F112's balance for the commutator part, and the second is the anticommutator part's balance, which until that morning was measured and not derived. What survives is

    asymmetry = Im z = Σ_σ 2·a_σ b_σ·4^N·(E_Z(σ) − E_Y(σ)) ,

which is 0 when bit_b(σ) is even (the two signs agree) and 4^(N+1)·(−1)^{#Z}·a_σ b_σ when it is odd (they are opposite).  ∎

## (f) Corollaries

1. **B = 0, i.e. Hermitian H.** Balanced: a bilinear form needs both arguments (gate G6).
2. **A carries only bit_b-even strings.** Balanced against every B (G6).
3. **B ∝ I.** Balanced against every A: the identity is bit_b-even (G6).
4. **Disjoint odd supports.** If the bit_b-odd content of A and of B sits on disjoint strings, balanced, by diagonality (G6).
5. **The sorting is one-directional.** bit_b-odd content in A is necessary and NOT sufficient; corollaries 3 and 4 are the counterexamples. This is the sharp form of the one-directional statement the F112 scope note recorded as measured on 2026-08-19.
6. **Identity shifts are invisible.** A ↦ A + cI and B ↦ B + cI leave the value unchanged (G10). Physically: the law is **blind to the gain/loss offset**, so balanced gain/loss and pure loss of the same imbalance read the same number.
7. **The textbook PT dimer is balanced at every coupling.** For H = J(XX + YY)/2 + i·g·(Z₀ − Z₁)/2 the value is exactly 0 at g = 0, 0.5, 0.99, 1.0, 1.5, 3.0, because XX and YY are bit_b-even; g = J = 1 is where this dimer's single-excitation block sits at its exceptional point (its two eigenvalues coalesce at 0 and go imaginary just above), and the asymmetry does not notice. Add an ordinary 0.4·Z₀ detuning and the value is exactly −12.8·g, a straight line with no feature at g = 1 (G11). Note that the detuned system has no exceptional point on the real g axis at all: its single-excitation gap runs 1.99, 1.80, 1.81, 2.55 across g = 0.5 … 1.5 and never closes. **The quantity is a detuning × gain-imbalance meter, not a PT meter**, which is the same lesson `docs/CAUGHT_ERRORS.md` records for the class-shaped claims this replaces.
8. **The transpose-paired lift is balanced for ANY K.** For K not Hermitian, split K = K₁ + iK₂ into Hermitian parts; the self-terms are real as above and the two cross terms are symmetric in the pair, so they cancel (G12). This makes the drain's chirality-freedom a theorem for every collapse operator, c†c being Hermitian in any case.
9. **The commutator of an arbitrary matrix is balanced.** For H = A + iB it is −i(𝕃−ℝ)_A + (𝕃−ℝ)_B, and its two cross terms cancel against each other identically (G12). This is the aggregate face of F112's non-Hermitian extension, not that theorem: F112 proves the PER-PAIR identity F(σ_α, σ_β) = 0 for every pair of strings, and diagonality means the cross-pair information is exactly what drops out here, so this route cannot reach it.
10. **F113.** With a_{Z_l} = ω_l/2 and b_{Z_l} = (γ_T1,l − γ_pump,l)/4, and #Z = 1 on every Z_l so the sign is −1,

        asymmetry = 4^(N+1)·Σ_l (−1)·(ω_l/2)·(γ_T1,l − γ_pump,l)/4 = (4^N/2)·Σ_l ω_l·(γ_pump,l − γ_T1,l) ,

    which is F113 exactly (G7, random per-site rates). **Scope of that identification:** F113's own theorem is stated for the FULL Lindbladian; this is its no-jump face. That the recycling term Σ_k c_k ⊗ c_k* leaves the value unchanged is a separate MEASURED fact (F113's entry, N = 2, 3, 4, σ⁻/σ⁺ family), not derived here, and it is false outside that family.

## (g) The dephase siblings: the X and Y conventions

*(Added 2026-08-19, the same day as the theorem above, closing the OpenArc `polarity_break_dephase_siblings`. Gate: [`simulations/f155_dephase_siblings.py`](../../simulations/f155_dephase_siblings.py), 10 gates, ~2 s.)*

Π above is the Z-dephasing palindrome. The X- and Y-dephasing conventions Π_X, Π_Y are the other two members of the Klein triple, and they have been available as a parameter since long before this arc (`framework.symmetry.pi_action`, `compute/RCPsiSquared.Core/Symmetry/PiOperator.cs`), and callers have been passing letters into `polarity_coordinates_from_L` for months: [`simulations/f112_klein_v4_cross_dephase_verify.py`](../../simulations/f112_klein_v4_cross_dephase_verify.py), the verifier of the Klein-V₄ proof cited twice below, runs X and Y and already implements the χ_ℓ half of the theorem below. What had never been asked in a non-Z convention is the **value** of the asymmetry for this generator; that script asks whether it **vanishes**, and for the full Lindbladian. (An earlier draft of this section claimed no caller had ever passed a letter, with the counterexample sitting in a file it linked. The narrower gap is real and is what the gate script needs: `polarity_coordinates` hard-wires the Z operator in its chain-facing entry points, and only `polarity_coordinates_from_L` accepts an injected Π.) **The answer is not two new laws.** Both siblings reduce to the theorem above, by two different mechanisms, and only one of them is even a change of law.

**Theorem (the letter as a parameter).** For dephasing letter ℓ ∈ {X, Y, Z},

    asymmetry_ℓ(G_H) = 4^(N+1) · Σ_{σ : χ_ℓ(σ) odd} (−1)^(#ℓ(σ)) · a_σ · b_σ ,

where χ_ℓ is the parity Π_ℓ² grades by (F38/F88a: bit_a for ℓ = X, bit_b for ℓ = Y and Z, which `pi_squared_eigenvalue` already carries) and #ℓ counts the occurrences of the dephasing letter **itself**. ℓ = Z returns the theorem of §(a).

**Which half of that is frame-free, because "the letter itself" is pretty and only half true.** The **support** χ_ℓ is convention-free. The **weight** is not: it is the dephasing letter only in the twisted pairing §(a) names. In the untwisted pairing the two per-letter signs exchange and the weights become (−1)^#Y for **both** Z and X and (−1)^#Z for Y, so there the weight letter is the dephasing letter for none of them (gate S3 measures the exchange in all three conventions). This is exactly the freedom §(a) prices for its own law, inherited unchanged, and the same three quantities stay convention-free: the magnitude, the zero-versus-nonzero verdict, and the ratio of signs between any two strings in the support.

*Proof.* Steps 1, 2 and 4 never mention the letter: Lemma 1 needs only that Ad_Π is order-4 and Frobenius-unitary, Π_ℓ = π_ℓ^⊗N holds by construction for all three letters, and Lemma 3 is Π-free. Only Lemma 2 is letter-specific, and recomputing it gives the three tables (gate S2, bit-exact):

| convention | I | X | Y | Z | support | weight |
|---|---|---|---|---|---|---|
| Π_Z | (+,+) | (+,+) | (−,+) | (+,−) | bit_b-odd | (−1)^#Z |
| Π_X | (+,+) | (+,−) | (−,+) | (+,+) | bit_a-odd | (−1)^#X |
| Π_Y | (+,+) | (+,+) | (+,−) | (−,+) | bit_b-odd | (−1)^#Y |

The support is where ε_L and ε_R disagree, and to see that this is generic rather than a fit to three rows, redo §(e) in letter-generic symbols: with s_L(σ) = Π ε_L, s_R(σ) = Π ε_R and c_σ = −i a_σ + b_σ, Lemma 3 leaves z = 4^N Σ_σ [conj(c_σ)²·s_R(σ) + c_σ²·s_L(σ)], whose imaginary part is 4^N Σ_σ 2 a_σ b_σ (s_R(σ) − s_L(σ)). Only the **difference** of the two signs survives, for any letter. So the support is the disagreement set and the weight is the letter carrying ε_R = −1, both read off the table. For ℓ = Z, s_L = (−1)^#Y and s_R = (−1)^#Z, which is §(e) as written there. ∎

**Y is not a sibling.** Π_Y = Π_Z⁻¹ exactly, at every N, which the repo already owns ([the factorisation proof](PROOF_PI_FACTORS_AS_R_TIMES_D.md), [the Klein-V₄ dephase-swap proof](PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md), `reflections/D_PI_Z_EQUALS_PI_Y.md`). Hence Ad_{Π_Y} = Ad_{Π_Z}³, and inverting an order-4 map exchanges its +i and −i eigenprojections. So for **every** superoperator M, with no hypothesis and no bilinear argument at all,

    asymmetry_Y(M) = −asymmetry_Z(M) .

Gate S4 checks this **bit-exactly** on fifteen arbitrary complex superoperators, none of generator form: the two squared norms are the same two floats subtracted in the opposite order. The Y row of the table above is that identity restricted to G_H, and the global minus is the same single sign the twisted-versus-untwisted pairing freedom of §(a) already carries. What this adds to the owned Π_Y = Π_Z⁻¹ is that the fact alone decides the asymmetry, without touching the generator.

**X is the one genuine sibling, and it is this theorem in a rotated frame.** A second and independent route reaches the X law without recomputing anything: for h = H^⊗N, the **same** Hadamard on **every** site (the X ↔ Z relabelling),

    asymmetry_X(A, B) = asymmetry_Z(hAh, hBh) ,

with no sign and no scale correction (gate S10, at the noise floor). It is a **route** and not a numerical coincidence, because Q_zx·Π_Z·Q_zx⁻¹ = Π_X with Q_zx the Hadamard conjugation, which §(d) of [the Klein-V₄ cross-dephase proof](PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md) already states, and Ad_{Q_zx} is Frobenius-unitary, so the two readings differ only in which frame the Pauli content is expressed. **Every site must be rotated:** a Hadamard on a proper subset gives a mixed-letter law, not the X one, and disagrees macroscopically (gate S10 checks that too). This is the same Hadamard transport that proof uses for the F112 siblings, and it completes gate G13 of the parent script: G13 measured that the value **moves** under a per-site Hadamard, and S10 says **where it moves to**, namely onto the X-convention reading of the unrotated system.

**The consequence worth quoting, and it is a SILENCE rather than a balance.** The configuration that makes this theorem subsume F113, a per-site Z-drive meeting amplitude damping, reads exactly 0.0 under Π_X, and the two ±i norms are **individually** and exactly zero (gate S8): there is no cancellation, there is nothing there. That is the repo's third answer beside balanced and broken, the one `is_structurally_degenerate` names, and reading "balanced" off it would confirm nothing.

It is also more general than F113. If every Pauli string in A and B is bit_a-**even**, then G_H lies in the +1 eigenspace of Ad_{Π_X}², so its anti part is the exact zero array and both halves vanish, at **every** N and every rate. Gate S8 checks that general statement on a bit_a-even but non-diagonal configuration (a YY bond against a single Z), not only on the F113 family. So the asymmetry is a property of the physics **paired with a palindrome convention**, which is the sharpest form of the frame warning in the next section.

**What this corrects.** The scope note this section replaces, its copy in [the F-registry](../ANALYTICAL_FORMULAS.md), the OpenArc that carried this item, and the parent gate script's own SCOPE block all described the siblings as arising "under the axis substitution bit_b ↔ bit_a", which lumps Y in with X. Their own cited source says the opposite (the Klein-V₄ cross-dephase proof: F112-Y shares the bit_b axis with Z and comes from the other route), and the table above settles it: **Y shares Z's support exactly and differs only by a global sign; X is the only one that moves the axis.** By D·Π_Z·D = Π_Y, which that same proof owns, the Y minus is not merely the analogue of the pairing freedom but literally the same operator D acting. Separately, gate S5 pins Π_ℓ² to one parity per letter, refuting the sentence "Π²_Y activates bit_a + bit_b simultaneously" that stood in the typed claim `compute/RCPsiSquared.Diagnostics/F87/DissipatorAxisSelectsPolarity.cs`, which had confused the Klein index of the **letter** Y, which is indeed (1,1), with the grading of Π_Y², which is bit_b alone. The claim did not invent that conflation: its anchor `hypotheses/THE_POLARITY_LAYER.md` carried it in the same words, in a list explicitly of polarity axes, and from there it reached a second typed claim and a live arc. All of them are repaired in the same change, and `docs/CAUGHT_ERRORS.md` records the episode.

**Still not claimed.** The X- and Y-dephase break of the **full Lindbladian**, recycling term included, is a different object and remains open; the recycling term's inertness is measured only for the σ⁻/σ⁺ family and is false outside it.

## (h) Scope, and what is not claimed

- **The palindrome convention of §(a) to §(f) is Z**, and §(g) closes the other two: Y is this theorem with a global minus (an exact consequence of Π_Y = Π_Z⁻¹, universal in M), and X is this theorem on the bit_a-odd support with weight (−1)^#X, equivalently this theorem in the Hadamard-rotated frame.
- **The value is not basis-invariant, and that is the scope fact a reader most needs.** It moves under a per-site Hadamard (an X ↔ Z relabelling) and under a local rotation (gate G13). The law reads the Pauli content of A and B in the frame the Z-palindrome convention fixes, so it must be quoted in that frame, exactly as F113's entry says its number is a lab-frame reading.
- **Which generator.** −i(Hρ − ρH†) is the un-normalised conditional generator of post-selected dynamics, and it is the operator PT and gain-loss models are usually written with. It is NOT the full Lindbladian: a physical gain medium pumps incoherently, and the recycling term it adds changes the asymmetry outside the σ⁻/σ⁺ family. It is also not trace-preserving; the state-dependent normalisation of post-selection is nonlinear and is not part of this object. Corollary 6 is what makes the un-normalised reading defensible, since the state-independent part of that normalisation is exactly an identity shift.
- **A and B Hermitian** is no restriction on H: every matrix splits that way, uniquely.
- **Nothing about spectra.** No claim about eigenvalues, dynamics, or PT-breaking thresholds; corollary 7 shows the quantity is blind to the threshold.

## (i) Gates

`simulations/f155_polarity_break_bilinear.py`, 20 gate instances, all passing, about 5 seconds. G1 Lemma 1; G2 the factorisation of Π; G3a and G3b the swap rule, per site and per string; G3c lift orthogonality; G4 diagonality; G5 the closed form on dense random Hermitian A and B at N = 2, 3, 4; G6 corollaries 1 to 4; G7 the F113 reduction; G8 a positive control; G9 the two pairings and their #Y ↔ #Z mirroring; G10 identity shifts; G11 the PT dimer through its exceptional point; G12 corollaries 8 and 9; G13 the basis dependence.

[`simulations/f155_dephase_siblings.py`](../../simulations/f155_dephase_siblings.py), 10 gate instances, all passing, about 2 seconds, for §(g). S1 Π_Z unitary, of order 4 and symmetric, with Π_Y its inverse entry for entry; S2 the three swap tables, measured and checked by exact set equality against the law as written; S3 the string form per convention and the untwisted spelling; S4 asymmetry_Y = −asymmetry_Z bit-exactly on arbitrary superoperators; S5 the Π² grading, one parity per letter; S6 the closed form across nine decades of input scale at N = 2, 3, 4; S7 diagonality, exactly; S8 the F113 reduction and its silence under Π_X, with B built from the collapse operators; S9 the positive controls with their off-support halves; S10 the Hadamard transport and its failure under a partial rotation.

The exact steps (G2, G3a, G3b, G3c, and S1 to S5 and S7) are compared exactly, because an exact route exists. The float comparisons are gated against eps·‖M‖²_F; the parent takes that error model from `f112_gain_loss_carrier_check.py`, and the sibling script states it on its own page and gates it as a law, requiring the ratio to stay put across nine decades rather than merely to fit under a ceiling.
