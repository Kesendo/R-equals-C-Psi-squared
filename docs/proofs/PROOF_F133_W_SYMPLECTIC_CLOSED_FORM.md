# PROOF F133: the symplectic closed form of the cofactor W

*2026-07-17. Resolves the open lead of [PROOF_F128_FLIP_SUM_FACTORIZATION](PROOF_F128_FLIP_SUM_FACTORIZATION.md) §5, which had left the cofactor W as the whole content of 𝔉 up to the explicit prefactor and named W's own closed form the natural next question. Consumes F128 and the classical Weyl denominator formula for Sp(12); nothing else.*

Before the machinery: what this feels like. F128 left one opaque object standing, the odd-projected cotangent bracket W, and everything about the cross form that was not the visible prefactor lived inside it. Today the inside is written out: six sines, two cosine Vandermondes, and one symmetric polynomial, a sum of 143 Weyl characters of Sp(12) whose integer coefficients never exceed 8, divided by the product of the 32 sheet sines. The 64 sheet vectors that have haunted this arc since the residue collapse turn out to be the weight system of the spin representation of so(13), and the characters that tame them live not in so(13) but in its Langlands dual. The cross form has no secrets left at the level of shape; the one thing still open is the law behind the 143 integers.

## §1 Definitions and the statement

All functions live on the six angles x = (a₁, a₂, a₃, b₁, b₂, b₃); write x_u, u = 1..6, for the flat list and c_u = cos x_u. From F127/F128 (definitions verbatim in [PROOF_F128_FLIP_SUM_FACTORIZATION §§1, 3](PROOF_F128_FLIP_SUM_FACTORIZATION.md); P̃ is §3's):

- e₁ = Σ cos a_i, f₁ = Σ cos b_j, s = (Σa + Σb)/2,
- V_a = Π_{i<j} sin((a_i − a_j)/2), V_b likewise, P = Π_{i,j} sin((a_i + b_j)/2), P̃ = Π_{i,j} sin((a_i − b_j)/2),
- 𝒪 = the projection onto the part odd in each of the six angles separately,
- 𝔉 = the committed cross form; F128: 𝔉 = −(e₁ − f₁)² · W with **W = 𝒪[cos s · cot s · V_a V_b / P]**, the cofactor this proof is about.

New objects:

- Δ̂ = Π_{u<v} sin((x_u − x_v)/2), the full six-angle half-angle Vandermonde (= V_a V_b P̃ in the ordering a₁..a₃, b₁..b₃),
- the **sheets**: L ∈ {±1}⁶ acting by (L∘x)_k = L_k x_k; canonical means L₁ = +1 (32 of them); s_L = (L·x)/2 the sheet half-angle, so s = s_{1⁶},
- **SP** = Π_{L canonical} sin(s_L), the product over all 32 canonical sheets,
- **V_c(a)** = Π_{i<j} (cos a_i − cos a_j), V_c(b) likewise; Vand_c = Π_{u<v} (c_u − c_v) the full cosine Vandermonde,
- **χ^{C₆}_λ** = the Weyl character of Sp(12) (type C₆) at z_u = e^{i x_u}, as a bialternant: with ρ = (6, 5, 4, 3, 2, 1) and μ = λ + ρ, χ^{C₆}_λ = A_μ / A_ρ where A_μ = det[z_u^{μ_j} − z_u^{−μ_j}]_{u,j}.

**Theorem (F133, the closed form).** As an identity of rational functions of the six angles,

  **W = −2⁹ · (Π_u sin x_u) · V_c(a) · V_c(b) · K(c) / SP(x),**

where K is the symmetric polynomial in (c₁, …, c₆)

  **K = 2⁻³⁰ · Σ_λ n_λ · χ^{C₆}_λ,**

a sum over exactly 143 partitions λ (≤ 6 parts, even |λ| ≤ 18, λ₁ + λ₂ ≤ 10), with integer coefficients n_λ, |n_λ| ≤ 8 (65 of the 143 are ±1, the largest single class). As a polynomial in c, K has monomial degrees 6 to 18 only; the low-degree characters (down to λ = ∅ with n_∅ = −1) cancel below degree 6. The coefficient table is committed at [simulations/results/f133_w_closed_form/chiC_coeffs.txt](../../simulations/results/f133_w_closed_form/chiC_coeffs.txt) (the same K in the monomial-symmetric basis, 190 dyadic-rational terms, at [m_coeffs.txt](../../simulations/results/f133_w_closed_form/m_coeffs.txt)).

**Corollary (the cross form, fully explicit).** With F128,

  𝔉 = +2⁹ · (e₁ − f₁)² · (Π_u sin x_u) · V_c(a) V_c(b) · K(c) / SP(x).

Every named factor is elementary; the only nontrivial content of 𝔉 is the single symmetric polynomial K, i.e. the 143 integers n_λ.

## §2 The reduction to one trigonometric polynomial

Three identity steps move W from a rational function with an odd projection outside to the odd projection of a polynomial.

**[R1] The escape.** Pairwise, sin((a_i + b_j)/2) · sin((a_i − b_j)/2) = (cos b_j − cos a_i)/2 (gate G1, exact dicts, all nine pairs; "dict" throughout = the F128 toolkit's sparse exponent-vector → integer-coefficient dictionary in the half-angle units t_u = e^{i x_u/2}, exact arithmetic over ℤ), so P·P̃ = 2⁻⁹ Π_{i,j}(cos b_j − cos a_i) is even in every angle and passes through 𝒪. Multiplying inside W by P̃/P̃ and using V_a V_b P̃ = Δ̂:

  W = 2⁹ · 𝒪[cos s · cot s · Δ̂] / Π_{i,j}(cos b_j − cos a_i).

**[R2] The sin-s lemma.** cos s · cot s = csc s − sin s, and

  **𝒪[sin s · Δ̂] ≡ 0** (gate G2, exact over ℤ: the 64-flip signed character sum annihilates the integer half-angle polynomial, zero monomials survive; the same mechanism as F128's flip lemma, one power of the sheet term lower).

So the core is the single object 𝒞 = 𝒪[Δ̂ / sin s] (the trigonometric core; the Lie type is always written C₆ or Sp(12) here, never bare), and W = 2⁹ · 𝒞 / Π(cos b_j − cos a_i).

**[R3] Clearing the pole.** SP is flip-invariant: a generator flip ε_u permutes the 64 sheet forms; re-canonicalizing (sin is odd, so sin(s_L) = −sin(s_{−L})) costs one sign per sheet whose first entry flips, and the signs multiply to +1 for every generator (for u = 1 all 32 sheets re-canonicalize, (−1)³² = +1; for u ≥ 2 none do; gate G3, integer bookkeeping). Hence SP passes through 𝒪, and since SP / sin s = Π_{L canonical, L ≠ 1⁶} sin(s_L) (31 factors),

  **𝒞 · SP = 𝒪[X], X = Δ̂ · Π_{31 sheets} sin(s_L),**

a trigonometric *polynomial* with integer coefficients in the half-angle units t_u = e^{i x_u/2}.

## §3 The alternant read-off: the 143 coefficients, derived

X is S₆-alternating: Δ̂ alternates (exact dict check under all 15 transpositions, gate G4), and the 31-sheet product is S₆-invariant (a transposition permutes the 31 non-trivial canonical sheets; the re-canonicalization signs count the canonical sheets whose first entry moves onto a −1, which is 0 or 16, always even; gate G4). Therefore 64·𝒪[X] is antisymmetric under the full hyperoctahedral Weyl group W(C₆) = S₆ ⋉ (ℤ/2)⁶.

A W(C₆)-antisymmetric integer Laurent polynomial is an integer combination of the symplectic alternants A_μ (μ strictly dominant: μ₁ > … > μ₆ > 0), and the coefficients are read off at the dominant exponents: distinct alternants have disjoint W(C₆)-orbits of exponent vectors, and each A_μ carries its own dominant monomial with coefficient exactly +1. In t-units (z = t², so A_μ sits at t-exponents ±2μ):

  **n_raw(λ) = Σ_{ε ∈ {±1}⁶} sgn(ε) · [X]_{ε∘2(λ+ρ)},**

no fit and no degree cap anywhere. The extraction is pure integer dictionary arithmetic, meet-in-the-middle: P₁ = Δ̂ · 15 sheet factors (590 016 monomials after cancellation), P₂ = the other 16 sheet factors (5 817 monomials), and each coefficient of X = P₁·P₂ is a convolution lookup. The feared 2³¹-monomial common product never materializes; cancellation keeps the halves small.

**Support is forced, then swept.** X's per-variable t-exponents lie in [−36, 36] (5 from Δ̂, 31 from the sheets; asserted on the computed dict in gate G5), so a contributing alternant needs 2(λ₁ + 6) ≤ 36, i.e. **λ₁ ≤ 12**. The complete sweep over all 18 564 dominant λ with ≤ 6 parts ≤ 12 (gate G5 with `--full`, ~9 min; the default ~25 s run checks the 557-candidate even-degree window only, the exhaustive claim needs `--full`) leaves exactly 143 nonzero values, none outside the committed table, none missing, none different. The even-degree-only support |λ| ≤ 18 and the pair cap λ₁ + λ₂ ≤ 10 are corollaries of the surviving support, not inputs (the survivors span |λ| = 0, 2, …, 18; the 6..18 degree range belongs to K as a *polynomial*, i.e. to the monomial-symmetric table, where the low-degree characters cancel). Everything satisfies n_raw(λ) = 2·n_λ; the factor 2 is absorbed into the normalization chase of §4.

The same 143 integers were found first by an independent route (GF(p) linear solve in the character basis over two 30-bit primes, 40 consistent overdetermination rows per prime, CRT + rational reconstruction); the read-off above re-derives every value. Gate G7 additionally pins the two committed tables (character basis and monomial-symmetric basis) against each other exactly at random GF(p) points.

## §4 The Weyl denominator and the assembly

The classical type-C₆ Weyl denominator formula, as an exact dict identity (gate G6, 46 080-monomial determinant expansion versus the product form):

  A_ρ = det[z_u^{7−j} − z_u^{−(7−j)}] = Π_u (z_u − z_u^{−1}) · Π_{u<v} (z_u + z_u^{−1} − z_v − z_v^{−1}) = (2i)⁶ · 2¹⁵ · (Π_u sin x_u) · Vand_c.

Now the constants, once. In dict units X carries (2i)⁴⁶ (15 Vandermonde factors + 31 sheet factors, each stored as t^v − t^{−v} = 2i sin), and the dict projector is 64·𝒪. So §3 says 64·(2i)⁴⁶·𝒪[X] = Σ_λ n_raw(λ)·A_{λ+ρ} = 2·Σ_λ n_λ·A_{λ+ρ}, i.e.

  𝒞 · SP = 𝒪[X] = 2⁻⁵·(2i)⁻⁴⁶ · Σ n_λ A_{λ+ρ}.

Dividing by (Π sin x_u)·Vand_c = A_ρ · (2i)⁻⁶·2⁻¹⁵ and using χ_λ = A_{λ+ρ}/A_ρ:

  K := 𝒞·SP / ((Π sin x_u)·Vand_c) = (2i)⁶·2¹⁵ / (2⁵·(2i)⁴⁶) · Σ n_λ χ_λ = **2⁻³⁰ · Σ n_λ χ^{C₆}_λ**

(i⁴⁰ = 1). The committed common denominator 2³⁰ of the character-basis fit is thereby derived, not observed. Finally, assemble: Vand_c = V_c(a)·V_c(b)·Π_{i,j}(c_{a_i} − c_{b_j}) and Π(c_{a_i} − c_{b_j}) = (−1)⁹·Π(cos b_j − cos a_i), so [R1]-[R3] give

  W = 2⁹·𝒞/Π(cos b_j − cos a_i) = 2⁹·(Π sin x_u)·V_c(a)V_c(b)·(−1)⁹·K/SP = −2⁹·(Π sin x_u)·V_c(a)V_c(b)·K/SP. ∎

Gate G8 corroborates the assembled identity numerically: W against the literal 64-flip evaluation of the F128 cofactor, and 𝔉 against the committed `cross_form`, rel. deviation ~10⁻¹⁰ at generic points.

## §5 The naming, and what stays open

The 64 sheet vectors (±1)⁶, read as half-angle forms L·x/2 = (±½, …, ±½)·x, are the weight system of the **spin representation of so(13) = B₆** (dimension 2⁶ = 64); the 32 canonical sheets are one choice of positive system (L₁ = +1, a set of ±-representatives), and SP is the spinor sine-product Π_{w > 0} sin⟨w, x⟩ of that representation. The character basis in which K collapses is type **C₆ = Sp(12), the Langlands dual of B₆**; the type-B₆ character basis was tried and is rank-deficient on this problem (it fails consistently), so the symplectic basis is forced, not chosen (a computed rank-deficiency, at the same code-trust grade as the rest). SP expanded in χ^{C₆} alone is not compact (1 096 terms, coefficients to ±41); the compactness lives only in the full ratio.

Open, deliberately: a closed form for the 143 integers n_λ themselves (nearest literature kin: King et al., arXiv:2303.00576, spin ↔ symplectic dual pairs with residue-periodic small-integer coefficients; Rains-Warnaar, arXiv:1506.02755, bounded Littlewood identities of the right index genre). Simple product ansätze over natural factor families are ruled out (5 316 exponent combinations, pointwise GF(p), zero constant ratios). Nothing about an n_λ law is claimed here.

## §6 Evidence grade

Every link of §§2-4 is exact: integer dict identities (G1, G2, G6), integer bookkeeping with structural bijections (G3, G4), and integer coefficient extraction with a forced, fully swept support window (G5). The 143 coefficients are derived twice by independent exact routes (GF(p)+CRT fit; alternant read-off). The one numeric link is corroboration, not load-bearing: G8 re-checks the assembled rational-function identity against the committed `cross_form` at generic points. The standing caveat is the arc's usual one: code trust (single-language implementations of each gate; the C# witness is the second implementation of the read-off core).

## §7 Verification index

| gate | proves | method |
|---|---|---|
| G1 | the P·P̃ escape, pairwise | 9 exact dict identities over ℤ |
| G2 | the sin-s lemma 𝒪[sin s·Δ̂] ≡ 0 | 64-flip signed character sum annihilates the integer polynomial |
| G3 | SP flip-invariance | re-canonicalization bookkeeping per generator flip, signs multiply to +1 |
| G4 | X is S₆-alternating | Δ̂ antisymmetric as exact dict; 31-sheet product invariant, even sign counts |
| G5 | the 143 n_λ, derived | meet-in-the-middle integer extraction; window λ₁ ≤ 12 forced by X's support; `--full` sweeps all 18 564 dominant λ |
| G6 | the C₆ Weyl denominator | 46 080-monomial determinant expansion == product form, exact over ℤ |
| G7 | the two committed tables agree | K evaluated both ways at exact GF(p) points, 2 primes |
| G8 | the assembled closed form | numeric pin vs `W_brute` and vs the committed `cross_form`, 10 generic points |

Run: `python simulations/f133_w_closed_form.py` (~25 s, exit 0 iff all pass); `--full` adds the ~9 min support sweep: [f133_w_closed_form.py](../../simulations/f133_w_closed_form.py).
