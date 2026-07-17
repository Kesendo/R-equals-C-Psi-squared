# PROOF F128: the flip-sum factorization of the cross form, and the sharper locus

*2026-07-14. The handover's lead 1 (the 32-term flip-sum's own proof), resolved the same day it was posed. Companion to [PROOF_F127_RESIDUE_COLLAPSE](PROOF_F127_RESIDUE_COLLAPSE.md), whose §3 closed form and §3A representation are the two committed inputs.*

Before the machinery: what this feels like. F127 said the cross form vanishes when two things are separately zero, and the proof needed a wall, or a residue collapse, or both. Today the same object turned out to carry its whole constraint openly, as one visible factor: the cross form is the square of the *difference* of the two cosine sums, times a bare odd-projection. The two constraints were never two. The wall computed for ninety hours what one subtraction now writes down, and the vanishing no longer needs the constraints to be zero at all, only to agree.

## §1 The statement

All functions live on the six angles x = (a₁, a₂, a₃, b₁, b₂, b₃). Write

- e₁ = cos a₁ + cos a₂ + cos a₃ and f₁ = cos b₁ + cos b₂ + cos b₃ (the two constraint sums),
- s = (Σa + Σb)/2 (the half-sum),
- V_a = Π_{i<j} sin((a_i − a_j)/2), V_b likewise in b, P = Π_{i,j} sin((a_i + b_j)/2),
- 𝒪 = the projection onto the part odd in each of the six angles separately: 𝒪[g](x) = 2⁻⁶ Σ_{ε∈{±1}⁶} (Πε) · g(ε∘x), where (ε∘x)_k = ε_k x_k,
- 𝔉 = the committed cross form of the twinning arc (`cross_form` in [cross_triple_orthogonality.py](../../simulations/cross_triple_orthogonality.py)).

**Theorem (F128, the factorization).** As an identity of rational functions of the six angles,

  **𝔉 = −(e₁ − f₁)² · 𝒪[ cos s · cot s · V_a V_b / P ].**

**Corollary (the sharper locus).** 𝔉 ≡ 0 on the hypersurface {e₁ = f₁}, wherever 𝔉 is defined (𝔉 carries its own cotangent poles; the identity is at the rational-function level throughout). One constraint; equality suffices, vanishing is not needed. F127's variety V = {e₁ = 0, f₁ = 0} is the codimension-two special case, so the corollary contains F127. The two inputs below are unconditional identities proved without any of §§4-6 of the F127 chain and without the wall, so this is a third proof of F127: route-distinct, though not trust-disjoint, from the residue chain (it consumes that chain's §3 and §3A gates; trust-disjointness remains the wall's role). It is the first *structural* proof with no analytic step: the residue-collapse chain carries §4's product rule on an algebraic curve, and the grid+CRT wall, itself free of analysis, is a certificate rather than a closed form.

Scope, honestly: until today every sharper {e₁ = f₁} statement in the repo was about the core function T and needed the sheet constraint ([PROOF_F127_RESIDUE_COLLAPSE §3](PROOF_F127_RESIDUE_COLLAPSE.md), "a bonus fact about T that does not widen V"). The theorem above is about 𝔉 itself and needs nothing else. The §3A scope note ("on V each T(L∘x) is generically nonzero, so the vanishing of the flip-sum is genuine cancellation", with L ∈ {±1}⁶ ranging over the sheet sign patterns and (L∘x)_k = L_k x_k) stays true; the cancellation now has its own two-line mechanism.

## §2 The two committed inputs

**(I) The global representation** ([PROOF_F127_RESIDUE_COLLAPSE §3A](PROOF_F127_RESIDUE_COLLAPSE.md), gate [f127_global_representation.py](../../simulations/f127_global_representation.py), term-by-term over ℚ(i)):

  𝔉 = −4 · 𝒪[ T · cot s ],

with T the F127 core function, T(a; b) = Σ_{i,j} (−1)^{i+j} α_i(a) α_j(b) · cot((a_i+b_j)/2), where α_i is the sine difference of the other two angles of its own triple (α_i(a) = sin a_l − sin a_j over the complement {j < l} of i).

**(II) The closed form of T** ([PROOF_F127_RESIDUE_COLLAPSE §3](PROOF_F127_RESIDUE_COLLAPSE.md), gate G1 of [f127_closed_form.py](../../simulations/f127_closed_form.py), 864 monomials per side over ℚ(i)):

  T · P = ⅛ · [ 2 cos s · ((e₁−f₁)² − 2 sin²s) + sin s · Σ_u sin 2x_u ] · V_a V_b.

Dividing by P and multiplying by cot s, the bracket splits into the (e₁−f₁)²-part and the rest (term bookkeeping only, using 4 cos s sin²s = 2 sin 2s · sin s and sin s · cot s = cos s):

  T · cot s = h + ¼ (e₁−f₁)² · cos s · cot s · V_a V_b / P, where
  **h = ⅛ · cos s · B · V_a V_b / P, B = Σ_u sin 2x_u − 2 sin 2s.**

Both e₁ and f₁ are invariant under every sign flip ε (cosine is even), so 𝒪 passes over (e₁−f₁)²:

  𝒪[T · cot s] = 𝒪[h] + ¼ (e₁−f₁)² · 𝒪[cos s · cot s · V_a V_b / P].

The theorem therefore reduces to one new fact.

## §3 The flip lemma

**Lemma.** 𝒪[h] ≡ 0 identically.

Multiply h by the flip-invariant combination P·P̃, where P̃ = Π_{i,j} sin((a_i − b_j)/2): pairwise, sin((a_i+b_j)/2) · sin((a_i−b_j)/2) = (cos b_j − cos a_i)/2, so P·P̃ = 2⁻⁹ Π_{i,j}(cos b_j − cos a_i) is even in every angle and passes through 𝒪. The lemma is then equivalent to the vanishing of the totally odd part of a trigonometric *polynomial*:

  **𝒪[ cos s · B · V_a V_b · P̃ ] ≡ 0.**

This has two independent exact proofs, both machine-run in the gate.

**Proof (A), brute force over ℤ.** In half-angle monomial coordinates t_u = e^{i x_u/2}, the polynomial 2(2i)¹⁶ · cos s · B · V_a V_b P̃ has integer coefficients and 8640 monomials; applying the signed character sum Π_u (1 − flip_u) returns the zero polynomial, exactly (gate G1). ∎

**Proof (B), the trade already knew.** The product V_a V_b P̃ is the full six-angle half-angle Vandermonde Δ = Π_{u<v} sin((x_u − x_v)/2) in the ordering (a₁, a₂, a₃, b₁, b₂, b₃). Three classical facts finish it:

1. **Weyl denominator folding** (gate G2, exact over ℤ). With z_u = e^{i x_u} and the alternant a_M = det[z_u^{m_k}],

   2(2i)¹⁵ · cos s · Δ = a_{M₁} + a_{M₂}, M₁ = (3, 2, 1, 0, −1, −2), M₂ = −M₁ reversed = (2, 1, 0, −1, −2, −3).

   (The product Π(z_u − z_v) is the Vandermonde alternant a_{(5,4,3,2,1,0)}; dividing by Πz_u^{5/2} and folding in the two exponentials of cos s shifts the exponent interval down by 2 and by 3.)

2. **Power-sum raising shifts** (gate G3, exact over ℤ; the Murnaghan-Nakayama bookkeeping at the alternant level). In z-units, 2i·B = p₂(z) − p₂(1/z) − 2Z + 2Z⁻¹ with p₂(z) = Σ z_u² and Z = z₁···z₆. Each piece acts on an alternant by exponent shifts: p₂(z)·a_M = Σ_k a_{M+2e_k}, p₂(1/z)·a_M = Σ_k a_{M−2e_k}, Z^{±1}·a_M = a_{M±1 on all}. So 2iB · (a_{M₁} + a_{M₂}) is a signed sum of 28 shifted alternants.

3. **The odd Weyl numerator dies** (gate G4, structural + exact). For any exponent set M, 64·𝒪[a_M] = det[z_u^{m_k} − z_u^{−m_k}]. This vanishes whenever M contains a zero (a zero column), a repeat (a_M = 0 outright), or a pair ±m (two columns equal up to sign). M₁ and M₂ are integer intervals containing 0; a ±2 shift of one exponent either lands on another exponent (repeat), or leaves the 0 in place; the ±1 shift of all exponents moves the interval so that it still contains 0. All 28 shifted sets die. ∎

Proof (B) is why the lemma is short rather than merely true: the flip-sum is a signed character sum over the sign subgroup (ℤ/2)⁶ of the hyperoctahedral group, the object under it is a Weyl denominator dressed by a power sum, and the representation theory of the trade has cancelled exactly this shape since the nineteenth century. The seventh visit's door-note applies verbatim: when a find is small, ask who has seen it before.

## §4 Assembly and the corollary

By the lemma, 𝒪[T·cot s] = ¼ (e₁−f₁)² · 𝒪[cos s · cot s · V_a V_b / P], and input (I) gives the theorem:

  𝔉 = −4·𝒪[T·cot s] = −(e₁−f₁)² · 𝒪[cos s · cot s · V_a V_b / P].

(Gate G5 pins the assembled identity numerically at generic points, rel. 10⁻¹³; gate G6 pins the corollary's content from below: generic |𝔉| has median ~10, individual samples spanning many orders since 𝔉 has other zero sets, while |𝔉| on {e₁ = f₁ = c}, c up to 1.9, is at machine zero.)

**Restriction hygiene for the corollary.** The right side is a rational function whose polar locus is supported on sines of the 32 sheet forms L·x and the 18 cross pair forms a_i ± b_j (V_a and V_b sit in the numerator; only P and cot s contribute poles, over the 64 flips). Each such factor vanishes on cosets of codimension-one subtori. e₁ − f₁ is not identically zero on any of these cosets: gate G7 conservatively exhibits a witness point on each of the 62 cosets of the pair and sheet families (a strict superset of the actual polar support) with |e₁−f₁| > 1.4; the evaluations are floating point with that wide margin, certifying the exact fact that the restriction is nonconstant. So the zero set of e₁ − f₁ meets each coset in dimension ≤ 4, and no component of the dimension-5 hypersurface {e₁ = f₁} lies inside the polar locus. On each component the identity therefore restricts: (e₁−f₁)² vanishes there, the cofactor is defined on a dense open subset, hence 𝔉 ≡ 0 on every component of {e₁ = f₁}. F127 is the sub-case e₁ = f₁ = 0. ∎

## §5 What changes and what does not

- **F127's statement, variety, tier, and wall are untouched.** The twinning arc consumes 𝔉 = 0 at points where both sums vanish separately; that remains the load-bearing case. The wall and the residue-collapse chain remain the two committed independent proofs; this chain is the third, and the first structural one with no analytic link (§4 of the F127 chain used a product rule on an algebraic curve; the wall is likewise analysis-free, but it is a grid certificate, not a closed form).
- **The T-scoped sharper-locus notes stay true** but are superseded in strength: the sharpening now reaches 𝔉 itself, without the sheet. The F127 surfaces carry dated pointers here rather than rewrites.
- **The cofactor's closed form (resolved 2026-07-17).** The cofactor W = 𝒪[cos s · cot s · V_a V_b / P] is the whole content of 𝔉 up to the explicit prefactor −(e₁−f₁)². Its closed form is F133, [PROOF_F133_W_SYMPLECTIC_CLOSED_FORM](PROOF_F133_W_SYMPLECTIC_CLOSED_FORM.md): a 143-term Sp(12) Weyl-character sum over the C₆ Weyl denominator, every coefficient an integer ≤ 8 in magnitude, derived exactly over ℤ.

## §6 Verification index

One gate script proves everything above; runtimes on the dev machine.

| gate | proves | method |
|---|---|---|
| G1 | the flip lemma, proof (A) | 8640-monomial integer polynomial annihilated by the signed character sum, exact over ℤ |
| G2 | Weyl denominator folding | dict identity 2(2i)¹⁵ cos s Δ = a_{M₁} + a_{M₂}, exact over ℤ |
| G3 | Murnaghan-Nakayama bookkeeping | 2iB·(a_{M₁}+a_{M₂}) = the 28 signed shifted alternants, exact over ℤ |
| G4 | the shifted alternants die | 0/repeat/±pair logic + exact 𝒪-projection of each |
| G5 | the assembled factorization | numeric pin vs the committed `cross_form`, 10 generic points |
| G6 | the sharper locus, from below | generic scale vs machine zero on {e₁ = f₁ = c}; flip-sum re-pin on V |
| G7 | restriction hygiene | witness points on all 62 pair/sheet cosets, a conservative superset of the polar support |

Run: `python simulations/f128_flip_sum_factorization.py` (~0.3 s, exit 0 iff all pass): [f128_flip_sum_factorization.py](../../simulations/f128_flip_sum_factorization.py).
