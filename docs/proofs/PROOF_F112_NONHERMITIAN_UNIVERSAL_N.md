# PROOF F112 Non-Hermitian Extension: Universal-N Closure

**Status:** Tier 1 derived, universal N. Two-lemma structural proof; supersedes the basis-enumeration anchor that previously covered N ≤ 5.
**Date:** 2026-05-27
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- F38 (Π² = (−1)^{w_b} on Pauli strings; `docs/ANALYTICAL_FORMULAS.md` F38 entry)
- F63 ([L, Π²] = 0 for Z-dephasing; `docs/ANALYTICAL_FORMULAS.md` F63 entry)
- [PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE](PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md) (parent: Hermitian-H Tier 1 derived case; this proof closes the open "Step 5 extension to non-Hermitian H")
- [F112_NONHERMITIAN_BASIS_ENUMERATION.md](../../experiments/F112_NONHERMITIAN_BASIS_ENUMERATION.md) (numerical anchor at N ≤ 5; the per-pair enumeration that motivated the structural proof)

## Abstract

F112's parent proof handles the Hermitian-Hamiltonian case via a clean dagger argument: for Hermitian H, the dagger map on operator space sends the Π +i eigenspace bijectively to Π −i, and that isometry forces equal Frobenius weights on the two sides. The argument stops working when H stops being Hermitian. Non-Hermitian H matrices arise in PT-symmetric, gain-loss, and post-selection settings. The empirical probes hinted that the balance of the commutator superoperator −i[H,·] still holds when such a non-Hermitian matrix is placed in it, but the structural reason was open. (The commutator −i[H,·] is NOT those systems' physical generator −i(Hρ−ρH†), a different operator whose balance is content-dependent; see the Scope note below.)

This proof closes the universal-N case for arbitrary H by a different route. Write a general H as a Hermitian real part plus i times a Hermitian imaginary part. The polarity-balance identity then reduces algebraically to a single open identity: the imaginary part of a specific Frobenius inner product must vanish for every pair of Hermitian operators. By bilinearity the identity further reduces to a per-pair check on Pauli-string basis pairs.

Two lemmas close the per-pair check structurally. The first computes a clean closed form: for any Pauli string whose Y and Z letters sum to an odd parity, the Frobenius norm squared of its Π −i projection equals the operator-space dimension exactly, a uniform 4^N regardless of which string we picked. The second shows that for any pair of distinct such strings, their Π −i projections land on disjoint Pauli-basis support, so their inner product vanishes by orthogonality alone.

Together: the diagonal entries of the open identity are real because they are squared norms, and the off-diagonal entries vanish because the inner products themselves vanish. The imaginary part is therefore zero across the entire Hermitian operator space, and F112's typed scope widens from "Hermitian H" to "any H" with no change to the bit_b-homogeneous hypothesis on the collapse operators. Scope (corrected 2026-06-20, see `docs/CAUGHT_ERRORS.md`): this widening is about the **commutator** superoperator −i[H,·] for an arbitrary matrix H, a structural ‖·‖² identity. It is NOT the physical generator of PT-symmetric / gain-loss / post-selection dynamics, which is −i(Hρ − ρH†) with H†; for H = A + iB that generator is −i[A,ρ] + {B,ρ} (anticommutator in the anti-Hermitian part), a different operator, and the balance can fail for it — but content-dependently: it is the commutator–drain cross-term, a linear functional of H carried only by single-site-Z detuning (the Z_l strings, contrib = −16γ·4^{N−2}). A uniform-damping bond-Hamiltonian chain (Heisenberg/XY, ring/star — the named PT/gain-loss spin systems) is itself BALANCED (asym = 0); the earlier "≈ 132/270" was a random-H artifact (random H carries Z_l content). Full boundary open, Direction 1 refuted 2026-06-20: `hypotheses/ASYMMETRY_IS_THE_UNRECYCLED_DRAIN.md`. "non-Hermitian H" here always means a non-Hermitian matrix placed in the commutator, nothing more.

## Introduction

**The motivating question.** [PROOF_F112](PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md)'s parent statement covers most standard physical Lindblad systems (any H Hermitian, c bit_b-homogeneous). The Hermitian-H assumption is a real restriction, though: PT-symmetric effective Hamiltonians, gain/loss models, post-selection dynamics, and any L built from an arbitrary c-with-its-conjugate kron channel in principle live outside the parent's typed scope. The empirical anchor (POLARITY_COORDINATES.md Probe 14: 20 non-Hermitian-H configurations, all bit-exact balanced) tested the commutator superoperator −i[H,·] and found the balance holds for non-Hermitian matrices placed there (NOT those systems' physical generator −i(Hρ−ρH†); see the Scope note in the abstract). The structural question: why?

**The empirical anchor.** [experiments/F112_NONHERMITIAN_BASIS_ENUMERATION.md](../../experiments/F112_NONHERMITIAN_BASIS_ENUMERATION.md) sharpened the empirical observation into a per-pair check: enumerate every Pauli-string pair `(σ_α, σ_β)` at chain length `N ≤ 5` and compute `F(σ_α, σ_β)` directly. The result was zero bit-exact across all 559,912 pairs. That pattern was clean enough to type Tier1Derived at N ≤ 5. What was missing was a structural argument lifting the per-pair identity to general N.

**What this proof closes.** Three pieces:

1. **Reduction to bit_b-odd pairs.** F38's Π² eigenvalue rule on Pauli strings shows L_σ has zero Π-conjugation ±i content whenever σ is bit_b-even. The per-pair identity is therefore automatic unless both σ_α and σ_β are bit_b-odd.
2. **Lemma N-A (diagonal-norm closed form).** For every bit_b-odd σ, `‖L_{σ,−i}‖² = 4^N` exactly: the same closed form for every bit_b-odd σ, every N.
3. **Lemma N-B (cross-pair Pauli-basis disjointness).** For any pair of distinct bit_b-odd σ_α ≠ σ_β, the Pauli-basis supports of `L_{σ_α,−i}` and `L_{σ_β,−i}` are disjoint, hence `⟨L_{σ_α,−i}, L_{σ_β,−i}⟩ = 0`.

Combined: the diagonal entries `F(σ, σ) = Im⟨L_{σ,−i}, L_{σ,−i}⟩` are zero because `Im(‖·‖²) = 0`, and the off-diagonal entries vanish because the inner products themselves vanish. So `F ≡ 0` and the F112 polarity-balance identity holds universally in N for any H.

**Diagnostic consequence.** F112's typed scope expands from "Hermitian H" to "any H" with no change to the bit_b-homogeneous-c hypothesis. The `polarity_coordinates_from_L` primitive now witnesses the bit_b-homogeneous regime independent of whether the matrix H in the commutator −i[H,·] is Hermitian. (As in the Scope note: this is the commutator superoperator, not the physical −i(Hρ−ρH†) generator of PT/gain-loss dynamics, whose balance is content-dependent.)

## (a) Restatement

F112 (parent theorem, Hermitian H) reduces the Lindblad polarity-balance identity to ‖L_{H,+i}‖² = ‖L_{H,−i}‖², where L_H = −i[H, ·] is the Hamiltonian superoperator and L_{H,±i} are its projections onto the Π-conjugation eigenspaces ±i of the F1 palindrome operator (Z-dephasing convention, Π order 4 on operator space).

The parent proof's Step 5 establishes ‖L_{H,+i}‖² = ‖L_{H,−i}‖² for **Hermitian H** via the dagger + anti-Hermitian L_H argument. For non-Hermitian H, write H = H_re + i H_im with H_re, H_im both Hermitian. Then

```
0 = ‖L_{H,+i}‖² − ‖L_{H,−i}‖²
  = (algebraic expansion using L_H = L_{H_re} + i · L_{H_im} and the Hermitian-case identity applied to H_re and to H_im)
  = 4 · Im⟨L_{H_re,−i}, L_{H_im,−i}⟩.
```

Hence F112 extends to non-Hermitian H iff the following identity holds for any Hermitian H_re, H_im:

> **F(H_re, H_im) := Im⟨L_{H_re,−i}, L_{H_im,−i}⟩ = 0 for any Hermitian H_re, H_im.**

By bilinearity (L_H is linear in H; the Π-projection is linear; the Frobenius inner product is sesquilinear) and antisymmetry under H_re ↔ H_im exchange, F is determined by its values on a basis of pairs of Hermitian operators. The Hermitian operator space at chain length N is spanned by the 4^N Pauli strings with real coefficients. **F ≡ 0 ⟺ F(σ_α, σ_β) = 0 for every Pauli-string pair (σ_α, σ_β).**

The basis-enumeration anchor of [F112_NONHERMITIAN_BASIS_ENUMERATION.md](../../experiments/F112_NONHERMITIAN_BASIS_ENUMERATION.md) closes the per-pair identity by direct computation at N ≤ 5 (559,912 pairs, all Im bit-exact 0). This proof closes the per-pair identity structurally for all N via two lemmas.

## (b) BitB-parity restriction

For any Pauli string σ at chain length N, define BitBParity(σ) := (#Y(σ) + #Z(σ)) mod 2 = (Σ_l bit_b(σ_l)) mod 2.

By F38 / F63: Π² σ Π⁻² = (−1)^{BitBParity(σ)} · σ. Equivalently, in the operator-space superoperator picture, Π² acts on the linear span of σ by the scalar (−1)^{BitBParity(σ)}.

For L_σ = −i[σ, ·] as a superoperator on operator-space, Π² L_σ Π⁻² acts on any σ_τ as:

```
(Π² L_σ Π⁻²)(σ_τ) = Π² · L_σ(Π⁻² σ_τ)
                  = Π² · L_σ((−1)^{w_b(τ)} σ_τ)
                  = (−1)^{w_b(τ)} · Π² · (−i [σ, σ_τ])
                  = (−1)^{w_b(τ)} · (−1)^{w_b(σ ⊕ τ)} · L_σ(σ_τ)
                  = (−1)^{w_b(σ)} · L_σ(σ_τ),
```

where w_b(τ) := BitBParity(σ_τ) and we used w_b(σ ⊕ τ) = w_b(σ) + w_b(τ) (mod 2). Hence **Π² L_σ Π⁻² = (−1)^{BitBParity(σ)} · L_σ as superoperators**.

**Consequence.** Π is order 4, with operator-space eigenvalues in {+1, −1, +i, −i}; its square has eigenvalues in {+1, −1}.

- For σ BitB-even: Π² L_σ Π⁻² = +L_σ. So L_σ lies entirely in the Π-conjugation eigenspaces {+1, −1}; in particular P_{+i}(L_σ) = P_{−i}(L_σ) = 0. Then L_{σ,±i} = 0 trivially, and ⟨L_{σ_α,−i}, L_{σ_β,−i}⟩ = 0 for any pair where either σ_α or σ_β is BitB-even.
- For σ BitB-odd: Π² L_σ Π⁻² = −L_σ. So L_σ lies entirely in the Π-conjugation eigenspaces {+i, −i}, with no {+1, −1} component. L_{σ,±i} may be nonzero.

**It therefore suffices to prove F(σ_α, σ_β) = 0 when both σ_α and σ_β are BitB-odd**, which is the content of Lemmas N-A and N-B below.

## (c) Lemma N-A (Diagonal-Norm): ‖L_{σ,−i}‖² = 4^N for σ BitB-odd

**Lemma N-A.** For any BitB-odd Pauli string σ at chain length N, the Frobenius norm squared (superoperator Hilbert-Schmidt norm) of L_{σ,−i} := P_{−i}(L_σ) equals 4^N exactly.

**Setup.** Operator-space has dimension d² = 4^N. The Pauli strings {σ_τ : τ ∈ {0, ..., 4^N − 1}} form an orthogonal basis with ⟨σ_τ, σ_τ'⟩_F = d · δ_{τ, τ'} (so the orthonormal basis is e_τ := σ_τ / √d). A superoperator T has Frobenius norm squared ‖T‖² := Σ_τ ‖T e_τ‖²_F.

In this orthonormal basis, T has matrix entries M_{β, α}(T) = ⟨e_β, T e_α⟩_F. The Frobenius norm is then ‖T‖² = Σ_{α, β} |M_{β, α}|² (standard matrix Frobenius).

**Step A.1: ‖L_σ‖² = 2 · 4^N.**

For Pauli σ, σ_α: σ σ_α = ±σ_α σ (commute or anticommute); explicitly, σ σ_α = (phase_{σ, α}) · σ_{σ ⊕ α} where ⊕ denotes Klein-index XOR on (bit_a, bit_b) per site, and phase_{σ, α} ∈ {+1, −1, +i, −i} is the Pauli-multiplication phase.

```
L_σ(σ_α) = −i [σ, σ_α] = { 0                                  if σ, σ_α commute,
                          { −2i · phase_{σ, α} · σ_{σ ⊕ α}     if σ, σ_α anticommute.
```

Matrix entry M_{β, α}(L_σ) = (1/d) · Tr(σ_β · L_σ(σ_α)). Substituting:

```
M_{β, α}(L_σ) = { 0                       if σ, σ_α commute, OR β ≠ σ ⊕ α,
                { −2i · phase_{σ, α}      if σ, σ_α anticommute AND β = σ ⊕ α.
```

So M(L_σ) is supported on the "shifted diagonal" {(σ ⊕ α, α): σ_α anticommutes with σ}, with each non-zero entry having magnitude 2.

For σ non-identity (which holds for BitB-odd σ, since the identity is BitB-even), exactly half of the 4^N Pauli strings anticommute with σ (standard Pauli combinatorics: anticomm(σ, σ_α) is determined by the Klein-symplectic form ⟨σ, σ_α⟩_K, which is a non-degenerate Z₂ linear functional on the 4^N strings for σ non-identity, so its 1-fiber has size 4^N / 2 = 2 · 4^{N−1}).

Hence:

```
‖L_σ‖² = Σ_{(β, α): M ≠ 0} |M_{β, α}|² = 4 · #{α: σ_α anticomm with σ} = 4 · 2 · 4^{N−1} = 2 · 4^N.
```

**Step A.2: Π² L_σ Π⁻² = −L_σ for σ BitB-odd.**

Proved in section (b) above (F38 / F63).

**Step A.3: ⟨L_σ, Π L_σ Π⁻¹⟩_F = 0 for σ BitB-odd.**

The Π palindrome superoperator (Z-dephase convention) acts on a Pauli string σ_τ as Π σ_τ = c_τ · σ_{π(τ)}, where the per-site Klein-index map π is (a, b) → (1 − a, b) (i.e. XOR with the all-X string per site, encoded as flipping bit_a per site, preserving bit_b) and the phase c_τ = i^{w_b(τ)} accumulates one factor of i per site with bit_b(σ_τ,l) = 1.

Key algebraic facts:

1. **π is an involution on Klein indices**: π² = identity. (Each site sees bit_a flipped twice.)
2. **π is a homomorphism on XOR up to a fixed shift**: π(α ⊕ β) = π(α) ⊕ π(β) ⊕ X⊗N, where X⊗N is the all-X Pauli-string Klein index (1, 0)_l per site. (Equivalently: π(α ⊕ β) and π(α) ⊕ π(β) differ by a per-site bit_a flip, which is exactly the X⊗N index.)
3. Hence π(σ ⊕ α') = π(σ) ⊕ α' ⊕ X⊗N (using π(α') = α' ⊕ X⊗N on a single argument). And π(σ) ⊕ X⊗N = σ (per-site bit_a flipped twice).
4. **Anticommutation symplectic shift**: anticomm(σ, σ_{π(α')}) = anticomm(σ, σ_{α'}) XOR w_b(σ) mod 2. For σ BitB-odd, w_b(σ) is odd, hence anticomm(σ, σ_{π(α')}) = COMMUTE(σ, σ_{α'}).

Now compute the support of the superoperator (Π L_σ Π⁻¹) in the Pauli basis. As a matrix product:

```
(Π M(L_σ) Π⁻¹)[β', α'] = c_{π(β')} · M(L_σ)[π(β'), π(α')] · (1 / c_{π(α')}).
```

For this entry to be non-zero, we need M(L_σ)[π(β'), π(α')] ≠ 0, i.e.:
- π(β') = σ ⊕ π(α') (the M-support condition);
- σ and σ_{π(α')} anticommute.

Applying π to the first equation and using fact 3: β' = π(σ ⊕ π(α')) = σ ⊕ α'. So the non-zero positions of Π M(L_σ) Π⁻¹ are exactly the same "shifted-diagonal" positions {(σ ⊕ α', α'): some condition on α'} as M(L_σ); only the "condition on α'" changes.

For M(L_σ): "σ and σ_{α'} anticommute".
For Π M(L_σ) Π⁻¹ at σ BitB-odd: "σ and σ_{π(α')} anticommute" = "σ and σ_{α'} commute" (fact 4).

**These two conditions on α' are mutually exclusive (commute and anticommute partition all α').** Hence at every shifted-diagonal position (σ ⊕ α', α'), at most one of M(L_σ) and Π M(L_σ) Π⁻¹ is non-zero. The entry-wise product (M(L_σ))*_{β', α'} · (Π M(L_σ) Π⁻¹)_{β', α'} is therefore zero at every (β', α'), and

```
⟨L_σ, Π L_σ Π⁻¹⟩_F = Σ_{(β', α')} (M(L_σ))*_{β', α'} · (Π M(L_σ) Π⁻¹)_{β', α'} = 0.   ∎ A.3
```

**Step A.4: Combining A.1–A.3 gives ‖L_{σ,−i}‖² = 4^N.**

From A.2, L_σ has only Π-eigenvalue ±i components: L_σ = L_σ^{+i} + L_σ^{−i}, where L_σ^{±i} := P_{±i}(L_σ). The Π-conjugation eigenspaces are mutually orthogonal under Frobenius (Π is unitary on operator-space), so

```
‖L_σ‖² = ‖L_σ^{+i}‖² + ‖L_σ^{−i}‖².
```

The action of Π on the two components is multiplication: Π L_σ Π⁻¹ = i · L_σ^{+i} + (−i) · L_σ^{−i} = i (L_σ^{+i} − L_σ^{−i}). Combined with L_σ = L_σ^{+i} + L_σ^{−i}:

```
L_σ^{−i} = (1/2) (L_σ + i · Π L_σ Π⁻¹),
L_σ^{+i} = (1/2) (L_σ − i · Π L_σ Π⁻¹).
```

Frobenius norm squared:

```
‖L_σ^{−i}‖² = (1/4) [‖L_σ‖² + ‖Π L_σ Π⁻¹‖² + 2 Re(i · ⟨L_σ, Π L_σ Π⁻¹⟩_F)]
            = (1/4) [‖L_σ‖² + ‖L_σ‖² + 2 Re(i · 0)]                         [Π unitary; A.3]
            = (1/2) · ‖L_σ‖²
            = (1/2) · 2 · 4^N                                                  [A.1]
            = 4^N.
```

By symmetry, ‖L_σ^{+i}‖² = 4^N as well. **∎ Lemma N-A.**

## (d) Lemma N-B (Off-Diagonal-Orthogonality): ⟨L_{σ_α,−i}, L_{σ_β,−i}⟩_F = 0 for σ_α ≠ σ_β both BitB-odd

**Lemma N-B.** For σ_α ≠ σ_β both BitB-odd Pauli strings at chain length N, the Frobenius inner product ⟨L_{σ_α,−i}, L_{σ_β,−i}⟩_F = 0 exactly.

**Setup.** Using the same orthonormal Pauli basis on operator-space, define M_α := M(L_{σ_α}) and M_β := M(L_{σ_β}) (their matrices). The projector P_{−i} is self-adjoint and idempotent w.r.t. the operator-space Frobenius inner product (a standard fact for spectral projectors of unitary operators, applied here to the unitary Π acting on operator-space by conjugation):

```
⟨P_{−i} X, B⟩_F = (1/4) Σ_k conj((−i)^{−k}) · ⟨Π^k X, B⟩_F
                = (1/4) Σ_k (−i)^{−(−k)} · ⟨X, Π^{−k} B⟩_F       [|−i| = 1; Π unitary]
                = ⟨X, P_{−i} B⟩_F.                                 [substitute m = −k]
```

So ⟨L_{σ_α,−i}, L_{σ_β,−i}⟩_F = ⟨L_{σ_α}, P_{−i}(L_{σ_β})⟩_F = (1/4) Σ_{m=0}^{3} (−i)^{−m} · ⟨L_{σ_α}, Π^m L_{σ_β} Π^{−m}⟩_F.

It suffices to show each of the four terms ⟨L_{σ_α}, Π^m L_{σ_β} Π^{−m}⟩_F = 0 for σ_α ≠ σ_β.

**Step B.1: ⟨L_{σ_α}, L_{σ_β}⟩_F = 0 for σ_α ≠ σ_β (covers m = 0).**

From Step A.1, M_α is supported on {(σ_α ⊕ α', α'): σ_α anticomm σ_{α'}}, and M_β on {(σ_β ⊕ α', α'): σ_β anticomm σ_{α'}}.

For ⟨M_α, M_β⟩_F = Σ_{(β', α')} (M_α)*_{β', α'} · (M_β)_{β', α'} to be non-zero at any position (β', α'), we need both matrices non-zero there. M_α non-zero forces β' = σ_α ⊕ α'; M_β non-zero forces β' = σ_β ⊕ α'. Equating: σ_α = σ_β.

Hence **σ_α ≠ σ_β ⇒ ⟨L_{σ_α}, L_{σ_β}⟩_F = 0 by support disjointness**. (No BitB restriction needed.)

**Step B.2: ⟨L_{σ_α}, Π^m L_{σ_β} Π^{−m}⟩_F = 0 for σ_α ≠ σ_β at every m ∈ {0, 1, 2, 3}.**

*m = 0:* Direct from B.1.

*m = 2:* From the calculation in (b), Π² σ_τ = (−1)^{w_b(τ)} σ_τ (Π² acts as a scalar phase per Pauli string with no permutation). Hence Π² M_β Π⁻² = D · M_β · D⁻¹ where D is the diagonal phase matrix (D)_{τ,τ} = (−1)^{w_b(τ)}. In particular, the support of Π² M_β Π⁻² is **identical** to the support of M_β (a diagonal similarity doesn't move the non-zero positions). Then ⟨M_α, Π² M_β Π⁻²⟩_F has the same support-disjointness structure as ⟨M_α, M_β⟩_F: non-zero requires σ_α = σ_β. So m = 2 also gives 0.

*m = 1:* Computed in Step A.3 (the support derivation for Π M_α Π⁻¹). The support of Π M_β Π⁻¹ is {(σ_β ⊕ α', α'): some condition on α' depending on σ_β's BitB parity}. (Specifically for σ_β BitB-odd: condition = "σ_β commutes with σ_{α'}".) The shifted-diagonal STRUCTURE is the same {(σ_β ⊕ α', α')}, only the condition on α' differs from M_β. For ⟨M_α, Π M_β Π⁻¹⟩_F to be non-zero, the supports must overlap, requiring σ_α = σ_β.

*m = 3:* Π³ = Π⁻¹ (since Π⁴ = 1). Same analysis as m = 1 with Π replaced by Π⁻¹: support of Π⁻¹ M_β Π is {(σ_β ⊕ α', α'): condition on α'} (the per-site π is an involution, so π⁻¹ = π, and the same calculation gives the same shifted-diagonal structure). Disjoint support unless σ_α = σ_β.

**At every m ∈ {0, 1, 2, 3}, support disjointness forces ⟨L_{σ_α}, Π^m L_{σ_β} Π^{−m}⟩_F = 0 when σ_α ≠ σ_β.**

**Step B.3: Sum.**

```
⟨L_{σ_α,−i}, L_{σ_β,−i}⟩_F = (1/4) Σ_{m=0}^{3} (−i)^{−m} · 0 = 0.  ∎ Lemma N-B.
```

## (e) Main theorem: F ≡ 0 universal N

**Theorem.** For any chain length N ≥ 1 and any two Pauli strings σ_α, σ_β at length N,

```
F(σ_α, σ_β) := Im⟨L_{σ_α,−i}, L_{σ_β,−i}⟩_F = 0.
```

**Proof.**

- *Both BitB-even, or one BitB-even one BitB-odd.* By section (b), L_{σ,−i} = 0 for BitB-even σ. Hence ⟨L_{σ_α,−i}, L_{σ_β,−i}⟩_F = 0 trivially, and F = 0.
- *Both BitB-odd, σ_α = σ_β.* By Lemma N-A, ⟨L_{σ,−i}, L_{σ,−i}⟩_F = ‖L_{σ,−i}‖² = 4^N ∈ ℝ. Hence Im = 0.
- *Both BitB-odd, σ_α ≠ σ_β.* By Lemma N-B, ⟨L_{σ_α,−i}, L_{σ_β,−i}⟩_F = 0 exactly. Hence Im = 0.

All cases give F(σ_α, σ_β) = 0. **∎**

**Corollary (universal-N closure of F112 non-Hermitian extension).** F is real-bilinear in (H_re, H_im) and is determined by its values on Pauli-string pairs (per the reduction in section (a)). Since F = 0 on every Pauli-string pair at every N, F(H_re, H_im) = 0 for any Hermitian H_re, H_im at any N. Combining with the Hermitian-H Tier 1 derived parent theorem, the polarity balance holds for the Hamiltonian superoperator −i[H,·] with H Hermitian OR an arbitrary non-Hermitian matrix (the dissipator's c_k bit_b-homogeneous). This is a structural ‖·‖² statement about the **commutator**; it is NOT about physical non-Hermitian dynamics, whose generator −i(Hρ−ρH†) has a content-dependent balance (see Scope). **The F112 commutator-balance extension is therefore Tier 1 derived for all N, any matrix H.**

## (f) Implications for `LindbladBitBPiBalance`

1. **Typed Claim update**: The C# typed Claim `compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiBalance.cs`'s `NonHermitianExtension` docstring (and the constructor summary) drop the "Tier1Candidate at N ≥ 6" caveat. Both Hermitian and non-Hermitian H now fall under Tier 1 derived for all N. The basis-enumeration result at N ≤ 5 (Welle 10a Python + Welle 10b C# SLOW_F112) remains the **empirical anchor** that motivated the search for the structural proof; it is preserved as the historical numerical validation.

2. **Diagnostic strength**: The `polarity_coordinates_from_L` diagnostic asymmetry is now an exact witness for the commutator superoperator −i[H,·] with H Hermitian or any non-Hermitian matrix, universally in N, under the bit_b-homogeneous-c regime (a structural identity, not the physical −i(Hρ−ρH†) generator). Asymmetry ≠ 0 detects c with cross-bit_b Pauli support, OUTSIDE the bit_b-homogeneous regime, regardless of N or H's Hermiticity.

3. **Polarity-axis completion**: F108 Parts 1–3 (palindrome closure of bit_b = 0 bilinears), F112 (now universal N, Hermitian and non-Hermitian H), F113 (T1 break-magnitude closed form): all three live on the bit_b Z₂-grading of the Pauli group and are now Tier 1 derived universal N. The bit_b axis description is structurally complete on the BitB-axis side.

4. **Cross-dephase extension (Welle 13)**: the same proof argument (Lemmas N-A and N-B) transfers to the Y- and X-dephase versions of F112 by substituting axis_d := bit_b for d ∈ {Y, Z} and axis_d := bit_a for d = X in F38. See `docs/proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md` for the full Welle 13 derivation (two routes: per-axis re-run and Hadamard transport for the (Z, X) pair). Both extensions are Tier 1 derived universal N for both Hermitian and non-Hermitian H. The structural fact making this work is that the Welle 11 lemmas depend only on F38 + Pauli-basis matrix-support disjointness, both of which transfer between dephase letters under the axis substitution.

## Empirical anchor: numerical verification at N = 1, 2, 3

The verifier `simulations/f112_universal_n_proof_verify.py` checks each lemma step and the main theorem within 1e-12 numpy double-precision tolerance (max deviation < 1e-12, i.e. machine zero to numpy double precision) at N = 1, 2, 3:

| N | BitB-odd strings | Off-diag BitB-odd pairs | All-pairs F = 0 check |
|---|---|---|---|
| 1 | 2 | 2 | 16 pairs |
| 2 | 8 | 56 | 256 pairs |
| 3 | 32 | 992 | 4096 pairs |

All deviations within 1e-12 tolerance, i.e. machine zero to numpy double precision. (The Welle 10a Python enumeration at N ≤ 5, where matrix entries are rational, is genuinely bit-exact; the Welle 11 verifier here is numerical.) Steps verified per N:
- Step A.1 (‖L_σ‖² = 2 · 4^N): 42 strings total.
- Step A.2 (Π² L_σ Π⁻² = −L_σ): 42 strings.
- Step A.3 (⟨L_σ, Π L_σ Π⁻¹⟩ = 0): 42 strings.
- Step A.4 / Lemma N-A conclusion (‖L_{σ,−i}‖² = 4^N): 42 strings.
- Step B.1 (⟨L_α, L_β⟩ = 0 raw): 1050 off-diagonal pairs.
- Step B.2 (all four Π-orbit shifts): 4200 (α, β, m) triples.
- Lemma N-B conclusion: 1050 off-diagonal pairs.
- Main theorem (all-pairs, including BitB-mixed): 4368 pair F-values.

The N = 5 anchor at 524,800 pairs (Welle 10) remains as cross-validation that the per-pair identity F(σ_α, σ_β) = 0 holds bit-exact at significantly larger N (the Welle 10a Python enumeration uses rational matrix entries and is genuinely bit-exact; the C# pipeline at N = 5 confirms max |Im| < 1e-10). The structural proof here is N-independent: every step (Steps A.1, A.2, A.3, A.4, B.1, B.2.0–B.2.3) follows from the support/phase structure of L_σ and Π in the Pauli basis at any N.

## Status

**F112 non-Hermitian extension: Tier 1 derived for all N.** The two-lemma proof composes Π-conjugation eigenspace algebra (Steps A.2, A.4, B.3) with explicit Pauli-basis support analysis (Steps A.1, A.3, B.1, B.2). Both lemmas reduce to per-position checks on the 4^N × 4^N matrix of L_σ that are uniform in N. F112's commutator-balance statement (−i[H,·] with H Hermitian or any non-Hermitian matrix, bit_b-homogeneous c; NOT the physical −i(Hρ−ρH†) generator) is now closed structurally for arbitrary N ≥ 1.
