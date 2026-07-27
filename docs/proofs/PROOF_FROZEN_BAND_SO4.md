# What the frozen band contains: two ladders, an invariant disagreement, and the seed rung in closed form

**Status:** Section 2 (both ladders, Σ-oddness as the exact condition for the SPIN ladder, and Lemma 2.5 deriving the off-diagonal floor from the corner's by injectivity) Tier 1 derived, with Corollary 2.4 recording that Σ-oddness is sufficient for that floor and NOT necessary for the off-diagonal band to exist; the matching off-diagonal EQUALITY stays measured, as it is the ceiling; Section 3 (the band as ⌊N/2⌋ copies of one SO(4) irrep) Tier 1 derived as a FLOOR, with the equality resting on the same measured depth Section 7 leaves open, and with j_spin ≤ 1 taken from the measured band statement of [XY_FROZEN_BAND](../../experiments/XY_FROZEN_BAND.md); Section 4 (F141, the invariant disagreement and the inheritance of a rung's spectrum from the rungs below) Tier 1 derived, with the criterion it yields explicitly one-directional and living on ker(ad_h) rather than on the full lowest-weight spaces; Section 5 (the large-J reduction, and every even order vanishing by transpose parity) Tier 1 derived; Section 6 (F143, the seed rung's spectrum in closed form) Tier 1 derived. The upper bound at rungs ℓ ≥ 2, which is what an all-N ceiling needs, is NOT proved here and Section 7 says exactly what is missing.
**Date:** 2026-07-27
**Authors:** Thomas Wicht, Claude (Anthropic)
**Statement:** On the open XY chain under Z-dephasing, the band of blocks carrying the frozen root λ = −4γ̄ **contains** a single representation-theoretic object: ⌊N/2⌋ copies of the irreducible representation of Spin(4) = SU(2) × SU(2) with (η-spin, spin) = (N/2 − 1, 1), which descends to SO(4) only at even N, which places ⌊N/2⌋ in every one of the 3(N − 1) band blocks and is the proved floor. That the band is nothing more than those copies is the per-N measured ceiling, open for all N and stated as such in Section 7. The two ladders that build it are the two commuting SU(2)s of the Hubbard model in disguise; the η one needs nothing of the Hamiltonian beyond number-conserving quadratic, and the spin one commutes exactly when the single-excitation matrix h is odd under the staggering Σ = diag((−1)^l). That condition is sufficient for the off-diagonal floor through the spin ladder and is not necessary for the off-diagonal band to be occupied, so it names a mechanism and does not replace the measured band-existence gate. The disagreement count is η-invariant, so it acts only on the multiplet-counting spaces, and the spectrum on rung p is the union of the lowest-weight spectra at rungs 0 ≤ ℓ ≤ min(p, N−p). At the seed rung the resulting operator is G = (1/M)(J + (I+R)/2) with M = N + 1, whose spectrum is {0 with multiplicity ⌊N/2⌋, 1/M with multiplicity ⌈N/2⌉ − 1, 1 simple}: the frozen kernel is exactly the chiral-odd sector, and the gap above it is 1/(N+1).
**Verification:** [`simulations/eta_ceiling_reduction.py`](../../simulations/eta_ceiling_reduction.py) (must print "eta-ceiling reduction gate: ALL GREEN", 170 checks, about a minute; `--deep` adds the larger rungs, runs a few minutes, and carries F143 to N = 40)
**Depends on:** [PROOF_R90_FROZEN_DIVISOR](PROOF_R90_FROZEN_DIVISOR.md) (the frozen root, its corner, and the proved floor), [XY_FROZEN_BAND](../../experiments/XY_FROZEN_BAND.md) (the band |p − q| ∈ {0, 2}, the η ladder, and the per-N ceiling this document tries to make N-free)

---

## What this means

The frozen band was found one block at a time: a diagonal line of blocks all carrying the same undisturbable decay rate, plus two lines beside it, and a ladder that climbs the diagonal one excitation at a time. This document says the three lines are not three facts. They are one object seen from three angles, and the object is an old one: the chain under dephasing is the Hubbard model wearing a disguise, and a Hubbard model has two independent SU(2) symmetries rather than one. The ladder already known is the first of them. The second is what steps sideways, off the diagonal onto the two neighbouring lines, and that is the whole reason the band is three lines wide and not one.

The second symmetry comes with a condition, and the condition is about the matrix rather than about its spectrum: the one-excitation matrix has to connect only sites of opposite parity. That is strictly stronger than having a spectrum symmetric about zero, which it implies and which does not imply it. What it buys is a mechanism: where it holds, the sideways ladder carries the corner's frozen modes onto the two neighbouring lines, and the floor there is the corner's floor for a reason rather than by measurement.

It does not, however, decide whether those lines are occupied. A star, one site joined to every other, fails the parity condition and fills them anyway. So the measured test the band note uses to predict occupancy is not superseded here, and the parity condition is one sufficient route to the floor rather than the gate itself. Getting this the wrong way round is easy, and this document did, until the star was asked.

There is a second thing the disguise gives, and it is the one that makes the remaining question small. Ask how much the two sides of a coherence disagree, count the sites where one has an excitation and the other does not, and that count does not notice the ladder at all: adding a matched pair to both sides adds no disagreement. So the whole dissipative half of the problem lives on the *count of multiplets* rather than on the multiplets themselves, and the question of how deep the freezing goes on a high rung reduces to the same question on the rungs below it. Everything the band does is then decided by what starts, and where.

At the bottom rung that decision can be written out completely. The operator that survives the reduction there is a two-line matrix built from the all-ones matrix and the mirror that swaps mode k with mode N+1−k, and reading its spectrum takes no more than noticing that the two commute. The frozen modes are exactly the combinations odd under that mirror, one per mirror pair, which is the ⌊N/2⌋ the whole arc has been counting, now in the mode picture rather than the site picture. And the next eigenvalue above them is 1/(N+1), which is not a chain constant but the length of the discrete sine transform the chain diagonalises under.

## 1. Setting

Open chain of N sites, H = (J/2)·Σ_b (X_b X_{b+1} + Y_b Y_{b+1}), so the single-excitation matrix h is the N×N hopping matrix with h_{ab} = J on neighbours. Site-resolved Z-dephasing with rates γ_l, mean γ̄, and

  L(ρ) = −i[H, ρ] + Σ_l γ_l (Z_l ρ Z_l − ρ).

Cells are |A⟩⟨B| with A and B site SETS; L preserves the joint-popcount block (p, q) = (|A|, |B|). The dephasing part is diagonal on cells with entry

  −2·Σ_{l ∈ AΔB} γ_l,

so the **disagreement set** AΔB is the only thing it sees. The frozen root is λ = −4γ̄, and the band is the set of blocks carrying it, |p − q| ∈ {0, 2} less the two blocks (0,0) and (N,N), as established in [XY_FROZEN_BAND](../../experiments/XY_FROZEN_BAND.md).

**The doubled picture.** Read the ket index as one fermion species and the bra index as a second: a cell |A⟩⟨B| becomes a state with ↑-occupation A and ↓-occupation B. Right multiplication transposes, so H acting from the right becomes h^T on the ↓ species, and

  −i[H, ·]  becomes  −i·K,  K = Σ_{ab} h_{ab} c†_{a↑}c_{b↑} − Σ_{ab} h_{ba} c†_{a↓}c_{b↓}.

Write **𝒦** for the disagreement count as an operator: 𝒦 = Σ_l n_{l↑}(1 − n_{l↓}), which on a cell reads |A \ B|, and write 𝒦' for its mirror |B \ A|. At uniform rates the dissipator is −2γ(𝒦 + 𝒦'). **On a diagonal block, and only there**, the two agree, |A \ B| = |B \ A| = |AΔB|/2, and the dissipator collapses to −4γ·𝒦. On an off-diagonal band block it does not: a frozen cell there has |A \ B| = 2 and |B \ A| = 0, so the dissipator is −4γ while −4γ·𝒦 would be −8γ. Everything below that evaluates the spectrum of 𝒦 does so on a diagonal rung, where the collapse is valid.

## 2. The two ladders

Two maps on cells, both raising:

  **Φ**(ρ) := Σ_l d†_l ρ d_l = Σ_l c†_{l↑}c†_{l↓},   (p, q) ↦ (p+1, q+1)
  **S⁺**(ρ) := Σ_l (−1)^l d†_l ρ d†_l = Σ_l (−1)^l c†_{l↑}c_{l↓},   (p, q) ↦ (p+1, q−1)

Φ is Yang's η-pairing without the staggering, S⁺ is the spin raising with it. The swap of which one carries the stagger is the partial particle-hole transformation that turns this Liouvillian into a Hubbard Hamiltonian: the ↓ species hops with the opposite sign here, and that sign is exactly what moves the stagger from one SU(2) to the other.

**Lemma 2.1 (both are blind to the watching, at any profile).** Φ and S⁺ commute with the dissipator for every rate profile γ_l, and not only for the uniform one.

*Proof.* The dissipator is diagonal on cells with entry −2Σ_{l ∈ AΔB} γ_l, so it suffices that both maps leave AΔB unchanged. Φ acts on a cell only at a site l with l ∉ A and l ∉ B, so l ∉ AΔB before, and afterwards l lies in both A and B, so l ∉ AΔB again; no other site moves. S⁺ acts only at a site l with l ∉ A and l ∈ B, so l ∈ B \ A ⊆ AΔB before, and afterwards l ∈ A \ B ⊆ AΔB; again no other site moves. ∎

**Lemma 2.2 (Φ and the turning).** [K, Φ] = 0 for every single-excitation matrix h, symmetric or not.

*Proof.* [Σ_{ab} h_{ab} c†_{a↑}c_{b↑}, Σ_l c†_{l↑}c†_{l↓}] = Σ_{a,l} h_{al} c†_{a↑}c†_{l↓}, and the ↓ term contributes −Σ_{a,l} h_{la} c†_{l↑}c†_{a↓}, which after renaming (a ↔ l) is −Σ_{a,l} h_{al} c†_{a↑}c†_{l↓}. The two cancel. ∎

This is the ladder already carried by [XY_FROZEN_BAND](../../experiments/XY_FROZEN_BAND.md); it is restated here because Lemma 2.3 is its sibling and reads best beside it. Note what the two lemmas together say: Φ commutes with the two halves of L **separately**, and for independent reasons.

**Lemma 2.3 (S⁺ and the turning).** Let Σ := diag((−1)^l). For real symmetric h,

  [K, S⁺] = Σ_{a,l} h_{al}·((−1)^a + (−1)^l)·c†_{a↑}c_{l↓},

which vanishes identically **iff h_{al} = 0 whenever a and l have the same parity**, that is, iff **Σ h Σ = −h**.

*Proof.* The ↑ term of K contributes Σ_{a,l} h_{al}(−1)^l c†_{a↑}c_{l↓}; the ↓ term contributes +Σ_{b,l} h_{bl}(−1)^l c†_{l↑}c_{b↓}, which after renaming is Σ_{a,l} h_{la}(−1)^a c†_{a↑}c_{l↓}. Adding and using h_{la} = h_{al} gives the displayed form. Each coefficient vanishes separately, and (−1)^a + (−1)^l is 0 for opposite parities and ±2 for equal ones. ∎

**Corollary 2.4 (what Σ-oddness does and does not decide).** Σ h Σ = −h implies that the spectrum of h is symmetric about zero, because Σ is an involution conjugating h to −h, and the converse fails: the **star**, one site joined to every other, is bipartite and so has a symmetric spectrum, but its two colours are centre-against-leaves rather than even-against-odd, so Σ h Σ ≠ −h and S⁺ does not commute there.

So Σ-oddness is **sufficient** for the off-diagonal floor through the spin ladder: where it holds, S⁺ carries the corner's frozen space onto the block (2,0) and the floor is inherited. It is **not necessary** for that block to be occupied. The star fills it regardless, carrying 2, 5 and 9 frozen modes at N = 5, 6 and 7 while failing Σ-oddness, and its (1,1) block carries far more than ⌊N/2⌋ besides, because its hopping matrix is highly degenerate.

Those three counts are for the **spin** sector, which is what L is built from, and a reader who reaches for free fermions will get different ones. The star's bonds are not adjacent in any site ordering, so the Jordan-Wigner strings do not cancel and the star is not a free-fermion model; pricing the same bonds by order instead of by occupation gives 4, 7 and 11 there, the Slater count 1 + C(N−2, 2). That is the arc's own pricing distinction turning up again, and it is worth stating because the free-fermion number is the one a quick recomputation produces. The measured band-existence test of [XY_FROZEN_BAND](../../experiments/XY_FROZEN_BAND.md) is therefore not superseded by anything here, and Σ-oddness is not "the seed gate"; it is the condition for one mechanism.

This paragraph replaces an earlier one that had it backwards, and the way it went wrong is worth keeping. Gate check V1(b) verifies that Σ-oddness and S⁺-commutation agree, which is Lemma 2.3, and the star was added to it as the case separating Σ-oddness from the symmetric spectrum. That much is real. But S⁺-commutation is not band occupancy, and the check was then read as evidence for a claim about occupancy that it never touched. **A verification can be green, and be about the wrong thing.** V1(c) now asks the star the occupancy question directly, and the answer is the one above.

**Lemma 2.5 (S⁺ does not annihilate a seed, so the off-diagonal floor is derived and not measured).** Assume Σ h Σ = −h. Then S⁺ is injective on the frozen space of the corner block, and consequently each of (2,0) and (0,2) carries at least ⌊N/2⌋ frozen modes, at every N.

*Proof.* Three steps, each exact.

1. **S⁺ is the commutator with Σ.** On the corner block write X for the N×N matrix of the cell coefficients. S⁺ lands in the block (2,0), which is Λ²(C^N), and for **symmetric** X its image has Λ²-coordinates (S⁺X)_{ab} = [Σ, X]_{ab}. (Gate V2(b) pins this identity at residual exactly zero; it is the statement that acting with two creation operators antisymmetrises, and the staggering supplies the Σ.)
2. **The seeds are Σ-odd.** Σ conjugates h to −h, so it carries the ε_k eigenvector to an eigenvector of −ε_k, giving Σ P_k Σ = P_{k̄} for the chiral partner k̄. Hence Σ (P_k − P_{k̄}) Σ = −(P_k − P_{k̄}), and the same holds for every real combination of seeds, since the relation is linear.

   Worth spelling out, because it is the step a reader is most likely to doubt: these seeds are frozen at **every** coupling, not only in the large-J limit, and they span the whole corner frozen space there. Two facts do it, both from Section 6. Each seed commutes with h, so the Hamiltonian part annihilates it; and each has zero diagonal, because chiral partners have identical squared amplitudes, so the rate part multiplies it by the single constant −4γ̄ that the two-disagreement cells carry. Neither fact mentions J. The corner frozen space is accordingly a space of SYMMETRIC matrices at every coupling, which is what steps 1 and 3 need, and it is not the generic non-symmetric space one might expect from a non-Hermitian generator.
3. **An odd symmetric matrix has a nonzero commutator.** For Σ X Σ = −X and Σ² = I one has X Σ = −Σ X, so [Σ, X] = Σ X − X Σ = 2 Σ X, which vanishes only if X = 0 because Σ is invertible.

So S⁺X ≠ 0 for every nonzero X in the seed space, the image has dimension ⌊N/2⌋, and it consists of frozen vectors because S⁺ commutes with L. ∎

This is what upgrades the off-diagonal floor. Commutation alone gives only that S⁺ maps the corner's frozen space **into** (2,0); a seed annihilated by S⁺ would be a spin singlet and would contribute nothing there. Lemma 2.5 rules that out for all N, where gate check V2 only ever measured it at the N its census reaches. What stays measured is the matching **equality**, that (2,0) carries no more than ⌊N/2⌋, which is the ceiling and is not this section's.

## 3. The band as one SO(4) multiplet

**Lemma 3.1 (the seeds are η-lowest-weight).** Every frozen vector of the corner block (1,1) is traceless, hence annihilated by Ψ = Φ†.

*Proof.* Ψ maps (1,1) to the one-dimensional block (0,0) and is there proportional to the trace. The corner's frozen space is spanned by the differences P_k − P_{k̄} of spectral projectors of h at chiral partner pairs (Section 6 identifies them; each has trace 1 − 1 = 0). ∎

**Proposition 3.2.** Given the band statement |p − q| ∈ {0, 2}, the frozen space of the whole band is ⌊N/2⌋ copies of the single SO(4) irreducible representation with

  (j_η, j_spin) = (N/2 − 1, 1),  of dimension (2j_η + 1)(2j_spin + 1) = 3(N − 1),

(at odd N the η-spin is half-integer, so strictly this is a Spin(4) = SU(2) × SU(2) representation that does not factor through SO(4) = Spin(4)/Z₂; the name follows standard Hubbard usage and the dimension count, which is what the argument uses, is the same either way)

one dimension in each of the 3(N − 1) band blocks.

*Proof.* By Lemma 3.1 each of the ⌊N/2⌋ corner seeds is η-lowest-weight; its η-weight there is (2·1 − N)/2 = 1 − N/2, so its multiplet has j_η = N/2 − 1 and spans the diagonal rungs p = 1 up to p = N − 1, which is exactly the diagonal line of the band. For the spin: by Lemma 2.5 no seed is annihilated by S⁺, so none is a spin singlet and j_spin ≥ 1 (gate check V2 additionally measures that the image is all of the frozen space of (2,0), which is the matching ceiling and is not needed for the floor); and j_spin ≥ 2 would put frozen weight on a block with |p − q| = 4, which the band statement excludes, giving j_spin = 1. Counting dimensions, ⌊N/2⌋ copies of a 3(N − 1)-dimensional irrep fill each band block with ⌊N/2⌋, which is the measured depth. ∎

Two halves of this lean on measurement and both are named. First, j_spin ≤ 1 is the band statement, which [XY_FROZEN_BAND](../../experiments/XY_FROZEN_BAND.md) establishes by a census plus a mechanism-unavailability argument, not by a proof of the null. Second, and this is the sharper one: what the argument derives is the **floor**, that ⌊N/2⌋ copies are present. That each band block carries **only** ⌊N/2⌋, which is what makes "is" rather than "contains" the right verb in the proposition, is the same per-N depth that Section 7 leaves open for all N; the proof imports it at the words "which is the measured depth". Read as a floor the proposition is derived; read as an equality it consumes the ceiling.

## 4. F141: the disagreement is η-invariant, and what that buys

**Theorem 4.1 (F141).** [𝒦, Φ] = 0 and [𝒦, Ψ] = 0 exactly, at every rung. (𝒦 is the profile-independent integer count, so the quantifier over rate profiles has nothing to bite on here; it belongs to the stronger Lemma 2.1, that the DISSIPATOR commutes with both ladders at any profile.)

*Proof.* Lemma 2.1 for Φ; Ψ is its adjoint and the disagreement count is self-adjoint. In the fermion language the same statement reads 𝒦 = n_↑ − D̂ with D̂ = Σ_l n_{l↑}n_{l↓} the double occupancy, and [n_↑, Φ] = Φ and [D̂, Φ] = Φ cancel: adding a matched pair at an empty site raises the particle number by one and the double occupancy by one, so it raises their difference by nothing. ∎

**Corollary 4.2 (the spectrum on a rung is inherited from below).** Since Φ, Ψ and 𝒦 commute in the stated way, the space decomposes into sl(2) multiplets and 𝒦 acts only on the multiplicity spaces. Writing LW_ℓ for the lowest-weight space at rung ℓ,

  **spec(𝒦 on rung p) = ⋃_{0 ≤ ℓ ≤ min(p, N−p)} spec(𝒦 on LW_ℓ),  with multiplicities.**

Both limits are load-bearing. The upper one is min(p, N−p) and **not** p: a multiplet seeded at rung ℓ has η-spin N/2 − ℓ and so spans rungs ℓ through N − ℓ, which does not reach rung p once p > N − ℓ, so the union over ℓ ≤ p overcounts on the upper half of the ladder. The lower one is 0 and **not** 1: LW_0 is the vacuum cell, whose multiplet is the identity direction climbing every rung, and it contributes one mode at 𝒦 = 0. Dropping it leaves the union short by exactly one dimension at every rung. Gate check V3(b) pins both limits rung by rung, upper half included; it reads them on V₀, which is where the reduction lives, and the limits are the same on the full rung.

**What this does and does not reduce.** It does not yet say the ceiling is a statement about lowest weights, and two things stand in the way, both of which Section 5 removes rather than this section.

- **"Frozen" is not "𝒦 = 1" in general.** At finite coupling the frozen condition is (−i·ad_h + 4γ(1 − 𝒦))v = 0, and 𝒦 does not commute with ad_h, so a frozen mode is not usually a 𝒦-eigenvector at all. The two conditions coincide only on ker(ad_h), which is where the large-J reduction of Section 5 puts the question.
- **On the full lowest-weight spaces the criterion is vacuous.** 𝒦 is diagonal on cells with integer entries and preserves ker Ψ, so the eigenvalue 1 is present there with large multiplicity at essentially every N, including every N at which the ceiling holds. Read on the full spaces, "no LW_ℓ carries 𝒦 = 1" is simply false, and false in a way that would prove nothing if it were true.

So the usable form of the reduction is the one Section 5 supplies, on **V₀ := ker(ad_h)**, and it runs one way only:

  **no LW_ℓ ∩ V₀ with 2 ≤ ℓ ≤ min(p, N−p) carries the eigenvalue 1  ⟹  the ceiling holds at rung p.**

The converse fails, and N = 5 is the counterexample: there the reduced lowest-weight space at ℓ = 2 does carry 𝒦 = 1, the first-order certificate is not full rank, and the ceiling holds anyway, closed at third order. That is the same one-directional shape as Proposition 5.1, which bounds the multiplicity from above and never from below.

Gate check V3(a) pins the commutators at machine zero on uniform and graded profiles; V3(b) pins the union identity.

## 5. The large-J reduction, and why the even orders are empty

Write the frozen question on a diagonal rung as a kernel. With h = J·h₀ and s = 1/J, and after dividing by J,

  **N(s) = A + s·B,  A = −i·ad_{h₀},  B = 4γ(1 − 𝒦),**

whose kernel for small s > 0 is the frozen space. A is anti-Hermitian and semisimple, so in the mode operator basis |u_a⟩⟨u_b| it is diagonal with entries −i(E_a − E_b), the projector P onto V₀ := ker A is orthogonal, and the reduced resolvent S is elementwise i/(E_a − E_b) off V₀.

**Proposition 5.1 (the reduction is an upper bound).** The number of eigenvalue branches of N(s) that vanish identically is at most the multiplicity of 0 in B₁ := P B P restricted to V₀. Because eigenvalue multiplicity is upper semicontinuous in the parameters and the proved floor holds at *every* J, a bound valid for all large J is a bound at generic J, and the floor then makes it an equality.

*Proof.* Standard Rayleigh-Schrödinger with A semisimple: the branches emanating from 0 are s·λ_j(B₁) + O(s²), so a branch identically zero forces λ_j(B₁) = 0. ∎

The semicontinuity step is the informal bridge, and it is not what carries the weight. The rigorous per-N form of this upper bound is the exact rank certificate of [ETA_CEILING_REDUCTION](../../experiments/ETA_CEILING_REDUCTION.md): a full column rank over GF(q) forces a full rank over ℚ, which settles the non-existence outright and needs no limiting argument. Section 5 supplies the object on which that rank is read; it does not close the bound by itself.

That is what collapses the cost. The block has dimension C(N,p)², the reduced object has dimension dim V₀ = Σ_E m_E², the number of pairs of p-subsets of modes with equal Slater energy. At N = 12 on the middle rung that is 3584 against 853776.

**Proposition 5.2 (every even order vanishes).** Let T be the transpose involution (a, b) ↦ (b, a) on the operator basis. Then B is T-even and S is T-odd, so the n-th order reduced operator has T-parity (−1)^{n−1}. The frozen vectors are symmetric matrices, hence of one T-parity, and a T-odd operator has zero matrix elements between two vectors of the same parity. Therefore **every even order vanishes identically on the surviving space**, and the next order after the first that can bite is the third.

*Proof.* In the mode basis the (a,b) row of B is R_{ab}[c,d] = (U^T diag(u_a) Dm diag(u_b) U)[c,d] with Dm the symmetric cell-level multiplier, so R_{ab}[d,c] = R_{ba}[c,d], which is T-evenness. S[c,d] = i/(E_c − E_d) is manifestly T-odd. The n-th order operator is a product of n factors B and n − 1 factors S. ∎

Gate check V6 measures both halves: the second order is the zero matrix at every N tested, and where the first order is loose the third is what closes it.

## 6. F143: the seed rung in closed form

Let M := N + 1 and v_a(l) = √(2/M)·sin(πal/M) be the modes of h, with ε_a = 2J·cos(πa/M).

**Lemma 6.1 (the reduced operator at the seed).** On rung p = 1 the reduced operator P D̂ P is, in the mode basis, the matrix

  G_{ac} = Σ_l v_a(l)²·v_c(l)² = (W^T W)_{ac},  W_{lk} := v_k(l)²,

and since W is symmetric, G = W². The frozen condition there is G = 0.

**Theorem 6.2 (F143).** Let R be the chiral involution a ↦ M − a on modes and J the all-ones N×N matrix. Then

  **G = (1/M)·(J + (I + R)/2),**

and therefore, since J and R commute,

  **spec(G) = { 0 with multiplicity ⌊N/2⌋,  1/M with multiplicity ⌈N/2⌉ − 1,  1 simple }.**

The kernel is exactly the R-odd sector, spanned by the chiral differences e_k − e_{M−k}, one per mirror pair. The gap above the frozen eigenvalue is **1/M = 1/(N+1)**, for N ≥ 3; at N = 2 the middle band is empty and the gap is 1.

*Proof.* Write W_{lk} = (2/M)sin²(πkl/M) = (1/M)(1 − cos(2πkl/M)); W is symmetric in k and l by inspection. Then

  G_{ac} = (1/M²)·Σ_{l=1}^{M−1} [1 − cos(2πal/M) − cos(2πcl/M) + cos(2πal/M)cos(2πcl/M)].

Use Σ_{l=1}^{M−1} cos(2πml/M) = M·[M | m] − 1. The first three terms give (M−1) + 1 + 1. The fourth is ½[(M·[M | a−c] − 1) + (M·[M | a+c] − 1)]. For a ≠ c and a + c ≠ M this is −1, giving G_{ac} = 1/M. For a + c = M it is M/2 − 1, giving 3/(2M). For a = c with 2a ≠ M it is likewise M/2 − 1, giving 3/(2M); for the self-paired mode a = M/2 both delta terms fire and the entry is 2/M. Those four cases are exactly the entries of (1/M)(J + (I+R)/2).

For the spectrum: R permutes coordinates and J is all-ones, so RJ = JR and the two are simultaneously diagonalisable. J has eigenvalue N on the all-ones vector, which is R-even, and 0 on its orthogonal complement. (I+R)/2 is the orthogonal projector onto the R-even sector. Hence G is (N + 1)/M = 1 on the all-ones vector, (0 + 1)/M = 1/M on the R-even vectors orthogonal to it, and 0 on the R-odd vectors. The R-odd dimension is the number of 2-cycles of R on {1, …, M−1}, which is ⌊N/2⌋ (the mode a = M/2 is fixed and exists exactly when N is odd), and the remaining R-even dimension beyond the all-ones vector is ⌈N/2⌉ − 1. ∎

Two readings worth keeping. First, the kernel is the chiral-odd sector, so the ⌊N/2⌋ frozen modes at the seed are the antisymmetric combinations of chiral mode pairs, one per pair, with the self-paired middle mode at odd N contributing nothing. That is the same count the site picture calls one frozen mode per balanced pair, and the same seat left empty that [PROOF_R90_FROZEN_DIVISOR](PROOF_R90_FROZEN_DIVISOR.md) counts as the middle. Second, the 1/(N+1) is not a property of the chain but of the transform length: M = N + 1 is the length of the discrete sine transform, and the gap is one over it.

Gate check V5 pins the identity and the whole spectrum with multiplicities, and the deep run carries it to N = 40.

## 7. What is not proved

**The upper bound at rungs ℓ ≥ 2.** By Corollary 4.2 the sufficient condition is that no LW_ℓ ∩ V₀ above the seed carries 𝒦 = 1. The measurement is sharper than that: the smallest value of 𝒦 on **LW_ℓ ∩ V₀** appears to be ℓ(N − ℓ)/(N + 1), exactly, which exceeds 1 for all N ≥ 6 and equals 1 at N = 5. The intersection with V₀ is load-bearing here as everywhere: on the full lowest-weight space the minimum is 0, as Corollary 4.2 already records, so the same sentence without it would contradict that corollary. That law is measured and not derived; [ETA_CEILING_REDUCTION](../../experiments/ETA_CEILING_REDUCTION.md) states it as a conjecture, records the evidence, gives the extremal vectors in closed form, and names the exceptional N. What this document proves is everything around it: the decomposition that makes the question a question about lowest weights, the reduction that makes it small, and the bottom rung in full.

**The off-diagonal equality.** Lemma 2.5 derives the floor ⌊N/2⌋ on (2,0) and (0,2). That those blocks carry no MORE than that is the ceiling on the off-diagonal lines, and it stays measured, as [XY_FROZEN_BAND](../../experiments/XY_FROZEN_BAND.md) states for them.

**The band statement itself.** Section 3 uses |p − q| ∈ {0, 2} as an input. That statement is measured, with a mechanism-unavailability argument for its edge, in the band note; nothing here strengthens it.
