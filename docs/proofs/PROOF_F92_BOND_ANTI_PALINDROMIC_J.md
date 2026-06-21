# PROOF F92: F71-anti-palindromic J spectral invariance (J-side Pi2-Z₄ twin of F91)

**Status:** Tier 1 derived (algebraic argument parallel to F91 + bit-exact empirical witness at N=4, 5)
**Date:** 2026-05-12
**Authors:** Thomas Wicht, Claude (Anthropic)
**Typed claim:** [`F92BondAntiPalindromicJSpectralInvariance.cs`](../../compute/RCPsiSquared.Core/SymmetryFamily/F92BondAntiPalindromicJSpectralInvariance.cs)

---

## Abstract

F71 is the spatial mirror b ↦ N−2−b acting on the bond couplings; F92 reads the anti-palindromic reshuffle of J as the parameter-side shadow of the Pi2 operator-side 90° turn. For the XY chain with inhomogeneous bond couplings J_b under uniform Z-dephasing, the F71-refined diagonal-block eigenvalue multiset is invariant not only when J is palindromic (J_b = J_{N−2−b}) but across the strictly larger anti-palindromic orbit J_b + J_{N−2−b} = 2·J_avg. The diagonal-block spectrum sees only the F71-symmetric pair-sums of J, never the anti-symmetric pair-differences; the full Liouvillian L does change (F71 breaks as an L-symmetry, the off-block Frobenius norm is nonzero, and the full-L decay spectrum genuinely moves), so the entire breaking is carried into the eigenvectors while only the diagonal-block eigenvalues stand still.

This makes F92 the J-axis twin of F91 (the same structural identity run on the dephasing rates γ) and the sibling of F93 (the h-detuning axis): one algebraic identity, instantiated on every parameter coordinate that feeds the Liouvillian. It sits strictly between F71-as-J-symmetry, which it weakens, and F1, which it strengthens, and its observable-side dual is F100, where the closure-breaking c₁/Q_peak coefficient is exactly odd in the anti-palindromic part J_anti that the diagonal-block spectrum ignores. The diagonal-block spectrum keeps the symmetric half of J; the geometry of the modes (and so the full-L spectrum) keeps the rest.

## Statement

For the chain XY + uniform Z-dephasing Liouvillian L on N qubits with **inhomogeneous bond couplings** J_b (b ∈ {0..N−2}, Hamiltonian H = (1/2) Σ_b J_b · (X_b X_{b+1} + Y_b Y_{b+1})), the **F71-refined diagonal-block eigenvalue multiset** is invariant under any J-distribution satisfying

    J_b + J_{N−2−b} = 2·J_avg = (2/(N−1))·Σ_b J_b   for all b ∈ {0..N−2}

i.e. **J is F71-anti-palindromic around its mean**. The full operator L itself generally changes (F71 broken as L-symmetry, off-block Frobenius in F71-refined basis nonzero), but the diagonal-block eigenvalues coincide; the breaking lives entirely in eigenvectors (F71-cross-blocks).

Strictly weaker than F71-as-J-symmetry (J_b = J_{N−2−b}); strictly stronger than F1 alone (which only requires Σγ_l invariant, and γ here is uniform anyway).

**Middle-bond constraint:** for N even (so N−1 bonds with odd count), the middle bond at b = (N−2)/2 is F71-fixed (b = N−2−b). Anti-palindromy then requires J_{(N−2)/2} = J_avg. For N odd (even bond count), all bonds pair up; no middle constraint.

## Parameter-side Klein V₄ (shadow of the Pi2 operator-side Z₄)

The Pi2-Foundation Z₄ (`NinetyDegreeMirrorMemoryClaim` in `Pi2KnowledgeBaseClaims.cs`) is genuinely cyclic of order 4 only on the **operator side**, where the quarter-turn is the literal i in F80's Spec(M) = ±2i·Spec(H) (i² = −1, i⁴ = 1). On the **parameter vector** J_b the realized structure is the **Klein four-group V₄ = Z₂×Z₂**, two commuting involutions, not a cyclic Z₄:

| V₄ element | Action on J_b | Order | Effect on L | Effect on F71-refined diagonal spectrum |
|---|---|---|---|---|
| e | J_b ↦ J_b | 1 | unchanged | unchanged |
| F71 (palindromic mirror) | J_b ↦ J_{N−2−b} | 2 | F71 holds as L-symmetry | unchanged |
| R₉₀ (anti-palindromic reshuffle) | J_b ↦ 2·J_avg − J_{N−2−b} | 2 | F71 broken | **unchanged** |
| F71∘R₉₀ (mean-reflection) | J_b ↦ 2·J_avg − J_b | 2 | F71 broken | unchanged |

R₉₀ is an **involution** (R₉₀² = e, NOT F71): it **preserves each F71-pair difference and reflects each pair-sum about 2·J_avg** (T_b ↦ 4·J_avg − T_b). The anti-palindromic orbit T_b = 2·J_avg is its fixed-point set (there R₉₀ acts as the identity). Diagonal-block matrix elements depend on the pair-sums only, so they are R₉₀-invariant on the orbit; cross-block elements depend on the pair-differences. R₉₀ is the order-2 parameter-side **shadow** of the operator-side 90° turn, not itself a quarter-turn; the genuine i⁴ = 1 lives on the operator side. Same structural identity as F91/F93; different parameter axis.

## Algebraic proof (parallel to PROOF_F91 § Algebraic proof)

### Setup

Same Liouville-space basis pairs |a⟩⟨b| labelled by computational basis indices (a, b) ∈ {0, 1}^N × {0, 1}^N. The F71-refined basis vectors are `|sym⟩ = (|a⟩⟨b| + |a'⟩⟨b'|)/√2` and `|antisym⟩ = (|a⟩⟨b| − |a'⟩⟨b'|)/√2` for F71-orbit-2 pairs, where (a', b') is the bit-reversed image.

### Step 1. Hamiltonian action on basis pair |a⟩⟨b|.

The XY chain Hamiltonian H = (1/2) Σ_b J_b · (X_b X_{b+1} + Y_b Y_{b+1}) acts on |a⟩⟨b| via [H, |a⟩⟨b|] = H|a⟩⟨b| − |a⟩⟨b|H. Each XY-bond b applies a swap |…01…⟩ ↔ |…10…⟩ on adjacent sites b, b+1 if the bits differ.

### Step 2. F71-action on bond index.

For the F71-image (a', b'), bit-reversal maps site index ℓ to N−1−ℓ. Therefore bond b (connecting sites b, b+1) maps to bond N−2−b (connecting sites N−1−b and N−2−b). The Hamiltonian's bond-index relabels:

    F71 · H(J) · F71 = H(F71(J))   where (F71(J))_b = J_{N−2−b}

### Step 3. F71-symmetric/anti-symmetric J-decomposition.

Since H is linear in J (H(J) = Σ_b J_b·h_b with h_b = (1/2)(X_bX_{b+1} + Y_bY_{b+1})), we can decompose:

    J_sym  := (J + F71(J)) / 2,     i.e. (J_sym)_b = T_b/2 = (J_b + J_{N−2−b})/2
    J_anti := (J − F71(J)) / 2,     i.e. (J_anti)_b = B_b/2 = (J_b − J_{N−2−b})/2

with J = J_sym + J_anti, hence H(J) = H(J_sym) + H(J_anti).

**J_sym is F71-palindromic by construction:** F71(J_sym) = (F71(J) + J)/2 = J_sym. Therefore H(J_sym) commutes with F71-conjugation: F71·H(J_sym)·F71 = H(J_sym), and the Liouvillian super-operator [H(J_sym), ·] is F71-block-diagonal in the F71-refined basis (same Step-5 argument as PROOF_F91 § Step 5 for F71-symmetric H).

**J_anti is F71-anti-palindromic** by construction: F71(J_anti) = −J_anti. Therefore F71·H(J_anti)·F71 = H(F71(J_anti)) = H(−J_anti) = −H(J_anti). The super-operator [H(J_anti), ·] anti-commutes with F71-conjugation, so in the F71-refined basis (which diagonalises F71-conjugation with eigenvalues ±1) it is **purely off-block**: it maps F71-even sub-sectors to F71-odd and vice versa, contributing only to ⟨sym|·|antisym⟩ and ⟨antisym|·|sym⟩ matrix elements, never to ⟨sym|·|sym⟩ or ⟨antisym|·|antisym⟩.

**Conclusion of Step 3:** ⟨sym|[H(J), ·]|sym⟩ = ⟨sym|[H(J_sym), ·]|sym⟩, depending on J only through (J_sym)_b = T_b/2. The same holds for ⟨antisym|·|antisym⟩. Cross-block matrix elements ⟨sym|·|antisym⟩ pick up the [H(J_anti), ·] contribution and depend on J only through (J_anti)_b = B_b/2.

The Z-dephasing dissipator is γ-uniform (γ_l = γ ∀l), so its contribution is identical across all γ-permutations (palindromic by triviality). It is a Π-commuting diagonal operator (eigenvalue −2γ·popcount(a⊕b) on |a⟩⟨b|, varying with the pair's Hamming distance), so it never mixes F71-even with F71-odd and preserves the F71-block structure.

### Step 4. Conclusion.

Diagonal-block matrix elements of L = −i[H, ·] + D_uniform-γ in the F71-refined basis are linear functionals of J depending only on the indexed F71-pair-sums T_b (equivalently, on J_sym as a vector / on H(J_sym)), not merely on their multiset: reassigning which site-pair carries which pair-sum changes the spectrum (a same-multiset, different-assignment pair at N=6 gives distinct diagonal-block spectra). The reshuffle R₉₀ : J_b ↦ 2·J_avg − J_{N−2−b} maps T_b ↦ 4·J_avg − T_b. Anti-palindromy T_b = 2·J_avg ∀b is the orbit (R₉₀'s fixed-point set) on which J_sym is the constant vector J_avg, so all J there give identical diagonal-block spectra equal to uniform J_avg (the multiset/assignment distinction collapses on a constant vector).

## Empirical witness

| N | uniform J = 1.0 | anti-palindromic J | full-L Frobenius leak (F71 off-block) |
|---|---|---|---|
| 4 | spectrum reference | bit-identical | nonzero (J_b ≠ J_{N−2−b} per pair) |
| 5 | spectrum reference | bit-identical | nonzero |

Verified by `F92BondAntiPalindromicJSpectralInvarianceTests.Spectrum_InvariantUnderAntiPalindromicJ_*` using the F71-refined diagonal-block spectrum (`F71MirrorBlockRefinement.ComputeSpectrumPerBlock`); the full-L spectrum (via `LiouvillianBlockSpectrum.ComputeSpectrumPerBlock`) does NOT match across anti-palindromic J profiles, consistent with F91's analogous caveat.

## Connection to F91

F91 (γ-axis) and F92 (J-axis) are two instances of the same Pi2 structural identity (operator-side Z₄, parameter-side V₄ shadow) applied to different parameter families. Both are corollaries of:

> Any L-parameter family P_α (α indexing the family) for which F71 conjugation acts as P ↦ F71(P), and for which the L-matrix elements in the F71-refined basis are linear in P, carries a **Klein V₄** on its parameter space, generated by the palindromic mirror F71 and the anti-palindromic involution R₉₀ : P ↦ 2·P_avg − F71(P) (two commuting order-2 maps); the F71-refined diagonal-block spectrum is invariant on the anti-palindromic orbit P_α + P_{F71(α)} = 2·P_avg ∀α (the fixed-point set of R₉₀). The genuine order-4 quarter-turn is the operator-side i⁴ = 1, of which R₉₀ is the order-2 parameter-side shadow.

This generalisation predicts F93 (h-detuning axis) automatically; whether further parameter axes apply depends only on whether F71 acts non-trivially on them.

## What this is NOT

- This is **not** a claim about the full-L spectrum. The physical Liouvillian decay rates do move on the anti-palindromic orbit: at fixed J_avg the full-L spectrum differs from the uniform-J reference (verified from below, max |Δλ| ≈ 5.8 at N=4, ≈ 8.1 at N=5). Only the F71-**compressed** diagonal blocks (P₊LP₊ ⊕ P₋LP₋, cross-terms discarded) stand still; once F71 is broken those blocks are a basis-relative quantity, not L's eigenvalues. (Same caveat as F91's "What this is NOT".)
- The genuine J-independent true-eigenvalue structure is the **joint-popcount (U(1)×U(1)) sectors**, which survive anti-palindromic J intact (J only enters popcount-preserving XY swaps); the F71 even/odd split is a within-popcount compression once F71 is broken.
- The right-sized physical reading is **first-order stationarity**: by F100's Hellmann-Feynman result δλ = 0 along J_anti, no eigenvalue moves to first order; the bit-exact all-orders invariance is a property of the compressed operator, not of L. See `./PROOF_F100_C1_QPEAK_MIRROR_J_PARITY.md`.
- This is **not** a Z₄ acting on the bond couplings: the parameter-side group is Klein V₄ (see above); the genuine order-4 quarter-turn is operator-side.

## Anchors

- Memory record: `project_anti_palindromic_is_ninety_degrees.md` (γ-side); F92 added to the inventory.
- Typed claim: `F92BondAntiPalindromicJSpectralInvariance` (Tier1Derived, two parents).
- Sister proof (γ-side): `./PROOF_F91_GAMMA_NINETY_DEGREES.md`.
- Future sister (h-side): `./PROOF_F93_DETUNING_ANTI_PALINDROMIC.md`.
- Observable-side twin: `./PROOF_F100_C1_QPEAK_MIRROR_J_PARITY.md` (F100). Same J_sym/J_anti split, read on the c₁/Q_peak bond-mirror observable instead of the diagonal-block spectrum: F92 says the spectrum depends only on J_sym, F100 says the c₁/Q_peak deviation depends only on J_anti and is exactly odd in it.
- Inventory: `docs/SYMMETRY_FAMILY_INVENTORY.md`.
- Reflection: `reflections/ON_THE_SYMMETRY_FAMILY.md`.
