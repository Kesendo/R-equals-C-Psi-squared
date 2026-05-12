# PROOF F92: F71-anti-palindromic J spectral invariance (J-side Pi2-Z₄ twin of F91)

**Status:** Tier 1 derived (algebraic argument parallel to F91 + bit-exact empirical witness at N=4, 5)
**Date:** 2026-05-12
**Authors:** Thomas Wicht, Claude (Anthropic)
**Typed claim:** [`F92BondAntiPalindromicJSpectralInvariance.cs`](../../compute/RCPsiSquared.Core/SymmetryFamily/F92BondAntiPalindromicJSpectralInvariance.cs)

---

## Statement

For the chain XY + uniform Z-dephasing Liouvillian L on N qubits with **inhomogeneous bond couplings** J_b (b ∈ {0..N−2}, Hamiltonian H = (1/2) Σ_b J_b · (X_b X_{b+1} + Y_b Y_{b+1})), the **F71-refined diagonal-block eigenvalue multiset** is invariant under any J-distribution satisfying

    J_b + J_{N−2−b} = 2·J_avg = (2/(N−1))·Σ_b J_b   for all b ∈ {0..N−2}

i.e. **J is F71-anti-palindromic around its mean**. The full operator L itself generally changes (F71 broken as L-symmetry, off-block Frobenius in F71-refined basis nonzero), but the diagonal-block eigenvalues coincide; the breaking lives entirely in eigenvectors (F71-cross-blocks).

Strictly weaker than F71-as-J-symmetry (J_b = J_{N−2−b}); strictly stronger than F1 alone (which only requires Σγ_l invariant, and γ here is uniform anyway).

## Pi2-Z₄ structure (parameter-side, J-axis)

The Pi2-Foundation Z₄ rotational axis (typed in `Pi2KnowledgeBaseClaims.cs` as `NinetyDegreeMirrorMemoryClaim`) has four group elements `{e, i, i², i³}` from `i⁴ = 1`. F91 typed the γ-axis Z₄; F92 types the J-axis Z₄.

| Z₄ element | Action on J_b | Effect on L | Effect on F71-refined diagonal spectrum |
|---|---|---|---|
| e | J_b ↦ J_b | unchanged | unchanged |
| i² (180°, F71-palindromic) | J_b ↦ J_{N−2−b} | F71 holds as L-symmetry | unchanged |
| i (90°, F71-anti-palindromic) | J_b ↦ 2·J_avg − J_{N−2−b} | F71 broken | **unchanged** |
| i³ (270°) | composition | F71 broken | unchanged |

The 90° rotation flips the F71-pair *difference* but preserves the F71-pair *sum*. Diagonal-block matrix elements depend on pair-sums only; cross-block on pair-differences. Same structural identity as F91; different parameter axis.

## Algebraic proof (parallel to PROOF_F91 § Algebraic proof)

### Setup

Same Liouville-space basis pairs |a⟩⟨b| labelled by computational basis indices (a, b) ∈ {0, 1}^N × {0, 1}^N. The F71-refined basis vectors are `|sym⟩ = (|a⟩⟨b| + |a'⟩⟨b'|)/√2` and `|antisym⟩ = (|a⟩⟨b| − |a'⟩⟨b'|)/√2` for F71-orbit-2 pairs, where (a', b') is the bit-reversed image.

### Step 1. Hamiltonian action on basis pair |a⟩⟨b|.

The XY chain Hamiltonian H = (1/2) Σ_b J_b · (X_b X_{b+1} + Y_b Y_{b+1}) acts on |a⟩⟨b| via [H, |a⟩⟨b|] = H|a⟩⟨b| − |a⟩⟨b|H. Each XY-bond b applies a swap |…01…⟩ ↔ |…10…⟩ on adjacent sites b, b+1 if the bits differ.

### Step 2. F71-action on bond index.

For the F71-image (a', b'), bit-reversal maps site index ℓ to N−1−ℓ. Therefore bond b (connecting sites b, b+1) maps to bond N−2−b (connecting sites N−1−b and N−2−b). The Hamiltonian's bond-index relabels:

    F71 · H(J) · F71 = H(F71(J))   where (F71(J))_b = J_{N−2−b}

### Step 3. F71-refined basis matrix elements.

Following the same calculation pattern as PROOF_F91 §§3–4:

⟨sym|H|sym⟩ involves contributions from H acting on |a⟩⟨b| and on |a'⟩⟨b'|, plus cross terms ⟨a,b|H|a',b'⟩ which vanish for non-overlapping (a, b) and (a', b') (or contribute symmetrically). After the F71-orbit averaging, **the H-contribution to ⟨sym|H|sym⟩ depends on J only through pair-sums** T_b := J_b + J_{N−2−b}.

The off-diagonal ⟨sym|H|antisym⟩ depends on pair-differences B_b := J_b − J_{N−2−b}.

The Z-dephasing dissipator is γ-uniform (γ_l = γ ∀l), so its contribution is identical across all γ-permutations (palindromic by triviality), and contributes the same diagonal shift to every F71-block.

### Step 4. Conclusion.

Diagonal-block matrix elements of L = −i[H, ·] + D_uniform-γ in the F71-refined basis are linear functionals of J depending only on the multiset of F71-pair-sums {T_b}. The 90°-rotation J_b ↦ 2·J_avg − J_{N−2−b} maps T_b ↦ 4·J_avg − T_b. Anti-palindromy T_b = 2·J_avg ∀b is the orbit on which all J give identical diagonal-block spectra equal to uniform J_avg.

## Empirical witness

| N | uniform J = 1.0 | anti-palindromic J | full-L Frobenius leak (F71 off-block) |
|---|---|---|---|
| 4 | spectrum reference | bit-identical | nonzero (J_b ≠ J_{N−2−b} per pair) |
| 5 | spectrum reference | bit-identical | nonzero |

Verified by `F92BondAntiPalindromicJSpectralInvarianceTests.Spectrum_InvariantUnderAntiPalindromicJ_*` using the F71-refined diagonal-block spectrum (`F71MirrorBlockRefinement.ComputeSpectrumPerBlock`); the full-L spectrum (via `LiouvillianBlockSpectrum.ComputeSpectrumPerBlock`) does NOT match across anti-palindromic J profiles, consistent with F91's analogous caveat.

## Connection to F91

F91 (γ-axis) and F92 (J-axis) are two instances of the same Pi2-Z₄ structural identity applied to different parameter families. Both are corollaries of:

> Any L-parameter family P_α (α indexing the family) for which F71 conjugation acts as P ↦ F71(P), and for which the L-matrix elements in the F71-refined basis are linear in P, has a Z₄ orbit `{e, i² = palindromic, i = anti-palindromic, i³}` on its parameter space, and the F71-refined diagonal-block spectrum is invariant on the anti-palindromic orbit (P_α + P_{F71(α)} = 2·P_avg ∀α).

This generalisation predicts F93 (h-detuning axis) automatically; whether further parameter axes apply depends only on whether F71 acts non-trivially on them.

## Anchors

- Memory record: `project_anti_palindromic_is_ninety_degrees.md` (γ-side); F92 added to the inventory.
- Typed claim: `F92BondAntiPalindromicJSpectralInvariance` (Tier1Derived, two parents).
- Sister proof (γ-side): `docs/proofs/PROOF_F91_GAMMA_NINETY_DEGREES.md`.
- Future sister (h-side): `docs/proofs/PROOF_F93_DETUNING_ANTI_PALINDROMIC.md`.
- Inventory: `docs/SYMMETRY_FAMILY_INVENTORY.md`.
- Reflection: `reflections/ON_THE_SYMMETRY_FAMILY.md`.
