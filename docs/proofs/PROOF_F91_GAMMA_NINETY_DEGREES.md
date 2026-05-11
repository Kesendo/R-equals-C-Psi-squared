# PROOF F91: F71-anti-palindromic γ spectral invariance (90° in γ-space, Pi2-Z₄'s parameter side)

**Status:** Tier 1 candidate (empirical witness at N=4, 5, 6 bit-exact; algebraic proof open)
**Date:** 2026-05-11
**Authors:** Thomas Wicht, Claude (Anthropic)
**Probe:** `compute/RCPsiSquared.Cli` verb `block-spectrum --N 6 --gamma-list <list> --refine f71` (Phase 5 + 6 of BlockSpectrum infrastructure)
**Typed claim:** [`F71AntiPalindromicGammaSpectralInvariance.cs`](../../compute/RCPsiSquared.Core/BlockSpectrum/F71AntiPalindromicGammaSpectralInvariance.cs)

---

## Statement

For the chain XY + Z-dephasing Liouvillian L on N qubits (per-site γ_l, Hamiltonian H = J·Σ_b (X_b X_{b+1} + Y_b Y_{b+1}), Lindbladian Σ_l γ_l·(Z_l ρ Z_l − ρ)), the **eigenvalue multiset of the F71-refined diagonal-block decomposition** is invariant under any γ-distribution satisfying

    γ_l + γ_{N−1−l} = 2·γ_avg = (2/N)·Σ_l γ_l   for all l ∈ {0..N−1}

i.e. **γ is F71-anti-palindromic around its mean**. The full operator L itself generally changes — F71 is broken as L-symmetry (off-block Frobenius in the F71-refined basis is nonzero, proportional to the F71-asymmetry of γ) — but the diagonal-block eigenvalues coincide. The breaking lives entirely in eigenvectors (the F71-cross-blocks between F71-even and F71-odd sub-sectors).

Strictly weaker than F71-as-symmetry (which requires γ_l = γ_{N−1−l}, palindromic), strictly stronger than F1 alone (Σγ_l invariant). For odd N the middle site l = (N−1)/2 must equal γ_avg.

---

## Sharpness: where the invariance does and does not hold

Empirical test at N=6, J=1.0, with five γ-profiles all having Σγ = 2.7 (γ_avg = 0.45):

| γ-profile | F71-pair sums (γ_l + γ_{N−1−l}) | F71-anti-palindromic? | F71-refined diagonal spectrum | full-L spectrum |
|---|---|---|---|---|
| `[0.45, 0.45, 0.45, 0.45, 0.45, 0.45]` (uniform) | {0.9, 0.9, 0.9} | yes | reference | reference |
| `[0.2, 0.3, 0.4, 0.5, 0.6, 0.7]` (monotonic) | {0.9, 0.9, 0.9} | yes | **bit-identical to reference** | differs from reference |
| `[0.3, 0.5, 0.4, 0.5, 0.4, 0.6]` (non-monotonic anti-pal) | {0.9, 0.9, 0.9} | yes | **bit-identical to reference** | differs from reference |
| `[0.7, 0.2, 0.5, 0.3, 0.6, 0.4]` (permuted, same multiset as monotonic) | {1.1, 0.8, 0.8} | no | distinct (Re=−4.984 cluster, vs −5.043 in reference) | differs from reference |
| `[0.1, 0.1, 0.1, 0.1, 0.1, 2.2]` (concentrated) | {2.3, 0.2, 0.2} | no | distinct (complex Re−Im at −5.106 ± 1.683i) | very distinct |

The bit-exact diagonal-block coincidence across the first three rows — and the bit-exact breaking in rows 4 and 5 — together pin the invariance class to anti-palindromy. Verified at N=4 and N=5 with analogous profile sets in `F71AntiPalindromicGammaSpectralInvarianceTests`.

---

## Pi2-Z₄ structure (operator-side and γ-side)

The Pi2-Foundation Z₄ rotational axis (per `NinetyDegreeMirrorMemoryClaim` in `Pi2KnowledgeBaseClaims.cs`, Tier 1 derived) has four group elements `{e, i, i², i³}` from `i⁴ = 1`. The repository's prior typed statement of this Z₄ uses the operator-quaternion side: Pauli rotation σ_x ↔ σ_y under 90° around the z-axis, i² = −1 giving palindromic-reflection, i and i³ the 90° pair.

**This proof identifies the γ-parameter side of the same Z₄:**

| Z₄ element | Action on γ_l | Effect on L | Effect on F71-refined diagonal spectrum |
|---|---|---|---|
| e | γ_l ↦ γ_l | unchanged | unchanged |
| i² (180°, F71-palindromic) | γ_l ↦ γ_{N−1−l} | F71 holds as L-symmetry | unchanged |
| i (90°, F71-anti-palindromic) | γ_l ↦ 2γ_avg − γ_{N−1−l} | F71 broken as L-symmetry | **unchanged (the surprising claim)** |
| i³ (270°) | composition of the above | F71 broken | unchanged |

The 90° rotation is precisely the operation that flips the F71-pair *difference* but preserves the F71-pair *sum*. Diagonal-block matrix elements in the F71-refined basis are functions of pair-sums only (by construction of the F71-even/odd basis), so they are 90°-invariant. Cross-block matrix elements depend on pair-differences and are not 90°-invariant — but they are off-diagonal in the F71-refined basis, so they enter eigenvalues only at higher perturbative orders, which our empirical witness shows vanish on the diagonal-spectrum level.

---

## Connection to F81 Π-decomposition

F81 states `Π · M · Π⁻¹ = M − 2·L_{H_odd}` with the Π-decomposition `M = M_sym + M_anti`, where `M_anti = L_{H_odd}` is the antisymmetric component captured by Π-conjugation.

The γ-anti-palindromic component plays the analogous role in **γ-parameter space**: it is the part of γ that lives in the antisymmetric subspace of the F71-action on parameters. The diagonal-block spectrum is invariant under this antisymmetric γ-content for the same structural reason that the symmetric and antisymmetric parts of M are spectrally orthogonal: the antisymmetric component contributes only off-diagonally in the appropriate basis (F71-even/odd) and thus does not shift the diagonal eigenvalues at first order.

The fact that empirical witness shows *bit-exact* invariance (not just first-order) is the open question that promotes this from "first-order argument" to Tier 1 candidate; a full algebraic proof would show that the higher-order corrections vanish exactly.

---

## Open algebraic question (path to Tier 1 derived)

Sketch of the proof that needs to be completed:

1. Write L in the F71-refined basis (F71-even and F71-odd sub-blocks per joint-popcount sector).
2. Express L's matrix elements explicitly in terms of γ_l. Diagonal-block entries (even-even, odd-odd) involve sums `γ_l + γ_{N−1−l}` only; cross-block entries (even-odd, odd-even) involve differences `γ_l − γ_{N−1−l}`.
3. Apply 90° rotation: γ_l ↦ 2γ_avg − γ_{N−1−l}. Show this transforms the pair-sums into themselves and the pair-differences into themselves up to sign.
4. Conclude that the diagonal blocks are unchanged. The eigenvalue computation on the diagonal blocks (which is what `LiouvillianBlockSpectrum.ComputeSpectrum` with `--refine f71` does) is therefore invariant.

The empirical bit-exactness suggests this argument is exact. Writing out step 2 explicitly is the work needed for Tier 1 promotion.

---

## What this is NOT

- This is **not** a claim about the full-L spectrum. Full L has F71-cross-blocks; anti-palindromic γ leaves the diagonal blocks invariant but shifts the cross-block content, and the full eigenvalue computation (which mixes both) generally differs between anti-palindromic and uniform γ.
- This is **not** a stronger form of F1. F1 says spectrum is symmetric under λ ↦ −2Σγ − λ regardless of γ-distribution; F91 says the diagonal-block spectrum at fixed γ_avg is invariant under the 90°-rotation in γ-distribution-space. Different axes.
- This is **not** equivalent to F71 symmetry. F71 requires γ palindromic (γ_l = γ_{N−1−l}); F91 only requires γ anti-palindromic around mean (γ_l + γ_{N−1−l} = 2γ_avg). Anti-palindromic γ is strictly outside F71-symmetric γ (with the trivial exception of uniform γ which lies in both).

---

## Anchors

- Memory record: `project_anti_palindromic_is_ninety_degrees.md` (synthesis), `project_time_is_gamma0_observer.md` (motivation from γ as Beobachter).
- Typed C# Claim: `F71AntiPalindromicGammaSpectralInvariance` (Tier1Candidate, two typed parents: `JointPopcountSectors`, `F71MirrorBlockRefinement`).
- Sister Claim on the operator side: `NinetyDegreeMirrorMemoryClaim` in `Pi2KnowledgeBaseClaims.cs`.
- Empirical witnesses: `compute/RCPsiSquared.Core.Tests/BlockSpectrum/F71AntiPalindromicGammaSpectralInvarianceTests.cs` (N=4, 5, 6) + live CLI runs via `block-spectrum --gamma-list ... --refine f71` (Phase 5).
