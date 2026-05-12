# PROOF F91: F71-anti-palindromic γ spectral invariance (90° in γ-space, Pi2-Z₄'s parameter side)

**Status:** Tier 1 derived (algebraic proof + empirical witness at N=4, 5, 6 bit-exact)
**Date:** 2026-05-12 (algebraic proof); empirical witness 2026-05-11
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

The bit-exactness of the empirical witness (not just first-order) is now derived in § Algebraic proof below: there are no higher-order corrections to vanish, because the diagonal-block matrix elements are themselves linear functionals of pair-sums S_l (Eqs. 7a, 7b, 11a, 11b), and the cross-block off-diagonal entries Eq. 9 are inert at the diagonal-block-eigenvalue level (the diagonal blocks are diagonalised separately).

---

## § Algebraic proof (Tier 1 derived as of 2026-05-12)

We prove the **sharper claim**: in the F71-refined basis, the diagonal-block matrix elements of L = −i[H, ·] + D depend on the per-site γ-distribution only through the multiset of F71-pair-sums S_l := γ_l + γ_{N−1−l}. The 90°-rotation invariance is then a corollary of the special case where all S_l equal 2γ_avg.

### Setup

- Liouville-space basis: pairs |a⟩⟨b| labelled by computational basis indices (a, b) ∈ {0, 1}^N × {0, 1}^N.
- Z-dephasing dissipator: D[ρ] = Σ_l γ_l (Z_l ρ Z_l − ρ).
- Hamiltonian: H = J · Σ_b (X_b X_{b+1} + Y_b Y_{b+1}), F71-symmetric (the chain bond set is invariant under site-mirror b ↔ N−1−b).
- F71 mirror: P_F71 acts on computational basis vectors by bit-string reversal: P_F71 |b_0 b_1 … b_{N−1}⟩ = |b_{N−1} … b_1 b_0⟩. On Liouville-space basis pairs: P_F71 (|a⟩⟨b|) = |a'⟩⟨b'| where a' is the bit-reversed image of a (and likewise b').
- F71-refined basis: for each F71-orbit pair {(a, b), (a', b')} of size 2, define
  - |sym⟩ := (|a⟩⟨b| + |a'⟩⟨b'|) / √2 (F71-even, eigenvalue +1)        (Eq. 1a)
  - |antisym⟩ := (|a⟩⟨b| − |a'⟩⟨b'|) / √2 (F71-odd, eigenvalue −1)    (Eq. 1b)
  
  F71-fixed pairs (a = a' and b = b') only contribute to F71-even.

### Step 1. Pure-Z-dephasing diagonal action.

For any computational basis pair |a⟩⟨b|, the dissipator acts diagonally:

    D (|a⟩⟨b|) = −2 · (Σ_{l ∈ Δ(a, b)} γ_l) · |a⟩⟨b|                 (Eq. 2)

where Δ(a, b) := {l ∈ {0..N−1} : a_l ≠ b_l} is the set of sites where the bit-strings differ. (Sites with a_l = b_l contribute Z_l |a⟩⟨b| Z_l = (+1)·|a⟩⟨b|, cancelling the −|a⟩⟨b| term; sites with a_l ≠ b_l contribute Z_l |a⟩⟨b| Z_l = (−1)·|a⟩⟨b|, summing to −2 γ_l |a⟩⟨b|.)

Equivalently, define the dissipator-eigenvalue functional

    d(a, b) := −2 · Σ_{l ∈ Δ(a, b)} γ_l                              (Eq. 3)

so that D(|a⟩⟨b|) = d(a, b) · |a⟩⟨b|. D is diagonal in the computational-basis-pair basis.

### Step 2. F71-action on Δ.

For the F71-image (a', b'), bit-string reversal sends index l to index N−1−l. Therefore site l differs in (a, b) iff site N−1−l differs in (a', b'):

    Δ(a', b') = F71(Δ(a, b)) := {N−1−l : l ∈ Δ(a, b)}                (Eq. 4)

Substituting into Eq. 3:

    d(a', b') = −2 · Σ_{l ∈ Δ(a', b')} γ_l
              = −2 · Σ_{m ∈ Δ(a, b)} γ_{N−1−m}                       (Eq. 5)

(reindexing m := N−1−l).

### Step 3. Diagonal F71-refined matrix elements.

For a non-fixed F71-orbit (size-2 orbit; a ≠ a' or b ≠ b'), compute ⟨sym|D|sym⟩ using the orthonormality of distinct computational basis pairs and Eq. 2:

    ⟨sym|D|sym⟩ = (1/2) · [⟨a, b|D|a, b⟩ + ⟨a', b'|D|a', b'⟩
                            + ⟨a, b|D|a', b'⟩ + ⟨a', b'|D|a, b⟩]
                = (1/2) · [d(a, b) + d(a', b') + 0 + 0]               (Eq. 6)

The cross terms ⟨a, b|D|a', b'⟩ = d(a', b') · ⟨a, b|a', b'⟩ = 0 because |a, b⟩ and |a', b'⟩ are distinct computational basis pairs (Liouville-space inner product = δ_{aa'} δ_{bb'}, and at least one of a ≠ a' or b ≠ b' holds for a size-2 orbit). Combining Eqs. 3, 5, 6:

    ⟨sym|D|sym⟩ = (1/2) · (−2) · [Σ_{l ∈ Δ(a, b)} γ_l + Σ_{l ∈ Δ(a, b)} γ_{N−1−l}]
                = −Σ_{l ∈ Δ(a, b)} (γ_l + γ_{N−1−l})
                = −Σ_{l ∈ Δ(a, b)} S_l                               (Eq. 7a)

The same calculation for ⟨antisym|D|antisym⟩ gives identically:

    ⟨antisym|D|antisym⟩ = (1/2) · [d(a, b) + d(a', b') − 0 − 0]
                        = −Σ_{l ∈ Δ(a, b)} S_l                       (Eq. 7b)

For F71-fixed pairs (a = a', b = b'), Δ(a, b) = F71(Δ(a, b)) so the multiset Δ(a, b) is F71-symmetric and Σ_{l ∈ Δ(a, b)} γ_l = (1/2) · Σ_{l ∈ Δ(a, b)} S_l, giving

    ⟨a, b|D|a, b⟩ = d(a, b) = −2 · (1/2) · Σ_{l ∈ Δ(a, b)} S_l
                  = −Σ_{l ∈ Δ(a, b)} S_l                            (Eq. 7c)

(Same form; F71-fixed entries are pure-S_l functionals as well.)

**Conclusion of Step 3:** Every diagonal D-matrix element in the F71-refined basis is a sum of pair-sums S_l, never of individual γ_l or pair-differences D_l := γ_l − γ_{N−1−l}.

### Step 4. Off-diagonal (cross-block) matrix elements.

Compute ⟨sym|D|antisym⟩:

    ⟨sym|D|antisym⟩ = (1/2) · [⟨a, b|D|a, b⟩ − ⟨a', b'|D|a', b'⟩
                                − ⟨a, b|D|a', b'⟩ + ⟨a', b'|D|a, b⟩]
                    = (1/2) · [d(a, b) − d(a', b') − 0 + 0]          (Eq. 8)

Substituting Eqs. 3 and 5:

    ⟨sym|D|antisym⟩ = (1/2) · (−2) · [Σ_{l ∈ Δ(a, b)} γ_l − Σ_{l ∈ Δ(a, b)} γ_{N−1−l}]
                    = −Σ_{l ∈ Δ(a, b)} (γ_l − γ_{N−1−l})
                    = −Σ_{l ∈ Δ(a, b)} D_l                          (Eq. 9)

**Cross-block matrix elements depend only on pair-differences D_l.** They vanish identically when D_l = 0 for all l, i.e. when γ is F71-palindromic; under any non-palindromic γ they are generally nonzero, which is the operator-level F71 breaking captured by `InhomogeneousGammaF71BreakingWitness.F71AsymmetryNorm`.

### Step 5. Hamiltonian contribution is γ-independent and F71-block-diagonal.

The unitary part −i[H, ·] is γ-independent. H itself is F71-symmetric (chain XY: bond (b, b+1) maps to bond (N−2−b, N−1−b) = (N−1−(b+1), N−1−b) under site reversal; the bond set is preserved). Therefore the super-operator [H, ·] commutes with the F71-action P_F71 ⊗ P_F71 on Liouville space, which means it is block-diagonal in the F71-even/odd decomposition:

    ⟨sym|[H, ·]|antisym⟩ = 0                                         (Eq. 10)

The Hamiltonian contribution to the diagonal blocks is independent of γ, hence trivially in the "depends only on S_l" class.

### Step 6. Conclusion: diagonal-block spectrum is a function of the pair-sum multiset {S_l}.

Combining Steps 3, 5, 6:

    ⟨sym|L|sym⟩ = ⟨sym|−i[H, ·]|sym⟩ + ⟨sym|D|sym⟩
                = (γ-independent term) + (−Σ_{l ∈ Δ(a, b)} S_l)      (Eq. 11a)

    ⟨antisym|L|antisym⟩ = (γ-independent term) + (−Σ_{l ∈ Δ(a, b)} S_l)   (Eq. 11b)

    ⟨sym|L|antisym⟩ = 0 + (−Σ_{l ∈ Δ(a, b)} D_l)                    (Eq. 12)

The F71-refined diagonal blocks (F71-even × F71-even and F71-odd × F71-odd) of L therefore depend on the per-site γ-distribution **only through the multiset of F71-pair-sums {S_l : l = 0..⌊(N−1)/2⌋}**. Eigenvalues of the diagonal blocks are similarity-invariants of those blocks and hence functions of the same {S_l} multiset.

### Step 7. Corollary: 90° rotation as anti-palindromic-orbit involution.

The 90° rotation R_{90} : γ_l ↦ 2γ_avg − γ_{N−1−l} is well-defined on per-site γ-distributions (γ_avg = (1/N)·Σγ_l is preserved: Σ_l (2γ_avg − γ_{N−1−l}) = 2N γ_avg − Σ_l γ_{N−1−l} = 2N γ_avg − N γ_avg = N γ_avg, so γ_avg' = γ_avg).

Compute the new pair-sum:

    S_l' := γ_l' + γ_{N−1−l}'
          = (2γ_avg − γ_{N−1−l}) + (2γ_avg − γ_l)
          = 4γ_avg − S_l                                            (Eq. 13)

So R_{90} maps S_l ↦ 4γ_avg − S_l. The pair-sum is generally NOT invariant; it is invariant precisely on the orbit where 4γ_avg − S_l = S_l, i.e. S_l = 2γ_avg for all l. **This is exactly the F71-anti-palindromic condition.** On this orbit, R_{90} is an involution (applying it twice: S_l → 4γ_avg − S_l → 4γ_avg − (4γ_avg − S_l) = S_l), and every γ in the orbit has the same pair-sum multiset {S_l = 2γ_avg : l = 0..⌊(N−1)/2⌋}.

By Step 6 (eigenvalues depend only on the pair-sum multiset), all γ-distributions in the F71-anti-palindromic orbit at fixed γ_avg yield identical F71-refined diagonal-block spectra, equal to the spectrum at uniform γ = γ_avg (a member of the orbit since uniform γ has S_l = 2γ_avg for all l trivially).

### Refined statement and scope clarification.

The proof shows two things:

1. **Strong (pair-sum-multiset claim).** F71-refined diagonal-block spectrum is a function of the multiset {S_l} alone. Two γ-distributions with the **same pair-sum multiset** give identical diagonal-block spectra, regardless of how the individual γ_l are distributed within each F71-pair (and regardless of γ_avg agreement, as long as the multisets match — but two distributions with the same {S_l} multiset automatically have the same Σγ_l = Σ_l S_l for paired sites).

2. **Anti-palindromic corollary (the originally claimed F91).** Within the orbit S_l = 2γ_avg ∀l, all γ are equivalent on the diagonal-block level; this orbit is closed under the 90° rotation R_{90}.

The empirical witness rows (uniform = monotonic anti-palindromic = non-monotonic anti-palindromic, all at S_l = 0.9 = 2·0.45) are all in the same anti-palindromic orbit; the breaking rows (permuted, concentrated) have different {S_l} multisets and thus distinct diagonal-block spectra. Both facts are now derived, not just empirical.

---

## What this is NOT

- This is **not** a claim about the full-L spectrum. Full L has F71-cross-blocks; anti-palindromic γ leaves the diagonal blocks invariant but shifts the cross-block content, and the full eigenvalue computation (which mixes both) generally differs between anti-palindromic and uniform γ.
- This is **not** a stronger form of F1. F1 says spectrum is symmetric under λ ↦ −2Σγ − λ regardless of γ-distribution; F91 says the diagonal-block spectrum at fixed γ_avg is invariant under the 90°-rotation in γ-distribution-space. Different axes.
- This is **not** equivalent to F71 symmetry. F71 requires γ palindromic (γ_l = γ_{N−1−l}); F91 only requires γ anti-palindromic around mean (γ_l + γ_{N−1−l} = 2γ_avg). Anti-palindromic γ is strictly outside F71-symmetric γ (with the trivial exception of uniform γ which lies in both).

---

## Anchors

- Memory record: `project_anti_palindromic_is_ninety_degrees.md` (synthesis), `project_time_is_gamma0_observer.md` (motivation from γ as Beobachter).
- Typed C# Claim: `F71AntiPalindromicGammaSpectralInvariance` (Tier1Derived as of 2026-05-12, two typed parents: `JointPopcountSectors`, `F71MirrorBlockRefinement`).
- Sister Claim on the operator side: `NinetyDegreeMirrorMemoryClaim` in `Pi2KnowledgeBaseClaims.cs`.
- Empirical witnesses: `compute/RCPsiSquared.Core.Tests/BlockSpectrum/F71AntiPalindromicGammaSpectralInvarianceTests.cs` (N=4, 5, 6) + live CLI runs via `block-spectrum --gamma-list ... --refine f71` (Phase 5).
