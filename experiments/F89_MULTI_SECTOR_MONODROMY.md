# The Multi-Sector Monodromy Census: the S₈ braid is confined at N=4 but spreads to a 12-sector diamond at N=5

**Status:** Tier 1 Candidate. The N=4 confinement (braid on the D₄ orbit, dense core braid-free) and the N=5 spread (a 12-sector joint-popcount diamond, two cross-fold-conjugate families sharing a byte-identical eigenvalue λ) are EXACT and gate-tested (the exact 12-set is pinned as a regression test). The mechanism, free-fermion / AT additivity, is CONFIRMED for the embedding (the diagonal-spectator byte-identity and the F89d cross-fold, both to machine precision) and for the membership rule; the reduced single-particle-level object whose defective coalescence gives the λ VALUE (and the codim-1-by-additivity theorem) is the remaining from-below build.

**Date:** 2026-07-01.

**Authors:** Thomas Wicht, Claude (Opus 4.8)

## What this is about

The F89 story found that the (1,2) = (SE, DE) coherence block of the XY (Δ=0) chain under Z-dephasing carries a Galois **braid**: as the coupling q = J/γ loops the complex plane, eight eigenvalues of the block's residual octic braid, and the braiding IS the Galois group S₈ (monodromy = Galois over the function field; [F89OcticMonodromyClaim](../compute/RCPsiSquared.Core/Symmetry/F89OcticMonodromyClaim.cs)). The braid-carrying branch points are the **defective** exceptional points (Jordan blocks, eigenvectors coalesced, a loud √-branch), as opposed to the **diabolic** degeneracy (semisimple, eigenvectors independent, a silent double-discriminant-zero) at q ≈ 0.659 ([F89Path3OcticEpClaim](../compute/RCPsiSquared.Core/Symmetry/F89Path3OcticEpClaim.cs)).

The Liouvillian L(q) decomposes into (N+1)² joint-popcount sectors (p, q̃), where p is the bra (column) excitation count and q̃ the ket (row) count. The (1,2) octic lives in one such sector. **The question:** is the S₈ braid LOCALIZED to the (1,2) sector (and its symmetry images), or is it SHARED across the joint-popcount lattice — do other sectors carry the same defective branch points at the same q-loci?

The answer is **N-dependent, and that is the finding.**

## The instrument: the AT-aware sector probe

The census probes EVERY joint-popcount sector (dim > 1) at the (1,2) octic's reference DEFECTIVE loci ([ReferenceDefectiveLoci.For(N)](../compute/RCPsiSquared.Diagnostics/Foundation/ReferenceDefectiveLoci.cs), the P₁₀ simple zeros, conjugation-closed, the silent diabolic excluded), and classifies each coalescence as defective (a shared braid) or semisimple (a shared node) via the artifact-free `EpCharacter` reader.

The naive global-closest-pair probe **fails at N ≥ 5**: many raw joint-popcount blocks carry a PERMANENT semisimple degeneracy (a gap ~1e-15 eigenvalue pair at every q, sitting on the Absorption-Theorem dissipative rate lines). That permanent pair is the globally-closest pair, so the probe locks onto it, reads it diabolic, and MASKS the residual braid-carrying √-EP (gap ~5e-4) at a different λ. The fix is [`SectorEpProbe.ProbeDefectiveAnywhere`](../compute/RCPsiSquared.Diagnostics/Foundation/SectorEpProbe.cs): single-linkage-cluster every near-coalescence below a gap threshold, characterize each cluster on its own Riesz contour, and report `HasDefective` iff ANY cluster is a genuine Jordan block. The permanent semisimple cluster is skipped; the residual √-EP is reached. (The masking gotcha and the cluster fix are documented in [`SectorEpProbe`](../compute/RCPsiSquared.Diagnostics/Foundation/SectorEpProbe.cs)'s `ProbeDefectiveAnywhere` docstring.)

## The N=4 verdict: CONFINED to the D₄ orbit

At N=4 (the degenerate self-fold anchor, where the cross-fold partner w_{N−2} = w2 = DE coincides with the reference) the braid is shared across exactly the four sectors

> **{(1,2), (2,1), (2,3), (3,2)}**

— the D₄ orbit of (1,2) under transpose (p ↔ q̃) and the bra/ket bit-flips — all reading Defective at all reference loci with an identical coalescence gap (~3.7e-5). It is CONFINED there: it does **not** reach the dense half-filled core (2,2), which carries only a semisimple degeneracy. This echoes the Door-C finding that dissipative chaos is a filling threshold: the braid lives on the dilute orbit, not the dense core.

**Why N=4 is confined — the same self-fold that runs through all of F89.** N=4 is not incidentally special here: it is the F89 path-3 anchor, and its cross-fold q̃ ↦ N−2 = 2 maps (1,2) → (1,2) **onto itself** (the degenerate partner = self case, [F89CrossFoldSimilarityClaim](../compute/RCPsiSquared.Core/Symmetry/F89CrossFoldSimilarityClaim.cs) — the within-block self-fold that put one diabolic on the real axis, and the same N=4→N=5 dichotomy [F89_PATH_K_DIABOLIC](F89_PATH_K_DIABOLIC.md) reads on the diabolics). Because the cross-fold is a self-map at N=4, the two families the spread splits into (below) COINCIDE — Family B folds back onto Family A — so there is no distinct family to spread into and no dense core to reach. The confinement is that degeneracy, read on the braid.

## The N=5 verdict: SPREADS to a 12-sector diamond

At N=5, the first generic length, the braid **reaches beyond the D₄ orbit** to a symmetric 12-sector diamond:

| p \ q̃ | 0 | 1 | 2 | 3 | 4 | 5 |
|--------|---|---|---|---|---|---|
| **1** | · | n | **B** | **B** | n | · |
| **2** | n | **B** | **B** | **B** | **B** | n |
| **3** | n | **B** | **B** | **B** | **B** | n |
| **4** | · | n | **B** | **B** | n | · |

**B** = braid (defective), **n** = node-only (a semisimple/AT near-coalescence, no braid), **·** = no coalescence. The 12 braid sectors are

> **(1,2) (1,3) (2,1) (2,2) (2,3) (2,4) (3,1) (3,2) (3,3) (3,4) (4,2) (4,3)** — INCLUDING the dense diagonal cores (2,2) and (3,3).

All 12 read Defective at all 24 reference loci with a byte-identical gap signature (min 8.35e-5, max 6.12e-4). The set splits into two **cross-fold-conjugate families of 6**, each carrying a byte-identical shared eigenvalue λ (the same branch point, to 6 digits across the whole family, not merely the same gap):

| Family | Re λ range (over the 24 loci) | Sectors |
|--------|-------------------------------|---------|
| **A** | [−6.158591, −2.211959] | (1,2) (2,1) (2,3) (3,2) (3,4) (4,3) |
| **B** | [−7.788041, −3.841409] | (1,3) (2,2) (2,4) (3,1) (3,3) (4,2) |

and the two families are related by the F89d cross-fold **exactly**: `λ_B(q̄) = −conj(λ_A(q)) − 2N` (N=5). For example at q = 0.470875 − 0.047493i, (1,2) has λ = −4.028447 + 2.375482i, and −conj(λ) − 10 = −5.971553 + 2.375482i, which is exactly the (2,2) eigenvalue at q̄ ([F89CrossFoldSimilarityClaim](../compute/RCPsiSquared.Core/Symmetry/F89CrossFoldSimilarityClaim.cs), machine-zero for all N and q).

So the verdict is **N-dependent: confined at N=4, shared/coincident at N=5** — via SHARED SPECTRAL CONTENT, the same defective eigenvalue living in many joint-popcount sectors, a symmetry broader than the naive Klein-four.

## The mechanism: free-fermion / AT additivity

The membership is a clean rule:

> **braid set = {(p, q̃) : |p − q̃| = 1, both popcounts ∈ [1, N−1]} ∪ its q̃ ↦ N−q̃ cross-fold image.**

Family A is exactly the |bra − ket| = 1 sectors that keep at least one and at most N−1 excitations on each side; the |Δ| = 1 edge sectors that touch popcount 0 or N (an empty or full side) carry no coalescence — they fall in the braid-free "neither" set, its |Δ| = 1 members being (0,1),(1,0),(4,5),(5,4). Family B is the cross-fold image, not a separate phenomenon.

This single rule reproduces BOTH verdicts, and pins the N-dependence to the N=4 self-fold: at N=4 the cross-fold q̃ ↦ N−q̃ maps Family A **onto itself** (e.g. (1,2)→(1,2), (2,1)→(2,3)), so the rule collapses to the four-sector orbit and no dense core is reached — the confined verdict. At N ≥ 5 the self-fold lifts (the partner (1, N−2) ≠ (1,2)), Family B separates from Family A and brings in the dense diagonal cores, and the diamond opens — the spread verdict. **The N-dependence of the verdict is exactly the N-dependence of whether the cross-fold is a self-map on the |Δ| = 1 orbit, which it is only at N=4.**

Two operations, both verified live to machine precision, generate the whole diamond from the single elementary EP of the (1,2) SE-DE "rung":

1. **Diagonal mode-spectator (additivity).** Adding one excitation to BOTH the bra and the ket in the same single-particle mode leaves the energy difference E_bra − E_ket, and the entire EP (λ, character, gap), invariant. So (1,2) ≡ (2,3) ≡ (3,4) byte-identically: at q = 0.620878 the (1,2) and (2,3) defective λ agree to 1.2e-14. The single-particle energies are the free-fermion Bloch band E_n = 4·cos(πn/(N+1)) = {3.464, 2, 0, −2, −3.464} for N=5 (the XY chain is free-fermion integrable under Jordan-Wigner; the two-excitation energies are exact sums of one-excitation energies — [XyJordanWignerModes](../compute/RCPsiSquared.Core/F86/JordanWigner/XyJordanWignerModes.cs), [DIABOLIC_BY_INTEGRABILITY](../hypotheses/DIABOLIC_BY_INTEGRABILITY.md)).

2. **The F89d cross-fold.** q̃ ↦ N−q̃ maps Family A to Family B with λ ↦ −λ̄ − 2N (an exact antiunitary similarity, [--root crossfold](../compute/RCPsiSquared.Diagnostics/Foundation/CrossFoldSimilarityWitness.cs)).

This is the free-fermion-kernel rule (**defective = single-particle coalescence; diabolic = frequency-difference coincidence**) operating across sectors: the branch point is a property of the |Δ| = 1 coherence rung, made popcount-translation-invariant by diagonal spectators and cross-fold-invariant by F89d, so the same defective eigenvalue lives in every diamond sector. This is the sector-resolved face of the [nonhermitian_diabolic_codimension](../compute/RCPsiSquared.Core/F86/LocalGlobalEpLink.cs) result (free-fermion additivity collapses the generically codim-3 coincidence to codim-1).

## What is open

The additivity EMBEDDING is confirmed — which sectors carry the braid, and that they share the eigenvalue. But the λ VALUE itself is **not** a bare single-particle sum. At the real-q reference loci the merged eigenvalue is real with `Re λ ≈ −4.62` (q = 0.620878) and `−3.79` (q = 1.077615) — NOT the −2·integer the Absorption Theorem quantizes eigenmode rates to. The defective EP sits OFF the AT rate lines (where the permanent semisimple degeneracies live): it is a **γ-driven** coalescence, so the dephasing is essential to the λ value, not just the free-fermion frequencies. The remaining from-below build is the reduced single-particle-level object whose defective coalescence gives that λ, and the codim-1-by-additivity theorem stated formally.

## Reproduce

```bash
# The live, browsable verdict (default N=4, the fast anchor):
dotnet run --project compute/RCPsiSquared.Cli -c Release -- inspect --root sectorbraid
# The generic spread verdict (N=5, ~2 min: runs the full census live):
dotnet run --project compute/RCPsiSquared.Cli -c Release -- inspect --root sectorbraid --N 5

# The gate-tested verdicts (the exact 12-set pinned; SLOW_MSM):
dotnet test compute/RCPsiSquared.Diagnostics.Tests -c Release \
  --filter "FullyQualifiedName~MultiSectorMonodromyCensusTests"
```

The census + probe live in [`compute/RCPsiSquared.Diagnostics/Foundation/`](../compute/RCPsiSquared.Diagnostics/Foundation/) (`MultiSectorMonodromyCensus`, `SectorEpProbe`, `ReferenceDefectiveLoci`, `SectorBraidWitness`); the exact N=5 12-set is the regression test `Census_N5_BraidSpreadsBeyondOrbit_ReachesDenseCore`.
