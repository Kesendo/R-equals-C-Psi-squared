# The Multi-Sector Monodromy Census: the S₈ braid is confined at N=4 but spreads to a 12-sector diamond at N=5

**Status:** Tier 1 Candidate. The N=4 confinement (braid on the D₄ orbit, dense core braid-free) and the N=5 spread (a 12-sector joint-popcount diamond, two cross-fold-conjugate families sharing a byte-identical eigenvalue λ) are EXACT and gate-tested (the exact 12-set is pinned as a regression test). The mechanism, free-fermion / AT additivity, is CONFIRMED for the embedding (the diagonal-spectator byte-identity and the F89d cross-fold, both to machine precision) and for the membership rule; the λ VALUE is grounded as the AT Theorem-2 rate of the mixed defective eigenvector (Re λ = −2·⟨n_diff⟩ at the real loci, machine-exact). The closed-form mixture ⟨n_diff⟩(q) is now RESOLVED via the free-fermion mode geometry: the mode-density overlap I(a,b) is exactly quantized (1/(N+1) or 3/(2(N+1))), N=4 pins ⟨n_diff⟩ = N/2 by a vanishing off-diagonal mixing (a second proof beside the Re λ = −N self-fold axis), and N ≥ 5 is non-radical because the function ⟨n_diff⟩(q) carries √-branch points at the S₈ exceptional points. The general-N membership and the codim-1-by-additivity theorem stay open.

**Date:** 2026-07-01.

**Authors:** Thomas Wicht, Claude (Opus 4.8)

## Setup (plain words, for a first-time reader)

The system is an **N-site spin chain**: N qubits in a row, neighbours coupled by the XY interaction (strength J, no ZZ term, so anisotropy Δ = 0), and each qubit independently **dephased** by local Z-noise (rate γ). This is an OPEN quantum system, so its density matrix ρ evolves by a Lindblad master equation whose generator is the **Liouvillian L**: a matrix acting on the space of operators (on ρ itself, dimension 4^N). The single dimensionless knob is **q = J/γ**, and the analysis continues q into the complex plane.

Because both the XY coupling and Z-dephasing conserve the number of up-spins, L splits into blocks labelled by two excitation counts. A **joint-popcount sector (p, q̃)** is the block of coherences |a⟩⟨b| whose bra index a has p up-spins (popcount p) and whose ket index b has q̃ up-spins. The **(1,2) = (SE, DE) coherence block** is thus single-excitation-bra against double-excitation-ket. (Notation warning: q̃, a popcount label, is NOT the coupling q.)

As q loops in the complex plane the eigenvalues of a block move and can collide. A **defective** collision (an *exceptional point*: a Jordan block, eigenvectors coalesced, a √-branch point) is genuine and "loud"; a **diabolic** collision (a plain repeated eigenvalue, eigenvectors still independent) is silent. Now trace q around a closed loop that encircles the defective points and returns to its start: the block's eigenvalues come back **permuted** among themselves. The set of all permutations reachable this way is the **monodromy**; for the (1,2) block's characteristic polynomial (degree 8 after removing a trivially-solvable, rate-locked part, "**the octic**") the monodromy is the full symmetric group **S₈** (every reshuffle of the 8 roots is achievable, and it equals the polynomial's Galois group). That S₈ structure is what we call "**the braid**." **This experiment asks whether the braid lives only in the (1,2) block or is shared across the whole (p, q̃) sector lattice.**

One fact used throughout is the **Absorption Theorem (AT)**: the decay rate of a coherence |a⟩⟨b| (the real part of its Liouvillian eigenvalue) is fixed purely by how many sites the two indices differ on, **Re λ = −2γ·⟨n_diff⟩** with n_diff = popcount(a⊕b) the Hamming distance. A pure eigenmode has an integer n_diff; a mixed (exceptional-point) eigenvector has a fractional average ⟨n_diff⟩, which is exactly the lever the λ-value result below turns on.

One symmetry recurs throughout: the **cross-fold** flips the ket count, q̃ ↦ N − q̃ (a particle-hole flip on the ket side, swapping its excited sites for the empty ones). It is an exact structural map between the (p, q̃) and (p, N − q̃) blocks that sends an eigenvalue λ to −λ̄ − 2N, so a defective point in one block forces a partner in the other. When N − q̃ = q̃ the fold sends a block to ITSELF (a **self-fold**); this happens for the (1,2) block exactly at N = 4, and that degeneracy is the reason the N = 4 case turns out special (it is why the braid confines there but spreads at N ≥ 5).

## What this is about

The F89 story found that the (1,2) = (SE, DE) coherence block of the XY (Δ=0) chain under Z-dephasing carries a Galois **braid**: as the coupling q = J/γ loops the complex plane, eight eigenvalues of the block's residual octic braid, and the braiding IS the Galois group S₈ (monodromy = Galois over the function field; [F89OcticMonodromyClaim](../compute/RCPsiSquared.Core/Symmetry/F89OcticMonodromyClaim.cs)). The braid-carrying branch points are the **defective** exceptional points (Jordan blocks, eigenvectors coalesced, a loud √-branch), as opposed to the **diabolic** degeneracy (semisimple, eigenvectors independent, a silent double-discriminant-zero) at q ≈ 0.659 ([F89Path3OcticEpClaim](../compute/RCPsiSquared.Core/Symmetry/F89Path3OcticEpClaim.cs)).

The Liouvillian L(q) decomposes into (N+1)² joint-popcount sectors (p, q̃), where p is the bra (column) excitation count and q̃ the ket (row) count. The (1,2) octic lives in one such sector. **The question:** is the S₈ braid LOCALIZED to the (1,2) sector (and its symmetry images), or is it SHARED across the joint-popcount lattice; do other sectors carry the same defective branch points at the same q-loci?

The answer is **N-dependent, and that is the finding.**

## The instrument: the AT-aware sector probe

The census probes EVERY joint-popcount sector (dim > 1) at the (1,2) octic's reference DEFECTIVE loci ([ReferenceDefectiveLoci.For(N)](../compute/RCPsiSquared.Diagnostics/Foundation/ReferenceDefectiveLoci.cs), the P₁₀ simple zeros, conjugation-closed, the silent diabolic excluded), and classifies each coalescence as defective (a shared braid) or semisimple (a shared node) via the artifact-free `EpCharacter` reader.

The naive global-closest-pair probe **fails at N ≥ 5**: many raw joint-popcount blocks carry a PERMANENT semisimple degeneracy (a gap ~1e-15 eigenvalue pair at every q, sitting on the Absorption-Theorem dissipative rate lines). That permanent pair is the globally-closest pair, so the probe locks onto it, reads it diabolic, and MASKS the residual braid-carrying √-EP (gap ~5e-4) at a different λ. The fix is [`SectorEpProbe.ProbeDefectiveAnywhere`](../compute/RCPsiSquared.Diagnostics/Foundation/SectorEpProbe.cs): single-linkage-cluster every near-coalescence below a gap threshold, characterize each cluster on its own Riesz contour, and report `HasDefective` iff ANY cluster is a genuine Jordan block. The permanent semisimple cluster is skipped; the residual √-EP is reached. (The masking gotcha and the cluster fix are documented in [`SectorEpProbe`](../compute/RCPsiSquared.Diagnostics/Foundation/SectorEpProbe.cs)'s `ProbeDefectiveAnywhere` docstring.)

## The N=4 verdict: CONFINED to the D₄ orbit

At N=4 (the degenerate self-fold anchor, where the cross-fold partner w_{N−2} = w2 = DE coincides with the reference) the braid is shared across exactly the four sectors

> **{(1,2), (2,1), (2,3), (3,2)}**

These are the D₄ orbit of (1,2) under transpose (p ↔ q̃) and the bra/ket bit-flips — the canonical operator-space mirror group [D₄ = ⟨R, D⟩, Π = R·D](../docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md) (F118 `MirrorGroupD4Claim`), block-resolved onto the (bra-weight, ket-weight) lattice — all reading Defective at all reference loci with an identical coalescence gap (~3.7e-5). It is CONFINED there: it does **not** reach the dense half-filled core (2,2), which carries only a semisimple degeneracy. This echoes the Door-C finding that dissipative chaos is a filling threshold: the braid lives on the dilute orbit, not the dense core.

**Why N=4 is confined, the same self-fold that runs through all of F89.** N=4 is not incidentally special here: it is the F89 path-3 anchor, and its cross-fold q̃ ↦ N−2 = 2 maps (1,2) → (1,2) **onto itself** (the degenerate partner = self case, [F89CrossFoldSimilarityClaim](../compute/RCPsiSquared.Core/Symmetry/F89CrossFoldSimilarityClaim.cs); the within-block self-fold that put one diabolic on the real axis, and the same N=4→N=5 dichotomy [F89_PATH_K_DIABOLIC](F89_PATH_K_DIABOLIC.md) reads on the diabolics). Because the cross-fold is a self-map at N=4, the two families the spread splits into (below) COINCIDE (Family B folds back onto Family A), so there is no distinct family to spread into and no dense core to reach. The confinement is that degeneracy, read on the braid.

## The N=5 verdict: SPREADS to a 12-sector diamond

At N=5, the first generic length, the braid **reaches beyond the D₄ orbit** to a symmetric 12-sector diamond:

| p \ q̃ | 0 | 1 | 2 | 3 | 4 | 5 |
|--------|---|---|---|---|---|---|
| **1** | · | n | **B** | **B** | n | · |
| **2** | n | **B** | **B** | **B** | **B** | n |
| **3** | n | **B** | **B** | **B** | **B** | n |
| **4** | · | n | **B** | **B** | n | · |

**B** = braid (defective), **n** = node-only (a semisimple/AT near-coalescence, no braid), **·** = no coalescence. The 12 braid sectors are

> **(1,2) (1,3) (2,1) (2,2) (2,3) (2,4) (3,1) (3,2) (3,3) (3,4) (4,2) (4,3)**, including the dense diagonal cores (2,2) and (3,3).

All 12 read Defective at all 24 reference loci with a byte-identical gap signature (min 8.35e-5, max 6.12e-4). The set splits into two **cross-fold-conjugate families of 6**, each carrying a byte-identical shared eigenvalue λ (at any fixed locus all six sectors of a family agree on λ to 6 digits, the same branch point, not merely the same gap; the Re λ *range* in the table below is the variation ACROSS the 24 loci, not a spread within a family at fixed q):

| Family | Re λ range (over the 24 loci) | Sectors |
|--------|-------------------------------|---------|
| **A** | [−6.158591, −2.211959] | (1,2) (2,1) (2,3) (3,2) (3,4) (4,3) |
| **B** | [−7.788041, −3.841409] | (1,3) (2,2) (2,4) (3,1) (3,3) (4,2) |

and the two families are related by the F89d cross-fold **exactly**: `λ_B(q̄) = −conj(λ_A(q)) − 2N` (N=5). For example at q = 0.470875 − 0.047493i, (1,2) has λ = −4.028447 + 2.375482i, and −conj(λ) − 10 = −5.971553 + 2.375482i, which is exactly the (2,2) eigenvalue at q̄ ([F89CrossFoldSimilarityClaim](../compute/RCPsiSquared.Core/Symmetry/F89CrossFoldSimilarityClaim.cs), machine-zero for all N and q).

So the verdict is **N-dependent: confined at N=4, shared/coincident at N=5**, via SHARED SPECTRAL CONTENT, the same defective eigenvalue living in many joint-popcount sectors, a symmetry broader than the naive Klein-four.

## The mechanism: free-fermion / AT additivity

The membership is a clean rule:

> **braid set = {(p, q̃) : |p − q̃| = 1, both popcounts ∈ [1, N−1]} ∪ its q̃ ↦ N−q̃ cross-fold image.**

Family A is exactly the |bra − ket| = 1 sectors that keep at least one and at most N−1 excitations on each side; the |Δ| = 1 edge sectors that touch popcount 0 or N (an empty or full side) carry no coalescence: they fall in the braid-free "neither" set, its |Δ| = 1 members being (0,1),(1,0),(4,5),(5,4). Family B is the cross-fold image, not a separate phenomenon.

This single rule reproduces BOTH verdicts, and pins the N-dependence to the N=4 self-fold: at N=4 the cross-fold q̃ ↦ N−q̃ maps Family A **onto itself** (e.g. (1,2)→(1,2), (2,1)→(2,3)), so the rule collapses to the four-sector orbit and no dense core is reached: the confined verdict. At N ≥ 5 the self-fold lifts (the partner (1, N−2) ≠ (1,2)), Family B separates from Family A and brings in the dense diagonal cores, and the diamond opens: the spread verdict. **The N-dependence of the verdict is exactly the N-dependence of whether the cross-fold is a self-map on the |Δ| = 1 orbit, which it is only at N=4.**

Two operations, both verified live to machine precision, generate the whole diamond from the single elementary EP of the (1,2) SE-DE "rung":

1. **Diagonal mode-spectator (additivity).** Adding one excitation to BOTH the bra and the ket in the same single-particle mode leaves the energy difference E_bra − E_ket, and the entire EP (λ, character, gap), invariant. So (1,2) ≡ (2,3) ≡ (3,4) byte-identically: at q = 0.620878 the (1,2) and (2,3) defective λ agree to 1.2e-14. The single-particle energies are the free-fermion Bloch band E_n = 4·cos(πn/(N+1)) = {3.464, 2, 0, −2, −3.464} for N=5 (the XY chain is free-fermion integrable under Jordan-Wigner; the two-excitation energies are exact sums of one-excitation energies; see [XyJordanWignerModes](../compute/RCPsiSquared.Core/F86/JordanWigner/XyJordanWignerModes.cs), [DIABOLIC_BY_INTEGRABILITY](../hypotheses/DIABOLIC_BY_INTEGRABILITY.md)).

2. **The F89d cross-fold.** q̃ ↦ N−q̃ maps Family A to Family B with λ ↦ −λ̄ − 2N (an exact antiunitary similarity, [--root crossfold](../compute/RCPsiSquared.Diagnostics/Foundation/CrossFoldSimilarityWitness.cs)).

This is the free-fermion-kernel rule (**defective = single-particle coalescence; diabolic = frequency-difference coincidence**) operating across sectors: the branch point is a property of the |Δ| = 1 coherence rung, made popcount-translation-invariant by diagonal spectators and cross-fold-invariant by F89d, so the same defective eigenvalue lives in every diamond sector. This is the sector-resolved face of the [nonhermitian_diabolic_codimension](../compute/RCPsiSquared.Core/F86/LocalGlobalEpLink.cs) result (free-fermion additivity collapses the generically codim-3 coincidence to codim-1).

## The λ value, from below: the AT rate of the mixed eigenvector

The (1,2) block is the linear pencil **L(q) = A + q·C**, with **A = −2·diag(n_diff)** the REAL Absorption-Theorem rates (n_diff = popcount(a⊕b), so n_diff ∈ {1,3} for (1,2): overlap rate 2, no-overlap rate 6) and **C the ANTI-HERMITIAN free-fermion hopping** (the −i[H,·] part). Because C is anti-Hermitian, v^H·C·v is purely imaginary, so at real q the coherent term q·v^H·C·v contributes only to Im λ, and for any right-eigenvector v:

> **Re λ = v^H·A·v = −2·⟨n_diff⟩_v** (exact).

At the real defective loci this is machine-exact for the coalescing eigenvector: q = 0.620878 gives Re λ = −4.6189 = −2·⟨n_diff⟩ with ⟨n_diff⟩ = 2.3094 (residual 7e-15); q = 1.077615 gives −3.7917 with ⟨n_diff⟩ = 1.8958 (residual 1e-15). So the λ VALUE is NOT off the Absorption Theorem: it is the AT Theorem-2 rate of the coalescing eigenvector, which is a q-tuned MIXTURE of the overlap (rate 2) and no-overlap (rate 6) eigenmodes. ⟨n_diff⟩ is a weighted average between 1 and 3, so the rate lies between the pure −2 and −6 lines (off the integer-quantized lines when the mixture is not 50/50; at N=4 the defective EP sits at the −4 midpoint, ⟨n_diff⟩ = 2, yet is still Jordan). The defectiveness is the eigenvector coalescence (geo < alg, a genuine Jordan block), SEPARATE from the rate. At complex q the coherent term also feeds Re λ (Im q mixes it in), so the clean identity is a real-locus, physical-coupling statement. Live at `inspect --root sectorbraid` (the "λ VALUE from below" node).

## The closed-form mixture ⟨n_diff⟩(q): the free-fermion mode geometry

The (1,2) block is a free-fermion **contact problem**, and that fixes the mixture ⟨n_diff⟩ = 3 − 2⟨Ô⟩ completely. Rewrite the pencil as

> **L(q) = −6·I + 4·Ô + q·Ĥ**,

a **single particle on one side of the coherence and two on the other**, hopping as free fermions (Ĥ, the anti-Hermitian part), with a **contact attraction Ô** that counts coincidences: Ô = 1 when the lone excitation sits on one of the two paired sites, since n_diff = 3 − 2·Ô on this block. (This is the same pencil as L(q) = A + q·C from the section above, only regrouped: A = −2·diag(n_diff) = −6·I + 4·Ô because n_diff = 3 − 2·Ô, and Ĥ = C. The mode-product basis below writes the single-excitation side as the ket |k⟩ and the pair as the bra ⟨k₁,k₂|; which physical side that is is a harmless relabeling, since n_diff and Ô are symmetric under bra↔ket.) Without Ô every eigenvalue has Re = −6 (⟨n_diff⟩ = 3, the particles never coincide); the whole mixture is generated by the contact term.

Decompose the defective eigenvector into the free-fermion **mode-product basis** |k⟩⟨k₁,k₂| (the Jordan-Wigner sine modes ψ_m(x) = √(2/(N+1))·sin(πm x/(N+1)), single-particle energies ε_m = 2cos(πm/(N+1)), the [XyJordanWignerModes](../compute/RCPsiSquared.Core/F86/JordanWigner/XyJordanWignerModes.cs) dispersion; for N=5 this is {√3, 1, 0, −1, −√3}, half the block-scale band 4cos(·) of the mechanism section, which absorbs the J = 2q coupling convention; [JwBlockBasis](../compute/RCPsiSquared.Core/F86/JordanWigner/JwBlockBasis.cs)). In that basis the coincidence expectation of a pure mode-product is the overlap of the ket-mode density with the bra-mode density,

> **⟨Ô⟩_diag(k; k₁,k₂) = I(k,k₁) + I(k,k₂),   I(a,b) = Σ_l ψ_a(l)²·ψ_b(l)²**,

and this overlap is **EXACTLY QUANTIZED** (via the Dirichlet sum Σ_{m=1}^N cos(2π c m/(N+1)) = −1 unless (N+1)|c):

> **I(a,b) = 1/(N+1)** generic;  **I(a,b) = 3/(2(N+1))** if a = b or chiral a + b = N+1.

So the diagonal part of ⟨Ô⟩ is always rational. The defective eigenvector, written in this basis as v = Σ_α c_α·|k⟩⟨k₁,k₂| (the c_α its mode-product amplitudes), is a **mixture within a δ-degenerate multiplet** (δ = ε_k − ε_{k₁} − ε_{k₂}; the contact term is off-diagonal in the mode basis and couples equal-δ triples, chirally ±-symmetric because the defective λ is real), so

> **⟨Ô⟩ = Σ_α |c_α|²·⟨Ô⟩_diag(α)  +  (off-diagonal δ-multiplet mixing).**

**The off-diagonal mixing is the closed-form boundary.** At **N=4 it is exactly zero** at every locus: the self-fold pairs the I = 3/(2(N+1)) chiral partners with equal weight and the cross-terms cancel, so ⟨Ô⟩ = 1/2 and ⟨n_diff⟩ = N/2 = 2 exactly. This is a **second, independent proof** of the N=4 pinning (the first is the eigenvalue-side Re λ = −N self-fold fixed axis above): the mode geometry shows the pinning as a *vanishing of the mixing*, not merely a fixed point. At **N ≥ 5** there is no self-fold, the mixing is nonzero (N=5: the δ=±1 multiplet at locus 1, mixing −0.105; the δ=±(√3−1) multiplet at locus 2, +0.190), and it shifts ⟨Ô⟩ off the rational, so the per-locus value is the non-radical real part of an S₈ octic root. Live at `inspect --root sectorbraid` (the mixture-table and mode-geometry nodes); the decomposition is [SectorBraidModeGeometry](../compute/RCPsiSquared.Diagnostics/Foundation/SectorBraidModeGeometry.cs).

**The function ⟨n_diff⟩(q) itself.** Tracking the defective branch as a smooth function of real q ([SectorBraidModeGeometry.SweepDefectiveBranch](../compute/RCPsiSquared.Diagnostics/Foundation/SectorBraidModeGeometry.cs); the AT identity Re λ = −2⟨n_diff⟩ is now gate-verified along the *whole* branch, real and complex λ, residual < 1e-9) shows ⟨n_diff⟩(q) is a smooth algebraic function with a **√-branch point exactly at each exceptional point**. For N=5 the branch is real below the exceptional point, ⟨n_diff⟩ falling from the q=0 quantized endpoint value 3 (low-q as 3 − c·q²) down through ~2.35, then steepening to the branch-point value 2.31 as q → 0.622 and the slope diverges; at q ≈ 0.622 it hits the √-branch point and continues as a complex-conjugate pair. There is no rational or radical global form: **the √-branch points at the S₈ exceptional points ARE the analytic obstruction.** The S₈ "wall" is not only the Galois group of the loci; it is the branch-point structure of the function ⟨n_diff⟩(q) itself, and a radical per-locus form exists exactly when the self-fold acts (N=4).

So the closed-form question is resolved: the mixture is the amplitude-weighted contraction ⟨Ô⟩ = v†Ô v / v†v (the right eigenvector v, Hermitian inner product, as in the λ-value section) of the quantized overlaps I(a,b) over the δ-multiplet eigenvector; it closes at N=4 (the self-fold zeroes the mixing → ⟨n_diff⟩ = N/2), and it is non-radical for N ≥ 5 because the function ⟨n_diff⟩(q) carries √-branch points at the S₈ exceptional points.

## What is open

The mixture ⟨n_diff⟩(q) is now understood (the quantized-overlap contraction above; closed at N=4, non-radical at N ≥ 5 by the √-branch-point structure). What remains: the **general-N membership** beyond N = 4, 5 (does the diamond grow with N as additivity predicts, verified at N=6?), and the **codim-1-by-additivity theorem** (why free-fermion additivity collapses the generically codim-3 coincidence to codim-1) as a stated Tier-1 theorem rather than a numerically-confirmed rule.

## Reproduce

```bash
# The live, browsable verdict (default N=4, the fast anchor):
dotnet run --project compute/RCPsiSquared.Cli -c Release -- inspect --root sectorbraid
# The generic spread verdict (N=5, ~2 min: runs the full census live):
dotnet run --project compute/RCPsiSquared.Cli -c Release -- inspect --root sectorbraid --N 5

# The gate-tested verdicts (the exact 12-set pinned; SLOW_MSM):
dotnet test compute/RCPsiSquared.Diagnostics.Tests -c Release \
  --filter "FullyQualifiedName~MultiSectorMonodromyCensusTests"

# The AT identity along the whole defective branch (the function thread; SLOW_MSM):
dotnet test compute/RCPsiSquared.Diagnostics.Tests -c Release \
  --filter "FullyQualifiedName~SectorBraidBranchSweep"
```

The census + probe + mode geometry live in [`compute/RCPsiSquared.Diagnostics/Foundation/`](../compute/RCPsiSquared.Diagnostics/Foundation/) (`MultiSectorMonodromyCensus`, `SectorEpProbe`, `ReferenceDefectiveLoci`, `SectorBraidWitness`, `SectorBraidModeGeometry`); the exact N=5 12-set is the regression test `Census_N5_BraidSpreadsBeyondOrbit_ReachesDenseCore`, and the AT-identity-along-the-branch gate is `Sweep_N5_DefectiveBranch_ATIdentity_AndDump`. The `inspect --root sectorbraid` witness now carries the ⟨n_diff⟩ mixture table (all real loci) and the free-fermion mode-geometry decomposition of the defective eigenvector.
