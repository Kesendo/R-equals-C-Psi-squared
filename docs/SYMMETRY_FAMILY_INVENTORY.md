# Symmetry Family Inventory: chain XY + Z-dephasing Liouvillian

**Purpose:** Systematic typed inventory of all discrete symmetries of the canonical chain XY + per-site Z-dephasing Liouvillian L = −i[H, ·] + Σ_l γ_l (Z_l ρ Z_l − ρ). Each symmetry is classified by axis (operator vs parameter) and effect (sector-splitting vs sector-pairing). New typed Claims live under `compute/RCPsiSquared.Core/SymmetryFamily/`.

**Scope of "discrete":** Z_n-style finite groups acting on the Liouville space or on the parameter coefficients. Excludes continuous U(1) (already covered by `JointPopcountSectors`).

## The five known symmetry-axes

| # | Symmetry | Axis | Effect | Status | Anchor |
|---|---|---|---|---|---|
| 1 | **U(1) × U(1)** per-side popcount | operator (continuous) | block-diagonal (N+1)² sectors | Tier1Derived | `JointPopcountSectors` |
| 2 | **F71** chain spatial mirror (Z₂) | operator | within-sector splitting (factor 2), chain only, palindromic γ only | Tier1Derived | `F71MirrorBlockRefinement` |
| 3 | **F91** γ-Z₄ anti-palindromic | parameter (γ_l) | spectral invariance on diagonal-block layer | Tier1Derived | `F71AntiPalindromicGammaSpectralInvariance` |
| 4 | **F92** J-Z₄ anti-palindromic (NEW) | parameter (J_b) | spectral invariance on diagonal-block layer | Tier1Derived | `F92BondAntiPalindromicJSpectralInvariance` |
| 5 | **F93** h-Z₄ anti-palindromic (NEW) | parameter (h_l detuning) | spectral invariance on diagonal-block layer | Tier1Derived | `F93DetuningAntiPalindromicSpectralInvariance` |

## Sector-pairing primitives (do not split; identify equal-spectrum sectors)

| # | Symmetry | Axis | Effect | Status | Anchor |
|---|---|---|---|---|---|
| 6 | **Z⊗N** Pauli-letter parity (Z₂) | operator | trivially redundant with joint-popcount-parity (typed for inventory completeness) | Tier1Derived | `ZGlobalMirrorRefinement` |
| 7 | **X⊗N** charge-conjugation (Z₂) | operator | sector-pairing (p_c, p_r) ↔ (N−p_c, N−p_r); halves number of eig-calls | Tier1Derived | `XGlobalChargeConjugationPairing` |

## Negative results (typed for completeness)

| # | Symmetry | Status | Anchor |
|---|---|---|---|
| 8 | **F71_col × F71_row** factor-4 split | Tier2Empirical (does NOT hold: Z⊗Z dissipator correlates both sides; only diagonal product survives) | `F71BilateralBlockRefinement` |

## Symmetries from other axes (typed elsewhere; not in SymmetryFamily but relevant)

- **K** (chiral / sublattice / AZ class BDI, KHK = −H): Hamilton-symmetry only; NOT a full L-symmetry under Z-dephasing → not a BlockSpectrum refinement. Typed as `ChiralK.cs` in `Symmetry/`.
- **Π** (palindrome master, F1): non-diagonal master operator; F1 spectrum mirror around −Σγ. Typed as `PiOperator.cs` in `Symmetry/` and as `F1` in `ANALYTICAL_FORMULAS.md`.
- **Pi2-Z₄** (operator-quaternion, NinetyDegreeMirrorMemory): operator-side of the same Z₄ that F91/F92/F93 inhabit on the parameter side. Typed in `Pi2KnowledgeBaseClaims.cs`.

## Combinatorial reach

For chain XY + uniform Z-dephasing, the simultaneous use of:
- joint-popcount: (N+1)² sectors
- F71 split: factor ~2
- X⊗N pairing: halves number of distinct eig-calls (paired sectors share spectrum)

gives at N=10 a max-block of `C(10,5)² / 2 / 2 ≈ 16k`, still beyond commodity hardware. F92/F93 do NOT add block-reduction (they describe spectral-invariance under parameter rotation, not sector splitting). Therefore N=10 push requires a **new** symmetry not yet found, or a different attack (matrix-free per-block eig, etc.) — the latter now realised, see the next section.

## N=10 push: the realised "different attack" (2026-05-16)

Phase 2 of the N=10 push landed five primitives under `compute/RCPsiSquared.Core/BlockSpectrum/JordanWigner/`:

| Primitive | Role |
|---|---|
| `JwSlaterPairBasis` | Basis transform U: maps the (p_c, p_r) computational-basis sector to the JW Slater-pair basis; L_H diagonal with eigenvalues `−i·(Σε(L) − Σε(K))`. |
| `JwSlaterPairLProjection` | Dense `U^†·L·U` + Slater-swap sparsity witness `nnz_off/row ≤ (1 + p_r·(N−p_r))·(1 + p_c·(N−p_c)) − 1`; density falls from ~37 % at N=4 to ~10 % at N=7, extrapolated ~1 % at N=10. Verifies the bound where dense `U^†LU` still fits. |
| `JwSlaterPairSparseLBuilder` | Direct sparse CSR construction of L_JW element-by-element from the analytic `Z_l ⊗ Z_l` action on Slater pairs — no dense `U^†LU` detour. At N=10 (p_c=5, p_r=5, sectorDim 63 504): ~14 M nnz, ~280 MB, ~5 s build. **This is the N=10 enabler.** |
| `JwSlaterPairArnoldiEig` | Managed non-symmetric Arnoldi with Modified Gram-Schmidt on the CSR matvec; returns top-k largest-magnitude Ritz values. Matvec parallelised via `Parallel.For` over rows. |
| `JwSlaterPairShiftInvertArnoldi` | Same outer Arnoldi on `(L − σI)^(−1)`; inner solver is an in-house BiCGStab on the same CSR matvec — no MathNet sparse conversion. Recovers eigenvalues nearest σ, including the steady state at λ = 0 (machine precision at N=10). |

Shared low-level complex-vector kernels (`NormSquared`, `ConjugateDot`, `AxpyInPlace`, `RandomNormalized`) live in `KrylovOps.cs` to keep the two Arnoldi sites lean.

The combinatorial table above is unchanged — joint-popcount + F71 + X⊗N still do not shrink the max block enough for dense LAPACK at N=10 — but the matrix-free Krylov path renders that limit moot for top-k / near-σ work.

**What this path does NOT yet give:**
- **Full N=10 spectrum** (all ~63 k eigenvalues). For F1 palindrome verification at N=10 the structurally cleaner route is Prosen third quantization (Phase 3 of the N=10 plan): under uniform γ, the chain XY + Z-dephasing Liouvillian reduces to a 2N × 2N Nambu spectral problem (Medvedyeva, Essler, Prosen 2016). Phase 3 has not yet been implemented.
- **Tight shift-invert convergence** at N=10. The in-house BiCGStab stagnates at ~866 inner iter per outer step (budget cap 1000) on the highly non-normal Liouvillian; the steady state and slow modes are still recovered cleanly, but preconditioning (Jacobi or ILU0) would cut wall time substantially.
- **Inhomogeneous γ spectrum** at N=10 is supported by the sparse path but not yet exercised in a science-facing pipeline; Prosen's reduction does not extend there.

## Cross-references

- `JointPopcountSectors`, `F71MirrorBlockRefinement`, `F71BilateralBlockRefinement`, `F71AntiPalindromicGammaSpectralInvariance`: see `compute/RCPsiSquared.Core/BlockSpectrum/`.
- `F1`, `F71`, `F91`, F92, F93: see `docs/ANALYTICAL_FORMULAS.md`.
- F91/F92/F93 algebraic proofs: `docs/proofs/PROOF_F91_GAMMA_NINETY_DEGREES.md`, `docs/proofs/PROOF_F92_BOND_ANTI_PALINDROMIC_J.md`, `docs/proofs/PROOF_F93_DETUNING_ANTI_PALINDROMIC.md`.
- N=10 push primitives: `compute/RCPsiSquared.Core/BlockSpectrum/JordanWigner/JwSlaterPairBasis.cs`, `JwSlaterPairLProjection.cs`, `JwSlaterPairSparseLBuilder.cs`, `JwSlaterPairArnoldiEig.cs`, `JwSlaterPairShiftInvertArnoldi.cs`, `KrylovOps.cs`.
- Synthesis: `reflections/ON_THE_SYMMETRY_FAMILY.md`, `reflections/ON_THE_NINETY_DEGREE_GAMMA.md`.
