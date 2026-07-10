# Symmetry Family Inventory: chain XY + Z-dephasing Liouvillian

**Purpose:** Systematic typed inventory of all discrete symmetries of the canonical chain XY + per-site Z-dephasing Liouvillian L = −i[H, ·] + Σ_l γ_l (Z_l ρ Z_l − ρ). Each symmetry is classified by axis (operator vs parameter) and effect (sector-splitting vs sector-pairing). New typed Claims live under `compute/RCPsiSquared.Core/SymmetryFamily/`.

**Scope of "discrete":** Z_n-style finite groups acting on the Liouville space or on the parameter coefficients. Excludes continuous U(1) (already covered by `JointPopcountSectors`).

## The five known symmetry-axes

| # | Symmetry | Axis | Effect | Status | Anchor |
|---|---|---|---|---|---|
| 1 | **U(1) × U(1)** per-side popcount | operator (continuous) | block-diagonal (N+1)² sectors | Tier1Derived | `JointPopcountSectors` |
| 2 | **F71** chain spatial mirror (Z₂) | operator | within-sector splitting (factor 2), chain only, palindromic γ only | Tier1Derived | `F71MirrorBlockRefinement` |
| 3 | **F91** γ anti-palindromic (parameter-side Klein V₄, shadow of the operator-side Z₄) | parameter (γ_l) | spectral invariance on diagonal-block layer | Tier1Derived | `F71AntiPalindromicGammaSpectralInvariance` |
| 4 | **F92** J anti-palindromic (parameter-side Klein V₄, shadow of the operator-side Z₄) (NEW) | parameter (J_b) | spectral invariance on diagonal-block layer | Tier1Derived | `F92BondAntiPalindromicJSpectralInvariance` |
| 5 | **F93** h anti-palindromic (parameter-side Klein V₄, shadow of the operator-side Z₄) (NEW) | parameter (h_l detuning) | spectral invariance on diagonal-block layer | Tier1Derived | `F93DetuningAntiPalindromicSpectralInvariance` |

## Sector-pairing primitives (do not split; identify equal-spectrum sectors)

| # | Symmetry | Axis | Effect | Status | Anchor |
|---|---|---|---|---|---|
| 6 | **Z⊗N** Pauli-letter parity (Z₂) = Π²_X (F61) | operator | trivially redundant with joint-popcount-parity (typed for inventory completeness) | Tier1Derived | `ZGlobalMirrorRefinement` |
| 7 | **X⊗N** charge-conjugation (Z₂) = Π² (F1²) | operator | sector-pairing (p_c, p_r) ↔ (N−p_c, N−p_r); halves number of eig-calls | Tier1Derived | `XGlobalChargeConjugationPairing` |

## Architecture aggregators and refinement Claims (inventory state; the finer YParity / BitA / BitB / Cubic3 refinement axes beyond the top-level five)

| # | Symmetry | Axis | Effect | Status | Anchor |
|---|---|---|---|---|---|
| 8 | PolarityCubeMap | aggregator | Z2Axis classification of all Pi²-Inheritance Claims | gap-detection via null BitATwin (Tier2Empirical, snapshot strength bounded by weakest contributing parent; cubic architecture itself is Tier1Derived) | `PolarityCubeMap` (added 2026-05-24) |
| 10 | YParityIndependenceAtK3 | YParity (term-level Z₂) | classifies Pauli terms by n_Y mod 2 | independent of Klein at k_body≥3 (Tier1Derived; F102) | `YParityIndependenceAtK3` (added 2026-05-24) |
| 11 | F87Z2CubedRefinementN4K3 | YParity (term-level Z₂) | refines F87 trichotomy at k=3 with five sub-statements (truly y_par=0-pure, hard diagonal 42:8 with Y-inversion, diagonal soft 13:13, mother soft y_par=1-pure, off-diagonal soft Pattern B + C) | empirical anchor at N=4 (294 pairs; Tier1Derived; F103). First derived class of `F87Z2CubedRefinementBase` | `F87Z2CubedRefinementN4K3` (added 2026-05-24) |
| 12 | F87Z2CubedRefinementN5K3 | YParity (term-level Z₂) | refines F87 trichotomy at N=5 k=3 with the same five sub-statements as F103; same 294-pair enumeration as F103 (N4K3); F85 N-stability lift test on the y_par axis | empirical anchor at N=5 (294 pairs; Tier1Derived; F105). F85 N-stability lift to y_par axis CONFIRMED bit-exactly: all 5 sub-statement counts at N=5 k=3 match F103's at N=4 k=3 identically. Second derived class of `F87Z2CubedRefinementBase` | `F87Z2CubedRefinementN5K3` (added 2026-05-24) |
| 13 | F87Z2CubedRefinementN4K4 | YParity (term-level Z₂) | refines F87 trichotomy at N=4 k=4; new 4248-pair enumeration; tests k-stability where F105 confirmed N-stability | empirical anchor at N=4 k=4 (4248 pairs; Tier1Derived; F106). Mixed k-stability: two purity statements HELD (truly y_par=0-pure 3924/0; mother soft y_par=1-pure (0, 300) ×3), Y-inversion HELD qualitatively; the k=3 42:8 hard ratio sharpened to 100% pure (228, 0), the 13:13 diagonal soft symmetry broke to asymmetric (300, 528). Off-diagonal: 3 cells stay Pattern C analog (0, 528), 3 cells became fully symmetric (528, 528). Third derived class of `F87Z2CubedRefinementBase` | `F87Z2CubedRefinementN4K4` (added 2026-05-24) |
| 14 | TrulyYParityZeroPurity | YParity (term-level Z₂) | truly classification forces y_par = 0 on every Pauli term under any single-letter dephase channel; closed-form corollary of F85's k-body truly criterion extended to X and Y dephase via per-dephase Π² eigenvalue + dissipator commutativity | Tier1Derived; F107. First DERIVED-not-EMPIRICAL Claim in the F87 Z₂³ refinement family; explains the y_par=1 truly count of zero observed empirically across F103/F105/F106 (4524 truly classifications total) | `TrulyYParityZeroPurity` (added 2026-05-24) |
| 15 | KleinEightCellClaim | Cubic3 (full Z₂³) | F88a + F102 8-cell Z₂³ decomposition (bit_a, bit_b, y_par); structural anchor for F87 Z₂³ refinement work | Tier1Derived; Stage 2b. First Cubic3-axis Claim. At fixed even k_body collapses to F88a's Klein 4-cell (y_par = bit_a XOR bit_b); at fixed odd k_body collapses to the opposite 4. Mixed-k_body enumerations populate all 8 | `KleinEightCellClaim` (added 2026-05-24) |
| 16 | MotherSoftYParityOnePurity | YParity (term-level Z₂) | mother sector Klein (0,0) soft pairs under any single-letter dephase channel are y_par=1 pure; all three dephase branches closed-form via the F108 Part 1+2+3 Π_5bilinear family | Tier1Derived (fully unconditional after F108 Part 1+2+3 closure 2026-05-25); F109. Sister to F107 on the y_par axis; together they pin the y_par signature of truly + mother-soft cells. Explains the y_par=0 mother-soft count of zero observed empirically across F103/F105/F106 (1026 mother-soft classifications) | `MotherSoftYParityOnePurity` (added 2026-05-24, fully unconditional 2026-05-25 after F108 Part 1+2+3 closure) |
| 17 | F108Part1Pi2EvenAlwaysPalindromic | BitB (Π²_Z = X⊗N) | every Π²_Z-even Hamiltonian H built from bilinears {XX, YY, YZ, ZY, ZZ} + Z-dephasing admits EXACT operator-level palindrome via the Π_5bilinear phase-variant operator (I→+X, X→−I, Y→+iZ, Z→−iY). Closes the F87-hardness question for the Π²-even half of Klein cells under Z-dephasing | Tier1Derived; F108 Part 1, closed 2026-05-25. First BitB-axis Claim of the post-F107 wave. BitATwin = Filled (typed ctor parent = F108 Part 2); Y-dephasing analog (F108 Part 3) has no covering Claim yet. Retroactively closes F109's Step 5 for Z-dephasing | `F108Part1Pi2EvenAlwaysPalindromic` + `Pi5BilinearOperator` (added 2026-05-25) |
| 18 | F108Part2Pi2XEvenAlwaysPalindromic | BitA (Π²_X = Z⊗N) | every Π²_X-even Hamiltonian H built from bilinears {ZZ, XX, XY, YX, YY} + X-dephasing admits EXACT operator-level palindrome via the X-dephasing variant of Π_5bilinear (I→+Z, Z→−I, X→−iY, Y→+iX). Closes the F87-hardness question for the Π²_X-even half of Klein cells under X-dephasing | Tier1Derived; F108 Part 2, closed 2026-05-25. BitA twin of F108 Part 1; second BitA-axis Claim (after F61). Pi5BilinearOperator extended with dephase parameter (Z, X, Y all supported). Closes F109's Step 5 X-dephasing branch closed-form | `F108Part2Pi2XEvenAlwaysPalindromic` + `Pi5BilinearOperator` (X-deph branch, added 2026-05-25) |
| 19 | F108Part3Pi2YEvenAlwaysPalindromic | BitB (Π²_Y = bit_b parity, same axis as Π²_Z) | every Π²_Y-even Hamiltonian H built from bilinears {XX, YY, YZ, ZY, ZZ} (same set as Π²_Z-even) + Y-dephasing admits EXACT operator-level palindrome via the Y-dephasing variant of Π_5bilinear (I→+X, X→−I, Y→−iZ, Z→+iY; only Y/Z 2-cycle phase differs from Part 1's +i to Y-deph's −i). Closes the F87-hardness question for Π²_Y-even cells under Y-dephasing | Tier1Derived; F108 Part 3, closed 2026-05-25. Y-dephasing sibling of F108 Part 1 (same BitB axis, same bilinear set, different dephase letter); BitBSpecific BitATwin slot (Y-deph has no meaningful bit_a analog). Reuses F108 Part 1's bilinear predicates. Closes F109's Step 5 Y-dephasing branch closed-form; F109 now fully unconditional across {Z, X, Y} | `F108Part3Pi2YEvenAlwaysPalindromic` + `Pi5BilinearOperator` (Y-deph branch, added 2026-05-25) |
| 20 | HardCellYInversionPattern | YParity (term-level Z₂) | F87-hard pairs only in diagonal Klein cells (Aspect A closed-form via F108 Part 1+2+3 + F107 + F109 + F87 dissipator-resonance) with Y-inversion structural pattern (dominant y_par in hard cell equals dephase letter's own y_par; Z/X→0, Y→1) and k-purity sharpening (42:8 @ k=3 → 228:0 @ k=4 with Y-inversion preserved) | Tier1Derived (promoted 2026-06-10); F110, added 2026-05-25. Seventh YParity-axis Claim; completes the y_par-axis classification of the F87 trichotomy together with F107 (truly) and F109 (mother soft). Aspect A closed-form, Aspect B+C derived via F103 §6/§7 (anchored at F103/F105/F106); the windowed hard-direction converse closed 2026-06-10 (WindowedConverseAllGammaClaim, no residual). F111 closes Aspect B at k=N=4 via the Pure-D Template Rule (Tier1Derived) | `HardCellYInversionPattern` (added 2026-05-25) |
| 21 | HardCellPureDTemplate | YParity (term-level Z₂) | at k=N=4 in diagonal Klein cell (D.BitA(), D.BitB()) for dephase D, pair is F87-hard iff at least one term is a "pure-D template" (length-4 string with only D and I letters). Per-cell decomposition: 36 Pure-Pure (HARD) + 192 Pure-Mixed (HARD) + 300 Mixed-Mixed (SOFT) = 528 pairs, 228 hard. Implies F110 Aspect B at k=4 as immediate corollary (pure-D templates have y_par=y_par(D) by construction) | Tier1Derived (promoted 2026-06-10); F111, added 2026-05-25. Eighth YParity-axis Claim, sharpens F110 Aspect B. Subclaim (a) pure-D=hard heuristically derived via dissipator-commute mechanism; subclaim (d) Mixed+Mixed=soft closed modulo M via PROOF_F103 §7.4 (2026-05-30); the hard-direction converse behind (a)/(c) closed 2026-06-10 (WindowedConverseAllGammaClaim, Pascal-Gram positivity F117, no residual), promoting to Tier1Derived. Empirical anchor: 1584 pair classifications across 3 dephases match the rule, zero exceptions | `HardCellPureDTemplate` (added 2026-05-25) |
| 22 | LindbladBitBPiBalance | BitB (Π² = (-1)^bit_b) | for any Lindblad-form L = -i[H, ·] + Σ_k γ_k · np.kron(c_k, c_k^*) with Hermitian H and each c_k bit_b-homogeneous (every Pauli string σ in c_k's expansion shares bit_b(σ) = (#Y + #Z) mod 2 = const), the `polarity_coordinates_from_L` decomposition of M = Π L Π⁻¹ + L + 2σ·I satisfies ‖M_plus_half‖² = ‖M_minus_half‖² bit-exactly. Structural identity behind the polarity_coordinates diagnostic; asymmetry ≠ 0 detects c with cross-bit_b Pauli support (outside the F108 closure regime) | Tier1Derived; F112, added 2026-05-26. Hermitian H is the typed scope (rigorous 5-step proof via F38 / F63 Π² eigenvalue formula + dagger / anti-Hermitian L_H argument); non-Hermitian H extension empirical only (20 random configs N=2, 3 bit-exact; Tier1Candidate documented in inspectables, not typed as a separate Claim). F87 ↔ F112 orthogonal axes on shared bit_b Z₂-grading: F87 lives in spec(L) palindromy, F112 lives in M_anti's Π +i / −i split. Ctor parent: F108 Part 1 (shared bit_b foundation). BitBSpecific BitATwin slot (intrinsically bit_b-axis theorem) | `LindbladBitBPiBalance` (added 2026-05-26) |
| 23 | LindbladBitBPiBreakMagnitude | BitB (Π² = (−1)^bit_b) | closed-form magnitude for the F112 polarity-asymmetry counterexample when F112's typed scope is violated by the canonical Z-drive × amplitude-damping interference: for Lindblad-form L with Hermitian H containing single-site Z-drives Σ_l (ω_l/2)·Z_l plus any bit_b-homogeneous additions (each F112-balanced) and dissipator c containing σ⁻_l rate γ_T1,l + σ⁺_l rate γ_pump,l per site (standard physics convention: σ⁻ lowering), asymmetry = (4^N / 2) · Σ_l ω_l · (γ_pump,l − γ_T1,l) bit-exactly. T1 cooling at positive ω gives negative asymmetry; pumping flips to positive; detailed balance cancels. Structural origin in [Z, σ⁻] = −2·σ⁻ (proportional to the non-Hermitian σ⁻ itself, producing Π-eigenspace ±i imbalance); X- and Y-drives give Hermitian commutators that remain F112-symmetric, hence contribute 0. Per-site additivity from same-site locality [Z_l, σ⁻_m] = −2·σ⁻_m · δ_{lm}. Helpers: `PredictAsymmetry(ω_l, γ_T1,l, γ_pump,l, N)`, `PredictAsymmetryUniform(ω, γ_T1, γ_pump, N)` | Tier1Derived at N=2, 3, 4 (bit-exact via constructive parameter sweep + per-site / cross-site / sign-flip / detailed-balance / non-uniform-rate verification, `simulations/f113_break_formula_derivation.py`); F113, added 2026-05-26. Tier1Derived for general N (2026-05-26): the (1/2)·4^N coefficient's rigorous algebraic derivation from Π-eigenspace structure is in `PROOF_F113_COEFFICIENT_DERIVATION.md`. Sister to F112 on the same bit_b axis: F112 covers in-scope balance (asymmetry = 0), F113 the out-of-scope counterexample magnitude; together complete polarity-axis description of the standard Lindblad family. Ctor parent: F112 (LindbladBitBPiBalance). BitBSpecific BitATwin slot (Z-axis single-site drives intrinsically bit_b, no meaningful bit_a analog). Hardware fingerprinting: Welle 2 f95 fit (ω=0.13, γ_T1≈0.001, N=2) gives F113-predicted −2.08e-3, matching Kingston bit-exact (experiments/F112_HARDWARE_LENS_KINGSTON.md) | `LindbladBitBPiBreakMagnitude` (added 2026-05-26) |
| 24 | F38BitAInvolutionInheritance | BitA (Π²_X = (−1)^{n_XY}) | Π²_X involution eigenvalue formula on Pauli strings, with the same half-half (2^{2N-1}/2^{2N-1}) eigenspace partition as F38's Π²_Z. Structural BitA twin of F38: by Hadamard conjugation Π_X = H^⊗N·Π_Z·H^⊗N, the n_XY parity counter mirrors the w_YZ counter exactly. Foundational for the bit_a-axis structure | Tier1Derived; Welle 7, added 2026-05-26. First TRIVIAL Z↔X mirror typed in the post-F108 era; fills the BitA twin slot of F38 (the most-cited Π² involution Claim) | `F38BitAInvolutionInheritance` (added 2026-05-26) |
| 25 | F39DetPiBitAInheritance | BitA (Π²_X = (−1)^{n_XY}) | det(Π_X) = (−1)^{N·4^{N−1}}: identical to det(Π_Z) via Hadamard conjugation invariance det(Π_X) = det(H^⊗N · Π_Z · H^⊗N) = det(Π_Z). Both equal +1 for N ≥ 2, −1 for N = 1 | Tier1Derived; Welle 7, added 2026-05-26. Second TRIVIAL Z↔X mirror; fills F39 BitA twin slot. Algebraic identity, no physics asymmetry | `F39DetPiBitAInheritance` (added 2026-05-26) |
| 26 | F63BitAReference | BitA (Π²_X = (−1)^{n_XY}) | lightweight bit_a-axis Claim establishing [L, Π²_X] = 0 as the structural twin of F63's [L, Π²_Z] = 0. The actual proof lives in F61BitAParityPi2Inheritance (which proves the bit_a-axis [L, Π²] commutation directly); this Claim provides the typed BitA twin slot for F63 without forming a F61 → F63 → F38 → F61 ctor cycle. F61 referenced in docstring + inspectables, not as ctor parent | Tier1Derived; Welle 7, added 2026-05-26. Lightweight cycle-breaking reference Claim; the cleanest way to expose F61's bit_a-axis content as F63's typed BitA twin | `F63BitAReference` (added 2026-05-26) |
| 27 | ZGlobalEigenstateMirrorBitAInheritance | BitA (Π²_X = (−1)^{n_XY}) | Z⊗N\|ψ⟩ = ±\|ψ⟩ ⟹ γ_X = ±1 ⟹ α_total = 0 on the X-Frobenius axis. Computational basis states \|0...0⟩, \|1...1⟩ are the canonical Z⊗N eigenstates with γ_X = ±1, mirror of \|+...+⟩, \|−...−⟩ (X⊗N eigenstates with γ_Z = ±1) from `XGlobalEigenstateMirrorPi2Inheritance` | Tier1Derived; Welle 7, added 2026-05-26. F99 0°-anchor closed on both Π² axes; fills the BitA twin slot of XGlobalEigenstateMirror via Hadamard rotation | `ZGlobalEigenstateMirrorBitAInheritance` (added 2026-05-26) |

## Negative results (typed for completeness)

| # | Symmetry | Status | Anchor |
|---|---|---|---|
| 9 | **F71_col × F71_row** factor-4 split | Tier2Empirical (does NOT hold: Z⊗Z dissipator correlates both sides; only diagonal product survives) | `F71BilateralBlockRefinement` |

## Symmetries from other axes (typed elsewhere; not in SymmetryFamily but relevant)

- **K** (chiral / sublattice / AZ class BDI, KHK = −H): Hamilton-symmetry only; NOT a full L-symmetry under Z-dephasing → not a BlockSpectrum refinement. Typed as `ChiralK.cs` in `Symmetry/`.
- **Π** (palindrome master, F1): non-diagonal master operator; F1 spectrum mirror around −Σγ. Typed as `PiOperator.cs` in `Symmetry/` and as `F1` in `ANALYTICAL_FORMULAS.md`. Its square Π² is the X⊗N charge-conjugation (#7 above, registered as F1²).
- **Pi2-Z₄** (operator-quaternion, NinetyDegreeMirrorMemory): the genuine order-4 Z₄ (`i⁴=1`, `i²=−1`) on the operator side. F91/F92/F93 inhabit its order-2 *shadow* on the parameter side, a Klein V₄ = Z₂×Z₂ (palindromic mirror F71 × anti-palindromic involution R₉₀), not the order-4 Z₄ itself. Typed in `Pi2KnowledgeBaseClaims.cs`. Its from-below companion at the real exceptional point: the eigenVECTOR-frame holonomy around the living N=9 defective seed rotates the coalescing-eigenvector frame 90° per loop (biorthogonal vᵀv gauge; frame-monodromy generator with eigenvalues ±i, so M₂ = −I, M₄ = +I, single-valued only after four loops = `i⁴=1`), the eigenVECTOR-phase partner of the eigenVALUE-swap `Monodromy.cs`. The mod-4 is the vᵀv-gauge reading; the gauge-invariant kernel is the order-2 eigenvalue swap. A noted correspondence to this algebraic Z₄, not a derived identity. Typed as `SeedHolonomyClaim` (Tier1Candidate); live via `inspect --root holonomy` (`SeedHolonomyWitness`).

## Combinatorial reach

For chain XY + uniform Z-dephasing, the simultaneous use of:
- joint-popcount: (N+1)² sectors
- F71 split: factor ~2
- X⊗N pairing: halves number of distinct eig-calls (paired sectors share spectrum)

gives at N=10 a max-block of `C(10,5)² / 2 / 2 ≈ 16k`, still beyond commodity hardware. F92/F93 do NOT add block-reduction (they describe spectral-invariance under parameter rotation, not sector splitting). Therefore N=10 push requires a **new** symmetry not yet found, or a different attack (matrix-free per-block eig, etc.); the latter now realised, see the next section.

## N=10 push: the realised "different attack" (2026-05-16 / 2026-05-17)

Phase 2 of the N=10 push landed six primitives under `compute/RCPsiSquared.Core/BlockSpectrum/JordanWigner/`:

| Primitive | Role |
|---|---|
| `JwSlaterPairBasis` | Basis transform U: maps the (p_c, p_r) computational-basis sector to the JW Slater-pair basis; L_H diagonal with eigenvalues `−i·(Σε(L) − Σε(K))`. |
| `JwSlaterPairLProjection` | Dense `U^†·L·U` + Slater-swap sparsity witness `nnz_off/row ≤ (1 + p_r·(N−p_r))·(1 + p_c·(N−p_c)) − 1`; density falls from ~37 % at N=4 to ~10 % at N=7, extrapolated ~1 % at N=10. Verifies the bound where dense `U^†LU` still fits. |
| `JwSlaterPairSparseLBuilder` | Direct sparse CSR construction of L_JW element-by-element from the analytic `Z_l ⊗ Z_l` action on Slater pairs, no dense `U^†LU` detour. At N=10 (p_c=5, p_r=5, sectorDim 63 504): ~14 M nnz, ~280 MB, ~3 s build. **This is the N=10 enabler.** |
| `JwSlaterPairArnoldiEig` | Managed non-symmetric Arnoldi with Modified Gram-Schmidt on the CSR matvec; returns top-k largest-magnitude Ritz values. Matvec parallelised via `Parallel.For` over rows. |
| `JwSlaterPairShiftInvertArnoldi` | Same outer Arnoldi on `(L − σI)^(−1)`; inner solver is a Jacobi-preconditioned in-house BiCGStab on the same CSR matvec, no MathNet sparse conversion. Recovers eigenvalues nearest σ, including the steady state at λ = 0 (machine precision at N=10). |
| `JwSlaterPairF1PalindromeProbe` | Two-shift F1 probe: runs ShiftInvertArnoldi at σ_slow ≈ 0 AND σ_fast = F1(σ_slow) = −σ_slow − 2·Σγ, reports per-pair residual on the F1 mirror map λ → −λ − 2·Σγ with the L-conjugation degree of freedom absorbed. At N=10 (5,5) the steady-state pair matches at machine precision (~3.6e-14), the two slowest decay-mode pairs at ~3e-7 / ~3e-6; deeper Ritz values may not pair (Arnoldi-depth limit). Informational, not a strict F1 prover. |

Shared low-level complex-vector kernels (`NormSquared`, `ConjugateDot`, `AxpyInPlace`, `RandomNormalized`) live in `KrylovOps.cs` to keep the Arnoldi sites lean.

The combinatorial table above is unchanged; joint-popcount + F71 + X⊗N still do not shrink the max block enough for dense LAPACK at N=10, but the matrix-free Krylov path renders that limit moot for top-k / near-σ work.

**What this path does NOT yet give:**
- **Full N=10 spectrum** (all ~63 k eigenvalues per max-block, ~1 M total). For F1 palindrome verification at N=10 as a structural theorem-check the cleaner route is Prosen third quantization (Phase 3 of the N=10 plan): under uniform γ, the chain XY + Z-dephasing Liouvillian reduces to a 2N × 2N Nambu spectral problem (Medvedyeva, Essler, Prosen 2016). Phase 3 has not yet been implemented.
- **Strict F1 pair-matching for deeper Ritz values.** `JwSlaterPairF1PalindromeProbe` recovers F1-paired sets reliably for the closest-to-σ eigenvalues (steady state, first two slow modes) but at finite Krylov dimension the deeper Ritz values returned by each side's shift-invert Arnoldi can land on different "approximately K-closest" subsets, breaking pair-matching at the K-th probe slot even though F1 holds at the theorem level. This is an Arnoldi-depth limit, not an F1 break.
- **Inhomogeneous γ spectrum** at N=10 is supported by the sparse path but not yet exercised in a science-facing pipeline; Prosen's reduction does not extend there.

**Phase 2 timing at N=10 (5,5)** (Jacobi-preconditioned BiCGStab, commodity hardware): sparse build ~3 s, two-shift probe with K=4 per end at numIter=20, inner tol 1e-8 ~8 s. Inner BiCGStab mean iter dropped from ~866 (pre-Jacobi) to ~22 (post-Jacobi commit `45e6a40`), a ~40× wall-clock improvement.

## Phase 3b: Klein 4-group internal symmetry on self-paired sectors (2026-05-17)

The X⊗N global charge-conjugation primitive `XGlobalChargeConjugationPairing` halves the
number of distinct sector eig-calls by pairing (p_c, p_r) ↔ (N−p_c, N−p_r). At even N
the sector (N/2, N/2) is X⊗N-self-paired, its pair *is itself*, so the inter-sector
halving cannot apply. But X⊗N is still a Z₂ symmetry of that sector's Liouville-space
basis, and combined with the F71 chain spatial mirror (also Z₂, commutes with X⊗N) it
generates the Klein four-group K = {1, F71, X⊗N, F71·X⊗N} acting *internally* on the
(N/2, N/2) basis. K splits the sector into 4 character sub-blocks (++, +-, -+, --).

This closes the gap between the inventory's advertised "≈ 16k max-block at N=10" (which
was implicitly assuming the Klein splitting) and the previously delivered 31752 dim from
F71 alone. The Klein refinement primitive
`compute/RCPsiSquared.Core/BlockSpectrum/KleinFourGroupSelfPairedRefinement.cs` builds:

| Klein character | (++) | (+-) | (-+) | (--) | sum |
|-----------------|------|------|------|------|-----|
| Sub-block dim at N=10 (5, 5) | 16132 | 15620 | 15620 | 16132 | 63504 ✓ |

Per sub-block element-wise construction (no full sector L materialised): ~19 s per sub-block
at N=10 (5, 5), mean nnz per row 10.7, **sparse**, 0.07 % density, ~170 k nnz vs 260 M
dense entries. Cross-validated at N=4 (2, 2), N=6 (3, 3), N=8 (4, 4): union of the 4
sub-block dense Evds matches the direct `PerBlockLiouvillianBuilder` sector Evd as a
multiset within 1e-9.

**Where the literature does not go.** Medvedyeva-Essler-Prosen (2016) reach N=10 via
the imaginary-U Bethe ansatz, analytic, doesn't need K splitting. For our
computational path, K on self-paired sectors is repo-specific (F71 spatial Z₂ ×
X⊗N charge Z₂ both typed primitives, combination unique to this codebase).

**Phase 3c (delivered).** The dense Evd on Klein sub-blocks at dim 16132 is blocked by
the .NET/MKL int32 array-size cap (4.2 GB > 2 GB managed marshaling limit), but the
sub-blocks are *sparse* (mean 10.86 nnz per row at N=10 (5, 5), 0.017 % density,
~170 k nnz per sub-block, ~7 MB sparse storage vs 4.2 GB dense). Phase 3c stores
each sub-block as CSR and extracts top-K slow modes via per-sub-block shift-invert
Arnoldi:

| Primitive | Role |
|---|---|
| `KleinFourGroupSelfPairedSparseLBuilder` | Direct sparse-CSR build of one Klein character sub-block, element-wise reconstruction of `KleinFourGroupSelfPairedRefinement.BuildSubBlockL`. At N=10 all 4 sub-blocks build in ~0.4 s total. |
| `SparseShiftInvertArnoldi` | Generic shift-invert Arnoldi on any CSR matrix (extracted from `JwSlaterPairShiftInvertArnoldi`'s algorithm). Inner BiCGStab with Jacobi preconditioning; outer Modified Gram-Schmidt Arnoldi. Reusable by both the JW Slater-pair pipeline AND the Klein sub-block pipeline (a single source of algorithmic truth). `KrylovOps` promoted from `internal` to `public` to share the complex-vector primitives. |

At N=10 (5, 5), top-4 slow modes per Klein sub-block (16 modes total) in ~22 s wall.
The slow modes match Phase 2 `JwSlaterPairF1PalindromeProbe` on the full sector for
the modes both runs sample (e.g. the steady state at (0, 0) and the (−0.200, 0)
slow-decay mode in ++; (−0.182, 0) in +-) AND extend coverage with Klein-character-
specific modes inaccessible to Phase 2 (e.g. (−0.336) in +-, (−0.313) in -+).

**Phase 3d (closed): correct preconditioner choice → strict eigenvalues.** The initial
Phase 3c run used Jacobi preconditioning by default and observed BiCGStab saturating
at the iteration cap (~1000 iter, vs Phase 2 JW's ~22 mean), Ritz values were
APPROXIMATE (not strict eigenvalues), a "trap" per Tom's principle that nothing is
imprecise in the quantum domain. Root cause diagnosed: Klein computational-basis
sub-blocks have diagonal = only −2γ·hamming (max magnitude 1.0; many rows exactly 0
for diagonal-projector orbits) and off-diagonal = ±iJ (magnitude 1.0). The matrix is
NOT diagonally dominant, so Jacobi preconditioner amplifies near-zero shifted diagonals
to ~1/σ ≈ 1000 at σ = (0, 0.001), destabilising BiCGStab. The JW basis used in Phase 2
has ε-dominant diagonals (Hamiltonian eigen-energy differences on the diagonal), making
Jacobi appropriate there. The fix is structural: `SparseShiftInvertArnoldi.Run`
exposes a `PreconditionerKind` parameter (default `Jacobi` for backward compatibility,
`Identity` for non-diagonally-dominant matrices). With `PreconditionerKind.Identity`
and numIter=80 at N=10 (5, 5), all 16 slow modes converge to machine precision
(1e-10 to 1e-13 distance from dense eigenvalues at small N), and the recovered slow
modes exactly reproduce Phase 2's outputs:

| Phase 2 slow mode (full sector) | Klein sub-block | Recovered (Phase 3d) |
|---|---|---|
| (0, 0) steady state | ++ | (0.0, 0.0) machine precision |
| (−0.182, 0) | +- | (−0.182, 0) × 4 (degenerate) |
| (−0.200, 0) | ++ | (−0.200, 0) × 2 (degenerate) |
| (−0.163, ±0.239) | -- | (−0.163, ±0.238) conjugate pair |

Plus Klein-character-specific modes inaccessible to Phase 2 (e.g. (−0.288) in ++,
(−0.287, ±0.069) in -+, (−0.322, ±0.070) in --). Wall time at N=10 (5, 5):
~65 s total for 16 slow modes across 4 sub-blocks (vs Phase 2's 8 s for 4 modes).
4× spectral coverage at ~8× wall time, all strict eigenvalues.

## Phase 3a: Prosen rapidities for the one-sided sectors (2026-05-17)

The (p_c = 0, p_r = m) sectors of chain XY + uniform Z-dephasing admit a closed-form
spectrum without any diagonalisation. The dissipator reduces to a constant `−2γm`
(Hamming distance of any m-fermion ket from the vacuum bra is exactly m), and `H|0⟩ = 0`
removes the unitary contribution from the bra side, so `L_(0,m) = −iH_m − 2γm·I` and

```
λ_S  =  Σ_{k∈S} β_k  =  −2γm − i·Σ_{k∈S} ε_k     for  S ⊆ {1..N},  |S| = m
```

where `β_k = −2γ − i·ε_k` are the **N Prosen rapidities** of the model and
`ε_k = 2J·cos(πk/(N+1))` is the OBC sine-mode dispersion. This is the simplest leaf of
Medvedyeva-Essler-Prosen's (2016) imaginary-Hubbard programme; the (m, m̃) sectors with
m, m̃ ≥ 1 require the imaginary-U Bethe ansatz (not yet implemented in the repo).

`compute/RCPsiSquared.Core/BlockSpectrum/Prosen/OneSidedSectorClosedForm.cs` returns the
rapidities + full sector spectrum as subset sums. Tests cross-validate against
`PerBlockLiouvillianBuilder.BuildBlockZ` dense Evd (multiset match within 1e-10 at
N=3..6), verify the rapidities match `XyJordanWignerModes.Dispersion`, confirm the X⊗N
charge-conjugation pair (p_c = N, p_r = N−m) has equal spectrum, and check the F1 mirror
predictions `−λ − 2·Σγ` of every (0, m) eigenvalue appear somewhere in the full per-block
spectrum (F1's Π-mediated mirror distributes across sectors, not a simple sector
permutation).

**Combinatorial coverage at N=10.** Summing C(10, m) for m = 0..N gives 2^N = 1024
eigenvalues per one-sided family. With sister sectors (m, 0) → complex conjugates of
(0, m), and X⊗N pairs (N, N−m) and (N−m, N) → equal spectra, four 1024-eigenvalue
families are accessible analytically: ~4 k of the ~1 M total N=10 eigenvalues, all in
closed form. Each (0, m) sector sits at Re λ = −2γm so spans the slow (m=0, steady
state at λ=0) to fast (m=N=10, λ=−2Σγ = −1.0 at γ=0.05) range uniformly.

## Cross-references

- `JointPopcountSectors`, `F71MirrorBlockRefinement`, `F71BilateralBlockRefinement`, `F71AntiPalindromicGammaSpectralInvariance`: see `compute/RCPsiSquared.Core/BlockSpectrum/`.
- `F1`, `F71`, `F91`, `F92`, `F93`, `F100`, `F101`: see `docs/ANALYTICAL_FORMULAS.md`.
- F91/F92/F93 algebraic proofs: `docs/proofs/PROOF_F91_GAMMA_NINETY_DEGREES.md`, `docs/proofs/PROOF_F92_BOND_ANTI_PALINDROMIC_J.md`, `docs/proofs/PROOF_F93_DETUNING_ANTI_PALINDROMIC.md`.
- F100 (observable-side twin of F92): `docs/proofs/PROOF_F100_C1_QPEAK_MIRROR_J_PARITY.md`. The c₁/Q_peak bond-mirror deviation is exactly odd in the F71-anti-palindromic J component; where F92 keeps J_anti out of the diagonal-block spectrum, F100 localises the entire bond-mirror deviation in J_anti. Typed as `C1QPeakMirrorJParity` in `compute/RCPsiSquared.Core/F71/`.
- F101 (observable-side twin of F91): `docs/proofs/PROOF_F101_C1_MIRROR_GAMMA_PARITY.md`. The c₁ bond-mirror deviation is exactly odd in the F71-anti-palindromic component of the per-site γ profile; where F91 keeps γ_anti out of the diagonal-block spectrum, F101 localises the entire bond-mirror deviation in γ_anti. c₁ only (the F86c Q_peak observable needs a scalar γ₀). Typed as `C1MirrorGammaParity` in `compute/RCPsiSquared.Core/F71/`.
- Phase 2 N=10 push primitives: `compute/RCPsiSquared.Core/BlockSpectrum/JordanWigner/JwSlaterPairBasis.cs`, `JwSlaterPairLProjection.cs`, `JwSlaterPairSparseLBuilder.cs`, `JwSlaterPairArnoldiEig.cs`, `JwSlaterPairShiftInvertArnoldi.cs`, `JwSlaterPairF1PalindromeProbe.cs`, `KrylovOps.cs`.
- Phase 3a Prosen leaf: `compute/RCPsiSquared.Core/BlockSpectrum/Prosen/OneSidedSectorClosedForm.cs`.
- Phase 3b Klein-4 refinement: `compute/RCPsiSquared.Core/BlockSpectrum/KleinFourGroupSelfPairedRefinement.cs`.
- Phase 3c sparse Klein + shift-invert: `compute/RCPsiSquared.Core/BlockSpectrum/KleinFourGroupSelfPairedSparseLBuilder.cs` + `compute/RCPsiSquared.Core/BlockSpectrum/SparseShiftInvertArnoldi.cs`.
- Synthesis: `reflections/ON_THE_SYMMETRY_FAMILY.md`, `reflections/ON_THE_NINETY_DEGREE_GAMMA.md`.
