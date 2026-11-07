# The N = 10 block-spectrum push: the one-sided closed form, the sparse Slater-pair route, and the Klein sub-blocks

<!-- STAGE 0 SWEEP RECORD (2026-11-07). This page holds the N = 10 computation that docs/SYMMETRY_FAMILY_INVENTORY.md carried until that page became an inventory of symmetries alone; the symmetries it uses are rows there.
docs/ANALYTICAL_FORMULAS.md: F1 (the mirror), F2b (the sine modes), F125 and F142 (the two ladders), F140 (the frozen divisor's corner block), F153 (the pinned boundary blocks), F171 (the grey rooms, the ladder images counted).
docs/proofs/: PROOF_F1_GENERAL_TOPOLOGY (the N = 9 chain through the native MklDirect path, the N = 10 frontier named there), PROOF_FROZEN_BAND_SO4 (the ladders, section 2), PROOF_R90_FROZEN_DIVISOR (section 9.2, the population-free half of (1, 1)), PROOF_CODIM1_BY_ADDITIVITY (Theorem B, the spectator), MIRROR_SYMMETRY_PROOF.
docs/proofs/derivations/: D10 is the Heisenberg chain's vacuum-coherence block |vac⟩⟨ψ| (Step 3, the graph Laplacian), whose XY counterpart is the adjoint of the closed form below at m = 1, with F2b's sine modes; F2's note draws that line, and F2b's body places the sine modes as the frequencies of the vacuum coherences.
docs/LITERATURE_REVIEW.md: Medvedyeva, Essler and Prosen (2016), the imaginary-interaction Hubbard map on any bipartite lattice and the Bethe ansatz in one dimension; Prosen (2012), the PT symmetry behind the palindrome; that the paper's Bethe equations are the periodic chain's is recorded in docs/CAUGHT_ERRORS.md (2026-11-04).
experiments/: SLOW_MODE_R_PARITY (reads JwBlockBasis, the c = 2 case the Slater-pair basis generalises), THE_PAIRING_CONDITION; no flight concerns N = 10.
fw.Confirmations: nothing at N = 10 (the longest hardware line is eight qubits, f129_standing_fringe_kingston_july2026).
OpenArcs registry: site_resolved_vacuum_block (the vacuum-coherence block derived from below three times without cross-citation, per the arc; this page routes the XY case to F2b and the Heisenberg case to D10); xxz_axis_handover (the interacting Bethe ansatz absent from the repo).
docs/GLOSSARY.md: the fermion section (Jordan-Wigner fermions, the spectator W, η-pairing).
docs/CAUGHT_ERRORS.md: 2026-11-01 (the N = 10 max-block after the splits is 26 460, and the 16 132-row dense solve is past MathNet's managed limit but not past the native path), 2026-11-04 (no Nambu reduction holds at every N; the Bethe ansatz is the tool for the interacting model).
Typed layer: Core/BlockSpectrum/Prosen, Core/BlockSpectrum/JordanWigner, the Klein builders in Core/BlockSpectrum, BlockSpectrumOpenQuestions (the open question on the two-sided sectors), Diagnostics GreyRoomLadderImagesClaim; the tests named below print the numbers quoted.
Adjacency: the closed form at m = 1 is the adjoint of the XY chain's vacuum-coherence block (F2b), the counterpart of D10's Heisenberg one; the Klein sub-block builder and the F1 general-topology page describe the same 64 GB block from two sides. -->

**What this is.** At N = 10 the largest joint-popcount sector of the chain XY + Z-dephasing Liouvillian, (5, 5), has C(10, 5)² = 63 504 rows, and its dense eigendecomposition needs 64 GB of native memory per block, which fits the 128 GB dev machine only one block at a time and exhausts the parallel budget of the per-block path that reached N = 9 ([the F1 general-topology proof](../docs/proofs/PROOF_F1_GENERAL_TOPOLOGY.md)). This page records what reaches N = 10 instead: the one-sided sectors in closed form, the slow modes of the full (5, 5) sector through a sparse Jordan-Wigner route, and the Klein split of that sector into four sparse sub-blocks with slow modes of their own. It ends with what the route does not give. The symmetries it uses are rows of [the symmetry family inventory](../docs/SYMMETRY_FAMILY_INVENTORY.md).

**Setting.** N = 10, open chain, H = (J/2)·Σ_b (X_b X_{b+1} + Y_b Y_{b+1}) with J = 1, uniform rates γ = 0.05, so σ = Σγ = 0.5 and the far end of the palindrome sits at −2σ = −1.0. Sectors are (p_c, p_r) = (bra popcount, ket popcount). The measurements below are printed by tests in `compute/RCPsiSquared.Core.Tests/BlockSpectrum/`, named at the number, from two Debug-build runs of 2026-11-07 on the dev machine; the counts and sizes (675, 4092, 26 460, 15 876, the memory figures) are arithmetic or come from the pages cited.

## 1. The one-sided sectors in closed form

In a sector (0, m) the bra is the vacuum. H annihilates it (XX + YY on |0…0⟩ is zero), so the Hamiltonian part acts on the ket alone, and the Hamming distance of any m-fermion ket from the vacuum bra is m, so the dissipator is the constant −2γm:

    L_(0,m) = −i·H_m − 2γm·I,
    λ_S = Σ_{k∈S} β_k = −2γm − i·Σ_{k∈S} ε_k    for S ⊆ {1..N}, |S| = m,

with ε_k = 2J·cos(πk/(N+1)) the open-chain sine-mode dispersion ([F2b](../docs/ANALYTICAL_FORMULAS.md#f2b)) and β_k = −2γ − i·ε_k the N Prosen rapidities. No diagonalisation. This is the simplest leaf of Medvedyeva, Essler and Prosen's imaginary-interaction Hubbard programme ([literature review](../docs/LITERATURE_REVIEW.md)).

`compute/RCPsiSquared.Core/BlockSpectrum/Prosen/OneSidedSectorClosedForm.cs` returns the rapidities and the sector spectrum as subset sums. `OneSidedSectorClosedFormTests` cross-validates the spectrum against the dense block eigendecomposition of `PerBlockLiouvillianBuilder.BuildBlockZ` at N = 3 to 6 (multiset match within 1e-10), the rapidities against `XyJordanWignerModes.Dispersion`, the X⊗N partner (N, N−m) for equal spectra, and that the F1 images −λ − 2σ of every (0, m) eigenvalue appear in the union of the dense per-block spectra at (N, m) = (3, 1), (4, 2), (5, 2); Π's sector map places them in (N−m, 0), Π sending |a⟩⟨b| to |b⟩⟨ā|.

At m = 1 the closed form is the block of the single-excitation coherences |ψ⟩⟨vac|, N eigenvalues −2γ − i·ε_k; its adjoint |vac⟩⟨ψ_k|, the page's (1, 0), carries the conjugate values and is the block whose frequencies F2b names, the sine modes. The Heisenberg chain's |vac⟩⟨ψ| block, with the degree diagonal its ZZ term adds to the hopping, the graph Laplacian, is [D10](../docs/proofs/derivations/D10_W1_DISPERSION.md)'s and [F2](../docs/ANALYTICAL_FORMULAS.md#f2)'s; F2's note draws the line between the two, and F2b's body places its sine modes as the frequencies of the |vac⟩⟨ψ_k| coherences.

**Coverage at N = 10.** Σ_m C(10, m) = 2^10 = 1024 eigenvalues per one-sided family; the sister family (m, 0) carries the complex conjugates, and the X⊗N partners (N, N−m) and (N−m, N) the same spectra, so four families of 1024 over 40 sectors, sharing the four corner sectors, give 4092 eigenvalues in closed form, counted with sector multiplicity (the X⊗N partners repeat the values, and so does the conjugate family, each (0, m) multiset being closed under conjugation by the chiral pairing ε_{N+1−k} = −ε_k), out of 4^10 = 1 048 576. Each (0, m) sits at Re λ = −2γm = −0.1·m, from the steady state at λ = 0 to the far end −1.0 at m = 10.

**Pieces of the two-sided sectors that close without the Bethe ansatz.** The images of the boundary blocks (k, 0) under the two ladders of [PROOF_FROZEN_BAND_SO4](../docs/proofs/PROOF_FROZEN_BAND_SO4.md) §2, Φ = W and the transpose S⁻ of S⁺, which lie in the grey rooms, of dimension at least C(N, k) at −2γk in every block whose Hamming window holds k ([F171](../docs/ANALYTICAL_FORMULAS.md#f171), `inspect --root greyroom`); and the half of (1, 1) that holds no populations, which at every rate profile is, up to the shift −4γ̄ (γ̄ the mean rate), the operator ρ' ↦ aρ' + ρ'a on Λ²(ℂᴺ), a the single-particle matrix the proof calls m and ρ' the block state times the stagger, its eigenvalues the pair sums of a's eigenvalues minus 4γ̄ ([PROOF_R90_FROZEN_DIVISOR](../docs/proofs/PROOF_R90_FROZEN_DIVISOR.md) §9.2).

## 2. The sparse Slater-pair route to the full (5, 5) sector

Six primitives under `compute/RCPsiSquared.Core/BlockSpectrum/JordanWigner/`, with the complex-vector kernels (`NormSquared`, `ConjugateDot`, `AxpyInPlace`, `RandomNormalized`) in `KrylovOps.cs`:

| Primitive | Role |
|---|---|
| `JwSlaterPairBasis` | The basis change U from the computational-basis sector (p_c, p_r) to the Jordan-Wigner Slater-pair basis, in which L_H is diagonal with eigenvalues −i·(Σε(L) − Σε(K)). |
| `JwSlaterPairLProjection` | The dense U†·L·U where it still fits, and the Slater-swap sparsity witness: the dissipator, quartic in the fermions, acts on each side by a single orbital swap, the support being the Cartesian product of the two sides' swaps, so nnz_off per row ≤ (1 + p_r·(N−p_r))·(1 + p_c·(N−p_c)) − 1, which at N = 10, (5, 5) is 675 of 63 504, about 1.1 %. |
| `JwSlaterPairSparseLBuilder` | The direct sparse CSR build of L in the Slater-pair basis from the analytic action of Z_l ⊗ Z_l on Slater pairs, without the dense detour. At N = 10, (5, 5): 13 993 504 nonzeros, at most 386 per row, built in about six seconds (one Debug-build run, 2026-11-07). This is what makes N = 10 reachable. |
| `JwSlaterPairArnoldiEig` | Managed non-symmetric Arnoldi with modified Gram-Schmidt on the CSR matvec (rows in parallel), returning the largest-magnitude Ritz values. |
| `JwSlaterPairShiftInvertArnoldi` | The same Arnoldi on (L − shift·I)⁻¹, the inner solves by a Jacobi-preconditioned BiCGStab on the same matvec; eigenvalues nearest the shift, the steady state at λ = 0 among them for a shift near 0. In this basis the diagonal carries the energy differences −i·(Σε(L) − Σε(K)) beside the dissipator's part, and the Jacobi-preconditioned inner solves converge in about 22 iterations (section 3 says where Jacobi fails). |
| `JwSlaterPairF1PalindromeProbe` | The two-shift F1 probe: shift-invert runs at a slow shift near 0 and at the fast shift, its mirror −(slow shift) − 2σ, and per slow value the residual, the least distance from the mirror −λ − 2σ of the value or of its conjugate to any recovered fast value. |

The probe at N = 10, (5, 5), four Ritz values per end, twenty Arnoldi steps, inner tolerance 1e-8, the inner BiCGStab averaging 22 iterations per solve, the whole probe about 24 s in the same run (`JwSlaterPairF1PalindromeProbeTests.Build_AtN10HalfFilling_RunsAndRecoveredEigenvaluesAreBounded`):

| Slow end, as printed | Predicted mirror −λ − 2σ | Fast end, as printed, in slot order | Residual, the least distance to any fast value |
|---|---|---|---|
| (−0.000000, −0.000000), the steady state | (−1.000000, 0.000000) | (−1.000000, −0.000000) | 3.58e-14 |
| (−0.181763, −0.000000) | (−0.818237, 0.000000) | (−0.818237, 0.000000) | 2.78e-07 |
| (−0.199997, −0.000000) | (−0.800003, 0.000000) | (−0.800007, 0.000000) | 3.62e-06 |
| (−0.162240, 0.238781) | (−0.837760, −0.238781) | (−0.732455, 0.000558) | 0.24 (unpaired) |

The steady state pairs with the far end at residual 3.58e-14; the next two pair at 2.78e-7 and 3.62e-6, the depth the twenty-step Krylov space reaches; the fourth slot does not pair: its residual 0.24 is the distance to the second fast value, and the fast end's own fourth value, (−0.732455, 0.000558), is a twenty-step Ritz value whose mirror near −0.2675 is none of the sixteen sub-block values of section 3, nothing here showing it to be an eigenvalue. The probe is informational, not a proof of the pairing: at finite Krylov dimension the two ends can land on different sets of four nearest values, so a deeper slot may fail to pair even though F1 holds at the theorem level. The slow values are twenty-step Ritz values: section 3's eighty-step sub-block run places −0.199997 at −0.200000, the fourth slow value (−0.162240, 0.238781) at (−0.162803, 0.237879), 1.06e-3 away, and −0.181763 within 2e-5 to 4e-5 of the four +− values.

## 3. The Klein split of (5, 5) and the slow modes per sub-block

X⊗N pairs (p_c, p_r) with (N−p_c, N−p_r) and so pairs (5, 5) with itself; the pairing cannot halve it. But X⊗N is still a Z₂ on the sector's basis, it commutes with the F71 site reflection, and the Klein group K = {1, F71, X⊗N, F71·X⊗N} acts inside the sector and splits it into four character sub-blocks:

| Klein character | ++ | +− | −+ | −− | sum |
|---|---|---|---|---|---|
| Sub-block dimension at N = 10, (5, 5) | 16 132 | 15 620 | 15 620 | 16 132 | 63 504 |

`KleinFourGroupSelfPairedRefinement.cs` builds the orbits and the sub-blocks; `KleinFourGroupSelfPairedRefinementTests` asserts these dimensions and, at N = 4, 6, 8, that the union of the four sub-block spectra equals the sector's dense spectrum as a multiset within 1e-9. The largest block of the whole N = 10 problem is then (4, 5) at 26 460 rows ([the inventory](../docs/SYMMETRY_FAMILY_INVENTORY.md), the N = 10 section).

**Sparse, in the computational basis.** `KleinFourGroupSelfPairedSparseLBuilder.cs` builds one sub-block directly as CSR. At N = 10 the four sub-blocks hold 175 261, 169 645, 169 681 and 175 297 nonzeros, a mean of 10.86 to 10.87 per row and at most 18, built together in under a second (`Build_AtN10_AllFourSubBlocksHaveExpectedSparsity`), against 16 132² ≈ 260 M dense entries; the CSR is a few megabytes where the dense sub-block is 4.2 GB (16 132² × 16 B), past MathNet's 2 GB managed marshaling limit. The native MklDirect path that solved the N = 9 chain's largest block (15 876 rows) passes that limit and has not been run on these sub-blocks; the test that would do it, `SLOW_Build_AtN10_FullSectorSpectrumViaFourSubBlockEvd_F1PalindromeHolds`, is skipped.

**The preconditioner.** `SparseShiftInvertArnoldi.cs` is a generic shift-invert Arnoldi on any CSR matrix with a `PreconditionerKind` parameter; the Slater-pair route still runs its own `JwSlaterPairShiftInvertArnoldi`, and the two share only `KrylovOps`. In the computational basis a Klein sub-block has diagonal −2γ·(Hamming distance), many rows exactly zero, and off-diagonal entries of magnitude comparable to J: not diagonally dominant. A Jacobi preconditioner on (L − 0.001i·I) then amplifies the near-zero diagonals by about 1/0.001 = 1000 and the inner BiCGStab stalls at its iteration cap, which leaves approximate Ritz values. `PreconditionerKind.Identity` keeps the natural conditioning. That the eighty-step Ritz values are then eigenvalues is witnessed at N = 6 (++ and −−) and N = 8 (++) with inner tolerance 1e-12: each recovered value lies within 1e-6 of an eigenvalue of the dense sub-block (`Build_SparseShiftInvertWithIdentityPrecond_RecoversGenuineEigenvalues`; the distances printed on 2026-11-07 run from 1e-16 to 4e-14, and 7.7e-8 for the fourth value of the N = 8 ++ sub-block). At N = 10 there is no dense spectrum to compare with: the values are Ritz values, the inner solves ran to a relative residual of 1e-8 or to the cap of 1000 iterations (the printout gives the mean counts, 570 to 728, not the maximum), and no Arnoldi residual is printed. The Slater-pair basis keeps Jacobi, its diagonal being the energy differences.

**The sixteen slow modes at N = 10** (`Build_AtN10_SparseShiftInvertExtractsSlowModesPerSubBlock`, shift 0.001i, four Ritz values per sub-block, eighty Arnoldi steps, Identity preconditioner, inner tolerance 1e-8 and at most 1000 inner iterations; one Debug-build run of 2026-11-07, the four sub-blocks built in under a second together and the four solves taking 13.6 minutes in all with the inner BiCGStab averaging 570 to 728 iterations per solve, a reading of the cost and not a benchmark):

| Sub-block | Recovered λ, as printed (six decimals) |
|---|---|
| ++ | (0.000000, 0.000000), the steady state; (−0.200000, 0.000000) twice; (−0.287688, −0.000001) |
| +− | (−0.181723, −0.000000), (−0.181782, −0.000000), (−0.181797, 0.000000), (−0.181802, −0.000000) |
| −+ | (−0.286576, ±0.069437), (−0.290573, ±0.067426) |
| −− | (−0.162803, ±0.237879), (−0.322187, ±0.070162) |

Two readings of this table. The four +− values lie within 1e-4 of −0.1818, and the full-sector probe's −0.181763 lies 2e-5 to 4e-5 from each; whether they are one fourfold eigenvalue or four close ones is not decided by this route, since at N = 10 there is no dense spectrum to compare with. And the modes the two routes both sample agree to 3e-6 on the real value, −0.200000 against the probe's −0.199997, and to 1.06e-3 on the complex pair, (−0.162803, 0.237879) against the probe's (−0.162240, 0.238781), the steady state agreeing at 1e-14. A single sweep on the full sector returns four values and the sub-block runs four each, sixteen in all; the printouts show proximity and nothing more: the probe's −0.181763 lies near the four +− values and its −0.199997 near the two −0.200000 values of ++, and nothing printed gives a multiplicity.

## What this route does not give

- **The full N = 10 spectrum**, all ~63 k eigenvalues of (5, 5) and ~1 M in all. The palindrome itself needs no spectrum: Π, and Prosen's PT relation, prove it at every N. The route the literature offers is Medvedyeva, Essler and Prosen (2016): the dephased tight-binding model on any bipartite lattice, and through the Jordan-Wigner transform the XY chain with open ends, becomes a Hubbard model with imaginary interaction, Bethe-ansatz integrable in one dimension. In the Jordan-Wigner fermions Z_l = 1 − 2n_l is quadratic, so Z_l ρ Z_l is quartic, and no reduction quadratic in those fermions holds at every N: the even-parity sector of a quadratic generator is point-symmetric about its mean, every even subset of the 2N rapidities having an even complement, and the Liouvillian's is not at N = 3 and 5 (J = 1, γ = 0.3, [`simulations/even_sector_point_symmetry.py`](../simulations/even_sector_point_symmetry.py)); at other odd N this is not proven, and at even N the test is silent. The paper writes the Bethe equations for the periodic chain; the open chain becomes the imaginary-interaction Hubbard chain with open ends, whose Bethe equations have not been taken up here, and completeness at imaginary interaction is not established here. The open question is kept in `BlockSpectrumOpenQuestions.cs`.
- **The F1 pairing of the deeper Ritz values** of the two-shift probe (section 2).
- **Inhomogeneous γ at N = 10**: the sparse builders take a profile; the Klein builder's rates are indexed by bit, bit l being site N−1−l (a trap the first ledger entry of 2026-08-20 records), which is harmless exactly where the sub-blocks are blocks of L at all, on rates and couplings symmetric under the site reflection, since the F71 factor of K needs that symmetry of H and of the rates; no science-facing run has used a profile, and the Bethe solution assumes uniform γ and J.
- **The dense eigendecomposition of a Klein sub-block** through the native path: not run.

## Where things live

- Closed form: `compute/RCPsiSquared.Core/BlockSpectrum/Prosen/OneSidedSectorClosedForm.cs`, tests `compute/RCPsiSquared.Core.Tests/BlockSpectrum/Prosen/`.
- Slater-pair route: `compute/RCPsiSquared.Core/BlockSpectrum/JordanWigner/` (the six primitives and `KrylovOps.cs`), tests `compute/RCPsiSquared.Core.Tests/BlockSpectrum/JordanWigner/`; the modes themselves in `compute/RCPsiSquared.Core/F86/JordanWigner/XyJordanWignerModes.cs`.
- Klein route: `compute/RCPsiSquared.Core/BlockSpectrum/KleinFourGroupSelfPairedRefinement.cs`, `KleinFourGroupSelfPairedSparseLBuilder.cs`, `SparseShiftInvertArnoldi.cs`, tests `compute/RCPsiSquared.Core.Tests/BlockSpectrum/KleinFourGroupSelfPaired*.cs`.
- The symmetries used: [SYMMETRY_FAMILY_INVENTORY.md](../docs/SYMMETRY_FAMILY_INVENTORY.md); the N = 9 dense path and the 64 GB block: [PROOF_F1_GENERAL_TOPOLOGY.md](../docs/proofs/PROOF_F1_GENERAL_TOPOLOGY.md).
