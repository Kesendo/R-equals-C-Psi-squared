# PROOF F93: F71-anti-palindromic h spectral invariance (h-detuning Pi2-Z₄ twin of F91/F92)

**Status:** Tier 1 derived (algebraic argument parallel to F91/F92 + bit-exact empirical witness at N=4, 5)
**Date:** 2026-05-12
**Authors:** Thomas Wicht, Claude (Anthropic)
**Typed claim:** [`F93DetuningAntiPalindromicSpectralInvariance.cs`](../../compute/RCPsiSquared.Core/SymmetryFamily/F93DetuningAntiPalindromicSpectralInvariance.cs)

---

## Statement

For the chain XY + per-site Z-detuning + uniform Z-dephasing Liouvillian L on N qubits (Hamiltonian H = (1/2) Σ_b J_b · (X_b X_{b+1} + Y_b Y_{b+1}) + Σ_l h_l · Z_l with uniform J, inhomogeneous h), the **F71-refined diagonal-block eigenvalue multiset** is invariant under any h-distribution satisfying

    h_l + h_{N−1−l} = 2·h_avg = (2/N)·Σ_l h_l   for all l

i.e. **h is F71-anti-palindromic around its mean**.

**Scope restriction:** only **longitudinal** h_l Z_l detuning is in scope. Transverse detuning h_l X_l or h_l Y_l flips popcount and breaks the joint-popcount conservation that the entire BlockSpectrum framework rests on; F93 does not apply there.

## Pi2-Z₄ structure (parameter-side, h-axis)

Same Z₄ as F91 (γ-axis) and F92 (J-axis). The 90°-rotation h_l ↦ 2·h_avg − h_{N−1−l} preserves the anti-palindromic orbit pair-sum h_l + h_{N−1−l} = 2·h_avg.

## Algebraic mechanism (parallel to F91/F92)

The h_l Z_l detuning contributes to L diagonally in computational basis (Z_l is diagonal):

    [h_l Z_l, |a⟩⟨b|] = h_l (z_l(a) − z_l(b)) · |a⟩⟨b|

where z_l(x) = +1 if bit l of x is 0, else −1. The contribution is purely diagonal in Liouville-space.

In the F71-refined basis: diagonal-block matrix elements ⟨sym|h_l Z_l|sym⟩ pick up the average over the F71-orbit pair (a, b) and (a', b') = (F71(a), F71(b)). The h-coefficient that enters is h_l + h_{N−1−l} (since F71 maps site l ↔ N−1−l, and z_l(a') = z_{N−1−l}(a) in the orbit average).

Therefore **diagonal-block matrix elements depend on h only through pair-sums** h_l + h_{N−1−l}; cross-block entries depend on pair-differences h_l − h_{N−1−l}.

The XY chain Hamiltonian (uniform J) is F71-symmetric per Step 5 of PROOF_F91; it contributes the same to every F71-block independent of h. The uniform Z-dephasing contributes the same diagonal shift to every F71-block independent of h.

Anti-palindromy h_l + h_{N−1−l} = 2·h_avg is the orbit on which all h give identical diagonal-block spectra equal to uniform h_avg.

## Empirical witness

| N | uniform h = 0.4 | anti-palindromic h | full-L Frobenius leak (F71 off-block) |
|---|---|---|---|
| 4 | spectrum reference | bit-identical | nonzero (h_l ≠ h_{N-1-l} per pair) |
| 5 | spectrum reference | bit-identical | nonzero |

Verified by `F93DetuningAntiPalindromicSpectralInvarianceTests.Spectrum_InvariantUnderAntiPalindromicH_*` using the F71-refined diagonal-block spectrum (`F71MirrorBlockRefinement.ComputeSpectrumPerBlock`); the full-L spectrum (via `LiouvillianBlockSpectrum.ComputeSpectrumPerBlock`) does NOT match across anti-palindromic h profiles, consistent with F91/F92's analogous caveats.

## The three-axis family is complete (so far)

F91 + F92 + F93 cover the three parameter axes that:
1. Appear in chain XY+Z-deph+h_l Z_l Hamiltonian/dissipator at the bilinear-coupling level
2. Preserve joint-popcount conservation (no transverse fields, no T1 amplitude damping)
3. Have non-trivial F71 action (each axis transforms under bit-reversal of the site index)

Further parameter axes (e.g. T1 rates if T1 is added per site; ZZ-couplings if Heisenberg chain replaces XY chain) would extend the family. The structural identity is generic; each new axis is a new typed Claim by the same template.

## Anchors

- Memory record: F91 + F92 + F93 close the three-axis family for chain XY+Z-deph+h_l Z_l
- Typed claim: `F93DetuningAntiPalindromicSpectralInvariance` (Tier1Derived, two parents)
- Sister proofs: F91 `docs/proofs/PROOF_F91_GAMMA_NINETY_DEGREES.md`, F92 `docs/proofs/PROOF_F92_BOND_ANTI_PALINDROMIC_J.md`
- Inventory: `docs/SYMMETRY_FAMILY_INVENTORY.md`
- Reflection: `reflections/ON_THE_SYMMETRY_FAMILY.md`
