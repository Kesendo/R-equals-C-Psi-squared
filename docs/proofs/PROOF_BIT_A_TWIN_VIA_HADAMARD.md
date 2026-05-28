# PROOF: bit_a Twins via the X to Z Hadamard Duality (consolidation theorem)

**Status:** Tier1Derived consolidation. Collects the per-formula bit_a transports already proven in Welle 12 to 15 into one general classification theorem for the BitATwin slots. No new heavy machinery; it names the precondition under which a BitB claim automatically owns its bit_a twin, and the exception classes where it does not.
**Date:** 2026-05-28
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Builds on:**
- [PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md](PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md): the operator-space Klein-V₄ element `Q_zx` with `Q_zx · Π_Z · Q_zx⁻¹ = Π_X` (universal N, bit-exact).
- [PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md](PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md): the worked Lindblad transport ("F112 is one identity, not three") and the lift caveat.
- [PROOF_F108_KLEIN_V4_EQUIVALENCE.md](PROOF_F108_KLEIN_V4_EQUIVALENCE.md): the cautionary bespoke-operator case (Π_5bilinear does not transport via the operator-space element alone).
- [PROOF_ABSORPTION_THEOREM.md](PROOF_ABSORPTION_THEOREM.md): Re(λ) = −2γ⟨popcount⟩, H-independent for any Hermitian H (the dynamics enabler).

## Abstract

The bit_a-twin slots of the cubic Z₂³ polarity architecture (`BitATwinClassification`) had 52 BitB claims marked NeedsDerivation. They are not 52 independent problems. The bit_a axis (Π²_X = Z⊗N) is the X to Z image of the bit_b axis (Π²_Z = X⊗N), and that image is realized by a proven duality. This note states the general theorem: a BitB claim whose content reduces to the spectrum, eigenspaces, or operator-identities of Π / L (or to the Absorption-Theorem popcount reading) automatically owns its bit_a twin by the X to Z duality, and is classified `CoveredByHadamardDuality`. Two exception classes keep `NeedsDerivation` (bespoke operators) or `BitBSpecific` (no twin exists).

## Statement

Let C be a BitB `IZ2AxisClaim` (a statement about an N-qubit system under site-local Z-dephasing and any Hermitian Hamiltonian H, classified on the Π²_Z = X⊗N axis). Let C' be its bit_a image: the same statement with Π_Z replaced by Π_X and Z-dephasing replaced by X-dephasing.

**Theorem.** C' holds, with identical structure, whenever C's content is one of:

1. **operator-space** (the spectrum, eigenvalues, eigenspace dimensions, or an operator identity of Π_Z, or a Frobenius-norm statement on the 4^N Pauli space), or
2. **Lindblad-spectral** (a statement about Re(λ) of the Liouvillian L, or about the Z-dephasing dissipator and its popcount reading).

In case 1, C' is the conjugate of C by the operator-space Klein element Q_zx (`Q_zx · Π_Z · Q_zx⁻¹ = Π_X`, PROOF_KLEIN_V4). In case 2, C' is the conjugate of C by the genuine Hilbert-space global Hadamard `U_H = H^⊗N` (`H X H = Z` per site), which maps Z-dephasing to X-dephasing.

## Proof

**Case 1 (operator-space).** Q_zx is an order-2 element of the Klein-V₄ subgroup acting on the 4^N Pauli operator space, with `Q_zx · Π_Z · Q_zx⁻¹ = Π_X`, proven bit-exact at universal N in PROOF_KLEIN_V4 (per-site 4×4 identity tensor-powered by the Kronecker mixed-product property). Conjugation by Q_zx is a similarity transform, so it preserves spectra, eigenvalue multiplicities, eigenspace dimensions, and operator identities verbatim. Any C of operator-space type therefore transports to C' by replacing Π_Z with its Q_zx-conjugate Π_X. The F38 / F39 / F63 operator-identity triple are the already-typed instances (F38BitA: Π²_X = (−1)^{n_XY}, half-half eigenspace split identical to F38 by this argument).

**Case 2 (Lindblad-spectral).** By the Absorption Theorem, Re(λ_k) = ⟨v_k|Herm(L)|v_k⟩ depends only on the dissipator: Herm(L) is the pure Z-dephasing dissipator, diagonal with entries −2γ·popcount(i XOR j), and the Hamiltonian part is anti-Hermitian for any Hermitian H, so it contributes nothing to Re(λ). The global Hadamard U_H = H^⊗N is a Hilbert-space unitary that maps each jump operator Z_l to X_l, hence Z-dephasing to X-dephasing, and popcount in the Z-basis to popcount in the X-basis. Conjugating L by U_H therefore carries the bit_b spectral statement to the bit_a one. F112-X is the worked instance ("F112 is one identity, not three", PROOF_F112).

## Preconditions and exception classes

- **Collapsible (status `CoveredByHadamardDuality`):** C reduces to operator-space or Lindblad-spectral content as above. This covers the Absorption-descendant family (F33, F50, F55, F64, F65, F66, F67, F68, F74, F89), the operator-identity and parity claims, and the spectral / mirror claims. Each is a corollary of this theorem; no bespoke per-claim derivation is owed.

- **Lift caveat:** of the Klein-V₄ elements, only {I, U_H (Hadamard)} lifts to a Hilbert-space unitary that preserves Lindblad form. D (Z to Y) and Q_yx are operator-space-only and break Lindblad form (PROOF_F112). So a Lindblad-dynamics claim transports only along the Z to X leg, not to the Y axis. This bounds case 2 to the X to Z duality.

- **Bespoke-operator residue (keeps `NeedsDerivation`):** claims built on operators other than the canonical Π_d, e.g. Π_5bilinear (F108 family), are not transported by the operator-space Q_zx alone; they need a deeper Hilbert-space Hadamard at the spin-algebra level (PROOF_F108). These are flagged for individual attention, not covered here.

- **bit_b-specific (status `BitBSpecific`):** claims whose content is intrinsically tied to Z-dephasing or to non-dephasing dissipation (amplitude damping F1T1 / F82 / F84, and the Y-dephase / break-magnitude claims F91 / F93 / F108Part3 / F113 / F112-Y) have no bit_a twin to mirror. F61's break conditions show that the σ⁺/σ⁻ jumps of amplitude damping flip bit_a, so the symmetry does not exist.

## Consequence for the map

The 52 NeedsDerivation slots split into the collapsible majority (this theorem, `CoveredByHadamardDuality`) and the small bespoke-operator residue (genuine open work). Together with the 7 already-typed Filled twins and the 8 BitBSpecific slots, the bit_a side of the cubic architecture is closed in principle: one duality theorem plus a short residue list, not 52 independent derivations.

## Anchors

- [PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md](PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md), [PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md](PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md), [PROOF_F108_KLEIN_V4_EQUIVALENCE.md](PROOF_F108_KLEIN_V4_EQUIVALENCE.md), [PROOF_ABSORPTION_THEOREM.md](PROOF_ABSORPTION_THEOREM.md).
- Typed: `BitATwinClassification.CoveredByHadamardDuality` (the status this theorem licenses on the collapsible BitB claims; counted by `PolarityCubeMap`).
- Map: [docs/PI2KB_INHERITANCE_MAP.md](../PI2KB_INHERITANCE_MAP.md) (Z2Axis classification section).
- Already-typed instances: `F38BitAInvolutionInheritance.cs`, `F39DetPiBitAInheritance.cs`, `F63BitAReference.cs`, `F108Part2Pi2XEvenAlwaysPalindromic.cs`, `LindbladBitAPiBalance.cs`.
