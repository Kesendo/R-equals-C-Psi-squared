# PROOF F108 Part 3: Π²_Y-Even Bilinears Always Admit an Exact Palindrome Operator under Y-Dephasing

**Status:** Tier 1 derived (closed-form via Y-dephasing variant of Π_5bilinear + F1-style algebra; Y-dephasing sibling of F108 Part 1).
**Corollary of Part 1:** Part 3 also follows from Part 1, either by D on the mirror (D · Π_5b(Z) · D = Π_5b(Y), the bilinear set D-invariant, the dissipator identity re-checked for Y) or by the quarter turn about X, which carries the Lindbladian and fixes Π_5b(Z) = (−1)^N · Π_5b(Y)⁻¹; see `docs/proofs/PROOF_F108_KLEIN_V4_EQUIVALENCE.md`. The direct proof below is the canonical Part 3 derivation.
**Date:** 2026-05-25 (direct proof); 2026-05-27 (corollary of Part 1).
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [F108 Part 1](PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md) (F108 Part 1, Z-dephasing on the same BitB axis; this Part 3 mirrors its proof structure with the Y-dephase-appropriate phase choice)
- [F108 Klein-V₄ equivalence](PROOF_F108_KLEIN_V4_EQUIVALENCE.md) (Part 3 from Part 1 by D on the mirror, the bilinear set being D-invariant and the dissipator identity re-checked for Y, or by the quarter turn about X on the Lindbladian.)
- [F85 k-body generalization](PROOF_F85_KBODY_GENERALIZATION.md) (Z-dephasing k-body truly criterion)
- [F107](PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md) (Y-dephasing truly criterion: #Y even AND #Z even, identical to Z-deph since Π_Y = Π_Z⁻¹)
- [Palindromic Symmetry Beyond Heisenberg](../../experiments/NON_HEISENBERG_PALINDROME.md) (Π-family taxonomy; F108 Part 3 sits in the same P1 family as Part 1, with the Y/Z 2-cycle phase variant)
- [`compute/RCPsiSquared.Core/Symmetry/Pi5BilinearOperator.cs`](../../compute/RCPsiSquared.Core/Symmetry/Pi5BilinearOperator.cs) (Π_5bilinear builder, Y-deph variant)

## Abstract

Part 3 closes the third axis of the F108 trinity: Π²_Y-even bilinears under Y-dephasing. The bilinear set turns out to be identical to Part 1's, because Y and Z share the bit_b axis under the Π² eigenvalue rule. The same five bilinears {XX, YY, YZ, ZY, ZZ} come out even on both axes. The difference is in the conjugating operator: the Y-dephasing variant of Π_5bilinear has the same letter permutation as the Z-dephasing version (I↔X, Y↔Z) but flips the Y/Z 2-cycle phase from +i to −i.

The proof carries over from Part 1 once the phase variant is identified. The anti-commutation argument with the Π²-even bilinears works the same way; the dissipator-side identity holds because the letter swap exchanges Y-dephasing's undamped pair {I, Y} with its damped pair {X, Z}. The operator-level palindrome holds bit-exactly for any Y-dephasing-axis Π²-even bilinear Hamiltonian on any sites with any per-site Y-dephasing rates.

The structural consequence is the third leg of the F87-hardness collapse: no pair of Π²_Y-even bilinears can be F87-hard under Y-dephasing. Together with Parts 1 and 2, the F108 trinity rules out hardness on Π²-even bilinears across all three single-letter dephase channels.

Part 3 also follows from Part 1 ([F108 Klein-V₄ equivalence](PROOF_F108_KLEIN_V4_EQUIVALENCE.md)): either D, which maps Π_5bilinear(Z) to Π_5bilinear(Y) and fixes the bilinear set, the dissipator identity being checked for Y separately, or the quarter turn about X, which carries the Part-1 Lindbladian to the Part-3 one and fixes Π_5bilinear(Z). The direct proof is preserved as the canonical Part 3 construction.

**Statement (Theorem F108 Part 3):** For any Hamiltonian H built as a linear combination of Π²_Y-even 2-site bilinears {XX, YY, YZ, ZY, ZZ} on N sites with arbitrary real bond coefficients, and Y-dephasing on every site with arbitrary per-site rates γ_l, there exists a per-site Liouville-space operator Π_5bilinear (Y-deph variant) such that

  Π_5bilinear · L · Π_5bilinear⁻¹ = −L − 2σ·I exactly, where σ = Σ_l γ_l.

In particular, spec(L) is palindromic around −σ, hence no pair of these bilinears can be F87-hard under Y-dephasing.

This is the Y-dephasing sibling of F108 Part 1 (same BitB axis, same bilinear set, different dephase letter). Together with Part 1 (Z-deph) and Part 2 (X-deph, BitA twin) the three parts cover the F108 Π²-even palindrome family across all three single-letter dephase channels.

## Π²_Y-even bilinear set is identical to Π²_Z-even

Per `PiOperator.SquaredEigenvalue`, Π²_Y and Π²_Z both compute (−1)^bit_b (the Y and Z dephase letters share the bit_b classification axis; only X-dephasing flips to bit_a). The Π²_Y-even (bit_b = 0) 2-site bilinear set therefore equals the Π²_Z-even set: {XX, YY, YZ, ZY, ZZ}. The truly criterion under Y-dephasing is #Y even AND #Z even, identical to the Z-deph criterion (F107).

## The Π_5bilinear operator (Y-dephasing variant)

Per-site Liouville-space automorphism with action on the four Pauli labels:

  I → +1 · X,    X → −1 · I,    Y → −i · Z,    Z → +i · Y.

Same I↔X, Y↔Z permutation as the Z-dephasing variant (since Y-deph shares Π²_Z's permutation structure per `PiOperator`'s Y branch); the only difference is the Y/Z 2-cycle phase: −i (Y-deph) versus +i (Z-deph). The I↔X 2-cycle phases (I → +1, X → −1) are identical to the Z-deph variant.

Key per-site facts:

1. **M is a Liouville-space automorphism, not a Hilbert-space conjugation.** Same subtlety as F108 Parts 1 and 2.
2. **M² = diag(−1, −1, +1, +1) on {I, X, Y, Z}.** Identical sign pattern to F108 Part 1's M², because the I↔X 2-cycle product is unchanged and the Y↔Z 2-cycle product (−i)·(+i) = +1 matches the Z-deph variant's (+i)·(−i) = +1. So M⁴ = I and M is order-4.
3. **Π_5bilinear (Y-deph) is unitary on the d²-dim Liouville space.** Signed permutation matrix; each column has one non-zero entry of unit modulus.

## Proof

### Step 1: anti-commutation with every Π²_Y-even bilinear

Let Q = M^⊗N be the full N-site Π_5bilinear (Y-deph variant). For every Π²_Y-even 2-body bilinear B ∈ {XX, YY, YZ, ZY, ZZ}, the commutator superoperator [B, ·] anti-commutes with Q:

  {Q, [B, ·]} = Q · [B, ·] + [B, ·] · Q = 0.

Verified bit-exactly (residual = 0 at machine precision) at the 2-qubit level. The 4 Π²_Y-odd 2-body bilinears {XY, XZ, YX, ZX} produce residual = 8.00 (clean separation).

Why this lifts from Part 1: the per-site M (Y-deph variant) differs from the Z-deph variant only in the Y/Z 2-cycle phase. The anti-commutation calculation for Π²-even bilinears (which contain even numbers of Y and Z together) is sign-pattern-driven, not phase-magnitude-driven. The flipped phase on the Y/Z 2-cycle leaves the {Q, [B, ·]} = 0 algebra intact for the bit_b-even bilinear set. The 2-qubit anti-commutation result transfers to N qubits exactly as in Part 1.

**Consequence for the Hamiltonian part of L:** L_H = −i [H, ·]. For H = Σ_b α_b B_b a sum of Π²_Y-even bilinears with α_b ∈ ℝ,

  Q · L_H · Q⁻¹ = −L_H.

### Step 2: per-site identity for the Y-dephasing dissipator

The Lindblad Y-dephasing dissipator on site l is

  D[Y_l] · ρ = γ_l · (Y_l · ρ · Y_l − ρ).

In vec basis: D[Y_l] = γ_l · (Y_l ⊗ Y_l* − I_{d²}). Per site, conjugation by the single-site M satisfies

  M · D[Y] · M⁻¹ = −D[Y] − 2γ · I_4.

Verified bit-exactly at the 1-qubit level (residual = 0). The mechanism is the diagonal-permutation argument from F108 Parts 1 and 2. In the {I, X, Y, Z} Pauli basis the Y-dephasing dissipator is

  D[Y]_pauli = γ · diag(0, −2, 0, −2)

(zeros on the {I, Y} commuting sector; −2γ on the {X, Z} anti-commuting sector). M's (I↔X, Y↔Z) per-site swap permutes the diagonal entries by 2-cycle (each phase of M meets its own inverse in M⁻¹):

  M · D[Y]_pauli · M⁻¹ = γ · diag(−2, 0, −2, 0) = −D[Y]_pauli − 2γ · I_4.

The identity transfers from the Pauli basis to the standard vec basis by the unitary change-of-basis T.

**Consequence for the dissipator part of L:** L_D = Σ_l D[Y_l]; M acts as a per-site product Q = M^⊗N, so Q · L_D · Q⁻¹ = −L_D − 2σ · I_{d²}.

### Step 3: combining Hamiltonian and dissipator

  Q · L · Q⁻¹ = Q · L_H · Q⁻¹ + Q · L_D · Q⁻¹ = −L_H − L_D − 2σ · I = −L − 2σ · I.

Bit-exact for every H in the Π²_Y-even bilinear family + Y-dephasing on every site.

### Step 4: spectral palindrome and F108 Part 3 corollary

  spec(L) = spec(Q · L · Q⁻¹) = spec(−L − 2σ · I) = {−λ − 2σ : λ ∈ spec(L)},

so spec(L) is palindromic around −σ.

**F87 corollary:** A pair is F87-hard under Y-dephasing iff spec(L) breaks palindromy. Since spec(L) is palindromic for every H of the family (truly or non-truly), no pair of Π²_Y-even bilinears can be F87-hard under Y-dephasing. ∎

## Empirical verification

Bit-exact residual ‖Π_5bilinear (Y-deph) · L · Π⁻¹ + L + 2σ · I‖_F = 0 at machine precision, across:

| Setup | N range | residual |
|-------|---------|----------|
| All 9 pure-Π²_Y-even non-truly pairs (YZ, ZY, XX+YZ, XX+ZY, YY+YZ, YY+ZY, YZ+ZY, YZ+ZZ, ZY+ZZ) | N = 3, 4, 5 | 0 |
| 15 random non-uniform-J instances on Π²_Y-even bilinear family (5 trials × N ∈ {3, 4, 5}) | N = 3, 4, 5 | 0 |
| Pure D[Y]^⊗N dissipator (no Hamiltonian) | N = 1, 3, 4, 5 | 0 |

Reproduction: [`simulations/f108_part3_y_dephasing_scan.py`](../../simulations/f108_part3_y_dephasing_scan.py); C# tests in [`compute/RCPsiSquared.Core.Tests/Symmetry/F108Part3Pi2YEvenAlwaysPalindromicTests.cs`](../../compute/RCPsiSquared.Core.Tests/Symmetry/F108Part3Pi2YEvenAlwaysPalindromicTests.cs).

## Significance

F108 Part 3 completes the F108 Π²-even palindrome family across all three single-letter dephase channels:

- **F108 Part 1** (BitB axis, Z-deph, Tier 1 derived, 2026-05-25)
- **F108 Part 2** (BitA axis, X-deph, Tier 1 derived, 2026-05-25)
- **F108 Part 3** (BitB axis, Y-deph, Tier 1 derived, 2026-05-25 THIS PROOF)

Each of the three variants reaches exactly the Π²-D-even strings of even weight (an even number of non-identity letters), D its own dephase letter; the Y-deph variant is Π_Y ∘ Ad_{Y^⊗N} ([Part 1](PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md), Significance). The mother cell's non-truly strings, which [F109](PROOF_F109_MOTHER_SOFT_Y_PARITY_ONE_PURITY.md) is about, have odd weight; their softness is the colouring's (F109 Step 5).

The structural pattern transfers cleanly across the three parts via the per-site permutation + phase-flip algebra:
- Part 1 (Z-deph) and Part 3 (Y-deph) share the I↔X, Y↔Z permutation (both bit_b axis); they differ only in Y/Z 2-cycle phase (+i vs −i) to match each dephase letter's canonical Π phase.
- Part 2 (X-deph) uses the I↔Z, X↔Y permutation (bit_a axis), with analogous back-arrow phase flips.
- The dissipator-side proof (diagonal permutation in the Pauli basis) transfers identically across all three parts.

## Sibling y_par-axis claims

F107 (TrulyYParityZeroPurity), F109 (MotherSoftYParityOnePurity), F110 (HardCellYInversionPattern) and F111 (HardCellPureDTemplate), all Tier1Derived, sit on the y_par axis; F108 sits on the BitB axis (Parts 1 and 3) and the BitA axis (Part 2).

## Open

- k ≥ 5 empirical confirmation of F103/F106 pattern stability beyond N=4, for the soft, truly and off-diagonal cells (the hard diagonal cell is read at k = 4, 5 by `DiagonalCellComplementConnectionTests`).
- Hardware QPU confirmation at k ≥ 3 (no F87 QPU confirmations exist beyond Marrakesh k=2).

∎
