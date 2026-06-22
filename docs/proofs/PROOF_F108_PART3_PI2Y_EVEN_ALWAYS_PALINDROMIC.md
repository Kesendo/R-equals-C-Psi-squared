# PROOF F108 Part 3: Π²_Y-Even Hamiltonians Always Admit an Exact Palindrome Operator under Y-Dephasing

**Status:** Tier 1 derived (closed-form via Y-dephasing variant of Π_5bilinear + F1-style algebra; Y-dephasing sibling of F108 Part 1).
**Klein-V₄ corollary:** Welle 14 (2026-05-27) showed Part 3 follows from Part 1 by operator-space D-conjugation (D · Π_5b(Z) · D = Π_5b(Y) bit-exact); see `docs/proofs/PROOF_F108_KLEIN_V4_EQUIVALENCE.md`. The direct proof below is the canonical Part 3 derivation; the Klein-V₄ corollary makes the equivariance explicit.
**Date:** 2026-05-25 (direct proof); 2026-05-27 (Klein-V₄ corollary added).
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md](PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md) (F108 Part 1, Z-dephasing on the same BitB axis; this Part 3 mirrors its proof structure with the Y-dephase-appropriate phase choice)
- [PROOF_F108_KLEIN_V4_EQUIVALENCE.md](PROOF_F108_KLEIN_V4_EQUIVALENCE.md) (Welle 14: Part 3 as a direct Klein-V₄ corollary of Part 1 via the D operator. The bilinear set is D-invariant; D · Π_5b(Z) · D = Π_5b(Y) bit-exact, so the Part 1 proof transfers verbatim with d=Z→Y substitution.)
- [PROOF_F85_KBODY_GENERALIZATION.md](PROOF_F85_KBODY_GENERALIZATION.md) (Z-dephasing k-body truly criterion)
- [PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md](PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md) (Y-dephasing truly criterion: #Y even AND #Z even, identical to Z-deph since the Y↔Z 2-cycle fixes the {#Y, #Z} pair setwise)
- [NON_HEISENBERG_PALINDROME](../../experiments/NON_HEISENBERG_PALINDROME.md) (Π-family taxonomy; F108 Part 3 sits in the same P1 family as Part 1, with the Y/Z 2-cycle phase variant)
- [`compute/RCPsiSquared.Core/Symmetry/Pi5BilinearOperator.cs`](../../compute/RCPsiSquared.Core/Symmetry/Pi5BilinearOperator.cs) (Π_5bilinear builder, Y-deph variant)

## Abstract

Part 3 closes the third axis of the F108 trinity: Π²_Y-even bilinears under Y-dephasing. The bilinear set turns out to be identical to Part 1's, because Y and Z share the bit_b axis under the Π² eigenvalue rule. The same five bilinears {XX, YY, YZ, ZY, ZZ} come out even on both axes. The difference is in the conjugating operator: the Y-dephasing variant of Π_5bilinear has the same letter permutation as the Z-dephasing version (I↔X, Y↔Z) but flips the Y/Z 2-cycle phase from +i to −i.

The proof carries over from Part 1 once the phase variant is identified. The anti-commutation argument with the Π²-even bilinears works the same way; the dissipator-side identity works because Y-dephasing commutes with itself and the phase change does not affect the dissipator's per-site structure. The operator-level palindrome holds bit-exactly for any Y-dephasing-axis Π²-even bilinear Hamiltonian on any sites with any per-site Y-dephasing rates.

The structural consequence is the third leg of the F87-hardness collapse: no Π²_Y-even Pauli pair can be F87-hard under Y-dephasing. Together with Parts 1 and 2, the F108 trinity rules out hardness on Π²-even bilinears across all three single-letter dephase channels.

Welle 14 later showed that Part 3 also follows from Part 1 as a D-conjugation corollary on operator space (the Welle 12 Task 1 diagonal involution D maps Π_5bilinear(Z) to Π_5bilinear(Y) bit-exactly). The direct proof is preserved as the canonical Part 3 construction; the D-conjugation route makes the equivariance with Part 1 explicit, and lets the F108 trinity be read as a single result with three Klein-V₄ images rather than three independent proofs.

**Statement (Theorem F108 Part 3):** For any Hamiltonian H built as a linear combination of Π²_Y-even 2-site bilinears {XX, YY, YZ, ZY, ZZ} on N sites with arbitrary real bond coefficients, and Y-dephasing on every site with arbitrary per-site rates γ_l, there exists a per-site Liouville-space operator Π_5bilinear (Y-deph variant) such that

  Π_5bilinear · L · Π_5bilinear⁻¹ = −L − 2σ·I exactly, where σ = Σ_l γ_l.

In particular, spec(L) is palindromic around −σ, hence no pure-Π²_Y-even Pauli pair can be F87-hard under Y-dephasing.

This is the Y-dephasing sibling of F108 Part 1 (same BitB axis, same bilinear set, different dephase letter). Together with Part 1 (Z-deph) and Part 2 (X-deph, BitA twin) the three parts cover the F108 Π²-even palindrome family across all three single-letter dephase channels. F109's Step 5 is now closed-form across the full {Z, X, Y} dephasing set, promoting F109 to fully unconditional Tier1Derived.

## Π²_Y-even bilinear set is identical to Π²_Z-even

Per `PiOperator.SquaredEigenvalue`, Π²_Y and Π²_Z both compute (−1)^bit_b (the Y and Z dephase letters share the bit_b classification axis; only X-dephasing flips to bit_a). The Π²_Y-even (bit_b = 0) 2-site bilinear set therefore equals the Π²_Z-even set: {XX, YY, YZ, ZY, ZZ}. F85's truly criterion under Y-dephasing is #Y even AND #Z even, identical to the Z-deph criterion per F107.

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

(zeros on the {I, Y} commuting sector; −2γ on the {X, Z} anti-commuting sector). M's (I↔X, Y↔Z) per-site swap permutes the diagonal entries by 2-cycle (phase factors cancel pairwise: +1·−1 on I↔X, −i·+i on Y↔Z):

  M · D[Y]_pauli · M⁻¹ = γ · diag(−2, 0, −2, 0) = −D[Y]_pauli − 2γ · I_4.

The identity transfers from the Pauli basis to the standard vec basis by the unitary change-of-basis T.

**Consequence for the dissipator part of L:** L_D = Σ_l D[Y_l]; M acts as a per-site product Q = M^⊗N, so Q · L_D · Q⁻¹ = −L_D − 2σ · I_{d²}.

### Step 3: combining Hamiltonian and dissipator

  Q · L · Q⁻¹ = Q · L_H · Q⁻¹ + Q · L_D · Q⁻¹ = −L_H − L_D − 2σ · I = −L − 2σ · I.

Bit-exact for every Π²_Y-even H + Y-dephasing on every site.

### Step 4: spectral palindrome and F108 Part 3 corollary

  spec(L) = spec(Q · L · Q⁻¹) = spec(−L − 2σ · I) = {−λ − 2σ : λ ∈ spec(L)},

so spec(L) is palindromic around −σ.

**F87 corollary:** A Π²_Y-even Pauli pair is F87-hard under Y-dephasing iff spec(L) breaks palindromy. Since spec(L) is palindromic for every Π²_Y-even H (truly or non-truly), no Π²_Y-even pair can be F87-hard under Y-dephasing. ∎

## Empirical verification

Bit-exact residual ‖Π_5bilinear (Y-deph) · L · Π⁻¹ + L + 2σ · I‖_F = 0 at machine precision, across:

| Setup | N range | residual |
|-------|---------|----------|
| All 9 pure-Π²_Y-even non-truly pairs (YZ, ZY, XX+YZ, XX+ZY, YY+YZ, YY+ZY, YZ+ZY, YZ+ZZ, ZY+ZZ) | N = 3, 4, 5 | 0 |
| 15 random non-uniform-J instances on Π²_Y-even bilinear family (5 trials × N ∈ {3, 4, 5}) | N = 3, 4, 5 | 0 |
| Pure D[Y]^⊗N dissipator (no Hamiltonian) | N = 1, 3, 4, 5 | 0 |

Reproduction: [`simulations/_f108_part3_y_dephasing_scan.py`](../../simulations/_f108_part3_y_dephasing_scan.py); C# tests in [`compute/RCPsiSquared.Core.Tests/Symmetry/F108Part3Pi2YEvenAlwaysPalindromicTests.cs`](../../compute/RCPsiSquared.Core.Tests/Symmetry/F108Part3Pi2YEvenAlwaysPalindromicTests.cs).

## Significance

F108 Part 3 completes the F108 Π²-even palindrome family across all three single-letter dephase channels:

- **F108 Part 1** (BitB axis, Z-deph, Tier 1 derived, 2026-05-25)
- **F108 Part 2** (BitA axis, X-deph, Tier 1 derived, 2026-05-25)
- **F108 Part 3** (BitB axis, Y-deph, Tier 1 derived, 2026-05-25 THIS PROOF)

**F109 promotion:** F109's Step 5 (Π²-even non-truly ⟹ soft) is now closed-form across all three dephase letters. F109 promotes from "Tier1Derived for Z- and X-dephasing, Y empirically anchored" to **fully unconditional Tier1Derived**.

The structural pattern transfers cleanly across the three parts via the per-site permutation + phase-flip algebra:
- Part 1 (Z-deph) and Part 3 (Y-deph) share the I↔X, Y↔Z permutation (both bit_b axis); they differ only in Y/Z 2-cycle phase (+i vs −i) to match each dephase letter's canonical Π phase.
- Part 2 (X-deph) uses the I↔Z, X↔Y permutation (bit_a axis), with analogous back-arrow phase flips.
- The dissipator-side proof (diagonal permutation in the Pauli basis) transfers identically across all three parts.

## Sibling y_par-axis claims

Closed 2026-05-25: F108 Part 1+2+3 (Π²-even palindrome family, Tier1Derived); F109 (MotherSoftYParityOnePurity, Tier1Derived unconditional); F110 (HardCellYInversionPattern, Tier1Derived since 2026-06-10). Together F107+F109+F110 close the y_par-axis F87 trichotomy classification.

## Open

- Closed-form derivation of F110 Aspect C exact ratios (42:8 at k=3, 228:0 at k=4) per Pauli-letter combinatorics. F103 Section 5 explicitly lists as open.
- k ≥ 5 empirical confirmation of F103/F106 pattern stability beyond N=4.
- Hardware QPU confirmation at k ≥ 3 (no F87 QPU confirmations exist beyond Marrakesh k=2).

∎
