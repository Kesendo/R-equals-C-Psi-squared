# PROOF F108 Part 2: Π²_X-Even Hamiltonians Always Admit an Exact Palindrome Operator under X-Dephasing

**Status:** Tier 1 derived (closed-form via X-dephasing variant of Π_5bilinear + F1-style algebra; BitA twin of F108 Part 1).
**Klein-V₄ corollary:** Welle 14 (2026-05-27) showed Part 2 also follows from Part 1 by Hilbert-space Hadamard transport (`docs/proofs/PROOF_F108_KLEIN_V4_EQUIVALENCE.md`); the direct proof below is the canonical Π_5b(X) version and is preserved here.
**Date:** 2026-05-25 (direct proof); 2026-05-27 (Klein-V₄ corollary added).
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [F108 Part 1](PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md) (F108 Part 1, BitB-axis sibling under Z-dephasing; this Part 2 mirrors its proof structure exactly)
- [F108 Klein-V₄ equivalence](PROOF_F108_KLEIN_V4_EQUIVALENCE.md) (Welle 14: Part 2 as a Klein-V₄ corollary of Part 1 via Hadamard transport. The proof below is the direct canonical-Π_5b(X) version; the Klein-V₄ corollary establishes existence of a Π_5b-family palindrome operator for L_X via a different representative U_op · Π_5b(Z) · U_op^†.)
- [F85 k-body generalization](PROOF_F85_KBODY_GENERALIZATION.md) (Z-dephasing k-body truly criterion)
- [F107](PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md) (X-dephasing truly criterion: #X even AND #Y even, derived via Π letter-cycle transport from F85)
- [Palindromic Symmetry Beyond Heisenberg](../../experiments/NON_HEISENBERG_PALINDROME.md) (Π-family classification: P1, P4, alternating, continuous per-site, all local)
- [`compute/RCPsiSquared.Core/Symmetry/Pi5BilinearOperator.cs`](../../compute/RCPsiSquared.Core/Symmetry/Pi5BilinearOperator.cs) (Π_5bilinear builder, X-deph variant)

## Abstract

Part 1 closes the palindrome for Π²-even bilinears under Z-dephasing. Part 2 closes the same kind of statement on the orthogonal axis: Π²_X-even bilinears under X-dephasing. The bilinear set changes because X-dephasing classifies by the bit_a parity instead of bit_b; the five bilinears that come out even under X-dephasing's Π² rule are {ZZ, XX, XY, YX, YY} rather than Part 1's {XX, YY, YZ, ZY, ZZ}. Three bilinears (XX, YY, ZZ) are common to both axes; the other two swap.

The proof structure mirrors Part 1 exactly. There is an X-dephasing-specific phase variant of the Π_5bilinear operator, with its own per-site letter permutation (I↔Z, X↔Y, picking out the Π_X letter swap) and its own two-phase-flip pattern. With this operator in hand, the operator-level palindrome identity holds bit-exactly for any X-dephasing-axis Π²-even bilinear Hamiltonian on any sites with any per-site X-dephasing rates.

The structural consequence is symmetric to Part 1: no Π²_X-even Pauli pair can be F87-hard under X-dephasing, because the spectrum is palindromic by construction. F87-hardness on the X-axis is confined to Π²_X-odd or mixed-parity content, exactly mirroring the Z-axis story.

Welle 14 later showed that Part 2 also follows from Part 1 as a Hadamard-transport corollary (the Hilbert-space Hadamard rotates the spin algebra to exchange X and Z, turning a Z-dephasing system into an X-dephasing one). The direct proof here, written by the F1 algebra applied to the X-axis Π_5bilinear variant, is preserved because it is the canonical construction; the Klein-V₄ corollary route is the higher-level perspective.

**Statement (Theorem F108 Part 2):** For any Hamiltonian H built as a linear combination of Π²_X-even 2-site bilinears {ZZ, XX, XY, YX, YY} on N sites with arbitrary real bond coefficients, and X-dephasing on every site with arbitrary per-site rates γ_l, there exists a per-site Liouville-space operator Π_5bilinear (X-deph variant) such that

  Π_5bilinear · L · Π_5bilinear⁻¹ = −L − 2σ·I exactly, where σ = Σ_l γ_l.

In particular, spec(L) is palindromic around −σ, hence no pure-Π²_X-even Pauli pair (truly or non-truly) can be F87-hard under X-dephasing.

This is the BitA-axis twin of F108 Part 1; together they cover the Z- and X-dephasing branches of the F108 Π²-even palindrome family. The Y-dephasing analog is F108 Part 3 ([`PROOF_F108_PART3_PI2Y_EVEN_ALWAYS_PALINDROMIC`](PROOF_F108_PART3_PI2Y_EVEN_ALWAYS_PALINDROMIC.md), Tier 1 derived 2026-05-25 via the Y-deph variant of Π_5bilinear, same I↔X / Y↔Z permutation as Part 1 with Y-deph's −i phase convention), completing the Z/X/Y trio and promoting F109 to fully unconditional Tier1Derived.

## The Π_5bilinear operator (X-dephasing variant)

Per-site Liouville-space automorphism with action on the four Pauli labels:

  I → +1 · Z,    Z → −1 · I,    X → −i · Y,    Y → +i · X.

In the 4×4 label-basis matrix form on {I, X, Y, Z}:

```
        I   X   Y   Z
   I  [ 0   0   0  -1 ]
   X  [ 0   0  +i   0 ]
   Y  [ 0  -i   0   0 ]
   Z  [ 1   0   0   0 ]
```

Same I↔Z, X↔Y permutation as the canonical X-dephasing Π (per `PiOperator` with `PauliLetter.X`), with two phase flips relative to the canonical choice (Z → +I, Y → −iX): the Z→I and Y→X back-arrows carry sign −1 and +i respectively (canonical: +1 and −i).

Key per-site facts:

1. **M is a Liouville-space automorphism, not a Hilbert-space conjugation.** Same subtlety as F108 Part 1's Π_5bilinear.
2. **M² = diag(−1, +1, +1, −1) on {I, X, Y, Z}.** So M⁴ = I and M is order-4. The {I, Z} 2-cycle squares to −1 (the immune-pair under X-dephasing); the {X, Y} 2-cycle squares to +1 (the damped-pair).
3. **Π_5bilinear is unitary on the d²-dim Liouville space.** Each column has one non-zero entry of unit modulus; columns and rows are pairwise orthogonal.

The sign-pattern of M² is structurally the mirror of F108 Part 1's M² = diag(−1, −1, +1, +1) under the Z↔X label swap, matching the bit_a vs bit_b roles in the two dephasing pictures.

## Proof

### Step 1: anti-commutation with every Π²_X-even bilinear

Let Q = M^⊗N be the full N-site Π_5bilinear operator (X-deph variant) on the 4^N-dim Pauli basis. For every Π²_X-even 2-body bilinear B ∈ {ZZ, XX, XY, YX, YY}, the commutator superoperator [B, ·] anti-commutes with Q:

  {Q, [B, ·]} = Q · [B, ·] + [B, ·] · Q = 0.

This is verified bit-exactly (residual = 0 at machine precision) at the 2-qubit level for the 16×16 superoperator. The 4 Π²_X-odd 2-body bilinears {XZ, YZ, ZX, ZY} produce residual = 8.00 (clean separation; Π_5bilinear (X-deph variant) does NOT anti-commute with them, as expected).

The extension from 2 sites to N sites follows the same argument as F108 Part 1 Step 1: B is a 2-body operator, the commutator [B_(l,l+1), ·] acts non-trivially only on the (l, l+1) factor, M acts on each identity factor as a permutation with overall sign ±1 that cancels in the {·, ·} bracket.

**Consequence for the Hamiltonian part of L:** L_H = −i [H, ·]. For H = Σ_b α_b B_b a sum of Π²_X-even bilinears (each B_b in the set above with coefficient α_b ∈ ℝ),

  Q · L_H · Q⁻¹ = Σ_b α_b · (−i) · (−[B_b, ·]) = −L_H.

### Step 2: per-site identity for the X-dephasing dissipator

The Lindblad X-dephasing dissipator on site l is

  D[X_l] · ρ = γ_l · (X_l · ρ · X_l − ρ).

In vec basis: D[X_l] = γ_l · (X_l ⊗ X_l* − I_{d²}). Per site, conjugation by the single-site M satisfies

  M · D[X] · M⁻¹ = −D[X] − 2γ · I_4.

Verified bit-exactly at the 1-qubit level (residual = 0). The mechanism is a diagonal permutation in the Pauli basis, exactly mirroring F108 Part 1 Step 2 under the Z↔X swap. The single-qubit X-dephasing dissipator in the {I, X, Y, Z} Pauli basis is

  D[X]_pauli = γ · diag(0, 0, −2, −2)

(zeros on the {I, X} commuting sector, −2γ on the {Y, Z} anti-commuting sector). M is the per-site signed permutation with permutation (I↔Z, X↔Y) and phases; the conjugation M · D · M⁻¹ for a diagonal D in this basis permutes the diagonal entries by the underlying letter permutation (the phase factors cancel pairwise on each 2-cycle: +1 · −1 on I↔Z, −i · +i on X↔Y). Applying the swap (I↔Z, X↔Y) to diag(0, 0, −2, −2) yields

  M · D[X]_pauli · M⁻¹ = γ · diag(−2, −2, 0, 0) = −D[X]_pauli − 2γ · I_4.

The identity transfers from the Pauli basis to the standard vec basis by the unitary change-of-basis T, since both sides of the identity are unchanged by similarity.

**Consequence for the dissipator part of L:** L_D = Σ_l D[X_l]. M acts as a per-site product Q = M^⊗N, so

  Q · L_D · Q⁻¹ = Σ_l (−D[X_l] − 2γ_l · I_{d²}) = −L_D − 2σ · I_{d²},

where σ = Σ_l γ_l.

### Step 3: combining Hamiltonian and dissipator

  Q · L · Q⁻¹ = Q · L_H · Q⁻¹ + Q · L_D · Q⁻¹ = −L_H − L_D − 2σ · I = −L − 2σ · I.

Bit-exact for every H in the Π²_X-even bilinear family + X-dephasing on every site.

### Step 4: spectral palindrome and F108 Part 2 corollary

From Q · L · Q⁻¹ = −L − 2σ · I and unitarity of Q:

  spec(L) = spec(Q · L · Q⁻¹) = spec(−L − 2σ · I) = {−λ − 2σ : λ ∈ spec(L)}.

So spec(L) is palindromic around −σ.

**F87 corollary:** A Π²_X-even Pauli pair is F87-hard under X-dephasing iff spec(L) breaks palindromy. Since spec(L) is palindromic for every Π²_X-even H (truly or non-truly), no Π²_X-even pair can be F87-hard under X-dephasing. ∎

## Empirical verification

Bit-exact residual ‖Π_5bilinear (X-deph) · L · Π⁻¹ + L + 2σ · I‖_F = 0 at machine precision, across:

| Setup | N range | residual |
|-------|---------|----------|
| All 9 pure-Π²_X-even non-truly pairs (single-bilinear XY/YX + two-term combinations) | N = 3, 4, 5 | 0 |
| 15 random non-uniform-J instances on Π²_X-even bilinear family (5 trials × N ∈ {3, 4, 5}) | N = 3, 4, 5 | 0 |
| Pure D[X]^⊗N dissipator (no Hamiltonian) | N = 1, 3, 4, 5 | 0 |

Reproduction: [`simulations/f108_part2_x_dephasing_scan.py`](../../simulations/f108_part2_x_dephasing_scan.py); C# tests in [`compute/RCPsiSquared.Core.Tests/Symmetry/F108Part2Pi2XEvenAlwaysPalindromicTests.cs`](../../compute/RCPsiSquared.Core.Tests/Symmetry/F108Part2Pi2XEvenAlwaysPalindromicTests.cs).

## Significance

F108 Part 2 completes the BitA twin of F108 Part 1; together they close the F108 Π²-even hardness question across Z- and X-dephasing:

- **F108 Part 1** (BitB axis, Tier 1 derived, 2026-05-25): no Π²_Z-even pair is F87-hard under Z-dephasing.
- **F108 Part 2** (BitA axis, Tier 1 derived, 2026-05-25 THIS PROOF): no Π²_X-even pair is F87-hard under X-dephasing.
- **F109** (Tier 1 derived): mother sector Klein (0, 0) soft ⟹ y_par = 1. After F108 Part 1+2, the Z- and X-dephasing branches of F109 Step 5 are both closed-form; only the Y-dephasing branch remains empirically anchored.

The proof's structural pattern transfers cleanly from F108 Part 1 by the bit_a ↔ bit_b mirror: per-site permutation (I↔X, Y↔Z) under Z-deph maps to (I↔Z, X↔Y) under X-deph; M² sign-pattern diag(−1, −1, +1, +1) maps to diag(−1, +1, +1, −1); D[Z]_pauli diagonal (0, −2, −2, 0) maps to D[X]_pauli diagonal (0, 0, −2, −2). The diagonal-permutation mechanism in Step 2 transfers identically.

## Sibling y_par-axis claims

Closed 2026-05-25: F108 Part 1+2+3 (Π²-even palindrome family, Tier1Derived); F109 (MotherSoftYParityOnePurity, Tier1Derived unconditional); F110 (HardCellYInversionPattern, Tier1Derived since 2026-06-10). Together F107+F109+F110 close the y_par-axis F87 trichotomy classification.

## Open

- Closed-form derivation of F110 Aspect C exact ratios (42:8 at k=3, 228:0 at k=4) per Pauli-letter combinatorics. F103 Section 5 explicitly lists as open.
- k ≥ 5 empirical confirmation of F103/F106 pattern stability beyond N=4.
- Hardware QPU confirmation at k ≥ 3 (no F87 QPU confirmations exist beyond Marrakesh k=2).

∎
