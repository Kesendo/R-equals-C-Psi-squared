# PROOF F108 Part 2: Π²_X-Even Hamiltonians Always Admit an Exact Palindrome Operator under X-Dephasing

**Status:** Tier 1 derived (closed-form via X-dephasing variant of Π_5bilinear + F1-style algebra; BitA twin of F108 Part 1)
**Date:** 2026-05-25
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md](PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md) (F108 Part 1, BitB-axis sibling under Z-dephasing; this Part 2 mirrors its proof structure exactly)
- [PROOF_F85_KBODY_GENERALIZATION.md](PROOF_F85_KBODY_GENERALIZATION.md) (truly criterion for k-body Pauli terms; X-dephase truly = #X even AND #Y even)
- [NON_HEISENBERG_PALINDROME](../../experiments/NON_HEISENBERG_PALINDROME.md) (Π-family classification: P1, P4, alternating, non-local)
- [`compute/RCPsiSquared.Core/Symmetry/Pi5BilinearOperator.cs`](../../compute/RCPsiSquared.Core/Symmetry/Pi5BilinearOperator.cs) (Π_5bilinear builder, X-deph variant)

**Statement (Theorem F108 Part 2):** For any Hamiltonian H built as a linear combination of Π²_X-even 2-site bilinears {ZZ, XX, XY, YX, YY} on N sites with arbitrary real bond coefficients, and X-dephasing on every site with arbitrary per-site rates γ_l, there exists a per-site Liouville-space operator Π_5bilinear (X-deph variant) such that

  Π_5bilinear · L · Π_5bilinear⁻¹ = −L − 2σ·I exactly, where σ = Σ_l γ_l.

In particular, spec(L) is palindromic around −σ, hence no pure-Π²_X-even Pauli pair (truly or non-truly) can be F87-hard under X-dephasing.

This is the BitA-axis twin of F108 Part 1; together they cover the Z- and X-dephasing branches of the F108 Π²-even palindrome family. The Y-dephasing analog (F108 Part 3) is open.

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
| 5 random non-uniform-J instances on Π²_X-even bilinear family | N = 3, 4, 5 | 0 |
| Pure D[X]^⊗N dissipator (no Hamiltonian) | N = 1, 3, 4, 5 | 0 |

Reproduction: [`simulations/_f108_part2_x_dephasing_scan.py`](../../simulations/_f108_part2_x_dephasing_scan.py); C# tests in [`compute/RCPsiSquared.Core.Tests/Symmetry/F108Part2Pi2XEvenAlwaysPalindromicTests.cs`](../../compute/RCPsiSquared.Core.Tests/Symmetry/F108Part2Pi2XEvenAlwaysPalindromicTests.cs).

## Significance

F108 Part 2 completes the BitA twin of F108 Part 1; together they close the F108 Π²-even hardness question across Z- and X-dephasing:

- **F108 Part 1** (BitB axis, Tier 1 derived, 2026-05-25): no Π²_Z-even pair is F87-hard under Z-dephasing.
- **F108 Part 2** (BitA axis, Tier 1 derived, 2026-05-25 THIS PROOF): no Π²_X-even pair is F87-hard under X-dephasing.
- **F109** (Tier 1 derived): mother sector Klein (0, 0) soft ⟹ y_par = 1. After F108 Part 1+2, the Z- and X-dephasing branches of F109 Step 5 are both closed-form; only the Y-dephasing branch remains empirically anchored.

The proof's structural pattern transfers cleanly from F108 Part 1 by the bit_a ↔ bit_b mirror: per-site permutation (I↔X, Y↔Z) under Z-deph maps to (I↔Z, X↔Y) under X-deph; M² sign-pattern diag(−1, −1, +1, +1) maps to diag(−1, +1, +1, −1); D[Z]_pauli diagonal (0, −2, −2, 0) maps to D[X]_pauli diagonal (0, 0, −2, −2). The diagonal-permutation mechanism in Step 2 transfers identically.

## Open

- **F108 Part 3** (Y-dephasing analog): no covering Claim yet. Would derive the Π_5bilinear variant for Y-dephasing on the Π²_Y-even bilinear set + D[Y] per-site identity. The structural pattern should transfer the same way (per-site permutation under Y-deph is (I↔X, Y↔Z) per `PiOperator`, same as Z-deph but with different phase factors). When F108 Part 3 lands, F109 will be closed-form across all three dephase letters.
- **F110** (hard cells y_par-pure with Y-inversion): higher difficulty, per-dephase-letter algebra on the F87-hard pair-set.

∎
