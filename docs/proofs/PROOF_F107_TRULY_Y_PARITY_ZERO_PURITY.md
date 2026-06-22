# PROOF F107: F87 Truly Classification Forces y_par = 0 (All Dephase Letters)

**Status:** Tier 1 derived (closed-form corollary of F85's Z-dephasing truly criterion + Π letter-cycle transport to X- and Y-dephasing via F108)
**Date:** 2026-05-24
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [PROOF_F85_KBODY_GENERALIZATION.md](PROOF_F85_KBODY_GENERALIZATION.md) (k-body truly criterion under Z-dephasing: #Y even AND #Z even, Step 3)
- [PROOF_F108_PART2_PI2X_EVEN_ALWAYS_PALINDROMIC.md](PROOF_F108_PART2_PI2X_EVEN_ALWAYS_PALINDROMIC.md) (X-dephasing Π_X palindrome closure)
- [PROOF_F108_PART3_PI2Y_EVEN_ALWAYS_PALINDROMIC.md](PROOF_F108_PART3_PI2Y_EVEN_ALWAYS_PALINDROMIC.md) (Y-dephasing Π_Y palindrome closure)
- [MIRROR_SYMMETRY_PROOF.md](MIRROR_SYMMETRY_PROOF.md) (Π definitions per dephase letter)
- [`compute/RCPsiSquared.Core/Symmetry/PiOperator.cs`](../../compute/RCPsiSquared.Core/Symmetry/PiOperator.cs) (per-dephase Π² eigenvalue rules)
- [`compute/RCPsiSquared.Core/Symmetry/TrulyYParityZeroPurity.cs`](../../compute/RCPsiSquared.Core/Symmetry/TrulyYParityZeroPurity.cs) (`TrulyCriterionHolds`: per-dephase truly criterion encoded in C#)

## Abstract

F102 surfaced y-parity as a real third axis of the polarity cube. F103, F105, F106 then mapped out empirically how the F87 trichotomy splits along that axis, and one pattern stood out clearly across every anchor: F87-truly classifications always land on y-parity zero, never y-parity one. The question is whether this purity is an accident of the specific (N, k) regimes tested, or a structural truth that survives at any chain length and any body count under any of the three dephase letters.

The answer is the structural one. Under each of the three dephase letters, the F87-truly criterion forces y-parity to vanish: every truly Pauli term has an even number of Y letters, full stop. The proof imports F85's structural Z-dephasing truly criterion (a Pauli term contributes M=0 iff it has #Y and #Z both even) and transports it to the X- and Y-dephasing cases via the Π letter-cycle permutations that F108 establishes. All three criteria include "#Y even" as a sub-condition; the rest is bookkeeping.

The diagnostic upshot is that y-parity zero is universal in F87-truly classification. A measured truly-class Pauli pair that carries y-parity one would be either a hardware bug, a non-standard dephase channel, or evidence of a missing classification axis that the polarity cube does not yet capture. The 4524 empirically observed truly classifications across F103+F105+F106 all sit on y-parity zero, which is what F107 closes by closed-form.

The companion proof F109 closes the other purity statement of the trichotomy: mother soft is y-parity one. F110 explores the harder cells (which are not purity-classified but instead carry an inversion pattern). Together F107 + F109 + F110 + F111 give the full y-parity-axis classification of the F87 trichotomy.

**Statement (Theorem F107):** Under any single-letter dephase channel (Z, X, or Y), if a Pauli term σ_α is classified as truly by the F87 trichotomy, then y_par(σ_α) = (#Y in α) mod 2 = 0.

By extension, for any y_par-homogeneous Pauli pair classified as truly, the pair's shared y_par value is 0. Empirically confirmed across F103/F105/F106 anchors (N=4 k=3, N=5 k=3, N=4 k=4): across the three F103/F105/F106 anchor regimes, zero truly classifications carry y_par=1.

## Proof

### Step 1: per-dephase truly criterion (F85 + Π letter-cycle transport)

[PROOF_F85](PROOF_F85_KBODY_GENERALIZATION.md) Step 3 establishes the **Z-dephasing** truly criterion structurally: a Pauli term σ_α contributes M = 0 under Z-dephasing iff

    #Y(α) even  AND  #Z(α) even.

(F85's derivation combines Π's bond-mirror action with the Z-dephasing dissipator's Z-commutativity; the criterion is body-count-independent and bit-exact-verified at k=2, 3, 4 there.)

The X- and Y-dephasing truly criteria follow by **Π letter-cycle transport**. Each canonical Π for dephase letter D is built so its per-letter action permutes the four Pauli letters in the Klein-Vierergruppe pattern fixed by D, and D commutes with the corresponding dissipator. F85's structural derivation does not single out the letter Z; it relies on the pair (Π's letter permutation, dissipator's commuting letter) being matched, which is true for each D ∈ {Z, X, Y} after substituting D into both:

- **Z-dephasing:** Π_Z permutes (I ↔ X, Y ↔ Z); the Π²-odd letter pair is {Y, Z}; criterion = "#Y even AND #Z even" (F85 Step 3 directly).
- **X-dephasing:** Π_X permutes (I ↔ Z, X ↔ Y); the Π²-odd letter pair is {X, Y}; applying F85's derivation with the substitution Z ↔ X gives criterion = "#X even AND #Y even". The palindrome closure under Π_X is the [F108 Part 2](PROOF_F108_PART2_PI2X_EVEN_ALWAYS_PALINDROMIC.md) result; without that closure the truly criterion would be ill-posed at the operator level.
- **Y-dephasing:** Π_Y permutes (I ↔ X, Y ↔ Z), the same per-letter swap as Π_Z (phases differ; see [F108 Part 3](PROOF_F108_PART3_PI2Y_EVEN_ALWAYS_PALINDROMIC.md)); the Π²-odd letter pair is again {Y, Z}; criterion = "#Y even AND #Z even", identical to Z-dephasing because the Y ↔ Z 2-cycle fixes the pair {#Y, #Z} setwise.

The dependence on dephase letter D therefore gives:

| Dephase letter D | Π²-odd letter pair | Truly criterion |
|------------------|---------------------|------------------|
| Z                | {Y, Z}              | #Y even AND #Z even |
| X                | {X, Y}              | #X even AND #Y even |
| Y                | {Y, Z}              | #Y even AND #Z even |

These match `TrulyYParityZeroPurity.TrulyCriterionHolds` branch-for-branch in the typed claim, and are bit-exact-verified at N=3–5, k=3 and k=4 by the F103/F105/F106 anchors.

### Step 2: every truly criterion includes "#Y even"

Reading the criterion column:

- Z-dephasing truly: #Y even AND #Z even, so #Y even.
- X-dephasing truly: #X even AND #Y even, so #Y even.
- Y-dephasing truly: #Y even AND #Z even, so #Y even.

All three dephase letters force #Y even as a sub-condition.

### Step 3: #Y even ⟹ y_par = 0

By definition, y_par(σ) = (#Y in σ) mod 2. #Y even ⟹ y_par = 0.

### Step 4: y_par-homogeneous pair corollary

A Klein-homogeneous + y_par-homogeneous pair (term1, term2) has shared y_par value y_par(term1) = y_par(term2). The pair is truly iff both terms individually satisfy the truly criterion (per F85, M decomposes term-by-term in the linear regime). Each truly term has y_par = 0; hence pair y_par = 0.

∎

## Empirical confirmation

| Anchor | Cells with non-zero truly counts | Total truly | y_par=1 truly |
|--------|----------------------------------|-------------|----------------|
| F103 (N=4 k=3) | 6 of 12 (Klein × dephase) | 300 | 0 |
| F105 (N=5 k=3) | 6 of 12 (Klein × dephase) | 300 | 0 |
| F106 (N=4 k=4) | 9 of 12 (Klein × dephase) | 3924 | 0 |

Total: 4524 truly classifications observed across the three F103/F105/F106 anchor regimes (each regime is a specific (N, k) point); zero have y_par=1. F107 explains this bit-exactly as a closed-form corollary.

## Cross-letter empirical spot-check at Klein (1,0) Y-dephase, N=4 k=3

Klein (1,0) requires bit_a = #X+#Y odd, bit_b = #Y+#Z even. Y-dephase truly criterion adds #Y even AND #Z even.

From #Y even + #X+#Y odd: #X odd.
From #Y even + #Y+#Z even: #Z even.

So Y-dephase Klein (1,0) truly terms have #X odd, #Y even, #Z even. At k=3 letter sequence: enumerate constrained tuples:

- #X=1, #Y=0, #Z=0 (k_body=1): XII, IXI, IIX → 3 sequences
- #X=1, #Y=2, #Z=0 (k_body=3): XYY, YXY, YYX → 3 sequences
- #X=1, #Y=0, #Z=2 (k_body=3): XZZ, ZXZ, ZZX → 3 sequences
- #X=3, #Y=0, #Z=0 (k_body=3): XXX → 1 sequence

Total: 10 terms. Klein-homogeneous + y_par-homogeneous (all y_par=0) unordered pairs with self-pairs: 10·11/2 = 55. **Matches F103 empirical: 55 y_par=0 truly pairs at Klein (1,0) Y-dephase.** ✓

## Sibling y_par-axis claims

Closed 2026-05-25: F108 Part 1+2+3 (Π²-even palindrome family, Tier1Derived); F109 (MotherSoftYParityOnePurity, Tier1Derived unconditional); F110 (HardCellYInversionPattern, Tier1Derived since 2026-06-10). Together F107+F109+F110 close the y_par-axis F87 trichotomy classification.

## Open

- Closed-form derivation of F110 Aspect C exact ratios. The **42:8 at k=3 is answered** (2026-05-29) by the diagonal-cell hardness rule in [PROOF_F103](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) §6; the **228:0 at k=4** is the sibling case via F111's pure-D template rule. Atomic sub-rules verified-not-yet-palindrome-proven.
- k ≥ 5 empirical confirmation of F103/F106 pattern stability beyond N=4.
- Hardware QPU confirmation at k ≥ 3 (no F87 QPU confirmations exist beyond Marrakesh k=2).
