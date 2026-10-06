# PROOF F107: F87 Truly Classification Forces y_par = 0 (All Dephase Letters)

**Status:** Tier 1 derived (the Z criterion from Π_Z = R·D at every body count, carried to X and Y by the exact transports between the canonical mirrors)
**Date:** 2026-05-24
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [F85 k-body generalization](PROOF_F85_KBODY_GENERALIZATION.md) (the Z-dephasing truly criterion #Y even AND #Z even, stated and verified at k = 2, 3, 4)
- [the Π factorization](PROOF_PI_FACTORS_AS_R_TIMES_D.md) (§4(d) and (f): the Z criterion from Π_Z = R·D at every body count; §4(a): Π_Y = Π_Z⁻¹)
- [the Klein-V₄ dephase swaps](PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md) (Q_zx·Π_Z·Q_zx⁻¹ = Π_X, Q_zx the conjugation by the Hadamard)
- [the Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md) (the canonical Π_Z and its palindrome)
- [`compute/RCPsiSquared.Core/Symmetry/PiOperator.cs`](../../compute/RCPsiSquared.Core/Symmetry/PiOperator.cs) (the canonical Π_Z, Π_X, Π_Y and their Π² eigenvalue rules)
- [`compute/RCPsiSquared.Core/Symmetry/TrulyYParityZeroPurity.cs`](../../compute/RCPsiSquared.Core/Symmetry/TrulyYParityZeroPurity.cs) (`TrulyCriterionHolds`: per-dephase truly criterion encoded in C#)

## Abstract

F102 surfaced y-parity as a real third axis of the polarity cube. F103, F105, F106 then mapped out empirically how the F87 trichotomy splits along that axis, and one pattern stood out clearly across every anchor: F87-truly classifications always land on y-parity zero, never y-parity one. The question is whether this purity is an accident of the specific (N, k) regimes tested, or a structural truth that survives at any chain length and any body count under any of the three dephase letters.

The answer is the structural one. Under each of the three dephase letters, the F87-truly criterion forces y-parity to vanish: every truly Pauli term has an even number of Y letters. The proof takes the Z-dephasing criterion (a Pauli term contributes M = 0 iff it has #Y and #Z both even) from the factorization Π_Z = R·D, which derives it at every body count, and carries it to X- and Y-dephasing by the exact transports between the canonical mirrors: Π_Y is the inverse of Π_Z and Π_X its Hadamard conjugate. All three criteria include "#Y even" as a sub-condition; the rest is bookkeeping. The criteria are those of the canonical palindromizers, their letters (X or Z) and their phases together. Against other palindromizers of the same dissipator truly strings can carry y-parity one: the quarter-turned mirror of Z- or X-dephasing flips by Y and asks for #X and #Z even, which puts the truly strings of the Y cell at y-parity one, and F108 Part 1's Π_5bilinear keeps the letter X with other phases and puts those of the X cell there ([F103 §8](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md)).

The diagnostic upshot is that y-parity zero is universal in F87-truly classification, which is defined against the canonical palindromizers. A measured truly-class Pauli pair that carries y-parity one would be either a hardware bug, a non-standard dephase channel, a palindromizer other than the canonical one, or evidence of a missing classification axis that the polarity cube does not yet capture. The 4524 empirically observed truly classifications across F103+F105+F106 all sit on y-parity zero, which is what F107 closes by closed-form.

The companion proof F109 closes the other purity statement of the trichotomy: mother-soft pairs of one y-parity are y-parity one. F110 explores the harder cells (which are not purity-classified but instead carry an inversion pattern). Together F107 + F109 + F110 + F111 give the y-parity signature of the truly pairs, the Mother cell's soft pairs and the hard pairs; the soft pairs of the other cells are F103's (§3.3, §8).

**Statement (Theorem F107):** Under any single-letter dephase channel (Z, X, or Y), if a Pauli term σ_α is classified as truly by the F87 trichotomy, then y_par(σ_α) = (#Y in α) mod 2 = 0.

By extension, for any y_par-homogeneous Pauli pair classified as truly, the pair's shared y_par value is 0. Empirically confirmed across F103/F105/F106 anchors (N=4 k=3, N=5 k=3, N=4 k=4): across the three F103/F105/F106 anchor regimes, zero truly classifications carry y_par=1.

## Proof

### Step 1: per-dephase truly criterion

For **Z-dephasing** the criterion is derived at every body count in [the Π factorization](PROOF_PI_FACTORS_AS_R_TIMES_D.md) §4(d) and §4(f). There Π_Z = R·D, where D is the transpose ρ ↦ ρᵀ (that proof's name for it), which multiplies L_σ by (−1)^(n_Y+1) (F114) and R (right multiplication by X^⊗N) keeps L_σ when n_Y + n_Z is even and turns it into an anticommutator otherwise. So Π_Z·L_σ·Π_Z⁻¹ = −L_σ exactly when

    #Y(σ) even  AND  #Z(σ) even,

the criterion [F85](PROOF_F85_KBODY_GENERALIZATION.md) states and verifies at k = 2, 3, 4. The canonical mirrors of the other two letters are exact transports of Π_Z, phases included:

- **Y-dephasing:** Π_Y = Π_Z⁻¹ (§4(a) there). An operator and its inverse flip L_σ for the same strings, so the criterion is the one for Z. That Π_Y palindromizes the Y dissipator, so that M is again the Hamiltonian part alone, comes from the quarter turn R_x(π/2)^⊗N, which carries Π_Z to Π_Y and L_Z to L_Y ([the Klein-V₄ dephase swaps](PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md), implications).
- **X-dephasing:** Π_X = Q_zx·Π_Z·Q_zx⁻¹ ([the Klein-V₄ dephase swaps](PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md)), with Q_zx the conjugation ρ ↦ U·ρ·U† by U = U_H^⊗N, U_H the Hadamard. It carries L_σ to L_{UσU†} and exchanges X and Z letter by letter (Y changes sign), so σ is truly against Π_X exactly when UσU† is truly against Π_Z: #X even AND #Y even.

| Dephase letter | Its canonical mirror flips by | Truly criterion |
|----------------|-------------------------------|------------------|
| Z                | X            | #Y even AND #Z even |
| X                | Z            | #X even AND #Y even |
| Y                | X            | #Y even AND #Z even |

The two letters that must be even are the two other than the letter the mirror flips by. Per site that is the whole mechanism: each canonical mirror carries left multiplication by a letter to plus or minus right multiplication by the same letter and back, with the minus signs on exactly those two letters, one on each side (Π_Z: Z on the left, Y on the right; Π_X: X on the left, Y on the right; Π_Y: Y on the left, Z on the right). For a string both sides are products over sites, and M_σ = 0 exactly when both sign products are +1. The criteria belong to the canonical mirrors, their letters and their phases together: a mirror with the same letter swap and other phases has other truly strings ([F103 §8](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md)).

These match `TrulyYParityZeroPurity.TrulyCriterionHolds` branch for branch in the typed claim. [`f107_f110_route_gate.py`](../../simulations/f107_f110_route_gate.py) checks the per-site signs of the three canonical mirrors exactly and compares the criterion with the directly computed M_σ on every non-identity string at N = 4, body counts 1 to 4.

### Step 2: every truly criterion includes "#Y even"

Reading the criterion column:

- Z-dephasing truly: #Y even AND #Z even, so #Y even.
- X-dephasing truly: #X even AND #Y even, so #Y even.
- Y-dephasing truly: #Y even AND #Z even, so #Y even.

All three dephase letters force #Y even as a sub-condition.

### Step 3: #Y even ⟹ y_par = 0

By definition, y_par(σ) = (#Y in σ) mod 2. #Y even ⟹ y_par = 0.

### Step 4: y_par-homogeneous pair corollary

A Klein-homogeneous + y_par-homogeneous pair (term1, term2) has shared y_par value y_par(term1) = y_par(term2). The pair is truly iff both terms individually satisfy the truly criterion: M is linear in H, and the part a string σ contributes sends each string S to a multiple of the one string σ·S, so different strings never cancel ([F103 §8](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md)). Each truly term has y_par = 0; hence pair y_par = 0.

∎

## Empirical confirmation

| Anchor | Cells with non-zero truly counts | Total truly | y_par=1 truly |
|--------|----------------------------------|-------------|----------------|
| F103 (N=4 k=3) | 6 of 12 (Klein × dephase) | 300 | 0 |
| F105 (N=5 k=3) | 6 of 12 (Klein × dephase) | 300 | 0 |
| F106 (N=4 k=4) | 6 of 12 (Klein × dephase) | 3924 | 0 |

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

F109 (MotherSoftYParityOnePurity, Tier1Derived); F110 (HardCellYInversionPattern, Tier1Derived); F111 (HardCellPureDTemplate, Tier1Derived).

## Open

- k ≥ 5 empirical confirmation of F103/F106 pattern stability beyond N=4, for the soft, truly and off-diagonal cells (the hard diagonal cell is read at k = 4, 5 by `DiagonalCellComplementConnectionTests`).
- Hardware QPU confirmation at k ≥ 3 (no F87 QPU confirmations exist beyond Marrakesh k=2).
