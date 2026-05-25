# PROOF F110: F87-Hard Cells Exhibit Y-Inversion Pattern

**Status:** Tier 1 Candidate (Aspect A closed-form via F108 Part 1+2+3 + F87 dissipator-resonance; Aspect B Y-inversion and Aspect C k-purity sharpening empirically anchored, closed-form derivation of 42:8 / 228:0 ratios open)
**Date:** 2026-05-25
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md](PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md)
- [PROOF_F108_PART2_PI2X_EVEN_ALWAYS_PALINDROMIC.md](PROOF_F108_PART2_PI2X_EVEN_ALWAYS_PALINDROMIC.md)
- [PROOF_F108_PART3_PI2Y_EVEN_ALWAYS_PALINDROMIC.md](PROOF_F108_PART3_PI2Y_EVEN_ALWAYS_PALINDROMIC.md)
- [PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md](PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md)
- [PROOF_F109_MOTHER_SOFT_Y_PARITY_ONE_PURITY.md](PROOF_F109_MOTHER_SOFT_Y_PARITY_ONE_PURITY.md)
- F87 dissipator-resonance law (anchored at N=4 k=3 over 294 pairs; documented in `compute/RCPsiSquared.Core/Lindblad/PauliDephasingDissipator.cs` docstring)

## 1. Statement

**Aspect A (closed-form):** For any dephase letter D ∈ {Z, X, Y}, F87-hard Pauli pairs appear only in the diagonal Klein cell, the cell whose Klein index matches the dephase letter's own Klein index: Z → (0, 1), X → (1, 0), Y → (1, 1).

**Aspect B (empirical Y-inversion):** Within each diagonal hard cell, the dominant y_par equals y_par(dephase letter). Concretely:
- Z-deph + Klein (0, 1) hard: dominantly y_par = 0
- X-deph + Klein (1, 0) hard: dominantly y_par = 0
- Y-deph + Klein (1, 1) hard: dominantly y_par = 1 (Y-INVERSION)

**Aspect C (empirical k-purity sharpening):**
- k = 3, N = 4 (F103 anchor): 42:8 biased split per diagonal cell
- k = 3, N = 5 (F105 anchor): identical 42:8 (N-stable per F85)
- k = 4, N = 4 (F106 anchor): 228:0 fully pure with Y-inversion preserved

## 2. Proof of Aspect A (closed-form)

Partition all 4 Klein cells into three classes per dephase letter D:

(i) **Π²-D-even cells.** Per F108 Part 1 (D = Z), Part 2 (D = X), Part 3 (D = Y), every Π²-D-even Hamiltonian admits an EXACT operator-level palindrome via the matching Π_5bilinear variant. Hence spec(L) is palindromic in these cells, hence no pair in a Π²-D-even cell can be F87-hard.

For each D, the Π²-D-even cells are:
- D = Z: Klein (0, 0) and (1, 0) (both have bit_b = 0)
- D = X: Klein (0, 0) and (0, 1) (both have bit_a = 0)
- D = Y: Klein (0, 0) and (1, 0) (both have bit_b = 0; same as Z since Y and Z share bit_b parity per `PiOperator.SquaredEigenvalue`)

(ii) **Mother sector Klein (0, 0).** Per F107 (truly y_par = 0) + F109 (mother soft y_par = 1), Klein (0, 0) under each dephase letter contains only truly + soft classifications, never hard. (Klein (0, 0) is also Π²-D-even for all three D, so this is subsumed by (i) but documented separately because F109's structural derivation is independent of F108.)

(iii) **The remaining diagonal cell.** By exclusion, hard can appear only in the single Klein cell that is NOT Π²-D-even and NOT the Mother sector. For each D this is exactly the cell whose Klein index matches D's own Klein index:
- D = Z (Klein index (0, 1)): the remaining cell is (0, 1)
- D = X (Klein index (1, 0)): the remaining cell is (1, 0)
- D = Y (Klein index (1, 1)): the remaining cell is (1, 1)

This is the F87 dissipator-resonance law itself (separately anchored at N = 4 k = 3 over 294 pairs in `PauliDephasingDissipator.cs`). The closed-form derivation here shows it as a corollary of F108 Part 1+2+3 + F107 + F109. ∎

## 3. Aspect B + C (empirical)

From the F103/F105/F106 frozen count tables:

| Anchor | Klein cell | Dephase | Hard split (y_par=0, y_par=1) |
|---|---|---|---|
| F103 N=4 k=3 | (0, 1) | Z | (42, 8), dominantly y_par=0 |
| F103 N=4 k=3 | (1, 0) | X | (42, 8), dominantly y_par=0 |
| F103 N=4 k=3 | (1, 1) | Y | (8, 42), DOMINANTLY y_par=1 (Y-INVERSION) |
| F105 N=5 k=3 | (0, 1) | Z | (42, 8) (N-stable) |
| F105 N=5 k=3 | (1, 0) | X | (42, 8) (N-stable) |
| F105 N=5 k=3 | (1, 1) | Y | (8, 42) (Y-inversion N-stable) |
| F106 N=4 k=4 | (0, 1) | Z | (228, 0), fully pure y_par=0 |
| F106 N=4 k=4 | (1, 0) | X | (228, 0), fully pure y_par=0 |
| F106 N=4 k=4 | (1, 1) | Y | (0, 228), fully pure y_par=1 (Y-INVERSION) |

**Structural reading of Aspect B:** the dephase letter enters the dissipator as a single-letter "preferred" content; in the diagonal hard cell, the y_par favored by the dephase letter's own Y-content dominates. The Y-letter carries y_par = 1, which inverts the otherwise-y_par = 0-preferred pattern.

**Aspect C:** the asymmetry sharpens with k_body. At k = 3 the split is biased (84% : 16%); at k = 4 the split is fully pure (100% : 0%). Closed-form derivation of the exact ratios per Pauli-letter combinatorics is open (F103 Section 5).

## 4. Empirical verification

Bit-exact verification via the `HardCellYInversionPatternEnumerationTests` SLOW_F110 trait at all three anchors (k=3 N=4, k=3 N=5, k=4 N=4). The C# test class uses `Z2HomogeneousKBodyEnumeration.Enumerate(k)` + `PauliPairTrichotomy.Classify(chain, terms, dephase)` to re-classify every pair and assert per-cell hard counts match F110's expected split with Y-inversion. 6/6 SLOW_F110 tests pass.

## 5. Significance

F110 completes the y_par-axis classification of the F87 trichotomy:

- **F107 (Tier1Derived):** truly classifications have y_par = 0 across all dephase letters and all Klein cells.
- **F109 (Tier1Derived):** mother sector Klein (0, 0) soft classifications have y_par = 1 across all dephase letters.
- **F110 (Tier1Candidate, THIS PROOF):** F87-hard classifications appear only in the diagonal Klein cell, with dominant y_par equal to the dephase letter's own y_par (Y-inversion).

Together F107 + F109 + F110 characterize the y_par-axis structure of every F87 trichotomy class. The remaining open work is exclusively Aspect C: the closed-form derivation of the exact 42:8 (k=3) and 228:0 (k=4) ratios.

## 6. Open

- Closed-form derivation of the 42:8 (k=3) and 228:0 (k=4) hard split ratios per Pauli-letter combinatorics. F103 Section 5 explicitly listed as open.
- k ≥ 5 empirical confirmation: F106 anchors k=4 only at N=4. Predictions for k=5 are unverified.
- Hardware QPU confirmation at k ≥ 3: no F87 QPU confirmations exist beyond Marrakesh k=2.

∎
