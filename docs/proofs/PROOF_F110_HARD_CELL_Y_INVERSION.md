# PROOF F110: F87-Hard Cells Exhibit Y-Inversion Pattern

**Status:** Tier 1 Candidate (Aspect A closed-form via F108 Part 1+2+3 + F107 + F109 + F87 dissipator-resonance; Aspect B Y-inversion closed-form at k = N = 4 via sibling F111 Pure-D Template Rule, and at k = 3 via the F103 §6 diagonal-cell rule (atomics verified, palindrome-proof pending); Aspect C k-purity 42:8 ratio at k = 3 likewise derived by the F103 §6 rule)
**Date:** 2026-05-25
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md](PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md)
- [PROOF_F108_PART2_PI2X_EVEN_ALWAYS_PALINDROMIC.md](PROOF_F108_PART2_PI2X_EVEN_ALWAYS_PALINDROMIC.md)
- [PROOF_F108_PART3_PI2Y_EVEN_ALWAYS_PALINDROMIC.md](PROOF_F108_PART3_PI2Y_EVEN_ALWAYS_PALINDROMIC.md)
- [PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md](PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md)
- [PROOF_F109_MOTHER_SOFT_Y_PARITY_ONE_PURITY.md](PROOF_F109_MOTHER_SOFT_Y_PARITY_ONE_PURITY.md)
- F87 dissipator-resonance law (Tier1Derived, `compute/RCPsiSquared.Diagnostics/F87/DissipatorResonanceLaw.cs`, anchored at N=4 k=3 over 294 pairs)

## Abstract

F107 and F109 closed the two clean purity statements of the F87 trichotomy on the y-parity axis: truly is always y-parity zero, mother-soft is always y-parity one. The third class, F87-hard, is where things get more interesting. F110 maps out three aspects of its structure.

The first aspect is the cleanest. F87-hard Pauli pairs appear only in one specific Klein cell per dephase letter, the cell whose Klein index matches the dephase letter itself: Z-hardness lives in Klein (0,1), X-hardness in Klein (1,0), Y-hardness in Klein (1,1). This is closed-form, a direct corollary of F108 (Π²-even bilinears never produce hardness) plus F107 (truly is purity-zero) plus F109 (mother soft is purity-one) plus the F87 dissipator-resonance law that selects the diagonal cell from among the remaining candidates.

The second aspect is the Y-inversion observation. Within each diagonal hard cell, the dominant y-parity equals the y-parity of the dephase letter. For Z- and X-dephasing the diagonal is dominantly y-parity zero (matching Z and X both being y-parity-zero letters); for Y-dephasing the diagonal flips to dominantly y-parity one (matching Y being a y-parity-one letter). At k = N = 4 this dominance is bit-exactly pure (228:0 split per cell), closed-form via the sibling Pure-D Template Rule (F111). At k = 3 the dominance is empirical and not yet closed-form (the 42:8 split).

The third aspect is the k-dependent sharpening. At k = 3 the hard cells split 42:8 with the dominant y-parity carrying 84% of the weight; at k = 4 the same cells go fully pure (100% on the dominant side). The pattern is a sharpening, not a re-shaping. The exact 42:8 ratio at k = 3 is derived (2026-05-29) by the diagonal-cell hardness rule, [PROOF_F103](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) §6; the atomic sub-rules remain verified-not-yet-palindrome-proven.

The diagnostic upshot is that y-parity completes the F87 trichotomy classification: truly = y-parity-zero, mother-soft = y-parity-one, hard-on-diagonal = y-parity-of-the-dephase-letter (dominantly at k=3, purely at k=N=4). Outside the diagonal cell, hardness simply does not occur. The cube has full structure.

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

(iii) **The remaining Π²-D-odd non-mother cells.** After (i) + (ii), the remaining-as-possibly-hard cells per dephase D are TWO: the Π²-D-odd non-mother cells (the two Klein cells whose Π²-D parity is odd; the mother (0,0) is even and already removed by (i)/(ii)).
- D = Z: {(0, 1), (1, 1)}  (Π²-Z-odd means bit_b = 1)
- D = X: {(1, 0), (1, 1)}  (Π²-X-odd means bit_a = 1)
- D = Y: {(0, 1), (1, 1)}  (Π²-Y-odd means bit_b = 1, same axis as Z)

(iv) **F87 dissipator-resonance selects one cell.** The F87 dissipator-resonance law (separately Tier1Derived, encoded in `compute/RCPsiSquared.Diagnostics/F87/DissipatorResonanceLaw.cs`, anchored at N = 4 k = 3 over 294 pairs) selects ONE of the two remaining cells per dephase: Z → (0, 1), X → (1, 0), Y → (1, 1). Combining (i) + (ii) + (iv) gives Aspect A: hard appears only in the diagonal Klein cell.

Aspect A is derived as a corollary of F108 Part 1+2+3 + F107 + F109 + F87 dissipator-resonance law. ∎

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

**Aspect B at k = N = 4 (closed-form, Tier1Candidate):** The sibling Claim F111 (HardCellPureDTemplate, 2026-05-25) sharpens Aspect B at k = N = 4: a pair (P, Q) in the diagonal cell is F87-hard iff at least one of P, Q is a "pure-D template" (length-4 string with only D and I letters). Pure-D templates have y_par = y_par(D) by construction, so the F106 N = 4 k = 4 228:0 split follows immediately. See [PROOF_F111_HARD_CELL_PURE_D_TEMPLATE.md](PROOF_F111_HARD_CELL_PURE_D_TEMPLATE.md). F111 is Tier1Candidate (not Tier1Derived) because subclaim (d) Mixed+Mixed = soft at k=N=4 lacks an operator-level closed-form. At k = 3 the 42:8 dominance remains empirical: F111's Pure-D Template Rule is anchored at k = N = 4 and does not transport down to k = 3 as a 1:1 structural correspondence (the F103 enumeration at k_body=3 admits pure-D letter-sequences only as the single all-D string per diagonal cell, far short of the 8 pure-D templates the k = 4 rule relies on, so the 36 + 192 + 0 decomposition does not reproduce the F103 50-pair hard count).

**Aspect C:** the asymmetry sharpens with k_body. At k = 3 the split is biased (84% : 16%); at k = 4 the split is fully pure (100% : 0%, closed-form via F111 at the k=N=4 anchor). The exact 42:8 ratio at k = 3 is derived by the F103 §6 diagonal-cell rule (atomics verified, palindrome-proof pending).

## 4. Empirical verification

Bit-exact verification via the `HardCellYInversionPatternEnumerationTests` SLOW_F110 trait at all three anchors (k=3 N=4, k=3 N=5, k=4 N=4). The C# test class uses `Z2HomogeneousKBodyEnumeration.Enumerate(k)` + `PauliPairTrichotomy.Classify(chain, terms, dephase)` to re-classify every pair and assert per-cell hard counts match F110's expected split with Y-inversion. 7/7 SLOW_F110 tests pass (k=3 N=4 counts + diagonal-only, k=3 N=5 N-stability, k=4 N=4 counts + diagonal-only, dominant-y_par structural reading at k=3 and k=4).

## 5. Significance

F110 completes the y_par-axis classification of the F87 trichotomy:

- **F107 (Tier1Derived):** truly classifications have y_par = 0 across all dephase letters and all Klein cells.
- **F109 (Tier1Derived):** mother sector Klein (0, 0) soft classifications have y_par = 1 across all dephase letters.
- **F110 (Tier1Candidate, THIS PROOF):** F87-hard classifications appear only in the diagonal Klein cell, with dominant y_par equal to the dephase letter's own y_par (Y-inversion).

Together F107 + F109 + F110 characterize the dominant y_par signature of every F87 trichotomy class. The remaining open work (now that the 42:8 k=3 ratio is derived by the F103 §6 rule, atomics pending a palindrome proof) is the closed-form completion of F111 subclaim (d) Mixed+Mixed = soft, which is what gates F111 (and hence the k=4 228:0 split) from Tier1Candidate to Tier1Derived.

## 6. Open

- Closed-form derivation of the 42:8 (k=3) hard split ratio. **ANSWERED 2026-05-29** by the diagonal-cell hardness rule in [PROOF_F103](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) §6 (all-diagonal templates + single-diagonal adjacency; Y-inversion forced by the templates' y_par; verified N=4,5). The atomic sub-rules remain verified-not-yet-palindrome-proven. (The k = 4 228:0 ratio is closed-form Tier1Candidate via F111, with subclaim (d) Mixed+Mixed = soft as the remaining mechanism gap.)
- k ≥ 5 empirical confirmation: F106 anchors k=4 only at N=4. Predictions for k=5 are unverified.
- Hardware QPU confirmation at k ≥ 3: no F87 QPU confirmations exist beyond Marrakesh k=2.

∎
