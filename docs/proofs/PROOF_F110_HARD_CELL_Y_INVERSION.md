# PROOF F110: F87-Hard Cells Exhibit Y-Inversion Pattern

**Status:** Tier 1 derived (Aspect A by the colouring of the three non-diagonal Klein cells, with F158; Aspects B and C by the F103 §6 counting rule at k = 3 and the F111 pure-D template rule at k = N = 4, both resting on the windowed all-γ converse, WindowedConverseAllGammaClaim with the Pascal-Gram positivity of F117)
**Date:** 2026-05-25
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [F103 §8](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) (every non-diagonal Klein cell is coloured) and §6-§7 (the counting rule and the bipartite mechanism)
- [PROOF_PALINDROME_TWO_END_COUNT.md](PROOF_PALINDROME_TWO_END_COUNT.md) (F158 §(e): an invertible operator that commutes with H and anticommutes with every jump reflects L to −L† − 2σ, and the spectrum pairs about −σ)
- [The palindrome as a colouring](../../experiments/THE_PALINDROME_AS_A_COLOURING.md) (a lit Pauli string commuting with H is such an operator)
- [F111](PROOF_F111_HARD_CELL_PURE_D_TEMPLATE.md) (Aspect B at k = N = 4, derived there inside the diagonal cell that Aspect A leaves)
- [F107](PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md), [F109](PROOF_F109_MOTHER_SOFT_Y_PARITY_ONE_PURITY.md) (the two purity statements beside this one)

## Abstract

F107 and F109 closed the two clean purity statements of the F87 trichotomy on the y-parity axis: truly is always y-parity zero, and a mother-soft pair of one y-parity is always y-parity one. The third class, F87-hard, is where things get more interesting. F110 maps out three aspects of its structure.

The first aspect is the cleanest. F87-hard Pauli pairs whose two terms share a Klein cell appear only in one specific cell per dephase letter, the cell whose Klein index matches the dephase letter itself: Z-hardness lives in Klein (0,1), X-hardness in Klein (1,0), Y-hardness in Klein (1,1). This holds at every body count, every coupling and every per-site rate, because each of the other three cells is coloured: some letter string that anticommutes with every jump commutes with every string of the cell, and by F158 that makes the spectrum palindromic. Only the diagonal cell's strings anticommute with both such letter strings.

The second aspect is the Y-inversion observation. Within each diagonal hard cell, the dominant y-parity equals the y-parity of the dephase letter. For Z- and X-dephasing the diagonal is dominantly y-parity zero (matching Z and X both being y-parity-zero letters); for Y-dephasing the diagonal flips to dominantly y-parity one (matching Y being a y-parity-one letter). At k = N = 4 this dominance is bit-exactly pure (228:0 split per cell), closed-form via the sibling Pure-D Template Rule (F111). At k = 3 the dominance is derived by the F103 §6 counting rule (the 42:8 split).

The third aspect is the k-dependent sharpening. At k = 3 the hard cells split 42:8 with the dominant y-parity carrying 84% of the weight; at k = 4 the same cells go fully pure (100% on the dominant side). The pattern is a sharpening, not a re-shaping. The exact 42:8 ratio at k = 3 is derived (2026-05-29) by the diagonal-cell hardness rule, [F103](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) §6, whose hard direction is the windowed all-γ converse (§3).

The diagnostic upshot is that y-parity reads the truly pairs, the Mother cell's soft pairs and the hard pairs: truly = y-parity-zero, mother-soft of one y-parity = y-parity-one, hard-on-diagonal = y-parity-of-the-dephase-letter (dominantly at k=3, purely at k=N=4); the soft pairs of the other cells are F103's (§3.3, §8). Outside the diagonal cell, hardness does not occur for pairs whose terms share a Klein cell.

## 1. Statement

**Aspect A (closed-form):** For any dephase letter D ∈ {Z, X, Y}, F87-hard Klein-homogeneous Pauli pairs (both terms in one Klein cell, as in the F103/F105/F106 enumerations) appear only in the diagonal Klein cell, the cell whose Klein index matches the dephase letter's own Klein index: Z → (0, 1), X → (1, 0), Y → (1, 1).

**Aspect B (Y-inversion, derived):** Within each diagonal hard cell, the dominant y_par equals y_par(dephase letter). Concretely:
- Z-deph + Klein (0, 1) hard: dominantly y_par = 0
- X-deph + Klein (1, 0) hard: dominantly y_par = 0
- Y-deph + Klein (1, 1) hard: dominantly y_par = 1 (Y-INVERSION)

**Aspect C (k-purity sharpening, derived via F103 §6/§7):**
- k = 3, N = 4 (F103 anchor): 42:8 biased split per diagonal cell
- k = 3, N = 5 (F105 anchor): identical 42:8 (N-stable from N = 4 per F103 §6; F85's own N-stability is the per-term Π²-class, a different cut of the word, see `experiments/SOFTNESS_IS_N_DEPENDENT.md`)
- k = 4, N = 4 (F106 anchor): 228:0 fully pure with Y-inversion preserved

## 2. Proof of Aspect A

Every string has a Klein letter K, the product of its letters up to phase: I for (0, 0), X for (1, 0), Z for (0, 1), Y for (1, 1). A string commutes with X^⊗N when bit_b = 0, with Z^⊗N when bit_a = 0 and with Y^⊗N when bit_a = bit_b, so a string of Klein letter K ≠ I commutes with K^⊗N and anticommutes with the other two letter strings, and a string of Klein letter I commutes with all three. Under dephasing by the letter D, a letter string a^⊗N anticommutes with every jump D_l exactly when a ∉ {I, D}: the two lit letters.

- **K ∉ {I, D}**, the two off-diagonal cells: K itself is lit, and K^⊗N commutes with every string of the cell.
- **K = I**, the Mother sector: its strings commute with both lit letter strings.

In these three cells a lit letter string F commutes with every Hamiltonian H of the cell and anticommutes with every jump. Right multiplication ρ ↦ ρ·F then keeps −i[H, ·] and reflects the dissipator, R·L·R⁻¹ = −L† − 2σ ([F109](PROOF_F109_MOTHER_SOFT_Y_PARITY_ONE_PURITY.md) Step 5 writes the two lines out), and by the sufficiency step of [F158](PROOF_PALINDROME_TWO_END_COUNT.md) §(e) and hermiticity preservation the spectrum pairs about −σ, for every coupling and every per-site rate. No pair of these cells is F87-hard. This is the colouring of [the palindrome as a colouring](../../experiments/THE_PALINDROME_AS_A_COLOURING.md); in the two off-diagonal cells it completes to the operator identity W·L·W⁻¹ = −L − 2σ with W(ρ) = D^⊗N·ρ·K^⊗N·D^⊗N ([F103 §8](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md)).

- **K = D**, the diagonal cell: its strings anticommute with both lit letter strings, and no letter string colours it. This is where hard pairs occur (§3).

So hard Klein-homogeneous pairs appear only in the diagonal Klein cell: Z → (0, 1), X → (1, 0), Y → (1, 1). Of the three coloured cells, the Mother sector and the cell of the canonical mirror's flip letter are the Π²-D-even ones (Π_D² is the turn by that letter string), and the third is Π²-D-odd (in [F103 §8](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) the flip letter's cell is Pattern C and the third Pattern B). The F87 dissipator-resonance law (`compute/RCPsiSquared.Diagnostics/F87/DissipatorResonanceLaw.cs`) is the census of the same statement at N = 4, k = 3: 50 hard pairs of 76 in the matched cell and none in the other three, under each letter. [`f107_f110_route_gate.py`](../../simulations/f107_f110_route_gate.py) checks the colouring identity exactly at N = 4 on random Hamiltonians of every non-diagonal cell, any body count and random per-site rates, and checks that no lit letter string colours the diagonal cell. ∎

## 3. Aspect B + C (derived via F103 §6/§7, anchored by F103/F105/F106)

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

The three F105 rows read "N-stable" from N = 4 upward, and F103 §7's criterion carries
them through N = 8. At N = 3 the same three cells read (34, 0), (34, 0) and (0, 34): the
split with F103 §6's adjacency rule (b) absent, since a single window closes no odd cycle
([F105 §5](PROOF_F105_F87_Z2_CUBED_REFINEMENT_N5K3.md)). The Y-inversion survives that
step, because it comes from the templates' Y content rather than from the adjacency half.

**Structural reading of Aspect B:** the dephase letter enters the dissipator as a single-letter "preferred" content; in the diagonal hard cell, the y_par favored by the dephase letter's own Y-content dominates. The Y-letter carries y_par = 1, which inverts the otherwise-y_par = 0-preferred pattern.

**Aspect B at k = N = 4 (closed-form, Tier1Derived):** The sibling Claim F111 (HardCellPureDTemplate, 2026-05-25, Tier1Derived since 2026-06-10) sharpens Aspect B at k = N = 4: a pair (P, Q) in the diagonal cell is F87-hard iff at least one of P, Q is a "pure-D template" (length-4 string with only D and I letters). Pure-D templates have y_par = y_par(D) by construction, so the F106 N = 4 k = 4 228:0 split follows immediately. See [F111](PROOF_F111_HARD_CELL_PURE_D_TEMPLATE.md). (F111 was promoted to Tier1Derived once subclaim (d) Mixed+Mixed = soft closed modulo M via PROOF_F103 §7.4 and the hard-direction converse closed via WindowedConverseAllGammaClaim.) At k = 3 the 42:8 dominance is derived instead by the F103 §6 counting rule, not by F111's k = 4 rule: F111's Pure-D Template Rule is anchored at k = N = 4 and does not transport down to k = 3 as a 1:1 structural correspondence (the F103 enumeration at k_body=3 admits pure-D letter-sequences only as the single all-D string per diagonal cell, far short of the 8 pure-D templates the k = 4 rule relies on, so the 36 + 192 + 0 decomposition does not reproduce the F103 50-pair hard count).

**Aspect C:** the asymmetry sharpens with k_body. At k = 3 the split is biased (84% : 16%); at k = 4 the split is fully pure (100% : 0%, closed-form via F111 at the k=N=4 anchor). The exact 42:8 ratio at k = 3 is derived by the F103 §6 diagonal-cell rule; the windowed hard-direction converse it relied on closed 2026-06-10 (WindowedConverseAllGammaClaim, no residual).

## 4. Empirical verification

Bit-exact verification via the `HardCellYInversionPatternEnumerationTests` SLOW_F110 trait at four anchors (k=3 N=3, k=3 N=4, k=3 N=5, k=4 N=4). The C# test class uses `Z2HomogeneousKBodyEnumeration.Enumerate(k)` + `PauliPairTrichotomy.Classify(chain, terms, dephase)` to re-classify every pair and assert per-cell hard counts match F110's expected split with Y-inversion. 9/9 SLOW_F110 tests pass (k=3 N=4 counts + diagonal-only, k=3 N=5 N-stability, k=3 N=3 below the floor + the three records that hold there, k=4 N=4 counts + diagonal-only, dominant-y_par structural reading at k=3 and k=4).

## 5. Significance

F110 is the third of the y_par-axis statements on the F87 trichotomy:

- **F107 (Tier1Derived):** truly classifications have y_par = 0 across all dephase letters and all Klein cells.
- **F109 (Tier1Derived):** mother sector Klein (0, 0) soft pairs of one y-parity have y_par = 1 across all dephase letters.
- **F110 (THIS PROOF; Tier1Derived since 2026-06-10):** F87-hard Klein-homogeneous pairs appear only in the diagonal Klein cell, with dominant y_par equal to the dephase letter's own y_par (Y-inversion).

Together F107 + F109 + F110 give the y_par signature of the truly pairs, the Mother cell's soft pairs and the hard pairs; the soft pairs of the other cells are F103's (§3.3, §8). F111 subclaim (d) Mixed+Mixed = soft is closed modulo M via PROOF_F103 §7.4, and the hard-direction converse is the windowed all-γ theorem (WindowedConverseAllGammaClaim, Pascal-Gram positivity F117, no residual); both F110 and F111 are Tier1Derived.

## 6. Open

- Closed-form derivation of the 42:8 (k=3) hard split ratio. **ANSWERED 2026-05-29** by the diagonal-cell hardness rule in [F103](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) §6 (all-diagonal templates + single-diagonal adjacency; Y-inversion forced by the templates' y_par; verified N=4,5). The windowed hard-direction converse the atomic sub-rules relied on closed 2026-06-10 (WindowedConverseAllGammaClaim, no residual). (The k = 4 228:0 ratio is closed-form via F111, Tier1Derived since 2026-06-10; subclaim (d) Mixed+Mixed = soft closed modulo M via PROOF_F103 §7.4.)
- k ≥ 5 empirical confirmation: F106 anchors k=4 only at N=4. Predictions for k=5 are unverified.
- Hardware QPU confirmation at k ≥ 3: no F87 QPU confirmations exist beyond Marrakesh k=2.

∎
