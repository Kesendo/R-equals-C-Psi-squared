# PROOF F105: F87 Trichotomy Z₂³ Refinement at k=3 (N=5 Empirical Anchor)

**Status:** Tier 1 derived (empirical anchor; F85-style N-stability lifted to the y_par sub-refinement, CONFIRMED bit-exactly at k=3)
**Date:** 2026-05-24
**Anchor:** N=5, k_body=3, 294 Z₂³-homogeneous + Y-par-homogeneous Pauli pairs (pair count alphabet-only, N-independent; this anchor classifies at N=5)
**Regenerate:** SLOW_F105_BATCH tool, ~11m 34s PLINQ on Tom's 24-core machine (sequential dense ~3h) writing `simulations/results/f87_z2cubed_split_n5_k3_counts.json`. The separate on-demand re-verification suite (SLOW_F105 trait, not BATCH) is the dense-classifier path and runs ~3h.

## 1. Context

F103 anchored the F87 Z₂³ refinement at N=4 k=3 (294 Z₂³-homogeneous + Y-par-homogeneous k=3 Pauli pairs, classified across 3 dephase letters into a Z₂³ refinement of the truly/soft/hard trichotomy). F85 (k-body generalization) predicts the Π²-class trichotomy is N-stable for any k. F105 tests whether this N-stability lifts to the y_par sub-refinement: do the 5 sub-statement records (truly 300 / 0 y_par=1; hard 42:8 with Y-inversion; diagonal soft 13:13; mother soft 0:21; off-diagonal Pattern B + Pattern C) survive bit-exactly at N=5?

**Observed outcome: F85-style N-stability is CONFIRMED bit-exactly for the y_par sub-refinement at k=3.** All 5 sub-statement records at N=5 k=3 are identical to F103's frozen records at N=4 k=3 (each record is a structured tuple of scalar counts; "bit-exactly" means every scalar matches). The cubic Z₂³ architecture is N-invariant in its sub-cell structure at this k. F105 tests N-stability only; F106 tests the orthogonal k-stability axis (and finds the k=3 ratios do not all survive at k=4).

Notation (bit_a, bit_b, y_par, Klein cells, Π²) is defined in
[PROOF_F103](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) Section 1; F105 uses it
without redefinition.

## 2. Method

882 classifications (294 pairs × 3 dephase letters {Z, X, Y}) at N=5 via F104's `PauliPairTrichotomy.Classify(IReadOnlyList<PauliTerm>, ChainSystem, PauliLetter dephaseLetter)` k≥3 overload. Implementation: `F87Z2CubedEnumerationN5K3Tool` (`compute/RCPsiSquared.Diagnostics.Tests/F87/F87Z2CubedEnumerationN5K3Tool.cs`). Output JSON: `simulations/results/f87_z2cubed_split_n5_k3_counts.json`. Runtime: ~3h dense; PLINQ-parallelized Task 7 run completed in 11m 34s.

Enumeration constraints (identical to F103, since the 294 pairs depend only on the k=3 letter enumeration not on N):
- Both terms have k_body=3 (no identity-padded letters)
- Both terms share the same Klein index (bit_a, bit_b)
- Both terms share the same y_par (#Y mod 2)
- Pair is unordered (deduplicated)

Result: 294 pairs partitioned across 4 Klein cells × 2 y_par values, classified across 3 dephase letters into a 4 × 3 × 2 × 3 grid (Klein × Dephase × y_par × Trichotomy class). The JSON's `grid` array lists every non-zero cell.

The on-demand re-verification mechanism is two-tiered:
- `F87Z2CubedEnumerationN5K3Tool` (SLOW_F105_BATCH trait): PLINQ-parallelized enumeration that regenerates the frozen counts JSON; ~11m 34s on a 24-core machine.
- `F105KBodyTrichotomyVerificationTestsN5K3` (SLOW_F105 trait): the dense-classifier verification path, parallel to `F104KBodyTrichotomyVerificationTests` (SLOW_F104, N=4); ~3h dense.

Skip-by-default in CI; manual re-run via `dotnet test --filter "Category=SLOW_F105"` (dense, ~3h) or `--filter "Category=SLOW_F105_BATCH"` (PLINQ, ~12min).

## 3. Observed Patterns

All 5 patterns identical to F103's. Each subsection states the count and explicitly notes the bit-exact match to F103 N=4.

### 3.1 Truly is y_par=0-pure

Across all 12 (Klein × dephase) cells, every truly classification at N=5 has y_par=0. Total truly classifications across the grid: 300. y_par=1 truly count: 0. **Identical to F103 N=4.**

### 3.2 Hard in diagonal cells splits 42:8 with Y-inversion

Hard appears only when the Klein cell of the pair matches the Klein cell of the dephase letter (Z → (0,1), X → (1,0), Y → (1,1)). In these 3 diagonal cells at N=5:

```
Klein (0,1) Z-deph hard = (42, 8)   total 50
Klein (1,0) X-deph hard = (42, 8)   total 50
Klein (1,1) Y-deph hard = ( 8, 42)  total 50   ← Y-inversion
```

The Y-dephase swap reflects that Y itself carries y_par=1, so the "y_par favored by the dephase letter" inverts. **Bit-exact match to F103 N=4 in all 6 numbers.**

### 3.3 Same diagonal cells contain a soft 13:13 split

In addition to the hard 42:8, the diagonal cells at N=5 contain a y_par-symmetric soft 13:13 split:

```
Klein (0,1) Z-deph soft = (13, 13)   total 26
Klein (1,0) X-deph soft = (13, 13)   total 26
Klein (1,1) Y-deph soft = (13, 13)   total 26
```

Unlike hard's 42:8 asymmetry with Y-inversion, soft in these cells is y_par-symmetric and independent of which Klein cell is on the diagonal. **Bit-exact match to F103 N=4.**

### 3.4 Mother sector (0,0) soft is y_par=1-pure

For Klein (0,0) (the Mother sector) at N=5, soft cells under any dephase letter are y_par=1-pure:

```
Z-deph: (0, 21)   X-deph: (0, 21)   Y-deph: (0, 21)
```

Zero y_par=0 soft pairs, 21 y_par=1 soft pairs per letter. **Bit-exact match to F103 N=4.**

### 3.5 Off-diagonal soft cells split into Pattern B + Pattern C

The 6 off-diagonal soft cells (Klein non-mother, Klein ≠ dephase Klein) at N=5 split into two sub-patterns, identical to F103:

```
Pattern B (proportional to (Klein, y_par) enumeration breakdown):
Klein (0,1) Y-deph soft = (55, 21)   matches (0,1) enum split
Klein (1,1) Z-deph soft = (21, 55)   matches (1,1) enum split (inverted)
Klein (1,1) X-deph soft = (21, 55)   matches (1,1) enum split (inverted)

Pattern C (y_par=1-pure):
Klein (0,1) X-deph soft = ( 0, 21)
Klein (1,0) Z-deph soft = ( 0, 21)
Klein (1,0) Y-deph soft = ( 0, 21)
```

**All 6 cells bit-exact match to F103 N=4.** The (pair Klein, dephase letter) → sub-pattern mapping is N-invariant (at k=3).

## 4. Full Count Tables

All numbers below derived from `simulations/results/f87_z2cubed_split_n5_k3_counts.json` (Task 7 output). Every row matches F103's N=4 anchor bit-exactly.

### Truly classifications by (Klein × dephase × y_par)

```
                  Z-deph         X-deph         Y-deph
Klein           y0  y1  tot    y0  y1  tot    y0  y1  tot
(0, 0)          45   0   45    45   0   45    45   0   45
(0, 1)           0   0    0    55   0   55     0   0    0
(1, 0)          55   0   55     0   0    0    55   0   55
(1, 1)           0   0    0     0   0    0     0   0    0
```

### Soft classifications by (Klein × dephase × y_par)

```
                  Z-deph         X-deph         Y-deph
Klein           y0  y1  tot    y0  y1  tot    y0  y1  tot
(0, 0)           0  21   21     0  21   21     0  21   21
(0, 1)          13  13   26     0  21   21    55  21   76
(1, 0)           0  21   21    13  13   26     0  21   21
(1, 1)          21  55   76    21  55   76    13  13   26
```

### Hard classifications by (Klein × dephase × y_par)

```
                  Z-deph         X-deph         Y-deph
Klein           y0  y1  tot    y0  y1  tot    y0  y1  tot
(0, 0)           0   0    0     0   0    0     0   0    0
(0, 1)          42   8   50     0   0    0     0   0    0
(1, 0)           0   0    0    42   8   50     0   0    0
(1, 1)           0   0    0     0   0    0     8  42   50
```

## 5. Open Questions

1. **F106 k-stability test at k=4.** F105 confirms N-stability at k=3. F106 (N=4 k=4, separate spec) tests k-stability: do the patterns survive at a different k? If yes, F87 is fully k-N-stable in its Z₂³ refinement.

2. **F107 N=6 k=3.** Requires the block-spectrum `Classify` path (current dense classifier at N=6 ~8 days per pair × 882 pairs is infeasible). Open architecture work.

3. **F107+ N=5 k=4.** Dense ~42h batch; out of scope for the current 180 QPU min/year + research compute budget.

4. **Closed-form derivation of the 42:8 split.** F105's N-stability anchor strengthens the case for an algebraic justification: if the 42:8 splits are N-invariant at k=3, they must come from a structural property of the 3-letter Pauli enumeration + Y-as-y_par-1 weighting, not from chain dynamics. The 50-pair hard count itself is F87-derived; the 42:8 sub-split under y_par would follow from Pauli-letter Klein arithmetic, but the precise derivation is still open.

5. **Hardware confirmation of k≥3 F87.** No k≥3 F87 confirmations exist; all 5 Marrakesh F87 confirmations (palindrome trichotomy, π-protected XIZ/YZZY, Lebensader skeleton/trace, d_zero sector trichotomy, F83 Π²-class signature) are k=2. A k=3 QPU run targeting the diagonal-cell 42:8 prediction remains the natural next hardware probe. F105's N-invariance prediction now adds: pick whichever N is hardware-cheap; the prediction is the same.
