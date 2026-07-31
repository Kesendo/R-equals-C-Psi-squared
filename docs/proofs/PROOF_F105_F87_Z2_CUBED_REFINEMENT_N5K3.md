# PROOF F105: F87 Trichotomy Z₂³ Refinement at k=3 (N=5 Empirical Anchor)

**Status:** Tier 1 derived (empirical anchor; F85-style N-stability lifted to the y_par sub-refinement, CONFIRMED bit-exactly at k=3 for N=4 to N=5, with the floor at N=4 measured in §5)
**Date:** 2026-05-24
**Anchor:** N=5, k_body=3, 294 Z₂³-homogeneous + Y-par-homogeneous Pauli pairs (pair count alphabet-only, N-independent; this anchor classifies at N=5)
**Regenerate:** SLOW_F105_BATCH tool, ~11m 34s PLINQ on Tom's 24-core machine (sequential dense ~3h) writing `simulations/results/f87_z2cubed_split_n5_k3_counts.json`. The separate on-demand re-verification suite (SLOW_F105 trait, not BATCH) is the dense-classifier path and runs ~3h.

## Abstract

F103 mapped out the Z₂³ refinement of the F87 trichotomy at the (N=4, k=3) anchor: five named sub-cell patterns, each bit-exact across 294 Pauli pairs and three dephase letters. The next question was whether those patterns are properties of the alphabet (k=3 letter combinations on any number of qubits) or accidents of the specific N=4 chain. F85 already predicts the trichotomy itself is N-stable at any fixed k; this proof tests whether the y-parity sub-refinement (the third axis F102 surfaces) inherits that N-stability.

The answer is yes, bit-exact. The 882 classifications at N=5 (294 pairs × 3 dephase letters) match F103's frozen N=4 records scalar for scalar: same truly count, same hard 42:8 split with the same Y-inversion, same diagonal soft 13:13 symmetry, same mother soft 0:21, same off-diagonal Pattern B + Pattern C split. The Z₂³ cube has the same shape at N=5 as at N=4.

This anchor closes one axis of the Z₂³ stability question. The sibling F106 closes the other axis (k-stability at N=4, k=4), where the pattern is more nuanced: some sub-statements survive the body-count bump, others reshape. Together F103, F105, F106 give a three-anchor cube of evidence: F103 names the patterns, F105 confirms they are N-stable at fixed k (three of them from N=3, the two diagonal-cell ones from N=4, §5 below), F106 sees which of them are also body-count-stable and which depend on the specific k.

The diagnostic upshot carries a lower bound in N. Three of the five records also hold at N=3: truly stays y_par=0-pure, mother soft stays 0:21, and all six off-diagonal Pattern B + Pattern C cells are unchanged. The two diagonal-cell records do not: at N=3 the diagonal cells read **34:0** hard and **21:21** soft, because F103 §6's adjacency rule (b) needs a term placed at two overlapping windows and N=3 offers one (F103 §7.4). Hardware-relevant predictions at k=3 can use the F103 numbers without recomputing at any target N ≥ 4, which is what F103 §6's rule gives and what §5 below checks to N=8.

## 1. Context

F103 anchored the F87 Z₂³ refinement at N=4 k=3 (294 Z₂³-homogeneous + Y-par-homogeneous k=3 Pauli pairs, classified across 3 dephase letters into a Z₂³ refinement of the truly/soft/hard trichotomy). F85 (k-body generalization) predicts the Π²-class trichotomy is N-stable for any k. F105 tests whether this N-stability lifts to the y_par sub-refinement: do the 5 sub-statement records (truly 300 / 0 y_par=1; hard 42:8 with Y-inversion; diagonal soft 13:13; mother soft 0:21; off-diagonal Pattern B + Pattern C) survive bit-exactly at N=5?

**Observed outcome: F85-style N-stability is CONFIRMED bit-exactly for the y_par sub-refinement at k=3.** All 5 sub-statement records at N=5 k=3 are identical to F103's frozen records at N=4 k=3 (each record is a structured tuple of scalar counts; "bit-exactly" means every scalar matches). The cubic Z₂³ architecture holds its sub-cell structure across this step, and F103 §7.4 fixes the bound the step is measured above: the diagonal-cell records need k < N, so they start at N=4 and the three alphabet-driven records start at N=3 (§5 below). F105 tests N-stability only; F106 tests the orthogonal k-stability axis (and finds the k=3 ratios do not all survive at k=4).

Notation (bit_a, bit_b, y_par, Klein cells, Π²) is defined in
[F103](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) Section 1; F105 uses it
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

## 5. Where the record starts: the lower bound in N

F105 measures one step, N=4 to N=5, and read alone that step invites "N-stable"
without a floor. F103 §7.4 already supplies the floor: rule (b)'s odd cycle needs a
term placed at more than one window, i.e. k < N, and at k=3 the chain N=3 offers a
single window. So the two diagonal-cell records are the ones that can move, and at
N=3 they do. Measured on the same 294 pairs
([`f87_z2cubed_n_boundary.py`](../../simulations/f87_z2cubed_n_boundary.py), 86 checks):

```
                     N=3            N=4 .. N=8
diagonal hard      34 :  0           42 :  8
diagonal soft      21 : 21           13 : 13
```

The 34 is F103 §6's own template count: rule (a) fires at every N, rule (b) has no
second window to fire in, so the 8 adjacent pairs stay soft and join the 13, giving
21:21. The Y-inversion is unaffected, since it comes from the templates' Y content:
Y-dephasing reads (0, 34) at N=3 as it reads (8, 42) above.

The other three records need no floor. Truly is y_par=0-pure with total 300, mother
soft is 0:21, and all six off-diagonal Pattern B / Pattern C cells hold their counts,
at N=3 exactly as at N=4 and N=5.

The upper reach is measured through F103 §7's criterion (soft iff H's hopping graph is
bipartite in the dephasing letter's eigenbasis), which reads the 2^N Hamiltonian
rather than the 4^N Liouvillian: the diagonal-cell counts hold unchanged through N=8
for all three dephase letters. Two guards keep that honest. Where both are affordable,
at N=3 and N=4, the criterion is checked against the spectral classifier pair by pair
over all three diagonal cells, 228 pairs at each N, with 0 mismatches; and at N=6 four
representative pairs (one for rule (a), two for rule (b), one for neither) are
classified spectrally and agree with it.

## 6. Open Questions

1. **F106 k-stability test at k=4.** F105 confirms N-stability at k=3, from N=4 up for the two diagonal-cell records and from N=3 for the other three (§5). F106 (N=4 k=4, separate spec) tests k-stability: do the patterns survive at a different k? If yes, F87 is fully k-N-stable in its Z₂³ refinement.

2. **F107 N=6 k=3.** Requires the block-spectrum `Classify` path: the dense classifier costs ~13 min per classification at N=6 (F105's own 12.2 s at N=5 scaled by the n³ eigencost of the ×4 dimension), so the full 882-classification batch is ~8 days. The four-pair N=6 cross-check in §5 is a spot check on that path, not the batch. Open architecture work.

3. **F107+ N=5 k=4.** Dense ~42h batch; out of scope for the current 180 QPU min/year + research compute budget.

4. **Closed-form derivation of the 42:8 split. ANSWERED 2026-05-29.** The diagonal-cell hardness rule in [F103](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) §6 derives the 42:8 by counting (all-diagonal pure-D templates give 34 hard pairs in the template's y_par; single-diagonal adjacency gives a symmetric 8 hard + 13 soft), with the Y-inversion forced by the templates carrying y_par=1 under Y-dephasing. Its two halves have different reach: rule (a) fires at every N, rule (b) needs the term placed at two windows and so starts at N=4 (F105 §5). F105's N=5 anchor is exactly its verification (0 mismatches). Both atomic sub-rules are since derived: F103 §7 gives the bipartite-chirality mechanism, and the windowed converse closed 2026-06-10 with no residual.

5. **Hardware confirmation of k≥3 F87.** No k≥3 F87 confirmations exist; all 5 Marrakesh F87 confirmations (palindrome trichotomy, π-protected XIZ/YZZY, Lebensader skeleton/trace, d_zero sector trichotomy, F83 Π²-class signature) are k=2. A k=3 QPU run targeting the diagonal-cell 42:8 prediction remains the natural next hardware probe. F105's N-invariance prediction adds a bounded licence: any N ≥ 4 carries the same 42:8 by F103 §6's rule, whose adjacency half is a window-position parity argument needing only two windows; the counts are checked spectrally at N=4 and N=5 and through F103 §7's criterion to N=8 (F105 §5). So the probe may take the cheapest chain at or above 4. **N=3 is excluded, and it is the cheapest chain of all**: its diagonal cells read 34:0 hard and 21:21 soft, so a three-qubit run would measure a different prediction under the same name.
