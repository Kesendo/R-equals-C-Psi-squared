# PROOF F106: F87 Trichotomy Z₂³ Refinement at k=4 (N=4 Empirical Anchor)

**Status:** Tier 1 derived (empirical anchor; mixed k-stability outcome at k=4)
**Date:** 2026-05-24
**Anchor:** N=4, k_body=4, 4248 Z₂³-homogeneous Pauli pairs
**Regenerate:** SLOW_F106_BATCH tool (~2-3min PLINQ on Tom's 24-core machine, ~40min sequential) writing `simulations/results/f87_z2cubed_split_n4_k4_counts.json`

## 1. Context

F103 anchored the F87 Z₂³ refinement at N=4 k=3 (294 Z₂³-homogeneous + Y-par-homogeneous k=3 Pauli pairs). F105 confirmed F85's N-stability lift at the y_par axis: identical counts at N=5 k=3 to F103. F106 tests the orthogonal axis (k-stability) at N=4 k=4. F85 does NOT predict k-stability of the y_par sub-refinement; specifically, the Klein (0,0) enum balance shifts from 45/21 at k=3 to 780/300 at k=4, and the structural worry was that "mother soft is y_par=1-pure" would break.

**Observed outcome (mixed):** Three patterns SURVIVED in some form at k=4; two patterns BROKE structurally. Concretely:

- Truly y_par=0-purity: **HELD bit-structurally** (every truly classification at k=4 has y_par=0; total 3924, y_par=1 count 0).
- Mother soft y_par=1-purity: **HELD bit-structurally** (every mother-soft classification still has y_par=1; counts moved from (0, 21) at k=3 to (0, 300) at k=4 across all 3 dephase letters, but purity is intact).
- Hard diagonal Y-inversion: **HELD qualitatively** (Y-dephase still inverts the y_par assignment relative to Z/X-dephase).
- Hard diagonal 42:8 ratio: **BROKE** (the mixed-y_par 42:8 split sharpened to fully polarized 228:0 / 0:228, i.e. 100% pure per cell rather than 84%:16% mixed).
- Diagonal soft 13:13 symmetry: **BROKE** (the per-diagonal-cell y_par symmetry became asymmetric 300:528 for Z/X and 528:300 for Y).

The k=4 enumeration is genuinely different (4248 vs 294 pairs; 12744 vs 882 classifications); some F103 sub-statements were y_par-axis invariants that genuinely lift across k (the two purity statements survive bit-structurally), while others were k-specific ratios that re-scale with the enum balance shift.

## 2. Method

12744 classifications (4248 pairs × 3 dephase letters {Z, X, Y}) at N=4 via F104's `PauliPairTrichotomy.Classify(IReadOnlyList<PauliTerm>, ChainSystem, PauliLetter dephaseLetter)` k≥3 overload. Implementation: `F87Z2CubedEnumerationN4K4Tool` (`compute/RCPsiSquared.Diagnostics.Tests/F87/F87Z2CubedEnumerationN4K4Tool.cs`). Output JSON: `simulations/results/f87_z2cubed_split_n4_k4_counts.json`. Runtime: 3m 59s PLINQ on 24 cores.

Enumeration constraints (k_body now 4, not 3; pair count is 4248):
- Both terms have k_body=4 (no identity-padded letters; total non-identity 4^4 - 1 = 255 letter strings per term before pair filter)
- Both terms share the same Klein index (bit_a, bit_b)
- Both terms share the same y_par (#Y mod 2)
- Pair is unordered (deduplicated)

Shared enumerator: `Z2HomogeneousKBodyEnumeration.Enumerate(int k)` (replaces 3 previously duplicated k=3-only copies in F104 fixture, F105 fixture, F105 tool). At k=4 the base-4 indexing materializes all 256 k-letter sequences once; filtered Cartesian product produces 4248 unordered pairs.

The on-demand re-verification mechanism lives in `F106KBodyTrichotomyVerificationTestsN4K4` (SLOW_F106 trait). Skip-by-default in CI; manual re-run via `dotnet test --filter "Category=SLOW_F106"` takes ~2-3min PLINQ.

## 3. Observed Patterns

### 3.1 Truly is y_par=0-pure (HELD)

Across all 6 (Klein × dephase) cells where truly appears at k=4, every truly classification has y_par=0. Total truly classifications across the grid: 3924. y_par=1 truly count: 0.

```
Klein (0,0) X-deph truly y_par=0 = 780
Klein (0,0) Y-deph truly y_par=0 = 780
Klein (0,0) Z-deph truly y_par=0 = 780
Klein (0,1) X-deph truly y_par=0 = 528    (off-diagonal: dephase=X, Klein=(0,1))
Klein (1,0) Y-deph truly y_par=0 = 528    (off-diagonal: dephase=Y, Klein=(1,0))
Klein (1,0) Z-deph truly y_par=0 = 528    (off-diagonal: dephase=Z, Klein=(1,0))
```

Total = 3 × 780 + 3 × 528 = 3924. Y_par=0-purity from k=3 is preserved structurally at k=4.

### 3.2 Hard in diagonal cells: 42:8 broke, Y-inversion held (MIXED)

Hard appears only when the Klein cell of the pair matches the Klein cell of the dephase letter (Z → (0,1), X → (1,0), Y → (1,1)). In these 3 diagonal cells at k=4:

```
Klein (0,1) Z-deph hard = (228, 0)    total 228   pure y_par=0
Klein (1,0) X-deph hard = (228, 0)    total 228   pure y_par=0
Klein (1,1) Y-deph hard = (0, 228)    total 228   pure y_par=1   ← Y-inversion
```

The qualitative Y-inversion structure HELD (Y still inverts to the opposite y_par). The 42:8 mixed-y_par split at k=3 BROKE: each cell at k=4 is fully polarized to a single y_par value (228:0 or 0:228). The total per cell (228) replaces k=3's (42 + 8 = 50).

### 3.3 Diagonal cells contain soft: 13:13 broke (BROKE)

In the same 3 diagonal cells at k=4, soft has an asymmetric split (vs k=3's symmetric 13:13):

```
Klein (0,1) Z-deph soft = (300, 528)    total 828
Klein (1,0) X-deph soft = (300, 528)    total 828
Klein (1,1) Y-deph soft = (528, 300)    total 828   ← Y-inverted asymmetric
```

The y_par symmetry from k=3 (each cell was (13, 13)) BROKE: now asymmetric, with the same Y-inversion as the hard split. Total per cell (828) replaces k=3's (13 + 13 = 26).

### 3.4 Mother (Klein (0,0)) soft: y_par=1-purity held (HELD)

Klein (0,0) under each of the 3 dephase letters:

```
Klein (0,0) Z-deph soft = (0, 300)    pure y_par=1
Klein (0,0) X-deph soft = (0, 300)    pure y_par=1
Klein (0,0) Y-deph soft = (0, 300)    pure y_par=1
```

This was the structurally-worried-about case (enum balance shifted from 45/21 to 780/300, so the framework prediction in the F106 plan was that this would break). It DID NOT break: mother soft remains y_par=1-pure at k=4 across all 3 letters. The counts re-scale (k=3 had (0, 21), k=4 has (0, 300)), but the purity is preserved.

### 3.5 Off-diagonal soft: 6 cells, two-tier structure (MIXED)

6 off-diagonal cells (Klein non-(0,0) and Klein != dephase Klein):

```
Klein (0,1) X-deph soft = (0, 528)      Pattern C analog: pure y_par=1
Klein (0,1) Y-deph soft = (528, 528)    NEW: fully y_par-symmetric
Klein (1,0) Y-deph soft = (0, 528)      Pattern C analog: pure y_par=1
Klein (1,0) Z-deph soft = (0, 528)      Pattern C analog: pure y_par=1
Klein (1,1) X-deph soft = (528, 528)    NEW: fully y_par-symmetric
Klein (1,1) Z-deph soft = (528, 528)    NEW: fully y_par-symmetric
```

Three cells (those associated with X-dephase at Klein (0,1) and Y/Z-dephase at Klein (1,0)) match F105's Pattern C (y_par=1-pure (0, x), here x=528). The other three cells (which at k=3 were Pattern B asymmetric (55, 21) or (21, 55)) became FULLY SYMMETRIC (528, 528) at k=4.

The Pattern C / Pattern B split from k=3 is preserved as a 3+3 partition; only the asymmetric Pattern B values changed shape (from (55, 21)/(21, 55) to symmetric (528, 528)).

## 4. Full Count Tables

The 4 Klein × 3 dephase × 2 y_par × 3 trichotomy class = 72-cell grid is in `simulations/results/f87_z2cubed_split_n4_k4_counts.json` and reproduced in the frozen records of `F87Z2CubedRefinementN4K4.cs`. The JSON omits zero cells; only 27 cells have non-zero counts (12744 classifications total).

Grand totals reconciliation:
- Truly total = 3924 (780 × 3 mother + 528 × 3 off-diagonal-y_par=0)
- Hard total = 684 (228 × 3 diagonals)
- Soft total = 8136 (mother 900 + diagonal 2484 + off-diagonal 4752)
- Grand total = 3924 + 684 + 8136 = 12744 ✓

## 5. Open Questions

- Closed-form derivation of the (228, 0) / (0, 228) hard diagonal at k=4, and the (300, 528) / (528, 300) diagonal soft split (analogous to the long-standing F103 42:8 closed-form open question; the new ratios may be more tractable due to their polarization)
- Closed-form for the off-diagonal Pattern B → fully symmetric transition at k=4 (the cells that were asymmetric (55, 21)/(21, 55) at k=3 became (528, 528))
- N>4 at k=4 (~42h dense per batch at N=5; impractical without block-spectrum Classify)
- k>4 (impractical; pair count grows by factor of ~16 per +1 to k; 67968 at k=5 ahead of the 4248 at k=4)
- Hardware k≥3 F87 confirmation (still open from F103/F105)
