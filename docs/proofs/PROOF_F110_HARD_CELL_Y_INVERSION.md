# PROOF F110: F87-Hard Pairs Only in Diagonal Klein Cells with Y-Inversion

**Status:** Tier1Candidate (Aspect A closed-form, Aspect B + C empirically anchored)
**Date:** 2026-05-25 (placeholder created during Task 5 registration; full proof to be written in Task 7)
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md](PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md) (closes Z-dephasing Π²-even cells)
- [PROOF_F108_PART2_PI2X_EVEN_ALWAYS_PALINDROMIC.md](PROOF_F108_PART2_PI2X_EVEN_ALWAYS_PALINDROMIC.md) (closes X-dephasing Π²-even cells)
- [PROOF_F108_PART3_PI2Y_EVEN_ALWAYS_PALINDROMIC.md](PROOF_F108_PART3_PI2Y_EVEN_ALWAYS_PALINDROMIC.md) (closes Y-dephasing Π²-even cells)
- [PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md](PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md) (Mother sector truly y_par = 0)
- [PROOF_F109_MOTHER_SOFT_Y_PARITY_ONE_PURITY.md](PROOF_F109_MOTHER_SOFT_Y_PARITY_ONE_PURITY.md) (Mother sector soft y_par = 1)

**Statement (Theorem F110, three aspects):**

**Aspect A (closed-form):** Under any single-letter dephase channel D ∈ {X, Y, Z}, every F87-hard Pauli pair lives in the diagonal Klein cell whose Klein index matches the dephase letter:
- Z-dephase: hard only in Klein (0, 1)
- X-dephase: hard only in Klein (1, 0)
- Y-dephase: hard only in Klein (1, 1)

All other cells admit no F87-hard pairs (forced by F108 Part 1+2+3 closing Π²-D-even cells and F107+F109 closing the Mother sector Klein (0, 0)).

**Aspect B (empirical Y-inversion):** Within the diagonal Klein cell, the dominant y_par equals y_par(dephase letter):
- Z-dephase, X-dephase: dominant y_par = 0
- Y-dephase: dominant y_par = 1 (inverted)

**Aspect C (empirical k-purity sharpening):** The dominance is k-dependent:
- k = 3, N = 4 (F103): 42:8 biased per diagonal cell
- k = 3, N = 5 (F105): identical 42:8 (N-stable)
- k = 4, N = 4 (F106): 228:0 fully pure (Y-inversion preserved)

Closed-form derivation of the exact 42:8 and 228:0 ratios is open (F103 Section 5).

## Proof (placeholder for Task 7)

Detailed proof of Aspect A by exclusion across F108 Part 1+2+3 + F107 + F109; structural reading of Aspect B; empirical anchors for Aspect C will be written in Task 7 of the F110 implementation plan.
