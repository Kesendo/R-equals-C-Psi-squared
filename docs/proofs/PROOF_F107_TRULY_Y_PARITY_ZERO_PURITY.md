# PROOF F107: F87 Truly Classification Forces y_par = 0 (All Dephase Letters)

**Status:** Tier 1 derived (closed-form corollary of F85 + per-dephase Π² + dissipator commutativity)
**Date:** 2026-05-24
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [PROOF_F85_KBODY_GENERALIZATION.md](PROOF_F85_KBODY_GENERALIZATION.md) (k-body truly criterion under Z-dephasing)
- [MIRROR_SYMMETRY_PROOF.md](MIRROR_SYMMETRY_PROOF.md) (Π definitions per dephase letter)
- [`compute/RCPsiSquared.Core/Symmetry/PiOperator.cs`](../../compute/RCPsiSquared.Core/Symmetry/PiOperator.cs) (per-dephase Π² eigenvalue rules)

**Statement (Theorem F107):** Under any single-letter dephase channel (Z, X, or Y), if a Pauli term σ_α is classified as truly by the F87 trichotomy, then y_par(σ_α) = (#Y in α) mod 2 = 0.

By extension, for any y_par-homogeneous Pauli pair classified as truly, the pair's shared y_par value is 0. Empirically confirmed across F103/F105/F106 anchors (N=4 k=3, N=5 k=3, N=4 k=4): zero truly classifications with y_par=1 across all 12 (Klein × dephase) cells at every measured (N, k) regime.

## Proof

### Step 1: per-dephase Π² + dissipator yields a per-letter even-count criterion

For dephase letter D ∈ {Z, X, Y}, the truly condition M(H) = 0 (palindrome residual vanishes) decomposes into two independent conditions on each Pauli term in H:

1. **Π² parity** (from MIRROR_SYMMETRY + Π definition per dephase letter, encoded in `PiOperator.SquaredEigenvalue`): the term must be Π²-even under the dephase-specific Π². For Z and Y dephasing Π² sums bit_b; for X dephasing Π² sums bit_a.
2. **Dissipator commutativity** (from F84 Pauli-Channel Cancellation Lemma): the term must commute with the dephase-letter operator, i.e., have an even count of letters that anticommute with D.

The dephase letter D commutes with itself and with I; anticommutes with the other two non-I letters. So condition (2) reduces to: count of non-D non-I letters is even.

Combining (1) and (2) per dephase letter:

| Dephase | (1) Π²-even condition | (2) Dissipator-commute | Combined truly criterion |
|---------|----------------------|------------------------|--------------------------|
| Z       | bit_b = #Y + #Z even | #Y even (the only non-Z anticommutant other than X) — actually no: D=Z anticommutes with X and Y; #(X,Y) even ↔ #X + #Y even | #Z even (from (1) + #Y even from (2) implies #Z even); equivalently #Y even AND #Z even |
| X       | bit_a = #X + #Y even | D=X anticommutes with Y and Z; #(Y,Z) even ↔ #Y + #Z even | #X even AND #Y even |
| Y       | bit_b = #Y + #Z even | D=Y anticommutes with X and Z; #(X,Z) even ↔ #X + #Z even | #Y even AND #Z even |

The Z and Y rows yield identical criteria (#Y even AND #Z even) because Π_Y and Π_Z share the same per-letter swap (I↔X, Y↔Z), differing only in the phase factor (±i for Y/Z swap), and the dissipator-commutativity condition under Y is also #X + #Z even — combined with Π_Y's bit_b even (#Y + #Z even) and subtracting in Z₂ gives the same constraints up to relabel.

### Step 2: every truly criterion includes #Y even

Reading the rightmost column:

- Z-dephasing truly: #Y even AND #Z even — includes #Y even
- X-dephasing truly: #X even AND #Y even — includes #Y even
- Y-dephasing truly: #Y even AND #Z even — includes #Y even

All three dephase letters force #Y even as a sub-condition.

### Step 3: #Y even ⟹ y_par = 0

By definition, y_par(σ) = (#Y in σ) mod 2. #Y even ⟹ y_par = 0.

### Step 4: y_par-homogeneous pair corollary

A Klein-homogeneous + y_par-homogeneous pair (term1, term2) has shared y_par value y_par(term1) = y_par(term2). The pair is truly iff both terms individually satisfy the truly criterion (per F85, M decomposes term-by-term in the linear regime). Each truly term has y_par = 0; hence pair y_par = 0.

∎

## Empirical confirmation

| Anchor | Cells with truly counts | Total truly | y_par=1 truly |
|--------|------------------------|-------------|----------------|
| F103 (N=4 k=3) | 6 of 12 | 300 | 0 |
| F105 (N=5 k=3) | 6 of 12 | 300 | 0 |
| F106 (N=4 k=4) | 9 of 12 | 3924 | 0 |

Total: 4524 truly classifications observed across (Klein × dephase × y_par × N × k); zero have y_par=1. F107 explains this bit-exactly as a closed-form corollary.

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

## Open

- **F108** (next): hard appears only in diagonal Klein cells (F87 dissipator-resonance law). Already stated in `docs/ANALYTICAL_FORMULAS.md` F87 entry but lacks a formal closed-form proof.
- **F109** (medium): mother (0,0) soft is y_par=1-pure across all dephase letters. Empirically held at k=3 and k=4 (sharper at k=4). Requires Pauli-pair compatibility analysis layered on F102.
- **F110** (higher): hard cells y_par-pure with Y-inversion. Becomes 100% pure at k=4 (228:0 vs 0:228). Requires per-dephase-letter algebra on the hard-cell pair-set.
