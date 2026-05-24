# PROOF F109: Mother Sector Soft is y_par = 1 Pure (All Dephase Letters)

**Status:** Tier 1 derived (closed-form modulo one named empirical dependency)
**Date:** 2026-05-24
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md](PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md) (per-dephase truly criteria)
- [PROOF_F85_KBODY_GENERALIZATION.md](PROOF_F85_KBODY_GENERALIZATION.md) (k-body truly criterion under Z-dephasing)
- **Open: Π²-even-non-truly is never hard (F108 Part 1, deferred).** Empirically verified across F103/F105/F106 (5346+ pairs, zero hard observed in any Π²-even cell). F109 inherits this dependency at Step 5.

**Statement (Theorem F109):** Under any single-letter dephase channel (Z, X, or Y), every Pauli pair classified as soft and located in the Mother sector Klein (0, 0) has shared y_par = 1.

Empirically confirmed across F103 (mother soft (0, 21) ×3 dephase letters), F105 (same), F106 (mother soft (0, 300) ×3 — sharpened at k=4). Total: 1026 mother-soft classifications, all y_par=1, zero y_par=0.

## Proof

### Step 1: Klein (0, 0) forces all three (#X, #Y, #Z) to share the same parity

Klein (0, 0) means bit_a = 0 AND bit_b = 0, i.e., #X + #Y ≡ 0 (mod 2) AND #Y + #Z ≡ 0 (mod 2).

From #X + #Y even AND #Y + #Z even: subtracting gives #X − #Z ≡ 0 (mod 2), so #X and #Z have the same parity. Combined with #X + #Y even: #X and #Y same parity. Hence **#X, #Y, #Z all share the same parity** (all even or all odd).

### Step 2: Per-dephase F107 truly criteria collapse on Klein (0, 0) to "all three even"

Per F107, the per-dephase F87 truly criterion is:
- Z-dephase: #Y even AND #Z even
- X-dephase: #X even AND #Y even
- Y-dephase: #Y even AND #Z even

Combined with Step 1's same-parity constraint, each dephase's truly criterion forces all three to be even:
- Z-dephase truly at Klein (0, 0): #Y even AND #Z even AND (all same parity) ⟹ all three even.
- X-dephase truly: #X even AND #Y even AND (all same parity) ⟹ all even.
- Y-dephase truly: #Y even AND #Z even AND (all same parity) ⟹ all even.

So Klein (0, 0) truly under any dephase = all #X, #Y, #Z even. These terms have y_par = #Y mod 2 = 0 (consistent with F107).

### Step 3: Klein (0, 0) non-truly = all three counts odd

Step 1's same-parity constraint gives two cases:
- All three even ⟹ truly (Step 2).
- All three odd ⟹ NOT truly under any dephase.

So Klein (0, 0) non-truly = #X, #Y, #Z all odd.

### Step 4: Klein (0, 0) is Π²-even under every dephase letter

The Π² eigenvalue per dephase (per PiOperator.SquaredEigenvalue):
- Z-dephase: Π²_Z parity = bit_b = 0 ⟹ Π²-EVEN
- X-dephase: Π²_X parity = bit_a = 0 ⟹ Π²-EVEN
- Y-dephase: Π²_Y parity = bit_b = 0 ⟹ Π²-EVEN

Klein (0, 0) is the only cell that is Π²-even under all three dephase letters simultaneously.

### Step 5 (depends on the open Π²-even-soft observation): Klein (0, 0) non-truly is SOFT

Π²-even non-truly pairs have M_anti = L_{H_odd} = 0 (per F81 / F85: no Π²-odd Hamiltonian content). Empirically across F103 (N=4 k=3), F105 (N=5 k=3), F106 (N=4 k=4): every Π²-even non-truly pair is classified soft, zero hard observed across 5346+ pairs.

The closed-form proof of "Π²-even non-truly ⟹ soft" requires a block-restricted palindrome lemma on the Π²-eigenspace decomposition of L (Stage F108 Part 1; currently open). F109 inherits this dependency. When F108 Part 1 closes, F109 becomes fully closed-form Tier1Derived.

Under this empirically-overwhelming hypothesis: Klein (0, 0) non-truly pairs are SOFT (not hard).

### Step 6: All Klein (0, 0) soft pairs have y_par = 1

Combining Steps 3 and 5: Klein (0, 0) soft pairs have terms with #X, #Y, #Z all odd. y_par = #Y mod 2 = 1 for each term. For y_par-homogeneous pairs (both terms share y_par): pair y_par = 1.

∎

## Empirical confirmation

| Anchor | Mother soft cells | (y_par=0, y_par=1) per cell | Match Step 6? |
|--------|--------------------|------------------------------|---------------|
| F103 (N=4 k=3) | 3 (Z, X, Y dephase) | each (0, 21) | ✓ |
| F105 (N=5 k=3) | 3 | each (0, 21) | ✓ |
| F106 (N=4 k=4) | 3 | each (0, 300) | ✓ |

Total: 1026 mother-soft classifications, all y_par=1, zero y_par=0. F109 explains this bit-exactly.

## Cross-letter spot-check: enumerate Klein (0, 0) non-truly k=3 terms

Per Step 3, Klein (0, 0) non-truly k=3 terms have #X, #Y, #Z all odd, summing to ≤ 3 (since k=3 letter sequence length). The only triple with all-odd and sum ≤ 3 is (1, 1, 1), so terms are permutations of XYZ with one of each non-I letter.

Number of such terms: 3! = 6 (XYZ, XZY, YXZ, YZX, ZXY, ZYX).

Number of y_par-homogeneous unordered pairs (including self-pairs): 6·7/2 = 21. **Matches F103/F105 mother soft = (0, 21) per dephase letter exactly.** ✓

For k=4 letter sequences (N=4 enumeration): Klein (0, 0) non-truly = #X, #Y, #Z all odd, sum ≤ 4. Only (1, 1, 1) with #I = 1: 4!/(1!1!1!1!) = 24 letter sequences. Unordered pairs with self: 24·25/2 = 300. **Matches F106 mother soft = (0, 300) per dephase letter exactly.** ✓

## Open

- **F108 Part 1** (the Π²-even-soft step): the block-restricted palindrome lemma on Π²-eigenspace decomposition of L. F109's full Tier1Derived status pending.
- **F110** (hard cells y_par-pure with Y-inversion): the harder dephase-letter-specific analysis on the diagonal-Klein cells. F107 + F109 do not directly attack F110.
