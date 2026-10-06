# PROOF F109: Mother Sector Soft is y_par = 1 Pure for y_par-Homogeneous Pairs (All Dephase Letters)

**Status:** Tier 1 derived (Steps 1-4 and 6 count letters; Step 5 is the colouring of the mother cell, with F158)
**Date:** 2026-05-24
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md](PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md) (per-dephase truly criteria)
- [PROOF_F85_KBODY_GENERALIZATION.md](PROOF_F85_KBODY_GENERALIZATION.md) (k-body truly criterion under Z-dephasing)
- [PROOF_PALINDROME_TWO_END_COUNT.md](PROOF_PALINDROME_TWO_END_COUNT.md) (F158 §(e): an invertible operator that commutes with H and anticommutes with every jump reflects L to −L† − 2σ, and the spectrum pairs about −σ)
- [The palindrome as a colouring](../../experiments/THE_PALINDROME_AS_A_COLOURING.md) (a lit Pauli string commuting with H is such an operator)
- [PROOF_F103 §8](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) (every non-diagonal Klein cell is coloured)

## Abstract

F107 closed one of the two purity statements that the empirical Z₂³ refinement made obvious: F87-truly is always y-parity zero. F109 closes the other one. In the Mother sector (the Klein (0,0) cell where bit_a and bit_b both vanish), a soft pair of one y-parity is always of y-parity ONE, never zero. The two purity statements are mirror images on the y-parity axis: truly forces an even number of Y letters on every term, and a Mother-cell term that is not truly carries an odd number.

The proof is six steps, with a clean structural story. The Mother sector forces all three Pauli-letter counts (#X, #Y, #Z) to share the same parity (Klein (0,0) means bit_a and bit_b both even, which forces the parities to match). Combined with the per-dephase truly criterion (from F107), the all-even case is exactly truly, leaving the all-odd case as everything that is not truly in the Mother sector. The all-odd case has #Y odd, hence y-parity one.

What remains is whether the all-odd non-truly case is F87-soft rather than F87-hard, and the answer is a colouring. A Mother-sector string commutes with every letter string X^⊗N, Y^⊗N, Z^⊗N. For dephase letter D, the letter string of the canonical mirror's flip letter (X for Z- and Y-dephasing, Z for X-dephasing) anticommutes with every jump, so right multiplication by it leaves the Hamiltonian part of L alone and reflects the dissipator, and F158 turns that into a spectrum palindromic about −σ. The all-odd case is soft at every body count, every coupling and every per-site rate.

The diagnostic upshot is that the Mother sector's strings separate cleanly on the y-parity axis, the truly ones at zero and the others at one, so its y_par-homogeneous truly and soft pairs are pure in both directions. The 1026 mother-soft classifications observed across F103+F105+F106 (63 + 63 + 900 across three dephase letters and three (N, k) anchors) all sit on y-parity one with zero exceptions. F107 + F109 together close both purity statements of the y-parity-axis classification; F110 and F111 handle the harder cells where the pattern is more nuanced.

**Statement (Theorem F109):** Under any single-letter dephase channel (Z, X, or Y), every y_par-homogeneous Pauli pair classified as soft and located in the Mother sector Klein (0, 0) has shared y_par = 1.

Empirically confirmed across F103 (mother soft (0, 21) ×3 dephase letters), F105 (same), F106 (mother soft (0, 300) ×3, sharpened at k=4). Total: 1026 mother-soft classifications, all y_par=1, zero y_par=0.

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

### Step 5: Klein (0, 0) non-truly is SOFT, at every body count

A string commutes with X^⊗N when #Y + #Z is even, with Z^⊗N when #X + #Y is even and with Y^⊗N when #X + #Z is even. In Klein (0, 0) all three counts share one parity (Step 1), so every Mother-sector string commutes with all three letter strings, and so does every Hamiltonian H built from them.

Let A be the letter the canonical mirror flips by: X for Z- and Y-dephasing, Z for X-dephasing. A anticommutes with the dephase letter, so F = A^⊗N anticommutes with every jump. Right multiplication R(ρ) = ρ·F then leaves the Hamiltonian part alone, −i[H, ρF] = (−i[H, ρ])·F since F commutes with H, and reflects each dephasing term, γ_l(D_l·ρF·D_l − ρF) = −γ_l(D_l·ρ·D_l + ρ)·F since F anticommutes with D_l. Together

    R · L · R⁻¹ = −L† − 2σ,    σ = Σ_l γ_l.

This is the colouring of [the palindrome as a colouring](../../experiments/THE_PALINDROME_AS_A_COLOURING.md), a lit string commuting with H, and the sufficiency step of [F158](PROOF_PALINDROME_TWO_END_COUNT.md) §(e) turns it into a palindrome: spec(L) = {−λ̄ − 2σ : λ ∈ spec(L)}, and since L preserves hermiticity its spectrum is closed under conjugation, so spec(L) = {−λ − 2σ : λ ∈ spec(L)}. That holds for every Hamiltonian of the cell, every coupling and every per-site rate, so Klein (0, 0) holds no hard pair, and its non-truly pairs are soft. For Z-dephasing the map R is the factor R of Π_Z = R·D in [the Π factorization](PROOF_PI_FACTORS_AS_R_TIMES_D.md) §3 (D there is the transpose), which keeps the Hamiltonian commutator exactly when H commutes with X^⊗N. [F103 §8](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) settles the other non-diagonal Klein cells the same way.

No mirror of the F108 family reaches these strings: its three per-site mirrors carry L_σ to −L_σ for exactly the Π²-D-even strings σ of even weight (an even number of non-identity letters), and an all-odd string has odd weight. The canonical Π_D has M ≠ 0 on them by Step 3. [`f107_f110_route_gate.py`](../../simulations/f107_f110_route_gate.py) checks R·L·R⁻¹ = −L† − 2σ exactly at N = 4 on random Π²-D-even Hamiltonians of any body count, mother strings among them, with random per-site rates, and checks that none of the 24 mother non-truly strings at N = 4 is reached by an F108 mirror.

### Step 6: All Klein (0, 0) soft pairs have y_par = 1

A soft pair is not truly, so at least one of its terms is non-truly, hence all-odd by Step 3, with y_par = #Y mod 2 = 1. The pair is y_par-homogeneous, so its other term has y_par = 1 as well, and in Klein (0, 0) that makes it all-odd too: pair y_par = 1. Step 5 adds that every non-truly pair of the cell is soft, so the cell's soft pairs are exactly its non-truly pairs.

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

## Sibling y_par-axis claims

F107 (TrulyYParityZeroPurity, Tier1Derived); F110 (HardCellYInversionPattern, Tier1Derived); F111 (HardCellPureDTemplate, Tier1Derived).

## Open

- k ≥ 5 empirical confirmation of F103/F106 pattern stability beyond N=4, for the soft, truly and off-diagonal cells (the hard diagonal cell is read at k = 4, 5 by `DiagonalCellComplementConnectionTests`).
- Hardware QPU confirmation at k ≥ 3 (no F87 QPU confirmations exist beyond Marrakesh k=2).
