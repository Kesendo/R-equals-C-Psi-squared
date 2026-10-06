# PROOF F111: Hard Cell Pure-D Template Rule (Tier1Derived, promoted 2026-06-10)

**Status:** Tier 1 derived (promoted 2026-06-10; typed `HardCellPureDTemplate`, Tier1Derived). Empirical anchor across 3 dephase letters at N=4 k=4. **Subclaim (d) Mixed+Mixed = soft closed modulo M = −2i(H⊗I) on 2026-05-30** via the bipartite-chirality mechanism of [PROOF_F103 §7.4](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) (the operator is a *second* mirror, the chiral K composed after Π, not the better Π the three Task-1 paths sought); subclaim (b) follows likewise. The **hard-direction converse behind subclaims (a)/(c) closed 2026-06-10** (WindowedConverseAllGammaClaim, Pascal-Gram positivity F117, no residual), which is what promoted F111 to Tier1Derived. See the Open section.
**Date:** 2026-05-25
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [PROOF_F110_HARD_CELL_Y_INVERSION.md](PROOF_F110_HARD_CELL_Y_INVERSION.md) (Aspect A, which leaves the diagonal cell as the place of the hardness of Klein-homogeneous pairs; F111 sharpens F110 Aspect B there and implies it as corollary at k = N = 4)
- [PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md](PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md) (per-dephase truly criterion; pure-D template's #Y count derivation)
- [PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md](PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md) (Π_5bilinear for Z-deph; failed candidate for off-y_par palindromization in Task 1)
- [PROOF_F108_PART2_PI2X_EVEN_ALWAYS_PALINDROMIC.md](PROOF_F108_PART2_PI2X_EVEN_ALWAYS_PALINDROMIC.md) (Π_5bilinear for X-deph; likewise a failed candidate, Task 1 Path 2)
- [PROOF_F108_PART3_PI2Y_EVEN_ALWAYS_PALINDROMIC.md](PROOF_F108_PART3_PI2Y_EVEN_ALWAYS_PALINDROMIC.md) (Π_5bilinear for Y-deph; likewise a failed candidate, Task 1 Path 2)
- Theorem 2 of the [complement connection](PROOF_PALINDROME_COMPLEMENT_CONNECTION.md) (the exact palindrome criterion that decides (a) to (d) and the rule beyond k = N = 4)
- F106 N = 4 k = 4 empirical anchor (compute/RCPsiSquared.Core/Symmetry/F87Z2CubedRefinementN4K4.cs); the F87 dissipator-resonance law (compute/RCPsiSquared.Diagnostics/F87/DissipatorResonanceLaw.cs) is the N = 4, k = 3 census of where the hardness of Klein-homogeneous pairs lives

## Abstract

F110 left two pieces open. The first was the Y-inversion ratio at k = 3 (the empirical 42:8 split in the diagonal hard cells). The second was the sharpened version of the same observation at k = 4, where the split goes fully pure (228:0). F111 closes the k = 4 case with a structural rule and proves the corresponding F110 Aspect B as an immediate corollary.

The rule is the "Pure-D Template Rule". At k = N = 4 in the diagonal Klein cell for dephase letter D, a Pauli pair is F87-hard if and only if at least one of the two terms is a "pure-D template": a length-4 Pauli string built only from the letter D and the identity, with no other non-identity letter mixed in. So the eight strings with one or three D letters, {DIII, IDII, IIDI, IIID, DDDI, DDID, DIDD, IDDD} (an even D count would leave the diagonal cell), classifies the entire hard-pair population at k = N = 4 in the diagonal cell of D.

The empirical numerics fit cleanly. 8 pure-D templates × 9/2 self-pair counts = 36 pure-pure pairs, all hard. 8 pure-D templates × 24 mixed templates = 192 pure-mixed pairs, all hard. The remaining 24 × 25/2 = 300 mixed-mixed pairs are all soft. 36 + 192 + 0 = 228, matching the F106 anchor exactly across all three dephase letters. The F110 Aspect B corollary then drops out: pure-D templates carry y-parity equal to the y-parity of D itself, so every hard pair inherits that y-parity, which is the 228:0 split with Y-inversion.

F111 shipped Tier 1 Candidate on 2026-05-25 because the "Mixed+Mixed pair is soft" half of the rule (the 0 in 36+192+0) had only an empirical anchor across 300 pairs and no operator-level construction: three derivation paths were attempted and exhausted (per-site M^N tensor-product search, F108 Π_5bilinear extended action, Q_V × Π composition). It was promoted to Tier 1 Derived on 2026-06-10: subclaim (d) closed modulo M = −2i(H⊗I) via the chiral-K route (PROOF_F103 §7.4, 2026-05-30), and the hard-direction converse behind subclaims (a)/(c) closed via the windowed all-γ theorem (WindowedConverseAllGammaClaim, Pascal-Gram positivity F117, no residual). The similarity transformation that realizes the spectrum-level palindromy is that second mirror (the chiral K composed after Π), not the single better Π the Task-1 paths sought.

The diagnostic upshot is that F111 sharpens F110 Aspect B at the specific k = N = 4 anchor: from "Y-inversion empirical 228:0" to a structural per-pair rule. With F107, F109 and F110 it gives the y-parity signature of the truly, Mother-soft and hard pairs; the soft pairs of the other cells are F103's (§3.3, §8).

## Statement (Theorem F111)

At k = N = 4 in the diagonal Klein cell (D.BitA(), D.BitB()) for dephase letter D ∈ {Z, X, Y}, a Pauli pair (P, Q) is F87-hard if and only if at least one of P, Q is a "pure-D template" (a length-4 Pauli string composed only of D and I letters, no other non-I Pauli letter).

**Corollary (F110 Aspect B at k = 4):** Pure-D templates have y_par = y_par(D) by construction (templates contain only D and I; only Y has #Y = 1 mod 2 of itself, so y_par(pure-Z) = y_par(pure-X) = 0 and y_par(pure-Y) = 1). Therefore at k = N = 4 in the diagonal cell, every F87-hard pair has y_par(pair) = y_par(D). This is the F106 N = 4 k = 4 228:0 split across all 3 dephase letters.

Originally Tier1Candidate, not Tier1Derived: subclaim (d) (Mixed, Mixed) = soft at k = N = 4 lacked an operator-level closed-form derivation (see Section 3 below). Since then: (d) closed modulo M via PROOF_F103 §7.4 (2026-05-30), and the hard-direction converse behind (a)/(c) closed 2026-06-10 (WindowedConverseAllGammaClaim, Pascal-Gram positivity F117, no residual), promoting F111 to Tier1Derived.

## Empirical anchor

F106 N = 4 k = 4 enumeration (per `compute/RCPsiSquared.Core/Symmetry/F87Z2CubedRefinementN4K4.cs` frozen record and `simulations/results/f87_z2cubed_split_n4_k4_counts.json`):

| Klein cell | Dephase D | Hard count (y_par=0, y_par=1) | Pure-D Template Rule prediction |
|------------|-----------|-------------------------------|----------------------------------|
| (0, 1)     | Z         | (228, 0)                      | 228 hard (36 Pure-Pure + 192 Pure-Mixed); ✓ |
| (1, 0)     | X         | (228, 0)                      | 228 hard; ✓ |
| (1, 1)     | Y         | (0, 228)                      | 228 hard (Y-inversion); ✓ |

Independent verification: `simulations/f111_pair_off_ypar_verify.py` classifies all 528 pairs per on-y_par sub-cell of the diagonal cell × 3 dephases = 1584 pair classifications, all matching the Pure-D Template Rule with zero exceptions. (Each diagonal cell decomposes into an on-y_par sub-cell of 528 pairs and an off-y_par sub-cell of 528 pairs; the off-y_par sub-cell contains zero pure-D templates by Step 2 below, and the rule predicts zero hard pairs there, also empirically verified, so the rule covers the full diagonal cell.)

## Structural decomposition (per on-y_par sub-cell of the diagonal cell at k = N = 4)

The on-y_par sub-cell of the diagonal cell contains 32 length-4 templates (templates with the same y_par as D and the matching Klein index). The 528 = 32·33/2 unordered pairs (with self) on these templates decompose by template-membership:

| Pair class | Count formula | Count | F87 status |
|------------|---------------|-------|------------|
| Pure-Pure (both terms pure-D) | 8·9/2 | 36 | HARD (per subclaim a) |
| Pure-Mixed (one pure-D, one mixed) | 8·24 | 192 | HARD (per subclaim c) |
| Mixed-Mixed (both terms mixed) | 24·25/2 | 300 | SOFT (per subclaim d) |
| **Total** | | **528** | 36 + 192 = 228 hard |

The 36 + 192 + 0 = 228 count matches F106 exactly across all 3 dephase letters. The off-y_par sub-cell of the diagonal cell (another 528 pairs) contains no pure-D templates by construction (pure-D template's y_par equals y_par(D) by Step 2 below); under the Pure-D Template Rule it contributes zero hard pairs, also matching F106's empirical zero in the off-y_par half of the (228, 0) split.

Counting the Pure-D templates per on-y_par sub-cell: at k = 4 in the diagonal cell with bit_b (Z, Y) = 1 or bit_a (X) = 1, pure-D templates have either #D = 1, #I = 3 (4 placements) or #D = 3, #I = 1 (4 placements). Total: 8 pure-D templates in the on-y_par sub-cell, 0 in the off-y_par sub-cell.

## Proof

### Step 1: Pure-D templates lie in the diagonal cell

A pure-D template at k = 4 has #D = #D, #I = 4 − #D, and all other counts zero. The Klein index is computed from (bit_a, bit_b) = (#X + #Y mod 2, #Y + #Z mod 2).

For D = Z (pure-Z template): #X = #Y = 0, #Z = #D. Klein index = (0, #Z mod 2). The diagonal cell for Z-deph is (0, 1), so pure-Z templates with #Z odd are in the diagonal cell. At k = 4 the odd #Z values are 1 and 3. Count: 4 (placements of single Z) + 4 (placements of 3 Z's) = 8.

For D = X: analogous, with bit_a flipped instead of bit_b. Diagonal cell for X-deph is (1, 0); pure-X templates have Klein (#X mod 2, 0). With #X odd: 8 templates.

For D = Y: pure-Y templates have #Y = #D, #X = #Z = 0. Klein index = (#Y mod 2, #Y mod 2). For #Y odd: Klein (1, 1) = Y-deph diagonal. 8 templates.

In all three cases, the count of pure-D templates per diagonal cell at k = N = 4 is 8 (4 with one D + 4 with three D's). ✓

### Step 2: Pure-D template's y_par equals y_par(D)

Pure-D templates have #Y = #D (if D = Y) or #Y = 0 (if D = Z or X).

- D = Z: #Y = 0 ⟹ y_par = 0 = y_par(Z) (since y_par(Z) = Z.BitA() AND Z.BitB() = 0 AND 1 = 0). ✓
- D = X: #Y = 0 ⟹ y_par = 0 = y_par(X) (= 1 AND 0 = 0). ✓
- D = Y: #Y = #D = 1 or 3 (odd) ⟹ y_par = 1 = y_par(Y) (= 1 AND 1 = 1). ✓

Therefore: pure-D templates have y_par = y_par(D), and the Pure-D Template Rule directly implies F110 Aspect B at k = N = 4. ✓

### Step 3: Subclaim (a): pure-D single-term H is F87-hard

Pure-D template H is built from D and I letters only. The dephase letter D commutes with itself: [D, D] = 0 per-site. Therefore the Z-dephasing dissipator's per-site action D[D_l] · σ = γ_l · (D_l σ D_l − σ) vanishes when σ contains only D and I at site l (since D_l commutes with D and with I, so D_l σ D_l = σ).

For a pure-D template H, every Pauli string in the Hamiltonian commutes with every per-site D. Therefore [D[D_l], L_H] = 0 for all l. The Lindbladian decomposes as L = L_H + L_D with [L_H, L_D] = 0.

Both summands carry their own reflection symmetry. L_H is the commutator superoperator (−i[H, ·]) of a Hermitian H, so its spectrum is {−i(E_k − E_j)}, pure imaginary and symmetric about 0. L_D is symmetric about −σ: in the Pauli basis it is diagonal with eigenvalue −2·Σ_{l ∈ A} γ_l on the string whose D-anticommuting set is A, and the set-complement A ↦ {1..N} \ A is an involution on the basis carrying that eigenvalue to −2σ + 2·Σ_{l ∈ A} γ_l. This holds for arbitrary site-dependent γ_l, not only uniform γ.

So the failure is **not** in either summand. It is that no single involution realizes both reflections at once on the shared eigenbasis: L_H and L_D commute, but the complement map that reflects L_D's spectrum about −σ does not preserve L_H's eigenvalue on the same basis vector, so the reflected pair is not in the spectrum of the sum. The (multiset) sum of the two spectra is therefore the wrong object: for commuting operators the eigenvalues of L = L_H + L_D are the sums taken **along shared eigenvectors**, not the free sum of the two lists.

That failure has a single-eigenvalue witness. Every Liouvillian carries the steady eigenvalue 0, and the reflection λ ↦ −2σ − λ sends it to the anti-steady eigenvalue −2σ, so a palindromic spectrum must contain −2σ. Measured over the full k = N = 4 diagonal cell for all three dephase letters, with site-dependent γ = (0.31, 0.17, 0.53, 0.11) as well as uniform γ: **−2σ is absent for all 8 pure-D templates and present for all 24 mixed ones**, in every cell, and the split agrees eigenvalue-for-eigenvalue with the independent Π-residual classifier's 8 hard / 24 soft. Hence the F87 spectrum-pairing condition fails for pure-D H, and pure-D H is F87-hard.

Theorem 2 of the [complement connection](PROOF_PALINDROME_COMPLEMENT_CONNECTION.md) gives the *why* from the letter structure. Rotated to Z dephasing, a pure-D template is diagonal with an odd number of Z letters: it has no hopping edges, every bitstring is its own component, and its diagonal flips sign under the complement at every bitstring, so no component is good, the far count dim ker(L + 2σ) is 0, and −2σ is absent. A mixed template has no diagonal; X⊗N anticommutes with it, and its single flip mask makes the hopping graph bipartite, so every component is good. The hard direction is closed independently by the windowed all-γ converse (WindowedConverseAllGammaClaim, Pascal-Gram positivity F117, no residual, 2026-06-10), which is what promoted F111 to Tier1Derived; with Theorem 2 this step is a second, exact route at full support. Both the classification and the −2σ witness are computed by `simulations/f111_spec_palindrome_single_term.py`.

### Step 4: Subclaims (b), (c), (d)

**Subclaim (b):** Mixed single-term H (contains a non-D non-I letter) at k = N = 4 in the diagonal cell is F87-soft. Verified for all 24 mixed templates per diagonal cell; the mechanism is Theorem 2 as in Step 3 (one flip mask, a bipartite hopping graph), and PROOF_F103 §7.4's chiral K.

**Subclaim (c):** Pair (Pure-D, Mixed) H at k = N = 4 is F87-hard. Verified for all 192 Pure-Mixed pairs per diagonal cell; the mechanism is condition (ii) of Theorem 2 (see Open): the pure-D template's diagonal flips under the complement and the mixed template has none.

**Subclaim (d):** Pair (Mixed, Mixed) H at k = N = 4 is F87-soft. Verified for all 300 Mixed-Mixed pairs per diagonal cell, and derived by Theorem 2 of the complement connection: X⊗N anticommutes with every Mixed string of the cell, so H_x̄ȳ = −H_xy, condition (ii) holds (no diagonal) and (iii) asks for a bipartite hopping graph; at full support each template is one flip mask, and two nonzero masks always admit a linear φ equal to 1 on both (PROOF_F103 §7.4's chiral K), so the graph is bipartite and the pair palindromic; soft rather than truly because no diagonal-cell string is truly (its Y/Z count is odd) and a pair is truly only when both terms are ([F107](PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md)). For the y-parity-homogeneous pairs counted here the operator that carries it acts site by site: ρ ↦ Π(ρ)·K at y_par 0 and ρ ↦ K·Π(ρ) at y_par 1, with the chiral K = ⊗Z_l^{φ_l} and φ depending on the site and on the pair (PROOF_F103 §7.4). A y-parity-mismatched pair is palindromic as well, by the same graph. The three operator searches below did not reach it: a site-uniform M^⊗N (Path 1: 512 phase variants × 2 dissipator-valid permutations per dephase, no winners), F108's Π_5bilinear (Path 2: residual 32 uniformly on off-y_par single-term H) and global Q_V × Π compositions (Path 3: no zero-residual hits). See `simulations/f111_path1_operator_search.py` (Path 1), `simulations/f111_path2_pi5bilinear_test.py` (Path 2) and `simulations/f111_combined_operator_search.py` (Path 3) for the verification scripts.

The spectrum-level palindromy IS realized by a similarity S·L·S⁻¹ = −L − 2σ: with Z-dephasing the system is in F158's class, where a palindromic spectrum makes the two end counts agree and [the two-end count](PROOF_PALINDROME_TWO_END_COUNT.md) §(f8) turns that into a similarity of the whole Jordan structure; for the y-parity-homogeneous Mixed+Mixed pairs one such S is the per-site map above.

### Step 5: F111 statement follows from subclaims (a)-(d)

On the on-y_par sub-cell of the diagonal cell, pair (P, Q) is F87-hard iff at least one of P, Q is a pure-D template:
- Both pure-D: hard, by condition (ii) of Theorem 2 (see Open).
- One pure-D, one mixed: hard via (c).
- Both mixed: SOFT via (d), so NOT hard.

On the off-y_par sub-cell of the diagonal cell, no template is pure-D (Step 2: pure-D ⟹ y_par = y_par(D), which excludes the off-y_par sub-cell by construction), so the "at least one pure-D" condition is vacuously false everywhere; the rule predicts zero hard pairs in the off-y_par sub-cell, matching F106's empirical zero in the off-y_par half of the (228, 0) split.

Combining both sub-cells, the Pure-D Template Rule covers the full diagonal cell. Combined with Step 2 (pure-D ⟹ y_par = y_par(D)), F110 Aspect B at k = N = 4 follows: every F87-hard pair has y_par = y_par(D). ∎

## Path 1 / Path 2 / Path 3 derivation attempts (Task 1, BLOCKED)

The original Task 1 goal was to derive F111 as Tier1Derived via closed-form palindromization of the off-y_par sub-sector. Three paths attempted:

**Path 1 (per-site M⊗N tensor product):** Brute force over candidate per-site M operators that anti-commute with all k = 4 off-y_par(D) bilinears in the diagonal cell AND conjugate the D-dephasing dissipator correctly. Scan: 512 phase variants × 2 dissipator-valid letter permutations per dephase. Result: **zero winners.** None of the 512 scanned site-uniform M^⊗N achieves operator-level palindrome on any off-y_par single-term H.

**Path 2 (existing F108 Π_5bilinear extended action):** Tested Pi_5bilinear (Z, X, Y variants) on all 32 templates in the diagonal cell per dephase, classified by y_par. Result: **residual = 32 uniformly** on both off-y_par AND on-y_par. Pi_5bilinear is built for the Π²-D-even bilinears and gives residual = 0 on the Π²-D-even strings of even weight (an even number of non-identity letters); the diagonal cell is Π²-D-ODD, and there the operator gives no useful cancellation.

**Path 3 (Q_V × Π composition):** Q_V = Hilbert-space conjugation by V ∈ {X⊗N, Y⊗N, Z⊗N}. Found universal H-flipping V per dephase (V·H·V⁻¹ = −H universally on off-y_par H). Composition with canonical Π or Pi_5bilinear tested: **zero hits.** The required property "Π commutes with [H, ·] for off-y_par H AND gives −L_D − 2σI on dissipator" has no solution among standard Π operators.

Conclusion: none of the three paths found the palindromizer of the off-y_par sub-sector; it is not a better Π but a second mirror after Π, the chiral K of PROOF_F103 §7.4. Tier1Candidate shipped with the empirical anchor + Pure-D Template Rule + subclaim (d) as open work. (2026-05-30: (d) closed modulo M via the chiral route, PROOF_F103 §7.4; 2026-06-10: the converse closed, F111 Tier1Derived.)

Verification scripts: `simulations/f111_path1_operator_search.py`, `simulations/f111_path2_pi5bilinear_test.py`, `simulations/f111_combined_operator_search.py`. Logs: `simulations/results/f111_*.txt`.

## Significance

F111 sharpens F110 Aspect B from "empirical Y-inversion at k = 4 (228:0)" to a structural rule (Pure-D Template Rule) that implies the Y-inversion as immediate corollary AND provides additional structural content (the 36 + 192 + 0 decomposition; condition (ii) of Theorem 2 of the complement connection for subclaim (a)).

The rule predicts the F87-hard set in the diagonal cell EXACTLY at k = N = 4 across all 3 dephase letters with zero exceptions. Subclaim (d) (Mixed + Mixed = soft) is closed modulo M via [PROOF_F103 §7.4](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md), and the hard-direction converse behind subclaim (c) (a lifted diagonal ⟹ hard) closed 2026-06-10 (WindowedConverseAllGammaClaim, no residual), completing F111's Tier1Derived promotion.

The rule's structural origin is the turn by X⊗N on the diagonal: a pure-D template carries an odd number of D letters, so after the rotation to Z dephasing its diagonal flips under the complement, which condition (ii) of Theorem 2 of the complement connection forbids.

## Sibling y_par-axis claims and the cross-axis neighbours

F107 (truly ⟹ y_par=0, Tier1Derived); F109 (y_par-homogeneous mother soft ⟹ y_par=1, Tier1Derived); F110 (HardCellYInversionPattern, Tier1Derived).

The 8 YParity-axis Claims (per `IZ2AxisClaim.Z2Axis == Z2Axis.YParity`) are: F102 (YParityIndependenceAtK3), F103 (F87Z2CubedRefinementN4K3), F105 (F87Z2CubedRefinementN5K3), F106 (F87Z2CubedRefinementN4K4), F107 (TrulyYParityZeroPurity), F109 (MotherSoftYParityOnePurity), F110 (HardCellYInversionPattern), and F111 (HardCellPureDTemplate, this proof). Together they form the YParity-axis classification of the F87 trichotomy.

F108 Part 1+2+3 (Π²-even palindrome family, Tier1Derived 2026-05-25) are **not** YParity-axis sisters: per their `Z2Axis` declarations they live on the BitB axis (Parts 1 and 3) and the BitA axis (Part 2). The diagonal-cell scope of F110 and F111 comes from the colouring of the other three Klein cells ([F110 §2](PROOF_F110_HARD_CELL_Y_INVERSION.md)); F108's mirrors reach the Π²-D-even strings of even weight.

## Open

- **Subclaim (d) closed-form derivation , CLOSED modulo M (2026-05-30), and directly by Theorem 2 (Step 4).** Pair (Mixed, Mixed) at k = N = 4 is F87-soft. The operator-level mechanism is [PROOF_F103 §7.4](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md): in the dephasing basis each Mixed term is a single bit-flip mask, so a Mixed+Mixed pair at full support (k=N) has only two flip generators; two nonzero 𝔽₂ vectors always admit a linear φ with φ(both)=1, and the chiral K = diag((−1)^φ) , a *second* mirror after Π, not the better Π the three Task-1 paths sought , palindromizes the spectrum. So Mixed+Mixed ⟹ bipartite ⟹ soft, modulo the F80 one-sidedness M = −2i(H⊗I) (bit-exact). Subclaim (b) (single Mixed term ⟹ soft) follows the same way (|S| = 1). The converse behind subclaim (c) (a lifted diagonal ⟹ hard) closed 2026-06-10 (WindowedConverseAllGammaClaim, Pascal-Gram positivity F117, no residual), promoting F111 to Tier1Derived.
- **Pure-D Template Rule beyond k = N = 4. Answered** by Theorem 2 of the [complement connection](PROOF_PALINDROME_COMPLEMENT_CONNECTION.md), which decides each pair exactly on a graph of 2^N bitstrings (`DiagonalCellComplementConnectionTests`, Z dephasing). At full support the rule holds at every k (its "if" half below and subclaim (d)'s argument, which needs only two masks), measured at k = 5: 2056 hard pairs, every one with a pure-D template. Its "if" half holds at every N: a pure-D template of the diagonal cell has an odd number of D letters, so after the rotation that takes D to Z (the sign-free letter permutations X → Z → Y → X for X dephasing and X → Y → Z → X for Y, which carry y_par over for X and invert it for Y) its windows give H a diagonal that the complement flips. At the all-zero bitstring every pure-Z window contributes +1, so H_xx there is the number of pure-Z windows and H_x̄x̄ its negative; condition (ii) of Theorem 2 fails, and a mixed template of the cell has no diagonal to compensate. Below full support its "only if" half gives way: hard pairs of two mixed templates appear (192 at k = 4, N = 5; 448 at N = 6, 7, 8), above the window floor N ≥ 2k − 2 exactly [F115](../ANALYTICAL_FORMULAS.md#f115)'s valuation criterion. The counts are read under Z dephasing.
- **Hardware QPU confirmation at k ≥ 3:** open (no F87 QPU confirmations exist beyond Marrakesh k = 2).

∎ (Tier1Derived since 2026-06-10; subclaim (d) closed modulo M via F103 §7.4, 2026-05-30; the (a)/(c) hard-direction converse closed via the windowed all-γ theorem, Pascal-Gram positivity F117)
