# PROOF F108 Klein-V₄ Equivalence: Parts 2, 3 as Klein-V₄ Corollaries of Part 1

**Status:** Tier 1 derived universal N via two complementary Klein-V₄ routes (D-conjugation for Part 1 ↔ Part 3, Hadamard transport for Part 1 ↔ Part 2). Welle 14.
**Date:** 2026-05-27 (Welle 14)
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md](PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md) (the base claim under Z-dephasing)
- [PROOF_F108_PART2_PI2X_EVEN_ALWAYS_PALINDROMIC.md](PROOF_F108_PART2_PI2X_EVEN_ALWAYS_PALINDROMIC.md) (X-dephasing, BitA twin)
- [PROOF_F108_PART3_PI2Y_EVEN_ALWAYS_PALINDROMIC.md](PROOF_F108_PART3_PI2Y_EVEN_ALWAYS_PALINDROMIC.md) (Y-dephasing, BitB sibling)
- [PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md](PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md) (Welle 12: Klein-V₄ subgroup {I, D, H, Q_zx} of U(4^N))
- [PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md](PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md) (Welle 13: Route 1 / Route 2 pattern for F112)
- Verifier `simulations/_f108_klein_v4_equivalence_verify.py` + log `simulations/results/f108_klein_v4_equivalence_verify.txt` (Welle 14)

## (a) Question

F108 Parts 1, 2, 3 (closed 2026-05-25) state the operator-level palindrome identity Π_5b · L · Π_5b⁻¹ = −L − 2σ·I for Π²-even bilinear H + Z, X, Y dephasing respectively. The PROOF_KLEIN_V4 doc (Welle 12, §Implications point 4) and the PROOF_F112_CROSS_DEPHASE doc (Welle 13, §(g)) both flagged the question:

> Are F108 Parts 2 and 3 Klein-V₄ corollaries of Part 1, or do they require independent proofs?

This proof closes the question. **Both Part 2 and Part 3 are honest Klein-V₄ corollaries of Part 1, but via different mechanisms.** Part 3 follows from Part 1 by operator-space D-conjugation; Part 2 follows by Hilbert-space Hadamard transport. The Klein-V₄ Q_zx swap (Z↔X on the operator-space side, as for canonical Π_d) does NOT swap Π_5b(Z) ↔ Π_5b(X); the operator-space and Hilbert-space mechanisms are genuinely distinct here.

## (b) Statement

Let Π_5b(d) denote the Π_5bilinear operator for dephase letter d ∈ {Z, X, Y} as built in `Pi5BilinearOperator.BuildFull(N, d)`. Let F108-d be the per-d statement:

  Π_5b(d) · L_d · Π_5b(d)⁻¹ = −L_d − 2σ · I    (operator-level palindrome)

where L_d is the Lindbladian for any Π²_d-even bilinear Hamiltonian H and d-dephasing on every site, σ = Σ_l γ_l.

**Theorem (Klein-V₄ equivalence).**
1. **Z↔Y (Part 1 ↔ Part 3) via D:** D · Π_5b(Z) · D = Π_5b(Y) exactly, where D = ⊗_l diag(1, 1, 1, −1) on basis (I, X, Z, Y) is the Welle-12 Klein-V₄ diagonal involution. Together with the fact that the Π²-even bilinear set is fixed under D (Part 1 set = Part 3 set = {XX, YY, YZ, ZY, ZZ}), Part 3 follows from Part 1 by D-equivariance of the proof.
2. **Z↔X (Part 1 ↔ Part 2) via Hadamard:** Let U := U_H^⊗N (Hilbert-space Hadamard ⊗N) and U_op := U ⊗ U^* (operator-space lift, equal to Q_zx of Welle 12 up to a basis-convention permutation). Then U_op · L_Z(H_1) · U_op^† = L_X(U H_1 U^†) for any Part-1 Hamiltonian H_1, and the Hadamard rotation U_H per-letter sends X↔Z, Y↔−Y, I→I, which sends the Part-1 bilinear set bijectively to the Part-2 bilinear set. Hence "L_X admits a palindrome operator in the Π_5b family" follows from Part 1; the canonical Π_5b(X) is one such operator.
3. **Operator-space Q_zx and H do NOT swap Π_5b(Z) ↔ Π_5b(X) ↔ Π_5b(Y):** Q_zx · Π_5b(Z) · Q_zx ≠ ±Π_5b(X), H · Π_5b(Y) · H ≠ ±Π_5b(X) at the operator level (gap = 2.0 in Frobenius distance at all N tested). This negative result shows that Π_5b is a genuinely distinct family from the canonical Π_d palindrome operators; Klein-V₄ acts on Π_d via the full {D, Q_zx, H} but acts on Π_5b(Z) ↔ Π_5b(Y) via D only.

All three claims are Tier 1 derived universal N via per-site Kronecker factorization (each identity reduces to a 4×4 single-site check) plus the framework's `Pi5BilinearOperator` definition.

## (c) Proof — Part 1 ↔ Part 3 via D

### (c.1) D intertwines Π_5b(Z) and Π_5b(Y) per site

The Welle 12 D operator is per-site D_l = diag(1, 1, 1, −1) on basis (I, X, Z, Y); the −1 sits on the Y entry (letter index 3 in the a + 2·b packing, with Y = (1, 1)).

Π_5b(Z) per-site action: I → +X, X → −I, Y → +iZ, Z → −iY.
Π_5b(Y) per-site action: I → +X, X → −I, Y → −iZ, Z → +iY.

The two differ only in the Y/Z 2-cycle phase: +i ↔ −i. Concretely, write the per-site 4×4 matrices on basis (I, X, Z, Y):

```
              col=I  col=X  col=Z  col=Y
M_Z = π_Z:   [  0    -1     0      0   ]   (row=I: X→−I)
              [  1     0     0      0   ]   (row=X: I→+X)
              [  0     0     0     +i   ]   (row=Z: Y→+iZ)
              [  0     0    -i     0   ]   (row=Y: Z→−iY)

M_Y = π_Y:   [  0    -1     0      0   ]   (row=I: X→−I)
              [  1     0     0      0   ]   (row=X: I→+X)
              [  0     0     0     -i   ]   (row=Z: Y→−iZ)
              [  0     0    +i     0   ]   (row=Y: Z→+iY)
```

Note M_Y differs from M_Z only in the bottom-right 2×2 block: the (Z, Y) and (Y, Z) entries flip sign.

**Per-site D-conjugation:** D = diag(1, 1, 1, −1). Left-multiplication by D negates row Y (the last row); right-multiplication by D^{−1} = D negates column Y (the last column). For M_Z:

- Row Y (entries (0, 0, −i, 0)) → (0, 0, +i, 0) after left-D (sign flip on row 3).
- Then column Y (last column of intermediate) gets sign-flipped on the entries (0, 0, +i, 0): no change since column Y is (0, 0, +i, 0) and the row-Y entry is in row Y; after right-D, the column Y entries (rows I, X, Z, Y) become (0, 0, +i, 0) sign-flipped at row Y → (0, 0, +i, 0). But wait the (Z, Y) entry in original M_Z is +i in column Y row Z: right-D negates column Y, so (Z, Y) entry flips to −i.

Let me be explicit. Original M_Z entries (row, col):
- (I, X) = −1
- (X, I) = +1
- (Z, Y) = +i
- (Y, Z) = −i
- All others zero.

After D · M_Z (left-multiply, negate row Y): the entry (Y, Z) flips sign from −i to +i. All other entries unchanged.

After D · M_Z · D (right-multiply by D, negate column Y): the entry (Z, Y) flips sign from +i to −i. Entry (Y, Z) is in column Z, unchanged.

Resulting matrix:
- (I, X) = −1   (unchanged)
- (X, I) = +1   (unchanged)
- (Z, Y) = +i → −i   (sign flipped by right-D)
- (Y, Z) = −i → +i   (sign flipped by left-D)
- All others zero.

This is precisely M_Y. ∎

### (c.2) N-site lift by per-site tensor power

Both Π_5b(Z) = ⊗_l M_Z and Π_5b(Y) = ⊗_l M_Y are pure tensor powers (by definition; `BuildFullUncached` does per-site action then flat-index reassembly). D = ⊗_l D_l is also a pure tensor power.

By the Kronecker mixed-product property:

  D · Π_5b(Z) · D = (⊗_l D_l) · (⊗_l M_Z) · (⊗_l D_l) = ⊗_l (D_l · M_Z · D_l) = ⊗_l M_Y = Π_5b(Y).

Bit-exact universal N. ∎

### (c.3) The Π²-even bilinear set is D-invariant

The bilinear sets are {XX, YY, YZ, ZY, ZZ} for both Part 1 and Part 3 (per the F108 Part 3 doc, "Π²_Y-even bilinear set is identical to Π²_Z-even"). D = ⊗_l diag(1, 1, 1, −1) acts diagonally on Pauli strings with sign (−1)^n_Y, so for a 2-site bilinear B = σ_a ⊗ σ_b the conjugation D · B · D = (−1)^{n_Y(σ_a) + n_Y(σ_b)} · B. For B ∈ {XX, ZZ, YY, YZ, ZY}:

| B  | n_Y | sign |
|----|-----|------|
| XX | 0   | +1   |
| YY | 2   | +1   |
| YZ | 1   | −1   |
| ZY | 1   | −1   |
| ZZ | 0   | +1   |

So D maps the set {XX, YY, YZ, ZY, ZZ} bijectively onto itself, with sign flips on {YZ, ZY} (which appear with both signs in the sum). The set is D-invariant.

### (c.4) F108-Y from F108-Z via D-equivariance

The Part 1 proof has two pillars:
- (A) Anti-commutation {Π_5b(Z)^⊗N, [B, ·]} = 0 for every B ∈ {XX, YY, YZ, ZY, ZZ}.
- (B) Per-site dissipator identity M_Z · D[Z_l] · M_Z⁻¹ = −D[Z_l] − 2γ_l · I.

Both pillars transfer to Part 3 via:
- (A) The bilinear set is D-invariant (c.3). The anti-commutation is a 4^N-dim matrix identity {Q, C} = 0 where Q = Π_5b(Z) and C = [B, ·]. Conjugating both sides by D gives {D·Q·D, C'} = 0 where C' = D · [B, ·] · D⁻¹. Both Π_5b(Z) ↦ D·Q·D = Π_5b(Y) (by c.2) and [B, ·] ↦ D · [B, ·] · D⁻¹ = [D · B · D, ·] for any B (the commutator superoperator transforms covariantly under inner automorphisms). Combined with D · B · D = ±B per c.3, the anti-commutation identity transfers to Part 3.
- (B) The per-site M_Z · D[Z] · M_Z⁻¹ identity uses the (I↔X, Y↔Z) per-site permutation in M and the diagonal-permutation argument on D[Z]_pauli = γ · diag(0, −2, −2, 0). For M_Y, the per-site permutation is the SAME (I↔X, Y↔Z); only the Y/Z 2-cycle phase differs. In the diagonal-permutation argument, phase factors cancel pairwise on each 2-cycle (Y/Z phases are (+i, −i) for M_Z and (−i, +i) for M_Y; the products +i · −i = +1 and −i · +i = +1 are identical). Hence the dissipator identity transfers from M_Z · D[Z] · M_Z⁻¹ to M_Y · D[Y] · M_Y⁻¹.

Combining (A) + (B) gives the F108-Y palindrome identity. ∎

## (d) Proof — Part 1 ↔ Part 2 via Hadamard transport

### (d.1) The Hilbert-space Hadamard rotates Z-deph into X-deph

Let U := U_H^⊗N with U_H = (1/√2)[[1, 1], [1, −1]]. Per-site, U_H rotates Pauli operators:

  U_H · X · U_H = Z,    U_H · Z · U_H = X,    U_H · Y · U_H = −Y,    U_H · I · U_H = I.

Let U_op := U ⊗ U^*, the operator-space lift of conjugation by U: vec(U ρ U^†) = U_op · vec(ρ). Then for any Hamiltonian H_1 and any dephasing operator c_l (single-site Pauli on site l):

  U_op · L_Z(H_1) · U_op^† = U_op · [−i [H_1, ·] + Σ_l D[Z_l]] · U_op^†
                            = −i [U H_1 U^†, ·] + Σ_l D[U Z_l U^†]
                            = −i [U H_1 U^†, ·] + Σ_l D[X_l]
                            = L_X(U H_1 U^†).

Bit-exact (Lindblad-form is unitarily covariant; this is the operator-space lift of ρ → U ρ U^†). Verified numerically: ‖U_op · L_Z · U_op^† − L_X(rotated H)‖ < 1e-13 at N = 1, 2, 3.

### (d.2) Hadamard rotates Part-1 bilinear set into Part-2 bilinear set

Apply the per-letter Hadamard map (X↔Z, Y→−Y, I→I) to each Part-1 bilinear:

| Part 1 | Hadamard image | sign |
|--------|----------------|------|
| XX     | ZZ             | +1   |
| YY     | YY             | +1   |
| YZ     | YX             | −1   |
| ZY     | XY             | −1   |
| ZZ     | XX             | +1   |

The image set is {XX, XY, YX, YY, ZZ} — bit-exact equal to Part 2's set {ZZ, XX, XY, YX, YY}. So every Part-1 H is bijectively mapped (with coefficient sign flips on YZ/ZY → YX/XY) to a Part-2 H, and vice versa.

### (d.3) F108-X from F108-Z via Hadamard

Take any Part-1 H_1 + Z-dephasing on every site. By F108 Part 1:

  Π_5b(Z) · L_Z(H_1) · Π_5b(Z)⁻¹ = −L_Z(H_1) − 2σ · I.    (*)

Conjugate (*) by U_op:

  U_op Π_5b(Z) U_op^† · U_op L_Z(H_1) U_op^† · U_op Π_5b(Z)⁻¹ U_op^†
    = −U_op L_Z(H_1) U_op^† − 2σ · I.

Using (d.1): U_op · L_Z(H_1) · U_op^† = L_X(H_2) where H_2 := U H_1 U^† ∈ Part-2 class. Let Π̃ := U_op · Π_5b(Z) · U_op^†. The above becomes

  Π̃ · L_X(H_2) · Π̃⁻¹ = −L_X(H_2) − 2σ · I.    (**)

So L_X(H_2) admits the palindrome operator Π̃. Π̃ is a unitary signed permutation in the Π_5bilinear family (a U_op-rotated relative of Π_5b(Z)). It is NOT equal to the canonical Π_5b(X): we verified ‖Π̃ − Π_5b(X)‖_max = 2.0 (bit-exact non-zero) at N = 1, 2, 3. But both Π̃ AND Π_5b(X) satisfy the F108 palindrome identity for L_X(H_2): the verifier directly checks F108(L_X, Π̃) = 0 and F108(L_X, Π_5b(X)) = 0 simultaneously at N = 1, 2, 3 (both residuals < 1e-14).

Hence Part 2's claim (the existence of Π_5b(X) as a palindrome operator) follows from Part 1 in two ways:
- **Existence-version:** "L_X(H_2) admits a Π_5b-family palindrome operator" follows by Hadamard transport (using Π̃).
- **Specific-Π̃ version:** "Π_5b(X) is a palindrome operator for L_X(H_2)" requires the direct per-axis check (Route 1; the Part 2 proof). Π_5b(X) is a CANONICAL CHOICE; Π̃ is an equally valid representative of the Π_5b family for L_X.

The substantive content (existence of an exact-palindrome operator in the Π_5b family for L_X) is a Klein-V₄ corollary of Part 1 via Hadamard. ∎

### (d.4) The negative result: Q_zx does NOT swap Π_5b(Z) ↔ Π_5b(X) at the operator level

Welle 12 proved Q_zx · Π_Z · Q_zx = Π_X for the canonical Π_d operators. The natural question is whether Q_zx · Π_5b(Z) · Q_zx = Π_5b(X) also holds.

**No.** Numerical check at N = 1, 2, 3: ‖Q_zx · Π_5b(Z) · Q_zx − Π_5b(X)‖_max = 2.0 (and ‖Q_zx · Π_5b(Z) · Q_zx + Π_5b(X)‖_max = 2.0 also). The same holds for H · Π_5b(Y) · H vs Π_5b(X). So the Klein-V₄ subgroup action on Π_d does NOT extend straightforwardly to the Π_5bilinear family.

Why does D work but Q_zx fail? D is the IDENTITY on (I, X, Z) and only negates Y. Π_5b(Z) and Π_5b(Y) differ ONLY in the Y/Z 2-cycle phase; the I↔X part of the permutation is identical. D's localized action on Y is enough to flip the Y/Z phase. Q_zx instead permutes the BASIS INDICES X↔Z; this restructures the I↔X and Y↔Z 2-cycles of Π_5b(Z) into the I↔Z and X↔Y 2-cycles of Π_5b(X). The result is a Π_5b-family operator (the U_op-transport Π̃ from d.3) but NOT the canonical Π_5b(X) chosen for Part 2.

This is the precise honesty distinction: Klein-V₄ acts on canonical Π_d via the full subgroup {D, Q_zx, H}; on Π_5bilinear it acts via D only (Z↔Y), and the X-deph case enters via the Hilbert-space Hadamard at the level of the Lindbladian.

## (e) Verification

The numerical verifier [`simulations/_f108_klein_v4_equivalence_verify.py`](../../simulations/_f108_klein_v4_equivalence_verify.py) confirms all claims at N = 1, 2, 3 with output saved to [`simulations/results/f108_klein_v4_equivalence_verify.txt`](../../simulations/results/f108_klein_v4_equivalence_verify.txt):

| Claim | N=1 | N=2 | N=3 |
|---|---|---|---|
| F108 Part 1 residual (5 random Hamiltonians per N) | n/a | 0 | 3.8e-16 |
| F108 Part 2 residual (5 random Hamiltonians) | n/a | 0 | 0 |
| F108 Part 3 residual (5 random Hamiltonians) | n/a | 0 | 0 |
| D · Π_5b(Z) · D − Π_5b(Y) | 0 | 0 | 0 |
| Q_zx · Π_5b(Z) · Q_zx − Π_5b(X) (NEGATIVE) | 2.0 | 2.0 | 2.0 |
| H · Π_5b(Y) · H − Π_5b(X) (NEGATIVE) | 2.0 | 2.0 | 2.0 |
| U_op · L_Z · U_op^† − L_X(rotated H) | 1.6e-16 | 4.5e-15 | 2.2e-14 |
| F108 palindrome of L_X via U_op·Π_5b(Z)·U_op^† | 1.9e-49 | 3.3e-16 | 5.2e-15 |
| F108 palindrome of L_X via canonical Π_5b(X) | 0 | 3.3e-16 | 4.8e-15 |
| Anti-commutation {Π_5b^⊗2, [B, ·]} = 0 for all parts | n/a | 0 (all bilinears) | n/a |
| Dissipator identity M · D[d] · M⁻¹ = −D[d] − 2γ·I | 0 (all parts) | n/a | n/a |
| Bilinear-set bijection Part 1 → Part 2 under per-letter Hadamard | True | n/a | n/a |
| Bilinear-set fixity Part 1 = Part 3 under D | True | n/a | n/a |

## (f) Implications

1. **Consolidation IS possible, with full honesty.** F108 Parts 2 and 3 are Klein-V₄ corollaries of Part 1 via complementary mechanisms:
   - Part 1 → Part 3 by operator-space D-conjugation (bit_b-axis-preserving).
   - Part 1 → Part 2 by Hilbert-space Hadamard transport (bit_a-axis ↔ bit_b-axis swap).
   The Welle-13 Route 1 / Route 2 dichotomy applies cleanly here: Route 1 (per-axis re-run made explicit by D) for Part 3; Route 2 (Hadamard transport) for Part 2.

2. **Operator-level Klein-V₄ on Π_5b is PARTIAL, not full.** Only the {I, D} subgroup acts on Π_5bilinear as expected (D swaps Z↔Y phase variant). The Q_zx and H involutions do NOT swap Π_5b(Z) ↔ Π_5b(X) or Π_5b(Y) ↔ Π_5b(X) at the operator level. The X-deph case needs the Hilbert-space Hadamard mechanism (Route 2). This refines the Welle-12 PROOF_KLEIN_V4 doc's §Implications point 4 (the conjecture had asked all three to be Klein-V₄ equivalent without distinguishing the mechanism).

3. **The three Tier1Derived typed Claims are KEPT separate.** Each Part has its own typed Claim wired into the registry (`F108Part1Pi2EvenAlwaysPalindromic`, `F108Part2Pi2XEvenAlwaysPalindromic`, `F108Part3Pi2YEvenAlwaysPalindromic`); they hold distinct integration edges (Part 1 is the BitB twin of Part 2; Part 3 sits on the same BitB axis as Part 1 with BitATwinStatus = BitBSpecific). Keeping them as siblings, with cross-references to this proof, preserves the typed-knowledge integration while making the equivalence visible. The PROOF_KLEIN_V4 doc's "F108 cross-dephasing (conjecture, not yet proven)" caveat (§Implications point 4) is now closed positively.

4. **For the Klein-V₄ on canonical Π_d:** D · Π_Z · D = Π_Y exactly (Welle 12 Task 1). For Π_5bilinear: D · Π_5b(Z) · D = Π_5b(Y) ALSO exactly. So D acts the same way on both Π families on the Z↔Y axis. This is consistent with the Y/Z phase-flip structure shared by both Π_d and Π_5b on the bit_b axis.

5. **The Π_5b family is not unique.** For a Part-2 L_X, both the canonical Π_5b(X) and the U_op-transport U_op · Π_5b(Z) · U_op^† achieve the F108 palindrome. There may be infinitely many other Π_5b-family representatives. The Welle-14 verifier explicitly demonstrates this non-uniqueness at N = 1, 2, 3.

## (g) Open questions and follow-ups

- **Characterization of the full Π_5bilinear orbit under U_op for arbitrary unitary U** (not just Hadamard): which subgroup of U(2)^⊗N acting via U_op on the operator space transports valid F108 palindrome operators? The Hadamard case shows the unit-cell subgroup; whether a larger subgroup works is open.

- **Part 2's BitA-twin status:** F108 Part 2 sits on BitA axis. Could the bit_a-axis version of D (a Z↔X swap analog) intertwine Π_5b(Z) ↔ Π_5b(X) directly on operator space? Welle 14 shows the canonical Q_zx fails. A search for a different operator-space involution that intertwines Π_5b(Z) ↔ Π_5b(X) directly (analogous to D for Z↔Y) is open. The natural ansatz would be an operator-space involution sitting on the bit_a axis with the Y phase fixed but the X/Z phases flipped; we have not formally enumerated.

- **F109's status after Welle 14:** F109's Step 5 closure across {Z, X, Y} is now DOUBLY anchored: each Part separately (independent proofs, 2026-05-25), and now via Klein-V₄ corollary from Part 1 (this proof). F109's full unconditional Tier1Derived status is reinforced; the cross-dephase coverage is now structurally redundant rather than three independent results.

## (h) Conclusion

F108 Parts 2 and 3 are Klein-V₄ corollaries of Part 1 via:
- **Part 3 ↔ Part 1 via D** (operator-space Z↔Y intertwiner; bit_b-axis-preserving). The proof is D-equivariant.
- **Part 2 ↔ Part 1 via Hadamard U_op = U_H^⊗N ⊗ (U_H^⊗N)^*** (Hilbert-space lift; bit_a ↔ bit_b swap on the bilinear set). Existence of an F108 palindrome operator for L_X follows; the canonical Π_5b(X) is a specific representative.
- **Operator-space Q_zx, H do NOT swap Π_5b(Z) ↔ Π_5b(X)** (negative result; the action of Klein-V₄ on Π_5b is partial). The Welle-13 Route 1 / Route 2 distinction matters here as it did for F112.

The three typed Claims remain separate in the registry to preserve their independent integration edges, but cross-reference this Klein-V₄ equivalence proof. The PROOF_KLEIN_V4 doc's conjecture about F108 cross-dephasing equivalence is now closed positively, with the precise mechanism distinction made explicit.

∎
