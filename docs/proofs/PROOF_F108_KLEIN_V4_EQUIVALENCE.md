# PROOF F108 Klein-V₄ Equivalence: Parts 2, 3 as Corollaries of Part 1

**Status:** Tier 1 derived universal N. The Klein-V₄ acts simply transitively on the four oriented Π_5b mirrors, and Parts 2 and 3 follow from Part 1 by moves that carry the Lindbladian: the Hadamard or H for Part 2, the quarter turn about X for Part 3. Welle 14.
**Date:** 2026-05-27 (Welle 14)
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:**
- [F108 Part 1](PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md) (the base claim under Z-dephasing)
- [F108 Part 2](PROOF_F108_PART2_PI2X_EVEN_ALWAYS_PALINDROMIC.md) (X-dephasing, BitA twin)
- [F108 Part 3](PROOF_F108_PART3_PI2Y_EVEN_ALWAYS_PALINDROMIC.md) (Y-dephasing, BitB sibling)
- [Klein-V₄ dephase-letter swaps on operator space](PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md) (Welle 12: the Klein-V₄ {I, D, H, Q_zx} of U(4^N))
- [F112 cross-dephase extension](PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md) (Welle 13: the quarter turn about X carries L_Z to L_Y)
- Gate [`simulations/f108_klein_transport_gate.py`](../../simulations/f108_klein_transport_gate.py) (the identities listed in §(e), exact); the verifier `simulations/f108_klein_v4_equivalence_verify.py` with its log `simulations/results/f108_klein_v4_equivalence_verify.txt` (residuals at N = 1, 2, 3)

**Swept before writing:** the three Part proofs and their typed claims (`F108Part1Pi2EvenAlwaysPalindromic` and siblings), the registry entries F108, F114 and F155, [the Klein-V₄ proof](PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md) (its precision note holds the canonical regular action on {Π_X, Π_X⁻¹, Π_Y, Π_Z}), [the F112 cross-dephase proof](PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md) (the quarter turn about X carries L_Z to L_Y), the gate [`f107_f110_route_gate.py`](../../simulations/f107_f110_route_gate.py) G1 (Π_5b(Z) = Π_Z ∘ Ad_{Y^⊗N}, Π_5b(Y) = Π_Y ∘ Ad_{Y^⊗N}, Π_5b(X) = Π_X⁻¹ ∘ Ad_{Y^⊗N}, the V₄ commuting with Ad_{Y^⊗N}, and its swaps of the three variants), `docs/CAUGHT_ERRORS.md` (the entry of 2026-10-04 lists this proof's earlier readings as false; item (vi) of 2026-09-05 the Route-1 necessity in the two Klein proofs), the OpenArcs registry (nothing on F108's Klein action), the typed claim `Pi2KleinV4DephaseSwapGroup` (it cites this proof and still carries the reading this proof replaces), the Diagnostics layer (F108 and Π_5b appear only in passing, in the LindbladBitA/BitB witnesses), F118 and `MirrorGroupD4Claim` (Π_Z = R·D, through which the registry reads Π_5b = Π_Z ∘ Ad_{Y^⊗N}; nothing there on the Klein action), and `cube_old_questions_gate.py` G11 (the quarter turn on the Part-1 strings). The action on the four oriented mirrors is G1's carried over; what this proof adds is reading it as the two orientations of each mirror, and the Lindbladian moves.

## Abstract

F108 came in three flavors, one per dephasing axis. Each says: take a Hamiltonian built from a specific "even" two-site Pauli set, add dephasing along one axis (Z, X, or Y), and an operator Π_5bilinear closes the open-system dynamics into an exact palindrome. The three Parts were proved independently, sharing a skeleton. This proof asks whether Parts 2 and 3 are images of Part 1.

They are, and the reason is one small structure. A mirror can be read in two orientations: if A · L · A⁻¹ = −L − 2σ·I, then A⁻¹ does the same, and so does −A. Per site the three Π_5b maps and one more operator are the four oriented mirrors M_Z, M_X, M_Y = −M_Z⁻¹ and −M_X⁻¹, and the Klein-V₄ of operator-space moves acts on these four simply transitively: D (the transpose, Y ↦ −Y per site) swaps M_Z with M_Y and M_X with −M_X⁻¹, H (the X↔Z swap fixing Y) swaps M_Z with M_X, and Q_zx (the Hadamard's action, X↔Z with Y ↦ −Y) swaps M_X with M_Y. So Π_5b(Y) = (−1)^N · Π_5b(Z)⁻¹ is Π_5b(Z) read the other way, and every move lands on a canonical mirror up to orientation and sign.

What a move does to the Lindbladian decides which Part it reaches. Q_zx is the operator-space form of the Hilbert-space Hadamard and carries L_Z(H₁) to L_X(U·H₁·U†). H = Q_zx·D is not a Hilbert-space unitary but carries L_Z(H₁) to L_X(−U·H₁ᵀ·U†). D carries L_Z(H₁) to L_Z(−H₁ᵀ), the same letter. The step from Z to Y on the Lindbladian is the quarter turn about X, a Clifford outside the V₄ ([F112](PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md)), which fixes Π_5b(Z).

So Part 2 follows from Part 1 by the Hadamard, landing on (−1)^N · Π_5b(X)⁻¹, and by H, landing on Π_5b(X) itself; Part 3 follows from Part 1 by the quarter turn, landing on Π_5b(Z) = (−1)^N · Π_5b(Y)⁻¹. D proves Part 3 too, on the mirror, with the dissipator pillar checked again for Y.

## Introduction

**The motivating question.** F108 Parts 1, 2, 3 gave three independent proofs of the Π_5bilinear palindrome identity, one per dephase letter. After the Klein-V₄ group on operator space ([Klein-V₄ dephase-letter swaps on operator space](PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md)) was used to transport F112 between dephase letters ([F112 cross-dephase extension](PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md)), the natural follow-up was whether the same moves collapse F108 Parts 1/2/3 into one identity plus two corollaries.

**What this proof closes.**

1. **The mirrors.** The four oriented mirrors and the simply transitive Klein-V₄ action on them, each a 4×4 identity per site lifted by tensor power. This is the Π_5b counterpart of the canonical family's Π_Y = Π_Z⁻¹ ([F155](../ANALYTICAL_FORMULAS.md)) and of the regular action on {Π_X, Π_X⁻¹, Π_Y, Π_Z} in the Klein-V₄ proof.
2. **The Lindbladian.** Q_zx and H both carry Z-dephasing to X-dephasing, Q_zx as a Hilbert-space unitary and H with the Hamiltonian sent to minus its transpose; D keeps the letter; the quarter turn about X carries Z to Y.
3. **The corollaries.** Part 2 from Part 1 by the Hadamard or by H; Part 3 from Part 1 by the quarter turn, and by D with the dissipator pillar re-checked.

## (a) Question

F108 Parts 1, 2, 3 (closed 2026-05-25) state the operator-level palindrome identity Π_5b · L · Π_5b⁻¹ = −L − 2σ·I for Π²-even bilinear H + Z, X, Y dephasing respectively. The question: are Parts 2 and 3 images of Part 1 under moves that carry the mirror and the Lindbladian, or do they need proofs of their own?

This proof answers it: both are corollaries of Part 1, Part 2 by the Hadamard or by H, Part 3 by the quarter turn about X or by D with one pillar re-checked.

## (b) Statement

Let Π_5b(d) denote the Π_5bilinear operator for dephase letter d ∈ {Z, X, Y} as built in `Pi5BilinearOperator.BuildFull(N, d)`. Let F108-d be the per-d statement:

  Π_5b(d) · L_d · Π_5b(d)⁻¹ = −L_d − 2σ · I    (operator-level palindrome)

where L_d(H₁) ρ = −i[H₁, ρ] + Σ_l γ_l (d_l ρ d_l − ρ) for any Π²_d-even bilinear Hamiltonian H₁ (real coefficients) and d-dephasing on every site, σ = Σ_l γ_l. The bilinear sets are {XX, YY, YZ, ZY, ZZ} for Parts 1 and 3 and {ZZ, XX, XY, YX, YY} for Part 2.

**Theorem.**
1. **The four oriented mirrors.** The reason is [`f107_f110_route_gate.py`](../../simulations/f107_f110_route_gate.py) G1: Π_5b(Z), Π_5b(Y), Π_5b(X) are Π_Z, Π_Y, Π_X⁻¹ composed with Ad_{Y^⊗N}, and the V₄ commutes with Ad_{Y^⊗N}, so its regular action on the canonical {Π_Z, Π_Y, Π_X, Π_X⁻¹} carries over. Concretely, with M_d the per-site map of Π_5b(d), D = ⊗_l diag(1, 1, 1, −1) on the basis (I, X, Z, Y), H the per-site X↔Z swap fixing I and Y (the Klein-V₄ proof's Q_yx), and Q_zx the per-site X↔Z swap with Y ↦ −Y: M_Y = −M_Z⁻¹, and the Klein-V₄ permutes {M_Z, M_X, M_Y, −M_X⁻¹} simply transitively, D as (M_Z, M_Y)(M_X, −M_X⁻¹), H as (M_Z, M_X)(M_Y, −M_X⁻¹), Q_zx as (M_X, M_Y)(M_Z, −M_X⁻¹). D·H = Q_zx. In particular Π_5b(Y) = (−1)^N · Π_5b(Z)⁻¹, and each element sends the registry variant it does not swap to (−1)^N · Π_5b(X)⁻¹, the operator the registry writes Π_X ∘ Ad_{Y^⊗N}.
2. **Z↔X (Part 1 ↔ Part 2).** Let U := U_H^⊗N and U_op := U ⊗ U^*; on the Pauli basis U_op is Q_zx. Then U_op · L_Z(H₁) · U_op^† = L_X(U H₁ U^†), the Hadamard sends the Part-1 bilinear set onto the Part-2 set, and the transported mirror is (−1)^N · Π_5b(X)⁻¹, the canonical mirror in the other orientation. And H carries L_Z(H₁) to L_X(−U·H₁ᵀ·U†) while carrying Π_5b(Z) to Π_5b(X) itself.
3. **Z↔Y (Part 1 ↔ Part 3).** The quarter turn about X, R = (e^{iπ/4·X})^⊗N, carries L_Z(H₁) to L_Y(R·H₁·R†), maps the Part-1 set onto itself and fixes Π_5b(Z) = (−1)^N · Π_5b(Y)⁻¹. And D carries Π_5b(Z) to Π_5b(Y) and fixes the bilinear set, but keeps every dissipator on its letter, so on that route the dissipator pillar is re-checked for Y (§(c.4)).

All operator identities are exact at universal N via per-site Kronecker factorization, plus the framework's `Pi5BilinearOperator` definition; the Lindbladian moves follow from covariance and the per-site actions. What is gated where is listed in §(e).

## (c) Proof: Part 1 ↔ Part 3

### (c.1) D intertwines Π_5b(Z) and Π_5b(Y) per site

The Welle 12 D operator is per-site D_l = diag(1, 1, 1, −1) on basis (I, X, Z, Y); the −1 sits on the Y entry. On operators D is the transpose, Pᵀ = (−1)^{n_Y}·P for a Pauli string P.

Π_5b(Z) per-site action: I → +X, X → −I, Y → +iZ, Z → −iY.
Π_5b(Y) per-site action: I → +X, X → −I, Y → −iZ, Z → +iY.

The two differ only in the Y/Z 2-cycle phase: +i ↔ −i. On basis (I, X, Z, Y):

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

M_Y differs from M_Z only in the (Z, Y) and (Y, Z) entries. D · M_Z · D negates row Y and column Y: the (Y, Z) entry −i becomes +i (row Y), the (Z, Y) entry +i becomes −i (column Y), and the (I, X) and (X, I) entries are untouched. The result is M_Y. ∎

For §(d), the X variant: Π_5b(X) per-site action I → +Z, Z → −I, X → −iY, Y → +iX.

### (c.2) N-site lift by per-site tensor power

Π_5b(Z) = ⊗_l M_Z, Π_5b(Y) = ⊗_l M_Y and D = ⊗_l D_l are tensor powers, so by the mixed-product property

  D · Π_5b(Z) · D = ⊗_l (D_l · M_Z · D_l) = ⊗_l M_Y = Π_5b(Y),

exact at universal N. ∎

### (c.3) The Π²-even bilinear set is D-invariant

For a 2-site bilinear B = σ_a ⊗ σ_b, D(B) = Bᵀ = (−1)^{n_Y(σ_a) + n_Y(σ_b)} · B. For B ∈ {XX, ZZ, YY, YZ, ZY}:

| B  | n_Y | sign |
|----|-----|------|
| XX | 0   | +1   |
| YY | 2   | +1   |
| YZ | 1   | −1   |
| ZY | 1   | −1   |
| ZZ | 0   | +1   |

So D maps the set {XX, YY, YZ, ZY, ZZ} onto itself, with sign flips on {YZ, ZY}, which appear with both signs in the span.

### (c.4) F108-Y from F108-Z

The Part 1 proof has two pillars:
- (A) Anti-commutation {Π_5b(Z), [B, ·]} = 0 for every B ∈ {XX, YY, YZ, ZY, ZZ}.
- (B) Per-site dissipator identity M_Z · D[Z_l] · M_Z⁻¹ = −D[Z_l] − 2γ_l · I.

(A) transfers by conjugation with D. D is the transpose, so D · [B, ·] · D = [ ·, Bᵀ] = −[Bᵀ, ·] (F114), and with Bᵀ = ±B from (c.3) this is ∓[B, ·]. Conjugating {Π_5b(Z), [B, ·]} = 0 by D gives {Π_5b(Y), ∓[B, ·]} = 0, and the sign does not matter for an anticommutator.

(B) does not transfer by D, because D keeps every single-letter dissipator on its own letter: conjugating the dissipator D[Z_l] by the transpose returns it unchanged. It is re-checked for Y directly. The per-site identity uses only the letter permutation of M (I↔X, Y↔Z) acting on the diagonal D[Z]_pauli = γ · diag(0, −2, −2, 0) on the basis (I, X, Y, Z); each phase of M meets its own inverse in M⁻¹. M_Y has the same permutation, and (I↔X, Y↔Z) exchanges Y-dephasing's undamped pair {I, Y} with its damped pair {X, Z}, as it exchanges {I, Z} with {X, Y} for Z-dephasing. Hence M_Y · D[Y] · M_Y⁻¹ = −D[Y] − 2γ · I.

(A) + (B) give the F108-Y palindrome identity. ∎

### (c.5) The Lindbladian route: the quarter turn about X

D keeps the letter on the Lindbladian: D · L_Z(H₁) · D = L_Z(−H₁ᵀ), a Z-dephasing Lindbladian again. The step to Y-dephasing is made by a Hilbert-space Clifford outside the V₄, R = (e^{iπ/4·X})^⊗N, which fixes X and exchanges Y and Z up to one sign. Lindblad form is unitarily covariant, so R_op · L_Z(H₁) · R_op^† = L_Y(R H₁ R^†) with R_op = R ⊗ R^* ([F112](PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md) uses a quarter turn about X, written R_x(π/2), for F112-Y, where it carries the canonical Π_Z to Π_Y = Π_Z⁻¹; on Π_5b(Z) it acts trivially, see below). On the Part-1 set it acts as XX ↦ XX, YY ↔ ZZ, YZ ↦ −ZY, ZY ↦ −YZ, the set onto itself. And it fixes the mirror: R_op · Π_5b(Z) · R_op^† = Π_5b(Z). Conjugating F108-Z by R_op therefore gives Π_5b(Z) · L_Y(H₃) · Π_5b(Z)⁻¹ = −L_Y(H₃) − 2σ·I for every Part-3 Hamiltonian H₃. Since Π_5b(Z) = (−1)^N · Π_5b(Y)⁻¹, this is Part 3 with its canonical mirror read in the other orientation. ∎

## (d) Proof: Part 1 ↔ Part 2

### (d.1) The Hilbert-space Hadamard rotates Z-deph into X-deph

Let U := U_H^⊗N with U_H = (1/√2)[[1, 1], [1, −1]]. Per site, U_H rotates Pauli operators:

  U_H · X · U_H = Z,    U_H · Z · U_H = X,    U_H · Y · U_H = −Y,    U_H · I · U_H = I.

Let U_op := U ⊗ U^*, the operator-space lift of conjugation by U: vec(U ρ U^†) = U_op · vec(ρ), with row-major vec throughout (in column-major it would read U^* ⊗ U). On the Pauli basis it acts per site as X↔Z, Y ↦ −Y, which is Q_zx. For any Hamiltonian H_1 and any single-site Pauli dephasing:

  U_op · L_Z(H_1) · U_op^† = −i [U H_1 U^†, ·] + Σ_l D[U Z_l U^†] = −i [U H_1 U^†, ·] + Σ_l D[X_l] = L_X(U H_1 U^†).

Lindblad form is unitarily covariant; this is the operator-space lift of ρ → U ρ U^†.

### (d.2) Hadamard rotates Part-1 bilinear set into Part-2 bilinear set

Apply the per-letter Hadamard map (X↔Z, Y→−Y, I→I) to each Part-1 bilinear:

| Part 1 | Hadamard image | sign |
|--------|----------------|------|
| XX     | ZZ             | +1   |
| YY     | YY             | +1   |
| YZ     | YX             | −1   |
| ZY     | XY             | −1   |
| ZZ     | XX             | +1   |

The image set is {XX, XY, YX, YY, ZZ}, Part 2's set. So every Part-1 H maps (with coefficient sign flips on YZ/ZY → YX/XY) to a Part-2 H, and back.

### (d.3) F108-X from F108-Z via the Hadamard

Take any Part-1 H_1 + Z-dephasing on every site. By F108 Part 1:

  Π_5b(Z) · L_Z(H_1) · Π_5b(Z)⁻¹ = −L_Z(H_1) − 2σ · I.    (*)

Conjugating (*) by U_op and using (d.1), with H_2 := U H_1 U^† in the Part-2 class and Π̃ := U_op · Π_5b(Z) · U_op^†:

  Π̃ · L_X(H_2) · Π̃⁻¹ = −L_X(H_2) − 2σ · I.    (**)

So L_X(H_2) admits the palindrome operator Π̃. Per site Q_zx · M_Z · Q_zx = −M_X⁻¹, so Π̃ = (−1)^N · Π_5b(X)⁻¹: the canonical mirror of Part 2 in the other orientation, which palindromizes whenever Π_5b(X) does. It is not ±Π_5b(X) itself (largest entry of the difference 2.0 at N = 1, 2, 3).

### (d.4) F108-X from F108-Z via H, and the Klein action

H = Q_zx · D, the per-site X↔Z swap fixing I and Y (the Klein-V₄ proof's Q_yx). It is not a Hilbert-space unitary: a swap of X and Z that fixes Y is not a Pauli automorphism, since Y = iXZ forces any unitary X↔Z swap to send Y ↦ −Y; on a Hermitian ρ it acts as the Hadamard composed with complex conjugation. It keeps Lindblad form because D does, sending L_Z(H₁) to L_Z(−H₁ᵀ), and Q_zx is conjugation by a Hilbert-space unitary, under which Lindblad form is covariant. For Part 2 it does both jobs:

- **The mirror:** H · Π_5b(Z) · H = Π_5b(X), a 4×4 identity per site lifted by tensor power.
- **The Lindbladian:** H · L_Z(H₁) · H⁻¹ = L_X(−U·H₁ᵀ·U†). D sends L_Z(H₁) to L_Z(−H₁ᵀ), and Q_zx then carries Z-dephasing to X-dephasing; −U·H₁ᵀ·U† runs over the Part-2 class as H₁ runs over the Part-1 class, H₁ᵀ being ±H₁ string by string (§(c.3)).

Conjugating (*) by H gives Part 2 with Π_5b(X) itself.

The three elements and the four oriented mirrors: D swaps M_Z with M_Y and M_X with −M_X⁻¹, H swaps M_Z with M_X and M_Y with −M_X⁻¹, Q_zx swaps M_X with M_Y and M_Z with −M_X⁻¹. Each non-identity element acts as two disjoint swaps with no fixed point, so the group of order four acts freely on the four mirrors and hence simply transitively; D·H = Q_zx. For the canonical Π_d the same group acts simply transitively on {Π_X, Π_X⁻¹, Π_Y, Π_Z}, with D serving Z↔Y, Q_zx serving Z↔X and H = Q_yx serving Y↔X ([the Klein-V₄ proof](PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md)); on Π_5b the roles of H and Q_zx are exchanged.

## (e) Verification

[`simulations/f108_klein_transport_gate.py`](../../simulations/f108_klein_transport_gate.py) checks exactly: the four oriented mirrors, the simply transitive action and D·H = Q_zx, per site (T1); that the label maps of T1 are the same operators as the superoperators of T2 at N = 3 (the transpose, the Hadamard conjugation and their product); the four Lindbladian moves at N = 3 on a Hamiltonian with every two-site string on some bond and two fields, integer couplings, with the unnormalised Cliffords X + Z and I + iX so that every comparison is == (T2); and the mirror statements at N = 3, with wrong mirrors as controls (T3); T2 and T3 use rates 1, 2, 3; the orientation signs (−1)^N are checked at N = 2, 3, 4 on the tensor powers. The quarter turn's per-site letter map (X ↦ X, Y ↦ −Z, Z ↦ Y), from which its action on the Part-1 strings follows, is gated here too (T2); [`cube_old_questions_gate.py`](../../simulations/cube_old_questions_gate.py) G11 checks that it carries the Part-1 set onto itself. The factorization Π_5b(d) = Π_d^{±1} ∘ Ad_{Y^⊗N} and the V₄ commuting with Ad_{Y^⊗N} are [`f107_f110_route_gate.py`](../../simulations/f107_f110_route_gate.py) G1. The Hadamard's bijection of the bilinear sets and the per-site dissipator identities are in the verifier's log below.

The verifier [`simulations/f108_klein_v4_equivalence_verify.py`](../../simulations/f108_klein_v4_equivalence_verify.py) records float residuals at N = 1, 2, 3 (the exact evidence for the identities is the gate above) in [`simulations/results/f108_klein_v4_equivalence_verify.txt`](../../simulations/results/f108_klein_v4_equivalence_verify.txt):

| Claim | N=1 | N=2 | N=3 |
|---|---|---|---|
| F108 Part 1 residual (5 random Hamiltonians per N) | n/a | 0 | 3.8e-16 |
| F108 Part 2 residual (5 random Hamiltonians) | n/a | 0 | 0 |
| F108 Part 3 residual (5 random Hamiltonians) | n/a | 0 | 0 |
| D · Π_5b(Z) · D − Π_5b(Y) | 0 | 0 | 0 |
| Q_zx · Π_5b(Z) · Q_zx − Π_5b(X) | 2.0 | 2.0 | 2.0 |
| H · Π_5b(Y) · H − Π_5b(X) | 2.0 | 2.0 | 2.0 |
| H · Π_5b(Z) · H − Π_5b(X) (`f108_bita_d_search.py`) | 0 | 0 | 0 |
| U_op · L_Z · U_op^† − L_X(rotated H) | 1.6e-16 | 4.5e-15 | 2.2e-14 |
| F108 palindrome of L_X via U_op·Π_5b(Z)·U_op^† | 1.9e-49 | 3.3e-16 | 5.2e-15 |
| F108 palindrome of L_X via canonical Π_5b(X) | 0 | 3.3e-16 | 4.8e-15 |
| Anti-commutation {Π_5b^⊗2, [B, ·]} = 0 for all parts | n/a | 0 (all bilinears) | n/a |
| Dissipator identity M · D[d] · M⁻¹ = −D[d] − 2γ·I | 0 (all parts) | n/a | n/a |
| Bilinear-set bijection Part 1 → Part 2 under per-letter Hadamard | True | n/a | n/a |

The log's Step 6b says the two Lindbladian sides are related only by a non-Lindblad similarity transformation; that is false, since D · L_Z(H₁) · D = L_Z(−H₁ᵀ) is of Lindblad form (gate T2). The table's H row comes from `f108_bita_d_search.py`, not from the log.

## (f) Implications

1. **Both Parts are corollaries of Part 1.** Part 2 by the Hadamard (landing on (−1)^N · Π_5b(X)⁻¹) or by H (landing on Π_5b(X)); Part 3 by the quarter turn about X (landing on Π_5b(Z) = (−1)^N · Π_5b(Y)⁻¹), or by D on the mirror with the dissipator pillar re-checked.

2. **Two jobs, two answers.** On the mirrors the Klein-V₄ acts simply transitively on four oriented operators. On the Lindbladian Q_zx and H move Z to X, D keeps the letter, and Z to Y needs the quarter turn, a Clifford outside the V₄. Whether a move is a Hilbert-space unitary decides neither job by itself: H is not one and carries both.

3. **The three Tier1Derived typed Claims are KEPT separate.** Each Part has its own typed Claim (`F108Part1Pi2EvenAlwaysPalindromic`, `F108Part2Pi2XEvenAlwaysPalindromic`, `F108Part3Pi2YEvenAlwaysPalindromic`) with distinct integration edges (Part 2 fills Part 1's BitA-twin slot; Part 3 sits on the same BitB axis as Part 1 with BitATwinStatus = BitBSpecific). Keeping them as siblings, with cross-references to this proof, preserves the typed-knowledge integration while making the equivalence visible.

4. **The same shape in both families.** For the canonical operators Π_Y = Π_Z⁻¹ (F155) and D · Π_Z · D = Π_Y; for the bilinear family Π_5b(Y) = (−1)^N · Π_5b(Z)⁻¹ and D · Π_5b(Z) · D = Π_5b(Y). In both, the Y mirror is the Z mirror read the other way.

## (g) Open questions and follow-ups

- **The orbit under local unitaries.** Read the question on the letter cube, and extend F108 by one step: call a setting a pair (P, Q) of distinct letters, jumps along P on every site of a chain of N ≥ 2 sites and a Hamiltonian in the span of the two-site strings that commute with Q^⊗N, every real coefficient allowed. Part 1 is (Z, X), Part 2 is (X, Z), Part 3 is (Y, X); the other three pairs are the extension. A product U = u₁ ⊗ … ⊗ u_N, acting on operators as U ⊗ U*, carries a setting to a setting when the Lindbladian of every Hamiltonian of the first becomes the Lindbladian of a Hamiltonian of the second. At H = 0 the dissipator alone must go to one along a single letter, and the Lindblad form of a traceless jump is unique, so every u_l sends Z to ± one common letter P′ and u_l = C_l · e^{iθ_l Z} with C_l a Clifford. On a bond (l, l+1) the rotation turns XX into c_l c_{l+1}·XX − c_l s_{l+1}·XY − s_l c_{l+1}·YX + s_l s_{l+1}·YY, with c = cos 2θ and s = sin 2θ, and each C sends X and Y onto the two letters other than P′. The image must commute with Q′ ⊗ Q′, and of the four strings over those two letters only the two repeated ones do. Whichever way the two sites order the letters, two of the four coefficients must vanish, a cos·sin pair or the cos·cos and sin·sin pair, and either way 2θ_l and 2θ_{l+1} are multiples of π/2: u_l is itself a Clifford. On every bond the image of YZ is u_l(Y) ⊗ (±P′), which commutes with Q′ ⊗ Q′ only if u_l(Y) is ± the third letter R′, and ZY gives the same for u_{l+1}; so on a chain every site carries the same letter permutation (Z ↦ P′, Y ↦ R′, X ↦ Q′), with its Pauli freedom per site, and fixing the setting leaves the local Pauli group. The subgroup that carries settings to settings is therefore the local Clifford group with one letter permutation on all sites, up to phases, of order 6·4^N; Part 1's orbit is all six settings, and in each the transported operator (U ⊗ U*) · Π_5b(Z) · (U ⊗ U*)^† palindromizes the image by covariance. Counted over all 24³ site-dependent local Cliffords at N = 3: six settings, 64 ways each, one permutation per hit, the stabilizer the 4³ Paulis. Read with the three named Parts alone, only the local Paulis map that set to itself, and the maps from Part 1 to Part 2 and Part 3 (the letter swaps X ↔ Z and Y ↔ Z, each with its Pauli freedom) do not form a group. Pieces of the orbit were already in the repo: [F112](PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md)'s quarter turn about X carries Z-dephasing to Y-dephasing, the move from (Z, X) to (Y, X); [F103 §8](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) carries the two flip letters of one dephase letter into each other by the quarter turn about it, the move from (Z, X) to (Z, Y); and the system [Non-Heisenberg Palindrome](../../experiments/NON_HEISENBERG_PALINDROME.md)'s P4 family palindromizes is the setting (Z, Y). New here are the subgroup, its stabilizer, and that nothing continuous joins it. Gate: [`cube_old_questions_gate.py`](../../simulations/cube_old_questions_gate.py) G11 (the Clifford count, the stabilizer, the quarter turn about X carrying Part 1's strings onto Part 3's, and exact F158 end counts for the three settings no Part names; the step from continuous rotations to Cliffords is the argument above).

## (h) Conclusion

F108 Parts 2 and 3 are corollaries of Part 1:
- **The mirrors:** per site M_Y = −M_Z⁻¹, and the Klein-V₄ permutes the four oriented mirrors {M_Z, M_X, M_Y, −M_X⁻¹} simply transitively, D·H = Q_zx.
- **Part 2** by the Hadamard U_op = Q_zx, which carries the Lindbladian and lands on (−1)^N · Π_5b(X)⁻¹, and by H, which carries L_Z(H₁) to L_X(−U·H₁ᵀ·U†) and Π_5b(Z) to Π_5b(X).
- **Part 3** by the quarter turn about X, which carries the Lindbladian and fixes Π_5b(Z) = (−1)^N · Π_5b(Y)⁻¹, and by D on the mirror with the dissipator pillar re-checked.

The three typed Claims remain separate in the registry to preserve their independent integration edges and cross-reference this proof.

∎
