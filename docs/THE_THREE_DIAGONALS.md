# THE THREE DIAGONALS: the dephasing diagonal is one of three (Q_X, Q_Y, Q_Z), one basis-S₃ orbit

<!-- Keywords: three dephasing diagonals Q_X Q_Y Q_Z, the one diagonal is one of three,
basis-S3 orbit, single-qubit Clifford basis change Hadamard R_x, Y transpose minus sign
Y^T = -Y, spectrum vs operator gate, build the operator not the spectrum, mirror group
S3 semidirect D4, dephasing diagonal Q = sum_l kron(P_l, P_l transpose), disagreement count,
ThreeDephasingDiagonalsOrbitClaim, inspect root diagonal, R=CPsi2 -->

**Status:** Synthesis. The result is Tier 1 derived, typed as `ThreeDephasingDiagonalsOrbitClaim`, and recomputed live at `inspect --root diagonal` (the "orbit" node). Verified bit-exact at N = 2, 3 (orbit size, the basis moves) and N = 2..4 (same spectrum).
**Date:** 2026-06-15
**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Builds on:**
- [PROOF_ABSORPTION_THEOREM.md](proofs/PROOF_ABSORPTION_THEOREM.md) §4.7: the Z diagonal Q_Z and L_D = γ·(Q_Z − N·I).
- [reflections/ON_THE_ONE_DIAGONAL.md](../reflections/ON_THE_ONE_DIAGONAL.md): one diagonal, read three ways (rate / mirror / judge), the D₄ side.
- [reflections/D_PI_Z_EQUALS_PI_Y.md](../reflections/D_PI_Z_EQUALS_PI_Y.md): D·Π_Z·D = Π_Y, the three-fold at the palindromizer level.
- [PROOF_PI_FACTORS_AS_R_TIMES_D.md](proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md) §5: the S₃ ⋉ D₄ completion, named open.

---

## What this document is about

When an environment keeps "looking at" a quantum system along one direction, asking each spin, over and over, "are you up or down along Z?", that watching (called dephasing) damages one particular ledger of the system. We call that ledger the diagonal. Almost everything this project has proved about how open quantum systems decay was, underneath, the structure of this one ledger.

The point of this document is that the ledger is not special: it is one of three identical siblings, one for each direction the environment could have watched (X, Y, or Z). And they are not three separate things. They are the same object seen from three angles, and a simple change of viewpoint, like turning a die to bring a different face up, carries one into the next. Because they are one shape rotated, they share what matters: the very same set of decay rates.

One trap is worth telling, because the way it was caught is a small lesson in how to check things. The Y sibling hides a minus sign (the Y matrix flips sign when you transpose it; X and Z do not). Forget that sign and the three stop lining up. Here is the trap: a check that only compares their fingerprints, the bare list of decay rates, passes anyway, because flipping the sign of a symmetric list leaves it unchanged. The fingerprint is blind to the sign. So the three siblings had to be built as actual operators and matched one to one, not read off from their spectra. Build the thing; do not trust the shadow.

These three diagonals sit inside a larger web of symmetries (two distinct three-folds, woven together), and one piece of that web is still open. The map at the end of the document says exactly where.

---

## Abstract

Z-dephasing of a spin chain touches exactly one diagonal of the Liouvillian: Q_Z = Σ_l Z_l ⊗ Z_l, whose integer levels are the disagreement count k = popcount(i ⊕ j) and whose rate is Re λ = −2γk (the Absorption Theorem, §4.7). That diagonal is one of **three**: Q_X, Q_Y, Q_Z, one per dephasing axis, and {Q_X, Q_Y, Q_Z} is exactly **one orbit** of the single-qubit Clifford basis-change S₃. The three are conjugate, hence co-spectral, so "the one diagonal" is provably one face of a three-fold. The Y diagonal carries a sign (Yᵀ = −Y), and that sign is what closes the orbit; the lesson of how it was found is the reason this is generated as three operators, not three spectra. The three diagonals (a basis-S₃) and the three readings within one diagonal (the mirror group D₄) are two distinct three-folds, semidirectly coupled: the structure is S₃ ⋉ D₄.

## §1 The object: one diagonal per dephasing axis

Dephasing in light P (for P = Z: the environment repeatedly asking each qubit "are you up or down along Z?", which erases the off-diagonal coherences between states of different Z-content) sends ρ ↦ Σ_l P_l ρ P_l. On the d²-dimensional coherence space (d = 2^N; row-stacking vec, kron(A, B): ρ ↦ A ρ Bᵀ) this is the diagonal generator

  **Q_P = Σ_l kron(P_l, P_lᵀ).**

For Z this is the canonical §4.7 object (of the Absorption Theorem), Q_Z = Σ_l Z_l ⊗ Z_l, diagonal with

  Q_Z(|i⟩⟨j|) = N − 2·popcount(i ⊕ j) = N − 2k,   L_D = γ·(Q_Z − N·I),   Re λ = −2γk,

the disagreement-count rung ladder. Here k = popcount(i ⊕ j) is the **disagreement count**: for a coherence |i⟩⟨j| (row i, column j of the density matrix), it is the number of qubit positions where the bit-strings i and j differ. X, Y, Z are the standard single-qubit Pauli matrices. The X diagonal is built the same way (Xᵀ = X); the Y diagonal is the meeting point of the transpose, because alone among the three Y is antisymmetric, **Yᵀ = −Y**, so

  **Q_Y = Σ_l kron(Y_l, −Y_l) = −Σ_l Y_l ⊗ Y_l.**

The minus sign is not cosmetic. It is exactly what lets the three diagonals form one orbit, and forgetting it is the trap of §3.

## §2 One of three: the basis-S₃ orbit

The three diagonals are one orbit of the single-qubit Clifford basis-change

  **S₃ = ⟨h_zx, h_yz⟩,**   h_zx = Ad_{H^⊗N} (the Hadamard, Z↔X),   h_yz = Ad_{R_x(π/2)^⊗N} (Z↔Y),

each acting by conjugation, Ad_U(Q) = U·Q·U†:

  **h_zx · Q_Z · h_zx⁻¹ = Q_X**   (dev 4.4·10⁻¹⁶ at N = 3),
  **h_yz · Q_Z · h_yz⁻¹ = Q_Y**   (dev 8.9·10⁻¹⁶ at N = 3).

The orbit of Q_Z under ⟨h_zx, h_yz⟩ has size exactly **3** (verified N = 2, 3): it carries Q_Z to the three dephasing diagonals and no further. Conjugate operators have the same eigenvalues (they are co-spectral), so

  **spec(Q_X) = spec(Q_Y) = spec(Q_Z)**   (verified N = 2..4).

The basis-change permutes the three letters {X, Y, Z}; the three diagonals are its orbit. "The one diagonal" of `ON_THE_ONE_DIAGONAL` is one face of this three-fold.

## §3 The lesson: generate the operators, not the spectra

The orbit was found by a physics-first gate, and the way it fired is worth keeping. With the naive Q_Y = +Σ Y_l ⊗ Y_l (forgetting Yᵀ = −Y), the orbit does **not** close: applying h_yz gives −Σ kron(Y, Y), one sign flip from the naive +Σ kron(Y, Y). Yet a **spectrum-only** check passes anyway, because the spectrum of Q_Y is symmetric about 0, so +Q_Y and −Q_Y are co-spectral. A test that compared eigenvalues would have reported success while the operator identity was off by a sign.

So the gate separated **spectrum** from **operator**: the witness builds all three Q_X, Q_Y, Q_Z as operators and checks the conjugation identity h·Q·h⁻¹ = Q′ bit-exactly, rather than merely comparing eigenvalues. The physical content lives in the operator (the Y-transpose); the spectrum is blind to it. That is why the three diagonals are *generated*, not inferred from their spectra.

## §4 The structure: S₃ ⋉ D₄

There are two distinct three-folds in the diagonal story, and they are not the same S₃:

- the **basis-S₃** of this note, which permutes the three **diagonals** {Q_X, Q_Y, Q_Z};
- the **mirror group D₄ = ⟨R, D⟩**, which moves **one** diagonal three ways into its three readings (rate = D fixes Q, mirror = R reflects Q, judge = the F87 truly cell); see `ON_THE_ONE_DIAGONAL` and `PROOF_PI_FACTORS_AS_R_TIMES_D`.

They are semidirectly coupled, each basis move commuting with one mirror generator but not the other:

  [h_zx, D] = 0  but  [h_zx, R] ≠ 0,    [h_yz, R] = 0  but  [h_yz, D] ≠ 0,

so the structure is **S₃ ⋉ D₄** (a semidirect product: the basis-S₃ and the mirror-D₄ lock together, the basis moves acting non-trivially on the mirrors).

One subtlety pins the two apart. **D (the transpose) FIXES every diagonal**, D·Q·D = +Q (it is the rate reading); it does **not** permute them. The basis-S₃ permutes the diagonals; D does not. The proof's phrase "D = the Z↔Y swap" lives on the **palindromizer** Π (a symmetry operator from the referenced proofs that exchanges the slow and fast decay modes, the spectral palindrome), where D·Π_Z·D = Π_Y (`D_PI_Z_EQUALS_PI_Y`), not on the diagonal Q. The Π-level three-fold and the Q-level three-fold are parallel but distinct realizations of the same letter permutation.

## §5 The open completion

The D₄ core ⟨R, D⟩ is typed (`MirrorGroupD4Claim`); the full mirror group of the palindrome family is **S₃-letter-action ⋉ D₄** (order 48), which `PROOF_PI_FACTORS_AS_R_TIMES_D` §5 names open. The Z↔Y transposition is already D (inside D₄); the X↔Z and X↔Y swaps need the X↔Z basis permutation (the dephase-letter swaps Q_zx, Q_yx of `Pi2KleinV4DephaseSwapGroup`) and sit outside ⟨R, D⟩. The basis-change permutes the three dephasing directions {X, Y, Z}; separately, the mirror group can permute the three palindromizers {Π_X, Π_Y, Π_Z}. Both look like S₃ (the symmetric group on three things), but they could be independent structures. The open physics question is whether the dephase-letter-swap S₃ (acting on the palindromizers Π) is the **same** abstract S₃ as the basis-change S₃ of this note (acting on the diagonals Q). Tracked as the open arc `linear_s3_mirror_completion` (`inspect --root arcs`).

## Where it lives

- **Typed:** `ThreeDephasingDiagonalsOrbitClaim` (`compute/RCPsiSquared.Core/Symmetry/`, Tier 1 derived, dual parents `MirrorGroupD4Claim` + `AbsorptionTheoremClaim`, the physics edge that welds the mirror-group and absorption clusters).
- **Live:** the project's C# layer recomputes Q_X, Q_Y, Q_Z and the orbit on demand (`DiagonalWitness`, the "orbit" node; for codebase users, the CLI command `inspect --root diagonal`).
- **Verifier:** `simulations/one_diagonal_mirror_group.py` (self-validating: the Y-transpose gate, the basis-S₃ orbit, the semidirect structure).
