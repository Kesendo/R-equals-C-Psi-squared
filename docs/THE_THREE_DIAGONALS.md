# THE THREE DIAGONALS: the dephasing diagonal is one of three (Q_X, Q_Y, Q_Z), one basis-S₃ orbit

<!-- Keywords: three dephasing diagonals Q_X Q_Y Q_Z, the one diagonal is one of three,
basis-S3 orbit, single-qubit Clifford basis change Hadamard R_x, Y transpose minus sign
Y^T = -Y, spectrum vs operator gate, build the operator not the spectrum, mirror group,
letter S3 does not normalize D4, no order-48 completion, closure 96 times 2^N,
dephasing diagonal Q = sum_l kron(P_l, P_l transpose), disagreement count,
ThreeDephasingDiagonalsOrbitClaim, inspect root diagonal, R=CPsi2 -->

**Status:** Synthesis. The result is Tier 1 derived, typed as `ThreeDephasingDiagonalsOrbitClaim`, and recomputed live at `inspect --root diagonal` (the "orbit" node). Verified at N = 2, 3 (orbit size, the basis moves) and N = 2..4 (same spectrum). These are float comparisons at tolerance 10⁻¹⁰, not bit-exact ones: both conjugators carry a 1/√2.
**Date:** 2026-06-15
**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Builds on:**
- [PROOF_ABSORPTION_THEOREM.md](proofs/PROOF_ABSORPTION_THEOREM.md) §4.7: the Z diagonal Q_Z and L_D = γ·(Q_Z − N·I).
- [reflections/ON_THE_ONE_DIAGONAL.md](../reflections/ON_THE_ONE_DIAGONAL.md): one diagonal, read three ways (rate / mirror / judge). The three readings are there in words; the D₄ that carries them is in the R·D proof below, not in that reflection.
- [reflections/D_PI_Z_EQUALS_PI_Y.md](../reflections/D_PI_Z_EQUALS_PI_Y.md): D·Π_Z·D = Π_Y, the three-fold at the palindromizer level.
- [PROOF_PI_FACTORS_AS_R_TIMES_D.md](proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md) §5: the letter-S₃ side of the mirror group, and its resolution (the expected order-48 completion does not exist).

---

## What this document is about

When light keeps arriving at every spin of a quantum system, separating each place that holds one way of being on one face and the other on the other, and passing through wherever the two faces agree, that separating (called dephasing) damages one particular ledger of the system. We call that ledger the diagonal. Almost everything this project has proved about how open quantum systems decay was, underneath, the structure of this one ledger.

The point of this document is that the ledger is not special: it is one of three identical siblings, one for each axis a system can be separated along (X, Y, or Z). And they are not three separate things. They are the same object seen from three angles, and a simple change of viewpoint, like turning a die to bring a different face up, carries one into the next. Because they are one shape rotated, they share what matters: the very same set of decay rates.

One trap is worth telling, because the way it was caught is a small lesson in how to check things. The Y sibling hides a minus sign (the Y matrix flips sign when you transpose it; X and Z do not). Forget that sign and the three stop lining up. Here is the trap: a check that only compares their fingerprints, the bare list of decay rates, passes anyway, because flipping the sign of a symmetric list leaves it unchanged. The fingerprint is blind to the sign. So the three siblings had to be built as actual operators and matched one to one, not read off from their spectra. Build the thing; do not trust the shadow.

These three diagonals sit inside a larger web of symmetries, two distinct three-folds. The expectation was that the two would lock into a single larger group; they do not, and the map at the end of the document says exactly where that expectation broke.

---

## Abstract

Z-dephasing of a spin chain touches exactly one diagonal of the Liouvillian: Q_Z = Σ_l Z_l ⊗ Z_l, whose integer levels are the disagreement count k = popcount(i ⊕ j), making L_D diagonal with entries −2γk (the Absorption Theorem, §4.7). That diagonal is one of **three**: Q_X, Q_Y, Q_Z, one per dephasing axis, and {Q_X, Q_Y, Q_Z} is exactly **one orbit**, of size three, under the single-qubit letter moves. The three are conjugate, hence co-spectral, so "the one diagonal" is provably one face of a three-fold. The Y diagonal carries a sign (Yᵀ = −Y), and that sign is what closes the orbit; the lesson of how it was found is the reason this is generated as three operators, not three spectra. The three diagonals (permuted by the letter-S₃) and the three readings within one diagonal (the mirror group D₄) are two distinct three-folds, and they do **not** assemble into a semidirect product: the letter moves do not normalise D₄ (§4).

## §1 The object: one diagonal per dephasing axis

Dephasing in letter P (for P = Z: light arriving at each qubit again and again, separating along Z, erasing the off-diagonal coherences between states that differ anywhere along Z) sends ρ ↦ Σ_l P_l ρ P_l. On the d²-dimensional coherence space (d = 2^N; row-stacking vec, kron(A, B): ρ ↦ A ρ Bᵀ) this is the generator

  **Q_P = Σ_l kron(P_l, P_lᵀ).**

All three are diagonal in the Pauli-string basis, which is what "the diagonal" names throughout; in the computational coherence basis just written down, only Q_Z is (Q_X and Q_Y carry off-diagonal entries of modulus 1 there). For Z the two bases agree, and this is the canonical §4.7 object (of the Absorption Theorem), Q_Z = Σ_l Z_l ⊗ Z_l, diagonal with

  Q_Z(|i⟩⟨j|) = N − 2·popcount(i ⊕ j) = N − 2k,   L_D = γ·(Q_Z − N·I),   L_D(|i⟩⟨j|) = −2γk·|i⟩⟨j|,

the disagreement-count rung ladder. Those are the dissipator's own levels. They are the decay rates of L = L_H + L_D only where the Hamiltonian does not mix weight sectors; in general the Absorption Theorem gives Re λ = −2γ⟨n_XY⟩_v, an eigenvector-weighted average that need not be an integer multiple of 2γ. Here k = popcount(i ⊕ j) is the **disagreement count**: for a coherence |i⟩⟨j| (row i, column j of the density matrix), it is the number of qubit positions where the bit-strings i and j differ. X, Y, Z are the standard single-qubit Pauli matrices. The X diagonal is built the same way (Xᵀ = X); the Y diagonal is the meeting point of the transpose, because alone among the three Y is antisymmetric, **Yᵀ = −Y**, so

  **Q_Y = Σ_l kron(Y_l, −Y_l) = −Σ_l Y_l ⊗ Y_l.**

The minus sign is not cosmetic. It is exactly what lets the three diagonals form one orbit, and forgetting it is the trap of §3.

## §2 One of three: the letter orbit

The three diagonals are one orbit of the single-qubit basis moves

  h_zx = Ad_{H^⊗N} (the Hadamard, Z↔X),   h_yz = Ad_{R_x(π/2)^⊗N} (Z↔Y),

each acting by conjugation, Ad_U(Q) = U·Q·U†:

  **h_zx · Q_Z · h_zx⁻¹ = Q_X**,
  **h_yz · Q_Z · h_yz⁻¹ = Q_Y**   (both machine-zero, below 10⁻¹⁴ at N = 3; the exact residual is BLAS-dependent, so no digits are pinned here).

The group these two generate is **not** S₃ and it is worth saying why, because the name has been wrong here. R_x(π/2) is a quarter-turn on the letters (Y → Z → −Y → −Z), so h_yz has order 4 and ⟨h_zx, h_yz⟩ has order **24**, the single-qubit Clifford group modulo phases. The genuine letter-S₃ is generated by the **involutive** transposition Cliffords, ⟨h_zx, t_yz⟩ with t_yz = Ad of (Y + Z)/√2, and has order 6. Both groups induce the same action here: the orbit of Q_Z has size exactly **3** under either (verified N = 2, 3, 4), carrying Q_Z to the three dephasing diagonals and no further, and both realise the full S₃ on those three. The order-6 group acts faithfully; the order-24 one acts through its quotient by a kernel of order 4, because Q_P is quadratic in P_l and so cannot see a letter move that only changes signs. Conjugate operators have the same eigenvalues (they are co-spectral), so

  **spec(Q_X) = spec(Q_Y) = spec(Q_Z)**   (verified N = 2..4).

The letter moves permute {X, Y, Z}; the three diagonals are their orbit. "The one diagonal" of `ON_THE_ONE_DIAGONAL` is one face of this three-fold.

## §3 The lesson: generate the operators, not the spectra

The orbit was found by a physics-first gate, and the way it fired is worth keeping. With the naive Q_Y = +Σ Y_l ⊗ Y_l (forgetting Yᵀ = −Y), the orbit does **not** close: applying h_yz gives −Σ kron(Y, Y), one sign flip from the naive +Σ kron(Y, Y). Yet a **spectrum-only** check passes anyway, because the spectrum of Q_Y is symmetric about 0, so +Q_Y and −Q_Y are co-spectral. A test that compared eigenvalues would have reported success while the operator identity was off by a sign.

So the gate separated **spectrum** from **operator**: the witness builds all three Q_X, Q_Y, Q_Z as operators and checks the conjugation identity h·Q·h⁻¹ = Q′ entry by entry (to a 10⁻¹⁰ tolerance; the 1/√2 in the conjugators rules out an exact comparison), rather than merely comparing eigenvalues. The physical content lives in the operator (the Y-transpose); the spectrum is blind to it. That is why the three diagonals are *generated*, not inferred from their spectra.

## §4 The structure: two three-folds that do not lock into one group

There are two distinct three-folds in the diagonal story, and they are not the same S₃:

- the **letter-S₃** of this note, which permutes the three **diagonals** {Q_X, Q_Y, Q_Z};
- the **mirror group D₄ = ⟨R, D⟩**, which moves **one** diagonal three ways into its three readings (rate = D fixes Q, mirror = R reflects Q, judge = the F87 truly cell); see `ON_THE_ONE_DIAGONAL` and `PROOF_PI_FACTORS_AS_R_TIMES_D`.

Each letter move commutes with one mirror generator but not the other:

  [h_zx, D] = 0  but  [h_zx, R] ≠ 0,    [h_yz, R] = 0  but  [h_yz, D] ≠ 0.

That pattern is a fact about the two factors, and it is **not** evidence of a semidirect product; relations of exactly this shape hold in any group containing both. The semidirect product S₃ ⋉ D₄ would need the letter moves to normalise D₄, and they do not:

  **h_zx · R · h_zx⁻¹ = the one-sided multiplication by Z^⊗N, which is outside the eight elements of ⟨R, D⟩.**

The coherence-space closure ⟨R, D, h_zx, t_yz⟩ has order 96·2^N (384 at N = 2, 768 at N = 3), not the 48 an S₃ ⋉ D₄ would have. This is settled in `PROOF_PI_FACTORS_AS_R_TIMES_D` §5 ("Resolution of the S₃ side", 2026-06-15, with its own gate at N = 2, 3): the abstract S₃ ⋉ D₄ has no faithful finite realisation here. What survives, and is the content of this note, is the orbit: three diagonals, one letter-S₃ action, one shared spectrum.

One subtlety pins the two apart. **D (the transpose) FIXES every diagonal**, D·Q·D = +Q (it is the rate reading); it does **not** permute them. The basis-S₃ permutes the diagonals; D does not. The proof's phrase "D = the Z↔Y swap" lives on the **palindromizer** Π (a symmetry operator from the referenced proofs that exchanges the slow and fast decay modes, the spectral palindrome), where D·Π_Z·D = Π_Y (`D_PI_Z_EQUALS_PI_Y`), not on the diagonal Q. The Π-level three-fold and the Q-level three-fold are parallel but distinct realizations of the same letter permutation.

## §5 The completion, and how it closed

The D₄ core ⟨R, D⟩ is typed (`MirrorGroupD4Claim`). The completion by the letter moves was expected to be the order-48 **S₃-letter-action ⋉ D₄**, and `PROOF_PI_FACTORS_AS_R_TIMES_D` §5 named it open on the day this note was written. It closed the same day, and the answer is the one §4 states: the linear letter-S₃ exists as superoperators, but it does not normalise D₄, so the order-48 group is not there to be built.

The question this section used to pose closed with it. It asked whether the dephase-letter-swap S₃ acting on the palindromizers Π is the **same** abstract S₃ as the letter-S₃ of this note acting on the diagonals Q. There were never two S₃'s. One letter-S₃ acts in both places, and it acts differently because the two objects have different sensitivities: it acts faithfully on the diagonals, which are quadratic in the letters and therefore blind to sign, and it carries Clifford signs on the palindromizers, giving a larger orbit there. The same asymmetry is what §4's closing subtlety records, where D fixes every diagonal but swaps Π_Z with Π_Y.

The arc `linear_s3_mirror_completion` is retired (`inspect --root arcs`); its resolution note is the anchor.

## Where it lives

- **Typed:** `ThreeDephasingDiagonalsOrbitClaim` (`compute/RCPsiSquared.Core/Symmetry/`, Tier 1 derived, dual parents `MirrorGroupD4Claim` + `AbsorptionTheoremClaim`, the physics edge that welds the mirror-group and absorption clusters).
- **Live:** the project's C# layer recomputes Q_X, Q_Y, Q_Z and the orbit on demand (`DiagonalWitness`, the "orbit" node; for codebase users, the CLI command `inspect --root diagonal`).
- **Verifier:** `simulations/one_diagonal_mirror_group.py` (self-validating: the Y-transpose gate, the letter orbit, the commutator pattern of §4). The non-existence of the order-48 completion has its own gate, `simulations/linear_s3_mirror_closure.py`, at the R·D proof.
