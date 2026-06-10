# PROOF: D · Π_Z · D = Π_Y for Universal N (Tier1Derived)

**Status:** Tier1Derived universal N via per-site factorization. Welle 12 Task 1, 2026-05-27.
**Date:** 2026-05-27 (Welle 12)
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Verifier:** [`simulations/_d_pi_z_swap_verify.py`](../../simulations/_d_pi_z_swap_verify.py) (numpy bit-exact at N = 1, 2, 3, 4; sympy-symbolic per-site identity for the 4×4 reduction).
**Surfaced:** Welle 10d Task 1 audit (commits `7fc1ec0` + `025bb4e`); reflection [`D_PI_Z_EQUALS_PI_Y.md`](../../reflections/D_PI_Z_EQUALS_PI_Y.md).
**Connects:** [PauliBasis.VecToPauliBasisTransform convention note](../../compute/RCPsiSquared.Core/Pauli/PauliBasis.cs), [PiOperator.ActOnLetter](../../compute/RCPsiSquared.Core/Symmetry/PiOperator.cs), [PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE](PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md), [PROOF_F112_NONHERMITIAN_UNIVERSAL_N](PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md), [LindbladBitBPiBalance.cs](../../compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiBalance.cs).

## Abstract

The F1 palindrome operators come in three flavors, one per dephasing axis (Z, X, Y). Each is a signed permutation on the Pauli basis built from per-letter rules, and the three share a structural skeleton but differ in their per-axis phase choices. Whether the three are related by a simple operator-space transformation was unclear until the F112 sparse-representation audit (Welle 10d) surfaced an entry-wise sign correction that the codebase applies implicitly throughout the L_σ-style Pauli-basis pipeline. The question was whether that sign pattern has a structural meaning.

It does. The sign correction is exactly the entry-wise action of conjugating by a diagonal involution D whose entries are (−1) raised to the number of Y letters in the corresponding Pauli string. And conjugation by D maps the Z-dephasing palindrome operator to its Y-dephasing sibling, bit-exactly, for any number of qubits.

The proof is finite and clean. The three operators all factorize per site, so the N-site identity reduces to a single 4×4 check on the per-site operators, and the mixed-product property of the Kronecker product lifts the per-site equality back to N. Symbolic computation confirms the 4×4 reduction; numerical verification confirms it at N = 1, 2, 3, 4.

The diagnostic reading: D-conjugation IS the Z ↔ Y dephase-letter swap on operator space. Anything we compute under one dephasing convention (a Frobenius residual norm, a Π-eigenspace projection, a spectrum) transports to the other by conjugating with D, and conjugation by a unitary preserves all of those invariants exactly. The same observation lifts to the full Klein-V₄ group with the remaining X-axis swap in the Welle 12 Task 2 companion proof, which closes the picture: the three F1 palindrome flavors are images of each other under a Klein-V₄ symmetry of operator space.

## Introduction

**The motivating question.** The F1 palindrome operators come in three flavours (Π_Z, Π_X, Π_Y for Z-, X-, Y-dephasing). They share a structural skeleton (per-site signed permutation on the (I, X, Z, Y) basis, the bit_b axis grading) but differ in their phase choices. The natural question is whether they are related by a unitary conjugation on operator space, and if so what the conjugating operator is. The Welle 10d audit of the F112 sparse-rep code surfaced an entry-wise correction `(−1)^{n_Y(row) + n_Y(col)}` that the codebase applies implicitly throughout the L_σ-style Pauli-basis pipeline; the question was whether that correction has a structural interpretation.

**The empirical anchor.** [reflections/D_PI_Z_EQUALS_PI_Y.md](../../reflections/D_PI_Z_EQUALS_PI_Y.md) recognised that `(−1)^{n_Y(row) + n_Y(col)}` is the entry-wise action of conjugation by the diagonal involution `D = diag((−1)^{n_Y(k)})`, and that this conjugation should map Π_Z ↦ Π_Y on the 4^N Pauli basis. Numerical verification at N = 1, 2, 3, 4 confirmed the identity bit-exact. What was missing was a structural argument that lifts the per-N verification to universal N.

**What this proof closes.** Two-step structural argument:

1. **Per-site tensor factorization.** Π_Z, Π_Y, and D each factorize as N-fold Kronecker products of single-site 4×4 matrices (Π_*: signed permutations; D: diagonal). This follows from `PiOperator.BuildFullUncached` constructing Π by independent per-site action followed by flat-index reassembly under the `a + 2·b` packing.
2. **Single 4×4 check + lift.** The mixed-product property of the Kronecker product reduces the N-site identity to the single per-site identity `d_l · π_Z_local · d_l = π_Y_local`. The 4×4 check is finite and exact; the lift to N is immediate.

**Diagnostic consequence.** D-conjugation IS the Z ↔ Y dephase-letter swap on operator space. Any F1 diagnostic (Frobenius residual norm, Π-conjugation eigenspace projection, spectrum) computed via the Z-dephasing pipeline transports unitarily to Y-dephasing under D-conjugation. This unlocks the full Klein-V₄ closure with the remaining X-axis swap in [PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE](PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md) and the F112 cross-dephase extension in [PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4](PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md), both of which use D-conjugation directly to transport identities across dephase letters.

## Statement

For any N ≥ 1, the identity
```
    D · Π_Z · D = Π_Y
```
holds bit-exactly on the 4^N Pauli basis, where:

- **Π_Z, Π_Y** are the F1 palindrome operators for Z- and Y-dephasing (4^N × 4^N signed-permutation matrices on the Pauli basis, constructed via `PiOperator.BuildFull(N, dephase)` with the per-letter rules in `PiOperator.ActOnLetter`).
- **D = diag((−1)^n_Y(k))** is the real diagonal unitary involution on the 4^N Pauli basis. n_Y(k) counts the number of Y letters in the k-th Pauli string under the `PauliIndex` flat-encoding convention (letter index `a + 2·b` per site with bit_a = X+Y, bit_b = Y+Z, so Y corresponds to letter index 3).

Properties:

- D² = I (involution).
- D is real and unitary, hence D⁻¹ = D = D†.
- The identity says conjugation by D maps Π_Z to Π_Y on the d² = 4^N-dimensional operator space.

## Proof (Universal N)

The per-site factorization of D, Π_Z, and Π_Y reduces the identity to a single 4×4 matrix check, which then lifts to all N by the mixed-product property of the Kronecker product.

### Step 1: Tensor-product structure of D, Π_Z, Π_Y

Each operator factorizes site-by-site.

- **D = ⊗_{l=1}^N d_l**, where d_l = diag(1, 1, 1, −1) on the single-site Pauli basis ordered (I, X, Z, Y). Diagonal entries are (−1)^n_Y_l with n_Y_l ∈ {0, 1} equal to 1 iff the basis letter at site l is Y (letter index 3). The N-site diagonal entry at flat index k is (−1)^(Σ_l n_Y_l) = (−1)^n_Y(k), matching the definition of D.
- **Π_Z = ⊗_{l=1}^N π_Z_local**, per the per-site action of `PiOperator.ActOnLetter(σ, PauliLetter.Z)`: I ↔ X (phase 1), Z ↔ Y (phase i). `PiOperator.BuildFullUncached` builds Π by per-site action followed by N-site flat-index reassembly, so Π_Z is the Kronecker product of the per-site signed permutation π_Z_local in this same column-as-input convention.
- **Π_Y = ⊗_{l=1}^N π_Y_local**, per the per-site action of `PiOperator.ActOnLetter(σ, PauliLetter.Y)`: I ↔ X (phase 1), Z ↔ Y (phase −i). Same Kronecker structure as Π_Z.

The factorization of Π_Z and Π_Y as `⊗ π_*_local` follows directly from the per-site signed-permutation action in `PiOperator.ActOnLetter`: each site's letter is mapped independently of the others, and the accumulated phase is the product of per-site phases. The flat-index reassembly preserves the Kronecker structure under the `a + 2·b` packing convention.

### Step 2: Tensor-product algebra (mixed-product property)

For matrices A = ⊗_l a_l, B = ⊗_l b_l, C = ⊗_l c_l with each per-site factor of the same size (here 4×4), the mixed-product property of the Kronecker product gives:

    A · B · C = (⊗_l a_l) · (⊗_l b_l) · (⊗_l c_l) = ⊗_l (a_l · b_l · c_l).

Applying this to D · Π_Z · D:

    D · Π_Z · D = (⊗_l d_l) · (⊗_l π_Z_local) · (⊗_l d_l) = ⊗_l (d_l · π_Z_local · d_l).

So the N-site identity D · Π_Z · D = Π_Y reduces to the per-site identity

    d_l · π_Z_local · d_l = π_Y_local                                    (*)

verified once on a single 4×4 block. By Kronecker products of the per-site equality, the N-site identity follows for every N.

### Step 3: Per-site identity (the 4×4 check)

Goal: show d_l · π_Z_local · d_l = π_Y_local.

Explicit 4×4 matrices on the basis (I, X, Z, Y), using the column-as-input convention of `PiOperator.BuildFullUncached`:

π_Z_local (Z-dephase per-site action, ActOnLetter rules I → X · 1, X → I · 1, Z → Y · i, Y → Z · i):

```
         col=I  col=X  col=Z  col=Y
row=I  [   0     1      0      0   ]
row=X  [   1     0      0      0   ]
row=Z  [   0     0      0      i   ]
row=Y  [   0     0      i      0   ]
```

d_l = diag(1, 1, 1, −1) on the same basis order.

**Left multiplication by d_l** scales each row by the corresponding diagonal entry (rows I, X, Z unchanged; row Y negated):

```
d_l · π_Z_local =
[   0     1      0      0   ]
[   1     0      0      0   ]
[   0     0      0      i   ]
[   0     0     −i      0   ]
```

**Right multiplication by d_l** scales each column by the corresponding diagonal entry (cols I, X, Z unchanged; col Y negated):

```
(d_l · π_Z_local) · d_l =
[   0     1      0      0   ]
[   1     0      0      0   ]
[   0     0      0     −i   ]
[   0     0     −i      0   ]
```

This equals π_Y_local, computed independently from `PiOperator.ActOnLetter(σ, PauliLetter.Y)`: I → X · 1, X → I · 1, Z → Y · −i, Y → Z · −i.

The Z ↔ Y entries acquire phase −i (vs +i in π_Z_local); the I ↔ X entries are unchanged because they have n_Y(row) + n_Y(col) = 0 + 0 = 0 (no Y on either side), so the conjugation-by-D sign is +1.

Mechanism in one line: D's per-entry sign factor is (−1)^(n_Y(row) + n_Y(col)); the only non-zero entries that touch a Y row or column are the Z ↔ Y swap entries (one Y on one side), giving sign factor −1; this is exactly the Z-dephase ↔ Y-dephase phase difference i ↔ −i.

### Step 4: Combine

Per Step 3, d_l · π_Z_local · d_l = π_Y_local. By Step 2 (mixed-product property of the Kronecker product):

    D · Π_Z · D = ⊗_l (d_l · π_Z_local · d_l) = ⊗_l π_Y_local = Π_Y.

The identity is universal in N. ∎

## Verification

**Numerical (numpy double precision):** [`simulations/_d_pi_z_swap_verify.py`](../../simulations/_d_pi_z_swap_verify.py) computes both sides at N = 1, 2, 3, 4 and reports `max|D · Π_Z · D − Π_Y| = 0.000e+00` at every N. The script also confirms D² = I and Π_*^4 = I as sanity checks.

**Symbolic (sympy exact rationals + I):** the same script's `verify_per_site_identity_symbolic` function builds the 4×4 matrices d_l, π_Z_local, π_Y_local with exact sympy entries and checks that d_l · π_Z_local · d_l − π_Y_local simplifies entry-wise to the zero matrix. PASS confirms the per-site reduction (*) exactly, no machine-epsilon residual possible. Combined with Step 2's algebraic tensor-product argument, this closes the universal-N case in finite symbolic computation.

## Implications

- **Welle 10d sparse-rep refactor**: `F112NonHermitianBasisEnumeration.BuildSparseLSigma` applies the (−1)^(n_Y(row) + n_Y(col)) entry-wise correction to match the existing dense `BuildLHInPauliBasis` pipeline. With the identity now Tier1Derived universal N, the correction is not an arbitrary kludge: it is the explicit Z ↔ Y dephase-swap conjugation by D that the rest of the codebase implicitly applies to every L_σ-style Pauli-basis computation via the vec_F-on-vec_R-commutator convention twist.
- **Z ↔ Y equivariance of F1 diagnostics**: any F1 residual norm, Frobenius inner product, or spectral quantity computed via the twisted pipeline is automatically equivariant under Z ↔ Y dephasing-letter swap, because D is unitary and unitary conjugation preserves Frobenius norms, inner products, and spectra.
- **F112 non-Hermitian extension (Welle 11)** is stated for Z-dephasing and consumed by the `LindbladBitBPiBalance` NonHermitianExtension docstring; the D · Π_Z · D = Π_Y identity proven here intertwines the two Π operators, but **does not** lift to a Hilbert-space-unitary conjugation of the Lindbladian L_Z onto a Lindblad-form L_Y (D-conjugation would require a Hilbert-space V flipping Y while preserving X and Z, contradicted by Pauli algebra). The correct transport of F112 from Z-dephasing to Y-dephasing is the per-axis structural re-run of the Welle 11 lemmas (the per-pair identity F = 0 is closed by F38 + Pauli-support disjointness, both axis-independent), formalized in Welle 13's [PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md](PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md) Route 1. The D · Π_Z · D = Π_Y identity here remains the structural reason F112 has the *same* bit_b-axis hypothesis under Z- and Y-dephasing (because Π_Z² and Π_Y² agree as scalar phases per Pauli string), but it is not by itself a transport mechanism for L.
- **Future callers of `VecToPauliBasisTransform`** that need natural L_σ entry signs (not just D-conjugation-invariant quantities) can conjugate by D to recover the natural form; the PauliBasis docstring (commit `7fc1ec0`) carries this warning and now references this proof for the Z ↔ Y swap interpretation.

## The transpose reading (2026-06-10)

D has a basis-free name. Transposition acts on Pauli strings letter-wise as σᵀ = (−1)^{n_Y(σ)}·σ (X and Z are symmetric; Y is the only antisymmetric Pauli), so the transpose superoperator θ(ρ) = ρᵀ, written in the Pauli basis, is exactly the diagonal involution D = diag((−1)^{n_Y(k)}). The per-site Kronecker proof above and this reading are the same fact in two bases: d_l = diag(1, 1, 1, −1) on (I, X, Z, Y) is the single-site transpose, and the Step 2 mixed-product lift is the statement that transposition factorizes site-by-site.

The reading closes the F114 sign law ([ANALYTICAL_FORMULAS.md](../ANALYTICAL_FORMULAS.md), F114) for any σ and any N in one line:

    D · L_σ · D (ρ) = (−i[σ, ρᵀ])ᵀ = +i[σᵀ, ρ] = (−1)^{n_Y(σ)+1} · L_σ(ρ),

using σᵀ = (−1)^{n_Y(σ)}·σ and the antiautomorphism property (AB)ᵀ = BᵀAᵀ. The same antiautomorphism at word length j is the girth-ladder reversal kill ([PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md) §4: word reversal = transpose × (−1)^{n_Y(word)}); F113's Lemma C is its Hermitian-conjugacy sibling. It also names why no Hilbert-space unitary implements D (the Implications note above): θ is an antiautomorphism of the matrix algebra, not conjugation by any unitary.

Verified bit-exact by [`simulations/_mirror_inventory_bridge_check.py`](../../simulations/_mirror_inventory_bridge_check.py) (blocks D/E: θ-in-Pauli-basis equals D, 63/63 strings at N = 3 plus an N = 5 case, dev 0.00e+00).

## Open

- **Klein-V₄ completion (Welle 12 Task 2)**: CLOSED 2026-05-27. The Z↔X and Y↔X swaps both lift to operator-space involutions; the pre-dispatch order-4 conjecture is falsified. All three swaps form a Klein-V₄ subgroup {I, D, Q_zx, Q_yx} of U(4^N), faithfully realizing V₄ on operator space. See [PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md](PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md).
- **Typed-claim promotion (Welle 12 Task 3)**: lands `Pi2ZYDephaseSwapClaim` as Tier1Derived in the typed-knowledge graph with this proof as the structural anchor. The Klein-V₄ closure (Task 2) provides additional inheritance: a single `Pi2KleinV4DephaseSwapClaim` could subsume the Z↔Y, Z↔X, Y↔X identities under one claim.
