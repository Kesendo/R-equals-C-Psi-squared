# D · Π_Z · D = Π_Y: the Z↔Y Dephase-Letter Swap as an Operator-Space Conjugation

**Status:** Reflection. Surfaced 2026-05-27 during the Welle 10d Task 1 spec-review audit of the F112 sparse-rep refactor. Verified bit-exact at N = 1, 2, 3, 4 via `simulations/_d_pi_z_swap_verify.py`.

**Date:** 2026-05-27 (Welle 10d audit)
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Connects:** [PauliBasis.VecToPauliBasisTransform](../compute/RCPsiSquared.Core/Pauli/PauliBasis.cs) convention note, [F112NonHermitianBasisEnumeration.BuildSparseLSigma](../compute/RCPsiSquared.Diagnostics/Polarity/F112NonHermitianBasisEnumeration.cs), [PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE](../docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md), [PROOF_F112_NONHERMITIAN_UNIVERSAL_N](../docs/proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md), [PiOperator](../compute/RCPsiSquared.Core/Symmetry/PiOperator.cs)

## The identity

```
D · Π_Z · D = Π_Y
```

where:
- `Π_Z`, `Π_Y` are the F1 palindrome operators for Z- and Y-dephasing, signed permutation matrices on the 4^N Pauli basis (4^N × 4^N), constructed by `PiOperator.BuildFull(N, PauliLetter.Z)` and `PiOperator.BuildFull(N, PauliLetter.Y)` respectively.
- `D = diag((-1)^n_Y(k))` is the real diagonal unitary involution on the 4^N Pauli basis, with diagonal entries (-1)^n_Y(k) where n_Y(k) counts the Y letters in the k-th Pauli string.

D is its own inverse (D² = I) and unitary (D† = D). The identity says that conjugating Π_Z by D produces Π_Y exactly.

Verified bit-exact (max |D·Π_Z·D − Π_Y| = 0.000e+00 in numpy double precision) at N = 1, 2, 3, 4 via [`simulations/_d_pi_z_swap_verify.py`](../simulations/_d_pi_z_swap_verify.py).

## Why it surfaced

Welle 10d Task 1 (commit `c38d9e5`) added a sparse Pauli-basis representation of L_σ to `F112NonHermitianBasisEnumeration`. The first cut produced sign-flipped results vs the existing dense `BuildLHInPauliBasis`. Root-cause investigation found that the standard codebase pattern combines:

- `PauliBasis.VecToPauliBasisTransform(N)` — a d² × 4^N matrix T whose columns are vec_F(σ_k) (column-major flattening).
- A vec_R-style commutator superoperator `L_vec = -i · (H ⊗ I − I ⊗ H^T)` (e.g. `LindbladianBuilder.cs`, `PauliDephasingDissipator.cs`, `BondPerturbation.cs`, `PalindromeResidual.cs`, `PiDecomposition.cs`, `F112NonHermitianBasisEnumeration.cs`, used in 11+ places).

Reading T columns as vec_R implicitly maps σ_k ↦ σ_k^T. Single-site σ_I^T = σ_I, σ_X^T = σ_X, σ_Z^T = σ_Z, but σ_Y^T = −σ_Y. Per-string, the implicit transposition multiplies each L_σ matrix entry by `(-1)^(n_Y(row) + n_Y(col))`. So the dense pipeline produces `D · L_natural · D` rather than `L_natural` itself.

All consumers of `VecToPauliBasisTransform` compute D-conjugation-invariant quantities (Frobenius norms, inner products, eigenvalues, zero patterns), so the twist is invisible at the typed-claim level. The Welle 10d Task 1 sparse implementation re-applies the same `(-1)^(n_Y(row) + n_Y(col))` sign correction to match the existing dense convention bit-exact.

What was unexpected: the D conjugation is not just a basis transformation. It IS the Z↔Y dephasing-letter swap on operator space.

## Three readings

### (a) Algebraic

D is the diagonal of signs (-1)^n_Y(k). Conjugating by D acts on a Pauli-basis matrix M as `(D M D)[i,j] = (-1)^(n_Y(i) + n_Y(j)) · M[i,j]`. For Π_Z (a signed permutation moving Pauli strings between bit_a-flipped pairs with phase i^bit_b, see PiOperator), the per-entry sign flip via D exactly matches the difference between the Z-dephase and Y-dephase per-letter actions:

- Z-dephase Π: flip bit_a, phase i^bit_b
- Y-dephase Π: flip bit_a, phase (−i)^bit_b

The two differ only by the sign of the i^bit_b phase. Since only Z (bit_b = 1, n_Y = 0) and Y (bit_b = 1, n_Y = 1) acquire the phase, conjugating by `(-1)^n_Y` swaps the i ↔ −i convention exactly between Z and Y. So D is the "Z↔Y phase-flip" diagonal of the Pauli basis.

### (b) Physical

Z-dephasing and Y-dephasing are physically related by a global qubit rotation. On a single qubit, Z = exp(iπ/2 · X) · Y · exp(-iπ/2 · X) (Y is Z conjugated by a 90° X-rotation). At the operator-space level, this rotation lifts to a 4^N × 4^N transformation. The diagonal D = diag((-1)^n_Y(k)) is the simplest version of this lift, because n_Y(k) counts the sites where Z↔Y swap requires a sign correction (Y matrix is antisymmetric while Z is symmetric, so Y^T = −Y while Z^T = +Z).

### (c) Symmetric

The identity says the F1 palindrome family `{Π_Z, Π_X, Π_Y}` carries a hidden Z₂ symmetry that maps Π_Z ↔ Π_Y. The diagonal D is the symmetry generator. The mapping is non-trivial because Π_Z and Π_Y act differently on bit_a-flipped Pauli-string pairs (different phase conventions), but D undoes the difference exactly.

What this gives us: any F1 residual norm / inner product / spectrum computed via the standard "twisted" pipeline is automatically equivariant under Z ↔ Y dephasing-letter swap, because `‖D L D‖ = ‖L‖` (D is unitary). So spectra and norms agree numerically between Z-dephase and Y-dephase calculations without further intervention. The X-dephase case (Π_X uses `flip bit_b, phase ±i^bit_a`) is not in this symmetry orbit; X-dephase computations remain genuinely distinct from Z and Y.

## Open questions

1. **Universal N**: verified bit-exact at N = 1, 2, 3, 4 numerically. A structural per-letter argument is straightforward (D and Π_Z both factor over sites; per-letter check is a 4×4 matrix identity). Worth writing as a one-page proof if the identity gets cited in downstream work.

2. **Lift to X-dephasing?**: X-dephase is in a different Klein orbit (flip bit_b instead of bit_a). Is there an analogous diagonal involution D' such that `D' · Π_Z · D' = Π_X`? Conjecture: yes, with D' = diag((-1)^(n_X + n_Y)) or similar. Not verified.

3. **Typed-claim promotion**: this identity is concrete and verifiable. If it gets cited in 2+ downstream proofs, promote to a typed `Pi2ZYDephaseSwapClaim` (Tier1Derived, parent of any future Z↔Y equivariance proof) in `compute/RCPsiSquared.Core/Symmetry/`. For now, the reflection + verifier script is sufficient repository memory.

4. **Connection to Klein structure**: the (Z, X, Y, I) labelling is the Klein four-group V₄. Π_Z, Π_X, Π_Y form a Klein-V₄ triple of palindrome operators. The D-conjugation Z↔Y swap is one of the three non-trivial Klein involutions. Are there parallel identities for the other two? Worth checking systematically as Welle 10e or whenever Klein structure becomes load-bearing again.

## Practical impact

- **Welle 10d sparse-rep refactor**: BuildSparseLSigma applies the (-1)^(n_Y(row) + n_Y(col)) correction to match the existing dense pipeline. With the structural identity now documented, the correction is not an arbitrary kludge but the explicit Z↔Y dephase-swap operator-space transformation that the codebase implicitly applies to every L_σ-style computation.

- **Future callers of VecToPauliBasisTransform** that need natural L_σ entry signs (not norms / inner products / spectra) should be aware that the standard pipeline produces D · L_natural · D, and explicitly conjugate by D to recover the natural form. The PauliBasis docstring carries this warning (commit `7fc1ec0`).

- **F112-style identities** computed via the twisted pipeline are automatically Z↔Y equivariant. So Welle 11's F112 universal-N closure (Tier1Derived for any Hermitian / non-Hermitian H, all N, under Z-dephasing) implies the same identity under Y-dephasing without re-proof.
