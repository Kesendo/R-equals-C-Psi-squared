# D · Π_Z · D = Π_Y: the Z↔Y Dephase-Letter Swap as an Operator-Space Conjugation

**Status:** Reflection. Surfaced 2026-05-27 during the Welle 10d Task 1 spec-review audit of the F112 sparse-rep refactor. Verified bit-exact at N = 1, 2, 3, 4 via `simulations/d_pi_z_swap_verify.py`.

**Date:** 2026-05-27 (Welle 10d audit)
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Connects:** [PauliBasis.VecToPauliBasisTransform](../compute/RCPsiSquared.Core/Pauli/PauliBasis.cs) convention note, [F112NonHermitianBasisEnumeration.BuildSparseLSigma](../compute/RCPsiSquared.Diagnostics/Polarity/F112NonHermitianBasisEnumeration.cs), [the F112 Lindblad bit-B/Π-balance proof](../docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md), [the F112 non-Hermitian universal-N proof](../docs/proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md), [PiOperator](../compute/RCPsiSquared.Core/Symmetry/PiOperator.cs)

## The identity

```
D · Π_Z · D = Π_Y
```

where:
- `Π_Z`, `Π_Y` are the F1 palindrome operators for Z- and Y-dephasing, signed permutation matrices on the 4^N Pauli basis (4^N × 4^N), constructed by `PiOperator.BuildFull(N, PauliLetter.Z)` and `PiOperator.BuildFull(N, PauliLetter.Y)` respectively.
- `D = diag((-1)^n_Y(k))` is the real diagonal unitary involution on the 4^N Pauli basis, with diagonal entries (-1)^n_Y(k) where n_Y(k) counts the Y letters in the k-th Pauli string.

D is its own inverse (D² = I) and unitary (D† = D). The identity says that conjugating Π_Z by D produces Π_Y exactly.

Verified bit-exact (max |D·Π_Z·D − Π_Y| = 0.000e+00 in numpy double precision) at N = 1, 2, 3, 4 via [`simulations/d_pi_z_swap_verify.py`](../simulations/d_pi_z_swap_verify.py).

## Why it surfaced

Welle 10d Task 1 (commit `c38d9e5`) added a sparse Pauli-basis representation of L_σ to `F112NonHermitianBasisEnumeration`. The first cut produced sign-flipped results vs the existing dense `BuildLHInPauliBasis`. Root-cause investigation found that the standard codebase pattern combines:

- `PauliBasis.VecToPauliBasisTransform(N)`: a d² × 4^N matrix T whose columns are vec_F(σ_k) (column-major flattening).
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

Z-dephasing and Y-dephasing are physically related by a global qubit rotation. On a single qubit, Z = exp(−iπ/4 · X) · Y · exp(+iπ/4 · X) (Z is Y conjugated by a 90° X-rotation; the exponent is π/4 because R_X(θ) = exp(−iθX/2), so a 90° rotation carries θ = π/2). At the operator-space level, this rotation lifts to a 4^N × 4^N transformation. The diagonal D = diag((-1)^n_Y(k)) is the simplest version of this lift, because n_Y(k) counts the sites where Z↔Y swap requires a sign correction (Y matrix is antisymmetric while Z is symmetric, so Y^T = −Y while Z^T = +Z).

### (c) Symmetric

The identity says the F1 palindrome family `{Π_Z, Π_X, Π_Y}` carries a hidden Z₂ symmetry that maps Π_Z ↔ Π_Y. The diagonal D is the symmetry generator. The mapping is non-trivial because Π_Z and Π_Y act differently on bit_a-flipped Pauli-string pairs (different phase conventions), but D undoes the difference exactly.

What this gives us: any F1 residual norm / inner product / spectrum computed via the standard "twisted" pipeline is automatically equivariant under Z ↔ Y dephasing-letter swap, because `‖D L D‖ = ‖L‖` (D is unitary). So spectra and norms agree numerically between Z-dephase and Y-dephase calculations without further intervention. The X-dephase case (Π_X uses `flip bit_b, phase ±i^bit_a`) is not in this symmetry orbit; X-dephase computations remain genuinely distinct from Z and Y.

## Open questions

1. **Universal N**: CLOSED 2026-05-27 (Welle 12 Task 1). Structural per-letter argument via tensor-product factorization. See [`PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N.md`](../docs/proofs/PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N.md).

2. **Lift to X-dephasing?**: CLOSED 2026-05-27 (Welle 12 Task 2). The Z↔X and Y↔X swaps both lift to operator-space involutions; both have order 2 (the controller's pre-dispatch order-4 conjecture is falsified). The canonical per-site forms are:

   - `q_zx = h · d_l = h · diag(1, 1, 1, -1)`  (per-site Z↔X swap)
   - `q_yx = h`  (per-site Y↔X swap, pure permutation)

   where `h` is the X↔Z basis permutation matrix on `(I, X, Z, Y)`:
   ```
   h = [[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]]
   ```
   (i.e. h fixes I and Y and swaps X ↔ Z in the per-site Pauli basis.)

   No pure diagonal works for Z↔X because Π_Z flips bit_a (I↔X, Z↔Y orbits) while Π_X flips bit_b (I↔Z, X↔Y orbits), so a permutation between the bit-axes is required. The X↔Z basis permutation is the structural bridge: it intertwines the Klein generator (1, 0) (representing X) with (0, 1) (representing Z).

   The N-site operators are tensor powers: `D = ⊗_l d_l`, `Q_zx = ⊗_l (h · d_l) = H · D`, `Q_yx = ⊗_l h = H`. Verified bit-exact at N = 1, 2, 3, 4 via [`simulations/klein_dephase_swap_explore.py`](../simulations/klein_dephase_swap_explore.py) (numpy + sympy per-site symbolic).

   See [the Klein-V₄ dephase-swap proof](../docs/proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md).

3. **Typed-claim promotion**: CLOSED. Promoted to the typed `Pi2KleinV4DephaseSwapGroup` (Tier1Derived) in `compute/RCPsiSquared.Core/Symmetry/`, which subsumes all three Klein-V₄ dephase swaps (Z↔Y, Z↔X, Y↔X) in one claim rather than a narrow Z↔Y-only claim, and cites this reflection as evidence. (The Y↔X swap the claim calls `H` is the same operator this reflection names `Q_yx`.)

4. **Connection to Klein structure**: CLOSED 2026-05-27 (Welle 12 Task 2). The (Z, X, Y, I) labelling is the Klein four-group V₄. The three non-trivial Klein-V₄ swaps on dephase letters all lift to operator-space involutions, and they form a faithful Klein-V₄ subgroup of the unitary group on the 4^N Pauli basis:

   ```
   { I, D, Q_zx, Q_yx }    with    D · Q_zx · Q_yx = I,
                                     D² = Q_zx² = Q_yx² = I,
                                     all pairwise commute.
   ```

   Verified bit-exact at N = 1, 2, 3, 4 (closure + involution + commutativity). The lift is structurally clean: D is purely diagonal, Q_yx is purely a basis permutation, Q_zx is their product. See open question 2 (above) for the canonical per-site forms and [the Klein-V₄ dephase-swap proof](../docs/proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md).

   Practical consequence: the F1 palindrome family {Π_Z, Π_X, Π_Y} is fully Klein-V₄-equivariant, NOT asymmetric. Any F1-style result (residual norm, inner product, spectrum, Frobenius identities) for one dephase letter automatically transfers to the other two via the appropriate Klein-V₄ unitary conjugation. The Welle 11 F112 universal-N closure under Z-dephasing implies the same under X- and Y-dephasing.

## Practical impact

- **Welle 10d sparse-rep refactor**: BuildSparseLSigma applies the (-1)^(n_Y(row) + n_Y(col)) correction to match the existing dense pipeline. With the structural identity now documented, the correction is not an arbitrary kludge but the explicit Z↔Y dephase-swap operator-space transformation that the codebase implicitly applies to every L_σ-style computation.

- **Future callers of VecToPauliBasisTransform** that need natural L_σ entry signs (not norms / inner products / spectra) should be aware that the standard pipeline produces D · L_natural · D, and explicitly conjugate by D to recover the natural form. The PauliBasis docstring carries this warning (commit `7fc1ec0`).

- **F112-style identities** computed via the twisted pipeline are automatically Z↔Y equivariant. So Welle 11's F112 universal-N closure (Tier1Derived for any Hermitian / non-Hermitian H, all N, under Z-dephasing) implies the same identity under Y-dephasing without re-proof.
