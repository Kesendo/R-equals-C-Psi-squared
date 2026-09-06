# PROOF: F1 palindrome residual closed form under depolarizing noise

**Status:** Tier 1 derived. Closed-form Frobenius norm of the F1 residual M for the depolarizing channel, derived analytically from the per-site action of Π and the depolarizing dissipator; verified to machine precision at N = 2, 3, 4, 5 across uniform and non-uniform γ.
**Date:** 2026-05-18
**Authors:** Thomas Wicht, Claude (Opus 4.7)

## Abstract

The F1 palindrome theorem closes for pure Z-dephasing. Every other physical noise channel breaks the closure to some degree, and each break leaves a structural fingerprint in the residual M. T1 amplitude damping leaves one fingerprint (the closed form proven in the sibling document); depolarizing noise leaves another, and this proof writes it down.

The centered depolarizing fingerprint is cleaner than T1's: it is purely local and graph-independent. The uncentered sum ΠLΠ⁻¹+L does contain a cooperative 16(Σγ)² term, but that term is exactly the squared diagonal mean. The F1 shift +2Σγ·I removes it, leaving only (16/9)Σγ_l².

A second structural surprise is that the depolarizing residual is purely symmetric: the algebraic Π-conjugation identity gives M_anti = 0. T1 by contrast has a nonzero M_anti carrying its σ⁻ off-diagonal content, but it also has a substantial symmetric component. Depolarizing noise is diagonal in the Pauli basis (each Pauli letter goes to itself with a scaling factor), and the diagonal is preserved by Π, so its symmetric / anti-symmetric split is trivial.

This proof owns the pure-depolarizing block and its orthogonality to the Hamiltonian and Z-dephasing blocks. It does **not** make the T1 and depolarizing norms additive: the unshifted T1 residual and the F1-centered depolarizing residual have a nonzero Frobenius cross-term. Consequently the three separate closed forms are not, by themselves, an inversion dictionary for a mixed T1-plus-depolarizing channel, and M_anti alone does not recover all T1 content.

## Statement

Let L_depol be the Lindbladian for the per-site depolarizing channel with rates {γ_l}, no Hamiltonian and no Z-dephasing:

    L_depol(ρ) = Σ_l γ_l · (1/3) · Σ_{P ∈ {X,Y,Z}} (P_l ρ P_l − ρ).

Let Π be the F1 palindrome operator for Z-dephasing (see [the Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md)). Define the centered F1 residual

    M_F1 := Π · L_depol · Π⁻¹ + L_depol + 2(Σ_l γ_l)·I.

In the framework's orthonormal Pauli-string basis (the basis used by `palindrome_residual` in [`framework/lindblad.py`](../../simulations/framework/lindblad.py)):

**Theorem (depolarizing block closed form).** The pure-depol residual satisfies

    ‖M_F1(depol)‖²_F = 4^(N−1) · (16/9) · Σ_l γ²_l.

Three structural facts follow immediately and are verified to machine precision:

1. **Per-site only.** No graph-parameter (B, D2) dependence: the centered residual scales only with `Σγ²`; its cross-site term is zero.
2. **Hamiltonian-independent.** L_depol does not couple to H, and the depol-block is Frobenius-orthogonal to the H-block (same proof structure as T1 in [F1 residual under T1 damping](PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md) Step 6).
3. **Π²-decomposition is trivial.** M is Pauli-basis-diagonal, so Π·M·Π⁻¹ = M exactly (not just Π²·M·Π²⁻¹ = M). Hence M_anti = 0 and ‖M_sym‖² = ‖M‖². Contrast T1 where M_anti = D_{T1, odd} carries the F82/F84 amplitude-damping content; depol has no σ⁻-style off-diagonal Pauli-basis channel.

## Conventions

- **Pauli letters** are indexed (a, b) ∈ {(0,0), (1,0), (0,1), (1,1)} = (I, X, Z, Y) following the framework's Klein-Vierergruppe convention ([`framework/pauli.py`](../../simulations/framework/pauli.py)). `bit_a` is the X/Y indicator, `bit_b` is the Y/Z indicator.
- **Pauli-string basis on N sites** is the 4^N orthonormal basis {σ_α} with the inner product ⟨σ_α | σ_β⟩ = δ_{αβ}; equivalently `Tr(σ_α^† σ_β) / 2^N = δ_{αβ}`. Operators on operator-space (4^N × 4^N "super-operator" matrices) inherit this orthonormality via the Hilbert-Schmidt product `⟨A, B⟩ = Tr(A^† B)`. The framework's `palindrome_residual` enforces this via the transform `L_pauli = M_basis^† · L_vec · M_basis / 2^N` in [`framework/lindblad.py`](../../simulations/framework/lindblad.py).
- **Π** is the F1 palindrome operator for Z-dephasing, acting per site on Pauli letters as
      I ↔ X (phase +1),    Y ↔ Z (phase +i).
  See [the Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md) (the "Conjugation Operator Π" section). On the full Pauli string Π is the tensor product of per-site rules with phases multiplied. Π is unitary, order-4: Π² is diagonal with eigenvalue (−1)^{Σ_l bit_b(α_l)} on the string σ_α (see [F81](PROOF_F81_PI_CONJUGATION_OF_M.md) Step 1).
- **Depolarizing dissipator** uses the standard isotropic-Pauli convention. The per-site channel is
      D_{depol, l}(ρ) = (γ_l / 3) · Σ_{P ∈ {X, Y, Z}} (P_l ρ P_l − ρ).
  Equivalently in jump-operator form (used by `lindbladian_general`): three Lindblad operators per site, c_{l, P} = √(γ_l/3) · P_l for P ∈ {X, Y, Z}; D[c]ρ = c ρ c† − ½{c†c, ρ} reduces to (γ_l/3) · (P_l ρ P_l − ρ) since each P is Hermitian and squares to I. Same normalization as F5 in [Why Depolarizing Noise Breaks the Palindrome](../../experiments/DEPOLARIZING_PALINDROME.md).

## Empirical anchor (motivation, observed earlier)

[F5](../ANALYTICAL_FORMULAS.md#f5-depolarizing-error-tier-1-proven) records the extreme pair-sum shortfall `error = (2/3) · Σ_l γ_l` for the depolarizing channel, Hamiltonian-independent and linear in γ ([Why Depolarizing Noise Breaks the Palindrome](../../experiments/DEPOLARIZING_PALINDROME.md)). In the centered diagonal residual this is its spectral norm. The Frobenius norm derived here is a different diagnostic of the same split, not a trace projection of the bare residual. The earlier `F1OpenQuestions` item "depolarizing noise: residual scaling" is closed by this document, promoting the entry to the Tier-1-derived [`F1DepolResidualClosedForm`](../../compute/RCPsiSquared.Core/F1/F1DepolResidualClosedForm.cs) claim on `F1KnowledgeBase`.

## Proof

### Step 1: Depolarizing dissipator's per-site Pauli-basis matrix

In the single-qubit orthonormal Pauli basis (I, X, Y, Z), the depolarizing dissipator D_{depol, local} is a 4 × 4 matrix with `D_{αβ} = (1/2) Tr(σ_α · D_{depol, local}(σ_β))`. Direct computation: since P · σ_α · P − σ_α = −2 σ_α when {P, σ_α} = 0 and equals 0 when [P, σ_α] = 0, and exactly 2 of the 3 letters {X, Y, Z} anti-commute with each non-identity Pauli, we get

|       | I        | X        | Y        | Z        |
|-------|----------|----------|----------|----------|
| **I** | 0        | 0        | 0        | 0        |
| **X** | 0        | −4γ/3    | 0        | 0        |
| **Y** | 0        | 0        | −4γ/3    | 0        |
| **Z** | 0        | 0        | 0        | −4γ/3    |

(Take γ = 1 for the closed-form derivation; the formula will be γ²-homogeneous in the per-site contribution and γ_l γ_{l′}-bilinear in the cross-site contribution.)

The depolarizing dissipator is **Pauli-basis-diagonal**: each non-identity Pauli decays uniformly at rate 4γ/3, the identity is immune. This is the canonical "isotropic Pauli channel" structure ([Why Depolarizing Noise Breaks the Palindrome](../../experiments/DEPOLARIZING_PALINDROME.md) §2: "the 1:3 split").

### Step 2: Π conjugation on the per-site Pauli-basis matrix

Π acts on the single-letter Pauli basis as a 4 × 4 signed permutation:

    Π = ⎡ 0  1  0   0 ⎤
        ⎢ 1  0  0   0 ⎥
        ⎢ 0  0  0   i ⎥
        ⎣ 0  0  i   0 ⎦

(rows = target letter, columns = source letter, in the order I, X, Y, Z). Π is unitary, Π² = diag(+1, +1, −1, −1) (the bit_b-parity sign), and Π⁻¹ = Π†.

Conjugating the diagonal D_{depol, local} (with γ = 1) by Π permutes its eigenvalues by Π's underlying permutation: I ↔ X swaps entries `D_II ↔ D_XX`, Y ↔ Z swaps `D_YY ↔ D_ZZ`. The +i phases on Y ↔ Z cancel in conjugation (since `(+i) · (+i)* = 1`). Hence

    Π · D_depol · Π⁻¹ = diag(D_XX, D_II, D_ZZ, D_YY) = diag(−4/3, 0, −4/3, −4/3).

### Step 3: Bare per-site kernel and its F1 centering

Summing the two diagonal matrices (with γ = 1):

    M_l^bare = diag(−4/3, −4/3, −8/3, −8/3).

Four non-zero diagonal entries; no off-diagonal entries (contrast T1, which has off-diagonal (Y, X) and (Z, I) entries). Frobenius norm squared:

    ‖M_l^bare‖²_F = (4/3)² · 2 + (8/3)² · 2 = **160/9**.

Diagonal trace:

    tr(M_l^bare) = **−8**,    |tr(M_l^bare)|² = **64**.

The F1 shift contributes +2γ_l I on this site. At γ_l=1 the centered local kernel is

    M_l = M_l^bare + 2I = diag(+2/3, +2/3, −2/3, −2/3),

so `tr(M_l)=0` and `‖M_l‖²_F=4·(2/3)²=16/9`. The two values remain distinct, hence the palindrome is still broken; centering removes only the mean.

### Step 4: Multi-site assembly via per-site action

For the full N-qubit problem with γ_l per site, the multi-qubit M_depol is

    M_depol = Σ_l γ_l · M_l^{(super)}

where M_l^{(super)} acts as M_l on the Pauli-basis coordinate of site l and as the 4 × 4 identity on every other site:

    (M_l^{(super)})_{α β} = (M_l)_{α_l, β_l} · ∏_{k ≠ l} δ_{α_k, β_k}.

This per-site / identity-elsewhere structure follows because the depolarizing dissipator on site l acts only on the site-l qubit and the identity Pauli on every other site is fixed under the partial-trace projection that defines the per-site Pauli-basis component.

Equivalently in tensor notation: M_l^{(super)} = I_4^{⊗l} ⊗ M_l ⊗ I_4^{⊗(N−l−1)}.

### Step 5: Frobenius norm of the centered residual

For tensor-product matrices in the Pauli basis,

    tr( (M_l^{(super)})^† · M_{l′}^{(super)} ) = tr_{site l}(M_l^† · X_l) · tr_{site l′}(X_{l′} · M_{l′}) · ∏_{k ≠ l, l′} tr_{site k}(I_4)

where X_l, X_{l′} are the corresponding identity factors of the other operator at site l, l′. Concretely:

- **l = l′ (same site):** all other sites contribute tr(I_4) = 4, giving
      tr( (M_l^{(super)})^† · M_l^{(super)} ) = **(16/9) · 4^(N−1)**.
- **l ≠ l′ (different sites):** `tr(M_l)=0`, so
      tr( (M_l^{(super)})^† · M_{l′}^{(super)} ) = **0**.

Substituting into the Frobenius norm of M_depol:

    ‖M_F1,depol‖²_F = **4^(N−1) · (16/9) · Σ_l γ²_l**.    ∎

For comparison only, assembling the bare kernels before centering gives `4^(N−1)[(16/9)Σγ_l²+16(Σγ_l)²]`. The F1 shift removes the second term exactly.

### Step 6: Π²-decomposition is trivial (M_anti = 0)

The centered per-site M_l = diag(+2/3, +2/3, −2/3, −2/3) is **Pauli-basis-diagonal**. The multi-site M_depol = Σ_l γ_l · M_l^{(super)} is therefore also diagonal in the Pauli-string basis. Π swaps equal entries within the (I,X) and (Y,Z) pairs, so conjugation acts as the identity.

For our centered M_l, the diagonal pairs (I, X) both carry value +2/3, and the pairs (Y, Z) both carry value −2/3. Hence

    Π · M_l · Π⁻¹ = M_l    (exactly, per site)
    ⟹ Π · M_depol · Π⁻¹ = M_depol    (exactly, multi-site).

Therefore the Π²-orthogonal Pythagorean split of M_depol degenerates:

    M_anti(depol) = (M − Π·M·Π⁻¹) / 2 = 0,
    M_sym(depol)  = (M + Π·M·Π⁻¹) / 2 = M,
    ‖M_anti(depol)‖² = 0,
    ‖M_sym(depol)‖²  = ‖M_F1(depol)‖² = 4^(N−1) · (16/9)·Σγ².

This is the structural distinction from [F1 residual under T1 damping](PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md) Step 7: T1's M_l has off-diagonal entries (the (Z, I) channel from σ⁻ amplitude damping plus the (Y, X) entry produced by Π conjugation), giving a non-trivial Π²-anti-symmetric piece ‖M_anti(T1)‖² = 4^(N−1)·Σγ² that maps onto F82's D_{T1, odd}. Depolarizing has no such off-diagonal Pauli-basis channel: every Pauli decays into itself, so the per-site D_l is diagonal, Π conjugation merely permutes (already-equal) diagonal values, and the Π²-anti-symmetric piece vanishes identically.

A companion typed claim `F1DepolResidualPi2Decomposition` would be 150 lines for the single fact M_anti = 0. The point is inlined as an `ExtraChildren` node on the parent claim instead; see [`F1DepolResidualClosedForm`](../../compute/RCPsiSquared.Core/F1/F1DepolResidualClosedForm.cs).

### Step 7: F1 centering and relationship to F5

The F1 residual is conventionally written `M := Π·L·Π⁻¹ + L + 2σ·I` with σ chosen to absorb the "diagonal background" produced by the dissipator. For Z-dephasing σ = Σγ works exactly. `L_Z` itself is weight-dependent—on a Pauli string α its eigenvalue is `−2Σ_l γ_l bit_a(α_l)`—but Π complements that local bit, so `Π·L_Z·Π⁻¹ + L_Z = −2Σγ·I` uniformly and the shift cancels the sum.

For depolarizing, no scalar can eliminate the two-level split, but the F1 choice **σ=Σγ** removes its mean. The centered one-site entries are ±2γ/3, the minimum-norm scalar centering; the residual remains nonzero.

Numerical confirmation at N = 3, uniform γ = 0.1: ‖M‖² with σ = 0 is 23.893, while the centered F1 residual at σ = Σγ has norm squared 0.853333. The typed predictor uses the latter convention.

**F5 relation.** F5 records `error = (2/3)Σγ`, the extreme pair-sum shortfall. For the centered diagonal residual this is its spectral norm. `F1DepolResidualClosedForm` instead measures the squared Frobenius norm, `4^(N−1)(16/9)Σγ_l²`. They are distinct norms of the same 2:2-versus-1:3 obstruction.

## Verification

[`simulations/f1_depol_residual_verify.py`](../../simulations/f1_depol_residual_verify.py) verifies the closed form in seven sections:

1. **F1 sanity.** ‖M‖² for pure Z-dephasing is at machine precision (confirms the framework's Π is the right one).
2. **Pure depol numerical fit.** At N = 2, 3, 4, 5 (both uniform γ = 0.1 and non-uniform γ = [0.05, 0.10, ..., 0.05·N]) the centered fit is (a, b) = (16/9, 0); the prediction matches to ~10⁻¹³.
3. **Orthogonality H ⊥ depol for truly H.** Heisenberg H has ‖M‖² = 0; adding depol gives exactly ‖M(depol)‖² with zero cross-term.
4. **Orthogonality with Z-dephasing (H + Z + depol).** Adding Z-dephasing to the H + depol setup leaves ‖M‖² = ‖M(depol)‖² (since Z and H both contribute 0 to ‖M‖²).
5. **Orthogonality with soft (Π²-odd) H (XY+YX).** Soft H gives non-zero ‖M(H)‖² (the F49 closed-form prediction); adding depol gives exactly ‖M(H)‖² + ‖M(depol)‖² with cross-term ~10⁻¹³.
6. **Per-site kernel.** Displays the bare kernel (`tr=-8`, norm squared `160/9`) and its F1-centered form (`tr=0`, norm squared `16/9`), which kills cross-site terms.
7. **Π²-trivial split.** Asserts `‖M − Π·M·Π⁻¹‖_F < 1e-13` and `‖M − Π²·M·Π²⁻¹‖_F < 1e-13` at N = 2, 3, 4: M_anti(depol) = 0 exactly.

All verifications pass at machine precision. Summary of section-2 numerics:

| N | uniform γ = 0.1 | non-uniform [0.05·(k+1)] | fitted (a, b) |
|---|-----------------|--------------------------|----------------|
| 2 | obs 0.142222 = pred 0.142222 | obs 0.088889 = pred 0.088889 | (1.777778, 0.000000) |
| 3 | obs 0.853333 = pred 0.853333 | obs 0.995556 = pred 0.995556 | (1.777778, 0.000000) |
| 4 | obs 4.551111 = pred 4.551111 | obs 8.533333 = pred 8.533333 | (1.777778, 0.000000) |
| 5 | obs 22.755556 = pred 22.755556 | obs 62.577778 = pred 62.577778 | (1.777778, 0.000000) |

(1.777778 = 16/9 to 6 decimals.)

## Diagnostic interpretation

The closed form makes the centered F1 depol-block residual a quantitative, **Hamiltonian-independent and topology-independent** diagnostic for depolarizing content:

- **Pure-depol inversion (uniform γ).** `γ = √(‖M_F1‖²_F / [4^(N−1)·(16/9)·N])`.
- **Non-uniform limitation.** One centered Frobenius reading determines Σ_lγ_l², not the individual profile or Σ_lγ_l; those require independent measurements.
- **N-scaling.** At uniform γ the centered squared norm is `4^(N−1)·(16/9)·N·γ²`; there is no cross-site term.
- **No graph dependence.** Depolarizing is per-site only; the centered residual scales with `Σγ²`. Contrast the H-block, whose scaling can depend on graph data.

## Cross-references

### Repository entries

- **F1 palindrome equation** ([`docs/ANALYTICAL_FORMULAS.md` F1](../ANALYTICAL_FORMULAS.md#f1-palindrome-equation-tier-1-proven), [the Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md)): the underlying Π·L·Π⁻¹ + L + 2Σγ·I = 0 identity for Z-dephasing.
- **F5 depolarizing error** ([`docs/ANALYTICAL_FORMULAS.md` F5](../ANALYTICAL_FORMULAS.md), [Why Depolarizing Noise Breaks the Palindrome](../../experiments/DEPOLARIZING_PALINDROME.md)): the extreme pair-sum shortfall `(2/3)Σγ`, equivalently the spectral norm of the centered diagonal residual.
- **F49 Frobenius residual scaling** ([`docs/ANALYTICAL_FORMULAS.md` F49](../ANALYTICAL_FORMULAS.md#f49-cross-term-formula-tier-1-proven), [the Cross-Term Formula proof](PROOF_CROSS_TERM_FORMULA.md)): companion closed form for the Hamiltonian block.
- **F1 T1-residual closed form** ([F1 residual under T1 damping](PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md)): sibling closed form for amplitude damping; contrasts with depol via the Π²-decomposition non-triviality (T1's M_anti = D_{T1, odd}, depol's M_anti = 0).

### Typed claims

- [`compute/RCPsiSquared.Core/F1/F1DepolResidualClosedForm.cs`](../../compute/RCPsiSquared.Core/F1/F1DepolResidualClosedForm.cs): Tier-1-derived typed claim for this closed form, registered on `F1KnowledgeBase`. Replaces the earlier `F1OpenQuestions` item "depolarizing noise: residual scaling", which was closed by this proof on 2026-05-18.
- [`compute/RCPsiSquared.Core/F1/F1T1ResidualClosedForm.cs`](../../compute/RCPsiSquared.Core/F1/F1T1ResidualClosedForm.cs): sibling T1 closed-form claim.
- [`compute/RCPsiSquared.Core/Symmetry/F5DepolarizingErrorPi2Inheritance.cs`](../../compute/RCPsiSquared.Core/Symmetry/F5DepolarizingErrorPi2Inheritance.cs): F5 scalar error claim with Pi2-Foundation inheritance.

### Scripts

- [`simulations/f1_depol_residual_verify.py`](../../simulations/f1_depol_residual_verify.py): the verification script for this proof.

### Memory

- `project_palindrome_frobenius_scaling`: recorded the per-class dissipator Frobenius scaling pattern on 2026-04-29; depol arm closed 2026-05-18 by this proof.
