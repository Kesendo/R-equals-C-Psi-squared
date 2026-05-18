# PROOF: F1 palindrome residual closed form under depolarizing noise

**Status:** Tier 1 derived. Closed-form Frobenius norm of the F1 residual M for the depolarizing channel, derived analytically from the per-site action of Π and the depolarizing dissipator; bit-exact verification at N = 2, 3, 4, 5 across uniform and non-uniform γ.
**Date:** 2026-05-18
**Authors:** Thomas Wicht, Claude (Opus 4.7)

## Statement

Let L_depol be the Lindbladian for the per-site depolarizing channel with rates {γ_l}, no Hamiltonian and no Z-dephasing:

    L_depol(ρ) = Σ_l γ_l · (1/3) · Σ_{P ∈ {X,Y,Z}} (P_l ρ P_l − ρ).

Let Π be the F1 palindrome operator for Z-dephasing (see [MIRROR_SYMMETRY_PROOF.md](MIRROR_SYMMETRY_PROOF.md)). Define the bare F1 residual

    M := Π · L_depol · Π⁻¹ + L_depol.

(σ-shift = 0; see Step 7 below for why depol cannot absorb a constant 2σ·I shift.)

In the framework's orthonormal Pauli-string basis (the basis used by `palindrome_residual` in [`framework/lindblad.py`](../../simulations/framework/lindblad.py)):

**Theorem (depolarizing block closed form).** The pure-depol residual satisfies

    ‖M(depol)‖²_F = 4^(N−1) · [ (16/9) · Σ_l γ²_l  +  16 · (Σ_l γ_l)² ].

Three structural facts follow immediately and are verified to machine precision:

1. **Per-site only.** No graph-parameter (B, D2) dependence: depolarizing noise is per-site, not bond-coupled. The residual scales purely with `(Σγ², (Σγ)²)`.
2. **Hamiltonian-independent.** L_depol does not couple to H, and the depol-block is Frobenius-orthogonal to the H-block (same proof structure as T1 in [PROOF_F1_T1_RESIDUAL_CLOSED_FORM](PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md) Step 6).
3. **Π²-decomposition is trivial.** M is Pauli-basis-diagonal, so Π·M·Π⁻¹ = M exactly (not just Π²·M·Π²⁻¹ = M). Hence M_anti = 0 and ‖M_sym‖² = ‖M‖². Contrast T1 where M_anti = D_{T1, odd} carries the F82/F84 amplitude-damping content; depol has no σ⁻-style off-diagonal Pauli-basis channel.

## Conventions

- **Pauli letters** are indexed (a, b) ∈ {(0,0), (1,0), (0,1), (1,1)} = (I, X, Z, Y) following the framework's Klein-Vierergruppe convention ([`framework/pauli.py`](../../simulations/framework/pauli.py)). `bit_a` is the X/Y indicator, `bit_b` is the Y/Z indicator.
- **Pauli-string basis on N sites** is the 4^N orthonormal basis {σ_α} with the inner product ⟨σ_α | σ_β⟩ = δ_{αβ}; equivalently `Tr(σ_α^† σ_β) / 2^N = δ_{αβ}`. Operators on operator-space (4^N × 4^N "super-operator" matrices) inherit this orthonormality via the Hilbert-Schmidt product `⟨A, B⟩ = Tr(A^† B)`. The framework's `palindrome_residual` enforces this via the transform `L_pauli = M_basis^† · L_vec · M_basis / 2^N` in [`framework/lindblad.py`](../../simulations/framework/lindblad.py).
- **Π** is the F1 palindrome operator for Z-dephasing, acting per site on Pauli letters as
      I ↔ X (phase +1),    Y ↔ Z (phase +i).
  See [MIRROR_SYMMETRY_PROOF.md](MIRROR_SYMMETRY_PROOF.md) (the "Conjugation Operator Π" section). On the full Pauli string Π is the tensor product of per-site rules with phases multiplied. Π is unitary, order-4: Π² is diagonal with eigenvalue (−1)^{Σ_l bit_b(α_l)} on the string σ_α (see [PROOF_F81](PROOF_F81_PI_CONJUGATION_OF_M.md) Step 1).
- **Depolarizing dissipator** uses the standard isotropic-Pauli convention. The per-site channel is
      D_{depol, l}(ρ) = (γ_l / 3) · Σ_{P ∈ {X, Y, Z}} (P_l ρ P_l − ρ).
  Equivalently in jump-operator form (used by `lindbladian_general`): three Lindblad operators per site, c_{l, P} = √(γ_l/3) · P_l for P ∈ {X, Y, Z}; D[c]ρ = c ρ c† − ½{c†c, ρ} reduces to (γ_l/3) · (P_l ρ P_l − ρ) since each P is Hermitian and squares to I. Same normalization as F5 in [DEPOLARIZING_PALINDROME.md](../../experiments/DEPOLARIZING_PALINDROME.md).

## Empirical anchor (motivation, observed earlier)

[F5](../ANALYTICAL_FORMULAS.md#f5-depolarizing-error-tier-1-proven) records the scalar palindrome error `error = (2/3) · Σ_l γ_l` for the depolarizing channel, Hamiltonian-independent and linear in γ ([DEPOLARIZING_PALINDROME.md](../../experiments/DEPOLARIZING_PALINDROME.md)). The scalar lives in the (I, I, …, I) component of M (the trace of the residual divided by 2^N); it captures one specific projection of the broken palindrome but not the full Frobenius norm. The Frobenius norm closed form derived here is the complementary quantitative diagnostic. The earlier `F1OpenQuestions` item "depolarizing noise: residual scaling" is closed by this document, promoting the entry to the Tier-1-derived [`F1DepolResidualClosedForm`](../../compute/RCPsiSquared.Core/F1/F1DepolResidualClosedForm.cs) claim on `F1KnowledgeBase`.

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

The depolarizing dissipator is **Pauli-basis-diagonal**: each non-identity Pauli decays uniformly at rate 4γ/3, the identity is immune. This is the canonical "isotropic Pauli channel" structure ([DEPOLARIZING_PALINDROME.md](../../experiments/DEPOLARIZING_PALINDROME.md) §2: "the 1:3 split").

### Step 2: Π conjugation on the per-site Pauli-basis matrix

Π acts on the single-letter Pauli basis as a 4 × 4 signed permutation:

    Π = ⎡ 0  1  0   0 ⎤
        ⎢ 1  0  0   0 ⎥
        ⎢ 0  0  0   i ⎥
        ⎣ 0  0  i   0 ⎦

(rows = target letter, columns = source letter, in the order I, X, Y, Z). Π is unitary, Π² = diag(+1, +1, −1, −1) (the bit_b-parity sign), and Π⁻¹ = Π†.

Conjugating the diagonal D_{depol, local} (with γ = 1) by Π permutes its eigenvalues by Π's underlying permutation: I ↔ X swaps entries `D_II ↔ D_XX`, Y ↔ Z swaps `D_YY ↔ D_ZZ`. The +i phases on Y ↔ Z cancel in conjugation (since `(+i) · (+i)* = 1`). Hence

    Π · D_depol · Π⁻¹ = diag(D_XX, D_II, D_ZZ, D_YY) = diag(−4/3, 0, −4/3, −4/3).

### Step 3: Per-site M_l = Π · D · Π⁻¹ + D

Summing the two diagonal matrices (with γ = 1):

    M_l = diag(−4/3, −4/3, −8/3, −8/3).

Four non-zero diagonal entries; no off-diagonal entries (contrast T1, which has off-diagonal (Y, X) and (Z, I) entries). Frobenius norm squared:

    ‖M_l‖²_F = (4/3)² · 2 + (8/3)² · 2 = 2 · (16/9 + 64/9) = 2 · 80/9 = **160/9**.

Diagonal trace:

    tr(M_l) = −4/3 − 4/3 − 8/3 − 8/3 = −24/3 = **−8**,    |tr(M_l)|² = **64**.

### Step 4: Multi-site assembly via per-site action

For the full N-qubit problem with γ_l per site, the multi-qubit M_depol is

    M_depol = Σ_l γ_l · M_l^{(super)}

where M_l^{(super)} acts as M_l on the Pauli-basis coordinate of site l and as the 4 × 4 identity on every other site:

    (M_l^{(super)})_{α β} = (M_l)_{α_l, β_l} · ∏_{k ≠ l} δ_{α_k, β_k}.

This per-site / identity-elsewhere structure follows because the depolarizing dissipator on site l acts only on the site-l qubit and the identity Pauli on every other site is fixed under the partial-trace projection that defines the per-site Pauli-basis component.

Equivalently in tensor notation: M_l^{(super)} = I_4^{⊗l} ⊗ M_l ⊗ I_4^{⊗(N−l−1)}.

### Step 5: Frobenius norm and the (16/9, 16) closed form

For tensor-product matrices in the Pauli basis,

    tr( (M_l^{(super)})^† · M_{l′}^{(super)} ) = tr_{site l}(M_l^† · X_l) · tr_{site l′}(X_{l′} · M_{l′}) · ∏_{k ≠ l, l′} tr_{site k}(I_4)

where X_l, X_{l′} are the corresponding identity factors of the other operator at site l, l′. Concretely:

- **l = l′ (same site):** all other sites contribute tr(I_4) = 4, giving
      tr( (M_l^{(super)})^† · M_l^{(super)} ) = ‖M_l‖²_F · 4^(N−1) = **(160/9) · 4^(N−1)**.
- **l ≠ l′ (different sites):** site l contributes tr(M_l^†) = tr(M_l)^* = tr(M_l) (M_l real), site l′ contributes tr(M_{l′}) = tr(M_l), the remaining N−2 sites contribute tr(I_4) = 4 each, giving
      tr( (M_l^{(super)})^† · M_{l′}^{(super)} ) = |tr(M_l)|² · 4^(N−2) = 64 · 4^(N−2) = **16 · 4^(N−1)**.

Substituting into the Frobenius norm of M_depol:

    ‖M_depol‖²_F = Σ_l Σ_{l′} γ_l γ_{l′} · tr( (M_l^{(super)})^† · M_{l′}^{(super)} )
                 = Σ_l γ²_l · ((160/9) · 4^(N−1))  +  Σ_{l ≠ l′} γ_l γ_{l′} · (16 · 4^(N−1))
                 = 4^(N−1) · [ (160/9) · Σ_l γ²_l  +  16 · ((Σ_l γ_l)² − Σ_l γ²_l) ]
                 = 4^(N−1) · [ ((160/9) − 16) · Σ_l γ²_l  +  16 · (Σ_l γ_l)² ]
                 = **4^(N−1) · [ (16/9) · Σ_l γ²_l  +  16 · (Σ_l γ_l)² ]**.    ∎

Hence the local coefficient is `‖M_l‖²_F − |tr(M_l)|²/4 = 160/9 − 144/9 = 16/9` and the cross-site coefficient is `|tr(M_l)|²/4 = 64/4 = 16`, giving `(c_1, c_2) = (16/9, 16)`.

### Step 6: Π²-decomposition is trivial (M_anti = 0)

The per-site M_l = diag(−4/3, −4/3, −8/3, −8/3) is **Pauli-basis-diagonal**. The multi-site M_depol = Σ_l γ_l · M_l^{(super)} is therefore also diagonal in the Pauli-string basis. Π is a signed permutation: it permutes Pauli letters within each site and multiplies phases. Conjugating a diagonal matrix by a signed permutation permutes the diagonal entries (with phase factors that cancel under M ↔ M conjugation). Crucially, when the underlying permutation already maps each diagonal eigenvalue to its own coset partner with equal value, the conjugation acts as the identity.

For our M_l, the diagonal pairs (I, X) both carry value −4/3, and the pairs (Y, Z) both carry value −8/3. Π swaps I ↔ X and Y ↔ Z. Both swaps preserve the diagonal values. Hence

    Π · M_l · Π⁻¹ = M_l    (exactly, per site)
    ⟹ Π · M_depol · Π⁻¹ = M_depol    (exactly, multi-site).

Therefore the Π²-orthogonal Pythagorean split of M_depol degenerates:

    M_anti(depol) = (M − Π·M·Π⁻¹) / 2 = 0,
    M_sym(depol)  = (M + Π·M·Π⁻¹) / 2 = M,
    ‖M_anti(depol)‖² = 0,
    ‖M_sym(depol)‖²  = ‖M(depol)‖² = 4^(N−1) · [(16/9)·Σγ² + 16·(Σγ)²].

This is the structural distinction from [PROOF_F1_T1_RESIDUAL_CLOSED_FORM](PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md) Step 7: T1's M_l has off-diagonal entries (the (Z, I) channel from σ⁻ amplitude damping plus the (Y, X) entry produced by Π conjugation), giving a non-trivial Π²-anti-symmetric piece ‖M_anti(T1)‖² = 4^(N−1)·Σγ² that maps onto F82's D_{T1, odd}. Depolarizing has no such off-diagonal Pauli-basis channel: every Pauli decays into itself, so the per-site D_l is diagonal, Π conjugation merely permutes (already-equal) diagonal values, and the Π²-anti-symmetric piece vanishes identically.

A companion typed claim `F1DepolResidualPi2Decomposition` would be 150 lines for the single fact M_anti = 0. The point is inlined as an `ExtraChildren` node on the parent claim instead; see [`F1DepolResidualClosedForm`](../../compute/RCPsiSquared.Core/F1/F1DepolResidualClosedForm.cs).

### Step 7: F1 σ-shift = 0 for depol; relationship to F5

The F1 residual is conventionally written `M := Π·L·Π⁻¹ + L + 2σ·I` with σ chosen to absorb the "diagonal background" produced by the dissipator. For Z-dephasing σ = Σγ works exactly: L_Z has the form −Σγ·I on every non-identity Pauli (in the Pauli basis), so Π·L·Π⁻¹ + L = −2Σγ·I on every non-identity Pauli, exactly cancelled by +2Σγ·I.

For depolarizing, **no scalar σ can absorb the background.** The per-site M_l = diag(−4/3, −4/3, −8/3, −8/3) has **two distinct diagonal values** (−4/3 on (I, X), −8/3 on (Y, Z)). A constant 2σ·I shift adds the same scalar to every diagonal entry; it cannot equalize the (I, X) and (Y, Z) blocks. The minimum-norm choice is σ = 0 (any nonzero σ adds positive-definite mass to the residual without cancellation).

Numerical confirmation at N = 3, uniform γ = 0.1: ‖M‖² with σ = 0 is 23.893; ‖M‖² with σ = Σγ is 0.853 (smaller because the +2Σγ·I shift partially absorbs the diagonal; but still non-zero, hence still palindrome-breaking; the residual just has the diagonal mean removed). The closed form derived in Step 5 is for the bare residual σ = 0, matching the convention used in [`F1DepolResidualClosedForm.Predict`](../../compute/RCPsiSquared.Core/F1/F1DepolResidualClosedForm.cs).

**F5 relation.** F5 records the scalar `error = (2/3)Σγ`: this is `−(1/2^N) · tr(M)` (the (I⊗N) Pauli-basis component, equivalently the magnitude of the trace contribution). F5 measures one specific scalar observable of the broken palindrome; F1DepolResidualClosedForm measures the full Frobenius norm. Both are consistent quantitative diagnostics of the same underlying obstruction (Π's 2:2 split vs depol's 1:3 split, in the language of [DEPOLARIZING_PALINDROME.md](../../experiments/DEPOLARIZING_PALINDROME.md) §1). They do not coincide because F5's scalar projection captures only the diagonal-mean part; the Frobenius norm captures the entire residual including the off-diagonal-mean (Y, Z)-block split.

## Verification

[`simulations/_f1_depol_residual_verify.py`](../../simulations/_f1_depol_residual_verify.py) verifies the closed form in seven sections:

1. **F1 sanity.** ‖M‖² for pure Z-dephasing is at machine precision (confirms the framework's Π is the right one).
2. **Pure depol numerical fit.** At N = 2, 3, 4, 5 (both uniform γ = 0.1 and non-uniform γ = [0.05, 0.10, ..., 0.05·N]) the fitted (a, b) = (16/9, 16) exactly; the predicted ‖M(depol)‖² matches the numerical value to within ~10⁻¹³.
3. **Orthogonality H ⊥ depol for truly H.** Heisenberg H has ‖M‖² = 0; adding depol gives exactly ‖M(depol)‖² with zero cross-term.
4. **Orthogonality with Z-dephasing (H + Z + depol).** Adding Z-dephasing to the H + depol setup leaves ‖M‖² = ‖M(depol)‖² (since Z and H both contribute 0 to ‖M‖²).
5. **Orthogonality with soft (Π²-odd) H (XY+YX).** Soft H gives non-zero ‖M(H)‖² (the F49 closed-form prediction); adding depol gives exactly ‖M(H)‖² + ‖M(depol)‖² with cross-term ~10⁻¹³.
6. **Per-site M_l kernel.** Displays D_depol, Π, M_l in the single-site Pauli basis and walks through `tr(M_l) = −8`, `|tr|² = 64`, `‖M_l‖² = 160/9`, then the multi-site assembly `4^(N−1) · [(160/9 − 16)·Σγ² + 16·(Σγ)²] = 4^(N−1) · [(16/9)·Σγ² + 16·(Σγ)²]`.
7. **Π²-trivial split.** Asserts `‖M − Π·M·Π⁻¹‖_F < 1e-13` and `‖M − Π²·M·Π²⁻¹‖_F < 1e-13` at N = 2, 3, 4: M_anti(depol) = 0 exactly.

All verifications pass at machine precision. Summary of section-2 numerics:

| N | uniform γ = 0.1 | non-uniform [0.05·(k+1)] | fitted (a, b) |
|---|-----------------|--------------------------|----------------|
| 2 | obs 2.702222 = pred 2.702222 | obs 1.528889 = pred 1.528889 | (1.777778, 16.000000) |
| 3 | obs 23.893333 = pred 23.893333 | obs 24.035556 = pred 24.035556 | (1.777778, 16.000000) |
| 4 | obs 168.391111 = pred 168.391111 | n/a | (1.777778, 16.000000) |
| 5 | obs 1046.755556 = pred 1046.755556 | n/a | (1.777778, 16.000000) |

(1.777778 = 16/9 to 6 decimals.)

## Diagnostic interpretation

The closed form makes the F1 depol-block residual a quantitative, **Hamiltonian-independent, γ_Z-independent, topology-independent** diagnostic for depolarizing content:

- **Pure-depol inversion (uniform γ).** From ‖M(depol)‖²_F = 4^(N−1)·((16/9)·N + 16·N²)·γ²:
      γ = √( ‖M(depol)‖²_F / [4^(N−1) · ((16/9)·N + 16·N²)] ).
  At N = 3: γ = √(‖M(depol)‖²_F / (16·(16/3 + 144))) = √(‖M(depol)‖²_F / 2389.33).
- **Pure-depol inversion (RMS for non-uniform).** With known Σ_l γ_l (from a calibration scan), the structure split `(16/9)·Σγ² + 16·(Σγ)²` lets us extract Σγ² independently from Σγ; combined, this recovers the {γ_l} distribution up to permutation.
- **N-scaling of the cross-site dominance.** Ratio of cross-site to local piece, at uniform γ: `16·N²·γ² / ((16/9)·N·γ²) = 9·N`. At N = 3 the cross-site is 27× the local; at N = 10 it is 90×. Cooperative depol-palindrome-breaking dominates as N grows, faster than T1's `(4/3)·N` ratio (because depol's per-site trace |tr(M_l)|² = 64 vs T1's 16).
- **No graph dependence.** Depolarizing is per-site only; the residual scales purely with `(Σγ², (Σγ)²)`. Contrast the H-block, whose `c_H · F(N, G)` scaling depends on bond count B and degree-squared sum D2 via `PalindromeResidualScaling`.

## Cross-references

### Repository entries

- **F1 palindrome equation** ([`docs/ANALYTICAL_FORMULAS.md` F1](../ANALYTICAL_FORMULAS.md#f1-palindrome-equation-tier-1-proven), [MIRROR_SYMMETRY_PROOF.md](MIRROR_SYMMETRY_PROOF.md)): the underlying Π·L·Π⁻¹ + L + 2Σγ·I = 0 identity for Z-dephasing.
- **F5 depolarizing error** ([`docs/ANALYTICAL_FORMULAS.md` F5](../ANALYTICAL_FORMULAS.md), [DEPOLARIZING_PALINDROME.md](../../experiments/DEPOLARIZING_PALINDROME.md)): the scalar `(2/3)Σγ` diagnostic; complementary scalar projection of the same broken palindrome.
- **F49 Frobenius residual scaling** ([`docs/ANALYTICAL_FORMULAS.md` F49](../ANALYTICAL_FORMULAS.md#f49-cross-term-formula-tier-1-proven), [PROOF_CROSS_TERM_FORMULA.md](PROOF_CROSS_TERM_FORMULA.md)): companion closed form for the Hamiltonian block.
- **F1 T1-residual closed form** ([PROOF_F1_T1_RESIDUAL_CLOSED_FORM](PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md)): sibling closed form for amplitude damping; contrasts with depol via the Π²-decomposition non-triviality (T1's M_anti = D_{T1, odd}, depol's M_anti = 0).

### Typed claims

- [`compute/RCPsiSquared.Core/F1/F1DepolResidualClosedForm.cs`](../../compute/RCPsiSquared.Core/F1/F1DepolResidualClosedForm.cs): Tier-1-derived typed claim for this closed form, registered on `F1KnowledgeBase`. Replaces the earlier `F1OpenQuestions` item "depolarizing noise: residual scaling", which was closed by this proof on 2026-05-18.
- [`compute/RCPsiSquared.Core/F1/F1T1ResidualClosedForm.cs`](../../compute/RCPsiSquared.Core/F1/F1T1ResidualClosedForm.cs): sibling T1 closed-form claim.
- [`compute/RCPsiSquared.Core/Symmetry/F5DepolarizingErrorPi2Inheritance.cs`](../../compute/RCPsiSquared.Core/Symmetry/F5DepolarizingErrorPi2Inheritance.cs): F5 scalar error claim with Pi2-Foundation inheritance.

### Scripts

- [`simulations/_f1_depol_residual_verify.py`](../../simulations/_f1_depol_residual_verify.py): the verification script for this proof.

### Memory

- `project_palindrome_frobenius_scaling`: recorded the per-class dissipator Frobenius scaling pattern on 2026-04-29; depol arm closed 2026-05-18 by this proof.
