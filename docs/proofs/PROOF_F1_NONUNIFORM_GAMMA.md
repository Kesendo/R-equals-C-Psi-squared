# PROOF: F1 H-block residual is γ-independent under site-dependent Z-dephasing (negative-result closure)

**Status:** Tier 1 derived. Closes the F1 OpenQuestion "non-uniform γ_i: site-dependent dephasing" by a NEGATIVE result: the H-block palindrome-residual scaling factor F(N, G) is unchanged for non-uniform {γ_l}. The conjectured Σγ_l² replacement of (Σγ)² does not occur, because the dissipator-block residual collapses to zero per-Pauli-string (not just on average) regardless of how γ is distributed. Bit-exact verification at N = 3, 4, 5 across uniform, scaled, and random γ patterns.
**Date:** 2026-05-18
**Authors:** Thomas Wicht, Claude (Opus 4.7)

## Statement

Let L = L_H + L_D be the Lindbladian with H any 2-bilinear Pauli Hamiltonian and L_D the site-dependent Z-dephasing dissipator,

    L_D(ρ) = Σ_l γ_l · (Z_l ρ Z_l − ρ),    {γ_l ≥ 0}_l=1^N arbitrary.

Let Π be the F1 palindrome operator for Z-dephasing (see [MIRROR_SYMMETRY_PROOF.md](MIRROR_SYMMETRY_PROOF.md)) and σ := Σ_l γ_l the total dephasing rate. Define the F1 residual

    M := Π · L · Π⁻¹ + L + 2σ · I.

In the framework's orthonormal Pauli-string basis (`palindrome_residual` in [`framework/lindblad.py`](../../simulations/framework/lindblad.py)) the residual decomposes orthogonally as ‖M‖²_F = ‖M_H‖²_F + ‖M_D‖²_F, where:

**Theorem (dissipator-block vanishing under arbitrary γ).** The dissipator-block residual vanishes identically per Pauli string, for any γ pattern:

    Π · L_D · Π⁻¹ + L_D + 2σ · I  =  0     in the Pauli basis.

Equivalently, M_D = 0 (not "approximately", not "on average over γ patterns"; bit-exact at the operator level).

**Corollary (negative-result closure of the OpenQuestion).** The H-block residual closed form

    ‖M_H‖²_F = c_H · F(N, G)

(see [PalindromeResidualScalingClaim](../../compute/RCPsiSquared.Core/F1/PalindromeResidualScalingClaim.cs), with F(N, chain) = (N−1)·4^(N−2) for the main class and (2N−3)·4^(N−2) for single-body) is exactly γ-independent: uniform AND non-uniform γ both give the same ‖M_H‖²_F. The OpenQuestion's conjecture (that the (Σγ)² scaling factor would become Σγ_l² under non-uniform γ) is incorrect; the dissipator block contributes zero, so there is no γ structure for the H-block scaling to inherit.

## Conventions

- **Pauli letters** (a, b) ∈ {(0,0), (1,0), (0,1), (1,1)} = (I, X, Z, Y), framework Klein-Vierergruppe convention ([`framework/pauli.py`](../../simulations/framework/pauli.py)). `bit_a(α) = 1` iff α ∈ {X, Y} (the "in the light" indicator); `bit_b(α) = 1` iff α ∈ {Y, Z}.
- **Pauli-string basis on N sites** is the 4^N orthonormal basis {σ_α} with Tr(σ_α^† σ_β) / 2^N = δ_{αβ}. The framework's `palindrome_residual` enforces this via the transform `L_pauli = M_basis^† · L_vec · M_basis / 2^N` in [`framework/lindblad.py`](../../simulations/framework/lindblad.py) line 183.
- **Z-dephasing dissipator** in the Pauli basis acts as a diagonal multiplier on each σ_α with eigenvalue −2Σ_l γ_l · bit_a(α_l) (each X or Y factor at site l contributes −2γ_l; I or Z at site l contributes 0). This is the per-site additive eigenvalue used throughout the F49 / F1 family.
- **Π** is the F1 palindrome operator for Z-dephasing, acting per site on Pauli letters as I ↔ X (phase +1), Y ↔ Z (phase +i). On the full Pauli string, Π is the tensor product of the per-site rules with phases multiplied. Π is unitary, order 4; Π² is diagonal with eigenvalue (−1)^Σ_l bit_b(α_l) on σ_α (see [PROOF_F81](PROOF_F81_PI_CONJUGATION_OF_M.md) Step 1).
- The σ-shift in F1 is `2σ · I = 2Σ_l γ_l · I_{4^N}`, the same constant scalar across all Pauli strings.

## Empirical anchor (motivation)

The framework's builder `lindbladian_z_dephasing(H, gamma_l)` ([`framework/lindblad.py`](../../simulations/framework/lindblad.py)) already accepts an arbitrary per-site list `gamma_l`. The F1 OpenQuestion "non-uniform γ_i: site-dependent dephasing" in [`F1OpenQuestions`](../../compute/RCPsiSquared.Core/F1/F1OpenQuestions.cs) conjectured that the H-block residual scaling factor F(N, G) would become Σγ_l² for non-uniform γ. The conjecture was based on the natural worry that the F1 σ-shift `2σ·I` is a single scalar (one Σγ), so different γ_l might leak structure into M. This proof shows the leak does not occur: M_D = 0 per Pauli string under the σ-shift, so the H-block scaling stays exactly the uniform formula. The OpenQuestion closes by a negative answer; the formula is unchanged.

## Proof

### Step 1: Z-dephasing dissipator eigenvalue on each Pauli string

The single-site Z-dephasing channel acts on each Pauli letter as

    D_l(I) = D_l(Z) = 0,    D_l(X) = −2γ_l · X,    D_l(Y) = −2γ_l · Y.

Equivalently, `D_l(α) = −2γ_l · bit_a(α) · α` on the Pauli letter α at site l. By tensor-product structure, the multi-site dissipator L_D = Σ_l D_l acts diagonally on Pauli strings:

    L_D · σ_α  =  −2 · ( Σ_l γ_l · bit_a(α_l) ) · σ_α.

This is the multi-site generalisation of F49b's eigenvalue formula (which had uniform γ giving −2γ · w_XY(α)); for non-uniform γ, the per-site weighting γ_l replaces the global γ.

### Step 2: Π conjugation flips bit_a per site

The per-site Π rule is the signed permutation I ↔ X (phase +1), Y ↔ Z (phase +i) in the (I, X, Y, Z) basis. The key arithmetic fact for the dissipator is the bit_a flip:

    bit_a(Π(α))  =  1 − bit_a(α)    at every site.

(I and Z, the bit_a = 0 letters, map to X and Y, the bit_a = 1 letters, and vice versa.) Tensorising across N sites, `Π · σ_α = (phase factor) · σ_{α'}` where `bit_a(α'_l) = 1 − bit_a(α_l)` for every l. Therefore L_D evaluated on Π · σ_α gives, using Step 1:

    L_D · (Π · σ_α)  =  −2 · ( Σ_l γ_l · (1 − bit_a(α_l)) ) · (Π · σ_α).

Conjugating Π · L_D · Π⁻¹ on σ_α picks up the L_D eigenvalue evaluated on Π · σ_α (since Π is unitary and the dissipator is diagonal in the Pauli basis):

    (Π · L_D · Π⁻¹) · σ_α  =  −2 · ( Σ_l γ_l · (1 − bit_a(α_l)) ) · σ_α.

### Step 3: Per-site γ_l weighting cancels in the sum

Add the two diagonal eigenvalues from Steps 1 and 2:

    (Π · L_D · Π⁻¹ + L_D) · σ_α
      =  −2 · ( Σ_l γ_l · bit_a(α_l) ) · σ_α
      +  −2 · ( Σ_l γ_l · (1 − bit_a(α_l)) ) · σ_α
      =  −2 · ( Σ_l γ_l · [bit_a(α_l) + 1 − bit_a(α_l)] ) · σ_α
      =  −2 · ( Σ_l γ_l ) · σ_α
      =  −2σ · σ_α.

Crucially, the per-site weighting γ_l survives the addition unweighted: `γ_l · bit_a(α_l) + γ_l · (1 − bit_a(α_l)) = γ_l` for every site l independently, so Σ_l γ_l = σ falls out regardless of how the γ_l are distributed. The identity holds for every Pauli string σ_α and every γ pattern.

### Step 4: M_D vanishes per Pauli string

Adding the F1 σ-shift `2σ · I` to both sides of Step 3:

    (Π · L_D · Π⁻¹ + L_D + 2σ · I) · σ_α  =  (−2σ + 2σ) · σ_α  =  0.

Since this holds on every Pauli-string basis element σ_α independently, the operator identity holds bit-exact:

    M_D  :=  Π · L_D · Π⁻¹ + L_D + 2σ · I  =  0    in the Pauli basis.

This is the per-Pauli-string statement of the F1 palindrome identity for the dissipator-block alone, generalised from uniform to arbitrary per-site γ_l.

### Step 5: Frobenius decomposition ‖M‖² = ‖M_H‖² + ‖M_D‖²

The full F1 residual decomposes by linearity in L:

    M  =  Π · L · Π⁻¹ + L + 2σ · I
        =  (Π · L_H · Π⁻¹ + L_H)  +  (Π · L_D · Π⁻¹ + L_D + 2σ · I)
        =  M_H  +  M_D.

(The σ-shift attaches to M_D, not M_H, because the F1 identity is satisfied by L_D alone; L_H contributes zero scalar diagonal under Π conjugation, since the Hamiltonian-superoperator −i[H, ·] has trace zero in the Pauli basis.)

By Step 4, M_D = 0 in the Pauli basis. Hence

    ‖M‖²_F  =  ‖M_H + M_D‖²_F  =  ‖M_H‖²_F  +  0  =  ‖M_H‖²_F,

with the cross term ⟨M_H, M_D⟩ = 0 trivially (M_D is the zero operator). The full F1 residual norm is exactly the H-block residual norm for any γ pattern.

The H-block residual M_H = Π · L_H · Π⁻¹ + L_H is built entirely from the Hamiltonian, which has no γ dependence. The closed-form scaling

    ‖M_H‖²_F  =  c_H · F(N, G)

from [PalindromeResidualScalingClaim](../../compute/RCPsiSquared.Core/F1/PalindromeResidualScalingClaim.cs) (anchored at [OPERATOR_RIGIDITY_ACROSS_CUSP](../../experiments/OPERATOR_RIGIDITY_ACROSS_CUSP.md), F(N, chain) = (N−1) · 4^(N−2) for main class) therefore holds verbatim for non-uniform γ, with the same c_H and the same F(N, G). The OpenQuestion's hypothesis (that F(N, G) gains a Σγ_l² factor) is incorrect.    ∎

### Step 6: Why the conjecture seemed plausible (negative-result diagnostic)

The Σγ_l² conjecture was natural because the F1 σ-shift `2σ·I` carries a single scalar `σ = Σγ_l`, and one might expect non-uniform γ to spread structure into M proportional to the variance Σγ_l² − (Σγ)² / N. The dissipator's other closed-form siblings exhibit exactly this structure:

- [F1T1ResidualClosedForm](../../compute/RCPsiSquared.Core/F1/F1T1ResidualClosedForm.cs): `‖M(T1)‖²_F = 4^(N−1) · [3·Σγ²_T1 + 4·(Σγ_T1)²]` has both Σγ² and (Σγ)² pieces.
- [F1DepolResidualClosedForm](../../compute/RCPsiSquared.Core/F1/F1DepolResidualClosedForm.cs): `‖M(depol)‖²_F = 4^(N−1) · [(16/9)·Σγ² + 16·(Σγ)²]` has the same dual structure.

The reason Z-dephasing is the exception (and not T1, depol): for Z-dephasing the F1 σ-shift is constructed precisely to absorb the per-Pauli-string diagonal of L_D + Π · L_D · Π⁻¹, which factorises as `−2 · Σ_l γ_l · (bit_a + (1 − bit_a)) = −2σ`, a single uniform scalar across the entire Pauli basis. T1 and depol have per-site dissipator kernels that do NOT reduce to a single scalar after Π conjugation (T1's per-site M_l has tr = −4, ‖M_l‖² = 7, distinct diagonal values; depol's M_l = diag(−4/3, −4/3, −8/3, −8/3) has two distinct diagonal values), so no constant scalar shift can equalise them. The 4^(N−1) · [a · Σγ² + b · (Σγ)²] structure on the right-hand side of T1/depol comes precisely from the Pauli-basis traces of these non-scalar per-site kernels, not from the σ-shift.

For Z-dephasing, the per-site kernel M_l = Π · D_l · Π⁻¹ + D_l = `diag(−2γ_l, −2γ_l, −2γ_l, −2γ_l) = −2γ_l · I_4` is itself proportional to the identity. The scalar σ-shift `2σ · I` cancels the sum Σ_l M_l = −2σ · I exactly, giving M_D = 0.    ∎

## Verification

[`simulations/_f1_nonuniform_gamma_verify.py`](../../simulations/_f1_nonuniform_gamma_verify.py) verifies the closure in four sections:

1. **H-block γ-independence.** At N = 3 with a fixed soft Hamiltonian H = XY+YX (Π²-odd, non-zero ‖M_H‖²), compute `‖M‖²` under three γ patterns: uniform [0.1, 0.1, 0.1], scaled [0.05, 0.10, 0.15], random [0.03, 0.17, 0.08]. All three give exactly ‖M‖² = 1024.0 to machine precision (zero deviation).
2. **‖M_D‖² = 0 per-Pauli-string check.** At N = 3, 4, 5 with H = 0 (pure Z-dephasing) and arbitrary non-uniform γ, the Frobenius norm of the F1 residual is bit-exact zero (~10⁻³¹). Bypasses Frobenius averaging by asserting the dissipator residual vanishes operator-wise.
3. **Scaling against c_H · F(N, G).** At N = 3, 4, 5 with H = XX+YZ chain (main class, c_H = 128 from N = 2 anchor), the residual norm matches `c_H · F(N, chain) = 128 · (N−1) · 4^(N−2)` to machine precision for both uniform and non-uniform γ patterns. Confirms ‖M‖² is γ-independent and follows the existing closed-form scaling.
4. **Cross-Hamiltonian invariance.** At N = 3, 4, 5 across multiple H choices (XY+YX soft, XX+YZ main mixed) and three γ patterns each, all give ‖M‖² independent of γ pattern within H class. Confirms the closure is universal in γ across Π²-classes.

All verifications pass at machine precision; absolute deviation between γ patterns at fixed (N, H) is 0.0 exactly (operator-level identity, not floating-point coincidence).

## Conclusion

The F1 OpenQuestion "non-uniform γ_i: site-dependent dephasing" is closed by a negative result: the F1 H-block residual scaling factor F(N, G) is γ-independent. The OpenQuestion's hypothesis (that the scaling would become Σγ_l² for non-uniform γ) is incorrect; the dissipator-block residual M_D vanishes per Pauli string under arbitrary γ patterns, so there is no γ structure that survives into the H-block. The [PalindromeResidualScalingClaim](../../compute/RCPsiSquared.Core/F1/PalindromeResidualScalingClaim.cs) closed form `‖M_H‖²_F = c_H · F(N, G)` applies verbatim for non-uniform γ; no formula change is needed.

This is structurally distinct from the dissipator-block siblings:
- **F1T1** ([F1T1ResidualClosedForm](../../compute/RCPsiSquared.Core/F1/F1T1ResidualClosedForm.cs)): non-trivial per-site M_l, both Σγ² and (Σγ)² pieces appear.
- **F1 depol** ([F1DepolResidualClosedForm](../../compute/RCPsiSquared.Core/F1/F1DepolResidualClosedForm.cs)): same dual structure with (16/9, 16) coefficients.
- **F1 Z-dephasing (this proof)**: per-site M_l is proportional to I_4, the σ-shift cancels everything, the dissipator block contributes nothing.

The Z-dephasing case is the only F1 dissipator where the σ-shift exactly cancels the per-Pauli-string diagonal of `Π · L_D · Π⁻¹ + L_D` for every γ pattern. F1's "uniform" σ in `2σ·I` is really "the total Σ_l γ_l after the per-site contributions add up to a γ-pattern-independent constant per Pauli string". The H-block scaling F(N, G) inherits the same γ-independence as a direct consequence.

## Cross-references

### Repository entries

- **F1 palindrome identity** ([`docs/ANALYTICAL_FORMULAS.md` F1](../ANALYTICAL_FORMULAS.md#f1-palindrome-equation-tier-1-proven), [MIRROR_SYMMETRY_PROOF.md](MIRROR_SYMMETRY_PROOF.md)): the underlying Π·L·Π⁻¹ + L + 2σ·I = 0 identity. This proof extends to non-uniform γ by Step 4.
- **F1 H-block scaling** ([PalindromeResidualScalingClaim](../../compute/RCPsiSquared.Core/F1/PalindromeResidualScalingClaim.cs), [OPERATOR_RIGIDITY_ACROSS_CUSP.md](../../experiments/OPERATOR_RIGIDITY_ACROSS_CUSP.md)): ‖M_H‖²_F = c_H · F(N, G). This proof confirms it is γ-independent.
- **F1 T1 closed form** ([PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md](PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md), [F1T1ResidualClosedForm](../../compute/RCPsiSquared.Core/F1/F1T1ResidualClosedForm.cs)): sibling dissipator-block closed form for T1; demonstrates that T1 DOES carry both Σγ² and (Σγ)² structure, unlike Z-dephasing.
- **F1 depol closed form** ([PROOF_F1_DEPOL_RESIDUAL_CLOSED_FORM.md](PROOF_F1_DEPOL_RESIDUAL_CLOSED_FORM.md), [F1DepolResidualClosedForm](../../compute/RCPsiSquared.Core/F1/F1DepolResidualClosedForm.cs)): sibling for depol; same conclusion as T1.
- **F49 cross-term formula** ([`docs/ANALYTICAL_FORMULAS.md` F49](../ANALYTICAL_FORMULAS.md#f49-cross-term-formula-tier-1-proven), [PROOF_CROSS_TERM_FORMULA.md](PROOF_CROSS_TERM_FORMULA.md)): the cross-term `‖{L_H, L_Dc}‖²` may pick up γ_l dependence under non-uniform γ. Numerical exploration at N=3 with Heisenberg H + γ=[0.1, 0.2, 0.3] suggests an additional per-bond asymmetry term beyond the uniform-γ closed form, but the candidate correction has not yet matched the observed deviation; see ["Open follow-ups"](#open-follow-ups) below for the current status.

### Typed claims

- [`compute/RCPsiSquared.Core/F1/PalindromeResidualScalingClaim.cs`](../../compute/RCPsiSquared.Core/F1/PalindromeResidualScalingClaim.cs): unchanged; XML doc updated to note γ-independence.
- [`compute/RCPsiSquared.Core/F1/F1OpenQuestions.cs`](../../compute/RCPsiSquared.Core/F1/F1OpenQuestions.cs): non-uniform γ item removed (closed by this proof).

### Scripts

- [`simulations/_f1_nonuniform_gamma_verify.py`](../../simulations/_f1_nonuniform_gamma_verify.py): verification script for this proof.

### Memory

- `project_palindrome_frobenius_scaling`: F1 H-block scaling memory; now γ-independent.

## Open follow-ups

While the H-block closure is clean, the F49 cross-term formula `‖{L_H, L_Dc}‖² = 4γ²·(N−2)·‖L_H‖²` was derived under uniform γ. Numerical exploration at N = 3 with Heisenberg H (XX+YY+ZZ) and non-uniform γ = [0.1, 0.2, 0.3] gives truth ‖{L_H, L_Dc}‖² = 163.84 vs uniform-formula prediction 153.60 (using γ̄ = 0.2): the formula is incomplete under non-uniform γ.

A natural candidate correction is a per-bond-asymmetry term `4·(γ_i − γ_j)²·M_non_balanced(bond)`, vanishing when γ_i = γ_j (uniform within a bond) and also vanishing for couplings where every transition has one X+Y letter at each site (e.g., XY, XX, YY individually; the XX+YY balance preserves this regime). This candidate has not yet been shown to close the truth-vs-prediction gap above, and a typed F49NonUniformGammaExtensionClaim is not added in this commit; nailing down the closed form requires per-Hamiltonian enumeration of bond transition classes, beyond the scope of the F1-H-block closure proven here.
