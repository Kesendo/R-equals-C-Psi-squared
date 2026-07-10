# PROOF: F1 palindrome residual closed form under T1 amplitude damping

**Status:** Tier 1 derived. Closed-form Frobenius norm of the F1 residual M for the T1 amplitude-damping block, derived analytically from the per-site action of Π and the T1 dissipator; bit-exact verification at N = 2, 3, 4, 5 across uniform and non-uniform γ_T1.
**Date:** 2026-05-18
**Authors:** Thomas Wicht, Claude (Opus 4.7)

## Abstract

The F1 palindrome theorem says that under pure Z-dephasing, the residual operator M (which measures how far the Lindbladian falls from its Π-conjugation image) closes to zero for any palindromic Hamiltonian. T1 amplitude damping breaks the closure: it is the canonical noise process that the Z-dephasing F1 symmetry does not see. The natural follow-up question is how the residual grows under T1, and what the closed form of its Frobenius norm looks like as a function of the T1 rates.

This proof gives the answer. The pure-T1 contribution to the residual norm splits into two structural pieces: a local part proportional to the sum of squared T1 rates, and a cooperative part proportional to the square of the sum of T1 rates. The cooperative part is the interesting one. It grows quadratically with the number of sites when the rates are uniform, signaling that T1 dissipators on different sites do not contribute independently in the residual norm; they overlap. The local part counts each site once; the cooperative part counts every pair once.

Two more observations make the formula self-contained. The pure-T1 contribution does not depend on the Hamiltonian (the T1 block is Frobenius-orthogonal to the H block in M), and it does not depend on the Z-dephasing rate (Z-dephasing absorbs into the +2σ·I shift that defines M itself). So the formula gives a clean separation of T1 noise from the rest of the system: change H, leave the T1 part alone; change Z-dephasing rate, leave the T1 part alone.

Combined with the F49 closed form for the H block (proven separately), this gives a complete Frobenius dictionary for the residual norm under standard Z + T1 channels: a closed-form sum of three pieces, each from a structurally distinct source. The diagnostic upshot is that any deviation between predicted and measured residual norm flags either a Hamiltonian outside F49's bilinear scope, a non-standard T1 channel, or the presence of additional noise (depolarizing, dephasing on other axes) not captured by Z + T1 alone.

## Statement

Let L = L_H + L_Z + L_T1 be the Lindbladian for a 2-bilinear Pauli Hamiltonian H, single-site Z-dephasing with rates {γ_l}, and single-site T1 amplitude damping with rates {γ^{T1}_l}. Let Π be the F1 palindrome operator for Z-dephasing (see [the Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md)). Define the F1 residual

    M := Π · L · Π⁻¹ + L + 2 · Σ_l γ_l · I.

In the framework's orthonormal Pauli-string basis (the basis used by `palindrome_residual` in [`framework/lindblad.py`](../../simulations/framework/lindblad.py); see also [F49](../ANALYTICAL_FORMULAS.md#f49-cross-term-formula-tier-1-proven) and [the F49 Frobenius dictionary in `project_palindrome_frobenius_scaling`](../../experiments/OPERATOR_RIGIDITY_ACROSS_CUSP.md)), the residual norm decomposes orthogonally as

    ‖M‖²_F = 2^(N+2) · n_YZ · ‖H‖²_F  +  4^(N−1) · [ 3 · Σ_l γ^{T1}²_l  +  4 · (Σ_l γ^{T1}_l)² ].

The H part is the F49 closed form ([`docs/ANALYTICAL_FORMULAS.md` F49](../ANALYTICAL_FORMULAS.md#f49-cross-term-formula-tier-1-proven), proven separately). This document closes the T1 part:

**Theorem (T1 block closed form).** The pure-T1 residual contribution satisfies

    ‖M(T1)‖²_F = 4^(N−1) · [ 3 · Σ_l γ^{T1}²_l  +  4 · (Σ_l γ^{T1}_l)² ].

Three structural facts follow immediately and are verified to machine precision:

1. **Hamiltonian-independent.** Pure-T1 contribution does not depend on H; H-block and T1-block are Frobenius-orthogonal in M.
2. **γ_Z-independent.** Z-dephasing is absorbed by the +2Σγ·I shift; pure-T1 contribution does not depend on {γ_l}.
3. **Two structural pieces.** `3·Σγ²` is the per-site (local) contribution; `4·(Σγ)²` is the cross-site (cooperative) contribution that grows quadratically in N when γ_T1 is uniform.

## Conventions

- **Pauli letters** are indexed (a, b) ∈ {(0,0), (1,0), (0,1), (1,1)} = (I, X, Z, Y) following the framework's Klein-Vierergruppe convention ([`framework/pauli.py`](../../simulations/framework/pauli.py)). `bit_a` is the X/Y indicator, `bit_b` is the Y/Z indicator.
- **Pauli-string basis on N sites** is the 4^N orthonormal basis {σ_α} with the inner product ⟨σ_α | σ_β⟩ = δ_{αβ}; equivalently `Tr(σ_α^† σ_β) / 2^N = δ_{αβ}`. Operators on operator-space (4^N × 4^N "super-operator" matrices) inherit this orthonormality via the Hilbert-Schmidt product `⟨A, B⟩ = Tr(A^† B)`. The framework's `palindrome_residual` enforces this via the transform `L_pauli = M_basis^† · L_vec · M_basis / 2^N` in [`framework/lindblad.py`](../../simulations/framework/lindblad.py) line 183.
- **Π** is the F1 palindrome operator for Z-dephasing, acting per site on Pauli letters as
      I ↔ X (phase +1),    Y ↔ Z (phase +i).
  See [the Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md) (the "Conjugation Operator Π" section). On the full Pauli string Π is the tensor product of per-site rules with phases multiplied. Π is unitary, order-4: Π² is diagonal with eigenvalue (−1)^{Σ_l bit_b(α_l)} on the string σ_α (see [F81](PROOF_F81_PI_CONJUGATION_OF_M.md) Step 1).
- **T1 dissipator** uses the lowering convention σ⁻ = (X + iY)/2 = [[0, 1], [0, 0]] (taking |1⟩ → |0⟩), σ⁺ = (σ⁻)†. The per-site channel is
      D_{T1, l}(ρ) = γ^{T1}_l · [σ⁻_l ρ σ⁺_l − ½ {σ⁺_l σ⁻_l, ρ}].
  See [`framework/lindblad.py`](../../simulations/framework/lindblad.py) `lindbladian_z_plus_t1`. Same convention as [F82](PROOF_F82_T1_DISSIPATOR_CORRECTION.md) Step 3.

## Empirical anchor (motivation, observed 2026-04-29)

Memory [`project_palindrome_frobenius_scaling`](../../) recorded the per-class dissipator scaling

    ‖M(L_T1)‖²_F = 4^(N−1) · [ c_1 · Σ_l γ²_T1_l + c_2 · (Σ_l γ_T1_l)² ]

with (c_1, c_2) = (3, 4) for σ⁻ amplitude damping, tabulated in the framework as `HARDWARE_DISSIPATORS['T1']` in [`framework/lindblad.py`](../../simulations/framework/lindblad.py). The entry was previously registered as an F1 open question ("T1 amplitude damping: full closed form"); verified bit-exact at N = 3..6, but analytical derivation was open. This document derives the (3, 4) pair from first principles, promoting the entry to the Tier-1-derived [`F1T1ResidualClosedForm`](../../compute/RCPsiSquared.Core/F1/F1T1ResidualClosedForm.cs) claim on `F1KnowledgeBase`.

## Proof

### Step 1: T1 dissipator's per-site Pauli-basis matrix

In the single-qubit orthonormal Pauli basis (I, X, Y, Z), the T1 dissipator D_{T1, local} is a 4 × 4 matrix with `D_{αβ} = (1/2) Tr(σ_α · D_{T1, local}(σ_β))`. Direct computation (also given in [F82](PROOF_F82_T1_DISSIPATOR_CORRECTION.md) Step 3) yields

|       | I        | X        | Y        | Z        |
|-------|----------|----------|----------|----------|
| **I** | 0        | 0        | 0        | 0        |
| **X** | 0        | −γ/2     | 0        | 0        |
| **Y** | 0        | 0        | −γ/2     | 0        |
| **Z** | +γ       | 0        | 0        | −γ       |

(Take γ = 1 for the closed-form derivation; the formula will be γ²-homogeneous in the per-site contribution and γ_l γ_{l′}-bilinear in the cross-site contribution.)

### Step 2: Π conjugation on the per-site Pauli-basis matrix

Π acts on the single-letter Pauli basis as a 4 × 4 signed permutation:

    Π = ⎡ 0  1  0   0 ⎤
        ⎢ 1  0  0   0 ⎥
        ⎢ 0  0  0   i ⎥
        ⎣ 0  0  i   0 ⎦

(rows = target letter, columns = source letter, in the order I, X, Y, Z). Π is unitary, Π² = diag(+1, +1, −1, −1) (the bit_b-parity sign), and Π⁻¹ = Π†.

Conjugating D_{T1, local} (with γ = 1) by Π gives

    Π · D · Π⁻¹ = ⎡ −1/2   0     0     0    ⎤
                  ⎢  0     0     0     0    ⎥
                  ⎢  0     i    −1     0    ⎥
                  ⎣  0     0     0    −1/2  ⎦

(The (Z, I) entry of D, value +1, maps to the (I, X) channel under Π but the X-column of the Y→iZ rule sends it back into the (X, ·) block; full computation in [`simulations/f1_t1_residual_verify.py`](../../simulations/f1_t1_residual_verify.py) Section 6.)

### Step 3: Per-site M_l = Π · D · Π⁻¹ + D

Summing the two matrices (with γ = 1):

    M_l = ⎡ −1/2   0      0      0    ⎤
          ⎢  0    −1/2    0      0    ⎥
          ⎢  0     i     −3/2    0    ⎥
          ⎣ +1     0      0     −3/2  ⎦

Six non-zero entries with squared magnitudes summing to

    ‖M_l‖²_F = (1/2)² + (1/2)² + 1² + (3/2)² + 1² + (3/2)² = 0.25 + 0.25 + 1 + 2.25 + 1 + 2.25 = **7**.

The diagonal trace

    tr(M_l) = −1/2 − 1/2 − 3/2 − 3/2 = **−4**,    |tr(M_l)|² = **16**.

(The +2Σγ·I shift in the definition of M is zero for the pure-T1 case since γ_l = 0, so M_l is just Π·D_l·Π⁻¹ + D_l.)

### Step 4: Multi-site assembly via per-site action

For the full N-qubit problem with γ^{T1}_l per site, the multi-qubit M_T1 is

    M_T1 = Σ_l γ^{T1}_l · M_l^{(super)}

where M_l^{(super)} acts as M_l on the Pauli-basis coordinate of site l and as the 4 × 4 identity on every other site:

    (M_l^{(super)})_{α β} = (M_l)_{α_l, β_l} · ∏_{k ≠ l} δ_{α_k, β_k}.

This per-site / identity-elsewhere structure follows because the dissipator on site l acts only on the site-l qubit and the identity Pauli on every other site is fixed under the partial-trace projection that defines the per-site Pauli-basis component.

Equivalently in tensor notation: M_l^{(super)} = I_4^{⊗l} ⊗ M_l ⊗ I_4^{⊗(N−l−1)}.

### Step 5: Frobenius norm and the (3, 4) closed form

For tensor-product matrices in the Pauli basis,

    tr( (M_l^{(super)})^† · M_{l′}^{(super)} ) = tr_{site l}( M_l^† · X_l ) · tr_{site l′}( X_{l′} · M_{l′} ) · ∏_{k ≠ l, l′} tr_{site k}(I_4)

where X_l, X_{l′} are the corresponding identity factors of the other operator at site l, l′. Concretely:

- **l = l′ (same site):** all other sites contribute tr(I_4) = 4, giving
      tr( (M_l^{(super)})^† · M_l^{(super)} ) = ‖M_l‖²_F · 4^(N−1) = **7 · 4^(N−1)**.
- **l ≠ l′ (different sites):** site l contributes tr(M_l^†) = tr(M_l)^* (since M_l is its own conjugate for our M_l matrix), site l′ contributes tr(M_{l′}) = tr(M_l), the remaining N−2 sites contribute tr(I_4) = 4 each, giving
      tr( (M_l^{(super)})^† · M_{l′}^{(super)} ) = |tr(M_l)|² · 4^(N−2) = 16 · 4^(N−2) = **4 · 4^(N−1)**.

Substituting into the Frobenius norm of M_T1:

    ‖M_T1‖²_F = Σ_l Σ_{l′} γ_l γ_{l′} · tr( (M_l^{(super)})^† · M_{l′}^{(super)} )
              = Σ_l γ²_l · (7 · 4^(N−1))  +  Σ_{l ≠ l′} γ_l γ_{l′} · (4 · 4^(N−1))
              = 4^(N−1) · [ 7 · Σ_l γ²_l  +  4 · ((Σ_l γ_l)² − Σ_l γ²_l) ]
              = 4^(N−1) · [ (7 − 4) · Σ_l γ²_l  +  4 · (Σ_l γ_l)² ]
              = **4^(N−1) · [ 3 · Σ_l γ²_l  +  4 · (Σ_l γ_l)² ]**.    ∎

Hence (c_1, c_2) = (3, 4) is derived from `(‖M_l‖²_F − |tr(M_l)|²) = 7 − 4 = 3` and `|tr(M_l)|² / 4 = 16 / 4 = 4`. The 4 is the spurious diagonal that would arise if all sites contributed identically; the 3 is the genuine single-site residual after subtracting the "uniform background" picked up by the cross-site sum.

### Step 6: Orthogonality of the T1 block to H and Z blocks

The cross-terms between L_T1 and L_H or L_Z in the residual decomposition vanish: in the Pauli basis L_T1 has its only non-zero Π-asymmetric entry at the (Z, I) site-l slot (Pauli-string transitions from `…I…` to `…Z…` at site l), while L_H supports Pauli-string transitions consistent with two-body commutators (no transitions of weight 0 → 1) and L_Z is diagonal in the Pauli basis and Π²-symmetric (see [F81](PROOF_F81_PI_CONJUGATION_OF_M.md) Step 4). The supports are disjoint, so the cross-terms in `Tr(M_T1^† · M_H)`, `Tr(M_T1^† · M_Z)` vanish identically.

This is the same orthogonality principle that gives F49 / F49b / F49c the additive cross-term structure. See [`docs/ANALYTICAL_FORMULAS.md` F49](../ANALYTICAL_FORMULAS.md#f49-cross-term-formula-tier-1-proven) and [F81](PROOF_F81_PI_CONJUGATION_OF_M.md) Step 6 for the formal Pythagorean splits.

### Step 7: Connection to F82

[F82](PROOF_F82_T1_DISSIPATOR_CORRECTION.md) isolates the **Π²-anti-symmetric** part of L_T1 (the (Z, I) site-l entries in Pauli basis) and computes

    ‖D_{T1, odd}‖²_F = 4^(N−1) · Σ_l γ²_T1_l.

The Π-anti-symmetric part of M splits the same way. For pure T1 (no H), F82's identity Π · M · Π⁻¹ = M − 2 · D_{T1, odd} gives

    M_anti = (M − Π·M·Π⁻¹) / 2 = D_{T1, odd},

so ‖M_anti‖²_F = ‖D_{T1, odd}‖²_F = 4^(N−1) · Σ_l γ²_T1_l. By Pythagoras and the (3, 4) total derived above,

    ‖M_sym‖²_F = ‖M‖²_F − ‖M_anti‖²_F
              = 4^(N−1) · [ 3·Σγ² + 4·(Σγ)² ] − 4^(N−1) · Σγ²
              = 4^(N−1) · [ 2·Σγ² + 4·(Σγ)² ].

At N = 3, uniform γ_T1 = 0.1: ‖M‖² = 7.2, ‖M_anti‖² = 0.48, ‖M_sym‖² = 6.72 (verified numerically in section 6 of the verification script alongside the F82 closed form). The (3, 4) total contains both the Π²-symmetric contributions from the (X, X), (Y, Y), (Z, Z) diagonal entries plus the (Y, X) off-diagonal of M_l (which the |tr(M_l)|² · 4^(N−2) cross-site sum picks up) AND the Π²-anti-symmetric (Z, I) entries that F82 measures in isolation. F82's diagnostic uses only the latter; the F1 residual sees both.

## Verification

[`simulations/f1_t1_residual_verify.py`](../../simulations/f1_t1_residual_verify.py) verifies the closed form in six sections:

1. **F1 sanity.** ‖M‖² for pure Z-dephasing and for Heisenberg + Z-dephasing is at machine precision (∼ 6·10⁻³¹ at N = 3): confirms the framework's Π is the right one and that ‖M‖² is genuinely orthogonal in the F1 blocks.
2. **Pure T1 numerical fit.** At N = 2, 3, 4, 5 (both uniform γ_T1 = 0.1 and non-uniform γ_T1 = [0.05, 0.10, ..., 0.05·N]) the fitted (a, b) = (3.000000, 4.000000) exactly; the predicted ‖M(T1)‖² matches the numerical value to within ~10⁻¹³ (the floating-point limit at these problem sizes).
3. **Orthogonality H ⊥ T1 for truly H.** Heisenberg H has ‖M‖² = 0; adding T1 gives exactly ‖M(T1)‖² with zero cross-term.
4. **Orthogonality with Z-dephasing.** Adding Z-dephasing to the H+T1 setup leaves ‖M‖² = ‖M(T1)‖² (since Z and H both contribute 0).
5. **Orthogonality with soft (Π²-odd) H (XY+YX).** Soft H gives non-zero ‖M(H)‖² = 1024 at N=3 (the F49 formula's prediction); adding T1 gives exactly ‖M(H)‖² + ‖M(T1)‖² = 1031.2, with cross-term ~10⁻¹³.
6. **Per-site M_l kernel.** Displays D_T1, Π, M_l in the single-site Pauli basis and walks through `tr(M_l) = −4`, `|tr|² = 16`, `‖M_l‖² = 7`, then the multi-site assembly `4^(N−1) · [(7 − 4)·Σγ² + 4·(Σγ)²]`.

All verifications pass at machine precision. Summary of section-2 numerics:

| N | uniform γ_T1 = 0.1 | non-uniform [0.05·k+0.05] | fitted (a, b) |
|---|--------------------|---------------------------|----------------|
| 2 | obs 0.880000 = pred 0.880000 | obs 0.510000 = pred 0.510000 | (3.000000, 4.000000) |
| 3 | obs 7.200000 = pred 7.200000 | obs 7.440000 = pred 7.440000 | (3.000000, 4.000000) |
| 4 | obs 48.640000 = pred 48.640000 | obs 78.400000 = pred 78.400000 | (3.000000, 4.000000) |
| 5 | obs 294.400000 = pred 294.400000 | obs 681.600000 = pred 681.600000 | (3.000000, 4.000000) |

## Diagnostic interpretation

The closed form makes the F1 T1-block residual a quantitative, **Hamiltonian-independent, γ_Z-independent** diagnostic for T1 content:

- **Pure-T1 inversion (uniform γ_T1).** From ‖M(T1)‖²_F = 4^(N−1)·(3N + 4N²)·γ²_T1:
      γ_T1 = √( ‖M(T1)‖²_F / [4^(N−1) · (3N + 4N²)] ).
  At N = 3: γ_T1 = √(‖M(T1)‖²_F / 720).
- **Pure-T1 inversion (RMS for non-uniform).** With known Σ_l γ_T1_l (e.g., from a calibration scan), the structure split `3·Σγ² + 4·(Σγ)²` lets us extract Σγ² independently from Σγ; combined, this recovers the {γ_T1_l} distribution up to permutation.
- **N-scaling of the cross-site dominance.** Ratio of cross-site to local piece, at uniform γ_T1: `4·N²·γ² / (3·N·γ²) = (4/3)·N`. At N = 3 the cross-site is `4×` the local; at N = 10 it is `13×`. Cooperative T1 palindrome-breaking dominates as N grows.

## Cross-references

### Repository entries

- **F1 palindrome equation** ([`docs/ANALYTICAL_FORMULAS.md` F1](../ANALYTICAL_FORMULAS.md#f1-palindrome-equation-tier-1-proven), [the Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md)): the underlying Π·L·Π⁻¹ + L + 2Σγ·I = 0 identity for Z-dephasing.
- **F49 Frobenius residual scaling** ([`docs/ANALYTICAL_FORMULAS.md` F49](../ANALYTICAL_FORMULAS.md#f49-cross-term-formula-tier-1-proven), [the Cross-Term Formula proof](PROOF_CROSS_TERM_FORMULA.md)): companion closed form for the Hamiltonian block.
- **F82 T1 dissipator correction** ([`docs/ANALYTICAL_FORMULAS.md` F82](../ANALYTICAL_FORMULAS.md#f82-pi-conjugation-of-m-under-t1-amplitude-damping-tier-1-proven), [the F82 T1 dissipator correction proof](PROOF_F82_T1_DISSIPATOR_CORRECTION.md)): isolates the Π²-anti-symmetric piece ‖D_{T1, odd}‖_F = √(Σγ²) · 2^(N−1); related but different quantity.
- **F84 Amplitude damping (thermal)** ([`docs/ANALYTICAL_FORMULAS.md` F84](../ANALYTICAL_FORMULAS.md#f84-pi-conjugation-of-m-under-thermal-amplitude-damping-tier-1-proven), [the F84 amplitude damping proof](PROOF_F84_AMPLITUDE_DAMPING.md)): F82's thermal generalization.

### Typed claims

- [`compute/RCPsiSquared.Core/F1/F1T1ResidualClosedForm.cs`](../../compute/RCPsiSquared.Core/F1/F1T1ResidualClosedForm.cs): Tier-1-derived typed claim for this closed form, registered on `F1KnowledgeBase`. Replaces the earlier `F1OpenQuestions` item 2 entry, which was closed by this proof on 2026-05-18.
- [`compute/RCPsiSquared.Core/F1/PalindromeResidualScalingClaim.cs`](../../compute/RCPsiSquared.Core/F1/PalindromeResidualScalingClaim.cs): companion claim covering the Hamiltonian (non-truly) part of ‖M‖²_F.
- [`framework.HARDWARE_DISSIPATORS['T1']`](../../simulations/framework/lindblad.py): already records `{'c1': 3.0, 'c2': 4.0}`; this proof derives the entries.

### Scripts

- [`simulations/f1_t1_residual_verify.py`](../../simulations/f1_t1_residual_verify.py): the verification script for this proof.

### Memory

- `project_palindrome_frobenius_scaling`: recorded the empirical (3, 4) anchor on 2026-04-29; analytical derivation thread closed 2026-05-18 by this proof.
