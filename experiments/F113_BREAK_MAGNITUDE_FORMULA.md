# F113: Closed-Form for the F112 Counterexample Asymmetry

**Status:** Tier 1 derived for N=2, 3, 4 via constructive parameter sweep. Closed-form formula for the polarity-asymmetry break magnitude when F112's typed scope is violated by the canonical Z-drive × amplitude-damping interference.
**Date:** 2026-05-26
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Script:** [`simulations/f113_break_formula_derivation.py`](../simulations/f113_break_formula_derivation.py)
**Builds on:** F112 ([the F112 proof](../docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md), [LindbladBitBPiBalance](../compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiBalance.cs)) + the structural counterexample, first noticed while fitting the Kingston f95 data ([F112 hardware lens on Kingston](F112_HARDWARE_LENS_KINGSTON.md)), whose protocol supplies the Z-drive

## Theorem (F113)

For a Lindblad-form Liouvillian L = -i[H, ·] + Σ_k γ_k · D[c_k] with:

- Hermitian H that may include single-site Z-drives of the form Σ_l (ω_l/2)·Z_l
- Dissipator c_k that may include σ⁻_l (amplitude damping) at rate γ_T1,l per site and σ⁺_l (pumping) at rate γ_pump,l per site
- Any other bit_b-homogeneous additions (Z-dephasing, ZZ/XX/YY/XY bilinear bonds, single-site X-drives, single-site Y-drives), which contribute 0 individually by F112

the F112 polarity-coordinate asymmetry has the closed form:

    asymmetry := ‖M_plus_half‖² − ‖M_minus_half‖²
              = (4^N / 2) · Σ_l ω_l · (γ_pump,l − γ_T1,l)

Verified at N=2, 3, 4 across multiple parameter samples, to a relative deviation below 1e-12, which is floating-point noise in the Frobenius norms rather than an exact zero. The sign of the asymmetry depends on a convention, spelled out under Sign convention below; the magnitude does not.

## Empirical anchor (constructive verification)

`simulations/f113_break_formula_derivation.py` runs the isolation cases, the three univariate scans and the multivariate fit below. The N-scaling, per-site, cross-site and non-uniform tables that follow were computed the same way but are not printed by it; they are reproduced by evaluating `polarity_coordinates_from_hc` at the stated parameters.

**1. Univariate scaling (each parameter independently):**

| swept parameter | observed scaling | R² |
|---|---|---|
| ω (fix γ_T1=0.001, γ_Z=0.005) | asym ∝ ω¹ | 1.000000 |
| γ_T1 (fix ω=0.13, γ_Z=0.005) | asym ∝ γ_T1¹ | 1.000000 |
| γ_Z (fix ω=0.13, γ_T1=0.001) | asym ∝ γ_Z⁰ (independent) | 1.000000 |

**2. Multivariate fit on 60 random (ω, γ_T1, γ_Z) samples at N=2:**

```
log|asym| = 2.7726 + 1.0000·log(ω) + 1.0000·log(γ_T1) + 0.0000·log(γ_Z)
```

Implied constant `exp(2.7726) = 16.000` bit-exact (std 0.000000 across 60 samples).

**3. N-scaling (per-site contribution structure):**

Reported as asym / (ω · γ_T1), so the minus is the cooling sign of γ_pump − γ_T1.

| N | predicted coefficient −(N/2)·4^N | measured (5 random samples) | max deviation |
|---|---|---|---|
| 2 | −16.0 | −16.000000 each | 7.46e-14 |
| 3 | −96.0 | −96.000000 each | 3.07e-12 |
| 4 | −512.0 | −512.000000 each | 3.18e-12 |

**4. Per-site decomposition (single-site Z-drive on q_l + σ⁻ on q_l only):**

| N | site l | observed asym / (ω · γ_T1) | predicted −(1/2)·4^N |
|---|---|---|---|
| 3 | 0 | −32.0 | −32.0 |
| 3 | 1 | −32.0 | −32.0 |
| 3 | 2 | −32.0 | −32.0 |

**5. Cross-site (Z-drive on q_a, σ⁻ on q_b, a ≠ b):** asym = 0.0 bit-exact. The break is local: only same-site (Z-drive_l, σ⁻_l) pairs contribute.

**6. Non-uniform rates:** for ω_l = (0.05, 0.1, 0.2), γ_T1,l = (0.001, 0.002, 0.003) at N=3:
- Σ_l 0.5·4^N · ω_l · (0 − γ_T1,l) = −0.027200 (formula prediction)
- Measured asymmetry = −0.027200 (ratio 1.000000)

## Scope (what doesn't contribute to F113)

The following H/c structures give F112 asymmetry = 0 bit-exact and therefore do NOT contribute to the F113 break magnitude:

- Single-site X-drive Σ_l h_x,l·X_l
- Single-site Y-drive Σ_l h_y,l·Y_l
- Bilinear bond Hamiltonians (ZZ, XX, YY, XY, YX, XZ, ZX, YZ, ZY) on any bond
- Z-dephasing dissipator γ_Z·Σ_l D[Z_l] (in F112's typed scope)
- σ⁻ + σ⁺ at equal rate (T1 cooling + pumping in detailed balance: γ_T1 = γ_pump → contributions cancel)

The formula is therefore additive in only two channels: single-site Z-drives crossed with same-site σ⁻ minus same-site σ⁺.

## Sign convention

Two separate conventions sit on this number, and both have to be named.

- **Which pairing builds the Liouvillian.** The pipeline reads a row-stack L against the order='F' Pauli transform. That pairing is mismatched, and the mismatch is exactly equivalent to conjugating Π. Either consistent pairing (column-stack L with order='F', or row-stack L with order='C') returns +2.08e-3 where this one returns −2.08e-3. Both Π and conj(Π) palindromize the Pauli-basis Liouvillian, which is exactly real, so neither is privileged and the direction is a choice, not a measurement. Everything below is stated in the pipeline's pairing.
- **Which operator σ⁻ names.** With σ⁻ = (X+iY)/2 = [[0,1],[0,0]], taking |1⟩ → |0⟩, cooling gives a negative asymmetry; (X−iY)/2 flips it.

Within one pairing the relative statements are firm: Sign(asymmetry) = −Sign(ω_l · γ_T1,l) per site (verified by reversing ω_l), and σ⁺ contributes opposite to σ⁻ at equal magnitude.

## What selects a Hamiltonian

The asymmetry reads H only through the single-site moment Tr(Z_l H). At N=2 with ω = 0.13 and γ_T1 = 0.001:

| H | Tr(Z₀H) | asymmetry |
|---|---|---|
| Z-drive | 0.26 | −2.08e-3 |
| Z-drive + XX bond | 0.26 | −2.08e-3 |
| Z-drive + ZZ bond | 0.26 | −2.08e-3 |
| Z-drive + Y-drive | 0.26 | −2.08e-3 |
| X-drive, Y-drive, ZZ, XY+YX, XX+YY alone | 0 | 0 |
| Z on site 0 only | 0.26 / 0 | −1.04e-3 |

Three different additions leave the value bit-identical. That is not a coincidence to be checked case by case: F112 kills the H-only part of the quadratic form for any Hermitian H, and a σ⁻ dissipator alone gives 0, so what survives is bilinear in (H, D) and hence LINEAR in H. Which linear functional it is, the table then shows and a wider sweep confirms: over random Hermitian H drawn on all 4^N Pauli strings at N=2 and N=3, with per-site rates all distinct, the asymmetry matches (4^N/2)·Σ_l [Tr(Z_l H)/2^(N−1)]·(γ_pump,l − γ_T1,l) to a relative 2e-12.

The mechanism behind the selector is the commutator algebra, for σ⁻ = (X+iY)/2 = [[0,1],[0,0]]:

- [Z, σ⁻] = +2·σ⁻, proportional to σ⁻ itself
- [X, σ⁻] = −Z
- [Y, σ⁻] = −i·Z

Only the first stays in the σ⁻ direction, and σ⁻ is the non-Hermitian object carrying the Π-eigenvalue +i / −i imbalance; X and Y drives leave it. Under σ⁻ = (X−iY)/2 the first two flip and [Y, σ⁻] = −i·Z is unchanged. Note that Hermiticity of the commutator is not the criterion: [Y, σ⁻] is anti-Hermitian and still gives zero.

This is the structural origin of F113's restriction to Z-drives.

## Implications

- **Hardware fingerprinting.** Any hardware protocol that combines a single-site Z-drive (deliberate or as a Stark shift) with amplitude damping will exhibit measurable F112 asymmetry. The asymmetry magnitude directly extracts ω · (γ_pump − γ_T1) / (predictable structural factor). With ω = 0.13 on both sites, γ_T1 ≈ 0.001 per μs, γ_pump = 0, N = 2, the formula gives 16 · 0.13 · (0 − 0.001) = −2.08e-3. The f95 angle-steering protocol drives ONE qubit of its pair, so at its parameters the value is −1.04e-3. No hardware anchor exists for it yet. The [Kingston survey](F112_HARDWARE_LENS_KINGSTON.md) reports nonzero asymmetry only on its Z-drive runs, but the single-site Z there is handed to the fitter as a known drive rather than fitted, so the split follows from the model set rather than from the chip; and the same data need a per-qubit detuning, itself a single-site Z, that those models do not carry. A survey of that kind could anchor a magnitude in any case, never a direction: it publishes |asymmetry| / ‖M‖², which is identical for σ⁻ and σ⁺.
- **Calibration tool, with one caveat that runs deep.** Inverted, F113 gives γ_T1,l · ω_l from the asymmetry. But the asymmetry is not an observable: it is a functional of a FITTED Liouvillian, and anyone holding a fitted L already holds γ_T1. The inversion is a consistency check on the fit, exact by construction, and it adds no calibration information the fit did not have. Could become a per-site T1-extraction protocol when the drive parameters are known.
- **F112 typed-scope sharpening.** F112's typed Tier1Derived covers Hermitian H + bit_b-homogeneous c, giving asymmetry = 0. F112's empirical envelope was loosely "bit_b-mixed c also balances", refuted by Welle 2. F113 provides the exact closed-form for the regime where the envelope breaks; together F112 + F113 give a complete picture of the polarity-axis behavior across the standard Lindblad family.

## Universal-N status

**Tier 1 derived for general N** (Welle 4, 2026-05-26): the rigorous derivation of the (1/2) · 4^N coefficient is in [`docs/proofs/PROOF_F113_COEFFICIENT_DERIVATION.md`](../docs/proofs/PROOF_F113_COEFFICIENT_DERIVATION.md). The structural decomposition is

    (1/2) · 4^N  =  4  ·  4^(N-1)  ·  (1/2)

- **factor 4**: cross-term reduction `asymmetry = 4 · Re⟨L_H,+i, L_T1,+i⟩` (from `‖L,+i‖² − ‖L,−i‖²` expansion + F112 typed + F112 non-Hermitian extension cancellations + cross-term equal-magnitude-opposite-sign relation).
- **factor 4^(N-1)**: N−1 spectator-site identity factors. Each spectator site contributes `⟨I_4, I_4⟩ = Tr(I_4) = 4` to the Frobenius inner product on tensor products. This is the operator-space dimension entering through the local Pauli dimension 4 per spectator site.
- **factor 1/2**: explicit single-site N=1 inner product `⟨(L_H,1)_{+i}, (L_T1,1)_{+i}⟩ = −ωγ/2`, derived via sympy from the explicit 4×4 matrices.

The proof has 8 steps (single-site sympy + tensor factorization of Π per site + per-site additivity of the inner product + sum over driven sites) and three lemmas (Lemma A and B from F112; Lemma C new: `L_T1` is real-valued in Pauli basis since `D[c]ρ` preserves Hermiticity for Hermitian ρ). One specific Frobenius equality in Lemma C step 5 is verified bit-exact at N = 1, 2, 3, 4, 5 but not yet algebraically closed from the support pattern; this is documented as a structural exercise and does not block the universal-N status given the bit-exact anchor across N = 1–5.

Verified numerically at N = 5 (524k+ Pauli-string pairs would be needed for full basis enumeration, but the proof reduces it to a single-site inner product × tensor factorization, so N = 5 verification is cheap).

## Reproduction

```
python -X utf8 simulations/f113_break_formula_derivation.py
```

Runs in under 30 seconds; produces the multivariate fit + N-scaling + per-site / cross-site / variant tests inline.

## Related

- [the F112 proof](../docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md): parent theorem, F113 lives in F112's "broken empirical envelope" regime.
- [F112 hardware lens on Kingston](F112_HARDWARE_LENS_KINGSTON.md): hardware discovery of the counterexample; F113 now gives its closed-form magnitude.
- [LindbladBitBPiBalance.cs](../compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiBalance.cs): F112's typed Claim; F113 informs the empirical-envelope-counterexample story documented in its inspectables.
- [LindbladBitBPiBalanceWitness.cs](../compute/RCPsiSquared.Diagnostics/Polarity/LindbladBitBPiBalanceWitness.cs): the StandardSet's 5th witness (`Zdrive_with_T1_envelope_BROKEN`) is a specific instance of F113 at ω=0.13, γ_T1=0.001, N=2.
- F112 (Tier1Derived Hermitian H + bit_b-homog c → asymmetry = 0) and F113 (closed-form for the bit_b-mixed-c break magnitude) together cover the standard-Lindblad-family polarity behavior completely.
