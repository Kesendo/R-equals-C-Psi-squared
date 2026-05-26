# F113: Closed-Form for the F112 Counterexample Asymmetry

**Status:** Tier 1 derived for N=2, 3, 4 via constructive parameter sweep. Closed-form formula for the polarity-asymmetry break magnitude when F112's typed scope is violated by the canonical Z-drive × amplitude-damping interference.
**Date:** 2026-05-26
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Script:** [`simulations/_f113_break_formula_derivation.py`](../simulations/_f113_break_formula_derivation.py)
**Builds on:** F112 ([PROOF_F112](../docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md), [LindbladBitBPiBalance](../compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiBalance.cs)) + the structural counterexample discovered in Welle 2 ([F112_HARDWARE_LENS_KINGSTON.md](F112_HARDWARE_LENS_KINGSTON.md))

## Theorem (F113)

For a Lindblad-form Liouvillian L = -i[H, ·] + Σ_k γ_k · D[c_k] with:

- Hermitian H that may include single-site Z-drives of the form Σ_l (ω_l/2)·Z_l
- Dissipator c_k that may include σ⁻_l (amplitude damping) at rate γ_T1,l per site and σ⁺_l (pumping) at rate γ_pump,l per site
- Any other bit_b-homogeneous additions (Z-dephasing, ZZ/XX/YY/XY bilinear bonds, single-site X-drives, single-site Y-drives), which contribute 0 individually by F112

the F112 polarity-coordinate asymmetry has the closed form:

    asymmetry := ‖M_plus_half‖² − ‖M_minus_half‖²
              = (4^N / 2) · Σ_l ω_l · (γ_T1,l − γ_pump,l)

Bit-exactly. Verified at N=2, 3, 4 across multiple parameter samples (relative deviation < 1e-12, dominated by floating-point precision).

## Empirical anchor (constructive verification)

`simulations/_f113_break_formula_derivation.py` runs three classes of test:

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

| N | predicted coefficient (N/2)·4^N | measured (5 random samples) | max deviation |
|---|---|---|---|
| 2 | 16.0 | 16.000000 each | 7.46e-14 |
| 3 | 96.0 | 96.000000 each | 3.07e-12 |
| 4 | 512.0 | 512.000000 each | 3.18e-12 |

**4. Per-site decomposition (single-site Z-drive on q_l + σ⁻ on q_l only):**

| N | site l | observed asym / (ω · γ_T1) | predicted (1/2)·4^N |
|---|---|---|---|
| 3 | 0 | 32.0 | 32.0 |
| 3 | 1 | 32.0 | 32.0 |
| 3 | 2 | 32.0 | 32.0 |

**5. Cross-site (Z-drive on q_a, σ⁻ on q_b, a ≠ b):** asym = 0.0 bit-exact. The break is local: only same-site (Z-drive_l, σ⁻_l) pairs contribute.

**6. Non-uniform rates:** for ω_l = (0.05, 0.1, 0.2), γ_T1,l = (0.001, 0.002, 0.003) at N=3:
- Σ_l 0.5·4^N · ω_l · γ_T1,l = 0.027200 (formula prediction)
- Measured asymmetry = 0.027200 (ratio 1.000000)

## Scope (what doesn't contribute to F113)

The following H/c structures give F112 asymmetry = 0 bit-exact and therefore do NOT contribute to the F113 break magnitude:

- Single-site X-drive Σ_l h_x,l·X_l
- Single-site Y-drive Σ_l h_y,l·Y_l
- Bilinear bond Hamiltonians (ZZ, XX, YY, XY, YX, XZ, ZX, YZ, ZY) on any bond
- Z-dephasing dissipator γ_Z·Σ_l D[Z_l] (in F112's typed scope)
- σ⁻ + σ⁺ at equal rate (T1 cooling + pumping in detailed balance: γ_T1 = γ_pump → contributions cancel)

The formula is therefore additive in only two channels: single-site Z-drives crossed with same-site σ⁻ minus same-site σ⁺.

## Sign convention

- Sign(asymmetry) = Sign(ω_l) for each l (verified by reversing ω_l: asymmetry flips sign exact).
- σ⁺ contributes opposite sign to σ⁻ (asymmetry = (4^N/2) · Σ_l ω_l · γ_T1,l with σ⁻ alone; same magnitude with sign flipped if c = σ⁺ instead).

## Why X-drive and Y-drive give zero

For [Pauli_letter, σ⁻] commutators:
- [Z, σ⁻] = -2·σ⁻ (proportional to σ⁻ itself, non-Hermitian)
- [X, σ⁻] = Z (Hermitian)
- [Y, σ⁻] = i·Z (imaginary Hermitian)

The F112 break requires non-Hermitian Π-eigenspace coupling between H and c. Only Z-drive produces this structure: its commutator with σ⁻ is proportional to the non-Hermitian σ⁻ itself, which carries Π-eigenvalue +i / −i imbalance. X and Y drives produce Hermitian commutators that are F112-symmetric.

This is the structural origin of F113's restriction to Z-drives.

## Implications

- **Hardware fingerprinting.** Any hardware protocol that combines a single-site Z-drive (deliberate or as a Stark shift) with amplitude damping will exhibit measurable F112 asymmetry. The asymmetry magnitude directly extracts γ_T1 · ω / (predictable structural factor). At the f95 angle-steering parameters (ω = 0.13, γ_T1 ≈ 0.001 per μs, N=2 effective), predicted asymmetry is 16 · 0.13 · 0.001 = 2.08e-3, matching the Welle 2 hardware-fit value bit-exact.
- **Calibration tool.** Inverted: given measured F112 asymmetry on hardware-effective L, F113 directly gives γ_T1,l · ω_l from the measurement. Could become a per-site T1-extraction protocol when the drive parameters are known.
- **F112 typed-scope sharpening.** F112's typed Tier1Derived covers Hermitian H + bit_b-homogeneous c, giving asymmetry = 0. F112's empirical envelope was loosely "bit_b-mixed c also balances", refuted by Welle 2. F113 provides the exact closed-form for the regime where the envelope breaks; together F112 + F113 give a complete picture of the polarity-axis behavior across the standard Lindblad family.

## Universal-N status

Verified bit-exact at N=2, 3, 4. The structural origin (commutator algebra [Z, σ⁻] = -2σ⁻ + the 4^N · (1/2) per-site factor) is N-universal in form: the 4^N prefactor matches the operator-space dimension 4^N = d²; the per-site additivity matches the locality of [Z_l, σ⁻_m] = -2σ⁻_m · δ_{lm}. A rigorous N-universal proof would derive the (1/2) · 4^N coefficient algebraically from the Π-eigenspace structure.

Status: **Tier 1 derived at N = 2, 3, 4** via constructive parameter sweep + per-site / cross-site / sign / variant verification. **Tier 1 candidate for general N** (form is universal-N-suggestive but algebraic proof open).

## Reproduction

```
python -X utf8 simulations/_f113_break_formula_derivation.py
```

Runs in under 30 seconds; produces the multivariate fit + N-scaling + per-site / cross-site / variant tests inline.

## Related

- [F112 PROOF](../docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md): parent theorem, F113 lives in F112's "broken empirical envelope" regime.
- [F112_HARDWARE_LENS_KINGSTON.md](F112_HARDWARE_LENS_KINGSTON.md): hardware discovery of the counterexample; F113 now gives its closed-form magnitude.
- [LindbladBitBPiBalance.cs](../compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiBalance.cs): F112's typed Claim; F113 informs the empirical-envelope-counterexample story documented in its inspectables.
- [LindbladBitBPiBalanceWitness.cs](../compute/RCPsiSquared.Diagnostics/Polarity/LindbladBitBPiBalanceWitness.cs): the StandardSet's 5th witness (`Zdrive_with_T1_envelope_BROKEN`) is a specific instance of F113 at ω=0.13, γ_T1=0.001, N=2.
- F112 (Tier1Derived Hermitian H + bit_b-homog c → asymmetry = 0) and F113 (closed-form for the bit_b-mixed-c break magnitude) together cover the standard-Lindblad-family polarity behavior completely.
