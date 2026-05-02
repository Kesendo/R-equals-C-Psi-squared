# F-Formula Crosswalk: `RCPsiSquared.Diagnostics` ↔ `docs/ANALYTICAL_FORMULAS.md` + `docs/proofs/`

Each Diagnostics module is a **verkettung** (composition) of `RCPsiSquared.Core` primitives
into a specific F-anchored reading. This map points each C# class to its F-formula in
[`docs/ANALYTICAL_FORMULAS.md`](../../docs/ANALYTICAL_FORMULAS.md), its proof document in
[`docs/proofs/`](../../docs/proofs/), and any hardware-confirmation entry in
[`Core.Confirmations.ConfirmationsRegistry`](../RCPsiSquared.Core/Confirmations/ConfirmationsRegistry.cs).

For Core-side F-formulas (Pauli, Lindblad, Symmetry, etc.) see
[`RCPsiSquared.Core/F_FORMULA_CROSSWALK.md`](../RCPsiSquared.Core/F_FORMULA_CROSSWALK.md).

---

## F49 / F85 — Frobenius scaling

| C# | F-formula | Proof | Confirmation |
|----|-----------|-------|--------------|
| `F49.FrobeniusScaling.PredictNormSquared` | F49 closed-form ‖M‖²_F = c_H · F(N, G) | OPERATOR_RIGIDITY_ACROSS_CUSP.md | — |
| `F49.FrobeniusScaling.PredictNormSquaredFromTerms` | F85 k-body generalisation via Π²-class | [`PROOF_F85_KBODY_GENERALIZATION.md`](../../docs/proofs/PROOF_F85_KBODY_GENERALIZATION.md) | — |

## F77 — Pauli-pair trichotomy

| C# | F-formula | Proof | Confirmation |
|----|-----------|-------|--------------|
| `F77.PauliPairTrichotomy.Classify` | F77 truly/soft/hard classifier | docs/ANALYTICAL_FORMULAS.md F77 entry | `palindrome_trichotomy` (Marrakesh 2026-04-26) |

## F80 — Bloch sign-walk

| C# | F-formula | Proof | Confirmation |
|----|-----------|-------|--------------|
| `F80.BlochSignWalk.PredictMSpectrumImaginaryParts` | Spec(M) = ±2i · Spec(H_non-truly) for chain Π²-odd 2-body | [`PROOF_F80_BLOCH_SIGNWALK.md`](../../docs/proofs/PROOF_F80_BLOCH_SIGNWALK.md) | — |

## F81 — Π-decomposition (Tier 1)

| C# | F-formula | Proof | Confirmation |
|----|-----------|-------|--------------|
| `F81.PiDecomposition.Decompose` | Tier-1: M_anti = L_{H_odd}, ‖M‖² = ‖M_sym‖² + ‖M_anti‖² | [`PROOF_F81_PI_CONJUGATION_OF_M.md`](../../docs/proofs/PROOF_F81_PI_CONJUGATION_OF_M.md) | `f83_pi2_class_signature_marrakesh` (April 30) |

## F82 — T1 dissipator correction

| C# | F-formula | Proof | Confirmation |
|----|-----------|-------|--------------|
| `F82.T1DissipatorPrediction.PredictViolation` | ‖D_{T1, odd}‖_F = √(Σγ²)·2^(N−1) | [`PROOF_F82_T1_DISSIPATOR_CORRECTION.md`](../../docs/proofs/PROOF_F82_T1_DISSIPATOR_CORRECTION.md) | `palindrome_trichotomy` T1-amplification refutation |
| `F82.T1DissipatorPrediction.EstimateRmsT1FromViolation` | inverse: γ_T1, RMS = f81_violation / (√N · 2^(N−1)) | same | — |

## F83 — Π-decomposition prediction (anti-fraction)

| C# | F-formula | Proof | Confirmation |
|----|-----------|-------|--------------|
| `F83.PiDecompositionPrediction.Predict` | F81/F85 ‖M‖²/‖M_anti‖²/‖M_sym‖² closed forms | [`PROOF_F83_PI_DECOMPOSITION_RATIO.md`](../../docs/proofs/PROOF_F83_PI_DECOMPOSITION_RATIO.md) + [`PROOF_F85_KBODY_GENERALIZATION.md`](../../docs/proofs/PROOF_F85_KBODY_GENERALIZATION.md) | `f83_pi2_class_signature_marrakesh` |
| `F83.PiDecompositionPrediction.AntiFraction` | anti-fraction = 1/(2 + 4·r) | same | — |

## F84 — Amplitude damping (thermal)

| C# | F-formula | Proof | Confirmation |
|----|-----------|-------|--------------|
| `F84.AmplitudeDampingPrediction.PredictViolation` | ‖D_{AmplDamp, odd}‖_F = √(Σ(γ_↓ − γ_↑)²)·2^(N−1) | [`PROOF_F84_AMPLITUDE_DAMPING.md`](../../docs/proofs/PROOF_F84_AMPLITUDE_DAMPING.md) | — |
| `F84.AmplitudeDampingPrediction.EstimateRmsNetCoolingFromViolation` | inverse: |Δγ|_RMS from violation | same | — |

## DZero — d=0 substrate (kernel of L)

| C# | F-formula | Proof | Confirmation |
|----|-----------|-------|--------------|
| `DZero.StationaryModes.Compute` | F4 sector projectors as kernel modes (N+1 dim) | docs/ANALYTICAL_FORMULAS.md F4 | — |
| `DZero.DZeroDecomposition.Decompose` | ρ = ρ_d0 + ρ_d2 via kernel projector; F4 + THE_POLARITY_LAYER.md | F4 derivation | — |
| `DZero.SectorPopulations.FromDensityMatrix` | p_n = Tr(P_n · ρ); Hamming-weight binning | F4 sector projectors | `d_zero_sector_trichotomy_marrakesh` (May 1) |

## Polarity — +0/0/−0 axis reading

| C# | F-formula | Proof | Confirmation |
|----|-----------|-------|--------------|
| `Polarity.PolarityDiagnostic.FromDensityMatrix` | per-site Bloch ⟨X⟩,⟨Y⟩,⟨Z⟩ — no specific F-formula | [`docs/THE_POLARITY_LAYER.md`](../../docs/THE_POLARITY_LAYER.md), [`docs/PRIMORDIAL_QUBIT.md`](../../docs/PRIMORDIAL_QUBIT.md) §9 | — |

## Ptf — Perspectival Time Field

| C# | Source | Hypothesis | Confirmation |
|----|--------|------------|--------------|
| `Ptf.PerturbationMatrixElements.Compute` | first-order PT ⟨W_s|V_L|M_{s'}⟩ matrix elements | [`hypotheses/PERSPECTIVAL_TIME_FIELD.md`](../../hypotheses/PERSPECTIVAL_TIME_FIELD.md) (Tier 2 after EQ-014) | — |
| `Ptf.PerturbationMatrixElements.EigenvectorShift` | δM_s = Σ_{s'≠s} [⟨W_{s'}|V_L|M_s⟩ / (λ_s − λ_{s'})] · M_{s'} | same | — |

---

## How to read a Diagnostics file without guessing

1. Open the C# file (e.g. `F77/PauliPairTrichotomy.cs`).
2. The class XML doc explicitly names the F-formula (e.g. "F77") and links to:
   - `docs/ANALYTICAL_FORMULAS.md` (the formula entry)
   - The proof in `docs/proofs/` if one exists
   - The Core primitive(s) it composes (`PauliHamiltonian`, `PalindromeResidual`, etc.)
   - A `ConfirmationsRegistry` entry name if hardware-confirmed
3. This crosswalk is the consolidated index — every Diagnostics class is one row.
4. Tests under `RCPsiSquared.Diagnostics.Tests/` carry the same F-anchor in xUnit `[Fact]` names
   (e.g. `XXplusYY_IsTruly_UnderZDephasing`).

## When to add a new entry

If you create a new diagnostic that:
- Implements an F-formula → add a row above with the F-entry + proof pointer.
- Composes Core machinery without a new F-formula → add a row noting the source docs and
  the Core primitives it chains.
- Is hardware-confirmed → cross-reference the `ConfirmationsRegistry` entry name.
