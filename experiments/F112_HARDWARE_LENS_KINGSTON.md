# F112 Hardware Lens on Kingston: 4 Datasets, 11 Pair-Runs, 55 Fitted Liouvillians

**Status:** F112 (`LindbladBitBPiBalance`) diagnostic applied to four existing IBM Kingston Tier-A datasets, no new QPU spend. The first reading (2026-05-26 block_cpsi single-pair) was extended same day to three more datasets, surfacing the first structural counterexample to the probes 1-14 broader empirical envelope (not to F112's typed Tier1Derived theorem, which remains intact 33/33 in-scope fits).
**Date:** 2026-05-26
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Data:**
- `data/ibm_block_cpsi_saturation_may2026/...20260508T032749Z.json` (Kingston q13-q14, 5 t-points, idle)
- `data/ibm_cusp_slowing_april2026/...20260416_212042.json` (Kingston A_mid q124-q125 + B_high, 2 pair-runs × 6 t-points, idle)
- `data/ibm_chain_gamma0_april2026/...20260419_110200.json` (Kingston Q12-Q19 chain, 4 pair-runs × 9 t-points, XY+YX bond Hamiltonian via Trotter, J = 1.838 rad/μs)
- `data/ibm_f95_angle_steering_may2026/...omega0.130_...json` and `...omega0.250_...json` (Kingston q82-q83 + q13-q14, 2 omega × 2 pairs × 6 t-points, with applied RZ Z-drive omega·Σ Z_l)
**Scripts:**
- [`simulations/_f112_block_cpsi_analysis.py`](../simulations/_f112_block_cpsi_analysis.py) (original single-dataset block_cpsi)
- [`simulations/_f112_hardware_lens_multi.py`](../simulations/_f112_hardware_lens_multi.py) (multi-dataset extension, Welle 2)
**Lens:** F112 polarity-asymmetry diagnostic on fitted effective Liouvillians

## Setup

The chosen dataset is the only one in the repo (per the F112-data inventory of 2026-05-26) that combines two qubit pairs of full ρ tomography with a multi-time-point delay sweep, making it the only Tier-A reconstructable trajectory for which F112's `polarity_coordinates_from_L` diagnostic can be evaluated on a fitted effective L without new QPU spend.

The protocol:
- Initial state: (|D_0⟩ + |D_1⟩) / √2 = (|00⟩ + (|01⟩+|10⟩)/√2) / √2 on Kingston qubits 13, 14
- Pure-decoherence idle (no applied Hamiltonian during the delay)
- 5 t-points: 0, 120, 240, 360, 480 μs
- 16 Pauli expectations per t-point (full 2-qubit tomography)
- Calibrated T2_min = 480 μs, γ_eff = 1/T2 ≈ 0.00208 per μs
- Documented anomaly in the dataset README: hardware C_block decays ~1.72× faster than pure-T2 predicts

The F112 question on this anomaly: does the hardware-effective noise that explains the trajectory sit inside F112's typed Tier1Derived scope (Hermitian H + bit_b-homogeneous c)?

## Method

For each of five candidate noise models, fit the model's parameters to the 5-point ρ(t) trajectory via Frobenius² least-squares (Nelder-Mead), then compute the F112 polarity asymmetry on the fitted L:

| Model | Parameters | F112 scope |
|---|---|---|
| `pure_Z` | γ_Z per qubit | in scope (Z is bit_b-homogeneous) |
| `Z + T1` | + σ⁻ amplitude damping per qubit | outside scope (σ⁻ = (X − iY)/2 bit_b-mixed) |
| `Z + ZZ` | + ZZ-crosstalk Hamiltonian | in scope (ZZ is bit_b-homogeneous, Hermitian) |
| `Z + T1 + ZZ` | combined | outside scope (T1 component) |
| `Z + h_y` | + single-site Y transverse field | in scope (Y is bit_b-homogeneous, Hermitian) |

For each fitted L, compute `polarity_coordinates_from_L(L_pauli, N=2, σ=Σγ)` and read `asymmetry`. F112 predicts asymmetry = 0 bit-exact for any model in its typed scope.

## Result

```
Model                       fit RMS   in F112 scope       F112 asym       F112 rel asym
----------------------------------------------------------------------------------------
pure_Z                       0.4500       YES         +0.000000e+00       0.0000e+00
Z_plus_T1                    0.4996        no         +0.000000e+00       0.0000e+00
Z_plus_ZZ                    0.4058       YES         +0.000000e+00       0.0000e+00
Z_plus_T1_plus_ZZ            0.4582        no         +0.000000e+00       0.0000e+00
Z_plus_hy                    0.3317       YES         +0.000000e+00       0.0000e+00
```

Two threads, one structural, one fit-quality:

**Structural (F112 asymmetry).** Every fitted model gives asymmetry = 0 bit-exact, including the two outside-scope cases (`Z + T1` and `Z + T1 + ZZ` with σ⁻ as bit_b-mixed collapse operator). This reproduces the empirical observation from probes 1–14: the standard Lindblad construction channel always produces asymmetry = 0, even for c that fall outside F112's typed Tier1Derived scope. The σ⁻ T1 case sits in F112's Tier1Candidate non-Hermitian-extension envelope; the empirical bit-exactness here is one more instance consistent with the open Im⟨L_{H_re,-i}, L_{H_im,-i}⟩ = 0 identity holding.

**Fit-quality (which model explains the data).** The best-fit model is `Z + h_y` (RMS = 0.332), a 26% improvement over the pure-Z baseline (RMS = 0.450). The fitted transverse field is h_y ≈ -0.00144 per μs, comparable to the T2 dephasing rate γ_z ≈ 0.0021 per μs. The `Z + ZZ` model is intermediate (RMS = 0.406, 10% improvement), with a small fitted ZZ coupling. The two T1-containing models do not improve on pure-Z (Nelder-Mead pushes the σ⁻ rates to ≈ 0); adding T1 does not help fit the data.

The interpretation: the dominant beyond-pure-T2 channel on Kingston q13–q14 during this idle protocol is consistent with a small effective transverse Y drift (a low-frequency Larmor-like detuning), not amplitude damping. This is a Hermitian, bit_b-homogeneous H addition, which F112 classifies as fully inside its Tier1Derived scope.

## What F112 added to the existing reading

The block-CΨ saturation dataset's own README documents the 1.72× hardware-over-T2 decay-rate anomaly. Without F112, that anomaly is a number without a structural classifier. F112 adds two pieces of information:

1. **A bit_b-axis classifier on the dominant extra channel.** Of the candidate noise models that materially improve the fit, the best one (Z + h_y) is bit_b-homogeneous on both H and c sides, locating the anomaly inside F112's typed scope. The two amplitude-damping models that would have placed the anomaly outside the typed scope are decisively rejected by the fit (T1 contributes negligibly).

2. **An empirical confirmation that the standard Lindblad construction preserves polarity balance across all five candidate models, including the two outside the typed scope.** The 240-config empirical envelope from probes 7–12 and the F87-orthogonality probe at N=3 are extended by 5 more datapoints from a different system (N=2 instead of N=3, different bit-arrangement of Bell-like Dicke initial state, hardware-derived rather than analytically constructed L). Every standard Lindblad channel observed empirically balanced; the non-Hermitian Tier1Candidate envelope of F112 (open Im⟨L_{H_re,-i}, L_{H_im,-i}⟩ = 0 identity) holds on this hardware-derived L too.

## Limits of this reading

- **N=2, single qubit pair.** F112's typed scope is N-universal but this reading is at the smallest non-trivial N. Other Tier-A datasets in the repo (cusp_slowing April26, chain_gamma0 April19, f95_angle_steering May16) provide additional N=2 trajectories that could be analyzed identically, multiplying the empirical anchor count.
- **5 t-points constrain only ~2-4 effective parameters tightly.** The model selection ranks the 5 candidates by RMS but does not give confidence intervals. A higher-density time grid (e.g., 9-15 points) would tighten the h_y estimate and make T1 contribution falsifiable rather than just dispreferred.
- **F112's diagnostic is structural, not predictive.** The asymmetry = 0 reading does not say "the hardware noise is Tier1Derived F112"; it says "every standard Lindblad model fitting the data preserves the polarity balance, including ones outside the typed scope." Distinguishing typed-scope models from out-of-scope models requires reading the model's c structure, not the asymmetry itself.
- **Fitted dephasing rates exceed T2 calibration substantially.** Even the best-fit Z + h_y model has γ_Z ≈ 0.015 per μs on qubit 14, ~7× the T2 calibration value. This is the "1.72× anomaly" multiplied by additional structural mismatch between the model family and the true hardware noise. Further model classes (mixed dephase letters, non-Lindblad terms) would need to be tested to close this gap.

## Welle 2 extension: three more Kingston datasets (2026-05-26)

The block_cpsi reading above was the first hardware application of F112 and operated on a single dataset (5 t-points, 1 qubit pair, idle protocol). Same-day extension covered the three remaining Tier-A datasets identified by the F112-data inventory:

| Dataset | Pairs | t-points / pair | Applied H during evolution | bit_b structure of H |
|---|---|---|---|---|
| block_cpsi (May 8) | q13-q14 (1) | 5 | none (idle) | trivial |
| cusp_slowing (Apr 16) | q124-q125 (A_mid), B_high (2) | 6 | none (idle) | trivial |
| chain_gamma0 (Apr 19) | Q12-Q13, Q13-Q14, Q14-Q15, Q15-Q19 (4) | 9 | J · (XX + YY)/2 bond bilinear (Trotter) | bit_b = 0 (Π²-Z-even) |
| f95 angle_steering (May 16) | A_mid q82-q83 + B_high q13-q14 × 2 omega values (4) | 6 | omega · Σ Z_l single-site Z drive | bit_b = 1 (Z-axis) |

Total: 11 pair-runs × 5 candidate noise models = 55 fitted Liouvillians.

### Aggregate result

```
F112 typed Tier1Derived scope (Hermitian H + bit_b-homog c):
  in-scope BALANCED:  33  (F112 prediction met, all bit-exact)
  in-scope BROKEN:    0   (no counterexamples to the typed theorem)

Out-of-scope (bit_b-mixed c via σ⁻ T1):
  out-scope BALANCED: 15  (broader empirical envelope holds)
  out-scope BROKEN:    7  (broader envelope counterexamples)

Max observed relative asymmetry: 4.24e-03 (f95 B_high q13-q14 omega=0.13, Z_plus_T1_plus_ZZ fit)
```

### The 7 counterexamples are structurally specific

Every BROKEN fit is in the f95 angle_steering dataset, where the protocol applies a single-site Z-drive omega · Σ Z_l during evolution. With idle protocols (block_cpsi, cusp_slowing) and the XY-bond chain Hamiltonian (chain_gamma0), the σ⁻ T1 fits gave asymmetry = 0 bit-exact regardless of T1 magnitude. With the Z-drive H, σ⁻ T1 gives non-zero asymmetry.

Synthetic isolation confirms it is structural, not optimization artifact (verified independently inline in `simulations/_f112_hardware_lens_multi.py`'s sanity loop):

```
Case A: Z-drive H + Z-deph             rel asym = 0.000        (in F112 scope: BALANCED)
Case B: idle H + σ⁻ T1                 rel asym = 0.000        (out of scope, still BALANCED)
Case C: Z-drive H + σ⁻ T1              rel asym = 3.85e-03     (out of scope, BROKEN bit-exact)
Case D: idle H + Z-deph + σ⁻ T1        rel asym = 0.000        (out of scope, still BALANCED)
Case E: Z-drive H + Z-deph + σ⁻ T1     rel asym = 3.85e-03     (out of scope, BROKEN bit-exact)
```

The asymmetry magnitude is identical in cases C and E, confirming the Z-deph term is irrelevant to the breaking (it stays in F112-scope and contributes nothing to the +i / −i Π-eigenvalue content). The breaker is the cross-bit_b interference between the Z-drive H (bit_b=1) and the σ⁻ collapse operator (bit_b ∈ {0, 1}).

### What this sharpens

The probes 1-14 reflection (`reflections/POLARITY_COORDINATES.md`) noted that bit_b-mixed σ⁻ T1 was empirically observed to preserve balance across every probe tested. F112's typed theorem covers only bit_b-homogeneous c; the broader observation was that bit_b-mixed c "happens to" also balance in every tested case. This was the empirical envelope behind the Tier1Candidate non-Hermitian extension scope (although technically the non-Hermitian Candidate concerns non-Hermitian H, not bit_b-mixed c).

The f95 hardware reading is the first structural counterexample to the bit_b-mixed-c-preserves-balance reading. The sharpened statement: bit_b-mixed c preserves balance **only when the Hermitian H is trivial or bit_b-homogeneous through a bilinear (XX, YY, etc.) channel that does not produce cross-bit_b matrix elements jointly with σ⁻**. A single-site bit_b-homogeneous H (omega·Z) DOES produce such cross-bit_b elements via [Z, σ⁻] = -2σ⁻ × bit_b shifts.

This does not threaten F112's typed Tier1Derived theorem. It does refine the broader empirical envelope and contains useful information for future predictions: **any hardware protocol that combines a Z-drive (or other bit_b-homogeneous single-site H) with measurable T1 amplitude damping will produce a non-zero F112 polarity asymmetry, of order γ_T1 · omega_Z / Σ_rates²**. The 4.24e-3 magnitude is consistent with this scale (omega = 0.13 / μs, γ_T1 ≈ 0.0005 / μs from the fit, Σ ≈ 0.13 / μs → expected scale 0.13 · 0.0005 / 0.017 ≈ 4e-3).

### Best-fit model distribution

Across the 11 pair-runs:
- Z + h_y (transverse Y drift): 7 pair-runs best fit
- Z + T1 + ZZ (combined channels): 3 pair-runs (all f95)
- Z + ZZ: 1 pair-run

The 7-of-11 dominance of `Z + h_y` is consistent across idle datasets and is the central within-F112-scope finding: most Kingston pairs are well-modeled by Z-dephasing plus a small effective transverse Y drift, both bit_b-homogeneous. The 3 f95 pair-runs that prefer T1+ZZ over the alternatives sit in a regime where the structural counterexample applies; F112's reading sharpens those into honest "out-of-scope" classifications rather than mis-attributing them to the typed-scope envelope.

## Connections

- **F112** (`docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md`, `docs/ANALYTICAL_FORMULAS.md` F112 entry, `compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiBalance.cs`): the typed Tier1Derived theorem and Claim. This experiment is the first hardware application.
- **Block-CΨ saturation lens** (`experiments/THEOREM_1_AND_2_BLOCK_CPSI_SATURATION.md`, Confirmations entry `block_cpsi_saturation_kingston_may2026`): the original lens through which this dataset was analyzed. F112 adds an independent structural reading.
- **F87 trichotomy** (orthogonal axis on shared bit_b grading): F87 was not applicable here (the protocol has no bilinear Hamiltonian during delay), so this reading is F112-only. A future dataset with non-trivial H during evolution could give the joint (F87 class, F112 asymmetry) reading.
- **Probes 1–14** (`reflections/POLARITY_COORDINATES.md` Closure section): the discovery arc that led to F112. This hardware reading extends the empirical envelope to a hardware-derived L.
- **`feedback_qpu_conservative`** (memory): F112's first hardware reading was produced from existing 2026-05-08 data, no new QPU spend, per the conservative-QPU discipline.

## Reproduction

```
python -X utf8 simulations/_f112_block_cpsi_analysis.py
```

Runs in under 5 seconds on a desktop. No QPU access required.
