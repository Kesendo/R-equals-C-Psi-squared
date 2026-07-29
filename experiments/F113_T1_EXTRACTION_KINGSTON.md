# F113 T1-Extraction on Kingston f95 Angle-Steering Data

**Status:** First application of F113 closed form as a hardware diagnostic: invert the F113 formula to extract a γ_T1 reading from the polarity-asymmetry measurement and compare to device-calibrated 1/T1. Self-consistency confirmed bit-exact; the fit-vs-calibration gap bounds how much non-T1 noise the minimal model is absorbing.
**Date:** 2026-05-26
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Script:** [`simulations/f113_t1_extraction_kingston.py`](../simulations/f113_t1_extraction_kingston.py)
**Data:** [`data/ibm_f95_angle_steering_may2026/`](../data/ibm_f95_angle_steering_may2026/) (2 omega × 2 pair-runs × 6 t-points each, Kingston, 2026-05-16)
**Builds on:** F113 ([`F113_BREAK_MAGNITUDE_FORMULA.md`](F113_BREAK_MAGNITUDE_FORMULA.md), [`PROOF_F113_COEFFICIENT_DERIVATION.md`](../docs/proofs/PROOF_F113_COEFFICIENT_DERIVATION.md))

## Idea

F113 (Welle 4 closure) gives the polarity asymmetry of a Z-drive + amplitude-damping Lindblad in closed form:

    asymmetry = (4^N / 2) · Σ_l ω_l · (γ_pump,l − γ_T1,l)

The inverse: given a measured polarity asymmetry, known drive ω, and γ_pump = 0, extract a γ_T1 estimate:

    γ_T1_F113 = −asymmetry / ((N/2) · 4^N · ω)

One driven site, so the sum has a single term: at N=2 this is `γ_T1_F113 = −asymmetry / (8 · ω)`.

The question this experiment answers: when applied to real hardware data, does this F113-inversion produce a γ_T1 value consistent with (a) the fit-direct γ_T1 (self-consistency check) and (b) the device-calibrated 1/T1 (the physical-T1 reading)?

## Pipeline per pair-run

1. **Load trajectory.** 6 ρ(t) snapshots per pair-run from f95 angle-steering dataset (full 2-qubit tomography reconstructed from 4×4 rho2_real + rho2_imag).
2. **Fit minimal Z + σ⁻ T1 Lindblad model** with the applied Z-drive H = (ω/2)·Z on the STEERING QUBIT (the protocol injects RZ on one qubit of the pair, and the dataset's own predicted crossing angles match φ₀ − ω·t, not φ₀ − 2ω·t), and (γ_z, γ_T1) as free parameters.
3. **Compute the F112 polarity asymmetry** of the fitted L via `polarity_coordinates_from_hc`.
4. **Invert F113.** Apply γ_T1_F113 = −asymmetry / (8 · ω) at N=2.
5. **Compare** the three γ_T1 readings: fit-direct, F113-inverted, device-calibrated.

## Result

| Pair-run | γ_z fit | γ_T1 fit | γ_T1 F113 | F113/fit | γ_T1 calib | fit/calib | RMS |
|---|---:|---:|---:|---:|---:|---:|---:|
| ω=0.13 A_mid q82-q83 | 0.0163 | 0.00577 | 0.00577 | **1.000000** | 0.00507 | 1.14 | 0.30 |
| ω=0.13 B_high q13-q14 | 0.0123 | 0.00615 | 0.00615 | **1.000000** | 0.00464 | 1.33 | 0.36 |
| ω=0.25 A_mid q82-q83 | 0.0155 | 0.00746 | 0.00746 | **1.000000** | 0.00507 | 1.47 | 0.30 |
| ω=0.25 B_high q13-q14 | 0.0114 | 0.00564 | 0.00564 | **1.000000** | 0.00464 | 1.22 | 0.37 |

(γ values in per-μs; γ_T1 calib = mean of 1/T1 across the two qubits in the pair as reported in the dataset metadata.)

## Reading

### (1) F113 self-consistency: confirmed bit-exact

`F113/fit = 1.000000` to all decimals across all 4 pair-runs. The F113 closed form is a faithful re-derivation of the fitted Lindblad's polarity asymmetry; inversion recovers the same γ_T1 that produced the asymmetry. This is the structural self-consistency check the experiment was designed to perform, and it passes.

### (2) Fit above calibration, and the ω-trend is not readable

`fit/calib` ranges from 1.14 to 1.47. The two ω=0.25 runs do not bracket the two ω=0.13 runs (1.47 and 1.22 against 1.14 and 1.33), so the pair-to-pair spread is as large as any drive dependence, and this data cannot support a "stronger drive, more excess" reading.

And the excess is not resolved everywhere. A parametric bootstrap at the run's own shot budget (binomial on each of the 15 Pauli expectations, 2048 shots per basis, 30 replicas of the detuned fit) gives σ(γ_T1) = 0.0012 on A_mid ω=0.13, so its calibration value sits **0.54σ** from the fit. That run shows no excess at all once the shot noise is counted; only the larger ratios are worth reading, and readout error, unmitigated here at 1-3% per qubit, sits on top of that as an unquantified systematic.

The direction of the excess, where it is real, is expected and structural: the minimal Z + σ⁻ T1 model has two parameters, so any noise outside pure Z-dephasing gets absorbed by the σ⁻ channel, the only other one available.

The obvious suspect for that absorbed noise is a missing per-qubit detuning, since this chip's coherences rotate. Adding one free Z coefficient per qubit (12-start Nelder-Mead) settles it, and the answer is that γ_T1 is robust:

| pair-run | γ_T1 (Z+T1) | γ_T1 (+detuning) | γ_z | RMS | Σδ | ω + 2Σδ |
|---|---:|---:|---|---|---:|---:|
| ω=0.13 A_mid | 0.00577 | 0.00571 | 0.0163 → 0.0045 | 0.297 → 0.052 | −0.0408 | 0.048 |
| ω=0.13 B_high | 0.00615 | 0.00617 | 0.0123 → 0.0016 | 0.361 → 0.037 | −0.0298 | 0.070 |
| ω=0.25 A_mid | 0.00746 | 0.00743 | 0.0155 → 0.0044 | 0.303 → 0.066 | −0.0399 | 0.170 |
| ω=0.25 B_high | 0.00564 | 0.00564 | 0.0114 → 0.0014 | 0.366 → 0.043 | −0.0290 | 0.192 |

RMS improves by 4.6 to 9.7× and γ_z falls further, while γ_T1 moves by at most 1.1%. Σδ is consistent across all four runs at −0.029 to −0.041 per μs, i.e. the chip's own detuning opposes the applied drive. The SPLIT between the two sites is a different matter. It is not unobservable, as one might expect from the Bell coherence turning at the sum: the measured initial state carries single-excitation coherences (ρ₀₁ = 0.002 − 0.016i, ρ₀₂ = 0.016 − 0.045i) that turn at 2δ₁ and 2δ₀ separately, and least squares does resolve them, with a sharp minimum (SSE 0.0136 at δ₀−δ₁ = −0.097, rising to 0.065 at ±0.3). But on B_high that landscape has a second, aliased minimum at a split near ±0.8 rad/μs, which no idle Kingston qubit carries. So Σδ is the robust readout here and the split is not. So the detuning is what γ_z was absorbing, not what γ_T1 was, and the fit-versus-calibration excess survives the repair intact.

The last column carries a warning of its own. The effective drive ω + 2Σδ is what a detuning-bearing model would hand to the F113 inversion, and at ω = 0.13 the chip cancels about 63% of what was applied on A_mid and 46% on B_high. So the inversion below reads the drive the experimenter set, not the drive the qubits felt, and the gap between them is not small.

### (3) The dephasing rates absorb the missing terms

The right comparison is not 1/T2. A `D[Z]` channel at rate γ_z decays coherences at 2γ_z, and 1/T2 = 2γ_z + 1/(2T1), so the Lindblad target is γ_z = (1/T2 − 1/(2T1))/2: 0.0012 and 0.0023 per μs for A_mid, 0.0004 and 0.0007 for B_high. The fitted 0.011 to 0.016 sits 7 to 30× above THAT, and the RMS around 0.3 says the same thing: the minimal model is still underfit, and γ_z is where the unmodelled structure lands. The F113-inversion is unaffected by that, because it reads the asymmetry of whatever L the fit returns; it is a faithful inversion of the model, not an independent measurement of the chip.

### (4) Sharpened F113-as-diagnostic interpretation

The downstream useful reading of this experiment:

> The F113-extracted γ_T1 is an **effective-T1 number that equates all bit_b-mixed broken-balance noise to the σ⁻ T1 channel**. It is at least stable against the channel most likely to have distorted it: adding a per-qubit detuning moves it by under 3%. That rules out the detuning, not the rest, so the gap to device-calibrated 1/T1 stays a bound on the absorbed non-T1 noise rather than a reading of the chip.

F113 makes this conversion **structural**: any polarity-asymmetry measurement on a Z-drive Bell-state protocol gives a single-number diagnostic that combines T1 with all other bit_b-mixed channels into one effective rate. This is useful as a quick hardware health check that doesn't require a full multi-channel noise-model fit.

## Limits

- **Minimal model**: the Z+T1 fit captures only two noise parameters. For richer noise models the F113 inversion still works structurally but the extracted γ_T1 absorbs more channels, making the "effective T1" reading drift further from the isolated-qubit T1.
- **No detuning term in the shipped model**: adding one leaves γ_T1 within 3% (see above), but it is itself a single-site Z, so it enters F113's ω_l as well as the fit. A detuning-bearing model therefore needs the inversion rewritten with ω_l → ω_l + 2δ_l before its asymmetry can be read the same way.
- **Per-pair, not per-qubit**: F113 gives one γ_T1 number per pair-run (assuming uniform-site rates). Per-qubit T1 calibration is more granular than what polarity asymmetry can resolve at this scope.
- **One trajectory class only (Z-drive + Bell)**: F113 inversion as T1-extraction works specifically for the H = (ω/2)·Z_drive + σ⁻ T1 noise family, with the drive on the steering qubit as here. Other Hamiltonian + dissipator combinations would need different inversion formulas (or to use F87 + F112 + F113 jointly as in Welle 5.B).

## Connection to existing readings

- **Multi-model analysis** ([`F112_HARDWARE_LENS_KINGSTON.md`](F112_HARDWARE_LENS_KINGSTON.md)): that survey identifies no channel, for reasons that partly apply here too, since the two share the missing detuning term. It is not independent evidence for anything in this document. The fit here is the minimal Z+T1 model on purpose, to make the F113 inversion clean.
- **F84 amplitude-damping correction**: F84 captures the F81-axis projection of amplitude damping. F113 here gives the polarity-axis projection: the same σ⁻ non-Hermiticity that breaks F81 also breaks F112, with magnitude given by F113.
- **Calibrated-T1 vs experiment-effective-T1 question**: this experiment makes the distinction concrete. Future hardware analyses can use F113 inversion as a quick check: "what γ_T1 does this protocol's polarity-asymmetry imply, and how does it compare to standalone characterization?"

## Reproduction

```
python -X utf8 simulations/f113_t1_extraction_kingston.py
```

Runs in ~15 seconds; produces the per-pair-run γ_T1 comparison table inline.
