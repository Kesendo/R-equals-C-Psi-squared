# F113 T1-Extraction on Kingston f95 Angle-Steering Data

**Status:** First application of F113 closed form as a hardware diagnostic: invert the F113 formula to extract a γ_T1 reading from the polarity-asymmetry measurement and compare to device-calibrated 1/T1. Self-consistency confirmed to machine precision; the fit-vs-calibration gap is READ as absorbed non-T1 noise, but it is not separable from epoch drift in the comparand (Reading (4)) (§ Limits: admitting a detuning leaves the fitted rate but not the extracted one).
**Date:** 2026-05-26
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Script:** [`simulations/f113_t1_extraction_kingston.py`](../simulations/f113_t1_extraction_kingston.py)
**Data:** [`data/ibm_f95_angle_steering_may2026/`](../data/ibm_f95_angle_steering_may2026/) (2 omega × 2 pair-runs × 6 t-points each, Kingston, 2026-05-16)
**Builds on:** F113 ([`F113_BREAK_MAGNITUDE_FORMULA.md`](F113_BREAK_MAGNITUDE_FORMULA.md), [`PROOF_F113_COEFFICIENT_DERIVATION.md`](../docs/proofs/PROOF_F113_COEFFICIENT_DERIVATION.md))

## Idea

F113 (Welle 4 closure) gives the polarity asymmetry of a Z-drive + amplitude-damping Lindblad in closed form:

    asymmetry = (4^N / 2) · Σ_l ω_l · (γ_pump,l − γ_T1,l)

The inverse: given a measured polarity asymmetry, known drive ω, and γ_pump = 0, extract a γ_T1 estimate:

    γ_T1_F113 = −asymmetry / ((1/2) · 4^N · ω)

One driven site, so the sum has a single term: at N=2 this is `γ_T1_F113 = −asymmetry / (8 · ω)`.
The all-sites-driven prefactor is a different number, `(N/2) · 4^N` (the
`PredictAsymmetryUniform` helper), which at N=2 is 16; folding the two together is the
notation the F113 registry entry corrects at its own hardware-application paragraph.

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

### (1) F113 self-consistency: confirmed to machine precision

`F113/fit = 1.000000` to the six decimals the table prints, across all 4 pair-runs (at full precision the four ratios sit within 4e-15 of 1, i.e. machine precision rather than bit-identical). The F113 closed form is a faithful re-derivation of the fitted Lindblad's polarity asymmetry; inversion recovers the same γ_T1 that produced the asymmetry. This is the structural self-consistency check the experiment was designed to perform, and it passes.

### (2) Fit above calibration, and the ω-trend is not readable

`fit/calib` ranges from 1.14 to 1.47. The two ω=0.25 runs do not bracket the two ω=0.13 runs (1.47 and 1.22 against 1.14 and 1.33), so the pair-to-pair spread is as large as any drive dependence, and this data cannot support a "stronger drive, more excess" reading.

And the excess is not resolved everywhere. A parametric bootstrap at the run's own shot budget (binomial on each of the 15 Pauli expectations, 2048 shots per basis, 30 replicas of the detuned fit) gives σ(γ_T1) = 0.0012 on A_mid ω=0.13, so its calibration value sits **0.54σ** from that same detuned fit (against the shipped fit in the table it would be 0.59σ). Note the 15 expectations come from 9 measured bases, three per basis, so they are not 15 independent binomials and the σ is approximate. That run shows no excess at all once the shot noise is counted, so only the larger ratios are worth reading. How much larger is not established here: the bootstrap was run for this one pair-run, the other three were never given a σ, and at 30 replicas the σ itself carries about 13% of its own scatter. Read the ordering, not a significance. Readout error sits on top as an unquantified systematic; it was not corrected here. IBM's device-class figure for Heron r2 is 1-3% per qubit, and the nearest Kingston calibration snapshot in the repo, 2026-05-05, eleven days before the run, puts this experiment's own four qubits at 0.49% (q13), 1.66% (q83) and 2.64% (q14), with q82 an outlier at 38.94% (0.51% in the June snapshot). So the class figure is not a conservative bound for this pair, and A_mid in particular is not covered by it.

The direction of the excess, where it is real, is expected and structural: the minimal Z + σ⁻ T1 model has two parameters, so any noise outside pure Z-dephasing gets absorbed by the σ⁻ channel, the only other one available.

The obvious suspect for that absorbed noise is a missing per-qubit detuning, since this chip's coherences rotate. Adding one free Z coefficient per qubit (12-start Nelder-Mead) settles it, and the answer is that γ_T1 is robust:

| pair-run | γ_T1 (Z+T1) | γ_T1 (+detuning) | γ_z | RMS | Σδ | ω + 2Σδ |
|---|---:|---:|---|---|---:|---:|
| ω=0.13 A_mid | 0.00577 | 0.00571 | 0.0163 → 0.0045 | 0.297 → 0.052 | −0.0408 | 0.048 |
| ω=0.13 B_high | 0.00615 | 0.00617 | 0.0123 → 0.0016 | 0.361 → 0.037 | −0.0298 | 0.070 |
| ω=0.25 A_mid | 0.00746 | 0.00743 | 0.0155 → 0.0044 | 0.303 → 0.066 | −0.0399 | 0.170 |
| ω=0.25 B_high | 0.00564 | 0.00564 | 0.0114 → 0.0014 | 0.366 → 0.043 | −0.0290 | 0.192 |

RMS improves by 4.6 to 9.7× and γ_z falls further, while γ_T1 moves by at most 1.1%. Σδ is consistent across three of the four runs at −0.029 to −0.041 per μs, i.e. the chip's own detuning opposes the applied drive. The exception is ω=0.25 B_high, where the shipped fit lands at Σδ = −0.029 but a wide multistart finds a lower minimum near Σδ = −0.87, in the same τ-grid alias band as the split (the shipped `fit_z_t1_detuned` reaches it at 30 and 60 starts with seed 2 but not with seeds 0 or 1, so which basin a run reports is seed-dependent). γ_T1 is identical to six digits in both (0.00564435) and γ_z agrees to four (0.0013569 vs 0.0013562), so nothing that reads a RATE below is affected. What is affected is anything reading ω+2Σδ: in the deeper basin the effective drive changes sign, so this row's contribution to the extracted-rate range in § Limits belongs to the shipped basin only. The SPLIT between the two sites is a different matter. It is not unobservable, as one might expect from the Bell coherence turning at the sum: the measured initial state carries single-excitation coherences (ρ₀₁ = 0.002 − 0.016i, ρ₀₂ = 0.016 − 0.045i) that turn at 2δ₁ and 2δ₀ separately, and least squares does resolve them, with a sharp minimum (SSE 0.0136 at δ₀−δ₁ = −0.097, rising to 0.065 at ±0.3). But on B_high that landscape has a second minimum, at individual detunings near ±0.8 rad/μs (a split δ₀−δ₁ of about 1.57, so read the ±0.8 as the δ_l and not as the split defined above). It is an aliasing artifact of the τ grid: B_high's minimum spacing is 3.71 μs, and π/Δt = 0.847 rad/μs is the scale it sits at, within a few percent. So Σδ is the robust readout here and the split is not. So the detuning is what γ_z was absorbing, not what γ_T1 was, and the fit-versus-calibration excess survives the repair intact.

The last column carries a warning of its own. The effective drive ω + 2Σδ is what a detuning-bearing model would hand to the F113 inversion, and at ω = 0.13 the chip cancels about 63% of what was applied on A_mid and 46% on B_high. So the inversion below reads the drive the experimenter set, not the drive the qubits felt, and the gap between them is not small.

### (3) The dephasing rates absorb the missing terms

The right comparison is not 1/T2. A `D[Z]` channel at rate γ_z decays coherences at 2γ_z, and 1/T2 = 2γ_z + 1/(2T1), so the Lindblad target is γ_z = (1/T2 − 1/(2T1))/2: 0.0012 and 0.0023 per μs for A_mid, 0.0004 and 0.0007 for B_high. The fitted 0.011 to 0.016 sits 6.9 to 30.7× above THAT, and the RMS around 0.3 says the same thing: the minimal model is still underfit, and γ_z is where the unmodelled structure lands. The F113-inversion is unaffected by that, because it reads the asymmetry of whatever L the fit returns; it is a faithful inversion of the model, not an independent measurement of the chip.

### (4) Sharpened F113-as-diagnostic interpretation

The downstream useful reading of this experiment:

> The F113-extracted γ_T1 is an **effective-T1 number that equates all bit_b-mixed broken-balance noise to the σ⁻ T1 channel**. The channel most likely to have distorted it is a per-qubit detuning, and admitting one moves the FIT's γ_T1 by under 3%. It does not leave the extracted number where it was: the inversion divides by the drive the experimenter set, while the detuned model's own effective drive is ω + 2Σδ, so the extracted rate moves from 0.00577 to 0.00213 on A_mid ω=0.13, a factor 2.7, and by 23 to 63% across the four runs. So what survives the detuning repair is the fitted rate, not the reading, and any detuning-bearing model has to rewrite the inversion with ω_l → ω_l + 2δ_l first. Nor is the gap to device-calibrated 1/T1 a bound on the absorbed non-T1 noise, because the comparand moves at least as much as the gap: the metadata's T1 is IBM's published backend calibration (`cal_snapshot_start` and `cal_snapshot_end` are bit-identical, so nothing brackets the run), and q13, which appears both here and in the F120 flight, is recorded at 232.9 μs in this dataset and 415.8 μs in that one, a factor 1.79 against a measured excess of 1.14 to 1.47; F120 also caught q13 at ≈197 μs inside a single morning. Read the excess as unresolved rather than as a bound.

F113 makes this conversion **structural**: any polarity-asymmetry measurement on a Z-drive Bell-state protocol gives a single-number diagnostic that combines T1 with all other bit_b-mixed channels into one effective rate. This is useful as a quick hardware health check that doesn't require a full multi-channel noise-model fit.

## Limits

- **Minimal model**: the Z+T1 fit captures only two noise parameters. For richer noise models the F113 inversion still works structurally but the extracted γ_T1 absorbs more channels, making the "effective T1" reading drift further from the isolated-qubit T1.
- **No detuning term in the shipped model**: adding one leaves γ_T1 within 3% (see above), but it is itself a single-site Z, so it enters F113's ω_l as well as the fit. A detuning-bearing model therefore needs the inversion rewritten with ω_l → ω_l + 2δ_l before its asymmetry can be read the same way.
- **One number per pair-run, and the limit is the FIT's not the asymmetry's**: the shipped model imposes a single γ_T1 on both qubits, so one number is all there is to extract. The asymmetry itself is finer than that: its weights are the ω_l, and in the shipped model the drive sits on one qubit (`DRIVE_SITE = 0`) with the other's ω exactly 0, so the undriven qubit's γ_T1 does not enter: on the A_mid ω=0.13 run the asymmetry is −6.0028e−3 whatever the undriven rate is, over a sweep from 0.001 to 0.5. That weighting belongs to the model, not to the chip. §2 above shows the model's ω is not the drive the qubits felt, and any per-qubit detuning sends ω_l → ω_l + 2δ_l, which gives the second qubit a nonzero weight; so this is a statement about which rate the fit's asymmetry is sensitive to, not an attribution to a named physical qubit. The table's `γ_T1 calib` column stays the pair mean, the fit's own comparand. The within-pair spread in 1/T1 is 30 % on A_mid (q82/q83 = 174.4/227.2 μs) and 16 % on B_high (q13/q14 = 232.9/200.6 μs), so a different comparand would move the `fit/calib` ratios; none is folded in here.
- **One trajectory class only (Z-drive + Bell)**: F113 inversion as T1-extraction works specifically for the H = (ω/2)·Z_drive + σ⁻ T1 noise family, with the drive on the steering qubit as here. Other Hamiltonian + dissipator combinations would need different inversion formulas (or to use F87 + F112 + F113 jointly as in Welle 5.B).

## Connection to existing readings

- **Multi-model analysis** ([`F112_HARDWARE_LENS_KINGSTON.md`](F112_HARDWARE_LENS_KINGSTON.md)): that survey identifies no channel, for reasons that partly apply here too, since the two share the missing detuning term. It is not independent evidence for anything in this document. The fit here is the minimal Z+T1 model on purpose, to make the F113 inversion clean.
- **F84 amplitude-damping correction**: F84 captures the F81-axis projection of amplitude damping. F113 here gives the polarity-axis projection: the same σ⁻ non-Hermiticity that breaks F81 also breaks F112, with magnitude given by F113.
- **Calibrated-T1 vs experiment-effective-T1 question**: this experiment makes the distinction concrete. Future hardware analyses can use F113 inversion as a quick check: "what γ_T1 does this protocol's polarity-asymmetry imply, and how does it compare to standalone characterization?"

## Reproduction

```
python -X utf8 simulations/f113_t1_extraction_kingston.py
```

Runs in ~33 seconds; produces the per-pair-run γ_T1 comparison table inline.
