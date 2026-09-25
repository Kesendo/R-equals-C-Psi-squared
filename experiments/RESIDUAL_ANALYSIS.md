# Residual Analysis: Late-Time Coherence Anomaly in IBM Torino Qubit 52

<!-- Keywords: IBM Torino qubit 52 residual coherence anomaly, late-time
off-diagonal directional consistency, Monte Carlo shot-noise null, static
late-time phase, measurement offset SPAM readout asymmetry basis pulse leak,
Re equals minus Im, selection-sensitive tail slope, TLS coupling hypothesis,
plus minus preparation control, R=CPsi2 residual analysis -->

**Status:** Tier 2, exploratory. The excess over the shot-noise null is real (9.5 σ); its direction is a static offset with equal ⟨X⟩ and ⟨Y⟩ parts, the pattern of a measurement (SPAM) offset, which a readout asymmetry and a basis-pulse leak of ⟨Z⟩ both leave and the record does not rank; its mechanism is open until a |+⟩/|−⟩ control and a readout-assignment measurement run.
**Date:** 2026-02-09
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [IBM Quantum Tomography](IBM_QUANTUM_TOMOGRAPHY.md), [Fixed Point Shadow](FIXED_POINT_SHADOW.md)

> The correct response to this data is not belief or disbelief. It is: measure again.

---

## What this document is about

After measuring a qubit's decoherence on IBM hardware, the raw data was
re-examined for structure in the noise. Coherence that should have faded into
shot noise had not: its late mean sits 9.5 standard deviations above a
shot-noise null, and in all 17 late samples it points into the same quadrant.
Three explanations were proposed: a calibration artifact, coupling to a material
defect, or an external coherent drive. The record itself says more than it was
first asked: the late component does not rotate, it does not point where the ¼
crossing's fixed point pointed (and through the crossing itself the coherence kept
turning), and it carries equal offsets in ⟨X⟩ and ⟨Y⟩, the pattern of a
measurement offset. The calibration artifact leads. The
March 2026 follow-up on two other qubits ([Fixed Point Shadow](FIXED_POINT_SHADOW.md))
measured their detuning and closed its window before their coherence had decayed
to a residual of this size, so it could not see one. The controls that would narrow the
mechanism have not been run. The document is kept as a case study in careful
anomaly investigation.

---

## Abstract

Reanalysis of IBM Torino qubit 52 tomography data found late-time coherence
(t/T₂(echo) ≥ 1.25) above all 10,000 runs of a shot-noise Monte Carlo null
(late mean 0.01852 against 0.00859 ± 0.00105, 9.5 σ), with 17 of 17 samples at
t/T₂(echo) ≥ 1 in the fourth quadrant (Re+ / Im−). A tail slope is positive but
selection-sensitive (p = 0.053 at the ≥ 1.5 cut, 0.014 to 0.033 at the others),
and the reported correlation with the distance to ¼ (r = −0.9955) reuses |ρ₀₁|
on both axes. Three hypotheses were proposed: SPAM error (State Preparation And
Measurement, systematic bias from imperfect calibration), TLS coupling
(Two-Level System, a material defect in the chip substrate that can exchange
energy with the qubit), and an external coherent drive; the companion page added
the ¼ boundary itself as a source. The record closes that last reading as it was
proposed: at the crossing the coherence pointed at +79° and +89° and kept turning
for another 150 μs, and the late direction, −48.4° ± 5.5°, lies 4.4 standard
errors from the −24° the fixed point pointed. The late phase is static
(−0.040°/μs, p = 0.18) and the late ⟨X⟩ and ⟨Y⟩ offsets are equal within 0.7 to
1.3 σ across the cuts, the pattern of a measurement offset: a readout asymmetry of
about 2.4% leaves it, and so does an amplitude error in the basis-change pulse,
which leaks ⟨Z⟩ into both readings. The record does not rank the two; the first
sample's ⟨Y⟩ stands against the asymmetry, the extrapolated ⟨Z⟩, with the relaxation
time fitted, for it. SPAM is
the leading candidate. The March 2026 follow-up on Q80 and Q102 measured their detuning
and could not see a residual of this size. A |+⟩/|−⟩ preparation test and a measured
readout-assignment matrix would narrow the mechanism.

## Discovery Date
2026-02-09 (same-day reanalysis of tomography data)

## Context

After completing the [IBM Quantum Tomography](IBM_QUANTUM_TOMOGRAPHY.md)
experiment, the raw density matrices were reanalysed for structure in the
residuals. The original experiment measured C·Ψ crossing the ¼ boundary and
compared the crossing with the generalized crossing equation, a qualitative
comparison on the same record. This analysis asks a different question: **Is
there structure in the data beyond the expected exponential decay?**

The time axis throughout is t/T₂(echo), the delay over the calibration's
Hahn-echo value T₂(echo) = 298.2 μs, as stored in the record. The qubit's own
free decay runs on T₂* = 110.7 μs, so t/T₂(echo) = 1 is t = 2.7 T₂*.

## Finding 1: Excess Late-Time Coherence

### Null Hypothesis Test

A Monte Carlo simulation was constructed:
- 10,000 synthetic tomography experiments
- Parameters: T₁ = 221.2 μs, T₂* = 110.2 μs (the recorded input; the record's own free fits give 110.1 and 110.7 μs), 8192 shots, 25 delay points, and the ideal amplitude ½ at t = 0
- Shot noise modelled as binomial sampling in the tomography bases
- One random phase per run (an unknown detuning phase)

The original script is not in the repository. The null rebuilds from this
description (the rebuild runs as a gate in
[`test_torino_shadow_records.py`](../simulations/tests/test_torino_shadow_records.py)),
and it lands on the recorded numbers. It lands there with the ideal amplitude ½;
with the record's fitted amplitude 0.48 the null's late mean drops to 0.00849,
which only widens the separation.

| Metric | IBM Torino (real) | Null (mean ± std) |
|--------|------------------|----------------------|
| Mean \|ρ₀₁\| for t/T₂(echo) ≥ 1.25 (15 samples, t > 370 μs) | 0.01852 | 0.00859 ± 0.00105 rebuilt (0.00861 ± 0.00105 recorded) |
| Runs reaching the hardware value | 0 of 10,000 | z = 9.5 |
| Samples above the per-sample 99th percentile | 9 of 25 (two early, seven late; the one at t/T₂(echo) = 2.125 sits at the percentile's own sampling edge) | 0.25 expected |

The strongest anomaly is at t/T₂(echo) = 2.625, where |ρ₀₁| = 0.0363, 2.2× the
null's 99th percentile there (0.0167 rebuilt, 0.0166 recorded).

The excess does not hang on the decay constant: with T₂* = 100, 110, 120 or
130 μs in the null, none of 4,000 runs reaches the hardware value. Only a
decay constant of 150 μs, far off the fitted 110 μs, lets a handful through.

### Important Caveat

The null model is simple: exponential decay plus shot noise. It does not include
TLS coupling, 1/f noise, readout error, or non-Markovian dynamics. Any of these
could close the gap. The point is not that the effect is unexplainable, but that
it is **quantifiably present** and deserves investigation.

An offset of the late samples' own size (δ = 0.024, their mean ⟨X⟩ and ⟨Y⟩) added
to the same null reproduces the late mean (0.0186 ± 0.0028 against the measured
0.0185). That holds by construction for any static offset of that size: it says the
excess is this mean offset and no extra scatter, not what made it. Finding 2 shows
why a measurement offset is the natural candidate, and why the record cannot say
which form.

## Finding 2: Directional Coherence (17/17 Sign Consistency)

This is the most striking anomaly.

For all 17 data points at t/T₂(echo) ≥ 1.0:
- **Re(ρ₀₁) > 0 in 17/17 measurements**
- **Im(ρ₀₁) < 0 in 17/17 measurements**

The residual coherence always points into the fourth quadrant of the complex plane.

| Component | Mean | Std | All same sign? |
|-----------|------|-----|---------------|
| Re(ρ₀₁) | +0.01137 | 0.00499 | Yes (17/17 positive) |
| Im(ρ₀₁) | −0.01279 | 0.00805 | Yes (17/17 negative) |

This is remarkable because decoherence should leave a residual of random phase: at
t = 3 T₂(echo) the prepared coherence is e^(−8) of its start, and what is left
should scatter around zero. Instead the residual has a fixed direction, −48.4° for
the mean. (A count of 34 independent fair signs gives (½)³⁴ ≈ 6 × 10⁻¹¹; that is
the arithmetic of coins, not a probability for anything the hardware does.)

The two components say more together than apart. In Pauli terms the late samples
carry ⟨X⟩ = 2 Re ρ₀₁ = 0.0227 ± 0.0025 and ⟨Y⟩ = −2 Im ρ₀₁ = 0.0256 ± 0.0040, 9.1 σ
and 6.4 σ from zero. They are equal within the noise at every cut:

| Cut | Samples | ⟨X⟩ − ⟨Y⟩ | Mean direction |
|-----|---------|-----------|----------------|
| ≥ 1.0 | 17 | −0.0028 ± 0.0041 (0.7 σ) | −48.4° ± 5.5° |
| ≥ 1.25 | 15 | −0.0045 ± 0.0044 (1.0 σ) | −50.2° ± 5.8° |
| ≥ 1.5 | 13 | −0.0063 ± 0.0049 (1.3 σ) | −51.9° ± 5.6° |

Re ρ₀₁ = −Im ρ₀₁ is the pattern of a measurement offset, and two simple forms leave
it. A readout asymmetry δ = P(0|1) − P(1|0), read out in any basis, adds +δ to the
measured expectation, so ρ₀₁ = (⟨X⟩ − i⟨Y⟩)/2 gains (δ − iδ)/2, at −45°, in the
fourth quadrant for the usual δ > 0. An amplitude error in the basis-change pulse
leaks ⟨Z⟩ into the X and the Y reading with one coefficient (the phase gate that
makes one basis change out of the other commutes with Z), again at −45° for a pulse
that falls short. The stored readings are raw count ratios, n/8192, with no
mitigation, so either would reach ρ unchanged. Fits of the late samples alone, where
a simple model fits (a constant, or a constant plus the decaying rotation; χ²/dof
1.5 to 2), put the offset between −49° and −56°, each within one and a half
standard errors of −45°. A fit of the whole record, a detuned decay plus a constant
complex offset, cannot weigh it: its early part misfits by χ²/dof = 34, and scaled
by that misfit the offset's Δχ² = 117 for two parameters is not significant
(p = 0.19).

The record does not rank the two forms. On the 34 late readings the constant
asymmetry (δ = 0.024) leaves χ² = 50.5 and the leak (0.039·⟨Z⟩) leaves 45.4, for 33
degrees of freedom each, a lean toward the leak (a likelihood ratio of about 13) and
not a verdict; the late ⟨X⟩ and ⟨Y⟩ slopes on ⟨Z⟩, 0.04 ± 0.04 and 0.15 ± 0.06, do not
separate a leak (0.039) from an asymmetry (0) either. At t = 0 an ideal |+⟩ has ⟨Y⟩ = 0, so the asymmetry
should read ⟨Y⟩ = +0.024; the record reads −0.018 ± 0.011, 3.8 σ below, which
inside the asymmetry reading needs a preparation or pre-rotation phase error of
about 2.6°, while the leak expects ⟨Y⟩ ≈ 0 there (1.6 σ). The relaxation curve of
⟨Z⟩, extrapolated to t = 0 with its relaxation time fitted (241 μs), reads
+0.021 ± 0.008, as the asymmetry predicts (+0.006 ± 0.007 with the calibration's T₁
held fixed, a fit worse by Δχ² ≈ 14), and a preparation made with the same short
pulse would start at ⟨Z⟩ = +0.039, 3.0 σ above the first sample's +0.005 ± 0.011.
The −45° pattern names the family, not the member.

The same late component appears in [Absorption Theorem on IBM Hardware](IBM_ABSORPTION_THEOREM.md)
§4 as the slow tail of the coherence envelope (amplitude 0.013, "readout offset
(most likely)"); a decay plus a constant puts the constant at 0.0123 ± 0.0051.

## Finding 3: A Mild Rise, Cut-Sensitive

A linear fit to |ρ₀₁| for t/T₂(echo) ≥ 1.5 shows a **positive slope**:

```
Slope = +0.00819 per T₂(echo) unit   (2.75 × 10⁻⁵ per μs)
```

It is not a clean trend. The amplitudes go up and down, and the verdict depends on
where the tail is cut:

| Cut | Samples | Slope per T₂(echo) | Two-sided p | 95% interval |
|-----|---------|--------------------|-------------|--------------|
| ≥ 1.0 | 17 | +0.00699 | 0.014 | [+0.00164, +0.01234] |
| ≥ 1.25 | 15 | +0.00878 | 0.015 | [+0.00203, +0.01554] |
| ≥ 1.5 | 13 | +0.00819 | 0.053 | [−0.00013, +0.01651] |
| > 1.5 | 12 | +0.01041 | 0.033 | [+0.00103, +0.01978] |

A rising off-diagonal element is also no violation of Markovian dynamics: a
coherent drive, or a stationary state that holds coherence, raises one under a
perfectly memoryless Lindblad generator. And a measurement-side error that leaks
the ground-state population into the transverse readings (Hypothesis 1 below) rises
with that population. So the slope witnesses neither memory nor anything against
a calibration artifact.

## Finding 4: Qubit Detuning and Phase Structure

A frequency detuning was extracted from the early high-SNR phase data, the first
five samples: their phases +1.1°, −62.4°, −146.7°, +79.4°, +89.0° unwrap to a slope
of −0.0357 rad/μs, −5.68 kHz (a detuning representative δ_f = +5.68 kHz under the
convention ρ₀₁ ∝ e^(−iδωt)), a period of about 176 μs. The samples sit on a
37.28 μs grid (T₂(echo)/8), so the frequency is known only modulo 26.82 kHz; the
value quoted is the representative nearest zero.

The rotation is not steady. The step-to-step turns run −63.5°, −84.3°, −133.9°,
+9.6°, −137.0°, −98.9°, −113.1°, −142.8° through t = 300 μs, and after derotation
residual phase remains at intermediate times (t/T₂(echo) = 0.375 to 1.375): the
frequency wanders, or the signal has more than one component.

Then it stops. From t/T₂(echo) = 1 on, the phase stands still: the 17 late phases
drift by −0.040° per μs (p = 0.18), where a detuning of the early size would turn
them by −76° every sample. Quantised to 30° bins, 9 of the 15 phases at
t/T₂(echo) ≥ 1.25 fall at −60° and 4 at −30°. So the late direction is not the
detuning frozen into place: a detuning moves, and this component does not. Nor was
it set at the ¼ crossing: the crossing falls at 114.7 μs, between the samples at
111.8 and 149.1 μs, where the phase reads +79.4° and +89.0°, and it keeps turning
through −48.1°, −146.9° and +100.0° until 261 μs. The early detuned coherence decays
away and a static offset is what remains.

## Finding 5: Spacing Between Late-Time Peaks

Local maxima in late-time |ρ₀₁| occur at:

| Peak | t/T₂(echo) | \|ρ₀₁\| |
|------|-------|---------|
| 1 | 1.500 | 0.01964 |
| 2 | 1.750 | 0.01234 |
| 3 | 2.000 | 0.02210 |
| 4 | 2.625 | 0.03629 |

Peak spacings: 0.250, 0.250 and 0.625 T₂(echo). They looked like a pattern (two
equal gaps of T₂/4, and a third exactly 2.5 times the first). They are the sampling
grid: every sample sits on a multiple of 0.125 T₂(echo), so every spacing is one,
and the peaks of a noisy tail on that grid carry no revival period.

## A Separate Simulator Fixture

The data folder also holds a simulator run of the same protocol
(`simulator_test_20260209_125106.json`), made with different parameters
(T₁ = 200 μs, T₂ = 150 μs, its time axis t/150 μs). Its late samples
(t/T₂ ≥ 1.25 on its own axis) average 0.0726, well above the hardware's 0.0185:
that model holds more late coherence than the chip. It stores magnitudes only, so it
says nothing about direction, and with other parameters it is no matched control;
the shot-noise null of Finding 1 is the comparison this analysis rests on.

## Three Hypotheses

### H1: Systematic SPAM Error (Most Conservative)

Tomography calibration has a small, constant error. When the true signal is zero,
the fixed offset dominates. This would explain the directional consistency.

**Why it leads.** Both simple measurement-side forms put equal offsets into ⟨X⟩
and ⟨Y⟩, at −45°: a readout asymmetry adds +δ to every reading, and an amplitude
error in the basis-change pulse leaks ⟨Z⟩ into both. The late samples carry equal
offsets within 0.7 to 1.3 σ (Finding 2). Which form it is, the record does not say:
each meets one tension of 3 to 4 σ in the first sample (the leak only if the
preparation shares the pulse), and the leak, which climbs with ⟨Z⟩ from 0.005 to
0.71 as the qubit relaxes, would also give the mild rise of Finding 3.

**Test:** Measure with |+⟩ and |−⟩ initial states. A SPAM offset is independent of the
initial state; a coherence carried from the preparation flips sign. An external drive
locked in phase would keep its sign too, so the test narrows the mechanism to one of
two families, and a measured readout-assignment matrix identifies the readout part.
Preparing |0⟩ and |1⟩ and reading each in X and Y separates the two SPAM forms
directly: a readout asymmetry reads +δ for both, a ⟨Z⟩ leak flips sign with ⟨Z⟩.

### H2: Two-Level System Coupling (Known Physics)

A TLS defect in the substrate is coupled to qubit 52 and feeds coherence back into
the system. TLS coupling is well documented in superconducting qubits and can
produce coherence revivals with a preferred phase.

**Status:** Not excluded, and nothing in the record requires it. A TLS that drifts
over hours looks static across a sub-millisecond delay window, so a fixed direction
does not argue against it; but the late component shows no rotation and no revival,
and a measurement offset explains its shape without a defect.

**Test:** Repeat the measurement on qubit 52 days later, and measure neighbouring
qubits: TLS effects are local to individual qubits, and a changed direction on its
own would not tell a TLS from calibration drift.

### H3: External Coherent Coupling (Speculative)

An external source is driving coherence into the qubit at a fixed phase. This would
explain excess coherence, fixed direction and a rising trend with one mechanism.

**Scope:** Such sources exist and are mundane (microwave leakage or crosstalk, chip
and package modes), so the refrigerator temperature does not exclude one; none was
identified or fitted here. Its test below needs late residuals on several qubits, and
the March run, which closed its window before the other qubits' coherence had
decayed, did not measure any.

**Test:** Multi-qubit correlation. If independent qubits on the same chip show
correlated excess coherence with the same directional signature, it cannot be local
TLS or individual SPAM.

## March 2026 Hardware Test Plan

Based on these findings, the March run should include:

1. **Reproduce:** Same qubit 52, same protocol. Does the 17/17 directional signature persist?
2. **SPAM discrimination:** Run with |+⟩ and |−⟩ initial states. Does the excess flip sign?
3. **Multi-qubit:** At least 5 qubits (ideally 10+). Check Re/Im sign consistency independently for each.
4. **Cross-correlation:** Do different qubits show correlated excess coherence?
5. **Extended time range:** Push to 5×T₂ or beyond. Does the rising trend continue?

If H1 (SPAM) survives: The anomaly is a calibration artifact. Document and move on.
If H1 is killed and H2 survives: Non-Markovian dynamics worth studying but conventional.
If both H1 and H2 are killed: We have a genuine open question.

The March run on 2026-03-09 carried out item 3 on two of the five qubits it asks
for, Q80 and Q102 ([Fixed Point Shadow](FIXED_POINT_SHADOW.md)), with its late
window timed on each qubit's echo time. That window closed while their coherence was
still twice the size of Q52's residual: the run measured each qubit's detuning (Q102
turns at tens of kHz, Q80 drifts slowly) and could not see a residual. Items 1, 2, 4
and 5 have not been run.

## Raw Numbers for Reference

Late-time off-diagonal elements (t/T₂(echo) ≥ 1.0):

```
t/T₂(echo)   Re(ρ₀₁)      Im(ρ₀₁)      |ρ₀₁|
1.000   +0.011597    -0.010742    0.015808
1.125   +0.013794    -0.004883    0.014633
1.250   +0.002930    -0.000977    0.003088
1.375   +0.015747    -0.010376    0.018858
1.500   +0.012939    -0.014771    0.019637
1.625   +0.007568    -0.007568    0.010703
1.750   +0.012329    -0.000488    0.012339
1.875   +0.002319    -0.011719    0.011946
2.000   +0.006104    -0.021240    0.022100
2.125   +0.011841    -0.012939    0.017540
2.250   +0.012817    -0.007080    0.014643
2.375   +0.007080    -0.020630    0.021811
2.500   +0.017578    -0.020508    0.027010
2.625   +0.021240    -0.029419    0.036285
2.750   +0.013184    -0.003906    0.013750
2.875   +0.016479    -0.016968    0.023653
3.000   +0.007690    -0.023193    0.024435
```

## Connection to R = CΨ² Framework

This analysis began as a search for structure beyond the C·Ψ = ¼ crossing, and the
excess lives after the crossing. Its link to the boundary is algebraic only: the saved
analysis computes Ψ = 2|ρ₀₁| and C·Ψ = C·2|ρ₀₁|, so the near-perfect correlation of
the residual with the distance to ¼ (r = −0.9955 over the 17 late samples) is one
quantity on both axes. And on this record the crossing leaves nothing the shadow
proposal predicted: the coherence pointed at +79° and +89° either side of it and kept
turning for another 150 μs, and the late direction lies 4.4 standard errors from the
fixed point's. The residual carries no information about the boundary that this
record can show; what it carries is information about the qubit and its measurement.

## Epistemic Status

- **Established:** Excess coherence beyond the shot-noise null (0 of 10,000 runs, 9.5 σ; the null rebuilt from its recorded description, robust to its decay constant)
- **Established:** Directional consistency Re+/Im− in 17/17 late samples, with a static phase and ⟨X⟩ and ⟨Y⟩ offsets equal within 0.7 to 1.3 σ across the cuts
- **Descriptive:** A mild rise in the late tail, significant or not depending on the cut
- **Leading candidate:** A measurement (SPAM) offset, the family [Absorption Theorem on IBM Hardware](IBM_ABSORPTION_THEOREM.md) §4 reads in the same tail; within it a readout asymmetry and a basis-pulse leak of ⟨Z⟩ fit the record about equally, each with one first-sample tension
- **Unknown:** Whether it is SPAM, TLS, or something else, until a |+⟩/|−⟩ control and a readout-assignment measurement run
- **Not claimed:** A non-Markovian witness or any property of the ¼ boundary

The numbers this page rests on are recomputed from the raw record by
[`test_torino_shadow_records.py`](../simulations/tests/test_torino_shadow_records.py).
