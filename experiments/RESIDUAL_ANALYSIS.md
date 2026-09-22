<!-- QUARTER-CURRENT -->
# Finite residual analysis for IBM Torino qubit 52

Current reading: the saved residuals, phase directions, and null calculations
belong to the named Q52 dataset and preprocessing choices.  They motivate noise
and drift hypotheses but do not establish a scalar-boundary mechanism or an
ontological change.

**Status:** Historical exploratory analysis; finite record retained, Q52 mechanism open

The completed Q80/Q102 comparison rejected a universal-boundary reading but did
not resolve the Q52 mechanism. Its Q80 phase-compatible fits do not resolve
the Q52 magnitude excess. The current result is [Q52 Residual Record](FIXED_POINT_SHADOW.md).

<!-- QUARTER-HISTORICAL -->
**Historical record:** the post-crossing interpretation below is retained with
the finite measurements that prompted it.

# Residual Analysis: Late-Time Coherence Anomaly in IBM Torino Qubit 52

<!-- Keywords: IBM Torino qubit 52 residual coherence anomaly, late-time
off-diagonal directional consistency, Monte Carlo null hypothesis p<0.0001,
rising coherence trend post-decoherence, fourth quadrant phase clustering,
SPAM error TLS coupling hypotheses, revival peak spacing T2/4, non-Markovian
coherence revival, R=CPsi2 residual analysis -->

> **Restoration note (March 14, 2026):** Originally written 2026-02-09, deleted March 12,
> restored March 14. The finite analysis is retained; the later cross-qubit comparison
> rejected a universal-boundary reading but did not resolve the Q52 mechanism
> (see [Q52 Residual Record](FIXED_POINT_SHADOW.md)).

**Status:** Historical exploratory analysis; finite record retained, Q52 mechanism open
**Date:** 2026-02-09
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [IBM Quantum Tomography](IBM_QUANTUM_TOMOGRAPHY.md), [Q52 Residual Record](FIXED_POINT_SHADOW.md)

---

## What this document is about

After measuring a qubit's decoherence on IBM hardware, the raw data was
re-examined for structure in the residuals. The finite record has a common
late-time quadrant and exceeds a recorded model described as exponential decay,
binomial shot sampling, and one random phase per synthetic run.
The positive tail slope and boundary-distance correlation depend on the
chosen preprocessing and do not identify a mechanism. The March 2026
comparison rejected a universal direction; its Q80 phase-compatible fits
do not resolve the Q52 magnitude excess. The dated analysis below is kept
as a case study in hypothesis generation, not as a current causal verdict.

---

## Abstract

Reanalysis of IBM Torino qubit 52 tomography data found a finite late-time
excess over a 10,000-draw recorded model combining exponential decay, binomial
shot sampling, and one random phase per synthetic run, with 17/17 rows
at `t/T2_echo >= 1` in one quadrant. The 13-row `>= 1.5` tail fit has slope
`+0.00819/T2_echo`, two-sided p = 0.0531, and a 95% interval that crosses
zero. The reported `r = -0.9955` boundary-distance correlation reuses
`|rho_01|` through `C*Psi` and is not independent boundary evidence. No
Q52-fitted time-dependent detuning/drift, SPAM, TLS, or memory alternative is
compared by the retained null description. The
March comparison found a different common direction on Q80 and no common
direction on Q102, rejecting the former universal reading without closing
the Q52 mechanism.

## Discovery Date
2026-02-09 (same-day reanalysis of tomography data)

## Status
**Exploratory.** Anomaly detected and quantified. Cause unknown. Three hypotheses proposed, none confirmed. March 2026 hardware run designed to discriminate.

## Context

After completing the IBM Quantum Tomography experiment (see IBM_QUANTUM_TOMOGRAPHY.md), the raw density matrices were reanalysed for structure in the residuals. The original experiment measured a C·Ψ crossing and compared it with the generalized crossing equation; that same-record comparison is qualitative crossing evidence, not an independent precision confirmation. This analysis asks a different question: **Is there structure in the data beyond the recorded exponential-decay null?**

## Finding 1: Excess relative to the recorded 10,000-draw null

### Null Hypothesis Test

A Monte Carlo simulation was constructed:
- 10,000 synthetic tomography experiments
- Recorded null inputs: T₁ = 221.2 μs, T₂* = 110.2 μs, 8192 shots, 25 delay points
- Shot noise modelled as binomial sampling across 3 tomography bases (X, Y, Z)
- Random phase per run (detuning)

The null producer is not retained in this repository. Its ensemble counts,
spread, and percentile threshold below are therefore historical recorded values;
only the hardware statistic and its selection are recomputed here from the saved
tomography rows.

**Result:** The saved IBM Torino late-time statistic (`t/T2_echo >= 1.25`, 15 stored rows) exceeded all 10,000 draws of this particular null.

| Metric | IBM Torino (real) | Recorded null ensemble (mean ± std) |
|--------|------------------|----------------------|
| Mean \|ρ₀₁\| for t/T2_echo >= 1.25 (15 rows) | 0.01852 | 0.00861 ± 0.00105 |
| Recorded exceedances | 0 of 10,000 | - |
| Points above 99th percentile | 9 of 25 | 0.25 expected |

At `t/T2_echo = 2.625`, the saved hardware row has `|rho_01| = 0.036285`.
The historical null table labels that value 2.2x its recorded 99th-percentile
threshold of 0.0166; without the producer, that percentile is not rerun here.

### Important Caveat

The recorded null combines exponential decay, binomial shot sampling, and one random phase per synthetic run. It does not include a Q52-fitted time-dependent detuning or drift, SPAM, TLS, memory, or other calibrated hardware alternative. Zero exceedances in this finite ensemble therefore describe separation from this null, not a p-value against all relevant hardware models.

## Finding 2: Directional Coherence (17/17 Sign Consistency)

This is the most striking anomaly.

For all 17 data points at `t/T2_echo >= 1.0`:
- **Re(ρ₀₁) > 0 in 17/17 measurements**
- **Im(ρ₀₁) < 0 in 17/17 measurements**

The residual coherence always points into the fourth quadrant of the complex plane.

| Component | Mean | Std | All same sign? |
|-----------|------|-----|---------------|
| Re(ρ₀₁) | +0.01137 | 0.00499 | Yes (17/17 positive) |
| Im(ρ₀₁) | -0.01279 | 0.00805 | Yes (17/17 negative) |

A naive calculation with 34 independent, equiprobable signs gives `(1/2)^34`, but the tomography coordinates share one acquisition and the relevant hardware null is not an independent-sign model. The durable observation is the finite 17/17 sign record, not that naive number as a calibrated p-value.

## Finding 3: Positive slope on one selected tail

A linear fit to \|ρ₀₁\| for `t/T2_echo >= 1.5` shows a **positive slope**:

```
Slope = +0.00819 per T2_echo unit
```

For the 13 rows selected by `t/T2_echo >= 1.5`, the two-sided slope p-value is 0.0531 and the 95% interval crosses zero. The amplitudes are non-monotone, and the slope changes when the endpoint convention changes. Markovian GKSL dynamics with coherent Hamiltonian evolution can also create or transiently increase a chosen off-diagonal element, so this finite slope is not a non-Markovianity witness.

## Finding 4: Qubit Detuning and Phase Structure

Using the first five saved rows and `np.unwrap` on phases recomputed from the raw complex ρ₀₁ coordinates, the fitted phase slope is `m = -0.035696 rad/us` (`-5.681 kHz`). Under `rho_01 ~ exp(-i delta_omega t)`, the selected near-zero detuning representative is therefore `delta_f = +5.681 kHz`. Those five delays lie on a `37.280936 us` grid, so frequency is identified only modulo `26.823361 kHz`; the five-row and unwrapping choices are part of this descriptive extraction, and no absolute detuning is identified without an additional physical prior.

After derotation, residual phase remains at intermediate times (`t/T2_echo = 0.375` to `1.375`). That is compatible with nonlinear drift or additional components but does not identify either mechanism.

Late-time phases quantised to 30° bins show strong clustering:
- -60°: 9 of 15 values (60%, naive uniform six-bin reference ~17%)
- -30°: 4 of 15
- 0°: 1 of 15
- -90°: 1 of 15

## Descriptive Reading 5: Spacing Between Selected Late-Time Local Maxima

Local maxima in late-time \|ρ₀₁\| occur at:

| Peak | t/T2_echo | \|ρ₀₁\| |
|------|-------|---------|
| 1 | 1.500 | 0.01964 |
| 2 | 1.750 | 0.01234 |
| 3 | 2.000 | 0.02210 |
| 4 | 2.625 | 0.03629 |

Peak spacings on the sampled grid are 0.250×T2_echo, 0.250×T2_echo, and 0.625×T2_echo. These selected local maxima do not establish a revival period.

## Separate Simulator Fixture

The retained fixture identifies itself only as `SIMULATOR_TEST`; no saved producer,
backend, or noise configuration establishes how it was generated. The comparison
below therefore uses only its saved analysis rows.

The saved February simulator fixture and hardware record each contain 15 rows at
their stored normalized cutoff `>= 1.25` (hardware `delay/T2_echo`; fixture
`delay/T2_parameter`). Recomputed from the saved `populations.rho_01_abs` rows, their
late-time means are respectively 0.072615 and 0.018520. This is a cross-fixture
comparison, not a sigma-significance statement or a matched causal control.

The retained null description and the separate fixture answer narrower questions. The recorded ensemble reports separation from its exponential-decay, binomial-shot, one-phase-per-run null; the separate fixture has larger late-time coherence. Neither comparison identifies the Q52 mechanism or calibrates the probability of the 17/17 directional record under the relevant hardware alternatives.

## Three Hypotheses

### H1: Systematic SPAM Error (Most Conservative)

A constant angular calibration error could dominate the estimate when the true signal is small and could produce directional consistency. This model was proposed, not fitted.

**Open point:** A constant-only offset was not fitted jointly with drift and tomography-systematic terms. The selected positive tail slope does not by itself exclude SPAM.

**Test:** Measure with \|+⟩ and \|−⟩ initial states. Under the proposed constant-offset model the SPAM term does not flip with the preparation, whereas the prepared coherence does.

### H2: Two-Level System Coupling (Known Physics)

A substrate TLS coupled to qubit 52 could feed coherence back into the system. This is a candidate mechanism, not an identification from the saved record.

**Scope:** Hours-to-days drift would ordinarily look nearly static across the sampled `0-894.742 us` evolution-time window. The acquisition order and wall-clock duration are not recorded here, so neither can be inferred from the maximum delay. The fixed direction is not surprising and is not diagnostic of TLS. A direction change days later would likewise not identify TLS, because detuning and calibration drift can also rotate the phase.

**Control:** Repeat across acquisition times and nearby qubits, log calibration drift, and compare an explicit TLS model against detuning/SPAM alternatives on held-out complex data. Persistence, drift, or locality alone would narrow candidates but would not identify one.

### H3: External Coherent Coupling (Speculative)

An external coherent source was considered as a speculative way to generate phase-directed coherence. No such model was fitted to this record.

**Scope:** This record identified and fitted no external coherent source or coupling path. Microwave leakage/crosstalk and chip or package modes are mundane coherent alternatives, and refrigerator temperature alone does not exclude a coherent drive; the saved record does not distinguish them.

**Test:** Multi-qubit correlation. A repeatable cross-qubit signature would narrow local TLS and individual-SPAM explanations, but would still require explicit controls before a mechanism claim.

## March 2026 Hardware Test Plan

Based on these findings, the March run should include:

1. **Reproduce:** Same qubit 52, same protocol. Does the 17/17 directional signature persist?
2. **SPAM discrimination:** Run with \|+⟩ and \|−⟩ initial states. Does the excess flip sign?
3. **Multi-qubit:** At least 5 qubits (ideally 10+). Check Re/Im sign consistency independently for each.
4. **Cross-correlation:** Do different qubits show correlated excess coherence?
5. **Extended time range:** Push to 5×T2_echo or beyond. Does the rising trend continue?

If a fitted SPAM model accounts for the record: document the calibration explanation.
If SPAM is disfavored and a TLS model survives its controls: study that conventional non-Markovian candidate.
If both are disfavored: retain the mechanism as an open question rather than assigning a cause.

## Raw Numbers for Reference

Late-time off-diagonal elements (`t/T2_echo >= 1.0`):

```
t/T2_echo   Re(ρ₀₁)      Im(ρ₀₁)      |ρ₀₁|
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

This analysis began as a search for structure after the measured C·Ψ = ¼ crossing. The scalar crossing does not identify a quantum-to-classical transition or a source for the residual. Because `|ρ₀₁|` also enters C·Ψ, proximity-to-boundary and residual-amplitude correlations from this same record are not independent boundary evidence.

## Epistemic Status

- **Recorded:** Zero exceedances in 10,000 draws of the named null; broader hardware alternatives were not tested
- **Recorded:** Directional consistency Re+/Im- in 17/17 sampled late points
- **Measured on the selected 13-row tail:** Positive slope, p = 0.0531, interval crossing zero; non-monotone and cut-sensitive
- **Unknown:** Whether this is SPAM, TLS, or something else
- **Not claimed:** A non-Markovian witness, a universal boundary mechanism, or a specific Q52 cause

The correct response to this data is not belief or disbelief. It is: measure again.
