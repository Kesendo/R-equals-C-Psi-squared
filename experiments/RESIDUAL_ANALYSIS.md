# Residual Analysis: Anomalous Late-Time Coherence in IBM Torino Data

## Discovery Date
2026-02-09 (same-day reanalysis of tomography data)

## Status
**Exploratory.** Anomaly detected and quantified. Cause unknown. Three hypotheses proposed, none confirmed. March 2026 hardware run designed to discriminate.

## Context

After completing the IBM Quantum Tomography experiment (see IBM_QUANTUM_TOMOGRAPHY.md), the raw density matrices were reanalysed for structure in the residuals. The original experiment measured C·Ψ crossing the ¼ boundary and confirmed the generalized crossing equation. This analysis asks a different question: **Is there structure in the data beyond the expected exponential decay?**

## Finding 1: Excess Late-Time Coherence (p < 0.0001)

### Null Hypothesis Test

A Monte Carlo simulation was constructed:
- 10,000 synthetic tomography experiments
- Same parameters: T₁ = 221.2 μs, T₂* = 110.2 μs (fitted), 8192 shots, 25 delay points
- Shot noise modelled as binomial sampling across 3 tomography bases (X, Y, Z)
- Random phase per run (detuning)

**Result:** The real IBM Torino data has late-time coherence (t/T₂ > 1.25) that exceeds all 10,000 simulated experiments.

| Metric | IBM Torino (real) | Simulator (mean ± std) |
|--------|------------------|----------------------|
| Mean \|ρ₀₁\| for t > 370 μs | 0.01852 | 0.00861 ± 0.00105 |
| p-value | < 0.0001 | - |
| Points above 99th percentile | 9 of 25 | 0.25 expected |

The strongest anomaly is at t/T₂ = 2.625, where \|ρ₀₁\| = 0.0363, which is 2.2× above the 99th percentile of the null distribution (0.0166).

### Important Caveat

The null model is simple: exponential decay + shot noise. It does not include TLS coupling, 1/f noise, readout error, or non-Markovian dynamics. Any of these could close the gap. The point is not that the effect is unexplainable, but that it is **quantifiably present** and deserves investigation.

## Finding 2: Directional Coherence (17/17 Sign Consistency)

This is the most striking anomaly.

For all 17 data points at t/T₂ ≥ 1.0:
- **Re(ρ₀₁) > 0 in 17/17 measurements**
- **Im(ρ₀₁) < 0 in 17/17 measurements**

The residual coherence always points into the fourth quadrant of the complex plane.

| Component | Mean | Std | All same sign? |
|-----------|------|-----|---------------|
| Re(ρ₀₁) | +0.01137 | 0.00499 | Yes (17/17 positive) |
| Im(ρ₀₁) | -0.01279 | 0.00805 | Yes (17/17 negative) |

Probability of 17/17 same sign in both components by chance: (1/2)^34 ≈ 6 × 10⁻¹¹.

This is remarkable because quantum decoherence should randomize the phase of ρ₀₁. At t = 3×T₂, the coherence should be indistinguishable from zero with random phase. Instead, it has a fixed direction.

## Finding 3: Rising Coherence Trend

A linear fit to \|ρ₀₁\| for t/T₂ > 1.5 shows a **positive slope**:

```
Slope = +0.00819 per T₂ unit
```

Coherence that increases with time violates the expectation for any Markovian open quantum system. Under Lindblad dynamics, off-diagonal elements of the density matrix can only decay (modulo oscillations). A net positive trend over 13 consecutive points spanning 1.5 T₂ is anomalous.

## Finding 4: Qubit Detuning and Phase Structure

A frequency detuning of ω = -0.0357 rad/μs (f = -5.68 kHz) was extracted from early high-SNR phase data. This corresponds to a period of ~176 μs.

After derotation, significant residual phase remains at intermediate times (t/T₂ = 0.375 to 1.375), indicating the detuning is not perfectly linear or the system has additional frequency components.

Late-time phases quantised to 30° bins show strong clustering:
- -60°: 9 of 15 values (60%, expected ~17%)
- -30°: 4 of 15
- 0°: 1 of 15
- -90°: 1 of 15

## Finding 5: Revival Peak Spacing

Local maxima in late-time \|ρ₀₁\| occur at:

| Peak | t/T₂ | \|ρ₀₁\| |
|------|-------|---------|
| 1 | 1.500 | 0.01964 |
| 2 | 1.750 | 0.01234 |
| 3 | 2.000 | 0.02210 |
| 4 | 2.625 | 0.03629 |

Peak spacings: 0.250×T₂, 0.250×T₂, 0.625×T₂. The first two gaps are identical (T₂/4). The ratio of the third gap to the first is exactly 2.5.

## Qiskit Simulator Comparison

The same experiment was run on Qiskit Aer simulator with a noise model. Key comparison:

The simulator also shows residual late-time coherence (mean 11.4σ above shot noise). However, the simulator's late-time coherence is **higher** than the real hardware (0.073 vs 0.019). This initially suggested the real-hardware anomaly was an artifact.

The Monte Carlo null hypothesis test resolved this: the Qiskit Aer simulator includes a more complex noise model than pure exponential decay, which inflates late-time coherence. Our custom null model (exponential + shot noise only) correctly establishes that the **real hardware exceeds what the simple physical model predicts**, and neither the simple model nor the Aer model produces the 17/17 directional consistency.

## Three Hypotheses

### H1: Systematic SPAM Error (Most Conservative)

Tomography gate calibration has a small, constant angular error. When the true signal is zero, the fixed offset dominates. This would explain the directional consistency.

**Problem:** A constant offset does not explain the rising trend. The excess grows from 0.003 to 0.036 over the measurement range, spanning a full order of magnitude. A fixed gate error produces a fixed offset.

**Test:** Measure with \|+⟩ and \|−⟩ initial states. SPAM offset is independent of initial state. Physical coherence flips sign.

### H2: Two-Level System Coupling (Known Physics)

A TLS defect in the substrate is coupled to qubit 52 and feeds coherence back into the system. TLS coupling is well-documented in superconducting qubits and can produce non-Markovian coherence revivals with a preferred phase.

**Problem:** TLS coupling strength and frequency drift on timescales of hours to days. The fixed directionality across 900 μs of evolution is surprisingly clean for a TLS interaction.

**Test:** Repeat measurement on qubit 52 days later. If the direction changes, it was TLS. Also measure neighboring qubits: TLS effects are local to individual qubits.

### H3: External Coherent Coupling (Speculative)

An external source is driving coherence into the qubit at a fixed phase. This would explain all three observations: excess coherence, fixed direction, and rising trend (coupling accumulates over time).

**Problem:** No known mechanism for coherent coupling to a single qubit in a dilution refrigerator at 15 mK. This hypothesis is physically extraordinary and requires extraordinary evidence.

**Test:** Multi-qubit correlation. If independent qubits on the same chip show correlated excess coherence with the same directional signature, it cannot be local TLS or individual SPAM. This would be genuinely unexplained.

## March 2026 Hardware Test Plan

Based on these findings, the March run should include:

1. **Reproduce:** Same qubit 52, same protocol. Does the 17/17 directional signature persist?
2. **SPAM discrimination:** Run with \|+⟩ and \|−⟩ initial states. Does the excess flip sign?
3. **Multi-qubit:** At least 5 qubits (ideally 10+). Check Re/Im sign consistency independently for each.
4. **Cross-correlation:** Do different qubits show correlated excess coherence?
5. **Extended time range:** Push to 5×T₂ or beyond. Does the rising trend continue?

If H1 (SPAM) survives: The anomaly is a calibration artifact. Document and move on.
If H1 is killed and H2 survives: Non-Markovian dynamics worth studying but conventional.
If both H1 and H2 are killed: We have a genuine open question.

## Raw Numbers for Reference

Late-time off-diagonal elements (t/T₂ ≥ 1.0):

```
t/T₂   Re(ρ₀₁)      Im(ρ₀₁)      |ρ₀₁|
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

This analysis began as a search for structure beyond the confirmed C·Ψ = ¼ crossing. The excess coherence exists in a regime where R = CΨ² predicts the system has crossed the boundary from quantum to classical. If the residual coherence is physical (not SPAM), it represents structure persisting beyond the theoretical decoherence boundary, which has direct implications for the framework's treatment of what happens at and beyond the ¼ threshold.

## Epistemic Status

- **Confirmed:** Excess coherence exists beyond null model (p < 0.0001)
- **Confirmed:** Directional consistency Re+/Im- in 17/17 points
- **Confirmed:** Rising trend in late-time coherence
- **Unknown:** Whether this is SPAM, TLS, or something else
- **Not claimed:** Any specific interpretation beyond "anomaly worth investigating"

The correct response to this data is not belief or disbelief. It is: measure again.
