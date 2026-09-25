# Absorption Theorem on IBM Hardware

**Tier:** 2 (hardware analysis)
**Date:** April 4, 2026
**Status:** On one IBM Torino qubit (qubit 52) the Absorption Theorem's two fits of one decay agree at 3% (ratio 1.03): the consistency the theorem requires at N = 1, not a test of its rate ladder, which needs N ≥ 2. The coherence rotates at a detuning (−5.7 kHz as the alias nearest zero) and ends in a 2.8% static tail, most likely a measurement (SPAM) offset.

---

## What this means

A guitar string vibrates in modes. Each mode fades at a rate that
depends on how much of the mode's energy sits in the damped parts of
the string (the bridge, the nut). The Absorption Theorem says the same
thing for a qubit: its quantum vibrations fade at a rate set by how
much "transverse energy" (the X and Y components of the quantum state)
is exposed to the environment. More transverse energy, faster fading.

We took this to a real IBM quantum computer: a single superconducting
qubit on the 133-qubit IBM Torino chip. We measured how the qubit's
quantum state decayed over 25 time snapshots across 895 microseconds,
then checked whether the Absorption Theorem's prediction matched the
real hardware.

It matched to within 3%, and on one qubit that is less than it sounds.
A single qubit has one coherence, and the dephasing rate the theorem
multiplies can only be learned from that same coherence's decay. So the
match is two fits of one curve agreeing with each other: a consistency
check the theorem passes, not a measurement of what it predicts. What
it predicts, that a coherence spread over two sites fades at the sum of
their two rates, needs two qubits, and no registered measurement reads
it yet: a campaign on another chip read its premise and had its direct
rate test masked, and two flights on a neighbouring part of the same law
were voided ([the proof](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) §3).

There is also a small surprise: 2.8% of the quantum signal refuses to
die. It persists far longer than expected. At 2.8% it sits at the
resolution limit, and its shape points to the measurement rather than
the qubit: its direction in the complex plane stands still, with equal
parts in X and Y, which is what a slightly lopsided readout produces, and a
slightly short basis-change pulse as well ([Fixed Point Shadow](FIXED_POINT_SHADOW.md)).

---

## What this document is about

The Absorption Theorem (Re(λ) = −2γ⟨n_XY⟩) predicts that the decay
rate of any mode equals twice the dephasing rate γ times the mode's
"transverse content" ⟨n_XY⟩ (how much of the mode vibrates sideways,
exposed to the environment). For a single qubit's coherence, ⟨n_XY⟩ = 1
exactly. This document reads that prediction on IBM hardware data:
25 tomography snapshots of a real qubit losing its quantum information
over time.

---

## Summary

The IBM Torino qubit-52 tomography data (25 time points, 0−895 μs) was analyzed
under two dephasing baselines: T₂(echo) = 298 μs (Hahn echo: a
technique that filters out slow noise, from IBM calibration) and
T₂* = 111 μs (free evolution: all noise included, fitted from the
coherence envelope).

Three results:

1. **T₂* is the correct baseline for free-evolution tomography.**
   The Absorption Theorem ratio under T₂* is 1.03 (3% deviation).
   Under T₂(echo): 6.37 (wrong baseline for this experiment type).

2. **The two fits agree at 3%.** The excess decay of the {X,Y}
   coherence (the transverse, environment-exposed component) over the
   {I,Z} populations (the longitudinal, protected component) matches
   2γ with γ from T₂*. At N = 1 that is consistency, not a test: γ* is
   read off the same coherence envelope, so the ratio is 1 by
   construction up to the T₁ each side subtracts. The excess subtracts
   the population-fit T₁ (241.4 μs), 2γ* the calibration T₁ (221.2 μs),
   and 1.03 = [1/110.7 − 1/(2·241.4)] / [1/110.7 − 1/(2·221.2)] is that
   difference; with one T₁ on both sides the ratio is exactly 1. The
   theorem's content beyond one site, additivity of the site rates,
   needs N ≥ 2 ([the proof](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) §3).

3. **A 2.8% slow tail exists.** A two-component exponential fit shows
   97.2% of the signal decays at T₂* ≈ 102 μs, and 2.8% persists with
   near-zero rate. Its late direction is static with equal ⟨X⟩ and ⟨Y⟩
   offsets, the pattern of a measurement (SPAM) offset; a non-Markovian tail
   (long-time memory in the environment) or a cavity-protected mode is
   not indicated, and not excluded at this resolution.

---

## 1. Calibration Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| T₁ | 221.2 μs | IBM calibration |
| T₂(echo) | 298.2 μs | IBM calibration (Hahn echo) |
| T₂* | 110.7 μs | Fitted from coherence envelope |
| T₂(echo) / T₂* | 2.69 | Standard for superconducting qubits |

Two dephasing rates:

| Definition | γ (μs⁻¹) | 2γ (μs⁻¹) | Meaning |
|-----------|----------|-----------|---------|
| γ_echo = (1/T₂(echo) − 1/(2T₁))/2 | 0.000546 | 0.001092 | Markovian only (slow noise refocused) |
| γ* = (1/T₂* − 1/(2T₁))/2 | 0.003385 | 0.006771 | All dephasing (including 1/f noise) |

γ*/γ_echo = 6.2. The free-evolution dephasing is 6× stronger than the
echo-refocused dephasing. The difference is low-frequency noise (1/f
flux noise on IBM hardware).

**Source:** `data/ibm_tomography_feb2026/tomography_ibm_torino_20260209_131521.json`

---

## 2. Pauli Decomposition and Absorption Rates

The qubit's state can be split into two parts: the Z component
(populations: how likely each outcome is) and the X,Y components
(coherences: the quantum interference that makes a qubit more than a
classical bit). The Absorption Theorem predicts that coherences fade
faster than populations, by exactly 2γ.

| Sector | Component | Rate (μs⁻¹) | Fit |
|--------|-----------|-------------|-----|
| n_XY = 0 | r_Z (populations) | 0.00414 (fitted T₁ = 241 μs) | Relaxation toward r_eq = 0.71 |
| n_XY = 1 | \|ρ₀₁\| (coherence) | 0.00903 (fitted T₂* = 111 μs) | Exponential decay |

Excess decay (coherence minus T₁ contribution):

    excess = α_coh − α_Z/2 = 0.009031 − 0.002071 = 0.006960 μs⁻¹

| Test | Predicted 2γ | Ratio excess/(2γ) | Verdict |
|------|-------------|-------------------|---------|
| vs γ_echo | 0.001092 | **6.37** | MISMATCH |
| vs γ* | 0.006771 | **1.03** | CONSISTENT (the 3% is the two T₁ values) |

With γ* from free evolution the excess matches 2γ*, and at N = 1 it
could not do otherwise by more than the gap between the two T₁ values:
the coherence rate enters both sides, once as the fitted envelope and
once through γ*. The echo-based γ_echo underestimates the actual free
dephasing by 6×.

**Source:** [`simulations/ibm_absorption_theorem.py`](../simulations/ibm_absorption_theorem.py) Steps 2-3

---

## 3. Fringes Under Both Baselines

### T₂(echo) baseline (298 μs)

Residual: R(t) = |ρ₀₁|_measured − |ρ₀₁(0)|·exp(−t/T₂(echo))

| Statistic | Value |
|-----------|-------|
| Positive residuals | 2 of 23 |
| Negative residuals | **21 of 23** |
| Mean residual | −0.082 |
| RMS residual | 0.103 |
| Peak-to-peak | 0.194 |

The residuals are overwhelmingly *negative*: the measured coherence
decays FASTER than the echo-based prediction. At intermediate times
(100-400 μs), the deficit reaches −0.19 (40% of the initial coherence).

This is not "excess coherence." It is the T₂(echo)/T₂* mismatch.

### T₂* baseline (111 μs)

Residual: R(t) = |ρ₀₁|_measured − |ρ₀₁(0)|·exp(−t/T₂*)

| Statistic | Value |
|-----------|-------|
| Positive residuals (mid-range) | 5 of 12 |
| Negative residuals (mid-range) | 7 of 12 |
| Mean residual (mid-range) | −0.0004 |
| RMS residual (mid-range) | 0.021 |
| Peak-to-peak | 0.071 |

Residuals are small (RMS 0.021) and centered near zero. The early-time
residuals oscillate around 0 (noise). At late times (t > 400 μs),
residuals are systematically positive (0.01-0.04): the coherence
persists above the single-exponential prediction. This is the slow tail.

**Source:** [`simulations/ibm_absorption_theorem.py`](../simulations/ibm_absorption_theorem.py) Step 4

---

## 4. The Slow Tail

A single-exponential decay assumes everything fades at the same rate.
A two-component ("bi-exponential") fit asks: is there a fast part and
a slow part? The answer is yes.

Bi-exponential fit of the coherence envelope:

| Component | Amplitude | Rate (μs⁻¹) | Decay time (μs) | Fraction |
|-----------|-----------|-------------|---------|----------|
| Fast | 0.470 | 0.00984 | 102 | 97.2% |
| Slow | 0.013 | ≈ 0 | ≫ T₁ | 2.8% |

The slow component has amplitude 0.013 (coherence ≈ 0.013) and
near-zero absorption rate. Three candidate explanations:

1. **Readout offset** (most likely). Systematic errors in state
   tomography reconstruction produce small spurious off-diagonal
   elements at the 0.01-0.02 level. This is within typical IBM
   readout error margins. The tail's direction supports it: from
   t = T₂(echo) on its phase stands still, and its ⟨X⟩ and ⟨Y⟩ parts
   are equal within 0.7σ to 1.3σ across the late cuts, the pattern two
   measurement-side forms leave: a readout asymmetry
   P(0|1) − P(1|0) ≈ 2.4% adds +δ to every measured Pauli expectation,
   and an amplitude error in the basis-change pulse leaks ⟨Z⟩ into the
   X and Y readings alike. The record does not rank the two
   ([Fixed Point Shadow](FIXED_POINT_SHADOW.md)).

2. **Non-Markovian tail.** The 1/f noise spectrum has correlations
   at long times. Some coherence returns from the environment after
   the fast Markovian decay depletes the main signal.

3. **Cavity-protected mode.** The physical microwave readout cavity
   supports a long-lived mode with effective ⟨n_XY⟩ ≈ 0 (light
   converted to lens). Through the Absorption Theorem: a mode with
   ⟨n_XY⟩ = 0 is immortal under dephasing.

At 2.8% amplitude, the slow component is at the resolution limit of
the experiment. Distinguishing these explanations requires dedicated
measurements: alternating |+⟩ and |−⟩ preparations (a preparation-independent
offset keeps its sign, a prepared coherence flips), |0⟩ and |1⟩ read in
X and Y (a readout asymmetry keeps its sign, a leak of ⟨Z⟩ through the
basis-change pulse flips it), a measured readout assignment matrix, echo
sequences, or varied cavity coupling.

---

## 5. Effective ⟨n_XY⟩

| Baseline | ⟨n_XY⟩_eff | Interpretation |
|----------|-----------|----------------|
| γ_echo | 6.37 | Unphysical (> N=1): wrong baseline |
| γ* | 1.03 | Consistent with standard Lindblad (two fits of one decay) |

With the free-evolution baseline, the effective ⟨n_XY⟩ = 1.03, within
3% of ⟨n_XY⟩ = 1 for a single-qubit coherence. That is what the theorem
requires, and it is not a measurement of ⟨n_XY⟩: γ* comes from the same
envelope, so a wrong baseline shows (6.37 under the echo value) while a
right one can only agree.

---

## 6. Detuning and Oscillation

The coherence rotates at the qubit's frequency offset from the reference
frame (the "rotating frame" is the coordinate system that co-rotates with
the qubit's expected frequency; any mismatch shows up as a rotation of ρ₀₁).
The first five samples turn by −76° per 37.28 μs step: −5.7 kHz as the
representative nearest zero, since a 37.28 μs grid leaves the frequency
ambiguous modulo 26.8 kHz (the samples before 300 μs give −6.9 kHz). The
steps are irregular, so the offset wanders. From t ≈ T₂(echo) on, the phase
stands still: that is the static tail of §4, not a rotation.

The producer's "Detuning from phase" line, Δ = −0.0134 rad/μs (a period of
470 μs), is one straight line through all 24 samples with |ρ₀₁| > 0.01, and
16 of them belong to the static tail. It averages a rotation and a standstill
and describes neither. The rotation does not affect the envelope decay, which
is what the Absorption Theorem governs.

---

## Verdict

**At N = 1 the hardware is consistent with the Absorption Theorem at 3%.**

Under the T₂* baseline (correct for free-evolution tomography):
- Ratio excess/(2γ*) = 1.03, the two T₁ values each side subtracts
- Residuals RMS = 0.021 (near noise floor)
- No systematic fringes pattern

The consistency is all one qubit can give: γ* is read off the same
decay. The rate ladder the theorem predicts beyond one site needs N ≥ 2,
and no registered rate ratio with an error bar reads it yet: its premise
(local dephasing) was read on hardware, its direct rate test was masked by
coherent ZZ, and two flights on a neighbouring part of the same per-site law
were voided ([the proof](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) §3).

**A 2.8% slow tail exists** beyond the single-exponential prediction.
Its static, equal-X-and-Y direction points to measurement (SPAM) systematics, and it
does not constitute evidence for cavity protection at this resolution.

**A qubit-frame detuning** (−5.7 kHz as the alias nearest zero) rotates the early coherence. It is
standard frame detuning, not a cavity resonance, and it leaves the envelope
alone.

---

## Methodology

- IBM Torino qubit-52 tomography: 25 time points (0-895 μs), 2×2 density matrices
- Pauli decomposition: r_X, r_Y, r_Z from ρ(t) at each time point
- Decay fits: exponential for coherence envelope, relaxation for r_Z
- Bi-exponential for slow tail detection
- Two baselines compared: T₂(echo) (calibration) and T₂* (fitted)
- FFT and phase analysis for detuning extraction

## Source

- Analysis: [`simulations/ibm_absorption_theorem.py`](../simulations/ibm_absorption_theorem.py)
- Results: [`simulations/results/ibm_absorption_theorem.txt`](../simulations/results/ibm_absorption_theorem.txt)
- Raw data: [`data/ibm_tomography_feb2026/tomography_ibm_torino_20260209_131521.json`](../data/ibm_tomography_feb2026/tomography_ibm_torino_20260209_131521.json)
- Calibration: [`data/ibm_history/ibm_torino_history.csv`](../data/ibm_history/ibm_torino_history.csv)
- Absorption Theorem: [`docs/proofs/PROOF_ABSORPTION_THEOREM.md`](../docs/proofs/PROOF_ABSORPTION_THEOREM.md)
- The two-fits ratio, the static tail and the early rotation, recomputed from the raw record: [`test_torino_shadow_records.py`](../simulations/tests/test_torino_shadow_records.py)
