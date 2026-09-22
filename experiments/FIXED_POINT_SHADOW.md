<!-- QUARTER-CURRENT -->
# Q52 Residual Record and the Failed Universal Shadow Interpretation

<!-- Keywords: IBM Torino Q52 late-time coherence, finite residual record,
Q80 Q102 cross-qubit comparison, phase-compatible fit, open hardware
mechanism, C*Psi quarter interpretation audit -->

**Status:** Finite Q52 record retained; universal-boundary/non-Markovian-witness interpretation closed; Q52 late-time excess mechanism open
**Date:** 2026-02-09 record; 2026-03-09 cross-qubit comparison
**Scope:** Measurement record and hypothesis audit. No typed Claim, Witness, or hardware Confirmation owns the Q52 late-time residual/magnitude-excess mechanism.
**Related:** [IBM Quantum Tomography](IBM_QUANTUM_TOMOGRAPHY.md), [Finite Residual Analysis](RESIDUAL_ANALYSIS.md), [Predictions and Interpretive Questions](../docs/PREDICTIONS.md)

## Current result

The Q52 tomography remains a finite hardware record with directional late-time coherence and an excess over one narrow null model. It does not establish a physical fixed-point remnant, a special boundary process, or a cause for the excess magnitude.

Only the universal-boundary/non-Markovian-witness interpretation is closed.

Detuning is the preferred explanation for the phase component.

The Q52 late-time excess mechanism remains unresolved absent a Q52-specific fit/control.

The March Q80/Q102 comparison is load-bearing only against a universal reading. The Q80 fit is compatible with detuning, but it neither fits Q52 nor supplies an independent causal diagnosis.

## The Q52 record

The saved record contains 25 tomography points at 8,192 shots. Its `delay_over_T2` coordinate uses the IBM calibration value `T2_echo = 298.247 us`; it is not the fitted free-induction value `T2* = 110.7 us`. The late slope is reported only on that echo-normalized axis.

At or beyond t/T2_echo = 1.0, all 17 sampled points have Re(rho_01) > 0 and Im(rho_01) < 0.

For the 13 rows selected by t/T2_echo >= 1.5, a least-squares line through |rho_01| has slope +0.00819 per T2_echo.

For that 13-row fit, the correlation coefficient is 0.5469, the two-sided p-value is 0.0531, and the 95% slope interval `[-0.00013, +0.01651]` crosses zero. Changing the cut to strict `t/T2_echo > 1.5` leaves 12 rows and changes the slope to `+0.01041`; the tail statement is selection-sensitive.

Those 13 amplitudes are non-monotone, so the finite positive slope is not evidence of a growing asymptote.

The corresponding slope is `2.7457e-5/us` on the echo-normalized axis. Exponential damping also does not switch off after a fixed number of time constants.

Because |rho_01| also enters C*Psi, the reported r = -0.9955 correlation with distance from 1/4 is descriptive, not independent boundary evidence.

In the saved analysis, `Psi = 2|rho_01|` and `C*Psi = C 2|rho_01|` exactly. Correlating `|rho_01|` with `0.25 - C*Psi` therefore reuses the same measured quantity on both axes.

The saved analysis reports zero exceedances in 10,000 draws of a model described as exponential decay, binomial shot sampling, and one random phase per synthetic run.

That description is not a comparison against a Q52-fitted time-dependent detuning or drift, SPAM, TLS, or memory model. This document and its retained provenance links do not identify the producer for that null calculation, so the recorded `p < 0.0001` cannot carry a more specific null than the retained method description. The nominal `4^-17 = 5.82e-11` sign probability likewise assumes independent, uniformly distributed phases and a pre-specified quadrant. It is not a mechanism probability.

A comparison between the phase of an algebraic root `R-` and the phase of `rho_01` is not a dynamical result. The roots belong to the chosen scalar recurrence; no transport law in this record makes `rho_01` track one of them or preserve its approach direction.

The separate Q52 quarter crossing remains a qualitative reconstructed-trajectory record. Its precision and same-record prediction scope are owned by the cockpit validation, not by the residual analysis here.

## What the March comparison establishes

The March run retained ten tomography rows for each of Q80 and Q102: two reference rows and eight sampled late rows.

Q80 has eight of eight sampled late points in quadrant 1; Q102 has mixed signs and quadrants across its eight sampled late points.

The cross-qubit comparison rejects the former universal-boundary reading; it does not identify the Q52 excess mechanism.

Q102 supplies no consistent direction in its eight sampled late points; that observation does not identify the cause as shot noise.

Q102 also has eight negative real residuals and seven stored significant-excess flags, so a mixed phase pattern cannot be promoted to a complete null diagnosis. Q80's phases are directional but not monotone. No retained producer supports a fast-rotator/slow-drifter or shared hardware-skeleton mechanism.

The checked simulator record reports two simulated qubits with no directional match. The cockpit has a matched simulator record for Q80 but none for Q102. These are finite controls, not a universal simulator theorem.

## The Q80 fit, scoped

All Q80 delays lie on a `5.412785 us` grid, so complex samples identify frequency only modulo `184.747761 kHz`.

The quoted Q80 frequencies are near-zero representatives of those alias classes, not identified absolute detunings.

The exploratory Q80 script fits the fixed-T2 phase-line model

`rho_fit(t) = A0 exp(-t/T2) exp(i(phi_0 + m t))`, with `A0 = |rho_01(0)|`.

Under the Hamiltonian convention `rho_01(t) proportional to exp(-i delta_omega t)`, the physical detuning parameter is `delta_omega = -m`. The selected near-zero phase-slope representative is +1.27 kHz; under the stated `exp(-i delta_omega t)` convention its detuning representative is `delta_f = -1.27 kHz`.

It reports two different same-record comparisons:

| Q80 comparison | Fitted quantities | Mean complex error | Nested comparator | Scope |
|---|---:|---:|---:|---|
| Hahn-T2 fixed phase line | slope and intercept on eight late rows | 0.0356 | intercept-only 0.0508 | both scored on nine nonzero-time rows |
| Free complex fit | `T_eff`, `delta_omega` on nine rows | 0.0138 | fitted no-detuning envelope 0.0487 | both fitted and scored on the same nine rows |

The fixed-T2 phase-line curve has mean complex error 0.0356 versus 0.0508 for the intercept-only comparator, a 1.4x in-sample error ratio. Both curves use the same `A0`, Hahn `T2`, phase-fit rows, and scoring rows; the line adds only the fitted slope.

The free Q80 fit gives T_eff = 23.25 us and the near-zero alias representative delta_f = -2.58 kHz, with mean complex error 0.0138 versus 0.0487 for its fitted no-detuning envelope comparator, a 3.5x in-sample error ratio.

The fitted Q80 models are same-record, in-sample comparisons; neither is a Q52 fit or a held-out prediction.

An ideal static Z detuning rotates `rho_01` but leaves `|rho_01|` unchanged. It can explain a phase component without, by itself, explaining the Q52 magnitude excess. The free-fit `T_eff = 23.25 us` lies below the March 9 Hahn value `T2_echo = 27.06 us`. The record's `10.83 us` field is `T2_echo/2.5`, a scheduling proxy, not a measured Ramsey time; the March 18 Ramsey `T2* = 17.36 us` belongs to a different run and was not an input to this fit. None of these time comparisons identifies a microscopic cause.

## Ownership and open status

The cross-term formula does not identify a Q52 hardware mechanism.

Its proved scope concerns dephasing cross terms for named multi-site Hamiltonian couplings. A local single-qubit `delta_omega Z/2` term commuting with Z dephasing is compatible algebra, but that observation is not a Q52 parameter estimate or causal test.

Detuning is the preferred explanation for the phase component.

The Q52 late-time excess mechanism remains unresolved absent a Q52-specific fit/control.

OQ-033 and OQ-098 remain open.

The typed hardware registry owns the qualitative Q52 crossing and the same-record absorption-ratio comparison. It does not own a residual-coherence mechanism. No live Witness or OpenArc closes that gap.

This record does not establish:

- a physical remnant of the scalar fixed point;
- a special quarter-boundary source or universal scar;
- a demonstrated non-Markovian revival or TLS cause;
- a non-decaying or growing late-time asymptote;
- a Q52 causal diagnosis imported from the Q80 fit;
- a shared three-qubit-simulation and one-qubit-hardware mechanism.

## Discriminating controls

A Q52 mechanism claim needs a Q52 experiment designed to distinguish the alternatives, not another fit to Q80. A useful minimum set is:

1. repeat the complex Q52 trajectory with acquisition order randomized and calibration drift logged;
2. measure same-session Ramsey and Hahn-echo twins so detuning and envelope decay are independently constrained;
3. sweep a known positive and negative detuning and test the predicted phase reversal while checking whether the magnitude excess changes;
4. alternate `|+>` and `|->` preparations to separate a physical sign reversal from a fixed SPAM offset;
5. pre-register fit and hold-out times, then compare detuning, offset/SPAM, colored-noise, and TLS-compatible models on held-out complex data;
6. repeat across time and neighboring qubits to distinguish persistent calibration bias from local drift or a moving defect.

Until such a control exists, the correct result is the finite record plus an open mechanism.

## Reproducibility

- Q52 raw tomography: [Q52 tomography JSON](../data/ibm_tomography_feb2026/tomography_ibm_torino_20260209_131521.json)
- March Q80/Q102 record: [March combined JSON](../data/ibm_shadow_march2026/shadow_hardware_combined_20260309_181852.json)
- March simulator record: [`shadow_simulate_20260309_181709.json`](../data/ibm_shadow_march2026/shadow_simulate_20260309_181709.json)
- Q80-only exploratory fit: [Q80 in-sample fit script](../simulations/shadow_ibm_retrodict.py)
- Current crossing and simulator scope: [cockpit output](../simulations/results/cockpit_validation.txt)
- Historical residual calculations and hypotheses: [Finite Residual Analysis](RESIDUAL_ANALYSIS.md)
- Current summary owner: [Predictions and Interpretive Questions](../docs/PREDICTIONS.md)
- Open questions: [OQ-033 and OQ-098](../review/OPEN_QUESTIONS_INDEX.md)
