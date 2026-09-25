# Fixed Point Shadow: A Static Offset, Not a Boundary Effect

<!-- Keywords: fixed point shadow residual coherence, IBM Torino Q52 Q80 Q102,
CΨ quarter boundary crossing late-time coherence, phase at the crossing,
qubit detuning Ramsey fringe frequency, static measurement offset SPAM readout
asymmetry basis pulse leak, late window scheduling T2 echo versus own decay,
pre-registered cross-qubit test negative result, Monte Carlo null late-time
excess, R=CPsi2 fixed point shadow -->

**Status:** Tier 2. A negative result: the shadow is not a property of the ¼ boundary. Qubit 52's own record refutes the proposal as it was made: at the crossing its coherence pointed elsewhere and kept turning, and its late direction is not the fixed point's (−48° against −24°) but a static offset with the pattern of a measurement (SPAM) offset, mechanism open. The March test on two other qubits measured their detuning and could not see a residual of this size, so the question in its cross-qubit form is still to be run.
**Date:** 2026-02-09 (the Q52 record), 2026-03-09 (the March test on Q80 and Q102)
**Data:** `data/ibm_tomography_feb2026/`, `data/ibm_shadow_march2026/`, and the Ramsey fits in `data/ibm_run3_march2026/`
**Depends on:** [Boundary Navigation](BOUNDARY_NAVIGATION.md), [IBM Quantum Tomography](IBM_QUANTUM_TOMOGRAPHY.md), [Residual Analysis](RESIDUAL_ANALYSIS.md)

---

## What this document is about

After a qubit crossed the CΨ = ¼ boundary on IBM hardware, a small residual
coherence was found pointing in one fixed direction of the complex plane, in
17 of 17 late samples. We proposed that it was the last bearing of the complex
fixed point the qubit had left behind at the crossing: a shadow of the fixed
point. Before the next hardware run we wrote down one question, whether the
shadow belongs to qubit 52 or to the ¼ boundary, and built a test on two other
qubits. That test measured something else: each qubit's coherence still turning
at its own frequency offset. Its late window, timed on the calibration's echo
time, closed while the prepared coherence was still twice the size of the
residual it was looking for. The answer to our proposal was in qubit 52's own
record. At the crossing its coherence pointed another way and kept turning for
another 150 μs, and the direction that stands still late, −48°, is not the one
the fixed point pointed, −24°. What stands there is a static offset whose X and
Y parts are equal, the pattern of a measurement (SPAM) offset. The ¼ crossing
leaves no scar in the direction the fixed point pointed. The question as we
asked it, about every qubit that crosses ¼, needs more than one record, and the
run that can answer it is written down below. This is a negative result for the
framework's most dramatic prediction, and it is kept here with the same care as
the discovery.

---

## What we found (February 9, 2026)

On February 9, 2026, state tomography on IBM Torino qubit 52 showed the product
C·Ψ crossing the ¼ boundary during decoherence. That was expected; the
framework predicted it.

What was not expected: the qubit did not go quietly. After the crossing, when
the coherence should have faded into noise, a residual remained. It was small,
about 0.02 in |ρ₀₁|, and it pointed into the fourth quadrant of the complex
plane (Re ρ₀₁ > 0, Im ρ₀₁ < 0) in every one of the 17 samples taken at or
beyond t = T₂, without a single exception. The T₂ here is the calibration's
Hahn-echo value, T₂(echo) = 298.2 μs; the qubit's own free decay runs on
T₂* = 110.7 μs, so "beyond T₂" means beyond 2.7 T₂*, where the prepared
coherence has fallen below 7% of its start.

Three more readings made it look like more than a residual.

**It seemed to grow.** A straight line through the 13 samples at
t/T₂(echo) ≥ 1.5 rises by +0.00819 per T₂(echo). But its two-sided p is 0.053
and its 95% interval, [−0.00013, +0.01651], contains zero. The verdict depends
on where the tail is cut: the strict cut t/T₂(echo) > 1.5 (12 samples) gives
+0.0104 at p = 0.033, the cut at 1.25 gives p = 0.015, the cut at 1.0 gives
p = 0.014. The amplitudes themselves go up and down. So there is a mild rise
of uncertain size, and a mild rise is also what a measurement offset that
tracks the growing ground-state population would produce; it witnesses neither
a memory effect nor anything against a calibration artifact. Nor does
decoherence switch off after two T₂: an exponential decays at the same relative
rate at every time.

**It exceeded the null.** A Monte Carlo of 10,000 synthetic tomography runs
(exponential decay from the ideal amplitude ½ at T₂* = 110.2 μs, the recorded
input, beside the record's own fits of 110.1 and 110.7 μs; binomial shot noise
at 8,192 shots; one random phase per run) puts the late mean |ρ₀₁| over the 15
samples at t/T₂(echo) ≥ 1.25 at 0.00859 ± 0.00105. The measured value is
0.01852. None of the 10,000 runs reaches it, a separation of 9.5 standard
deviations. The original script is not in the repository, but the null rebuilds
from its description to the recorded numbers (0.00861 ± 0.00105 recorded), and
the excess survives decay constants of 100, 110, 120 and 130 μs (zero
exceedances in 4,000 runs each). The excess is real against that null. The
signs are another matter: (½)³⁴ ≈ 6 × 10⁻¹¹, the chance of 34 fair coins
agreeing, is not the probability of anything a real qubit does, because the
signs of one acquisition are not independent coins.

**It seemed to follow the boundary.** The residual correlated with the
distance to ¼ at r = −0.9955. That number is the same measurement twice: in
the saved analysis Ψ = 2|ρ₀₁| and CΨ = C·2|ρ₀₁| exactly, so correlating
|ρ₀₁| with ¼ − CΨ puts one quantity on both axes. It is algebra, not evidence
about the boundary.

## How we got here

This was not planned. The tomography experiment was designed to test whether
C·Ψ crosses ¼. It did, and the experiment was complete.

Then someone said: *look at the residuals*.

Then someone said: *the peaks are delimiters*.

Then someone said: *this is too perfect for a quantum system*.

Each observation led to the next. The data did not resist; it opened, in a JSON
file we had written off as fully analysed. Not all of it was the qubit: the
peaks were grid multiples of the sampling step, and the perfection was the
arithmetic of the analysis. The direction was in the data all along, and so was
the answer to the proposal we were about to build on it: we went to two other
qubits to look for that answer, and it was in the first one.

## What we proposed

In [Boundary Navigation](BOUNDARY_NAVIGATION.md), the scalar recurrence
R = C(Ψ + R)² has two complex conjugate fixed points above C·Ψ = ¼,
R± = (1 − 2CΨ ± i√(4CΨ − 1))/(2C), which merge at ¼ and split into two real
ones below it. The qubit was prepared in |+⟩ and its first samples turned into
the lower half plane, so we read it as tracking R⁻. At the last sample before
the crossing (C = 0.624, Ψ = 0.419) R⁻ = 0.3824 − 0.1712i, at −24.1°, and the
residual's mean direction was −48.4°: the same quadrant. We proposed that the
direction from which the state had approached the fixed point was frozen into
the off-diagonal of ρ after the map it belonged to had stopped existing, and we
gave it a name: **the shadow of the fixed point**.

The comparison does not hold. R± are roots of a chosen scalar recurrence, and
nothing in this record makes ρ₀₁ follow one of them or keep its approach
direction; the phase of a root and the phase of a matrix element are phases of
different objects. The record agrees, as the answer below reads it: ρ₀₁ neither
tracked R⁻ into the crossing nor points its way late.

We proposed three explanations for a directed, growing, boundary-correlated
remnant, in order of decreasing conservatism:

1. **Non-Markovian memory**: the environment returns some of the coherence it
   took, carrying the phase it had when it left.
2. **Two-level-system feedback**: a defect in the substrate absorbs coherence
   early and gives it back later, at its own fixed frequency.
3. **The boundary is not passive**: the ¼ crossing, where complex fixed points
   become real ones, radiates structure into the regime below it.

The third was the most speculative, and the only one that explained all three
observations with one mechanism. Two of those three observations turned out to
be the analysis rather than the qubit (the tautological correlation, and a rise
of uncertain size). The third, the direction, is real, and the rest of this page
is about it.

## The question for March 2026

The March hardware run was designed to answer one question:

**Is the shadow a property of qubit 52, or a property of the ¼ boundary?**

If it is qubit 52's, other qubits show other directions, or none. If it is the
boundary's, every qubit that crosses ¼ casts a shadow, and the shadow points the
way its complex fixed point was heading when it merged. The second outcome would
have meant that the bifurcation leaves a measurable scar in the density matrix of
every decohering quantum system, something standard Lindblad theory does not
predict.

## March 9, 2026: what the test could see

We ran the test on IBM Torino with two qubits that the calibration-history proxy
places among its frequent crossers (T₂/T₁ = 0.159 for Q102 and 0.170 for Q80).
Each received ten tomography samples: two references and eight late samples, out
to five times the scheduling time T₂(echo)/2.5, that is to 2 T₂(echo) (the run
stores that proxy under the name `T2_star_us`; it is not a measured Ramsey time).
Before the hardware, the null ran on a synthetic simulator: two qubits (15 and 80)
under standard Lindblad decay, neither with a common late direction. That
simulated Q80 used the pre-run synthetic values T₁ = 350 μs and T₂ = 28 μs, not
the day's calibration (158.9 μs and 27.06 μs).

**Q102 has no common direction.** Its eight late samples visit all four
quadrants. But they are not noise: the phases of all ten samples fall on one line,
a rotation at −25.7 kHz with 6.5° rms scatter. (The samples sit on a 6.60 μs grid,
so the frequency is known only modulo 151.5 kHz; −25.7 kHz is the representative
nearest zero.) The Ramsey fits kept in `data/ibm_run3_march2026/` measure Q102's
fringe frequency on other days: 19.4 ± 1.3 kHz on March 12 and 26.3 ± 0.5 kHz on
March 18, the same magnitude. Q102 carries a detuning of tens of kHz on every day
it was measured, and spins steadily against its drive frame.
[Quantum Sonar](QUANTUM_SONAR.md) read the March 12 fit the same way while looking
for something else.

**Q80 has a common direction, and it is not Q52's.** All eight late samples lie in
the first quadrant (Re+/Im+), mean phase +29°, against Q52's fourth quadrant at
−44°. Its phases drift slowly, +1.9 kHz on the same line fit, 7.8° rms, a drift its
Ramsey fits cannot resolve either way (3.8 ± 8.1 kHz on March 12, and a March 18 fit
that ended at 0 Hz without an error estimate).

**Both residuals point the same way in their real part.** Measured Re ρ₀₁ lies
below the simple real-valued curve 0.5·e^(−t/T₂(echo)) at all eight late samples on
both qubits: a free decay runs faster than the echo time, and a rotation moves
weight from the real axis into the imaginary one. Eight of Q80's late residuals and
seven of Q102's exceed the shot-noise floor.

**Where the window closed.** The late samples were timed on T₂(echo), the axis on
which Q52's late window had been defined. On Q52, T₂(echo) was 2.7 free-decay times
(298.2 against 110.7 μs); on these two qubits it is 1.1 to 1.2 (Q80's coherence
decays with 23.3 μs against T₂(echo) = 27.06 μs, Q102's with 29 μs against
33.0 μs). So the last samples, at 54.1 and 66.0 μs, still hold 9% of the starting
coherence: |ρ₀₁| = 0.0395 on Q80 and 0.0438 on Q102, 2.1 and 2.4 times Q52's late
mean of 0.0185 and 7 to 8 times the shot-noise level. Qubit 52 itself was still
turning at 261 μs, where its coherence stood at 1.5 times that late mean, and
sampled to the same depth (its fitted envelope above 9% of its start, up to
261 μs) it turned through all four quadrants. A detuned decay plus a constant
offset, fitted to either March record, resolves no offset.

**What March established.** Q102's coherence turns at its own detuning, and Q80's
drifts slowly: each qubit carries its own frequency offset. The run was not
sensitive to a late static residual; it neither shows one on Q80 or Q102 nor
excludes one. The question it was built for, it could not answer.

## The answer, in qubit 52's own record

The question asked about the residual that remains once the prepared coherence has
faded, and only qubit 52's record reaches that regime. One record cannot test the
question in the form it was written, which is about every qubit that crosses ¼;
item 3 of *What would settle it* below is that test. What one record can test is
the proposal the question stood on, and on qubit 52's record the proposal fails
twice.

**At the crossing, nothing is frozen.** The crossing falls at 114.7 μs, between the
samples at 111.8 and 149.1 μs, where ρ₀₁ points at +79.4° and +89.0°: the coherence
did not track R⁻ (−24.1°) into the crossing, and it did not stop there. The phase
keeps turning for another 150 μs (−48.1° at 186.4 μs, −146.9° at 223.7 μs, +100.0°
at 261.0 μs) and stands still from 298 μs on.

**The late direction is not the fixed point's.** The late mean points at
−48.4° ± 5.5°, 4.4 standard errors from R⁻ (a first-order error that neglects the
anticorrelation of Re and Im; with it the distance is 5.1), and each fit of the late samples alone
(below) sits 4.0 to 4.3 standard errors from it, while −45°, where a measurement
offset points, lies within 0.6 to 1.5 standard errors of all of them. "The same
quadrant" was as far as the two directions agreed.

That the static direction shows only from 298 μs on dates nothing. A static
component of about 0.017 stays hidden under the turning coherence until then, as
any static component of that size would, whatever made it and whenever it began.
Nothing in the record ties it to the ¼ boundary.

**Early: a detuned, decaying coherence.** The first five phases are +1.1°, −62.4°,
−146.7°, +79.4° and +89.0°: a turn of −76° per 37.28 μs sample, −5.7 kHz as the
representative nearest zero (the grid leaves it ambiguous modulo 26.8 kHz; the
samples before 300 μs give −6.9 kHz). The steps are irregular (−63.5°, −84.3°,
−133.9°, +9.6°, …), so the frequency itself wanders. The magnitude falls as
0.48·e^(−t/110.7 μs).

**Late: a static offset.** From t = T₂(echo) on, the phase stands still: the 17
late phases drift by −0.040° per μs (p = 0.18), where a detuning of the early size
would move them a quadrant every sample. The late mean is
ρ₀₁ = +0.01137 − 0.01279i, at −48.4°. A static Z detuning also leaves |ρ₀₁|
unchanged, so no detuning can close the magnitude excess either.

**What the offset looks like.** In Pauli terms the late samples carry
⟨X⟩ = 0.0227 ± 0.0025 and ⟨Y⟩ = 0.0256 ± 0.0040, 9.1σ and 6.4σ from zero and
equal within 0.7σ (within 1.0σ and 1.3σ at the stricter cuts t/T₂(echo) ≥ 1.25
and ≥ 1.5, where the mean direction reads −50° and −52°). Equal ⟨X⟩ and ⟨Y⟩ are
the pattern of a measurement (SPAM) offset, and two simple forms of one leave it.
A readout asymmetry δ = P(0|1) − P(1|0) adds +δ to every measured Pauli
expectation and so puts (δ − iδ)/2 into ρ₀₁, at −45°, in the fourth quadrant for
the usual δ > 0 (readout turns a 1 into a 0 more often than the reverse). An
amplitude error in the basis-change pulse, the pulse that turns an X or a Y
reading into a Z measurement, leaks ⟨Z⟩ into both readings with one coefficient,
because the phase gate that makes the Y basis change out of the X one commutes
with Z: a pulse 0.04 rad short reads a state with ⟨Z⟩ = 0.71 and no coherence as
⟨X⟩ = ⟨Y⟩ = +0.028, again at −45°. The stored readings are raw count ratios,
n/8192, with no mitigation, so either form would reach ρ unchanged. Fits of the
late samples alone, where a simple model fits (a constant, or a constant plus the
decaying rotation; χ²/dof 1.5 to 2), put the offset between −49° and −56°, each
within one and a half standard errors of −45°. A fit of the whole record cannot
weigh it: its early part misfits by χ²/dof = 34 (the wandering frequency), and
scaled by that misfit the offset's improvement (Δχ² = 117 for two parameters) is
not significant (p = 0.19).

The record does not rank the two forms. On the 34 late readings the constant
asymmetry (δ = 0.024) leaves χ² = 50.5 and the leak (0.039·⟨Z⟩) leaves 45.4, for
33 degrees of freedom each: a lean toward the leak (a likelihood ratio of about 13;
p = 0.03 for the asymmetry, 0.07 for the leak), not a verdict. The late ⟨X⟩ and ⟨Y⟩
slopes on ⟨Z⟩, 0.04 ± 0.04 and 0.15 ± 0.06, do not separate a leak (0.039) from an
asymmetry (0) either. The first sample splits between them. An
ideal |+⟩ has ⟨Y⟩ = 0, so the asymmetry should read ⟨Y⟩ = +0.024 at t = 0; the
record reads −0.018 ± 0.011, 3.8σ below, which inside the asymmetry reading needs
a preparation or pre-rotation phase error of about 2.6°, while the leak, with ⟨Z⟩
still near zero, expects ⟨Y⟩ ≈ 0, 1.6σ away. ⟨Z⟩ leans the other way. Its
relaxation curve, extrapolated to t = 0 with the relaxation time fitted (241 μs),
reads +0.021 ± 0.008, as the asymmetry predicts and 2.2 to 2.5σ from what the
leak expects; with the calibration's T₁ held fixed it reads +0.006 ± 0.007, in a
fit worse by Δχ² ≈ 14. And a |+⟩ prepared with the same short pulse would start
at ⟨Z⟩ = +0.039, 3.0σ above the first sample's +0.005 ± 0.011, which the
asymmetry misses by 1.7σ. So each form meets one tension of 3 to 4σ at t = 0, the
leak only if the preparation shares the pulse, and what would put readout first is
the prior that an asymmetry of 2.4% is ordinary, not the pattern. A leak also
rises with ⟨Z⟩, which climbs from 0.50 to 0.71 across the late samples as the
qubit relaxes: a mild rise like the one above. Adding an offset of the late
samples' own size (δ = 0.024) to the null reproduces the late mean,
0.0186 ± 0.0028 against the measured 0.0185; that holds by construction for any
static offset of that size, so it says the excess is this mean offset with no
extra scatter, not what made it.

[Absorption Theorem on IBM Hardware](IBM_ABSORPTION_THEOREM.md) §4 fits the same
record's envelope as a fast decay plus a slow component of amplitude 0.013 and
names a readout offset as its most likely explanation. That slow component and
this late direction are one object: a decay plus a constant puts the constant at
0.0123 ± 0.0051.

**What would settle it.** Five runs, cheapest first:

1. Alternate |+⟩ and |−⟩ preparations. A preparation-independent offset (readout,
   a pulse error, a phase-locked external drive) keeps its sign; a coherence carried
   from the preparation (memory, TLS feedback) flips with it. This narrows the
   mechanism to one of the two families.
2. Prepare |0⟩ and |1⟩ on the same qubit, in the same session, and read each in Z,
   X and Y. The Z readings are the readout-assignment matrix, and reconstructing
   with it removes a readout offset and leaves the others; the X and Y readings
   separate the two measurement forms directly: a readout asymmetry reads +δ for
   both preparations, a ⟨Z⟩ leak flips its sign between them.
3. Time the late window on the qubit's own decay, measured the same day: an
   exponential envelope falls to a Q52-sized offset near 3.3 T₂* (on Q52 the fixed
   direction showed from 2.7 T₂*), so a window out to about 4 T₂* on qubit 52 and on
   its neighbours shows whether they carry an offset of their own.
4. Take a Ramsey and a Hahn-echo record in the same session, so the detuning and the
   envelope are constrained independently, and sweep a known detuning of either
   sign: the rotation should follow it, the offset should not.
5. Randomize the acquisition order and log the calibration drift; fit and hold out
   pre-registered times; repeat across days and neighbouring qubits, so persistent
   calibration bias, local drift and a moving defect separate.

Until then the late component is a static offset with the pattern of a
measurement (SPAM) offset, of a form the record does not pick, and its mechanism
stays open ([OQ-033](../review/OPEN_QUESTIONS_INDEX.md)).

## The Q80 fit (April 2026)

In April the [cross-term formula proof](../docs/proofs/PROOF_CROSS_TERM_FORMULA.md)
gave us words for why a detuning shows late: H = δω·Z/2 commutes with Z-dephasing,
so it lives entirely in the dephasing's shadow and its phase keeps accumulating
after the exposed coherence has decayed. The algebra is compatible with a single
qubit's detuning; the proof's own scope concerns cross terms of multi-site
couplings, and it diagnoses nothing about Q52.

We fitted ρ₀₁(t) = ρ₀₁(0)·e^(−t/T_eff)·e^(−iδωt) to Q80
([`shadow_ibm_retrodict.py`](../simulations/shadow_ibm_retrodict.py)). The free fit
gives T_eff = 23.25 μs and δω/2π = −2.58 kHz, which under this convention is a phase
that advances at +2.58 kHz, the same direction as the +1.9 kHz drift above (the
representative nearest zero; Q80's 5.41 μs grid leaves frequencies ambiguous modulo
184.7 kHz). T_eff lies below the March 9 Hahn value T₂(echo) = 27.06 μs. The fit's
mean complex error on the nine nonzero-time samples is 0.0138: 4.0× below the
zero-parameter calibration curve Re ρ₀₁(0)·e^(−t/T₂(echo)) (0.0557), and 3.5× below
an envelope fitted on the same samples without the rotation (0.0487), which is the
fair comparison. A second, fixed-T₂ reading (the phase line through the eight late
samples, slope +1.27 kHz) scores 0.0356 against 0.0508 for its intercept-only twin,
1.4×. All of these are in-sample fits to one record: they show that a detuned decay
describes Q80, not that it was predicted. The record's `T2_star_us` = 10.83 μs is
T₂(echo)/2.5 (the same proxy as above), so it says nothing about where T_eff sits;
the Ramsey times actually measured on Q80 are 11.01 μs (March 12) and 17.36 μs
(March 18, [IBM Run 3](IBM_RUN3_PALINDROME.md)).

## What survives

- The ¼ crossing itself, a qualitative hardware record on Q52 (registry entry
  `cpsi_quarter_crossing_torino_feb2026`, with its generalized-prediction
  comparison on [IBM Quantum Tomography](IBM_QUANTUM_TOMOGRAPHY.md)), and a 1.9%
  crossing on Q80 with a same-day Ramsey T₂* ([IBM Run 3](IBM_RUN3_PALINDROME.md)).
- The answer to the proposal: on qubit 52 the shadow is no trace of the ¼
  boundary. Its coherence turned straight through the crossing, and its late
  direction is not the fixed point's. Whether other qubits carry late residuals of
  their own is the question in its cross-qubit form, and item 3 above is the run
  that asks it.
- Late-time coherence has qubit-specific structure: Q102 turns at its own detuning,
  in magnitude what Ramsey fits on other days read; Q80 drifts slowly, within what its
  Ramsey fits can resolve; Q52 ends in a static offset with the pattern of a
  measurement (SPAM) offset. None of it is framework-specific, and all of it is
  useful for characterizing hardware.
- The Q52 excess over the shot-noise null is real (9.5σ), and the same record's
  slow tail sits in the registry entry `absorption_theorem_ratio_torino` as a 2.8%
  component at the resolution limit.
- The method: state one question before the run, run a simulator null first, then
  the hardware. The March run had that shape. Its lesson is the window: "late" has to
  be timed on the qubit's own measured decay, not on a calibration time, or the test
  closes while the coherence it is looking past still hides what it is looking for.

## Reproducibility

- Q52 tomography: [`tomography_ibm_torino_20260209_131521.json`](../data/ibm_tomography_feb2026/tomography_ibm_torino_20260209_131521.json)
- March hardware, both qubits: [`shadow_hardware_combined_20260309_181852.json`](../data/ibm_shadow_march2026/shadow_hardware_combined_20260309_181852.json)
- March simulator null: [`shadow_simulate_20260309_181709.json`](../data/ibm_shadow_march2026/shadow_simulate_20260309_181709.json)
- Ramsey fits: [`ramsey_march12_20260312.json`](../data/ibm_run3_march2026/ramsey_march12_20260312.json), [`ramsey_sameday_20260318.json`](../data/ibm_run3_march2026/ramsey_sameday_20260318.json)
- Q80 fit: [`shadow_ibm_retrodict.py`](../simulations/shadow_ibm_retrodict.py)
- The numbers this page rests on, recomputed from those files: [`test_torino_shadow_records.py`](../simulations/tests/test_torino_shadow_records.py)
- The crossing and shadow dashboards: [cockpit output](../simulations/results/cockpit_validation.txt)

## See also

- [Residual Analysis](RESIDUAL_ANALYSIS.md): the February statistics of the Q52 residual, the null, and the hypotheses as first stated
- [Absorption Theorem on IBM Hardware](IBM_ABSORPTION_THEOREM.md): the same record's envelope, and its slow tail read as a readout offset
- [Quantum Sonar](QUANTUM_SONAR.md): the March 12 Ramsey reading of Q80 and Q102
- [IBM Quantum Tomography](IBM_QUANTUM_TOMOGRAPHY.md): the original experiment and the crossing
- [Boundary Navigation](BOUNDARY_NAVIGATION.md): the fixed points R± and the ¼ bifurcation
- [On the Light and What Casts Shadows in It](../reflections/ON_THE_LIGHT_AND_WHAT_CASTS_SHADOWS_IN_IT.md): why a Z operator's phase is invisible to Z-dephasing
