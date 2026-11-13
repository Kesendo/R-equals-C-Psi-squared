# IBM Marrakesh F83 4-Hamiltonian Signature Test, April 2026

Live hardware verification of F83's structural prediction that the four
F87 Π²-classes (truly / pure Π²-odd / pure Π²-even-non-truly / mixed)
produce **operationally distinguishable Pauli-expectation patterns** on
N=3 chains, read through the F-toolkit lens.

Run on `ibm_marrakesh` (Heron r2) on path **[4, 5, 6]** (top-ranked
by 2026-04-30T16:25Z calibration, score 0.0162, rank 1 of 223
3-qubit chains).

## Files

| File | Description |
|------|-------------|
| `f83_signature_ibm_marrakesh_20260430_190035.json` | Raw counts and reconstructed 2-qubit Pauli expectations for **4 Hamiltonians × 9 tomography bases** on (q0=4, q2=6); 36 circuits, 4096 shots/circuit. |

## Experiment summary

- **Backend:** ibm_marrakesh (Heron r2). Same backend as April 26 soft_break and other April 2026 runs.
- **Date:** 2026-04-30, ~19:00 UTC.
- **Job ID:** `d7pol1e7g7gs73cf7j90`.
- **Path:** [4, 5, 6]. Calibration: T1 = 206 / 144 / 351 μs, T2 = 184 / 121 / 151 μs, RO err = 0.42 / 0.38 / 0.37 %, CZ err = 0.633 % (4-5) / 0.287 % (5-6).
- **Initial state:** \|+−+⟩ X-Néel (Hadamard on q0, q2; Hadamard+Z on q1).
- **Four Hamiltonians:**
  - `truly_unbroken` = XX+YY (truly; M = 0 idealized)
  - `pi2_odd_pure` = XY+YX (pure Π²-odd; F83 anti-fraction = 1/2)
  - `pi2_even_nontruly` = YZ+ZY (pure Π²-even non-truly; F83 anti-fraction = 0)
  - `mixed_anti_one_sixth` = XY+YZ (equal-Frobenius mixed; F83 anti-fraction = 1/6)
- **Evolution:** uniform J = 1, t = 0.8, n_trotter = 3 (first-order Trotter via Qiskit `PauliEvolutionGate`).
- **Tomography:** 9 Pauli bases on (q0, q2), 4096 shots/basis. Total 36 circuits.
- **QPU cost:** ~3 minutes wall-clock.

## Reading through the F-toolkit lens

Each Π²-class is its own diagnostic lens. The framework tells us in
advance what each class reveals:

| Π²-class | F-toolkit reading | Idealized M | What it shows |
|----------|--------------------|-------------|----------------|
| `truly_unbroken` (XX+YY) | F87 truly: M = 0; F84 says environmental σ⁻ amplitude damping would break the ⟨Z⟩-conservation the Hamiltonian stipulates | 0 | The environment, cleanly: truly, Z-dephasing and the state respect X^⊗3, so a Z-only model holds ⟨Z₀⟩ = ⟨Z₂⟩ = 0, and an F82/F84 σ⁻ contribution would show there. |
| `pi2_odd_pure` (XY+YX) | F87 pure-Π²-odd; F81 says M_anti = L_{H_odd} dominates; F83 anti = 1/2 | M_anti dominant | The dynamics-generator side of M. F83 anti-fraction is the structural prediction. |
| `pi2_even_nontruly` (YZ+ZY) | F87 pure-Π²-even-non-truly; F81 says M_sym dominates; F83 anti = 0 | M_sym dominant | The mirror-echo side of M (memory-only, no drive). |
| `mixed_anti_one_sixth` (XY+YZ) | F83 anti = 1/6 with r = 1 (equal Frobenius mix) | M_anti and M_sym in fixed ratio | The continuous tunability of F83 at one decimal point. |

This guides what we measure and how to read it. The two pure-class
categories are the F83 anchors. truly is the lens for environmental
F82/F84. mixed tests F83 at a non-trivial ratio.

## The prediction is the circuit the hardware ran

The Trotter + γ_Z prediction (`fw.predict_signature_table`) runs the
same first-order Trotter step as the hardware: the runner
`run_soft_break.py` builds one `SparsePauliOp` with the bond loop outside
the term loop, and Qiskit's default Lie-Trotter synthesis applies the
terms in that list order. The framework's step is gated against the
synthesised circuit (`gate.definition`) in
[`test_predict_signature_table.py`](../../simulations/framework/tests/workflows/test_predict_signature_table.py).
On N = 3 with two terms A + B the two orders differ by one swap: terms
outside runs A₀₁ A₁₂ B₀₁ B₁₂, bonds outside runs A₀₁ B₀₁ A₁₂ B₁₂, so the
order is visible exactly when B on bond (0,1) fails to commute with A on
bond (1,2), at O(δt²) with δt = 0.267. They meet on the shared site 1:
truly (Y against X) and mixed (Z against X) anticommute there, pi2_odd_pure
(X against X) and pi2_even_nontruly (Y against Y) commute, and for those
two the order is invisible. The data choose the same order on their own:
over the 28 datapoints of this run the bonds-outside prediction fits at
RMS 0.024, the terms-outside one at best 0.127 (on April 26, over 45
datapoints, 0.039 against 0.081).

## F83 anti-fraction discrimination: CONFIRMED

The four F87 Π²-classes are operationally distinguishable on hardware
via Pauli-expectation pattern. Each category has at least one
unique-fingerprint observable that separates it from all three other
categories at >>10σ (statistical error 0.0156 at 4096 shots):

- `truly_unbroken`: ⟨Y₀Z₂⟩ = +0.67 (others all near 0); ⟨X₀Z₂⟩ ≈ 0
- `pi2_odd_pure`: ⟨X₀Z₂⟩ = −0.85 (uniquely large negative)
- `pi2_even_nontruly`: ⟨X₀X₂⟩ = +0.92 (uniquely large positive)
- `mixed_anti_one_sixth`: ⟨Z₀X₂⟩ = −0.72, ⟨X₀Z₂⟩ = +0.15 (sign-flip on Z,X axis vs pi2_odd)

Quantitative match to the Trotter + γ_Z model at the path-fit γ_Z = 0.03
(the minimum of the 28-datapoint sweep below):

| Category | RMS over 7 key Paulis | Reading |
|----------|------------------------|---------|
| `truly_unbroken` | **0.012** | M = 0 lens: on the prediction over these 7 two-body Paulis; the single-site window is read below. |
| `pi2_odd_pure` | **0.036** | F83 anchor at anti = 1/2: on the prediction. |
| `pi2_even_nontruly` | **0.014** | F83 anchor at anti = 0: on the prediction. |
| `mixed_anti_one_sixth` | **0.024** | F83 at anti = 1/6: on the prediction. |

All four classes land on the framework's prediction to within 0.04 RMS
over 7 observables, a class RMS between 0.8 and 2.3 shot-noise σ. The
fit has one parameter, γ_Z, taken from these same 28 points, and single
observables still reach 4.3σ (pi2_odd_pure ⟨Y₀Z₂⟩ = +0.067 where the
model gives 0), 3.4σ (pi2_odd_pure ⟨Z₀Y₂⟩ = +0.053 against 0) and 3.1σ
(mixed ⟨Z₀X₂⟩ = −0.721 against −0.769); what
carries those we do not know. On 2026-10-29 the pi2_odd_pure circuits
ran again on the same path, gate for gate (with an explicit schedule and
q1 read as well), and read ⟨Y₀Z₂⟩ nonzero again; a fresh compile of the
same logical circuit read it near zero. Every fit arm of that flight
carries an odd-Y excess, and no single one of nineteen physical letters
accounts for it
([the reality-window flight](../../experiments/REALITY_WINDOW_HARDWARE_PREDICTION.md)).
Away from these forbidden readings, the Trotter + γ_Z model that carries F83's
structural classes holds on hardware in every class, not only in the two
pure anchors; the anti-fraction r itself is not measured (see below).

## The X^⊗3 window: quiet on truly, open on pi2_even

Because truly XX+YY has M = 0 idealized, it is the lens for the
environmental contribution F82/F84 describes:

- F82 closed form: ‖D_T1_odd‖_F = γ_T1 · √N · 2^(N-1). T1 amplitude damping leaks into M_anti.
- F84: among the single-qubit amplitude-damping and Pauli channels, only σ⁻ and σ⁺ reach M_anti; in general a traceless jump does exactly when its Pauli parts mix bit_b parity (the F82 proof, "Other dissipators"). D[Z], D[X], D[Y] are Π²-symmetric and contribute zero. The T1 channel decomposes as σ⁻ (cooling) + σ⁺ (heating); their imbalance at T → 0 is what reaches M_anti.

Which observable carries it is a symmetry question. XX+YY commutes with
X^⊗3, Z-dephasing respects X^⊗3, and |+, −, +⟩ is an X^⊗3 eigenstate, so
any Z-only model keeps the single-site ⟨Z⟩ and ⟨Y⟩ at zero exactly (both
anticommute with X^⊗3). σ⁻ breaks X^⊗3 and pulls ⟨Z⟩ positive. ⟨Z₀Z₂⟩
commutes with X^⊗3 and is nearly blind to σ⁻ on truly:

| truly, γ_Z = 0.03 | γ_T1 = 0 | γ_T1 = 0.046 | γ_T1 = 0.2 | HW |
|---|---|---|---|---|
| ⟨Z₀⟩ | 0 | +0.037 | +0.152 | +0.043 |
| ⟨Z₂⟩ | 0 | +0.040 | +0.162 | +0.012 |
| ⟨Z₀Z₂⟩ | +0.231 | +0.221 | +0.208 | +0.215 |

(σ⁻ per site at a uniform rate γ_T1 inside each Trotter step,
`marrakesh_t1_amplification_test.predict_trotter`; 0.046 is the scale the
F88b baseline derives from the 2026-04-25 Marrakesh calibration, not a fit
to this path.)

On truly, ⟨Z₀⟩ reads 2.75σ above zero and ⟨Z₂⟩ 0.8σ in the one-basis
estimator (pooled over three bases, below, both sit near 3σ). The readout asymmetry of the 2026-04-30 calibration (P(0|1) − P(1|0) =
0.0059 on Q4, 0.0039 on Q6) accounts for +0.006 and +0.004 of that, not
for the rest. Calibration-scale σ⁻ would fit as well as none: whether
there is a T1 contribution on truly we do not know.

YZ+ZY commutes with X^⊗3 too, and there the window is open. Its
single-site readings, all exactly 0 in a Z-only model:

| pi2_even_nontruly | HW | σ from 0 | σ⁻ at γ_T1 = 0.046 |
|---|---|---|---|
| ⟨Z₀⟩ | +0.019 | 1.2 | +0.011 |
| ⟨Z₂⟩ | +0.063 | 4.0 | +0.007 |
| ⟨Y₀⟩ | +0.048 | 3.1 | 0 |
| ⟨Y₂⟩ | +0.045 | 2.9 | 0 |

Something breaks X^⊗3 on this class at 4σ, and it is not calibration-scale
σ⁻: σ⁻ leaves ⟨Y⟩ at zero, and gives ⟨Z₂⟩ a tenth of what is measured.
Nor is it the Trotter circuit: YZ and ZY each commute with X^⊗3, so every
Trotter factor keeps it, as the continuous evolution does
([the palindrome as a colouring](../../experiments/THE_PALINDROME_AS_A_COLOURING.md),
its paragraph "Against each term"); the native rz, sx and cz that compile a
factor do not commute with X^⊗3 one by one, so a coherent error inside a
factor can break it. Stochastic Pauli noise between the factors cannot
(depolarizing, dephasing along any axis, incoherent crosstalk): a Pauli
channel commutes with conjugation by X^⊗3, and so does a preparation error
that mixes eigenstates of X^⊗3. Inside a factor it is carried through native
gates that need not commute with X^⊗3; in the flown circuits the only rz that
are no Clifford are the factors' own rotations and one small pair in the
preparation, so this enters at most at the size of that pair. What breaks
the window is coherent and not a Pauli channel, or not unital, as σ⁻ is.
Pooled over the three bases that read each letter (the table above takes one
basis each), the four readings are ⟨Z₀⟩ +0.023, ⟨Z₂⟩ +0.057, ⟨Y₀⟩ +0.039,
⟨Y₂⟩ +0.042 at a binomial σ ≈ 0.009. The three bases scatter more than that:
over the 24 single-site letters of truly and YZ+ZY on this run and the
Kingston run below, their spread is 74 on 48 degrees of freedom, about 1.24
times the binomial σ, which leaves ⟨Z₂⟩ near 5σ and truly's pooled ⟨Z₀⟩
near 3σ.

The break is not this chip's alone. Kingston on 2026-05-05 flew the same
four classes on path [43, 56, 63]
(`data/ibm_soft_break_april2026/soft_break_ibm_kingston_20260505_102806.json`),
and pooled the same way YZ+ZY reads ⟨Z₀⟩ +0.061 (6.8σ binomial), ⟨Y₂⟩ +0.048
(5.3σ), ⟨Y₀⟩ +0.022 and ⟨Z₂⟩ +0.009, while truly reads nothing beyond 2σ
there; its ⟨Y₀X₂⟩ reads −0.022 where the continuous evolution makes it equal
to ⟨Y₀⟩, so its ⟨Y⟩ rests mainly on ⟨Y₂⟩. The large ⟨Z⟩ sits on q2 on Marrakesh and on q0 on Kingston, which
points at something local to a qubit.

Which letters are lit says what kind of term it is. On YZ+ZY from |+, −, +⟩,
a small Z field on the end qubits lights only ⟨Y⟩, ⟨X₀Y₂⟩ and ⟨Y₀X₂⟩; a ZZ
coupling only ⟨Y₀Z₂⟩ and ⟨Z₀Y₂⟩; a turn about Y only ⟨Z⟩, ⟨X₀Z₂⟩ and
⟨Z₀X₂⟩ (continuous evolution, checked numerically). Both runs light the
first and the third group, and ZZ is silent (⟨Y₀Z₂⟩ +0.001, ⟨Z₀Y₂⟩ −0.006
here). In the circuit a detuning is not that Z field: between the pulses the
frames are turned, and the fit below shows a detuning lighting mostly ⟨Z⟩.
For ⟨Y⟩ there is a lever: at the end of the circuit the end qubits of YZ+ZY
keep ⟨X⟩ ≈ 0.9, truly's about 0.14, so one small turn about Z shows some
seven times larger on YZ+ZY, and truly's ⟨Y⟩ stays near zero. For ⟨Z⟩ the
lever does not account for the classes: truly's ⟨Z₀⟩ here is as large as
YZ+ZY's, as calibration-scale σ⁻ would make it (the table above). A tilt
common to every sx and x pulse on a qubit cannot carry it: it is a turn of
that qubit's drive frame, the identity on the logical circuit
([the driver's frame](../../experiments/DRIVE_FRAME_HARDWARE_PREDICTION.md)
§1), and the replay gives it machine zero. A tilt of the sx pulses against
the x pulses can, 0.06 rad on q4 and q6 moving the forbidden readings that
carry no X letter by up to 0.09, and so can the pulses' amplitude, 1 % on q4 by up to 0.05 and on
q6 by 0.014 (arm B of the pre-registration linked below, under the
October snapshot; [`x3_window_power.py`](../../simulations/x3_window_power.py)). On
the flown circuits, replayed with the day's T1, T2 and readout as in
[REALITY_WINDOW_HARDWARE_PREDICTION](../../experiments/REALITY_WINDOW_HARDWARE_PREDICTION.md)
§2, against the 18 readings these two classes forbid (16 by X^⊗3, two of
YZ+ZY's by the parity of n_Y), a detuning of the two end qubits inside its
bound takes the sum of squared residuals from 94 to 70 (76 with the gates as
soon as possible); it gives YZ+ZY's ⟨Z₂⟩ +0.034 against +0.057 but ⟨Y₀⟩
+0.014 and ⟨Y₂⟩ +0.003 against +0.039 and +0.042, and leaves truly's ⟨Z⟩ at
zero. All nineteen letters of that document's §5, each held inside its
bound, bring it to 3.6 (4.2), with more letters than readings and seven of
them at their bounds, which decides nothing; that document's own fit of all
four classes reached its small sum only with a detuning far beyond its
bound. The readings share shots,
some strongly (⟨Z₀⟩ and ⟨Z₀X₂⟩ come from the same counts), so these sums
are no χ² of a known distribution
([`f83_xn_window.py`](../../simulations/f83_xn_window.py) →
[`f83_xn_window.txt`](../../simulations/results/f83_xn_window.txt)). What
breaks it we do not know. One time point cannot separate the candidates, and
a scan over the circuit's length does not either: more Trotter steps add
gates and time together, and over a few µs against T1 of 144 to 351 µs
damping grows as linearly as the rest. A sign twin from |−, +, −⟩ and a
twirl of the cz gates do not separate the candidates on a fresh compile
of these circuits either: the twin changes the sign of some single-site readings
and keeps others (under a phase of 0.05 after every cz ⟨Z₂⟩ goes from
−0.34 to +0.11 while ⟨Z₀⟩ goes from −0.155 to −0.136), and a twirl suppresses
single-qubit letters as well
([`x3_window_arms_scout.py`](../../simulations/x3_window_arms_scout.py)).
What is pre-registered instead is
[the X^⊗3 window under a turned frame](../../experiments/X3_WINDOW_HARDWARE_PREDICTION.md):
these circuits read with one qubit's drive frame turned, which Z letters and
own-frame letters do not see and a foreign transverse letter does, beside
the idle detuning of each qubit measured directly.

⟨Z₀Z₂⟩ separates the classes only under σ⁻: truly and pi2_odd_pure predict
the same +0.231 in a Z-only model, but σ⁻ moves pi2_odd's (0.193 at
γ_T1 = 0.046) more than truly's (0.221). The measured +0.215 and +0.210 sit
about 1σ from either hypothesis.

## The Y/Z asymmetry on truly is the circuit's

The truly Hamiltonian XX+YY and the |+, −, +⟩ initial state are both
symmetric under the reflection q0 ↔ q2; the Trotter circuit is not. Each
step runs XX₀₁ YY₀₁ XX₁₂ YY₁₂; its mirror image XX₁₂ YY₁₂ XX₀₁ YY₀₁ is a
different product, because the cross-bond pairs (YY₀₁ with XX₁₂, XX₀₁
with YY₁₂) anticommute on the shared site.
(With the terms outside, XX₀₁ XX₁₂ YY₀₁ YY₁₂ reflects into itself, since
factors of one term commute, and the prediction would be symmetric.) The
prediction already carries the asymmetry:

| | Pred (γ_Z = 0.03) | HW |
|---|---|---|
| ⟨Y₀Z₂⟩ | +0.659 | +0.670 |
| ⟨Z₀Y₂⟩ | +0.192 | +0.185 |

The measured asymmetry is the predicted one to within 0.02, about one
shot-noise σ. A contribution from per-qubit T2 differences on the path
(T2 = 184 μs on Q4 against 151 μs on Q6) is bounded by that, not needed.

## Residual phenomenology: γ_Z_eff differs between the two runs

A post-hoc γ_Z sweep (`simulations/f83_gamma_z_sweep.py`) over the
28-datapoint set of this run (4 classes × 7 Paulis), and over the 14
comparable datapoints of the April 26 run (truly and XY+YX):

| Path | Date | Best γ_Z_eff | RMS at best | Anchor XY+YX ⟨X₀Z₂⟩ pred | HW |
|------|------|---------------|-------------|---------------------------|-----|
| [4, 5, 6] | 2026-04-30 | 0.03 | 0.024 | −0.835 | −0.849 |
| [48, 49, 50] | 2026-04-26 | 0.10 | 0.037 | −0.723 | −0.711 |

The two runs fit different effective dephasing rates, by a factor of
about three, and both fits are good. They differ in path, day and job, so
the difference is not isolated to the path. γ_Z_eff is a phenomenological
parameter that absorbs coherent gate errors, transpilation overhead,
idle decoherence and crosstalk; it does not reduce to a single hardware
time constant (Q5 has the shortest T2 on [4,5,6] but the path fits the
*lower* γ_Z_eff). What sets the difference between these two runs we do
not know.

## What this run does NOT establish

- That F83's anti-fraction r = ‖H_even_nontruly‖² / ‖H_odd‖² is invertible
  from data. The prediction is structural (closed-form on H), not fitted
  from observables. We measured signatures, not r.
- Anything quantitative about F82/F84. truly's single-site ⟨Z₀⟩, ⟨Z₂⟩
  sit at most near 3σ, and the X^⊗3 break on pi2_even_nontruly is not
  calibration-scale σ⁻ and not the Trotter circuit; what opens it is not
  known (the X^⊗3 window above).
- That the discrimination scales to N > 3. Open.
- That arbitrary mixed Hamiltonians (other than 1/2, 0, 1/6) follow the
  same pattern. Three F83-discrete cases tested; the continuous tunability
  of r remains a structural prediction, not a hardware-mapped one.

## Reproducing the run

```bash
cd "D:\Entwicklung\Projekte\.NET Projekte\AIEvolution\AIEvolution.UI\experiments\ibm_quantum_tomography"
python run_soft_break.py --hardware --backend ibm_marrakesh --path 4,5,6 --shots 4096
```

The run-time on Marrakesh queue + execution: ~3 minutes total.

## Reading

- [PROOF_F83_PI_DECOMPOSITION_RATIO](../../docs/proofs/PROOF_F83_PI_DECOMPOSITION_RATIO.md): the closed-form anti-fraction.
- [PROOF_F82_T1_DISSIPATOR_CORRECTION](../../docs/proofs/PROOF_F82_T1_DISSIPATOR_CORRECTION.md): the T1-into-M_anti closed form; on this run truly's single-site ⟨Z⟩ sits at most near 3σ.
- [PROOF_F84_AMPLITUDE_DAMPING](../../docs/proofs/PROOF_F84_AMPLITUDE_DAMPING.md): among the single-qubit amplitude-damping and Pauli channels, only σ⁻/σ⁺ reach M_anti.
- [ON_THE_RESIDUAL](../../reflections/ON_THE_RESIDUAL.md): consolidating reflection on F80–F85.
- [ibm_soft_break_april2026/](../ibm_soft_break_april2026/): the April 26 anchor run on path [48, 49, 50] that this test extends.
- [`f83_signature_predictions.py`](../../simulations/f83_signature_predictions.py): closed-form Trotter+γ_Z=0.1 predictions per category.
- [`f83_aer_preflight.py`](../../simulations/f83_aer_preflight.py): Aer noise simulation with [4, 5, 6] calibration values; matched Trotter to 0.001 on the anchor observable; both sat 0.13 from the hardware there, the gap the path's lower γ_Z_eff closes.
- [`f83_hy_field_check.py`](../../simulations/f83_hy_field_check.py): tested h_y = 0.05 (zn_mirror April 29 finding) as alternative explanation; rejected (h_y rotates Z↔X, doesn't amplify Y).
- [`f83_gamma_z_sweep.py`](../../simulations/f83_gamma_z_sweep.py): per-path γ_Z_eff fit and per-observable residual diagnostic.
- [`2qubit_dissipator_exploration.py`](../../simulations/2qubit_dissipator_exploration.py) (commit ba8e861): partial F86, single-bond closed form for σ-channels; overlap structure open.
