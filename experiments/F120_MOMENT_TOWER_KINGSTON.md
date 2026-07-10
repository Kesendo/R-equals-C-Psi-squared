# F120 on hardware: the moment tower read by the chip's own damping (ibm_kingston, 2026-06-11)

**Date:** 2026-06-11 (calibration 06:33Z, main run 07:39Z, arbiter 07:55Z)
**Machine:** ibm_kingston (Heron r2), qubits q149 / q13 / q9, no two-qubit gates anywhere
**Jobs:** d8l6c7rqv2lc73863acg (Arm A), d8l6c832d42s73cb16a0 (Arm B), d8l6h03nn5bs738rmrug (T1 arbiter)
**Data:** [data/ibm_moment_tower_june2026/](../data/ibm_moment_tower_june2026/) (main JSON + arbiter JSON + the calibration snapshot)
**Theory:** [the moment-tower pump-channel proof](../docs/proofs/PROOF_MOMENT_TOWER_PUMP_CHANNEL.md) (F120); pipeline script `run_moment_tower.py` (AIEvolution, imports the repo framework and calls `fw.moment_tower` / `fw.predict_pump_slope` at startup)
**Status:** the structural law is confirmed on hardware. The rate layer's first reading ("q13 violates pump ≤ Γ") was **corrected the same day** by the prep-conditioned re-analysis (§ Correction below): in-situ, every qubit satisfies the bound; the real finding is minute-scale T1 telegraphing on q13 and q9, and the discovery that the protocol measures pump and Γ simultaneously, self-arbitrating.

## What was asked

[F120](../docs/ANALYTICAL_FORMULAS.md) says the device needs no instrument to read the girth ladder's deg-1 tower beyond its own amplitude damping: from the maximally mixed state, d/dt⟨A⟩ = (1/d)·Σ_l Δγ_l·Tr(A Z_l), so the energy-moment polynomials ⟨H_p^j⟩ respond to the chip's damping rung by rung, and the first responding rung is the girth. The protocol is the cleanest we have ever sent to a QPU: X gates to prepare the eight computational basis states (their average is I/d exactly), a delay, single-qubit basis rotations, measure. No entangling gate exists in any circuit, so the three qubits were chosen by quality alone: q149 (T1 = 447.5 μs at calibration), q13 (415.8 μs, the same physical qubit as the [F112 Kingston lens](F112_HARDWARE_LENS_KINGSTON.md)), and q9 (101.5 μs, the deliberate short-T1 contrast).

The measured polynomial is the girth-2 witness H_p = X₀ + X₀Z₁ + 0.7·X₁X₂ with the exact identity H_p² = 2.49·I + 2·Z₁ + 1.4·X₀X₁X₂ and tower t₁ ≡ 0, t₂ = [0, 16, 0]: the first rung must stay silent while the second fires, and the firing reads the *middle* qubit's pump. Two arms permute which physical qubit sits in the middle (Arm A: q13; Arm B: q9), making the site resolution a parameter-free ratio. Three measurement settings cover both rungs; τ ∈ {0, 25, 50, 75, 100, 150} μs; 4096 shots; 288 circuits; dynamical decoupling explicitly disabled (DD's X pulses during the idle would invert the very populations the pump builds).

## Layer one: the structural law, confirmed

- **The double null held.** slope⟨H_p⟩ = +2.4·10⁻⁴/μs (z = +1.47) in Arm A and **−6.8·10⁻⁶/μs (z = −0.04)** in Arm B: the first rung is silent on real hardware, through real idle-time parasitics, exactly as the evolution-blindness predicts (static ZZ and Stark shifts are Z-flavored and cannot reach the slope). The X-strings stayed dark (|⟨XXX⟩| ≤ 0.02 at every τ).
- **The second rung fired, and it is purely the predicted pump.** The operator identity ⟨H_p²⟩ = 2.49 + 2⟨Z₁⟩ + 1.4⟨XXX⟩ holds row by row in the measured table (τ = 100, Arm A: 2.49 + 2·0.4147 + 1.4·0.0186 = 3.345 vs 3.3455 measured): the H²-moment grows as exactly twice the middle site's pump curve, which is the law with j = 2 in action.
- **Site tracking works.** The firing followed the middle-qubit identity across the arms, and the per-qubit early-window pump slopes reproduce across the two arms at **0.3% (q9), 1.9% (q13), 5.7% (q149)**, measured in different chain roles, in different jobs. The girth read from the data is 2, hence m\* = 5 for this witness, exactly the certificate F120 promised.

## Layer two: the rates, and the finding

The parameter-free rate predictions used the morning calibration (Δγ_l = 1/T1_l), and they did not match: Arm A fired ~30% high, Arm B ~20% low, ratio 1.88 instead of 4.10. The first suspect was stale calibration, so a five-circuit standard T1 arbiter ran immediately on the same qubits (|1⟩, delay, measure): **q149 = 424.6 μs, q13 = 430.3 μs, q9 = 99.9 μs**, all consistent with the morning values. No drift. The deviation is physics, not bookkeeping:

| qubit | pump slope (Arm A / B, μs⁻¹) | arbiter 1/T1 | ratio (= z∞ if the noise were pure amplitude damping) |
|---|---|---|---|
| q149 | 2.33·10⁻³ / 2.19·10⁻³ | 2.36·10⁻³ | **0.99 / 0.93, textbook** |
| q13 | 3.03·10⁻³ / 3.09·10⁻³ | 2.32·10⁻³ | **1.30 / 1.33, violation** |
| q9 | 5.79·10⁻³ / 5.78·10⁻³ | 1.00·10⁻² | **0.58, deficit** |

Within *any* two-level Lindblad model with amplitude damping (γ↓, γ↑), dephasing, and unitary terms, the pump slope is γ↓ − γ↑ and the |1⟩-decay rate is γ↓ + γ↑, so **pump ≤ Γ_T1 strictly**, with equality only at zero temperature. The protocol is therefore, unplanned, a *model test*, and the three qubits return three different verdicts:

- **q149 passes**: the textbook qubit, pump/Γ = 0.93-0.99, the law and the noise model both hold.
- **q13 violates the inequality at 4-6σ**, reproducibly (1.9% across arms): its Z-pump runs 30% *faster* than its relaxation rate, which no thermal correction can produce (z∞ ≤ 1). Something non-unital beyond two-level amplitude damping pumps q13 toward |0⟩: candidate mechanisms are leakage through the third transmon level, a TLS, or correlated channels. Consistent with this, q13's Arm-A pump curve shows a telegraph-like jump between τ = 75 and 100 μs (0.23 → 0.41) that Arm B, minutes earlier, does not have: the extra channel appears to switch.
- **q9's deficit (0.58) would mean 21% thermal population if read as temperature, which is absurd**; the arbiter resolves it differently: q9's |1⟩-decay is itself non-exponential (local rate 13.4·10⁻³ → 7.7·10⁻³ μs⁻¹ across the window), so "1/T1" is not a well-defined number for this qubit (calibrated T2 = 28.5 μs ≪ T1 already flagged it). The pump channel and the |1⟩-decay weight the same non-exponential bath differently, and neither is wrong.

The honest summary of layer two: **F120's slope reads Tr(Z_l·D(I)), the device's true non-unital pump vector, and on this chip the true pump vector disagrees with calibrated amplitude damping on two of three qubits.** Equating the pump with 1/T1 was the textbook-model assumption, and the chip declines it, per qubit, with signatures (an inequality violation, a non-exponential arbiter, a telegraph switch) that point at specific physics beyond the model. This is the same shape as the [F112 Kingston reading](F112_HARDWARE_LENS_KINGSTON.md), where the lens found a transverse-field anomaly: the channel works, and what it reads through is the gap between the chip and its datasheet.

## Correction, same day: the violation was an epoch artifact, and the protocol is its own arbiter

The model-test table above compares the run's pump slopes against the arbiter's Γ, measured 16 minutes later. That comparison contains a hidden variable, and the run's own data exposes it. The pump protocol prepares all eight basis states, so conditioning each qubit's ⟨Z⟩(τ) on its own preparation bit splits the same circuits into a |0⟩-branch and a |1⟩-branch: s₀ measures −2γ↑ (heating only; in any two-level model ⟨Z⟩ cannot rise from the ground state), s₁ measures the decay, pump = (s₁+s₀)/2 and **Γ = (s₁−s₀)/2 in-situ, from the same shots, the same minutes**. The inequality pump ≤ Γ is equivalent to s₀ ≤ 0. (One bit-ordering gotcha for reproducers: the saved count keys are little-endian; the mapping is verified bit-exactly against the unconditioned curves, e.g. q13's mixed ⟨Z₁⟩ at τ = 50 equals (0.9947 + (−0.6737))/2 = 0.1605 exactly. Re-analysis: [`simulations/f120_prep_split_reanalysis.py`](../simulations/f120_prep_split_reanalysis.py).)

The in-situ table (early window τ ≤ 75 μs, readout-corrected):

| qubit | pump (Arm A / B, μs⁻¹) | Γ in-situ (A / B) | pump/Γ (A / B) | s₀ |
|---|---|---|---|---|
| q149 | 2.33·10⁻³ / 2.19·10⁻³ | 2.38·10⁻³ / 2.25·10⁻³ | 0.979 / 0.972 | < 0 ✓ |
| q13 | 3.03·10⁻³ / 3.10·10⁻³ | **3.14·10⁻³ / 3.21·10⁻³** | **0.965 / 0.966** | < 0 ✓ |
| q9 | 5.79·10⁻³ / 5.78·10⁻³ | **5.83·10⁻³ / 5.80·10⁻³** | **0.994 / 0.996** | < 0 ✓ |

**No violation, anywhere.** Every qubit sits 1-3% below the bound, and the margin is itself a measurement: −s₀/2 = γ↑ reads the per-qubit thermal excitation directly (q13 ≈ 1.8%, q149 ≈ 0.7-2%, q9 ≈ 0.3%, realistic transmon numbers). What the cross-epoch comparison actually detected is **minute-scale T1 telegraphing**: during the run q13's Γ was 3.14-3.21·10⁻³ (T1 ≈ 315 μs) and by the arbiter it was 2.32·10⁻³ (T1 = 430 μs, matching the morning calibration); q9 switched the other way (in-run T1 ≈ 172 μs, arbiter early-window ≈ 75 μs); q149 held still, which is why it looked "textbook" in the flawed comparison. All three earlier verdicts (violation, deficit, pass) collapse into one statement: **two of three qubits telegraph their T1 between epochs, the two-level model holds within each epoch, and the epoch was the hidden variable.** q13's mid-run jump (τ = 75 → 100, Arm A, localized in the |1⟩-branch) is the telegraph caught inside a single arm.

The upgrade this leaves behind: the pump protocol never needed the separate arbiter. Prep-conditioning makes it **self-arbitrating** (pump, Γ, and γ↑ from one set of circuits, epoch-matched by construction), and the model test in its sharp in-situ form is simply s₀ ≤ 0 per qubit per block. The original table above is kept as the honest record of how the wrong reading arose.

## The telegraph chase (same day, 12:02Z): no qubit is special, time is

The follow-up run (`run_moment_tower_telegraph.py`, job d8l8f7r2d42s73cb3q7g; a first submission d8l78s32d42s73cb2b00 sat RUNNING for 90 minutes with zero quantum seconds, an open-plan fair-share stall, and was cancelled; the Batch-mode resubmission ran cleanly) made the epoch the measured axis: 16 blocks inside one job, each block reading s₀ and s₁ per qubit from two preparations (|000⟩, |111⟩) over τ ∈ {0, 20, 40, 70} μs, 8 circuits per block, self-arbitrating by construction.

**Within the job (~75 s execution): everything still.** All three qubits' Γ(b) traces are flat at the shot-noise floor (block-to-block sd 0.15-0.21·10⁻³ against floors 0.19-0.63·10⁻³), no switch was caught in the act, and the in-situ model test s₀ ≤ 0 passed in 47 of 48 block-qubit cells (the one excursion, q149 block 3 at +2σ, is the expected false-positive rate). Per-block thermometry stayed at the few-tenths-of-a-percent level.

**Across the day's three epochs: every qubit moved, including the control.**

| qubit | 07:39 in-situ Γ (μs⁻¹) | 07:55 arbiter | 12:02 in-situ Γ | bias-corrected T1 now |
|---|---|---|---|---|
| q149 | 2.33 / 2.19·10⁻³ | 2.36·10⁻³ | **3.10·10⁻³** | ≈ 285 μs (was ≈ 430) |
| q13 | 3.14 / 3.21·10⁻³ | 2.32·10⁻³ | **4.31·10⁻³** | ≈ 197 μs (a third state: 430 ↔ 315 ↔ ~200) |
| q9 | 5.83 / 5.80·10⁻³ | 1.00·10⁻² | **6.95·10⁻³** | ≈ 108 μs (back near the arbiter state) |

(The 07:39 and 12:02 columns share the same window-matched linear estimator and compare directly; the arbiter column is a global exponential fit and is only roughly comparable. The bias correction inverts the finite-window factor calibrated in the script's simulate mode.)

The chase's verdict closes the q13 investigation: **q13 was never special.** Its T1 has now been seen in three states across five hours, but so has q149's, the qubit the morning's comparison crowned "textbook-stable" on the strength of two epochs that happened to agree. On this device the T1 landscape breathes everywhere, by 1.5-2×, on timescales between minutes and hours, while within any ~minute window the two-level model holds cleanly (flat Γ, s₀ ≤ 0, thermometry steady). The hidden variable of the original "violation" was time, and the prep-split protocol is the right instrument precisely because it is epoch-immune: pump, Γ, and γ↑ from the same shots, every time. Data: `moment_tower_telegraph_ibm_kingston_20260611T100206Z.json` in [data/ibm_moment_tower_june2026/](../data/ibm_moment_tower_june2026/).

## What this confirms, and what it opens

Confirmed (registered as `f120_moment_tower_kingston_june2026` in the Confirmations registries): the pump-slope law's structure on hardware: the double null (z = 1.47 / 0.04), the row-exact H² identity, the girth-2 pattern (rung 1 silent, rung 2 firing), site tracking across arms, per-qubit pump rates reproducible to 0.3-5.7%, and, after the same-day correction, the **in-situ pump ≤ Γ bound holding on all three qubits at 1-3% margins that read the per-qubit thermal population**. The deg-1 hardness rung of a programmed Hamiltonian is, as of this run, a quantity a quantum computer measures about itself by doing nothing but decaying.

Opened, then partly answered the same day: the time-resolved chase RAN (the telegraph section above): no switch inside a ~75 s window, every qubit moved across the day's epochs, so the breathing timescale is bounded between ~minutes and hours and a longer-span monitor (blocks spread over tens of minutes, e.g. across several Batch jobs) is the remaining instrument to catch a switch in the act. Still open: the prep-conditioned **self-arbitrating form of the protocol** as the standing cheap noise diagnostic (pump, Γ, γ↑ from one circuit set, no tomography, no separate T1 experiment, the in-situ model test being s₀ ≤ 0 per block); the F113 complementarity prediction, that a deliberately injected X/Y-drive shifts the pump *curvature* while a Z-drive cannot, as the second hardware arm with a validated baseline; and the d = 3 response-chain reader for the higher moment classes.

## Honest fences

The girth certificate is one-sided by design (a silent deg-1 tower proves nothing about softness; the deg-5 control IIXY+ZXZY lives in the verifier, not on hardware yet). The exponential full-curve fits are degenerate when the window covers < 50% of saturation (q13's z∞ railed in Arm A); the early-window linear slopes are the conditioned readout and carry the conclusions, and they carry a finite-window bias common to all blocks and arms (calibrated in the telegraph script's simulate mode and inverted only for the bias-corrected T1 column), so **rates compare only between same-estimator, same-window readings**; the arbiter's global exponential fit is a different estimator and only roughly comparable, which is part of how the original misread arose. The F113 static twin was reported per arm (−4³·slope⟨H⟩ = −1.5·10⁻² and +4.4·10⁻⁴ against a prediction of 0) and is consistent with the null at the same significance as the slopes themselves. Readout correction uses the morning assignment matrices (errors ≤ 1%; slopes are τ-differences and largely immune, and no readout systematic was needed to explain anything once the epochs were matched). The telegraph chase's per-block fits rest on 4 τ-points × 2048 shots, resolving Γ steps of ≳ 0.5·10⁻³ μs⁻¹ per block and nothing finer; its ~75 s execution window bounds within-job switching only from below, and the day's three epochs are three samples, not a statistics of the breathing. And the morning's verdict "q149 stable" is itself fenced: it rested on two epochs that happened to agree, and the 12:02Z epoch retired it.
