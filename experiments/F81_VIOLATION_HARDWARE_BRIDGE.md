# F81 Violation: the Hardware Bridge (from operator diagnostic to measurable number)

**Status:** The f81_violation discriminator, until now a post-fit operator diagnostic computed from a Liouvillian we build in simulation, is grounded to a measurable object and read out on existing hardware data for the first time. Three legs, all zero QPU: (1) the method demo on the F113-fitted Kingston Lindbladians, with the tautology of a parameterized fit made explicit; (2) the from-below grounding: the violation IS 2^(N−1) times the RMS velocity with which the maximally mixed state polarizes, per site the net cooling flux γ↓−γ↑; (3) the first non-tautological readout on the Marrakesh price_pair T1-leg data via the asymptote recipe a = z∞/T1, no Lindblad fit anywhere. The missing half (a |0⟩ heating leg that separates γ↓ from γ↑ and settles the z∞ attribution) is simulator-validated and pre-registered (same day, below; artifacts in `data/ibm_heating_leg_july2026/`); the hardware run waits for sane queues and Tom's go.
**Date:** 2026-07-05
**Authors:** Thomas Wicht, Claude (Fable 5)
**Scripts:** [`simulations/f81_violation_on_f113_fits.py`](../simulations/f81_violation_on_f113_fits.py) (leg 1), [`simulations/f81_identity_velocity_grounding.py`](../simulations/f81_identity_velocity_grounding.py) (leg 2), [`simulations/f84_net_cooling_readout_marrakesh.py`](../simulations/f84_net_cooling_readout_marrakesh.py) (leg 3)
**Data:** [`data/ibm_f95_angle_steering_may2026/`](../data/ibm_f95_angle_steering_may2026/) (Kingston, 2026-05-16), [`data/ibm_price_pair_july2026/`](../data/ibm_price_pair_july2026/) (Marrakesh, 2026-07-04, runs 1+2)
**Builds on:** F81/F82/F84 ([`PROOF_F81_PI_CONJUGATION_OF_M.md`](../docs/proofs/PROOF_F81_PI_CONJUGATION_OF_M.md), [`PROOF_F82_T1_DISSIPATOR_CORRECTION.md`](../docs/proofs/PROOF_F82_T1_DISSIPATOR_CORRECTION.md), [`PROOF_F84_AMPLITUDE_DAMPING.md`](../docs/proofs/PROOF_F84_AMPLITUDE_DAMPING.md)), F113 ([`F113_T1_EXTRACTION_KINGSTON.md`](F113_T1_EXTRACTION_KINGSTON.md)), the price_pair campaign ([`PRICE_PAIR_HARDWARE_PREDICTION.md`](PRICE_PAIR_HARDWARE_PREDICTION.md))

## The problem

The theorem chain F81 → F82 → F84 makes `f81_violation = ‖M_anti − L_{H_odd}‖_F` a one-scalar readout of the population-inverting (non-unital) noise content: every Pauli channel (X/Y/Z dephasing, depolarizing, correlated ZZ-dephasing) contributes exactly zero; only σ⁻/σ⁺ content survives, with the closed form `√(Σ_l (γ↓−γ↑)²_l) · 2^(N−1)`, Hamiltonian-independent and γ_z-independent. But the discriminator's input is the full 4^N×4^N superoperator L. In simulation we build L; on hardware nobody hands it to you. This document is the bridge: what the number IS physically, and how to measure it without building L.

**The tautology trap** (the reason a naive "compute it on a fitted L" is empty): if the fit is Z+T1-parameterized, the violation of the fitted L equals the F82 closed form of the fitted γ_T1 by construction; you learn nothing beyond the fit parameter. Leg 1 demonstrates this exactly; leg 3 escapes it.

## Leg 1: the method demo on the F113 Kingston fits (and the tautology, made visible)

The F113 pipeline ([`f113_t1_extraction_kingston.py`](../simulations/f113_t1_extraction_kingston.py)) fit a minimal Z+T1 Lindblad model (known Z-drive ω, free γ_z and γ_T1) to the four f95 pair-run trajectories. We refit (values reproduce the F113 table to ≤0.05%) and compute the violation on each fitted L directly.

A technical note that matters for reuse: the f95 Hamiltonian is the 1-body drive H = (ω/2)·Σ_l Z_l, which `pi_decompose_M`'s terms interface cannot express (1-letter terms fall through its 2-body/k-body branches). The script therefore hand-rolls the decomposition from the same framework primitives (`palindrome_residual`, `build_pi_full`, the Pauli-basis transform), with H_odd projected from below via X^⊗N conjugation, and GATES the hand-rolled path against `pi_decompose_M` on a terms-expressible case first: agreement is bit-exact (diff 0.0).

| Pair-run | γ_z fit | γ_T1 fit | violation (/µs) | violation / F82(γ_T1 fit) | γ_T1,RMS inverse | T1 implied | γ_T1 calib | overshoot |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| ω=0.13 A_mid q82-q83 | 0.0980 | 0.00574 | 1.623e-2 | **1.000000000000** | 0.00574 | 174.3 µs | 0.00507 | 1.13 |
| ω=0.13 B_high q13-q14 | 0.1402 | 0.00616 | 1.742e-2 | **1.000000000000** | 0.00616 | 162.4 µs | 0.00464 | 1.33 |
| ω=0.25 A_mid q82-q83 | 0.3361 | 0.00722 | 2.041e-2 | **1.000000000000** | 0.00722 | 138.6 µs | 0.00507 | 1.42 |
| ω=0.25 B_high q13-q14 | 3.5237 | 0.00564 | 1.595e-2 | **1.000000000000** | 0.00564 | 177.4 µs | 0.00464 | 1.22 |

Readings:

1. **The ratio column is the tautology, confirmed to all digits.** On a Z+T1-parameterized fit the violation is the closed form of the fitted γ_T1. The pipeline (hardware fit → violation → inverse readout) runs end-to-end and is exact; as a measurement it carries no information beyond γ_T1_fit.
2. **The physics that survives is the overshoot column:** the violation-implied γ_T1,RMS exceeds the independent calibration 1/T1 by 1.13-1.42, exactly the F113 underfit signature (the 2-parameter model dumps all non-T1 noise into the σ⁻ channel). The violation of a parameterized fit reads an EFFECTIVE non-unital content, an upper bound on the true σ± rates.
3. One honesty note on the F113 narrative: the overshoot is ω-monotone on A_mid (1.13 → 1.42) but not on B_high (1.33 → 1.22, the run whose γ_z fit degenerated to 3.5/µs). "Higher drive → more absorbed noise" holds as a trend, not per-pair.

## Leg 2: the grounding (what the number IS)

Strip the labels and the violation has a state-space meaning. The maximally mixed state I/2^N is the fixed point of every unital channel: all Pauli watching fixes it. Non-unital content gives the identity a velocity. [`f81_identity_velocity_grounding.py`](../simulations/f81_identity_velocity_grounding.py) confirms, by state propagation (no Pauli-basis entry inspection anywhere):

1. **v_l := d⟨Z_l⟩/dt at ρ = I/2^N equals γ↓,l − γ↑,l**, and on all four fitted Kingston L's `2^(N−1)·√(Σ_l v_l²)` reproduces the leg-1 violations to all printed digits. The discriminator IS the RMS identity-escape velocity, in the framework's Pauli-basis packaging.
2. **Temperature-independence at the identity (F84's vacuum reading), exact:** for γ↓ = 0.007+c, γ↑ = c the velocity stays 0.0070000000 for every c up to 0.05 (total rates 8× the net). The symmetric thermal traffic cancels; the velocity reads only the spontaneous (vacuum) component.
3. **The measurement recipe:** on the Bloch equation dz/dt = a − b·z (a = γ↓−γ↑, b = γ↓+γ↑ = 1/T1), the affine constant is a = z∞·b = z∞/T1. Both factors come out of the STANDARD T1 leg: the decay rate and the asymptote. Recovered from an exact ⟨Z⟩(t) trajectory to 8.7e-19.
4. **Support check** (the "single-entry lead" of [`PROOF_F82_T1_DISSIPATOR_CORRECTION.md`](../docs/proofs/PROOF_F82_T1_DISSIPATOR_CORRECTION.md) Step 3, confirmed on the fitted hardware L): at N=1 the difference M_anti − L_{H_odd} has exactly one entry, D[Z ← I] = γ_T1; at N=2 exactly the 8 site-l I→Z transfer entries (2 sites × 4 spectator strings), each of magnitude γ_T1. Nothing else.

Two corollaries worth stating precisely:

- **The number is not 1/T1.** It is the net one-way flux γ↓−γ↑ = z∞/T1, the piece that survives when the up and down thermal traffic cancels. "T1-vs-T2 discriminator" is the cold-qubit approximation (z∞ ≈ 1); leg 3 shows the distinction is 10-40% on real hardware, not academic.
- **ZZ-blindness on all three levels:** coherent ZZ is Π²-even Hamiltonian content (never enters M_anti); ZZ-dephasing is a Pauli channel (cancels exactly, F84 lemma); and in the recipe, diagonal Hamiltonians do not move populations, so the T1-leg trajectory is unaffected. The always-on ≈ −3.9 kHz ZZ that the price_pair campaign decoded is invisible to this readout, structurally.

## Leg 3: the first non-tautological readout (Marrakesh, existing data, zero QPU)

The price_pair campaign's Block B ([`PRICE_PAIR_HARDWARE_PREDICTION.md`](PRICE_PAIR_HARDWARE_PREDICTION.md), 2026-07-04) stored readout-mitigated ⟨Z_l⟩(t) trajectories of the |1⟩ free decay: 10 delays, 0-320 µs (≈ one full T1), 8192 shots, two independent 3-qubit lines. That is exactly the R2 recipe's input. No Lindblad channel fit anywhere: per qubit, a free-asymptote exponential fit gives (z∞, b), and a = z∞·b is the measured identity-escape velocity. [`f84_net_cooling_readout_marrakesh.py`](../simulations/f84_net_cooling_readout_marrakesh.py):

| Line | Qubit | z∞ | b (/µs) | a = z∞·b (/µs) | naive 1/T1 | χ²/dof |
|---|---|---:|---:|---:|---:|---:|
| run 1 [2,3,4] | q2 | +0.901 ± 0.082 | 0.00316 | **0.00285 ± 0.00008** | 0.00316 | 0.39 |
| | q3 | +0.866 ± 0.075 | 0.00333 | **0.00288 ± 0.00008** | 0.00333 | 0.87 |
| | q4 | +0.941 ± 0.031 | 0.00541 | **0.00509 ± 0.00005** | 0.00541 | 0.59 |
| run 2 [93,94,95] | q93 | +0.697 ± 0.051 | 0.00408 | 0.00284 ± 0.00007 | 0.00408 | 0.78 |
| | q94 | +0.576 ± 0.017 | 0.00778 | (0.00448 ± 0.00005) | 0.00778 | **12.49** |
| | q95 | +0.754 ± 0.050 | 0.00413 | 0.00311 ± 0.00006 | 0.00413 | 1.32 |

Packaged (N=3 line, factor 2^(N−1) = 4):

- **run 1: f81_violation(local σ± family) = 0.02603 ± 0.00024 /µs**, identity-escape RMS velocity 0.00376/µs (T1-equivalent 266 µs). The naive all-z∞=1 packaging would read 0.02839: T1 relaxometry without the asymptote overestimates the palindrome-breaking content by 8-14% on this line.
- run 2: 0.02460 ± 0.00024 /µs, RMS velocity 0.00355/µs, **flagged**: see below.

Readings:

1. **This is the bridge.** The number that was an operator diagnostic is now two decay parameters per qubit, from a standard relaxometry leg with asymptote fitting. It is non-tautological in the operative sense: no Lindblad channel set was chosen anywhere; the inputs are directly measured decay constants, and the packaging is a theorem (F84), not a fit.
2. **The recipe self-diagnoses.** q94's single-exponential fit fails loudly (χ²/dof = 12.5; early decay ≈ 3× faster than late, the anomalous-T1 qubit of run 2). Its z∞ is a fit parameter, not a measurement, so the run-2 line number is flagged rather than trusted. A clean χ² is part of the readout.
3. **The asymptotes sit measurably below +1** on clean fits: z∞ = 0.87-0.94 (run 1) and 0.70-0.75 (run 2 clean qubits), i.e. steady-state excited populations of 3-7% and 12-15%. The stored per-qubit readout errors (0.2-1.3%, already mitigated) are an order of magnitude too small to produce this, so it is real steady-state physics of the idle qubit (bath heating γ↑ > 0) or a slow systematic in the leg; one leg alone cannot separate the two. This is exactly the F84 distinction at work: the discriminator wants the net flux, and on this device the net flux is NOT 1/T1.
4. Cross-line consistency: the two lines agree on the RMS velocity (0.00376 vs 0.00355 /µs) far better than on their T1 profiles (185-316 vs 129-245 µs). One session, one device; noted as an observation, not a claim.

## Limits

- **The R2 recipe reads the local σ± family only.** It is a lower bound on the full violation: correlated or non-local Π²-odd channels live on other entries of M_anti − L_{H_odd} and need the free-form route (R1: an unconstrained Lindblad fit on an idle tomography grid with tens of time points; the f95 6-point trajectories are too thin, the F113 lesson).
- **The recipe needs a diagonal idle Hamiltonian** (detuning + ZZ are fine; an exchange term would mix site populations and break the per-site Bloch form). Superconducting idle lines satisfy this; other platforms must check.
- **The z∞ attribution is open.** Genuine bath heating and a slow leg systematic both depress the asymptote; the discriminating experiment is the |0⟩ heating leg (below), not run.
- One session per line, asymptotes extrapolated from t_max ≈ 1·T1; this is a method demonstration with honest error bars, not a device characterization.

## The missing half: the heating leg (simulator-validated 2026-07-05, pre-registered, not yet flown)

Prepare |0⟩ instead of |1⟩ beside the |1⟩ leg, same 10-delay grid to 320 µs, both legs in ONE job, measure ⟨Z⟩: under a thermal bath both legs relax to the SAME z∞ with the SAME rate b, so the meeting test |z∞_down − z∞_up| ≤ 2σ separates bath heating from a leg systematic, and the joint fit (shared z∞ and b, free per-leg z0) separates γ↓ = b·(1+z∞)/2 and γ↑ = b·(1−z∞)/2 per qubit. Runner: `run_heating_leg.py` in the external tomography pipeline beside `run_price_pair.py`; 22 circuits (2 legs × 10 delays + 2 readout calibrations) × 8192 shots ≈ 1.4 QPU min (anchor: 76 circuits ≈ 4.7 min).

**Simulator validation** (qiskit-aer 0.17.2: the actual circuits through AerSimulator with `RelaxationNoisePass` thermal relaxation incl. excited-state population on the delays, per-qubit asymmetric readout error, tensor-product mitigation; artifacts + full numbers in [`data/ibm_heating_leg_july2026/`](../data/ibm_heating_leg_july2026/)):

1. Noiseless parity: legs flat, all rates ≈ 0, violation 0.00000 ± 0.00010; the meeting test correctly reports N.A. where no decay is resolved (model-free guard).
2. Warm scenario (run-2's Block-B reading taken at face value: T1 = 245/129/242 µs, p_exc = 0.152/0.212/0.123): MEET on all three qubits, every planted (γ↓, γ↑, a) recovered within 2σ, violation 0.02469 ± 0.00017 /µs, matching the run-2 hardware reading 0.02460.
3. Cold scenario (p_exc = 0.01): MEET, planted recovered, γ↑ resolved at ≈ 3σ even at 1% excited population; violation 0.03794 ± 0.00023 /µs (higher than warm: heating reduces the net flux, exactly F84).
4. Planted TLS control (q1's T1 a 50/50 mixture of 60/250 µs, the q94 signature): the qubit trips BOTH tripwires (meeting test SPLIT and joint χ²/dof = 11.05 against the flag threshold 4); the clean qubits stay clean.
5. 20-seed scatter calibration: quoted fit errors are slightly conservative (scatter/σ ratios 0.59-0.86), the safe side of the price_pair fit-σ lesson.

**Pre-registered predictions for the hardware run** (line [93, 94, 95] preferred, same as price_pair run 2; stated before any shot):

- **P1 (meeting):** every χ²-clean qubit MEETs at 2σ. A SPLIT on a χ²-clean qubit = leg systematic; that qubit's R2 packaging is withheld and γ↑ is quoted from the up-leg alone.
- **P2 (the attribution):** if run 2's z∞ readings were bath physics, γ↑ lands at ≈ 0.0005-0.0017 /µs per qubit (≥ 10σ above zero at the validated error bars ± 0.00002-0.00004); if they were a leg systematic, γ↑ < 0.0002 /µs and the up-leg stays flat within ≈ 3% absolute. The two hypotheses are separated by ≫ 5σ either way; either answer closes the attribution.
- **P3 (q94):** if its non-exponentiality is a stable device property, the χ² flag (> 4) and/or SPLIT recurs on q94; a clean q94 means the 2026-07-04 shape was transient (TLS switching), and its two-leg numbers replace the flagged Block-B reading.
- **Instrument conditions, declared now:** χ²/dof > 4 on any fit = that qubit's numbers flagged, not trusted; unresolved decay = meeting test N.A.; fresh same-session calibration required (the stale-calibration trap, 61% vs 1.9%).

Confirmations-grade on a P2 pass. Gated on Tom's go per the standing QPU rule. **Queue note (2026-07-05):** all three available IBM backends showed > 30k pending jobs (suspected stuck external submit loop); nothing was submitted, this session ran simulator-only.

## Where this lands (the S4 socket)

This bridge is the missing piece named in the outbound arc (`outbound_label_adapters`, S4: the watching → quantum noise spectroscopy). What the adapter can now hand over, in their objects: a one-scalar, closed-form readout of the non-unital noise content, exactly zero for every unital (Pauli) channel including correlated dephasing, equal to the RMS spontaneous-emission (net cooling) rate, temperature-independent, Hamiltonian-blind, ZZ-blind, and measurable by standard T1 relaxometry with asymptote fitting. The measured Marrakesh number above is the demonstration. The adapter itself (genre: the three sisters in [`docs/outbound/`](../docs/outbound/)) is the next step of the arc, not this document.

## Reproduction

```
python -X utf8 simulations/f81_violation_on_f113_fits.py        # leg 1, ~30 s (4 Nelder-Mead refits)
python -X utf8 simulations/f81_identity_velocity_grounding.py   # leg 2, ~5 s
python -X utf8 simulations/f84_net_cooling_readout_marrakesh.py # leg 3, ~2 s
```
