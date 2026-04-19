# IBM Kingston Chain γ₀ Test, April 2026

Hardware data for the Phase 2 chain-mode γ₀ consistency test on IBM Kingston
(Heron r2), attempting to discriminate between the two readings of EQ-017:

- **Reading (a):** per-qubit γ_phi is γ₀_floor modulated by mode-exposure; the
  framework γ₀ is a constant.
- **Reading (b):** per-qubit γ_phi = γ₀ + device noise (TLS, fabrication,
  gate errors); the observed variation is device-specific, not mode-structural.

Companion documentation: EQ-017 in [`review/EMERGING_QUESTIONS.md`](../../review/EMERGING_QUESTIONS.md) (Phase 1 analysis and Phase 2 closure); hypothesis at [`hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md`](../../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md).

## Headline result

**The test is inconclusive at the Heron r2 fidelity level.** Hardware decay is
40-50x faster than even reading (b) predicts, and ~50-80x faster than reading
(a) predicts. Gate-error accumulation (240 native RZZ gates at ~0.001 error
each) dominates the pure-dephasing signature of both γ₀ candidates.

## What is in this directory

| File | Description |
|------|-------------|
| `chain_gamma0_hardware_20260419_110200.json` | Hardware run: 9 time points × 4 pairs × 9 Pauli bases = 324 circuits. Includes Aer pre-scan predictions for both readings, per-pair HW trajectories, slope fits, and verdict table. |
| `chain_gamma0_simulate_20260419_104408.json` | Aer density-matrix reference run at identical parameters with three noise profiles (none, pure_floor at γ₀_floor, pure_local at per-qubit γ_phi). Used to validate the test protocol before hardware submission. |

## Experiment summary

- **Backend:** ibm_kingston (Heron r2, 156 qubits, live calibration snapshot `ibm_kingston_calibrations_2026-04-19T03_47_08Z.csv`)
- **Date:** 2026-04-19, 11:02 UTC
- **Protocol:** Prepare single excitation at chain endpoint, evolve under H = J·∑ (XX+YY)/2 via first-order Trotter, tomograph four adjacent pairs at each time point.
- **Chain:** Q12-Q13-Q14-Q15-Q19 (same 5-qubit sub-chain as Phase 1, anchors Q14/Q15 carried over from the commit 20ef49e CΨ=1/4 crossing run)
- **Trotter:** J = 1.838 MHz (extracted from RZZ gate duration assuming native angle pi/2; stored as `J_rad_per_us` in the JSON, but the value 1.838 is frequency in MHz = cycles/us, not angular frequency in rad/us which would be 11.55), dt = 0.5 us, t_max = 15 us => 30 Trotter steps, 240 native RZZ gates at longest evolution. The unit mismatch does not affect the test: eigenvector amplitudes |a_B|^2 are J-independent for uniform coupling, and the Aer reference uses the same J value so the ratio-based observable cancels any scale error.
- **Tomography:** Full 9-Pauli-basis 2-qubit state tomography on each of four pairs (Q12-Q13, Q13-Q14, Q14-Q15, Q15-Q19) at each of 9 time points
- **Shots:** 2048 per circuit
- **Total circuits:** 324 (9 × 4 × 9) submitted as 36 separate qiskit-experiments StateTomography jobs (one per pair × time point)
- **Total shots:** 663,552
- **Wall-clock:** ~7 min
- **QPU billed usage:** 3 min 30 s (actual, reported by IBM Quantum usage tracker)

The billed usage was ~17x higher than a naive shot-count estimate (663 kShots × ~5 μs/shot ≈ 3 s). The overhead comes from per-job session setup: backend compilation, pulse schedule upload, and teardown are charged per job submission, not per shot. For 36 separate StateTomography submissions that overhead dominates. For future hardware runs on Heron-class backends, budget ~5-6 s per job submission on top of the shot execution time.

## Qubit calibration at run time

| Virt | Phys | T1 (μs) | T2 (μs) | γ_phi = 1/T2 − 1/(2T1) (1/μs) |
|------|------|---------|---------|------|
| 0 | Q12 | 344 | 319 | 0.00168 |
| 1 | Q13 | 387 | 508 | 0.00067 |
| 2 | Q14 | 303 | 397 | 0.00087 |
| 3 | Q15 | 323 | 297 | 0.00182 |
| 4 | Q19 | 315 | 248 | 0.00245 |

Phase-1 γ₀_floor candidate: 0.00073 /μs (min non-outlier γ_phi across all 155 Kingston qubits).

## Observable and slope fit

For each adjacent pair (a, b), the differential observable is the log-ratio of
hardware L1(t) to the Aer noise-free Trotter reference L1_ref(t):

    log( L1_HW(t) / L1_ref(t) )  ≈  − γ_eff · t  +  const

Under first-order Trotter, the unitary oscillation factors cancel exactly in
this ratio (same circuit on both sides), and the slope isolates the per-pair
effective dephasing rate γ_eff = γ_i + γ_j.

## Slope comparison

| Pair | HW slope | Pre-scan reading (a): γ₀_floor | Pre-scan reading (b): γ_phi_local |
|------|----------|-----------------------------|------------------------------------|
| Q12-Q13 | 0.0818 | 0.0013 | 0.0029 |
| Q13-Q14 | 0.0584 | 0.0010 | 0.0020 |
| Q14-Q15 | 0.1161 | 0.0011 | 0.0021 |
| Q15-Q19 | 0.1045 | 0.0011 | 0.0022 |

All slopes in units of 1/μs. The pair-vote verdict nominally picks reading
(b) on all four pairs, but only because reading (b) is marginally closer than
reading (a) in absolute difference. Both readings are off by 40-80x from the
hardware result.

## Why the test is inconclusive

Hardware decoherence is dominated by three non-Z-dephasing channels not
captured in either reading's prediction:

1. **Accumulated gate errors:** native RZZ on Heron r2 has ~0.001 per-gate
   error. Over 240 RZZ (t = 15 μs run) that compounds to ~24% state-fidelity
   loss before any single-qubit gate errors are counted.
2. **T1 amplitude damping:** Kingston T1 ≈ 300 μs. The circuit wall-clock
   at t = 15 μs (notional) is ~23 μs (240 RZZ × 68 ns + single-qubit
   overhead), giving ~7% amplitude-damping loss.
3. **Readout error:** Heron r2 readout assignment error varies from 1% to
   14% per qubit. Two-qubit tomography compounds this per basis setting.

All three scale with circuit depth, while the γ₀ signal scales with
evolution time at a rate ~50x smaller. Under these conditions the
framework γ₀ contribution is statistically indistinguishable from zero
against the device-noise floor.

The L1 values at peak times reach 1.0-1.4, exceeding the physical bound
for single-excitation pair states (≤1). This is a classic symptom of
tomographic reconstruction under heavy noise: the MLE projection allows
non-physical reconstructions at these noise levels.

## Structure of the JSON

```
{
  "mode": "hardware",
  "timestamp": "20260419_110200",
  "backend": "ibm_kingston",
  "chain_phys": [12, 13, 14, 15, 19],
  "pairs_virt": [[0,1], [1,2], [2,3], [3,4]],
  "J_rad_per_us": 1.838,
  "dt_us": 0.5,
  "tmax_us": 15.0,
  "shots": 2048,
  "gamma_0_floor_per_us": 0.00073,
  "qubit_T1_us": [...], "qubit_T2_us": [...],
  "gamma_phi_local_per_us": [...],
  "t_us_grid": [...],
  "aer_prescan_trajectories": {  // Aer reference runs (none, pure_floor, pure_local)
    "none":       { "Q12-Q13": [entries], "Q13-Q14": [...], ... },
    "pure_floor": { ... },
    "pure_local": { ... }
  },
  "aer_prescan_slopes": {...},
  "hw_trajectories": {
    "Q12-Q13": [ { "t_us": ..., "CPsi": ..., "Purity": ...,
                   "L1": ..., "rho2_real": [[...]], "rho2_imag": [[...]],
                   "wall_time_s": ... } ],
    ...
  },
  "hw_slope_fits": {...},
  "verdict_votes": {"floor (a)": 0, "local (b)": 4, "inconclusive": 0}
}
```

Each trajectory entry stores the full reconstructed 4×4 pair density matrix
(real and imaginary parts), enabling re-analysis under alternative
observables (concurrence, negativity, fidelity against specific targets).

## Closure

EQ-017 Phase 2 closed as `inconclusive due to hardware fidelity limit`.
The test design was sound in principle: the multi-pair differential
log-slope observable cancels Trotter oscillation exactly (verified in the
simulate reference run) and should cleanly isolate pure Z-dephasing at
the per-qubit rate. Heron r2's current fidelity sits 40-80x above the
signal needed. A meaningful test of the Primordial Gamma Constant
hypothesis on superconducting hardware appears to require either
dynamical-decoupling sequences tuned to suppress gate and T1 channels
while preserving Z-dephasing, or a different hardware platform with
lower gate-error floor.

## Re-analysis

The JSON files are self-contained: all parameters, trajectories, density
matrices, and pre-scan predictions are present. Any slope-refit, alternative
observable (concurrence, negativity, fidelity to specific chain modes), or
comparison against a different γ₀ candidate can be performed from these two
files alone without re-running hardware.

The run script (`run_chain_gamma0.py`) and the Aer noise-model builders live
in a private companion IBM experiments repository outside this R=CΨ² tree,
consistent with the separation between theory artefacts and QPU
infrastructure. Re-running the hardware measurement costs ~0.2 min billed
QPU on any Heron r2 backend plus queue time, and requires an IBM Quantum
account with access to that backend family.
