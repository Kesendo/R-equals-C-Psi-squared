# Heating-Leg Protocol: Simulator Validation, July 2026

Simulator validation AND the hardware run of the two-legged T1 protocol (cooling leg |111⟩ + heating leg |000⟩, one job). The validation ran on the IBM simulator (qiskit-aer 0.17.2, AerSimulator + `RelaxationNoisePass` with per-qubit T1/T2/excited-state population on the delay instructions + per-qubit asymmetric readout error + tensor-product mitigation) during the morning's > 30k-queue outage; the hardware run flew the same day when the Kingston queue emptied (0 pending, fresh calibration 08:28:50Z, Tom's go). Pre-registration, amendment, and the HARDWARE RECORD with verdicts live in [`experiments/F81_VIOLATION_HARDWARE_BRIDGE.md`](../../experiments/F81_VIOLATION_HARDWARE_BRIDGE.md).

Purpose: the protocol is the "missing half" of the f81 hardware bridge. The price_pair Block B (|1⟩ decay) found asymptotes measurably below +1; one leg cannot attribute this (bath heating γ↑ > 0 vs a slow leg systematic). Both legs relax to the SAME z∞ with the SAME rate under a thermal bath, so the pre-registered meeting test |z∞_down − z∞_up| ≤ 2σ separates the hypotheses, and the joint fit separates γ↓ = b(1+z∞)/2 from γ↑ = b(1−z∞)/2 per qubit.

## Runner

`run_heating_leg.py` in the external tomography pipeline (`D:\Entwicklung\Projekte\.NET Projekte\AIEvolution\AIEvolution.UI\experiments\ibm_quantum_tomography\`, beside `run_price_pair.py`). 22 circuits: 2 legs × 10 delays (0-320 µs, the price_pair Block-B grid) + 2 readout calibrations; 8192 shots. Hardware cost anchored to the price_pair campaign (76 circuits ≈ 4.7 QPU min): ≈ 1.4 QPU min.

## What is in this directory

| File | Mode | Result |
|------|------|--------|
| `heating_leg_aer_parity_*.json` | actual circuits, noiseless Aer | legs flat at ±1; all rates ≈ 0; meeting test correctly reports N.A. (decay unresolved); violation 0.00000 ± 0.00010 |
| `heating_leg_aer_warm_*.json` | thermal noise, warm scenario (run-2-like: T1 = 245/129/242 µs, p_exc = 0.152/0.212/0.123) | MEET on all 3 qubits; all planted (γ↓, γ↑, a) recovered within 2σ; violation 0.02469 ± 0.00017 /µs (the run-2 hardware reading was 0.02460) |
| `heating_leg_aer_cold_*.json` | thermal noise, cold scenario (same T1/T2, p_exc = 0.01) | MEET on all 3; planted recovered within 2σ; 1% excited population separates γ↑ from zero at 1.6-3.2σ per qubit; violation 0.03794 ± 0.00023 /µs (higher than warm: heating reduces the net flux, F84) |
| `heating_leg_aer_tls_*.json` | q1's T1 planted as a 50/50 mixture of 60/250 µs (the q94 signature) | q1 trips BOTH tripwires: meeting test SPLIT (+ joint χ²/dof = 11.05, flag threshold 4); q0/q2 stay clean (MEET, no flag) |
| `heating_leg_seeds_*.json` | 20-seed scatter calibration, warm | seed-to-seed scatter of a vs mean fit σ: ratios 0.59/0.86/0.80, i.e. the quoted errors are slightly conservative (the safe side of the price_pair fit-σ lesson) |
| `heating_leg_ibm_kingston_20260705_105648.json` | **HARDWARE** (ibm_kingston, job d951mhkql68s73ca3u0g, path [82, 83, 13], ≈ 1.6 QPU min) | **SPLIT on all three qubits** (13.7σ/5.1σ/6.6σ): up-legs flat at +0.98-1.00 (γ↑ = 2-7·10⁻⁵ /µs, p_th ≲ 1%, cold bath); down-legs still rising at 320 µs, non-exponential (pinned-asymptote two-rate mixture accepted, single exponential rejected). Attribution closed: the sub-unity down-leg asymptote is a leg systematic, not thermal population. Verdicts + post-hoc analysis in the bridge doc's HARDWARE RECORD; typed Confirmation `f84_heating_leg_attribution_kingston_july2026`. |

## Reading

The instrument does what the pre-registration needs: it separates warm from cold baths at ≫ 5σ (γ↑ errors ± 0.00002-0.00004 /µs vs hypothesis gap ≥ 0.0005 /µs), the meeting test passes exactly when the physics is a bath and reports N.A. when no decay is resolved, and a fluctuating-T1 qubit cannot silently launder into a clean number: it is flagged twice (SPLIT + χ²). Fit errors are validated conservative by the seed sweep.

The hardware pre-registration (predictions P1-P3, verdict bands, instrument conditions) lives in the bridge document, added 2026-07-05 after this validation.
