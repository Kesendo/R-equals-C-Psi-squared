# IBM single-excitation population handover (Kingston, 2026-05-31)

`ibm_ep_onset_may2026` is the historical run identifier, not the current physical verdict.

Single-excitation coherent walk on a 3-site chain (Q13-Q14-Q15) of ibm_kingston, the hardware
reading of [`experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md`](../../experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md)
(a finite-time population record beside the ideal model's separate 1/N stationary-sector theorem).

- Prepare |100>, Trotter gates RXX(J dt) RYY(J dt), corresponding to H = (J/2)·(XX+YY)
  per bond (runner parameter J = 1.5 rad/us, dt = 0.5 us),
  read per-site populations <n_l> (Z-basis, no tomography), 4096 shots, 21 time points to 20 us.
- Runner: `run_ep_onset.py` (--simulate / --hardware / --analyze) in the IBM pipeline.
- Chain verified on same-day calibration: 13-14-15 a connected path, all operational,
  T2 = 302 / 235 / 176 us.

## Result (Job d8dr7dfd0j8c73f4man0)

The excitation sloshes 0 → 2 → 1 → 0 with a ~3 us period; the site-0 revival
falls 0.84 → 0.43 over 15 us. At 20 us the three measured Z-basis marginals are
0.339 / 0.426 / 0.338 (sum 1.103). They have moved closer together, but are not a normalized
one-excitation distribution and do not establish convergence or an asymptotic 1/N fixed point.
transfer = 0.773, swing = 0.687.

This is one finite-time population trajectory on a real chip, populations only. The ideal
number-conserving connected-chain model's 1/N stationary-sector limit is a separate theorem.

**Honest caveat.** The revivals fade faster than the T1+T2-only simulation predicted (0.43 vs
~0.83 at t = 18 us). That discrepancy is consistent with accumulated circuit/gate cost, but
the stored record contains no calibrated gate/readout/leakage model that identifies a unique cause.

## Part B — injected-dephasing population scan (Job d8drjbfd0j8c73f4mobg)

The historical runner scanned injected dephasing (random-Z twirl, K=16 instances). Its random
phase variance σ²=2Γdt makes coherences decay as exp(−Γt), so its labels are Q_label=J/Γ.
The repository's canonical Lindblad jump √γ Z makes coherences decay as exp(−2γt): γ=Γ/2 and
**Q_Lindblad = 2 Q_label**. The revival (max <n_0> for t >= 2 us) is compared with the 1/N reference level:

| Q_label | Q_Lindblad | revival | population reading |
|---|---|---|---|
| 0.5 | 1 | 0.30 | near 1/N |
| 1.0 | 2 | 0.36 | near 1/N |
| 1.5 | 3 | 0.34 | near 1/N |
| 2.5 | 5 | 0.49 | larger return |
| 5.0 | 10 | 0.56 | larger return |
| 20 | 40 | 0.70 | larger return |

The sampled population handover lies between Q_label=1.5 and 2.5, equivalently canonical
Q_Lindblad=3 and 5. Its spectral character remains open: populations do not certify critical
damping, an EP, eigenmode coalescence, or Jordan structure. The result tracks the
validated exact-statevector twirl simulation (0.28 → 0.84 across the same scan). At Q_label=20
the measured 0.703 is below the exact-twirl simulation 0.842 (raw 0.8417853730254796); the difference is consistent with accumulated
circuit/gate cost, but is not uniquely attributed by these records. The 1/N line is a reference,
not a measured floor; the sampled change is not a certified onset or threshold.

## Files
- `ep_onset_hardware_ibm_kingston_20260531_060943.json` , Part A hardware (coherent walk, per-site populations, job id)
- `ep_onset_simulate_20260531_060202.json` , Part A Aer simulate (ideal + Kingston T1+T2), same-day calibration
- `ep_onset_hardware_ep_ibm_kingston_20260531_064022.json` — Part B hardware (historically named twirl scan, revival vs Q_label, job id)
- `ep_onset_simulate_twirl_20260531_063048.json` — Part B twirl validation (K=16, exact statevector), confirms the injected-dephasing population curve
