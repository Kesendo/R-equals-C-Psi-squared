# IBM EP Onset (Kingston, 2026-05-31)

Single-excitation coherent walk on a 3-site chain (Q13-Q14-Q15) of ibm_kingston, the hardware
reading of [`experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md`](../../experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md)
(the post-EP flow: memory sloshing → forgetting → 1/N target).

- Prepare |100>, Trotter evolution of H = J·(XX+YY) per bond (J = 1.5 rad/us, dt = 0.5 us),
  read per-site populations <n_l> (Z-basis, no tomography), 4096 shots, 21 time points to 20 us.
- Runner: `run_ep_onset.py` (--simulate / --hardware / --analyze) in the IBM pipeline.
- Chain verified on same-day calibration: 13-14-15 a connected path, all operational,
  T2 = 302 / 235 / 176 us.

## Result (Job d8dr7dfd0j8c73f4man0)

The excitation sloshes 0 → 2 → 1 → 0 with a ~3 us period (the reborn memory); the site-0 revival
fades monotonically 0.84 → 0.43 over 15 us (the forgetting); the populations converge toward
1/3 = 1/N at 20 us (the target: 0.34 / 0.43 / 0.34). transfer = 0.773, swing = 0.687.

The whole arc of the experiment doc (memory born → fades → 1/N fixed point) in one population
trajectory, on a real chip, populations only.

**Honest caveat.** The revivals fade faster than the T1+T2 simulate predicted (0.43 vs ~0.83 at
t = 18 us). The cause is two-qubit gate error (~160 RZZ gates by 20 us at ~0.5% each), not
dephasing; the decay envelope is Trotterization-limited (~9 us), not T2-limited (~200 us). The
sloshing and the flow into 1/N are clean; only the decay rate is gate-cost, not physics.

## Part B , the EP onset (Job d8drjbfd0j8c73f4mobg)

Scanning injected dephasing (random-Z twirl, K=16 instances) to push Q = J/γ down through
Q_EP ~ 1.5. The revival (max <n_0> for t >= 2 us, the memory's return) collapses to the
equipartition floor:

| Q | revival | regime |
|---|---|---|
| 0.5 | 0.30 | overdamped (forgotten, ~1/3) |
| 1.0 | 0.36 | overdamped |
| 1.5 | 0.34 | ~EP |
| 2.5 | 0.49 | sloshing lifts off |
| 5.0 | 0.56 | memory returning |
| 20 | 0.70 | memory present |

The memory switches on as Q crosses ~1.5 to 2.5: the EP onset on a real chip. Tracks the
validated twirl simulate (0.31 → 0.84 across the same scan); the high-Q side is suppressed
(0.70 vs 0.84 at Q=20) by gate error, but the floor and the onset are clean. The RZ twirl gates
are virtual on IBM (error-free); only the RXX/RYY carry gate cost.

## Files
- `ep_onset_hardware_ibm_kingston_20260531_060943.json` , Part A hardware (coherent walk, per-site populations, job id)
- `ep_onset_simulate_20260531_060202.json` , Part A Aer simulate (ideal + Kingston T1+T2), same-day calibration
- `ep_onset_hardware_ep_ibm_kingston_20260531_064022.json` , Part B hardware (EP-onset twirl scan, revival vs Q, job id)
- `ep_onset_simulate_twirl_20260531_063048.json` , Part B twirl validation (K=16, exact statevector), confirms the twirl reproduces the dephasing EP onset
