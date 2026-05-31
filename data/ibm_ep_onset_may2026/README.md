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

## Files
- `ep_onset_hardware_ibm_kingston_20260531_060943.json` , hardware run (per-site populations, metrics, job id)
- `ep_onset_simulate_20260531_060202.json` , Aer simulate (ideal + Kingston T1+T2 walk, plus the EP-onset dephasing scan), same-day calibration
