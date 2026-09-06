# Gamma Sets a Dissipative Scale; It Does Not Define Time

<!-- Keywords: gamma dephasing timescale, parameter time, finite trajectory
comparison, tau gamma t scaling, no experienced-time ontology -->

**Status:** Computed finite trajectory comparisons (Tier 2); no time ontology
**Date:** March 22, 2026; current-state rewrite September 4, 2026
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Evidence:** [two_qubits_no_noise.py](../simulations/two_qubits_no_noise.py),
[gamma_unit_scaling_gate.py](../simulations/gamma_unit_scaling_gate.py)

---

## Result

In the declared Lindblad models, `t` is the evolution parameter, `J` sets a
coherent rate, and `gamma` sets a dissipative rate. The dimensionless
trajectory depends on ratios such as `Q = J/gamma`; changing gamma while
holding J fixed changes both the scale and the trajectory's shape. Gamma is
therefore neither time itself nor a sufficient description of a clock.

This document makes no claim about experienced time, a fundamental arrow of
time, or the microscopic origin of gamma.

## What was actually computed

The finite two-qubit comparison used two preparations:

- At `gamma = 0`, the tested Bell-plus preparation is stationary under the
  chosen Hamiltonian, while the tested `|01>` preparation evolves unitarily
  and recurrently.
- At positive gamma, those same finite trajectories show damping in the
  declared Markovian dephasing model.
- Matching only `tau = gamma*t` while holding J fixed does not collapse all
  observables, because this also changes `Q = J/gamma`.
- Under the joint rescaling `L(lambda J, lambda gamma) = lambda L(J, gamma)`,
  with Q held fixed and time rescaled inversely, the compared trajectories
  collapse to numerical precision.

These observations distinguish stationary, recurrent, and damped trajectories
in a small declared model. They are not necessary or sufficient conditions
for experience, memory, or a thermodynamic arrow.

## What the old argument got wrong

Earlier versions promoted the finite comparison into the statements "gamma is
time" and "without gamma there is no time." Neither follows. Hamiltonian
dynamics still uses `t` at `gamma = 0`; a stationary preparation says something
about that preparation, not about whether time exists; and dissipative damping
does not by itself define experience.

The former producer
[gamma_is_time_proof.py](../simulations/gamma_is_time_proof.py) and its
[committed output](../simulations/results/gamma_is_time_proof.txt) are retained
as withdrawal-safe tombstones. They do not supply trajectory evidence.

## Precise language

| Statement | Status |
|---|---|
| Gamma sets a dissipative timescale in the declared generator | Established by the model definition |
| The tested positive-gamma trajectories damp | Finite numerical observation |
| The tested zero-gamma trajectories are stationary or recurrent | Finite numerical observation for two preparations |
| `tau = gamma*t` universally determines the trajectory | False when other dimensionless ratios change |
| Gamma defines experienced time or its origin | Not established |
| A nonzero dissipator identifies its microscopic source | Not established; see the [Incompleteness Proof](proofs/INCOMPLETENESS_PROOF.md) |

## Reproduction

Run:

```bash
python simulations/two_qubits_no_noise.py
python simulations/gamma_unit_scaling_gate.py
```

Interpret the resulting numbers only within their stated system, preparation,
time window, and observable set.
