# Gamma Sets a Dissipative Scale; It Does Not Define Time

<!-- Keywords: gamma dephasing timescale, parameter time, finite trajectory
comparison, tau gamma t scaling, no experienced-time ontology -->

**Status:** Computed finite trajectory comparisons (Tier 2); no time ontology
**Date:** March 22, 2026
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Evidence:** [gamma_is_time_proof.py](../simulations/gamma_is_time_proof.py)
([output](../simulations/results/gamma_is_time_proof.txt)),
[two_qubits_no_noise.py](../simulations/two_qubits_no_noise.py),
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

## What gamma does to the trajectory

Two preparations on the N=2 Heisenberg chain, gamma = 0 against gamma = 0.05,
over t in [0, 50]:

| Preparation | gamma | trace distance returns | CPsi crossings down/up | \|\|drho/dt\|\| at t=50 |
|---|---|---|---|---|
| `\|01>` | 0.00 | yes | 63/64 | 2.809609 |
| `\|01>` | 0.05 | no | 2/2 | 0.019149 |
| Bell+ | 0.00 | yes | 0/0 | 0.000000 |
| Bell+ | 0.05 | no | 1/0 | 0.000006 |

At gamma = 0 the `|01>` trajectory is exactly recurrent. `H` has spectrum
`{-3J, J}`, a single gap `4J`, so the period is `pi/(2J)` and the returns over
this window are set by that period rather than by the sampling: a coarser grid
resolves fewer of the same dips. At gamma = 0.05 the returns are gone and the
state has all but stopped. What the dephasing removes is the recurrence. The
Bell-plus preparation is stationary under this Hamiltonian at gamma = 0, so it
reports on that preparation rather than on whether anything is moving.

## The tau collapse, and what it actually measures

The generator is homogeneous: `L(J, gamma)*t = tau * L(J/gamma, 1)` with
`tau = gamma*t`. So every observable is a function of `tau` **and** of
`Q = J/gamma`, and a collapse test has to be run twice or it says nothing.
Spread is the largest gap between two curves at equal tau, range is how far the
observable itself travels, and collapse is called only under 5% of the range,
so an observable cannot pass by standing still.

Holding `J = 1` and sweeping gamma from 0.01 to 0.20, which sweeps `Q` from 100
to 5:

| Observable | spread | range | spread/range | collapses? |
|---|---|---|---|---|
| `S(rho_A)` | 0.789752 | 1.000000 | 0.7898 | no |
| `Tr(rho^2)` | 0.057319 | 0.491263 | 0.1167 | no |
| `CPsi` | 0.259139 | 0.307348 | 0.8431 | no |
| Concurrence | 0.861100 | 0.959735 | 0.8972 | no |

Holding `Q = 20` instead and sweeping the same gammas, with `J` moving along:

| Observable | spread | range | spread/range | collapses? |
|---|---|---|---|---|
| `S(rho_A)` | 0.000000 | 0.999986 | 0.0000 | yes |
| `Tr(rho^2)` | 0.000000 | 0.490770 | 0.0000 | yes |
| `CPsi` | 0.000000 | 0.307348 | 0.0000 | yes |
| Concurrence | 0.000000 | 0.959735 | 0.0000 | yes |

Both arms come out as the identity requires, and the producer raises if either
one does not. At fixed `Q` the curves land on each other to machine precision:
`tau = gamma*t` **is** this generator's own time. At fixed `J` they do not, and
the reason is not that tau is the wrong clock. Holding `J` while sweeping gamma
moves `Q`, so those five runs are five different systems compared at matched
tau, and what the failure measures is the second knob.

The left table alone, read as "irreversible observables do not scale with tau",
mistakes a change of system for a failure of the time variable. Read together,
the two tables say something narrower and firmer: a trajectory here needs two
numbers, and gamma supplies one of them. Fixing the clock only alongside a
second knob is not what "gamma is time" claims.

One observable had to be dropped from these tables rather than reported.
`arg(rho_01)` looks like a phase and is not: `rho_01` links popcount 0 to
popcount 1, and both the Hamiltonian and the Z-dephasing conserve popcount, so
for this preparation the element is identically zero at every gamma and every
`t`. Its neighbour `arg(rho_12)`, the `|01><10|` coherence this state does
populate, is no better as a table row: it takes three values, so its spread
equals its range by construction and it can only ever print a failure. Both are
observables whose verdict is fixed before the run.

## Precise language

| Statement | Status |
|---|---|
| Gamma sets a dissipative timescale in the declared generator | Established by the model definition |
| The tested positive-gamma trajectories damp | Finite numerical observation |
| The tested zero-gamma trajectories are stationary or recurrent | Finite numerical observation for two preparations |
| `tau = gamma*t` alone determines the trajectory at fixed J | False; measured to fail on all four observables above |
| `tau = gamma*t` determines the trajectory at fixed `Q = J/gamma` | True; measured to collapse exactly |
| Gamma defines experienced time or its origin | Not established |
| A nonzero dissipator identifies its microscopic source | Not established; see the [Incompleteness Proof](proofs/INCOMPLETENESS_PROOF.md) |

Hamiltonian dynamics still uses `t` at gamma = 0, a stationary preparation says
something about that preparation rather than about whether time exists, and
dissipative damping does not by itself define experience. The numbers on this
page are consistent with all three of those, and with nothing stronger.

## Reproduction

```bash
python simulations/gamma_is_time_proof.py
python simulations/two_qubits_no_noise.py
python simulations/gamma_unit_scaling_gate.py
```

Interpret the resulting numbers only within their stated system, preparation,
time window, and observable set.
