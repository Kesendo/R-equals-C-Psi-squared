<!-- F14-CURRENT -->

# K-Dosimetry: finite purity-threshold readings and an unresolved discrepancy

**Status:** Unregenerated legacy calculation; the intermediate-rate discrepancy is unresolved.
**Date:** April 4, 2026.
**Source:** [k_dosimetry.py](../simulations/k_dosimetry.py),
[stored output](../simulations/results/k_dosimetry.txt).

The stored J = 1 calculation reports a 61.816% intermediate-rate deviation in
K = γt for |++⟩ at target purity 0.26. This is a historical numerical
discrepancy, not a certified physical reciprocity failure. With isotropic
two-qubit coupling and equal local Z-dephasing the stated swap-symmetric
preparation remains Hamiltonian-dead throughout its trajectory.

That distinction matters: F14 concerns a fixed readout and a fixed
dimensionless evolution book. A general fixed-J state can depend on J/γ,
but this particular symmetric preparation does not supply that counterexample.
The exact product-trajectory explanation is given in
[Observer-Gravity Bridge](OBSERVER_GRAVITY_BRIDGE.md). The legacy source
pairs independently computed eigensystems by column without proving their
correspondence. This repair leaves that reconstruction unchanged and does not
assign the discrepancy to a diagnosed mechanism.

The target here is **purity 0.26**, not CΨ = ¼. A purity-target time
must not be relabelled as a scalar quarter crossing or a measurement event.

<!-- F14-HISTORICAL -->

**Historical record:** The following tables retain the finite calculations.
They are not regenerated in this label repair; the source writes its result
at module scope and includes larger computations.

## Preparations at N = 2, γ = 0.05, purity target 0.26

| State | Initial purity | t_cross | K = γ × t |
|---|---|---|---|
| \|++⟩ | 1.000 | 19.7 | 0.983 |
| Bell+ | 1.000 | (not reached) | -- |
| \|01⟩ | 1.000 | (not reached) | -- |

The recorded Bell+ and |01⟩ arms did not reach the target. This is a
finite record, not a general classification of preparations by absorbed dose.

## The reported gamma sweep

| γ | t_cross | K = γ × t | Deviation |
|---|---------|-----------|-----------|
| 0.001 | 980.7 | 0.9807 | 0.00% |
| 0.01 | 101.2 | 1.0117 | 3.2% |
| 0.05 | 31.7 | **1.5870** | **61.8%** |
| 0.10 | 9.80 | 0.9805 | 0.03% |
| 0.50 | 1.96 | 0.9807 | 0.01% |
| 1.00 | 0.98 | 0.9805 | 0.03% |
| 5.00 | 0.20 | 0.9805 | 0.02% |

The displayed intermediate deviation is 61.8%, or about 62%; the stored
six-digit output reports 61.816%. Its overall K spread is 58.326%.
The sampled extreme-gamma rows cluster near K ≈ 0.98, but neither that
agreement nor the intermediate discrepancy certifies the legacy propagator.
The proposed explanation by competing Hamiltonian oscillations does not
apply to this swap-symmetric trajectory.

## Purity-target sweep

| Target purity | t_cross | K = γ × t |
|---|---|---|
| 0.900 | 0.54 | 0.027 |
| 0.750 | 1.56 | 0.078 |
| 0.500 | 4.41 | 0.220 |
| 0.400 | 6.64 | 0.332 |
| 0.300 | 13.5 | 0.677 |
| 0.260 | 19.9 | 0.994 |
| 0.255 | 25.3 | 1.265 |
| 0.251 | 32.2 | 1.612 |

These are the reported threshold readings, with distinct targets. Purity 0.25
is the maximally mixed value for two qubits; it is not interchangeable with
the algebraic CΨ quarter boundary.

## Uniform-rate bookkeeping across N

| N | K_qubit = γ × t | K_system = Σγ × t | Ratio |
|---|---|---|---|
| 2 | 0.983 | 1.966 | 2.0 |
| 3 | 0.894 | 2.681 | 3.0 |
| 4 | 0.913 | 3.653 | 4.0 |
| 5 | 0.670 | 3.349 | 5.0 |

K_system = (Σγ)t = Nγt = N·K_qubit is an arithmetic identity for a
uniform rate profile. The ratio N cannot demonstrate independent absorption.
The stated target also changes with N, as 1/d + 0.01.

## The sacrifice profile

| Profile | Σγ | t_cross | K_total = Σγ × t |
|---|---|---|---|
| Uniform | 0.200 | 18.3 | 3.653 |
| Sacrifice | 0.200 | 22.6 | 4.517 |

The stored sacrifice/uniform K_total ratio is 1.2367, about 24% higher.
These rows change total exposure to the selected threshold; they do not
demonstrate preserved dose, independent shares, or a cause for an MI benefit.
The [Beer-Lambert comparison](BEER_LAMBERT_BREAKDOWN.md) is a separate
coupled-system experiment.

<!-- F14-INTERPRETIVE -->

**Interpretive invitation, not a result:** A photographer can trade aperture
against shutter time. That makes dose an attractive picture for γt: bright
and brief, dim and long, a scene slowly appearing on film. The tables invited
us to ask whether a quantum state has a comparable sensitivity or response
curve. The picture remains worth exploring, but no photon absorption model
or physical crystallization event is defined by these purity thresholds.

<!-- F14-CURRENT -->

## Reproduction boundary

Read the [stored output](../simulations/results/k_dosimetry.txt) together with
its current-reading wrapper. The script is syntax-checked only here; its
legacy numerical reconstruction remains an open diagnostic task.
The independently specified fixed-J negative control lives in
[Gamma-Time Distinction](../docs/GAMMA_TIME_DISTINCTION.md).
