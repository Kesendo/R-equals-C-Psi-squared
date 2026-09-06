# N=3 Direct Pauli-Observable Time Traces

<!-- Keywords: N=3 observable trace, Pauli expectation, Heisenberg dephasing,
conditional standing wave -->

**Status:** finite direct-propagation reading; no mode-weight or spatial-wave verdict
**Script:** [standing_wave_analysis.py](../simulations/standing_wave_analysis.py)

## What was repaired

The earlier report treated `sum |c_k|^2` in a right-eigenvector expansion as a
fraction of state weight. A non-normal Liouvillian has no orthonormal right
eigenbasis, and within degenerate eigenspaces the coefficients depend on the
chosen basis. Those percentages and every conclusion built from them are
withdrawn.

The replacement computes density-matrix propagation directly and reports
expectation values `Tr(P rho(t))` for named Pauli observables. These readings
do not depend on a choice of Liouvillian eigenvectors.

## Finite setup

- `N=3` open Heisenberg chain, `J=1`
- local Z-dephasing `gamma=0.05`
- four initial states: GHZ, W, Bell(0,1), and `|+++>`
- seven observables: `ZZZ`, `IYY`, `XXZ`, `ZXX`, `YYI`, `XZX`, `YIY`
- 101 samples on `t in [0,10]`

The committed output gives, for every state-observable row, the initial value,
sampled minimum, sampled maximum, and half-range. The half-range is only a
coordinate-free description of the sampled expectation trace. It is not a
spectral weight and not a standing-wave certificate.

## Scope fence

This run does not test spatial counter-propagation, balanced excitation of an
F1 partner pair, or the relative phase needed for a stationary spatial
pattern. It therefore supports statements about these direct observable
traces only. The exact F1 partner map remains a separate algebraic result.

## Reproduction

- [producer](../simulations/standing_wave_analysis.py)
- [committed output](../simulations/results/standing_wave_analysis.txt)
- [conditional standing-wave account](../docs/STANDING_WAVE_THEORY.md)
