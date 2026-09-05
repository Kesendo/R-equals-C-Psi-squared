# N=3 Oscillation and Pauli-Fingerprint Analysis

<!-- Keywords: N=3 oscillation census, Pauli fingerprint, conditional standing
wave, Heisenberg dephasing, state Hamiltonian cross table -->

**Status:** Computational reading at `N=3`, 6 Hamiltonians and 8 initial states
**Date:** March 19, 2026
**Script:** [standing_wave_analysis.py](../simulations/standing_wave_analysis.py)

---

## Scope

This experiment decomposes specified `N=3` evolutions into Liouvillian modes
and reports which Pauli observables have oscillatory weight. It does not prove
that every palindromic eigenvalue pair is a pair of spatially
counter-propagating waves.

F1 supplies the linear spectral transport
`lambda -> -lambda-2 Sigma_gamma`, or `mu -> -mu` after centering. Complex
conjugation is a separate spectral closure. A standing-wave interpretation is
available only for diagonalizable or semisimple pairs on the centered
imaginary axis, after opposite spatial propagation and the relevant
excitation/readout overlaps have been independently established. When
`Re(mu) != 0`, relative envelope drift remains; defective blocks add Jordan
terms. The calculations below are therefore described directly as oscillation
and Pauli-fingerprint measurements.

## Setup and spectral reading

The producer uses `N=3`, `J=1`, local Z-dephasing `gamma=0.05`, and a full
64-dimensional Liouvillian eigendecomposition. It reports 32 palindromic
pairs with maximum matching error `1.44e-14`.

| Category | Pairs | Frequencies |
|---|---:|---|
| Steady-XOR | 4 | centered `mu=+/-0.15`, real |
| Decay | 8 | centered `mu=+/-0.05`, real |
| Fundamental | 8 | `omega approx 2J` |
| Second harmonic | 4 | `omega approx 4J` |
| Third harmonic | 8 | `omega approx 6J` |

In this run, the measured `|Re(mu)| approx 0.017` for the oscillatory bands is
roughly 100--350 times smaller than their frequencies. This finite numerical
ratio describes slow envelope drift in the tested model; it is not a proof of
a stationary spatial pattern.

## Pauli-observable results

For the tested Hamiltonian-by-state grid, `ZZZ` has zero reported oscillatory
weight. Oscillatory readings occur in Pauli strings including `IYY`, `XXZ`,
`ZXX`, `YYI`, `XZX`, and `YIY`. “Node” and “antinode” may be used as shorthand
for zero and nonzero oscillatory observable weight in this finite table, not as
a universal spatial-wave theorem.

| | Heisenberg | XY | Ising | DM | XXZ | Heisenberg+DM |
|---|---:|---:|---:|---:|---:|---:|
| GHZ | 0% | 0% | 0% | 0% | 0% | 0% |
| W | 0% | 5.6% | 44.4% | 50% | 1.3% | 10.4% |
| Bell | 48.6% | 40.6% | 50% | 40.6% | 65.5% | 65.3% |
| `|+++>` | 0% | 40.6% | 62.5% | 40.6% | 38.3% | 43.4% |

The percentage is the producer's fraction of state weight in modes it labels
oscillatory. Thus, within this exact `N=3` grid:

- GHZ has 0% in all six columns.
- Bell has 40.6--65.5% in all six columns.
- W varies from 0% to 50%, demonstrating Hamiltonian dependence.
- `ZZZ` has zero oscillatory weight for all eight tested states and all six
  Hamiltonians.

These are finite-grid empirical statements. They do not license “for every
Hamiltonian,” “for every initial state,” or an extension to `N>=4`.

## State fingerprints under the Heisenberg run

| State | Oscillatory weight | Top Pauli oscillator | Nonzero Pauli entries |
|---|---:|---|---:|
| `|010>` | 44.4% | `ZIZ` | 26 |
| `|+-+>` | 44.5% | `IXI` | 26 |
| Bell(0,1) | 48.6% | `XXZ` | 26 |
| Bell + `|+++>` | 20.2% | `YYI` | 48 |
| Bell + W | 28.8% | `YYX` | 50 |
| GHZ | 0% | none | 0 |
| W | 0% | none | 0 |

The table establishes dependence on both preparation and Hamiltonian. A
physical standing wave would additionally require evidence that the relevant
mode pair propagates in opposite spatial directions and that the preparation
and observable excite both members with the required relation.

## Reproducibility

- [standing_wave_analysis.py](../simulations/standing_wave_analysis.py): full
  eigendecomposition, state decomposition, and Pauli fingerprints
- [standing_wave_analysis.txt](../simulations/results/standing_wave_analysis.txt):
  committed output
- [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md): F1 partner map
- [Conditional standing-wave account](../docs/STANDING_WAVE_THEORY.md)
