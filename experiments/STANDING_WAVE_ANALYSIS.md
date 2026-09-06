# N=3 Direct Pauli-Observable Time Traces and Centred Spectral Census

<!-- Keywords: N=3 observable trace, Pauli expectation, Heisenberg dephasing,
conditional standing wave, centred spectral census -->

**Status:** finite direct-propagation reading plus an eigenvalue census; no mode-weight or spatial-wave verdict
**Script:** [standing_wave_analysis.py](../simulations/standing_wave_analysis.py)

## Why the reading is coordinate-free

A non-normal Liouvillian has no orthonormal right eigenbasis, and within
degenerate eigenspaces the expansion coefficients depend on the chosen basis.
So `sum |c_k|^2` over right eigenvectors is not a fraction of state weight and
carries no reading at all.

Two things do survive that, and this page reports both. Density-matrix
propagation gives expectation values `Tr(P rho(t))` directly, with no
eigenvectors anywhere. And the eigenvalues themselves are invariant under any
change of eigenvector basis and under degeneracy, so the spectrum can be
counted even where its eigenvectors cannot be weighed.

## The centred spectral census

Centring on `mu = lambda + Sigma_gamma` turns F1 into `mu -> -mu`, and at
`N=3`, `gamma=0.05` the 64 Liouvillian eigenvalues give:

| reading | value |
|---|---|
| F1 pairing | 64/64 matched, 32 pairs, 0 unmatched |
| worst pairing residual | `2.8e-14` |
| purely imaginary `mu` (standing-wave candidates) | **0** |
| purely real `mu` | 24 |
| mixed decay and oscillation | 40 |

The zero is the load-bearing entry, and it is a measurement rather than a
threshold effect: the nearest mode to the imaginary axis sits `0.016658` away,
which is `1.1e+10` times the eigensolver's noise floor `256*eps*||L||`. Nothing
here is close enough to the axis for the count to depend on where a line is
drawn.

Read against the three conditions in the
[conditional standing-wave account](../docs/STANDING_WAVE_THEORY.md), the first
one already fails at `N=3`, though not for want of shared envelopes: the 40
mixed modes form 20 conjugate pairs with a common real part each, six distinct
values in all (`±0.05`, `±0.0167243`, `±0.0166584`). What condition 1 requires
is a **flat** shared envelope, `Re mu = 0`, and no pair has one. Spatial
counter-propagation and the phase relation, the second and third conditions,
are not reached.

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
pattern. It supports statements about these direct observable traces and about
the eigenvalue counts above, and about nothing further. The census counts modes;
it does not choose a basis of modes, and in a defective block it would not:
F1 transports a whole generalized eigenspace with its Jordan-chain data, and a
multiplicity count says nothing about that. The exact F1 partner map remains a
separate algebraic result.

## Reproduction

- [producer](../simulations/standing_wave_analysis.py)
- [committed output](../simulations/results/standing_wave_analysis.txt)
- [conditional standing-wave account](../docs/STANDING_WAVE_THEORY.md)
