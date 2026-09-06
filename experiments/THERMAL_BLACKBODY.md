# Thermal-Channel Numerical Census

<!-- Keywords: thermal Lindblad channel, amplitude damping, oscillating
eigenvalue count, n_bar, numerical tolerance -->

**Status:** finite `N=4` eigensolver census; no blackbody, EP, or wave mechanism
**Verification:** [`thermal_blackbody.py`](../simulations/thermal_blackbody.py),
[`thermal_ep_analysis.py`](../simulations/thermal_ep_analysis.py)

## Channel and endpoint

The tested generator adds local emission and absorption to the Heisenberg
chain with Z-dephasing:

```text
sigma- = |0><1|  (emission),    rate gamma_T (n_bar+1)
sigma+ = |1><0|  (absorption),  rate gamma_T n_bar.
```

The spontaneous-emission term is therefore present at `n_bar=0`. The zero row
is not the pure-Z-dephasing model and must not be compared to it as though only
thermal absorption had been switched off. The producer checks this explicitly:
`||L(0)||_F > 0` for a one-site emission-only control, while
`||L(1e-9)-L(0)||_F = 1.22e-10` verifies continuity near the endpoint.

## Finite census

For `N=4`, `J=1`, `gamma_Z=gamma_T=0.05`, the producer calls an eigenvalue
oscillating when `|Im lambda| > 1e-6`:

| `n_bar` | Oscillating count | Fraction | Mean decay | `Q_max` |
|---:|---:|---:|---:|---:|
| 0 | 210 | 82.031% | 0.300000 | 54.627 |
| 0.001 | 210 | 82.031% | 0.300200 | 54.544 |
| 0.01 | 210 | 82.031% | 0.302000 | 53.803 |
| 0.1 | 210 | 82.031% | 0.320000 | 47.397 |
| 0.5 | 210 | 82.031% | 0.400000 | 31.232 |
| 1 | 212 | 82.812% | 0.500000 | 22.167 |
| 2 | 210 | 82.031% | 0.700000 | 14.555 |
| 5 | 212 | 82.812% | 1.300000 | 7.558 |
| 10 | 212 | 82.812% | 2.300000 | 4.235 |
| 50 | 206 | 80.469% | 10.300000 | 0.986 |

These are numerical counts at the stated tolerance. They show that the count
changes on the tested grid; they do not establish an invariant fraction.

## Transition brackets, not exceptional points

The companion script samples `n_bar` and reports intervals across which the
tolerance-defined count changes. It calls them transition brackets. A collision
or a real-to-complex change is not enough to certify an exceptional point:
defectiveness, geometric-multiplicity loss, or Jordan-chain growth must be
tested independently. No such gate is present here, so this document makes no
EP claim.

## The blackbody reading, asked directly

The scope argument below says these data are the wrong kind of object to test a
radiation law. The producer also asks the question anyway, so the reading fails
on its own terms and not only on a definition.

**Planck fit to the frequency distribution.** Binned mode frequencies at
`nbar = 0.1, 0.5, 1, 5, 10`, fitted with the 1D Planck density
`omega/(exp(omega/T)-1)`. The best any fit reaches is `R^2 = 0.0480`, against
the `R^2 >= 0.90` a curve must clear to be called a fit; the five values are
0.0382, 0.0382, 0.0466, 0.0480, 0.0270. A bare exponential does equally badly
(0.0253 to 0.0446), and the widest gap between the two is 0.0034, far smaller
than the 0.8520 that separates the better curve from being a fit at all, so
ranking one above the other reads a preference out of noise.

**Absorption scaling.** The excess mean decay rate over the cold-bath baseline
is not quartic (Stefan-Boltzmann) or quadratic (a 1D Planck law). It is exactly
linear:

    excess = N * gamma_T * nbar

read off the trace rather than the spectrum, so the check needs no eigensolver
and no tolerance. At a dyadic rate `gamma_T = 0.25` the deviation is `0.0` to
the bit; at the census rate `0.05` it is `8.9e-16`, which is the sqrt-then-square
roundtrip in the jump operators rather than the law. The blackbody reading
predicts a power the generator does not carry.

## What is not inferred

The spectrum of this finite spin-chain Lindbladian is not a sampled radiation
intensity distribution. Consequently these data do not test Planck's law or
Stefan-Boltzmann scaling, do not measure a blackbody phase transition, and do
not show thermal photons activating or destroying physical standing waves.
Likewise, the decrease of the numerical `Q_max` is a spectral observation in
this model, not a mechanism for IBM hardware or a statement about biological
systems.

## Reproduction

- [census producer](../simulations/thermal_blackbody.py)
- [committed census](../simulations/results/thermal_blackbody.txt)
- [transition-bracket producer](../simulations/thermal_ep_analysis.py)
- [F52 registry entry](../docs/ANALYTICAL_FORMULAS.md#f52-thermal-oscillating-count-census-tier-2-numerical-n4)
