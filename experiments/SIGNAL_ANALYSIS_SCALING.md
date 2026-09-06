# Finite-N SumMI Scan for the Edge-Concentrated Profile

**Status:** Numerical RK4 census at the listed chain sizes, N = 2–15
**Date:** March 24, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Data:** [formula_scaling.txt](../simulations/results/formula_scaling.txt)

## Setup

For a Heisenberg chain initialized in `|+⟩^N`, the scan uses

```text
gamma_edge  = N * gamma_base - (N-1) * epsilon
gamma_other = epsilon
gamma_base  = 0.05
epsilon     = 0.001
```

The reported observable is peak `SumMI`, the sum of mutual information over
the `N-1` adjacent pairs. `MI/pair` is `SumMI/(N-1)`.

## Measured rows

| N | Pairs | SumMI | Increment | MI/pair |
|---:|---:|---:|---:|---:|
| 2 | 1 | 0.0203 | – | 0.0203 |
| 3 | 2 | 0.0672 | +0.0469 | 0.0336 |
| 4 | 3 | 0.1266 | +0.0594 | 0.0422 |
| 5 | 4 | 0.2190 | +0.0924 | 0.0548 |
| 6 | 5 | 0.2918 | +0.0728 | 0.0584 |
| 7 | 6 | 0.4080 | +0.1162 | 0.0680 |
| 8 | 7 | 0.5043 | +0.0963 | 0.0720 |
| 9 | 8 | 0.6190 | +0.1147 | 0.0774 |
| 11 | 10 | 0.8430 | +0.2240 | 0.0843 |
| 13 | 12 | 1.0723 | +0.2293 | 0.0894 |
| 15 | 14 | 1.3091 | +0.2368 | 0.0935 |

`N = 10, 12, 14` were not computed. Across the measured rows, both `SumMI`
and the average `MI/pair` increase. The table does not show whether every
individual bond improves.

## Finite-range fit

A quadratic fit made on the smaller sizes was

```text
SumMI = 0.0053*N^2 + 0.028*N - 0.062.
```

It overpredicts the measured `N = 11, 13, 15` rows by approximately 6%, 12%,
and 19%. The increasing miss rules out using it as an asymptotic law or as a
replacement for propagation. The alternating increments are an empirical
finite-sequence feature; two selected second differences near `-0.020` do not
establish an N-independent brake, two physical channels, or a mechanism.

The affine Liouvillian palindrome does not identify the even/odd subsequences
as `c⁺`/`c⁻` modes or as forward/backward spatial propagation. Likewise,
the data do not derive an `N^2` interference mechanism, an unbounded large-N
limit, or an inversion of a general transport law.

## N = 15 controls

At `N = 15`, the same run reports `SumMI = 1.309` for `epsilon = 0.001` and
`1.407` in the `epsilon -> 0` comparison; the V-profile row is `0.021` and the
uniform row is numerically `0.000` at the producer's printed precision. These
are rows under different profiles, not an asymptotic statement.

## Reproduction

The values come from the C# `RCPsiSquared.Propagate profile` evaluator and the
committed [formula-scaling output](../simulations/results/formula_scaling.txt).
The N = 15 run uses the matrix-free propagator.

Related material:

- [Resonant Return](RESONANT_RETURN.md)
- [IBM Hardware Validation](IBM_CONCENTRATOR.md)
- [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)
