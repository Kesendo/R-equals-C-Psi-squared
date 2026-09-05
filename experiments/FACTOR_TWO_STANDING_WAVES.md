# Palindromic Pair Census and the Conditional Standing-Wave Reading

<!-- Keywords: palindromic eigenvalue pair census, decay-rate sum, conditional
standing wave, centered spectrum, Pi linear transport, R=CPsi2 -->

**Status:** Pair census confirmed for the reported finite systems; standing-wave
language is conditional
**Date:** April 4, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Verification:** [`simulations/factor_two_standing_waves.py`](../simulations/factor_two_standing_waves.py)

---

## What is established

For the dephasing spin systems in scope, the F1 palindromizer obeys

```text
Pi L Pi^-1 = -L - 2 Sigma_gamma I.
```

Thus `L v = lambda v` implies

```text
L(Pi v) = (-lambda - 2 Sigma_gamma)(Pi v).
```

In centered coordinates `mu = lambda + Sigma_gamma`, the linear partner map is
`mu -> -mu`. It reverses both real and imaginary parts. Lindbladian spectra are
also closed under complex conjugation for a separate reason: Hermiticity
preservation. Composing that closure with F1 gives the frequency-preserving
spectral representative `-conj(lambda)-2 Sigma_gamma`; it is not the vector
produced by `Pi v` alone.

Both representatives have the same complementary decay-rate relation:

```text
d + d_partner = 2 Sigma_gamma,    d = -Re(lambda).
```

This rate sum is the theorem-level content. Calling its members forward and
backward waves, or calling their sum a standing wave, requires more.

## When a standing-wave reading is licensed

A physical standing wave needs independently established opposite spatial
propagation, compatible amplitudes/phases, and an excitation and observable
that see both members. For a diagonalizable pair on the centered imaginary
axis, `mu = +/- i omega`, whose eigenvectors have independently been shown to
propagate oppositely, suitable even/odd superpositions can then give the usual
stationary spatial pattern. The spectral palindrome alone supplies neither
spatial propagation nor that preparation/readout condition.

If `Re(mu) != 0`, the two centered factors also have relative growth and decay,
so they are not a stationary equal-envelope wave pair. If a block is defective,
its evolution contains Jordan-polynomial factors such as `t^k exp(lambda t)`;
an eigenvalue-pair count does not remove them. Self-pairing on the centered
line likewise means only that the eigenvalue lies on the fixed locus, not that
the mode is automatically a node of a spatial standing wave.

## Finite pair census

The producer enumerates 21,840 eigenvalues for `N=2,...,7` and finds the
following complete palindromic accounting:

| N | Total | Distinct partner members | Self-paired on `Re(lambda)=-Sigma_gamma` | Accounted for |
|---|---:|---:|---:|---:|
| 2 | 16 | 6 | 10 | 100% |
| 3 | 64 | 64 | 0 | 100% |
| 4 | 256 | 104 | 152 | 100% |
| 5 | 1,024 | 1,024 | 0 | 100% |
| 6 | 4,096 | 2,260 | 1,836 | 100% |
| 7 | 16,384 | 16,384 | 0 | 100% |

Equivalently, the census contains 9,921 distinct unordered pairs and 1,998
self-paired eigenvalues. These are multiplicity counts, not counts of physical
traveling or standing waves.

The measured mean decay over each complete spectrum is `Sigma_gamma`:

| N | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---:|---:|---:|---:|---:|---:|
| Mean decay | 0.100000 | 0.150000 | 0.200000 | 0.250000 | 0.300000 | 0.350000 |
| `Sigma_gamma` | 0.10 | 0.15 | 0.20 | 0.25 | 0.30 | 0.35 |

For a pair, the member closer to `Re(lambda)=0` is the **slow** decay and the
member closer to `Re(lambda)=-2 Sigma_gamma` is the **fast** decay. Therefore

```text
d_slow + d_fast = 2 Sigma_gamma.
```

The earlier inverse naming is not used here.

## Topology sweep

The same finite pairing census was run on chain, star, and ring topologies at
`N=3,4,5`; all spectra were fully accounted for. At `N=4`:

| N=4 | Chain | Star | Ring |
|---|---:|---:|---:|
| Self-paired eigenvalues | 152 | 126 | 156 |
| Distinct unordered pairs | 52 | 65 | 50 |
| Accounted for | 100% | 100% | 100% |

Topology changes the spectrum and fixed-locus multiplicity. This table tests
the F1 pairing in those models; it does not test spatial counter-propagation.

## Factor two

The factor two is not a universal lifetime ratio between two classes of modes.
Within an F1 spectrum that reaches both endpoints, it is the full decay range
`0...2 Sigma_gamma` divided by its center `Sigma_gamma`, and pairwise it is the
sum `d_slow+d_fast=2 Sigma_gamma`. Calling this a cavity round trip is an
optional optical reading, not an additional dynamical theorem.

## Reproduction

- Script: [`simulations/factor_two_standing_waves.py`](../simulations/factor_two_standing_waves.py)
- Output: [`simulations/results/factor_two_standing_waves.txt`](../simulations/results/factor_two_standing_waves.txt)
- Algebraic source: [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)
