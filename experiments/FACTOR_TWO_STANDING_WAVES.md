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
an eigenvalue-pair count does not remove them. The linear F1 fixed locus is the
single complex point `lambda=-Sigma_gamma`, not the whole centered vertical
line. F1 maps generalized eigenspaces (and their Jordan-chain lengths) at
`lambda` to those at `-lambda-2 Sigma_gamma`; it does not select individual
vectors inside a degenerate eigenspace.

## Finite pair census

The producer consumes the committed
`simulations/results/rmt_eigenvalues_N{2..7}.csv` spectra: 21,840 eigenvalues
for `N=2,...,7`. Each CSV contains only the columns `Re` and `Im`; it does not
embed a backend, source revision, invocation, or timestamp. The matching
regeneration route in the current C# source is
`dotnet run -c Release --project compute/RCPsiSquared.Compute -- rmt chain`,
whose source constants are chain topology, `J=1.0`, and uniform `gamma=0.05`.
That current-source route is reproducibility metadata, not metadata recovered
from inside the historical CSV artifacts.

At primary matching tolerance `1e-8`, repeating the complete census at
`1e-6`, `1e-8`, and `1e-10` gives the same multiplicities; these are
tolerance-stable eigensolver counts, not exact-arithmetic multiplicities. There
are two different complete orbit censuses, because there are two different involutions.
For the **linear F1 map** `lambda -> -lambda-2 Sigma_gamma`:

| N | Total | Unordered two-member F1 orbits | F1-fixed multiplicity at `lambda=-Sigma_gamma` | Accounted for |
|---|---:|---:|---:|---:|
| 2 | 16 | 6 | 4 | 100% |
| 3 | 64 | 32 | 0 | 100% |
| 4 | 256 | 121 | 14 | 100% |
| 5 | 1,024 | 512 | 0 | 100% |
| 6 | 4,096 | 2,040 | 16 | 100% |
| 7 | 16,384 | 8,192 | 0 | 100% |

Thus the linear census contains 10,903 unordered two-member orbits and 34
fixed eigenvalues, counted with algebraic multiplicity.

For the **composite map**
`lambda -> -conj(lambda)-2 Sigma_gamma`, obtained by following conjugate closure
with F1, the fixed locus is the full line `Re(lambda)=-Sigma_gamma`:

| N | Total | Unordered two-member composite orbits | Composite-fixed multiplicity on `Re(lambda)=-Sigma_gamma` | Accounted for |
|---|---:|---:|---:|---:|
| 2 | 16 | 3 | 10 | 100% |
| 3 | 64 | 32 | 0 | 100% |
| 4 | 256 | 52 | 152 | 100% |
| 5 | 1,024 | 512 | 0 | 100% |
| 6 | 4,096 | 1,130 | 1,836 | 100% |
| 7 | 16,384 | 8,192 | 0 | 100% |

This second census contains 9,921 unordered two-member orbits and 1,998 fixed
eigenvalues. Neither census counts physical traveling or standing waves. The
composite census preserves `Im(lambda)` and is useful for grouping equal
frequencies, but only the linear census is the orbit structure transported by
`Pi`.

The producer also emits the measured mean decay over each complete spectrum;
at the printed precision it is `Sigma_gamma`:

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

The producer's topology negative control is restricted to chain, star, and
ring at `N=4`; all three spectra are fully accounted for:

| Linear F1 census at N=4 | Chain | Star | Ring |
|---|---:|---:|---:|
| F1-fixed multiplicity at `lambda=-Sigma_gamma` | 14 | 16 | 24 |
| Unordered two-member F1 orbits | 121 | 120 | 116 |
| Accounted for | 100% | 100% | 100% |

Topology changes the spectrum and exact-point multiplicity. This table tests
the linear F1 orbit accounting in those models; it does not test spatial
counter-propagation.

## Factor two

The factor two is not a universal lifetime ratio between two classes of modes.
Within an F1 spectrum that reaches both endpoints, it is the full decay range
`0...2 Sigma_gamma` divided by its center `Sigma_gamma`, and pairwise it is the
sum `d_slow+d_fast=2 Sigma_gamma`. Calling this a cavity round trip is an
optional optical reading, not an additional dynamical theorem.

## Reproduction

- Script: [`simulations/factor_two_standing_waves.py`](../simulations/factor_two_standing_waves.py)
- Output: [`simulations/results/factor_two_standing_waves.txt`](../simulations/results/factor_two_standing_waves.txt)
- Inputs: committed `simulations/results/rmt_eigenvalues_N{2..7}.csv` files
  (tab-separated `Re`, `Im` only; no embedded run metadata)
- Current regeneration command: `dotnet run -c Release --project
  compute/RCPsiSquared.Compute -- rmt chain` (`J=1.0`, uniform `gamma=0.05`)
- Algebraic source: [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)
