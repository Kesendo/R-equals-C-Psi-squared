# Pi as a Centered Spectral Mirror

<!-- Keywords: Pi palindromizer, centered spectral mirror, structural time
reversal caveat, populations coherences, conditional standing wave -->

**Status:** Algebraic identity proven; physical time-reversal and standing-wave
readings require additional hypotheses
**Date:** March 19, 2026
**Verification:** [`simulations/pi_time_reversal_verify.py`](../simulations/pi_time_reversal_verify.py)

---

## Exact statement

For the F1 model family, the linear palindromizer satisfies

```text
Pi L Pi^-1 = -L - 2 Sigma_gamma I.
```

Writing `L_c = L + Sigma_gamma I` gives

```text
Pi L_c Pi^-1 = -L_c,
Pi exp(L_c t) Pi^-1 = exp(-L_c t).
```

Consequently, if `L v = lambda v`, then

```text
L(Pi v) = (-lambda-2 Sigma_gamma)(Pi v),
mu = lambda+Sigma_gamma  ->  -mu.
```

This is an exact structural mirror of the centered generator. It is not, by
itself, physical time reversal: it does not reverse the dissipative envelope,
construct a realizable reversed trajectory, supply an antiunitary operation,
or establish opposite spatial propagation.

Complex conjugation must be kept separate. Hermiticity preservation gives
conjugate closure of the spectrum; composing it with F1 yields the
frequency-preserving spectral representative
`-conj(lambda)-2 Sigma_gamma`. Pi acting on `v` alone yields the linear partner
`-lambda-2 Sigma_gamma`, whose imaginary part has the opposite sign.

## Pauli action

In one gauge for Z-dephasing, Pi acts locally as

```text
I -> X,    X -> I,    Y -> iZ,    Z -> iY.
```

It exchanges the all-site classes `{I,Z}` and `{X,Y}` and maps total XY
weight `k` to `N-k`. The first class commutes with the pure-dephasing
dissipator; the second is maximally exposed when every site carries `X` or
`Y`. Mixed Pauli strings need not be purely populations or purely
coherences, so those words should not be applied to the entire weight ladder
without qualification.

The exact rate statement is that a partner below the center is the **slow**
decay and its partner above the center is the **fast** decay:

```text
d_slow + d_fast = 2 Sigma_gamma.
```

## What the verification establishes

At `N=3`, the producer verifies 32/32 Pi transports into the appropriate
partner eigenspace, with maximum residual `2.68e-13`. The XY-weight sum is
`3` to maximum deviation `8.88e-16`, and `Pi(ZZZ)=-i YYY` in this gauge.
These are operator and eigenspace facts. Degeneracy means the invariant
partner subspace, rather than an arbitrarily chosen eigenvector, is the robust
object.

The separate `N=3` oscillation producer reports zero oscillatory weight for
`ZZZ` and nonzero weights for specified XX/YY-containing observables across a
finite 6-Hamiltonian by 8-state grid. Those readings remain in
[N=3 Oscillation and Pauli-Fingerprint Analysis](STANDING_WAVE_ANALYSIS.md);
they do not follow from `mu -> -mu` alone.

## Conditional standing-wave reading

For a diagonalizable or semisimple centered pair `mu=+/- i omega`, if its
eigenvectors have independently established opposite spatial propagation and
an excitation and observable address both with suitable amplitudes and phases,
even/odd combinations can form a physical standing pattern. None of those
spatial and operational conditions is supplied merely by the palindrome.

When `Re(mu) != 0`, the centered pair has relative growth/decay. For defective
blocks, Jordan chains produce factors `t^k exp(lambda t)`. In either case the
generic “forward plus backward equals standing wave” inference fails.

Past/future, decided/undecided, and present-as-interference remain philosophical
labels. They are not consequences of the F1 identity.

## Relation to time reversal

The phrase “time reversal” is retained only as a search label for this record.
The technically accurate name is **centered structural mirror**. The identity
relates two linear evolutions after removing a scalar envelope; it does not run
the physical Lindblad semigroup backward. See
[Time Irreversibility Exclusion](../docs/proofs/TIME_IRREVERSIBILITY_EXCLUSION.md)
and [On Two Times](../reflections/ON_TWO_TIMES.md).

## Reproducibility

- Script: [`simulations/pi_time_reversal_verify.py`](../simulations/pi_time_reversal_verify.py)
- Output: [`simulations/results/pi_time_reversal_verify.txt`](../simulations/results/pi_time_reversal_verify.txt)
- Proof: [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)
- Conditional wave account: [Standing Wave Theory](../docs/STANDING_WAVE_THEORY.md)
