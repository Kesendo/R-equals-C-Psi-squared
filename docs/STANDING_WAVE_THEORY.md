# Standing Waves: Conditions Beyond the Spectral Palindrome

<!-- Keywords: standing wave conditions, palindromic pair, opposite spatial
propagation, semisimple imaginary pair, Jordan caveat, Pi centered mirror -->

**Status:** F1 algebra proven; standing-wave interpretation conditional;
philosophical analogy separated
**Authors:** Thomas Wicht, Claude (Anthropic)

---

## The exact algebra

The F1 palindromizer for the dephasing spin family satisfies

```text
Pi L Pi^-1 = -L - 2 Sigma_gamma I.
```

It follows linearly that

```text
lambda -> -lambda-2 Sigma_gamma,
mu=lambda+Sigma_gamma -> -mu.
```

This establishes a centered spectral reflection and the complementary decay
sum `d+d_partner=2 Sigma_gamma`. It does not establish that either eigenmode is
a spatial traveling wave.

Hermiticity preservation separately closes the Lindbladian spectrum under
complex conjugation. Composing conjugation with F1 gives the
frequency-preserving spectral representative
`-conj(lambda)-2 Sigma_gamma`. That composite operation must not be confused
with the vector transport `v -> Pi v`, which yields `-lambda-2 Sigma_gamma`.

## What a physical standing wave requires

For two modes to support the familiar standing-wave construction, the analysis
must additionally establish all of the following:

1. The relevant centered pair is diagonalizable or semisimple and lies on the
   imaginary axis, `mu=+/- i omega`, so the two members share an envelope.
2. Their eigenvectors represent opposite **spatial** propagation under a
   specified geometry and observable. The sign reversal of a complex
   eigenvalue is not evidence of spatial direction.
3. The preparation excites both members, and the readout sees their coherent
   combination with the required amplitude and phase relation.

Only under these conditions can even/odd combinations reproduce the stationary
spatial pattern of two counter-propagating waves. If `Re(mu) != 0`, relative
growth and decay remain in the centered frame. If the block is defective,
Jordan chains add polynomial factors `t^k exp(lambda t)`. Neither case is
captured by the two-sinusoid guitar-string formula.

A self-paired eigenvalue on `Re(lambda)=-Sigma_gamma` is merely on the spectral
fixed locus. It is not automatically a spatial node or standing wave.

## What the N=3 computation measured

The finite calculation in
[N=3 Oscillation and Pauli-Fingerprint Analysis](../experiments/STANDING_WAVE_ANALYSIS.md)
uses 6 Hamiltonians and 8 initial states. In that grid:

- `ZZZ` has zero reported oscillatory weight.
- Bell preparations have 40.6--65.5% oscillatory weight across the six tested
  Hamiltonians.
- GHZ has 0% across those same six columns.
- W varies from 0% to 50%, showing that the reading depends on both state and
  Hamiltonian.

The Heisenberg run includes oscillatory Pauli observables at frequencies near
`2J`, `4J`, and `6J`. These are real finite-model observations. “Node” and
“antinode” are acceptable shorthand inside that measured observable table, but
not universal labels derived from F1.

## What Pi means here

Pi exchanges local Pauli classes `{I,Z}` and `{X,Y}` and maps total XY weight
`k` to `N-k`. This explains the complementary dephasing exposure. It does not
turn classical populations into a physical past or quantum coherences into a
physical future. Those temporal words are an interpretive analogy.

Likewise, the centered propagator identity

```text
Pi exp((L+Sigma_gamma I)t) Pi^-1
    = exp(-(L+Sigma_gamma I)t)
```

is a structural mirror after a scalar envelope is removed, not a general
physical time-reversal operation and not a reversed Lindblad trajectory.

## Interpretive layer

The earlier picture of future and past meeting at a mirror can still be read as
philosophy: a standing wave as “present,” fixed observable components as
“memory,” and oscillatory components as “possibility.” Nothing in F1 proves
those identifications, a consciousness mechanism, retrocausality, or a stable
reality created by observation.

The analogy becomes a physics claim only in a concrete model that passes the
three standing-wave conditions above. Otherwise the supported statement is the
spectral palindrome and its rate sum.

## Sources

- [Mirror Symmetry Proof](proofs/MIRROR_SYMMETRY_PROOF.md): exact F1 identity
- [Pi as a Centered Spectral Mirror](../experiments/PI_AS_TIME_REVERSAL.md):
  operator action and `N=3` eigenspace verification
- [N=3 Oscillation and Pauli-Fingerprint Analysis](../experiments/STANDING_WAVE_ANALYSIS.md):
  finite empirical observable table
- [Palindromic Pair Census](../experiments/FACTOR_TWO_STANDING_WAVES.md):
  finite multiplicity census and rate sums
