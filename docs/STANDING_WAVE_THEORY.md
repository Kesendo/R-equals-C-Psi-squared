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

Conditions 1 and 2 are properties of the generator; condition 3 is a property
of the preparation and the readout. So a standing wave here is a state and
Hamiltonian property jointly, and not a property of either alone: the same
Hamiltonian rings for one preparation and sits still for another, and the same
preparation rings under one Hamiltonian and not under the next.

Only under these conditions can even/odd combinations reproduce the stationary
spatial pattern of two counter-propagating waves. If `Re(mu) != 0`, relative
growth and decay remain in the centered frame. If the block is defective,
Jordan chains add polynomial factors `t^k exp(lambda t)`. Neither case is
captured by the two-sinusoid guitar-string formula.

The fixed locus depends on the map. Linear F1 self-pairing requires the exact
complex point `lambda=-Sigma_gamma`; the whole line
`Re(lambda)=-Sigma_gamma` is fixed only by the composite
`lambda -> -conj(lambda)-2 Sigma_gamma`. Neither condition makes a mode a
spatial node or standing wave. In a defective block, F1 transports the whole
generalized eigenspace and its Jordan-chain data; a multiplicity census alone
does not choose a basis of modes.

For the committed `N=2,...,7` spectra, a numerical multiplicity-aware census
has 10,903
unordered two-member orbits and 34 fixed eigenvalues. The composite census has
9,921 unordered two-member orbits and 1,998 eigenvalues on its fixed line. Both
counts consume occurrences with multiplicity and are stable across matching
tolerances `1e-6`, `1e-8`, and `1e-10`; they are not exact-arithmetic counts,
and neither is a wave count.

## What the N=3 computation now measures

The repaired finite calculation in
[N=3 Direct Pauli-Observable Time Traces](../experiments/STANDING_WAVE_ANALYSIS.md)
propagates four specified states in the `N=3` Heisenberg chain and reports
sampled ranges of seven direct Pauli expectations. It deliberately does not
report `sum |c_k|^2` as eigenmode “state weight”: for a non-normal generator
that quantity is not invariant, and degeneracy makes it basis-dependent.

The direct table is an observable-time-trace result only. Its nonzero sampled
half-ranges show variation of the named expectation values; they do not by
themselves establish a spatial node, antinode, or standing wave.

The same producer counts the spectrum, which the eigenvector defect never
reached, and the count settles condition 1 at this size. Centred on
`mu = lambda + Sigma_gamma`, the 64 eigenvalues pair 64/64 into 32 F1 pairs at a
worst residual of `2.8e-14`, and **none of them is purely imaginary**: 24 are
purely real and 40 carry both decay and oscillation. The nearest mode to the
imaginary axis sits `0.016658` from it, about `10^10` times the eigensolver's
noise floor, so this is measured rather than thresholded.

Shared envelopes are not what is missing. The 40 mixed modes form 20 conjugate
pairs whose members have a common real part, drawn from just six values
(`±0.05`, `±0.0167243`, `±0.0166584`), so each pair shares an envelope exactly.
Condition 1 asks for more than a shared envelope: it asks for a **flat** one,
`Re mu = 0`, so that neither member grows against the other. That is what is
absent at `N=3`, and conditions 2 and 3 are never reached.

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
