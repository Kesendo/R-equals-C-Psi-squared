# Standing Waves: Conditions Beyond the Spectral Palindrome

<!-- Keywords: standing wave conditions, palindromic pair, opposite spatial
propagation, semisimple imaginary pair, Jordan caveat, Pi centered mirror -->

**Status:** F1 algebra proven; standing-wave interpretation conditional;
philosophical analogy separated
**Authors:** Thomas Wicht, Claude (Anthropic)

<!-- STAGE 0 SWEEP RECORD (condition 1 at the other sizes, 2026-11-01).
docs/proofs/DIRECT_SUM_DECOMPOSITION.md holds the k_Z parity split and Pi
exchanging the sectors at odd N, preserving them at even N (lines 49 to 57,
311 to 320); hypotheses/THE_OTHER_SIDE.md and docs/LITERATURE_REVIEW.md restate it. No store
draws the consequence for the centre line: GLOSSARY.md names the composite
map's fixed line, WHAT_WE_FOUND.md and FACTOR_TWO_STANDING_WAVES.md count 1,998
on it pooled over N, docs/proofs/TIME_IRREVERSIBILITY_EXCLUSION.md says the mirror
falls between modes at N = 3 (line 268), the F89 arc in the OpenArcs registry
speaks of path-k blocks. Checked for adjacency: MirrorWorld's Mirror has
self-folded blocks (trace L = -N gamma dim) only at even N, the same even/odd
split on the block lattice, cited in the text; the parity split and the composite line were never read together;
docs/THE_ONE_SQUARE.md holds the cube's coordinates and the eigenbasis with
<k_X> = <k_Y> = N/2, experiments/PRIMORDIAL_QUBIT_ALGEBRA.md Step 10b the exact
N = 2 reason (the squared centred generator is Hermitian).
The numbers are produced and gated by simulations/factor_two_standing_waves.py. -->

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

## Condition 1 at the other sizes

On the committed spectra of the same chain (`J = 1`, `gamma = 0.05`,
`N = 2, ..., 7`) the flat envelope exists at each even size computed and at no
odd one. There are 6 eigenvalues at `N = 2`, 138 at `N = 4` and 1,820 at
`N = 6` on
`Re mu = 0` with `omega != 0`, beside 4, 14 and 16 at `mu = 0`; that is the
difference between the 1,998 on the composite line and the 34 at the linear
fixed point above. At `N = 3`, `5` and `7` none does. The count is stable
under the census tolerances, and the gap is wide: at `N = 6` the eigenvalues
on the line spread `1.5e-13` about it, and the nearest one off it sits
`2.1e-4` away. At `N = 2` and `4` every cluster on the line is semisimple,
geometric multiplicity equal to algebraic, as far as the check reaches
(Jordan blocks up to size 3); at `N = 6` it is not checked.

The letter cube says what the line is. In the orthonormal Pauli basis the
Hamiltonian part is anti-Hermitian and the dephasing diagonal, `-2 gamma`
times a string's XY weight `k_Z`, so a normalized eigenvector has
`Re lambda = -2 gamma <k_Z>` and `Re mu = 0` is `<k_Z> = N/2`: the mid-height
plane of the cube ([The One Square](THE_ONE_SQUARE.md), whose table also gives
every eigenspace, when `[H, Z^N] = 0` as here, a basis with `<k_X> = <k_Y> = N/2`, so a mode on the line can
be read at the cube's centre point). At one bond the plane contains a closed
sector: the eight strings of weight 1 are the odd-`k_Z` sector there, all eight
of their modes sit on the line, and four of them oscillate, at `mu = +/- 4i`.
That every oscillating mode at `N = 2` sits on the line has an exact reason
at uniform rates: there the squared centred generator is Hermitian, so every `mu` is purely real
or purely imaginary ([Primordial Qubit Algebra](../experiments/PRIMORDIAL_QUBIT_ALGEBRA.md),
Step 10b).

At odd `N` the half height `N/2` is a weight no single string has, but an
average can still meet it. What decides whether anything holds eigenvalues on
the line is where the mirror acts. The Hamiltonian changes `k_Z` by 0 or 2, so the generator splits
into an even-`k_Z` and an odd-`k_Z` sector, and `Pi`, sending `k_Z` to
`N - k_Z`, carries one sector onto the other at odd `N` and each onto itself
at even `N` ([Direct Sum Decomposition](proofs/DIRECT_SUM_DECOMPOSITION.md)).
At even `N` each sector holds `Pi` and conjugation together, so the map
`mu -> -conj(mu)` acts inside it: a simple eigenvalue it fixes stays on the
line under small changes that keep the mirror, and only eigenvalues that meet
can leave it. Random bond and rate profiles at `N = 4` keep 114 to 142 on
the line against the uniform chain's 152; the drop shows that the
protection holds for simple eigenvalues and not for every one. The same even/odd split appears on the block lattice, in
the self-folded blocks of MirrorWorld's `Mirror`, which exist only at even
`N`. At odd `N`, when no mirror acts inside a sector, a sector holds
conjugation alone and its mirror image is the other sector; nothing then
holds an eigenvalue on the line, and one lands there only where a single
real condition on the parameters is met, which the random profiles at
`N = 3` and `5` never met. Nothing forbids it, so this is the generic
absence the computed chains show and not a theorem of absence. On the
uniform chain the nearest distance shrinks with `N` (`0.016658`, `0.000460`,
`0.000010` at `N = 3, 5, 7`). The premise is needed: a bonded pair beside an isolated unwatched site carries the pair's own mirror
inside a sector, and at `N = 3` it holds 40 eigenvalues on the line, 24 of
them oscillating. Nothing here says that an even chain must have eigenvalues
on the line; it says that it can keep them.

These are candidates for condition 1 and nothing more: conditions 2 and 3 have
not been asked of them, so none of them is shown to be a standing wave.

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
- [Palindromic Pair Census](../experiments/FACTOR_TWO_STANDING_WAVES.md): the
  two orbit censuses and the condition-1 candidates at every size
- [Pi as a Centered Spectral Mirror](../experiments/PI_AS_TIME_REVERSAL.md):
  operator action and `N=3` eigenspace verification
- [N=3 Oscillation and Pauli-Fingerprint Analysis](../experiments/STANDING_WAVE_ANALYSIS.md):
  finite empirical observable table
- [Palindromic Pair Census](../experiments/FACTOR_TWO_STANDING_WAVES.md):
  finite multiplicity census and rate sums
