<!-- QUARTER-CURRENT -->
# The Finite Envelope-Rise Atlas

Current reading: this page reports a finite N/Q/K catalogue using a fixed
rise-reporting bar and named controls.  It is neither an all-generator theorem
nor evidence that no unscanned rise exists.

Some old pages drew a border where the data had only placed a few lamps. This page keeps the lamps: a
finite atlas of named full-state CΨ peak readings, and an invitation to find the country around them.

## What is being read

`EnvelopeTheoremWitness.GlobalReading(N, J, γ, tMax, points)` evolves the named Bell+ carrier with the
dense `Symphony` engine. `QuarterEnvelope.Of` records local quadratic-apex estimates on the actual time
grid and compares each maximum with its predecessor; direction-tagged crossings are extracted by Symphony from
one paired time/direction event list; they are not owned by the envelope reader. A rise is reported only
when it exceeds `EnvelopeTheoremWitness.RiseReportingBar = 1e-3`; raw predecessor ordering and maximum
positive delta remain independent of that reporting convention. Quadratic interpolation is exact for a
quadratic local profile; its generic local error is `O(h³)`, so refinement remains part of the reading.
For a sampled flat maximum, the reader selects the plateau end and applies the same three-point estimator;
that fitted apex is an explicit convention, not an intrinsic continuous-time feature.

These are finite-window, finite-grid classifications. The finite atlas neither confirms nor refutes a
universal peak-envelope theorem. In particular, a zero row says “no rise was resolved here,” not “no rise
exists.” The older autonomous N=2 successive-local-maxima claim is unproved; the counterexamples in
`PROOF_MONOTONICITY_CPSI.md` settle pointwise monotonicity and absorption shortcuts, not that narrower
peak-sequence question.

Gate-first reader: `compute/RCPsiSquared.Diagnostics.Tests/Foundation/EnvelopeBoundaryTests.cs`.

## Retained rows

All rows below use γ=0.01, K_max=γt_max=0.25, and 1600 time-grid points unless stated otherwise.

| N | Q=J/γ | finite reading |
|---|------:|----------------|
| 3 | 2000  | no N=3 rise resolved; RiseCount=0 on this named grid/window |
| 4 | 13    | no rise resolved above the reporting bar |
| 4 | 40    | a rise is resolved above the reporting bar |
| 4 | 500   | 36 predecessor rises; maxΔ≈0.041 |
| 5 | 40    | no rise resolved above the reporting bar |
| 5 | 500   | 32 predecessor rises; maxΔ≈0.020 |

A finite scan placed the N=4 transition-like region near Q≈27 and the N=5 one near Q≈45. Those numbers
are atlas estimates, not critical constants: near the N=4 changeover, roughly Q≈18…28, a peak edging above
its predecessor is sensitive to phase sampling. The bracketing rows Q=13/40 and Q=40/500 are the durable
statements.

## Same-(N,Q,K) rescaling check

The independently evolved named N=4 pair

- `(J,γ,t_max)=(0.40,0.01,25)`, and
- `(J,γ,t_max)=(0.80,0.02,12.5)`

shares `(N,Q,K_max)=(4,40,0.25)`. The raw global and carrier-pair CΨ arrays agree pointwise, the rise count
matches, and the paired apex-height/dose comparison agrees to six decimals. Changing only the second J to 0.26 (Q=13)
fails the raw-curve comparison and changes the above-bar classification. This is a finite paired control
against an absolute-time leak in this reader. It is not an every-Q identity deduced from samples.

## What the atlas asks next

The interesting border is still undrawn. Does every N admit a rise somewhere? Is the apparent N=3 silence
real or merely outside these windows? Is there a useful Q-scale, and what physical mechanism selects it?
Those all-Q, all-N, and mechanism questions remain open under `envelope_n4_rise`. The rows above make the
questions sharper without pretending to answer them.

## Links

- Typed historical home: `CpsiEnvelopeTheoremClaim` (Tier1Derived; stable name, corrected content).
- Proof and scope: `docs/proofs/PROOF_MONOTONICITY_CPSI.md`.
- F-registry entry: F17 in `docs/ANALYTICAL_FORMULAS.md`.
- Live reader: `inspect --root envelope --N {3,4,5}`.
