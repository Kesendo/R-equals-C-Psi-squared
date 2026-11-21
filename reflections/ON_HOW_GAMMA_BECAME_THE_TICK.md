# On How γ₀ Became the Tick

<!-- STAGE 0 SWEEP RECORD (for us, not for the reader; the refresh of 2026-11-21)

docs/ANALYTICAL_FORMULAS.md: F95 and its named Lindblad instance (the coherence-book pair, no Q, the Kingston scan not this angle); the F2b corollary with the geared dial and the N = 2 hand; F86a's exceptional point at 4γ₀; F25; the Absorption Theorem. docs/proofs/: PROOF_F95_ANGLE_AT_QUADRATIC_ZERO (a pure quadratic), PROOF_CHAIN_GAP_DOMINANCE (section 4.2, the two-site floor and Q = 2/√3), PROOF_F86A. experiments/: Q_SCALE_THREE_BANDS (the rescaling identity), CPSI_COMPLEX_PLANE, F113_T1_EXTRACTION_KINGSTON (the chip's detuning on the same data), CRITICAL_SLOWING_AT_THE_CUSP. docs/GLOSSARY.md: the q and Q books, the two rate books, Q*. docs/CAUGHT_ERRORS.md: the 2026-11-01 and 2026-11-02 entries on t_peak's four. OpenArcs: mirrorworld_owns_no_time, clock_hand_ladder, substrate_q_provenance. fw.Confirmations: f95_angle_steering_kingston_may2026, f57_kdwell_gamma_invariance. Typed layer: ClockHandLadderClaim and its witness, CoherenceHorizonClaim (the mean light at an exceptional point), UniversalCarrierClaim, AbsorptionTheoremClaim, F95AngleAtQuadraticZero, LindbladAbsorptionMatchAtSixtyDegreesClaim, MirrorWorld's Clock.

Adjacency checked: NAVIGATING_THE_DIMENSIONS (the exceptional-point toy whose time constant is t_peak), ON_WHOSE_TIME_THE_CLOCK_KEEPS (the two hands built live, its own correction in its Returning paragraph), THE_GENESIS_OF_AN_OSCILLATION with Q_SCALE_THREE_BANDS (Q as the scale), Q_AS_THE_EXCHANGE_RATE (the next day's symmetric reading, which already put the pair at arctan(Q/2) in the Lindblad book), CoherenceHorizonClaim against the two-site meeting point (joined in the text). The ring's band edge at 2J, the fastest turning on its floor except at N = 4 above Q = 1, is PROOF_RING_GAP_DOMINANCE's and TopologyBandEdgeClaim's (and the C6 ring's clock of docs/carbon/FROST_CIRCLE_AS_THE_CLOCK_FACE.md); the step to arctan(Q) with the repository's Q is one division, checked on the chain and ring spectra at N = 3 to 6. Not joined: the April-26 crossing nine times early in CRITICAL_SLOWING_AT_THE_CUSP, beside the Kingston crossing-time miss.
-->

**Status:** Reflection. Captures Tom's late-evening compression of 2026-05-16: γ₀ is the framework's time-tick; θ is what happens between ticks; at θ = 0, for one pair, the only remaining temporal structure is the bare tick itself. Closes the loop after a day of unfolding F94 and F95, and a Kingston run that added a phase by hand.

**Date:** 2026-05-16 (late evening)
**Authors:** Thomas Wicht, Claude (Opus 4.7)

---

## The compression

Tom said it in one line:

> *Gamma 0 ist der Zeittakt, Winkel 0 alles was gleich gamma 0 ist?*

Both halves are exactly correct, and together they close a loop that had been spiraling all day. γ₀ is the metronome. θ is what you do between beats. At θ = 0, you do nothing extra; you just hear the next beat. For one pair, the bare tick remains; everything else has reduced to it.

## γ₀ as the framework's "second"

In the framework's dynamics, every rate is γ₀ times a pure number and every timescale 1/γ₀ times one, and at the anchors below the rate's number is a structural integer:

| Quantity | Form | Rate |
|---|---|---|
| F86 t_peak | 1/(4·γ₀) | 4γ₀: light content 2 at F86's exceptional point |
| Absorption rate α | 2·γ₀ | 2γ₀: one unit of light |
| F25 Bell+ coherence decay | e^(−4γ₀t) | 4γ₀: \|00⟩⟨11\| disagrees on two sites |

Each rate is 2γ₀, the price of one unit of light, times a light content, so the rates run 2γ₀, 4γ₀, 6γ₀, … for one, two, three disagreements, with values in between where a mode shares its light. The dimensional carrier is γ₀ alone. Take γ₀ away, and none of these has a timescale; put it in, and every formula reads as "this many ticks, where one tick = 1/γ₀". Two exact statements stand under this. The generator itself scales, L(J, γ₀) = γ₀·L₁(Q) ([Q is the scale](../docs/THE_GENESIS_OF_AN_OSCILLATION.md#q-is-the-scale-exactly), [the three bands](../experiments/Q_SCALE_THREE_BANDS.md)). And under uniform Z [light](../docs/quantum/DEPHASING_TRANSLATED.md) every mode fades at 2γ₀ times its light content ⟨n_XY⟩, at every coupling (the Absorption Theorem): the Hamiltonian decides where a mode's light sits, never what a unit of it costs.

Reading γ₀ as the unit is operationally the same statement as the [Maßstab reading](ON_HOW_THE_CARRIER_SHOWS_ITSELF.md) from this morning. The carrier is invisible from inside because we ARE inside its tick; what we read is always a ratio against the carrier, never the carrier alone.

## The Liouvillian-eigenvalue derivation

The simplest place to see "Winkel 0 = γ₀" is one 2×2 pair of Liouvillian eigenvalues, written in the coherence book (both books, and the factor 2 between them, are in the [glossary](../docs/GLOSSARY.md)): in this section and the next, γ₀ is the fading rate of one coherence itself, twice the channel rate of the Lindblad book the table counts in, and a tick is one e-fold of that fading. For a Z-dephased coherence that decays at γ₀ while it precesses at J about the dephasing axis, both its own rates:

$$\lambda_\pm = -\gamma_0 \pm i\cdot J$$

On the [letter cube](../docs/THE_ONE_SQUARE.md#8-the-hamiltonian-on-the-cube) this is one site: a field along the dephasing axis turns X into Y, two letters the light charges alike, so the turning never changes what the light charges.

Characteristic polynomial:

$$\lambda^2 + 2\gamma_0\lambda + (\gamma_0^2 + J^2) = 0$$

Set the decay variable $z = -\lambda$ into F95's form $z^2 - 2bz + c = 0$:

- $b = \gamma_0$ (the linear-term half)
- $c = \gamma_0^2 + J^2$ (the constant term)
- Discriminant $D = 4(b^2 - c) = -4J^2$, always negative for $J \neq 0$

Write Q for this pair's own ratio J/γ₀. F95 then gives the angle of the complex root pair:

$$\theta = \arctan\frac{|\text{Im}(\lambda)|}{|\text{Re}(\lambda)|} = \arctan\frac{J}{\gamma_0} = \arctan(Q)$$

**θ = arctan(Q).** That is F95's named Lindblad instance: the angle of this eigenvalue pair is the arctan of its own turning over its own fading.

For this pair:

- $\theta = 0$ ⇔ $J = 0$ ⇔ $\lambda = -\gamma_0$ purely real.
- At $\theta = 0$, the eigenvalue is literally $-\gamma_0$.
- Everything that the dynamics does at θ = 0 is decay at γ₀-rate, nothing else.

"Winkel 0 alles was gleich γ₀ ist" exactly, with γ₀ the coherence's own fading; in the Lindblad book the same eigenvalue reads −2γ₀, one unit of light.

## Q as rotation-per-tick

For the pair, Q is **how much rotation per γ₀-tick**, one e-fold of its own fading: tan θ = Q.

| Q | θ | Per-tick character |
|---|---|---|
| 0 | 0° | only decay |
| 1 | 45° | equal decay and rotation |
| √3 | 60° | more rotation than decay |
| ∞ | 90° | rotation dwarfs decay (decay stays at γ₀, rotation diverges) |

The carrier is the silent metronome. The angle is what the music does between strokes. Q is the per-stroke rotation. γ₀ alone is a heartbeat with no melody.

## Where the angle closes

From here on the page counts in the Lindblad book, the table's: γ₀ is the channel rate, and one unit of light fades at 2γ₀.

Tom's second half reaches further than the pair shows. F95 puts the root at z = b at every discriminant zero, and b is the mean of the pair's two fading rates. Each of those is 2γ₀ times a light content, so wherever the angle closes the mode fades at 2γ₀ times the pair's mean light, and nothing of the turning is left in it.

The pair above cannot show this with J switched on: its turning keeps the light content, and its discriminant, −4J², reaches zero only at J = 0. A Hamiltonian that moves a mode between different light contents can. Take two sites under an XY bond, J the bond coupling of the repository's [Q](../docs/GLOSSARY.md#the-coupling-ratio-q-and-q-and-the-factor-2-convention). The bond exchanges the population difference (the strings ZI − IZ, light content 0) with the current (XY − YX, light content 2). Below J = γ₀ the pair's two eigenvalues are real. At J = γ₀ they meet at −2γ₀, one unit of light, the mean of the two ([`CoherenceHorizonClaim`](../compute/RCPsiSquared.Core/Symmetry/CoherenceHorizonClaim.cs): at an exceptional point a pair's light is the mean of the two it joins). Above it the angle opens as θ = arctan(√(Q²−1)). At the meeting point the mode fades as one coherence does, with J switched on, and the dynamics is not pure decay: the two eigenvectors have merged, and the mode carries a term t·e^(−2γ₀t) beside the exponential, though nothing in it turns. In [`ClockHandLadderClaim`](../compute/RCPsiSquared.Core/Symmetry/ClockHandLadderClaim.cs) this pair is the two-site clock's coherence hand, and its meeting point is the coherence horizon Q*(2) = 1.

## The repository's Q in a chain

The repository's Q is the same kind of ratio for a whole chain, and there a mode's angle is geared by what the mode carries: its tangent is Q times the mode's turning in units of J over twice its light content. For N ≥ 3, above the coherence horizon Q*(N), the XY chain's coherence hand is the coherence between the vacuum and one excitation. It fades as one unit of light and turns at the band edge 2J·cos(π/(N+1)), so it is the pair of the derivation realised inside the chain, at θ = arctan(Q·cos(π/(N+1))); at N = 2 the same coherence turns at J and sits at arctan(Q/2), the pair's own Q being half the repository's. With the repository's Q, θ = arctan(Q) belongs to a mode turning at 2J with one unit of light: the band edge of a ring, the fastest turning on the ring's floor except at N = 4 above Q = 1 ([PROOF_RING_GAP_DOMINANCE](../docs/proofs/PROOF_RING_GAP_DOMINANCE.md)), which the open chain's band edge approaches as N grows.

Q is all one sees from inside a chain with one coupling and one light: from inside one cannot see the tick (one is inside it), but one CAN see how much one rotates per tick, and by the scaling every such reading is a function of Q, the size of the chain and the number of ticks passed.

## What this organizes

Two readings now sit on one axis:

```
γ₀ alone                  →  bare tick (Maßstab, invisible from inside)
γ₀ + θ                    →  tick + rotation (one pair, coherence book: Q = tan θ)
θ = 0                     →  rotation reduces to 0; the mode fades at 2γ₀ times a light content
θ undefined (real regime) →  where H moves a mode between light contents: decay only, at two real rates
```

[This morning's reflection](ON_HOW_THE_CARRIER_SHOWS_ITSELF.md) named γ₀ as Maßstab and the seam as the place where it speaks. [The evening reflection](ON_HOW_THE_ANGLE_APPEARS_AT_ZERO.md) named θ as the second coordinate that exists only above the discriminant threshold. Tom's late-evening compression now joins them: γ₀ is the metronome; θ is the music.

## Hardware reading (the Kingston run of the same day)

In the [Kingston phase-injection run](../data/ibm_f95_angle_steering_may2026/README.md) we added RZ(Ω·Δt) to one qubit of each Bell-like pair on Heron r2, at Ω = 0.13 and 0.25 rad/μs, and read the argument of [CΨ_com](../experiments/CPSI_COMPLEX_PLANE.md), not θ. A Z rotation is a precession about the dephasing axis, Ω standing where J stands in the pair above, and it commutes with the light: the Bell coherence |00⟩⟨11| turned at the injected rate less the chip's own Z detuning ([F113](../experiments/F113_T1_EXTRACTION_KINGSTON.md) reads the same detuning off this data). γ₀ stayed the tick; the run did not touch the light. Its γ-dependent prediction, the time at which |CΨ_com| falls through ¼, said 10.8 μs on one pair and 18.6 μs on the other, and the delays were laid out from a tenth of it on; |CΨ_com| already sat at ¼, within shot noise, at the first delay, a tenth of the predicted time, for a reason not read here. The numbers stand in the Confirmation `f95_angle_steering_kingston_may2026`.

## Coda

> *γ₀ is the metronome.*
> *θ is what we do between beats.*
> *Q = tan θ is how far we rotate per beat.*
> *At θ = 0, we hear only the bare tick.*
> *At θ undefined, even the rhythm dissolves and only decay remains.*

> Live now: we still do not know what time is; that, Tom says, is why there are two clocks. The
> two hands are a typed object, `ClockHandLadderClaim`, read live by `inspect --root clock`: a Takt
> hand for the fading, 2γ₀ for one unit of light, and a coherence hand for the turning, with the
> angle between them. MirrorWorld's `Clock`, given N, owns both hands and that angle. Where the
> angle closes the mode fades at 2γ₀ times its mean light, and a chain gears the angle by its modes
> ("Where the angle closes" and "The repository's Q in a chain"). Closed forms: [`ANALYTICAL_FORMULAS.md`](../docs/ANALYTICAL_FORMULAS.md),
> F2b corollary.

---

**Anchors:**

- Universal Carrier typed claim: [`UniversalCarrierClaim.cs`](../compute/RCPsiSquared.Core/Symmetry/UniversalCarrierClaim.cs)
- F95 typed claim: [`F95AngleAtQuadraticZero.cs`](../compute/RCPsiSquared.Core/Symmetry/F95AngleAtQuadraticZero.cs)
- F95 proof: [`PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md`](../docs/proofs/PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md)
- F86 t_peak inheritance: [`F86TPeakPi2Inheritance.cs`](../compute/RCPsiSquared.Core/Symmetry/F86TPeakPi2Inheritance.cs)
- Absorption theorem: [`AbsorptionTheoremClaim.cs`](../compute/RCPsiSquared.Core/Symmetry/AbsorptionTheoremClaim.cs) (every mode fades at 2γ₀ times its light content; the integers that multiply γ₀)
- The two clocks: [`ClockHandLadderClaim.cs`](../compute/RCPsiSquared.Core/Symmetry/ClockHandLadderClaim.cs), live as `inspect --root clock`, and MirrorWorld's [`Clock.cs`](../compute/MirrorWorld/Clock.cs) given N; the horizon and the meeting at a pair's mean light: [`CoherenceHorizonClaim.cs`](../compute/RCPsiSquared.Core/Symmetry/CoherenceHorizonClaim.cs); the two-site floor in detail: [PROOF_CHAIN_GAP_DOMINANCE](../docs/proofs/PROOF_CHAIN_GAP_DOMINANCE.md) §4.2
- Kingston phase-injection run: [`data/ibm_f95_angle_steering_may2026/`](../data/ibm_f95_angle_steering_may2026/README.md) raw Heron r2 data + run writeup; `Confirmations.lookup('f95_angle_steering_kingston_may2026')` holds the numbers; [`F113_T1_EXTRACTION_KINGSTON.md`](../experiments/F113_T1_EXTRACTION_KINGSTON.md) reads the chip's own detuning off the same data
- Companion reflections, all today:
  - [`ON_HOW_THE_CARRIER_SHOWS_ITSELF.md`](ON_HOW_THE_CARRIER_SHOWS_ITSELF.md): γ₀ as Maßstab, the seams as extraction routes
  - [`ON_HOW_THE_ANGLE_APPEARS_AT_ZERO.md`](ON_HOW_THE_ANGLE_APPEARS_AT_ZERO.md): θ as the second coordinate above the discriminant zero
  - [`ON_HOW_FOUR_THIRDS_APPEARED.md`](ON_HOW_FOUR_THIRDS_APPEARED.md): F94's leading (4/3)·Q²·K³, the day's other closed form, at the per-outcome layer
- Readings of the same clock: [`Q_AS_THE_EXCHANGE_RATE.md`](../hypotheses/Q_AS_THE_EXCHANGE_RATE.md) (the next day: Q as the exchange rate between the two clocks, against taking γ₀ as the primary clock); [`ON_WHOSE_TIME_THE_CLOCK_KEEPS.md`](ON_WHOSE_TIME_THE_CLOCK_KEEPS.md) (the two hands built as a voice on the live system); [`NAVIGATING_THE_DIMENSIONS.md`](../docs/NAVIGATING_THE_DIMENSIONS.md) (the exceptional-point axis: the angle born from zero, at F86's meeting whose time constant is t_peak = 1/(4γ₀); it reads a real pair as θ = 0, where F95 leaves the angle undefined); [`THE_GENESIS_OF_AN_OSCILLATION.md`](../docs/THE_GENESIS_OF_AN_OSCILLATION.md) and [`Q_SCALE_THREE_BANDS.md`](../experiments/Q_SCALE_THREE_BANDS.md) (Q as the scale, γ₀ its unit)
- The inside reading: [Primordial Qubit §9](../hypotheses/PRIMORDIAL_QUBIT.md#9-operational-test-from-the-inside-april-15-2026) (the inside observer reads Q only)
- Canonical CAPS anchor docs this reflection extends:
  - [`MIRROR_THEORY.md`](../MIRROR_THEORY.md): the standing-wave reading; γ₀ is the tick of the standing wave we are
  - [`THE_BRIDGE_WAS_ALWAYS_OPEN.md`](../docs/THE_BRIDGE_WAS_ALWAYS_OPEN.md): γ as the source of time's arrow; the K = γ·t pure-number identity, constant within one book (F14), reads today as "K counts γ₀-ticks, not seconds"
  - [`ON_TWO_TIMES.md`](ON_TWO_TIMES.md): γ₀-time vs felt-time; today's tick reading sharpens the γ₀-time side (the metronome that flows but is not directly felt)
  - [`ON_THE_NINETY_DEGREE_GAMMA.md`](ON_THE_NINETY_DEGREE_GAMMA.md): the same γ₀ read through its distribution over the sites (F91)
  - [`PRIMORDIAL_GAMMA_CONSTANT.md`](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md): γ₀ as framework constant; the tick reading takes that constant as its unit

---

*Tom and Claude, 2026-05-16 (late evening). Tom's one-line compression closed the loop after a day spent unfolding F94 and F95. The metronome was there in every formula all along; today we heard it as a tick.*
