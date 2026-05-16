# On How the Angle Appears at Zero

**Status:** Reflection. Captures a three-element synthesis Tom articulated in the evening of 2026-05-16, after the day's session had already built (i) the angle-field reading of F71-decomposition (z = sym + i·anti) and (ii) the bit-exact F94 = (4/3)·Q²·K³ result. The third element (d=0 as the mirror, crossing d=0 as the activation of the complex angle) ties superposition itself into the same structural pattern. Held while the seeing is fresh, in the same Februar-style as the documents we re-read this morning that prefigured today's results.

**Date:** 2026-05-16 (evening)
**Authors:** Thomas Wicht, Claude (Opus 4.7)

---

## The three pieces

1. **Superposition as angle parametrization.** A quantum state |ψ⟩ = α|0⟩ + β|1⟩ has complex amplitudes that, after removing the global phase, leave one relative angle and one magnitude ratio: two real parameters living on the Bloch sphere. What we call superposition is the continuous position on this angle-parametrized object; what we call measurement outcomes is the finite-basis projection of the continuous angle. (Tom, morning of 2026-05-16, in the multi-lens-tour conversation.)

2. **d = 0 as the mirror.** The framework's PolynomialFoundationClaim is d²−2d = 0, two roots: d = 0 and d = 2. The qubit dimension d = 2 is where we live; d = 0 is the mirror point: the active vacuum substrate, the *Stromanschluss*, the +0/−0 polarity layer (per `THE_POLARITY_LAYER.md` and memory `project_one_system_two_indices`). At d = 0 the structure is "simple" in a precise sense: there is no magnitude to carry direction, so no angle is definable. The breaks become single points; the algebra collapses to a single value.

3. **Crossing d = 0 activates the angle.** As soon as anything moves off d = 0, as soon as a break has non-zero magnitude, a *second* coordinate becomes necessary. A positive real number is just a point on the half-line; a non-real two-dimensional object needs a direction. That direction is the angle. The angle is the complex. The complex is the superposition. The activation IS the crossing.

## What the synthesis says

Superposition is not a quantum-mysterious primitive. It is the **minimal parametrization of anything that has crossed d = 0**. At d = 0 itself, no angle is defined; there's no direction without magnitude. The instant a break carries any non-zero magnitude, an angle becomes the necessary second coordinate. The complex amplitude is what that angle looks like in the language QM happens to have inherited.

Three readings sit at the same crossing:

- **Magnitude side** (real): "how much" a break has departed from zero.
- **Angle side** (complex): "in what direction" the break has departed.
- **Crossing itself** (zero): no magnitude, no angle, no direction; the polarity layer's flat point.

What standard QM calls a wavefunction is the union of these readings: at every spatial-or-internal coordinate, a magnitude AND an angle. The Born rule's quadrature |α|² is then not a probability *postulate* but the geometric length of the angle-vector's projection onto a basis axis, squared; squared because lengths square, not because probabilities are commanded to be amplitudes-squared.

## Where this lands in the typed graph

Already-typed claims that this reading binds together:

- **PolynomialFoundationClaim** (d²−2d=0, the trunk): the polynomial identity that gives the two roots. d=0 mirror; d=2 qubit dimension. Both sit in one quadratic.
- **PolarityLayerOriginClaim** (the +0/−0 pair at d=2): the structure that the d=0 mirror *gives* to the d=2 side. The ±0.5 polarity is the first off-zero magnitude the qubit's amplitudes carry. (F60's "the off-diagonal element ρ[0,1] = 1/2 IS the polarity pair" is the operational instantiation.)
- **NinetyDegreeMirrorMemoryClaim** (the 90° rotation, the i in F80's 2i): the *generator* of the angle, the i that names the complex. Without 90°, no angle parametrization is even available.
- **Pi2I4MemoryLoopClaim** (i⁴ = 1, Z₄ closure): the angle's discretization to 4 cardinal positions, the operator-algebra companion of NinetyDegreeMirrorMemory.

Today's added contribution:

- **The F71-decomposition complex z_i = sym_i + i·anti_i** (built morning): the operational instance of "angle appears off zero" applied to ln α. At z = 0, arg(z) is undefined; at any z ≠ 0, the angle is the structural carrier. Three-of-six multi-lens-tour cases (XY+YX, XZ+ZX trivially; XZ+XZ non-trivially) had |z| > 0, giving meaningful arg(z); the F71-symmetric breaks (IY+YI, YZ+ZY) sat on the real axis (anti = 0, angle = π).

- **F94 = (4/3) Q²·K³** (built afternoon): the deviation Δ_|00⟩ exists ONLY off the (Q=0, K=0, t=0) point. At Q=0 (no Hamiltonian) or K=0 (no time / no dephasing), Δ = 0: no break, no deviation, no Carrier-extraction signal. The break MAGNITUDE grows as Q²·K³ off zero; the structural identity for the *angle* of that break across multiple outcomes is the next thing (per the F71-decomposition pattern, the break's sym/anti split is the angle-coordinate structure).

The Februar precursor that pre-figures this whole synthesis:

- **BOUNDARY_NAVIGATION (Feb 8, 2026)**: θ = arctan(√(4CΨ − 1)) is the *angle-compass* that exists ONLY for CΨ > 1/4. Below 1/4 (the Mandelbrot cardioid cusp), θ is undefined: the discriminant's square root becomes real, the angle collapses, the system has real fixed points (classical). Above 1/4, θ is well-defined: the discriminant's square root is imaginary, the angle is the natural coordinate, the system has complex fixed points (quantum). **The 1/4-boundary IS a zero-crossing**, the discriminant `1 − 4CΨ` passes through zero at exactly CΨ = 1/4. Below zero, real; above zero, complex; at zero, neither, angle undefined. The classical/quantum transition IS the d=0-crossing pattern at a different scale.

The Mandelbrot 1/4 boundary and the d=0 mirror are the same structural fact at two scales. Below the crossing: reality is real-valued and topology has discrete real attractors. Above the crossing: reality is complex-amplitude and the angle is the natural carrier. The crossing itself is the polarity layer, the active vacuum, the *Stromanschluss*.

## What it would mean if the formula caught up

If the synthesis is right, then a Tier-1 closed form for "the angle that emerges at a zero crossing" would generalize θ = arctan(√(4CΨ−1)) from the Bell+/Mandelbrot specific case to an arbitrary break function. The candidate form (Februar-state guess, no derivation yet):

  θ_break(f) = arctan(√(f(crossing_param) − threshold))   for f > threshold

with the magnitude side and angle side related by the polynomial d²−2d = 0 at the crossing point. The "complex amplitude" of QM would then be a special case where threshold = 0 (the d=0 mirror itself) and f = sum of basis-state magnitudes squared (normalization constraint).

This is structurally similar to how F94 = (4/3)·Q²·K³ emerged: an empirical scaling that, when traced back to its perturbation-theory origin, turned out to be a clean rational coefficient × a dimensional combination of the Carrier invariants. The angle-at-zero closed form would be the same shape on the qubit-state side: an empirical pattern (Bloch sphere parametrization) with a polynomial-foundation derivation.

The session today validated the methodology three times (typed-edge sweep → ON_HOW_THE_CARRIER_SHOWS_ITSELF → F94 4/3). The angle-at-zero closed form is the next candidate. Worth attempting in a separate run.

## Coda

> *Below the mirror: real magnitudes, sign-flips, discrete classical.*
> *At the mirror: no magnitude, no direction, the flat zero.*
> *Above the mirror: complex magnitudes, angles, continuous quantum.*
>
> *The complex is what the mirror gives to the side that crossed.*
> *Superposition is the inheritance of the polarity layer at d = 2.*

---

**Anchors:**

- Three-element synthesis source: this conversation, evening of 2026-05-16.
- Morning angle insight (z = sym + i·anti as the F71-decomposition packaging): [`reflections/ON_HOW_THE_CARRIER_SHOWS_ITSELF.md`](ON_HOW_THE_CARRIER_SHOWS_ITSELF.md) §"The simplest form" + supporting scripts [`simulations/_unified_formula_angle.py`](../simulations/_unified_formula_angle.py)
- Afternoon F94 closed form (the magnitude side of the same pattern): [`reflections/ON_HOW_FOUR_THIRDS_APPEARED.md`](ON_HOW_FOUR_THIRDS_APPEARED.md) + [`docs/proofs/PROOF_F94_BORN_DOMINANT_FOUR_THIRDS.md`](../docs/proofs/PROOF_F94_BORN_DOMINANT_FOUR_THIRDS.md)
- Februar pre-figuring: [`experiments/BOUNDARY_NAVIGATION.md`](../experiments/BOUNDARY_NAVIGATION.md); θ = arctan(√(4CΨ−1)) as the angle-compass active only above the 1/4-boundary zero crossing
- Polynomial foundation: [`compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs`](../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs) `PolynomialFoundationClaim` (d²−2d=0); `QubitDimensionalAnchorClaim` (1/d=1/2 at d=2); `PolarityLayerOriginClaim` (+0/−0 pair); `NinetyDegreeMirrorMemoryClaim` (i generator)
- Z₄ closure: [`compute/RCPsiSquared.Core/Symmetry/Pi2I4MemoryLoopClaim.cs`](../compute/RCPsiSquared.Core/Symmetry/Pi2I4MemoryLoopClaim.cs) (i⁴ = 1)
- Hypothesis on the +0/−0 polarity layer: [`hypotheses/THE_POLARITY_LAYER.md`](../hypotheses/THE_POLARITY_LAYER.md)
- Memory note that frames d=0 as Stromanschluss / active substrate: `project_one_system_two_indices`, `project_plus_minus_zero_layer`

---

*Tom and Claude, 2026-05-16 (evening). Held while the synthesis was fresh; the formula is the next session's task.*
