# The Flow Between Two Singularities

**Date:** 2026-05-30
**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Status:** SEEN, not yet understood. Parked for future us. The numbers below are
verified (machine precision / bit-exact); the *interpretation* is deliberately left open.
We showed what happens; we do not yet grasp why. This document holds the structure so a
later session can pick it up without re-deriving, and resists the temptation to sink into
guessing while the numbers are in front of us but not yet readable.

**Probes:** [`the_flow_endpoints.py`](../simulations/the_flow_endpoints.py),
[`at_the_target.py`](../simulations/at_the_target.py),
[`post_ep_dynamics_4d.py`](../simulations/post_ep_dynamics_4d.py),
[`f86_ep_through_the_clock.py`](../simulations/f86_ep_through_the_clock.py),
[`after_the_collapse.py`](../simulations/after_the_collapse.py),
[`memory_fate_across_the_takt.py`](../simulations/memory_fate_across_the_takt.py)

---

## What we set out to see

F86's exceptional point, read through the clock (the Takt/Rotation reading on
MirrorSystem), gives the dynamics a *birth*: forgetting (real decay) collapses into the
EP, and above it a rotating mode opens (the memory, F80's imaginary axis). See
[F86_EP_THROUGH_THE_CLOCK](F86_EP_THROUGH_THE_CLOCK.md). The story does not stop there. We
asked what happens *after*: where does the reborn rotating mode go? That is a dynamical
question, and the static spectrum collapses it. It is only readable in the loop (the RK4
trajectory), so we computed it.

## What the loop showed

**1. The post-EP dynamics is a 4D object (Q × site × t).** Below the rotation onset
(small Q ~ J/γ₀) the single excitation just diffuses and decays (forgetting). Above it,
the excitation sloshes site to site, a wave that propagates and reflects (remembering).
Edge and bulk sites differ. The dimensionless knob is Q alone.

**2. The target is 1/N, the equipartitioned state.** Every trajectory, at every Q, ends at
⟨n_site⟩ = 1/N (verified N=3,4,5,6 → 0.3333, 0.2500, 0.2000, 0.1667). This is the
fully-forgotten, uniform single-excitation state. It equals ¼ at N=4 *by coincidence*; it
is NOT the framework's universal CΨ=¼ fold (that is N-independent coherence, a different
quantity; the steady coherence here is 0).

**3. Two singularities bracket the flow, and they are different types** (verified in
`the_flow_endpoints.py`, N=3..6):

| | the **EP** (birth) | the **target** (death) |
|---|---|---|
| where | parameter space, at Q_EP | state space, the fixed point |
| type | **DEFECTIVE** (Jordan block, rank-deficient, eigenvectors coalesce, Petermann K → ∞) | **SIMPLE** (λ=0, geometric = algebraic multiplicity) |
| L | rank-deficient at that Q | singular (det L = 0, the kernel) |
| memory | is born | dies (forgotten) |

The target sits in the kernel of L: `L · vec(1/N uniform) = 0` bit-exact. The kernel is
**(N+1)-dimensional**, one fixed point per particle-number sector; the single-excitation
trajectory lands on its own sector's fixed point.

## The structural picture (seen, not interpreted)

The flow runs from a defective singularity (the EP, where the memory is born) to a simple
singularity (the target, the λ=0 fixed point, where the memory dies). The EP is a *pinch*;
the target is a *fixed point*. The journey between them is the loop, and the loop has no
closed form: that is exactly where the dynamics lives. So the clean, closed-form structure
sits at the two singular *endpoints* (the EP's defectiveness; the target's kernel), and the
un-closed-form dynamics is the *middle*. This reframes the long-standing "F86 has no closed
form" not as a gap but as a location: the closed form ends where the loop begins.

## At the target: a point, and then?

The target is one point. What happens there, and how it goes on, is computed in
[`at_the_target.py`](../simulations/at_the_target.py).

**It is a sink, reached only asymptotically.** Every non-kernel mode of L has Re λ < 0, so
all nearby states flow in (max non-kernel Re λ = −2 at N=4, Q=2). But the approach is
exponential: |⟨n⟩ − 1/N| halves at a fixed rate and reaches 0 only in the limit. The point
is approached forever, never landed in finite time.

**The future is already present, and conserved.** Decompose ρ(0) (a localized excitation)
in L's eigenmodes: the λ=0 component already equals the 1/N target, bit-exact (‖difference‖
= 1.4e-15), at t=0. Its coefficient is exp(0·t) = 1, so it never grows and never decays. The
target is not built up over time; it is the invariant part of the state, present from the
start at full strength. What changes is everything else: the non-kernel (transient) part,
norm 0.866 at t=0, fades to zero and *reveals* the target that was always already there. The
loop does not carry us to the target; it dissolves the present to expose it.

**How it goes on: the point's role flips across the mirror.** The 1/N point stays the fixed
point for *every* net dephasing Σγ (`L · vec(1/N) = 0` throughout); only its stability
changes:

| Σγ | max non-kernel Re λ | role |
|---|---|---|
| > 0 | < 0 | **sink** (the flow falls in; the end of history) |
| = 0 | = 0 | **neutral** (eternal oscillation; the mirror) |
| < 0 | > 0 | **source** (the flow is pushed away; the Hopf, runaway) |

The point does not move; the arrows around it invert. And there, at the dynamical layer,
Tom's {−, 0, +} reappears: not as the operator polarity triple (it did not fit there), but
as the *sign of the spectral gap*, the three fates of the fixed point (sink −, mirror 0,
source +) across the Σγ axis of [ZERO_IS_THE_MIRROR](../hypotheses/ZERO_IS_THE_MIRROR.md).

## Tier-5 reading: the target as the future (our motor)

*Below this point is interpretation, in the spirit of [THE_OTHER_SIDE](../hypotheses/THE_OTHER_SIDE.md)
§14: the structure above is computed (Tier 1-2); the reading here is Tier-5. We keep it
because, as Tom put it on 2026-05-30, it is our motor and drive, not a truth-claim.*

The target has the shape of the future as we know it: a point you approach but never reach
(the moment you would arrive, it is the present, not the future), the thing the arrow of
time points at (the sink). And the computation says something stranger and lovelier, the
future already exists. It is the conserved λ=0 mode, present in ρ(0) at full strength from
t=0, invariant. We do not travel to it; the present fades and uncovers it. In the timeless
spectrum (the Compute engine) all modes, including the future, are there at once; the loop
(Propagate) is only the unfolding.

This is the computed face of the oldest wild thing in the project: "Both sides exist
simultaneously", R = CΨ² and Ψ = √(R/C), "neither first, both simultaneous"
([THE_OTHER_SIDE](../hypotheses/THE_OTHER_SIDE.md) §3, §17). What was poetry with eigenvalues
is now a conserved fixed point we can touch. Two honest seams stay at the label, not the
structure: today's "future" is the *decided* equilibrium (the kernel), where THE_OTHER_SIDE
called the future the *undecided* coherences, opposite labels but the same simultaneity; and
the mirror-world reading itself stays Tier-5, where that document drew its own line. We
returned with new sight and the wild thing sharpened into something we could compute. We did
not invent it; we learned to see the structure that was already there.

## What we do NOT yet understand (for future us)

- Why two singularities, and why of *different* types (one defective, one simple)?
- What does the journey between them *mean*, beyond "relaxation to a fixed point"?
- Does the EP → target flow carry structure we have not named (a conserved quantity along
  the way, a geometry of the approach)?
- The kernel is (N+1)-dimensional, one fixed point per number sector. Does that per-sector
  structure connect to anything, or is it just number conservation?
- **Tom's polarity thread (2026-05-30), recorded as an open viewpoint.** The {−0, 0, +0}
  polarity layer at d=0 and its d=2 reading {−½, 0, +½} (the
  [PolarityLayerOrigin](../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs)
  0.5-shift, deepened in [POLARITY_COORDINATES](../reflections/POLARITY_COORDINATES.md) and
  [THE_POLARITY_LAYER](../hypotheses/THE_POLARITY_LAYER.md)) may touch this problem at the
  **target**: the target (1/N, uniform) is the *unpolarized* state, and the flow to it is a
  *depolarization*, the localized excitation (polarized, "alive") relaxing back to the
  uniform 0 (the d=0 layer, the "forgotten"). So the {−½, 0, +½} → {−0, 0, +0} collapse
  reads as the return to the unpolarized layer. The honest seam: there appear to be *two
  different zeros*. The **polarity-0** (unpolarized = the target, reached at fixed Σγ>0) and
  the **Σγ-mirror-0** of [ZERO_IS_THE_MIRROR](../hypotheses/ZERO_IS_THE_MIRROR.md) (Σγ=0, no
  decay, the eternal standing wave, reached by *sliding* Σγ) are orthogonal axes, not the
  same point. And the EP does not sit on the polarity triple at all, it is a
  dynamics-layer singularity. Whether the operator-level polarity triple (the residual M's
  Π +i/0/−i content) and the state-level polarization (localized vs uniform) are the same
  structure seen from two sides is itself open. Pursue as a viewpoint; do not force a single
  three-point map.
- The continuation past the target is now partly seen as the role-flip above
  (sink/neutral/source); the memory's three fates along Σγ are in
  [`memory_fate_across_the_takt.py`](../simulations/memory_fate_across_the_takt.py) and the
  [FRAGILE_BRIDGE](../hypotheses/FRAGILE_BRIDGE.md). What stays open is whether the three
  fates and the EP's birth are one structure.

## For future us

Start by re-running the six probes above. The verified facts are the two singularity types,
the 1/N target, the future-already-present (the conserved λ=0 mode), and the sink/neutral/source
role-flip. The open work is the *reading*: the polarity-0 / Σγ-mirror-0 seam, the
meaning of the journey, and whether the EP and the target are two faces of one structure or
genuinely two. Seeing was enough for today; the grasp is the next session's.

---

*Seen 2026-05-30. The numbers are on the table; reading them is still ahead.*
