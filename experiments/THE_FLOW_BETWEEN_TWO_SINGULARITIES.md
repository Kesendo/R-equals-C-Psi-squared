# The Flow Between Two Singularities

**Date:** 2026-05-30
**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Status:** SEEN, not yet understood. Parked for future us. The numbers below are
verified (machine precision / bit-exact); the *interpretation* is deliberately left open.
We showed what happens; we do not yet grasp why. This document holds the structure so a
later session can pick it up without re-deriving, and resists the temptation to sink into
guessing while the numbers are in front of us but not yet readable.

**Probes:** [`the_flow_endpoints.py`](../simulations/the_flow_endpoints.py),
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
- The continuation past the target: sliding the net dephasing Σγ from Nγ₀ (target fades)
  through 0 (eternal, the mirror) into gain (the Hopf, runaway), the memory's three fates,
  is in [`memory_fate_across_the_takt.py`](../simulations/memory_fate_across_the_takt.py) and
  the [FRAGILE_BRIDGE](../hypotheses/FRAGILE_BRIDGE.md). The local EP only points there.

## For future us

Start by re-running the five probes above. The verified facts are the two singularity types
and the 1/N target. The open work is the *reading*: the polarity-0 / Σγ-mirror-0 seam, the
meaning of the journey, and whether the EP and the target are two faces of one structure or
genuinely two. Seeing was enough for today; the grasp is the next session's.

---

*Seen 2026-05-30. The numbers are on the table; reading them is still ahead.*
