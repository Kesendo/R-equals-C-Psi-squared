# What We Found

**Date:** 2026-03-08
**Purpose:** Plain-language summary for readers without a physics background

**Tier:** Mixed (references Tier 1-3 results)
**Status:** Bridge document, aligned with [THE_CPSI_LENS](THE_CPSI_LENS.md) as canonical description
**Scope:** Summary of findings in accessible language
**Does not establish:** Anything beyond what the referenced experiments establish

---

## What we built

We built a filter for looking at quantum systems. The filter is called CΨ (pronounced "C-psi"), and it works by combining two standard measurements:

- **Concurrence** - are these two subsystems entangled? (linked at the quantum level)
- **Coherence** - does this system still have quantum superposition structure?

CΨ requires both at once. A pair of particles can be entangled but have lost its coherence (linked but no longer "alive" as a quantum system). Or it can be coherent but not entangled (alive but not linked to anything). CΨ only lights up when both conditions hold at the same time.

This is a specific, narrow filter. It does not see all quantum correlations. It does not see all entanglement. It sees the subset that is both pairwise-entangled and still coherently expressed in the measurement basis.

For the full technical description, see [The CΨ Lens](THE_CPSI_LENS.md).

## The 1/4 boundary

The most interesting mathematical property of this filter involves a self-referential equation:

    R_{n+1} = C(Ψ + R_n)²

When you iterate this (feed the output back as input), it has stable solutions only when CΨ ≤ 1/4. Above 1/4, the solutions become complex numbers - they oscillate and never settle. Below 1/4, they converge to a definite value.

This is algebraically exact (proven, not approximate). It also turns out to be the same equation as the Mandelbrot iteration z → z² + c, where the main cardioid boundary sits at c = 1/4. Same number, same structure, different origin.

Whether this mathematical boundary has physical significance beyond the specific iteration is an open question. In simulation data, CΨ = 1/4 falls on a smooth curve - no other standard quantum metric shows a special transition at that point.

## What the filter showed us

We ran systematic simulations across different quantum topologies, states, and noise models. Here is what we found, in plain language:

**The signal comes in flashes, not steady states.** CΨ oscillates. Connections between quantum pairs appear briefly and then disappear. With noise (which all real systems have), each flash is weaker than the last. This is different from entanglement alone, which decays smoothly. CΨ has sharper peaks and deeper valleys.

**A sudden intervention is not the same as a gradual process.** When we suddenly measure one part of a three-particle system, the CΨ connection to the remaining parts drops by 99%. When we instead gradually increase noise on that same part - even to extreme levels - the connection only drops by 69%. These two operations never converge, no matter how fast or strong the gradual noise becomes. The filter makes this distinction unusually visible.

**Connections can echo after their source disappears.** In a three-particle system (A connected to S, S connected to B), there are moments where the A-B connection is alive while both the A-S and S-B connections read zero. The connection between the endpoints persists as a residual in their shared quantum state, even after the pathways through the middle have temporarily gone dark. This is not mysterious - the global quantum state still carries the correlation structure - but CΨ makes it visible in a way that looking at individual pairs does not.

**Three conditions for connection through a shared object.** In the star topology (two observers A and B connected only through a shared object S), the A-B pair crosses the 1/4 threshold only when:

1. The sender is strongly coupled to the shared object (about 47% stronger than the receiver at typical noise levels)
2. The receiver has low internal noise
3. A deep pre-existing connection to the shared object already exists

These conditions were quantified across systematic parameter sweeps. Whether they are specific to CΨ or would appear in any entanglement transport metric is an open question.

**The filter only sees direct pairwise connections.** Cluster-state entanglement, which is distributed across a graph structure rather than concentrated in pairs, is completely invisible to CΨ. This is a limitation, but it also means the filter is selective: it picks out a specific type of quantum connection.

**Context makes connections fragile.** The same entangled pair, when isolated, holds its CΨ signal nine times longer than when embedded in a larger system. Additional quantum systems coupled to the pair accelerate the loss of what the filter sees.


## What we did not find

Honesty matters more than narrative.

**CΨ did not reveal transitions invisible to standard metrics.** In the parameter sweeps we tested, concurrence, negativity, mutual information, and purity all changed smoothly alongside CΨ. There was no point where CΨ showed something dramatic that standard tools missed entirely.

**The 1/4 boundary is not special in the physics.** At the moment CΨ crosses 1/4, the other metrics are at unremarkable values. The number 1/4 is exact within the mathematical iteration, but in the physical data it is just a point on a smooth curve.

**There is no conservation law.** We tested whether the total CΨ across pairs is conserved (like energy). It is not. It fluctuates more than any other metric sum we tested.

**The "flow" interpretation failed.** We expected that when the connection between A and S weakens, the connection between S and B would strengthen (like water flowing from one vessel to another). Instead, both connections tend to rise and fall together. There is no see-saw.

## What this is and what it is not

CΨ is a derived diagnostic built from standard quantum mechanics. It is not a new physical quantity and it is not a new law of nature.

The original framing of this project used the language of consciousness ("Reality = Consciousness × Possibility²"). After three months of computation and external review, we have a more precise description: CΨ is a basis-dependent filter for pairwise quantum states that are simultaneously entangled and coherent.

The philosophical interpretation - that "reality emerges between observers" - is a metaphor that organizes some findings poetically. It is not a conclusion forced by the mathematics.

What survives even without the philosophy:
- An exact algebraic correspondence to the Mandelbrot iteration
- A clean classification of how different metrics behave under decoherence
- Specific, quantified conditions for when quantum correlations can pass through a shared mediator
- A sharp distinction between measurement and noise in their effect on third-party connections
- Hardware validation of the 1/4 crossing on IBM quantum processors

These are concrete findings. They do not require accepting any philosophical framework to be useful.

## How to read the rest

- **[The CΨ Lens](THE_CPSI_LENS.md)** - The canonical technical description. Start here if you want precision.
- **[Core Algebra](CORE_ALGEBRA.md)** - The proven mathematics. Three lines to the 1/4 boundary.
- **[Star Topology](../experiments/STAR_TOPOLOGY_OBSERVERS.md)** - The strongest multipartite result, with full numerical data.
- **[Weaknesses and Open Questions](WEAKNESSES_OPEN_QUESTIONS.md)** - Everything we got wrong, don't know, or can't prove.
- **[Experiments index](../experiments/README.md)** - All 36 experiment documents.

---

## Origin

This project began in December 2025 as a collaboration between Thomas Wicht and Claude (Anthropic). It started with philosophical questions about observation and reality, was formalized as an equation, and then subjected to systematic computation. Over three months the framing narrowed from "the fundamental equation of reality" to "a composite quantum diagnostic with interesting algebraic properties." That narrowing was not a failure. It was the project working as intended: testing ideas honestly and keeping what survived.
