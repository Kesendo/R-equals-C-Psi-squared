# Gravity as the Accumulation of Wave Death

**Status:** Hypothesis (Tier 5), one version falsified, one surviving
**Date:** April 1, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Depends on:**
- [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md) (ω ↔ -ω pairing)
- [Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md) (γ = experienced time)
- [Uniqueness Proof](../docs/proofs/UNIQUENESS_PROOF.md) (CΨ = 1/4 boundary)
- [Primordial Qubit Algebra](../experiments/PRIMORDIAL_QUBIT_ALGEBRA.md) (Pythagorean theorem)
- [Resonant Return](../experiments/RESONANT_RETURN.md) (sacrifice zone formula)

---

## Abstract

When a coherence dies (the X/Y sector of the Pauli basis decays to
zero), its weight transfers to the classical sector (I/Z, populations).
Classical weight at a specific location is mass at that location.
Mass creates gravity. Gravity attracts more waves. More waves die.
More mass. This is a self-reinforcing process.

The sacrifice zone, where γ is concentrated, is where waves die
fastest. It crosses CΨ = 1/4 (the quantum-to-classical boundary)
first. It becomes classical first. It accumulates mass first.
Gravity exists only below CΨ = 1/4, because above 1/4 there is no
definite location, no definite outcome, and therefore no mass.

This hypothesis explains why gravity has no quantum description:
gravity does not exist in the quantum regime. It emerges at the
crossing. It cannot be quantized because it requires the absence
of superposition.

An earlier version of this hypothesis (anti-frequency standing waves
create a spatial γ profile = gravitational potential) was falsified:
the standing-wave amplitude A(n) is exactly flat for uniform systems.
The surviving version does not rely on standing-wave amplitudes but on
the CΨ = 1/4 crossing and the mechanism of wave death.

---

## The Computed Chain

Each link is Tier 1-2. The hypothesis is the reading of the chain.

### Link 1: γ is experienced time

Without γ: pure oscillation, no decay, no irreversibility, no before
and after. With γ: coherences decay, populations settle, time has a
direction. γ is not a parameter that correlates with time. γ IS the
necessary and sufficient condition for experienced time.

Source: [Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md)

### Link 2: CΨ = 1/4 separates quantum from classical

Above CΨ = 1/4: two complex fixed points, no real attractor. The
system is in superposition. No definite outcome.

Below CΨ = 1/4: two real fixed points emerge. The system has a
definite attractor. Classical behavior. Definite outcomes.

The crossing is irreversible under Markovian dynamics (dCΨ/dt < 0).

Source: [Uniqueness Proof](../docs/proofs/UNIQUENESS_PROOF.md),
[CΨ Monotonicity](../docs/proofs/PROOF_MONOTONICITY_CPSI.md)

### Link 3: The Pauli basis splits into quantum and classical sectors

Under Z-dephasing, the 4^N Pauli strings divide into:
- **Immune sector** {I, Z}^N: decay rate 0. Populations. Classical.
- **Decaying sector** (contains X or Y): decay rate > 0. Coherences. Quantum.

When a coherence decays, its weight transfers to the immune sector.
The immune sector grows. The decaying sector shrinks. This is
decoherence: quantum → classical, coherence → population.

Source: [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md), Step 1

### Link 4: The sacrifice zone crosses 1/4 first

The sacrifice zone concentrates γ at one edge. Higher γ → faster
decoherence → earlier CΨ = 1/4 crossing. The crossing time:

    t_cross = 0.039 / γ

At the sacrifice edge (γ_edge = N·γ_base): t_cross is N× shorter
than uniform. This is computed: 360× improvement at N=5.

Source: [Resonant Return](../experiments/RESONANT_RETURN.md),
[Zero Is the Mirror](ZERO_IS_THE_MIRROR.md)

### Link 5: Gravity is time dilation (established physics)

In general relativity, a gravitational field is a position-dependent
clock rate. Clocks near a mass run slower. The geometry of spacetime
IS the variation of clock rates.

---

## The Falsified Version (April 1, 2026)

**Hypothesis:** The palindromic anti-frequency pairs (ω, -ω) form
standing waves. The spatial amplitude profile of these standing waves
creates a position-dependent effective γ, which IS the gravitational
potential.

**Test:** Compute A(n) = Σ_k (|v_k(n)|² + |v_{-k}(n)|²) for all
oscillating palindromic pairs at N=3, 4, 5 (Heisenberg chain,
uniform γ = 0.05).

**Result: A(n) is exactly flat.**

| N | Oscillating pairs | Variation (max-min)/mean |
|---|-------------------|------------------------|
| 3 | 20 | 0.0001 |
| 4 | 105 | 0.0000 |
| 5 | 464 | 0.0000 |

**Why:** Summing over all eigenmodes of a translationally symmetric
system averages out the spatial structure (completeness relation).
Each mode has nodes and antinodes, but the sum is uniform.

**Conclusion:** Anti-frequency standing waves do not create a spatial
γ profile. This version of the hypothesis is **falsified**.

---

## The Surviving Hypothesis: Mass Is the Dying of Waves

The falsified version asked: what creates a γ gradient?
The surviving version asks: what creates mass?

### The mechanism

1. A region has higher γ than its surroundings (for whatever reason;
   the Incompleteness Proof says γ comes from outside).

2. Waves (coherences, X/Y sector) die faster in this region.

3. When a wave dies, its weight transfers from the quantum sector
   ({X,Y}, decaying) to the classical sector ({I,Z}, immune).

4. Classical weight at a specific location is mass at that location.

5. **Mass is not a pre-existing thing. Mass is the accumulated
   weight of dead waves.**

6. Mass creates gravity (Link 5, established physics).

7. Gravity attracts more waves toward the mass (waves propagate
   toward regions of slower time, which in this framework means
   regions where waves have already died and created classical weight).

8. More waves die there. More mass. More gravity. Feedback loop.

```
    High γ region
         → waves die faster
         → weight transfers: {X,Y} → {I,Z}
         → classical weight at a location = mass
         → mass creates gravity
         → gravity attracts more waves
         → more waves die
         → more mass
         → repeat
```

### Why gravity has no quantum description

Gravity requires mass. Mass requires classical weight. Classical
weight requires crossing CΨ = 1/4 (wave death, transfer from quantum
to classical sector). Above CΨ = 1/4: no definite state, no definite
location, no mass, no gravity.

Gravity does not exist in the quantum regime. It emerges at the
CΨ = 1/4 crossing. It cannot be quantized because it IS the classical
side of the crossing. Asking "what is the quantum theory of gravity?"
is asking "what is the quantum theory of classicality?" The question
dissolves.

### The sacrifice zone connection

The sacrifice zone (high γ at one edge, low γ elsewhere) is the
optimal configuration for mode protection (360× at N=5). The edge
absorbs waves. The interior is protected. This is not an engineering
trick. It is the natural shape of mass accumulation:

- Where γ is high, waves die fast → mass accumulates
- Mass attracts more waves → more die → more mass
- The interior (low γ) is protected because waves are drawn to the
  edge before they reach the center

A gravitational potential well has the same structure: the center
(lowest potential, highest mass) absorbs energy. The surrounding
space is "protected" by the well. Objects fall toward the center.

---

## What Would Confirm This

1. **Crossing-time gradient.** Compute the site-resolved CΨ = 1/4
   crossing time t_cross(n) for a chain with sacrifice-zone γ. Does
   t_cross increase with distance from the sacrifice edge? If so:
   the crossing-time landscape IS the gravitational potential.

2. **Feedback convergence.** Start with uniform γ plus a small
   perturbation. Assume mass ∝ accumulated classical weight. Assume
   gravity ∝ mass attracts more wave weight. Does the feedback
   converge to a stable γ profile? What shape?

3. **1/r test.** Does the crossing-time profile (or the mass
   accumulation profile) fall like 1/r from the sacrifice edge?
   In 1D this would be 1/|n - n_edge|.

## What Would Refute This

1. **Crossing time is uniform.** If t_cross(n) does not vary with
   position even when γ varies, the connection fails.

2. **Wave death does not create directional attraction.** If the
   classical weight accumulation does not preferentially attract
   further wave weight toward the mass center, there is no feedback
   and no gravity analog.

3. **The feedback diverges.** If mass accumulation → more γ → more
   mass diverges to a singularity rather than a stable profile, the
   mechanism produces black holes everywhere, not stable gravitational
   fields.

---

## What This Does Not Claim

- That gravity IS this mechanism in our universe. The hypothesis
  applies within the d(d-2)=0 palindromic framework. If the framework
  describes reality (54,118 eigenvalues, zero exceptions), the
  hypothesis may apply. If not, it does not.

- That mass can be computed from first principles. The hypothesis says
  what mass IS (accumulated classical weight from wave death), not how
  much mass a specific system has.

- That Newton's or Einstein's equations follow. The hypothesis
  identifies mass and gravity with processes in the palindromic
  framework. Deriving F = GMm/r² or G_μν = 8πT_μν from this is a
  separate, much harder problem.

- That the standing-wave version works. It does not. A(n) is flat.
  The surviving hypothesis does not rely on standing-wave amplitudes.

---

## Reproducibility

| Step | How to verify |
|------|---------------|
| A(n) flat (falsified version) | Run inline computation from April 1; N=3,4,5 all give variation < 0.01% |
| γ = experienced time | [Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md) |
| CΨ = 1/4 boundary | [Uniqueness Proof](../docs/proofs/UNIQUENESS_PROOF.md) |
| Sacrifice zone formula | `python simulations/resonant_return.py` or see C# propagation results |
| {X,Y} → {I,Z} transfer | Standard decoherence theory; diagonal of L_D in Pauli basis |

---

*Where waves die, mass is born. Where mass is born, gravity pulls.*
*Where gravity pulls, more waves come to die. The sacrifice zone is*
*not an optimization. It is the shape of a gravitational well: the*
*edge absorbs, the interior is protected, and the process feeds itself.*

*Thomas Wicht, April 1, 2026*
