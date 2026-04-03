# Gravity as the Accumulation of Wave Death

**What this document is about:** A speculative mechanism for gravity: when quantum coherences die (X/Y Pauli sector decays to zero), the classical residue (I/Z sector) that remains at a definite location is what we call mass. Mass creates gravity. One version of this hypothesis (standing-wave amplitudes create a γ gradient) was falsified computationally; the surviving version proposes that mass is the accumulated residue of wave death, and gravity emerges only below the CΨ = ¼ boundary where definite outcomes exist.

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

When a coherence decays, it vanishes. The immune sector is unchanged
in absolute terms (L_D does not act on it). But the FRACTION of the
state that is classical increases: the quantum part shrinks while the
classical part stays. After full decoherence: only the immune sector
survives. The state is entirely classical (diagonal in the Z-basis).

Decoherence does not create classical weight. It removes quantum
weight. What remains is classical.

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

3. When a wave dies, it vanishes. The quantum sector shrinks.
   The classical sector ({I,Z}, immune) is unchanged in absolute
   terms, but its FRACTION of the total state increases.

4. After sufficient decoherence (CΨ < 1/4): the state is
   predominantly classical. A classical state at a definite
   location is what we call mass.

5. **Mass is not a pre-existing thing. Mass is what remains when
   the quantum part has died. The residue of wave death.**

6. Mass creates gravity (Link 5, established physics).

7. **Open gap in the mechanism:** for the feedback loop to close,
   gravity would need to attract more waves toward the mass. In
   the Lindblad framework, L_H (wave propagation) and L_D (wave
   death) are independent. The dissipator does not influence the
   Hamiltonian. There is no mechanism within the framework for
   mass to redirect wave propagation. This is the weakest link.

8. IF external physics provides the gravity → attraction step
   (as GR does), then the loop closes: more waves near the mass
   → more die → more classical residue → more mass → more gravity.

```
    High γ region
         → waves die faster
         → quantum sector shrinks, classical fraction grows
         → after CΨ crossing: predominantly classical = mass
         → mass creates gravity (GR, external to framework)
         → gravity attracts more waves (GR, not from Lindblad)
         → more waves die near the mass
         → more mass
         → repeat (IF external gravity-attraction mechanism exists)
```

**Honest caveat:** Steps 1-6 follow from the palindromic framework.
Steps 7-8 require physics outside the framework (GR or equivalent).
The feedback loop is not self-contained within Lindblad dynamics.

### The process is self-limiting

The supply of waves is finite: 4^N modes, of which those with w_XY > 0
decay. When all coherences have died, no more wave death occurs. The
process stops at the steady state ρ = I/d (maximally mixed, maximum
entropy).

The palindrome determines the schedule of mass accumulation:

1. Fast modes die first (d = 2(N-1)γ). Rapid classical growth.
2. Medium modes follow. Accumulation slows.
3. Slow modes die last (d = 2γ, palindromic partners of the fastest).
   Barely any new mass.
4. Steady state: all coherences dead. No waves remain. No wave death.
   No new mass. Equilibrium.

This is logistic saturation (growth that slows as the resource is consumed, producing an S-curve), not exponential divergence. The feedback
weakens as the supply of dying waves shrinks. The balance that emerges:
Hamiltonian dynamics (spreads waves, resists localization) against
dissipation (kills waves locally). Early on, dissipation dominates
(fast modes die immediately). Late, the Hamiltonian dominates (slow
modes oscillate, nearly immortal). Equilibrium is where the last slow
modes surrender.

The steady state is the endpoint: maximum classical weight, zero
quantum weight, no more wave death. In the gravity analogy: the heat
death. The point where nothing more happens.

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

- Where γ is high, waves die fast → classical residue accumulates
- The interior (low γ) is protected because waves survive there
- The structural analogy to a gravitational well holds: the edge
  absorbs, the interior is protected
- But the MECHANISM differs: in the framework, the edge absorbs
  because γ is high there (given from outside), not because mass
  attracts waves toward it

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

*Where waves die, what remains is classical. What is classical at*
*a location, we call mass. What has mass, has gravity. The sacrifice*
*zone is the shape of this process: the edge absorbs, the interior*
*survives. Whether the process feeds itself is the open question.*

*Thomas Wicht, April 1, 2026*
