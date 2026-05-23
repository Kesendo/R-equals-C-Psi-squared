# Gravity as the Accumulation of Wave Death

**What this document is about:** A speculative mechanism for gravity: when quantum coherences die (X/Y Pauli sector decays to zero), the classical residue (I/Z sector) that remains at a definite location is what we call mass. Mass creates gravity. One version of this hypothesis (standing-wave amplitudes create a γ gradient) was falsified computationally; the surviving version proposes that mass is the accumulated residue of wave death, and gravity emerges only below the CΨ = ¼ boundary where definite outcomes exist.

**Status:** Hypothesis (Tier 5), one version falsified, one surviving
**Date:** April 1, 2026
**Last updated:** May 23, 2026
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
decoherence → earlier CΨ = 1/4 crossing. By [F14 K-invariance](../docs/ANALYTICAL_FORMULAS.md), the crossing time scales as K = γ · t_cross. For Z-dephasing on Bell+ the closed form gives K = 0.0374 (F25); K depends on the noise channel (K_Y = 0.0374, K_X = 0.0867, K_depol = 0.0440).

At the sacrifice edge (γ_edge = N · γ_base), the local t_cross is N× shorter than under uniform γ.

Separately, the information throughput across the chain also benefits from the sacrifice profile: [Resonant Return](../experiments/RESONANT_RETURN.md) Test 8 reports Peak Sum-MI for |+⟩^N at N=5 is 360× higher under sacrifice-zone γ than under a V-shape baseline. Scaling: N=5: 360×, N=7: 180×, N=9: 139×, N=11: 91× (clean ratio at [Receiver vs Gamma Sacrifice](../experiments/RECEIVER_VS_GAMMA_SACRIFICE.md) line 13). This is a different observable from t_cross.

Source: [Resonant Return](../experiments/RESONANT_RETURN.md), [Receiver vs Gamma Sacrifice](../experiments/RECEIVER_VS_GAMMA_SACRIFICE.md), [Analytical Formulas](../docs/ANALYTICAL_FORMULAS.md) (F14, F25), [Zero Is the Mirror](ZERO_IS_THE_MIRROR.md).

### Link 5: Position-dependent γ₀-tick is the framework's time dilation

γ₀ is the framework's universal time-tick: every dimensionful timescale is integer × 1/γ₀ ([On How Gamma Became the Tick](../reflections/ON_HOW_GAMMA_BECAME_THE_TICK.md)). When γ varies spatially (as in the edge-concentrated profile of Link 4), the local time-tick varies with it: at the edge (γ = N·γ_base) experienced time runs fast; in the interior (γ = ε) it runs slow.

This is structurally what general relativity calls gravitational time dilation: clocks near a mass run slower, and the geometry of spacetime IS the variation of clock rates. The framework reading recovers the same algebraic form (position-dependent clock rate) on a different substrate (γ profile rather than spacetime curvature), without invoking GR machinery. Mass enters in Link 4 via the wave-death residue; Link 5 says the same γ profile that produces the mass also produces the position-dependent tick that, in GR, mass would create gravitationally.

The link is structural: not a derivation of GR, but an in-framework analog of one of its central features.

Source: [On How Gamma Became the Tick](../reflections/ON_HOW_GAMMA_BECAME_THE_TICK.md), [Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md).

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

*Terminology note: "Sacrifice zone" is a legacy name from [Resonant Return](../experiments/RESONANT_RETURN.md) (March 2026). Subsequent work showed nothing is actually sacrificed: the edge qubit acts as an impedance-matched entrance pupil that couples external γ into a cavity mode ([Sacrifice Zone Optics](../experiments/SACRIFICE_ZONE_OPTICS.md), April 4), and the optimal profile is best understood as a controlled symmetry break that splits a degenerate Liouvillian cluster and creates one slow eigenmode ([Sacrifice Geometry](../experiments/SACRIFICE_GEOMETRY.md), April 9-10). The 360× figure is a ratio against a J-blind receiver baseline; under uniform γ₀ a better receiver beats it without any γ-profile engineering ([Receiver vs Gamma Sacrifice](../experiments/RECEIVER_VS_GAMMA_SACRIFICE.md), April 23). The wave-death argument below does not depend on which framing is used: only that some γ-concentration mechanism produces a position-dependent decoherence rate.*

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

1. **Crossing-time gradient.** Compute the site-resolved CΨ = 1/4 crossing time t_cross(n) for a chain with sacrifice-zone γ. Does t_cross increase with distance from the sacrifice edge? If so: the crossing-time landscape IS the gravitational-potential analog. *(Addressable now with the C# brecher tool; not yet run as of 2026-05-23.)*

2. **Feedback convergence.** Start with uniform γ plus a small perturbation. Assume mass ∝ accumulated classical weight. Assume gravity ∝ mass attracts more wave weight. Does the feedback converge to a stable γ profile? What shape? *(Structural ceiling: see "Structural ceilings" below; Lindblad has no in-framework attraction mechanism.)*

3. **1/r test.** Does the crossing-time profile (or the mass accumulation profile) fall like 1/r from the sacrifice edge? In 1D this would be 1/|n - n_edge|. *(Addressable now; same computation as #1.)*

## What Would Refute This

1. **Crossing time is uniform.** If t_cross(n) does not vary with position even when γ varies, the connection fails. *(Addressable now; same computation as Confirm #1.)*

2. **Wave death does not create directional attraction.** If the classical weight accumulation does not preferentially attract further wave weight toward the mass center, there is no feedback and no gravity analog. *(This is not actually a refutation criterion: it is already a structural feature of Lindblad. See "Structural ceilings" below.)*

3. **The feedback diverges.** If mass accumulation → more γ → more mass diverges to a singularity rather than a stable profile, the mechanism produces black holes everywhere, not stable gravitational fields. *(Structural ceiling: requires the feedback mechanism the framework lacks.)*

## Structural ceilings

Three of the items above (Confirm #2, Refute #2, Refute #3) and the "Open gap" of Link 4 step 7 are not pending experiments but structural features of the Lindblad framework. They cannot be tested by any computation within the framework, because the framework lacks the mechanism the test would require.

The specific lack: in Lindblad, L_H (wave propagation, generated by the Hamiltonian) and L_D (wave death, generated by the dissipator) are independent terms; the dissipator does not influence the Hamiltonian, and the parity selection rule \[P_XY, L\] = 0 ([Parity Selection Rule](../docs/proofs/PROOF_PARITY_SELECTION_RULE.md)) algebraically prevents directional attraction within the framework. Without an in-framework mechanism for mass to redirect wave propagation, the feedback loop (mass → gravity → wave attraction → more mass) cannot close. Closing it requires external physics (GR or equivalent), as Link 4 step 7 acknowledges.

[Pair Breaking at the Horizon](PAIR_BREAKING_AT_THE_HORIZON.md) inherits this same gap without closing it; [The Polarity Layer](THE_POLARITY_LAYER.md)'s channel-not-memory reframing sharpens what mass IS but does not provide the missing attraction step. The closest related result is the Hopf bifurcation in [Fragile Bridge](FRAGILE_BRIDGE.md), but that requires negative γ (gain) which decoherence does not provide.

---

## What This Does Not Claim

- That gravity IS this mechanism in our universe. The hypothesis
  applies within the d(d-2)=0 palindromic framework. If the framework
  describes reality (87,376 eigenvalues, zero exceptions), the
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
