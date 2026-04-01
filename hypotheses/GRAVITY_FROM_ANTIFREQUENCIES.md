# Gravity from Anti-Frequencies

**Status:** Hypothesis (Tier 5), grounded in Tier 1-2 palindrome results
**Date:** April 1, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Depends on:**
- [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md) (ω ↔ -ω pairing)
- [Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md) (γ = experienced time)
- [Time Irreversibility Exclusion](../docs/proofs/TIME_IRREVERSIBILITY_EXCLUSION.md) (cross term)
- [Pi as Time Reversal](../experiments/PI_AS_TIME_REVERSAL.md) (standing waves)
- [Zero Is the Mirror](ZERO_IS_THE_MIRROR.md) (Σγ = 0 as origin)

---

## Abstract

The palindromic spectrum pairs every oscillation frequency ω with its
negation -ω. These pairs form standing waves. In a spatially extended
system, the standing-wave pattern creates a position-dependent
effective dephasing profile: a γ landscape. Since γ is the source of
experienced time (Incompleteness Proof), a spatial variation of γ is
indistinguishable from gravitational time dilation. The hypothesis:
gravity is not an external force acting on the palindromic system.
Gravity IS the spatial γ profile induced by the anti-frequency
standing waves. The anti-frequencies do not describe gravity. They
generate it.

---

## The Computed Chain

Each link is Tier 1-2. The hypothesis is the reading of the chain.

### Link 1: Every frequency has an anti-frequency

The palindrome equation Π·L·Π⁻¹ = -L - 2Σγ·I implies: if λ = -d + iω
is an eigenvalue, then -(λ + 2Σγ) = (d - 2Σγ) - iω is also an
eigenvalue. The frequency mirrors: ω ↔ -ω. Every oscillation has a
counter-oscillation. This is proven for all N, all topologies, all
XXZ couplings under Z-dephasing.

### Link 2: The pairs form standing waves

In the rescaled frame (removing the uniform decay envelope e^{-Σγt}),
the palindromic pair at ±ω produces:

    e^{+iωt} + e^{-iωt} = 2cos(ωt)

A standing wave. Not a traveling wave. The mode does not propagate; it
oscillates in place. The node structure of cos(ωt) is fixed in space
(for a chain: determined by the dispersion relation ω_k = 4J sin²(πk/2N)).

Source: [Pi as Time Reversal](../experiments/PI_AS_TIME_REVERSAL.md)

### Link 3: γ is experienced time

The Incompleteness Proof establishes: without γ, no decay, no
irreversibility, no before and after. γ is not a parameter that
correlates with time. γ is the necessary and sufficient condition for
experienced time. Remove γ and t loses its meaning.

The crossing time formula: t_cross × γ = const. Time multiplied by
the rate of time gives a pure number.

Source: [Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md),
Corollary 2

### Link 4: Gravity is time dilation

In general relativity, a gravitational field is equivalent to a
position-dependent clock rate. A clock at lower gravitational potential
runs slower. The Schwarzschild metric:

    dτ² = (1 - 2GM/rc²) dt²

The factor (1 - 2GM/rc²) is the local time dilation. Clocks at
different r tick at different rates. This is not an analogy. This is
the definition of gravity in GR: the geometry of spacetime IS the
variation of clock rates.

---

## The Hypothesis

If γ is experienced time (Link 3), then a spatial variation of γ is
a spatial variation of clock rates. If gravity is a spatial variation
of clock rates (Link 4), then:

**A spatial γ profile IS gravity.**

Not analogous to gravity. Not a model of gravity. Identical to it.
A region with higher γ has faster time (faster decoherence, faster
crossing of CΨ = 1/4). A region with lower γ has slower time. The
difference between them is what we call a gravitational field.

The standing waves from the anti-frequency pairs (Links 1-2) create
a spatial modulation of the effective dynamics. In a chain of N qubits,
the mode at wave number k has spatial profile cos(πkn/N) at site n.
The superposition of all palindromic standing waves creates a spatially
varying envelope.

The hypothesis: this envelope, summed over all palindromic pairs, acts
as an effective γ profile. Regions where constructive interference
concentrates the standing-wave amplitude have higher effective γ
(faster time). Regions of destructive interference have lower effective
γ (slower time). The γ landscape IS the gravitational potential.

---

## What Would Confirm This

1. **Standing-wave envelope = γ profile.** Compute the spatial
   amplitude profile A(n) = Σ_k |c_k|² cos²(πkn/N) for a chain of
   N qubits with specific initial state. Show that A(n) creates an
   effective position-dependent decoherence rate. This is computable.

2. **The profile matches 1/r.** For a 1D chain, the gravitational
   analog would be a monotonically decreasing profile from a "mass"
   at one end. If the standing-wave envelope produces a 1/r-like
   profile for appropriate initial conditions, the connection becomes
   quantitative.

3. **γ dilation matches Schwarzschild.** The time dilation factor
   dτ/dt = √(1 - 2GM/rc²) should correspond to γ(r)/γ_∞ for some
   mapping of chain position to radial coordinate.

## What Would Refute This

1. **The standing-wave envelope is flat.** If all palindromic pairs
   contribute equally at all sites, there is no spatial γ profile and
   no gravitational analog.

2. **The profile does not match any gravitational potential.** If A(n)
   is oscillatory, random, or incompatible with 1/r or any known
   metric, the hypothesis fails.

3. **γ variation does not produce the correct equations of motion.**
   Gravity is not just time dilation; it is geodesic motion in curved
   spacetime. A γ profile that produces time dilation but not geodesic
   deviation would be incomplete.

---

## First Computation: A(n) is Flat (April 1, 2026)

N=3, 4, 5 Heisenberg chain, γ = 0.05 uniform. Computed:

1. All palindromic eigenvalue pairs (ω_k, -ω_k)
2. The right eigenvectors v_k(n) at each site n
3. The standing-wave amplitude profile:
   A(n) = Σ_k (|v_k(n)|² + |v_{-k}(n)|²) for all oscillating pairs

**Result: A(n) is exactly flat.**

| N | Oscillating pairs | Variation (max-min)/mean | Edge/Center |
|---|-------------------|------------------------|-------------|
| 3 | 20 | 0.0001 | 1.0001 |
| 4 | 105 | 0.0000 | 1.0000 |
| 5 | 464 | 0.0000 | 1.0000 |

Each individual palindromic pair has spatial structure (nodes and
antinodes). But when summed over ALL pairs, the structure cancels.
The total standing-wave envelope is uniform. No γ landscape. No
spatial variation. No gravitational profile.

**The hypothesis in its direct form does not work for uniform
systems.** The anti-frequencies create standing waves, but the sum
of all standing waves is flat.

### Why it is flat

The uniform Heisenberg chain with uniform γ has approximate
translational symmetry (exact for ring topology, approximate for
open chain with small edge effects). Summing over all eigenmodes
of a translationally symmetric system averages out the spatial
structure. This is the completeness relation: Σ_k |v_k(n)|² is
constant for a complete orthonormal set.

### What could survive

The direct sum over all modes is flat. But:

1. **Non-uniform γ breaks the symmetry.** The sacrifice-zone result
   shows that concentrating γ on one edge creates dramatic effects
   (139-360x improvement). A pre-existing γ gradient would break
   translational symmetry and create a non-trivial A(n). The question
   becomes: does the palindromic structure amplify an initial γ
   perturbation into a self-consistent profile? This is a bootstrap
   question.

2. **State-dependent weighting.** The sum over ALL modes weights them
   equally. A specific initial state excites specific modes. The
   standing-wave profile for a Bell state at one edge would not be
   flat: it would be peaked where the state lives. The "gravitational
   field of a localized excitation" might have structure even when the
   total does not.

3. **The sacrifice zone itself.** The result that edge sacrifice
   outperforms uniform dephasing by 100x+ means the palindromic
   structure IS sensitive to spatial γ profiles, even though the
   summed standing waves are flat. The sensitivity is in the
   eigenvalue response (rates change with γ profile), not in the
   mode amplitudes (which sum to flat).

---

## The Deeper Reading

The Incompleteness Proof says: γ comes from outside. There is no
internal source.

This hypothesis says: γ does not come from outside uniformly. It comes
from the standing-wave structure of the palindrome itself. The
anti-frequencies create a spatial pattern. The pattern is the
gravitational field. Gravity is not imposed on the system from outside.
Gravity emerges from the palindromic spectrum when the anti-frequencies
interfere spatially.

If this is correct, then the "outside" of the Incompleteness Proof is
not a separate entity. It is the global pattern of the palindrome,
viewed from the local perspective of a single site. Each site
experiences γ differently because the standing waves concentrate
differently at different positions. The variation is gravity. The
source is the palindrome. The palindrome is the system itself.

Nothing comes from outside. Everything comes from the pattern.

---

*The anti-frequencies are not a mirror image of the frequencies.*
*They are the other leg of the standing wave. And the standing wave*
*is not flat. It has nodes and antinodes. The antinodes are where*
*time runs fast. The nodes are where time runs slow. The difference*
*between them is what holds us to the ground.*

*Thomas Wicht, April 1, 2026*
