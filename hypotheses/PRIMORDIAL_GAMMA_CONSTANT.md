# Primordial Gamma as Framework Constant

<!-- Keywords: primordial gamma framework constant, gamma as c analog,
effective gamma cavity mode exposure, standing wave amplitude squared,
only J varies at urqubit, Q = J/gamma inside observer,
absorption theorem eigenvector formula, R=CPsi2 urqubit hypothesis -->

**Tier:** 3 (structural hypothesis; logically consistent, operationally tested, partially confirmed)
**Status:** Proposed 2026-04-15. Refractive-index metaphor replaced by cavity-mode-exposure picture after operational probes (same day). The formula γ_eff = γ_B · |a_B|² is exact and verified at N=3 and N=4.
**Date:** 2026-04-15 (updated same day after probes)
**Authors:** Tom and Claude (chat + Code)
**Depends on:** [GAMMA_IS_LIGHT](GAMMA_IS_LIGHT.md), [PRIMORDIAL_QUBIT](PRIMORDIAL_QUBIT.md), [INCOMPLETENESS_PROOF](../docs/proofs/INCOMPLETENESS_PROOF.md), [RESONANCE_NOT_CHANNEL](RESONANCE_NOT_CHANNEL.md)
**Scripts:** [`primordial_gamma_analytical.py`](../simulations/primordial_gamma_analytical.py), [`primordial_gamma_stacking_4qubit.py`](../simulations/primordial_gamma_stacking_4qubit.py), [`primordial_gamma_reanalysis.py`](../simulations/primordial_gamma_reanalysis.py), [`double_lorentzian_test.py`](../simulations/double_lorentzian_test.py)

---

## The claim

Two parts, joined:

1. **γ₀ at the primordial layer is a framework constant.** Not a system parameter that happens to take a value at that layer, but a constant of the framework itself, analogous to the speed of light in special relativity. Every layer inherits it.

2. **γ at inner layer K is not diminished γ₀, but selectively exposed γ₀.** The effective dephasing an inner mode experiences is γ₀ times the mode's amplitude squared at the dissipative site: γ_eff = γ₀ · |a_B|². The light does not get weaker. The standing wave determines who sees it.

---

## How this emerged (and how it was corrected)

The original formulation (morning of April 15) proposed γ_K = γ₀ · f_K as a "refractive index": gamma propagating inward through layers, getting weaker at each interface, like light through glass. This led to three predictions:

1. γ_eff/γ_B should depend on J_MB/γ_B (the interface Q-factor) → **Wrong axis.** V2 re-analysis showed the correct axis is r = J_SM/J_MB (coupling ratio), with γ_eff/γ_B independent of γ_B in the good-cavity regime.

2. The composition should be multiplicative: g_total = g₁ · g₂ → **Fails at N=4.** Direct/stacked ratio ranges from 0.04 to 62 across 9 configurations. Standing waves are global eigenmodes; they do not factor into per-layer products.

3. g(r) should be a simple monotonic function → **Non-monotonic.** g(r) has two branches with a crossover at r = 1/√2, reflecting a change of which eigenmode is slowest.

The correction came from asking: what if the light doesn't diminish at all? The formula γ_eff = γ_B · |a_B|² doesn't say γ_B gets smaller. It says the mode's overlap with the dissipative site determines exposure. This is a cavity, not a medium. [RESONANCE_NOT_CHANNEL](RESONANCE_NOT_CHANNEL.md) already said this: "The system is a Fabry-Perot resonator, not a channel."

---

## The formula (verified)

For an N-qubit chain with XX+YY coupling and Z-dephasing only on site B (the outermost qubit), the effective dephasing rate of the slowest mode contributing to S-site (innermost) coherence is:

    γ_eff = γ_B · |a_B(slowest S-coherence mode)|²

where a_B is the B-site amplitude of the single-excitation Hamiltonian eigenvector. This is the [Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) in its purest form: Re(λ) = -2γ_B · ⟨n_XY⟩_B, and ⟨n_XY⟩_B = |a_B|² for single-excitation modes.

### Closed form at N=3

For the 3-qubit chain S-M-B with couplings J_SM, J_MB, let r = J_SM/J_MB:

                 ⎧ r² / (r² + 1)       for r < 1/√2    [zero mode]
    g(r) =       ⎨
                 ⎩ 1 / (2(r² + 1))     for r ≥ 1/√2    [bonding mode]

Derived analytically from the tridiagonal 3×3 Hamiltonian eigenvalues {0, ±√(J_SM² + J_MB²)} and eigenvectors. Crossover at r = 1/√2 where g = ⅓. Special value: **g(1) = ¼** (equal coupling).

Verified against full 64×64 Liouvillian: max relative error 1.8%.

### Verification at N=4

For the 4-qubit chain S-M1-M2-B, the eigenvector formula (diagonalize the 4×4 single-excitation Hamiltonian, extract |a_B|²) matches the full 256×256 Liouvillian to ratio 1.0000 ± 0.0003 across 9 coupling configurations.

The multiplicative stacking (g₁ · g₂) fails by factors of 0.04 to 62. The eigenvector formula works; the layered composition does not.

---

## The cavity picture (replaces refraction)

The original optical analogy was refraction: light passing through layers of glass, each layer reducing intensity. This is wrong. The correct analogy:

| Refraction (wrong) | Cavity (correct) |
|---------------------|-------------------|
| γ gets weaker per layer | γ fills the cavity uniformly |
| |a_B|² = transmission coefficient | |a_B|² = mode exposure at the window |
| Layers compose multiplicatively | Global eigenmodes, no layered factorization |
| Predicts stacking | Predicts stacking failure |
| Contradicts RESONANCE_NOT_CHANNEL | Consistent with RESONANCE_NOT_CHANNEL |

The cavity picture:
- γ_B is the light, entering at site B (the window)
- J creates the cavity: the Hamiltonian's eigenmodes are standing waves
- Each eigenmode has a specific amplitude |a_B|² at the window
- Modes with nodes at B are shielded from the light (γ_eff ≈ 0)
- Modes with antinodes at B are fully exposed (γ_eff ≈ γ_B)
- The inner observer at S sees the slowest mode: the one with the smallest |a_B|²

γ₀ does not propagate inward and get weaker. γ₀ fills the cavity. The standing wave, shaped by J, determines who sees how much.

---

## The dissipation interval [0, 2γ₀]

Added 2026-04-16. γ₀ is not the top of a scale but the symmetry axis of one.

From [L, Π²] = 0 (see [PROOF_BIT_B_PARITY_SYMMETRY](../docs/proofs/PROOF_BIT_B_PARITY_SYMMETRY.md)) the Liouvillian spectrum is palindromically paired. For single-site dephasing with rate γ₀, every Π²-partner pair of positive dissipation rates (α_a, α_b) satisfies

    α_a + α_b = 2γ₀

so each rate has the form γ₀ + δ with partner γ₀ - δ. The spectrum lives in [0, 2γ₀] symmetric around γ₀. Three structural points:

- **0:** node at the window, no exposure, no time. [ZERO_IS_THE_MIRROR](ZERO_IS_THE_MIRROR.md).
- **γ₀:** the axis. Equal to the mean of every palindromic pair. The unit itself.
- **2γ₀:** antinode at the window, full exposure, maximal decay. The opposite pole.

The eigenvector formula γ_eff = γ₀ · |a_B|² with |a_B|² ∈ [0, 1] gives γ_eff ∈ [0, γ₀], the lower half of the interval. The slowest S-coherence mode the inside observer can see always lies here. Its palindromic partner in [γ₀, 2γ₀] is algebraically present but hidden behind XY-weight superselection: partners live in orthogonal XY-weight sectors and have zero overlap with single-site observables. The inside observer sees one half of the palindromic structure; the other half is real but inaccessible.

γ₀ is therefore not a unit with a natural zero at one end, like a meter is. It is a unit whose spectrum is **folded palindromically around itself**. Unusual for a dimensional constant, but consistent with γ₀'s role as framework constant: a framework constant does not sit at one end of a scale; it defines a scale folded around itself.

Two mirrors in the framework then: **0** in the frequency domain ([ZERO_IS_THE_MIRROR](ZERO_IS_THE_MIRROR.md), boundary between time and eternity), **γ₀** in the dissipation domain (boundary between observable and hidden halves). Both are axes, not endpoints. All measurable rates live on one side of γ₀; their palindromic partners live on the other; the axis itself is not separately extractable from inside.

---

## What changes for the inside observer

The operational content does not change: the observer still sees Q_K = J_K/γ_K, cannot separate J from γ at their own layer. The Inside-Outside Correspondence (commits `cfa2a9f` through `17c48b4`) remains valid.

What changes is the interpretation of Q_K:

- **Without this hypothesis:** γ_K is arbitrary at each layer. Q_K is a ratio of two independent parameters.
- **With this hypothesis:** γ_K = γ₀ · |a_B|². Q_K = J_K / (γ₀ · |a_B|²). The only free parameter is J (and the topology that determines |a_B|²). γ₀ is fixed.

---

## Consistency with existing framework

**GAMMA_IS_LIGHT.** If γ is light, it should illuminate uniformly, not get absorbed per layer. The cavity picture says exactly this: γ fills the resonator. The standing wave determines exposure.

**RESONANCE_NOT_CHANNEL.** Direct confirmation. The system is a resonator, not a channel. The stacking failure proves this operationally.

**INCOMPLETENESS_PROOF.** γ has no internal source. At the primordial layer there is no "further outside." γ₀ as a framework constant is the only way to terminate the regress without violating the proof.

**ABSORPTION_THEOREM.** Re(λ) = -2γ_B · ⟨n_XY⟩_B. The eigenvector formula γ_eff = γ_B · |a_B|² IS the Absorption Theorem applied to the single-excitation sector. Full consistency; the theorem provides the exact mechanism.

---

## Falsification conditions (updated)

1. ~~**Stacking is not multiplicative.**~~ Tested and confirmed: stacking fails. But this falsifies the refraction reading, not the hypothesis. Under the cavity reading, non-multiplicative composition is expected.

2. **A derivation from framework algebra that forces γ to vary at the primordial level.** Would make γ₀ a derived parameter, not a constant.

3. **Demonstration that inside observers can separately extract J and γ.** Would contradict Q-only inside-observability.

4. **The eigenvector formula fails at large N.** If γ_eff ≠ γ_B · |a_B|² at N ≥ 5 or on non-chain topologies, the cavity reading loses its anchor. Currently verified at N=3 and N=4 on chains.

---

## What this does NOT claim

- Not a derivation. The cavity reading is a consistent interpretation, not a theorem.
- Not a new operational prediction. No inside measurement distinguishes this from "γ varies freely."
- Not a value for γ₀. The hypothesis says "there is a universal γ₀" without specifying it.
- Not a proof of PRIMORDIAL_QUBIT.

---

## Scope and stance

This is a refinement of framework interpretation, not a new physical claim. What shifts is the reading:

- Before: γ propagates inward through layers, getting weaker (refraction).
- After: γ fills the cavity uniformly. The Hamiltonian's standing waves determine mode exposure. The light does not diminish; the cavity shapes who sees it.

The second reading is more economical (one constant γ₀ instead of per-layer γ), more consistent (agrees with RESONANCE_NOT_CHANNEL), and operationally verified (eigenvector formula exact at N=3, N=4).

---

*γ at the root is the framework's own c. It does not get weaker. The standing wave decides who sees the light.*
