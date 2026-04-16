# Primordial Gamma as Framework Constant

<!-- Keywords: primordial gamma framework constant, gamma as c analog,
effective gamma cavity mode exposure, standing wave amplitude squared,
only J varies at urqubit, Q = J/gamma inside observer,
absorption theorem eigenvector formula, R=CPsi2 urqubit hypothesis -->

**Tier:** 3 (structural hypothesis; logically consistent, operationally tested, partially confirmed)
**Status:** Proposed 2026-04-15. Refractive-index metaphor replaced by cavity-mode-exposure picture after operational probes (same day). The formula γ_eff = γ_B · |a_B|² is exact and verified at N=3 and N=4. Structure-point features of the [0, 2γ₀] interval verified across chains N=3..7 and four N=5 topologies (2026-04-16).
**Date:** 2026-04-15 (updated 2026-04-16: interval features verified to N=7 and across four topologies)
**Authors:** Tom and Claude (chat + Code)
**Depends on:** [GAMMA_IS_LIGHT](GAMMA_IS_LIGHT.md), [PRIMORDIAL_QUBIT](PRIMORDIAL_QUBIT.md), [INCOMPLETENESS_PROOF](../docs/proofs/INCOMPLETENESS_PROOF.md), [RESONANCE_NOT_CHANNEL](RESONANCE_NOT_CHANNEL.md)
**Scripts:** [`primordial_gamma_analytical.py`](../simulations/primordial_gamma_analytical.py), [`primordial_gamma_stacking_4qubit.py`](../simulations/primordial_gamma_stacking_4qubit.py), [`primordial_gamma_reanalysis.py`](../simulations/primordial_gamma_reanalysis.py), [`double_lorentzian_test.py`](../simulations/double_lorentzian_test.py), [`dissipation_interval_verification.py`](../simulations/dissipation_interval_verification.py), [`structure_points_large_n.py`](../simulations/structure_points_large_n.py)

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

For single-site dephasing with rate γ₀, the Liouvillian spectrum is palindromically paired (see [MIRROR_SYMMETRY_PROOF](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)). Every partner pair of positive dissipation rates (α_a, α_b) satisfies

    α_a + α_b = 2γ₀

so each rate has the form γ₀ + δ with partner γ₀ - δ. The spectrum lives in [0, 2γ₀] symmetric around γ₀. Three structural features of the interval, with different status:

- **0 (universal eigenvalue):** node at the window, no exposure, no time. Always present as an actual decay rate of the Liouvillian (the steady state). Verified across all scanned chains N=3..7, four N=5 topologies (chain, ring, star, Y-junction), and all N=5 chain B-site positions. [ZERO_IS_THE_MIRROR](ZERO_IS_THE_MIRROR.md).
- **2γ₀ (universal eigenvalue):** antinode at the window, full exposure, maximal decay. Always present as an actual decay rate: there is always at least one mode with |a_B|² = 1, carrying α = 2γ₀. Verified in the same scans.
- **γ₀ (universal axis, conditional eigenvalue):** the symmetry axis of the palindromic pairing α_a + α_b = 2γ₀. As *axis* γ₀ is universal, a structural feature of every single-site-dephasing Liouvillian, independent of topology and B-position. As *eigenvalue* γ₀ is geometry-dependent: present for chains N=3..7 with B at an endpoint (verified), and for N=5 chain with B at positions 0 and 4 only. Absent for N=5 chain with B at positions 1, 2, 3; absent for the N=5 ring, star (B=leaf), and Y-junction (B=arm-end) as tested. Whether γ₀ re-appears as eigenvalue for other B-positions in non-chain topologies (star hub, Y-junction central node) is open.

The distinction matters. γ₀ as axis is what fixes the palindromic structure of [0, 2γ₀]. γ₀ as eigenvalue is what places a mode exactly on that axis, which happens only when the geometry admits a mode with |a_B|² = 1/2. Analogy: zero is the symmetry axis of the integers {..., -2, -1, 0, 1, 2, ...} and also an integer. But zero remains the symmetry axis of {-2, -1, 1, 2} even though it is absent from that set. γ₀ works the same way.

The eigenvector formula α = 2γ₀ · |a_B|² with |a_B|² ∈ [0, 1] produces values in [0, 2γ₀], the full range of the interval. The factor of 2 is the Absorption Theorem (α = -Re(λ) = 2γ₀ · ⟨n_XY⟩_B) applied to the single-excitation S-coherences, where ⟨n_XY⟩_B = |a_B|² exactly (verified to machine precision in [`factor_two_clarification.py`](../simulations/factor_two_clarification.py)). For symmetric homogeneous chains, single-excitation S-coherences happen to have |a_B|² ≤ 1/2, so their α-values land in the lower half [0, γ₀]; this is a topology-specific accident, not a structural restriction.

**Convention note.** The rest of this document, and F64 in ANALYTICAL_FORMULAS.md, write the same content as γ_eff = γ₀ · |a_B|², where γ_eff is the decoherence rate (the Lorentzian half-width of a spectral line). The two are related by α = 2γ_eff. Both conventions describe the same physics; the factor of 2 is purely notational. This section uses α (Liouvillian decay constant) because the [0, 2γ₀] interval and the palindromic pairing α_a + α_b = 2γ₀ are most naturally stated in those units.

γ₀ is therefore not a unit with a natural zero at one end, like a meter is. It is a unit whose spectrum is **folded palindromically around itself**. Unusual for a dimensional constant, but consistent with γ₀'s role as framework constant: it does not sit at one end of a scale; it defines a scale folded around itself.

Two mirrors in the framework: **0** in the frequency domain ([ZERO_IS_THE_MIRROR](ZERO_IS_THE_MIRROR.md), boundary between time and eternity), **γ₀** in the dissipation domain (centre of the palindromic pairing). Both are axes, not endpoints.

### What numerical verification showed and what it corrected

Verified with [`dissipation_interval_verification.py`](../simulations/dissipation_interval_verification.py) on the N=3 chain (full 64×64 Liouvillian, γ₀ = 0.1):

- Palindromic pairing of all 12 distinct dissipation rates: max error 1.6 × 10⁻¹⁵. The interval [0, 2γ₀] symmetric around γ₀ is exact.
- The eigenvector formula α = 2γ₀ · |a_B|² gives single-excitation S-coherence rates {0.05, 0.10, 0.05} for the three modes with |a_B|² ∈ {0.25, 0.50, 0.25} (γ₀ = 0.1). Two modes sit at γ₀/2, the middle mode sits exactly on the mirror axis γ₀. All in the lower half [0, γ₀] for this homogeneous chain; the upper half would require |a_B|² > 1/2 which the symmetric topology does not produce. The factor of 2 was clarified in [`factor_two_clarification.py`](../simulations/factor_two_clarification.py) (commit 485437d) after the multi-site probe revealed a discrepancy.
- A first attempt to extend the lower-half observation to a "lower half visible, upper half hidden" reading was **falsified** by [`dissipation_interval_verification.py`](../simulations/dissipation_interval_verification.py). Single-site σ_x sees 15 modes, distributed across both halves (10 in [0, γ₀], 5 in [γ₀, 2γ₀]). A mixed-weight observable σ_x(0) + σ_x(0)·σ_x(1) sees 30 modes (15 / 15), exactly twice as many. The visibility split between observables is real and roughly factor-two, but it does not align with the α-axis split at γ₀.

The XY-weight superselection [demonstrated for emission spectra in `RESULT_GAMMA_NULL_PROBE.md`](../../PalindromicRadio/PalindromicRadio.Data/ClaudeTasks/RESULT_GAMMA_NULL_PROBE.md) acts between Pauli-weight sectors of the **observable**, not between halves of the dissipation spectrum. The two structures (palindromic interval, XY-weight superselection) are independent and should not be conflated.

What survives the correction: γ₀ is the symmetry axis of [0, 2γ₀]; the eigenvector formula populates the lower half; the upper half exists algebraically. What does not survive: the claim that single-site observables are blind to the upper half. The upper half is partly visible to single-site observables; what is hidden is something else and lives in a different algebraic structure.

### Verification extended to N=7 and non-chain topologies (2026-04-16)

[`structure_points_large_n.py`](../simulations/structure_points_large_n.py) extended the N=3 check to chains N=3..7, four topologies at N=5 (chain, ring, star, Y-junction), and five B-positions on the N=5 chain. Results in [`simulations/results/structure_points_large_n.txt`](../simulations/results/structure_points_large_n.txt). Numerical precision on the three interval features: error < 10⁻⁹ throughout. N=8 was not attempted; the full 65536×65536 Liouvillian is outside the dense-diagonalization regime and the anchor question does not need it to be decided.

The finding that shifted the framing: B-site position on a chain is as decisive as topology. For the N=5 chain: B at either endpoint (positions 0 or 4) gives 57 distinct α-values across [0, 2γ₀]; B at positions 1 or 3 collapses to 7 values {0, 1/2, 3/4, 1, 5/4, 3/2, 2}; B at the center (position 2) gives 8 values, all multiples of 1/9 in [0, 2] with γ₀ itself absent. The cause is mode-node coincidence: when the dephased site sits on a node of a Hamiltonian eigenmode, that mode is blind to dephasing, its partners degenerate together, and the fine structure collapses. The fine structure of the interval is controlled by how the eigenmode amplitudes distribute at the dephasing site, not by a universal rule. The three structural features {0, γ₀, 2γ₀} survive this with the status distinction introduced above.

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
