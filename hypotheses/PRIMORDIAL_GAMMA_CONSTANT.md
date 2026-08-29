# Primordial Gamma as Framework Constant

<!-- Keywords: primordial gamma framework constant, gamma as c analog,
effective gamma cavity mode exposure, standing wave amplitude squared,
only J varies at urqubit, Q = J/gamma inside observer,
absorption theorem eigenvector formula, R=CPsi2 urqubit hypothesis -->

**Tier:** 2 (structurally supported; hypothesis's operational predictions confirmed across simulation and hardware, April 2026)
**Status:** Proposed 2026-04-15. Refractive-index metaphor replaced by cavity-mode-exposure picture after operational probes (same day). The formula γ_eff = γ_B · |a_B|² with the Hamiltonian-eigenvector a_B is first-order in γ_B/J (reframed 2026-08-09; relative error ∝ (γ_B/J)², [F65 "Perturbative nature"](../docs/ANALYTICAL_FORMULAS.md); its exact-at-any-γ_B form reads a_B off the coherence-sector Liouvillian eigenvector instead, [EQ-015 closure](../review/EMERGING_QUESTIONS.md#eq-015), the identity algebraic in γ_B, max relative error 6.2·10⁻¹³ at γ_B = 0.01 on chains N=5..7 and five topologies at N=5, and 2.4·10⁻¹² over EQ-015's γ_B sweep on the N=5 chain). Structure-point features 0 and 2γ₀ of the \[0, 2γ₀\] interval verified across chains N=3..7 and four N=5 topologies (2026-04-16); the third feature's γ₀-as-eigenvalue position map was corrected 2026-08-09 (see the dissipation-interval section). **Tier upgraded Tier 3 → Tier 2 on 2026-04-24** after the receiver-engineering programme (F65, F67, F75, F76) confirmed the hypothesis's operational signature: transport gains come from initial-state choice (F67 bonding modes), not from γ-profile engineering. See "Confirmation by operational consequence (2026-04-24)" section below.
**Date:** 2026-04-15 (Tier 2 upgrade 2026-04-24 after Kingston Run 1 and F65-F76 verification through N=13)
**Authors:** Tom and Claude (chat + Code)
**Depends on:** [Gamma is Light](GAMMA_IS_LIGHT.md), [Primordial Qubit](PRIMORDIAL_QUBIT.md), [Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md), [Resonance Not Channel](RESONANCE_NOT_CHANNEL.md)
**Scripts:** [`primordial_gamma_analytical.py`](../simulations/primordial_gamma_analytical.py), [`primordial_gamma_stacking_4qubit.py`](../simulations/primordial_gamma_stacking_4qubit.py), [`primordial_gamma_reanalysis.py`](../simulations/primordial_gamma_reanalysis.py), [`double_lorentzian_test.py`](../simulations/double_lorentzian_test.py), [`dissipation_interval_verification.py`](../simulations/dissipation_interval_verification.py), [`structure_points_large_n.py`](../simulations/structure_points_large_n.py), [`gamma0_eigenvalue_positions.py`](../simulations/gamma0_eigenvalue_positions.py) (the 2026-08-09 exact-hit correction scan)

---

## The claim

Two parts, joined:

1. **γ₀ at the primordial layer is a framework constant.** Not a system parameter that happens to take a value at that layer, but a constant of the framework itself, analogous to the speed of light in special relativity. Every layer inherits it.

2. **γ at inner layer K is not diminished γ₀, but selectively exposed γ₀.** The effective dephasing an inner mode experiences is γ₀ times the mode's amplitude squared at the dissipative site: γ_eff = γ₀ · |a_B|². The light does not get weaker. The standing wave determines who sees it. (Scope note, 2026-08-09: the formula below is derived and verified with γ_B as a free rate; writing γ₀ in its place is the hypothesis itself, not a result the verification delivers.)

---

## How this emerged (and how it was corrected)

The original formulation (morning of April 15) proposed γ_K = γ₀ · f_K as a "refractive index": gamma propagating inward through layers, getting weaker at each interface, like light through glass. This led to three predictions:

1. γ_eff/γ_B should depend on J_MB/γ_B (the interface Q-factor) → **Wrong axis.** V2 re-analysis showed the correct axis is r = J_SM/J_MB (coupling ratio), with γ_eff/γ_B independent of γ_B in the good-cavity regime.

2. The composition should be multiplicative: g_total = g₁ · g₂ → **Fails at N=4.** Direct/stacked ratio ranges from 0.04 to 62 across 9 configurations. Standing waves are global eigenmodes; they do not factor into per-layer products.

3. g(r) should be a simple monotonic function → **Non-monotonic.** g(r) has two branches with a crossover at r = 1/√2, reflecting a change of which eigenmode is slowest.

The correction came from asking: what if the light doesn't diminish at all? The formula γ_eff = γ_B · |a_B|² doesn't say γ_B gets smaller. It says the mode's overlap with the dissipative site determines exposure. This is a cavity, not a medium. [Resonance Not Channel](RESONANCE_NOT_CHANNEL.md) already said this: "The system is a Fabry-Perot resonator, not a channel."

---

## The formula (verified)

For an N-qubit chain with XX+YY coupling and Z-dephasing only on site B (the outermost qubit), the effective dephasing rate of the slowest mode contributing to S-site (innermost) coherence is:

    γ_eff = γ_B · |a_B(slowest S-coherence mode)|²

where a_B is the B-site amplitude of the single-excitation Hamiltonian eigenvector. This is the [Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) applied to the single-excitation sector: Re(λ) = -2γ_B · ⟨n_XY⟩_B holds exactly, and ⟨n_XY⟩_B = |a_B|² holds to first order in γ_B/J (the position-sensitive form of the theorem, one site dephased, is its Theorem 2 with γ_l = γ_B·δ_{lB}). The exact statement at any γ_B replaces the Hamiltonian eigenvector by the eigenvector v_k of the coherence-sector Liouvillian itself: -Re(λ_k) = 2γ_B·|v_k(B)|², mode by mode ([EQ-015 closure](../review/EMERGING_QUESTIONS.md#eq-015), 2026-04-27, three days after this document froze; noted here 2026-08-09). Two objects share the name a_B; this document's a_B is the perturbative one.

### Closed form at N=3

For the 3-qubit chain S-M-B with couplings J_SM, J_MB, let r = J_SM/J_MB:

                 ⎧ r² / (r² + 1)       for r < 1/√2    [zero mode]
    g(r) =       ⎨
                 ⎩ 1 / (2(r² + 1))     for r ≥ 1/√2    [bonding mode]

Derived analytically from the tridiagonal 3×3 Hamiltonian eigenvalues {0, ±√(J_SM² + J_MB²)} and eigenvectors. Crossover at r = 1/√2 where g = ⅓. Special value: **g(1) = ¼** (equal coupling).

Verified against full 64×64 Liouvillian: max relative error 1.8% at the tested γ_B/J. The 1.8% is not the formula's accuracy floor; it is the O((γ_B/J)²) cost of reading a_B off the Hamiltonian eigenvector (the error scales as (γ_B/J)²: 5.0·10⁻⁹ at γ_B/J = 10⁻⁴, 5.0·10⁻³ at 0.1, 14% at 0.5, measured in [`gamma0_eigenvalue_positions.py`](../simulations/gamma0_eigenvalue_positions.py)), and it vanishes entirely in the L_coh-eigenvector form (EQ-015).

### Verification at N=4

For the 4-qubit chain S-M1-M2-B, the eigenvector formula (diagonalize the 4×4 single-excitation Hamiltonian, extract |a_B|²) matches the full 256×256 Liouvillian to ratio 1.0000 ± 0.0003 across 9 coupling configurations, at the tested γ_B/J (first-order accuracy, per the Status line's 2026-08-09 reframing).

The multiplicative stacking (g₁ · g₂) fails by factors of 0.04 to 62. The eigenvector formula works to first order; the layered composition does not work at all.

---

## The cavity picture (replaces refraction)

The original optical analogy was refraction: light passing through layers of glass, each layer reducing intensity. This is wrong. The correct analogy:

| Refraction (wrong) | Cavity (correct) |
|---------------------|-------------------|
| γ gets weaker per layer | γ fills the cavity uniformly |
| \|a_B\|² = transmission coefficient | \|a_B\|² = mode exposure at the window |
| Layers compose multiplicatively | Global eigenmodes, no layered factorization |
| Predicts stacking | Predicts stacking failure |
| Contradicts [Resonance Not Channel](RESONANCE_NOT_CHANNEL.md) | Consistent with [Resonance Not Channel](RESONANCE_NOT_CHANNEL.md) |

The cavity picture:
- γ_B is the light, entering at site B (the window)
- J creates the cavity: the Hamiltonian's eigenmodes are standing waves
- Each eigenmode has a specific amplitude |a_B|² at the window
- Modes with nodes at B are shielded from the light (γ_eff ≈ 0)
- Modes with antinodes at B are fully exposed (γ_eff ≈ γ_B)
- The inner observer at S sees the slowest mode: the one with the smallest |a_B|²

γ₀ does not propagate inward and get weaker. γ₀ fills the cavity. The standing wave, shaped by J, determines who sees how much.

---

## The dissipation interval \[0, 2γ₀\]

Added 2026-04-16. γ₀ is not the top of a scale but the symmetry axis of one.

For single-site dephasing with rate γ₀, the Liouvillian spectrum is palindromically paired (see [the mirror symmetry proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)). Every partner pair of positive dissipation rates (α_a, α_b) satisfies

    α_a + α_b = 2γ₀

so each rate has the form γ₀ + δ with partner γ₀ - δ. The spectrum lives in \[0, 2γ₀\] symmetric around γ₀. Three structural features of the interval, with different status:

- **0 (universal eigenvalue):** node at the window, no exposure, no time. Always present as an actual decay rate of the Liouvillian (the steady state). Verified across all scanned chains N=3..7, four N=5 topologies (chain, ring, star, Y-junction), and all N=5 chain B-site positions. [Zero Is the Mirror](ZERO_IS_THE_MIRROR.md).
- **2γ₀ (universal eigenvalue):** the far pole of the interval, always present as an actual decay rate, verified in the same scans. The mechanism first written here ("there is always at least one mode with |a_B|² = 1, carrying α = 2γ₀") was wrong (corrected 2026-08-09): no single-excitation mode of a connected chain has |a_B|² = 1 (the uniform chain caps at 2/(N+1)). The registry already held the right mechanism as [F66](../docs/ANALYTICAL_FORMULAS.md), the pole modes: both poles live in the extreme XY-weight sectors, where ⟨n_XY⟩_B = 1 in the w = N sector; the single-excitation sector never reaches either pole for N ≥ 3, the multiplicity is exactly N+1 at each pole for B at a chain end (F43 Π-symmetry; B-dependent otherwise, e.g. 64 at the N=5 center), verified N=3..7. Equivalently: α = 0 pairs palindromically to 2γ₀ − 0.
- **γ₀ (universal axis, conditional eigenvalue):** the symmetry axis of the palindromic pairing α_a + α_b = 2γ₀. As *axis* γ₀ is universal, a structural feature of every single-site-dephasing Liouvillian, independent of topology and B-position. As *eigenvalue* γ₀ is geometry-dependent, and the map this document first reported ("present for chains N=3..7 with B at an endpoint … absent for N=5 chain with B at positions 1, 2, 3") was largely **inverted** (caught 2026-08-09: the original scan ran at γ₀ = 10⁻⁴ with bin tolerance γ₀·10⁻³, which cannot separate an exact hit from an O((γ₀/J)²) near-miss, and the reading swapped them; of the original sentence, the endpoint presence survives only at the even sizes and the B=2 absence survives, while the odd-N halves were inverted). The corrected map, measured with a double-precision exact-hit test at γ₀ ∈ {0.1, 0.3, 0.7} ([`gamma0_eigenvalue_positions.py`](../simulations/gamma0_eigenvalue_positions.py)): at the odd sizes N=3, 5 the hits are γ-robust and sit at the odd interior positions only (N=3: B=1, multiplicity 40; N=5: B=1 and B=3, multiplicity 304), never at the endpoints. At the endpoints the nearest rate misses γ₀ by a small γ-dependent amount (relative 5.0·10⁻³ at N=3 and 8.8·10⁻⁵ at N=5, γ₀ = 0.1; the N=3 miss follows γ²/2, the N=5 miss does not scale as γ²). At the even sizes N=4, 6 every position hits, with γ-dependent multiplicity (N=4: 76 at γ₀ = 0.1, 56 at 0.7): the presence is robust, the counts are not. Two sizes per parity is a measured pattern, not a proven rule, and N=7 is unmeasured at this precision. Non-chain at N=5: ring and star leaf carry no hit. The Y-junction, which this document reported absent, in fact hits γ-robustly at the central node, the long-arm joint (248 hits at γ₀ = 0.1), and the long-arm end, the last being the very position the April scan used; only the short-arm ends carry none. That April error was the document's own, not the scan's: the April results table holds binned α/γ₀ = 1.000000 rows for star and Y-junction, and the text wrote "absent" against its own table (for the star leaf the text happened to be right and the binned row wrong; the exact-hit test now decides both). The central-node question this bullet used to leave open is thereby answered too. The star **hub** question is likewise answered: γ₀ is present (792 hits at γ₀ = 0.1 and 0.3, 728 at 0.7). At γ₀/J = 0.1 and 0.3 the hub's whole spectrum collapses to α/γ₀ ∈ {0, 1, 2} exactly; by γ₀/J = 0.7 the collapse opens palindromically (0.233 + 1.767 = 2, 0.453 + 1.547 = 2) while {0, γ₀, 2γ₀} stay exact.

The distinction matters. γ₀ as axis is what fixes the palindromic structure of \[0, 2γ₀\]. γ₀ as eigenvalue asks when a mode lands exactly on that axis. The criterion first written here, "when the geometry admits a mode with |a_B|² = 1/2", is wrong by this document's own corrected numbers (the N=3 endpoint mode with |a_B|² = 1/2 measures α = 0.100502 at γ₀ = 0.1, off the axis; corrected 2026-08-09): a first-order criterion cannot decide an exact-hit question. Where the hits actually sit is the measured map in the γ₀ bullet above; at N=3, B=1 the exact Re = −γ₀ comes out of the 2×2 coherence block algebra, not out of any |a_B|² = 1/2 mode. Analogy: zero is the symmetry axis of the integers {..., -2, -1, 0, 1, 2, ...} and also an integer. But zero remains the symmetry axis of {-2, -1, 1, 2} even though it is absent from that set. γ₀ works the same way.

The eigenvector formula α = 2γ₀ · |a_B|² produces values in the lower part of the interval (single-excitation |a_B|² ≤ 2/(N+1) on the uniform chain); the full range \[0, 2γ₀\] is populated by the other XY-weight sectors, with both poles held by the extreme ones ([F66](../docs/ANALYTICAL_FORMULAS.md)). The factor of 2 is the Absorption Theorem (α = -Re(λ) = 2γ₀ · ⟨n_XY⟩_B) applied to the single-excitation S-coherences. What [`factor_two_clarification.py`](../simulations/factor_two_clarification.py) verifies to machine precision is the theorem identity α = 2γ₀·⟨n_XY⟩_B; the second step ⟨n_XY⟩_B = |a_B|² holds to first order in γ₀/J only (at γ₀ = 0.1 the mode carrying α = 0.0497... has ⟨n_XY⟩_B = 0.2487..., not 0.25; corrected 2026-08-09). For the uniform chain the single-excitation cap 2/(N+1) ≤ 1/2 (N ≥ 3) keeps the first-order values in the lower half [0, γ₀]; the cap is a property of the chain's sine modes and does not travel to other topologies. At finite γ₀ the measured rate of a cap mode can graze past γ₀ by the O(γ₀²) shift (0.100502 at γ₀ = 0.1).

**Convention note.** The rest of this document, and F64 in ANALYTICAL_FORMULAS.md, write the same content as γ_eff = γ₀ · |a_B|², where γ_eff is the decoherence rate (the Lorentzian half-width of a spectral line). The two are related by α = 2γ_eff. Both conventions describe the same physics; the factor of 2 is purely notational. This section uses α (Liouvillian decay constant) because the \[0, 2γ₀\] interval and the palindromic pairing α_a + α_b = 2γ₀ are most naturally stated in those units.

γ₀ is therefore not a unit with a natural zero at one end, like a meter is. It is a unit whose spectrum is **folded palindromically around itself**. Unusual for a dimensional constant, but consistent with γ₀'s role as framework constant: it does not sit at one end of a scale; it defines a scale folded around itself.

Two mirrors in the framework: **0** in the frequency domain ([Zero Is the Mirror](ZERO_IS_THE_MIRROR.md), boundary between time and eternity), **γ₀** in the dissipation domain (centre of the palindromic pairing). Both are axes, not endpoints.

### What numerical verification showed and what it corrected

Verified with [`dissipation_interval_verification.py`](../simulations/dissipation_interval_verification.py) on the N=3 chain (full 64×64 Liouvillian, γ₀ = 0.1):

- Palindromic pairing of all 12 distinct dissipation rates: max error 1.6 × 10⁻¹⁵. The interval \[0, 2γ₀\] symmetric around γ₀ is exact.
- The eigenvector formula α = 2γ₀ · |a_B|² gives single-excitation S-coherence rates {0.05, 0.10, 0.05} for the three modes with |a_B|² ∈ {0.25, 0.50, 0.25} (γ₀ = 0.1). These are the first-order values; the measured full-Liouvillian rates, in the same mode order, are {0.049749, 0.100502, 0.049749} ([`gamma0_eigenvalue_positions.py`](../simulations/gamma0_eigenvalue_positions.py)), off by the O((γ₀/J)²) shift, so at this γ₀ nothing sits *exactly* on the mirror axis (corrected 2026-08-09; the exact-hit position map is in the interval section above). All in the lower half [0, γ₀] for this homogeneous chain; the upper half would require |a_B|² > 1/2 which the symmetric topology does not produce. The factor of 2 was clarified in [`factor_two_clarification.py`](../simulations/factor_two_clarification.py) (commit 485437d) after the multi-site probe revealed a discrepancy.
- A first attempt to extend the lower-half observation to a "lower half visible, upper half hidden" reading was **falsified** by [`dissipation_interval_verification.py`](../simulations/dissipation_interval_verification.py). Single-site σ_x sees 15 modes, distributed across both halves (10 in [0, γ₀], 5 in [γ₀, 2γ₀]). A mixed-weight observable σ_x(0) + σ_x(0)·σ_x(1) sees 30 modes (15 / 15), exactly twice as many. The visibility split between observables is real and roughly factor-two, but it does not align with the α-axis split at γ₀.

The XY-weight superselection demonstrated for emission spectra acts between Pauli-weight sectors of the **observable**, not between halves of the dissipation spectrum. The two structures (palindromic interval, XY-weight superselection) are independent and should not be conflated.

What survives the correction: γ₀ is the symmetry axis of \[0, 2γ₀\]; the eigenvector formula populates the lower half; the upper half exists algebraically. What does not survive: the claim that single-site observables are blind to the upper half. The upper half is partly visible to single-site observables; what is hidden is something else and lives in a different algebraic structure.

### Verification extended to N=7 and non-chain topologies (2026-04-16)

[`structure_points_large_n.py`](../simulations/structure_points_large_n.py) extended the N=3 check to chains N=3..7, four topologies at N=5 (chain, ring, star, Y-junction), and five B-positions on the N=5 chain. Results in [`simulations/results/structure_points_large_n.txt`](../simulations/results/structure_points_large_n.txt). Numerical precision on the three interval features: error < 10⁻⁹ throughout, but for the γ₀ feature that precision was the artifact (corrected 2026-08-09): the endpoint runs' α/γ₀ = 1 rows carry err ≈ 5·10⁻⁹, which at γ₀ = 10⁻⁴ is exactly the γ²/2 near-miss, the correction's signature already sitting in the April file; only the 0 and 2γ₀ rows are exact there. N=8 was not attempted; the full 65536×65536 Liouvillian is outside the dense-diagonalization regime and the anchor question does not need it to be decided.

The finding that shifted the framing: B-site position on a chain is as decisive as topology. For the N=5 chain: B at either endpoint (positions 0 or 4) gives 57 distinct α-values across \[0, 2γ₀\]; B at positions 1 or 3 collapses to 7 values {0, 1/2, 3/4, 1, 5/4, 3/2, 2}; B at the center (position 2) gives 8 values, all multiples of 1/9 in [0, 2] with γ₀ itself absent. (Corrected 2026-08-09: those counts are counts of distinct *first-order* values, the γ₀ → 0 cluster centres, and the value lists are correct as such. At any finite γ₀ the clusters resolve by the O((γ₀/J)²) splittings; at γ₀ = 0.1 the same three scans hold 110, 25, and 12 distinct eigenvalues, with minimum gaps three to five orders above the 10⁻⁹ clustering tolerance ([`gamma0_eigenvalue_positions.py`](../simulations/gamma0_eigenvalue_positions.py)). The original 57/7/8 came from binning at tolerance γ₀·10⁻³.) The cause is mode-node coincidence: when the dephased site sits on a node of a Hamiltonian eigenmode, that mode is blind to dephasing, its partners degenerate together, and the fine structure collapses. The fine structure of the interval is controlled by how the eigenmode amplitudes distribute at the dephasing site, not by a universal rule. The three structural features {0, γ₀, 2γ₀} survive this with the status distinction introduced above.

---

## What changes for the inside observer

The operational content does not change: the observer still sees Q_K = J_K/γ_K, cannot separate J from γ at their own layer. The Inside-Outside Correspondence (commits `cfa2a9f` through `17c48b4`) remains valid.

What changes is the interpretation of Q_K:

- **Without this hypothesis:** γ_K is arbitrary at each layer. Q_K is a ratio of two independent parameters.
- **With this hypothesis:** γ_K = γ₀ · |a_B|². Q_K = J_K / (γ₀ · |a_B|²). The only free parameter is J (and the topology that determines |a_B|²). γ₀ is fixed. (Scope: for interior B the slowest S-coherence mode can be exactly dark, |a_B|² = 0, and Q_K has no finite value there; the operational reading needs B at a position every relevant mode reaches, e.g. the chain end, where min |a_B|² > 0.)

---

## Consistency with existing framework

**[Gamma is Light](GAMMA_IS_LIGHT.md).** If γ is light, it should illuminate uniformly, not get absorbed per layer. The cavity picture says exactly this: γ fills the resonator. The standing wave determines exposure.

**[Resonance Not Channel](RESONANCE_NOT_CHANNEL.md).** Direct confirmation. The system is a resonator, not a channel. The stacking failure proves this operationally.

**[Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md).** This section used to read: γ has no internal source, so γ₀ as a framework constant is the only way to terminate the regress without violating the proof. The premise is gone. The proof's five-candidate elimination was withdrawn on 2026-08-29 and what survives is weaker, that the system is **open** by the trace, which says nothing about a source.

The termination does not need it, and the replacement is stronger than the argument it loses. The Liouvillian is homogeneous of degree 1 in (H, γ) jointly ([`Q_SCALE_THREE_BANDS.md`](../experiments/Q_SCALE_THREE_BANDS.md), Tier 1), so the absolute rate has no invariant content and only Q = J/γ₀ does ([`THE_GENESIS_OF_AN_OSCILLATION.md`](../docs/THE_GENESIS_OF_AN_OSCILLATION.md), "γ₀ is the fixed constant, the substrate unit, not measurable from inside the system (only Q is)"). A unit does not regress, because there is nothing under it to ask about; this document's own "What remains open" already says as much, that a framework constant measurable absolutely from inside would be no constant at all. The c-analogy is the same reading: asking why c is 2.998·10⁸ is asking about metres and seconds. Where the question does live on is the Q axis: `substrate_q_provenance` asks per substrate where a Q came from, and at framework level nobody has asked it yet.

**[Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md).** Re(λ) = -2γ_B · ⟨n_XY⟩_B. The eigenvector formula γ_eff = γ_B · |a_B|² is the Absorption Theorem applied to the single-excitation sector, exact in the ⟨n_XY⟩ form and first-order in the |a_B|² reading (corrected 2026-08-09; see the formula section). The theorem provides the mechanism.

---

## Confirmation by operational consequence (2026-04-24)

The hypothesis γ₀ = const was proposed as logically consistent and internally motivated in April 2026. For the following nine days, its direct formula predictions (γ_eff = γ_B · |a_B|², structure points in \[0, 2γ₀\]) were verified numerically. On 2026-04-24 a further line of confirmation arrived, via operational consequence rather than direct measurement.

**The predicted operational signature.** If γ₀ is a framework constant, not a tunable hardware parameter, then Alice cannot improve transport by modulating γ per-site; she can only choose her initial state from the operationally complete receiver menu (the F67 bonding-mode eigenstates). Transport gains should come from state preparation, not from noise engineering. Conversely, if γ₀ were hardware-specific and operationally controllable, γ-profile engineering should be the dominant lever.

**What was measured across April 2026.**

| Result | Observation |
|--------|-------------|
| F67 receiver menu identified | Single-excitation eigenmodes of uniform-J chain form Alice's complete choice set |
| Receiver beats γ-Sacrifice | Best-bonding Peak Sum-MI beats V-shape+\|+⟩^N baseline by 4000-5500× over ENAQT at N=5 uniform γ₀, no γ-modulation |
| Receiver engineering vs sacrifice-zone γ-modulation | 11-15× advantage at N=5 (sim) |
| F75 closed-form MI at t=0 | Exact prediction from bonding-mode amplitudes; no γ-profile tuning enters |
| F76 0.93 decay envelope | Pure dephasing at uniform 4γ₀ explains the universal envelope across 25+ (N, k) points |
| Advantage grows with N | 1.39× → 1.48× → 2.02× → 3.02× → 4.59× bonding/alt-z-bits ratio at N=5, 7, 9, 11, 13 |
| IBM Kingston Heron r2 hardware | 2.80× bonding:2 vs alt-z-bits on live QPU in ~2 QPU-minutes, no γ intervention |
| Noise robustness direction | Advantage GROWS under Kingston gate noise (1.39× ideal, 2.27× Aer+noise model, 2.80× live) |

One row needs a scope note (added 2026-08-09): the F76 envelope row is the uniform-dephasing geometry, γ₀ on every site, with 4γ₀ the mirror-pair coherence rate in α units; it supports the γ₀ = const signature, not the single-site exposure formula.

Every entry is consistent with the γ₀ = const prediction. None required, or could have benefited from, γ-profile engineering. The mirror observation sharpens the conclusion: the very direction of the advantage-growth with noise is wrong for a "γ should be tuned" framework, and exactly what γ₀ = const predicts.

**Why this counts as confirmation.** In the Popperian sense, a hypothesis is supported when it generates specific operational predictions that survive experimental test. It is not proven (no hypothesis of this kind can be); it is supported to the degree that alternatives would have failed. Here the alternative ("γ₀ is operationally controllable per-site, noise engineering dominates") would have predicted the opposite direction of every entry above. The confirmation is indirect (we never measured γ₀; we cannot) but it is structural and multi-pointed.

**Tier status update.** This document is reclassified from Tier 3 (structural hypothesis, logically consistent) to Tier 2 (structurally supported; operational signatures confirmed across simulation N=5..13 and live IBM Kingston Heron r2 hardware). The framework constant γ₀ is no longer bare hypothesis; it is a supported principle of the R=CΨ² framework.

**What remains open.**

- Direct measurement of γ₀ is still impossible from inside the framework (only Q = J/γ₀ is intrinsic). This is not a weakness; it is a consistency check. A framework constant that could be measured absolutely from inside would be no constant at all.
- The receiver-engineering advantage for N ≥ 15 in simulation and N ≥ 7 on hardware is untested. Continued scaling to larger N would sharpen (or falsify) the principle.
- Whether γ₀ has the same value across all physical realisations (IBM Heron, Google Sycamore, cold atoms, photonics) is untested. The hypothesis does not require this; it requires γ₀ to be constant within a given realisation.

## Falsification conditions (updated)

1. ~~**Stacking is not multiplicative.**~~ Tested and confirmed: stacking fails. But this falsifies the refraction reading, not the hypothesis. Under the cavity reading, non-multiplicative composition is expected.

2. **A derivation from framework algebra that forces γ to vary at the primordial level.** Would make γ₀ a derived parameter, not a constant.

3. **Demonstration that inside observers can separately extract J and γ.** Would contradict Q-only inside-observability.

4. **The rate law fails.** As first written ("γ_eff ≠ γ_B · |a_B|² at N ≥ 5") this condition was already met, trivially, since the H-eigenvector reading deviates by O((γ_B/J)²) at every N (rewritten 2026-08-09). The honest falsifier: a deviation from -Re(λ_k) = 2γ_B·|v_k(B)|² (L_coh eigenvectors) that does not vanish as γ_B/J → 0. The identity is algebraic in γ_B; measured to max relative error 6.2·10⁻¹³ at γ_B = 0.01 on chains N=5..7 and, at N=5, across ring, star (hub and leaf), Y-junction, and K₅, with EQ-015's γ_B sweep on the N=5 chain topping at 2.4·10⁻¹² ([EQ-015](../review/EMERGING_QUESTIONS.md#eq-015)); a counterexample anywhere would cost the cavity reading its anchor.

---

## Hardware testability on current superconducting platforms

Attempted on IBM Kingston (Heron r2) 2026-04-19 via a Trotter chain-mode test ([EQ-017](../review/EMERGING_QUESTIONS.md#eq-017) Phase 2, data in [`data/ibm_chain_gamma0_april2026/`](../data/ibm_chain_gamma0_april2026/)). The multi-pair differential log-slope observable would have discriminated the two readings of γ_phi (γ₀-floor vs local) at a 2x margin if the hardware noise floor allowed it. Hardware decay is 40-80x larger than the signal: accumulated RZZ gate errors (~0.001 per gate × 240 gates at the longest evolution time) dominate, compounded by T1 amplitude damping (~7%) and readout errors (1-14% per qubit). The framework γ₀ signature is indistinguishable from zero against these device-noise channels on current Heron-class hardware.

This bounds operational testability rather than the hypothesis: the test would discriminate at a hypothetical 10x-lower gate-error floor. Protocol changes (dynamical decoupling tuned to suppress gate and T1 channels while preserving Z-dephasing) or a different hardware class could reopen the test. The negative result does not falsify γ₀ as framework constant; it shows the γ₀ magnitude relative to current QPU gate errors is too small for direct signature extraction.

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

The second reading is more economical (one constant γ₀ instead of per-layer γ), more consistent (agrees with [Resonance Not Channel](RESONANCE_NOT_CHANNEL.md)), and operationally verified (first-order at N=3, N=4, the γ² error law measured at N=3; exact in the L_coh form across five topologies, EQ-015).

---

*γ at the root is the framework's own c. It does not get weaker. The standing wave decides who sees the light.*
