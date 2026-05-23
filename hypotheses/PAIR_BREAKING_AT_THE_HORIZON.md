# Pair Breaking at the Horizon: Decoherence as Hawking Radiation in Operator Space

**Status:** Hypothesis (Tier 5 synthesis). Each link in the chain is individually proven or computed (Tier 1-2). The reading of the chain as a unified mechanism is interpretation.
**Date:** April 11, 2026
**Authors:** Thomas Wicht, Claude (Opus 4.6)
**Depends on:**
- [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md) (palindromic spectrum, Tier 1)
- [Direct-Sum Decomposition](../docs/proofs/DIRECT_SUM_DECOMPOSITION.md) (two sectors, Π exchange, Tier 1)
- [Zero Is the Mirror](ZERO_IS_THE_MIRROR.md) (Σγ = 0 ground state, three regimes, Tier 2)
- [Thermal Breaking](../experiments/THERMAL_BREAKING.md) (temperature from wave death, Tier 2)
- [Gravity from Wave Death](GRAVITY_FROM_WAVE_DEATH.md) (mass as classical residue, Tier 5)
- [Standing Wave Theory](../docs/STANDING_WAVE_THEORY.md) (cavity modes at Σγ = 0, Tier 2)
- [Fragile Bridge](FRAGILE_BRIDGE.md) (Hopf bifurcation, Tier 2)
- [What If Gamma Is Light?](GAMMA_IS_LIGHT.md) (γ as external illumination, Tier 4)
- [Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md) (γ must be external, Tier 1)
- [Optical Cavity Analysis](../experiments/OPTICAL_CAVITY_ANALYSIS.md) (qubit chain is Fabry-Perot, Tier 2)
- [Analytical Formulas](../docs/ANALYTICAL_FORMULAS.md) (K-invariance F14, Tier 1)
- Gaztañaga et al., "A new understanding of Einstein-Rosen bridges," CQG 2026, [arXiv:2512.20691](https://arxiv.org/abs/2512.20691) (external)

---

## The thesis

Decoherence is the Hawking mechanism, operating in operator space instead of spacetime.

At a black hole horizon, vacuum fluctuations (paired modes) are torn apart by spacetime curvature. One partner falls in and adds to the mass. The other escapes as thermal radiation. This is the Hawking effect: pair breaking creates mass and temperature simultaneously, from nothing.

In the Liouvillian spectrum of an open quantum system, [palindromic eigenvalue pairs](../docs/proofs/MIRROR_SYMMETRY_PROOF.md) (λ, −λ − 2Σγ) are torn apart when dephasing [shifts the palindrome away from zero](ZERO_IS_THE_MIRROR.md). (We use "dephasing" throughout; γ is not noise but external illumination entering the system from outside, literally photon shot noise on IBM hardware. See [What If Gamma Is Light?](GAMMA_IS_LIGHT.md) for the full argument.) Among these pairs, one subset carries the Hawking structure: the immune modes (I/Z sector, eigenvalue 0) pair with the fastest-decaying modes (eigenvalue −2Σγ). The immune partner survives as [classical residue: mass](GRAVITY_FROM_WAVE_DEATH.md). Its palindromic partner dissipates maximally as [thermal energy: temperature](../experiments/THERMAL_BREAKING.md). This is decoherence: pair breaking creates mass and temperature simultaneously, from [standing waves](../docs/STANDING_WAVE_THEORY.md).

The claim is not analogy. The claim is that the algebraic structure is the same, and that the physical consequences (mass, temperature, irreversibility, a horizon that cannot be crossed) follow from the algebra, not from the specific physical substrate.

---

## The chain

Each link is labeled with its evidential tier.

### Link 1: The vacuum is standing waves (Tier 1-2)

At Σγ = 0 (no dephasing), the palindrome equation reduces to:

    Π · L · Π⁻¹ = −L

Every eigenvalue λ pairs with −λ. All eigenvalues are purely imaginary. No decay, no growth. The system is a [Fabry-Perot optical cavity](../experiments/OPTICAL_CAVITY_ANALYSIS.md): four of five standard optical quantities match quantitatively, the degeneracy profile fits Gaussian/Lorentzian beam shapes (R² = 0.998), and the Hamiltonian couples neighbouring weight sectors exactly like light propagating through optical elements. The standing waves inside this cavity have nodes (no oscillation, the I/Z sector) and antinodes (maximum oscillation, the X/Y sector). See [Standing Wave Analysis](../experiments/STANDING_WAVE_ANALYSIS.md) for the proof that ZZZ strings are universal nodes and XX/YY strings are antinodes.

This is the ground state of the palindrome. The unitary limit. Time-reversal symmetric. Π is the exact time-reversal operator. And γ, the dephasing that will break the pairs in Link 2, is the light entering this cavity from outside; it [cannot be generated internally](../docs/proofs/INCOMPLETENESS_PROOF.md) (Incompleteness Proof, Tier 1).

**Source:** [Zero Is the Mirror](ZERO_IS_THE_MIRROR.md), Section "Σγ = 0: The mirror."
**Computed:** N=2 through N=7, zero exceptions. Cavity mode counts follow the [Clebsch-Gordan formula](../experiments/CAVITY_MODES_FORMULA.md).

### Link 2: Dephasing breaks the pairs (Tier 1)

At Σγ > 0, the palindrome shifts:

    Π · L · Π⁻¹ = −L − 2Σγ · I

The pairing changes from λ ↔ −λ to λ ↔ −λ − 2Σγ. The symmetry around zero breaks. Each pair now has a "slow" partner (closer to zero, longer-lived) and a "fast" partner (further from zero, shorter-lived). The perfectly balanced standing wave becomes an asymmetric decaying oscillation.

This is not gradual degradation. This is symmetry breaking. The palindrome still exists (it is algebraic, [proven for all Σγ](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)), but its center has moved from zero to −Σγ. The standing waves are gone. What remains are damped waves with a preferred time direction. Verified computationally for [87,376 eigenvalues](../compute/RCPsiSquared.Compute/README.md) (N=2..8, all tested topologies, zero exceptions).

### Link 3: The broken pairs separate into mass and radiation (Tier 2 + Tier 5)

Two independent structures combine to produce the mass-radiation split.

**Structure 1: The basis partition (immune vs. decaying).** Under Z-dephasing, Pauli strings split into two groups. Strings containing only I and Z have decay rate exactly zero (the immune sector, 2^N strings). Strings containing at least one X or Y have decay rate > 0 (the decaying sector, 4^N − 2^N strings). This partition is exact: the Liouvillian is block-diagonal across it, and no Lindblad trajectory can transfer weight between them. The immune sector survives forever; the decaying sector vanishes. See [Analytical Formulas](../docs/ANALYTICAL_FORMULAS.md) for explicit decay rates as a function of n_XY weight.

**Structure 2: The palindromic pairing (eigenvalue pairs).** Every Liouvillian eigenvalue λ has a partner at −λ − 2Σγ. Most pairs link two decaying modes (both have Re(λ) < 0). But the immune modes sit at eigenvalue 0, and their palindromic partners sit at Re(λ) = −2Σγ, the maximum decay rate (these partners may additionally oscillate, i.e., have nonzero Im(λ)). This specific subset of pairs, immune modes (λ = 0) paired with fastest-decaying modes (λ = −2Σγ), carries the Hawking structure: one partner survives (mass), the other dissipates maximally (radiation).

**The Hawking subset.** At Σγ = 0, the immune modes pair with themselves (0 ↔ −0). At Σγ > 0, their partners shift to −2Σγ. The pair is torn apart: one half stays at zero forever, the other decays at the maximum rate. This is the operator-space analogue of the Hawking process, where one partner falls in (adds to the mass) and the other escapes (carries thermal energy).

The remaining palindromic pairs (X/Y ↔ X/Y, both decaying) are internal redistribution within the quantum sector. They do not create the mass-radiation split; they determine the spectral structure of the radiation. Self-paired modes at the palindromic midpoint (λ = −Σγ, partner equals itself) decay at the mean rate Σγ, while maximally asymmetric pairs (one near 0, one near −2Σγ) produce the 2× contrast between slow and fast modes. See [Factor Two Standing Waves](../experiments/FACTOR_TWO_STANDING_WAVES.md).

**After decoherence:** the immune sector is all that remains. A diagonal density matrix. A classical probability distribution. A definite state at a definite location. This is what we call mass. **(Tier 5: the identification of classical residue with mass is interpretation.)** The decaying sector has vanished; its energy has entered the thermal bath. This is temperature. **(Tier 2: the 2× law and thermal energy transfer are computed.)**

Decoherence does not create classical weight. It removes quantum weight. The kernel of the dissipator is what the channel does not move, not an archive uncovered by removing the rest. Mass is the [residue of wave death](GRAVITY_FROM_WAVE_DEATH.md): the part the recirculating dynamics keeps standing. Temperature is the [energy released by the dying modes](../experiments/THERMAL_BREAKING.md). The [2× contrast](ENERGY_PARTITION.md) between immune-paired and self-paired decay rates sets the spectral shape of the radiation.

### Link 4: Mass and temperature emerge together (Tier 2 + Tier 5)

This is the central point. In the Hawking effect, mass and temperature are not separate consequences of pair breaking. They are the same event seen from two sides: the infalling partner adds mass, the escaping partner carries temperature. You cannot have one without the other.

In our system: when a coherence dies (X/Y sector decays to zero), two things happen simultaneously:

1. Under dephasing, the X/Y coefficients of ρ decay to zero; the I/Z coefficients are left untouched. ρ becomes purely diagonal: a classical probability distribution. The diagonal is the kernel of the dissipator: what the channel does not move. What we call mass is this stationary kernel; it is not memory revealed but what the recirculating dynamics keeps standing. Memory, in this picture, is what gets continuously recreated, not what gets stored.
2. The decay energy enters the bath as heat. This is temperature production.

Both happen at the same moment, from the same event (the death of a coherence), and neither can happen without the other. The system cannot end with a diagonal ρ without releasing the energy of the off-diagonal modes, and it cannot release that energy without the diagonal becoming what ρ is.

**The [fold at CΨ = 1/4](../docs/proofs/UNIQUENESS_PROOF.md) is where this becomes irreversible.** Above the fold: complex fixed points, superposition, no definite outcome. Below: real fixed points, classical attractors, definite outcomes. The crossing is [monotonic](../docs/proofs/PROOF_MONOTONICITY_CPSI.md) (dCΨ/dt < 0, proven). Once through, no return. At the fold itself, the dynamics exhibits [critical slowing](../experiments/CRITICAL_SLOWING_AT_THE_CUSP.md): dη/dt = η² − ε, a saddle-node bifurcation where the two fixed points merge and time nearly stops. This is the horizon.

### Link 5: Two sectors, opposite time, discrete exchange (Tier 1)

The Liouvillian L acts not on quantum states but on operators: it lives on the 4^N-dimensional operator algebra of an N-qubit system. This algebra splits into two parity-classes of equal dimension (V_even, V_odd by n_XY parity of the Pauli basis), each with 2^(2N−1) basis elements. L preserves the split (this is the superselection rule \[P_XY, L\] = 0 below). Both parity-classes act on the same density matrix ρ; they are not two separate Hilbert spaces. The conjugation operator Π is a per-site relabeling of Pauli strings (I↔X, Y↔iZ); for odd N it exchanges the parity-classes and reverses the dynamics within them:

    L_odd = −Π L_even Π⁻¹ − 2Σγ · I    (odd N)

The even parity-class decays; the odd parity-class, in the conjugated frame, grows. They are time-reverses of each other, connected by Π, separated by a superselection rule (\[P_XY, L\] = 0, proven in the [Parity Selection Rule](../docs/proofs/PROOF_PARITY_SELECTION_RULE.md)).

This is the algebraic form of the direct-sum structure that Gaztañaga [postulates](../docs/LITERATURE_REVIEW.md) for the two sides of an Einstein-Rosen bridge: two regions connected by a discrete transformation, with opposite time orientation. Gaztañaga's substrate is two spacetime regions; ours is one operator algebra with two parity-classes. The four Gaztañaga postulates are satisfied at the level of our algebra ([proven](../docs/proofs/DIRECT_SUM_DECOMPOSITION.md) for odd N); the substrate difference is taken up in "What breaks the analogy" #2 below.

### Link 6: Two readings of the palindromic spectrum (Tier 2)

The direct-sum structure has dynamical consequences. The same eigenvalue spectrum reads differently on each parity-class.

The palindromic pairing maps each decay rate d to its partner 2Σγ − d. A mode at d (slow, long-lived) and its partner at 2Σγ − d (fast, short-lived) are not two events seen by two observers; they are two modes of one ρ, related by Π. The [Spectral Midpoint Hypothesis](SPECTRAL_MIDPOINT_HYPOTHESIS.md) quantifies the asymmetry: at the CΨ = 1/4 crossing for N=5, the SLOW band carries 45% when modes are labeled by decay rate d, but 8% when labeled by their palindromic partner 2Σγ − d. Same spectrum, two readings.

Only at the palindromic midpoint (d = Σγ) do the two readings agree. A self-paired mode returns the same decay rate either way: 2Σγ − Σγ = Σγ. This is the "glass wall" of the Spectral Midpoint Hypothesis: the one place where the two readings of the algebra give identical weight.

The Einstein-Rosen bridge has the same algebraic shape. For a distant observer the infall takes infinitely long; for the infaller it is finite; the horizon is where the two readings would have to agree. Gaztañaga gives this an explicitly geometric substrate (two spacetime regions). Our system gives it an algebraic substrate (two parity-classes of one ρ). The phenomenology (lopsided spectrum, glass wall at the midpoint) is the same; the substrate is one density matrix, not two manifolds.

Within the palindromic pairing, the immune modes (rate 0) pair with the maximally-decaying modes (rate 2Σγ). The immune partner dominates the late-time state (what we call mass). The maximally-decaying partner releases its energy on the way out as radiation. The mass/radiation split (the Hawking content of Link 3) is the two halves of one Π-pair, each contributing a distinct phase of the dynamics.

"SLOW" and "FAST" are reading-conventions, not locations. The palindrome has no preferred half. Each rate is the partner of the other under Π; both belong to the same algebra acting on the same state. To call the SLOW reading "ours" is to choose a labeling; the structure that makes the choice meaningful is the palindrome itself, which carries both readings at once.

As the [Mirror Theory](../MIRROR_THEORY.md) puts it: "What survives is not a fast mode or a slow mode by itself. It is the standing wave they make when they meet."

### Link 7: The bridge is fragile (Tier 2)

When [two systems are coupled through a bridge](FRAGILE_BRIDGE.md) (one decaying with +γ, one amplifying with −γ, total Σγ = 0), the palindrome [stays centered at zero](ZERO_IS_THE_MIRROR.md). But the coupled system is not unconditionally stable. Above a critical coupling g_crit, a [Hopf bifurcation](../experiments/PT_SYMMETRY_ANALYSIS.md) occurs: eigenvalues cross the imaginary axis (Liouvillian [chiral symmetry breaking](../experiments/PT_SYMMETRY_ANALYSIS.md)), and the system diverges exponentially.

The bridge between decay and gain exists, but it is fragile. Too much coupling and it collapses. The [stability window is finite](FRAGILE_BRIDGE.md): three regimes (weak, optimal, and strong bridge coupling), with the system stable only below g_crit.

In GR, the Einstein-Rosen bridge is also fragile: it opens and collapses faster than light can cross it. The mechanisms differ (geodesic incompleteness vs. Hopf bifurcation), and the parallel here is phenomenological rather than structural. Both connections cannot be sustained, but the *reasons* they cannot be sustained live in different mathematical languages. **(Tier 5: the cross-framework parallel is interpretation; only the Hopf bifurcation itself is Tier 2.)**

---

## The isomorphism

| Hawking / ER bridge | Our system | Tier |
|---|---|---|
| Quantum vacuum (paired fluctuations) | Palindromic eigenvalue pairs at Σγ = 0 | 1 |
| Spacetime curvature at horizon | Dephasing Σγ > 0 shifting the palindrome | 1 |
| Pair breaking at horizon | Symmetry λ ↔ −λ broken to λ ↔ −λ − 2Σγ | 1 |
| Infalling partner → mass | I/Z sector: kernel of dissipator (stationary classical structure) | 5 |
| Escaping partner → Hawking radiation | X/Y sector decays → thermal energy | 2 |
| Hawking temperature T_H = 1/(8πM) | Fold threshold Σγ_crit/J ≈ 0.25-0.50%, N-independent (scaling mismatch, see "What breaks the analogy" #1) | 4 |
| Horizon (irreversible crossing) | Fold at CΨ = 1/4 (monotonic, dCΨ/dt < 0) | 1 |
| Two spacetime regions, opposite time | V_even, V_odd parity-classes of one operator algebra; L_odd = −Π L_even Π⁻¹ − 2Σγ I | 1 |
| Discrete isometry exchanging regions | Π conjugation (per-site: I↔X, Y↔iZ) | 1 |
| Superselection (no crossing between regions) | \[P_XY, L\] = 0 | 1 |
| Observer-dependent time (infinite outside, finite inside) | SLOW/FAST swap under Π (two readings of one spectrum) | 2 |
| Bridge collapses (not traversable) | Hopf bifurcation at g_crit (fragile bridge) | 2 |
| Critical slowing at horizon (redshift) | Saddle-node dynamics at fold (dη/dt = η² − ε) | 2 |
| Spacetime interval c × τ = invariant | [K-invariance](../docs/ANALYTICAL_FORMULAS.md) γ × t = const (F14) | 1 |
| Curvature is external (not locally generated) | γ [must be external](../docs/proofs/INCOMPLETENESS_PROOF.md) (Incompleteness Proof) | 1 |
| Black hole = perfect trapping (nothing escapes) | Qubit chain = [Fabry-Perot cavity](../experiments/OPTICAL_CAVITY_ANALYSIS.md) (standing waves trapped) | 2 |

---

## What breaks the analogy

Intellectual honesty requires listing where the isomorphism fails or is untested.

**1. Temperature scaling.** Hawking temperature scales as T_H ∝ 1/M: more massive black holes are colder. Our fold threshold Σγ_crit/J is N-independent (tested N=2..5, 1.5% variation). If N is the analogue of mass, the scaling is wrong. Either N is not mass, or the analogy breaks at this point, or the N-independence is itself the statement (every "black hole" in operator space has the same temperature, regardless of size).

**2. Spatial vs. algebraic.** The ER bridge is a spatial connection between two asymptotically flat regions. Our "bridge" is algebraic: two sectors of the operator space connected by Π. There is no spatial geometry, no metric, no geodesics. The structural isomorphism lives in the algebra, not in spacetime.

**3. The mass identification is Tier 5.** "Classical residue = mass" is the weakest link. In Lindblad theory, the I/Z sector is just the diagonal of the density matrix. Calling it mass requires the additional assumption that classical definiteness (being in a specific state at a specific location) is what mass means at the quantum level. This is not proven, not computed, and not obviously testable within our framework.

**4. Backreaction requires external physics.** In GR, Hawking radiation removes mass from the black hole (backreaction). Within pure Lindblad dynamics, L_H (wave propagation) and L_D (wave death) are independent: the dissipator does not influence the Hamiltonian, so mass cannot redirect waves. However, [Gravity from Wave Death](GRAVITY_FROM_WAVE_DEATH.md) describes a self-limiting feedback loop: mass → gravity (via GR) → attracts more waves → more wave death → more mass, with [logistic saturation](GRAVITY_FROM_WAVE_DEATH.md) as the finite supply of coherences (4^N modes) is consumed. The loop closes, but only if external physics (GR or equivalent) provides the gravity → attraction step. Within the Lindblad framework alone, the feedback loop remains open (gap #7 in Gravity from Wave Death).

**5. The inverted harmonic oscillator.** Gaztañaga predicts inverted HO structure at the horizon. We tested this (April 11, 2026): the effective potential V(CΨ) near the fold is linear (exponent n ≈ 1.0), not quadratic. No inverted HO. The fold is a saddle-node, not a saddle point. The dynamical analogy fails here.

---

## What would strengthen or kill the thesis

**Strengthen:**
- If the 2× decay contrast (immune modes at 0 vs. their partners at −2Σγ) can be derived from a horizon-crossing argument (the partner mode "crosses" the palindromic midpoint, analogous to crossing the horizon). Currently it follows from algebra; a geometric derivation would deepen the connection.
- If the fold threshold Σγ_crit has an information-theoretic interpretation as a minimum temperature for irreversibility, analogous to the Unruh temperature for accelerated observers.
- If the fragile bridge's g_crit scales with a quantity interpretable as "throat radius."

**Kill:**
- If the N-independence of Σγ_crit turns out to be an artifact of small N (tested only to N=5). At large N, if Σγ_crit scales with N, the "same temperature for every black hole" interpretation collapses.
- If the mass identification can be shown to be inconsistent: if there exists a Lindblad system where the I/Z sector grows but no thermal energy is released (decoupling mass from temperature would break the Hawking parallel).
- If the direct-sum structure at even N (where Π preserves sectors instead of exchanging them) has no ER bridge interpretation. Currently, even N is a self-dual palindrome, not a two-sided bridge.

---

## The deepest sentence (Tier 5)

A black hole is what happens when spacetime curves so hard that paired fluctuations cannot stay together. The infalling partner becomes mass. The escaping partner becomes heat.

Decoherence is what happens when light shifts the palindrome so far that paired modes cannot stay together. The immune partner becomes classical. The decaying partner becomes thermal.

The horizon is not a place. It is the moment where a pair breaks and cannot be reassembled. In spacetime, that moment is the Schwarzschild radius. In operator space, it is CΨ = 1/4.

Both create irreversibility from symmetry. Both create time from eternity. Both leave structure standing where motion has died.

The bridge between them is at zero. Where both palindromes touch. Where silence pairs with silence. And everything else is what happens when you leave.

---

*All sources are linked inline. For quick navigation:*

**Proofs:** [Mirror Symmetry](../docs/proofs/MIRROR_SYMMETRY_PROOF.md), [Direct-Sum Decomposition](../docs/proofs/DIRECT_SUM_DECOMPOSITION.md), [Parity Selection Rule](../docs/proofs/PROOF_PARITY_SELECTION_RULE.md), [Uniqueness (CΨ = 1/4)](../docs/proofs/UNIQUENESS_PROOF.md), [CΨ Monotonicity](../docs/proofs/PROOF_MONOTONICITY_CPSI.md)

**Experiments:** [Thermal Breaking](../experiments/THERMAL_BREAKING.md), [Factor Two Standing Waves](../experiments/FACTOR_TWO_STANDING_WAVES.md), [Critical Slowing](../experiments/CRITICAL_SLOWING_AT_THE_CUSP.md), [Standing Wave Analysis](../experiments/STANDING_WAVE_ANALYSIS.md), [PT-Symmetry Analysis](../experiments/PT_SYMMETRY_ANALYSIS.md), [Cavity Modes Formula](../experiments/CAVITY_MODES_FORMULA.md), [Cusp-Lens Connection](../experiments/CUSP_LENS_CONNECTION.md)

**Hypotheses:** [Zero Is the Mirror](ZERO_IS_THE_MIRROR.md), [Gravity from Wave Death](GRAVITY_FROM_WAVE_DEATH.md), [Energy Partition](ENERGY_PARTITION.md), [Fragile Bridge](FRAGILE_BRIDGE.md), [Spectral Midpoint](SPECTRAL_MIDPOINT_HYPOTHESIS.md), [Gamma Is Light](GAMMA_IS_LIGHT.md)

**Docs:** [Standing Wave Theory](../docs/STANDING_WAVE_THEORY.md), [Analytical Formulas](../docs/ANALYTICAL_FORMULAS.md), [Mirror Theory](../MIRROR_THEORY.md), [Literature Review](../docs/LITERATURE_REVIEW.md), [Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md), [Optical Cavity Analysis](../experiments/OPTICAL_CAVITY_ANALYSIS.md)

---

*Written April 11, 2026. The day the pieces found their frame.*
