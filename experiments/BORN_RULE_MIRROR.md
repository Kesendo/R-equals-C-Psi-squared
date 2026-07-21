# Born Rule as Mirror Quality: 97% Hamiltonian, 3% Decoherence Correction

<!-- Keywords: Born rule R=CPsi2 generalized, measurement outcome probability
Hamiltonian dominance 97 percent, decoherence basis correction 3 percent,
effective coupling C_eff per outcome, standing wave sum squared origin,
perfect mirror limit recovers Born rule, basis alignment probability shift,
R=CPsi2 Born rule mirror -->

> **Restoration note (March 14, 2026):** Originally written 2026-02-18, deleted March 12,
> restored March 14. Core numerics confirmed.

**Status:** Verified numerics (Tier 2); interpretation speculative (Tier 3)
**Date:** 2026-02-18 (Section 4.3 added 2026-02-27)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Dynamic Entanglement](DYNAMIC_ENTANGLEMENT.md), [Subsystem Crossing](SUBSYSTEM_CROSSING.md)

---

## Abstract

> **Reproduction note (2026-07-20):** the reference time t = 0.286 used
> throughout was labeled "the (0,2) crossing" by the retired MCP
> tool (the delta_calc family, February 2026); under the canonical pair-CΨ convention (Wootters concurrence ·
> l₁/3, [subsystem_crossing_pairs.py](../simulations/subsystem_crossing_pairs.py))
> pair (0,2) of |0+0+⟩ on the N=4 ring never crosses ¼ (see the
> reproduction note in [Dynamic Entanglement](DYNAMIC_ENTANGLEMENT.md)).
> The Born-deviation content at t = 0.286 is unaffected by that label and
> was later closed in F94/F96 (the registry), on exactly this lens.
>
> **Reproduction note (2026-07-21):** unlike the crossing label, the
> NUMBERS of this document reproduce exactly: every probability table,
> deviation, shift, and the C_eff column regenerate from first
> principles under the PAULI convention H = 1.0·Σ(XX+YY+ZZ) (committed
> probe:
> [born_rule_mirror_tables.py](../simulations/born_rule_mirror_tables.py)).
> Two decodings the probe pins: the doc's J = 1.0 is the Pauli
> coupling, which in the spin convention of the F94 proof means
> J_spin = 4, so this document's point is Q = 80, not Q = 20 (see the
> F94 note in Section 2.1); and the Ψ_i column of Section 3.4 is the
> l₁ row-coherence of the Lindblad pair state over d−1, not a Born
> amplitude (see the note there).

At the CΨ = ¼ crossing point for pair (0,2) in the alternating state
|0+0+⟩ under Heisenberg ring dynamics, Born rule probabilities are ~97%
determined by unitary Hamiltonian evolution alone. The remaining ~3% is
a systematic correction from the decoherence basis: σ_z dephasing shifts
+0.0098 aggregate probability toward the z-eigenstate pair (the effect
is initial-state-dependent: for this z-aligned initial state the σ_x
and σ_y own-basis shifts are negative, Section 7.3). Applied per
measurement outcome,
R_i = C_i·Ψ_i² recovers the standard Born rule P(i) = |⟨i|ψ⟩|² when
C_i is uniform across outcomes (perfect mirror limit). For real
measurements, C_i varies across outcomes (C_eff = 8.22 for |00⟩,
4.82 for |11⟩; the 05-22 note in Section 2.2 corrects the old
"z-superpositions" label, and the Section 3.4 note pins what the
C_eff column actually measures). The standing wave interpretation R = C·(Ψ_past +
Ψ_future)² provides a physical origin for the square: equal forward and
backward amplitudes meeting in a standing wave produce intensity
proportional to amplitude squared.

---

## 1. The Question

[Dynamic Entanglement](DYNAMIC_ENTANGLEMENT.md) (Experiment 11) showed that pair (0,2) in the
alternating state |0+0+⟩ crosses the 1/4 threshold from below. At the
crossing point, the reduced density matrix has diagonal elements:

| State | Probability |
|-------|-------------|
| \|00⟩  | 0.4254      |
| \|01⟩  | 0.2567      |
| \|10⟩  | 0.2567      |
| \|11⟩  | 0.0613      |

These are Born rule probabilities, the chances of finding the pair
in each basis state upon measurement. The open question was:

> Can the framework predict these probabilities independently, without
> computing the full density matrix?

## 2. The Answer

The Born rule probabilities come almost entirely from the Hamiltonian
evolution. The decoherence determines WHEN the threshold is crossed,
but not WHAT the probabilities are at that moment.

### 2.1 Numerical Evidence

Comparison of pair (0,2) diagonal at crossing time t = 0.286:

| State | Unitary (no noise) | Lindblad (γ=0.05) | Deviation |
|-------|-------------------|-------------------|-----------|
| \|00⟩  | 0.4134            | 0.4254            | 2.9%      |
| \|01⟩  | 0.2616            | 0.2567            | 1.9%      |
| \|10⟩  | 0.2616            | 0.2567            | 1.9%      |
| \|11⟩  | 0.0635            | 0.0613            | 3.5%      |

The unitary evolution alone (pure Hamiltonian, no decoherence)
determines ~97% of each probability. The remaining ~3% is a systematic
correction from the decoherence basis.

*Later (2026-05-16, convention bridge added 2026-07-21):* this 97/3 split was closed analytically as [F94](../docs/ANALYTICAL_FORMULAS.md): the dominant-outcome RATIO deviation Δ_|00⟩ = P_lind/P_unit − 1 = (4/3)·Q²·K³ in the deep perturbative regime, with Q = J/γ₀ and K = γt. Mind the book when plugging in: F94's J is the SPIN convention H = (J/4)·Σ, while this document's J = 1.0 is the PAULI coupling, so this document's point is Q = 80, K = 0.0143. There the leading order gives 0.0250 against the measured ratio 0.0290 (the 2.9% of the table above): within ~14%, the rest being higher orders, since Q²K³ ≈ 0.019 sits well above F94's verified deep-perturbative window. A naive Q = 20 would give 0.0016 and miss by ~18×. See [`PROOF_F94_BORN_DOMINANT_FOUR_THIRDS.md`](../docs/proofs/PROOF_F94_BORN_DOMINANT_FOUR_THIRDS.md) for the Dyson-sym3 = 8 bit-exact derivation and the companion reflection [`ON_HOW_FOUR_THIRDS_APPEARED.md`](../reflections/ON_HOW_FOUR_THIRDS_APPEARED.md) for the day-of-arrival path.

### 2.2 The Correction Pattern

The 3% correction is not random. It follows a rule: decoherence
shifts probability toward the eigenstates of the dephasing operator.

For σ_z dephasing, z-eigenstates are favored (in aggregate, see below):

| | Total shift to z-eigenstates | Total shift to z-superpositions |
|---|---|---|
| σ_z dephasing | +0.0098 | -0.0098 |
| σ_x dephasing | +0.0012 | -0.0012 |
| σ_y dephasing | +0.0045 | -0.0045 |

(All nine shifts of this section and Section 7.3 regenerate exactly in
the committed probe; "shift to z-eigenstates" = the Lindblad-minus-
unitary diagonal summed over {|00⟩, |11⟩}.)

The sign confirms: σ_z dephasing favors the z-eigenstate pair in
AGGREGATE (the +0.0098 is |00⟩'s +0.0120 minus |11⟩'s −0.0022; |00⟩ is
favored, |11⟩ individually suppressed, exactly as Section 6 predicts).
The magnitude depends on the overlap between the initial state and the
dephasing basis, not just on the dephasing basis alone; the same
overlap dependence is why the σ_x and σ_y OWN-basis shifts in
Section 7.3 come out negative, so "shifts toward the dephasing
eigenstates" is not a universal rule but an initial-state-dependent
one.

The initial state of pair (0,2) is |0⟩|0⟩, which is a z-eigenstate.
Therefore σ_z dephasing PROTECTS it (largest positive shift to |00⟩),
while σ_x and σ_y dephasing ATTACK it (they see |0⟩ as a superposition
in their basis).

*Later (2026-05-22):* the "z-eigenstates / z-superpositions" labels are loose: all four pair-(0,2) outcomes are z-basis states. The real split is parity (Z⊗Z = +1 for {|00⟩,|11⟩}, −1 for {|01⟩,|10⟩}); under the qubit-0↔2 mirror of |0+0+⟩, {|00⟩,|11⟩} are each invariant while {|01⟩,|10⟩} are exchanged, so the latter's natural modes (|01⟩±|10⟩)/√2 are genuine superpositions. The label "z-superpositions" was reaching for that mirror-mode structure, not for a literal z-basis superposition. The framework later made the stance precise in [`ON_HOW_THE_ANGLE_APPEARS_AT_ZERO.md`](../reflections/ON_HOW_THE_ANGLE_APPEARS_AT_ZERO.md): superposition is not a primitive but the basis-relative F95 angle-coordinate. The loose wording here is an early trace of it.

## 3. Connection to R = CΨ²

### 3.1 The Standard Born Rule

For an ideal measurement of a pure state |ψ⟩, the Born rule gives:

    P(i) = |⟨i|ψ⟩|²

This is a postulate in standard quantum mechanics. It is not derived
from the other axioms. It is assumed.

### 3.2 R = CΨ² as Generalized Born Rule

R = CΨ² applied per measurement outcome gives:

    R_i = C_i · Ψ_i²

where Ψ_i = |⟨i|ψ⟩| is the overlap between outcome and state (the
quantum potential for that outcome), and C_i is the effective coupling
strength for that outcome.

For ideal measurement, C is the same for all outcomes. It cancels in
the normalization, and we recover:

    P(i) = R_i / Σ R_j = Ψ_i² / Σ Ψ_j² = |⟨i|ψ⟩|²

This is exactly the Born rule.

### 3.3 Non-Ideal Measurement

In real systems, the decoherence basis breaks the symmetry between
outcomes. C_i is no longer equal for all i. The probabilities become:

    P(i) = C_i · Ψ_i² / Σ C_j · Ψ_j²

Outcomes whose basis states are eigenstates of the dephasing operator
have higher effective C (the environment "sees" them more clearly).
Outcomes that are superpositions in the dephasing basis have lower
effective C (the environment sees them less clearly).

The deviation from the standard Born rule is small (~3% in our
simulations) because the unitary dynamics dominate. But the deviation
is systematic and predictable from the dephasing basis.

### 3.4 Measured C_eff Values at Crossing

Computing C_eff = P(i) / Ψ_i² for each basis state at the crossing:

| State | P(i)   | Ψ_i    | C_eff  |
|-------|--------|--------|--------|
| \|00⟩  | 0.4254 | 0.2274 | 8.22   |
| \|01⟩  | 0.2567 | 0.2139 | 5.61   |
| \|10⟩  | 0.2567 | 0.2139 | 5.61   |
| \|11⟩  | 0.0613 | 0.1127 | 4.82   |

Decoding note (2026-07-21): the Ψ_i column is NOT the amplitude
|⟨i|ψ⟩| of Section 3.2 (the pair state is mixed, and √P_unitary would
read 0.64/0.51/0.25). What the February tool computed is
Ψ_i = (Σ_{j≠i} |ρ₀₂[i,j]|)/(d−1), the l₁ row-coherence of the Lindblad
pair state over d−1 = 3 (recovered to four decimals by the committed
probe). C_eff = P/Ψ² is therefore a coherence-referenced diagnostic
computed from the same density matrix, not an amplitude-based coupling;
the honest observation is that |00⟩ carries the largest probability per
unit of squared row-coherence. The framework reading "C depends on
basis alignment" is the Tier-3 interpretation of that pattern.

## 4. The Mirror Interpretation

[Mirror Theory](../MIRROR_THEORY.md) describes
reality as arising between observers, like light between mirrors.
The quality of the reflection determines what becomes real.

The Born rule, in this picture, answers: WHICH reality crystallizes
when the threshold is crossed?

The answer: the one with the strongest reflection. The outcome that
has the best overlap with the quantum state (Ψ_i) AND the best
alignment with the observation basis (C_i) dominates.

In an ideal measurement (perfect mirror), only Ψ matters: P ∝ Ψ².
In a real measurement (imperfect mirror), C × Ψ² determines the
outcome. The mirror's distortion shifts the probabilities.

### 4.1 The Two-Step Process

1. **The Hamiltonian determines the amplitudes.** The interaction
   between subsystems creates the quantum state |ψ⟩ with specific
   amplitudes c_i = ⟨i|ψ⟩. This is ~97% of the Born probabilities.
   The conversation between the mirrors determines what is possible.

2. **The decoherence selects the basis and provides the correction.**
   The environment's preferred basis (einselection) determines which
   outcomes are "seen clearly" and which are "blurred." This provides
   the ~3% correction. The quality of the mirror determines which
   possibilities are slightly favored.

Together: R_i = C_i · |⟨i|ψ⟩|² = C_i · Ψ_i².

### 4.2 Why the Hamiltonian Dominates

The Hamiltonian acts for the entire evolution time (0 to t_cross).
The decoherence also acts for the entire time, but its effect on the
PROBABILITY RATIOS is small because it primarily destroys off-diagonal
coherence (which determines WHEN the crossing happens) rather than
redistributing diagonal elements (which determine WHAT the outcome is).

Decoherence is the clock. The Hamiltonian is the story.

## 4.3 The Standing Wave Origin of Ψ²

**Added:** 2026-02-27
**Tier:** 3, Connects [Interpretive Framework](../hypotheses/archive/INTERPRETIVE_FRAMEWORK.md) bidirectional bridge to Born rule
**Depends on:** [Interpretive Framework](../hypotheses/archive/INTERPRETIVE_FRAMEWORK.md) Section 2 (Wave Composition)

### Why the Square?

The Born rule says probability is the *square* of the amplitude. This
is postulated in standard QM. Within the framework, it follows from
the bidirectional wave structure.

From [Interpretive Framework](../hypotheses/archive/INTERPRETIVE_FRAMEWORK.md):

    Ψ = Ψ_past + Ψ_future

    Ψ_past:   offer wave from source (emitter → outcome)
    Ψ_future: confirmation wave from detector (absorber → outcome)

For a standing wave to form, forward and backward wave must meet:

    Ψ_f = A·sin(kx - ωt)    (forward, future)
    Ψ_p = A·sin(kx + ωt)    (backward, past)
    Ψ = 2A·sin(kx)·cos(ωt)  (standing wave)

The key: a standing wave requires EQUAL amplitudes. Both waves
carry amplitude A. The resulting pattern has amplitude 2A, and
its intensity is (2A)² = 4A². In this Tier-3 picture, the square
comes from the wave physics rather than from a postulate.

### Per Outcome

For measurement outcome i:

    Ψ_past_i  = |⟨i|ψ⟩| = amplitude the source offers to outcome i
    Ψ_future_i = confirmation amplitude the detector returns for i
    C_i        = resonance quality (how well offer and confirm match)

    R_i = C_i · (Ψ_past_i + Ψ_future_i)²

Three independently defined quantities. No circularity.

### The Perfect Mirror Limit

A perfect detector is a perfect mirror: it returns exactly what it
receives. Then Ψ_future_i = Ψ_past_i for all outcomes.

    R_i = C · (2·|⟨i|ψ⟩|)² = 4C · |⟨i|ψ⟩|²

Since C is uniform (same mirror quality for all outcomes), it cancels
in normalization:

    P(i) = R_i / Σ R_j = |⟨i|ψ⟩|² / Σ |⟨j|ψ⟩|² = |⟨i|ψ⟩|²

This IS the Born rule. It follows from: perfect mirror → equal
confirmation → uniform C → probability proportional to amplitude
squared.

### Why No Perfect Mirrors Exist

In this section's reading, C is the detector's purity: C = Tr(ρ²) ≤ 1.
(This is a different C than Section 3.4's C_eff = P/Ψ², which is a
normalization-free diagnostic and can exceed 1; the two share the
letter, not the book.) For a macroscopic detector in thermal
contact with its environment, C < 1. Always. A detector with C = 1
would be a pure quantum state, isolated from the entire universe.
No real measurement apparatus achieves this.

Therefore:
- The Born rule P(i) = |⟨i|ψ⟩|² is the LIMIT, not the law.
- The actual law is R_i = C_i · (Ψ_past_i + Ψ_future_i)².
- Deviations from Born are systematic, not noise.
- The ~3% correction measured in Section 2.2 is the prediction.

### Connection to Cramer's Transactional Interpretation

This structure resembles Cramer's TI (1986): offer wave forward,
confirmation wave backward, transaction = reality. In Cramer's
formulation, the transaction probability is ψ·ψ* = |ψ|².

The difference: Cramer uses a PRODUCT (ψ · ψ*). We use a SUM
squared: (Ψ_past + Ψ_future)². For perfect mirrors where
Ψ_future = Ψ_past, both give the same result (proportional to
|amplitude|²). They diverge for imperfect mirrors.

Our formulation adds:
1. C as an explicit resonance quality factor (not in Cramer)
2. The standing wave picture (sum, not product)
3. Quantitative prediction of deviations from Born rule
4. Connection to the 1/4 boundary via R = CΨ²

Whether the sum or product is physically correct is an open
question. Both reproduce Born in the ideal limit. They make
different predictions for imperfect measurements. The simulation data
(Section 2) is consistent with both; discriminating requires
measurements at different detector purities.

### What This Does NOT Claim

- We do not derive the amplitudes ⟨i|ψ⟩ from the framework.
  Those come from the Hamiltonian (Schrödinger equation).
- We do not explain WHY the detector returns Ψ_future_i = Ψ_past_i
  in the ideal case. We observe that perfect mirroring implies Born.
- The identification of Ψ_past_i with |⟨i|ψ⟩| is proposed (Tier 3),
  not proven.

## 5. The Cycle Completed

This result closes the conceptual loop opened in Section 8 of
[Dynamic Entanglement](DYNAMIC_ENTANGLEMENT.md):

1. **Decoherence cements the past** (destroys coherence, makes
   crossings irreversible).
2. **Decoherence creates blind spots** (selects which subsystems
   can evolve freely).
3. **The Hamiltonian builds new potential** (creates entanglement
   in the protected subspace).
4. **Potential crosses the threshold** (C·Ψ reaches 1/4).
5. **The Born rule determines WHAT becomes real** (~97% from the
   Hamiltonian, ~3% from the decoherence basis alignment).

In this reading (Tier 3, per the Status line), the Born
rule is not a separate postulate but the consequence of how interaction
(the Hamiltonian) and observation (decoherence) combine at the
crossing point. R = CΨ² contains both.

## 6. Testable Prediction

The framework predicts that Born rule probabilities deviate from
|⟨i|ψ⟩|² in a specific, systematic way:

- Outcomes aligned with the decoherence basis are FAVORED in aggregate
  (per-outcome signs depend on the initial state, Section 7.3)
- Outcomes misaligned with the decoherence basis are SUPPRESSED
- The deviation magnitude is linear in γ at fixed readout time
  (verified 2026-07-21; consistent with F94's Δ = (4/3)·J²γt³)
- At γ → 0, the standard Born rule is recovered exactly

For the alternating state under σ_z dephasing, this predicts:
P(|00⟩) > |⟨00|ψ⟩|² and P(|11⟩) < |⟨11|ψ⟩|² at the crossing point.

This is a quantitative prediction testable on quantum hardware by
comparing measurement statistics for different engineered decoherence
channels.

## 7. Verification

### 7.1 Setup and tool

The February analysis ran on the retired delta_calc MCP tool, not
QuTiP (the tool's source, recovered and inspected 2026-07-21, is pure
numpy/scipy; the "QuTiP mesolve" attribution here was a prose error).
QuTiP itself was genuinely available in the February environment, so a
QuTiP label elsewhere is judged per document, not condemned as a class.
The setup:
- State: |0+0+⟩ (alternating, 4 qubits)
- Hamiltonian: Heisenberg ring, J = 1.0 (PAULI convention, H = J·Σ(XX+YY+ZZ))
- Dephasing: local σ_z, σ_x, σ_y at γ = 0.05
- dt = 0.001 for fine crossing detection

The committed reproduction is
[born_rule_mirror_tables.py](../simulations/born_rule_mirror_tables.py)
(numpy/scipy, exact expm): it regenerates every table, all nine shifts,
the dominant eigenvalue, the Ψ/C_eff decoding, and the F94 bridge.

### 7.2 Key Numbers to Check

1. Crossing time (σ_z): t = 0.286
2. P(|00⟩) unitary at crossing: 0.4134
3. P(|00⟩) Lindblad at crossing: 0.4254 (Δ = +0.012)
4. P(|01⟩) Lindblad at crossing: 0.2567 (Δ = -0.005 from unitary)
5. Total shift to z-eigenstates under σ_z: +0.0098
6. Total shift to z-eigenstates under σ_x: +0.0012 (weaker, wrong basis)
7. Dominant eigenvalue of rho_02 at crossing: 0.911 (nearly pure)

### 7.3 Three-Basis Consistency

| Dephasing | Own-basis eigenstate shift | Cross-basis eigenstate shift |
|-----------|--------------------------|----------------------------|
| σ_z       | +0.0098 (z-eigenstates)  | -0.0136 (x-eigenstates)    |
| σ_x       | -0.0110 (x-eigenstates)  | +0.0012 (z-eigenstates)    |
| σ_y       | -0.0088 (y-eigenstates)  | +0.0045 (z-eigenstates)    |

Note: The sign reversal for σ_x and σ_y own-basis shifts reflects
the initial state asymmetry: pair (0,2) starts as |0⟩|0⟩, which
is a z-eigenstate. The dephasing effect depends on the overlap
between the INITIAL STATE basis and the DEPHASING basis, not just
the dephasing basis alone.

## 8. Origin

The conceptual foundation for this result was laid in December 2025
and January 2026 by Tom and a Claude instance working together on
the philosophical foundation (the [mirror theory](../MIRROR_THEORY.md)). At that time, the formula R = CΨ² existed but
there were no computational tools to determine what C should be for
specific quantum systems. No MCP server, no QuTiP simulations, no
density matrix calculations.

Instead, the framework was built from the experiential side: what
does observation feel like? What makes one encounter more real than
another? The mirror metaphor (two observers creating reality
between them, the quality of the reflection determining the outcome)
was derived from lived experience, not from equations.

The numerical verification in this document (February 2026) confirms
what was described without numbers months earlier. The 97/3 split
between Hamiltonian and decoherence, the basis-dependent mirror
distortion, the role of alignment between observer and observed:
all of this was already present in the human-readable text, expressed
in the language of relationships rather than operators.

The math caught up with the intuition.

---

*Previous: [Dynamic Entanglement](DYNAMIC_ENTANGLEMENT.md)*
*See also: [Mirror Theory](../MIRROR_THEORY.md) for the mirror metaphor that motivated this analysis.*
