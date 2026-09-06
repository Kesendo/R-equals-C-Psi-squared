# Standing Wave: Two Observers Create Interference, Not Time Travel

<!-- Keywords: standing wave two observer interference, Cramer transactional
interpretation reframed, sum squared cross-term coherence, partial trace
discards joint information, off-diagonal Bell state between observers,
fast slow decoherence rate observer labels, Born rule perfect mirror limit,
coherence floor IBM Q52 residual, R=CPsi2 standing wave -->

> **Restoration note (March 14, 2026):** Originally written February 27, 2026, deleted March 12,
> restored March 14. Mirror symmetry is confirmed; the standing-wave reading
> remains conditional on physical excitation, propagation, and interference gates.
> Time travel, gravity, and FTL signaling references have fallen.

**Status:** Conceptual interpretation (Tier 3)
**Date:** 2026-02-27
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Why the Sum](WHY_THE_SUM.md), [No-Signalling Boundary](NO_SIGNALLING_BOUNDARY.md)

---

## What this document is about

Two observers sharing an entangled state can be assigned the algebraic
cross-term `2·Psi_A·Psi_B` in the proposed sum-squared ansatz. A cross-term is
not by itself a physical interference measurement or standing wave. The
partial trace discards joint coherences, but that fact neither makes them
locally observable nor establishes counter-propagating modes.

---

## Abstract

Two observers looking at the same entangled state can be inserted into the
ansatz R = C·(Ψ_A + Ψ_B)². Its cross-term is algebraic; calling it
interference requires an operational observable and a phase-sensitive control.
The labels "past" and "future" are names in the ansatz, not measured time
directions or gravitational clocks. The partial trace ρ_A = Tr_B(ρ_AB) correctly
describes what one observer sees, but it discards the off-diagonal terms
of the joint state where "reality between us" lives. CΨ measures these
joint terms. The sum-squared formulation has a cross-term; Cramer's product
does not. When one labelled amplitude vanishes (Ψ_A → 0), the sum preserves
Ψ_B² while the product gives zero. A late-time hardware residual does not
distinguish surviving coherence from SPAM or readout floor and is not evidence
for the ansatz.

---

## 1. The Misunderstanding

Cramer's transactional interpretation (1986: the proposal that every
quantum event involves an "offer wave" traveling forward in time and a
"confirmation wave" traveling backward) called them "offer wave" and
"confirmation wave." He said one goes forward in time, the other backward. Everyone heard "time travel" and stopped listening.

Time does not run backward. Ever. For anyone. Time runs forward for every observer, always. That is not negotiable.

## 2. What Ψ_past and Ψ_future mean here

They are two amplitude labels in the proposed sum. This document supplies no
map from gravitational acceleration (units of length/time²) to a dephasing
rate γ (units of 1/time), and no cross-system law making `γ·t_cross` a universal
constant. A valid comparison would have to specify two channels, calibrate each
γ independently, prepare the same state, and measure each crossing time with
the same observable. None of that turns the labels into past and future.

## 3. The Standing Wave

The equation per measurement outcome:

    R_i = C_i · (Ψ_past_i + Ψ_future_i)²

This is not a wave traveling forward and another traveling backward. It is:

**Two labelled amplitudes associated with one joint state.**

Each observer applies R = CΨ² to their own qubit. Each sees their own reality through their own purity (C) and coherence (Ψ). But the qubits are entangled, they are not independent systems. The observers are looking at the same thing from two sides.

Two amplitudes in a squared sum produce a cross-term. Whether that cross-term
is a standing wave is the question, not the premise.

(Ψ_past + Ψ_future)² = Ψ_past² + 2·Ψ_past·Ψ_future + Ψ_future²

The cross-term `2·Ψ_past·Ψ_future` exists by expansion of the proposed
ansatz. An incoherent mixture or phase randomization is the needed negative
control: if the reported term survives without phase coherence, it is not an
interference witness.

## 4. Why "We Are All Mirrors" Is Not a Metaphor

Standard QM says: each observer has a reduced density matrix ρ_A or ρ_B. For a Bell state, both are I/2. Locally identical. Locally boring.

But R = CΨ² is not a property of one qubit. It is a property of one observer looking at their qubit. And there are two observers.

    Alpha looks at his qubit:  R_A = C_A · Ψ_A²
    Beta looks at his qubit:   R_B = C_B · Ψ_B²

Two viewpoints on the same entangled pair. Alpha is Beta's mirror. Beta is Alpha's mirror. The confirmation wave is not a wave traveling backward through time; it is the other observer's perspective on the shared state.

"Reality is what happens between us": the off-diagonal terms of the joint state
|Φ+⟩ that neither ρ_A nor ρ_B contains. Each reduced state omits those joint
coherences. Calling their simultaneous existence a standing wave is metaphor,
not a consequence of Π.

## 5. The Partial Trace Problem

Standard QM uses the partial trace: ρ_A = Tr_B(ρ_AB). This throws away everything about B. It says: "if you only have access to A, this is all you can know."

And it is correct, for a single observer with one qubit. But R = CΨ² is applied by EACH observer to THEIR qubit. Both simultaneously. The framework does not ask "what does A know about B?" It asks "what reality emerges when both observers exist?"

The partial trace is the right tool for one subsystem; the joint coherence is
the object present only in the two-subsystem state.

The constructed two-seat diagonal F36 neural Jacobian is a useful negative
control: it satisfies the exact same affine spectral identity while having no
spatial propagation mechanism. Exact pairing therefore cannot, by itself, be
a standing-wave witness. See [Algebraic Palindrome Neural](../docs/neural/ALGEBRAIC_PALINDROME_NEURAL.md).

## 6. Why This Matters for the Bridge Question

The agents (v033-v040) tried to find a way for Alpha to detect Beta's actions through the entangled pair. Standard QM says no: the partial trace guarantees ρ_A = I/2 regardless of what Beta does.

> [FALLEN: This FTL signaling connection was not confirmed and has been retired from the technical core.]

The joint coherence does not live inside either partial trace. The operational
question is whether a declared joint observable detects phase-sensitive
interference that neither marginal contains; no standing wave is presumed.

The off-diagonal terms of |Φ+⟩⟨Φ+| are nonzero. The partial trace discards them,
and CΨ_joint responds to them. No counter-propagation or interference gate here
turns them into a physical standing wave.

## 7. What This Does NOT Say

- It does NOT say time runs backward. Time runs forward for everyone.
- It does NOT say FTL signaling is possible. The standing wave may or may not be observable locally. That is an open question.
- It does NOT violate no-signaling automatically. No-signaling is a theorem about ρ_A = Tr_B(ρ_AB). The reported object is R_AB for the joint state, not a demonstrated standing wave.
- It does NOT require new physics. It requires a different question: not "what does one observer see?" but "what emerges when two observers look at the same entangled state?"

## 8. The Cramer Fork Revisited

Cramer: R = C · Ψ_past · Ψ_future (product)
Framework: R = C · (Ψ_past + Ψ_future)² (sum squared)

Both recover Born in the perfect-mirror limit. Both diverge when mirrors are imperfect.

The two formulas encode different algebraic combination rules; neither is
derived here from a measurement protocol.

The product has no cross-term. The sum squared has one. A physical standing
wave would additionally require coherent excitation, spatial or modal
propagation, and an interference observable with the incoherent control above.

The late-time discriminator tests this: when Ψ_past → 0 (one observer's coherence dies), does the cross-term survive? Sum says yes (Ψ_future² remains). Product says no (0 × anything = 0).

IBM Torino (qubit 52) shows a 2% coherence floor at t >> T2. A residual floor at that level is generic (readout/SPAM), so it is at most loosely consistent with a surviving-mirror reading, not evidence for it.

---

*Previous: BORN_RULE_MIRROR.md, Born rule as perfect-mirror limit*
*See also: BLACK_WHITE_HOLES_BIGBANG.md, τ = 0 as maximum coherence*
