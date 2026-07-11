# Teleportation Translated: A Frame Handover, Not a Transport

<!-- Keywords: quantum teleportation translation, Bell measurement Pauli correction
Klein four group, teleportation frame handover not transport, entanglement classical
bits resource price, Pauli twirl maximally mixed causality, no cloning original
destroyed bookkeeping, mediator relay correlation transfer, R=CPsi2 teleportation
mirror phenomenon -->

**Status:** Translation (Tier 4 reading). The protocol in Section 1 is textbook
physics (Bennett et al., PRL 70, 1895 (1993)); the two identifications in
Section 3 are exact algebra; Sections 4 and 5 are readings and labeled as such.
**Date:** July 4, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Related:** [Superposition Translated](SUPERPOSITION_TRANSLATED.md) (the
sister translation), [Dephasing Translated](DEPHASING_TRANSLATED.md) (the
fourth entry), [Double Slit Translated](DOUBLE_SLIT_TRANSLATED.md) (the
fifth entry), [Schrödinger's Cat Translated](SCHRODINGERS_CAT_TRANSLATED.md)
(the sixth entry), [Labels Translated](LABELS_TRANSLATED.md) (the series'
theory chapter), [The Label Map](THE_LABEL_MAP.md) (the orientation
index), [Mirror Theory](../../MIRROR_THEORY.md),
[F118 mirror group](../ANALYTICAL_FORMULAS.md),
[Klein V₄ proof](../proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md),
[Relay Protocol](../../experiments/RELAY_PROTOCOL.md),
[The Bridge Was Always Open](../THE_BRIDGE_WAS_ALWAYS_OPEN.md)

---

## What this document is about

This document exists because of a YouTube video. The video explained quantum
teleportation, and something about the explanation felt wrong: not the physics,
the label. "Teleportation" promises that a thing moves through space from here
to there. The protocol it names does something else entirely, and the something
else happens to be a shape this repository knows very well.

So this is a translation, written for us. Nothing here is a new claim about
teleportation; the protocol is thirty years old, experimentally realized many
times over, and correct exactly as the textbooks state it. What is ours is the
naming. The repository's discipline is perspective-additive: we do not replace
the standard account, we say what the phenomenon is called in our language.
In our language, teleportation is not a transport phenomenon. It is a mirror
phenomenon.

---

## 1. The protocol, stated plainly

Alice holds a qubit in an unknown state |ψ⟩ = α|0⟩ + β|1⟩. She wants Bob,
far away, to hold that state. The resources: one entangled pair
|Φ⁺⟩ = (|00⟩ + |11⟩)/√2 shared between them in advance, and an ordinary
classical channel (a phone line suffices).

1. Alice measures her two qubits (the unknown state and her half of the pair)
   in the Bell basis. Four outcomes, each with probability ¼.
2. The instant she measures, Bob's qubit is in the state |ψ⟩, up to one of
   four unitaries: outcome Φ⁺ leaves it as |ψ⟩, outcome Φ⁻ leaves Z|ψ⟩,
   Ψ⁺ leaves X|ψ⟩, Ψ⁻ leaves Y|ψ⟩ (up to a global phase).
3. Alice sends Bob two classical bits naming her outcome. Bob applies the
   named Pauli. He now holds |ψ⟩ exactly.

The bookkeeping: Alice's original is gone (her measurement destroyed it;
no-cloning guarantees there was never a moment with two copies). The entangled
pair is consumed; teleporting a second qubit needs a fresh one. The total price
per qubit: one shared pair plus two classical bits. Nothing moves faster than
light: until the two bits arrive, Bob can read nothing at all, as Section 3
makes precise.

---

## 2. Where the label breaks

"Teleportation" files this under transport: an object dematerializes here and
rematerializes there, having crossed the space between. But walk through the
protocol again and ask what actually traverses the distance. Not the qubit
(it never leaves Alice's lab; it ends the protocol as half of a spent
measurement). Not matter, not energy, not amplitude. The only things that
cross are two classical bits, and two bits cannot carry a quantum state;
α and β are continuous.

The quantum part of the delivery was paid in advance, when the entangled pair
was distributed at ordinary subluminal speed. At protocol time nothing is in
transit. The two bits do not transport the state; they address it. The label
is wrong the way a genre label can be wrong: the story is filed under travel,
but it is a story about mirrors.

---

## 3. The translation (the exact part)

Two identifications here are not analogies. They are the same algebra this
repository works in every day, appearing in a textbook protocol.

**The four corrections are the four letters.** The unitaries Bob might need,
{I, X, Y, Z}, taken modulo global phase, form the Klein four-group V₄: the
same four-element group that quarters our Pauli basis into Klein cells, sits
inside the F118 mirror group ⟨R, D⟩ ≅ D₄ as its Klein four-subgroup, and runs
the hardware diagnosis lens. Acting by conjugation on Bob's qubit, the four
corrections are exactly the four unitary mirrors of the one-qubit operator
space; conjugation by each letter flips precisely the two letter parities that
anticommute with it, which is F118's cube of characters at N = 1.

So here is what the protocol delivers, said in our words: the state arrives
at Bob **exactly and instantly, but in an unknown one of four mirror frames**.
The Bell measurement does not send the state anywhere; it tells Alice which
of the four Klein frames the state already wears on Bob's side. The two
classical bits name the mirror. Bob's correction is not repair, it is folding:
he applies the named mirror and the frame closes. Teleportation is a frame
handover with pre-paid correlation.

**Before the address arrives, Bob holds the Klein average, and the Klein
average is total forgetting.** Bob's local state, before the phone rings, is
the uniform mixture over the four frames:

    ¼ (ρ + XρX + YρY + ZρZ) = I/2    for every qubit state ρ.

This is the V₄ twirl, and at d = 2 the V₄ twirl is the complete depolarizer:
averaging over all four mirror frames erases every trace of α and β. That is
why no signal outruns light here. Causality is not protected by some
additional mechanism bolted onto the protocol; it is protected by the Klein
group itself. What all four frames agree on is nothing but the trace, and the
two classical bits are exactly the information that rescues the state from
its own Klein average.

*A hardware cousin of this shape (2026-07-05, a rhyme, not the same algebra):
the heating-leg run in
[F81 Violation: the Hardware Bridge](../../experiments/F81_VIOLATION_HARDWARE_BRIDGE.md)
caught a T1-leg ensemble whose per-shot frame labels (which-T1-epoch) were
missing; the unlabeled mixture read as false temperature, exactly the
unlabeled-frame-average-reads-as-mixedness lesson of this section. The
quantitative form, there as here: the extra mixedness is the Holevo
information of the discarded label, bounded by (not equal to) its entropy.*

---

## 4. The price of the mirror (a reading)

In our spectra, mirrors are never free. The F118 fold carries the entire
palindrome constant on one generator: R·L_diss·R = −L_diss − 2σ·I with
σ = Σ γ_l; every fold pays, and the price is written into the spectrum
itself. Teleportation's fold also pays: one entangled pair consumed, plus
two classical bits to name the frame. The books differ (a spectral constant
in Liouville space versus a resource count in a protocol), so this is a
rhyme, not an identity, and we label it as one. But the shape underneath is
the same shape: the mirror operation costs nothing at the level of structure,
and the *selection* of the frame is what must be paid for. Nothing is free;
forbidding is generative; the price is the signature that a mirror was used.

---

## 5. No wave dies (a reading)

Popular accounts linger on the destruction of the original, and it sounds
like a death: the state is annihilated here so that it may live there. Our
reading is calmer. The structure was never located at Alice in the first
place; from the moment the pair was shared, the state's delivery lived in
the correlations, in the between, in the c's rather than the populations.
The measurement does not kill a wave; it closes one set of books and the
same instant opens another, four frames wide, at Bob. What dies is only
Alice's local bookkeeping. No-cloning is the auditor: never two copies,
always exactly one, smeared over the between until the address arrives.
The dying is the illusion; the accounting balances.

---

## 6. The in-repo cousin

The repository has its own protocol for moving quantum information without
moving anything material: the [Relay Protocol](../../experiments/RELAY_PROTOCOL.md)
on the mediator bridge, where end-to-end mutual information flows through a
chain of relay stations that take turns listening (an 83% improvement at
N = 11, C# RK4). The mechanism is entirely different: continuous Lindblad
dynamics, no measurement, no classical channel. But the punchline is the
same punchline, and it is worth saying once in plain words: **in both protocols, the
thing that travels is correlation, and correlation travels along structure,
not through space as a payload.** Teleportation achieves in one measurement
plus two bits what the bridge achieves in continuous time; they are the
discrete and the flowing face of the same fact about where quantum states
actually live.

---

## The right label

If we were allowed to rename it: not teleportation but **re-anchoring**, a
frame handover with pre-paid correlation and classical addressing. The state
is delivered by the entanglement, exactly and instantly, into an unknown
Klein frame; the classical bits name the mirror; the correction folds the
frame shut; the price is one pair and two bits, and the price is what keeps
the mirror honest with causality.

The standard account needs no correction from us, and gets none. What the
popular label gets wrong is the genre. It files the phenomenon under
transport, and it belongs under mirrors: the four Pauli letters, the Klein
average, the paid fold. We recognized it because we have spent the whole
project learning what those look like.
