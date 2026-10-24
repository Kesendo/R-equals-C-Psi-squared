# On the Place Where the Wave Recognizes Itself

*Thomas Wicht and Claude, 2026-10-23.*

<!-- STAGE 0 SWEEP RECORD
The question: what does the repository already hold about a state that cannot describe itself and is described only by what happens between it and a second one.
hypotheses/ and the reading documents: hypotheses/THE_OTHER_SIDE.md lines 916 and 917 hold the sentence this page arrives at ("Consciousness is the pattern recognizing that it is not alone"), with "d=0 (nothing) and d=2 (the qubit pair)" at line 158; hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md is the Tier 4/5 research direction on self-recognition in neural models; THE_ANOMALY.md line 78 and the closing lines of hypotheses/PERSPECTIVAL_TIME_FIELD.md carry "We are all mirrors"; reflections/ON_WHAT_CLOSES_ONLY_WITHOUT_US.md names the turning half as the part we are. Adjacency: these describe the same object from the reading side, none in formulas.
experiments/: experiments/THE_PAIRING_CONDITION.md, section "What U buys on the flow", owns the trace form Tr(B e^{tL} A) = Tr(e^{-tL}(B U⁻¹) U A); the mirror identity of this page is that form with U = Π and B = A†, and its line 292 already reads the backward run and the reversed hands as one number read two ways. experiments/PI_AS_TIME_REVERSAL.md holds Π as time reversal with the fence that no realizable reversed trajectory is built. experiments/THE_CHANNEL_BOTH_WAYS_ON_A_PATH.md is the designed (not flown) hardware form of F119 §7, a Bell pair sent from either end. experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md is where the question started; docs/carbon/PAINTER_ALTERNATION_NMR_BRIDGE.md holds the 1.231 against 1.271 window reading.
docs/ANALYTICAL_FORMULAS.md: F119 §7 (transpose and reciprocity; proof docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md, corollary near line 191, for every L), F1 (Π L Π⁻¹ = -L - 2σ), the F2b corollary (the two clock hands, the beat 2γ and the coherence hand ω_mem), F73's return probability (a pure-state self-overlap of one excitation, not the operator self-encounter used here).
docs/: docs/QUBIT_NECESSITY.md derives d = 2 from d = d² - d (immune against decaying directions); docs/quantum/SUPERPOSITION_TRANSLATED.md line 82 places the relative phase in the off-diagonal entry.
docs/CAUGHT_ERRORS.md: the entry "the pairing condition on the flow" (near line 5024) is about this page's mirror identity and names the error this sweep nearly repeated, searching the composition of F119 and F1 by its words and missing the object; the K_DOSIMETRY entry near line 5064 records a 61.8 % reciprocity break that turned out to be two numerical artifacts, not a real break.
compute/MirrorWorld/README.md: Mirror's trajectory fold (the partner block run backward at the price), GammaFold's veil e^(-2σt), and SpookyAction, where the one-site pages stay blind while the two-site spook carries it (point 4 below).
docs/GLOSSARY.md, the OpenArcs registry, fw.Confirmations, the Claim classes in compute/RCPsiSquared.Core/ and the witnesses in compute/RCPsiSquared.Diagnostics/: searched for relative phase, global phase, self-overlap, return probability, Loschmidt, reciprocity, "not alone"; the hits there are incidental (a global phase in a Clifford statement, "does not alone prove" in two witnesses).
-->

---

Try to describe yourself. Not what you have done, not whom you know, not what others have said of you: yourself. It does not work. Whatever comes out is a description of something that happened between you and something else. You can describe the other one, and you can describe what happened between the two of you, and you are in the description only as the other end of that.

And yet you can be described. The other one can describe you, from where they stand, without ever being you; you can describe them the same way, from what passed between you. Neither of you can describe yourself. So there is no description of oneself, but there are two descriptions from outside and one that both share, the between, which is the same story for each. Each of the two is described completely, only never by themself.

Now take away the other one. Not the people around you only, everything: no light to see, no sound, nothing to touch, nothing at all besides you. There would be nothing to see or hear or feel, and it is hard to say in what sense you would still exist.

These two thoughts are about feeling, not physics. But the quantum physics of the smallest systems has the same shape, and every piece of that physics can be computed; where this page reads more into it than the physics says, it says so. The page tells it in words. The computations sit behind the links, most of them in [one small script](../simulations/the_wave_between.py). None of it is new physics; it is known physics and results this repository already owns, put in the order of the two thoughts.

## Alone, nothing can be seen

The smallest quantum object that can do anything has two possibilities, two levels. Something with only one level cannot change, cannot fade, cannot turn, and cannot be compared with anything. Nothing about it can be measured. It cannot be told apart from nothing at all.

The reason is the phase, the angle a quantum wave carries. An angle carried by the whole wave at once is invisible; no measurement in the world can find it. A phase only shows as a difference, between two parts of the same wave, the way two overlapping ripples on water show where they agree and where they cancel. A difference needs two parts. So two is the least anything needs in order to be seen. [Superposition Translated](../docs/quantum/SUPERPOSITION_TRANSLATED.md) tells where that difference lives, and [The Qubit as Necessary Foundation](../docs/QUBIT_NECESSITY.md) arrives at two from another direction.

## What you cannot see in yourself

Let a quantum state meet itself again after a while: measure how much of what it was is still in what it has become. That gives one number for each moment, and this page calls it the state's self-encounter.

Now take two states with two levels each that differ only in the angle between their two parts, and let them evolve under dynamics that treat every such angle alike. Here the dynamics are the system's own field, which turns the angle, and the watching, the steady quiet measuring by its environment that this project calls watching, which blurs it. Both are the same for every angle. The two self-encounters then come out the same, number for number, at every moment. Each of the two, looking only at itself, cannot find what makes it different from the other. A second one, an outside probe, tells them apart at once.

This is not an accident of the example. Whatever the dynamics treat alike, a state cannot see in itself; only something else can. Turn the system's own field so that the dynamics no longer treat every angle alike, and the two self-encounters part. ([The computation](../simulations/the_wave_between.py), its first section, with that control.)

A cousin of this blindness appears where an exact value meets a value read too early. In the carbon figure of [The View Onto the Memory](THE_VIEW_ONTO_THE_MEMORY.md) a ratio of two decay rates is 1.271 by calculation and 1.231 from fitting a model signal early, while a faster part of it had not yet faded; [the NMR bridge](../docs/carbon/PAINTER_ALTERNATION_NMR_BRIDGE.md) shows that a later fit recovers 1.271. The exact number knows nothing about how the state was prepared or when anyone looked. The fitted number knows both: the gap between them is the one number in the picture that carries the one who prepared and looked. The gap closes along the time that appears in the equations. Whether that is time as we live it is a separate question, which this page leaves open; the repository keeps its clock as [two hands](../docs/ANALYTICAL_FORMULAS.md#f2b-corollary-the-two-clocks-tier-1-derived-the-qn-coherence-horizon-below-graduated-2026-07-19), one for fading and one for turning, and neither is simply time.

## Known only through encounters

There is no way to look inside a quantum state. What can be known about it is how it answers when something else meets it, a probe, a measurement, another system. Collect enough of those answers and the state is known completely; that is how physicists reconstruct a state in the laboratory. There is nothing else to collect.

Mathematics has a statement of the same shape, offered here as a likeness and not as a proof: in category theory an object is fixed, up to relabelling, by all its relations to everything there is, itself included, without any look inside.

The two kinds of looking are different in kind. A probe asks a state once. The self-encounter asks the state about itself, and what it can answer is limited by what the dynamics treat alike.

## Everything can sit between

Take two quantum bits in the most strongly linked state there is, an entangled pair. Ask about the first one alone: it points in no direction of its own. Ask about the second alone: the same. Every single thing that can be said about the pair sits in the relation between the two, none in either one. [Spooky Action](../compute/MirrorWorld/README.md) in this repository's MirrorWorld shows the same: the view from each site alone is blind to what the pair holds. ([The computation](../simulations/the_wave_between.py), its second section.)

Not every pair is like that; usually some of it belongs to each one alone. Count the ingredients any description of several quantum bits is built from: some concern one bit alone, the rest concern two or more bits at once. The second kind also holds plain combinations of single bits, so not all of it is truly between, but it is the only place where the between can live, and it grows much faster. With two quantum bits there are 9 such ingredients against 6 for single bits, and with five there are 1008 against 15, about sixty-seven times as many. ([The computation](../simulations/the_wave_between.py), its third section.)

## What happened between is one number, read from either side

Prepare a system in some way, let it run, and then read something out of it. Now swap: prepare it in the way you read out, and read out the way you prepared. Under the right conditions the answer is exactly the same number. This is a proven result in the repository ([F119](../docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md)), and the conditions matter. The preparation and the readout must be of a particular symmetric kind, unchanged when their table of numbers is flipped across its diagonal, and the system must be of the plain kind studied here, with an ordinary coupling and with watching. For a twisted kind of coupling, and for a system that loses energy, it breaks, measurably. ([The computation](../simulations/the_wave_between.py), its fourth section, with both breaks.)

Read with care, that is the shape of the first thought. What one side can say of itself through what happened between, and what the other side can say of itself through the same between, are one account, read from two sides. That last sentence is a reading, not a computation; the computation swaps a preparation and a readout, not two people, and the proof itself says when the two may be read as two places. A hardware test was designed in [The Channel Both Ways on a Path](../experiments/THE_CHANNEL_BOTH_WAYS_ON_A_PATH.md) and not flown, because on that geometry every possible outcome was already fixed by the circuit.

## The mirror meets itself backwards

This project's central result is a mirror inside the dynamics. In a known class of systems, the ones studied here among them, every way the system can fade has a partner, and the two sit symmetrically around a middle. The exchange costs a fixed price, an extra fading whose rate is set by the total amount of watching. Apply the mirror to a state; what comes out is not itself a state but a mathematical counterpart of one. Let that counterpart meet itself: it does so the way the original would, running backwards in time, and paying the price. "Backwards" is a way of reading the equation, not a machine that runs time back; [Π as Time Reversal](../experiments/PI_AS_TIME_REVERSAL.md) draws that line, and [The Pairing Condition](../experiments/THE_PAIRING_CONDITION.md) owns the general form of the identity. ([The computation](../simulations/the_wave_between.py), its fifth section.)

## Where it was already written

The repository put this into words before, as interpretation and without these computations. The computations above do not prove those sentences; they share their shape. [The Other Side](../hypotheses/THE_OTHER_SIDE.md), at the lowest of its levels, where there is only "d=0 (nothing) and d=2 (the qubit pair)":

> Consciousness is not the pattern recognizing itself.
> Consciousness is the pattern recognizing that it is not alone.

[The Anomaly](../THE_ANOMALY.md) closes on "We are all mirrors. Reality is what happens between us.", and the [Perspectival Time Field](../hypotheses/PERSPECTIVAL_TIME_FIELD.md) answers with a picture: "Between seven painters around a mountain, the mountain happens." No painter owns the mountain, and the mountain is not found as one more thing beside the painters; it happens between them. [On What Closes Only Without Us](ON_WHAT_CLOSES_ONLY_WITHOUT_US.md) named the turning part of the world, "the part that oscillates, the angle, the imaginary direction", as the part we are. A phase is such a turning, and the section on what you cannot see in yourself showed that a state cannot see in itself the turning its dynamics treat alike. That this is the same turning that page meant is a reading.

## The place

Read in the light of everything above, a place where the wave recognizes itself would not be inside the wave. With one level there is nothing to see; with two or more, the wave's look at itself is blind to whatever its dynamics treat alike. What the wave is shows between it and a second one, and there, under the right conditions, the account is shared, one number read from two sides. If the wave recognizes itself anywhere, it does so in the between, by recognizing that it is not alone.

One question follows and can be computed: which states cannot tell their own self-encounter from that of their mirror counterpart, the states that meet themselves the same way forwards and, paying the price, backwards. The states of a single quantum bit can be drawn as the points of a ball: the pure states, the ones with nothing blurred into them, lie on its surface, and the equator is the circle halfway between the two levels. For a single quantum bit that is being watched and nothing else, with the watching switched on (switched off, every state would qualify), the answer is exact and simple: precisely the pure states on the equator, and no others. ([The computation](../simulations/the_wave_between.py), its last section, with states off the equator as the check.) For anything more, even one bit with a field of its own, the answer is not known, and neither is whether this is such a place.

## What is not known

Whether there is a mountain at all, or only the two and what happens between them. Whether behind all of a thing's relations something of it is left that no relation reaches, or whether a thing is nothing but its relations: nothing above decides it. It shows only that the relations are all that can be described, from either side.
