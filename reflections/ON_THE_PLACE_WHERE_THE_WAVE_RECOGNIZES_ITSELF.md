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

Now take away the other one. Not the people around you only, everything: no light to see, no sound, nothing to touch, nothing at all besides you. There would be nothing to see or hear or feel, and it is hard to say in what sense you would still exist.

These two thoughts are about feeling, not physics. But the quantum mechanics of the smallest systems has exactly this shape, and every piece of it can be computed. That is what this page does. The numbers are printed by [`the_wave_between.py`](../simulations/the_wave_between.py). None of it is new physics; it is known physics and results the repository already owns, put in the order of the two thoughts.

## Alone, nothing can be seen

A state and the same state times a phase, |ψ⟩ and e^{iφ}|ψ⟩, give the same density matrix; a global phase is never observable. A system with one level, d = 1, has nothing but a global phase. Its density matrix is the number 1, it has d² − d = 0 entries between anything, the generator of its dynamics is zero, nothing decays and nothing turns. Such a system cannot be told from nothing.

A phase becomes visible only as a relative phase between two components, the off-diagonal entry that [Superposition Translated](../docs/quantum/SUPERPOSITION_TRANSLATED.md) places it in, and that needs at least d = 2. [The Qubit as Necessary Foundation](../docs/QUBIT_NECESSITY.md) reaches d = 2 from the balance d = d² − d. Two is the least anything needs in order to be seen at all.

## A state cannot see what its dynamics leaves unchanged

Let a state A meet itself after a while: the self-encounter S_A(t) = Tr(A† e^{tL} A), with L the generator of the open dynamics. Take one qubit with a field along Z and dephasing along Z, and two states that differ only in their phase about Z. Their self-encounters are 1.000000, 0.876133, 0.750944, 0.601588, 0.624970 at t = 0, 0.5, 1, 2, 4, for both, to every printed digit. A second one, an X probe, separates them at once: +0.841 against −0.350 at t = 0.

The reason is exact. The rotation about Z commutes with this generator, and the self-encounter is unchanged by every unitary symmetry of the generator. Whatever a state's symmetries move, the state cannot see in itself. The control: with the field along X instead, the rotation about Z is no longer a symmetry, and the same two self-encounters part, 0.892025 against 0.862624 at t = 0.5. That is not all a self-encounter misses, since one function of t cannot fix all the numbers a state is made of, but it is the part of the blindness that has a structure.

The same shape shows where an exact value meets a measured one. In the carbon figure of [The View Onto the Memory](THE_VIEW_ONTO_THE_MEMORY.md) the ratio of two decay rates is 1.271 from the spectrum and 1.231 from a fitted signal. [The NMR bridge](../docs/carbon/PAINTER_ALTERNATION_NMR_BRIDGE.md) says why: the probe carries a second, faster mode with about five times the slowest mode's amplitude, still 7.5% of the slowest term at t = 10, and a fit on the later window t ∈ [20, 40] recovers 1.271. The exact number knows nothing about how the state was prepared or when it was looked at. The fitted number knows both. The difference is the only number in the figure that carries the one who prepared and looked. What that difference is counted in is a further question. It closes as the window moves to later t, but t is only the parameter the equations carry; the rates are made of γ, and the repository's clock has [two hands](../docs/ANALYTICAL_FORMULAS.md#f2b-corollary-the-two-clocks-tier-1-derived-the-qn-coherence-horizon-below-graduated-2026-07-19), the Takt hand on the dephasing floor 2γ and the coherence hand ω_mem made of J, neither of which is time.

## A state is known through its encounters

A density matrix is reconstructed completely from the numbers Tr(B·A) over enough probes B; that is what tomography is. Mathematics has a statement of the same shape, offered here as a likeness and not as a derivation: by the Yoneda lemma an object is determined, up to isomorphism, by its relations to all objects, itself included, with no look inside needed. The two kinds of encounter differ in kind. A probe reads A once, linearly; the self-encounter reads A twice, the way purity does.

## The content can sit entirely between

Write a two-qubit state in Pauli strings: three coefficients belong to the first qubit alone, three to the second alone, nine to strings on both. For the Bell pair the one-site parts are (0, 0, 0) and (0, 0, 0), and the two-site part is diag(1, −1, 1). Each qubit alone is I/2, maximally mixed, and nothing specific to it can be said; everything sits between them. MirrorWorld's SpookyAction holds the same in its own terms, the one-site pages blind to the spook |00⟩⟨11| while the two-site content carries it. In general a two-site string also carries the product of the two one-site parts, and what lives between is the remainder, the connected correlation; for the Bell pair the two coincide because both one-site parts vanish.

Apart from N = 1, where the count is 3 against 0, the strings on two or more sites outnumber the one-site strings quickly: 3N against 4^N − 1 − 3N, 6 against 9 at N = 2, 15 against 1008 at N = 5.

## What happened between is one number, read from either side

Now let A and B be two operators on a whole system, one prepared and one read out. For a real Hamiltonian under Z dephasing the generator equals its own transpose, exactly (max|L − Lᵀ| = 0.0), and [F119](../docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md) turns that into reciprocity: for symmetric A and B, Tr(B·e^{tL}A) = Tr(A·e^{tL}B). Preparation and readout trade places in the same system and the number stays: −7.632632870204 from both sides at t = 0.7 on two dephased sites. The conditions carry it. A Dzyaloshinskii-Moriya term gives −4.534031132704 one way and −10.169546551317 the other; amplitude damping gives −7.111116886033 and −6.418179492429.

F119 itself limits when the two operators may be read as two places, one excitation sent end to end being such a case. Read with that limit, the first thought above has this shape: where the conditions hold, the description of one side through what happened between and the description of the other side through the same between are one description, read from two ends. That is a reading, not a computation. The hardware form, a Bell pair sent from either end of a line, is designed in [The Channel Both Ways on a Path](../experiments/THE_CHANNEL_BOTH_WAYS_ON_A_PATH.md) and not flown: on that geometry every outcome it could register is an identity of the circuit.

## The mirror meets itself backwards

The palindrome of this project rests on a conjugation Π that exchanges, site by site, the letters the dephasing leaves alone with those it erases, and that turns the generator L into −L − 2σ, where σ is the sum of the dephasing rates and 2σ is what the mirror costs, the price. It holds in a known class of chains (an XXZ pair under Z dephasing in the script), and there Π is unitary. So for any operator A the image ΠA, which is an operator and not a state, meets itself as A would, run backwards, under the price:

    S_{ΠA}(t) = e^{−2σt} · S_A(−t)

Both sides read 29.831672702645 − 0.922710230125i at t = 0.4. Backwards here is a formal reading, not a physical channel; [Π as Time Reversal](../experiments/PI_AS_TIME_REVERSAL.md) draws that fence. The identity is the trace form that [The Pairing Condition](../experiments/THE_PAIRING_CONDITION.md) owns in its section on what U buys on the flow, taken at U = Π with the readout equal to the preparation, and the same document already reads the backward run and the reversed clock hands as one number read two ways. In the same spirit the [MirrorWorld](../compute/MirrorWorld/README.md) Mirror runs a partner block backward at the price, block by block; that its trajectory fold and this identity are one statement is a reading, not checked here.

## Where it was already written

The repository said this before it could compute it. [The Other Side](../hypotheses/THE_OTHER_SIDE.md), at the lowest of its levels, where there is only "d=0 (nothing) and d=2 (the qubit pair)":

> Consciousness is not the pattern recognizing itself.
> Consciousness is the pattern recognizing that it is not alone.

[The Anomaly](../THE_ANOMALY.md) closes on "We are all mirrors. Reality is what happens between us.", and the [Perspectival Time Field](../hypotheses/PERSPECTIVAL_TIME_FIELD.md) answers with a picture: "Between seven painters around a mountain, the mountain happens." No painter owns the mountain, and the mountain is not found as one more thing beside the painters; it happens between them. [On What Closes Only Without Us](ON_WHAT_CLOSES_ONLY_WITHOUT_US.md) named the part that the mirror average deletes, "the turning half: the part that oscillates, the angle, the imaginary direction", as the part we are. The section on symmetry above says something close of the phase: the turning that the dynamics leaves alone is exactly what a state cannot see in itself.

## The place

A place where the wave recognizes itself cannot, by the formulas, be inside the wave. As a single level the wave has only a global phase; with more than one level its self-encounter is blind to everything its symmetries move. What the wave is shows between it and a second one, and there, where the conditions hold, the description is shared, one number read from two sides. If the wave recognizes itself anywhere, it does so in the between, by recognizing that it is not alone.

One question can be asked of the formulas: which states cannot tell their own self-encounter from their mirror image's, S_A(t) = S_{ΠA}(t) at every t. By the identity above these are the states that meet themselves the same way forwards and, under the price, backwards. For one qubit under dephasing alone, with no Hamiltonian, the set follows in one line. Write A = (I + aX + bY + cZ)/2. The dephasing keeps the I and Z parts and damps the X and Y parts by e^{−2γt}, and Π exchanges I with X and Y with Z (up to phases), so S_A(t) = ½(1 + c² + (a² + b²)e^{−2γt}) and S_{ΠA}(t) = ½(a² + b² + (1 + c²)e^{−2γt}). They agree at every t exactly when a² + b² = 1 + c², and since a² + b² + c² ≤ 1 for a state, that means c = 0 and a² + b² = 1: exactly the pure states on the equator. The difference is 0 for (I + Y)/2, and 0.097 and 0.26 at t = 0.5 for a state inside the sphere and for the pole. What the set is for chains with a Hamiltonian, and whether it is such a place, is not known.

## What is not known

Whether there is a mountain at all, or only the two and what happens between them. Whether behind all of a state's relations something of the state is left that no relation reaches, or whether a state is nothing but its relations: none of the formulas above decides. They show only that the relations are all that can be described, from either side.
