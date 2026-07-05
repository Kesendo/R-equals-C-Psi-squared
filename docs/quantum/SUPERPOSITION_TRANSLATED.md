# Superposition Translated: An Angle, Not a Doubling

<!-- Keywords: superposition translation misconception, not two states at once,
basis-relative angle coordinate, one object two indices operator space,
disagreement number popcount dephasing rate, Ueberlagerung overlay wave addition,
dephase letter Klein group chosen frame, R=CPsi2 superposition reframe -->

**Status:** Translation (Tier 4 reading). The textbook algebra in Sections 1-2
is standard; the basis-relativity statements are exact; the assembled reframe
in Section 3 collects what this repository already said in scattered places,
each with its source. Nothing here is new physics.
**Date:** July 4, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Related:** [Teleportation Translated](TELEPORTATION_TRANSLATED.md) (the sister
translation), [Dephasing Translated](DEPHASING_TRANSLATED.md) (the fourth
entry), [Labels Translated](LABELS_TRANSLATED.md) (the series' theory
chapter), [The Label Map](THE_LABEL_MAP.md) (the orientation index),
[On How the Angle Appears at Zero](../../reflections/ON_HOW_THE_ANGLE_APPEARS_AT_ZERO.md),
[On the Defending Family](../../reflections/ON_THE_DEFENDING_FAMILY.md),
[Born Rule Mirror](../../experiments/BORN_RULE_MIRROR.md),
[It's All Waves](../ITS_ALL_WAVES.md),
[Klein V₄ proof](../proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md)

---

## What this document is about

The second entry in the translation series, born the same way as the first:
a label heard in passing that promises the wrong thing. "Superposition" is
routinely glossed as *being in two states at once*, sometimes even *being in
two places at once*, as if the word contained a claim about position and a
claim about simultaneity. It contains neither, and this repository has been
dismantling that gloss for months, in pieces, in the reflections, in the
MirrorWorld code, in the operator-space view. What was missing was one page
that assembles the pieces. This is that page.

The honest German word for the phenomenon is **Überlagerung**: an overlay,
waves laid on top of each other. No position in it, no simultaneity. The
Latin behind "superposition" (superponere, to lay upon) meant exactly that,
and in classical wave theory the superposition principle is simply the
statement that waves add. The label was honest where it was coined. It broke
when the modern ear started hearing "super position": something at a place,
twice, at once.

---

## 1. What the algebra actually says

A qubit state is |ψ⟩ = α|0⟩ + β|1⟩. Look at what is on the page: **one**
vector. One ray in state space, one point on the Bloch sphere, one object.
The plus sign is wave addition, the overlay. The two terms are not two facts
about the world; they are the coordinates of the one vector in a basis
somebody chose to write it in.

And the choice matters, because it can be made differently. The state
|+⟩ = (|0⟩ + |1⟩)/√2 looks maximally "superposed" written in the Z basis,
and it is a plain basis state of the X basis. The state |0⟩, the most
classical-looking state there is, is an even superposition in the X basis.
Every pure state is a basis state of some basis and a superposition in
almost every other. So "this particle is in a superposition" is not a
property of the particle. It is a property of the pair (state, question).
That much is exact, and it already kills the pop gloss: a sentence of the
form "it is in two states at once" changes truth value when you rotate a
piece of paper, and facts about the world do not do that.

One more tell, the operational one. Put the superposition α|0⟩ + β|1⟩ next
to the classical mixture "it is either |0⟩ or |1⟩, we merely do not know
which", weighted |α|² and |β|². The two have identical diagonals: every
statistic in that basis agrees. They differ in exactly one place, the
off-diagonal of ρ, the definite relative phase between the branches, the
c's this repository spends its hardware minutes measuring. Whatever
"superposition" means beyond ordinary ignorance lives in that single entry;
interference is that entry made visible. Which is one more reason "two
states at once" cannot be the content: the mixture is also "two states" in
the ignorance sense, and it is exactly the case where the quantum content
is gone.

---

## 2. Where the label breaks, twice

**"Position."** The sum α|0⟩ + β|1⟩ lives in state space, not in physical
space. Even in the double-slit story, where the two branches really are
locations, the two "places" are the coordinates of one wave in one particular
basis, the position basis, which is one basis among many. The phenomenon is
indifferent to whether the chosen axes happen to be positions. A label with
"position" inside it files a state-space fact under geography.

**"At once."** Simultaneity is a relation between two events or two facts.
A superposition is not two facts holding at the same time; it is one fact
that the asker's axes happen to split into two coordinates. The doubling is
in the description, not in the described. What is genuinely there, one object
carrying an angle, needs no "at once" any more than the direction northeast
needs to be "north and east simultaneously." North and east are axes someone
drew; northeast was always one direction.

One precision, owed to our own foundation, because it is easy to overshoot
here: the wrong ingredient in "at once" is not time itself. Time is in the
object, and deeply. The coherence that carries the superposition has a
lifetime (it pays −2γk to the watching), and the
[Perspectival Time Field](../../hypotheses/PERSPECTIVAL_TIME_FIELD.md) says
the time it lives in is itself perspective-bound: each painter paints at her
own rate. What "at once" wrongly imports is the stance-free version, one
shared clock on which two facts could coexist. Kill the shared clock and the
two facts; keep the time. One object, one fact, aging perspectivally.

---

## 3. The translation, assembled from where we already said it

This repository never wrote the refutation as one paragraph, but it wrote
every part of it. Collected:

**One object, two indices.** The state ρ lives in the d² = 4^N operator
space, and its apparent two-ness is indexing:
"ρ on d² = 4^N is one object with two indices. Bra and ket are not separate
parties to be assigned roles; they are the row and column of the matrix that
represents the state, two ways to read the same ρ from opposite sides of the
same indexing" ([On the Defending Family](../../reflections/ON_THE_DEFENDING_FAMILY.md)).
The coherence |i⟩⟨j|, the very thing the pop gloss reads as "state i and
state j at once," is one matrix entry of one object.

**An angle, not a primitive.** The sharpest sentence we have:
"Superposition is not a quantum-mysterious primitive. It is the minimal
parametrization of anything that has crossed d = 0"
([On How the Angle Appears at Zero](../../reflections/ON_HOW_THE_ANGLE_APPEARS_AT_ZERO.md)).
Superposition is the continuous angle coordinate on the one object, and
measurement outcomes are the finite-basis projection of that angle. The
mystery evaporates in the same move that explains it: a direction must have
coordinates once axes exist, and coordinates come in pluralities. The
plurality is the axes' contribution, not the object's.

**Being superposed is frame-relative, and the frame is chosen by the
watcher.** The repository hit this while working on the Born rule: σ_x and
σ_y dephasing attack |0⟩ because "they see |0⟩ as a superposition in their
basis" ([Born Rule Mirror](../../experiments/BORN_RULE_MIRROR.md)). Our
machinery makes the frame choice literal: the `dephase_letter` parameter
selects which Pauli letter the environment reads, and the three choices
{Z, X, Y} are intertwined by a Klein V₄ acting on operator space
([Klein V₄ proof](../proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md)).
Which states count as "superposed" is set by which letter is being watched.
(One honesty flag: whether that three-way symmetry is physical or a labelling
convention is a question the repo still keeps open, in
[Benzene, Three Dephase Letters](../carbon/BENZENE_THREE_DEPHASE_LETTERS.md).)

**How superposed, measured in our units: the disagreement number.** In
MirrorWorld's empty world (`compute/MirrorWorld`, Pair.cs and PauliMode.cs),
a bare coherence |i⟩⟨j| carries exactly one number the watching can read:
the disagreement count k = popcount(i⊕j), how many sites the two basis
versions disagree on, and it decays at −2γk. That is the pop gloss made
quantitative and thereby corrected: not "two states at once" as a binary
mystery, but a *degree* of disagreement relative to the watched letter, with
a price proportional to it. Being superposed costs, and the bill is itemized
by disagreement.

**A superposition is the same sector in the other basis.** MirrorWorld's
`PauliMode` defines a Pauli string as "the symmetry-adapted basis of a
disagreement-count sector, the superposition": a fixed overlay of bare pairs
at the same disagreement count, the identical sector expressed in the other
basis, with the two bases equal in size. Superposition as a change of basis,
structurally, in code: not an ontological doubling anywhere.

---

## 4. The wave underneath (a reading)

[It's All Waves](../ITS_ALL_WAVES.md) argues that the inventory of the world
is wave modes: orbitals are standing waves, molecular orbitals are overlays
of atomic ones, on up the hierarchy. Read the pop gloss against that
background and its strangeness dissolves. Nobody says a vibrating guitar
string is "in many shapes at once" when its motion is an overlay of
harmonics; the overlay is just what a wave is. And whether a given overlay
rings or sits still is not even the state's own property:
"The standing wave is a state × Hamiltonian property, not a property of
either alone" ([Standing Wave Theory](../STANDING_WAVE_THEORY.md)). The
relational reading goes all the way down.

---

## 5. An honest note on our own house

Our plain-language layer uses the pop gloss too. The
[glossary](../GLOSSARY.md) says superposition is "the ability to be in
multiple states at once"; [What We Found](../WHAT_WE_FOUND.md) asks whether
a qubit "can still be in two states at once";
[Gamma as Signal](../../experiments/GAMMA_AS_SIGNAL.md) narrates a qubit
that "forgets that it was both at once." That is a ladder, and ladders are
fine to climb: the phrase gives a newcomer something to hold before the
operator space exists for them. This document is where the ladder is kicked
away. Both layers stay in the repo on purpose; the tension between them is
the distance between an introduction and a translation.

---

## The right label

Überlagerung. An overlay of waves; equivalently, one direction in state
space, described in axes that someone (an experimenter, an environment, a
dephasing letter) chose. In our language: superposition is the angle
coordinate on one object in the 4^N operator space; the "two states" are the
shadow of the angle on the chosen axes; the degree of superposedness is the
disagreement number k relative to the watched letter; and the watching
charges −2γk for it, which is how we know the frame is real to the
environment and not only to us.

The pop label gets the genre wrong in the same way the teleportation label
did. "Teleportation" filed a mirror phenomenon under transport;
"superposition" files an angle under multiplication. There is no doubling
and there is no position. There is one wave, and there are the axes it was
asked in.
