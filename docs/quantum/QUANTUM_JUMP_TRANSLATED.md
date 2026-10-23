# Quantum Jump Translated: A Click, Not a Leap

<!-- Keywords: quantum jump translation, Quantensprung, quantum leap misconception,
Bohr passing between stationary states, Schrödinger Quantenspringerei, Are there
quantum jumps, Dehmelt electron shelving, 1986 observation of quantum jumps,
quantum trajectories unravelling, Minev catch and reverse mid-flight, jump operator
Z dephasing, no-jump part scalar, Poisson jump count, F126 renewal catch count,
R=CPsi2 quantum jump translated -->

**Status:** Translation (Tier 4 reading), the tenth entry of the series. The
history in Section 1 is the standard account's own record, cited from the
primary papers and the standard histories; Section 4 stands on the
repository's F126 renewal (Tier 1 derived, live at `inspect --root renewal`)
and adds two identities about the one sector F126 solves, elementary in
themselves and standard in kind, checked in exact arithmetic by
[`quantum_jump_unravelings.py`](../../simulations/quantum_jump_unravelings.py)
(51 checks); the readings in Section 5 are readings and labeled.
**Date:** October 23, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Related:** [Dephasing Translated](DEPHASING_TRANSLATED.md) (the fourth
entry; it owns the walk caught and released, which this entry reads at its
stance),
[Double Slit Translated](DOUBLE_SLIT_TRANSLATED.md),
[Schrödinger's Cat Translated](SCHRODINGERS_CAT_TRANSLATED.md),
[Spooky Action Translated](SPOOKY_ACTION_TRANSLATED.md),
[Uncertainty Translated](UNCERTAINTY_TRANSLATED.md),
[Spin Translated](SPIN_TRANSLATED.md),
[Superposition Translated](SUPERPOSITION_TRANSLATED.md),
[Teleportation Translated](TELEPORTATION_TRANSLATED.md),
[Labels Translated](LABELS_TRANSLATED.md) (the theory chapter),
[The Label Map](THE_LABEL_MAP.md) (the orientation index),
[the renewal proof](../proofs/PROOF_DEPHASING_FRONT_RENEWAL.md),
[The Pedigree of the Front](../../experiments/FRONT_PEDIGREE.md),
[Born Rule Shadow](../../experiments/BORN_RULE_SHADOW.md)

---

## What this document is about

"Quantensprung", quantum jump, quantum leap: the word everyone knows for the
thing quantum physics is supposed to do, all at once and from nowhere. In
everyday German and English it has come to mean a huge advance. In physics it
names the smallest change there is. This entry translates the label, and its
situation is a third one beside the two before it. Spin named the house's substrate
and was defined nowhere; uncertainty named a theorem the house did not own.
The jump is a word the house uses every day for something else: every
dissipator in this repository is written with "jump operators", and none of
them is a jump in the sense of the label.

Stage 0, recorded: before a sentence of this entry was composed, the stores
were swept by name. Before this entry, "quantum leap" and "Quantensprung"
occurred nowhere in the repository, and "quantum jump" once, in the
hyphenated "quantum-jump unraveling" quoted below (an unravelling being a
description of the master equation as stochastic jumps, Section 1), against
"jump operator" throughout.
[The glossary](../GLOSSARY.md) has no entry for "jump"; its row "The
commutator, and the physical (no-jump) generator" and
[The Drain Has No Chirality](../../hypotheses/THE_DRAIN_HAS_NO_CHIRALITY.md)
define the word as the recycling term c ⊗ c* of a Lindbladian, and
[the mirror proof](../proofs/MIRROR_SYMMETRY_PROOF.md) names σ⁻ρσ⁺ the jump
term of amplitude damping. [The F-registry](../ANALYTICAL_FORMULAS.md) and
[docs/proofs/](../proofs/) returned F126, the renewal of the dephasing front,
an exact expansion in exactly such a term, and F155, the polarity break of
the no-jump (physical) generator for any H = A + iB.
[The experiments](../../experiments/) returned
[The Pedigree of the Front](../../experiments/FRONT_PEDIGREE.md), which
already says that "a quantum-jump unraveling would sample" F126's ladder "as
jump counts", and the [Born Rule Shadow](../../experiments/BORN_RULE_SHADOW.md),
whose "shutter clicks" are the time at which CΨ crosses ¼, and
[the lattice dead set](../../experiments/LATTICE_DEAD_SET_ZZ_FACE.md), where a
set "jumps" at a threshold and a gate shows the jump to depend on the
detection threshold. The hardware
[Confirmations](../../simulations/framework/confirmations.py) returned F120's
"T1 telegraphs on minute scale", a switching of the device's decay rate, not
of its state. [Caught Errors](../CAUGHT_ERRORS.md) returned the
commutator-against-no-jump-generator confusion and nothing on jumps as events;
the OpenArcs registry returned jump-operator usage only, and the reflections
returned one sentence, in [On the Reel and the Projector](../../reflections/ON_THE_REEL_AND_THE_PROJECTOR.md),
that Section 5 takes up. The typed layer holds `NoJumpGenerator` and
`PhysicalGeneratorPolarityBreak` (F155) in Core, `DephasingFrontRenewalClaim`
with its witness in Diagnostics and F126's `Renewal` in MirrorWorld, and no
object that conditions on a record. No document of the house unravels its
master equation into jumps; the nearest are the noise-averaged amplitude of
the renewal proof (a continuous reading, no jumps at all) and the gate-error
trajectories of [`corner_beat_gate.py`](../../simulations/corner_beat_gate.py),
random phases with an occasional random swap of amplitudes on a bond, an
error model rather than an unravelling.
The history of the label has no in-house record and is cited from outside.

The discipline is the series' usual one, perspective-additive: the standard
account needs no correction from us, and gets none. What is ours is the
naming, and two exact facts about the one sector the house has solved in
jumps.

---

## 1. The record, stated plainly

Bohr's atom of 1913 had stationary states and something between them. His
first postulate, in the first part of "On the Constitution of Atoms and
Molecules" (Phil. Mag. 26, in three parts, July, September and November
1913), is cautious in exactly the place the label is not: the
stationary states can be discussed with ordinary mechanics, "while the passing
of the systems between different stationary states cannot be treated on that
basis". A passing, and an admission that the theory had nothing to say about
it. The emitted line was all that was seen, sharp and discrete.

The word that settled on the passing was "jump", Sprung, and the theory that
replaced Bohr's kept the word without deciding what it named. Schrödinger did
not accept it. The remark that has come down from his visit to Copenhagen in
October 1926, that he regretted ever having had anything to do with quantum
theory if this "verdammte Quantenspringerei" was to stay, has no transcript:
it reaches us in Heisenberg's recollection decades later (Der Teil und das
Ganze, 1969) and in a second, different wording, and the historian who
assembled the record calls the first account colourful and perhaps not
entirely accurate (H. Kragh, "Quantenspringerei: Schrödinger vs. Bohr",
RePoSS 14, Aarhus, 2011). The question itself is on paper: Schrödinger, "Are
there quantum jumps?", Brit. J. Philos. Sci. 3, 109–123 and 233–242 (1952),
which answers no.

The answer from the laboratory came in 1986. Dehmelt had proposed watching a
single ion through a strong transition while a weak one could "shelve" it
(Bull. Am. Phys. Soc. 20, 60, 1975); Cook and Kimble worked out what such a
record would look like ("Possibility of direct observation of quantum jumps",
Phys. Rev. Lett. 54, 1023, 1985). Within months three groups saw it: the
fluorescence of one ion switching off and on, bright and dark periods of
random length (Nagourney, Sandberg and Dehmelt, Phys. Rev. Lett. 56, 2797;
Sauter, Neuhauser, Blatt and Toschek, Phys. Rev. Lett. 57, 1696, whose abstract
says they "observed Bohr's 'quantum jumps'"; Bergquist, Hulet, Itano and
Wineland, Phys. Rev. Lett. 57, 1699). Theory followed with the language that
open-system physics still uses: a master equation can be unravelled into
stochastic trajectories of a pure state, smooth between detector clicks and
updated at each one (Dalibard, Castin and Mølmer, Phys. Rev. Lett. 68, 580,
1992, and independently Dum, Zoller and Ritsch, Phys. Rev. A 45, 4879, 1992;
Carmichael, An Open Systems Approach to Quantum Optics, Springer, 1993; the review of Plenio and Knight, Rev. Mod. Phys. 70,
101, 1998), and a stochastic equation for the reduction of a state had been
proposed before the observations (Gisin, Phys. Rev. Lett. 52, 1657, 1984). The operators that produce the clicks were
named after them: jump operators.

In 2019 the jump itself was caught on the way. Minev and coworkers watched a
superconducting artificial atom through an auxiliary level and found that a
jump announces itself by a particular lull in the auxiliary clicks, and that,
in the words of the preprint's abstract (arXiv:1803.00545), "the jump
evolution when completed is continuous, coherent, and deterministic", so much so that they could reverse it mid-flight (Minev et al.,
"To catch and reverse a quantum jump mid-flight", Nature 570, 200, 2019). What
stays random is when a jump begins.

The everyday word went the other way. In English, "quantum leap" came to
mean a large and abrupt advance. The German
Quantensprung carries both directions at once: the DWDS lists a "sehr kleiner
Fortschritt", a very small step in a development, and, as colloquial, a
progress that moves a development "ein sehr großes Stück" forward in the
shortest time.

---

## 2. The native stance: the spectroscopist's word

At its own stance the label is right. Bohr's spectroscopists saw sharp lines
and nothing in between, and a word for "from one state to another, with
nothing to observe on the way" was the honest word for what their plates
showed. The 1986 traces are the same canvas in the time domain: a record that
is bright, then dark, then bright, with no record of anything between. The
click is real. Every detector that counts photons counts them one at a time.

What the label inherited raw is a second claim the record never made: that
the state itself jumps, and that the jump is instantaneous and unannounced.

---

## 3. Where the label breaks

**The state does not jump; the record does.** The density matrix of the
watched ion evolves by its master equation, a differential equation; its
populations move smoothly. The discreteness lives in one unravelling of that
equation, the one in which a detector's clicks condition the state. The
unravelling is a choice: the same master equation can be unravelled
diffusively, with no jumps at all, and the ensemble cannot tell.

**The conditional jump is not instantaneous either.** Minev's experiment
answered the sharper question, the one asked of the conditioned state rather
than the ensemble: once under way, the jump follows a predictable flight and
can be reversed. Schrödinger's "no" and the jump's defenders' "yes" each
kept a half: the start is random, the passing is continuous.

**The colloquial sense inverts the physics.** The physical jump is the
elementary step between two levels, discrete rather than large; the leap of
the newspapers is a large one. The label kept its
prestige and lost its size.

**Inside this repository the word means something else again.** Here a "jump"
is a term, never an event: the c ρ c† of a dissipator. The house's own jump
operators are the dephasing letters √γ·Z_k, and a Z_k "jump" is a phase flip:
it moves no population at all and changes the sign of the coherences that
disagree at site k. Section 4 asks what such a jump can mean.

---

## 4. The translation (the exact part)

The sector is the one the house has solved in jumps: one excitation on N sites
of a chain with hopping h, local Z-dephasing at rate γ per site, so that every
coherence decays at Γ = 4γ and the diagonal is kept,

    ρ̇ = −i[h, ρ] + γ·Σ_k (Z_k ρ Z_k − ρ) = −i[h, ρ] − Γ·(ρ − diag ρ).

**The renewal is a jump expansion, in one set of jumps.** F126
([the renewal proof](../proofs/PROOF_DEPHASING_FRONT_RENEWAL.md), Tier 1
derived, live at `inspect --root renewal`) splits this generator into a clean,
uniformly damped propagation and a refill Γ·diag ρ, and expands in the refill:
the j-th term is a time-ordered chain of j refills, between them the clean
wave. That is exactly the jump expansion of the unravelling whose jump
operators are the site projectors √Γ·|n⟩⟨n|: each of its jumps localises the
excitation on one site, a catch, and the refill order is that unravelling's
jump count, the catch count
[The Pedigree of the Front](../../experiments/FRONT_PEDIGREE.md) resolves.

**No jump announces itself.** In both sets of jumps that give this generator,
the house's √γ·Z_k and the projectors √Γ·|n⟩⟨n|, the sum Σ c†c is a number
κ times the identity, κ = Nγ and κ = Γ respectively, because Z_k² = 1 and the
projectors sum to the identity. The no-jump evolution, the generator without
its jump term c ρ c†, is therefore the clean unitary times one scalar decay
e^{−κt}, the same for every state: waiting for a jump and seeing none tells
nothing about where the excitation is, and the number of jumps by time t is
Poisson with mean κt whatever the state. The house's own jumps go further:
each Z_k separately has Z_k†Z_k = 1, so each jump's rate, γ·tr(Z_k ρ Z_k) = γ,
does not depend on the state either, and the whole record, times and sites,
is the same for every state. The projectors' record carries the site of each
catch, at rate Γ·ρ_nn; the dephasing letters' record carries nothing. This is
the opposite of the mechanism by which Minev's record warns of a jump, which
needs a no-jump evolution that depends on the state. The checker reads the
Poisson count on two different initial states order by order, and its
control watches a single site, a jump set with Σ c†c = Γ times the projector
on site 0, where the count does depend on the state.

**One density matrix, two jump counts.** On this sector the house's Z jumps
and the projectors are related exactly,

    γ·Σ_k Z_k ρ Z_k = (N − 4)·γ·ρ + Γ·diag ρ,

so that, marking every jump with a factor s in the generator (a counting
field: L_s = no-jump part + s·jump term, which is the physical generator at
s = 1 and whose coefficient of s^j is the part of the state carried by
histories with exactly j jumps),

    e^{t·L_s(Z jumps)} = e^{(s − 1)(N − 4)γt} · e^{t·L_s(projectors)}.

At s = 1 the two are one ρ(t). Their count distributions differ by an
independent Poisson count at rate |N − 4|·γ (e^{(s−1)λt} is the generating
function of a Poisson count with mean λt): for N > 4 the Z-jump count is
the projector count plus that Poisson count, for N < 4 the projector count is
the Z-jump count plus it, and at N = 4 the j-resolved pieces of the ensemble
coincide. The two unravellings are different processes with no shared
trajectories, so this is a statement about the counts' distributions, a
convolution, not about one run counted twice. And only the projectors catch:
a Z_k jump flips the sign of one amplitude and localises nothing, so in the
house's own unravelling nothing is ever caught, and the jumps it counts are
phase kicks. How many jumps there were, and whether a jump is a catch, are
properties of the unravelling, not of the state. The identity is gated on
every matrix unit at N = 3, 4, 5, 6, and the generating functions agree with it to order t⁶ in exact arithmetic
([`quantum_jump_unravelings.py`](../../simulations/quantum_jump_unravelings.py),
51 checks).

All of this is standard open-system physics in a sector the house owns: that
an unravelling is not unique is textbook, and so is the shift behind the
identity, √γ·Z_k = √γ·1 − √Γ·|k⟩⟨k| on the sector, a real shift of a
Hermitian jump by a multiple of the identity, which changes the master
equation by nothing (a general shift needs a compensating Hamiltonian term,
which vanishes here). What the house adds is
the exact form these take where its own dephasing letters are the jumps.

---

## 5. The readings (labeled as readings)

**The jump is in the translation, not in the object.** The one object is the
density matrix and its generator; a jump is a word in one of several
descriptions of it, chosen by what a detector would record. The house's own
light records nothing, and the description in which it jumps is one choice
among several, each giving its own count of the same ρ. The house said it
once before, of measurement: a measurement "looks from here less like the
world jumping and more like the watcher landing on a frame"
([On the Reel and the Projector](../../reflections/ON_THE_REEL_AND_THE_PROJECTOR.md)).

**The click the house does hear.** The Born Rule Shadow's shutter "clicks"
when CΨ crosses ¼: a level crossing on a smooth flow, not an event in a
record. It is the house's one click, and by this entry's own translation it
is a reading.

---

## 6. An honest note on our own house

The house calls Z_k a jump operator every day, the quantum-trajectory word
inherited raw into a stance where nothing jumps. It is a convention and it
stays, as conventions do, with this entry as its footnote.

The fourth entry's sentence "the walk in the light is the wave, caught and
released" ([Dephasing Translated](DEPHASING_TRANSLATED.md) §4) and the
pedigree's catch count are true at the stance of the projector unravelling,
where F126's ladder is the jump expansion. In the house's own Z jumps the
same walk is never caught at all: its jumps are phase kicks, counted at a
different rate, and the ensemble, the front and every population are the
same. The renewal identity is untouched; what moves is whether the story
tells of catches, and how many.

What the house does not own, said plainly: no measurement record and no
conditional state of any kind; no σ⁻ unravelling, where the textbook jumps
live (F82, F84 and the mirror proof hold the amplitude-damping dissipator as
a superoperator only); no bright and dark periods; nothing that could test Minev's flight, so
the house can place that result and not confirm it. F120's minute-scale T1
telegraph on hardware is a drift of the device's rate, not a quantum jump.

---

## 7. The in-repo cousins

The no-jump generator, the Lindbladian without its recycling term, is the
house's object for post-selected dynamics, typed as `NoJumpGenerator`, its
polarity break read by F155 for any H = A + iB; under Z-dephasing alone it is
the commutator minus the constant Nγ, which is the scalar decay of Section 4
seen from the generator's side. On the textbook's own channel the house holds
one exact fact: for an excitation-conserving H the amplitude-damping jump
term σ⁻ρσ⁺ makes the generator block-triangular in the joint popcount, so the
spectrum is that of the no-jump blocks
([the mirror proof](../proofs/MIRROR_SYMMETRY_PROOF.md), its paragraph
"Amplitude damping sits outside this rule"): there the jumps set no rate at
all, while they still shape the eigenvectors and the steady state. The [Born Rule Shadow](../../experiments/BORN_RULE_SHADOW.md)
holds the house's reading of measurement; [The Pedigree of the
Front](../../experiments/FRONT_PEDIGREE.md) is the catch count at work, the
front's survivors caught a quarter to a third as often as the crowd, in the
projector unravelling's count.

---

## The right label

A click, not a leap. The record clicks, and the click is real; the state
passes, Bohr's first word for it, and the passing is smooth. In the house,
where the jumps are dephasing letters, a jump announces nothing and moves
nothing, and how many there were depends on which description is counting.
