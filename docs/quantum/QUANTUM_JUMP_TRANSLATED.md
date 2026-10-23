# Quantum Jump Translated: A Click, Not a Leap

<!-- Keywords: quantum jump translation, Quantensprung, quantum leap misconception,
Bohr passing between stationary states, Schrödinger Quantenspringerei, Are there
quantum jumps, Dehmelt electron shelving, 1986 observation of quantum jumps,
quantum trajectories unravelling, Minev catch and reverse mid-flight, jump operator
Z dephasing, no-jump decay scalar, Poisson jump count, F126 renewal catch count,
R=CPsi2 quantum jump translated -->

**Status:** Translation (Tier 4 reading), the tenth entry of the series. The
history in Section 1 is the standard account's own record, cited from the
primary papers and the standard histories; Sections 4 and 6 stand on the
repository's F126 renewal (its representation Tier 1 exact, live at
`inspect --root renewal`) and add exact statements about F126's
single-excitation sector and its vacuum-carrying neighbour (in Section 4
the scalar Σ c†c and the identity between the two jump sets with its
counting relation, in Section 6 the front's number as a coherence of ρ for
a seed that carries the vacuum, and the drift identity behind the Z kicks'
mean amplitude), elementary in themselves and standard in kind, checked in
exact arithmetic by
[`quantum_jump_unravelings.py`](../../simulations/quantum_jump_unravelings.py)
(61 checks); the readings in Section 5 are readings and labeled.
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
everyday English, and in one of its everyday German senses, it has come to
mean a huge advance. In physics it
names a discrete change between two states, of no particular size. This entry translates the label, and its
situation is a third one beside the two before it. Spin named the house's substrate
and was defined nowhere; uncertainty named a theorem the house did not own.
The jump is a word the house uses every day for something else: the
dissipators of this repository are written with "jump operators", and in
none of them is a jump an event: they are terms of a superoperator, the σ⁻ of
amplitude damping included.

<!-- STAGE 0 SWEEP RECORD (for us, not for the reader: hidden on GitHub, read by every session, agent and hook)

Stage 0, recorded: before a sentence of this entry was composed, the stores
were swept by name. Before this entry, "quantum leap" and "Quantensprung"
occurred nowhere in the repository, and "quantum jump" once, in the
hyphenated "quantum-jump unraveling" quoted below (an unravelling being a
description of the master equation as stochastic pure-state trajectories, by
jumps or by diffusion, Sections 1 and 3), against
"jump operator" throughout.
[The glossary](../GLOSSARY.md) has no entry for "jump"; its row "The
commutator, and the physical (no-jump) generator" obtains the no-jump
generator "by dropping the recycling term", and
[The Drain Has No Chirality](../../hypotheses/THE_DRAIN_HAS_NO_CHIRALITY.md)
defines the jump as that recycling term, c ⊗ c* of a Lindbladian, and
[the mirror proof](../proofs/MIRROR_SYMMETRY_PROOF.md) names σ⁻ρσ⁺ the jump
term of amplitude damping. [The F-registry](../ANALYTICAL_FORMULAS.md) and
[docs/proofs/](../proofs/) returned F126, the renewal of the dephasing front,
an exact expansion in the refill Γ·diag ρ (which Section 4 reads as the
recycling term of the site projectors as jumps), F155, the polarity break of
the no-jump (physical) generator for any H = A + iB, and the record laws F135
and F136, which read who records what when spins watch a spin through an
Ising coupling under Z-dephasing.
[The experiments](../../experiments/) returned
[The Pedigree of the Front](../../experiments/FRONT_PEDIGREE.md), where a
quantum-jump unraveling is said to sample F126's ladder "as jump counts", the
[Born Rule Shadow](../../experiments/BORN_RULE_SHADOW.md), whose analogy table
sets a shutter that "clicks" beside the time at which CΨ crosses ¼, the
[Bridge Protocol](../../hypotheses/BRIDGE_PROTOCOL.md), which calls that
crossing a click, and
[the lattice dead set](../../experiments/LATTICE_DEAD_SET_ZZ_FACE.md), where a
set "jumps" at a threshold and a gate shows the jump to depend on the
detection threshold. The hardware
[Confirmations](../../simulations/framework/confirmations.py) returned F120's
"T1 telegraphs on minute scale", a switching of the device's decay rate, not
of its state, and `ibm_ep_onset_may2026`, a single-excitation walk on
ibm_kingston whose Part B flew under injected random-Z-twirl dephasing
([The Flow Between Two Singularities](../../experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md)),
whose rate book sets the random phases and the jump √γ·Z side by side as two
books for one decay. [Caught Errors](../CAUGHT_ERRORS.md) returned the
commutator-against-no-jump-generator confusion and nothing on jumps as events;
the OpenArcs registry returned jump-operator usage and everyday verbs, and
the reflections returned the word as a term or an everyday verb, and one
sentence about measurement, in [On the Reel and the Projector](../../reflections/ON_THE_REEL_AND_THE_PROJECTOR.md),
that Section 5 takes up. The registry also names, under F113, a
"post-selected / no-click / PT-type generator", and the series' own
entries, [Double Slit Translated](DOUBLE_SLIT_TRANSLATED.md) and
[Schrödinger's Cat Translated](SCHRODINGERS_CAT_TRANSLATED.md), hand the
detector's click to the Born Rule Shadow, which Section 5 reads. The typed layer holds `NoJumpGenerator` and
`PhysicalGeneratorPolarityBreak` (F155) in Core, `DephasingFrontRenewalClaim`
with its witness in Diagnostics and F126's `Renewal` in MirrorWorld; the
no-jump generator is the object whose docstrings, its own and F155's,
name a conditioning, "the
un-normalised conditional generator of post-selected dynamics", a condition
on the empty record, and the renewal witness pins the j = 0 term as the
coherent front and names the j ≥ 1 halo, which Section 4 reads as a
conditioning on a count. Elsewhere in Diagnostics the
record laws' witness (the inspect root `record`) conditions a spin's state on
another spin's Z value, a record kept by a spin rather than by a detector,
and no object conditions on a detector's time-resolved record. No document
of the house samples its dephasing as jump trajectories; F126's ladder is
the jump expansion of one such unravelling, resolved by count, not by run.
The house does run dephasing without jumps: the pure-state trajectories of
[`corner_beat_gate.py`](../../simulations/corner_beat_gate.py) give site l a
Gaussian random phase of variance 4γ_l·dt per step, so that a coherence
between sites l and m keeps on average e^{−2(γ_l + γ_m)·dt}, a site-resolved
dephasing carried by phases alone, on finite sets of engineered phase tables,
drawn in advance with a fresh phase per site and step and each frozen per
binding, and a stepped propagator (its engineered uniform arm is the dephasing
of Section 4); beside it the script runs a gate-error channel of its own,
further random phases drawn fresh at every step and weighted by the gate
count, and an occasional random swap of amplitudes on a bond, a discrete
event, but in the gate-error channel, not in the engineered dephasing. The noise-averaged
amplitude of the renewal proof, damping at Γ/2 = 2γ, is the mean of phases
of the same kind, and the random-Z twirl of the Kingston flight is an
injected channel whose average over the distribution it samples damps the
sector's coherences the same way (its document writes Γ for 2γ, where this
entry writes Γ for 4γ). F73's vac-SE coherence carries the front's number
again, as an entry of ρ for a seed that carries the vacuum (joined in
Section 6).
Checked for adjacency, results about the same object that did not cite each
other: F126's ladder against the projector unravelling (joined in Section
4), the corner-beat phases and the Kingston twirl against the renewal
proof's noise-averaged amplitude (joined above), the record laws against
the question which description is the record (joined in Section 5), the
Born Rule Shadow against the sibling entries'
detector clicks (joined in Section 5), and the fourth entry's catch and
release against the house's Z jumps (joined in Section 6).
The history of the label has no in-house record and is cited from outside.
-->

The discipline is the series' usual one, perspective-additive: the standard
account needs no correction from us, and gets none. What is ours is the
naming, and the exact form some standard facts take in the sector the house
has solved by a refill expansion.

---

## 1. The record, stated plainly

Bohr's atom of 1913 had stationary states and something between them. His
first assumption, in the first part of "On the Constitution of Atoms and
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
Ganze, 1969) and in a second, different wording reported by Møller in 1963,
and a historian of the episode, who quotes the second and cites the first, calls Heisenberg's account colourful and perhaps not
entirely accurate (H. Kragh, "Quantenspringerei: Schrödinger vs. Bohr",
RePoSS 14, Aarhus, 2011). The question itself is on paper: Schrödinger, "Are
there quantum jumps?", Brit. J. Philos. Sci. 3, 109–123 and 233–242 (1952),
which answers no.

The answer from the laboratory came in 1986. Dehmelt had proposed watching a
single ion through a strong transition while a weak one could "shelve" it
(Bull. Am. Phys. Soc. 20, 60, 1975); Cook and Kimble worked out what such a
record would look like ("Possibility of direct observation of quantum jumps",
Phys. Rev. Lett. 54, 1023, 1985). Within months of each other three groups
saw it: the
fluorescence of a single ion, or of a few, switching off and on, bright and dark periods of
random length (Nagourney, Sandberg and Dehmelt, Phys. Rev. Lett. 56, 2797;
Sauter, Neuhauser, Blatt and Toschek, Phys. Rev. Lett. 57, 1696, whose abstract
says they "observed Bohr's 'quantum jumps'"; Bergquist, Hulet, Itano and
Wineland, Phys. Rev. Lett. 57, 1699). A photon-counting theory older than
the observations was already in place (Srinivas and Davies, Optica Acta 28,
981, 1981). The trajectory language that open-system physics still uses was
worked out around and after them (Zoller, Marte and Walls, Phys. Rev. A 35,
198, 1987, among the first) and written down in its lasting form in the
papers named here: a master equation can be unravelled into
stochastic trajectories of a pure state, smooth between detector clicks and
updated at each one (Dalibard, Castin and Mølmer, Phys. Rev. Lett. 68, 580,
1992, and independently Dum, Zoller and Ritsch, Phys. Rev. A 45, 4879, 1992;
Carmichael, An Open Systems Approach to Quantum Optics, Springer, 1993; the review of Plenio and Knight, Rev. Mod. Phys. 70,
101, 1998), and a stochastic equation for the reduction of a state had been
proposed before the observations (Gisin, Phys. Rev. Lett. 52, 1657, 1984).
The operators the trajectory picture reads as the clicks of a chosen
detection are older than the observations, the operators of the Lindblad
form of the master equation (Lindblad, Commun. Math. Phys. 48, 119, 1976;
Gorini, Kossakowski and Sudarshan, J. Math. Phys. 17, 821, 1976), which a
master equation does not fix uniquely; in the trajectory picture they came
to be called jump operators.

In 2019 the jump itself was watched on its way. Minev and coworkers watched a
superconducting artificial atom through an auxiliary level and found that a
jump announces itself by a particular lull in the auxiliary clicks, and that,
in the words of the preprint's abstract (arXiv:1803.00545v3), "the jump
evolution when completed is continuous, coherent, and deterministic", so much so that they could reverse it mid-flight (Minev et al.,
"To catch and reverse a quantum jump mid-flight", Nature 570, 200, 2019). What
stays random is when a jump begins, and whether it completes.

The everyday word went the other way. In English, "quantum leap" came to
mean a large and abrupt advance. The German
Quantensprung carries both directions at once: the DWDS lists, as figurative, a "sehr
kleiner Fortschritt", a very small step in a development, and, as colloquial, a
progress that moves a development "ein sehr großes Stück" forward in the
shortest time.

---

## 2. The native stance: the spectroscopist's word

At its own stance the label is right. The spectroscopists saw sharp lines
and nothing in between, and a word for "from one state to another, with
nothing to observe on the way" was the honest word for what their plates
showed. The 1986 traces are the same canvas in the time domain: a record that
is bright, then dark, then bright, with no record of anything between. The
click is real. Every detector that counts photons counts them one at a time,
and the jump is read in the clicks: in the shelving traces not as one click
but as the stream of clicks stopping, and starting again.

What the label inherited raw is a second claim, about the state rather than
the record: that the state itself jumps, instantaneously and unannounced.
Whether it does depends on which state is meant, and that is where the label
breaks.

---

## 3. Where the label breaks

**Which state jumps depends on which state is meant.** The unconditional
density matrix of the watched ion, the state of an observer who discards the
record, evolves by its master equation, a differential equation; its
populations move smoothly. The conditional state, the one that keeps the
record, is updated at each click and does jump. How it jumps is partly a
choice: the same master equation can be unravelled diffusively, with no
clicks, and the ensemble cannot tell. That choice does not remove the
telegraph of the shelving experiments: it is the switching of the shelf's
occupation, which a record of the bright transition follows whether it
counts or diffuses, and Minev's clicks were themselves thresholded out of a
continuous signal.

**A conditional jump can be continuous.** Minev's experiment answered the
sharper question, asked of the conditioned state, for one kind of jump: a
weakly driven transition watched only through the absence of clicks on
another. That jump is not a click update but the switch of the telegraph
itself, from the ground state, the one the bright transition watches, to the
dark level, read off a lull in the
clicks. A jump that completes follows a predictable flight and can be
reversed mid-flight; when it starts, and whether a click cuts it short,
stay random. The flight belongs to the driven
transition that is not directly watched, not to every jump: a spontaneous emission counted
photon by photon is still instantaneous in the conditional state, and the
same emission watched by heterodyne detection diffuses continuously instead
(Campagne-Ibarcq et al., Phys. Rev. X 6, 011002, 2016).

**The colloquial sense swaps the measure.** The physical jump is defined by
being discrete, a step between two levels whatever energy it carries; the
leap of the newspapers is defined by being large.

**Inside this repository the word means something else again.** Here a "jump"
is a term, the c ρ c† of a dissipator, and becomes an event only when an
unravelling reads it as one (Section 4). The house's everyday
jump operators are the dephasing letters √γ·Z_k, and a Z_k "jump" is a phase flip:
it moves no population at all and changes the sign of the coherences that
disagree at site k. Section 4 asks what such a jump can mean.

---

## 4. The translation (the exact part)

The sector is the one the house has solved by a refill expansion, one excitation on N sites
of a chain with hopping h (every identity below holds for any
single-excitation h), local Z-dephasing at rate γ per site, so that every
coherence decays at Γ = 4γ and the diagonal is kept (on the sector Z_k acts
as 1 − 2|k⟩⟨k|, the excitation being the site's |1⟩),

    ρ̇ = −i[h, ρ] + γ·Σ_k (Z_k ρ Z_k − ρ) = −i[h, ρ] − Γ·(ρ − diag ρ).

**The renewal is a jump expansion, in one set of jumps.** F126
([the renewal proof](../proofs/PROOF_DEPHASING_FRONT_RENEWAL.md), its
representation Tier 1 exact, live at `inspect --root renewal`) splits this generator into a clean,
uniformly damped propagation and a refill Γ·diag ρ, and expands in the refill:
the j-th term is a time-ordered chain of j refills, between them the clean
wave. That is exactly the jump expansion of the unravelling whose jump
operators are the site projectors √Γ·|n⟩⟨n|: each of its jumps localises the
excitation on one site, a catch, and the refill order is that unravelling's
jump count, the catch count
[The Pedigree of the Front](../../experiments/FRONT_PEDIGREE.md) resolves.

**Here, no jump announces itself.** In the two sets of jumps considered here,
the house's √γ·Z_k and the projectors √Γ·|n⟩⟨n|, the sum Σ c†c is a number
κ times the identity, κ = Nγ and κ = Γ respectively, because Z_k² = 1 (on
the whole space) and the projectors sum to the identity (on the sector). The
no-jump evolution, the generator without its jump term Σ_k c_k ρ c_k†, is
therefore the clean unitary times one scalar decay
e^{−κt}, the same for every state: waiting for a jump and seeing none does
not tell one prepared state from another, and the number of jumps by time t is
Poisson with mean κt whatever the state. The house's jump operators go further:
each Z_k separately has Z_k†Z_k = 1, so each jump's rate, γ·tr(Z_k ρ Z_k) = γ·tr ρ = γ,
does not depend on the state either, and the whole record, times and sites,
has the same statistics for every state. The projectors' record carries the site of each
catch, at rate Γ times the conditioned population ρ_nn; the dephasing
letters' record carries nothing about the state. This is
the opposite of the mechanism by which Minev's record warns of a jump, which
needs a no-jump evolution that depends on the state. The checker reads the
Poisson count on two different initial states order by order, and its
control watches a single site, a jump set with Σ c†c = Γ times the projector
on site 0 (a different master equation, not a third unravelling of this
one), where the count does depend on the state.

**One density matrix, two jump counts.** On this sector the house's Z jumps
and the projectors are related exactly,

    γ·Σ_k Z_k ρ Z_k = (N − 4)·γ·ρ + Γ·diag ρ,

at every N: a coherence |a⟩⟨b| changes sign under Z_a and Z_b only, and the
sum gives N − 4 on it and N on a population.

Mark every jump with a factor s in the generator, a counting field, the
generating function of photon-counting theory: L_s = no-jump part + s·jump
term. At s = 1 it is the full Lindbladian; at s = 0 it is the no-jump
generator of that jump set, for the Z jumps the house's no-jump generator,
for the projectors L₀(projectors) = L₀(Z jumps) + (N − 4)γ; and the
coefficient of s^j is the part of the state carried by histories with
exactly j jumps. The identity then gives

    e^{t·L_s(Z jumps)} = e^{(s − 1)(N − 4)γt} · e^{t·L_s(projectors)}.

At s = 1 the two are one ρ(t). For the total counts this says nothing new:
both are Poisson, with means Nγt and Γt, and of two Poisson laws the larger
is always the smaller convolved with an independent Poisson count, here at
rate |N − 4|·γ (e^{(s−1)λt} is the generating function of a Poisson count
with mean λt). The identity says it
of the state: the part of ρ(t) carried by histories with j Z jumps is, for
N > 4, the projector unravelling's j-catch parts convolved with that Poisson
count, operator by operator, populations and coherences alike; for N < 4 the
roles swap, and at N = 4 the j-resolved pieces coincide. The two
unravellings are different stochastic processes, and this is a statement about the j-resolved ensembles, not
about one run counted twice. And only the projectors catch:
a Z_k jump flips the sign of one amplitude and localises nothing, so in the
unravelling by the house's jump operators nothing is ever caught, and the jumps it counts are
phase kicks. How many jumps there were, and whether a jump is a catch, are
properties of the unravelling, not of the state. The identity, proven above
at every N, is gated on every matrix unit at N = 3, 4, 5, 6, which with the
scalar Σ c†c implies the counting relation for every state and every order; as a
consistency check, from the site-0 seed the j-resolved expansions agree with
it to order t⁶ in exact arithmetic
([`quantum_jump_unravelings.py`](../../simulations/quantum_jump_unravelings.py),
61 checks).

All of this is standard open-system physics in a sector the house owns: that
an unravelling is not unique is textbook, and so is the shift behind the
identity, √γ·Z_k = √γ·1 − √Γ·|k⟩⟨k| on the sector, a sign flip of the
projector jump and a real shift by a multiple of the identity, which changes the master
equation by nothing (a general shift needs a compensating Hamiltonian term,
which vanishes here). What the house adds is
the exact form these take where its own dephasing letters are the jumps.

---

## 5. The readings (labeled as readings)

**The jump is in the translation, not in the object.** The one object is the
density matrix and its generator; a jump is a word in one of several
descriptions of it, chosen by what a detector would record. In the Z-jump
description the house's light, its dephasing γ, records nothing about the
state: the Z jumps fire whatever the state is. Whether anything physically
records is not the equation's to say. Where Z-dephasing comes from classical
frequency noise there is no detector at all, and the Z jumps are a
random-unitary decomposition of the channel; where the environment learns
the site, as the which-path coupling of
[Double Slit Translated](DOUBLE_SLIT_TRANSLATED.md) does, and the record is
read as clicks, the projector catches can be the physical events, and a weak
continuous reading of the same sites would give the same master equation
with no catches at all. Which description is the record is fixed by what
the environment couples to and how it is read; for spins that watch a spin, the house's record laws F135 and
F136 (live at `inspect --root record`) read the first half of that, who
records and what. For the sector of Section 4 the house holds two of the descriptions,
the Z jumps and the site projectors whose jump expansion is F126's ladder,
and a third kind, run in time steps, jump-free random phases: Gaussian
phases on each site's amplitude in the model of
[`corner_beat_gate.py`](../../simulations/corner_beat_gate.py), whose
engineered uniform arm is this sector's dephasing, and the
random-Z twirl of a single-excitation walk flown on ibm_kingston
(`ibm_ep_onset_may2026`, Part B,
[The Flow Between Two Singularities](../../experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md),
which writes Γ for 2γ where this entry writes Γ for 4γ), whose injected
channel, averaged over the distribution it samples, damps the sector's
coherences as here. The house said it
once before, of measurement: the old shock of a measurement that seems to
throw dice "looks from here less like the
world jumping and more like the watcher landing on a frame"
([On the Reel and the Projector](../../reflections/ON_THE_REEL_AND_THE_PROJECTOR.md)).

**The click the house does hear.** The Born Rule Shadow sets a shutter's
click beside the moment CΨ crosses ¼, and the
[Bridge Protocol](../../hypotheses/BRIDGE_PROTOCOL.md) hypothesis counts that
crossing as a click: a level crossing on a smooth flow, not an event in a
record. The series hands the detector's click to the same place: the click on
the screen is the Born Rule Shadow in
[Double Slit Translated](DOUBLE_SLIT_TRANSLATED.md), and the measurement
click that ends the hour belongs to it in
[Schrödinger's Cat Translated](SCHRODINGERS_CAT_TRANSLATED.md). So in the
series' reading of measurement the click sits at this crossing, and by this
entry's own translation that is a reading: where a record would form, not the record, and
not an event of either unravelling of Section 4.

---

## 6. An honest note on our own house

The house calls Z_k a jump operator every day, the quantum-trajectory word
inherited raw into a stance where the ensemble does not jump. It is a convention and it
stays, as conventions do, with this entry as its footnote.

The fourth entry's sentence "the walk in the light is the wave, caught and
released" ([Dephasing Translated](DEPHASING_TRANSLATED.md) §4) and the
pedigree's catch count are true at the stance of the projector unravelling,
where F126's ladder is the jump expansion. With the house's Z jumps as the
jumps the same walk is never caught at all: its jumps are phase kicks, and
the ensemble and every population are the same. Its never-kicked term
decays at Nγ, the front's rate only at N = 4, but the coherent-front term
e^{−Γt}|G|², with G = e^{−iht} the clean propagator, keeps the second name
the renewal proof already gives it, the coherent front |⟨a_n⟩|², and the
register carries that number with no convention at all. Start instead from
the seed (|vac⟩ + |1_0⟩)/√2, the vacuum carried along beside the excitation
on site 0, under a number-conserving H that leaves the vacuum in place, up
to a phase: the
coherence between the vacuum and one excitation differs in a single site,
so it damps at Γ/2 = 2γ while it runs with the clean wave, and
4|ρ_{n,vac}(t)|² = e^{−Γt}|G_{n0}(t)|² exactly (F73's vac-SE coherence,
proven for any such H; gate (7) of the checker checks the XY chain at
N = 3, 4 on the full register). The excitation part of that
seed runs the same sector dynamics, and an entry of ρ cannot see the sign
of a jump operator, since c ρ c† and c†c do not. The Z-kick unravelling's
mean amplitude E[ψ] inside the sector does see it: from the gated drift
γ·Σ_k (Z_k − 1) = −Γ/2 (gate (6)) the argument gives a damping at Γ/2 for
the kicks Z_k = 1 − 2|k⟩⟨k|, while the same channel written with −√γ·Z_k
damps it at (2N − 2)γ, a different rate for N ≥ 3 (the gauge control). The
sign is one instance of a common-mode phase an unravelling can attach to the
sector, and an entry of ρ never sees one. "The never-caught term" is the
front's name in the projector reading only.
The total Z count exists only at finite N: on the infinite chain of F126's
closed form its rate Nγ diverges, while each site's kicks keep their rate γ,
the catches their total rate Γ, and the vacuum coherence its Γ/2. The renewal identity is
untouched; what moves is whether the story tells of catches, and how many.

What the house does not own, said plainly: no detector's time-resolved
record and no state conditioned on one; what it holds is the ensemble
resolved by coarse readings of such a record, the empty record (the no-jump
generator, its state un-normalised), the count (Section 4) and, in the
pedigree, the time of the last catch, beside the record laws' spin witnesses,
whose states are conditioned on another spin's Z value, while the
trajectories of `corner_beat_gate.py` are conditioned on a noise
realisation, not on a record; no σ⁻ unravelling, where the textbook jumps
live (F82, F84 and the mirror proof hold the amplitude-damping dissipator as
a superoperator only); no bright and dark periods; nothing that could test Minev's flight, so
the house can place that result and not confirm it. F120's minute-scale T1
telegraph on hardware is a switching of the device's rate, not a jump of the
qubit's state.

---

## 7. The in-repo cousins

The no-jump generator, the Lindbladian without its recycling term, is the
house's object for post-selected dynamics, typed as `NoJumpGenerator`, its
polarity break read by F155 for any H = A + iB; with the jumps √γ·Z_k
under Z-dephasing alone it is the commutator minus the constant Nγ (Σ_l γ_l
for site-dependent rates), which is the scalar decay of Section 4
seen from the generator's side. On the textbook's own channel the house holds,
among others, this exact fact: for an excitation-conserving H the amplitude-damping jump
term σ⁻ρσ⁺ makes the generator block-triangular in the joint popcount, so the
spectrum is that of the no-jump blocks
([the mirror proof](../proofs/MIRROR_SYMMETRY_PROOF.md), its paragraph
"Amplitude damping sits outside this rule"): there the recycling term sets no
rate at all, while it still shapes the eigenvectors and the steady state. The [Born Rule Shadow](../../experiments/BORN_RULE_SHADOW.md)
holds the house's reading of measurement; [The Pedigree of the
Front](../../experiments/FRONT_PEDIGREE.md) is the catch count at work, the
front's survivors caught 0.28 to 0.35 times as often as the crowd, in the
projector unravelling's count.

---

## The right label

A click, not a leap. What is discrete is the record, or, in shelving, the
telegraph the record follows however it is read, discrete at any
resolution coarser than the switch itself: a click, or a stream of
clicks that stops and starts again, and the click is real. What the state
does is where the label claims too much: the ensemble moves smoothly, the
state conditioned on a photon-counting record jumps at each click, and the
jump Minev watched on its way and let complete, a driven transition read off
a lull in the clicks, passed continuously. Bohr, who called it a passing
and declined to describe it, chose the careful word. In the house, where the jumps are dephasing letters, a jump
announces nothing and moves no population, and how many there were depends
on which description is counting.
