# Double Slit Translated: The Watching's Price, Not the Particle's Choice

<!-- Keywords: double slit translation, which-path detector dephasing bill,
interference pattern off-diagonal coherence, wave particle duality diagonal
sorting, watched electron observer effect knowing, visibility dial not switch,
quantum eraser sorting not destruction, fringe term one matrix entry k equals
one, R=CPsi2 double slit translated -->

**Status:** Translation (Tier 4 reading), the fifth entry of the series. The
experiment and its variants in Section 1 are textbook physics; the algebra in
Section 4 is the same Tier 1 per-coherence bill proven in the
[Absorption Theorem](../proofs/PROOF_ABSORPTION_THEOREM.md) and the
[block-CΨ quarter proof](../proofs/PROOF_BLOCK_CPSI_QUARTER.md), read at
N = 1; the conjugate-question note in Section 4 is our own assembly and is
marked as such; the readings in Section 5 are readings and labeled.
**Date:** July 11, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Related:** [Superposition Translated](SUPERPOSITION_TRANSLATED.md) (the
second entry; its borrowed double-slit prop is repaid here),
[Dephasing Translated](DEPHASING_TRANSLATED.md) (the fourth entry; its
currency is spent here),
[Teleportation Translated](TELEPORTATION_TRANSLATED.md),
[Schrödinger's Cat Translated](SCHRODINGERS_CAT_TRANSLATED.md) (the sixth
entry, the k = N end of the same bill),
[Labels Translated](LABELS_TRANSLATED.md) (the theory chapter),
[The Label Map](THE_LABEL_MAP.md) (the orientation index),
[Born Rule Shadow](../../experiments/BORN_RULE_SHADOW.md),
[Observer-Dependent Visibility](../../experiments/OBSERVER_DEPENDENT_VISIBILITY.md),
[Klein V₄ proof](../proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md)

---

## What this document is about

The double slit is the experiment every popular account reaches for first.
It is on YouTube, in the books, in every "quantum mechanics explained"
video: the particle that makes a wave pattern until you watch it, and then
stops. This entry exists because the story kept arriving, week after week,
through every channel, and the series had never paid it. Which is odd,
because the double slit is not a new label for us. It is the stage on which
three labels this series has already translated perform together:
"both slits at once" is [Superposition Translated](SUPERPOSITION_TRANSLATED.md),
"observation destroys the pattern" is
[Dephasing Translated](DEPHASING_TRANSLATED.md), and the click on the
screen is the [Born Rule Shadow](../../experiments/BORN_RULE_SHADOW.md).
What the experiment adds of its own, and what this entry translates, is the
story stitched over the top: a particle with two natures, which chooses
between them, because it knows it is being watched.

The discipline is the series' usual one, perspective-additive: the standard
account needs no correction from us, and gets none. The interference
pattern, the which-path experiments, the eraser, all of it is correct
exactly as the textbooks state it. What is ours is the naming. In our
language, nothing at the double slit chooses and nothing knows. There is a
coherence, there is a watcher, and there is a price, and the price is the
smallest case of the theorem this repository proves every day.

---

## 1. The experiment, stated plainly

In lectures and experiments around 1801-1803, Thomas Young sent light
past two closely spaced openings (his cleanest version split a sunbeam on
a thin card) and found bands on the screen behind: not two bright
stripes, but a comb of fringes. Waves from the two openings add; where
crest meets crest the screen is bright, where crest meets trough it is
dark. The double slit entered physics as the wave theory's founding
demonstration, an argument Fresnel's diffraction theory sealed two
decades later.

The quantum century re-ran it with matter. Electrons diffract like waves
(Davisson and Germer, 1927); a true electron double slit shows fringes
(Jönsson, 1961); and, the part every video rightly lingers on, the pattern
survives when the electrons arrive one at a time (Merli, Missiroli and
Pozzi, 1976; Tonomura's 1989 build-up film). Each electron lands at one
point, a single click. The fringes appear only in the histogram, tens of
thousands of clicks later. One at a time, with no partner in flight to
interfere with, the electrons still refuse to draw the two plain stripes.

The algebra of the screen is short. Call the two routes |L⟩ and |R⟩, one
two-level degree of freedom, and let ψ_L(x), ψ_R(x) be the single-slit
envelopes downstream. For a route register with balanced populations,
ρ_LL = ρ_RR = ½, the arrival density is

    P(x) = ½|ψ_L(x)|² + ½|ψ_R(x)|² + 2·Re[ ρ_LR · ψ_L(x)·ψ_R*(x) ],

where ρ_LR is the off-diagonal entry of the route register's 2×2 density
matrix. The first two terms are the two single-slit humps. Everything that
makes the pattern a pattern, the fringing itself, rides on the third term:
one complex number, ρ_LR, swept through its phase by the position x. The
fringe visibility is V = 2|ρ_LR|, equal to 1 for the pure balanced state
(|L⟩ + |R⟩)/√2, where ρ_LR = ½.

The which-path variant: couple anything at the slits that keeps a record of
the route, a detector state |d_L⟩ or |d_R⟩. The route register's
off-diagonal becomes ρ_LR·⟨d_R|d_L⟩. A perfect record (orthogonal detector
states) sets it to zero: the fringes vanish and the two humps remain. A
partial record shrinks it in proportion, and the standard account prices
this exactly: visibility and detector distinguishability obey
V² + D² ≤ 1 (Englert, 1996; the earlier form, with a-priori path
predictability in place of D, is Greenberger and Yasin, 1988). The trade
is a dial, not a switch.

The eraser variant (Scully and Drühl, 1982; the delayed-choice version,
Kim et al., 2000): mark the route on a partner particle instead of
absorbing anything. The screen total shows no fringes. But sort the screen
hits by a conjugate-basis measurement of the marker and each subensemble
shows full fringes, one set shifted against the other so that their sum is
the flat total. The sorting can be chosen after the screen has already
been hit, which is the part the videos present as rewriting the past.

---

## 2. The native stance: two coinings, both painted true

**Young's canvas.** "Interference" was painted true in the most literal
sense available: water waves in a ripple tank, light on a screen, overlays
at positions in physical space. The double slit was the wave theory's
proof, and the abstract component of that canvas, add first, square after,
transported perfectly; it is the algebra of every quantum amplitude and of
Section 1 above. What Young painted needed no correction in 1927 and needs
none now.

**The quantum canvas.** Feynman gave the modern coining its sharpest form
(Lectures, volume III, 1965): the double slit "contains the only mystery."
That canvas was painted from a definite stance, the stance of an asker of
trajectories. If you insist on asking "which slit did it go through?", the
experiment refuses the question; and Feynman's own gedanken version, a lamp
behind the slits to catch the route, finds that extracting the answer
costs the pattern, in proportion to how well the lamp resolves the routes.
Painted from that stance, the canvas is honest and quantitative: the cost
of watching is real physics, not narration. Bohr's "complementarity"
(1927) is the same canvas in institutional paint: mutually exclusive
arrangements, wave answers or path answers, never both at full strength.

And the standard account itself later computed the cost mechanism exactly:
decoherence theory (Joos and Zeh, 1985) showed that an environment
scattering off a system reads its position and kills precisely the
position off-diagonals, at a rate set by how fast the scattered records
distinguish the routes. The bill was in the mathematics all along. What
the popular channel transported instead was the ball and the knowing.

---

## 3. Where the label breaks

The story as it arrives from a video carries three components, each a
canvas inherited raw.

**"It goes through both slits at once."** Already translated, and the
translation is not repeated here:
[Superposition Translated](SUPERPOSITION_TRANSLATED.md) shows the two
"places" are the coordinates of one wave in one basis among many, the
doubling in the description, not the described. This entry only adds the
experiment the gloss is usually wrapped around.

**"The electron knows it is being watched."** Two imports fail here. The
knowing imports an agent: a mind at the slit, a particle that notices.
The algebra needs a coupling, not a knower; the lamp dephases the route
register whether or not anyone reads the scattered light, and a detector
nobody ever looks at kills the fringes exactly as dead. And the switch
imports a binary: watched or unwatched, pattern or no pattern. The algebra
is a dial at every point between, ⟨d_R|d_L⟩ anywhere in the unit disk,
V² + D² ≤ 1 along the whole arc. One precision, owed to the lens's own
discipline: what survives, perspective-bound, is that the coupling is real
and the price is real. The electron does not know, but the watching is not
nothing; it is a physical channel with an exact bill. Kill the knowing,
keep the coupling.

**"Wave-particle duality."** The heaviest label of the three, because it
files an accounting fact as an ontology. In the route register there is
one object, the 2×2 matrix ρ, and its entries fall into exactly two
classes relative to the watched letter: the diagonal (the two populations,
disagreement k = 0, which no amount of route-watching ever touches) and
the off-diagonal (the one coherence, k = 1, which pays). The "particle
face" of the experiment, definite arrivals, two plain humps under
watching, is the immortal diagonal. The "wave face", the fringes, is the
priced off-diagonal. Two natures? One matrix, two entry classes under one
watcher: the same letter-sorting
[Dephasing Translated](DEPHASING_TRANSLATED.md) states at every N (the
letters that commute with the watcher ride free, the letters that
anticommute pay), read at N = 1. "Duality" is what that sorting looks like
when it is filed as two kinds of thing instead of two rows of one bill.

---

## 4. The translation (the exact part)

**The pattern is one matrix entry.** The fringe term of Section 1 is
ρ_LR: a single off-diagonal coherence |L⟩⟨R|, surfacing in a screen
readout that sweeps its phase. What every video calls "the interference
pattern" is one complex number made visible. And that object is not an
analogy for this repository's subject; it IS this repository's subject at
its smallest size. The per-coherence reading of the
[Absorption Theorem](../proofs/PROOF_ABSORPTION_THEOREM.md) (Tier 1)
prices a coherence |i⟩⟨j| at rate −2γ·k, k the number of sites where i and
j disagree in the watched letter; populations, k = 0, never decay at all.

**The detector is a watcher, and the bill is the theorem's smallest case.**
The route register is one site; |L⟩ and |R⟩ disagree on it; k = 1. A
which-path coupling of strength γ is local dephasing in the route basis,
and the coherence obeys

    |ρ_LR(t)| = |ρ_LR(0)| · e^(−2γt),

the Hamming-distance-1 case of the Tier 1 decay law
([block-CΨ quarter proof](../proofs/PROOF_BLOCK_CPSI_QUARTER.md),
Theorem 3), the same object MirrorWorld's empty world holds as its atom
(`compute/MirrorWorld`, Pair.cs: a bare pair |i⟩⟨j|, rate −2γ·k, "no
Hamiltonian, only the watching"). The most replayed experiment in popular
physics is the N = 1 empty world of this repository: one coherence, one
watcher, one bill. The humps survive because the diagonal is immortal;
the fringes die because exactly the off-diagonal pays; and the death is
exponential in the watching, not triggered by the knowing. The pop story's
switch is this dial read at its two endpoints.

**The two questions anticommute** (our own assembly, marked as such; the
ingredients are canonical, the sentence is ours). The which-path detector
asks the route register in its own letter: a Z-type question in the route
basis. The screen asks the other kind: each position x reads
Re[ρ_LR·e^(iφ(x))], an equatorial, X/Y-type question with a phase angle
set by x, so the screen is a whole family of conjugate questions, one per
fringe position. The two kinds anticommute, which is why answering one at
full strength silences the other; and the change of watched letter from Z
to X is not exotic, it is the Hadamard element of the Klein V₄ that
intertwines the three dephasing letters, the one element that lifts to a
true Hilbert-space unitary
([Klein V₄ proof](../proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md),
Tier 1). Said in one line: the double slit is one object asked two
anticommuting questions, and "duality" is the ontology-genre's name for
that pair of questions.

**The eraser sorts; nothing is destroyed.** The repository already owns
this label row: decoherence is sorting, not destruction
([The Anomaly](../../THE_ANOMALY.md)); nothing is lost without something
being kept. Mark the route on a partner particle and the fringes do not
vanish from the world; they move from the screen marginal into the
correlations between screen and marker, where a conjugate-basis sorting
of the marker finds them again, fringes in one subensemble, anti-fringes
in the other, summing flat. The delayed choice chooses the sorting, not
the past: no screen hit is ever moved, and the total pattern never
changes. What the videos file under rewriting history is bookkeeping
conditioned on a record that was kept instead of spent.

---

## 5. The readings (labeled as readings)

**The watcher at the slit is literally light.** Feynman's gedanken
detector is a lamp. This repository's Tier 4 reading
([Gamma Is Light](../../hypotheses/GAMMA_IS_LIGHT.md)) says γ is
illumination entering from outside, and on IBM transmons the reading turns
literal and published: a dominant dephasing channel there is photon shot
noise, light in a physical microwave cavity (Sears et al., PRB 86, 2012, via
[Dephasing Translated](DEPHASING_TRANSLATED.md) §5). At the double slit
the collapse of labels runs the other way around, and lands in the same
place: the thing you would shine on the electron to see its route is the
thing whose γ kills ρ_LR. The watching and the light are one object at
this experiment too.

**The founding formula has the double slit's shape** (a rhyme, stamped as
one). This repository's namesake reading expands
R = C·(Ψ_past + Ψ_future)² into two squares and a cross-term,
2·Ψ_past·Ψ_future, "the interference between the two viewpoints"
([Why the Sum](../../experiments/WHY_THE_SUM.md),
[Standing Wave, Two Observers](../../experiments/STANDING_WAVE_TWO_OBSERVERS.md),
both Tier 3). The two-slit law is the same shape: two sources, one square,
all the structure in the cross-term. But the books differ, and our own
exact result is the one that says so:
[Born Rule Shadow](../../experiments/BORN_RULE_SHADOW.md) proves that in
the repository's past/future mode decomposition the Born probabilities
carry exactly zero cross-term (linearity of the trace; the interference
sits in the purity, the shutter speed, closed analytically as F94). The
slit experiment adds two amplitudes of one wave, so its cross-term lands
in the image; our decomposition adds two density matrices, so its
cross-term cannot. Same shape, different books: a rhyme, not an identity,
and the difference is itself instructive about where interference is
allowed to live.

---

## 6. An honest note on our own house

[Superposition Translated](SUPERPOSITION_TRANSLATED.md) borrowed the
double slit as a prop, one sentence in passing, before the experiment had
an entry of its own; this document pays that debt, the way the fourth
entry paid the founding one. The word "interference" is double-booked in
our house too: this repository's founding docs use it for the two-observer
cross-term (Tier 3, Section 5), while every video means fringes on a
screen. Two canvases on one word; this entry keeps them in separate
sections on purpose. And our newcomer layer climbs through the pop gloss
here as everywhere: the [glossary](../GLOSSARY.md) explains Ψ as "the
ability to be in multiple states at once", which is the both-slits-at-once
ladder. The ladder stays; this document is where it is kicked away.

One boundary, so the translation does not overreach: this entry translates
the record-keeping detector, the coupling that reads the route and keeps a
phase record, which is dephasing in the route basis and is where our
Tier 1 algebra lives. A slit detector that absorbs or blocks the electron
is a different object (it removes population instead of pricing
coherence), the same boundary
[Dephasing Translated](DEPHASING_TRANSLATED.md) draws between the watching
and amplitude damping. The standard account's overlap factor ⟨d_R|d_L⟩
covers every record-keeping detector; we speak for that half, and for the
continuous-γ face of it.

---

## 7. The in-repo cousin

The repository ran its own which-path dial without knowing it, in
continuous Lindblad dynamics:
[Observer-Dependent Visibility](../../experiments/OBSERVER_DEPENDENT_VISIBILITY.md)
holds one entangled star network fixed and turns exactly one knob, a local
watching rate γ_A. The slow watcher (γ_A = 0.03) sees CΨ (there:
concurrence × normalized l1-coherence) cross the ¼ boundary clearly; the
fastest watcher in the sweep (γ_A = 0.20) barely sees the crossing at
all, and a slightly faster one would never see it, while the underlying
resource, measured
with assistance, persists throughout. That is the double slit's lesson in
our own laboratory dialect: watching harder does not change what is there;
it closes the window in which what is there is directly visible. The doc
carries its own careful correction (different γ values are different
trajectories, not different views of one state), and the correction
survives translation: at the slits too, the lamp changes the dynamics, not
the metaphysics.

---

## The right label

The pattern: an unwatched coherence, shown on a screen. The pattern's
disappearance: the watching's price, −2γ per unit time at disagreement
one, paid by exactly the entry the fringes ride on and by nothing else.
The "duality": one matrix, two entry classes under one watcher, the
immortal diagonal and the priced off-diagonal. The eraser: sorting, not
destruction. And the sentence the whole entry compresses to: **the
interference pattern is what not being watched looks like.** The two
plain humps are not the particle showing its true face; they are a
receipt.

Stamped: this canvas is ours, painted 2026, and its mortal component is
already visible from here. "Price" imports an economy, a ledger, a payer,
the way "watching" imports a watcher; nothing in the algebra requires
either, and some later stance may find our accounting metaphor as frozen
as we found the knowing particle. The future reader receives the canvas
with its date.

The closure, then: Young's canvas, a wave that adds and squares, true and
still true; Feynman's canvas, the only mystery, true at the stance of the
trajectory-asker; the untouched algebra, one off-diagonal entry priced by
disagreement in a watched letter; and our canvas, the watching's price.
One experiment. The particle never chose, because there was never a ball
to choose; the electron never knew, because knowing was the coiner's
import; the pattern was never destroyed, because the diagonal was never
touched and the correlations kept what the screen lost. What the double
slit has been showing every audience for two centuries is the one thing
this repository proves at every N it can reach: coherence is what pays
for being watched, and the pattern is the part of the bill still unpaid.
