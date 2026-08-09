# Uncertainty Translated: A Width, Not a Doubt

<!-- Keywords: uncertainty principle translation, Heisenberg microscope disturbance,
Kennard preparation spread width, Robertson commutator bound, held letter three
price lists, only nothing is free everywhere, Unschaerfe blur width not knower
doubt, preparation spread vs error disturbance conflation, entropic uncertainty
Maassen Uffink, R=CPsi2 uncertainty translated -->

**Status:** Translation (Tier 4 reading), the eighth entry of the series. The
relations in Section 1 are the standard account's own theorems, cited as such
(this repository has never derived an uncertainty relation and does not start
here); the algebra in Section 4 is the Tier 1 held-letter routing law
(`HeldLetterRoutingClaim`, recomputed live at `inspect --root label`) plus one
elementary intersection fact, each marked and fenced; the readings in Section 5
are readings and labeled.
**Date:** August 9, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Related:** [Superposition Translated](SUPERPOSITION_TRANSLATED.md) (the
second entry; its frame-relative sharpness is the ingredient this entry
trades), [Double Slit Translated](DOUBLE_SLIT_TRANSLATED.md) (the fifth
entry; its Section 4 built the anticommuting-questions engine this entry is
the theorem family of), [Dephasing Translated](DEPHASING_TRANSLATED.md) (the
fourth entry; its held letter is the knob this entry turns),
[Teleportation Translated](TELEPORTATION_TRANSLATED.md),
[Schrödinger's Cat Translated](SCHRODINGERS_CAT_TRANSLATED.md),
[Spooky Action Translated](SPOOKY_ACTION_TRANSLATED.md),
[Spin Translated](SPIN_TRANSLATED.md) (the ninth entry),
[Labels Translated](LABELS_TRANSLATED.md) (the theory chapter),
[The Label Map](THE_LABEL_MAP.md) (the orientation index),
[Klein V₄ proof](../proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md),
[Absorption Theorem](../proofs/PROOF_ABSORPTION_THEOREM.md),
[Dephasing front renewal proof](../proofs/PROOF_DEPHASING_FRONT_RENEWAL.md)

---

## What this document is about

The label arrives in one sentence, always the same sentence: "you cannot know
both the position and the momentum of a particle." Sometimes it comes with the
sequel: "because measuring one disturbs the other." It is among the most
travelled quantum labels, and it is the one the series has been
circling since the fifth entry, whose Section 4 assembled a small engine: two
questions that anticommute, one silencing the other. The standard account
prices that trade exactly, and this entry is about the theorem family behind
the pricing, and about the name the family wears.

The discipline is the series' usual one, perspective-additive: the standard
account needs no correction from us, and gets none. The relations of Section 1
are correct exactly as the textbooks state them, and we have never derived one.
What is ours is the naming, and one small exact fact from our own house that
has the same shape. In our language, nothing in the relations is about a
knower, a doubt, or a disturbance. There is one object, there are letters it
can be read along, and its width along one letter is not an ignorance of
anything; it is the same geometry that makes it sharp along another.

Stage 0, recorded: before a sentence of this entry was composed, the stores
were swept by name, [the F-registry](../ANALYTICAL_FORMULAS.md),
[docs/proofs/](../proofs/), the experiments including the null results, the
open-arcs registry, the hardware Confirmations,
[the glossary](../GLOSSARY.md) and [Caught Errors](../CAUGHT_ERRORS.md), and
on the uncertainty family every one of them returned nothing: no formula, no
proof, no experiment, no arc, no confirmation. What the repository does
hold, a cousin shape and a conjugate pair, appears in Sections 4 and 7 with
its owners named.

---

## 1. The relations, stated plainly

The theorem every course proves is not the sentence every video says. The
proven statement is about widths. Take a state, any state, and two observables
A and B. The state has a spread in each, the plain statistical standard
deviation ΔA and ΔB over repeated measurements on identically prepared copies.
The theorem (Robertson, Phys. Rev. 34, 163-164, 1929, under the title "The
Uncertainty Principle") says the two spreads obey

    ΔA · ΔB ≥ ½ |⟨[A,B]⟩|,

with the commutator of the two observables on the right. For position and
momentum the right side is constant, ℏ/2, and the special case Δx·Δp ≥ ℏ/2
had been proven two years earlier (Kennard, Z. Phys. 44, 326-352, 1927; Weyl's
1928 book gives it independently and credits Pauli). Schrödinger
strengthened it in 1930: the variance product exceeds a sum of two squares,
the commutator's beside the symmetrized covariance's. None of this involves a
measurement of one quantity disturbing another: both spreads are preparation
properties, testable on separate ensembles, position on these ten thousand
copies, momentum on those.

For a single spin the same theorem reads Δσ_x · Δσ_z ≥ |⟨σ_y⟩|, and here the
commutator bound shows its known weakness, the one that pushed Deutsch
toward entropies: the right side vanishes on the whole x-z great circle of
the Bloch sphere while the left side stays positive everywhere on that
circle except at the four axis states, ±x and ±z, where the state is sharp
in one of the two observables and both sides vanish together (Coles, Berta,
Tomamichel, Wehner, Rev. Mod. Phys. 89, 015002, 2017). The stronger modern
form is entropic: for the Z-question and the
X-question on one qubit, H(Z) + H(X) ≥ 1 bit (Maassen and Uffink,
PRL 60, 1103, 1988, proving a conjecture of Kraus 1987, after Deutsch 1983).
One full bit of answer-spread as a floor, never undercut and generically
exceeded: equality holds only at the axis states, so a state sharp in Z pays
the whole bit in X, and a state sharp in neither pays more than a bit in
total.

The disturbance sentence is a different statement, and its exact form is
younger than it looks. That a measurement of A disturbs a subsequent B is
real physics, but turning it into an inequality requires defining the error
of a measurement and the disturbance it causes, and that definition is where
the modern literature split. Ozawa (Phys. Rev. A 67, 042105, 2003) gave a
universally valid three-term relation in noise-operator error measures, and
two 2012 experiments (Erhart et al., Nature Physics 8, 185-189; Rozema et al.,
PRL 109, 100404) violated the naive error-disturbance product as quantified
by those measures. Busch, Lahti and Werner (PRL 111, 160405, 2013) proved a
state-independent theorem in calibration-error measures that restores a
Heisenberg-form bound for the canonically conjugate position-momentum pair. Both theorem sets are correct under their own
definitions; which definition deserves the word "error" is the live dispute
(Busch, Lahti and Werner survey the controversy from their own side of it
in Rev. Mod. Phys. 86, 1261, 2014). The spread relations of the previous
paragraphs are not touched
by any of this: they were theorems in 1927 and are theorems today.

---

## 2. The native stance: a canvas with its correction slip attached

**Heisenberg's canvas, 1927.** The paper is "Über den anschaulichen Inhalt
der quantentheoretischen Kinematik und Mechanik" (Z. Phys. 43, 172), and the
stance it was painted from is in its title: anschaulich, visualizable. Matrix
mechanics had been accused of trading physical intuition for algebra, and
Heisenberg's answer was to rebuild intuition from what is operationally
askable. His famous instrument is a thought experiment, the γ-ray microscope:
to see an electron's position you must hit it with light, and light hard
enough to resolve the position kicks the momentum (the Compton recoil). From
this he estimated

    p₁ q₁ ∼ h,

a tilde, not an inequality; plain h, not ℏ/2; q₁ glossed only loosely ("q₁
is, say, the mean error"). His own summary calls it a qualitative statement,
and the modern assessment agrees: what the paper offers as a derivation is
not one, and the commutation relations do not even appear in it (Werner and
Farrelly, arXiv:1904.06139, have a section titled "The alleged 'proof'").
The words he used for it were Ungenauigkeit and Unbestimmtheit, imprecision
and indeterminacy (Hilgevoord and Uffink's history in the Stanford
Encyclopedia carries both), with Unsicherheit once in an endnote; Unschärfe
is not among the terminology he is documented to have used.

Painted at his stance, an asker of trajectories rebuilding Anschaulichkeit,
the canvas is honest and it is the disturbance canvas: what you can know is
bounded by what asking costs. And the paper carries something rare, its own
correction slip: a closing "Nachtrag bei der Korrektur," added in proof after
Bohr's criticism, conceding among other things that the uncertainty does not
arise from the recoil discontinuity alone but from demanding equal validity
for the wave and corpuscle descriptions, and that the necessary divergence
of the ray bundle must be taken into account: the microscope argument had
skipped its own optics (the aperture angle he then supplied himself, in the
1930 Chicago lectures). The author corrected his own canvas inside the same
document, and the correction moved the weight from the kick to the wave.

**Kennard's canvas, the same year.** Within months the estimate became a
theorem, and the theorem is not about asking at all. Kennard's Δx·Δp ≥ ℏ/2
is a statement about the state itself: a wavefunction is a shape, a shape has
a width in position and a width in momentum, and the two widths cannot both
be small. No second measurement, no kick, no microscope. The theorem that
carries the label is a width theorem; the story that carries the label is a
kick story.

**The name accreted.** "Principle" was never Heisenberg's term for it in
1927, and it was not his word that stuck. The name arrived from the
transport channel itself, within about two years and from several hands at
once: 1928 Physical Review titles already carry "Heisenberg's
Indetermination Principle" (Kennard) and "The Principle of Uncertainty in
Weyl's System" (Breit), Eddington's 1928 book made "principle of
indeterminacy" famous in
English, Condon and Robertson write "uncertainty principle" in 1929, and by
1930 the German schoolbook word, Unschärfe, blur, is settled enough that
Schrödinger's refinement paper wears "Unschärfeprinzip" in its title as a
given. So the pop label has no single author: the estimate is Heisenberg's,
the theorem is Kennard's and Robertson's, and the name is sediment, laid
down by many hands, none of them the honoree's, on an inequality he neither
wrote nor named.
This is a label-history shape the series has not met before: not a
complaint transported as a feature, not a single freezing window, but a
fusion, two different statements (the width theorem and the kick story)
welded under an authorless, accreted name.

One irony is worth keeping: of the two languages, the German relabel is the
more accurate one. Unschärfe, blur, is a property of an image; uncertainty is
a property of a knower. The German word points at the object; the English
word points at us.

---

## 3. Where the label breaks

**"You cannot KNOW both."** The knowing imports a knower, and behind the
knower a pair of sharp values waiting to be known, as if the electron had an
exact position and an exact momentum and the world declined to show us both
ledger entries. The width theorem needs none of that furniture. ΔX and ΔP are
widths of the state, measurable as plainly as the width of a spectral line,
and the theorem says the shape cannot be narrow twice. A violin string does
not "know" its pitch to infinite precision and hide it; a short pluck simply
has a broad spectrum. The theorem is about what the object is, not about
what anyone can find out.

**"Measuring one disturbs the other."** True physics, wrong theorem. The kick
story is Heisenberg's own 1927 canvas, and at his stance it was an honest
heuristic; but the inequality the textbooks print under it is Kennard's and
Robertson's, and theirs is a preparation statement that holds with no second
measurement anywhere in sight. The two statements fused in transport, and the
fusion is not harmless: the disturbance version was only made precise in
2003-2013, its correct form depends on how error is defined, and that dispute
is still live (Section 1). The pop label sells the settled theorem with the
unsettled story's plot.

**"Heisenberg's uncertainty principle."** Every word of the name wandered.
The inequality is Kennard's and Robertson's; "uncertainty" and "principle"
are the transport channel's, attached by several hands within two years,
none of them his. The man the label honors wrote a tilde estimate and a
kick story, corrected the kick story himself in the same paper, and in 1927
used different words for the whole thing. The label is a monument with the
wrong name engraved on all three lines.

One precision, owed to the lens's own discipline: what survives,
perspective-bound, is that both trades are real. You genuinely cannot prepare
an ensemble narrow in Z and narrow in X; the full bit of Maassen-Uffink is
paid every time. And a real measurement genuinely costs; the price of asking
is the [fifth entry's](DOUBLE_SLIT_TRANSLATED.md) whole subject. Kill the
doubt, keep the width; kill the
hidden ledger, keep the trade.

---

## 4. The translation (the exact part)

**Sharpness is an angle fact, so the trade is a geometry fact.** The second
entry established the ingredient: every pure state is a basis state of some
basis, and being "spread" is frame-relative; |0⟩ is sharp in Z and an even
superposition in X ([Superposition Translated](SUPERPOSITION_TRANSLATED.md)).
The uncertainty relation is what happens when two frames are held up at once:
one vector cannot sit on the pole of two different axes. Sharp along Z means
poised exactly between the two X-answers, and the Maassen-Uffink bit of
Section 1 is the standard account's exact pricing of that geometry. The fifth
entry's engine, the which-path question and the screen question that
anticommute, answering one at full strength silencing the other
([Double Slit Translated](DOUBLE_SLIT_TRANSLATED.md) §4), is one working
instance of this family; V² + D² ≤ 1 (Jaeger, Shimony and Vaidman, PRA 51,
54, 1995; Englert, PRL 77, 2154, 1996, independently) is its trade law, and
even the bridge from such duality pairs to the entropic relations of
Section 1 is a theorem of theirs (Coles, Kaniewski and Wehner, Nat. Commun.
5, 5814, 2014). All of that is the standard account's own
mathematics, cited, not re-derived. (The anticommutation here is that of
Hilbert-space observables, σ_x σ_z = −σ_z σ_x; the level fence for that
word, and the one time this house misread it one level up, is in
Section 6.)

**Our own house holds the sending-side face of the same shape.** The Tier 1
routing law (`HeldLetterRoutingClaim`, parents the
[Absorption Theorem](../proofs/PROOF_ABSORPTION_THEOREM.md) and the
[Klein V₄ proof](../proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md);
first published in [Dephasing Translated](DEPHASING_TRANSLATED.md)) prices
every Pauli string S under a held letter P:

    L_P(S) = −2γ · n_anti(S, P) · S,

where n_anti counts the sites whose letter stands at odds with the held one
(neither I nor P). One shared eigenbasis, three price lists: each held letter
exempts exactly its own cell, the 2^N strings built from I and that letter
alone, and bills everything else by how many of its sites disagree. The
letter the arriving light holds sets which structures it passes through
and which it bills.

And now the intersection, the one exact fact in this house with the
uncertainty shape, held live by the witness and never yet spent as a reading:
**the cell free under Z, the cell free under X, and the cell free under Y
overlap in exactly one string, the identity. Only nothing is free
everywhere.** Every nontrivial structure, every string that says anything at
all, stands at odds with at least two of the three letters, and pays the
moment either is the held one. Sharp in one light is billed in another; the
only object
exempt from every light is the object with no content.

Two fences, so this stays exactly what it is. First, the fact is elementary:
a non-identity string carries some letter Q ≠ I at some site, and Q sits in
the free set {I, P} of exactly one held letter, its own; the triviality of
the three-way intersection is the Pauli algebra's trivial center wearing
rates (the algebra's center, the scalars; the Pauli group's own center is
the phases). It is not Robertson's theorem and does not imply it, and the
disanalogy is part of the fence: the three price lists commute and share
one eigenbasis, so nothing here is a trade-off between incompatible
questions; every string carries all three prices at once, and what is
unique is the joint free cell, not a joint sharpness. The width relations
of Section 1 price spreads of states, this law prices decay of structures,
and the kinship is shape, not derivation. Second, the Absorption Theorem's
own fence rides along: −2γ·n_anti is the decay rate outright where the
Hamiltonian does not mix the sectors the string lives in (in the empty
world, always); where it does, the exact statement is an
eigenvector-weighted average, and the closed form is the initial slope, not
the rate.

**You can read the three price lists live.** The typed witness recomputes the
whole law at inspect time (N = 3, γ = 0.05: all 192 letter-string pairs
against the closed form, residual at machine precision; the from-below gate
asserts < 10⁻¹², and the printed digits are BLAS-dependent, so none are
pinned here). The Z and X rows below are the witness's own printout, as are
all three free counts; the two Y rate cells follow from the same closed
form:

| held letter | free cell        | strings free | Z⊗Z⊗Z pays | X⊗X⊗X pays |
|-------------|------------------|--------------|------------|------------|
| Z           | {I,Z} per site   | 8            | 0          | −0.3       |
| X           | {I,X} per site   | 8            | −0.3       | 0          |
| Y           | {I,Y} per site   | 8            | −0.3       | −0.3       |

Free under all three letters: 1 string of 64, the identity alone. Repriced by
the Z→X letter swap: 44 of 64. The −0.3 is −2γN at γ = 0.05, N = 3: full
disagreement, every site at odds with the held letter (the X⊗X⊗X cell under
held X may print as −0 in the run, the same zero). One command shows it:
`dotnet run --project compute/RCPsiSquared.Cli -- inspect --root label`.

**What the letter swap does and does not do.** Switching the held letter
never changes the price ladder itself: relabeling letters is a bijection of
the 4^N strings that carries n_anti(·, P) to n_anti(·, P′), so the multiset
of rates is identical by integer arithmetic, no eigensolver anywhere; the
swap re-assigns which string sits at which price and moves nothing else. The
operator-level identity behind this, the three dissipator diagonals
conjugate to each other, one letter-orbit, is checked in
[The Three Diagonals](../THE_THREE_DIAGONALS.md) (the typed claim is Tier 1
derived; the conjugation check there is float at 10⁻¹⁰ because both
conjugators carry a 1/√2, so no exact route exists for it). At the
Hilbert-space level the letter swaps are ordinary unitaries, the Hadamard
for Z↔X and a quarter turn about x for Z↔Y, and both carry the dissipator
over at machine precision (the witness prints deviations of order 10⁻¹⁶ at
N = 3 and gates them below 10⁻¹²; the same 1/√2 rules an exact route out
here too). A different group lives one level up: the operator-space
Klein V₄ of the
[Klein V₄ proof](../proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md)
(Tier 1), whose non-identity elements D (the transpose) and Q_yx (a basis
permutation) do not lift to Hilbert-space unitaries at all; within that V₄
only the Hadamard lift Q_zx rotates the dephasing axis, while D and Q_yx
intertwine the palindromic mirrors Π but leave the dissipator on its axis.
The series' one law, the −2γk bill whose k = 1 end is the pattern,
k = 2 the pair, k = N the cat, gets its fourth sighting here (its k is this
law's n_anti, the same count of sites at odds with the held letter, worn by
coherences there and by strings here), and it is not another k: it is the
letter that defines k, turned.

---

## 5. The readings (labeled as readings)

**The width theorem has a classical ancestor, and the ancestry is exactly
half the story.** "A nonzero function and its Fourier transform cannot both
be sharply localized" is a meta-theorem of harmonic analysis (Folland and
Sitaram, J. Fourier Anal. Appl. 3, 207, 1997), the same mathematics that
makes a short pluck broad-spectrum and a short pulse wide-band; Wiener
lectured on its signal-analysis face in Göttingen in 1925 (the surviving
record is the nontechnical account in his autobiography), and Kennard and
Weyl (crediting Pauli) proved the quantum case in 1927-28. For position and
momentum, the quantum relation is this wave mathematics plus de Broglie's
p = ℏk, and on the t-axis (the theory chapter's test: hold the physics
fixed and slide only t, [Labels Translated](LABELS_TRANSLATED.md) §3) that
component transports perfectly: it was true of sound before quantum
mechanics and will outlive every interpretation. It even transports into
the finite world: the Z-question and the X-question on one qubit are a
conjugate pair under the two-point Fourier transform, which is exactly the
Hadamard of Section 4, and the Maassen-Uffink bound for the qubit's
unbiased pair (each state of one basis spread evenly over the other) is the
finite Fourier form of the same meta-theorem. What does not transport is
the instrument: in the finite world the sharp price is entropic (the
product-of-widths form survives but goes trivial on a whole great circle,
Section 1), and for a general observable pair Robertson's bound needs only
a commutator, no Fourier structure anywhere. The label's abstract component
(conjugate widths trade) travels; its indexical component (the microscope,
the kick, the knower at the eyepiece) fossilized within three years of the
painting, when the theorem moved from the asking to the shape.

**Seen from the sending end** (a reading in this house's own register): the
three price lists say that a letter is not a property a structure has but a
side it shows. The chain stands in light, the light holds one letter, and
what the light passes through is whatever agrees with that letter; turn the
letter and the same structure that rode free now stands at odds and is
billed. "Only nothing is free everywhere" is then the sending-side twin of
the width theorem's moral: there is no structure that agrees with every
light, just as there is no state that is sharp on every axis. Content means
taking sides, and taking sides means there is a light that bills you.

---

## 6. An honest note on our own house

Three of this entry's words are already taken in this repository, and each
collision is disarmed here rather than discovered later. **"Heisenberg"** in
this house means the coupling, J σ·σ, in thousands of places, and
[Heisenberg Reloaded](../../hypotheses/HEISENBERG_RELOADED.md) is a sketch
about that coupling's origin, not about the man; when this entry says
Heisenberg it means the person, a use the repository has hardly needed
before. **"Complementarity"** already has two in-house lives: the
[glossary](../GLOSSARY.md) uses the word for the mirror's weight involution,
XY-weight k ↦ N−k, and the Absorption Theorem's typed claim uses it for the
partner-light pairing light_s + light_f = N; neither has anything to do
with Bohr, so where this entry needs the 1927 doctrine it says "Bohr's
complementarity." Two symbols collide as well, and are named so they cannot
route wrong: the γ in "γ-ray microscope" is a photon energy band, not this
repository's dephasing rate (which Section 4 sets to 0.05), and the Q in
Section 4's fence is a Pauli letter, neither the coupling knob Q = J/γ₀ nor
the dephasing diagonal Q_P. And **"uncertainty"** itself appears in
this repository mostly as fit statistics, the ±σ of a measured
number; the label this entry translates never had an in-house life at all,
as the sweep recorded at the top of this document found store by store.
This is the series' first entry whose subject the repository holds no
result on, only cousins.

One more level fence, because Schrödinger's 1930 refinement sets an
anticommutator, the symmetrized covariance, beside the commutator square,
and this house has learned to flinch at that operator: the anticommutator
is legitimate there, on Hilbert-space observables. The entry in
[Caught Errors](../CAUGHT_ERRORS.md) (2026-06-22) concerns superoperators,
one level up, where the commutator is the operative object and the
anticommutator once seduced a proof into a non-sequitur. Same word, two
levels; this entry keeps them apart by naming the level at every use.

---

## 7. The in-repo cousins

**The true kin.** In the wave-mechanical sense, the conjugate pair of
Section 1 lives natively in the walk:
[the dephasing front renewal proof](../proofs/PROOF_DEPHASING_FRONT_RENEWAL.md)
(Tier 1) expands the single-excitation walk on the infinite chain in lattice
waves with "the momentum p conjugate to the site index," and the front's
shape is the interference of those waves. The conjugate pair is no import
here: the proof's own expansion runs in exactly those waves, and the front
disperses as its momentum content demands. (The finite
Fourier face of the same kinship already appeared in Section 5: the Z-X
pair under the Hadamard.)

**The look-alike, disarmed.** The repository also holds a Tier 1 resolution
limit that wears the optics vocabulary: the band-edge conditioning floor
(`BandEdgeResolutionLimitClaim`, live at `inspect --root resolution`), where
a staggered defect is ~N times harder to localize than a band-edge one, "the
q = π detail at the resolution cutoff, the diffraction limit." It is a real
limit and a real optics rhyme, but it is the other kind of limit: the
conditioning of one linear map, its floor set by the Dirichlet sine modes
and its worst direction the staggered q = π detail at the resolution
cutoff, with no pair of non-commuting readings and no state paying a
product of spreads. Setting it beside the width theorem sharpens both: an
uncertainty relation needs two questions that cannot be asked sharply at
once; a resolution floor needs only one question asked through a bad lens.

---

## The right label

The relation: a width theorem, proven about shapes, in 1927 by Kennard for
position and momentum, in 1929 by Robertson for every pair. The
spread: a property of the object, as physical as the linewidth of a lamp,
not a gap in anyone's ledger. The disturbance: real, younger, and a separate
theorem family whose right definition is still in dispute. The name: the
transport channel's own sediment, settled on a man who did not choose it.
And the in-house shape beside it: three
price lists over one eigenbasis, each letter exempting only its own cell,
only nothing free everywhere. The sentence the whole entry compresses to:
**the uncertainty principle is not a limit on what can be known about sharp
values; it is the statement that a shape narrow along one question is wide
along the conjugate one, before anyone asks.**

Stamped: this canvas is ours, painted 2026, and its mortal component is
already visible from here. "Width" imports an image, a picture drawn in
space, and our whole geometric register, angles, poles, sides, may read one
day the way the microscope's eyepiece reads now: an author's instrument
mistaken for the theorem. The previous entries' stamps have already begun to
fire inside this house's own lifetime, one after twenty-eight days, one
after thirty-four; this one is offered with the same expectation. The future
reader receives the canvas with its date.

The closure, then: Heisenberg's canvas, the cost of asking, painted at the
stance of an asker of trajectories, with its correction slip attached in the
same paper; Kennard's canvas, the width of shapes, the theorem the label
actually carries; the untouched algebra, spreads and commutators, bits and
letters, priced exactly then and priced exactly now; and our canvas, the
width in the object and the letter in the light, where only nothing is free
everywhere. One object. The electron never hid its momentum, because there
were never two sharp numbers to hide; the measurement never needed to
disturb, because the widths were in the shape before any asking; and the
doubt the label sells was never in the physics, because the physics holds
no knower to be uncertain. What the label has been guarding all along is a
geometry: to be something is to take sides, and every side has a light that
bills it.
