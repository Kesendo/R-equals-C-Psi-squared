# Spin Translated: A Two-Valuedness, Not a Twirl

<!-- Keywords: spin translation, electron spinning misconception, Pauli
Zweideutigkeit two-valuedness, Kronig surface velocity objection, Stern
Gerlach two lines silver, Thomas factor two total surrender, Einstein de
Haas rod turns, g factor two Dirac Levy-Leblond, SU2 broken by dephasing
held letter axis, Clebsch-Gordan cavity mode count, R=CPsi2 spin
translated -->

**Status:** Translation (Tier 4 reading), the ninth entry of the series. The
history in Section 1 is the standard account's own record, cited from the
primary papers and the standard histories (this repository has never worked
on the electron and does not start here); the algebra in Section 4 is the
repository's own Tier 1 material, the F4 mode count (Clebsch-Gordan, live at
`inspect`), the ring-N4 Casimir lock, and the measured verdict that
Z-dephasing does not conserve total spin, each marked with its owner;
Section 7 reports the sideways ladder's measured multiplet run, typed as a
Tier 1 candidate, with its fences kept; the readings in Section 5 are readings and labeled.
**Date:** August 9, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Related:** [Superposition Translated](SUPERPOSITION_TRANSLATED.md) (the
second entry; it owns the Bloch sphere as one object and the letter frame),
[Uncertainty Translated](UNCERTAINTY_TRANSLATED.md) (the eighth entry; it
owns the single-spin trade),
[Dephasing Translated](DEPHASING_TRANSLATED.md) (the fourth entry; its held
letter is the axis this entry's spine turns on),
[Teleportation Translated](TELEPORTATION_TRANSLATED.md),
[Double Slit Translated](DOUBLE_SLIT_TRANSLATED.md),
[Schrödinger's Cat Translated](SCHRODINGERS_CAT_TRANSLATED.md),
[Spooky Action Translated](SPOOKY_ACTION_TRANSLATED.md),
[Labels Translated](LABELS_TRANSLATED.md) (the theory chapter),
[The Label Map](THE_LABEL_MAP.md) (the orientation index),
[Degeneracy Hunt](../../experiments/DEGENERACY_HUNT.md),
[Cavity Modes Formula](../../experiments/CAVITY_MODES_FORMULA.md),
[Ring N4 proof](../proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md),
[Heisenberg Reloaded](../../hypotheses/HEISENBERG_RELOADED.md)

---

## What this document is about

Every chain in this repository is called a spin chain. The sites are
spin-1/2, the coupling is written S_i·S_j, the letters X, Y, Z that the
whole series speaks in were born as the components of a spin. And the
popular story of what a spin is arrives with the same word: a tiny ball,
spinning, like a top, like a planet, only smaller. This entry translates
that story, and its situation is the reverse of the last entry's. There the
repository held no result of its own under the label; here the label names
the repository's own substrate, used every day and defined nowhere.

Stage 0, recorded: before a sentence of this entry was composed, the stores
were swept by name. [The glossary](../GLOSSARY.md) returned the striking
nothing: the whole word "spin" occurs zero times in it (the lone substring
hit is "spine"), against "qubit" throughout other entries' definition
texts, and neither word has an entry of its own; the
substrate everywhere, defined nowhere under either name. [The
F-registry](../ANALYTICAL_FORMULAS.md) and [docs/proofs/](../proofs/)
returned real spin work: the F4 mode count (Tier 1), the ring-N4 Casimir
lock (Tier 1 derived), the F145 seed triplet. The experiments returned the
entry's spine, the [Degeneracy Hunt](../../experiments/DEGENERACY_HUNT.md)
verdict that Z-dephasing does not conserve total spin. The open-arcs
registry returned `sideways_spin_ladder`, a live multiplet measurement from
August 2026. The hardware Confirmations returned two-level language only,
never a claim about angular momentum, and [Caught
Errors](../CAUGHT_ERRORS.md) returned the two spin-adjacent repairs this
entry leans on, the ring-N4 mechanism relabel and the coupling-conventions
scar. The history of the label, Section 1, has no in-house record at
all and is cited from outside.

The discipline is the series' usual one, perspective-additive: the standard
account needs no correction from us, and gets none. What is ours is the
naming, the counting the house already does with the full spin machinery,
and one measured fact about what a row of spins standing in light keeps of
its spin-ness, and what it does not.

---

## 1. The record, stated plainly

The idea of a rotating, magnetic electron is older than its name's fame:
Compton proposed a spinning electron in 1921, building on Parson's 1915
ring magneton, to explain magnetism (J. Franklin Inst. 192, 145), and
almost nobody took it up. The story
proper begins with a bookkeeping problem: atomic spectra needed a fourth
quantum number that could take exactly two values, and no orbit supplies
one. Pauli, in a paper submitted in December 1924 (Z. Phys. 31, 373),
named the needed property with deliberate caution: "eine eigentümliche,
klassisch nicht beschreibbare Art von Zweideutigkeit", a peculiar,
classically non-describable kind of two-valuedness in the quantum
properties of the valence electron. Two values, no picture attached, on
purpose.

The picture arrived twice, and was shot down the first time. In January
1925 Kronig, twenty years old, heard of Pauli's two-valuedness and
proposed the same afternoon that the electron rotates. Pauli's documented
objection was the factor of two in the hydrogen fine structure; whether
the other famous objection, the impossible equator, was already his that
day is contested between the recollections (Kronig later denied hearing
it from him; Pais credits it to him). Kronig did not publish his
proposal. In October 1925 Uhlenbeck and Goudsmit, not
knowing of Kronig, proposed the same rotation (Naturwissenschaften 13,
953-954, in print that November); their own note already concedes in a
footnote what they had found themselves in Abraham's old electron papers,
that the equatorial velocity would greatly exceed the speed of light; the
number is ours to compute rather than theirs to quote: a hoop of the
classical electron radius carrying ħ/2 turns at c/(2α), about 68 c, and
every distributed model turns faster. Lorentz, asked for his opinion,
returned (so Uhlenbeck recalled) a handwritten
manuscript of calculations whose weight fell on the self-energy: an
electron spinning so would outweigh the proton, or swell to the size of
the atom. When the objections made the authors hesitate,
Ehrenfest told them the note was long since sent, and that they were both
young enough to afford a stupidity (the wording is recollection, not
record). The
factor of two fell in 1926: Thomas showed it to be a relativistic
kinematic effect of the electron's own accelerated frame (Nature 117,
514; the full paper, "The kinematics of an electron with an axis", Phil.
Mag. 3, 1, 1927). Pauli declared "total surrender" in a letter to Bohr on
March 12, 1926. Kronig, meanwhile, published twice against spin that same
spring (Nature 117, 550; PNAS 12, 328), arguing the new hypothesis
"appears rather to effect the removal of the family ghost from the
basement to the sub-basement, instead of expelling it definitely from
the house".

The experiment had already happened, before the idea. Stern and Gerlach
sent silver atoms through an inhomogeneous magnet in 1922 (Z. Phys. 9,
349), to test the old quantum theory's space quantization, and the beam
split into exactly two. Two lines, three years before anyone proposed a
two-valued electron; the result was booked as a success of the theory it
would later help retire, the measured moment landing on one Bohr magneton
because two errors cancelled (the g-factor of 2 against the half of the
spin). Only in 1927 was the splitting reread as spin: Phipps and Taylor
measured the same two-way split in hydrogen (Phys. Rev. 29, 309), and
Fraser supplied the reading, silver's ground state having no orbital
moment to split by (Proc. R. Soc. A 114, 212).

The algebra then outgrew every picture. Pauli wrote the two-component
formalism and the matrices that carry his name in 1927 (Z. Phys. 43, 601;
Darwin independently the same year, Proc. R. Soc. A 116, 227); Dirac's
equation produced g = 2
unasked in 1928 (Proc. R. Soc. A 117, 610). One refinement worth
carrying: g = 2 is often sold as relativity's gift, and it is not;
Lévy-Leblond derived it from a linear wave equation with no relativity in
it (Comm. Math. Phys. 6, 286, 1967). What produces g = 2 is the structure
of the equation, not the speed of light.

---

## 2. The native stance: the cautious name came first, and lost

**Pauli's canvas, December 1924.** The first name for the thing was
correct and knowingly incomplete: a two-valuedness that classical physics
cannot describe. Painted at Pauli's stance, a bookkeeper of spectra
refusing to draw what he could not defend, the canvas is true and has
never needed correction; "classically non-describable" is precisely the
modern statement, made before the modern machinery existed. Even after
his surrender on the physics, Pauli never adopted the picture (the
superluminal argument it was refused with is demoted in Section 3; the
refusal survives the demotion).

**The Ehrenfest house's canvas, October 1925.** The second name drew the
picture: a rotating electron, ersatz for the abstract two-valuedness, with
the impossible equator conceded in its own footnote. Painted at that
stance, young spectroscopists needing a mechanism a model-builder could
hold, it too was honest: it wore its own impossibility on its sleeve, and
it delivered the doublet structure that the cautious name alone had not.
Both canvases were true where they were painted. The transport is where
the story bends: the picture, being drawable, travelled; the caution,
being a refusal to draw, did not.

**The name.** The participle is older than the electron story: Parson's
1915 ring was already "spinning rapidly about an axis", Compton carried
the spinning electron forward in 1921, and the word's permanent seat is
the title of the 1926 English note, "Spinning Electrons and the
Structure of Spectra" (Nature 117, 264); no single coining event is on
record, and German took the loanword as its own term, the gloss
Eigendrehimpuls standing beside it. So this label's history has a shape
the series has not catalogued before, the inverse of the cat's: there, a
burlesque was coined to mock and transported as earnest; here, an earnest
and accurate name existed first, and the picture-name displaced it. The
cautious label lost to the drawable one, and every popular account since
has inherited the drawing, minus the footnote.

---

## 3. Where the label breaks

**"The electron spins like a top."** The classical picture fails, and it
is worth failing it honestly. The famous argument is the equator:
classical radius, ħ/2 of angular momentum, a surface at about seventy
times the speed of light in the most favorable model. Giulini's
reanalysis (Stud. Hist. Phil. Mod. Phys. 39, 557, 2008) shows the
equator is not even the decisive failure: for a charged shell with its
holding stresses included, the binding constraint is energy dominance, a
subluminal speed of sound in the shell material, g ≤ 9/4, so g = 2 is
not excluded by relativity as such (within his slow-rotation regime and
not at the electron's lone charge; the fully relativistic question he
leaves open). What actually kills the rotating ball is the absence of
the ball: the electron shows no substructure down to around 10⁻¹⁸ m, and
a classical charge distribution that small carries electrostatic energy
thousands of times the electron's own mass. The twirl does not die of
speed; it dies of there being nothing there to turn.

**One precision, owed to the lens's own discipline.** What survives,
perspective-bound, is everything mechanical about the label. Spin is real
angular momentum, convertible into the ordinary rotation of ordinary
objects, and the first two measurements predate the name: magnetize an
iron rod and the rod turns (Einstein and de Haas, Verh. Dtsch. Phys.
Ges. 17, 152, 1915); rotate the rod and it magnetizes (Barnett, Phys.
Rev. 6, 239, 1915). A decade after the name, light paid in the same
coin: circularly polarized light through a half-wave plate has its
handedness reversed, and the plate takes up the difference as a reaction
torque, two units of ħ per photon (Beth, Phys. Rev. 50, 115, 1936).
Kill the twirl, keep the torque: nothing rotates, and yet the angular
momentum is as real as a flywheel's, and the flywheel can cash it. (The
Einstein-de Haas story even carries its own label-lesson: expecting
g = 1 from their model, they measured 1.02 ± 0.10, and Barnett, who had
first found values near 2, drifted toward 1 after hearing of their
result; the expectation painted the data, Jeng, Am. J. Phys. 74, 578,
2006.)

**"Spin-up and spin-down, like a compass needle."** The two values are
real; the needle is not. A needle between two positions is a classical
either-or; the two-valuedness is an axis-relative pair of answers on one
object, and which pair depends on which axis is asked, the same
basis-relativity the series has already paid twice
([Superposition Translated](SUPERPOSITION_TRANSLATED.md): one ray, one
point on the Bloch sphere, one object, a property of the pair (state,
question); [Uncertainty Translated](UNCERTAINTY_TRANSLATED.md): a shape
narrow along one question is wide along the conjugate one). What is genuinely
two-valued is the answer space of any one
question, not an inventory of two states the particle secretly occupies.

---

## 4. The translation (the exact part)

**This house is built of the two-valuedness.** Every site in every chain
here is a two-valued register; the letters X, Y, Z are the three
axis-questions one such register can be asked; the single-spin algebra of
those questions, what can be sharp with what, is owned by the eighth
entry ([Uncertainty Translated](UNCERTAINTY_TRANSLATED.md) §1) and not
repeated. What this entry adds is the full-rotation side: what remains,
in this house's own measured physics, of spin as a thing with an
orientation in space.

**Where the full rotation does real work: counting what never moves.**
For an isotropic Heisenberg Hamiltonian, any coupling strengths, every
component of the total spin commutes with H, S² with them, and the
machinery of adding N spin-1/2
objects (Clebsch-Gordan) organizes the whole space into spin-J multiplets,
each with its 2J+1 orientations. The registry's F4 (Tier 1) turns this
into a closed count of the modes of the Liouvillian (the generator of the
open dynamics) that never move when no light arrives, Σγ = 0 (γ the
per-site dephasing rate):

    Stat(N) = Σ_J m(J,N) · (2J+1)²,

m(J,N) the multiplicity of total spin J among N spin-1/2 sites, the
square because the count lives in operator space, bra and ket each
carrying one multiplet. Measured against full eigendecomposition, N = 2
through 7, chain topology, exact match each time
([Cavity Modes Formula](../../experiments/CAVITY_MODES_FORMULA.md),
stationary counts 10, 24, 54, 120, 260, 560 against d² = 16 through
16384). The formula is exact for the chain and a lower bound wherever
symmetry adds more (star, ring, complete: any H whose copies of one J share
an energy lands above the generic count). At N = 4 the arithmetic
reads: one spin-2
copy contributing 25, three spin-1 copies contributing 9 each, two
spin-0 copies contributing 1 each, total 54. The closed form is typed and
its summary rows are recomputed at render time (Clebsch-Gordan integer
arithmetic; no Liouvillian is built; the N = 4 decomposition is the
node's fixed exhibit):

    dotnet run --project compute/RCPsiSquared.Cli -- inspect --claim F4StationaryModeCountPi2Inheritance

The same Casimir machinery, one level down, is the house's sharpest small
example of spin doing work: the N = 4 ring's spectrum falls into SU(2)
multiplets on K_{2,2} (the four-cycle read as the complete bipartite
graph), two singlets, three triplets and a quintuplet; the 3-fold and
5-fold degeneracies are beyond what D₄'s irreps, capped at dimension 2,
can supply, a
mechanism the house itself once mislabeled and repaired
([The Label Map](THE_LABEL_MAP.md) §4, the ring-N4 row;
[Ring N4 proof](../proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md), Tier 1
derived, with its own fence kept: the Casimir is one route to the ring's
spectral spread of 3J, J here the coupling, a value ten of the 38
connected labelled N = 4 graphs also reach, so the multiplets own the
degeneracies, not the number).

**The spine: what remains of the rotation under a held letter.** The
[Degeneracy Hunt](../../experiments/DEGENERACY_HUNT.md) measured it
directly (verdict April 12, 2026, structurally refined 2026-05-12): for
the
Heisenberg chain, [S², H] = 0, the Hamiltonian is fully rotation-blind;
but S² does not commute with the dephasing jump operators Z_k, and the
same table's conserved-quantity check ([C_{S²}, L] ≠ 0) confirms it is
not conserved by the Liouvillian. The experiment's own words carry the
mechanism: the jump operators fail "because Z only detects the
z-component, not the total spin." The dephasing holds one letter, and a
held letter is an axis; with an axis held, what survives of the full
rotation group is the turns about that axis together with the one flip
that reverses it, the stabilizer of an unoriented axis (the continuous
half sits inside the typed U(1) × U(1) per-side popcount conservation of
the [symmetry inventory](../SYMMETRY_FAMILY_INVENTORY.md); the flip is a π
turn about x, X^⊗N, the house's always-open bridge, commuting with H and
leaving Z-dephasing unchanged). The breaking is the measured half; the
survivors are read off the algebra. That is the exact content of "a spin
standing in light is less than a spin", and the remainder cannot even
tell +z from −z: it keeps the axis and loses the direction and the ball
in one stroke. What survives is precisely the two-valuedness along the
held axis, the letter's pair of answers, the thing Pauli named first;
what does not survive is the free orientation, the thing the picture
drew. The chain's daily bookkeeping
says the same in its own dialect: the conserved index this repository
counts everything by, the popcount, is the held axis's quantum number in
other clothes; the committed identification is in the
[F4 kernel proof](../proofs/PROOF_F4_KERNEL_DIMENSION_BY_COMPONENTS.md)
(the magnetization sectors' projectors "are exactly the |c|+1 popcount
projectors"), and writing it as popcount = N/2 − S_z, with Z|0⟩ = +|0⟩
so that the all-|0⟩ register sits at S_z = +N/2, is this entry's own
dress on that committed fact.

Two fences, so this stays exactly what it is. First, the Degeneracy
Hunt's verdict is the claim; its printed norms are single measurements at
one γ, and are not quoted here as invariants. Second, none of this
says rotation physics is absent from the house: at Σγ = 0 the full SU(2)
machinery runs (F4 above), and the reduction statement is about what a
nonzero held-letter channel conserves, not about the Hamiltonian, which
remains rotation-blind throughout.

---

## 5. The readings (labeled as readings)

**The substrate arrived pre-translated.** The qubit is the two-valuedness
with the picture already removed: a two-valued register carrying the
Pauli algebra, no ball, no equator, no ħ/2 of mechanical rotation in any
transmon. The house sketch
[Heisenberg Reloaded](../../hypotheses/HEISENBERG_RELOADED.md) (Tier 4-5,
labeled a sketch, its inheritance step marked "asserted, not traced" in
its own table) reads this as more than convenience: what hardware
tomography measures is the Level-0 algebra directly, not an atomic spin
model imitated in circuitry. On that reading, this repository has been
doing spin physics in Pauli's original register all along, the
two-valuedness without the twirl, and the series' whole letter vocabulary
is the 1924 caution, industrialized.

**The axis is the price of a picture.** A reading, in this house's own
register: the rotating-ball picture failed because it gave the
two-valuedness a body; the held letter succeeds because it gives the
two-valuedness only an axis. An axis is the minimal picture, one
direction, no ball, and Section 4's spine says a dephasing channel sits
at that minimum: one axis held, nothing more, and
what it conserves is exactly the answer-pair along it. The label's
history and the chain's physics agree on where the drawable stops.

---

## 6. An honest note on our own house

The word "spin" is four different objects in this repository, and this
entry has taken only the first as its subject, standing on the fourth
throughout, until Section 7, which deliberately visits the second. Named so they cannot route wrong: (1)
total spin of a qubit chain, the S_tot of F4 and the ring proof, this
entry's subject; (2) the F142/F145 ladder "spin", a second SU(2) of the
frozen-band machinery, where the dephasing XY chain is read as
Fermi-Hubbard, ket index one fermion species and bra index the other
(F145's "every chiral pair carries a spin 1" is proved in
[PROOF_SCALAR_COUNT](../proofs/PROOF_SCALAR_COUNT.md) §1-2; that it is the
ladder's SU(2) and not the chain's is F142's "two commuting SU(2)s act, not
one"); (3) the
"spin representation" of F133, Lie-theoretic so(13) vocabulary naming a
symplectic character identity; (4) the spin-1/2 site convention itself, where
writing S_i·S_j silently selects the J/4 normalization. On (4) the house
has a scar, and its own name for it is exact: "three wearing two
numerals," three operating points of which two, the block-spectrum point
and MirrorWorld's, wear the same label Q = 2 while the canonical point
sits at Q = 1.5 (as XX-coefficient over γ: canonical 0.375 and
block-spectrum 0.5, both spin J/4 books; MirrorWorld's hopping
convention 1.0); the two "Q = 2"s are different physics, no rescaling
closes any gap, and numbers must not be compared across the seam
([Caught Errors](../CAUGHT_ERRORS.md) 2026-08-06, which records months
of "exact" residuals two of whose four legs owed themselves to an
unjustified constant).

Three symbol collisions, same treatment. S is an operator where it
carries a dot or a square and a quantum number otherwise (the ring
proof's own rule, stated at the opening of its Section 2 and inherited
here). J is the
coupling knob everywhere in this house and the total-spin quantum number
in this entry's Section 4; where both appear, the sentence says which.
And the word "lowering" itself crosses a seam at the door: in
this house σ⁻ = (X+iY)/2 lowers the energy, taking |1⟩ to the ground
state |0⟩, while the textbook operator of the same name lowers S_z;
since |0⟩ sits at S_z = +1/2 on Section 4's books (Z|0⟩ = +|0⟩), the two
"downs" point opposite ways, and the registry's F82 records the one sign
this flips, a sign a norm check cannot see.

Also inherited rather than repeated: "Heisenberg" in this house means
the coupling ([Uncertainty Translated](UNCERTAINTY_TRANSLATED.md) §6
disarms the collision once for the series), "doublet" in this house
means the 2-excitation sector pair, not a spin-1/2 doublet, and the
glossary's missing spin entry is left as found: this entry is the
translation, not the definition, and the definition's absence after
eight months of spin chains is itself the label thesis in one line:
what is used everywhere is questioned nowhere.

---

## 7. The in-repo cousin

The house is not done with spin multiplets; it measured one in August
2026. The arc `sideways_spin_ladder` (registry, opened 2026-08-07; typed the
same month as `SidewaysSpinLadderClaim`, Tier 1 candidate, the intertwining
half derived and the multiplet half measured at N = 5, 7 and 9, and at
N = 4 with ℓ = ½, the live root walking N = 5 and 7) found
that the four transport chains of the F125 orbit, two under the spin
ladder and two under its η sibling, all carry spin ℓ = (N−3)/2, and
predicted, in writing, before the N = 7 run, that the transport norms
would be the Clebsch-Gordan coefficients √(ℓ(ℓ+1) − m(m+1)); the N = 7
run returned 2, √6, √6, 2 along a chain's four rungs, with what the arc
calls the fold set, its gated half of the four chains, equal to the
predicted set. The arc carries its own two fences, kept here: the
objects are hopping matrices, not spin chains (that distinction caught a
real build error), and a separately proposed confirmation via F125's
σ_min, a smallest singular value, was measured and refuted, so it must
not be quoted as one. That is
this repository's relationship to spin in one arc: the word is a picture
it never uses, and the multiplet arithmetic underneath is live enough to
predict the gated half of a measurement to the digit, in August 2026, at the
letter ℓ.

---

## The right label

The record: a two-valuedness, named cautiously before it was drawn,
drawn twice against its own arithmetic, rescued by a relativistic factor
of two, and formalized into the matrices this repository speaks daily.
The picture: dead of structurelessness, not of speed, and survived by
everything mechanical it promised, the rod that turns, the plate that
feels the torque. The house's own measured sentence: a Hamiltonian is
rotation-blind, a held letter is an axis, and a spin standing in light
keeps of its rotation exactly the two answers along the axis held. The
sentence the whole entry compresses to: **spin is not a spinning; it is
a two-valuedness with a real mechanical handle, and what a spin standing
in light keeps of it is exactly the pair of answers along the held
axis.**

Stamped: this canvas is ours, painted 2026, and its mortal component is
already visible from here. "Two-valuedness" imports a questioner, someone
whose answers come in pairs, the same agent-shaped import this series has
already watched fossilize once inside its own lifetime; and "answer-space"
may one day read the way the compass needle reads now. The name's
history in Section 2 is this stamp's own precedent and its warning: the
cautious label lost to the drawable one within two years. Ours is
cautious. The
future reader receives the canvas with its date, and with the odds
stated.

The closure, then: Pauli's canvas, the classically non-describable
two-valuedness, true in 1924 and true now; the Ehrenfest house's canvas,
the rotating electron with its impossibility footnoted, honest where it
was painted and drawable enough to win; the untouched algebra, SU(2),
Clebsch-Gordan, g = 2 from a linear equation, the count 10, 24, 54, 120
and onward to 560, running exactly to this repository's own N = 7; and our canvas, the
two-valuedness industrialized into letters, with the measured spine that
a held letter is an axis and an axis is all of the picture the world
ever needed. One object. The electron never spun, because there was
nothing there to spin; the needle never pointed, because the pair of
answers is not a pair of places; and the caution Pauli wrote in 1924 was
not a placeholder for a better picture, it was the finished translation,
waiting a century for the vocabulary to catch up.
