# Where the ceiling actually fails: the couplings at which one more mode freezes

**Status:** the measuring, deliberately unnamed. The existence is proved. The count on the
singlet space of the block (2,2) is exact through N = 8 and measured beyond, and it is the
full rung-2 count provably through N = 6, where the whole block was enumerated. The complete
count, over all rungs, is exact only at N = 5, where the whole band was enumerated at every
coupling at once. No F number, because
no law is known for any of the counts and because from N = 6 on no small closed form appears
for the couplings themselves, so the registry of closed forms is the wrong home for either.
(Not proved: that no radical expression exists. What is shown, at N = 6, is that the minimal
polynomial of the smallest coupling has degree 12 and that no small relation appears, which is
a weaker and sufficient reason; at N = 7 and 8 not even that is read.)
What is proved and what is measured are separated below and should stay that way.
**Date:** 2026-07-28
**Verification:** [`simulations/exceptional_couplings.py`](../simulations/exceptional_couplings.py)
(must print "exceptional couplings gate: ALL GREEN"; the run states its own check count, so no
number is carried here to drift). The default run reaches N = 6; `--deep` adds N = 7, the
whole-block enumeration at N = 6 and a numeric read of the N = 7 rung; `--slow` adds the two
exact counts that cost about half an hour, the rung-3 count at N = 7 and the rung-2 count at
N = 8. What remains a session measurement,
in no mode of the gate, is the whole band at N = 7, the rung-3 count at N = 8, and the
continuation of the rung-2 count past N = 8; those are recorded here as observations.
**Grew out of:** the step from large J to generic J. Section 5 of
[PROOF_FROZEN_BAND_SO4](../docs/proofs/PROOF_FROZEN_BAND_SO4.md) marked it asserted and not
argued and its Section 8 listed it first among the open items; it now carries it as
Proposition 5.3, and what that proposition leaves finite is this note.

This note is the companion measuring, not a standalone document: the band, the frozen root,
the two ladders Φ and S⁺ with their lowest-weight spaces, V₀ = ker(ad_h) and the transform
length M = N + 1 are all defined in
[PROOF_FROZEN_BAND_SO4](../docs/proofs/PROOF_FROZEN_BAND_SO4.md), Section 1, and are used here
in exactly that sense. Two of its words are spent twice in this repository and its Section 1
says which sense it means: **block** is the joint-popcount block (p, q), and **seed** is the
bottom rung of a ladder with the frozen modes that start there, never MirrorWorld's `Seed`,
which is a coupling-independent count and has nothing to do with this. A third word joins
them below: **field** is always the algebraic one, a number field such as ℚ(cos(π/M)), never
MirrorWorld's `Field`, which is the empty world running. Two shapes are used
here throughout: the **side lines** are the band's two off-diagonal lines |p − q| = 2, against
the diagonal line |p − q| = 0, and the **corner** is the block (1,1), where the ladders start.

Two more things to carry before the machinery. **Every coupling quoted in this note is in
units of the watching.** C is γ̄ times an integer matrix and A₀ carries no γ̄ at all, so
C + i·J·A₀ = γ̄·(C₁ + i·(J/γ̄)·A₀): the exceptional set is a set of RATIOS J/γ̄, and every number
tabulated below is that ratio at γ̄ = 1. At γ̄ = 2 the N = 5 point sits at J² = 6 rather than at
3/2, exactly (gate block E2). And **frozen** means the rate is undisturbable by the coupling,
not that the mode persists: the frozen root is λ = −4γ̄, so these modes decay like the rest of
the band, and what is frozen is their rate against J and not their amplitude against time.

## What this is about, before any of the machinery

The frozen band carries ⌊N/2⌋ undisturbable modes per block. That the count is at least
⌊N/2⌋ is proved at every coupling; that it is no more than ⌊N/2⌋ was proved only in the
large-J limit, and the missing step was expected to be a formality: a bound that holds for
all large J *ought* to hold at generic J, with the floor closing it to an equality. The
question nobody had asked is what the word "generic" was hiding. It is hiding something.
There are real couplings, at ordinary values, at which one more mode freezes.

## The finding

At the uniform watching point, on the diagonal blocks of the band, the frozen count is
⌊N/2⌋ for all but finitely many couplings, and the finite set is **not empty**. At an
exceptional coupling the count is one more, and exactly one more (rather than at least one
more) wherever the whole block has been enumerated, which is N = 5 and N = 6.

Two things make it worth a note rather than a footnote:

1. The extra mode is identified. At the couplings the block (2,2) can see, it is
   η-lowest-weight at rung ℓ = 2, so it occupies the diagonal rungs 2 through N − 2 and
   nothing else, and it is a **spin singlet**, where the band's own modes carry j_spin = 1.
   That block sees the rung-2 couplings and no others, because a block (p, p) carries only
   the rungs ℓ ≤ min(p, N − p); higher rungs carry couplings of their own, counted in
   [the rung the block cannot see](#the-rung-the-block-cannot-see) below. Whichever rung
   they sit on, the extra modes enlarge the diagonal line and never step sideways: the band
   statement |p − q| ∈ {0, 2} is untouched, and so is the count on the two side lines. Read as
   physics rather than as ladders, that says what the extra freezing IS: a coherence between
   two states of the same excitation number, never one bridging different particle numbers. That
   last sentence is exact over all couplings at N = 5, read at one coupling at N = 6, and
   read numerically over the whole band at N = 7; the sections below keep those apart.
2. The bridge that was open becomes elementary once the failure is known. Nullity of a
   matrix pencil is constant and minimal off a finite set, and larger on it, and the floor holds
   everywhere, so a **single** coupling with nullity ≤ ⌊N/2⌋ settles the generic value.
   Proposition 5.1 with F144 supplies one for every N ≥ 6 in the limit, and an exact GF(p)
   rank supplies a concrete finite one, J = 1, at every N the rank reaches. No semicontinuity
   argument is needed, and none would have delivered more, because "all but finitely many" is
   the true statement and not a technical hedge.

## Proved, and how

The existence, the counts, and everything stated for N = 5 are exact arithmetic over ℚ, with
no floating point anywhere; at N = 6 the whole-block count is exact too, under `--deep`. Two
things are numerical and are marked where they appear: the confirmation at N = 5 and 6, and
at N = 7 under `--deep`, that a root of the singlet polynomial makes the whole block singular,
which the exact enumeration supersedes at N = 5 and 6 and on which no N = 5 statement rests,
and the reading of WHICH blocks are raised at N = 6,
both at the rung-2 couplings and at the rung-3 ones of the section further down. Everything
past N = 7, and the whole band at N = 7, is numerical throughout and is kept in its own two
sections.

**The reformulation that makes it rational.** The frozen condition on a block is
(C + i·J·A₀)v = 0 with C the integer rate diagonal shifted by the frozen root and A₀ the
integer Hamiltonian part. Writing v = x + i·J·w turns it into

    C x − J²·A₀ w = 0,     A₀ x + C w = 0,

in which only J² appears and every entry is rational. So at a rational J² the whole kernel
question is a rank over ℚ.

**N = 5, by rank.** At γ uniform and **J² = 3/2** the blocks (2,2) and (3,3) carry 3 frozen
modes where ⌊5/2⌋ = 2, as an exact rank over ℚ; at J² = 1 and J² = 2 they carry 2. Every
other band block carries 2 at all three couplings, side lines included. Gate block E3.

**Every N, by determinant.** The extra mode lives in V := ker Ψ ∩ ker S⁺ ∩ ker S⁻ inside the
block, the η-lowest-weight spin singlets, where Ψ = Φ† lowers the η ladder that Φ raises and
S⁻ = (S⁺)† lowers the spin ladder. Singlet is meant in the sense of that second ladder, the
one the Hubbard disguise supplies, not of the chain's own spin. (V is not V₀ = ker(ad_h), the space the large-J
reduction of [ETA_CEILING_REDUCTION](ETA_CEILING_REDUCTION.md) lives on; the two are different
objects and the neighbouring letters are worth guarding against.) All three maps have integer matrices, and both
halves of L restrict to V exactly (Lemmas 2.1, 2.2, 2.3 of the proof), so

    q(z) := det(C_V + z·A_V) ∈ ℚ[z],   z = i·J

is an exact polynomial, obtained by interpolation from exact rational determinants. It is
even in z, and for a reason rather than by measurement: the chain is bipartite, so with the
site parity G_p = diag(∏(−1)^l over the occupied sites) the involution P = G_p ⊗ G_q fixes C
and flips A₀ entry for entry, whence det(C + z·A₀) = det(C − z·A₀) on every block, side lines
included. So q(i·J) is real, and its real roots are the exceptional couplings, counted and
isolated by Sturm sequences. Gate blocks E4 and E5, the second of which also checks that at
every exact root the FULL block, which the subspace never assumed anything about, is singular
beyond the floor and so carries AT LEAST ⌊N/2⌋ + 1. That the count is exactly one more is the
next section's, by simplicity of the roots.

The counts in the table are counts on V, so they are the couplings of RUNG 2 and not every
exceptional coupling; the rungs above contribute their own, and are counted two sections down.

| N | dim block | dim V | deg q | rung-2 J > 0 | smallest |
|---|-----------|-------|-------|--------------|----------|
| 5 | 100 | 35 | 30 | 3 | J² = 3/2, i.e. J = √6/2 |
| 6 | 225 | 84 | 78 | 6 | 0.749042443688 |
| 7 | 441 | 168 | 160 | 11 | 0.952056678 |
| 8 | 784 | 300 | 288 | 15 | 0.790268421 |

At N = 5 the smallest is the quadratic J² = 3/2 and the other two are the roots of
121J⁴ − 640J² + 832, that is J² = (320 ± 24√3)/121; the √3 there is the LARGEST mode energy in
units of J at M = 6, the one at a = 1, since ε_a = 2J·cos(πa/M). That the couplings should live in the field of
the mode energies is a reading N = 5 invites and N = 6 refuses, which
[beyond the exact reach](#beyond-the-exact-reach-and-the-laws-that-died) settles. At N = 6 the six roots sit on two irreducible
factors of degree 12 each, four on one and two on the other, and the smallest of the six,
J₀ = 0.74904244368842257949, has minimal polynomial

    64827J¹² − 1185408J¹⁰ − 912576J⁸ + 50083328J⁶ + 35151872J⁴ − 218365952J² + 102760448,

degree 12, which is why no small closed form exists and why a PSLQ search finds nothing.

**Why the extra mode is independent of the floor's ⌊N/2⌋ seeds**, with no spin argument
needed: a corner seed v is traceless, so Ψv = 0 and ΨΦv = [Ψ, Φ]v, which is a nonzero
multiple of v because the η-weight at the corner is 1 − N/2 ≠ 0, and the same holds further
up with the multiplet's own ladder coefficients. So no nonzero seed combination is
η-lowest-weight at rung p, while all of V is, and the two spaces meet only in zero.

## The whole block, not just the subspace

The singlet space V is where the extra mode was FOUND, and a count read there could in
principle miss exceptional couplings living outside it. It does not, and that is now exact
rather than sampled. The fixed part of every kernel is ker C ∩ ker A₀, a rational subspace of
dimension ⌊N/2⌋ carrying the floor's seeds; since M(z) = C + z·A₀ is complex symmetric, a
rational congruence putting a basis of it first makes the first ⌊N/2⌋ rows AND columns
vanish, and the determinant of what remains is again an exact polynomial, this time of the
WHOLE block with nothing assumed. Reading it gives:

- **N = 5, the entire band, at every coupling at once:** three exceptional couplings, on the
  blocks (2,2) and (3,3) and on no other, side lines and corner included. So at N = 5 the
  absence on the side lines is a statement about all J and not about the couplings sampled.
- **N = 6, the block (2,2):** six, which is the singlet count.
- **Every positive root is simple**, wherever the whole block has been enumerated, so the
  count at an exceptional coupling is EXACTLY ⌊N/2⌋ + 1 there: a kernel of dimension
  ⌊N/2⌋ + k would force the determinant to have a root of multiplicity k. That is the whole
  band at N = 5 and the block (2,2) at N = 6. At N = 7 simplicity is not tested anywhere, so
  only "at least one more" holds.

The two lists then coincide with no numerical comparison anywhere. At N = 5 the gate checks
that by DIVISIBILITY of the exact polynomials rather than by comparing digits; at N = 6 the
two exact counts are computed separately and agree, which with the same subset argument is
what closes it there. A root of the
singlet
polynomial gives a nonzero v in V with M(J₀)v = 0, and V meets ker C ∩ ker A₀ only in zero
(the seed argument above), so J₀ is a root of the whole block's polynomial too. The singlet
roots are therefore a subset, and equal counts make the two sets equal.

## The rung the block cannot see

Every count in the sections above was read on the block (2,2), and at N = 6 and beyond that
block cannot answer the question it was asked. A rung-ℓ lowest weight and its η multiplet occupy the diagonal blocks ℓ through
N − ℓ, so a block (p, p) carries the rungs ℓ ≤ min(p, N − p) and no others. At N = 6 the
block (2,2) therefore sees rung 2 alone, by construction, and so does its chiral partner
(4,4), for which min(4, 2) = 2 as well.

The block (3,3) sees rung 3 too, and it carries **ten** couplings: the six of rung 2, digit
for digit, and four more. The four are exact. The η-lowest-weight singlets of (3,3) form a
space of dimension 70, both halves of L restrict to it, and Sturm on the exact determinant
polynomial counts exactly four positive roots,

    0.804151383,  0.969363535,  1.405807973,  2.473408936,

none of which is one of the six: the gcd of the two exact polynomials carries no positive
root at all. At each of the four the nullity rises on (3,3) and on no other block of the
band, which is what a rung-3 lowest weight at N = 6 must do, its η-spin being N/2 − 3 = 0 so
that its multiplet is a single block, and which a rung-2 mode cannot do, those occupying
2 through N − 2. Gate block E9.

So higher rungs do contribute, and the numbers the sections above report are the **rung-2
count** rather than the count. At N = 7 the whole band was enumerated numerically, every
block of it: the chiral partners (2,2) and (5,5) carry the same eleven, the partners (3,3)
and (4,4) carry the same twenty-four, containing those eleven, and every side line carries
none at any coupling. The complete count per N, which lives in the middle block, then reads

| N | rung 2 | rung 3 | complete |
|---|--------|--------|----------|
| 5 | 3 | (no rung 3) | 3 |
| 6 | 6 | 4 | 10 |
| 7 | 11 | 13 | 24 |
| 8 | 15 | 25 | not read |

Every cell of the two rung columns is an exact Sturm count on an η-lowest-weight singlet
space, of dimension 70 for rung 3 at N = 6 and 294 for rung 3 at N = 7, except the rung-3
twenty-five at N = 8, which is measured; the rung-2 fifteen beside it is exact.

The complete column is a different kind of number in its first row than in the two below it.
At N = 5 it is exact and it is everything: the whole band was enumerated over ℚ at every
coupling. At N = 6 and N = 7 it is the whole-block count of the middle block read numerically,
and what the exact rung counts add up to, 6 + 4 and 11 + 13, meets it. That the sum is
everything, that no exceptional mode outside the singlet spaces hides in those blocks, is what
the agreement supports and does not prove. Every rung is visible in the middle block, since a
rung-ℓ multiplet spans the blocks ℓ through N − ℓ; at N = 8 that block is (4,4) of dimension
4900 and rung 4 lives nowhere else, so the complete count there is not read at all, and the
twenty-five is only what the block (3,3) adds over the rung-2 fifteen.

## Beyond the exact reach, and the laws that died

The exact route runs out after N = 8, where dim V reaches 300 and the determinant is
interpolated from 301 exact rational determinants, one per node; the dimensions grow 35, 84,
168, 300 and the cost with their fourth power. The pencil is still readable well past that, so the rung-2
count continues, by two detectors that share only the two builders: the finite spectrum of the
deflated pencil, and a sign bracket on the determinant, which is a real even polynomial in J
so that a simple root is a sign change. At N = 8, where the exact count is also available,
all three agree on fifteen, the third being a count of the determinant's sign changes over the
whole range on a fine grid, which uses no pencil at all.

| N | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
|---|---|---|---|---|---|----|----|
| rung-2 count | 3 | 6 | 11 | 15 | 23 | 29 | 38 |

exact through N = 8, measured from N = 9. Exact means the count on the singlet space, which
is the rung-2 count itself only where the whole block has been enumerated, at N = 5 and N = 6;
from N = 7 on the identification rests on the side lines staying empty and is not proved. No law is visible and several natural ones die:
(N − 4)² + 2 and ⌊((N − 2)² − 3)/2⌋ both fail first at N = 8, C(N − 2, 2) fails at N = 7, and
the quadratic fitted through the three odd values predicts 39 at N = 11 where the measurement
gives 38. Three points fixing a three-parameter family is a prediction and never evidence,
which is the only use those fits were put to.

One reading of the couplings themselves also dies here. At N = 5 all three squared couplings
lie in ℚ(cos(π/6)) = ℚ(√3), the field of the mode energies ε_a = 2J·cos(πa/M), and the √3 in
the section above invites reading that as structure. At N = 6 not one of the six lies in the
mode field: each squared coupling is algebraic of degree 6 over ℚ while ℚ(cos(π/7)) has
degree 3, and 6 does not divide 3. The two irreducible factors do sit differently beside it, which is worth having.
The one carrying four roots splits over the mode field into quadratics, so each of those four
generates a field CONTAINING it; the one carrying two stays irreducible of degree 6 there, so
those two generate a field meeting it in ℚ alone. Either way the containment runs the wrong
way for the reading N = 5 invited. N = 5 is a coincidence of the smallest case, which is the N
this arc has had to hold apart before.

**A trap for whoever reads the pencil next.** Near J = 0 both detectors report nonsense, and
not because of tolerances: C is diagonal and vanishes on every pair whose disagreement is 2,
so z = 0 is a root of the determinant with high multiplicity, and its numerical image smears
over a radius of about 10⁻⁴, which is the fourth root of the machine epsilon and so the size a
root with Jordan blocks of size four smears to. The candidate
at 2.4·10⁻⁴ at N = 8 is that smear and not a coupling; counting it turns the fifteen into a
sixteen and restores a law that is not there. The tempting exact test, comparing the nullity
of A₀ compressed to ker C against the floor, is not a test: that surplus is large at every N
including N = 5 and N = 6, where the exact count says there is no small root. Singularity
there raises the multiplicity of the root AT zero and says nothing about a root beside it.

## Measured, and named as such

- The rung-2 count, 3, 6, 11, 15 at N = 5 through 8, is exact per N but has no law. From
  N = 7 on it is the count on the SINGLET space; the whole block has not been enumerated
  there over ℚ.
- At **N = 7** only the singlet list is exact. The whole band there has been enumerated
  numerically, block by block and at every coupling at once, which is what says the side lines
  stay empty and what gives the twenty-four of the middle blocks; over ℚ none of that is run.
- The continuation of the rung-2 count past N = 8, and the rung-3 count at N = 8, are numeric
  throughout. Two detectors agree on every value, and at N = 8, where the rung-2 count is also
  exact, all three agree there. That is evidence and not a proof. The rung-3 count at N = 7 is
  not in this list: it is exact, by the same Sturm route, under `--slow`.
- WHICH blocks are raised at N = 6 is read numerically, at one coupling. The corresponding
  N = 5 statement is exact and over all J. The N = 6 COUNT is exact under `--deep`, where the
  whole block (2,2) is enumerated over ℚ; only the block-by-block reading there is numerical.
- Whether the exceptional points are semisimple or defective is not read here. The
  neighbouring question is answered in the other direction by
  [PROOF_R90_FROZEN_DIVISOR](../docs/proofs/PROOF_R90_FROZEN_DIVISOR.md) §9, and see the word
  warning below.

## The word that is already spent

**Exceptional coupling** is already a repo term, and it means something adjacent but not the
same. In `PROOF_R90_FROZEN_DIVISOR` §9 it is the corner block on the R90 locus with a
non-uniform profile, where the ALGEBRAIC multiplicity gains a rung and the point turns out
to be defective while the geometric count stays at ⌊N/2⌋. Here it is the uniform point, the
rungs ℓ ≥ 2, and the jump is GEOMETRIC. Both are roots of the same kind of tightness
criterion, so the word fits; what must never be dropped is which block and which
multiplicity. The same document also proves that at the uniform point the corner has no
exceptional coupling at all, which is exactly what the scans show: the corner block (1,1)
carries ⌊N/2⌋ at every coupling tested, and the extra mode appears first at rung 2.

## What is open

1. **A law for either count.** The rung-2 count reads 3, 6, 11, 15, 23, 29, 38 at N = 5
   through 11 and the complete count 3, 10, 24 at N = 5, 6, 7. Nothing is known about how
   either grows, and the candidates tried are dead ones, listed above.
2. **How far up do the rungs go?** Rung 3 carries couplings of its own, four at N = 6 and
   thirteen at N = 7. No rung above 3 has been read at any N, and nothing yet says whether
   every rung ℓ ≤ ⌊N/2⌋ carries some or whether the contribution stops.
3. **Semisimple or defective** at the exceptional point, and whether the answer differs from
   §9's, where the corner's exceptional points are defective.
4. **The off-diagonal side lines.** At N = 5 they carry no exceptional coupling at ANY J,
   exactly. At N = 7 the same is now read numerically over the whole band. If the extra mode
   is always a singlet the absence is forced, but the forcing has not been written down.
5. **What it means.** A coupling at which the chain freezes one more mode than it generically
   can is a resonance between the turning and the watching, and the note has not asked what
   the resonance is between.

## How it was found, and the two traps on the way

The first pass scanned the coupling on a grid and tracked the (⌊N/2⌋+1)-th singular value. It
found one dip, and only because a grid point happened to land within 1.7·10⁻⁴ of a root. **A
grid scan under-detects here by construction**, since the crossing is transversal and a root
is a point. The honest tool is the pencil's finite eigenvalues, which give every candidate at
once, then bisection to confirm each, and only then the exact polynomial. Read the other way
round, this is why the exceptional couplings were not noticed earlier: nothing in the arc ever
scanned J, because every statement in it was designed to be coupling-free.

The second trap was the identification. The first exceptional coupling found, at N = 6, was
pinned to more than thirty digits and handed to an integer-relation search, which found
nothing. That was
correct and it was not a failure of precision: the number is algebraic of degree 12. The
lesson is that "no closed form appears" is a finding when the exact minimal polynomial is in
hand, and an absence of evidence when only digits are.

The gate builds every object from scratch, and its first block checks the block generator
against the same block cut out of the dense 4^N Liouvillian, entry for entry, so that the
agreement with [XY_FROZEN_BAND](XY_FROZEN_BAND.md) and
[ETA_CEILING_REDUCTION](ETA_CEILING_REDUCTION.md) is evidence rather than a shared bug.
