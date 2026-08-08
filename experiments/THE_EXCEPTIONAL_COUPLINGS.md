# Where the ceiling actually fails: the couplings at which one more mode freezes

**Status:** the measuring, deliberately unnamed. The existence is proved. The count on the
singlet space of the block (2,2) is exact through N = 8 and measured beyond, and it is the
full rung-2 count provably through N = 6, where the whole block was enumerated. The complete
count, over all rungs, is read through N = 8 and exact only at N = 5, where the whole band was
enumerated at every coupling at once; every numeric BLOCK count anywhere below is a verified
LOWER bound, and a difference of two of them is bounded in neither direction, for the reason
given in
[beyond the exact reach](#beyond-the-exact-reach-and-the-laws-that-died). No F number, because
no law is known for any of the counts and because from N = 6 on no small closed form appears
for the couplings themselves, so the registry of closed forms is the wrong home for either.
(Not proved: that no radical expression exists. What is shown, at N = 6, is that the minimal
polynomial of the smallest coupling has degree 12 and that no small relation appears, which is
a weaker and sufficient reason; at N = 7 and 8 not even that is read.)
What is proved and what is measured are separated below and should stay that way.
**Date:** 2026-07-28
**Verification:** [`simulations/exceptional_couplings.py`](../simulations/exceptional_couplings.py)
(must print "exceptional couplings gate: ALL GREEN"; the run states its own check count, so no
number is carried here to drift). Under `--rungs` that criterion covers the numbers in the
table two ways, and deliberately by two separate checks: one fails if a block count came out
BELOW what is recorded here, which is the detector regressing, and one fails if it came out
ABOVE, which is not a regression at all but means this note is out of date and the new root
wants writing up. A single one-sided check would have let the second case pass green. The default run reaches N = 6; `--deep` adds N = 7, the
whole-block enumeration at N = 6 and a numeric read of the N = 7 rung; `--slow` adds the two
exact counts that cost about half an hour, the rung-3 count at N = 7 and the rung-2 count at
N = 8; `--rungs` adds the numeric block counts at N = 8 up to the middle block, where the
polynomial is out of reach and the pencil is not, and turns the other two flags on with it,
about an hour and a half in all. It needs them: it validates a numeric detector against counts
this note proves exactly, and it reads every one of those expected values out of the exact
blocks of the same run rather than carrying them as literals, so a stale one cannot survive
in the validation. What remains a session
measurement, in no mode of the gate, is the whole band at N = 7, the continuation of the rung-2
count past N = 8, the grid detector, the mode-field factorisations at N = 6 and the dead
small-J criterion; those are recorded here as observations.
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
which is a coupling-independent count and has nothing to do with this. A **pair** in the
sections on the numeric route is always the ±z pair a single root arrives as in the pencil
spectrum, never MirrorWorld's `Pair`, which is a bare coherence |i⟩⟨j|. A fourth word joins
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
   last sentence is exact over all couplings at N = 5, read at five couplings at N = 6, and
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
- **Every positive root is simple**, wherever the whole block has been enumerated EXACTLY,
  over ℚ, so the count at an exceptional coupling is EXACTLY ⌊N/2⌋ + 1 there: a kernel of
  dimension ⌊N/2⌋ + k would force the determinant to have a root of multiplicity k. That is
  the whole band at N = 5 and the block (2,2) at N = 6, and nowhere else. At N = 7 and N = 8
  simplicity is not tested anywhere, so only "at least one more" holds. The numeric block
  reads there do not extend it: a sign change proves a root of ODD order, which is not the
  same as a simple one.

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

The block (3,3) sees rung 3 too, and it carries **at least ten** couplings: the six of rung 2,
digit for digit, and four more. Both lists are exact, so the ten are certain; what is not
proved is that there is no eleventh, since the whole block (3,3) was never enumerated over ℚ
at N = 6. It was read numerically, and the read returns those ten and nothing else, value for
value: the gate requires agreement to 10⁻⁹ relative and the run reads 4·10⁻¹⁵. That is the
validation the numeric detector is put through before it is trusted at N = 8. The four are exact. The η-lowest-weight singlets of (3,3) form a
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
none at any coupling. At N = 8 the three diagonal blocks (2,2), (3,3) and (4,4), of dimensions
784, 3136 and 4900, were read the same way, and their root lists nest value for value; that is
gate block E11 under `--rungs`, and only those three blocks were read there, no side line and
no partner. The complete count per N, which lives in the middle block, then reads

| N | rung 2 | rung 3 | rung 4 | complete |
|---|--------|--------|--------|----------|
| 5 | 3 | (no rung 3) | (no rung 4) | 3 |
| 6 | 6 | 4 | (no rung 4) | 10 |
| 7 | 11 | 13 | (no rung 4) | 24 |
| 8 | 15 | 25 | 12 | 52 |

The empty cells are empty by construction and not for want of measuring: a rung-ℓ multiplet
spans the blocks ℓ through N − ℓ, so a rung needs ℓ ≤ ⌊N/2⌋ to exist at all, and rung 4 first
becomes possible at N = 8.

Every cell of the three rung columns is an exact Sturm count on an η-lowest-weight singlet
space, of dimension 70 for rung 3 at N = 6 and 294 for rung 3 at N = 7, except the rung-3
twenty-five and the rung-4 twelve at N = 8, which are measured; the rung-2 fifteen beside them
is exact. The twenty-five has been reached twice, and only one of the two is in the gate: the
session measurement on the singlet space of (3,3) at N = 8, and the difference of the numeric
(3,3) and (2,2) counts, which is what the gate computes. Their agreement checks the singlet
identification at that block, since the two spaces are different; it is not a second detection,
the pencil being the same in both. Inside the gate alone the twenty-five is the subtraction and
nothing else, which is why the gate prints it rather than asserting it as a reproduction.

**What the two differences at N = 8 do and do not license.** The three block counts nest, so
40 − 15 counts what (3,3) adds over (2,2) and 52 − 40 counts what (4,4) adds over (3,3). Each
accepted root is certified on its own by a sign change of the determinant, and nothing
certifies that no root was MISSED, so each block count is a lower bound on its block. **A
difference of two lower bounds is bounded in neither direction**, and that is the whole of the
caution here: if the (2,2) list is short the true rung 3 is smaller than 25, and if the (3,3)
list is short it is larger. What the measurement supports with nothing further assumed is only
this: at least 25 certified roots of (3,3) lie outside (2,2)'s list, and at least 12 of (4,4)
lie outside (3,3)'s. Add the assumption that the SMALLER list is complete and the difference
becomes a lower bound on the rung. Equality, that is calling 25 and 12 the rung counts, needs
two further things: both lists complete, and the two rungs sharing no coupling. The second is
easy to overlook and is not implied by the first. A coupling exceptional for a rung-2 mode AND
a rung-3 one would sit in (2,2)'s list, be subtracted away, and leave the true rung-3 count
above 25. The note checks exactly that at N = 6 and at N = 7, by taking the gcd of the two
exact polynomials and finding no positive root, N = 7 being the last N where both polynomials
exist; at N = 8 there is no exact polynomial to take a gcd of, and nothing checks it.

**The twelve is the one number in the table that has never been read a second way at all.** The
twenty-five has the session measurement beside it, sharing the detector; the twelve has nothing
beside it. The subtraction that produces it is itself validated wherever an exact rung count
exists to validate it against:
10 − 6 = 4 at N = 6 and 24 − 11 = 13 at N = 7 both reproduce the exact Sturm count on the
rung-3 singlet space, and the gate runs both. Taken alone each tests completeness and
disjointness together, so a failure could not say which; the gcd above separates them at both
N, leaving completeness as what the agreement actually tests. No such check exists for rung 4 at N = 8 and none can: the rung-4 multiplet there
has η-spin 0, so it lives in (4,4) alone and there is no larger block to read it against. The
one route that would give the twelve a second reading is a Sturm count on the rung-4 singlet
space of (4,4), the same machinery one rung up; it has not been run, and its dimension is what
would decide whether it can be.

The complete column is a different kind of number in each of its four rows, and the gradient
runs one way. At N = 5 it is exact and it is everything: the whole band was enumerated over ℚ
at every coupling. At N = 6 and N = 7 it is the middle block read numerically, and what the two
EXACT rung counts add up to, 6 + 4 and 11 + 13, meets it independently. At N = 8 there is no
such meeting at all, because only the rung-2 fifteen is exact there: 15 + 25 = 40 is the (3,3)
count, and nothing outside the (4,4) read itself says anything about the 52. That the sum is
everything, that no exceptional mode outside the singlet spaces hides in those blocks, is what
the agreement supports and does not prove.

The word **complete** in that column means the rung inventory is closed, and at N = 8 it is.
Every rung is visible in the middle block, since min(⌊N/2⌋, N − ⌊N/2⌋) = ⌊N/2⌋. The rungs that
could carry anything run ℓ = 2 through ⌊N/2⌋ = 4: rung 0 is the vacuum cell, where C acts as
+4γ̄ and A₀ as zero, so it is never frozen at any coupling, and rung 1 is the corner, which is
excluded by a theorem rather than by a scan, `PROOF_R90_FROZEN_DIVISOR` §7 making the
determinant a J-monomial at the uniform point so that the multiplicity is exactly ⌊N/2⌋ for
every J ≠ 0 (see [the word that is already
spent](#the-word-that-is-already-spent)). So 15 + 25 + 12 is all of it. Closed as an inventory
of rungs, which is a statement about which rungs can exist and not about how well any of them
was counted: the fifteen is exact on the singlet space, and the twenty-five and the twelve
carry the caution of the paragraph above.

## Beyond the exact reach, and the laws that died

The exact route runs out after N = 8, where dim V reaches 300 and the determinant is
interpolated from 301 exact rational determinants, one per node; the dimensions grow 35, 84,
168, 300 and the cost with their fourth power. The pencil is still readable well past that, so the rung-2
count continues. The route has two parts and they are not the same kind of thing. **Detection**
is the finite spectrum of the deflated pencil, which offers every candidate at once.
**Certification** is a sign bracket on the determinant, which is a real even polynomial in J so
that a root of odd order is a sign change; it is applied only to candidates the pencil already
proposed, so it can confirm one and can never find one. The independent second detector is a
different tool again: a count of the determinant's sign changes over the whole range on a fine
grid, which uses no pencil at all. At N = 8, where the exact count is also available, the
pencil, the grid and the exact polynomial agree on fifteen. The grid is not a safety net for
what the pencil misses, since a grid under-detects in its own way, two roots inside one cell
showing as no sign change at all.

**Every numeric BLOCK count in this note is a verified lower bound, and the asymmetry is
structural.** The word block is doing work: the rung counts 25 and 12 are differences of two
such bounds and are not bounds themselves, as the section on the table sets out. The sign
bracket certifies each ACCEPTED root, since a sign change across a
window narrower than the gap to the next root proves an odd-order root inside it; nothing
certifies that no root was MISSED, because detection is by a pencil spectrum and a genuine
root can be pushed off the real axis by a nearby cluster. So the numbers below, the numeric
block counts at N = 7 and N = 8, and the complete counts that rest on them, are bounds from
one side only. This matters most where a law is being fitted: a curve through counts that may
be truncated can fit perfectly because the truncation delivered the number the curve wanted.

Being a bound even in that one direction is not free, and exactly one premise stands under it.
A count is too high if one root is certified twice, and two windows that both flip can only be
one root twice if they both contain it; if they both contain ρ then J₂(1 − ε) ≤ ρ ≤ J₁(1 + ε),
so their gap is at most ε(J₁ + J₂). **Disjoint windows therefore contain distinct roots and the
count is a lower bound**, and the gate measures the disjointness rather than assuming it. It is
sufficient and not equivalent, which is worth saying since the tempting stronger claim is
false: two overlapping windows each holding a root of its own outside the overlap also count
correctly, so an overlap would mean the premise is no longer measured, not that the bound has
failed. Sufficient is all the count needs. The margin is the tight number of the whole
construction: at N = 8 the smallest gap between consecutive accepted roots is 2.8 times the
distance ε·J the two windows reach toward each other, which is **1.4 times their combined
width**. A factor, not decades, and the one to watch at N = 9, where ε is a fixed constant and
the roots are getting denser.

Two things sit beside that premise without being it. The **±z parity**: the pencil carries
every root as a pair and the two copies are merged before counting. A pair that fails to merge
does happen, and at exactly the block that matters, the z = 0 smear arriving at (4,4) as two
unmerged candidates, which the gate's record now shows as merged 1 where all 52 accepted roots
show merged 2. But a split pair that survived would be caught by the disjointness above, and a
LOST copy is under-detection, which a lower bound tolerates. So parity is detector health, not
a second premise, and the note does not lean on it. And the **small-J cut**, which is read
rather than trusted: the threshold 0.05 says nothing on its own, but the smallest accepted root
sits 2065 times above the largest candidate rejected as smear, so nothing there was decided by
the threshold. Those two split copies would in any case not have been certified: their windows
sit at 2.6·10⁻⁴ ± 1e-5 relative and do not contain 0, so the determinant does not change sign
across either.

| N | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
|---|---|---|---|---|---|----|----|
| rung-2 count | 3 | 6 | 11 | 15 | 23 | 29 | 38 |

exact through N = 8, measured from N = 9. Exact means the count on the singlet space, which
is the rung-2 count itself only where the whole block has been enumerated, at N = 5 and N = 6;
from N = 7 on it rests on the numeric whole-block read returning exactly the singlet count, at
N = 7 and at N = 8 alike, and is not proved. (The empty side lines are not what supports it:
the forced implication runs from singlet-ness to the emptiness, not back.) No law is visible
and several natural ones die:
(N − 4)² + 2 and ⌊((N − 2)² − 3)/2⌋ both fail first at N = 8, C(N − 2, 2) fails at N = 7, and
the quadratic fitted through the three odd values predicts 39 at N = 11 where the measurement
gives 38. Three points fixing a three-parameter family is a prediction and never evidence,
which is the only use those fits were put to.

The complete count, now that it reads 3, 10, 24, 52 at N = 5 through 8, has one live
candidate: those four satisfy a(N + 1) = 2·a(N) + 4, that is a(N) = 7·2^(N−5) − 4, exactly.
It is a candidate and stays one, for three reasons stacked the same way. Four points against a
two-parameter family is two degrees of confirmation. This arc has already buried three laws
that fitted every point they were fitted to, listed just above. And three of the four points
are lower bounds, so the exactness of the fit is a reason for suspicion rather than for
comfort. The value it predicts at N = 9 is 108, on the block (4,4) of dimension 15876, which
cubic scaling puts near twenty hours and which has not been run.

One reading of the couplings themselves also dies here. At N = 5 all three squared couplings
lie in ℚ(cos(π/6)) = ℚ(√3), the field of the mode energies ε_a = 2J·cos(πa/M), and the √3 in
the section above invites reading that as structure. At N = 6 not one of the six lies in the
mode field: each squared coupling is algebraic of degree 6 over ℚ while ℚ(cos(π/7)) has
degree 3, and 6 does not divide 3. The two irreducible factors do sit differently beside it, which is worth having.
The one carrying four roots splits over the mode field into quadratics in J², so each of those
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
sixteen and restores a law that is not there. The gate rejects it by name rather than silently
and prints each rejection of the three blocks it measures: at N = 8 that is one candidate on
(2,2), one on (3,3) and two on (4,4), all of them between 2.4 and 2.9·10⁻⁴. Every rejection in
the run is one of these; at the smaller N the detector produces no candidate down there at all.
The tempting exact test, comparing the nullity
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
  The twenty-four is now inside the gate as well, under `--rungs`, where it also serves as the
  detector's validation; the side-line statement over the whole band there remains a session
  measurement.
- The continuation of the rung-2 count past N = 8 is numeric throughout. The pencil and the
  grid, which detect independently, agree on every value, and at N = 8, where the rung-2 count
  is also exact, all three agree there. That is evidence and not a proof. The rung-3 count at N = 7 is not in this list: it is exact,
  by the same Sturm route, under `--slow`.
- At **N = 8** the block counts 40 and 52, and with them the rung-3 twenty-five, the rung-4
  twelve and the complete fifty-two, are numeric. They are gated under `--rungs` rather than
  left to a session, and inside that gate the detector is first made to reproduce counts this
  note proves exactly elsewhere, including the fifteen at the very N it then measures, though
  on the block BELOW the two being measured, 784 against 3136 and 4900, and conditioning is the
  stated risk. That row also compares a whole-block numeric count against a singlet exact one,
  so it carries the open identification with it. Each
  accepted root is certified and each BLOCK count is a lower bound; the two differences are
  not bounded in either direction, for the reason given with the table.
- WHICH blocks are raised at N = 6 is read numerically at five couplings, and over the whole
  band at each: one rung-2 coupling in gate block E6 and all four rung-3 ones in E9(c). At the
  other six E9 reads one block only, whether (3,3) is raised, which is a smaller statement.
  The corresponding
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
   through 11 and the complete count 3, 10, 24, 52 at N = 5 through 8. For the rung-2 count
   nothing is known and every candidate tried is dead, listed above. For the complete count
   there is one live candidate, 7·2^(N−5) − 4, which is fitted to four points of which three
   are lower bounds and is therefore not evidence of anything yet; the section above says why,
   and the number that would test it is 108 at N = 9.
2. **How far up do the rungs go?** Rung 3 carries couplings of its own, four at N = 6 and
   thirteen at N = 7, and rung 4 carries twelve at N = 8, the first N at which it can exist.
   Rung 5 is unread and needs N = 10. So no rung yet reached is empty, and nothing says
   whether every rung ℓ ≤ ⌊N/2⌋ carries some or whether the contribution stops.
3. **Semisimple or defective** at the exceptional point, and whether the answer differs from
   §9's, where the corner's exceptional points are defective.
4. **The off-diagonal side lines.** At N = 5 they carry no exceptional coupling at ANY J,
   exactly. At N = 7 the same is now read numerically over the whole band. At N = 8 they were
   not read at all, which is why the complete count there rests on three diagonal blocks where
   the N = 5 one rests on the whole band. If the extra mode is always a singlet the absence is
   forced, but the forcing has not been written down. There is a proposed argument that would
   close the N = 8 gap without measuring anything, and it is recorded here as proposed and not
   as established: S± commutes with L and shifts p − q by 2, so every frozen multiplet carries
   INTEGER spin and therefore has a nonzero zero-weight component on a diagonal block, at a
   rung the middle block sees. If that holds, no exceptional coupling of the band can be
   invisible to (4,4), and completeness never needed the side lines to be empty, only that they
   add nothing new. Worth checking before it is used; the same shape of gift has been wrong in
   this arc before.
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
