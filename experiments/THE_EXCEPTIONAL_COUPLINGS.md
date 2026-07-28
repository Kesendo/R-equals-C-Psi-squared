# Where the ceiling actually fails: the couplings at which one more mode freezes

**Status:** the measuring, deliberately unnamed. The existence and the counts are exact; no
F number, because the count 3, 6, 11 has no law yet and because from N = 6 on no small closed
form appears for the couplings themselves, so the registry of closed forms is the wrong home
for either.
(Not proved: that no radical expression exists. What is shown is that the minimal polynomial
has degree 12 and that no small relation appears, which is a weaker and sufficient reason.)
What is proved and what is measured are separated below and should stay that way.
**Date:** 2026-07-28
**Verification:** [`simulations/exceptional_couplings.py`](../simulations/exceptional_couplings.py)
(must print "exceptional couplings gate: ALL GREEN", 35 checks, about a minute; `--deep` adds
N = 7 and the whole-block enumeration at N = 6, 41 checks in about nine minutes)
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
which is a coupling-independent count and has nothing to do with this.

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
exceptional coupling the count is one more, and exactly one more wherever the whole block has
been enumerated, which is N = 5 and N = 6.

Two things make it worth a note rather than a footnote:

1. The extra mode is identified. It is η-lowest-weight at rung ℓ = 2, so it occupies the
   diagonal rungs 2 through N − 2 and nothing else, and it is a **spin singlet**, where the
   band's own modes carry j_spin = 1. So it enlarges the diagonal line and never steps
   sideways: the band statement |p − q| ∈ {0, 2} is untouched, and so is the count on the two
   side lines. That last sentence is exact over all couplings at N = 5, read at one coupling
   at N = 6, and unmeasured at N = 7; the sections below keep those apart.
2. The bridge that was open becomes elementary once the failure is known. Nullity of a
   matrix pencil is constant off a finite set and minimal there, and the floor holds
   everywhere, so a **single** coupling with nullity ≤ ⌊N/2⌋ settles the generic value.
   Proposition 5.1 with F144 supplies one for every N ≥ 6 in the limit, and an exact GF(p)
   rank supplies a concrete finite one, J = 1, at every N the rank reaches. No semicontinuity
   argument is needed, and none would have delivered more, because "all but finitely many" is
   the true statement and not a technical hedge.

## Proved, and how

The existence, the counts, and everything stated for N = 5 are exact arithmetic over ℚ, with
no floating point anywhere; at N = 6 the whole-block count is exact too, under `--deep`. Two
things are numerical and are marked where they appear: the default run's confirmation at
N = 6 and 7 that a root of the singlet polynomial makes the whole block singular, which the
exact enumeration supersedes at N = 6, and the reading of WHICH blocks are raised at N = 6.

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
block, the η-lowest-weight spin singlets. (V is not V₀ = ker(ad_h), the space the large-J
reduction of [ETA_CEILING_REDUCTION](ETA_CEILING_REDUCTION.md) lives on; the two are different
objects and the neighbouring letters are worth guarding against.) All three maps have integer matrices, and both
halves of L restrict to V exactly (Lemmas 2.1, 2.2, 2.3 of the proof), so

    q(z) := det(C_V + z·A_V) ∈ ℚ[z],   z = i·J

is an exact polynomial, obtained by interpolation from exact rational determinants. It is
even in z, so q(i·J) is real, and its real roots are the exceptional couplings, counted and
isolated by Sturm sequences. Gate blocks E4 and E5, the second of which also checks that at
every exact root the FULL block, which the subspace never assumed anything about, is singular
beyond the floor and so carries AT LEAST ⌊N/2⌋ + 1. That the count is exactly one more is the
next section's, by simplicity of the roots.

| N | dim block | dim V | deg q | exceptional J > 0 | smallest |
|---|-----------|-------|-------|-------------------|----------|
| 5 | 100 | 35 | 30 | 3 | J² = 3/2, i.e. J = √6/2 |
| 6 | 225 | 84 | 78 | 6 | 0.749042443688 |
| 7 | 441 | 168 | 160 | 11 | 0.952056678 |

At N = 5 the smallest is the quadratic J² = 3/2 and the other two are the roots of
121J⁴ − 640J² + 832, that is J² = (320 ± 24√3)/121; the √3 there is the top mode energy in
units of J at M = 6, since ε_a = 2J·cos(πa/M). At N = 6 the six roots sit on two irreducible
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

The two lists then coincide with no numerical comparison anywhere, and the gate checks the
same by DIVISIBILITY of the exact polynomials rather than by comparing digits. A root of the
singlet
polynomial gives a nonzero v in V with M(J₀)v = 0, and V meets ker C ∩ ker A₀ only in zero
(the seed argument above), so J₀ is a root of the whole block's polynomial too. The singlet
roots are therefore a subset, and equal counts make the two sets equal.

## Measured, and named as such

- The count of exceptional couplings, 3, 6, 11 at N = 5, 6, 7, is exact per N but has no law.
  At N = 7 it is the count on the SINGLET space; the whole block has not been enumerated there.
- At **N = 7** only the singlet list is exact; the whole-block enumeration has not been run
  there, so completeness at N = 7, and with it "every exceptional mode is an ℓ = 2 singlet"
  as a general statement, remains measured. Nothing says a higher rung cannot contribute at
  some larger N.
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

1. **A law for the count** 3, 6, 11. Nothing is known about how it grows.
2. **Do higher rungs ever contribute?** Measured absent at N = 5, 6, 7.
3. **Semisimple or defective** at the exceptional point, and whether the answer differs from
   §9's, where the corner's exceptional points are defective.
4. **The off-diagonal side lines.** At N = 5 they carry no exceptional coupling at ANY J,
   exactly. At N = 6 that is read at one coupling only. If the extra mode is always a singlet
   the absence is forced, but the forcing has not been written down.
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
