# What the R₉₀ locus buys, and why two results are standing on it

*2026-08-14, repaired 2026-08-15. [F153](../docs/ANALYTICAL_FORMULAS.md) puts a pinned block
(min(p,q) = 0 or max(p,q) = N) entirely on its decay floor at UNIFORM γ, its first fence; under
a profile a block of dimension above one leaves the floor, and what an involution then pairs its
cells about is the trace constant, a different number.
[F140](../docs/ANALYTICAL_FORMULAS.md) freezes an eigenvalue in the corner block by counting
the rooms an involution leaves short. Both stand on the same anti-palindromic γ locus, the page
says they are different objects, and nothing said why one locus carries both. This note answers
that. Almost every ingredient is already the repo's, and six review rounds were needed to find
that out; what the note contributes is one identification, of F153's size-class condition with
the defect F140 carries on its diagonal cells, and it is an identification of the defect MATRIX
and not of a constant. The observation that the two vanish on the zero-mean stratum is NOT ours,
it is PROOF_R90_FROZEN_DIVISOR §5's; what is ours there is the reach of the measurement. Most of
what it says about F153 is scope rather than error, with TWO exceptions it also repaired in
`docs/ANALYTICAL_FORMULAS.md`: F153's entry gave the wrong REASON for the longitudinal field
leaving the criterion alone, a copy of which this note had inherited, and its Δ = 0 / Δ = 0.5
saturation reading was missing the clause saying it is OFF the locus, which is the very omission
the surrounding paragraphs were written to repair.*

## What this is about

Two of the project's registered laws stand on the same special lighting of
the chain, the mirror-balanced profile in which each site and its mirror
partner together receive the same total light
([Dephasing Translated](../docs/quantum/DEPHASING_TRANSLATED.md)). One law
concerns the families of patterns in which one of the two connected spin
arrangements is completely empty or completely full, and says where an uneven
lighting lets their decay rates sit; the other freezes eigenvalues in a different family, the one
whose patterns connect two arrangements holding a single excitation each, at
a rate that no strength of coupling can move. The registry treated them as
different objects, and nothing said why one lighting carries both. This note
answers that: they are not two tenants of one address but two readings of one
and the same act of reflection, taken by two genuinely different instruments. For a chain whose bonds also read the
same from either end, reflecting it pairs the decay rates of a family about a
centre, and the centre is set by how many sites the two arrangements a
pattern connects disagree at. Where a family holds patterns of two different
such counts, no single centre fits both; the mismatch vanishes only where the
lighting averages out to zero, which for ordinary light means no light at
all. The mismatch that remains is a
concrete leftover, a matrix, and the note's one contribution is the
identification: this leftover is the very defect standing behind the freezing
law, the same object matrix for matrix, not merely the same size. Along the way the note repairs two errors
in the registry entry it examined, and withdraws half of a mechanism it had
drafted, keeping the half that a later note went on to explain.

## What the repo already held, store by store

The sweep ran before the derivation and twice more under review, and the misses are the useful
part. The first pass searched for "an antilinear involution of a block" and missed the stores
that own it, because the repo files that object as the **cross-fold** and the **self-folded
block**. The third pass found that the one sentence the note still claimed as new is in F153's
own registered text. The label was the whole difficulty each time.

`docs/ANALYTICAL_FORMULAS.md` returned **F89d** (Tier 1 derived), the identity this note uses,
typed as `F89CrossFoldSimilarityClaim` with the live witness `inspect --root crossfold`:

    L₍1,N−2₎(q̄, Δ)  =  −P · conj(L₍1,2₎(q, Δ)) · Pᵀ  −  2N·I

(that display is the registry's, at unit γ where σ = Σ_l γ_l = N, and its q is the COUPLING J/γ while the
subscripts are popcounts; everywhere else in this note q is a popcount, with one flagged
exception in the §7(c) paragraph below, which borrows that proof's letters) with P the
bra-complement permutation, machine zero for N = 4..9 at every coupling (real and complex)
and every Δ, because the Δ·ZZ frequency is bit-flip-even. Its own carrier derives the shift
from the σ form, R·L_diss·R = −L_diss − 2σ (σ = Σ_l γ_l throughout, γ̄ = σ/N its mean), docking
onto `F1PalindromeIdentity`, which holds for Z-dephasing "uniform or site-dependent".

It returned **F153**, and F153 states source 1's map and its profile generality outright:
"the bitwise COMPLEMENT of the free configuration stays inside the block and sends
rate ↦ Σγ − rate = 2·mean(rate), for ANY profile". It is repeated in the OpenArcs entry and in
`PinnedBlockFloorClaim`'s doc-comment, and it is already gated on a non-uniform off-locus
profile by `PinnedBlockFloorClaimTests` ("the complement source is blind to the bond profile").
F153 also carries the general trace constant, −2γ̄·(p + q − 2pq/N), from which the palindrome
axis of any index-N/2 block follows in one line as −Nγ̄ = −σ, for every p and any profile.
It returned **F140**, whose mechanism is a room count of the cell mirror τQ and which says
outright that no symmetry is behind it.

`docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md` returned §7, the fold lattice. Consequence **(b)**:
"a block fixed by a fold leg (q̃ = N/2 or p = N/2, even N) has an s-symmetric spectrum",
s(λ) = −λ − 2N, with no pinning hypothesis, gated by `BlockLatticeFoldGroupTests` including the
N = 6 self-folded line. Consequence **(c)** is the Δ-scoping: at Δ ≠ 0 its ingredient (iv), the
bipartite gauge, fails off the two fixed lines, and (c)'s own conclusion is that the fold legs then
connect the COUPLING q with −q instead, still holomorphically, with q-evenness surviving exactly on
the fixed lines (this paragraph borrows §7's letters, where q is the coupling and q̃ the popcount;
it is the one place this note's popcount rule for q does not apply). What this note's gate measures,
that at a FIXED coupling the linear-with-gauge form fails at Δ ≠ 0, is the fixed-coupling shadow of
that statement, not (c)'s wording; the antilinear form's survival at every
Δ is the gate's own measurement. §6 owns the rate window the size classes live in.

`docs/proofs/PROOF_R90_FROZEN_DIVISOR.md` returned τQ as `(a,b) ↦ (R(b), R(a))`, the
fixed-count argument, the statement that the (1,1)-type corner block's spectrum ("the corner"
below) is not palindromic about the root,
the defect `τQ(L + 4γ̄)τQ = −(L + 4γ̄) + 8γ̄·P_D` localised on the diagonal cells, and in §5 the
size-class criterion already applied to blocks: "the center depends on the size |S|. The corner
block has a single off-diagonal size class (|S| = 2), so one center serves; the block (2,2) has
classes |S| ∈ {0, 2, 4} and admits no single center." §5 also carries the γ̄ = 0 fence that §3
and §6 name as a hypothesis and that this note needed and at first dropped. And the proof's own
gate ledger already gates one half of the linear/antilinear split this note measures: "the mirror
identity with defect 8γ̄·P_D at machine zero for N = 4, 5, 6, and that the antilinear variant
fails (the mirror is linear)". The τQ half is committed there; what is new below is the
reflection's half.

`docs/proofs/PROOF_F91_GAMMA_NINETY_DEGREES.md` returned R₉₀ as an involution on the γ VECTOR
whose fixed-point set is the locus, with the warning that on that set it acts as the identity,
so "R₉₀ is not what carries the invariance"; what carries it is within-pair redistribution at
fixed indexed pair-sums.

`experiments/` returned [XY_FROZEN_BAND](XY_FROZEN_BAND.md) (2026-07-25), which had already
sighted the antilinear involution on the pinned block (0,2): "The involution that makes M odd
there is antilinear: reversal of the bra index composed with complex conjugation ... The
oddness is exact, entry for entry", with M defined in that file as L_block + 4γ̄, so the
constant is there too. On (0,2) the two-leg reflection collapses to bra-only reversal, so the
sighting cannot distinguish the one-leg and two-leg maps, and it states no spectral palindrome.
It also holds the argument that keeps the two results apart, which this note does not overturn:
"an antilinear involution has equal-dimensional real forms, so it yields no shortage and cannot
be the counter". That paragraph closed with "Hand it over as a fact, not as a lead."

The OpenArcs registry returned `f_registry_meets_the_typed_layer` step (3), which is this
material: both sources, the palindrome-not-flatness reading, and the disowning of F140 ("The
shared thing is the locus, and the locus is F91's"). It also returned, found only under review,
the retired arc `two_coast_classifier_repair`, which had already litigated the longitudinal-field
correction this note repairs in F153: "'a longitudinal Z field commutes with the Z-dephasing and
so cannot touch the decay rates at all' is false as stated", with
`docs/proofs/MIRROR_SYMMETRY_PROOF.md` (lines 229-235) having it right, "uniformity doing the
work, not U(1)". The repo owned the correction's SHAPE before this note re-derived it; what the
note adds is only the scalar-dissipator scoping on the blocks. It returned `f91_scope_fences`, and the
credit there took two goes to get right. An earlier draft named the arc as the SOURCE of the
identity-on-the-locus statement; the arc's six listed findings are all repairs to F91's own
prose and none of them is that statement, whose source is PROOF_F91:68, cited above where it is
used. But the arc does carry the fact and does lean on it: its NextStep kills a proposed wording
as "VACUOUS on the locus (R₉₀ is the identity there, finding (2) of this arc's own list)". So the
arc is a second reader of PROOF_F91:68 rather than its origin, and the repair that removed the
credit outright went one step too far. The typed layer returned the
two carriers, whose distinctive parents differ, F153 hanging off the Absorption Theorem and
F140 off F91, but which meet at `JointPopcountSectors`. `compute/MirrorWorld/Mirror.cs`
returned the fold lattice adopted at a SCALAR γ on the XY chain.
`BranchLocusPalindromeWitness` returned the reading in as many words: "the antilinear mirror
λ ↦ −λ̄ − 2σ (reflect Re about −σ, preserve Im)".

`docs/GLOSSARY.md` returned NOTHING on R₉₀, the locus, τQ, the frozen divisor, or the pinned
block. `fw.Confirmations` returned nothing; no flight touches this. `docs/CAUGHT_ERRORS.md`
returned the 2026-06-20 entry correcting the F91 family's group to Klein V₄, and nothing on
F140 against F153.

## The identity, written as a self-map

A block (p, q) carries an antilinear involution about a REAL constant c when a permutation m of
its cells satisfies, entry for entry,

    L[m i, m j] = −conj(L[i, j]) − c·δ_ij.                                    (*)

This is `Pᵀ L P = −conj(L) − c·I` with P a permutation matrix, a similarity, so
spec(L) = −conj(spec(L)) − c with algebraic multiplicities, needing no normality: these blocks
are decidedly not normal. It needs no diagonalisability either, but that is a weaker remark than
it looks and an earlier draft made it a stronger one by writing that the blocks lack both. They
do not lack diagonalisability. At N = 4, Δ = 1, locus γ, J = 1.3 the eigenvector-MATRIX condition
numbers κ(V) of the blocks run from 6 to 25 (the per-eigenvalue condition numbers, ranging over
individual eigenvalues rather than blocks, from 1.0 to 10.4),
and defective couplings are isolated points: for the corner block PROOF_R90_FROZEN_DIVISOR §9
catalogues them, and other blocks carry their own, e.g. the (1,2) real EP that
PROOF_CODIM1_BY_ADDITIVITY §7 records. With c real the pairing is Re λ ↦ −Re λ − c and
Im λ unchanged. The map m need not be an involution for any of this; involutivity only makes the
antiunitary v ↦ P·conj(v) square to the identity, and nothing below leans on it except the closing
real-forms remark, which does need the square (both sources supply it, being involutions). Source 1 below is
F89d's cross-fold in the case where the fold's image block IS the block, i.e. §7(b)'s self-folded
block; source 2 is not F89d's map at all, being the site reflection with its own two hypotheses.

Two maps supply it, and they are F153's two sources, and **they do not carry the same hypothesis
on H**. **Source 1**: if p = N/2 or q = N/2, the bitwise complement of that leg maps the block to
itself, with c = 2σ, σ = Σ_l γ_l, at any γ profile and any bond profile. **Source 2**: if γ is
anti-palindromic **and H is itself reflection-symmetric**, the site reflection on both legs, with
c = 4γ̄·|S| on a single-size-class block, which on a pinned block is c = 4γ̄·|p − q|.

The bond hypothesis on source 2 is not decoration and it is not inherited from anywhere: a bond
list that is not a palindrome breaks the identity in proportion to its asymmetry, with no
threshold: for bonds [1, 1+ε, 1+2ε] at N = 4 the residual is 4ε exactly, measured from ε = 0.25
down to 10⁻⁶ (the committed gate pins one asymmetric list per N; the ε sweep is a separate
measurement). Measured on the locus at Δ = 1 with bonds [1, 1.25, 1.5, …], the residual on the
pinned blocks is 1.0, 2.0, 3.0 at N = 4, 5, 6, while source 1 on the same chains stays exactly
0.0 at N = 4 and 6 (vacuously at N = 5, which has no index-N/2 block). F153 gates this distinction on the pinned blocks
(`SourceTwoNeedsAReflectionSymmetricH_SourceOneDoesNot`); what is new here is only that the
hypothesis must be carried into every statement below, the γ̄ = 0 stratum included, because that
stratum is this note's own result and does not inherit F153's fence.

**One fence worth naming, and it says what (*) is NOT.** Source 1 is the global complement X^N
on one leg, so it needs [H, X^N] = 0, and a LONGITUDINAL FIELD breaks it: measured at N = 4 on
the block (1,2), the residual goes from exactly 0 to 1.75 for a generic field and is still
1.00 for a reflection-symmetric one. F153's criterion meanwhile SURVIVES that field exactly;
measured at uniform γ on the same chain, every pinned block stays on −2γ|p − q| to machine zero
while the residual reads 1.75. The reason is not that the field is diagonal in the computational
basis, which an earlier draft gave and which does not follow: adding an imaginary diagonal to a
NON-NORMAL matrix moves real parts in general. Nor is it, as the FIRST repair of that draft said,
that the Hamiltonian part of L is −i·ad_H and therefore anti-Hermitian, so that no
block-preserving H can reach Re λ. That is the same overreach in the other direction and it is
false: at uniform γ and Δ = 0 a longitudinal field takes the N = 4 block (2,2) from Re range
[−4.0000, 0.0000] to [−3.8680, 0.0000], and (1,2) from [−3.0000, −1.0000] to
[−2.9464, −1.0000]. What is true is the SCOPED version, and its scope is F153's criterion case,
the pinned block at uniform γ, a strict subset of the entry's subject:
on a block whose dissipator is a real SCALAR, which at uniform γ is exactly a pinned block, L is
c·Id plus an anti-Hermitian part, and THEN no H reaches Re λ. The same field leaves the pinned
(0,2) block at −2.0000 unmoved. F153's own entry makes the point from the other side, its
criterion surviving a Peierls phase that is neither real nor diagonal-symmetric. So on the h axis the involution breaks and the conclusion does not, which
means (*) is a sufficient condition and not the mechanism of the criterion there. Both halves
are in the gate, since a control that measured only the residual would support the opposite
reading.

## Why source 2 needs one size class, and where F140's defect comes from

On the locus the reflection sends a cell's disagreement sum γ(S) to `2γ̄·|S| − γ(S)`, S the
disagreement set: an affine reflection whose CENTRE, **−2γ̄·|S| in λ**, is fixed by the cell's
size class. The factor is worth care, because the rate and λ are two variables and λ = −2·rate:
γ̄·|S| is the centre in the RATE, and the λ centre is −c/2 = −2γ̄·|S|, which is what
PROOF_R90_FROZEN_DIVISOR §5 writes and what the |S| = 2 line below reproduces. So a single c
serves a block exactly when all its cells carry the same |S|, **provided γ̄ ≠ 0**. The window is `{|p − q|, …, min(p + q, 2N − p − q)}` in steps of 2 and holds a single class exactly
when `min(p, q) = 0 or max(p, q) = N`; both the window and that equivalence are the repo's, in
three places that cross-cite each other so they cannot drift (PROOF_CODIM1_BY_ADDITIVITY §6, the
typed `DisagreementWindow`, and `WindowShellLemmaTests`).

**The identification, and it is a defect MATRIX and not a defect number.** F140's (1,1)-type
block has window {0, 2}, two classes. The |S| = 2 cells, the off-diagonal ones, force the
constant c = 8γ̄, whose axis −4γ̄ is F140's frozen root; the |S| = 0 cells, the diagonal ones,
force c = 0. One constant cannot serve both. The mismatch is not merely the same NUMBER as the
even defect `8γ̄·P_D` that PROOF_R90_FROZEN_DIVISOR §2 localises on D: it is the same defect
MATRIX, reached by a second and antilinear identity. Measured on the corner at N = 4, 5, 6 and
J ∈ {0.7, 1.3}, the whole residual of the plain-c version sits on the |S| = 0 cells at exactly
8γ̄, is exactly 0.0 on every entry off the matrix diagonal (i ≠ j), and on the matrix-diagonal
entries of the |S| = 2 cells is zero within the same error model the rest of this gate carries (0.0 at N = 4, one ZZ rounding at N = 5 and 6;
calling that column "exactly 0.0" was a repair's own overstatement), so that

    L[m i, m j] = −conj(L[i, j]) − 8γ̄·δ_ij + 8γ̄·(P_D)_ij

holds entry for entry with the same P_D. §2's identity is LINEAR and this one is ANTILINEAR, and
on that two-class corner each map satisfies exactly one of the two and fails the other, which is
what makes this an identification of the defects rather than a rhyme between two constants. Half
of that split is already committed, PROOF_R90_FROZEN_DIVISOR's gate ledger gating that τQ's
antilinear variant fails ("the mirror is linear"); the reflection's half, the SAME defect matrix
reached antilinearly by a second map, is this note's.

It is the DEFECT that is identified and not the tax: the tax of §3.1 is `dim D₋ = ⌊N/2⌋`, a
count of τQ's 2-cycles on D, and the size classes say that there IS a defect on D, never how
large the subtraction is. This note does not reach that count.

**A word on "served", because it is exact as arithmetic and unsupported as geometry.** The |S| = 2
cells forcing c = 8γ̄ is a statement about constants, and entry-wise on their submatrix it is
exactly 0.0. It is not a statement that a spectrum sits about the axis −4γ̄ there: that submatrix
is not L-invariant, the coupling to the diagonal cells reaching 2.6 = 2J at N = 4, J = 1.3. §2 says the same
thing from its own side, the corner spectrum being palindromic about no constant at all.

**And the two vanish together, which is §5's observation and our measurement.**
PROOF_R90_FROZEN_DIVISOR §5 already states it: "at γ̄ = 0 it is zero for every size class at
once, so no block needs a center of its own and the argument of Section 3 runs on **every**
block", and §5 checks it entry-wise itself. What is added here is reach rather than insight.
§5 runs two mirrors at N = 3, 4, 5, τQ on the blocks with q = p and its two-sided-X^N-bridged
partner on q = N − p; measured here with the bare
reflection on EVERY block and over Δ, (*) holds with c = 0 on all 21 blocks of dimension above
one at N = 4, all 32 at N = 5 and all 45 at N = 6, at Δ = 0 and 1, residual exactly 0.0 at the
gate's dyadic coupling J = 1.25 (at a non-dyadic coupling the stratum carries the same one
ZZ-accumulation rounding as everywhere else in this gate, 2⁻⁵⁰ at J = 1.3, an independent
builder confirming both readings), including the blocks (1,1), (1,2) and (2,2) that are controls
everywhere else. That reach
inherits source 2's bond hypothesis and does not get it from §5: with a non-palindromic bond list
the stratum breaks too, at 2.0 and 4.0 for N = 4 and 5. γ̄ = 0 is also exactly where F140's own
count changes, the diagonal cells having no tax to charge and the count becoming N rather than
⌊N/2⌋. The counterexample to a careless "exactly when" and the corroboration of the
identification are one fact seen twice.

That stratum carries a fence, and it is a physical one rather than a bookkeeping one: γ̄ = 0
together with γ ≥ 0 forces γ ≡ 0, so on this locus, where γ_{R(l)} = −γ_l pairs the sites,
any nonzero zero-mean profile puts gain on one site of every nonzero pair, up to ⌊N/2⌋ sites
(the centre site at odd N carrying zero)
and the generator is no longer completely positive. What is read there is algebra about the
same matrix, not a statement about an open system, and F140's own zero-mean stratum carries the
same caveat.

So the two results are two readings of the same reflection action rather than two tenants of
one locus. F153 reads the linear algebra, needs the action to be a single reflection of the
whole block, and is therefore confined to the single-class blocks at γ̄ ≠ 0. F140 reads an
INDEX, needs oddness on one cell set only, and carries the defect on the other, which is what
lets it produce a count where F153 produces a symmetry.

**What this does not do.** It does not make the two involutions the same map. τQ is the
transpose dressed with the reflection, `(a,b) ↦ (R(b), R(a))`; source 2's map is the bare
reflection, `(a,b) ↦ (R(a), R(b))`. On the rate diagonal the transpose is invisible, a
disagreement set being symmetric in its two legs, and both act identically there, an identity of
the arithmetic rather than a measurement. They part in three places, and an earlier draft named only the
smallest. The transpose separates them on the hop part; at Δ ≠ 0 it also moves the ZZ IMAGINARY
diagonal, by 5.2 against the reflection's 0.0 on the N = 4 corner. And the primary difference is
not a part of L at all: **τQ's identity is linear and the reflection's is antilinear**, and on
that corner, with both written about the same 8γ̄·P_D defect, each map satisfies exactly one of
the two and fails the other. And XY_FROZEN_BAND's argument stands: an antilinear involution has
equal-dimensional real forms, so (*) can never be what counts F140's rooms.

## Two scope notes on F153, and the error checking them turned up

**Source 1 is stated inside the pinned universe.** F153 introduces it through "the block's FREE
INDEX, the one of p, q that is not pinned to 0 or N", and free index is defined only for pinned
blocks. That is the entry's subject and not a mistake; the complement argument needs only one
index equal to N/2, which is what §7(b) states block-generally. Measured here entry-wise on the
enlarged set, both legs, at Δ = 0, 1, 2, on and off the locus: (*) holds for every block with an
index N/2, nine of dimension above one at N = 4 and thirteen at N = 6, of which five and nine
respectively are not pinned.

**What extends is the PALINDROME, not the flat set.** An earlier draft derived this from (*)
together with the one-sided Absorption bound, saying that symmetric about the bound and bounded
on one side forces every eigenvalue onto it. That argument is wrong twice over and both errors
are instructive, because they are the merge of *flat* and *on the floor* that F153's own entry
spends a paragraph forbidding, committed in a note about F153.

It is wrong first because the criterion needs no involution at all. At UNIFORM γ a pinned block
has one size class, so its dissipator is the real scalar c = −2γ|p − q| and L = c·Id + i·S with
S exactly Hermitian, which puts every Re λ at c with neither the locus nor (*) in the argument.
That is F153's own derivation, and (*) is one sufficient route to the same place, not the
mechanism. It is wrong second because the table below is measured on a PROFILE, where the
palindrome axis and the floor are two different numbers: on the very locus profile used here the
axis of the (0,3) block sits at −3.000 while its floor −2·min(rate) is −0.750, so the spectrum
is not symmetric about the bound and could not be forced onto it.

What the PINNED blocks do is F153's own result and is Δ-robust in the LABEL, not in the onset:
past the collapse onset they flatten onto the trace constant at every Δ measured, and the onset
itself is strongly Δ-dependent, already complete at J = 1.25 at Δ = 0 (the spread is 1.15 at
J = 1, so an onset exists there too) yet past J = 10³ at Δ = 2 (the (0,3) spread still reads 4.3
at J = 100 there). What the NON-PINNED index-N/2 blocks do, (3,3) excepted (below), turns out to depend on Δ;
this note reports the measurement and gates it, and the mechanism, found after this note
closed, is [THE_SPREAD_IS_A_RESONANCE](THE_SPREAD_IS_A_RESONANCE.md)'s.
A second draft did claim a mechanism, that a multi-class block saturates on the interval spanned
by its size-class centres −2γ̄·|S| and never flattens, and Δ = 2 refutes the second half of it.
The four measured Δ show THREE behaviours, not two; the successor note shows flat is the
generic case and the other two are resonances, this Δ set happening to contain all three
resonant points. Measured at N = 6 on the locus, at J = 10⁴ and
converged (identical at 10⁵ and 10⁶; that convergence check is a separate measurement, the
committed gate's ladder stopping at 10⁴), spread with the Re-range in brackets:

    block    window          Δ = 0            Δ = 0.5           Δ = 1            Δ = 2
    (0,3)    [3]     0.000 [−3.00]     0.000 [−3.00]    0.000 [−3.00]    0.000 [−3.00]
    (1,3)    [2,4]   2.000 [−2,−4]     0.449            2.000 [−2,−4]    0.000 [−3.00]
    (2,3)    [1,3,5] 4.000 [−1,−5]     0.566            4.000 [−1,−5]    0.000 [−3.00]
    (3,3)    [0,…,6] 6.000 [0,−6]      6.000 [0,−6]     6.000 [0,−6]     6.000 [0,−6]

The gate reads all seven (p,3) rows; the omitted (4,3), (5,3) and (6,3) repeat (2,3), (1,3) and
(0,3) value for value in its output.

**Read the J column before reading the Δ row.** A ladder stopping at J = 100 reads Δ = 2 as
broken when it has merely not converged: the (0,3) spread there is 4.300 at J = 1.25, 4.287 at
J = 100, and flat at the eigensolver's backward error (~10⁻¹⁰) only from J = 10⁴. The gate's
ladder reaches 10⁴ for that reason, and an
earlier version of it failed twenty-eight rows by stopping short. Δ = 0.5 is NOT a slow
version of the same thing: 0.449 and 0.566 are stable to five digits from 10⁴ to 10⁶, so that
row is a third behaviour and not a transient. The table's VALUES are γ̄-proportional and its
READING is γ̄-independent: at γ̄ = 0.75 every value is scaled by exactly 3/2 and every label is
unchanged, so the bracketed Re-ranges above are the γ̄ = ½ ones.

Three readings, and each is asserted only where it is measured. Where a block SATURATES, it
saturates on the interval spanned by its size-class centres, endpoints included: at γ̄ = 0.5 the
(1,3) ends land on −2.0000 and −4.0000, its centres exactly, and at γ̄ = 0.75 on −3.0000 and
−6.0000. That is what makes it a size-class reading rather than a width that happens to fit, and
γ̄ is swept because an earlier draft read the width as "twice the window half-width", which is
the width in the RATE variable and agrees with the λ one only at γ̄ = ½, where the whole table
happened to be run. The saturation mechanism lives two notes downstream, in
[THE_ENDPOINTS_ARE_A_DENSITY_LAW](THE_ENDPOINTS_ARE_A_DENSITY_LAW.md): wherever the colliding
eigenspaces carry reflection-symmetric compressed site densities, the locus compression is
−2γ̄·Π N_XY Π (N_XY the disagreement-count operator of F122) and the centres are its Rayleigh
bounds, attained by pure-size-class vectors; at Δ = 0 and ½ parity forces that on every
colliding eigenspace of the blocks in that note's census, and at Δ = 1 a few parity-mixed
eigenspaces break the density symmetry in the INTERIOR (measured on (1,3) and (2,4), the
second not a row of this table) and stay inside the centres. Where a block FLATTENS it sits on the trace constant −2γ̄(p + q − 2pq/N),
which is −3.00 for every one of these blocks. And where it does NEITHER, at Δ = 0.5, the note
has no name for what it does beyond the number: the spectrum sits about the trace constant with
a residual spread that is neither zero nor the size-class width; what selects it is the
successor note's answer (a small resonance, rank-2 collisions with an exact double level). And (3,3), the only block with BOTH indices at N/2, saturates at every Δ measured.

**Why Δ selects among the three was this note's sharpest open question, and
[THE_SPREAD_IS_A_RESONANCE](THE_SPREAD_IS_A_RESONANCE.md) answers it**: there is no selection.
Flat is the generic case, exact and profile-free (an X^N density cancellation on the
half-filling leg), and the spread is a RESONANCE of the sector difference spectrum, this
four-point table having sampled the three resonant Δ and one generic one. The observation
this paragraph used to close on, that the exempt block (3,3) is exactly the one where source 1
applies on BOTH legs, is the p = q = N/2 special case of the successor's wider fact: any p = q
gives the standing i = j collision at every Δ.

A second connection suggested itself and turned out to be false, which is worth recording
because it was one sentence away from being written as support. F153's entry names Δ = 0 and
Δ = 0.5 as the values where the pinned blocks' span SATURATES instead of collapsing, and those
are two of the four Δ measured here, with Δ = 0.5 carrying this table's third behaviour. The
rhyme is not one. F153's reading is OFF the locus and this table is ON it, and on the locus at
J = 10⁴ every pinned block collapses at every Δ measured, including 0 and 0.5, on both the dyadic
ramp and linspace(0.1, 1, N) (the linspace half is a separate measurement, not in the committed
gate, whose profiles are the locus, off-locus and zero-mean triple). Checking that is what turned
up the missing clause in F153's entry,
now repaired there together with an N clause the Δ = 0.5 half also needs.

So F153's flat set is exactly right and
unchanged, and its gate `TheGenericFlatSet_IsTheFreeIndexAtHalfFilling_AndNotAParity` is
correct as written: both sides of its equality range over pinned blocks, so nothing here can
falsify it. An earlier draft of this note claimed that assert was stale. It is not, and the
claim is withdrawn rather than carried forward as work.

## The gate

`simulations/two_sources_gate.py`, writing `simulations/results/two_sources_gate.txt`. It pins
(*) with `== 0.0` on dyadic γ; sweeps N = 4 and 6 for source 1 (odd N has no index-N/2 block at
all) at Δ = 0, 1, 2, and N = 4, 5, 6 for source 2 at Δ = 0, 1, the latter over a locus profile,
an off-locus profile and the zero-mean stratum; tests both legs of source 1; and runs the
Δ-behaviour table at N = 6 over Δ ∈ {0, 0.5, 1, 2}, J ∈ {1.25, 100, 10⁴} and γ̄ ∈ {0.5, 0.75},
each row's label computed and compared with the committed one.

The controls. Source 2 off the locus and on the non-pinned blocks (1,1), (1,2), (2,2) must
break, and there the gate scans c and prints the SMALLEST residual, so the row shows that no
constant serves the block UNDER THIS MAP rather than only that the claimed one does not. The
scoping matters and an earlier draft dropped it: the scan is over constants at a fixed
permutation, so it cannot speak for other permutations, and the map-independent statement needs
the spectrum instead of the entries. The scan
window is derived from the candidate constants and a boundary hit fails the gate, an earlier
fixed window having sat below 2σ so that every printed minimum was its own edge. Two values of
γ̄ are swept, since at a single γ̄ the constant 4γ̄·|p − q| cannot be told from one with the γ̄
dropped, and the saturation table is swept over γ̄ too, which it was not until a review round
found that its law was a γ̄ = ½ reading. Source 2's other hypothesis, that H itself is
reflection-symmetric, has its own control now: a non-palindromic BOND LIST must break the
identity and does, at 1.0, 2.0, 3.0 for N = 4, 5, 6 on the pinned blocks and at 2.0 and 4.0 on
the zero-mean stratum, while source 1 on the same chains stays exactly 0.0 where it exists
(N = 4 and 6; N = 5 has no index-N/2 block). Until that row
existed the entire source-2 section ran at uniform J and could not see the hypothesis it stands
on. Source 1's control is the longitudinal field, in both of its halves. The
linear-with-gauge reading, which is `Mirror`'s leg, is gated in both directions, exactly 0.0 at
Δ = 0 and required to FAIL at Δ ≠ 0, which is §7(c) measured; that failure is gauge-independent,
the whole residual sitting on the diagonal where the gauge squares to one. τQ and the bare
reflection agree exactly on the rate diagonal by arithmetic, a disagreement set being symmetric
in its two legs; that is an IDENTITY the gate states rather than measures (an earlier version
dressed it as a gated row whose failure branch no input could reach), and it is what lets the
two results share one reflection action. What the gate measures is two of the three places they
part, the ZZ imaginary diagonal at Δ ≠ 0 and the linear-against-antilinear split (the third, the
hop part, is where the transpose acts and is stated, not gated); the split is gated in both
directions so that a map satisfying BOTH identities would fail the gate rather than pass it. The F140 defect is gated as a MATRIX and not as a number: the
plain-c residual is required to be exactly 8γ̄ on the diagonal cells and within the error model
everywhere else, so a residual of the right SIZE in the wrong PLACE would not pass.

One trap the gate caught in its own construction and now documents: an arithmetic ramp is
automatically anti-palindromic, so linspace and every dyadic ramp lie ON the locus and cannot
serve as an off-locus profile. Every profile's class is asserted rather than commented.

The spectral consequence is checked against a dense eigendecomposition as a ratio to the
block's spectral radius. It is asserted only where it has POWER, and there is one family where
it has none: on the blocks (0, N/2), (N/2, 0), (N/2, N) and (N, N/2) the two constants coincide identically, since
|p − q| = N/2 makes 4γ̄·|p − q| = 2Nγ̄ = 2σ at any profile, so source 1 supplies the pairing
off the locus as well and those rows would read as carried whatever the reflection did. They
are marked in the output and excluded from the spectral assertion; their entry-wise verdict
still counts. The eigenvalue sort keys on rounded (Re, Im) jointly, since a naive complex sort
scrambles near-tied real parts and fabricates a large defect.

**One residual is bounded and then read.** Source 2 at Δ ≠ 0 carries one rounding of the ZZ
diagonal, 2⁻⁵¹ at N = 5 and 2⁻⁵⁰ at N = 6 at the gated couplings (J = 2.7 reaches 2⁻⁵⁰ and
2⁻⁴⁹, printed in the gate's own rounding ladder). The mechanism is the ORDER in which that diagonal
accumulates over bonds, not the coupling's UNREPRESENTABILITY: J = π and J = 1.1 are as
unrepresentable as 1.3 yet give exactly 0.0, while 1.3, 0.9 and 2.7 do not, and N = 4 is exactly
0.0 at every coupling tried, which the gate computes rather than asserts in prose. So the residual is a deterministic
function of an input the identity's value does not contain and vanishes for settings of it,
which is this repo's test for reading a rounding rather than gating it. It is nonetheless
bounded by 4·eps·max|L| and fires above that: an earlier version gated only the dyadic
couplings and printed the rest, and an injected O(1) break at the non-dyadic coupling passed
unnoticed.

## Which modes sit on the axis below the threshold, and why it is not a new law

*Added 2026-08-22, from the class-B fencing pass, which needed to know what a fence may say about
this block under a profile. Recorded as a clause of [F153](../docs/ANALYTICAL_FORMULAS.md#f153),
NOT as a number of its own: a draft of this section was minted as "F156" and withdrawn the same
hour, because everything it rests on is committed and the entry would have been a corollary in a
narrower special case than the lemma it quotes.*

The section above gives (\*) and reads off, for real c, that the pairing is `Re λ ↦ −Re λ − c` with
**Im λ unchanged**, with algebraic multiplicities and without normality. F153 adds that the resulting
symmetry about the trace constant holds at EVERY real coupling, while flatness, all eigenvalues ON
the axis, sets in only above one. What was not written down is the step between those two, and it is
one line:

> an eigenvalue whose imaginary part is shared with no other DISTINCT eigenvalue of the block is its
> own image under the pairing, hence sits on the axis exactly, threshold or no threshold.

Flatness is then simply the case where every mode is Im-lone. On the (0,1) block the axis is `−2γ̄`.

**This is an owned argument read on a new object, not a new argument.**
[`F89BranchLocusPalindromeClaim`](../compute/RCPsiSquared.Core/Symmetry/F89BranchLocusPalindromeClaim.cs)
makes exactly this step for the antiunitary `λ ↦ −λ̄ − 2σ`: "reflect Re about −σ, preserve Im, so
every EP lies on the line or in a mirror pair across it, **no orphan**". Same step, different pairing.
Whoever needs it elsewhere should look there first.

**Measured** on (0,1) at N = 7 on the locus, Δ = 0, J real: three of seven modes are Im-lone and on
the axis to 10⁻¹⁵ while the block's span is 0.87, staying there through 0.66, 0.49, 0.12 as J rises,
until every mode is lone and the span is zero. At Δ = 1, N = 7, J = 1 five of seven are lone and on
the axis to 4·10⁻¹⁵ at a span of 1.08. Gate:
[`band_edge_profile_fence_gate.py`](../simulations/band_edge_profile_fence_gate.py).

**Four fences, each of which a draft of this section got wrong.** The coupling must be **real**: at
N = 7 with `J = 1 + 0.3i` the identity residual is 1.20 and the "lone" modes miss the axis by 1.19,
which is F153's own REAL-coupling fence and not a new one. The identity's form is (\*)'s
`Pᵀ L P = −conj(L) − c·Id`, **not** `R M R = −M† − c`: those agree only because a real hopping makes M
symmetric, and on a Peierls chain at N = 6 the dagger form misses by 1.6 where (\*) holds at
2·10⁻¹⁶. Below the threshold the statement is often **vacuous** rather than false, at J = 0 always,
the block being real diagonal with every Im equal. And the hypotheses are **sufficient, not
necessary**: an N = 4 ring whose rates leave the locus (R₉₀ defect 0.25) fails (\*) at 0.50 and still
has every mode exactly on the axis, to 10⁻¹⁶, so a failed identity forbids nothing.

**A sharpening of source 2's own hypothesis, found by the gate rather than by reading.** Source 2 asks
for `H` reflection-symmetric. What (\*) actually needs is `P A P = Aᵀ`, and the two are not the same:
a flux ring with hopping phase 0.7 has `‖R H R − H‖ = 1.29` and yet `‖P A P − Aᵀ‖ = 0` exactly, so
(\*) holds for it at 10⁻¹⁶ and every mode is on the axis. For a real symmetric hopping the weaker
condition reduces to source 2's, which is why the difference never showed. A review round read the
flux ring as a counterexample to the identity; it had measured the dagger form, and the ring is
INSIDE the theorem once (\*) is used.

**The loneness predicate needs a tolerance, and the axis residual is a law rather than a number.** At
N = 4, γ = (0.9, 0.65, 0.35, 0.1), J = 0.3 two genuinely degenerate Im values differ by 7·10⁻¹⁶ in
float, so a strict inequality declares all four modes lone and the claim appears to fail by 0.304. The
test needs a tolerance at or above `eps·|λ|`, and the residual on the axis scales as `O(1)·eps·|λ|`, so
a fixed absolute bound breaks at large J rather than at large N.

**One direction, and the converse fails only trivially.** Degenerate Im permits a spread and does not
force one. The witnesses found are at uniform γ, where the diagonal is already the scalar `−2γ·Id` and
the conclusion holds for an elementary reason that never touches the pairing (a star, `N − 2` modes at
`Im = 0`). Whether a genuinely non-uniform on-locus profile can be pinned with a degenerate Im
spectrum is **open**: a local search of 48 000 draws over chain, ring and complete at N = 4..7 and four
couplings found none, and that search is a scout rather than a committed gate.

Not to be merged with two neighbours. It is not the route this note rejects above, which argued from
the palindrome plus the one-sided Absorption bound and merged *flat* with *on the floor*; nothing here
mentions the floor, and on the profile measured the axis sits 0.875 away from it. And it is not
[THE_SPREAD_IS_A_RESONANCE](THE_SPREAD_IS_A_RESONANCE.md)'s `X^N` cancellation, which needs an index at
`N/2` that (0,1) never has for `N > 2`.

## What is still open

The bond profile on the enlarged block set. F153 gates that source 2 needs H
reflection-symmetric while source 1 does not
(`PinnedBlockFloorClaimTests.SourceTwoNeedsAReflectionSymmetricH_SourceOneDoesNot`), so the
question is closed on the pinned blocks, and the gate here now carries its own bond control on
the pinned blocks and on the zero-mean stratum. What is still not measured is the non-pinned
index-N/2 blocks under a bond profile: source 1 is blind to it by the argument and by the
control, but the reflection's behaviour there is untested.

And the count. This note says why F140's argument may carry a defect on its second size class
and F153's may not; it does not derive ⌊N/2⌋, which stays where PROOF_R90_FROZEN_DIVISOR §3.1
has it, as fixed points on the off-diagonal cells minus 2-cycles on the diagonal ones.
