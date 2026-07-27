# What reaches the ceiling: the measuring, in the order it happened

This is the measurement side of **F145 and F146**, proved in
[PROOF_SCALAR_COUNT](../docs/proofs/PROOF_SCALAR_COUNT.md). It is kept in the order
the afternoon actually went, including the two readings that were written down as if
true and killed within the hour, because that is where the obvious interpretation
breaks and the next person deserves to see it.

One word does double duty below and the repo spends it elsewhere: **block** on its
own is the joint-popcount block (p, q), while **2-block** and **3-block** always mean
the connected maximizers, which turn out to be the metric and the volume of SO(3).

The title says ceiling and F144 says floor, and they are the same states: the
disagreement 𝒦 and the double occupancy D̂ add to ℓ on a diagonal rung, so the floor
of what the watching charges is the ceiling of what survives it.

The result, in one line: the states attaining the F144 disagreement floor on ℓ chiral
pairs are the rotation-invariant couplings of one spin-1 per pair, so there are
C(⌊N/2⌋, ℓ)·R_ℓ of them with R the Riordan number, and they are products of two
connected blocks. The note was written and carried UNNAMED for a day, on the rule that
we cannot name what we do not yet know; the name came when the mechanism did.

**Date:** 2026-07-27, in three passes: the count in the afternoon, the reach in ℓ and
the resonance criterion in the evening, the mechanism and the proof last.
**Grew out of:** [ETA_CEILING_REDUCTION](ETA_CEILING_REDUCTION.md) open item 1, the
attainment. F144 ([PROOF_FROZEN_BAND_SO4](../docs/proofs/PROOF_FROZEN_BAND_SO4.md)
Section 7) proves the floor 𝒦 ≥ ℓ(N − ℓ)/(N + 1) on LW_ℓ ∩ V₀. That the floor is
ATTAINED, and by how many states, is what this note measures.
**Reproduction:** [`simulations/scalar_count.py`](../simulations/scalar_count.py),
67 checks in about a minute and 78 under `--deep`, which adds the rungs ℓ = 9 and
ℓ = 10, the M = 44 prediction, AND a wider range of N in the checks W2, W3, W5 and W6.
It carries the load-bearing numbers: the count at every rung, the blocks against the
two invariants, the relations against the classical ones, the triplet sector, and the
resonances. It does NOT carry every number quoted below. Section 5's cell counts,
Section 8's cutting table, the ℓ = 2 sweep to M = 84, the 235 triples against the
solved block, and the ℓ = 9 relation comparison all come from the exploration; an
independent audit recomputed them.

---

## The question

By F144 the disagreement on a lowest-weight, turning-blind state cannot fall below
ℓ(N − ℓ)/(N + 1). Gate check V12 of `simulations/eta_ceiling_reduction.py` says the
bound is attained on every rung it reaches, and that the saturating space is exactly
what Corollary 7.3 describes. What it also says, and what this note starts from, is
that the DIMENSION of that space is C(⌊N/2⌋, ℓ) only up to ℓ = 3 and larger from
ℓ = 4 on. The committed note calls that open and leaves it. So: how many states reach
the ceiling, and what are they made of?

## 1. The count factorizes, and the pair choice is one factor

Write p := ⌊N/2⌋ for the number of chiral pairs {a, M − a}. Measured dimensions of
the saturating space:

| N | p | ℓ = 2 | ℓ = 3 | ℓ = 4 | ℓ = 5 | ℓ = 6 |
|---|---|---|---|---|---|---|
| 6, 7 | 3 | 3 | 1 | | | |
| 8, 9 | 4 | 6 | 4 | 3 | | |
| 10 | 5 | 10 | 10 | 15 | 6 | |
| 12 | 6 | | | 45 | 36 | 15 |
| 13 | 6 | | | 45 | | |

Divide by C(p, ℓ) and the table collapses to one row:

  **dim = C(⌊N/2⌋, ℓ) · m(ℓ),  m(ℓ) = 1, 1, 3, 6, 15 for ℓ = 2, 3, 4, 5, 6.**

The second factor is free of N: m(4) = 3 reads the same at p = 4, 5 and 6, and
m(5) = 6 at p = 5 and 6. So WHICH chiral pairs a maximizer uses is a free choice,
and HOW MANY maximizers live on a given choice is a pure function of the rung.

Two supporting facts, both from the support of the space in the mode-cell basis:

- **Every maximizer touches exactly ℓ chiral pairs.** Not fewer, not more.
- With ℓ pairs and ℓ particles the occupancies are forced: if d pairs are doubly
  occupied then exactly d are empty and ℓ − 2d are singly occupied, since the pairs
  touched number ℓ. A full pair carries Slater energy zero, so V₀ couples precisely
  the states that differ in WHICH of the remaining 2d pairs are full. The observed
  occupancy patterns are exactly the admissible d: at ℓ = 4 the patterns (1,1,1,1),
  (2,1,1,0) and (2,2,0,0) all appear, at ℓ = 5 only (2,1,1,1,0) and (2,2,1,0,0), the
  all-singly-occupied pattern being absent at odd ℓ.

## 2. What they are made of: two connected blocks and products

On a fixed choice of ℓ pairs the maximizers are spanned by PRODUCTS of connected
blocks, and connected blocks come in exactly two sizes.

- **The 2-block** is the closed-form ℓ = 2 maximizer that
  [ETA_CEILING_REDUCTION](ETA_CEILING_REDUCTION.md) already carries, unchanged, six
  terms with entries in {1, 2}.
- **The 3-block** has to exist because three pairs admit no splitting into blocks of
  size ≥ 2, and at ℓ = 3 it is the only maximizer. It is now in closed form too:
  twelve terms, every entry ±1, and the cell pattern is one sentence. **One chiral
  pair sits FULL on the left, another sits FULL on the right, and a single mode of
  the third pair is shared by both.** Writing the three pairs in the order of their
  smaller mode, P_i = {a_i, ā_i} with a_i < ā_i,

      T₃ = Σ_{i ≠ j} Σ_{m ∈ P_k} sgn(j − i)·σ(m)·|P_i ∪ {m}⟩⟨P_j ∪ {m}|,

  where k is the remaining index, σ(a_k) = +1, σ(ā_k) = −1, and the cells are read
  in the sorted-mode fermionic convention. Six ordered (i, j) times the two modes of
  P_k give the twelve terms. The 3-block is transpose-ANTIsymmetric where the
  2-block is symmetric, and it carries no diagonal cell at all. Checked against the
  SOLVED block, the one-dimensional exact rational nullspace on the triple's six
  modes, on all 235 triples at N = 6, 7, 8, 9, 10, 12, 13, 15, 16, 18, without
  exception. That check belongs to the exploration; what the gate checks is the same
  closed form against the invariant volume, on 83 blocks by default and 258 under
  `--deep`.

Everything above that is a product, and the products are enough:

- At ℓ = 4 the three products of two 2-blocks (the three pairings of four pairs) are
  linearly INDEPENDENT and span the whole 3-dimensional space at N = 8. No connected
  4-block is needed.
- At ℓ = 5 all ten products (2-block on a pair of pairs, 3-block on the complementary
  triple) lie in the space and span it, rank 6 of dimension 6. No connected 5-block
  is needed.

Fermionic products are taken with the merge sign of sorting the two mode sets.

## 3. The count is the Riordan sequence, and the test that says so

1, 1, 3, 6, 15 continues, if it is the Riordan sequence A005043, with 36 and 91. The
obvious competitor is the double factorial (ℓ − 1)!! = 1, 3, 15, 105 for even ℓ,
which counts the perfect matchings of the ℓ pairs and agrees with Riordan up to
ℓ = 6. **They part at ℓ = 8: 105 against 91.** So the rank of the matching products
at ℓ = 8 decides which sequence the multiplicity follows, and it can be read without
touching a Liouvillian block, since the products are built from the 2-block alone.

- ℓ = 2, 4, 6: 1, 3, 15 matching products, exact rank 1, 3, 15. Independent.
- **ℓ = 8: 105 matching products, exact rational rank 91** (fraction-free Bareiss on
  the integer Gram matrix). Fourteen relations appear.
- ℓ = 10: 945 matching products, rank 603 over two large primes, matching
  R_10 = 603. A rank mod q can only drop, so this is a lower bound. Read at N = 20,
  which Section 9 now knows to be RESONANT, so it is a statement about the products
  and not about that space; Section 4 reads ℓ = 10 again at N = 21.

So the multiplicity per pair-choice follows the Riordan numbers, and the first
relation among the products appears at ℓ = 8 rather than at ℓ = 4, which is worth
noticing: whatever produces the relations is not a local identity on four pairs.

## 4. The extrapolation, removed: the dimension of the TRUE space to ℓ = 10

What Section 3 measures is a rank of products. That it still equals the DIMENSION of
the maximizer space is what this section had to leave open, and does not any more.
The dimension is squeezed between two computations that are one-sided in directions
that compose, and they meet at every rung:

- **from below**, the block products, each verified over the integers to satisfy
  every Corollary 7.3 condition, with a rank mod p. A rank mod p can only drop, so
  it is a lower bound on the dimension;
- **from above**, the NULLITY of the whole condition system, also mod p. A nullity
  mod p can only grow, so it is an upper bound.

| ℓ | read at N | products built | lower | upper | dimension | R_ℓ | (ℓ−1)!! |
|---|---|---|---|---|---|---|---|
| 6 | 12 | 25 | 15 | 15 | **15** | 15 | 15 |
| 7 | 15 | 105 | 36 | 36 | **36** | 36 | . |
| 8 | 16 | 385 | 91 | 91 | **91** | 91 | 105 |
| 9 | 18 | 1540 | 232 | 232 | **232** | 232 | . |
| 10 | 21 | 7245 | 603 | 603 | **603** | 603 | 945 |

So at ℓ = 8 the true dimension is **91 and not 105**: the sequence is Riordan, read
on the space itself rather than on a rank of products. Two further things fall out
of the same table.

- **The products SPAN at every rung tested**, since the lower bound reaches the
  upper one. The block model of Section 2 is therefore complete to ℓ = 10, where
  before it was verified to ℓ = 5. At ℓ = 10 that means 7245 products of rank 603:
  the model is nowhere near a basis, and it is still exactly the space.
- **The pair-choice factorization survives past ℓ = 6.** At N = 16 with p = 8 the
  rung ℓ = 7 has dimension 288 = C(8, 7)·36 exactly, so C(⌊N/2⌋, ℓ)·m(ℓ) is not an
  artifact of the small rungs.

Every number in the table is read at a NON-resonant N; see Section 9, which now
carries a criterion for which N those are.

## 5. Why this is cheap, and what makes it exact

The whole computation rests on one structural fact. Every Corollary 7.3 condition is
built from X(a,b): (A,B) ↦ (A\a, B\b), which is INJECTIVE on the cells it does not
annihilate, so every ROW of X(a,b) carries exactly ONE entry. Hence the conditions
are not linear algebra at all until the last step:

- **X(a,ā)v = 0 deletes coordinates.** Every cell with a in A and ā in B is zero.
  At N = 16, ℓ = 8 this is 97444 cells of V₀ down to 12870.
- **(X(a,a) + X(ā,ā))v = 0 pairs the survivors**, two cells and a sign per equation,
  which is a signed union-find. 12870 cells down to 1107 free variables. A sign
  clash around a cycle sets a whole component to zero.
- **Only the spin ladder S⁻ is left**, a sparse system on those few variables, and
  its rank is what the primes are for.

Two guards were run rather than assumed. The lowest-weight condition Ψv = 0 needs no
separate imposition, because Ψ = Σ_a X(a,a) and the pairing makes the sum vanish;
adding the Ψ rows anyway changes no nullity, which the gate measures at every rung to
ℓ = 7. And the rank is computed by two independent paths, a sparse exact elimination
at a prime near 2³¹ and a dense one in float64 at a prime near 10⁶; the gate runs both
and compares them to ℓ = 6, and during the exploration they were also run against each
other at ℓ = 9, where both returned 232.

## 6. Why Riordan: the maximizers ARE the multilinear invariants of ℓ vectors

The count, the two block sizes, and the relations all come from one place, and it is
classical. Each chiral pair carries a spin-1 triplet, and a maximizer on ℓ pairs is
an SO(3)-invariant multilinear form in those ℓ triplets.

**The triplet.** Write P = {a, ā} for a chiral pair. Three cell patterns survive the
two chiral conditions on that pair:

    u₊(P) = both modes of P in A,   u₋(P) = both modes in B,
    u₀(P) = (a in A and in B) − (ā in A and in B).

That is exactly the 4 → 3 collapse the reduction of Section 5 performs: the deletion
kills the mixed cell a|ā and the signed pairing glues a|a to ā|ā.

**The two blocks are the two invariants of SO(3), in that basis.** Verified cell by
cell over the integers on all 258 blocks at N = 6 to 16:

    2-block(P, Q) = u₀(P)u₀(Q) − 2[u₊(P)u₋(Q) + u₋(P)u₊(Q)]      exactly,
    3-block(P, Q, R) = −ε(P, Q, R)                                exactly,

where ε is the antisymmetric form that is nonzero exactly when the three m-values are
+1, 0, −1 in some order, with the permutation sign. So the 2-block is the invariant
METRIC and the 3-block is the invariant VOLUME, in a basis whose ± components carry
weight 2. In that normalization the classical syzygy reads ε·ε = −¼·det(Gram), so the
3-block is (i/2) times the ordinary triple product, and a product carrying t
three-blocks is off from the ordinary invariant by (i/2)^t.

**The relations are the classical ones, tested as a prediction.** Evaluate the same
partitions as products of dot products and triple scalar products of ℓ random vectors
in ℝ³ and read the relations among them. Then compare:

| ℓ | products | rank ours / classical | relations ours / classical | joint span |
|---|---|---|---|---|
| 4 | 3 | 3 / 3 | 0 / 0 | . |
| 5 | 10 | 6 / 6 | 4 / 4 | 4 |
| 6 | 25 | 15 / 15 | 10 / 10 | 10 |
| 7 | 105 | 36 / 36 | 69 / 69 | 69 |
| 8 | 385 | 91 / 91 | 294 / 294 | 294 |
| 9 | 1540 | 232 / 232 | 1308 / 1308 | 1308 |

The joint span equals each kernel, so the two relation spaces are THE SAME subspace,
not merely spaces of the same size. The one free constant, the weight (−4)^{t/2} that
compensates the (i/2)^t above, was fitted at ℓ = 6 alone (all ten relations returned
−4 to ten digits) and then used unchanged at ℓ = 7, 8, 9, where it is a prediction.
The gate carries this comparison to ℓ = 7 by default and ℓ = 8 under `--deep`; the
ℓ = 9 row is from the exploration.

**What this answers.** The three open items were one item.

- **Why Riordan.** R_ℓ is the dimension of the SO(3) invariants of ℓ vectors, that is
  the singlet count of the ℓ-fold tensor of the spin-1 representation. The count was
  never a coincidence with a singlet number; it IS one, once the singlets are taken
  multilinearly, one triplet per chiral pair.
- **Why the blocks stop at three.** First fundamental theorem for SO(3): the
  invariants of vectors are generated by the dot products and the triple products.
  There is nothing of degree four to find.
- **What the relations are.** Second fundamental theorem: the syzygies, whose
  multilinear part starts with the vanishing 4 × 4 Gram determinant. That is why the
  first relation among matchings appears at ℓ = 8 and not at ℓ = 4: four vectors in
  three dimensions need four more slots to state their dependence.

**The group is the arc's own ladder, not an isomorphic copy.** The three patterns sit
in the blocks (2,0), (1,1) and (0,2), which is where 𝔖± moves, and the ladder acts on
them as a spin-1 multiplet:

    𝔖⁻ u₊ = ∓u₀,   𝔖⁻ u₀ = ∓2u₋,   𝔖⁻ u₋ = 0,   and 𝔖⁺ the mirror image,

checked on every chiral pair at N = 6, 8, 9, 10, 12, 15. The coefficients 1 and 2 rather
than √2 and √2 are the same basis normalization that puts the weight 2 on the ±
components of the metric; their product is 2, which is the spin-1 invariant. So the
SO(3) of this section IS the F142 SU(2), one triplet per chiral pair, and the whole
identification is group-backed rather than numerological.

**Where the argument stood when this was written**, and it is worth keeping because
it is what pointed at the last step: the blocks ARE the invariants (exact, cell by
cell), and the invariants of ℓ vectors are generated in degrees 2 and 3 with the
classical syzygies (theorem), so the only measured link left was that the products
SPAN. The next two subsections close it for three families of N, and
[PROOF_SCALAR_COUNT](../docs/proofs/PROOF_SCALAR_COUNT.md) carries the finished
argument.

### Where the spanning comes from, and where it breaks

The argument for spanning is a chain of three steps,

    (A) every maximizer lives in the TRIPLET sector, one pattern per chiral pair,
    (B) inside that sector the conditions say exactly "singlet",
    (C) the singlets are spanned by products of the metric and the volume,

of which (B) is representation theory and (C) is the first fundamental theorem. Only
(A) was ever measured, and it is now measured sharply enough to see its cause. Call a
cell a TRIPLET cell when every chiral pair contributes u₊, u₋, one of the two u₀
patterns, or nothing.

**At a chiral-only rung the restriction removes NOTHING.** After the chiral deletion,
every surviving cell of V₀ is already a triplet cell. The counts before and after the
restriction are equal at every rung tested: 80 at (N, ℓ) = (8, 3), 350 at (10, 4),
1050 at (12, 4), 924 at (12, 6), 3432 at (15, 7) and 12870 at (16, 8).
So (A) is not a further condition to impose; it is what the deletion and the energy
classes have already left. And the reason is arithmetic, not dynamics: a chiral pair
sums to zero energy, so swapping one full pair for another is the ONLY way two
different mode sets can carry the same Slater energy, unless that N has extra cosine
coincidences.

**At a resonant N the restriction removes exactly the surplus**, and what is left is
the contract:

| N | ℓ | cells | dim | triplet cells | triplet dim | C(⌊N/2⌋,ℓ)·R_ℓ |
|---|---|---|---|---|---|---|
| 11 | 2 | 72 | 14 | 60 | 10 | 10 |
| 11 | 3 | 280 | 22 | 200 | 10 | 10 |
| 11 | 4 | 490 | 27 | 350 | 15 | 15 |
| 14 | 2 | 150 | 29 | 126 | 21 | 21 |
| 14 | 4 | 3570 | 213 | 2450 | 105 | 105 |
| 17 | 2 | 192 | 36 | 168 | 28 | 28 |
| 20 | 2 | 282 | 49 | 270 | 45 | 45 |
| 20 | 3 | 3000 | 246 | 2400 | 120 | 120 |

The triplet part is the law's value in every case, and the whole surplus sits on the
non-triplet cells, which exist only because that N carries the extra coincidences of
Section 9. **So the extra maximizers at a resonant N are a different animal, not more
of the same**, and the law is not violated there so much as accompanied.

That moves the open item. "Why do the products span" is no longer a blind measured
step; it is the statement that at a chiral-only rung the only ℓ-subset energy
coincidences compatible with the deletion are the chiral ones. That is a question
about vanishing sums of cosines, which is exactly the ground
[F89_SEED_EXISTENCE_REDUCTION](F89_SEED_EXISTENCE_REDUCTION.md) already stands on.

### And for a prime M it is a theorem

The cosine statement is not open everywhere. Written as a vanishing sum it is the
sibling of [PROOF_F129_LEVEL_COLLISION_LAW](../docs/proofs/PROOF_F129_LEVEL_COLLISION_LAW.md)
§2, in the same idiom and with the same imported lemmas. Two ℓ-subsets A, B of equal
Slater energy give, on their DIFFERENCES X = A \ B and Y = B \ A, since the shared
modes cancel,

    Σ_{x ∈ X} (ζ^x + ζ^{−x}) − Σ_{y ∈ Y} (ζ^y + ζ^{−y}) = 0,   ζ = e^{iπ/M},

and moving the minus signs through ζ^M = −1 makes it a vanishing sum of 2M-th roots
of unity with nonnegative coefficients, of weight 4·|X|, whose **exponents avoid 0 and
M** (each term is ±x or ±y + M with 1 ≤ x, y ≤ M − 1) and whose terms are pairwise
DISTINCT precisely because of the deletion, which forbids a mode in A whose chiral
partner is in B.

Now let M = p be prime. Then 2M = 2p has two prime divisors, so by (LL2) every
minimal piece is a rotated 2-cycle or a rotated p-cycle.

- A rotated **p-cycle** is a coset of μ_p: its exponents step by 2, so they are ALL
  residues of one parity class mod 2p. The even class contains 0 and the odd class
  contains p, and the sum has neither. **No p-cycle can occur.**
- A rotated **2-cycle** is an antipodal pair {ζ^a, ζ^{a+M}}. In our exponent set the
  only antipodal pairs available are a mode together with its chiral partner INSIDE X
  or inside Y: the cross cases are excluded, one by disjointness and one by the
  deletion.

So the whole sum decomposes into chiral pairs, which says exactly that the two
DIFFERENCES are unions of whole chiral pairs; the shared modes never entered the
argument. **For M prime there is no non-chiral coincidence at any rung**, step (A) holds, and with (B) and (C) the maximizer count and the spanning are
theorems there rather than measurements.

Checked against the search anyway, since a clean argument can still be a wrong one:
CLEAN to the rung cap at M = 7, 11, 13, 17, 19, 23 and to rung 4 or 5 at
M = 29, 31, 37, 41, 47, with no exception. The family is N = p − 1, which contains
three of this note's own readings: ℓ = 6 at N = 12, ℓ = 8 at N = 16, ℓ = 9 at N = 18.

### M = 2p closes too, and the deletion is what closes it

Take M = 2p, so 2M = 4p, still two prime divisors, so (LL2) still leaves only 2-cycles
and p-cycles. Two things change: the mode M/2 = p is now SELF-chiral and carries zero
energy, and a p-cycle is now a full residue class **mod 4** rather than a parity
class. The classes 0 and 2 are excluded as before (they contain the exponents 0 and
M = 2p, which the sum never has), so a p-cycle would have to fill class 1 or class 3
completely, each of which has exactly p members.

Count what can fill them. Only the ODD modes contribute to the odd classes, there are
exactly p of them (1, 3, …, 2p−1), and each one used contributes exactly ONE exponent
to class 1 (and one to class 3). So a p-cycle forces every odd mode to be used exactly
once, and then a short calculation pins the sides. Writing m̄ = 2p − m,

    (m in X) and (m̄ in Y) give the SAME class-1 exponent, and so do (m in Y) and (m̄ in X).

The exponents are pairwise distinct, by the deletion, so each class-1 exponent is hit
once; covering both exponents of a chiral pair {m, m̄} therefore forces one choice from
each colliding group, and the only combinations left put **m and m̄ on the same side**.
The even modes cannot enter a p-cycle at all, so they pair off in 2-cycles, which is
the same conclusion. Either way X and Y are unions of whole chiral pairs.

A round of review then found the argument's own gap and made it shorter. The mode
M/2 = p is its OWN chiral partner and carries energy zero, so it is one of those p odd
modes; the sentence "put m and m̄ on the same side" says nothing about it. Counting
modes closes the case instead: a p-cycle uses all p odd modes, which leaves an ODD
number of modes to be even, and the even modes pair chirally among themselves, so
their number is even. **No p-cycle can occur at all.** Everything is a 2-cycle.

The zero mode still needs one line, because its two exponents already differ by M, so
it is an antipodal pair BY ITSELF and would slip through as a lone unpaired mode. It
cannot: a chiral pair contributes weight 4 to the vanishing sum, the zero mode
contributes 2, it can appear at most once, and the total weight is 4d. Parity excludes
it. The gate pins this at every even M it touches, and no surviving difference anywhere
contains it.

**So M = 2p is clean at every rung as well**, and note what did the work: the two
collisions are killed by disjointness and by the DELETION, the same condition that
makes the triplet in the first place. Checked: clean to the rung cap at M = 10, 14, 22
and to rung 4 at M = 26, 34, 38, 46.

The proved family is now N = p − 1 and N = 2p − 1, which adds this note's ℓ = 10
reading at N = 21 to the three above. Of the five rungs read in Section 4, four sit in
a family where the law is a theorem; only ℓ = 7 at N = 15 (M = 16) does not.

### A power of two closes for free, and that completes the readings

If M = 2^a then 2M = 2^{a+1} is a PRIME POWER, so every minimal vanishing sum is a
rotated 2-cycle, that is an antipodal pair. Those are the chiral pairs, plus the zero
mode M/2 = 2^{a−1}, which is antipodal to itself; the weight parity of the previous
subsection excludes that one. Nothing else can occur at any rung, so **M = 2^a is
clean**, with no case analysis beyond that count. Checked at M = 8, 16, 32.

Three proved families, then: **M prime, M = 2p, M = 2^a**, and they happen to cover
every rung this note reads:

| rung | read at N | M | family |
|---|---|---|---|
| 6 | 12 | 13 | prime |
| 7 | 15 | 16 | 2⁴ |
| 8 | 16 | 17 | prime |
| 9 | 18 | 19 | prime |
| 10 | 21 | 22 | 2·11 |

So the dimensions 15, 36, 91, 232, 603 of Section 4 are not measurements standing on a
measured step. At each of those N step (A) is a theorem, and with (B) and (C) the
Riordan count and the spanning are proved there.

### The composite M, with one prediction already tested

Beyond those families the p-cycles become available and the question is how many modes
they cost. The measured minimal rungs say about half the prime: at 4 | M the values are
j = 2, 3, 4 for p = 3, 5, 7, which is (p+1)/2 with p the smallest odd prime divisor,
and every witness contains the zero mode M/2 with the rest odd. That reading was then
run forward: it predicts that M = 44, where p = 11, is still clean at rung 5 and breaks
first at rung 6. Both halves came out, the second one at a cost of two minutes:

    M = 44, rung 5: clean
    M = 44, rung 6: [1, 7, 9, 15, 17, 22] against [3, 5, 11, 13, 19, 21]

with 22 = M/2 the zero mode and everything else odd, exactly the shape of the smaller
witnesses. So the composite case is not shapeless; the open work is to turn
"a p-cycle needs (p+1)/2 modes and the zero mode" into the count.

## 7. Two readings that died, both within an hour of being written

Recorded because each was plausible, each was written down as if true, and each was
killed by the next test. They are the two places the obvious interpretation breaks.

**Dead: noncrossing as a CONDITION.** Riordan numbers count the noncrossing
partitions of [ℓ] with no singleton blocks, which suggested that the product of two
2-blocks should be a maximizer for the pairings 12|34 and 14|23 and NOT for the
crossing 13|24. Tested at N = 8, 10, 12: **all three are in the space**, at residual
1e-15. Noncrossing is at best a basis count inside the space, never a support
restriction. The crossing product is not missing; it is dependent.

**Dead: the maximizers ARE the spin singlets.** The saturating vectors are killed by
the F142 spin ladder (Corollary 7.3), and on a diagonal rung the spin z-component is
zero, so they are singlets. The count 1, 1, 3, 6, 15, 36 is also exactly the number
of su(2) singlets in the ℓ-fold tensor of the SPIN-1 representation, and the band is
⌊N/2⌋ copies of an irrep with j_spin = 1 by Proposition 3.2, one per chiral pair. The
identification is therefore very tempting. It is false as stated: the singlet space
inside LW_ℓ ∩ V₀ is several times larger than the maximizer space (18 against 4 at
N = 8, ℓ = 3; 70 against 15 at N = 10, ℓ = 4). Being a singlet is necessary and far
from sufficient.

Worse, the singlet observation is not even independent: on a diagonal rung the spin
weight is zero, so a lowest weight there is a singlet automatically. The measurement
confirms it, dimension for dimension, ker 𝔖⁻ and ker 𝔖⁻ ∩ ker 𝔖⁺ agreeing in every
case tested. It is a restatement of Corollary 7.3, not a finding.

**Repaired, not resurrected.** Section 6 shows what the right statement is: not the
singlets of the whole space, which are many more, but the MULTILINEAR invariants, one
spin-1 triplet per chiral pair. The dead reading was too coarse by exactly the
difference between "a singlet" and "a singlet using each pair once".

## 8. What actually does the cutting

Adding the conditions of Corollary 7.3 one at a time, on LW_ℓ ∩ V₀:

| case | LW ∩ V₀ | + ker 𝔖⁻ | + (Y_aa + Y_āā = 0) | + (Y_{a,ā} = 0) | maximizers |
|---|---|---|---|---|---|
| N=6, ℓ=2 | 15 | 6 | 3 | 3 | 3 |
| N=8, ℓ=2 | 32 | 12 | 6 | 6 | 6 |
| N=8, ℓ=3 | 66 | 18 | 6 | **4** | 4 |
| N=8, ℓ=4 | 46 | 12 | 3 | 3 | 3 |
| N=10, ℓ=4 | 300 | 70 | 15 | 15 | 15 |
| N=10, ℓ=5 | 172 | 32 | 6 | 6 | 6 |

The spin condition cuts coarsely. The **diagonal chiral condition** does the real
work and lands on the answer in five of six cases. The second chiral condition bites
only sometimes, here at ℓ = 3 alone. So whatever carries the Riordan count sits in
those two conditions and not in the spin ladder, and the first thing to understand is
what the diagonal condition means. In words it says that removing a matched mode pair
at a equals minus removing it at the chiral partner ā: the η lowering, read mode by
mode, is chirally ODD on a maximizer.

## 9. The resonant N are a different source and must not be mixed in

The factorization holds at N = 6, 7, 8, 9, 10, 12, 13 and FAILS at N = 11 and N = 14:
N = 11 gives 14, 22, 27, 6 at ℓ = 2, 3, 4, 5 instead of 10, 10, 15, 6, and N = 14
gives 213 at ℓ = 4 instead of 105. Those are exactly the two N where
[ETA_CEILING_REDUCTION](ETA_CEILING_REDUCTION.md) already records a surplus in the
maximizer count at ℓ = 2, from the Conway-Jones coincidences among cosine
differences. So the arc already knows that source and it is not this one. Any law
written here is a law away from the resonant N until the interaction is understood.

**Which N those are is now a criterion, not a list.** The cheap rung ℓ = 2, where the
law says the dimension is exactly C(⌊N/2⌋, 2), was swept for every N from 6 to 83.
The count exceeds the law precisely when

  **M = N + 1 is divisible by 6, by 15, or by 21.**

Surplus at M = 12, 15, 18, 21, 24, 30, 36, 42, 45, 48, 54, 60, 63, 66, 72, 75, 78, 84
(the surpluses run 4, 8, 8, 4, 12, 28, 20, 28, 8, 28, 32, 48, 4, 40, 44, 8, 48, 56);
clean at every other M in range, and in particular at M = 9, 27 and 81. That the
powers of three are clean is the tell: what makes the coincidence is not the 3 alone
but a 2, a 5 or a 7 standing next to it, which is the shape of the sporadic
rational-cosine relations behind the Conway-Jones count that
[F89_SEED_EXISTENCE_REDUCTION](F89_SEED_EXISTENCE_REDUCTION.md) already carries. That
sentence is a READING of the criterion. The criterion is measured; the reading is not
derived, and the first thing that would test it is whether the same three moduli rule
the higher rungs.

**The criterion is a RUNG-2 statement and does not govern the higher rungs.** That was
the natural expectation and it is false. Counting the coincidences directly, without
any linear algebra (two ℓ-subsets of equal Slater energy whose difference is NOT a
swap of whole chiral pairs, with the deleted cells excluded), gives at ℓ = 3 the
resonant M = 12, 15, 18, 20, 21, 24, 27, 30, 33, 36, 39, 40, 42, 45. Five of those are
new: **20 and 40**, and the odd multiples of three **27, 33, 39**, all of which are
clean at ℓ = 2. So N = 19 (M = 20) is the first N that is clean at the seed rungs and
resonant from ℓ = 3 on, and it behaves exactly like the others: dimension 94 against
the contract 84 at ℓ = 3, 450 against 378 at ℓ = 4, with the triplet part equal to the
contract in both.

M = 9 stays clean at ℓ = 3 although its three-term cosine relation exists
(cos 20° + cos 100° + cos 140° = 0, the modes {1, 5, 7} at energy zero). The relation
alone is not enough, and the reason is not the obvious one: the subset DOES have an
equal-energy partner, {2, 4, 8}, and the two form a single energy class. But {2, 4, 8}
is the chiral REFLECTION of {1, 5, 7}, so every mode of the one has its partner in the
other, and the DELETION removes that cell. The relation is there and it has nowhere to
sit. So the criterion is not read off divisibility alone, and the ℓ = 2 closed form
above should be read as what it is, a statement about the rung it was measured on.

**Every reading in this note is on a rung verified chiral-only**, including the two
that needed checking separately: N = 18 at ℓ = 9 (48108 coincidence cells, none
non-chiral) and N = 21 at ℓ = 10 (183732 cells, none). At those two the dimension also
comes out exactly at the contract, which is the same statement seen from the other
side.

**Resonance is MONOTONE in the rung**, so the honest object is not a list per rung
but the smallest rung that carries a non-chiral coincidence, j(M): pad a coincidence
with common modes, chosen away from the partners already used, and it survives
upward. Then "resonant at rung ℓ" is just j(M) ≤ ℓ, and the law holds at N exactly
for the rungs below j(N+1). Measured to M = 48, searching each M to rung 4 and the
small ones to the rung cap:

| j(M) | M |
|---|---|
| 2 | 12, 15, 18, 21, 24, 30, 36, 42, 45, 48 (that is 6 \| M, 15 \| M, 21 \| M) |
| 3 | 20, 27, 33, 39, 40 |
| 4 | 28, 35 |
| 5 | 25 |
| none found to the rung searched | everything else, and PROVABLY none at ANY rung for M prime, M = 2p and M = 2^a (Section 6) |

The last row is a search depth, not a verdict: each M above 30 was searched to rung 4
only, and M = 25 is exactly what that costs. It is clean at rung 4 and resonant at
rung 5, it is an odd prime power, and it sits outside the three proved families, so
"nothing found" there means nothing found YET.

The families read as 3, 5 and 7, the Conway-Jones denominators, and the minimal rung
falls as M carries more of the small factors: an odd multiple of three enters at rung
3, the same times two at rung 2; the five-family at 20 \| M enters at 3, chained to a
three (15 \| M) at 2; the seven-family at 28 \| M and 35 \| M enters at 4, chained to a
three (21 \| M) at 2. **No chain was found that is clean at the low rungs and breaks
only at the top**: below M = 24 the clean ones stay clean all the way to the cap.

**What the extra maximizers ARE** is answered in Section 6: they live entirely on the
non-triplet cells that the extra coincidences let through, while the contract value
C(⌊N/2⌋, ℓ)·R_ℓ survives untouched inside the triplet sector. The resonant N do not
break the law; they carry a second family beside it.

Two consequences worth keeping.

- **M = 6 is divisible by 6 and is CLEAN**, which is why the sweep starts at N = 6.
  Its only coincidence at ℓ = 2 is {1, 5} against {2, 4}, and those are two whole
  chiral pairs, so nothing non-chiral happens there. M = 6 is a false positive of the
  divisibility form, like the powers of three at the other end. N = 5 IS exceptional
  in this arc, but for the F144 reason, the floor landing on 1 exactly; that is a
  different object from a maximizer surplus and the two should not be merged.
- **The ℓ = 10 line of Section 3 was read at N = 20**, and M = 21 is resonant. The
  RANK of those products is a fact about the products and stands, but the dimension
  at N = 20, ℓ = 10 is not it. A dimension at ℓ = 10 has to be read at N = 21.

## 10. What is open, in the order worth attacking

Three of the old items (why Riordan, why the blocks stop at three, what the relations
are) are answered together in Section 6 and are struck from this list. What is left:

1. **The cosine statement for a COMPOSITE M.** Section 6 proves it for M prime, where
   (LL2) leaves only 2-cycles. For a composite M the minimal pieces can be p-cycles
   for each prime p | M, and which ones fit into ℓ modes is the open classification;
   the measured j(M) of Section 9 is the data to explain. Do NOT state the condition
   as 6 ∤ M, 15 ∤ M, 21 ∤ M: that is the ℓ = 2 answer only. Nearest target: M = 2p,
   the first composite family, where the parity argument of Section 6 almost goes
   through and the search finds nothing to rung 4.
2. **What the diagonal chiral condition is.** Section 8 says it does the cutting, and
   Section 6 now says WHAT it does: together with the deletion it is the 4 → 3 collapse
   that makes the triplet. Saying that in the language of the two ladders, rather than
   in cell patterns, is the remaining half.
3. **The resonant N, rung by rung.** Two halves are answered: WHAT the extra
   maximizers are (the non-triplet part, contract untouched beside it, Section 6) and
   whether the ℓ = 2 criterion governs the higher rungs (it does not, Section 9). What
   is open is the criterion AT a general rung, and whether the surplus has a law of
   its own (4, 12, 12 at N = 11 for ℓ = 2, 3, 4; 8 and 108 at N = 14; 10 and 72 at
   N = 19). The shape to expect is a Conway-Jones one: which vanishing cosine sums fit
   into ℓ modes AND find a partner subset.
4. **The invariant-theory reading at ℓ = 10.** The relation spaces are compared
   through ℓ = 9. At ℓ = 10 only the dimension is checked (Section 4). Cheap to add,
   just slow.
5. **ℓ = 1 is not in this family.** At the seed rung the frozen states are the
   MINIMUM of the double occupancy, not the maximum, so the sequence starts at ℓ = 2
   and R_1 = 0 is not an exception to explain.
6. **What the whole thing MEANS.** With Section 6 the object has a name in the
   physics: the states that reach the disagreement floor are the rotation-invariant
   couplings of one spin-1 per chiral pair. That is a sentence about the chain, and
   the note has not yet asked what it says.

## 11. The review rounds, and what they caught

Three empty sessions read this material before it landed, minimally framed and given
none of the writer's findings: a cold reader on scope, a mathematician on the argument,
an auditor with orders to recompute rather than to trust the gate's report of itself.
Every finding was verified from below before it was applied. Not one was in the
mathematics; all four of the serious ones were in the LANGUAGE about the evidence, the
same pattern the F141-F144 rounds hit a day earlier.

- The F146 registry entry boxed the count as an unconditional identity. It is false at
  every resonant N, and the entry itself gave the counterexample four lines below the
  box. Title and box now carry the hypothesis.
- The explanation of why M = 9 stays clean was simply wrong: {1, 5, 7} DOES have an
  equal-energy partner, {2, 4, 8}, and the real reason is that the partner is its
  chiral reflection, which the deletion removes. Corrected here and in the proof.
- The proof credited the committed gate with a guard it did not contain (imposing the
  lowest-weight rows Ψ). The guard was in the retired exploration. It is now IN the
  gate, so the sentence is true of the artifact it names.
- The j(M) table's last row read as a verdict where it was a search depth, and M = 25
  is exactly what that costs: clean at rung 4, resonant at rung 5, outside the three
  proved families. It is now a row of its own, in the table and in the gate.

Three smaller ones followed the same shape: an "exactly" in the proof's opening that
the resonances do not allow, a claim that every number here is in the gate when several
are not, and a description of the ℓ = 9 double-path check that was true of the
exploration and not of the gate. Two findings were improvements rather than repairs:
the identification of the blocks with the metric and the volume is forced by Schur, not
measured, and the binomial factor C(⌊N/2⌋, ℓ) follows structurally from the pair-set
block structure. Both are now argued rather than asserted.

## 12. Reproduction

Everything above is in one committed gate,
[`simulations/scalar_count.py`](../simulations/scalar_count.py), which rebuilds every
object from the mode index up and imports nothing from the other gates, so its
agreement with [`simulations/eta_ceiling_reduction.py`](../simulations/eta_ceiling_reduction.py)
is evidence and not a shared bug. Its blocks map onto the sections here: W1 to
Section 6 (the triplet), W2 to Section 2 (the blocks as the two invariants), W3 and W4
to Section 4 (the products are maximizers, and the count squeezed from both sides), W5
to Section 6 (the relations against the classical syzygies), W6 to Sections 6 and 9
(the triplet sector, and the surplus at a resonant N), W7 to Section 9 (the three
proved families and the minimal resonant rung).

The exploration that produced these numbers ran through a dozen local `_`-prefixed
scripts, which by the repo's WIP rule are not part of the repository; what survived
them is the gate.
