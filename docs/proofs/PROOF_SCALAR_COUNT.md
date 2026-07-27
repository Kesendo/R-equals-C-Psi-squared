# The scalar count: what reaches the disagreement floor, and how much of it there is

**Status:** F145 proved at every N; F146 proved at every rung for M = N + 1 prime,
M = 2p and M = 2^a, and measured elsewhere. Gate:
[`simulations/scalar_count.py`](../../simulations/scalar_count.py), 67 checks in
about a minute; `--deep` reaches 78, adding the rungs ℓ = 9 and ℓ = 10, the M = 44
prediction and a wider range of N in the checks W2, W3, W5 and W6,
and takes well over an hour.

Before the machinery, what this feels like. [F144](../ANALYTICAL_FORMULAS.md) says the
watching charges every state on a rung at least ℓ(N − ℓ)/(N + 1), and leaves open who
pays exactly that. The answer turns out not to be a list of vectors but an old object
in a new place: **each chiral pair of standing waves carries a spin 1, and the states
that pay the floor are the rotation-invariant couplings of those spins.** That last
sentence has one condition on N, stated as step (A) in Section 4 and proved there for
three families; on the cosine resonances of Section 7 those couplings are still all
there and a second family stands beside them.

So the count is a singlet count, the building blocks are the two invariants any three
dimensions have, a dot product and a volume, and the relations between them are the
classical ones. Nothing here is new mathematics; what is new is that this is where the
chain keeps its longest-lived memories.

---

## Symbols

| symbol | meaning |
|---|---|
| N, M | chain length and the transform length M := N + 1 |
| h | the single-excitation hopping matrix, eigenmodes v_k(l) ∝ sin(πkl/M), energies ε_k = 2J cos(πk/M) |
| k̄ | the chiral partner M − k of a mode, with ε_k̄ = −ε_k |
| cell (A, B) | \|A⟩⟨B\| with A, B sets of MODES; on a diagonal rung \|A\| = \|B\| = ℓ |
| V₀ | ker(ad_h), spanned by the cells whose two mode sets carry equal Slater energy |
| Φ, Ψ | the eta ladder and its adjoint; LW_ℓ := ker Ψ at rung ℓ |
| 𝔖⁺, 𝔖⁻ | the F142 SPIN ladder, (p,q) → (p±1, q∓1) |
| 𝒦, D̂ | the disagreement count and the double occupancy, 𝒦 = ℓ − D̂ on a diagonal rung |
| X(a,b) | the map (A,B) ↦ (A\a, B\b) with the fermionic signs, zero unless a ∈ A and b ∈ B |
| R_ℓ | the Riordan number (A005043): 1, 0, 1, 1, 3, 6, 15, 36, 91, 232, 603 for ℓ = 0..10 |
| p | in "M = 2p" and "j(M)" an odd prime; elsewhere ⌊N/2⌋, the number of chiral pairs. The two never appear in one formula |
| block | bare, the repo's joint-popcount block (p, q), as in "the blocks (2,0), (1,1), (0,2)". With a prefix, **2-block** and **3-block** are the connected maximizers of this document, which Section 3 identifies with the metric and the volume. The two senses are unrelated and both are needed here |

**The maximizer space** M(N, ℓ) is the set of v ∈ LW_ℓ ∩ V₀ attaining the F144 floor,
equivalently the top eigenspace of D̂ there. By F144 (Corollary 7.3 of
[PROOF_FROZEN_BAND_SO4](PROOF_FROZEN_BAND_SO4.md)) it is exactly the set of v ∈ V₀ with

    (i)   X(a, ā) v = 0 for every mode a,
    (ii)  (X(a,a) + X(ā,ā)) v = 0 for every mode a,
    (iii) 𝔖⁻ v = 0,

whenever that set is nonzero. The lowest-weight condition needs no separate
imposition: Ψ = Σ_a X(a,a), and since a ↦ ā is a bijection of the modes,
Ψ = ½ Σ_a [X(a,a) + X(ā,ā)], which (ii) annihilates. Gate W4 measures it as well, by
imposing the Ψ rows and finding the nullity unchanged.

---

## 1. The collapse: two conditions, and only three patterns survive per pair

Every one of the conditions above is built from X(a,b), and X(a,b) is INJECTIVE on the
cells it does not annihilate: distinct (A,B) with a ∈ A, b ∈ B go to distinct
(A\a, B\b). So every ROW of X(a,b) carries AT MOST one entry, and the two
chiral conditions are not linear algebra at all:

- **(i) deletes coordinates.** X(a,ā)v = 0 says v vanishes on every cell with a in A
  and ā in B.
- **(ii) pairs the survivors.** Its row at (A₂, B₂) has exactly the two entries
  (A₂ ∪ a, B₂ ∪ a) and (A₂ ∪ ā, B₂ ∪ ā), so it glues two coordinates with a sign, and
  a sign clash around a cycle sets a component to zero.

For one chiral pair P = {a, ā} the cells that survive (i) contribute one of

    u₊(P) = both modes of P in A,        u₋(P) = both modes in B,
    u₀(P) = (a in A and in B) − (ā in A and in B),        or nothing at all,

the last line being what (ii) leaves of the two mixed patterns. **Four patterns become
three.** The mixed cell a|ā is deleted by (i); a|a and ā|ā are glued by (ii).

## 2. F145: the three patterns are a spin-1 multiplet of the F142 ladder

The three patterns live in the blocks (2,0), (1,1) and (0,2), which is exactly where
𝔖± moves, and the ladder acts on them as one spin-1:

    𝔖⁻ u₊ = ∓u₀,    𝔖⁻ u₀ = ∓2u₋,    𝔖⁻ u₋ = 0,

with 𝔖⁺ the mirror image and the per-pair sign immaterial. The coefficients are 1 and
2 rather than √2 and √2 because u₀ is normalized as a difference of two cells and u±
as single ones; their PRODUCT is 2, which is the basis-free spin-1 invariant, and the
same normalization is what puts the weight 2 on the ± terms of the metric in Section 3.

*Proof.* Direct: 𝔖⁻ moves a mode from the bra to the ket at the chiral partner index.
On u₊ the two available moves take the modes a and ā to a and ā on the other side,
which is the u₀ combination with the sign the fermionic reordering gives; on u₀ the
single available move on each term produces u₋ twice, hence the 2; on u₋ there is no
mode left in the bra. □ (Gate W1, at N = 6, 8, 9, 10, 12, 15, every chiral pair.)

**What u₀ is.** It is the F143 **seed** P_k − P_k̄ written in the cell basis: the
difference of the spectral projectors at chiral partners. In the site basis that
object has zero diagonal everywhere, lives only on site pairs at odd distance, and
commutes with h at every coupling. So the middle of the triplet is the frozen mode of
the seed rung, and u± are the coherences that create and annihilate the same mirror
pair out of the vacuum. **Seed, birth, death.**

## 3. F146, first half: the two connected blocks ARE the two invariants of SO(3)

Two spaces need separate names, because they sit one step apart. The **triplet-cell
space** on a set S of chiral pairs is the span of the cells in which every pair
contributes one of the patterns of Section 1 with a|a and ā|ā still separate;
condition (ii) still cuts inside it, and it is what step (A) and gate W6 talk about.
What (ii) leaves is the **sector** T(S), the span of the products of one of u₊, u₀,
u₋ per pair. On T(S) the ladder acts as a sum of one-pair operators, hence as the
tensor product action, so

    T(S) ≅ V₁^{⊗|S|}

as a representation of the F142 SU(2), with V₁ the spin-1.

The two connected blocks of the maximizer family are then the two invariants that
three dimensions have, exactly and cell by cell:

    2-block(P, Q)    = u₀(P)u₀(Q) − 2[u₊(P)u₋(Q) + u₋(P)u₊(Q)]        (the METRIC)
    3-block(P, Q, R) = −ε(P, Q, R)                                     (the VOLUME)

where ε is the antisymmetric form supported where the three m-values are +1, 0, −1 in
some order, with the permutation sign. The 2-block is the six-term closed form
[ETA_CEILING_REDUCTION](../../experiments/ETA_CEILING_REDUCTION.md) already carried;
the 3-block is twelve terms with entries ±1, one pair full in the bra, another full in
the ket, a shared mode of the third, and it is transpose-ANTIsymmetric where the
2-block is symmetric. (Gate W2: 83 blocks in the default run, all 258 at
N = 6..16 under `--deep`.)

**Neither identification rests on the finite range of that gate.** On two pairs the
space of ladder-annihilated vectors is one-dimensional, and so is Inv(V₁ ⊗ V₁); the
same holds on three pairs, R₂ = R₃ = 1 being exactly that statement. By Schur the
block therefore IS the metric, respectively the volume, up to a scalar at every N.
What W2 adds is the scalar, and that it is 1 and −1 in the normalization above.

In this normalization the classical syzygy reads ε·ε = −¼·det(Gram), so the 3-block is
(i/2) times the ordinary triple product, which is the one constant the comparison in
Section 6 needs.

## 4. F146, second half: the count

**Theorem.** Let N be such that step (A) below holds at rung ℓ. Then

    dim M(N, ℓ) = C(⌊N/2⌋, ℓ) · R_ℓ,

and M(N, ℓ) is spanned by the products of 2- and 3-blocks over the partitions of the
chosen ℓ chiral pairs into blocks of size 2 and 3.

The proof is three steps, of which only the first depends on N.

### (A) Every maximizer lies in the triplet sector

Equivalently: if two ℓ-subsets A, B of modes carry the same Slater energy and no mode
of A has its chiral partner in B, then A \ B and B \ A are unions of whole chiral
pairs, **and neither contains the zero mode** M/2 when M is even. Call such a rung
**chiral-only**. The zero mode needs saying separately because it is its own chiral
partner, so a lone copy of it would satisfy "union of chiral pairs" on a reading the
count C(⌊N/2⌋, ℓ) does not use: that count is over the ⌊N/2⌋ genuine PAIRS, and the
zero mode is not one of them.

**The weight bookkeeping, once, for all three cases.** Write d := |A \ B| = |B \ A|.
Each mode contributes two exponents, so the vanishing sum has weight 4d. A chiral pair
sitting inside one difference contributes exactly two antipodal 2-cycles, weight 4.
The zero mode, if it occurs, contributes ONE antipodal 2-cycle by itself, since its two
exponents M/2 and 3M/2 already differ by M: weight 2. And it can occur at most once,
the two differences being disjoint. So if every minimal piece is a 2-cycle, the weight
reads 4d = 4·(number of chiral pairs) + 2·(0 or 1), which forces the second term to be
0. **The zero mode is excluded by parity**, and what is left is chiral pairs.

Writing ζ = e^{iπ/M}, a primitive 2M-th root of unity, and moving the minus signs
through ζ^M = −1, the energy equation becomes a vanishing sum of 2M-th roots of unity
with nonnegative coefficients,

    Σ_{x ∈ A\B} (ζ^x + ζ^{−x}) + Σ_{y ∈ B\A} (ζ^{y+M} + ζ^{M−y}) = 0,

of weight 4·|A \ B|, whose exponents **avoid 0 and M** and whose terms are **pairwise
distinct**, the latter precisely because a mode of A whose partner lies in B is
forbidden. This is the same reduction that
[PROOF_F129_LEVEL_COLLISION_LAW](PROOF_F129_LEVEL_COLLISION_LAW.md) §2 performs for
level collisions, and it imports the same two facts (Lam-Leung, *On vanishing sums of
roots of unity*, J. Algebra 224 (2000)): every vanishing sum with nonnegative
coefficients splits into minimal ones, and a minimal vanishing sum of m-th roots with
m = p^a q^b is a rotated p-cycle or q-cycle.

**Case M prime.** Then 2M = 2p has two prime divisors, so every minimal piece is a
rotated 2-cycle or a rotated p-cycle. A rotated p-cycle steps by 2, so its exponents
are ALL residues of one parity class mod 2p; the even class contains 0 and the odd
class contains p, and the sum has neither, so no p-cycle occurs. A rotated 2-cycle is
an antipodal pair {ζ^u, ζ^{u+M}}; the only antipodal pairs available are a mode
together with its chiral partner inside A \ B or inside B \ A, since the cross cases
are excluded by disjointness and by the deletion. Hence both differences are unions of
chiral pairs. **Every rung is chiral-only.** □

**Case M = 2p.** Now 2M = 4p, still two prime divisors, and a p-cycle is a full
residue class mod 4. Classes 0 and 2 contain the exponents 0 and M = 2p, which are
absent, so a p-cycle would have to fill class 1 or class 3, each of which has p
members. Only the ODD modes contribute to those classes, there are exactly p of them
(1, 3, …, 2p−1), and each one used contributes exactly one exponent to each odd class.
So a p-cycle would force every odd mode to be used exactly once, leaving 2d − p modes
to be even. But the even modes are 2, 4, …, 2p−2, they pair chirally among themselves
with no self-paired member, so they occur in even number; and 2d − p is odd, p being
odd. **No p-cycle can occur.** Every piece is then a 2-cycle, and the bookkeeping above
gives chiral pairs and no zero mode. □

**Case M = 2^a.** Then 2M is a prime power, so every minimal piece is a rotated
2-cycle, that is an antipodal pair. Here the zero mode M/2 = 2^{a−1} exists and is its
own antipodal pair, so "antipodal pair = chiral pair" is not quite the whole story;
the weight bookkeeping above supplies the rest and excludes it. **Every rung is
chiral-only**, still with no case analysis beyond the parity count. □

Outside these three families (A) is measured rather than proved; Section 7 says what
is measured and where it fails.

### (B) Inside the triplet sector the conditions say exactly "singlet"

By Section 1 the conditions (i) and (ii) carry the triplet-cell space down to the
sector T(S): (i) deletes, (ii) glues a|a to ā|ā into u₀. So on T(S) nothing of (i) or
(ii) is left to impose, and what remains is (iii), 𝔖⁻v = 0.

By Section 3 that sector is V₁^{⊗ℓ} as a representation of the spin SU(2). Its weight
is (|A| − |B|)/2, which is 0 on a diagonal rung. A vector killed by the lowering
operator is a lowest-weight vector of some irreducible summand V_j and has weight −j;
weight 0 forces j = 0. So the vectors satisfying (iii) inside the sector are exactly
the **singlets** of V₁^{⊗ℓ}, and conversely every singlet satisfies it. □

### (C) The singlets are the multilinear SO(3) invariants

so(3) ≅ su(2) with V₁ the vector representation, so the singlets of V₁^{⊗ℓ} are the
multilinear SO(3)-invariant functions of ℓ vectors in three dimensions. Two classical
facts finish the count.

- **First fundamental theorem for SO(3).** The invariants of vectors are generated by
  the inner products and the triple scalar products. Multilinearly in ℓ vectors that
  means: spanned by the products over partitions of {1..ℓ} into blocks of size 2 and 3,
  which by Section 3 are exactly the products of 2- and 3-blocks. **This is why the
  connected blocks stop at three: there is nothing of degree four to find.**
- **The dimension.** dim Inv(V₁^{⊗ℓ}) = R_ℓ, the Riordan number.

The factor C(⌊N/2⌋, ℓ) is the choice of which ℓ chiral pairs are touched, and the
choices do not mix. Count modes: if a triplet cell has p₊, p₀, p₋ pairs in the
patterns u₊, u₀, u₋ then |A| = 2p₊ + p₀ and |B| = 2p₋ + p₀, so p₊ = p₋ and the number
of pairs touched is p₊ + p₀ + p₋ = ℓ exactly. Every operator in the conditions keeps
that set: 𝔖⁻ moves a mode from the bra to the ket WITHIN its own pair, and the
X-conditions act inside one pair as well. So the system is block diagonal over the
C(⌊N/2⌋, ℓ) pair-sets, each block contributing R_ℓ on disjoint support, and the
binomial is structural rather than measured; gate W4 corroborates it by building the
products over EVERY choice rather than one choice times a binomial. □

**What the relations are.** The second fundamental theorem: the syzygies among dot and
triple products, whose multilinear part begins with the vanishing 4 × 4 Gram
determinant. That is why the first relation among the pure-2-block products appears at
ℓ = 8 and not at ℓ = 4: four vectors in three dimensions need four more slots to state
their dependence.

## 5. What the theorem gives, rung by rung

| ℓ | read at N | M | family | products | dim |
|---|---|---|---|---|---|
| 6 | 12 | 13 | prime | 25 | 15 |
| 7 | 15 | 16 | 2⁴ | 105 | 36 |
| 8 | 16 | 17 | prime | 385 | 91 |
| 9 | 18 | 19 | prime | 1540 | 232 |
| 10 | 21 | 22 | 2·11 | 7245 | 603 |

Each of those N lies in one of the three proved families, so these are theorems and
not measurements; the gate squeezes each of them anyway between a product rank from
below and a condition nullity from above, both mod p and one-sided in directions that
compose, the rungs ℓ = 9 and ℓ = 10 under `--deep`. At ℓ = 8 the two competitors part: the perfect matchings would give 105 and
the true dimension is 91.

## 6. The measured side

Three things in this document are measured and not derived, and they are stated here
rather than left to be inferred.

1. **Step (A) outside the three families.** See Section 7.
2. **The constant in the relations comparison.** The classical relations transfer to
   our products once each product carrying t three-blocks is weighted by (−4)^{t/2},
   which is (i/2)^{−t}. The constant was fitted at ℓ = 6, where all ten relations
   return −4, and then held fixed as a prediction at ℓ = 7, 8, 9; the relation
   subspaces coincide, not merely their dimensions. Gate W5 carries this to ℓ = 7 by
   default and ℓ = 8 under `--deep`; the ℓ = 9 comparison is from the exploration.
3. **The spanning at ℓ = 10.** Verified as a dimension, with the relation subspace
   compared only through ℓ = 9.

## 7. Where (A) fails, and what happens there

A non-chiral coincidence at rung j survives to every higher rung (pad both sides with
common modes chosen away from the partners already used), so resonance is MONOTONE and
the object is j(M), the smallest rung carrying one. Measured to M = 48:

| j(M) | M |
|---|---|
| 2 | 12, 15, 18, 21, 24, 30, 36, 42, 45, 48, that is 6 \| M or 15 \| M or 21 \| M |
| 3 | 20, 27, 33, 39, 40 |
| 4 | 28, 35 |
| 5 | 25 |
| none found to the rung searched | everything else, and provably none at any rung for the three families of Section 4 |

The last row records a search depth. Each M above 30 was searched to rung 4, and
M = 25 shows what that costs: clean at rung 4, resonant at rung 5, an odd prime power
outside the three proved families.

The families read as the Conway-Jones denominators 3, 5 and 7, and the minimal rung
falls as M carries more of the small factors. At 4 \| M the measured values are
j = 2, 3, 4 for the smallest odd prime divisor p = 3, 5, 7, which is (p+1)/2, and every
witness contains the zero mode M/2 with the rest odd. Run forward, that predicts
M = 44 clean at rung 5 and breaking first at rung 6; both halves hold, with witness
[1, 7, 9, 15, 17, 22] against [3, 5, 11, 13, 19, 21] (gate W7 under `--deep`).

**At a resonant N the law is not violated, it is accompanied.** The triplet part of the
maximizer space still has exactly the dimension C(⌊N/2⌋, ℓ)·R_ℓ, and the whole surplus
sits on the non-triplet cells that the extra coincidences let through: 14 against 10 at
N = 11, ℓ = 2; 213 against 105 at N = 14, ℓ = 4; 94 against 84 at N = 19, ℓ = 3 (gate
W6). So the extra maximizers are a second family beside the scalars, not more of them.

## 8. What is not proved

- **(A) for a general composite M.** The classification of which vanishing cosine sums
  fit into ℓ modes AND find a partner subset. M = 9 shows why divisibility alone is not
  the criterion: its three-term relation exists (cos 20° + cos 100° + cos 140° = 0,
  that is the modes {1, 5, 7} at energy zero), and M = 9 is still clean at that rung.
  The reason is not that the relation has no partner, it does: {2, 4, 8} carries the
  same energy and the two form a single energy class. It is that {2, 4, 8} is the
  chiral REFLECTION of {1, 5, 7}, so every mode of the one has its partner in the
  other and the DELETION removes that cell. The relation is there and it has nowhere
  to sit.
- **A law for the surplus** at the resonant N. The numbers are 4, 12, 12 at N = 11 for
  ℓ = 2, 3, 4; 8 and 108 at N = 14; 10 and 72 at N = 19.
- **Anything about ℓ = 1.** At the seed rung the frozen states are the MINIMUM of the
  double occupancy rather than the maximum, so the family starts at ℓ = 2 and R_1 = 0
  is not an exception to explain.

---

**Measured in** [WHAT_REACHES_THE_CEILING](../../experiments/WHAT_REACHES_THE_CEILING.md),
which carries the measurements in the order they were made, including the two readings
that died on the way.
**Gate:** [`simulations/scalar_count.py`](../../simulations/scalar_count.py).
**Parent:** [PROOF_FROZEN_BAND_SO4](PROOF_FROZEN_BAND_SO4.md) Section 7 (F144), whose
open attainment question this answers.
