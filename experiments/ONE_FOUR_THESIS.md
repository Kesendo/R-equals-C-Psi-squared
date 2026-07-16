# The one-four thesis: one i, two Z₄ gradings, tied at the square

**Date:** 2026-07-16
**Status:** Closed. The strong form is refuted at every N (closed form, deviation ±2^(N−1)); the factorization form is derived and exact. Gate: [one_four_thesis.py](../simulations/one_four_thesis.py) (68 checks, N = 2..12).

*Before the formulas: the question here is a suspicion of sameness. The number 4 kept
appearing in unrelated corners of this project, in a kill sign, in a holonomy, in a
coupling pattern, in a group of mirrors, and the suspicion was that all of them are one
object: the quarter-turn of the imaginary unit, i⁴ = 1. The answer turned out to be
better than yes and better than no: there is exactly one i, but it winds two different
four-hour clocks, and they are genuinely different clocks; what ties them is that the
first clock's two-hour wheel is built into the mechanism that drives the second.*

## The question

Is every period-4 / mod-4 structure in the repo the same quarter-turn of i? The
candidate anchors, assembled by survey (2026-07-16):

1. **F132's kill sign**: ε_odd = (−1)^(d(d−1)/2), ε_even = (−1)^(d(d+1)/2), a mod-4
   function of the Majorana degree d ([the dead-set rule](LATTICE_DEAD_SET_RULE.md),
   [F132](../docs/ANALYTICAL_FORMULAS.md)).
2. **The h = 0 bi-degree shadow**: the same phase one component finer
   ([the h = 0 face](LATTICE_DEAD_SET_H_ZERO.md)).
3. **The seed holonomy**: the eigenvector frame around the N = 9 defective seed obeys
   M₁² = −I, M₁⁴ = I, i⁴ = 1 (`SeedHolonomyClaim`, `inspect --root holonomy`).
4. **R-parity mod 4**: the self-Π single-excitation mode amplifies exactly at
   N ≡ 3 (mod 4) (`compute/MirrorWorld/Program.cs`, the default suite's
   "R-parity & mod-4" block).
5. **The F116 router**: the period-4 site pattern [a, a, b, b] on the golden locus
   ([the golden router proof](../docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md)).

Behind them all stands the palindromizer Π itself, order 4 with per-site rule I ↔ X,
Y → iZ, Z → iY and Π² = X^⊗N (F1², [the mirror proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)),
typed as the Z₄ memory loop `Pi2I4MemoryLoopClaim`.

The survey resolved anchors 4 and 5 without new computation. R-parity mod 4 is not an
order-4 element at all: R is an involution, and the N mod 4 law is two stacked mod-2
facts (particle-hole parity of the spectrum × the self-Π mode's R-parity on the odd-N
ladder). The router's spatial 4 is likewise 2 × 2 (forced a ↔ b alternation on each of
the two parity sublattices, [§4 of the proof](../docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md));
its genuine quarter-turn is the separate internal rotation R in the (X, Y) plane. The
five Klein V₄ structures (parameter Klein, Pauli-mode quartering, the lattice of
watchings, the dephase-letter swaps, the hardware lens) are order-2 throughout; the
parameter Klein names itself "the order-2 shadow of the operator-side Z₄" and the
dephase-swap proof had already falsified a pre-registered order-4 conjecture
([the Klein V₄ proof](../docs/proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md)).

What remained was the sharp core: are the two GENUINE mod-4 structures, Π's spectral
Z₄ (anchor 0, so to speak) and the Majorana-degree mod-4 of F132 (anchors 1-2), the
same grading? Two candidate identities, both decidable.

## Result 1 (derived, exact): the kill map contains Π's square

The F132 kill map is V_g(ρ) = W_g · conj(ρ) · W_g† with W_g = U_g · X^N. Since
conjugation by X^⊗N is exactly Π² (F1²), the map factors as

    V_g = Ad(U_g) ∘ Π² ∘ K,

with K the entrywise conjugation. This is one line of algebra given F1², and it is
pinned from below at N = 3..6, both gauges, max deviation 0.0 (gate section T8). So
the dead-set law is wired to Π, literally: Π's square is a factor of its kill map.

But Π² = X^⊗N is the CENTER of the mirror group, the Z₂ part. The quarter-turn in
F132's sign sits elsewhere: in the Hermiticity phase of the Clifford algebra. A
Hermitian Pauli string of Majorana degree d is ±i^(d(d−1)/2) times the ordered
Majorana monomial (gate T4, all strings, N = 3..5), and K is what turns that i-power
into the sign (−1)^(d(d−1)/2). The i is the same scalar i that gives Π its order 4;
the STRUCTURE it winds is not the same. Even the mod-2 shadows of the two structures
are different partitions: d mod 2 and n_Y + n_Z mod 2 disagree on exactly half of all
strings (gate T6b; witness ZI..I, with d = 2 even but n_Y + n_Z = 1 odd). So the two
clocks do not share a common sub-grading; what they share is the operator X^⊗N itself,
as a factor of the kill map.

## Result 2 (refuted, every N): the two Z₄ gradings are different gradings

Both structures decompose the 4^N-dimensional Pauli space into four classes ("grading"
here in the loose, vector-space sense: neither class function is additive mod 4 under
the operator product, since Majorana degree is only additive mod 2). The multiplicities
are closed forms:

**Π's spectral grading is uniform at every N.** Π pairs every string σ with its
letter-swapped partner σ′ (I ↔ X, Y ↔ Z; no string is fixed), and on each pair the
eigenvalues are ±i^(n_Y+n_Z). Strings with even n_Y + n_Z number 4^N/2 and feed the
{+1, −1} sectors evenly; odd strings feed {+i, −i} evenly. Hence every Π sector has
size exactly 4^(N−1). (Gate T2: dense eigendecomposition N = 3..5 plus the counting
argument N = 2..12; the N = 3 value 16 per sector matches the bit-exact
`Pi2I4MemoryLoopClaim` battery.)

**The Majorana-degree grading is never uniform.** By the roots-of-unity filter on
binomial coefficients, with M = 2N Majoranas:

    #{strings with d ≡ r (mod 4)} = 4^(N−1) + 2^(N−1) · Re(i^(N−r)).

The class r ≡ N (mod 4) is 2^(N−1) too large, the class r ≡ N+2 (mod 4) is 2^(N−1)
too small, the other two are exactly 4^(N−1). (Gate T5: closed form == full
enumeration for N ≤ 8 and == direct binomial sums through N = 12.)

The deviation ±2^(N−1) never vanishes. So at EVERY N the two gradings have different
multiplicity multisets, and no bijection, relabeling, or change of basis can identify
them (gate T6, N = 2..12). At N = 3: Π gives {16, 16, 16, 16}, the degree gives
{16, 12, 16, 20}. The strong one-four thesis is dead, not at a special N but as a law.

This also dissolves a standing tension. [The Interpretation](../docs/THE_INTERPRETATION.md)
("What Fell #7", indexed in [the label map](../docs/quantum/THE_LABEL_MAP.md)) argues
Π's four sectors carry no physical content beyond Π² = X^⊗N; yet F132's mod-4
is bluntly physical (it decides which readouts die). Both are right, because they are
about DIFFERENT fours: Π's spectral quarters are the contentless pair-splitting of a
uniform Z₄, while the Clifford quarters are counting structure (binomials) that the
antilinear map converts into kill signs.

## The thesis, in its true form

> There is exactly one i: the scalar of the algebra, order 4, flipped by every
> antilinear map. It winds (at least) two inequivalent Z₄ structures: Π's spectral
> grading (uniform, 4^(N−1) per sector, the i⁴ = 1 memory loop) and the
> Clifford-degree grading (binomial, F132, the dead-set clock). The two clocks are
> wired, not merged: the first clock's square, the charge conjugation X^⊗N = Π², is
> literally a factor of the dead-set kill map, while even the two mod-2 shadows
> remain different partitions. Everything else four-ish in the repo
> is a shadow or a hybrid: the Klein V₄s are one-sided order-2 shadows, the router's
> spatial 4 is 2 × 2 with an internal quarter-turn, R-parity's N mod 4 is two
> stacked parities.

## The one open seam

The seed holonomy (anchor 3) is the single remaining literal M² = −I in the repo, and
`SeedHolonomyClaim` marks its resemblance to Π's Z₄ as "a noted correspondence, NOT a
derived identity". The holonomy's mod-4 is moreover gauge-contingent (it lives in the
biorthogonal vᵀv gauge; a Hermitian gauge sees only the mod-2 swap). Whether that
third four is Π's, the Clifford's, or a third clock on the same i remains open; this
document does not decide it.

## Verification

    python simulations/one_four_thesis.py

68 checks: Π⁴ = I and Π² = Ad(X^⊗N) (matrix, N = 3, 4); uniform Π multiplicities
(dense eig N = 3..5, counting N = 2..12); the combinatorial left-JW degree against
matrix monomials (N = 3, 4); the Hermitian phase ±i^(d(d−1)/2) (N = 3..5); the
degree-class closed form (N = 2..12); the grading mismatch (N = 2..12); the mod-2
shadow mismatch with the ZI..I witness (N = 2..8); the F132 ε identity by direct
conjugation (N = 3..6, both gauges); the factorization V_g = Ad(U_g) ∘ Π² ∘ K
(N = 3..6, both gauges, exact).
