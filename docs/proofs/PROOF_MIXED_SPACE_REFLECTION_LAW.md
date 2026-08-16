# PROOF: the mixed-space reflection law, or why the parity halves cannot hear each other

*2026-08-15. Closes the first open item of
[THE_ENDPOINTS_ARE_A_DENSITY_LAW](../../experiments/THE_ENDPOINTS_ARE_A_DENSITY_LAW.md)
(the arc `compressed_density_laws`, NextStep 1): the 91 parity-mixed ad_H
eigenspaces at N = 8, Δ = 0 satisfy C_l = 0 although the symmetry argument leaves
their cross block free. That was a LAW, measured and underived. It is now a
theorem, and a stronger one than was measured: not only the odd combination
N_l − N_{N−1−l} but EVERY compressed site density Π_Ω N_l Π_Ω separately is
parity-block-diagonal. The tool the arc named, F129's collision machinery, is the
right one, one weight level down: what decides is not the triple level map but
the PAIR-sum map of the free-fermion comb, and at M = 9 that map is injective up
to the chiral zeros. Gate:
[`simulations/mixed_space_reflection_gate.py`](../../simulations/mixed_space_reflection_gate.py),
VERDICT green 2026-08-15.*

## What this is about

An open chain whose couplings read the same from either end has a site mirror
down its middle: run the sites left to right or right to left and the physics
is the same. The chain's basic patterns, coherences connecting two spin
arrangements as the two faces of one object
([Superposition Translated](../quantum/SUPERPOSITION_TRANSLATED.md)),
therefore come in two families, mirror-even and mirror-odd.
Patterns that oscillate at exactly the same frequency share one room, and
inside such a room the symmetry alone does not keep a site-by-site question
from coupling an even pattern to an odd one; the question being whether the
pattern's two faces disagree at this particular site. Measurement had
found the lopsided combination of these questions, site minus mirror site,
silent in every such room. This proof explains why, and proves more than was
measured: every site's question separately respects the divide. The setting
is the free-hopping chain with equal couplings on every bond, where the
question reduces to arithmetic among the chain's standing-wave energies: coupling the two families would need two pairs of those energies to
add up to the same value, and the only coincidences of that kind the ladder of
cosines allows turn out, on inspection of their signs, to connect even to even
and odd to odd, never across. So each site's question sounds entirely within
one family, and the mirror halves cannot hear each other. The protection is
proved outright for one chain length, forced by the same argument for an
infinite family of others, and checked for several more; the proof also maps
where it ends, the lengths at which new coincidences become
available, the doors through which mixing could enter. Away from the
free-hopping setting the argument gives nothing, and at the fully symmetric
setting the law is measured to break; the document says so.

## What the repo already holds

The sweep, by store. `docs/ANALYTICAL_FORMULAS.md`: [F2b](../ANALYTICAL_FORMULAS.md#f2b)
owns the OBC dispersion and sine modes this proof stands on; F129/F130 own the
six-cosine collision machinery and its Lam-Leung imports; F146 owns the
chiral-only rung notion and names its failure "the cosine resonances of F129's
kind"; F154 is the entry this proof upgrades: its mixed-space clause was
registered 2026-08-15 as measured, not derived, and now reads DERIVED with a
pointer here (the co-landing edit).
`docs/proofs/`: [PROOF_F129_LEVEL_COLLISION_LAW](PROOF_F129_LEVEL_COLLISION_LAW.md)
§2-§3 owns the reduction of equal cosine sums to vanishing root sums and the
(LL1)-(LL3) source pins, its (LL2) is exactly the two-prime decomposition
fact the pair lemma below needs at order 18, and its §4 imports the
Poonen-Rubinstein classification of minimal vanishing sums of weight ≤ 12,
which already closes the repeated-exponent subtlety the lemma meets (height 1,
all coefficients +1) and would settle the weight-8 case directly;
[PROOF_SCALAR_COUNT](PROOF_SCALAR_COUNT.md) §4 proves every rung chiral-only
for ℓ-subsets at M prime, M = 2p, M = 2^a under its deletion hypothesis,
measures the rung-2 doors 6|M, 15|M, 21|M (§7; M = 6 itself is F146's
recorded exception), and its §8 is the nearest
prior statement about this very gap ("M = 9 is still clean at that rung");
none of that covers the unrestricted multiset statement at M = 9 = 3² that
the pair lemma below proves in the same frame;
[PROOF_C1_MIRROR_SYMMETRY](PROOF_C1_MIRROR_SYMMETRY.md) (Ingredient 2) owns
R ψ_k = (−1)^(k+1) ψ_k (measured and derived in
`review/OBC_SINE_BASIS_FINDINGS.md`).
`experiments/`: the endpoint note holds the measured law, the census (91 mixed
spaces: 20 of 59 on (1,3), 21 of 83 on (1,4), 50 of 143 on (2,4)) and the open
item; `experiments/SLOW_MODE_R_PARITY.md` holds the operator-space product rule
"R-parity of σ_(a,b) is the product ε_a·ε_b" and the diagonal action of R on
Bogoliubov modes; `experiments/F89_SEED_EXISTENCE_REDUCTION.md` Lemma 3 holds the
mode-triple reflection sign, the p = 3 case of the product rule proved below in
general. `fw.Confirmations`: nothing (no hardware content). `docs/GLOSSARY.md`:
nothing for the mixed/scalar split. `docs/CAUGHT_ERRORS.md`: the standing shapes
"generalizes verbatim is worth a gate every time" and the F101 scope-qualifier
catch, both applied here (the general-M section states exactly which M each
argument covers). The OpenArcs registry: `compressed_density_laws` NextStep 1 is
this question, verbatim. No store derives the law.

## §1 Setting and statement

The uniform open XXZ chain in the Pauli book, H = J·Σ_bonds(XX + YY + Δ·ZZ), at
**Δ = 0** and **uniform J**. M = N + 1 is the comb size; the single-particle
modes are F2b's sine modes ψ_k(z) = √(2/M)·sin(kπ(z+1)/M), k = 1..M−1, with
energies ε_k = 4J·e_k in this book, e_k = cos(kπ/M) the normalized comb (F2b's
book carries J/2 and reads 2J·e_k; every collision predicate below is
factor-blind, a factor-2 book remark of the same kind as PROOF_F130 §1's). A coherence block
(p, q) holds the cells |a⟩⟨b| with popcount(a) = p, popcount(b) = q; ad_H acts
on it; R is the site reflection (a, b) ↦ (rev a, rev b); N_l the diagonal
indicator that a cell disagrees at site l; Π_Ω, comp, scalar/mixed as in the
endpoint note's vocabulary block.

**Theorem (the mixed-space reflection law, derived).** At Δ = 0, uniform J,
N = 8: on every ad_H eigenspace Ω of every coherence block, and for every site
l separately,

    Π_Ω N_l Π_Ω  has zero matrix elements between the R-parity halves of Ω.

Consequently every real combination Σ_l c_l Π_Ω N_l Π_Ω is parity-block-diagonal,
and C_l = Π_Ω (N_l − N_{N−1-l}) Π_Ω vanishes on every mixed Ω: the cross block by
this theorem, the two diagonal blocks by the R-odd conjugation argument the
endpoint note already proved. The N = 8, Δ = 0 clause of
[F154](../ANALYTICAL_FORMULAS.md#f154) is thereby DERIVED, and the saturation on
the size-class centres at that (N, Δ) no longer rests on a measurement.

The proof has four steps: the Slater pair eigenbasis (§2), the parity product
rule (§3), the one-body selection rule and the collision equation it leaves
(§4), and the pair-sum lemma at M = 9 with the parity bookkeeping that closes it
(§5). §6 verifies, §7 states what the same argument gives at other N and where
it stops.

## §2 The Slater pair eigenbasis

At Δ = 0 the chain is free-fermionic (Jordan-Wigner; the route the repo already
owns: `docs/SYMMETRY_FAMILY_INVENTORY.md`, `JwSlaterPairBasis` in
Core/BlockSpectrum/JordanWigner, PROOF_SCALAR_COUNT's Symbols block). The p-particle sector
is diagonalized by the Slater states |A⟩ = b†_{k₁}···b†_{k_p}|0⟩ over p-subsets
A ⊂ {1..M−1}, E_A = Σ_{k∈A} ε_k, and the block operator
ad_H = H^(p) ⊗ I − I ⊗ H^(q) (H real symmetric, so the bra-side transpose is
itself) is diagonalized by the pair states |A⟩⟨B|, eigenvalue determined by
E_A − E_B. Every eigenspace Ω is a span of pair states, and distinct pair states
are orthonormal, so the parity halves of Ω can be read off pair state by pair
state.

## §3 The parity product rule

**Lemma.** On the p-particle sector, R is diagonal in the Slater basis:
R |A⟩ = s_p · ε(A) · |A⟩ with ε(A) = Π_{k∈A} (−1)^(k+1) and s_p = (−1)^(p(p−1)/2)
a sector constant. On a block cell therefore
R-parity(|A⟩⟨B|) = s_p s_q · ε(A) ε(B), and the RELATIVE parity of two cells of
one block is ε(A)ε(B) against ε(A')ε(B').

*Proof.* R on spin states is the bare bit reversal, |a⟩ ↦ |rev a⟩, no sign. Let
R_f be the fermionic one-body reflection, R_f c†_z R_f⁻¹ = c†_{N−1−z},
R_f|0⟩ = |0⟩. On a fermionized basis state c†_{z₁}···c†_{z_p}|0⟩ (sites
ascending), R_f produces the descending product on the reflected sites;
reordering p factors costs (−1)^(p(p−1)/2), so R = s_p·R_f on the sector with
s_p = (−1)^(p(p−1)/2). On modes, R_f b†_k R_f⁻¹ = Σ_z ψ_k(N−1−z) c†_z =
(−1)^(k+1) b†_k (the wavefunction form is PROOF_C1_MIRROR_SYMMETRY,
Ingredient 2; the operator form is SLOW_MODE_R_PARITY's "R is diagonal on
Bogoliubov modes"), and R_f fixes |0⟩, so R_f|A⟩ = ε(A)|A⟩ with no
reordering. The cell statement is the bra-ket product (the operator-space
product rule of SLOW_MODE_R_PARITY, now with its Slater form); the bra-side
conjugation is vacuous because the sine modes and the eigenvalues ±1 are
real; s_p s_q is one constant for the whole block and cancels from every
relative parity. ∎

This is the general-p form of F89_SEED_EXISTENCE_REDUCTION's Lemma 3 (which
proved the p = 3 determinant case); nothing else in the repo states it, and the
endpoint gate measured parities numerically instead.

## §4 One-body reach, and the equation a connection must solve

N_l on the block is (one-body) ⊗ (one-body): as a superoperator,
N_l = n̂_l ⊗ (1 − n̂_l) + (1 − n̂_l) ⊗ n̂_l with n̂_l = c†_l c_l, so a matrix
element between pair states factorizes into one-body density elements,

    ⟨A'| n̂_l |A⟩ · ⟨B'| 1 − n̂_l |B⟩  +  ⟨A'| 1 − n̂_l |A⟩ · ⟨B'| n̂_l |B⟩.

Since n̂_l = Σ ψ_j(l) ψ_k(l) b†_j b_k moves at most one mode, the element
vanishes unless A' differs from A by at most one orbital AND B' from B by at
most one orbital (the one-body selection rule, Slater-Condon). For two pair
states in the SAME eigenspace Ω this leaves exactly three shapes:

- **(0,0)** A' = A, B' = B: the diagonal, parity-preserving trivially.
- **(1,0)** A' = A∖{i}∪{m}, B' = B: equal eigenvalue forces e_i = e_m, dead
  because the comb is simple (cos is strictly decreasing on (0, π)); (0,1)
  symmetrically.
- **(1,1)** A' = A∖{i}∪{m}, B' = B∖{j}∪{n}, i ≠ m, j ≠ n: equal eigenvalue
  forces e_i − e_m = e_j − e_n, i.e. the **pair-sum collision**

      e_i + e_n = e_m + e_j.

So everything rests on which pair-sum collisions the comb allows.

## §5 The pair-sum lemma at M = 9, and the parity bookkeeping

**Lemma (pair-sum injectivity up to the chiral zeros).** Let M = 9,
e_k = cos(kπ/9), k = 1..8. If e_x + e_y = e_u + e_v with multisets
{x, y} ≠ {u, v}, then x + y = 9 = u + v (both sides are chiral pairs, both sums
zero). In particular the pair-sum map on 2-multisets is injective away from the
single chiral-zero group {1,8}, {2,7}, {3,6}, {4,5}.

*Proof.* If the multisets share an element, cancel it: the remainder is a
single-cosine equality, dead by simplicity. So assume {x, y} ∩ {u, v} = ∅. In
roots, with ζ = e^(iπ/9) of order 18 and −1 = ζ⁹,

    W = ζ^x + ζ^(−x) + ζ^y + ζ^(−y) + ζ^(9+u) + ζ^(9−u) + ζ^(9+v) + ζ^(9−v) = 0,

a vanishing sum of weight 8 of 18th roots with nonnegative coefficients (the
PROOF_F129 §2 reduction, at pair weight; exponents may repeat, e.g. when
x + u = 9, which is a coefficient 2, not a cancellation, and repeated terms
are within the decomposition machinery: PROOF_F129 §4's imported weight ≤ 12
classification is height-1, closing this subtlety there in the same words). The order is
18 = 2·3², two primes, so by (LL2) of PROOF_F129 §3 every minimal piece of any
decomposition of W is a rotated 2-cycle R₂ = {ζ^c, ζ^(c+9)} or a rotated 3-cycle
R₃. The weight 8 is not divisible by 3, so every decomposition (8 = 2+3+3 or
2+2+2+2) contains at least one R₂: two exponents of W differing by 9 mod 18.
Enumerating the pair types of the exponent multiset
{±x, ±y, 9±u, 9±v}: x with −y differ by 9 iff x + y = 9; (9+u) with (9−v)
differ by 9 iff u + v = 9; every other type demands x − y ≡ 9 or 2x ≡ 9
(impossible in range or by parity), x ≡ u mod 18 (the excluded shared element),
or x + u ≡ 0 mod 18 (impossible, 2 ≤ x + u ≤ 16). So x + y = 9 or u + v = 9;
either one makes its side sum to zero, hence the other side sums to zero too,
and e_u = −e_v = e_(9−v) forces u = 9 − v by simplicity. (The repetition case
x = y with u ≠ v runs the same way and ends in e_x = 0, which the odd comb
does not carry: no zero mode at M = 9; with u = v as well, both live clauses
die at once and W has no R₂ at all, dead one step earlier.) ∎

**Parity bookkeeping.** A (1,1) connection multiplies the cell parity ε(A)ε(B)
by (−1)^(i+m) · (−1)^(j+n). The lemma leaves exactly two solution shapes, and
both are even:

- **trivial** {i, n} = {m, j}, i.e. i = j and n = m (i = m is excluded):
  the factor is (−1)^(2(i+m)) = +1;
- **chiral** i + n = 9 = m + j: the factor is (−1)^(i+m+j+n) = (−1)^18 = +1.

The (0,0) and dead (1,0)/(0,1) shapes preserve parity trivially. So no matrix
element of N_l connects cells of opposite relative parity inside any Ω:
Π_Ω N_l Π_Ω is parity-block-diagonal, for every l, on every eigenspace of every
block. That is the Theorem. ∎

Two remarks the measured picture asked for. First, the generic R-odd diagonal
contrast of the endpoint gate (0.119) is not touched by this theorem: a generic
diagonal on cells is a MULTI-body operator, outside the one-body reach of §4, so
nothing forces its cross block, and indeed it does not vanish; the theorem
explains why the site densities are special. Second, the mixed eigenspaces
themselves track the same arithmetic one rung up: gate section (4) buckets
Slater pairs in mode space across the frontier and finds mixed collisions
absent at M = 6, 7, 8, 10, 11, exactly the endpoint counts (20, 21, 50) at
M = 9, and abundant at M = 12 (155, 307, 609), so on the tested frontier they
fire exactly at 3|M with M > 6, the same threshold as F129's door (M = 6
carries the 3 but stays empty, its only collisions being chiral and
parity-even, F146's recorded exception one rung up; a derivation of the
mixed counts from F129's families is not attempted here). The collisions that MIX parity are
there, in one-body terms far apart; what the theorem says is that every
collision within one-body REACH preserves parity. Mixing and coupling live at
different Hamming distances in mode space.

## §6 Verification

One gate,
[`simulations/mixed_space_reflection_gate.py`](../../simulations/mixed_space_reflection_gate.py)
(output `simulations/results/mixed_space_reflection_gate.txt`, VERDICT last
line), four sections. (1) The FORCED combs, exact: the pair-sum census of
{1..M−1} bucketed over GF(p) for two primes p ≡ 1 (mod 2M) per M
(distinctness mod one prime is proof of exact distinctness, the
LevelCollisionCensus soundness direction; a mod-p single-particle degeneracy
returns a sentinel that fails the comparison rather than passing), gated ==
the chiral-zero group alone at M = 6, 9, 10, 11, 13, 14, 16, 25, 27; the
M = 9 row is the lemma of §5. (2) The DOOR combs: a non-chiral colliding
group present and a named parity-odd instance exact at M = 12, 15, 18, 21,
plus the 6|M family identity of §7 at every admissible k for
M = 12, 18, 24, 30 and the empty admissible range at M = 6. (3) The
sharpened prediction at N = 8, Δ = 0 on (1,3), (1,4), (2,4), (2,3) and the
population-carrying (4,4): per-l cross-parity spectral norms of comp(N_l) on
every parity-mixed eigenspace (a scalar one has no cross block), and a
random combination Σ c_l N_l, all gated at 64·eps·‖A‖/gap (measured worst
ratios 0.001-0.011 across 207 mixed spaces, 91 of them the endpoint trio's);
the generic multi-body diagonal contrast gated from below (0.113-0.165; this
gate's contrast object is a plain random diagonal, the endpoint gate's 0.119
an R-odd-projected one, two different generics of similar size); and the
smallest true difference-spectrum gap gated far above the grouping tolerance
(measured 4.46·10⁻² against gtol = 10⁻⁹, so grouping cannot merge distinct
eigenspaces). (4) The frontier census in mode space, exact over two primes:
Slater pairs bucketed by E_A − E_B with their parities, mixed-bucket counts
zero at M = 6, 7, 8, 10, 11, exactly (20, 21, 50) at M = 9, populated at
M = 12.

## §7 Other N: where the same argument forces the law, and where it stops

The proof used Δ = 0, uniform J, and two facts about the comb: simplicity,
and pair-sum injectivity up to parity-even solutions. The chiral solutions
are parity-even at EVERY M (their factor is (−1)^(2M)), so what decides is
whether the comb carries any NON-chiral pair-sum coincidence; the subset
face of that input is F146's chiral-only notion at rung 2, while the
MULTISET face is what §5 proves at M = 9 and the gate certifies per comb.
Gate sections (1), (2), (4) certify each claim below exactly.

- **M an odd prime power.** The §5 argument runs unchanged after 9 → M,
  18 → 2M (order 2M = 2p^a, two primes, so (LL2) still yields R₂/R_p
  pieces; the branch 2x ≡ M stays dead because M is odd; for p ≥ 5 the
  weight-8 decompositions are all-R₂, for p = 3 at least one R₂ survives
  the count): rung-2 is chiral-only and the law is FORCED. Beyond M = 9 this
  covers every odd prime M (N = 10 and N = 12 among them; for M = p,
  PROOF_SCALAR_COUNT §4 proves the subset form under its deletion
  hypothesis) and the higher odd prime powers M = 25, 27, 49, ...
  (N = 24, 26, 48, ...), gated at M = 11, 13, 25, 27.
- **M = 2p and M = 2^a.** Forced at the gated combs M = 10, 14, 16, where
  the unrestricted census is exact; for the full families the forced status
  is EXPECTED, not derived here. The reason a citation does not close them:
  PROOF_SCALAR_COUNT §4's chiral-only is a statement about ℓ-SUBSETS under
  its deletion hypothesis, while the shapes of §4 above are MULTISETS (m and
  j sit on opposite sides of the cell) and include configurations outside
  that hypothesis, e.g. the doubled zero mode {1, M−1} ~ {M/2, M/2}
  (parity-even, index sum 2M, joining the chiral-zero group) and, in
  principle, cross-chiral shapes {x, y} ~ {M−x, v} whose parity is not
  fixed; the gated censuses show none of the latter at the three combs.
- **The doors 6|M (M > 6), 15|M, 21|M.** Exactly PROOF_SCALAR_COUNT §7's
  measured rung-2 doors, and each is a PARITY door: the forcing fails with
  an explicit parity-odd coincidence, gated exactly at M = 12, 15, 18, 21
  ({1,9} ~ {5,6}, {1,5} ~ {3,4}, {2,14} ~ {8,9}, {3,9} ~ {6,7}, index sums
  odd). At 6|M the mechanism is a family:
  e_(k−M/3) + e_(k+M/3) = e_k + e_(M/2), the 60° identity plus the zero
  mode, holding at every admissible M/3 < k < 2M/3, k ≠ M/2, with parity
  factor (−1)^(3k+M/2), parity-odd for at least one admissible k at every
  gated M (2 of 2, 2 of 4, 4 of 6, 4 of 8 at M = 12, 18, 24, 30; the
  excluded k = M/2 breaks the naive alternation); M = 6 itself is F146's
  recorded divisible exception,
  the family's admissible range being empty there. At 15|M and 21|M there is
  no zero mode (M odd) and the door has a different shape, exhibited by the
  gated instances. First chain instances: N = 11, 14, 20.
- **The frontier, in mode space.** Gate section (4): mixed difference-buckets
  on (1,3), (1,4), (2,4) are absent at M = 6, 7, 8, 10, 11, exactly
  (20, 21, 50) at M = 9, and populated at M = 12 (155, 307, 609): on this
  frontier, mixed collisions fire exactly at 3|M with M > 6 (the M = 6
  exception again), and N = 11 holds BOTH
  ingredients, abundant mixed spaces and a parity-odd rung-2 coincidence.
  What fails there is the FORCING: this proof makes no claim that the law
  breaks at N = 11, only that its protection ends; whether a parity-odd
  coincidence lands inside a mixed ad_H eigenspace with a nonvanishing
  one-body element is the arc's open NextStep 4, now with a concrete
  candidate and a well-posed first N.
- **Δ ≠ 0.** Nothing here survives: the argument is Slater-additive through
  and through (the resonance note's own fence on F129), and the measured
  Δ = 1 breaks at N = 6 show the conclusion itself is false there.
- **Non-uniform J.** The comb is gone. For palindromic J the theorem layer
  ([H, R] = 0) stands and the endpoint gate measured zero mixed collisions,
  so there is nothing for this proof to add; for non-palindromic J even the
  theorem layer is void.
