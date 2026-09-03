# PROOF: The Crack and the End-Pair Anisotropy Carry One Blindness Locus, and at Odd N One Count

**Registry:** [F162](../ANALYTICAL_FORMULAS.md#f162), for §(i) alone; §(a) to §(h) carry none.
§(a) to §(f) are harvested from [F157](../ANALYTICAL_FORMULAS.md) and
[F160](../ANALYTICAL_FORMULAS.md) and add no closed form those two do not already carry. §(g) does
add closed forms, though fewer than it might look: the crack's count at u = +1, ⌊(N−1)/2⌋ at every
seat, is the companion page's gate B3 and is committed, and what §(g) adds there is the PROOF;
what is new is the other end, ⌈(N−1)/2⌉ at u = −1, and the centre seat's (N−1)/2 at every coupling
on both axes. What §(g) adds no closed form for is the crack's count BETWEEN the ring
ends; §(h) closes that at even N as 2·#{k odd : Δ_k = t}, and at odd N away from the
reflection-fixed centre seat it is F157's own #{k}, that seat paying (N−1)/2 at every coupling by
§(g). What neither closed was the two sector halves-resultants of §(g)'s Corollary, and §(i)
closes those: the sign is (−1)^binom(p+1,2) in the fold coordinate alone, each halves-resultant is
its own constant times the product over the non-pole roots of its own comb, and what is factorised
is F157's generator up to a second sign given in closed form. That is what F162 registers,
and it registers §(i) and nothing else. Whether the §(a) to §(f) harvest earns a number is still
left open rather than decided here.
**Status:** Tier 1 derived. Lemma 1 is a matrix identity at every N and every complex u.
Lemma 2a is the fence-free Cramer argument of [The Seat That
Cuts](../../experiments/THE_SEAT_THAT_CUTS.md) §7 read per eigenvalue; Lemma 2b, the sector split,
is this file's. Lemma 3 is Lemma J of the node-lemma proof. Lemma 4 is this file's own
boundary-system argument and is deliberately independent of F160, whose simplicity clause is
fenced to u ≥ 0 while half of the u axis here is negative. Lemmas 6 and 7, the split of the struck
characteristic polynomial and the count decomposition it carries, are this file's, as is the
strengthening of Lemma 5 from sets to counts. Lemma 8 is a TRANSPORT and not a new sign: the
reflection parity (−1)^{k+1} of a mode is owned in six places for the chain's own comb, and what
is this file's is that the blind mode of an END-DETUNED chain is still a single sine, so the law
survives both the change of modulus to N_node and the loss of uniformity, and the index it then
puts on §(g)'s two summands. Lemmas 9, 10 and 11 and the four corollaries of §(i) are this file's;
the Cassini identity Lemma 9 turns on is `PROOF_CRACKED_RING_EXACT_CURVE` §(b)'s and is cited
there. Both Theorems are exact at every N ≥ 3 and both
parities. **Date:** 2026-09-02, sections (h) and (i) 2026-09-03. **Authors:** Thomas Wicht, Claude (Opus 5). **Script:**
[`simulations/blind_seat_two_axes_proof.py`](../../simulations/blind_seat_two_axes_proof.py), 84
checks under 81 labels in eight blocks (L2a fires twice, L2b three times), exact in sympy, about
two minutes measured quiet under sympy 1.14.0.
Every resultant in block W that carries a LAW is a Sylvester determinant the file builds itself,
sympy's own routine appearing there only where it is the object under test; W0 to W0d say why. Run
committed at
[`blind_seat_two_axes_proof_run.txt`](../../simulations/results/blind_seat_two_axes/blind_seat_two_axes_proof_run.txt).
**Typed:** §(i) alone, as
[`BlindSeatSectorFactorisationClaim`](../../compute/RCPsiSquared.Core/Symmetry/BlindSeatSectorFactorisationClaim.cs),
Tier 1 derived, parents `SeatCutBlindnessClaim` (F157, the locus it factorises) and
`CrackedRingExactCurveClaim` (F160, the Cassini step); live at
`dotnet run --project compute/RCPsiSquared.Cli -- inspect --claim BlindSeatSectorFactorisationClaim`.
It recomputes Lemmas 9, 10 and 11 and Corollaries 10a, 10b and 11b at inspect time, and
deliberately not on this script's road:
resultants by the Euclidean remainder sequence rather than as Sylvester determinants, one integer
knob at a time with an interpolation whose held-back samples make an underestimated knob-degree
red, and the sector parity as the exact remainders S_{n+1} + 1 ≡ 0 (mod β_E) and
S_{n+1} − 1 ≡ 0 (mod β_O) rather than as readings in the node field. On that second road it
reproduces four of block W's pinned populations: the 72 seats, the fold profile
{1:22, 2:18, 3:14, 4:10, 5:6, 6:2}, the 20 sector readings where a pole drops the knob-degree, and
the 16 that carry a repeated factor. What it does NOT carry is what §(i) itself marks as read
rather than derived, W5's identification of the sector halves and the step from a multiplicity to
Corollary 8a's b_E, together with Corollary 11a's seat-versus-fold reading, the comparison
against F157's four committed rows with its pinned per-row scales, the K1B_SIGNS convention
reading and W0 to W0d, all of which stay in this script.
**Closes:** the two questions carried by open item 1 of [The Blind Seat on the
Road](../../experiments/THE_BLIND_SEAT_ON_THE_ROAD.md). Its §(c) was an observation with an
identified mechanism at odd N = 5..17 and at 44 of 66 even-N seats; that page states the relation,
and §(a) to §(f) here state why it holds, why it fails at the other 22, and why it holds at every
N rather than at the ones that were reached. The item's second question asked whether a seat is
blind by the same AMOUNT at a shared point, and recorded that the construction there cannot see
it; §(g) answers it, and shows that the construction can see it after all, doubled. §(g) then
generated one more question on that page, which sector each of F157's node indices k belongs
to, and §(h) closes that too, so the even-N count is closed-form. §(g) and §(h) then left two
constants unidentified between them, the ±1 of §(g)'s Corollary and the rational constant of
§(h)'s closing paragraph, both carried as open by the arc `the_forced_and_the_met`; §(i) closes
both, and finds on the way that §(g)'s ±1 had no convention-free answer to give.

## Statement

Let A be the adjacency matrix of the open N-site path, the single-excitation XY block, and put

    P = E₀,ₙ₋₁ + Eₙ₋₁,₀        the wrap bond, F160's road, u = 0 open chain, u = 1 ring
    D = E₀₀ + Eₙ₋₁,ₙ₋₁         the end-pair diagonal, F157's anisotropy axis

For a seat j write **locus(j)** for the set of real knob values at which seat j is blind, that is
the real root set of Res_x(χ(H), χ(H_j)) with H_j the matrix H with row and column j struck.
Let E(j) and O(j) be the corresponding sets computed inside the reflection-even and
reflection-odd sector of A + tD, and let −S denote {−s : s ∈ S}. Then

> **Theorem.** u-locus(j) = E(j) ∪ (−O(j)) ∪ {+1, −1} and Δ-locus(j) = E(j) ∪ O(j), at every
> N ≥ 3 and every seat. Hence
>
> - at **odd N**: u-locus(j) = Δ-locus(j) ∪ {+1, −1}, at every seat, every odd N;
> - at **even N**: u-locus(j) = E(j) ∪ {+1, −1} and Δ-locus(j) = E(j) ∪ (−E(j)), so the same
>   relation holds **exactly when** E(j) ∪ {+1, −1} is negation-closed.

The odd-N line is the companion page's §(c) headline, there measured at N = 5..17. The even-N line
is a criterion rather than the relation, and it is that page's gate C8, which measured *every
break is a non-closure* over N = 6..16 and called §(c) *"an observation with an identified
mechanism, not a theorem"*, having no mechanism for the converse; that page now cites this file
for it. It also names the breaking seats: over N = 6..16 the criterion predicts 22 seats, and they
are the 22 the companion gate found.

A locus is a set, and F157's blind count at one of its points is a multiplicity, so the Theorem
leaves open whether a seat blind at the same knob value on both axes is blind there by the same
AMOUNT. Writing b_E(j; t) and b_O(j; t) for the number of eigenvalues in each sector whose
eigenvector vanishes at j, §(g) settles that too:

> **Theorem (count).** At every seat j and every real t,
>
> - at **odd N and |t| ≠ 1**: blind(j; H_u(t)) = blind(j; H_Δ(t));
> - at **even N and |t| ≠ 1**: blind(j; H_u(t)) = 2·b_E(j; t) while
>   blind(j; H_Δ(t)) = b_E(j; t) + b_E(j; −t);
> - at **t = ±1**, both parities: blind(j; H_u(t)) is ⌊(N−1)/2⌋ at u = +1 and ⌈(N−1)/2⌉ at
>   u = −1, at every seat, whatever the chain pays there.

So the two axes carry one count wherever they carry one locus, over everything measured, but the
REASON is one only at odd N and off the ring ends; at even N the agreement is a coincidence of two
sector counts, and where it fails the locus fails with it. At even N the crack pays the even
sector TWICE where the chain pays each sector once, and at the ring ends the crack is blind at
every seat at once. The parity that split the loci splits the counts the same way, through the
same commutation ΣR = (−1)^{N−1}RΣ.

The two summands b_E and b_O are still measurements at that point. §(h) puts F157's own index
on them: the blind mode at node index k is R-even for odd k and R-odd for even k, so
at every seat with N_node ≥ 2, b_E counts the odd k landing on a value and b_O the even ones, and
the even-N count becomes 2·#{k odd : Δ_k = t} off the ring ends, a closed form.

## What the repo already held, store by store

Swept twice. The second sweep, 2026-09-03 and for §(h) alone, is recorded inside that
section rather than here, because what it returned changed that section from a discovery into a
transport: the sign law it needed is owned in six places, and only the modulus is new. The
first, for §(a) to §(g), was 2026-09-02 by three agents, one per primitive (the sector reduction, the blindness
criterion, level crossings between symmetry sectors), each over the named stores rather than over
the arc's vocabulary, then corrected by three empty review rounds which moved four claims here
against the author. What was drafted as a new lemma turned out to be owned in folded form; what
was drafted as the hard step turned out to be a corollary; a blanket "nothing" over three stores
turned out to repeat the exact defect `docs/CAUGHT_ERRORS.md` logged for the companion page the
same day; and Lemma 4's original citation turned out not to reach the negative half of its own
axis. All four are recorded rather than quietly absorbed.

- **`docs/proofs/PROOF_CRACKED_RING_EXACT_CURVE.md`** (F160), Corollary B: **owns the folded
  half of Lemma 1.** *"R commutes with H(u) and H splits into an R-even and an R-odd block.
  Folding the ring along R, each block is a tridiagonal matrix with every off-diagonal entry
  nonzero ... and a diagonal carrying ±u at the crack end"*, with its helper `folded_blocks`
  (whose body is literally `e[0, 0] = u; o[0, 0] = -u`) gated symbolically in u at N = 3..9 as
  gate P9 of
  [`cracked_ring_exact_curve_proof.py`](../../simulations/cracked_ring_exact_curve_proof.py). So
  Lemma 1 is not a new lemma. It is one comparison line: the same fold applied to the
  anisotropic OPEN chain puts the SAME +t in both folded blocks, which is the crack's even block
  at u = +t and its odd block at u = −t. That file's simplicity clause also gives a route to
  Lemma 4, and §(d) explains why this file does not take it.
- **`experiments/THE_SEAT_THAT_CUTS.md` §7**: owns Lemma 2a in all but the phrasing, fence-free
  and for arbitrary real symmetric H. *"By Cramer's rule the (j, j) entry of the adjugate of
  xI − H IS the characteristic polynomial of the principal submatrix, so
  χ(H_cut)/χ(H) = [(xI − H)⁻¹]_{jj} = Σ_λ ‖P_λ e_j‖²/(x − λ) ... An eigenvalue the seat projects
  onto keeps one factor fewer in the numerator than in χ(H); one the seat misses keeps at least
  as many, and the gcd takes the minimum."* That sentence is the per-eigenvalue bookkeeping
  Lemma 2a needs; the lemma is its restatement, not a derivation.
- **`docs/proofs/PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md`**: owns the same counting as
  Corollary C, `blind(j) = Σ_λ (m_λ − s_λ)` with `s_λ ∈ {0, 1}`, and owns **Lemma J**, the
  simplicity of an unreduced Jacobi spectrum, which is the whole of Lemma 3. Its own per-λ
  eigenvector-vanishing statement (J3) is fenced to the two-halves Jacobi case, which is why
  Lemma 2a is taken from §7 above and not from here.
- **`docs/proofs/PROOF_COLLISION_GAP_ODD_ORDERS.md`** (F161) §(b), the lemma the companion page
  cites as *"Lemma B"*: the wrap bond in the mode basis, `⟨ψ_k|V|ψ_l⟩ = a_k·a_l·(η_k + η_l)`,
  rank one per sector. That is Lemma 1 seen in the mode basis, and it is what let the companion
  page write *"two end-pair objects differing only in that factor"*. That page has since taken
  the last step itself, citing this file.
- **`docs/ANALYTICAL_FORMULAS.md`**: F157 owns the Δ axis in closed form (Δ_k, the resultant
  packaging P_j, the live witness `inspect --root blindlocus`). F160's entry owns the u axis curve
  and, correcting a draft of this bullet that said otherwise, **does carry the reflection
  sectors**: *"its simplicity clause (G = 2AB, the two reflection sectors, each a nonvanishing
  prefactor times an unreduced Jacobi block's characteristic polynomial, and a Bézout identity
  forbidding a common zero unless u² = 1: the spectrum is simple for every u ≥ 0 except u = 1)"*,
  and it carries the fence this file crosses, *"Not the blind seat: THE_SEAT_THAT_CUTS's open item
  asks for a detuned bond under a seat cut, a different object"*. F161's entry names the sectors
  too. What no entry carried, before the commit that landed §(a) to §(f), was a statement relating
  the two PERTURBATIONS; F157's entry carries one now.
- **`compute/MirrorWorld/BlindSeat.cs`**: owns the Δ-locus as a resultant at a named non-centre
  seat, *"at N = 9 seat 1 that resultant factors as `128*D*(D-1)*(D+1)*(D^2-3)`, so the seat is
  blind at Delta = 0, +-1 and +-sqrt(3) and nowhere else"*, and it owns §(a)'s third convention
  too, the centre seat's two principal submatrices being *"conjugate by the chain reflection, so
  they carry the same characteristic polynomial and share every root, at any Delta and any
  palindromic profile (their resultant is identically zero)"*. That is §(a)'s object in the
  small. `Crack.cs` names F160's factorization in passing and does no fold.
- **`compute/RCPsiSquared.Diagnostics/Foundation/DiabolicReflectionParityWitness.cs`**: owns the
  SHAPE of Lemma 5 on a different object, the (SE,DE) coherence block: *"At EVEN nBlock the two
  sectors have equal dimension ... At ODD nBlock the reflection-fixed central site makes
  dim(R-even) − dim(R-odd) = (nBlock−1)/2 ≠ 0; the dimension mismatch forbids the cross-pairing,
  forcing the antiunitarity to act within each sector."* Same two operators in kind, an
  antiunitary and a site reflection, and the same preserve-versus-swap dichotomy by the parity
  of N. Lemma 5 is that dichotomy on the single-excitation chain; it is a second instance, not a
  new shape.
- **`docs/GLOSSARY.md`** holds the blind seat itself, under the headword *The blind seat (F157,
  August 2026)*, with `blind(j) = N − dim Krylov(e_j)` and both gcd laws: §(a)'s vocabulary is
  already written there. What it has no headword for is reflection sector, end pair, crack or
  comb.
- **`reflections/ON_LEAVING_THE_CIRCLE.md`** holds Lemma 5's dichotomy in prose, *"in an even
  chain the two halves are each other's reflection ... An odd chain has a middle seat that
  reflection cannot move"*, about the unmirrorable seat rather than about a locus.
- **`hypotheses/DIABOLIC_BY_INTEGRABILITY.md`** turns on the same Δ axis and records the
  opposite case to §(d)'s, *"No symmetry puts the two modes in different sectors so they cannot
  couple: the site-reflection R commutes but both modes are R = +1"*: a degeneracy INSIDE one
  sector, which is exactly what Lemma 3 excludes here.
- **Nothing** from `recovered/` and from `fw.Confirmations`, both checked entry by entry; the
  nearest confirmation is 24, on the n = 9 comb, a neighbour and not a hit.
- **A near-miss worth naming, because it looks like a hit:** `experiments/CROSSING_TAXONOMY.md`,
  `SUBSYSTEM_CROSSING.md` and every "crossing" row in the glossary are the CΨ boundary crossing,
  not a level crossing. The repo owns no general statement of when levels in different
  reflection sectors cross.

**The COUNT was swept separately, on 2026-09-02 by three further agents**, one per store, before
§(g) was drafted, because it is a different primitive from the locus and the first sweep had not
looked for it.

- **`docs/ANALYTICAL_FORMULAS.md`**: F157 owns the count on the Δ axis and states the rule §(g)
  needs, *"the multiplicity of a root of P_j IS the blind count at that Δ"*, with the worked
  double root at N = 11 seat 2 where counting distinct roots would break the XY law. What it did
  not own, before this file, was any count on the u axis or any N_node analogue there; the F157
  entry carries §(g)'s reading now. F160 and F161 own no blindness count at all.
- **The typed layer**: `SeatCutBlindnessClaim` (F157) carries the count at the two endpoints in
  closed form, as `BlindHeisenberg` and `BlindXy`; `SeatBlindnessDeltaLocusWitness` carries it
  along the Δ axis, as `BlindAtRational` and `BlindAtAlgebraic`, the latter explicitly a
  multiplicity by exact integer polynomial division. Every matrix builder in both is an OPEN
  chain, and **no C# member anywhere computes a blind count as a function of u**. `Crack.cs` owns
  `road` and `departures` and no Krylov space.
- **The operational layer**: `blind_seat_two_axes_proof.py` and `blind_seat_on_the_road.py` both
  define a count on a numeric matrix, and both had read it on the u axis. Where they had read it
  matters. Between the ends they had read it at every seat once, at the single generic u = 1/3 of
  gate C3b, which lies on no NON-CENTRE seat's locus, and at six couplings at the reflection-fixed
  CENTRE seat, where the count cannot move. Seats that CAN lose their blindness had been counted
  where they are blind, at the ring end u = +1 by that page's gate B3 and at u = 0 by its A1 and
  B4. What had never been done, on either page, is to compare a count on one axis with a count on
  the other.
- **`docs/proofs/PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md`** owns Corollary C,
  `blind(j) = Σ_λ (m_λ − s_λ)`, which is the whole arithmetic Lemma 7 needs, and its Lemma J is
  what makes each sector spectrum simple. Lemma 7 is that corollary read one reflection sector
  down, exactly as Lemma 2b reads the kernel.
- **The OpenArcs registry**, which is the store this file's own arc lives in and the one
  `CLAUDE.md` names by path: `the_forced_and_the_met` held the prior §(g) overturns, *"whether
  they carry the same COUNT at a shared point is untouched, the resultant being a multiple of
  F157's P_j SQUARED."* That sentence WAS the fourth copy of the P_j-squared claim and the reason
  the question was live; this change rewrites it, and the three on the companion page with it.
- **`docs/GLOSSARY.md`** owns the count's definition, `blind(j) = N − dim Krylov(e_j)`, under the
  headword *The blind seat*, and no comparison of counts across two axes; it is a hit on the
  primitive and not on the question.
- **`experiments/`, null results included**: the question itself is from
  `THE_BLIND_SEAT_ON_THE_ROAD`, `THE_SEAT_THAT_CUTS` owns the Cramer bookkeeping and the
  Jacobi-end statement §(g)'s centre-seat corollary rests on, and `THE_BLIND_SITE` owns the
  direct-sum form of the struck matrix. No page there compares a count across two axes.
- **`docs/CAUGHT_ERRORS.md`**: its 2026-09-02 entries bear on the method rather than the result.
  The third holds the two sympy traps and the controls that could not pass, all four items of it
  about this arc, and they are respected below.
- **Nothing** on a count comparison across the two axes, from `hypotheses/`, `reflections/`,
  `recovered/`, `review/`, the derivations D01 to D10, or the hardware Confirmations registries in
  either language. And **nothing anywhere derives** the P_j-squared claim, which four passages
  state in near-identical words: three on the companion page (its gate-C1 paragraph, its scope
  bullet, its old open item) and one in the arc. §(g)'s last corollary is the first derivation of
  it, and it corrects the conclusion those passages drew from it.

## (a) The object, and three conventions the numbers depend on

Two one-parameter families on the same N sites, both real symmetric, both commuting with the
chain reflection R: l ↦ N − 1 − l:

    H_u(t) = A + t·P        the cracked ring
    H_Δ(t) = A + t·D        the end-pair-anisotropic open chain

The scale is the companion page's, hop 1 and a bare t, which differs from F157's committed
convention (hop 2, 2Δ per end, a Δ(N−5)·I shift) by a common scaling and a common shift. Neither
moves a coincidence between χ(H) and χ(H_j), and the companion gate's C1 checks that against
F157's committed table rather than asserting it.

A **seat** j ∈ {0, …, N−1} is where the watching sits and blind(j) is F157's count, defined in
`docs/GLOSSARY.md`. The **locus** is the real root set of Res_x(χ(H), χ(H_j)) in t, a SET; F157's
count at a given point of that set is a multiplicity, out of scope until §(g), which is where
it is met. Three conventions,
because without them not one number below can be reproduced:

1. **How E(j) and O(j) are computed.** Both are the same strike-and-resultant construction
   applied INSIDE a sector: the sector's matrix in the pair basis {e_p + η·e_{N−1−p}}, with the
   row and column at the fold coordinate min(j, N−1−j) struck. A sector vector's amplitude at site
   j is
   ± that pair coordinate, the sign being −1 only in the odd sector at j > N−1−j, and vanishing
   is unaffected by a sign, so the criterion transports.
2. **The sector block is not symmetric.** The pair basis is orthogonal but not orthonormal, so
   the block is the conjugate of a real symmetric matrix by the positive diagonal
   diag(√‖·‖). That conjugation changes neither χ nor which coordinate of an eigenvector
   vanishes, which is all Lemma 2a reads, so Lemma 2a applies to it; but the block itself does
   not satisfy Lemma 2a's stated hypothesis, and saying so is not optional.
3. **The reflection-fixed centre seat of an odd chain is a degenerate case of the definition.**
   There every R-odd eigenvector vanishes (ψ_j = −ψ_j), so O(j) is all of ℝ and
   Res_x(χ(H), χ(H_j)) is identically zero: the locus is ℝ on both axes and is not the root set
   of anything. The Theorem holds there because both sides are ℝ, not because a polynomial was
   compared. The companion page counts 7 such seats in its odd range, the script carries them as an
   explicit `'ALL'` branch, and `BlindSeat.cs` already records the identically vanishing
   resultant.

## (b) Lemma 1, the sector reduction

> **Lemma 1.** On the R-even subspace P acts as D; on the R-odd subspace P acts as −D. Hence
> H_u(t)|even = H_Δ(t)|even and H_u(t)|odd = H_Δ(−t)|odd, as matrices, at every N and every t.

*Proof.* Both P and D commute with R, so both preserve the two sectors and the restrictions are
defined. Let ψ be R-even, so ψ_{N−1} = ψ_0. Then (Pψ)_0 = ψ_{N−1} = ψ_0 = (Dψ)_0 and
(Pψ)_{N−1} = ψ_0 = ψ_{N−1} = (Dψ)_{N−1}, and both vanish at every other site. So Pψ = Dψ. Let ψ
be R-odd, so ψ_{N−1} = −ψ_0; the same two lines give (Pψ)_0 = −ψ_0 = −(Dψ)_0 and
(Pψ)_{N−1} = ψ_0 = −ψ_{N−1} = −(Dψ)_{N−1}, so Pψ = −Dψ. ∎

The content is one substitution, ψ_{N−1} = ±ψ_0, and it is why the two axes can meet at all:
**the cracked ring is, sector by sector, the anisotropic open chain, at +u and at −u.** At odd N
the reflection-fixed site is untouched, neither P nor D having any support there.

Gate L1a checks the matrix identity symbolically in u at N = 3..14 in both sectors, exact zero
matrix. L1b is its anti-vacuity partner on the same code path: the WRONG sign in the odd sector
disagrees at every N, so L1a is not comparing a matrix with itself. L1c checks the mechanism
rather than a difference of values: R commutes with both P and D and does NOT commute with a
one-end diagonal, both directions asserted at every N.

## (c) Lemma 2, the locus decomposition, and where the striking goes

The obstacle the open item ran into is real and must be met head on: F157's count is defined by
striking row and column j, and **striking a non-central seat destroys the reflection**, so
Lemma 1 does not transfer to the struck matrix. It does not have to, because the striking can be
removed from the argument before the sectors are used.

> **Lemma 2a.** For real symmetric H and any seat j, λ is a common root of χ(H) and χ(H_j) if
> and only if some nonzero eigenvector of H at λ vanishes at j.

*Proof.* This is `THE_SEAT_THAT_CUTS` §7 read one eigenvalue at a time, and is repeated here
only because the two directions are used separately below. Write m_λ = dim E_λ and
s_λ = dim span{P_λ e_j} ∈ {0, 1}; the eigenvectors at λ vanishing at j are E_λ ∩ e_j^⊥, of
dimension m_λ − s_λ, because P_λ e_j spans the only direction of E_λ not orthogonal to e_j.

Should such an eigenvector ψ exist, delete its j-th entry to get ψ′ ≠ 0; for i ≠ j,
(H_j ψ′)_i = (Hψ)_i − H_{ij}ψ_j = λψ_i, so λ ∈ spec(H_j) and λ is a common root.

Conversely let λ be a common root and suppose no eigenvector at λ vanishes at j, so m_λ = s_λ,
which with s_λ ≤ 1 and λ ∈ spec(H) forces m_λ = s_λ = 1. Then the residue of the resolvent entry
at λ is ‖P_λ e_j‖² ≠ 0, so by §7's Cramer identity χ(H_j)/χ(H) has a genuine pole at λ: λ is a
root of χ(H) to one order higher than of χ(H_j), and with m_λ = 1 that means λ is not a root of
χ(H_j) at all, contradicting the assumption. ∎

Hermiticity is what the argument uses, and it already gives spec(H) ⊂ ℝ, so that a real root set
is the whole of the common-root set. What real symmetry adds here is only that the loci are read
as real root sets of polynomials with rational coefficients, which is what makes them comparable
factor by factor; both families have it.

> **Lemma 2b.** If H commutes with R, write E(j) for the knob values at which some nonzero
> R-even eigenvector vanishes at j, O(j) the same for R-odd, and **C** for the knob values at
> which H has any eigenvalue of multiplicity at least two. Then locus(j) = E(j) ∪ O(j) ∪ C, for
> every seat, and C does not depend on the seat.

*Proof.* For ⊇: E(j) and O(j) lie in locus(j) by Lemma 2a; and if dim ker(H − λ) ≥ 2 then that
kernel, a subspace of ℝ^N of dimension at least two, meets the hyperplane {ψ_j = 0} in a nonzero
vector, so by Lemma 2a again the value lies in locus(j), whatever j is.

For ⊆: let ψ ∈ E_λ be nonzero with ψ_j = 0 and write ψ = ψ⁺ + ψ⁻ along the R-invariant splitting
E_λ = E_λ⁺ ⊕ E_λ⁻. Should one part vanish, the other is a sector eigenvector vanishing at j and
the value lies in E(j) or O(j). Should both be nonzero, dim E_λ ≥ 2 and the value lies in C. ∎

**C carries no seat.** That is the force of the third term: a degeneracy blinds every seat at
once, so it appears in every seat's locus and in none of them for a reason having to do with that
seat. It is also the missing entry in the companion page's taxonomy, which splits a vanishing
mode into *forced* (a reflection fixes the seat) and *accidental* (a sine lands on a multiple of
π) and has no room for this. Stated on its own:

> **Remark.** For any Hermitian H, any λ with m_λ ≥ 2 is blind at EVERY seat.

It is one line from `THE_SEAT_THAT_CUTS` §7 (an eigenvalue of multiplicity m keeps at least
m − 1 factors in the gcd whether the seat projects onto it or not) and the sweep found it drawn
nowhere. It is not the same statement as C: C is the set of knob values where such a λ exists,
and the Remark is why that set enters every seat's locus.

Gate L2a checks the decomposition on both axes, 60 seats each over N = 4..11. Its three
mutations L2b drop C, O and E from the union in turn and are gated against the literal counts
54, 10 and 14 of the 60 crack seats.

## (d) The two degeneracy sets, which is where the axes differ

> **Lemma 3.** C_Δ = ∅ at every N and every real t.

*Proof.* H_Δ(t) is tridiagonal with every off-diagonal entry equal to 1, an unreduced Jacobi
matrix, whose spectrum is simple (Lemma J of `PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA`).
Unreducedness is a fact about the zero pattern and t
never touches an off-diagonal, so this holds at every t, including the |t| at which levels leave
the band. ∎

> **Lemma 4.** C_u = {+1, −1} at every N ≥ 3.

F160's Corollary B offers a route here, through its Bézout identity, and this file deliberately
does not take it: that clause is stated *"for every u ≥ 0 with u ≠ 1"*, and half of this file's u
axis is negative, so the citation would not reach the point (u = −1) at which the statement has
its content. The argument below is self-contained, needs no fold, no band restriction and no
prefactor relation, and holds at every real u.

*Proof.* Fix a value u of the knob and let H_u(u)ψ = λψ. The interior rows are the recursion
ψ_{l−1} + ψ_{l+1} = λψ_l for
l = 1, …, N−2, whose solution is ψ_l = ψ_1·U_{l−1} − ψ_0·U_{l−2} with U_n = U_n(λ/2) the
Chebyshev polynomials of the second kind and the conventional U_{−1} = 0, U_{−2} = −1.
Substituting into the two boundary rows leaves a 2×2 homogeneous system in (ψ_0, ψ_1),

    ⎡ −(λ + u·U_{N−3})     1 + u·U_{N−2} ⎤
    ⎣    U_{N−2} + u         −U_{N−1}    ⎦

and since an eigenvector is determined by (ψ_0, ψ_1), dim ker(H_u − λ) ≥ 2 holds exactly when all
four entries vanish. Now

    u·(1 + u·U_{N−2}) − (U_{N−2} + u)  =  U_{N−2}·(u² − 1)

so a common zero of the two anti-diagonal entries, (0,1) and (1,0), forces u² = 1 or
U_{N−2} = 0. On the branch
U_{N−2} = 0 the same two entries read 1 and u, so the first is never zero and the branch is
empty. Hence C_u ⊆ {+1, −1}.

For the converse, at λ = 2cos(mπ/N) with 1 ≤ m ≤ N−1 one has U_{N−1} = 0 and
U_{N−2} = sin((N−1)mπ/N)/sin(mπ/N) = (−1)^{m+1}, and U_{N−3} = (−1)^{m+1}·λ. All four entries
then vanish for u = −U_{N−2}, that is at u = +1 for even m and at u = −1 for odd m. Both
occur for every N ≥ 3, in counts ⌊(N−1)/2⌋ and ⌈(N−1)/2⌉. ∎

**This is the entire difference between the two axes.** The open chain is Jacobi; the ring is not,
because the wrap entry lifts it out of that form. What the companion page carried as an
unexplained appendix *"∪ {+1, −1}"* is the ring degeneracy, and it appears at every seat because C
carries no seat.

**Which multiplicity, and why the gate may certify it.** Lemma 4 bounds the GEOMETRIC
multiplicity, dim ker, and its forward direction is a polynomial identity that holds over ℂ. The
gates read the ALGEBRAIC multiplicity, off disc_x(χ). For real t both families are real symmetric
and the two coincide, which is what entitles a discriminant to certify a kernel statement; over ℂ
they part, and at N = 3, u = 2√2·i the discriminant vanishes with dim ker = 1, an exceptional
point. So the reality of the knob is named here rather than assumed. A free by-product of the
same 2×2: dim ker(H_u − λ) ≤ 2 at every N and every u, since an eigenvector is determined by
(ψ_0, ψ_1).

Gates C1 and C2 read both degeneracy sets off the same solved object, the real root set of
disc_x(χ(H(t))), at N = 3..14, so the empty answer and the nonempty one are not different code
paths; each is the other's control. B0 is the tie the rest of the block needs: the 2×2 is built
from the Chebyshev ansatz and never touches H_u, so a slip in either boundary row would go
unseen; det(the 2×2) = χ(H_u), exactly and symbolically, pins both rows, the ansatz and the
kernel isomorphism against the real matrix at once. B1 checks the displayed identity, B1b the
empty branch, B3 the converse in the boundary system itself, by exact reduction modulo U_{N−1},
which is also the only gate that reads entry (0,0); and B2 checks the same converse from the
outside as an exact multiplicity count in the discriminant, 2⌊(N−1)/2⌋ at u = +1 and
2⌈(N−1)/2⌉ at u = −1, the expectations combinatorial in N rather than read off the measurement.
B2 gates attainment only; the *exactly* in C_u = {±1} is C2's. C3 and C3b read the Remark on the
ring at N = 5..11 against a generic u = 1/3 on the same code path.

## (e) Parity, which is the only place N's parity enters

Let Σ = diag((−1)^l). Then ΣAΣ = −A, ΣDΣ = D and ΣPΣ = (−1)^{N−1}P, so

    Σ H_Δ(t) Σ = −H_Δ(−t)     at every N
    Σ H_u(t) Σ = −H_u(−t)     at odd N only

which is the companion gate's C5b and C5 in that order, the crack's being C5 and holding at odd N
only; the second identity is what that page calls the staggering identity. If Hψ = λψ then
H(−t)(Σψ) = −λ(Σψ), and Σ preserves |ψ_j|, so blindness at seat j at t implies blindness at seat j
at −t: **every Δ-locus is negation-closed, at every N.** The arc `the_forced_and_the_met` already
records this as the reason Δ = −1 was not worth probing.

The sector-level statement is the one this file adds, and it is a single commutation:

    Σ R = (−1)^{N−1} R Σ

> **Lemma 5.** Read on H_Δ, whose staggering identity holds at every N: at odd N, Σ preserves
> each reflection sector, so E(j) and O(j) are EACH negation-closed; at even N, Σ swaps the
> sectors, so O(j) = −E(j).

*Proof.* From the commutation applied to the map ψ ↦ Σψ of Lemma 2b's sector eigenvectors,
together with Σ H_Δ(t) Σ = −H_Δ(−t). Naming the family is not decoration: the crack's own
staggering identity holds at ODD N only, and read on H_u at even N the second clause is false,
measured at 8 of the 24 seats of N = 6, 8, 10, where O_u(j) = E_u(j) instead. §(f) applies the
lemma to H_Δ alone, so the Theorem is untouched. ∎

The same preserve-versus-swap dichotomy, by the same parity of N and between the same kinds of
operator, is already live in the repo on a different object: see
`DiabolicReflectionParityWitness` in the sweep above.

Gates P1 (odd N ∈ {5, 7, 9, 11, 13}), P2 (even N ∈ {6, …, 16}) and P3, which is P1's fence: at
even N the sectors are NOT separately closed, at 30 of the seats tested, so P1 is a statement
about odd N and not a statement that happens to be checked there.

## (f) The theorem

By Lemma 1 the crack's even sector at u IS the anisotropy's even sector at u, and the crack's
odd sector at u IS the anisotropy's odd sector at −u. Eigenvector vanishing at a site is a
property of the sector alone, so

    E_u(j) = E_Δ(j)        O_u(j) = −O_Δ(j)

and with Lemmas 2b, 3 and 4,

    u-locus(j) = E_Δ(j) ∪ (−O_Δ(j)) ∪ {+1, −1}
    Δ-locus(j) = E_Δ(j) ∪ O_Δ(j)

**Odd N.** By Lemma 5, −O_Δ = O_Δ, so u-locus = E_Δ ∪ O_Δ ∪ {±1} = Δ-locus ∪ {±1}, at every
seat and every odd N. The companion page's §(c) headline is proved, and proved past the range
N = 5..17 in which it was measured.

**Even N.** By Lemma 5, O_Δ = −E_Δ, so u-locus = E_Δ ∪ {±1} and Δ-locus = E_Δ ∪ (−E_Δ). The
relation u-locus = Δ-locus ∪ {±1} is therefore equivalent to −E_Δ ⊆ E_Δ ∪ {±1}, and since {±1}
is itself negation-closed that is equivalent to E_Δ ∪ {±1} being negation-closed, that is to the
**u-locus** being negation-closed. ∎

Two things follow that the companion page could not state.

1. **Its gate C8 is a consequence, not a measured biconditional.** That page proves one direction (a
   non-closed u-locus cannot equal a negation-closed right-hand side) and has no mechanism for
   the converse. The converse is Lemma 5 at even N.
2. **The breaking seats are named.** Applying the criterion over N = 6..16 predicts 22 seats,
   and they are exactly the 22 the companion gate reports. Gate T4 checks this against that
   gate's list taken as a literal; nothing in this file derives the expectation from a
   measurement.

The ±1 in the criterion is load-bearing rather than decorative. Dropping it predicts 30 seats
instead of 22, the eight extra being (6,1), (6,4), (10,2), (10,7), (12,4), (12,7), (14,3),
(14,10), where E_Δ is not negation-closed on its own and the ring ends restore the closure. Gate
T5 asserts exactly that difference, against the eight seats as a literal. It is the same clause
under which the companion page's gate C9 finds its two exceptions at N = 12 seats 1 and 10, where
the Δ-locus contains ±1 itself.

Gates T1 (a consistency check, see the table), T2 (odd N, 45 seats), T3, T4 and T5.

## (g) The count, which the loci leave open

The Theorem compares two SETS, and F157's blind count at a point of a locus is a MULTIPLICITY.
Nothing so far forbids a seat from being blind at the same t on both axes and blind there by a
different amount, and the companion page's construction cannot decide it. The same lemmas settle
it, once the seat's own weight is split the way Lemma 2b split the kernel.

> **Lemma 6.** Let H be real symmetric and commute with R, write E and O for its restrictions to
> the two sectors and jr = min(j, N−1−j) for the seat's fold coordinate. At a seat not fixed by R,
>
>     χ(H_j) = ½·( χ(E_jr)·χ(O) + χ(E)·χ(O_jr) )
>
> and at the reflection-fixed centre seat c of an odd chain, χ(H_c) = χ(E_c)·χ(O).

*Proof.* The vectors f_p^± = (e_p ± e_{N−1−p})/√2 for p < N/2, together with e_c in the even
sector at odd N, are orthonormal and span the two sectors. A seat not fixed by R has
e_j = (f_jr^+ ± f_jr^−)/√2, the sign being + at j < N−1−j and − at j > N−1−j. Write P_λ^± for
the projector onto E_λ intersected with that sector; then
‖P_λ e_j‖² = ‖P_λ^+ e_j‖² + ‖P_λ^− e_j‖², and each term is ½·|ψ_jr|² for that sector's unit
eigenvector at λ, the coordinate taken in the f-basis, or 0 if the sector does not carry λ. That
covers an eigenvalue carried by BOTH sectors, which is the case Lemma 7 turns on. Now §7's Cramer
identity, applied to H and then inside each sector separately, gives

    χ(H_j)/χ(H) = Σ_λ ‖P_λ e_j‖²/(x − λ) = ½·[ χ(E_jr)/χ(E) + χ(O_jr)/χ(O) ]

and multiplying by χ(H) = χ(E)·χ(O) is the display. At the centre seat e_c = f_c lies wholly in
the even sector, so the odd term is absent and the ½ with it. ∎

The ½ is the whole content, and it is the counting statement of §(a)'s first convention: **a seat
the reflection does not fix sits half in each sector.** The blocks the gate computes are §(a)
convention 2's conjugates, and neither χ(E) nor χ(E_jr) is moved by that conjugation, so the
display holds for them as written.

> **Lemma 7.** Let H be real symmetric and commute with R, and let each sector block have a
> simple spectrum. Write b_E(j) for the number of eigenvalues of E whose eigenvector vanishes at
> j, and b_O(j) likewise. Then
>
>     blind(j) = b_E(j) + b_O(j) + #{ λ ∈ spec(E) ∩ spec(O) : NEITHER sector eigenvector
>                                     vanishes at j }
>
> and in particular blind(j) = b_E(j) + b_O(j) at every knob value off C.

*Proof.* blind(j) = Σ_λ (m_λ − s_λ) with s_λ ∈ {0, 1} is Corollary C of the node-lemma proof,
and m_λ − s_λ = dim(E_λ ∩ e_j^⊥) is the dimension Lemma 2a already computes; that corollary is
stated for any Hermitian H, with no simplicity assumption, so it is available at every knob value
on both families and at the ring ends included. The simplicity hypothesis holds for both
families: folding either along R gives two tridiagonal blocks with every off-diagonal entry
nonzero, which is F160's Corollary B for the crack and the same fold applied to the anisotropy
for the chain, and the knob never touches an off-diagonal. Those blocks are unreduced Jacobi
matrices up to §(a) convention 2's positive-diagonal conjugation, which moves neither χ nor which
coordinate vanishes, so Lemma J of the node-lemma proof applies and each sector spectrum is
simple at every real t. Hence m_λ ∈ {1, 2}. An eigenvalue
carried by one sector only contributes 1 − s_λ, which is 1 exactly when its eigenvector vanishes
at j: those are the ones b_E and b_O count. A shared eigenvalue has m_λ = 2, and s_λ = 0 exactly
when both sector eigenvectors vanish at j, so it contributes 1 + eo, writing e and o for the two
indicators; against the e + o it contributes to b_E + b_O that is a surplus of
1 + eo − e − o = (1 − e)(1 − o), which is 1 exactly when neither vanishes. Summing over λ gives
the display, and off C the last sum is over the empty set. ∎

**The third term never fires, and that is a theorem rather than a range.** A shared level is
missed by the seat only when BOTH its sector eigenvectors vanish there. Off the ring ends neither
family has a shared level at all, and at the ring ends the two are the cosine and the sine of one
mode about the reflection centre, whose squares sum to 1, so no site is a node of both. Hence
#{both vanish} = 0 on both families and the correction is #shared − #{even vanishes} −
#{odd vanishes}. Gate K2c reads it, over 544 readings carrying 552 shared levels, because a term
that is always zero can carry a sign error forever: flipping its sign leaves every other check in
the file green.

**The correction is not a nuisance term; it is where the two axes must part.** By Lemma 3 the
chain has no shared eigenvalue at any t, so on the Δ axis the count is additive without
qualification. By Lemma 4 the crack has them at u = ±1 and nowhere else. So the third term of
Lemma 2b, which put {±1} into every seat's u-locus, puts a correction into
the u-count at those two points and at no other, though not at every seat: it vanishes at 12 of
the 136 seat-ends over N = 5..12. Where it is 0 is not simply the centre
seat: at the reflection-fixed centre of an odd chain it vanishes because every shared level
already has its odd eigenvector vanishing there, which the centre-seat corollary below turns on,
and it also vanishes at four non-centre seat-ends over N = 5..12, (6,1), (6,4), (10,2) and (10,7)
at u = +1, where each shared level has exactly one of its two sector eigenvectors vanishing. Those
four are exactly the seat-ends gate K4c reports as the ones where the even-N display accidentally
survives, and gate K2b pins the count they belong to.

> **Lemma 5′.** At odd N, b_E(j; t) = b_E(j; −t) and b_O(j; t) = b_O(j; −t). At even N,
> b_O(j; −t) = b_E(j; t) and b_O(j; t) = b_E(j; −t).

*Proof.* Lemma 5's map is ψ ↦ Σψ. It is a linear bijection from the eigenvectors of H_Δ(t) onto
those of H_Δ(−t), at the negated eigenvalue; it is diagonal, so it preserves vanishing at j; and
by ΣR = (−1)^{N−1}RΣ it preserves each sector at odd N and swaps them at even N. A bijection
carries the number and not only the emptiness, which is all this adds to Lemma 5. ∎

> **Theorem (count).** At every seat j and every real t,
>
> - at **odd N, and |t| ≠ 1**: blind(j; H_u(t)) = blind(j; H_Δ(t)). The two axes carry the same
>   points AND the same count.
> - at **even N, and |t| ≠ 1**: blind(j; H_u(t)) = 2·b_E(j; t) while
>   blind(j; H_Δ(t)) = b_E(j; t) + b_E(j; −t), so the two agree exactly where the even sector
>   pays alike at t and at −t.
> - at **t = ±1**, at either parity: the crack pays ⌊(N−1)/2⌋ at u = +1 and ⌈(N−1)/2⌉ at
>   u = −1, at EVERY seat, while the chain pays what P_j says there.
>
> The fence on the first two bullets is not a formality. At even N and t = ±1 the crack pays
> Lemma 4's count at every seat while 2·b_E(j; t) is 0 at all but a handful, so the even-N
> display is false there and not merely unproved; gate K4c reads both sides at both ring ends
> and finds them equal at 4 of 72 even-N seat-ends.

*Proof.* By Lemma 1 the crack's even sector at t is the chain's even sector at t, and the crack's
odd sector at t is the chain's odd sector at −t; vanishing at a site is a property of the sector
alone, so b_{E_u}(j; t) = b_{E_Δ}(j; t) and b_{O_u}(j; t) = b_{O_Δ}(j; −t). By Lemma 3, C_Δ = ∅,
so Lemma 7 on the chain reads blind_Δ = b_{E_Δ}(t) + b_{O_Δ}(t) at every t. By Lemma 4,
C_u = {±1}, so off those two points Lemma 7 on the crack reads
blind_u = b_{E_Δ}(t) + b_{O_Δ}(−t).

At odd N Lemma 5′ gives b_{O_Δ}(−t) = b_{O_Δ}(t) and the two expressions coincide. At even N it
gives b_{O_Δ}(−t) = b_{E_Δ}(t) and b_{O_Δ}(t) = b_{E_Δ}(−t), which are the two displays.

At t = ±1 the wrap bond is ±1, that is a periodic or an antiperiodic boundary condition, and the
levels are 2cos(θ_m) with θ_m = 2πm/N at u = +1 and θ_m = (2m+1)π/N at u = −1. Both readings are
F160's road at its two ends, and gate K5c decides them by an exact integer identity rather than by
a trigonometric one or a simplifier: χ(H_u(±1)) = 2·T_N(x/2) ∓ 2 at N = 3..14, T_N the Chebyshev
polynomial of the first kind, whose level sets T_N = ±1 are exactly the two combs.

R here is l ↦ N − 1 − l, so the modes adapted to it are the ones centred on the reflection
centre c = (N−1)/2, namely cos(θ_m(l − c)) and sin(θ_m(l − c)); the plain cos(2πml/N) is not an
R-eigenvector and naming it would prove nothing. The reflection puts one member of each pair in
each sector, so every degenerate level is shared, and the levels the pairing misses carry a
nowhere vanishing eigenvector: at u = +1 those are m = 0 and, at even N, m = N/2; at u = −1 there
are none at even N and one at odd N, the alternating vector at λ = −2. Hence b_E + b_O counts only
shared levels. And no site is a node of both the cosine and the sine of one mode, since their
squares sum to 1, so no shared level has both sector eigenvectors vanishing. Lemma 7 then gives
blind_u = #{shared levels}, which is Lemma 4's own count, ⌊(N−1)/2⌋ at u = +1 and ⌈(N−1)/2⌉ at
u = −1, and it does not depend on the seat. ∎

**So the answer is neither yes nor no, and the parity that split the loci splits the counts the
same way.** At odd N the crack pays the chain's count exactly, at every seat and every knob value
except the two ring ends; the shared points really are shared readings. At even N the crack pays
the EVEN sector twice where the chain pays the two sectors once each, and the two agree only when
that sector happens to pay alike at t and at −t.

At even N that condition is not independent of the locus criterion, and one direction is a
theorem. Should a seat's two loci differ, the symmetric difference off the ring ends is (−E) \ E,
so a point s of it has b_E(s) = 0 and b_E(−s) ≥ 1: the crack pays 2·b_E(s) = 0 there while the
chain pays b_E(s) + b_E(−s) ≥ 1, and the counts part at that seat too. The converse is not proved,
and it is worth saying which way the evidence runs: over N = 6..14 there is no seat at all
carrying the same locus and a different count, so within that range the counts part at exactly the
seats where the loci break and nowhere else, which is 14 of the committed 22 (the other 8 lie at N
= 16, outside the range gate K4d reads). K4d reads both halves.

At the two ring ends the crack is blind at every seat at once, and there the two axes part almost
everywhere: they agree at the reflection-fixed centre seat of an odd chain, where both pay (N−1)/2
by the corollary below. That half is proved. That they agree at NO other seat is measured, by gate
K5b over N = 5..11 at both ends, and is not proved here. In particular the agreement fails at
seats whose Δ-locus does not reach ±1 at all, where the chain simply pays nothing and the crack
pays Lemma 4's count.

> **Corollary (the centre seat).** At odd N the reflection-fixed centre seat pays exactly
> (N−1)/2, at every real t, on both axes.

*Proof.* Every R-odd eigenvector vanishes at the fixed site, so b_O is the whole odd sector,
(N−1)/2. In the folded even block the centre is the last coordinate, an END of that Jacobi chain,
and no eigenvector of an unreduced Jacobi matrix vanishes at an end; so b_E = 0. That step is not
this file's: it is (J1) of `PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA`, *"no eigenvector of H has two
consecutive zero entries; in particular no eigenvector vanishes at the first or last site"*, and
`THE_SEAT_THAT_CUTS` states it in as many words and calls it a theorem rather than a measurement.
Off C Lemma 7 then gives (N−1)/2 directly, and at u = ±1 every shared level already has its odd
eigenvector vanishing at the seat, so the correction term is 0 there too. ∎

That closes the companion page's gate B2 from the other side. That page could establish
blind(centre) ≥ (N−1)/2 by the reflection and had to MEASURE the reverse inequality at six
rational couplings; the end-of-a-Jacobi-chain argument gives it at every coupling, rational or
not, and on both axes.

> **Corollary (the multiplicity is in the resultant, doubled).** At a seat that R does not fix,
> and off C, the order of vanishing of Res_x(χ(H), χ(H_j)) at a locus point is exactly
> 2·blind(j).

*Proof.* Lemma 6 and the multiplicativity of the resultant give

    Res_x(χ(H), χ(H_j)) = 2^{−N}·(−1)^{dE·dO}·R_E·R_O·Res_x(χ(E), χ(O))²

with dE, dO the sector dimensions and R_S = Res_x(χ(S), χ(S_jr)) the sector resultants: each
factor is ∏ over the roots of one sector's characteristic polynomial of Lemma 6's right-hand
side, at which one of the two summands vanishes. (The seat must not be R-fixed: at the centre
seat of an odd chain both sides are identically zero, by §(a)'s third convention, and an order
of vanishing is not defined.)

Inside a sector the seat DISCONNECTS the tridiagonal block into halves L and R′, so
χ(S_jr) = χ(L)·χ(R′) and, by the node lemma, a blind λ is a common root of the two halves. Write
Q_S = Res_x(χ(L), χ(R′)) for the sector's halves-resultant; then R_S is a t-free rational
constant times Q_S², gate K1c. **The two Q_S multiply to F157's own polynomial**: F157's
definition route is the resultant of the two halves the seat cuts the FULL chain into, and on the
chain that equals ±Q_E·Q_O, gate K1b, which pins the ratio to ±1 and pins the split, 22 seats
carrying −1 and 10 carrying +1 over N = 4..10. WHICH sign is §(i)'s Corollary 10a, once an
argument order is named, and that split is a reading of TWO orders at once, sympy's and this
gate's own, which §(i) separates. The
same identity does not survive on the ring. What breaks it is not that striking leaves one path,
though it does: for every interior seat of the ring BOTH principal submatrices omit one end of
the wrap bond, so both are knob-free while the two sectors' halves still carry the knob, and the
ratio stops being a constant. It fails at 20 of the same 32 seats and survives at the other 12,
so gate K1b2 pins that count rather than asserting a universal.

The order of vanishing is the number of coinciding root pairs only under transversality, and here
it is: the fold puts t in one diagonal entry alone, so χ(R′) is t-free, and dα_i/dt = ±|v_i(0)|²,
the sign being − in the crack's odd sector where the fold puts −t at coordinate 0; either way it
is nonzero by (J1), because coordinate 0 is an end of the unreduced block. Off C the cross-sector
factor does not vanish, so the order at a locus point is the order in R_E·R_O, which by Lemma 7
is twice the count. ∎

On the chain the two displays compose into the statement four passages of this repository make and
none derives: Res_x(χ(H), χ(H_j)) is a t-free rational constant times F157's P_j SQUARED. What
makes the constant a constant is that Res_x(χ(E), χ(O)) is itself t-free on the chain, and Lemma 3
is not enough for that: Lemma 3 is a statement about real t, while a polynomial in t with no REAL
roots need not be constant. The complex statement is one line and does not need Lemma 3. H_Δ(t) is
tridiagonal with every off-diagonal entry 1 at every complex t, and the recursion of (J2)
determines an eigenvector from its first entry, so every eigenvalue has geometric multiplicity at
most 1 over ℂ; a value of t at which the two sectors shared an eigenvalue would give H a
two-dimensional eigenspace. So they share none at any complex t, the resultant has no roots at
all, and it is a constant, which gate K1e pins by VALUE as ±2^⌊N/2⌋ at N = 3..14. The ring has no
such argument, being no longer tridiagonal, and there the same factor carries the ring ends. Gate
K1e reads the composed statement, with K1e2 as its control. On the ring the same product carries
the cross-sector factor as well, and that factor is where the ring ends live.

**The companion page's reading of its own obstruction was one step short.** It records that the
resultant it forms is a rational multiple of F157's P_j squared and concludes that the
construction cannot see the multiplicity; what the display above says is WHERE the squaring
happens, namely inside each sector and because striking disconnects a Jacobi block, and that the
count is therefore present DOUBLED rather than absent. Halving is exact off C. What `_squarefree`
discards there is recoverable, and at the ring ends it is mostly recoverable too, which is worth
saying because the fence invites the opposite guess: gate K6b finds halving exact at 76 of the 84
ring-end seat-ends over N = 5..10, failing at exactly the eight seat-ends the gate names as a
literal, which are exactly the ones where b_E + b_O > 0, that is where some shared level has a
sector eigenvector vanishing at the seat. The cross-sector factor vanishes at all 84, so it is
not what breaks the halving; the correction term of Lemma 7 is.

Gates K0 to K7, twenty-eight checks: the split and its two controls, the resultant factorization
with the bridge to F157's P_j and its two next-nearest-entry controls, the field arithmetic
against two independent oracles, the decomposition and its correction term, the odd-N theorem
with two controls, the even-N law with its fence and its locus-versus-count reading, the ring
ends, the doubling, and the centre seat. The table in §(j) says what would redden each.

## (h) Which k sits in which sector, and the index it puts on b_E and b_O

§(g) writes the blind count as b_E + b_O and leaves the two summands as things a sector block is
asked for. F157 writes the same count as a number of node indices k. The two readings meet in one
line. The line itself is not new: it is this repository's own reflection-parity law, read on a
different comb. What IS new is that it survives the reading. The six places that own the law all
read an eigenvector of the UNIFORM chain, and the matrix here has both ends detuned; that the
blind mode of an end-detuned chain is still a SINGLE sine, on the seat's own node modulus, is
what this section establishes and what puts F157's index on §(g)'s two summands.

**What the sweep returned, and it moved this section from a discovery to a transport.** The sign
law is owned in at least six places, none of them about F157:
[PROOF_COLLISION_GAP_ODD_ORDERS](PROOF_COLLISION_GAP_ODD_ORDERS.md) §(b) states it as *"which
sector one is in is decided by the parity of k"*;
[SLOW_MODE_R_PARITY](../../experiments/SLOW_MODE_R_PARITY.md) derives it, *"b_k is R-even when k
is odd, R-odd when k is even"*; [PROOF_K_PARTNERSHIP](PROOF_K_PARTNERSHIP.md) writes
Rψ_k = (−1)^{k+1}ψ_k and credits F71; [HANDSHAKE_GEOMETRY](../../hypotheses/HANDSHAKE_GEOMETRY.md) carries it with a
carrier-parity generalization (−1)^{k−c}; `Diagnostics/Ptf/DefectReadingEquivarianceClaim.cs`
types it Tier 1; and `compute/MirrorWorld/Formulas.cs` holds it as a member,
`F71_ReflectionParity(k)`, which is the bare parity with no modulus in it at all. Every one of
the other five reads the UNIFORM chain's own Dirichlet comb, whose modulus is N+1. What
returned nothing was the assignment for F157's k, and it said so in writing in three
places: [ANALYTICAL_FORMULAS](../ANALYTICAL_FORMULAS.md) F157, the arc
`the_forced_and_the_met` in `OpenArcsRegistry.cs`, and [The Blind Seat on the
Road](../../experiments/THE_BLIND_SEAT_ON_THE_ROAD.md), the last adding *"A parity of k is the
obvious guess and is not checked here"*. **This paragraph records what the sweep found, not what
those three now say**: the same change that added this section rewrote all three to carry the
answer, so the sentence quoted above no longer exists in the tree. `fw.Confirmations` (all 24
entries), `docs/GLOSSARY.md` and
[PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA](PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md) returned nothing
on the sector ASSIGNMENT; the last of those does read the count one sector down, which is what
§(g) uses. The gap is an index change and not a labelling oversight: F157's k rides the
seat-dependent node modulus N_node = |N − 1 − 2j|, which is at most N − 1 and so never the chain's
own N+1. That is arithmetic and no gate certifies it; what the gates read is everything downstream
of it. `docs/CAUGHT_ERRORS.md` returned no entry on the assignment and two standing traps that sit
on this path, and both are fenced below: the Ad_R grading read as Hilbert-space parity, and the
vacuous anti-vacuity check that the 2026-09-02 entry logs for block K and that this block shipped
twice before three rounds took it out.

**One convention question, answered here so no reader has to ask it twice.** §(a) warns that this
file's book (hop 1, a bare t) differs from F157's (hop 2, 2Δ per end, a Δ(N−5)·I shift) by a
common scaling and a common shift, so a reader is right to look for a conversion. There is none to
make. Subtracting the shift and dividing by 2 turns F157's half-block into this file's exactly, and
Δ_k is a RATIO of two sines, so neither operation moves it: F157's Δ_k is a value of this file's t
unchanged, and every display below uses it that way.

> **Lemma 8.** Write n = N_node = |N − 1 − 2j|, θ_k = kπ/n and, with F157,
> Δ_k = sin((j+1)θ_k)/sin(jθ_k). Let n ≥ 2 and k ∈ 1..n−1 be a non-pole index, that is
> sin(jθ_k) ≠ 0. Then
>
>     ψ_l = sin((j − l)·θ_k),      l = 0 .. N−1
>
> is an eigenvector of A + Δ_k·D at the eigenvalue 2cos θ_k, it vanishes at seat j and is not
> the zero vector, and
>
>     ψ_{N−1−l} = (−1)^{k+1}·ψ_l.
>
> So the blind mode at node index k lies in the R-EVEN sector when k is odd and in the R-ODD
> sector when k is even.

*Proof.* Three claims, and the first is F157's own derivation read forward rather than backward.

The interior rows are the Chebyshev recursion sin((j−l+1)θ) + sin((j−l−1)θ) = 2cos θ·sin((j−l)θ),
which holds at every angle, so rows 1 to N−2 cost nothing (and at j = 0 or j = N−1 there is no
interior row j, but there sin(jθ_k) = 0 for every k, so every index is a pole and the lemma
says nothing). Row 0 reads Δ·ψ_0 + ψ_1 = 2cos θ·ψ_0,
which collapses to Δ·sin(jθ) = sin((j+1)θ); with ψ_{N−1} = −sin(qθ) and ψ_{N−2} = −sin((q−1)θ)
for q = N−1−j, row N−1 collapses the same way to Δ·sin(qθ) = sin((q+1)θ). Those are exactly the
two half-block conditions F157 derives separately and then matches, so both end rows hold at once
at F157's own θ_k, and by F157's own branch analysis of cos((m+1)θ) = cos((1−m)θ) only there; the
common value of the two ratios is Δ_k. Row j is one of the
interior rows and is carried by the same recursion, which is the whole of the second claim: the
two half-solutions F157 matches are ONE sine, and the relative sign the matching needs is the
oddness of the sine rather than a second constant. ψ_j = sin 0 = 0, and ψ_0 = sin(jθ_k) ≠ 0 off the pole, so the
vector is blind at j and not the zero vector.

For the parity put m = j − q = 2j − N + 1, so |m| = n and mθ_k = ±kπ. Then

    ψ_{N−1−l} = sin((m − (j − l))·θ_k) = sin(mθ_k)·cos((j−l)θ_k) − cos(mθ_k)·sin((j−l)θ_k)

and sin(mθ_k) = sin(±kπ) = 0 while cos(mθ_k) = cos(kπ) = (−1)^k, which leaves
ψ_{N−1−l} = −(−1)^k·ψ_l. ∎

The mechanism is sin(kπ − φ) = (−1)^{k+1}·sin φ, which is the same identity that gives the law on
the chain's own comb; only the angle that closes at kπ has changed, from (N+1)θ to N_node·θ. That
is why the sign comes out the same on a modulus that is never the same. It is not quite the whole
of the step, and the section's own word TRANSPORT names only half of it: every one of the six
places above reads an eigenvector of the UNIFORM chain, while the matrix here is A + Δ_k·D, whose
two ends are detuned. What is not obvious in advance is that the blind eigenvector of an
end-detuned chain is still a SINGLE sine, and that is what the first two claims of the proof
establish and what gate S1 reads against the family object.

> **Corollary 8a.** At every seat with N_node ≥ 2 and every real Δ,
>
>     b_E(j; Δ) = #{ k odd,  1 ≤ k ≤ N_node−1, non-pole : Δ_k = Δ }
>     b_O(j; Δ) = #{ k even, 1 ≤ k ≤ N_node−1, non-pole : Δ_k = Δ }

*Proof.* It is a squeeze, and the direction Lemma 8 gives on its own is only one of the two.
Distinct k in 1..N_node−1 give distinct θ_k in (0, π), hence distinct eigenvalues 2cos θ_k, which
is F157's own multiplicity clause; each sector block has a simple spectrum, which Lemma 7 ASSUMES
and its proof discharges for both families, so distinct k contribute distinct levels of the
sector Lemma 8 assigns them to. That gives

    b_E(j; Δ) ≥ #{k odd : Δ_k = Δ}     and     b_O(j; Δ) ≥ #{k even : Δ_k = Δ}.

The reverse is not Lemma 8's to give: nothing in it says a blind level of a sector has to be one
of the 2cos θ_k. What says so is that the two sides have equal SUMS. By Lemma 3 the chain's
degeneracy set is empty, so Lemma 7 reads blind(j; Δ) = b_E + b_O; and F157's locus is exact, the
seat blind at Δ exactly at the Δ_k, with blind(Δ) = #{k : Δ_k = Δ}. So the two sums agree, and
two inequalities of the same sign under one equality of sums cannot be strict. ∎

**Remark, and it is Lemma 5 arrived at from the index.** It sits after Corollary 8a and not
before it, because the clauses below are about the SETS E(j) and O(j) of §(a), and it is 8a and
not Lemma 8 that identifies those with the odd-k and even-k Δ values; Lemma 8 alone would give
inclusions and not closure. Applying sin(aπ − x) = (−1)^{a+1}·sin x
to the numerator and the denominator of Δ_{N_node−k} leaves exactly one sign, so

    Δ_{N_node−k} = −Δ_k,

and k is a pole exactly when N_node − k is, both conditions being N_node | jk. Now
N_node ≡ N−1 (mod 2). At ODD N the modulus is EVEN, so k and N_node−k share a parity and, by
Lemma 8 with Corollary 8a, a sector: each sector's locus is negation-closed. At EVEN N the
modulus is ODD, so they differ in parity and swap sectors: O = −E. Those are Lemma 5's two
clauses at every seat with N_node ≥ 2, reached from the index rather than from the commutation
ΣR = (−1)^{N−1}RΣ, and there the two routes agree. They do not agree everywhere and the index
route is the weaker: at the odd-N centre seat N_node = 0 and there is no index at all, while
Lemma 5 still speaks. §(a)'s third convention has the sector the other way round there: O(j) is
all of ℝ, every R-odd eigenvector vanishing at the fixed site, and E(j) is empty by (J1). Neither
is anything the index could have said. Gate S2c reads the two identities the argument rests on,
Δ_{N_node−k} = −Δ_k and the shared pole condition; what follows them, that k and N_node−k share a parity exactly when N
is odd, is an identity of the integers and is not gated.

> **Corollary 8b, the even-N count in closed form.** At even N and |t| ≠ 1,
>
>     blind(j; H_u(t)) = 2·#{ k odd, 1 ≤ k ≤ N_node−1, non-pole : Δ_k = t }

*Proof.* The Theorem (count)'s second bullet reads blind(j; H_u(t)) = 2·b_E(j; t); Corollary 8a
supplies b_E where N_node ≥ 2, and by Lemma 1 the crack's even sector at t IS the chain's, so the
Theorem's b_E and Corollary 8a's are one number. At even N the node modulus |N − 1 − 2j| is ODD,
so it is never 0 and the only case left out is N_node = 1, where the index range 1..N_node−1 is
EMPTY: F157's locus is empty there and both sides are 0. That last case is argued and not gated,
S4 skipping the seats it covers. ∎

That is what open item 1 of the companion page asked for. The odd-N side needs nothing new: there
the crack pays what the chain pays, and away from the reflection-fixed centre seat the chain's
count is F157's own #{k}. AT that seat N_node = 0, so #{k} = 0 while the seat is in fact blind at
every coupling by (N−1)/2; §(g)'s centre-seat corollary covers it and §(h) does not speak there.

**Three fences, and each is a confusion this repository has already paid for once.**

- **Not the chiral K.** The partner map ψ_k ↦ ψ_{N+1−k} is a different operator on the same
  modes, and [PROOF_K_PARTNERSHIP](PROOF_K_PARTNERSHIP.md) keeps the two apart. The node index
  has its own involution k ↦ N_node − k, and the Remark above says what that one does; it is not
  the chiral partner map and it does not act on the same modulus.
- **Not Ad_R.** The reflection's grading on OPERATOR space is not this Hilbert-space parity.
  [The Blind Site](../../experiments/THE_BLIND_SITE.md) records that confusion as caught twice,
  and every ψ here is a single-excitation vector.
- **Not the superseded reflection-parity reading of the blind LAW.**
  `Core/Symmetry/SeatCutBlindnessClaim.cs` records a *"reflection-parity reading"* as a rival to
  the divisor law that was superseded, and warns that it agrees with the divisor law at every
  seat of a prime chain, so N = 5 and N = 11 cannot tell them apart. That reading is about the
  SEAT index; Lemma 8 is about the MODE index, and it is a statement about which sector a mode
  lies in rather than about whether a seat is blind. The checks that read a seat against a sector
  are S3 and S4, whose ranges are N = 5..14 and even N = 6..14, so the composites that
  discriminate, 6, 8, 9, 10, 12 and 14, are swept rather than assumed away. S1, S2, S2b, S2c and
  S2d reach N = 20, but they read modes and never a seat count, so the wider range is not what carries
  this fence.

Gates S0 to S4, eight checks. S1 and S2b read the mode and its parity at 812 (N, seat, k) modes
over N = 5..20, exactly, as polynomial identities in the field ℚ(2cos(π/N_node)); S3 reads
Corollary 8a against the SECTOR BLOCKS, which is block K's own route through a gcd over
ℚ[t]/(μ), μ the Δ value's minimal polynomial, at 196 (seat, Δ) readings over N = 5..14. That is
the point of S3: Corollary 8a's right-hand side counts k, while the left-hand side is computed by
a route that never mentions a mode or a k. S2 reads the two ingredients of the parity apart, since
the parity is their product; S2d is the block's one CONTROL, asserting that the same construction
on the modulus N_node+1 carries no parity at all, so no mutation of the objects tests it and what
does is feeding it the right modulus; S4 reads Corollary 8b against `blind_at` on the full crack
matrix. The exactness discipline is block K's: the node angle never reaches a simplifier, every
value lives in an integer minimal polynomial's quotient, and the field itself, the one every
reduction happens modulo, is verified by S0 rather than trusted.

**Block S was mutated nine times**, by hand and with no committed artifact: the parity exponent
k+1 to k; the odd and even roles swapped in the k split; the sweep narrowed from N = 5..20 to
5..13; the node modulus |N−1−2j| to |N−2j|; the anisotropy in `H_aniso` moved from (0, N−1) to
(0, N−2); `min_poly_of` returning a non-minimal annihilator; S2d's control fed the RIGHT modulus;
the node field replaced by the one at N_node+1; and S3's assertion that its eight balanced
readings sit at Δ = 0 retargeted to Δ = 1, which reddens S3, since a clause that only ever agrees
is the defect this block already shipped three times. Every check reddens under at least one. S2d
is the block's only control, so a mutation of the objects is not what tests it; the seventh,
feeding it the right modulus, is, and it reddens S2d alone. Three of the nine are there because an
earlier build of this block was GREEN under them: the third, when four checks still had a bare
sweep range; the fifth, when S1 wrote its own rows instead of reading `H_aniso`; and the seventh,
which tests a control that exists only because the dead fence it replaced could not fire. That
earlier build also had a ninth check asserting that the SWAPPED sector assignment is refuted,
which is a strict consequence of S3 and could never redden alone; it had a tenth reading arithmetic
about |N−1−2j| as though it certified the separation of two combs; and S1 wrote its own
tridiagonal rows instead of reading `H_aniso`, so moving the anisotropy left it green. The first
two are gone, the third is the fifth mutation above. THREE dead clauses went with them, ψ_j = 0
and ψ ≠ 0 inside S1 and "the opposite sign at none" inside S2b: the first two are true of any
input by how the mode is built, and the third is excluded by the pole test, so all three were
red-makers no input could reach. The swap clause took two passes to remove rather than one. It
began as a check of its own, was folded into S3's predicate when that check turned out to be S3's
own consequence, and was still inert there for the same reason; what refutes the swap is S3's
equality together with the pinned 188 readings that have #odd ≠ #even, and that is where it now
lives.

**What this opened, and §(i) is that step.** §(g)'s Corollary factors F157's own definition
polynomial as P_j = ±Q_E·Q_O, the two sectors' halves-resultants, and says in as many words that
WHICH sign is not a law it identifies. Corollary 8a names the roots of each factor: Q_E should be
the odd-k product and Q_O the even-k product, up to a rational constant. That is not checked here
and it is not a corollary of anything above, since a factorization of a polynomial into two
factors is not determined by the root sets alone without an argument about multiplicities.
§(i) supplies both: the constant in closed form, and the multiplicity by reading Lemma 8 on the
middle route, a uniform path, rather than on the end-detuned chain.

## (i) The two constants, and the convention the sign lives in

§(g) and §(h) leave two constants unidentified between them. This section closes both, and they
turn out to be one object: a single congruence produces the product, and what is left over is a
resultant that evaluates to a sign.

**Notation, since this section computes rather than counts.** Write S_m for sin(mθ)/sin θ as a
monic integer polynomial in x = 2cos θ, which is `S_cheb` in the gate and the S of §(h)'s block;
in F157's letter S_m(x) = U_{m−1}(x/2), so F157's own P_j(Δ) = Res_x(U_{N_node−1}, Δ·U_{j−1} − U_j)
reads here as Res(S_n, Δ·S_j − S_{j+1}). Put α_p = S_{p+1} − t·S_p, the characteristic polynomial
of the p-site path carrying the knob at coordinate 0. Throughout, p = jr = min(j, N−1−j) is §(g)'s
fold coordinate and n = N_node = |N−1−2j| is §(h)'s node modulus. The knob is this file's t
throughout; §(h)'s convention paragraph settles that F157's Δ is that same t with no conversion to
make, and Corollary 11a below settles the remaining half of it, that F157's SEAT-indexed Δ_k is
the FOLD-indexed one this section computes with. Q_E and Q_O are §(g)'s: Q_S = Res(χ(L), χ(R′))
for the two halves the struck seat cuts sector S's block into. **Every resultant that carries a
LAW here is the Sylvester determinant lc(f)^{deg g}·∏g(α), built in the gate as `res_std`;
sympy's own routine appears in block W only where it is the object under test.** That distinction
is the section's subject, not its housekeeping: a resultant is antisymmetric up to
(−1)^{deg f·deg g}, so a sign law has to name an argument order, and `sympy.resultant` does not
keep the one it is given: it agrees with the Sylvester determinant when deg f ≥ deg g and is
(−1)^{deg f·deg g} times it when deg f < deg g (gates W0 to W0c). At equal degrees it is
order-dependent, and correctly so, gate W0d, which is why the rule is stated by degree and never
by symmetry.

**What the repo held, store by store.** [`docs/ANALYTICAL_FORMULAS.md`](../ANALYTICAL_FORMULAS.md):
F157 owns P_j as a DEFINITION, together with Δ_k, the pole rule and the multiplicity clause; it
states no leading coefficient, no sign and no factorisation, and its four worked rows are printed
primitive with a positive leading coefficient, which is the normalisation that discards exactly
the number wanted here. F160 owns the ring's quantisation curve in the same U_n(x/2)
normalisation. F139 owns a polynomial it also writes S_m, in the same 2cos normalisation but
indexed one higher (F139's S_m(2cosθ) = sin((m+1)θ)/sinθ), and owns the product-to-sum rule, not
the addition formula. [`docs/proofs/`](.):
[PROOF_CRACKED_RING_EXACT_CURVE](PROOF_CRACKED_RING_EXACT_CURVE.md) §(b) owns the Cassini identity
S_{p+1}² − S_p·S_{p+2} = 1, as its own gate P2, and until this change owned it alone; Lemma 9
below uses it and cites it rather than repeating it. The Chebyshev ADDITION formula
S_{a+b} = S_a·S_{b+1} − S_{a−1}·S_b is in no proof of this repository under any spelling, so
Lemma 9's other half is written out here. Five other files in `docs/proofs/` form a resultant:
`PROOF_CODIM1_BY_ADDITIVITY`, `PROOF_EDGE_BLOCK_DEFECTIVE_UNDER_PROFILE` and
`PROOF_F103_F87_Z2_CUBED_REFINEMENT` are F89 fold certificates and divisibility arguments, and
`PROOF_F127_RESIDUE_COLLAPSE` and `PROOF_F87_WINDOWED_MONOMIAL_CONVERSE` are resultant-and-Sturm
questions on a different parameter; none states a closed form for a resultant.
[`experiments/`](../../experiments/), null results included: nothing on either constant;
[The Blind Seat on the Road](../../experiments/THE_BLIND_SEAT_ON_THE_ROAD.md) meets the rational
multiple, records that its `_squarefree` discards it, and writes *"No gate here reads it."*
[`OpenArcsRegistry.cs`](../../compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs): the arc
`the_forced_and_the_met` carried both as open and as two DIFFERENT indeterminacies, which was the
right reading of what was then known and is the reading this section overturns: one congruence
answers both.
`fw.Confirmations`: nothing, and nothing here is proposed for a flight.
[`docs/GLOSSARY.md`](../GLOSSARY.md): nothing, as §(a) already records.
[`docs/CAUGHT_ERRORS.md`](../CAUGHT_ERRORS.md): no entry on the unresolved sign, but the
2026-09-02 entry records that "exactly −1" was written once and repaired, so a constant sign is a
shape this page has already shipped. `compute/`: `Cyclotomy.cs` holds integers and no polynomial
algebra, `Crack.cs` holds a Chebyshev routine and no resultant, and `BlindSeat.cs` carries the
N = 9 seat 1 locus un-normalised as 128·Δ(Δ−1)(Δ+1)(Δ²−3). That 128 is in F157's hop-2 book, not
this file's, and this section does not identify it.

**Lemma 9 (the congruence).** For p ≥ 1 and n ≥ 1, S_p·α_{N−1−p} ≡ −S_n (mod α_p).

*Proof.* α_p is monic in x, so the reduction is a division in ℤ[t][x]. Write q = N−1−p, so that
q − p = n. The addition formula at (a, b) = (n, p+1) and at (n, p) gives S_{q+1} = S_n·S_{p+2} −
S_{n−1}·S_{p+1} and S_q = S_n·S_{p+1} − S_{n−1}·S_p, hence

    α_q = S_{q+1} − t·S_q = S_n·(S_{p+2} − t·S_{p+1}) − S_{n−1}·(S_{p+1} − t·S_p),

and the second summand is S_{n−1}·α_p, which vanishes modulo α_p. Multiply the first by S_p and
use t·S_p ≡ S_{p+1}: it becomes S_n·(S_p·S_{p+2} − S_{p+1}²), and Cassini says that bracket
is −1. ∎

The congruence is the whole mechanism. S_n is the characteristic polynomial of the MIDDLE ROUTE,
the uniform stretch the seat leaves standing between its own two mirror images, and the two
sector LEFT halves are both literally α_p, because the reflection's sign only reaches the far end
of the folded block and the near end is where the knob sits. So Q_E·Q_O = Res(α_p, S_n), and the
only thing between that and F157's polynomial is the factor S_p the lemma had to introduce.

**Lemma 10 (the common factor).** For p ≥ 1, Res(α_p, S_p) = (−1)^binom(p,2), with no t in it.

*Proof.* S_p is monic of degree p−1 and α_p monic of degree p, so p(p−1) is even and the two
orders agree; take the second, over the roots of S_p. Those are y_k = 2cos(kπ/p) for
k = 1..p−1, simple and exhausting S_p. There S_p vanishes, so α_p(y_k) = S_{p+1}(y_k), and the
recursion S_{p+1} + S_{p−1} = x·S_p gives S_{p+1}(y_k) = −S_{p−1}(y_k). Finally
S_{p−1}(y_k) = sin((p−1)kπ/p)/sin(kπ/p) = (−1)^{k+1}, so the product over k is
(−1)^{p−1}·(−1)^{Σ(k+1)} = (−1)^{(p−1)(p+4)/2}. That exponent differs from binom(p,2) by 2(p−1),
which is even at every integer p. ∎

**Corollary 10a (the sign, identified).** At every interior seat the reflection does not fix,

    Res(α_p, α_{N−1−p})  =  (−1)^binom(p+1,2) · Q_E · Q_O,      p = min(j, N−1−j),

with every resultant taken FOLD HALF FIRST, in the two sectors as well as outside.

*Proof.* Lemma 9 gives S_p·α_{N−1−p} ≡ −S_n (mod α_p); α_p is monic, so Res(α_p, ·) depends on
its second argument modulo α_p only and is multiplicative in it. Hence Res(α_p, S_p)·Res(α_p,
α_{N−1−p}) = Res(α_p, −S_n) = (−1)^p·Res(α_p, S_n) = (−1)^p·Q_E·Q_O, the power being deg α_p = p,
and Lemma 10 turns the left factor into (−1)^binom(p,2). The exponent is binom(p,2) + p =
binom(p+1,2). ∎

**The sign is a function of the fold coordinate and of nothing else.** N does not appear in it,
and that is an algebraic statement about the display above, not a reading of a sweep. One caution
about what gate W3 can see: swapping the two arguments of every resultant in the block leaves W3
green, the induced factor being (−1)^{p(p−1)} = +1, so the ORIENTATION is read by W0 and by W6's
degree line and not by the sign law's own check. It is worth
saying because the reading §(g) pinned looks N-dependent: it is the same law seen through a
seat-indexed argument order.

**Corollary 10b (the factorisation, with its constant).** Let β_S be the sector characteristic
polynomial of the middle route, so that β_E·β_O = S_n. Split β_S against S_p into its pole part
g_S = gcd(β_S, S_p) and its pole-free part h_S = β_S/g_S, and put n_S = deg h_S, r_S = deg β_S.
Then

    Q_S  =  c_S · ∏ over the non-pole roots x_k of β_S of (t − Δ_k),
    c_S  =  (−1)^(p·r_S + n_S) · Res(h_S, S_p) · Res(g_S, S_{p+1}),

one factor per NON-POLE root of β_S, a pole root contributing a constant instead, so a value
shared by several appears with that multiplicity and deg_t Q_S falls below deg β_S at 20 of the
144 readings, the smallest being N = 7 seat 2 where β_E's one root is a pole and Q_E is a
constant. In
particular deg_t Q_S = n_S and lc_t Q_S = c_S, and c_S is a nonzero integer.

*Proof.* Q_S = Res(α_p, β_S) = (−1)^{p·r_S}·Res(β_S, α_p), and β_S is monic, so the second factor
is ∏ over the roots x_k of β_S of α_p(x_k). At a non-pole root α_p(x_k) = −S_p(x_k)·(t − Δ_k),
with Δ_k = S_{p+1}(x_k)/S_p(x_k); at a pole root S_p(x_k) = 0 and the factor is the constant
S_{p+1}(x_k). Collecting the n_S signs and reading the two remaining products as resultants
against h_S and g_S gives the display. Both factors of c_S are integers because g_S and h_S are
monic integer polynomials, β_S dividing the monic S_n; and c_S ≠ 0 because β_S divides S_n, which
is squarefree, so h_S is coprime to S_p, while gcd(S_p, S_{p+1}) = 1 by the recursion. ∎

**Corollary 10b is about roots; which roots they are is Lemma 8's, and the gate reads the two
apart.** The identity of Corollary 10b holds for any monic SQUAREFREE β in
place of β_S, so it carries no information about node indices at all, and gate W6, which reads it,
is blind to the parity. Squarefreeness is where the pole split has to capture the full
multiplicity; β_S has it because β_S divides S_n, and without it the identity fails, gate W6c,
which draws squarefree and non-squarefree β through the same door and pins both outcomes. What supplies the parity is §(h)'s Lemma 8 read on the middle route, which
is a UNIFORM path and hence F71's own comb rather than the end-detuned one §(h) had to transport
it to: β_E carries the odd node indices and β_O the even ones, decided in the node field and not
asserted. Two gates carry that jointly and neither alone: W5b reads the closed forms against the
node comb and W5 reads the sector blocks against the closed forms, so swapping the parity inside
the block reddens W5 only, while swapping it in the closed forms reddens both. Only
with W5b does the factorisation say what Corollary 10b's display means, and only then does it
compose with §(h): the multiplicity of a value Δ in Q_E is #{k odd, non-pole : Δ_k = Δ}, which is
Corollary 8a's b_E(j; Δ) unchanged, and the same for Q_O and b_O. That composition is a corollary
of Lemma 8 and W5b; no gate measures it against a b_E.

**The factors repeat, and simplicity would contradict F157.** At N = 10 seat 2 the odd indices
k = 1 and k = 3 give one and the same Δ, so Q_E = −(t−1)²; at N = 14 seat 3 the multiplicity is
three. F157's committed row at N = 11 seat 2, 3Δ⁴ − 4Δ², has 0 as a double root carried by k = 2
and k = 4, both even and so both landing in Q_O, which is a repeated factor of one halves-resultant
and not merely of the product. An earlier build of this section wrote "every factor simple" and
its gate read only the degree and the leading coefficient, which cannot see a repeated factor;
that is recorded in `docs/CAUGHT_ERRORS.md` for 2026-09-03.

**Lemma 11 (the node identity).** At every root x_k = 2cos(kπ/n) of S_n, for n ≥ 2 and every
p ≥ 0 under the sine convention S_{−1} = −1,

    S_{p+n}(x_k) = (−1)^k · S_p(x_k).

*Proof.* The addition formula at (a, b) = (p, n) gives S_{p+n} = S_p·S_{n+1} − S_{p−1}·S_n. At x_k
the second term vanishes and S_{n+1}(x_k) = sin(kπ + θ_k)/sin θ_k = (−1)^k. Those two readings are
exactly what gate S2 pins for §(h). ∎

**Corollary 11a.** F157's Δ_k, written with the SEAT index as sin((j+1)kπ/n)/sin(jkπ/n), is the
fold-indexed S_{p+1}(x_k)/S_p(x_k) used above, and the two agree on which k is a pole. At a seat
with j ≤ N−1−j the two are the same expression. At the mirror seats j = p + n, and Lemma 11 gives
S_{j+1}(x_k) = (−1)^k·S_{p+1}(x_k) and S_j(x_k) = (−1)^k·S_p(x_k), so the factor cancels in the
ratio and the vanishing of the denominator is unmoved. Gate W5c. ∎

**Corollary 11b (what is factorised is F157's own generator).** With G_j = Res(S_n, Δ·S_j −
S_{j+1}), which is F157's definition route in this section's letters,

    Res(α_p, α_{N−1−p})  =  (−1)^e · G_j,     e = (n−1)(p+1) + p + binom(p,2) + [j > N−1−j]·binom(n,2).

*Proof.* G_j = Res(S_n, −α_j) = (−1)^{n−1}·(−1)^{(n−1)j}·Res(α_j, S_n), the first power because
deg S_n = n−1 and the second by antisymmetry. When j = p, Corollary 10a's chain gives
Res(α_p, S_n) = (−1)^{p + binom(p,2)}·Res(α_p, α_{N−1−p}), which is the first three terms of e.
When j > N−1−j then j = p + n, and Lemma 11 gives α_j(x_k) = (−1)^k·α_p(x_k) at every root of S_n,
so Res(S_n, α_j) = ∏(−1)^k·Res(S_n, α_p) = (−1)^binom(n,2)·Res(S_n, α_p), which is the fourth
term. Two bookkeepings are free and are named rather than performed: the antisymmetry factors
(−1)^{(n−1)j} and (−1)^{(n−1)p} differ by (−1)^{(n−1)n}, which is +1, and e's first term may be
written with j or with p for the same reason. Gate W7c, at all 72 seats; W7b reads the same thing against F157's four committed rows,
where the sign is what F157's primitive-and-positive normalisation discards. ∎

**Composition.** lc_t(P_j) = (−1)^binom(p+1,2)·c_E·c_O and deg_t P_j = n_E + n_O. Gate W7 reads
it, and reddens only where W3 or W6's degree-and-leading-coefficient line does: given those two it
is a theorem rather than an independent measurement, which is what it is there for.

**What §(g)'s literal is a reading of.** The identity P_j = ±Q_E·Q_O and everything squared that
rests on it are untouched by the convention question, since a per-argument sign of ±1 cancels in
|ratio| = 1 and in K1c's and K1e's squares. What the committed literal K1B_SIGNS = {−1: 22,
+1: 10} is, on the other hand, is a reading of TWO choices at once: sympy's ordering, and K1b's own
seat order, which puts the left half of the chain first whether or not it is the fold half. Gate W4
reproduces the committed split from here and pins the three readings apart; W4b names the seats
rather than counting them, and counts the two causes apart over its own sweep: over N = 4..12 the
convention alone moves six of fifty seats and the seat order alone moves six. Over N = 4..10 the
convention flips (5,1), (7,1), (9,3) and (9,7) while leaving the TOTAL at 22:10, so the committed
literal is accidentally convention-stable in exactly the range it is pinned over, and is not so at
N = 4..12, where the same comparison reads 34:16 against 32:18. A total is what hid this; the seat
lists are what the gate pins now.

**Four fences, and each is a place an earlier build of this section was wrong.**

1. **The R-fixed centre seat is out**, because there the seat cuts the chain into two matrices
   the reflection exchanges, which therefore share a characteristic polynomial, so the outer
   resultant vanishes identically in t and there is no ratio to carry a sign. Gate W3c reads that
   at the five odd chains in range, and what such a reading can do is worth saying: the vanishing
   is forced by R-symmetry alone, so no mutation of an R-symmetric family reddens W3c, and it
   reddens only through its own population pin or through a broken resultant. It is a READ that
   decides a word, not evidence. §(a)'s third convention already covers that seat and §(g)'s centre-seat
   corollary already pays its count; Corollary 10a does not speak there. The end seats j = 0 and
   j = N−1 are out of the gate's sweep as well; there α_p is the empty determinant and both sides
   of Corollary 10a are 1, so the statement is true and says nothing.
2. **"Fold half first" is a statement about the sectors too**, not only about the outer resultant.
   Taking Q_S with the right half first breaks the law at 12 of 50 seats over N = 4..12, gate W8b.
3. **"Per non-pole root", never "simple"**, as above. The population that shows it begins at
   N = 10: over N = 4..9 no reading in this BLOCK carries a repeated factor at all, so a sweep narrowed below
   N = 10 would make the false clause unfalsifiable. Gate W6b pins the set of readings where the
   repetition occurs, and empties, and reddens, if the sweep is narrowed.
4. **What is thin at the largest fold coordinate is the gate's coverage, not the claim.**
   Corollary 10a is algebraic in p and holds at every N; nothing about it rests on the sweep. But
   over N = 4..14 the seats per p run 22, 18, 14, 10, 6, 2, so a reader checking the law by
   reading the run sees p = 6 only twice. Gate W3b pins that profile, so a narrowed sweep reddens
   rather than quietly thinning the one place the reading is thin.

**Mutation record.** Block W carries twenty-three checks, and they fall into three kinds.
**Six assert that something BREAKS**, W0b, W0c, W0d, W6c, W8 and W8b, so each is GREEN in a
healthy block and the input that reddens it is the one that would repair what it fences, like the
six §(j) already carries. **Two redden only through their own pin or a broken resultant**, W3b
whose profile no arithmetic defect can move, and W3c whose content is forced by the reflection
(fence 1). **The remaining fifteen redden under a mutation of the OBJECTS**, and those run were: the congruence's sign, both exponents of Lemma 10 and Corollary 10a,
Lemma 11's (−1)^k, the constant's n_S term, the middle route's parity swap, W5b's parity
condition, W5c's fold index, W7c's mirror-side term, the routing of the sector order, the routing
of the outer order, the coefficient blocks of the Sylvester matrix, the swap rule stated the other
way round, each pinned literal in turn, and, for the pinned checks, each literal in turn and the
sweep narrowed to N = 4..9. Three are worth naming. Swapping the middle route's parity reddens W5 and W5b and leaves W6 GREEN, which is
what shows W6 reads roots and not indices, and is why W5b exists at all. Dropping the sign in the
monic product reddens W6 on the PRODUCT line alone, leaving its degree and leading-coefficient
reading green, which is what shows the earlier check was blind rather than merely weak. Feeding
W8's control the R-SYMMETRIC family through the same door reddens W8, which is what shows the
control's door is the one the real family uses. Two attempted mutations turned out not to mutate the objects at all, one inserting a dead loop and
one turning a check into a tautology; both were reported as dead checks before they were re-read,
and both are in the ledger rather than in this list.

## (j) The gates, in one table

| gate | what would make it red |
|---|---|
| L1a | the sector reduction failing at one N or one sector; exact zero matrix, symbolic in u |
| L1b | the identity holding with the WRONG odd-sector sign, which would make L1a vacuous |
| L1c | R failing to commute with P or D, or commuting with a one-end diagonal |
| L2a | one seat where locus(j) ≠ E ∪ O ∪ C, either axis |
| L2b ×3 | a dropped summand noticed by a number of seats other than 54, 10, 14 |
| C1 | one N or one Δ at which the anisotropy's spectrum degenerates |
| C2 | the crack's degeneracy set differing from {+1, −1} at one N |
| C3 | one ring seat blind by less than the number of degenerate pairs |
| C3b | every seat blind at a generic u, which would make C3 a property of the code path |
| B0 | det(the 2×2) ≠ χ(H_u) at one N, i.e. the Chebyshev reduction not describing the real matrix |
| B1 | the boundary identity failing at one N, which would break Lemma 4's ⊆ direction |
| B1b | the U_{N−2} = 0 branch turning out to be nonempty |
| B2 | a ring end not attained, or attained with a multiplicity other than the predicted one. It gates ATTAINMENT only; the *exactly* in C_u = {±1} is C2's |
| B3 | one of the four entries failing to vanish at u = −U_{N−2}; the only gate that reads entry (0,0) |
| P1 | one odd-N seat whose sector locus is not negation-closed |
| P2 | one even-N seat where O ≠ −E |
| P3 | the even-N sectors turning out separately closed after all |
| T1 | **a consistency check throughout, not independent evidence.** Its even leg cannot fail at all: E_u = E_Δ is L1a, the two blocks being the same matrix. Its odd leg tests that the locus of a family at −t is the negation of its locus, i.e. that `real_factors` and `negate` commute with t ↦ −t |
| T2 | one odd-N seat where u-locus ≠ Δ-locus ∪ {±1} |
| T3 | one even-N seat where the two-term forms fail |
| T4 | a predicted break set differing from the committed 22 |
| T5 | the ring-ends clause turning out not to change the criterion, or changing it at seats other than the eight named |
| K0 | one seat where χ(H) is not χ(E)·χ(O), or the struck polynomial is not the sector split; exact zero polynomial, symbolic in t and x, population pinned at 98 |
| K0b | one summand of the split turning out to suffice on its own at a non-centre seat, or the even one failing to suffice at a centre seat. Two verdicts through one door |
| K0c | the one-end diagonal, which does not commute with R, satisfying the split anyway through the same door |
| K1 | the resultant factorization failing, its constant not being exactly 2^(−N), or the sign law being wrong at one N |
| K1b | on the CHAIN, F157's definition polynomial not being a constant times the product of the two SECTOR halves' resultants |
| K1b2 | that identity holding on the RING at more than the 12 seats where it survives; there both principal submatrices are knob-free while the sector halves are not, which is what breaks it |
| K1c | a sector resultant not being a t-free constant times its own halves' resultant SQUARED |
| K1d | that square surviving one next-nearest entry in the sector block, which would mean K1c reads the construction and not the disconnection |
| K1e | on the CHAIN, Res(χ(H), χ(H_j)) not being a t-free constant times F157's P_j SQUARED. **This is the claim four passages of this repository state and none derives** |
| K1e2 | that constant surviving one next-nearest entry in the chain, which would mean K1e reads the diagonal's position rather than the disconnection |
| K2a | the Q[t]/(μ) arithmetic disagreeing with the committed rational `blind_count`, or with sympy's own gcd over ℚ(α) at two quadratic and one QUARTIC μ; the quartic is what reaches the modular inverse's polynomial branch, which no quadratic does |
| K2 | one locus point off C where blind(j) ≠ b_E + b_O; population pinned at 96 |
| K2b | the correction at a ring end differing from #{shared levels neither sector eigenvector misses}, or biting at a number of seat-ends other than the pinned one |
| K2c | the both-vanish term of that correction turning out nonzero anywhere, or the loop that reads it walking past no shared level at all, which would make "always zero" a statement about an empty set |
| K3 | one odd-N locus point off the ring ends where the two axes pay different counts. **It is an equality between two readings of one routine, so a defect that moves BOTH axes alike survives it, a uniform miscount included**; K3c is the control against the other failure |
| K3b | the two axes AGREEING at u = +1 at the four seats whose two values are derived from F157's committed P_j and generator and from Lemma 4, never from this run. Its ring column repeats what K5 asserts over a wider range; its content is the four chain numbers |
| K3c | the R-breaking one-end family tracking the crack anywhere off the knob value 0, where all three families are one matrix |
| K3d | Lemma 5′ read on the counts: an odd-N sector count that is not negation-invariant, or an even-N pair that does not cross over. K3 uses this lemma and does not test it |
| K4 | one even-N locus point where the ring is not twice the even sector, or the chain is not b_E(t) + b_E(−t). **It reads `locus_mus`, which removes the ring ends**, so the fence on the even-N bullet is K4c's business and not K4's |
| K4b | the even-N counts not parting at the four literal triples, two of which have the ring paying 0 where the chain pays 1; a point that left the locus reddens it rather than crashing it |
| K4c | the `\|t\| ≠ 1` fence turning out not to be load-bearing, that is blind_u equalling 2·b_E at the ring ends at more than the four of 72 even-N seat-ends named |
| K4d | an even-N seat where the counts part and the loci do not, or a break seat where they do not part; the first half is the direction §(g) proves, the second is measured over N = 6..14 |
| K5 | one seat of the ring at u = ±1 paying other than Lemma 4's degenerate-pair count |
| K5b | the two axes agreeing at a ring end anywhere but the odd-N centre seat, or failing to agree there |
| K5c | the ring ends not carrying the periodic and the ANTIperiodic comb; an exact integer identity against 2·T_N(x/2) ∓ 2, which is what the count Theorem's third bullet rests on |
| K6 | one locus point off C where the resultant's order of vanishing is not exactly twice the count |
| K6b | halving turning out to fail at the ring ends anywhere but the eight named seat-ends, or at one of those turning out not to have b_E + b_O > 0 |
| K7 | the centre seat paying other than (N−1)/2 at one of six knob values, two of them irrational and two of them the ring ends |
| S0 | the node field Q(2cos(π/n)) failing to be irreducible, or to carry the degree φ(2n)/2, or to divide U_{n−1}. Every reduction in the block happens modulo this polynomial, so a wrong one from sympy's algebraic machinery would silence the block rather than redden it |
| S1 | one row of `H_aniso(N, Δ_k)`ψ = 2cos θ_k·ψ not vanishing in the field, or a population other than the pinned 812 modes and 412 poles. It reads the FAMILY OBJECT: an earlier build wrote the rows out by hand and stayed green when the anisotropy moved to (0, N−2). It asserts nothing about ψ_j = 0 or ψ ≠ 0, which are true of any input by construction |
| S2 | the node angle not closing: S_{N_node}(x_k) ≠ 0, or T_{N_node}(x_k) ≠ (−1)^k, or a population other than the pinned 1224 angles. The parity is the product of the two, so this reads them apart, where a compensating pair of sign errors would be invisible |
| S2b | one site where ψ_{N−1−l} ≠ (−1)^{k+1}ψ_l, or a population other than the pinned 812. An earlier build added "or the OPPOSITE sign also passing", which cannot happen: both signs at once force ψ = 0, which the pole test excludes. S2d is the live control instead |
| S2c | Δ_{N_node−k} ≠ −Δ_k, or the pole condition not being shared by k and N_node−k, or a population other than the pinned 812. Those two readings are the Remark's content; its third clause, that k and N_node−k share a parity exactly when N is odd, was dropped from the check because it is an identity of the integers |
| S2d | **a control: it asserts that something BREAKS.** The same mode built on the modulus N_node+1 carrying either reflection sign at even one of the pinned 851 constructions, or a population other than 851. So no mutation of the OBJECTS tests it; feeding it the right modulus does, and does redden it. It is what makes S2b a reading of the node angle rather than of the shape of the formula |
| S3 | one (seat, Δ) reading where the sector blocks' b_E and b_O are not the odd-k and even-k counts, or one of the 8 balanced readings sitting anywhere but Δ = 0, or a minimal polynomial that does not annihilate its element or is not irreducible, or a sweep shrunk below the pinned 196 readings and 188 asymmetric ones |
| S4 | one even-N locus point off the ring ends where the crack does not pay 2·#{k odd on that Δ}, read against `blind_at` on the full N × N crack matrix, which knows nothing about sectors or modes; or a population other than the pinned 36. A locus point carrying no node index predicts 0 and is read, not skipped |
| W0 | res_std disagreeing with the DEFINITION lc(f)^deg(g)·∏g(α) at one of 60 random pairs, or the pairs losing their non-monic and degree-0 members, pinned at 31 and 15. Without those the lc(f)^deg(g) factor is 1 throughout and both early returns of the routine go untested, though production reaches them |
| W0b | sympy.resultant agreeing with the Sylvester determinant everywhere, or differing at a population other than the pinned 5 of 60 |
| W0c | the difference being anything but the swap: one pair where sympy is not res_std times (−1)^{deg f·deg g} exactly when deg f < deg g, or a population other than the 33 pairs with a nonzero resultant of which 4 sit at EQUAL degrees, where the rule's other branch asserts plain agreement |
| W1 | one seat where S_p·α_{N−1−p} + S_{N_node} is not divisible by α_p; an exact remainder in ℤ[t][x], population pinned at 72 |
| W1b | one (n, p, k) where S_{p+n}(x_k) ≠ (−1)^k·S_p(x_k), or a population other than the pinned 1260 of which 1014 are not the trivial 0 = 0 at S_p(x_k) = 0. It is Lemma 11, and Corollaries 11a and 11b both stand on it |
| W2 | Res(α_p, S_p) carrying a t, or differing from (−1)^binom(p,2) at one p |
| W3 | one non-centre seat where the ratio is not (−1)^binom(p+1,2). N appears nowhere in that exponent, so a law that secretly depended on N reddens here |
| W3b | the seats-per-fold-coordinate profile differing from the pinned one. W3's own total already reddens on a narrowed sweep; what this row adds is WHICH p lost its seats, the largest being the one a reader can least check by eye |
| W3c | its own population pin, or a broken resultant. **Not a mutation of any R-symmetric family**: the vanishing at the centre seat is forced by the reflection, so this row is a READ that decides a word, not evidence |
| W4 | K1b's committed `K1B_SIGNS` not being reproduced from here, or the three readings not splitting as pinned |
| W4b | the flipped seats differing from the two pinned LISTS, or the two causes not moving six of fifty seats each over N = 4..12. It names seats and not counts where it can, because over N = 4..10 the two totals agree while four seats do not |
| W5 | one sector right half whose characteristic polynomial is not the middle route's comb of that parity, or a population other than the pinned 144 of which 110 have a nonempty right half rather than the trivial 1 = 1. With W5b it carries the parity jointly: swapping it inside the sector block reddens W5 alone, swapping it in the closed forms reddens both |
| W5b | one (N, seat, sector, k) where β_E fails to vanish at x_k for odd k, or β_O for even k, or a population other than the pinned 500. **This is the parity W6 cannot see**, W6's identity holding for any monic SQUAREFREE β; swapping the parity in the closed forms reddens W5 and W5b and leaves W6 green, while swapping it inside the sector block reddens W5 alone, so the two carry it jointly |
| W5c | one non-pole k where F157's seat-indexed Δ_k differs from the fold-indexed one, or a pole the two disagree about, or a population other than the pinned 113. It is read only at the 36 seats where the two are different expressions |
| W6 | one reading where Q_S is not c_S times the product over the NON-POLE roots of β_S, or the pinned 20 readings where a pole root leaves deg_t Q_S below deg β_S. It reads the product and then degree and leading coefficient; the identity holds for any monic SQUAREFREE β (W6c), so this row is about roots and never about which node indices they are |
| W0d | sympy turning out symmetric at equal degrees, or a population other than the pinned 12 asymmetric pairs of 200 with 55 at equal degrees. **A control**: W0c compares one fixed order and cannot see this, so the sentence that says it has its own reading |
| W6c | **a control.** W6's identity surviving a non-squarefree β, or the two draw populations differing from the pinned 54 squarefree with 0 failures and 26 non-squarefree with 4. It is what makes "any monic SQUAREFREE β" a measured hypothesis rather than a convenient word |
| W6b | the repeated-factor population differing from the pinned set, which a sweep narrowed below N = 10 empties |
| W7 | one seat where lc_t(P_j) ≠ (−1)^binom(p+1,2)·c_E·c_O, or where deg_t P_j ≠ n_E + n_O. It reddens only where W3 or W6's degree line does: given those two it is a theorem, and it is here so that the composition is written down and not inferred |
| W7b | the polynomial built here differing from F157's four committed rows by anything but the pinned sign |
| W7c | one seat where what this section factorises is not F157's own generator Res(S_n, Δ·S_j − S_{j+1}) times Corollary 11b's sign. Dropping the mirror-side term reddens it at half the mirror seats |
| W8 | **a control: it asserts that something BREAKS.** The R-breaking one-end family keeping the law at more than the three pinned seats. It is fed through the same door the real family uses, and feeding it the R-SYMMETRIC family reddens it |
| W8b | **a control.** Q_S taken right-half-first keeping the law at other than the pinned 12 of 50 seats, which would mean "fold half first" need not be said of the sectors |

Three mutations of the objects were run by hand against an earlier build of this gate and all
three went red: the wrong odd-sector sign, the crack fed through C1's door, and the criterion
without the ring-ends clause. None of the three is a hand mutation any longer. The first is gate
L1b, the second is what C2 asserts through C1's own door, and the third is gate T5.

Block K was mutated the same way, twelve times, by hand and with no committed artifact, against
the objects and against the field arithmetic: the fold coordinate, the split's coefficient, Lemma
5′'s sign, the gcd degree, the modular inverse's polynomial branch, the sector halves, the odd
sector's pair basis, the centre seat's forced count, the fence read at a generic knob, the ring
end read at u = 1/3 instead of u = 1, the wrap bond moved to (0, N−2), and the anisotropy moved to
the pair (0, N−2). Every check in the block reddens under at least one of them except K0c, K1d and
K1e2. Those three, and K1b2 and K3c beside them, are the block's controls: each asserts that
something BREAKS, so each is GREEN in a healthy block like every other check here, and red under
the input that would repair what it fences; a mutation of the object is not what tests them.
(This clause read "red exactly when the block is healthy" until 2026-09-03, which the committed
run contradicts, all five being PASS; §(h)'s S2d row would have inherited the same wording.)

Three of the twelve did not stay hand mutations, and each names a check that a whole class of
defect walks past. K3 is an equality between two readings of one routine, so a uniform miscount
leaves it green; K3c is the standing control that fires on what it cannot see. K1c and K1e are
statements about a tridiagonal matrix rather than about these families, so moving the
perturbation leaves them green; K1d and K1e2 break tridiagonality through the same door. An
earlier build of this block had a fourth, an anti-vacuity partner to K0 that asserted the split
with coefficient 1 must fail: given K0 that is equivalent to the struck characteristic polynomial
being identically zero, so no input could ever redden it. K0b now reads whether either summand
suffices alone, which the centre seat answers one way and every other seat the other.

## (k) Scope and fences

- **Multiplicity is in scope from §(g) on.** §(a) to §(f) compare SETS, with one exception that
  is not a comparison: gate C3 reads a COUNT on the ring, as a lower bound, to check the Remark
  of §(c). Nothing before §(g) compares two counts; §(g) compares COUNTS and is where F157's
  multiplicity is met. They are kept apart on purpose, because a locus is a set and a count is
  not read off a set: the odd-N equality of counts needs Lemma 5′, a bijection, where the odd-N
  equality of loci needed only Lemma 5. Where the two answers coincide is itself a result rather
  than an assumption: at even N one direction is proved in §(g) (loci parting forces counts
  parting) and the converse is measured over N = 6..14, where no seat carries the same locus and
  a different count.
- **What §(g) does NOT claim.** It compares the count of the crack at u with the count of the
  end-pair anisotropy at Δ = u, and nothing else. It says nothing about the count under any
  other perturbation, nothing about the Heisenberg book, and nothing about how the count varies
  ALONG either axis; F157's Δ_k comb remains the only closed form for that, and this file adds
  no closed form for the crack's count between the ring ends.
- **What §(i) does NOT claim.** It is about the anisotropy axis and the chain. It says nothing
  about the crack, where §(g) already records that the K1b identity fails at 20 of the same 32
  seats; nothing about F157's hop-2 normalisation, so the un-normalised 128 that
  [`BlindSeat.cs`](../../compute/MirrorWorld/BlindSeat.cs) carries at N = 9 seat 1 is not
  identified here; and nothing about which sign a reader would get from a third argument order.
  A sign law is a statement about a NAMED order, and §(i) names one.
- **F157's standing fence is not lifted here.** *"A Δ is NOT the detuned bond that The Seat That
  Cuts leaves open; do not report one as the other."* Lemma 1 identifies no operators: P is a
  bond and D a diagonal, they sit in different positions on the full space, and what coincides
  is their RESTRICTIONS to each reflection sector, at opposite signs of the knob in the odd one.
  Everything here reports loci, never one perturbation as the other.
- **The XY book only.** Every H here is the single-excitation XY block, F160's book. F157's
  Heisenberg law is on modulus N, a different comb, and is not on this road.
- **The theorem is about the two named families**, the wrap bond and the end-pair diagonal. The
  companion page measured four interior perturbations and found none that carries the whole
  Δ-locus, the best reaching 44 of 48 seats by containment. Which interior perturbations do, and
  at which seats, is untouched here. Lemma 1 says why the end pair is the one the crack meets,
  and says nothing about why a symmetric interior diagonal pair comes close.
- **"The end pair is special" is still refuted** for the locus, as the companion page found, and
  nothing here reinstates it. Lemma 1 is about the crack and the anisotropy sharing an end pair,
  not about the end pair being the only site pair that can carry a locus.
- **What is this file's own**: the comparison line in §(b), Lemma 2b, Lemma 4, Lemma 5's sector
  clause, the Remark in §(c), and the whole of §(g), that is Lemma 6, Lemma 7, the strengthening
  of Lemma 5 to counts, both corollaries and the reading of where the squaring happens; within
  §(g) the centre seat's Jacobi-end step is (J1) of the node-lemma proof and is cited there. In
  §(h) this file's own are the TRANSPORT of the sign law to the node modulus, that is the single
  sine on an end-detuned chain, together with the Remark and Corollaries 8a and 8b. The sign law
  (−1)^{k+1} itself is F71's and is owned in six places named in §(h); Δ_k, the pole rule and the
  multiplicity clause are F157's. In §(i) this file's own are Lemma 9, Lemma 10, Lemma 11 and
  Corollaries 10a, 10b, 11a and 11b, and the reading of what block K's committed sign literal is
  a reading OF; Lemma 11 stands on the two node readings gate S2 pins for §(h), and the parity
  the factorisation needs is §(h)'s Lemma 8 read on a uniform path, not a new sign. The
  Cassini identity Lemma 9 turns on is `PROOF_CRACKED_RING_EXACT_CURVE` §(b)'s and is cited, not
  repeated; the Chebyshev addition formula is standard and is written out in Lemma 9 because no
  proof in this repository carries it.
  Lemma 1's folded half is F160's, Lemma 2a is
  `THE_SEAT_THAT_CUTS` §7's, Lemma 3 is Lemma J of the node-lemma proof, and the Σ identities are
  the companion page's.
- **Ordinary words here that are typed objects elsewhere, and none of them is meant.** A **pair** on
  this page is the site pair {p, N−1−p} the reflection joins, never
  [`Pair.cs`](../../compute/MirrorWorld/Pair.cs)'s bare coherence |i⟩⟨j| with its rate −2γk; the
  **end pair** is F157's {0, N−1}, which the companion page already fences against F140's
  "corner block". The **divisor law** §(h)'s third fence names is F157's gcd law on the seat
  index, never [`Divisor.cs`](../../compute/MirrorWorld/Divisor.cs)'s frozen divisor of F140,
  which lives on the R₉₀ locus and has no seat in it. A **block** here is the matrix of H
  restricted to one reflection sector, never
  [`Block.cs`](../../compute/MirrorWorld/Block.cs)'s joint-popcount block (p, q). A **sector**
  here is always a reflection sector: the whole file sits inside one excitation number, so no
  popcount grading is in play anywhere in it, and neither
  [`JointPopcountSectors.cs`](../../compute/RCPsiSquared.Core/BlockSpectrum/JointPopcountSectors.cs)
  nor its witnesses are meant; "block K", "block S" and "block W" are the one place a **block**
  here is a lettered group of gates rather than a sector's matrix. A **comb** is F157's node comb, the Δ_k indexed by
  k ∈ 1..N_node−1, and in §(i) also the middle route's own Dirichlet comb of x-values; never
  [`Cyclotomy.cs`](../../compute/MirrorWorld/Cyclotomy.cs)'s two combs, which are turn fractions
  and not Δ values. A **pole** is F157's pole INDEX, a k at which S_p(x_k) = 0 so that Δ_k is not
  defined and the degree drops; never
  [`Cat.cs`](../../compute/MirrorWorld/Cat.cs)'s two immortal poles. The **folded half** is
  F160's geometric fold along R and the **fold coordinate** is min(j, N−1−j); neither is
  [`Mirror.cs`](../../compute/MirrorWorld/Mirror.cs)'s f_P and f_Q, which pay λ → −λ − 2Nγ, and
  no spectrum is folded here. §(i) calls Lemma 10's object **the common factor** and never a
  divisor, because F139 uses "Chebyshev divisor" in the same polynomial-algebra sense one file
  over; The uniform stretch between a seat's two **mirror images** is
  named **the middle route** so that it has a name at all; a mirror image here is a site's partner
  under the chain reflection R and nothing more, never
  [`Mirror.cs`](../../compute/MirrorWorld/Mirror.cs)'s block-lattice group of eight, which pays
  λ → −λ − 2Nγ, and never Π; §(i) does also call it a path, in the
  ordinary graph sense the Statement uses for the open N-site path, and what is never meant is
  F89's `path-k` topology or
  [`Cyclotomy.cs`](../../compute/MirrorWorld/Cyclotomy.cs)'s `PathOrders`. No γ
  and no Liouvillian appears anywhere in this file.
- **Bare-letter lemma names are ambiguous in this neighbourhood** and this file numbers its own
  for that reason. Three of the cited proofs have a Corollary C, two of them a Corollary B, and
  two a Lemma B, all of different objects, `PROOF_COLLISION_GAP_ODD_ORDERS.md`'s Lemma B being one
  this page quotes; `PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md` adds a fourth object under the same
  letter. Every citation here names the file.
- **No hardware claim.** Nothing here is proposed for a flight.

## Where it came from

The companion page landed on 2026-09-02 with §(c) as an observation and its open item 1 carrying
two questions, the first asking for a proof of the shape
u-locus = (even-sector Δ-locus) ∪ −(odd-sector Δ-locus). That shape is correct
and is Lemma 1 plus Lemma 2b. What the item did not anticipate is that the same lemmas settle the
even-N side as well, converting its own gate C8 from a measured biconditional into a consequence,
and that the appendix {+1, −1} is not an appendix but the third term of the decomposition.

The item's second question was written the same day and in the opposite spirit, as a limitation
rather than a question: that page recorded of its own construction that it *"cannot see it (the
resultant is a multiple of P_j squared), so a different one is needed"*. No different construction
was needed. What was needed was to split the
seat the way Lemma 2b had already split the kernel, and the squaring then turns out to be the
answer rather than the obstacle: it is the sector's own strike disconnecting a Jacobi block, so
the count is in the resultant doubled. The item's own diagnosis was right about the polynomial
and wrong about what follows from it.

Related: [The Blind Seat on the Road](../../experiments/THE_BLIND_SEAT_ON_THE_ROAD.md) ·
[PROOF_CRACKED_RING_EXACT_CURVE](PROOF_CRACKED_RING_EXACT_CURVE.md) ·
[PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA](PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md) ·
[PROOF_COLLISION_GAP_ODD_ORDERS](PROOF_COLLISION_GAP_ODD_ORDERS.md) ·
[The Seat That Cuts](../../experiments/THE_SEAT_THAT_CUTS.md) ·
[The Cracked Bell](../../experiments/THE_CRACKED_BELL.md)
