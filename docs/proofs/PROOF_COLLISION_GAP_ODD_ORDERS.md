# PROOF: Each Computed Order of a Collision Gap Reads One Comb Under Integer Multipliers, and Which Orders Vanish Is Decided by a gcd

**Registry:** [F161](../ANALYTICAL_FORMULAS.md#f161).
**Status:** Tier 1 derived. Theorem A (the u-series of a chain-end level to fifth order)
is a consequence of F160's polynomial identity and holds for every N ≥ 3 and every mode;
it is symbolic in n, with the sign η treated as the number ±1 and both values checked
separately. Lemma B (the wrap bond is rank one inside each reflection sector) is exact
for every N and is the reason only the ODD orders carry η. Corollary C (the multiplier
form) is a change of trigonometric basis and carries no hypothesis. Theorem D (the
Galois kill) holds at every ODD n for every collision pair, at every rung j with
gcd(j, n) = 1, and its converse direction is an equivalence, which is what makes the
negative control sharp. Theorem E (the ROT3 rung lemma) is proved for a parity-uniform
triple whose doubled-label set is a union of two order-3 cosets, and only the direction
3 ∤ j is forced: at 3 | j the value is a cosine that may or may not vanish, and both
outcomes occur in this census. Theorem F (where the ladder stops) is unconditional, from
F129's own firing condition. Corollary G (the second order) starts from Theorem D on the
other ladder and then sharpens the criterion: what a multiplier must be invertible on is
every minimal PIECE of the pair's vanishing sum, not the whole ring, and under that reading
Theorems D and E are one criterion at its two ends. The sufficient direction is a proof; the converse at m = 3 is checked on all 2558 pairs of the
census, against piece signatures that are themselves decided at two primes rather than exactly,
and identifies them as one named F129 family. **What is READ and not derived:** that all twelve even-n non-mirror
standing pairs of the n ≤ 30 census have both triples of the ROT3 shape AND share a parity
class. Theorem E turns the shape into the vanishing of c₃ and the parity match into the
vanishing of c₁; nothing here says a standing pair must have either.
**Date:** 2026-09-02.
**Authors:** Thomas Wicht, Claude (Opus 5).
**Script:** [`simulations/collision_gap_odd_orders.py`](../../simulations/collision_gap_odd_orders.py)
(fifteen gate blocks L1 to L14, about fourteen seconds measured quiet; `--fast` SKIPS L4 and the
two larger moduli of L8, so it passes without the end-to-end anchor and without most of the
negative control, and buys almost nothing. Every gate is exact except L4, an error-model law on
an eigensolver, where no exact route to an eigenvalue exists;
run committed at
[`collision_gap_odd_orders_run.txt`](../../simulations/results/collision_gap_odd_orders/collision_gap_odd_orders_run.txt)).
It imports the exact cyclotomic ring and the census from
[`f129_level_collision_law.py`](../../simulations/f129_level_collision_law.py) and
[`comb_road_f129.py`](../../simulations/comb_road_f129.py) rather than restating either.

**Builds on:**
- [PROOF_CRACKED_RING_EXACT_CURVE](PROOF_CRACKED_RING_EXACT_CURVE.md) Theorem A for the
  polynomial det(x·I − H(u)) = U_N(x/2) − u²·U_{N−2}(x/2) − 2u, which is the whole input,
  and Theorem G for the first order, which is this file's Theorem A read at m = 1.
- [THE_COMB_ON_THE_ROAD](../../experiments/THE_COMB_ON_THE_ROAD.md), whose census this
  file explains rather than repeats: its gate R3b is Theorem D at j = 1, its R4c is
  §(f)'s term-by-term argument for the Θ-mirror pairs, and the 212 pairs it recorded as reported
and not explained are what this file answers. That page is edited by the same change that lands this one, so its
  sentence is no longer there to quote; what it said is what is answered here.
- [PROOF_F129_LEVEL_COLLISION_LAW](PROOF_F129_LEVEL_COLLISION_LAW.md) for the firing
  condition 3 | n or 10 | n, which Theorem F uses and does not reprove.
- [F89_SEED_EXISTENCE_REDUCTION](../../experiments/F89_SEED_EXISTENCE_REDUCTION.md) for
  the Conway-Jones classification of vanishing cosine triples into TRIV, ROT3 and PENT (the
  first two unions of cosets of a prime-order subgroup, PENT not), and for the general
  Galois criterion *"every Galois automorphism σ_j: ζ ↦ ζ^j
  (j coprime to 2(N+1))"*. Theorem E's ROT3 is that file's family, cited by its name.
- [PROOF_K_PARTNERSHIP](PROOF_K_PARTNERSHIP.md) for R·ψ_k = (−1)^{k+1}·ψ_k, the chain
  reflection's action on the sine modes, which is Lemma B's only input beyond arithmetic.

## Notation, and one trap before anything else

**n is the comb modulus, n = N + 1**, as in the whole F129 stack
(`f129_level_collision_law.py`, `LevelCollisionCensus.cs`, `LevelCollision.Ncomb`). In the
crack CODE the same letter is the SITE COUNT: `Crack.cs` opens its own workings with
`int n = N;`. The crack PROOF never uses lowercase n for the site count, writing that as capital N throughout
and spending lowercase n on a free Chebyshev degree and once on a site occupation n(j), so the
trap is between this file and `Crack.cs`, not between this file and that one. Every formula below is in the comb
book; a formula carried into the code without changing the letter is off by one.

N sites, nearest-neighbour hopping 1, single-excitation N × N block. The wrap bond between sites
0 and N − 1 carries strength u, so H(u) = H_chain + u·V with V = |0⟩⟨N−1| + |N−1⟩⟨0|. At u = 0
the eigenvectors are the sine modes ψ_k(i) = √(2/n)·sin((i+1)·kπ/n), i = 0..N−1 (F2b), with
levels E_k = 2·cos θ_k, θ_k = kπ/n, k = 1..N, labelled by DESCENDING energy. η_k = (−1)^{k+1}
and a_k = ψ_k(0) = √(2/n)·sin θ_k, the endpoint amplitude.

For a triple τ ⊂ {1..n−1} and an integer j,

    M_j(τ)   = Σ_{k∈τ} cos(j·kπ/n),          X_2j(τ) = Σ_{k∈τ} η_k·cos(2j·kπ/n).

Δ always means τ minus σ. A COLLISION pair is two distinct clean triples with ΔM_1 = 0
(F129's level equality; *clean* is F129's condition k_i + k_j ≠ n). The GAP is
g(u) = Σ_{k∈τ} E_k(u) − Σ_{k∈σ} E_k(u), and c_m is its u^m coefficient. That letter is spent
in F160 on a different quantity, its split correction c_m = ½ − 1/(N·sin²k_m), where m indexes
a MODE and not an order; the two never meet on this page, but the letter is the same and the
reader should know it. *Rung* means one
value of j, and NOT the glossary's rung, which is a rate (2γ, 2Nγ) or a Q threshold, nor
the arc `sideways_spin_ladder`'s unit of checking. *Ladder* means the sequence of evaluations
M_{n+2j} of this one comb, and nothing else: not F142's spin ladder, not
`Pi2DyadicLadderClaim`'s dyadic ladder, not `ClockHandLadderClaim`'s. *Road* is F160's u axis. *Θ-mirror pair* means σ = n − τ and is
[PROOF_ZETA2_ANTI_PROTECTION](PROOF_ZETA2_ANTI_PROTECTION.md)'s term, used in its sense. A
*pair* here is two triples, THE_COMB_ON_THE_ROAD's fence, not MirrorWorld's `Pair` (a bare
coherence) and not F75's mirror pair (two sites).

Several letters do double duty and are worth naming once rather than renaming. A and B are the
two Chebyshev polynomials here and the two reflection-sector factors of
PROOF_CRACKED_RING_EXACT_CURVE's G = 2AB. P is the characteristic polynomial in §(a) and the
±label set in §(f). m is an expansion order in c_m, d_m and a multiplier in *the map k ↦ k·m mod
2n*. s is sin θ_k from §(c) onward and nothing else; the perturbation scale in §(b) is written
out as 2ηu and the angle equation's is v = ηu. j is a rung index everywhere except in Lemma B's
*the chain reflection j ↦ N − 1 − j*, where it is a site. a is the endpoint amplitude in §(b)
and a coset label in §(f). D is three objects and none of them is this file's: D_m in §(b) is
the η-free part of d_m, PROOF_CRACKED_RING_EXACT_CURVE's D_n is its path determinant and its
Theorem D is the join, and THE_COMB_ON_THE_ROAD's D is what this file calls c₁. All seven theorem letters here overlap
PROOF_CRACKED_RING_EXACT_CURVE's, so *Theorem A* is qualified by its file wherever both are in
one sentence. Two more, unavoidable and local: R is the chain reflection in §(b) and R_k the
sector resolvent three lines later; σ is the second triple throughout and σ_m a Galois
automorphism in §(e).

## What was learned

Switching on the wrap bond makes each order of a level's motion re-read the chain's own mode
comb at a higher harmonic, and a level collision, being a rational identity among roots of
unity, is invisible to every harmonic that merely permutes that comb. Five orders are computed
here; that the pattern is the general shape of the series is not proved.

Which harmonics merely permute is a gcd, so the orders that vanish are decided by arithmetic and
not by the physics. For a pair that STANDS at first order the surviving harmonics run out until
the first one that fails to permute, which is 3; and at an ODD modulus 3 divides n only because F129 needed it there for the collision to exist
at all, so the arithmetic that builds the coincidence is the same arithmetic that ends its
protection. For the other 2335 pairs the gap opens at first order
and there is nothing to protect. The arc `the_forced_and_the_met` had two of these clauses before this file, and only the first
order: *the coincidence met, its derivative forced*. This file adds the third, that the ORDER at
which the forcing runs out is forced too, and Corollary G then takes the last coincidence the arc still called *met* and shows it forced,
after which the word is retired from the sentence as a vocabulary decision, recorded in this
date's error log: on this object *met* has never once survived being looked at, so it is written
as *no forcing yet found*, a claim about us. The three-clause version now in the arc is this
change's own wording and not prior work; an earlier draft of this paragraph quoted it back as
though it were, having taken a review lens's citation of the staged file for a citation of the
repository.

## What the repo already held, store by store

Swept 2026-09-02 by three agents, one per layer, before any of the algebra below was
written down, and the finding that mattered most was a negative: the pieces are all
committed and the step between them is not.

- **[`docs/ANALYTICAL_FORMULAS.md`](../ANALYTICAL_FORMULAS.md).** F160 holds the road and
  its curve; its entry carries the second-order coefficient of the SPLIT at the ring end,
  c_m = ½ − 1/(N·sin²k_m), which is a different expansion point and a different object
  from anything here (fenced again in Scope). F65 holds the endpoint rate comb
  α_k/γ₀ = (4/n)·sin²θ_k and, in its Niven face, writes it on the DOUBLED angle
  (2/n)(1 − cos 2θ_k); that combination is exactly d₁/η below, the X₀ and X₂ rungs together,
  and F65 is where the doubled angle entered the arc. The sign η = (−1)^{k+1} is F71's, carried in F75's bonding-mode specialisation and in
  PROOF_K_PARTNERSHIP's closing section. F77 owns the identity
  sin⁴ = (3 − 4cos 2x + cos 4x)/8, in its entry in those words, but spends it on a sum over
  SITES at a fixed mode (an entropy sum), not over modes; the identity is cited here, the
  object is not the same. F129 holds the law under test and its firing condition. F123 is the Re-side sibling of the velocity in the MOVE and not in the form,
THE_COMB_ON_THE_ROAD's own qualifier: a Hellmann-Feynman derivative in a bond's knob, but of a
decay rate and on a different object. **Nothing** in the registry holds a
  second- or third-order coefficient of a single level at the CHAIN end, and nothing
  anywhere pairs cos(4kπ/n) with a level shift.
- **[`docs/proofs/`](.).** PROOF_CRACKED_RING_EXACT_CURVE's Theorem E is a genuine
  second-order elimination, but in δ = 1 − u about the RING end and for a pair split.
  Theorem G is the first order at the chain end, by two routes. Corollary B's simplicity
  clause already owns the reflection split of the crack, as the real-space factorization
  G = 2AB into an R-even and an R-odd Jacobi block; Lemma B below is the same fold seen in
  the mode basis, and says in one line what that factorization says in the site basis.
  PROOF_K_PARTNERSHIP carries R·ψ_k = (−1)^{k+1}·ψ_k in its closing section and credits it
  to F71; the derivation home is PROOF_C1_MIRROR_SYMMETRY §2.
  PROOF_MIRROR_ORDER_SORTING (F131) is the repo's general home for *which orders may enter a
response*: its §1 sorts them by a parity into four cells, of which one is *generic (no
constraint)*, and its Theorem A is one column of that table. This file's first version did not
cite the file at all; a review lens supplied it. The eleven Θ-mirror pairs are an instance of
its EVEN cell, with M = K and the scan parameter u. The other 212 get nothing from it: they fall
under the cell that imposes no constraint, and they can fall nowhere else, because F131's
carrier is an involution and grades orders by ℤ/2 at most, while *orders 1 and 3 vanish and 2
and 5 do not* is no ℤ/2 grading. The carrier here is a Galois group, (ℤ/2n)*, which can leave an
arbitrary set of orders standing. That is the sharpest available statement of what this file
adds, and the limit is F131's by construction, not a defect of it.
  PROOF_ZETA2_ANTI_PROTECTION owns an all-orders IDENTITY on the ζ axis,
  θ_τ(ζ) = −θ_ν(−ζ), whose consequence splits by order and is its §4 heading, the title naming the even half only: the ODD
  orders cancel in the mirror pair's difference (the protection) and the EVEN ones double
  in it (the anti-protection). It is the shape §(f)'s Θ-mirror argument reproduces on the u axis, with u for ζ; that argument
is term-by-term and is NOT a case of Theorem E, one of the eleven pairs failing Theorem E's
hypothesis. **Nothing** in `docs/proofs/` holds a third-order coefficient
  of a level in a bond knob, and **nothing** holds any of the five comb laws under a bond
  detuning.
- **[`experiments/`](../../experiments/), the prior flights included; no null result in that
store bears on this object, which is a statement about a store that was opened and not a claim
that none exist.** THE_COMB_ON_THE_ROAD is the census this file explains; before the change that lands this file,
its record read the 212 as reported and not explained, and its R4 measured, by 40-digit decade
ratios,
  that exactly 60 pairs have a vanishing second-order coefficient. Corollary G reproduces
  that 60 by an exact route that never sees a float, which is the strongest independent
  check in this file and was not looked for. The per-modulus collision counts that census
  column rests on are themselves a committed closed form,
  [`f129_family_inventory.py`](../../simulations/f129_family_inventory.py) and F129's family
  inventory, so the 1313 at n = 30 is not this file's measurement either. F89_SEED_EXISTENCE_REDUCTION classifies the vanishing triples into TRIV, ROT3 and PENT, the first two unions of cosets of a
  prime-order subgroup and PENT not, and states the general Galois criterion this file's Theorem D uses.
  COUPLING_DEFECT_WALK_TIME_STEP carries a Taylor series in the same bond knob, of the
  Wigner delay and not of a level. DICKE_VS_ENDPOINT_PROBE_JW_N11 item 5 holds the F71 mirror sign rule for a bond matrix element
in the sine-mode basis, and THE_COMB_ON_THE_ROAD carries the DIAGONAL element of this very V,
⟨ψ_k|V|ψ_k⟩ = 2ψ_k(0)ψ_k(N−1), with the reflection substitution and the sign already named. What
**nothing** holds is the OFF-diagonal element and the rank-one-per-sector consequence, which is
Lemma B.
- **[`hypotheses/`](../../hypotheses/).** PERSPECTIVAL_TIME_FIELD §3.2 (Tier 2) carries
  the first-order bond overlap 2ψ_k(b)ψ_k(b+1), instantiated at the TERMINAL bond (0,1),
  which shares site 0 with the wrap bond, and at the central bond (3,4) as its control.
  Nothing beyond first order OF A LEVEL there; the file does reach second order, on a
  rescaling sum, where it reports a partial second-order cancellation.
- **[`fw.Confirmations`](../../simulations/framework/confirmations.py) and the C#
  `ConfirmationsRegistry`.** Twenty-four entries in each. Confirmation 24,
  `f129_standing_fringe_kingston_july2026`, IS this census's n = 9 pair (1,5,7) ~ (2,4,8),
  flown 2026-07-15 with ZZ and Z-detuning as its perturbations. It is a sighting of the
  collision on the comb, on a different axis; nothing here proposes a flight and nothing
  here is claimed to be confirmed by it. Of the other twenty-three, nothing ring-shaped
  and nothing bond-defect-shaped.
- **The OpenArcs registry
  ([`OpenArcsRegistry.cs`](../../compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs)),
  `the_forced_and_the_met`.** Its NextStep (1) records the comb test as RUN and quotes its
  own result including *c3 = 0 for all 223*; this file is that item's explanation and not
  a second run of it. Its fence *"THE_SEAT_THAT_CUTS's own open item asks for a DETUNED
  BOND rather than a Delta, so do not report one as the other"* is respected: the bond here
  is the ring's wrap bond under a level law. `sideways_spin_ladder`, `f89_galois_open_doors`
  and `compressed_density_laws` were opened by name. All three turn a coupling and none
  expands a level in one: the first runs its 220 checks at NON-UNIFORM J_b and γ_l but reads
  a defective locus in q = J/γ, the second is Galois on the relaxation spectrum in q with its
  threads B and D unstarted, the third carries a bond control on the R₉₀ locus and a PAIR-sum
  injectivity argument.
- **[`docs/GLOSSARY.md`](../GLOSSARY.md) and the words.** No *multiplier*, and that one is
  free. *Rung* is NOT free and the first version of this bullet said it was: the glossary
  spends it four times on the rate ladder (*"the ladder's top rungs (2Nγ, and F3's 2(N−1)γ
  …)"*, *"between the coherence-horizon rungs Q\*(2) = 1 and Q\*(3) = √2"*, *"keeps its
  minimum at or above its rung at every coupling"*, *"the 2γ rung exists at every
  coupling"*), where a rung is a RATE or a Q threshold, and the arc `sideways_spin_ladder`
  counts its 220 checks in rungs as well. The fence is in Notation. Same for *ladder*,
  *road*, *pair* and *Θ-mirror pair*; THE_COMB_ON_THE_ROAD's own fence list is inherited
  unchanged.
- **[`docs/CAUGHT_ERRORS.md`](../CAUGHT_ERRORS.md).** The FIRST entry of 2026-09-02 is
  this object's own previous round; the second is this file's, written after it. The trap CLAUDE.md names, a law concluded from an input that could not break it, was walked into
again while deriving this file and is recorded where it happened, in §(f) and §(k) here: the
ROT3 converse was read off n = 30, a modulus that cannot break it, and n = 24 does. This date's
second entry, which this change adds, carries the twenty-six shapes the review rounds found.

## (a) The object

Along the road the single-excitation block is H(u) = H_chain + u·V. F160's Theorem A gives
its characteristic polynomial as an exact polynomial in x and u,

    P(x, u) = det(x·I − H(u)) = A(x) − u²·B(x) − 2u,     A = U_N(x/2),  B = U_{N−2}(x/2),

with U the Chebyshev polynomial of the second kind in the 2cos normalization. Two features
of P are used below and are worth naming before they are: the u-dependence enters only
through u² and through the single term −2u, and that term is a constant in x.

At u = 0 the levels are the comb E_k, simple, and by Corollary B of the same file the spectrum
stays simple at every u ≥ 0 with u ≠ 1, so each level continues analytically and the labels k =
1..N by descending energy carry along the road, the single point u = 1 excepted. The expansion below is at u = 0.

## (b) Lemma B: the wrap bond is rank one inside each reflection sector

Let R be the chain reflection j ↦ N − 1 − j. It maps every chain bond to a chain bond and
the wrap bond to itself, so it commutes with H(u) at every u, and on the sine modes
R·ψ_k = η_k·ψ_k (PROOF_K_PARTNERSHIP, cited). In particular ψ_k(N−1) = η_k·ψ_k(0) = η_k·a_k.
Then

    ⟨ψ_k|V|ψ_l⟩ = ψ_k(0)·ψ_l(N−1) + ψ_k(N−1)·ψ_l(0) = a_k·a_l·(η_k + η_l).

**Two consequences.** V vanishes whenever η_k = −η_l, so it never connects modes of
opposite reflection parity. And inside the sector where the common value is η, it is
2η·a·aᵀ: rank ONE, the outer square of the endpoint amplitudes. So on each sector the perturbation is u·V = 2ηu·a·aᵀ, and every Rayleigh-Schrödinger coefficient
in that scale is built from the sector's amplitudes and gaps alone, η entering only as the
factor 2ηu.
Which sector one is in is decided by the parity of k, so this alone gives d_m = η^m·(a function
OF THE SECTOR); the first non-trivial one is the SECTOR RESOLVENT R_k = Σ_{l ≠ k, l ≡ k (mod 2)} a_l²/(E_k −
E_l), the sum over the mode's OWN sector; the other sector does not enter, V having no matrix
element into it. R_k is not an independent claim here: the second-order Rayleigh-Schrödinger
coefficient of a rank-one perturbation of scale 2ηu is d₂ = (2η)²·a_k²·R_k = 4·a_k²·R_k, and
with d₂ from L2 and a_k² = (2/n)·sin²θ_k this reads back as R_k = (n−3)·cos θ_k/(2n), the same
closed form in either sector.

That the function is the same for both sectors AT EVERY ORDER, and not only in the five computed
below, has a one-line proof that does not go through the sectors at all. Put x = 2cos φ with φ =
θ_k + ε in the characteristic equation multiplied by sin φ, sin(nφ) − u²·sin((n−2)φ) − 2u·sin φ
= 0. Using sin(nθ_k) = 0 and cos(nθ_k) = −η and dividing by −η, it becomes

sin(nε) − v²·sin(nε − 2φ) + 2v·sin φ = 0, v = η·u,

in which η appears NOWHERE ELSE. The left side vanishes at (ε, v) = (0, 0) and its ε-derivative
there is n ≠ 0, so the implicit function theorem gives a unique analytic branch ε(v; θ_k, n)
with ε(0) = 0. It depends on v, θ_k and n, and on η only through v; hence so does δ = 2cos(θ_k +
ε) − 2cos θ_k, and d_m = η^m·D_m(θ_k, n) for every m with D_m free of η. Gated at L2b. Therefore

    d_m = η^m · (an η-free function of θ_k and n),

which is why the odd orders carry the reflection sign and the even orders do not. Gate L2b
checks the matrix element exactly for N = 4..9 and checks the η-power structure on the
five coefficients derived below.

This is the mode-basis face of the fold that PROOF_CRACKED_RING_EXACT_CURVE's Corollary B
uses in the site basis, where the same reflection splits the curve as G = 2AB. The rank is
two overall, as that file states, because there are two sectors; it is one in each.

## (c) Theorem A: the u-series of a chain-end level, to fifth order

Write x = E_k + δ and solve P(E_k + δ, u) = 0 order by order in u. Every derivative of A
and B needed is a closed form at a comb root, obtained by working in θ with
x = 2cos θ, so that d/dx = −(1/(2 sin θ))·d/dθ, on

    A(θ) = sin(nθ)/sin θ,     B(θ) = sin((n−2)θ)/sin θ,

and then substituting the two facts that define the root: sin(nθ_k) = 0 and
cos(nθ_k) = (−1)^k = −η_k. The first three are

    A(E_k) = 0,     A′(E_k) = n·η/(2 sin²θ),     B(E_k) = 2η·cos θ,

and the rest are produced by the script rather than typed here; A″, A‴, A⁗, A⁗′ and B′, B″,
B‴ enter only through the elimination. Solving to fifth order gives, with s = sin θ_k and c = cos θ_k,

    d₁ = η·(4/n)·s²
    d₂ = (4(n−3)/n²)·s²·c
    d₃ = (4(n−2)(n−4)/(3n³))·η·s²·(4c² − 1)
    d₅ = η·F·[ (5/4)(n−1)·cos 2θ + (3n²/2 − 8n + 8)·cos 4θ − (3/4)(2n−3)(n−3)·cos 6θ ],
         F = 4(n−6)(n−2)/(15 n⁵)

with d₄ computed and not displayed: it carries no η and plays no part in an odd order, but
it is needed to reach d₅ and the script keeps it. d₁ is F160's Theorem G, recovered here as
the first term of the same series, and its magnitude is F65's α_k/γ₀ with F71's mirror sign, the one F75's bonding-mode specialisation
carries.

**Gates.** L1 pins P against the actual matrix for N = 3..9, exactly, so the polynomial is
this file's own and not a quotation. L2 pins d₁, d₂, d₃ and d₅ symbolically in n, at η = +1
and η = −1 separately, because η is a sign and not a symbol that sympy may square at will.
L4 is the end-to-end check that the series is the real matrix's: at N = 14 and N = 19 the
residual of the five-term series against eigenvalues of H(u) computed at 60 digits (the
coefficients themselves at 50) is divided by u⁶ at FOUR values of u, 1e-3 to 1e-6, giving q(u) → d₆. The error model is q(u) =
d₆ + c·u + O(u²), so successive differences fall by one decade and the RATIO of differences
approaches 1/10 from its own O(u). Four values give two such ratios, and the gate is that the second is nearer 1/10 than the first and that both are within 0.002 of it,
which is the model converging rather than one quotient landing in a window. The three cases are
(N, k) = (14, 4), (14, 9) and (19, 7), and each moves a full decade closer: 0.099983508 to
0.099998346, 0.09995787 to 0.09999578, 0.099953513 to 0.099995348. That is the model itself,
not a window a single quotient has to sit inside. There is no exact route to an eigenvalue, so this is one of the two gates on the page that are
not exact; the other is L13, whose piece signatures come from a decomposition that decides
vanishing at two primes (see Scope).

## (d) Corollary C: each computed coefficient is a signed combination on a multiplier ladder

Rewrite each d in the cosine basis: s² = (1 − cos 2θ)/2, s²c = (cos θ − cos 3θ)/4, and
s²(4c² − 1) = (cos 2θ − cos 4θ)/2. Then

    d₁ = (2/n)·η·(1 − cos 2θ)
    d₂ = ((n−3)/n²)·(cos θ − cos 3θ)
    d₃ = (2(n−2)(n−4)/(3n³))·η·(cos 2θ − cos 4θ)

and d₅ is already in that basis above. Summing over a triple and differencing gives, for a
collision pair,

    c₁ = (2/n)·(ΔX₀ − ΔX₂)
    c₂ = ((n−3)/n²)·(ΔM₁ − ΔM₃) = −((n−3)/n²)·ΔM₃          since ΔM₁ = 0 IS the collision
    c₃ = (2(n−2)(n−4)/(3n³))·(ΔX₂ − ΔX₄)
    c₅ = F·[ (5/4)(n−1)·ΔX₂ + (3n²/2 − 8n + 8)·ΔX₄ − (3/4)(2n−3)(n−3)·ΔX₆ ]

Through third order each of these is a DIFFERENCE OF TWO neighbouring rungs. c₅ is not: it
carries three rungs with unequal coefficients, and nothing here says the pattern of two
continues past c₃. What does hold in every case computed is that the rungs are consecutive
and that the coefficients sum to zero.

Two ladders, and whether they are the same one is decided by n. The EVEN orders sit on M with
ODD multipliers 1, 3, 5, …, and the ODD orders on X with even ones. That split is a theorem at every order, not a pattern in the five computed, and the proof is in
the angle equation of §(b) rather than in the comb. Write E(θ, ε, v) = sin(nε) − v²·sin(nε − 2φ)
+ 2v·sin φ with φ = θ + ε. Two substitutions send E to −E and therefore fix its zero set:

(θ, ε, v) ↦ (π − θ, −ε, −v) and (θ, ε) ↦ (−θ, −ε).

By the uniqueness of the analytic branch ε(v; θ) with ε(0) = 0, the first gives ε(−v; π − θ) =
−ε(v; θ), hence δ(−v; π − θ) = −δ(v; θ), hence **D_m(π − θ) = (−1)^{m+1}·D_m(θ) as functions of
a free θ**; the second gives D_m(−θ) = D_m(θ), so D_m is even in θ and a cosine polynomial in
the first place, which the conclusion needs and which nothing else on this page supplies. A term
cos(rθ) in D_m therefore requires r ≡ m + 1 (mod 2). Both substitutions are gated at L2b,
symbolically in n, and neither uses the comb points, the chiral K, the parity of N, or any bound
on the degree of D_m. The chiral K reaches the same conclusion through the labels, but only AT
the n − 1 comb points, from which the statement about the FUNCTION would need exactly the degree
bound §(k) leaves open; the first version of this paragraph took that route and a review lens
found the gap. X is then the comb
read under a SHIFTED multiplier. The shift is η itself: η_k = −cos(n·kπ/n) is the comb read at multiplier
n, so attaching the reflection sign to a harmonic just translates its index by n. At ODD n
that keeps it odd, both ladders are rungs of the one odd-multiplier ladder, and both obey
gcd(m, n) = 1; at even n the translation throws the odd orders onto the EVEN multipliers,
where no automorphism lives at all. That is why the parity of the modulus decides which half
of the series can be killed, in one line rather than a table.

    η_k·cos(2j·kπ/n) = −cos(k(n+2j)π/n),    so   X_2j(τ) = −M_{n+2j}(τ).

That identity is the comb page's R3b move at j = 1, written for every j. Note also that d₅ carries no X₀ part, visible in the form displayed above and pinned by L2,
and that its three coefficients sum to zero together with that absent fourth one, gated at
L3; the second is just d₅ → 0 as θ → 0, where the mode has no endpoint amplitude.

## (e) Theorem D: the Galois kill, and its converse

The map k ↦ k·m mod 2n induces the Galois automorphism σ_m: ζ_2n ↦ ζ_2n^m of ℚ(ζ_2n)
exactly when gcd(m, 2n) = 1 (the forward direction is F89_SEED_EXISTENCE_REDUCTION Step 4, *every Galois automorphism σ_j: ζ ↦
ζ^j (j coprime to 2(N+1))*, cited; the equivalence is the standard description of Gal(ℚ(ζ_m)/ℚ) as
(ℤ/m)* and is not that file's). For the X ladder
m = n + 2j, and

    gcd(n + 2j, 2n) = gcd(j, n)     for every odd n,          gated at L7 for n = 9..59

because at odd n the number n + 2j is odd, so the factor 2 never contributes, and
gcd(n + 2j, n) = gcd(2j, n) = gcd(j, n). For the M ladder m = 2j + 1 is odd already and the
criterion is gcd(2j + 1, n) = 1.

**Statement.** Let (τ, σ) be a collision pair at odd n, so ΔM₁ = 0. For every rung j with
gcd(j, n) = 1,

    ΔX_2j = −ΔM_{n+2j} = −σ_{n+2j}(ΔM₁) = 0.

At odd n both j = 1 and j = 2 always satisfy gcd(j, n) = 1. Hence **ΔX₂ = ΔX₄ = 0 for every
collision pair at odd n**, and therefore

    c₁ = (2/n)·ΔX₀ = (4/n)·(o_τ − o_σ)      o = the number of odd labels, since Σ_τ η = 2o − 3
    c₃ = 0

The first line is THE_COMB_ON_THE_ROAD's gate R3b, recovered as the j = 1 rung. The second
is new, and it is not a property of the pairs that stand: **c₃ = 0 holds for all 627
collision pairs at odd n, standing or separating.** Gated at L6.

**The converse, which is what makes the control sharp.** σ_{n+2} is an automorphism of ℚ(ζ_2n),
hence injective, so ΔX₂ = 0 forces ΔM₁ = 0. At odd n, therefore, ΔX₂ = 0 **if and only if**
the pair collides. Gate L8 is the negative control the arc asks for, a configuration where
the law must not fire: over all 496 clean-triple pairs at n = 9, all 39060 at n = 15, and
40000 sampled at each of n = 21 and n = 27, ΔX₂ vanishes on the collisions and on nothing
else.

## (f) Theorem E: the ROT3 rung lemma, which carries twelve of the twenty-three even-n pairs

At even n no member of the X ladder is a Galois image of ΔM₁: the multiplier n + 2j is even,
so gcd(n + 2j, 2n) ≥ 2 at every rung, and the collision hypothesis never reaches the X ladder
at all. Automorphisms still act WITHIN that ladder, carrying ΔX₂ to ΔX_2j whenever
gcd(j, 2n) = 1; what is missing at even n is the first step, from the collision to any X.
What carries the even-n vanishing instead is the shape of the triples.

**Definition.** Let 3 | n. A triple τ is *doubled-label ROT3* if its ±label set
P = {±k mod n : k ∈ τ}, taken mod n, has six elements and is the union of two cosets of the
order-3 subgroup ⟨n/3⟩ ≤ ℤ/n. This is Conway-Jones' ROT3 shape, F89_SEED_EXISTENCE_REDUCTION's family, read on the DOUBLED
labels. `Seed.cs`'s `TripleFamily.Rot3` is the same shape on a narrower domain, the triples
whose own level sum vanishes, and the two do not coincide here: (1,2,4) at n = 9 is
doubled-label ROT3 with a level sum of 3.7588, outside that classifier entirely.

**The two cosets are C and −C.** P is closed under negation and has six elements, none of
them 0. A coset fixed by negation would contain a fixed point of x ↦ −x, an involution on a set of odd
size 3 having one; and the fixed points in ℤ/n are 0 alone at odd n, 0 and n/2 at even n; 0 is not in P and n/2 lies in at most one coset, so negation cannot fix both cosets; and a permutation of a two-element set that fixes one fixes
the other, so it fixes neither and swaps them. Write them C = a + ⟨n/3⟩ and −C.

**Statement.** Let τ be parity-uniform (all labels of one parity, so η_k = η_τ throughout)
and doubled-label ROT3, with coset label a. Writing ω = ζ_n,

    2·X_2j(τ) = η_τ · Σ_{e∈P} ω^{j·e}.

Multiplication by j sends the coset a + ⟨n/3⟩ to j·a + ⟨j·n/3⟩. If 3 ∤ j then
⟨j·n/3⟩ = ⟨n/3⟩, the image is again a full coset of the order-3 subgroup, and a full coset
of a nontrivial cyclic group sums to zero. Both cosets do, so

    **X_2j(τ) = 0 for every j with 3 ∤ j, forced.**

If 3 | j then j·n/3 ≡ 0 mod n, the coset collapses to one point taken three times, and

    X_2j(τ) = 3·η_τ·cos(2πja/n),

which is a specific cosine, not a forced zero. **BREAK-INPUT, and it is not hypothetical:**
at n = 24 the triple (1, 7, 9) is clean, parity-uniform and ROT3 with cosets {1,9,17} and
{7,15,23}, so a = 1, and that cosine is cos(π/4) ≠ 0 at j = 3 and cos(π/2) = 0 at j = 6. So the
lemma has exactly one forced direction, and a version reading *X_2j = 0 if and only if
3 ∤ j* is false. It survives a census at n = 30, where the collapse cosine never vanishes, and dies at n = 24 and
at n = 12; n = 30 alone could not have broken it, which is the point. Gate L10 pins the forced direction on every parity-uniform ROT3 triple
at every 3 | n ≤ 30 and every rung j ≤ 9 with 3 ∤ j, pins the collapse formula exactly, and
exhibits both outcomes at 3 | j.

**Consequence, and what it does NOT give.** If both triples of a pair are parity-uniform
doubled-label ROT3, then ΔX₂ = ΔX₄ = 0 term by term, so

    c₃ = 0     and     c₁ = (2/n)·ΔX₀ = (4/n)·(o_τ − o_σ).

The first order is left standing by the shape alone: ΔX₀ is not constrained by it, and c₁
vanishes only when the two triples carry the same odd-label count, which at even n means the
same parity class. That is not a technicality. Fifty-eight collision pairs of this census satisfy the hypothesis in full, both triples
parity-uniform and doubled-label ROT3, and have c₁ ≠ 0, thirty of them at
odd n and twenty-eight at even n, so the shape fails to give c₁ = 0 on both sides of the parity.
At the smallest modulus that carries one, n = 9, the pair is (1,5,7) ~ (2,4,8): two zero-sum ROT3
triples of opposite parity class, ΔX₀ = 6, c₁ = 4/3, and it is the pair Confirmation 24 flew, named in the sweep record above for another reason. The smallest |c₁| over the 58 is 2/5, at n = 30, and the reason is that all 58 have |o_τ − o_σ|
= 3, so |c₁| = 12/n and the smallest is simply the largest modulus. The twelve that stand do share a parity class. The shape is what makes c₃
vanish; the parity match is what makes c₁ vanish.

**Where the hypothesis is free, and where it is not.** At an EVEN n parity-uniformity comes
with the shape and is not a second assumption: ROT3 needs 3 | n, an even n gives 6 | n, so
n/3 is even, every coset {a, a + n/3, a + 2n/3} is parity-homogeneous, and so is −C since n
is even. Gated at L9 on all 80 ROT3 triples at n = 12, 18, 24, 30. At an ODD n the hypothesis
is real and load-bearing: (1,2,4) at n = 9 is clean and doubled-label ROT3 with mixed parity,
and its X₂ is not zero.

**Census fact, READ.** All twelve even-n non-mirror standing pairs at n ≤ 30 (four at
n = 24, eight at n = 30) have both triples doubled-label ROT3 and share a parity class,
gated at L9. The lemma turns the shape into the vanishing of c₃ and the parity match turns it
into the vanishing of c₁; neither says a standing pair must have either, and nothing here
decides whether a pair of another shape can stand.

**The eleven Θ-mirror pairs, on the same ladder but NOT by Theorem E.** Their argument is the
one line below and uses no ROT3 hypothesis, which matters because one of the eleven, (6,18,20) ~
(10,12,24) at n = 30, is not doubled-label ROT3 at all. If σ = n − τ then at even n both
η_{n−k} = η_k and cos(2j(n−k)π/n) = cos(2j·kπ/n), so X_2j(σ) = X_2j(τ) term by term for
EVERY j. Through fifth order that is every odd coefficient of the gap, by Corollary C. The
all-orders statement needs a route that does not pass through the ladder, and the comb page
has one: K·H(u)·K = −H(−u) at even n makes the gap S_τ(u) + S_τ(−u), even in u by
construction. That is [PROOF_ZETA2_ANTI_PROTECTION](PROOF_ZETA2_ANTI_PROTECTION.md)'s protection and the comb
page's gate R4c. The two routes divide cleanly: the ladder reaches every odd coefficient through
the fifth WITHOUT the K identity, and the K identity is what carries it to all orders. They
agree where they overlap.
Gate L9 pins ΔX_2j = 0 for those eleven at every rung j ≤ 3, and L11 records that among the
223 standing pairs, ΔX₆ = 0 holds for exactly those eleven.

## (g) Theorem F: where the ladder stops, and why it is F129's own condition

F129 fires only at 3 | n ≥ 9 or 10 | n ≥ 20 (PROOF_F129_LEVEL_COLLISION_LAW, cited). Since
10 | n forces n even, the second clause admits no odd modulus at all, and **every ODD firing
modulus has 3 | n**. That is one line of arithmetic on F129's condition and it carries no
gate: a gate whose expectation is read off the same condition it tests cannot fail, and the
first version of this file carried exactly that gate.

The rungs of the X ladder are j = 0, 1, 2, …, and j = 0 is never an automorphism: gcd(n + 0, 2n)
= n. That rung is exactly the one that leaves c₁ standing, ΔX₀ being the odd-label count, and it
is why 427 of the 627 odd-n collision pairs separate at first order. Among the rungs j ≥ 1, at
odd n the first that is not an automorphism is j = 3, and at even n the first that collapses a
ROT3 coset is also j = 3. So for a pair that STANDS at first order, where the j = 0 rung is
empty too, X₆ is the first rung that can survive.
By Corollary C, X₆ first appears at FIFTH order, with coefficient
−(3/4)(2n−3)(n−3)·F, whose factors (2n−3), (n−3), (n−6) and (n−2) are all positive for
n ≥ 9, so it vanishes for no firing modulus; spot-gated at L3 for n = 9..60, with positivity
as the reason, and note that F itself does vanish at n = 6, outside the firing range.
Therefore for a pair that
stands, where ΔX₂ = ΔX₄ = 0,

    g(u) = c₂·u² + c₄·u⁴ + c₅·u⁵ + O(u⁶),      c₅ = −(3/4)(2n−3)(n−3)·F·ΔX₆,

so the ODD part of the gap begins at u⁵ exactly when ΔX₆ ≠ 0. That the exponent is 5 and not
7 is therefore not a theorem about all standing pairs but a decidable question about each
one, and the census decides it exactly: ΔX₆ ≠ 0 for all 212 non-mirror standing pairs (200 at odd n and 12 at even n) and = 0 for
exactly the eleven Θ-mirror ones (L11), at n ≤ 30. What IS a theorem is the vanishing below it, c₃ = 0, and the reason X₆ is the first rung that
can survive at all. One qualification, since it is easy to lose: at odd n that theorem is
unconditional, covering all 627 pairs; at even n it reaches the twelve through Theorem E, whose
hypothesis, the ROT3 shape, is a census observation and not derived, so for those twelve the
vanishing is a theorem ON a read hypothesis. The eleven Θ-mirror pairs need neither. Together they
account for THE_COMB_ON_THE_ROAD's measured *"the odd part starts at u⁵ (decade ratio within
15% of 1e5 for every one)"*, exactly rather than at 40 digits, and give the coefficient.

The reason deserves one sentence on its own. The multiplier that breaks the ladder is 3,
and 3 is the prime F129 needs in order to have a collision at an odd modulus at all. The
same arithmetic that lets the coincidence exist is the arithmetic that stops the road from
dissolving it any faster.

## (h) Corollary G: the second order, and the criterion that was too crude

c₂ = −((n−3)/n²)·ΔM₃, and the first reading of it is the one Theorem D gives: the M ladder's
criterion is gcd(3, 2n) = gcd(3, n) = 1, i.e. 3 ∤ n, which among the firing moduli is the
10 | n family with 3 ∤ n, smallest member n = 20. There k ↦ 3k is an automorphism,
ΔM₃ = σ₃(ΔM₁) = 0, and **c₂ = 0 for every collision pair**. Census, gated at L12: exactly 60
of the 2558 pairs have c₂ = 0, all 20 at n = 20 and 40 of the 1313 at n = 30.

The first version of this section called those 40 *met*, a coincidence at a modulus where the
kill does not apply. They are not. **The criterion gcd(m, 2n) = 1 is sufficient and too crude,
and the sharp one is local to the pieces.**

A collision reduces to a vanishing sum of roots of unity, twelve of them or eight once a pair's
shared labels cancel (F129's §2 reduction), and such a sum decomposes into MINIMAL vanishing
pieces, meaning vanishing subsets with no proper vanishing subset. The family vocabulary below
is F129's family inventory's, not this file's: a *piece* is one such subset, its *ratio-order*
is defined two paragraphs down, R_p is a full set of p-th roots summed to zero, the *zero mode*
is the self-antipodal pair at n/2, and a family's *door* is the divisibility of n that admits
it. What matters about a piece is NOT that it is a coset of a prime-order subgroup. That is
false in general and false here: the inventory's own family signatures carry pieces of ratio-order 30, 42, 66, 70 and 210, of
which 30 and 42 occur inside this census, and the Conway-Jones and Poonen-Rubinstein classifications exist
precisely to describe the minimal sums that are not prime cosets, which is why this file's own
Builds-on bullet says PENT is not one. What matters is the piece's RATIO-ORDER, and that is what
`f129_family_inventory.piece_decomposition` computes and reports.

Take a minimal piece S with exponent set E ⊂ ℤ/2n, fix e₀ ∈ E, and put g = gcd({e − e₀ : e ∈ E}
∪ {2n}) and o = 2n/g. Every exponent of S lies in the single coset e₀ + ⟨g⟩, so S = ζ^{e₀}·w
where w is a sum of powers of ζ^g, a primitive o-th root of unity; since ζ^{e₀} is a unit and S
= 0, also w = 0 in ℤ[ζ_o]. Now read the same formal sum with every exponent multiplied by m: it
is ζ^{m·e₀} times the sum w with ζ^g replaced by ζ^{mg}. That is again a primitive o-th root
exactly when gcd(m, o) = 1, and then the bracket is σ(w) for some σ ∈ Gal(ℚ(ζ_o)/ℚ), hence zero.
Only a piece whose RATIO-ORDER shares a factor with m can survive. So:

> **The criterion.** If SOME minimal tiling of the pair's vanishing sum has every
> RATIO-ORDER coprime to m, then reading the sum under the multiplier m leaves it zero.
> For ODD m that reads back as ΔM_m = 0; the pair's minus sign is carried by ζ^n = −1, so
> at even m the same reading gives the SUM of the two triples' evaluations and not their
> difference. Global invertibility, gcd(m, 2n) = 1, is the special case where every
> tiling qualifies, since every ratio-order divides 2n. The statement does not care which
> modulus the pieces are written in: a ratio-order is intrinsic to the piece.

At m = 3 the census says the converse holds as well: **c₂ = 0 exactly when some minimal tiling
is entirely 3-free, on all 2558 pairs** (L13). The wording matters and the first version of it
was wrong: a minimal tiling is NOT unique, so *no piece has a 3-divisible ratio-order* is not a
property of the sum at all. The pair (15,23,25) ~ (17,19,29) at n = 30 tiles both as family C's
3-free ((2,2),(5,5),(5,5)) and as ((6,30),(6,30)), which is not; twenty of the sixty do. *Some
tiling* is the form the sufficient direction actually needs, and it is decomposition-free. Gated
as a break-input at L13. And the pairs with no 3-divisible piece are not scattered: inside this census they are one named
family, F129's **family C**, zero mode plus two R₅ pieces, door 10 | n, whose committed
closed-form count 2(n − 10) gives 20 at n = 20 and 40 at n = 30, total 60. So the 40 are as
forced as the 20; what differs is only that at n = 20 the whole modulus is 3-free while at n =
30 it is the PIECES that are. Gated at L13, which enumerates the vanishing subsets exactly in ℤ[ζ_2n] and keeps the
inventory's greedy two-prime decomposition only as a cross-check; the two agree on every pair.

**C is not the only 3-free family, and that gives the file its one prediction past the census.**
Reading the same committed table against the same criterion, exactly two of its thirteen
families have every ratio-order coprime to 3: C, and **L**, the zero mode with a single piece of
ratio-order 70, whose door is 70 | n and whose count is the constant 20. Since 70 = 2·5·7
carries no 3, the sufficient direction proved above applies to L unchanged, so

#{pairs with c₂ = 0} = 2(n − 10)·[10 | n] + 20·[70 | n],

each term carrying its own family's door, without which the first would read −2 at n = 9. It gives 20 at n = 20 and 40 at n = 30, matching the census, and then 140 at n = 70 and 420 at n
= 210, both firing moduli and both outside it. Only the second tests anything local: 3 ∤ 70, so
at n = 70 global invertibility already gives c₂ = 0 for every collision there and 140 is simply
the total. At n = 210 the modulus carries a 3 and the count is a claim the global criterion
cannot make. Past n = 30 the equality is a LOWER BOUND carrying
the proved direction only: every C and L pair has c₂ = 0, and whether any other pair there does
is the converse, which this file checks at n ≤ 30 and nowhere else. At n
= 210 the modulus is not 3-free at all (3 | 210), so that case is the local criterion doing work
no global one can. Gated at L13 as a statement about the committed table; the counts themselves
are not checked here and are the natural next thing to check. THE_COMB_ON_THE_ROAD's gate R4 counts the same 60 through
40-digit decade ratios, by a route with no arithmetic in it: two routes, one number, and
neither written with the other in view.

The same gate settles the complement exactly: no pair that stands at first order has c₂ = 0,
so *every one of the 223 leaves at second order* stops being a 40-digit reading and becomes a
statement about ΔM₃ in ℤ[ζ_2n].

**Theorems D and E run on this one mechanism, applied to two different sums.** Theorem D
decomposes the PAIR's vanishing sum in ℤ[ζ_2n] and concludes ΔM_m = 0; Theorem E decomposes a
SINGLE triple's doubled-label sum in ℤ[ζ_n] and concludes X_2j(τ) = 0, which the collision
hypothesis does not supply and the shape does. The argument is the same in both, word for word
with 2n replaced by n: Theorem D is the case where m is invertible on every piece at once,
Theorem E the case of pieces of ratio-order 3, where *3 ∤ j* is exactly *the multiplier is
coprime to their ratio-order*. What looked like two mechanisms, one for odd n and one for even, is one mechanism
read at two resolutions. The non-uniqueness of a minimal tiling is therefore not a worry but the reason the criterion is
quantified the way it is: σ_m(W) = 0 is a property of W, the hypothesis asks only that ONE
tiling qualify, and both sides are then decomposition-free. What the census reads in the
inventory's single greedy tiling is the sufficient hypothesis in its easiest form; the
existential is the statement.

One consequence worth keeping, because it collapses two sections of this file into one
event. At an odd firing modulus 3 | n, so the multiplier 3 fails on the M ladder and the
multiplier n + 6 fails on the X ladder for the same reason, gcd(6, n) = 3. **The second order cannot be killed and the fifth cannot be killed, and it is the same prime
that fails in both.** Whether the fifth order actually is nonzero is then the per-pair question
§(g) fences, ΔX₆ ≠ 0; what is one fact seen from two sides is the failure of the kill, not the
nonvanishing.

Note the complementarity, which is the compact form of the whole file. At odd n the X ladder's
MIDDLE rungs are automorphisms, which empties c₃ while leaving c₁ on the j = 0 rung and c₅ on j
= 3, and the M ladder's rung m = 3, the one c₂ needs, is not an automorphism since 3 | n (others are:
at n = 9, gcd(5, 18) = 1, so ΔM₅ = 0 on every collision pair there). At even n it is the other way
round: no rung of the X ladder is reached from the collision, and the M ladder can fire. The
parity of n decides which ladder has automorphisms; it never kills a whole half of the series.

## (i) What the sign of u is, and what c₁ = c₃ = 0 therefore says

A reading, gated only in its algebra (L14), and worth having because without it the theorem
is bookkeeping.

Flipping one bond of a loop flips the product of the hoppings around the loop, and no
diagonal gauge undoes that: u > 0 and u < 0 are half a flux quantum apart. F160's own curve
says so in two lines. At u = +1, G = 2·sin k·(cos Nk − 1), whose interior zeros are the periodic comb 2πm/N; at u =
−1, G = 2·sin k·(cos Nk + 1), whose interior zeros are the ANTI-periodic comb (2m+1)π/N. Two caveats the sibling proof carries and this reading inherits. At u = ±1 every such zero is a
DOUBLE root, cos Nk = ±1 being an extremum. And the open interval misses a band-edge level: at u
= +1 that is k = 0 and, at even N, k = π; at u = −1 the point k = 0 is not a level at all and
the missing one is k = π at ODD N. G vanishes at k = 0 and k = π for free at every u, which is
why the interval is open. Those two combs are already a committed object:
[PROOF_RING_GAP_DOMINANCE](PROOF_RING_GAP_DOMINANCE.md) identifies them as the two
Jordan-Wigner fermion-parity sectors of the ring. What is added here is that the u axis carries them, one on each side of the chain. Note the
range: §(a) rests the level labels on simplicity at u ≥ 0, and the sibling proof fences its own
departure count the same way, so u = −1 is read here only through the curve G, which is a
polynomial identity at every real u, and not through any labelled level.

So the EVEN part of a level in u is its response to how strongly the loop is closed, and the
ODD part is its response to WHICH parity sector. Read back onto the theorem, and only for a pair that ALSO stands at first order, since c₁ is the
leading odd coefficient and is nonzero for the other 2335: **c₁ = c₃ = 0 says such a collision
cannot tell the ring's two parity sectors apart below fifth order**, and for the eleven Θ-mirror
pairs it can never tell them apart at all. A pair with c₁ ≠ 0 distinguishes them at first order,
and 2335 of the 2558 do. §(f)'s term-by-term identity X_2j(σ) = X_2j(τ) proves it through fifth order, by Corollary C,
and the K identity proves it at all orders (see §(f), where the two routes are divided and where
it is said why neither is Theorem E); K·H(u)·K = −H(−u) is itself the statement *energy
reflection and the π flux are the same operation here*.

## (j) The gates, in one table

| Gate | What it pins | Kind |
|------|--------------|------|
| L1 | det(x·I − H(u)) = U_N − u²U_{N−2} − 2u, N = 3..9 | exact symbolic, against the matrix |
| L2 | d₁, d₂, d₃, d₅ in closed form, symbolic in n, both signs of η | exact symbolic |
| L2b | ⟨ψ_k\|V\|ψ_l⟩ = a_k a_l(η_k + η_l) and V rank one per sector, exact at N = 4..9; the angle equation, d_m = η^m·(η-free), and the two symmetries of that equation that give the parity split and the evenness of D_m at every m, symbolic in n | exact, part concrete in N |
| L3 | the multiplier form of d₁, d₂, d₃ and the vanishing sum of d₅'s four coefficients, symbolic; d₅'s X₆ coefficient nonzero, spot-checked at n = 9..60 | exact, part concrete in n |
| L4 | the five-term series against 60-digit eigenvalues at four values of u: the two difference ratios of q = residual/u⁶ both sit within 0.002 of 1/10 and the finer is nearer | error model, eigensolver |
| L5 | the census 2558 / 2335 / 223, plus a regression pin of ΔX₀ − ΔX₂ against the committed `n_times_D_vec` (the same quantity retyped, not a second route) | exact, ℤ[ζ_2n] |
| L6 | odd n: ΔX₂ = ΔX₄ = 0 on the whole population, which is all 627 pairs, so c₃ = 0 on every one; and ΔX₀ against its integer form 2(o_τ − o_σ) | exact, ℤ[ζ_2n] |
| L7 | odd n: ΔX_2j = 0 at every coprime rung j ≤ 8; the gcd identity at every odd n from 9 to 59 | exact, ℤ[ζ_2n] |
| L8 | NEGATIVE CONTROL: ΔX₂ = 0 exactly on the collisions; every pair at n = 9 and 15, 40000 sampled at each of n = 21 and 27 | exact, ℤ[ζ_2n] |
| L9 | even n: 23 standing, 11 Θ-mirror with ΔX_2j = 0, the 12 both-triples ROT3 and parity-matched, and parity-uniformity forced on all 80 ROT3 triples | exact, ℤ[ζ_2n] |
| L10 | the ROT3 rung lemma forced at 3 ∤ j; the collapse formula; the BREAK-INPUT at 3 \| j | exact, ℤ[ζ_2n] |
| L11 | ΔX₆ = 0 on exactly the 11 Θ-mirror standing pairs, so ΔX₆ ≠ 0 on the other 212 | exact, ℤ[ζ_2n] |
| L12 | ΔM₃ = 0 on exactly 60 pairs, all 20 at n = 20 and 40 at n = 30, and on no standing pair (the formula itself is pinned at L3) | exact, ℤ[ζ_2n] |
| L13 | c₂ = 0 exactly when SOME minimal tiling is entirely 3-free, on all 2558, the tiling
lattice enumerated exactly and the non-uniqueness pinned as a break-input; the 60 are F129's
family C at 2(n−10); C is not alone, family L being 3-free too | exact, ℤ[ζ_2n] |
| L14 | u = +1 and u = −1 give the periodic and anti-periodic ring combs, and G vanishes at all 29 + 34 interior points of both, N = 4..12 | exact, part concrete in N |

## (k) Scope and fences

- **The expansion point is u = 0, the chain end.** Like F160's Theorem G and unlike its
  Theorem E, which expands in δ = 1 − u about the ring end and is about a pair SPLIT. The
  two are different objects on different combs (modulus N there, n = N + 1 here) and
  nothing in this file transfers to that one.
- **Five orders, not a general form.** Theorem A gives d₁ to d₅. The PARITY of the multipliers is
a theorem at every m (§(d), from the chiral K); what is NOT proved is the RANGE, that an odd
d_{2m+1} reaches no further than X_{2m+2}, which is only visible in d₁, d₃ and d₅, and the fact
that d₃ and d₅ carry no X₀ while d₁ does, which is read off the two cases. The ladder's general
shape is open in that sense.
- **Theorem E has one forced direction only.** 3 ∤ j forces the vanishing; 3 | j does not
  force anything, and n = 24 breaks the converse. A reading that treats the lemma as an equivalence is wrong: n = 24 and n = 12 both break it, and
this file's own first version made that error against n = 30, a modulus that could not.
- **The ROT3 shape of the twelve is READ, at n ≤ 30.** The consequence is a theorem, the
  hypothesis is a census observation. Whether a standing pair at even n must be ROT3, and
  whether a mixed-parity pair can stand at all, are the comb page's open items and stay
  open.
- **Every count on this page is at n ≤ 30**, the census's range, except those read off the
committed family table (the thirteen families, the ratio-orders 30, 42, 66, 70 and 210, and
§(h)'s 140 and 420): 2558, 2335, 223, 627, 427, 23,
11, 212, 60 with its 20 and 40, 12 with its 4 and 8, 1313, the 58 counterexamples of §(f) with
their 30 and 28, the 80 ROT3 triples, and L8's populations. Those exceptions are read off the committed family table; everything else in the list is
measured by the gate named beside it in §(j). The theorems (D at every odd n, E at every
  3 | n, F at every firing n) do not depend on that bound; no count says anything about a
  larger modulus.
- **The piece criterion has one proved direction.** *Some tiling with every ratio-order coprime to
m ⟹ the reading is zero* is a proof. The converse at m = 3 is a CHECK, not a proof, and it is
exact: L13 enumerates the vanishing subsets of each pair in ℤ[ζ_2n], takes the minimal ones and
asks whether an exact cover by 3-free pieces exists, on all 2558 pairs, with the inventory's
greedy two-prime tiling kept beside it as a cross-check. Nothing here proves the converse at any
other m or past n = 30.
- **§(i) is a reading.** Its algebra is gated (L14) and the two combs are
  PROOF_RING_GAP_DOMINANCE's; the sentence about what a level collision can and cannot
  distinguish is an interpretation of the theorem, not a second theorem.
- **No hardware claim.** Confirmation 24 is cited as a sighting of the n = 9 collision on
  the comb, on the ζ axis. Nothing here is confirmed by it and nothing here proposes a
  flight.
- **The letter n.** Comb modulus throughout, n = N + 1. See Notation.
