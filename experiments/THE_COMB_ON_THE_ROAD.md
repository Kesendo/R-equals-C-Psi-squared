# The Comb on the Road: F129's collisions under F160's crack

**Date:** 2026-09-02. **Authors:** Thomas Wicht, Claude (Fable 5.1). **Arc:** `the_forced_and_the_met`,
NextStep (1), *"RUN THE COMB TEST: take the cheapest comb law and detune the comb, solving rather than
sampling."* **Instrument:** [F160](../docs/ANALYTICAL_FORMULAS.md), the cracked ring's exact curve, whose
road u = J′/J on one bond leaves the open chain's cosine comb at u = 0 and reaches the ring's at u = 1
([PROOF_CRACKED_RING_EXACT_CURVE](../docs/proofs/PROOF_CRACKED_RING_EXACT_CURVE.md)). **Law under
test:** [F129](../docs/ANALYTICAL_FORMULAS.md), the level-collision law, the cheapest of the five comb
laws the arc names (F145/F146 counted as one) and the only one of them that transports to a detuned
spectrum at all (below). **Gate:** [`simulations/comb_road_f129.py`](../simulations/comb_road_f129.py),
nine gates (R1 counts two) and one reading, run committed at
[`comb_road_f129_run.txt`](../simulations/results/comb_road/comb_road_f129_run.txt).

**Notation, once.** N is the ring's site count and n = N + 1 is F129's modulus, so the chain comb at
u = 0 is E_k = 2cos(kπ/n), k = 1..N, in units of J. S(τ) = Σ cos(k_iπ/n) is F129's level of a triple, a
sum of three cosines; S_τ(u) = Σ_{k∈τ} E_k(u) is the road's sum of the three LEVELS, so S_τ(0) = 2·S(τ),
and every gap and every table entry below is in the level book, twice F129's. ζ_2n is a primitive 2n-th
root of unity; the ζ of the ζ² anti-protection law, quoted below, is that law's ZZ knob and a different
letter's worth of object. D(τ, σ) is the first-order gap defined below, not the proof's path
determinant D_n; D carries a 1/n, so n·D is the integer object in ℤ[ζ_2n]. α_k/γ₀ is F65's endpoint
dissipation rate of chain level k in units of γ₀, a γ-free number, (4/n)·sin²(kπ/n). *Velocity* here is
dE/du along the road, a parameter derivative, not F2b's group velocity dE/dk on the same chain. *Stands*
means D = 0 on the u axis; Confirmation 24's *standing fringe* is the ζ axis of the Floquet step, where
the flown pair stands, and here it separates.

**In one paragraph.** Every one of the 2558 exact collision pairs of F129 at the nine firing moduli
n ≤ 30 separates as the road leaves the comb. 2335 separate at first order in u, and that is a theorem
per pair: the first-order gap is an element of ℚ(ζ_2n), decided exactly and found nonzero. The other 223 stand at first order (the same element is exactly zero) and every one of them leaves
at second order; that was read here at 40 digits with a decade law and is exact since
[F161](../docs/ANALYTICAL_FORMULAS.md), which writes c₂ = −((n−3)/n²)·ΔM₃ and finds ΔM₃ ≠ 0 on
every one of the 223 in ℤ[ζ_2n]. At u = 1/2 no former collision is closed, to forty digits.
So the arc's prediction, *the comb family dissolves together*, holds for this law, and it holds with more
structure than the prediction had. The first-order motion of every chain level along the road is the
signed endpoint rate comb of F65, dE_k/du = (−1)^{k+1}·α_k/γ₀ (Theorem G of the proof), so a
collision stands at first order exactly when the SIGNED F65 sums of its two triples agree. At odd n
that condition reduces, by a Galois automorphism the collision itself is invariant under, to one
integer: the two triples carry the same number of odd labels, and the separation speed takes only the
values 0, 4/n, 8/n, 12/n. At even n every pair that stands, in this census, is a collision between two
ZERO-SUM triples whose doubled sums agree as well, 23 of them, and 11 of those are Θ-mirror pairs
σ = n − τ (the ζ² proof's word), for which the gap is even in u to all orders by the chiral K, the same
linear half of the mirror that protects the flown F129 pair at first order in ZZ (the ζ² anti-protection
law). A Θ-mirror collision pair stands at first order exactly when n is even, which is why the flown
n = 9 pair separates on the road although it is protected against ZZ. Each collision is a coincidence with no forcing yet found; the rate at which the road dissolves
it is forced, and since F161 so is the order at which that forcing runs out.

## What the repo already holds

The sweep was run 2026-09-02 by three agents over the markdown layer, the typed layer and the five
comb laws, store by store, and re-run by two record auditors and a reader whose corrections are folded
in: a first version of this section reported `fw.Confirmations` empty when it holds the flown
confirmation on this census's own first row, listed the repo's bond detunings as exhaustive while
omitting F130's, and did not know that the census's first column is a committed closed form.

- **[`docs/ANALYTICAL_FORMULAS.md`](../docs/ANALYTICAL_FORMULAS.md).** F129 states the law: for
  distinct CLEAN triples τ, σ ⊂ {1..n−1} (no internal pair summing to n), S(τ) = S(σ) with
  S(τ) = Σ cos(k_iπ/n) forces 3|n or 10|n, with converses at every 3|n ≥ 9 and 10|n ≥ 20; away from
  both families the level map is injective on clean triples, certified to n = 210; the collision
  COUNTS per firing n are the thirteen closed forms of F129's family inventory
  ([F129_FAMILY_INVENTORY](F129_FAMILY_INVENTORY.md), typed as `CollisionFamilyInventory`, gate
  `f129_family_inventory.py`), and the census column below reproduces them row for row (run with
  `--fast` on 2026-09-02: 1, 25, 127, 162, 20, 255, 411, 244, 1313). F160 holds the road, and its
  sweep record carried, until this page, the sentence *"none of those has been run on this road"*; the
  same change that lands this page rewrites it. F65 holds the endpoint rates
  α_k/γ₀ = (4/(N+1))·sin²(kπ/(N+1)) for the dephased endpoint, which enter this page not as a law
  under test but as the road's velocity (the phrase *rate comb* for that set is this page's, F65 says
  *the dissipation rates*); its Niven face already writes them on the DOUBLED angle,
  (2/(N+1))·(1 − cos(2kπ/(N+1))), which is the angle this page's even-n coincidence lives on. F75 names
  the sign: the mirror sign η = (−1)^{k+1}. F123 is the Re-side sibling of the velocity in the MOVE and
  not in the form: the same Hellmann-Feynman derivative in one bond's knob, there giving the squared
  difference of the mode's site occupations across that bond, (n(j) − n(j+1))², which vanishes at the
  ends where Theorem G's signed endpoint product lives. F89's resonance count has no registry formula
  line of its own (it lives in F89_SEED_EXISTENCE_REDUCTION Pieces 2-3); F145/F146 are the spin-1
  triplet and the scalar count C(⌊N/2⌋, ℓ)·R_ℓ; F144 is the disagreement floor ℓ(N−ℓ)/(N+1) with its
  N = 5 exception. The formula names for the velocity's magnitude are F64's |a_B|² and F65's α_k; the
  Absorption Theorem supplies the light_l(v) they are instances of.
- **[`docs/proofs/`](../docs/proofs/).** [PROOF_F129_LEVEL_COLLISION_LAW](../docs/proofs/PROOF_F129_LEVEL_COLLISION_LAW.md)
  §1 defines the level as a sum of three cosines and §2 reduces a collision to a vanishing sum of comb
  cosines, then of roots of unity; its mechanism (Lam-Leung, Poonen-Rubinstein) has no input off the
  comb, so what this page tests is the empirical face, injectivity.
  [PROOF_ZETA2_ANTI_PROTECTION](../docs/proofs/PROOF_ZETA2_ANTI_PROTECTION.md) (F131 Theorem B) holds
  F129's own n = 9 pair (1,5,7) ~ (2,4,8) under site-dependent ZZ and Z-detuning to all orders: the
  Θ-mirror pair, ν = 9 − τ, is protected at first order and pushed apart at second, *"twice one
  branch's shift"*. That is an F129 collision perturbed off the comb, by site fields and not by a bond,
  and its mechanism reappears on this road (below), with its two fences quoted where it does.
  [PROOF_SCALAR_COUNT](../docs/proofs/PROOF_SCALAR_COUNT.md) owns F145/F146 and
  [PROOF_FROZEN_BAND_SO4](../docs/proofs/PROOF_FROZEN_BAND_SO4.md) §7 owns F144; neither takes a level
  list as input (below). [PROOF_ABSORPTION_THEOREM](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) names
  the velocity's magnitude up to a factor two: the light at a site, light_l(v) = Σ_x Δ_l(x)|v_x|²/Σ|v_x|²,
  is F64's |a_B|² at the dephased endpoint, and the velocity's magnitude α_k/γ₀ = 2|a_B|² is twice it.
  [PROOF_F92_BOND_ANTI_PALINDROMIC_J](../docs/proofs/PROOF_F92_BOND_ANTI_PALINDROMIC_J.md) and
  [PROOF_F93_DETUNING_ANTI_PALINDROMIC](../docs/proofs/PROOF_F93_DETUNING_ANTI_PALINDROMIC.md) own
  bond and field detunings of the chain on the Liouvillian side (the parameter Klein), not of the comb's
  collisions. [PROOF_DEPHASING_FRONT_RENEWAL](../docs/proofs/PROOF_DEPHASING_FRONT_RENEWAL.md) stands
    on the cracked ring since 2026-08-31 (its Status line), a proof on the object and not about
  it. [PROOF_K_PARTNERSHIP](../docs/proofs/PROOF_K_PARTNERSHIP.md) Lemma 1 owns KHK = −H for
  bipartite nearest-neighbour hopping with non-uniform J, the identity the mirror-pair
  argument below stands on, and holds it under a bond detuning: a symmetry of H, not one of
  the five laws' arithmetic. No proof holds any of the five laws' arithmetic under a BOND
  detuning of the comb.
- **[`experiments/`](.)**, including the null results and the prior flights.
  [IBM_F129_RAMSEY_FRINGE](IBM_F129_RAMSEY_FRINGE.md) flew the n = 9 pair on 2026-07-15 (Confirmation
  24, below) with single-qubit Z-detuning and ZZ as its perturbations: F129's levels detuned by fields,
  at the flown pair only. The nearest bond detuning of a collision object is F130's, not one of the
  five: `f130_time_domain_decoupling.py` gate G4, *"the broken comb: a defect bond (delta = 0.15) destroys
  the eigenmultiplet property"*, on a level-collision multiplet whose pair 2 + 10 = 12 = n is not clean
  (`simulations/endpoint_density_gate.py` also builds chains from per-bond vectors and reads level
  coincidences of ad_H eigenspaces off them, a different object again, no F129 in it).
  Every other detuning in the tree is on other laws or other objects: F160/Crack/Warble, F157/BlindSeat, the two
  gap-dominance Scope fences, F155/F112/F113, XY_FROZEN_BAND, F92/F93, and F142's ladder break-inputs
  (THE_EXCEPTIONAL_COUPLINGS, `eta_ladder_breakinput.py`). Three near-misses, so they are not mistaken
  for prior work: `eta_ceiling_reduction.py`'s `h_chain` accepts a bond profile but every call with one
  sits in F142's ladder block, never in F144's floor blocks V10-V12; PROOF_FROZEN_BAND_SO4's *"all but
  finitely many couplings"* is the scalar J of a pencil, the comb untouched; and F89_SEED_EXISTENCE_REDUCTION
  discusses the RING for its count and records it as a break (*"a ring would let particles pass around
  the wrap"*; `Seed.cs`: *"on the RING it is an artifact"*), a change of graph rather than a walk along
  u. **None of the five laws had been evaluated at any u ≠ 0, 1**, as of HEAD `1a65a9c`. Two neighbours
  of Theorem G's MOVE, found by the reader: F89_SEED_EXISTENCE_REDUCTION Step 4 already runs a Galois
  automorphism σ_j: ζ ↦ ζ^j on this comb (as a norm squeeze; the single multiplier k ↦ k(n+2) below is
  new, the technique is not), and F86c's per-bond Hellmann-Feynman response K_b (PROOF_F100 §Q_peak)
  is the same derivative in a bond's knob on a Liouvillian observable. THE_BLIND_SITE §5 holds the
  rank-one member of the crack's determinant family; the crack is the rank-two member.
- **[`hypotheses/`](../hypotheses/).** [PERSPECTIVAL_TIME_FIELD](../hypotheses/PERSPECTIVAL_TIME_FIELD.md)
  §3.2   (Tier 2) carries Theorem G's eigenvector route one bond over: *"A = 2 ψ_k(b) ψ_k(b+1) is the (real)
  bond overlap of the single-excitation mode at the defect bond (b, b+1)."* What Theorem G adds is the closed
  form at the wrap bond, the product of the two END amplitudes being the signed F65 rate.
- **[`fw.Confirmations`](../simulations/framework/confirmations.py) and the C# `ConfirmationsRegistry`.**
  `f129_standing_fringe_kingston_july2026` (Confirmation 24, `ibm_kingston`, 2026-07-15, job
  `d9br4vmg26ic73dgbgk0`): *"the exact clean-clean level collision (1,5,7)~(2,4,8) of the n = 9 comb"*,
  flown as a two-branch 3-magnon cat and found standing to its computed drift. That pair is the single
  exact collision the census below finds at n = 9, so the road's first row has a hardware sighting, on
  the comb. Nothing ring-shaped.
- **The OpenArcs registry** (`compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs`),
  `the_forced_and_the_met`. NextStep (1) is quoted in the header; its own fence, *"THE_SEAT_THAT_CUTS's
  own open item asks for a DETUNED BOND rather than a Delta, so do not report one as the other"*, is
  respected: this page's bond is the wrap bond of a ring under F129's level sums, not a blind seat.
  The arc names the laws slightly more narrowly than "F65" and "F144": *F65's Niven root* and *F144's
  N = 5 exception*.   Of its three method words, one is used here in its sense: the balanced-pair check (R2) is what
  the arc calls a NEGATIVE CONTROL (a configuration where the law must not fire); deciding n·D in
  ℤ[ζ_2n] instead of sampling u is the arc's own first trap avoided, the F157 Δ-locus move on the u axis.
- **The code.** `compute/MirrorWorld/LevelCollision.cs` (`CensusOf(int n)`), `Seed.cs`
  (`VanishingTriples()`), `simulations/scalar_count.py` (`energy_classes(N, l)`),
  `simulations/f129_level_collision_law.py` (`level_vec(n, t)`): every one builds the comb from n
  internally and none accepts a level array. `Crack.Levels(points)` returns the road's N levels at any
  rational u with no eigensolver; F129's census and the crack's road had never been composed. This
  page reuses the F129 script's exact cyclotomic layer (`root_sum_vec`, `level_vec`,
  `collision_groups_mod_p`, integer vectors in ℤ[x]/Φ_2n) rather than restating it, and
  `Cyclotomy.PathRateOrders` is the doubled angle as an integer object. `compute/RCPsiSquared.Core/Numerics/LevelCollisionCensus.cs` is the
  typed census of F129's collisions by the same exact route, live to n ≤ 60; like the others it builds
  the comb from n.   Two committed scripts already write Theorem G's Hellmann-Feynman element on the
  chain's sine modes, one bond in: `simulations/handshake_rk_first_principles.py` (*"d eps_1/dJ_b =
  <psi_1|V_b|psi_1> = 2 psi_1(b)psi_1(b+1)"*) and `simulations/handshake_rk_block.py`; the wrap-bond
  form and its identity with F65's rate are new, the move is not.
- **[`docs/GLOSSARY.md`](../docs/GLOSSARY.md) and the words.** No *comb*, no *crack*. *Road* is spent in
  F89's arc in a path sense F160's sweep record judged shared rather than colliding. *Survivor* is
  MirrorWorld's `Survivor`, the slowest mode, and F123's noun, so this page says *a pair that stands at
  first order* and *left to second order* instead (the verb *survive* stays a plain verb here, as in F65's
  and F160's entries). *Quantization condition* is this road's own name for G(k) = 0, so the discrete set
  of speeds below is not called quantized; *lattice* is MirrorWorld's `Lattice`, the bridged lattice of
  worlds, and the README's word for `Cyclotomy`, so it is not used for that set either. *Octave* is a
  factor-two bin in δ on the sibling page, so it is not used here; the page says *the doubled angle* and
  *twice the modulus*. *Field* is MirrorWorld's empty world and *pair* MirrorWorld's bare coherence; here a
  field is a site term in H and a pair is two triples. *Θ-mirror pair* is
  PROOF_ZETA2_ANTI_PROTECTION's term for σ = n − τ and is used in its sense, always with the Θ: F75's
  *mirror pair* is two SITES (ℓ, N−1−ℓ), a different object.
- **[`docs/CAUGHT_ERRORS.md`](../docs/CAUGHT_ERRORS.md).** Of the eight entries of 2026-08-31 and
  2026-09-01, four are the crack's own rounds (the landing round, the exactness round, the five review
  rounds, the F160 mint); nothing on the comb laws. The round on this page is the FIRST entry of 2026-09-02; the second, from the round that
landed F161, carries this page's own repairs.

## Which of the five laws can be run on a detuned spectrum at all

The arc's cut says the five laws *"all stand on ONE object, the spectrum {2cos(k*pi/(N+1))}"*. That is
true of what they are about and not of what they CONSUME, and the difference decides what "run it on
the road" can mean.

| Law | What it consumes | Transports to an arbitrary level list? |
|-----|------------------|----------------------------------------|
| F129 | levels as sums of three, the index condition CLEAN (k_i + k_j ≠ n), a conclusion in n (roots of unity) | the equality predicate yes, with CLEAN carried by the label k; the mechanism no |
| F89's resonance count | the predicates λ_a + λ_b = λ_c and λ_a + λ_b + λ_c = 0, tied to a NULLITY whose proof needs the open chain's ordering sectors | the predicates yes; the nullity identity no, and the road is a ring at every u > 0 |
| F145 | the chiral involution ā = M − a and the ladder coefficients 1 and 2 | no: there is no M on the road |
| F146 | Slater-energy coincidence classes (yes) against C(⌊N/2⌋, ℓ)·R_ℓ (index arithmetic) | the hypothesis yes, the conclusion no |
| F65's Niven root | an eigenvector amplitude at the dephased ENDPOINT | no: the ring has no endpoint (the rationality question about the road's own values would be well posed, and is not this page's) |
| F144's N = 5 exception | ℓ(N−ℓ) = N+1, Diophantine in (N, ℓ), and the chiral pairing | no |

So the cheapest comb law is F129 and it is the only one of the five whose central object, a
coincidence among sums of levels, moves to a detuned spectrum with nothing but a labelling
convention. The convention is: label the road's N levels k = 1..N by descending energy, which is the
chain comb's own order at u = 0. It agrees with the continuation of the labels along the whole road
0 ≤ u < 1 because no two levels ever cross there, and since 2026-09-02 that is a theorem, not the
measurement the sibling page had fenced (*"nothing here proved it for every N and δ"*, E2d(i)'s sweep
over 1050 (N, δ) points): the curve factors into the two reflection sectors through the crack,
G = 2AB, each factor a nonvanishing prefactor times the characteristic polynomial of an unreduced
Jacobi block, and a common zero of A and B forces u² = 1 (Corollary B of the proof, its simplicity
clause, found by a review lens on this page's first version and repaired by another). Everything
decided exactly on this page is decided at u = 0 in any case, where the labels are the comb's. F65 did
enter, from the other side: it is not run on the road, the road runs on it.

## The velocity, and what decides first order

Theorem G of the proof, restated for this page. Along the road H(u) = H_chain + u·V with
V = |0⟩⟨N−1| + |N−1⟩⟨0|, the chain spectrum is simple, and

    dE_k/du |_{u=0} = ⟨ψ_k|V|ψ_k⟩ = 2ψ_k(0)ψ_k(N−1) = (−1)^{k+1}·(4/n)·sin²(kπ/n),        n = N + 1.

In the repo's own words: the magnitude is F65's α_k/γ₀, twice the Absorption Theorem's light at the
dephased endpoint (light_B(v) = F64's |a_B|² = (2/n)·sin²(kπ/n)), and the sign is F75's mirror sign
η = (−1)^{k+1}, F71's reflection parity of the mode, not the chiral K, which swaps k ↔ n − k rather
than fixing k. The dissipator prices the end amplitude squared; the wrap bond prices it times the
other end; on a reflection-symmetric chain the two differ by the mirror's sign and nothing else. So
for a collision S(τ) = S(σ) on the comb, the first-order gap on the road is u·D with

    D(τ, σ) = Σ_τ (−1)^{k+1}·w_k − Σ_σ (−1)^{k+1}·w_k,        w_k = (4/n)·sin²(kπ/n) = (2/n)·(1 − cos(2kπ/n)),

an element of ℚ(ζ_2n); n·D lies in ℤ[ζ_2n] and is decided EXACTLY in the F129 script's own layer
(integer vectors modulo Φ_2n), and the two verdicts have different standing: **D ≠ 0 is a theorem** that the pair separates as
the road leaves the comb (an analytic gap with a nonzero derivative is nonzero on a punctured
neighbourhood of u = 0); **D = 0 leaves the pair to second order**, where this page reads rather than derives; F161 derives
it afterwards, and the reading held. Gates R1a (sympy, n = 5, 6, 7, 9) and R1b (finite differences at every firing n, worst
1.1e-9) hold the velocity; the proof's gate P7 has it from the polynomial's derivative as well, exactly.

**The negative control the CLEAN condition points at (gate R2).** A balanced pair k + k′ = n, which
CLEAN excludes, sums to zero on the comb. Its first-order derivative on the road is
w_k·[(−1)^{k+1} + (−1)^{n−k+1}]: exactly 0 at odd n (even N, where the cracked ring is still bipartite
and E ↔ −E holds at every u, so the pair's sum stays zero for the whole road) and nonzero at even n
(odd N, where the crack breaks bipartiteness). Exact at every n ≤ 30. The labels' chiral pairing is a
symmetry-forced coincidence at even N and a met one at odd N, and the road tells them apart at first
order, before any collision is looked at. This is the arc's negative control: the configuration where
the law must not fire, and does not.

## The census (gates R3, R3b, R3c)

Every exact collision pair at every firing n ≤ 30, enumerated by the F129 script's mod-p groups and
refuted exactly in ℤ[x]/Φ_2n; D decided exactly for each.

| n | N | exact collision pairs | separate at first order (D ≠ 0, a theorem) | stand at first order (D = 0), left to second |
|---|---|---|---|---|
| 9 | 8 | 1 | 1 | 0 |
| 12 | 11 | 25 | 24 | 1 |
| 15 | 14 | 127 | 107 | 20 |
| 18 | 17 | 162 | 160 | 2 |
| 20 | 19 | 20 | 20 | 0 |
| 21 | 20 | 255 | 195 | 60 |
| 24 | 23 | 411 | 404 | 7 |
| 27 | 26 | 244 | 124 | 120 |
| 30 | 29 | 1313 | 1300 | 13 |
| total | | 2558 | 2335 | 223 |

The n = 9 row is the flown pair, (1,5,7) ~ (2,4,8), Confirmation 24. On the road it separates at first
order: τ is all-odd and σ all-even, so the signs make the difference a sum,
D = (2/9)·[Σ_τ(1 − cos(2kπ/9)) + Σ_σ(1 − cos(2kπ/9))] = (2/9)·(3 + 3) = 4/3, both doubled sums being the
n = 9 zero sums themselves.

**At odd n the pairs that stand obey a theorem (gate R3b).** (−1)^{k+1}·cos(2kπ/n) = −cos(k(n+2)π/n),
and for odd n the map k ↦ k(n+2) mod 2n is a Galois automorphism of ℚ(ζ_2n) (gcd(n+2, 2n) = 1). So
the signed cosine part of D is minus the Galois conjugate of S(τ) − S(σ) = 0 and vanishes for EVERY
collision; what is left of D is its constant part,

    D(τ, σ) = (4/n)·(o_τ − o_σ),        o = the number of ODD labels in the triple,

which holds exactly on all 627 collision pairs at n = 9, 15, 21, 27, both as an identity of
cyclotomic vectors and as the set equality *the pairs with D = 0 are the equal-odd-count pairs*, both
inclusions. Two things follow. A collision at odd n stands at first order exactly when its two
triples carry the same number of odd labels, an integer condition with no cosine in it. And the
separation speed can only take the values 0, 4/n, 8/n, 12/n, in units of J per unit of u, and all
four occur (|o_τ − o_σ| over the pairs: n = 9 {3: 1}; n = 15 {0: 20, 1: 73, 2: 12, 3: 22}; n = 21
{0: 60, 1: 72, 2: 60, 3: 63}; n = 27 {0: 120, 3: 124}, so at n = 27 every collision either stands or
separates at the full 12/n): the collision is a coincidence with no forcing yet found, and the rate at which the road
dissolves it is forced, by the same Galois structure that made the collision possible. The odd-label
count is the counting form of the character k ↦ (−1)^k, which the chiral involution k ↦ n − k (F145's
ā = M − a) flips at odd n and preserves at even n; that flip is what freezes R2's balanced pair at odd n,
and the same character is what R3b counts.

**At even n the pairs that stand are read, and their vanishing is then derived (gate R3c).** k ↦ k(n+2) is not an
automorphism there (n + 2 is even). All 23 pairs found (n = 12: 1; 18: 2; 24: 7; 30: 13) are
collisions between two ZERO-SUM triples, S(τ) = S(σ) = 0 exactly, which is F89's resonance object
(`Seed.cs`'s `VanishingTriples()`): the count F89 builds on it does not transport, the object did.
For the 11 Θ-mirror pairs among them the zero sum is forced, not met: at u = 0 the chain is
bipartite, E_{n−k} = −E_k, so S_{n−τ} = −S_τ and a collision S_τ = S_{n−τ} is a zero sum (the same
holds for the 11 Θ-mirror collision pairs at odd n, printed by gate R3, none of which stands). For the
other 12 it is read. Each of the 23 is parity-uniform of one class (both triples all-odd or both
all-even), where D reduces to ±(2/n)·(Σ_σ cos(2kπ/n) − Σ_τ cos(2kπ/n)), so D = 0 exactly when the
DOUBLED sums agree: for the 12 non-mirror pairs a coincidence on F65's own doubled angle, at twice the modulus, whose
SHAPE is named by [F161](../docs/ANALYTICAL_FORMULAS.md): both triples are Conway-Jones ROT3 in
the doubled labels, a union of two order-3 cosets, and once the shape is there the vanishing is
forced at every rung not divisible by 3, the quadrupled angle included, so what is met here is the
shape and not each rung separately; for the 11 Θ-mirror pairs one line, since w_{n−k} = w_k
and at even n (−1)^{(n−k)+1} = (−1)^{k+1}, so the signed weight is invariant under k ↦ n − k and
D(τ, n−τ) = 0 term by term. In 22 of the
23 the doubled sums are both zero (the n = 18 pairs (2,10,14) ~ (4,8,16) and (1,11,13) ~ (5,7,17) are
the n = 9 zero-sum triples read at twice the modulus); in one, (6,18,20) ~ (10,12,24) at n = 30, they
are equal and nonzero. Eleven of the 23 are Θ-mirror pairs, σ = n − τ, the ζ² proof's kind. Whether a
mixed-parity pair can stand at even n is not decided here: none occurred to n = 30, and the page says
no more than that.

**The Θ-mirror pairs' parity corollary, in one line.** The signed weight (−1)^{k+1}·w_k is invariant
under k ↦ n − k exactly when n is even (w_{n−k} = w_k always; the sign flips at odd n), so for a
Θ-mirror pair σ = n − τ the first-order gap is D = 0 at even n and D = 2·Σ_τ(−1)^{k+1}w_k = (4/n)·(2o_τ − 3) ≠ 0
at odd n (the Galois form of the same statement, o_σ = 3 − o_τ there). A Θ-mirror collision pair
therefore stands at first order exactly when n is even, and the chiral K below is needed not for this
but for what only it gives, the evenness of the gap in u to ALL orders. The flown pair lives at n = 9
and separates at 4/3 per unit of u.

## Second order, and the finite-u reading (gates R4, R4b, R4c, R5)

**The derivative is a law, not a number (R4), on every pair.** For all 2335 pairs with D ≠ 0 the
read gap [S_τ(u) − S_σ(u)]/u, from 40-digit eigenvalues at u = 1e-2, 1e-3, 1e-4, converges to the
exact D with a residual whose leading power is 1 or 2: at the finer decade pair (1e-3 → 1e-4) the
ratio is within 15% of 10 for 2275 pairs and within 15% of 100 for 60 (the pairs whose second-order
coefficient vanishes), and the finer pair's ratio is never farther from a pure power than the coarser
one's (a first version of this sentence, of the gate's label and of the script's docstring had the
inequality the other way round, while the code had it right). Four pairs at n = 30 show why the class has to be read at the finer pair, and what they really are
is exact: their ΔM₃ vanishes in ℤ[ζ_2n], so by [F161](../docs/ANALYTICAL_FORMULAS.md) their c₂ is
not small but ZERO, and what still competes at u = 1e-2 is the u³ term against the u⁴ one, which
is why the coarse ratio sits between 56 and 140 instead of at 100. Three earlier
versions of this gate were wrong in three ways the record now carries: one selected its pairs with a
float `D == 0.0` and let exact zeros through as separating (five of the 40 pairs then sampled; on the
full population 173 of the 223 exact zeros have a nonzero float D, a count gate R3 now prints), while
the exact test stood three lines above; one demanded the class at the coarse pair and failed on those four; one tested the
convergence with the inequality inverted. The law is *the residual is O(u)*, and O(u²) is O(u).

**Every pair that stands leaves at second order (R4b).** For all 223, [S_τ(u) − S_σ(u)]/u² converges
to a nonzero c₂ (the smallest |c₂| is 3.27e-2, far above the 40-digit floor), with a residual falling two decades per decade in every one of the 446 ratios: no u³ term in any
of them. Both halves are exact since F161: c₂ ≠ 0 is ΔM₃ ≠ 0, and the absent u³ term is ΔX₂ = ΔX₄.

**Why there is no u³ term (R4c for eleven of them, F161 for all of them).** K = diag((−1)^site), the chiral K of
`ChiralKClaim` and of PROOF_K_PARTNERSHIP's Lemma 1 (KHK = −H for bipartite nearest-neighbour
hopping, non-uniform J included), flips every chain bond and, at odd N (even n), leaves the wrap
bond alone, so K·H(u)·K = −H(−u) and, with the levels labelled by descending order,
E_k(u) = −E_{n−k}(−u) (gate R4c reads at u = ±1e-2 and ±1e-3, the negative side labelled by descending
energy as well, which this identity itself justifies). For a
Θ-mirror pair σ = n − τ the gap is therefore S_τ(u) − S_σ(u) = S_τ(u) + S_τ(−u), EVEN in u to all
orders: every odd coefficient vanishes, and c₂(gap) = 2·c₂(S_τ), twice one branch's own shift. That is
the shape of the ζ² anti-protection law of
[PROOF_ZETA2_ANTI_PROTECTION](../docs/proofs/PROOF_ZETA2_ANTI_PROTECTION.md) with u for ζ: the same
knob-flip identity, the same conclusion, a bond instead of site fields, and the linear half K of that
law's antiunitary Θ = T·K doing the work here. Two fences of that proof, before it is quoted further.
Its first, *"Θ-mirror pairs ONLY"* (F129's entry shortens it to *"Mirror pairs ONLY"*), this page
shares. Its second reads *"Open chain, NN hopping … A ring
or beyond-NN hop reintroduces parity strings and the statements need revisiting"*; the road is a ring
at every u > 0, and the fence does not bite here because that proof's exact statements are about a
3-magnon COMPOUND built through the Jordan-Wigner map, while this page never forms one: it sums three
eigenvalues of the N×N single-excitation block, on which K acts without strings. And the parity is
complementary: the ζ² law's flown case is n = 9, odd, where a Θ-mirror pair separates on this road at
first order; the road's protected Θ-mirror pairs are at even n. Gate R4c reads the odd part of the gap
for all 223: for the 11 Θ-mirror pairs it sits at the 40-digit floor (worst 2.3e-41), a theorem measured;
for the other 212, none of which is a Θ-mirror pair, at both parities of n, the odd part starts at u⁵
(decade ratio within 15% of 1e5 for every one). So c₃ = 0 for every pair that stands at first order, and since 2026-09-02 that is a theorem for
all 223, unconditional at odd n and resting at even n on the ROT3 shape of R3c, which is read at
n ≤ 30 and not derived. The exponent 5 is settled too, though not as a theorem about all of them: it is an exact
decision per pair, ΔX₆ ≠ 0, which holds for all 212 here and fails for exactly the 11 Θ-mirror
pairs, whose gap is even in u throughout.
[F161](../docs/ANALYTICAL_FORMULAS.md) carries the level's motion to fifth order and finds each
of the five computed coefficients of the gap on a MULTIPLIER LADDER: the odd orders read the comb under
k ↦ k(n+2j), an object it writes X_2j = −M_{n+2j}, and that map is a Galois automorphism of
ℚ(ζ_2n) exactly when gcd(n+2j, 2n) = 1. At odd n the criterion is gcd(j, n) = 1, which j = 1
and j = 2 always meet, so ΔX₂ = ΔX₄ = 0 and c₃ = 0 for EVERY collision pair there, standing
or separating, the j = 1 rung being R3b above. At even n no rung is an automorphism and the
twelve are carried instead by the ROT3 shape of R3c. In both cases j = 3 is the first surviving rung, for two different reasons: at odd n because every
firing modulus there has 3|n, F129 firing only at 3|n or 10|n and 10|n forcing n even, so the
multiplier stops being an automorphism; at even n because 3|j is what collapses the ROT3 coset. X₆
first enters at FIFTH order, which is the u⁵ measured here.
The coefficient is c₅ = −(3/4)(2n−3)(n−3)·(4(n−6)(n−2)/(15n⁵))·ΔX₆. See
[PROOF_COLLISION_GAP_ODD_ORDERS](../docs/proofs/PROOF_COLLISION_GAP_ODD_ORDERS.md).

**At u = 1/2, every former collision is open (R5).** The smallest gap between ANY two distinct
clean-triple sums, over all pairs, and the range of the former collision pairs' gaps, at 40 digits:

| n | clean triples | smallest gap between distinct sums | former collision pairs' gaps |
|---|---|---|---|
| 9 | 32 | 3.27e-2 | 0.614 |
| 12 | 120 | 1.70e-3 | 0.062 to 0.415 |
| 15 | 280 | 2.03e-3 | 0.042 to 0.424 |
| 18 | 560 | 5.52e-5 | 1.23e-3 to 0.377 |
| 20 | 816 | 9.61e-5 | 1.78e-4 to 0.295 |
| 21 | 960 | 6.34e-5 | 0.020 to 0.321 |
| 24 | 1540 | 1.73e-6 | 4.14e-4 to 0.290 |
| 27 | 2288 | 1.49e-5 | 0.011 to 0.256 |
| 30 | 3276 | 4.86e-7 | 1.44e-6 to 0.260 |

(The n = 9 smallest gap, 0.0327233, and the smallest |c₂| above, 0.03269, are two different
quantities that happen to print alike.) This is a READING and is fenced as one: a gap of 4.9e-7 at 40
digits is not zero, but nothing here proves that no exact coincidence exists among the road's
algebraic levels at u = 1/2. The exact certificate this page holds is first order: 2335 pairs separate
by a theorem, and 11 pairs are even in u by a theorem. What it held at second order was a decade law on all 223, and F161 turned that into ΔM₃ ≠ 0 in
ℤ[ζ_2n], exactly. What it holds at finite u is forty digits.

## What this says for the arc

The arc predicted that *the comb family dissolves together*, and for the one law of the five that can
be asked the question, it does: no F129 collision at n ≤ 30 remains at any order read. But the dissolution is
not one event. It is staged, and the staging is itself a comb law: at odd n, first order is decided by
the parity count of the labels (forced, by Galois), and the speed is an integer multiple of 4/n; at even n, first order is decided by whether the collision persists on the doubled angle, which
F161 then traces to the ROT3 shape of the doubled labels, so what is met there is the shape and
not the persistence; and the Θ-mirror pairs among those are even in u by the chiral K (forced). Two senses of *forced*
run through that sentence and they should be kept apart: forced by a symmetry, hence robust (R2's
balanced pairs at even N, R4c's Θ-mirror pairs), and forced to a discrete value (R3b's set 4ℤ/n),
the memory law's sense; here the Galois structure supplies both. The four other laws are not
refuted and not confirmed by this page; they are shown not to be questions a detuned level list can
answer, each for its own reason in the table above (their conclusions do not transport; F89's
predicates and F146's hypothesis do), and F65 turned out to be the road's velocity rather than a law
to test on it. Method: the condition was solved (D in ℤ[ζ_2n]) and not sampled, the arc's first trap;
the balanced pairs are its negative control; the whole page is the F157 Δ-locus move on the u axis.

## Fences

- **The book.** The XY adjacency book of PROOF_RING_GAP_DOMINANCE; n = N + 1 is F129's modulus and
  the chain comb is 2cos(kπ/n) in units of J. The ring end u = 1 is not F129's comb (its levels come
  in the pairs F122 names as *k ↔ −k*, this page's m ↔ N−m) and is not examined here.
- **What is exact and what is read.** D and its vanishing: exact (cyclotomic vectors, the F129
  script's layer). The odd-n law: a theorem at every odd n (the Galois automorphism k ↦ k(n+2)),
  checked exactly on all 627 pairs to n = 30. The Θ-mirror pairs' evenness in u: a   theorem, measured at the floor. The even-n pairs' zero-sum and doubled-sum facts: forced for the 11
  Θ-mirror pairs (bipartiteness at u = 0; the k ↦ n − k invariance of the signed weight), read exactly
  for the other 12, whose ROT3 shape F161 turns into a theorem at every rung not divisible by
  3; the shape itself is a census observation at n ≤ 30 and no more. Second order: no longer a reading at all. F161 writes c₂ = −((n−3)/n²)·ΔM₃ and settles both halves
in ℤ[ζ_2n], the 223 that do not vanish and the 60 that do, and it identifies the 60 as one F129
family, C, at its own committed closed form 2(n−10); this page's R4 decade classes are reproduced
by a route with no float in it. u = 1/2: a 40-digit reading, not a certificate and not a decade law, since R5 reads a gap table
and no power.
The
  labelling at finite u: a convention resting on a proven simplicity (0 ≤ u < 1).
- **The orders past the first come from F161**, not from this page. Theorem G is a statement at
  the chain end, and "separates" means as
  the road leaves the comb; what a finite u does is the R5 reading. Nothing here says how the ring
  end's own coincidences behave as u leaves 1; that is F160's split, already solved, and a different comb.
- **Not the blind seat.** The arc's own fence. This page's bond is the ring's wrap bond under a level
  law, not a seat cut under a bond deformation.
- **The modulus.** The census runs to n = 30. The odd-n first-order law is a theorem at every odd n;
  everything counted, read or listed (223, 23, 11 and 11, 212, 173, the tables) is n ≤ 30 and says
  nothing about a standing pair at a larger modulus.
- **No hardware claim.** The flown F129 pair is cited as what it is, a sighting of the collision on
  the comb; nothing here proposes a flight.

## Anchors

- Gate: [`simulations/comb_road_f129.py`](../simulations/comb_road_f129.py) (R1 the velocity two ways;
  R2 the balanced-pair negative control; R3 the census; R3b the odd-n theorem on every pair; R3c the
  even-n pairs read exactly, zero-sum and mirror flags; R4 the derivative's decade law on every
  separating pair; R4b the second order on every standing pair; R4c the parity of the gap in u on
  every standing pair; R5 the finite-u reading, not a gate), imports the exact layer of
  [`f129_level_collision_law.py`](../simulations/f129_level_collision_law.py); run committed at
  [`comb_road_f129_run.txt`](../simulations/results/comb_road/comb_road_f129_run.txt). Cross-check of the
  census column: [`f129_family_inventory.py`](../simulations/f129_family_inventory.py) `--fast`.
- Proof of the orders past the first: [PROOF_COLLISION_GAP_ODD_ORDERS](../docs/proofs/PROOF_COLLISION_GAP_ODD_ORDERS.md),
  gate [`simulations/collision_gap_odd_orders.py`](../simulations/collision_gap_odd_orders.py),
  which reproduces this page's census from the same committed layer before explaining it.
- Proof: [PROOF_CRACKED_RING_EXACT_CURVE](../docs/proofs/PROOF_CRACKED_RING_EXACT_CURVE.md) Theorem G
  (the velocity from the polynomial and from the eigenvectors) and Corollary B (the simplicity the
  labelling rests on); [PROOF_ZETA2_ANTI_PROTECTION](../docs/proofs/PROOF_ZETA2_ANTI_PROTECTION.md)
  (the Θ-mirror pairs' evenness, on the ζ axis).
- Sibling: [The Cracked Bell](THE_CRACKED_BELL.md), where the road was found and where its
  own §E is the crack's spectral reading; it lists this page among its results, and this
  line is the edge back.
- Registry: [F161](../docs/ANALYTICAL_FORMULAS.md) (the orders past the first, this page's
  census explained), [F160](../docs/ANALYTICAL_FORMULAS.md) (the road), [F129](../docs/ANALYTICAL_FORMULAS.md)
  (the law under test; its entry carries a pointer to this page), [F65](../docs/ANALYTICAL_FORMULAS.md)
  (the velocity's magnitude, and the doubled angle), F75 (the sign), F123 (the Re-side sibling in the
  move), F130 (the nearest bond detuning of a collision object).
- Hardware: Confirmation 24, `f129_standing_fringe_kingston_july2026`
  ([IBM_F129_RAMSEY_FRINGE](IBM_F129_RAMSEY_FRINGE.md)), the n = 9 pair on the comb.
- Typed: `CrackedRingExactCurveClaim` (`compute/RCPsiSquared.Core/Symmetry/`, node *the velocity at the
  chain end is F65's comb, signed*); the road's levels live in `compute/MirrorWorld/Crack.cs`.
- Arc: `the_forced_and_the_met` NextStep (1), which this page answers for F129 and leaves open for
  nothing else it named, the other four being shown not to transport.
- Record: [`docs/CAUGHT_ERRORS.md`](../docs/CAUGHT_ERRORS.md), both entries of 2026-09-02.
