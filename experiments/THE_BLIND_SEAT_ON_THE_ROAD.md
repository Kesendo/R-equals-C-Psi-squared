# The Blind Seat on the Road: which half of F157's count survives a boundary knob

**Date:** 2026-09-02. **Authors:** Thomas Wicht, Claude (Opus 5).
**Arc:** [`the_forced_and_the_met`](../compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs).
Its item (1), *"RUN THE COMB TEST"*, the arc now records as done — but with a qualifier this page
keeps: *"item (1) is DONE for the one law that could be asked and is not a verdict on the other
four."* It was closed on 2026-09-02 by [The Comb on the Road](THE_COMB_ON_THE_ROAD.md) and F161; the
live items are (2) and (3). Inside
item (3) sits the trap this page obeys: *"THE_SEAT_THAT_CUTS's own open item asks for a DETUNED BOND
rather than a Delta, so do not report one as the other."* That is a fence, not a commission; this
page runs the detuned bond and keeps the two apart.
**Laws:** [F157](../docs/ANALYTICAL_FORMULAS.md), the seat-cut blindness laws and their Δ locus, and
[F161](../docs/ANALYTICAL_FORMULAS.md), which reads each computed order of a collision gap as one
comb under an integer multiplier and decides by a gcd which orders vanish
([PROOF_COLLISION_GAP_ODD_ORDERS](../docs/proofs/PROOF_COLLISION_GAP_ODD_ORDERS.md) — the file name
says "odd orders", but §(a) below leans on the M ladder, which carries the even ones).
**Instrument:** [F160](../docs/ANALYTICAL_FORMULAS.md), the cracked ring's exact curve
([PROOF_CRACKED_RING_EXACT_CURVE](../docs/proofs/PROOF_CRACKED_RING_EXACT_CURVE.md)), whose road
u = J′/J on the wrap bond leaves the open chain at u = 0 and reaches the ring at u = 1. J itself
never enters: every H below is in units of J.
**Gate:** `python simulations/blind_seat_on_the_road.py`
([source](../simulations/blind_seat_on_the_road.py)), 25 checks under 20 labels in three blocks
(A4 fires three times, C1 four), about 21 seconds measured quiet under **sympy 1.14.0**, which is
the version the numbers below were produced with and matters because the third trap is a
sympy-behaviour trap with two instances. sympy is not in the dependency line in `CLAUDE.md`. Run at
[`blind_seat_on_the_road_run.txt`](../simulations/results/blind_seat_on_the_road/blind_seat_on_the_road_run.txt).
**Registry:** unregistered. §(c) was an exactly solved relation over a finite range of N with its
mechanism identified but not derived, and
[PROOF_BLIND_SEAT_TWO_AXES](../docs/proofs/PROOF_BLIND_SEAT_TWO_AXES.md) derives it: the underlying
decomposition holds at every N, §(c)'s relation is a theorem at every odd N, and at even N it
becomes a criterion which turns this page's gate C8 into a consequence and names its 22 breaking
seats.

**Notation, once.** N is the chain's site count; **n = N + 1** is the *comb modulus*, the letter
F129, F161 and this page share. A **comb** is the arithmetic progression of angles a chain's
standing waves are forced onto — kπ/n for the open N-site chain — and the repo also uses the word
for the spectrum {2cos(kπ/n)} it produces and for a law whose content is a coincidence among those
cosines. It has no glossary headword. A **seat** j ∈ {0, …, N−1} is where the watching sits;
**blind(j)** is F157's *count*, N minus the dimension of the Krylov space of e_j, equivalently
deg gcd(χ(H), χ(H_j)) where H_j is H with row and column j struck; a seat is **sighted** when
blind(j) = 0. The **aliasing degree** of a multiplier m on ℤ/mod is mod divided by the size of the
image of multiplication by m (written μ_m); it equals gcd(m, mod) and is never computed that way
below. **N_node = |N−1−2j|** is F157's node modulus for seat j. A **locus** is the set of real
parameter values at which a given seat is blind; it is the real root set of Res_x(χ(H), χ(H_j)) in
that parameter, which the page writes **t** generically, **u** for the crack and **Δ** for the
anisotropy. The two sites {0, N−1} are the **end pair**, not "the corner" (F140 owns "corner block"
for a different object). This page writes H = A + Δ·(E₀₀ + E_{N−1,N−1}); F157's committed convention
writes hop 2, 2Δ per end and a Δ(N−5)·I shift, which differ by a common scaling and a common shift,
neither of which moves a coincidence between χ(H) and χ(H_j) — gate C1 checks this rather than
asserting it.

## What this is about

Put a light on one seat of a chain of spins and ask what it cannot see. The chain's standing waves
carry everything that moves, and a wave that is exactly zero at the lit seat is never touched by the
light. Counting those waves is the seat's blindness, and on 2026-08-24 we learned which seats are
blind: the answer was a divisor.

This page asks which seats keep that blindness when the chain stops being perfect.

A wave can vanish at a seat for two reasons. It can be **forced**: the chain looks the same from
both ends, so every wave antisymmetric about the middle must vanish at the middle seat, whatever the
chain is made of. Or it can be **accidental**: on the even chain the waves are sines, and for
certain lengths and seats a sine lands on zero by arithmetic alone. On the perfect chain the divisor
counts both without distinguishing them.

So we spoiled the chain, two ways — once by joining the two ends with a weak bond, once by putting
an energy offset on those same two end sites — and watched. The middle seat of an odd chain stayed
blind at every setting of either knob, and by the same amount, both provably. Every other blind
seat lost its blindness at once and got it back only at isolated, exceptional settings. Then the
part we were not looking for: over the lengths we could reach, the exceptional settings of the two
quite different knobs are the same ones, apart from the two settings at which the chain closes
into a ring. That holds at every odd length we tried and at two thirds of the even ones. We
already had a closed formula for one knob's settings. Where it holds, it is the formula for the
other's as well.

## What the repo already held, store by store

Swept 2026-09-02 by three agents, adversarially by a fourth, and corrected by three empty review
rounds. **Three** entries below were called "nothing" by an earlier sweep of this page and are not
nothing; all three are recorded rather than quietly filled, because `CAUGHT_ERRORS.md` had logged
that exact shape on a neighbouring page earlier the same day, and this page then repeated it three
times.

- **[`docs/ANALYTICAL_FORMULAS.md`](../docs/ANALYTICAL_FORMULAS.md), F157 — the load-bearing
  prior.** It owns the whole Δ axis in closed form: the node-count reduction
  gcd(j+1, N_node) = gcd(j+1, N+1); the generating formula
  Δ_k = sin((j+1)kπ/N_node)/sin(jkπ/N_node) over those **k ∈ 1..N_node−1** whose denominator does
  not vanish — F157 calls that range *"not decoration"* — with the locus empty exactly when
  N_node | j; its resultant packaging
  P_j(Δ) = Res_x(U_{N_node−1}(x), Δ·U_{j−1}(x) − U_j(x)); a four-row table of worked seats that the
  entry itself fences as *"an ILLUSTRATION and not the content"*; the
  statement that *the multiplicity* of a root of P_j is the blind count at that Δ; and the live
  witness `dotnet run --project compute/RCPsiSquared.Cli -- inspect --root blindlocus`
  (`SeatBlindnessDeltaLocusWitness`). Gates C1 and C2 pin against the table and against the
  generating formula. It also states F157's fence-free identity as *"the identity unconditional for
  real symmetric H"*, by Cramer's rule — every locus here is a resultant of those two
  characteristic polynomials and inherits it, so nothing below needs the unreduced-Jacobi
  hypothesis that `PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA` places on a *different* identity,
  deg gcd(χ_L, χ_R). That proof draws the same line itself: *"the two-halves phrasing goes wrong
  there while the fence-free Cramer identity of The Seat That Cuts §7 does not"*.
- **F157's standing fence, on §(c).** *"A Δ is NOT the detuned bond that
  [The Seat That Cuts](THE_SEAT_THAT_CUTS.md) leaves open; do not report one as the other."* §(c)
  asserts that two *loci* coincide as sets. The fence still holds and is worth stating exactly,
  because [PROOF_BLIND_SEAT_TWO_AXES](../docs/proofs/PROOF_BLIND_SEAT_TWO_AXES.md) does identify
  something: not the two operators, which sit in different positions on the full space, one a bond
  and one a diagonal, but their RESTRICTIONS to each reflection sector, and there at opposite
  signs of the knob. A Δ is still not a detuned bond.
- **[`docs/proofs/`](../docs/proofs/).** `PROOF_COLLISION_GAP_ODD_ORDERS` §(e) supplies the
  M-ladder criterion, *"For the M ladder m = 2j + 1 is odd already and the criterion is
  gcd(2j + 1, n) = 1"* — a remark in that section, not Theorem D, whose Statement is about the X
  ladder. §(h) supplies the fence §(a) obeys: *"The criterion gcd(m, 2n) = 1 is sufficient and too
  crude, and the sharp one is local to the pieces."* §(b), Lemma B, *"the wrap bond is rank one
  inside each reflection sector"*, is why a crack and a diagonal can share a locus at all, and is
  quoted in §(c) rather than reproved. When this page was written nothing in `docs/proofs/` carried
  §(c), and `PROOF_CRACKED_RING_EXACT_CURVE` said so of itself: *"Not the blind seat:
  THE_SEAT_THAT_CUTS's open item asks for a detuned bond under a seat cut, a different object."*
  That is the object this page built, and
  [PROOF_BLIND_SEAT_TWO_AXES](../docs/proofs/PROOF_BLIND_SEAT_TWO_AXES.md) now carries §(c);
  the fence is still exact, since what that proof identifies is sector restrictions and not the
  two perturbations.
- **[`experiments/`](../experiments/).** `THE_SEAT_THAT_CUTS` §4 states the node condition
  m(j+1) ≡ 0 (mod n) and gives one derivation with three verification routes; its §7 owns the
  Cramer identity. `THE_COMB_ON_THE_ROAD` names F157 three times, twice as a *method* precedent and
  once as an exclusion, never as the same arithmetic. `XY_FROZEN_BAND` is quoted in the fences.
- **[`fw.Confirmations`](../simulations/framework/confirmations.py) — NOT nothing, and an earlier
  sweep of this page said it was.** Confirmation 24, `f129_standing_fringe_kingston_july2026`, is
  the F129 census's n = 9 collision pair (1,5,7) ~ (2,4,8), flown 2026-07-15 on ibm_kingston — the
  same comb §(a) makes an implication about. What it measured is a Ramsey fringe phase slope
  against Trotter step count at fixed θ, with ζ entering only as a pre-registered systematic
  budget term; it is on neither the u axis nor the Δ axis, and nothing on this page is confirmed by
  it, the blindness count never having been flown. But the store is not empty, and
  `CAUGHT_ERRORS.md` recorded this identical miss on `THE_COMB_ON_THE_ROAD` earlier the same day.
- **[`reflections/`](../reflections/) — NOT nothing either.**
  `ON_LEAVING_THE_CIRCLE.md`, an eleven-visit reflection on whether F89 and PTF are connected,
  carries in one bolded sentence — *"The seat that cannot be mirrored is how the circle touches the
  line"* — the reading of the unmirrorable middle seat of an odd chain: the object §(b)
  calls forced. It speaks of "the discrete comb" and "a comb of cosines" further down. It carries no
  count and no locus, so it changes nothing computed here, but it is the reading that already
  existed for §(b)'s survivor.
- **The open-arcs ledger.** `the_forced_and_the_met` sorts our laws by track record into ones whose
  arithmetic is *"a COUNTING argument … rather than a coincidence among cosines"* and comb laws, and
  predicts that *"the comb family dissolves together"*. §(b) revises where F157 sits and leaves the
  prediction untouched. The slug `blind_site_divisor_law`, which an earlier draft of this header
  named as a second arc, **is not an arc**: it is a private memory slug, the trap the memory index
  already records for `w1_dispersion_anchor`.
- **[`compute/MirrorWorld/`](../compute/MirrorWorld/).** `Crack.cs` is what sent this page down the
  road: a class-header comment recording that the world held the two combs as a switch *"three times
  over (Cyclotomy.Own, Divisor's clock modulus, BlindSeat's two gcd laws) and as a road zero
  times"*. `Cyclotomy.cs` owns *"the order of a turn fraction"*, n/gcd(j, n); `BlindSeat.cs` calls
  `Cyclotomy.Gcd` but not that order.
- **[`docs/GLOSSARY.md`](../docs/GLOSSARY.md) — NOT nothing, and this is the third store an
  earlier sweep of this page wrongly called empty.** It carries a section headed *"The blind seat
  (F157, August 2026)"* defining seat, the blind subspace, blind(j) = N − dim Krylov(e_j), F157's
  Breaks-for fence and both uniform gcd laws — that is the notation block above, already written.
  What is genuinely absent is the **comb** vocabulary: zero occurrences of "comb", no headword for
  a comb multiplier or a node mode, and "Galois" once, in an F89 row about a different object. That
  absence is why the notation block defines "comb" and should not have been extended to the rest.
- **[`hypotheses/`](../hypotheses/) and [`recovered/`](../recovered/): nothing on the loci.**
  `hypotheses/DIABOLIC_BY_INTEGRABILITY.md` turns a Δ anisotropy on this same chain family and is
  cited by the arc's own Origin as method prior art, but it carries no blindness count and no
  locus.
- **[`docs/CAUGHT_ERRORS.md`](../docs/CAUGHT_ERRORS.md).** Its 2026-08-30 fourth entry records that
  the Δ-locus test bench was stated wrongly and replaced by the closed form the same day. Its
  2026-09-01 fourth entry names the failure this page's first sweep then repeated: *"the sweep
  searched the OBJECT (a cracked ring) and not the METHOD (an integer polynomial identity, a
  Chebyshev polynomial)"*, naming F157's locus polynomial as what was missed.

## (a) On the uniform chain the count is an aliasing degree

F157's node condition is m(j+1) ≡ 0 (mod n), from ψ_m(j) ∝ sin(πm(j+1)/n). Read as a map rather than
a congruence its solutions are the kernel of μ_{j+1} on ℤ/n, so

    blind(j) + 1  =  aliasing degree of μ_{j+1} on ℤ/n,

gate A1 at every seat of N = 3..30, the left side from F157's *definition* (N minus the seat's
Krylov rank, exact rational elimination) and the right by enumerating an image; no gcd is called on
either side. 158 seats alias, so the identity is not 1 = 1, and reading the aliasing on a wrong
modulus disagrees 154, 313 and 253 times (2n, N, 2N) on that same range.

**This is the whole count, not one part of it.** At the centre seat of an odd chain the formula
returns gcd(n/2, n) − 1 = (N−1)/2, which is exactly the forced count §(b) attributes to the
reflection. The arithmetic finds the reflection-odd modes too; it does not know they are forced.
Forced and accidental is therefore not a decomposition of this number but a statement about which
seats keep their count when the chain is spoiled, and §(b) is where it is measured.

This much is a relabelling: the registry writes the same count out as *"N_node | (j+1)k ⟺
(N_node/h) | k with h = gcd(j+1, N_node)"*. What it buys is a consequence. μ_m is a Galois
automorphism of ℚ(ζ_2n) iff gcd(m, 2n) = 1, which for **odd** m is the sighting condition
gcd(m, n) = 1 — gate A2, over N = 3..30, with 39 of the odd-multiplier seats not sighted. For even m
the criterion is simply silent (gcd(m, 2n) ≥ 2 always, so μ_m is never an automorphism) while 105
such seats are sighted anyway; gate A3 exhibits that, because A2's restriction is a fence and not a
convenience. Given A1, A2 and A3 are one-line facts about integers; A1 is the check with a matrix on
one side.

F161's §(e) needs that automorphism on this same n. Its M ladder carries m = 2j+1, so:

> **For a collision pair at odd n = N+1: seat 2j of the open chain sighted ⟹ F161's M-ladder
> rung j is killed.**

The hypothesis is not decoration: §(e)'s kill transports ΔM₁ = 0, which *is* the collision, so
with no pair there is nothing to transport and nothing is killed. One direction only — §(h) calls the criterion *"sufficient and too crude"*, and F161's census has 40
pairs at n = 30 with c₂ = 0 whose multiplier is no automorphism, forced by the sharp local criterion
instead. Three things the sentence does not carry on its face: F161's j is a **rung index** and this
page's j is a **seat**, one letter for two objects; the seat 2j must exist, so j ≤ (N−1)/2, and past
that the multiplier needs reducing mod 2n first; and nothing here computes a collision pair, a ΔM or
anything in ℤ[ζ_2n] — the F161 side is quoted from its proof, and a gate that earned more would
import that census.

## (b) Which seats keep the count

A mode vanishes at a seat for one of two reasons.

**Forced.** If a reflection of H fixes the seat, every reflection-odd mode vanishes there. This is
representation theory: no length, no arithmetic, no coupling enters, so nothing detunes it away.

**Accidental.** The sine lands on a multiple of π by arithmetic. This is §(a)'s reading of the same
number, and it is a comb coincidence.

The knob separates them, and it separates *seats*. Gate B1: at odd N = 5..13 the centre seat's
resultant vanishes **identically** on both axes — a statement about the polynomial, not about
sampled points, though the reflection argument already forces it and B1 is therefore a consistency
check on the construction rather than evidence. Gate B2 adds what B1 does not say, that the
centre's *value* stays (N−1)/2 at six rational couplings on both axes, over odd N = 5..11; the
reflection argument already gives ≥ (N−1)/2 as a theorem, so what B2 tests is the ≤ direction,
which
[PROOF_BLIND_SEAT_TWO_AXES](../docs/proofs/PROOF_BLIND_SEAT_TWO_AXES.md) §(g) then gives at every
coupling rather than at six: in the folded even block the centre is an END of a Jacobi chain, and
an eigenvector there cannot vanish. Gate B4 is the anti-vacuity partner: all 16
accidentally blind seats of N = 5..11 have a finite locus **on the u axis**, so they do move; the Δ
leg falls out of C2 instead. Gate B5 records that at
the centre the whole count is the forced one, which given A1 is a line of integer arithmetic rather
than independent evidence.

The ring end is the same law again rather than a second phenomenon: a ring has a reflection through
each site, so every seat is forced there, and gate B3 finds the constant at every seat at **both**
parities over N = 5..13, (N−1)/2 at odd N and (N−2)/2 at even.

**Where this leaves the ledger's cut.** The arc places F157 among the laws that have survived
detuning, and that is right; what survives is the forced count plus a definition that outlives the
comb. The tempting next step — "the real cut is comb laws with a counting continuation versus
without" — is a two-case promotion from one law, and §(c) shows it is also the wrong cut. The arc's
prediction is untouched by this page.

## (c) Where the accidental blindness goes

Off the comb the accidental blindness survives on a finite set: the real roots of
Res_x(χ(H), χ(H_j)) in the knob. On the Δ axis that set is **already owned**. Gate C1 checks this
page's construction against all four rows of F157's committed table — N = 9 seat 1 → Δ⁵ − 4Δ³ + 3Δ;
N = 9 seat 2 → 2Δ² − 1; N = 11 seat 1 → Δ⁷ − 6Δ⁵ + 10Δ³ − 4Δ; N = 11 seat 2 → 3Δ⁴ − 4Δ² — as **real
root sets**, which is all it can check: the resultant this page forms is a rational multiple of
F157's P_j *squared*, so the multiplicity that F157 says carries the blind count at each Δ is
present doubled, and `_squarefree` discards it before any comparison. No gate here reads it.
Why it is doubled is [PROOF_BLIND_SEAT_TWO_AXES](../docs/proofs/PROOF_BLIND_SEAT_TWO_AXES.md)
§(g): the squaring happens inside each reflection sector, because striking a seat disconnects
that sector's tridiagonal block into two halves and a blind eigenvalue is a common root of both.
The multiple THIS page discards is the one §(g)'s gate K1e reads, and K1e asserts only that it is a t-free nonzero constant, not which. What §(i) and
F162 add is one level down: each sector's halves-resultant is an integer constant times the
product of (t − Δ_k) over the non-pole node indices of its parity, and the two multiply to
Res(α_p, α_{N−1−p}) times (−1)^binom(p+1,2) in the fold coordinate p, that object being F157's own
generator up to a second sign §(i) also gives in closed form. So the constants inside the
factorisation are numbers now, whether or not this page's own normalisation still throws one away.
So the count is recoverable from the same resultant by halving, exactly, off the ring ends and
away from the forced centre seat. The typed home of that level is
[`BlindSeatSectorFactorisationClaim`](../compute/RCPsiSquared.Core/Symmetry/BlindSeatSectorFactorisationClaim.cs)
(`inspect --claim BlindSeatSectorFactorisationClaim`), which recomputes it from (N, j) alone.
At the ring ends the two axes part at every seat but that one,
which is measured over N = 5..11 rather than proved. Gate C2 then checks F157's *generating* formula Δ_k at all 52
non-centre seats of N = 5..11,
of which 26 have a nonempty value set and 26 are empty (the 14 end seats by F157's degenerate case,
the rest by N_node | j); it certifies that the locus is the set of *real conjugates* of the Δ_k,
marginally weaker than the set claim. Since its `P is None` branch never fires, C2 also establishes
what B4 leaves open on the Δ axis: no non-centre seat of N = 5..11 is identically blind there.

That closed form is itself a comb: at seat 1, Δ_k = sin(2kπ/N_node)/sin(kπ/N_node) = 2cos(kπ/N_node).
F157's fence-free continuation does **not** step outside comb country; it moves the comb from the
spectrum onto the parameter axis and changes the modulus from n to the seat-dependent N_node. That
is a sharper question than the one §(b) declines to promote: not *does the law have a counting
continuation*, but *which modulus is its comb on, and does the detuning move that modulus or destroy
it*.

This page's own result is that the crack's axis carries the same set:

> **At odd N = 5, 7, …, 17, u-locus(j) = Δ-locus(j) ∪ {+1, −1}, at every seat.**

Gate C3: all 77 seats. The relation has content at the 48 whose Δ-locus is nonempty; 22 more read
u-locus = the ring ends alone; and the 7 forced centres are identically blind on both axes and have
no locus polynomial at all. The ±1
are the ring end, where the spectrum degenerates; they are not exclusive to the u axis, since the
Δ-locus contains ±1 at N = 9 seats 1 and 7. Because F157 generates the Δ side in closed form, this
hands the u axis one it did not have.

The measured range is this page's; the statement is not bounded by it.
[PROOF_BLIND_SEAT_TWO_AXES](../docs/proofs/PROOF_BLIND_SEAT_TWO_AXES.md) proves it at every odd N,
and reads the ±1 as the third term of a decomposition rather than as an appendix: where the two
reflection sectors share an eigenvalue the kernel is two-dimensional and therefore contains a
vector vanishing at any single site, so a degeneracy blinds every seat at once. That is why the
ring ends appear at every seat, and it is a third entry beside the *forced* and *accidental* of
§(b) above.

**What is not true: that the end pair is what does it.** Two earlier versions of this page claimed
so, and got the question wrong in opposite directions. The first asked whether an interior
perturbation satisfies the *crack's* relation, locus = Δ-locus ∪ {±1} — but ±1 is the u axis's ring
end and a diagonal has none, so every control failed before its locus was compared. The second
scored the controls on **equality** with the Δ-locus while scoring the crack on the relation *with*
the ring ends: two predicates, two denominators.

Asked the same way of everything, on the 48 seats with a nonempty Δ-locus, the predicate that means
*carries the whole Δ-locus* is containment:

| perturbation | locus = Δ-locus | Δ-locus ⊆ locus |
|---|---|---|
| the crack (wrap bond) | 8 of 48 | **48 of 48** |
| an interior bond, asymmetric | 0 of 48 | 36 of 48 |
| a symmetric interior bond pair | 0 of 48 | 36 of 48 |
| a symmetric interior diagonal pair (1, N−2) | 12 of 48 | 44 of 48 |
| a symmetric interior diagonal pair (2, N−3) | 16 of 46 | 22 of 46 |

The table is a **read**, not a gate. On **equality** the crack loses, and must: it always carries
the ring ends as well. On containment it is alone at 48 of 48 — which is C3 restated rather than a
second measurement, since C3's equality already implies it — and gate C7, the one assertion here
that carries new information, is that none of the four interior perturbations tried reaches that.
The best of them is at 44, so the margin is four seats, not a gulf. Which interior
perturbations carry the locus, and at which seats, is an open question this page did not have before
the retraction.

**Why the two axes can meet.** By Lemma B the wrap bond is rank one inside each reflection sector,
with a sector-dependent sign η; the anisotropy is a diagonal on the same two sites and is likewise
rank one per sector, ⟨ψ_k|D|ψ_l⟩ = a_k a_l(1 + η_k η_l) against the crack's 2η·a aᵀ, where a_k is
mode k's amplitude on the end pair and η_k = ±1 its reflection parity. Two end-pair objects
differing only in that factor.

One step further and the two stop merely resembling each other. On the R-even subspace
ψ_{N−1} = ψ_0, so the wrap bond acts there exactly as the end-pair diagonal, and on the R-odd
subspace exactly as its negative: the cracked ring IS the anisotropic open chain, sector by sector,
at +u and at −u. `PROOF_CRACKED_RING_EXACT_CURVE`'s Corollary B already folds the crack that way;
what was missing was the line that folds the anisotropy the same way and compares. That line, and
what follows from it at both parities, is
[PROOF_BLIND_SEAT_TWO_AXES](../docs/proofs/PROOF_BLIND_SEAT_TWO_AXES.md).

**Why odd N, and what parity actually decides.** Write Σ = diag((−1)^l). Gate C5 verifies as exact
matrix algebra that Σ·H_u·Σ = −H_u(−t) at odd N and fails at every even N, the wrap entry picking up
(−1)^{N−1}; C5b that it holds for the anisotropy at every N. Since {±1} is negation-closed too, the
right-hand side of C3 is negation-closed at every N, so a u-locus that is not cannot equal it:
**non-closure forces a break, and that half is a theorem no check can report false.** What gate C8
measures is the converse, and over even N = 6..16 every break is a non-closure — the same 22 seats,
N = 8 seats 1 and 6; N = 10 seats 1 and 8; N = 12 seats 1, 2, 3, 8, 9, 10; N = 14 seats 1, 2, 11,
12; N = 16 seats 1, 2, 3, 4, 11, 12, 13, 14, eleven reflection orbits, with no break of any other
kind.

Parity therefore does not partition the behaviour. The relation still **holds at 44 of the 66
even-N seats**. What this page's own identity buys at odd N is weaker than forcing: it makes the
u-locus negation-closed, which is *necessary* for the relation and never sufficient. The sufficient
half is a sharper reading of the same Σ, one sector down.
[PROOF_BLIND_SEAT_TWO_AXES](../docs/proofs/PROOF_BLIND_SEAT_TWO_AXES.md) reads ΣR = (−1)^{N−1}RΣ:
at odd N, Σ preserves each reflection sector, so each sector's locus is negation-closed **on its
own** and the relation follows; at even N it SWAPS them, so the odd-sector locus is the negation of
the even one and the relation holds exactly when the u-locus is negation-closed. C8's biconditional
is that statement, so the converse it measures is a consequence rather than a coincidence, and the
22 breaking seats are predicted rather than found.

Seen from the data side alone, the sector-split shape is visible at a break but not universal, and
that is what gate C9 *asserts*. The rate beside it is a read: u minus the ring ends takes exactly
one member of each ± pair of the Δ-locus at **20 of the 22 breaks**, failing at N = 12 seats 1 and
10, where the Δ-locus contains ±1 itself and the u-locus carries those two only as ring ends. The
proof's criterion carries the same clause and the same two exceptions, so the 20-of-22 is the shape
of that clause and not a defect in the reading. An earlier version of this page asserted the shape
at every break and gated nothing.

## Scope and fences

- **§(c) is a relation between two loci, not an identification of two objects.**
- **This page measures §(c) at odd N = 5..17 and at 44 of the 66 even-N seats of N = 6..16,
  breaking at the other 22.** Gate C4 is required to report breaks, so the measured scope is
  earned rather than decorative. The statement itself is not bounded by that range:
  [PROOF_BLIND_SEAT_TWO_AXES](../docs/proofs/PROOF_BLIND_SEAT_TWO_AXES.md) proves it at every odd
  N, gives the even-N criterion, and derives the 22 seats. This page's staggering identity supplies
  only the *necessary* half; the sufficient half is the same Σ read one sector down.
- **Every LOCUS here is a SET, and no gate here compares two counts.** Gates B2, B3 and B5 do
  read a count, at the centre seat and at the ring end; what this page never does is compare the
  count on one axis with the count on the other, and the resultant it forms is a rational multiple
  of F157's P_j *squared*, so the multiplicity survives doubled and `_squarefree` discards it
  before any comparison. What the counts do is in
  [PROOF_BLIND_SEAT_TWO_AXES](../docs/proofs/PROOF_BLIND_SEAT_TWO_AXES.md) §(g), and it is not
  the same answer as for the loci. Off the ring ends |u| = 1 the two axes carry one count at odd N,
  while at even N the crack pays the even sector twice and the chain pays each sector once; AT the
  ring ends the crack is blind at every seat at once, at either parity. Where the
  loci break the counts break with them, provably; whether the converse holds is measured over
  N = 6..14 and not proved.
- **The u ≥ 0 clause this page respects is the departure COUNT's, and it lives in
  `PROOF_CRACKED_RING_EXACT_CURVE` rather than in F160** (F160 has a u ≥ 0 clause of its own, on
  simplicity): *"its departure COUNT is stated for u ≥ 0, the object's range (at negative u
  the count reads differently and is not this file's)"*. Every locus here is the full real root set
  and contains negative values; no departure count is used, but the negative half of each locus is
  this page's own object and carries none of that file's readings.
- **"The end pair is special" is refuted for this locus and is not claimed here.** Independently,
  `XY_FROZEN_BAND` records a refutation in the same direction inside a positive law. Correcting an
  earlier reading of its own that blamed the diagonal, it finds an XY chain keeps its diagonal band
  *"for both a diagonal on the end sites and one on the inner sites"* — end against inner is the
  control arm of that refutation, which is why it bears on this page.
- **The centre-seat argument is the reflection-odd one, not the two-block one.** Striking the centre
  of a cracked ring leaves one path, so the two-blocks-conjugate-by-reversal argument is a u = 0
  argument. What transfers is that i ↦ 2c − i (mod N) is a symmetry of H_u at every u fixing only
  the centre.
- **§(a)'s aliasing reading is a relabelling of a committed derivation**; only the F161 implication
  is new, and it is one-directional and bounded to j ≤ (N−1)/2.
- **The XY book only.** Every H here is the single-excitation XY block, F160's book. F157's
  Heisenberg law, on modulus N, is a different comb and is not on this road.
- **No hardware claim.** Confirmation 24 stands on the same comb as §(a)'s implication but is on the
  ζ axis and confirms nothing here; no flight is proposed.

## Three traps

Kept because a next session will walk into them; the rest of this page's drafting history goes to
`docs/CAUGHT_ERRORS.md` in the change that lands this page, not here.

**The three moduli.** F157's Heisenberg law lives on modulus N, its XY law on N+1, F161 on n = N+1.
F161 matches the Heisenberg law in its *multiplier* and the XY law in its *modulus*, and neither in
both. A join asserted on two matching gcd's is asserted on the wrong one.

**A sufficient condition reads like an equivalence.** "Killed exactly when" survived three readings
of this page before §(h) was opened. The Galois criterion is one-way wherever it is invoked, and the
proof says so in bold.

**A convenience function must not judge an exact question, and the failures are silent.** Twice, in
this gate. `sp.roots` returns an irreducible cubic's roots in *casus irreducibilis* form, where
`.is_real` is `None`, not `True`; a truthiness filter dropped them, and every even-N count on an
earlier version of this page was wrong — the break list read 4 seats where it is 22, and omitted
N = 10 entirely. Then `sp.simplify` failed to reduce P(2cos(3π/7)) to 0 and reported a true root as
a counterexample. Both are now decided by polynomial arithmetic: squarefree parts, Sturm counting,
minimal-polynomial divisibility, and exact expansion in place of the remaining `simplify` zero
tests.

## What this opens

1. **Which sector each of F157's k belongs to. CLOSED 2026-09-03.** The two questions this item
   carried are both answered in [PROOF_BLIND_SEAT_TWO_AXES](../docs/proofs/PROOF_BLIND_SEAT_TWO_AXES.md): §(a) to
   §(f) prove the shape conjectured here, and §(g) settles the count, which this page had fenced
   off as unreachable. It is reachable, and the answer splits by parity: at odd N the crack pays
   the chain's count exactly, off the ring ends; at even N it pays the EVEN sector twice where the
   chain pays each sector once. That left one closed form, and §(h) of the same proof supplies it, a
   day later. F157 writes the Δ-locus as Δ_k over k ∈ 1..N_node−1 and reads the blind count as
   the number of k landing on a value; §(g) reads the same count as b_E + b_O. The mode at the
   node angle is ψ_l = sin((j−l)θ_k) and its reflection parity is (−1)^{k+1}, so at every seat with
   N_node ≥ 2, **b_E counts the ODD k on a value and b_O the even ones**, and the even-N count
   closes to 2·#{k odd : Δ_k = t} off the ring ends. The odd-N centre seat has no index at all
   and stays §(g)'s (N−1)/2 at every coupling. **The guess named here, a parity of k, was right**, and
   right as stated, per index: odd k lie in the R-even sector, even k in the R-odd one. What was
   missing was the proof, and the fact that the sign law is not ours to find. It is F71's,
   owned in six places for the chain's own comb, and what §(h) adds is that it
   survives the change of modulus to the seat's N_node, where the matrix is no longer uniform.
2. **Which interior perturbations carry the locus?** New, and generated by the retraction above.
   Read on EQUALITY, a reflection-symmetric interior diagonal pair reproduces the Δ-locus at 12 of
   48 seats at (1, N−2)
   and 16 of 46 at (2, N−3), an interior bond at none. Which seats, and why those, is unasked.
3. **Ask each comb law which modulus its comb is on.** F157 survives detuning because its comb moves
   from the spectrum to the parameter, on the seat-dependent N_node. For F89's seed resonance,
   F145/F146, F65's Niven root and F144's exception the question is whether a parameter-side comb
   exists at all. What an answer looks like is F157's Δ_k: a closed form in the parameter whose
   values are cosines of a rational angle.
4. **A candidate lemma, with a floor of two instances.** "A fixed point of an involution cannot be
   moved by anything commuting with it" underlies §(b)'s forced count and
   [`Seed.cs`](../compute/MirrorWorld/Seed.cs)'s N−1 forced seeds at odd N against 0 at even —
   where that 0 is a forced *lower bound*, as `Seed.cs` says itself, not an exclusion. A
   review proposed three further instances in Core; one of them,
   `BetaExoticPerNExclusionClaim`, does not carry the fixed-point count attributed to it (its parent
   `SeedExistenceCountingClaim` does), and the other two are unchecked. The sweep is unmade.
5. **F157's Heisenberg law has no road.** Its comb is the ring's, modulus N, and F160's road ends
   there. Whether the two gcd laws are the two ends of one *count* rather than of one switch is this
   page's question, generated by `Crack.cs`'s observation that the pair had been held as a switch
   "and as a road zero times"; that comment frames no question of its own.

## Where it came from

The blind seat is [The Seat That Cuts](THE_SEAT_THAT_CUTS.md) and
[The Blind Site](THE_BLIND_SITE.md), F157, whose Δ-locus closed form is the load-bearing prior here.
The road is [The Cracked Bell](THE_CRACKED_BELL.md) and its proof, F160. The ladder §(a) meets is
[The Comb on the Road](THE_COMB_ON_THE_ROAD.md) and F161. The reading §(b)'s survivor already had is
[On Leaving the Circle](../reflections/ON_LEAVING_THE_CIRCLE.md). The sentence that sent this page
down the road is a class-header comment in
[`compute/MirrorWorld/Crack.cs`](../compute/MirrorWorld/Crack.cs).
