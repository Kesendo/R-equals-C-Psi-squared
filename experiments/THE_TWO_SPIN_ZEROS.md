# The two spin zeros

**2026-08-16.** [THE_ENDPOINTS_ARE_A_DENSITY_LAW](THE_ENDPOINTS_ARE_A_DENSITY_LAW.md)
measured, at Δ = 1 and N = 6, that mixedness is necessary but not sufficient for
breaking the reflection law C_l = 0: on (1,3) two parity-mixed eigenspaces break
and two hold at the floor. The arc `compressed_density_laws` carried the question
as the first half of its NextStep (3): what separates them? The answer is a
rule built from two Clebsch–Gordan zeros. Its hold direction is a theorem, its
break direction is genericity, and the rule is verified 60/60 and
element-exactly on ALL mixed eigenspaces of all fifteen N = 6 blocks (48
breaks, 12 holds), and it comes with a second finding the question never asked
for: every breaking residual measures as the square root of a rational
number, and four of the forty-eight breaks sit BELOW the 0.05 threshold
the earlier census classified with, so the census here classifies by the
eigensolver floor instead. (The name says "spin zeros" and not "two zeros"
deliberately: "quotient of two zeros" already means something else in this
repo, the SILENT third answer of F112 in the polarity witnesses.)

Gate: [`simulations/mixed_break_rule_gate.py`](../simulations/mixed_break_rule_gate.py),
writing `simulations/results/mixed_break_rule_gate.txt`, VERDICT last line, all
green 2026-08-16. Every number quoted here is from that gate or from a named
committed predecessor.

## What this is about

A coherence is one object with two faces, two spin arrangements read from
opposite sides, and the outside light separates the faces wherever they
disagree and passes through wherever they agree
([Dephasing Translated](../docs/quantum/DEPHASING_TRANSLATED.md)). The patterns
that oscillate at one shared frequency share one room, and the earlier notes
showed that a room whose patterns mix the two mirror families, the mirror-even
and the mirror-odd, is the only place where a pattern's weight can sit
unevenly between a site and its mirror partner. But measured at the fully
symmetric coupling, some of those mixed rooms still kept the even weight. This
note gives the missing piece, and it is about what a single site can and
cannot report. Because the fully symmetric coupling preserves the chain's
total spin, every pattern carries a spin label, and a single site's occupation
is a one-step instrument on that label: between two patterns whose total spin
disagrees by two steps or more, the site reports nothing. And in the
half-filled room, where up and down are in perfect balance, the site also
reports nothing between any two patterns of EQUAL total spin, because what it
would report is proportional to the imbalance, which is zero; the same
balance pins each pattern's own reading to exactly one half. A mixed room can
lose the even weight only where some cross-family pair of patterns is
reportable on both of its sides at once, and in every measured room that
suffices. Where every such pair contains a silence, the room provably keeps
the even weight, and the half-filled balance turns out to protect exactly one
direction of the frequency ladder. The break sizes are the products of the
two reports: exact fractions under a square root, never noise.
([Spin Translated](../docs/quantum/SPIN_TRANSLATED.md) is the door to which
spin is meant here: the chain's total spin, the label the fully symmetric
coupling preserves, not the η-ladder of the sideways arc.)

## What the repo already holds

The sweep, by store. `docs/ANALYTICAL_FORMULAS.md`: F154 owns the saturation
law; its containment covers the Δ = 1, N = 6 breaking spaces too (measured,
with slack), and what its entry fences is the DERIVATION status, a theorem
only on the all-scalar territory; F4 owns the Clebsch–Gordan stationary-mode
count and F50's max-spin sub-clause the Dicke transitions, the same arithmetic
in other objects; no entry states either zero here. `docs/proofs/`:
[PROOF_MIXED_SPACE_REFLECTION_LAW](../docs/proofs/PROOF_MIXED_SPACE_REFLECTION_LAW.md)
is the Δ = 0 sibling, the same question answered by a free-fermion pair-sum
lemma, a different engine on the same object;
[PROOF_UNIFORM_LAW](../docs/proofs/PROOF_UNIFORM_LAW.md) lives at ω = 0 behind
its own fence ("no degenerate eigenspace at ω ≠ 0 meets a class purely",
itself gated at N = 5..8, not derived), and this note lives exactly on the
other side of that fence, at ω ≠ 0;
[PROOF_CODIM1_BY_ADDITIVITY](../docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md) §7
holds the X^N object at block-index level, the Klein full flip and the
dressed antiunitary one-sided fold (its clause (c) scoping the latter to
Δ = 0), while the undressed one-sided cell fold F: (a,b) → (a, complement b)
is the endpoint note's fold section and the arc's ParkedAt; neither is the
face used here (this note uses the multiplet-parity face, X^N acting inside
the half-filled sector). `experiments/`:
[THE_ENDPOINTS_ARE_A_DENSITY_LAW](THE_ENDPOINTS_ARE_A_DENSITY_LAW.md) holds
the measured break/hold census on (1,3) and (2,4) with the numbers 0.129/0.258
this note derives; [THE_SPREAD_IS_A_RESONANCE](THE_SPREAD_IS_A_RESONANCE.md)
holds the fourth behaviour and, importantly, already owns the diagonal half of
zero (b): its X^N argument gives n_l = 1/2 exactly for every NON-degenerate
eigenvector at half filling. The OpenArcs registry: `compressed_density_laws`
NextStep (3) is this question verbatim; `sideways_spin_ladder` is the other
live Clebsch–Gordan arc, transport norms √(ℓ(ℓ+1) − m(m+1)) measured at
N = 7, the same su(2) arithmetic but a different algebra (the η-pairing
SU(2)), so the kinship is arithmetic, not identity.
`docs/GLOSSARY.md`: nothing for these objects. `fw.Confirmations` (24
entries): nothing. `docs/CAUGHT_ERRORS.md`: no CG bookkeeping entry, but the
lesson about threshold-defined quantities stands there, and this note walks into
that territory: four of the breaks below sit where a chosen 0.05 fence would
misread them, and only a floor classification shows them.

## The objects

Fix the open uniform chain at Δ = 1 (Pauli book H = Σ_bonds XX + YY + ZZ), a
coherence block (p, q) with the convention p ≤ q, and the site reflection R.
Because H is SU(2)-invariant and reflection-symmetric, each magnetization
sector p carries an orthonormal eigenbasis with sharp labels (E, S, r):
energy, total spin, reflection parity. At N = 6 every such label is
nondegenerate in every sector 1..5 (gate section A), so everything below is
basis-unique up to sign. The ad_H eigenspace at ω ≠ 0 on the block is spanned
by multiplet dyads |i⟩⟨j| with E_i − F_j = ω (E_i the ket-side, F_j the
bra-side energy), and a dyad's parity under R is
the product r_i · r_j; a mixed eigenspace is one carrying both parities. The
compressed site density between dyads factors into one-sector legs: writing
n_l for the site occupation and ⟨i|n_l|i'⟩ for a leg,

    ⟨(i,j)| N_l |(i',j')⟩ = −2 ⟨i|n_l|i'⟩ ⟨j|n_l|j'⟩          (i ≠ i', j ≠ j')
    ⟨(i,j)| N_l |(i,j')⟩  = ⟨j|n_l|j'⟩ (1 − 2⟨i|n_l|i⟩)        (i = i', j ≠ j')

## The rule, and where the zeros come from

A leg ⟨i|n_l|i'⟩ is DEAD iff

- **(a) the triangle zero:** |S_i − S_i'| ≥ 2. Off the diagonal only the Z_l
  half of n_l contributes (the identity half needs i = i'), and Z_l is the
  q = 0 component of a rank-1 spherical tensor (the site spin), so the
  Wigner–Eckart triangle kills the element.
- **(b) the half-filling zero:** the sector is half filling (magnetization
  m = 0) and S_i = S_i'. The Clebsch–Gordan factor in the fixed-m sector is
  ⟨S m; 1 0 | S m⟩ = m/√(S(S+1)), zero at m = 0. This zero, and only this
  one, has a second route: Z_l anticommutes with X^N, X^N maps the m = 0
  sector to itself with eigenvalue i^N (−1)^S on a spin-S state, which is
  (−1)^(S+1) at N = 6, so same-S states share the parity and an X^N-odd
  operator cannot connect them (that route kills every even ΔS at m = 0,
  more than the rule needs). The same zero makes every diagonal exact,
  ⟨i|n_l|i⟩ = 1/2, which kills the i = i' factor (1 − 2⟨i|n_l|i⟩).

**The rule:** a mixed eigenspace breaks C_l = 0 iff some cross-parity dyad
pair has BOTH legs alive. The two directions carry different weight. The HOLD
direction is a theorem: R-conjugation makes the C_l element of a cross-parity
dyad pair equal to twice its comp(N_l) element and cancels every same-parity
pair identically, so if each cross-parity pair contains a dead leg, C_l = 0
exactly. The BREAK direction is genericity: "alive" means "not forced to zero
by the two symmetries", a leg alive by the rule could still vanish at every
common site, and the census answer is that in all 60 mixed spaces none does.
Verified 60/60 over ALL fifteen N = 6 blocks p ≤ q, and element-exactly:
every rule-dead compressed element measures at the eigensolver floor (worst
1.1e−16), every hold is a floor-hold (worst 5.6e−16), and the dyad route is
bridged per eigenspace against an independent cells-route compression (gate
section B). Both zeros are exact in exact arithmetic, so the float gate
carries an error model, not a chosen number: eigh backward error,
eps·‖H_sector‖/gap, with the classification living in the empty band of
nearly fourteen decades between the holds and the smallest break.

## The census, and the sizes nobody asked for

All 60 mixed spaces, with resC the worst spectral norm of the compressed
C_l = comp(N_l − N_{N−1−l}) and each break given as 2160·resC² (an integer,
gated; the twelve holds, and the four sub-fence breaks discussed below, in
bold):

| block | ω = −4 | ω = −2 | ω = +2 | ω = +4 | m = 0 side |
|-------|--------|--------|--------|--------|------------|
| (1,1) | 20 | 80 | 80 | 20 | none |
| (1,2) | 32 | 128 | 20 | **5** | none |
| (1,3) | 36 | 144 | **hold** | **hold** | bra, protects ω > 0 |
| (1,4) | 32 | 128 | 20 | **5** | none |
| (1,5) | 20 | 80 | 80 | 20 | none |
| (2,2) | 8 | 32 | 32 | 8 | none |
| (2,3) | 9 | 36 | **hold** | **hold** | bra, protects ω > 0 |
| (2,4) | 8 | 32 | 32 | 8 | none |
| (2,5) | **5** | 20 | 128 | 32 | none |
| (3,3) | **hold** | **hold** | **hold** | **hold** | both, protects everything |
| (3,4) | **hold** | **hold** | 36 | 9 | ket, protects ω < 0 |
| (3,5) | **hold** | **hold** | 144 | 36 | ket, protects ω < 0 |
| (4,4) | 8 | 32 | 32 | 8 | none |
| (4,5) | **5** | 20 | 128 | 32 | none |
| (5,5) | 20 | 80 | 80 | 20 | none |

Reading the table: resC = √(k/2160), so (1,3)'s breaks are √(1/60) = 0.129099
and √(1/15) = 0.258199, exactly the endpoint note's 0.129/0.258. The sizes
are a measurement, not a derivation: the legs entering these particular mixed
spaces measure as square roots of rationals and the breaking elements are
their products, but aliveness does not imply rationality (generic alive legs
in the same sectors measure with irrational squares, a contrast checked off
the gate), and the common denominator 2160 is a fit. The four bold **5**
entries are what a threshold would have hidden: (1,2) and (1,4) at ω = +4,
(2,5) and (4,5) at ω = −4 break at
resC = √(1/432) = 0.048113, thirteen decades above the floor and just below
the 0.05 fence the earlier census used. No committed claim is affected: the
earlier census's blocks, (1,3) and (2,4), carry no break in that band, and
had the fence been walked onto these four blocks it would have misread them;
the honest fragility on the old blocks is that (2,4)'s smallest break,
0.060858, clears 0.05 by only 22%. Classification here is therefore by the
floor: holds ≤ 64 × floor, breaks ≥ 0.04, the band between gated empty.

Two structural echoes, stated as observations, not gated: the sizes depend
only on the pair (min(p, N−p), min(q, N−q)) up to the direction of ω, so the
fifteen blocks fall into six classes (the maps behind that reduction are the
per-side spin flip p → N−p at fixed ω and the adjoint (p,q) → (q,p) with
ω → −ω, whose composition pairs the table's mirror rows); and within each
block the dim-3 break size k is 4× the dim-2 break size of the same ω side,
which is 2× in resC.

## The protection map

Each half-filled SIDE of a block protects one ω sign: the bra side (q = N/2)
protects ω > 0, so (1,3) and (2,3) hold there; the ket side (p = N/2)
protects ω < 0, so (3,4) and (3,5) hold there; (3,3) carries both sides and
holds at all four ω; and the ten blocks without an m = 0 sector hold
nowhere (gate section C). The sign asymmetry is the rule plus sector
spectroscopy: at m = 0 the half-filling zero leaves only ΔS = 1 legs alive,
and in all 60 measured spaces the only alive S-pair at m = 0 is (2,3), whose
S = 3 member is the ferromagnetic multiplet, the global energy maximum E = 5
(this "every alive m = 0 leg touches S = N/2" fact is gated). A dyad carrying
it therefore has ω = E_i − 5 ≤ 0 when the m = 0 sector is the bra side and
ω = 5 − F_j ≥ 0 when it is the ket side: the protected sign is the one those
dyads cannot reach. (That no other ΔS = 1 pair appears in the mixed spaces is
census, not derivation.) This is the same boundary the resonance note's
fourth behaviour drew, an index at N/2 as the precondition for the X^N
cancellation; here it returns one level up, as the multiplet parity.

## What is ours and what was owned

The measured break/hold map on (1,3) and (2,4), with 0.129/0.258 and the
"mixedness is necessary but not sufficient" sentence, is the endpoint note's.
The diagonal half of zero (b), n_l = 1/2 at half filling, is the resonance
note's, proven there by the same X^N argument for every non-degenerate
eigenvector. What this note adds: the off-diagonal zeros (the triangle, and
the m = 0 same-S zero that needs no nondegeneracy hypothesis), the
factorization into legs, the hold-direction theorem, the 60-space census
with the floor classification and the four sub-fence breaks, the measured
√rational sizes, and the protection map with its mechanism.

## Bond control and the frontier

The 60 mixed spaces are a uniform-chain resonance, as the Δ = 0 collisions at
N = 8 were: a fixed dyadic palindromic bond profile and a seeded six-profile
random dyadic palindromic sweep both leave zero mixed spaces on all fifteen
blocks, with [R, H] = 0 asserted exactly per profile (gate section E). And at
Δ = 1 there are no mixed spaces at ω ≠ 0 at any N in 3..8 except N = 6, on
all blocks p ≤ q per size (N = 7: 21 blocks, N = 8: 28): among N ≤ 8 the
phenomenon lives at N = 6 alone. Both statements are measured fences, not
derivations; no claim is made past the tested sizes.

## What stays open

The second half of NextStep (3) is untouched: why the breaking spaces'
compressed spectra stay INSIDE the size-class-centre interval on the locus,
containment measured with slack, not derived. The break direction of the rule
is genericity, not a theorem: nothing forbids an alive cross-parity pair
whose legs miss each other on every site, the census just contains none. The
√rational sizes and the 2160 are measured, unexplained. And whether the two
zeros organize the N = 11 frontier the arc's NextStep (4) opens (mixed spaces
are abundant there at Δ = 0) is not asked here: this note is Δ = 1 only.
