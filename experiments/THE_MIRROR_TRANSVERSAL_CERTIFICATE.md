# The mirror-transversal certificate

**2026-08-16.** The second half of arc `compressed_density_laws` NextStep (3),
open since the endpoint note: on the mirror-balanced locus
γ_l + γ_{N−1−l} = 2γ̄, the eigenspaces that BREAK the reflection law C_l = 0
at Δ = 1, N = 6 were measured to keep their compressed dissipator spectra
INSIDE the size-class-centre interval [−2γ̄·s_max, −2γ̄·s_min], with slack, on
one tested interior profile, and nothing derived it. This note closes it for
EVERY physical locus profile at once: the physical locus set is a box whose
corners are the mirror transversals (one site of each mirror pair lit at
2γ̄, its partner dark), extreme eigenvalues of an affine Hermitian family
are attained at corners, and at a corner the compression splits along the
mirror parity into a slider-blind even part and an odd part whose size a
closed form controls. The result is the reduction DERIVED plus a finite
certificate, closed-form on all 48 breaking spaces (exact in exact
arithmetic, gated at the eigensolver floor), and the certificate is sharp:
on the six dim-3 spaces of the corner blocks (1,1), (1,5), (5,5) it holds
with exact equality (1/3 = 1/3), a transversal eigenvalue reaching an
interval edge precisely, so the containment cannot be improved; the slack
the endpoint note measured is the margin the closed form leaves everywhere
else.

Gate: [`simulations/locus_containment_gate.py`](../simulations/locus_containment_gate.py),
writing `simulations/results/locus_containment_gate.txt`, VERDICT last line,
all green 2026-08-16. Every number quoted here is from that gate or a named
committed predecessor; the breaking census is
[THE_TWO_SPIN_ZEROS](THE_TWO_SPIN_ZEROS.md)' (48 breaking spaces, sitting on
14 of the fifteen-block census's blocks; (3,3) carries only holds).

## What this is about

The earlier notes established where the decay rates of same-frequency
pattern rooms must land when the lighting is mirror-balanced: inside the
window spanned by the smallest and largest disagreement count, because a
room that keeps its weight mirror-evenly cannot tell balanced lighting from
uniform lighting. The rooms that mix the two mirror families can lose the
even weight ([THE_TWO_SPIN_ZEROS](THE_TWO_SPIN_ZEROS.md) decides which do),
and for those rooms the argument breaks; yet measured, their rates stayed
inside the window too. This note explains that, and the explanation has the
shape of a walk to the corners. Balanced lighting is a family with a slider
per mirror pair: slide light from a site to its partner and the balance is
kept. The most extreme balanced lightings sit at the ends of the sliders,
where one site of each pair carries all the pair's light and its partner
stands fully dark (an image this note introduces; the repo had no word for
it). The worst rate over the whole family always leans to the ends of the
sliders, never to their middles, so the worst case over ALL balanced
lightings is one of these corner lightings, and at N = 6, with three
sliders, there are only eight of them. At a corner, what a room pays
splits into two halves: the half that counts disagreements evenly across
each pair, which is blind to the sliders altogether, and the uneven half,
which is exactly the object the two spin zeros left alive. The two halves
do not simply add: they combine crosswise, like the two legs of a right
triangle, and what must fit into the margin the even half leaves to the
window's edge is the triangle's long side. The certificate is that check,
room by room and corner by corner; on the corner blocks the long side
uses the margin exactly, to the last drop, and holds.

## What the repo already holds

The sweep, by store. `docs/proofs/`:
[PROOF_CODIM1_BY_ADDITIVITY](../docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md)
§6 owns the rate WINDOW (Bendixson: real parts inside the range of the cell
rates, exact for every eigenvalue), the window-edge lemma (an eigenvalue AT
an edge forces a joint eigenvector), and the window combinatorics, all at
uniform γ; it contains no parameter-extremality argument, and reading its
window under a profile as the cell-rate range, wider than the centre
interval certified here, is this note's extrapolation, not §6's own
statement. [PROOF_F101](../docs/proofs/PROOF_F101_C1_MIRROR_GAMMA_PARITY.md)
is this note's direct ancestor: it splits the profile γ = γ_sym + γ_anti
along the site mirror and proves an exact parity in γ_anti; the split below
is the same decomposition applied to the compression.
[PROOF_MIXED_SPACE_REFLECTION_LAW](../docs/proofs/PROOF_MIXED_SPACE_REFLECTION_LAW.md)
contains no γ at all (its object is ad_H). `docs/ANALYTICAL_FORMULAS.md`:
F154 owns the containment where C_l = 0 (theorem, any locus profile) and
names this note's territory in its own words, "the breaking interiors
measured to stay inside"; F153's pinning needs uniform γ (bar its four
one-cell blocks); F91 owns the
locus and the sharper indexed pair-sum law, answered below rather than
contradicted; F140's frozen divisor is the locus neighbour with the same
pair bookkeeping; F64 holds the single-lit-site atom (Z-dephasing on one
site B, any graph topology). `experiments/`:
[RESONANT_RETURN](RESONANT_RETURN.md) holds the repo's only prior
corner-of-the-γ-set idea: it poses the concentrator as a constrained
optimization over γ profiles and asserts the optimum is "a delta function
at the boundary" (which the document's own restatement places at the
chain's edge qubit, a spatial concentration), explicitly leaving "a
global-optimum proof over
all γ profiles" open; this note's reduction is the first such argument
carried through, on a different set (the locus box, corners =
transversals, not single-site spikes).
[PROOF_RING_N4_DIHEDRAL_LOCK](../docs/proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md)
holds the repo's only computed partially-zero per-site profile, as a GATED
identity (γ = (2,0,0,0) matches the single-jump-operator spectrum), with
the sentence that a zero rate removes that site's jump operator
altogether: a corner lighting is exactly that, still completely positive
(rates 0 and 2γ̄, both ≥ 0), a boundary point of the closed box, with
⌊N/2⌋ sites entirely undephased. The OpenArcs registry:
`compressed_density_laws` carries the question verbatim;
`site_resolved_vacuum_block` already walks the convex-combination /
Bendixson move on the (0,1) block and records the refutation "a profile
opens the floor line" (false in general), which this note does not
re-walk. `docs/GLOSSARY.md` and `docs/quantum/THE_LABEL_MAP.md`: no entry
for a dephasing profile, the locus, or any lit/dark image; the nearest
house word is the Concentrator (one site absorbs the whole budget), and a
mirror transversal is its paired cousin, one lit site PER PAIR.
`fw.Confirmations`: nothing. `docs/CAUGHT_ERRORS.md`: the 2026-08-06
lesson (a fitted number in the vocabulary of a law) and the
threshold-defined-quantity lesson both apply and are obeyed below: every
fence here is an enumerated finite minimum stated as such, a gated
exactness, or a closed form.

## The reduction: corners of the locus box

comp(D) = −2 Σ_l γ_l comp(N_l) is affine in the profile, and the
compressions comp(N_l) are profile-independent because the eigenspace Ω is
built from ad_H alone; without that the argument would collapse, so it is
said. The physical locus set {γ_l ≥ 0, γ_l + γ_{N−1−l} = 2γ̄} is a product
of ⌊N/2⌋ segments, one slider per mirror pair (positivity of both partners
forces γ_l ≤ 2γ̄), and its corners are the 2^⌊N/2⌋ mirror transversals A:
one site per pair at 2γ̄, its partner at 0. For a Hermitian affine family,
λ_max is convex in the parameters and λ_min is concave (the variational
characterization), so over the box λ_max peaks at a corner and λ_min
bottoms at a corner: two-sided containment for every physical locus
profile is EQUIVALENT to containment at the transversals. At even N, this
note's case, comp(D) = −4γ̄·comp(N_A) there; at odd N the middle site is
its own mirror, stays pinned at γ̄, and adds the fixed term
−2γ̄·comp(N_mid) to every corner, so the reduction survives verbatim while
the corner form carries the middle term. The step is an argument, not a
computation, and it is block-general; gate section (B) cross-checks its
consequence on 20 random locus profiles per space (worst overshoot beyond
the transversal hull 0.0, tightest distance to the hull 4.4e−3). Two fences,
stated against their neighbours: the convexity used is of the EXTREME
eigenvalues of a Hermitian family only, not of per-mode rate curves, which
[LIGHT_DOSE_RESPONSE](LIGHT_DOSE_RESPONSE.md) (superseding
[GAMMA_AS_BINDING](GAMMA_AS_BINDING.md)'s per-sector reading) measured to
be genuinely nonlinear in the global scale; and off the locus box (γ̄ = 0
gain profiles, the two-sources territory) nothing here applies.

## The split: F91's parities, compressed

Let R̃ be the compressed site reflection on a breaking eigenspace Ω (an
involution, since [ad_H, R] = 0 and R preserves Ω; in the multiplet-dyad
basis it is the diagonal of dyad parities). For a transversal A,
R̃·comp(N_A)·R̃ = comp(N_{A^c}) and N_A + N_{A^c} = N_XY, so

    comp(N_A) = ( comp(N_XY) + Σ_{l<N/2} σ_l C_l ) / 2,

σ_l = ±1 the transversal's choices, C_l = comp(N_l − N_{N−1−l}) the
two-spin-zeros object. The even part comp(N_XY) commutes with R̃ and every
C_l anticommutes with it (both gated), so the C_l are supported on
cross-parity dyad pairs only. This answers F91's tension rather than
contradicting it: F91's sharper law makes the F71-refined diagonal blocks
of L depend only on the indexed pair-sums γ_l + γ_{N−1−l}, with the
pair-difference dependence confined to the F71-cross-blocks, and F91
insists its diagonal-block quantities are not decay rates of anything; the
carrier here is a different object, the compression onto an ad_H
eigenspace, and it obeys the same grammar: the pair-sum part is
−2γ̄·comp(N_XY), blind to the sliders as F91 demands on the locus, and the
slider dependence enters only through the C_l, nonzero precisely because a
breaking Ω carries cross-parity content. Measured structure of the even
part, gated in section (C): comp(N_XY) is DIAGONAL in the multiplet-dyad
eigenbasis on every breaking space, with /24-rational entries (5/3, 11/6,
7/3, 29/12, 8/3, 65/24, 3, 10/3, 79/24, 11/3, 43/12, 25/6, 13/3, the
sorted set gated); it is scalar on the 24 dim-3 spaces and the 4 dim-2
spaces of the m = 0 blocks, so the diagonality statement has content
exactly on the remaining 20 dim-2 spaces, where the two diagonal entries
genuinely differ, with gaps exactly 1/24, 1/12, 1/6 (gated as integer
gaps on the 24ths, no float fence). The diagonality and the rationals are
measured facts, not derived. Two pieces of the interval arithmetic, by
contrast, ARE exact combinatorics. First, a (p,q) cell has
popcount(a⊕b) ≡ p+q (mod 2), so the block's class set is the ladder
{|p−q|, |p−q|+2, ..., min(p+q, 2N−p−q)}, giving s_min = |p−q| and
s_max = min(p+q, 2N−p−q) with no eigenvector, no mass sum and no
threshold; off-ladder classes are empty at block level (no discriminating
power there), and the measured content is that every breaking space
carries mass on EVERY ladder class (gated ≥ 0.1, measured minimum 0.119):
the support is the full ladder, endpoints live, non-contiguous by parity.
Second, for any block with an index at N/2 the ladder is centred at 3
identically ((3−p + p+3)/2 on the q = N/2 blocks, (q−3 + 9−q)/2 on the
p = N/2 blocks), so the interval CENTRE is derived; that the measured
scalar m̄ sits exactly there, m̄ = 3, is a property of comp(N_XY) and is
gated as a measurement. Where comp(N_XY) = m̄·Id, the spectrum of
comp(N_A) is symmetric about m̄/2, because the odd part anticommutes with
the involution R̃.

## The certificate

For every one of the 48 breaking spaces and every one of its 8
transversals, spec(comp N_A) ⊂ [s_min/2, s_max/2], gated at the eigensolver
floor (section D), and the whole check is one closed form. With
δ = half the gap of the comp(N_XY) diagonal (zero on the 28 scalar
spaces), w = max_σ ‖Σ_l σ_l C_l‖ the adversarial odd norm, and
bound = min(s_max − m̄, m̄ − s_min):

    √(δ² + w²) ≤ bound,

with the measured worst transversal slack EQUAL to (bound − √(δ² + w²))/2
on every space (gated at the floor): the closed form IS the spectrum, not
an estimate of it. On the scalar spaces it reduces to w ≤ bound; on the
dim-2 spaces it is the exact 2×2 eigenvalue, and it has to be, because the
cruder Weyl bound δ + w FAILS on the six corner dim-2 spaces (0.275783
against the bound 0.25, gated as failing) while the exact form clears it
(√(1/144 + 1/27) = 0.209718, slack 0.020141). Both ingredients are
rational under a square: 2160·w² is an integer on every space (the same
measured denominator as the two-spin-zeros sizes, a fit, not derived) and
2304·δ² is an integer, both gated. The tight set is exactly the six dim-3
spaces of the corner blocks (1,1), (1,5), (5,5) at ω = ±2: there w = 1/3
equals the bound and a transversal eigenvalue reaches an interval edge
exactly, the UPPER edge s_max/2 on (1,1) and (5,5) and the LOWER edge
s_min/2 on (1,5) (edge sides gated), the same shape as §6's window-edge
lemma, an endpoint reached structurally rather than by accident. Because
slack = (bound − √(δ² + w²))/2 exactly, tight versus non-tight is an
exact split, not a band: the six corner dim-2 spaces sit at slack
0.020141 and the remaining 36 at ≥ 0.122515, enumerated finite values.
Together with the theorem side (C_l = 0 forces containment wherever the
space holds, F154, any locus profile), the containment of the compressed
dissipator at Δ = 1, N = 6 is now settled for every physical locus
profile on every ω ≠ 0 mixed eigenspace: the reduction and the split
derived, the closed-form inequality finitely certified, exact in exact
arithmetic and gated at the floor.

## What is ours and what was owned

The containment-with-slack measurement on one interior profile is the
endpoint note's; the Bendixson window and the window-edge tightness shape
are PROOF_CODIM1 §6's; the γ_sym/γ_anti split is PROOF_F101's; the
corner-of-the-γ-set idea appears once before as RESONANT_RETURN's
undischarged optimization sketch; the C_l and the breaking census are
THE_TWO_SPIN_ZEROS'. What this note adds: the locus box and its mirror
transversals as the reduction that makes a profile-family statement
finite, the discharged convexity argument, the compressed parity split
with its measured diagonal rational even part, the exact ladder form of
the interval (s_min = |p−q|, s_max = min(p+q, 2N−p−q), the m = 0 centre
derived), the closed-form certificate √(δ² + w²) ≤ bound with slack
equality and exact tightness on the corner class, and the closing of
NextStep (3).

## What stays open

The diagonality of comp(N_XY) in the dyad eigenbasis, its /24-rational
entries, and the 2160/2304 denominators are measured, underived. The
corner tightness invites the window-edge treatment (which structural
eigenvector reaches the edge) but that reading is not written. Everything
here is Δ = 1, N = 6, uniform J, open chain, on the physical locus box;
the N = 11 frontier of NextStep (4) and the F154 typed carrier of (5)
remain the arc's open items.
