# The endpoints are a density law

**2026-08-15.** [THE_SPREAD_IS_A_RESONANCE](THE_SPREAD_IS_A_RESONANCE.md) left one
question at the top of its open list: why do the full resonances saturate on EXACTLY
the size-class-centre interval, endpoints attained? The answer has two layers, and the
first is a theorem. On every degenerate ad_H eigenspace where the chain reflection
acts as a scalar ±1, the compressed site density is reflection-symmetric BY SYMMETRY,
and the R₉₀ locus pairs γ_l + γ_{N−1−l} = 2γ̄ exactly against that, so containment on
the size-class centres and endpoint attainment follow with no input beyond [H, R] = 0
and the locus. The parity census then turns the theorem into the full answer for
almost every measured case: at Δ = 0 every colliding eigenspace on every census
block is scalar up to N = 7, and at Δ = ½ up to N = 8. The second layer is a law with real content: the
parity-MIXED eigenspaces that first appear at N = 8, Δ = 0 STILL carry the
reflection-symmetric compressed density, at the eigensolver floor, where symmetry
alone permits a cross block of order one; and at Δ = 1 those mixed spaces are exactly
where the law can break. The SU(2) resonance adds its own second law: the compressed
density on every pure size class is exactly UNIFORM, s/N per site, blind to the whole
γ profile.

Gate: `simulations/high_q_selection_gate.py`'s sibling
[`simulations/endpoint_density_gate.py`](../simulations/endpoint_density_gate.py),
writing `simulations/results/endpoint_density_gate.txt`, VERDICT last line, all green
2026-08-15. Every number quoted here is from that gate or from a named committed
predecessor. Registered 2026-08-15: the saturation law is
[F154](../docs/ANALYTICAL_FORMULAS.md#f154), and the scope-extended Mechanism went into
[F122](../docs/ANALYTICAL_FORMULAS.md#f122)'s own entry rather than a second number.

## What the repo already holds

The sweep, by store. `docs/ANALYTICAL_FORMULAS.md`: F122's registered Mechanism is
this compression at UNIFORM γ (high-Q degenerate PT; the literal P_Ω N_XY P_Ω form is
PROOF_STRUCTURAL_CEILING's §1); this note is its lift onto a γ profile. F153 owns the
involution source the theorem layer lifts (γ anti-palindromic, source 2, hypothesis H
reflection-symmetric) and also owns the fence the uniform law must not be misread
against: F153's fence that one size class alone does not suffice is ENTRY-WISE (the
cell rates differ within a class under a profile), while the uniformity here is a
statement about the COMPRESSION onto an eigenspace, a different object; both are
true at once.
F91 owns the locus, F140 the corner-block fence, and F144 is the genre precedent
(extremes of a compressed disagreement on ker ad_h, derived by an SU(2) ladder).
`docs/proofs/`: PROOF_R90_FROZEN_DIVISOR owns the |S|-resolved reflection reading
(on the rates, rate ↦ 2γ̄|S| − rate; in λ, "the center is −2γ̄|S|"); PROOF_ABSORPTION_THEOREM Theorem 2
owns the per-site split D = −2 Σ_l γ_l N_l that makes the site-resolved question
well-posed (its light_l is a Rayleigh quotient on one eigenvector, with a trace-form
corollary on slow subspaces of L; the compressed-operator statement here is new);
PROOF_FROZEN_BAND_SO4 §6 is the nearest relative, the compressed double-occupancy
D̂ solved in closed form on the ad_h kernel of the one-excitation rung at Δ = 0
(G = (1/M)(𝟏𝟏ᵀ + (I+R)/2)); its D̂ is the AGREEMENT (the rung-1 disagreement is
𝒦 = 1 − D̂) and its R is the chiral involution ON MODES, not the site reflection
this note is about, so the two must not be merged; PROOF_CODIM1_BY_ADDITIVITY §6
owns the Bendixson rate window over CELL rates (at uniform γ it coincides with the
centre interval; under a profile the centre interval is properly inside it, and the
content here is not the width but that the endpoints are the CENTRES), and its §7
fold lattice owns clause (b), the s-symmetric self-folded block, which this note
uses for the off-locus overshoot equality below.
`experiments/`: the predecessor's saturation table and its necessity/sufficiency
fence; [WHAT_THE_R90_LOCUS_BUYS](WHAT_THE_R90_LOCUS_BUYS.md) measured the endpoints
(its gate is `simulations/two_sources_gate.py`) and withdrew the "never flattens"
HALF of a drafted mechanism, while its first half, saturation on the centre interval,
stayed measured and unexplained and is what this note derives;
`review/OBC_SINE_BASIS_FINDINGS.md` owns the single-particle reflection parities
R ψ_k = (−1)^(k+1) ψ_k; SURVIVOR_FLIP_AND_REFLECTION_ODD's reflection-ODD density is
a signed zero-mean fluctuation of one survivor mode, not an occupation, and is not
this object; F89_MULTI_SECTOR_MONODROMY owns the quantized density overlap I(a,b) at
uniform γ, aggregated, never site-resolved. The OpenArcs registry: three arcs hold
thread content, `xxz_axis_handover` (the Δ axis), `topology_band_edge` (retired;
the registry's own carrier of the uniform-γ compression rule) and
`f_registry_meets_the_typed_layer` (F153's clauses); `site_resolved_vacuum_block` is
the genre neighbour for per-site density readings and carries the caution that
biorthogonal right-eigenvector densities are not probabilities (the compressions here
are on ORTHONORMAL ad_H eigenbases, the safe route); no arc is dedicated to this
thread, and the change that lands this note opens one (`compressed_density_laws`).
`fw.Confirmations`: nothing (no hardware content). `docs/GLOSSARY.md`: no entries for
the disagreement set or size classes, consistent with the predecessors' sweeps.
Nothing in any store states either law or the parity theorem.

## The objects

Vocabulary as in the predecessor's head. Additions: **N_l** is the diagonal indicator
that a cell's ket and bra disagree at site l, so the dissipator's site-resolved face
is D = −2 Σ_l γ_l N_l and F122's size operator is N_XY = Σ_l N_l (eigenvalue
|S| = popcount(a⊕b) on a cell). **R** is the site reflection, acting on cells as
(a, b) ↦ (rev a, rev b); it commutes with the block's ad_H exactly when H is
reflection-symmetric (palindromic bonds), and then preserves every degenerate
eigenspace **Ω** of ad_H. Π_Ω is the orthonormal projector onto Ω, comp(X) = Π_Ω X Π_Ω,
and **Ω_s** = Ω ∩ span of the class-s cells. An eigenspace is **scalar** when R
restricted to it is ±Id (all its R-eigenvalues one sign) and **mixed** otherwise.
Compression is on an orthonormal basis throughout; nothing here touches the
biorthogonal projectors of the non-normal L.

## Layer 1: the theorem, and the parity census that makes it the answer

Let C_l = Π_Ω (N_l − N_{N−1−l}) Π_Ω, the compressed odd density. C_l is R-odd
(R N_l R = N_{N−1−l}), so on a scalar Ω conjugation gives C_l = −C_l = 0: **on every
scalar eigenspace the reflection-symmetric compressed density is forced by
[H, R] = 0 alone.** On a mixed Ω the same argument kills only the two diagonal parity
blocks of C_l; the cross block is unconstrained by symmetry.

The census (gate section A, spectral-norm residuals against the error model
res ≤ 64·eps·‖A‖/gap per eigenspace): 30 case rows over N = 4, 5, 6, 7, 8,
Δ ∈ {0, ½}, eight blocks at N = 6, up to 143 eigenspaces per block and eigenspace
dimensions up to 72. **On every census row at Δ = 0 up to N = 7, and on every
Δ = ½ row including N = 8, every colliding eigenspace is scalar**, so on that
territory the strong law IS the theorem (the scoping matters: at Δ = 1, N = 6
mixed eigenspaces DO exist, four per block tested, the subject of Layer 2 below).
The mixed Δ = 0 eigenspaces first appear at N = 8: 20 of 59 on (1,3), 21 of 83 on
(1,4), 50 of 143 on (2,4).

## Layer 2: the two laws, measured

**The reflection law beyond the symmetry (Δ = 0, the free-fermion resonance).** On
all 91 mixed eigenspaces at N = 8, C_l = 0 still holds at the floor (max spectral
norm 6·10⁻¹⁵), where symmetry permits a cross block of order one: a generic R-odd
diagonal compressed on the same eigenspaces has norm 0.119. This is the note's
measured content at Δ = 0, and it is deliberately scoped: it rests on three blocks
at one N and one Δ, because that is where mixed collisions exist at all.

**The uniform law (Δ = 1, the SU(2) resonance).** On every Ω_s:

    comp(N_l) = (s/N) · Id        for every site l.

Gated on all eight N = 6 blocks and on (1,3), (1,4) at N = 8. Blind to the γ profile
by construction. At Δ = 0 this law is FALSE (max residual 0.226, gated from below):
the two resonances carry genuinely different density structure, mirrored vs uniform.

**Where the reflection law breaks, the parity census says why.** At Δ = 1, N = 6
every eigenspace that breaks C_l = 0 is MIXED, and every scalar one sits at the floor
(forced): on (1,3) the breaks are ω = −4 (dim 2, (1+,1−), res 0.129) and ω = −2
(dim 3, (1+,2−), res 0.258), while the mixed spaces at ω = +2 (dim 4, (2+,2−)) and
ω = +4 (dim 3, (1+,2−)) HOLD at the floor, so mixedness is necessary but not
sufficient; on (2,4) all four mixed spaces break (ω = ±2, ±4). The ω = 0 carriers of
the extremes are scalar on both blocks, and at N = 8 (1,3) and (2,4) every Δ = 1
eigenspace is scalar, so there the reflection law is entirely the theorem.

## The consequence: saturation derived

Sum rule plus C_l = 0 is the whole mechanism. On the locus γ_l + γ_{N−1−l} = 2γ̄,
pairing the sites in Σ_l γ_l comp(N_l) gives γ̄ comp(N_XY), so on every eigenspace
where C_l = 0 (forced or lawful)

    Π_Ω D Π_Ω = −2γ̄ · Π_Ω N_XY Π_Ω.

Gated as an identity (ratios ≤ 0.003 of the model floor) on (1,3) and (2,4) at
Δ = 0, (1,3) at Δ = ½, both at Δ = 1 on the law-abiding spaces, and (1,4) at N = 8.
Three consequences, each gated:

1. **Containment.** By Rayleigh on comp(N_XY), the compressed spectrum of Ω lies in
   [−2γ̄·s_max(Ω), −2γ̄·s_min(Ω)], inside the block's size-class-centre interval.
   Gated for EVERY multi-dim eigenspace of the six locus cases, the Δ = 1 breaking
   spaces included (their spectra sit inside with slack; that containment is
   measured, not derived).
2. **Attainment.** A pure-class vector v ∈ Ω_s satisfies N_XY v = s v, hence
   comp(N_XY) v = s v, so each nonempty extreme class pins an endpoint ON its
   centre; both endpoints gated where the extreme classes have pure vectors. Block
   level: (1,3) at N = 6 reads range [−4, −2] at Δ = 0 AND Δ = 1, and (1,4) at
   N = 8 reads [−5, −3] at Δ = 0, spread 2.0000 (the predecessor reads the same
   range at Δ = ½), all to the model floor. One honesty clause: those two headline blocks are SELF-FOLDED (an index
   at N/2), where the fold below makes the two endpoints one fact, so the
   genuinely two-sided attainment case among the locus rows is (2,4) at N = 6
   (block range [−6, −2] gated to the floor), whose lower centre −2γ̄·N is in turn
   an arithmetic endpoint (the global minimum of the rate diagonal for any
   profile); its upper endpoint is the one the density law earns.
3. **Collapse is the single-class case.** On a pinned block every cell has
   |S| = |p−q|, N_XY is scalar there, so the same identity forces
   Π D Π = −2γ̄|p−q| · Id: the locus collapse the predecessor measured on (0,3)
   "although the same collisions are present" is this law, not an exemption from it.

The derivation status is worth stating plainly, because it is stronger than the
first draft of this note knew. **Wherever every colliding eigenspace is scalar, the
locus containment is DERIVED** (theorem + locus pairing), and attainment follows
from the measured pure-vector census on top. On the blocks in the census, that
territory is: all Δ = 0 rows through N = 7, all Δ = ½ rows through N = 8, the
Δ = 1 rows at N = 8 ((1,3) and (2,4), gate section C), and the extreme-carrying ω = 0
spaces at Δ = 1, N = 6 (with the breaking spaces' interior containment measured,
not derived). At N = 8, Δ = 0 the saturation rests on the measured mixed-space law
and is exact but not yet derived.

And the sizes of the resonances stop being accidents. At Δ = ½, N = 6, NO colliding
eigenspace of (1,3) contains a pure-class vector (gated: zero of six), so nothing
can attain the centres and the resonance stays small (the predecessor's 0.4494),
while (2,4) at the SAME Δ carries pure vectors of both its extreme classes in its one
colliding eigenspace and saturates fully, range [−6, −2], spread 4.0000 (gated in
section D): same N, same Δ, opposite outcome, decided by the census alone; the
(1,4)/(2,4) magnitude contrast at N = 8 has the same shape, one dim-16 ω = 0
eigenspace on one block against rank-2 spaces on the other (the predecessor's
2.0000 vs 0.4467, measured there; that the dim-16 space is the one carrying
pure-class vectors is a local-scout reading, not gated). The interior of a mixed-class
eigenspace is comp(N_XY) arithmetic: on the dim-12 ω = 0 space of (1,3) at Δ = 0,
7·spec(comp N_XY) = {14, 20, 22, 28}, each threefold, integrality defect 4·10⁻¹⁴, so
the compressed rates on the locus are −2γ̄ · {2, 20/7, 22/7, 4}, a RATIONAL ladder
between the centres.

One bridge sentence the compression owes the reader: everything above lives at
J = ∞. The identification of spec(Π_Ω D Π_Ω) with the block's actual Re-spectrum at
strong coupling is F122's high-Q degenerate PT, and the predecessor's gate section
(A) established the truncation as O(1/J²); the saturation numbers it measured on
its J ladder (10³ to 10⁶) carry that error, not zero.

## The fences

**The bond hypothesis is inherited, and the control separates what exists to
separate.** The theorem needs H reflection-symmetric (F153's source-2 hypothesis).
At Δ = 0 the bipartite chiral symmetry of the hopping matrix (ε ↦ −ε for any
bonds) keeps difference collisions alive for ANY bond profile, but the MIXED
collisions are a resonance of the uniform chain: the fixed palindromic J and a
ten-profile random palindromic sweep at N = 8 all leave zero mixed eigenspaces
(only the forced layer remains, and it holds at the floor), and a non-palindromic
J breaks the law outright (0.456) on the same collision dimensions, R being no
symmetry there. Gated. So the mixed-space content is, so far, a uniform-chain
statement, and the non-palindromic arm fences the hypothesis rather than the
cross-block content.

**Off the locus the consequence is gone, and the two ends are not two facts.** The
laws are profile-free statements about H's eigenspaces; the locus enters only in
the pairing step. Off it, the Δ = 0 extremes of (1,3) leave the centre interval by
0.3307 on each end and the Δ = 1 extremes by 0.0051 on each end (gated from
below), and the equality of the two ends is EXACT, not approximate: (1,3) is a
self-folded block (bra index N/2), the ket/bra fold F: (a, b) ↦ (a, b̄) is an
exact involution that commutes with ad_H (from [H, X^N] = 0) and sends
F D F = −D − 2σ (entry-wise 0.0 in the gate, the dyadic profile making the float
sums exact), so EVERY compressed spectrum on such a block is symmetric about −σ
for any profile, and since s_min + s_max = N on a self-folded block, the two
overshoot ends coincide. The sibling statement is PROOF_CODIM1_BY_ADDITIVITY §7's
clause (b), the s-symmetric self-folded block; two differences keep the objects
apart. That clause speaks of the FULL block Liouvillian at uniform γ, where this
paragraph needs only the compressed D (F flips D exactly but not the whole L);
and the proof's fold leg dresses the bra complement with the bipartite staggering
𝒟, which is what scopes it to Δ = 0 (its clause (c)), while the F here is the
undressed one-sided X^N of MirrorWorld's Lattice, needing only [H, X^N] = 0,
hence any Δ and any profile. The first draft of this note called the equality
"measured structure with no explanation", which violated the house no-rounding
rule twice (it read an exact zero as "four digits" and missed the owned
mechanism). The genuinely two-sided
off-locus case is (2,4): its lower end never moves (the arithmetic endpoint), its
upper end escapes by 0.3307. The Δ = 1 overshoot is small because the uniform law
survives any profile: comp(D) on Ω_s is the centre −2γ̄s for any γ, so by
interlacing the extremes reach AT LEAST the centres, and only the excess beyond is
unprotected.

## What is ours and what was owned

Owned: the per-site split of D (Absorption Theorem 2), the size operator and the
aggregate compression (F122, PROOF_STRUCTURAL_CEILING), the |S|-resolved locus
reflection and its centre (PROOF_R90_FROZEN_DIVISOR, F153, the R₉₀ note's
`two_sources_gate.py`), the kernel-compressed double occupancy with the chiral
involution inside at rung 1 (PROOF_FROZEN_BAND_SO4 §6), the cell-rate window and the
self-folded fold symmetry (PROOF_CODIM1_BY_ADDITIVITY §6 and §7 (b)), the
single-particle reflection parities (review/OBC_SINE_BASIS_FINDINGS). Ours: the parity theorem on scalar eigenspaces and
the census that shows it covers everything below the N = 8, Δ = 0 frontier; the
mixed-space reflection law there (with the generic-R-odd contrast); the uniform law
at Δ = 1 and its failure at Δ = 0; the identity Π D Π = −2γ̄ Π N_XY Π on the locus
with containment-on-centres and attainment; the parity characterization of the Δ = 1
breaks (mixed necessary, not sufficient); the explanation of
collapse/saturation/small-resonance sizes as one mechanism; and the rational
interior.

## What is still open

Why the mixed-space reflection law holds at Δ = 0, N = 8: the cross matrix elements
between the parity halves vanish for a reason the one-body selection rules alone do
not supply (candidate collisions e_i − e_j = e_m − e_n with mixed parity products
need to be absent or to cancel; F129's collision law is the natural tool and is not
applied here). Why the uniform law holds at Δ = 1 (an SU(2) ladder argument in
F144's genre is the obvious candidate). What separates the Δ = 1, N = 6 mixed
spaces that break from the two on (1,3) that hold, and why the breaking spaces'
compressed spectra stay INSIDE the centre interval on the locus (measured with
slack, underived). And whether mixed collisions and the law survive at N ≥ 10,
where F129-type cosine collisions grow (the census stops at N = 8).

## The gate

`simulations/endpoint_density_gate.py`, five sections. (A) the parity census and
theorem layer (30 case rows over N = 4..8, one vacuous), scalar-side residuals
gated per eigenspace against res ≤ 64·eps·‖A‖/gap in the spectral norm (a
max-entry residual would depend on LAPACK's arbitrary basis inside a degenerate
space and is not reproducible; the norms are), plus the R-involution check and the
mixed counts (0 through N = 7 and at Δ = ½, 91 at N = 8, Δ = 0). (B) the
mixed-space law with the generic-R-odd contrast (0.119 vs floor) and the bond
control including the ten-profile random palindromic sweep, breaks gated from
below at 0.05. (C) the Δ = 1 parity characterization (breaks ⊂ mixed with ω, dim
and parity printed, ω = 0 scalar, N = 8 all scalar), the uniform-law census, and
the Δ = 0 uniform-law failure gated from below at 0.1. (D) the locus identity
(spectral norm), containment for EVERY eigenspace, attainment at both endpoints,
the zero-pure-vector count of the (1,3) Δ = ½ row gated exactly (its (2,4) sibling
saturates fully, [−6, −2], gated in the same section), block-level ranges, and
the dim-12 rational interior at ω = 0 (integer multiset compared exactly after
rounding, defect gated). (E) the fold identities on (1,3) compared == 0.0, both
off-locus overshoots with their difference gated at the fold floor, and the
(2,4) asymmetric pair (arithmetic lower end, escaping upper end). The eigenvalue
grouping tolerance sits in the difference spectrum's measured gap void (the
predecessor gate's gtol note). Measured-not-gated, complete: the predecessor's
0.4494 / 0.4467 / 0.5325 / 2.0000-at-(1,4) values (its gate and page), and the
exploratory per-eigenspace tables behind the prose (local scouts; everything
load-bearing is in the gate).
