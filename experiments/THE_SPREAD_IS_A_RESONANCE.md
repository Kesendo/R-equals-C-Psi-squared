# The spread is a resonance, and the Δ table had no three behaviours

*2026-08-15. [WHAT_THE_R90_LOCUS_BUYS](WHAT_THE_R90_LOCUS_BUYS.md) closed with one open
question it called its sharpest: on the non-pinned index-N/2 blocks, the strong-coupling
Re-spectrum saturates on the size-class centres at Δ = 0 and 1, flattens onto the trace
constant at Δ = 2, and does neither at Δ = 0.5, and nothing said what selects. This note
answers it, and the answer dissolves the question's premise: there is no selection among
three behaviours. On the index-N/2 blocks the question lives on, FLAT is the generic case,
exact and profile-free, and the spread is a RESONANCE: there it REQUIRES a collision of
the difference spectrum of the two Hamiltonian sectors (necessary on those blocks because
the X^N cancellation pins every collision-free mode, and not sufficient: on the R₉₀ locus
the pinned blocks collapse at Δ = 0 although collisions are present, which is F153's
involution content seen from this side; on blocks with NO index at N/2 the cancellation is
absent and a spread needs no collision at all, the fourth behaviour below). The four-point grid Δ ∈ {0, 0.5, 1, 2} had sampled the three resonant
points and one generic one. The projection rule itself is the repo's (F122's registered
mechanism); what is new is its scope, the cancellation that makes flat exact, the resonance
catalogue with an exact ℤ-certificate at Δ = 0.5, the derivation of F153's N clause from
that catalogue, and a fourth behaviour on the blocks nobody had measured.*

## What this is about

The chain's coupling has a tunable character, the anisotropy dial. A
predecessor note sampled four settings of that dial and saw what looked like
three different behaviours of the strong-coupling decay-rate spectrum: spread
out at two settings, collapsed onto a single flat line at a third, neither at
the fourth, with no rule saying which setting gets which. This note dissolves
the question: there is no selection rule, because there were never three
behaviours. The families of patterns in question are those where one of the
two spin arrangements the pattern connects has exactly half of its sites
excited, and on them flat is what a generic setting of the dial produces,
exactly and regardless of how the light, the outside dephasing that gives
every pattern its decay rate
([Dephasing Translated](../docs/quantum/DEPHASING_TRANSLATED.md)), falls
across the sites: the pattern's half-filled face is carried by a wave that
spreads its weight perfectly evenly across the sites, and that evenness wipes
the other face's leanings out of the average. A spread is the
exception and needs a resonance, two of the underlying energy ladders
producing coinciding differences, clocks ticking in step. The four settings
the predecessor happened to sample were the three structural resonant ones
and a single generic one, so the illusion of a rich selection rule was the
grid's, not the physics' (the scan establishes what is generic, it does not
exclude further isolated resonances; and one family carries such a
coincidence at every setting and is never in the generic case). The note names the source of each resonance,
certifies the subtlest one in exact integer arithmetic, and finds a genuinely
fourth behaviour on the families nobody had measured, where the half-filling
evenness is absent and a spread needs no coincidence at all.

## The system

System throughout: open chain, N sites, Pauli book H = J·Σ_bonds(XX + YY + Δ·ZZ), per-site
Z-dephasing γ_l, σ = Σγ_l, γ̄ = σ/N. A **coherence block** (p, q) is the invariant set of
cells |a⟩⟨b| with ket popcount p and bra popcount q; on it L(J) = J·A + D with A = −i·ad_H
at unit coupling (anti-Hermitian, carrying the hops and the Δ·ZZ frequency) and D the real
diagonal dissipator, D_cell = −2·Σ_{l: a_l ≠ b_l} γ_l. The **locus profile** here is
γ = (1, 2, 3, 13, 14, 15)/16 (anti-palindromic, γ̄ = ½), the **off-locus profile**
γ_l = (l+1)²/16; both are the profiles of this note's gate and of `two_sources_gate.py`
(F153's C# theory uses a different off-locus profile, `OffLocusProfile`, and gives
different values, as its entry says). A cell's **size class** |S| is popcount(a⊕b), the
number of sites where ket and bra disagree.

## What the repo already held, store by store

The sweep ran before the derivation (a dispatched scout over the named stores, its findings
verified against the files), and its headline is that the RULE was already registered and
the QUESTION already filed open, so the work below is a scope extension plus a mechanism,
not a new rule.

`docs/ANALYTICAL_FORMULAS.md` returned **F122**, whose registered *Mechanism* is this
note's projection rule: "high-Q degenerate PT: the decay rates are the eigenvalues
of N_XY (diagonal in the coherence basis, entry hamming(a,b)) block-diagonalized by the
ad_H eigenspaces"; the XY (Δ = 0), uniform-γ, slowest-mode scope is
PROOF_STRUCTURAL_CEILING's, the entry naming no Hamiltonian itself. It returned **F153**
with the Δ saturation readings and their N clause as measured facts without a mechanism,
and **F152**, whose (0,1)-block identity is an all-J statement available because on that
block "vacuum and one excitation disagree in exactly one bit"; F152's own fence says the
Laplacian form needs |Δ| = 1, a warning that the ZZ term acts inside the sector. It
returned **F146**, which already owns the Δ = 0 collisions in closed form: the chiral-only
rungs and their failure "on the cosine resonances of F129's kind" are exactly the
free-fermion difference-spectrum degeneracies of the catalogue below (F129 itself is NOT
cited for the Δ = 0.5 doublet: its level map presupposes Slater additivity, which the ZZ
term destroys at Δ ≠ 0).

`docs/proofs/` returned **PROOF_STRUCTURAL_CEILING §1**, which states the rule with a
numerical gate ("for each ad_H-eigenspace Ω the rates are 2γ · the eigenvalues of
P_Ω N_XY P_Ω", reproducing the full-L g2 at Q = 10³ to O(1/Q)), again XY, uniform γ,
slowest mode, and which already ties the ceiling to the DEGENERACY of the single-particle
spectrum. It also returned the negative result not to re-walk: the universal 4/(m+1) law in
the degeneracy m fails on the ring. **PROOF_FROZEN_BAND_SO4** owns the same first-order
step from the other side (its §5, the large-J reduction), with §4's warning that at finite
coupling the frozen condition's 𝒦 does not commute with ad_h, so the two conditions meet
only on ker(ad_h). **PROOF_CODIM1_BY_ADDITIVITY §6** owns the rate window every union below
can be read against (written there at uniform γ; the profile form is F153's γ-fence
closure, and the gate's pins are sharper than the window, which is not re-checked here). **MIRROR_SYMMETRY_PROOF** (line 928) states [X^⊗N, H_XXZ] = 0 outright, so the
cancellation's one ingredient is docs/proofs' as well as MirrorWorld's. The γ-fence
closure of F153 (Herm(L_block) = −2·diag(rate), exact) is the from-below constraint on
every prediction here: the union's trace and window are fixed before any perturbation
theory runs.

`experiments/` returned [WHAT_THE_R90_LOCUS_BUYS](WHAT_THE_R90_LOCUS_BUYS.md), the
immediate predecessor: the Δ table this note explains, the three named behaviours, the
warning that Δ = 0.5 is not a transient, the J-ladder trap, and the question as it stood
at commit b7490d0 ("Why Δ selects among the three is not answered here"; the predecessor's
paragraph now points here instead, repaired in the same change as this note). The committed implementation of
the rule is `simulations/topology_ceiling_rep_derivation.py` (Stage 0 an actual gate),
XY only. Three more experiments bear directly and were review-round finds, not first-pass
ones: [ETA_CEILING_REDUCTION](ETA_CEILING_REDUCTION.md) already computes the ω = 0
equal-energy classes exactly (Φ_2M divisibility), i.e. the difference-spectrum collision
counted and certified for its own rung question; [XY_FROZEN_BAND](XY_FROZEN_BAND.md)
carries an unrelated-looking N = 4 arithmetic ("N = 4 is a different reading … the frozen
root and its fold partner are the same number exactly at N = 4"), a DIFFERENT mechanism
for a DIFFERENT N = 4 specialness, named here so nobody conflates it with the N clause
below, which is a degeneracy count and not that coincidence; and
[THE_EXCEPTIONAL_COUPLINGS](THE_EXCEPTIONAL_COUPLINGS.md) is the standing answer to the
first-order-exactness scope question in the open list ("the finite set is not empty").

The OpenArcs registry returned three arcs, and the first sweep named only one of them (a
review round supplied the other two, which is a Stage-0 miss worth recording).
**`xxz_axis_handover`** owns the Δ axis: its object is the slowest mode's Δ* handover,
living in the diagonal population block (p, p) at p = ⌈N/2⌉, and its γ → 0 secular
reduction (a rate matrix over the half-filling eigenbasis whose diagonal is
−4·Σ_k Var(n_k)) is this note's density overlap at m = n, a cousin statement on the other
perturbative axis. **`topology_band_edge`** is the arc PROOF_STRUCTURAL_CEILING resolved,
i.e. the registry's own carrier of the rule, with the 4/(m+1) null beside it.
**`f_registry_meets_the_typed_layer`** holds, dated the same day as this note, the Δ = 0
strong-coupling re-flattening of the pinned blocks and F153's N clause; the note's
catalogue must not be read against that entry's Δ = 0 sentence, which is the pinned
on-locus case (collision present, spread closed), not this note's non-pinned one.
`docs/GLOSSARY.md` returned nothing for the strong-coupling projection.
`docs/CAUGHT_ERRORS.md` returned the two Δ traps this note must respect: the W1Dispersion
over-scoping (the clean band is Δ = 1 only; at Δ = 0 and 0.5 the sector spectrum is
completely different) and the Δ-inversion entry (Δ ≠ 0 flips the free-fermion diabolic to
defective), plus a recorded history of fabricated Δ = 0.5 digits, which is why every
GATED number below comes from the committed gate and the gate section closes with the
list of what is measured instead, and the recorded uniform-γ ⟹ profile
scope-extension failure (:583), which is this note's central move and the reason the
cancellation is gated on both profiles rather than argued once. `fw.Confirmations` returned nothing; no flight
touches this. MirrorWorld returned [H, X^N] = 0 including ZZ (the Lattice's standing
exact fact), which is the cancellation's one ingredient, and no J → ∞ rule of its own.

## The rule, scope-extended

F122's mechanism, run beyond its registration: at J → ∞ the Re-spectrum of the block is

    ⋃_ω spec(Π_ω D Π_ω),    Π_ω the eigenprojectors of A = −i·ad_H,

degenerate perturbation in 1/J. Gated here on XXZ (Δ ∈ {0, 0.5, 1, 2}), a γ PROFILE, the
FULL spectrum, and the four distinct (p,3) blocks at N = 6 (the mirrors (4,3), (5,3),
(6,3) repeat value for value, the predecessor's record), against an error model with two terms, and
the first is a decade STRONGER than first-order perturbation usually buys: the truncation
of the Re-spectrum is **O(1/J²), not O(1/J)**, because A is anti-Hermitian: with ω_k the
eigenvalues OF A (purely imaginary) and D Hermitian in A's orthonormal eigenbasis, the 1/J
correction Σ_m |D_mk|²/(J(ω_k − ω_m)) has a real numerator over imaginary denominators and
is purely imaginary; on a degenerate eigenspace the same holds for the effective
second-order operator Σ_m Π D|m⟩⟨m|D Π/(J(ω − ω_m)), Hermitian-PSD numerators over
imaginary denominators making it anti-Hermitian. The real parts get no 1/J term at all,
degenerate or not. Measured: the deviation
shrinks by ×100 per decade of J wherever the truncation dominates (gated as a ratio in
[30, 300] from J = 10³ to 10⁴). The second term is the dense REFERENCE degrading, not the
rule: the eigensolver's backward error scales as eps·J·‖A‖ and overtakes the truncation
between 10⁴ and 10⁶ (at 10⁶ the reference is only good to ~10⁻⁸ and the deviations sit
exactly there; measured prefactor 1 to 11, gated at 32). Below the rule's validity scale
J* ~ ‖D‖/(min difference-gap) the comparison is meaningless and the gate prints rather
than gates (the Δ = 2 rows at J = 10³; the predecessor's J-ladder trap is this scale
measured from the other side). This 1/J² fact retroactively explains a number the repo
already owned: PROOF_STRUCTURAL_CEILING's own gate agrees "≈ 10⁻⁶ at Q = 10³", which is
1/Q², and that proof's prose labels it O(1/Q); the label is repaired alongside this note. The
γ̄-proportionality the predecessor's table carries (every value scaled by exactly 3/2 at
γ̄ = 0.75) is this rule's free corollary, gated exactly: A depends only on H, D is linear
in γ, and pred(1.5γ) = 1.5·pred(γ) to 10⁻¹².

## Flat is generic, and the three behaviours were a grid artifact

On a 24-point Δ grid the predicted spread of (1,3) and (2,3) is zero except at exactly
three points:

    Δ:        0      0.1 … 0.4    0.5      0.6 … 0.9    1.0     1.1 … 4.0
    (1,3):    2.0000    0         0.4494      0         2.0000     0
    (2,3):    4.0000    0         0.5663      0         4.0000     0

The predecessor's four measured Δ were the three resonant points plus one generic one, so
its trichotomy FLAT@trace / sat@centres / NEITHER is not three regimes of a selection; it
is one generic behaviour plus a resonance catalogue, and NEITHER is simply a small
resonance where the other two are large ones. Two scope clauses. The grid establishes
genericity, not exhaustiveness: a collision condition is one equation in Δ, so further
isolated resonant Δ off this grid are not excluded, only the three structural ones at
Δ ≥ 0 are named. And the catalogue is Δ ↦ −Δ symmetric: the sublattice rotation
U = Π_{l odd} Z_l gives U·H(Δ)·U = −H(−Δ), so the sector spectra negate, the collision
structure transfers, and a scan over [−2, 4] finds exactly the mirror set
{−1, −0.5, 0, 0.5, 1} (the gate spot-checks −0.25 flat, −0.5 at 0.4494, −1 at 2.0000).

## Why flat is exact, and why it needs no locus

On a block (p, N/2) the eigenvectors of A are the product cells ψ_i ⊗ conj(φ_j), ψ from
sector p and φ from the half-filling sector, with eigenvalue −i(E_i − F_j). The chain
Hamiltonian commutes with the global spin flip X^N exactly (entry-wise zero; the repo's
standing fact), and X^N maps the half-filling sector to itself, so every NON-degenerate φ
satisfies |φ(a)|² = |φ(ā)|² and its site density is exactly n_l = ½. In the pair average

    ⟨rate⟩ = Σ_l γ_l·(m_l + n_l − 2·m_l·n_l),

n = ½ cancels the other sector's density COMPLETELY: ⟨rate⟩ = σ/2, so the mode sits at
exactly −σ, for ANY γ profile, with no locus anywhere in the argument. Gated on the blocks
(0,3), (1,3), (2,3), on both the locus and the off-locus profile, at three generic Δ:
spread and distance from −σ both below 10⁻⁶. (−σ is the trace constant
−2γ̄(p + q − 2pq/N) at q = N/2, so this agrees with the predecessor's flat value and
extends it off the locus. "Rate" here is the dissipator's argument Σγ_l over disagreeing
sites; the cell's Re λ is −2·rate, hence the factor between σ/2 and σ.) The one hypothesis
is SIMPLICITY OF THE DIFFERENCE EIGENVALUE E_i − F_j, not of φ itself, and the asymmetry
is worth stating: a simple φ_j makes Π D Π equal −σ·I on the entire family
{ψ_i ⊗ conj(φ_j)}, so degeneracy on the ψ side alone is harmless; conversely at Δ = 1
every φ is simple and the spread still opens, because the collisions there are
cross-sector. Where the DIFFERENCE is simple, the mode is exact; where it collides, the
resonances enter.

## The resonance catalogue

On a block with an index at N/2, spread requires a COLLISION of the difference spectrum
{E_i − F_j}: the cancellation pins every collision-free mode at −σ, so only a colliding
eigenspace CAN move (on blocks without an index at N/2 no collision is needed, the fourth
behaviour below; the (1,2) block at Δ = 2 spreads 0.9982 with a difference spectrum that
is entirely simple). A collision makes
mixing possible, not compulsory: on the R₉₀ locus at Δ = 0 the pinned (0,3) collapses
although the same collisions are present (Π D Π stays scalar on every colliding eigenspace
there), while off the locus the same block saturates at 0.5325. Necessity is the
mechanism's; sufficiency is the profile's and the block's. The measured sources:

**Δ = 0, free fermions.** The XX sector spectra are sums of single-particle levels and the
difference spectrum collides massively (eigenspace dimensions up to 12 on (1,3) at N = 6;
measured, not gated). Full saturation on the size-class-centre interval on the non-pinned
locus blocks. F146 already owns this collision structure in closed form, its chiral-only
rungs failing exactly "on the cosine resonances of F129's kind".

**Δ = 1, SU(2), and the collision is CROSS-sector.** The half-filling sector itself is
non-degenerate at N = 6, which refutes the first guess that within-sector degeneracy does
the work. What collides is the sectors against each other: every sector-1 state belongs to
a total-spin multiplet with S ≥ |S_z| = 2, whose S_z = 0 member lies in sector 3 at the
same energy, and the same argument runs from sector 2 (|S_z| = 1). The gate pins both
legs: all 6 sector-1 levels and all 15 sector-2 levels reappear exactly in sector 3 at
Δ = 1, none of either at Δ = 2. Shared levels make shared differences, and the resonance
is again full.

**Δ = 0.5, one exact double level at N = 6, and a different face at N = 8.** The N = 6
half-filling sector carries exactly one degenerate pair, at E = 5J/2, and the degeneracy is
EXACT: the characteristic polynomial of 2H (an integer matrix, built in exact ints)
satisfies p(5) = 0, p′(5) = 0, p″(5) ≠ 0 in exact integer arithmetic (Faddeev-LeVerrier
over ℤ in the gate). Every sector-1 level pairs with that doublet into a rank-2 difference
eigenspace: six of them on (1,3) and fifteen on (2,3) (one per sector-2 level), and those
small resonances are the whole NEITHER row, 0.4494 and 0.5663. At N = 8 the same Δ wears
the OTHER mechanism: the eight half-filling doublets are exactly the eight sector-1 levels
appearing in sector 4 (shared-level census 0 / 0 / 8 at N = 4 / 6 / 8, gated), so there the
Δ = 0.5 resonance is cross-sector like Δ = 1's, and the clean split "Δ = 1 is cross-sector,
Δ = 0.5 is a within-sector doublet" is an N = 6 reading, not a law. Where the N = 6 double
level comes from as algebra is left open below.

**(3,3), the standing collision.** On a block with p = q the difference spectrum contains
E_i − F_i = 0 for every i, an eigenspace of dimension 22 / 20 / 20 at Δ = 0.5 / 1 / 2
(pinned in the gate together with the commutant identity Σm_i²; the bare bound ≥ 20 holds
for ANY Hermitian H and is deliberately not the check). That makes the predecessor's
exempt block mechanical: (3,3) is never in the generic case, at any Δ.

## F153's N clause, derived

F153's Δ = 0.5 saturation reading carries an N clause (collapse at N = 4, saturation at
N = 6 and 8, measured on the (0, N/2) family on the off-locus profile) that was measured
and not explained. On that family the derivation is sharper than anywhere else, because
p = 0 makes the ket sector one-dimensional: the difference spectrum is the half-filling
spectrum negated and shifted (E_vac − F_j), so its collisions are exactly the half-filling
degeneracies, and the Δ = 0.5 double-level count

    N = 4: 0      N = 6: 1      N = 8: 8

(gated) is exactly collapse / saturation at 0.6456 / saturation at 1.2903. On this
off-locus profile at N = 4, 6, 8 the count matches the clause binarily: zero doublets on
the row that collapses, doublets on the rows that saturate. What the count does NOT give
is the magnitude, and an N = 8 pair is the counterexample kept on the page, measured on
the LOCUS profile (1, 2, 3, 4, 12, 13, 14, 15)/16, γ̄ = ½, and not in the gate: (1,4) and
(2,4) see the same eight doublets, and (1,4) saturates on its size-class centres (spread
2.0000, carried by one dim-16 eigenspace at ω = 0) while (2,4) reads 0.4467 (rank-2
spaces only). Which eigenspaces collide,
and how large they are, sets the value; the count sets only whether anything happens at
all. The rule also reproduces the family's Δ = 0 row, 0.2236 / 0.5325 / 1.0756 at
N = 4 / 6 / 8 (gated; F153 quotes the first two of this row, five of the six across both
rows, and 1.0756 matches the dense N = 8 measurement taken 2026-08-15).

## The fourth behaviour, predicted and then measured

The cancellation needs one leg in the half-filling sector. On a block with BOTH indices off
N/2 it is absent, and the rule predicts a fourth behaviour the predecessor's table could
not contain: a stable band whose width is the spread of the density overlaps
Σγ_l(m + n − 2mn) over eigenvector pairs, whose MEAN is the trace constant (the trace is
preserved; the range midpoint is not it) and which is pinned neither to it nor to the
size-class centres. Measured on (1,2) at
N = 6, Δ = 2, locus profile: spread 0.9982, range [−2.8510, −1.8528] about the trace
constant −7/3, identical at J = 10⁴ and 10⁶, predicted by the rule to 10⁻⁸. No committed
measurement had touched these blocks at strong coupling.

## What is ours and what was owned

Owned: the projection rule (F122, PROOF_STRUCTURAL_CEILING §1, gated at Δ = 0), the
[H, X^N] = 0 fact (MIRROR_SYMMETRY_PROOF and MirrorWorld's Lattice), the rate window
(PROOF_CODIM1 §6), the Δ = 0 collision structure (F146), every number in the
predecessor's Δ table. Ours: the scope extension of the rule to XXZ + profile + full
spectrum + all blocks, with the 1/J² truncation identified, derived and gated as a
scaling; the X^N density cancellation and its consequence that flat is exact,
profile-free and generic; the identification of the spread as a difference-spectrum
resonance (necessary, with the sufficiency boundary named) and the catalogue above,
including the N = 6/N = 8 mechanism split at Δ = 0.5; the exact ℤ-certificate of the
Δ = 0.5 double level; the derivation of F153's N clause as a binary resonance count on
the (0, N/2) family; and the fourth behaviour with its numbers. The name "high-Q degenerate
PT" is reused from F122 deliberately; "structural ceiling" is NOT used here because that
name is the slowest-mode consequence, and "coalescence" is avoided as EP-loaded.

## What is still open

Why the full resonances saturate on EXACTLY the size-class-centre interval is answered
in [THE_ENDPOINTS_ARE_A_DENSITY_LAW](THE_ENDPOINTS_ARE_A_DENSITY_LAW.md): on every
degenerate eigenspace where the chain reflection acts as a scalar, the compressed site
density is reflection-symmetric by symmetry alone, and the parity census shows that
covers every colliding eigenspace of its census blocks up to N = 7; on the locus
Π D Π = −2γ̄·Π N_XY Π follows, and the centres are Rayleigh bounds attained by the
pure-size-class vectors. Open there in turn: the parity-mixed eigenspaces that first
appear at N = 8, Δ = 0, where the same vanishing is measured but not forced. The origin of
the exact 5/2 double level at Δ = ½ (an exact rational double root in an otherwise
high-degree spectrum; whether it is the boundary of a family in N, and whether it connects
to the XXZ combinatorial point, is unasked here). The scope of the 1/J² truncation: the
gate shows the scaling, not a bound, and the repo's standing answer to "holds at all large
J ⟹ holds at generic J" is THE_EXCEPTIONAL_COUPLINGS' finite nonempty exceptional set, so
isolated couplings where the limit misbehaves are expected, not excluded. And the fourth
behaviour's band edges as a closed form.

## The gate

`simulations/high_q_selection_gate.py`, writing
`simulations/results/high_q_selection_gate.txt`; the last line is the VERDICT line.
Five sections mirror this note: (A) the rule with the two-term error model (the 1/J²
ratio gated in [30, 300] per decade, the reference floor at 32·eps·J·‖A‖, the below-J*
rows printed and labeled rather than gated); (B) the 24-point grid with the three
resonant values pinned, the Δ ↦ −Δ mirror spot-checked, and the γ̄-proportionality
compared exactly; (C) [H, X^N] compared == 0.0, the half-filling densities at Δ = 2 and
0.5 (the doublet exempted, the simple levels at ½ to 10⁻¹⁰), and flat-at-−σ on three
blocks and both profiles; (D) the catalogue: SU(2) shared-level counts, the Δ = 0.5
multiplicity HISTOGRAMS (so a triple cannot masquerade as a double), the shared-level
census 0/0/8, the exact ℤ charpoly certificate on an integer matrix built in exact ints
and compared == to the float route, and the (3,3) commutant dimension as the exact
equality Σm_i² (the bound ≥ 20 holds for any H and is deliberately not the check);
(E) the (0, N/2) family at Δ = 0 and 0.5 across N = 4, 6, 8, and the fourth behaviour.
One scope sentence the gate also carries: section (A) certifies the rule against a dense
reference built by the same block builder, so it is self-consistency of the projection,
and the physics is pinned by the hard-coded constants of (B) to (E). The
measured-not-gated list, complete: the eigenspace-dimension profiles, the (1,4)/(2,4)
magnitude contrast, the 0.01-step resonance scan over Δ ∈ [−2, 4] (the gate carries its
three spot checks), and the dense N = 8 Δ = 0 reading behind 1.0756 (the gate pins the
prediction; the dense reference at N = 8 was a one-off measurement).
