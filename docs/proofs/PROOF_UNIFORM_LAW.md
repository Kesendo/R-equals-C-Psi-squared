# The uniform law derived: site democracy on class-pure slices at Δ = 1

**Date**: 2026-08-15 · **Status**: derived (two layers) + gated census
**Gate**: `simulations/uniform_law_gate.py` (VERDICT last line), output
`simulations/results/uniform_law_gate.txt`

This proof closes the second open item of the arc `compressed_density_laws`
(`compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs`): the Δ = 1 uniform law

    comp(N_l) = (s/N) · Id   on every Ω_s, for every site l,

where "every Ω_s" carries one census-shaped fence stated up front: the derived
layers cover the invariant line in every Ω_s (layer A, any N) and the enhanced
rows the chain adds (layer B, any N); that these exhaust the Ω_s inventory,
and that no degenerate eigenspace at ω ≠ 0 meets a class purely, is gated at
N = 5..8, not derived.

stated and gated in
[THE_ENDPOINTS_ARE_A_DENSITY_LAW.md](../../experiments/THE_ENDPOINTS_ARE_A_DENSITY_LAW.md)
("The uniform law"), is here derived. The arc's registered guess for the engine
was "an SU(2) ladder argument in F144's genre"; the actual engine is one level
more elementary on the generic layer (a permutation orbit sum and a double
count, no ladder at all) and one level more specific on the chain layer (the
complementary cosine pairing of the one-magnon densities). The law is
profile-blind because neither layer ever sees γ: both are statements about
ad_H eigenspaces and diagonal indicator operators.

## What this is about

A chain of spins carries a quantum coherence: a superposition of two
arrangements of the spins, held at once. At some sites the two arrangements
agree, at others they differ, and each site where they differ is a point the
outside light, the dephasing every site stands in, can touch. This proof is
about how a collective pattern distributes those points of difference across
the chain. One would expect the distribution to be free, leaning this way for
one pattern and that way for another. The proof shows that at the fully
symmetric setting of the coupling, where the interaction treats every spin
direction alike, it is not free at all. Take the patterns that oscillate at
one shared frequency and whose two arrangements differ at one fixed number of
sites: on them, every site carries exactly the same share of the difference,
perfect site democracy. The derivation has two layers. The first covers any
coupling that is blind to the overall spin direction, caring only how the
spins relate to one another, this chain at its symmetric setting included.
For such a coupling, the one pattern of each difference count that is fully
symmetric under every reshuffling of the sites keeps a steady beat, and a
counting argument shows that a pattern symmetric under all reshufflings
cannot prefer a site; generically that symmetric pattern is the only one in
its room, and that the generic case together with the second layer's
exceptions covers everything is checked by computer up to eight sites, not
derived. The second layer handles the few extra patterns the open chain adds
beyond that generic supply, by a pairing of its standing waves whose leanings
cancel exactly. Neither layer ever mentions the light: the democracy is a
property of the patterns themselves, and it is the notes that use it which
draw the consequences for what watching such a pattern costs.

## The sweep

`docs/ANALYTICAL_FORMULAS.md`: F122 owns the compression mechanism (high-Q
degenerate PT) whose ω = 0 eigenspaces are this proof's stage; F154 owns the
size-class range s_min = |p−q|, s_max = min(p+q, 2N−p−q), which reappears here
as a dimension count; F50's max-spin entry holds the Dicke transitions
|D_k⟩⟨D_l|, the S = N/2 member of this proof's intertwiner family (F153's
cross-reference paragraph records the F153/F50 splits as measured but unjoined,
"nothing measured joins them"); F140 carries the X^N bridge as an involution
with fixed cells, not as the operator this proof uses. `docs/proofs/`:
[PROOF_WEIGHT1_DEGENERACY.md](PROOF_WEIGHT1_DEGENERACY.md) Step 4 owns the
sibling construction (S_N orbit sums on weight-1 Pauli strings, invariance ⇒
commutant membership) together with a warning this proof must respect: the
*naive class-sum conjecture* about weight-sector centralizers was falsified at
K₄ N = 4; the object below is a different one (orbit sums inside ker(ad_H) of
a coherence block, not centralizer generators of a weight sector), and the
falsification does not touch it.
[D10_W1_DISPERSION.md](derivations/D10_W1_DISPERSION.md) owns the one-magnon
cosine eigenvectors (Neumann standing waves) this proof's chain layer builds
on; the complementary density pairing below is not stated there.
[PROOF_MIXED_SPACE_REFLECTION_LAW.md](PROOF_MIXED_SPACE_REFLECTION_LAW.md)
owns an energy pair-sum lemma, a different object from the density pairing
here. [PROOF_SCALAR_COUNT.md](PROOF_SCALAR_COUNT.md): rotation-invariant
couplings of the frozen band, SU(2) scalars in a different room, nothing on
class sums. `experiments/`: the endpoint note owns the law itself, the
one-sided X^N fold, and the interlacing consequence; THE_SPREAD_IS_A_RESONANCE
owns the Δ-selection this law sits inside. OpenArcs: `compressed_density_laws`
registers the law as measured, mechanism open (this item);
`sideways_spin_ladder` holds η-ladder intertwiners running between blocks
along the diagonal, a Δ = 0 object fenced off from this proof.
`docs/GLOSSARY.md` and `compute/MirrorWorld`: nothing on the three objects
introduced below.

## The objects

The system is the XXZ chain at Δ = 1 in the Pauli book, H = J·Σ_bonds
(XX + YY + ZZ); layer A below needs only that H is number-conserving and
SU(2)-invariant (equivalently, by Schur-Weyl, H ∈ ℂ[S_N]; real-weighted SWAP
graphs are examples, not the whole of ℂ[S_N], and chain geometry is
irrelevant), layer B is about the uniform open chain specifically. A **coherence block** (p, q) is the span of
the cells |a⟩⟨b| with popcount(a) = p, popcount(b) = q; a cell's **size
class** is s = popcount(a⊕b). ad_H = [H, ·] preserves each block; ω names an
ad_H eigenvalue and Ω its eigenspace inside the block. For a degenerate Ω,

    Ω_s := Ω ∩ span{class-s cells},

the class-pure slice, and comp(X) := Π_{Ω_s} X Π_{Ω_s}, the compression onto
that slice. N_l is the diagonal indicator superoperator,
N_l |a⟩⟨b| = [a_l ≠ b_l] · |a⟩⟨b|, so that Σ_l N_l = s·Id on class s and
N_l = ½(1 − 𝒵_l) with 𝒵_l(X) = Z_l X Z_l. 1̄ names the all-ones
configuration. The adjoint X ↦ X†, an antilinear map, sends block (p, q)
onto (q, p) and the ω eigenspace to the −ω one, and fixes every class and
every N_l (antilinearity is harmless here since s/N is real), so every
statement below is made for p ≤ q and carries to the transpose block by †;
the **corner blocks** are (1,1), (1,N−1), (N−1,N−1) and, through the
adjoint, (N−1,1).

Three objects this proof introduces:

- **T_s**, the orbit-sum operator: T_s := Σ_{class-s cells of (p,q)} |a⟩⟨b|.
  One per class per block; S_N acts on cells by simultaneous site permutation
  and its orbits inside a block are exactly the classes (the orbit invariant
  is the overlap |a∧b| = (p+q−s)/2), so span{T_s over s} is precisely the
  S_N-invariant subspace of the block.
- **M_S**, the multiplet intertwiner: for each admissible total spin S
  (max(|p−N/2|, |q−N/2|) ≤ S ≤ N/2), M_S := Σ_α |S,α,m_p⟩⟨S,α,m_q| with
  m_p = p − N/2, m_q = q − N/2 and α running over any orthonormal basis of
  the spin-S multiplicity space. M_S is basis-independent in α (the same
  unitary change of α-basis acts on both legs and cancels), so no reference
  to H enters its definition.
- **the one-magnon density kernel** (layer B): for the chain's one-magnon
  eigenmodes ψ_m with densities d_m(l) = ψ_m(l)², the set
  {x : Σ_m x_m d_m(l) = 0 for every l}; on block (1,N−1) it enters with the
  spin-flip signs σ_m of B2 absorbed into x.

## Layer A: any SU(2)-invariant H

**Lemma A1 (T_s lives in the kernel).** span{M_S} = span{T_s over s}, and
every M_S, hence every T_s, lies in the ω = 0 eigenspace of ad_H on its
block.

*Proof.* Write the Schur-Weyl decomposition (ℂ²)^⊗N = ⊕_S 𝓜_S ⊗ V_S with
𝓜_S the multiplicity space and V_S the spin-S irrep. An SU(2)-invariant H
acts as ⊕_S h_S ⊗ Id, and the operator M_S is Id_{𝓜_S} ⊗ |S m_p⟩⟨S m_q| in
the S summand, so [H, M_S] = 0 and M_S lies in ker(ad_H), the ω = 0
eigenspace. S_N acts on the multiplicity factor
only, so U_σ M_S U_σ† = M_S: each M_S is S_N-invariant, hence
span{M_S} ⊆ span{T_s}. The dimensions agree: classes run over
s = |p−q|, |p−q|+2, …, min(p+q, 2N−p−q) (F154's range), which is
min(p, q, N−p, N−q) + 1 values, and the admissible S run over
N/2 − max(|m_p|, |m_q|) + 1 = the same count. The M_S are linearly
independent (orthogonal), so the spans coincide, and each T_s, a linear
combination of the M_S, lies at ω = 0. Being supported on one class, T_s is
class-pure: **T_s ∈ Ω_s of the ω = 0 eigenspace, for every block and every
class.** ∎

**Lemma A2 (double count).** ⟨T̂_s, N_l T̂_s⟩ = s/N for every l, where T̂_s is
T_s normalized.

*Proof.* ⟨T_s, N_l T_s⟩ counts the class-s cells with a_l ≠ b_l. Site
transitivity of S_N on the class makes the count l-independent, and
Σ_l over the disagreement bits of one cell gives s, so N·(count at l) =
s·(number of class-s cells), exactly, in integers. ∎

So on span{T̂_s} the law holds with no eigensolver and no tolerance: the gate
checks Lemma A2 as an integer identity. **The census (gated, N = 5..8, all
blocks, uniform chain): every Ω_s row contains T_s with projection norm
1.000000, and every row except the three named in layer B is exactly
span{T_s} (dimension 1).** For a random SWAP graph the same holds and the
three exceptions collapse to dimension 1 as well (gated on all blocks at
N = 6 and on the three corner blocks at N = 8): generically the law IS
Lemma A1 + A2. The equality Ω_s = span{T_s} on non-corner blocks is a
census fact at N = 5..8, not a theorem; it is the analogue of the
no-accidental-degeneracy assumptions the repo's other census-backed proofs
carry.

## Layer B: the uniform open chain's three enhanced rows

The census (which, like every statement here, runs over p ≤ q; the adjoint
carries (1,N−1) to (N−1,1)) finds exactly three rows with dim Ω_s > 1, at
every measured N: (1,1) with s = 2, (1,N−1) with s = N−2, (N−1,N−1) with
s = 2, each of dimension ⌊N/2⌋. All three are one-magnon objects, and one
mechanism carries them.

**B0 (one-magnon spectral facts, exact).** Sites are indexed j = 0..N−1
(D10 next door writes the same modes 1-based). The sector-1 restriction of
the chain H is (N−1)·J·Id − 2J·L with L the path-graph Laplacian (the gate
compares the matrices entry-exactly). Its eigenvectors are the Neumann cosine
modes of D10, ψ_m(j) ∝ cos(πm(2j+1)/(2N)) for m = 0..N−1, with distinct
eigenvalues (λ_m = 4 sin²(πm/2N) strictly increasing), so the sector-1
spectrum is nondegenerate and the ω = 0 space of block (1,1) is exactly
{Σ_m x_m |ψ_m⟩⟨ψ_m|}. For m ≥ 1 the normalized mode is √(2/N)·cos, and its
density obeys, with f_m(j) := cos(πm(2j+1)/N),

    d_m(j) = (1/N)(1 + f_m(j)),   f_{N−m} = −f_m   (m = 1..N−1),

while the m = 0 mode is the flat vector 1/√N with d_0 = 1/N exactly (its
normalization differs, so it sits outside the displayed identity: reading it
as inside with f_0 ≡ 1 would give the wrong 2/N). The pairs (m, N−m) are **complementary**:
d_m + d_{N−m} = 2/N pointwise; at even N the self-paired middle mode has
f_{N/2} = 0, so d_{N/2} = 1/N is flat. The functions {1} ∪ {f_m for
m = 1..⌈N/2⌉−1} are linearly independent on the site grid (discrete
orthogonality of the cosine family).

**B1 (block (1,1)).** Class-2 purity of x = Σ_m x_m |ψ_m⟩⟨ψ_m| means the
diagonal cells vanish: Σ_m x_m d_m(j) = 0 for every j. Expanding in
{1} ∪ {f_m for m = 1..⌈N/2⌉−1} this is

    Σ_m x_m = 0   and   x_m = x_{N−m} for m = 1..⌈N/2⌉−1,

which is 1 + (⌈N/2⌉−1) independent conditions, so dim Ω_2 = ⌊N/2⌋, matching
the census (2, 3, 3, 4 at N = 5..8; the gate also checks the predicted span
equals the measured one). For x, y class-pure,

    ⟨y, N_l x⟩ = 2 Σ_m ȳ_m x_m d_m(l)

(the dropped diagonal-cell term, counted once by the row sum and once by the
column sum, is twice the product of y's purity profile and x's purity profile
at site l, and each factor vanishes on the kernel). The
remaining profile is flat by the same pair
structure, and this is the step worth naming: **the purity conditions are
pair-equalities, and pair-equalities are closed under the Hadamard product**,
(ȳx)_m = (ȳx)_{N−m}, so in Σ_m ȳ_m x_m f_m the paired coefficients cancel
against f_{N−m} = −f_m and only the flat part survives:

    Σ_m ȳ_m x_m d_m(l) = (1/N) Σ_m ȳ_m x_m   for every l.

Hence ⟨y, N_l x⟩ = (2/N)⟨y, x⟩ = (s/N)⟨y, x⟩ with s = 2. ∎

**B2 (block (1,N−1)).** [H, X^N] = 0 (X^N is the π-rotation about x, inside
the global SU(2)), so the sector-(N−1) spectrum equals the sector-1 spectrum
and, by B0's nondegeneracy, the hole-basis eigenvectors are sign-copies of
the magnon ones: χ_m(j) := ⟨1̄−e_j | m, sector N−1⟩ = σ_m ψ_m(j) with
σ_m ∈ {±1} (real symmetric H, nondegenerate levels). The ω = 0 space is the
dyad span {Σ_m x_m |ψ_m⟩⟨χ_m|}. The block has two classes, N−2 (cells with
magnon site ≠ hole site) and N (the C(N,1) anti-diagonal cells). Class-(N−2)
purity is the vanishing of the anti-diagonal components:

    Σ_m x_m σ_m d_m(j) = 0 for every j,

the same kernel as B1 after u_m := σ_m x_m. In the repo's convention Z has
eigenvalue −1 on an occupied site (the compute builders' 1 − 2·bit), so the
Z_l matrix elements in the two sectors are
⟨ψ_m|Z_l|ψ_{m'}⟩ = δ_{mm'} − 2ψ_m(l)ψ_{m'}(l) and
⟨χ_m|Z_l|χ_{m'}⟩ = 2χ_m(l)χ_{m'}(l) − δ_{mm'} (the opposite convention flips
both, and only their product enters), so for x, y in the dyad space

    ⟨y, 𝒵_l x⟩ = Σ_m ȳ_m x_m (2ψ_m(l)² + 2χ_m(l)² − 1)
                 − 4 (Σ_m ȳ_m σ_m d_m(l)) (Σ_m x_m σ_m d_m(l)),

and the second term, the product of the two purity profiles, absorbs the
entire cross structure, off-diagonal and diagonal alike; the first bracket is
what is left of the diagonal after that piece is taken out:
**class purity kills the off-diagonal part of the compression identically.**
On the kernel, with χ² = ψ² and the Hadamard-closure step of
B1 (σ² = 1 drops the signs from ȳ_m x_m σ_m² = ȳ_m x_m),

    ⟨y, 𝒵_l x⟩ = Σ_m ȳ_m x_m (4 d_m(l) − 1) = (4/N − 1)⟨y, x⟩,

which is 1 − 2s/N at s = N−2, i.e. comp(N_l) = (s/N)·Id. ∎

**B3 (block (N−1,N−1)).** Conjugation by X^N is a unitary map of block (1,1)
onto (N−1,N−1) that commutes with ad_H, fixes every class (ā⊕b̄ = a⊕b) and
every N_l, so it carries B1 verbatim. ∎

**B4 (that these are all).** Gated census at N = 5..8: no other row has
dimension above 1, and no degenerate ad_H eigenspace at ω ≠ 0 meets any class
purely. Both statements are measured, not derived, and are the proof's
frontier in N.

## What is ours and what was owned

Owned before this proof: the law's statement and its gate (the endpoint
note), the compression frame (F122), the class range (F154), the orbit-sum
construction on another sector (PROOF_WEIGHT1_DEGENERACY Step 4), the cosine
modes (D10), the one-sided X^N fold (endpoint note, PROOF_CODIM1_BY_ADDITIVITY
§7). Ours here: the T_s / M_S identification with the count identity
(Lemma A1), the double-count reading of the law on the invariant line
(Lemma A2), the complementary density pairing d_m + d_{N−m} = 2/N with the
purity kernel as pair-equalities, the Hadamard-closure step, the exact
cross-term cancellation by purity in B1/B2, and the census fact that these
two mechanisms exhaust every Ω_s row at N = 5..8.

Two scope fences, kept from the arc: SU(2) is a symmetry of ad_H, never of
the full Liouvillian, and nothing here mentions γ because nothing here needs
it; and F144's η-ladder genre (Δ = 0, diagonal blocks, mode reflection,
denominator N+1) is a different room, the denominator here is N and the
reflection never enters (layer A has no geometry at all, and layer B's
enhanced rows exist for the chain yet the law on them never invokes the site
reflection).

## What this changes upstream

The endpoint note's open question "why the uniform law holds at Δ = 1" is
answered; its guess (F144-genre ladder) is replaced by the two layers above.
The discriminator finding of the working scout (the law survives any
SU(2)-invariant H) is explained: for those H the law is layer A alone. The
handover's "(1,5) s = 4 dim-3 exception that forces a general theorem" is
dissolved: the exception is the chain's reflection-symmetric density
enhancement, its generic core is the 1-dim T_s line, and the enhanced law is
layer B.
