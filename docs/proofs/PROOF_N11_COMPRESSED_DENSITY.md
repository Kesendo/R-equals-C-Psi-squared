# N=11: the fixed-frequency mirror contrast speaks, while the `(1,1)` interval holds

*2026-09-23. A scoped sequel to [the mixed-space reflection law](PROOF_MIXED_SPACE_REFLECTION_LAW.md) and [F154](../ANALYTICAL_FORMULAS.md#f154). The exact physical-cell certificate is `simulations/n11_compressed_density_gate.py`; the typed reading is `CompressedDensityLocusClaim` with `CompressedDensityN11Witness` (`inspect --root compresseddensity`), which recomputes the same radicals in exact Q(√2, √3) arithmetic, all 55 N = 11 frequency rooms included. The measured readings beyond N = 11 are produced by `simulations/n11_compressed_density_census.py` into `simulations/results/n11_compressed_density_census.txt`, sections (A) to (E), each with its error model.*

## What this is about

Consider an eleven-site open XY spin chain under local Z-dephasing: each site
is watched at its own rate. A `(1,1)` coherence |ψ_a⟩⟨ψ_b| connects two
one-excitation standing waves. Two such coherences can turn at exactly the
same frequency but react with opposite signs to the chain's site mirror.
Below N = 11 the `(1,1)` block holds no such pair at all. Here the question
at the left end does connect them, while the mirrored question
at the right end answers with the opposite sign. A balanced but uneven
watching profile therefore retains a trace of *where* it watches: the old
shortcut that replaced it by uniform watching stops working. Yet the rate
of every complete fixed-frequency `(1,1)` compression stays inside the
familiar interval. That protection comes from the physical disagreement cost and the
positive weight on cells whose two faces coincide, rather than from the
mirror-cancellation shortcut. The interval's two ends sit in the
zero-frequency room, where every coherence is mirror-even and the shortcut
is a theorem.

## Abstract

For the uniform open XY chain at N = 11, two one-excitation coherences in one
Hamiltonian frequency space have opposite site-reflection parity. Their
projected left-site minus right-site disagreement has the exact matrix
element ⟨Y,C₀X⟩ = −√2/72. An explicit nonnegative mirror-balanced profile
therefore falsifies F154's conditional compressed-density identity when its
C_l = 0 premise is dropped. Separately, for any Hermitian, simple-spectrum,
reflection-symmetric one-excitation Hamiltonian and any nonnegative
mirror-balanced profile, the complete `(1,1)` frequency compression has the form
D_Ω = −4γ̄ I_Ω + 4P_Ω T P_Ω and spectrum in [−4γ̄,0]; this interval theorem is
the new part. Both ends of the interval are attained for every such Hamiltonian,
in the zero-frequency room, and there they belong to results the repo already
held: that room is a scalar-parity room where F154's identity is a theorem, its
compression is −4γ̄(I − G) with G the Gram of squared mode amplitudes, in closed
form on the uniform XY chain by [PROOF_R90_FROZEN_DIVISOR](PROOF_R90_FROZEN_DIVISOR.md)
Lemma 5 (re-derived as [F143](../ANALYTICAL_FORMULAS.md#f143)),
and for real h the lower end is [F140](../ANALYTICAL_FORMULAS.md#f140)'s frozen
root, an exact Liouvillian eigenvalue at every coupling.

Here mirror-balanced means γ_l+γ_{N−1−l}=2γ̄; C_l is the projected
disagreement at site l minus that at its mirror site. F154's conditional
identity replaces the uneven profile by its mean rate only when every
C_l vanishes on the frequency space.

## What the repo already holds

The sweep, by store, and what each returned.
`docs/ANALYTICAL_FORMULAS.md`: F154 is the entry this proof extends, with its
parity theorem (C_l = 0 on every scalar-parity eigenspace) and its pure-vector
attainment clause; [F122](../ANALYTICAL_FORMULAS.md#f122) owns the
strong-coupling compression; [F91](../ANALYTICAL_FORMULAS.md#f91) the locus;
F140 owns −4γ̄ as an exact eigenvalue of the `(1,1)` corner block at every J,
with multiplicity at least ⌊N/2⌋, for any real symmetric reflection-invariant h
and any real locus profile; F143 owns the zero-frequency reduced operator at
uniform rate on the XY chain, G = (1/M)(𝟏𝟏ᵀ + (I+R_χ)/2), M = N + 1, R_χ the
chiral mode flip k ↦ M − k, with spectrum {0 ×⌊N/2⌋, 1/M ×(⌈N/2⌉−1), 1 ×1}
and the chiral differences as kernel, and its own "Anticipated by" paragraph
names the earlier owner below; F144
owns the chiral transpose (a,b) ↦ (M−b, M−a); F145 names P_k − P_(M−k) the seed;
F146 records the N = 11, ℓ = 2 resonant rung (14 against 10).
`docs/proofs/`: [PROOF_MIXED_SPACE_REFLECTION_LAW](PROOF_MIXED_SPACE_REFLECTION_LAW.md)
names N = 11, 14, 20 as the first parity doors and the M = 12 coincidence
e₁ + e₉ = e₅ + e₆; [PROOF_R90_FROZEN_DIVISOR](PROOF_R90_FROZEN_DIVISOR.md)
§1-3 is F140's theorem, its Lemma 5 (2026-07-22, five days before F143)
derives the site-indexed law B Bᵀ = (1 − 1/M)𝟙𝟙ᵀ/N + (I + R)/(2M), B_ak = u_k(a)²,
for both open chains (XY at M = N + 1, Heisenberg at M = N) and with it the
spectrum {1, 1/M ×(⌈N/2⌉−1), 0 ×⌊N/2⌋} of the mode-indexed G = BᵀB; at
M = N + 1 the law is F143's G itself, and its Lemma 2 holds the recentred rate operator of the
corner block with defect 8γ̄·P_D, of which the decomposition D = −2(Γ⊗I + I⊗Γ) + 4T
below is the compressed sibling; [PROOF_FROZEN_BAND_SO4](PROOF_FROZEN_BAND_SO4.md)
§5-6 holds the large-J reduction G = P D̂ P and "the corner's frozen space is
spanned by the differences P_k − P_k̄"; [PROOF_CODIM1_BY_ADDITIVITY](PROOF_CODIM1_BY_ADDITIVITY.md)
§6 puts the N = 5 chiral pairs at the window bottom; [PROOF_SCALAR_COUNT](PROOF_SCALAR_COUNT.md)
§7 measures the rung-2 resonances at 6|M, 15|M and 21|M;
[PROOF_UNIFORM_LAW](PROOF_UNIFORM_LAW.md) holds the Δ = 1 second law.
`experiments/`: [THE_ENDPOINTS_ARE_A_DENSITY_LAW](../../experiments/THE_ENDPOINTS_ARE_A_DENSITY_LAW.md),
[THE_TWO_SPIN_ZEROS](../../experiments/THE_TWO_SPIN_ZEROS.md),
[THE_MIRROR_TRANSVERSAL_CERTIFICATE](../../experiments/THE_MIRROR_TRANSVERSAL_CERTIFICATE.md),
[XY_FROZEN_BAND](../../experiments/XY_FROZEN_BAND.md) (−4γ̄ along the XY band,
`(1,1)` included, at depth ⌊N/2⌋), [WHAT_THE_R90_LOCUS_BUYS](../../experiments/WHAT_THE_R90_LOCUS_BUYS.md)
(the corner's 8γ̄·P_D); no null result and no flight on this object.
`docs/CAUGHT_ERRORS.md`: its 2026-08-04 entry records F143 itself as a narrower
re-derivation of a lemma the repo already held, Lemma 5 above, which is the
shape to guard against here; the N = 11 Seed record there is the `(1,2)` pencil of F89, a
different object from this `(1,1)` room. `docs/GLOSSARY.md`: nothing on the
locus, the compressed density or the frozen divisor. OpenArcs: the arc
`compressed_density_laws`; the F-number census arc, which records F154's typed
carrier and nothing on the endpoints; the F89 seed arcs are the `(1,2)` object
again.
Confirmations: nothing. The typed layer: `FrozenDivisorClaim` (F140) with
`FrozenDivisorWitness` (`inspect --root divisor`), `SeedRungGramClaim` (F143)
with `SeedRungGramWitness` (`inspect --root seedrung`), `PinnedBlockFloorClaim`
(F153) with `PinnedBlockFloorWitness`, `AbsorptionTheoremClaim`,
`F71AntiPalindromicGammaSpectralInvariance` and `JointPopcountSectors`;
MirrorWorld's `Divisor` counts F140's frozen modes by exact GF(p) ranks and
holds nothing on the compressed rooms.

Checked for adjacency, pairs about one object that did not cite each other: F140
and the lower endpoint (the endpoint is F140's root, below); Lemma 5 and F143
against the zero-frequency room (on the uniform chains the room's compression is
their Gram on every locus profile, below); F154's own parity theorem and the endpoint room (a
scalar-parity room, below); F144's chiral transpose and the contact matrix of
the N = 11 room (V_(1,6) = V_(6,11) and V_(3,7) = V_(5,9) are chiral-transpose
pairs, below); F146 and PROOF_SCALAR_COUNT §7 against the N = 11 door (the same
coincidence e₁ + e₉ = e₅ + e₆; the `(1,1)` rooms that mix parity appear exactly
at the rung-2 resonant M of that section's table, measured to N = 30); and
PROOF_R90_FROZEN_DIVISOR Lemma 2 against the decomposition of D.

## The two questions

In the Pauli convention, take the uniform open XY chain at N = 11,
H = JΣ(XX+YY), J > 0. Its one-excitation modes are
ψ_k(z) = √(2/12) sin(kπ(z+1)/12), with ε_k = 4J cos(kπ/12).
The `(1,1)` operator cells |a⟩⟨b| form an orthonormal basis. A site-disagreement
operator N_l multiplies |a⟩⟨b| by 1 if exactly one of a,b equals l and by 0
otherwise. The size operator is N_XY = Σ_l N_l: its physical `(1,1)` cells
have disagreement class 0 when a = b and class 2 when a ≠ b, whose uniform
rate centres are 0 and −4γ̄. Write P_Ω for the orthogonal projector onto a
**complete** fixed frequency space of ad_H and C_l = P_Ω(N_l−N_{10−l})P_Ω.

F154's projection identity on a mirror-balanced profile assumes C_l = 0.
The first question is whether that premise persists at N = 11. The second is
whether the `(1,1)` compressed dissipator can nevertheless leave the interval
of its two size-class centres, [−4γ̄, 0]. They have different answers.

## An exact failure of the projection identity

Let X = |ψ₁⟩⟨ψ₆| and Y = |ψ₅⟩⟨ψ₉|. The cosine identity
cos(π/12)−cos(6π/12) = cos(5π/12)−cos(9π/12) places both at frequency
ω = 4J cos(π/12). Reflection gives them opposite parities:
R X R = −X and R Y R = +Y. Neither an accidental frequency near-match nor a
generic diagonal contrast is involved.

For distinct mode indices on both legs, the physical-cell action gives

    ⟨ψ_c|⟨ψ_d| N_l |ψ_a⟩|ψ_b⟩
      = −2 ψ_c(l)ψ_a(l)ψ_d(l)ψ_b(l),

where the double-ket notation represents the Hilbert-Schmidt dyad basis and
all sine modes are real. At the left endpoint this yields

    ⟨Y,N₀X⟩ = −2 ψ₅(0)ψ₁(0)ψ₉(0)ψ₆(0) = −√2/144.

Reflection changes the sign at the right endpoint, so
⟨Y,N₁₀X⟩ = +√2/144 and **⟨Y,C₀X⟩ = −√2/72 ≠ 0**.
The same-frequency space is four-dimensional, with dyads ordered
(1,6), (3,7), (5,9), (6,11). Direct cell summation gives

    C₀ = −√2/72 · [[0,1,1,0], [1,0,0,1],
                     [1,0,0,1], [0,1,1,0]],

with rank 2 and spectrum {−√2/36, 0, 0, +√2/36}. The exact gate constructs
that matrix from the 121 physical cells, checks the complete frequency
membership, and fails under a wrong cell action. A same-parity pair with a
nonzero individual N₀ connection gives a separate zero-C₀ control.

The same formula decides every other `(1,1)` door at once. Since
ψ_k(N−1−z) = (−1)^(k+1) ψ_k(z), the mirror site carries the same product times
(−1)^(a+b+c+d), so ⟨y,C₀x⟩ = −2ψ_cψ_aψ_dψ_b(0)·[1 − (−1)^(a+b+c+d)], which is
−4ψ_cψ_aψ_dψ_b(0) for every opposite-parity pair x = (a,b), y = (c,d). No sine
mode vanishes at site 0, and two distinct dyads in one frequency room differ on
both legs (equal bras and equal frequencies force equal kets). So **every
parity-mixed `(1,1)` room fires, at every N**. Where such rooms exist is
arithmetic: measured on N = 2..30 (census section (B), frequencies grouped at
40 digits), at N = 11, 14, 17, 20, 23 and 29 (6, 12, 6, 6, 12 and 14 rooms, all
firing), and nowhere else, in particular not at N = 8, 10 or 12.

The failure is physically active. Put γ₀ = 2g, γ₁₀ = 0, every other γ_l = g,
with g > 0. This is a nonnegative mirror-balanced profile with γ̄ = g. As
D = −2Σ_l γ_l N_l, its compression is

    P_Ω D P_Ω = −2g P_Ω N_XY P_Ω − 2g C₀.

The uniform term has zero X-to-Y entry, while
⟨Y,D X⟩ = +g√2/36. Thus F154's **conditional** identity
P_Ω D P_Ω = −2γ̄ P_Ω N_XY P_Ω is false on this physical locus profile.
At the uniform profile its correction vanishes. The result decides this
N=11 `(1,1)` premise, not every profile or every block at N=11.

## A separate interval theorem on the `(1,1)` block

The interval survives for a reason independent of C_l. Let h be a Hermitian
one-excitation Hamiltonian with **simple** eigenvalues and site-reflection
symmetry. Let γ_l ≥ 0 and γ_l+γ_{N−1−l} = 2γ̄. This includes the uniform open
XY chain above, at any N and nonzero J. On the `(1,1)` operator space set
Γ = diag(γ_l) and T = Σ_l γ_l |ll⟩⟨ll|. The physical disagreement law is

    D = −2(Γ⊗I + I⊗Γ) + 4T,

where T is positive semidefinite. Also D is diagonal in the physical cells
with entries −2Σ_l γ_l N_l ≤ 0.

Use the complete orthonormal dyads |ψ_a⟩⟨ψ_b| in a fixed frequency space Ω.
An off-diagonal matrix element of Γ⊗I requires the bra modes to agree; equal
frequencies and simple one-body energies then force the ket modes to agree as
well. The other one-sided term is identical with the legs exchanged. Each
eigenmode has reflection-even density, so mirror balance gives
⟨ψ_a|Γ|ψ_a⟩ = γ̄. Therefore, as an **operator identity on Ω**,

    P_Ω(Γ⊗I + I⊗Γ)P_Ω = 2γ̄ I_Ω,
    D_Ω := P_Ω D P_Ω = −4γ̄ I_Ω + 4P_Ω T P_Ω.

The second line and T ≥ 0 give D_Ω ≥ −4γ̄ I_Ω. Compressing D ≤ 0 gives
D_Ω ≤ 0. Hence **spec(D_Ω) ⊂ [−4γ̄, 0]** for every complete `(1,1)` frequency
space under these hypotheses, including the N=11 space above where C₀ ≠ 0.
The two lines above use mirror balance only; the positivity γ_l ≥ 0 enters
through T ≥ 0 and D ≤ 0. This theorem proves containment, not endpoint
attainment or F154's projection identity; containment is a statement about the
compression, and at small J the `(1,1)` Liouvillian can have rates below −4γ̄
(below).

The same one-sided cancellation identifies the contrast's physical carrier.
Write P_ll = |ll⟩⟨ll| in the doubled one-excitation cell space and
m = N−1−l. The compression of
(n_l−n_m)⊗I + I⊗(n_l−n_m) is zero, because its matrix is diagonal in the
dyad frequency basis and each diagonal mode density is reflection-even.
Consequently

    C_l = −2P_Ω(P_ll−P_mm)P_Ω.

This is the **agreement** or double-occupancy projector already used at
uniform rate and zero frequency in
[PROOF_FROZEN_BAND_SO4](PROOF_FROZEN_BAND_SO4.md) §6 (F143); here it is
rate-weighted and compressed at any frequency. Its off-diagonal matrix
elements are signed coherent overlaps on the physical population cells,
not probabilities. The interval theorem needs none of that closed form;
the endpoints below are where it returns.

## The endpoints: F154's scalar room, Lemma 5's Gram, F140's root

The zero-frequency room of a simple h is Ω₀ = span{P_k}, P_k = |ψ_k⟩⟨ψ_k|.
Every P_k is reflection-even, so R acts on Ω₀ as the scalar +1, and F154's
parity theorem gives C_l = 0 there. F154's identity is therefore a theorem on
Ω₀, for every simple reflection-symmetric h, and the endpoint question on this
room is F154's own attainment question. It has an exact answer. Write
W_lk = |ψ_k(l)|², whose rows l and N−1−l agree because each mode density is
reflection-even. In the P_k coordinates, P_Ω₀ T P_Ω₀ = Wᵀ Γ W, and pairing each
site with its mirror turns this into γ̄ WᵀW on every locus profile. So

    D_Ω₀ = −4γ̄ (I − G),   G = WᵀW,

**profile-blind**: the compression hears the mean rate and nothing else, and
it needs no sign on the rates. Both ends follow.

- **Upper, 0.** I = Σ_k P_k is the one-excitation identity, diagonal in the
  physical cells, so D I = 0. In coordinates this is G𝟏 = 𝟏, W being doubly
  stochastic.
- **Lower, −4γ̄.** At most ⌈N/2⌉ rows of W differ, so rank W ≤ ⌈N/2⌉ and
  ker G = ker W has dimension at least ⌊N/2⌋; on it D_Ω₀ = −4γ̄. On random real
  and complex reflection-symmetric h, N = 2..11 (census section (C)), the
  multiplicity reads ⌊N/2⌋ and 0 reads simple, the next eigenvalue in each case
  at least 9·10¹² rounding units away.

On the uniform open XY chain G is Lemma 5's law at M = N + 1, in F143's form
G = (1/M)(𝟏𝟏ᵀ + (I+R_χ)/2), so the whole room is known in closed form:
spec D_Ω₀ = {−4γ̄ ×⌊N/2⌋, −4γ̄(1 − 1/M) ×(⌈N/2⌉−1), 0 ×1}. On the isotropic
Heisenberg chain Lemma 5's law is the site-indexed B Bᵀ, which the mode-indexed
G = BᵀB matches only in spectrum; that is enough for the room, whose spectrum is
the same closed form with M = N (the two square matrices B Bᵀ and BᵀB share
their spectrum), while the matrix identity is an XY statement. At N = 11 on XY that is
{−4γ̄ ×5, −(11/3)γ̄ ×5, 0 ×1}; the witness checks D_Ω₀ = −4(I − G) entry by
entry, exactly, at the balanced (2,1,…,1,0) and at the signed
(5,−2,1,…,1,4,−3), both with γ̄ = 1. On XY the kernel is spanned by the chiral differences
Q = P_k − P_(M−k): since ψ_(M−k)(z) = (−1)^z ψ_k(z), Q has zero physical
diagonal, TQ = 0, and D_Ω₀ Q = −4γ̄ Q. This is F145's seed and the chiral kernel
PROOF_FROZEN_BAND_SO4 §6 names; the chiral density cancellation was recorded for
the N = 5 flip image in [PROOF_CODIM1_BY_ADDITIVITY](PROOF_CODIM1_BY_ADDITIVITY.md) §6.

Both endpoint values are also exact eigenvalues of the finite-J `(1,1)`
Liouvillian block. The sector projector I commutes with every number-conserving
H and is diagonal under dephasing, so the generator annihilates it at every J.
For real h, F140 puts −4γ̄ in the block's spectrum with multiplicity at least
⌊N/2⌋ at every J; at N = 11 on the balanced profile the kernel of L + 4γ̄ is
exactly 5 at J = 3/10, 1/2, 1, 37/10 and 25, F140's lower bound met by an upper
bound from exact ranks modulo two primes (census section (D)). At J = 3/10 the
sixth singular value is 1.2·10⁻¹⁰, inside any rounding tolerance and still
nonzero, following F140's distance ladder J^(2d) with d = 10: the outermost pair
departs last. What is a compression statement is the
containment, not the ends.

**Room by room.** The chiral transpose τ(a,b) = (M−b, M−a) of F144 keeps a
dyad's frequency, since ε_(M−k) = −ε_k, and keeps its contact amplitudes
V_{z,(a,b)} = ψ_a(z)ψ_b(z), by the chiral sign law above. For each non-fixed
orbit {x, τx}, the vector e_x − e_τx is annihilated by V and so by T_Ω, and it sits
at −4γ̄. Hence in every room, for every balanced profile,

    mult(−4γ̄) ≥ #non-fixed chiral-transpose orbits.

Equality holds exactly in all 55 N = 11 rooms at the balanced profile, where the
lower endpoint is attained in 45 of them (the witness, by exact ranks), and it
is measured in every room at N = 2..14 (census section (A), three random locus
profiles each, the next eigenvalue at least 10¹³ rounding units away). The four-dimensional N = 11 room above
is the example: V_(1,6) = V_(6,11) and V_(3,7) = V_(5,9) are its two
chiral-transpose pairs, V has rank 2, and −4γ̄ has multiplicity two. For the
explicit profile (2,1,...,1,0), direct physical-cell compression gives

    spec(D_Ω) = {−4, −4, −10/3−√2/18, −10/3+√2/18},

with characteristic polynomial (x+4)²(x² + 20x/3 + 1799/162), so the upper
endpoint 0 is absent from that room. At N = 11 on the balanced profile 0 lies
in the zero-frequency room only (the witness), and so it reads at every tested
profile of section (A).

Three hypotheses have teeth. Off the locus, γ₀ = 1 and all other rates zero
make the one-sided Γ expectations vary between modes; on the same N=11
frequency space, ⟨Y,DY⟩ = −(22+5√3)/72 < −4/11, so Rayleigh forces a
compressed eigenvalue below the centre interval calculated with γ̄ = 1/11.
The displayed diagonal is a Rayleigh value, **not** asserted to be an
eigenvalue. Without a simple spectrum the one-sided cancellation fails: at
N = 6, two mirror-image dimers on the sites (0,1) and (4,5) with coupling 1 and
a middle dimer on (2,3) with coupling 0.3 give a reflection-symmetric h with
levels ±1 doubly degenerate, and on the locus profile (2,2,1,1,0,0) its complete
four-dimensional frequency rooms at ω = ±0.7 and ±1.3 reach −6 < −4γ̄ = −4
(census section (E)). At finite J the interval need not hold for the whole `(1,1)`
Liouvillian: for balanced rates (2,2,1,1,1,1,1,1,1,0,0), the physical cell
|0⟩⟨1| has dissipative rate −8 at J = 0, below −4. The characteristic roots of
the finite-dimensional block depend continuously on J, so an eigenvalue with
real part below −4 persists for all sufficiently small J > 0; census section
(D) reads the lowest real part at −7.995, −7.923 and −6.639 at J = 0.05, 0.2 and
1. F122 describes the
strong-coupling compression limit; it does not promote this interval to
arbitrary J.

## Scope of the closure

The N=11 example closes the first parity door left by the mixed-space proof:
a parity-odd equal-frequency dyad connection does reach the **physical** site
density, and on the `(1,1)` block every larger-N door fires the same way. The
`(1,1)` interval has an independent all-N proof under simple
reflection-symmetric one-body dynamics and physical mirror-balanced rates. Its
two ends are attained for every such h in the zero-frequency room, where they
are F154's, Lemma 5's and F143's on the uniform XY chain and, for real h, F140's. What stays open: the parity doors
of the higher excitation blocks, the upper endpoint room by room at profiles
that leave pairs dark, and the M₀ arithmetic of the N = 6 transversal
certificate.
