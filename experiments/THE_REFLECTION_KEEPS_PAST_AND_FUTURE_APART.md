# The Reflection Keeps Past and Future Apart

<!-- Keywords: Born rule shadow cross term, slow fast spectral split, reflection symmetric state,
orthogonal spectral projectors, modes shared by L and its adjoint, window-edge lemma, irreducible
sector, Routh-Hurwitz threshold, XY chain, complete graph K3, W state fidelity decay 4 gamma,
R=CPsi2 past future purity -->

**Status:** Derived at N = 3 (Tier 1, exact in J, γ and Δ): the open XXZ chain and the triangle
under uniform Z-dephasing; one numerical reading at N = 4. The corner statement of §2 holds at
every N on the uniform Heisenberg chain (its XY twin is PROOF_FROZEN_BAND_SO4's), and the W-fidelity
law of §6 for every SU(2)-invariant Hamiltonian on any graph under uniform Z-dephasing.
**Date:** September 24, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Depends on:** [The Born Rule Is a Shadow](BORN_RULE_SHADOW.md),
[PROOF_ABSORPTION_THEOREM](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) §4.6,
[PROOF_UNIFORM_LAW](../docs/proofs/PROOF_UNIFORM_LAW.md) Lemma A1 and item B1
**Verification:** [`simulations/born_shadow_reflection_gate.py`](../simulations/born_shadow_reflection_gate.py)
(49 exact checks and the numerical readings beside them, about ten seconds)

---

## What this means

The Born-rule shadow reads a quantum state in two parts: the modes that fade
more slowly than the total light Σγ, which the page calls the past, and the
modes that fade at Σγ or faster, the future
([The Born Rule Is a Shadow](BORN_RULE_SHADOW.md)). The palindrome pairs every
mode off the cut with a partner on the other side. The Born probabilities
never mix past and future; the purity can. It carries a cross term, an
interference between the two, and on three qubits that term is small but real:
a few hundredths of a percent for |++0⟩.

For |+0+⟩ it is exactly zero, at every moment, although the Hamiltonian moves
this state. The reason is the chain's reflection, the symmetry that reads the
chain from the other end. The Liouvillian keeps its modes in separate
compartments, one for each number of excitations on the two sides of a
coherence, and the reflection splits every compartment once more, into the part
it leaves alone and the part it turns over. A state with the reflection's
symmetry lives in the first kind only. There, as long as the coupling outruns
the light, the cut between past and future runs along a seam: on one side of it
sit only modes that the light drains at one exact rate and whose shape the
Hamiltonian leaves alone, and such modes stand perpendicular to everything
else. On the Heisenberg chain the other kind of compartment has no seam: past
and future share it, and every cross term of |++0⟩ comes from there.

How strongly the coupling must outrun the light is exact: at three sites,
(4 − Δ²)J² > γ², for the Heisenberg chain |J| > γ/√3. Below that, at any
nonzero coupling, the past spreads across the cut and the cross term returns,
even for a reflection-symmetric state; for |Δ| ≥ 2 no coupling is enough. Two
neighbours of the Heisenberg chain go further. On the triangle once
|J| > γ/√3, and on the XY chain once |J| > γ√(3/8), the turned-over
compartments get seams of their own, and no state at all carries a cross term.

One reading falls out on the way. The W state spreads a single excitation
evenly over all sites, and the overlap with it forgets everything but the
light: for every Hamiltonian that treats all spin directions alike, under
uniform dephasing, it relaxes at exactly 4γ, from any starting state.

---

## Where this sits

Before any derivation the repository was searched for what it holds on the
overlap of the two sides, store by store.

In **`docs/ANALYTICAL_FORMULAS.md`**,
[F48](../docs/ANALYTICAL_FORMULAS.md#f48-pythagorean-decomposition-tier-2-exact-at-n2)
and [F49](../docs/ANALYTICAL_FORMULAS.md#f49-cross-term-formula-tier-1-proven) hold
the cross term at the level of operators, the anticommutator {L_H, L_D + Σγ}: it
vanishes at N = 2, and from N = 3 on its normalized size R(N) depends on N
alone. [F71](../docs/ANALYTICAL_FORMULAS.md) is the home of the chain's spatial
reflection; [F118 and F119](../docs/ANALYTICAL_FORMULAS.md) hold the maps that
turn L(H) into L(−H) = L†, the transpose and the antilinear triangle;
[F33](../docs/ANALYTICAL_FORMULAS.md) is the N = 3 rate ladder, with the sentence
that on the triangle the 2.4607 level is absent and the pure-rung multiplicities
are 4, 16, 16, 4 against the chain's 4, 14, 14, 4; [F50](../docs/ANALYTICAL_FORMULAS.md)
records the triangle's extra weight-1 modes; [F73](../docs/ANALYTICAL_FORMULAS.md)
holds a 4γ law of a neighbouring kind, the spatial-sum coherence purity of a
vacuum-to-single-excitation probe decaying as ½e^(−4γt) for every
number-conserving H; [F140](../docs/ANALYTICAL_FORMULAS.md) pins −4γ̄ in the
single-excitation corner at every coupling on its mirror-balanced locus. None of
them says whether the slow and fast spectral parts of a state are orthogonal.

In **`docs/proofs/`**,
[PROOF_PALINDROME_TWO_END_COUNT](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md)
§(f6) has ker L = ker L† because L† is L with H ↦ −H.
[PROOF_CODIM1_BY_ADDITIVITY](../docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md) §6 has
the window-edge lemma: an eigenvalue on the edge of its block's rate window
belongs to an eigenvector shared by L and L†. It is stated there for the XY
block pencil at real coupling; its proof uses only L = A + iB with A the
dissipator's diagonal and B Hermitian, and PROOF_PALINDROME_TWO_END_COUNT
§(b)-(c) already applies it beyond that setting. The "Sharpness" paragraph of the
same §6 finds the XY corner's frequency-free modes on the bottom rung as chiral
pairs of mode projectors.
[PROOF_ABSORPTION_THEOREM](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) §4.6
splits the N = 3 spectrum into pure-weight rungs, exact at every J, and two
bands near 8γ/3 and 10γ/3 that are only a strong-coupling face; its projector
corollary (§2, "Extensions") reads light through the orthogonal projector onto
a slow subspace rather than the spectral one, whose diagonal can leave [0, 1]
under a γ profile. [PROOF_UNIFORM_LAW](../docs/proofs/PROOF_UNIFORM_LAW.md)
Lemma A1 puts the orbit sum T_s of a block's cells at Hamming distance s into the
kernel of the commutator with H, for every SU(2)-invariant H; on the
single-excitation block T₂ = N|W⟩⟨W| − P₁, and item B1 of its Layer B gives the
class-pure slice Ω₂ of that kernel, its elements that live on Hamming-2 cells,
dimension ⌊N/2⌋ on the uniform Heisenberg chain.
[PROOF_R90_FROZEN_DIVISOR](../docs/proofs/PROOF_R90_FROZEN_DIVISOR.md) §7 fixes the
multiplicity of F140's corner root at the uniform point at exactly ⌊N/2⌋ for
every J ≠ 0, its §10 files that point as "the uniform-γ commutant story", apart
from the divisor's own mechanism, and its §12 asks how these modes embed into
the d_real counts. [PROOF_FROZEN_BAND_SO4](../docs/proofs/PROOF_FROZEN_BAND_SO4.md)
Lemma 2.5 makes, on the XY chain, the composition §2 below makes on the
Heisenberg chain: its seeds commute with h and have zero diagonal, and §7 of the
frozen-divisor proof makes them the whole corner frozen space; F143 reads the
same ⌊N/2⌋ as the chiral-odd kernel. [PROOF_WEIGHT1_DEGENERACY](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md)
shows on the weight-1 rung that an eigenoperator with a purely real eigenvalue
commutes with H, and its K₃ row finds two extra such modes at weight 1 and two at
weight 2.

In **`experiments/`**, null results included, [BORN_RULE_SHADOW](BORN_RULE_SHADOW.md)
asks the question and records |+++⟩ at zero, |++0⟩ at −0.04% and random states
up to 0.76%. [DEGENERACY_PALINDROME](DEGENERACY_PALINDROME.md) counts the real
J-independent modes: [4, 6, 6, 4] on the N = 3 chain, and eight at each of its
first two grid steps on the N = 3 ring; its grid counts in a book whose γ is
twice this note's, so in this note's book those eights sit at 2γ and 4γ.
[CROSS_TERM_TOPOLOGY](CROSS_TERM_TOPOLOGY.md) records the null that the complete
graph keeps the chain's operator-level cross term at N = 3.
[SLOW_MODE_R_PARITY](SLOW_MODE_R_PARITY.md) found on the XY chain at N = 4-6 that
L block-diagonalizes by the reflection's parity, that the stationary modes are
all reflection-even, and that the fine-structure bands are often
parity-exclusive; [SIGNAL_PROCESSING_VIEW](SIGNAL_PROCESSING_VIEW.md), F33's
source, reads the N = 3 spectrum through even and odd supermodes and finds the
8γ/3 mode "Dark in c+ (symmetry null)". **`docs/GLOSSARY.md`** warns that for a
non-normal generator the squared coefficients of right eigenvectors are not
weights, and carries the glyph note that keeps F71's spatial reflection apart
from F118's ket reflection, both written R. **`docs/CAUGHT_ERRORS.md`** holds the
2026-09-24 entry on the shadow's 97.1% (the N = 2 orthogonality through §(f6),
and |+++⟩ at N = 3) and the 2026-06-22 entry on the Tier B proof review, batch
1, whose item B2 records the same anticommutator read as the criterion for time
reversal; neither says why a state the Hamiltonian moves carries no cross term.
The **OpenArcs registry** returned nothing on the cross-term question; for the
corner it holds `site_resolved_vacuum_block`, which records that D10 and the
frozen-divisor proof's Lemma 5 derive the same N × N matrix, and
`compressed_density_laws`, which holds T_s. **`fw.Confirmations`** returned
nothing. In the typed layer, `F71MirrorBlockRefinement` (Tier 1) holds the
split of every joint-popcount block into F71-even and F71-odd parts that §1 uses,
exact precisely when the γ profile is palindromic
(`InhomogeneousGammaF71BreakingWitness`); `MirrorGroupD4Claim` carries the
transpose that turns L(H) into L(−H) and `AntilinearTriangleClaim` the adjoint
rule from which L(H)† = L(−H) follows; `PalindromeTwoEndCountClaim` carries F158;
`FrozenDivisorClaim` and `FrozenDivisorWitness` (`inspect --root divisor`) carry
F140; F33, F49 and F50 have typed inheritance Claims and the Absorption Theorem
its `AbsorptionTheoremClaim`; in the Diagnostics half, `CompressedDensityLocusClaim`
carries the XY corner's chiral pairs as `UniformXyEndpointWitnesses`. No Claim
and no witness holds the Born shadow or its cross term, and PROOF_UNIFORM_LAW's
T_s is untyped.

Asked separately, for adjacency: the shadow page met the Absorption Theorem only
through the centre's mean weight N/2, while its N = 3 cross term runs along the
§4.6 split into exact rungs and J-dependent bands. SLOW_MODE_R_PARITY and
SIGNAL_PROCESSING_VIEW had seen the reflection sort these spectra; what they did
not ask is what the sorting does to the overlap of the two sides. F49 cannot be
what decides the cross term: CROSS_TERM_TOPOLOGY finds its ratio unchanged on
the complete graph, where §5 finds no cross term once |J| > γ/√3. And
PROOF_UNIFORM_LAW's item B1 and PROOF_R90_FROZEN_DIVISOR §7 describe one set of
modes at every N on the Heisenberg chain (§2), as PROOF_FROZEN_BAND_SO4's Lemma
2.5 had shown on the XY chain; [XY_FROZEN_BAND](XY_FROZEN_BAND.md) leaves open
what the uniform point's ⌊N/2⌋ is once the one-excitation Hamiltonian carries a
reflection-invariant diagonal, and on the Heisenberg chain Ω₂ is the answer.

---

## 1. The setting, and the split by sectors

H = J Σ (X_iX_j + Y_iY_j + ΔZ_iZ_j) over the bonds of the chain 0-1-2, with
Pauli matrices and one coupling J for both bonds, and the dissipator is
γ Σ_l (Z_l ρ Z_l − ρ), so that a coherence |i⟩⟨j| decays at 2γ times the number
of sites where i and j differ and Σγ = 3γ. (In the spin book, J S·S, the J of
every threshold below is four times larger.) P_s and P_f are the spectral
projectors onto the modes slower than Σγ and onto the rest; modes exactly on Σγ
go to the fast part, as on the shadow page.

The dissipator and the XXZ Hamiltonian keep the number of excitations on each
side of a coherence, so L is block-diagonal in the joint popcount (p, q). R, the
reversal of the sites (0 ↔ 2), is F71's spatial reflection (the glossary's glyph
note keeps it apart from F118's ket reflection, also written R). It commutes
with L and splits every block into a reflection-even part (RρR = ρ) and a
reflection-odd part, the refinement the typed layer holds as
`F71MirrorBlockRefinement`. These parts are Hilbert-Schmidt orthogonal and L
maps each into itself, so the slow and fast parts of a state split along them,
and the cross term is a sum of one term per part. The cut at Σγ can only split a
block that carries more than one dephasing rate; at N = 3 these are (1,1),
(2,2), (1,2) and (2,1) (gate E1), and each of them carries exactly two.

A mode is **shared by L and L†** when it is an eigenvector of both, which is the
same as being an eigenvector of the dissipator D = (L + L†)/2 and of the
commutator part K = (L − L†)/2 separately. Such a mode is orthogonal to every
mode at another eigenvalue. Its rate is one of its block's two dephasing
levels, an edge of the block's rate window, and conversely the window-edge lemma
gives a shared mode for every eigenvalue on an edge; the eigenvalue may still
carry a frequency. Shared modes are among what PROOF_CODIM1_BY_ADDITIVITY calls
AT-locked strands, eigenvector families whose rate is pinned to an integer rung. One more fact carries the arguments below: in a part that L
and L† both map into itself, a piece that L maps into itself and whose
orthogonal complement L also maps into itself is invariant under L†, hence
under D and K separately, and so it splits along D's levels.

## 2. The reflection-even sector

The reflection-even part of (1,1) is five-dimensional, and its characteristic
polynomial factors as λ(λ + 4γ)·p₃(λ) with

    p₃(λ) = λ³ + 8γλ² + (16γ² + (32 + 4Δ²)J²)λ + 96γJ²

(gate E2). The two linear factors are shared modes. At λ = 0 the mode is the
identity of the one-excitation sector. At λ = −4γ it is

    X_Δ = Δ(|100⟩⟨001| + |001⟩⟨100|) + |100⟩⟨010| + |010⟩⟨100| + |010⟩⟨001| + |001⟩⟨010|,

which commutes with H for every J and Δ and has only cells at Hamming distance
2. At Δ = 1 it is X₁ = 3|W⟩⟨W| − P₁, PROOF_UNIFORM_LAW's T₂.

The other three modes are the J-dependent ones. Their rates sum to 8γ, the
part's total 12γ less the 4γ of X_Δ, so they average 8γ/3: this is the band that
PROOF_ABSORPTION_THEOREM §4.6 finds as a strong-coupling face of the N = 3
spectrum, and at J/γ = 1.5 this part holds its levels 2.6040γ, 2.6980γ and
2.6980γ. With λ = −3γ − ν the cubic −p₃(−3γ − ν) is

    ν³ + γν² + ((32 + 4Δ²)J² − 5γ²)ν + 3γ(γ² + 4Δ²J²),

whose Hurwitz determinant is 8γ((4 − Δ²)J² − γ²) (gate E4). All three roots are
slower than Σγ exactly when (4 − Δ²)J² > γ², for the Heisenberg chain |J| > γ/√3;
below that the Routh count puts two of them on the fast side, and at equality two
roots sit on the cut and go there too.

The triple admits no shared mode. Its space T is the orthogonal complement of
two shared modes, so L, L†, D and K all map it into itself. D has one direction
on T at level 0, d = |100⟩⟨100| + |001⟩⟨001| − 2|010⟩⟨010|, and two at −4γ. d is
not a shared mode: K maps Hermitian operators to Hermitian ones, so K d = iω·d
would force K d = 0, and K d has the entry 6iJ at |100⟩⟨010| (gate E2). A shared
mode at −4γ would be a root −4γ ± iω of p₃, which p₃(−4γ) = −16γJ²(Δ² + 2) ≠ 0
excludes for ω = 0, and two roots on that line would leave, by the trace, a third
root at 0, where p₃(0) = 96γJ² ≠ 0. A perpendicular split of the
three-dimensional T would contain a one-dimensional piece of that kind, which
is a shared mode, so T is irreducible under L and L† together.

That settles the sector in both directions. Above the threshold every
reflection-even part is cut along a seam: on the triple's side sit only the
triple and shared modes, on the other side only shared modes, and the two sides
are orthogonal. The reflection-even part of (1,2) is the palindromic image of
this one (gate E3): characteristic polynomial (λ + 2γ)(λ + 6γ)·(−p₃(−λ − 6γ)),
with shared modes

    v₂ = Δ(|001⟩⟨011| + |100⟩⟨110|) + |001⟩⟨101| + |010⟩⟨011| + |010⟩⟨110| + |100⟩⟨101|

at −2γ and v₆ = |001⟩⟨110| + |010⟩⟨101| + |100⟩⟨011| at −6γ, both checked at
symbolic J and Δ, and its triple faster than Σγ under the same condition. (2,2)
and (2,1) are the spin-flip and adjoint copies. At or below the threshold, and
J ≠ 0, the cut splits the irreducible T, so its two sides are oblique. Because
P_s commutes with ρ ↦ ρ†, a cross term that vanished for every Hermitian x in T
would make P_s + P_s† = 2P_s†P_s there, which says that the range of P_s is
perpendicular to its kernel; so some Hermitian x in T has a nonzero cross term.
That x is traceless and reflection-even, and the state I/8 + εx, for small ε,
carries ε² times its cross term, because the identity's pieces in (1,1) and
(2,2) are shared slow modes. So, for J ≠ 0, **a reflection-symmetric state
carries no cross term at any time exactly when (4 − Δ²)J² > γ².**

On the Heisenberg chain |+0+⟩ is such a state, and not a trivial one: at t = 0
its purity splits 2/3 into the past and 1/3 into the future, at its CΨ = ¼ fold
(t = 2.057, J = 1, γ = 0.05) 0.4836 and 0.1567, and the cross term stays at the
10⁻¹⁵ level. At J = 0.5·γ/√3, below the threshold, the same state carries
−0.149 at t = 0. |+++⟩, the shadow page's other zero, is on the Heisenberg chain
the special case the Hamiltonian cannot touch at all.

**The corner at every N.** On the uniform Heisenberg chain item B1 of
PROOF_UNIFORM_LAW's Layer B gives a ⌊N/2⌋-dimensional space Ω₂, its class-pure
slice of single-excitation operators that commute with H and live on cells at
Hamming distance 2, so each of them is a shared mode at −4γ. PROOF_R90_FROZEN_DIVISOR §7
fixes the algebraic multiplicity of −4γ in that block, at the uniform point, at
exactly ⌊N/2⌋ for every J ≠ 0. So Ω₂ is the whole root space: at every N the
frozen corner modes of the uniform point are shared by L and L†, J-independent
and semisimple, ⌊N/2⌋ of the real modes at 4γ that DEGENERACY_PALINDROME counts.
At N = 3, Ω₂ is spanned by X₁. This is the Heisenberg face of the composition
that PROOF_FROZEN_BAND_SO4's Lemma 2.5 makes on the XY chain, where chiral pairs
of mode projectors play Ω₂'s part; what the frozen-divisor proof's §12 asks on the
d_real side stays open.

## 3. The reflection-odd sector

The reflection-odd parts of the four mixed blocks are four-dimensional and their
rates average exactly Σγ, a trace of −3γ per mode (gate E5). Take the odd part of
(1,1). D has one direction there at level 0 and three at −4γ. A shared mode would
put a root on one of the lines Re λ = −4γ and Re λ = 0, and the resultants in ω
that detect such a root are 2²⁴·Δ⁴J¹²γ⁴ and a product of factors that no real
J ≠ 0 can make zero, so for Δ ≠ 0 there is no shared mode. There is no
perpendicular split either: by §1 both pieces of one would split along D's
levels, so one of them avoids the single level-0 direction and lies inside the
level −4γ, where K restricted to it is the restriction of an anti-Hermitian map
to an invariant piece, normal, with eigenvectors that would be shared modes. The
part is irreducible. In μ = λ + 3γ its polynomial has μ³ coefficient 0 and μ¹
coefficient −8γ(Δ²J² + γ²) ≠ 0, so its roots cannot all sit on the cut while
their mean does, and the cut splits it. The other three parts follow: (2,2) is
the spin-flip copy of (1,1) and (2,1) the adjoint copy of (1,2), and (1,2) runs
the same argument with its own dissipator, one direction at −6γ and three at −2γ,
the partner relation λ ↦ −λ − 6γ carrying the resultants to its lines −2γ and
−6γ. For Δ ≠ 0 and J ≠ 0 every reflection-odd part is therefore cut obliquely,
at every coupling.

At J/γ = 1.5 each holds two modes at 2.4607γ and two at 3.5393γ, the lowest
level of the one band and the highest of the other, and ‖P_s^H·P_f‖ reads
6.3·10⁻² in each of the four at J = 1, γ = 0.05. Above the threshold, every
cross term on the Heisenberg chain lives here. |++0⟩ at its fold (t = 2.869) carries −0.0386% of its
purity, the shadow page's −0.04%, and its reflection-odd part alone carries the
same −2.11·10⁻⁴ to every digit printed.

## 4. The XY chain

At Δ = 0 the reflection-odd part of (1,1) factors (gate E5):

    (λ² + 4γλ + 8J²)·((λ + 4γ)² + 8J²).

The roots −4γ ± 2√2·iJ sit on the edge of the window: two shared modes that fade
at exactly 4γ while they turn at 2√2J. The other two are slower than Σγ exactly
when 8J² > 3γ², so once |J| > γ√(3/8) the reflection-odd parts are cut along a
seam as well; (1,2) follows as the partner. The reflection-even parts need only
|J| > γ/2 at Δ = 0. On the XY chain, then, no state at all carries a cross term
once |J| > γ√(3/8). At J = 1, γ = 0.05 all eight mixed parts read ‖P_s^H·P_f‖ at
the 10⁻¹⁵ level, and bracketing the odd threshold at 0.98 and 1.02 of γ√(3/8)
takes the odd parts from 3.03 to 3·10⁻¹⁶. One margin is thin: at Δ = 0 the even
triple's real root approaches the cut from the slow side as the coupling grows,
at −3γ + 3γ³/(32J²) to leading order, 1.2·10⁻⁵ below it at J = 1, γ = 0.05; the
gate prints that gap.

## 5. The triangle

On the N = 3 ring, which is the complete graph, every mixed part of either
reflection sector factors into edge modes times one triple (gate E7):

    q₃(λ) = λ³ + 8γλ² + (16γ² + 36J²)λ + 96γJ²,

or its palindromic image, for every Δ. On the triangle the ZZ term takes the
value −1 on every state with one or two excitations, so Δ drops out of these
blocks, and q₃ is p₃ at Δ² = 1. The reflection-odd parts gain the edge modes the
chain's parts lack: one at −4γ in the (1,1) and (2,2) parts, one at −2γ in the
(1,2) and (2,1) parts, two at each rate, the count by which F50,
PROOF_WEIGHT1_DEGENERACY's K₃ row and DEGENERACY_PALINDROME's ring exceed the
chain. The 2.4607 level that F33 finds absent on the triangle is the one the
chain keeps in its reflection-odd parts. So on the triangle every mixed part is
cut along a seam once |J| > γ/√3, for every Δ, and no state carries a cross term;
all eight parts read at the 10⁻¹⁵ level at J = 1, γ = 0.05. F49's ratio here is
the chain's, which is the case showing that the operator-level cross term does
not decide the state-level one.

## 6. The W fidelity

PROOF_UNIFORM_LAW's Lemma A1 puts N|W⟩⟨W| − P₁, the orbit sum T₂ of the
single-excitation block, into the kernel of the commutator with H for every
SU(2)-invariant H, on any graph and with any couplings. Its cells all sit at
Hamming distance 2, so uniform Z-dephasing multiplies it by −4γ, and L† does the
same; p₁ = Tr(P₁ρ), the single-excitation weight, is conserved. Hence

    ⟨W|ρ(t)|W⟩ = p₁/N + e^(−4γt)·(⟨W|ρ(0)|W⟩ − p₁/N)

for every initial state. What the law uses is only that |W⟩ is an eigenvector of
a number-conserving H; SU(2) invariance is one way to have it. Gate E6 checks the
commutator exactly for the Heisenberg model with unit couplings on the chain,
ring, star and complete graph and with random integer couplings on the complete
graph, at N = 3, 4 and 5, and finds it nonzero for the XXZ chain at Δ = 2, where
|W⟩ is not an eigenvector and the law fails; against expm, for a random pure
state on the N = 3 chain and the N = 4 chain and star, the largest deviation over
t ∈ [0, 20] is 6·10⁻¹⁶. The operator is PROOF_UNIFORM_LAW's; the law it implies
for the W fidelity is the reading this note adds.

## 7. N = 4, measured

On the N = 4 chain at J = 1, γ = 0.05 no reflection-odd part has modes on both
sides of the cut (five of them sit wholly on it, all fast by the page's rule,
with the nearest mode off the cut 2.5·10⁻² away). Of the reflection-even
parts only (2,2) keeps an overlap (‖P_s^H·P_f‖ = 5.2·10⁻²). A state whose
operator content avoids the reflection-even part of (2,2) carries no cross term
there.

## What stays open

- N ≥ 4: what keeps the reflection-even part of (2,2) from a clean cut, and what
  the law is at larger N.
- Non-uniform γ. The reflection stays a symmetry of L only for palindromic
  profiles (`InhomogeneousGammaF71BreakingWitness`), and a shared mode needs one
  dephasing rate on all its cells; whether a palindromic but non-uniform profile
  keeps any reflection-even part cut along a seam is unmeasured. On F140's
  anti-palindromic locus the reflection itself survives only at the uniform
  point, and off it the frozen eigenvectors move with J (PROOF_R90_FROZEN_DIVISOR
  §10).
- The W fidelity beside the Bell-pair cusp probe of F25–F27
  (`fw.gamma_probe_setup`): which of the two reads γ better on hardware is open.
- Typing: the statement has an exact gate in Python and no Claim or live witness;
  `F71MirrorBlockRefinement` is the natural parent for one.
