# Which Bonds Pair Under Cooling

<!-- Keywords: amplitude damping, cooling alone, T1, palindrome, spectral pairing, two-site bond, loss and gain
exchange, turns about Z, chain reversal, odd power sum, mixed letter, ring, Sextant, Theodolite, F137 -->

**Status:** a classification on the open chain with one bond repeated, 2026-10-26, and the same questions on
rings. The carriers are derived below and checked exactly, the turned and reversed ones each with a control that
fails, in [`simulations/cooled_bonds_pairing.py`](../simulations/cooled_bonds_pairing.py) (93 checks, ALL PASS; 103
with --deep, which adds the exact polynomials at larger N that the closed forms for every N rest on, and one more point modulo a prime); the breaks are exact where the seventh power sum below
decides them, on the chain's other profiles where the ninth does, for a mixed letter where the fifth does, and on the settings named; the rest is a double-precision
survey or open, and the page says which is which.

<!-- STAGE 0 SWEEP RECORD (for us, not for the reader: hidden on GitHub, read by every session, agent and hook)

Swept: docs/ANALYTICAL_FORMULAS.md F137 (S_XY for every H commuting with X⊗N and Y⊗N, XYZ bonds with any couplings;
S_IZ and the loss and gain criterion for H commuting with Z⊗N; the turns about Z on trees; a site permutation
composed with turns; the conjugation with a sign on bipartite graphs; the exact break of YZ − ZY on the open chain
of three, Tr((L + c/2)⁵) = 240 under cooling alone, read there as the formal setting at zero splitting; its open
sentence on Hamiltonians that do not commute with Z⊗N), F82, F84, F158 (the size of the Π violation; Pauli jumps
only), F1 (its DM clause on bipartite graphs); docs/proofs/MIRROR_SYMMETRY_PROOF.md's T1 paragraph (sufficient
conditions, no necessity), PROOF_F84, PROOF_F155 (nothing on pairing); experiments/ THERMAL_BREAKING,
KMS_DETAILED_BALANCE (scope rows only), NON_HEISENBERG_PALINDROME and THE_PALINDROME_AS_A_COLOURING (the bond census
under Z dephasing and the frame theorem, none under cooling), THE_CUBE_SITE_BY_SITE (finding 1, the element onto
−L − c at ratio zero; the elements here are onto −L† − c); docs/THE_ONE_SQUARE.md §9 (amplitude damping's arrow and
the turns that survive damping, F63 and the Z⊗N partnership); review/OPEN_QUESTIONS_INDEX.md (OQ-168, detailed
balance; no question on bonds); the OpenArcs registry (one measured cooling plus Z reading, no classification arc);
docs/CAUGHT_ERRORS.md (the scoping of "T1 alone keeps the palindrome", with its counterexamples: transverse fields,
transverse DM bonds, longitudinal profiles); the typed layer (F137 untyped; DirectSumDecompositionClaim repeats its
scope); MirrorWorld SextantTests (S_XY under loss; DM under loss witnessed through the theodolite) and TheodoliteTests; the
scripts cube_old_questions_gate.py G7, cube_moves_gate.py M10/M11, pauli_weight_conjugation.py's 18 configurations.
No store held a census of two-site bonds under cooling, and none claimed necessity outside the class of F137's
criterion. Adjacency: every carrier used here is one of F137's general forms; F137 states them for XXZ bonds with
axial DM terms (w = 0 below) and for H commuting with X⊗N and Y⊗N. New here: the split of J into the parts that see
the difference and the sum of two site angles, which says exactly where turns exist; the one rotation that joins
Jᵀ and F·J·F for every real J, which makes the chain's reversal a carrier wherever the rates allow; the census; and
the seventh power sum in closed form. For that last part the sweep went to the odd power sums the repo holds in closed
form: docs/proofs/PROOF_MOMENT_TOWER_PUMP_CHANNEL.md (the odd moments under Z dephasing, odd in γ by a word reversal,
the class P_{5,3} = 20·4^N·[...] a sum over weight-three Z strings and sites, the same counting of local words with
another symmetry doing the work); its sibling PROOF_F87_WINDOWED_MONOMIAL_CONVERSE §4 with the girth dichotomy
m* = 2ℓ + deg (the first odd power sum that fires set by the shortest odd closed word and the γ degree, scripts
f87_girth_dichotomy.py, typed as GirthLadder in compute/RCPsiSquared.Diagnostics/F87 and the claims
WindowedConverseAllGammaClaim and MomentTowerPumpChannelClaim), the same shape as the phase count here, with the
shortest turn-invariant word odd under H ↦ H̄ in place of the odd closed walk;
THE_MIRROR_UNDER_A_SHARED_FIELD (low odd sums vanishing while the spectrum does not pair, the K3 pair), the neighbour
of the ninth-sum remark; F117 (the odd power sums as the trace face of hardness), F137's single values (240 at the fifth for
YZ − ZY, −656·sin Φ on the damped triangle, a seventh of −263277/400 there with X dephasing). None is a closed form
under cooling alone, and none reads a rate profile against its reverse. For the mixed letters and the rings
the same stores again: F137's open sentence on Hamiltonians that do not commute with Z⊗N and its YZ − ZY value 240,
which the closed form for the mixed letters reproduces (the page says so), F87's positivity of the first odd sum as a sum of squares (PROOF_F87_WINDOWED_MONOMIAL_CONVERSE §5,
the same shape as S here, under Z dephasing), and no store with a power sum on a ring under cooling alone (F137's damped triangle and its 4 ring carry heating
too). Adjacency checked: F137's XXZ bonds with axial DM terms pair at one ratio on any graph, ratio zero included,
so w = 0 bonds pair on rings at any flux, consistent with the ring forms here vanishing at w = 0. For the ninth cumulant on the chain: no store holds a
cumulant of a Liouvillian spectrum or a power sum read as a sum over windows of a rate profile; the nearest,
PROOF_MOMENT_TOWER_PUMP_CHANNEL's P_{5,3} and its census values of p₉, are uniform rates under Z dephasing, the same
counting of local words with no profile to read, so adjacent in method and not in object.
-->

Put cooling σ⁻ = |0⟩⟨1| on every site of an open chain, each site at its own rate, and nothing else: no heating,
no dephasing, no field. Couple neighbours by one bond repeated along the chain, a sum of two-letter Pauli terms.
When is the spectrum of the generator symmetric about −c/2, c the sum of the rates?
[F137](../docs/ANALYTICAL_FORMULAS.md#f137) answers yes for every H commuting with X⊗N and Y⊗N (the element
S_XY: ρ ↦ X⊗N·ρ·Y⊗N, XYZ bonds with any couplings), derives the cases of XXZ bonds with axial DM terms, shows one
bond that breaks (YZ − ZY), and leaves the rest open. Here is the rest, for this chain.

**Which setting this is.** XXZ bonds with axial DM terms conserve the excitation number and are read, as in F137,
in the frame rotating with the qubits. Every other bond below does not conserve it (an anisotropy XX − YY or
XY + YX, a mixed letter), so no rotating frame removes the qubits' frequency, and cooling alone on such a bond is
F137's formal setting at zero splitting: the generator read as written. The physical bath at zero splitting has
equal cooling and heating rates, and there every H commuting with Z⊗N pairs (F137). The classification below is
of the formal generator.

## The class and the criterion

A bond with a mixed letter, X or Y on one site and Z on the other, does not commute with Z⊗N. Every bond without
one does: J_ab·σ_a⊗σ_b over a, b ∈ {X, Y} plus Δ·ZZ, the real 2×2 matrix J carrying the X and Y letters. For this
class F137 holds a criterion: the shift S_IZ: ρ ↦ ρ·Z⊗N carries L onto the reflection of the Liouvillian with
loss and gain exchanged, so the spectrum pairs exactly when exchanging them changes nothing. The spin flip X⊗N
exchanges them back and turns H into X⊗N·H·X⊗N, which in this class, with real couplings, is H̄, every Y letter's
sign reversed. So the spectrum pairs exactly when spec L(H) = spec L(H̄), and one answer is a unitary V that keeps
the dissipator, with V·H̄·V† = H; the element is then Ad_V∘S_XY, onto −L† − c (L† the adjoint in the
Hilbert–Schmidt inner product).

Write J = s·1 + a·ε + T, with ε = [[0, 1], [−1, 0]] and T traceless and symmetric: s is the trace part, a the
antisymmetric part (the axial DM term, the flux of the hopping), T the anisotropy; set z = s − i·a and
w = T₁₁ + i·T₁₂. The conjugation H ↦ H̄ sends J to F·J·F, F = diag(1, −1), that is z to z̄ and w to w̄. The unitaries
that keep σ⁻ up to a phase on every site are the turns about Z, which rotate a site's X and Y letters by an angle
θ_l (X ↦ cos θ_l·X + sin θ_l·Y); they send the bond (l, l+1) to R(θ_l)·J·R(θ_{l+1})ᵀ with
R(θ) = [[cos θ, −sin θ], [sin θ, cos θ]], which turns z by θ_l − θ_{l+1} and w by θ_l + θ_{l+1}: the two parts
see the difference and the sum of the angles. ZZ is untouched by both.

## The answer

| bond on the chain (rows 1 to 5 without a mixed letter), cooling at any rates per site | pairs | carrier |
|---|---|---|
| Δ = 0, any J | always | F137's conjugation with a sign: σ⁻ is real, so spec L(H) = spec L(−H̄) for every H, and Z on every other site sends −H̄ to H̄ |
| Δ ≠ 0, a = 0 (J symmetric), s = 0 (J traceless, z = 0 included) or w = 0 (XXZ + DM) | always | turns about Z: Ad_V∘S_XY |
| Δ ≠ 0, N = 2 | always | turns about Z |
| Δ ≠ 0, N ≥ 3, s, a, w all nonzero, rates the chain's reversal keeps | yes | the reversal times one uniform turn, then S_XY |
| Δ ≠ 0, N ≥ 3, s, a, w all nonzero, other rates | no where γ₁γ₂ ≠ γ_{N−1}γ_N, at every N, and at N = 3 wherever γ₁ ≠ γ₃; with equal end products wherever the ninth cumulant below is not zero, and on the two settings computed where it is | none |
| a mixed letter (XZ, ZX, YZ, ZY) anywhere in the bond | no where the fifth power sum below is not zero, at every N: always without anisotropy (w = 0) at rates positive on every site; at a site without cooling it can pair (XZ at (1, 0) and (1, 0, 1)); elsewhere no on every setting computed | none |

**The turns.** Turns undo the conjugation on a bond when z·e^{i(θ_l − θ_{l+1})} = z̄ and w·e^{i(θ_l + θ_{l+1})} = w̄.
With z and w both nonzero the difference and the sum of the angles are fixed on every bond, and two neighbouring
bonds then give θ_l − θ_{l+2} = 0 from the sums and 2(θ_l − θ_{l+1}) from the differences, so 2·arg(z̄/z) vanishes
modulo 2π: z² is real, a = 0 or s = 0. Where it is, angles alternating between two values do it; where z or w is
zero only one of the two conditions is there and a turn from site to site meets it. One bond fixes the sum and the
difference of its two free angles independently, so N = 2 always pairs. Since the turns are all the site maps that
keep σ⁻ up to a phase, this is exactly where a product of them carries the pairing. F137's tree turn,
θ_b = θ_a − arg(J + iD_z), is the turn W that makes H real; V here is its square, W⁻² (the angle differences are
2·arg z), and the condition above extends it from w = 0 to the whole class. The Theodolite reads these elements off the
letters, its site factor [[0, x], [1, 0]] carrying the turn x: a real X/Y part (the letters XX and YY only) needs no turn
(S_XY itself, F137's element for H commuting with X⊗N and Y⊗N), an imaginary one (XY and YX only) alternates by a
half turn (x changing sign from site to site),
and a symmetric J turns every site alike (at XX + 2XY + 2YX + ZZ, x = ±(−1 + 4i)/√17 on every site).

**The reversal.** Where the turns cannot, the chain's reversal P can, an instance of F137's site permutation
composed with turns. It sends J to Jᵀ, and Jᵀ and F·J·F have the same trace, the same antisymmetric part and
traceless symmetric parts of one length, so one rotation R takes the first to the second, R·Jᵀ·Rᵀ = F·J·F, for every
real J. With U the uniform turn by R, U·P·H·P†·U† = H̄, and where the rates are symmetric about the chain's middle
P keeps the dissipator, so V = (U·P)† carries it. Exact at N = 3 and 4 on the rates (1, 2, 1), (1, 1, 1),
(1, 2, 2, 1) and on J = [[3, −1], [2, 5]] with −2·ZZ at (2, 7, 2), the carrier failing at (2, 1, 1) and (2, 7, 3); the first
three are XX + 2XY + ZZ, where the rotation is by the phase of −1 + 2i, cos = −1/√5, sin = 2/√5.

**The breaks.** The odd centred power sums Tr((L + c/2)^k), k odd, all vanish wherever the spectrum pairs, so one
that does not vanish rules out every element. With Δ, s, a and w all nonzero, the first that is not zero is the
seventh wherever γ₁γ₂ ≠ γ_{N−1}γ_N; the next section says why and gives it in closed form (section D of the gate asserts it at XX + 2XY + ZZ at (2, 1, 1), (1, 2, 3) and, N = 4, (1, 2, 3, 4), and at
2XX + 3XY − YX + 2YY + ZZ at (1, 2, 3)). A bond with a mixed letter breaks at the fifth wherever its closed form below is not
zero (section D asserts it at XZ and at ZX + XX). The survey, in double precision, with one random integer bond in ±{1, 2, 3} for each of the 511 subsets
of the nine letters, read as paired when every eigenvalue's mirror lies within 10⁻⁶, every bond with a mixed letter
misses its mirror (480 of 480 at each of six rate profiles, N = 2, 3, 4, by at least 0.06; the closed form decides
most of them, and those where anisotropy cancels the fifth are survey only), every bond of the class
pairs at N = 2 and at rates the reversal keeps ((1, 1, 1), (1, 2, 1), (1, 2, 2, 1)), and at (1, 2, 3) and
(1, 2, 3, 5) the bonds drawn with s, a, w all nonzero break, 8 of 8 (one of the nine subsets with ZZ drew a
vanishing part, and both profiles have unequal end products, so the closed form below requires these breaks), while
the other 23 pair.

## Where the pairing has to fail

**The reversal decides.** In the class, the spin flip exchanges loss and gain back and turns H into H̄, and U·P
carries H̄ to H while moving each site's rate to the mirror site; so L(H̄) at the rates γ is, up to that unitary, L(H)
at the reversed rates Pγ. F137's criterion then reads: for every J and Δ, the spectrum pairs exactly when reversing
the rate profile leaves the spectrum unchanged. Write p_k(γ) = Tr((L + c/2)^k). The power sums are real (L keeps Hermiticity,
so its spectrum is closed under conjugation), and S_IZ gives p_k at γ as (−1)^k times p_k of the loss and gain
exchange (F137: every odd centred trace changes sign under the exchange), which is p_k at Pγ. So p_k(Pγ) = p_k(γ) for even k, and no even power sum can tell a profile from its
reverse, while p_k(Pγ) = −p_k(γ) for odd k.

**Why the seventh.** Each factor L + c/2 is linear in the couplings and the rates together, so p_k is homogeneous
of degree k in them. With every bond's couplings taken as their own variables, each odd power sum is odd under
H ↦ H̄ (the step above, which needs no reversal); since σ⁻ is real, complex conjugation gives the spectrum of L(H)
as the conjugate of that of L(−H̄), so p_k(H) = p_k(−H̄) and, with the first, p_k is odd in H as a whole; it is even
in J alone (Z on every other site flips J and keeps Δ), so odd in Δ; and it is unchanged by the turns about Z. A turn gives z
on the bond (l, l+1) the phase θ_l − θ_{l+1} and w the phase θ_l + θ_{l+1}. On one bond the only products these
phases leave alone are |z|² and |w|², which H ↦ H̄ keeps; the first product it changes, with the phases cancelling,
has four factors on two neighbouring bonds, z·w̄ on the first and z·w on the second. For a uniform bond that is
z²·|w|², and its part that H ↦ H̄ reverses is Im(z²)·|w|² = −2s·a·|w|². With an odd power of Δ beside these four J factors the
degree in H is at least five and odd; the rates enter at least once (a commutator's spectrum pairs about 0 on its own) and, the total degree
being odd, an even number of times. So p₁, p₃ and p₅ vanish at every N, the first that can break is the seventh,
and it needs two bonds: N = 2 always pairs.

**What the seventh sees.** Read p_k/4^N as the k-th moment of the spectrum. Since p₁, p₃ and p₅ vanish, the
seventh moment equals the seventh cumulant, and the cumulants add over connected clusters of local terms (the
linked-cluster theorem: terms on disjoint sites commute, the trace of the exponential factorizes, and its logarithm
keeps no cross term). A cluster that contributes carries the four J factors on two
neighbouring bonds, Δ on a bond touching them and two rates on its sites: at most four sites. So p₇ is 4^N times a
sum over windows of three and four neighbouring sites, each window's polynomial in its own rates the same wherever it
sits; the chain of three fixes the three-site window, the chain of four then fixes the four-site one (it comes out
zero), and that gives

  Tr((L + c/2)⁷) = −7·2¹⁴·4^(N−3)·Δ·s·a·|w|²·(γ₁γ₂ − γ_{N−1}γ_N).

Each three sites l, l+1, l+2 contribute γ_{l+1}(γ_l − γ_{l+2}), and along the chain these telescope: the seventh
sum sees only the two ends. Section E of the gate interpolates the polynomial in the rates modulo primes at N = 3
and 4 for three bonds (and at N = 2 for one, where it is zero), compares it exactly, and checks the closed form at N = 5 and 6 at one random
profile for each of two of the bonds, modulo a prime.

At N = 3 the sum is a multiple of γ₂(γ₁ − γ₃): with s, a, w and Δ nonzero and γ₂ ≠ 0 the spectrum pairs exactly
when γ₁ = γ₃. From N = 4 on it breaks wherever γ₁γ₂ ≠ γ_{N−1}γ_N. On the profiles with equal end products that the
reversal moves, and at N = 3 with no cooling in the middle, the seventh is silent and the ninth has to decide.

**What the ninth sees.** Off the seventh's zeros the ninth moment carries a product of the seventh and second
cumulants, p₉ = 4^N·κ₉ + 36·p₇·p₂/4^N (no other product survives, the third and fifth cumulants being zero), and
the ninth cumulant κ₉ is again a sum over windows. The phase count fixes how the couplings enter. A product of
J factors that H ↦ H̄ reverses has at least four of them, so with p₉ odd in H the couplings enter five or seven
times and the rates four or two (nine would leave no rate, and a commutator alone pairs). With five, the only product is the seventh's, Δ·Im(z²)·|w|²; with seven, the
products are that one times Δ², |z|² or |w|², and Δ·Im(z⁴)·|w|², whose reversed part is again a multiple of
s·a·|w|²·Δ. So the whole of κ₉ is Δ·s·a·|w|² times a polynomial in the rates whose quartic part is the same for
every bond and whose quadratic part is linear in Δ², s², a² and |w|². A window with four rates spans at most four
sites and is the same for every bond, so one bond fixes it. A window with two spans at most six: its J part lies in
the span of Δ², s², a² and |w|² up to five sites, and on six only the four J factors with Δ on each of the three
other bonds reach, so there it is Δ² alone. Exact polynomials at N = 3 and 4 for six bonds (four with independent
parts fix the quadratic windows, two more agree), at N = 5 for four and at N = 6 for one therefore fix every window,
those of five and six sites come out zero, and

  κ₉ = 144·Δ·s·a·|w|²·[Q(γ) + E(γ₁, γ₂, γ₃) − E(γ_N, γ_{N−1}, γ_{N−2}) + 16·(e(γ₁, γ₂, γ₃) − e(γ_N, γ_{N−1}, γ_{N−2}))],

  Q(γ) = Σ_{l=1}^{N−1} 2γ_lγ_{l+1}(γ_l² − γ_{l+1}²) + Σ_{l=1}^{N−2} γ_lγ_{l+2}·[11γ_{l+1}(γ_l − γ_{l+2}) − 7(γ_l² − γ_{l+2}²)],

  E(x, y, u) = xy·(38x² − 3xy + 40y²) − 34x²yu,

  e(x, y, u) = −6(s² + a² + |w|²)·x² + (16s² + 40a² − 40|w|²)·y² − (111Δ² + 226s² + 374a² + 66|w|²)·xy
               − (23Δ² + 33s² + 181a² − 61|w|²)·yu − (7Δ² + 196s² − 28a² − 28|w|²)·xu,

from N = 3 on (section H of the gate, N = 5 and 6 with --deep; also at one random bond with every part nonzero and a profile
with unequal end products at N = 4 to 6 modulo a prime, and 7 with --deep). Every part changes
sign under the reversal, as it must. The quadratic part sits at the ends only, since a quadratic sum over windows
that the reversal reverses telescopes, the argument of the rings. The quartic part does not: Q is a sum over the
whole chain, each bond and each pair of next neighbours contributing a term odd under its own reversal, and it does
not telescope. The seventh reads the profile at the two ends only; the ninth reads it everywhere.

Two things follow. At N = 3 with no cooling in the middle,

  κ₉ = −144·Δ·s·a·|w|²·(γ₁² − γ₃²)·(7γ₁γ₃ + 96·(s² + a² + |w|²)),

and the last factor is positive at rates not negative, so for every bond with Δ, s, a, w nonzero the ninth breaks
wherever γ₁ ≠ γ₃. With the seventh, at N = 3 the spectrum pairs exactly when γ₁ = γ₃, whatever the middle rate,
which closes N = 3 (for XX + 2XY + ZZ at (1, 0, 2) the gate also has the ninth not zero, exactly). From N = 4 on the profiles
with equal end products break wherever κ₉ ≠ 0, which for a given bond leaves a surface of one dimension less, and
the surface is not empty off the reversal's profiles: the closed form picks −6XX − 2XY + YX + YY − 2ZZ at
(1, 6, 3, 2) and −6XX − 6XY + 2YX + 2YY + 3ZZ at (1, 8, 1, 8), both with κ₉ = 0, and there the seventh and ninth are
exactly zero and the eleventh is not (in the gate; XX + 2XY + ZZ at (1, 6, 2, 3), where κ₉ ≠ 0, has the ninth
not zero, also in the gate).

## The bonds with a mixed letter

Write the mixed part of the bond in complex form, μ·(σ⊗Z) + ν·(Z⊗σ) with μ = c_XZ + i·c_YZ and ν = c_ZX + i·c_ZY
(c the coefficients). A turn about Z gives μ on the bond (l, l+1) the phase θ_l and ν the phase θ_{l+1}. These
letters anticommute with Z⊗N, so S_IZ no longer applies and nothing makes the odd power sums odd under H ↦ H̄; the
turns and complex conjugation remain. The phase count now finds smaller invariants: on one bond μ̄·ν·z and μ·ν·w̄,
across two bonds ν_b·μ̄_{b+1}. A half turn on every site flips the mixed letters and keeps the others, so they enter
every power sum with an even degree. For every bond of the nine letters p₃ vanishes (exact at N = 2 and 3, where its
clusters, at most two bond terms and one rate, all fit), so the first that can break is the fifth, and

  Tr((L + c/2)⁵) = 10·4^(N−2)·[3·S − 16·Im(μνw̄)·Σ_b γ_bγ_{b+1}],

  S = Σ_{m=2}^{N−1} γ_m·|μ·γ_{m+1} + ν·γ_{m−1}|² + |μ|²·γ₁γ₂² + |ν|²·γ_{N−1}²·γ_N

    = Σ_b (|μ|²·γ_bγ_{b+1}² + |ν|²·γ_b²γ_{b+1}) + 2·Re(μν̄)·Σ_m γ_{m−1}γ_mγ_{m+1},

the second line its windows: one per bond and one per three neighbouring sites.

The fifth moment is the fifth cumulant here (p₁ and p₃ vanish), a contributing cluster has five terms with at least
one rate and so spans at most five sites, and the polynomial in all nine coefficients and the rates is exact at
N = 2 to 5 (section F of the gate, N ≥ 3 with --deep; at N = 3 to 6 also at random points modulo a prime). Of the
class's parts only the anisotropy w enters, once, through the invariant μ·ν·w̄; that μ̄·ν·z, Δ, s and a do not
enter is what the exact polynomials say, not something the phase count forces.

F137's break of YZ − ZY on three sites at rate 1 is an instance: μ = i, ν = −i, w = 0, so the square vanishes,
the end terms give S = 2, and the fifth sum is 10·4·3·2 = 240.

Without anisotropy S alone remains, a sum of squares weighted by the rates plus two end terms that are positive as
soon as μ or ν is not zero and the rates are. So at rates positive on every site no bond with a mixed letter and
w = 0 pairs, at any N (a site without cooling can empty S, and then the spectrum can pair: XZ alone at rates (1, 0) and at (1, 0, 1)
pairs, exact by the characteristic polynomial, while (1, 1, 0) does not): the X/Y part
absent, or XXZ with an axial DM term. With anisotropy the cross term can cancel the squares: XZ + ZX − 3XY − 3YX at
N = 2 and rates (8, 8) has p₅ = 0, and the seventh is not zero there (exact, in the gate). A mixed letter with
anisotropy is decided by the fifth wherever 3·S ≠ 16·Im(μνw̄)·Σ_b γ_bγ_{b+1}, and elsewhere by a higher one on the
settings computed.

## Rings

On a ring every site has two neighbours and the chain's ends are gone. The phase count carries over with one
change. Every factor z, z̄, w or w̄ of a bond gives each of the bond's two sites a phase ±θ, and at each site the
factors of its two bonds must cancel, so two neighbouring bonds carry numbers of factors of the same parity. On a
chain the end sites force every number even; round a ring they are all even or all odd, and all odd needs a factor
on every bond, N of them at least. So on a ring of N ≥ 4 sites, as on the chain, a product that H ↦ H̄ reverses
has at least four X/Y factors; p_k being odd in H, its degree in the couplings is odd and at least five, and the
rates enter at most k − 5 times.

The ring brings two symmetries the chain lacks. A rotation only relabels the rates, so p_k keeps every rotation of
the profile; and F137's criterion, read through the reversal as above, holds with a reflection of the ring in place
of the chain's reversal (the uniform turn that carries the reflected bond Jᵀ to F·J·F works on any ring), so every odd p_k changes sign
when the profile is reflected. A polynomial in the rates of degree two or less that keeps every rotation is a
constant plus a multiple of Σγ_l plus a quadratic form whose matrix depends only on the distance round the ring between
two sites, and each of these keeps the reflections too; changing sign as well, it is zero. So on every ring of four
or more sites p₃, p₅ and p₇ vanish at every rate profile, for every bond of the class repeated round the ring, and the ninth is the first
power sum that can read the reflection.

On the triangle one factor per bond closes round the ring. In p₇ a rate degree d leaves the couplings the odd
degree 7 − d; the parts with d ≤ 2 vanish as on the longer rings, and with d ≥ 6 too few factors remain, so d = 4
is the only one left, with exactly one factor on each bond and no Δ. The products that close are z³ and z·|w|² and
their conjugates (a w closes only against a w̄), so the part H ↦ H̄ reverses is α·Im z³ + β·Im z·|w|², with α and
β quartics in the rates, the same for every bond. The triangle's symmetries are all six permutations of its sites,
the reflections being the transpositions; Im z³ and Im z·|w|² being independent functions of J, α and β each
change sign under every transposition, and so each is a number times (γ₁ − γ₂)(γ₁ − γ₃)(γ₂ − γ₃)(γ₁ + γ₂ + γ₃). The seventh of the two bonds 3XX − XY + 2YX + 5YY − 2ZZ
and −2XX + XY + 3YX + YY + 3ZZ is zero as a polynomial in the rates on the triangle and on the ring of four (in
the gate), and their pairs (Im z³, Im z·|w|²) are independent, so both numbers are zero: the seventh vanishes on the
triangle for every bond of the class as well. At one random class bond on each ring of 3 to 6 sites p₃, p₅ and p₇ are zero
and p₉ is not (modulo a prime). Without the rotation the count still holds: with a different bond on each edge the
triangle's fifth is not zero, its three factors closing round it, while on the rings of four and five the fifth
stays zero and the seventh is not (in the gate). On the ring of four,

  Tr((L + c/2)⁹) = −9·2¹²·Δ·s·a·|w|²·(γ₁ − γ₃)(γ₂ − γ₄)·(Q(γ₁, γ₃) − Q(γ₂, γ₄)),  Q(x, y) = 5x² + xy + 5y²,

exact as a polynomial in the rates for two bonds. The first two factors are the reflections through opposite sites;
the reflections through opposite bonds swap the two pairs and make the last factor vanish, which it also does on a
surface that no reflection keeps, where the ninth is silent; at (5, 3, 9, 10) on that surface, Q(5, 9) = Q(3, 10)
= 575, the ninth is exactly zero and the eleventh is not (in the gate). On the triangle the ninth sum is a nonzero constant
times (γ₁ − γ₂)(γ₁ − γ₃)(γ₂ − γ₃)(γ₁ + γ₂ + γ₃) for the two bonds computed, the constant depending on J otherwise
than through Δ·s·a·|w|² (the two constants are not in that ratio; the triangle is not bipartite, so the sum need
not be odd in Δ); for those bonds the
triangle pairs exactly when two rates agree, which is exactly when a reflection keeps the profile.

The bonds with a mixed letter lose their end terms on a ring too: on the ring of four the fifth sum is the same
expression with S summed round the ring and without the two end terms (at one point modulo a prime; on the
triangle, where three bond terms already close round the ring, it differs, in the gate for one bond), so at equal rates and μ = −ν every square vanishes.
XZ − ZX on the ring of four at equal rates has p₅ = p₇ = 0 exactly, and the ninth is not zero: it does not pair
either. For the class on the surface where the ninth is silent and for XZ − ZX, where the fifth and seventh are, a
higher sum breaks on the settings computed; whether a ring with these bonds ever pairs where its first sums vanish
is open.

## What the instruments did

The Sextant decides by the same odd power sums, modulo a prime, so the exact nos are its kind of answer (computed
here in Python, sections D to G of the gate); its own search runs over uniform moves only, and for the yeses its
decision asks the Theodolite. That reads the turned yeses as exact elements off the letters (section C), and the
reversal after relabelling the sites (`Theodolite.ReadRelabelled`, the element E·P checked densely on the full
generator in MirrorWorld's TheodoliteTests): the reversal keeps every string's letter counts and so every height, and
turns the letters of a bond around. `sextant decide 3 XX:1,XY:2,ZZ:1 1,2,1 loss` reads the kept rates as witnessed and
`… 1,2,3 loss` as obstructed at the seventh trace. The cube gave the why: the
split of a bond into a part that sees the difference of two site angles and a part that sees their sum, and the
reversal's one rotation, which a computation then confirmed exactly; and the phase counts that make the seventh
power sum the first that can break for the class and the fifth for a mixed letter, which the exact polynomials then
gave in closed form.

## What stays open

Whether every bond with Δ, s, a, w all nonzero breaks on the profiles from N = 4 on where both the seventh and the
ninth cumulant vanish and the reversal does not keep the profile (the eleventh breaks on the two computed); whether a bond with a mixed letter and anisotropy can pair where its fifth power sum
vanishes, and which bonds with a mixed letter pair at sites without cooling (in the formal setting, the counterpart
of F137's open sentence, which asks about the thermal bath); rings: the surface Q(γ₁, γ₃) = Q(γ₂, γ₄) of the ring of four beyond the one point
computed, the triangle's ninth beyond the two bonds computed, the mixed letters on rings beyond the bonds computed, and longer rings; bonds that differ along the chain, other graphs (a ring adds a cycle to the angle
conditions, a star a centre to the reversal; a ring's turns are not examined here), fields, and a thermal bath at a ratio other than zero.
