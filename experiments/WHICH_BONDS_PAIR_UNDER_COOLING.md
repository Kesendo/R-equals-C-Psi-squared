# Which Bonds Pair Under Cooling

<!-- Keywords: amplitude damping, cooling alone, T1, palindrome, spectral pairing, two-site bond, loss and gain
exchange, turns about Z, chain reversal, odd power sum, Sextant, Theodolite, F137 -->

**Status:** a classification on the open chain with one bond repeated, 2026-10-26. The carriers are derived below
and checked exactly, the turned and reversed ones each with a control that fails, in
[`simulations/cooled_bonds_pairing.py`](../simulations/cooled_bonds_pairing.py) (27 checks, ALL PASS); the breaks
are exact on the settings named and a double-precision survey everywhere else, and the page says which is which.

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
scope); MirrorWorld SextantTests (S_XY under loss; DM under loss read as symmetric modulo p) and TheodoliteTests; the
scripts cube_old_questions_gate.py G7, cube_moves_gate.py M10/M11, pauli_weight_conjugation.py's 18 configurations.
No store held a census of two-site bonds under cooling, and none claimed necessity outside the class of F137's
criterion. Adjacency: every carrier used here is one of F137's general forms; F137 states them for XXZ bonds with
axial DM terms (w = 0 below) and for H commuting with X⊗N and Y⊗N. New here: the split of J into the parts that see
the difference and the sum of two site angles, which says exactly where turns exist; the one rotation that joins
Jᵀ and F·J·F for every real J, which makes the chain's reversal a carrier wherever the rates allow; and the census.
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
| Δ ≠ 0, N ≥ 3, s, a, w all nonzero, other rates | no, on every setting computed | none |
| a mixed letter (XZ, ZX, YZ, ZY) anywhere in the bond | no, on every setting computed | none |

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
that does not vanish rules out every element. With s, a and w all nonzero and rates the reversal moves, the first
that is not zero is the seventh on every setting computed (XX + 2XY + ZZ at (2, 1, 1), (1, 2, 3) and, N = 4, (1, 2, 3, 4);
2XX + 3XY − YX + 2YY + ZZ at (1, 2, 3)). A bond with a mixed letter breaks at the fifth (XZ; ZX + XX); F137's
YZ − ZY has 240 there. That these break everywhere outside the carriers is
a survey, not a proof: in double precision, with one random integer bond in ±{1, 2, 3} for each of the 511 subsets
of the nine letters, read as paired when every eigenvalue's mirror lies within 10⁻⁶, every bond with a mixed letter
misses its mirror (480 of 480 at each of six rate profiles, N = 2, 3, 4, by at least 0.06), every bond of the class
pairs at N = 2 and at rates the reversal keeps ((1, 1, 1), (1, 2, 1), (1, 2, 2, 1)), and at (1, 2, 3) and
(1, 2, 3, 5) the bonds drawn with s, a, w all nonzero break, 8 of 8 (one of the nine subsets with ZZ drew a
vanishing part), while the other 23 pair.

## What the instruments did

The Sextant decides by the same odd power sums, modulo a prime, so the exact nos are its kind of answer (computed
here in Python, section D of the gate); its search runs over uniform moves only and cannot see the turned yeses.
The Theodolite reads those as exact elements off the letters (section C), and is silent on the reversal, which no
product of site maps is. The cube gave the why: the
split of a bond into a part that sees the difference of two site angles and a part that sees their sum, and the
reversal's one rotation, which a computation then confirmed exactly.

## What stays open

Whether every bond with s, a, w all nonzero breaks at every N ≥ 3 and every rate profile the reversal does not
keep; whether a bond with a mixed letter can pair at all under cooling on some graph or rate profile (F137's open
sentence, for this chain); bonds that differ along the chain, other graphs (a ring adds a cycle to the angle
conditions, a star a centre to the reversal), fields, and a thermal bath at a ratio other than zero.
