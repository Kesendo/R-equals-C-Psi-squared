# The Cube, Site by Site

<!-- Keywords: letter cube instrument, site-local palindromizer, product reflector, steady-state weight, M_beta,
KMS detailed balance, temperature gradient, golden router uniform, Heisenberg DM spiral gauge, Groebner exact -->

**Status:** an instrument and three findings, 2026-10-20; finding 1's bound through Z⊗N 2026-10-30. The findings are
reproduced exactly, each with a control that fails, in
[`simulations/cube_site_by_site_findings.py`](../simulations/cube_site_by_site_findings.py) (56 checks, ALL PASS), and
the bound also lives in C# (`MirrorLocalityBoundWitness`, `inspect --root locality`). The least rank at two sites is
computed in [`simulations/kms_element_least_rank.py`](../simulations/kms_element_least_rank.py), exact where it says
so and numerical at the rank-5 point; its exact exclusion of rank 4 also lives in C# (`ThermalElementLeastRankWitness`,
`inspect --root leastrank`). Where a statement below rests on an argument or on a computation rather than on
a check, it says so.

<!-- STAGE 0 SWEEP RECORD (for us, not for the reader: hidden on GitHub, read by every session, agent and hook)

**What the repo already held (sweep).** Whether a product of one-site maps reflects the generator is the house's
per-site routing question: `KBodyPalindromeRouting` and `PalindromeSoftCertifier` (Diagnostics), the golden router
of [PROOF_CEILING_GOLDEN_ROUTER](../docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md) and its `GoldenRouterClaim`, the
router census of [SOFTNESS_IS_N_DEPENDENT](SOFTNESS_IS_N_DEPENDENT.md), the Q families of
[NON_HEISENBERG_PALINDROME](NON_HEISENBERG_PALINDROME.md), all under Z dephasing. Beyond it the repo held S_XY and
its carriers ([F137](../docs/ANALYTICAL_FORMULAS.md#f137)), the thermal element M_β of
[KMS_DETAILED_BALANCE](../docs/KMS_DETAILED_BALANCE.md) Question 4 with an open question beside it and in
[OQ-168](../review/OPEN_QUESTIONS_INDEX.md), and MirrorWorld's `Sextant`, whose search over the uniform one-site
moves finds S_XY in a thermal bath by itself. Swept besides: `docs/proofs/` (the complement connection, the
windowed converse, F103, F108 Part 2, MIRROR_SYMMETRY_PROOF), the Core claim graph (`AntilinearTriangleClaim`,
`F84ThermalAmplitudeDampingPi2Inheritance`, `PalindromeTwoEndCountClaim`), MirrorWorld's `AntilinearTriangle`,
the OpenArcs registry (`f138_converse_failures`, a monomial search never ported), `fw.Confirmations`
(`palindrome_trichotomy`), the glossary (the Π row's per-site families) and `docs/CAUGHT_ERRORS.md`. Checked for adjacency: the instrument's
support-set equations against `KBodyPalindromeRouting` (the same condition, continuous here); its identity row
against the router proof's identity column (§4); its weight against M_β's thermal weight and the transpose law of
`AntilinearTriangleClaim`; its turned maps against F137's turns about Z; its answer with fields against F137's
chained carriers, which compose into it. For the zero-ratio boundary of finding 1: F137's pure-cooling
paragraph (the pairing through the effective Hamiltonian and S_XY, no Π-shaped element), MIRROR_SYMMETRY_PROOF's T1
paragraph (no site-local Π of Π's own shape), KMS Question 4's pure-cooling sentence (the weight singular, no claim
on an element), F82 and F84 (the size of the violation, no element), the OpenArcs registry (nothing on M_β or ratio
zero), OQ-168 and `docs/CAUGHT_ERRORS.md` (finite ratios only), the pins `test_cube_reflector.py` (no cube move
for cooling alone, the machine's NONE on cooling-only chains), TheodoliteTests (no copy on the damped chain) and
`cube_old_questions_gate.py` G7 (the H_eff route), the Core amplitude-damping claims (violation size), Sextant's loss
tests (S_XY a reflection, four uniform reflections, no Π-reflection asserted); nothing contradicted it. Adjacency:
MIRROR_SYMMETRY_PROOF's T1 paragraph and KMS's pure-cooling sentence were about the same object and cited neither
this finding nor each other; both now point here.

For the Z-corner bound (2026-10-30, from KMS Question 4's open sentence): swept by name `docs/proofs/` (PROOF_PALINDROME_TWO_END_COUNT §(f8) holds the similarity L ~ −L − 2σ from the real structure,
stated for Hermitian involutive jumps and never applied to thermal ones; PROOF_ANTILINEAR_TRIANGLE §7 asks only after
the form of the weight; PROOF_CEILING_GOLDEN_ROUTER closes the F87 "6 → 4 → 2" count to zero, Hamiltonian
palindromizers under Z dephasing, no thermal element), `experiments/` (THE_PAIRING_CONDITION lines 157-165, the same
chain, and line 719 keeping EXISTS apart from EXPLAINS; PI_OPERATOR_ENTANGLEMENT, whose Schmidt rank 9 of an
eigenvector-pairing Π was the representative's, not the family's, and NON_HEISENBERG_PALINDROME line 241, the
Z-dephasing twin of the open rank question, where the least rank over the valid-Π space is 1 and the 9 was a search
floor; THERMAL_BREAKING, the X⊗N swap at line 71 and S_XY at line 341; CEILING_FOUR_NONLOCAL_CASES),
`docs/ANALYTICAL_FORMULAS.md` (F137's S_XY line stops at −L† − c; F158's entry names §(f8); F119 nothing),
the OpenArcs registry (`f138_converse_failures` line 4049 "EXISTENCE TURNS OUT TO BE FREE" and line 4098
"explaining nothing", its item FOUR at line 4286, "widen the S search past site-locality", never run;
`handshake_decoder`, whose location dictionary is a Gram matrix on bond space in the single-excitation sector and
not a Liouville-space metric), OQ-168, the Core claims
(`PalindromeTwoEndCountClaim` "which no gate here observes", `QuditProductMirrorCap`'s non-product partial
intertwiners of the dissipator alone, F121), the Diagnostics witnesses (none on a thermal element), MirrorWorld
(`EndCount` leaves §(f8)'s similarity outside by design; `Sextant` searches uniform one-site moves and reads the
pairing; `hidden_pi_search.py`'s eigenvector-pairing Q are all ill-conditioned, two by pseudo-inverse, residuals up
to 6.6e4), `fw.Confirmations`, the glossary and
`docs/CAUGHT_ERRORS.md` (line 5043, the real-structure repair). No tool solved E·L = (−L − c)·E on the full
superoperator space, and no store held a Schmidt-rank bound on an element. Adjacency: KMS Question 4's open sentence
and OQ-168 are about the same object as §(f8) and THE_PAIRING_CONDITION and cited neither; both now point here.

For the least rank (2026-10-30, three scouts by store). Markdown: the open item stood here, in KMS Question 4 and in
OQ-168, nowhere answered; THE_PAIRING_CONDITION lines 157-182 hold the chain L ~ −L† − 2σ ~ −L − 2σ ("the metric
relation is the law") under Z dephasing; PROOF_ANTILINEAR_TRIANGLE §7 holds the weight at one ratio and left open
whether a weight of another kind serves across a gradient (line 158); KMS line 452 writes M_β as Γ_β·T·S_XY;
PI_OPERATOR_ENTANGLEMENT line 46 and NON_HEISENBERG_PALINDROME lines 240-242 record a search floor once mistaken for a
least rank, which is why a rank floor seen only by failing searches (one was seen for the Hermiticity-preserving
elements) is not stated here; nothing on a charge-graded Schmidt rank. By name: `docs/ANALYTICAL_FORMULAS.md` holds
F137 and F158 §(f8) and no F-number on an element's locality; `fw.Confirmations` holds hardware entries only, none on
an element; the glossary has no Schmidt-rank entry; `docs/CAUGHT_ERRORS.md` has none on rank searches. Typed layer:
`MirrorLocalityBoundClaim` named the least rank open and has no child claims; the OpenArcs registry holds nothing on
it; the one Diagnostics witness that splits a Schmidt rank by charge is this finding's own port,
`ThermalElementLeastRankWitness`; F137 holds S_XY's reflection onto −L† − c at any rates and
that thermal jumps keep the coherence order which S_XY reverses; MirrorWorld's
`AntilinearTriangle.ThermalTransposeResiduals` is Γ_β at one ratio, Theodolite reads S_XY under a gradient and is
silent for Π's shape. Adjacency, checked by a third scout and read at the source: PROOF_ANTILINEAR_TRIANGLE line 158
asks after the same object, a weight across a gradient; it now says what this rate point shows and keeps the general
case open, and points here; the coherence order is THE_ONE_SQUARE §10's d (its m is another charge, 2η^z);
THE_CHANNEL_BOTH_WAYS_ON_A_PATH line 135 names XY − YX as the part that is not its own transpose, there in a generator
and not in a state, a rhyme rather than the same object. Nothing adjacent on the rank split by charge. -->

Nothing in the repo held a site map that carries a weight beside the move beyond Z dephasing.

## The instrument

The [letter cube](../docs/THE_ONE_SQUARE.md) reads a Pauli string by its letter counts, which is what makes it
valid at every N and what makes it lose the site, the side of a multiplication, the sign and the coefficient. The
instrument works one level beneath. Per site it takes one of the cube's 192 one-site moves, dressed with a weight
that is diagonal on the cells |i⟩⟨j| (a turn about Z, a thermal weight), before or after the move, and asks
exactly whether the product of these site maps carries the generator L onto −L − c (a **Π-reflection**, the point
reflection of the spectrum about −c/2) or onto −L† − c (a **reflection**, the line reflection across Re = −c/2),
with c = −Σ_l tr L_l / 2 the total rate, dephasing included. On a given graph the question splits into one equation
per site and one per support set of the Hamiltonian. The letters answer it. When every site factor is diagonal or
anti-diagonal, a sandwich ρ ↦ A·ρ^(T)·B with A and B products over the sites, each factor acts on the letters as a
pure turn: it multiplies σ⁺ by a power of its turn or swaps σ⁺ and σ⁻, and keeps or flips Z. Every equation then
becomes binomial, a product of turns equal to a constant read off the rates and the couplings: the site's heights and
arrows give one or two, each bond's letter table the rest, the left and the right factor apart. Such a system is
solved by integer elimination on the exponents and radicals of the constants, in closed form and without a Gröbner
basis, and a step per site in the same equations gives every N: a map the same on every site, turned from site to
site, or periodic, answers on every open chain, and on the rings whose length the period divides (for a turned map,
where the total turn closes). Where the letters are silent, the identity pre-check, the cube search and a general
solve over all invertible site maps (Gröbner bases over ℚ(i)) shake hands: no site-local reflector at all, or one
that is not such a sandwich. Code: [`framework/cube_tables.py`](../simulations/framework/cube_tables.py),
[`framework/diagnostics/cube_letters.py`](../simulations/framework/diagnostics/cube_letters.py) (`read`, `every_n`)
and [`framework/diagnostics/cube_reflector.py`](../simulations/framework/diagnostics/cube_reflector.py) (`answer`,
`every_n`), with their tests. The instrument is the Theodolite, and its C# twin,
[`compute/MirrorWorld/Theodolite.cs`](../compute/MirrorWorld/Theodolite.cs), reads the same letters as equations and
never as numbers: the turns satisfy x² = 3 rather than equal a root, and solvability is decided by integer
elimination alone (`Binomials`). Out of its reach by design: a relabelling of the sites,
a weight along X or Y, an element that is not a product over sites, a term on four sites.

The findings below came from running it on settings taken from the dephasing documents (F1, MIRROR_SYMMETRY_PROOF,
NON_HEISENBERG_PALINDROME, DEPOLARIZING_PALINDROME, F108 Part 2), the thermal ones (F137, THERMAL_BREAKING, KMS
Question 4, the rows M9 to M11 of `cube_moves_gate.py`) and the routing ones (the ceiling cases, the golden router,
SOFTNESS_IS_N_DEPENDENT, THE_PALINDROME_AS_A_COLOURING, the `SextantTests` palindromes, F138's SWAP rows). Those runs
are not part of the repo; what they found is reproduced in the script above.

## What it found

**1. The weight of the thermal element is the one-site steady state.** Take a site with cooling γ↓, heating γ↑ and
dephasing at rate γ⊥ along any axis a·X + b·Y of the XY plane. Its steady state has the population ratio
q² = (γ↑ + γ⊥)/(γ↓ + γ⊥), and the site map ρ ↦ w·Y·ρᵀ·X·w with w = diag(1, q) solves the site equation exactly
there and nowhere else (the residual is proportional to q²(γ↓ + γ⊥) − (γ↑ + γ⊥)). The other half is the identity
row: applied to the identity, a Π-reflection E forces E†(I) to be an eigenvector of L† (a left eigenvector of L)
at −c; with damping on the site the site equation fixes its factor there to Z + κ·I, κ = (γ↑ − γ↓)/(γ↑ + γ↓ + 2γ⊥)
= (q² − 1)/(q² + 1), and the product of these factors must commute with H, which across a hopping bond holds exactly
when the two κ, hence the two steady-state ratios, are equal. Together, for XXZ bonds with axial DM terms, without
fields or Z dephasing: **where the one-site steady-state ratio differs across a hopping bond between damped sites,
no site-local Π-reflection exists; where it is finite, nonzero and the same on every set of sites joined by hopping,
M_β with the steady state's weight is one,** carrying L onto −L − c with c = Σ(γ↓ + γ↑ + 2γ⊥). The bond half of the
sufficiency is KMS Question 4's proof with equal weights at both ends; without dephasing the steady state is the
bath's and the element is M_β as written. Checked: the site identity symbolically for every XY axis; the Heisenberg
chain with cooling and X dephasing at q² = 1/3, and a triangle with flux π/2, Δ = 7⁄10 and three unaligned axes
(3X + 4Y)/5, X, (5X − 12Y)/13 at q² = 1/4, both exactly on the full generator, q = 1 failing on each; raising one
heating rate on the triangle by 1⁄10 breaks its spectrum, Tr((L − μ)⁵) = 192. On the Heisenberg chain with a
temperature gradient the instrument's pre-check proves that no site-local Π-reflection exists, the identity row's
commutator being ‖[H, ⊗(Z_l + κ_l)]‖_F = 8√803/85, while the reflection −L† − c is found (S_XY). At a ratio of
zero or infinity (cooling or heating alone, without dephasing) the weight is singular, and the element is gone at
every bond: **a site with cooling alone and no field on it has a site-local Π-reflection while it is isolated, and
none as soon as any two-site term touches it, whatever the neighbour; heating alone is the same through X⊗N.** The
reason is in the heights. A product of site maps keeps the split of the generator by support, so the site equation
and each bond's equation hold apart. The cooled site's rates are 0 on |0⟩⟨0|, γ/2 on X and Y and γ on Z, so its
map must exchange the steady population and Z: E(|0⟩⟨0|) = (a/2)·Z, E(Z) = 2b·|0⟩⟨0|, X and Y among themselves.
The identity row and column are not what forbids it: they ask that |1…1⟩ be an eigenvector of H and that H commute
with Z⊗N, which a number-conserving bond such as the hopping satisfies; the bond's letters forbid it.
Write the bond as K_Z ⊗ Z + M ⊗ σ + M† ⊗ σ†, with σ = |0⟩⟨1| on the cooled site and K_Z, M traceless on the
neighbour (their trace parts would be a field on the cooled site). On A ⊗ |0⟩⟨0| the cooled site's populations see only the letter Z, and |0⟩⟨0| is an eigenletter of
multiplication by Z while Z is not: one side of the equation has a Z component there and the other an I component,
so a·[K_Z, E₁(A)] = 0 for every A, and K_Z = 0. On A ⊗ σ the populations give 2b·E₁(M†A) = [N, E₁(A)], with
N = β₁₁M† + β₂₁M and β the map on {σ, σ†}; for an invertible M the whole image of E₁ is then a commutator, hence
traceless, and E₁ is not invertible. A traceless M of rank one is |u⟩⟨v| with ⟨v|u⟩ = 0, which a unitary on the
neighbour, a turn about Z on the cooled site and a scale bring to the hopping |1⟩⟨0| ⊗ σ + h.c., and there an exact
solve with the neighbour's map left free has no invertible solution. So nothing invertible replaces M_β at ratio zero, whose limit there is singular (w → |0⟩⟨0|): at a
finite ratio a ZZ bond between two different ratios keeps the element (KMS Question 4's gate), at ratio zero not
even ZZ does. For the Hamiltonians F137 names the spectrum still pairs, through S_XY onto −L† − c and through its effective
Hamiltonian (a Z⊗X bond between two cooled sites does not pair); what
no product of site maps of any shape can do there is the Π-reflection. For a temperature gradient this settles the site-local half of OQ-168 and of the open
question in KMS Question 4; F137 scopes its sentence on the flux and the dephasing axes to S_XY, since for the
steady-state element the axes in the XY plane do not matter and the ratios do. Where fields break M_β, F137's
chained carrier gives an element of the same shape: S_XY, the thermal transpose, the conjugation and the turns W
about Z with W·(−H̄)·W† = H are each an identity between generators, so composed they are one site-local linear
map, ρ ↦ W·π_β^½·X⊗N·ρ·Y⊗N·π_β^½·W†, carrying L onto −L − c wherever the chain applies (π_β the bath's thermal
state). On F137's triangle with fields (3⁄10, −⅕, ½) at ratio ¼, Δ = 0 and flux ±π/2 (the hopping phase on the
bond (2, 0), W = Z₁) the instrument returns it as ρ ↦ π_β·Z₀Z₂·X⊗N·ρ·X⊗N·Z₁, the composite after the similarity
by Z⊗N·π_β^½, which commutes with L. Checked exactly there and on an open chain of four sites with hoppings 1,
3 + 4i and −5 + 12i; the same maps fail at Δ = 7⁄10 and at flux 0.

**The Z corner bounds how far from a product any element must sit.** Call Z⊗N the Z corner: the cube's vertex at
(0, N, N), every letter Z, the deepest cell under the baths. The site-local NONE across a temperature gradient has a
second route through it, and the route says more. Z⊗N is an eigenvector of L at −c for every H with [H, Z⊗N] = 0:
each site charges Z its height −(γ↓ + γ↑), the heights add, and the hopping and ZZ terms commute with Z⊗N (checked
entry for entry on the integer generator). An element E with E·L·E⁻¹ = −L − c sends an eigenvector at λ to one at
−λ − c, so E(Z⊗N) lies in the kernel of L, the steady state's line: the steady state is unique because σ⁻ and
σ⁺ on every bathed site have a trivial commutant. A map of operator Schmidt rank r across a cut turns a product
across that cut into an operator of rank at most r, and Z⊗N is a product across every cut. So across every cut
**Schmidt-rank(E) ≥ Schmidt-rank(ρ_ss)**: the element is at least as entangled as the state the bath holds the
chain in. At one temperature that state is the Gibbs product and the bound is 1, which M_β attains within its scope
(XXZ bonds with axial DM, no fields, a finite nonzero ratio). For XXZ bonds, two-body hopping, with a gradient along
a hopping bond the steady state is no product over the sites, and the reason is short. A product steady state is
diagonal site by site, by U(1) covariance. The dissipator keeps such a product diagonal and the hopping takes it
off the diagonal, so both parts vanish apart. Each factor is then the site's own thermal state, and the hopping
commutes with their product only at equal ratios along every bond. An invertible product map keeps the Schmidt
rank. So no product map over the sites is an element there, sandwich or not; the identity-row commutator above is
not needed for that conclusion. Across a single cut the bound is what the steady state's rank across that cut
says; that it is entangled across every cut, including a cut that crosses only an equal-ratio bond, is witnessed on
the chains below and not derived. What the letters show at N = 2 is a current, XY − YX on the bond. The ranks are
exact over ℚ (the steady state recovered from its residues modulo primes and verified by L·ρ_ss = 0, the rank by
elimination; sympy over ℚ(i) in the findings script, block G, and live in C# in `MirrorLocalityBoundWitness`,
`inspect --root locality`, claim `MirrorLocalityBoundClaim`): on the two-site
chain with γ↓ = (3, 2), γ↑ = (1, 1) and J = 1, ρ_ss = ¼·II + 5/56·IZ + 27/224·ZI + 67/1568·ZZ + 1/112·(XY − YX),
rank 4 of 4, the current ⟨i(σ⁺₀σ⁻₁ − σ⁻₀σ⁺₁)⟩ = −1/28; at N = 3 with γ↓ = (3, 2, 1), γ↑ = 1 and Δ = 1, rank 4 across
both cuts, and with γ↓ = (3, 3, 2), γ↑ = 1, rank 4 across the cut that crosses only the equal-ratio bond; at N = 4
with γ↓ = (6, 5, 4, 3), γ↑ = 2, rank 16 of 16 at the middle cut; at one temperature rank 1, and without the hop
rank 1. That some element exists is free wherever a reflection carries the pairing, as S_XY does for XXZ bonds at
any rates, since L is real in the Pauli basis and every matrix is similar to its transpose
([PROOF_PALINDROME_TWO_END_COUNT](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md) §(f8),
[THE_PAIRING_CONDITION](THE_PAIRING_CONDITION.md)); that settles the existence question of KMS Question 4 and
OQ-168 there, and leaves it open where the spectrum pairs and none of S_XY, M_β or the chain carries it. Whether an
entangled sandwich A·ρᵀ·B exists off the product line is not known; how far from a product an element must sit is
the next paragraph.

**At one two-site rate point no element that keeps the coherence order has rank 4, without ZZ or with it at Δ = 1, nor
one that reverses it, of which with ZZ there is none; one of rank 5 is found numerically.** Every element factors as E
= U·S_XY through F137's S_XY, ρ ↦ X⊗N·ρ·Y⊗N, which carries L onto −L† − c for XXZ bonds without fields or Z dephasing,
at any rates, so that U = E·S_XY⁻¹ solves L·U = U·L†. Since L commutes with ρ ↦ ρ†, L† = T·Lᵀ·T with T the transpose
in the computational basis, and Γ = U·T solves Lᵀ = Γ⁻¹·L·Γ, the relation of
[PROOF_ANTILINEAR_TRIANGLE](../docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md)'s proposition. Some such Γ exists for every
generator, every matrix being similar to its transpose; what S_XY adds is the correspondence E = Γ·T·S_XY, which keeps
the Schmidt rank. The proposition's argument still gives L·Γ(I) = Γ·Lᵀ(I) = 0, so Γ sends the identity onto ρ_ss. For
M_β at one ratio Γ is the thermal weight Γ_β, and Γ_β·T·S_XY is M_β up to the sign (−1)^N; the other elements there
have other Γ. On the two-site chain above, at its one rate point γ↓ = (3, 2), γ↑ = (1, 1), J = 1, and without ZZ, no
weight by ρ_ss serves: not from the left or from the right (exact), not split symmetrically as ρ_ss^½·X·ρ_ss^½
(numerical), each with the transpose or without it. How entangled must Γ, and with it E, be there? Take the coherence
order d = n_ket − n_bra of [THE_ONE_SQUARE](../docs/THE_ONE_SQUARE.md) §10, the excitation number's charge on the
cells, which the thermal jumps keep and S_XY reverses (F137), as does T, being ρ ↦ ρ† composed with complex
conjugation; so E keeps d exactly when Γ keeps it, and reverses d exactly when U keeps it. Write σ⁺ = |1⟩⟨0| and σ⁻ =
|0⟩⟨1| for the cells of order ±1 and n₀ = |0⟩⟨0|, n₁ = |1⟩⟨1| for the populations, and call the realigned matrix of a
two-site map the one indexed by its factors on the two sites, whose rank is the operator Schmidt rank across the cut.
For a map that keeps d the realigned matrix is block diagonal once the order shift δ of a factor on the first site, δ
= −2 … 2, is paired with the shift −δ on the second, with square blocks R_δ of sides 1, 4, 6, 4, 1, so its Schmidt
rank is the sum of the block ranks. The steady state splits the same way, 1 + 2 + 1: its populations have rank 2 at δ
= 0, and each half of the current, σ⁺⊗σ⁻ and σ⁻⊗σ⁺, has rank 1 at δ = ±1, with and without ZZ. E(Z⊗Z) ∝ ρ_ss for an E
that keeps d, and U(I) ∝ ρ_ss for the U of an E that reverses it (U(I) ∝ E(Z⊗Z), the Z corner again), puts that split
into the blocks, so rank 4 would need exactly R_{±2} = 0, rank R_{±1} = 1 and rank R₀ = 2. At that rate point the
system has no invertible solution: for the elements that keep d, without ZZ and with ZZ at Δ = 1 (the isotropic bond,
the one ZZ strength computed), and for those that reverse d, without ZZ (decided on their U). With ZZ no element
reverses d: at Δ = 1 the characteristic polynomial of L on the order-1 sector has no common factor with its complex
conjugate, so L·U = U·L† forces the order-±1 blocks of a U that keeps d to vanish (the U that keep d span 10
dimensions, all on the orders 0 and ±2, both exact), and every U·S_XY is singular; without ZZ that spectrum is closed
under conjugation. A Gröbner basis over ℚ(i) is {1} in every chart, with one Rabinowitsch variable per order sector
(z·det = 1, which forces that sector's block to be invertible) and charts over the two scalings that are symmetries,
the global scale and E ↦ ω^d·E for every nonzero complex ω. A third coordinate may not be normalised: a minor mixes
entries of different order sectors, so the system is not homogeneous sector by sector, and normalising a third
coordinate can fake an empty answer. Two controls: at one temperature, γ↓ = (3, 3), γ↑ = (1, 1) and c = 8, the exact
product element (on each site n₀ ↦ n₁, n₁ ↦ −3·n₀, σ⁺ ↦ σ⁺, σ⁻ ↦ −3·σ⁻) satisfies every equation with every sector
invertible; and at the gradient the equations of the sectors d ≥ 0 alone have solutions. The full system was not
decided on a case with solutions, which sympy did not finish. Rank 5 is reached without ZZ by an element that keeps d,
with block ranks 0, 1, 3, 1, 0. It is held to 80 digits, and its rank tails follow the digits kept, near 10⁻²⁰, 10⁻⁴⁰
and 10⁻⁸⁰, at a condition number of 10: a numerical point, not an exact one. So at that rate point the least rank
among the elements that keep d is 5 without ZZ (the 5 numerical) and at least 5 with ZZ; among those that reverse d it
is at least 5 without ZZ, and with ZZ there are none. For other ZZ strengths, for an element that mixes the orders, at
other rates and beyond two sites it is not known.

Why five, we cannot say in words. The cube was asked: one site and one bond give the hop's Jordan-Wigner step, which
pairs the four operators of order 1, σ⁺⊗I with Z⊗σ⁺ and σ⁺⊗Z with I⊗σ⁺, and the same at order −1; it does not say
where a fifth term sits, and the cube is silent there. What is known, in two classes. At the rank-5 element, which
keeps d, R₀ has rank 3, the rank of its population block (that block's singular values, scaled to its largest: 1,
0.31, 4·10⁻⁴, and a fourth that follows the digits kept, 4·10⁻⁸¹ at 80), so its generalized Schur complement vanishes
and the entries carrying the hop coherence |10⟩⟨01| = σ⁺⊗σ⁻, its mirror and the double coherences σ⁺⊗σ⁺ and σ⁻⊗σ⁻ onto
themselves are fixed by the population block and its cross entries. And without ZZ, on the maps U that keep d, which
belong to the elements that reverse d, the rank-4 equations of the sectors d ≥ 0 alone (the order-zero block, the
sector d = 1 and σ⁺⊗σ⁺) have, modulo the two scalings and with those sectors' blocks invertible, exactly six
solutions, the roots of three quadratics over ℚ(i), with a lex basis in shape position (every other variable a
polynomial in one); none of them extends to the full system, which is empty. The equations of the sectors d ≤ 0 have
the mirror solutions under ρ ↦ ρ†, which commutes with L, flips d and, applied on both sides, keeps a map linear and
its block ranks up to δ ↦ −δ. The computations of this and the previous paragraph are in
[`simulations/kms_element_least_rank.py`](../simulations/kms_element_least_rank.py): A exact, B exact over ℚ(i)
(its exclusion a second time in C#, each chart's Gröbner run over GF(p) replayed exactly over ℚ(i) to a nonzero
constant), C the rank-5 point as a law, D the six solutions exact and the population block at the rank-5 point as a law; the mirror and
the remark on normalising rest on the one-line arguments given with them. Keeping and reversing d do not exhaust the
elements: the solution spaces of the two classes share the maps of the order-zero sector onto itself and together span
28 of the 36 dimensions (20 of 28 with ZZ), the rest mixing the orders.

**2. The golden family routes with one map on every site at c = 0, and the golden bond with a uniform golden map.**
At the c = 0 station of the metallic line (windows XZY + YZX) the site map ρ ↦ (X + Y)·ρ·diag(1, −i), the same on
every site, solves the site and window equations; its identity column is the Pauli axis Y. So it routes every open
chain and every ring, and it fails at c = 3/2 and at the weights (0, 1, 2). The router proof shows uniform maps
empty for c ≠ 0 and names period-2 routers at c = 0; period 1 routes there as well. For the bond XX + XY + YX,
whose letter table [[1, 1], [1, 0]] is the golden window's, the uniform map ρ ↦ X·(g·ρ)ᵀ·X with g = φX + Y (φ the
golden ratio) is a Π-reflection on every graph whose bonds carry that table, and g = 2X + Y is not. The golden
frame is forced: applied to the identity, a uniform Π-reflection of that bond gives g⊗g, g = αX + βY its identity
column, which must commute with the bond, and it does only where α² − αβ − β² = 0; without the middle Z there is
nothing to force the period-4 rhythm. In letters φ is a half angle: the bond's σ⁺σ⁺ amplitude is 1 + 2i, the
element's diagonal factor turns by its phase arctan 2, and tan(½·arctan 2) = 1/φ. `PalindromeSoftCertifier` certifies the same bond through the chain's
reversal, a branch verified empirically and not derived ([SOFTNESS_IS_N_DEPENDENT](SOFTNESS_IS_N_DEPENDENT.md) holds
it to N = 6); the golden map is a product certificate, exact at every N. It is also a closed form, one map of
order four, for [CEILING_FOUR_NONLOCAL_CASES](CEILING_FOUR_NONLOCAL_CASES.md)'s XIX + XIY + YIX, whose next-nearest
bonds carry the same table and which the optimizer found routable; every uniform router of it, the optimizer's
included, has its identity column in the golden frame, so that local case meets the router family at c = 1.

**3. Heisenberg + DM is turned, not alternating.** Under Z dephasing a chain with Heisenberg bonds and a DM term
has no period-2 product of any invertible site maps that reflects it onto −L − c (checked at D/J = 1/2; the tied
equations of three sites reappear on every longer chain, so none exists at any N ≥ 3), while pure DM has one. What
reflects it is the canonical Π turned about Z by arctan(D/J) from site to site, the spiral gauge of the DM term
(F137's turns, θ_b = θ_a − arg(J + iD)); turned by twice that angle it fails. Read as a cell weight beside a move
that reverses turns, the instrument's output advances by twice the angle from site to site; the turn of the map is
arctan(D/J). For pure DM the angle is π/2, and since a turn by π only flips Π's sign, that is the period-2 map again.

## What stays open

Weights along X or Y (the same search in the Hadamard and quarter-turn frames); site relabellings; across a
temperature gradient, existence where the spectrum pairs with no reflection known and whether an entangled sandwich
A·ρᵀ·B exists off the product line; at the two-site rate point of finding 1, the least Schmidt rank of an element that
mixes the coherence orders, whether one that reverses them reaches 5 without ZZ, whether one that keeps them reaches 5
with ZZ, an exact element of rank 5 and why five in words; the least rank at ZZ strengths other than Δ = 1, at other
rates and beyond two sites; the per-site turn reported relative to the neighbour rather than as one representative's
absolute value; the instrument's name.
