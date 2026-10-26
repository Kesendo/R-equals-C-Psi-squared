# The Cube, Site by Site

<!-- Keywords: letter cube instrument, site-local palindromizer, product reflector, steady-state weight, M_beta,
KMS detailed balance, temperature gradient, golden router uniform, Heisenberg DM spiral gauge, Groebner exact -->

**Status:** an instrument and three findings, 2026-10-20. The findings are reproduced exactly, each with a control
that fails, in [`simulations/cube_site_by_site_findings.py`](../simulations/cube_site_by_site_findings.py) (44
checks, ALL PASS); where a statement below rests on an argument or on a computation rather than on a check, it
says so.

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
-->

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

Weights along X or Y (the same search in the Hadamard and quarter-turn frames); site relabellings; elements that are
not products over sites; the per-site turn reported
relative to the neighbour rather than as one representative's absolute value; the instrument's name.
