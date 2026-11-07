# The palindrome without a carrier of its own

**Tom + Claude, 2026-11-06.** Cool every site of the XY triangle at its own rate, with no heating and no
dephasing, put a flux π/2 through it and a field on each site. [F137](../docs/ANALYTICAL_FORMULAS.md#f137)
says its spectrum is palindromic: at one common positive ratio of heating to cooling on all sites a local
element carries the pairing, and since the characteristic polynomial depends polynomially on that ratio, the
pairing survives at ratio zero. The element does not, because it contains the thermal weight π_β^½, which is
singular there. This page asks what carries the palindrome at ratio zero, and finds three things.

1. **It is always carried, onto the triangle with the flux reversed.** A site-local antilinear map takes the
   generator at flux π/2 to minus the adjoint of the generator at flux −π/2, shifted by c. Together with the
   triangular shape of the generator this gives the palindrome directly, exactly, at every point.
2. **Almost everywhere it is also carried onto itself**, by a local map that turns the damping's arrow around
   followed by two non-local dressings that take the arrow away and put the turned one in.
3. **At special points it is not.** There the triangle at flux π/2 and the same triangle at flux −π/2 have
   the same spectrum but opposite Jordan structure: one has a Jordan block at an eigenvalue λ where the other
   has two eigenvectors, and the other way round at the mirror point μ = −c − λ. Then no invertible map,
   linear or antilinear, takes L to −L − c or to −L† − c. This is read numerically, by an error law at a point
   located to ninety digits, not in exact arithmetic.

Where the palindrome has no carrier of its own, the Jordan structure feels the direction of the flux. Why the
points sit where they sit, we do not know.

<!-- STAGE 0 SWEEP RECORD. Asked: whether a palindromic spectrum under cooling alone (sigma minus jumps, A squared not one) is carried by a similarity, and what carries it. Stores swept by name. docs/ANALYTICAL_FORMULAS.md: F137 derives the cooling-alone palindrome spectrally (L block-triangular in the joint popcount, the spectrum that of the diagonal blocks through H_eff, without fields), states the Jordan structure only for S_XY where H commutes with X^N and Y^N, has the one-site identity S_XY D S_XY^-1 = -D^dagger - c that the arrow identity here reads, has S_IZ onto -L'^dagger - c with every jump replaced by its adjoint (a sibling, not the same object: here only the dissipator is adjointed), carries ratio zero by the polynomial dependence on the rates and says the element does not survive there, speaks of evenness under flux reversal, which the flux-reversed carrier here makes a map, and in its bullet "A unitary V whose conjugation commutes with the dissipator" reverses the flux by a site permutation that keeps rates and bonds and negates the fields, related to the relabelling here; its bullet "The three chained, with fields" says the palindrome holds at cooling alone through the characteristic polynomial, the carrier being gone (docs/CAUGHT_ERRORS.md 2026-11-06, F137's chained carrier). F158's fence on T1 jumps points here and states the Jordan lift only for its own class. docs/THE_ONE_SQUARE.md section 9, the subsections "Which rows survive amplitude damping" and "The element S_XY and the thermal palindrome": the arrow I to Z, and a reflection has to turn it round (gate M9). docs/KMS_DETAILED_BALANCE.md: the singular weight of pure cooling, the H_eff route, nothing on Jordan structure. docs/proofs: PROOF_PALINDROME_TWO_END_COUNT section (f8) lifts multiset pairing to a similarity for A squared one only, and its Outside paragraph says nothing on the Jordan structure for A squared not one; PROOF_ANTILINEAR_TRIANGLE section 7 records that a cold bath has no thermal weight; PROOF_RING_FROZEN_BAND Lemma 1 gives the wrap sign s = (-1)^(p-1) used below; PROOF_COLLISION_GAP_ODD_ORDERS section (i) reads the two Jordan-Wigner parity sectors as half a flux quantum apart, without cooling and without an external flux; no proof states a Jordan structure under T1. experiments: THE_CUBE_SITE_BY_SITE finding 1 proves that no invertible product of site maps reflects onto -L - c at ratio zero once a bond touches a site with cooling alone and no field (whether its proof extends to a field on that site is not checked here); WHICH_BONDS_PAIR_UNDER_COOLING power sums only; THE_PAIRING_CONDITION the F158 similarity under dephasing only; PT_SYMMETRY_ANALYSIS Jordan blocks at dephasing exceptional points, unrelated; nothing else near Jordan and cooling. OpenArcs registry: nothing on T1 pairing, Jordan structure or a singular weight; f138_converse_failures is about dephasing. Typed layer: F137 untyped; PalindromeTwoEndCountClaim names the similarity as unobserved by its gates; MirrorLocalityBoundClaim says existence stays open where the spectrum pairs and no reflection is known, which this page answers at its points; CooledBondPowerSumsClaim and the Diagnostics witnesses ThermalElementLeastRankWitness and CooledBondPowerSumsWitness carry nothing on Jordan structure or ratio zero. Glossary: Jordan-Wigner and the wrap sign, nothing on Jordan structure under cooling. CAUGHT_ERRORS: the 2026-08-30 third entry, a criterion for a Jordan block built from the characteristic polynomial, which cannot see one, is the same blindness this page meets in F137's polynomial argument. Confirmations: no T1 pairing entry. Adjacency checked: F137's flux condition e^(2i Phi) = (-1)^l and PROOF_RING_FROZEN_BAND Lemma 1 with PROOF_COLLISION_GAP_ODD_ORDERS section (i) are about the same parity sectors and do not cite each other, joined below only as a reading; F137's polynomial limit and the CAUGHT_ERRORS 2026-08-30 entry are joined here; F137's cooling paragraphs and F158 section (f8) are joined through F158's fence; the resonance condition on sums of damped single-particle energies is a level-collision question like F129's (MirrorWorld LevelCollision), there on the undamped cosine comb of the chain, adjacent and not joined; P's action on the blocks, (p, q) to (N - p, N - q), is the composite of the two folds of MirrorWorld's Mirror, joined in the text. -->

## The setting

Three sites on a ring. The Hamiltonian is J(XX + YY) on the bonds (0, 1) and (1, 2), J(XY − YX) on the seam
(2, 0), which puts a flux π/2 through the triangle (−J(XY − YX) puts −π/2), and fields Σ h_l Z_l. Each site is
cooled at its own rate γ_l by the jump σ⁻ = |0⟩⟨1|, with Z|0⟩ = |0⟩. Write the generator as

    L = L_H + D,   D(ρ) = Σ_l γ_l (σ⁻_l ρ σ⁺_l − ½{n_l, ρ}),   n_l = σ⁺_l σ⁻_l,

and L(Φ) where the flux matters. The term σ⁻ρσ⁺ lowers ket and bra together, so L is block-triangular in the
joint popcount (p, q), the grading of MirrorWorld's [`Block`](../compute/MirrorWorld/Block.cs): every other
term keeps (p, q), and this one, the **feeding**, steps (p, q) → (p − 1, q − 1). The blocks with one p − q form
a line, and a **level** is one block of a line together with its eigenvalues. Call the block-diagonal part Λ
and the feeding F↓, so L = Λ + F↓. The spectrum of L is the spectrum of Λ, and it pairs about −Σγ_l/2, so
c = Σγ_l.

On the letter cube this is the reading of σ⁻ in
[THE_ONE_SQUARE](../docs/THE_ONE_SQUARE.md#9-the-moves-and-the-two-questions) §9 ("Which rows survive amplitude damping" and "The element S_XY and the thermal palindrome"): two transverse lights, a height, plus one **arrow** I → Z per site, and a reflection has to turn
that arrow round. The heights can be mirrored; the question of this page is the arrow.

## The local map turns the arrow around

The Hilbert-Schmidt adjoint of the cooling dissipator is

    D†(ρ) = Σ_l γ_l (σ⁺_l ρ σ⁻_l − ½{n_l, ρ}),

the same dissipator with its feeding turned around: the anticommutator, and with it every diagonal block, is
unchanged, and the arrow now climbs, (p, q) → (p + 1, q + 1). Call **L̃ = L_H + D† = Λ + F↑** the **turned
generator**, F↑ = F↓ᵀ (σ⁻ is real); it is L(−H)†, exactly (gate Part A).

F137's element S_XY: ρ ↦ X⊗N·ρ·Y⊗N satisfies the one-site identity S_XY·D·S_XY⁻¹ = −D† − c (THE_ONE_SQUARE
§9, gate M9). For a Hamiltonian that S_XY keeps, XX + YY bonds without fields, L_H = −L_H† makes S_XY a
reflection of the whole L onto −L† − c, the case F137 and THE_ONE_SQUARE gate. The triangle here is not such a
Hamiltonian: X⊗N negates the fields and the seam's XY − YX. A turn about Z on site 1 gives a map that negates
every term of H, so it treats the Hamiltonian part like Π (L_H ↦ −L_H) and the dissipator like the reflection
(D ↦ −D† − c):

    P·L·P⁻¹ = −L̃ − c,   P(ρ) = Z₁X⊗3·ρ·Y⊗3·Z₁.

F137's flux condition e^{2iΦ} = (−1)^ℓ, for a cycle of length ℓ with flux Φ, is what allows such a turn; on the
triangle ℓ = 3 and Φ = ±π/2. In block terms the reason the arrow turns is visible: X⊗N on the ket and Y⊗N on
the bra send the block (p, q) to (N − p, N − q), the composite of the two folds of MirrorWorld's
[`Mirror`](../compute/MirrorWorld/Mirror.cs), so a step down, (p, q) → (p − 1, q − 1), lands as a step up,
(N − p, N − q) → (N − p + 1, N − q + 1). The identity is checked exactly in integers at a rational point with
unequal fields and rates, and no such sandwich exists at flux 0 on the odd triangle (gate Part A).

The identity also proves the palindrome directly, at every point and every rate: L and L̃ have the same
diagonal blocks, so the same spectrum, and P gives spec L = −c − spec L̃ = −c − spec L.

It also gives a reading, not a proof, of why F137's thermal element is local at a positive ratio κ = γ↑/γ↓ and
the Π-shaped one is not at κ = 0. The weight π_β^½ multiplies each block (p, q) by κ^((p+q)/2); conjugating the
feeding by it multiplies the downhill feeding by κ and divides the uphill one by κ, which balances the two
arrows. At κ = 0 there is nothing left to balance. For a site with cooling alone and no field,
[THE_CUBE_SITE_BY_SITE](THE_CUBE_SITE_BY_SITE.md) finding 1 proves that no product of site maps reflects onto
−L − c once a bond touches it; whether that proof extends to a field on the cooled site is not checked here.

## The flux-reversed triangle

Complex conjugation of the state sends L(−H) to L(H̄), and H̄ is the triangle at the reversed flux, with the
same fields. With L̃ = L(−H)† this makes the site-local antilinear map

    ρ ↦ conj(Z₁X⊗3·ρ·Y⊗3·Z₁)   carry L(π/2) onto −L(−π/2)† − c,

exactly in integers at the rational point of Part A (gate Part A). So the palindrome is always carried, by a
local map, from the triangle onto the same triangle with its flux reversed. Relabelling sites 0 and 2 also
reverses the flux, so L(−π/2) is, after relabelling, exactly the triangle at flux π/2 with h₀ ↔ h₂ and
γ₀ ↔ γ₂, its **mirror twin** (gate Part A); the search and Part C below use the twin, because it keeps the seam where it is. On the
symmetric locus h₀ = h₂, γ₀ = γ₂ the relabelling maps the triangle onto its own flux reversal, the two are
similar, and the palindrome is carried onto the triangle itself. F137's bullet "A unitary V whose
conjugation commutes with the dissipator" reverses the flux by a site permutation too, there with rates the
permutation keeps and fields it carries into their negatives; a related case, not this one.

The two lines L̃ = L(−H)† and K·L(−H)·K = L(H̄) together give L̃(π/2) = L(−π/2)ᵀ, exactly (gate Part A): the
generator with the arrow turned is the transpose of the generator with the flux turned. A transpose keeps the
spectrum and the Jordan structure, so off the symmetric locus the question of a carrier of its own becomes:
is L similar to L̃, that is, do the triangle at π/2 and at −π/2 have the same Jordan structure?

## Why the Jordan structure usually follows

L = Λ + F↓ and L̃ = Λ + F↑ have the same diagonal blocks; one feeds down the lines p − q = const, the other up.
The lines do not talk to each other, so each line can be taken alone. If no two blocks on a line share an
eigenvalue, each of L and L̃ is similar to Λ by a unipotent block-triangular matrix (the Sylvester equations for
the off-diagonal blocks are solvable), so L ~ Λ ~ L̃ ~ −L − c: a carrier of its own, made of P and two
non-local unipotent dressings, one taking the downhill arrow away and one putting the uphill arrow in.

At a rational point with fields (1, −2, 3), rates (2, 5, 7) and J = 1, all fourteen pairs of blocks on a common
line have coprime characteristic polynomials, exactly over ℚ(i) (gate Part B). Coprimality is a resultant
that does not vanish, and a polynomial in the parameters that is nonzero at one point is nonzero on a dense
open set, so generic points carry their palindrome with the whole Jordan structure. Exceptional points
inside the levels do not change this by themselves. At rates (8, 8 ± 4√3) with a uniform field the
single-particle block has one triple eigenvalue with one eigenvector, and at fields (0, ½, ½) with rates
(6, 6 ± 2√7) one double eigenvalue with one eigenvector (both exact, gate Part E). At both points L has twelve
defective repeated factors modulo each of two primes in which i and the square root exist, and for every
repeated factor the kernel dimensions of its powers match those of its mirror factor, so the Jordan
structure is mirrored modulo both primes; at flux 0 the rates of either point give no palindrome at all
(gate Part E).

## Where two levels can meet

The XY ring under cooling is quadratic in Jordan-Wigner fermions. By
[PROOF_RING_FROZEN_BAND](../docs/proofs/PROOF_RING_FROZEN_BAND.md) Lemma 1 the sector of p particles is the
ring with the fixed sign s = (−1)^(p−1), so with the seam's flux Φ it sees the single-particle matrix

    K_p = h(Φ + (p − 1)π) + diag(−2h_l − iγ_l/2),

h the hopping matrix (entries 2J, the seam's carrying the phase), and the block (p, q) has the eigenvalues
−i(E_A − Ē_B), E_A a sum of p eigenvalues of K_p and E_B of q eigenvalues of K_q (gate Part D, in floating
point). At Φ = π/2 on the triangle the two parities see Φ = π/2 and Φ = −π/2, and K(−Φ) = K(Φ)ᵀ, so both
sectors have the same three single-particle energies a₁, a₂, a₃ (gate Part D): on this page F137's odd-ring
condition reads as the two fermion parities being transposed images of each other. This joins F137's flux
condition to Lemma 1 and to [PROOF_COLLISION_GAP_ODD_ORDERS](../docs/proofs/PROOF_COLLISION_GAP_ODD_ORDERS.md)
§(i) as a reading only.

The imaginary parts decide which levels can meet; this list is argued from them and not gated. The real part
of an eigenvalue of block (p, q) is the sum of the imaginary parts of the p + q single-particle energies
involved, each negative, and the three sum to −Σγ_l/2. The lines with p − q < 0 are the conjugates of those
with p − q > 0 (ρ ↦ ρ†), so the positive lines and the diagonal say everything.

- **Equal rates never resonate.** Then every single-particle energy has imaginary part −γ/2, and the real
  parts of neighbouring levels on a line differ by γ.
- **Levels two steps apart never meet, and neither do (3, 3) and (2, 2).** (3, 2) against (1, 0): −Σγ/2 plus
  two more negative terms against a single one above −Σγ/2. (2, 2) and (3, 3) against (0, 0): a strictly
  negative real part against zero. (3, 3) against (1, 1): −Σγ against a sum of two single-particle imaginary
  parts, above −Σγ. (3, 3) against (2, 2): −Σγ against −Σγ minus two negative terms.
- **(1, 1) with (0, 0) and (3, 1) with (2, 0) never meet**: both would need a single-particle imaginary part
  of the wrong sign.
- **(2, 1) with (1, 0)** can meet, with unequal rates, and so can its image (3, 2) with (2, 1) under P followed
  by ρ ↦ ρ†. A resonance asks a_i + a_j − ā_k = a_m. Its imaginary part is Im a_i + Im a_j + Im a_k = Im a_m,
  so m is not i or j (else two negative parts would sum to zero), hence m is the third index; and k is i or j
  (else Im a_i + Im a_j = 0). With k = i the condition is Re a_j = Re a_m and Im a_m = 2 Im a_i + Im a_j. Two
  real conditions; the points of this page have J of about 0.14 and 0.21 against rates between 0.4 and 4.2.
- **(2, 2) with (1, 1)** is not excluded by the imaginary parts; this page does not follow that line.

## The point where nothing carries it onto itself

At a resonance λ between a level A = (p, q) and the level B = (p − 1, q − 1) below it, with λ simple in each
and absent from the other blocks on that line, L has a Jordan chain at λ exactly when the feeding couples A's
right eigenvector to B's left one, η↓ = ⟨ℓ_B|F↓|r_A⟩ ≠ 0, and L̃ exactly when η↑ = ⟨ℓ_A|F↑|r_B⟩ ≠ 0 (the
script's c_down and c_up). If one vanishes and the other does not, L and L̃ differ at λ, and through P, L and
−L − c differ at the mirror point μ = −c − λ. Each coupling depends on how the eigenvectors are scaled; their
product does not, and the solver works with the product.

Reversing the flux exchanges the two couplings. Since L̃(π/2) = L(−π/2)ᵀ, and a transpose keeps the eigenvalues
and exchanges left and right eigenvectors, η↑ of the triangle at λ has the size of η↓ of its flux reversal at
the same λ, and the other way round. This is an identity; Part C checks one instance, on the twin. So the points come in
pairs, a triangle with η↑ = 0 and its flux reversal with η↓ = 0, with the same spectrum and opposite Jordan
structure. On the symmetric locus the two couplings are equal in size, and no such point sits there.

With seven parameters (three fields, three rates, J), the resonance (2, 1) ↔ (1, 0) and η↓·η↑ = 0 are four
real conditions. A seeded search from 120 random starts converged at 14 points with η↑ = 0 and at 4 with
η↓ = 0, and at every one of them the twin has the other coupling vanishing (`--search`); the split 14 to 4
reflects where the search started, not the physics. At one point with η↑ = 0, three parameters held and four
solved by Newton, the Jacobian of the four conditions is nonsingular (smallest singular value 0.014 against a
largest of 7.5), so near it they leave a smooth set of dimension three in the seven parameters, and at ninety
digits

    fields (−0.2569, −1.1562, −0.2672),  rates (3.6080, 0.4112, 1.5775),  J = 0.2086,

with a resonance gap of 6·10⁻⁹¹, |η↑| = 10⁻⁹⁰ and |η↓| = 0.095 (normalised). λ is simple in both blocks and
absent from the third block (3, 2) on the line, each by more than 0.1; the whole spectrum pairs about
−Σγ_l/2 = −2.7983, and μ is a double eigenvalue like λ (gate Part C). Whether L has a Jordan block at a double
eigenvalue is read from how a fixed random perturbation εE splits it, by a law and not by a threshold: a
Jordan block of size two splits like ε^½, a semisimple double eigenvalue like ε.

| ε | split at λ = −1.5378 − 0.4535i | split at μ = −4.0588 + 0.4535i |
|---|---|---|
| 10⁻²⁰ | 1.037·10⁻¹⁰ | 6.766·10⁻²⁰ |
| 10⁻²⁸ | 1.037·10⁻¹⁴ | 6.766·10⁻²⁸ |
| 10⁻³⁶ | 1.037·10⁻¹⁸ | 6.766·10⁻³⁶ |
| 10⁻⁴⁴ | 1.037·10⁻²² | 6.766·10⁻⁴⁴ |

The ratios to ε^½ at λ and to ε at μ are constant over twenty-four decades: L has a Jordan block at λ and two
eigenvectors at μ. As a control, at a resonance where both couplings are nonzero (|η↓| = 0.056 and
|η↑| = 0.031), both sides split like ε^½ with constant ratios, so the law tells the two cases apart (gate
Part C). That control's λ ≈ −2 lies 2·10⁻⁷ off the real axis, so its conjugate is 5·10⁻⁷ away; the splits,
at most 4·10⁻¹¹, stay far below that.

So at the exact point these digits approximate, L is not similar to −L − c: no invertible linear S has
S·L·S⁻¹ = −L − c. Nor does F137's reflection target help, since −L† − c is similar to −L − c (L is real in the
Hermitian Pauli basis, so L† ~ Lᵀ ~ L). No antilinear S does either: an antilinear carrier would match L's
structure at λ with its structure at −c − λ̄ = μ̄, and L commutes with ρ ↦ ρ†, so its structure at μ̄ is its
structure at μ, which is semisimple. What still carries the palindrome there is the flux-reversed map.

Read as dynamics: at flux π/2 a state with a component along the Jordan chain at λ decays partly as
t·e^{λt}, while at μ every state decays as a pure exponential; at flux −π/2 it is the other way round. The
spectrum pairs, the time dependence does not, and reversing the flux exchanges the two ends.

## What we do not understand

- **Why there.** We can say when an arrow coupling vanishes, not why: there is no reading yet, in the cube's
  words or the fermions', of what makes η↑ vanish at a resonance, beyond its exchange with η↓ under flux
  reversal.
- **Beyond this triangle.** Whether such points exist on other rings, at other N, with dephasing added, is
  open, and so is whether the flux-reversed carrier has a counterpart there.
- **An exact point.** The evidence for the opposite Jordan structure is the error law at a point located to
  ninety digits, not exact arithmetic. The exact algebraic point and an exact Jordan computation there are
  not done.
- **What it means.** A palindrome whose carrier always leads to the flux-reversed system and only sometimes
  back to itself contradicts no statement in the repository: F158's lift to a similarity is for jumps with
  A² = 1, and its fence says so. What it means that the Jordan structure can feel the direction of the flux
  while the spectrum does not, we do not know yet.

## Reproduce

```
python simulations/cooled_ring_arrow_mirror.py            # the gates, about two minutes
python simulations/cooled_ring_arrow_mirror.py --search   # the seeded search, about two and a half minutes
```

numpy, scipy, sympy and mpmath. Part A is exact in integers, Part B exact over ℚ(i), Part C the error law at
ninety digits with its control, Part D the free-fermion reading in floating point, Part E exact over ℚ(√3, i)
and ℚ(√7, i) for the single-particle blocks and modulo two primes for L.
