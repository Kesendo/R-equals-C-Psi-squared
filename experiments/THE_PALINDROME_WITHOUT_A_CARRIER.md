# The palindrome without a carrier

**Tom + Claude, 2026-11-06.** Under cooling alone, amplitude damping on every site with no heating and no
dephasing, the spectrum of an XY ring with longitudinal fields stays palindromic where
[F137](../docs/ANALYTICAL_FORMULAS.md#f137) says it does: the pairing holds at every positive ratio of heating
to cooling and is polynomial in the rates, so it holds at ratio zero too. The site-local element that carries
it at a positive ratio has a weight that is singular at ratio zero. This page asks what is left on the
smallest odd ring, the triangle at flux π/2: is there still a matrix S with S·L·S⁻¹ = −L − c, a carrier of any
shape, or only the eigenvalues? Almost everywhere a carrier exists, and it comes in two parts: a local mirror
that turns the dissipator's arrow around, and a non-local dressing that undoes the turned arrow. At special
points, where two levels the arrow links share an eigenvalue and the arrow chains them in L but the turned
arrow does not chain them in the mirror image, nothing carries it: the spectrum pairs, the Jordan structure
does not, and no linear or antilinear map takes L to −L − c. We located such points numerically and read the
Jordan structure at one of them by an error law; we do not yet understand why they sit where they sit.

<!-- STAGE 0 SWEEP RECORD. Asked: whether a palindromic spectrum under cooling alone (sigma minus jumps, A squared not one) is carried by a similarity, and what carries it. Stores swept by name. docs/ANALYTICAL_FORMULAS.md: F137 derives the cooling-alone palindrome spectrally (L block-triangular in the joint popcount, the spectrum that of the diagonal blocks through H_eff, around line 6740), states the Jordan structure only for S_XY where H commutes with X^N and Y^N (around line 6793), has the one-site identity S_XY D S_XY^-1 = -D^dagger - c (around line 6790) that the arrow identity here reads, has S_IZ onto -L'^dagger - c with every jump replaced by its adjoint (a sibling, not the same object: here only the dissipator is adjointed), carries ratio zero by the polynomial dependence on the rates and says the element does not survive there (lines 6878 to 6883), and its chained-carrier bullet says the palindrome holds at cooling alone through the characteristic polynomial, the carrier being gone (docs/CAUGHT_ERRORS.md 2026-11-06 records the two sentences that once read otherwise). F158's fence sends T1 jumps to F137 and states the Jordan lift only for its own class. docs/THE_ONE_SQUARE.md section 9, the amplitude-damping paragraph: the arrow I to Z, and a reflection has to turn it round (gate M9), which is the reading used here. docs/proofs: PROOF_PALINDROME_TWO_END_COUNT section (f8) lifts multiset pairing to a similarity for A squared one only, and its Outside paragraph says the criterion says nothing on the Jordan structure for A squared not one; PROOF_ANTILINEAR_TRIANGLE section 7 records that a cold bath has no thermal weight; PROOF_RING_FROZEN_BAND Lemma 1 gives the wrap sign s = (-1)^(p-1) used below; PROOF_COLLISION_GAP_ODD_ORDERS section (i) reads the two Jordan-Wigner parity sectors as half a flux quantum apart, without cooling and without an external flux; no proof states a Jordan structure under T1. experiments: THE_CUBE_SITE_BY_SITE finding 1 proves that no invertible product of site maps reflects onto -L - c at ratio zero once a bond touches a site with cooling alone and no field; this page has fields on every site, where that is not proved; KMS_DETAILED_BALANCE (docs) has the singular weight of pure cooling; WHICH_BONDS_PAIR_UNDER_COOLING power sums only; THE_PAIRING_CONDITION the F158 similarity under dephasing only; PT_SYMMETRY_ANALYSIS Jordan blocks at dephasing exceptional points, unrelated; nothing else near Jordan and cooling. OpenArcs registry: nothing on T1 pairing, Jordan structure or a singular weight; f138_converse_failures is about dephasing. Typed layer: F137 untyped; PalindromeTwoEndCountClaim names the similarity as unobserved by its gates; MirrorLocalityBoundClaim says existence stays open where the spectrum pairs and no reflection is known, which this page answers negatively at its points; CooledBondPowerSumsClaim and the Diagnostics witnesses ThermalElementLeastRankWitness and CooledBondPowerSumsWitness carry nothing on Jordan structure or ratio zero. Glossary: Jordan-Wigner and the wrap sign, nothing on Jordan structure under cooling. CAUGHT_ERRORS: the 2026-08-30 third entry, a criterion for a Jordan block built from the characteristic polynomial, which cannot see one, is the same blindness this page meets in F137's polynomial argument. Confirmations: no T1 pairing entry. Adjacency checked: F137's flux condition e^(2i Phi) = (-1)^l and PROOF_RING_FROZEN_BAND Lemma 1 with PROOF_COLLISION_GAP_ODD_ORDERS section (i) are about the same parity sectors and do not cite each other, joined below only as a reading; F137's polynomial limit and the CAUGHT_ERRORS 2026-08-30 entry are joined here; F137's cooling paragraphs and F158 section (f8) are joined through F158's fence. -->

## The setting

Three sites on a ring. The Hamiltonian is J(XX + YY) on the bonds (0, 1) and (1, 2), J(XY − YX) on the seam
(2, 0), which puts a flux π/2 through the triangle, and fields Σ h_l Z_l. Each site is cooled at its own rate
γ_l by the jump σ⁻ = |0⟩⟨1|, with Z|0⟩ = |0⟩. Write the generator as

    L = L_H + D,   D(ρ) = Σ_l γ_l (σ⁻_l ρ σ⁺_l − ½{n_l, ρ}),   n_l = σ⁺_l σ⁻_l.

The term σ⁻ρσ⁺ lowers ket and bra together, so L is block-triangular in the joint popcount (p, q): every
other term keeps (p, q), and this one, the **feeding**, steps (p, q) → (p − 1, q − 1). Call the block-diagonal
part Λ and the feeding F↓, so L = Λ + F↓. The spectrum is the union of the diagonal blocks' spectra, and F137
pairs it about −Σγ_l/2, so c = Σγ_l.

On the letter cube this is the reading of σ⁻ in
[THE_ONE_SQUARE](../docs/THE_ONE_SQUARE.md#9-the-moves-and-the-two-questions) §9 ("Which rows survive amplitude
damping"): two transverse lights, a height, plus one **arrow** I → Z per site, and a reflection has to turn
that arrow round. The heights can be mirrored; the question of this page is the arrow.

## The mirror turns the arrow around

F137's element S_XY: ρ ↦ X⊗N·ρ·Y⊗N satisfies the one-site identity S_XY·D·S_XY⁻¹ = −D† − c (THE_ONE_SQUARE
§9, gate M9). The Hilbert-Schmidt adjoint of the cooling dissipator is

    D†(ρ) = Σ_l γ_l (σ⁺_l ρ σ⁻_l − ½{n_l, ρ}),

the same dissipator with its feeding turned around: the anticommutator, and with it every diagonal block, is
unchanged, and the arrow now climbs, (p, q) → (p + 1, q + 1). Write L̃ = L_H + D† = Λ + F↑ with F↑ = F↓ᵀ
(σ⁻ is real). For a Hamiltonian that S_XY keeps, XX + YY bonds without fields, L_H = −L_H† makes S_XY a
reflection of the whole L onto −L† − c, the case F137 and THE_ONE_SQUARE gate. The triangle here is not such a
Hamiltonian: X⊗N negates the fields and the seam's XY − YX. A turn about Z on site 1, which F137's flux
condition e^{2iΦ} = (−1)^ℓ allows on the triangle at flux π/2, gives a map that negates every term of H, so
it treats the Hamiltonian part like Π (L_H ↦ −L_H) and the dissipator like the reflection (D ↦ −D† − c):

    P·L·P⁻¹ = −L̃ − c,   P(ρ) = Z₁X⊗3·ρ·Y⊗3·Z₁.

So the local mirror does reach a target, half Π and half reflection, and the target is the generator with the
arrow reversed. The identity is checked exactly in integers at a rational point with unequal fields and rates,
and no such sandwich exists at flux 0 on the odd triangle (gate Part A).

It also gives a reading, not a proof, of why the thermal element of F137 is local at a positive ratio κ =
γ↑/γ↓ and nothing local is known at κ = 0. The element's weight π_β^½ multiplies each block (p, q) by
κ^((p+q)/2); moved across the identity, it multiplies the downhill feeding by κ and divides the uphill one by
κ, which balances the two arrows. At κ = 0 there is nothing left to balance. For a site with cooling alone and
no field, [THE_CUBE_SITE_BY_SITE](THE_CUBE_SITE_BY_SITE.md) finding 1 proves that no product of site maps
reflects onto −L − c once a bond touches it; with fields on every site, as here, that is not proved, and at
the points below no map of any shape does.

The whole question of a carrier is now whether L and L̃ are similar.

## Why the Jordan structure usually follows

L = Λ + F↓ and L̃ = Λ + F↑ have the same diagonal blocks; one feeds down the lines p − q = const, the other up.
The lines do not talk to each other, so each line can be taken alone. If no two blocks on a line share an
eigenvalue, each of L and L̃ is similar to Λ by a unipotent block-triangular matrix (the Sylvester equations for
the off-diagonal blocks are solvable), so L ~ Λ ~ L̃ ~ −L − c: a carrier made of the local mirror P and two
non-local unipotent dressings, one taking the downhill arrow away and one putting the uphill arrow in.

At a rational point with fields (1, −2, 3), rates (2, 5, 7) and J = 1, all fourteen pairs of blocks on a common
line have coprime characteristic polynomials, exactly over ℚ(i) (gate Part B). Coprimality is a resultant
that does not vanish, and a polynomial in the parameters that is nonzero at one point is nonzero on a dense
open set, so generic points are carried with their whole Jordan structure. Exceptional points do not change
this by themselves: at a triple exceptional point of the single-particle block (rates (8, 8 ± 4√3), uniform
field) and at an exceptional point with unequal fields (0, ½, ½) and rates (6, 6 ± 2√7), L has twelve
defective repeated factors, Jordan blocks inside its levels, and every one is mirrored, exactly modulo two
primes in which i and the square root exist (gate Part E). Those blocks sit in Λ, which L and L̃ share.

## Where two levels can meet

The XY ring under cooling is quadratic in Jordan-Wigner fermions. By
[PROOF_RING_FROZEN_BAND](../docs/proofs/PROOF_RING_FROZEN_BAND.md) Lemma 1 the sector of p particles is the
ring with the fixed sign s = (−1)^(p−1), so with the seam's flux Φ it sees the single-particle matrix

    K_p = h(Φ + (p − 1)π) + diag(−2h_l − iγ_l/2),

h the hopping matrix (entries 2J, the seam's carrying the phase), and the block (p, q) has the eigenvalues
−i(E_A − Ē_B), E_A a sum of p eigenvalues of K_p and E_B of q eigenvalues of K_q (gate Part D, in floating
point). At Φ = π/2 on the triangle the two parities see Φ = π/2 and Φ = −π/2, and K(−Φ) = K(Φ)ᵀ, so both
sectors have the same single-particle spectrum: on this page F137's odd-ring condition reads as the two
fermion parities being transposed images of each other. This joins F137's flux condition to Lemma 1 and to
[PROOF_COLLISION_GAP_ODD_ORDERS](../docs/proofs/PROOF_COLLISION_GAP_ODD_ORDERS.md) §(i) as a reading only.

The imaginary parts decide which levels can meet. The real part of an eigenvalue of block (p, q) is the sum
of the imaginary parts of the p + q single-particle energies involved, each negative, and the three of a
sector sum to −Σγ_l/2.

- **Equal rates never resonate.** Then every single-particle energy has imaginary part −γ/2, and the real
  parts of neighbouring levels on a line differ by γ.
- **Levels two steps apart never meet.** (3, 2) against (1, 0): −Σγ/2 plus two more negative terms against a
  single one above −Σγ/2. (2, 2) and (3, 3) against (0, 0): a strictly negative real part against zero.
  (3, 3) against (1, 1): −Σγ against something above −Σγ/2.
- **(1, 1) with (0, 0) and (3, 1) with (2, 0) never meet**: both would need a single-particle imaginary part
  of the wrong sign.
- **(2, 1) with (1, 0)**, and its mirror (3, 2) with (2, 1), can meet, with unequal rates: for some labelling
  a_l, a_x, a_y of the single-particle energies, Re a_l = Re a_y and Im a_l = 2 Im a_x + Im a_y, and no other
  labelling works. Two real conditions; the points of this page have J of about 0.14 and 0.21 against rates
  between 0.4 and 4.2.
- **(2, 2) with (1, 1)** is not excluded by the imaginary parts; this page does not follow that line.

## The point where nothing carries it

At a resonance λ between a level A = (p, q) and the level B = (p − 1, q − 1) below it, with λ simple in each
and absent from the other blocks on that line, L has a Jordan chain at λ exactly when the feeding couples A's
right eigenvector to B's left one, c↓ = ⟨ℓ_B|F↓|r_A⟩ ≠ 0, and L̃ exactly when c↑ = ⟨ℓ_A|F↑|r_B⟩ ≠ 0. If one
vanishes and the other does not, L and L̃ differ at λ, and through P, L and −L − c differ at the mirror point
μ = −c − λ. Each coupling depends on how the eigenvectors are scaled; their product does not, and the solver
works with the product.

With seven parameters (three fields, three rates, J), the resonance (2, 1) ↔ (1, 0) and c↓·c↑ = 0 are four
real conditions. A seeded search from 120 random starts converged at 14 points with c↑ = 0 and at 4 with
c↓ = 0 (`--search`), so both arrows can fall silent; when c↓ = 0 the roles exchange and the Jordan block sits
at μ instead. At one point with c↑ = 0, three parameters held and four solved by Newton, the convergence is
quadratic, and at ninety digits

    fields (−0.2569, −1.1562, −0.2672),  rates (3.6080, 0.4112, 1.5775),  J = 0.2086,

with a resonance gap of 6·10⁻⁹¹, |c↑| = 10⁻⁹⁰ and |c↓| = 0.095 (normalised). λ is simple in both blocks, the
next eigenvalue 0.55 away, and absent from the third block (3, 2) on the line; the whole spectrum pairs about
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
eigenvectors at μ. As a control, at a resonance where both couplings are nonzero (λ ≈ −2, μ ≈ −5.2,
|c↓| = 0.056 and |c↑| = 0.031), both sides split like ε^½ with constant ratios, so the law tells the two cases
apart (gate Part C).

So at the exact point these digits approximate, L is not similar to −L − c, and no linear map carries the
palindrome. No antilinear one does either: an antilinear carrier would match L's structure at λ with its
structure at −c − λ̄ = μ̄, and L commutes with ρ ↦ ρ†, so its structure at μ̄ is its structure at μ, which is
semisimple. The eigenvalues pair, and nothing pairs them.

Read as dynamics: a mode at λ decays as t·e^{λt} in part, while its mirror partner at μ decays as a pure
exponential. The two ends of the spectrum are mirror images; the two ends of time are not.

## What we do not understand

- **Why there.** We can say when an arrow coupling vanishes, not why: there is no reading yet, in the cube's
  words or the fermions', of what makes ⟨ℓ_A|F↑|r_B⟩ vanish at a resonance while ⟨ℓ_B|F↓|r_A⟩ does not, or the
  other way round.
- **The size of the set.** The Newton steps converge quadratically with three parameters held, which suggests
  that the four conditions are independent near this point and leave a set of dimension three in the seven
  parameters. Whether such points exist on other rings, at other N, with dephasing added, is open.
- **An exact point.** The evidence is the error law at a point located to ninety digits, not exact arithmetic.
  The exact algebraic point and an exact Jordan computation there are not done.
- **What it means.** A palindrome with no carrier contradicts no statement in the repository: F158's lift to a
  similarity is for jumps with A² = 1, and its fence says so. What it means for the mirror that the spectrum
  can be symmetric while the dynamics is not, we do not know yet.

## Reproduce

```
python simulations/cooled_ring_arrow_mirror.py            # the gates, about two minutes
python simulations/cooled_ring_arrow_mirror.py --search   # the seeded search, about two and a half minutes
```

numpy, scipy, sympy and mpmath. Part A is exact in integers, Part B exact over ℚ(i), Part C the error law at
ninety digits with its control, Part D the free-fermion reading in floating point, Part E exact modulo two
primes.
