# The Double Root

*A place, not a number. What our ¼, our exceptional points and our folds have in
common, and how far the one test that tells them apart actually reaches.*

---

## Where this page comes from

The sweep, store by store. `docs/ANALYTICAL_FORMULAS.md` returned F95, F97 and
the F89 path-3 entry that reads a character off a discriminant's multiplicity.
`docs/proofs/` returned the F95 angle proof, the F97 cardioid proof,
`PROOF_ROADMAP_QUARTER_BOUNDARY` with its "One Word, Two Seams" disambiguation
and its fold normal form, `PROOF_F86_QPEAK` with the fence that three distinct
Q-thresholds live on one axis, `PROOF_F86A_EP_MECHANISM`, which names the Kato
lemma and the Puiseux exponent in one sentence, and `PROOF_CODIM1_BY_ADDITIVITY`,
which uses the order rule operationally. `experiments/` returned the rule itself
as a named lemma in `F89_SEED_EXISTENCE_REDUCTION` and the exponent table in
`F89_BETA_EXOTIC_GENERICITY`. The typed layer returned `EpCharacter` in
`compute/RCPsiSquared.Core/Numerics/`, its eig-based sibling `PhaseRigidity`
(corroborating, not load-bearing), and the roots `epcharacter`, `f89octic`,
`galoismonodromy` and `f89galois`. The open-arcs registry returned the standing
of the counting route and the gap the reviews named against the test's converse.
`docs/GLOSSARY.md` returned the unit trap this page's roster runs straight into,
and it is quoted below; neither it nor `docs/READING_GUIDE.md` carries a
defective-versus-diabolic entry.

`docs/CAUGHT_ERRORS.md` returned this page's parent:

> the EP-character trilogy is complete and internally coherent — F86a (near-EP,
> no coalescence, retracted) / coherence-horizon (genuine defective EP,
> confirmed) / F89-octic (diabolic crossing, corrected) — all three settled by
> the same artifact-free apparatus (`EpCharacter` + the discriminant
> order-of-vanishing), and each a genuinely different outcome

So the family was assembled before this page, for the spectrum, with its
verdicts, inside an errors ledger. What is new here is the family across BOTH
rooms, and a sharper account of how far the order test reaches.

---

## The object

Take a family of polynomials with one parameter, and watch two of its roots.
Almost everywhere they are two. On a thin set they are one. That thin set is the
**double root**, and it is a place in parameter space, not a value of anything.

The reason one number, ¼, keeps appearing is not that ¼ is a special number. It
is that we keep meeting the same place, and our quadratics are normalized so
that it sits at ¼.

Two rooms. **In parameter space** an equation in an unknown has two solutions
that merge as a knob is turned; nothing is an operator and there are no
eigenvectors. **In a spectrum** two eigenvalues of an operator meet, and the
eigenvectors may or may not meet with them. The two rooms are not the same room,
and `PROOF_ROADMAP_QUARTER_BOUNDARY` Layer 7 says so under *One Word, Two
Seams*: the repo uses "fold" and "cusp" for one pair in each room, and the two
spectral seams coincide only at N = 2.

---

## The test: read the ORDER of the zero

The discriminant vanishing says the roots met. **How** it vanishes says something
about the kind of meeting, without an eigensolver and without a magnitude a
finer grid could move.

The lemma is ours, and the fact under it is classical. It is stated in
[F89 Seed Existence Reduction](../experiments/F89_SEED_EXISTENCE_REDUCTION.md)
as a Lemma graded Tier 1 derived, resting on what it calls the classical Kato
fact: **a colliding pair with leading Puiseux exponent e makes the discriminant
vanish to order 2e.**

Read backwards it is a test, and the test is one-sided:

| Discriminant vanishes to | What it forces |
|---|---|
| a **simple** zero | exponent ½ exactly, so a square-root branch point and a 2×2 Jordan block. Defective. **Decided.** |
| order **≥ 2** | a diabolic crossing, an analytic-defective 2×2, a cubic branch point, or two coincident defective pairs. **Not decided.** |

**The one-sidedness is the whole discipline, and it is easy to lose in one
sentence.** A higher-order zero rules out the exponent-½ branch point. It does
**not** rule out defectiveness, because an analytic pair can be defective too:
the matrix [[0, 1], [t², 0]] has eigenvalues ±t, perfectly analytic,
discriminant 4t², a zero of order 2, and at t = 0 its geometric multiplicity is
1 against algebraic 2.

We do not need the toy. The `sideways_spin_ladder` arc measured the real thing
at the (1,3) block of N = 6: fifteen real-q defective loci, each with its
q\*, and "every observed disc zero EVEN-order". So the unscoped converse,
"defective ⟹ simple zero", is not merely underived: the registry calls it FALSE
in doubled settings. Carry the scope with it. The doubling there belongs to that
block's even character, and the registry reads the odd-character case, (1,2) at
N = 8, the other way, a generic codim-1 defective locus being a simple zero.

That is also why an even-order zero is not by itself a verdict of semisimple.
The claim class for the F89 path-3 octic states the grid-free chain "a double
zero in q ⟹ the eigenvalues cross linearly/analytically ⟹ semisimple" as its
proof, with the artifact-free discriminators listed as corroboration. The first
arrow is the lemma. The second is the step both the toy and the fifteen loci
break, so the load rests on the corroborators the claim actually lists:
geometric multiplicity 2 equal to the algebraic one, rank(L − λ·I) = n − 2, a
vanishing departure from normality, and no generalized eigenvector. That
departure is read against a relative tolerance of 1e-2 in `EpCharacter`, so it
is a measured smallness and not an exact zero.

**Two "orders" live in this section and they are different numbers.** The order
of the discriminant's zero is 1 for the ordinary defective case; the order of the
exceptional point counts how many eigenvalues meet, and is 2 in that same case.
The Kato lemma reads "a defective, order-2 exceptional point, at a simple zero"
in one sentence, and both numbers are correct.

**When two honest scans disagree about a character, the first hypothesis is two
objects, not one contradiction.** That is the ledger's lesson and it was paid
for: a min-gap scan read DIABOLIC where a witness read DEFECTIVE, and both were
right, because the pair at Re = −2γ and the pair at Re = −4γ are different pairs
in the same block. Check it is the SAME pair before calling anything a
contradiction. The ledger also names the measurable form of the order test, for
when no exact discriminant is available: the split-scaling exponent, ½ against 1.

---

## The roster

**Read the units before putting two rows on one axis.** `docs/GLOSSARY.md` keeps
the trap this table walks into: the octic's q ≈ 0.659 is 1.318 in carrier-clock
units, which lands between the horizon rungs Q\*(2) = 1 and Q\*(3) = √2, so
"convert (halve a Q, or double a q) before putting an octic q and a horizon Q on
one axis". The q in the horizon's dispersion below is a mode wavevector, another
object again; the glossary entry lists further uses of the same letter.

### In parameter space

| Site | The equation | What the repo records about the meeting |
|---|---|---|
| The CΨ = ¼ fold ([Uniqueness Proof](proofs/UNIQUENESS_PROOF.md)) | CR² + (2CΨ−1)R + CΨ² = 0, D = 1 − 4CΨ | [The Roadmap](proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md) Layer 6 puts it in normal form, "the fold catastrophe x² + a = 0 IS the recursion R = C(Ψ+R)² with a = 1−4CΨ", and [Critical Slowing](../experiments/CRITICAL_SLOWING_AT_THE_CUSP.md) carries the same content as a saddle-node, η± = ±√ε |
| [F95](proofs/PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md), the angle | z² − 2bz + c, D = 4(b² − c) | The root's doubleness is tabulated, "one degenerate real root z = b (double)". The order of D's zero in c is not stated |
| [F97](proofs/PROOF_F97_CARDIOID_HALF_FIXED_POINT.md), the cardioid cusp | z² − z + c, D = 1 − 4c, cusp at φ = 0 | The cardioid identity, verified to a max residual of 1.24e-16 over 1000 sampled φ; the order is not stated |
| The Rényi family ([Roadmap](proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md) Layer 6) | R = C_α(Ψ+R)^α | Recorded as a fold threshold in CΨ, "the value of CΨ at which its two real fixed points merge", one per α, state-independent only at α = 2 |
| **The near-miss:** [F116](proofs/PROOF_CEILING_GOLDEN_ROUTER.md), the metallic means | r² = cr + 1, D = c² + 4 | The router theorem is asserted for every **real** c, where D never vanishes, so the metallic line carries no fold on its own axis. The window lemma is stated more widely, "for all real, indeed all complex, c"; that D would vanish at c = ±2i is this page's arithmetic, not that proof's |

### In a spectrum

| Site | What meets | What is recorded about the order |
|---|---|---|
| The coherence horizon Q\*(N) ([proof](proofs/PROOF_COHERENCE_HORIZON_SLOPE.md)) | the {0,2}-coherence pair, tuned by Q = J/γ. The full-ladder dispersion is λ² + 8γλ + 4J²q², with q the mode wavevector, and q_min → π/N is where N enters that form. [Fold and Cusp](../experiments/FOLD_AND_CUSP_TWO_SEAMS.md) records that the double root is not pinned, drifting −2γ → −4γ with N, the short ladder λ² + 4γλ holding at N = 2, 3 | Named a **simple** zero in the [F89 contrast](ANALYTICAL_FORMULAS.md); defective, confirmed artifact-free |
| The F89 real seed ([reduction](../experiments/F89_SEED_EXISTENCE_REDUCTION.md)) | two eigenvalue branches of the (1,2) block's pencil, tuned by q | **Simple**, Puiseux ½, so defective by the lemma; the source calls the question closed modulo two remaining items, the non-simple-zero genericity among them |
| The F89 path-3 octic | two of eight roots at q ≈ 0.658983 | **Double**: the EP-condition (3q⁴+q²−1) enters disc(F₈) squared. Diabolic, by the discriminators and by monodromy, not by the order |
| The F89 branch points | pairs of the same eight roots, at complex q | **Simple** zeros of the squarefree layer; defective, the S₈ transpositions |
| The edge block under a γ profile ([proof](proofs/PROOF_EDGE_BLOCK_DEFECTIVE_UNDER_PROFILE.md)) | two eigenvalues of an unreduced tridiagonal, tuned by one site's rate | What is certified in exact ℚ is a **sign change**, which proves odd order; simplicity is the hypothesis, applied at the bracket ends. Defectiveness is forced separately: a non-derogatory matrix has no diabolic alternative |
| [F140](proofs/PROOF_R90_FROZEN_DIVISOR.md) §9 | the frozen root's block, tuned by J | No discriminant on either stratum: on the taxed one the couplings come off a cofactor, or off kernel dimensions of powers; on the zero-mean one the cofactor vanishes identically and they must be read from the characteristic polynomial. Defective: one Jordan block of size 2 at every exceptional coupling reached so far on the taxed stratum, and a size-3 block at one zero-mean profile, at a plus/minus pair of couplings, which the proof offers as showing the structure is not universally 2x2 rather than as a law |
| The cracked ring at u = 1 ([proof](proofs/PROOF_CRACKED_RING_EXACT_CURVE.md)) | the m ↔ N−m pair, tuned by the wrap detuning | Simplicity for every u ≥ 0 except 1 is the work of the proof, carried by a fold and a Bézout certificate; u < 0 is outside its declared range |
| The level collisions of [F129](proofs/PROOF_F129_LEVEL_COLLISION_LAW.md) / [F130](proofs/PROOF_F130_COLLISION_DECOUPLING.md) | two clean triples reaching one level, on the degenerate level-S eigenspace of the hop | Not read as an order. What is proved is that the pair does not couple, B(τ,σ) = 0 at every collision, so "no avoided crossing, ever"; the coincidence can still be shifted apart at second order, only never hybridized |
| **Not a member:** [the exceptional couplings](../experiments/THE_EXCEPTIONAL_COUPLINGS.md) | a mobile eigenvalue crosses an already-degenerate line and one more mode freezes | A pencil determinant, not a discriminant. Its roots are proved simple where the block was enumerated exactly, the whole band at N = 5 and one block at N = 6; at N = 7 and 8 simplicity is untested, and a sign change proves only odd order. Whether those points are semisimple or defective is an open item there |

**A name already spent.** "Exceptional coupling" means F140 §9's corner block on
the R₉₀ locus in one row and the uniform-point object in the other.
`THE_EXCEPTIONAL_COUPLINGS.md` has a section on exactly this, headed *The word
that is already spent*, and its rule is the one to carry: never drop which block
and which multiplicity.

---

## The one place the word is missing

F129 and F130 describe a repeated level and prove the degenerate pair does not
couple. The words "diabolic" and "semisimple" appear zero times in either proof
and zero times in `compute/MirrorWorld/LevelCollision.cs`, while they appear
across dozens of files elsewhere. The proofs say it in physics words instead.
Whether that object belongs in the spectral table above under a character label
is a question this page raises and does not settle.

---

## The question the roster raises

The order route is how a count gets closed exactly, and the live entry in the
arc `zeros_connecting_structure` says where that stands: the complex-q set "is
not claimed complete (Route B / the exact F_18 discriminant is unattempted, the
infeasibility read retired)". Unattempted, and no longer for the reason once
given.

## What this document does not claim

That the sites above are the *same* double root, or that one causes another.
They are the same KIND of place. Two quadratics sharing the shape of a
discriminant is an identity of the algebra and says nothing on its own about the
physics on either side; the repo's own fences are
`experiments/MANDELBROT_CONNECTION.md`, which declines to claim more than shared
structure, and `experiments/INFORMATION_GEOMETRY.md`, which looked for something
physical at ¼ and found the geometry finite and smooth there.

The value of the family is the test, and knowing where the test stops.
