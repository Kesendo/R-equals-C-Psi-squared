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
Q-thresholds live on one axis, and `PROOF_F86A_EP_MECHANISM`, which names the
Kato lemma and the Puiseux exponent in one sentence.
`PROOF_CODIM1_BY_ADDITIVITY` owns the twin-scalar condition this page's octic
paragraph leans on, and separately assigns the layers, "the simple layer
carries the √-branch loci, the doubled layer the silent diabolics". What it
does NOT do is read a character off an ORDER: its fold-resultant certificate
is deliberately character-agnostic, holding "defective or diabolic, either
parity" precisely because it never extracts the split. `experiments/` returned
the rule itself as a named lemma in `F89_SEED_EXISTENCE_REDUCTION` and the
exponent table in `F89_BETA_EXOTIC_GENERICITY`. The typed layer returned
`EpCharacter` in `compute/RCPsiSquared.Core/Numerics/`, its eig-based sibling
`PhaseRigidity` (corroborating, not load-bearing), and the roots
`epcharacter`, `f89octic`, `galoismonodromy` and `f89galois`. The open-arcs
registry returned the standing of the counting route and the gap the reviews
named against the test's converse. `docs/GLOSSARY.md` returned the unit trap
this page's roster runs straight into, and it is quoted below; neither it nor
`docs/READING_GUIDE.md` carries a defective-versus-diabolic entry.

`docs/CAUGHT_ERRORS.md` returned this page's parent:

> the EP-character trilogy is complete and internally coherent — F86a (near-EP,
> no coalescence, retracted) / coherence-horizon (genuine defective EP,
> confirmed) / F89-octic (diabolic crossing, corrected) — all three settled by
> the same artifact-free apparatus (`EpCharacter` + the discriminant
> order-of-vanishing), and each a genuinely different outcome

and, further down the same bullet, the fence on that apparatus, which belongs
on a page about how far a test reaches:

> ‖P‖ measures oblique embedding, NOT defective-vs-diabolic — a closed-block
> defective EP reads ‖P‖≈1 (Gate0a), only the eig-Petermann K diverges, which
> `EpCharacter` deliberately avoids

So the family was assembled before this page, for the spectrum, with its
verdicts, inside an errors ledger. What is new here is the family across BOTH
rooms, and a sharper account of how far the order test reaches.

---

## The object

Take a family of polynomials with one parameter, and watch two of its roots.
Almost everywhere they are two. On a thin set they are one. That thin set is
the
**double root**, and it is a place in parameter space, not a value of anything.

The reason one number, ¼, keeps appearing is not that ¼ is a special number.
It is that we keep meeting the same place, and our quadratics are normalized
so that it sits at ¼.

Two rooms. **In parameter space** an equation in an unknown has two solutions
that merge as a knob is turned; nothing is an operator and there are no
eigenvectors. **In a spectrum** two eigenvalues of an operator meet, and the
eigenvectors may or may not meet with them. The two rooms are not the same
room, and `PROOF_ROADMAP_QUARTER_BOUNDARY` Layer 7 says so under *One Word,
Two Seams*: the repo uses "fold" and "cusp" for one pair in each room, and the
two spectral seams coincide only at N = 2.

---

## The test: read the ORDER of the zero

The discriminant vanishing says the roots met. **How** it vanishes says
something about the kind of meeting, without an eigensolver and without a
magnitude a finer grid could move.

The lemma is ours, and the fact under it is classical. It is stated in [F89
Seed Existence Reduction](../experiments/F89_SEED_EXISTENCE_REDUCTION.md) as a
Lemma graded Tier 1 derived, resting on what it calls the classical Kato fact:
**a colliding pair with leading Puiseux exponent e makes the discriminant
vanish to order 2e.**

The lemma needs three things, and they hold everywhere on this page because
every knob here enters polynomially: the family depends analytically on ONE
parameter, "order" means order of vanishing at a point of that parameter's
axis, and "the discriminant" is disc_λ of the characteristic polynomial read
as a function of that parameter. Drop analyticity and the forward direction
goes: A(t) = [[0, √t], [√t, 0]] is continuous and real symmetric on t ≥ 0 with
char poly λ² − t, a SIMPLE zero of the discriminant, and A(0) = 0, perfectly
semisimple.

Read backwards it is a test, and the test is one-sided:

| Discriminant vanishes to | What it forces |
|---|---|
| a **simple** zero | exponent ½ exactly, so a square-root branch point and a 2×2 Jordan block. Defective. **Decided.** |
| order **≥ 2** | a diabolic crossing, an analytic-defective 2×2, a cubic branch point, or two coincident defective pairs. **Not decided.** |

**The one-sidedness is the whole discipline, and it is easy to lose in one
sentence.** A higher-order zero rules out the exponent-½ branch point. It does
**not** rule out defectiveness, because an analytic pair can be defective too:
the matrix [[0, 1], [t², 0]] has eigenvalues ±t, perfectly analytic,
discriminant 4t², a zero of order 2, and at t = 0 its geometric multiplicity
is 1 against algebraic 2.

We do not need the toy. The `sideways_spin_ladder` arc measured the real thing
at the (1,3) block of N = 6: fifteen real-q defective loci, each with its q\*,
and "every observed disc zero EVEN-order". That evenness is no longer only
observed: [Path-K Diabolic](../experiments/F89_PATH_K_DIABOLIC.md) derives it
from the composition identity disc_Λ(F) = ±4^m·f·disc_M(G)², which with f a
perfect square admits no multiplicity-1 or multiplicity-3 layer at all. So the
unscoped converse, "defective ⟹ simple zero", is not merely underived: the
registry calls it FALSE in doubled settings. Carry the scope with it. The
doubling there belongs to that block's even character. For the odd-character
case, (1,2) at N = 8, a generic codim-1 defective locus being a simple zero,
the registry records the reading and fences it in the same breath: "a
gitignored design note only, no committed store", and it went moot for that
block the same day, the certificate finding no real roots q ≠ 0 of the
RESIDUAL to read at all (the AT factor keeps its own exact degeneracies,
semisimple by construction and untouched by that statement); the stronger
committed form is that the (1,2) block at N = 8 carries no real-q defective EP
whatever. So the even-character doubling is what the axis actually shows, and
the other side of it is not yet a citeable reading.

That is also why an even-order zero is not by itself a verdict of semisimple.
The claim class for the F89 path-3 octic states the grid-free chain "a double
zero in q ⟹ the eigenvalues cross linearly/analytically ⟹ semisimple" as its
proof, with the artifact-free discriminators listed as corroboration.

**Neither arrow carries the chain by itself, and they fail for different reasons.**
The first arrow runs the lemma BACKWARDS, and backwards it needs a hypothesis
the chain does not state: that exactly TWO roots collide. Without it, disc(λ³
− t) = −27t² is a double zero whose branches are t^(1/3), no analytic crossing
at all. The fifteen loci break this arrow too, and by the mechanism the
registry gives: their zeros are even because two DISTINCT defective pairs
coincide, each still exponent ½, so the eigenvalues there are not analytic
either. The second arrow is what the toy breaks: ±t is as analytic as a
crossing gets and the matrix is defective. On the octic the first arrow does
hold, because "two of eight roots" is known separately; that scope is the
whole of its content and it is not in the chain.

[The Galois reading](../experiments/F89_PATH_K_GALOIS.md) already says the
first half, that the perfect square "is consistent with an analytic crossing
but is
**corroborating, not decisive, on its own**", and points at the second: the
semisimplicity is "established decisively by the scalar-λI restriction of the
octic onto the coalescing span".

**Taken as it stands that phrase is nearly circular**, since "the restriction is
λ·I, therefore semisimple" puts all the content into establishing scalarity
and supplies no route to it. The route is committed and it is a theorem, not a
docstring: [PROOF_CODIM1_BY_ADDITIVITY](proofs/PROOF_CODIM1_BY_ADDITIVITY.md)
Regime 2 owns the twin-scalar condition, and [Diabolic by
Integrability](../hypotheses/DIABOLIC_BY_INTEGRABILITY.md) holds the N = 4
mechanism, H-scalar by Slater additivity and D-scalar at the AT midpoint. That
source states this page's whole subject in one line: "the discriminant
double-zero is the **algebraic shadow** of this twin scalarity, not an
independent cause". The order sits downstream of the mechanism, which is why
reading a character off it can only ever be one-sided.

`EpCharacter` implements it and its docstring states it compactly: split the
pencil on the coalescing GENERALIZED eigenspace (the orthonormal range of the
Riesz projector, not the eigenspace, or the argument would be vacuous). Given
the H-half restriction scalar, which free-fermion additivity supplies, and a
coalescence, the D-half follows, because A₂ = M₂ − qC₂ is then Hermitian with
a repeated eigenvalue and a Hermitian 2×2 with a repeated eigenvalue IS
scalar. So at real q a diabolic is automatically twin-scalar. That is the
load-bearing step, and it is conditional on the H-half, not unconditional.

The artifact-free discriminators sit beside it as confirmation, and there are
fewer of them than the list suggests: at algebraic multiplicity 2, "geometric
multiplicity 2", "rank(L − λ·I) = n − 2" and "no generalized eigenvector" are
one condition written three ways, sharing one failure mode, a numerical rank
taken at a float q near an algebraic q\*. The independent one is the departure
from normality, and it is read against a relative tolerance of 1e-2 in
`EpCharacter` (relative with a floor, the normaliser being max(1, ‖A‖)), so it
is a measured smallness and not an exact zero. Nothing in this paragraph is
exact arithmetic; the exact half is the Hermitian-compression argument above
it.

**Two "orders" live in this section and they are different numbers.** The order
of the discriminant's zero is 1 for the ordinary defective case; the order of
the exceptional point counts how many eigenvalues meet, and is 2 in that same
case. The Kato lemma reads "a defective, order-2 exceptional point, at a
simple zero" in one sentence, and both numbers are correct.

**When two honest scans disagree about a character, the first hypothesis is two
objects, not one contradiction.** That is the ledger's lesson and it was paid
for: a min-gap scan read DIABOLIC where a witness read DEFECTIVE, and both
were right, because the pair at Re = −2γ and the pair at Re = −4γ are
different pairs in the same block. Check it is the SAME pair before calling
anything a contradiction. The ledger also names the measurable form of the
order test, for when no exact discriminant is available: the split-scaling
exponent, ½ against 1.

---

## The roster

**"Simple" is spent twice in these tables, and this page of all pages owes the
distinction.** A **simple zero** is an order-of-vanishing of a discriminant; a
**simple root** or **simple spectrum** means distinct, no collision at all. The
F89 rows use the first, the cracked-ring and exceptional-coupling rows the
second. The two are nearly opposite: a simple ZERO is where a collision
happens, a simple SPECTRUM is where none does.

**Read the units before putting two rows on one axis.** `docs/GLOSSARY.md` keeps
the trap this table walks into: the octic's q ≈ 0.659 is 1.318 in
carrier-clock units, which lands between the horizon rungs Q\*(2) = 1 and
Q\*(3) = √2, so "convert (halve a Q, or double a q) before putting an octic q
and a horizon Q on one axis". The q in the horizon's dispersion below is a
mode wavevector, another object again; the glossary entry lists further uses
of the same letter.

**The test crosses only one of these tables.** In parameter space there are no
eigenvectors, so "defective" and "diabolic" are not available verdicts there
at all; those rows are folds, and what the order reads is the fold's
degeneracy, not a Jordan structure. Read down the parameter-space table and
the three quadratics are one monic quadratic in three costumes, related by
affine reparametrisation: for x² + px + q the double root sits at q = p²/4,
and ¼ appears exactly when p = ±1. That is an identity of the algebra. The
content of that half of the roster is where each normalization comes from, not
that a discriminant exists.

### In parameter space

| Site | The equation | What the repo records about the meeting |
|---|---|---|
| The CΨ = ¼ fold ([Uniqueness Proof](proofs/UNIQUENESS_PROOF.md)) | CR² + (2CΨ−1)R + CΨ² = 0, D = 1 − 4CΨ exactly, the C²Ψ² cancelling with no hidden normalization; at C = 0 the equation is not quadratic and the row does not apply | [The Roadmap](proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md) Layer 6 puts it in normal form, "the fold catastrophe x² + a = 0 IS the recursion R = C(Ψ+R)²". Depressing the quadratic by x = R + (2CΨ−1)/(2C) and dividing by C gives that form exactly, with a = (4CΨ−1)/(4C²) = −D/(4C²) for every C ≠ 0; since 4C² > 0 the two real roots sit where a < 0, i.e. CΨ < ¼. The knob is D only up to the positive factor 4C², which on the theorem's domain 0 ≤ C ≤ 1 is 1 only at C = ½ (the roadmap also states the algebra for any real C, where C = −½ does it too), and [Critical Slowing](../experiments/CRITICAL_SLOWING_AT_THE_CUSP.md) carries the same content as a saddle-node, η± = ±√ε |
| [F95](proofs/PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md), the angle | z² − 2bz + c, D = 4(b² − c) | The root's doubleness is tabulated, "one degenerate real root z = b (double)". The order of D's zero in c is not stated |
| [F97](proofs/PROOF_F97_CARDIOID_HALF_FIXED_POINT.md), the cardioid cusp | z² − z + c, D = 1 − 4c, cusp at φ = 0 | The cardioid identity, verified to a max residual of 1.24e-16 over 1000 sampled φ; the order is not stated |
| The Rényi family ([Roadmap](proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md) Layer 6) | R = C_α(Ψ+R)^α | Recorded as a fold threshold in CΨ, "the value of CΨ at which its two real fixed points merge", one per α, state-independent only at α = 2 |
| **The near-miss:** [F116](proofs/PROOF_CEILING_GOLDEN_ROUTER.md), the metallic means | r² = cr + 1, D = c² + 4 | The theorem block is the c = 1 case; the extension to every **real** c is the metallic family section. There D never vanishes, so the metallic line carries no fold on its own axis. The window lemma is stated more widely, "for all real, indeed all complex, c"; that D would vanish at c = ±2i is this page's arithmetic, not that proof's |

### In a spectrum

| Site | What meets | What is recorded about the order |
|---|---|---|
| The coherence horizon Q\*(N) ([proof](proofs/PROOF_COHERENCE_HORIZON_SLOPE.md)) | the {0,2}-coherence pair, tuned by Q = J/γ. **Q\*, the EP, and not the handover Q_h just below it**, which is where the FULL Liouvillian's slowest mode stops oscillating; the two are one event at N = 2, 3 and separate from N = 4 (1.87874 against 1.87854), so this row is the EP's. The full-ladder dispersion is λ² + 8γλ + 4J²q², with q the mode wavevector, and q_min → π/N is where N enters that form. [Fold and Cusp](../experiments/FOLD_AND_CUSP_TWO_SEAMS.md) records that the double root is not pinned, drifting −2γ → −4γ with N, the short ladder λ² + 4γλ + cJ²q² (c = 4, 2) holding at N = 2, 3, whose double root gives Q* = 2/√c = 1, √2; the truncation λ² + 4γλ is shorthand for the λ-linear coefficient and has discriminant 16γ², never zero | Named a **simple** zero in the [F89 contrast](ANALYTICAL_FORMULAS.md); defective, confirmed artifact-free |
| The F89 real seed ([reduction](../experiments/F89_SEED_EXISTENCE_REDUCTION.md)) | two eigenvalue branches of the (1,2) block's pencil, tuned by q | **Simple**, Puiseux ½, so defective by the lemma; the source calls the question closed modulo two remaining items, the non-simple-zero genericity among them |
| The F89 path-3 octic | two of eight roots at q ≈ 0.658983 | **Double**: the EP-condition (3q⁴+q²−1) enters disc(F₈) squared. Diabolic, by the discriminators and by the Hermitian-compression argument, not by the order. The EIGENVALUE monodromy (a loop around q_EP returns the identity) is not a second route: trivial eigenvalue monodromy says the branches are single-valued, which is the even-order reading again, and the toy above has trivial monodromy while being defective. Eigenvector holonomy would be a second route and is a different object |
| The F89 branch points | pairs of the same eight roots, at real AND complex q | **Simple** zeros of the squarefree layer P₂₀; defective, the S₈ transpositions. P₂₀(q) = P₁₀(q²) is EVEN, so its roots come in ± pairs and **eight** of the twenty sit on the REAL axis: the certificate's four committed loci on q > 0 (q = 0.460, 0.854, 0.857, 1.738) and their negatives. The squared factor likewise contributes two real diabolics, ±0.659. So the real axis carries both characters, and the row above is one of the real-q sites rather than the only one. The near-twin pair 0.854 / 0.857, 0.003 apart, is what a 0.05-cell lasso reads as a single point |
| The edge block under a γ profile ([proof](proofs/PROOF_EDGE_BLOCK_DEFECTIVE_UNDER_PROFILE.md)) | two eigenvalues of an unreduced tridiagonal, tuned by one site's rate | What is certified in exact ℚ is a **sign change**, which proves odd order; simplicity is the hypothesis, applied at the bracket ends. Defectiveness is forced separately: a non-derogatory matrix has no diabolic alternative |
| [F140](proofs/PROOF_R90_FROZEN_DIVISOR.md) §9 | the frozen root's block, tuned by J | No discriminant on either stratum: on the taxed one the couplings come off a cofactor (kernel dimensions of powers give the Jordan structure at a coupling already found, not the coupling itself); on the zero-mean one the cofactor vanishes identically and they must be read from the characteristic polynomial. Defective: one Jordan block of size 2 at every exceptional coupling reached so far on the taxed stratum, and a size-3 block at one zero-mean profile, at a plus/minus pair of couplings, which the proof offers as showing the structure is not universally 2x2 rather than as a law |
| The cracked ring at u = 1 ([proof](proofs/PROOF_CRACKED_RING_EXACT_CURVE.md)) | the m ↔ N−m pair, tuned by the wrap detuning | Simplicity for every u ≥ 0 except 1 is the work of the proof, carried by a fold and a Bézout certificate; u < 0 is outside its declared range |
| The level collisions of [F129](proofs/PROOF_F129_LEVEL_COLLISION_LAW.md) / [F130](proofs/PROOF_F130_COLLISION_DECOUPLING.md) | two clean triples reaching one level, on the degenerate level-S eigenspace of the hop | Not read as an order. What is proved is that the pair does not couple, B(τ,σ) = 0 at every collision, so "no avoided crossing, ever"; the coincidence can still be shifted apart at second order, only never hybridized. F130 alone carries the shift, and it is DIAGONAL: the colliding triples differ in Σλ² over their own modes, λ_k = 2cos(kπ/n), 9.73 against 2.27, which is why a coincidence can move without the pair ever repelling |
| **Not a member:** [the exceptional couplings](../experiments/THE_EXCEPTIONAL_COUPLINGS.md) | a mobile eigenvalue crosses an already-degenerate line and one more mode freezes | A pencil determinant, not a discriminant. Its roots are proved simple where the block was enumerated exactly, the whole band at N = 5 and one block at N = 6; at N = 7 and 8 simplicity is untested, and a sign change proves only odd order. Whether those points are semisimple or defective is an open item there |

**A name already spent.** "Exceptional coupling" means F140 §9's corner block on
the R₉₀ locus in one row and the uniform-point object in the other.
`THE_EXCEPTIONAL_COUPLINGS.md` has a section on exactly this, headed *The word
that is already spent*, and its rule is the one to carry: never drop which
block and which multiplicity.

---

## The one place the word is missing

F129 and F130 describe a repeated level and prove the degenerate pair does not
couple. The words "diabolic" and "semisimple" appear zero times in either
proof and zero times in `compute/MirrorWorld/LevelCollision.cs`, while they
appear across dozens of files elsewhere. The proofs say it in physics words
instead. Whether that object belongs in the spectral table above under a
character label is a question this page raises and does not settle.

---

## The question the roster raises

The order route is how a count gets closed exactly, and the live entry in the
arc `zeros_connecting_structure` says where that stands: the complex-q set "is
not claimed complete (Route B / the exact F_18 discriminant is unattempted,
the infeasibility read retired; see the obstacle note in this registry)".
Unattempted, and no longer for the reason once given.

## What this document does not claim

That the sites above are the *same* double root, or that one causes another.
They are the same KIND of place. Two quadratics sharing the shape of a
discriminant is an identity of the algebra and says nothing on its own about
the physics on either side; the repo's own fences are
`experiments/MANDELBROT_CONNECTION.md`, which declines to claim more than
shared structure, and `experiments/INFORMATION_GEOMETRY.md`, which looked for
something physical at ¼ and found the geometry finite and smooth there.

The value of the family is the test, and knowing where the test stops.
