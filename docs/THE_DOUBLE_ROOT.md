# The Double Root

*A place, not a number. What our ¼, our exceptional points and our folds have in
common, and how far the one test that tells them apart actually reaches.*

---

## Where this page comes from

The sweep, store by store. `docs/ANALYTICAL_FORMULAS.md` returned F95, F97 and
the F89 path-3 entry, where discriminant multiplicity locates and orders a
degeneracy but pair isolation, the twin-scalar restriction and `EpCharacter`
decide its character.
`docs/proofs/` returned the F95 angle proof, the F97 cardioid proof,
`PROOF_ROADMAP_QUARTER_BOUNDARY` with its "One Word, Two Seams" disambiguation
and its fold normal form, `PROOF_F86_QPEAK` with the fence that three distinct
Q-thresholds live on one axis, and `PROOF_F86A_EP_MECHANISM`, which names the
Kato lemma and the Puiseux exponent in one sentence.
`PROOF_CODIM1_BY_ADDITIVITY` owns the twin-scalar condition this page's octic
paragraph leans on. Its exact factorization separates two layers. Because A₁
is squarefree and coprime to A₂, every nonzero A₁ root is a simple discriminant
zero; w=q² is locally invertible at q≠0, so the simple-zero lemma already forces
one Puiseux-½ defective EP2 there. The doubled A₂² layer carries order-two
silent-degeneracy candidates whose pair and local character are not fixed by
the order alone. The fold-resultant certificate itself remains
character-agnostic, holding "defective or diabolic, either parity" because it
never extracts that split. `experiments/` returned
the rule itself as a named lemma in `F89_SEED_EXISTENCE_REDUCTION` and the
exponent table in `F89_BETA_EXOTIC_GENERICITY`. The typed layer returned
`EpCharacter` in `compute/RCPsiSquared.Core/Numerics/`, its eig-based sibling
`PhaseRigidity` (corroborating, not load-bearing), and the roots
`epcharacter`, `f89octic`, `galoismonodromy` and `f89galois`. The open-arcs
registry returned the standing of the counting route and the gap the reviews
named against the test's converse. `docs/GLOSSARY.md` returned the unit trap
this page's roster runs straight into, and it is quoted below; neither it nor
`docs/READING_GUIDE.md` carries a defective-versus-diabolic entry.
The F160 owner returned the continuous Hermitian road through the clean-triple
level collision; F161 returned the collision-gap series and gcd criteria along
that road through fifth order. The Python and C# confirmation registries returned two adjacent sightings.
F25 follows a 19-point Bell+ trajectory through CΨ=¼ on Kingston (RMS 0.0097):
it confirms that physical trajectory, not the discriminant order and not an
EP/Jordan character. F129's stepped-Floquet fringe supports the hardware
existence reading of that level collision, again not its order or character.
The typed sweep also returned `TransitionBridgeF95SiblingClaim`, which already
owns the narrower CΨ-cusp/F86-EP siblinghood: two quadratics, two anchors, two
spaces and no hidden identity.

`docs/CAUGHT_ERRORS.md` returned this page's parent and, later, its correction.
The older entry grouped `EpCharacter` and discriminant order as though both
settled all three character verdicts. The September 6 correction records why
that grouping was the error: order can establish an isolated transverse
square-root split, but not semisimplicity. One fence from the older entry still
belongs here:

> ‖P‖ measures oblique embedding, NOT defective-vs-diabolic — a closed-block
> defective EP reads ‖P‖≈1 (Gate0a), only the eig-Petermann K diverges, which
> `EpCharacter` deliberately avoids

So the family was assembled before this page, for the spectrum, with its
verdicts, inside an errors ledger, and the cusp/EP sibling pair was already
typed. What this page adds is the broader roster across both rooms and a
sharper account of how far the order test reaches.

---

## The object

Take a family of polynomials with one parameter, and watch two of its roots.
Almost everywhere they are two. On a thin set they coincide. The repeated
root value is the **double root**; the thin parameter set is its
**double-root locus**. This page's "place, not a number" is the locus, not a
redefinition of the root.

The reason ¼ keeps appearing in the parameter-space rows is not that every
double-root locus is the same place. It is the normalization of those
quadratics; the spectral rows meet the same KIND of place at other values.

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
for one isolated colliding pair with leading Puiseux exponent e, while every
other root gap stays nonzero, **the discriminant vanishes to order 2e.** With
several simultaneous pairs, those contributions add. The analytic-family
background is classical Kato perturbation theory; a modern matrix statement
with the generic single-Jordan-block premise is given by
[Welters](https://arxiv.org/abs/0905.4051).

The test needs more than a polynomial-looking knob. The operator or polynomial
family must be analytic in ONE local parameter; persistent factors must first
be removed so that the discriminant is not identically zero; the colliding
pair must be isolated before its exponent is inferred from the TOTAL
discriminant; and the parameter must cross the locus transversely. Order is
preserved by a locally invertible analytic coordinate change u(t) with
u'(0) ≠ 0. A tangent or ramified path multiplies it: λ²−t has discriminant 4t,
whereas the same locus read through t=s² gives λ²−s² and discriminant 4s².

F97 is the repository's own warning. In the coefficient c,
D(c)=1−4c has a simple zero at c=¼. Along the cardioid parameter
c(φ)=½e^(iφ)−¼e^(2iφ), however,
D(φ)=(e^(iφ)−1)² has order 2 at φ=0 because the curve is tangent to the
locus. That order-2 reading does not turn the underlying transverse fold into
a different singularity.

Drop analyticity of the matrix family and even a simple discriminant zero no
longer fixes the Jordan character: A(t) = [[0, √t], [√t, 0]] is continuous and
real symmetric on t ≥ 0 with char poly λ² − t, a simple zero of the
discriminant, and A(0) = 0, perfectly semisimple.

Read backwards it is a test, and the test is one-sided:

| Discriminant vanishes to | What it forces |
|---|---|
| a **simple** zero | For an analytic polynomial family and one isolated double root: exponent ½ exactly. If this is the local characteristic factor of an analytic matrix/operator and the collision has algebraic multiplicity 2, the nonanalytic split forces one 2×2 Jordan block. In parameter space there is no Jordan verdict. |
| order **≥ 2** | A diabolic crossing, an analytic-defective 2×2, a cubic branch point, several coincident defective pairs, or merely a tangent/ramified scan. **Not decided.** |

**The one-sidedness is the whole discipline, and it is easy to lose in one
sentence.** For one isolated pair in a transverse analytic coordinate, a
higher-order zero rules out exponent ½ in that coordinate. The total
discriminant does not say this until isolation and transversality are proved.
Even then, higher order does **not** rule out defectiveness, because an
analytic pair can be defective too:
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
Two different failures have to stay visible. Without pair isolation,
disc(λ³−t)=−27t² has a double zero but cubic branches t^(1/3); the fifteen
sideways-ladder loci are even for another reason, two distinct square-root
defective pairs coinciding. Even after an isolated pair has been shown to
cross analytically, the toy [[0,1],[t²,0]] shows that analyticity does not
decide its Jordan character. For the octic, pair isolation supplies the
analytic-crossing reading; it still does not supply semisimplicity.

[The Galois reading](../experiments/F89_PATH_K_GALOIS.md) already says the
first half, that the perfect square "is consistent with an analytic crossing
but is
**corroborating, not decisive, on its own**", and points at the second: the
semisimplicity is "established decisively by the scalar-λI restriction of the
octic onto the coalescing span".

The scalar restriction carries content only once scalarity has been
established. [PROOF_CODIM1_BY_ADDITIVITY](proofs/PROOF_CODIM1_BY_ADDITIVITY.md)
owns the conditional twin-scalar lemma: the residual regime is semisimple
where both compressed halves are scalar; this is proved at the N=4 octic point
and checked on the stated N=5 loci. [Diabolic by
Integrability](../hypotheses/DIABOLIC_BY_INTEGRABILITY.md) gives the Tier-2
physical explanation for N=4, H-scalar by Slater additivity and D-scalar at the
AT midpoint. Its useful sentence is correctly causal: the discriminant
double-zero is the **algebraic shadow** of twin scalarity, not an independent
cause. The Tier-1 character verdict and the Tier-2 explanation must not be
collapsed into one grade.

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
multiplicity 2", "rank(L − λ·I) = n − 2", "no generalized eigenvector" and a
scalar 2×2 compression are equivalent mathematical conditions. Departure from
normality is a differently conditioned numerical counter-reading from that
same Riesz compression, not an independent condition. `EpCharacter` reads it
against a relative tolerance of 1e-2 (with max(1, ‖A‖) as normaliser), so it is
a measured smallness and not an exact zero. The exact half is the
Hermitian-compression argument above it.

**Two "orders" live in this section and they are different numbers.** The order
of the discriminant's zero is 1 for the ordinary defective case. The order of
an exceptional point is the length of its longest coalescing eigenvector chain,
equivalently the size of the largest Jordan block in the isolated collision;
it is 2 for an EP2. It is not the total algebraic multiplicity: two simultaneous
EP2s do not become an EP4 merely because four eigenvalues meet at one parameter.
The Kato lemma reads "a defective EP2 at a simple discriminant zero" in one
sentence, and both numbers are correct.

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
degeneracy, not a Jordan structure. The quadratic rows share the local fold
normal form x²+a=0, but not one global quarter-normalized monic polynomial.
F97 and F95's b=½ specialization use x²+px+q with p=−1, so q=¼ at the double root; general F95 has p=−2b and locus c=b². In the CΨ recursion,
y=CR and s=CΨ instead give y²+(2s−1)y+s²=0; its discriminant vanishes at
s=¼ although its coefficients there are p=−½ and q=1/16. The content of this
half of the roster is where each normalization and scan direction comes from,
not merely that a discriminant exists.

### In parameter space

| Site | The equation | What the repo records about the meeting |
|---|---|---|
| The CΨ = ¼ fold ([Uniqueness Proof](proofs/UNIQUENESS_PROOF.md)) | CR² + (2CΨ−1)R + CΨ² = 0, D = 1 − 4CΨ exactly, the C²Ψ² cancelling with no hidden normalization; at C = 0 the equation is not quadratic and the row does not apply | [The Roadmap](proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md) Layer 6 puts it in normal form, "the fold catastrophe x² + a = 0 IS the recursion R = C(Ψ+R)²". Depressing the quadratic by x = R + (2CΨ−1)/(2C) and dividing by C gives that form exactly, with a = (4CΨ−1)/(4C²) = −D/(4C²) for every C ≠ 0; since 4C² > 0 the two real roots sit where a < 0, i.e. CΨ < ¼. In the local coordinate s=CΨ the zero is simple; a physical-time scan inherits that order only when ds/dt ≠ 0 at the crossing. [Critical Slowing](../experiments/CRITICAL_SLOWING_AT_THE_CUSP.md) carries the transverse content as η± = ±√ε |
| [F95](proofs/PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md), the angle | z² − 2bz + c, D = 4(b² − c) | At fixed b, scanning c transversely gives a simple zero at c=b² and the tabulated double root z=b. A different path (b(t),c(t)) needs its own order |
| [F97](proofs/PROOF_F97_CARDIOID_HALF_FIXED_POINT.md), the cardioid cusp | z² − z + c, D = 1 − 4c, cusp at φ = 0 | Simple in the transverse coefficient c. Along the cardioid itself, D(φ)=(e^(iφ)−1)² has order 2: the curve is tangent at the cusp. The 1.24e-16 numerical residual verifies the cardioid identity, not either order statement |
| The Rényi family ([Roadmap](proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md) Layer 6) | R = C_α(Ψ+R)^α | Recorded on the selected positive real branch as a fold threshold in CΨ, one per α, state-independent only at α = 2. For general α this row is locally analytic on that branch, not a polynomial family |
| **The near-miss:** [F116](proofs/PROOF_CEILING_GOLDEN_ROUTER.md), the metallic means | r² = cr + 1, D = c² + 4 | The theorem block is the c = 1 case; the extension to every **real** c is the metallic family section. There D never vanishes, so the metallic line carries no fold on its own axis. The window lemma is stated more widely, "for all real, indeed all complex, c"; that D would vanish at c = ±2i is this page's arithmetic, not that proof's |

### In a spectrum

| Site | What meets | What is recorded about the order |
|---|---|---|
| The coherence horizon Q\*(N) ([proof](proofs/PROOF_COHERENCE_HORIZON_SLOPE.md)) | the {0,2}-coherence pair of the open XY/Haken-Strobl single-excitation block under uniform local Z-dephasing, tuned by carrier Q = J/γ. **Q\*, the EP, and not the handover Q_h just below it**, which is where the FULL Liouvillian's slowest mode stops oscillating; the two are one event at N = 2, 3 and separate from N = 4 (1.87874 against 1.87854), so this row is the EP's. The large-N/long-wave full-ladder dispersion is λ² + 8γλ + 4J²q², with q the mode wavevector and q_min → π/N; it owns the asymptotic slope, not an exact finite-N quadratic. [Fold and Cusp](../experiments/FOLD_AND_CUSP_TWO_SEAMS.md) records that the double root drifts −2γ → −4γ with N. The exact short ladders λ² + 4γλ + cJ² (c = 4, 2) hold at N = 2, 3 and give Q\* = 2/√c = 1, √2 | Named a **simple** zero in the [F89 contrast](ANALYTICAL_FORMULAS.md); the live artifact-free defectiveness witness is gated at N = 2..5, while the asymptotic dispersion has its separate derivation |
| The F89 real seed ([reduction](../experiments/F89_SEED_EXISTENCE_REDUCTION.md)) | two branches of the open nearest-neighbour XY chain's uniform-Z (1,2) joint-popcount pencil, tuned by doubled-book q=J/γ | At N=5,7,9 the individually certified seeds are simple, Puiseux-½ defective EP2s. For arbitrary odd N only the nullity-surplus COUNT is a theorem; the literal finite-q drop and √-type genericity remain open outside the certified sizes (N=11 already exhibits third-order lift-off in the count bookkeeping) |
| The F89 path-3 octic | two of eight roots at q ≈ 0.658983 | **Double**: the degeneracy-locus factor (3q⁴+q²−1) enters disc(F₈) squared. Diabolic, by the discriminators and by the Hermitian-compression argument, not by the order. The EIGENVALUE monodromy (a loop around q_EP returns the identity) is not a second route: trivial eigenvalue monodromy says the branches are single-valued, which is the even-order reading again, and the toy above has trivial monodromy while being defective. Eigenvector holonomy would be a second route and is a different object |
| The F89 branch points | pairs of the same eight roots, at real AND complex q | **Simple** zeros of the squarefree layer P₂₀; defective, the S₈ transpositions. P₂₀(q) = P₁₀(q²) is EVEN, so its roots come in ± pairs and **eight** of the twenty sit on the REAL axis: the certificate's four committed loci on q > 0 (q = 0.460, 0.854, 0.857, 1.738) and their negatives. The squared factor likewise contributes two real diabolics, ±0.659. So the real axis carries both characters, and the row above is one of the real-q sites rather than the only one. The near-twin pair 0.854 / 0.857, 0.003 apart, is what a 0.05-cell lasso reads as a single point |
| The edge block under a γ profile ([proof](proofs/PROOF_EDGE_BLOCK_DEFECTIVE_UNDER_PROFILE.md)) | two eigenvalues of the open-XY (0,1) tridiagonal block, tuned by one site's rate at fixed J | Exact discriminant signs prove an odd TOTAL multiplicity of real-pair transitions inside the bracket. Only an isolated unique zero can therefore be assigned odd order. Defectiveness is forced separately: a non-derogatory matrix has no diabolic alternative |
| [F140](proofs/PROOF_R90_FROZEN_DIVISOR.md) §9 | the frozen root's block, tuned by J | No discriminant on either stratum: on the taxed one the couplings come off a cofactor (kernel dimensions of powers give the Jordan structure at a coupling already found, not the coupling itself); on the zero-mean one the cofactor vanishes identically and they must be read from the characteristic polynomial. Defective: one Jordan block of size 2 at every exceptional coupling reached so far on the taxed stratum, and a size-3 block at one zero-mean profile, at a plus/minus pair of couplings, which the proof offers as showing the structure is not universally 2x2 rather than as a law |
| The cracked ring at u = 1 ([proof](proofs/PROOF_CRACKED_RING_EXACT_CURVE.md)) | the m ↔ N−m pair of the real-symmetric single-excitation XY adjacency Hamiltonian, tuned by wrap detuning u; γ is absent | Simplicity for every u ≥ 0 except 1 is the work of the proof, carried by a fold and a Bézout certificate; u < 0 is outside its declared range. At u=1 the Hermitian ring degeneracy is semisimple: independent eigenvectors, not a Liouvillian EP and not noise-induced |
| **Not a member:** [the exceptional couplings](../experiments/THE_EXCEPTIONAL_COUPLINGS.md) | a mobile eigenvalue crosses an already-degenerate line and one more mode freezes | A pencil determinant, not a discriminant. Its roots are proved simple where the block was enumerated exactly, the whole band at N = 5 and one block at N = 6; at N = 7 and 8 simplicity is untested, and a sign change proves only odd order. Whether those points are semisimple or defective is an open item there |

**A name already spent.** "Exceptional coupling" means F140 §9's corner block on
the R₉₀ locus in one row and the uniform-point object in the other.
`THE_EXCEPTIONAL_COUPLINGS.md` has a section on exactly this, headed *The word
that is already spent*, and its rule is the one to carry: never drop which
block and which multiplicity.

---

## A related fixed degeneracy outside the test

F129 is an arithmetic equality between clean triples at discrete comb modulus
n, not an analytic one-knob family, so it does not belong in the order-test
table. F130 proves second-order non-mixing; diagonal self-shifts may still
separate the equal levels. F160 supplies the missing continuous road: its
one-particle H₁(u)=H_chain+uV is a real-analytic Hermitian matrix family, so
Rellich's theorem supplies local real-analytic eigenvalue parametrizations;
their three-mode sums inherit that regularity. The F129 equality is a repeated level of the induced three-fermion
exterior-power Hamiltonian, which is Hermitian for real u; at u=0 it is
therefore semisimple. F161 reads the collision-gap series and gcd criteria
along this road through fifth order; it does not supply the road or the
character. The mathematical character is settled; whether the repository
spends the house word "diabolic" on this Hermitian collision is only a naming
choice. The hardware confirmation sees the stepped-Floquet near-degeneracy and
uses verified θ³ scaling to infer the continuum collision; it does not measure
the exact equality, discriminant order or character directly.

---

## The question the roster raises

The order route is how a count gets closed exactly, and the live owner is the
arc `diabolic_over_higher_n`. Route B concerns the exact discriminants of the
path-4 N=5, (1,2) residuals: F_18 in R-even (q-degree 274) and the distinct
F_17 in R-odd (q-degree 246), with no M-descent shortcut because this block is
not fold-fixed. The exact full-degree factorization has already landed: in
w=q² each sector obeys `disc_Λ(F_res)=C·w^v·A₁·A₂²`, with squarefree coprime
layers; at N=5 both A₁ and A₂ are irreducible over ℚ. The q-degree layers are
[56,32] in R-even and [56,26] in R-odd, with v_q=154/138. Route B's complete
root/locus inventory, including certified complex-root isolation, is now given
below. The remaining work is pair isolation and a character gate on the doubled
A₂² loci. The A₁ loci need no second character gate: squarefreeness and
coprimality make their discriminant zeros simple, and the simple-zero lemma
already certifies Puiseux-½ defective EP2 character.

That first remaining Route B step is closed exactly at N=5. Exact rational
root isolation gives the complete inventory below; each tuple is
`(negative-real, positive-real, nonreal)` in w. Because every listed root is
nonzero, `w=q²` lifts each one to two q-loci.

| R parity | layer | exact w-root inventory | w roots | induced q-loci | character |
|---|---:|---:|---:|---:|---|
| even | A₁ | (0, 2, 26) | 28 | 56 | all Puiseux-½ defective EP2 |
| odd | A₁ | (0, 2, 26) | 28 | 56 | all Puiseux-½ defective EP2 |
| even | A₂ | (6, 0, 10) | 16 | 32 | locally undecided |
| odd | A₂ | (6, 1, 6) | 13 | 26 | locally undecided except the gated positive-real diabolic |

Thus A₁ is complete at **56 w roots / 112 q-loci**, all certified EP2. A₂ is
complete as a locus inventory at **29 w roots / 58 q-loci**, but its local
character classification is the next Route B step. The certificate is the
default N=5 run of `simulations/o2b_gcd_certificate.py`; it gates the four count
tuples in addition to the earlier exact factorization and real-positive anchors.

## What this document does not claim

That the sites above are the *same* double root, or that one causes another.
They are the same KIND of place. Two quadratics sharing the shape of a
discriminant is an identity of the algebra and says nothing on its own about
the physics on either side; the repo's own fences are
`experiments/MANDELBROT_CONNECTION.md`, which declines to claim more than
shared structure, and `experiments/INFORMATION_GEOMETRY.md`, which looked for
something physical at ¼ and found the geometry finite and smooth there.

The value of the family is the test, and knowing where the test stops.
