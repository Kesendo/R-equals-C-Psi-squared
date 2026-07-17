# The β-exotic genericity: is each count-dropping seed a √-type EP2?

*Extracted 2026-07-13 from [F89_SEED_EXISTENCE_REDUCTION.md](F89_SEED_EXISTENCE_REDUCTION.md) as its own
concern (the F89 corpus de-monolith, step 1). That note proves the seed COUNT, r(0⁺) − r(∞) = N − 1;
this one asks whether each of those count-dropping seeds is a generic **√-type** defective exceptional
point (an EP2, geometric multiplicity 1, Puiseux exponent ½, s₆ ≠ 0) rather than a non-generic order-≥3
structure. It is settled at N = 5, 7 AND 9, all unconditional over ℚ (N = 5 and 7 since
2026-07-16: the two named N = 7 premises are both discharged, the base-polynomial CRT grade by the
a-priori ℓ1 bound and the mod-p layer identification by reconstruct-and-verify, see the two
premise-discharge subsections; N = 9 since 2026-07-17, the same engine at scale with both
mechanisms built in), both R-parities, and is **open for all N**.
Shared objects live in the parent note and the C# core: the (1,2)-block pencil L(q) = A + qC and its
recorded real defective seeds (`WeightCoherenceBlock`, `RealDefectiveSeeds`); the typed claim is
`BetaExoticPerNExclusionClaim` (Tier 1 derived, scoped); the verifiers are the `o2b_*` scripts (see
Reproduce, below).*

**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Registry:** F89h in [ANALYTICAL_FORMULAS.md](../docs/ANALYTICAL_FORMULAS.md) (letter minted 2026-07-13, the de-monolith)

## The β-exotic, sharpened: a count-dropping seed is √-type iff s₆ ≠ 0 (open, 2026-07-08; taxonomy refined 2026-07-09)

Pieces 1-3 of the parent [F89_SEED_EXISTENCE_REDUCTION.md](F89_SEED_EXISTENCE_REDUCTION.md) prove the
count drops by exactly N − 1. That a drop happens does not yet say each drop is a
**√-type** defective exceptional point (an EP2, geometric multiplicity 1, Puiseux exponent ½) rather
than one of the two non-generic order-3 structures: the
**β-exotic**, a real-to-complex transition through a formally semisimple, higher-order-degenerate point
(normal form β(s) = [[0, s], [s², 0]], eigenvalues ±s^{3/2}), or its *defective twin*, the
exponent-3/2 EP2 (geometric multiplicity 1 with vanishing leading branch coefficient; the "Three
attacks" section below is where it enters the story). Excluding both for all odd N is the last
item. This section records how far it reduces, and why the natural proof route does not finish it.

**The reduction (sound, and sharper than the discriminant's Galois group).** For any eigenvector r of
the pencil L(q) = D + qC (D the −2/−6 dephasing diagonal, C = iK the coherent hop, K real symmetric, so
L is complex-symmetric, Lᵀ = L), from Cr = (1/q)(λ − D)r:

> rᵀCr = (1/q)[(λ + 2)·rᵀr + 4·s₆],   s₆ := Σ_{i on the −6 rung} (r_i)²   (transpose sense).

At a coalescence the coalescing eigenvector is self-orthogonal (rᵀr → 0, the complex-symmetric EP
signature), so **rᵀCr → (4/q\*)·s₆**. And rᵀCr is exactly the branch-coefficient numerator: for a
defective coalescence of Jordan size m the branch coefficient μ (the splitting is ∼ μ·(q − q\*)^{1/m})
obeys μ^m = rᵀCr / (rᵀr_{m−1}). An
*isolated* count-dropping seed (a drop of exactly 2) is a real-to-complex transition, hence an
*odd*-order zero of the relevant real discriminant factor (two seeds coinciding at one q\* would drop
the count by 4 at an even-order zero; that case is handled separately in the certificate section).
Therefore, at an isolated forced seed, and under H1 below (algebraic multiplicity exactly 2),

> **s₆ ≠ 0 ⟺ the zero is simple (order 1) ⟺ a √-type defective EP2** (exponent ½; s₆ read at
> geometric multiplicity 1). **The zero has order ≥ 3 ⟺ the seed is not √-type**: at geometric
> multiplicity 1 that is the defective exponent-3/2 twin, with s₆ = 0, and at geometric multiplicity 2
> the semisimple β-exotic, where s₆ has no invariant value; both die per N by the same
> disc-multiplicity bound.

The s₆-reading is scoped to geometric multiplicity 1, where the coalescing eigenvector, hence s₆,
is well-defined up to the gauge. At a geometric-multiplicity-2 point (the semisimple β-exotic) s₆
depends on the eigenvector choice and carries no invariant content; there the classification runs
through the kernel dimension itself, which is why the direct test further down is defectiveness, not
s₆.

So the order-3 exclusion (β-exotic and twin alike) **reduces to: s₆ ≠ 0 at every forced seed, for all
odd N.** This is a genuine
sharpening of the dead discriminant-Galois route (which died because a full symmetric Galois group does
not force a squarefree discriminant: the group sees which sheets swap, never the ½-vs-3/2 *exponent*).
The quantity s₆ sees exactly that exponent.

**Caveat (H1), stated honestly.** The step "s₆ = 0 ⟺ an order-≥3 transition" assumes the seed has algebraic
multiplicity exactly 2 (no Jordan block of size ≥ 3; for the defective sub-case that is the m = 2
well-definedness rᵀr₁ ≠ 0 of the branch-coefficient formula above, with r₁ the Jordan generalized vector). The count-dropping structure H1 must rule out is a **Jordan
block of size ≥ 4** (disc order ≥ 3, so caught by the same multiplicity bound as the β-exotic); the
smallest, a 4×4 (exponent ¼), drops the count 2 real → 0 and is nonetheless not a defective EP2. (Even
Jordan sizes drop the count, odd sizes ≥ 5 do not; H1 excludes all sizes ≥ 3 regardless.) A 3×3 Jordan
is *not* such a confounder:
it is a cubic branch point (exponent ⅓, disc order 2, one real branch and one conjugate pair on both
sides of q\*), so it does not drop the real count and cannot be a forced seed at all (see the
multiplicity table in the certificate section below). Algebraic-multiplicity-2 holds numerically at
every seed through N = 11 but is not itself proved for all odd N; it is a premise separate from
s₆ ≠ 0, and the sign law bears only on the latter. (At N = 5 and N = 7 H1 is no longer numerical: it
is a corollary of the certificate below, "The β-exotic is excluded at N = 5 and N = 7"; at N = 9 it
follows the same way from the gcd-certificate section's proved layer identity, max multiplicity 2
plus the simple-zero lemma on A₁, with the same per-N reality check. For all odd N
it remains open, alongside the reduction of this section.)

**Numerics (grounding for the identity, not evidence for all N).** At the 17 census seeds
N = 5, 7, 9 (the census loci; the sign-law sweep below counts 19 once the re-entrant coalescences
beyond the census scan window are included, one each at N = 7 and N = 9) the anchored self-orthogonal eigenvector (a Takagi-minimal vᵀv extraction at the refined q\*, i.e. the
direction in the colliding 2-plane that minimizes |vᵀv|, the isotropic one; a naive
min-gap finder drifts to the −6-rung edge crossings at λ = −6, which are not self-orthogonal and not the
seeds) gives |s₆| ∈ [0.025, 0.38] (magnitudes; the sign of s₆ is a gauge convention, see the trap at the
end of this section), bounded away from zero, with rᵀCr = (4/q\*)s₆ holding to 0.5%. This
confirms the identity and the reduction machinery only; it is the same evidence class as the N ≤ 11
census and carries no all-N weight on its own.

**Why the natural proof route does not close it.** Splitting the eigen-equation by rung yields transpose
relations (α = γ = −β, with α = r₂ᵀK₂₂r₂, β = r₂ᵀK₂₆r₆, γ = r₆ᵀK₆₆r₆) and Hermitian relations (the
positive norms with (λ + 2)‖r₂‖² + (λ + 6)‖r₆‖² = 0, forcing λ ∈ (−6, −2)). These moment relations do
not force s₆ ≠ 0: an explicit in-class matrix (complex-symmetric and pseudo-Hermitian, two rungs) is a
clean defective 2×2 (geometric multiplicity 1, algebraic 2, rᵀr₁ = 1) with s₂ = s₆ = 0 satisfying all of
them (K = [[0,2,0,0],[2,0,−4,0],[0,−4,0,−2],[0,0,−2,0]], D = diag(−2, −2, −6, −6),
T = diag(1, −1, 1, −1), r = (1, i, 1, −i), at q = 1, λ = −4). It is an
*analytic-defective* 2×2 in the certificate section's sense, not the exponent-½ EP2 of the reduction:
s₆ = 0 forces rᵀCr = 0, so μ² = 0 and its two branches stay analytic (no √-split); such a point keeps the
real-eigenvalue count unchanged across q\* (this witness holds a complex-conjugate pair that only touches
the axis at q\* = 1, so the real count is 0 on both sides), does not drop the count, and is thereby barred
from a forced seed, which is exactly why its s₆ = 0 is consistent with the count-dropping biconditional
above. Three
structural reasons: the transpose quantity s₆ is decoupled from every positivity-bearing relation; the
pseudo-Hermitian relation L†T = TL collapses to (λ\* − λ)r†Tr = 0, vacuous on the real axis; so any real
proof must re-enter the exact spectral arithmetic of the two closed-form blocks, K₂₂ (the N − 1 paths of
Piece 1) and K₆₆ = −3H₃ (Piece 3). The precise missing ingredient is a *definiteness*: a coercivity from
the −3H₃ spectrum, or a sign-rigidity from the path structure, that pins the indefinite, opposite-Krein-
sign colliding pair off the s₆ = 0 locus. It is unidentified.

**A named candidate for that ingredient: the class-imbalance sign law (2026-07-09).** The
"sign-rigidity from the path structure" just called unidentified now has a concrete, gauge-invariant
candidate. Split each rung by the bipartite sign into class **E** (T = +1) and class **O** (T = −1);
the rung's **class imbalance** |E_rung| − |O_rung| is (computed exactly at each N) **+(N − 1)** on the
−2 rung, equal to nullity(K₂₂), and **−3(N − 1)/2** on the −6 rung. The conjecture is that the
gauge-invariant Krein index κ_rung := v†·T·P_rung·v / ‖v‖² (P_rung the rung projector; the trap
paragraph at the end of this section holds the gauge discussion) tracks that imbalance,

> **sign(κ_rung) = sign(|E_rung| − |O_rung|) at every forced seed,**  i.e. κ₋₂ > 0 and κ₋₆ < 0.

Because κ₋₂ + κ₋₆ = v†Tv/‖v‖² = 0 at a defective seed (the antilinear gauge again, valid where the
eigenspace is one-dimensional so that 𝒜v = c·v fixes conj(v) = c·Tv), the two signs are one
statement, and at such a seed κ₋₂ > 0 is the *signed* form of s₆ ≠ 0: κ₋₆ < 0 ⟹ κ₋₆ ≠ 0, and the
same gauge gives |κ₋₆| = |s₆|/‖v‖², so κ₋₆ ≠ 0 ⟺ s₆ ≠ 0 *there*. The sign is strictly stronger than
that nonvanishing, stated in the *gauge-invariant* κ (sidestepping the dead "s₆ > 0", whose sign is a
branch convention). What this does, and does not, say about the β-exotic globally is spelled out
below; it is not a standalone exclusion.

Tested properly, it holds. The extraction is the whole difficulty, and it has four traps, not the two
above: a min-gap finder drifts onto the diabolic real crossings; |vᵀv| ≈ 0 certifies nothing on its
own; a stray near-real crossing masquerades as a seed unless one also demands v†Tv ≈ 0 (a genuine
seed is an 𝒜-eigenvector, so this must hold); and at N ≥ 9 the seed sits in a dense complex thicket.
Locating every real↔complex transition as a change in the complex-pair count (which diabolic
crossings do not cause), gating on the eigenvalue gap **and** the 1-dim-kernel signature
σ₂(L − λ\*I) = O(1) **and** v†Tv/‖v‖² < 10⁻⁵, and splitting by R-parity to thin the thicket (R = the
site reflection; the R-sectors are defined in the certificate section below), the law
holds at **every** defective coalescence at **N = 5, 7, 9**, both parities: the count matches the
forced net (N − 1)/2 pairs at each (2, 3, 4 pairs; 4, 7, 8 coalescences counting the re-entrant
complex→real ones), κ₋₂ ≳ 0.025 (minimum 0.0253) across all 19, and the sweep finds **zero**
semisimple (2-dim-kernel) transitions, corroborating (numerically, not a proof) the β-exotic
exclusion at N = 5, 7 by a route orthogonal to the disc-multiplicity certificate. Verifier:
`simulations/o2b_krein_sign_law.py` (self-asserting at N = 5, 7 with no argument; pass `5 7 9` to
include N = 9, or a bare `9` to run N = 9 alone).

**The gauge geometry, made explicit (2026-07-09).** The antilinear gauge conj(r) = T r puts the real
and imaginary parts of the coalescing vector on opposite classes. Write r = x + i y with x, y real;
conj(r) = T r reads x − i y = T x + i T y, hence **T x = x and T y = −y**: x is supported entirely on
class E (the +1 eigenspace of T) and y entirely on class O. Every E/O cross term then vanishes on the
disjoint supports, so rᵀr = ‖x‖² − ‖y‖² and

> κ₋₂ = (‖P₋₂ x‖² − ‖P₋₂ y‖²) / (‖x‖² + ‖y‖²),   with ‖x‖² = ‖y‖² forced by rᵀr = 0.

So the sign law κ₋₂ > 0 is exactly: the class-E real part carries more −2-rung weight than the class-O
imaginary part. This is the weight inequality below, now in closed geometric form (and κ₋₂ + κ₋₆ = 0
is just ‖x‖² = ‖y‖²).

**The mechanism, from below, and where the band-wide reading breaks (2026-07-09, corrected same day by the
five-discipline lens sweep).** The candidate above
(the seed eigenvector "inherits" E-dominance on the −2 rung from ker(K₂₂) ⊂ class E) explains the
*births* but not the whole law. Measured seed by seed, the −2-rung part of x overlaps ker(K₂₂) heavily
at the born seeds (≈ 0.90–0.95 of its weight at three of the four N = 5 seeds) but only weakly at the
*re-entrant* complex→real seeds (0.10 at the N = 5 seed q\* = 5.61, 0.14 at the N = 7 seed q\* = 11.36;
these q\* are on the unit-hop axis of the sweep, twice the octic q\* used in the trap at the end of this
section), where κ₋₂ > 0 holds nonetheless. So ker-proximity is not the universal cause: a re-entrant seed
satisfies the law with its −2-part mostly outside the kernel.

A five-discipline lens sweep (2026-07-09) then killed the natural next guess and sharpened what remains.
The guess was that κ₋₂ > 0 is a band-wide weight fact, "class-E |v|²-weight on the −2 rung exceeds class-O
weight for the real band eigenvectors with λ ∈ (−6, −2)." **That reading is false.** Sweeping every real
band eigenvector, κ₋₂ < 0 on a substantial fraction (order 20–40 % depending on the sampling),
concentrated in a band-center window (λ roughly −5 to −3.6, drifting lower at larger N) with a small tail
toward the lower edge λ = −6. The reason is a theorem the linear-algebra view
had not invoked: a *pure* bipartite hop H = [[0, B],[Bᵀ, 0]] has κ = 0 on every nonzero-energy eigenvector
(Bw = εu and Bᵀu = εw give ‖u‖² = ‖w‖²), so a pure hop carries no *net* sublattice imbalance (v†Tv = 0)
on its nonzero modes; any rung-resolved split is balanced, κ₋₂ = −κ₋₆. A net, sign-definite κ₋₂ is
therefore not something K supplies on its own: it needs the on-site D, which is also what creates the real
band (L = iqK alone is anti-Hermitian, its spectrum pure-imaginary). This is why the state-count imbalance
+(N − 1) does not imply the weight inequality, and why ker-inheritance (a zero-mode fact) was incomplete.

What survives is a *conditional*, verified exactly. The **separating law**: κ₋₂ < 0 ⟹ v†Tv < 0, zero
exceptions across the full N = 5, 7 band sweep (thousands of eigenvectors). Its contrapositive is the
**free half** of the target: every *positive-type* band eigenvector (v†Tv ≥ 0) has κ₋₂ ≥ 0. Because a real
pair leaves the axis only by a Krein collision of opposite-type branches (the pseudo-Hermitian collision
requirement, MacKay/Krein), the positive-type colliding branch has κ₋₂ ≥ 0 automatically; strict
positivity there (= s₆ ≠ 0) is *not* free, and it, together with the *negative-type* branch whose
Krein-negativity must stay confined to the −6 rung, is the whole remaining content.

At a defective seed the statement is EP-free, gauge-free, and real-symmetric. The gauge geometry above
(x = Re v ∈ class E, y = Im v ∈ class O, ‖x‖ = ‖y‖) plus the sum rule
‖P₂x‖²/‖x‖² + ‖P₂y‖²/‖y‖² = (λ + 6)/2 (a two-line theorem, not a measurement: adding
xᵀ(D − λ)x = q·xᵀKy and yᵀ(D − λ)y = −q·xᵀKy gives 4(‖P₂x‖² + ‖P₂y‖²) = (λ + 6)(‖x‖² + ‖y‖²), and
‖x‖ = ‖y‖ at a seed does the rest (the undivided identity checked to 1e-15 at all 11 seeds N = 5, 7;
the /2 form inherits the ~|vᵀv| accuracy of the located seed) collapses κ₋₂ to a
single closed form, reached independently by two of the lenses:

> κ₋₂ = ‖P₂x‖²/‖x‖² − (λ + 6)/4 = q·xᵀKy / (2‖v‖²).

So the sign law is exactly one real-symmetric inequality, **f₂ := ‖P₂x‖²/‖x‖² > (λ + 6)/4**, with x the
class-E solution of the real self-consistent eigenproblem [D − q²K(λ − D)⁻¹K] x = λx. Its named mechanism
is a **Semenoff sublattice mass**: integrating out the −6 rung dresses the −2 rung with a self-energy
Σ = −q²K₂₆(λ + 6 − iqK₆₆)⁻¹K₆₂ whose even-in-K₆₆ part is a chiral-symmetry-*breaking* on-site mass, the
only source of the −2-rung E/O asymmetry the pure-hop theorem allows. Its leading sign is a Lieb
dimensional bias: E₂ couples to the −6 *majority* sublattice O₆, O₂ to the *minority* E₆, an edge-count
asymmetry of +4(N − 1) > 0 (measured +16, +24, +32 at N = 5, 7, 9), and with λ + 6 > 0 across the band it
polarizes the in-gap vector onto E₂. Proving the mass keeps that sign after resumming the full resolvent
over λ ∈ (−6, −2), through the closed K₆₆ = −3H₃ spectrum, is the coercivity that remains. (The band-wide
numerics and the separating law are the same anchored-sweep class as the committed
`simulations/o2b_krein_sign_law.py`; the sum rule and closed form were verified independently this session.)
*[Update, later the same day: the coercivity as just stated is falsified per unit norm; the mechanism is
not attraction dominance. The mass identities survive and turn into a proximity statement, and the whole
target is re-based one section below ("Three attacks on the sign law"). The leading-order Lieb bias
+4(N − 1) remains true and remains the leading order; what fails is only its persistence after
resummation, seed by seed.]*

**What proving the sign law would, and would not, do.** It bears only on the *defective* seeds: the
identity κ₋₆ = c·s₆/‖v‖² holds where the eigenspace is one-dimensional, so κ₋₂ > 0 is the
gauge-invariant, signed reading of s₆ ≠ 0 *there*. It does **not** by itself exclude the β-exotic: at
a (hypothetical) β-exotic the eigenspace is two-dimensional, κ₋₂ is not the gauge-linked quantity, and
κ₋₂ > 0 would not certify s₆ ≠ 0. (At the defective exponent-3/2 twin the gauge relation *does* hold
and forces κ₋₂ = 0, so a measured κ₋₂ > 0 at a gate-confirmed geometric-multiplicity-1 seed certifies
√-type there.) The direct β-exotic test is the defectiveness itself, geometric
multiplicity 1 (the 1-dim-kernel signature σ₂ = O(1)), and the same sweep corroborates
it numerically at every forced seed at N = 5, 7, 9 (zero semisimple points; at N = 9 there is no
certificate to corroborate, so the sweep is the only evidence there). So the sign law is a *gauge-invariant,
signed sharpening* of the s₆ ≠ 0 reading at the defective seeds, plus a from-below mechanism candidate
for why they sit off s₆ = 0, not a standalone O2b closer (O2b: the arc ledger's label for this all-N
s₆ ≠ 0 item). Three chain lengths confirm the law; they do
not prove it (the arc's standing rule). What is new is that the missing definiteness now has a name, a
gauge-invariant form, and a sharper (EP-free, real-symmetric, seed-local) target to attack,
**f₂ := ‖P₂x‖²/‖x‖² > (λ + 6)/4**, with both the band-wide weight reading and the ker-inheritance
candidate ruled out as its route.

**Two ways it can still advance.** An exact *per-N* certificate (now taken, at N = 5, 7 and 9: the
section "The β-exotic is excluded at N = 5 and N = 7" below, and at N = 9 the gcd-certificate
section's proved layer identity), or the all-N *sign-half* prize: prove the
class-imbalance sign law, i.e. the definiteness ingredient just named (attacked and re-based the same
day: the section "Three attacks on the sign law" below, where the sign law's proved reach and its open
remainder are separated).

The reduction was put through three independent empty reviews (a mathematician, a strategy pre-mortem, an
all-N-validity audit) before landing: the sharpened reduction is their agreed survivor, the
moment-relation proof route their agreed casualty (the in-class witness above is the concrete kill).

**One trap, from the numerics above.** The sign of s₆ is **not** gauge-invariant. The colliding
eigenvector is fixed only up to scale, and the antilinear symmetry 𝒜 = T∘conj (T the metric of the
pencil section) has two branches, 𝒜r = ±r, exchanged by r ↦ i·r; s₆ flips sign between them. The Krein
index of the −6 rung, κ₋₆ := r†T P₋₆ r / ‖r‖² (P₋₆ the rung projector; the subscript is the *rung*, not
the reflection parity R), does not: it is unchanged by any rescaling r ↦ c·r. So "s₆ ∈ [0.025, 0.38]"
records the *magnitudes*, and the interval's positivity is a convention of the extraction, not content:
the same physical mode at the N = 5 seed q\* = 0.643037 gives s₆ = −0.140184 in the other branch (an
extraction-dependent number; it needs the Takagi-minimal vᵀv anchor, since a naive min-gap finder drifts
to the λ = −6 rung edge). **Only |s₆| ≠ 0, equivalently κ₋₆ ≠ 0, is the statement.** Do not try to prove
s₆ > 0.

## Three attacks on the sign law: a reduction to a regular band statement, and the O2b nonvanishing goes rational (2026-07-09)

Three independent attacks were run in parallel on the seed-local inequality f₂ > (λ + 6)/4: the
Feshbach/Semenoff coercivity route the previous section named, a Krein-branch route through the collision
itself, and an exact-certificate hunt at small N. One route died at a precisely identified inequality,
and the two others each moved the target: the sign law now reduces to a *regular, EP-free, band-wide*
conditional, and the O2b nonvanishing itself acquired an integer-polynomial form. Unless noted
otherwise, "verified" below means recomputed and asserted at all 11 defective seeds of N = 5 and N = 7
(both R-parity sectors, unit-hop axis, gauge conj(v) = Tv, v = x + iy), at the tolerances stated next. The
identities, seed values, flip witnesses, sign patterns, kernel facts and Q-pencil claims of this
section are recomputed and asserted
by the committed verifier `simulations/o2b_three_attacks_audit.py` (self-asserting at N = 5 plus the
N = 7 kernel facts with no argument; `7` adds the N = 7 seed sweep, `b7`/`b9` the N = 7 and N = 9 band
sweeps). Precision bookkeeping: the branch identities measure ~1e-14 and the cell identities ~1e-11
on the band sweeps, grid-dependent (asserted at 1e-10 and 1e-9); the seed-local and finite-difference checks run at
their methods' tolerances (the ratio identity 1e-6, κ₋₂ = −S₆/S_T 1e-4, the Q-pencil derivative 1e-3;
each stated in the verifier). Band-sweep *sample counts* are grid-dependent and quoted as measured
in-session.

**The Semenoff coercivity is falsified per unit norm; what survives is a proximity statement.** Write
R₆ := ((λ + 6)² + q²K₆₆²)⁻¹ ≻ 0 for the −6-rung resolvent kernel (R₆ is a matrix on the −6 rung, not
the spatial reflection R of the certificate section), Σ_even := −q²(λ + 6)·K₂₆R₆K₆₂ for the mass, and
W(λ, q) := |λ + 2| + Σ_even for the dressed binding on the −2 rung (class-diagonal blocks W_E, W_O), and
let α := x₂ᵀ(−Σ_even)x₂/‖x₂‖² and β := y₂ᵀ(−Σ_even)y₂/‖y₂‖² be the mass felt by the two class
components. The exact two-component system W x₂ = qK̃y₂, W y₂ = −qK̃ᵀx₂ (K̃ = K₂₂ − q²K₂₆K₆₆R₆K₆₂ the
dressed hop) gives the balance identity x₂ᵀWx₂ = q·x₂ᵀK̃y₂ = −y₂ᵀWy₂ and the **ratio identity**

> ‖x₂‖²/‖y₂‖² = (β − |λ+2|) / (|λ+2| − α),

so f₂ > (λ + 6)/4 ⟺ **||λ+2| − α| < ||λ+2| − β|**: the class-E mass value sits *closer to the binding
depth* than the class-O one. That is a proximity statement, not a dominance statement, and dominance is
what fails: at 3 of the 11 seeds the binding order flips to α > |λ+2| > β (witnesses, unit-hop axis:
N = 5 q\* = 1.286074 with α = 1.8622, |λ+2| = 1.8196, β = 1.7474; N = 7 q\* = 1.077392 and 1.447107),
yet κ₋₂ > 0 holds through the other branch of the ratio identity. Both W_E and W_O are indefinite at
every seed, so no per-class definiteness exists, and the +4(N − 1) trace bias does not act pointwise.
The surviving margins are thin (0.019 vs 0.027 at the N = 7 seed q\* = 1.077392), so any proof of the
proximity form must be exact, not perturbative.

**The branch identity, and the reduction that is the section's centerpiece.** For *every* simple real
band eigenvector v of L(q) (not only at seeds), with σ := v†Tv/‖v‖² its normalized Krein form value
(its sign is the branch's Krein type) and λ̇ := dλ/dq the
branch slope (exact from the two real equations (D − λ)x = qKy, (D − λ)y = −qKx together with the slope
identity λ̇ = −2·xᵀKy/vᵀv; found independently by two of the attacks in algebraically identical forms):

> **κ₋₂ = σ·[(λ + 6) − q·λ̇]/4,  κ₋₆ = σ·[q·λ̇ − (λ + 2)]/4.**

Now let a defective seed be approached along its two colliding real branches: at an *exit* (this arc's
"birth": the two real eigenvalues leave as a complex-conjugate pair, the count-drop)
λ± = λ\* ± μ·(q\* − q)^{1/2} + O(|q − q\*|), and at a re-entrant *entry* (a "death" of the complex pair)
the same with (q − q\*)^{1/2}. Mind what H1 does **not** give: under H1 the two branches carry odd
half-powers, λ± = λ\* ± c₁δ^{1/2} ± c₃δ^{3/2} + …, and a count-drop forces only *some* odd coefficient
nonzero, not the leading one. Call a seed **√-type** when μ = c₁ ≠ 0; since μ² = rᵀCr/rᵀr₁ =
(4/q\*)·s₆/rᵀr₁ (rᵀr₁ ≠ 0 being exactly H1), √-type ⟺ **s₆ ≠ 0** there; so demanding μ ≠ 0 would
silently assume the very nonvanishing O2b must prove (a defective exponent-3/2 point, disc order 3,
is H1-compatible; per N it is excluded only by the disc-multiplicity certificate, which kills every
order-≥3 locus). At a √-type seed the slopes q·λ̇ diverge to ±∞ while κ₋₂ (|κ₋₂| ≤ 1 by definition)
converges to a common value κ₋₂\* on both branches, the eigenvectors converging at a
geometric-multiplicity-1 point, so the divergent factors in the identity must cancel. The
**separating law** of the previous
section (κ₋₂ < 0 ⟹ v†Tv < 0 for simple real band eigenvectors) enters through its contrapositive: on a
positive-type branch κ₋₂ ≥ 0, so q·λ̇ ≤ λ + 6 < 4 is bounded above. The branch whose slope diverges to
+∞ therefore cannot be positive-type, so the positive-type branch takes the −∞ slope, falling into
every exit and away from every entry, and the identity evaluates to κ₋₂\* = q\*·|μ|·|s|/8 =
|s₆|/‖r‖² > 0, where s is the √-coefficient of the Krein form along the branch,
σ = s·|q − q\*|^{1/2} + O(|q − q\*|), with s = ±2μ·rᵀr₁/‖r‖² (r₁ the Jordan generalized vector); the two
colliding branches carry opposite Krein types automatically (σ± = ±s·|q − q\*|^{1/2} to leading order).
At a defective seed that is *not* √-type the same limit gives κ₋₂\* = 0 consistently. So the separating
law supplies exactly the **sign** (it forbids the branch κ₋₂\* = −|s₆|/‖r‖²), and the **strictness is
the nonvanishing s₆ ≠ 0 itself**, which nothing in this argument supplies: that is target 2 below.
One scope premise, flagged rather than proven: the reduction needs λ\* strictly inside (−6, −2) so the
colliding branches are band eigenvectors within the separating law's stated scope; the Hermitian
relation (λ + 2)‖r₂‖² + (λ + 6)‖r₆‖² = 0 gives that whenever both rung components are nonzero, a
boundary seed would need single-rung support (r₂ ∈ ker(K₆₂) or r₆ ∈ ker(K₂₆)), and none of the 19
measured seeds is boundary.
Hence, modulo H1 (the standing premise, untouched):

> **The band-wide separating law ⟹ κ₋₂ ≥ 0 at every isolated defective seed, with strict positivity
> exactly at the √-type seeds (s₆ ≠ 0), for every odd N.**

(Isolated: one colliding pair at q\*. Two coincident EP2s at one q\*, the drop-by-4 case, put geometric
multiplicity 2 at the point and the two-branch argument does not apply as written; the sharpened
section defers that case to the certificate, and so does this reduction.)

(All 19 defective seeds measured at N = 5, 7, 9 are √-type, κ₋₂ ≥ 0.025.)

The *sign half* of the seed-local inequality f₂ > (λ + 6)/4 is thereby upgraded to a *regular*
statement with no coalescence, no gauge, and no exceptional point in it; the nonvanishing half stays
with s₆. The deflation of the previous section carries
over verbatim: this bears on the defective seeds only (at a hypothetical β-exotic the argument yields
κ₋₂\* = 0, consistent with the biconditional: σ → 0 on both branches there, the branch eigenvectors
converging to a common limit in the β normal form), so it is a route to the sign law's sign, not by
itself an O2b closer. The separating law itself is now verified at **N = 5, 7, 9** (0 violations in 1290, 2012,
2654 simple real band samples, both R-parity sectors, q ∈ (0, 20]; a different sweep than the previous
section's "thousands of eigenvectors" (5110, as the ledger records), hence the different counts, and
the committed verifier's band legs sweep the full block with no R-split; the
minima of κ₋₂ over positive-type
samples were 0.085, 0.037, 0.051; sample minima, not band minima: a denser N = 5 sweep reaches 0.073),
three chain lengths and not a proof. Its sharpest *sufficient* form, where any proof must
consume the class imbalance: with a := ‖(Kx)₂‖², b := ‖(Kx)₆‖², a′ := ‖(Ky)₂‖², b′ := ‖(Ky)₆‖² (the
images in the four class-rung cells; Kx ∈ O, Ky ∈ E), one has ‖x₂‖² − ‖y₂‖² = q²(a′ − a)/(λ + 2)² and
σ‖v‖² = q²[(a′ − a)/(λ + 2)² + (b′ − b)/(λ + 6)²], and the candidate reads **b′ ≥ b ⟹ a′ ≥ a**
(equivalently κ₋₆ ≥ 0 ⟹ κ₋₂ ≥ 0). Mind the logic: this is *strictly stronger* than the separating law,
not equivalent to it (the law permits a′ < a with b′ ≥ b whenever the (λ + 2)⁻²-weighted deficit
outweighs the (λ + 6)⁻²-weighted surplus, making σ‖v‖² negative), but it implies the law, it holds
with zero violations on the same sweeps
(and on an independent audit sweep), and proving it closes the sign law just the same. The E↔O-swapped
statement is false (a seed has σ = 0 with ‖x₂‖ > ‖y₂‖), so the imbalance data must enter; newly measured
structure for that purpose (measured N = 5, 7; the −6-rung and whole-K nullities are the resonance-free
*baseline*, not an all-N law (this note's own Piece-2 count is 21 ≠ 15 at N = 11), while
nullity(K₂₂) = N − 1 is the Piece-1 theorem, its class-E localization the endpoints-in-E add-on
asserted at the swept N): ker(K₆₆) lies entirely in class O with baseline
nullity 3(N − 1)/2, and ker(K) has baseline nullity 3(N − 1)/2 splitting as (N − 1)/2 in E plus
N − 1 in O.

**The certificate goes rational: κ₋₂ as a ratio of two integer polynomials.** At every defective seed
adj(λI − L) = c·rrᵀ with c ≠ 0 (rank-1 adjugate at geometric multiplicity 1), and in the antilinear
gauge r†r = rᵀTr, whence

> **κ₋₂ = −S₆(λ, q)/S_T(λ, q)**,  S₆ := tr(P₋₆·adj(λI − L)),  S_T := tr(T·adj(λI − L)),

both *real polynomials in (λ, q²) with integer coefficients* (realness and evenness are a two-line
consequence of T-pseudo-reality and the q ↦ −q conjugation, machine-checked at generic points;
integrality holds because D, K, T, P₋₆ have integer entries; the *value* κ₋₂ = −S₆/S_T is asserted at
all 11 seeds by the verifier's SVD-product route, and an independent in-session resolvent-limit check
agreed at the four N = 5 seeds). The consequences are
locus-independent. First, **S₆ ≠ 0 at a count-drop point ⟹ geometric multiplicity 1 (the β-exotic is
excluded there) and s₆ ≠ 0**: S₆ ≠ 0 forces adj ≠ 0, hence a 1-dimensional kernel, hence c ≠ 0 and
s₆ = S₆/c ≠ 0. The whole O2b target becomes the nonvanishing of one integer polynomial on the seed
variety (the H1 caveat stands exactly as before: this does not bound the algebraic multiplicity).
Second, the sign law itself reads **S₆·S_T < 0 on the seed locus** (verified at all 11 seeds; the
*product* is convention-independent, the individual signs are not: on the full block, seeds sorted by
q\*, the N = 5 signs are (+,−,−,+) against (−,+,+,−), asserted by the verifier, while per-R-sector
adjugates give (+,+,−,+) against (−,−,+,−), recorded but not asserted; the co-sector characteristic
factor flips both together). A word of caution on
per-N use: the locus S₆ = 0 necessarily *contains* every point of geometric multiplicity ≥ 2 (adj ≡ 0
there), i.e. all diabolic crossings **and** any hypothetical β-exotic, so an elimination certificate on
{F_res = ∂_Λ F_res = S₆ = 0, q real > 0} cannot by itself tell the two apart (F_res, ∂_Λ and the
R-sectors are defined in the certificate section, next); the honest per-N
exclusion remains the disc-multiplicity certificate of that next section. **(This caution is
retired 2026-07-10 per-N: the gcd certificate section further below separates the diabolics onto
the doubled disc layer A₂ and proves geometric multiplicity 1 directly on the simple layer,
exactly at N = 5 and, at the time of writing, modulo two named premises at N = 7, both
discharged 2026-07-16, so exactly at N = 7 too; kept here so the paragraph reads as it
stood when written.)** What the rational form
adds is the *all-N shape*: one integer polynomial whose nonvanishing at the count-dropping loci is
O2b, handing the question the vertex-deleted characteristic polynomials of the full pencil L (whose
diagonal blocks are the closed K₂₂-path and H₃ structures) in place of an eigenvector. A supporting
hint, a search rather than an exclusion: a PSLQ hunt on Newton-refined N = 5 loci (46 digits) surfaced no
low-complexity minimal polynomial for q\*², λ\*, κ₋₂; the values are algebraic *by construction* (degree
bounded by deg disc), and 46 digits cannot certify absence anywhere near the searched degree-24,
height-10¹⁴ box, so this says only that no low-complexity closed form presented itself, consistent with
this arc's forced/free law; a certificate must be identical in (λ, q), which is what the rational form
provides.

**A second gauge-free equivalent, banked for a global attack (cashed and bounded the same arc, next
paragraph).** The Hermitian **Krein pencil**
Q(λ, q) := T(D − λ) + iq·TK (Hermitian for real λ, q since TK is real antisymmetric) has real
eigencurves ν(λ); real eigenvalues of L are their zeros, the Krein sign of a branch is −ν′(λ), and a
seed is a *tangential, T-neutral touching* ν = ν′ = 0. Exactly there ∂ν/∂q = −(4/q)·κ₋₂ (verified at
all 11 seeds), so the sign law is: **every T-neutral touching of 0 by an eigencurve of Q descends in
q.** The exits (complex pair born, the count-drops) open downward (ν″ := ∂²ν/∂λ² ∈ [−0.61, −0.35]) and
the re-entrant entries open upward (ν″ ∈ [+2.56, +3.72]) at the 11 seeds, exactly the pattern the
descending law forces: a downward-opening maximum descending through 0 loses its two real zeros (an
exit), an upward-opening minimum descending through 0 gains two (an entry). This is an equivalence, not a
proof, and a counting argument alone cannot close it (nothing in the net count forbids a compensated
ascending/descending pair); its value is that the target lives on one Hermitian family with no
eigenvector-weight language left in it.

**The banked handle, cashed and honestly bounded (2026-07-16 night): ν″ is the sign
characteristic, and the sign characteristic is blind to O2b.** The touching curvature has a closed
second-order (Rellich) form: with v the unit kernel vector of L − λ\* and w the minimum-norm
generalized vector ((L − λ\*)w = v), Q′ = −T and Q″ = 0 give ν″(λ\*) = −2·v†Tw/‖v‖². The number
v†Tw is real, and it is exactly the **Gohberg-Lancaster-Rodman sign characteristic** of the size-2
Jordan block of the T-selfadjoint L (the root-subspace Gram [[v†Tv, v†Tw], [w†Tv, w†Tw]] has
determinant −|v†Tw|², and the GLR canonical form makes that Gram nondegenerate for a real
eigenvalue). Consequence, all N: **at any genuine size-2 Jordan block the λ-touching is quadratic,
never flatter, forced by structure** (a semisimple β-exotic would instead make ν = 0 a double
eigenvalue of Q, and any higher Jordan structure would flatten the touching, ν″ = 0). Read as a
certificate it runs the other way: computing ν″ ≠ 0 together with ν = 0 simple and the T-neutral
touching certifies "Jordan size exactly 2" at that seed, from the local pair (v, w) alone, cheaper
than a full spectral decomposition (still a per-seed computation; the all-N content is the forcing
direction, not a proof of H1). What ν″ does **not**
certify is s₆: the O2b scalar is the orthogonal jet ∂ν/∂q = i·v†TKv/‖v‖² = −(4/q\*)·κ₋₂ (the
relation above, now verified three ways at every N = 5 seed: formula, finite difference, and the
banked identity; the first two agree exactly and to 1e−10, the finite differences to 1e−4), and
the two jets are independent, ν″/κ₋₂ spanning −4.33 to +40.82
with a sign flip across the four N = 5 seeds. The moment-relation kill witness of the "Why the
natural proof route does not close it" section separates them outright: it is a genuine T-selfadjoint size-2 Jordan block with
ν″ = −0.5 ≠ 0 and s₆ = κ₋₂ = ∂ν/∂q = 0 exactly. So the Krein pencil *restates* κ₋₂ as the q-jet
of the touching and supplies no new lever for the nonvanishing; the measured exit/entry ν″ signs
are a defectiveness/sign-characteristic diagnostic, not evidence toward s₆ ≠ 0. The route "prove
s₆ ≠ 0 for all N via GLR sign-characteristic theory" joins the do-not-retry list with that precise
reason. Gate: [the Krein touching dictionary](../simulations/o2b_krein_touching_dictionary.py)
(T-selfadjointness and Hermiticity exact 0.0, ν = 0 simple, T-neutrality < 1e−13, the ν″ and
∂ν/∂q formulas against finite differences at all four forced N = 5 seeds, the ratio spread, and
the witness, all asserted).

**New dead ends, measured (these join the do-not-retry list).** Per-unit-norm mass dominance (the three
flipped seeds above). Band-wide monotonicity σ·λ̇ < 0, equivalently xᵀKy > 0 band-wide: false at
452/2310 samples (N = 5) and 860/3622 (N = 7), with σ among violators reaching ≈ −0.1 (sweep-dependent),
so no uniform σ-gap exists and the quantitative strengthening
"κ₋₂ < 0 ⟹ σ ≤ −c with an N-free constant c > 0" is dead too. sign(S₆) alone as an invariant: false,
the rank-1 factor c changes sign from seed to seed. Sign-definiteness of S₆ on the strip
λ ∈ (−6, −2), q > 0: false (mixed signs on a 60×60 grid in-session; the verifier's grid is 30×30). The
band-wide J-cone law (J := iTK Hermitian,
"v†Jv < 0 for all real band eigenvectors"): false at 111/600 (N = 5) and 166/864 (N = 7), all violations
on T-negative vectors; it is collision-local only (pointwise, v†Jv < 0 is the *same* condition as
σ·λ̇ < 0, so these two dead ends are one falsification counted on two in-session sweeps, not two
independent ones). The laws and identities above are asserted by the committed
`simulations/o2b_three_attacks_audit.py`; the violation counts are grid statistics of the in-session
sweeps.

**Where this leaves O2b.** Two named targets, both sharper than what this section started with: prove
the separating law (or its stronger cell form b′ ≥ b ⟹ a′ ≥ a; either closes the *sign half* (κ₋₂
never negative, strict wherever s₆ ≠ 0) at all defective seeds, all odd N, modulo H1), and prove
S₆ ≠ 0 at the count-dropping loci (the O2b nonvanishing itself, which is also what would make the sign
strict everywhere; mind the per-N caution above: the S₆ = 0 locus contains the diabolics, so this is
an all-N shape, not a ready per-N certificate; **that caution is retired 2026-07-10**: the gcd
certificate section below separates the diabolics onto A₂ and delivers the per-N certificate at
N = 5 and 7). The two are complementary and only together complete
κ₋₂ > 0: the first is the sign, the second the exclusion.

## The β-exotic is excluded at N = 5 and N = 7, both R-parities (2026-07-09)

The per-N certificate is taken, and it is cheaper than the grind it was expected to be: it reads a
**multiplicity**, not a Galois group. It settles one chain length at a time, and so far it has
settled two.

**Three names this section needs** (the "Three attacks" section above already borrowed F_res and the
R-sectors; their definitions live here). The spatial reflection R (site i ↦ N+1−i)
commutes with both A (the dephasing diagonal, written D in the two sections above) and C, so it splits
the block into two invariant sectors, **R-even** and **R-odd**
(dimensions 26 and 24 at N = 5). Within each sector the characteristic polynomial factors as
χ = **AT** · **F_res**: the **AT factor** (AT = Absorption Theorem) collects the rate-confined strands,
whose eigenvalues run q-linearly, λ = r₀ + iκq, with the decay rate locked to one value of the dephasing
diagonal; they are the Absorption-Theorem-locked modes of
[F89_PATH_K_DIABOLIC.md](F89_PATH_K_DIABOLIC.md), and **F_res**, the *residual* factor, carries
everything else, degree 18 (R-even) and 17 (R-odd) at N = 5. Λ denotes the spectral variable of these
polynomials, λ an individual eigenvalue branch. Discarding AT loses nothing, for a reason worth stating
once: the AT subspace is *q-independent* and invariant for every q, so A and K each preserve it, and
each being Hermitian preserves its orthogonal complement too. The splitting is therefore
block-**diagonal**, not merely block-triangular, so no Jordan chain can run across the seam between AT
and F_res where disc_Λ(F_res) would fail to see it. (Verified from below at N = 5: all four forced
seeds lie wholly in F_res, two in each R-parity sector.)

**What the discriminant's multiplicity knows.** At a branch locus q\* of F_res, the order of vanishing
of disc_Λ(F_res) reads the Puiseux exponent of the branches that meet there:

| local structure at q\* | branches | ord_{q\*} disc |
|---|---|---|
| √-type defective EP2 (the seeds) | λ± = λ₀ ± c·(q − q\*)^{1/2} | 1 |
| diabolic (semisimple) crossing | analytic, cross transversally | 2 |
| cubic branch point (3×3 Jordan) | three sheets, exponent ⅓ | 2 |
| **the β-exotic** (semisimple) | λ± = ±s^{3/2} | **3** |
| **its defective exponent-3/2 twin** (geom mult 1) | λ± = λ₀ ± c·(q − q\*)^{3/2} | **3** |

A 2-cycle with exponent e gives (λ₁ − λ₂)² ~ (q − q\*)^{2e}, so e = ½ ⟹ 1 and e = 3/2 ⟹ 3. The cubic
row is the Kato count of the Lemma above, not this formula: a k-cycle branch point makes the discriminant
vanish to order k − 1, so a 3-cycle gives 2. Every
*other* colliding pair at the same q\* contributes a **non-negative** order, so a coincidental extra
collision can only raise the total (3 + 1, 3 + 2), never mask it. Hence:

> **disc_Λ(F_res) has no root of multiplicity ≥ 3 off q = 0  ⟹  no count-dropping order-≥3 structure
> on that block: no β-exotic and no defective exponent-3/2 twin** (count-drops provably live in F_res,
> two paragraphs below; disc_Λ(F_res) is blind to AT-internal degeneracies, which cannot drop the count).

This is deliberately weaker than squarefreeness: the diabolic loci are genuine *double* roots and must
survive. That weakness is exactly why the dead Galois route is not needed here. A Galois group sees
which sheets swap; it never sees the ½-versus-3/2 exponent. A multiplicity does.

**The certificate.** `FoldResultantCertificate` (built for a different statement, the remainder-R1
gcd) already computes D(q) = disc_Λ(F_res)(q) mod p by interpolation (this D is the engine's name for
the discriminant, not the dephasing diagonal and not Piece 3's degenerate-resonance count) and splits
it into squarefree
layers. The certified readings, both parities and both chain lengths:

| N | sector | residual degree | layers | deg_q D (bound) | v_q(D) | primes / bound | max multiplicity |
|---|---|---|---|---|---|---|---|
| 5 | R-odd | 17 | [56, 26] | 246 (254) | 138 | 139 / 138 | **2** |
| 5 | R-even | 18 | [56, 32] | 274 (286) | 154 | 156 / 155 | **2** |
| 7 | R-odd | 52 | [222, 390] | 2508 (2556) | 1506 | 1327 / 1326 | **2** |
| 7 | R-even | 53 | [228, 420] | 2612 (2672) | 1544 | 1377 / 1376 | **2** |

Two layers everywhere, never a third. The simple roots are the √-branch (defective) loci; the forced
real seeds are among them, together with the complex defective loci this note does not track. The double
roots are *not* identified here: by the table above, an order-2 zero is a diabolic crossing, an
analytic-defective 2×2, a cubic branch point, or two coincident defective pairs. The theorem does not
need to know which. It needs only that no layer of multiplicity 3 exists.

The lift to an exact statement over ℚ(i) is **one-way**: reduction modulo a prime ideal π_p that does
not divide lc_q(D) is a degree-preserving ring homomorphism, so a factor (q − α)^m survives and
distinct roots can only *merge*. Therefore max-mult(D mod p) ≥ max-mult(D), and reading multiplicity 2
at a single certified prime proves max-mult(D) ≤ 2 exactly. The prime must be good at **both ends of the
q-axis**: attaining the true deg_q D (no root escaped to q = ∞) and the true q-valuation (no nonzero root
collapsed onto q = 0, where it would be stripped away with the q-power). A Hadamard bound on ‖D‖, the
same device the engine already uses for the resultant's degree, bounds how many of the sampled primes
can lose the degree, and separately how many can raise the valuation; sampling past that bound certifies
each of the two true values. It does *not*, by itself, promise that one prime attains both. That prime
is searched for among the primes the run samples anyway (139 / 156 at N = 5 on the D-only path,
1327 / 1377 at N = 7; the earlier full-path N = 5 run sampled 242 / 256, which is why the arc ledger
carries those larger counts), and
the certificate **fails closed** if none does. In all four cases the very first sampled prime attains both,
and further primes reproduce the same layer degrees exactly (a second and an eighth at N = 5, a fourth
at N = 7; the gates `TheLayerReading_DoesNotDependOnWhichPrime` and its N = 7 sibling), so the verdict
does not hang on one reduction.

**What it closes.** The β-exotic, and with it the defective exponent-3/2 twin (the whole order-≥3
class), is excluded at N = 5 and at N = 7, on both R-parity sectors,
unconditionally: the multiplicity bound is the whole argument, and it needs nothing else.

**What it closes as a corollary: H1.** The subsidiary premise H1 (each forced seed has algebraic
multiplicity exactly 2, no 3×3 Jordan) follows wherever the certificate runs, by **elimination, not by
counting the order of the zero**. Do not run the tempting shortcut "a count-drop changes the discriminant's sign, so the zero
has odd order, so with max multiplicity 2 the order is 1". That step silently assumes exactly one
conjugate pair is born at q\*; two forced seeds coinciding at one q\* would drop the real count by 4 at
an order-2, sign-preserving zero, and the shortcut breaks. The argument that holds needs no order parity:

Maximum multiplicity 2 already excludes the β-exotic and its defective exponent-3/2 twin (ord 3 each),
every Jordan block of size ≥ 4 (ord ≥ 3),
and every semisimple degeneracy of three or more branches (ord ≥ 3). What can still sit at a branch
locus is a √-type defective EP2 (ord 1), a diabolic crossing (ord 2), an analytic-defective 2×2 (ord 2:
a Jordan block whose two branches stay analytic, as in [[0,1],[0,s]]), a cubic branch point (ord 2), or
a coincidence of these. Of that list, **only the √-type EP2 changes the real-eigenvalue count**: the two analytic
cases keep the real-count unchanged across q\* (either two real branches throughout, or a
complex-conjugate pair that only touches the axis at q\*), and a cubic branch point keeps one real branch
and one conjugate pair on both sides of q\*. Better still, the bound forbids the coincidences that would muddy this: an EP2 sharing
its locus with any of the ord-2 structures would give ord ≥ 3. So a locus that carries a count-drop
carries √-type EP2s and **nothing else** (one, or two coincident), each a defective 2×2 Jordan block,
algebraic multiplicity exactly 2. That is H1.

One thing this needs is that the count-drops live in F_res at all, and not in the AT factor where
disc_Λ(F_res) is blind. They do, for every N and with no check: an AT strand is q-linear with purely
imaginary slope, λ = r₀ + iκq, so it is real at every q > 0 when κ = 0 and complex at every q > 0 when
κ ≠ 0. The AT factor's real count is constant on q > 0. It cannot drop. Every count-drop is a branch
locus of F_res.

The other thing it needs is that F_res has **real** coefficients, so that "real branch" and
"conjugate pair" mean anything. Two steps, and the second is the one to watch. The first sharpens the
pencil section's statement, and the two versions are easy to confuse: **whole-block** self-conjugacy
holds at every N, but F_res's roots come from a single R-parity sector, so what is needed is **per-sector**
self-conjugacy. That is what odd N buys. At odd N the bipartite
metric T commutes with the reflection R (the site-sum parity is preserved because 3(N+1) is even), so T
restricts to each R-parity sector and T L T = L† holds there: **each sector's spectrum is
self-conjugate**. The AT spectrum inside that sector is separately self-conjugate, because the AT strand
slopes κ come in ± pairs (the chiral pairing λ_{N+1−k} = −λ_k, `ChiralKClaim`). Hence the roots of
F_res, the sector's spectrum with AT's removed, are self-conjugate, and the monic F_res is real. Both
steps are checked from below at N = 5 and N = 7, and at N = 9 by the gcd certificate's runtime asserts
(AT is asserted even in u with integer coefficients, exactly the ± slope pairing, and the exact
division leaves F_res monic in λ with integer coefficients, so F_res is real by construction there);
the second step rests on the chiral pairing being what fixes
the AT slopes, which we have not derived in general.

**How N = 7 was reached, and why not by "the same call".** `CertifyComplete` proves the remainder-R1 gcd
as well, and R1's resultant runs against the corner block (p_c+1, p_c+1), whose dimension
C(N, (N+1)/2+1)² is 25 at N = 5 but **441** at N = 7: deg_q R jumps from 422 to roughly 53·441 ≈ 23000
interpolation nodes per prime, times the primes the Mignotte lift needs. A first N = 7 probe ran past 28
minutes without finishing. The β-exclusion requires none of that; it reads only D, whose size is
about resDeg·(resDeg − 1) ≈ 2600 (that crude product slightly overestimates; the table's "(bound)"
column records the engine's exact per-run degree bounds, 2556/2672 at N = 7, which the fail-closed
logic uses). So N = 7 runs through a separate entry point, `CertifyDiscMultiplicity`,
which drops the corner, the resultant and the Mignotte lift (that bound certifies the *gcd*, a statement
this one does not make). It is pinned where both can run: at N = 5 the two paths agree number for
number, and the D-only one is about ten times faster. N = 7 then costs 2.4 and 2.7 minutes per parity.

**What it does not close.** Not N = 9: the block there is 324-dimensional, and its exact bivariate
ℤ[i][q] charpoly is a different engineering problem, not a longer run of this one (the gcd
certificate section below reaches N = 9 by a different route, the proved layer identity,
2026-07-17). And not any N at
once. This is a per-N certificate, not a law; it retires chain lengths one at a time. The all-N item,
s₆ ≠ 0, stays open, and nothing here bears on it.

Verified from below at all four N = 5 seeds: each lies wholly in F_res (no AT component), is defective,
and has splitting exponent 0.500.

Gates. Fast (Categories `FOLDRESULTANT` and `DISCMULT`, seconds):
`DiscHasNoMultiplicityThreeRoot_ExcludesTheBetaExotic`, the agreement test
`N5_DiscOnlyPath_ReproducesTheFullCertificate_ExactlyAndFaster`, and the independent-prime falsifier
`TheLayerReading_DoesNotDependOnWhichPrime`. Slow (Category `SLOW_DISCMULT_N7`, about five minutes):
`N7_DiscOnlyPath_ExcludesTheBetaExotic` and its two-prime falsifier. Typed as
`BetaExoticPerNExclusionClaim` with the live witness `inspect --root betaexotic` (N = 5 by default;
`--N 7` pays the minutes).

## The gcd certificate: the per-N nonvanishing closes, exactly at N = 5, 7 AND 9 (2026-07-10; both N = 7 premises discharged 2026-07-16, see the two subsections inside; N = 9 landed 2026-07-17)

The disc-multiplicity certificate above bounds an *exponent*: after it, no count-drop at N = 5 or
N = 7 can be an order-≥3 point. What it never said is the O2b scalar itself, s₆ ≠ 0 at the seeds; a
multiplicity bound reads the branching, not the weight. That gap closes now, rationally, by one gcd
per sector. The computation was then adversarially reconfirmed by an independent route (independent
pencil builders, Faddeev-LeVerrier adjugate recursion in place of vertex-deleted characteristic
polynomials, a second factorization library; zero discrepancies at N = 5), and extended to N = 7 by
the mod-p/CRT engine described below. On 2026-07-17 the same engine, with the disc sweeps fanned
out over a multiprocessing pool, extended everything to N = 9 (sector dims 164/160, F_res degrees
116/115); at N = 9 the layer discharge is not only the premise-killer but the *only* route to max
multiplicity 2, since the DISCMULT engine never reached N = 9.

Everything lives per R-sector over ℤ[λ, w], w = q² (this note's unit-hop axis), with the "Three
attacks" section's two integer polynomials S₆ = tr(P₋₆·adj(λI − L)) and S_T = tr(T·adj(λI − L)).
Four facts per sector, at each of the three N (grades per fact and per N in the bookkeeping paragraph below):

1. **S₆ splits as strands × one core.** S₆ = (integer content) · (prefactors) · G with G irreducible
   over ℚ, and every prefactor either an AT-strand factor or a linear rung factor (λ+2), (λ+6). At
   N = 5: S₆ᴱ = 4·(w + (λ+2)²)·Gᴱ with Gᴱ of degree (23, 11); S₆ᴼ = 2·(λ+2)(λ+6)·(3w + (λ+2)²)·Gᴼ,
   Gᴼ of degree (19, 8). At N = 7: S₆ᴱ = 2·(λ+2)(λ+6)²·N₂·N₂′·Gᴱ, Gᴱ of degree (67, 31);
   S₆ᴼ = N₄·N₄′·Gᴼ, Gᴼ of degree (63, 31) (Nₖ the AT strand factors of that sector). At N = 9:
   S₆ᴱ = 4·N₄·N₄′²·Gᴱ, Gᴱ of degree (151, 75); S₆ᴼ = 4·(λ+2)(λ+6)³·N₄·N₄′²·Gᴼ, Gᴼ of degree
   (143, 68) (both N₄, N₄′ of degree (4, 2), the second squared). Each prefactor
   is **sign-fixed on the strip** λ ∈ (−6, −2), w > 0, and provably so by shape: every nonlinear AT
   factor is a Galois product of ((λ − r₀)² + κ²w) terms with κ real (spectrum of a symmetric
   restriction of K), strictly positive for w > 0, and the linear factors do not vanish on the open
   strip. So on the strip, S₆ ≠ 0 ⟺ G ≠ 0: the whole nonvanishing lives in one irreducible bivariate
   polynomial per sector.
2. **The discriminant layers, exactly.** disc_Λ(F_res) = c · w^v · A₁ · A₂² per sector (c an integer
   constant; A₁, A₂ primitive with positive leading coefficient), with **A₁ (the simple layer)
   irreducible over ℚ** and A₂ irreducible (A₂'s irreducibility exact at N = 5, evidence-grade at
   N = 7 and 9; nothing below uses it), and max multiplicity 2 off w = 0. Degrees
   in w, slashes in **E/O order** (mind: the DISCMULT table in the certificate section lists R-odd
   first): at N = 5, A₁ has degree 28 in both sectors (A₂: 16/13, v: 77/69); at N = 7, A₁ has degree
   114/111 (A₂: 210/195, v: 772/753); at N = 9, A₁ has degree 288/286 (A₂: 1072/1021, v: 4054/4020).
   The N = 5 computation is exact over ℤ, which upgrades the
   certificate section's fail-closed mod-p reading of "max multiplicity 2" to an exact equality
   there; the N = 7 layer degrees reproduce the DISCMULT-certified table after the unit change
   (that table is in q, this list in w = q², so every degree doubles: 2·114 = 228, 2·210 = 420,
   2·772 = 1544). At N = 9 there is no DISCMULT table to reproduce (that engine never reached
   N = 9); the proved identity itself delivers max multiplicity 2 there, see the second
   premise-discharge subsection.
3. **The gcds.** gcd(Res_Λ(F_res, S₆), A₁) = 1 in both sectors at all three N. And, closing a gap the
   first pass had not addressed: the **cross-sector** certificate gcd(Res_Λ(F_res_s, χ_other), A₁_s) = 1
   in both directions, at all three N (a λ\* shared with the *other* sector's spectrum would give
   full-block geometric multiplicity 2 despite sector S₆ ≠ 0; it cannot happen on the simple layer).
4. **No coincident pair.** The one order-2 *count-dropping* configuration the certificate section
   keeps alive, two coincident √-type EP2s at a single w\* (distinct λ\*, dropping the count by 4),
   would sit on A₂, where fact 3 says nothing. It is excluded by a subresultant certificate: two
   double λ-roots at one w₀ force psc₀(w₀) = psc₁(w₀) = 0 (psc₁ the first principal subresultant
   coefficient of (F_res, ∂_λF_res); psc₀ is the discriminant), and
   **gcd(disc_Λ(F_res)/w^v, psc₁/w^{v′}) = 1** in both sectors at all three N, so no w₀ ≠ 0, real or
   complex, carries more than one double λ-root. This fact is free of the layer premise at
   every N (psc₁ is a fixed integer determinant in F_res's coefficients, so its mod-p reduction
   commutes unconditionally, and the gcd lifts from a good prime attaining the discriminant's
   certified degree); at N = 5 the exact inventory says the same
   thing concretely (A₂'s real positive roots: none in R-even, exactly the known diabolic in R-odd,
   carrying exactly one real double λ).

The chain from the four facts to the theorem, each link named. F_res is **monic in λ** (checked),
so the resultant specializes exactly at every w° (no leading-coefficient degeneration):
Res_Λ(F_res, S₆)(w°) = ∏ᵢ S₆(λᵢ(w°), w°) over the F_res roots. A count-drop is a real double root of
F_res at w\* > 0. By the certificate section's taxonomy the non-coincident order-2 classes
(diabolic, analytic-defective 2×2, cubic branch point) drop no count, the coincident-EP2 pair is
excluded by fact 4, and order ≥ 3 is excluded by fact 2's max multiplicity 2; what remains has
discriminant order 1. So **every count-drop sits on the simple layer A₁** (and
w\* ≠ 0, A₁(0) ≠ 0 checked). Fact 3 then forbids Res_Λ(F_res, S₆)(w\*) = 0, hence S₆(λ\*, w\*) ≠ 0,
hence adj(λ\*I − L) ≠ 0: **sector geometric multiplicity 1**; the cross-sector gcd extends that to
the full block. At geometric multiplicity 1 the adjugate is rank one, adj = c·rrᵀ with c ≠ 0 in the
antilinear gauge, so s₆ = S₆/c ≠ 0. That rank-one gauge step is the single non-rational link in the
chain (standard linear algebra at geometric multiplicity 1; corroborated numerically at all 19
seeds, where κ₋₂ = −S₆/S_T from the exact polynomials matches the eigenvector-side values to the
1e−5 gate, printed at six decimals: N = 5: 0.314463, 0.140184, 0.351821, 0.062939; N = 7: 0.090208, 0.238066, 0.184075,
0.060903, 0.025288, 0.377212, 0.050276; N = 9: 0.064236, 0.190886, 0.146061, 0.090363, 0.139896,
0.379238, 0.049866, 0.040848 (the full defective census at N = 9: 8 coalescences, 4 per R-parity,
2 of them re-entrant complex→real ones; one convention for all three N: the per-N seed lists count
ALL defective coalescences including the re-entrant ones, 4/7/8 at N = 5/7/9 with 1/2/2 re-entrant,
while the forced NET drop is (N − 1)/2 conjugate pairs = 2/3/4);
S₆·S_T < 0 at every seed, the sign law's product form).
And at a simple discriminant zero this note's own lemma gives exactly two branches meeting in a
square-root point, so algebraic multiplicity 2 (H1) *follows* on A₁ rather than being assumed.
Hence:

> **At N = 5, 7 and 9, both R-parities: every forced count-drop is a √-type defective EP2 with
> s₆ ≠ 0. The per-N O2b nonvanishing holds there, eigensolver-free and exactly, at all three N**
> (the two premises the first N = 7 landing carried are both discharged 2026-07-16: the
> base-polynomial grade by the a-priori ℓ1 bound, the layer identification by
> reconstruct-and-verify; the two premise-discharge subsections below. N = 9 runs the same
> proof shape end to end, landed 2026-07-17).

Grades, in one line: N = 5 unconditional (exact ℤ end to end, base polynomials by exact integer
Faddeev-LeVerrier, inventory included); N = 7 unconditional over ℚ since 2026-07-16, both former
premises discharged: the CRT verification grade of the base polynomials (proof grade, first
subsection below) and the mod-p layer identification (proved, second subsection below); N = 9
unconditional over ℚ since 2026-07-17, by the same two mechanisms (never carrying the premises:
it landed with them already built in).

Certification bookkeeping, so nothing is over-read. Proved over ℚ outright: the N = 5 facts, exact ℤ
arithmetic end to end (independently reproduced by a second exact engine during the review round).
The N = 7 statements carried, until 2026-07-16, one **identification premise** that had to be named:
everything at N = 7 is computed mod p, and the identification "the mod-p multiplicity-1 layer = A₁
mod p" can fail at a bad prime (a mod-p merge of two A₁ roots moves degree from the simple to the
doubled layer while preserving the total degree, the valuation, and max multiplicity ≤ 2, so it is
invisible to the fail-closed DISCMULT logic, whose one-prime proof covers only the multiplicity
*bound*); multi-prime agreement was strong evidence with a named finite failure set, not a ℚ-proof.
**That premise is DISCHARGED** (second premise-discharge subsection below): A₁ and A₂ are now
reconstructed as exact primitive ℤ polynomials and disc = C·w^v·A₁·A₂² is proved over ℤ, so no
degree/shape-preserving prime can misbook a layer, and the committed verifier asserts the Yun
layers against the exact A₁/A₂ at every certificate prime, at every certified N. On that footing
the N = 7 and 9 gcds are rational: F_res is monic and p ∤ lc_λ(S₆), so the mod-p reduction of the
resultant is exact and Gauss's lemma lifts coprimality from one good prime; and the N = 7 and 9
irreducibility certificates (mod-p degree-partition subset-sum, fed by reductions of the proved
exact A₁, never by unverified Yun layers) prove irreducibility of A₁ itself. The
S₆-core irreducibility at N = 7 and 9 never needed the premise (the core comes from the exact CRT
polynomial). The finer split Res = c·w^m·A₂²·B with gcd(B, disc) = 1 is exact at N = 5 but
multi-prime evidence only at N = 7 and 9 (finitely many primes never lift divisibility), and nothing
above uses it. One accounting for the base polynomials, applied everywhere: at N = 5, χ, S₆, S_T are
computed by exact integer Faddeev-LeVerrier with the node count justified by the degree bound (no
CRT, no stability stop), so the N = 5 chain is exact ℤ end to end; at N = 7 they were, at the first
landing, CRT-with-stability, verification-grade rather than proof-grade. That grade was upgraded
2026-07-16 (subsection below): the CRT prime pool is now consumed until its product exceeds twice a
rigorous a-priori coefficient bound, which makes the symmetric-range reconstruction provably exact,
so every N = 7 statement, including fact 1's core irreducibility and the premise-free-of-layer
fact 4, now rests on proof-grade base polynomials; the layer-identification premise above was the
one N = 7 premise left, and it fell the same evening (second subsection below). N = 9 entered
(2026-07-17) with both mechanisms already in the engine, so it never carried either premise.

### Premise discharge I: the N = 7 (and 9) base polynomials are proof grade (2026-07-16; N = 9 on 2026-07-17)

The bound is the permanent/row-sum bound, and it is deliberately elementary. χ = det(λI − D − uK)
is a signed sum over permutations of entry products; the ℓ1 norm of coefficients (in the two
variables λ, u) is submultiplicative on products and subadditive on sums, so
‖χ‖₁ ≤ perm(𝒩) ≤ ∏_a rowsum_a, with 𝒩 the matrix of entry-wise ℓ1 norms
(rowsum_a = 1 + |D_a| + Σ_b |K_ab|) and the permanent of a nonnegative matrix bounded by its
row-sum product. Each adjugate diagonal entry is an (a,a) minor, so ‖adj_aa‖₁ ≤ ∏_{a′≠a} rowsum_{a′}
(deleting a column only lowers row sums), giving Σ-of-minor bounds for S₆, S₂, S_T. The u → w = u²
collection merges no terms (everything is even in u), so the same integers bound every (λ, w)
coefficient. `l1_bounds()` in [the gcd certificate](../simulations/o2b_gcd_certificate.py) computes
them exactly, the CRT loop refuses to close before ∏p > 2B, and the numbers come out comfortable:
at N = 7 the bounds are 2^251..2^254 (R-even, dim 75) and 2^239..2^242 (R-odd, dim 72), met at 11
25-bit primes with ∏p = 2^274; at N = 9 they are 2^563..2^567 (R-even, dim 164) and 2^547..2^551
(R-odd, dim 160), met at 23 primes with ∏p = 2^574; at N = 5 the exact-ℤ path validates the bound
from below (max|coeff| ≤ B asserted for all four polynomials, both sectors). With ∏p > 2B the symmetric-range
CRT reconstruction is exact by construction, no stability assumption left; the stability stop and
the fresh-prime extra-node checks stay in as independent cross-checks. What this does NOT touch:
the layer-identification premise (a factorization statement about disc_Λ(F_res), downstream of but
not reducible to exact base polynomials); its natural discharge is a reconstruct-and-verify pass,
built and landed the same evening:

### Premise discharge II: the layer identification is proved (2026-07-16, evening)

Reconstruct, then verify; the reconstruction may be heuristic because the verification is a
complete gate. **Stage L (lift):** at each good prime the Yun squarefree layers of disc_Λ(F_res)
mod p are monic images of the true layers; CRT across primes plus per-coefficient rational
reconstruction (half-extended Euclid, |num|·|den| < M/2) recovers the monic ℚ-coefficients, and
clearing denominators gives primitive ℤ candidates A₁, A₂ with lc > 0 (stability under 2 extra
primes as a cross-check, not the proof). **Stage V (the proof):** (i) B_D, a rigorous ℓ1 bound on
disc's coefficients, by the same permanent/row-sum argument one level up: disc = det of the
(2m−1)×(2m−1) Sylvester matrix of (F_res, ∂_λF_res), so ‖disc‖₁ ≤ (Σ_dl ‖F_dl‖₁)^{m−1} ·
(Σ_dl dl·‖F_dl‖₁)^m from the exact F_res coefficients; (ii) C from the leading-coefficient ratio
lc(disc mod p)/(lc(A₁)·lc(A₂)²) by CRT at degree-stable primes with ∏p > 2·B_D, node-independent,
never fitted to the identity; (iii) the identity **disc = C·w^v·A₁·A₂²** coefficient-wise mod
every prime of a pool with ∏p > B_D + B_R, B_R = |C|·‖A₁‖₁·‖A₂‖₁² a rigorous bound on the right
side, so the difference, integer-bounded by B_D + B_R and divisible by the pool product, is zero
over ℤ; (iv) SHAPE: A₁, A₂ squarefree, coprime, and nonzero at w = 0, each lifted from one
degree-preserving prime by Gauss (asserted at 3). Uniqueness of the squarefree decomposition then
forces: **the multiplicity-1 layer of disc off w = 0 IS A₁ and the multiplicity-2 layer IS A₂
over ℚ**: the identification is a theorem, and the verifier now asserts the mod-p Yun layers
against the exact A₁/A₂ at every certificate prime, at every certified N; the A₁-irreducibility
certificate
is fed by reductions of the exact A₁ (a fresh unverified Yun layer could otherwise smuggle the
premise back in, a review catch). Numbers: N = 7 R-even B_D = 2^14200, pool 569 primes
(∏ = 2^14224), C 442 digits; R-odd B_D = 2^13723, pool 549 primes, C 408 digits; ~11-12 min per
sector. At N = 9 the same stages run at scale (2026-07-17): R-even B_D = 2^69908, Stage L stable
at 504 primes, identity pool 2797 primes with ∏ = 2^69922 > B_D + B_R, C 931 digits, sign +;
R-odd B_D = 2^68845, Stage L 480 primes, identity pool 2754 primes (∏ = 2^68847), C 956 digits,
sign −; the disc-mod-p sweeps are fanned out over a multiprocessing
pool (~1-1.5 h per sector at 22 workers; the per-prime math is byte-identical to the serial path).
At N = 9 this identity is doing double duty: **max multiplicity 2 off w = 0 at N = 9 is itself a
corollary** (the proved shape C·w^v·A₁·A₂² with A₁, A₂ squarefree and coprime has no
multiplicity-3 layer), which the DISCMULT engine could never reach at N = 9; the localization
half of the argument is thereby served without DISCMULT.
Validated from below at N = 5: the reconstructed (C, A₁, A₂) equal the exact-ℤ inventory
exactly, and max|coeff(disc)| = 2^430/2^372 ≤ B_D = 2^1548/2^1393. What this deliberately does
NOT upgrade: the A₂²-split of Res_Λ(F_res, S₆) stays multi-prime evidence (nothing uses it), and
A₂'s own irreducibility at N = 7 and 9 stays evidence-grade (nothing uses it either).

What this hands the all-N question: **the seed-bearing layer does not fragment.** A₁ is irreducible
over ℚ in both sectors at all three N (six data points: hypothesis (i) below). If A₁ is
irreducible at every odd N, then
gcd(Res_Λ(F_res, S₆), A₁) is 1 or A₁, so the *resultant's* nonvanishing at one point of the simple
layer decides all of it: the **simple-layer half** of O2b reduces to (i) A₁-irreducibility for all
odd N and (ii) Res_Λ(F_res, S₆) ≠ 0 at one accessible w° with A₁(w°) = 0, i.e. S₆ ≠ 0 at **every**
F_res branch over that one fiber, not merely at the defective one (S₆ could vanish at a spectator
branch while the seed's s₆ ≠ 0 still holds, so a one-branch evaluation would neither force gcd = 1
nor would gcd = A₁ falsify O2b). The *localization* half, "every count-drop sits on the simple
layer", stays per-N: it consumed max multiplicity 2 (DISCMULT at N = 5 and 7; at N = 9 the proved
layer identity delivers it directly, the paragraph above), the psc₁ coincident-pair exclusion,
the cross-sector gcd, and per-sector reality of F_res,
none of which the two hypotheses subsume. That two-hypothesis shape for the simple-layer half
is new; the previous all-N form (the "Three attacks" caution) had no way to separate the diabolics
inside {S₆ = 0} from the seeds, and here the separation is the layer structure itself: the diabolics
sit on A₂ (at N = 5, R-even has no real positive A₂ root at all, R-odd exactly one, the known
diabolic at w = 5.100831, where S₆ does vanish, as adj ≡ 0 forces), the seeds on A₁, and the four
real positive A₁ roots at N = 5 are exactly the four seeds: the inventory is closed (asserted
exactly over ℤ by the committed verifier).

One numerical trap, recorded for reuse: S₆'s integer coefficients cancel heavily at seed loci
(8 to 18 digits, seed- and measure-dependent), so a locus refined only to 1e−8 can flip the
*apparent* sign of S₆; the seed must be pinned to
~1e−20 by explicit 2D Newton on (F_res, ∂_λF_res) with the exact Jacobian at ≥ 60 digits before
reading any polynomial value there. mpmath's `findroot` absolute tolerance is useless at these
coefficient scales.

Committed verifier: `simulations/o2b_gcd_certificate.py` (self-asserting at N = 5, ~10 s; `7` adds
N = 7, ~30 min since the layer discharge runs inside it; `9` adds N = 9, ~3.4 h measured with the
disc sweeps on a 22-worker pool; it recomputes the polynomials, the
factorizations, the exact layers A₁/A₂ with the full reconstruct-and-verify proof, all three gcd
families including the psc₁ coincident-pair certificate, the strand positivity symbolically, the
exact N = 5 disc/inventory over ℤ, and the seed evaluations from scratch; its docstring carries the
same grade ladder as the bookkeeping paragraph above).

## The cell law, tightened: decomposition, sharp constants, and where the band tightness comes from (2026-07-10)

The "Three attacks" section left the sign half on two band-wide candidates, the separating law and
the stronger cell form b′ ≥ b ⟹ a′ ≥ a. A dedicated session sharpened both. Notation as there
(v = x + iy in the gauge conj(v) = Tv, x on class E, y on class O; a = ‖(Kx)₂‖², b = ‖(Kx)₆‖²,
a′ = ‖(Ky)₂‖², b′ = ‖(Ky)₆‖²); write Δ₂ := ‖x₂‖² − ‖y₂‖² and Δ₆ := ‖x₆‖² − ‖y₆‖², so that
κ₋₂‖v‖² = Δ₂, κ₋₆‖v‖² = Δ₆ (rung-wise proportionality, exact: (λ+2)x₂ = −q(Ky)₂ and its three
siblings), and the cell form says the quadrant (Δ₆ ≥ 0, Δ₂ < 0) is empty.

- **The decomposition (proved; here are the four lines).** Cell form ⟺ separating law
  (σ ≥ 0 ⟹ κ₋₂ ≥ 0) **and** its rung-6 dual (σ < 0 ⟹ κ₋₆ < 0 strictly). Proof, with
  σ‖v‖² = Δ₂ + Δ₆: (⟸) if Δ₆ ≥ 0 then σ < 0 is excluded by the dual, so σ ≥ 0, so Δ₂ ≥ 0 by the
  separating law. (⟹) If σ ≥ 0 and Δ₂ < 0 then Δ₆ = σ‖v‖² − Δ₂ > 0, which the cell form forbids
  alongside Δ₂ < 0; if σ < 0 and Δ₆ ≥ 0 then the cell form gives Δ₂ ≥ 0, so σ‖v‖² ≥ 0,
  contradiction. ∎ The dual is exactly the "strictly stronger"
  content; it was tested standalone: 0 violations in 3575/3615/1036 negative-type samples at
  N = 5/7/9. The cell form itself now stands at **N = 5, 7, 9, 11** (0 violations in 2906, 2218,
  1172, 894 + an independent 208-sample N = 11 audit; N = 11 is the first *resonant* datum, in the
  resonant-N section's sense).
  Practical note for re-runs: the λ = −2 band edge carries pure ker-K modes whose four cells are all
  ~1e−30; without a 1e−16 floor on the cell norms they read as phantom violations.
- **Sharp constants with a rational equality family.** On positive-type vectors the separating law
  tightens to **κ₋₂ ≥ 2σ/(N+1)** (measured minima of κ₋₂/σ: 0.33648, 0.25104, 0.20101 at
  N = 5, 7, 9 against 1/3, 1/4, 1/5, approached from above as q → ∞); equivalently, on
  positive-type branches with σ > 0, qλ̇ ≤ λ − λ∞ᴱ (at σ = 0 the κ-form is vacuous while the slope
  form still constrains; quote the κ-form as the law). The equality family is exact and rational: at
  q → ∞ the *in-band, mixed-rung* all-class-E ker(K) modes sit at λ∞ᴱ = −(6N−2)/(N+1) with
  (κ₋₂, κ₋₆, σ) = (2/(N+1), (N−1)/(N+1), 1) exactly, the in-band all-class-O modes at −6(N−1)/(N+1)
  with (−3/(N+1), −(N−2)/(N+1), −1) (verified rationally at N = 5, 7, 9; the *pure-rung* ker(K)
  modes sit at the band edges −2 and −6 instead, the same family the first bullet's practical note
  floors out). The dual's candidate sharpening, **scoped to negative-type vectors (σ < 0):
  κ₋₆ ≤ σ/3**, a strengthened negativity there (unscoped it is false: on positive type the equality
  family itself has κ₋₆ = (N−1)/(N+1) > σ/3), has thin positive margins (~1.3e−3 on the densest
  sweep, growing in N); its extremizers look generic (no closed form surfaced; the forced/free law
  again).
- **A fixed-form reformulation, banked for a global attack.** z := x + y satisfies the *real* affine
  pencil (D + qKT − λ)z = 0, and a′ − a = zᵀ(KTP₂K)z, b′ − b = zᵀ(KTP₆K)z with KTP₂K, KTP₆K fixed
  integer symmetric matrices: the entire law is two fixed quadratic forms on the null vectors of a
  one-knob real family. (The pointwise S-procedure on the corresponding forms fails, measured:
  relaxing the kernel constraint to its quadratic form is too weak, so any operator proof must
  consume the full kernel structure.)
- **The moment algebra closes; one scalar is free.** With p₂ := x₂ᵀKy₂, p₆ := x₆ᵀKy₆,
  c₂₆ := x₂ᵀKy₆, c₆₂ := x₆ᵀKy₂ (all four bilinears of the rung-split eigen-equations), the exact
  system gives: p₂ and p₆ both lie in the interval [−c₂₆, −c₆₂] (endpoint order −c₂₆ ≤ −c₆₂,
  forced by λ + 2 < 0), the normalized positions being the
  shares e₂ := ‖x₂‖²/(‖x₂‖² + ‖y₂‖²) and o₆ := ‖y₆‖²/(‖x₆‖² + ‖y₆‖²), and the cell form forbids
  exactly p₂ strictly below the midpoint m := −(c₂₆ + c₆₂)/2 while p₆ is at or below it (the p₆ = m
  boundary, κ₋₆ = 0, is where the rung-6 dual's strictness lives). The dead end xᵀKy ≥ 0 is
  *exactly* the odds-ratio strengthening ‖x₂‖‖y₆‖ ≥ ‖y₂‖‖x₆‖ (e₂ + o₆ ≥ 1) of the true law, which
  locates precisely why it failed. The 4×4 moment system has rank 3: one scalar is free, so no
  moment identity can decide the law; a proof must consume K's internal structure. A new reduction
  with clean margins: the law follows from [m > 0 ⟹ Δ₂ < 0 ∧ Δ₆ < 0] (an independent measured law,
  0 violations, margins ≥ 0.02 at N = 5..11) together with the m ≤ 0 core
  [σ ≥ 0 ⟹ p₂ ≥ m] ∧ [σ < 0 ⟹ p₆ > m], whose margins ~0.002 are strangely N-stable across
  N = 5..11 (unexplained; the one place the law is nearly tight away from seeds).
- **Edge-birth theorem (proved, second-order degenerate perturbation theory; the derivation, so the
  label carries something checkable).** A level leaving the −2 rung has first-order shifts given by
  the eigenvalues of the compression of qC = iqK to the rung, purely imaginary unless zero (the
  pencil-section fact), so only rung-kernel modes stay real past q = 0. For ξ ∈ ker(K₂₂) ⊂ E with
  K₆₂ξ ≠ 0 (the pure ker(K) directions have K₆₂ξ = 0 and stay at the edge λ = −2, never entering
  the band): x = ξ + O(q²), y₆ = (q/4)K₆₂ξ + O(q³), y₂ = O(q³), λ = −2 − (q²/4)‖K₆₂ξ‖², giving
  Δ₂ = ‖ξ‖² + O(q²) > 0 and Δ₆ = −(q²/16)‖K₆₂ξ‖² < 0 strictly. For ξ ∈ ker(K₆₆) ∩ O with
  K₂₆ξ ≠ 0: y = ξ + O(q²), x₂ = (q/4)K₂₆ξ + O(q³), λ = −6 + (q²/4)‖K₂₆ξ‖², giving Δ₂ > 0,
  Δ₆ = −‖ξ‖² + O(q²) < 0. So every real branch
  entering the band from an edge at q → 0⁺ comes from a rung-kernel mode (proved), and the two
  baseline families that actually enter, ξ ∈ ker(K₂₂) with K₆₂ξ ≠ 0 (class E) and ξ ∈ ker(K₆₆) ∩ O
  with K₂₆ξ ≠ 0, enter with (Δ₂ > 0, Δ₆ < 0)
  strictly (proved): **at every N where ker(K₆₆) is pure class O and ker(K₂₂) pure class E, all
  real edge births start in the allowed corner, as a theorem.** Both purity conditions are inputs:
  a hypothetical ξ ∈ ker(K₂₂) ∩ O would enter with the mirrored signs, i.e. forbidden, so the
  ker(K₂₂)-in-E localization (the endpoints-in-E add-on, asserted at the swept N, now including
  10E+0O at N = 11 and 16E+0O at N = 17) carries load here too. A birth from ξ ∈ ker(K₆₆) ∩ E with
  K₂₆ξ ≠ 0 would enter in the
  forbidden corner; that no such birth materializes at resonant N is the *measured* twinning input
  (N = 11, 17; open in general), so there the allowed-corner statement is theorem plus measurement,
  not theorem. The resonant-N section of [F89_SEED_EXISTENCE_REDUCTION.md](F89_SEED_EXISTENCE_REDUCTION.md) holds that piece.
- **Where the band tightness lives.** The corner (κ₋₂, κ₋₆) → (0, 0) is never approached: the
  band-wide floor of max(|κ₋₂|, |κ₋₆|) equals the minimal seed value κ₋₂\* = |s₆|/‖r‖² (0.0629 at
  N = 5, 0.0253 at N = 7, both attained at *entry* seeds), and far from seeds the floor is an order
  of magnitude higher. So no separate band mechanism binds: the band laws' global tightness IS the
  seed-local O2b quantity, one more way the sign half and the nonvanishing half are the same flesh.
- **New dead ends, measured** (join the do-not-retry list): the swap (p₂ < 0 ∧ p₆ < 0 ⟹ σ < 0)
  is false (359 counterexamples at N = 5); pointwise kernel-refuge (kernel share of x₂ forcing the
  class-E share ≥ ½) is false (shares 0.003..0.86 on both sides); any uniform gap
  max(e₂, o₆) ≥ ½ + c dies (sample minimum 0.539 at N = 7, N-decreasing); rung-restricted Loewner
  definiteness of KTP_rK is absent (mixed inertia).

These are in-session measurements, two independent recomputations where marked; the committed
verifier for this block's predecessors remains `simulations/o2b_three_attacks_audit.py`, and the
resonant-N facts in [F89_SEED_EXISTENCE_REDUCTION.md](F89_SEED_EXISTENCE_REDUCTION.md) are flagged for a committed witness of their own if they become load-bearing.



## Reproduce

```bash
python simulations/seed_existence_nullity_check.py     # the seed object + the count (parent note)
python simulations/o2b_krein_sign_law.py               # the class-imbalance sign law, N = 5, 7 (add 9)
python simulations/o2b_three_attacks_audit.py          # the "Three attacks" section (N=5; also 7/b7/b9)
python simulations/o2b_gcd_certificate.py              # the gcd certificate section (N=5; add 7 and/or 9)
```

## Status

The count r(0⁺) − r(∞) = N − 1 is a theorem (the parent
[F89_SEED_EXISTENCE_REDUCTION.md](F89_SEED_EXISTENCE_REDUCTION.md), Pieces 1-3). What this note owes the
seed-existence *conclusion* is that each count-dropping seed is genuinely √-type:

**The codim-2 β-exotic (OPEN for all N; CLOSED EXACTLY at N = 5, 7 and 9, both parities,
   2026-07-09/10, N = 9 on 2026-07-17):** a count-dropping transition is a √-type defective EP2 (exponent ½) unless it is a
   non-generic order-≥3 point: the formally semisimple β-exotic (normal form β(s) = [[0, s], [s², 0]],
   eigenvalues ±s^{3/2}), its defective exponent-3/2 twin (geometric multiplicity 1, s₆ = 0), or the
   higher-order relatives (Jordan size ≥ 4, higher odd exponents), all caught by the same bound.
   - *At N = 5 and N = 7 it is settled*, unconditionally, by the certified disc-multiplicity reading
     (section "The β-exotic is excluded at N = 5 and N = 7" above): max root multiplicity 2 off q = 0,
     and a β-exotic would need 3. The subsidiary premise H1 (algebraic multiplicity exactly 2, no
     Jordan block of size ≥ 3) **follows there, with one checked premise**: the argument needs F_res
     real, whose per-sector reality rests on the chiral pairing being what fixes the AT slopes,
     checked at N = 5 and N = 7 rather than derived in general, and at N = 9 by the gcd
     certificate's runtime asserts (the certificate section's reality
     paragraph is the canonical wording). Cite H1 with that qualifier. N = 9 is out of reach by the
     DISCMULT route; there the max-multiplicity-2 bound comes from the proved layer identity
     instead (the second premise-discharge subsection).
     **2026-07-10, strictly more**: the gcd certificate (its own section above) proves the
     nonvanishing itself, s₆ ≠ 0 at every forced count-drop
     (gcd(Res_Λ(F_res, S₆), A₁) = 1 plus the cross-sector gcd): unconditionally at N = 5 (exact ℤ),
     unconditionally at N = 7 since 2026-07-16 (both premises of that section's bookkeeping
     paragraph discharged: the base-polynomial CRT grade by the a-priori ℓ1 bound, the mod-p layer
     identification by reconstruct-and-verify; the two premise-discharge subsections), and
     unconditionally at N = 9 since 2026-07-17 (same engine, premises built in from the start).
     So at those N the O2b scalar is not
     just β-excluded but positively settled: every forced count-drop is a √-type defective EP2.
   - *For all N it is open*, and reduced to the scalar statement **s₆ ≠ 0 at every forced seed for all
     odd N** (section "The β-exotic, sharpened" above; read via the polynomial S₆ = tr(P₋₆·adj(λI−L))
     where the eigenspace could be two-dimensional, since s₆ itself is defined only at geometric
     multiplicity 1), which excludes the β-exotic only under the
     separate premise H1 (algebraic multiplicity exactly 2, only numerical for N ≥ 11). It is a genuine
     sharpening of the failed discriminant-Galois route. The moment-relation proof route (transpose +
     Hermitian relations from the rung-split eigen-equation) is proved insufficient by an explicit
     in-class witness. The remaining need, a definiteness ingredient from the K₆₆ = −3H₃ spectrum or
     the N − 1 path structure, now has a **named candidate**: the class-imbalance sign law
     sign(κ_rung) = sign(|E_rung| − |O_rung|), equivalently κ₋₂ > 0, which holds at every defective
     coalescence at N = 5, 7, 9 (both parities; verifier `simulations/o2b_krein_sign_law.py`) but is
     not proved. It is the gauge-invariant, signed reading of s₆ ≠ 0 at the *defective* seeds (the
     gauge relation κ₋₆ = c·s₆/‖v‖² holds only at geometric multiplicity 1), and a from-below
     mechanism candidate, not a standalone β-exotic exclusion: that still needs the defectiveness
     (geometric multiplicity 1) itself, which the certificate supplies at N = 5, 7, 9 and the σ₂
     sweep corroborates numerically at the same N. The 2026-07-09 section "Three attacks on the sign law"
     re-based both halves: the sign law's *sign half* reduces (modulo H1) to the band-wide *separating
     law* (giving κ₋₂ ≥ 0 at every defective seed, strict exactly where s₆ ≠ 0), and the nonvanishing
     acquires the integer-polynomial form κ₋₂ = −S₆/S_T with S₆ ≠ 0 ⟹ geometric
     multiplicity 1 ∧ s₆ ≠ 0; the two named targets live there (verifier
     `simulations/o2b_three_attacks_audit.py`). **2026-07-10 re-base of the all-N shape**: the gcd
     certificate section replaces the old caution ("{S₆ = 0} contains the diabolics, no separation")
     with a measured separation, the diabolics on A₂ and the seeds on A₁, and reduces the
     *simple-layer half* of all-N O2b to
     two uniform hypotheses: **A₁ irreducible over ℚ for every odd N** (true at N = 5, 7, 9, both
     sectors: six data points) **plus Res_Λ(F_res, S₆) ≠ 0 at one accessible simple-layer fiber**
     (S₆ ≠ 0 at every
     branch over it, not merely at the defective one); the localization half (every count-drop on
     A₁) keeps its per-N inputs (max-mult 2 by DISCMULT at N = 5, 7 and by the layer identity at
     N = 9, the psc₁ leg, the cross-sector gcd, F_res
     reality). The sign half's
     band-wide targets are tightened in the cell-law section (the separating law with sharp constant
     κ₋₂ ≥ 2σ/(N+1), its rung-6 dual, the edge-birth theorem, and the twinning protection at
     resonant N, which is also the cheapest falsification target).
