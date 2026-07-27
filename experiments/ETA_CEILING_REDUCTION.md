# The ceiling, reduced: an exact certificate on the multiplet count, and the law that closed it

*2026-07-27. [XY_FROZEN_BAND](XY_FROZEN_BAND.md) left the ceiling settled one N at a time, by an exact rank read on the block itself, and said what was missing was the argument that does all N at once. This note is the measuring that found the law that argument turned out to be. It moves the per-N certificate off the block and onto an object smaller by two orders of magnitude, so a rung that took hours on the block takes seconds and rungs that were out of reach become routine; it states the scope of that carefully below, because the two certificates do not certify the same thing, the block-level one settling a whole band at N = 9 and this one settling a rung at a time. It turns the certificate into a NON-existence question, which is the direction a rank over a finite field can actually settle and the earlier direction could not. And it states the closed form the runs kept reproducing, together with the vectors that attain it and the reason exactly one N escapes. The inequality half of that closed form is now proved for every N, in [PROOF_FROZEN_BAND_SO4](../docs/proofs/PROOF_FROZEN_BAND_SO4.md) Section 7, so the per-N runs below are a check rather than the argument. The structural half, the decomposition that makes all of this a question about where multiplets start, is proved in the same document; this note stays the measuring.*

## What is being measured, and on what

By [F141](../docs/ANALYTICAL_FORMULAS.md) the disagreement count 𝒦 commutes with the η-ladder, so it acts only on the multiplicity spaces of the sl(2), and the spectrum on a diagonal rung p is the union of the lowest-weight spectra at rungs ℓ ≤ min(p, N−p).

**The symbols, in one place.** A cell |A⟩⟨B| of the open XY chain has ket and bra site sets A and B; **𝒦** is the disagreement count |A \ B|, and **D̂** the double occupancy |A ∩ B|, so on a diagonal rung 𝒦 = ℓ − D̂ and "min 𝒦" and "max D̂" are one statement. **Φ** = Σ_l d†_l ρ d_l is the η ladder, (p,q) ↦ (p+1,q+1), and **Ψ := Φ†** removes such a pair; **LW_ℓ** is the lowest-weight space of that ladder at rung ℓ, that is ker Ψ, never of the spin ladder beside it. **M := N + 1** is the length of the discrete sine transform the open chain diagonalises under, with modes v_a(l) = √(2/M)·sin(πal/M) and energies 2J·cos(πa/M); it is neither a popcount nor a prime. **ad_h** is the map X ↦ [h, X] for h the single-excitation hopping matrix, and **m_E** the number of ℓ-subsets of modes at Slater energy E. All of these are defined and used in [PROOF_FROZEN_BAND_SO4](../docs/proofs/PROOF_FROZEN_BAND_SO4.md); they are repeated here because this note is read on its own. One word is worth pinning: **seed** here is the bottom rung of a ladder and the frozen modes that start there, not MirrorWorld's `Seed`, which is the within-block self-dual seed [F89](../docs/ANALYTICAL_FORMULAS.md) holds as a count. Neither statement depends on the other.

That alone is not yet the ceiling, and the restriction that turns it into one is **not** a convenience. At finite coupling "frozen" is not "𝒦 = 1": the frozen condition mixes 𝒦 with the hopping, and 𝒦 does not commute with the hopping. The two agree only on

  **V₀ := ker(ad_h)**,

spanned by the pairs of ℓ-subsets of modes with equal Slater energy, which is exactly where the large-J reduction of [PROOF_FROZEN_BAND_SO4](../docs/proofs/PROOF_FROZEN_BAND_SO4.md) Section 5 puts the question. On the FULL lowest-weight spaces the same sentence would be vacuous: 𝒦 is diagonal with integer entries there and preserves ker Ψ, so the eigenvalue 1 sits in it with large multiplicity at essentially every N, the ones where the ceiling holds included. So the sentence to test is

**no LW_ℓ ∩ V₀ with 2 ≤ ℓ ≤ min(p, N−p) carries the eigenvalue 𝒦 = 1,**

and it is **sufficient, not necessary**: it implies the ceiling at rung p, and the converse fails at N = 5, where the reduced space does carry the eigenvalue and the ceiling holds anyway. Everything below tests the sufficient condition, and the one N where it fails is handled separately.

The reduction is also what makes it cheap. The block of dimension C(N,ℓ)² gives way to V₀ of dimension Σ_E m_E², and the gap between the two is why any of this runs:

| N | ℓ | block C(N,ℓ)² | dim V₀ |
|---|---|---|---|
| 9 | 4 | 15 876 | 252 |
| 12 | 6 | 853 776 | 3 584 |
| 16 | 8 | 165 636 900 | 97 444 |

Two things had to be got right before any of it could be trusted. The equal-energy classes are computed **exactly**, not to a tolerance: 2cos(kπ/M) = ζ^k + ζ^(−k) for a primitive 2M-th root of unity, so a Slater energy is an integer polynomial in ζ and two energies agree exactly when their difference is divisible by the cyclotomic polynomial Φ_2M. Sums of cosines get arbitrarily close without being equal, and one wrongly merged class would corrupt dim V₀ and everything after it. And the whole reduced operator is **integer**: expanding the four sine modes and summing the geometric series gives

  2M·Σ_l v_a v_b v_c v_d = Σ over the eight sign patterns of (∏σ)·[2M divides σ_a a + σ_b b + σ_c c + σ_d d],

zero when a + b + c + d is odd. The eight signs sum to zero, so no constant survives, and Z := 2M·D̂ is an integer matrix. Everything below is arithmetic on Z.

## The certificate, and the direction that makes it work

The ceiling is a **non**-existence: 𝒦 = 1 must not occur on a lowest-weight space. Writing it out, a violating vector would satisfy (Z − 2M(ℓ−1))v = 0 and Ψv = 0 at once. So

  **LW_ℓ ∩ V₀ carries no 𝒦 = 1 ⟺ the stacked integer matrix [ Z − 2M(ℓ−1)·I ; Ψ ] has full column rank,**

and that in turn IMPLIES the ceiling at rung ℓ, one way and not the other. The biconditional is between the rank and the absence of the eigenvalue; between the absence and the ceiling there is only the implication, as N = 5 shows by holding the ceiling while failing the rank.

A rank over GF(q) can only **drop** relative to the rank over ℚ. So full rank mod q proves full rank over ℚ, which proves the non-existence. That is the direction the ceiling needs, and it is cheap. The opposite question, whether some value **is** an eigenvalue, cannot be settled mod q for the same reason and has to be done over ℚ; both directions are used below, each where it works.

What it reads, over two large primes (no square root of −1 is needed, the matrices being integer rather than Gaussian-integer). **Four of these six rows need `--deep`**: the plain run caps the reduced dimension at 700 and stops at N = 13, so it reaches only the N = 5 and N = 9 rows below, and the seconds are one machine's, recorded for the ratio between rows rather than as a reproducible number.

| N | ℓ | dim V₀ | rank | full | seconds |
|---|---|---|---|---|---|
| 5 | 2 | 12 | 11 | **no** | 0.0 |
| 9 | 4 | 252 | 252 | yes | 0.1 |
| 11 | 5 | 2 100 | 2 100 | yes | 35.1 |
| 12 | 4 | 1 425 | 1 425 | yes | 12.5 |
| 13 | 4 | 1 885 | 1 885 | yes | 32.1 |
| 16 | 3 | 1 232 | 1 232 | yes | 6.4 |

Every rung tested is certified except one, and the exception is N = 5, which the reduction had already flagged from the float side.

**Scope, stated plainly, because the comparison is easy to overstate.** Three separate restrictions, and none of them is decoration.

- **The certificate is per RUNG, not per N.** Certifying rung ℓ rules out 𝒦 = 1 on LW_ℓ ∩ V₀. The ceiling at a band block p needs that for every ℓ ≤ min(p, N−p). So a band block is certified only when all of its rungs are.
- **It certifies the reduced first-order operator**, which by [PROOF_FROZEN_BAND_SO4](../docs/proofs/PROOF_FROZEN_BAND_SO4.md) Proposition 5.1 bounds the true multiplicity from above **at large J**. With the proved floor holding at every coupling, that makes the multiplicity exactly ⌊N/2⌋ there, on the rungs that are certified. Carrying it to generic J is the further upper-semicontinuity step named after that proposition, which is asserted and not argued, so "at generic coupling" is a stronger sentence than anything the certificate delivers by itself.
- **The reach is not what a single headline number suggests.** The block-level certificate of [XY_FROZEN_BAND](XY_FROZEN_BAND.md) stops at N = 9, takes about three and a half hours there, and settles the WHOLE band. The committed gate here settles **every** rung only to N = 11; at N = 12 and 13 it reaches ℓ ≤ 4, and at N = 14, 15 and 16 only ℓ ≤ 3, because the size cap binds. The middle rungs above those, up to dim V₀ = 97 444 at N = 16, ℓ = 8, are not certified by it at all. So "further than N = 9" is true rung by rung and false band by band, and the honest summary is that the cost per rung collapsed, not that the certified region simply grew.

## The law that replaced it

  **min 𝒦 on LW_ℓ ∩ V₀ = ℓ(N − ℓ)/(N + 1),  for ℓ ≥ 2.**

This note first stated the law as a conjecture with the evidence below it. The half of it the certificate needs, the inequality **≥**, is now proved for every N and every ℓ, in [PROOF_FROZEN_BAND_SO4](../docs/proofs/PROOF_FROZEN_BAND_SO4.md) Section 7 ([F144](../docs/ANALYTICAL_FORMULAS.md)): the double occupancy on a lowest-weight vector of V₀ obeys an exact identity, M·D̂ = ℓ(ℓ+1) minus three squares, whose leading term is the norm of the vector under the SPIN ladder. So the per-rung rank certificate below is superseded as an argument, though not as a check: what it computed one rung at a time now follows in one line for all rungs, and the table of runs stays as the record that the two agree. The **equality** half, that the bound is attained, is what remains measured: at ℓ = 2 the maximizers are the closed-form vectors below and gate V9 confirms them per N, and for ℓ ≥ 3 attainment is the evidence in this section and nothing more.

The intersection with V₀ is not decoration and dropping it makes the statement false: on the full lowest-weight space the minimum is 0, since 𝒦 is diagonal there with integer entries and preserves ker Ψ. Every law, bound and conjecture in this note is on **LW_ℓ ∩ V₀**, which is what the certificate reads and what the gate computes. Equivalently, in terms of the double occupancy, max D̂ on LW_ℓ ∩ V₀ = ℓ(ℓ+1)/(N+1). The evidence is exact rational agreement, never approximate: complete rung sweeps over every ℓ from 2 to ⌊N/2⌋ for **N = 4 through 14**, eleven consecutive N, and partial sweeps beyond, the largest single read being N = 19 at ℓ = 8 on a reduced space of dimension 921 852. Values run 4/5, 1, 8/7, 9/7, 5/4, 3/2, 4/3, 5/3, 16/9, 7/5, 9/5, 2, 16/11, 21/11, 24/11, 25/11, 3/2, 2, 7/3, 5/2, 20/13, 27/13, 32/13, 35/13, 36/13, 11/7, 15/7, 18/7, 20/7, 3, 8/5, and on. Zero misses. The committed gate reproduces the law on every rung to N = 11, on ℓ ≤ 4 at N = 12 and 13, and on ℓ ≤ 3 at N = 14, 15 and 16, the size cap binding above that. Everything past that, including the complete sweeps at N = 14 and the reads to N = 19, was measured in the session with the same sparse mode-basis construction and is recorded here as an observation rather than as gate output.

The form is worth reading. The uncorrelated expectation of 𝒦 for ℓ particles of each species spread over N sites would be ℓ(N − ℓ)/**N**. The law is that value with N replaced by the transform length N + 1, the same substitution that turns the seed gap into 1/(N+1) in [F143](../docs/ANALYTICAL_FORMULAS.md). Lowest-weight states are the maximally decorrelated ones, measured against the transform rather than against the sites.

**What the law buys.** The ceiling needs min 𝒦 > 1, that is ℓ(N − ℓ) > N + 1. On 2 ≤ ℓ ≤ ⌊N/2⌋ the product ℓ(N − ℓ) increases with ℓ, so the binding case is ℓ = 2 and the condition is 2(N − 2) > N + 1, that is **N ≥ 6**. One argument, every N, every rung. Only the inequality is used here, which is the proved half. What it then gives, and what it does not, is worth keeping apart. It gives the depth exactly ⌊N/2⌋ **in the large-J regime** at every N ≥ 6, because the multiplicity of 𝒦 = 1 on a rung becomes writable rather than measurable: the seed rung contributes ⌊N/2⌋ and no rung above it contributes anything, so [PROOF_FROZEN_BAND_SO4](../docs/proofs/PROOF_FROZEN_BAND_SO4.md) Proposition 5.1 caps the frozen space there, and the proved floor meets the cap. It does not give the same at generic coupling: that passes through the semicontinuity step of Section 5 of that document, which is not proved, or through the certificate below where the certificate reaches. The per-N computation left the argument; the bridge did not.

## Why exactly one N escapes, in one line

Equality ℓ(N − ℓ) = N + 1 rearranges to

  N = (ℓ² + 1)/(ℓ − 1) = ℓ + 1 + 2/(ℓ − 1),

which is an integer only when ℓ − 1 divides 2, so only for ℓ = 2 and ℓ = 3. The first gives N = 5; the second gives N = 5 as well, where ℓ = 3 exceeds ⌊N/2⌋ = 2 and no lowest-weight space exists. **N = 5 at ℓ = 2 is the unique solution.** So the single exception the certificate finds is not an accident of small N but the only place two integers can meet, and the resonance the reduction sees from the float side, the certificate sees as a rank drop, and the arithmetic sees as a divisibility, are the same fact three times.

At N = 5 the first order is loose by exactly one and the **third** order closes it, back to the proved floor, with a clean gap of 1/24 (the gap itself is a session measurement; the gate reads the resulting dimension, not the gap). The second order cannot help: by [PROOF_FROZEN_BAND_SO4](../docs/proofs/PROOF_FROZEN_BAND_SO4.md) Section 5 every even order vanishes identically on the surviving space, for a transpose-parity reason, so the next order after the first that can bite is the third.

At N = 4 the law gives min 𝒦 = 4/5 < 1, so the sufficient condition does not apply; the ceiling holds there anyway because the spectrum on that lowest-weight space is {4/5, 6/5, 6/5, 2} and simply misses 1. That case is settled by inspection, not by the criterion.

## The vectors that attain it

At ℓ = 2 the extremal vectors are known in closed form, and there is nothing approximate about them. Take two **chiral pairs** of modes, {p, s} and {q, r} with p + s = q + r = M and p < q < r < s, and form

  v = |pq⟩⟨pq| − |pr⟩⟨pr| − |qs⟩⟨qs| + |rs⟩⟨rs| − 2·(|ps⟩⟨qr| + |qr⟩⟨ps|).

Then **Z v = 12 v and Ψ v = 0**, exactly, in integer arithmetic, at every N from 5 to 16, which is D̂ v = (6/M)·v = (ℓ(ℓ+1)/M)·v on a lowest-weight vector.

One convention note, because a re-implementation will trip on it otherwise. The ±1 and ±2 above are written in the **Jordan-Wigner sign gauge** the V₀ basis carries, the one the gate's `sign_remove` applies to both the operator and the basis. In a raw unsigned matrix-unit convention only some choices of the two chiral pairs satisfy the relation, and the formula looks wrong. Everything gauge-invariant is unaffected: the eigenvalue, the multiplicity C(⌊N/2⌋, 2), and the surplus at the resonant N are convention-free. Six entries, values in {1, 2}, one vector per pair of chiral pairs. That is the lower half of the closed form, uniform in N and needing no computation per N; the upper half, that nothing exceeds it, is what stays open.

The count of these vectors is C(⌊N/2⌋, 2), and it matches the exact rational nullity of [Z − 2ℓ(ℓ+1)I ; Ψ] at every N tested except two: **N = 11 gives 14 against 10, and N = 14 gives 29 against 21.** That surplus is a claim that vectors EXIST, so it needs the rank over ℚ and not over a prime, where a rank can only drop; gate check V9 computes it fraction-free (N = 10 and N = 11 by default, N = 12 and N = 14 under `--deep`).

**Two thresholds, not one, and they are easy to conflate.** The certificate stacks Ψ under Z − 2M(ℓ−1)·I, the FROZEN value 𝒦 = 1. The surplus count stacks it under Z − 2ℓ(ℓ+1)·I, the MAXIMUM of D̂. Those are different matrices at every N but one: they coincide exactly when 2M(ℓ−1) = 2ℓ(ℓ+1), which at ℓ = 2 is M = 6, that is N = 5. The coincidence of the two thresholds is the resonance, seen a fourth way. Those are M = 12 and M = 15, and 3 | 12, 3 | 15, 15 | 15 are exactly the divisibilities carrying the [3|n] and [15|n] terms of the Conway-Jones resonance count that [F89_SEED_EXISTENCE_REDUCTION](F89_SEED_EXISTENCE_REDUCTION.md) already classifies. The extra vectors are that arithmetic showing up at a new place, not a new phenomenon.

## Two traps on the way, both of which passed a smaller check

Recorded because each survived a run and was caught only by a test pointing the other way.

**Divisibility is not equality.** The eight conditions in the overlap formula are "2M divides σ·(a,b,c,d)", not "σ·(a,b,c,d) = 0". Since a + b + c − d runs up to 3M − 4 it can hit 2M, and a − b − c − d runs down to 4 − 3M and can hit −2M. Read literally, with all eight conditions as "the signed sum vanishes", the formula is wrong on **226 of the 1296 quadruples at N = 6**. The form that actually shipped and had to be caught was subtler and still wrong on **80**: it had the all-plus condition right, since a + b + c + d can only ever reach 2M and equality is correct there, and the other seven wrong. That is why it survived a look. What caught it was the certificate returning "no full rank" on **every** rung, which would have broken a ceiling measured forty times over. A result that agreed with expectation would not have exposed it.

**An eigensolver's labelling is not the analytic one.** `eigh` returns eigenvectors by ascending eigenvalue and with arbitrary signs, while ε_k = 2cos(kπ/M) decreases in k. Checking an analytic formula against column a − 1 therefore compares two different modes, and the check fails for a reason that has nothing to do with the formula. Both paths are individually safe, because each is internally consistent and the mismatch is the chiral relabelling composed with a per-mode gauge, neither of which moves a spectrum; mixing them is not. The fix that made both trustworthy was to tie the integer path to the float path by their **spectra**, which agree to 5e−15, rather than by their matrices, which do not.

## What is open

1. **The attainment, and its count.** min 𝒦 ≥ ℓ(N − ℓ)/(N + 1) is proved (Section 7 of the proof document, F144). That it is an EQUALITY needs a saturating vector at each ℓ; at ℓ = 2 the closed-form family below supplies one and is checked per N, and above that it is measured. The equality is not needed for the ceiling, so this is a question about the law rather than about the band. Gate check V12 reads two things that point at where an argument would start. The bound IS attained on every rung it reaches, and the saturating space is exactly what Corollary 7.3 describes, the vectors annihilated by BOTH ladders that also satisfy the two chiral defect conditions: all three quantities vanish on it to machine zero. And the DIMENSION of that space is C(⌊N/2⌋, ℓ) only up to ℓ = 3; from ℓ = 4 it is strictly larger (3 against 1 at N = 8 and N = 9, 15 against 5 at N = 10, ℓ = 4, and 6 against 1 at N = 10, ℓ = 5). So the ℓ = 2 count below does not extend by the obvious guess, and what the extra maximizers are made of is open.
2. **The larger ℓ at larger N.** The closed form is verified on every rung to N = 14 and on some rungs to N = 19. Past that the binding case ℓ = 2 is the one that matters and is the cheapest to extend.
3. **A C# witness.** By the repo's fifth cockpit rule this certificate is a Python verifier that has outlived its session, so it is a witness waiting to be ported; the exact GF(q) rank machinery it needs already exists in `Divisor.RankModP`, and only an integer builder for the reduced object is missing.
