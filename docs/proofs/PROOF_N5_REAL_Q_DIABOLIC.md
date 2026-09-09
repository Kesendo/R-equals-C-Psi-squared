# The N=5 real-q diabolic is semisimple, by Galois conjugacy

**Tier 1, derived.** The A₂ locus at q = ±1.129 250 970 874 767 1, λ = −4.791 960 365 179 641 0
carries a **semisimple** double eigenvalue. The proof is a field argument over two
facts the repo already certifies exactly. The proof itself uses no eigensolver and
no tolerance; the gate beside it is numerical, and says so.

**The repo classifies this locus and does not certify it, and it says so
consistently.** The verdict is gated in
`compute/RCPsiSquared.Diagnostics.Tests/Foundation/RouteBA2PositiveAnchorTests.cs`,
at these very numbers (`WAnchor = 5.10083102`, `LambdaAnchor = -4.79196037`, both q
lifts, R-odd), asserting `EpKind.Diabolic` with `Algebraic == 2` and
`Geometric == 2`. That gate is floating point throughout: `distances[1] < 1e-4`,
relative departure `< 1e-6`, dephase departure `< 1e-6`. `docs/CAUGHT_ERRORS.md`
(2026-09-06, the same day) caught an inference at this locus, that an order-two
discriminant zero gives semisimplicity, and closed with the sentence that
"`EpCharacter` on the isolated local operator or exact rank must determine its
**geometric multiplicity**", and its "Current boundary" paragraph then states that
"the positive-real R-odd A2 root and 16 nonreal w roots (34 q-loci together) remain
**character-uncertified**. Route B remains open until all 29 w roots / 58 q-loci
reconcile" (lines 2675-2680).

That is not in tension with the stores that call the locus classified. The two words
carry different weight here, and `docs/THE_DOUBLE_ROOT.md:348-351` uses both in one
breath about the same 34 loci: "all 34 stable EpCharacter readings are Diabolic,
alg=geo=2, on the full sector at three isolating radii. **No exact fallback was
used; the artifact exactRankCertificates array is empty.**" CLASSIFIED means a stable
numerical `EpCharacter` reading; CERTIFIED means an exact one. Every store agrees
that this locus has the first and lacks the second. What follows supplies the second,
for the 2 real-q loci of those 34; the other 32 keep only the reading.

The locus: the open N=5 XY chain under uniform local Z dephasing, in the (1,2)
coherence block, R-odd sector, at Δ = 0 and γ = 1 per site. Both q and λ are real,
so this is an operating point one can set, not an analytic continuation. It is the
only such point in the Route-B corpus at these two N: at N=5, 2 of the 58 q-loci
are strictly real and they are ±the same point; at N=6, counting the stored atlas
by `tReal == 0`, **0 of 266** are, with 118 at purely imaginary q where the block
is Hermitian and 148 generic complex.

## The sweep, stores by name

`experiments/F89_PATH_K_DIABOLIC.md` holds this locus, and it holds **both** things.
Lines 126 and 137 flag a gap: "This does not settle the positive-real R-odd A2
locus at N=5", "These R-even readings do not decide the R-odd A2 character at
N=5". The gap is one of scope, because the reflection-parity mechanism that
grounds the odd-N real-q diabolics was read entirely on the **R-even** sector and
this locus is R-odd. Lines 3 and 283 then record the answer that arrived by the
other route; line 283 has it in full, "both q lifts of the positive-real R-odd A2
root are Diabolic, alg=geo=2, by stable full-sector EpCharacter readings", and line
3 the shorter form, "classifies both lifts of the positive-real R-odd A2 root as
diabolic". The same file supplies the
defective control used below, the R-odd real EP at q ≈ 2.804888, λ ≈ −4.4882
(line 146).

`docs/THE_DOUBLE_ROOT.md` holds the N=5 layer certificate (both layers irreducible)
and, in its N=5 character table at lines 332-337 and the prose at 348-351, the same
verdict. `compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs` carries the
mechanism under `diabolic_over_higher_n` and its walls, and the verdict again at
line 274, "the positive R-odd A2 pair is already classified Diabolic with
alg=geo=2". The same arc entry states the gap this document closes, at lines 58-59:
"No exact fallback was used; the artifact exactRankCertificates array is empty.
PSC1 supplies algebraic pair uniqueness and S1 supplies the repeated-lambda seed;
**neither determines Jordan character.**"

`simulations/results/route_b_a2_n5.json` holds the inventory: 29 A₂(w) roots, 58
q-loci, exact isolation boxes. Its `exactRankCertificates` list is **empty**, and
the root records carry no character field at all; the N=5 verdict is produced at
run time by `RouteBA2CharacterClassifier`, not stored.

`simulations/o2b_gcd_certificate.py` supplies both inputs of the proof, and
supplies them **executed** over ℤ rather than asserted: line 1336 factors the exact
A₂ layer and asserts irreducibility over ℚ, lines 1342-1343 assert the exact root
inventory (R-odd: 6 negative real, 1 positive real, 6 nonreal). Its item 3b
excludes, at every nonzero w real or complex, both a second double λ-root and a
triple root. One reading note on the irreducibility assert, since the proof leans
on it: a single entry in `factor_list` would in general also admit A₂ = f², but A₂
is the Yun multiplicity-2 layer of a squarefree decomposition and is squarefree by
construction, so the single entry does mean irreducible.

`docs/ANALYTICAL_FORMULAS.md` holds F163, and its proof document
`docs/proofs/PROOF_ROUTE_B_N6_UNFOLDING.md` carries the conjugacy argument this one
adapts: §2 certifies the 133 even loci from one Hermitian conjugate, §5 the odd
sector. It is in the direct-t parameter, where the factor already has ODD degree 133
and therefore hands the argument a real root for free. That free real root is exactly
what the N=5 lift below has to buy. Neither document makes an N=5 character claim. F150 and F151 hold the
fold-checkerboard and gauge criteria that decide when a sector discriminant is
real; neither decides a Jordan type.

The typed layer, both halves, plus the tests. `DiabolicReflectionParityWitness`
(`inspect --root diabolicparity`) is live and grounds the odd-N mechanism, but reads
the **R-even** sector across N=5..9, so it does not reach this locus. Diagnostics
also carries `RouteBA2Inventory` and `RouteBA2CharacterClassifier`. No `Claim` and
no `IInspectable` owns this character; the C# gate named at the top does, and it is
a test rather than either half of the typed layer, which is why an earlier draft of
this sweep missed it.

**Step 2 below is prior work, and saying so strengthens it.** The Hermitian anchors
are not asserted here for the first time: `experiments/F89_PATH_K_DIABOLIC.md:145`
already reads the R-odd negative-w family as "semisimple by Hermiticity exactly as on
the R-even side (the block is real-symmetric at Re q = 0)", and
`docs/CAUGHT_ERRORS.md:2676-2677` states the same with the generator written out,
"M(q)=D/2+iqK/2 there, with D real diagonal and K real symmetric". What is new is not
the anchors but that a single Galois orbit connects them to the physical point.

Two further stores, swept and returning a different object.
`experiments/F89_PATH_K_GALOIS.md` is the repo's Galois file, but it certifies that
the path-3 octic's group is S₈, a negative result about radicals, and carries no
conjugacy transfer; its "isolate before DDF" degree-pollution principle is, however,
the same care the χ = AT·F_res remark below takes.
`docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md` §8 is the standing semisimplicity
instrument at this block, but it is Theorem A's two regimes of silence, a different
mechanism.

`docs/GLOSSARY.md` returns nothing on A₂ loci or Route B, but carries the q versus
Q factor-2 trap at line 339 that the gate's anchor had to pin. `fw.Confirmations`
(24 entries) returns one EP entry, `ibm_ep_onset_may2026`, an IBM Kingston
population trajectory from 2026-05-31; it is not a reading of this locus. No
hardware measurement of this point exists.

What is new is one fact and its consequence: that A₂_O(−t²) is irreducible over ℚ,
so that the 26 q-lifts of the 13 A₂_O roots form a **single** Galois orbit, and that
this transfers semisimplicity from the Hermitian anchors to the physical point. The
number 26 itself is not new; `docs/THE_DOUBLE_ROOT.md:337` has carried the R-odd row
as 13 roots and 26 loci since the inventory closed. Nor is the argument's shape: it
is the N=6 conjugacy argument, applied one lift step lower down.

## What the block is

The (1,2) coherence block of the open N-site XY chain under uniform Z dephasing has
dimension N·C(N,2) = 50 at N=5. The site reflection R commutes with the generator,
exactly, and splits it into an R-even sector of dimension 26 and an R-odd sector of
dimension 24. Write t = i·q_unitHop, where the inventory's parameter is
w = q_unitHop² and q_unitHop = 2·q_physical, so that

    w = −t².

The whole argument runs **inside the R-odd sector**, where

    L_O(t) = D_O + t·K_O

with D_O the rational dephasing diagonal (entries −2 and −6 at γ = 1) and K_O the
rational hopping, integer in the ±1 orbit basis.

**Why the sector and not the full block.** Not for Hermiticity: the full 50×50 block
is already real symmetric at real t, measured `max|L − Lᵀ| = 0.0`. Nor because the
even sector fails: K_E's asymmetry in the ±1 orbit basis is only the orbit-size
normalisation (E mixes column norms 1 and √2 where O is uniformly √2), and in an
orthonormal basis **both** restrictions are exactly symmetric, measured 0.0. The
reason is that
**A₂_O is a per-sector object**: the discriminant layer the whole argument is built
on is the R-odd sector's, and there is no full-block polynomial whose irreducibility
would carry the same 13 roots. The sector is where the certified input lives, and
K_O being integer and symmetric there is what lets steps 2 and 3 consume it
directly.

## The theorem

**Let w₀ be the positive real root of the R-odd A₂ layer at N=5. The double
eigenvalue of L_O at w₀ is semisimple.**

*Proof.* Three steps.

**1. One Galois orbit.** A₂_O is irreducible over ℚ of degree 13, so
[ℚ(w₀) : ℚ] = 13, and since w₀ is a real number, ℚ(w₀) is a subfield of ℝ. The
corresponding parameter is t₀ = √(−w₀), which is **nonzero** purely imaginary
because w₀ > 0 and A₂_O does not vanish at w = 0. No nonzero purely imaginary
number lies in a subfield of ℝ, so t₀ ∉ ℚ(w₀) and [ℚ(t₀) : ℚ(w₀)] = 2, giving
[ℚ(t₀) : ℚ] = 26. The minimal polynomial of t₀ therefore has degree 26 and divides
A₂_O(−t²), which also has degree 26, so the two agree up to a constant: A₂_O(−t²)
is irreducible over ℚ and its 26 roots form a **single** Galois orbit. The orbit
accounts for every root: 1 positive w gives 2 imaginary t, 6 negative w give 12
real t, 6 nonreal w give 12 nonreal t.

**2. The orbit meets the Hermitian axis.** A₂_O has 6 negative real roots. For each,
−w > 0, so the corresponding t are **real**, and there L_O = D_O + t·K_O is a real
symmetric matrix, hence Hermitian, hence geometric multiplicity equals algebraic
multiplicity for every one of its eigenvalues. Twelve of the orbit's 26 members are
such points. What step 1 supplies is not the parity of 13 but the fact that the
single orbit contains **both** signs of w: the positive one is the physical point,
the negative ones are the Hermitian anchors, and both are in the certified root
inventory.

**3. Transfer.** Let M be the Galois closure of ℚ(t₀, λ₀) and E ⊆ M the splitting
field of A₂_O(−t²), the containment holding because M contains t₀ and, being normal,
every conjugate of it. Since E/ℚ is normal, the restriction Gal(M/ℚ) → Gal(E/ℚ) is
surjective, so the automorphism of E carrying t₀ to a real member of the orbit
lifts to some σ ∈ Gal(M/ℚ). Because L_O = D_O + t·K_O with D_O and K_O over ℚ,
applying σ entry-wise gives σ(L_O(t₀)) = L_O(σt₀), and a field automorphism applied
entry-wise preserves rank, so

    geo(t₀, λ₀) = geo(σt₀, σλ₀).

At σt₀ the matrix is Hermitian by step 2, so geo(σt₀, σλ₀) = alg(σt₀, σλ₀). And σ
permutes the roots of the characteristic polynomial, which lies in ℚ[λ, t],
preserving their multiplicities, so alg(σt₀, σλ₀) = alg(t₀, λ₀). Chaining,

    geo(t₀, λ₀) = alg(t₀, λ₀),

which is semisimplicity at w₀. ∎

The chain never needs the *value* of the multiplicity, only that the two agree, and
that is what makes it safe: the characteristic polynomial of the sector is
χ = AT·F_res, and the certificate of `o2b_gcd_certificate.py` item 3b bounds
repeated roots of F_res only, while AT's non-constant factors have real roots that
sweep the spectrum at negative w. Asserting "algebraic multiplicity 2" at the
transferred point would have leaned on a certificate that does not reach there.
Read at t₀ instead, the value follows as a separate remark: AT does not vanish at
(λ₀, w₀), and item 3b gives F_res exactly one double root, so alg(t₀, λ₀) = 2, in
agreement with the C# gate's reading.

## What is held, and what is not

Held: the geometric multiplicity at the positive real R-odd A₂ root of N=5, for the
uniform chain at γ = 1 per site, Δ = 0.

Not held. Nothing about the other 57 N=5 q-loci beyond what the Hermitian axis
already gives, nothing about the R-even sector's scan, nothing about any N ≥ 6. No
ε radius and no perturbative statement; F163's unfolding is a separate object.
Nothing at Δ ≠ 0, where the ZZ term breaks the free-fermion structure and the A₂
layer is a different polynomial. No metrological claim and no hardware proposal.

The remark that alg = 2 leans on a third input beyond the two below, namely that AT
does not vanish at (λ₀, w₀). That one is read off the repo's own description of the
rungs rather than a factorisation: `experiments/F89_PATH_K_DIABOLIC.md:143` gives the
R-odd AT strands as λ = r₀ + 2isq with r₀ ∈ {−2, −6}, which at real q is non-real
unless s = 0, and at s = 0 is −2 or −6. λ₀ is real and is neither. The
semisimplicity statement itself does not use this input at all.

Two dependencies, named plainly because the argument rests on them. The
irreducibility of A₂_O over ℚ and its negative-real root count are taken from
`o2b_gcd_certificate.py`, not re-derived here; they are executed there in exact
integer arithmetic, and the gate re-runs the inventory shape against the committed
artifact, but the irreducibility itself is that file's result. And the exclusion of
a second double root or a triple root at every nonzero w is that file's item 3b.

## The gate

`simulations/n5_a2_rodd_character_gate.py`, five steps. The proof is exact; this
gate is not, and it does not fill the artifact's empty `exactRankCertificates`.

**Anchor.** The inventory shape the proof needs is checked first: A₂_O of degree
13, split 6 negative / 1 positive / 6 nonreal. Then the q/J book, pinned by
requiring exactly one of three scalings to put the stored 49-digit λ₀ in the R-odd
spectrum to 1e-12; it pins to J = q_unitHop with |Δλ| = 1.78e-15, and the R-even
sector is checked to **not** carry λ₀ (distance 1.16), so the parity label is
measured and not assumed. The anchor pins the label only; the factor-2 convention
is pinned by the source step, where a wrong halving turns all six roots red.

**Structure.** Exact reads, all `== 0.0`: L is Hermitian at real t; R commutes
with L **and** is the reflection the two projectors are built from, so the sector
split is the one the proof means and not merely some commuting permutation; and
K_O is exactly symmetric and exactly integer in the ±1 orbit
basis, `max|K_O − K_Oᵀ| = 0.0` and `max|K_O − round(K_O)| = 0.0`, on a K_O whose
largest entry is 1.0, which is the nonzero paired read. Integrality is the property
the proof consumes, and the restriction is taken as the exact halving ½·U_Oᵀ K U_O
rather than by a least-squares solve, which would have measured only rounding.
The physical-axis Hermiticity companion returns 1.4270, which is 2·MAG by
construction, so it reports that the hopping is present rather than discriminating
anything, and is labelled that way.

A fourth check in this step is a **construction** check, on the commutator sign,
and it is there because no spectral check can do the job. Flipping the bra-side
sign conjugates the entire R-odd spectrum, and the A₂ inventory is closed under
conjugation (its nonreal roots come in conjugate pairs), so every anchor available
reads identically either way. The gate measures the blindness rather than asserting
it: at the anchor coupling the flipped build differs from the correct one by
**4.5170** in the full block, yet its R-odd spectrum is the conjugate of the correct
one to residual **0.0**, and λ₀ being real, the three singular values differ by
**exactly 0.0**. This is a real blindness of the numerical apparatus, recorded
rather than papered over.

**Source.** All six negative-real w roots read nullity 2, σ₁ and σ₂ between 8.3e-18
and 1.4e-16 against σ₃ between 6.56e-3 and 1.78e-2.

**Conclusion.** At w₀ the read is nullity 2, σ = 1.71e-17, 7.53e-17, 8.58e-02. Two
controls run on the same code path in the same sector: the defective EP reads
nullity **1** (σ₂ = 1.96e-02), and a generic λ 0.137 away reads nullity **0**. The
repo records the defective point to six digits only, so the gate refines it, and
certifies that the refined point really is defective by the local law: the pair
splits as C·√|J − J*|, measured exponent **0.4997**, where a semisimple coalescence
would return 1.0. The assert keeps a residual-gap bound beside the exponent, but it
is the exponent that discriminates.

**Error model.** The thresholds rest on a measured law rather than on a chosen
number. Detuning the coupling by δJ/J and re-reading at w₀:

| δJ/J | σ₂ at w₀ | σ₂ / (δJ/J) | σ₂ at the defective point |
|---|---|---|---|
| 1e-12 | 5.8395e-14 | 5.8395e-02 | 1.9550e-02 |
| 1e-10 | 5.8393e-12 | 5.8393e-02 | 1.9550e-02 |
| 1e-08 | 5.8393e-10 | 5.8393e-02 | 1.9550e-02 |
| 1e-06 | 5.8392e-08 | 5.8392e-02 | 1.9550e-02 |
| 1e-04 | 5.8385e-06 | 5.8385e-02 | 1.9549e-02 |

The ratio is flat to 1.000 18 across eight decades, so σ₂ = 5.84e-2 · (δJ/J) at w₀ under a uniform detune,
while at the defective point σ₂ stays pinned at the coupling scale and does not
move at all. (the constant c depends on the detune direction; the flatness does not). At zero
detune σ₂ = 7.53e-17 corresponds to δJ/J = 1.29e-15, about 6.6 units in the last
place of J, so the read sits on the rounding floor of the eigen path rather than at
a tuned threshold. Two further things the gate measures rather than restates: across
the same eight decades σ₂ at the defective point stays in [1.9549e-02, 1.9550e-02],
a swing of 1.000 07. That read says NOT semisimple; it is the nullity-1 read above,
not this one, that says defective;
and the margin is not uniform, the σ₃ side having a computed factor of 6.56 of
headroom rather than the fifteen orders the w₀ read alone suggests.

## Links

- [F89 path-k diabolic census](../../experiments/F89_PATH_K_DIABOLIC.md), which
  holds the locus, the R-even scope caveats, and the numerical verdict
- [The double root](../THE_DOUBLE_ROOT.md), the N=5 layer certificate and character
  table
- [Analytical formulas](../ANALYTICAL_FORMULAS.md), F163 for the N=6 conjugacy
  argument this one adapts, F150 and F151 for the reality criteria
- [Caught errors](../CAUGHT_ERRORS.md), 2026-09-06, the inference this replaces
