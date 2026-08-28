# The pairing condition

**2026-08-26.** The dephasing palindrome has one condition, and it is not the
two clauses F138 was minted with. On every grid tested, the characteristic
polynomial of the Lindbladian satisfies p(x) ≡ p(−x − 2σ) exactly when a single
operator exists:

> There is an invertible U with **[U, H] = 0** and **U A_l U⁻¹ = −A_l** for every
> dephasing site l.

**What U is.** Not only an algebraic device. The two conditions give
`L(U) = −2σ·U` outright, so U is an exact eigenmode of the Liouvillian at −2σ:
the mirror partner of the steady state at 0, sitting at the far end of the very
axis the palindrome folds about. Read that way the criterion says the spectrum
pairs about −σ exactly when the turning has a conserved quantity that every
dephasing site charges in full.

Scored against the exact GF(p) verdict in both directions, false negatives as
loud as false positives: **0 and 0** on the five full N=3 grids, on the whole
bond-letter axis where F138's converse fails hardest, under non-uniform J and γ,
past F138's ceiling of two dephasing axes, and at N=4 on the path and the ring.

**Since 2026-08-28 the criterion is a theorem, registered as
[F158](../docs/ANALYTICAL_FORMULAS.md), and it is no longer an existence
statement.** U exists exactly when two nullities agree, `dim ker L =
dim ker(L + 2σ)`, and both directions of that are proved in
[PROOF_PALINDROME_TWO_END_COUNT.md](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md).
The far-end eigenspace IS the set of admissible U (an equality-case argument in
the Hilbert-Schmidt norm), the near-end kernel is the commutant, both eigenvalues
are semisimple, and an invertible element exists exactly when the two dimensions
match. Everything below stands as written and was the route; three of its fences
have moved and each says so where it stands.

**Read the left side precisely.** What is measured is the identity of
characteristic polynomials, that is the eigenvalue **multiset** with
multiplicities. It is one notch stronger than the set statement the palindrome
is usually written as, and one notch weaker than F1's operator identity, which
also sees the Jordan structure. `F1PalindromeIdentity` makes exactly this
distinction and scopes it to F138.

---

## What the repo already held, store by store

This is the Stage-0 record, and it is the reason this page is small.

**The field name was already here, and this page is the reconciliation of two
sentences we had written and not joined.**
[KMS_DETAILED_BALANCE.md](../docs/KMS_DETAILED_BALANCE.md) records that Π is
closest to **Q₋, anti-pseudo-Hermiticity**, which asks for an operator with
Q·L·Q⁻¹ = −L†, and notes *"but our condition has L, not L†"*. The one-sided
ℒ_U below is precisely such a Q₋ operator, and the −L† it produces is the step
that was missing. And
[SYMMETRY_CENSUS.md](SYMMETRY_CENSUS.md) carries **both** conditions of the
criterion in the wild, months old: *"UHU† = H … UL_kU† = −L_k … but the
dissipator is quadratic in L_k, so the overall Liouvillian is invariant."* That
sentence is correct and it is the reason nothing followed from it: under
**two-sided** conjugation the sign cancels and U is a symmetry. Under **one-sided**
multiplication it does not cancel, and the same U reflects. The census wrote
down Π² and read it as invariance, which is the same reading a session repeated
on 2026-08-26 before this page existed.

The **open arc** `f138_converse_failures`
([OpenArcsRegistry.cs](../compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs))
asked since 2026-08-03 for a reflector **S** with S L S⁻¹ = −L − 2σ, named the
row to use, and recorded that the bare site reflection was refuted from below.
It did not ask for U; the criterion below is a different object, and **as the
arc stood** it contained nothing of the form [U,H] = 0. It does now, because
this work put it there, and the arc file is modified in the same change that
carries this page: a reader checking at HEAD will find [U,H] = 0 there twice and
should not read that as this sentence being wrong. What this page owes the arc
is the question and the row, not the answer.

**[ANALYTICAL_FORMULAS.md](../docs/ANALYTICAL_FORMULAS.md) F138** held the law
being tested, its converse withdrawn on 2026-08-03, and every count quoted here:
22 and 104 at two bond letters, 776 / 732 / 520 at one, 78 of 21,952 at the full
bond with equal end magnitudes, 78 per exchangeable pair on K₃ and 234 for all
three. It also holds the open item *a derivation of the two-term proviso*.

**[PROOF_PI_FACTORS_AS_R_TIMES_D.md](../docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md)**
(F118) held the factorization Π_Z = R·D and the division of labour that turns
out to be the general one: the transpose leg flips the Hamiltonian part, the
one-sided multiplication reflects the dissipator and carries the entire −2Σγ. It
asserts in a parenthesis, unproved, that any palindromizer must exchange the
dephasing-dead letters with the lit ones.
**[PROOF_QUDIT_PARTIAL_PALINDROME.md](../docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md)**
(F121) proves a class-swap requirement, but only for per-site **product** mirrors
intertwining the **dissipator** palindrome, and its own section 7 records that
non-product intertwiners beat that cap. Precedent for the shape, neither one the
criterion.

**[F1PalindromeIdentity.cs](../compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs)**
held the spectrum-versus-operator distinction, already scoped to F138 and
already saying *the char-poly identity, multiplicities included*.
**[ON_THE_SOFT_BREAK.md](../reflections/ON_THE_SOFT_BREAK.md)** held the repo's
name for that gap and a predictor of the same kind for a neighbouring family.
**[THE_PALINDROME_CLASSIFIER.md](THE_PALINDROME_CLASSIFIER.md)** held the honest
statement that asking whether some operator exists is a search, and thresholded
a spectral quantity anyway.

**[CAUGHT_ERRORS.md](../docs/CAUGHT_ERRORS.md)** returned the lesson the
exhaustive-search section below re-derives independently, dated 2026-08-06:
*a search built from the known instances measures the search, not the codebase.*
It has no entry on conjugator searches as such.

Stores that returned **nothing**: `docs/GLOSSARY.md` has no headword for the
criterion, the class swap, or a weak symmetry; `fw.Confirmations` has nothing,
being a hardware registry; `hypotheses/` and `recovered/` returned nothing.

**And one store answered that had not been asked**: three scouts of 2026-08-04,
gitignored under the WIP rule and therefore invisible to every sweep, had
already built most of a search for S. They left no artifact and their only
surviving verdict was one sentence in a docstring. A local scout is not the repo
answering.

---

## The criterion, and why it is decidable rather than searchable

Write the generator as

    L(ρ) = −i[H, ρ] + Σ_l γ_l (A_l ρ A_l − ρ),   σ = Σ_l γ_l

with every A_l a Pauli letter on one site, so A_l² = 1 and Σ_l c_l† c_l = σ·1
automatically. That premise is doing real work and is named again under Scope.

Both conditions on U are **linear**: [U, H] = 0 and U A_l + A_l U = 0 are linear
systems in the entries of U. The admissible U therefore form a subspace whose
dimension is a rank over GF(p), and the predicate is *the subspace contains an
invertible element*. No catalogue of candidate operators, no group to enumerate,
no eigensolver.

That is not a stylistic choice. A catalogue version of the same predicate, U
ranging over graph automorphisms times Pauli strings, agreed everywhere except
on rows whose two field letters differ, which no Pauli can rotate into one
another; the linear solve resolves all of them. A search space is a hypothesis
about the answer, and this repo already owns the lesson: *a search built from the
known instances measures the search, not the codebase.*

### What U buys, in two steps, and only the first is unconditional

With [U,H] = 0 and U A_l U⁻¹ = −A_l, one-sided multiplication ℒ_U(ρ) = Uρ gives

    ℒ_U L ℒ_U⁻¹ = −L† − 2σ

directly. Only the LEFT factor is conjugated: the commutator term becomes
U H U⁻¹ ρ − ρ H = Hρ − ρH and is untouched, the dissipator's U A U⁻¹ ρ A becomes
−AρA, and the surviving −σρ against the +σρ of −L† is where the −2σ comes from.
This is the Q₋ shape KMS_DETAILED_BALANCE names.

**The step from −L† to the palindrome is not the transpose, and it is not
modular.** spec(−L† − 2σ) = −conj(spec L) − 2σ, a reflection in the vertical line
Re = −σ. It collapses onto the palindrome λ ↦ −λ − 2σ because **L preserves
hermiticity**, so its characteristic polynomial is real and its spectrum is
conjugation-closed. That is the whole sufficiency argument, it holds with no
condition on H, and it can never be seen in GF(p), where there is no conjugation.
It is stated here and gated nowhere, deliberately.

When H is additionally real symmetric in the computational basis, the transpose
upgrades the metric relation to a **similarity**:

    Π = ℒ_U ∘ T,   Π(ρ) = U ρᵀ,   Π L Π⁻¹ = −L − 2σ

and Π² = Ad_U, a commuting operator. **A Y field breaks the precondition**, since
Yᵀ = −Y, and 37 of the 64 field patterns on the scored grid carry one. There Π
is not a reflector: gate 4 exhibits a row where U exists, the metric relation
holds, Π misses by 128 entries, and the palindrome holds regardless. The
similarity is the special case; the metric relation is the law.

That Π² = Ad_U is commuting is an identity, not a measurement: Π = (U⊗I)·T and
T(U⊗I)T = I⊗U, so Π² = U⊗U for any invertible U at all. The **reading** is what
matters and it cost this arc three weeks: **a commuting weak symmetry is
frequently the square of the reflector, not the wrong kind of object.**

---

## Why an exhaustive search returned nothing, and why that is consistent

An exhaustive search over monomial S in the Pauli-string basis, with the
permutation part site-local, found none on the exception rows. That computation
is not part of this page's gate and is not evidence here; what matters is a
structural fact that needs no search. The reflector Π = ℒ_U ∘ T with
U = SWAP₀₂·Z₀Z₁Z₂ is **not monomial in the Pauli basis**: a SWAP is a Clifford,
so it permutes Pauli strings under *conjugation*, while multiplied on one side it
turns a string into a sum of strings.

This also settles an apparent contradiction with the arc, which records that *the
site reversal dressed with signed axis permutations cannot be S at any γ > 0 on
any row*, since the identity string is a null vector of L and every algebra
automorphism fixes it. That argument is correct and Π escapes it twice over: ℒ_U
is a one-sided multiplication rather than an automorphism, and it moves the
identity, Π(1) = U, which is exactly the eigenmode at −2σ. The arc's argument
bounds the monomial family; the answer was leaning on its wall.

---

## What the exceptions were

F138's exceptions at equal field magnitudes are the rows where **U is allowed to
use the graph**. On P₃ with the two end magnitudes equal, U = SWAP₀₂·Z⊗³ works
because the site reversal and the global z-rotation each flip the field
h(X₀ − X₂) and the two flips cancel; neither leg commutes with H alone. That the
exceptions track the graph's automorphisms is measurable on the row: breaking the
reversal by J = (0.43, 0.55) empties the subspace and the palindrome goes with it,
restoring it with J = (0.43, 0.43) brings both back.

F138's clauses are the same criterion read under a restriction:

| F138 clause | what it is, read through the criterion |
|---|---|
| at most two dephasing axes per component | not a property of the jumps at all: three axes on three different sites have a common anticommutant, and the criterion carries no ceiling. What binds is [U,H] = 0 alone, and the gate scores rows with three distinct axes exactly |
| the field has one common axis, orthogonal to every dephasing axis | U must commute with the field. Restricted to U with no site-permutation part, that forces the global π-rotation about the field axis, and orthogonality to the dephasing letters is what makes it anticommute with them |
| the two-term proviso | measured, not asserted: at a single-letter bond H keeps enough symmetry that a U survives where the three-letter bond kills it. The gate scores all seven bond sets |

The demonstration above is one 6-family row. That U exists for all 78 is the
gate's aggregate FN = 0; no operator is exhibited for the other 72, and the
mechanism there is not shown.

---

## Evidence

Gate suite: [`simulations/f138_pairing_condition.py`](../simulations/f138_pairing_condition.py)
→ [`f138_pairing_condition.txt`](../simulations/results/f138_pairing_condition.txt),
committed in this change. It runs on the exact GF(p) kernel
[`f138_exact_palindrome_test.py`](../simulations/f138_exact_palindrome_test.py)
and the committed generator in
[`f138_clause_two_sweep.py`](../simulations/f138_clause_two_sweep.py), and
re-implements the Hamiltonian rather than importing it, so gates 2 to 4 would
catch a disagreement on the named row.

**The named row, assembled once.** N = 3, graph P₃ with edges (0,1) and (1,2),
each bond J(XX+YY+ZZ) at J = 1, dephasing along X on site 1 alone at γ = 1/20,
an on-site X field on sites 0 and 2 of equal magnitude 30/100 with signs
(+, +, −), site 1 fieldless. Its U is SWAP₀₂·Z₀Z₁Z₂, and the anticommutation
that makes it work is Z₁ against A = X₁.

All **42 gates pass**.

| what is scored | rows | false positives | false negatives |
|---|---:|---:|---:|
| P₃, end magnitudes equal | 21,952 | 0 | 0 |
| K₃, all magnitudes equal | 21,952 | 0 | 0 |
| bond plus isolated site, all equal | 21,952 | 0 | 0 |
| P₃, committed tuple (30, 22, 41) | 21,952 | 0 | 0 |
| K₃, committed tuple | 21,952 | 0 | 0 |
| the seven bond-letter sets, 4,096 each | 28,672 | 0 | 0 |
| non-uniform J per bond and γ per site | 1,200 | 0 | 0 |
| two distinct dephasing axes | 626 | 0 | 0 |
| three distinct dephasing axes, past the F138 ceiling | 123 | 0 | 0 |
| N=4 path, equal and generic magnitudes | 240 | 0 | 0 |
| N=4 ring, equal magnitudes | 120 | 0 | 0 |
| N=4 path, non-uniform J | 120 | 0 | 0 |
| **total** | **140,861** | **0** | **0** |

Rows whose admissible subspace was nonempty yet yielded only singular draws,
reported in full because they were where the residual doubt sat: 174, 270 and
1,332 on the three equal-magnitude grids, none at the committed tuple, 18 at
each of the three single-letter bonds, 9 in the non-uniform J and γ row, and
none anywhere else. Two corrections. **These rows are confirmed nonempty at ONE
prime, not at three**: `admissible_multi` returns as soon as a prime yields no
invertible draw, so the later primes are never consulted on them, and the
sentence that stood here said otherwise. **And the doubt is gone**: on every
such row `dim 𝒲 < dim 𝒩`, and
[the two-end proof](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md) §(d) shows
that a strict inequality is a proof that no invertible element exists. What was
sampled is now decided.

**The two verdicts are not alike.** The palindrome verdict is proved on a break
and certified over three primes on a hold. The criterion verdict is **exact when
the admissible subspace is empty**, because a nullity read mod p can only come
out too large; a nonempty subspace can be an artifact of the modulus and is
confirmed at all three primes before a row counts as admissible. Invertibility
inside a nonempty subspace is decided by sampling, with a degree-2^N determinant:
below 10⁻²⁴ at N = 3 and below 5·10⁻²⁴ at N = 4 for three independent singular
draws on a subspace that holds an invertible element. Rows returning
nonempty-but-singular are counted and printed, and the residual doubt there
points at a hidden **false positive**, which is the load-bearing direction.

**One prediction the criterion makes before any of this is run.** Neither
condition on U mentions γ. The criterion is γ-blind by construction, so it
predicts that the palindrome cannot depend on the γ-profile at all. The
non-uniform-γ rows are that prediction under test, not another axis.

---

## Scope, and what would falsify it

Tested: N = 3 and N = 4; chain, ring, triangle, and a disconnected bond plus
isolated site; all seven bond-letter sets; uniform and non-uniform rational J and
γ; one, two and three distinct dephasing axes; Pauli-letter jump operators only.

**Not tested.** An anisotropy axis: every bond carries weight 1 per letter and
there is no Δ. N beyond 4 for the criterion as scored on this page.

**Off-axis dephasing has since been both proved and tested, and it was never a
new axis.** The theorem asks of a jump operator only that it be Hermitian and
square to the identity, which any `n̂·σ⃗` at a unit direction satisfies, and so
does any multi-site Pauli string. Both are scored in
[`f138_rank_criterion.py`](../simulations/f138_rank_criterion.py), 1,692 rows
with rational directions and string jumps, FP = 0 and FN = 0, and the sibling
arc `f138_clause_two_sweep`'s question about clause 2's word *orthogonal* is
answered there: orthogonality of directions is not what the criterion asks for.

**And one axis where the criterion is not merely untested but mis-centred.**
Non-Pauli jump operators break Σ c†c = σ·1, and the palindrome does not simply
die there: **F137** records that T1 alone keeps it, about **−σ/2** rather than −σ.
So a criterion phrased "pairs about −σ" is answering the wrong question in that
regime, and the premise is doing more work than a footnote. F82 to F84 own the
separate matter of Π failing there.

**Necessity rested on 0 false negatives, and now it does not.** When this page
was written, sufficiency was a short calculation and necessity was measured on
the grids above and nowhere derived. It is derived in
[PROOF_PALINDROME_TWO_END_COUNT.md](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md)
§(e), and the route is not this page's: not through the operator at all, but
through the two ends of the axis. L is unital, so 0 is an eigenvalue; the
palindrome carries its multiplicity to −2σ; both eigenvalues are semisimple; so
the two kernels have equal dimension, and equal dimension is what forces an
invertible U. The census is now an illustration rather than the evidence.

What would falsify it: one row, at any setting, where an invertible U exists and
the characteristic polynomial does not pair.

---

## What this leaves open, and what it makes stale

**The typed layer.** Nothing here is citable from the repository until the
construction lands as a claim with a live witness. The gate script is committed;
the witness is not written.

**The arc `f138_converse_failures` now needs a pass** (and got one, in this same
change). Its ParkedAt says the named candidate space is empty and that existence
is free where diagonalizability holds, and it grades its own evidence there more
carefully than this sentence first did: *"CERTIFIED ONLY: the named row measures
30 of 64 … diagonalizability rests on g(L) == 0 mod p at three primes, which is
not a decision … treat this half as unfinished rather than as measured."* So it
is not "every row tested", and the arc says so itself. Both remain true, and both are now beside the point:
the operator exists, it is structured, and it was found by a linear solve scored
in both directions rather than inferred from cospectrality. Its finding (5) files
V = R·U as *not the answer, for a reason of kind*; V is Ad_U, which is Π², so
that is the corrected reading and not the current one.

**F138 carried one stale sentence, and it is retired in this same change**:
*"No operator explaining the exceptions has been exhibited, and the first
candidate for one, the bare site reflection, was refuted from below."* Its first
clause was false on the 6-family, and F138's row now carries two exhibited
operators instead, with the reason the old candidate failed: U is a SUM of two
Pauli strings, neither a single string nor a signed site permutation, so the
refutation was right about its candidate and was never evidence that nothing
exists. The row keeps the separation the arc prescribed, EXISTS from EXPLAINS.
F138's headline, that the law is an implication whose converse is open, stays
exactly as it is; the criterion is not F138's two clauses.

**One surface this page does not fix, and one that is fixed in this change.**
[MIRROR_SYMMETRY_PROOF.md](../docs/proofs/MIRROR_SYMMETRY_PROOF.md) said the
palindrome holds *exactly when* and carried no converse fence, although F138
names its Scope paragraphs as its Proof anchor; the word is now *when*, with the
converse failure and its counts beside it, and F138's Proof field is narrowed to
say those paragraphs carry the sufficient direction as a measured census. Still
unfixed:
[PROOF_F111_HARD_CELL_PURE_D_TEMPLATE.md](../docs/proofs/PROOF_F111_HARD_CELL_PURE_D_TEMPLATE.md)
asserts in an aside that spectrum-level palindromy is realized by some similarity
with existence guaranteed by the palindromic spectrum, which is free only where
the generator is diagonalizable.

**The neighbouring predictor.** The bipartite criterion in
[ON_THE_SOFT_BREAK.md](../reflections/ON_THE_SOFT_BREAK.md), with its windowed
converse in [PROOF_F103_F87_Z2_CUBED_REFINEMENT.md](../docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md),
answers a question of the same shape for a different family and has never been
run against these rows.

---

## How this was found, because the route is part of the result

The derivation was not wrong anywhere. It established, exactly, that V = R·U
commutes with L rather than reflecting it, and filed V as the wrong kind of
object. The step that was missing was one question: what does the field call the
relation between an operator that commutes with a generator and one that reflects
its spectrum. The answer, that the first is a weak symmetry and that squaring a
reflector always produces one, turns V from a dead end into the square of the
answer. And the field name was not even outside: two of our own files carried it,
unjoined.

That question is now a lens, `naming-what-we-already-own`, in the local
superpowers fork. It points INWARD and is fenced against outward searches: the
name we needed was in two of our own files, and joining them is the move. Its
acceptance test is that a name handing over no theorem, no method and no new
question did not fire, and a renaming is not a finding.
