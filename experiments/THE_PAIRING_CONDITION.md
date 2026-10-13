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

**The criterion is a theorem, registered as
[F158](../docs/ANALYTICAL_FORMULAS.md#f158).** U exists exactly when two nullities agree, `dim ker L =
dim ker(L + 2σ)`, and both directions of that are proved in
[PROOF_PALINDROME_TWO_END_COUNT.md](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md).
The far-end eigenspace IS the set of admissible U (an equality-case argument in
the Hilbert-Schmidt norm), the near-end kernel is the commutant, both eigenvalues
are semisimple, and an invertible element exists exactly when the two dimensions
match. Everything below is the route; where a fence has moved, the passage says
so.

**Read the left side precisely.** What is measured is the identity of
characteristic polynomials, that is the eigenvalue **multiset** with
multiplicities. It is one notch stronger than the set statement the palindrome
is usually written as. Inside this page's class it is not weaker than a
similarity: where the multiset pairs, [F158](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md)
§(f8) gives L ~ −L − 2σ with the whole Jordan structure. What F1 adds is one
particular operator, Π, unitary and in closed form.

---

## What the repo already held, store by store

This is the Stage-0 record, and it is the reason this page is small.

**The field name was already here, and this page is the reconciliation of two
sentences we had written and not joined.**
[KMS_DETAILED_BALANCE.md](../docs/KMS_DETAILED_BALANCE.md) lists
**Q₋, anti-pseudo-Hermiticity**, by name among the symmetries the Sá-Ribeiro-Prosen classification uses; in that
classification it is the symmetry that asks, on the trace-shifted generator
L_c = L + σ, for an operator with Q·L_c·Q⁻¹ = −L_c†, while Π's own relation has −L_c where Q₋ has −L_c†. The one-sided
ℒ_U below has that shape (whether it meets the classification's own conditions
is not checked here), and the −L† it produces is the step
that was missing. And
[SYMMETRY_CENSUS.md](SYMMETRY_CENSUS.md) carries **both** conditions of the
criterion in the wild, months old: *"UHU† = H … UL_kU† = −L_k … but the
dissipator is quadratic in L_k, so the overall Liouvillian is invariant."* That
sentence is correct and it is the reason nothing followed from it: under
**two-sided** conjugation the sign cancels and U is a symmetry. Under **one-sided**
multiplication it does not cancel, and the same U reflects. The census wrote
down Ad_U, the two-sided map Π² equals for an orthogonal U, and read it as
invariance, which is the same reading a session repeated
on 2026-08-26 before this page existed.

The **open arc** `f138_converse_failures`
([OpenArcsRegistry.cs](../compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs))
asked since 2026-08-03 for a reflector **S** with S L S⁻¹ = −L − 2σ, named the
row to use, and recorded that the bare site reflection was refuted from below.
It did not ask for U; the criterion below is a different object. What this page
owes the arc is the question and the row, not the answer.

**[ANALYTICAL_FORMULAS.md](../docs/ANALYTICAL_FORMULAS.md) F138** held the law
being tested, its converse withdrawn on 2026-08-03, and every count quoted here:
22 and 104 at two bond letters, 776 / 732 / 520 at one, 78 of 21,952 at the full
bond with equal end magnitudes, 78 per exchangeable pair on K₃ and 234 for all
three. Its two-term proviso is derived in
[the palindrome as a colouring](THE_PALINDROME_AS_A_COLOURING.md): one-letter bonds of one letter leave the colour free to change from site to site.

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

with every A_l a Pauli letter on one site on this page's grids. The argument
uses only that each A_l is Hermitian with A_l² = 1, so Σ_l c_l† c_l = σ·1; that
premise is doing real work and is named again under Scope.

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

**The step from −L† to the palindrome is not the state transpose ρ ↦ ρᵀ.**
spec(−L† − 2σ) = −conj(spec L) − 2σ, a reflection in the vertical line
Re = −σ. It collapses onto the palindrome λ ↦ −λ − 2σ because **L preserves
hermiticity**: in the Pauli-string basis L is a real matrix, so L† is its matrix
transpose, and every square matrix is similar to its transpose over any field.
Hence L ~ −L† − 2σ ~ −L − 2σ, a similarity and not only an equal spectrum, as
[the two-end count §(f8)](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md)
derives for the Jordan structure. That is the whole sufficiency argument, and it
holds with no condition on H beyond hermiticity. No step needs a conjugation, so
it runs over GF(p): gate 4, on a row with a Y field where the explicit similarity
below fails, checks modulo p the metric relation, that L† is the Pauli-basis
transpose of L (in the matrix-unit basis it is not, for any non-scalar H, which
is what makes the check one that can fail), and the palindrome the chain
concludes; the step M ~ Mᵀ is a theorem and is not computed.

When H is additionally real symmetric in the computational basis, the transpose
upgrades the metric relation to an **explicit similarity**:

    Π = ℒ_U ∘ T,   Π(ρ) = U ρᵀ,   Π L Π⁻¹ = −L − 2σ

and Π²(ρ) = U ρ Uᵀ, a commuting operator, which is Ad_U when U is orthogonal
(X^⊗N and the named row's U are). **A Y field breaks the precondition**, since
Yᵀ = −Y, and 37 of the 64 field patterns on the scored grid carry one. There Π
is not a reflector: gate 4 exhibits a row where U exists, the metric relation
holds, Π misses by 128 entries, and the palindrome holds regardless. The
explicit Π is the special case; the metric relation is the law.

That Π² is the two-sided map ρ ↦ U ρ Uᵀ is an identity, not a measurement: Π = (U⊗I)·T and
T(U⊗I)T = I⊗U, so Π² = U⊗U for any invertible U at all. The **reading** is what
matters and it cost this arc three weeks: **a commuting weak symmetry is
frequently the square of the reflector, not the wrong kind of object.**

### What U buys on the flow

The repository already holds the relation on the flow once (the gamma fold),
two similarity cousins of it (the centred Π fold, the block fold), and beside it
Krein metrics of its own and a twin of ℒ_U. This subsection joins them and adds three things to
the repository: the trace form for every admissible U, written down as a
corollary of the relation; the reading of the conservation law it implies
between two forward runs, a law standard for any pseudo-Hermitian generator,
with ℒ_U in the metric's place; and ℒ_U set beside the repository's other Krein
metrics.

- **The trajectory, for U = X^⊗N.** [The gamma fold](GAMMA_FOLD_PAIR_OF_MIRRORS.md)
  docks onto the lattice as X^⊗N ρ(t) = e^{−2σt}·[gain evolution of X^⊗N ρ(0)]
  (`GammaFold.ReadThroughVeil`, gated in `GammaFoldTests`). For a Hermitian H and
  Hermitian jumps the gain generator is exactly −L†, so that
  identity is ℒ_U e^{tL} ℒ_U⁻¹ = e^{−2σt}·e^{−tL†} (the gain generator
  L(−γ) = L_H − L_D, L_D the dissipator, negates the dephasing: a formal
  generator, not a channel). The relation itself, as an
  operator identity, is in [MIRROR_SYMMETRY_PROOF](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)
  and [the two-end count](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md).
- **The centred fold.** [Π as time reversal](PI_AS_TIME_REVERSAL.md) writes
  L_c = L + σ and, where F1's similarity holds, Π e^{tL_c} Π⁻¹ = e^{−tL_c} with
  μ ↦ −μ for μ = λ + σ, and says what this is not: it does not reverse the
  dissipative envelope or build a realizable reversed trajectory.
- **The block fold.** MirrorWorld's `Mirror.TrajectoryFold` runs the same
  reading block by block, in its own scope (the XY chain, uniform γ): the bra
  fold with the bipartite gauge carries a block (p, q) to its partner (p, N − q),
  a similarity to −L − 2σ there, and the folded trajectory, weighted by the
  price, solves ẇ = −L w on the partner block, "the
  partner world running backward, at the price" in the words of its code
  comment (the fold lattice itself is
  [PROOF_CODIM1_BY_ADDITIVITY](../docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md) §7).
- **Krein metrics and a twin.** [The weight-one degeneracy](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md),
  Lemma 1: on a diagonal block (p, p) (ket popcount p, bra popcount p) of a
  real, number-conserving H under uniform Z dephasing (its argument, τDτ = D,
  carries over to a profile), the transpose τ (the T above)
  makes the block self-adjoint for an indefinite form, so every mode with a
  complex eigenvalue is τ-neutral. [The fragile bridge](../hypotheses/FRAGILE_BRIDGE.md)
  reads its threshold as a Krein collision under the chain reflection G, and
  with [the PT-symmetry analysis](PT_SYMMETRY_ANALYSIS.md) links G to the right
  multiplication R(ρ) = ρ·X^⊗N (here R is this map, not the site reversal R
  of V = R·U under "What this leaves open") through their S = G∘R (not the reflector S of
  the arc above), which commutes with all of L
  when H commutes with X^⊗N and is reflection-even and the rate profile is
  reflection-odd,
  γ_{N−1−l} = −γ_l (so σ = 0, with rates of both signs, a formal generator
  rather than a dephasing channel), and maps block (p, q) to (p, N − q).
  [The β-exotic genericity study](F89_BETA_EXOTIC_GENERICITY.md) holds a
  bipartite Krein metric of the same kind. R is the R of F118's Π_Z = R·D ([the Π factorization](../docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md),
  MirrorWorld's `MirrorGroup`, which calls it the ket reflection; this page, like
  the lattice, counts right multiplication on the bra); it reflects L to
  −L† − 2σ, as F109 records, by the calculation of
  ℒ_U done from the other side. The sources hold R as a reflector, not as a
  form; read as a form, ⟨r, r·X^⊗N⟩ = ⟨r†, X^⊗N r†⟩, R gives a mode the sign
  that ℒ_U gives its conjugate mode, so it adds no sign of its own.

What follows uses the Hilbert-Schmidt product ⟨A, B⟩ = Tr(A†B), writes L† for
the adjoint of L in it (the generator of the Heisenberg picture), and keeps the
name L_c = L + σ. The relation reads ℒ_U L_c ℒ_U⁻¹ = −L_c†. Each mode of L sits
at λ = −α + iω, its decay rate α (the fade) and its angular frequency ω (the
turning), the radial and the angular hand of
[the F95 clock](../compute/MirrorWorld/Clock.cs); in the centred frame it sits
at μ = λ + σ = (σ − α) + iω. The **centre** Re μ = 0 holds the modes that fade at
exactly σ = Σγ. It is the cut that [the Born-rule shadow](BORN_RULE_SHADOW.md)
and [the reflection keeps past and future apart](THE_REFLECTION_KEEPS_PAST_AND_FUTURE_APART.md)
draw and read as the line between past (slower than Σγ) and future (Σγ or
faster): an interpretive label of those pages, not a time ordering, and one
that puts the centre itself on the future side.

**The trace form, for any admissible U.** Pair e^{tL_c} = ℒ_U⁻¹ e^{−tL_c†} ℒ_U
with a readout B, move the adjoint onto the other factor, and use that L
preserves hermiticity, (e^{sL_c}A)† = e^{sL_c}(A†):

    Tr(B · e^{tL_c}(A)) = Tr(e^{−tL_c}(B U⁻¹) · U A).

Here A and B are any operators, not the jumps A_l. The left side is an experiment (prepare A, run forward, read B) times e^{σt}.
The right side is not one. B·U⁻¹ and U·A are in general not Hermitian (for
U = X^⊗N they are the lattice's two one-sided readings, X^⊗N applied on the bra
side of B and on the ket side of A), and e^{−tL_c} is the inverse of a forward
run, not a channel. The same number also has a form with forward runs only,
Tr(UBU⁻¹ · e^{tL_c}(UAU⁻¹)), by the two-sided weak symmetry Ad_U, the
cancellation the census recorded. Where B·U⁻¹ is an eigenoperator of L the
backward run can be evaluated: B = U gives Tr(U·e^{tL}(A)) = e^{−2σt}·Tr(UA),
the far-end decay of the linear sibling below. Taken over all A, B and t the
trace form is the relation itself. Nor is it a
reciprocity in the sense of
[F119 §7](../docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md#7-one-level-up-the-generators-transpose-and-the-two-ends-of-a-channel),
which runs forward on both sides and exchanges physical operators for a real H.
What it says, in standard words: the centred Heisenberg generator is minus the
centred Schrödinger generator dressed by ℒ_U. That is the shape of anti-pseudo-Hermiticity,
Q₋, one of the symmetries the Sá-Ribeiro-Prosen classification is built from,
which [KMS and detailed balance](../docs/KMS_DETAILED_BALANCE.md) lists by name
and which acts on the trace-shifted generator L_c. Whether a given U meets that
classification's own conditions on its symmetries is not checked here.

Mode by mode, a term e^{μt} on the left is carried on the right by the partner
mode at −μ. The two hands are reversed by two different maps, the split that
Π as time reversal draws for Π and this page draws above for the spectrum.
ℒ_U alone carries a mode of L_c at μ to an eigenoperator of L_c† at −μ, the dual
of the mode at −μ̄: the fade mirrored about σ, the turning kept. Hermiticity
preservation alone sends μ to μ̄: the turning reversed. The trace form uses
both. Whether its factor e^{μt} is read (a reading, either way) as the partner's two
hands reversed at the same t, or as the partner's exponent run at −t (the block fold's "the
partner world running backward"), is one number read two ways; nothing in the
mathematics chooses. [On two times](../reflections/ON_TWO_TIMES.md) holds both
in one paragraph, both partners evolving forward in γ₀-time and Π running the
centred flow backward on the mirrored state, and reads the look-back as
structural, not temporal.

For U = X^⊗N and Z dephasing the relation can also be written with F119's
transpose Lᵀ (the transpose of the generator's matrix in the matrix-unit basis)
and with F1's Π_Z(ρ) = ρᵀ·X^⊗N, the inverse of the Π above:
Π_Z Lᵀ Π_Z⁻¹ = −L − 2σ. This is the relation rewritten, through Lᵀ = conj(L†),
conjugation by the state transpose turning ℒ_U into R, and hermiticity
preservation once more. It
holds where neither of its two halves does, the Π_Z similarity and Lᵀ = L, as
the page already found for the similarity alone. Since Lᵀ has the spectrum of L,
it is for U = X^⊗N a similarity route to the palindrome, checked exactly below. The largest entry of each
residual, Z dephasing on every site, the same at N = 3 and N = 4 (the couplings
and term strengths are the gate's):

| system | Π_Z L Π_Z⁻¹ + L + 2σ | Lᵀ − L | Π_Z Lᵀ Π_Z⁻¹ + L + 2σ |
|---|---:|---:|---:|
| Heisenberg chain, uneven J and γ | 0 | 0 | 0 |
| the chain + a Y₀Z₁ term | 1 | 1 | 0 |
| the chain + an X₀Y₁Z₂ term | 2 | 2 | 0 |
| the chain + a Dzyaloshinskii-Moriya term XY − YX on every bond | 4 | 4 | 4 |
| the chain + a Z field on one site | 1.5 | 0 | 1.5 |
| a ring with a quarter flux | 4 | 4 | 4 |
| the chain, with σ⁻ added on one site | 0.375 | 0.25 | 0.375 |

The Z field keeps F119's transpose and breaks the rest; the other three controls
break all three.

**A conserved pairing between two forward runs.** Run two states forward with
the same t. Their U-weighted product

    ⟨ρ_B, ρ_A⟩_U = ⟨ρ_B, U ρ_A⟩ = Tr(ρ_B(t)† · U · ρ_A(t))

is constant under L_c and decays as e^{−2σt} under L; its derivative is
⟨ρ_B, (L_c† + ℒ_U L_c ℒ_U⁻¹)·Uρ_A⟩ = 0. This is the standard conservation law of
a pseudo-Hermitian generator, here iL_c with ℒ_U as its metric (a Hermitian
metric when U is Hermitian, and then always indefinite, since such a U
anticommutes with a jump and is traceless; the one-line calculation holds for
any admissible U). Under L both runs are physical (the product then decays as
e^{−2σt}), and the product is not a single-copy
expectation value but a two-copy one: for a unitary U and Hermitian states,
Tr(ρ_B U ρ_A) = Tr[(U ⊗ 1)·SWAP·(ρ_B ⊗ ρ_A)], the expectation of a unitary on
the two runs side by side (the standard swap trick), which Hadamard tests read. As a meter of σ it adds
nothing to the linear sibling below, which carries the same e^{−2σt} in one
copy; its content is structural, the pairing of modes it imposes. For
eigenmodes r_j, r_k of L_c, (μ_k + μ̄_j)·⟨r_j, r_k⟩_U = 0: the product pairs a
mode only with the mode at −μ̄, its partner across the centre turning the same
way, and a mode off the centre is neutral, ⟨r, r⟩_U = 0. For a Hermitian U with
U² = 1, X^⊗N among them, the product is an indefinite Hermitian form, the Krein
form of iL_c, and with P± = (1 ± U)/2 acting on the ket,

    ⟨r, r⟩_U = ‖P₊r‖² − ‖P₋r‖²,   lean(r) = ⟨r, r⟩_U / ‖r‖² ∈ [−1, 1],

the weight on the U-even side minus the weight on the U-odd side, the shape of
Lemma 1's normalized κ(v) with ℒ_U in place of τ. So, for a diagonalizable L_c,
⟨ρ, ρ⟩_U is a sum of terms that pair a mode slower than σ with one faster, and
of terms within the centre. The reflection page reads the two sides as past and
future, and reads the purity's cross term between them with the positive form;
whether the indefinite one bears on its first open item (what keeps the
reflection-even part of (2, 2) from a clean cut at N ≥ 4) is not checked here.

**Two forms, one for each hand: a reading.** What is proved: τ (Lemma 1, proved
there and not measured here: a real, number-conserving H, the diagonal blocks)
forces neutrality on every mode that turns, ω ≠ 0. ⟨·, ·⟩_U (U commuting with H;
the jumps Hermitian, squaring to 1, each anticommuting with U) forces it on
every mode off the centre, α ≠ σ, and pairs a block with its image under ℒ_U,
which for U = X^⊗N sends ket popcount p to N − p. Both give signs on one block
only where both exist and meet: for U = X^⊗N and a real H, the block
(N/2, N/2) at even N. With a complex, number-conserving H that commutes with
X^⊗N, such as the chain plus i(σ⁺₀σ⁻₁ − σ⁻₀σ⁺₁)Z₂ (not gated), the U form
exists and τ does not. On (N/2, N/2) with a real H, a turning mode on the centre is τ-neutral and can carry
a U sign; a non-turning mode off the centre can carry a τ sign and is U-neutral;
a turning mode off the centre is neutral for both; a non-turning mode on the
centre can carry both. The gate below evidences the U side, on (2, 2) at N = 4. The two are not
symmetric: hermiticity preservation reverses the turning for every L, and τ only
adds signs to it, while the fade is mirrored only where an admissible U exists.
The reading, on that block: measure the fade from the palindrome's centre σ
rather than from zero, and each hand of the clock has its own form, τ for the
turning (whose mirror, ω ↦ −ω, is hermiticity's) and ℒ_U for the fade (form and
mirror at once); a hand is off its line only on modes its form leaves neutral.

**How modes leave the centre.** On the block (2, 2) at N = 4, with
U = X^⊗N, H = J(XX + YY) or J(XX + YY + ZZ) on each bond, J = 1, and the same γ
on every site (σ = 4γ). Physically, a mode on the centre fades at exactly σ at
each γ for as long as it stays there. At a collision, leaving is the analogue of the
unbroken-to-broken transition of the standard passive PT picture, two rates
splitting to σ ∓ δ. On the XY chains most departures, counted at γ = 10⁻⁶, are
consistent with the standard thresholdless breaking of a degenerate eigenvalue of
mixed sign; on the
Heisenberg chain 4 modes leave at once and 4 more at the first collision. Later
collisions, and returns to the centre at strong dephasing, are not tracked
here. On
a degenerate eigenvalue only the inertia (n₊, n₋) of the form on the whole
eigenspace is invariant, and by Krein's rule a small perturbation releases at
most 2·min(n₊, n₋) of its modes; a definite eigenspace releases none. At γ = 0
every mode of the block is on the centre, and at γ = 10⁻⁶ the first-order
splitting of degenerate frequencies (first order in γ by degenerate perturbation
theory; the gate counts at this one γ) has moved 22 modes
off it on the uniform XY chain (Krein's bound from the γ = 0 eigenspaces: 24), 4 on the Heisenberg
chain (bound 6) and 22 on the uneven XY chain (bound 24). Among the modes that
stay, the first collision detected takes two of opposite lean off the centre at
+ω, with their complex conjugates at −ω, and their lean closes as
√(1 − γ/γ*), consistent with the generic square root at a second-order
exceptional point that the refinement of γ* assumes (the last column averages
the two):

| chain at N = 4, J = 1 unless stated | γ* (from the squared splitting) | the two that leave at +ω, just before (sign of the lean) | lean² / (1 − γ/γ*) → |
|---|---:|---|---:|
| XY, uniform | 0.29961298 | ω = 4.1546 (−), 4.1548 (+) | 1.7726 |
| Heisenberg, uniform | 0.06358514 | ω = 1.2642 (−), 1.2643 (+) | 2.0640 |
| XY, J = (1, 0.6, 1.4) | 1.56493622 | ω = 3.6442 (−), 3.6455 (+) | 2.7808 |

The collision rule is Krein's (Krein, Gelfand and Lidskii), and the square-root
closure is the generic behaviour of the signature at a second-order exceptional
point; the repository reads both already, with other metrics, in the β-exotic
genericity study, in the fragile bridge and in the PT-symmetry analysis (the
signatures ±0.0158 at δ = −10⁻⁴, shrinking to zero at the collision). Whether
ℒ_U, G, τ and the bipartite metric give the same signs where more than one
exists is not checked here.

**The world and its one-sided reading.** For U = X^⊗N, ℒ_U ρ is the reading
X^⊗N ρ that [the bridged lattice](../compute/MirrorWorld/Lattice.cs) runs beside
the world ρ, and ⟨ρ, ρ⟩_U = Tr(X^⊗N ρ²) is their overlap:

    Tr(ρ(t) · X^⊗N · ρ(t)) = e^{−2σt} · Tr(ρ(0) · X^⊗N · ρ(0))

for every state, every coupling with [X^⊗N, H] = 0 and every profile of local Z
dephasing (independent sites; under a shared field the price is no longer one
number, [the mirror under a shared field](THE_MIRROR_UNDER_A_SHARED_FIELD.md)).
It has a linear sibling the page already owns: U is the far-end eigenmode,
L(U) = −2σU and likewise L†(U) = −2σU, so Tr(X^⊗N ρ(t)) decays as e^{−2σt} for
every state too. On the cat pair cos θ·|0…0⟩ + sin θ·|1…1⟩ both are
2·cos θ·sin θ·e^{−2σt} (the textbook decay of a GHZ coherence at 2σ, the sum of
the sites' coherence rates 2γ_l), twice the spook term cos θ·sin θ·e^{−2σt} of
[the lattice opening law](LATTICE_OPENING_LAW.md). For every pure initial state
the two coincide at all t: they are equal at t = 0 and obey the same law. The
overlap differs from the one-copy parity ⟨X^⊗N⟩ only for a mixed initial state,
and then only through its value at t = 0, a static two-copy property of the
prepared state. That law's opening, an entry-wise maximum, is J-blind because the
cat sector is H-dead; the overlap is J-blind for every state. Its
Hilbert-Schmidt companion is ‖ρ − X^⊗N ρ‖² = 2‖ρ‖²·(1 − lean(ρ)): zero at lean
+1, largest at lean −1 (the anti-cat (|0…0⟩ − |1…1⟩)/√2, for one). The
lattice's own relation between its two one-sided readings, X^⊗N ρ and
ρ·X^⊗N, is the adjoint (`Restless.cs`, at `SeedRaw`), not a form.

**Swept for this subsection.** `docs/ANALYTICAL_FORMULAS.md`: F119 holds the
forward reciprocity and names the weight-one Krein form among its neighbours;
F1, F138 and F158 hold the one-sided relation, F109 its right-hand twin
R·L·R⁻¹ = −L† − 2σ; F50 the Krein count, F118 the factorization Π = R·D.
`docs/proofs/`: Lemma 1 (τ); the two-end count and MIRROR_SYMMETRY_PROOF the
operator relation; the Π factorization (R); the fold-lattice lemma of
PROOF_CODIM1_BY_ADDITIVITY §7 (the block fold); the antilinear triangle §7 leaves
open whether θ∘L∘θ = L\* fills a known symmetry slot. `experiments/`: the gamma
fold (the trajectory face), Π as time reversal (the centred fold, and the split
of the two hands for Π), the PT-symmetry analysis (Π∘† and the commuting S),
the β-exotic genericity study (the Krein collision and its square root), the
lattice opening law (the cat case), the Born-rule shadow and the reflection
page (the slow and fast sides of σ), the mirror under a shared field (the same
one-sided computation with the price an operator). `hypotheses/`: the fragile
bridge (G and R); [the pattern recognizes itself](../hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md)
compares expm(J·t)·Q with e^{−2st}·Q·expm(−J·t) for a neural Jacobian, a fold
with the envelope and no exchange of ends. `reflections/`: nothing on
indefinite forms; [on two times](../reflections/ON_TWO_TIMES.md) on Π as a
structural look-back rather than time reversal. `docs/`: KMS and detailed balance (Q₋ among the SRP classes,
L† against Π), and [the shifted order-4 chiral symmetry](../docs/outbound/SHIFTED_ORDER4_CHIRAL_SYMMETRY.md),
which places Π against the same SRP set. The typed layer:
`PalindromeTwoEndCountClaim` carries the relation, `AntilinearTriangleClaim`
F119, and `ExceptionalCouplingSetClaim` and `ExceptionalCouplingWitness` in
Diagnostics the inertia identity and Krein debt built on Lemma 1;
`LocalGlobalEpLink` in Core types the fragile bridge's collision;
`compute/MirrorWorld/` carries `Mirror`'s trajectory fold (F89's antiunitary
cross-block fold, in Core's `F89BranchLocusPalindromeClaim`, is a different
map), GammaFold's
cross-dock, `MirrorGroup`'s R and the lattice. The OpenArcs registry holds the
undressed one-sided fold F D F = −D − 2σ on the dissipator alone (the arc on
the self-folded (1, 3)) and other Krein signs; `fw.Confirmations`,
`simulations/framework/`, `compute/RCPsiSquared.Cli/`, the glossary and the
caught errors returned nothing on the trace form, the U-weighted product, or ℒ_U as a
Krein metric. Checked for adjacency: the gamma fold's cross-dock, the block
fold, Π as time reversal and Lemma 1 are about the same object as this subsection and cite
neither each other nor this page on it; R is ℒ_U's right-hand twin, as F109
already uses it; the product
splits along the reflection page's cut; the lattice opening law's spook is the
product and its linear sibling on one state.

Gate: [`simulations/pairing_condition_flow.py`](../simulations/pairing_condition_flow.py)
→ [`pairing_condition_flow.txt`](../simulations/results/pairing_condition_flow.txt).
The relation ℒ_U L_c ℒ_U⁻¹ = −L_c† is compared exactly at N = 3 and 4: residual
0.0, because ℒ_U for these U permutes matrix units up to sign (for X^⊗N·2^{N̂},
up to powers of two) and the only
floating sums on the diagonal, rates and energy differences, are dyadic here or (the named row,
U = SWAP₀₂·Z₀Z₁Z₂ at γ = 1/20) a single jump; a non-dyadic profile of four rates
leaves 2·10⁻¹⁶ on the diagonal, a rounding of the rate sums, printed and not
gated. An admissible U that is neither unitary nor Hermitian, X^⊗N·2^{N̂},
chosen so that its inverse is dyadic too, holds exactly in both
orders of evaluation, and its Π² = (ρ ↦ UρUᵀ) commutes with L exactly while
differing from Ad_U.
The trace form (N = 3, 4), the U-weighted product and the overlap (N = 4) are
consequences of that exact relation, which carries the evidence; the gate also
runs them through matrix exponentials as a consistency check, beside the forward-only rewrite, which holds by Ad_U alone (for
any U commuting with H with UAU⁻¹ = ±A for every jump, so it tests nothing
one-sided; only the backward form does), compared against a
heuristic rounding scale rather than an error law, a holding row passing within
ten times it: for the unitary U the holding rows sit two to four decades under
it and the controls more than ten decades above it; for X^⊗N·2^{N̂²} the scale
carries the factor ‖U‖·‖U⁻¹‖ and the holding rows sit six to nine decades
under it, its controls eight to eleven above. The neutral
modes sit within ten times the heuristic scale ε·dim·cond(V), with a Z-field
control that leans 0.90; the departures are compared against Krein's
bound; the colliding pair is tracked across γ*, located first on a 1.3× grid,
then by bisection, then from the vanishing of the squared splitting. That
refinement assumes the square-root splitting of a second-order exceptional
point, so the lean's square root is checked against it, two observables of one
structure, not measured independently (with the bisection's γ* instead, the
fourth point already fails the convergence check on the two uniform chains): the series lean²/(1 − γ/γ*), at four
points over three decades, converges geometrically, each step at most a third
of the last (which rules out exponents below ½ and slightly above it), to a
limit that dominates its correction and is at least a tenth of the series'
first value (which rules out the larger exponents). The trace form and the forward-only rewrite are also run with
U = X^⊗N·2^{N̂²}, whose square is not a multiple of the identity, so that U and
U⁻¹ cannot be exchanged unseen. The two-copy formula is checked numerically, its wrong order of
the copies as the control. The table's zeros are compared exactly.

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
| the field has one common axis, orthogonal to every dephasing axis | U must commute with the field. For U a single Pauli string that puts the field's letter (or I) on the field's site, and a bond of two or three letters carries one letter across the component, so U is the global π-rotation about the field axis, which anticommutes with the dephasing letters because they are orthogonal to it; one-letter bonds of one shared letter let the letter change from site to site ([the palindrome as a colouring](THE_PALINDROME_AS_A_COLOURING.md)) |
| the two-term proviso | derived in [the palindrome as a colouring](THE_PALINDROME_AS_A_COLOURING.md) for one bond letter per component: a single-letter bond only puts its two ends in one of two letter classes, so with one bond letter a single Pauli string U can change its letter from site to site (bonds of different single letters can leave a site no letter), where a bond of two or three letters forces one letter per component. The gate scores all seven bond sets |

The demonstration above is one of the six rows among the 78 whose field is
anti-invariant under the end-site swap. All 78 carry the same kind of U,
SWAP₀₂·(n·σ)⊗(n·σ)⊗(n·σ): the site reversal composed with the π rotation about an
axis n that carries the one end field onto the other, fixes the middle field and
is orthogonal to the middle jump. Here n = Z; on the other 72 the axis is unique up to
sign, a letter or the bisector of two letters ([F138](../docs/ANALYTICAL_FORMULAS.md#f138), gated row by
row in [`EndCountTests`](../compute/MirrorWorld.Tests/EndCountTests.cs)).

---

## Evidence

Gate suite: [`simulations/f138_pairing_condition.py`](../simulations/f138_pairing_condition.py)
→ [`f138_pairing_condition.txt`](../simulations/results/f138_pairing_condition.txt).
It runs on the exact GF(p) kernel
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
none anywhere else. **These rows are confirmed nonempty at one prime, not
at three**: `admissible_multi` returns as soon as a prime yields no invertible
draw, so the later primes are never consulted on them. **None of them hides an
invertible element**: every such row breaks, so by the theorem `dim 𝒲 < dim 𝒩`
there (inferred, not computed row by row), and
[the two-end proof](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md) §(d) shows
that a strict inequality is a proof that no invertible element exists, so what
the draws sampled the count decides.

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

The subsection on the flow has its own gate,
[`pairing_condition_flow.py`](../simulations/pairing_condition_flow.py), named
there.

---

## Scope, and what would falsify it

Tested: N = 3 and N = 4; chain, ring, triangle, and a disconnected bond plus
isolated site; all seven bond-letter sets; uniform and non-uniform rational J and
γ; one, two and three distinct dephasing axes; single-site Pauli-letter jumps (off-axis directions and Pauli-string jumps are scored below).

**Not tested.** An anisotropy axis: every bond carries weight 1 per letter and
there is no Δ. N beyond 4 for the criterion as scored on this page.

**Off-axis dephasing is proved and tested, and it was never a new axis.** The theorem asks of a jump operator only that it be Hermitian and
square to the identity, which any `n̂·σ⃗` at a unit direction satisfies, and so
does any multi-site Pauli string. Both are scored in
[`f138_rank_criterion.py`](../simulations/f138_rank_criterion.py), 1,692 rows
with rational directions and string jumps, FP = 0 and FN = 0, and the sibling
arc `f138_clause_two_sweep`'s question about clause 2's word *orthogonal* is
answered there: orthogonality of directions is not what the criterion asks for.

**And one axis where the criterion is not merely untested but mis-centred.**
Jump operators with A² ≠ 1 break Σ c†c = σ·1, and the palindrome does not simply
die there: **F137** records that T1 alone keeps it, about **−σ/2** rather than −σ.
So a criterion phrased "pairs about −σ" is answering the wrong question in that
regime, and the premise is doing more work than a footnote. F82 to F84 own the
separate matter of Π failing there. Z dephasing with correlated site fields sits
on the same axis, with Hermitian jumps whose Σ c†c is not a number: the same one-sided
computation with U = X^⊗N, which anticommutes with every jump built linearly from the Z_k
(H commuting with it), premise dropped, gives −L† − {Σ c†c, ·}, so the shift 2σ becomes
the operator {Σ c†c, ·}. When the off-diagonal Kossakowski entries (the fields'
covariances) of a positive semidefinite Kossakowski matrix over the Z_k sum to a positive
number and H keeps the all-0 and all-1 states as
eigenstates, no single centre serves the whole spectrum.
Under a field shared alike by every pair, with a popcount-conserving H that
commutes with X^⊗N, the block (p, q) (ket popcount p, bra popcount q) pairs with
(q, N − p): ℒ_U sends it to (N − p, q), and ρ ↦ ρ† turns that into (q, N − p), each
pair about a centre of its own
([The Mirror under a Shared Field](THE_MIRROR_UNDER_A_SHARED_FIELD.md)).

**Necessity is derived, not only measured.** Sufficiency is a short
calculation; necessity is proved in
[PROOF_PALINDROME_TWO_END_COUNT.md](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md)
§(d) and §(e), and the route is not this page's: not through the operator at all, but
through the two ends of the axis. L is unital, so 0 is an eigenvalue; the
palindrome carries its multiplicity to −2σ; both eigenvalues are semisimple; so
the two kernels have equal dimension, and equal dimension is what forces an
invertible U. The census illustrates the theorem; it is not its evidence.

What would falsify it: one row, at any setting, where an invertible U exists and
the characteristic polynomial does not pair.

---

## What this leaves open

**On the flow.** Whether the U-weighted form, the chain reflection's, τ and the
bipartite metric give the same signs where more than one exists, and whether the indefinite form
bears on the reflection page's cut at N ≥ 4. The flow subsection's gate is
Python; the typed layer states the relation in `PalindromeTwoEndCountClaim`'s
comment and checks it nowhere entry by entry, so the gate is a witness waiting
to be ported.

**The typed layer.** The witness is `PalindromeTwoEndCountWitness`, live at
`inspect --root twoend` beside F158.

**The arc `f138_converse_failures`** says the named candidate space is empty and
that existence is free where diagonalizability holds, and grades that evidence as
*"CERTIFIED ONLY: the named row measures 30 of 64 … diagonalizability rests on
g(L) == 0 mod p at three primes, which is not a decision"*. Both stay true and
neither is needed: the operator exists, it is structured, and it is found by a
linear solve scored in both directions rather than inferred from cospectrality.
The arc's finding (5) files V = R·U as not the answer, for a reason of kind; V is
Ad_U, which is Π² for that row's orthogonal U.

**F138 exhibits operators on both kinds of exception**: the site reversal composed with a π rotation, SWAP₀₂·(n·σ)^⊗3, on all 78 of its equal-magnitude exceptions, and on a named family of its clause-1 exceptions a SUM of two Pauli strings,
neither a single string nor a signed site permutation. The bare site reflection,
a signed site permutation, is not a reflector on those rows, which says nothing
against the existence of one. The row keeps the separation the arc prescribed,
EXISTS from EXPLAINS.
F138 stays an implication whose converse is false, measured in its letter-field sweeps to fail at fewer
bond letters and at coincident field magnitudes; the criterion is not F138's two
clauses.

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
