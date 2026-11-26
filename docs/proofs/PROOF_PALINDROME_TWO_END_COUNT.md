# PROOF: The Palindrome Is a Count of the Two Ends

**Registry:** [F158](../ANALYTICAL_FORMULAS.md).
**Status:** Tier 1 derived. The theorem holds for any Hermitian H, at least one
jump operator, all of them Hermitian and squaring to the identity, any strictly
positive rate profile, any topology, any N, and any finite dimension. What the
CRITERION decides is the eigenvalue multiset; sufficiency in fact delivers more,
a full similarity `L ~ −L − 2σ`, so where the criterion holds the entire Jordan
structure is reflection-symmetric (§(f8)). Two fences are load bearing rather
than decorative. The first is Hermitian jumps with `A² = 1`: the hermiticity makes the
dissipator self-adjoint, which sufficiency spends, and `A² = 1` then makes each jump unitary,
which the three lemmas and necessity spend (F137 recentres the palindrome under T1, which is not
unitary, where it survives; unitary jumps whose dissipator is self-adjoint keep the theorem with the same proof, F172).
The second is `γ_l > 0` with at least one jump present (a zero rate drops a condition and
can move the verdict; with no jump at all the equivalence holds trivially, and
§(f5)'s odd-d corollary, Lemma 3's 𝒩 ∩ 𝒲 = 0 and §(f9) fail). Even TOTAL dimension d is not a fence but a
consequence: with a jump present no palindrome exists at odd d (§(f5)), and the
theorem holds there with both sides false. It concerns the total dimension, not the
local one: on a qutrit tensored with a qubit, d = 6, the row H = 0, A = 1₃ ⊗ Z pairs,
with U = 1₃ ⊗ X. The rate fence is gated (10b) and the odd-d inequality behind the
corollary (10c); `A² = 1` is the boundary of the stated class, and only F172's gate crosses it, on the unitary side. The gates
are of three kinds and only the first is exact in the strict sense: exact over
ℚ(i) with Fraction arithmetic on named rows, exact over GF(p) at scale where a
nullity can read too large, and a float route whose thresholds are gated on a
measured separation.
**Date:** 2026-08-28
**Authors:** Thomas Wicht, Claude (Opus 5)
**Script:** [`simulations/f138_rank_criterion.py`](../../simulations/f138_rank_criterion.py)
→ [`f138_rank_criterion.txt`](../../simulations/results/f138_rank_criterion.txt)
**Builds on:**
- [The pairing condition](../../experiments/THE_PAIRING_CONDITION.md), whose
  criterion this sharpens and whose sufficiency calculation it consumes
  unchanged; this file supplies the derivation of necessity.
- [F1](../ANALYTICAL_FORMULAS.md) for the palindrome, and
  [F138](../ANALYTICAL_FORMULAS.md#f138) for the law whose converse is false.
- [PROOF_CODIM1_BY_ADDITIVITY](PROOF_CODIM1_BY_ADDITIVITY.md) §6, whose
  window-edge lemma owns Lemma 2 in a stronger form, and which the arc ledger
  instructs be cited rather than re-derived.
- [PROOF_F103_F87_Z2_CUBED_REFINEMENT](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md)
  §7.5 for the forcing step in a narrower setting, and §7.12 for the three rows
  gated here.
- [Degeneracy Palindrome](../../experiments/DEGENERACY_PALINDROME.md) Result 2
  for both counts of the canonical case.

<!-- STAGE 0 SWEEP RECORD (for us, not for the reader: hidden on GitHub, read by every session, agent and hook)

What the repo already held, store by store

The sweep ran through two agents and a hand pass, and three adversarial
reviewers checked it against the stores at source.

- **[`docs/ANALYTICAL_FORMULAS.md`](../ANALYTICAL_FORMULAS.md).** F1 holds
  `Π·L·Π⁻¹ = −L − 2Σγ·I` with Π unitary and Π² = X^⊗N. **F4 holds one of the two
  counts**, in three separate statements that must not be fused: its Σγ = 0
  bullet says the kernel is the full commutant of H, while the Wedderburn form
  `dim ker = Σ mᵢ²` appears only in the one-seat, γ > 0, single-popcount-sector
  bullet for the algebra ⟨H_w, n_seat⟩. F138 holds the law under test, with two
  operators exhibited on its exceptions. F137 fences the premise rather than answering
  it: T1 jumps keep a palindrome for the Hamiltonians it names and recentre it, and F137 states that centre
  **exactly**, as trace(L)/dim. Nothing anywhere on the two counts as a
  **criterion**.
- **[`docs/proofs/`](.), and the sharpest hit first.**
  [PROOF_CODIM1_BY_ADDITIVITY](PROOF_CODIM1_BY_ADDITIVITY.md) §6's **window-edge
  lemma** owns Lemma 2 below in a stronger form: *"If an eigenvalue of L sits on
  an edge of the block's rate window … then its eigenvector v is a joint
  eigenvector of A and B, and λ is semisimple."* That is the same
  numerical-range argument at **every** window edge rather than at the two
  extremes, and it also delivers most of Lemma 1's reverse inclusion. The arc
  `site_resolved_vacuum_block` carries a standing instruction about the
  rate-window paragraph this lemma sits in, *"Cite it rather than re-deriving it a third time"*, so Lemma 2 is
  written out below for readability and is **not** claimed as new. That arc
  also carries the γ-profile
  bracket §(b)'s corollary sits inside, and a minting caution that applies here
  directly: *"Do not mint it; the F156 withdrawal of the same morning is the
  precedent."*
  [PROOF_F103_F87_Z2_CUBED_REFINEMENT](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md)
  §7.5 owns the forcing step: *"the equation Σ_l Z_l A Z_l = −N·A forces, per
  site, Z_l A Z_l = −A"*, with A drawn from H's commutant, which is this file's
  pair of conditions exactly. Its scope is the whole difference: §7.5 runs on
  the first-order ω = 0 block of the dissipator restricted to the commutant, for
  F87's windowed diagonal cell; Lemma 1 runs on the full L at finite γ. §7.12
  then bounds its own reach and **does not leave the restoring operator
  unexhibited**: the
  sentence after the one worth quoting reads *"that operator is no longer a
  mystery: it is the hidden-Q routing, a per-site Q from the P1/P4 families,
  which `TwoTermPalindromeRouting` classifies bit-exactly for 2-term pairs"*,
  and those C# files are on disk. Its three rows are gated here as an
  **agreement test against a route the repo already owns**.
  [PROOF_STAR_OPTICAL_CONFOCAL_SATURATION](PROOF_STAR_OPTICAL_CONFOCAL_SATURATION.md)
  (registered as **F147**) owns the same equality-case move on the imaginary
  axis, with the rigidity written out. [PROOF_ABSORPTION_THEOREM](PROOF_ABSORPTION_THEOREM.md)
  holds `Re(λ) = −2 Σ_l γ_l·light_l(v)` and reads one of its extremes in its own
  opening
  (*"eigenmodes of pure {X, Y}^⊗N content die fastest (rate 2Nγ)"*). It also records
  `(X_k + Z_k)/√2`, a member of this file's class, as holding only in the frame
  turned 45° about y, so nothing here leans on it.
  [MIRROR_SYMMETRY_PROOF](MIRROR_SYMMETRY_PROOF.md) owns the first link, *"the
  palindrome then forces its partner −2Σγᵢ to be an eigenvalue too"*, from trace
  preservation. Its Scope states F138's clauses with *when*, not *exactly when*,
  and the converse failure and its counts beside them.
  [PROOF_F111_HARD_CELL_PURE_D_TEMPLATE](PROOF_F111_HARD_CELL_PURE_D_TEMPLATE.md)
  names the far end **the anti-steady eigenvalue** and pairs a mechanism with a
  **measured** witness (there is no live C# witness for it); the why is answered
  there through Theorem 2 of the complement connection, the route of §(f10), and
  §(f4) translates the same question to the letters.
  [PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA](PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md)
  computes a commutant by hand and reads it off in **Wedderburn** form; the
  words Schur and double commutant do not appear in it.
  [PROOF_ASYMPTOTIC_SECTOR_PROJECTION](PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md)
  Step 2 owns the fixed-point algebra as decoherence-free subalgebra ∩
  commutant, which is the frame Lemma 1's kernel end lands in.
- **[`experiments/`](../../experiments/), and the join is narrow.** [DEGENERACY_PALINDROME](../../experiments/DEGENERACY_PALINDROME.md)
  Result 2 carries **both** counts, on **one page**, and states the bijection
  outright: *"At Re = 0: these are the N + 1 conserved quantities … At Re = −Nγ:
  these are the N + 1 XOR sector modes, the fastest-decaying … conserved
  quantities at Re = 0 bijectively to the XOR modes at Re = −Nγ."* So the repo
  held the two numbers together, on one page, and read the bijection as a
  consequence of Π. What is new is not the pairing of the two
  counts but that their **equality decides the palindrome in general**, off the
  canonical case where Π is available.
  **And the two γ books must not be lumped, which is the repo's own documented
  trap** (fenced in [GLOSSARY](../GLOSSARY.md) and
  [CAUGHT_ERRORS](../CAUGHT_ERRORS.md)):
  DEGENERACY_PALINDROME writes the far end as −Nγ at γ = 0.1 per site, while
  [XOR_SPACE](../../experiments/XOR_SPACE.md) writes it as −2Σγ in the canonical
  Lindblad book. Those are the same point (γ_there = 2γ_canonical; at N = 3 and
  γ = 0.05 the GLOSSARY row gives 2Σγ = 0.30), and this file works throughout in the canonical
  book, in which the far end is −2σ.
  [SYMMETRY_CENSUS](../../experiments/SYMMETRY_CENSUS.md) carries both operator
  conditions months old and reads them as invariance, correctly for the two-sided
  conjugation; the pairing-condition page adds that one-sided multiplication by
  the same X^⊗N reflects.
- **[`docs/GLOSSARY.md`](../GLOSSARY.md)** carries the far-end count measured
  three times and named only as modes: *"NONE at all under a generic Hermitian
  H, whose maximum rate falls short of 2Σγ entirely (measured at N = 3,
  γ = 0.05: 4, 8 and 0 modes …)"*. Those three integers count modes on the
  line Re λ = −2σ. A count of modes on that line can include oscillating ones,
  which §(f1) keeps out of the kernel, so such a count is `dim ker(L + 2σ)` only
  where no mode on the line oscillates; Lemma 2 says only that at the point −2σ
  itself the algebraic and geometric multiplicities agree.
  No headword for the criterion.
- **[The OpenArcs registry](../../compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs),
  and TWO arcs answered.**
  `site_resolved_vacuum_block` (opened 2026-08-02, still Open) holds the
  citation instruction quoted above and the numerical-range bracket
  `Re λ ∈ [−2σ, 0]` as a bound; `the_gate_that_does_not_gate` (opened
  2026-08-22) holds the commutant reading of a kernel dimension, which this file
  reaches through the markdown layer instead. And
  arc `f138_converse_failures` asked for an S
  with `S L S⁻¹ = −L − 2σ` and recorded that the
  bare site reflection was refuted from below. Its caution needs both of its
  conjuncts: an unstructured S exists wherever L is diagonalizable **and the
  characteristic polynomials already agree**, i.e. wherever the palindrome
  already holds, so free existence is not an alternative explanation of anything.
  §(d) is why the answer here is not free in any case: it is a dimension count,
  it is often unequal, and where it is equal the operator is forced.
- **[`docs/CAUGHT_ERRORS.md`](../CAUGHT_ERRORS.md)** returned *"a search built
  from the known instances measures the search, not the codebase"*, dated
  2026-08-06 and written about a grep sweep for retired code rather than about
  operator searches. It transfers by analogy, which is how the pairing-condition
  page uses it, and this file has no search in it at all.
- **[`compute/MirrorWorld/`](../../compute/MirrorWorld/).** −2σ appears
  mostly as a **shift constant**: GammaFold adds `−2σ·Id` per cell and wears
  `e^(−2σt)` as a veil; Mirror's price 2Nγ is a fold's affine cocycle between
  two different blocks. Three objects read the far end itself: EndCount counts
  both ends of this theorem modulo two primes and returns far elements, Stitch
  builds far elements X^⊗N·D off the flat sections, and Carrier holds one row's
  single far element; `Cat`'s `|0…0⟩⟨1…1|` at −2Nγ lives at H = 0.
- **Returned nothing, or less than it seems.** `recovered/` returned nothing.
  `reflections/` did NOT: `ON_TWO_TIMES.md` carries the far end (*"Near λ = 0 the
  memory is long; at λ = −2Σγ it is shortest"*), while this
  file's own companion page cites `ON_THE_SOFT_BREAK.md` twice.
  `fw.Confirmations` and the C# `ConfirmationsRegistry` carry nothing on the
  two-end count. And `simulations/framework/diagnostics/f158_end_count.py`
  counts both ends exactly over ℚ(i) on the Pauli strings, beside
  `far_element_shape.py` and `f158_odd_word.py`, while
  `diagnostics/d_zero.py::stationary_modes` returns a kernel
  dimension and pins N+1, which is §(f2)'s near-end count computed and running,
  by the float route F157's registry entry warns against. And
  `hypotheses/THE_OTHER_SIDE.md` Q6 is **not** open. It is marked *"ANSWERED 2026-06-01 (Klein routing): the mapping is now
  built and verified bit-exactly (Q·L·Q⁻¹ = −L−2Σγ·I to ‖·‖ ≤ 10⁻¹¹,
  N=3,4,5)"*, which is the repo constructively exhibiting a palindromizer for a
  family of palindromic cases.
-->

**And one page in exactly this territory:** [DEPOLARIZING_PALINDROME](../../experiments/DEPOLARIZING_PALINDROME.md)
§8 states an iff in exactly this territory, *"the Liouvillian palindrome holds under dephasing
noise if and only if the noise has at most two Pauli axes"*, and it is the
SAME-SITE reading, where F138 clause 1 and the pairing page answer the
across-sites one (three axes on three different sites have a common
anticommutant and the criterion carries no ceiling). Different objects, and
the same-site one is where an iff can bite. That law is also **conditional on
H**, as that page's abstract states (its §8 theorem is stated for Heisenberg
coupling on any graph): the
criterion says the U for two-axis noise is the global product of the MISSING
letter, so H may carry a field along that letter and no other. Under Z + X noise
on every site with U = ΠY, a Y field commutes with U and an X or a Z field does
not, since conjugation by Y flips both. The Y-field row is the one nobody would
guess, and it is a scope fence this theorem owes that one.

**What is new, and it is narrow.** Not the
forcing step (F103 §7.5, in a smaller setting). Not the equality-case technique
(F147 on the imaginary axis, PROOF_CODIM1's window-edge lemma on this one). Not
semisimplicity, which PROOF_CODIM1 owns in a stronger form. Not either count
(F4, and DEGENERACY_PALINDROME which holds both together). New here: that the
far end is a **kernel of the same kind** as the near end, with the two spaces
differing by one sign; that their two dimensions are equal **exactly** when the
palindromizer exists, which turns an existence question into a rank comparison;
and therefore the necessity direction, which the pairing-condition page could
only measure.

## Abstract

Write the generator with Hermitian, unitary jumps,

    L(ρ) = −i[H, ρ] + Σ_l γ_l (A_l ρ A_l − ρ),   A_l† = A_l,  A_l² = 1,
    γ_l > 0,   σ = Σ_l γ_l .

Its kernel is the commutant of the algebra generated by H and the jumps. Its
eigenspace at the opposite end of the palindrome axis, −2σ, is the same object
with one sign flipped: the operators that commute with H and ANTIcommute with
every jump. Both eigenvalues are semisimple, so both geometric dimensions are
also algebraic multiplicities, and a palindromic spectrum forces the two to
agree. Conversely, when they agree the anticommutant contains an invertible
element, and one-sided multiplication by it reflects the spectrum. Hence

> **The spectrum multiset of L is closed under λ ↦ −λ − 2σ if and only if
> `dim ker L = dim ker(L + 2σ)`.**

There is no operator to find and no subspace to sample: the criterion is two
nullities compared, each an exact rank. On the canonical Heisenberg chain under
Z-dephasing both counts are N+1, which is the repo's own founding pair of
numbers; what the theorem adds is that their equality is not a consequence of
the palindrome but the same statement.

## (a) Setting, and what each hypothesis is doing

Let ℋ = ℂ^d, H = H† on ℋ, and A_1 … A_m Hermitian with A_l² = 1 (so each A_l is
also unitary). Let γ_l > 0 and σ = Σ_l γ_l > 0. On the Hilbert-Schmidt space
(ℬ(ℋ), ⟨A, B⟩ = Tr(A†B)) define

    L(ρ) = −i[H, ρ] + Σ_l γ_l (A_l ρ A_l − ρ) .

This is a Lindbladian with c_l = √γ_l·A_l, since c_l† c_l = γ_l·1 makes the
anticommutator term −γ_l ρ. Two consequences are used and neither is optional:
L preserves hermiticity, and L is **unital**, L(1) = 0.

Write 𝒜 for the unital *-algebra generated by {H, A_1 … A_m}; it is a *-algebra
because every generator is Hermitian. Write

    𝒩 = {X : [H, X] = 0 and A_l X A_l = X for every l}    (the commutant 𝒜′)
    𝒲 = {W : [H, W] = 0 and A_l W A_l = −W for every l}   (the anticommutant)

Since A_l² = 1, `A_l X A_l = ±X` is the same statement as `A_l X = ±X A_l`.

In Buča and Prosen's terms ([arXiv:1203.0943](https://arxiv.org/abs/1203.0943)),
the unitaries in 𝒩 are exactly the strong symmetries of this generator, the
unitaries U with [U, H] = 0 and [U, A_l] = 0 for every jump. A unitary U in 𝒲
commutes with H and carries every jump onto its negative, U·A_l·U† = −A_l, so
the conjugation ρ ↦ UρU† keeps each dissipator and commutes with L: a weak
symmetry, and not a strong one, since [U, A_l] = 2·U·A_l ≠ 0. Not every
unitary weak symmetry lies in 𝒩 or 𝒲: one that keeps some jumps and reverses
others, or exchanges two jumps of equal rate, lies in neither. [The colouring page](../../experiments/THE_PALINDROME_AS_A_COLOURING.md)
reads a colouring, a Pauli string in 𝒲, this way.

**Two fences, both load bearing, and one consequence.** The rate fence is gated (10b);
Hermitian jumps with `A² = 1` are the boundary of the stated class, and the theorem rests on two hypotheses: the hermiticity makes the dissipator self-adjoint, which sufficiency spends, and `A² = 1` then makes each jump unitary, which the three lemmas and necessity spend; even total d follows (below).

**γ_l > 0.** Not a formality. A rate set to zero removes that site's condition
from both 𝒩 and 𝒲, and the verdict can move: on a ZZ bond with an X field on site
0 and X-dephasing on both sites, the spectrum is broken while both rates are on
and palindromic the moment site 0 stops being watched (gate 10b). It need not move: with
two jumps both equal to Z, switching one off drops no condition, since both
impose the same one, and the verdict stays whatever H is. So the
boundary of the positive orthant is where the verdict CAN change, and "the
palindrome does not depend on γ" in §(f3) means inside the open orthant.
Negative rates are outside the class altogether, which is worth saying because
[PROOF_R90_FROZEN_DIVISOR](PROOF_R90_FROZEN_DIVISOR.md) works a zero-mean
stratum where σ = 0 with jumps present; that is a different object. σ = 0 with
positive rates would need **no jumps at all**, which the class excludes; there 𝒩
and 𝒲 coincide, the two ends of the axis are one point, the spectrum of −i·ad_H
is symmetric about 0 for any Hermitian H, and the equivalence holds trivially.
Without a jump §(f5)'s odd-d corollary fails, and so do Lemma 3's 𝒩 ∩ 𝒲 = 0 and §(f9).

**A_l² = 1** is where the stated class ends. F137 records that under T1 jumps, which
are not unitary, the palindrome survives, for the Hamiltonians it names, about a different centre, which F137
states exactly as trace(L)/dim, so a criterion phrased about −2σ is answering a
different question there. The theorem reaches past that boundary on the unitary side: for unitary
jumps whose dissipator is self-adjoint, Σγ_l A_lρA_l† = Σγ_l A_l†ρA_l, it holds by this proof with
A_lρA_l† in place of A_lρA_l, and without the self-adjointness only its necessity is guaranteed,
and where the counts agree what remains is one spectral equality, L against L with every jump replaced by its adjoint
([the clock mirror](PROOF_QUDIT_CLOCK_MIRROR.md) §5, F172). What Hermitian jumps with `A² = 1` buy, positively, is broader than
the letters: any `n̂·σ⃗` at a unit direction qualifies, as does any Pauli STRING,
as does a full depolarizing site (X, Y and Z each satisfy it). All three
families are gated in §(g), and the last of them is F1's canonical break, which the criterion
gets right from outside the sample.

**Even total d.** Not assumed, but derived: §(f5) shows that with at least one jump
no palindrome is possible at odd total d, so the theorem is true there with both
sides false; gate 10c checks the inequality dim 𝒲 < dim 𝒩 behind it at d = 3 and 5.

## (b) Lemma 1: the two ends are the commutant and the anticommutant

> **Lemma 1.** `ker L = 𝒩` and `ker(L + 2σ) = 𝒲`.

**The easy inclusions are pure algebra and hold over any field.** If X ∈ 𝒩 then
L(X) = 0 + Σ γ_l(X − X) = 0. If W ∈ 𝒲 then
L(W) = 0 + Σ γ_l(−W − W) = −2σ·W.

**The reverse inclusions are equality cases, and they need ℂ.** Let
L(W) = −2σ·W. Take the Hilbert-Schmidt inner product with W:

    ⟨W, −i[H,W]⟩ + Σ_l γ_l ( ⟨W, A_l W A_l⟩ − ‖W‖² ) = −2σ‖W‖² .

The first term is purely imaginary: ⟨W, −i[H,W]⟩ = −i·Tr((WW† − W†W)H), and the
trace of a product of two Hermitian operators is real. Taking real parts and
using Σγ_l = σ,

    Σ_l γ_l · Re⟨W, A_l W A_l⟩  =  −σ‖W‖² .

Each A_l is unitary, so ‖A_l W A_l‖ = ‖W‖ and Cauchy-Schwarz gives
Re⟨W, A_l W A_l⟩ ≥ −‖W‖² term by term. The weights γ_l are strictly positive and
sum to σ, so a sum of terms each bounded below by −‖W‖² can equal −σ‖W‖² only if
**every** term attains its bound. Equality in Cauchy-Schwarz between two vectors
of equal norm forces

    A_l W A_l = −W    for every l ,

and substituting that back into L(W) = −2σW leaves −i[H, W] = 0. So W ∈ 𝒲.

The kernel end is the same four lines with one sign changed: L(X) = 0 gives
Σ γ_l Re⟨X, A_l X A_l⟩ = +σ‖X‖², every term is bounded ABOVE by ‖X‖², so every
term attains it, A_l X A_l = X, and then [H, X] = 0. ∎

The kernel half is the standard fixed-point statement for a unital semigroup and
the repo already works in that frame
([PROOF_ASYMPTOTIC_SECTOR_PROJECTION](PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md)
Step 2, F4). The far half is the one that was not written down, and it is the
same argument.

**One corollary, free from the same computation and used in §(f1):**

    Re⟨X, L X⟩ = Σ_l γ_l ( Re⟨X, A_l X A_l⟩ − ‖X‖² ) ∈ [−2σ‖X‖², 0] ,

so the numerical range of L lies in the strip −2σ ≤ Re λ ≤ 0, and both ends of
the palindrome axis are edges of it. That is the setting
[PROOF_CODIM1_BY_ADDITIVITY](PROOF_CODIM1_BY_ADDITIVITY.md) §6's window-edge
lemma already works in.

**A second corollary, for Pauli-string jumps.** When every A_l is a Pauli string
(up to sign), conjugation by A_l sends each Pauli string P to +P or −P, the minus
sign exactly when P and A_l anticommute. Expanding an operator in the string
basis, A_l X A_l = X for every l holds exactly when X is supported on the strings
that commute with every jump, and A_l W A_l = −W for every l exactly when W is
supported on the strings that anticommute with every jump. So

    𝒩 = ker ad_H on span{P : P commutes with every A_l} ,
    𝒲 = ker ad_H on span{P : P anticommutes with every A_l} .

Commutation between the jumps is not used. When H is written as a sum of
distinct strings T with coefficients h_T, [T, P] is 0 or 2·T·P, and for a fixed P
distinct terms give distinct strings T·P, so ad_H/2i on either span is a matrix whose
every entry is ±h_T: an integer matrix when the coefficients are integers, which
a rational H becomes after scaling by its common denominator, and scaling moves
neither kernel. The second span is empty when no string anticommutes with every
jump, and then 𝒲 = 0. This is the route the live witness `inspect --root
twoendstrings` (`PalindromeStringSpanWitness`) and
[`compute/MirrorWorld/EndCount.cs`](../../compute/MirrorWorld/EndCount.cs)
count by, with the ranks taken modulo two primes, and the route of
[`simulations/anticommuting_sum_census.py`](../../simulations/anticommuting_sum_census.py).

## (c) Lemma 2: both ends are semisimple

> **Lemma 2.** The eigenvalues 0 and −2σ of L carry no Jordan blocks, so at both
> ends the algebraic multiplicity equals the geometric one:
> `alg(0) = dim 𝒩` and `alg(−2σ) = dim 𝒲`.

**This is not new**, and the repo says where it lives: PROOF_CODIM1's window-edge
lemma gives semisimplicity at every edge of a rate window, and by §(b)'s
corollary 0 and −2σ are the two ends of the strip that bounds this generator's
spectrum, and they are edges of it exactly where they are eigenvalues, the only
case in which semisimplicity there says anything. It is written
out because the chain below consumes it, not because it is being claimed.

Let M = L + 2σ. By §(b)'s corollary L is **dissipative** (Re⟨X, LX⟩ ≤ 0) and M
is **accretive** (Re⟨X, MX⟩ ≥ 0). Let M W = V with M V = 0, i.e. W is a Jordan
vector of M above 0. For every real t,

    0 ≤ Re⟨W + tV, M(W + tV)⟩ = Re⟨W, MW⟩ + t‖V‖² ,

because ⟨V, MW⟩ = ⟨V, V⟩ = ‖V‖², ⟨W, MV⟩ = 0 and ⟨V, MV⟩ = 0. Letting t → −∞
forces ‖V‖² = 0, so V = 0 and ker M² = ker M. A Jordan chain of any length would
contain such a pair at its foot, so there is none of any length. The same
argument with −L in place of M gives ker L² = ker L. ∎

## (d) Lemma 3: an invertible element exists exactly when the counts agree

> **Lemma 3.** `dim 𝒲 ≤ dim 𝒩` always, and 𝒲 contains an invertible element
> if and only if `dim 𝒲 = dim 𝒩`.

**(⟹), and it is three lines.** Suppose U ∈ 𝒲 is invertible. If X ∈ 𝒩 then
XU ∈ 𝒲, since [H, XU] = 0 and A_l(XU)A_l = (A_lXA_l)(A_lUA_l) = X·(−U). If
W ∈ 𝒲 then WU⁻¹ ∈ 𝒩, since A_lU A_l = −U gives A_lU⁻¹A_l = −U⁻¹ and the two
signs cancel. So 𝒲 = 𝒩·U, and multiplication by an invertible U is a linear
isomorphism: dim 𝒲 = dim 𝒩.

**(⟸), and the trick is to grade the OTHER algebra.** The obvious move is to
grade 𝒜 by the parity of the A-letter count and hope the sum is direct, which it
need not be, and that route costs a case split and a recursion into corners.
Grade

    ℬ = 𝒩 ⊕ 𝒲

instead. The four products behave
because the signs simply multiply: 𝒩𝒩 ⊆ 𝒩, 𝒩𝒲 ⊆ 𝒲, 𝒲𝒩 ⊆ 𝒲, 𝒲𝒲 ⊆ 𝒩. Both
summands are closed under the adjoint, since H is Hermitian and A_l X A_l = ±X
survives it, and 1 ∈ 𝒩. And 𝒩 ∩ 𝒲 = 0 as soon as **at least one jump is
present**, because X = −X there. So ℬ is a unital ℤ₂-graded *-subalgebra of
M_d(ℂ), hence **semisimple** (a finite-dimensional *-subalgebra of M_d(ℂ) has no
nonzero nilpotent *-ideal, since x*x nilpotent forces x = 0; the step is
standard and is named because Lemma 3 is where it is spent), and

    β(X + W) = X − W

is an order-2 *-automorphism of ℬ, with no case analysis anywhere.

β permutes the simple blocks of ℬ as an involution. A SWAPPED pair M_n ⊕ M_n
contributes n² to each part. A FIXED block M_n has β = Ad(u) with u² = 1 after
rescaling, and ±1 multiplicities p + q = n, contributing p² + q² to the even
part and 2pq to the odd one. Hence

    dim 𝒲 = Σ_pairs n² + Σ_fixed 2p_k q_k  ≤  Σ_pairs n² + Σ_fixed (p_k² + q_k²)
           = dim 𝒩 ,

with equality exactly when p_k = q_k in every fixed block. And there an
invertible odd element exists blockwise: [[0,1],[1,0]] in the u-eigenbasis on a
fixed block, and (x, −φ(x)) with x invertible on a swapped pair, where β acts as
(x, y) ↦ (φ⁻¹(y), φ(x)) for an isomorphism φ of the two copies. ∎

That is the whole of Lemma 3, both halves, with no ideal, no corner and no
recursion, and it is what turns the pairing-condition page's sampled predicate
into a decided one: on every row that page reports as nonempty-but-singular the
palindrome breaks, so by the theorem `dim 𝒲 < dim 𝒩` there (inferred, not computed
row by row), and the inequality **proves** that no invertible element is hiding
there.

## (e) The theorem

> **Theorem.** With the hypotheses of §(a), the spectrum multiset of L is closed
> under λ ↦ −λ − 2σ if and only if `dim ker L = dim ker(L + 2σ)`.

**Necessity.** L is unital, so 0 ∈ spec(L) with the identity in its kernel. The
reflection carries the eigenvalue 0 to −2σ with the same multiplicity, so
alg(0) = alg(−2σ). By Lemma 2 both are geometric, so dim ker L = dim ker(L+2σ).

**Sufficiency.** By Lemmas 1 and 3, equal dimensions give an invertible U with
[U, H] = 0 and U A_l U⁻¹ = −A_l. One-sided multiplication ℒ_U(ρ) = Uρ then gives

    ℒ_U L ℒ_U⁻¹ = −L† − 2σ

directly, the computation carried out in
[The pairing condition](../../experiments/THE_PAIRING_CONDITION.md): only the
left factor is conjugated, the commutator term is untouched, the dissipator's
U A U⁻¹ ρ A becomes −AρA, and the surviving −σρ against the +σρ of −L† is where
the −2σ comes from. Within §(a)'s class, where the Hermitian jumps make the dissipator self-adjoint, it needs U invertible and nothing more; for a unitary
involution U it is the PT condition of Prosen (2012), see
[the literature review](../LITERATURE_REVIEW.md#pt-symmetric-liouvillians-prosen-2012). Similar operators
have equal spectra with multiplicity, so
spec(L) = spec(−L† − 2σ) = {−λ̄ − 2σ : λ ∈ spec L}. And L preserves hermiticity,
so its characteristic polynomial is real and spec(L) is closed under
conjugation, which collapses that to spec(L) = {−λ − 2σ}. ∎

The conjugation-closure step is stated, not gated, and cannot be: GF(p) has no
conjugation. It is one line and it is where the reflection about the vertical
line Re = −σ becomes the palindrome about the point −σ.

**On the characteristic-polynomial form.** With D = d² the degree,
p(−x − 2σ) = (−1)^D p(x) whenever the multiset is reflection-closed, so the
identity the repo's exact kernel tests, `p(x) = p(−x − 2σ)`, is the multiset
statement **for even d**, which is every qubit register. At odd d the sign flips
and the identity is unsatisfiable, which is consistent with §(f5): at odd d, with
at least one jump, the multiset statement fails too.

## (f) What follows

**(f1) The kernel always dominates, and the far space is not "the fastest
modes".** `dim ker(L + 2σ) ≤ dim ker L` for every H, every unitary-Hermitian
jump set and every positive profile. This is Lemma 3's inequality and it holds
with no palindrome anywhere in sight.

The tempting name for `ker(L + 2σ)` is *the fastest-decaying modes*, and it is
wrong twice. First, when 𝒲 = 0 the set is empty while the spectrum still has a
leftmost point: a [GLOSSARY](../GLOSSARY.md) row measures a generic Hermitian H
whose maximum rate falls short of 2σ entirely. Second, and less obviously, even
when the value −2σ IS attained on the line Re = −2σ, `ker(L + 2σ)` need not
contain those modes: §(b)'s argument applies to any eigenvector with
Re λ = −2σ and forces the anticommutation, but it leaves −i[H,W] = iθW with
θ ≠ 0 available, so an oscillating mode can sit at the left edge of the strip
without lying in the kernel at the point −2σ. The smallest example is one
qubit with H = Z and A = Z: `ker(L + 2σ) = 0` while two eigenvalues sit at
Re = −2γ, at −2γ ± 2i. The honest name is **the non-oscillating modes at the
left edge**, and the object the theorem is about is the eigenvalue −2σ, not the
line through it.

**(f2) The canonical chain.** For the Heisenberg chain under Z-dephasing on
every site, 𝒩 is spanned by the N+1 magnetization-sector projectors and 𝒲 by
the N+1 XOR-sector modes. Both counts are on one page,
[DEGENERACY_PALINDROME](../../experiments/DEGENERACY_PALINDROME.md) Result 2,
which states the bijection between them outright and reads it as a consequence
of Π. The theorem inverts that reading: the bijection is not downstream of the
palindrome, it **is** the palindrome, and it decides cases where Π is not
available.

**(f3) The criterion cannot see γ, inside the open orthant.** Neither defining
condition of 𝒩 or 𝒲 mentions γ_l, so the theorem predicts, before any run, that
the palindrome cannot depend on the rate profile at all as long as every rate is
strictly positive. The non-uniform block gates that more sharply than it looks:
`row_spaces` takes the bond couplings and the letters and DROPS the rate
profile, so on those 600 rows the criterion is computed with no γ information
whatever while the palindrome side receives γ = (3, 7, 11)/100, and they agree
on all 600. What it does not say: the SPECTRUM depends on γ throughout,
the centre −σ moves with it, and at the boundary γ_l = 0 the verdict can
move, because a zero rate can drop a condition (gate 10b). It is the pairing, and
only inside the open orthant, that is blind.

**(f4) The anti-steady mode, and F103's three rows.**
In [PROOF_F111](PROOF_F111_HARD_CELL_PURE_D_TEMPLATE.md) the pure-D templates
lose the anti-steady mode; by Lemma 1 they lose it exactly when no operator
anticommutes with every jump while commuting with H, which translates the
question from the spectrum to the letters, and F111 reads the why on the letters
through Theorem 2 of the complement connection, §(f10)'s route.
[F103](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) §7.12's three rows, soft with
non-bipartite basis-state graphs and therefore beyond any 2-colouring, are
decided correctly by the count at N = 3, 4 and 5. The operator the criterion
returns there is worth naming precisely: 𝒲 is **one-dimensional** and spanned by a single Pauli string,
Y^⊗N for XX+XZ and XX+XZ+ZX and X^⊗N for YY+YZ, and the X^⊗N of that last
family is F1's own Π². For the two 2-term rows this is an agreement test against
a routing the repo classifies bit-exactly; for XX+XZ+ZX, past the 2-term routing
§7.12 names as what remains, the count supplies the element Y^⊗N itself.

**(f5) Every jump must be traceless, and at odd d none can be.** An invertible
U with U A_l U⁻¹ = −A_l makes A_l and −A_l similar, so the ±1 eigenvalue
multiplicities of each A_l must agree: **every jump operator is traceless**. That
is a necessary condition which is H-free, costs one trace, and is strictly
stronger than a parity count on d, because it bites at even d too:
A = diag(1, 1, 1, −1) at d = 4 kills the palindrome for every Hamiltonian, and
nothing in the gated scope exercises it, every jump at even d there (n̂·σ⃗, Pauli
strings, depolarizing sites) being automatically traceless. Since a traceless
involution needs an even dimension, odd d admits no such A_l at all, so at odd d
**with at least one jump** the spectrum never pairs; gate 10c checks the inequality
dim 𝒲 < dim 𝒩 that carries this at d = 3 and d = 5, a check of Lemma 3 at those
dimensions rather than of the spectrum.

The hypothesis is not a formality. With NO jump, 𝒩 and 𝒲 are the same space, and at d = 3 with
H = diag(1, 2, 5)/7 the multiset IS reflection-closed, exactly. Note also that
the characteristic-polynomial FORM proves nothing here: p(x) = p(−x − 2σ) is
unsatisfiable for every monic polynomial of odd degree, whatever the physics.

And the condition generalises for free, because U commutes with H as well:
**every word in {H, A_1 … A_m} with an odd number of A-letters is traceless**.
That is an O(d³) pre-filter with no Liouvillian in it, and the converse-side
companion to §(d)'s rank comparison. The pre-filter is built, in exact
arithmetic: `fw.odd_word_obstruction`
([`simulations/framework/diagnostics/f158_odd_word.py`](../../simulations/framework/diagnostics/f158_odd_word.py))
finds the shortest such word with a nonzero trace, multiplying Pauli strings over
the Gaussian rationals, so a nonzero trace is a certificate and not a reading. On
the N = 3 census of [the palindrome as a colouring](../../experiments/THE_PALINDROME_AS_A_COLOURING.md)
(36,288 rows) it never fires on the 9,609 rows this theorem's two ranks call
palindromic, and it certifies 24,671 of the 26,679 broken rows with one jump letter
(Tr(H^k·A), k ≤ 4) and 25,831 with up to three, each by an exact trace, so those
verdicts need no rank at all
([`simulations/f158_odd_word_prefilter.py`](../../simulations/f158_odd_word_prefilter.py)).
For a qubit chain dephased along Z on site l the shortest word is already
Tr(H·Z_l) = 2^N·δ_l, δ_l the coefficient of Z_l in H: a longitudinal field on a
dephased site rules the full palindrome out, a uniform one included, while the
palindrome of the decay rates alone can survive it (Heisenberg under a uniform
field). A finite budget of traceless words decides nothing; over every word the
condition is also sufficient, since the traces of all words determine a
finite-dimensional *-representation up to unitary equivalence, so (H, A) and
(H, −A) would be unitarily equivalent. The live witness `inspect --root
twoendstrings` reads the one-jump words Tr(H^k·A) for k ≤ 4, exact over the
Gaussian integers, beside the string-span counts of the second corollary to
Lemma 1; MirrorWorld's [`EndCount`](../../compute/MirrorWorld/EndCount.cs)
carries the words and the order of `fw.odd_word_obstruction`, the three-jump words
included on request (its default reads one jump letter per word, where the
Python's also reads three).

**(f6) The two ends are both what they are called.** For this class L† is L with
H ↦ −H, so `ker L† = ker L = 𝒩`: the conserved quantities and the steady modes
are the same subspace here, which is what makes the join in (f2) safe. The
distinction still has teeth in general and the repo keeps it
([PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA](PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md):
*"That kernel is a complex vector space; the manifold of stationary STATES
inside it, positive and of unit trace, has one dimension fewer"*), and
everything computed here is the kernel.

**(f7) The pairing page's necessity rests on no census.** That page scores
140,861 rows, and its necessity, palindrome ⟹ an invertible U, does not lean on
them: §(e) gives the equal counts and §(d)'s graded algebra ℬ = 𝒩 ⊕ 𝒲 turns
equal counts into an invertible element by algebra alone.

**(f8) Sufficiency gives the Jordan structure too, and with it the arc's S.**
The criterion is stated for the multiset, but the sufficiency argument produces
a SIMILARITY, ℒ_U L ℒ_U⁻¹ = −L† − 2σ. In the Hermitian basis L is a real matrix,
so L ~ conj(L); and M ~ Mᵀ for every square M; so L† ~ conj(L) ~ L. Chaining,

    L  ~  −L − 2σ ,

so where the criterion holds the ENTIRE Jordan structure is reflection-symmetric
and not merely the eigenvalue multiset. That also answers the request the arc
`f138_converse_failures` has carried since 2026-08-03, for an S with
S L S⁻¹ = −L − 2σ: such an S exists whenever the two counts agree. What §(e)
exhibits is the one-sided U giving −L† − 2σ; the S is a derivation on top of it
and is not gated here.

**(f9) Noise in the system's own eigenbasis never pairs.** If every jump is a
FUNCTION of the Hamiltonian, A_l = f_l(H), then 𝒲 = 0 and the spectrum never
pairs, at any rate profile and whatever H is. One line: W ∈ 𝒲 commutes with H,
hence with f_l(H) = A_l, so A_l W A_l = A_l² W = +W, which contradicts
A_l W A_l = −W unless W = 0. Meanwhile 𝒩 is as large as H's own commutant, so
the two counts are as far apart as they get. A = sign(H − c), with c not an
eigenvalue of H, which is Hermitian and squares to 1, is such a jump squarely INSIDE the class.

The physical reading is worth the sentence, because it is the sharpest thing the
criterion says without any computation at all: **a channel that watches in the
system's own energy basis cannot produce a palindrome**. The pairing needs the
watching to be transverse to the turning in the strong sense of anticommuting
with it, and a function of H is the exact opposite of that.

**(f10) With every site dephased, the criterion is a flat connection.** When every site carries exactly one single-site Pauli jump, turning each jump into Z makes the dark span the diagonal matrices and the lit span X^⊗N times them, and the far end becomes the flat sections of the connection d_y / d_x = H_xy / H_x̄ȳ on H's hopping graph over bitstrings: the palindrome holds exactly when every component carries one. For Heisenberg bonds on a connected graph under uniform Z, with fields along letters, it follows that the palindrome holds exactly when no field lies along Z and the fields do not use both X and Y, F138's converse holding in that class; under mixed axes with Heisenberg bonds on a connected graph, exactly when two axes occur and every field lies along the third letter, the near end one-dimensional; and with the dephased sites each under one single-site jump and every undephased site keeping a letter H conserves, both ends split into sums over pairs of sectors of the same scalar problem; and under Heisenberg bonds on a connected graph with letter fields and exactly one undephased site, the palindrome holds exactly when a colouring exists. Proof and gates: [PROOF_PALINDROME_COMPLEMENT_CONNECTION](PROOF_PALINDROME_COMPLEMENT_CONNECTION.md).

**(f11) Where the counts differ, a Hermitian far element splits off a part that pairs.**

<!-- STAGE 0 SWEEP RECORD for (f11) (for us, not for the reader)

The sweep went to every store by name through one scout and by hand, by the primitive (a nonzero singular anticommutant; an invariant piece that pairs inside a broken spectrum; a map with a kernel against the reflected generator) and by the words that might name it, and two rounds of fresh readers checked it at source. This proof: gate 11 builds rows strictly between the counts as H₁ ⊕ H₂ with disjoint spectra, and gate 12's one strict row is that construction at d = 32; §(g) counts 212 strict rows in the sweep; Lemma 3's display gives dim 𝒩 and dim 𝒲 over its fixed blocks and swapped pairs, the difference being Σ (p_k − q_k)² over the fixed blocks; §(f5) reads the whole palindrome as (H, A) and (H, −A) unitarily equivalent; §(e) multiplies on the left, ℒ_U(ρ) = Uρ, where step 3 multiplies on the right; §(h) says the criterion says nothing about the Jordan structure where the counts differ. Nothing in it read what a singular far element does. The registry: F22 puts every X^⊗N·P_k at −2Σγ under [H, X^⊗N] = 0 and [H, Σ_l Z_l] = 0, the GHZ coherence alone wherever |0…0⟩ and |1…1⟩ are eigenstates of H, at −2Σγ − i(E_{0…0} − E_{1…1}), and under a uniform field shifts X^⊗N·P_k by 2ih(N − 2k), so that at even N the element X^⊗N·P_{N/2} stays at −2Σγ: a far element whose range is the half-filling shell (a reviewer's exact end count, Heisenberg chain N = 4, (5, 1), rank 6), the instance behind HeisenbergColourLawClaim's remark that a uniform field along Z at even N leaves one non-invertible far element, neither citing the other; F167's invariant hulls pair about their own rungs, which at J₀ = 0 equal σ, so there a part pairs about −σ inside a spectrum F158 breaks with 𝒲 = 0 (the jumps leave one lit string, Q, and F167's entry has [H, Q] ≠ 0), a part no far element reaches, the room the Corollary's "at least" leaves; THE_SODIUM_CHAIN leaves open whether one statement holds its hulls with F157's blind subspace, the undamped space and F89's AT lock, naming F108's hull beside them; F166's blockwise centres sit in §(h)'s paragraph; F121's qudit partial palindrome is another object, and CAUGHT_ERRORS (2026-11-07) records the word misused once. The proofs: the complement connection's Theorem 2 gives one far dimension per good component, its Theorems 4 and 7 write both ends as sums over pairs of charge sectors, with a charge chosen by hand where step 4's comes from the far element, and its Theorem 6 (a) leaves a far end at a singular chord (CAUGHT_ERRORS 2026-10-02); PROOF_ONE_DEPHASED_SEAT writes the invariant decomposition and the two ends for S (step 7 of the section with a planar field on d) and for the double singlet at N = 5, holds at N = 4 a singular far element on a row that pairs (live in DephasedSeatFlipBlockWitness), and names the row at J = (−1, 1, −1) that the gate reads; PROOF_R90_FROZEN_DIVISOR §9.1 holds a chiral sector symmetric about the frozen root, and its index count (MirrorWorld's Divisor) is a count of the same kind as V^⊥'s, neither citing the other on this point. The experiments: THE_PAIRING_CONDITION counts rows whose admissible subspace yielded only singular draws and closes them by this theorem, none hiding an invertible element; the colouring page holds rows one mode short at the far end and, on rows that pair, sector pieces that are singular far elements whose sum is the colour, this split read on a row that pairs, and the polar part of a far element; XOR_SPACE reads the XOR sector as the far end of the axis; LATTICE_OPENING_LAW reads the cat sector as H-dead; THE_SODIUM_NUCLEUS_TWO_MIRRORS holds a row with near 4 and far 2; SOFTNESS_IS_N_DEPENDENT reads a strict row at N = 6 (near 8, far 4) through Theorem 2; THE_MIRROR_UNDER_A_SHARED_FIELD, with its twin in MIRROR_SYMMETRY_PROOF, pairs the block (p, q) with (q, N − p) about a centre moved by h(2p − N) under a uniform field, so at even N the blocks with a popcount N/2 keep −σ, the same split on F22's uniform-field row by another route. The OpenArcs registry: the retired arc benzene_center_tier_upgrade measures uniform-field rows pairing in part (the C4 ring 182 of 256, the 4-chain 157 of 256, the C3 ring none), above the Corollary's bound of 156 on the first two (a reviewer's exact end count: (5, 1), rank 6) and at its 0 on the third; nothing else on singular far elements; the retired golden-router arc's singular strata (PROOF_CEILING_GOLDEN_ROUTER §5) are singular solutions of W·L = (−L − 2σ)·W, maps with a kernel against the reflected generator, not far elements. hypotheses/ and recovered/: nothing; reflections/ON_THE_WAVE_RECOGNIZING_THAT_IT_IS_NOT_ALONE holds, in its sweep comment, far elements that are sums of strings, the Carrier, the Stitch's carrier and F22's far sector. fw.Confirmations: nothing. The glossary: no headword for far element or anticommutant, while its rows on the XOR space and modes and on the colouring hold the object under those names. CAUGHT_ERRORS also: 2026-10-08 (only rows strictly between can give a false negative), 2026-10-16 (a left multiplication remembered as the right one), and 2026-11-07's R∘Φ = S⁺∘R, a reflection with a kernel. The typed layer: PalindromeTwoEndCountClaim and its witnesses count the two ends; TrianglePalindromeLocusWitness names S's far element singular; DephasedSeatFlipBlockWitness holds the N = 4 element and names the double singlet's; SodiumChainWitness holds the hulls; in MirrorWorld, Sextant.Intertwines reads a map with a kernel and keeps it out of Decide (step 3's map is one, which it may read on the cat; not run), Lattice holds both one-sided readings by X^N and GammaFold docks the left one, EndCount.FarElement returns an element where the far count is one and IntertwinerRank counts the far end as intertwiners, Carrier holds one row's single far element and CarrierTests decline the S row. The literature: Buča and Prosen 2012 (arXiv:1203.0943) define a strong symmetry as a unitary S with [S, H] = 0 and [S, L_m] = 0 for every jump and split the generator into its invariant blocks; KMS_DETAILED_BALANCE and LITERATURE_REVIEW state both kinds in the paper's sense, MIRROR_SYMMETRY_PROOF reads X^⊗N, which commutes with the XXZ Hamiltonian, under Z jumps as weak and not strong, the colouring page reads a colouring as weak and not strong and the product of two as strong, THE_SODIUM_CHAIN reads ⊗ZZ_k, lit to jumps, as weak where H commutes with it, and Z_N_PARTNERSHIP, which THE_ONE_SQUARE §9's table cites, uses the paper's sense. Joined in the text: §(d), §(e), §(f5), gates 11 and 12, F22's GHZ coherence, the seat proof's S and double singlet, the corner shift of THE_ONE_SQUARE §9, Buča and Prosen's strong symmetry. Adjacent and not joined: F22's uniform-field element with HeisenbergColourLawClaim's remark, the shared-field page's blocks at popcount N/2 and the benzene arc's counts, F167's hulls and the sodium chain's open question, F140 §9.1 and the Divisor's index, Theorems 2, 4 and 7, the colouring page's rows one mode short and its sector pieces, the sodium row, the Lattice's one-sided readings and LATTICE_OPENING_LAW's cat, Sextant.Intertwines and R∘Φ, the golden router's singular strata. Words: V for the range, since R is the right multiplication by X^⊗N of the Π factorization; block is the joint-popcount block and corner is in use for many objects, so the parts are the system on V, the system on V^⊥ and the two cross parts, and the gate's pieces are the operators whose bra lies in a given range; the noun watching is retired (THE_LABEL_MAP), the bits here are dephased.
-->

Call the elements of 𝒲 far elements and those of 𝒩 near elements. Lemma 3 leaves no invertible far element on a row with 0 < dim 𝒲 < dim 𝒩, and the theorem then says the multiset breaks. The far elements that do exist still say where part of the spectrum pairs. 𝒲 is closed under the adjoint, H and every A_l being Hermitian, so it is spanned by Hermitian elements; let W ∈ 𝒲 be Hermitian, V its range and r = dim V.

1. *V and V^⊥ are invariant under H and under every jump.* H·W·v = W·H·v and A_l·W·v = −W·A_l·v, so H and each A_l map V into itself, and being Hermitian they keep V^⊥ too. Equivalently, the projector P_V onto V commutes with H and with every jump, so every unitary e^{iφP_V} lies in 𝒩 and is a strong symmetry in the sense of Buča and Prosen (§(a)), and for e^{iφ} ≠ 1 the four parts below are the invariant subspaces into which it splits the generator. ℋ = V ⊕ V^⊥ splits L into four L-invariant parts, mutually orthogonal: the system on V, the system on V^⊥, and two cross parts, the operators with ket in V^⊥ and bra in V and those with ket in V and bra in V^⊥.
2. *The system on V pairs.* On V, W is invertible (it is Hermitian with range V), commutes with H restricted to V and anticommutes with every A_l restricted to V, which are Hermitian and square to the identity there. The theorem applies to the system on V, and its r² eigenvalues pair about −σ.
3. *The two cross parts pair together.* For every operator ρ,

       L(ρW) = [(−L† − 2σ)(ρ)]·W ,

   since the commutator term passes through W and A_lρWA_l = −(A_lρA_l)W. §(e) makes the same computation on the left, with an invertible element; here the element multiplies on the right and is not inverted, which is the corner shift ρ ↦ ρ·F of the table of moves in [THE_ONE_SQUARE](../THE_ONE_SQUARE.md) §9 with F a far element, a sum of strings that may have a kernel, in place of one Pauli string. Right multiplication by W carries the cross part with bra in V one to one onto itself, so there L is similar to −L† − 2σ; the parts being orthogonal and invariant, L† restricted to that part is the adjoint of L restricted to it, and the spectrum of the part is closed under λ ↦ −λ̄ − 2σ. The other cross part carries the conjugate spectrum, L(ρ†) = L(ρ)†, so the two together are closed under λ ↦ −λ − 2σ: 2r(d − r) eigenvalues more.
4. *The difference of the counts sits in V^⊥.* Each end of L is the sum of the ends of the four parts. On V the two agree (step 2). On a cross part, multiplication by W on the side of V takes near elements to far ones and far ones to near ones (near times far is far and far times far is near, §(d)), injectively both ways since W is invertible on V, so the two ends of a cross part have one dimension. So dim 𝒩 − dim 𝒲 is the same for the whole and for the system on V^⊥.
5. *At maximal rank nothing far is left.* If W has maximal rank among the Hermitian elements of 𝒲, the system on V^⊥ has no far element: a Hermitian one, W′, extended by zero on V, would be a far element of the whole, and W + W′ a Hermitian element of 𝒲 of rank r + rank W′.

> **Corollary.** Let V be the range of a Hermitian element of 𝒲 of maximal rank r (r is fixed by H and the jumps, V need not be). Then V and V^⊥ are invariant under H and every jump, the system on V pairs about −σ, so do the two cross parts together, and the system on V^⊥ has no far element and carries the whole difference dim 𝒩 − dim 𝒲. Exactly d² − (d − r)² eigenvalues lie in these parts that pair, so at least that many pair, and the spectrum pairs exactly when V = ℋ.

"At least", because the system on V^⊥ breaks only as a whole: its near end holds its own identity and its far end is empty, while some of its eigenvalues may still pair among themselves (on the gated swap-symmetric row four of them do). In §(f5)'s terms, which read the whole palindrome as (H, A) and (H, −A) being unitarily equivalent, the polar part W·|W|⁻¹ on V, itself a far element since |W| is near, is such an equivalence on V, and on V^⊥ none exists. [PROOF_ONE_DEPHASED_SEAT](PROOF_ONE_DEPHASED_SEAT.md), on a Heisenberg triangle of sites d, u and v with the jump Z on d alone and its classes named A, S and so on, writes the invariant decomposition and the two ends for its own rows: on the singlet class S outside A, whose singular far element (f̂_d·σ⃗_d) ⊗ P_s stands against a near end holding 1 and 1_d ⊗ P_s, H keeps ℂ²_d ⊗ (the singlet of u and v) and its complement ([the section with a planar field on d](PROOF_ONE_DEPHASED_SEAT.md#the-triangle-with-a-planar-field-on-d), step 7), and on the double singlet at N = 5 with the fields equal within each pair both ends split between a plane and its complement. Gate 11, and gate 12's one strict row, build the split as H₁ ⊕ H₂ with disjoint spectra, which leaves the cross parts without near or far elements; a general row need not have that, and on the gated row with z fields −1 and +1 on u and v, no field on d, two of the row's four far elements lie in the cross parts, one in each.

On the cube the step that does the work is the third. Under Z dephasing a cell |x⟩⟨y| is charged 2γ_l for every dephased bit l on which x and y differ. Right multiplication by a far element, which anticommutes with every Z_l, exchanges on each dephased bit the cells where ket and bra agree with those where they differ, so the price p becomes 2σ − p, while the commutator does not notice, the element commuting with H. An invertible element moves every bra, which is the whole palindrome; one with a kernel moves the bras in its range. The cat shows it plainly: wherever H keeps |0…0⟩ and |1…1⟩ as eigenstates of one energy ([F22](../ANALYTICAL_FORMULAS.md#f22-ghz-maximum-disagreement-support-tier-2-verified-n2-5), its GHZ coherence), W = |0…0⟩⟨1…1| + |1…1⟩⟨0…0| anticommutes with every Z_l and commutes with H, V = span{|0…0⟩, |1…1⟩}, and H on V is a multiple of the identity. A cell |x⟩⟨0…0| pays for the dephased ones of x, the cell |x⟩⟨1…1| for its dephased zeros, and their prices add up to 2σ. Two further kinds on the triangle with one dephased seat: the singlet class S, whose far element outside A is X_d ⊗ P_s with f_d along x, has V = ℂ²_d ⊗ (the singlet), with H on V a field perpendicular to the jump's axis plus a multiple of the identity; and on the gated swap-symmetric row with fields along z, V is spanned by |0⟩_d|11⟩ and |1⟩_d(|01⟩ + |10⟩)/√2 with equal diagonal entries of H in that basis, coupled only through d's hop.

Gated in [`far_element_mirror_gate.py`](../../simulations/far_element_mirror_gate.py) → [`far_element_mirror_gate.txt`](../../simulations/results/far_element_mirror_gate.txt), exact over ℚ(i) on eight rows strictly between the counts: the triangle with one dephased seat (the cat at fields along z summing to zero, a swap-symmetric row with fields along z whose two-dimensional V is resonant, S as [the section with fields off the plane](PROOF_ONE_DEPHASED_SEAT.md#the-triangle-with-fields-off-the-plane) reads it, with f_u = f_v carrying z parts and a field on d along x, and z fields −1 and +1 on u and v at couplings (−1, 1, −1) with no field on d), an XY chain of three with every site dephased at unequal rates and fields along z summing to zero, three qubits in four parts with two of them carrying an invertible anticommuting element, and two qubits in two parts, the construction of gates 11 and 12, with one jump and with two at unequal rates. For every far element read (ranks 1, 2 and 4, its range equal to that of W† or not): the invariance of its range and of the range of W†, the identity of step 3 on every matrix unit, and the line reflection between the operators whose bra lies in its range and those whose bra lies in the range of W†, as an identity of characteristic polynomials; neither the invariance nor the identity uses that W is Hermitian, and right multiplication by W carries the first of those two pieces one to one onto the second. For a Hermitian far element of maximal rank among those read: the system on V pairs, the two cross parts pair together, and the system on V^⊥ has the whole difference and no far element; the near and far counts of each cross part, and how many eigenvalues of the system on V^⊥ pair among themselves, are printed. That the characteristic polynomials of the parts multiply to that of L is checked too, but follows from the invariance and is no independent evidence. The controls go through the same doors: a basis line kept by every jump but not by H fails the invariance test on every row, W + 1 breaks the identity on every row, the pieces of the singular near element the gate picks do not reflect (four rows), the operators whose bra lies in the complement of a far element's range do not reflect onto themselves (every row), on the two rows with several jumps at unequal rates the reflection about the largest single rate instead of σ fails, the point reflection fails on 19 of the 22 far elements read, and on the four-part row a Hermitian far element of rank 2 leaves a far element in the system on its complement. Argued, not gated: steps 1 to 5 as statements about every row, and the cube's reading, the forms of V and of H on V in the paragraph above among it. Not live in the typed layer.

## (g) What is gated, and how

Most of the above is checked from below in
[`simulations/f138_rank_criterion.py`](../../simulations/f138_rank_criterion.py),
and the exceptions are named rather than covered by a blanket: §(f1)'s naming
argument, §(f6)'s `ker L† = ker L`, and §(d)'s graded-algebra step (the
Wedderburn decomposition of ℬ = 𝒩 ⊕ 𝒲) are reasoned and not gated; the
step's conclusion, invertible element ⟺ equal counts, is gated on 5,000 rows.
Also not gated, each checked by hand or by a scout only: §(f9) (noise that is a
function of H), the qutrit ⊗ qubit row at d = 6, §(f5)'s `diag(1, 1, 1, −1)` at
d = 4, the DEPOLARIZING Y-field row, and the sufficiency of §(f5)'s condition
over every word.
§(f5)'s odd-d rows vary H while holding A fixed at `diag(1, …, 1, −1)`.
The gate is built in **three** layers, deliberately, and the layers encode the
conditions differently (the exact layer as `A⊗Aᵀ ± 1`, the modular layer as a
commutator or anticommutator, the float layer from scratch), so a coding error
in one does not hide in the others.

- **Exact over ℚ(i), Fraction arithmetic, `== 0` and no tolerance:** Lemma 1 in
  both directions on 7 named rows, by taking the exact kernel and demanding that
  EVERY basis element satisfy the two operator conditions; and §(b)'s
  inequalities as exact rational comparisons. This layer covers N = 3, uniform
  γ = 1/20 and single-site Pauli letters only; no off-axis direction, no Pauli
  string, no N = 4 and no non-uniform profile is ever checked exactly.
- **Exact over GF(p):** semisimplicity as `rank(M²) = rank(M)`; Lemma 3's module
  identity `𝒩U = 𝒲`; and the scoring against the palindrome. **Prime counts
  differ by block and the difference is not cosmetic**: the palindrome side of
  the 13,540-row block uses the committed three-prime verdict of
  [`f138_clause_two_sweep.py`](../../simulations/f138_clause_two_sweep.py),
  while every nullity in the scoring, and both sides of the 1,875 rows in gates
  7 and 8, run at one prime. A nullity read mod p can only come out too
  **large** (reduction is a ring map, so a rank can only come out too small), so
  both nullities are upper bounds and a mod-p equality is evidence, not proof.
  The direction is USUALLY the safe one: a bad prime inflating ONE nullity
  either leaves the verdict as it was or flips it, and a flip makes a gate FAIL
  rather than pass. The exception is a prime
  inflating BOTH sides compensatingly, which would pass wrongly, and that is the
  residue the three-prime blocks reduce and the one-prime blocks do not.
- **Float, sharing no construction code with any of it:** gate 9 rebuilds the
  whole object in dense complex numbers, calling this file's `row_spaces` only as
  the other side of its cross-check, and compares its two nullities row by row
  against the modular route, its verdict only against its own nullities.
  Its two thresholds are gated on measured separations rather than chosen.

**Sampling has not vanished from the script**, and §(d) does not claim it has.
Three elements remain and they do NOT all fail safe. `invertible_in` draws
random elements and tests a determinant, load-bearing in the gates that check
Lemma 3 itself (gates 3 and 3b); there the direction is safe, since a missed
invertible element makes a gate fail rather than pass. Gates 10c and 11 also
report what it finds, and there a missed element would pass, so their "no
invertible element" is the search agreeing with Lemma 3, which the deterministic
`dim 𝒲 < dim 𝒩` they check already implies. `gate2b` draws six random vectors per row
for the two inequalities. And the palindrome verdict itself is a random-point
determinant test, in the committed three-prime route and in the single-prime
one alike, which carries the one-sidedness its own source module states in a
line this file should carry too: *"False is a PROOF of a break; True is a
certificate."* A wrongly-True palindrome on a row where the criterion also says
yes records no false positive, so that one can make a gate pass rather than
fail. The probability is negligible (degree at most 1024 against p near 10⁹, six
points, three primes on the largest block) but the direction is the unsafe one
and saying so is the point. A second way to a wrong True is a bad prime, at which
p(x) − p(−x−2σ) vanishes identically mod p although it does not over ℚ(i); more
points cannot catch that, only more primes, and the 1,875 rows of gates 7 and 8
are read at one. What §(d) removes is the sampling from the
CRITERION, which is now two ranks.

The gate scores dim 𝒩 and dim 𝒲 from the operator equations, the commutant and
the jump-signed anticommutant; that these are dim ker L and dim ker(L + 2σ) is
Lemma 1, checked exactly on gate 1's seven rows and in floats on gate 9's 885.

**All 131 gates pass.** The criterion is scored against the palindrome on
**15,415 rows**, in both directions, of which **2,596 hold and 12,819 break**
(a false positive could have occurred on any holding row; a false negative
only on a breaking row with 0 < dim 𝒲 < dim 𝒩, the 212 strict rows counted
below, since on the others the identity in 𝒩 forces the counts apart and the
unital L cannot pair without an element of 𝒲):

| what is scored | rows | holds | FP | FN |
|---|---:|---:|---:|---:|
| P₃ / K₃ / bond+isolate / P₃ generic, N=3 letter grids | 10,000 | 1,503 | 0 | 0 |
| three single- and two-letter bond sets | 2,700 | 832 | 0 | 0 |
| non-uniform J per bond and γ per site | 600 | 77 | 0 | 0 |
| N=4 path and ring | 240 | 6 | 0 | 0 |
| off-axis dephasing, `n̂·σ⃗` at rational directions | 792 | 36 | 0 | 0 |
| multi-site Pauli-string jumps, one and two of them | 900 | 47 | 0 | 0 |
| all 1-, 2- and 3-term bond words (F87 territory) | 174 | 86 | 0 | 0 |
| F103 §7.12's three rows, N = 3, 4, 5 | 9 | 9 | 0 | 0 |
| **total** | **15,415** | **2,596** | **0** | **0** |

**217** of the rows have no jump at all (the letter grids include the empty
dephasing pattern): they sit outside the theorem's class, where σ = 0 and
𝒩 = 𝒲, and hold trivially, so within the class the count is 15,198 rows of
which 2,379 hold. The tally gate counts them.

On **212** of those rows `0 < dim 𝒲 < dim 𝒩`, which is where the criterion's
equality says more than "the anticommutant is nonempty"; see the third caveat
below.

Beside the scoring, on rows outside the 15,415 tally:

| what is checked | scale | result |
|---|---|---|
| Lemma 1, both inclusions, exactly over ℚ(i) | 7 named rows | every kernel basis element satisfies both operator conditions, `== 0` |
| Lemma 2, semisimplicity at both ends | 7 named rows × 3 primes, 250 grid rows | `rank(M²) = rank(M)` everywhere; 33 of the 250 carry a kernel at −2σ, 29 of them with a jump (on the other 4, −2σ = 0), so the check is not vacuous |
| §(b)'s inequalities, exact rationals | 4 named rows × 6 random vectors | dissipative at 0, accretive at −2σ (sampled vectors, not a per-row proof) |
| Lemma 3 (⟹), `𝒩·U = 𝒲` | 5 named rows | spans exactly, dimensions equal; the other 2 named rows have no invertible element and assert the contrapositive of (⟸), `dim 𝒲 < dim 𝒩`, instead |
| Lemma 3, invertible element ⟺ equal counts | 5,000 rows | 0 mismatches |
| §(f1), `dim 𝒲 ≤ dim 𝒩` | 3,600 rows | holds, largest observed gap 12 |
| §(f2), the canonical chain | N = 2, 3, 4 | both counts N+1, palindrome holds |
| §(f5), odd d | d = 3, 5, 24 random Hermitian H and 4 multiples of 1 (where dim 𝒲 > 0) | `dim 𝒲 < dim 𝒩` always, so by Lemma 3 no invertible element (a sampled search agrees) |
| gate 10, depolarizing and the rate fence | N = 2, 3 | criterion tracks the palindrome through both |
| gate 11, `0 < dim 𝒲 < dim 𝒩` by construction | 8 built rows at d = 4 | strictly between on all eight (so no invertible element; a sampled search agrees), criterion and spectrum agree |
| gate 12, N = 5 on two routes | 14 built rows, 1 strict row and 8 random rows at d = 32 | 10 pair with agreeing counts, 13 break with differing counts, the strict row among them with (dim 𝒩, dim 𝒲) = (192, 64); the float route's verdicts and nullities match the modular ones on all 23; the moved centre and the wrong shift call every pairing row broken |

**And one route whose construction shares no code with any of it.** The failure
this file could not otherwise rule out is a shared helper making the two sides
agree by construction. Gate 9 rebuilds the entire object in dense complex
numbers: its own Pauli matrices, its own Liouvillian, an eigensolver, and an
OPTIMAL MATCHING of the spectrum against its reflection. It makes **two**
comparisons and they are not the same one: its verdict against its own two
nullities, which is internal to the float route, and its two nullities against
the GF(p) eliminations row by row, which is the one that crosses. On 885 rows,
142 palindromic and 743 broken, both come back clean, **0** disagreements and
**0** nullity mismatches. The second of those is also direct evidence that the
mod-p nullities are not inflated on those rows, which is the fence two
paragraphs up answered from below rather than argued. Note what the gate does
import: `row_spaces`, and only as the other SIDE of the cross-check, since
without it there is nothing modular to compare against. Both of its thresholds
are defended as laws rather than chosen: the
spectral cut separates 5.2e−14 from 2.1e−3 (**10.6 decades**) and the
singular-value cut separates 1.2e−15 from 9.5e−4 (**11.9 decades**), and the
gates are on the separations.

That gate carries a warning in its own docstring, and it is here too because it
is the kind of error that arrives while checking for errors: **sorting** both
spectra and comparing them elementwise reports every row broken, the canonical
chain included. A lexicographic sort is not a matching: two spectra can be the
same multiset while sorting differently, because near-degenerate real parts break
the tie in opposite orders.

**Three caveats that belong here rather than in a footnote.**

The scored blocks are **samples, not full grids**: `grid3()` holds 21,952 letter
configurations and the four N=3 blocks each score the same 2,500 of them, so
"10,000 rows" is 2,500 configurations under four settings.

Every block carries both verdicts and the run prints the counts, but two carry
them very unevenly: the per-site mixed off-axis directions block holds on 1 of
its 192 rows, so it is effectively scoring the false-positive side alone, and
the nine F103 rows are all palindromic, so they test necessity only.

And the sharpest one: on most rows the criterion's equality is **0 versus
nonzero**, not a comparison of two positive integers. Rows with
`0 < dim 𝒲 < dim 𝒩`, the only rows where equality says more than "the
anticommutant is nonempty", are a minority: **212 of the 15,415 scored rows**,
concentrated where the algebra is reducible (145 on the disconnected graph, 28
on K₃, 20 on P₃ with equal ends, 6 and 5 at the single-letter XX and ZZ bonds, 5
with a Pauli string as the jump, 3 under non-uniform couplings and rates, and
none at all on the generic-magnitude or N=4 blocks). The distinction is exercised, and on every one of those rows the
criterion is right, but the headline row count overstates how much of the sweep
tests the equality rather than nonemptiness. **Gate 11 is the answer to that**:
eight rows are BUILT to sit strictly between, by taking H = H₁ ⊕ H₂ with
disjoint spectra where one block carries an invertible anticommuting element and
the other carries none, so 𝒲 is nonzero and entirely singular by construction.
The criterion says BROKEN on all eight and the spectrum agrees.

**Gate 12 is N = 5.** Through gate 7, N = 5 is reached by F103's three
rows, all three with agreeing counts and a pairing spectrum. Which rows can
fail is worth one sentence before the counts. On a row with dim 𝒲 = 0 nothing
can: the identity lies in 𝒩, so the counts differ, and L is unital, so 0 is an
eigenvalue and a pairing would need an eigenvector at −2σ, an element of 𝒲.
Such a row shows only that the two routes agree. A row can fail in two ways:
with agreeing counts, where a break would be a false positive, and with
0 < dim 𝒲 < dim 𝒩, where a pairing would be a false negative. Gate 12 has
both kinds at d = 32. Twenty-three rows, fourteen built, one strict and eight
random, each read on the modular route (the two nullities as 1024 − rank over
GF(p), the palindrome at two evaluation points) and on the float route (the
spectrum matched optimally against its reflection, the nullities from singular
values), the four readings compared row by row: ten rows pair with agreeing
counts, thirteen break with differing counts, and the two routes agree on
every verdict and every nullity. The built rows are the ones F138 names,
carried to N = 5: the canonical chain with and without a field, X and Y
fields, alternating signs, two dephasing axes, three axes on Heisenberg
(broken) and on the Ising bond (pairing), the field along the dephasing axis,
the mixed field of clause 2's boundary, and the two-letter term XX + XY. The
strict row is gate 11's construction at d = 32, H = (1 + Z₄)X₀ + (1 − Z₄)(Z₀ + 5)
with the one jump Z₀: two blocks of disjoint spectra split by Z₄, the first
carrying X₀, which anticommutes with the jump, and the second Z₀ + 5, which
commutes with it, so 𝒲 is nonzero and entirely singular, (dim 𝒩, dim 𝒲) =
(192, 64), the counts differ and the spectrum breaks; a pairing there would
need a Jordan chain at −2σ, so this is the one N = 5 row on which Lemma 2 is
exercised. The count is derived, not read: H = 2P₊X₀ + 2P₋(Z₀ + 5) with
P± = (1 ± Z₄)/2 has four eigenspaces of dimension 8, so its commutant has
dimension 4·64; Ad(Z₀) swaps the two X₀ eigenspaces, giving 64 to 𝒩 and 64 to
𝒲, and fixes each Z₀ eigenspace of the other block, giving 128 to 𝒩 and nothing
to 𝒲. Every other breaking row of the gate has dim 𝒲 = 0, and the random
rows of the seed are not strict. Two controls through the same door: the matching
with the centre moved to −1.1σ and the modular palindrome at the wrong shift
call all ten pairing rows broken. The thresholds are read as separations: the
matching distance separates 7.1e−13 from 3.8e−4 (**8.7 decades**), the
singular values 1.1e−14 from 3.8e−4 (**10.6 decades**), the smallest break
being the bisector control below, whose two field magnitudes differ by a
sixth, the largest zero singular value the strict row's.

The row that matters most there is **the bisector at N = 5**: one Z jump on
the middle site, an X field on site 0 and a Y field on site 4 of equal
magnitude. No single Pauli string lies in the far kernel (the 1,023
non-identity strings tried; the non-identity strings commuting with the
Heisenberg chain are X^⊗5, Y^⊗5 and Z^⊗5, and none of them commutes with both
fields; the same search returns X^⊗5 on the canonical chain with an X field,
its positive control), yet the two counts agree at 1 and the spectrum pairs.
The element is the one of F138's U family, the site reflection composed with
the half-turn about the bisector (X + Y)/√2, which carries the site-0 X field
onto the site-4 Y field and negates the middle Z, exhibited and checked
exactly (U scaled to Gaussian-integer entries, [U, H] = 0 and U Z₂ = −Z₂ U
compared with == 0; with the site-4 field doubled the commutator is nonzero):
the "if" of the theorem beyond any colouring, at N = 5 on the chain with four
undephased sites, outside the sweep of
[the colouring page](../../experiments/THE_PALINDROME_AS_A_COLOURING.md),
whose five-site rows with two undephased sites find no symmetry element on
the chain, 36 on the star and 60 on K₅, and which states the in-place rule on
the five-site ring: every dephased site must stay fixed. Its control with the magnitudes unequal
breaks, as at N = 3. And the same row with the Z jumps on
sites 1, 2 and 3 breaks as well, its far kernel zero on both routes; the
U-family element is not in it because the reflection carries Z₁ to −Z₃
(U Z₁ = −Z₃ U exactly, while U Z₂ = −Z₂ U still holds), so U negates the jump
set and not every jump in place, and 𝒲 asks for every jump in place, A_l W A_l
= −W for every l (Lemma 1), which is the condition MirrorWorld's
`EndCount.SymmetryElement` builds its element under: an involution fixing
every site a jump acts on, each jump going to minus itself. At equal rates on
the permuted jumps the two-sided Ad_U still commutes with L and is a copy; the
one-sided action is a reflection only where the permutation fixes every
dephased site. With the jumps on sites 1 and 3 alone the far kernel is zero as
well.

One flag the gate carries is a tautology and is gated as one, so that it cannot
be read as evidence: under Z-dephasing "the operator is not diagonal" holds for
every element of 𝒲, since A W A = −W already forces every diagonal entry to
vanish.

## (h) Scope, and what would falsify it

**Proved for:** any Hermitian H; any nonempty finite set of Hermitian jump operators
squaring to the identity, single-site letters, off-axis directions `n̂·σ⃗`,
multi-site Pauli strings and full depolarizing sites alike, and, by the same proof with A_lρA_l† for A_lρA_l, any unitary jumps whose dissipator is self-adjoint ([the clock mirror](PROOF_QUDIT_CLOCK_MIRROR.md) §5); any strictly
positive rate profile, uniform or not; any topology, which is not an independent
axis but a special case of "any Hermitian H"; any N; any finite dimension, with
the odd-d case true because both sides are false (§(f5)).

**Gated for:** d = 2^N with N ≤ 5, plus d = 3 and d = 5 for §(f5). N = 3
carries both directions on the letter grids, and N = 2 through gate 11's strict
rows; at N = 4 the 240 path and ring rows have no strict row, so only the
false-positive side is exercised there; at N = 5 gate 12 carries ten
rows with agreeing counts, where a break would be a false positive, and one
strict row 0 < dim 𝒲 < dim 𝒩, where a pairing would be a false negative, each
read on two routes (§(g)); its twelve other broken rows have dim 𝒲 = 0, where
neither side can fail. Everything outside that range rests on the proof
alone.

**Outside:** jump operators that are not unitary, where the question changes (F137 recentres the palindrome under T1 where it survives); unitary jumps whose dissipator is not self-adjoint, where only necessity is guaranteed and, where the counts agree, the palindrome holds exactly when replacing every jump by its adjoint keeps the spectrum (F172's lone clock jump); and rates that are zero or negative. Where the two counts differ
the multiset is broken and the criterion says nothing about the Jordan
structure, while §(f11) reads, wherever a far element exists, a part that still pairs: the system on the range of a
Hermitian far element of maximal rank and the two cross parts together. Where they agree, §(f8) gives the Jordan structure in full. Z dephasing with
correlated site fields is outside as well, not unitary: its jumps c_μ are built linearly
from the Z_k, c_μ = Σ_k b_μk Z_k (b real) with Kossakowski matrix a_kl = Σ_μ b_μk b_μl, and such a
jump squares to a multiple of 1 only when it is a multiple of a single Z_k. There F137's only
candidate centre, −σ with σ = tr a, fails whenever the positive semidefinite a has covariances
summing to a positive number and H keeps the all-0 and all-1 states as eigenstates: then no
single centre serves the spectrum. For an H that conserves the
magnetization and commutes with X^⊗N, the sector of coherences |i⟩⟨j| with ket and bra
popcounts (p, q) still pairs with the sector (q, N − p), about −(K_p + K_q)/2, wherever
K = Σ_μ c_μ†c_μ takes one value K_p on popcount p and one value K_q on popcount q, which every
sector does under one covariance common to every pair
([The Mirror under a Shared Field](../../experiments/THE_MIRROR_UNDER_A_SHARED_FIELD.md)).

**What would falsify it:** one configuration where the two nullities agree and
the spectrum does not pair, or one where it pairs and they differ. Either would
break a specific lemma, and the gate reports which side it fell on.
