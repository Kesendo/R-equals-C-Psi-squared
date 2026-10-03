# PROOF: The Palindrome Is a Count of the Two Ends

**Registry:** [F158](../ANALYTICAL_FORMULAS.md).
**Status:** Tier 1 derived. The theorem holds for any Hermitian H, at least one
jump operator, all of them Hermitian and squaring to the identity, any strictly
positive rate profile, any topology, any N, and any finite dimension. What the
CRITERION decides is the eigenvalue multiset; sufficiency in fact delivers more,
a full similarity `L ~ −L − 2σ`, so where the criterion holds the entire Jordan
structure is reflection-symmetric (§(f8)). Two fences are load bearing rather
than decorative: `A² = 1` (F137 recentres the palindrome for non-unitary jumps)
and `γ_l > 0` with at least one jump present (a zero rate drops a condition and
can move the verdict; with no jump at all the equivalence holds trivially, and
§(f5)'s odd-d corollary, Lemma 3's 𝒩 ∩ 𝒲 = 0 and §(f9) fail). Even TOTAL dimension d is not a fence but a
consequence: with a jump present no palindrome exists at odd d (§(f5)), and the
theorem holds there with both sides false. It concerns the total dimension, not the
local one: on a qutrit tensored with a qubit, d = 6, the row H = 0, A = 1₃ ⊗ Z pairs,
with U = 1₃ ⊗ X. The rate fence is gated (10b) and the odd-d inequality behind the
corollary (10c); `A² = 1` is the boundary of the class and no gate crosses it. The gates
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

## What the repo already held, store by store

The sweep ran through two agents and a hand pass, and three adversarial
reviewers checked it against the stores at source.

- **[`docs/ANALYTICAL_FORMULAS.md`](../ANALYTICAL_FORMULAS.md).** F1 holds
  `Π·L·Π⁻¹ = −L − 2Σγ·I` with Π unitary and Π² = X^⊗N. **F4 holds one of the two
  counts**, in three separate statements that must not be fused: its Σγ = 0
  bullet says the kernel is the full commutant of H, while the Wedderburn form
  `dim ker = Σ mᵢ²` appears only in the one-seat, γ > 0, single-popcount-sector
  bullet for the algebra ⟨H_w, n_seat⟩. F138 holds the law under test, with two
  operators exhibited on its exceptions. F137 fences the premise rather than answering
  it: T1 jumps keep a palindrome and recentre it, and F137 states that centre
  **exactly**, as trace(L)/dim. Nothing anywhere on the two counts as a
  **criterion**.
- **[`docs/proofs/`](.), and the sharpest hit first.**
  [PROOF_CODIM1_BY_ADDITIVITY](PROOF_CODIM1_BY_ADDITIVITY.md) §6's **window-edge
  lemma** owns Lemma 2 below in a stronger form: *"If an eigenvalue of L sits on
  an edge of the block's rate window … then its eigenvector v is a joint
  eigenvector of A and B, and λ is semisimple."* That is the same
  numerical-range argument at **every** window edge rather than at the two
  extremes, and it also delivers most of Lemma 1's reverse inclusion. The arc
  `site_resolved_vacuum_block` carries a standing instruction about exactly this
  paragraph, *"Cite it rather than re-deriving it a third time"*, so Lemma 2 is
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
  `(X_k + Z_k)/√2` as breaking its own reading, which is a member of this file's
  class, so nothing here leans on it.
  [MIRROR_SYMMETRY_PROOF](MIRROR_SYMMETRY_PROOF.md) owns the first link, *"the
  palindrome then forces its partner −2Σγᵢ to be an eigenvalue too"*, from trace
  preservation. Its Scope states F138's clauses with *when*, not *exactly when*,
  and the converse failure and its counts beside them.
  [PROOF_F111_HARD_CELL_PURE_D_TEMPLATE](PROOF_F111_HARD_CELL_PURE_D_TEMPLATE.md)
  names the far end **the anti-steady eigenvalue** and pairs a mechanism with a
  **measured** witness (there is no live C# witness for it), and states as open
  *"why the pure-D templates are exactly the ones that lose the anti-steady
  mode"*: §(f4) translates that question rather than closing it.
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
  Lindblad book at γ = 0.05. Those are the same point (γ_there = 2γ_canonical,
  and both give −0.30 at N = 3), and this file works throughout in the canonical
  book, in which the far end is −2σ.
  [SYMMETRY_CENSUS](../../experiments/SYMMETRY_CENSUS.md) carries both operator
  conditions months old and reads them as invariance, which the pairing-condition
  page already corrects.
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
  everywhere as a **shift constant**: GammaFold adds `−2σ·Id` per cell and wears
  `e^(−2σt)` as a veil; Mirror's price 2Nγ is a fold's affine cocycle between
  two different blocks. The one genuine eigenvector at the far end is `Cat`'s
  `|0…0⟩⟨1…1|` at −2Nγ, and it lives at H = 0.
- **Returned nothing, or less than it seems.** `recovered/` returned nothing.
  `reflections/` did NOT: `ON_TWO_TIMES.md` carries the far end (*"Near λ = 0 the
  memory is long; at λ = −2Σγ it is shortest"*), while this
  file's own companion page cites `ON_THE_SOFT_BREAK.md` twice.
  `fw.Confirmations` and the C# `ConfirmationsRegistry` carry the F1 identity
  verbatim in `lebensader_skeleton_trace_decoupling`; what they do not carry is
  anything on the two-end count. And `simulations/framework/` builds no
  commutant, but `diagnostics/d_zero.py::stationary_modes` returns a kernel
  dimension and pins N+1, which is §(f2)'s near-end count computed and running,
  by the float route F157's registry entry warns against. And
  `hypotheses/THE_OTHER_SIDE.md` Q6 is **not** open. It is marked *"ANSWERED 2026-06-01 (Klein routing): the mapping is now
  built and verified bit-exactly (Q·L·Q⁻¹ = −L−2Σγ·I to ‖·‖ ≤ 10⁻¹¹,
  N=3,4,5)"*, which is the repo constructively exhibiting a palindromizer for a
  family of palindromic cases.

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

**Two fences, both load bearing, and one consequence.** The rate fence is gated (10b);
`A² = 1` is the boundary of the class; even total d follows (below).

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

**A_l² = 1** is where the class ends. F137 records that under T1 jumps, which
are not unitary, the palindrome survives about a different centre, which F137
states exactly as trace(L)/dim, so a criterion phrased about −2σ is answering a
different question there. What the hypothesis buys, positively, is broader than
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
the −2σ comes from. It needs U invertible and nothing more. Similar operators
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
leftmost point: the GLOSSARY row quoted above measures a generic Hermitian H
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
[PROOF_F111](PROOF_F111_HARD_CELL_PURE_D_TEMPLATE.md) asks why exactly the
pure-D templates lose the anti-steady mode; by Lemma 1 they lose it exactly when
no operator anticommutes with every jump while commuting with H, which
translates the question from the spectrum to the letters without answering it.
[F103](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) §7.12's three rows, soft with
non-bipartite basis-state graphs and therefore beyond any 2-colouring, are
decided correctly by the count at N = 3, 4 and 5. The operator the criterion
returns there is worth naming precisely: 𝒲 is **one-dimensional** and spanned by a single Pauli string,
Y^⊗N for XX+XZ and XX+XZ+ZX and X^⊗N for YY+YZ, and the X^⊗N of that last
family is F1's own Π². So this is an agreement test against a mechanism the repo
already owns and classifies bit-exactly, not a supply of something §7.12 was
missing.

**(f5) Every jump must be traceless, and at odd d none can be.** An invertible
U with U A_l U⁻¹ = −A_l makes A_l and −A_l similar, so the ±1 eigenvalue
multiplicities of each A_l must agree: **every jump operator is traceless**. That
is a necessary condition which is H-free, costs one trace, and is strictly
stronger than a parity count on d, because it bites at even d too:
A = diag(1, 1, 1, −1) at d = 4 kills the palindrome for every Hamiltonian, and
nothing in the gated scope exercises it, every jump family there (n̂·σ⃗, Pauli
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
  the other side of its cross-check, and compares
  both its verdict and its two nullities row by row against the modular route.
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

**All 103 gates pass.** The criterion is scored against the palindrome on
**15,415 rows**, in both directions, of which **2,596 hold and 12,819 break**:

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

Beside the scoring, and on rows where the palindrome is not consulted at all:

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

One flag the gate carries is a tautology and is gated as one, so that it cannot
be read as evidence: under Z-dephasing "the operator is not diagonal" holds for
every element of 𝒲, since A W A = −W already forces every diagonal entry to
vanish.

## (h) Scope, and what would falsify it

**Proved for:** any Hermitian H; any nonempty finite set of Hermitian jump operators
squaring to the identity, single-site letters, off-axis directions `n̂·σ⃗`,
multi-site Pauli strings and full depolarizing sites alike; any strictly
positive rate profile, uniform or not; any topology, which is not an independent
axis but a special case of "any Hermitian H"; any N; any finite dimension, with
the odd-d case true because both sides are false (§(f5)).

**Gated for:** d = 2^N with N ≤ 5, plus d = 3 and d = 5 for §(f5). N = 5 is
thin and deserves its caveat: it is reached only by F103's three rows, all of
them palindromic, so the false-positive direction is untested there. N ≤ 4
carries both directions. Everything outside that range rests on the proof
alone.

**Outside:** jump operators with A² ≠ 1, where F137 recentres the palindrome and
the question changes; and rates that are zero or negative. Where the two counts differ
the multiset is broken, and the criterion says nothing about the Jordan
structure there; where they agree, §(f8) gives it in full.

**What would falsify it:** one configuration where the two nullities agree and
the spectrum does not pair, or one where it pairs and they differ. Either would
break a specific lemma, and the gate reports which side it fell on.
