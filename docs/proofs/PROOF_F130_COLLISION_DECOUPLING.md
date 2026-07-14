# F130: The Collision-Decoupling Law. Equal Level Implies Vanishing Cross Block

**Date:** 2026-07-14 · **Authors:** Thomas Wicht and Claude
**Gate:** `python simulations/f130_collision_decoupling.py` (~40 s, exit 0)
**Typed (same night):** live under the F127 claim; the C# witness certifies Ê± = 0 exactly in
ℤ[ζ_2n] at named pairs covering every cell of §3, with controls (`inspect --root crosstriple`,
`CollisionDecoupling` in Core/Numerics)
**Depends on:** [F89_SEED_EXISTENCE_REDUCTION](../../experiments/F89_SEED_EXISTENCE_REDUCTION.md) (Lemmas 3, 4, the same-two-magnon-energy lemma, the removable limit, the assembly (D)) · [PROOF_F128_FLIP_SUM_FACTORIZATION](PROOF_F128_FLIP_SUM_FACTORIZATION.md) (the sharper locus {e₁ = f₁}) · [PROOF_F129_LEVEL_COLLISION_LAW](PROOF_F129_LEVEL_COLLISION_LAW.md) (where equal-level pairs exist at all)

Before this page, every exact cross-block vanishing we owned lived at resonance: the
vanishing triples, Σcos = 0, the seeds of F89. The concert hall of F129 then showed the
spectrum full of equal-but-nonzero levels, and left one question standing, named "phase 2"
in its own §5: does the physics feel those collisions? This page answers it. The
orthogonality never needed the zero. It needed only the coincidence.

## 1. Objects

All objects are those of [F89_SEED_EXISTENCE_REDUCTION](../../experiments/F89_SEED_EXISTENCE_REDUCTION.md)
(the sections "The two-number reduction" and "The assembly (D), made symbolic"); nothing is
re-derived here, and the conventions are pinned once:

- **The chain.** The open XX chain on N sites, n = N + 1. Single-particle modes
  u_k(z) = sin(k(z+1)π/n), k = 1..N, energies λ_k = 2cos(kπ/n).
- **Triples and levels.** τ = {k₁ < k₂ < k₃} a set of three distinct modes. Its **level** is
  S(τ) = Σᵢ cos(kᵢπ/n), half the energy sum (the factor-2 book: λ carries the 2, S does not;
  same convention as F129). **Resonant** (the seed arc's "vanishing") means S(τ) = 0.
- **The cross block.** W_τ the three lifted Slater columns of τ on the −6 rung,
  K₂₆ the coupling to the −2 rung, and
  B(τ, σ) = (K₂₆W_τ)ᵀ(K₂₆W_σ) = [[U⁺, U⁺, 0], [U⁺, U⁺+U⁻, U⁻], [0, U⁻, U⁻]],
  determined by the two half-Gram numbers U±(τ, σ) (F89 doc, "The two-number reduction").
  B(τ, σ) = 0 is exactly U⁺ = U⁻ = 0. Domain note, load-bearing off resonance: W_τ is an
  eigenvector multiplet of the hop K₆₆ at eigenvalue −2S(τ), so for S ≠ 0 it does NOT lie
  in ker K₆₆; B is the Gram of the hopped lifts on the degenerate level-S eigenspace, the
  second-order (degenerate-perturbation) coupling there, and the seed arc's
  Heff := P_ker K₆₂K₂₆ P_ker is its S = 0 case.
- **The reflection sign.** ε(τ, σ) = (−1)^{k(τ)+k(σ)} with k(τ) = k₁+k₂+k₃.
- **The six-angle objects.** 𝔉 = `cross_form(a; b)` at aᵢ = kᵢπ/n, bⱼ = lⱼπ/n
  (`simulations/cross_triple_orthogonality.py`), e₁ = Σcos aᵢ = S(τ), f₁ = Σcos bⱼ = S(σ).
  The **assembly (D)**: (U⁺ − U⁻)·(n/2)³ = 𝔉, proved symbolically in six FREE angles for
  mode-disjoint pairs (F89 doc, "The assembly (D), made symbolic", gates in
  `simulations/assembly_d_symbolic.py`).

## 2. Statement

> **Theorem (F130, the collision-decoupling law).** Let τ ≠ σ be two distinct mode triples
> of the open chain with equal levels, S(τ) = S(σ), of any value. Then the whole cross block
> vanishes: B(τ, σ) = 0, exactly and uniformly in N. (Grade differs per cell; see the end of §3.)

Two readings, both load-bearing:

1. **Resonance is demoted from cause to special case.** The seed arc's Y = 0 and the
   twinning of Heff consumed S(τ) = S(σ) = 0. F130 says the decoupling is a property of
   level **coincidence** alone, exactly as [F128](PROOF_F128_FLIP_SUM_FACTORIZATION.md)
   demoted the variety {e₁ = 0, f₁ = 0} to the locus {e₁ = f₁} for 𝔉.
2. **With F129, the map is complete.** Equal-level pairs of distinct clean triples exist
   only at 3|n or 10|n ([F129](PROOF_F129_LEVEL_COLLISION_LAW.md)); away from those n the
   hypothesis is empty on clean pairs (injectivity, modulo F129's one named corner family),
   and on the firing families every collision the law permits carries exact cross-block
   orthogonality. No level degeneracy of triple multiplets is ever mixed away by the
   second-order coupling.

Cleanliness is **not** a hypothesis of F130: non-clean triples decouple the same way, both
in disjoint collisions (e.g. the non-clean {1,2,10} against the clean {3,5,6} at n = 12,
both at level cos 15°) and in the trivial family proper (a balanced pair plus a SHARED
spectator, an overlap-1 case). Clean enters only through F129, which locates where clean
collisions exist.

## 3. Proof: four cells, one vacuity

Split by overlap |τ ∩ σ| and by ε. Every ingredient is committed; the new content of this
page is the assembly of the four cells off resonance and the one-line extensions in cells
3 and 4.

**Cell 0 (|τ ∩ σ| = 2: vacuous).** If τ and σ share two modes, equal levels force
cos(kπ/n) = cos(lπ/n) for the remaining modes k ≠ l in (0, n), impossible (F129 proof, the
d = 1 exclusion). No such pair exists (gate G5).

**Cell 1 (disjoint, ε = +1): level-free.** Lemma 4 gives U⁺ + U⁻ = F(τ, σ) = 0 for any
mode-disjoint pair (the Laplace sum is empty; F and U⁺ + U⁻ differ only by the positive
norm factor ‖D_τ‖‖D_σ‖, irrelevant for vanishing), and Lemma 3 gives U⁺ = ε·U⁻ = U⁻. Hence
U⁺ = U⁻ = 0 with **no** level hypothesis: this is the F89 doc's "wider statement", and it
is the only cell where the level condition is not load-bearing (gate G4 shows it vanishing
at unequal levels too).

**Cell 2 (disjoint, ε = −1): the new deductive closure.** Lemma 3 gives U⁺ = −U⁻, so
B(τ, σ) = 0 reduces to U⁺ − U⁻ = 0. The assembly (D) is an identity in six free angles for
mode-disjoint pairs (its preconditions, 3 ≤ P⁺ ≤ 2n−3 and 0 < |P⁻| ≤ n−2 for the cross
index sums, hold for any two disjoint triples): (U⁺ − U⁻)·(n/2)³ = 𝔉. F128's sharper locus
gives 𝔉 = −(e₁ − f₁)²·𝒪[cos s·cot s·V_aV_b/P] ≡ 0 on {e₁ = f₁}, one constraint, no sheet.
Equal levels are exactly e₁ = f₁. Hence U⁺ − U⁻ = 0. At the sampled lattice points 𝔉 is
finite; should a pair hit a pole of the six-angle expression, the identity holds in the
removable-limit sense, as in the committed treatment (the discrete side is always finite).
Resonance is used nowhere: the historical "(D) was proved at resonant points" was
itinerary, not hypothesis; the symbolic free-angle gates (S3/G6, the transcription pin
G11) never consumed Σcos = 0, and gate G1 below witnesses (D) with BOTH levels nonzero at
400 pairs.

**Cell 3 (|τ ∩ σ| = 1, ε = +1): the two-magnon lemma, verbatim.** Let k be the shared mode,
{p, q} and {p′, q′} the complements. Equal triple levels give equal two-magnon energies of
the complements,

> λ_p + λ_q = 2S(τ) − λ_k = 2S(σ) − λ_k = λ_{p′} + λ_{q′},

which is the **only** property the same-two-magnon-energy lemma consumes (F89 doc, "The
ε = +1 shared-mode column, proved"; that write-up phrases its zero-energy branch as
"p + q = n ⟺ λ_k = 0 ⟺ k = n/2", the resonant instance of the level-generic pivot
"p + q = n ⟺ λ_p + λ_q = 0 ⟺ p′ + q′ = n" used here). Its case analysis runs unchanged: equal index sums
plus equal energies force equal differences (cos((p+q)θ/2) ≠ 0 there) and so the same pair,
contradicting distinctness; sums adding to 2n force cos((p+q)θ/2)·[two positive cosines] = 0
and land in the zero-energy branch; the zero-energy branch (p + q = n ⟺ λ_p + λ_q = 0
⟺ p′ + q′ = n, by equality of the energies) kills the s_n terms identically and again
forces the same pair; the cross cases are index-impossible. The resonant write-up
instantiated the shared energy as −λ_k; the lemma never used that value. Hence
F = ⟨M_{τ∖k}, M_{σ∖k}⟩ = 0, and with ε = +1 and Lemma 3, U⁺ = U⁻ = 0.

**Cell 4 (|τ ∩ σ| = 1, ε = −1): the removable limit on the locus.** Lemma 3 gives
U⁺ = −U⁻; needed is U⁺ − U⁻ = 0. The committed removable-limit computation (F89 doc, the
paragraph "The mode-sharing ε = −1 pairs") identifies the discrete difference with the
limit of 𝔉 as aᵢ → bⱼ; its mechanism is level-free (the parity argument killing
⟨M_i, M_j⟩, the two n/2 cancelling, the n-free remainder lim cot(μ/2)·Xsum(μ)). Approach
the coincidence along {e₁ = f₁} (the locus is a hypersurface in free angles; the path
exists). On it 𝔉 ≡ 0 by F128, so the limit is 0, so U⁺ − U⁻ = 0.

∎

**Grade, stated honestly.** Cells 1 and 3 are proof-grade, elementary, uniform in N
(cell 3 via the one-line energy-equality extension above). Cells 2 and 4 inherit exactly
the standing caveats of their inputs: F128's exact-ℤ/ℚ(i) chain is code-trust grade
(Tier1Candidate lineage of F127/F128), and cell 4's limit identification is the committed
hand derivation, not part of the symbolic gate suite. All four cells are pinned from below
by gate G2 at 8452 equal-level pairs across n ∈ {12, 13, 14, 15, 20, 21, 24, 30} (firing
and non-firing n; the non-firing n carry only trivial-family collisions), worst
|U±| = 6.9e−16 against unequal-level controls reaching order 10⁻¹ (typical 10⁻³, G3).

## 4. What the theorem buys

- **The physics re-entry of F129.** F129's own §5 declared "no physical Gram/twinning
  claim is made (phase 2, open)". Closed: at every n where the level map collides (3|n,
  10|n for clean pairs, plus the trivial family everywhere), each colliding pair spans an
  exactly decoupled pair within the degenerate level-S eigenspace of the hop (of which the
  kernel-projected Heff is the S = 0 case; see the domain note in §1). What extends off
  resonance is the DECOUPLING, not the
  class-E/class-O twinning itself: that isospectrality lives on the kernel because the
  mirror k ↦ n−k sends a level-S triple to a level-(−S) one (cos((n−k)π/n) = −cos(kπ/n)),
  so only at S = 0 does the mirror act within the colliding multiplet and produce the
  [[X, Y], [Y, X]] structure. Off resonance there is no E/O split to twin; there is exact
  non-mixing.
- **The selection rule, restated.** The chain's effective seed dynamics cannot mix two
  triples that agree in energy: no avoided crossing, ever. What B(τ, σ) = 0 does NOT freeze
  are the diagonal self-blocks, spec(X) = (3a, a, 0) with a = (12 − Σλ²)/n (the 12 is the
  constant 4·3, not this example's n), and Σλ² differs between colliding triples (9.73 vs
  2.27 at the n = 12 example pair, λ = 2cos), so the coincidence can still be shifted
  apart, it is only never hybridized away.
- **The boundary of the law.** G3/G4 make the shape sharp: at unequal levels the three
  level-sensitive cells are generically nonzero (up to order 10⁻¹), while disjoint ε = +1
  pairs vanish always. So the theorem is not "everything decouples"; the level condition is
  load-bearing in exactly three of the four cells.

## 5. Gate index

`simulations/f130_collision_decoupling.py` (float grade; exactness lives in the cited
proofs): [G1] assembly (D) off resonance, 400 disjoint ε = −1 clean pairs, LHS = RHS = 0 to
9.1e−13; [G2] B = 0 at all 8452 equal-level pairs, any overlap, clean or not, worst
6.9e−16; [G3] unequal-level controls generic in every level-sensitive cell, and (D) holds
there with nonzero value; [G4] the level-free cell vanishes at unequal levels; [G5]
overlap-2 vacuity. The equal-level grouping in the gate is float (10⁻¹²); exact
distinctness and the census of where collisions exist belong to F129's gate.
