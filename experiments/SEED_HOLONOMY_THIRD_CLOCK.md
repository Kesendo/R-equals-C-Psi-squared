# The seed holonomy is the third clock: generic EP2 geometry, not Π's Z₄

**Date:** 2026-07-16
**Status:** Closed. The one seam [the one-four thesis](ONE_FOUR_THESIS.md) left open is decided: the i⁴ = 1 eigenvector-frame holonomy around the defective seed is generic exceptional-point geometry, present identically in systems with no palindrome, chain, or Clifford structure. Gate: [seed_holonomy_generic.py](../simulations/seed_holonomy_generic.py) (24 checks); C# pins in `EigenvectorHolonomyTests` (the disguised Π-free EP2 + the Lᵀ = L precondition).

*Before the formulas: the one-four thesis ended with one four undecided. The frame of
eigenvectors circling our N = 9 seed turns 90° per loop and comes home after four, i⁴ = 1,
and nobody had shown whose four that is: the mirror's, the fermions', or its own. The
answer is: its own. Any two waves merging at an exceptional point turn this way, in any
system, with or without our mirrors. It is the third clock, and all three run on the
same i.*

## The question

[SeedHolonomyClaim](../compute/RCPsiSquared.Core/Symmetry/SeedHolonomyClaim.cs) records
that transporting the two coalescing right eigenvectors around the (1,2)-block defective
seed (in the biorthogonal vᵀv gauge natural to a complex-symmetric matrix) gives a frame
map M₁ with eigenvalues ±i: M₂ = −I, M₄ = +I, single-valued only after four loops. Until
today the claim marked its resemblance to Π's algebraic Z₄ as "a noted correspondence,
NOT a derived identity". [The one-four thesis](ONE_FOUR_THESIS.md) sharpened the
alternatives: Π's spectral Z₄, the Clifford-degree mod-4, or a third clock on the same i.

The discriminating test is the same move that refuted the strong one-four thesis: if a
system WITHOUT Π shows the identical holonomy, the holonomy cannot be Π's.

## The derivation (why every complex-symmetric EP2 turns this way)

Setup: a family M(q) of complex-symmetric matrices (Mᵀ = M) with two eigenvalues
coalescing at an order-2 exceptional point q* (a square-root branch point: the
eigenvalue gap scales as √|q − q*|). Near q*, write δ = q − q*. The two eigenvectors
behave as v± ≈ v₀ ± √δ · w + O(δ), where v₀ is the coalesced eigenvector, which at the
EP is SELF-ORTHOGONAL: v₀ᵀv₀ = 0 (the complex-symmetric EP signature). That kills the
O(1) term of the self-overlap, leaving vᵀ±v± ≈ ±c√δ → 0 with c = 2·v₀ᵀw. The vᵀv gauge
normalizes ṽ± = v±/√(vᵀ±v±).

One loop sends δ → δ·e^(2πi), so √δ → −√δ. Two things happen at once:

1. **The sheets swap** (the eigenvalue monodromy, order 2): v₊ and v₋ exchange roles.
2. **The normalizer winds a quarter-turn**: √(vᵀv) ~ (c√δ)^(1/2) ~ δ^(1/4) picks up
   e^(iπ/2) per loop.

Together, with the ±1 sign continuation fixed by continuity, the normalized frame maps
as ṽ₁ → ṽ₂, ṽ₂ → −ṽ₁: the matrix M₁ = [[0, 1], [−1, 0]] (a representative; which
off-diagonal carries the −1 is a column-ordering convention, the invariants are
eigenvalues ±i, det 1, M₁² = −I), order 4. This
is the Heiss-type state exchange (ψ₁, ψ₂) → (ψ₂, −ψ₁) of standard EP2 physics. Nothing
in the argument mentions Π, a chain, or a Majorana: the only inputs are the order-2
branch point and the complex-symmetric self-orthogonality that makes the vᵀv normalizer
degenerate. The mod-4 is the branch structure of the FOURTH root δ^(1/4), i.e. the same
scalar i, winding its own clock.

The gauge-contingency corollary is visible in the same formula: a Hermitian gauge
normalizes by vᴴv, which does NOT vanish at the EP, so no δ^(1/4) appears; only the
mod-2 sheet swap survives. The mod-4 is real but lives in the biorthogonal reading,
exactly as the claim has always stated.

## The from-below controls (the gate, 24 checks)

1. **Random Π-free pencils** (the decider): five random complex-symmetric 6-dim pencils
   A + qC, no palindrome, no chain, no Clifford structure. EPs located by minimizing the
   eigenvalue gap in the complex q-plane (residual gap ~1e-8; since gap ~ √|q − q*|, that
   is |q − q*| ~ 1e-16, machine-exact), order-2 pinned by the gap slope 0.500. All five:
   M₁ eigenvalues exactly ±i, M₂ = −I and M₄ = +I to machine zero, span residuals ~1e-15.
2. **A deterministic disguised EP2** (no optimizer): B(s) = blockdiag([[s, 1], [1, −s]], S₄)
   has its EP2 at s = i exactly; conjugating by a random complex-orthogonal Q (QᵀQ = I,
   complex Givens rotations) preserves complex symmetry and defectiveness while mixing the
   toy across all six dimensions, with spectators 1.1 away from the EP eigenvalue. Same
   signature, deviations 0.0. This construction is also the C# pin
   (`DisguisedPiFreeEP2_SameMod4Holonomy_TheThirdClockIsGeneric`), closing the gap the
   repo's own docstring flagged: the committed 2×2 toy was not discriminating (a bare
   2×2 has a preferred basis and no spectators).
3. **Gauge contingency**: the SAME loop transported in the Hermitian (continuity-phase)
   gauge shows M₂ ≈ +I (deviation 0.04, slow spectator phase drift) and never −I: mod-2
   only, as the derivation predicts.
4. **ε-stability**: ε ∈ {5e-4, 1e-3, 2e-3}, unchanged.
5. **The diabolic control** stays where it was: a diabolic crossing (independent
   eigenvectors, no √-branch) gives M₁ ≈ I, not ±i
   (`DiabolicCrossing_IsNotTheMod4Loop_FrameDoesNotRotate`).

Additionally the vᵀv gauge's precondition on OUR block, Lᵀ = L, established in prose in
[the β-exotic genericity note](F89_BETA_EXOTIC_GENERICITY.md) (real dephasing diagonal +
hop C = iK, K real symmetric), is now pinned from below at N = 5 and N = 9
(`WeightCoherenceBlock_IsComplexSymmetric_TheVTvGaugePrecondition`).

## What is decided, and what is not

**Decided:** the holonomy Z₄ is not Π's Z₄ and not the Clifford-degree mod-4. A Π-free
system carries it identically, so no derivation from Π is possible or needed: it is the
third clock, wound by the same one i (the fourth root δ^(1/4) at a square-root EP in the
self-orthogonal gauge). [SeedHolonomyClaim](../compute/RCPsiSquared.Core/Symmetry/SeedHolonomyClaim.cs)'s
correspondence node is amended accordingly; the one-four inventory now reads: one i,
THREE clocks (Π's spectral Z₄, the Clifford degree, the EP holonomy), pairwise distinct,
the first two wired at X^⊗N = Π², the third free-standing wherever two eigenvectors merge.

**Not decided (unchanged, correctly open elsewhere):** whether OUR seed is a generic
√-type EP2 at every odd N is the β-exotic item of
[the genericity note](F89_BETA_EXOTIC_GENERICITY.md) (settled N = 5 unconditionally,
N = 7 modulo the layer-identification premise; the second premise was discharged
2026-07-16). The third-clock statement is conditional on EP2-ness and
does not touch that question; the claim stays Tier 1 candidate for exactly that reason.

## Verification

    python simulations/seed_holonomy_generic.py
    dotnet test compute/RCPsiSquared.Core.Tests --filter "FullyQualifiedName~EigenvectorHolonomy"

24 Python checks (five random pencils × three certificates each, the disguised EP2 × four,
the Hermitian-gauge contingency, three ε values, the spectator guard); 5 C# tests (the
three prior pins + the disguised Π-free EP2 + the Lᵀ = L precondition).
