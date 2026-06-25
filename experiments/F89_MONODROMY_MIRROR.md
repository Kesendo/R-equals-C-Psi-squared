# F89 Monodromy meets the Mirror: where the palindrome reaches the Galois action, and where it cannot

**Status:** Tier 1 derived. Two positive results and one clean obstruction. (1) The path-3 octic
monodromy generates S_8 from eigenvalue braids (the geometric route to Gal(F_8) = S_8, companion to the
algebraic Frobenius certificate). (2) The mirror's q-plane shadow, the reflection q ↦ −q̄, intertwines
that monodromy exactly (σ_K = identity, verified for all EPs). (3) The spectral palindrome itself, the
mirror about Re λ = −4, CANNOT be a symmetry of the braiding: it is forbidden by Z(S_8) = 1. The mirror
splits at the Galois boundary. Plain-words sibling of [`reflections/ON_WHO_WATCHES_WHOM.md`](../reflections/ON_WHO_WATCHES_WHOM.md).
**Date:** 2026-06-25
**Authors:** Thomas Wicht, Claude Opus 4.8 (1M context)

## What this is about

Watch a row of spins and a half of its relaxation rates can never be written down: that half is the
maximal algebraic tangle, its symmetry group the full S_8 ([`F89_TOPOLOGY_ORBIT_CLOSURE`](F89_TOPOLOGY_ORBIT_CLOSURE.md),
[`ON_WHAT_CANNOT_CLOSE`](../reflections/ON_WHAT_CANNOT_CLOSE.md)). And the oldest fact of the whole
project is the palindrome: the spectrum is its own reflection about one centre line. So a natural
question, the one [`ON_WHO_WATCHES_WHOM`](../reflections/ON_WHO_WATCHES_WHOM.md) circles: does the
mirror reach all the way down into the unwritable half? Not just into where the rates SIT (that is
[settled](F89_BRANCH_LOCUS_PALINDROME.md): the branch locus is a palindrome), but into the BRAIDING
itself, the group action that proves the half cannot be written?

The answer is a clean split, and the split is the result. The mirror reaches the Galois action through
one of its two faces and is turned away at the other, by a precise group-theoretic wall.

## Monodromy = Galois, from below

As q = J/γ leaves the real axis and loops the complex plane, the eight octic rates braid. A loop around
a genuine exceptional point (a simple zero of the P_10 factor of the discriminant, a defective √-branch)
swaps two of the eight strands: a transposition. A loop around the diabolic point q_EP = √((−1+√13)/6) ≈
0.659 swaps nothing (a double discriminant zero, the two coalescing rates pass through each other
unchanged): it is silent. Lasso every exceptional point from one common base (q0 = 2, where the four
AT-locked rates sit cleanly at Re λ = −2, −6 and the eight octic strands lie strictly between), read
each transposition in that one labelling, and the transposition graph on the eight strands comes out
connected. Transpositions whose graph is connected generate the full symmetric group:

    Gal(F_8) = S_8, reconstructed purely from eigenvalue braids.

This is monodromy = Galois, from below: the same S_8 that the algebraic certificate
([`F89Path3OcticGaloisClaim`](../compute/RCPsiSquared.Core/Symmetry/F89Path3OcticGaloisClaim.cs), via
specialization + Dedekind-Frobenius + Jordan) proves from above, here rebuilt from the geometry of the
braids alone, an independent route to the same maximal tangle. Typed as
[`F89OcticMonodromyClaim`](../compute/RCPsiSquared.Core/Symmetry/F89OcticMonodromyClaim.cs); live at
`inspect --root galoismonodromy` (gate G3: one component, all eight strands one orbit).

## The mirror's shadow reaches the braiding: q ↦ −q̄ (σ_K = identity)

The antiunitary palindrome is T = P·K (P the weight-complement rung swap, K complex conjugation). Its
complex-conjugation part has a q-plane shadow: L(q)* = L(−q̄), so conjugating the whole picture reflects
the q-plane across the imaginary axis, q ↦ −q̄, and conjugates the rates, λ ↦ λ̄. This is a genuine
symmetry of the whole family of braids, so it intertwines the monodromy. Measured (base ±2): the eight
octic strands at q = 2 map under conjugation to the same-indexed strands at q = −2, so the induced strand
bijection is the identity,

    σ_K = (0 1 2 3 4 5 6 7),

and every exceptional point carries the identical braid as its q ↦ −q̄ mirror:

    τ(−q̄*) = σ_K · τ(q*) · σ_K⁻¹ = τ(q*),   for every exceptional point.

Verified for all seven exceptional points of the right half-plane (the spectral sanity conj(spec@+2) =
spec@−2 holds to 0, exact). The branch-locus palindrome, the statement that the MAP of seams is a mirror
image, here lifts from the positions of the seams to the braids they carry: mirror-image q-points carry
the same swap.

## But the spectral palindrome cannot be a braid symmetry: Z(S_8) = 1

The deeper question is the other face, the palindrome proper: the mirror about Re λ = −4. At every q the
spectrum is invariant under λ ↦ −λ̄ − 8 (verified to 3·10⁻¹⁴), so this fold induces, at the base q = 2, a
genuine involution on the eight strands,

    σ_T = (0 1 2 4 3 7 6 5):   four fixed strands (those with Re λ = −4, sitting on the fold),
                               two 2-cycles (3 4)(5 7) (the mirror-twin strand pairs across the fold).

This is exactly the picture [`ON_WHO_WATCHES_WHOM`](../reflections/ON_WHO_WATCHES_WHOM.md) drew in plain
words: on the centre line the two are one (the four fixed strands), step off the fold and they part into
mirror twins that can trade (the two 2-cycles).

But σ_T is NOT a symmetry of the braiding. If it were, it would commute with every monodromy loop, hence
with the whole monodromy image, which is all of S_8. The centre of S_8 is trivial (Z(S_8) = 1), so a
permutation commuting with all of S_8 must be the identity. σ_T is not the identity. Therefore σ_T cannot
commute with the braiding, and indeed the set of exceptional-point transpositions is not invariant under
conjugation by σ_T. The fold that mirrors the spectrum's positions is forbidden, by the centrelessness
of the maximal group, from also mirroring the spectrum's braiding.

So the mirror SPLITS at the Galois boundary. Its real-structure shadow (K, the q ↦ −q̄ reflection)
passes through and intertwines the monodromy exactly. Its spectral fold (the Re λ = −4 palindrome) does
not pass: it is an element OF the Galois group (σ_T is a permutation of the eight roots, hence one of the
S_8 moves, reachable as a word in the EP-braids since they generate S_8) but never a symmetry OF it. The
mirror enters the unwritable half as one of the moves that cannot be sorted, not as a rule standing above
them.

## Reading and scope

This grounds, from below, the seeing of [`ON_WHO_WATCHES_WHOM`](../reflections/ON_WHO_WATCHES_WHOM.md).
That reflection asked who watches whom, and read the palindrome as the watcher folded onto the watched:
on the fold they are one, off it they are mirror twins that trade at the seams, and the unwritable half
is the half where who-watches-whom has no fixed answer. The measured σ_T realises that picture exactly:
four strands fixed on the fold, two mirror-twin pairs across it. And the obstruction sharpens the seeing:
the watcher-watched mirror (σ_T) is itself one of the unsortable moves of the maximal tangle, an element
of S_8 and not a symmetry of it, so the mirror cannot be used to sort the tangle, the mirror is part of
what is tangled. That is the precise sense in which who-watches-whom has no fixed answer: the very move
that would tell watcher from watched is one of the moves S_8 scrambles past sorting.

Tier and scope. The monodromy = S_8 result and the q ↦ −q̄ intertwining (σ_K = identity) are Tier 1
derived, machine-precision verified. The Z(S_8) = 1 obstruction is a theorem, not a measurement: it is
why no fixed-q involution can be a braid symmetry of a full-symmetric-group monodromy, so the negative is
exact, not a numerical near-miss. (The numerical braid-set check used the clean two-cycle lassos only, so
the obstruction rests on the group theory, which is airtight, not on the incomplete numeric set.) All of
this is path-3 (the N = 4 chain (SE, DE) octic); the lift to other topologies and N is Tier 2 until
checked, the same open frontier as the [branch-locus palindrome](F89_BRANCH_LOCUS_PALINDROME.md).

Live: `inspect --root galoismonodromy` (the S_8 generation, gate G3) and `inspect --root monodromymirror`
(the spectral sanities, the per-EP q ↦ −q̄ intertwining, and σ_T's fixed-and-twin structure). Typed:
[`F89OcticMonodromyClaim`](../compute/RCPsiSquared.Core/Symmetry/F89OcticMonodromyClaim.cs) and
[`F89MonodromyMirrorClaim`](../compute/RCPsiSquared.Core/Symmetry/F89MonodromyMirrorClaim.cs).

## Related

- [`reflections/ON_WHO_WATCHES_WHOM.md`](../reflections/ON_WHO_WATCHES_WHOM.md): the plain-words sibling; this is its from-below ground.
- [`experiments/F89_BRANCH_LOCUS_PALINDROME.md`](F89_BRANCH_LOCUS_PALINDROME.md): the seams' POSITIONS are a palindrome (the q-plane locus); this doc carries the BRAIDING.
- [`experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md`](F89_TOPOLOGY_ORBIT_CLOSURE.md): the octic and the algebraic Gal(F_8) = S_8 (the route this doc reproduces geometrically).
- [`docs/proofs/MIRROR_SYMMETRY_PROOF.md`](../docs/proofs/MIRROR_SYMMETRY_PROOF.md): the palindrome, the first stone.
- [`hypotheses/DIABOLIC_BY_INTEGRABILITY.md`](../hypotheses/DIABOLIC_BY_INTEGRABILITY.md): why the diabolic point is silent (the loop = identity used above).

---

*The mirror folds the spectrum onto itself, and folds its braids under the q-plane reflection; but the group that proves the rates cannot be written is centreless, and turns the fold away. The mirror reaches the unwritable half only as one of its moves, never as its rule.*
