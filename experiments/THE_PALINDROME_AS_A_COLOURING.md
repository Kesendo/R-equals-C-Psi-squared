# The palindrome as a colouring

**Status:** Tier 1 where it restates [F158](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md) (a Pauli string in the far kernel carries the palindrome) and where it derives F138's clauses from the colouring rule, which is exact; the census counts rows at N = 3 on three graphs whose palindrome verdicts come from singular-value ranks separated by twelve decades, not from an exact rank certificate, gated by [`simulations/f138_palindrome_colouring.py`](../simulations/f138_palindrome_colouring.py) → [`f138_palindrome_colouring.txt`](../simulations/results/f138_palindrome_colouring.txt).
**Date:** 2026-09-29
**Authors:** Thomas Wicht, Claude (Anthropic, Opus 5.5)

F138 states the dephasing palindrome's boundary as two clauses: at most two dephasing axes per component, and a field orthogonal to all of them, with a proviso that a single-letter bond such as ZZ tolerates all three axes. It measures them sufficient and knows their converse is false. F158 decides the palindrome exactly through its far kernel

  𝒲 = {W : [H, W] = 0 and A W A = −W for every jump A},

which holds an invertible element exactly when the palindrome holds. This page looks for the simplest such element, a single Pauli string, and finds that looking for one is a **colouring of the graph**: every site takes one Pauli letter, its colour, from the letters its jumps light (those that anticommute with every dephasing letter of the site), and every bond and every field of H has to accept the colours it touches. F138's two clauses are the colourings with one colour per component, and its proviso is the freedom a one-letter bond leaves. The colouring explains most of F138's exceptions at a one-letter bond and none at a two-letter bond, and it shows where the proviso overreaches: on a site carrying two axes.

The word is used once already. `PalindromeSoftCertifier`'s `LinearSiteColoring` 2-colours the site graph: a per-site string K, Z on one sublattice, with K·H·K = −H, a string that **anticommutes** with H; F103 §7.12 treats that site-graph test as a proxy for a diagonal operator on the basis-state graph. The colouring here chooses a string that **commutes** with H and anticommutes with the jumps. The two are siblings, one for the Hamiltonian's sign flip and one for the dissipator's.

## What the repo already held

Two scouts went through the stores by name before this page was written, and two review rounds checked the record against its sources.

- **docs/proofs/.** [The two-end count](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md) (F158) decides the palindrome by dim ker L = dim ker(L + 2σ) and builds the palindromizer from an invertible element of 𝒲; its §(f4) finds 𝒲 "one-dimensional and spanned by a single Pauli string" on F103's three rows, and its sweep notes, reading the depolarizing experiment, that "the U for two-axis noise is the global product of the MISSING letter, so H may carry a field along that letter and no other". [MIRROR_SYMMETRY_PROOF](../docs/proofs/MIRROR_SYMMETRY_PROOF.md) closes its per-site Z, X, Z row with "U = Y⊗Y⊗Y commutes with the Heisenberg H and anticommutes with every jump", and already fences the proviso on one side: "Three axes stacked on a single site is depolarizing there, and breaks the Ising bond too." [The Π factorization](../docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md) writes the canonical palindromizer as R·D with R(ρ) = ρ·X^⊗N, a one-sided multiplication by a Pauli string. F103 §7.12 names the per-site hidden-Q routing and carries the 2-colouring above. F121 §6 bounds per-site product mirrors of the dissipator.
- **docs/ANALYTICAL_FORMULAS.md.** F138 with its clauses, its measured exception counts (22 / 104 / 0 at a two-letter bond on P₃ / a bond plus an isolated site / K₃, 776 / 520 / 732 at one letter, 78 at coincident field magnitudes), its two exhibited families, U = SWAP₀₂·Z₀Z₁Z₂ and the sums U = YYZ + ZXX, XXZ + ZYY; F158; F5, whose depolarizing witness is "a global Pauli string with no identity letter that commutes with H".
- **experiments/.** [The pairing condition](THE_PAIRING_CONDITION.md) reads F138's clauses as its criterion "read under a restriction", and reports that "an exhaustive search over monomial S in the Pauli-string basis, with the permutation part site-local, found none on the exception rows". It also names F118's division of labour, the one-sided multiplication carrying −2Σγ, as the general one. [The depolarizing palindrome](DEPOLARIZING_PALINDROME.md) holds the same-site reading of the axis count.
- **CAUGHT_ERRORS.** The entry that scoped the two-axis law: "with X and Z noise on every site the far end is Y^⊗N alone, so the palindrome holds exactly when H commutes with it".
- **The OpenArcs registry.** `f138_converse_failures`: every dressing of the named kind is a monomial matrix in the Pauli basis and the exception rows it examined return none; whether F158 should become F138's Proof anchor, dropping the clauses to a corollary, is left open. `f138_clause_two_sweep`: staging only.
- **The typed layer.** `PalindromeTwoEndCountClaim` and its witness; `PalindromeSoftCertifier` with `LinearSiteColoring`; the routing classes `TwoTermPalindromeRouting` and `KBodyPalindromeRouting`. Nothing on single strings in 𝒲.
- **hypotheses/.** THE_OTHER_SIDE: the uniform per-site routers are the dephase-letter palindromizers.
- **The glossary, MirrorWorld, reflections, fw.Confirmations.** Nothing on this.

Checked for adjacency. F5's commuting string and F158's string in 𝒲 are the same kind of object at the two edges of the rate window; F5 cites F158's zero-error end and F158 cites F5, but neither names the other's string. The per-site routers are never phrased as a string in 𝒲. F103's 2-colouring and this colouring are siblings that do not cite each other. No store states that F138's clauses are the one-colour solutions of a colouring, or counts how many exceptions a single string explains. That is what this page adds; the sufficiency it uses is F158's.

## The colouring

Take H a sum of distinct Pauli strings with nonzero coefficients (bonds P⊗P with P ∈ {X, Y, Z}, and one-site fields) and one Pauli letter per dephasing jump. Look for a single Pauli string F in 𝒲.

- **Anticommuting with the jumps** is per site: the letter of F at site l must anticommute with every dephasing letter of site l, i.e. be lit by all of them. At a site dephased along one letter two letters qualify, along two letters one, along all three none; at an undephased site any letter does, the identity included. These are the site's **colours**.
- **Commuting with H** is per term. A Pauli string T commutes with F or anticommutes with it, and [T, F] = 2·T·F in the second case. Distinct terms give distinct products T·F, so nothing cancels: [H, F] = 0 exactly when every term commutes with F.
- **A bond P⊗P** commutes with F when both or neither of its sites anticommute with P.
- **A field P at site l** commutes with F when the colour there is I or P.

So an admissible colouring is exactly a Pauli string in 𝒲, and by F158 it carries the palindrome: the one-sided multiplication ρ ↦ F·ρ gives L ~ −L† − 2σ, and F158 §(f8) gives L ~ −L − 2σ. The one-sided factor R of the canonical Π_Z = R·D is the colouring F = X^⊗N.

**What the rule forces, bond by bond.** A letter's three answers to "do you anticommute with X, with Y, with Z" always contain an even number of yeses (I: none; X: Y and Z; Y: X and Z; Z: X and Y), so any two of them name the letter.

- A **Heisenberg** bond, XX + YY + ZZ, and any **two-letter** bond, asks two of those questions of both ends and gets the same answers, so it forces the same colour on its two ends. A connected component then needs one colour, lit on every site: at most two dephasing axes in it, **F138's clause 1**.
- A **field** fixes the colour of its site to its own letter (or I), and the colour must be lit: with one colour per component the field lies along it, orthogonal to every dephasing axis, **F138's clause 2**.
- A **one-letter bond** P⊗P asks only one question, so it only puts both ends in the same class, {I, P} or the other two letters. The second class is open at every site dephased along a single axis, whichever the axis, so the colours can change from site to site and three axes coexist. That is **F138's two-term proviso**, derived: the colour freedom of a one-letter bond.

F138's clauses 1 and 2 are the colourings with one colour per component; the proviso is where a one-letter bond lets the colour change.

## The census

N = 3 on the path P₃, the triangle K₃ and a bond plus an isolated site; every nonempty set of bond letters from {XX, YY, ZZ}; F138's own domain of at most one dephasing axis per site; every pattern of dephasing and of fields (4⁶ = 4096 rows per graph and bond set), fields at F138's committed magnitudes (0.30, 0.22, 0.41), all positive, the grid F138's own counts were measured on. The palindrome is read from F158's count, the two nullities by singular values, separated by twelve decades.

| per bond set | P₃ | K₃ | bond + isolated site |
|---|---:|---:|---:|
| one-letter bond: palindromic rows the clauses reject | 776 | 732 | 520 |
| of these, with an admissible colouring | **714** | **714** | **468** |
| two-letter bond: palindromic rows the clauses reject | 22 | 0 | 104 |
| of these, with an admissible colouring | **0** | – | **0** |
| three-letter bond | 0 | 0 | 0 |

The three one-letter sets and the three two-letter sets give the same numbers each. The rejected-palindrome counts are F138's own, reproduced, which anchors the census on the registry. No row has an admissible colouring without the palindrome.

At a one-letter bond the colouring explains 92 % of F138's exceptions on P₃, 98 % on K₃ and 90 % on the bond with an isolated site: those palindromes have a single string as palindromizer, and the clauses were only too coarse to see it. At a two-letter bond it explains none, and cannot: a two-letter bond forces one colour per component just as Heisenberg does, so its colourings are F138's. Those exceptions, and the 62, 18 and 52 one-letter rows left over, have palindromizers built from an invertible element of 𝒲 that is no single string, like the sums YYZ + ZXX on F138's clause-1 rows.

**Fields on a three-axis component.** The clauses admit a three-axis component only without a field, since no direction is orthogonal to all three axes. The palindrome admits more: on P₃ with the ZZ bond and axes X, Z, Y, an X field on the Z-dephased middle site pairs, because X is that site's colour (the colouring is YXX), while a uniform X field and a Z field do not (gate stage D). The rows F138 quotes for this case, N = 3 under ZZ with axes X, Y, Z and a field along one letter on every site, pair 64 of 64 eigenvalues with no field, none under X or Y and 8 under Z; each field row is broken, and the colouring predicts all four, since a field along P on every site meets the site dephased along P, which has no colour P.

## Where the proviso overreaches

Read as written, the proviso lets a component whose bonds carry one letter hold three axes in any arrangement. Three axes on one site were already excluded: F138 calls the depolarizing channel clause 1 failing, and MIRROR_SYMMETRY_PROOF says the stacked site breaks the Ising bond too. A site with **two** axes was excluded nowhere, and there the proviso is wrong. Such a site has a single colour, the third letter. When that letter is the bond's own, a one-letter bond asks it (or I) of each neighbour, and a neighbour dephased along it has neither. Over every pattern with up to three axes per site, with the clauses read in that fair way, they accept 882 rows whose palindrome fails, summed over the three graphs and the seven bond sets; in every one a two-axis site sits in a bonded component with one bond letter and three axes, and none has an admissible colouring. F138 measured one axis per site, where its sufficiency never broke. The colouring carries the right scope by construction: a site's colours are the letters all its jumps light.

## Scope, and what would falsify it

The colouring rule and the derivation of F138's clauses hold for any graph and any N, for H a sum of distinct Pauli-string bonds P⊗P and one-site fields with single-letter dephasing jumps. The census is N = 3, three graphs, unit bond weights, one field-magnitude tuple, positive field signs. Coincident field magnitudes, F138's second failure, add palindromes by the graph's own symmetry (U = SWAP₀₂·Z₀Z₁Z₂, not a string), and are not in this census.

What would falsify it: a row with an admissible colouring and no palindrome (that would contradict F158), or a rejected-palindrome count differing from F138's recorded ones on the same grid.

## Open

- **The rows beyond the colouring.** 62 / 18 / 52 at a one-letter bond and every two-letter exception: their 𝒲 holds an invertible element but no single string. One of them, worked out exactly (gate stage E): P₃ with H = a(XX + YY) on the first bond and b(XX + YY) on the second, jumps X, Z, Y on the three sites. Both bonds force one colour and the three axes leave none, yet U = a·YYZ + b·ZXX commutes with H, anticommutes with every jump and squares to (a² + b²)·I, so it carries the palindrome at every positive rate; neither string commutes with H alone, their two commutators cancel, and at a = b = 1 it spans 𝒲. The equal-weight U is F138's; the weighted form was handed over by a second session. How a sum of strings does the exchange a single string does, and whether a coarser rule (pairs of letters, rotated local colours, the graph's own symmetry) catches these rows, is the next question.
- **F138's Proof field.** The colouring derives the clauses; whether F158 becomes their Proof anchor is the arc's open question, and the colouring is the corollary they would drop to.
- **The typed layer.** The colouring is a Python gate only; its C# home would be beside `PalindromeTwoEndCountWitness`, next to `LinearSiteColoring`.
