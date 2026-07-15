# The F129 family inventory: every collision counted, family by family

*2026-07-15. The morning after the corner closed. [PROOF_F129_LEVEL_COLLISION_LAW](../docs/proofs/PROOF_F129_LEVEL_COLLISION_LAW.md) told us WHERE collisions live (3|n or 10|n, nowhere else, unconditionally); this inventory says HOW MANY, and the answer is thirteen families, each with an exact closed form. The family LIST is complete at every n (a one-paragraph consequence of the corner-closure import, below); the COUNTS are VERIFIED on every firing n ≤ 140 plus the capstones n = 150 and n = 210, NOT derived (the ROT3-multiplicity precedent). Gate: [f129_family_inventory.py](../simulations/f129_family_inventory.py).*

Before the machinery: what this feels like. Yesterday's law said the concert hall has exactly two doors. Today we stood in the doorways and counted who walks through each, and the count is not a mess: it is thirteen clean families, and the SHAPE of each family's count function tells you the shape of its key. Families built from freely rotatable pieces grow quadratically or linearly with the hall; families built from one rigid sporadic piece are CONSTANT, the same 60 or 20 or 100 pairs at every tested n their door opens.

## The classification rule

Every collision pair (τ, σ) reduces (proof §2) to a vanishing sum W of weight 8 (shared mode, d = 2) or 12 (disjoint, d = 3) over μ_2n. W splits into minimal vanishing pieces, and every piece of weight ≤ 12 is on the complete Poonen-Rubinstein/CDK list (arXiv:2008.11268, Table 1): the same import that closed the corner. The family of a pair = the multiset of (weight, ratio order) of its pieces under a DETERMINISTIC decomposition: first-fit, with minimal pieces ordered by size then by subset index (an artifact of the fixed root enumeration order), committing to the first piece whose remainder is still tileable and backtracking only on failure. A W admitting several tilings is filed under the one first-fit finds; the rule is deterministic but not canonical. Corollary for the honest grade: the closed forms were FIT to this rule's output, so gate I2 verifies that the counts follow the forms, not that a label means what its name says. The asymmetry to hold on to: the PAIRS are exact (byte-equal ℤ[ζ_2n] level vectors, the committed census machinery); only the LABELS ride the mod-p layer (piece vanishing tested mod two independent primes ≡ 1 (mod 2n), false-positive odds ~1/p², no false negatives).

## The inventory

k = ⌊(n−9)/6⌋ for A, k = (n−12)/6 for B. "Door" = the divisibility that admits the family (each piece's ratio order must divide 2n).

| | pieces | d | door | count(n) | shape |
|---|---|---|---|---|---|
| A | R₃ + R₃ + R₃ + R₃ | 3 | 3\|n | (20k+1)(k+1) | quadratic |
| B | zero + R₃ + R₃ | 2 | 6\|n | 12(3k+2)(k+1) | quadratic |
| C | zero + R₅-conjugate-pair | 3 | 10\|n | 2(n−10) | linear (the pentagon) |
| D | (R₅:3R₃), weight 8 | 2 | 15\|n | 12(n−9) | linear |
| E | R₃ + R₃ + (R₅:R₃) | 3 | 15\|n | (20n−264)/3 odd n, (20n−324)/3 even | linear, parity-split |
| F | (R₅:R₃) + conjugate | 3 | 15\|n | 10n−149 odd n, 10n−275 even | linear, parity-split |
| G | zero + (R₅:R₃) | 2 | 30\|n | 6(n−8) | linear |
| H | (R₇:R₃), weight 8 | 2 | 21\|n | 6(n−9) | linear |
| I | (R₇:5R₃), weight 12 | 3 | 21\|n | 60 | constant |
| J | zero + (R₇:3R₃) | 3 | 42\|n | 60 | constant |
| K | (R₁₁:R₃), weight 12 | 3 | 33\|n | 20 | constant |
| L | zero + (R₇:R₅) | 3 | 70\|n | 20 | constant |
| M | weight-12, order 210 | 3 | 105\|n | 100 | constant |

## Why exactly thirteen

The LIST (not the counts) is forced, at every n, by what the corner closure already imported. W has weight 8 or 12, at most one zero mode (proof §2), no real root, and is conjugation-symmetric; every minimal piece of weight ≤ 12 is on the complete CDK list. An odd-weight piece can never be self-conjugate (conjugation on an odd multiset has a fixed point, a real root), so odd pieces, the prime cycles AND the weight-9/11 sporadics alike, come in conjugate pairs of equal weight: 7+7, 9+9, 11+11 all overshoot, and a 5+7 pairing is impossible (conjugates have equal weight). R₂ beyond the zero mode is an antipodal pair, excluded by cleanliness; and no weight-4 minimal sum exists (Table 1 has no weight-4 row), so {4,4} and {8,4} are empty slots. That leaves, as piece-weight partitions: for weight 8, {8}, {2,3,3}, {2,6}; for weight 12, {12}, {2,10}, {6,6}, {3,3,6}, {3,3,3,3}, {2,5,5}. Fanning each slot out by the admissible CDK orders (weight 6: order 30 only; weight 8: 30 or 42; weight 10: 42 or 70; weight 12: 42, 66, or 210) gives exactly the thirteen signatures of the table, and no fourteenth can ever fire, at any n. What stays verified-not-derived is only how MANY pairs each signature collects.

Readings, from the table:

- **The degree counts the freedom.** Two independently rotatable R₃-conjugate-pairs → quadratic (A, B). One free rotation against a fixed partner → linear (C through H). A single rigid weight-12 sporadic, pinned by conjugation symmetry up to a finite orbit → constant (I, K, M); the same for a rigid weight-10 sporadic chained to the zero mode (J, L).
- **The onsets are IN the formulas.** A(6) = 0, B(6) = 0, C(10) = 0: the F129 thresholds (3|n from 9, 10|n from 20) are not extra conditions, they are zeros of the inventory. (Narrative: those n sit outside the gated range.)
- **The sub-law, refined.** Every d = 2 family (B, D, G, H) has a door divisible by 3: the F129 sub-law (overlap-1 forces 3|n), now visible piece by piece: a weight-8 W is zero + R₃-pair (B), zero + (R₅:R₃) (G), or a single weight-8 sporadic (D, H), and every route carries the 3: the order-30 and order-42 pieces by their orders, the R₃-pair by itself.
- **The corner-closure's second mechanism is real.** Family L, zero mode + (R₇:R₅), is exactly the one assembly the corner closure singled out as possible at 70|n. It walks: 20 pairs at n = 70, 20 at n = 210.
- **Odd-n overlap-1 needs a sporadic.** At odd n there is no zero mode, so a weight-8 W must be a single piece: that is why ov1 collisions at odd n exist exactly at 15|n (D) and 21|n (H) and nowhere else; the empirical zeros at n = 9, 27, 33, 39, 51, ... are the absence of a weight-8 sporadic whose door divides them.

## Verification

The committed gate [f129_family_inventory.py](../simulations/f129_family_inventory.py) re-derives the whole table from scratch: exact collision pairs from the F129 census machinery, greedy decomposition per pair, then asserts (I1) every pair decomposes, (I2) every family count equals its closed form, (I3) the thirteen families partition the census, (I4) a family appears iff its door divides n. Default run: every firing n ≤ 140 (~15 min). `--deep` adds the capstones n = 150 and n = 210, where seven and twelve families fire simultaneously; both pass, with four of the five constant families co-firing at 210 (K needs 33|n). `--fast` runs n ≤ 66 (~1 min); the n = 210 point dominates `--deep` (~40 min total). The n = 210 capstone total (72,269 pairs) and every sub-count, and later the n = 132 and n = 140 points, were PREDICTED from the smaller-n fits before their runs (this file's authoring session; every pre-registered number hit).

The MirrorWorld run mode `collision N` (`LevelCollision.cs`) prints the same census totals (pairs/disj/ov1) live and independently over GF(p); its totals for n ≤ 60 match the family sums row for row (checked in this file's authoring session).

## What is NOT claimed

- **No derivation of the counts.** The closed forms are exact fits verified on 55 firing n including the two capstones; the counting argument (enumerate conjugation-symmetric, antipodal-free piece assemblies × valid A/B splits per W) is set up by the proof's §2 reduction but not carried out. That derivation is the natural next step; until then the counts are verified-grade, like the ROT3 multiplicity before it. The degree-counts-freedom reading above is an observation on the fits, not a theorem.
- **The constants are thin where the doors are rare.** L is tested at three n only (70, 140, 210) and M at TWO (105, 210); "constant" is the shape of two or three points plus the rigidity reading, not a dense verification. K has four points (33, 66, 99, 132), I seven, J four.
- **Family M is not sub-classified, and no piece NAME is re-verified.** The classifier records (weight, ratio order) only; the CDK type names in the table are annotations by that signature. M = 100 may bundle CDK's three weight-12 order-210 types; and family F's name is loose at its one level-0 member per firing n, whose two weight-6 pieces are SELF-conjugate rather than a conjugate pair.
- **Level-0 pairs are included.** F129 asks for equal levels; equal-to-zero qualifies, so the F127 dark-fringe world is inside these counts (at n = 30: 28 of family A's 244, 4 of D's 252, 12 of E's 92, and F's single degenerate member; recomputed in this file's authoring session).
- **The parity splits in E and F are recorded, not explained** (the even-n deficits −20 and −126 relative to the odd-n lines await the counting argument).
- **The labels depend on the first-fit rule.** A pair whose W admits several minimal decompositions is counted once, under the first-fit tiling. Any derivation must match THIS rule (or replace it with a canonical invariant and re-verify).

## Links

The law: [PROOF_F129_LEVEL_COLLISION_LAW](../docs/proofs/PROOF_F129_LEVEL_COLLISION_LAW.md) (§5 names this inventory as the open follow-up; closed at verified grade by this file). Registry: [ANALYTICAL_FORMULAS F129](../docs/ANALYTICAL_FORMULAS.md). The piece list: arXiv:2008.11268, Table 1 (the corner-closure import). The three-cosine ancestor: the ROT3/PENT multiplicities of the seed arc (`Seed.cs`, [F89_SEED_EXISTENCE_REDUCTION](F89_SEED_EXISTENCE_REDUCTION.md)).
