# PROOF F129 addendum: the family-inventory counts, derived

*2026-07-15, the same day the inventory landed. [experiments/F129_FAMILY_INVENTORY](../../experiments/F129_FAMILY_INVENTORY.md) recorded thirteen families with closed-form counts at verified grade; this document derives the counts. The derivation was carried out in four parallel strands (A+B, the d=2 trio, E+F, the constants), each independently verified against the exact census; residual code-trust steps are collected in §8, honestly. One discovery en route: the bundled family M splits as 40 + 0 + 60 over the three CDK order-210 types, and the middle type can NEVER fire (§6).*

## §1 The counting frame

Fix n, m = 2n, ζ = e^{iπ/n}. From the reduction of [PROOF_F129](PROOF_F129_LEVEL_COLLISION_LAW.md) §2: a disjoint collision pair {τ, σ} ↔ (E, A) where E = τ ⊎ (n−σ) is a 6-element multiset in (0, n), antipodal-free (no distinct pair e + e′ = n), n/2 at most once, and W = Σ_E (ζ^e + ζ^{−e}) = 0; the split A ⊂ E (a 3-sub-multiset) recovers τ = A, σ = n−(E∖A). An overlap-1 pair ↔ (E, A, s) with |E| = 4, A a 2-sub-multiset, and s the shared mode.

**Lemma 1 (automatic validity).** For antipodal-free E, every split yields τ ∩ σ = ∅ and τ ≠ σ: a shared or self-mirror element would force an antipodal pair in E. So a split fails only by REPEATED elements landing three-in-a-triple.

**Lemma 2 (split counts).** Generic E (all elements distinct): C(6,3) = 20 valid splits (d=3) or C(4,2) = 6 (d=2). Exactly one element doubled: the two copies must split across A and E∖A, leaving 6 (d=3). Fully doubled E = 2·{a,b,c}: exactly 1 split, the self-mirror pair σ = n−τ, which sits at level 0 (S = −S). Height 1 of all minimal pieces (the CDK import of PROOF_F129 §4) bounds multiplicities by 2 and makes supp(E) = |E| for single-piece and zero-mode-plus-piece families.

**Lemma 3 (the flip).** E ↦ n−E is exponent-shift by n (W ↦ −W); it is the second description of the same unordered pair (it swaps the roles of τ and σ). E is flip-symmetric only if fully-doubled degenerates permit it; for the families below the flip is fixed-point-free on valid configurations (checked per family), so

  #pairs = Σ over flip-orbits of (valid splits) [d=3],  #pairs = Σ over flip-orbits of 6·(#valid s) [d=2].

**Lemma 4 (the s-count, d=2).** For a fixed E and split, the shared mode s ranges over {1..n−1} minus the exclusion set X = supp(E) ∪ supp(n−E) (distinctness excludes the four private modes = one description's support; cleanliness of both triples excludes their antipodes = the other support), and X is split-independent. |supp(E)| = 4 (height 1); the two supports intersect exactly in n/2 when n/2 ∈ E (the zero-mode families) and are disjoint otherwise. Hence

  #valid s = n − 1 − |X| = **n − 8** (zero-mode families B, G),  **n − 9** (n/2-free families D, H).

Every surviving s gives a valid pair, distinct (split, s) give distinct pairs, and the pair's own reduction returns E; so per flip-orbit the pair count is exactly 6·(#s).

## §2 Families A and B: the R₃ label arithmetic

An R₃ piece is a coset {c, c+m/3, c+2m/3}, label c mod q, q = m/3 = 2n/3; conjugation is c ↦ −c, the flip shifts labels by q/2.

**Family A** (W = two conjugate R₃-coset-pairs {±c₁}, {±c₂}, labels allowed to coincide). Write t = n/3, so q = 2t and the flip on label classes is c ↦ t−c. Exclusions: real roots c ∈ {0, t}; self-antipodal class 2c ≡ t (solvable only for t even, at c = t/2). Good classes: G = 2(k+1) in BOTH parities of t, with k = ⌊(n−9)/6⌋; this parity collapse IS the floor in the formula. Cross exclusion between the two classes: exactly c₁ + c₂ = t, and these are precisely the flip-FIXED configurations, so removing them leaves the flip free. Count: doubled configurations (c₁ = c₂, E fully doubled, the level-0 mirror pairs) = G, in k+1 flip-orbits of 1 split each; generic = C(G,2) − (k+1) = 2k(k+1), in k(k+1) flip-orbits of 20. Total

  **A(n) = (k+1)·1 + k(k+1)·20 = (20k+1)(k+1).**

**Family B** (W = zero mode + one conjugate R₃-pair; d = 2). The zero mode occupies exponents {m/4, 3m/4}, whose R₃ labels are exactly the quarter-labels {q/4, 3q/4}; those are simultaneously the self-antipodal labels, so the exclusion sets coincide and the valid labels are q − 4 = 4(k+1), k = (n−12)/6; conjugation halves to 2(k+1) W's, the flip (fixed only at the excluded quarter-labels) halves to **k+1 orbits**. With Lemma 4:

  **B(n) = (k+1) · 6 · (n−8) = 12(3k+2)(k+1).**

The A and B strands verified EXACT SET EQUALITY of the generated pairs against the census (A: n ≤ 105; B: n ≤ 90).

## §3 Family C: the pentagon (the template)

W = zero mode + R₅-conjugate-coset-pair, label c mod q₅ = m/5 = 4j at n = 10j. Exclusions c ∈ {0, q₅/2} (real roots) and 2c ≡ q₅/2 (the antipodal condition, coinciding with the coset containing n/2): q₅ − 4 valid labels, (q₅−4)/2 conjugate classes, flip c ↦ q₅/2 − c fixed-point-free, all 20 splits valid:

  **C(n) = (q₅−4)/2 · 20 / 2 = 2(n−10).**

## §4 Families D, G, H: the weight-8 sporadics

W is a single rotated minimal sum (D: (R₅:3R₃), order 30; H: (R₇:R₃), order 42) or zero mode + (R₅:R₃) (G, order 30). A minimal sum of order M sits in μ_m as k·Q + b (k = m/M); conjugation-symmetry of the whole W forces 2b ≡ k·c₀ (mod m), where −Q = Q + c₀. Each pattern has trivial rotational stabilizer and even c₀, so the congruence has EXACTLY 2 solutions at every n on the door, differing by n: one flip-orbit per rotation class, at every n. Rotation classes are necklaces of expanded positions on the R_p skeleton: (R₅:3R₃) = 2 marks on a 5-cycle = 2 necklaces (gap 1 or 2) → **2 orbits**; (R₇:R₃) = 1 mark → **1 orbit**; (R₅:R₃) = 1 mark → **1 orbit**. With Lemma 4:

  **D(n) = 2·6·(n−9) = 12(n−9),  H(n) = 6(n−9),  G(n) = 6(n−8).**

## §5 Families E and F: the order-30 piece and the parity splits

Write u = n/15 = j, and P(s) = the rotation by ζ^s of the weight-6 order-30 piece; its label s is rigid (trivial stabilizer), conj P(s) = P(−s), so P(s) is self-conjugate ⟺ s ∈ {0, m/2}. All internal pattern differences are multiples of u: labels off the u-lattice are automatically valid ("generic").

**Family E** (W = R₃-conjugate-pair + one self-conjugate P; canonically s = 0, the flip identifying s = 0 ↔ m/2). R₃ label c ∈ (0, 5u) after conjugation; endpoint real roots exclude c ∈ {0, 5u}; the multiples of u are the special labels: c ∈ {2u, 4u} hit the R₅ shell and produce ONE doubled element (**2 six-split orbits, every n**), c ∈ {u, 3u} produce a forbidden antipodal pair (invalid), and c = 5u/2 (solvable only for u EVEN) makes the R₃ pair self-antipodal (invalid). Generic orbits: 5u − 5 (u odd), 5u − 6 (u even); the even-n deficit of the formula is exactly the one extra invalid label:

  **E(n) = 20·(5u−5) + 12 = (20n−264)/3** (n odd), **20·(5u−6) + 12 = (20n−324)/3** (n even).

**Family F** (W = P(s₁) + P(s₂), conjugation-closed). Both self-conjugate: {0,0} ~ {m/2,m/2} under the flip = **the 1 fully-doubled orbit** (the level-0 self-mirror member, its two pieces individually self-conjugate: the inventory's name "(R₅:R₃)+conjugate" is loose exactly there); {0, m/2} carries an antipodal pair, invalid. Conjugate pairs s₂ = −s₁: the label group ⟨−s, m/2−s⟩ has order 4, generic labels (off the u-lattice for u odd, off the u/2-lattice for u even) sit in free orbits: generic = 30(u−1)/4 resp. 30(u−2)/4 = **15⌊(j−1)/2⌋** in both parities. At u even, exactly four special labels s = (u/2)·{1,7,11,13} produce a single clean overlap (one doubled element): **4 six-split orbits, even n only**. One further special label, s = m/12 (u even), reassembles into zero mode + R₅-pair: it belongs to family C, and the first-fit classifier files it there (the "impostor"; it never enters F's tally).

  **F(n) = 1 + 20·15·(j−1)/2 = 10n − 149** (n odd), **1 + 24 + 20·15·(j−2)/2 = 10n − 275** (n even).

## §6 The constants I, J, K, L, M: the fixed-vertex mechanism, and M = 40 + 0 + 60

A single minimal piece on an ODD skeleton R_P (P = 7, 11): conjugation must act as a reflection of the regular P-gon, and every reflection axis of an odd gon passes through EXACTLY ONE vertex (F1). A bare term at the fixed vertex would sit at exponent 0 or m/2, a real root; so the fixed vertex always carries a fan (F2). All remaining freedom is the assignment of substitution branches to the (P−1)/2 conjugate vertex-pairs, a finite count with no n in it; the two admissible axis positions (exponent 0 vs m/2) are the flip pair. Hence #orbits = #branch-assignments, constant on the door:

- **I** (R₇:5R₃): the 2 bare vertices must form one conjugate pair, 3 choices → **3 orbits, 60 pairs**.
- **J** (zero + (R₇:3R₃)): the expanded set = fixed vertex + one pair, 3 choices → **3 orbits, 60**.
- **K** (R₁₁:R₃): the single fan is forced onto the fixed vertex → **1 orbit, 20**.
- **L** (zero + (R₇:R₅)): forced → **1 orbit, 20** (the corner-closure's second mechanism).
- **M**, the three CDK order-210 types: **(R₇:(R₅:2R₃))**: the 2 further-fanned R₅-subterms must be a conjugation-closed pair, {1,4} or {2,3}, 2 choices → 2 orbits, **40 pairs**. **(R₇:2R₃,R₅)**: the R₅ takes the fixed vertex, the two R₃'s one pair, 3 choices → 3 orbits, **60 pairs**. **(R₇:R₃,(R₅:R₃))**: TWO branches of distinct types would each need a fixed vertex, but there is only one (F1) → **0 orbits, impossible at every n** (verified: 0 of 5040 pieces at n = 105 are conjugation-symmetric).

  **M = 40 + 0 + 60 = 100**, and the bundle caveat of the inventory resolves: only two of the three CDK types ever fire.

## §7 Assembly

Every family formula of [experiments/F129_FAMILY_INVENTORY](../../experiments/F129_FAMILY_INVENTORY.md) is reproduced: A from two free coset labels (§2); B from one free label times the shared-mode freedom (§2); C, E, F from one free label line (§§3, 5); D, G, H from zero free labels times the shared-mode freedom (§4); the constants from the fixed-vertex assignment count (§6). The degrees-count-freedom reading of the inventory is now DERIVED, in the precise form: the polynomial degree equals the number of free coset labels PLUS one for the shared mode s of a d = 2 family (each free label and the s-line contribute one factor linear in n; a rigid assignment contributes a constant).

## §8 What remains code-trust (honesty)

Each strand verified its derivation against the exact census with independent generators (A/B: exact SET equality of pairs, n ≤ 105/90; D/G/H: all doors to 105; E/F: all firing n ≤ 150; constants: all doors incl. 210). The steps that are verified-exhaustively rather than proved uniformly in n:

1. **Height 1** of every piece (no repeated term) is the CDK import, as in PROOF_F129 §4; not re-proved.
2. **Necklace off-lattice rigidity** (§4): that no ℤ/m shift identifies expansions across necklace classes is checked exhaustively per pattern, not argued abstractly; likewise "no exotic embedding outside the scaled μ_M coset" leans on the census.
3. **The E/F special-label assignments** (§5: {2u,4u} doubled vs {u,3u} invalid; the F set (u/2)·{1,7,11,13}; the single impostor s = m/12) are congruence-derived and verified to n = 150, but "no further coincidence at large u" is not written as a uniform lemma.
4. **The constants' fingerprint constancy** (§6) is structural (F1/F2 + branch counting) with the "exactly 2 axis rotations" step confirmed by enumeration; a standalone rotation-counting lemma was not written.
5. **Lemma 3's flip fixed-point-freeness** is checked per family (A: c = t/2 excluded; B: quarter-labels excluded; the d=2 trio: E never flip-symmetric; E/F: the special labels absorb the fixed points), not proved by one uniform argument; **Lemma 4's** sufficiency legs (every surviving s valid, injectivity, reduction returns E) are proved from antipodal-freeness in the D/G/H strand and re-checked per family.
6. The strands' verification code lives in the session scratchpad (WIP policy); the committed, re-runnable verification is the census gate [f129_family_inventory.py](../../simulations/f129_family_inventory.py): its I2 asserts the thirteen totals to n = 210, and its I5 (added with this document) reconstructs the M sub-classification from committed substitution recipes. Beyond the totals and I5, the orbit/split factorizations quoted here are re-derivable from the census (the pairs are exact) but not pinned by a committed script.

## §9 Verification index

| claim | committed check |
|---|---|
| all thirteen totals, n ≤ 140 + capstones 150/210 | `f129_family_inventory.py` gates I1-I4 |
| the orbit/split factorizations quoted here | re-derivable from the census (the pairs are exact; group by the flip-quotiented E-multiset) |
| M = 40 + 0 + 60 at every 105\|n in range | gate I5: the census M W-sets equal EXACTLY the admissible rotations built from the committed substitution recipes (4 M1 + 6 M3 sets, completeness asserted loudly) |
| the closed-form SUMS (total + d-split) vs the exact census, live | the typed surface (2026-07-15): `CollisionFamilyInventory` (Core/Numerics) on the crosstriple witness, every n ≤ 60 + the n = 70 capstone, exact ℤ[ζ_2n]; the n = 105 sum 8858 pinned in `--filter CrossTripleOrthogonality` |

The counting derivation closes the "derivation open" item of the inventory; the M sub-classification and the impostor are new content beyond it.
