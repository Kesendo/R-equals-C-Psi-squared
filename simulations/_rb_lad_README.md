# Ladder round: the finest 44-symmetric sub-object is the FULL SUM (negative, decisive)

Scratch, local-only (`_rb_lad{1..6}.py`). Run from repo root, `PYTHONIOENCODING=utf-8`.
Exact over ℤ (`Xc, Yc, Fk, sds, TAILS` from `_rb_summed.py`). Builds on the FROZEN LADDER of
`_rb_kt_README.md` §2 (C1, 35/35): `F_k(m) = Σ_{ε∈{±1}⁵, r≥0} σ(ε)·Yc(m+2r+1, τ(ε,r))`,
`τ(ε,r)=sorted_desc(ε∘tail_k)+(2r+1)·𝟙₅`, `σ(ε)=sg(ε)·ss(ε)` both m-independent; tail_k=(2k+10,8,6,4,2).

**The atom.** For each ε define the comb `C_ε(m) := σ(ε)·Σ_{r≥0} Yc(m+2r+1, τ(ε,r))`, so
`F_k(m) = Σ_ε C_ε(m)` (reconstruction 35/35, `_rb_lad1`). Defect `D_ε(m) := C_ε(m) − C_ε(44−m)`
on the live window `m∈[12+2k, 32−2k]` (k=0..4; k=5 vacuous). A block B is 44-symmetric ⟺
`Σ_{ε∈B} D_ε(m) = 0` ∀m.

## Verdict: ONLY-FULL-SUM

The **finest 44-symmetric block is the entire sum over all 32 ε** (equivalently all (ε,r)). No
per-ε block, no natural ε-grouping, no (ε,r) sub-layer (r-parity, single H_r), and no *small*
zero-sum subset is individually 44-symmetric. The two load-bearing structural facts:

1. **The reflection defect is genuinely low-rank** — the 32 vectors `D_ε` span only
   `d = 3, 2, 2, 1, 1` dimensions for `k = 0,1,2,3,4` (windows of size 11,9,7,5,3; `_rb_lad2/3`).
   So the whole in-window identity `F_k(m)=F_k(44−m)` collapses to **`d` independent integer
   sum-rules** `Σ_ε a_ε^{(i)} = 0`, where `a_ε ∈ ℤ^d` are the defect coordinates.
2. **Those coordinates `a_ε` are arithmetically generic in ε** — 32 distinct large integers
   (per coordinate) with **no closed form / no ε-invariant structure** (`_rb_lad3` full tables).
   Hence the sum-rule is a *numeric* coincidence over the full 32-set, not a structural cancellation
   that factors through any grouping of ε. The germ we are hunting is these `d` scalar sum-rules
   themselves (a theta/affine-character statement, exactly as `_rb_bnd`/`_rb_kt`/`_rb_level` place it),
   NOT a partition of the index set.

## (a) Per-ε — DEAD

`_rb_lad1`: number of ε with `D_ε ≡ 0` in-window:

| k | window | per-ε symmetric | defective |
|---|---|---|---|
| 0 | [12..32] | 0/32 | 32/32 |
| 1 | [14..30] | 0/32 | 32/32 |
| 2 | [16..28] | 0/32 | 32/32 |
| 3 | [18..26] | 0/32 | 32/32 |
| 4 | [20,22,24] | 2/32 | 30/32 |

The 2 at k=4 are `+++++`,`++++-` gone zero only because the window is the 3-point stub `[20,22,24]`
(rank-1); they are not symmetric at any larger window. Witness (k=0, ε=`+++++`):
`C(12)−C(32) = 37 ≠ 0`.

## (b) Coarser / finer groupings — ALL fail (excluded partitions, with witnesses)

Every natural ε-grouping has ≥1 non-zero-sum block (`_rb_lad2`, all k):

| partition (of the 32 ε) | blocks | result |
|---|---|---|
| `prod(ε)=sg` | 2 | FAILS — +1 block D-sum `[-21,0,-39,0,-20,0,20,0,39,0,21]` (k=0) |
| `σ(ε)=sg·ss` | 2 | FAILS (2 bad) |
| number of minus signs | 6 | FAILS (6 bad) |
| `ε ↔ −ε` orbits | 16 | FAILS (16 bad) — total negation is NOT a symmetry of the comb |
| sign of head (2k+10 entry) | 2 | FAILS (2 bad; the lone k=4 "pass" is the rank-1 stub) |
| tail-multiset τ up to permutation | 1 | = full sum only (all 32 share |·|-multiset {2k+10,8,6,4,2}); the signed multisets are all distinct ⇒ this grouping is *either* the full sum *or* per-ε, nothing between |

**(ε,r)-level groupings also fail** (`_rb_lad4`):
- r-parity: `Σ_ε C_ε^{even-r}` and `Σ_ε C_ε^{odd-r}` are each NOT 44-symmetric (k=0: 6/11 mismatches; every k).
- single H_r layer `H_r(p)=Σ_ε σ(ε)Yc(p,τ(ε,r))`: not palindromic about 22 for any r (k=0, r=2..7 all False; the centres wander, cf `_rb_kt3`).

**Smallest proper zero-sum subset of the 32 combs** (`_rb_lad4/5`):

| k | smallest zero-sum subset | nature |
|---|---|---|
| 0 | **> 6** | none |
| 1 | **> 6** | none |
| 2 | **5** (exactly one) | large circuit, no shared ε-structure |
| 3 | 2 (one pair) | **accidental** `a_ε=−a_{ε'}` in the rank-1 stub |
| 4 | 2 (five pairs) | **accidental** (rank-1 stub) |

The size-2 pairs at k=3,4 cross every ε-invariant, so they are numeric accidents, not structure:
- k=3: `++--+ (sg+1) ↔ +-+-- (sg−1)` — different sg, different #minus (2 vs 3).
- k=4: e.g. `+++-- (sg+1,σ−1) ↔ +-+-+ (sg+1,σ+1)`; `+-+-+ ↔ -+-+-` — cross sg, σ and #minus.

So even the *finest arithmetic* zero-sum block is ≥5 (k=2) or >6 (k=0,1), and where it drops to 2
(k=3,4) it is a rank-1 coincidence pairing structurally unrelated ε. No structural sub-block exists.

## (c) Profiles

Per-ε comb profile `c_ε(m) := Σ_r Yc(m+2r+1, T(ε)+(2r+1))` (`σ` stripped), k=0 (`_rb_lad5`):
- **None is palindromic** — palindrome-centre search returns None for all 32 ε (about any centre, either sign).
- Support intervals *vary with ε*: `m∈[−4,24]` (6 nz) up to `[−4,32]` (9 nz); the lower/upper edge and the
  count both move with the sign pattern. So the combs are not translates of one profile.
- Frozen tails `T(ε)=sorted_desc(ε∘tail_k)` are **all 32 distinct** (distinct abs values ⇒ ε ↔ T(ε) is a
  bijection); every `T(ε)` carries the same |·|-multiset {10,8,6,4,2}, i.e. all 32 lie in one W(C₅)
  sign-orbit, but `Yc` is only **type-A (S₅) antisymmetric on the tail, NOT sign-flip symmetric** (an
  isolated coordinate flip is not a symmetry — `_rb_bnd` L3). Hence the sign-orbit does NOT force the
  profiles to coincide up to sign: the only exact relation is total negation `Yc(−ν)=−Yc(ν)`, and even
  that fails to relate `c_ε` to `c_{−ε}` cleanly (the head coordinate `m+2r+1` is not negated; supports of
  `c_ε` vs `c_{−ε}` differ, e.g. `+++++`:[0,12,16,20,24] vs `-----`:[0,4,8,12,16,20,24,28]). No profile-level
  palindrome to exploit.

## (d) Structural characterization of the (only) symmetric block

The only 44-symmetric block is `{all 32 ε}`. Its content is the rank-`d` sum-rule `Σ_ε a_ε = 0`
(`a_ε ∈ ℤ^d`, `d=3,2,2,1,1`). The germ of the identity is therefore **not** "what do the ε in a block
share" (there is no proper block) but **why the `d`-vector of defect functionals annihilates the whole
32-set** — a single affine/theta statement on the deleted-top-weight discriminant, matching the
`_rb_level` pin (character `(−,+)`, wall 11, step 22) and the `_rb_kt` handoff to route 1
(A₁₂⁽²⁾ Macdonald / Weyl-Kac). The reduction of the |window|-many equations to just `d` of them
(`5,4,3,2,1` naive reflection constraints collapse to `3,2,2,1,1`) is the one new quantitative asset here:
the reflection has **depth d**, the outer near-wall constraints being automatically satisfied by support
thinning.

## (e) Domain fence — respected

`_rb_lad6`, three-row tail `(2k+10, 2l+8, 6, 4, 2)`, reflect `m↦44−m` on the clean window:
- **Full sum HOLDS** at l=0,1,3,4 (0 mismatches, every reachable λ₂).
- **Full sum BREAKS** at l=2 (λ₂=2: 4/4; λ₂=3: 2/3; λ₂=4: 2/2) — e.g. `k2=2,l=2,m=16`: `1 ≠ 0`.
  Full-sum defect vector `(1,0,1,0,−1,0,−1)`, defect rank 2, per-ε symmetric 1/32.

So the only block we found (the full sum) correctly breaks at the ρ4-tail resonance l=2, exactly
reproducing `_rb_kt6` C4 / `_rb_schur` through the C_ε machinery. The block is not an artifact: there is
nothing finer to break, and the full sum itself fails at l=2 as the fence demands.

## Files
- `_rb_lad1.py` — comb `C_ε`, reconstruction (35/35), per-ε defect (a).
- `_rb_lad2.py` — structured-partition tests + defect rank (3,2,2,1,1).
- `_rb_lad3.py` — defect coordinates `a_ε` (full tables, all distinct, no closed form).
- `_rb_lad4.py` — r-parity split, single H_r layer, minimal zero-sum subset search (≤4).
- `_rb_lad5.py` — accidental-pair diagnosis (k=3,4), circuit size 5,6 push, per-ε profiles (c).
- `_rb_lad6.py` — three-row fence (e): full sum holds l≠2, breaks l=2.

## Net
Per-ε: DEAD (0/32). Finest structural level: **ONLY-FULL-SUM**. Excluded, with witnesses: per-ε,
`prod(ε)`, `σ(ε)`, #minus, `ε↔−ε`, sign-of-head, tail-multiset, r-parity, single-H_r; smallest arithmetic
zero-sum block ≥5 (k=2) / >6 (k=0,1), only accidental size-2 in the rank-1 stubs (k=3,4). The two facts
that matter: (1) the reflection defect is only **rank 3,2,2,1,1** (identity has depth d, not window-width);
(2) the defect coordinates `a_ε` are **arithmetically generic in ε** (no grouping factors the cancellation).
The residual is the affine/theta node, unchanged in kind from the level/kt/bnd rounds.
