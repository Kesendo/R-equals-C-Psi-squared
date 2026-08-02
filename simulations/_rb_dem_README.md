# The l=2 break, demoted: it is not a resonance, it is support truncation + coefficient asymmetry

Scratch, local-only (`_rb_dem{1..8}.py`). Run from repo root, `PYTHONIOENCODING=utf-8`,
exact over ℤ. Context: `_rb_schur_README.md` (the law: `n_raw` invariant under the affine
reflection `s₀: μ₁ ↦ 22−μ₁`, μ = λ+ρ, ρ = (6,5,4,3,2,1); three-row domain l ∈ {0,1,3} hold,
l = 2 breaks 8/16), `_rb_level_README.md` (the pin: s₀ at wall μ₁ = 11, step 22, character
`(−,+)`), and the F130 TEMPLATE (a resonance-bound law re-proved as the codim-drop of a wider
identity). Machinery: `f133_w_closed_form.build_halves` + `coeff_mu` (the alternant read-off as
a function of any integer μ; `coeff_mu(λ+ρ) == n_raw(λ) == 2·n_λ`, sanity-gated 143/143).

## Verdict

**PARTIAL / NULL.** The "l = 2 resonance" is DEMOTED, but the wider law is **not** a clean
codimension-1 identity. Every coordinate/collision/wall/modular/multiplicative predicate — the
task's (a)–(f) and 311 auto-generated atoms — **fails** to separate break from hold. The break
does have one sound geometric *sufficient* condition (support truncation), but its deepest face
is a pure symplectic-character coefficient asymmetry with **no coordinate shadow**.

## The correct object: the s₀-fixed TAIL (not l, not μ₃)

s₀ moves μ₁ only, holding the tail τ = (μ₂,μ₃,μ₄,μ₅,μ₆). So `n_λ = n_{s₀λ}` is a statement about
each tail: fix τ, ask whether `f_τ(μ₁) := n(μ₁,τ)` is symmetric about μ₁ = 11. A reflection pair
is genuine only when **both** μ₁ and 22−μ₁ are dominant (both > μ₂); non-dominant μ₁ are Weyl
images, not weights, and must be excluded (this was the first correction — including them falsely
flags two-row tails as breaks).

The schur-README domain is the **s₀-closed alcove strip** λ₁+λ₂ ≤ 10 (⟺ the partner 10−j ≥ k is
also dominant); it reproduces the counts 36/25/16/9 for l = 0..3 and the l = 2 8/16 break exactly
(`_rb_dem3.py`). The tail reframing reproduces the schur result on the nose: for **three-row**
tails τ = (μ₂,μ₃,3,2,1), **μ₃ = 6 ⟺ break** (0 exceptions: holds have μ₃ ∈ {4,5,7,8}, breaks have
μ₃ = 6 at μ₂ ∈ {7,8,9}). This confirms the tail is the right object.

## The extended domain table (the headline)

Per-weight, within-strip (`_rb_dem3.py`, 672 weights): **BREAK 124, HOLD 60, VAC 488.**
Reorganized by the s₀-fixed tail, μ₂ ≤ 11 (`_rb_dem4/5.py`, dominant-pairs-only):

| class | count | meaning |
|---|---|---|
| **break-tails** | **54** | `f_τ` asymmetric about 11 for some dominant pair |
| non-trivial hold-tails | 16 | ≥1 non-center dominant pair, all symmetric |
| trivial hold-tails | 9 | μ₂ ≥ 10: only the center μ₁ = 11 is a dominant pair |
| VAC | 383 | `f_τ ≡ 0` on the dominant window |

Three-row l-slice, within strip (extends the known l ∈ {0,1,3} hold / l = 2 break): l = 4 and
l = 5 **also hold** (0 breaks) — the resonance is l = 2 only *within three-row*, but that is a
projection artifact (below). Four/five/six-row breaks are pervasive (measured up to 6 rows).

**`μ₃ = 6` / `l = 2` is REFUTED beyond three rows.** Explicit counterexamples: the weight
λ = (4,3,2,2,1) has μ₃ = 6 and **holds**; the tail τ = (8,6,5,3,1) has μ₃ = 6 and is a **hold-tail**.
Same-tail accidental matches also occur: tail (6,5,4,2,1) holds on the pair (μ₁ = 7 ↔ 15) but
breaks on (μ₁ = 9 ↔ 13) — so a per-weight "hold" can be an accidental coincidence inside a
break-tail. The l = 2 fact is the three-row shadow of a much larger, non-resonance-shaped locus.

## The predicate hunt: NULL for geometry

Exhaustive search on the tail (`_rb_dem5/6.py`), 311 atoms — coordinate values/≤/≥, adjacent
and non-adjacent gaps, `μ_a = 2μ_b` (multiplicative), all pair-sums, mod-2 and mod-3 of every
coordinate and of Στ, support-size, `f(11) ≠ 0`:

- single-atom exact separators: **0**
- 2-conjunction exact separators: **0**
- 3-conjunction exact separators: **0**

The task's specific candidates all die: (a) coord = 11 — impossible in-strip (all coords ≤ μ₂ ≤ 10);
(b) pair-sum = 22 — impossible in-strip (μ₁+μ₂ ≤ 21); (c) `22−μ₁ ∈ tail` — 0/54; (d) `2μ₃ = 12`
i.e. μ₃ = 6 — refuted above; (e) reflected-vector repeat / frozen-tail hit — 0/54; (f) any
`μ_a = 2μ_b` — flags all 60 holds too. **The break is not a coordinate-wall collision.**

## The mechanism: three tiers (`_rb_dem7/8.py`)

All 16 non-trivial holds have **symmetric support endpoints** (lo + hi = 22). The 54 breaks split:

- **T1 = 39 — support truncated** (lo + hi ≠ 22). The sp(12) character support in μ₁ is cut on
  the high side below the reflection of the low end, so s₀ maps a supported weight to an
  unsupported one. → **SOUND SUFFICIENT PREDICATE:** *support of `f_τ` asymmetric about 11 ⟹ break*
  (39/54 breaks, **0** false positives on all 25 holds).
- **T2 = 1 — support-set gap** (endpoints symmetric, the set of nonzero μ₁ is not 11-symmetric).
- **T3 = 14 — pure coefficient asymmetry.** Support **set** identical under μ₁ ↦ 22−μ₁; only the
  integer coefficients differ (e.g. τ = (7,5,4,3,2): f(8) = +2 vs f(14) = −2). On the matched
  subset (**14 T3-breaks vs 16 holds, byte-identical symmetric support sets**) the exhaustive
  single/2/3-conjunction search again returns **0 separators**. These 14 are representation-theory
  with no geometry.

Defects `g_τ(μ₁) = f_τ(μ₁) − f_τ(22−μ₁)` are always **even**, `|g| ∈ {2,4,6,8,12,14,16}`,
concentrated on the near-center pairs (μ₁ ∈ {7,8,9,10} ↔ {15,14,13,12}).

## The F131 transfer: NOT found

Probe (`_rb_dem7.py`): does the defect `g` equal ±`f` at a box-modified (King-modification) tail?
Only **7/27** loose matches, no consistent σ. No `n(μ) ± n(σμ) = 0` doubling law was found. The
law does not visibly "double/transfer"; it simply fails where the character is asymmetric.

## The demotion statement (F130-style)

`n_λ = n_{s₀λ}` is **not** a resonance-gated law with l = 2 as the isolated breaking case. It is:
*the tail-fixed slice `f_τ(μ₁)` of the sp(12) = C₆ coefficient of S = X/A_ρ is symmetric about the
affine wall μ₁ = 11.* This holds for a tail iff **(i, geometric, sufficient)** the character support
in μ₁ is symmetric about 11 **and (ii, non-geometric)** the coefficients on that support are
symmetric. The l = 2 break is merely the first three-row tail where (i) fails; it is a special case
of pervasive support-truncation (T1) plus deep coefficient asymmetry (T3) — resonance demoted from
*cause* to *first visible instance*, exactly the F130 shape.

Where this **departs** from F130: F130's wider law was a clean codimension-1 identity. Here the
wider law's residual — the T3 coefficient asymmetry (14 tails) — has **no coordinate predicate**
(NULL under a 311-atom exhaustive search on a matched subset). The honest wider statement is
representation-theoretic: S = X/A_ρ is a *virtual/truncated* affine object, s₀-symmetric only on
the interior of its support and only where its truncation to the level does not perturb the
coefficients. Proving *that* is the affine-character-truncation problem the schur/level READMEs
already name; this campaign rules out the coordinate-geometry shortcut and pins the exact break
inventory (54 tails, 3 tiers) it must reproduce.

## Files
- `_rb_dem1.py` — build_halves + pickle cache + `coeff_mu(μ)` (any integer μ).
- `_rb_dem2.py` — first full map (over-wide; superseded by the alcove-strip restriction).
- `_rb_dem3.py` — within-strip per-weight map (672; reproduces 36/25/16/9, l=2 8/16) → `_rb_dem_strip.json`.
- `_rb_dem4.py` — tail reframing + f-value cache → `_rb_dem_fvals.json` (the expensive step, ~4 min).
- `_rb_dem5.py` — exhaustive predicate search + defect table.
- `_rb_dem6.py` — 311-atom search (multiplicative/mod-3/affine) + (k,k,0,0,0) mod-3 probe.
- `_rb_dem7.py` — truncation-vs-interior split + F131 transfer probe.
- `_rb_dem8.py` — three-tier decomposition + matched-subset separation (14 vs 16 → 0).
- Caches: `_rb_dem_cache.pkl` (halves; self-generated, safe), `_rb_dem_fvals.json`, `_rb_dem_strip.json`, `_rb_dem_map.txt`.
