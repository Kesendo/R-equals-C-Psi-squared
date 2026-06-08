# The F115 Obstruction-Size Distribution: the d=0 reduction and the mirror (x=1) odd-weight view

**Date:** 2026-06-08
**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Status:** Tier 2 (computed, in progress). The d=0 reduction, the mirror condition `a(1)+b(1)=1`, and the
Δ-organization are bit-exact on the grids below; the saturated per-Δ distribution and the Δ-bucket counts
are OPEN.
**Builds on:** [F115 / ANALYTICAL_FORMULAS.md](../docs/ANALYTICAL_FORMULAS.md) (the windowed-hardness GF(2)[x]
theory), [PROOF_F103 §7.7-§7.9](../docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md), and the lens of
[ZERO_IS_THE_MIRROR.md](../hypotheses/ZERO_IS_THE_MIRROR.md) (x=1 = the mirror / DC / Perron point).
**Scripts:** [`simulations/_f87_middle_distribution.py`](../simulations/_f87_middle_distribution.py),
[`simulations/_f87_oddweight_filtration.py`](../simulations/_f87_oddweight_filtration.py).

---

## The question

F115 settles, in closed form, almost everything about the F87 diagonal-cell hardness obstruction: the
maximum size `min(2W−1, 2k−3−2d)`, the total hard count `A203241 = (4^(k-1)−3·2^(k-1)+2)/3`, the d-layered
count `2^(d-1)·B(k-d)` with `B(k)=(4^k−12k+8)/18`, and the size-3 (triangle) sub-count. The one piece left
open is the **full per-size distribution**: how many hard pairs have an obstruction of size exactly s for
the *middle* sizes 3 < s < max. F115 records that this "stays window-dependent" because the obstruction is
the minimum ODD weight of a quasi-cyclic code, and a sparse multiplier can cancel a relation below the
gcd-formula popcount. This experiment opens that thread.

### Setup and vocabulary (self-contained)

A diagonal-cell pair is two nonzero even-popcount k-bit masks `p1, p2` (their X/Y flip positions). On an
N-site chain there are `W = N−k+1` windows; the pair's edge set is the shift family
`{x^w·p1, x^w·p2 : 0 ≤ w < W}`. The **obstruction** is the size of the smallest ODD subset of that set that
XORs to 0 (the minimal odd 𝔽₂-relation), 0 if none. Reading a k-bit mask as a GF(2)[x] polynomial (bit j =
x^j), with `v(·)` the (1+x)-adic valuation:

- **hard ⟺ `v(p1) ≠ v(p2)`** (the F115 one-number criterion);
- `g = gcd(p1,p2)`, `d = g_rest_degree` = degree of the non-(1+x) part of g; reduced generators
  `a = p2/g`, `b = p1/g` (coprime);
- `Δ = |v(p1) − v(p2)|` = the **(1+x)-valuation gap** (the two masks' "distance from the mirror" x=1).

---

## Finding 1: the distribution collapses to the d=0 layer (the d-reduction)

The whole per-size distribution factors through the d=0 layer:

> **full-dist(k, W) = Σ_d 2^(d-1) · D₀(k−d, W)**,  where D₀(k′, W) is the d=0 distribution at body k′.

Equivalently, the layer-d distribution equals `2^(d-1)` times the d=0 distribution at body `k−d`, at the
same window count W. Verified bit-exact for k = 4..7, all layers (`_f87_middle_distribution.py` part 2). The
test compares at fixed W, so it holds in the cancellation regime too: the d-structure is fully understood,
and only **D₀** remains.

This extends the already-closed d-layered *count* `2^(d-1)·B(k-d)` and *size cap* `2(k-d)-3` to the full
*distribution*: the same `2^(d-1)` × body-(k-d) reduction governs every size, not just the count and the cap.

---

## Finding 2: the QC-code (s-multiplier) view is bit-exact

For a d=0 pair, any windowed odd relation factors (a,b coprime) as `q_p1 = a·s`, `q_p2 = b·s`, so the
obstruction is

> the **minimum ODD weight of `{(a·s, b·s) : deg s ≤ W−1−max(deg a, deg b)}`** (a 2-generator quasi-cyclic /
> terminated rate-1/2 convolutional code). The `s = 1` term is the gcd-formula size `popcount(a)+popcount(b)`;
> any `s ≠ 1` that lowers the weight is the cancellation.

Validated against the direct windowed minimal-odd-relation search: 0 mismatches, k = 4,5,6, W ∈ {k−1, k, 2k}
(`_f87_oddweight_filtration.py` part A).

---

## Finding 3: the odd-weight condition IS the mirror statement `a(1)+b(1)=1`

For **every** d=0 hard pair (verified k = 5,6,7, no exceptions):

> **`a(1) + b(1) = 1`** (in GF(2)): exactly one reduced generator vanishes at the mirror point x=1.

Because a,b are coprime and the pair is hard, `{v(a), v(b)} = {0, Δ}`: the generator with v ≥ 1 has even
popcount (vanishes at x=1, `p(1)=0`); the one with v = 0 has odd popcount (`p(1)=1`). So **hardness is the two
reduced generators having opposite vanishing at the mirror**, and this is exactly what guarantees the code's
odd-weight coset is non-empty (every hard pair *has* an odd obstruction). This is the literal
"zero is the mirror" reading of F115 hardness: the verdict lives at x=1, the DC / Perron / mirror mode.

---

## Finding 4: Δ (the mirror gap) organizes the d=0 layer

The d=0 layer splits cleanly by Δ, and the buckets sum to B(k):

| k | Δ=1 | Δ=2 | Δ=3 | Δ=4 | Δ=5 | Σ = B(k) |
|---|-----|-----|-----|-----|-----|----------|
| 5 | 30  | 16  | 8   |     |     | 54  |
| 6 | 116 | 60  | 32  | 16  |     | 224 |
| 7 | 458 | 232 | 120 | 64  | 32  | 906 |

The high-Δ tail is clean: count(Δ = k−2, k) = `2^(k-2)` (8, 16, 32) and count(Δ = k−3, k) = `2^(k-1)`
(16, 32, 64). The low-Δ buckets deviate from powers of two (e.g. k=7: 458, 232, 120 vs 512, 256, 128); their
closed form is an OPEN sub-thread (a pure counting problem). Within each Δ the distribution is W-stable for
k ≤ 5 (no cancellation).

---

## Finding 5: the cancellation is a quantized −2 cascade, only k ≥ 6 (the MacWilliams middle, localized)

The window-dependence of the middle is entirely the cancellation, and it has a sharp shape:

- k ≤ 5: **no cancellation** (gcd-formula size == actual obstruction; the distribution saturates at W = k−1).
- k ≥ 6: the actual obstruction drops below the gcd-formula in **steps of −2** as W grows, saturating at
  large W. Cancellation events (at W = 2k):
  - k=6: `7→5 (×14), 9→7 (×3)`
  - k=7: `7→5 (×50), 9→7 (×50), 9→5 (×8, a double −4), 11→9 (×4)`

Drops are always even (the weight stays odd). At W → ∞ the distribution saturates to the **code's minimum
odd-weight distribution**, a clean coding-theory invariant (no window truncation). Saturated d=0 totals:

| k | saturated D₀ (W→∞) |
|---|--------------------|
| 4 | 3:9, 5:3 |
| 5 | 3:24, 5:26, 7:4 |
| 6 | 3:50, 5:131, 7:41, 9:2 |

So the once-vague "window-dependent middle" is now: **clean Δ-buckets (distance from the mirror) × a −2
cancellation cascade that appears only at k ≥ 6 and saturates at W → ∞**. The cascade is a MacWilliams-type
weight redistribution inside the odd coset, and the x=1 handle (`a(1)+b(1)=1`) is in hand.

---

## Status ledger

- **Bit-exact (solid):** the d-reduction (Finding 1); the QC-code view (Finding 2); `a(1)+b(1)=1` (Finding 3);
  the Δ-buckets summing to B(k) and the high-Δ powers `2^(k-2)`, `2^(k-1)` (Finding 4); the cancellation
  events and the saturated totals (Finding 5).
- **Pattern-observed (not yet derived):** the low-Δ bucket counts; the exact saturated per-Δ distribution; a
  closed form / generating function for the −2 cascade.

## Open / next

1. **The saturated per-Δ distribution** (W → ∞, the code's minimum-odd-weight distribution) as a closed form
   or generating function: the MacWilliams target, now with `a(1)+b(1)=1` and Δ as the axes.
2. **The Δ-bucket counts** (the low-Δ deviation from `2^j`): a separate, lighter counting problem.

Either path closes the last open piece of F115. The mirror lens (x=1) reduced the problem from "the
window-dependent middle is messy" to two sharp, named sub-questions.
