# The F115 Obstruction-Size Distribution: the d=0 reduction and the mirror (x=1) odd-weight view

**Date:** 2026-06-08
**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Status:** Tier 2 (computed, in progress). The d=0 reduction, the mirror condition `a(1)+b(1)=1`, and the
Δ-organization are bit-exact on the grids below; the **Δ-bucket count closed form**
(`n(Δ,k)=2^(Δ−1)·a(k−Δ+1)`) is **derived** (three-bijection proof of the doubling recurrence + the proven
`Σ_Δ n=B(k)`); the per-size **floor (size 3) is closed** `(k−1)²(k−2)/2`, the monomial column and repunit
ceiling are closed, and the per-size **middle (size ≥ 5) is located as a number-theoretic hard core** (weighted
coprime polynomial pairs in GF(2)[x]), not closed (Finding 6).
**Builds on:** [F115 / ANALYTICAL_FORMULAS.md](../docs/ANALYTICAL_FORMULAS.md) (the windowed-hardness GF(2)[x]
theory), [PROOF_F103 §7.7-§7.9](../docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md), and the lens of
[ZERO_IS_THE_MIRROR.md](../hypotheses/ZERO_IS_THE_MIRROR.md) (x=1 = the mirror / DC / Perron point).
**Scripts:** [`simulations/_f87_middle_distribution.py`](../simulations/_f87_middle_distribution.py),
[`simulations/_f87_oddweight_filtration.py`](../simulations/_f87_oddweight_filtration.py),
[`simulations/f87_per_delta_saturated.py`](../simulations/f87_per_delta_saturated.py) (the bucket-count
derivation + the per-size table),
[`simulations/f87_size_cells.py`](../simulations/f87_size_cells.py) (the per-size kernel: floor, monomial
column, repunit ceiling, the located hard core).

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

## Finding 4: the Δ-bucket counts close (the Δ-refinement of B(k))

The d=0 layer splits cleanly by Δ, the buckets sum to B(k), and the **per-bucket count is closed form**
(derived; verified bit-exact k = 4..9, `f87_per_delta_saturated.py`):

| k | Δ=1  | Δ=2  | Δ=3  | Δ=4 | Δ=5 | Δ=6 | Δ=7 | Σ = B(k) |
|---|------|------|------|-----|-----|-----|-----|----------|
| 4 | 8    | 4    |      |     |     |     |     | 12   |
| 5 | 30   | 16   | 8    |     |     |     |     | 54   |
| 6 | 116  | 60   | 32   | 16  |     |     |     | 224  |
| 7 | 458  | 232  | 120  | 64  | 32  |     |     | 906  |
| 8 | 1824 | 916  | 464  | 240 | 128 | 64  |     | 3636 |
| 9 | 7286 | 3648 | 1832 | 928 | 480 | 256 | 128 | 14558 |

The whole table is generated by one recurrence and one seed:

> **`n(Δ+1, k+1) = 2·n(Δ, k)`** (every bucket at fixed distance-from-top doubles with k), and with the
> proven `Σ_Δ n(Δ,k) = B(k)` this telescopes to the Δ=1 seed
> **`a(k) := n(Δ=1, k) = B(k) − 2·B(k−1) = (4^(k−1) + 6k − 16)/9`**, so
> **`n(Δ, k) = 2^(Δ−1) · a(k−Δ+1)`**.

The earlier "low-Δ buckets deviate from powers of two" puzzle dissolves: the high-Δ buckets only *look* like
clean `2^(k-2)`, `2^(k-1)` because the seed values `a(3)=2` and `a(4)=8` happen to be powers of two; the true
factor is `a(m) ≈ 4^(m−1)/9`, which is not. This is a viewpoint-rotation on the already-closed `B(k)`, the
Δ-refinement of it, not new content.

### Derivation of the doubling recurrence

The recurrence `n(Δ+1, k+1) = 2·n(Δ, k)` is not just observed — it is derived. Write a d=0 pair via its
reduced generators: `g = gcd = (1+x)^m` (pure power, m ≥ 1 since both masks have even popcount), `u = `
odd-valuation generator (`v(u)=0`), `w = ` the `v(w)=Δ` generator, coprime. Counting the valid powers m
under the body-k width bound gives

> `n(Δ, k) = Σ_D c_Δ(D)·max(0, k−1−D)`,  `c_Δ(D) = #{coprime (u,w) : v(u)=0, v(w)=Δ, max(deg u, deg w)=D}`.

The functions `{max(0,k−1−D)}_D` are linearly independent in k, so `n(Δ+1,k+1)=2n(Δ,k)` for all k **iff**
`c_{Δ+1}(D+1) = 2·c_Δ(D)` for all D. Split `c_Δ(D) = E + U + W` by which generator carries the top degree
(E = tie `deg u = deg w = D`; U = `deg u = D > deg w`; W = `deg w = D > deg u`). Three bijections close it:

- **(M) gap shift** `(u,w) ↦ (u, (1+x)w)`: a bijection {gap Δ} ↔ {gap Δ+1} (it preserves coprimality
  because `(1+x)∤u`). Pulling `c_{Δ+1}(D+1)` back through M and sorting by degree gives
  `c_{Δ+1}(D+1) = U_Δ(D+1) + E_Δ(D) + W_Δ(D)`.
- **(R) reflection** `(u,w) ↦ (u+w, w)`: when `deg u = deg w = D` the top terms cancel, so `deg(u+w) < D`,
  and `u+w` keeps `v=0` (odd+even popcount) and `gcd(u+w,w)=gcd(u,w)=1`. This is an involution swapping
  tie ↔ w-dominant, so **`E_Δ(D) = W_Δ(D)`**.
- **(P) per-w doubling**: for a fixed `w` with `(1+x)|w` (true, since `v(w)=Δ≥1`), the coprimality
  `gcd(u,w)=1` already forces `v(u)=0`, and the count of `u` of degree exactly `e` coprime to `w` is
  `2^(e−deg w)·φ(w)` (φ = units mod w) — it **doubles** with each step in `e`. Summing over w gives
  `U_Δ(D+1) = 2·(E_Δ(D) + U_Δ(D))`.

Combine: `c_{Δ+1}(D+1) = U_Δ(D+1) + E_Δ(D) + W_Δ(D) = 2(E_Δ(D)+U_Δ(D)) + 2E_Δ(D) = 4E_Δ(D) + 2U_Δ(D)
= 2·c_Δ(D)`. ∎ With the proven `Σ_Δ n(Δ,k) = B(k)`, the recurrence telescopes to `a(k) = B(k) − 2B(k−1)`,
so the bucket-count closed form is **derived**, not merely fitted. (All three ingredients + the assembled
lemma are checked in `verify_doubling_lemma`.)

Within each Δ the per-*size* distribution is W-stable for k ≤ 5 (no cancellation); the size split inside a
bucket is Finding 5.

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

## Finding 6: the MacWilliams kernel — floor closed, ceiling characterized, the hard core located

The saturated obstruction depends only on the reduced pair `(a,b)` (not on `m` or the window once W is large),
so the whole per-size distribution factors the same way the count did:

> `saturated_dist(k)[s] = Σ_D c(D,s)·max(0, k−1−D)`,  `c(D,s) = #{coprime (a,b): v(a)=0, v(b)≥1,
> max(deg a, deg b)=D, minweight(a,b)=s}`,  where **`minweight(a,b)` = the odd-restricted free distance of the
> convolutional code `⟨(a,b)⟩`** = `min_t` odd`[wt(a·t)+wt(b·t)]`.

So `c(D,s)` is the per-size analog of the now-closed `c_Δ(D)`. The deep look (`f87_size_cells.py`,
cap = D+4, cap-stability verified at D=5,6) closes everything that *can* close and locates the rest:

**The floor (size 3) is closed.** `minweight=3 ⟺ the v=0 generator is a monomial x^j and the v=Δ generator
has popcount 2` (forward: `t=1` gives `1+2=3`; converse: `minweight=3` ⟹ some `wt(a·t)=1` ⟹ `a·t=x^i` ⟹ `a`
is a monomial — only `x` is irreducible with that factorization — and `b` can't be the monomial since
`v(b)≥1` ⟹ even popcount). This forces **Δ to be a power of two** (a weight-2 `b=x^p(1+x^r)` has
`v=2^{v₂(r)}`), matching the data (size 3 present at Δ=1,2,4, absent at Δ=3,5). Hence

> `c(D,3) = 3D−1`,  and the d0 size-3 total `T₃(k) = Σ_{D=1}^{k−2}(3D−1)(k−1−D) = (k−1)²(k−2)/2`.

(Cross-check: via the d-reduction `T_full(k) = T₃(k) + Σ_{d≥1} 2^{d-1}T₃(k−d)` reproduces the catalogued full
triangle count `TriangleHardMaskCount`: 11, 37, … at k=4,5.)

**The monomial column is closed (polynomial).** For a fixed even weight β, `#{(x^j, weight-β b) coprime,
max deg D}` is a polynomial in D of degree `β−1` (leading difference `β+1`: β=2→3, β=4→5, β=6→7). A monomial
generator makes coprimality trivial (`gcd(x^j, ·)`), and the count stays polynomial. **This is exactly why
size 3 closes**: `3 = 1+2` is the only popcount split, and it has a monomial factor.

**The ceiling (max size 2D+1) is the repunit pair.** For D ≥ 4 the count is exactly **2**: one generator is
the full repunit `R_D = 1+x+…+x^D` (all ones), the other is `R_{D-1}` or `x·R_{D-1}` — the densest pairs, no
cancellation. (D=2,3 are small edge cases with 3, 4.) These are the same repunits that mark the Door-1 syzygy
extreme.

**The hard core is located, not closed.** The irregularity enters at the *first* popcount split with **both
popcounts ≥ 2** — `(3,2)` at size 5 — and that cell is genuinely **not polynomial in D**
(`3,16,51,114,215,348,556,822,1162,…`, erratic differences). It is honest GF(2)[x] weighted coprimality: a
weight-3 `a` is coprime to `b=x^p(1+x^r)` iff coprime to `1+x^r`, which depends on the factorization of
`1+x^r` (number-theoretic). The gcd-formula layer is already non-polynomial there (`f(D,5), f(D,7)`); the −2
convolutional cancellation (k ≥ 6) sits on top. So the MacWilliams middle resists closed form for a concrete
reason: **closed-formness holds exactly while some popcount split carries a monomial; from `(3,2)` on, the
count *is* the distribution of weighted coprime polynomial pairs.**

**Ported to C# (2026-06-08).** The floor closed forms live in `WindowedObstructionScan`
(`TriangleReducedPairCountByMaxDegree` = 3D−1, `TriangleHardCountBaseD0` = (k−1)²(k−2)/2,
`TriangleHardCountByGRestDegree` the d-layering, `Repunit`), recorded on `WindowedHardnessClaim`; the
F81↔F115 orthogonality (the next section) is the typed `AntiFractionObstructionOrthogonalityClaim`
(Tier1Derived, parents `F83AntiFractionPi2Inheritance` + `WindowedHardnessClaim`).

---

## Connection to the bit_b axis (F81/F83): orthogonal coordinates of M's polarity cube

This obstruction theory and the F81/F83 anti-fraction are the two bits of the *same* object — the residual
superoperator M (`PalindromeResidual`) — read on the two axes of its Klein polarity cube (`PolarityCubeMap`,
`Z2Axis`). The pieces are already typed in C#; what follows is only their assembly.

- **bit_b axis (F81/F83/F112).** `PiDecomposition` splits `M = M_sym + M_anti`; the anti-fraction
  `‖M_anti‖²/‖M‖² = 1/(2+4r)` (`F83AntiFractionPi2Inheritance`), `r = ‖H_even_nontruly‖²/‖H_odd‖²`.
  Continuous, hardware-measurable (F83, Marrakesh/Kingston).
- **bit_a axis (F87/F115).** `BitADyadicGrade` makes precise that the Klein `bit_a = popcount(mask) mod 2`
  is the `v=0` bottom rung of a `(1+x)`-valuation tower `v = 0…k−1`; the obstruction theory above is the
  graded refinement of that axis (verified bit-exact against `WindowedObstructionScan.ValuationAtOnePlusX`).
  Discrete, the hard/soft verdict + obstruction size.

**They meet at bit_b = 1 (the diagonal cell), and there the bit_b readout is degenerate.** Diagonal-cell
Mixed terms are pure Π²-odd, so `r = 0`, and F83's anti-fraction sits at its *maximum* ½ exactly
(`F83AntiFractionPi2Inheritance.MaximumAntiFraction`, the F81 Step-8 50/50). So the anti-fraction is pinned
at ½ across the *entire* F115 cell — identical for hard and soft pairs and every obstruction size. The bit_b
magnitude is blind precisely where the bit_a valuation carries all the information: **F81/F83 and F115 are
orthogonal coordinates of M's cube, and F115 is the strictly finer probe.** The link is structural (one M,
one shared bit_b = 1 gate), not quantitative — and it follows from an already-typed result, no new computation.

**Two breaks, one visible to the anti-fraction.** Per `DissipatorResonanceLaw`, hardness lives *inside* the
dissipator's Klein cell (the bit_a obstruction, F115); the transverse-field Brecher breaks it from *outside*
the cell (F89, predicted from F83). The anti-fraction sees the outside-cell Brecher, not the in-cell
obstruction — the two are the two poles of dissipator-letter resonance. (If we ever want this orthogonality
typed rather than recognized, its natural C# home is the open `BitATwin` slot that `PolarityCubeMap` already
tracks for `F83AntiFractionPi2Inheritance`.)

---

## Status ledger

- **Derived (proven):** the Δ-bucket **count** closed form `n(Δ,k) = 2^(Δ−1)·a(k−Δ+1)`,
  `a(k)=B(k)−2B(k−1)=(4^(k−1)+6k−16)/9` (Finding 4, three-bijection proof of the doubling recurrence); the
  size-3 floor `c(D,3)=3D−1`, `T₃(k)=(k−1)²(k−2)/2` with its monomial×weight-2 characterization (Finding 6).
- **Bit-exact (solid):** the d-reduction (Finding 1); the QC-code view (Finding 2); `a(1)+b(1)=1` (Finding 3);
  the cancellation events and the saturated totals (Finding 5); the monomial-column polynomiality (degree β−1)
  and the repunit ceiling (count 2, D≥4) (Finding 6).
- **Located, provably not closed:** the saturated per-size *middle* (size ≥ 5). It is the distribution of
  weighted coprime polynomial pairs in GF(2)[x], non-polynomial from the `(3,2)` cell on (Finding 6).

## Open / next

F115 is now closed except for one piece, and that piece is **identified as number-theoretic, not merely
unfinished**: the saturated middle (size ≥ 5) *is* the count of weighted coprime polynomial pairs in GF(2)[x]
(plus the −2 convolutional cancellation on top). The boundary cells are all closed — the size-3 floor exactly,
the monomial column as polynomials of degree β−1, the repunit ceiling at count 2 — and the closed-form
frontier is sharp: **a size cell is polynomial in D iff some popcount split carries a monomial; the first
irregular cell is `(3,2)` at size 5.**

A generating-function attack (the weight enumerator of `⟨(a,b)⟩` summed over the coprime family, or a
Möbius-over-gcd inclusion-exclusion on the weighted pair count) is the natural next tool if a closed form for
the middle is wanted; but the result here is that the middle's content is exactly a catalogued-hard object,
which is itself the answer to "why the MacWilliams middle resists."
