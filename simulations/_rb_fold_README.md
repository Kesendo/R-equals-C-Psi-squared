# Fold round: the affine (translation-carrying) symmetry search — bankable NULL, both faces

Scratch, local-only (`_rb_fold*.py`). Run from repo root, `PYTHONIOENCODING=utf-8`. All arithmetic
exact over ℤ. Gates: `python simulations/_rb_fold1.py` (Xc, ~1 min, `ALL GREEN`) and
`python simulations/_rb_fold2.py` (Φ, ~2 min, `ALL GREEN`). Builds on `_rb_summed.py` (`Xc`) and
`f133_w_closed_form.build_halves`.

Context: the two-row law `F_k(m) = F_k(44−m)` is proved as a finite object (`_rb_bnd`, five
palindromic `P_k`) but its **structural** (affine-character, generalizing) reason stays blocked. Prior
rounds pinned only the LINEAR symmetry of the coefficient array (`S₆ × ℤ₂(total-neg)`, `_rb_bnd` L3)
and showed a coordinate-1-only affine reflection is a symmetry for no center. The **untested space**
this round closes: affine maps `ι(v) = A·v + b` on the FULL 6-dim exponent lattice with `A` a signed
permutation AND `b` a genuine nonzero integer translation (possibly on tail coordinates), i.e. the
repo's own fold-lattice template (relabel `A` + gauge sign + a shift `b` carried on a subset;
`PROOF_CODIM1_BY_ADDITIVITY §7`, `Mirror.cs`).

## Headline — NULL, by a centroid obstruction that is GENERAL

There is **no translation-carrying signed-permutation affine symmetry** of either the 6-dim array
`Xc` or its 2-variable contraction `Φ`. A one-line reason kills the whole space and every gauge
refinement of it:

> **Centroid theorem.** `Xc(−v) = Xc(v)` EXACTLY (total negation multiplies each of the `15 + 31 = 46`
> antisymmetric sheet factors of `X = Dhat·∏₃₁` by `−1`, net `(−1)⁴⁶ = +1`). So the finite support
> `S = supp(Xc)` satisfies `S = −S`, hence `Σ_{v∈S} v = 0`. If `ι(v) = A·v + b` is a symmetry
> (`Xc(A·v+b) = ±Xc(v)` ∀v, any nonvanishing gauge), then `ι` bijects `S` onto `S`, so
> `Σ_{v∈S} ι(v) = A·0 + |S|·b = |S|·b` must equal `Σ_{v∈S} v = 0`. Therefore **`b = 0`**.

The same identity holds for `Φ` (`Φ(−e) = Φ(e)` exact, the tail negation absorbed into the `W(C₄)`
sum), so `Φ`'s translation is forced to `0` too. The affine group collapses to the linear group in
both cases. The `s₀`-translation (`m ↦ 44−m`, `e₁ ↦ 44−e₁`) that the two-row law needs is therefore
**not a global symmetry of the array** — it lives only in-window, exactly the obstruction `_rb_bnd`
L3/L4, `_rb_level`, `_rb_kt`, `_rb_lad`, and `_rb_def` all isolate. The residual stays a
theta / affine-C₆ character identity; no relabel-and-shift fold reaches it.

## Task 1 — Xc, the 6-dim array (`_rb_fold1.py`)

Support `B` = the F_k read-off region (all `sds([m, ε∘tail_k]) + c·𝟙`, `k=0..5`, `m∈[0,44]`,
`ε∈{±1}⁵`, `c∈[−8,8]`), **5523** nonzero exponent vectors.

- **Centroid / total-neg even**: `Xc(−v) = Xc(v)`, `5523/5523` (algebraic `(−1)⁴⁶ = +1`).
- **Linear group (b = 0), exact.** Since `Xc(P·diag(ε)·v) = sgn(P)·Xc(diag(ε)·v)` (`S₆`-antisym), a
  signed permutation is a symmetry iff its sign-pattern `diag(ε)` is. Testing all **64** patterns:
  exactly `ε = (+,+,+,+,+,+)` (identity, sign `+1`) and `ε = (−,…,−)` (total negation, sign `+1`) —
  flip-counts admitting a symmetry are `{0, 6}` only. **A partial/single flip is never a symmetry**
  (reconfirms `_rb_bnd` L3). Group order `720 · 2 = 1440 = S₆ × ℤ₂(total-neg)`.
- **The affine scan, b ≠ 0**: `8192` `(A,b)` pairs — the four prioritized directions
  `𝟙`, `e₁`, `(1,1,0,0,0,0)`, `(1,−1,0,0,0,0)`, each `t ∈ [−8,8]∖{0}`, × all `64` sign patterns × two
  permutations (identity, transposition) — **0 hits**. Consistent with the centroid theorem, which
  covers all `b` (not just the box).
- **Gauge extension**: a per-coordinate gauge `g(v) = ∏ g_i(v_i)` cannot rescue `b ≠ 0` — `g` is
  nonvanishing, so `ι` still bijects support to support and `|S|·b = 0` still forces `b = 0`.

**Complete affine symmetry group of Xc = `S₆ × ℤ₂`, order 1440, `b ≡ 0`.** This extends `_rb_bnd`'s
linear result to the full affine space: the translation direction is empty.

## Task 2 — Φ, the 2-variable contraction (`_rb_fold2.py`)

`Φ(e₁,e₂) = Σ_{w∈W(C₄)} sgn(w)·Xc(e₁,e₂,w(tail))`, `tail=(8,6,4,2)`. **Rebuilt on the full even
lattice** `e₁,e₂ ∈ [−36,36]` → **356** nonzero points.

> **Truncation correction (real).** The committed scratch `_rb_phi_dict.pkl` was built on the box
> `|e₂| ≤ 28` and holds only **340** points; **16** genuine support points with `|e₂| > 28` were
> clipped. Those 16 are exactly the "16 swap-anti exceptions" the earlier `_rb_phi_analyze` probe
> reported — an artifact of the box, not a real asymmetry. The read-off `nread` samples `Φ` only at
> `e₂ = 2(k+5) ∈ [10,20]`, well inside the box, so **the committed `n_λ` / two-row results are
> unaffected**; only the pkl-level global-symmetry probes were corrupted. Use the full-lattice build.

- **Centroid / total-neg even**: `Φ(−e) = Φ(e)`, `356/356` (exact).
- **Linear group (b = 0), order 4.** Of the 8 signed `2×2` permutations, exactly four are symmetries:
  identity `(+)`, **swap `(e₁,e₂)↦(e₂,e₁)` with sign `−1`**, total negation `(+)`, and swap∘negation
  `(−1)`. This is `S₂(rows) × ℤ₂(total-neg) ≅ V₄` (Klein four; the 2-row shadow of `Xc`'s
  `S₆ × ℤ₂`), sign character = swap parity, verified closed. The **swap-antisymmetry is exact on the
  full support** (`Φ(e₂,e₁) = −Φ(e₁,e₂)`, forced by `Xc`'s `S₆`-antisym transposition); the old pkl's
  16 "exceptions" were the truncated points.
- **The affine scan, b ≠ 0**: `672` `(A,b)` pairs — box `[−8,8]²` (even) × 8 linear parts, **plus the
  `s₀` / anti-period translations `e₁ ± 22`, `e₁ ± 44`** — **0 hits**.
- **The `s₀` reflection is window-only.** `e₁ ↦ 44−e₁` at fixed `e₂` fails **346/356** globally
  (sign `+`) — not a symmetry of `Φ`. Yet the in-window read-off reflection `nread(l₁,k) =
  nread(10−l₁,k)` holds **42/42**. So the "border constancy" the brief flagged **is** a translation
  statement, but an **in-window** one, not a global affine symmetry — the same split as `F_k` itself.

**Complete affine symmetry group of Φ = `S₂ × ℤ₂`, order 4, `b ≡ 0`.**

## Tasks 3 / 4 — the derivation vs. the bankable negative

No translation-carrying symmetry exists (tasks 1, 2 both NULL), so task 3 (deriving the law from
such a map) is vacuous. Per task 4, the searched spaces, exactly:

| object | linear group (b=0) | translations tested (all NULL) | general kill |
|---|---|---|---|
| `Xc` (6-dim) | `S₆ × ℤ₂`, order **1440** | 4 prioritized dirs × `[−8,8]` × 64 signs × 2 perms = **8192** `(A,b)` | centroid: `b=0` ∀ signed-perm `A`, ∀ gauge |
| `Φ` (2-var) | `S₂ × ℤ₂`, order **4** | box `[−8,8]²` × 8 D₄ parts + `e₁±22,±44` = **672** `(A,b)` | centroid: `b=0` ∀ signed-perm `A`, ∀ gauge |

The centroid theorem makes the negative bankable beyond the tested boxes: for **any** invertible
linear `A` (`A·0 = 0`) and **any** nonvanishing gauge, a symmetry of an origin-symmetric finite
array forces `b = 0`. The fold-lattice template needs a linear part that does NOT fix the origin
(a genuine affine reflection whose center is a support point, not the centroid) — and no signed
permutation supplies one while the array's mass is centered at 0. The `s₀` reflection has its center
at `22`, off the centroid `0`; it can only ever be a **window** identity, which is precisely the
open theta / affine-C₆ statement, unchanged in kind.

## Files
- `_rb_fold1.py` — Xc: support build, centroid/total-neg, 64-pattern linear group (1440),
  b≠0 prioritized scan (8192, empty), gauge note.
- `_rb_fold2.py` — Φ: full-lattice rebuild (356, +16 vs the truncated pkl), centroid, D₄ linear
  group (order 4, swap-antisym exact), b≠0 scan incl. `e₁±22,±44` (672, empty), `s₀` global-vs-window.
