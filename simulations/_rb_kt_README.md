# Route 2 (Koike-Terada / King modification rule): the finite half is the rule,
# the affine half is not

Scratch, local-only (`_rb_kt{1..6}.py`). Run from repo root, `PYTHONIOENCODING=utf-8`.
Exact over ℤ (dicts `Xc, Yc, Fk, sds` from `_rb_summed.py`; `build_halves`/`n_raw` from
`f133_w_closed_form.py`). Context: the open two-row law `n_(j,k)=n_(10−j,k)`
(`μ₁ ↦ 22−μ₁`), pinned by `_rb_schur_README.md` (0/36), `_rb_summed_README.md` (finite
half proved, affine translation-44 the residual), `_rb_level_README.md` (wall 11, step 22,
character `(−,+)`, coordinate peg `l+h∨=11`).

**Verdict up front.** Route 2 does NOT close the law. The Koike-Terada / King
modification rule is a **finite W(C₆) operation**; applied to the read-off it delivers
**exactly the oddness** `F_k(−m)=−F_k(m)` — the finite-Weyl half the summed-README already
PROVED. The reflection `μ₁ ↦ 22−μ₁ = (oddness about 0) ∘ (translation by 44)`; the residual
**translation-by-44 is the affine coroot-lattice translation, NOT a finite W(C₆) element**,
so no modification rule can produce it. The re-indexing hope of task (c) is REFUTED at the
term level in all coordinates (Xc, Yc, `(m,c)`): the two sides carry different `|Yc|`-
multisets with equal sums. The one positive contribution is a new structural reduction
(§2, C1) that pins the residual to a **first-coordinate-only summed identity with a frozen
tail-skeleton**; that residual is affine and belongs to route 1 (A₁₂⁽²⁾ Macdonald /
Weyl-Kac), not to route 2.

---

## §1 The modification rule, stated (task a)

For the symplectic type C₆, the alternant read-off is the W(C₆)-antisymmetrization. For any
integer vector `ν ∈ ℤ⁶`, with `A[ν] := Σ_{w∈W(C₆)} sgn(w)·t^{w·ν}` (`W(C₆)=S₆⋉{±1}⁶` acting
by signed permutation):

- **wall rule:** `A[ν] = 0` if some `ν_i = 0`, or `|ν_i| = |ν_j|` for `i≠j` (ν on a wall);
- **reduction rule:** otherwise `A[ν] = sgn(w_ν)·A[ν⁺]`, where `w_ν ∈ W(C₆)` is the unique
  signed permutation carrying `ν` into the strictly dominant chamber
  `ν⁺₁ > ν⁺₂ > … > ν⁺₆ > 0`, and `sgn(w_ν) = (sign-flip parity)·(sort permutation parity)`.

At the coefficient level (X is W(C₆)-antisymmetric, §3 of PROOF_F133): `[X]_ν =
sgn(w_ν)·[X]_{ν⁺]}` and `= 0` on the walls. This IS the code: `sds(v)` sorts by value and
returns `None` on a repeat (wall) or `(sorted, psign)` otherwise, and the `ε∈{±1}` loop
supplies the sign flips. The read-off `n_raw(λ) = Σ_ε sgn(ε)[X]_{ε∘2(λ+ρ)}` is nothing but
this rule executed at `μ = λ+ρ`.

**Crucial genre fact.** Every case of the rule is a FINITE-Weyl move. Its only reflection
hyperplanes are the finite ones `ν_i = ±ν_j`, `ν_i = 0`. In the one-variable read-off
`Φ_k(p) := F_k(k,p)` (first coordinate `p`, tail frozen) the sole finite reflection centre
is `p = 0`, giving `Φ_k(−p) = −Φ_k(p)`. There is no finite reflection at `p = 22`
(verified, §4 C3). The King over-length ("n-modification", border-strip) rule is likewise
finite and does not apply here (`λ` has `≤ 6` rows).

## §2 The ladder bookkeeping (task b)

**Telescope (verified `_rb_kt1.py` V1, 12/12 forward + 12/12 backward).**
`Xc(ν) = Σ_{r≥0} Yc(ν+(2r+1)𝟙) = −Σ_{r≥0} Yc(ν−(2r+1)𝟙)`, `𝟙=(1,…,1)`. Each rung adds the
same odd amount to ALL six coordinates. This is the Laurent expansion of `1/sin s`
(`s = s_{1⁶}`); the geometric ratio `M^{−2}` (`M = t^𝟙`) is the affine coroot translation by
`2` in every coordinate.

**Shared-tail reduction (the new structural fact, `_rb_kt6.py` C1, 35/35 exact).** Feeding
the forward telescope into `F_k(m) = Σ_ε σ(ε)·Xc(sds([m, ε∘tail]))` and using that `m`
exceeds the tail-max (so `m` is the strict first coordinate) gives

```
F_k(m) = Σ_{ε∈{±1}^5}  σ(ε) · Σ_{r≥0}  Yc( m+2r+1 , τ(ε,r) ),
         τ(ε,r) = sorted_desc(ε∘tail_k) + (2r+1)·𝟙₅ ,   σ(ε) = sg(ε)·ss(ε).
```

Both `τ` and `σ` are **independent of `m`**. So the reflection `m ↦ 44−m` moves ONLY the
first coordinate `m+2r+1 ↦ (44−m)+2r+1`; the weighted 5-entry tail-ladder
`{(σ(ε), τ(ε,r))}` is FROZEN and shared identically by the two partners. This was seen
first raw in `_rb_kt2.py`: at `k=0`, `m=16` vs `m'=28`, the term-lists have identical tail
patterns and identical weights `w`, differing only in the first coordinate (`17` vs `29`,
offset `m'−m=12`).

**Where the rungs re-enter, and the signs.**
- `σ(ε,r)=σ(ε)` (sign-of-sort of the 5 flipped tail entries), r-independent.
- First rung is `r=0` at first coordinate `m+1`; rung `r` sits at `m+2r+1`. Rungs survive
  while inside `Yc`'s coordinate wall (`≈ 37`), so the SHORT ladder is the LARGE-`m`
  partner (few rungs) and the LONG ladder is the SMALL-`m` partner (many rungs). E.g.
  `k=0 (16,28)`: `m=16` has rungs `r=2..5`, `m'=28` has `r=2,3` only (`_rb_kt1.py`,
  `_rb_kt3.py`).
- After the tail-antisymmetrization the first two rungs cancel: `T_0=T_1=0`; the ladder
  starts at `r=2` (observed all pairs).
- The row functions `H_r(p) := Σ_ε σ(ε)·Yc(p, τ(ε,r))` have a LOCAL palindrome centre at
  `p = 2r+1` (= the shift) plus a one-sided decaying tail toward the wall; the centres WANDER
  with `r`, so there is no single per-`r` palindrome (`_rb_kt3.py`, palindrome test NONE for
  every `r`). This is the summed-README's "the H^(r) centres wander", reconfirmed.

## §3 The re-indexing claim (task c): PARTIAL, then REFUTED

**Finite half — PROVED, and it IS the modification rule.** The finite W(C₆) reduction of
§1 gives, on the one-variable read-off, `Φ_k(−p) = −Φ_k(p)` and nothing else (only centre
`0`). `_rb_kt4.py`: `Φ_k(p)` is odd about `0` for all `k` (exact, all `p`). This equals the
summed-README's L1·L2·L3. So the modification rule reproduces the already-closed half.

**Affine half — the modification rule cannot reach it.** `μ₁ ↦ 22−μ₁` equals oddness
(centre 0) composed with translation by `44 = 2·22 = 2·2(l+h∨)` (level-README pin
`l+h∨=11`). A coroot-lattice translation is not in the finite `W(C₆)`, so it is outside the
modification rule's reach by genre. Confirmed structurally: `_rb_kt6.py` C3 finds the ONLY
finite reflection centre of `Φ_k` is `0` (no centre at `22`), every `k`.

**Re-indexing at the term level — REFUTED (witnesses).**
- Naive first-coordinate modification `Xc(P₁, tail) ↦ ±Xc(44−P₁, tail)` (`_rb_kt5.py`) sends
  the LOAD-BEARING (large) coefficients onto ZEROS. Witness `k=0, p=12→32`: e.g.
  `Xc(12,(10,4,-2,-6,-8)) = −9272` reflects to `Xc(32,(10,4,-2,-6,-8)) = 0`; of the 32
  contributing terms only 3 reflect to a nonzero `(=±1)` value. No term correspondence.
- Multiset test (`_rb_kt6.py` C2, all 6 live pairs): the two partners' sums agree, but the
  multisets of `|Yc|` present DIFFER (e.g. `k=0 (12,32)`: 123 vs 7 nonzero `Yc`-terms, only
  the value `1` shared; `k=2 (20,24)`: 80 vs 53 terms). So no bijection pairing equal
  `Yc`-values can exist. The identity is **intrinsically summed** — reconfirmed from the
  modification-rule side, in the Yc coordinates, matching the summed-README `(m,c)` result.
- Consequence: the two sides do NOT march through "the same E-coefficients in reflected
  order". The ladders differ in length AND in the `Yc`-multiset they carry. Task (c)'s
  re-indexing is dead.

**The residual, sharpened (the one thing route 2 adds).** By §2 the whole two-row law is
equivalent to the single first-coordinate summed identity

```
Σ_{ε,r} σ(ε) Yc( m+2r+1 , τ(ε,r) )  =  Σ_{ε,r} σ(ε) Yc( 44−m+2r+1 , τ(ε,r) )   (in window),
```

with the tail-ladder `{(σ(ε), τ(ε,r))}` frozen. The tail is fully reduced by the finite
modification rule (done); what remains is a first-coordinate affine reflection of `Yc`
summed against that frozen ladder, i.e. the anti-period `Φ(p+44)=−Φ(p)` on the resolvent
tail. That is the affine (Weyl-Kac) node, route 1's object, not a modification rule.

**Why "22 = 21+1" is an AFFINE pivot, not a modification pivot.** The wall `μ₁+μ₂ ≤ 21` is a
support bound (level-README). The reflection centre `22 = wall+1 = 2(l+h∨)` carries the `+1 =
ρ`/imaginary-root shift of an affine node. A finite modification pivot would sit on a finite
wall `μ_i=±μ_j`; `22` sits on the affine wall. The smell was right that `22` is `wall+1`, but
that `+1` is affine, which reinforces the placement in route 1.

## §4 Domain fence (required check): l = 2 breaks

`_rb_kt6.py` C4, three-row `λ=(λ₁, λ₂, l)` with `λ₂≥l`, reflect `λ₁` (`μ₁↦22−μ₁`), restricted
to the clean symmetric window `m ∈ [2λ₂+12, 32−2λ₂]` where both partners exceed the tail-max:

- `l = 0,1,3,4`: reflection HOLDS (0 mismatches) at every reachable `λ₂`.
- `l = 2`: reflection BREAKS at every reachable `λ₂` (`λ₂=2: 4/4`, `λ₂=3: 2/3`, `λ₂=4: 2/2`).

This reproduces `_rb_schur_README.md`'s `l∈{0,1,3}` hold / `l=2` (ρ₄-tail resonance) break
via the first-coordinate read-off, so the fence is respected: the shared-tail reduction (§2)
is an identity that holds for all rows, while the reflection itself carries the `l=2` break,
exactly as it must.

## §5 Numeric checks (every one, with count)

| check | file | statement | result |
|---|---|---|---|
| V0 | `_rb_kt1` | live reflection pairs `F_k(m)=F_k(44−m)` | 6/6 |
| V1 | `_rb_kt1` | forward / backward telescope | 12/12 · 12/12 |
| V2 | `_rb_kt1` | ZY: full W(C₆)-antisym of `Yc` ≡ 0 | 20/20 |
| V3 | `_rb_kt1` | ladder `F_k(m)=Σ_r T_r` | 35/35 |
| C1 | `_rb_kt6` | shared-tail reduction (τ,σ m-independent) | 35/35 |
| C2 | `_rb_kt6` | term-bijection refute (sums equal, `|Yc|`-multisets differ) | 6/6 pairs |
| C3 | `_rb_kt6` | oddness holds; only finite reflection centre = 0 | 5/5 k |
| C4 | `_rb_kt6` | domain fence: l=2 breaks, else holds | as §4 |
| step5 | `_rb_kt5` | naive first-coord modification lands on zeros | witnessed |

## §6 Files

- `_rb_kt1.py` — base facts V0-V3 + ladder tables.
- `_rb_kt2.py` — flat `Yc`-term dictionaries per partner (shared-tail discovery).
- `_rb_kt3.py` — first-coordinate slices `H_r(p)`, wandering-centre palindrome test.
- `_rb_kt4.py` — the global `Φ_k(p)` (tiny support, oddness, reflection window, defect).
- `_rb_kt5.py` — Xc-level term dump + naive first-coord modification (refutation witness).
- `_rb_kt6.py` — consolidation C1-C4 (reduction, term-bijection refute, finite/affine, fence).

## §7 Net

Route 2 = the finite modification rule = the already-proved oddness; the open half is the
affine translation-44, which the modification rule cannot produce. The re-indexing of task
(c) is refuted term-by-term in every coordinate; the sole new asset is the shared-tail
first-coordinate reduction (§2), which pins the residual to a single-variable summed affine
identity and hands it cleanly to route 1 (A₁₂⁽²⁾ Macdonald / Weyl-Kac denominator, the
`sin(s_{1⁶})`-deletion resolvent). No claim here is un-computed; the summed-not-bijective
constraint held throughout (all reductions are at the summed / frozen-ladder level, never a
`(m,c)` term map).
