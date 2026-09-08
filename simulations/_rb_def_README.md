# The completion picture: F_k = Θ_k − W_k (theta minus its beyond-wall shadow)

Scratch, local-only (`_rb_def*.py`). Run from repo root, `PYTHONIOENCODING=utf-8`. All arithmetic
exact over ℤ. Consolidated gate: `python simulations/_rb_def_gate.py` (~1 min build, prints
`ALL GREEN`). Builds on `_rb_summed.py` (`Fk`, exact ℤ) and `f133_w_closed_form.build_halves`.

Context: the two-row reflection law `F_k(m) = F_k(44−m)` holds **in-window** `m ∈ [12+2k, 32−2k]`
(midpoint 22), live `k = 0..4`; the level-pin fixed the character as the mixed `(−,+)` dihedral
sector (`s: m↦−m` odd, `s₀: m↦44−m` even), translation `44ℤ` carrying `−1`; the `_rb_bnd` round
proved the finite law by five explicit palindromic `P_k` and showed the 44-relation is **not** an
array symmetry (it FAILS off-window, `26/225` global violations). This round tests the HYPOTHESIS
that the true theta-like object is a corrected array `F̂_k` completing the `(−,+)` structure
globally, with the correction supported exactly on the off-window defect.

## Headline — DECOMPOSITION HOLDS

The picture is confirmed, literally. Define the finite theta

> **Θ_k(t) = (t²² − t⁻²²)·P_k(t)**, with `P_k` the known window closed forms
> (`P_0 = −(1+q²)(1+q⁸)/q^{-10}` … `_rb_bnd3`), support = window ∪ −window only.

Then, exact over ℤ and gated 0/N in-wall (`m ∈ [−33,33]`, one anti-period; `F_k` support `|m|≤32`):

1. **Θ_k is an honest `(−,+)` theta**, GLOBALLY: its periodic extension `F̂_k` (the unique
   `(−,+)`-equivariant lift of the window seed) is odd about 0, even about 22, anti-period 44 with
   **0/4050** violations, and `F̂_k = Θ_k` on one period (**0/522**).
2. **F_k = Θ_k − W_k**, reconstruction **0/402**, with `W_k = Θ_k − F_k = −D_k`.
3. **D_k = F_k − Θ_k is the defect**: identically 0 in-window (**0/N**), supported strictly
   off-window in-wall (**0/N**), and equal to the affine-factorization defect = `R_k`.
4. **W_k is the WALL-SHADOW**: every `W_k` support point `m` has its `s₀`-mirror `44−m` strictly
   beyond the wall, where `F_k = 0` (verified all k, all points). So `F_k` = (palindromic theta) −
   (the part whose mirror falls past the wall) — exactly "theta minus its beyond-wall shadow".
5. **Three-way tie exact** (**0/18**): `Θ_k = S22·P_k`, `P_k = Θ_k/S22`, `F̂_k = Θ_k` on one period.
6. **Fence PASSES**: at three-row `l = 2` the completion FAILS (§(e)) — the picture is not an artifact.

## (a) The defect array D_k (exact ℤ, in-wall)

`D_k(m) = F_k(m) − Θ_k(m)`. In-window `D_k ≡ 0`; off-window support (both signs, `m`-units):

| k | window `[lo,hi]` | D_k support (m : value) |
|---|---|---|
| 0 | [12,32] | `{}` (defect-free) |
| 1 | [14,30] | `{−10:−1, 10:+1}` |
| 2 | [16,28] | `{−12:+1, 12:−1}` |
| 3 | [18,26] | `{−10:−1, 10:+1}` |
| 4 | [20,24] | `{−16:−1, 16:+1}` |
| 5 | [22,22] | `{−18:+1, −14:+3, 14:−3, 18:−1}` (vacuous window) |

The completed object `F̂_k` (infinite periodic lift) satisfies the full `(−,+)` structure with
0 violations; the raw `F_k` fails the window-reflection `30/402` (witness `k=0, m=−32:
F(−32)=1, F(44−(−32))=F(76)=0`). The completion is exactly the `(−,+)` symmetrization forced from
the window values; where the propagation lands off-window in-wall (below the lower edge / beyond the
wall) it predicts 0, and the difference from the actual `F_k` **is** `D_k`.

## (b) D_k ⟷ R_k, exact (signs + units)

`D_k` **is** the affine-factorization defect `R_k` (same computed object). Two exact renderings:

- **Script unit** `n = (36−m)/2 = 18 − m/2` (both signs of m → n):
  `R_1={13:+1,23:−1}`, `R_2={12:−1,24:+1}`, `R_3={13:+1,23:−1}`, `R_4={10:+1,26:−1}`,
  `R_5={9:−1,11:−3,25:+3,27:+1}`, `R_0={}`.
- **The READMEs' stated R_k strings** (`R₁=x−x¹³`, `R₂=x²−x¹²`, `R₃=x−x¹³`, `R₄=x¹⁰−x⁴`,
  `R₅=3x¹¹+x⁹−x⁵−3x³`) are reproduced **exactly, all k** (`_rb_def2` gate ALL MATCH) by

  > **R_k^stated(n) = (−1)^k · [ D_k folded to n ≡ (36−m)/2 (mod 22) ]**.

  The `mod 22` folds the negative-m partner into the `[0,22)` fundamental domain (the `44 = 2·22`
  anti-period); the `(−1)^k` is the `P_k`-orientation sign. This is the precise
  units/sign correspondence the arc asked for — no residual ambiguity.

## (c) The decomposition F_k = Θ_k − W_k, with W_k the wall-shadow

`W_k = Θ_k − F_k` (= `−D_k`), reconstruction `F_k = Θ_k − W_k` holds **0/402**. Each `W_k` point is a
wall-shadow — its `s₀`-mirror lies past the wall where `F_k` vanishes:

| k | W_k (m : value) | mirror check (44−m outside window, F_k(44−m)=0) |
|---|---|---|
| 1 | `{−10:+1, 10:−1}` | `54, 34` ✓ |
| 2 | `{−12:−1, 12:+1}` | `56, 32` ✓ |
| 3 | `{−10:+1, 10:−1}` | `54, 34` ✓ |
| 4 | `{−16:+1, 16:−1}` | `60, 28` ✓ |
| 5 | `{−18:−1,−14:−3,14:+3,18:+1}` | `62, 58, 30, 26` ✓ |

So the WHOLE array is exactly **theta minus its beyond-wall shadow**, both parts explicit finite
ℤ-arrays. `Θ_k` is the palindromic core (the honest theta characteristic); `W_k` is precisely the
below-window content whose reflection would require support past the proven wall, so it cannot be
palindromized into the window.

## (d) Θ_k generating form and the three-way tie

`Θ_k = S22·P_k` with `P_k` the known factored window polynomials, recentered about 22 (`b = m−22`),
each self-reciprocal `P_k(1/t)=P_k(t)`:

```
P_0 = −(t⁴+1)(t¹⁶+1)/t¹⁰            P_1 = (t−1)²(t+1)²(t²+1)²(t⁴+1)²/t⁸
P_2 = 3(t⁴+1)/t²                    P_3 = −(t⁸+3t⁴+1)/t⁴
P_4 = (t⁴+1)/t²                     P_5 = 0
```
(In the `q = t`, `j = (m−12)/2` coordinate of `_rb_bnd3` these read `P_0=−(1+q²)(1+q⁸)`,
`P_1=q(1−q⁴)²`, `P_2=3q⁴(1+q²)`, `P_3=−q³(1+3q²+q⁴)`, `P_4=q⁴(1+q²)`.)
Tie (gate 0/18): `Θ_k/S22 = P_k` exactly (S22 divides Θ_k), and `F̂_k = Θ_k` on one period. So
`P_k` (the affine-factorization quotient), `Θ_k` (the finite theta), and `F̂_k` (the `(−,+)` lift)
are one object in three coordinates.

**On the "why".** The mechanism "beyond-wall shadow = actual defect" is **not independently
derived** by the completion: `Θ_k = S22·P_k` has NO support past the wall (its support is exactly
the window), so there is no separate beyond-wall continuation to match. The statement "`W_k` is a
wall-shadow" is therefore **equivalent to**, not a proof of, the in-window palindromy itself — i.e.
the same deleted-top-weight theta/affine-C₆ identity the level-pin and `_rb_bnd` rounds isolated,
now written with BOTH sides as explicit finite ℤ-arrays.

## (e) Fence: the completion FAILS at three-row l = 2

Same test as a harness (`_rb_def3`): per k, does a reflection center `c` admit a split into a
palindromic CORE (≥1 genuine reflected pair) plus DEFECT points that are all wall-shadows? A valid
picture needs a **single uniform c across the live rows**.

- **Positive control, two-row `l=0`**: clean centers `k0={22,30}, k1={22}, k2={22}, k3={22},
  k4={20,22}`, `k5=vacuous`; **intersection = {22}** — uniform, as the theory demands. ✓
- **Fence, three-row `l=2`**: `k0={22}, k2=NONE, k3={18}, k4=NONE, k5={16}` (`k1` empty). Rows
  `k=2,4` carry **substantial 4-point support yet admit NO clean completion center at all**; the
  rows that do complete **disagree** (22 / 18 / 16). **No uniform center exists.** ✓ FENCE PASSES.

The `l=2` completion is not merely shifted — it is structurally impossible for half the rows. The
two-row picture is therefore genuine, not an artifact of the construction.

## Verdict

**DECOMPOSITION HOLDS.** Exact statement: `F_k = Θ_k − W_k` with `Θ_k = (t²²−t⁻²²)·P_k` the honest
`(−,+)` theta (palindromic core) and `W_k` the beyond-wall shadow, both explicit finite ℤ-arrays;
gated 0/N on structure, reconstruction, shadow, tie, and the stated-`R_k` correspondence; fenced at
`l=2`.

**What it implies for the proof burden.** The completion **reframes but does not discharge** the
residual. It converts the open two-row identity into a single, fully explicit equation —
`W_k = Θ_k − F_k` is a wall-shadow, equivalently `F_k` is the wall-truncation of the explicit theta
`Θ_k` — with every term on both sides a known integer. This is a concrete target (an identity
between two given finite arrays), not a structural mystery, and it is **consistent with `_rb_bnd`
L3/L4**: it furnishes NO new array symmetry (the 44 is still not a symmetry of the coefficient
array; `Θ_k` is built from the window values, not from a symmetry of `X`). The surviving open core
is unchanged in kind — the deleted-top-weight theta/affine-C₆ character identity — but now stated
with both the theta part and the shadow part pinned as explicit closed forms, which is the useful
handle for the `n_λ` hunt.

## Files
- `_rb_def1.py` — D_k construction, F̂ global-structure gate, in-window vanishing, three-way tie (d).
- `_rb_def2.py` — stated-R_k correspondence gate (b'), first l=2 support/reflection probe.
- `_rb_def3.py` — sharpened fence (e): completion harness, two-row control (uniform 22) vs l=2 (fails).
- `_rb_def_gate.py` — consolidated `ALL GREEN` gate (a, b', c, d, e).
