# Attack B: the half-characteristic theta floor — the [½,½] pin, the dead single product,
# the NEW telescope law, and the closed resonance side-question

Scratch, local-only (2026-07-21, `_rb_theta{1..5}.py`). Run from repo root,
`PYTHONIOENCODING=utf-8`; scouts 1/3 need only the dem f-cache + the committed two-row
rows (cross-checked against each other before use); scouts 2/4/5 are sympy-exact on the
committed closed forms / the committed 143-term `chiC_coeffs.txt`. Context: the F134 plan's
Attack B (memory `f127_closed_form_find`), the Attack-A deliverable
`_rb_planA_sidebyside.md` (the ρ-condition, the 3-part fence), `_rb_def_README.md`
(Θ-decomposition), `_rb_level_README.md` (the (−,+) pin).

## Verdict up front

1. **Step 1 DONE — the characteristic is pinned: [a,b] = [½,½], i.e. θ₁, the odd one.**
2. **Step 2 — a SINGLE product identity is dead** (N(q,z) irreducible, rank 4 ⟹ ≥4 theta
   summands); the object is 8 sparse integers on a near-bijective (b,k) inventory.
3. **NEW LAW (the round's positive find): the two-row telescope** `N(q,q) = −1`, i.e.
   `Σ_{j+k=d} n_(j,k) = −δ_{d,0}` — a SECOND exact law living on exactly the two-row
   sector, transverse to the reflection, whose three-row defect is ODD about d = 11.
4. **Side-question CLOSED: the l=2 break is NOT a Conway–Jones angle resonance**
   (full-atlas scan, zero non-elementary vanishing sums).

## T1 — the spectral pin (`_rb_theta1.py`, gates: cross-check exact, DFT 40 dps vs sympy-exact)

The (−,+) sector on the μ₁-lattice (odd about 0, even about 11, anti-period 22) is spanned
EXACTLY by the half-integer frequencies `b_r(μ₁) = sin((2r+1)πμ₁/22)`, r = 0..10
(orthogonality on ℤ/44 verified to 6e-39). The lift F̂_k = Θ_k expands as

> **F̂_k(2μ₁) = (1/11) Σ_r A_{k,r} · sin((2r+1)πμ₁/22)**, with EXACT amplitudes
> **A_{k,r} = (−1)^r · P̃_k(u_r)**, u_r = 2cos((2r+1)π/22),
> P̃₀ = −u(u⁴−4u²+2) = −C₂C₈ | P̃₁ = u²(u²−4) | P̃₂ = 3u | P̃₃ = −(u²+1) | P̃₄ = u.

(55 pairs verified numeric-vs-exact to 1.3e-39; reconstruction on all 44 points exact.)
Half-integer frequency = characteristic a = ½; the alternating sign (−1)^r = b = ½. That is
literally the z-series shape of **θ₁ = θ[½,½]**, the ODD characteristic — the unique one of
the four classical thetas vanishing identically at the origin, the algebraic twin of the
deleted 1⁶ sheet (Tom's −½ | 0 | +½, now with the classical label attached). The lift is a
truncated θ₁ in z = μ₁/22 whose q-power profile is replaced by the finite profile P̃_k(u_r).

**The central Niven frequency r = 5 (angle π/2, u = 0): A_{k,5} = (0,0,0,1,0)** — every
column dies at the center EXCEPT k = 3, and P̃₃ = −(u²+1) is also the one column that is
NOT a C/S-product. k = 3 is the resonant column of the spectrum.

Frequency-axis parity (u ↦ −u under r ↦ 10−r): A symmetric for k odd, antisymmetric for
k even — the u-parity of P̃_k is (−1)^{k+1}.

## T2 — the product hunt (`_rb_theta2.py`, `_rb_theta4.py`, all sympy-exact)

- `N(q,z) = Σ_k P_k(q) z^k` is **irreducible over ℤ**; likewise under z ↦ −z, q-reciprocal
  (q¹⁰N(1/q,z) = N — the law itself), (q,z)-swap, and every regauge z ↦ q^s w, s ∈ [−4,4].
- Coefficient matrix n_(j,k) (11×5) has **rank 4** ⟹ no single theta product (rank 1); any
  theta-decomposition needs ≥ 4 summands. Only column relation: P₂ = 3P₄.
- Multiplier collapse (truncated-geometric signature): none of 17 candidate binomials
  (q±z, 1±qz, q²±z, …) makes (mult)·N factor or sparsify. The best compressions
  (q+z)·N and (1+qz)·N: 17 terms vs 14 — nothing.
- The sparse inventory (b-units, C_a = t^a+t^{−a}): k=0: −C₁₀−C₆ | k=1: C₈−2 | k=2: 3C₂ |
  k=3: −C₄−3 | k=4: C₂. Eight integers total; each frequency b ∈ {2,…,10} owned by exactly
  one k except b=2 (k=2,4, ratio 3); constants only at odd k. The **window-edge value** is
  (−1)^{k+1} for k ≠ 2 and 0 at k = 2 — the edge dies exactly once, at k = 2.

## T3 — the NEW telescope law (`_rb_theta4.py` T4c find, `_rb_theta5.py` reach)

> **N(q,q) = −1 exactly**, i.e. the anti-diagonal sum rule
> **Σ_{j+k=d} n_(j,k) = −δ_{d,0}** (verified d = 0..10 by hand and by machine).

Reach (committed 143-term table):
- restricted to rowcount ≤ 2: **HOLDS, 0 exceptions** — a genuinely TWO-ROW law;
- restricted to rowcount ≤ 3: FAILS at d ∈ {8,10,12,14} with defects **{−1, +2, −2, +1}**;
- the full 143-table: fails broadly (grand sum = 7, not −1).

**The three-row telescope defect is ODD about d = 11**: defect(d) = −defect(22−d) on all
four points. The s₀ center 11 reappears in the |λ|-grading of a DIFFERENT law. State of
affairs now: TWO exact laws (the reflection n\_(j,k) = n\_(10−j,k) and the telescope
Σ_{j+k=d} n = −δ_{d0}) live on exactly the two-row sector, both break at three rows, and
both breaks carry the center 11. The seam identity, whatever it is, should imply BOTH — a
much sharper target than the reflection alone. (A telescope is exactly what a truncated
geometric/theta object does at its ratio point; this is the first structure the truncated-θ₁
picture PREDICTED and the data then confirmed.)

Signed variants do not telescope ((−1)^{|λ|}: 7; (−1)^rows: −5; (−1)^{λ₂}: 55).

## T4 — the resonance side-question: CLOSED NEGATIVE (`_rb_theta3.py`)

Spectral form of a slice palindrome: f_τ symmetric about 11 ⟺ its antisymmetric part
h(b) = f(11+b) − f(11−b) has vanishing sine sums H_τ(ρ) = Σ_b h(b)·sin(ρπb/22). Scan of
ALL 54 break tails (tier census reproduced T1=39/T2=1/T3=14), all odd ρ ∈ [1,21], 40 dps,
candidates exact-verified:

- **No break anywhere in the atlas has a vanishing sum at any ρ ≠ 11.**
- Every ρ = 11 zero is ELEMENTARY: sin(11πb/22) = sin(πb/2) is rational (0, ±1) — parity
  zeros for even-offset defects, plus exactly two alternating-sum zeros (h(1) = h(3):
  tails (7,6,3,2,1) — the famous μ₂=7 l=2 break — and (7,6,5,3,2), both T1).
- The minimal T3 instance (4,4,2)↔(6,4,2) has single-pair defect h = {1: −6}: its spectrum
  is −6·sin(ρπ/22), nonzero at every lattice frequency — resonance-free by inspection.

So the l=2 break does NOT read as an angle resonance in the Conway–Jones genre; the
repo's native vanishing-cosine machinery has nothing to bite on here. Clean negative; the
plan's framing-note side-question is closed.

## Fence status (the Attack-A 3-part fence + the new material)

No candidate identity was minted this round, so nothing to fence yet. The fence GREW: a
future candidate must now reproduce (1) the 14 live reflection pairs, (2) truncation breaks
at μ₂ = 7,8, (3) the value-only break at (4,4,2), **(4) the telescope Σ_{j+k=d} n = −δ_{d0}
on two rows, and (5) the odd-about-11 three-row telescope defect {−1,+2,−2,+1}.**

## Hand-off (next moves inside Attack B)

- Hunt the **4-summand θ₁ decomposition** the rank forces (natural candidate gradings:
  the frequency inventory's b↔k near-bijection; the two window-edge diagonals).
- Try to PROVE the telescope law from the frozen-ladder/resolvent machinery (it is a
  single-variable identity Σ_k q^k P_k(q) = −1 over the committed closed forms — possibly
  elementary, and any proof mechanism is a candidate mechanism for the reflection too).
- The two-laws-one-structure angle: look for the object whose truncation yields BOTH laws
  with center 11 (the θ₁ ratio-point + reflection pair).

## ROUND 2 (same day): the telescope's anatomy — the corner interlock

Setting: g(μ₁,μ₂) := n_ext((μ₁,μ₂,4,3,2,1)) with n_ext(μ) = Σ_{ε∈{±1}⁶} sgn(ε)·Xc(ε∘2μ),
the alternant read-off at ARBITRARY integer μ (W(C₂)-antisymmetric in (μ₁,μ₂):
transposition and both sign flips, verified). Dominant-chamber support = EXACTLY the
committed 14-site two-row table, nothing outside the strip (`_rb_theta8.py`).

**The find (both families, exact):**

> **G(q,q) = G(q,1/q) = −2·q¹¹** where G(q,w) = Σ_{dominant} g·q^{μ₁}w^{μ₂}.
> Sum diagonals: Σ_{μ₁+μ₂ = c} g = −2δ_{c,11}. Difference diagonals:
> Σ_{μ₁−μ₂ = e} g = −2δ_{e,11}. (Rows/columns do NOT telescope — good contrast.)

**CORRECTED in round 3** (the first write-up wrongly placed both defects at one site):
the two defect lines are carried by the TWO ENDS of the extreme reflection pair — the
sum-defect line c = 11 contains ONLY the site (6,5) = λ(0,0) (value −2 = n_raw(0,0));
the diff-defect line e = 11 contains ONLY the wall corner (16,5) = λ(10,0) (value −2 =
n_raw(10,0)). Equality of the two defects ⟺ the extreme pair n\_(0,0) = n\_(10,0). The
sign flip μ₂ ↦ −μ₂ exchanges the two diagonal families and lands the corner's Weyl shadow
(16,−5) on the c = 11 line, balancing it in the extended plane (`_rb_theta7.py`: for
d ≥ 1 the sum-line support is exactly dominant chain ∪ mirror chain; only d = 0 carries
that extra pair).

**The corner interlock.** The telescope's two exceptions ARE the reflection law's extreme
pair, one end per diagonal family. The two laws are not merely parallel (same domain,
same center 11); they interlock at the wall. Also: the difference-diagonal e = 7 is
entirely EMPTY (a support hole: n\_(6,0) = n\_(7,1) = n\_(8,2) = n\_(9,3) = 0).

**Dead this round:** the naive wedge mechanism. The antisymmetric extension ĝ on [5,16]²
has **rank 8** (not 2 = one wedge h∧h′, not 4 = two wedges): the diagonal balance does not
come from a small-rank alternant of one-variable sequences. (A rank-2 wedge of palindromic
sequences would have implied the reflection outright — too good to be true, and it isn't.)

**Hand-off after round 2:** the seam identity's target statement is now
"the dominant two-row slice of S is diagonal-balanced in both W(C₂) families, with the
single defect −2q¹¹ at the wall corner, and s₀-palindromic in rows". One object, three
interlocking facts, all anchored at 11. The θ₁-characteristic (round 1) says WHICH modular
object should produce this; the corner says WHERE its boundary term lives. Natural next
probes: (i) the same anatomy one tail-step out (tail (5,3,2,1)-type / three-row slices:
does the corner interlock persist and where does l=2 sit in it); (ii) the Yc-resolvent
applied to the DIAGONAL generating function (the 𝟙-shift moves c by 12 — check whether
G(q,q) = −2q¹¹ follows from the telescope Xc = Σ_r Yc(·+(2r+1)𝟙) plus Yc's own diagonal
collapse); (iii) King §8's two-row ψ_{q,r} laws specialized to the corner pair.

## ROUND 3 (same day): the refinement — per-l anatomy and the rung resolution

**(i) One tail-step out (`_rb_theta9.py`, tails (μ₃,3,2,1), μ₃ = 4..8 ⟺ l = 0..4):**

| l | sum-diagonal defects (c: Σ) | diff-diagonal defects (e: Σ) |
|---|---|---|
| 0 | {11: −2} | {11: −2} |
| 1 | {16: −2, 18: −4, 20: −2} | {2: −2, 4: −4, 6: −2} |
| 2 | {15: +2, 17: +2, 19: +12} | {1: +2, 3: +8, 5: +6} |
| 3 | {18: −6, 20: −6} | {2: −6, 4: −6} |
| 4 | {19: +2, 21: +2} | {1: +2, 3: +2} |

- **The diagonal balance is STRICTLY two-row and INDEPENDENT of the reflection**: at
  l = 1 the reflection holds but both telescopes break broadly. The two laws share the
  domain boundary but not the mechanism gate — a very discriminating fence fact.
- **Pooled control**: converting per-l sum-defects to |λ| and summing reproduces T5c's
  rowcount ≤ 3 defect {0:−1, 8:−1, 10:+2, 12:−2, 14:+1} exactly, including a cross-l
  cancellation at |λ| = 6 (l=1's −1 against l=2's +1) invisible in the pooled view.
- **l = 2 is again the odd one**: it is the unique l where the sum- and diff-defect
  multisets DISAGREE ({2,2,12} vs {2,8,6}; every other l has matching multisets). A new,
  law-internal l=2 signature.

**(ii) Rung resolution along the Yc-resolvent (`_rb_theta10.py`):** writing
n_ext = Σ_r T_r with T_r(μ) = Σ_ε sgn(ε)·Yc(ε∘2μ + (2r+1)𝟙₆) (reconstruction gated 0/18):
individual rungs carry diagonal defects up to ±736 (rungs r = 3..9 live), canceling across
r to exactly {11: −2} in BOTH families — **ONLY-FULL-SUM again**, the same verdict the
old law's ladder produced. Per-rung sum-diagonal defects alternate in c; carrier profiles:
λ(0,0): {3:−5, 4:9, 5:−19, 6:26, 7:−14, 8:1}, λ(10,0): {3:−3, 4:19, 5:−44, 6:44, 7:−18},
each summing to −2. No rung-level mechanism exists.

**Fence part (6), new:** a candidate seam identity must make the diagonal telescopes
break at l = 1 WHILE the reflection holds there — any mechanism deriving both laws from
one symmetry of the slice is thereby wrong; the joint object must gate them differently.

## ROUND 4 (same day): the sgn-kernel lemma, the genre naming, universality refuted,
## and the recursive palindrome

**(a) The sgn-kernel lemma (exact algebra, no computation needed).** For the
W(C₂)-antisymmetric extension ĝ, the dominant diagonal sums are odd-kernel line sums:

> D⁺_c = ½·Σ_{a+b=c} sgn(a−b)·ĝ(a,b)   and   D⁻_e = ½·Σ_{a−b=e} sgn(a+b)·ĝ(a,b)

(each line's involution — transposition resp. anti-transposition — has sign −1 and only
zero-value fixed points, so the dominant restriction IS the sgn kernel). The laws say
these odd-kernel sums collapse to single boundary monomials −2q¹¹. Theta-type object ×
half-lattice sgn kernel with a single boundary term: that is the signature of the
**partial-theta / Appell–Lerch genre** (Zwegers' μ: two elliptic variables, antisymmetric,
modular up to an explicit boundary/shadow correction). GENRE HYPOTHESIS, banked: the seam
object is a truncated Appell–Lerch / partial theta on the (B₆,C₆) seam; the two −2q¹¹
defects (at the extreme-pair ends) are its shadow terms. This refines round 1's [½,½] pin
(θ₁ is exactly the theta whose Appell–Lerch companions carry μ-type shadows).

**(b) Pair-universality REFUTED (`_rb_theta11.py`).** The pair (μ₂,μ₃) with staircase
complement {μ₁} ∪ (3,2,1) does NOT balance: μ₁ = 12 breaks broadly (4 sum-defect lines),
μ₁ = 13..16 have single-site support (trivially defective), μ₁ = 17 empty. So the
diagonal balance is NOT a pair/staircase phenomenon of the alternant; it is irreducibly a
**two-free-rows** phenomenon — the same exclusivity as the reflection. Both laws live on
exactly the two-row sector; the ρ_{C₄}-tail (4,3,2,1) is the two-row sector's complement
(the balance held exactly there in `_rb_theta9.py` — the ρ/half-sum condition again).

**(c) The recursive palindrome (from the scout-9 defect tables, hand-verified).** The
per-l telescope-defect profiles are THEMSELVES palindromic for every live l ≠ 2, in both
families: l=1: (−2,−4,−2) about c*=18 / e*=4; l=3: (−6,−6) about 19 / 3; l=4: (2,2)
about 20 / 2. The centers pair as **c* + e* = 22 in every case** ⟺ the defect profile's
μ₁-centroid is pinned at the WALL 11. At l = 2 both palindromy ((2,2,12) / (2,8,6)) and
the centroid (21.75) break. So the reflection symmetry RECURSES into the defects of the
second law — with the same center and the same lone resonance. The structure is
self-similar one level down, and l=2 kills it at both levels.

**Still open from the round-3 hand-off:** King §8 ψ_{q,r} at the corner pair (needs the
paper's φ/ψ/a-definitions instantiated; the kw round's verdict "genre match, no closed
identity" stands meanwhile).

## ROUND 5 (same day): the determinacy audit — laws + three mirror anchors = everything

`_rb_theta12.py` (sympy-exact linear algebra, all gates green). Unknowns: the 21
parity-allowed dominant strip cells x\_(j,k). Constraint families, all discovered and
verified earlier: row palindromy (9), sum-diagonal balance d = 2..10 (5), diff-diagonal
balance e = 0..8 (5) [the two defect cells (0,0), (10,0) exempt], plus the natural extra
facts: column relation x\_(j,2) = 3x\_(j,4), the e = 7 support holes, the k = 5 vacuity.

> **Result: rank 18 of 21 — the laws pin the two-row sector to a 3-dimensional residual
> space, and the three residual directions are resolved by the MIRROR-FIXED anchors**
> n\_(0,0) = −1 (empty weight), n\_(5,1) = −2, n\_(5,3) = −3 (the s₀-fixed center values
> of the odd-k rows = exactly the constants of the C-inventory). Anchor matrix det = −2;
> coefficients (−3, −1, 1); reconstruction exact on all 21 cells.

**THE DETERMINACY STATEMENT.** {F134 row palindromy + both diagonal balances + column
relation + support holes + vacuity} ⊕ {the three values on the mirror-fixed locus} =
the ENTIRE two-row sector. The seam identity's remaining quantitative content is exactly
three integers, −1, −2, −3, all living ON the mirror (the s₀-fixed points and the empty
weight). Any candidate mechanism that produces the law families and evaluates correctly
at the fixed-point locus produces the whole table — the hunt's target is now minimal and
explicit. (Interpretive echo, label only: the laws are the pair of mirrors; what they
cannot generate is the seed sitting on the mirror itself — the fixed-point data.)

## ROUND 6 (same day): the anchor hunt — the wall column, the half-shift ladder, and
## the Λ^{h∨} address

**(a) The spectral consistency loop (hand algebra, closed but circular).** The anchors are
the wall values F_k(22) = (1/11)Σ_r P̃_k(u_r); power sums over the odd cosine lattice are
central binomials (Σ_r u_r^p = 11·C(p,p/2) for even p ≤ 10, else 0), giving P̃₁ ↦ 6−8 = −2,
P̃₃ ↦ −(2+1) = −3 ✓ — equivalently just CT(P_k). Consistent, but reads the table.

**(b) The wall inventory (`_rb_theta13.py`, one off-by-one in the tail map fixed:
ρ-tail is (5,4,3,2,1)).** The wall function w(τ) = f_τ(11) has 34 nonzero entries =
EXACTLY the λ₁ = 5 column of the committed table (cross-check 34/34, w = 2n). Values in
[−10,12] (n: [−5,6]). **Parity selection: every wall entry has |ν| ODD** (ν = λ minus the
leading 5) — opposite to the |λ|-even parity of the table itself.

**(c) The rank ladder carries the half-shift (`_rb_theta14.py`).** The rank-5 analog
X₅ = Δ̂₅·∏₁₅ has every coordinate in 4+15 = 19 factors ⟹ ALL exponents odd: X₅ lives on
the HALF-SHIFTED (spin/metaplectic) weight lattice, μ ∈ (ℤ+½)⁵. The rank ladder
alternates integer (even rank, 36 factors even) / half-integer (odd rank): the wall
column's odd-|ν| parity is the ½-shift made visible. (Tom's −½ | 0 | +½ as a RANK
phenomenon.) **The naive descent is REFUTED on data**: with either half-shift convention
the rank-5 read-off table is tiny (12 resp. 1 entries, all ±2) and matches neither the
wall support (34) nor constant ratios (−2, 3, 6, 4 at the four overlaps) — the wall
column is NOT the rank-5 table. (Convention caveat: the metaplectic read-off normalization
was guessed; but no convention fixes a 34-vs-12 support mismatch.)

**(d) THE DERIVED ADDRESS: the wall is the Λ^{h∨} slice.** Collecting t₁-degree in
X₆ = Δ̂₆·∏₃₁ (elementary algebra, no conjecture): every one of the 36 t₁-carrying factors
(5 differences + 31 sheets) contributes t₁^{±1}, so the coefficient of t₁^m factors as
Δ̂₅ × (a twisted e_{(36−m)/2}-sum over the 36-letter system {−e_v, v=2..6} ∪
{L′ ∈ {±1}⁵ ∖ 1⁵} = the weights of spin⊕vector of the x′-side, top deleted). At the wall
m = 22: **(36−22)/2 = 7 = h∨(C₆)** — the wall column is the co-Coxeter-depth slice, and
the three anchors are coefficients of a twisted **Λ^{h∨(C₆)}** character of the
spin⊕vector letter system paired against Δ̂₅. First Lie-theoretic address for the
anchor data; the h∨ that pins the level (11 = ℓ+h∨) reappears as the exterior-power
depth of the mirror-fixed seed. VERIFICATION of the slice identity + the Λ⁷ readout =
round 7's opening move.

## ROUND 7 (same day): the Λ^{h∨} slice VERIFIED, the anchors' grade profile, and the
## deletion cascade to the rank floor

**(a) The slice identity is exact (`_rb_theta15.py`, all gates green).** With letters
a_i (differences → t_v^{−1}, sheets → t′^{L′}): [t₁^{±22}]X₆ = −Δ̂₅·t′^{∓2𝟙}·e₇(a^{∓2});
verified against Xc on 50 sampled exponents (0 mismatches) and the WHOLE wall column
reproduced through it: **34/34**. The a^{−2} letter system = {2e_v} ∪ {2M: M ∈ {±1}⁵ ∖
{−1⁵}} = the (doubled) weights of vector ⊕ spin of the x′-side, bottom spin weight
deleted.

**(b) The anchors' grade profile** (contribution of Λ^j V ⊗ Λ^{7−j} Spin): anchor −2
(ν=(1)): {j=0: −4, 1: −2, 2: +4, 3: 0, 4: −2, 5: 0} summing to −4 = 2n; anchor −3
(ν=(3)): {0: −6, 1: −2, 2: +2, rest 0} summing to −6. Spread over grades — the seed is
not a single exterior-sector object.

**(c) THE DELETION CASCADE (`_rb_theta16.py`).** The FULL (bottom-restored, 37-letter)
system has ZERO wall content at Λ⁷ and again at Λ⁶ — so the wall falls twice:

> wall = −readoff[Δ̂₅·t^{−2𝟙}·e₇(36)] = +readoff[Δ̂₅·t^{−4𝟙}·e₆(36)]
>      = −readoff[Δ̂₅·t^{−6𝟙}·e₅(36)], each exact **34/34**,
> and the cascade FLOORS at Λ⁵: there the full 37-letter system carries nonzero content
> at every wall ν (34/34) — descent below k = 5 is not exact.

**Reading: 7 = h∨(C₆) → 6 → 5 = n (the rank).** The mirror-fixed seed has three exact
dresses (Λ^{h∨}, Λ^{h∨−1}, Λ^{rank} at bottom-shifts −2𝟙, −4𝟙, −6𝟙), and the content
irreducibly sits at the RANK exterior power. The rank-6 pattern "the full object has zero
content, everything lives in the deletion" (Φ_Y ≡ 0, kw round) recurses at the finite
Λ-level — twice — and stops exactly at Λ^n. The anchor hunt's literature door is now:
decompose Λ⁵(spin ⊕ vector of so(11)-type letters, bottom deleted) against Sp(10) with
the −6𝟙 (= −3× bottom) twist.

## ROUND 8 (same day): the wall character in closed form, and the obstruction anatomy

`_rb_theta17.py` (full spectral tabulations over ν₁ ≤ 13, ≤ 5 rows; split gate OK).

**(a) THE CLOSED FORM (the milestone).** The Λ⁵ floor object's Sp(10)-spectrum has
support EXACTLY on the 34 wall entries — no extras anywhere:

> **Θ_wall := Σ_ν n\_(5,ν)·χ^{Sp(10)}_ν = −½ × [C₅-expansion of Δ̂₅·t′^{−6𝟙}·e₅(V⊕Δ∖bottom)]**

(readoff = −2n on all 34, zero on every other dominant ν). The ENTIRE mirror column —
both anchors −2, −3 included — is now a single explicit finite object: the rank-th
elementary symmetric function of the 36 vector⊕spin letters, triple-bottom-twisted, read
against the A₄ half-angle Vandermonde. First closed form the seed data has ever had,
fully decoupled from the 16.3M-term X. (The third anchor −1 = n\_(0,0) is the telescope
boundary, interlocked separately.)

**(b) The obstruction anatomy.** O5 = expansion of Δ̂₅t^{−6𝟙}e₅(37 letters) and the
Λ⁴ remainder L4 = expansion of Δ̂₅t^{−8𝟙}e₄(36): both have the SAME 68-point support
(all |ν| ODD, all ν₁ ≤ 5), values all even, and **off the wall they cancel exactly**
(O5 = L4 on the 34 non-wall ν; 68 = 34 wall + 34 cancellation partners; split
F5 = O5 − L4 gated). So the descent below Λ⁵ fails not by absence but by a matched
pair: the full-system content and the Λ⁴ term agree everywhere except on the wall,
where their difference IS the wall character. Naming O5 (recognition as a known
Sp(10)/so(11) decomposition) stays open — the concrete literature question is now:
what is Λ⁵(spin ⊕ vector) of the so(11)-letter system as an Sp(10) virtual character
with the −3·bottom twist?

## ROUND 9 (2026-07-21 next session): THE ANCHORS HAND-DERIVED — the closed form
## collapses into Racah blocks of the spin cube

Goal (handover item (a)): derive n\_(5,1) = −2 and n\_(5,3) = −3 from the round-8
closed form Θ_wall = −½·[Δ̂₅·t′^{−6𝟙}·e₅(V⊕Δ∖bottom)] by hand. Result: the whole
closed form COLLAPSES into a finite block decomposition over the FULL spin cube,
via three hand-provable lemmas (scouts `_rb_theta{18..22}.py`; every step gated).

**Units.** Work in x = t² (the doubled theta17 lattice halved): Δ̂₅ = a_δ(x), the
plain A₄ alternant at δ = (2,1,0,−1,−2); letters e_v (vector) and M ∈ {±1}⁵∖{−𝟙}
(spin31); shift x^{−3𝟙}; readoff = W(C₅)-antisymmetrized coefficient at μ = ν+ρ₅.
UNIT TRAP (cost one wrong round): Δ̂₅ built from t-offsets e_u−e_v is a_δ only in
x-units; mixing lattices makes the Pieri check fail spuriously.

**L1 (Pieri-prefix, 3-line proof).** a_δ(x)·e_a(x₁..x₅) = a_{ζ_a} with ζ_a = δ +
(1^a, 0^{5−a}) — exactly ONE term, because δ's entries are consecutive integers:
adding 1 to a non-prefix subset collides two entries and kills the alternant.
(Gate: scout 21 Q1, OK for a = 0..5.)

**L2 (bottom recursion).** e_b(spin31) = Σ_j (−1)^j x^{−j𝟙} e_{b−j}(spin32); at the
coefficient level c31_b(w) = Σ_j (−1)^j c32_{b−j}(w + j𝟙) — the shift is **+j𝟙**
(SIGN TRAP: with −j𝟙 one gets a fake near-miss: every block except (0,0) vanishes
and the "anchors" come out −3/−4 — one off. Scout 20 is that wrong turn, kept.)

**L3 (flip absorption + stabilizer kill).** c32_m is W(C₅)-invariant (full-cube
flip/permutation symmetry), so the double sum over (ε, σ) absorbs into ONE signed
sum over w ∈ W(C₅):

> **n\_(5,ν) = −Σ_{a,j} (−1)^j·I(a,j)**, I(a,j) = Σ_{w∈W(C₅)} sgn(w)·c32_m(μ − w(ζ̃)),
> ζ̃ = ζ_a − (3+j)𝟙, m = 5−a−j;

and I(a,j) = 0 whenever ζ̃ has a zero entry or two entries of equal absolute value
(an odd reflection stabilizes ζ̃; pair w ↔ w·refl). The zero-entry kill fires
exactly when 3+j ∈ ζ_a — in particular every a ≥ 1, j = 0 block dies (ζ_a's top
entry is 3). Support windows (|arg_i| ≤ m, arg_i ≡ m mod 2) kill most of the rest.

**Gate (scout 21): the block formula reproduces the ENTIRE wall column 34/34.**

**The anchor evaluations (scout 22 prints the hand-checkable inventories):**

| ν | live blocks (a,j) | I values | contributions −(−1)^j·I | total |
|---|---|---|---|---|
| (1) | (0,0),(0,1),(0,2),(1,1),(2,1),(4,1) | −3,−7,−2,−1,+2,−1 | +3,−7,+2,−1,+2,−1 | **−2** |
| (3) | (0,0),(0,1),(0,2),(1,1),(2,1) | −4,−10,−3,−1,+1 | +4,−10,+3,−1,+1 | **−3** |

The (4,1) block of ν=(1) is literally ONE term: m = 0 forces μ = w(ζ̃) with
|ζ̃| = (1,2,3,4,6) = the multiset of μ = (6,4,3,2,1); the unique w is the
reversal (even) with all five signs flipped (sgn = −1), so I = −1 by inspection.
The m ≤ 3 blocks have 1–6 W-terms each; the two m=5 head blocks have 61/37 terms
in 12/9 |arg|-classes (scout 22 lists them; the c₅ class values 752, 355, 162, 70,
64, 28, 27, 10, 4, 3, 1 are cube-subset counts, each hand-derivable by the
coordinate-5 recursion c₅(5,v) = c₅ of the 4-cube).

**Reading (the structural yield).** The head block (a=0, j=0) has ζ̃ =
(−1,−2,−3,−4,−5) = −(6𝟙−ρ₅): it is the honest Racah/Brauer alternating sum, i.e.
**the Sp(10) decomposition multiplicity of Λ⁵(spin cube)**. The j ≥ 1 blocks are
the same sums at 𝟙-shifted arguments (the bottom-deletion echo: "false-level"
readings, the partial-theta flavor again); the a ≥ 1 blocks are the vector-letter
Pieri corrections. So the wall column = Λ^rank(Δ)-multiplicities PLUS a finite
tail of 𝟙-shifted and prefix-shifted corrections. HONEST QUANTIFIER (checked over
all 34): the head alone explains 0/34 — head and tail are CO-LEADING, mostly
opposite-signed (e.g. ν=(1): head +3 vs n = −2). The decomposition is exact and
finite, not head-dominated; the O5-naming question (open (b)) must name the tail
blocks too, not just recognize the head.

## ROUND 10 (same session): THE TAIL NAMED + THE LOOP REDUCED TO ONE FUNCTIONAL —
## the universal block formula, the Core/P split, and the even-part annihilation

Goal (handover items (b) + (c)): name the shifted tail; close the loop to the laws.
Scouts `_rb_theta{23..26}.py`, all gated.

**(1) The universal block formula (the whole two-row table, not just the wall).**
The general slice is elementary: [t₁^{2μ₁}]X₆ = (−1)^d·Δ̂₅·t′^{−2𝟙}·e_d(a^{−2}),
d = 18−μ₁ (gated vs Xc at μ₁ = 6, 9, 13, 16: 0 mismatches). The round-9 collapse
then gives, with κ = 1+j (twist x^{−𝟙} in x-units):

> n\_(λ₁,ν) = (−1)^d Σ_{a,κ} (−1)^{κ−1} Σ_{w∈W(C₅)} sgn(w)·c32_{m}(μ − w(ζ_a − κ𝟙)),
> m = d+1−a−κ, d = 12−λ₁ — verified on the ENTIRE committed two-row table (14/14
> live sites). The (a,κ,m) blocks are dress-independent (e₇/e₆/e₅ give the same
> blocks: the deletion cascade is manifest).

**(2) THE NAME.** Dominantizing ζ_a − κ𝟙 gives ξ = ρ₅ + (κ−3)𝟙 − (0^{5−a},1^a), so
each block is an honest tensor multiplicity with highest weight η(a,κ) =
((κ−3)^{5−a}, (κ−4)^a) = (κ−4)ω₅ + ω_{5−a}:

> **n\_(λ₁,ν) = Σ_{a,κ} (−1)^κ·[V_{(κ−4)ω₅+ω_{5−a}} ⊗ Λ^{m}(Δ) : V_ν]·(−1)^{d+1}**,
> Δ = the 32-weight spin cube. The ω₅-twist ladder = the false-level/partial-theta
> corrections; the Euler-characteristic shape (Σ (−1)^j twist^j ⊗ Λ^{b−j}) is a
> twisted Koszul complex on the spinor orbit. Dead blocks (κ ∈ ζ_a) = non-dominant
> η — the stabilizer kill and the dominance wall are the same fact.

**(3) Resummations that do NOT split the law** (two more ONLY-FULL-SUM verdicts):
per-(a,κ) one-row palindromy FAILS (scout 23 Q3); the r-diagonal (a+j = r, the
virtual system Ξ = {e_v} ⊖ {−𝟙}, e_r(Ξ) = Σ_{a+j=r}(−1)^j e_a(V)x^{−j𝟙}) kills
r = 0,1 identically but per-r palindromy FAILS too (scout 24). The reflection
refuses every per-block split tried so far — consistent with the arc's rung/comb/
constituent history.

**(4) THE CORE/P SPLIT (the loop's new floor; scout 25/26).** Λ36 = Core ⊔ P with
**Core = cube ∖ {±𝟙}** (30 letters, negation-closed, ∏ = x⁰) and **P = {e₁..e₅, 𝟙}**
(6 letters). Channels T_i(g) := R⁻[e_i(P)·e_g(Core)]. Two EXACT mirrors, both
3-line hand proofs:
- **(M1) Core mirror**: e_g(C) = e_{30−g}(C) as polynomials ⟹ T_i(g) = T_i(30−g);
- **(M2) P-skew**: e_{6−i}(P) = x^{2𝟙}e_i(P)(x⁻¹) + readoff inversion-oddness ⟹
  **T_{6−i} = −T_i**; in particular T₃ ≡ 0, channels i = 0,1,2 only.
Then c_d = Σ_{i≤2}[T_i(d−i) − T_i(d−6+i)] (gate 14/14), the trivial s-oddness is
automatic, and with S_i(g) := T_i(g) + T_i(8−g) (= T_i(g) + T_i(g+22) via M1 — the
22-TRANSLATION LIVES IN CORE DEPTH):

> **c_d − c_{14−d} = A(ν,d) := Σ_{i≤2}[S_i(d−i) − S_i(d−6+i)]** exactly, and
> **c_d = ½·Σ_{i≤2}[O_i(d−i) − O_i(d−6+i)]** with O_i(g) = T_i(g) − T_i(8−g)
> (gate 14/14): the table lives entirely in the ODD-about-4 channel parts.

**(5) THE SEAM STATEMENT (what remains to prove).** A(ν,d) = 0 on the F134 domain;
measured support of A (scout 26 G4): one-row k=0,1,3 all zero; k=2 only at the
out-of-window/non-partition edge d ∈ {2,12}; k=5 (the vacuous slice) nonzero;
l=1 (ν=(k,1)) ALL ZERO (the reflection's l=1 hold, now a channel fact); l=2
(ν=(2,2),(3,2)) NONZERO in-window (the famous breaks); l=3/l=4 nonzero only at
the d ∈ {2,12} edge. The whole break atlas is the support of ONE functional on
the even-about-4 channel parts of a negation-closed 30-letter core read through
a 6-letter window. The mechanism question is now: why does the P-window
functional annihilate the even parts exactly on ≤1-free-row ν (+ the l=1 slice)?

Warnings: (i) the earlier "V_i-defect" formulation (inline, before scout 26) had
a pairing error (16+2i ≠ 22); the correct pairing hits SUMS S_i, not differences;
(ii) T_i(g) for g > 16 must be folded by (M1) before lookup, e_g tables were only
built to g = 16.

## ROUND 11 (same session): THE SEAM IDENTITY — the wall is a Chebyshev divisor,
## the law is a remainder-degree bound, the committed P_k are the quotients

The annihilation step, executed to the end. Scouts `_rb_theta{27,28,29}.py`.

**(1) The cosine matrix (scout 27).** The Core = 15 negation pairs; per pair the
e-series factors as (1+z²) + z·C_M (empty/double carry weight zero). Hence exactly
E_C(z) = Σ_t z^t(1+z²)^{15−t} e_t({C_M}) and

> c_d = Σ_{i,t} ψ_{i,t}·Bin(15−t, (d−i−t)/2),  ψ_{i,t} := R⁻[e_i(P)·e_t({C_M})]

— the whole depth dependence is a binomial transform of a FINITE matrix of
cosine-pair readoffs (gate 14/14; skew ψ_{6−i} = −ψ_i, ψ₃ = 0).

**(2) The one-variable collapse (scout 28).** With y = z+z⁻¹: z^t(1+z²)^{15−t} =
z^{15}y^{15−t}, so φ_i(y) := Σ_t ψ_{i,t}y^{15−t} = R⁻[e_i(P)·∏_M(y+C_M)] and after
skew-folding (Chebyshev: z^{j}−z^{−j} = (z−z⁻¹)S_{j−1}(y)):

> **c(z) = −z^{18}·(z−z⁻¹)·Φ(y)**, Φ := (y²−1)φ₀ + yφ₁ + φ₂ (deg ≤ 15, G1)

— the trivial s-mirror is the factor (z−z⁻¹); the entire two-row object is ONE
polynomial in the cosine variable. In the S-basis (S_m(2cosθ) = sin((m+1)θ)/sinθ):
n\_(λ₁,k) = (−1)^d·b_{μ₁−1} (a-priori table gate: 27/27 strip cells incl. zeros).

**(3) THE SEAM IDENTITY (scout 29, all gates green).** Divide by S₁₀(y) =
sin(11θ)/sinθ = ∏_{r=1}^{10}(y − 2cos(rπ/11)) — the minimal polynomial of the
Niven 11-lattice, and note (z−z⁻¹)·S₁₀(y) = z¹¹−z⁻¹¹ = the committed wall factor:

> **Φ_k(y) = S₁₀(y)·Q_k(y) + R_k(y)** with (verified by exact division)
> Q₀ = S₁−S₅ = −P₀ | Q₁ = 2S₀+S₂−S₄ = −P₁ | Q₂ = 3S₁ = P₂ | Q₃ = 2S₀+S₂ = −P₃ |
> Q₄ = S₁ = P₄ | Q₅ = 0 = P₅  (the committed column polynomials ARE the quotients)
> and **deg R_k ≤ 4+k** for every k (0: −1(zero!), 1: 4, 2: 5, 3: 4, 4: 7, 5: 8).

**Why this IS the law:** S₁₀·S_j = S_{10+j}+S_{8+j}+…+S_{10−j} has S-coefficients
symmetric about 10 on [5,15]; R_k feeds only m ≤ 4+k, below the k-window's low
edge m = 5+k. So deg R_k ≤ 4+k ⟹ b_{10+u} = b_{10−u} on the window ⟹ F134.
The Θ-shadow decomposition (round 2C, "equivalent, not a proof") is now the
polynomial division: Θ_k ↔ S₁₀Q_k, the beyond-wall shadow W_k ↔ R_k.

**The fence lands exactly (G4):** l=1 remainders deg 3 (the l=1 hold), l=2
remainders deg 8–9 (OVERFLOW = the famous breaks), k=0 remainder ZERO
(Φ₀ = S₁₀·(S₁−S₅), a clean polynomial identity).

**Epistemic grade:** the chain F133-letters → ψ → Φ → division is a-priori (the
committed table enters only as gate); every reduction step is a hand-provable
lemma (pair factorization, z↦y, skew fold, Chebyshev algebra, Pieri, bottom
recursion, slice identity) and the end is a finite exact-ℤ division — the
F127-wall certificate class. F134 is thereby DERIVED from F133: the seam
identity the arc was hunting. STILL OPEN (beauty, not validity): a structural
derivation of deg R_k ≤ 4+k (why the shadow is small) without performing the
division; and the θ₁/Appell–Lerch reading of the remainder tower.

## Files
- `_rb_theta1.py` — the spectral pin: basis, exact amplitudes, [½,½], central frequency.
- `_rb_theta2.py` — factor/rank/column-inventory (single product dead).
- `_rb_theta3.py` — full-atlas Conway–Jones resonance scan (closed negative).
- `_rb_theta4.py` — multiplier/regauge probes; the T4c diagonal find N(q,q) = −1.
- `_rb_theta5.py` — the telescope's reach: two-row exact, three-row defect odd about 11.
- `_rb_theta6.py` — Φ anti-diagonal anatomy (+ self-generated pkl cache `_rb_theta6_phi.pkl`);
  lesson: Φ point values ≠ n (the read-off is the W(C₂)-antisymmetrization ψ).
- `_rb_theta7.py` — line anatomy at Xc level: d ≥ 1 clean, d = 0 carries the corner shadow.
- `_rb_theta8.py` — the 2D map: both diagonal telescopes, the corner, rank-8 wedge kill.
- `_rb_theta9.py` — round 3 (i): per-l diagonal anatomy, the l=1 independence, l=2 multiset
  mismatch, pooled T5c reproduction.
- `_rb_theta10.py` — round 3 (ii): rung resolution, ONLY-FULL-SUM, carrier profiles.
- `_rb_theta11.py` — round 4 (b): pair-(μ₂,μ₃) balance test, universality refuted.
- `_rb_theta12.py` — round 5: the determinacy audit + the anchor gate (laws + (−1,−2,−3)
  on the mirror-fixed locus = the whole two-row sector).
- `_rb_theta13.py` — round 6 (b): the wall inventory (34 entries, |ν| odd, 34/34 vs table).
- `_rb_theta14.py` — round 6 (c): the rank-5 build (X₅ half-shifted lattice), naive descent
  refuted.
- `_rb_theta15.py` — round 7 (a,b): Λ⁷ slice verified (Xc samples + wall 34/34), grade
  profiles, full-system zero.
- `_rb_theta16.py` — round 7 (c): the deletion cascade Λ⁷ → Λ⁶ → Λ⁵, floor at the rank.
- `_rb_theta17.py` — round 8: Θ_wall closed form (support exactly 34), obstruction/remainder
  spectra (shared 68-support, off-wall cancellation).
- `_rb_theta18.py` — round 9: raw (σ,ε)-subset reduction of the anchors (verifies −2/−3;
  267 signed terms, no collapse yet).
- `_rb_theta19.py` — round 9: the π-form (δ+3𝟙 = ρ₅) + the bottom-flip involution
  (correct but weak: 6382 survivors).
- `_rb_theta20.py` — round 9, the KEPT WRONG TURN: κ = 3−j sign error makes every block
  except (0,0) vanish and lands one off (−3/−4); also the x-vs-t unit trap in Pieri.
- `_rb_theta21.py` — round 9: the clean collapse (L1 Pieri-prefix, L2 +j𝟙 recursion,
  L3 W(C₅) Racah blocks + stabilizer kill); wall column 34/34; anchor block tables.
- `_rb_theta22.py` — round 9: hand-checkable per-block solution inventories for both
  anchors (|arg|-classes, net signs, c-values).
- `_rb_theta23.py` — round 10: general slice + universal block formula (whole two-row
  table 14/14); per-(a,κ) palindromy refuted.
- `_rb_theta24.py` — round 10: r-diagonal / virtual-Ξ resummation (r = 0,1 vanish);
  per-r palindromy refuted.
- `_rb_theta25.py` — round 10: the Core/P split, mirrors (M1)/(M2), channel tables;
  per-channel C22 refuted.
- `_rb_theta26.py` — round 10 GATE: the seam form (channels reconstruct the table,
  odd-part-only reconstruction, A(ν,d) support = the break atlas).
- `_rb_theta27.py` — round 11: the cosine matrix ψ (Core pair-factorization, binomial
  transform, 14/14).
- `_rb_theta28.py` — round 11: the one-variable collapse Φ(y), S-basis, window
  palindromy, Φ mod S₁₀ remainders.
- `_rb_theta29.py` — round 11 GATE: the seam identity end-to-end (deg lemma, a-priori
  table 27/27, division Q_k = ±P_k with deg R_k ≤ 4+k, l=2 remainder overflow).
