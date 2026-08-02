# The minimal crystal (k=4): the 32 integers, opened three depths -- OPAQUE, and why

Scratch, local-only (`_rb_cry{1..5}.py`). Run from repo root, `PYTHONIOENCODING=utf-8`.
Exact over ℤ/ℚ throughout (`Xc, Yc, Fk` from `_rb_summed`; `P1, P2` from
`f133_w_closed_form.build_halves`; `CHI_TAB` = the committed F133 chiC table).
Continues `_rb_lad` (the ONLY-FULL-SUM verdict) with the microscope that round did not use.

## TL;DR

- **The crystal IS a single F133 coefficient reflection.** Telescoping collapses the comb:
  `C_ε(m) = sig_ε · Xc(m, T_ε)` exactly (all six ladder coordinates shift together), so at
  k=4 `a_ε = sig_ε·[Xc(20,T_ε) − Xc(24,T_ε)]`, and `F_4(m) = +n_(m/2−6, 4)` — the two-row
  Sp(12) character coefficient of the F133 table. The whole 32-integer sum-rule
  `Σ_ε a_ε = 0` is the innermost cell **n_(4,4) = n_(6,4) = 1**. (Verified against the
  committed `chiC_coeffs.txt`, ratio `F_k/n = +1` on all 34 live cells, `_rb_cry3` (A).)
- **Microscope: OPAQUE at every constituent depth, and provably so.** The vanishing never
  becomes a shared-constituent column-cancellation. It is an irreducible **bilinear** value
  identity among independent alternant coefficients — exactly the arc's theta/affine residual.
- **Moment lens: a faithful restatement, not a mechanism.** Odd moments
  `M_p^(k) = 2^p Σ_j (j−5)^p n_(j,k)` vanish `⟺` the coefficient row is palindromic about
  `j=5`. Same obstruction as pointwise (per-ε first-coordinate moment is even in `b` for 0/32).
- **Fence respected:** both tools report ASYMMETRIC at the ρ4-resonance `l=2` for every
  reachable `k`, SYMMETRIC at `l≠2` (incl. l=3,4). `_rb_cry5`.

## Task 1 — the 32 integers `a_ε = sig_ε·[Xc(20,T_ε) − Xc(24,T_ε)]` (`_rb_cry1`)

`Σ_ε a_ε = 0` (= `F_4(20) − F_4(24) = 1 − 1`). Table (ε as +/- signs; `sg=∏ε`; `sig=sg·ss`):

| ε | a_ε | factor | ε | a_ε | factor |
|---|---|---|---|---|---|
| +++++ | 0 | 0 | -++++ | -11 | -11 |
| ++++- | 0 | 0 | -+++- | 25 | +5² |
| +++-+ | -2 | -2 | -++-+ | 97 | +97 |
| +++-- | 8 | +2³ | -++-- | -175 | -5²·7 |
| ++-++ | 2 | +2 | -+-++ | -21 | -3·7 |
| ++-+- | -3 | -3 | -+-+- | 8 | +2³ |
| ++--+ | -37 | -37 | -+--+ | 118 | +2·59 |
| ++--- | 37 | +37 | -+--- | -47 | -47 |
| +-+++ | 7 | +7 | --+++ | -182 | -2·7·13 |
| +-++- | -32 | -2⁵ | --++- | 305 | +5·61 |
| +-+-+ | -8 | -2³ | --+-+ | 85 | +5·17 |
| +-+-- | 34 | +2·17 | --+-- | -135 | -3³·5 |
| +--++ | 104 | +2³·13 | ---++ | -178 | -2·89 |
| +--+- | -86 | -2·43 | ---+- | 63 | +3²·7 |
| +---+ | -87 | -3·29 | ----+ | 69 | +3·23 |
| +---- | 63 | +3²·7 | ----- | -21 | -3·7 |

Range `[-182, 305]`, 28 distinct values, 2 zeros (`+++++`, `++++-`).

- **± pairing at the integer (multiset) level: NO** (asked explicitly). The multiset
  `{a_ε}` is *not* invariant under `v ↦ −v`: e.g. `305, 104, 118, 97, 85, 69, 63(×2)` have no
  negatives present; only `{±2, ±8, ±37}` happen to pair. So `_rb_lad`'s "arithmetically
  generic" survives the value-level test too — there is no hidden ±symmetry it missed.
- **Natural subgroup partial sums** (all-32 = 0): `sg=+1 → −36 / sg=−1 → +36`;
  `sig=+1 → +16 / −1 → −16`; last-sign `e4=+1 → −44 / −1 → +44`. Clean ± halves — but these
  are just the balanced-splits of one integer, no finer sum-rule.
- **Head-sign `e0` splits into 0/0** (each 16-block sums to 0). Genuine but **k=4-only**: at
  k=3 the same split gives ±5 (`_rb_cry1` note + direct check). A rank-1-stub coincidence
  (window = 3 points), not the germ — matching `_rb_lad`'s dismissal, now pinned as stub-only.
- **mod small primes**: no congruence structure (mod 2: 14 even/18 odd; mod 3/5/7/11/13 all
  spread). Generic.

**Most surprising arithmetic fact.** These 30 nonzero integers, ranging to ±305 with genuinely
generic factorizations, are the flip-sheet decomposition of the **most trivial identity in the
whole arc: `1 = 1`** (`F_4(20)=F_4(24)=1`). The "minimal crystal" bottoms out at two equal 1's
in the F133 Sp(12) character table; the fireworks above are one coefficient, computed 32 ways.

## Task 2 — the microscope: OPAQUE at every constituent depth (`_rb_cry2`, `_rb_cry3` (B), `_rb_cry4`)

Three depths opened; the vanishing is visible at none of them, for a structural reason.

- **Depth 1 — Xc targets under `⟨S₆-alternation, total negation⟩`** (`_rb_cry2`). Xc obeys
  `Xc(Pν)=sgn(P)Xc(ν)` and `Xc(−ν)=Xc(ν)` (proved, re-checked here). Canonicalizing the 64
  raw constituents `(head, T_ε)`, `head∈{20,24}`, gives **64 distinct orbits — no two merge**.
  Every column is a singleton ⇒ zero structural cancellation; the residual
  `Σ (coeffsum)·Xc = 0` is a sum of 60 nonzero, symmetry-unrelated integers.
- **Why no deeper depth can help** (`_rb_cry3` (B)). At the atomic monomial level
  `Xc(ν)=[X]_ν = Σ_{p1+p2=ν} P1[p1]P2[p2]`, a matching `(p1,p2)` has a **fixed sum = its
  target**. The 64 targets are pairwise distinct, so no matching is shared between two
  constituents: the `ε × matching` incidence matrix is **block-diagonal**. A block-diagonal
  system has no cross-column cancellation at *any* refinement of the monomial index. This is a
  proof that the constituent microscope cannot expose the vanishing — and it *explains*
  `_rb_lad`'s ONLY-FULL-SUM: there is provably no shared sub-object to group.
- **Depth 2 — the shared half-columns** (`_rb_cry4`). The one index that *is* shared is the
  meet-in-the-middle half (`X=P1·P2`, `|P1|=590016`, `|P2|=5817`). Grouping `Σ_ε a_ε` by
  P1-monomial `q` (coeff `sig_ε(P2[(20,T)−q]−P2[(24,T)−q])`) and by P2-monomial `p`
  (symmetrically): only **8.8 %** of P1-columns and **5.8 %** of P2-columns cancel
  structurally. The overwhelming residual needs the actual *values* of the other half — the
  identity is genuinely **bilinear** (`Σ coeffsum·P1[q] = 0` holds only after multiplying by
  P1). Half-depth = OPAQUE.

**Verdict:** STRUCTURED-AT-DEPTH = *none*. The crystal's vanishing is an irreducible bilinear
identity among independent Sp(12) alternant coefficients. This is a sharp negative asset: the
germ is **not** any partition/grouping/constituent-cancellation (all now excluded with proof),
confirming the arc's placement of the residual at the theta/affine-character node.

## Task 3 — the moment lens: faithful restatement (`_rb_cry3` (C,D))

Because `F_k(m) = n_(m/2−6, k)` and `m−22 = 2(j−5)` (`j=m/2−6`),
`M_p^(k) = Σ_m (m−22)^p F_k(m) = 2^p Σ_j (j−5)^p n_(j,k)`. All odd moments vanish (verified,
all k), **⟺ the two-row coefficient row `n_(·,k)` is palindromic about `j=5`** — which *is*
the theorem. The lens does not go one tier up:

- the per-ε first-coordinate moment `Σ_b b^p Xc(22+b, T_ε)` is **not** even in `b` (that is
  the per-ε defect, 0/32 in `_rb_lad`); only the ε-sum vanishes — the *same* obstruction as
  the pointwise law, so the odd-moment vanishing does **not** factor through an E-cosine parity.
- the even moments are geometric (`k=1: 2·8^p`, `k=2: 6·2^p`, ...) **only because each row's
  support is a single symmetric pair `{j,10−j}` with equal entries (+ optional center `j=5`)**:
  `k=0:(0,2,8,10)`, `k=1:(1,·,9)+c5`, `k=2:(4,6)`, `k=3:(3,·,7)+c5`, `k=4:(4,6)`. The
  palindrome is manifest in support+values — but *that manifestation is the theorem*, not a
  proof of it.

**Verdict:** the moment lens is a change of basis that faithfully mirrors the reflection; it
exposes the sparse symmetric-pair support (a genuine descriptive fact) but yields no mechanism.

## Task 4 — the fence (`_rb_cry5`), three-row tail `(2k+10, 2l+8, 6, 4, 2)`

Both descriptive tools re-run at every reachable `(k,l)`; the ρ4-resonance `l=2` must break.

| k | l=0 | l=1 | l=2 | l=3 | l=4 |
|---|---|---|---|---|---|
| 2 | SYM | SYM | **ASYM** (M₁=−8) | — | — |
| 3 | SYM | SYM | **ASYM** (M₁=12) | SYM | — |
| 4 | SYM | SYM | **ASYM** (M₁=−6) | SYM | SYM |

Odd moments are exactly 0 at `l≠2` and nonzero at `l=2`; the coefficient row is palindromic
`l≠2`, non-palindromic `l=2` (e.g. `k=2,l=2` row `(1,0,4,0,3,0,0)`, `F(20)=4≠3=F(24)`). The
tools **fail to certify l=2** and correctly certify l=3,4. Reproduces `_rb_lad6` / `_rb_kt6 C4`
through the crystal lens. Fence respected — no tool here proves the false statement.

## Files
- `_rb_cry1.py` — the 32 integers + factorizations + ± multiset test + subgroup sums + mod primes.
- `_rb_cry2.py` — Xc-target microscope under ⟨S₆, total-neg⟩ (64 distinct orbits).
- `_rb_cry3.py` — crystal = `n_(j,k)` identification (A); block-diagonal obstruction (B);
  two-row table + reflection (C); odd/even moments (D).
- `_rb_cry4.py` — the shared half-column microscope (P1 8.8 %, P2 5.8 % structural; bilinear residual).
- `_rb_cry5.py` — the fence: moment lens + row-palindrome at three-row `(k,l)`, break at l=2.

## Net
The minimal crystal is `n_(4,4)=n_(6,4)=1`, a single Sp(12) character-coefficient reflection.
Opened three depths: the vanishing is **provably not** a constituent-cancellation (block-diagonal
incidence at the Xc/monomial level; only ~6–9 % column-cancellation at the shared half-depth) —
it is an irreducible **bilinear** identity. The moment lens faithfully restates it as row-palindromy
without lifting to a parity mechanism. Everything breaks correctly at l=2. Residual unchanged in
kind: the affine C₆ / twisted-theta node the level/kt/bnd/lad rounds already named.
