# Round KW: the super / dual-pair face of the spinor discriminant X

Scratch, local-only (`_rb_kw*.py`). Run from repo root, `PYTHONIOENCODING=utf-8`. All coefficient
arithmetic exact over ℚ. Gate: `python simulations/_rb_kw1.py` (product-form test) plus the two inline
table probes in this note. Context: `_rb_mac_README.md` (the A₁₂⁽²⁾/BC₆ refutation), `_rb_level_README.md`,
`_rb_bnd_README.md`, PROOF_F133 §5.

The object: `X = Δ̂·∏₃₁`, Δ̂ = the type-A₅ half-angle Vandermonde `∏_{u<v} sin((x_u−x_v)/2)`, ∏₃₁ = the
product over 31 of the 32 canonical spin sheets `sin(L·x/2)` (`L ∈ {±1}⁶`, `L₁=+1`, TOP weight `1⁶`
deleted). `S = X / A_ρ^{C₆} = Σ n_λ χ^{Sp(12)}_λ`, 143 integer terms. The 64 sheets `(±½)⁶` are the weight
system of the spin representation of `so(13)=B₆`; the character basis is its Langlands dual `Sp(12)=C₆`.
Open law: `n_(λ₁,λ₂) = n_(10−λ₁,λ₂)` (two rows), `μ₁ ↦ 22−μ₁`, character `(−,+)`, breaks at three rows.

## Verdict up front (per sub-face)

- **(a) affine Lie superalgebra super-denominator (osp/spo, strange series, exceptional): NO-MATCH.**
  Every finite/affine Lie superalgebra whose even part carries `so(13)⊕sp(12)` has PAIRWISE odd roots
  `±ε_i±δ_j` (rank ≤ 2), never the six-fold spinor sheets. The one genre where B-type spinor sheets ARE
  fermionic roots is the exceptional F(4), and only for `B₃` (`so(7)`), with no `B₆` analogue.
- **(a′) specialized Weyl numerator of a larger algebra restricted to the B₆ torus: NO-MATCH.**
  The spin weights are not a root system (not closed under reflection); the 46-factor (or 32-sheet) count
  matches no rank-6/7 positive-root system; no |1|-graded simple algebra has `B₆⊕ℂ` with the 64-dim spin
  rep as `g₁`. X is genuinely a *minuscule-representation discriminant*, not a Weyl denominator.
- **(b) King dual pairs (arXiv:2303.00576), §7-8: PARTIAL — GENRE MATCH, no closed identity.**
  King's §8 is exactly the spin↔symplectic dual-pair machinery the F133 proof §5 named, and his coefficient
  laws (8.9) are reflection-with-sign of the SAME genre as our `n_(j,k)=n_(10−j,k)` with `(−,+)` character.
  But the direct dual-Cauchy object is a PRODUCT `∏_i Q(c_i)`; our K is provably NOT product-form
  (`_rb_kw1.py`), and our support is a STAIRCASE `λ₁+λ₂≤10`, not King's rectangle nor Rains-Warnaar's box.
- **(c) numeric fence: PASS.** The two-row reflection holds on ≤ 2 rows (`0/10`), breaks at three rows
  (`6/16`, the required l=2 break), wall `λ₁+λ₂≤10` exact (`0` violations).

Net: the super/dual-pair face refines round 1's naming. X is a **deleted-top-weight spinor discriminant on
the (B₆,C₆) Langlands seam**; the dual-pair genre (King §8, Hasegawa [Has]) supplies the reflection-with-sign
coefficient symmetry but no single closed identity, because our object carries the extra A₅ Vandermonde,
the deleted sheet, and a B-type staircase support glued to a C-type character basis.

## (a) Affine Lie superalgebra super-denominators — refuted with reasons

**The candidate.** Kac–Wakimoto super-denominators famously produce windowed theta identities with
defect/atypicality corrections and mixed signs (searched and confirmed: Kac–Wakimoto 1994; recent
orthosymplectic instance Matsusaka–Suzuki, arXiv:2502.06449, "Denominator identity for the affine Lie
superalgebra spo(2m,2m+1)"). The even part of `spo(2m,2m+1)=osp(2m+1|2m)` is exactly `so(2m+1)⊕sp(2m) =
B_m⊕C_m`; for `m=6` that is `so(13)⊕sp(12)` — precisely our spin/symplectic pairing. Natural hypothesis:
X is a finite specialization of this super-denominator.

**Refuted (structural, decisive).** From the source (Matsusaka–Suzuki §2, fetched), the finite root system
of `spo(2m,2m+1)` is:
- even roots `Φ₀`: `±(ε_i−ε_j), ±(ε_i+ε_j), ±2ε_p` (the `C_m` = sp part) and `±(δ_i−δ_j), ±(δ_i+δ_j), ±δ_e`
  (the `B_m` = so part);
- **odd roots `Φ₁`: `±(ε_i−δ_j), ±(ε_r+δ_s), ±ε_p`** — all PAIRWISE `ε±δ` (rank ≤ 2), plus single `±ε_p`.

The super-denominator `R = ∏_{Φ₀⁺}(1−e^{−α}) / ∏_{Φ₁⁺}(1+e^{−α})` therefore has a fermionic factor over
`ε±δ` pairs, NOT over the six-fold sheets `L·x/2 = (±½,…,±½)·x`. This is the **identical obstruction as the
A₁₂⁽²⁾ refutation R2** (`_rb_mac_README.md`): no spinor factor exists in the denominator, so the deleted
`1⁶` sheet has no home there. A super-denominator NUMERATOR is a single alternant (Weyl–Kac character), not
a 46-fold product. Verdict: NO-MATCH.

**The strange / exceptional series — refuted, and the one near-miss named.** The ONLY Lie superalgebra whose
FERMIONIC roots are literally B-type spinor sheets is the exceptional **F(4)** (confirmed, arXiv:1309.0418 and
the ABC monograph): even part `B₃⊕A₁ = so(7)⊕sl(2)`, odd roots `{½(±ε₁±ε₂±ε₃±δ)}` = the `so(7)=B₃` spinor
sheets `(±½,±½,±½)` tensored with `±½δ`. So spinor sheets as fermionic roots occur **only at n=3**. The finite
exceptional list is `F(4), G(3), D(2,1;α)` — there is NO `B₆` (n=6) analogue, hence no super home for our
32-sheet product. This is a clean, complete refutation, not a "not found."

## (a′) Larger-algebra Weyl-denominator restriction — refuted with reasons

The promising-and-cheap angle from the brief: is `∏_{L>0} 2sin(L·x/2)` (the spinor product SP) a
specialization of the Weyl denominator of a LARGER algebra (D₇, B₇, E₇) restricted to the B₆ torus?

- **The spin weights are not roots.** `(±½)⁶` is not closed under its own reflections: reflecting
  `v=(½,½,½,½,½,−½)` in `w=(½,…,½)` gives `v−(4/3)w`, not a `±½` vector. So the 64 sheets are not a
  sub-root-system of anything; they are minuscule WEIGHTS (the spin rep of `B₆` is minuscule — `W(B₆)` acts
  transitively on `(±½)⁶`), and weights of a minuscule rep are never the roots of a larger simple algebra
  here.
- **Factor counts match no root system.** X = 15 (A₅ differences) + 31 (sheets) = **46** factors; the full
  32-sheet version is 47. No rank-6 system has 46/47 positive roots (A₆:21, D₆:30, B₆/C₆/E₆:36), and no
  rank-7 system either (A₇:28, D₇:42, B₇/C₇:49, E₇:63). Same as R1 of `_rb_mac_README.md`, extended to
  rank 7. The D₇→B₆ folding gives B₆ ROOTS not spin weights (the brief's own note, confirmed).
- **No graded home.** A |1|-graded simple algebra `g = g₋₁⊕g₀⊕g₁` with `g₀ ⊇ B₆` and `g₁` = the 64-dim
  spin rep would need `dim g = 78+1+64+64 = 207`; no simple Lie algebra has this dimension (E₇=133, B₇=105,
  D₈=120, E₈=248). Equivalently, in a Dynkin diagram the spin node of `B₆` (`o-o-o-o-o⇒o`) sits behind the
  double bond, and no finite-type diagram attaches a further node there. So the spin rep is not the `g₁` of
  any minuscule/contact grading. NO-MATCH.
- **What SP actually is (derived, `_rb_kw` inline).** The spinor sine-product trivialises only at rank 2:
  `SP₂ = −½·(cos x₁ − cos x₂)` (a single cosine Vandermonde factor, exact). For `n≥3` SP is NOT any power
  `const·Vand_c^p` (`p≤3` tested, non-constant). So SP is a genuine spinor discriminant with no elementary
  root-system closed form — consistent with PROOF_F133 §5's "SP in χ^{C₆} is not compact (1096 terms)."
  The literature term "spinor discriminant" is not standard (searched; no hits) — the name is ours,
  descriptive, not a citation.

## (b) King dual pairs, arXiv:2303.00576 §7-8 — the genre match, exact

King, *Generating functions for some series of characters of classical Lie groups* (math.CO, 2023). Abstract
punchline, verbatim: "An alternative approach is then based on dual pairs of symplectic and/or orthogonal
groups. A byproduct of this approach is that expansions in terms of spin orthogonal group characters can
always be recovered from non-spin cases." This IS our (B₆-spin ↔ C₆-symplectic) recovery.

**§8, the dual-pair identities (verbatim eq. numbers).** The basic spin character of an even orthogonal
group, restricted to a Howe dual pair, expands over rectangular partitions `λ ⊆ (m^n)`:
- (8.2) `∏_{i=1}^n ∏_{j=1}^m (x_i+x̄_i+y_j+ȳ_j) = Σ_{λ⊆(m^n)} ch^{Sp(2n)}_λ(x,x)·ch^{Sp(2m)}_{λ̃}(y,y)`;
- (8.6) the same LHS times the spin factor `∏_i(x_i^{1/2}+x_i^{-1/2})` `= Σ_λ ch^{SO(2n+1)}_{Δ+λ}(x,x,1)·
  ch^{Sp(2m)}_{λ̃}(y,y)` — the one that NATIVELY carries `SO(2n+1)` (n=6 → `SO(13)`) SPIN characters with
  `Sp(2m)` dual coefficients.

**The reflection-with-sign coefficient laws (8.9), verbatim.** These are the "spin recovered from non-spin"
coincidences:
- `φ^{SO(2n+1)}_r(a½,1) = φ^{Sp(2n)}_r(a½−1,1)` (a SHIFT — the deletion genre);
- `φ^{O(2n)}_r(a½,1) = (−1)^r φ^{SO(2n+1)}_r(−a½+1,1)` (a REFLECTION `a½↦−a½+1` with sign `(−1)^r`);
- (two-row) `ψ^{O(2n)}_{q,r}(a½,a3/2,1) = (−1)^q ψ^{SO(2n+1)}_{q,r}(a½−a3/2+1,−a3/2+1,1)`.

The two-row `ψ_{q,r}` reflection-with-sign is structurally OUR genre: a reflection of the highest-weight
parameter carrying a parity sign `(−1)^q` — the same shape as `n_(j,k)=n_(10−j,k)` with character `(−,+)`.
King's coefficients `n_λ = ch^{Sp(2m)}_{λ̃}(y₀)` at fixed dual `y₀` are integers, bounded, and reflect under
rectangle complementation `λ↦λ̃` — precisely the "residue-periodic small-integer coefficients" the F133 proof
§5 flagged. The deepest structural kin is King's ref **[Has] = Hasegawa, PRIMS 25 (1989), "Spin Module
Versions of Weyl's Reciprocity Theorem for Classical Kac–Moody Lie Algebras"** (spin-module dual pairs,
Jimbo–Miwa branching duality, affine-character modular duality) — the mechanism under King §8.

**But not a literal King identity — two decisive obstructions:**

1. **K is NOT product-form (`_rb_kw1.py`, exact over ℚ).** A dual-Cauchy object at fixed `y` is a product
   `∏_i Q(c_i)` (with `X_i = x_i+x̄_i = 2c_i`, `Q(X)=∏_j(X+Y_j)`). Two exact failures:
   - necessary: K starts at monomial degree 6, so product form needs `Q(c)=c·R(c)` ⇒ `K` divisible by
     `e₆=∏c_i` ⇒ `K=0` when any `c_i=0`. FALSE: `K|_{c₁=0}` is generically nonzero (6/6 trials).
   - decisive: separability `K(u;A)K(v;B)=K(v;A)K(u;B)` FALSE; the ratio `K(c₁;A)/K(c₁;B)` is NOT constant
     in `c₁`. So K is not `∏_i Q(c_i)` for any single-variable `Q`.
   This is the fingerprint of the DELETED top sheet — and the deletion is even more radical than
   "product minus a correction": the FULL 32-sheet object `Δ̂·∏₃₂` has **identically zero C₆ content**
   (`_rb_kw` inline: `n_raw = 0` for all λ tested, reconfirming `Φ_Y≡0` of `_rb_endgame_README.md` via the
   C₆ read-off). So ALL of X's symplectic content is the *residue* of dividing out the top sheet `sin s`;
   there is no nonzero product-form parent to correct. The dual-Cauchy product mechanism is therefore dead
   for both the 31-sheet object (K not separable) and the 32-sheet object (content ≡ 0).
2. **Support is a STAIRCASE, not a rectangle/box.** Committed table structure (`chiC_coeffs.txt`, 143 terms):
   rows `{0:1,1:3,2:10,3:16,4:37,5:22,6:54}`, `max λ₁=10`, `max λ₂=5`, `max|λ|=18`, wall `λ₁+λ₂≤10` exact
   (0 violations). King's dual-Cauchy support is RECTANGULAR `λ⊆(m^n)`; Rains–Warnaar bounded-Littlewood
   (arXiv:1506.02755, the other F133 §5 kin) uses a BOX. Our `λ₁+λ₂≤10` is neither: it is the alcove of a
   LONG highest root `θ=ε₁+ε₂` (which gives `⟨λ,θ^∨⟩=λ₁+λ₂≤ℓ`), i.e. a `B`/`D`-type affine level wall — on
   the SPIN (B₆) side — while the character basis is `C₆`. The staircase-vs-box gap the memory already
   recorded (`f127_closed_form_find`, "R-W box ≠ our staircase") is here explained: the support wall is
   B-type, the characters are C-type; the object lives on the Langlands seam, which is why no single
   classical (rectangle or box) dual-pair identity closes it.

## (c) Numeric fence (exact, on the committed n_λ; the same object as live F_k via `n_raw/2`, 0/55 in `_rb_bnd`)

- reflection `n_(λ₁,rest)=n_(10−λ₁,rest)`: 1-row `0/3`, **2-row `0/10`** (holds), **3-row `6/16`** (BREAKS —
  the mandatory l=2 break), 4-row `24/37`. Any candidate that did NOT break at 3 rows would be wrong; ours
  breaks, King's `ψ_{q,r}` is inherently the `m=2` (two-row) sector and does not nest to `m=3`, consistent.
- wall `λ₁+λ₂≤10`: 0 violations.
- product-form (the one concrete King template): REFUTED as above.
- No surviving CLOSED identity remained to pin further on the live `F_k` pairs; the genre match (8.9) is a
  relation-between-groups, not a single self-contained identity, so there is nothing beyond the reflection
  law itself (already `0` mismatches, `_rb_bnd`/`_rb_endgame`) to check on `F_k`.

## The single most promising surviving thread

The **coefficient-level reflection-with-sign genre of King §8 (8.9) / Hasegawa spin-module dual pairs**, read
NOT as a product/dual-Cauchy identity (that mechanism is dead both ways above) but as a special-value
symmetry of the deleted-sheet RESIDUE. The object is the residue `X = (Δ̂·∏₃₂)/sin s` of a parent that has
identically zero C₆ content; its `n_λ` reflect as `n_(λ₁,λ₂)=n_(10−λ₁,λ₂)` with parity character `(−,+)`,
the exact shape of King's `ψ^{O(2n)}_{q,r} = (−1)^q ψ^{SO(2n+1)}_{q,r}(reflected)` and of the spin↔non-spin
SHIFT `a½↦a½−1` (8.9). The obstruction to a CLOSED identity is the (B₆,C₆) Langlands seam: a B-type
staircase support (`λ₁+λ₂≤10`, long root `θ=ε₁+ε₂`) glued to a C-type character basis, which no rectangle
(King dual-Cauchy) nor box (Rains–Warnaar) identity spans — precisely the theta/affine-C₆ residual the
`_rb_level`/`_rb_bnd` rounds isolated. The one concrete forward move consistent with all the data: pursue
the reflection as a modular/theta functional equation of the residue generating function on the seam (the
Matsusaka–Suzuki spo(2m,2m+1) affine super-denominator, arXiv:2502.06449, is the nearest theta machinery,
even though its finite denominator is not X). Everything product-form / root-system in the super/dual-pair
face is refuted with a reason; only this theta/coefficient-symmetry thread survives, and it has no cheap
finite closure — it is a special-function statement, matching every prior round's conclusion.

## Files
- `_rb_kw1.py` — K product-form test (e₆-divisibility + separability), exact over ℚ. Both FALSE.
- inline probes in this note — table structure (rows/λ₁/λ₂/wall), reflection by row count (l=2 break),
  spinor-product SP vs Vand_c^p (trivial only at n=2), full 32-sheet C₆ content ≡ 0 (n_raw=0 all λ tested).
- Literature: King arXiv:2303.00576 §7-8 (eqs 8.1-8.9); Matsusaka–Suzuki arXiv:2502.06449 (spo root data);
  F(4) roots arXiv:1309.0418; Hasegawa, PRIMS 25 (1989) [Has]; Rains–Warnaar arXiv:1506.02755 (box kin).
