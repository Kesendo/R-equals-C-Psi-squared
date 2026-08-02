# Summed-identity round: the affine reflection splits (finite negation PROVED × affine translation)

Scratch, local-only. Gates: `_rb_summed.py` (g_{k,m}(q) generating functions),
`_rb_summed2.py` (negation / period / defect decomposition), `_rb_summed3.py` (parity + anti-period
window), `_rb_summed4/5.py` (witness harvest + term-level symmetries), `_rb_summed6.py` (consolidated
gate). Run from repo root, `PYTHONIOENCODING=utf-8`. All arithmetic exact over ℤ (dict `Xc`). Reuses
`_rb_endgame*.py` definitions (`Xc, Yc, Fk, Gc`).

## The result of this round
The two-row reflection `F_k(m) = F_k(44−m)` (⇔ n_(j,k)=n_(10−j,k)) is an **affine Weyl reflection**
`μ₁ ↦ 22−μ₁` (m = 2μ₁, so m ↦ 44−m). Every affine reflection factors as **finite reflection ∘
translation**. This round PROVES the finite half and isolates the residual as the pure affine
coroot translation.

### PROVED this round (finite-Weyl half), all exact over ℤ, 0 exceptions
- **L1 — Xc is EVEN under total negation.** `Xc(−ν)=Xc(ν)` (0/8586 nonzero witnesses).
  Proof: Yc is ODD (Y=Δ̂·∏₃₂, 32 sheets ⇒ total sign flip gives (−1)^{15+32}·... = −1; confirmed
  Yc(7,5,3,1,−1,−3)=−2824 ↦ +2824), and the two-directional telescope
  `Xc(ν)=Σ_{r≥0}Yc(ν+(2r+1)𝟙) = −Σ_{r≥0}Yc(ν−(2r+1)𝟙)` (GATE 2) turns Yc-odd into Xc-even in one line.
- **L2 — term-level c/m negation.** `G_k(m,c)=G_k(−m,−c)` (0/954).
- **L3 — two-sided telescope vanishing.** `Σ_{c odd ∈ ℤ} G_k(m,c)=0` (0/354), so
  `Σ_{c<0}G_k(m,c) = −F_k(m)` (forward telescope = −backward telescope).
- **⇒ F_k(m) = −F_k(−m)** (the content is ODD in m): 0/726, all k. Clean proof:
  `F_k(−m) = Σ_{c<0}G_k(m,c)` [L3] `= Σ_{c>0}G_k(−m,c)` [L2] `= F_k(−m)`… i.e.
  `−F_k(m) = Σ_{c<0}G_k(m,c) = Σ_{c>0}G_k(−m,c) = F_k(−m)`. Finite-Weyl half discharged.

### The residual (the affine half), precisely named
Via the proven negation, the reflection is EQUIVALENT to the **pure affine coroot translation**
`F_k(m+44) = −F_k(m)` (anti-periodicity of period 44 = 2·22). Verified: reflection in-window ⇔
anti-period on the negated window `[−(32−2k), −(12+2k)]`, both True all k. This is the "affine C₆
wall reflection" the resumption point conjectured — now with the finite reflection stripped off, so
the open core is a single translation, not a folding.

### Two routes to the affine half, DOCUMENTED DEAD this round
- **Term-level (shifted-c bijection):** `G_k(m+44,c)` is empty (m+44 past the coord-1 support wall of
  Yc), so there is no per-term image; the translation is intrinsically summed. (Matches the endgame
  round's GATE 4/6 refutation from the other side.)
- **q-cyclotomic period:** `d(q)=g_{k,m}(q)−g_{k,44−m}(q)` has `d(1)=0` (the theorem) but
  `d(q)/(q²−1)` carries **no** cyclotomic/period-22 structure. Example k=0,m=16:
  `d = −5q⁵+17q⁷−16q⁹+4q¹¹ = −q⁵(q²−1)(2q²−1)(2q²−5)`, roots q²=1/2,5/2 — not roots of unity.
  So lever (1) (roots-of-unity specialisation) is negative.

### Net status
The affine reflection = (finite negation: **PROVED**, L1·L2·L3) ∘ (affine coroot translation-44:
**open**, but now the sole residual, with the term-level and q-cyclotomic sub-routes closed as
negative). The wall/defect law is unchanged: the translation holds while both partners lie inside
Yc's coord-1 support wall 32−2k, and fails exactly at the R_k defect (one partner past the wall).

The sharpest surviving thread: prove `F_k(m+44)=−F_k(m)` in-wall. It is an internal-window summed
identity of the deleted-sheet series; the two established handles are (i) the exact wall boundary law
and (ii) the geometric-series resolvent `X = Y·2i·Σ_r M^{−(2r+1)}` whose ratio M^{−2} IS the affine
coroot translation by 2 in every coordinate — 22 steps of it give the 44. Closing it needs the
in-wall truncation of that 22-step translation to reflect rather than vanish; neither a bijection nor
a q-period, so a genuinely summed / boundary argument on the resolvent tail.
