# Endgame round: the shifted deleted-sheet series (germ for the proof doc)

Scratch, local-only. Gates: `_rb_endgame2.py` (Xc/Yc, content_full/tail, GATE 0/1),
`_rb_endgame3.py` (telescope + survivor c, GATE 2/3), `_rb_endgame4/5/6.py` (pairing tests, GATE 4-7 + profiles).
Run from repo root, `PYTHONIOENCODING=utf-8`. All arithmetic exact over ℤ (dict Xc).

## Definitions
- `Xc(ν) = [X]_ν`, X = Δ̂·∏₃₁ (31-sheet product, 1⁶ deleted). `Yc(ν) := Xc(ν−𝟙) − Xc(ν+𝟙) = 2i[Y]`,
  Y = Δ̂·∏₃₂ (full 32-sheet product), 𝟙 = (1,1,1,1,1,1).
- `content_tail(cf,tail,m)` = Σ_{ε∈{±1}⁵} sgn(ε)·sds-sign·cf(sort[m, ε∘tail]): the C₅ read-off,
  slot-1 = slice var m (unflipped), slots 2-6 = the 5-tail flipped.
- `F_k(m) = content_tail(Xc,(2k+10,8,6,4,2),m)`; the two-row law is **F_k(m) = F_k(44−m)** in-window.
- `G_k(m,c) = Σ_ε sgn(ε)·sds-sign·Yc(sort[m,ε∘tail] + c·𝟙)`: the full-sheet content shifted uniformly by c.

## Proven exact over ℤ (this round)
- **GATE 0/1 — the full sheet contributes nothing.** Φ_Y (full C₆ content of Y) = 0 on all two-row
  points (0 nonzero); the C₆ content of X is live (14 nonzero). STRONGER: the 1-variable **slice**
  content `FY_k(m) ≡ 0` for all m, all k. So the reflection is a *pure deleted-sheet effect* — nothing
  survives from the full-32-sheet product at any slice.
- **GATE 2 — telescoping inverse (exact identity).** `Xc(ν) = Σ_{r≥0} Yc(ν+(2r+1)𝟙)` (0/6 probes).
  Proof: Σ_r [Xc(ν+2r𝟙) − Xc(ν+(2r+2)𝟙)] telescopes; Xc has finite support. No 2i ambiguity (both integer dicts).
- **GATE 3 — finite series + survivor set.** `F_k(m) = Σ_{c odd} G_k(m,c)` (0/81 mismatches, all k),
  with `G_k(m,0)=FY_k(m)=0` and **survivor support c ∈ {5,7,9,11,13,15,17,19}** for EVERY k
  (c=1,3 vanish identically; series truncates at 19). The set is symmetric about c=12.
- **GATE 5/7 — the reflection + defect law.** `F_k(m)=F_k(44−m)` in-window (0 fails, all k). The below-window
  defect points are exactly the R_k support, and each has its partner 44−m strictly outside the wall 32−2k
  (verified: e.g. k=0 defect m∈{−32,−28,−16,−12} → partners {76,72,60,56} all past wall 32). So
  **out-of-wall partner ⇔ unpaired ⇔ R_k defect**, exactly the SERIES agent's "leaves the wall."

## Refuted exact over ℤ (this round) — the clean pairing is DEAD
- **GATE 4/6 — no term-by-term involution.** The natural hypothesis `G_k(m,c) = G_k(44−m, 24−c)`
  (24 = 2·center) FAILS, even restricted to in-wall partners (34 mismatches k=0, down to 6 at k=5),
  and every signed/shifted variant fails too.
- **GATE 6 profiles — WHY (decisive).** At a reflected pair the c-profiles are structurally DISJOINT yet
  equal-sum. k=1: m=14 → [9,−34,36,−10] (spread c=5..11) vs m=30 → [1,0,0,0] (only c=5); both sum to F=1.
  k=2: m=20 → [−6,21,−12] vs m=24 → [−3,7,−1]; both sum to 3. Near the wall the profile collapses to a
  single c; near the window's lower edge it spreads. **The reflection reorganizes the entire series; it is
  intrinsically a summed identity, not a per-term one.** This confirms, now constructively, the parent's
  "structural obstruction": no character-multiplication / shift operator realizes the first-coordinate
  reflection while fixing the tail (M=∏t_u shifts all six coords — GATE 6 is that fact made explicit at the
  coefficient level).

## Where 44 and the wall are pinned (derived, not fitted)
- 44 = 2·22, the μ₁-reflection node (μ̂ = 2μ, wall λ₁+λ₂≤10 ⇒ μ₁+μ₂≤21, affine node 11 → 22 → 44).
- The c-survivor SET {5,…,19} is **k-independent** (GATE 3); the k-dependent wall 32−2k enters not through
  the c-cutoff but through which (m,c) cells land in Yc's support — the defect points (GATE 7) are exactly
  those whose reflected partner 44−m falls past 32−2k.
- 24 = 2·12 = the c-survivor center (symmetry of the survivor SET, though not of individual terms).

## Net proof status (honest)
- **Φ_Y ≡ 0 and its slice-strengthening FY_k ≡ 0: PROVED** (parent ALGEBRA sub-lemmas (a)/(b)/(c) +
  this round's GATE 0/1). The reflection is therefore a **pure deleted-sheet phenomenon**, expressed as a
  finite shifted series with an explicit survivor set — all exact over ℤ.
- **F_k(m)=F_k(44−m): established as a FACT to very high assurance** (independent exact-ℤ verifications:
  the slice content directly, the finite series, the in-wall/defect decomposition) but **NOT reduced to a
  closed term-level involution.** This round *refuted* the last natural candidate mechanism (the c-pairing)
  with an explicit constructive counterexample, so the residual is now sharply named: the two-row reflection
  is a genuinely summed (non-involutive) identity of the deleted-sheet series. Any future closure must be a
  summed/generating-function argument (e.g. a rationality/period-22 statement on Σ_c G_k(m,c)·q^c), NOT a
  bijection on (m,c). The "one covariance lemma" of the resumption point is thereby *transformed*, not
  discharged: Φ_Y-covariance is proved; the surviving open piece is the summed reflection, and the pairing
  route to it is now a documented dead end.
