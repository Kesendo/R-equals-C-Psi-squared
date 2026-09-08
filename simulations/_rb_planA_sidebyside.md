# Plan A: the side-by-side table (the nine canvases × the Berg, read in each canvas's units)

Scratch, local-only (2026-07-21). The PTF-as-method attack of the F134 plan (memory
`f127_closed_form_find`, "THE PLAN" block): rounds 1-3 hunted the mechanism INSIDE single
canvases and proved it is not there; this round paints the Berg quantities side by side in
every canvas's units and asks the three switching questions. Scouts:
`_rb_planA_tiers.py` (tier map of the l=2 slice, runs off `_rb_dem_fvals.json`, seconds) and
`_rb_planA_crystal_pair.py` (the crystal hold vs its T3 shadow, 32 sheets side by side,
builds the X dict ~2 min). Both green 2026-07-21; the committed F134 gate and the
`_pair_mirror_scout` re-verified green the same session.

## 0. The unit maps (all affine; the "α's" of the PTF reading)

One coordinate on the reflection axis, eight dresses. μ = λ+ρ, ρ = (6,5,4,3,2,1).

| canvas unit | from μ₁ | center | reflection s₀ | translation step | proven finite mirror s |
|---|---|---|---|---|---|
| μ₁ (character/schur, dem) | μ₁ | **11** | μ₁ ↦ 22−μ₁ | 22 | μ₁ ↦ −μ₁ |
| j = λ₁ (n-table) | μ₁−6 | **5** | j ↦ 10−j | 22 | j ↦ −j−12 |
| m (m-slice, ladder, Θ) | 2μ₁ | **22** | m ↦ 44−m | 44 | m ↦ −m |
| b (Θ recentered) | 2μ₁−22 | **0** | b ↦ −b | 44 | b ↦ −b−44 |
| n (script unit) | 18−μ₁ | **7** | n ↦ 14−n | −22 | n ↦ 36−n |
| c (c-ladder survivors) | — | **12** | (induced) | — | c ↦ −c (with m ↦ −m) |
| e₁ (Φ) | = m | **22** | e₁ ↦ 44−e₁ | 44 | e₁ ↦ −e₁ |
| σ-axis (spectral, home) | σ ↔ 11 | **σ** | anti-turn | price 2σ | γ ↦ −γ |

**First side-by-side seeing (feeds Attack B).** "Center = half the translation step" is TRUE
exactly in the ρ-shifted coordinates (μ₁: 11 = 22/2; m: 22 = 44/2; b: 0; e₁) and FALSE in
every unshifted label (j: 5 ≠ 11; n: 7 ≠ 11; c: 12). The coordinates in which the two-mirror
geometry (mirror at 0 + mirror at half-step) is visible are precisely the ones carrying
ρ = the half-sum. Tom's −½ | 0 | +½ is not one more observation about the object; it is the
CONDITION on the coordinate system for the dihedral pair to be visible at all. Theta
characteristics (Attack B) are exactly the formalism whose variables are ρ-shifted.

## 1. The table: Berg quantity × canvas

Canvases: (1) m-slice F_k(m) · (2) frozen ladder (ε,r) · (3) Φ contraction (e₁,e₂) ·
(4) character table n\_(j,k) · (5) Θ decomposition (b-units) · (6) twisted-Macdonald genre ·
(7) k=4 crystal (32 integers) · (8) array Xc (6-dim lattice) · (9) spectral/home (γ-axis).

| Berg quantity | 1 m-slice | 2 ladder | 3 Φ | 4 character | 5 Θ | 6 Macdonald | 7 crystal | 8 array | 9 spectral |
|---|---|---|---|---|---|---|---|---|---|
| wall / center 11 | window midpoint 22 = (edge+wall)/2, k-indep | **BLIND per-ε** (H_r centers wander at p=2r+1) | e₁=22 in-strip only | j=5 | b=0 by construction | ℓ+h∨ = 4+7 (C₆ peg); reads −2 in A₁₂⁽²⁾ | midpoint of {20,24} | **PROVABLY BLIND** (centroid 0, center off-centroid) | σ (the total watching) |
| anti-period (step −1) | F(m+44)=−F(m), in-wall only | none per rung | e₁±44: 0 hits globally | **BLIND** (finite table, translation leaves it) | GLOBAL on the lift F̂ (0/4050) — but constructed from the window | native genre (coroot lattice 2ℤ⁶, twisted (−,+)) | invisible | **PROVABLY ABSENT** (centroid theorem) | **exact & global**: veil e^{−2σt}, all t, machine zero |
| window [6+k,16−k] | [12+2k, 32−2k], wall = Yc coord-1 support wall | rung survival (short ladder = large m) | read-off strip e₂=2(k+5) | dominance strip λ₁+λ₂ ≤ 10 (B-side wall) | support = window exactly | s₀-closed alcove strip | 3-point stub [20,24] | support truncation of Xc | none yet (Attack C design) |
| C_k / P_k closed forms | five factored palindromic P_k(q) | **BLIND** (a_ε generic, no closed form) | border-orbit constancy | sparse symmetric pairs {j,10−j}, equal values | P₀=−(C₁₀+C₆) … P₅=0 | refuted as denominator | rows ARE the theorem | none | none yet |
| Θ = theta − wall-shadow | F = Θ−W (0/402) | — | — | symmetrized minus truncated | native | theta-with-characteristic target | — | absent | shadow analog OPEN |
| defects R_k | off-window support | defect rank d = 3,2,2,1,1 (depth 9) | — | unpaired below-window points | W_k wall-shadow, mirror past wall | — | k=4 defect-free | 26/225 off-window violations | OPEN |
| break atlas (l=2; T1/T2/T3) | three-row read-off breaks 8/16 | full-sum defect rank 2, per-ε 1/32 | — | **native**: f_τ slices; T1=39/T2=1/T3=14 | no uniform completion center (22/18/16) | must-break guard | M₁ ≠ 0 at l=2 | — | middle-seat odd/even (design) |
| crystal 1=1 | F₄(20)=F₄(24)=1 | 32 combs, ONLY-FULL-SUM | — | n\_(4,4)=n\_(6,4)=1 | — | — | native: 32 integers to ±305 | 64 distinct orbits, block-diagonal | — |
| counting facts (31−9=22, 14=2h∨) | — | — | — | — | — | level-pin breadcrumb | — | sheet counts | σ ↔ ℓ+h∨ (the level = the total watching) |

## 2. Switching question (i): the blind spots, and what they localize

- **Canvas 8 (array) is blind to the center BY THEOREM** (centroid obstruction): the one
  canvas holding the full object cannot even state where the second mirror stands.
- **Canvas 2 (ladder) is blind per-ε** (centers wander); the center exists only for the sum.
- **Canvas 4 (character) is blind to the translation**: the finite table never shows the
  anti-period; it sees only the reflection's footprint.
- **Canvas 9 (spectral/home) is the ONLY canvas where the translation is a global exact
  operator**: ρ_anti(t) = e^{−2σt}·ρ_gain(t) at every t, machine zero. Canvas 5's globality
  (the F̂ lift) is BOUGHT (constructed from the window), not seen.

Localization: every array-side canvas either lacks the center, lacks the translation, or
holds both only in-window. The two canvases in which the full dihedral structure is
unbroken are the theta lift (5, by construction) and the home dynamics (9, by nature).
That is the plan's attacks B and C, now derived from the blind-spot map rather than assumed.

## 3. Switching question (ii): joint relations no single canvas states

- **The ρ-condition** (§0): the pair of mirrors is visible iff the coordinate carries the
  half-sum shift. Unshifted canvases see a palindrome; shifted canvases see the dihedral.
- **Resolvent ↔ veil.** Canvas 2's telescope X = Y·2i·Σ_r M^{−(2r+1)} is a geometric sum in
  the coroot translation M⁻²; canvas 9's veil e^{−2σt} is the exponential of the generator
  shift −2σ·Id. Side by side: the SAME translation that on the array side exists only as a
  divergent-looking summed resolvent (and holds only in-wall) is, on the home side, the
  exact global exponential of the price. The joint statement no single canvas makes:
  *the translation is exact wherever the object is untruncated; the deletion (1⁶ sheet) and
  the wall are one truncation, and the window law is what survives of an exact translation
  after truncating.* On the home side nothing is deleted (the bridge is always open), which
  is WHY the law is global there.
- **Consistency cross-check (banked):** the dem f-cache reproduces the cry fence moments
  exactly after clean-window restriction (M₁ = −8, 12, −6 at μ₂ = 7, 8, 9) — computed
  independently in `_rb_planA_tiers.py` (μ-units, full dominant support: −22, 20, −24;
  drop the below-window point and the cry values reappear). Two canvases, one number.

## 4. Switching question (iii): the l=2 break in all canvases — THE FIND

`_rb_planA_tiers.py` (control: reproduces the full atlas census T1=39/T2=1/T3=14 from the
cache exactly). The famous 8 three-row l=2 break weights = 4 broken pairs = 3 break-tails
τ = (μ₂,6,3,2,1), and they are **tier-mixed**:

| μ₂ | tier | broken pairs (μ₁) | in n-units | cry-canvas M₁ |
|---|---|---|---|---|
| 7 | **T1** (support truncated, endpoints (8,12)) | (8,14), (10,12) | — | −8 |
| 8 | **T1** (endpoints (9,11)) | (9,13) | — | 12 |
| 9 | **T3** (support set symmetric {10,12}) | (10,12) | n\_(4,4,2) = 4 ≠ 1 = n\_(6,4,2) | −6 |
| 10 | hold (trivial: only center dominant) | — | — | — |

So the l=2 fence was never homogeneous: two of its three tails are the geometric face
(truncation), one is the deep face (pure coefficient asymmetry) — **T3 reaches down into
the three-row slice**, at exactly one weight pair, λ = (4,4,2) ↔ (6,4,2).

**And it is the minimal crystal's own shadow.** Same tail family μ₂ = 9, same pair
(μ₁, 22−μ₁) = (10,12), same m-points (20,24):

- two-row (l=0): n\_(4,4) = n\_(6,4) = 1 — the minimal crystal, the cleanest HOLD;
- three-row (l=2): n\_(4,4,2) = 4 ≠ 1 = n\_(6,4,2) — the minimal T3 BREAK, both partners
  alive, support sets identical, only the values differ.

The mechanism question now has its smallest possible discriminating instance: whatever seam
identity proves 1 = 1 at (4,4)/(6,4) must fail at (4,4,2)/(6,4,2) *without* being able to
blame support truncation. Every candidate from Attack B gets fenced on THIS pair first.

**The crystal canvas of both, side by side** (`_rb_planA_crystal_pair.py`, 32 sheets):

- hold: Σ a_ε = 0, 30/32 nonzero, range [−182, 305];
- break: Σ a_ε = 3 (= 4−1), 26/32 nonzero, range [−76, 122];
- **sign-coherence**: at every ε where both are nonzero, sign(a_hold) = sign(a_break)
  (26/26); the break's zero set contains the hold's; no constant ratio (magnitudes 1.6×–7.4×).

Descriptive, not yet mechanism: the T3 break is a sign-coherent shrinkage of the hold's
flip-sheet decomposition — the l=2 resonance deforms every sheet's magnitude while
preserving every sheet's sign, and the exact cancellation is what dies. (Caveat: the two
tails differ in one coordinate, so adjacent-coefficient correlation is expected; the
26/26 sign match is stronger than that expectation but is not yet a law. If a candidate
identity emerges, test whether it explains the sign-coherence for free.)

## 5. Deliverable state / hand-off

- Table complete (§1), blind spots localize the mechanism to theta-lift + home dynamics
  (§2) — Attack B and C confirmed as the right next moves, now by exclusion.
- Two exact scouts banked: tier map (`_rb_planA_tiers.py`) and crystal pair
  (`_rb_planA_crystal_pair.py`).
- NEW facts for the atlas: (a) the famous 8 = 2×T1 + 1×T3 tails; (b) T3's minimal instance
  is three-row λ=(4,4,2)↔(6,4,2), the l=2 shadow of the minimal crystal; (c) sign-coherent
  shrinkage in the 32-sheet decomposition; (d) the ρ-condition (center = half-step iff the
  coordinate is ρ-shifted); (e) dem-cache ↔ cry-fence moment agreement (two canvases, one
  number).
- Fence discipline for Attack B, sharpened: a candidate theta identity must (1) hold on the
  14 live pairs, (2) break at the two T1 tails BY truncation, and (3) break at (4,4,2)
  WITHOUT truncation — condition (3) is new and is the hard one.
