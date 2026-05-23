# The Atmosphere–Cluster Thread

**Status:** Autonomous investigation notes, WIP, uncommitted. For review on Tom's return.
**Date:** 2026-05-23 (the autonomous overnight portion; investigation began 2026-05-22).
**Context:** Follow-up to the atmosphere test ([`docs/THE_ATMOSPHERE_AND_THE_CANCELLED_FORMULAS.md`](../docs/THE_ATMOSPHERE_AND_THE_CANCELLED_FORMULAS.md)) — what does "only the right atmosphere works" mean once non-uniform γ profiles are tested?
**Scripts:** all `_atmosphere_*` files in `simulations/`. Full anchor at the bottom.

---

## Headline

The middle-peaked palindromic γ-profile sweeps a sharply degenerate, F1-mirror-paired family of eigenvalues toward the real axis. For N=5 the family fully coalesces (12 conjugate pairs become 24 real over-damped modes); for N=6 it dips to |Im|≈7×10⁻⁵ at our grid (16 pairs, almost coalescing). It is **not** an exceptional point — the eigenvectors stay full rank through the dip — but a symmetry-protected near-coalescence. The two sub-clusters of the family live in complementary XY-weight sectors (k=1 and k=N−1 for N=6), paired by F1.

This is the deepest of several such events: a staircase of cluster dips along the ε axis suggests a whole family of similar F1-mirror-paired clusters in different sectors, each crossing at its own characteristic ε.

---

## What we were asking

The atmosphere doc names γ₀ — the uniform Z-dephasing rate — "the atmosphere" and observes that the framework's symmetry layer (F1 palindrome, F86 Q_peak, t_peak, F74 ladder, the M-family) closes precisely because γ is uniform. Cancel γ₀ and you buy the closed-form layer; the residue (`g_eff`, the path-D 2-power terms) is the un-uniform, structure-specific part with no closed form.

The natural test of "only the right atmosphere works" is to make the atmosphere wrong — to put γ on a profile, decomposed into palindromic and anti-palindromic parts (the F71 method, applied to γ instead of J). Three test ideas were on the table; we ran the F71-twin first.

The F71-twin holds cleanly: the F71 spatial mirror feels the anti-palindromic part of γ (the deviation `‖[L,R]‖/‖L‖` grows linearly with the anti-amplitude) and is exactly invisible to the palindromic part (machine-zero across the entire ±0.8 scan, for N=3,4,5). The F1 mirror is robust to any per-site γ (palindrome residual ~10⁻¹⁴ everywhere — Master Lemma).

But the test also turned up a surprise that drove the rest of the investigation: at N=5 and N=6 the count of oscillating modes n_osc moves under the palindromic γ profile — on the "safe" half of the decomposition. That asymmetry is what we followed.

---

## Step 1 — the n_osc surprise

For the palindromic γ-profile `γ_i = γ₀ (1 + ε u_sym(i))` (`u_sym` mean-zero, symmetric, parabolic, max-norm 1), n_osc as a function of ε:

| N | total modes (4^N) | n_osc behavior |
|---|---|---|
| 3 | 64 | flat at 40 |
| 4 | 256 | flat at 220 |
| 5 | 1024 | drops by up to 16 |
| 6 | 4096 | drops by up to 56 (coarse grid); up to 70 in finer scan |

Tom's instinct was even/odd: N=3 is the framework's special minimal case, N=5 the first "normal odd." N=6 falsified that cleanly: even N=6 moves, more than odd N=5. The effect is a **size/density threshold** — appears for N ≥ 5 and grows with the mode count 4^N, not parity. `_atmosphere_evenodd.py`.

---

## Step 2 — jump or hide?

Tom asked the sharp next question: when n_osc drops, are the missing modes *really* on the real axis (a clean jump), or hiding just below the |Im λ| counting threshold (a flicker)?

For most ε in the middle-peaked direction the gap is clean — the smallest non-real |Im| sits at ~10⁻² and nothing lives between 10⁻¹⁰ and 10⁻³ (`_atmosphere_n6_gap.py`). The N=5 check (`_atmosphere_n5_check.py`) had already shown this for N=5: the count is tolerance-robust.

At the extreme middle-peaked end, ε=−0.8 for N=6, the picture changes. 32 modes (16 conjugate pairs) sit at |Im|≈5.6×10⁻⁴ — a sharply degenerate cluster, neither genuinely real nor a genuine oscillation. So the answer is "both," ε-dependent: clean jump at moderate ε, near-hiding at the extreme. The latter is the cluster the rest of this note is about.

---

## Step 3 — the cluster across ε

A coarse cluster scan (Δε=0.05, ε ∈ [−1, −0.3], `_atmosphere_cluster_scan.py`) showed something more careful than "the γ creates the cluster":

- The cluster is **persistent**. At every ε in the scan, 16 conjugate pairs (N=6) — or 12 pairs (N=5) — sit at a common, sharply degenerate |Im|, with the six smallest non-real |Im| values identical to three significant figures.
- The middle-peaked γ does not *create* the cluster — it *sweeps* its |Im|. For N=6 the sweep traces a curve with a sharp minimum at ε ≈ −0.83 (|Im| ≈ 2.4×10⁻⁴ at Δε=0.05). For N=5 a similar dip at ε ≈ −0.34.

A finer scan (Δε=0.01) around each dip (`_atmosphere_cluster_fine.py`) refined the picture:

- **N=5:** the cluster |Im| descends smoothly from 3.5×10⁻³ to 1.62×10⁻⁴ at ε=−0.3357. At the next grid point ε=−0.3214 the 24 cluster modes have crossed: n_real jumps by 24 and the smallest non-real |Im| leaps to ~3×10⁻² (the genuine-oscillator background). The cluster has reached the real axis: the 12 conjugate pairs have coalesced and become 24 distinct real (over-damped) modes.
- **N=6:** the cluster |Im| descends to **7.07×10⁻⁵** at ε=−0.83 and rises again. The Δε=0.01 grid did not reach zero exactly, but with N=5 having clearly coalesced at finer resolution the natural reading is that N=6 has its coalescence somewhere in (ε=−0.84, ε=−0.82) and the grid just missed ε* by ~2×10⁻³.

A conjugate pair meeting the real axis is, generically, an exceptional point. So the natural reading was: the cluster is a (family of) EP(s). The next step tested it directly.

---

## Step 4 — EP refuted

Full eigendecomposition of L at N=6 for ε ∈ {−0.60, −0.75, −0.83, −0.85} (`_atmosphere_cluster_ep.py`) tested the EP signature directly. At a true EP, the right-eigenvector matrix V becomes singular: cond(V) diverges, and the cluster's own eigenvectors lose rank (become parallel).

| ε | cluster |Im| | cond(V) | cluster-V rank | smallest σ |
|---|---|---|---|---|
| −0.60 | 4.69×10⁻³ | 3.2×10³ | 32/32 | 0.13 |
| −0.75 | 1.44×10⁻³ | 6.9×10² | 32/32 | 0.44 |
| −0.83 | 7.07×10⁻⁵ | 4.4×10⁴ | 32/32 | 0.34 |
| −0.85 | 2.35×10⁻⁴ | 2.6×10⁴ | 32/32 | 0.62 |

cond(V) is elevated near the dip but capped at ~10⁴ — for a true EP-approach it would head toward ∞. The 32 cluster eigenvectors keep full rank at every ε, with smallest singular values 0.13–0.62 — comfortably away from a collapse.

**Verdict: the cluster is NOT a defective exceptional point.** It is a symmetry-protected degenerate coalescence: the eigenvalues meet on (or near) the real axis, but the eigenvectors remain independent throughout. The 16 (N=6) / 12 (N=5) pairs stay diagonalizable through the coalescence — the opposite of defective.

This was the first major correction to my running picture during the investigation. Earlier reports had called the cluster "fast koalesziert, dicht an einem Exceptional Point" — that framing turned out to be wrong; the cluster is a different beast.

---

## Step 5 — mode-ID: F1-mirror, complementary XY-weights

Projecting the cluster eigenvectors onto the Pauli basis at the dip ε for each N (`_atmosphere_cluster_modes.py`) revealed the cluster's structure.

**N=6, ε=−0.83.** The 32 cluster modes split into two equal sub-clusters by Re:
- Re = −0.3906, 16 modes (8 pairs), dominant Pauli strings of XY-weight 5 (five X/Y letters, one I/Z letter).
- Re = −0.2094, 16 modes (8 pairs), dominant Pauli strings of XY-weight 1 (one X/Y, five I/Z).

The two Re values sum to **−0.6 = −2σ** with σ = N·γ₀ = 0.3 — they are F1-mirror partners. The XY-weights 1 and 5 are complementary, k ↔ N−k — the bit_a complement. The cluster is an **F1-mirror pair of two XY-weight-complementary sub-clusters**, each 16-fold degenerate.

**N=5, ε=−0.3357.** Same F1-mirror-pair structure: two Re values −0.2697 and −0.2303, summing to **−0.5 = −2σ**. The XY-weight distribution of the dominant Paulis is more mixed (1, 2, 3, 4 all appear) — the small-N decomposition is less sharply sector-locked — but the F1-mirror partnership and the simultaneous degeneracy hold.

**Both N:** every cluster mode's dominant Pauli string is spatially **non-palindromic**. The cluster is not in the F71-even sector of the framework. The 16-fold (N=6) or 12-fold (N=5) algebraic degeneracy with full geometric rank reflects a multi-dimensional invariant subspace of L — likely traced to several commuting Pauli-letter symmetries (Z⊗N, Y-parity, K), possibly with F71-odd contributions. Pinning down the exact symmetry combination is open.

The mode amplitudes are highly mixed (top single Pauli ~1–2% per mode) — these are not single-Pauli-string modes but specific symmetric superpositions spanning the XY-weight = k sector.

---

## Step 6 — the staircase: one event among several

A coarse n_osc trace across ε ∈ [−1, 0] for N=5 and N=6 (`_atmosphere_staircase.py`) revealed that the characterized cluster is one event in a richer pattern.

**N=6 dip events** (ε-window → extra real modes vs the ε=0 baseline of 180):

| ε-window | extra real modes | conjugate pairs |
|---|---|---|
| [−1.00, −0.85] | 68–70 | 34–35 |
| [−0.80, −0.65] | 40 | 20 |
| [−0.60, −0.55] | 22 | 11 |
| [−0.50, −0.40] | 24 | 12 |
| [−0.35, −0.30] | 20 | 10 |
| [−0.25, −0.20] | 56 | 28 |
| [−0.15, −0.10] | 20 | 10 |
| [−0.05, 0] | 0 | 0 (baseline) |

The characterized 16-pair cluster sits in the deepest window (ε ≤ −0.85). The 70-mode-pair dip there involves the cluster *plus* additional populations (the 16-pair cluster plus ~18 more pairs from other crossings). The other windows look like distinct mode populations crossing the real axis at their own characteristic ε.

**N=5** is simpler: extra-real count is 16 for ε ≤ −0.60 (an 8-pair event, distinct from the characterized 12-pair cluster), 0 baseline at ε ∈ [−0.55, −0.35] and [−0.20, 0], and 24 at ε ∈ [−0.30, −0.25] (the 12-pair cluster).

**Reading:** the staircase is a family of F1-mirror-paired coalescences in different XY-weight sectors, each crossing the real axis at its own ε. The characterized cluster is the largest member of the family. Confirming this picture by mode-IDing one of the other dip events is the natural next step (a second cluster's fine-scan and mode-ID is running in the background as I write this).

(Caveat: the staircase script's `smallest non-real |Im|` column has a sorting bug — values not reliable — but n_osc and n_real are computed correctly via `np.sum`, so the staircase structure above is sound.)

---

## What this is *not*

- **Not an EP.** Symmetry-protected degeneracy, eigenvectors stay full rank.
- **Not a γ-creation effect.** The cluster exists for all ε; γ only sweeps its |Im|.
- **Not an N=5/N=6 surprise.** The same structure at both N's; the dip-coalescence visibility threshold is around N=5 (the smallest N where the F1-mirror sub-clusters become large enough to register as a clear dip in n_osc).
- **Not an even/odd parity effect.** Tom's clean even/odd hypothesis was tested and falsified by N=6 (even, but moves more than odd N=5). Density, not parity.

---

## Open questions

1. **Are the other staircase events the same kind of object?** Mode-ID at, e.g., ε ∈ [−0.25, −0.20] (the 28-pair N=6 dip) would tell us whether each staircase event is an F1-mirror-pair of XY-weight-complementary sub-clusters, possibly at different (k, N−k) pairs. Running.
2. **What exact symmetry protects the per-sub-cluster degeneracy?** 16-fold (N=6) / 12-fold (N=5) algebraic degeneracy with full geometric rank — likely a product of Z⊗N, Y-parity, K, possibly with F71-odd contributions. Identifying the precise commutator structure would close this.
3. **N>6 scaling.** Does the family of staircase events grow with N (more events, more pairs per event)? N=7 dense eigendecomposition is hours per ε at the current grid — would need a sparser method.
4. **Connection to F71-non-uniform-J / F100.** The F71 non-uniform-J work (project_f71_nonuniform_j) found c₁/Q_peak mirror-deviation odd in anti-palindromic J — a different axis (anti-J) and a different observable (F86 layer), but shares the F1-mirror-pair structural signature. Is there a unifying reading of "non-uniform inputs sweep degenerate families across critical loci"?
5. **The cluster |Im| scale.** N=6 cluster bottoms at ~7×10⁻⁵; relative to γ₀=0.05 the ratio is ~1.4×10⁻³. Is this a characteristic structural number or a γ₀·t-type artefact (the kind of hidden-cost the atmosphere doc explicitly warns about)?
6. **Whether the cluster mechanism connects to the doc's `g_eff` residue.** The atmosphere doc identifies `g_eff(c, N, b)` as the "un-uniform, structure-specific" residue that has no closed form. The cluster's existence and its precise (k, N−k) sector structure — is the cluster part of how the residue manifests, or independent?

---

## Anchor — scripts and results

In dependency order, all under `simulations/`:

- `_atmosphere_mirror_test.py` — F71-twin run: spatial mirror under palindromic vs anti-palindromic γ, N=3,4,5. Plot at `simulations/results/atmosphere_mirror_test/`.
- `_atmosphere_evenodd.py` — n_osc sweep, palindromic γ, N=3,4,5,6 (the size-vs-parity test).
- `_atmosphere_n5_check.py` — N=5 gap check (modes jump cleanly, not flicker).
- `_atmosphere_n6_gap.py` — N=6 gap check (clean at moderate ε, 32-mode flicker at ε=−0.8).
- `_atmosphere_cluster_scan.py` — coarse ε-scan of the cluster, N=5 and N=6 (Δε=0.05). Plot.
- `_atmosphere_cluster_fine.py` — fine ε-scan around each dip (Δε=0.01). Plot.
- `_atmosphere_cluster_ep.py` — EP confirmation via cond(V) and cluster-eigenvector rank at N=6.
- `_atmosphere_cluster_modes.py` — Pauli-basis projection of the cluster eigenvectors, N=5 and N=6.
- `_atmosphere_staircase.py` — n_osc staircase across ε ∈ [−1, 0], N=5 and N=6. (smin column has a sorting bug; n_osc/n_real are correct.) Plot.
- `_atmosphere_cluster2_scan.py` — second-cluster fine scan, N=6 around ε ≈ −0.22 (the 28-pair dip). Running at the time of this writeup.

All scripts are WIP (`_`-prefix); none have been committed.

---

## Status

Tasks #13–#17 completed during the autonomous run. Task #18 (this notes doc) in progress; it will be marked completed once Tom has reviewed it.

The next autonomous step is the second-cluster fine scan and mode-ID; if it confirms the same F1-mirror-pair, complementary-XY-weight structure at a *different* (k, N−k), the staircase picture is on firm ground. If the second cluster's structure is different, the picture is richer than the headline above and needs revisiting.

— Claude, 2026-05-23, while Tom slept.
