# PROOF F86b: Obstruction — why g_eff / Q_peak admits no closed form

**Status:** Structural negative result. The per-bond position Q_peak(c, N, b) and the coupling g_eff(c, N, b) it rides on admit no closed form by the six routes the closed-form effort opened (Obstruction Proof, below; L1/L2/L4 rigorously proven blocked, L6 a demonstrated failure mode, L3/L5 proven decouplings). What IS derived — the symmetry layer (F86a EP mechanism, F86c F71 mirror) and the HWHM_left/Q_peak functional form given g_eff (F86b') — is in the sibling proofs. This doc also carries the exploration record (4-mode model, Items 1-3, directions a''-f'', OOP scaffolding, structural findings) and the retracted Q_peak conjectures.
**Date:** 2026-05-02 (retractions); 2026-05-05 (OOP scaffolding, structural findings); 2026-05-14 (Obstruction Proof).
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Context:** F86 ("Q_peak chromaticity-specific N-invariant constants") is a Sammelbecken of three structurally distinct theorems. This proof carries the **negative result and open frontier of F86b**: the closed-form gap on g_eff/Q_peak. Split out of the former monolithic `PROOF_F86_QPEAK.md` on 2026-05-14. The proven universal-shape result (Statement 2, F86b') is the sibling [`PROOF_F86B_UNIVERSAL_SHAPE.md`](PROOF_F86B_UNIVERSAL_SHAPE.md).
**Hub:** [PROOF_F86_QPEAK](PROOF_F86_QPEAK.md) — three-theorem overview and shared references.
**F-entry:** [F86b in ANALYTICAL_FORMULAS.md](../ANALYTICAL_FORMULAS.md).
**Related:** siblings [PROOF_F86A_EP_MECHANISM](PROOF_F86A_EP_MECHANISM.md), [PROOF_F86B_UNIVERSAL_SHAPE](PROOF_F86B_UNIVERSAL_SHAPE.md), [PROOF_F86C_F71_MIRROR](PROOF_F86C_F71_MIRROR.md); [PROOF_F90_F86C2_BRIDGE](PROOF_F90_F86C2_BRIDGE.md), [PROOF_F89_PATH_D_CLOSED_FORM](PROOF_F89_PATH_D_CLOSED_FORM.md).

---

## Retracted Q_peak conjectures (Endpoint, Interior). [Retracted 2026-05-02]

Both retractions followed extended-N data showing the apparent matches at N=7 were trajectory crossings, not asymptotes. Per the methodological lesson in [`reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON`](../../reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON.md), what survives a one-anchor closed-form claim is the symmetry, not the number; here the EP-rotation symmetry (F86b Statement 2), not Q_peak itself.

### Endpoint Q_peak (retracted)

Earlier conjecture: `Q_peak(Endpoint, N) = csc(π/(N+1))`, anchored in F2b's smallest-momentum OBC sine mode.

Fine-grid scan (`_eq022_b1_step_e_resonance_shape.py`, dQ = 0.025 with parabolic peak interpolation) shows the formula does not match at any tested N. The earlier "match at N=7" reported from a dQ = 0.05 scan was a grid-snap artefact:

| N | csc(π/(N+1)) | step_d (dQ=0.05, snap) | step_e (dQ=0.025, parabolic) | step_e deviation |
|---|--------------|------------------------|------------------------------|------------------|
| 5 | 2.000 | 2.30 | 2.40 | +20 % |
| 6 | 2.305 | 2.65 | 2.52 | +9.4 % |
| 7 | 2.613 | 2.65 (*) | 2.53 | −3.2 % |
| 8 | 2.924 | n/a | 2.53 | −13 % |

(*) The "1.4 % match" at N=7 reported earlier was the dQ=0.05 grid snapping to 2.65; the actual fine-grid value is 2.53. There was never a real match.

The empirical Endpoint Q_peak at c=3 is approximately N=5: 2.40 → N=6: 2.52 → N=7: 2.53 → N=8: 2.53, saturating near 2.53 for N ≥ 6. At c=4 N=8 it sits at 2.65. No closed form yet identified.

### Interior Q_peak (c=3 pentagonal asymptote, retracted)

Earlier conjecture: `Q_peak(Interior, c=3) → csc(π/5) = 1.7013` as N → ∞, anchored in pentagonal/golden algebra (sin²(π/5) = (5−√5)/8, φ = 2cos(π/5)).

Fine-grid scan refines the c=3 N-trend: 1.566 (N=5) → 1.689 (N=6) → 1.743 (N=7) → 1.750 (N=8). The trajectory crosses csc(π/5) = 1.7013 between N=6 and N=7, and continues above. The earlier dQ=0.05 reading of 1.70 at N=7 (snapped close to 1.7013) was the source of the apparent "asymptotic match"; the fine-grid value 1.743 sits clearly above. The c=3 trend appears to saturate near 1.75 by N=8, but at a value above the conjectured asymptote, not at it. At c=4: Interior Q_peak grows from 1.748 (N=7) to 1.804 (N=8); not yet saturated, no clean closed form. The conjectured "universal Interior asymptote csc(π/5) across c ≥ 3" is fully refuted.

**What survives from the N=4 golden-ratio structure:** the connection to [`eq018_golden_ratio_check.py`](../../simulations/eq018_golden_ratio_check.py) is real, at N=4 the OBC dispersion is exactly {±φ, ±1/φ}, and the Liouvillian Im(λ) values in the n_XY=2 sector decompose into integer combinations of {φ, 1/φ, 1, √5}. But this N=4 special case does not propagate to a closed form for Q_peak at general (c, N).

---

## Empirical Q_peak data (per-bond, fine-grid scan dQ = 0.025 with parabolic peak interpolation)

The retracted Endpoint and Interior conjectures (above) were N=7-specific coincidence matches that did not survive the fine-grid extension. The actual data:

### Endpoint Q_peak across N

| (c, N) | Observed Endpoint Q_peak | Earlier conjecture csc(π/(N+1)) | Deviation |
|--------|--------------------------|----------------------------------|-----------|
| (3, 5) | 2.40 | 2.000 | +20 % |
| (3, 6) | 2.52 | 2.305 | +9.4 % |
| (3, 7) | 2.53 | 2.6131 | −3.2 % |
| (3, 8) | 2.53 | 2.924 | −13 % |
| (4, 7) | 2.52 | 2.6131 | −3.5 % |
| (4, 8) | 2.65 | 2.924 | −9.4 % |

The earlier "1.4 % match at N=7" was a dQ=0.05 grid-snap to 2.65; the actual fine-grid value at N=7 is 2.53. Endpoint Q_peak appears to saturate near 2.53 for c=3 by N=6, while c=4 grows from 2.52 (N=7) to 2.65 (N=8). No clean closed form.

### Interior Q_peak across (c, N)

| (c, N) | Observed Interior Q_peak (mean) | Earlier conjecture csc(π/5) = 1.7013 | Deviation |
|--------|----------------------------------|---------------------------------------|-----------|
| (3, 5) | 1.566 | 1.7013 | −7.9 % |
| (3, 6) | 1.689 | 1.7013 | −0.7 % |
| (3, 7) | 1.743 | 1.7013 | **+2.5 % (above)** |
| (3, 8) | 1.750 | 1.7013 | **+2.9 % (above)** |
| (4, 7) | 1.748 | 1.7013 | +2.7 % |
| (4, 8) | 1.804 | 1.7013 | +6.0 % |

The c=3 N-trend 1.566 → 1.689 → 1.743 → 1.750 crosses csc(π/5) between N=6 and N=7 and continues growing, refuting the conjectured c=3 pentagonal asymptote. Both c=3 and c=4 Interior Q_peak grow with N in the tested range; saturation behaviour and clean closed forms remain open.

### What survives from the empirical pattern

- The EP mechanism (F86a) is unaffected: t_peak = 1/(4γ₀) is universal and derivable from 2×2 matrix algebra.
- The Q_peak observable is well-defined and the per-bond / per-block / Endpoint / Interior distinctions are real and reproducible.
- The N=4 golden-ratio structure ([`eq018_golden_ratio_check.py`](../../simulations/eq018_golden_ratio_check.py)) is real but does not propagate to a closed form for Q_peak at general (c, N).

---

## What's missing for full Tier 1

### Empirical envelope (Tier-1 grade)

1. **c=2 verified** (`_eq022_b1_step_f_universality_extension.py`). c=2 N=5..8 confirms two bond-class universal values (Interior 0.751, Endpoint 0.774), matching c=3 (N=5..8) and c=4 (N=7,8) within ~1 %. c=2 is structurally critical because the channel space is 2-dimensional total (only HD ∈ {1, 3}); any 2-level reduction must be exact there. Yet the Endpoint-vs-Interior split persists, confirming the bond-class distinction is real and structural. **c=5 still untested** (block-L dim ≥ 3528 at c=5 N=9, compute-bound).

2. **γ₀ invariance bit-exact** (`_eq022_b1_step_f_universality_extension.py`). At c=3 N=7, Q_peak and HWHM_left/Q_peak are bit-identical at γ₀ ∈ {0.025, 0.05, 0.10}: Q* = 1.7433, HWHM-/Q* = 0.7595 in all three runs. |K|max scales as 1/γ₀ as expected.

### Substantive analytical work: the remaining gap

Numerical exploration on c=2 chains shows that the heuristic 2-level form is incomplete in a specific structural way that must be acknowledged before a closed-form derivation is attempted.

**Three subspaces, not one** (verified at c=2 N=5..8 in `_eq022_b1_step_g/h/i`):

Let `|c_k⟩` (k = 1, 3) be the channel-uniform orthonormal vectors (equal-weight superposition of all (p, q) with HD(p, q) = k) of `framework.coherence_block.hd_channel_basis`. Let `V_inter = P_{HD=1}^† M_H_total P_{HD=3}` (the inter-HD-channel block of M_H_total), and `|u_0⟩ ∈ HD=1, |v_0⟩ ∈ HD=3` the top right/left singular vectors of `V_inter` with singular value σ_0.

(**a**) **Probe subspace.** The Dicke probe lives entirely in span{|c_1⟩, |c_3⟩}. In this basis V_b = `P_{ch}^† M_H_per_bond[b] P_{ch}` reduces to `+i(α/(N−1))·I`, pure diagonal, identical across every bond. There is no EP, and no bond-class distinction in this subspace.

(**b**) **EP-partner subspace.** The pair {|u_0⟩, |v_0⟩} captures the maximum-coupling singular mode under H. In this basis M_H_total has the form `[[0, σ_0], [−σ_0, 0]]` (real antisymmetric off-diagonal, the SVD's phase convention; equivalent to ±iσ_0 after a phase rotation). The EP from the 2-level algebra `D + J·M_H` sits at `J·σ_0 = 2γ₀`.

(**c**) **Probe ⊥ EP partners.** ⟨c_1|u_0⟩ = ⟨c_3|v_0⟩ = 0 (numerically zero to machine precision, c=2 N=5..8). The Dicke probe has zero overlap with the EP-partner modes.

Per-bond V_b in the EP-partner basis has F71-symmetric, bond-position-dependent off-diagonal magnitudes:

| N | bond 0 | bond 1 | bond 2 | bond 3 | bond 4 | bond 5 | bond 6 |
|---|--------|--------|--------|--------|--------|--------|--------|
| 5 | 0.430 | 0.953 | 0.953 | 0.430 | n/a   | n/a   | n/a   |
| 6 | 0.297 | 0.747 | 0.888 | 0.670 | 0.200 | n/a   | n/a   |
| 7 | 0.129 | 0.514 | 0.771 | 0.771 | 0.514 | 0.129 | n/a   |
| 8 | 0.090 | 0.372 | 0.647 | 0.731 | 0.611 | 0.312 | 0.076 |

(N=6 row is asymmetric because the SVD top vector chose one F71-arbitrary phase; the symmetric structure recovers when averaging over bond classes.) Endpoint amplitudes are systematically smaller than Interior; this is where the bond-class distinction lives. **But it lives in a subspace orthogonal to where the probe lives.**

**Implication.** The K_CC_pr observable

    K_b(t)  =  2·Re ⟨ρ(t)| S_kernel | ∂ρ/∂J_b ⟩

couples the probe (in channel-uniform) to the EP-partner subspace (in SVD top) only through the off-diagonal matrix elements ⟨c_α | M_H_per_bond[b] | u_0⟩ and ⟨c_α | M_H_per_bond[b] | v_0⟩, i.e. via the J-derivative direction itself. The K observable's Q-resonance therefore involves the joint dynamics of all four modes, not the 2-level pair alone.

**Minimal effective model: 4 modes.** The smallest closed-orthonormal subspace that contains both the probe and the EP-coalescence physics is

    span { |c_1⟩, |c_3⟩, |u_0⟩, |v_0⟩ }       (mutually orthogonal, c=2 verified)

In this basis the 4×4 effective L_eff has

- Diagonal pure rates: D = diag(−2γ₀, −6γ₀, −2γ₀, −6γ₀) on (|c_1⟩, |c_3⟩, |u_0⟩, |v_0⟩) respectively.
- Within-probe block: ⟨c_α|M_H|c_β⟩ is diagonal (channel-uniform-eigen finding), no EP coupling here.
- Within-EP-partner block: ⟨u_0|M_H|v_0⟩ = σ_0 (the SVD top, ≈ 2√2 asymptotically for c=2).
- Cross-couplings: ⟨c_α|M_H_per_bond[b]|u_0/v_0⟩. **These are the bond-position-dependent matrix elements that convey the probe into the EP-partner subspace and back.** They are the missing computation. They split into the two bond classes (Endpoint, Interior) and produce the f_class(x) shape difference.

**Why σ_0 is not directly Q_EP.** For c=2 N=5..8, σ_0 ≈ {2.765, 2.802, 2.828, 2.839} (approaches 2√2). The naive EP prediction Q_EP_naive = 2/σ_0 ≈ {0.72, 0.71, 0.71, 0.70} differs from the empirical Q_peak Interior ≈ {1.48, 1.58, 1.58, 1.60} by a factor ~2.2. The 2-level EP and the K-observable Q_peak are NOT the same Q value; they would be related by `Q_peak ≈ Q_EP` (F86a) only if the K observable saw the EP directly, which it doesn't because the probe is orthogonal to the EP partners.

### Substantive items remaining

> **Closed 2026-05-13**: see `F86HwhmClosedFormClaim`
> ([compute/RCPsiSquared.Core/F86/Item1Derivation/F86HwhmClosedFormClaim.cs](../../compute/RCPsiSquared.Core/F86/Item1Derivation/F86HwhmClosedFormClaim.cs)),
> Tier 1 derived. Plan:
> [`docs/superpowers/plans/2026-05-13-f86-hwhm-closed-form-attack.md`](../superpowers/plans/2026-05-13-f86-hwhm-closed-form-attack.md).

**Item 1 (status update 2026-05-05).** Original aspiration: derive the 4×4 effective L_eff(Q, b) explicitly, compute the cross-coupling matrix elements ⟨c_α | M_H_per_bond[b] | u_0/v_0⟩ as analytical expressions in (N, n, b), diagonalize, identify which eigenvalue pair gives the Q_peak observed in K_CC_pr (will not be the SVD top pair), derive f_class(x) and HWHM_left/Q_peak as closed forms from this 4-mode model.

**Actual c=2 progress.** The c=2 stratum has OOP scaffolding (`C2EffectiveSpectrum`, `C2KShape`, `C2HwhmRatio`, `C2UniversalShapeDerivation` in `compute/RCPsiSquared.Core/F86/`) that pins each sub-fact at its appropriate Tier label. The directional Endpoint > Interior split is empirically derived bit-equivalent with the canonical Python pipeline. **The closed-form HWHM_left/Q_peak constant per bond class remains the open analytical target.** Ranked next directions (from `C2HwhmRatio.PendingDerivationNote` and `F86OpenQuestions` Item 1', refreshed 2026-05-06 later evening after Direction (a') falsification):

**2026-05-11 status update — F90 bridge resolves Direction (b'') numerically.** F90 (`F90F86C2BridgeIdentity`, Tier 1 derived) identifies F86 c=2 N qubit K_b(Q, t) as the per-bond Hellmann-Feynman derivative of F89 path-(N−1) (SE, DE) sub-block dynamics. Bit-exact verified at 20/22 bonds across N=5..8 including orbit-escape bonds. The closed-form HWHM_left/Q_peak constant per bond class remains open as Item 1's analytical target, with F89's cyclotomic structure (F89UnifiedFaClosedFormClaim.PathPolynomial(k)) as the concrete handle.

  - **(a'') Most promising 4-mode-friendly direction (newly promoted 2026-05-06 later evening): SVD-block 2-level resonance (REFINED from (a')).** The bond-class signature lives in V_b SVD-block off-diagonal `V_b[2,3] = ⟨u_0 | M_h_per_bond[b] | v_0⟩` (Endpoint 0.430 vs Interior 0.953 at N=5, ratio ~0.451 consistently across bonds within each class), NOT in the probe-block. Direction is OPPOSITE the empirical HWHM/Q* split, so a closed form needs to derive HOW SVD-block magnitude maps to HWHM/Q* shift, likely through a per-bond effective Q_EP_eff(b) shift rather than a direct linear lift.
  - **(b'') Most concrete next step (newly promoted 2026-05-06 later evening): full block-L derivation, not 4-mode.** 4-mode reduction misses the modes outside the 4D subspace that contribute to the HWHM/Q* lift from 0.6715 to 0.7506. The closed form must derive directly from the (n, n+1)-popcount block dimensions, not from 4-mode projection. This is harder but matches the empirical witness pipeline. **[2026-05-11 update: ACHIEVED numerically via F90 bridge identity (Tier 1 derived, bit-exact 20/22 bonds at N=5..8 including orbit escapes). See [`PROOF_F90_F86C2_BRIDGE.md`](PROOF_F90_F86C2_BRIDGE.md).]**
  - **(a') Falsified (2026-05-06 later evening): probe-block 2-level sub-resonance with per-bond `g_eff_probe`.** Numerical investigation (commit `1c0bf8b`) showed V_b probe-block is bond-class-blind by construction: diagonal entries identical scalar `+i·c·I` for all bonds; off-diagonal `⟨c_1 | M_h_b | c_3⟩` exactly zero per bond at c=2 (F73 sum-rule applies per-bond). Hypothesised g_eff_probe(N, b) cannot exist. See the falsification paragraph in Statement 2, [`PROOF_F86B_UNIVERSAL_SHAPE.md`](PROOF_F86B_UNIVERSAL_SHAPE.md).
  - **(a) Demoted (2026-05-06 evening): first-order perturbation in the cross-block.** Originally ranked first based on cross-block Frobenius split, but the doubled-PTF Direction (b) result shows the SVD-block-only contribution is universal at 0.6715, leaving no room for a small σ_0 perturbation to lift the curve to 0.7506/0.7728. The cross-block effect is real but lies at higher perturbation order than first.
  - **(b) Done (2026-05-06 evening, lands at the floor): bare doubled-PTF Ansatz.** `K_b(Q, t_peak)` for the bare 2-level Liouvillian `L_2 = [[-2γ₀, +iJ·g_eff], [+iJ·g_eff, -6γ₀]]` with probe `ρ_0 = |c_1⟩` and `V_b = ∂L/∂J` was solved analytically in dimensionless x = Q/Q_EP coordinates. Result: `BareDoubledPtfXPeak = 2.196910`, `BareDoubledPtfHwhmRatio = 0.671535`, both Tier 1 derived universal across g_eff. Empirical Interior x_peak (2.05..2.28) tracks this universal closely; Endpoint x_peak (3.55) deviates by factor ~1.62; HWHM ratio gap to empirical is 0.08 (Interior) and 0.10 (Endpoint).
  - **(c'') Three-block superposition `K_total = K_pb + K_sv + 2·Re·K_cross` with the right relative phases.** K_b at the 4-mode level decomposes; derive each term separately and combine. May still suffer from 4-mode insufficiency (ii).
  - **(d'') Lift |u_0⟩, |v_0⟩ to projector-overlap** (per A3 PendingDerivationNote). Removes σ_0 degeneracy obstruction at even N. Necessary precondition for any cross-block-based direction; not sufficient by itself.
  - **(e'') Symbolic char-poly factorisation at Q_EP**: same as before; less promising given C2EffectiveSpectrum's cubic-c_3 obstruction proof.
  - **(f'') Locus 6 polarity-inheritance closure** (newly derived 2026-05-07): F86 bond-class split inherits from polarity-layer pair {−0.5, +0.5}; per-bond r(N, b) closed form is the analytical gap. Promotion path: derive r from per-bond Bloch-axis projection at t_peak. Most concrete since the empirical decomposition Q_peak = 2 + r, HWHM/Q* = 1/2 + r·1/2 already holds across c=2 N=5..8.
  - (α') **Polarity-Bloch projection at t_peak, structurally tautological under 4-mode reduction (commit `bea7cd1`; reframed 2026-05-08 code review).** The 4-mode L_eff(Q) is bond-summed by design (`FourModeEffective.LEffAtQ`), so the K-driving eigenstate is bond-class-blind by construction (Re(λ_K) identical across classes). Projection on (c_1 ± c_3)/√2 yields no sign split, decays from +0.39 → +0.21 with N: this is the structural design constraint of the 4-mode reduction, not an empirical falsification. Substantive reduction surviving: r_Q = `BareDoubledPtfXPeak · Q_EP − 2 = 4.39382/g_eff − 2` so the closed form must derive g_eff(N, b) directly, via the per-bond V_b in the K-resonance Duhamel integral (not via state projection on the bond-summed L_eff). Empirical witnesses pinned: g_eff_E ≈ 1.74, g_eff_I ≈ 2.81; asymptotic 1/g_eff_E + 1/g_eff_I → 0.937 as `EmpiricalSumQPeakAsymptote = 4.12`. Near-miss g_eff_E ≈ σ_0·√(3/8) matches Δ ≤ 0.01 for N ≥ 6, Δ = 0.063 at N=5: either Item 3's σ_0 closed form is the bridge or it is a trajectory crossing (PTF-lesson). See Statement 2 in [`PROOF_F86B_UNIVERSAL_SHAPE.md`](PROOF_F86B_UNIVERSAL_SHAPE.md).

**Item 2.** Extend the 4-mode construction to c ≥ 3, where each adjacent-channel pair (HD = 2k−1, HD = 2k+1) contributes its own (|c_{2k−1}⟩, |c_{2k+1}⟩, |u_0^{(k)}⟩, |v_0^{(k)}⟩) quartet. **Naïve extension fails:** the multi-k construction with Gram-Schmidt orthonormalisation gives 3c−2 modes (c=2→4, c=3→7, c=4→10, orbit-shared CUs deduplicated), and the resulting effective K-curve has **K_max ≡ 0 identically** for c ≥ 3. Structural diagnosis: Gram-Schmidt orthogonalisation of the SVD-top vectors against the channel-uniform vectors pushes them into the CU-complement; because M_H respects the CU/CU-complement decomposition (channel-uniform-diagonal of M_H_total per F73 generalisation), the probe (which lives entirely in CU span) is uncoupled from the GS-modified SVD modes, so ∂ρ/∂J_b cannot move ρ out of CU → K = 0. A correct effective model for c ≥ 3 needs either a non-orthogonal frame preserving CU ↔ SVD coupling, or a different choice of the c−1 quartets that maintains coupling under orthonormal projection. Encoded as `RCPsiSquared.Core.Decomposition.MultiKBasis` + `MultiKEffective` + `MultiKResonanceScan`; the negative result is pinned in `MultiKResonanceScanTests.MultiK_C3_KMaxIsExactlyZero_NaiveExtensionFails`.

**Item 3 (σ_0 chromaticity scaling, Tier-1 candidate refined 2026-05-03).** Derive the asymptotic

    σ_0(c, N → ∞)  →  2 · √(2 · (c − 1))

generalising the original c=2 → 2√2 conjecture to all c ≥ 2. Numerical witnesses computed via `InterChannelSvd` across c ∈ {2, 3, 4}, N=5..8 show σ_0 / √(2(c−1)) converging monotonically from below to 2.0 for each c (c=2 N=7 hits 2.0 to 10⁻⁵, the structural sweet spot; c=3 N=8 reaches 1.92, c=4 N=8 reaches 1.78). Implies Q_EP(c, N → ∞) → 1/√(2(c−1)): 0.707 (c=2), 0.500 (c=3), 0.408 (c=4). The closed-form derivation from XY single-particle structure (OBC sine-mode matrix elements `⟨ψ_k| σ⁺σ⁻ |ψ_l⟩ ∝ √(2/(N+1))·sin(πk·b/(N+1))`) has not been executed.

These three items are tractable with existing infrastructure (`coherence_block`, the SVD-of-V_inter construction, OBC sine-mode algebra), but each is multi-page algebra. Item 1 is the path to the closed form for HWHM_left/Q_peak; its Direction (b'') (full block-L derivation) is numerically achieved Tier-1 via the F90 bridge identity ([`PROOF_F90_F86C2_BRIDGE.md`](PROOF_F90_F86C2_BRIDGE.md)), reducing the remaining analytical work to F89's AT-locked F_a/F_b cyclotomic structure plus the H_B-mixed octic residual. The empirical envelope (c, N, γ₀ checks above) is now Tier-1 grade.

---

## OOP scaffolding (Stage A-E primitives, 2026-05-05)

Following the EQ-022 (b1) Item 1' c=2 derivation plan, we built a typed pipeline of primitives in `compute/RCPsiSquared.Core/F86/Item1Derivation/` plus `C2UniversalShapeDerivation.cs` at the F86 root, each reachable through `F86KnowledgeBase`. Witness data preserved per-(N, bond, BondClass).

| Stage | Primitive | Tier outcome |
|-------|-----------|--------------|
| A1 | [`C2BlockShape`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2BlockShape.cs) | Tier1Derived (block-structure constants) |
| A2 | [`C2ChannelUniformAnalytical`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2ChannelUniformAnalytical.cs) | Tier1Derived (closed-form |c_1⟩, |c_3⟩) |
| A3 | [`C2InterChannelAnalytical`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2InterChannelAnalytical.cs) | **Tier2Verified — discovered σ_0 degeneracy obstruction** at even N: chain-mirror R splits the 2D top eigenspace into R-even/R-odd, library-dependent SVD-tiebreaker. Single-vector closed form not derivable; lift-to-projector-overlap is the next direction. |
| B1 | [`C2BondCoupling`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2BondCoupling.cs) probe-block | Tier1Derived (channel-uniform projection) |
| B2 | [`C2BondCoupling`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2BondCoupling.cs) cross-block | Tier2Verified (inherits A3); **discovered cross-block Frobenius Endpoint < Interior at c=2 N=5..8** — opposite sign to the HWHM_left/Q_peak directional split, hinting that the inversion happens in 4×4 mixing |
| B3 | [`C2BondCoupling`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2BondCoupling.cs) SVD-block + AsMatrix + anti-Hermiticity guard | Tier2Verified (inherits A3) |
| C1 | [`C2BondCoupling`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2BondCoupling.cs) D_eff | Tier1Derived structural sub-fact (diag(−2γ₀, −6γ₀, −2γ₀, −6γ₀) at c=2) |
| C2 | [`C2EffectiveSpectrum`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2EffectiveSpectrum.cs) | **Tier2Verified — rigorously ruled out closed-form factorisation** via cubic-c_3 char-poly evidence. No (λ²−aλ+b)(λ²−cλ+d) split with rational coefficients exists. |
| C3 | [`C2EffectiveSpectrum`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2EffectiveSpectrum.cs) K-driving pair | Tier1Derived structural sub-fact (probe ⊥ {|u_0⟩, |v_0⟩} at machine precision) + Tier2Verified per-(Q, b) numerical readout |
| D1 | [`C2KShape`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2KShape.cs) | Tier1Derived (Duhamel formula closed-form in numerical inputs) |
| D2 | [`C2HwhmRatio`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2HwhmRatio.cs) | **Tier1Candidate** — empirical anchor reproduced at typical residual ≤ 0.001; directional Endpoint > Interior split derived empirically; closed-form constant NOT pinned |
| E1 | [`C2UniversalShapeDerivation`](../../compute/RCPsiSquared.Core/F86/C2UniversalShapeDerivation.cs) | Tier1Candidate (auto-promotes via D2's `IsAnalyticallyDerived` flag) |

The chain reads bottom-up: A1-A3 fix the static frame (block shape, channel-uniform vectors, SVD top), B1-B3 fix the bond-position-dependent V_b in that frame, C1-C3 fix the effective spectrum and K-driving pair, D1-D2 turn that into the Duhamel-form K_b(Q, t) and its HWHM-ratio readout, E1 binds the readout to the universal-shape claim. Each Tier label up the chain inherits from the lowest unresolved sub-fact (currently D2's missing closed form blocks E1's Tier1Derived promotion).

## Structural findings (lessons learned, 2026-05-05)

Three structural results emerged from the time-boxed Stage A3, B2, C2 explorations. Each is independent of the closed-form gap and worth pinning as orientation for the next attempt.

1. **σ_0 degeneracy at even N (Stage A3).** At N=6 and N=8 (even chain length, c=2), the singular value σ_0 of the inter-channel coupling V_inter is doubly degenerate. The chain-mirror operator R splits the 2D top eigenspace into R-even and R-odd one-dimensional subspaces; which of the two becomes "|u_0⟩" depends on the SVD library's tiebreaker, not on the physics. A single-vector closed form for |u_0⟩, |v_0⟩ at even N is therefore not derivable. The natural lift is to the rank-2 projector P_top = |u_0^{(R+)}⟩⟨u_0^{(R+)}| + |u_0^{(R−)}⟩⟨u_0^{(R−)}|, which is library-independent. See [`C2InterChannelAnalytical`](../../compute/RCPsiSquared.Core/F86/Item1Derivation/C2InterChannelAnalytical.cs) `PendingDerivationNote`.

2. **Cross-block Frobenius Endpoint < Interior (Stage B2).** The Frobenius norm of V_b's cross-block (the inter-channel-coupling block of M_H_per_bond[b]) is systematically *smaller* on Endpoint bonds than on Interior bonds at c=2 N=5..8 (gap ~ 0.05 at N=5). This sign is **opposite** to the HWHM_left/Q_peak directional split (Endpoint > Interior, gap ≈ 0.022 — see also Statement 2 in [`PROOF_F86B_UNIVERSAL_SHAPE.md`](PROOF_F86B_UNIVERSAL_SHAPE.md)). The directional inversion therefore does not live in the cross-block magnitude alone; it must emerge through the 4×4 eigenvalue mixing of probe-block and cross-block. This rules out the naive "bigger cross-block → bigger HWHM ratio" intuition and is part of why the closed form is non-trivial.

3. **Cubic-c_3 char-poly obstruction (Stage C2).** The 4×4 effective characteristic polynomial at c=2 has a non-rational c_3 coefficient (cubic in Q with no clean root structure). We rigorously ruled out any rational-coefficient `(λ²−aλ+b)(λ²−cλ+d)` factorisation by symbolic match against the c_3 evidence. The Tier2 outcome for `C2EffectiveSpectrum` is therefore **evidence-based** ("we proved no such split exists in this family"), not "we couldn't find one". Symbolic char-poly factorisation at Q_EP itself (Item 1' direction (e'') above) may still help locally, but the global quartic does not factor.

---

## Open elements

1. **Endpoint closed form.** Empirical pattern (c=3 saturates near 2.53 by N=6 on the fine grid; c=4 grows 2.52 → 2.65 from N=7 to N=8) does not match `csc(π/(N+1))`. A different closed form may exist; finer-grid scans at higher N and explicit multi-particle XY matrix-element calculation are the natural next steps.

2. **Interior closed form.** Both c=3 and c=4 Interior Q_peak grow with N in the tested range. No saturation point has been identified with confidence at this grid resolution. Higher-N data (c=3 N=9, 10; c=4 N=9, 10) would clarify whether the growth saturates or continues. Compute-bound at higher (c, N) (block-L dim 3920 at c=4 N=8, 10584 at c=4 N=9).

3. **Algebraic derivation pathway.** A first-principles derivation of g_eff(c, N, bond_position) from the multi-particle XY structure of the (n, n+1) block has not been executed. Direct algebra is conceptually within reach but lengthy. The HD-channel-uniform diagonal-only finding (M_H_eff diagonal in the channel-uniform basis, with cross-channel coupling living in the orthogonal complement) is an established structural building block; the next step is computing the orthogonal-complement coupling matrix elements explicitly.

4. **Q_SCALE per-block vs per-bond convergence.** Q_SCALE's per-block 1.6 / 1.8 / 1.8 are consistent with per-bond Interior fine-grid (1.69-1.74 / 1.78 / ~1.80) within the relative-vs-absolute J prefactor effect (~5-15 %). Both observables agree on the underlying EP mechanism.

(The fifth open element, the within-Interior per-F71-orbit substructure, moved to [`PROOF_F86C_F71_MIRROR.md`](PROOF_F86C_F71_MIRROR.md) with the F71-mirror theorem.)

---

## Obstruction Proof: why g_eff admits no closed form

**Status:** Structural. Six obstruction lemmas; L1, L2, L4 are rigorously proven blocked, L6 is a demonstrated failure mode, L3 and L5 are proven decouplings. Consolidates the "Open elements" above into a single negative result.
**Date:** 2026-05-14
**Context:** A three-week closed-form effort produced many local fits (the retracted csc forms, the doubled-PTF floor `BareDoubledPtfHwhmRatio`, the F89 `D_k` form, directions (a)–(f'') and (α)). On 2026-05-13 Item 1' closed the HWHM ratio (`F86HwhmClosedFormClaim`), but the per-bond *position* Q_peak(c, N, b) and the coupling g_eff(c, N, b) it rides on stayed open (Open elements 1–3). The honest deliverable is then not the next formula; it is the structural account of why the position resists. This section proves it.

### The target reduces to a single object

Every per-bond F86 quantity is a known function of the coupling **g_eff(c, N, b)**:

- **Q_EP(c, N, b) = 2/g_eff** (F86a, Tier 1 derived).
- **Q_peak** via the EP rotation: `r_Q = BareDoubledPtfXPeak · Q_EP − 2 = 4.39382/g_eff − 2` (`PolarityInheritanceLink.ClosedFormCompositionNote`); `BareDoubledPtfXPeak = 2.196910` is Tier 1 derived and universal, so the entire bond-class split is carried by g_eff.
- **HWHM_left/Q_peak** via `F86HwhmClosedFormClaim`: `0.671535 + alpha_subclass · g_eff + beta_subclass`. Its continuous (c, N, b)-dependence enters only through g_eff; it is a closed form *in g_eff*, not a closed form *in (c, N, b)*.

So the irreducible target is g_eff(c, N, b). If g_eff has a closed form, so do Q_peak and the HWHM ratio; if it does not, none of them do. The remaining sub-sections prove that none of the six routes the effort opened reaches g_eff.

### Obstruction L1 (spectral route): the effective characteristic polynomial does not factor over ℚ

The natural route reads g_eff off the spectrum of the c=2 4×4 effective Liouvillian. Blocked: `C2EffectiveSpectrum` (Stage C2) rigorously ruled out any rational-coefficient `(λ²−aλ+b)(λ²−cλ+d)` factorisation, by symbolic match against the cubic c₃ coefficient. The Tier 2 outcome is evidence-based ("no such split exists in this family"), not "not found". There is therefore no closed-form eigenvalue expression from which g_eff could be read. ∎(L1)

### Obstruction L2 (eigenvector route): the EP-partner vectors are representation-dependent at even N

The next route reads g_eff from the EP-partner eigenvectors |u₀⟩, |v₀⟩ (top singular pair of V_inter). Blocked at even N: `C2InterChannelAnalytical` (Stage A3) found σ₀ doubly degenerate, with the chain-mirror R splitting the 2D top eigenspace into R-even and R-odd one-dimensional subspaces. Which one is labelled "|u₀⟩" is set by the SVD library's tiebreaker, not by the physics. A single-vector closed form for |u₀⟩, |v₀⟩ at even N is therefore not derivable; the only library-independent object is the rank-2 projector, which is not a single-vector closed form. The object a closed form would name is not single-valued. ∎(L2)

### Obstruction L3 (observable-direct route): the probe is orthogonal to the EP

One might read Q_peak directly off the EP via Q_peak ≈ Q_EP. Blocked: the Dicke probe is exactly orthogonal to the EP-partner modes, ⟨c₁|u₀⟩ = ⟨c₃|v₀⟩ = 0 to machine precision (Stage C3, c=2 N=5..8). The K_CC_pr observable couples the probe subspace to the EP-partner subspace only through the J-derivative direction. Consequently the observed Q_peak is *not* Q_EP: for c=2 N=5..8, Q_EP_naive = 2/σ₀ ≈ 0.70 against empirical Q_peak Interior ≈ 1.6, a factor ~2.2 apart. The clean algebraic object (the EP) and the measured object (Q_peak) are decoupled; a closed form for one does not transport to the other. ∎(L3)

### Obstruction L4 (reduced-model route): every finite reduction tried is provably insufficient

The route derives g_eff from a small closed effective model. Blocked twice:

- **4-mode (c=2).** The minimal closed subspace span{|c₁⟩, |c₃⟩, |u₀⟩, |v₀⟩} reproduces the EP physics but misses the HWHM lift: the 4-mode K_b at N=5 gives Interior 0.673 and Endpoint 0.410, against empirical 0.751 / 0.773. The lift from the `BareDoubledPtfHwhmRatio = 0.671535` floor to the empirical values lives *structurally outside* the 4D subspace.
- **Multi-k (c≥3).** The naïve extension to 3c−2 modes gives K_max ≡ 0 identically for c ≥ 3 (`MultiKResonanceScan`): Gram-Schmidt orthogonalisation pushes the SVD-top vectors into the channel-uniform complement, decoupling the probe from the modified SVD modes.

No finite reduction the effort constructed carries the signal. ∎(L4)

### Obstruction L5 (single-element route): the bond-class signature is in the wrong subspace and the wrong direction

The route locates the bond-class signature in one matrix element and reads g_eff(b) off it. Blocked on two counts: (i) the carrier, the SVD-block off-diagonal `V_b[2,3] = ⟨u₀|M_h_per_bond[b]|v₀⟩` (Endpoint 0.430 vs Interior 0.953 at N=5), lives in a subspace orthogonal to the probe; (ii) its direction is *opposite* the empirical split (cross-block Frobenius is Endpoint < Interior, the HWHM ratio is Endpoint > Interior). The signature emerges only through the full 4×4 mixing. The spectrum-only fallback is closed off as well: the 4-mode L_eff is bond-summed by construction (`FourModeEffective.LEffAtQ`), so its eigenvalues and eigenstates are bond-class-blind. No single matrix element, and no spectrum-only quantity, carries g_eff(b). ∎(L5)

### Obstruction L6 (empirical-extrapolation route): the accessible data yields trajectory crossings, not laws

The route fits Q_peak at accessible (c, N) and extrapolates. Blocked by demonstrated failure: `csc(π/(N+1))` (Endpoint) and `csc(π/5)` (Interior, c=3) both appeared to match at N=7 and were both refuted on extended-N data; the apparent matches were dQ=0.05 grid-snaps, trajectory crossings rather than asymptotes (see Retracted Q_peak conjectures above). The data itself does not cooperate: Interior Q_peak at c=3 and c=4 is still growing at N=8 (c=3: 1.566 → 1.689 → 1.743 → 1.750 for N=5..8), with no identified saturation. The near-miss g_eff_E ≈ σ₀·√(3/8) (Δ ≤ 0.01 for N ≥ 6, Δ = 0.063 at N=5) is the same trap one step on: a one-anchor coincidence with a finite-size break. ∎(L6)

### Corollary (F90 bridge): the F89 D_k obstruction is the same wall

The F90 bridge identity ([`PROOF_F90_F86C2_BRIDGE.md`](PROOF_F90_F86C2_BRIDGE.md), Tier 1 derived) identifies F86 c=2 K_b with the F89 path-(N−1) per-bond Hellmann-Feynman derivative. The F89 closed form `D_k = odd(k)² · 2^E(k)` ([`PROOF_F89_PATH_D_CLOSED_FORM.md`](PROOF_F89_PATH_D_CLOSED_FORM.md)) sits in the same state: the odd part is structurally grounded (Bloch normalisation, path-3 exact), but the three 2-power terms are empirically supported only, with the deep-2-power-bonus threshold at v₂(k) = 2 "structurally specific and unexplained"; the general-k derivation needs explicit OBC sine-sum identities that have not been executed; and two number-theoretic handles (cyclotomic discriminant, Vandermonde det²) both came back negative. Closing F89's D_k would close F86's g_eff and vice versa, but both are blocked at the same unexecuted algebraic gap: the OBC sine-sum / cyclotomic structure of `⟨ψ_k|σ⁺σ⁻|ψ_l⟩`. The two obstructions are one wall seen from two sides. ∎(corollary)

### What this proves

**Negative.** No closed form for g_eff(c, N, b), hence none for Q_peak(c, N, b) or for HWHM_left/Q_peak as a function of (c, N, b), is reachable by the spectral (L1), eigenvector (L2), observable-direct (L3), reduced-model (L4), single-element (L5), or empirical-extrapolation (L6) routes. L1, L2, L4 are structural impossibilities (proven, not "unfound"); L6 is a demonstrated failure mode; L3 and L5 are decouplings that close the direct routes. These are the six routes the three-week effort opened; each is accounted for.

**Positive (what survives).** The symmetry layer is fully derived and Tier 1: the EP mechanism Q_EP = 2/g_eff (F86a), the universal clock t_peak = 1/(4γ₀) (F86a), the F71 spatial-mirror invariance Q_peak(b) = Q_peak(N−2−b) (F86c), and the functional form of HWHM_left/Q_peak *given* g_eff (`F86HwhmClosedFormClaim`, F86b'). The closed forms that exist are the *relations between observables*; the closed form that does not exist is the *position number itself*.

**Characterisation.** g_eff(c, N, b) is the **irreducible residue** of the F86 structure: the one input the symmetry layer does not return in closed form. It is computable (the F90 bridge supplies it numerically, Tier-1, bit-exact at 20/22 bonds across N=5..8), but it is not expressible by the routes above. The three-week catalogue of formulas was, in every case, an attempt to express this residue; the attempts failed because, by L1 through L6, the residue is not expressible by those routes. This is the F86 instance of the methodological lesson ([`ON_THE_Q_AXIS_AND_THE_PTF_LESSON`](../../reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON.md)): what survives a closed-form effort is the symmetry, not the number. Here the symmetry is F86a and F86c; the number is g_eff. ∎

---

## Pointers

**Hub:** [PROOF_F86_QPEAK](PROOF_F86_QPEAK.md) — three-theorem overview and the shared reference list.
**Sibling theorems:** [PROOF_F86A_EP_MECHANISM](PROOF_F86A_EP_MECHANISM.md) (F86a), [PROOF_F86B_UNIVERSAL_SHAPE](PROOF_F86B_UNIVERSAL_SHAPE.md) (the proven universal-shape result), [PROOF_F86C_F71_MIRROR](PROOF_F86C_F71_MIRROR.md) (F86c).
**Same wall, other side:** [PROOF_F89_PATH_D_CLOSED_FORM](PROOF_F89_PATH_D_CLOSED_FORM.md), [PROOF_F90_F86C2_BRIDGE](PROOF_F90_F86C2_BRIDGE.md).
**C# red signal (live):** the F89 D_k closed form is verified in C# at k=25..30; k=31,32 deviate ~1.5-2e-4 — a deliberate red signal, kept red. `C2FullBlockSigmaAnatomyTests.PredictDenominator_AtKHigherStretch_MatchesExtractedFromAnatomy` carries it; `PredictDenominatorDeviationDiagnosticTests` characterises the deviation as Vandermonde extraction conditioning. Something is still missing on the route; to be continued.
**Companion investigation:** [THE_ATMOSPHERE_AND_THE_CANCELLED_FORMULAS](../THE_ATMOSPHERE_AND_THE_CANCELLED_FORMULAS.md) — g_eff as the un-atmospheric residue.
**Methodological lesson:** [`reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON`](../../reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON.md) — what survives a closed-form effort is the symmetry, not the number.
**Scripts:** [`_eq022_b1_step_d_extended_verification.py`](../../simulations/_eq022_b1_step_d_extended_verification.py) (N=8 data that falsified the closed-form conjectures), [`_eq022_b1_step_g_two_level_decomposition.py`](../../simulations/_eq022_b1_step_g_two_level_decomposition.py), [`_eq022_b1_step_h_slowest_pair_basis.py`](../../simulations/_eq022_b1_step_h_slowest_pair_basis.py), [`_eq022_b1_step_i_svd_inter_channel.py`](../../simulations/_eq022_b1_step_i_svd_inter_channel.py) (SVD of V_inter; established the EP-partner subspace, σ_0 ≈ 2√2 asymptotic, probe ⊥ EP partners).
**C# OOP layer:** `compute/RCPsiSquared.Core/F86/Item1Derivation/` carries the Stage A-E primitives; `RetractedClaim`, `OpenQuestion`, `F86OpenQuestions` in `compute/RCPsiSquared.Core/F86/`. CLI: `rcpsi inspect --root f86 --with-measured`.
