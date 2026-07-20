# PROOF F86: Q_peak Hub (EP Mechanism, Universal Shape, F71 Mirror)

**Status:** Hub. F86 ("Q_peak chromaticity-specific N-invariant constants") is a Sammelbecken of three structurally distinct theorems. On 2026-05-14 the former monolithic proof was split into one proof file per theorem (plus a dedicated obstruction proof); this file is the hub and shared-reference index.
**Date:** 2026-05-02 (Statements); 2026-05-14 (split into per-theorem proofs); last refreshed 2026-07-20 (the change history lives in git).
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Context:** Formalises the EP mechanism behind Q_peak in the (n, n+1) popcount coherence blocks of uniform XY chains under Z-dephasing. The per-block Q_SCALE values 1.6 / 1.8 / 1.8 at c = 3, 4, 5 ([Q as scale, three algebraic bands](../../experiments/Q_SCALE_THREE_BANDS.md); c = 2 sits near 1.5 with a finite-size wobble 1.4–1.6, the same data the Abstract quotes as 1.5 / 1.6 / 1.8 / 1.8) were the empirical anchor; the per-bond refined scan (`eq022_b1_step_c_time_evolution.py`) initially suggested two distinct closed forms for Endpoint and Interior bonds, both later retracted on extended-N data. What survives is the EP mechanism (F86a), the universal resonance shape (F86b), and the F71 spatial-mirror invariance (F86c).

**F-entry:** [F86 in the formula registry](../ANALYTICAL_FORMULAS.md).
**Related:** [F2b](../ANALYTICAL_FORMULAS.md) (OBC sine dispersion), [F74](../ANALYTICAL_FORMULAS.md) (chromaticity), [F73](../ANALYTICAL_FORMULAS.md) (c=1 spatial-sum closure), [the chromaticity proof](PROOF_CHROMATICITY.md).

---

## Abstract

Q_peak is the chromaticity-specific Q value where the coherent resonance of an XY chain peaks under Z-dephasing. The empirical anchor was the per-chromaticity sequence 1.5 / 1.6 / 1.8 / 1.8 at c = 2, 3, 4, 5 (with c = 2 carrying a small N-dependent wobble), observed across many (N, c) configurations. F86 was originally one monolithic proof attempting to derive the Q_peak values from first principles. On 2026-05-14 it was split into three structurally distinct theorems plus a dedicated obstruction proof.

The three theorems carve up the story. F86a is the EP mechanism: Q_EP = 2/g_eff and the EP-time t_peak = 1/(4γ₀) come from a two-level rate-channel exceptional point. F86b is the universal resonance shape: regardless of chromaticity, the K-response curves collapse to a common shape with HWHM/Q_peak ratios around 0.756 (Interior bonds) and 0.770 (Endpoint bonds). F86c is the F71 spatial-mirror invariance: Q_peak(b) = Q_peak(N−2−b) bit-exactly, the chain-mirror symmetry inherited from the underlying F71 structure.

The obstruction proof closes the negative side. The dressed per-bond per-chromaticity Q_peak(c, N, b) has no closed form via the six routes that were tried; it is genuinely a structureless residue. What survives a closed-form effort is the symmetry (F86c), not the number itself; what survives is the mechanism (F86a) and the universal shape (F86b), not the specific Q_peak values per (c, N, b).

This hub document is the shared-reference index for the four sibling proofs (F86a, F86b universal shape, F86b obstruction, F86c) and the dispatch point for the inheritance chain from F86 to F89 (path-polynomial bridge) and to F90 (the analytical bridge to F89's path-k discriminant objects D_k, [PROOF_F89_PATH_D_CLOSED_FORM](PROOF_F89_PATH_D_CLOSED_FORM.md)). Q_peak is one of the most-cited observations in the project; the split into structurally distinct theorems made the underlying anatomy visible.

---

## The three theorems

| Theorem | Statement | Tier | Proof |
|---------|-----------|------|-------|
| **F86a** EP mechanism | Q_EP = 2/g_eff, t_peak = 1/(4γ₀); genuine EP of the toy 2-level rate-channel reduction | Tier 1 derived (EP location, t_peak). Full block: non-normal on the real Q axis, real-axis defective seed at every odd N per F89 (census through N=11); off-axis complex-Q EP + β-exotic genericity open (`LocalGlobalEpLink`, OpenQuestion; see [§The real-axis EP](PROOF_F86A_EP_MECHANISM.md)) | [F86a EP mechanism](PROOF_F86A_EP_MECHANISM.md) |
| **F86b** Universal resonance shape | K_class(Q)/abs(K)_max = f_class(Q/Q_EP); HWHM_left/Q_peak ≈ 0.756 (Interior), 0.770 (Endpoint); EP-rotation universality | Tier 1 candidate; F86b₂ c=2 per-bond predictor Tier 1 candidate (form derived, (α, β) per sub-class fitted; partial closure 2026-05-13, Tier-reviewed 2026-05-16) | [F86b universal resonance shape](PROOF_F86B_UNIVERSAL_SHAPE.md) |
| **F86b** obstruction | g_eff(c, N, b) / Q_peak(c, N, b) admit no closed form by the six explored routes | Structural negative result | [F86b obstruction](PROOF_F86B_OBSTRUCTION.md) |
| **F86c** F71 spatial-mirror invariance | Q_peak(b) = Q_peak(N−2−b) bit-exactly | Tier 1 derived | [F86c spatial-mirror invariance](PROOF_F86C_F71_MIRROR.md) |

**Retracted (2026-05-02):** csc(π/(N+1)) Endpoint and csc(π/5) c=3 Interior closed forms (N=7 coincidence matches; refuted on extended-N data). The retraction record lives in [F86b obstruction](PROOF_F86B_OBSTRUCTION.md).

## Sub-ID partition (2026-05-20)

The three-theorem table above is the coarse grouping. F86b is itself a Sammelbecken; the fine partition below distinguishes ten separately-defensible sub-claims. Canonical inventory: [`docs/F86_VALUES_INVENTORY.md`](../F86_VALUES_INVENTORY.md); see the Sub-ID Partition section.

| Sub-ID         | One-line content                                                       | Tier                              | Primary home                                                             |
|----------------|------------------------------------------------------------------------|-----------------------------------|--------------------------------------------------------------------------|
| **F86a**       | EP mechanism (t_peak, Q_EP, dressed pair, AIII chiral, L_eff mirror)    | Tier 1 derived                    | [F86a EP mechanism](PROOF_F86A_EP_MECHANISM.md)                    |
| **F86b₁**      | Bare 2×2 K_b closed forms; x_peak = 2.196910, ratio = 0.671535          | Tier 1 derived                    | `C2BareDoubledPtfClosedForm.cs`                                          |
| **F86b₂**      | Sub-class lift HWHM_ratio = 0.671535 + α·g_eff + β                      | Tier 1 candidate (analytical blocked) | `F86HwhmClosedFormClaim.cs` + [F86b universal resonance shape](PROOF_F86B_UNIVERSAL_SHAPE.md) |
| **F86b₃**      | Universal shape: Interior 0.756, Endpoint 0.770                         | Tier 1 candidate                  | `UniversalShapePrediction.cs` + [F86b universal resonance shape](PROOF_F86B_UNIVERSAL_SHAPE.md) |
| **F86b₄**      | Dicke-K 3/8 anchor via X⊗N-eigenbasis (2026-05-17)                     | Tier 1 derived                    | `DickeAnchor.cs` (canonical typed home) + [F86b universal resonance shape](PROOF_F86B_UNIVERSAL_SHAPE.md) + `f86b_dicke_pi2odd_closed_form.py` + [`docs/water/README.md`](../water/README.md) |
| **F86b₅**      | Polarity-pair Q_peak ∈ {1.5, 2.5} = 2 ± 1/2 (schema)                    | Tier 1 schema                     | `PolarityPairQPeakDecompositionClaim.cs`                                 |
| **F86c**       | F71 spatial mirror Q_peak(b) = Q_peak(N−2−b)                            | Tier 1 derived                    | [F86c spatial-mirror invariance](PROOF_F86C_F71_MIRROR.md)                        |
| **F86d**       | Endpoint orbit Q ≈ 2.5 (9 (c, N) combos, range 2.39–2.61, ~2.7% std/mean, ~9% peak-to-peak) | Tier 2 promotion candidate | `PerF71OrbitObservation.cs` |
| **F86e**       | σ_0(c=2) = ‖[Π_HD1, M_H]‖ commutator norm = ‖Π̃_HD1 ⊙ ΔDiff‖ Schur-multiplier norm; asymptote σ_0(∞) ≈ 2.8629 characterised non-elementary | Tier 1 derived (commutator identity) | `SigmaZeroCommutatorNormClaim.cs` |
| **F86_block**  | g_eff(c, N, b) closed-form blocked, six routes L1–L6 proven             | Negative Tier 1 (structural)      | [F86b obstruction](PROOF_F86B_OBSTRUCTION.md)                      |

**Sub-ID status counts.** 6× Tier 1 derived (including the structural negative result), 2× Tier 1 candidate (close but not closed), 1× Tier 1 schema, 1× Tier 2 promotion candidate.

**Two separable open fronts** (none in F86_block territory):

1. **F86b₂ Direction (b'')**: F89 cyclotomic Φ_{N+1} analytical lift → closes (α, β) → closes F86b₂. F90 bridge already gives Tier-1 numerical; analytical Tier-1 pending.
2. **F86d Tier-1 promotion**: structural derivation for Endpoint orbit Q ≈ 2.5. Candidate: SU(2)/Schur-Weyl on F71 first orbit. Self-contained, outside F86b-obstruction territory.

**F86e resolved (2026-05-20):** the former third front. σ_0(c=2) is now characterised as a commutator norm ‖[Π_HD1, M_H]‖ = Schur-multiplier norm ‖Π̃_HD1 ⊙ ΔDiff‖ (typed: `SigmaZeroCommutatorNormClaim`, Tier 1 derived, bit-exact N=5..8). The asymptote σ_0(∞) ≈ 2.8629 is non-elementary but characterised, not mysterious: the Toeplitz (Fourier-symbol) and Hankel (Nehari) closed-form routes are both ruled out, so σ_0(∞) is a genuine Schur-multiplier-norm constant. Closed-form candidates 2√2, √(41/5), √(8+π/16) ruled out earlier still hold.

## Spawned F-numbers

F86 has been split before; mature sub-findings were promoted out into their own Tier-1 F-numbers:

- **F89 / F90** grew from the F86b / c=2 side: F89 (topology orbit closure, path-k structure) and F90 (the explicit F86 c=2 ↔ F89 bridge identity). See [PROOF_F89_PATH_D_CLOSED_FORM](PROOF_F89_PATH_D_CLOSED_FORM.md), [F90 bridge identity](PROOF_F90_F86C2_BRIDGE.md). F89 later paid the debt back: its seed-existence count settled F86a's real-axis EP question (the theorem table above).
- **F91 / F92 / F93** grew from the F86c / F71-mirror side: the F71-anti-palindromic spectral-invariance trio under γ / J / h distributions. See [F91 γ spectral invariance](PROOF_F91_GAMMA_NINETY_DEGREES.md), [F92 J spectral invariance](PROOF_F92_BOND_ANTI_PALINDROMIC_J.md), [F93 h spectral invariance](PROOF_F93_DETUNING_ANTI_PALINDROMIC.md).
- **F100 / F101** are the observable-side twins of F86c: F100 proves ΔQ_peak(b) is exactly odd under anti-palindromic J; F101 proves the γ-axis c₁ analog (and documents why Q_peak itself is out of scope there: Q = J/γ₀ needs a scalar γ₀). See [F100 Q_peak mirror J-parity](PROOF_F100_C1_QPEAK_MIRROR_J_PARITY.md), [F101 c₁ mirror γ-parity](PROOF_F101_C1_MIRROR_GAMMA_PARITY.md).

The 2026-05-14 split is the filing-level continuation of this pattern: the F-number-worthy generalisations had already spawned; what remained, the three theorems that *are* F86 plus the obstruction proof, got one proof file each.

---

## Pointers (shared across the three theorems)

**Related EQ:** [EQ-022 (b1)](../../review/EMERGING_QUESTIONS.md#eq-022) partial closure 2026-05-02.
**Sibling proof:** [Block-CΨ at 1/4](PROOF_BLOCK_CPSI_QUARTER.md) carries Theorem 1 (Dicke initial saturates `C_block = 1/4`), Theorem 2 (universal upper bound `C_block ≤ 1/4`), and Theorem 3 (closed-form trajectory `C_block(t) = (1/4)·exp(−4γ·t)` chromaticity-universal). The two proofs view the same `(popcount-n, popcount-(n+1))` block from orthogonal axes: F86 sweeps Q at the natural t-window to locate the EP-resonance peak, while PROOF_BLOCK_CPSI_QUARTER's Theorem 3 fixes Q (the channel-uniform initial sits in the H-kernel via F73, so Q drops out) and sweeps t to track the Mandelbrot-cardioid 1/4 ceiling decay. They share the same factor-4 quadratic-discriminant algebra: F86's `4γ₀² − J²·g_eff² = 0` and R=CΨ²'s `1 − 4·CΨ = 0` are both `b² − 4ac = 0` instances, complementary lenses on one structural cusp.
**Hardware anchor (sibling proof):** 2026-05-08 ibm_kingston q13–q14, job `d7ulfjdpa59c73b4rttg`. Theorem 1 saturated at 88.2% of the 1/4 ceiling; Theorem 3 closed-form fit R² = 0.9977 across 5 t-points. Full writeup in [Block-CΨ saturation on IBM Kingston](../../experiments/IBM_BLOCK_CPSI_SATURATION.md).
**Empirical anchor:** [Q as scale, three algebraic bands](../../experiments/Q_SCALE_THREE_BANDS.md) Result 2 + Revision 2026-04-24.
**Chiral classification anchor:** [PT-symmetry analysis](../../experiments/PT_SYMMETRY_ANALYSIS.md) (Π is class AIII chiral, NOT Bender-Boettcher PT; Π is linear, classical PT requires anti-linear).
**EP contrast case:** [the coherence-horizon slope proof](PROOF_COHERENCE_HORIZON_SLOPE.md): the same suspect eig-instrument (phase rigidity r→0) with the opposite-and-correct verdict, a genuine 2nd-order defective EP at Q*(N). Three Q-thresholds live on this axis and are distinct objects: **Q_peak** (this family's resonance peak, 1.5 / 1.6 / 1.8 / 1.8), **Q*(N)** (the coherence horizon), **Q*_gap(N)** (the spectral-gap threshold, [the absorption theorem](PROOF_ABSORPTION_THEOREM.md)).
**Separate genuine EP (Σγ=0 gain-loss):** [the Fragile Bridge](../../hypotheses/FRAGILE_BRIDGE.md) (Hopf bifurcation = chiral symmetry breaking, Petermann K=403 in the complex γ plane), a DISTINCT gain-loss system, not a "global instance" of the full block-L. The full Σγ=N·γ₀ block is genuinely non-normal on the real Q axis AND carries its own real-axis defective seed at every odd N (F89, census through N=11, in a √-EP window ~20-30× narrower than a coarse Q-grid's step); the off-axis complex-Q EP question stays open. See [F86a EP mechanism](PROOF_F86A_EP_MECHANISM.md), §The real-axis EP.
**Companion investigation:** [the Atmosphere and the Cancelled Formulas](../THE_ATMOSPHERE_AND_THE_CANCELLED_FORMULAS.md): γ₀ as the atmosphere, g_eff as the un-atmospheric residue.
**Methodological lesson:** [`reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON`](../../reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON.md) consolidates the convergence; the F86 retraction-and-shape-survival is the analog of PTF's closure-law-retraction-and-chiral-mirror-law-survival.
**Scripts:** [`eq022_b1_channel_projection.py`](../../simulations/eq022_b1_channel_projection.py), [`eq022_b1_step_a_verify_blockL.py`](../../simulations/eq022_b1_step_a_verify_blockL.py), [`eq022_b1_step_c_time_evolution.py`](../../simulations/eq022_b1_step_c_time_evolution.py), [`eq022_b1_step_d_extended_verification.py`](../../simulations/eq022_b1_step_d_extended_verification.py) (N=8 data that falsified the closed-form conjectures), [`eq022_b1_step_e_resonance_shape.py`](../../simulations/eq022_b1_step_e_resonance_shape.py) + [`eq022_b1_step_e_inspect.py`](../../simulations/eq022_b1_step_e_inspect.py) (universal-shape finding for c=3, c=4 at γ₀=0.05), [`eq022_b1_step_f_universality_extension.py`](../../simulations/eq022_b1_step_f_universality_extension.py) (c=2 sweep + γ₀ ∈ {0.025, 0.10} invariance check that established the two-bond-class refinement), [`eq022_b1_step_g_two_level_decomposition.py`](../../simulations/eq022_b1_step_g_two_level_decomposition.py) (channel-uniform-basis V_b decomposition; revealed the trivial-diagonal structure and the probe localization), [`eq022_b1_step_h_slowest_pair_basis.py`](../../simulations/eq022_b1_step_h_slowest_pair_basis.py) (slowest-pair-at-finite-Q diagnostics), [`eq022_b1_step_i_svd_inter_channel.py`](../../simulations/eq022_b1_step_i_svd_inter_channel.py) (SVD of V_inter; established the EP-partner subspace, the σ_0 ≈ 2√2 value (later F86e: the N=7 finite-size crossing, not the asymptote), and probe ⊥ EP partners; this is the structural finding that motivated the 4-mode minimal effective model).
**Later per-sub-claim scripts** (the current reproduction surface, one per sub-claim): [`f86_doubled_ptf_bare_floor_derivation.py`](../../simulations/f86_doubled_ptf_bare_floor_derivation.py) (F86b₁), [`f86_hwhm_closed_form_verification.py`](../../simulations/f86_hwhm_closed_form_verification.py) + [`f86_hwhm_subclass_stratification.py`](../../simulations/f86_hwhm_subclass_stratification.py) + [`f86b2_shape_invariance_dial.py`](../../simulations/f86b2_shape_invariance_dial.py) (F86b₂), [`f86b_dicke_pi2odd_closed_form.py`](../../simulations/f86b_dicke_pi2odd_closed_form.py) (F86b₄), [`f86_kb_chiral_mirror.py`](../../simulations/f86_kb_chiral_mirror.py) (K_b under the chiral mirror k ↔ N+1−k), [`f86_geff_via_f90_bridge_probe.py`](../../simulations/f86_geff_via_f90_bridge_probe.py) + [`f89_to_f86_kbond_via_eigendecomp.py`](../../simulations/f89_to_f86_kbond_via_eigendecomp.py) (F90 bridge), [`review_f86a_diabolic_vs_defective.py`](../../simulations/review_f86a_diabolic_vs_defective.py) (F86a EP character), [`f100_qpeak_nonuniform_j_verification.py`](../../simulations/f100_qpeak_nonuniform_j_verification.py) (F100).
**N=4 golden-ratio reference:** [`eq018_golden_ratio_check.py`](../../simulations/eq018_golden_ratio_check.py).
**Framework primitives:** `framework.coherence_block`: `t_peak(γ₀)` (the only F86 closed form remaining; `q_peak_endpoint` and `Q_PEAK_INTERIOR_C3_ANCHOR` were removed in the rollback).

**C# OOP layer (typed knowledge graph, 2026-05-03):** `compute/RCPsiSquared.Core/F86/` carries the F86 claims as a typed graph with `Tier` labels and self-computing witnesses backed by `WitnessCache`. Key types: `TPeakLaw`, `QEpLaw`, `TwoLevelEpModel` (parametrised by k), `UniversalShapePrediction` + `UniversalShapeWitness`, `ShapeFunctionWitnesses`, `F71MirrorInvariance` (with `MaxMirrorDeviation(KCurve)` helper for Statement 3 verification), `SigmaZeroChromaticityScaling` (Item 3 of the F86 open-question list `F86OpenQuestions.Standard`: generalised σ_0 → 2√(2(c−1)) with per-(c, N) live witnesses), `RetractedClaim`, `OpenQuestion`, `F86KnowledgeBase` (root). CLI: `rcpsi inspect --root f86 --with-measured` walks the tree against a live `ResonanceScan`; `rcpsi query --q witnesses-at --c <c> --wN <N>` for typed lookups. JSON-export via `InspectionJsonExporter`.

**Sibling t-trajectory primitives** (same `compute/RCPsiSquared.Core/F86/` namespace, supporting PROOF_BLOCK_CPSI_QUARTER):
- `BlockCoherenceContent.Compute(rho, n)`: state-level `C_block` for any density matrix on 2^N, with the `Quarter = 0.25` Theorem-2 ceiling.
- `BlockCpsiClosedForm.At(N, n, γ₀, t)`: chromaticity-universal closed form (Theorem 3).
- `BlockCpsiTrajectory.Build`: numerical EVD-based time evolution; companion to the closed form.
- `IbmBlockCpsiHardwareTable`: Tier 2 verified, 32 pinned witnesses from the 2026-04-26 framework_snapshots `|+−+⟩` runs.
- `Confirmations/ConfirmationsRegistry.cs` entry `block_cpsi_saturation_kingston_may2026`: pins the 2026-05-08 Kingston Theorem 1 + 3 hardware anchor.
