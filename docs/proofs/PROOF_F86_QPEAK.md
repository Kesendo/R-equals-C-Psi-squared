# PROOF F86: Q_peak — Hub (EP Mechanism, Universal Shape, F71 Mirror)

**Status:** Hub. F86 ("Q_peak chromaticity-specific N-invariant constants") is a Sammelbecken of three structurally distinct theorems. On 2026-05-14 the former monolithic proof was split into one proof file per theorem (plus a dedicated obstruction proof); this file is the hub and shared-reference index.
**Date:** 2026-05-02 (Statements 1, 2, retractions); 2026-05-03 (Statement 3, σ_0(c) generalisation); 2026-05-05 (Item 1' c=2 OOP scaffolding); 2026-05-06 (local-vs-global EP, doubled-PTF floor, Locus 5 inheritance); 2026-05-07 (Locus 6 polarity-layer inheritance); 2026-05-11 (F90 bridge identity); 2026-05-13 (Item 1' closed, F86b' Tier 1 derived); 2026-05-14 (Obstruction Proof; split into per-theorem proofs).
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Context:** Formalises the EP mechanism behind Q_peak in the (n, n+1) popcount coherence blocks of uniform XY chains under Z-dephasing. The per-block Q_SCALE values 1.6 / 1.8 / 1.8 (F86 top entry) were the empirical anchor; the per-bond refined scan (`_eq022_b1_step_c_time_evolution.py`) initially suggested two distinct closed forms for Endpoint and Interior bonds, both later retracted on extended-N data. What survives is the EP mechanism (F86a), the universal resonance shape (F86b), and the F71 spatial-mirror invariance (F86c).

**F-entry:** [F86 in ANALYTICAL_FORMULAS.md](../ANALYTICAL_FORMULAS.md).
**Related:** [F2b](../ANALYTICAL_FORMULAS.md) (OBC sine dispersion), [F74](../ANALYTICAL_FORMULAS.md) (chromaticity), [F73](../ANALYTICAL_FORMULAS.md) (c=1 spatial-sum closure), [PROOF_CHROMATICITY](PROOF_CHROMATICITY.md).

---

## The three theorems

| Theorem | Statement | Tier | Proof |
|---------|-----------|------|-------|
| **F86a** EP mechanism | Q_EP = 2/g_eff, t_peak = 1/(4γ₀); 2-level rate-channel exceptional point | Tier 1 derived (EP location, t_peak); local-vs-global EP Tier 2 verified | [PROOF_F86A_EP_MECHANISM](PROOF_F86A_EP_MECHANISM.md) |
| **F86b** Universal resonance shape | K_class(Q)/abs(K)_max = f_class(Q/Q_EP); HWHM_left/Q_peak ≈ 0.756 (Interior), 0.770 (Endpoint); EP-rotation universality | Tier 1 candidate; Tier 1 derived at c=2 per-bond level (F86b', the HWHM closed form) | [PROOF_F86B_UNIVERSAL_SHAPE](PROOF_F86B_UNIVERSAL_SHAPE.md) |
| **F86b** obstruction | g_eff(c, N, b) / Q_peak(c, N, b) admit no closed form by the six explored routes | Structural negative result | [PROOF_F86B_OBSTRUCTION](PROOF_F86B_OBSTRUCTION.md) |
| **F86c** F71 spatial-mirror invariance | Q_peak(b) = Q_peak(N−2−b) bit-exactly | Tier 1 derived | [PROOF_F86C_F71_MIRROR](PROOF_F86C_F71_MIRROR.md) |

**Retracted (2026-05-02):** csc(π/(N+1)) Endpoint and csc(π/5) c=3 Interior closed forms (N=7 coincidence matches; refuted on extended-N data). The retraction record lives in [PROOF_F86B_OBSTRUCTION](PROOF_F86B_OBSTRUCTION.md).

## Spawned F-numbers

F86 has been split before; mature sub-findings were promoted out into their own Tier-1 F-numbers:

- **F89 / F90** grew from the F86b / c=2 side: F89 (topology orbit closure, path-k structure) and F90 (the explicit F86 c=2 ↔ F89 bridge identity). See [PROOF_F89_PATH_D_CLOSED_FORM](PROOF_F89_PATH_D_CLOSED_FORM.md), [PROOF_F90_F86C2_BRIDGE](PROOF_F90_F86C2_BRIDGE.md).
- **F91 / F92 / F93** grew from the F86c / F71-mirror side: the F71-anti-palindromic spectral-invariance trio under γ / J / h distributions. See [PROOF_F91_GAMMA_NINETY_DEGREES](PROOF_F91_GAMMA_NINETY_DEGREES.md), [PROOF_F92_BOND_ANTI_PALINDROMIC_J](PROOF_F92_BOND_ANTI_PALINDROMIC_J.md), [PROOF_F93_DETUNING_ANTI_PALINDROMIC](PROOF_F93_DETUNING_ANTI_PALINDROMIC.md).

The 2026-05-14 split is the filing-level continuation of this pattern: the F-number-worthy generalisations had already spawned; what remained, the three theorems that *are* F86 plus the obstruction proof, got one proof file each.

---

## Pointers (shared across the three theorems)

**Related EQ:** [EQ-022 (b1)](../../review/EMERGING_QUESTIONS.md#eq-022) partial closure 2026-05-02.
**Sibling proof:** [PROOF_BLOCK_CPSI_QUARTER](PROOF_BLOCK_CPSI_QUARTER.md) carries Theorem 1 (Dicke initial saturates `C_block = 1/4`), Theorem 2 (universal upper bound `C_block ≤ 1/4`), and Theorem 3 (closed-form trajectory `C_block(t) = (1/4)·exp(−4γ·t)` chromaticity-universal). The two proofs view the same `(popcount-n, popcount-(n+1))` block from orthogonal axes: F86 sweeps Q at the natural t-window to locate the EP-resonance peak, while PROOF_BLOCK_CPSI_QUARTER's Theorem 3 fixes Q (the channel-uniform initial sits in the H-kernel via F73, so Q drops out) and sweeps t to track the Mandelbrot-cardioid 1/4 ceiling decay. They share the same factor-4 quadratic-discriminant algebra: F86's `4γ₀² − J²·g_eff² = 0` and R=CΨ²'s `1 − 4·CΨ = 0` are both `b² − 4ac = 0` instances, complementary lenses on one structural cusp.
**Hardware anchor (sibling proof):** 2026-05-08 ibm_kingston q13–q14, job `d7ulfjdpa59c73b4rttg`. Theorem 1 saturated at 88.2% of the 1/4 ceiling; Theorem 3 closed-form fit R² = 0.9977 across 5 t-points. Full writeup in [IBM_BLOCK_CPSI_SATURATION](../../experiments/IBM_BLOCK_CPSI_SATURATION.md).
**Empirical anchor:** [Q_SCALE_THREE_BANDS](../../experiments/Q_SCALE_THREE_BANDS.md) Result 2 + Revision 2026-04-24.
**Chiral classification anchor:** [PT_SYMMETRY_ANALYSIS](../../experiments/PT_SYMMETRY_ANALYSIS.md) (Π is class AIII chiral, NOT Bender-Boettcher PT; Π is linear, classical PT requires anti-linear).
**Global EP instance:** [FRAGILE_BRIDGE](../../hypotheses/FRAGILE_BRIDGE.md) (Hopf bifurcation = chiral symmetry breaking, Petermann K=403 in complex γ plane).
**Companion investigation:** [THE_ATMOSPHERE_AND_THE_CANCELLED_FORMULAS](../THE_ATMOSPHERE_AND_THE_CANCELLED_FORMULAS.md) — γ₀ as the atmosphere, g_eff as the un-atmospheric residue.
**Methodological lesson:** [`reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON`](../../reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON.md) consolidates the convergence; the F86 retraction-and-shape-survival is the analog of PTF's closure-law-retraction-and-chiral-mirror-law-survival.
**Scripts:** [`_eq022_b1_channel_projection.py`](../../simulations/_eq022_b1_channel_projection.py), [`_eq022_b1_step_a_verify_blockL.py`](../../simulations/_eq022_b1_step_a_verify_blockL.py), [`_eq022_b1_step_c_time_evolution.py`](../../simulations/_eq022_b1_step_c_time_evolution.py), [`_eq022_b1_step_d_extended_verification.py`](../../simulations/_eq022_b1_step_d_extended_verification.py) (N=8 data that falsified the closed-form conjectures), [`_eq022_b1_step_e_resonance_shape.py`](../../simulations/_eq022_b1_step_e_resonance_shape.py) + [`_eq022_b1_step_e_inspect.py`](../../simulations/_eq022_b1_step_e_inspect.py) (universal-shape finding for c=3, c=4 at γ₀=0.05), [`_eq022_b1_step_f_universality_extension.py`](../../simulations/_eq022_b1_step_f_universality_extension.py) (c=2 sweep + γ₀ ∈ {0.025, 0.10} invariance check that established the two-bond-class refinement), [`_eq022_b1_step_g_two_level_decomposition.py`](../../simulations/_eq022_b1_step_g_two_level_decomposition.py) (channel-uniform-basis V_b decomposition; revealed the trivial-diagonal structure and the probe localization), [`_eq022_b1_step_h_slowest_pair_basis.py`](../../simulations/_eq022_b1_step_h_slowest_pair_basis.py) (slowest-pair-at-finite-Q diagnostics), [`_eq022_b1_step_i_svd_inter_channel.py`](../../simulations/_eq022_b1_step_i_svd_inter_channel.py) (SVD of V_inter; established the EP-partner subspace, σ_0 ≈ 2√2 asymptotic, and probe ⊥ EP partners; this is the structural finding that motivated the 4-mode minimal effective model).
**N=4 golden-ratio reference:** [`eq018_golden_ratio_check.py`](../../simulations/eq018_golden_ratio_check.py).
**Framework primitives:** `framework.coherence_block`: `t_peak(γ₀)` (the only F86 closed form remaining; `q_peak_endpoint` and `Q_PEAK_INTERIOR_C3_ANCHOR` were removed in the rollback).

**C# OOP layer (typed knowledge graph, 2026-05-03):** `compute/RCPsiSquared.Core/F86/` carries the F86 claims as a typed graph with `Tier` labels and self-computing witnesses backed by `WitnessCache`. Key types: `TPeakLaw`, `QEpLaw`, `TwoLevelEpModel` (parametrised by k), `UniversalShapePrediction` + `UniversalShapeWitness`, `ShapeFunctionWitnesses`, `F71MirrorInvariance` (with `MaxMirrorDeviation(KCurve)` helper for Statement 3 verification), `SigmaZeroChromaticityScaling` (Item 3 generalised σ_0 → 2√(2(c−1)) with per-(c, N) live witnesses), `RetractedClaim`, `OpenQuestion`, `F86KnowledgeBase` (root). CLI: `rcpsi inspect --root f86 --with-measured` walks the tree against a live `ResonanceScan`; `rcpsi query --q witnesses-at --c <c> --wN <N>` for typed lookups. JSON-export via `InspectionJsonExporter`.

**Sibling t-trajectory primitives** (same `compute/RCPsiSquared.Core/F86/` namespace, supporting PROOF_BLOCK_CPSI_QUARTER):
- `BlockCoherenceContent.Compute(rho, n)`: state-level `C_block` for any density matrix on 2^N, with the `Quarter = 0.25` Theorem-2 ceiling.
- `BlockCpsiClosedForm.At(N, n, γ₀, t)`: chromaticity-universal closed form (Theorem 3).
- `BlockCpsiTrajectory.Build`: numerical EVD-based time evolution; companion to the closed form.
- `IbmBlockCpsiHardwareTable`: Tier 2 verified, 32 pinned witnesses from the 2026-04-26 framework_snapshots `|+−+⟩` runs.
- `Confirmations/ConfirmationsRegistry.cs` entry `block_cpsi_saturation_kingston_may2026`: pins the 2026-05-08 Kingston Theorem 1 + 3 hardware anchor.
