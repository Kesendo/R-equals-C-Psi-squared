# F86 Values Inventory

Structured inventory separating **structural anchors** (named, N-invariant or saturated, Tier 1 or Tier 2 with saturation evidence) from **per-(c, N, orbit) empirical observations** (frozen measurements, may revise under refined methodology).

The anchor list is what `Q_REGIME_ANCHORS.md` should target; the empirical observations are useful diagnostics but should NOT be promoted to anchor status without further structural derivation.

---

## Sub-ID Partition (2026-05-20)

F86 was conceptually a Sammelbecken of three theorems (F86a / F86b / F86c). The fine partition below distinguishes ten separately-defensible sub-claims, each with its own tier status, derivation status, and home file. The earlier three-way structure remains valid as a coarse grouping; the ten-way structure shows where the work actually lives.

| Sub-ID         | Content                                                                | Tier                              | Home                                                                                                                            |
|----------------|------------------------------------------------------------------------|-----------------------------------|---------------------------------------------------------------------------------------------------------------------------------|
| **F86a**       | EP mechanism: t_peak, Q_EP, dressed pair, AIII chiral, L_eff mirror axis | Tier1Derived                      | [`TPeakLaw.cs`](../compute/RCPsiSquared.Core/F86/TPeakLaw.cs) + [`QEpLaw.cs`](../compute/RCPsiSquared.Core/F86/QEpLaw.cs) + [`LEffMirrorAxisClaim.cs`](../compute/RCPsiSquared.Core/F86/LEffMirrorAxisClaim.cs) + [`PROOF_F86A_EP_MECHANISM.md`](proofs/PROOF_F86A_EP_MECHANISM.md) |
| **F86b₁**      | Bare 2×2 K_b closed forms (post-EP, pre-EP, EP-limit); x_peak = 2.196910, ratio = 0.671535 | Tier1Derived                      | [`C2BareDoubledPtfClosedForm.cs`](../compute/RCPsiSquared.Core/F86/Item1Derivation/C2BareDoubledPtfClosedForm.cs)                |
| **F86b₂**      | Sub-class lift: HWHM_ratio = 0.671535 + α_subclass·g_eff + β_subclass  | Tier1Candidate (analytical blocked) | [`F86HwhmClosedFormClaim.cs`](../compute/RCPsiSquared.Core/F86/Item1Derivation/F86HwhmClosedFormClaim.cs) + [`PROOF_F86B_UNIVERSAL_SHAPE.md`](proofs/PROOF_F86B_UNIVERSAL_SHAPE.md) |
| **F86b₃**      | Universal shape: Interior HWHM/Q = 0.756 ± 0.005, Endpoint = 0.770 ± 0.005 | Tier1Candidate                    | [`UniversalShapePrediction.cs`](../compute/RCPsiSquared.Core/F86/UniversalShapePrediction.cs) + [`UniversalShapeWitness.cs`](../compute/RCPsiSquared.Core/F86/UniversalShapeWitness.cs) |
| **F86b₄**      | Dicke-K 3/8 anchor via X⊗N-eigenbasis decomposition                    | Tier1Derived (2026-05-17)         | [`DickeAnchor.cs`](../compute/RCPsiSquared.Core/Symmetry/DickeAnchor.cs) (canonical typed home, `KIntermediate => 3.0/8.0`) + [`PROOF_F86B_UNIVERSAL_SHAPE.md`](proofs/PROOF_F86B_UNIVERSAL_SHAPE.md) (derivation) + [`f86b_dicke_pi2odd_closed_form.py`](../simulations/f86b_dicke_pi2odd_closed_form.py) (verification) + [`docs/water/README.md`](water/README.md) (downstream reading) |
| **F86b₅**      | Polarity-pair Q_peak ∈ {1.5, 2.5} = 2 ± 1/2 (schema-derived)           | Tier1Schema (value Tier2Verified) | [`PolarityPairQPeakDecompositionClaim.cs`](../compute/RCPsiSquared.Core/F86/PolarityPairQPeakDecompositionClaim.cs)              |
| **F86c**       | F71 spatial mirror: Q_peak(b) = Q_peak(N−2−b) bit-exact                | Tier1Derived                      | [`F71MirrorInvariance.cs`](../compute/RCPsiSquared.Core/F86/F71MirrorInvariance.cs) + [`PROOF_F86C_F71_MIRROR.md`](proofs/PROOF_F86C_F71_MIRROR.md) |
| **F86d**       | Endpoint orbit Q ≈ 2.5 (9 (c, N) combinations, range 2.39–2.61, ~2.7% std/mean around 2.547; ~9% peak-to-peak) | Tier2Empirical (promotion candidate) | [`PerF71OrbitObservation.cs`](../compute/RCPsiSquared.Core/F86/PerF71OrbitObservation.cs) |
| **F86e**       | σ_0(c=2) = ‖[Π_HD1, M_H]‖ (commutator norm) = ‖Π̃_HD1 ⊙ ΔDiff‖ (Schur-multiplier norm); asymptote σ_0(∞) ≈ 2.8629 characterised non-elementary (Toeplitz + Hankel symbol routes ruled out) | Tier1Derived (commutator identity) | [`SigmaZeroCommutatorNormClaim.cs`](../compute/RCPsiSquared.Core/F86/SigmaZeroCommutatorNormClaim.cs) + [`SigmaZeroChromaticityScaling.cs`](../compute/RCPsiSquared.Core/F86/SigmaZeroChromaticityScaling.cs) + [`F86OpenQuestions.cs`](../compute/RCPsiSquared.Core/F86/F86OpenQuestions.cs) (asymptote-value home) |
| **F86_block**  | Obstruction proof: g_eff(c, N, b) has no closed form (6 routes L1–L6 blocked) | Negative Tier1 (structural result) | [`PROOF_F86B_OBSTRUCTION.md`](proofs/PROOF_F86B_OBSTRUCTION.md)                                                                  |

### Sub-ID status counts

- **6× settled at Tier 1:** F86a, F86b₁, F86b₄, F86c, F86e (commutator identity) are Tier1Derived positive results; F86_block is a Tier 1 negative/structural-obstruction result (the closed form does not exist)
- **2× Tier1Candidate (close, not closed):** F86b₂, F86b₃
- **1× Tier1Schema:** F86b₅
- **1× Tier2 promotion candidate:** F86d

### Two separable open fronts

The remaining work lives in two structurally distinct directions; neither overlaps with the F86_block obstruction:

1. **F86b₂ Direction (b'')**: F89 cyclotomic Φ_{N+1} analytical lift → closes (α, β) → closes F86b₂. F90 bridge already gives Tier-1 numerical; analytical Tier-1 pending.
2. **F86d Tier-1 promotion**: structural derivation for Endpoint orbit Q ≈ 2.5. Candidate: SU(2)/Schur-Weyl on F71 first orbit.

**F86e resolved (2026-05-20):** the former third front. σ_0(c=2) is now characterised as a commutator norm ‖[Π_HD1, M_H]‖ = Schur-multiplier norm ‖Π̃_HD1 ⊙ ΔDiff‖ (typed: `SigmaZeroCommutatorNormClaim`, Tier1Derived, bit-exact N=5..8). The asymptote σ_0(∞) ≈ 2.8629 is non-elementary but characterised, not mysterious: the Toeplitz (Fourier-symbol) and Hankel (Nehari) closed-form routes are both ruled out, so σ_0(∞) is a genuine Schur-multiplier-norm constant. Closed-form candidates 2√2, √(41/5), √(8+π/16) ruled out earlier still hold.

### Spawned out (already promoted to own F-numbers)

- **F89** (Path polynomial pipeline, Chebyshev + BigRational, D_k bit-exact): Tier1Derived since 2026-05-15
- **F90** (F86 c=2 K_b ↔ F89 path-(N−1) Hellmann-Feynman bridge): Tier1 numerical, analytical pending
- **F91 / F92 / F93** (F71-anti-palindromic spectral invariance on γ / J / h)

---

## Seen again (2026-06-11): the fresh-eyes wave and the hardware anchor

Three things landed after this inventory's last update (2026-06-01) that bear on its rows; recorded here dated, anchors unchanged unless said so.

1. **Q_EP got its hardware anchor.** The EP-onset run (ibm_kingston q13-q14-q15, 2026-05-31, registered as `ibm_ep_onset_may2026` in both Confirmations registries): the single-excitation revival sits pinned at the 1/N = 1/3 equipartition floor for Q ≤ 1.5 and lifts off as Q crosses 1.5 → 2.5 (0.34 → 0.49 → 0.56 → 0.70), the rotation born at the F86a exceptional point, watched on a real chip with populations only. This puts a **measured** row under the Q_EP entries of the anchor table below, which until now carried only the idealized 2.0 and the c=2 wobble. Writeup: [THE_FLOW_BETWEEN_TWO_SINGULARITIES](../experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md).
2. **g_eff's irreducible content halved.** The K_b mode mirror (2026-06-10): the antiunitary T = Σ₁∘conj gives K_b(b; mode k) = K_b(b; mode N+1−k) pointwise in Q and t, so **g_eff(c, N, b)'s irreducible content lives in the half-band of {k, N+1−k} orbits** ([PROOF_F86B_OBSTRUCTION](proofs/PROOF_F86B_OBSTRUCTION.md) "one more surviving symmetry", [PROOF_PTF_CHIRAL_MIRROR_RATE_LAW](proofs/PROOF_PTF_CHIRAL_MIRROR_RATE_LAW.md) §8). A structural bonus from the same run: on c = 1 blocks K_b ≡ 0 identically (the block generator is e^{−2γt} × unitary), the structural reason the F86 resonance needs c ≥ 2. This does not revise any anchor value, but it halves the target of the F86_block obstruction and of Direction (b'') below: whatever closed form g_eff has or refuses to have, it is a function on mode *pairs*.
3. **The clock and the seam, pinned.** The EP clock's peak angle is γ₀-invariant and provably not the diagonal: θ_peak = 44.3646° (45° would need x = √5; x_peak = 2.196910 of F86b₁ gives the 0.6354° shortfall; `ExceptionalPointClockTests`). And the Absorption Theorem extensions (2026-06-10) recentre the dissipative seam as L_D = γ(Q − N·I) with the per-eigenmode Rayleigh split and the orthogonal-projector form ([PROOF_ABSORPTION_THEOREM](proofs/PROOF_ABSORPTION_THEOREM.md) Extensions), the structural floor under every Re-rate in this inventory. The slow-light surface gained its basis-free primitive (`SlowLightDistribution`, [compute/RCPsiSquared.Diagnostics/Ptf/](../compute/RCPsiSquared.Diagnostics/Ptf/SlowLightDistribution.cs)): light_l = Tr(Π_V·Δ_l)/dim V, the s\*-surface read without hand-picked endpoint pairs.

---

## Anchors (structurally-named, stable across the verified parameter range)

### A1. Universal F86 constants (Tier1Derived)

| Anchor                       | Value              | Source                          |
|------------------------------|--------------------|---------------------------------|
| **Q_EP = 2/g_eff** (formula) | g_eff-dependent    | [`QEpLaw.cs:14-27`](../compute/RCPsiSquared.Core/F86/QEpLaw.cs) |
| **t_peak = 1/(4γ₀)**          | 5 sim-units (γ₀=0.05) | TPeakLaw (via F86KnowledgeBase) |
| **Q_EP at g_eff=1**          | **2.0** (idealized) | direct from QEpLaw at g_eff=1   |

The Q_EP = 2/g_eff law is the structural anchor; the value 2.0 emerges only at the ideal g_eff=1 case.

### A2. Per-chromaticity Q_peak (Tier2Empirical, saturated N-invariant)

| c   | Q_peak | Saturation status                         | Source                       |
|-----|--------|-------------------------------------------|------------------------------|
| 2   | 1.5    | Tier1Candidate, wobble 1.4-1.6 across N=4..9 (finite-size sensitive) | `PerBlockQPeakClaim.cs:63` |
| 3   | **1.6** | Tier2Empirical, saturated across N=5..9   | `PerBlockQPeakClaim.cs:65`   |
| 4   | **1.8** | Tier2Empirical, saturated across N=7..9   | `PerBlockQPeakClaim.cs:66`   |
| 5   | **1.8** | Tier2Empirical (single point N=9, plateau at 1.8 not 2.0) | `PerBlockQPeakClaim.cs:67` |

Saturation pattern (1.5 c=2, 1.6 c=3, 1.8 c=4 and c=5) refuted the previous linear-extrapolation hypothesis Q_peak(c=5)=2.0 (40h N=9 Phase-2 run, 2026-04-24). The plateau at 1.8 is the Tier-2 structural fact.

c=2's "1.5" carries the wobble caveat: it's *not* clean Tier 1 like the c≥3 values. Treated as Tier 1 Candidate.

### A3. Endpoint orbit Q (Tier2Empirical candidate)

| Orbit                  | Q stable value         | (c, N) range tested       | Variation       |
|------------------------|------------------------|---------------------------|-----------------|
| **Endpoint (orbit 0)** | **≈ 2.5**             | (c=2..4, N=5..8); 9 of 12 box cells filled (c=4 only at N=7) | 2.39-2.61 (~2.7% std/mean around 2.547; ~9% peak-to-peak) |

From `PerF71OrbitObservation.cs:23-44`: across 9 different (c, N) combinations, Endpoint orbit Q sits in [2.39, 2.61]. Stability across c and N makes this a candidate anchor not currently in QBasisAnkers.

**Status**: Tier 2 Empirical candidate for promotion. Pending structural derivation.

---

## Empirical observations (not anchor-grade, may revise)

The following values are measurements from the typed C# claims. They are useful for diagnostic comparison but should NOT be treated as structural anchors. Each carries a methodological caveat that may revise the specific number.

### B1. Per-orbit Q_EP values (computed from x_peak)

`C2BlockCpsiQScan.cs:19-20`: Q_EP per orbit = Q_peak / x_peak with x_peak = 2.196910.

| Orbit (N=5 c=2)      | Q_EP value         |
|----------------------|--------------------|
| Endpoint (Orbit 0)   | 1.138              |
| Interior (Orbit 1)   | 0.674              |

**Caveat**: x_peak = 2.196910 is the **bare-doubled-PTF SVD-block-only floor**; full block-L empirical values lift 8-10% above. These per-orbit Q_EP values therefore may shift if the full-block correction is properly included. Treat as derivative observations, not as primary anchors.

### B2. Off-grid escape regime (orbit 1 plateau, central plateau)

From `PerF71OrbitObservation.cs:23-44`:

| (c, N) | Orbit type           | Q-value reached |
|--------|----------------------|-----------------|
| (2, 7) | Orbit 1 plateau      | ≈ 7.24          |
| (2, 8) | Orbit 1 plateau      | ≈ 8.07          |
| (2, 8) | Central orbit plateau | ≈ 16.79         |

**Caveat**: these are NOT EP peaks but broad high-Q plateaus where K stays elevated past Q_EP. Per `project_q_peak_ep_structure` (line 63), "the plateau is multi-mode coherent enhancement, NOT a Petermann-type EP collapse" (Petermann factor probe disconnects the plateau from EP physics). The Q values here mark where K = max within a [0.20, 6.00] (or higher) grid, but they sit on a broad plateau, not a sharp peak.

These are observation-of-multi-mode-plateau values, NOT structural Q-anchors. The methodological reading: these orbits "escape" the F86 EP-mediated anchor structure and live in a multi-mode-coherent regime that has its own physics not yet typed.

### B3. Per-bond-class HWHM ratios (Tier1Candidate empirical means)

`C2HwhmRatio.cs:11,59` (Endpoint/Interior empirical anchors at line 11, BareDoubledPtfHwhmRatio floor at line 59):

| Bond class | Empirical mean | SVD floor | Lift above floor |
|------------|----------------|-----------|------------------|
| Endpoint   | 0.7728         | 0.671535  | +0.101           |
| Interior   | 0.7506         | 0.671535  | +0.079           |

**Caveat**: Tier1Candidate. The floor 0.671535 IS derived (Tier1Derived, bare-doubled-PTF closed form). The lifts +0.08 to +0.10 are empirical observations; the analytical derivation is **open** (PROOF_F90_F86C2_BRIDGE.md Item 1'-followup, F89 AT-locked F_a/F_b structure). The specific empirical mean values (0.7728, 0.7506) may shift under refined Q-grid resolution or under closed-form derivation.

### B4. Per-sub-class fitted (α, β): polyfit, not derived

`F86HwhmClosedFormClaim.cs:51-56`. Form: HWHM_ratio = 0.671535 + α·g_eff + β, all 6 (α, β) pairs **polyfit on N=5..8**, not analytically derived:

| BondSubClass         | α          | β          |
|----------------------|------------|------------|
| Endpoint             | −0.129110  | +0.227413  |
| Flanking             | −0.094978  | +0.193098  |
| Mid                  | +0.056559  | +0.005165  |
| CentralSelfPaired    | +0.057439  |  0.000000  |
| Orbit2Escape         | +0.698446  | −0.086386  |
| CentralEscapeOrbit3  | −0.400854  |  0.000000  |

**Caveat**: All 12 numbers are polyfit on a small N-range. Analytical derivation is open. The specific values are explicitly NOT anchors and should not be treated as structural constants.

---

## Anchor summary (combined with Q_REGIME_ANCHORS)

The anchors merging Q-band edges, per-c Q_peak, Endpoint orbit, and the universal constants:

| Anchor              | Q value | Source                                  | Tier            |
|---------------------|---------|-----------------------------------------|-----------------|
| onset start         | 0.2     | Q-band lower bound                      | Tier2Empirical  |
| onset end           | 0.35    | Q-band upper bound                      | Tier2Empirical  |
| Balance             | 1.0     | J = γ₀                                  | Tier1Derived    |
| peak start          | 1.2     | Q-band lower bound                      | Tier2Empirical  |
| F86 Q_peak c=2      | 1.5     | PerBlockQPeak (wobble 1.4-1.6)          | Tier1Candidate  |
| F86 Q_peak c=3      | 1.6     | PerBlockQPeak (saturated)               | Tier2Empirical  |
| F86 Q_peak c=4, c=5 | 1.8     | PerBlockQPeak (saturated)               | Tier2Empirical  |
| Q_EP at g_eff=1     | 2.0     | QEpLaw idealized                        | Tier1Derived    |
| Q_EP onset (hardware) | ≈ 1.5 | `ibm_ep_onset_may2026`: revival pinned at 1/N for Q ≤ 1.5, liftoff above (2026-05-31, Kingston) | Tier2Verified (hardware) |
| Endpoint orbit Q    | 2.5     | PerF71OrbitObservation (stable, ~2% N-variation) | Tier2Empirical candidate |

Only the anchors above belong in QBasisAnkers and Q_REGIME_ANCHORS. The per-orbit Q_EP, off-grid escapes, HWHM lifts, and fitted (α, β) values are empirical diagnostics that may revise; they live in this inventory and the typed C# claims, not in the anchor map.

## Anchors and cross-refs

- **Post-EP downstream (added 2026-06-01, interpretive, NO anchor change):** the after-life of the
  F86a EP is read in [`THE_FLOW_BETWEEN_TWO_SINGULARITIES`](../experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md)
  (the single-excitation flow from the EP into the 1/N fixed point) and assembled in
  [`THE_VIEW_ONTO_THE_MEMORY`](../reflections/THE_VIEW_ONTO_THE_MEMORY.md) (the depth = light = rate
  axis, the even/odd parity rail, the bilinear ½/¼ currency, the EP read as where the longest-lived
  mode switches rung). These read the post-EP *dynamics*; they do NOT revise the Q_peak / Q_EP /
  HWHM anchors above, and the flow's even-mode crown-switch EP is a distinct object from the c=2
  coherence-block Q_EP this inventory tracks (kept separate on purpose, not merged). The honest
  caution lives in THE_FLOW's "Seen again 2026-05-31" section: the tempting EP coincidences (slowest
  rate = 2γ, ⟨n_XY⟩ = 1 at the EP) are uniform-specific, broken by a non-uniform γ-profile, so they
  are not anchors.
- High-level Q-anchor map: [`Q_REGIME_ANCHORS.md`](Q_REGIME_ANCHORS.md)
- C# F86 KnowledgeBase root: [`compute/RCPsiSquared.Core/F86/`](../compute/RCPsiSquared.Core/F86/)
- Sub-ID partition (this document, top section): F86a / F86b₁..b₅ / F86c / F86d / F86e / F86_block
- Proof hub (sub-ID navigation): [`PROOF_F86_QPEAK.md`](proofs/PROOF_F86_QPEAK.md)
- Per-sub-claim proofs: [`PROOF_F86A_EP_MECHANISM.md`](proofs/PROOF_F86A_EP_MECHANISM.md), [`PROOF_F86B_UNIVERSAL_SHAPE.md`](proofs/PROOF_F86B_UNIVERSAL_SHAPE.md), [`PROOF_F86B_OBSTRUCTION.md`](proofs/PROOF_F86B_OBSTRUCTION.md), [`PROOF_F86C_F71_MIRROR.md`](proofs/PROOF_F86C_F71_MIRROR.md)
- Open derivation gap (HWHM lifts + α/β closed forms): PROOF_F90_F86C2_BRIDGE.md Item 1'-followup; only Direction (b'') survives Item 1' falsification arc (2026-05-11)
- Memory: `project_q_peak_ep_structure`, `project_q_middle_structure`, `project_no_classicalization`, `feedback_f86_needs_agents`; post-EP downstream (2026-05-31): `project_birth_canal_and_light`, `project_q_as_lifetime_of_the_new`
