# F86 Values Inventory

Structured inventory separating **structural anchors** (named, N-invariant or saturated, Tier 1 or Tier 2 with saturation evidence) from **per-(c, N, orbit) empirical observations** (frozen measurements, may revise under refined methodology).

The anchor list is what `Q_REGIME_ANCHORS.md` should target; the empirical observations are useful diagnostics but should NOT be promoted to anchor status without further structural derivation.

---

## Sub-ID Partition (2026-05-20)

F86 was conceptually a Sammelbecken of three theorems (F86a / F86b / F86c). The fine partition below distinguishes ten separately-defensible sub-claims, each with its own tier status, derivation status, and home file. The earlier three-way structure remains valid as a coarse grouping; the ten-way structure shows where the work actually lives.

| Sub-ID         | Content                                                                | Tier                              | Home                                                                                                                            |
|----------------|------------------------------------------------------------------------|-----------------------------------|---------------------------------------------------------------------------------------------------------------------------------|
| **F86a**       | EP mechanism: t_peak, Q_EP, dressed pair, Π-centered P-type mirror axis | Tier1Derived                      | [`TPeakLaw.cs`](../compute/RCPsiSquared.Core/F86/TPeakLaw.cs) + [`QEpLaw.cs`](../compute/RCPsiSquared.Core/F86/QEpLaw.cs) + [`LEffMirrorAxisClaim.cs`](../compute/RCPsiSquared.Core/F86/LEffMirrorAxisClaim.cs) + [`PROOF_F86A_EP_MECHANISM.md`](proofs/PROOF_F86A_EP_MECHANISM.md) |
| **F86b₁**      | Bare 2×2 K_b closed forms (post-EP, pre-EP, EP-limit); envelope peak x_peak = 1.7141131714, HWHM_left/x_peak = 0.7264892771; the fixed-τ = 1 slice peaks at x_peak = 2.196910, ratio 0.671535 | Tier1Derived                      | [`C2BareDoubledPtfClosedForm.cs`](../compute/RCPsiSquared.Core/F86/Item1Derivation/C2BareDoubledPtfClosedForm.cs)                |
| **F86b₂**      | Sub-class lift above the envelope floor 0.7264892771, described as α_subclass·g_eff + β_subclass; not an exact law (Endpoint and Flanking miss a line by 1.28·10⁻³ and 4.2·10⁻⁴) | Tier2Empirical (analytical lift open) | [`F86HwhmClosedFormClaim.cs`](../compute/RCPsiSquared.Core/F86/Item1Derivation/F86HwhmClosedFormClaim.cs) + [`PROOF_F86B_UNIVERSAL_SHAPE.md`](proofs/PROOF_F86B_UNIVERSAL_SHAPE.md) |
| **F86b₃**      | Shape band: Interior HWHM_left/Q_peak 0.7523..0.7637, Endpoint 0.7710..0.7815 (the span of the read witnesses, c = 2, 3 at N = 5..8 and c = 4 at N = 7, 8) | Tier2Empirical (a record of the read span) | [`UniversalShapeBand.cs`](../compute/RCPsiSquared.Core/F86/UniversalShapeBand.cs) + [`UniversalShapePrediction.cs`](../compute/RCPsiSquared.Core/F86/UniversalShapePrediction.cs) + [`UniversalShapeWitness.cs`](../compute/RCPsiSquared.Core/F86/UniversalShapeWitness.cs) |
| **F86b₄**      | Dicke-K 3/8 anchor via X⊗N-eigenbasis decomposition                    | Tier1Derived (2026-05-17)         | [`DickeAnchor.cs`](../compute/RCPsiSquared.Core/Symmetry/DickeAnchor.cs) (canonical typed home, `KIntermediate => 3.0/8.0`) + [`PROOF_F86B_UNIVERSAL_SHAPE.md`](proofs/PROOF_F86B_UNIVERSAL_SHAPE.md) (derivation) + [`f86b_dicke_pi2odd_closed_form.py`](../simulations/f86b_dicke_pi2odd_closed_form.py) (verification) + [`docs/water/README.md`](water/README.md) (downstream reading) |
| **F86b₅**      | Polarity-pair Q_peak ∈ {1.5, 2.5} = 2 ± 1/2 (schema-derived)           | Tier1Schema (value Tier2Verified) | [`PolarityPairQPeakDecompositionClaim.cs`](../compute/RCPsiSquared.Core/F86/PolarityPairQPeakDecompositionClaim.cs)              |
| **F86c**       | F71 spatial mirror: Q_peak(b) = Q_peak(N−2−b) exactly (readings agree to ~10⁻¹⁰, up to 7·10⁻¹⁰ at c = 4, N = 8) | Tier1Derived                      | [`F71MirrorInvariance.cs`](../compute/RCPsiSquared.Core/F86/F71MirrorInvariance.cs) + [`PROOF_F86C_F71_MIRROR.md`](proofs/PROOF_F86C_F71_MIRROR.md) |
| **F86d**       | Endpoint orbit Q ≈ 2.5 (10 (c, N) combinations, range 2.393–2.591, population std/mean 2.2% around 2.547, peak-to-peak 7.8%) | Tier2Empirical (promotion candidate) | [`PerF71OrbitObservation.cs`](../compute/RCPsiSquared.Core/F86/PerF71OrbitObservation.cs) |
| **F86e**       | σ_0(c=2) = ‖[Π_HD1, M_H]‖ (commutator norm) = ‖Π̃_HD1 ⊙ ΔDiff‖ (Schur-multiplier norm); asymptote σ_0(∞) ≈ 2.8629 characterised non-elementary (Toeplitz + Hankel symbol routes ruled out) | Tier1Derived (commutator identity) | [`SigmaZeroCommutatorNormClaim.cs`](../compute/RCPsiSquared.Core/F86/SigmaZeroCommutatorNormClaim.cs) + [`SigmaZeroChromaticityScaling.cs`](../compute/RCPsiSquared.Core/F86/SigmaZeroChromaticityScaling.cs) + [`F86OpenQuestions.cs`](../compute/RCPsiSquared.Core/F86/F86OpenQuestions.cs) (asymptote-value home) |
| **F86_block**  | Obstruction proof: g_eff(c, N, b) has no closed form (6 routes L1–L6 blocked) | Negative Tier1 (structural result) | [`PROOF_F86B_OBSTRUCTION.md`](proofs/PROOF_F86B_OBSTRUCTION.md)                                                                  |

### Sub-ID status counts

- **6× settled at Tier 1:** F86a, F86b₁, F86b₄, F86c, F86e (commutator identity) are Tier1Derived positive results; F86_block is a Tier 1 negative/structural-obstruction result (the closed form does not exist)
- **1× Tier1Schema:** F86b₅
- **2× Tier2Empirical:** F86b₂ (the lift is observed; no exact law holds it), F86b₃ (the band records the read span; it cannot fail)
- **1× Tier2 promotion candidate:** F86d

### Two separable open fronts

The remaining work lives in two structurally distinct directions; neither overlaps with the F86_block obstruction:

1. **F86b₂ Direction (b'')**: derive the lift above the envelope floor from the F89 cyclotomic Φ_{N+1} structure. The anchors do not support a line in g_eff, so the target is the lift itself, not a pair (α, β). The F90 bridge is the numerical identity under it (all 22 bonds of N = 5..8 agree with the F89 route to 10⁻¹⁰); the analytical step is open.
2. **F86d Tier-1 promotion**: structural derivation for Endpoint orbit Q ≈ 2.5. Candidate: SU(2)/Schur-Weyl on F71 first orbit.

**F86e resolved (2026-05-20):** the former third front. σ_0(c=2) is now characterised as a commutator norm ‖[Π_HD1, M_H]‖ = Schur-multiplier norm ‖Π̃_HD1 ⊙ ΔDiff‖ (typed: `SigmaZeroCommutatorNormClaim`, Tier1Derived, bit-exact N=5..8). The asymptote σ_0(∞) ≈ 2.8629 is non-elementary but characterised, not mysterious: the Toeplitz (Fourier-symbol) and Hankel (Nehari) closed-form routes are both ruled out, so σ_0(∞) is a genuine Schur-multiplier-norm constant. Closed-form candidates 2√2, √(41/5), √(8+π/16) ruled out earlier still hold.

### Spawned out (already promoted to own F-numbers)

- **F89** (Path polynomial pipeline, Chebyshev + BigRational, D_k bit-exact): Tier1Derived since 2026-05-15
- **F90** (F86 c=2 K_b ↔ F89 path-(N−1) Hellmann-Feynman bridge): Tier 1 derived (operator-exact identity at N=5..8); the closed-form lift it points to is pending
- **F91 / F92 / F93** (F71-anti-palindromic spectral invariance on γ / J / h)

---

## Seen again (2026-06-11): the fresh-eyes wave and the hardware handover

Three things landed after this inventory's last update (2026-06-01) that bear on its rows; recorded here dated, anchors unchanged unless said so.

1. **The SE-walk population handover was measured on hardware.** The Kingston q13-q14-q15 run (2026-05-31, historical registry id `ibm_ep_onset_may2026`) measured populations only. Its random-Z runner used coherence-rate labels Q_label=J/Γ, while the repository Lindblad convention has γ=Γ/2, hence **Q_Lindblad = 2 Q_label**. The revival stays near 1/N through Q_label=1.5 and is larger at Q_label=2.5 (canonical Q_Lindblad 3→5). This population scan is not an F86-block EP anchor. The walk's own spectral transition is the (1,1)-block EP at the coherence horizon Q*(3) = √2 (`CoherenceHorizonClaim`, Tier 1 derived), which is Q_label ≈ 0.71, below the bracket; the handover itself is a probe-time crossover of the revival amplitude, not a spectral event. Writeup: [The Flow Between Two Singularities](../experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md).
2. **g_eff's irreducible content halved.** The K_b mode mirror (2026-06-10): the antiunitary T = Σ₁∘conj gives K_b(b; mode k) = K_b(b; mode N+1−k) pointwise in Q and t, so **g_eff(c, N, b)'s irreducible content lives in the half-band of {k, N+1−k} orbits** ([the F86b obstruction proof](proofs/PROOF_F86B_OBSTRUCTION.md) "one more surviving symmetry", [the chiral mirror rate law proof](proofs/PROOF_PTF_CHIRAL_MIRROR_RATE_LAW.md) §8). A structural bonus from the same run: on c = 1 blocks K_b ≡ 0 identically (the block generator is e^{−2γt} × unitary), the structural reason the F86 resonance needs c ≥ 2. This does not revise any anchor value, but it halves the target of the F86_block obstruction and of Direction (b'') below: whatever closed form g_eff has or refuses to have, it is a function on mode *pairs*.
3. **The clock and the seam, pinned.** The EP clock's peak angle is γ₀-invariant and provably not the diagonal: θ_peak = 44.3646° (45° would need x = √5; the fixed-τ = 1 peak x_peak = 2.196910 of F86b₁ gives the 0.6354° shortfall; `ExceptionalPointClockTests`). And the Absorption Theorem extensions (2026-06-10) recentre the dissipative seam as L_D = γ(Q − N·I) with the per-eigenmode Rayleigh split and the orthogonal-projector form ([the Absorption Theorem](proofs/PROOF_ABSORPTION_THEOREM.md) Extensions), the structural floor under every Re-rate in this inventory. The slow-light surface gained its basis-free primitive (`SlowLightDistribution`, [compute/RCPsiSquared.Diagnostics/Ptf/](../compute/RCPsiSquared.Diagnostics/Ptf/SlowLightDistribution.cs)): light_l = Tr(Π_V·Δ_l)/dim V, the s\*-surface read without hand-picked endpoint pairs.

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
| **Endpoint (orbit 0)** | **≈ 2.5**             | (c=2..4, N=5..8); 10 of 12 box cells filled (c=4 at N=7, 8) | 2.393–2.591 (population std/mean 2.2% around 2.547, peak-to-peak 7.8%) |

From `PerF71OrbitObservation.cs:50-67` (the envelope Q_peak at γ₀ = 0.05): across 10 different (c, N) combinations, Endpoint orbit Q sits in [2.393, 2.591]. Stability across c and N makes this a candidate anchor not currently in QBasisAnkers.

**Status**: Tier 2 Empirical candidate for promotion. Pending structural derivation.

---

## Empirical observations (not anchor-grade, may revise)

The following values are measurements from the typed C# claims. They are useful for diagnostic comparison but should NOT be treated as structural anchors. Each carries a methodological caveat that may revise the specific number.

### B1. Per-orbit Q_EP values (computed from x_peak)

`C2BlockCpsiQScan.cs:19-21`: Q_EP per orbit = Q_peak / x_peak with x_peak = 2.196910.

| Orbit (N=5 c=2)      | Q_EP value         |
|----------------------|--------------------|
| Endpoint (Orbit 0)   | 1.146              |
| Interior (Orbit 1)   | 0.675              |

**Caveat**: x_peak = 2.196910 is the bare doubled-PTF model's fixed-τ = 1 peak, while Q_peak is the full block's envelope peak, so the quotient is a conversion the predictor defines, not a measured EP position. The bare model's own envelope peaks at x = 1.7141. Treat as derivative observations, not as primary anchors.

### B2. Far maxima (c = 2, N ≥ 7)

From `PerF71OrbitObservation.cs:53-59`. Every orbit has its canonical Q_peak, the first local maximum of the envelope, inside [0.2, 4] (A3 and B3). Past a minimum after it, the envelope rises again, and for these orbits the maximum over a wider window lies higher than the canonical one:

| (c, N) | Orbit                    | Q of the maximum | window      | envelope there vs canonical |
|--------|--------------------------|------------------|-------------|-----------------------------|
| (2, 7) | Orbit 1 (bonds 1, 4)     | 7.255            | [0.2, 20]   | 0.149 vs 0.109              |
| (2, 8) | Orbit 1 (bonds 1, 5)     | 8.225            | [0.2, 20]   | 0.154 vs 0.090              |
| (2, 8) | Orbit 2 (bonds 2, 4)     | 38.21            | [0.2, 60]   | 0.093 vs 0.082              |
| (2, 8) | Central orbit (bond 3)   | 20.09            | [0.2, 60]   | 0.101 vs 0.084              |

**Caveat**: the curves keep further local maxima (c = 2 N = 8 orbit 1 has three on [0.2, 60]), so this reading is a property of the window and is always quoted with it. The two values beyond Q = 20 are read by the second route alone (`simulations/f86_envelope_second_route.py`). What sets the far maxima is not known; they are not Q-anchors.

### B3. Per-bond-class HWHM ratios (c = 2 class values)

The envelope HWHM_left/Q_peak of the class-mean curve at c = 2, against the bare model's envelope floor (`C2BareDoubledPtfClosedForm.EnvelopeHwhmLeftOverXPeak`):

| Bond class | N = 5  | N = 6  | N = 7  | N = 8  | floor        | lift above floor  |
|------------|--------|--------|--------|--------|--------------|-------------------|
| Endpoint   | 0.7766 | 0.7815 | 0.7813 | 0.7808 | 0.7264892771 | +0.050 to +0.055  |
| Interior   | 0.7523 | 0.7569 | 0.7591 | 0.7602 | 0.7264892771 | +0.026 to +0.034  |

**Caveat**: the floor is derived (Tier1Derived, the bare doubled-PTF closed form read on its envelope). The lifts are observed; the analytical derivation is **open** (PROOF_F90_F86C2_BRIDGE.md Item 1'-followup, F89 AT-locked F_a/F_b structure). The Interior value rises with N (0.7621 at N = 12). The values are the reading on the default Q grid; per bond at c = 2, N = 5..8 that reading is within 3.5·10⁻⁵ of the continuous envelope peak (the class-mean curves were not checked against a continuous search).

### B4. Per-sub-class (α, β): a line through the anchors

`F86HwhmAlphaExtraction.cs` (the live recompute; `F86HwhmClosedFormClaim` reads it, no hand-transcription). Form: HWHM_ratio = 0.7264892771 + α·g_eff + β, with g_eff the fit's own coordinate 4.39382/(Q_peak + 2), on the canonical peaks of c = 2, N = 5..8. Whether a line holds is decided per sub-class by how many chain lengths carry an anchor:

| BondSubClass         | anchors        | verdict |
|----------------------|----------------|---------|
| Endpoint             | four N         | overdetermined; misses a line by up to 1.28·10⁻³, so the linear form is not an exact law |
| Flanking             | four N         | overdetermined; misses a line by up to 4.2·10⁻⁴, likewise |
| Mid                  | two N          | two-point line, exact by construction, tests nothing |
| CentralSelfPaired    | two N          | two-point line, likewise |
| Orbit2Escape         | none           | classifies a far maximum (B2); no canonical anchor, the fit throws |
| CentralEscapeOrbit3  | none           | likewise |

**Caveat**: the misses are 12 to 36 times the anchors' offset from the continuous envelope peak, so they are not resolution. The (α, β) values are recomputed live and frozen nowhere; the analytical derivation of the lift from F89/F90 structure remains open.

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
| SE-walk handover (hardware) | Q_label 1.5→2.5; Q_Lindblad 3→5 | `ibm_ep_onset_may2026`: population return near 1/N through the lower sample and larger at the upper (2026-05-31, Kingston); a probe-time crossover above the walk's (1,1)-block EP Q*(3) = √2 | Tier2Verified hardware observation, not a Q_EP anchor |
| Endpoint orbit Q    | 2.5     | PerF71OrbitObservation (range 2.393–2.591) | Tier2Empirical candidate |

Only the anchors above belong in QBasisAnkers and Q_REGIME_ANCHORS. The per-orbit Q_EP, far maxima, HWHM lifts, and fitted (α, β) values are empirical diagnostics that may revise; they live in this inventory and the typed C# claims, not in the anchor map.

## Anchors and cross-refs

- **Population-flow cross-reference (no anchor change):**
  [`THE_FLOW_BETWEEN_TWO_SINGULARITIES`](../experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md)
  follows a fixed-number single-excitation state from its parameter-space singularity, the (1,1)-block
  EP at the coherence horizon Q*(N) (`CoherenceHorizonClaim`), to its 1/N fixed point. That flow EP is a
  distinct object from the c=2 coherence-block Q_EP this inventory tracks (kept separate on purpose,
  not merged), and the tempting coincidences at it (the pair at Re = −2γ, ⟨n_XY⟩ = 1) hold at N = 2, 3
  only (from N = 4 the coalescing pair sits below that floor), so they are not anchors. The flow's
  hardware scan straddles that EP between its two lowest samples; its rate-book-corrected population
  handover lies above it and is not a continuation of an F86 eigenmode. None of this revises the Q_peak / Q_EP / HWHM anchors above.
- **Birth-channel cross-reference (no anchor change):** F86's (n, n+1) blocks are the
  number-changing (Δn = 1) channel that
  [`STERILE_BIRTHCANAL_AND_THE_JUNCTION`](STERILE_BIRTHCANAL_AND_THE_JUNCTION.md) reads by its slowest
  rate. On the open chain at uniform γ₀ every one of these blocks holds N free-fermion modes
  P_n·c_k·P_(n+1) at Re λ = −2γ₀ (the rate floor, not this inventory's envelope floor 0.7264892771),
  on HD = 1 cells only and closed under the hopping, so they never reach the HD = 3 cells
  ([`PROOF_CHAIN_GAP_DOMINANCE`](proofs/PROOF_CHAIN_GAP_DOMINANCE.md); `inspect --root fermionfloor`).
  The Dicke probe is silent exactly on the ones odd under the chain reflection of the coherence and
  reads the even ones in every block checked; in the (1,2) block
  these carry the k-labels of F89's F_a set. The family's span is invariant under L and its adjoint
  alike, so no exceptional point of a block involves a mode of this family. The sterile verdict
  (the slowest rate is Q-flat) reads this layer wherever it is the slowest. The two are readings of one channel, not of two
  channels. F86's K at the (1,2) block splits exactly by this family into floor², rest² and their
  interference (`inspect --root floorsplit`; [`ON_WHAT_CLOSES_ALONE`](../reflections/ON_WHAT_CLOSES_ALONE.md)).
- High-level Q-anchor map: [`Q_REGIME_ANCHORS.md`](Q_REGIME_ANCHORS.md)
- C# F86 KnowledgeBase root: [`compute/RCPsiSquared.Core/F86/`](../compute/RCPsiSquared.Core/F86/)
- Sub-ID partition (this document, top section): F86a / F86b₁..b₅ / F86c / F86d / F86e / F86_block
- Proof hub (sub-ID navigation): [`PROOF_F86_QPEAK.md`](proofs/PROOF_F86_QPEAK.md)
- Per-sub-claim proofs: [`PROOF_F86A_EP_MECHANISM.md`](proofs/PROOF_F86A_EP_MECHANISM.md), [`PROOF_F86B_UNIVERSAL_SHAPE.md`](proofs/PROOF_F86B_UNIVERSAL_SHAPE.md), [`PROOF_F86B_OBSTRUCTION.md`](proofs/PROOF_F86B_OBSTRUCTION.md), [`PROOF_F86C_F71_MIRROR.md`](proofs/PROOF_F86C_F71_MIRROR.md)
- Open derivation gap (HWHM lifts + α/β closed forms): PROOF_F90_F86C2_BRIDGE.md Item 1'-followup; only Direction (b'') survives Item 1' falsification arc (2026-05-11)
- Memory: `project_q_peak_ep_structure`, `project_q_middle_structure`, `project_no_classicalization`, `feedback_f86_needs_agents`; population-flow follow-up (2026-05-31): `project_birth_canal_and_light`, `project_q_as_lifetime_of_the_new`
