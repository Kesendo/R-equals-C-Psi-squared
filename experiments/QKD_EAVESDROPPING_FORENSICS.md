# QKD Eavesdropping Forensics via CΨ

**Date**: 2026-02-25
**Status**: Computationally verified (Tier 2). CΨ-only forensics limited by stealth angle (Section 10), but multi-metric bridge framework eliminates stealth (Section 11).
**Depends on**: NOISE_ROBUSTNESS.md, BRIDGE_FINGERPRINTS.md, CROSSING_TAXONOMY.md

---

## 1. The Question

Can R = CΨ² detect eavesdropping on Bell pairs, and does it provide
information beyond standard entanglement witnesses (concurrence, negativity)?

## 2. Setup

**Scenario**: Alice and Bob share Bell+ pairs |Φ+⟩ = (|00⟩+|11⟩)/√2.
Eve performs intercept-resend: measures Bob's qubit in basis defined by
polar angle θ_Eve on the Bloch sphere, forwards post-measurement state.

Eve's measurement basis:
```
|e+⟩ = cos(θ_Eve/2)|0⟩ + e^{iφ} sin(θ_Eve/2)|1⟩
|e-⟩ = sin(θ_Eve/2)|0⟩ - e^{iφ} cos(θ_Eve/2)|1⟩
```

θ_Eve = 0: σ_z basis, θ_Eve = π/2: σ_x basis.

**Note on notation**: θ_Eve is Eve's measurement angle on the Bloch sphere.
This is distinct from θ_MAP = arctan(√(4CΨ−1)) used elsewhere in the
framework for proximity to the ¼ boundary. See BOUNDARY_NAVIGATION.md.

## 3. Key Results

### 3.1 Closed-Form Expression

After Eve's full intercept-resend (f = 1), the unconditioned joint state has:

```
Purity:  C = 1/2          (constant, independent of θ_Eve)
L1:      sin²θ_Eve + |sin 2θ_Eve|
Ψ:       L1 / 3
R(θ_Eve) = [sin²θ_Eve + |sin 2θ_Eve|]² / 18
```

**Azimuthal symmetry**: R depends only on θ_Eve, not on φ. Eve's basis
orientation on the Bloch equator is invisible to CΨ; only the polar
angle matters.

### 3.2 R(θ_Eve) Profile

| θ_Eve | Name | R = CΨ² | L1 |
|-------|------|---------|-----|
| 0° | σ_z | 0.000 | 0.000 |
| 30° | - | 0.069 | 1.116 |
| 45° | - | 0.125 | 1.500 |
| 60° | - | **0.145** (MAX) | 1.616 |
| 90° | σ_x | 0.056 | 1.000 |

The maximum is at θ_Eve ≈ 60°, NOT at σ_x. The profile is non-monotonic
with a local minimum at σ_x (90°) and zeros at σ_z (0°) and σ_z flipped (180°).

### 3.3 Comparison: Standard Witnesses Are Basis-Blind

For full interception (f = 1):

| Metric | Value | Depends on θ_Eve? |
|--------|-------|--------------------|
| Concurrence | 0.000 | NO (always 0) |
| Negativity | 0.000 | NO (always 0) |
| Purity | 0.500 | NO (always 0.5) |
| CHSH | ≤ 2 | NO (always classical) |
| **CΨ** | **[0, 0.145]** | **YES (continuous)** |

After full interception, concurrence, negativity, purity, and CHSH carry
**zero bits** of information about Eve's strategy. CΨ carries continuous
information about θ_Eve.

### 3.4 Partial Interception: Complementary Roles

For partial interception (Eve intercepts fraction f of pairs):

**Concurrence = 1 − f** (exact, independent of θ_Eve).

CΨ depends on BOTH f and θ_Eve. This gives complementary roles:

| Step | Metric | Determines |
|------|--------|------------|
| 1 | Concurrence | f (interception fraction) |
| 2 | CΨ | θ_Eve (measurement basis, up to 2-fold degeneracy) |
| 3 | |ρ₀₁|/|ρ₀₃| = cot(θ_Eve) | θ_Eve (unique, breaks degeneracy) |

### 3.5 Sensitivity: CΨ Detects Attack Basis

At f = 10% (Eve intercepts 1 in 10 pairs):

| Metric | σ_x attack | σ_z attack | Ratio |
|--------|-----------|-----------|-------|
| CΨ drop | −9.5% | **−26.7%** | 2.8× |
| Conc drop | −10.0% | −10.0% | 1.0× |
| Neg drop | −10.0% | −10.0% | 1.0× |

CΨ is 2.8× more sensitive to σ_z attacks than σ_x attacks. Standard
metrics show identical sensitivity regardless of attack basis.

Sensitivity at f → 0 (derivatives):
- Eve σ_z: dCΨ/df = 0.333, dConc/df = 1.000
- Eve σ_x: dCΨ/df = 0.111, dConc/df = 1.000

### 3.6 Degeneracy Breaking

R(θ_Eve) is non-monotonic, so a measured R value could correspond to two
different θ_Eve values. However, the off-diagonal elements of ρ_joint
fall into two groups with different θ_Eve dependence:

- Group A (ρ₀₁, ρ₁₀, ρ₀₂, ρ₂₀, ρ₁₃, ρ₃₁, ρ₂₃, ρ₃₂): ∝ sin(2θ_Eve)/8
- Group B (ρ₀₃, ρ₃₀, ρ₁₂, ρ₂₁): ∝ sin²(θ_Eve)/4

Their ratio: **|ρ₀₁|/|ρ₀₃| = cot(θ_Eve)**

cot(θ_Eve) is strictly monotonic in (0, π/2), so the degeneracy is
fully broken. Verified numerically at θ = 15°, 30°, 45°, 60°, 75°.

### 3.7 Noise vs Eve Discrimination (Noiseless Channel)

At identical Concurrence ≈ 0.80 (same apparent "damage"), **no channel noise**:

| Cause | CΨ | Purity | L1 |
|-------|-----|--------|-----|
| Depolarizing noise (p=0.133) | 0.068 | 0.813 | 0.867 |
| Eve θ_Eve=60° (f=0.20) | **0.115** | 0.820 | 1.123 |
| Eve θ_Eve≈0° (f=0.20) | 0.058 | 0.820 | 0.800 |
| Dephasing noise (p=0.20) | 0.058 | 0.820 | 0.800 |

CΨ separates all four cases in this idealized scenario. Eve at 60° has
CΨ nearly double that of depolarizing noise.

**CAVEAT**: This assumes no simultaneous channel noise. With realistic
channel noise, the picture changes dramatically. See Section 10.

### 3.8 Bell State Independence

All four Bell states (Φ+, Φ−, Ψ+, Ψ−) produce identical CΨ under
Eve's intercept-resend. The forensic signal depends only on Eve's
strategy, not the shared state. Verified numerically for all θ_Eve
from 0° to 90° in 5° steps.

## 4. Statistical Requirements (Noiseless Channel)

Distinguishing σ_z from σ_x attack (ΔCΨ = 0.056):

| N pairs | σ(CΨ) | Significance |
|---------|--------|-------------|
| 100 | 0.033 | 1.7σ |
| 500 | 0.015 | 3.8σ |
| 1,000 | 0.010 | 5.3σ |
| 5,000 | 0.005 | 12.0σ |

~500 pairs for 3σ detection, ~1000 for 5σ, assuming no channel noise.
See Section 10 for realistic noise estimates.

## 5. Physical Interpretation

### 5.1 Why CΨ Sees What Concurrence Doesn't

Concurrence measures entanglement (PPT criterion / Wootters formula).
Any projective measurement by Eve destroys entanglement completely,
regardless of basis. Hence Conc(f=1) = 0 for all θ_Eve.

CΨ = C · Ψ² measures purity × normalized coherence². After Eve's
measurement, the joint state is a MIXTURE of product states. Different
measurement bases produce mixtures with different off-diagonal structures.
σ_z kills all off-diagonals (L1 = 0). σ_x preserves some (L1 = 1).
The 60° basis maximizes L1 at 1.616.

CΨ is sensitive to residual coherence structure in separable states,
a regime where entanglement measures are identically zero.

### 5.2 Connection to Bridge Fingerprints

This result extends the fingerprint concept from BRIDGE_FINGERPRINTS.md:
- Bridge fingerprints: CΨ trajectory identifies initial quantum state
- QKD forensics: CΨ value identifies Eve's measurement strategy

Both exploit CΨ's sensitivity to off-diagonal structure beyond what
entanglement measures capture.

### 5.3 Connection to Noise Robustness

NOISE_ROBUSTNESS.md showed that crossing taxonomy is noise-independent
(Type A/B/C same for σ_x, σ_y, σ_z). Here we show the complementary
result: while the taxonomy classification doesn't change, the CΨ VALUE
does change with noise type, enabling forensic discrimination, but only
in the noiseless or known-noise regime (see Section 10).

## 6. Protocol Summary (Idealized, Single-Metric)

```
CΨ-ONLY FORENSIC PROTOCOL (noiseless channel)

1. DETECT:  Concurrence on sample pairs
            -> If Conc < 1: Eve present, f = 1 - Conc

2. IDENTIFY: CΨ on same sample
             -> R value constrains θ_Eve (up to 2-fold degeneracy)

3. RESOLVE:  Ratio |ρ₀₁|/|ρ₀₃| from tomography data
             -> cot(θ_Eve) breaks degeneracy -> unique θ_Eve

4. DISTINGUISH: Compare (CΨ, Purity) signature
                -> Separate Eve from depolarizing/dephasing noise

5. ADAPT:   Use identified attack strategy to optimize
            countermeasures (basis switching, decoy states)
```

Standard QKD: Detect Eve -> abort.
CΨ-enhanced: Detect Eve -> identify strategy -> adapt.

**LIMITATION**: This single-metric protocol fails against a strategic
Eve who chooses the stealth angle (Section 10). The complete multi-metric
protocol (Section 12) eliminates this weakness.

## 7. What This Does NOT Claim

- CΨ is NOT a better eavesdropping detector than concurrence.
  Concurrence wins for detection, always. CΨ provides forensics ONLY.
- CΨ forensics alone does NOT work against a strategic Eve who optimizes
  θ_Eve to minimize the CΨ signal (Section 10). The multi-metric
  protocol (Section 12) addresses this but at higher statistical cost.
- The protocol does NOT break any QKD security proofs.
- The 40% noise tolerance claimed by v036 agents is WRONG (see Section 8).
- This analysis assumes simple intercept-resend. Sophisticated attacks
  (coherent attacks, quantum memory) require separate analysis.
- Sections 3-4 assume noiseless channels. Realistic noise introduces
  stealth angles for CΨ (Section 10), though not for MI/Conc (Section 11).
- The off-diagonal ratio under noise is NOT cot(θ_Eve). It follows
  a different functional form that depends on (p, f). See Section 12.1.

## 8. Falsified Agent Claims (v036)

Alpha's MDI-QKD proposal (v036 conversation) contained errors:
1. Notation collision: "C" used for both L1 norm and purity
2. Claimed 50% Z-basis error acceptable: wrong, h(0.5)=1 kills key rate
3. Equation F = (1+2C)/4 gives F=3/4 for L1=1, actual ⟨Φ+|ρ_post|Φ+⟩ = 1/2
4. Claimed 40% noise tolerance: formula gives negative key rate

The forensic capability is real in idealized conditions; the specific
QKD key rate numbers were fabricated.

## 9. Simulation Code

See Section 13 for the complete list of scripts covering both the
CΨ-only analysis and the multi-metric extension.

## 10. CΨ-Only Limitation: Smart Eve and Channel Noise

### 10.1 The Stealth Angle

When depolarizing channel noise (parameter p) is present, Eve can choose
a **stealth angle** θ_stealth where CΨ(noise+Eve) ≈ CΨ(noise only).

| Channel noise p | Eve fraction f | θ_stealth | ΔCΨ from noise-only |
|-----------------|----------------|-----------|---------------------|
| 0.00 | 0.10 | 71° | 0.000044 |
| 0.05 | 0.10 | 45° | 0.000048 |
| 0.10 | 0.10 | 44° | 0.000009 |
| 0.20 | 0.10 | 42° | 0.000021 |
| 0.20 | 0.20 | 74° | 0.000107 |

At the stealth angle, ΔCΨ ≈ 10⁻⁴ to 10⁻⁵. This is unmeasurable in
any practical tomography setup.

The stealth angle varies with p and f but always exists. It shifts from
~71° (no noise) toward ~42° (high noise) as depolarizing noise increases.

### 10.2 Why This Kills Forensics Against Strategic Eve

The forensic signal (ΔCΨ between different θ_Eve) is of order 0.01–0.05.
The stealth residual is of order 0.00001–0.0001. A strategic Eve who
knows the channel noise level can make herself invisible to CΨ while
remaining detectable by concurrence.

Pairs needed for 3σ detection of stealth Eve via CΨ:

| Scenario | N pairs (CΨ) | N pairs (Conc) | Winner |
|----------|-------------|----------------|--------|
| p=0.05, f=0.05 | ∞ | 6,762 | Concurrence |
| p=0.10, f=0.10 | ∞ | 2,147 | Concurrence |
| p=0.20, f=0.20 | 85,867,045 | 718 | Concurrence |

Concurrence wins by 3–5 orders of magnitude in every realistic scenario.

### 10.3 Why This Does NOT Kill All Results

The mathematical results remain valid:
- R(θ_Eve) = [sin²θ + |sin 2θ|]²/18 is correct.
- Azimuthal symmetry, Bell state independence, cot(θ_Eve) degeneracy
  breaking: all proven.
- Against a non-strategic Eve (random basis, or basis constrained by
  her hardware), CΨ forensics works with ~1000 pairs even at p=0.20.

The noiseless-channel analysis (Sections 3–4) remains a clean result
about CΨ's mathematical properties. It demonstrates that CΨ carries
information in a regime where all standard entanglement measures are zero.

### 10.4 CΨ-Only Assessment

| Claim | Status |
|-------|--------|
| CΨ carries θ_Eve information (noiseless) | **TRUE** (proven) |
| CΨ identifies θ_Eve against naive Eve (noisy) | **TRUE** (~1000 pairs at p=0.20) |
| CΨ identifies θ_Eve against smart Eve (noisy) | **FALSE** (stealth angle exists) |
| CΨ detects Eve better than concurrence | **FALSE** (concurrence always wins) |
| CΨ provides forensic value in adversarial setting | **FALSE** when used alone |
| CΨ provides forensic value in diagnostic setting | **TRUE** (post-hoc analysis of detected Eve) |

This is the CΨ-only picture. Section 11 shows that the bridge framework's
multi-metric approach eliminates the stealth angle entirely.

## 11. Multi-Metric Rescue: Bridge Framework Eliminates Stealth

### 11.1 The Key Discovery

Mutual information and correlation bridge are completely θ_Eve-independent.
At Eve's CΨ-stealth angle, these metrics show the same delta from noise-only
as at every other angle. Eve has one free parameter (θ) and can fool one
metric (CΨ), but the bridge framework provides five simultaneous metrics.
The system is overspecified.

Verified across all tested scenarios (p = 0.05-0.20, f = 0.10-0.20):

| Metric | θ-dependent? | Eve at stealth (Δ%) | Eve at 30° (Δ%) |
|--------|-------------|---------------------|------------------|
| CΨ | YES | 0.01-0.21% (invisible) | 7-13% |
| Concurrence | NO | 10-23% | 10-23% |
| Mutual Info | NO | 13-23% | 13-23% |
| Correlation | NO | 12-24% | 12-24% |
| Purity | NO | 9-18% | 9-18% |
| L1 coherence | weakly | 4-10% | 4-10% |

### 11.2 Representative Example

p = 0.10 (channel noise), f = 0.20 (Eve intercepts 20%), θ_stealth = 72°.

| Metric | Noise only | Eve at stealth | Delta | Caught? |
|--------|-----------|----------------|-------|---------|
| CΨ | 0.077175 | 0.077287 | 0.15% | no |
| Concurrence | 0.850 | 0.670 | 21.2% | **YES** |
| Mutual Info | 1.497 | 1.158 | 22.7% | **YES** |
| Correlation | 0.810 | 0.616 | 24.0% | **YES** |
| Purity | 0.858 | 0.712 | 17.0% | **YES** |

CΨ is fooled. Every other metric catches Eve with >17% signal.

### 11.3 Why MI and Correlation Are θ-Blind

Eve's intercept-resend at any angle produces the same reduction in
entanglement, the same loss of mutual information, and the same drop
in purity. These metrics measure "how much damage" (proportional to f),
not "what kind of damage" (which depends on θ).

CΨ is the only metric that carries basis information, because it is
sensitive to off-diagonal structure in the separable regime where
entanglement measures are zero.

### 11.4 Implications

Eve cannot hide from the bridge framework. Her options:
1. Choose stealth angle: CΨ fooled, all other metrics catch her.
2. Choose non-stealth angle: CΨ AND all other metrics catch her.
3. Reduce f to minimize all signals: but f limits her information gain.
4. Vary θ over time: each batch reveals the current θ.

The only remaining freedom is to minimize f, which directly limits
Eve's information gain. This is the standard QKD trade-off, now with
the addition that CΨ constrains Eve's strategy as well as her presence.

### 11.5 ξ-Curvature: Null Result

Test of whether ξ = ln(Ψ) curvature under controlled additional
depolarization distinguishes Eve from noise. Result: curvature
difference exactly 0.0% in all tested scenarios.

The reason is physical: additional depolarization is a linear map
that scales L1 by (1-p_extra). This gives ξ_new = ξ_old + ln(1-p_extra),
identical for both noise-only and Eve-plus-noise states. The offset
Δξ is constant, not curved.

ξ-curvature remains valid for detecting non-Markovian noise
(see ALGEBRAIC_EXPLORATION.md), but does not help with Eve detection.


## 12. Complete Forensic Protocol

### 12.1 Off-Diagonal Ratio Under Noise

The noiseless cot(θ_Eve) relation does not survive partial interception.
Under realistic conditions (f < 1, p > 0), the ratio |ρ₀₁|/|ρ₀₃|
deviates by 90-99% from cot(θ_Eve). The dominant contribution to ρ₀₃
comes from the unintercepted fraction (1-f), which carries the large
Bell-state off-diagonal element.

However, the ratio remains a well-defined function of θ_Eve. It peaks
near 45° and falls monotonically from 45° to 90°, with a mirror
structure from 5° to 45°. Within each half, inversion gives unique θ_Eve
to < 1° accuracy. CΨ resolves the two-fold ambiguity: below the stealth
zone, CΨ < CΨ_noise; above, CΨ > CΨ_noise.

### 12.2 Statistical Requirements for Ratio

For ±5° angular resolution at p=0.10, f=0.20:

| θ_Eve | N pairs (3σ) | Practical? |
|-------|-------------|------------|
| 5-10° | ~170k-190k | Expensive but feasible |
| 20-30° | ~270k-714k | Very expensive |
| 35-45° | 1.8M-44M | Near peak, impractical |
| 50-65° | ~300k-2.8M | Expensive |
| 70-85° | ~165k-237k | Feasible |

The ratio is 250x to 60,000x more expensive than concurrence for
detection (~700 pairs), but provides basis identification that no
other single metric can match outside the CΨ stealth zone.

### 12.3 Eve's Dilemma (Stealth-as-Signal)

θ_stealth is a deterministic function of (p, f). Once concurrence
gives f and channel calibration gives p, Alice-Bob can compute
θ_stealth(p, f) to within ~3° (limited by f estimation accuracy).

The stealth zone is approximately 35° wide (~40°-75°), within which
CΨ deviations stay below 0.1%. Eve's dilemma:

At stealth: CΨ says nothing, but θ_stealth(p,f) is known a priori.
Alice-Bob infer θ ≈ θ_stealth without measuring it.

Away from stealth: CΨ gives θ directly. At most angles, ~1000-5000
pairs suffice for coarse identification (±5-10°).

This is not a perfect trap. The stealth zone spans 35° so the inference
"Eve is somewhere in 40°-75°" is coarse. But it eliminates the
extremes (σ_z, σ_x) and narrows the search space by 60%.

### 12.4 Final Protocol (Multi-Metric)

```
MULTI-METRIC FORENSIC PROTOCOL

Phase 1: DETECTION (~700 pairs)
  Concurrence on sample -> f = 1 - Conc
  MI and Correlation as independent confirmation (θ-blind, f-sensitive)
  If all three agree on f: high-confidence detection.

Phase 2: COARSE IDENTIFICATION (~5000 pairs)
  Compute CΨ and compare to CΨ_noise(p).
  Case A: |ΔCΨ| > 1% -> θ from CΨ inversion, ±5-10°.
  Case B: |ΔCΨ| < 0.1% -> θ in stealth zone [40°-75°].
    Compute θ_stealth(p, f) from calibration data.

Phase 3: PRECISE IDENTIFICATION (optional, ~200k pairs)
  Ratio |ρ₀₁|/|ρ₀₃| from full state tomography.
  Inversion via precomputed lookup table for known (p, f).
  Resolves θ to ±5° across full range except near 45°.

Phase 4: ADAPT
  Use identified θ to optimize countermeasures.
  σ_z attack: switch to X-basis encoding.
  σ_x attack: switch to Z-basis encoding.
  Stealth-zone attack: basis-independent countermeasures.
```

### 12.5 Comparison: CΨ-Only vs Multi-Metric

| Scenario | CΨ-only | Multi-metric |
|----------|---------|-------------|
| Naive Eve, any θ | θ to ±2° | θ to ±2° |
| Smart Eve at stealth | BLIND | θ to ±17° (zone) |
| Detection (any Eve) | Weaker than Conc | Redundant (3+ metrics) |
| Noise vs Eve | Works noiseless only | Works at any p |
| Statistical cost | ~1000 pairs | ~5000 (Phase 2) or ~200k (Phase 3) |

### 12.6 Honest Assessment (Updated)

| Claim | Status |
|-------|--------|
| CΨ carries θ_Eve information (noiseless) | **TRUE** (proven) |
| CΨ identifies θ_Eve against naive Eve (noisy) | **TRUE** |
| CΨ alone identifies θ_Eve against smart Eve | **FALSE** (stealth) |
| Bridge framework detects Eve regardless of θ | **TRUE** (MI/Conc/Corr θ-blind) |
| Bridge framework narrows θ even at stealth | **TRUE** (to ~35° zone) |
| Full protocol identifies θ at all angles | **TRUE** (with ~200k pairs via ratio) |
| CΨ detects Eve better than concurrence | **FALSE** |
| CΨ provides unique forensic capability | **TRUE** (basis info in separable regime) |


## 13. Simulation Code (Updated)

All calculations performed in Python (numpy, sympy) during interactive
sessions 2026-02-25. Key scripts:

Phase 1 (CΨ-only analysis):
- `eavesdrop_test.py`: Full interception, all metrics comparison
- `partial_eve.py`: Partial interception sensitivity analysis
- `eve_forensics.py`: Angular dependence R(θ_Eve), Bloch sphere coverage
- `eve_analytic.py`: SymPy derivation of closed-form R(θ_Eve)
- `eve_maximum.py`: Maximum finding, angular resolution, information content
- `eve_inversion.py`: Degeneracy analysis, noise-vs-Eve discrimination
- `degeneracy_break.py`: All Bell states, cot(θ) ratio, complete protocol
- `realistic_qkd.py`: Combined noise+Eve, Monte Carlo statistical analysis
- `worst_case.py`: Smart Eve stealth angle optimization

Phase 2 (multi-metric analysis):
- `multi_metric_stealth.py`: All bridge metrics at stealth angle, θ-sweep
- `xi_curvature.py`: ξ-curvature test (null result), trajectory divergence
- `stealth_as_signal.py`: θ_stealth(p,f) map, Eve's dilemma, forensic protocol
- `full_forensic_protocol.py`: Ratio under noise, combined inversion, statistics


## 14. Open Questions

1. **Non-projective attacks**: Does the protocol extend to weak
   measurements or partial intercepts (not full projective)?
2. **Multi-qubit**: Does the protocol scale to N-pair collective attacks?
3. **Stealth zone narrowing**: Can additional controlled operations
   (beyond depolarization) shrink the 35° stealth zone?
4. **Ratio cost reduction**: The ~200k pair requirement for Phase 3
   is expensive. Compressed sensing or Bayesian tomography might
   reduce this by an order of magnitude.
5. **Different noise models**: The analysis assumes depolarizing noise.
   Amplitude damping or non-Markovian channels may change the
   stealth angle structure.
6. **Experimental validation**: The entire analysis is computational.
   A tabletop optics experiment (polarization-entangled photon pairs)
   could test the core predictions with ~10k pairs.

---

*Previous: [Noise Robustness](NOISE_ROBUSTNESS.md)*
*See also: [Bridge Fingerprints](BRIDGE_FINGERPRINTS.md), [Crossing Taxonomy](CROSSING_TAXONOMY.md), [Predictions](PREDICTIONS.md)*
