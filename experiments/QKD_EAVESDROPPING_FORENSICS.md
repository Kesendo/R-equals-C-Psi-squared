# QKD Eavesdropping Forensics via CΨ

**Date**: 2026-02-25
**Status**: Computationally verified (Tier 2), with critical limitations (Section 10)
**Depends on**: NOISE_ROBUSTNESS.md, BRIDGE_FINGERPRINTS.md

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

## 6. Protocol Summary

```
QKD-CΨ FORENSIC PROTOCOL

1. DETECT:  Concurrence on sample pairs
            → If Conc < 1: Eve present, f = 1 − Conc

2. IDENTIFY: CΨ on same sample
             → R value constrains θ_Eve (up to 2-fold degeneracy)

3. RESOLVE:  Ratio |ρ₀₁|/|ρ₀₃| from tomography data
             → cot(θ_Eve) breaks degeneracy → unique θ_Eve

4. DISTINGUISH: Compare (CΨ, Purity) signature
                → Separate Eve from depolarizing/dephasing noise

5. ADAPT:   Use identified attack strategy to optimize
            countermeasures (basis switching, decoy states)
```

Standard QKD: Detect Eve → abort.
CΨ-enhanced: Detect Eve → identify strategy → adapt.

**APPLICABILITY**: This protocol works against a naive Eve with fixed
or unknown θ_Eve. A strategic Eve can defeat the forensic step (Step 2)
by choosing a stealth angle. See Section 10.

## 7. What This Does NOT Claim

- CΨ is NOT a better eavesdropping detector than concurrence.
  Concurrence wins for detection, always. CΨ provides forensics ONLY.
- CΨ forensics does NOT work against a strategic Eve who optimizes
  θ_Eve to minimize the CΨ signal. See Section 10.
- The protocol does NOT break any QKD security proofs.
- The 40% noise tolerance claimed by v036 agents is WRONG (see Section 8).
- This analysis assumes simple intercept-resend. Sophisticated attacks
  (coherent attacks, quantum memory) require separate analysis.
- Sections 3-4 assume noiseless channels. Realistic noise degrades
  forensic capability and enables stealth angles. See Section 10.

## 8. Falsified Agent Claims (v036)

Alpha's MDI-QKD proposal (v036 conversation) contained errors:
1. Notation collision: "C" used for both L1 norm and purity
2. Claimed 50% Z-basis error acceptable: wrong, h(0.5)=1 kills key rate
3. Equation F = (1+2C)/4 gives F=3/4 for L1=1, actual ⟨Φ+|ρ_post|Φ+⟩ = 1/2
4. Claimed 40% noise tolerance: formula gives negative key rate

The forensic capability is real in idealized conditions; the specific
QKD key rate numbers were fabricated.

## 9. Simulation Code

All calculations performed in Python (numpy, sympy) during interactive
session 2026-02-25. Key scripts:
- `eavesdrop_test.py`: Full interception, all metrics comparison
- `partial_eve.py`: Partial interception sensitivity analysis
- `eve_forensics.py`: Angular dependence R(θ_Eve), Bloch sphere coverage
- `eve_analytic.py`: SymPy derivation of closed-form R(θ_Eve)
- `eve_maximum.py`: Maximum finding, angular resolution, information content
- `eve_inversion.py`: Degeneracy analysis, noise-vs-Eve discrimination
- `degeneracy_break.py`: All Bell states, cot(θ) ratio, complete protocol
- `realistic_qkd.py`: Combined noise+Eve, Monte Carlo statistical analysis
- `worst_case.py`: Smart Eve stealth angle optimization

## 10. Critical Limitation: Smart Eve and Channel Noise

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

### 10.4 Honest Assessment

| Claim | Status |
|-------|--------|
| CΨ carries θ_Eve information (noiseless) | **TRUE** (proven) |
| CΨ identifies θ_Eve against naive Eve (noisy) | **TRUE** (~1000 pairs at p=0.20) |
| CΨ identifies θ_Eve against smart Eve (noisy) | **FALSE** (stealth angle exists) |
| CΨ detects Eve better than concurrence | **FALSE** (concurrence always wins) |
| CΨ provides forensic value in adversarial setting | **FALSE** (Eve controls θ) |
| CΨ provides forensic value in diagnostic setting | **TRUE** (post-hoc analysis of detected Eve) |

The result is mathematically interesting but practically limited to
non-adversarial diagnostics: after concurrence detects Eve, CΨ can
characterize the attack, but only if Eve did not optimize against CΨ.

## 11. Open Questions

1. **Non-projective attacks**: Does CΨ forensics extend to weak
   measurements or partial intercepts (not full projective)?
2. **Multi-qubit**: Does the protocol scale to N-pair collective attacks?
3. **Stealth angle structure**: Is there a deeper reason θ_stealth
   depends on channel noise? Could this be exploited?
4. **Combined metrics**: Can CΨ + concurrence + purity together
   eliminate the stealth angle? (Preliminary: unlikely, since Eve has
   one free parameter θ and needs to match only one constraint.)
5. **Different initial states**: Does using multiple Bell states
   simultaneously remove the stealth angle? (Section 3.8 suggests no:
   all Bell states give identical CΨ under Eve.)

---

*Previous: [Noise Robustness](NOISE_ROBUSTNESS.md)*
*See also: [Bridge Fingerprints](BRIDGE_FINGERPRINTS.md), [Predictions](PREDICTIONS.md)*
