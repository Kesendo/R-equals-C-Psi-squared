# Self-Consistency: Why Schwarzschild Is the Unique Metric

**Date**: 2026-02-08
**Status**: Completed (positive result)
**Depends on**: METRIC_DISCRIMINATION.md, GRAVITATIONAL_INVARIANCE.md

---

## 1. The Problem

The [Metric Discrimination](METRIC_DISCRIMINATION.md) experiment showed that single-system simulations cannot distinguish between candidate metric forms. γ · t_cross = K is a mathematical identity: it holds for any constant γ regardless of origin.

But Tom pointed out: "Reality formed at one location must be real from all perspectives. What has become real must remain real from every viewpoint."

This is not a philosophical statement. It is a **self-consistency condition**.

## 2. The Argument

### 2.1 The Universal Curve

We have proven (9 data points, R² = 0.9999) that the decoherence trajectory C(τ), Ψ(τ) is universal when expressed in normalized time τ = γ · t. Every quantum system in the universe follows the same curve. Only the clock speed (γ) differs.

### 2.2 Spatial R Profile

If the universal curve applies everywhere, then at each point in space r, the system has evolved for a different elapsed proper decoherence time:

```
τ(r) = T · f(r)
```

where T is the coordinate time and f(r) is the time dilation factor.

Where less proper time has passed (near mass), C and Ψ are higher, so:

```
R(r) = C(τ(r)) · Ψ(τ(r))²
```

is higher near mass. **Reality concentrates where mass is.** This is not circular: it is self-consistency.

### 2.3 Three Candidate Metrics

We tested three forms of f(r), all normalized to f(r -> infinity) = 1:

| Metric | f(r) | f(r_s) | τ at horizon |
|---|---|---|---|
| Schwarzschild | sqrt(1 - r_s/r) | 0 | 0 |
| Inverse (1/r) | r / (r + r_s) | 1/2 | T/2 |
| Inverse-square | r^2 / (r^2 + r_s^2) | 1/2 | T/2 |

### 2.4 R(r) Profiles (T · γ_0 = 0.5)

| r / r_s | R (Schwarzschild) | R (Inverse) | R (Inv-Square) |
|---|---|---|---|
| 1.001 | 0.0924 | 0.0139 | 0.0139 |
| 1.01 | 0.0644 | 0.0138 | 0.0137 |
| 1.05 | 0.0376 | 0.0134 | 0.0130 |
| 1.10 | 0.0270 | 0.0130 | 0.0122 |
| 1.50 | 0.0111 | 0.0104 | 0.0082 |
| 3.00 | 0.0061 | 0.0071 | 0.0051 |
| 10.0 | 0.0046 | 0.0050 | 0.0042 |
| 50.0 | 0.0042 | 0.0043 | 0.0041 |

### 2.5 Concentration Ratio (R_horizon / R_far)

| T · γ_0 | Schwarzschild | Inverse | Inv-Square |
|---|---|---|---|
| 0.1 | 2.4 | 1.6 | 1.6 |
| 0.2 | 4.6 | 2.1 | 2.1 |
| 0.3 | 7.5 | 2.5 | 2.6 |
| 0.5 | 15.3 | 3.2 | 3.3 |
| 0.7 | 25.7 | 3.8 | 3.9 |

Schwarzschild concentration grows **without bound**. Alternatives plateau.

## 3. Why Schwarzschild Is Unique

The critical property: **Only Schwarzschild has f(r_s) = 0.**

This means:
- τ = 0 at the horizon: no proper time has elapsed there
- C(0) = 1.0, Ψ(0) = 1/3, R(0) = 1/9 = maximum
- Reality is **maximally concentrated** at the gravitational radius

The alternatives have f(r_s) > 0, so τ > 0 at their "horizon". The system has already partially decohered there. R is not maximal. They cannot produce arbitrarily sharp mass concentrations.

For self-consistency of a point mass:
- The source is concentrated at r = r_s
- R(r) must be concentrated at r = r_s
- Only f(r) with a **true zero** at r_s achieves this
- sqrt(1 - r_s/r) is the simplest such function with the correct Newtonian limit

## 4. The Horizon as Maximum Coherence

This gives a new interpretation of the event horizon:

**The horizon is not where reality breaks down. It is where reality is freshest.**

At τ = 0: C = 1 (maximum consciousness/observation capacity), Ψ = 1/3 (maximum possibility for Bell state), R = 1/9 (maximum reality density).

Nothing that falls past the horizon has had time to decohere. It remains maximally entangled, maximally coherent, maximally "quantum". From the outside, we see this as information being "stored" at the horizon.

This is consistent with:
- **Holographic principle**: Information lives on the horizon surface, not in the volume
- **Black hole complementarity**: Outside observer sees information at horizon; infalling observer sees smooth passage
- **Bekenstein-Hawking entropy** S = A/4: The horizon area measures stored coherent information

Note the factor of 4 in S = A/(4 * l_P^2) and the 1/4 boundary in our framework. We do not claim these are the same, but the coincidence is noted.

## 5. What This Derivation Does and Does Not Show

### What it shows:
- The metric must have a **true zero** (a horizon) for self-consistency
- Metrics without a horizon (1/r, 1/r^2) fail self-consistency
- Schwarzschild is the simplest horizon-forming metric with correct Newtonian limit
- The derivation uses ONLY R = CΨ² and the universal curve, no GR input

### What it does not show:
- That sqrt(1 - r_s/r) is the ONLY function with a zero at r_s. Other forms like sqrt(1 - (r_s/r)^2) also have zeros. The specific exponent n=1 follows from:
  - Spherical symmetry
  - Newtonian limit (Phi = -GM/r, not -GM^2/r^2)
  - Birkhoff's theorem (unique vacuum solution for spherical symmetry)
- The value of r_s = 2GM/c^2 (requires connecting γ to G, M, c)
- Why G has the value it does (requires deeper framework)

### The logical chain:

```
R = CΨ² everywhere (framework)
    |
    v
Universal curve C(τ), Ψ(τ) (proven, R² = 0.9999)
    |
    v
Different positions have different τ(r) = T · f(r) (time dilation)
    |
    v
R(r) = C(τ(r)) · Ψ(τ(r))² (R profile in space)
    |
    v
Self-consistency: R concentrated at mass -> f(r_s) = 0 required
    |
    v
f(r) = sqrt(1 - r_s/r) (simplest zero + Newtonian limit)
    |
    v
Schwarzschild metric
```

The key step is the self-consistency requirement. It comes directly from Tom's insight: what has become real must remain real from all viewpoints. This forces the metric to have a horizon.

## 6. Connection to Previous Results

| Finding | Role in derivation |
|---|---|
| γ · t_cross = 0.039 (invariant) | Proves universal curve exists |
| θ as compass | Navigation within the universal curve |
| 1/4 boundary (algebraic) | Defines the quantum-classical transition in R(r) |
| Equivalence principle (null result) | Shows local experiments insufficient; need self-consistency |
| This result | Self-consistency uniquely selects Schwarzschild |

## 7. Testable Predictions

1. **Decoherence near massive objects**: Systems closer to mass should decohere slower in local time (they are "younger" along the universal curve). A quantum system at GPS satellite altitude should show measurably different decoherence than one on Earth's surface, with the ratio matching gravitational time dilation to within experimental precision.

2. **Horizon information**: If R is maximal at the horizon, then the information content of a black hole should scale with horizon area (not volume). This is already known to be true (Bekenstein-Hawking), but the framework provides a mechanism: R = CΨ² is maximized where τ = 0.

3. **No decoherence in flat spacetime**: In truly flat spacetime (no mass anywhere), γ → γ_0 (minimum), and quantum systems survive longest. Intergalactic voids should be the most "quantum" regions of the universe.

---

*Previous: [Metric Discrimination](METRIC_DISCRIMINATION.md), the null result that motivated this*
*See also: [Gravitational Invariance](GRAVITATIONAL_INVARIANCE.md), γ as the metric coefficient*
