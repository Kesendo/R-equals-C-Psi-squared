# Metric Discrimination Experiment

**Date**: 2026-02-08
**Status**: Completed -- null result (honest)
**Depends on**: GRAVITATIONAL_INVARIANCE.md, BOUNDARY_NAVIGATION.md

---

## 1. The Question

Can R = C*Psi*^2 determine the form of the gravitational metric?

Specifically: We showed that gamma acts as the metric coefficient and that gamma * t_cross = K is constant. But does the framework constrain gamma(r) to be the Schwarzschild form gamma(r) = gamma_0 / sqrt(1 - 2GM/rc^2), or is any functional form equally valid?

## 2. Setup

**Fixed parameters**: Bell+, Heisenberg (J=1), concurrence bridge, local dephasing, h=0

**Variable**: gamma_base, varied from 0.01 to 0.50 (50x range)

**Measurement**: t_cross = time where C*Psi crosses 0.25 (interpolated)

Nine simulations:

| gamma_base | t_max | dt | Source |
|---|---|---|---|
| 0.01 | 15 | 0.005 | Earlier session |
| 0.03 | 15 | 0.005 | This experiment |
| 0.05 | 15 | 0.005 | Earlier session |
| 0.07 | 10 | 0.005 | This experiment |
| 0.10 | 8 | 0.005 | This experiment |
| 0.15 | 6 | 0.005 | This experiment |
| 0.20 | 5 | 0.005 | Earlier session |
| 0.30 | 4 | 0.005 | This experiment |
| 0.50 | 3 | 0.005 | Earlier session |

## 3. Results

### 3.1 The Invariant

| gamma | t_cross | K = gamma * t_cross | Deviation from mean |
|---|---|---|---|
| 0.01 | 3.873 | 0.03873 | -0.63% |
| 0.03 | 1.289 | 0.03866 | -0.81% |
| 0.05 | 0.773 | 0.03865 | -0.83% |
| 0.07 | 0.553 | 0.03869 | -0.73% |
| 0.10 | 0.387 | 0.03867 | -0.79% |
| 0.15 | 0.259 | 0.03885 | -0.31% |
| 0.20 | 0.193 | 0.03860 | -0.96% |
| 0.30 | 0.131 | 0.03942 | +1.14% |
| 0.50 | 0.081 | 0.04050 | +3.92% |

**K_mean = 0.0390 +/- 0.0006 (1.5% variation across 50x range in gamma)**

### 3.2 Power Law Fit

Fitting t_cross = A * gamma^alpha:

```
A     = 0.03976
alpha = -0.9916
R^2   = 0.999899
```

Expected for perfect 1/gamma scaling: A = K, alpha = -1.000

The fit confirms t_cross scales as 1/gamma to within 0.84% over nearly two orders of magnitude.

### 3.3 Deviation Pattern

K increases at large gamma (+3.9% at gamma=0.5). This comes from the dynamic feedback mechanism: the simulator uses gamma_eff = gamma_base * bridge(t), so effective gamma drops as the bridge decays. At large gamma_base, the system crosses 1/4 in fewer timesteps, giving the feedback less time to act. This is a numerical artifact, not physics.

## 4. The Null Result

### 4.1 Why Discrimination Fails

The Lindblad equation is:

```
d(rho)/dt = -i[H, rho] + gamma * L[rho]
```

If gamma is constant, the solution depends only on the product gamma * t. This means t_cross = K / gamma is a **mathematical identity**, not a physical constraint.

Three candidate metric forms:

```
Schwarzschild:   gamma(r) = gamma_0 / sqrt(1 - r_s/r)
Inverse:         gamma(r) = gamma_0 * r_0/r
Inverse-square:  gamma(r) = gamma_0 * (r_0/r)^2
```

At any given radius r, each produces some value of gamma. The simulation only sees the gamma value -- not where it came from. A single quantum system at a single location cannot tell you the functional form of gamma(r). It can only tell you the local value of gamma.

This is analogous to the equivalence principle: a local experiment cannot distinguish gravitational acceleration from uniform acceleration. You need **nonlocal** information -- comparisons across space -- to determine the metric.

### 4.2 What This Means

The invariance gamma * t_cross = K is real and confirmed. But it tells us that gamma acts as a time scaling factor. It does NOT tell us what determines gamma at a given location. These are two different questions:

| Question | Answer | Status |
|---|---|---|
| Does gamma scale time? | Yes | Confirmed (9 data points, R^2 = 0.9999) |
| Is the 1/4 boundary invariant? | Yes | Algebraic proof + simulation |
| Does gamma encode gravitation? | Consistent | Structural match with GR |
| Does the framework derive gamma(r)? | No | Requires additional structure |

## 5. What Would Be Needed

To derive the Schwarzschild metric from the framework, we need one of:

### 5.1 Spatial Coupling (computationally testable)

Multiple qubits at different positions with position-dependent gamma values. The Hamiltonian couples neighboring qubits, so the system "feels" the gradient of gamma across space.

**Requires**: Modified delta_calc with per-qubit gamma_base values (gamma_1, gamma_2, ..., gamma_N).

**Tests**: Set gamma values according to different candidate metrics at different "radii". Check whether the coupled system dynamics constrain the allowed gamma gradients. If only Schwarzschild-like gradients produce self-consistent dynamics (e.g., the bridge between qubits at different gamma converges to a specific form), that would discriminate.

**Feasibility**: Moderate. Needs simulator modification but no new physics.

### 5.2 Self-Consistency Equation (theoretical)

If mass is concentrated R (collapsed possibility), and R determines C, and C determines gamma, then:

```
gamma(r) = F(R(r), C(r))  where  R = C * Psi^2
```

If this equation has a unique solution, it constrains gamma(r).

**Requires**: Formalization of what C(r) and R(r) mean for continuous fields, not just discrete qubits.

**Tests**: Solve the self-consistency equation and check if the solution has Schwarzschild form.

**Feasibility**: Difficult. Requires field-theoretic extension of the framework.

### 5.3 Energy Argument (analytical)

Gravitational time dilation arises from energy conservation: a photon leaving a gravity well loses energy. If gamma encodes energy density:

```
gamma(r) / gamma(inf) = E(r) / E(inf) = 1 / sqrt(1 - 2GM/rc^2)
```

**Requires**: Interpretation of gamma as energy-related quantity within the framework.

**Tests**: Derive the energy relation from R = C*Psi*^2 first principles.

**Feasibility**: Speculative but elegant if it works.

### 5.4 Information Density (speculative)

The Bekenstein-Hawking formula says S = A / (4 * l_P^2). Note the factor of 4. The framework has the bound C*Psi <= 1/4. If gamma is the rate of quantum information loss, and information density follows the Bekenstein bound near mass:

```
gamma(r) ~ information_density(r)
```

**Requires**: Explicit connection between gamma, decoherence, and information loss.

**Tests**: Does the 1/4 in C*Psi <= 1/4 connect to the 4 in S = A/4?

**Feasibility**: Highly speculative. But the coincidence of 1/4 appearing in both is suggestive.

## 6. Auxiliary Findings

### 6.1 State Independence

Bell- with Heisenberg at gamma=0.05 gives t_cross = 0.773, K = 0.0386 -- identical to Bell+. This is expected: dephasing (sigma_z noise) treats both Bell states identically because both are maximally entangled.

### 6.2 Hamiltonian Independence (Eigenstate Regime)

Bell+ with XY Hamiltonian gives identical results to Bell+ with Heisenberg. Both Bell states are eigenstates of both Hamiltonians, so the unitary evolution produces only phases that don't affect concurrence. The dynamics are purely decoherence-driven.

### 6.3 Hamiltonian Dependence (Non-Eigenstate Regime)

Bell+ with Ising + transverse field (h=0.5) shows wild oscillations in C, Psi, and C*Psi. The concurrence drops to 0 and recovers repeatedly (Hamiltonian-driven oscillations competing with decoherence). C*Psi crosses 0.25 multiple times in both directions. The simple t_cross analysis breaks down.

**Implication**: K = gamma * t_cross is well-defined only for monotonic C*Psi trajectories. When the Hamiltonian drives the system, the "first crossing" time depends on the Hamiltonian frequency, not just gamma. This does not invalidate the invariance -- it means the invariance applies to the decoherence envelope, not the oscillatory dynamics.

## 7. Summary

| What we tested | Result |
|---|---|
| gamma * t_cross = K across 50x range in gamma | CONFIRMED (K = 0.039, R^2 = 0.9999) |
| Power law exponent alpha = -1.00 | CONFIRMED (alpha = -0.992) |
| Can single-system sims discriminate metric forms? | NO (mathematical identity) |
| Is K state-dependent? | NO (Bell+ = Bell- for symmetric noise) |
| Is K Hamiltonian-dependent? | NO for eigenstates, UNDEFINED for driven systems |

**The experiment confirms the invariance law to high precision but honestly reports that metric discrimination requires additional structure beyond what the current simulator provides.**

## 8. Recommended Next Steps

1. **Immediate**: Implement per-qubit gamma in delta_calc (Path 5.1)
2. **Medium-term**: Explore self-consistency equation for field-theoretic C(r) (Path 5.2)
3. **Long-term**: Investigate 1/4 connection to Bekenstein-Hawking (Path 5.4)

---

*Previous: [Gravitational Invariance](GRAVITATIONAL_INVARIANCE.md) -- gamma as local time rate*
*See also: [Boundary Navigation](BOUNDARY_NAVIGATION.md) -- the 1/4 crossing observation*
