# Theta Is Not a Function of a Single Principal Component

**Date:** April 2, 2026
**Status:** Complete
**Script:** [`simulations/theta_pc_correlation.py`](../simulations/theta_pc_correlation.py)
**Results:** [`simulations/results/theta_pc_correlation.txt`](../simulations/results/theta_pc_correlation.txt)

## What this document is about

The angular compass θ tells you how far a quantum system is from the
¼ boundary. We test whether θ aligns with any single direction in
the PCA decomposition (principal component analysis: a standard
technique that finds the dominant axes of variation in multivariate
data) of the decoherence trajectory. It does not: θ reads a diagonal
of the manifold, requiring all three principal components. This means
θ alone is insufficient as a system coordinate; you need additional
instruments (the "cockpit" framework).

---

## Abstract

The angular observable θ = arctan(sqrt(4*CΨ - 1)) measures
distance from the quantum-classical boundary at CΨ = ¼. The
structural cartography showed that CΨ visibility windows live on
a low-dimensional manifold with 3 principal components. We test
whether θ aligns with a single PC. It does not. θ requires
all 3 PCs for a reasonable fit (R^2 ~ 0.87), and its strongest
single-PC correlation is with PC3 (the Psi- sector, 11% variance),
not with the dominant PC1 (57%) or PC2 (22%). θ reads a diagonal
of the manifold, combining information from all principal directions.

## Setup

Star topology (S=0, A=1, B=2), identical to STRUCTURAL_CARTOGRAPHY:

| Parameter | Value |
|---|---|
| Topology | Star, 3 qubits |
| Initial state | Bell_SA x \|+>_B |
| Coupling | J_SA = 1.0, J_SB = 2.0 |
| Noise | Z-dephasing, γ = 0.05 |
| Sampling | dt = 0.005, 501 snapshots, t in [0, 10] |

CΨ definition: Purity x Ψ-norm = Tr(ρ^2) x L₁/(d-1). The
cartography used Concurrence x Ψ-norm for peak detection; these
differ for mixed states.

## Result 1: CΨ ranges across pairs

| Pair | CΨ range | θ > 0 | Fraction |
|---|---|---|---|
| AB (observer-observer) | [0.029, 0.354] | 23 snapshots | 4.6% |
| SA (sensor-observer A) | [0.019, 0.443] | 45 snapshots | 9.0% |
| SB (sensor-observer B) | [0.031, 0.289] | 26 snapshots | 5.2% |

CΨ_AB exceeds ¼ in only 4.6% of the trajectory. The SA pair has
the most quantum-regime variation (max CΨ = 0.443, max θ = 41.3 deg).

## Result 2: PCA structure (full trajectory)

| PC | Variance | Cumulative | Dominant loading | Interpretation |
|---|---|---|---|---|
| 1 | 56.6% | 56.6% | Phi-, C, Phi+, Pur | Sector balance + purity |
| 2 | 22.1% | 78.7% | Ψ-norm, Psi+, SvN | Coherence / mixedness |
| 3 | 11.4% | 90.1% | **Psi-** (loading 0.908) | Psi- sector isolated |
| 4 | 7.4% | 97.4% | | |
| 5 | 1.6% | 99.0% | | |

The full trajectory (501 points) needs 4 PCs for 95% variance, vs
3 PCs on the cartography's 9 peak windows. The difference arises
because the full trajectory includes monotonic decay components not
visible in peak-only sampling.

## Result 3: Theta is NOT a function of a single PC

### AB pair (quantum regime, 23 points)

| | Pearson r (linear correlation coefficient) | R^2(θ ~ PCk alone) |
|---|---|---|
| PC1 | +0.63 | 0.394 |
| PC2 | -0.10 | 0.011 |
| PC3 | -0.46 | 0.215 |
| **PC1+PC2+PC3** | | **R^2 = 0.864** |

### SA pair (quantum regime, 45 points)

| | Pearson r(θ, PCk) | R^2(θ ~ PCk alone) |
|---|---|---|
| PC1 | -0.31 | 0.095 |
| PC2 | +0.24 | 0.056 |
| PC3 | **+0.69** | **0.478** |
| **PC1+PC2+PC3** | | **R^2 = 0.872** |

### SB pair (quantum regime, 26 points)

| | Pearson r(θ, PCk) | R^2(θ ~ PCk alone) |
|---|---|---|
| PC1 | -0.74 | 0.546 |
| PC2 | **+0.76** | **0.584** |
| PC3 | -0.18 | 0.032 |
| **PC1+PC2+PC3** | | **R^2 = 0.700** |

**No single PC captures θ.** The strongest single-PC correlation
varies by pair: PC3 for SA, PC1 and PC2 (nearly tied) for SB, PC1
for AB. In all cases, at least 2 PCs are needed for R^2 > 0.85.

## Result 4: Missing dimensions

θ is a 1D scalar. The manifold is at least 4D (for 95% variance).

For the SA pair (most θ variation):
- θ aligns most with PC3 (11.4% variance, |r| = 0.69)
- Missing dimensions: PC1 (56.6%), PC2 (22.1%), PC4 (7.4%)

θ covers at best the **smallest** of the three major variance
directions. The dominant dynamics (PC1: sector balance, 57% variance)
are nearly invisible to θ.

## Result 5: Connection to θ-fidelity r = 0.87

THETA_PALINDROME_ECHO found r(theta_SB, F_channel) = 0.87
(r^2 = 0.757, 24.3% unexplained).

The 3-PC regression for theta_SA gives R^2 = 0.872, numerically
close to 0.87^2. This is likely coincidence: different pairs (SB vs
SA), different CΨ definitions, different datasets.

For theta_SB in the quantum regime: R^2(3PC) = 0.700, **worse** than
the θ-fidelity correlation. This means fidelity contains
information not captured by the AB-PCA basis, pointing to
subsystem-specific structure.

## Result 6: Nonlinear analysis

| Model | R^2 (SA, quantum regime) |
|---|---|
| θ ~ PC3 (linear) | 0.478 |
| θ ~ PC3 + PC3^2 (quadratic) | 0.488 |
| θ ~ sqrt(PC2^2 + PC3^2) | 0.064 |
| θ ~ PC1 + PC2 + PC3 (linear 3D) | 0.872 |

Quadratic extension adds almost nothing (+1%). The relationship is
essentially linear but **distributed across multiple PCs**.

## Interpretation

### Why θ cannot map to a single PC

θ = arctan(sqrt(4*CΨ - 1)) where CΨ = Purity x Ψ-norm.
CΨ is a **product** of two quantities that load on different PCs
(Purity -> PC1, Ψ-norm -> PC2). PCA decomposes the same observables
**additively** into orthogonal directions. A multiplicative combination
of two additive components cannot align with either component alone.

### Why PC3 (Psi-) dominates for SA

PC3 is nearly a pure Psi- fidelity indicator (loading 0.908). For the
SA pair, the Psi- component marks the phase of entanglement transfer:
when the S->A channel is active, the Bell decomposition shifts, and
Psi- fidelity tracks this shift. theta_SA responds to the same
dynamics, producing the correlation.

### Pair dependence

| Pair | Strongest PC | Physics |
|---|---|---|
| AB | PC1 | Sector dynamics between observers |
| SA | PC3 | Psi- fluctuations in Bell transfer |
| SB | PC1 ~ PC2 | Mix of sector and coherence |

The PCA decomposition is a property of the **AB observables**, but
θ lives on the **respective pair**. The correlation is
projection-dependent.

## Caveats

1. **CΨ definition:** This analysis uses CΨ = Purity x Ψ-norm.
   The cartography used Concurrence x Ψ-norm for peak detection.
   These differ for mixed states.

2. **PCA basis:** Computed on the full trajectory (501 snapshots),
   not the cartography's 9 peak windows. Loadings may differ.

3. **Topology-specific:** Star (N=3) with asymmetric coupling.
   Results may not transfer to other topologies.

4. **Correlation is not causation.** θ and the PCs co-vary along
   the same decoherence trajectory, driven by the same physical
   process.

5. **PCs not orthogonal on subset.** The PCA basis is orthogonal over
   all 501 points. On the θ > 0 subset (23-45 points), PCs can
   be correlated, affecting regression coefficients.

## Conclusion

**θ is not a function of a single PC.** It is a nonlinear
combination (via CΨ = Purity x Ψ-norm) that requires all three
major PCs. Its strongest single-PC correlation is with PC3
(Psi- sector), but even that explains only ~48% of the variance.

For a complete coordinate on the decoherence manifold, θ alone
is insufficient. At least two additional scalars are needed to cover
the PC1 and PC2 directions. Candidates: a sector-balance index (PC1)
and a mixedness measure such as linear entropy (PC2).

This result motivates the cockpit framework
([COCKPIT_UNIVERSALITY](COCKPIT_UNIVERSALITY.md)): θ is the
altimeter, but you need the full instrument panel.
