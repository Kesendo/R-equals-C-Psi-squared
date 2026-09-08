# Entropy-Rate and Spectral-Pairing Audit

<!-- Keywords: von Neumann entropy rate palindromic Liouvillian,
spectral rate pairing, algebraic rate ratio, dephasing dynamics,
R=CPsi2 entropy audit -->

**Status:** Tier 2 computational audit. The reported results are spectral-algebraic
or trajectory-specific; no thermodynamic interpretation is established.
**Date:** April 1, 2026
**Script:** [entropy_production.py](../simulations/entropy_production.py)
**Data:** [entropy_production.txt](../simulations/results/entropy_production.txt)
**Depends on:**
- [KMS Detailed Balance](../docs/KMS_DETAILED_BALANCE.md) (distinguishes the shifted symmetry from KMS/detailed balance)
- [PT-Symmetry Analysis](PT_SYMMETRY_ANALYSIS.md) (shifted spectral symmetry; full irreducible class open)
- [Fragile Bridge](../hypotheses/FRAGILE_BRIDGE.md) (gain-loss dynamics)
- [Energy Partition](../hypotheses/ENERGY_PARTITION.md) (filtered-spectrum audit and F8 range/centre distinction)

---

## What this document is about

This document audits whether the computed palindromic rate pairing is enough
to support thermodynamic language. It is not. The calculation contains no
work protocol, heat current, bath temperature, detailed-balance measure, or
forward/reverse trajectory ensemble. What survives is an exact rate-sum
identity, a logarithmic identity obtained from that sum, one trajectory's von
Neumann entropy derivative, an ad hoc exponential average, and an independent
quarter-valued population identity.

---

## Abstract

Five computed observations survive after removing the thermodynamic reading:

1. **Decay-rate pairing confirmed** (a restatement of the palindrome):
   d_k + d_k' = 2Σγ to the reported numerical precision. The trajectory
   quantity dS/dt = -Tr(dρ/dt · ln ρ) is not decomposable per Liouvillian
   eigenvalue because ln ρ depends nonlinearly on the full state.

2. **Algebraic rate-ratio identity:**
   ln(d_fast / d_slow) = 2 artanh(Δd / (2Σγ)), with effective inverse
   linear slope approximately 1/Σγ near Δd = 0. This is a direct
   rearrangement of the palindromic equation, not an inverse temperature.

3. **Ad hoc exponential average:** the script reports
   ⟨exp(-Δd)⟩ ≈ 0.93 rather than 1. This is not a Jarzynski test: Δd is a
   rate difference, not dimensionless work, and no protocol ensemble or
   thermodynamic weighting was defined.

4. **Fragile-bridge trajectory:** at the sampled Σγ = 0 operating point,
   the reported trajectory remains pure and has S(t)=0. This calculation
   defines neither reservoir temperatures nor heat and work, so no Carnot
   statement follows.

5. **CΨ = 1/4 and Var = 1/4: independent coincidence.** The occupation
   variance is 1/4 at the CΨ crossing, but trivially: Bell states have
   ⟨n_k⟩ = 1/2, so f(1-f) = 1/4 at ALL times. CΨ = 1/4 is a coherence
   threshold; Var = 1/4 is a population identity. Same number, different
   mechanism.

---

## Phase 1: Decay Rate Pairing

For Heisenberg chains N=2-4, uniform γ = 0.05:

| N | Eigenvalues | Pairs | Rate sum mean | Rate sum std | 2Σγ |
|---|-------------|-------|--------------|-------------|------|
| 2 | 16 | 8 | 0.200000 | 5.6e-16 | 0.200 |
| 3 | 64 | 32 | 0.300000 | 2.3e-15 | 0.300 |
| 4 | 256 | 128 | 0.400000 | 5.6e-15 | 0.400 |

Exact to machine precision. This is the palindromic equation
λ_k + λ_k' = -2Σγ restated as a rate sum.

### Von Neumann entropy rate along the sampled trajectory

dS/dt = -Tr(L[ρ] ln ρ), where this expression is defined on the support
of ρ (or by the corresponding limiting prescription).

For the Bell+ initial state at N=2, the script's largest sampled dS/dt is
0.69 nats per unit time at its first point, t=0.01, followed by a decrease on
that 50-point grid. This is not a continuum-time peak: the initially pure
state develops small eigenvalues proportional to t, so dS/dt grows logarithmically as t approaches zero from above.
This is a state-dependent dynamical reading that depends on all mode
amplitudes through the nonlinear ln ρ term. It cannot be decomposed into
independent per-mode contributions, and the run does not establish a unique
stationary state for every invariant sector.

**Observation:** The |0...0⟩ initial state produces dS/dt ≈ 0 here because
it is both a Z-dephasing pointer state and an eigenstate of this Hamiltonian:
the dissipator and commutator therefore vanish on it. Z-diagonality alone is
not enough. For example, |01⟩⟨01| initially has no Z-basis coherence, but the
exchange Hamiltonian generates coherence and subsequent mixing. Entropy growth
in the full H+D evolution can therefore begin from a Z-diagonal state that is
not Hamiltonian-invariant.

---

## Phase 2: Ad Hoc Exponential Average

The script evaluates ⟨exp(-(d_fast - d_slow))⟩ over all palindromic pairs,
including the stationary/max-rate endpoint pairs with d_slow = 0. The
positive-rate restriction belongs to the separate logarithmic identity in
Phase 3.
This expression is retained as a numerical transform of the rate differences,
not as a fluctuation-theorem observable: the exponent has not been derived as
dimensionless work and the average has no forward/reverse protocol weights.

| N | All palindrome pairs | ⟨Δd⟩ | ⟨exp(-Δd)⟩ | Jarzynski? |
|---|----------------|------|------------|-----------|
| 2 | 8 | 0.075 | 0.932 | No |
| 3 | 32 | 0.096 | 0.912 | No |
| 4 | 128 | 0.067 | 0.940 | No |
| 5 | 512 | 0.079 | 0.927 | No |

The computed average is approximately 0.93 rather than 1. The palindromic
pairing constrains rate sums (= 2Σγ), but does not fix this exponential
transform of the rate differences.

The exponential average is below 1 because the rate distribution is
concentrated near Σγ (most pairs have small Δd), with a few outliers
at large Δd that contribute little to the exponential average.

---

## Phase 3: Algebraic Rate-Ratio Identity

For each positive-rate palindromic pair (d_fast, d_slow) with
d_fast + d_slow = 2Σγ, the script fits

    ln(d_fast / d_slow) = b_fit × (d_fast - d_slow)

| N | Pairs | fitted slope b_fit | linear coefficient 1/Σγ | Residual | Fit description |
|---|-------|-------|------|----------|-------------|
| 3 | 28 | 7.05 | 6.67 | 0.000 | Perfect |
| 4 | 123 | 5.43 | 5.00 | 0.014 | Excellent |
| 5 | 506 | 4.49 | 4.00 | 0.021 | Excellent |

The nonlinear relation is algebraically fixed by the rate sum:

The exact identity: d_fast + d_slow = 2Σγ implies

    ln(d_fast/d_slow) = ln((Σγ + Δd/2) / (Σγ - Δd/2)) = 2 artanh(Δd/(2Σγ))

To first order this is Δd/Σγ. The fitted slope differs from 1/Σγ because
the fit compresses the cubic and higher artanh terms into one coefficient.
Neither coefficient is an effective inverse temperature without a separately
defined thermodynamic ensemble and detailed-balance relation.

---

## Phase 4: Fragile-Bridge Trajectory

At the sampled Σγ = 0 balanced gain-loss operating point:

- the initial state is the four-qubit GHZ/cat coherence
  `(|0000⟩+|1111⟩)/sqrt(2)`, not a Bell pair;

- the reported stable-regime eigenvalues lie on the imaginary axis;
- the simulated trajectory remains pure and S(t) = 0 throughout.

These are spectral and dynamical statements about that model run. No bath
temperatures, heat currents, work strokes, or cyclic protocol were defined,
so the calculation supports no efficiency or heat-engine verdict.

---

## Phase 5: CΨ = 1/4 and Occupation Variance

For N=2 (Bell+):
- the sampled grid brackets the crossing between t=0.5124 and t=0.7636;
  root refinement gives CΨ=1/4 at t ≈ 0.747
- Occupation variance at crossing: Var = 0.250000 = 1/4

But this is **trivially true**: the Bell state has ⟨n_k⟩ = 1/2 per
qubit, so Var = f(1-f) = 0.5 × 0.5 = 0.25 at ALL times. In this run,
both Z-dephasing and the Hamiltonian preserve this Bell+ population: Bell+ is
an eigenstate of the two-qubit Heisenberg Hamiltonian. Z-dephasing alone does
not justify population constancy for a generic initial state under H+D. Here
⟨n_k⟩ = 1/2 is constant, so the variance is 1/4 regardless of CΨ and no
variance-peak time is defined.

CΨ = 1/4 is a COHERENCE threshold (off-diagonal elements).
Var = 1/4 is a POPULATION identity (diagonal elements for ⟨n⟩ = 1/2).
They share the value 1/4 but for independent reasons. No causal link.

---

## What This Establishes

1. **2Σγ as a spectral pair sum.** For the tested palindromic spectra,
   paired decay rates add to 2Σγ. This is the width/shift parameter in the
   algebraic reflection, not a thermodynamic production rate.

2. **A rate-ratio identity and its linearization.** The artanh formula and
   the coefficient 1/Σγ follow directly from d_fast+d_slow=2Σγ. They do
   not define temperature or detailed balance.

3. **Trajectory entropy is a separate nonlinear observable.** The N=2 Bell+
   run has the reported largest sampled dS/dt at the grid's first point, while the Z-basis pointer-state control
   has approximately zero entropy change. Neither value is a modewise
   consequence of the spectral pair sum.

4. **The quarter-valued quantities are independent.** At the tested Bell+
   trajectory, Var=1/4 is fixed by the populations while CΨ crosses 1/4
   through its coherence dynamics.

## What This Does NOT Establish

- A modal entropy-production or fluctuation theorem (a modewise entropy
  contribution is not defined by these calculations)
- Any reservoir temperature, heat current, work extraction, or Carnot
  efficiency
- A unique maximally mixed stationary state across the full block-decomposed
  Liouvillian
- A deep CΨ-Fermi connection (the numerical coincidence is trivial)
- Whether a time-dependent driving protocol could extract work from
  the palindromic structure (not tested)

---

## Reproducibility

| Component | Location |
|-----------|----------|
| Script | [`simulations/entropy_production.py`](../simulations/entropy_production.py) |
| Output | [`simulations/results/entropy_production.txt`](../simulations/results/entropy_production.txt) |

---

## References

- Esposito, M., Lindenberg, K., Van den Broeck, C. (2010). "Entropy
  production as correlation between system and reservoir." NJP 12, 013013.
- Jarzynski, C. (1997). "Nonequilibrium equality for free energy
  differences." PRL 78, 2690.
- Crooks, G.E. (1999). "Entropy production fluctuation theorem and
  the nonequilibrium work relation for free energy differences."
  PRE 60, 2721.
- Alicki, R. (1979). "The quantum open system as a model of the heat
  engine." J. Phys. A 12, L103.

---

*The surviving result is spectral algebra plus trajectory dynamics: an exact
pair-rate sum, its logarithmic rearrangement, and state-dependent entropy
change. Thermodynamic meaning would require additional physical definitions
and calculations not present in this experiment.*
