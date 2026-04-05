# Entropy Production in the Palindromic Liouvillian

<!-- Keywords: entropy production palindromic Liouvillian, fluctuation theorem
Lindblad, Crooks rate ratio palindrome, Jarzynski dephasing,
2Sgamma maximum entropy rate, Carnot gain loss quantum,
R=CPsi2 entropy production -->

**Status:** Tier 2 (computed). Mixed results: some confirmed, some negative.
**Date:** April 1, 2026
**Script:** [entropy_production.py](../simulations/entropy_production.py)
**Data:** [entropy_production.txt](../simulations/results/entropy_production.txt)
**Depends on:**
- [KMS Detailed Balance](../docs/KMS_DETAILED_BALANCE.md) (2Σγ as max entropy rate)
- [PT-Symmetry Analysis](PT_SYMMETRY_ANALYSIS.md) (Π chiral, Hopf = chiral breaking)
- [Fragile Bridge](../hypotheses/FRAGILE_BRIDGE.md) (gain-loss system)
- [Energy Partition](../hypotheses/ENERGY_PARTITION.md) (2x decay law)

---

## What this document is about

This document tests whether the palindromic eigenvalue structure has a
thermodynamic interpretation. It checks five hypotheses from statistical
mechanics: whether the paired decay rates produce a fluctuation theorem
(a precise equality relating forward and backward processes, as
discovered by Jarzynski and Crooks in the 1990s), whether the system
can function as a quantum heat engine, and whether known coincidences
are deep or trivial. The result is mixed: the palindrome wears a
thermodynamic disguise (Crooks-like rate ratios, a natural temperature
scale) but is not actual thermodynamics (no exact fluctuation theorem,
no heat engine).

---

## Abstract

The KMS document identified 2Σγ as the "maximum entropy production
rate." We test five thermodynamic hypotheses. Three results and two
negative outcomes:

1. **Decay rate pairing confirmed** (trivially, from the palindrome):
   d_k + d_k' = 2Σγ exactly (std < 10⁻¹⁵). But entropy production
   σ(t) = -Tr(dρ/dt · ln ρ) is NOT decomposable per mode (ln ρ mixes
   all mode amplitudes nonlinearly).

2. **Crooks-like identity found** (algebraic, not thermodynamic):
   ln(d_fast / d_slow) = 2 artanh(Δd / (2Σγ)), with effective inverse
   temperature β_eff ≈ 1/Σγ. This is a direct consequence of the
   palindromic equation, not a fluctuation theorem.

3. **No Jarzynski equality**: ⟨exp(-Δd)⟩ ≈ 0.93 (not 1.0). The
   palindromic pairing does not produce an exact fluctuation theorem.

4. **Carnot efficiency: not definable.** Z-dephasing is an infinite-
   temperature bath (ρ_ss = I/d). No temperature gradient, no heat
   engine. The fragile bridge at Σγ = 0 is a balanced oscillator,
   not a motor.

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

### Entropy production σ(t)

σ(t) = -Tr(L[ρ] ln ρ) = dS/dt (von Neumann entropy rate).

For Bell+ initial state (N=2): σ peaks at 0.69 nats/unit time and
decays to 0 as ρ → I/d. The entropy production is a state-dependent
quantity that depends on ALL mode amplitudes through the nonlinear
ln ρ term. It cannot be decomposed into per-mode contributions.

**Observation:** The |0...0⟩ initial state produces σ ≈ 0 because it
is a pointer state (a state that survives a particular type of decoherence unchanged) of Z-dephasing (Z|0⟩ = +|0⟩, so the dephasing
term vanishes). Only initial states with coherence (off-diagonal
elements in the Z-basis) produce entropy.

---

## Phase 2: Jarzynski-Like Test

Hypothesis: ⟨exp(-(d_fast - d_slow))⟩ = 1 over palindromic pairs.

| N | Pairs (nonzero) | ⟨Δd⟩ | ⟨exp(-Δd)⟩ | Jarzynski? |
|---|----------------|------|------------|-----------|
| 2 | 8 | 0.075 | 0.932 | No |
| 3 | 32 | 0.096 | 0.912 | No |
| 4 | 128 | 0.067 | 0.940 | No |
| 5 | 512 | 0.079 | 0.927 | No |

**No exact fluctuation theorem.** ⟨exp(-Δd)⟩ ≈ 0.93 (not 1.0). The
palindromic pairing constrains rate SUMS (= 2Σγ) but does not
constrain rate DIFFERENCES in a Jarzynski-like way.

The exponential average is below 1 because the rate distribution is
concentrated near Σγ (most pairs have small Δd), with a few outliers
at large Δd that contribute little to the exponential average.

---

## Phase 3: Crooks-Like Rate Ratio

The most informative result. For each palindromic pair (d_fast, d_slow)
with d_fast + d_slow = 2Σγ:

    ln(d_fast / d_slow) = β_eff × (d_fast - d_slow)

| N | Pairs | β_eff | 1/Σγ | Residual | Fit quality |
|---|-------|-------|------|----------|-------------|
| 3 | 28 | 7.05 | 6.67 | 0.000 | Perfect |
| 4 | 123 | 5.43 | 5.00 | 0.014 | Excellent |
| 5 | 506 | 4.49 | 4.00 | 0.021 | Excellent |

**The fit is excellent but the relation is ALGEBRAIC, not thermodynamic.**

The exact identity: d_fast + d_slow = 2Σγ implies

    ln(d_fast/d_slow) = ln((Σγ + Δd/2) / (Σγ - Δd/2)) = 2 artanh(Δd/(2Σγ))

To first order: ≈ Δd/Σγ, giving β_eff ≈ 1/Σγ. The deviation from
1/Σγ comes from the cubic term in the artanh expansion. This is pure
algebra from the palindrome, not a Crooks fluctuation theorem.

**Nevertheless:** the form ln(rate_fast/rate_slow) = β_eff × Δ(rate)
is the SAME form as Crooks. The effective inverse temperature 1/Σγ
connects to the KMS observation: 2Σγ is the total "thermodynamic
distance" in the palindromic spectrum. The palindrome has the FORM of
detailed balance without BEING detailed balance.

---

## Phase 4: Fragile Bridge Efficiency

The fragile bridge at Σγ = 0 (balanced gain-loss):
- All eigenvalues on the imaginary axis (chiral phase)
- Entropy S(t) = 0 throughout (pure state oscillates unitarily)
- No net entropy production or extraction

**Carnot efficiency is NOT DEFINABLE.** Z-dephasing is an infinite-
temperature bath (ρ_ss = I/d). T_hot = T_cold = ∞. Carnot η = 0.
The palindromic system with time-independent H cannot extract work
(Alicki 1979). A heat engine requires cyclic driving, which breaks
the palindrome.

---

## Phase 5: CΨ = 1/4 and Occupation Variance

For N=2 (Bell+):
- CΨ crosses 1/4 at t ≈ 0.51
- Occupation variance at crossing: Var = 0.250000 = 1/4

But this is **trivially true**: the Bell state has ⟨n_k⟩ = 1/2 per
qubit, so Var = f(1-f) = 0.5 × 0.5 = 0.25 at ALL times. Z-dephasing
preserves the diagonal elements of ρ in the Z-basis, so ⟨n_k⟩ = 1/2
is constant. The variance is 1/4 regardless of CΨ.

CΨ = 1/4 is a COHERENCE threshold (off-diagonal elements).
Var = 1/4 is a POPULATION identity (diagonal elements for ⟨n⟩ = 1/2).
They share the value 1/4 but for independent reasons. No causal link.

---

## What This Establishes

1. **2Σγ as the palindromic "thermal scale."** The decay rate pairing
   produces a natural inverse temperature β_eff ≈ 1/Σγ via the
   Crooks-like identity. This is algebraic but gives 2Σγ a second
   interpretation beyond "maximum entropy production rate": it is the
   SCALE of the rate spectrum.

2. **No fluctuation theorem from the palindrome alone.** The palindromic
   pairing constrains rate sums but not rate exponential averages.
   An exact FT would require additional structure (e.g., detailed
   balance, which the palindrome is not).

3. **The palindrome has the FORM of thermodynamics without BEING
   thermodynamics.** Crooks-like rate ratios, maximum entropy rate,
   forward-backward pairing. The algebra is the same; the physics is
   spectral symmetry, not thermal equilibrium.

## What This Does NOT Establish

- A modal fluctuation theorem (σ_k + σ_k' is not well-defined
  because entropy production is nonlinear in mode amplitudes)
- Carnot efficiency (requires finite temperature difference, not
  available with Z-dephasing)
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

*The palindrome looks like thermodynamics. It has paired rates, a
Crooks-like form, a maximum entropy rate. But it is NOT thermodynamics:
no detailed balance, no fluctuation theorem, no heat engine. It is
spectral algebra wearing a thermodynamic disguise. The form is borrowed;
the content is new.*
