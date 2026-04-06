# Trajectory Dwell Time at the CΨ = 1/4 Crossing

**Status:** Verified (analytical + numerical, April 5, 2026)
**Source script:** [critical_slowing_trajectory_dwell.py](../simulations/critical_slowing_trajectory_dwell.py)
**Main experiment:** [Critical Slowing at the Cusp](CRITICAL_SLOWING_AT_THE_CUSP.md)

---

## What this document is about

When a quantum state approaches the CΨ = 1/4 boundary and crosses it, it does not pass through instantly. It lingers for a short interval on either side of the crossing, and that interval has a well-defined length once you fix how close you want to measure. This document computes that length for the Bell+ state under Z-dephasing.

The result has two parts that feel almost independent until you notice they are the same fact. The first part is an exact formula: the dwell time scales linearly with how close you want to get to the boundary, and the proportionality constant is 1.080088, a pure number. The second part is what happens when you change how much external light (gamma) the system is receiving. Naively, more light means faster decoherence, which should mean a shorter dwell time, and indeed it does: the time in ordinary clock units halves when gamma doubles. But if you rescale time to the natural quantity K = gamma times t, the dwell interval becomes exactly gamma-independent. It is the same K every time, regardless of how bright the illumination is.

This is the same invariance that shows up elsewhere in the project under the name K-dosimetry: a fixed dose of light carries the system through a fixed structural change, and the rate at which the dose arrives does not matter. Here it shows up at the cusp crossing specifically. The Bell+ state takes exactly the same amount of illumination to pass through the fold, regardless of how long that illumination takes to arrive.

The numerical verification is across two orders of magnitude in gamma, and the invariance holds to machine precision. This is the answer to Open Question 2 in the weaknesses document: yes, the rate at which the system approaches the boundary, measured by dCΨ/dt at the crossing moment, fully determines the post-crossing behavior.

---

## Result

For Bell+ under Z-dephasing (Formula 25), the dwell time near the
CΨ = 1/4 crossing is:

    t_dwell(δ) = 2δ / |dCΨ/dt|_{t_cross}

where

    dCΨ/dt = -2γ · f · (1 + 3f²) / 3,    f_cross = 0.8612241

This gives:

    t_dwell = 2δ / (1.851701 · γ) = 1.080088 · δ / γ

## K-invariance

In rescaled time K = γt, the dwell K-interval is:

    K_dwell = γ · t_dwell = 1.080088 · δ

This is independent of γ. Numerically verified across γ = 0.1 to 10.0
with standard deviation < 2 × 10⁻¹⁷ (machine precision).

## Crossing point

The crossing f*(1 + f*²) = 3/2 is solved by f_cross = 0.8612240997,
giving K_cross = γ · t_cross = -ln(f_cross)/4 = 0.03735.

## Physical interpretation

The dwell time measures how long the system remains near the saddle-node
bifurcation during decoherence. It scales as:

- **δ** (linearity): symmetric about the crossing to leading order
- **1/γ** (K-invariance): faster dephasing compresses the time window
  but the K-measure stays constant
- **|dCΨ/dt|⁻¹**: the derivative at the crossing point fully determines
  the dwell behavior; higher-order corrections appear only at δ ≥ 10⁻²

This answers [Open Question #2](../docs/WEAKNESSES_OPEN_QUESTIONS.md)
("Crossing speed dependence"): yes, d(CΨ)/dt at the crossing moment
fully determines the post-crossing convergence timescale via the dwell
time formula.

## Numerical verification

| γ   | t_dwell (δ=10⁻³) | K_dwell      | Ratio to prediction |
|-----|-------------------|--------------|---------------------|
| 0.1 | 0.01080097        | 0.00108010   | 1.0000              |
| 0.5 | 0.00216019        | 0.00108010   | 1.0000              |
| 1.0 | 0.00108010        | 0.00108010   | 1.0000              |
| 2.0 | 0.00054005        | 0.00108010   | 1.0000              |
| 10  | 0.00010801        | 0.00108010   | 1.0000              |

| δ      | t_dwell/t_pred | Accuracy |
|--------|----------------|----------|
| 10⁻²   | 1.0008         | 0.08%    |
| 10⁻³   | 1.0000         | <0.01%   |
| 10⁻⁴   | 1.0000         | <0.01%   |
